//! Soufflerie 3D : Lattice Boltzmann D3Q19 sur GPU (wgpu). Portage de
//! `js/lbm3d.js` : le WGSL est engendré ici exactement comme dans
//! l'application web.
//!
//! Un pas = une passe de calcul : propagation « pull » et collision fusionnées.
//! Populations stockées décalées de leur poids, g = f − w, en demi-précision
//! (deux par mot de 32 bits, pack2x16float), ou en simple précision au choix.
//! Collision régularisée (projection d'ordre 2) avec Smagorinsky aux Reynolds
//! de maille élevés. Parois : rebond interpolé de Bouzidi (fractions q des
//! liens calculées sur le CPU), efforts par échange de quantité de mouvement
//! accumulés par maille pariétale. Entrée à vitesse imposée (équilibre),
//! sortie à gradient nul avec zone éponge, bords latéraux périodiques.

use std::fmt::Write as _;
use std::sync::mpsc;

use rayon::prelude::*;
use wgpu::util::DeviceExt;

use crate::bodies3d::Body3D;
use crate::cpu::Ramp;

pub const C3: [[i32; 3]; 19] = [
    [0, 0, 0],
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1],
    [1, 1, 0],
    [-1, -1, 0],
    [1, 0, 1],
    [-1, 0, -1],
    [0, 1, 1],
    [0, -1, -1],
    [1, -1, 0],
    [-1, 1, 0],
    [1, 0, -1],
    [-1, 0, 1],
    [0, 1, -1],
    [0, -1, 1],
];

pub(crate) fn w3(i: usize) -> f64 {
    let c = C3[i];
    match c[0].abs() + c[1].abs() + c[2].abs() {
        0 => 1.0 / 3.0,
        1 => 1.0 / 18.0,
        _ => 1.0 / 36.0,
    }
}

fn opp3(i: usize) -> usize {
    let c = C3[i];
    C3.iter().position(|d| d[0] == -c[0] && d[1] == -c[1] && d[2] == -c[2]).unwrap()
}

/// Populations en demi-précision rangées une par case de 16 bits
/// (`array<u16>`, plan de N cases par direction) au lieu de deux par mot de
/// 32 bits : chaque population tirée d'un voisin est une lecture de 2 octets,
/// contiguë d'un fil à l'autre (mot de 4 octets à moitié utile sinon). Mêmes
/// conversions, mêmes valeurs au bit près. Fixé à l'ouverture du périphérique
/// (`SHADER_I16`, VK_KHR_16bit_storage) ; SOUFFLERIE_U16=0 garde les mots.
static POP16: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

pub(crate) fn pop16() -> bool {
    POP16.load(std::sync::atomic::Ordering::Relaxed)
}

/// Demi-mot de chaque direction dans le rangement en mots de 32 bits (deux
/// populations par mot) : directions appariées tirées de mailles voisines en
/// x (même rangée, fils voisins), puis les autres. SOUFFLERIE_PAIRES=0 :
/// ordre des directions (ancien rangement).
pub(crate) fn hslot() -> [usize; 19] {
    const PAIRES: [usize; 19] = [18, 0, 1, 10, 11, 12, 13, 2, 4, 6, 8, 14, 15, 5, 3, 9, 7, 16, 17];
    if std::env::var("SOUFFLERIE_PAIRES").is_ok_and(|v| v == "0") {
        std::array::from_fn(|i| i)
    } else {
        PAIRES
    }
}

/// Direction rangée dans chaque demi-mot (None : remplissage).
pub(crate) fn hdir() -> [Option<usize>; 20] {
    let mut d = [None; 20];
    for (i, &h) in hslot().iter().enumerate() {
        d[h] = Some(i);
    }
    d
}

/// Constante WGSL `HS` (demi-mot de chaque direction).
pub(crate) fn hs_wgsl() -> String {
    format!("const HS = array<u32, 19>({});\n", hslot().map(|h| format!("{h}u")).join(", "))
}

/// En-tête WGSL des shaders qui accèdent aux populations.
pub(crate) fn pop16_enable(p16: bool) -> &'static str {
    if p16 {
        "enable wgpu_int16;
"
    } else {
        ""
    }
}

pub(crate) const SOLID: u32 = 0xffff_ffff;
/// Maille d'une grille parente recouverte par une grille fine (voir
/// `Gpu3D::mark_covered`).
pub(crate) const COVERED: u32 = 0xffff_fffd;

/// Mise au point (SOUFFLERIE_TEMPS=1) : durée d'une étape de la mise en place.
pub(crate) fn temps(label: &str, t0: std::time::Instant) {
    if std::env::var("SOUFFLERIE_TEMPS").is_ok() {
        eprintln!("  [mise en place] {label} : {:.3} s", t0.elapsed().as_secs_f64());
    }
}

/// Module de calcul sans les contrôles que wgpu ajoute par défaut au code
/// engendré : bornes de chaque accès aux tampons (index ramené dans le
/// tableau), divisions entières protégées (fonctions d'enrobage pour chaque
/// `%` et `/`) et compteurs de boucle. Les shaders de la soufflerie
/// n'indexent que dans leurs tableaux (index calculés sur les dimensions des
/// grilles) et ne divisent que par des dimensions non nulles : mêmes
/// résultats, au bit près, pour bien moins d'instructions par maille.
/// SOUFFLERIE_CONTROLES=1 rétablit les contrôles (diagnostic).
pub(crate) fn compute_module(device: &wgpu::Device, label: &str, code: String) -> wgpu::ShaderModule {
    compute_module_raw(device, label, wgsl_hook(label, code))
}

/// Mise au point : SOUFFLERIE_WGSL=dossier écrit chaque shader dans
/// `dossier/label.wgsl`, ou le relit s'il y est déjà (essais sans recompiler).
/// Collision : code d'avant la spécialisation par grille.
fn wgsl_hook(label: &str, code: String) -> String {
    let mut code = code;
    if let Ok(dir) = std::env::var("SOUFFLERIE_WGSL") {
        let path = std::path::Path::new(&dir).join(format!("{label}.wgsl"));
        match std::fs::read_to_string(&path) {
            Ok(c) => code = c,
            Err(_) => {
                let _ = std::fs::write(&path, &code);
            }
        }
    }
    code
}

fn compute_module_raw(device: &wgpu::Device, label: &str, code: String) -> wgpu::ShaderModule {
    let desc = wgpu::ShaderModuleDescriptor {
        label: Some(label),
        source: wgpu::ShaderSource::Wgsl(code.into()),
    };
    if std::env::var("SOUFFLERIE_CONTROLES").is_ok_and(|v| v == "1") {
        return device.create_shader_module(desc);
    }
    // SAFETY : voir ci-dessus ; les index et diviseurs sont valides par
    // construction (dimensions des grilles, tables de blocs et listes
    // calculées sur le processeur).
    unsafe { device.create_shader_module_trusted(desc, wgpu::ShaderRuntimeChecks::unchecked()) }
}

/// Littéral flottant WGSL.
pub(crate) fn f(v: f64) -> String {
    if v.fract() == 0.0 {
        format!("{v:.1}")
    } else {
        format!("{v:?}")
    }
}

/// Somme Σ k(i)·fv[i] (k entier).
fn sum_expr(k: impl Fn(usize) -> i32) -> String {
    let mut t: Vec<String> = Vec::new();
    for i in 0..19 {
        match k(i) {
            0 => {}
            1 => t.push(format!("fv[{i}]")),
            -1 => t.push(format!("-fv[{i}]")),
            n => t.push(format!("{n} * fv[{i}]")),
        }
    }
    if t.is_empty() {
        "0.0".into()
    } else {
        t.join(" + ").replace("+ -", "- ")
    }
}

const PARAMS_WGSL: &str = "
struct Params {
  nx: u32, ny: u32, nz: u32, N: u32,
  tau0: f32, tauSponge: f32, spongeStart: u32, uIn: f32,
  accel: f32, c18: f32, writeMacro: u32, les: u32,
  rows: u32, sa: u32, nu: f32, lesA: f32,
  xtr: f32, band: f32, essai: u32, essaiNu: f32,
  bnx: u32, bny: u32, bnz: u32, bcoord: u32,
};
";

/// Rang d'une maille et position d'un rang : grille boîte (P.bnx = 0) ou
/// stockage creux par blocs de 8³ (table `blk`, puis coordonnées des blocs).
const SPARSE_FNS: &str = "
fn id_of(x: u32, y: u32, z: u32) -> u32 {
  if (P.bnx == 0u) { return x + P.nx * (y + P.ny * z); }
  let s = blk[(x >> 3u) + P.bnx * ((y >> 3u) + P.bny * (z >> 3u))];
  return s * 512u + (x & 7u) + 8u * ((y & 7u) + 8u * (z & 7u));
}
fn pos_of(i: u32) -> vec3<i32> {
  if (P.bnx == 0u) { return vec3<i32>(i32(i % P.nx), i32((i / P.nx) % P.ny), i32(i / (P.nx * P.ny))); }
  let c = blk[P.bcoord + (i >> 9u)];
  let l = i & 511u;
  return vec3<i32>(i32((c & 1023u) * 8u + (l & 7u)), i32(((c >> 10u) & 1023u) * 8u + ((l >> 3u) & 7u)), i32((c >> 20u) * 8u + (l >> 6u)));
}
";

/// Zone laminaire imposée : abscisse locale rapportée à la corde de la
/// tranche d'envergure y (étendue en x des mailles du corps dans la tranche).
const LAM_FNS: &str = "
fn lam_xl(x: i32, y: i32) -> vec2<f32> {
  let e = slab[u32(y)];
  return vec2<f32>(f32(x) - e.x, e.y - e.x);
}
";

/// Spalart–Allmaras : viscosité turbulente et loi de Spalding (comme en 2D).
const SA_FNS: &str = "
const KAPPA = 0.41;
// Hauteur de raccord de la loi de paroi (mailles) : au-delà de la première
// couche de nœuds pariétaux pour l'interpolation trilinéaire.
const WM_H = 2.0;
const CV1_3 = 357.911;
fn nu_turb(v: f32) -> f32 {
  if (v <= 0.0) { return 0.0; }
  let chi = v / P.nu;
  let c3 = chi * chi * chi;
  return v * c3 / (c3 + CV1_3);
}
// u_τ de la loi de Spalding passant par (y, u_t) (Newton sur u⁺).
fn spalding_utau(ut: f32, y: f32) -> f32 {
  if (!(ut > 1e-12)) { return 0.0; }
  let e_kb = 0.12873490;
  let rey = y * ut / P.nu;
  var up = sqrt(rey);
  if (rey >= 30.0) { up = log(rey) / KAPPA; }
  for (var k = 0; k < 20; k++) {
    let ku = KAPPA * up;
    let e = exp(ku);
    let yy = up + e_kb * (e - 1.0 - ku - ku * ku / 2.0 - ku * ku * ku / 6.0);
    let dy = 1.0 + e_kb * KAPPA * (e - 1.0 - ku - ku * ku / 2.0);
    let h = up * yy - rey;
    let next = up - h / (yy + up * dy);
    up = select(0.5 * up, next, next > 0.0);
    if (abs(h) < 1e-6 * rey) { break; }
  }
  return ut / up;
}
// Gradient de pression longitudinal (essai 8388608) : loi de Townsend (1961)
// et Mellor (1966). Couche interne d'équilibre, contrainte totale
// τ(y) = u_τ² + P y (P = (1/ρ) dp/ds, cinématique, s le long du glissement)
// et longueur de mélange κ y : du/dy = √τ / (κ y). Intégrée et raccordée à la
// loi logarithmique là où P y ≪ u_τ², elle ajoute à la loi sans gradient
//   Δ(y) = (2/κ) [√τ − u_τ + u_τ ln(2 u_τ / (√τ + u_τ))],
// nul pour P = 0 (Spalding inchangé), ≈ P y / (2κ u_τ) aux faibles
// gradients, et (2/κ)√(P y) à u_τ = 0 (profil de Stratford au décollement).
// Gradient favorable : τ borné à 0 (au-delà, couche sans contrainte).
fn mellor_du(u: f32, pg: f32, y: f32) -> f32 {
  if (!(u > 1e-12)) { return 2.0 / KAPPA * sqrt(max(pg * y, 0.0)); }
  let r = max(pg * y / (u * u), -1.0);
  // a = √(1+r) − 1 sans soustraction ; Δκ/(2u) = a − ln(1 + a/2).
  let a = r / (sqrt(1.0 + r) + 1.0);
  var g = a - log(1.0 + 0.5 * a);
  if (abs(a) < 1e-2) { g = a * (0.5 + a * (0.125 - a / 24.0)); }
  return 2.0 / KAPPA * u * g;
}
// ∫₀^y max(τ_w + P y', 0) dy'.
fn tau_int(tw: f32, pg: f32, y: f32) -> f32 {
  if (pg < 0.0 && tw + pg * y < 0.0) { return tw * tw / (-2.0 * pg); }
  return tw * y + 0.5 * pg * y * y;
}
// u_τ de la loi avec gradient : racine de F(u) = u − u_τ,Spalding(u_t − Δ(u))
// par fausse position (Illinois) dans un encadrement. Gradient adverse :
// F(0) ≤ 0 ≤ F(u_τ,0) (Spalding sans gradient). Favorable : racine au-delà
// de u_τ,0, encadrée par 2 u_τ,0 sauf sous-couche visqueuse (borne gardée).
fn mellor_utau(ut: f32, y: f32, pg: f32) -> f32 {
  let u0 = spalding_utau(ut, y);
  if (pg == 0.0 || !(u0 > 1e-12)) { return u0; }
  var lo = 0.0;
  var hi = u0;
  if (pg < 0.0) { lo = u0; hi = 2.0 * u0; }
  var flo = lo - spalding_utau(max(0.0, ut - mellor_du(lo, pg, y)), y);
  var fhi = hi - spalding_utau(max(0.0, ut - mellor_du(hi, pg, y)), y);
  if (fhi <= 0.0) { return hi; }
  if (flo >= 0.0) { return lo; }
  var side = 0;
  for (var k = 0; k < 12; k++) {
    if (!(fhi - flo > 0.0)) { break; }
    let m = clamp((lo * fhi - hi * flo) / (fhi - flo), lo, hi);
    let fm = m - spalding_utau(max(0.0, ut - mellor_du(m, pg, y)), y);
    if (abs(fm) < 1e-5 * u0) { return m; }
    if (fm < 0.0) {
      lo = m; flo = fm;
      if (side == -1) { fhi *= 0.5; }
      side = -1;
    } else {
      hi = m; fhi = fm;
      if (side == 1) { flo *= 0.5; }
      side = 1;
    }
  }
  if (!(fhi - flo > 0.0)) { return lo; }
  return clamp((lo * fhi - hi * flo) / (fhi - flo), lo, hi);
}
";

/// Transport explicite de ν̃ (même schéma que sa.wgsl, en 3D) : advection
/// décentrée amont, diffusion, production et destruction ; paroi : ν̃ = κ u_τ y.
const SA3_WGSL: &str = "
@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read> macroBuf: array<vec4<f32>>;
@group(0) @binding(2) var<storage, read> wall: array<u32>;
@group(0) @binding(3) var<storage, read> wallGeo: array<vec4<f32>>;
@group(0) @binding(4) var<storage, read> utau: array<f32>;
@group(0) @binding(5) var<storage, read> dist: array<f32>;
@group(0) @binding(6) var<storage, read> nutIn: array<f32>;
@group(0) @binding(7) var<storage, read_write> nutOut: array<f32>;
@group(0) @binding(8) var<storage, read> slab: array<vec2<f32>>;
fn lam_xl(x: i32, y: i32) -> vec2<f32> {
  let e = slab[u32(y)];
  return vec2<f32>(f32(x) - e.x, e.y - e.x);
}

const CB1 = 0.1355;
const CB2 = 0.622;
const SIGMA = 0.6666667;
const KAPPA = 0.41;
// Hauteur de raccord de la loi de paroi (mailles) : au-delà de la première
// couche de nœuds pariétaux pour l'interpolation trilinéaire.
const WM_H = 2.0;
const CW2 = 0.3;
const CW3_6 = 64.0;
const CV1_3 = 357.911;
const NU_TILDE_MAX = 0.15;

@group(0) @binding(9) var<storage, read> blk: array<u32>;
@group(0) @binding(10) var<storage, read> plain: array<u32>;
SPARSE_FNS

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let idx = gid.x + gid.y * nwg.x * 64u;
  if (idx >= P.N) { return; }
  let pp = pos_of(idx);
  let x = pp.x;
  let y = pp.y;
  let z = pp.z;
  // Masque des mailles fluides ordinaires (un bit) : le type complet
  // (4 octets) n'est lu que pour les autres.
  var wk0 = 0u;
  if (((plain[idx >> 5u] >> (idx & 31u)) & 1u) == 0u) { wk0 = wall[idx]; }
  // Maille recouverte par une grille fine : ν̃ recopié de celle-ci après ses pas.
  if (wk0 == 0xfffffffdu) { return; }
  let nu = P.nu;
  var out: f32;
  if (wk0 == 0xffffffffu) {
    out = 0.0;
  } else if (wk0 == 0xfffffffeu) {
    // Maille fantôme : ν̃ venu de la grille parente.
    out = nutIn[idx];
  } else if (x == 0) {
    out = 3.0 * nu;
  } else if (x == i32(P.nx) - 1) {
    out = nutIn[idx - 1u];
  } else if (wk0 != 0u && wallGeo[wk0 - 1u].w != 0.0) {
    // Premier nœud : ν̃ = κ u_τ h, viscosité turbulente de la couche
    // pariétale [0, h] de la loi de paroi (hauteur de raccord fixe : la
    // valeur ne saute pas quand la paroi traverse une rangée de nœuds).
    let k = wk0 - 1u;
    out = min(NU_TILDE_MAX, KAPPA * utau[4u * k] * WM_H);
    // Diagnostic : ν̃ = κ u_τ y_n, distance du nœud (524288).
    if ((P.essai & 524288u) != 0u) { out = min(NU_TILDE_MAX, KAPPA * utau[4u * k] * max(0.1, abs(wallGeo[k].w))); }
    // Gradient de pression (8388608) : viscosité de la couche pariétale de
    // la loi de Mellor, κ h √(u_τ² + P h), écrite par la passe de paroi.
    if ((P.essai & 8388608u) != 0u && (P.essai & 16777216u) == 0u) {
      let nw = utau[arrayLength(&utau) / 6u * 5u + k];
      if (nw >= 0.0) { out = min(NU_TILDE_MAX, nw); }
    }
    if ((P.essai & 1024u) != 0u) { out = min(out, P.essaiNu); }
  } else {
    // Voisins (comme `cell`, sans modulo ; ici 0 < x < nx − 1) : y
    // périodique, z périodique ou borné (sol).
    let xc = u32(x);
    let yc = u32(y);
    let zc = u32(z);
    let ground = (P.rows & 2u) != 0u;
    let ym = select(yc - 1u, P.ny - 1u, y == 0);
    let yp = select(yc + 1u, 0u, yc + 1u == P.ny);
    let zm = select(zc - 1u, select(P.nz - 1u, 0u, ground), z == 0);
    let zp = select(zc + 1u, select(0u, zc, ground), zc + 1u == P.nz);
    let ixm = id_of(xc - 1u, yc, zc);
    let ixp = id_of(xc + 1u, yc, zc);
    let iym = id_of(xc, ym, zc);
    let iyp = id_of(xc, yp, zc);
    let izm = id_of(xc, yc, zm);
    let izp = id_of(xc, yc, zp);
    let v = nutIn[idx];
    let vxm = nutIn[ixm];
    let vxp = nutIn[ixp];
    let vym = nutIn[iym];
    let vyp = nutIn[iyp];
    let vzm = nutIn[izm];
    let vzp = nutIn[izp];
    let u = macroBuf[idx].xyz;
    var adv = select(u.x * (vxp - v), u.x * (v - vxm), u.x > 0.0);
    adv += select(u.y * (vyp - v), u.y * (v - vym), u.y > 0.0);
    adv += select(u.z * (vzp - v), u.z * (v - vzm), u.z > 0.0);
    if ((P.essai & 1048576u) != 0u) {
      // Diagnostic : advection de Lax–Wendroff (ordre 2, diffusion
      // numérique u_a²/2 au lieu de |u_a|/2 par direction).
      adv = 0.5 * (u.x * (vxp - vxm) + u.y * (vyp - vym) + u.z * (vzp - vzm))
        - 0.5 * (u.x * u.x * (vxp - 2.0 * v + vxm) + u.y * u.y * (vyp - 2.0 * v + vym) + u.z * u.z * (vzp - 2.0 * v + vzm));
    }
    var diff = (nu + 0.5 * (v + vxp)) * (vxp - v) - (nu + 0.5 * (v + vxm)) * (v - vxm)
      + (nu + 0.5 * (v + vyp)) * (vyp - v) - (nu + 0.5 * (v + vym)) * (v - vym)
      + (nu + 0.5 * (v + vzp)) * (vzp - v) - (nu + 0.5 * (v + vzm)) * (v - vzm);
    let g = 0.5 * vec3<f32>(vxp - vxm, vyp - vym, vzp - vzm);
    diff = (diff + CB2 * dot(g, g)) / SIGMA;
    let uxm = macroBuf[ixm].xyz;
    let uxp = macroBuf[ixp].xyz;
    let uym = macroBuf[iym].xyz;
    let uyp = macroBuf[iyp].xyz;
    let uzm = macroBuf[izm].xyz;
    let uzp = macroBuf[izp].xyz;
    let w = 0.5 * vec3<f32>(
      (uyp.z - uym.z) - (uzp.y - uzm.y),
      (uzp.x - uzm.x) - (uxp.z - uxm.z),
      (uxp.y - uxm.y) - (uyp.x - uym.x));
    let omega = length(w);
    let d = max(0.5, dist[idx]);
    let kd2 = KAPPA * KAPPA * d * d;
    let chi = v / nu;
    let c3 = chi * chi * chi;
    let fv1 = c3 / (c3 + CV1_3);
    let fv2 = 1.0 - chi / (1.0 + chi * fv1);
    let s_tilde = max(omega + v / kd2 * fv2, 0.3 * omega);
    var prod = CB1 * s_tilde * v;
    if (P.xtr > 0.0 && dist[idx] < P.band) {
      let xl = lam_xl(x, y);
      if (xl.y > 0.0 && xl.x < P.xtr * xl.y) { prod = 0.0; }
    }
    var r = 10.0;
    if (s_tilde > 1e-12) { r = min(10.0, v / (s_tilde * kd2)); }
    let r2 = r * r;
    let gg = r + CW2 * (r2 * r2 * r2 - r);
    let g2 = gg * gg;
    let fw = gg * pow((1.0 + CW3_6) / (g2 * g2 * g2 + CW3_6), 1.0 / 6.0);
    let cw1 = CB1 / (KAPPA * KAPPA) + (1.0 + CB2) / SIGMA;
    let nv = (v + prod + diff - adv) / (1.0 + cw1 * fw * v / (d * d));
    out = clamp(nv, 0.0, NU_TILDE_MAX);
    if ((P.essai & 2048u) != 0u) { out = min(out, P.essaiNu); }
  }
  nutOut[idx] = out;
}
";

/// Transport de ν̃ (ou de k, ε) des grilles boîtes à grands plans : mailles parcourues
/// par briques de rangées x entières, SA_BY rangées en y et SA_BZ en z
/// (briques voisines en y consécutives), au lieu de plans entiers. Au-delà
/// de SA_PLAN mailles par plan, les trois plans voisins d'une maille ne
/// tiennent plus dans le cache (2 Mo) et ν̃, u y étaient relus trois fois ;
/// en deçà, l'ordre des plans reste plus rapide. Même calcul par maille que
/// `SA3_WGSL` (ou `KE3_WGSL`, même en-tête). Essais (Ahmed, plan de 53 k mailles) : 4 × 16 1,455 → 0,82 ms,
/// 8 × 8 0,85, 4 × 4 0,93, 8 × 2 1,10 ; colonnes z en registres (un fil par
/// (x, y) et 4 à 64 plans) 0,93 à 1,21, et plus lentes sur les petits plans.
const SA_BY: u32 = 4;
const SA_BZ: u32 = 16;
const SA_PLAN: usize = 32768;

fn bricked(base: String, bricks: bool) -> String {
    if !bricks {
        return base;
    }
    let head = "  let idx = gid.x + gid.y * nwg.x * 64u;
  if (idx >= P.N) { return; }
  let pp = pos_of(idx);
  let x = pp.x;
  let y = pp.y;
  let z = pp.z;
";
    assert!(base.contains(head), "transport : en-tête");
    let (by, bz) = (SA_BY, SA_BZ);
    base.replace(
        head,
        &format!(
            "  let t = gid.x + gid.y * nwg.x * 64u;
  let r = t / P.nx;
  let x = i32(t % P.nx);
  let b = r / {bb}u;
  let w = r % {bb}u;
  let nby = (P.ny + {by}u - 1u) / {by}u;
  let yu = b % nby * {by}u + w % {by}u;
  let zu = b / nby * {bz}u + w / {by}u;
  if (yu >= P.ny || zu >= P.nz) {{ return; }}
  let y = i32(yu);
  let z = i32(zu);
  let idx = u32(x) + P.nx * (yu + P.ny * zu);
",
            bb = by * bz
        ),
    )
}

/// k-ε RNG (Yakhot et al. 1992) : constantes, viscosité turbulente et état
/// amont, communs au transport, à la collision et aux raccords.
const KE_FNS: &str = "
const CMU = 0.0845;
// Plafond de ν_t (réseau) : stabilité de la diffusion explicite (6 voisins).
const NUT_KE_MAX = 0.08;
// Amont : intensité Tu = 0,1 % et ν_t/ν = 0,2 (celui de ν̃ = 3ν en SA).
const TU_INF = 0.001;
const RNUT_INF = 0.2;
// ν_t = C_μ k²/ε : seule définition de ν_t du modèle RANS.
fn nu_ke(k: f32, e: f32) -> f32 {
  if (k <= 0.0 || e <= 0.0) { return 0.0; }
  return min(NUT_KE_MAX, CMU * k * k / e);
}
// VLES : filtre de maille de Johansen, Wu et Shyy (2004),
// F = min(1, C3 Δ ε/k^{3/2}), C3 = 1, Δ = 1 (pas du réseau local) : là où
// l'échelle des tourbillons porteurs k^{3/2}/ε dépasse la maille, seule la
// part non résolue est modélisée (ν_t → C_μ Δ √k). Hors couche limite
// attachée seulement : f_d de la DDES (Spalart et al. 2006, voir `ke3d`),
// ν_t = ν_RANS (1 − f_d (1 − F)) ; f_d = 0 garde le RANS en paroi.
fn nu_vles(k: f32, e: f32, fd: f32) -> f32 {
  let nr = nu_ke(k, e);
  if (nr <= 0.0) { return 0.0; }
  let f = min(1.0, e / (k * sqrt(k)));
  return nr * (1.0 - fd * (1.0 - f));
}
fn k_inf() -> f32 {
  let t = TU_INF * P.uIn;
  return 1.5 * t * t;
}
fn e_inf() -> f32 {
  let k = k_inf();
  return CMU * k * k / (RNUT_INF * P.nu);
}
";

/// Transport explicite de k et ε (k-ε RNG), même schéma que ν̃ (`SA3_WGSL`) :
/// advection décentrée amont, diffusion conservative, production ν_t S²,
/// destructions semi-implicites, terme R de RNG dans C_ε2*. Termes
/// d'entretien amont (Spalart et Rumsey 2007) : k et ε n'y décroissent pas
/// entre l'entrée et la maquette. Paroi : k = u_τ²/√C_μ, ε = u_τ³/(κ h),
/// ν_t = κ u_τ h comme la couche pariétale de la loi de paroi. Champs rangés
/// par paires (k, ε) dans les tampons de ν̃.
const KE3_WGSL: &str = "
@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read> macroBuf: array<vec4<f32>>;
@group(0) @binding(2) var<storage, read> wall: array<u32>;
@group(0) @binding(3) var<storage, read> wallGeo: array<vec4<f32>>;
@group(0) @binding(4) var<storage, read> utau: array<f32>;
@group(0) @binding(5) var<storage, read> dist: array<f32>;
@group(0) @binding(6) var<storage, read> keIn: array<f32>;
@group(0) @binding(7) var<storage, read_write> keOut: array<f32>;
@group(0) @binding(8) var<storage, read> slab: array<vec2<f32>>;
@group(0) @binding(9) var<storage, read> blk: array<u32>;
@group(0) @binding(10) var<storage, read> plain: array<u32>;
SPARSE_FNS
KE_FNS
fn lam_xl(x: i32, y: i32) -> vec2<f32> {
  let e = slab[u32(y)];
  return vec2<f32>(f32(x) - e.x, e.y - e.x);
}

const KAPPA = 0.41;
const WM_H = 2.0;
// Yakhot et al. (1992) : σ_k = σ_ε = 1/1,393, C_ε1, C_ε2, η₀, β.
const SIGK = 0.7194;
const SIGE = 0.7194;
const CE1 = 1.42;
const CE2 = 1.68;
const ETA0 = 4.38;
const BETA = 0.012;

fn ke(i: u32) -> vec2<f32> {
  return vec2<f32>(keIn[2u * i], keIn[2u * i + 1u]);
}
// VLES (remplacé à la compilation) : ν_t réduit rangé dans le dernier tiers.
const VLES = false;
// ν_t réduit aussi dans le transport de k et ε (production, diffusion) :
// essai écarté (README, « VLES »). Par défaut k et ε suivent le k-ε RNG
// complet ; seule la collision prend le ν_t filtré (le τ_turb de Chen et
// al. 2003 n'entre que dans l'équation de Boltzmann).
const VLES_TR = false;
// ν_t d'une maille au pas précédent (diffusion).
fn nt_of(i: u32) -> f32 {
  if (VLES && VLES_TR) { return keIn[arrayLength(&keIn) / 3u * 2u + i]; }
  let v = ke(i);
  return nu_ke(v.x, v.y);
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let idx = gid.x + gid.y * nwg.x * 64u;
  if (idx >= P.N) { return; }
  let pp = pos_of(idx);
  let x = pp.x;
  let y = pp.y;
  let z = pp.z;
  var wk0 = 0u;
  if (((plain[idx >> 5u] >> (idx & 31u)) & 1u) == 0u) { wk0 = wall[idx]; }
  if (wk0 == 0xfffffffdu) { return; }
  let nu = P.nu;
  let kin = k_inf();
  let ein = e_inf();
  var out: vec2<f32>;
  // VLES : ν_t réduit de la maille (< 0 : celui de k-ε, hors cœur).
  var vt_out = -1.0;
  if (wk0 == 0xffffffffu) {
    out = vec2<f32>(0.0);
  } else if (wk0 == 0xfffffffeu) {
    out = ke(idx);
  } else if (x == 0) {
    out = vec2<f32>(kin, ein);
  } else if (x == i32(P.nx) - 1) {
    out = ke(idx - 1u);
  } else if (wk0 != 0u && wallGeo[wk0 - 1u].w != 0.0) {
    // Premier nœud : équilibre de la zone logarithmique à la hauteur de
    // raccord h de la loi de paroi (ν_t = C_μ k²/ε = κ u_τ h).
    let ut = utau[4u * (wk0 - 1u)];
    out = vec2<f32>(ut * ut / sqrt(CMU), ut * ut * ut / (KAPPA * WM_H));
    out = max(out, vec2<f32>(kin, ein));
  } else {
    let xc = u32(x);
    let yc = u32(y);
    let zc = u32(z);
    let ground = (P.rows & 2u) != 0u;
    let ym = select(yc - 1u, P.ny - 1u, y == 0);
    let yp = select(yc + 1u, 0u, yc + 1u == P.ny);
    let zm = select(zc - 1u, select(P.nz - 1u, 0u, ground), z == 0);
    let zp = select(zc + 1u, select(0u, zc, ground), zc + 1u == P.nz);
    let ixm = id_of(xc - 1u, yc, zc);
    let ixp = id_of(xc + 1u, yc, zc);
    let iym = id_of(xc, ym, zc);
    let iyp = id_of(xc, yp, zc);
    let izm = id_of(xc, yc, zm);
    let izp = id_of(xc, yc, zp);
    let v = ke(idx);
    let vxm = ke(ixm);
    let vxp = ke(ixp);
    let vym = ke(iym);
    let vyp = ke(iyp);
    let vzm = ke(izm);
    let vzp = ke(izp);
    let u = macroBuf[idx].xyz;
    var adv = select(u.x * (vxp - v), u.x * (v - vxm), u.x > 0.0);
    adv += select(u.y * (vyp - v), u.y * (v - vym), u.y > 0.0);
    adv += select(u.z * (vzp - v), u.z * (v - vzm), u.z > 0.0);
    // Diffusion : (ν + ν_t/σ) aux faces, ν_t moyen des deux mailles.
    let n0 = nt_of(idx);
    let s = vec2<f32>(1.0 / SIGK, 1.0 / SIGE);
    let diff = (nu + 0.5 * (n0 + nt_of(ixp)) * s) * (vxp - v) - (nu + 0.5 * (n0 + nt_of(ixm)) * s) * (v - vxm)
      + (nu + 0.5 * (n0 + nt_of(iyp)) * s) * (vyp - v) - (nu + 0.5 * (n0 + nt_of(iym)) * s) * (v - vym)
      + (nu + 0.5 * (n0 + nt_of(izp)) * s) * (vzp - v) - (nu + 0.5 * (n0 + nt_of(izm)) * s) * (v - vzm);
    // Taux de déformation S = √(2 S_ij S_ij), différences centrées.
    let gx = 0.5 * (macroBuf[ixp].xyz - macroBuf[ixm].xyz);
    let gy = 0.5 * (macroBuf[iyp].xyz - macroBuf[iym].xyz);
    let gz = 0.5 * (macroBuf[izp].xyz - macroBuf[izm].xyz);
    let sxy = gx.y + gy.x;
    let sxz = gx.z + gz.x;
    let syz = gy.z + gz.y;
    let s2 = 2.0 * (gx.x * gx.x + gy.y * gy.y + gz.z * gz.z) + sxy * sxy + sxz * sxz + syz * syz;
    var nt = nu_ke(v.x, v.y);
    if (VLES) {
      // DDES : r_d = (ν_t + ν)/(κ² d² √(U_ij U_ij)) ≈ 1 dans la zone
      // logarithmique, → 0 loin des parois ; f_d = 1 − tanh((8 r_d)³).
      let gn = sqrt(dot(gx, gx) + dot(gy, gy) + dot(gz, gz));
      let dw = max(0.5, dist[idx]);
      let rd = (nt + nu) / max(KAPPA * KAPPA * dw * dw * gn, 1e-20);
      let a = 8.0 * rd;
      let fd = 1.0 - tanh(min(a * a * a, 20.0));
      vt_out = nu_vles(v.x, v.y, fd);
      if (VLES_TR) { nt = vt_out; }
    }
    var prod = nt * s2;
    if (P.xtr > 0.0 && dist[idx] < P.band) {
      let xl = lam_xl(x, y);
      if (xl.y > 0.0 && xl.x < P.xtr * xl.y) { prod = 0.0; }
    }
    let kk = max(v.x, 1e-20);
    let ee = max(v.y, 1e-30);
    // RNG : C_ε2* = C_ε2 + C_μ η³ (1 − η/η₀)/(1 + β η³), η = S k/ε.
    let eta = sqrt(s2) * kk / ee;
    let e3 = eta * eta * eta;
    let ce2s = CE2 + CMU * e3 * (1.0 - eta / ETA0) / (1.0 + BETA * e3);
    let tinv = ee / kk;
    // Entretien amont (Spalart et Rumsey 2007) : ε∞ dans l'équation de k,
    // C_ε2 ε∞²/k∞ dans celle de ε.
    var sk = 0.0;
    var se = 0.0;
    if (kin > 0.0) {
      sk = ein;
      se = CE2 * ein * ein / kin;
    }
    let kn = (v.x + prod + sk + diff.x - adv.x) / (1.0 + tinv);
    // Destruction implicite pour C_ε2* > 0, explicite (source) sinon.
    let en = (v.y + CE1 * tinv * prod + se + max(-ce2s, 0.0) * ee * tinv + diff.y - adv.y) / (1.0 + max(ce2s, 0.0) * tinv);
    out = max(vec2<f32>(kn, en), 1e-4 * vec2<f32>(kin, ein));
  }
  keOut[2u * idx] = out.x;
  keOut[2u * idx + 1u] = out.y;
  if (VLES) {
    if (vt_out < 0.0) { vt_out = nu_ke(out.x, out.y); }
    keOut[arrayLength(&keOut) / 3u * 2u + idx] = vt_out;
  }
}
";

/// Loi de paroi (Spalart–Allmaras) : une passe sur les seules mailles
/// pariétales, avant chaque pas de collision. Tout ce qu'elle lit vient du
/// pas précédent (vitesse des nœuds, u_τ et glissement relaxé) : la
/// collision lisait autrefois la vitesse des nœuds voisins pendant que la
/// même passe la réécrivait (relevé à WM_H), tantôt au pas précédent,
/// tantôt au pas courant selon l'ordre d'exécution des groupes ; résultats
/// non reproductibles d'un calcul à l'autre (≈ 10⁻³ sur les efforts
/// instantanés en régime transitoire). Elle écrit u_τ, le glissement relaxé
/// et la viscosité de la couche pariétale (< 0 : pas de loi de paroi), que
/// la collision relit ; le Newton de Spalding n'y fait plus attendre les
/// autres fils des groupes de la collision.
const WALL3_WGSL: &str = "
@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read> macroBuf: array<vec4<f32>>;
@group(0) @binding(2) var<storage, read> wall: array<u32>;
@group(0) @binding(3) var<storage, read> wallGeo: array<vec4<f32>>;
@group(0) @binding(4) var<storage, read_write> utau: array<f32>;
@group(0) @binding(5) var<storage, read> wdist: array<f32>;
@group(0) @binding(6) var<storage, read> slab: array<vec2<f32>>;
@group(0) @binding(7) var<storage, read> blk: array<u32>;
// Rang de chaque maille pariétale (dans l'ordre de leur numéro).
@group(0) @binding(8) var<storage, read> wcell: array<u32>;
SPARSE_FNS
LAM_FNS
SA_FNS
fn cell(xi: i32, yi: i32, zi: i32) -> u32 {
  let ny = i32(P.ny);
  let nz = i32(P.nz);
  let yy = (yi + ny) % ny;
  let zz = (zi + nz) % nz;
  let xx = clamp(xi, 0, i32(P.nx) - 1);
  return id_of(u32(xx), u32(yy), u32(zz));
}

// (u, ρ) interpolés (trilinéaire) en p ; ρ = 0 si un sommet est solide ou
// hors veine.
fn mac_at(p: vec3<f32>) -> vec4<f32> {
  let i0 = vec3<i32>(floor(p));
  let fr = p - vec3<f32>(i0);
  var r = vec4<f32>(0.0);
  for (var c = 0u; c < 8u; c++) {
    let o = vec3<i32>(i32(c & 1u), i32((c >> 1u) & 1u), i32((c >> 2u) & 1u));
    let q3 = i0 + o;
    if (q3.x < 0 || q3.x >= i32(P.nx) || q3.z < 0 || q3.z >= i32(P.nz)) { return vec4<f32>(0.0); }
    let ci = cell(q3.x, q3.y, q3.z);
    if (wall[ci] >= 0xfffffffeu) { return vec4<f32>(0.0); }
    let wv = select(1.0 - fr, fr, vec3<bool>(o.x == 1, o.y == 1, o.z == 1));
    r += wv.x * wv.y * wv.z * macroBuf[ci];
  }
  return r;
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let wk = gid.x + gid.y * nwg.x * 64u;
  if (wk >= arrayLength(&wcell)) { return; }
  let idx = wcell[wk];
  let hasWall = true;
  let pp = pos_of(idx);
  let x = pp.x;
  let y = pp.y;
  let z = pp.z;
  // Loi de paroi : glissement tangentiel tiré de la vitesse du pas précédent.
  var uw = vec3<f32>(0.0);
  // Viscosité turbulente de la couche pariétale (κ u_τ h), < 0 hors paroi.
  var nuw = -1.0;
  if (P.sa != 0u && hasWall) {
    let gw0 = wallGeo[wk];
    // Paroi volumique (surfels) : maille réceptrice sans lien vers la paroi,
    // hors loi de paroi (maille ordinaire pour la collision et ν̃).
    if (gw0.w == 0.0) {
      utau[4u * wk] = 0.0;
      utau[arrayLength(&utau) / 6u * 5u + wk] = -1.0;
      return;
    }
    // Arête vive (distance négative) : pas de loi de paroi, paroi adhérente.
    // Diagnostic : loi de paroi aussi aux arêtes vives (262144).
    let edge = gw0.w < 0.0 && (P.essai & 262144u) == 0u;
    let gw = vec4<f32>(gw0.xyz, abs(gw0.w));
    let u0 = macroBuf[idx].xyz;
    // Vitesse tangentielle du nœud et sa distance à la paroi.
    let tv0 = u0 - dot(u0, gw.xyz) * gw.xyz;
    let ut0 = length(tv0);
    let y0 = max(0.1, gw.w);
    var ut = ut0;
    var yw = y0;
    var th = tv0;
    // Gradient de pression cinématique le long du glissement (essai 8388608).
    var pg = 0.0;
    {
      // Loi de paroi alimentée à une hauteur fixe h = WM_H au-dessus du pied
      // de la normale (Kawai et Larsson 2012) : au premier nœud, l'erreur
      // numérique et de sous-maille décale le profil logarithmique, et la
      // distance du nœud à la paroi saute d'une maille quand la paroi
      // traverse une rangée de nœuds. Repli sur le nœud si un voisin de
      // l'interpolation est solide, ou si le point relevé n'est pas à la
      // hauteur h de la paroi la plus proche (bord mince, coin : le relevé
      // tomberait près d'une autre paroi).
      let ps = vec3<f32>(f32(x), f32(y), f32(z)) + (WM_H - gw.w) * gw.xyz;
      let i0 = vec3<i32>(floor(ps));
      let fr = ps - vec3<f32>(i0);
      var us3 = vec3<f32>(0.0);
      var ds = 0.0;
      // Masse volumique et son gradient (dérivée de l'interpolation
      // trilinéaire) au point relevé.
      var rh = 0.0;
      var gr = vec3<f32>(0.0);
      var ok = true;
      for (var c = 0u; c < 8u; c++) {
        let o = vec3<i32>(i32(c & 1u), i32((c >> 1u) & 1u), i32((c >> 2u) & 1u));
        let q3 = i0 + o;
        if (q3.x < 0 || q3.x >= i32(P.nx) || q3.z < 0 || q3.z >= i32(P.nz)) { ok = false; break; }
        let ci = cell(q3.x, q3.y, q3.z);
        if (wall[ci] >= 0xfffffffeu) { ok = false; break; }
        let ob = vec3<bool>(o.x == 1, o.y == 1, o.z == 1);
        let wv = select(1.0 - fr, fr, ob);
        let sg = select(vec3<f32>(-1.0), vec3<f32>(1.0), ob);
        let m4 = macroBuf[ci];
        us3 += wv.x * wv.y * wv.z * m4.xyz;
        ds += wv.x * wv.y * wv.z * wdist[ci];
        rh += wv.x * wv.y * wv.z * m4.w;
        gr += m4.w * vec3<f32>(sg.x * wv.y * wv.z, wv.x * sg.y * wv.z, wv.x * wv.y * sg.z);
      }
      if (ok && abs(ds - WM_H) < 0.5) {
        th = us3 - dot(us3, gw.xyz) * gw.xyz;
        ut = length(th);
        yw = WM_H;
        // P = (1/ρ) dp/ds le long de la vitesse tangentielle relevée,
        // p = ρ c_s² = ρ/3.
        if ((P.essai & 8388608u) != 0u && ut > 1e-12) {
          pg = dot(gr, th) / (3.0 * rh * ut);
          // Essai 536870912 : gradient filtré, différence centrée de ρ
          // interpolé à h, à ±D = 2 mailles le long de u_h (moyenne du
          // gradient sur 4 mailles le long de la paroi) : les fluctuations
          // de pression du réseau, d'une à deux mailles, s'y annulent. Point
          // non valide (solide, hors veine, autre paroi) : pas de gradient.
          if ((P.essai & 536870912u) != 0u) {
            let e = th / ut;
            let mp = mac_at(ps + 2.0 * e);
            let mm = mac_at(ps - 2.0 * e);
            let ok2 = mp.w > 0.0 && mm.w > 0.0;
            pg = select(0.0, (mp.w - mm.w) / (12.0 * rh), ok2);
            // Essai 1073741824 (avec 536870912) : gradient de la pression
            // totale à h, d(p/ρ + |u|²/2)/ds. Le long d'une ligne de courant
            // il vaut (1/ρ) ∂τ/∂y à h, la divergence de contrainte que le
            // fluide y subit réellement ; près de la paroi elle tend vers
            // dp/ds. Nul quand l'écoulement à h est non visqueux (couche
            // limite plus mince que h : loi d'équilibre retrouvée), il ne
            // corrige que lorsque la couche limite épaissie atteint h.
            if ((P.essai & 1073741824u) != 0u) {
              pg = select(0.0, pg + (dot(mp.xyz, mp.xyz) - dot(mm.xyz, mm.xyz)) / 8.0, ok2);
            }
          }
        }
      }
    }
    // Diagnostic 33554432 (avec 8388608) : u_τ de Spalding, le gradient ne
    // sert qu'au profil de la couche pariétale (compatibilité à la paroi,
    // ∂τ/∂y = dp/ds, exacte en y = 0).
    var ut2 = mellor_utau(ut, yw, select(pg, 0.0, (P.essai & 33554432u) != 0u));
    if (P.xtr > 0.0) {
      // Zone laminaire : Cf de Blasius, 0,664/√Re_x, plafonné au point d'arrêt.
      let xl = lam_xl(x, y);
      if (xl.y > 0.0 && xl.x < P.xtr * xl.y) {
        let rex = max(1.0, P.uIn * max(xl.x, 0.5) / P.nu);
        let cf = min(0.02, 0.664 / sqrt(rex));
        ut2 = ut * sqrt(0.5 * cf);
        pg = 0.0;
      }
    }
    utau[4u * wk] = ut2;
    // Couche pariétale [0, h] de viscosité uniforme ν_t = κ u_τ h (celle de
    // la loi de paroi à la hauteur de raccord) : profil linéaire de pente
    // τ_p / (ν + ν_t), et la paroi glisse à la vitesse que ce profil donne en
    // y = 0 depuis le nœud, u_n − pente · y_n. La vitesse extrapolée à la
    // paroi ne dépend plus de la place du nœud dans la maille (avec
    // u_n − u_τ²y_n/(ν + κu_τy_n), elle sautait d'environ u_τ ln(y)/κ d'un
    // nœud à l'autre quand la paroi traverse une rangée : marches de
    // pression, volets Gurney parasites, portance qui croît en affinant).
    nuw = KAPPA * ut2 * yw;
    // Diagnostic : viscosité de la couche pariétale κ u_τ y_n (distance du nœud).
    if ((P.essai & 256u) != 0u) { nuw = KAPPA * ut2 * y0; }
    // Gradient de pression : viscosité de longueur de mélange de la loi à
    // la hauteur de raccord, κ h √τ(h), τ(h) = u_τ² + P h (borné à 0).
    // Diagnostic : κ u_τ h gardé (16777216).
    if (pg != 0.0 && (P.essai & 16777216u) == 0u) { nuw = KAPPA * yw * sqrt(max(ut2 * ut2 + pg * yw, 0.0)); }
    if ((P.essai & 1024u) != 0u) { nuw = min(nuw, P.essaiNu); }
    // Couche pariétale : contrainte u_τ² + P y (bornée à 0 comme dans la
    // loi), profil u(y) = u_s + ∫₀^y τ / (ν + ν_t) : le glissement ne dépasse
    // jamais la vitesse du nœud (sinon la paroi entraînerait le fluide).
    let dpn = tau_int(ut2 * ut2, pg, y0) / (P.nu + nuw);
    let dph = tau_int(ut2 * ut2, pg, yw) / (P.nu + nuw);
    var usv = vec3<f32>(0.0);
    if (edge) {
      // Arête vive : l'écoulement y décolle et la loi de paroi (couche
      // limite attachée le long d'une face) ne s'applique pas ; sans elle,
      // la paroi presque glissante laisse croître la vitesse de contournement
      // de l'arête (singulière en non visqueux) jusqu'à la divergence.
      nuw = -1.0;
    } else if (ut0 >= 1e-12) {
      // Profil monotone entre le nœud et la hauteur h : le profil de la
      // couche pariétale, évalué à la hauteur du nœud, ne dépasse pas la
      // vitesse relevée à h (u_s ≤ u_h − pente · y_n). La borne arrête un
      // nœud qui s'emballerait (un nœud presque à la paroi glisse presque
      // librement et rien d'autre ne le retient). Elle ne suppose pas que la
      // couche pariétale linéaire atteint h : quand la couche limite est plus
      // mince que h (bord d'attaque, et partout aux résolutions courantes),
      // ou sur une paroi convexe où la vitesse potentielle décroît en
      // s'éloignant de la paroi, la borne u_h − pente · h freinait le premier
      // nœud sans frottement correspondant et épaississait la couche lente
      // (plaque plane, `plaque3d` : θ −12 %, perte de pression totale du nez
      // −30 %).
      var us = max(0.0, min(ut0, ut) - ut2 * ut2 * y0 / (P.nu + nuw));
      if (pg != 0.0) { us = max(0.0, min(ut0, ut) - dpn); }
      // Plancher : le profil de la couche pariétale ancré à h. Là où la paroi
      // croise une rangée de nœuds, le nouveau premier nœud (y_n ≪ 1, dans
      // le coin d'une marche : liens vers la paroi dessous et vers l'amont)
      // a la vitesse du glissement par l'interpolation de Bouzidi ; le
      // glissement tiré de sa propre vitesse n'est alors plus fixé que par
      // lui-même et dérive vers le bas (u_n ≈ 0,5 U contre 0,9 U aux voisins,
      // NACA 0012 à 4°, 128 mailles par corde) : la couche lente d'une maille
      // qui en résulte est advectée jusqu'au bord de fuite (perte de pression
      // totale du premier rang ×3 à ×5). Le profil relevé à h n'en dépend
      // pas ; là où le nœud donne plus (paroi convexe, couche plus mince que
      // h, gradient adverse), le glissement reste celui du nœud. Diagnostic :
      // sans plancher (4194304).
      if ((P.essai & 4194304u) == 0u) { us = max(us, ut - select(ut2 * ut2 * yw / (P.nu + nuw), dph, pg != 0.0)); }
      // Diagnostics : sans borne (4096) ; ancienne borne u_h − pente · h (8192).
      if ((P.essai & 4096u) != 0u) { us = max(0.0, ut0 - ut2 * ut2 * y0 / (P.nu + nuw)); }
      // Diagnostic : profil ancré à h seulement (131072).
      if ((P.essai & 131072u) != 0u) { us = max(0.0, ut - ut2 * ut2 * yw / (P.nu + nuw)); }
      if ((P.essai & 8192u) != 0u) { us = max(0.0, min(ut0 - ut2 * ut2 * y0 / (P.nu + nuw), ut - ut2 * ut2 * yw / (P.nu + nuw))); }
      usv = tv0 * (us / ut0);
      // Essai 67108864 : glissement signé, porté par la vitesse relevée à h.
      // Le frottement de la loi s'oppose à u_h ; le glissement va dans le sens
      // de u_h et la vitesse du nœud n'y entre que par sa composante le long
      // de u_h (négative si le nœud est à contre-courant) ; plancher : profil
      // ancré à h. Variante 134217728 : sans plancher, le glissement suit le
      // nœud et peut s'inverser avant u_h.
      if ((P.essai & 67108864u) != 0u && ut > 1e-12) {
        let e = th / ut;
        let a0 = dot(tv0, e);
        let g = ut2 * ut2 / (P.nu + nuw);
        var s = max(min(a0, ut) - g * y0, ut - g * yw);
        if ((P.essai & 134217728u) != 0u) { s = min(a0, ut) - g * y0; }
        usv = e * s;
      }
    }
    if ((P.essai & 16u) != 0u || ((P.essai & 512u) != 0u && abs(gw.z) < 0.999)) {
      // Diagnostic : glissement libre (paroi non visqueuse).
      usv = tv0;
      nuw = 0.0;
      // Diagnostic : glissement libre à la vitesse relevée à h (65536).
      if ((P.essai & 65536u) != 0u) { usv = th; }
    }
    // Contrainte imposée comme force au nœud (collision) : la paroi glisse
    // librement ; viscosité du nœud κ u_τ h gardée, ou nulle (32768).
    if ((P.essai & 16384u) != 0u && !edge) {
      usv = tv0;
      if ((P.essai & 32768u) != 0u) { nuw = 0.0; }
    }
    // Relaxation sur une dizaine de pas : une paroi qui suit à chaque pas la
    // vitesse de son propre nœud (glissement presque libre près de la paroi)
    // n'a pas de rappel et diverge ; en régime établi, sans effet (β = 0,1 ou
    // 0,02 : mêmes efforts à 0,3 % près).
    let prev = vec3<f32>(utau[4u * wk + 1u], utau[4u * wk + 2u], utau[4u * wk + 3u]);
    uw = prev + select(0.1, 1.0, (P.essai & 8u) != 0u) * (usv - prev);
    utau[4u * wk + 1u] = uw.x;
    utau[4u * wk + 2u] = uw.y;
    utau[4u * wk + 3u] = uw.z;
  }
  utau[arrayLength(&utau) / 6u * 5u + wk] = nuw;
}
";

/// Sol (bit 2 de `rows`) : population venant de sous le plancher (z = −½),
/// renvoyée par une paroi qui défile à la vitesse amont ; venant d'au-dessus
/// du plafond, réfléchie (paroi glissante).
fn ground_code(i: usize) -> String {
    let [cx, cy, cz] = C3[i];
    if cz == 0 {
        return String::new();
    }
    if cz > 0 {
        let j = opp3(i);
        format!(
            " else if ((P.rows & 2u) != 0u && z - 1 < 0) {{
      v = pop({j}u, idx) + {} * P.uIn;
    }}",
            f(6.0 * w3(i) * cx as f64)
        )
    } else {
        let m = C3.iter().position(|d| d[0] == cx && d[1] == cy && d[2] == -cz).unwrap();
        format!(
            " else if ((P.rows & 2u) != 0u && z + 1 >= i32(P.nz)) {{
      v = pop({m}u, idx);
    }}"
        )
    }
}

/// Collision. `walls` : passe des seules mailles pariétales (liste `wcell`),
/// avec les liens de paroi ; sinon passe des autres mailles, sans le code
/// des liens, recopié dans chacune des 18 directions, qui demandait des
/// registres à toutes les mailles (−24 % de durée sans lui).
/// Taille des groupes de la collision (cœur).
pub(crate) fn step_wg() -> u32 {
    std::env::var("SOUFFLERIE_WG").ok().and_then(|v| v.parse().ok()).unwrap_or(64)
}

fn wgsl_step(fp16: bool, p16: bool, walls: bool, surf: bool) -> String {
    // Paroi volumique (surfels) : la propagation des mailles réceptrices
    // passe par les surfels (`surfel3d`), sans liens de Bouzidi.
    let bouzidi = walls && !surf;
    let words = if fp16 { 10 } else { 19 };
    let half = if fp16 { 5 } else { 10 };
    let load = if p16 {
        "fn pop(i: u32, c: u32) -> f32 {
  if (i < 10u) { return unpack2x16float(u32(srcA[i * P.N + c])).x; }
  return unpack2x16float(u32(srcB[(i - 10u) * P.N + c])).x;
}".to_string()
    } else if fp16 {
        format!(
            "{hs}fn pop(i: u32, c: u32) -> f32 {{
  let h = HS[i];
  let w = h >> 1u;
  var v: vec2<f32>;
  if (w < {half}u) {{ v = unpack2x16float(srcA[w * P.N + c]); }} else {{ v = unpack2x16float(srcB[(w - {half}u) * P.N + c]); }}
  return select(v.x, v.y, (h & 1u) == 1u);
}}",
            hs = hs_wgsl()
        )
    } else {
        format!(
            "fn pop(i: u32, c: u32) -> f32 {{
  if (i < {half}u) {{ return bitcast<f32>(srcA[i * P.N + c]); }}
  return bitcast<f32>(srcB[(i - {half}u) * P.N + c]);
}}"
        )
    };
    let store = if p16 {
        "fn putH(i: u32, c: u32, a: f32) {
  let v = u16(pack2x16float(vec2<f32>(a, 0.0)) & 0xffffu);
  if (i < 10u) { dstA[i * P.N + c] = v; } else { dstB[(i - 10u) * P.N + c] = v; }
}".to_string()
    } else if fp16 {
        format!(
            "fn putWord(w: u32, c: u32, a: f32, b: f32) {{
  let v = pack2x16float(vec2<f32>(a, b));
  if (w < {half}u) {{ dstA[w * P.N + c] = v; }} else {{ dstB[(w - {half}u) * P.N + c] = v; }}
}}"
        )
    } else {
        format!(
            "fn putPop(i: u32, c: u32, a: f32) {{
  if (i < {half}u) {{ dstA[i * P.N + c] = bitcast<u32>(a); }} else {{ dstB[(i - {half}u) * P.N + c] = bitcast<u32>(a); }}
}}"
        )
    };
    // Propagation : pour chaque direction i, population venant de x − c_i.
    let mut pull = String::new();
    for i in 1..19 {
        let [cx, cy, cz] = C3[i];
        let j = opp3(i);
        let wallpart = if bouzidi {
            format!(
                "    var q = 0.0;
    q = qa[{im1}];
    if (q > 0.0) {{
      let fj = pop({j}u, idx);
      let qo = qa[{jm1}];
      // Paroi mobile (glissement de la loi de paroi) : 6 w_i ρ (c_i · u_w),
      // u_w ramenée dans le plan de la face que touche le lien (aux arêtes,
      // la normale de la maille n'est celle d'aucune face : sa composante
      // normale ferait souffler la paroi).
      var uwl = uw;
      if (P.sa != 0u) {{
        let nf = oct_decode(ca[{im1}]);
        // Transpiration : vitesse normale de la face (déplacement de la
        // couche limite) ; la paroi mobile de Bouzidi injecte ce débit.
        uwl = uw - dot(uw, nf) * nf + vn * nf;
      }}
      let corr = {w6} * dot(vec3<f32>({ccx}, {ccy}, {ccz}), uwl);
      if (q < 0.5 && qo == 0.0) {{
        let nb = cell(x + {cx}, y + {cy}, z + {cz});
        v = 2.0 * q * fj + (1.0 - 2.0 * q) * pop({j}u, nb) + corr;
      }} else if (q < 0.5) {{
        v = fj + corr;
      }} else {{
        v = (fj + (2.0 * q - 1.0) * pop({i}u, idx) + corr) / (2.0 * q);
      }}
      let mom = fj + v + {w2};
      fsum += vec3<f32>({mcx}, {mcy}, {mcz}) * mom - uwl * (fj - v);
",
                im1 = i - 1,
                jm1 = j - 1,
                i = i,
                j = j,
                cx = cx,
                cy = cy,
                cz = cz,
                w2 = f(2.0 * w3(i)),
                w6 = f(6.0 * w3(i)),
                ccx = f(cx as f64),
                ccy = f(cy as f64),
                ccz = f(cz as f64),
                mcx = f(-cx as f64),
                mcy = f(-cy as f64),
                mcz = f(-cz as f64),
            )
        } else {
            "    if (false) {\n".to_string()
        };
        // Surfels : part transmise de la maille amont (poids T/D) et parts
        // émises par les surfels (Γᵒᵘᵗ/P rangées après leurs efforts dans
        // `force`, poids w/D) ; poids −1 : aucune source, rebond simple.
        let post = if walls && surf {
            format!(
                "    {{
      let a = qa[{im1}];
      if (a < 0.0) {{
        v = pop({j}u, idx);
      }} else {{
        var s = 0.0;
        let rb = 18u * WALLN_C;
        let re = 36u * WALLN_C + 1u;
        let e1 = bitcast<u32>(qbuf[rb + wk * 18u + {i}u]);
        for (var e = bitcast<u32>(qbuf[rb + wk * 18u + {im1}u]); e < e1; e++) {{
          let k = bitcast<u32>(qbuf[re + 2u * e]);
          s += qbuf[re + 2u * e + 1u] * force[k >> 2u][k & 3u];
        }}
        v = select(0.0, a * v, a > 0.0) + s;
      }}
    }}
",
                im1 = i - 1,
                i = i,
                j = opp3(i),
            )
        } else {
            String::new()
        };
        let _ = write!(
            pull,
            "  {{
    var v: f32;
{wallpart}    }}{ground}{perx} else if (x - {cx} < 0) {{
      v = feqIn({i}u);
    }} else if (x - {cx} >= i32(P.nx)) {{
      v = pop({i}u, idx);
    }} else {{
      v = pop({i}u, nbr(u32(x - {cx}), {yv}, {zv}, r_{yv}_{zv}));
    }}
{post}    g[{i}] = v;
  }}
",
            yv = ["ym", "yc", "yp"][(1 - cy) as usize],
            zv = ["zm", "zc", "zp"][(1 - cz) as usize],
            ground = ground_code(i),
            // Veine périodique en x (bit 8 de `rows`, cas de validation sans
            // entrée ni sortie, voir `set_periodic_x`).
            perx = if cx == 0 {
                String::new()
            } else {
                format!(
                    " else if ((P.rows & 8u) != 0u && (x - {cx} < 0 || x - {cx} >= i32(P.nx))) {{
      v = pop({i}u, nbr(u32((x - {cx} + i32(P.nx)) % i32(P.nx)), {yv}, {zv}, r_{yv}_{zv}));
    }}",
                    yv = ["ym", "yc", "yp"][(1 - cy) as usize],
                    zv = ["zm", "zc", "zp"][(1 - cz) as usize],
                )
            },
        );
    }
    let dotc = |i: usize, a: usize, b: usize| C3[i][a] * C3[i][b];
    // Collision par paires de directions opposées (i, −i) : partie paire
    // (équilibre en ρ, u², et partie hors équilibre d'ordre 2) et impaire
    // (terme en c·u, ordre 3 de Hermite, force) calculées une fois pour les
    // deux, f_i = paire + impaire, f_−i = paire − impaire. Même formule ;
    // l'arrondi change (≈ 10⁻⁷ par pas), collision 10 % plus rapide.
    let mut coll = String::new();
    let mut absorb_code = String::new();
    for i in 0..19 {
        let j = opp3(i);
        let [cx, cy, cz] = C3[i];
        let (cxf, cyf, czf) = (cx as f64, cy as f64, cz as f64);
        let w = w3(i);
        let _ = write!(
            absorb_code,
            "    {{
      let ci = {cxs} * P.uIn;
      o[{i}] -= sig * (o[{i}] - {wf} * (3.0 * ci + 4.5 * ci * ci - 1.5 * P.uIn * P.uIn));
    }}
",
            wf = f(w),
            cxs = f(cxf),
        );
        if j < i {
            continue;
        }
        let comps: Vec<String> = [(cx, "u.x"), (cy, "u.y"), (cz, "u.z")]
            .iter()
            .filter(|(c, _)| *c != 0)
            .map(|(c, s)| if *c == 1 { s.to_string() } else { format!("-{s}") })
            .collect();
        let cu = if comps.is_empty() {
            "0.0".to_string()
        } else {
            comps.join(" + ").replace("+ -", "- ")
        };
        // Termes de coefficient nul omis (0 · x n'est pas simplifiable en IEEE).
        let qp = [
            (cxf * cxf - 1.0 / 3.0, "pxx"),
            (cyf * cyf - 1.0 / 3.0, "pyy"),
            (czf * czf - 1.0 / 3.0, "pzz"),
            (2.0 * cxf * cyf, "pxy"),
            (2.0 * cxf * czf, "pxz"),
            (2.0 * cyf * czf, "pyz"),
        ]
        .iter()
        .filter(|(c, _)| *c != 0.0)
        .map(|(c, v)| format!("{} * {v}", f(*c)))
        .collect::<Vec<_>>()
        .join(" + ");
        // Correction du défaut cubique (essai, voir `MACH3_BIT`) : contribution
        // des termes correctifs s_αα au moment d'ordre 2 après collision.
        let qs = [(cxf * cxf - 1.0 / 3.0, "sxx"), (cyf * cyf - 1.0 / 3.0, "syy"), (czf * czf - 1.0 / 3.0, "szz")]
            .iter()
            .map(|(c, v)| format!("{} * {v}", f(*c)))
            .collect::<Vec<_>>()
            .join(" + ");
        let qs = format!("{M3_OPEN} + {} * ({qs}){M3_CLOSE}", f(4.5 * w));
        if j == i {
            let _ = writeln!(coll, "  o[{i}] = {wf} * rho * (1.0 - usq) + k * {w45} * ({qp}){qs} - {wf};", wf = f(w), w45 = f(4.5 * w));
            continue;
        }
        let accel = match cx {
            0 => String::new(),
            1 => format!(" + {} * rho * (P.accel)", f(3.0 * w)),
            _ => format!(" + {} * rho * (-P.accel)", f(3.0 * w)),
        };
        // Hermite d'ordre 3 (composantes portées par D3Q19, sans xyz).
        let h = |a: f64, b: f64| a * (b * b - 1.0 / 3.0);
        let h3: Vec<String> = [
            (h(cyf, cxf), "axxy"),
            (h(cxf, cyf), "axyy"),
            (h(czf, cxf), "axxz"),
            (h(cxf, czf), "axzz"),
            (h(czf, cyf), "ayyz"),
            (h(cyf, czf), "ayzz"),
        ]
        .iter()
        .filter(|(c, _)| c.abs() > 1e-12)
        .map(|(c, v)| format!("{} * {v}", f(*c)))
        .collect();
        let h3 = if h3.is_empty() { "0.0".to_string() } else { h3.join(" + ") };
        let _ = write!(
            coll,
            "  {{
    let cu = {cu};
    let ev = {wf} * rho * (1.0 + 4.5 * cu * cu - usq) + k * {w45} * ({qp}){qs};
    let od = {w3f} * rho * cu + {w135} * ({h3}){accel};
    o[{i}] = ev + od - {wf};
    o[{j}] = ev - od - {wf};
  }}
",
            wf = f(w),
            w3f = f(3.0 * w),
            w45 = f(4.5 * w),
            w135 = f(13.5 * w),
        );
    }
    let out: String = if p16 {
        (0..19).map(|i| format!("  putH({i}u, idx, o[{i}]);
")).collect()
    } else if fp16 {
        let hd = hdir();
        let o = |h: usize| hd[h].map_or("0.0".to_string(), |i| format!("o[{i}]"));
        (0..words).map(|w| format!("  putWord({w}u, idx, {}, {});\n", o(2 * w), o(2 * w + 1))).collect()
    } else {
        (0..19).map(|i| format!("  putPop({i}u, idx, o[{i}]);\n")).collect()
    };
    let feq_in: String = (1..19)
        .map(|i| format!("if (i == {i}u) {{ cx = {}; w = {}; }}\n  ", f(C3[i][0] as f64), f(w3(i))))
        .collect();
    let fv: String = (0..19).map(|i| format!("fv[{i}] = g[{i}] + {}; ", f(w3(i)))).collect();
    let pt = if p16 { "u16" } else { "u32" };
    format!(
        "{en}{PARAMS_WGSL}
@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read> srcA: array<{pt}>;
@group(0) @binding(2) var<storage, read> srcB: array<{pt}>;
@group(0) @binding(3) var<storage, read_write> dstA: array<{pt}>;
@group(0) @binding(4) var<storage, read_write> dstB: array<{pt}>;
@group(0) @binding(5) var<storage, read> wall: array<u32>;
@group(0) @binding(6) var<storage, read> qbuf: array<f32>;
@group(0) @binding(7) var<storage, read_write> force: array<vec4<f32>>;
@group(0) @binding(8) var<storage, read_write> macroBuf: array<vec4<f32>>;
// Spalart–Allmaras : ν̃ par maille ; par maille pariétale, normale et distance
// du nœud à la paroi (x, y, z, d), et u_τ écrit pour la passe de transport.
@group(0) @binding(9) var<storage, read> nut: array<f32>;
@group(0) @binding(10) var<storage, read> wallGeo: array<vec4<f32>>;
@group(0) @binding(11) var<storage, read> utau: array<f32>;
// Distance à la paroi de chaque maille (amortissement de la sous-maille).
@group(0) @binding(12) var<storage, read> wdist: array<f32>;
// Étendue en x du corps par tranche y (zone laminaire).
@group(0) @binding(13) var<storage, read> slab: array<vec2<f32>>;
@group(0) @binding(14) var<storage, read> blk: array<u32>;
// Rang de chaque maille pariétale (passe des mailles pariétales), puis,
// sur une grille creuse, de chaque maille à calculer du cœur.
@group(0) @binding(15) var<storage, read> wcell: array<u32>;
// Mailles fluides ordinaires (un bit par maille ; voir `plain_mask`).
@group(0) @binding(16) var<storage, read> plain: array<u32>;
{SPARSE_FNS}
{lam_fns}
{sa_fns}
{ke_fns}
// Viscosité turbulente d'une maille : SA (ν̃) ou k-ε RNG (paires k, ε).
fn nut_cell(i: u32) -> f32 {{
  if (P.sa == 2u) {{ return nu_ke(nut[2u * i], nut[2u * i + 1u]); }}
  // VLES : ν_t réduit de la passe `ke3d`, dernier tiers du tampon.
  if (P.sa == 3u) {{ return nut[arrayLength(&nut) / 3u * 2u + i]; }}
  return nu_turb(nut[i]);
}}
{load}
{store}

fn cell(xi: i32, yi: i32, zi: i32) -> u32 {{
  let ny = i32(P.ny);
  let nz = i32(P.nz);
  let yy = (yi + ny) % ny;
  let zz = (zi + nz) % nz;
  let xx = clamp(xi, 0, i32(P.nx) - 1);
  return id_of(u32(xx), u32(yy), u32(zz));
}}

// Rang d'un voisin : grille boîte, x plus le début de sa rangée ; stockage
// par blocs, table des blocs.
fn nbr(xx: u32, yy: u32, zz: u32, row: u32) -> u32 {{
  if (P.bnx == 0u) {{ return xx + row; }}
  return id_of(xx, yy, zz);
}}

fn oct_decode(k: f32) -> vec3<f32> {{
  let iu = floor(k / 2048.0);
  let iv = k - iu * 2048.0;
  var u = iu / 2047.0 * 2.0 - 1.0;
  var v = iv / 2047.0 * 2.0 - 1.0;
  let z = 1.0 - abs(u) - abs(v);
  if (z < 0.0) {{
    let a = u;
    u = (1.0 - abs(v)) * sign(a);
    v = (1.0 - abs(a)) * sign(v);
  }}
  return normalize(vec3<f32>(u, v, z));
}}

fn absorb(x: i32, y: i32, z: i32) -> f32 {{
  // Aux Reynolds de maille modérés, la viscosité suffit à amortir les ondes ;
  // les bandes y ajouteraient un confinement (traînée de la sphère +1 % à Re = 20).
  if (P.les == 0u || (P.rows & 1u) != 0u) {{ return 0.0; }}
  let wy = max(2.0, f32(P.ny) / 8.0);
  let wz = max(2.0, f32(P.nz) / 8.0);
  let wx = max(2.0, f32(P.nx) / 24.0);
  let dy = f32(min(y, i32(P.ny) - 1 - y));
  let dz = f32(min(z, i32(P.nz) - 1 - z));
  var sz = 1.0 - dz / wz;
  if ((P.rows & 2u) != 0u) {{ sz = 1.0 - f32(i32(P.nz) - 1 - z) / wz; }}
  var sy = 1.0 - dy / wy;
  // Envergure périodique : la maquette traverse les bords y.
  if ((P.rows & 4u) != 0u) {{ sy = 0.0; }}
  var s = max(max(sy, sz), 1.0 - f32(x) / wx);
  if (u32(x) > P.spongeStart) {{
    s = max(s, f32(u32(x) - P.spongeStart) / f32(max(1u, P.nx - 1u - P.spongeStart)));
  }}
  s = clamp(s, 0.0, 1.0);
  return 0.05 * s * s;
}}

fn feqIn(i: u32) -> f32 {{
  var cx = 0.0;
  var w = 1.0 / 36.0;
  {feq_in}
  let cu = cx * P.uIn;
  return w * (3.0 * cu + 4.5 * cu * cu - 1.5 * P.uIn * P.uIn);
}}

@compute @workgroup_size({wg})
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {{
{head}
  // Solide, maille fantôme d'une grille fine (remplie par refine3d), ou
  // maille recouverte (recopiée de la grille fine, voir `mark_covered`).
  if (wk0 >= 0xfffffffdu) {{ return; }}
{skip}
  let wk = select(wk0 - 1u, 0xffffffffu, wk0 == 0u);
  let hasWall = wk != 0xffffffffu;
  let pp = pos_of(idx);
  let x = pp.x;
  let y = pp.y;
  let z = pp.z;
  // Rangées voisines, bords y et z périodiques (comme `cell`, sans modulo).
  let yc = u32(y);
  let zc = u32(z);
  let ym = select(yc - 1u, P.ny - 1u, y == 0);
  let yp = select(yc + 1u, 0u, yc + 1u == P.ny);
  let zm = select(zc - 1u, P.nz - 1u, z == 0);
  let zp = select(zc + 1u, 0u, zc + 1u == P.nz);
  // Grille boîte : début de chaque rangée voisine (rang = x + début).
  let sxy = P.nx * P.ny;
  let oym = ym * P.nx;
  let oyc = yc * P.nx;
  let oyp = yp * P.nx;
  let ozm = zm * sxy;
  let ozc = zc * sxy;
  let ozp = zp * sxy;
  let r_ym_zm = oym + ozm;
  let r_ym_zc = oym + ozc;
  let r_ym_zp = oym + ozp;
  let r_yc_zm = oyc + ozm;
  let r_yc_zc = oyc + ozc;
  let r_yc_zp = oyc + ozp;
  let r_yp_zm = oyp + ozm;
  let r_yp_zc = oyp + ozc;
  let r_yp_zp = oyp + ozp;
  var g: array<f32, 19>;
  var fsum = vec3<f32>(0.0);
  // Loi de paroi (passe `wall3d`, juste avant ce pas, sur le champ du pas
  // précédent) : glissement de la paroi et viscosité de la couche pariétale.
  var uw = vec3<f32>(0.0);
  // Viscosité turbulente de la couche pariétale (κ u_τ h), < 0 hors paroi.
  var nuw = -1.0;
  if (P.sa != 0u && hasWall) {{
    uw = vec3<f32>(utau[4u * wk + 1u], utau[4u * wk + 2u], utau[4u * wk + 3u]);
    nuw = utau[arrayLength(&utau) / 6u * 5u + wk];
  }}
  // Décambrage visqueux : vitesse de transpiration de la maille pariétale,
  // rangée après u_τ et le glissement (4 + 1 valeurs par maille).
  var vn = 0.0;
  if (P.sa != 0u && hasWall) {{
    vn = utau[arrayLength(&utau) / 6u * 4u + wk];
  }}
{prefetch}  g[0] = pop(0u, idx);
{pull}
  var fv: array<f32, 19>;
  {fv}
  let rho = {rho};
  let ir = 1.0 / rho;
  let u = vec3<f32>({jx}, {jy}, {jz}) * ir;
  let usq = 1.5 * dot(u, u);
  let pxx = {pxx} - rho * (1.0 / 3.0 + u.x * u.x);
  let pyy = {pyy} - rho * (1.0 / 3.0 + u.y * u.y);
  let pzz = {pzz} - rho * (1.0 / 3.0 + u.z * u.z);
  let pxy = {pxy} - rho * u.x * u.y;
  let pxz = {pxz} - rho * u.x * u.z;
  let pyz = {pyz} - rho * u.y * u.z;
  var tau = P.tau0;
  if (P.sa != 0u) {{
    // Nœud pariétal : la viscosité de la couche [0, h] supposée par le
    // glissement, pour que la contrainte transmise soit τ_p.
    var nuc = select(nut_cell(idx), nuw, nuw >= 0.0);
    // Diagnostics : ν_t de Spalart–Allmaras hors premier nœud retiré (1),
    // premier nœud à la viscosité de Spalart–Allmaras (4).
    if ((P.essai & 1u) != 0u && nuw < 0.0) {{ nuc = 0.0; }}
    if ((P.essai & 4u) != 0u && nuw >= 0.0) {{ nuc = nut_cell(idx); }}
    tau += 3.0 * nuc;
  }}
  if (u32(x) > P.spongeStart) {{
    let s = f32(u32(x) - P.spongeStart) / f32(max(1u, P.nx - 1u - P.spongeStart));
    tau += (P.tauSponge - P.tau0) * s * s;
  }}
  var om = 1.0 / tau;
  // Temps de relaxation sans Spalart–Allmaras (essai 2097152).
  var taub = P.tau0;
  if (u32(x) > P.spongeStart) {{
    let s = f32(u32(x) - P.spongeStart) / f32(max(1u, P.nx - 1u - P.spongeStart));
    taub += (P.tauSponge - P.tau0) * s * s;
  }}
  if (P.les == 1u && (P.essai & 32u) == 0u) {{
    let qn = sqrt(2.0 * (pxx * pxx + pyy * pyy + pzz * pzz + 2.0 * (pxy * pxy + pxz * pxz + pyz * pyz)));
    // Amortissement de Van Driest près de la paroi, (1 − e^(−d/A))², comme en
    // 2D, sans loi de paroi seulement : avec elle, la sous-maille stabilise les
    // premiers nœuds à τ ≈ 0,5 (aile à Re = 10⁶, L = 32 : divergence sinon).
    var damp = 1.0;
    if (P.sa == 0u || (P.essai & 2u) != 0u) {{ damp = 1.0 - exp(-wdist[idx] / P.lesA); }}
    var c18 = P.c18;
    // Diagnostic : constante de Smagorinsky divisée par 2 (64) ou par 4 (128).
    if ((P.essai & 64u) != 0u) {{ c18 *= 0.25; }}
    if ((P.essai & 128u) != 0u) {{ c18 *= 0.0625; }}
    om = 2.0 / (tau + sqrt(tau * tau + c18 * damp * damp * qn * ir));
    // Essai 2097152 : ν_t = max(ν_SA, ν_Smagorinsky) au lieu de leur somme.
    // Là où Spalart–Allmaras modélise toute la turbulence (couche limite,
    // sillage), la sous-maille la compterait deux fois. Avec Π_neq = −2ρτS/3,
    // la viscosité de Smagorinsky seule est celle de la formule de Hou au
    // temps de relaxation sans SA, et le maximum des deux est cohérent.
    if ((P.essai & 2097152u) != 0u) {{
      let taus = 0.5 * (taub + sqrt(taub * taub + c18 * damp * damp * qn * ir));
      om = 1.0 / max(tau, taus);
    }}
  }}
  let k = 1.0 - om;
{M3_OPEN}
  // Défaut cubique du D3Q19 (essai MACH3_BIT) : Σ f c_α³ = Σ f c_α, le
  // moment d'équilibre ρ u_α³ manque, et la contrainte normale vaut
  // 2ρν ∂_α u_α (1 − 4,5 u_α²) au lieu de 2ρν ∂_α u_α. Terme source
  // (Dellar 2014 ; Hosseini et al. 2020) ajouté au moment d'ordre 2 après
  // collision : s_αα = −(1 − ω/2) ∂_α(ρ u_α³) ≈ −(1 − ω/2) 3ρ u_α² ∂_α u_α,
  // le gradient tiré localement de la partie hors équilibre,
  // ∂_α u_α = −3ω Π¹_αα / (2ρ) (comme la correction galiléenne du
  // cumulant, Geier et al. 2015). Le défaut en ∂_γ(ρ u_x u_y u_z) des
  // contraintes de cisaillement (c_x c_y c_z = 0 sur D3Q19) reste.
  let cm3 = (1.0 - 0.5 * om) * 4.5 * om;
  let sxx = cm3 * u.x * u.x * pxx;
  let syy = cm3 * u.y * u.y * pyy;
  let szz = cm3 * u.z * u.z * pzz;
{M3_CLOSE}
  // Régularisation récursive (Malaspinas 2015) : moments d'ordre 3 tirés de
  // ρ, u et de la partie hors équilibre d'ordre 2 après collision.
  let bxx = k * pxx{M3_OPEN} + sxx{M3_CLOSE};
  let byy = k * pyy{M3_OPEN} + syy{M3_CLOSE};
  let bzz = k * pzz{M3_OPEN} + szz{M3_CLOSE};
  let bxy = k * pxy;
  let bxz = k * pxz;
  let byz = k * pyz;
  let axxy = rho * u.x * u.x * u.y + u.y * bxx + 2.0 * u.x * bxy;
  let axyy = rho * u.x * u.y * u.y + u.x * byy + 2.0 * u.y * bxy;
  let axxz = rho * u.x * u.x * u.z + u.z * bxx + 2.0 * u.x * bxz;
  let axzz = rho * u.x * u.z * u.z + u.x * bzz + 2.0 * u.z * bxz;
  let ayyz = rho * u.y * u.y * u.z + u.z * byy + 2.0 * u.y * byz;
  let ayzz = rho * u.y * u.z * u.z + u.y * bzz + 2.0 * u.z * byz;
  // Bandes absorbantes : rappel vers l'écoulement amont près des bords
  // latéraux (périodiques), de l'entrée et dans la zone éponge de sortie.
  let sig = absorb(x, y, z);
  var o: array<f32, 19>;
{coll}
  // Bandes absorbantes : rappel vers l'écoulement amont près des bords.
  if (sig > 0.0) {{
{absorb_code}  }}
{wforce}{out}{fstore}
  if (P.writeMacro == 1u || P.sa != 0u{surf_macro}) {{
    macroBuf[idx] = vec4<f32>(u, rho);
  }}
}}
",
        en = pop16_enable(p16),
        // Paroi volumique : vitesses de tous les pas pour le relevé des surfels.
        surf_macro = if surf { " || true" } else { "" },
        fstore = if surf {
            ""
        } else {
            "
  if (hasWall) {
    force[wk] += vec4<f32>(fsum, 0.0);
  }"
        },
        wforce = if bouzidi {
            let mut a = String::new();
            for i in 1..19 {
                let [cx, cy, cz] = C3[i];
                let _ = writeln!(
                    a,
                    "    if (qa[{k}] > 0.0) {{ aw += {w6} * abs(dot(vec3<f32>({x}, {y}, {z}), oct_decode(ca[{k}]))); }}",
                    k = i - 1,
                    w6 = f(6.0 * w3(i)),
                    x = f(cx as f64),
                    y = f(cy as f64),
                    z = f(cz as f64)
                );
            }
            let mut d = String::new();
            for i in 1..19 {
                let [cx, cy, cz] = C3[i];
                let _ = writeln!(d, "      o[{i}] += {} * dot(vec3<f32>({}, {}, {}), fw);", f(3.0 * w3(i)), f(cx as f64), f(cy as f64), f(cz as f64));
            }
            format!(
                "  // Contrainte pariétale imposée comme force (essai 16384) : la paroi
  // glisse librement (passe `wall3d`) et le nœud cède ρ u_τ² A le long de
  // sa vitesse tangentielle, A = 6 Σ w_i |c_i · n| sur ses liens de paroi
  // (aire exacte en moyenne pour un plan quelconque). Le corps reçoit
  // l'opposé : conservatif par construction.
  if ((P.essai & 16384u) != 0u && P.sa != 0u && nuw >= 0.0) {{
    var aw = 0.0;
{a}    let nw = wallGeo[wk].xyz;
    let ut3 = u - dot(u, nw) * nw;
    let utl = length(ut3);
    if (utl > 1e-9) {{
      let ts = utau[4u * wk];
      // Pas d'inversion de la vitesse tangentielle en un pas.
      let fm = min(rho * ts * ts * aw, 0.5 * rho * utl);
      let fw = -ut3 * (fm / utl);
{d}      fsum -= fw;
    }}
  }}
"
            )
        } else {
            String::new()
        },
        sa_fns = SA_FNS,
        ke_fns = KE_FNS,
        lam_fns = LAM_FNS,
        head = if walls {
            "  let kw = gid.x + gid.y * nwg.x * 64u;\n  if (kw >= WALLN_C) { return; }\n  let idx = wcell[kw];\n  let wk0 = wall[idx];"
        } else {
            "  let t = gid.x + gid.y * nwg.x * WGu;\n  var idx = t;\n  // Grille creuse : liste des mailles à calculer (après les pariétales).\n  if (ACTN_C > 0u) {\n    if (t >= ACTN_C) { return; }\n    idx = wcell[WALLN_C + t];\n  } else if (idx >= P.N) { return; }\n  // Liste : mailles fluides sans lien de paroi, type connu sans le lire.\n  var wk0 = 0u;\n  if (ACTN_C == 0u && ((plain[idx >> 5u] >> (idx & 31u)) & 1u) == 0u) { wk0 = wall[idx]; }"
        },
        wg = if walls { 64 } else { step_wg() },
        // Mailles pariétales : fractions q et codes des normales des 18 liens
        // lus d'un coup (contigus), au lieu d'un par direction derrière ses
        // branches (−27 % sur la collision pariétale).
        prefetch = if walls && surf {
            format!(
                "  let qa = array<f32, 18>({});\n",
                (0..18).map(|k| format!("qbuf[wk * 18u + {k}u]")).collect::<Vec<_>>().join(", ")
            )
        } else if walls {
            format!(
                "  let qa = array<f32, 18>({});
  let ch = arrayLength(&qbuf) / 2u + wk * 18u;
  let ca = array<f32, 18>({});
",
                (0..18).map(|k| format!("qbuf[wk * 18u + {k}u]")).collect::<Vec<_>>().join(", "),
                (0..18).map(|k| format!("qbuf[ch + {k}u]")).collect::<Vec<_>>().join(", ")
            )
        } else {
            String::new()
        },
        skip = if walls { "" } else { "  // Maille pariétale : passe à part.\n  if (wk0 != 0u) { return; }" },
        rho = sum_expr(|_| 1),
        jx = sum_expr(|i| C3[i][0]),
        jy = sum_expr(|i| C3[i][1]),
        jz = sum_expr(|i| C3[i][2]),
        pxx = sum_expr(|i| dotc(i, 0, 0)),
        pyy = sum_expr(|i| dotc(i, 1, 1)),
        pzz = sum_expr(|i| dotc(i, 2, 2)),
        pxy = sum_expr(|i| dotc(i, 0, 1)),
        pxz = sum_expr(|i| dotc(i, 0, 2)),
        pyz = sum_expr(|i| dotc(i, 1, 2)),
    )
}

/// Paroi volumique : contrainte tangentielle relevée au-dessus des surfels
/// (paroi à sa place, ordre 2), sauf `SOUFFLERIE_SURFELS_O1=1` (rebond et
/// loi de paroi de la maille pariétale, schéma d'ordre 1 d'origine).
pub fn surfels_o2() -> bool {
    std::env::var("SOUFFLERIE_SURFELS_O1").map_or(true, |v| v.is_empty() || v == "0")
}

/// Passe des surfels (paroi volumique, `surfel3d`), un fil par surfel, sur
/// l'état d'avant le pas (populations après collision) : rassemblement des
/// populations entrantes de chaque paire de directions, rebond, loi de paroi
/// (part spéculaire de la quantité de mouvement tangentielle), effort sur le
/// corps, et Γᵒᵘᵗ/P (populations décalées de leur poids) pour la propagation
/// des mailles réceptrices. Tampon `force` : effort du surfel s en s, puis
/// Γᵒᵘᵗ/P de ses paires en 4 n + 12 s + p (cases flottantes).
fn wgsl_surfel(fp16: bool, p16: bool, n: u32, o2: bool) -> String {
    let half = if fp16 { 5 } else { 10 };
    let load = if p16 {
        "fn pop(i: u32, c: u32) -> f32 {
  if (i < 10u) { return unpack2x16float(u32(srcA[i * P.N + c])).x; }
  return unpack2x16float(u32(srcB[(i - 10u) * P.N + c])).x;
}".to_string()
    } else if fp16 {
        format!(
            "{hs}fn pop(i: u32, c: u32) -> f32 {{
  let h = HS[i];
  let w = h >> 1u;
  var v: vec2<f32>;
  if (w < {half}u) {{ v = unpack2x16float(srcA[w * P.N + c]); }} else {{ v = unpack2x16float(srcB[(w - {half}u) * P.N + c]); }}
  return select(v.x, v.y, (h & 1u) == 1u);
}}",
            hs = hs_wgsl()
        )
    } else {
        format!(
            "fn pop(i: u32, c: u32) -> f32 {{
  if (i < {half}u) {{ return bitcast<f32>(srcA[i * P.N + c]); }}
  return bitcast<f32>(srcB[(i - {half}u) * P.N + c]);
}}"
        )
    };
    let mut gather = String::new();
    let mut corr = String::new();
    let mut fin = String::new();
    for p in 0..9 {
        let (a, b) = (2 * p + 1, 2 * p + 2);
        let [cx, cy, cz] = C3[a];
        let cv = format!("vec3<f32>({}, {}, {})", f(cx as f64), f(cy as f64), f(cz as f64));
        let w = f(w3(a));
        let _ = write!(
            gather,
            "  {{
    let c1 = {cv};
    let cn = dot(c1, nrm);
    // Entrante : c₁ si c₁·n < 0, sinon son opposée.
    let ii = select({a}u, {b}u, cn > 0.0);
    var sum = 0.0;
    let e1 = gst[9u * s + {p1}u];
    for (var e = gst[9u * s + {p}u]; e < e1; e++) {{
      sum += bitcast<f32>(gat[2u * e + 1u]) * (pop(ii, gat[2u * e]) + {w});
    }}
    gin[{p}] = sum;
    pv[{p}] = abs(cn) * area;
    cin[{p}] = select(c1, -c1, cn > 0.0);
    fin += cin[{p}] * sum;
    mw += pv[{p}] * {w};
    mi += sum;
  }}
",
            p1 = p + 1,
        );
        let _ = write!(
            corr,
            "        d[{p}] = pv[{p}] * {w} * (lam.x - dot(lam.yzw, cin[{p}]));
        if (gin[{p}] + d[{p}] < 0.0) {{ sc = min(sc, gin[{p}] / max(-d[{p}], 1e-30)); }}
"
        );
        let _ = write!(
            fin,
            "  gout[{p}] = gin[{p}] + sc * d[{p}];
  ff += cin[{p}] * (gin[{p}] + gout[{p}]);
  {{
    let inv = sgeo[8u * s + {row}u].{comp};
    if (inv > 0.0) {{
      let k = 4u * SURFN + 12u * s + {p}u;
      force[k >> 2u][k & 3u] = gout[{p}] * inv - {w};
    }}
  }}
",
            row = 5 + p / 4,
            comp = ["x", "y", "z", "w"][p % 4],
        );
    }
    // Paroi à sa place (ordre 2) : contrainte tangentielle tirée du relevé
    // de la vitesse au-dessus du surfel, à la place de celle du rebond.
    let o2_block = if o2 {
        format!(
            "  // Contrainte tangentielle imposée, relevée au-dessus du surfel
  // (`surfel3d::probe_fit`, nœuds entièrement fluides, profil nul à la
  // paroi) : en laminaire ρ ν (∂u/∂η) A, la pente à la paroi ; en
  // turbulent ρ u_τ² A, u_τ de la loi de Spalding à la hauteur WM_H, comme
  // pour le rebond de Bouzidi. Le rebond renvoie masse et quantité de
  // mouvement normale ; la perturbation δⱼ = Pⱼ wⱼ (λ₀ + Λ·cⱼ) ramène son
  // effort tangentiel 2 Σ cᵢᵗ Γᵢⁱⁿ à l'effort visé. La place de la paroi
  // n'est plus celle, mélangée, des voxels coupés, mais celle du relevé,
  // exacte. Arête vive (loi de paroi) : rebond seul.
  {{
    let p0 = pst[s];
    let p1 = pst[s + 1u];
    var edge = false;
    if (home != 0xffffffffu) {{ edge = wallGeo[home].w < 0.0; }}
    if (p1 > p0 && area > 0.0 && !(P.sa != 0u && edge)) {{
      var ga = vec3<f32>(0.0);
      var gh = vec3<f32>(0.0);
      for (var e = p0; e < p1; e++) {{
        let uu = macroBuf[prb[3u * e]].xyz;
        ga += bitcast<f32>(prb[3u * e + 1u]) * uu;
        gh += bitcast<f32>(prb[3u * e + 2u]) * uu;
      }}
      ga -= dot(ga, nrm) * nrm;
      gh -= dot(gh, nrm) * nrm;
      let rho = select(1.0, mi / mw, mw > 0.0);
      var ft = rho * P.nu * area * ga;
      if (P.sa != 0u) {{
        let uh = length(gh);
        let ut = spalding_utau(uh, WM_H);
        ft = select(vec3<f32>(0.0), gh * (rho * ut * ut * area / max(uh, 1e-30)), uh > 1e-20);
      }}
      let fb = 2.0 * (fin - dot(fin, nrm) * nrm);
      let dm = fb - ft;
      let lam = (dm.x * sgeo[8u * s + 2u] + dm.y * sgeo[8u * s + 3u] + dm.z * sgeo[8u * s + 4u]) / area;
      sc = 1.0;
{corr}      done = true;
    }}
  }}
"
        )
    } else {
        String::new()
    };
    let en = pop16_enable(p16);
    let pt = if p16 { "u16" } else { "u32" };
    let sa_fns = SA_FNS;
    format!(
        "{en}{PARAMS_WGSL}
@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read> srcA: array<{pt}>;
@group(0) @binding(2) var<storage, read> srcB: array<{pt}>;
@group(0) @binding(3) var<storage, read_write> force: array<vec4<f32>>;
// Par surfel : (n, A), (centre, maille pariétale de u_τ), les trois
// colonnes de M (`surfel3d::wall_law_matrix`), puis 1 / Σ des poids de
// dispersion de chaque paire (`SurfelSetup::scat_inv`, trois vec4).
@group(0) @binding(4) var<storage, read> sgeo: array<vec4<f32>>;
// Rassemblement : début des entrées de chaque paire, entrées (maille, poids).
@group(0) @binding(5) var<storage, read> gst: array<u32>;
@group(0) @binding(6) var<storage, read> gat: array<u32>;
@group(0) @binding(7) var<storage, read> utau: array<f32>;
@group(0) @binding(8) var<storage, read> wallGeo: array<vec4<f32>>;
// Relevé de la vitesse (`SurfelSetup::probe`) : vitesses du pas précédent,
// début des entrées de chaque surfel, entrées (maille, α, β).
@group(0) @binding(9) var<storage, read> macroBuf: array<vec4<f32>>;
@group(0) @binding(10) var<storage, read> pst: array<u32>;
@group(0) @binding(11) var<storage, read> prb: array<u32>;
const SURFN = {n}u;
{sa_fns}
{load}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {{
  let s = gid.x + gid.y * nwg.x * 64u;
  if (s >= SURFN) {{ return; }}
  let g0 = sgeo[8u * s];
  let nrm = g0.xyz;
  let area = g0.w;
  let home = bitcast<u32>(sgeo[8u * s + 1u].w);
  var gin: array<f32, 9>;
  var pv: array<f32, 9>;
  var cin: array<vec3<f32>, 9>;
  var fin = vec3<f32>(0.0);
  var mw = 0.0;
  var mi = 0.0;
{gather}
  // Rebond : Γᵒᵘᵗ = Γⁱⁿ dans la direction opposée de chaque paire. Loi de
  // paroi (maille pariétale de la loi, ni arête vive) : l'effort tangentiel
  // du rebond, 2 Σ cᵢᵗ Γᵢⁱⁿ, ramené à ρ u_τ² A (au plus celui du rebond) par
  // la perturbation δⱼ = Pⱼ wⱼ (λ₀ + Λ·cⱼ) qui rend la quantité de mouvement
  // tangentielle (part spéculaire) sans toucher masse ni pression.
  var d: array<f32, 9>;
  var sc = 0.0;
  var done = false;
{o2_block}  if (!done && P.sa != 0u && home != 0xffffffffu && area > 0.0) {{
    if (wallGeo[home].w > 0.0) {{
      let ut = utau[4u * home];
      let rho = select(1.0, mi / mw, mw > 0.0);
      let fb = 2.0 * (fin - dot(fin, nrm) * nrm);
      let fl = length(fb);
      if (fl > 1e-20) {{
        let frac = min(1.0, rho * ut * ut * area / fl);
        let dm = fb * (1.0 - frac);
        let lam = (dm.x * sgeo[8u * s + 2u] + dm.y * sgeo[8u * s + 3u] + dm.z * sgeo[8u * s + 4u]) / area;
        // Positivité des populations sortantes : perturbation réduite
        // d'autant (même direction, contrainte plus forte).
        sc = 1.0;
{corr}      }}
    }}
  }}
  var gout: array<f32, 9>;
  var ff = vec3<f32>(0.0);
{fin}  force[s] += vec4<f32>(ff, 0.0);
}}
"
    )
}

fn wgsl_init(fp16: bool, p16: bool) -> String {
    let words = if fp16 { 10 } else { 19 };
    let half = if fp16 { 5 } else { 10 };
    let eq: String = (0..19)
        .map(|i| {
            let c = f(C3[i][0] as f64);
            format!(
                "let e{i} = {} * (3.0 * ({c} * u0) + 4.5 * ({c} * u0) * ({c} * u0) - 1.5 * u0 * u0);\n  ",
                f(w3(i))
            )
        })
        .collect();
    let put: String = if p16 {
        (0..19)
            .map(|i| {
                let v = format!("u16(pack2x16float(vec2<f32>(e{i}, 0.0)) & 0xffffu)");
                if i < 10 {
                    format!("A[{i}u * P.N + idx] = {v};
  ")
                } else {
                    format!("B[{}u * P.N + idx] = {v};
  ", i - 10)
                }
            })
            .collect()
    } else if fp16 {
        let hd = hdir();
        let e = |h: usize| hd[h].map_or("0.0".to_string(), |i| format!("e{i}"));
        (0..words)
            .map(|w| {
                let v = format!("pack2x16float(vec2<f32>({}, {}))", e(2 * w), e(2 * w + 1));
                if w < half {
                    format!("A[{w}u * P.N + idx] = {v};\n  ")
                } else {
                    format!("B[{}u * P.N + idx] = {v};\n  ", w - half)
                }
            })
            .collect()
    } else {
        (0..19)
            .map(|i| {
                if i < half {
                    format!("A[{i}u * P.N + idx] = bitcast<u32>(e{i});\n  ")
                } else {
                    format!("B[{}u * P.N + idx] = bitcast<u32>(e{i});\n  ", i - half)
                }
            })
            .collect()
    };
    let pt = if p16 { "u16" } else { "u32" };
    let en = pop16_enable(p16);
    format!(
        "{en}{PARAMS_WGSL}
@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read_write> A: array<{pt}>;
@group(0) @binding(2) var<storage, read_write> B: array<{pt}>;
@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {{
  let idx = gid.x + gid.y * nwg.x * 64u;
  if (idx >= P.N) {{ return; }}
  let u0 = P.uIn;
  {eq}
  {put}
}}
"
    )
}

/// Coupe : grandeur sur un plan, calculée sur le GPU depuis le champ macroscopique.
const SLICE_WGSL: &str = "
struct SliceParams {
  nx: u32, ny: u32, nz: u32, axis: u32,
  pos: u32, field: u32, w: u32, h: u32,
  u: f32, l: f32, pad0: f32, pad1: f32,
};
@group(0) @binding(0) var<uniform> S: SliceParams;
@group(0) @binding(1) var<storage, read> macroBuf: array<vec4<f32>>;
@group(0) @binding(2) var<storage, read> wall: array<u32>;
@group(0) @binding(3) var<storage, read_write> outv: array<f32>;

fn id(x: i32, y: i32, z: i32) -> u32 {
  let ny = i32(S.ny);
  let nz = i32(S.nz);
  let xx = clamp(x, 0, i32(S.nx) - 1);
  let yy = (y + ny) % ny;
  let zz = (z + nz) % nz;
  return u32(xx) + S.nx * (u32(yy) + S.ny * u32(zz));
}

// Vitesse (nulle dans le corps : paroi immobile).
fn vel(x: i32, y: i32, z: i32) -> vec3<f32> {
  let c = id(x, y, z);
  if (wall[c] == 0xffffffffu) { return vec3<f32>(0.0); }
  return macroBuf[c].xyz;
}

@compute @workgroup_size(8, 8)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  if (gid.x >= S.w || gid.y >= S.h) { return; }
  var x: i32; var y: i32; var z: i32;
  let a = i32(gid.x);
  let b = i32(gid.y);
  let p = i32(S.pos);
  if (S.axis == 0u) { x = p; y = a; z = b; }
  else if (S.axis == 1u) { x = a; y = p; z = b; }
  else { x = a; y = b; z = p; }
  let o = gid.x + S.w * gid.y;
  if (wall[id(x, y, z)] == 0xffffffffu) { outv[o] = 1e30; return; }
  let u = vel(x, y, z);
  if (S.field == 0u) { outv[o] = length(u) / S.u; return; }
  // Gradient par différences centrées : G[i][j] = ∂u_i/∂x_j.
  let dx = 0.5 * (vel(x + 1, y, z) - vel(x - 1, y, z));
  let dy = 0.5 * (vel(x, y + 1, z) - vel(x, y - 1, z));
  let dz = 0.5 * (vel(x, y, z + 1) - vel(x, y, z - 1));
  let s = S.l / S.u;
  if (S.field == 1u) {
    // Composante de la vorticité normale au plan.
    var w: f32;
    if (S.axis == 0u) { w = dy.z - dz.y; }
    else if (S.axis == 1u) { w = dz.x - dx.z; }
    else { w = dx.y - dy.x; }
    outv[o] = w * s;
    return;
  }
  // Critère Q = ½ (‖Ω‖² − ‖S‖²).
  let sxx = dx.x; let syy = dy.y; let szz = dz.z;
  let sxy = 0.5 * (dy.x + dx.y); let sxz = 0.5 * (dz.x + dx.z); let syz = 0.5 * (dz.y + dy.z);
  let oxy = 0.5 * (dy.x - dx.y); let oxz = 0.5 * (dz.x - dx.z); let oyz = 0.5 * (dz.y - dy.z);
  let ss = sxx * sxx + syy * syy + szz * szz + 2.0 * (sxy * sxy + sxz * sxz + syz * syz);
  let oo = 2.0 * (oxy * oxy + oxz * oxz + oyz * oyz);
  outv[o] = 0.5 * (oo - ss) * s * s;
}
";

/// Dérivées du champ : vorticité ω (xyz) et critère Q (w), par maille.
const DERIVE_WGSL: &str = "
struct D { nx: u32, ny: u32, nz: u32, n: u32 };
@group(0) @binding(0) var<uniform> P: D;
@group(0) @binding(1) var<storage, read> macroBuf: array<vec4<f32>>;
@group(0) @binding(2) var<storage, read> wall: array<u32>;
@group(0) @binding(3) var<storage, read_write> deriv: array<vec4<f32>>;

fn id(p: vec3<i32>) -> u32 {
  let q = clamp(p, vec3<i32>(0), vec3<i32>(i32(P.nx) - 1, i32(P.ny) - 1, i32(P.nz) - 1));
  return u32(q.x) + P.nx * (u32(q.y) + P.ny * u32(q.z));
}
fn solid(p: vec3<i32>) -> bool { return wall[id(p)] == 0xffffffffu; }
fn U(p: vec3<i32>) -> vec3<f32> {
  if (solid(p)) { return vec3<f32>(0.0); }
  return macroBuf[id(p)].xyz;
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let i = gid.x + gid.y * nwg.x * 64u;
  if (i >= P.n) { return; }
  let p = vec3<i32>(i32(i % P.nx), i32((i / P.nx) % P.ny), i32(i / (P.nx * P.ny)));
  let ex = vec3<i32>(1, 0, 0);
  let ey = vec3<i32>(0, 1, 0);
  let ez = vec3<i32>(0, 0, 1);
  let dx = 0.5 * (U(p + ex) - U(p - ex));
  let dy = 0.5 * (U(p + ey) - U(p - ey));
  let dz = 0.5 * (U(p + ez) - U(p - ez));
  let w = vec3<f32>(dy.z - dz.y, dz.x - dx.z, dx.y - dy.x);
  let s2 = dx.x * dx.x + dy.y * dy.y + dz.z * dz.z
    + 0.5 * ((dx.y + dy.x) * (dx.y + dy.x) + (dx.z + dz.x) * (dx.z + dz.x) + (dy.z + dz.y) * (dy.z + dz.y));
  var q = 0.5 * (0.5 * dot(w, w) - s2);
  // Au contact du corps, les différences finies traversent la paroi : Q y est artificiel.
  // Exclusion sur deux mailles : la couche limite, mal résolue, y donne des Q parasites.
  var near = solid(p);
  for (var k = 1; k <= 2; k++) {
    let a = vec3<i32>(k, 0, 0);
    let b = vec3<i32>(0, k, 0);
    let c = vec3<i32>(0, 0, k);
    near = near || solid(p + a) || solid(p - a) || solid(p + b) || solid(p - b) || solid(p + c) || solid(p - c);
  }
  near = near || solid(p + ex + ey) || solid(p + ex - ey) || solid(p - ex + ey) || solid(p - ex - ey)
    || solid(p + ex + ez) || solid(p + ex - ez) || solid(p - ex + ez) || solid(p - ex - ez)
    || solid(p + ey + ez) || solid(p + ey - ez) || solid(p - ey + ez) || solid(p - ey - ez);
  if (near) {
    q = 0.0;
  }
  deriv[i] = vec4<f32>(w, q);
}
";

/// Image de la veine par lancer de rayons : maquette ombrée, plan de coupe,
/// isosurfaces du critère Q translucides teintées par la vitesse.
const RENDER_WGSL: &str = "
struct R {
  eye: vec4<f32>,
  fwd: vec4<f32>,
  right: vec4<f32>,
  up: vec4<f32>,
  dims: vec4<f32>,
  bg: vec4<f32>,
  body: vec4<f32>,
  plane: vec4<f32>,
  iso: vec4<f32>,
  size: vec4<f32>,
  sdf0: vec4<f32>,
  sdfn: vec4<f32>,
  vorg: vec4<f32>,
  vdims: vec4<f32>,
};
@group(0) @binding(0) var<uniform> V: R;
// Distance signée à la maquette (affichage) : origine et pas dans sdf0, taille
// dans sdfn.xyz, épaississement affiché dans sdfn.w.
@group(0) @binding(5) var<storage, read> sdf: array<f32>;
// Type de chaque maille du calcul, deux bits : 0 fluide, 1 solide, 2 fluide
// portant seulement des liens de paroi mince. Mode d'affichage de la
// maquette dans size.z : 0 lisse, 1 lisse avec la grille, 2 mailles du calcul.
@group(0) @binding(6) var<storage, read> vox: array<u32>;
@group(0) @binding(1) var<storage, read> macroBuf: array<vec4<f32>>;
@group(0) @binding(3) var<storage, read> deriv: array<vec4<f32>>;
@group(0) @binding(4) var<storage, read_write> img: array<u32>;

fn seqMap(t0: f32) -> vec3<f32> {
  let t = clamp(t0, 0.0, 1.0);
  let a = vec3<f32>(0.05, 0.11, 0.27);
  let b = vec3<f32>(0.10, 0.45, 0.62);
  let c = vec3<f32>(0.42, 0.78, 0.62);
  let d = vec3<f32>(0.99, 0.88, 0.35);
  if (t < 0.33) { return mix(a, b, t / 0.33); }
  if (t < 0.66) { return mix(b, c, (t - 0.33) / 0.33); }
  return mix(c, d, (t - 0.66) / 0.34);
}
fn divMap(t0: f32) -> vec3<f32> {
  let t = clamp(t0, -1.0, 1.0);
  // Valeur nulle : teinte du fond, éclaircie (thème clair) ou non (sombre).
  let mid = clamp(V.bg.rgb * 1.04 + vec3<f32>(0.01), vec3<f32>(0.0), vec3<f32>(1.0));
  if (t < 0.0) { return mix(mid, vec3<f32>(0.13, 0.36, 0.70), -t); }
  return mix(mid, vec3<f32>(0.78, 0.20, 0.15), t);
}

fn id(p: vec3<i32>) -> u32 {
  let n = vec3<i32>(V.dims.xyz);
  let q = clamp(p, vec3<i32>(0), n - vec3<i32>(1));
  return u32(q.x) + u32(n.x) * (u32(q.y) + u32(n.y) * u32(q.z));
}
// Mailles affichées (grille la plus fine qui entoure la maquette, ou la veine) :
// position dans cette grille d'un point en mailles de la veine.
fn toVox(p: vec3<f32>) -> vec3<f32> { return (p - V.vorg.xyz) * V.vorg.w; }
fn voxKind(c: vec3<i32>) -> u32 {
  let n = vec3<i32>(V.vdims.xyz);
  if (any(c < vec3<i32>(0)) || any(c >= n)) { return 0u; }
  let i = u32(c.x) + u32(n.x) * (u32(c.y) + u32(n.y) * u32(c.z));
  return (vox[i >> 4u] >> ((i & 15u) * 2u)) & 3u;
}

// Première maille non fluide le long du rayon (traversée d'Amanatides et Woo,
// mailles centrées sur les entiers), en mailles de la grille affichée :
// (t, axe de la face d'entrée, type) ; t en mailles de la veine.
fn voxHit(ov: vec3<f32>, dir: vec3<f32>, t0v: f32, t1v: f32) -> vec3<f32> {
  let s = V.vorg.w;
  let o = toVox(ov);
  // Rayon limité à la boîte de la grille affichée.
  let inv = 1.0 / dir;
  let ba = (vec3<f32>(-0.5) - o) * inv;
  let bb = (V.vdims.xyz - vec3<f32>(0.5) - o) * inv;
  let t0 = max(t0v * s, max(max(min(ba.x, bb.x), min(ba.y, bb.y)), min(ba.z, bb.z)));
  let t1 = min(t1v * s, min(min(max(ba.x, bb.x), max(ba.y, bb.y)), max(ba.z, bb.z)));
  if (t0 > t1) { return vec3<f32>(1e30, 0.0, 0.0); }
  let r = voxWalk(o, dir, t0, t1);
  return vec3<f32>(select(r.x / s, 1e30, r.x > 1e29), r.y, r.z);
}
fn voxWalk(o: vec3<f32>, dir: vec3<f32>, t0: f32, t1: f32) -> vec3<f32> {
  let p = o + dir * (t0 + 1e-4);
  var c = vec3<i32>(round(p));
  let st = vec3<i32>(sign(dir));
  var tmax = vec3<f32>(1e30);
  var tdel = vec3<f32>(1e30);
  for (var a = 0; a < 3; a++) {
    if (abs(dir[a]) > 1e-9) {
      let b = f32(c[a]) + 0.5 * f32(st[a]);
      tmax[a] = (b - o[a]) / dir[a];
      tdel[a] = 1.0 / abs(dir[a]);
    }
  }
  var t = t0;
  var axis = 0;
  for (var k = 0; k < 20000; k++) {
    let kd = voxKind(c);
    if (kd != 0u) { return vec3<f32>(t, f32(axis), f32(kd)); }
    if (tmax.x <= tmax.y && tmax.x <= tmax.z) { axis = 0; } else if (tmax.y <= tmax.z) { axis = 1; } else { axis = 2; }
    t = tmax[axis];
    if (t > t1) { break; }
    c[axis] += st[axis];
    tmax[axis] += tdel[axis];
  }
  return vec3<f32>(1e30, 0.0, 0.0);
}

fn sd(i: vec3<i32>) -> f32 {
  let n = vec3<i32>(V.sdfn.xyz);
  let q = clamp(i, vec3<i32>(0), n - vec3<i32>(1));
  return sdf[u32(q.x) + u32(n.x) * (u32(q.y) + u32(n.y) * u32(q.z))];
}
// Distance signée (mailles) ; hors de la grille, distance à sa boîte.
fn S(p: vec3<f32>) -> f32 {
  let g = (p - V.sdf0.xyz) / V.sdf0.w;
  let n = V.sdfn.xyz - vec3<f32>(1.0);
  let outside = max(max(-g, g - n), vec3<f32>(0.0));
  if (any(outside > vec3<f32>(0.0))) {
    return length(outside) * V.sdf0.w + 1.0;
  }
  let b = vec3<i32>(floor(g));
  let f = g - floor(g);
  let c00 = mix(sd(b), sd(b + vec3<i32>(1, 0, 0)), f.x);
  let c10 = mix(sd(b + vec3<i32>(0, 1, 0)), sd(b + vec3<i32>(1, 1, 0)), f.x);
  let c01 = mix(sd(b + vec3<i32>(0, 0, 1)), sd(b + vec3<i32>(1, 0, 1)), f.x);
  let c11 = mix(sd(b + vec3<i32>(0, 1, 1)), sd(b + vec3<i32>(1, 1, 1)), f.x);
  return mix(mix(c00, c10, f.y), mix(c01, c11, f.y), f.z);
}
fn D(p: vec3<f32>) -> vec4<f32> {
  let b = vec3<i32>(floor(p));
  let f = p - floor(p);
  let c00 = mix(deriv[id(b)], deriv[id(b + vec3<i32>(1, 0, 0))], f.x);
  let c10 = mix(deriv[id(b + vec3<i32>(0, 1, 0))], deriv[id(b + vec3<i32>(1, 1, 0))], f.x);
  let c01 = mix(deriv[id(b + vec3<i32>(0, 0, 1))], deriv[id(b + vec3<i32>(1, 0, 1))], f.x);
  let c11 = mix(deriv[id(b + vec3<i32>(0, 1, 1))], deriv[id(b + vec3<i32>(1, 1, 1))], f.x);
  return mix(mix(c00, c10, f.y), mix(c01, c11, f.y), f.z);
}
fn Q(p: vec3<f32>) -> f32 { return D(p).w; }
fn Mv(p: vec3<f32>) -> vec3<f32> {
  let b = vec3<i32>(floor(p));
  let f = p - floor(p);
  let c00 = mix(macroBuf[id(b)].xyz, macroBuf[id(b + vec3<i32>(1, 0, 0))].xyz, f.x);
  let c10 = mix(macroBuf[id(b + vec3<i32>(0, 1, 0))].xyz, macroBuf[id(b + vec3<i32>(1, 1, 0))].xyz, f.x);
  let c01 = mix(macroBuf[id(b + vec3<i32>(0, 0, 1))].xyz, macroBuf[id(b + vec3<i32>(1, 0, 1))].xyz, f.x);
  let c11 = mix(macroBuf[id(b + vec3<i32>(0, 1, 1))].xyz, macroBuf[id(b + vec3<i32>(1, 1, 1))].xyz, f.x);
  return mix(mix(c00, c10, f.y), mix(c01, c11, f.y), f.z);
}

@compute @workgroup_size(8, 8)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  let w = u32(V.size.x);
  let h = u32(V.size.y);
  if (gid.x >= w || gid.y >= h) { return; }
  let sx = (f32(gid.x) + 0.5) / f32(w) * 2.0 - 1.0;
  let sy = 1.0 - (f32(gid.y) + 0.5) / f32(h) * 2.0;
  let th = V.up.w;
  let o = V.eye.xyz;
  let dir = normalize(V.fwd.xyz + V.right.xyz * (sx * th * V.right.w) + V.up.xyz * (sy * th));
  let lo = vec3<f32>(0.0);
  let hi = V.dims.xyz - vec3<f32>(1.0);
  let inv = 1.0 / dir;
  let t0v = (lo - o) * inv;
  let t1v = (hi - o) * inv;
  let tmin = max(max(min(t0v.x, t1v.x), min(t0v.y, t1v.y)), min(t0v.z, t1v.z));
  let tmax = min(min(max(t0v.x, t1v.x), max(t0v.y, t1v.y)), max(t0v.z, t1v.z));
  var col = V.bg.rgb;
  let t0 = max(tmin, 0.0);
  if (t0 < tmax) {
    // Plan de coupe.
    let ax = i32(V.plane.x);
    var ts = 1e30;
    if (abs(dir[ax]) > 1e-6) {
      let t = (V.plane.y - o[ax]) / dir[ax];
      if (t >= t0 && t <= tmax) { ts = t; }
    }
    let tend = min(tmax, ts);
    let stp = V.iso.z;
    var t = t0;
    var tp = t0;
    var tb = 1e30;
    var ti = 1e30;
    var tip = 0.0;
    let level = V.iso.x;
    let isoOn = V.iso.y > 0.5;
    let mode = u32(V.size.z);
    var vh = vec3<f32>(1e30, 0.0, 0.0);
    if (mode == 2u) {
      vh = voxHit(o, dir, t0, tend);
      tb = vh.x;
    }
    for (var k = 0; k < 4000; k++) {
      if (t > min(tend, tb)) { break; }
      let p = o + dir * t;
      if (isoOn && ti > 1e29 && Q(p) > level) { ti = t; tip = tp; }
      let ds = select(S(p), 1e30, mode == 2u);
      if (ds < V.sdfn.w) {
        // Surface exacte de la forme : bisection sur la distance signée.
        var a = tp;
        var b = t;
        for (var j = 0; j < 10; j++) {
          let m = 0.5 * (a + b);
          if (S(o + dir * m) < V.sdfn.w) { b = m; } else { a = m; }
        }
        tb = b;
        break;
      }
      tp = t;
      // Pas réduit près de la maquette : le bord de fuite est plus mince qu'une maille.
      t += clamp(ds - V.sdfn.w, 0.05, stp);
    }
    if (tb < 1e29 && mode == 2u) {
      // Mailles du calcul : cubes ombrés par face, arêtes marquées ; liens de
      // paroi mince en orange.
      let p = o + dir * tb;
      let a = i32(vh.y);
      var n = vec3<f32>(0.0);
      n[a] = -sign(dir[a]);
      let lam = abs(dot(n, -dir));
      let key = max(0.0, dot(n, normalize(vec3<f32>(-0.4, -0.3, 0.85))));
      var base = V.body.rgb;
      if (vh.z > 1.5) { base = vec3<f32>(0.93, 0.55, 0.18); }
      col = base * (0.35 + 0.45 * lam + 0.25 * key);
      let fr = abs(fract(toVox(p) + vec3<f32>(0.5)) - vec3<f32>(0.5));
      var e = 1.0;
      for (var b = 0; b < 3; b++) {
        if (b != a) { e = min(e, 0.5 - fr[b]); }
      }
      // Arêtes d'environ un pixel quelle que soit la distance : largeur en
      // mailles affichées ≈ angle d'un pixel × distance.
      let px = 2.0 * V.up.w / V.size.y * tb * V.vorg.w;
      let ew = clamp(px * 0.8, 0.015, 0.12);
      // Mailles de moins de quatre pixels environ : arêtes estompées (moiré).
      let fade = 1.0 - smoothstep(0.15, 0.35, px);
      col = mix(col, col * 0.45, fade * (1.0 - smoothstep(ew * 0.5, ew, e)));
    } else if (tb < 1e29) {
      let p = o + dir * tb;
      let e = 0.5 * V.sdf0.w;
      let g = vec3<f32>(
        S(p + vec3<f32>(e, 0.0, 0.0)) - S(p - vec3<f32>(e, 0.0, 0.0)),
        S(p + vec3<f32>(0.0, e, 0.0)) - S(p - vec3<f32>(0.0, e, 0.0)),
        S(p + vec3<f32>(0.0, 0.0, e)) - S(p - vec3<f32>(0.0, 0.0, e)));
      let n = normalize(g + vec3<f32>(1e-9));
      let lam = abs(dot(n, -dir));
      let key = max(0.0, dot(n, normalize(vec3<f32>(-0.4, -0.3, 0.85))));
      col = V.body.rgb * (0.35 + 0.45 * lam + 0.25 * key);
      if (mode == 1u) {
        // Grille du calcul tracée sur la surface : plans x, y, z = k + ½.
        let fr = abs(fract(p + vec3<f32>(0.5)) - vec3<f32>(0.5));
        let e = 0.5 - max(max(fr.x, fr.y), fr.z);
        if (e < 0.04) { col = mix(col, V.bg.rgb, 0.5); }
      }
    } else if (ts < 1e29) {
      let p = o + dir * ts;
      let inBody = select(S(p) < V.sdfn.w, voxKind(vec3<i32>(round(toVox(p)))) != 0u, mode == 2u);
      if (inBody) {
        col = V.body.rgb * 0.6;
      } else if (V.plane.z < 0.5) {
        col = seqMap(length(Mv(p)) / (1.6 * V.plane.w));
      } else {
        let dv = D(p);
        var wn = dv.y;
        if (ax == 0) { wn = dv.x; } else if (ax == 2) { wn = dv.z; }
        col = divMap(wn * V.dims.w / V.plane.w / V.iso.w);
      }
    }
    if (ti < 1e29 && ti < tb) {
      var a = tip;
      var b = ti;
      for (var j = 0; j < 5; j++) {
        let m = 0.5 * (a + b);
        if (Q(o + dir * m) > level) { b = m; } else { a = m; }
      }
      let p = o + dir * b;
      let e = 1.0;
      let gq = vec3<f32>(
        Q(p + vec3<f32>(e, 0.0, 0.0)) - Q(p - vec3<f32>(e, 0.0, 0.0)),
        Q(p + vec3<f32>(0.0, e, 0.0)) - Q(p - vec3<f32>(0.0, e, 0.0)),
        Q(p + vec3<f32>(0.0, 0.0, e)) - Q(p - vec3<f32>(0.0, 0.0, e)));
      let n = -normalize(gq + vec3<f32>(1e-9));
      let lam = abs(dot(n, -dir));
      let base = seqMap(0.15 + 0.75 * length(Mv(p)) / (1.3 * V.plane.w));
      let c = base * (0.35 + 0.75 * lam) + vec3<f32>(0.25) * pow(lam, 24.0);
      col = mix(col, c, 0.88);
    }
  }
  img[gid.x + w * gid.y] = pack4x8unorm(vec4<f32>(col, 1.0));
}
";

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct RenderUniform {
    eye: [f32; 4],
    fwd: [f32; 4],
    right: [f32; 4],
    up: [f32; 4],
    dims: [f32; 4],
    bg: [f32; 4],
    body: [f32; 4],
    plane: [f32; 4],
    iso: [f32; 4],
    size: [f32; 4],
    sdf0: [f32; 4],
    sdfn: [f32; 4],
    /// Grille des mailles affichées : origine (mailles de la veine) et facteur
    /// de raffinement, puis dimensions.
    vorg: [f32; 4],
    vdims: [f32; 4],
}

/// Caméra orbitale autour d'un point (mailles) : azimut et élévation en
/// degrés, distance en mailles, champ vertical de 40°.
#[derive(Clone, Copy, Debug)]
pub struct Camera {
    pub eye: [f64; 3],
    pub fwd: [f64; 3],
    pub right: [f64; 3],
    pub up: [f64; 3],
    pub tan_half: f64,
}

impl Camera {
    pub fn orbit(target: [f64; 3], az_deg: f64, el_deg: f64, dist: f64) -> Camera {
        let (az, el) = (az_deg.to_radians(), el_deg.to_radians());
        let eye = [
            target[0] - dist * el.cos() * az.cos(),
            target[1] - dist * el.cos() * az.sin(),
            target[2] + dist * el.sin(),
        ];
        let n = |v: [f64; 3]| {
            let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt().max(1e-12);
            [v[0] / l, v[1] / l, v[2] / l]
        };
        let cross = |a: [f64; 3], b: [f64; 3]| [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
        let fwd = n([target[0] - eye[0], target[1] - eye[1], target[2] - eye[2]]);
        let right = n(cross(fwd, [0.0, 0.0, 1.0]));
        let up = cross(right, fwd);
        Camera {
            eye,
            fwd,
            right,
            up,
            tan_half: (20f64).to_radians().tan(),
        }
    }

    /// Projection d'un point (mailles) en coordonnées normalisées [−1, 1]
    /// (y vers le haut) pour un rapport largeur/hauteur donné ; None derrière l'œil.
    pub fn project(&self, p: [f64; 3], aspect: f64) -> Option<[f64; 2]> {
        let d = [p[0] - self.eye[0], p[1] - self.eye[1], p[2] - self.eye[2]];
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let z = dot(d, self.fwd);
        if z <= 1e-6 {
            return None;
        }
        Some([dot(d, self.right) / (z * self.tan_half * aspect), dot(d, self.up) / (z * self.tan_half)])
    }
}

/// Réglages d'une image 3D.
#[derive(Clone, Copy, Debug)]
pub struct RenderOptions {
    pub camera: Camera,
    pub width: usize,
    pub height: usize,
    /// Plan de coupe : axe (0 : x, 1 : y, 2 : z) et position (maille).
    pub plane_axis: usize,
    pub plane_pos: f64,
    /// Coupe en vorticité (sinon vitesse).
    pub vorticity: bool,
    /// Isosurfaces Q > iso_level·(U/L)².
    pub iso: bool,
    pub iso_level: f64,
    /// Longueur de référence (mailles).
    pub l: f64,
    pub background: [f32; 3],
    pub body: [f32; 3],
    /// Maquette : forme lisse, avec la grille, ou mailles du calcul.
    pub body_view: BodyView,
}

/// Représentation de la maquette dans la vue 3D.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[repr(u8)]
pub enum BodyView {
    /// Surface exacte de la forme.
    #[default]
    Smooth = 0,
    /// Surface exacte, grille du calcul tracée dessus.
    Grid = 1,
    /// Ce que voit le calcul : mailles solides, et mailles à liens de paroi mince.
    Cells = 2,
}

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
pub(crate) struct Params {
    nx: u32,
    ny: u32,
    nz: u32,
    n: u32,
    tau0: f32,
    tau_sponge: f32,
    sponge_start: u32,
    u_in: f32,
    accel: f32,
    c18: f32,
    write_macro: u32,
    les: u32,
    rows: u32,
    sa: u32,
    nu: f32,
    /// Longueur d'amortissement pariétal de la sous-maille (mailles).
    les_a: f32,
    /// Transition imposée à cette fraction de la corde locale (0 : turbulent).
    xtr: f32,
    /// Épaisseur (mailles) de la zone sans production autour d'une paroi laminaire.
    band: f32,
    /// Variantes de diagnostic (bits, voir `Gpu3D::set_essai`).
    essai: u32,
    /// Plafond de viscosité turbulente des diagnostics 1024 et 2048.
    essai_nu: f32,
    /// Stockage creux : blocs selon x, y, z (0 : grille boîte) et début des
    /// coordonnées des blocs dans la table.
    bnx: u32,
    bny: u32,
    bnz: u32,
    bcoord: u32,
}

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct SliceParams {
    nx: u32,
    ny: u32,
    nz: u32,
    axis: u32,
    pos: u32,
    field: u32,
    w: u32,
    h: u32,
    u: f32,
    l: f32,
    pad: [f32; 2],
}

/// Pas par soumission (paramètres en décalage dynamique).
pub(crate) const BATCH: usize = 64;
pub(crate) const SLOT: u64 = 256;
/// Soumissions de pas en vol au plus (`Gpu3D::run`).
const EN_VOL_MAX: usize = 8;

/// Grandeur affichée sur une coupe.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SliceField {
    /// |u|/U
    Speed,
    /// Composante de ω normale au plan, ·L/U
    Vorticity,
    /// Critère Q·(L/U)²
    Q,
}

/// Coupe : `w × h` valeurs, ligne j = coordonnée croissante ; 1e30 dans le corps.
#[derive(Clone, Debug, Default)]
pub struct Slice {
    pub w: usize,
    pub h: usize,
    pub data: Vec<f32>,
}

/// Efforts moyens depuis la relecture précédente (unités réseau).
#[derive(Clone, Copy, Debug, Default)]
pub struct Force3 {
    pub fx: f64,
    pub fy: f64,
    pub fz: f64,
    pub steps: u64,
}

/// Octets par maille sur le GPU (populations ping-pong, table pariétale, champ).
pub fn bytes_per_cell(fp16: bool) -> usize {
    (if fp16 { 2 * 10 * 4 } else { 2 * 19 * 4 }) + 4 + 16
}

/// Modèle de turbulence RANS (avec loi de paroi), ou aucun.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Turbulence {
    /// Sans modèle ni loi de paroi (sous-maille seule).
    Laminar,
    /// Spalart–Allmaras (défaut).
    Sa,
    /// k-ε RNG (Yakhot et al. 1992), voir `KE3_WGSL`.
    KEpsRng,
    /// VLES sur la base k-ε RNG : ν_t = C_μ k²/ε réduit par le filtre de
    /// maille de Johansen et al. (2004), min(1, Δ ε/k^{3/2}), hors de la
    /// couche limite attachée (protection f_d de la DDES, Spalart et al.
    /// 2006). Voir `KE_FNS` et la section « VLES » du README.
    Vles,
}

impl Turbulence {
    /// Modèle RANS des bancs : `SOUFFLERIE_RANS=keps` choisit k-ε RNG,
    /// `SOUFFLERIE_RANS=vles` sa réduction VLES,
    /// sinon Spalart–Allmaras.
    pub fn rans_from_env() -> Turbulence {
        match std::env::var("SOUFFLERIE_RANS").as_deref() {
            Ok("keps") => Turbulence::KEpsRng,
            Ok("vles") => Turbulence::Vles,
            _ => Turbulence::Sa,
        }
    }
}

/// `true` : Spalart–Allmaras, `false` : aucun modèle (appels existants).
impl From<bool> for Turbulence {
    fn from(sa: bool) -> Turbulence {
        if sa { Turbulence::Sa } else { Turbulence::Laminar }
    }
}

pub struct Gpu3D {
    pub(crate) device: wgpu::Device,
    pub(crate) queue: wgpu::Queue,
    adapter_name: String,
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub n: usize,
    pub fp16: bool,
    /// Rangement des populations (voir `POP16`).
    pub(crate) p16: bool,
    pub u_in: f32,
    pub tau0: f64,
    tau_sponge: f64,
    pub les: bool,
    sponge_start: u32,
    cs: f64,
    pub(crate) pops: [[wgpu::Buffer; 2]; 2],
    pub(crate) wall: wgpu::Buffer,
    /// Distance signée à la maquette (affichage) et sa grille (origine, pas ; taille).
    sdf: wgpu::Buffer,
    sdf_geom: ([f32; 4], [f32; 4]),
    /// Type des mailles pour l'affichage (2 bits par maille).
    vox: wgpu::Buffer,
    _qbuf: wgpu::Buffer,
    /// Spalart–Allmaras avec loi de paroi : ν̃ (deux tampons alternés),
    /// distance à la paroi, géométrie pariétale et u_τ.
    pub sa: bool,
    /// Modèle k-ε RNG au lieu de Spalart–Allmaras (`sa` vrai aussi) :
    /// tampons `nut` à deux valeurs (k, ε) par maille.
    pub keps: bool,
    /// VLES (`keps` vrai aussi) : ν_t réduit, rangé après les paires (k, ε)
    /// (dernier tiers des tampons `nut`), écrit par la passe `ke3d`.
    pub vles: bool,
    pub(crate) nu: f64,
    /// Transition imposée (fraction de la corde locale ; 0 : turbulent partout).
    xtr: f64,
    _slab: wgpu::Buffer,
    /// Mailles par longueur de référence (Reynolds de maille, sous-maille).
    l: usize,
    les_a: f64,
    pub(crate) nut: [wgpu::Buffer; 2],
    _dist: wgpu::Buffer,
    _wall_geo: wgpu::Buffer,
    _utau: wgpu::Buffer,
    pub(crate) sa_pipe: wgpu::ComputePipeline,
    pub(crate) sa_bind: Vec<wgpu::BindGroup>,
    /// Loi de paroi : passe sur les mailles pariétales (None sans SA).
    pub(crate) wall_pass: Option<WallPass>,
    /// ν̃ initial (3ν dans le fluide, 0 dans le corps), pour `reset`.
    nut0: Vec<f32>,
    force: wgpu::Buffer,
    force_read: wgpu::Buffer,
    pub(crate) params: wgpu::Buffer,
    pub(crate) macro_buf: wgpu::Buffer,
    pub(crate) step_pipe: wgpu::ComputePipeline,
    /// Collision des mailles pariétales (None sans paroi), et ses groupes.
    pub(crate) step_wall_pipe: Option<wgpu::ComputePipeline>,
    /// Disposition des pipelines de collision, et (sous-maille, bords,
    /// diagnostics) pour lesquels ils sont spécialisés.
    step_layout: wgpu::PipelineLayout,
    spec: [u32; 3],
    pub(crate) wall_disp: (u32, u32),
    pub(crate) disp_step: (u32, u32),
    /// Grille creuse : mailles du cœur de la collision (liste), 0 : toutes.
    n_active: u32,
    _step_list: wgpu::Buffer,
    /// Mailles fluides ordinaires, un bit par maille.
    plain: wgpu::Buffer,
    _wcell: wgpu::Buffer,
    init_pipe: wgpu::ComputePipeline,
    slice_pipe: wgpu::ComputePipeline,
    derive_pipe: wgpu::ComputePipeline,
    render_pipe: wgpu::ComputePipeline,
    /// ω et Q par maille (alloué au premier rendu 3D).
    deriv: Option<wgpu::Buffer>,
    /// Tampons du rendu réutilisés d'une image à l'autre : taille, image, relecture.
    render_bufs: Option<(usize, usize, wgpu::Buffer, wgpu::Buffer)>,
    pub(crate) bind: Vec<wgpu::BindGroup>,
    init_bind: Vec<wgpu::BindGroup>,
    pub(crate) disp: (u32, u32),
    /// Lancement du transport de ν̃ (k-ε).
    pub(crate) sa_disp: (u32, u32),
    pub(crate) cur: usize,
    pub time: u64,
    pub(crate) force_steps: u64,
    pub wall_count: usize,
    pub solid_count: usize,
    /// Mailles solides (affichage, pochoirs des raccords).
    pub(crate) solid: Vec<bool>,
    /// Couche creuse : mailles fantômes (non calculées ; vide pour une
    /// boîte). Les relevés sur la surface les évitent.
    pub(crate) ghost: Vec<bool>,
    /// Le GPU ne répond plus (lecture impossible, pilote réinitialisé).
    pub lost: bool,
    /// Grille fine imbriquée : enveloppe de mailles fantômes, ni entrée, ni
    /// sortie, ni bandes absorbantes.
    pub(crate) nested: bool,
    /// Sol mobile en z = −½ et plafond glissant (au lieu de bords périodiques en z).
    ground: bool,
    /// Maquette d'envergure périodique (plaque plane, profil infini) : pas
    /// de bandes absorbantes sur les bords latéraux y, que la maquette
    /// traverse (elles freineraient la couche limite de ses bords).
    pub periodic_span: bool,
    /// Veine périodique en x, sans entrée ni sortie (cas de validation).
    periodic_x: bool,
    /// Variantes de diagnostic de la loi de paroi (bits, `set_essai`).
    essai: u32,
    essai_nu: f32,
    /// Grilles fines autour de la maquette (raffinement local).
    pub(crate) refine: Option<Box<crate::refine3d::Refinement>>,
    /// Maille de chaque maille pariétale (dans l'ordre de leur numéro).
    wall_cells: Vec<u32>,
    /// Couche pariétale fine (grille partielle qui épouse la surface) :
    /// mailles fantômes à remplir depuis la grille parente (voisines d'une
    /// maille active), et mailles actives loin de toute maille fantôme (leur
    /// état est recopié vers la grille parente). None : grille boîte.
    pub(crate) band_fill: Option<Vec<u32>>,
    pub(crate) band_deep: Option<Vec<bool>>,
    /// Positions (indices de la boîte) des fantômes à remplir.
    pub(crate) band_fill_pos: Option<Vec<u32>>,
    /// Stockage creux par blocs (couche pariétale), et sa table pour les shaders.
    pub(crate) sparse: Option<std::sync::Arc<SparseMap>>,
    pub(crate) blk_buf: wgpu::Buffer,
    /// Plus petit cosinus entre les normales des faces touchées par les liens
    /// de chaque maille pariétale (arête vive : < 0,5).
    edge_cos: Vec<f32>,
    /// Part de la pression de référence (ρ₀ c_s²) dans l'effort de chaque
    /// maille pariétale : Σ −2 w_i c_i sur ses liens. Nulle en somme sur une
    /// maquette fermée, mais pas maille par maille ; retirée des efforts par
    /// maille pour n'y garder que la pression relative et le frottement.
    rest_force: Vec<[f32; 3]>,
    /// Maillage de la maquette dans les coordonnées de cette grille (efforts
    /// par intégration sur la surface).
    surface: Option<std::sync::Arc<crate::stl::MeshBody>>,
    /// Décambrage visqueux par transpiration (voir `enable_transpiration`).
    transp: Option<Box<Transpiration>>,
    /// Loi de paroi hors équilibre (TBLE, `set_paroi_tble`) et son état
    /// (profils 1D sous chaque nœud pariétal) ; None : Spalding (défaut).
    tble: Option<(crate::tble::ParoiTble, wgpu::Buffer)>,
    /// Paroi volumique par surfels (`Body3D::surfels`).
    pub(crate) surf: Option<Box<SurfPass>>,
    /// Accélération uniforme constante selon x (`set_accel`).
    accel_x: f32,
}

/// Passe des surfels (rassemblement, réflexion ; voir `surfel3d`).
pub(crate) struct SurfPass {
    pipe: wgpu::ComputePipeline,
    bind: Vec<wgpu::BindGroup>,
    disp: (u32, u32),
    /// Nombre de surfels, et leurs centres (efforts par position).
    n: usize,
    center: Vec<[f32; 3]>,
    /// Volume fluide de chaque maille réceptrice (bilan de masse).
    vol: Vec<f32>,
    stats: crate::surfel3d::Stats,
    gpu_bytes: usize,
    _bufs: [wgpu::Buffer; 5],
}

/// Passe de la loi de paroi (voir `WALL3_WGSL`).
pub(crate) struct WallPass {
    pipe: wgpu::ComputePipeline,
    bind: wgpu::BindGroup,
    disp: (u32, u32),
}

/// Masque des mailles fluides ordinaires (type 0 : ni paroi, ni solide, ni
/// fantôme, ni recouverte), un bit par maille : la collision des grilles
/// boîtes et le transport de ν̃ ne lisent le type complet (4 octets) que
/// pour les autres (ν̃ −10 %).
fn plain_mask(wall: &[u32]) -> Vec<u32> {
    let mut m = vec![0u32; wall.len().div_ceil(32).max(1)];
    for (i, &w) in wall.iter().enumerate() {
        if w == 0 {
            m[i >> 5] |= 1 << (i & 31);
        }
    }
    m
}

/// Bit de `set_essai` : correction du défaut cubique du D3Q19 (moment
/// d'équilibre ρ u_α³ manquant ; README, « Domaine de Mach »).
pub const MACH3_BIT: u32 = 1 << 28;
/// Bornes du code de cette correction dans la collision : retiré sans le
/// bit (code par défaut inchangé à l'octet près).
const M3_OPEN: &str = "/*M3[*/";
const M3_CLOSE: &str = "/*]M3*/";

fn mach3_strip(code: &str, keep: bool) -> String {
    if keep {
        return code.to_string();
    }
    let mut out = String::with_capacity(code.len());
    let mut rest = code;
    while let Some(a) = rest.find(M3_OPEN) {
        out.push_str(&rest[..a]);
        let b = rest[a..].find(M3_CLOSE).expect("borne M3 non fermée") + a;
        rest = &rest[b + M3_CLOSE.len()..];
    }
    out.push_str(rest);
    out
}

/// Code de la collision où les paramètres fixes d'une grille deviennent des
/// constantes (SA, sous-maille, type de bords, diagnostics, stockage) : le
/// compilateur retire leurs branches (collision 4 % plus rapide).
fn specialize(code: &str, c: [u32; 7]) -> String {
    let code = mach3_strip(code, c[3] & MACH3_BIT != 0);
    let body = code
        .replace("P.lesA", "P.\u{1}A")
        .replace("P.essaiNu", "P.\u{2}Nu")
        .replace("P.sa", "SA_C")
        .replace("P.les", "LES_C")
        .replace("P.rows", "ROWS_C")
        .replace("P.essai", "ESSAI_C")
        .replace("P.bnx", "BNX_C")
        .replace("P.\u{1}A", "P.lesA")
        .replace("P.\u{2}Nu", "P.essaiNu");
    // Directives `enable` en tête du module.
    let (en, body) = match body.strip_prefix("enable wgpu_int16;\n") {
        Some(b) => ("enable wgpu_int16;\n", b.to_string()),
        None => ("", body),
    };
    format!(
        "{en}const SA_C = {}u;\nconst LES_C = {}u;\nconst ROWS_C = {}u;\nconst ESSAI_C = {}u;\nconst BNX_C = {}u;\nconst WALLN_C = {}u;\nconst ACTN_C = {}u;\n{body}",
        c[0], c[1], c[2], c[3], c[4], c[5], c[6]
    )
}

/// Pipelines de la collision (cœur, et mailles pariétales s'il y en a) pour
/// les constantes `c` (voir `specialize`).
fn step_pipelines(
    device: &wgpu::Device,
    layout: &wgpu::PipelineLayout,
    fp16: bool,
    p16: bool,
    walls: bool,
    surf: bool,
    c: [u32; 7],
) -> (wgpu::ComputePipeline, Option<wgpu::ComputePipeline>) {
    let mk = |label: &str, code: String| {
        let module = compute_module_raw(device, label, specialize(&wgsl_hook(label, code), c));
        device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some(label),
            layout: Some(layout),
            module: &module,
            entry_point: Some("main"),
            compilation_options: Default::default(),
            cache: None,
        })
    };
    let bulk = mk("d3q19", wgsl_step(fp16, p16, false, surf).replace("WGu", &format!("{}u", step_wg())));
    let wall = walls.then(|| mk("d3q19_paroi", wgsl_step(fp16, p16, true, surf)));
    (bulk, wall)
}

impl Gpu3D {
    /// Pipelines de collision à jour des paramètres (sous-maille après
    /// `set_viscosity`, bords, diagnostics `set_essai`) : recompilés s'ils
    /// ont changé depuis.
    pub(crate) fn ensure_spec(&mut self) {
        let p = self.params(0.0, 0.0, false);
        let k = [p.les, p.rows, p.essai];
        if k == self.spec {
            return;
        }
        let (a, b) = step_pipelines(
            &self.device,
            &self.step_layout,
            self.fp16,
            self.p16,
            self.step_wall_pipe.is_some(),
            self.surf.is_some(),
            [p.sa, k[0], k[1], k[2], p.bnx, self.wall_count as u32, self.n_active],
        );
        self.step_pipe = a;
        self.step_wall_pipe = b;
        self.spec = k;
    }

    /// Un pas de collision (cœur, puis mailles pariétales), état `cur`,
    /// paramètres du créneau `slot`.
    pub(crate) fn encode_step(&self, pass: &mut wgpu::ComputePass, slot: u32) {
        // Surfels : rassemblement et réflexion sur l'état d'avant le pas.
        if let Some(sp) = &self.surf {
            pass.set_pipeline(&sp.pipe);
            pass.set_bind_group(0, &sp.bind[self.cur], &[slot]);
            pass.dispatch_workgroups(sp.disp.0, sp.disp.1, 1);
        }
        pass.set_pipeline(&self.step_pipe);
        pass.set_bind_group(0, &self.bind[self.cur], &[slot]);
        pass.dispatch_workgroups(self.disp_step.0, self.disp_step.1, 1);
        if let Some(p) = &self.step_wall_pipe {
            pass.set_pipeline(p);
            pass.dispatch_workgroups(self.wall_disp.0, self.wall_disp.1, 1);
        }
    }
}

impl WallPass {
    pub(crate) fn encode(&self, pass: &mut wgpu::ComputePass, slot: u32) {
        pass.set_pipeline(&self.pipe);
        pass.set_bind_group(0, &self.bind, &[slot]);
        pass.dispatch_workgroups(self.disp.0, self.disp.1, 1);
    }
}

/// Couche limite intégrale sur le maillage de la grille entière la plus
/// fine, et sa vitesse de transpiration par maille pariétale.
struct Transpiration {
    topo: crate::ibl3d::Topology,
    /// Par grille (None : la veine sans raffinement, Some(k) : niveau k de
    /// `refine`) : triangle de chaque maille pariétale et v_n imposée.
    maps: Vec<(Option<usize>, Vec<u32>, Vec<f32>)>,
    xtr: f64,
    /// Vitesse extérieure prise au bord de la couche déficitaire que le
    /// calcul résout déjà, et seul le déplacement qui lui manque imposé
    /// (δ* − δ*_calcul) ; sinon u_e à WM_H et δ* entier.
    resolved: bool,
}

/// Bilan d'un couplage de transpiration (unités de la grille la plus fine).
#[derive(Clone, Debug, Default)]
pub struct TranspirationReport {
    pub layer: crate::ibl3d::Layer3,
    /// Vitesse extérieure par triangle.
    pub ue: Vec<[f64; 3]>,
    /// Aires et centres des triangles (mailles de cette grille).
    pub area: Vec<f64>,
    pub centers: Vec<[f64; 3]>,
    /// Normales sortantes des triangles.
    pub normals: Vec<[f64; 3]>,
    /// Épaisseur de déplacement que le champ du calcul contient déjà,
    /// ∫₀^{h_e} (1 − u_t/u_e) dh, et hauteur h_e du bord de cette couche
    /// (premier maximum de u_t au-dessus de WM_H), par triangle.
    pub dstar_res: Vec<f64>,
    pub h_edge: Vec<f64>,
    /// v_n imposée par triangle (avant sous-relaxation).
    pub vn: Vec<f64>,
    /// Part de la surface où u_e n'a pu être relevée.
    pub missed: f64,
    /// v_n/u_∞ moyenne (aire) des parties soufflantes et maximale.
    pub vn_mean: f64,
    pub vn_max: f64,
    /// Part de la surface décollée (H plafonné).
    pub separated: f64,
    /// Débit de déplacement qui part dans le sillage, Σ v_n dA / (u_∞ · L²).
    pub outflow: f64,
    /// Mailles par longueur de référence sur cette grille.
    pub l: f64,
}

/// Efforts par maille pariétale, moyennés par pas, en unités de la veine :
/// position (mailles de la veine) et force.
/// Efforts intégrés sur la surface du maillage (unités de la veine) :
/// pression à la paroi, frottement de la loi de paroi (Spalding) sur la
/// vitesse tangentielle relevée à une maille de la paroi.
#[derive(Clone, Debug, Default)]
pub struct SurfaceForces {
    /// Pression relevée à une maille de la paroi (sans extrapolation).
    pub pressure: [f64; 3],
    /// Pression extrapolée linéairement jusqu'à la paroi (relevés à h et h + 1).
    pub pressure_wall: [f64; 3],
    pub friction: [f64; 3],
    /// Part de la surface où le relevé a échoué (point d'échantillonnage
    /// entouré de mailles solides).
    pub missed: f64,
    /// Par triangle : centre (mailles de la veine) et effort (pression à la
    /// paroi et frottement), pour le moment de tangage et la répartition en
    /// envergure.
    pub loads: Vec<([f64; 3], [f64; 3])>,
    /// Par triangle (même ordre que `loads`) : normale sortante et pression
    /// à la paroi (unités du réseau, p = (ρ − 1)/3).
    pub walls: Vec<([f64; 3], f64)>,
    /// Par triangle (même ordre que `loads`) : ordonnées y de ses trois
    /// sommets, rangées (mailles de la veine), pour répartir son effort sur
    /// son étendue en envergure (`span_distribution`).
    pub spans: Vec<[f64; 3]>,
    /// Par triangle (même ordre que `loads`) : son rang dans le maillage,
    /// la contrainte pariétale moyenne de la loi de paroi (vecteur, unités
    /// du réseau de la veine, ρu_τ² le long de la vitesse relevée) et la
    /// hauteur de raccord de la loi de paroi de la grille qui le relève
    /// (mailles de la veine).
    pub tri: Vec<usize>,
    pub shear: Vec<[f64; 3]>,
    pub wm_h: Vec<f64>,
}

/// Répartition en envergure d'efforts portés par des triangles : chaque
/// effort est réparti sur l'étendue en y de son triangle selon la densité
/// d'aire du triangle le long de y (en chapeau : croissante du sommet le plus
/// bas au sommet médian, décroissante ensuite), et non posé en son centre.
/// Posé au centre, un maillage à stations d'envergure régulières met une ou
/// deux rangées de triangles dans chaque tranche selon où tombent leurs
/// centres : creux et bosses parasites (aile AR 4 à stations de 0,087 m,
/// tranches de 1/32 m : aire par tranche variant de 1 à 1,5). Tranches
/// [y0 + k dy, y0 + (k + 1) dy), k < nb ; rend la somme par tranche.
pub fn span_distribution(spans: &[[f64; 3]], values: &[[f64; 3]], y0: f64, dy: f64, nb: usize) -> Vec<[f64; 3]> {
    let mut out = vec![[0.0f64; 3]; nb];
    // Part de l'aire du triangle (sommets en y : a ≤ b ≤ c) sous l'ordonnée y.
    let cdf = |t: &[f64; 3], y: f64| -> f64 {
        let [a, b, c] = *t;
        if y <= a {
            0.0
        } else if y >= c {
            1.0
        } else if y <= b {
            (y - a) * (y - a) / ((b - a) * (c - a))
        } else {
            1.0 - (c - y) * (c - y) / ((c - b) * (c - a))
        }
    };
    let bin = |y: f64| (((y - y0) / dy).floor().max(0.0) as usize).min(nb - 1);
    for (t, v) in spans.iter().zip(values) {
        if t[2] - t[0] < 1e-9 * dy {
            let k = bin(t[0]);
            for a in 0..3 {
                out[k][a] += v[a];
            }
            continue;
        }
        let (k0, k1) = (bin(t[0]), bin(t[2]));
        for k in k0..=k1 {
            let lo = y0 + k as f64 * dy;
            let f = cdf(t, lo + dy) - cdf(t, lo);
            for a in 0..3 {
                out[k][a] += f * v[a];
            }
        }
    }
    out
}

impl SurfaceForces {
    pub fn total(&self) -> [f64; 3] {
        [0, 1, 2].map(|a| self.pressure_wall[a] + self.friction[a])
    }
}

/// Traînée en champ lointain (voir `Gpu3D::wake_drag`), unités du réseau.
#[derive(Clone, Debug, Default)]
pub struct WakeDrag {
    /// Perte de pression totale dans le sillage visqueux (traînée de profil).
    pub loss_wake: f64,
    /// Perte de pression totale hors du sillage (traînée parasite numérique).
    pub loss_outside: f64,
    /// Énergie cinétique transversale (traînée induite).
    pub induced: f64,
    /// Débit aval / débit amont.
    pub mass_ratio: f64,
    /// Aire du sillage visqueux dans le plan aval (mailles de la veine).
    pub wake_area: f64,
}

/// Hauteurs (mailles) des deux relevés de pression extrapolés jusqu'à la paroi.
const PW_H1: f64 = 0.5;
const PW_H2: f64 = 1.0;

/// Boîtes des grilles des bords (mailles de la grille parente) : (coin bas,
/// coin haut, domaine propre en mailles de la grille parente). Bord
/// d'attaque : sommets des faces tournées vers l'amont (n_x < −0,5). Bord de
/// fuite : arêtes vives saillantes (`geom3d::feature_edges`, seuil de la loi
/// de paroi) orientées en envergure (|d_y| > 0,5 |d|), donc pas les contours
/// de saumon. Points rangés en y et regroupés en tronçons d'étendue en x et
/// z au plus w = max(0,15 L, 4 e) : une aile en flèche est couverte par une
/// suite de boîtes au lieu d'une seule. Marges : e = MARGIN + 4 mailles
/// partout, plus 10 % de L vers l'aval au bord d'attaque (pic de dépression
/// et début de la recompression) et 5 % de L de part et d'autre au bord de
/// fuite (fin de la recompression, départ du sillage). Domaines propres :
/// les tronçons voisins se partagent l'envergure au milieu de leur écart.
#[allow(clippy::type_complexity)]
fn edge_boxes(tris: &[crate::stl::Tri], l: f64, dims: [usize; 3], le: bool, te: bool) -> Result<Vec<([usize; 3], [usize; 3], [[f64; 3]; 2])>, String> {
    let e = (crate::refine3d::MARGIN3 + 4) as f64;
    let w = (0.15 * l).max(4.0 * e);
    let mut sets: Vec<(Vec<[f64; 3]>, [f64; 2])> = Vec::new();
    let nrm = |t: &crate::stl::Tri| {
        let e1 = [0, 1, 2].map(|a| t[1][a] - t[0][a]);
        let e2 = [0, 1, 2].map(|a| t[2][a] - t[0][a]);
        [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]]
    };
    if le {
        let vol: f64 = tris
            .iter()
            .map(|t| {
                let n = nrm(t);
                (0..3).map(|a| t[0][a] * n[a]).sum::<f64>()
            })
            .sum();
        let sgn = if vol >= 0.0 { 1.0 } else { -1.0 };
        let mut pts = Vec::new();
        for t in tris {
            let n = nrm(t);
            let m = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            if m > 1e-12 && sgn * n[0] / m < -0.5 {
                pts.extend_from_slice(t);
            }
        }
        if pts.is_empty() {
            return Err("raffinement du bord d'attaque : aucune face tournée vers l'amont".into());
        }
        // Marges en x : [amont, aval].
        sets.push((pts, [e, e.max(0.1 * l)]));
    }
    if te {
        let ix = crate::geom3d::Indexed::from_tris(tris);
        let mut pts = Vec::new();
        for f in crate::geom3d::feature_edges(&ix, crate::geom3d::FEATURE_DEG) {
            let d = [0, 1, 2].map(|a| f.b[a] - f.a[a]);
            let dl = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            if f.convex && dl > 0.0 && d[1].abs() > 0.5 * dl {
                // Points toutes les demi-mailles le long de l'arête.
                let n = (dl / 0.5).ceil().max(1.0) as usize;
                for q in 0..=n {
                    let t = q as f64 / n as f64;
                    pts.push([0, 1, 2].map(|a| f.a[a] + t * d[a]));
                }
            }
        }
        if pts.is_empty() {
            return Err("raffinement du bord de fuite : aucune arête vive en envergure".into());
        }
        let m = e.max(0.05 * l);
        sets.push((pts, [m, m]));
    }
    let band = (crate::refine3d::MARGIN3 + 2) as f64;
    let mut out = Vec::new();
    if std::env::var("SOUFFLERIE_BOITE_ENTIERE").is_ok() {
        // Diagnostic : une boîte qui entoure toute la maquette.
        let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
        for t in tris {
            for v in t {
                for q in 0..3 {
                    lo[q] = lo[q].min(v[q]);
                    hi[q] = hi[q].max(v[q]);
                }
            }
        }
        let blo = [0, 1, 2].map(|a| (lo[a] - e - band).floor().max(3.0) as usize);
        let bhi = [0, 1, 2].map(|a| ((hi[a] + e + band).ceil() as usize).min(dims[a] - 4));
        let core = [[0, 1, 2].map(|a| blo[a] as f64 + band), [0, 1, 2].map(|a| bhi[a] as f64 - band)];
        return Ok(vec![(blo, bhi, core)]);
    }
    for (mut pts, mx) in sets {
        pts.sort_by(|a, b| a[1].total_cmp(&b[1]));
        // Tronçons d'étendue en x et z au plus w.
        let mut chunks: Vec<(usize, usize)> = Vec::new();
        let mut start = 0;
        let (mut x0, mut x1, mut z0, mut z1) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
        for (k, p) in pts.iter().enumerate() {
            let (nx0, nx1, nz0, nz1) = (x0.min(p[0]), x1.max(p[0]), z0.min(p[2]), z1.max(p[2]));
            if k > start && (nx1 - nx0 > w || nz1 - nz0 > w) {
                chunks.push((start, k));
                start = k;
                (x0, x1, z0, z1) = (p[0], p[0], p[2], p[2]);
            } else {
                (x0, x1, z0, z1) = (nx0, nx1, nz0, nz1);
            }
        }
        chunks.push((start, pts.len()));
        let nc = chunks.len();
        for (c, &(a, b)) in chunks.iter().enumerate() {
            let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
            for p in &pts[a..b] {
                for q in 0..3 {
                    lo[q] = lo[q].min(p[q]);
                    hi[q] = hi[q].max(p[q]);
                }
            }
            // Domaine propre en y : jusqu'au milieu de l'écart avec les voisins.
            let ylo = if c == 0 { lo[1] - e } else { 0.5 * (pts[a - 1][1] + lo[1]) };
            let yhi = if c + 1 == nc { hi[1] + e } else { 0.5 * (hi[1] + pts[b][1]) };
            let blo = [
                (lo[0] - mx[0] - band).floor().max(3.0) as usize,
                (ylo - band - e).floor().max(3.0) as usize,
                (lo[2] - e - band).floor().max(3.0) as usize,
            ];
            let bhi = [
                ((hi[0] + mx[1] + band).ceil() as usize).min(dims[0] - 4),
                ((yhi + band + e).ceil() as usize).min(dims[1] - 4),
                ((hi[2] + e + band).ceil() as usize).min(dims[2] - 4),
            ];
            // Domaine propre : la boîte moins sa bande de bord ; en y, la part
            // d'envergure du tronçon.
            let core = [
                [blo[0] as f64 + band, ylo.max(blo[1] as f64 + band), blo[2] as f64 + band],
                [bhi[0] as f64 - band, yhi.min(bhi[1] as f64 - band), bhi[2] as f64 - band],
            ];
            out.push((blo, bhi, core));
        }
    }
    Ok(out)
}

/// Épaisseurs (mailles de chaque couche) de couches pariétales imbriquées
/// rendues compatibles : une couche fine d'épaisseur D est remplie, à D
/// mailles fines (D/2 parentes) de la paroi, depuis sa couche parente, avec
/// un pochoir de ±2 mailles parentes et la bande de raccord ; la couche
/// parente doit donc s'étendre au-delà : D_parent ≥ D/2 + 2 MARGIN + 3.
/// Sinon les fantômes de la couche fine sont remplis depuis des mailles que
/// la couche parente ne calcule pas (écoulement faux).
/// `--bande auto` : assez de couches pariétales pour mettre au moins 2
/// mailles dans le plus petit rayon de bord d'attaque des coupes (hors
/// saumons), au plus 4, au-dessus de `depth` grilles boîtes. Rend le nombre
/// de couches, ce rayon (mailles de la veine) et le plus grand coin entre
/// facettes au bord d'attaque (au-delà de 60°, le rayon mesuré est celui des
/// facettes). (0, ∞, 0) sans maillage.
pub fn auto_band_count(body: &Body3D, depth: usize) -> (usize, f64, f64) {
    let secs = body.mesh.as_ref().map(|m| crate::geom3d::sections(&m.tris, 20)).unwrap_or_default();
    let r = secs.iter().map(|s| s.le_radius).filter(|r| r.is_finite() && *r > 0.0).fold(f64::INFINITY, f64::min);
    let corner = secs.iter().map(|s| s.le_corner_deg).fold(0.0f64, f64::max);
    let mut n = 0;
    while r.is_finite() && n < 4 && r * ((1u64 << (depth + n)) as f64) < 2.0 {
        n += 1;
    }
    (n, r, corner)
}

/// Région d'une couche limitée aux bords (attaque, fuite, saumon) : mailles
/// à moins de `radius` (mailles de la couche) des arêtes de `set`.
#[derive(Clone, Debug)]
pub struct EdgeRegion {
    pub set: crate::geom3d::FeatureSet,
    pub radius: f64,
}

/// Région d'une couche limitée à une partie de la surface (raffinement de
/// paroi par zones, `raff_paroi`) : triangles retenus du maillage de la
/// maquette (même ordre que `Body3D::mesh`, que `Body3D::refined` garde).
/// Mailles actives : à moins de `band` mailles de ces triangles.
#[derive(Clone, Debug)]
pub struct ZoneRegion {
    pub tris: Vec<bool>,
}

/// Couche pariétale (grille creuse) : mailles actives à moins de `band`
/// mailles (de cette couche) de la paroi ; plus, si `wake` > 0, celles
/// balayées vers l'aval (+x) sur `wake` mailles depuis la dernière maille
/// de la couche de chaque rangée ; si `edges`, seulement près des arêtes ;
/// si `zone`, seulement près des triangles retenus.
#[derive(Clone, Debug, Default)]
pub struct LayerSpec {
    pub band: f64,
    pub wake: f64,
    pub edges: Option<std::sync::Arc<EdgeRegion>>,
    pub zone: Option<std::sync::Arc<ZoneRegion>>,
}

/// Plan de raffinement : grilles boîtes (emprise amont, aval, de côté en
/// longueurs de référence autour de la maquette, comme `EXTENT3`), puis
/// couches pariétales, ou grilles des bords (`le`, `te`, sans couches).
#[derive(Clone, Debug, Default)]
pub struct RefinePlan {
    pub boxes: Vec<(f64, f64, f64)>,
    pub layers: Vec<LayerSpec>,
    pub le: bool,
    pub te: bool,
}

pub fn nested_bands(bands: &[f64]) -> Vec<f64> {
    let mut b = bands.to_vec();
    let m = (2 * crate::refine3d::MARGIN3 + 3) as f64;
    for i in (0..b.len().saturating_sub(1)).rev() {
        b[i] = b[i].max((0.5 * b[i + 1] + m).ceil());
    }
    b
}

/// Côté des blocs du stockage creux (mailles).
pub(crate) const SB: usize = 8;
const SB3: usize = SB * SB * SB;

/// Stockage creux par blocs d'une couche pariétale : la maille (x, y, z)
/// du bloc de rang `slot` est rangée en slot · 512 + local. Le rang 0 est un
/// bloc vide (mailles fantômes) qui reçoit toute lecture hors des blocs.
#[derive(Clone, Debug)]
pub(crate) struct SparseMap {
    /// Nombre de blocs selon x, y, z (boîte entière).
    bd: [usize; 3],
    /// Rang de chaque bloc de la boîte (0 : absent).
    table: Vec<u32>,
    /// Coordonnées (en blocs) de chaque rang.
    coords: Vec<[u32; 3]>,
}

impl SparseMap {
    fn new(bd: [usize; 3], used: &[bool]) -> SparseMap {
        let mut table = vec![0u32; used.len()];
        let mut coords = vec![[u32::MAX; 3]];
        for (b, &u) in used.iter().enumerate() {
            if u {
                table[b] = coords.len() as u32;
                coords.push([(b % bd[0]) as u32, ((b / bd[0]) % bd[1]) as u32, (b / (bd[0] * bd[1])) as u32]);
            }
        }
        SparseMap { bd, table, coords }
    }

    pub(crate) fn slots(&self) -> usize {
        self.coords.len()
    }

    pub(crate) fn bd(&self) -> [usize; 3] {
        self.bd
    }

    pub(crate) fn table_len(&self) -> usize {
        self.table.len()
    }

    /// Rang de la maille (x, y, z) (dans le bloc vide si son bloc est absent).
    pub(crate) fn id(&self, x: usize, y: usize, z: usize) -> usize {
        let b = x / SB + self.bd[0] * (y / SB + self.bd[1] * (z / SB));
        let slot = self.table.get(b).copied().unwrap_or(0) as usize;
        slot * SB3 + (x % SB) + SB * ((y % SB) + SB * (z % SB))
    }

    /// Position de la maille de rang `i`, None dans le bloc vide ou hors de
    /// la boîte `dims`.
    pub(crate) fn pos(&self, i: usize, dims: [usize; 3]) -> Option<[usize; 3]> {
        let c = self.coords[i / SB3];
        if c[0] == u32::MAX {
            return None;
        }
        let l = i % SB3;
        let p = [c[0] as usize * SB + l % SB, c[1] as usize * SB + (l / SB) % SB, c[2] as usize * SB + l / (SB * SB)];
        (p[0] < dims[0] && p[1] < dims[1] && p[2] < dims[2]).then_some(p)
    }

    /// Données des shaders : table des blocs puis coordonnées tassées
    /// (x | y << 10 | z << 20) de chaque rang.
    fn gpu_table(&self) -> Vec<u32> {
        let mut v = self.table.clone();
        v.extend(self.coords.iter().map(|c| if c[0] == u32::MAX { 0 } else { c[0] | (c[1] << 10) | (c[2] << 20) }));
        v
    }
}

/// Tableaux par maille d'une grille, rangés comme sur le GPU.
struct CellSetup {
    /// Stockage creux par blocs (couche pariétale) ; None : boîte entière.
    sparse: Option<SparseMap>,
    /// Mailles stockées.
    ns: usize,
    /// Type de chaque maille : 0 fluide, k + 1 maille pariétale k, SOLID, GHOST.
    wall: Vec<u32>,
    /// Fractions q des liens de chaque maille pariétale (18 par maille).
    qs: Vec<f32>,
    wall_count: usize,
    /// Rang et position de chaque maille pariétale.
    wall_cells: Vec<u32>,
    wall_pos: Vec<[f64; 3]>,
    edge_cos: Vec<f32>,
    /// Mailles qui ne portent que des liens de paroi mince (affichage).
    thin: Vec<u32>,
    solid_count: usize,
    band_fill: Option<Vec<u32>>,
    band_fill_pos: Option<Vec<u32>>,
    band_deep: Option<Vec<bool>>,
    /// Distance à la paroi de chaque maille.
    dist: Vec<f32>,
    /// Étendue en x du corps par tranche y.
    slab: Vec<[f32; 2]>,
    /// Paroi volumique : tables des surfels, et mailles réceptrices de la
    /// loi de paroi (un lien vers la paroi).
    surf: Option<(crate::surfel3d::SurfelSetup, Vec<bool>)>,
}

/// Liens de paroi d'une maille fluide en `p` : fractions q, normales des
/// faces (codées), plus petit cosinus entre elles, et si un voisin est
/// solide. `solid_at` dit si le nœud voisin est solide.
fn wall_links(body: &Body3D, p: [usize; 3], solid_at: &dyn Fn(i64, i64, i64) -> bool) -> Option<([f32; 18], [f32; 18], f32, bool)> {
    let mut row = [0f32; 18];
    let mut nrow = [0f32; 18];
    let mut nraw: Vec<[f64; 3]> = Vec::new();
    let mut beside_solid = false;
    for i in 1..19 {
        let [cx, cy, cz] = C3[i];
        let pf = [p[0] as f64, p[1] as f64, p[2] as f64];
        let d = [-cx as f64, -cy as f64, -cz as f64];
        // Maillage : intersection exacte du lien avec les triangles, y
        // compris entre deux mailles fluides (paroi plus mince qu'une maille).
        let hit = body.mesh.as_ref().and_then(|m| m.segment_hit_normal(pf, d));
        let solid = solid_at(p[0] as i64 - cx as i64, p[1] as i64 - cy as i64, p[2] as i64 - cz as i64);
        beside_solid |= solid;
        let q = match (hit, solid) {
            (Some((t, n)), _) => Some(((t as f32).clamp(0.02, 1.0), Some(n))),
            (None, true) => Some((link_q(body, pf, d), None)),
            (None, false) => None,
        };
        if let Some((q, n)) = q {
            row[i - 1] = q;
            // Sans triangle : gradient de la distance juste devant la paroi.
            let n = n.unwrap_or_else(|| {
                let t = (q as f64 - 0.25).max(0.0);
                body.wall_normal([pf[0] + t * d[0], pf[1] + t * d[1], pf[2] + t * d[2]]).0
            });
            nrow[i - 1] = oct_encode(n);
            nraw.push(n);
        }
    }
    if nraw.is_empty() {
        return None;
    }
    let mut cmin = 1.0f64;
    for a in 0..nraw.len() {
        for b in a + 1..nraw.len() {
            cmin = cmin.min((0..3).map(|t| nraw[a][t] * nraw[b][t]).sum::<f64>());
        }
    }
    Some((row, nrow, cmin as f32, beside_solid))
}

/// Grille boîte (veine ou grille fine qui entoure la maquette).
fn dense_setup(body: &Body3D, nested: bool) -> CellSetup {
    let (nx, ny, nz) = (body.nx, body.ny, body.nz);
    let n = nx * ny * nz;
    let b = body.bounds;
    let idx = |x: usize, y: usize, z: usize| x + nx * (y + ny * z);
    let x0 = ((b.x0.floor() as i64) - 1).max(1) as usize;
    let x1 = ((b.x1.ceil() as i64) + 1).min(nx as i64 - 2) as usize;
    let y0 = ((b.y0.floor() as i64) - 1).max(0) as usize;
    let y1 = ((b.y1.ceil() as i64) + 1).min(ny as i64 - 1) as usize;
    let z0 = ((b.z0.floor() as i64) - 1).max(0) as usize;
    let z1 = ((b.z1.ceil() as i64) + 1).min(nz as i64 - 1) as usize;
    let mut wall = vec![0u32; n];
    let mut solid_count = 0;
    for z in z0..=z1 {
        for y in y0..=y1 {
            for x in x0..=x1 {
                if body.inside(x as f64, y as f64, z as f64) {
                    wall[idx(x, y, z)] = SOLID;
                    solid_count += 1;
                }
            }
        }
    }
    let mut qs: Vec<f32> = Vec::new();
    let mut qn: Vec<f32> = Vec::new();
    let mut edge_cos: Vec<f32> = Vec::new();
    let mut thin: Vec<u32> = Vec::new();
    let mut wall_cells: Vec<u32> = Vec::new();
    let mut wall_pos: Vec<[f64; 3]> = Vec::new();
    {
        // Voisin solide, bords latéraux périodiques (comme les tirages).
        let wall_ro = wall.clone();
        let solid_at = |x: i64, y: i64, z: i64| {
            let x = x.clamp(0, nx as i64 - 1) as usize;
            let y = y.rem_euclid(ny as i64) as usize;
            let z = z.rem_euclid(nz as i64) as usize;
            wall_ro[idx(x, y, z)] == SOLID
        };
        for z in z0.saturating_sub(1)..=(z1 + 1).min(nz - 1) {
            for y in y0.saturating_sub(1)..=(y1 + 1).min(ny - 1) {
                for x in (x0 - 1).max(1)..=(x1 + 1).min(nx - 2) {
                    let c = idx(x, y, z);
                    if wall[c] == SOLID {
                        continue;
                    }
                    if let Some((row, nrow, cmin, beside)) = wall_links(body, [x, y, z], &solid_at) {
                        wall_cells.push(c as u32);
                        wall_pos.push([x as f64, y as f64, z as f64]);
                        wall[c] = wall_cells.len() as u32;
                        qs.extend_from_slice(&row);
                        qn.extend_from_slice(&nrow);
                        edge_cos.push(cmin);
                        if !beside {
                            thin.push(c as u32);
                        }
                    }
                }
            }
        }
    }
    let wall_count = wall_cells.len();
    // Normales des liens rangées après les fractions q (même tampon).
    qs.extend_from_slice(&qn);
    let solid_mask: Vec<bool> = wall.iter().map(|&w| w == SOLID).collect();
    let dist = body.wall_distance(&solid_mask);
    if nested {
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    if x == 0 || y == 0 || z == 0 || x == nx - 1 || y == ny - 1 || z == nz - 1 {
                        wall[idx(x, y, z)] = crate::refine3d::GHOST;
                    }
                }
            }
        }
    }
    let slab = slab_of(ny, wall.iter().enumerate().filter(|(_, &w)| w != 0 && w < crate::refine3d::GHOST || w == SOLID).map(|(c, _)| [c % nx, (c / nx) % ny]));
    CellSetup {
        sparse: None,
        ns: n,
        wall,
        qs,
        wall_count,
        wall_cells,
        wall_pos,
        edge_cos,
        thin,
        solid_count,
        band_fill: None,
        band_fill_pos: None,
        band_deep: None,
        dist,
        slab,
        surf: None,
    }
}

/// Grille boîte à paroi volumique (surfels) : mailles solides = voxels sans
/// volume fluide (les voxels coupés dont le nœud est dans le corps sont
/// calculés), mailles « pariétales » = réceptrices des surfels.
fn surfel_setup(body: &Body3D, nested: bool) -> Result<CellSetup, String> {
    let mesh = body.mesh.as_ref().ok_or("paroi volumique (surfels) : maquette maillée seulement")?;
    if body.ground {
        return Err("paroi volumique (surfels) : sol mobile pas encore pris en charge".into());
    }
    let (nx, ny, nz) = (body.nx, body.ny, body.nz);
    let n = nx * ny * nz;
    let idx = |x: usize, y: usize, z: usize| x + nx * (y + ny * z);
    let tt = std::time::Instant::now();
    // Veine : bords y, z périodiques (x : la maquette n'en approche pas, ou
    // veine périodique de validation) ; grille fine : couronne fantôme.
    let st = crate::surfel3d::prepare(&mesh.tris, [nx, ny, nz], [!nested; 3], nested);
    temps("surfels", tt);
    if std::env::var("SOUFFLERIE_TEMPS").is_ok() {
        eprintln!("  [mise en place] surfels : {:?}", st.stats);
    }
    let mut wall = vec![0u32; n];
    let mut solid_count = 0;
    for c in 0..n {
        if st.vol[c] == 0.0 {
            wall[c] = SOLID;
            solid_count += 1;
        }
    }
    for (r, &c) in st.wall_cells.iter().enumerate() {
        wall[c as usize] = r as u32 + 1;
    }
    let wall_ro = wall.clone();
    let solid_at = |x: i64, y: i64, z: i64| {
        let x = x.clamp(0, nx as i64 - 1) as usize;
        let y = y.rem_euclid(ny as i64) as usize;
        let z = z.rem_euclid(nz as i64) as usize;
        wall_ro[idx(x, y, z)] == SOLID
    };
    let pos: Vec<[usize; 3]> = st.wall_cells.iter().map(|&c| { let c = c as usize; [c % nx, (c / nx) % ny, c / (nx * ny)] }).collect();
    // Liens vers la paroi (loi de paroi, arêtes vives) : mêmes règles que
    // pour le rebond de Bouzidi.
    let links: Vec<Option<f32>> = pos.par_iter().map(|&p| wall_links(body, p, &solid_at).map(|l| l.2)).collect();
    let law: Vec<bool> = links.iter().map(|l| l.is_some()).collect();
    let edge_cos: Vec<f32> = links.iter().map(|l| l.unwrap_or(1.0)).collect();
    let wall_count = st.wall_cells.len();
    let wall_cells = st.wall_cells.clone();
    let wall_pos: Vec<[f64; 3]> = pos.iter().map(|p| p.map(|v| v as f64)).collect();
    // Tampon des liens : poids T/D, débuts des entrées, entrées (case de
    // Γᵒᵘᵗ/P dans `force` : après les efforts des surfels, 12 cases par surfel).
    let ns_ = st.area.len() as u32;
    let mut qs: Vec<f32> = st.recv_a.clone();
    qs.extend(st.recv_start.iter().map(|&v| f32::from_bits(v)));
    for &(k, w) in &st.recv {
        let (sf, p) = (k / 9, k % 9);
        qs.push(f32::from_bits(4 * ns_ + 12 * sf + p));
        qs.push(w);
    }
    let solid_mask: Vec<bool> = wall.iter().map(|&w| w == SOLID).collect();
    let dist = body.wall_distance(&solid_mask);
    if nested {
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    if x == 0 || y == 0 || z == 0 || x == nx - 1 || y == ny - 1 || z == nz - 1 {
                        wall[idx(x, y, z)] = crate::refine3d::GHOST;
                    }
                }
            }
        }
    }
    let slab = slab_of(ny, wall.iter().enumerate().filter(|(_, &w)| w != 0 && w < crate::refine3d::GHOST || w == SOLID).map(|(c, _)| [c % nx, (c / nx) % ny]));
    Ok(CellSetup {
        sparse: None,
        ns: n,
        wall,
        qs,
        wall_count,
        wall_cells,
        wall_pos,
        edge_cos,
        thin: Vec::new(),
        solid_count,
        band_fill: None,
        band_fill_pos: None,
        band_deep: None,
        dist,
        slab,
        surf: Some((st, law)),
    })
}

/// Étendue en x du corps (mailles solides et pariétales) dans chaque tranche y.
fn slab_of(ny: usize, cells: impl Iterator<Item = [usize; 2]>) -> Vec<[f32; 2]> {
    let mut slab = vec![[1.0f32, 0.0f32]; ny];
    for [x, y] in cells {
        let x = x as f32;
        if slab[y][0] > slab[y][1] {
            slab[y] = [x, x];
        } else {
            slab[y] = [slab[y][0].min(x), slab[y][1].max(x)];
        }
    }
    slab
}

/// Couche pariétale d'épaisseur `band` (mailles) en stockage creux : seuls
/// les blocs de 8³ mailles près de la surface sont examinés (distance au
/// maillage du centre du bloc), puis, maille par maille : solide, fluide à
/// moins de `band` de la paroi (active), ou fantôme. Mémoire et calcul de la
/// mise en place proportionnels à la surface, pas au volume de la boîte.
fn band_setup(body: &Body3D, spec: &LayerSpec) -> Result<CellSetup, String> {
    let band = spec.band;
    let mesh = body.mesh.clone().ok_or("couche pariétale : maquette maillée seulement")?;
    // Région étendue (sillage, bords) : chemins ajoutés ; sans eux, mêmes
    // opérations qu'avant, au bit près.
    let extended = spec.wake > 0.0 || spec.edges.is_some() || spec.zone.is_some();
    // Zone : maillage des seuls triangles retenus (repère de cette couche).
    let zmesh = match spec.zone.as_ref() {
        Some(z) => {
            if z.tris.len() != mesh.tris.len() {
                return Err("couche par zones : masque de triangles d'une autre maquette".into());
            }
            let t: Vec<crate::stl::Tri> = mesh.tris.iter().zip(&z.tris).filter(|(_, &k)| k).map(|(t, _)| *t).collect();
            if t.is_empty() {
                return Err("couche par zones : aucun triangle retenu".into());
            }
            Some(crate::stl::MeshBody::new(t, false))
        }
        None => None,
    };
    let near_edge = |p: [f64; 3], slack: f64| -> bool {
        spec.edges.as_ref().is_none_or(|e| e.set.distance(p, e.radius + slack + 1.0) <= e.radius + slack)
            && zmesh.as_ref().is_none_or(|z| z.distance(p, band + slack + 1.0) <= band + slack)
    };
    let (nx, ny, nz) = (body.nx, body.ny, body.nz);
    let dims = [nx, ny, nz];
    let bd = [nx.div_ceil(SB), ny.div_ceil(SB), nz.div_ceil(SB)];
    // Blocs : centre à moins de band + demi-diagonale (√3 · 4 ≈ 6,9) + 2 de
    // la surface (de quoi loger les fantômes voisins des mailles actives).
    let tb = std::time::Instant::now();
    let reach = band + 0.5 * (3.0f64).sqrt() * SB as f64 + 2.0;
    let used: Vec<bool> = (0..bd[0] * bd[1] * bd[2])
        .into_par_iter()
        .map(|b| {
            let c = [b % bd[0], (b / bd[0]) % bd[1], b / (bd[0] * bd[1])].map(|v| (v * SB) as f64 + 0.5 * (SB - 1) as f64);
            let bb = body.bounds;
            let lo = [bb.x0, bb.y0, bb.z0];
            let hi = [bb.x1, bb.y1, bb.z1];
            if (0..3).any(|a| c[a] < lo[a] - reach || c[a] > hi[a] + reach) {
                return false;
            }
            mesh.distance(c, reach + 1.0) <= reach && near_edge(c, reach - band)
        })
        .collect();
    let map = SparseMap::new(bd, &used);
    temps("couche : blocs", tb);
    // Type et distance de chaque maille gardée.
    let ghost = crate::refine3d::GHOST;
    // Distance de chaque maille bornée par celle du centre de son bloc
    // (inégalité triangulaire) : au-delà du plafond sans recherche, sinon
    // recherche arrêtée plus tôt ; même valeur au bit
    // (`MeshBody::distance_bounded`, `distance_beyond`).
    let cap = band + 3.0;
    let hd = 0.5 * (3.0f64).sqrt() * SB as f64;
    let cells: Vec<(u32, f32)> = (0..map.slots())
        .into_par_iter()
        .flat_map_iter(|slot| {
            let c = map.coords[slot].map(|v| v as f64 * SB as f64 + 0.5 * (SB - 1) as f64);
            let mut dc: Option<f64> = None;
            (0..SB3)
                .map(|l| {
                    let i = slot * SB3 + l;
                    let Some(p) = map.pos(i, dims) else { return (ghost, 1e3) };
                    if p[0] == 0 || p[1] == 0 || p[2] == 0 || p[0] == nx - 1 || p[1] == ny - 1 || p[2] == nz - 1 {
                        return (ghost, 1e3);
                    }
                    let pf = p.map(|v| v as f64);
                    if body.inside(pf[0], pf[1], pf[2]) {
                        return (SOLID, 0.0);
                    }
                    let dc = *dc.get_or_insert_with(|| mesh.distance(c, cap + hd + 1.0));
                    let r = ((pf[0] - c[0]).powi(2) + (pf[1] - c[1]).powi(2) + (pf[2] - c[2]).powi(2)).sqrt();
                    let d = if dc - r > cap + 1e-6 { mesh.distance_beyond(pf, cap) } else { mesh.distance_bounded(pf, cap, dc + r) };
                    (if d > band || !near_edge(pf, 0.0) { ghost } else { 0 }, d as f32)
                })
                .collect::<Vec<_>>()
        })
        .collect();
    temps("couche : types et distances", tb);
    // Sillage : dans chaque rangée (y, z), les mailles à moins de `wake`
    // en aval de la dernière maille active deviennent actives (blocs
    // ajoutés au besoin), avec leur distance à la paroi (Spalart–Allmaras).
    let (map, cells) = if spec.wake > 0.0 {
        let w = spec.wake;
        let mut last = vec![-1i64; ny * nz];
        for (i, c) in cells.iter().enumerate() {
            if c.0 == 0 {
                if let Some(p) = map.pos(i, dims) {
                    let r = &mut last[p[1] + ny * p[2]];
                    *r = (*r).max(p[0] as i64);
                }
            }
        }
        let in_wake = |p: [usize; 3]| {
            let x0 = last[p[1] + ny * p[2]];
            x0 >= 0 && (p[0] as i64) > x0 && (p[0] as f64) <= x0 as f64 + w
        };
        let mut used_w = vec![false; bd[0] * bd[1] * bd[2]];
        for c in map.coords.iter().skip(1) {
            used_w[c[0] as usize + bd[0] * (c[1] as usize + bd[1] * c[2] as usize)] = true;
        }
        for z in 0..nz {
            for y in 0..ny {
                let x0 = last[y + ny * z];
                if x0 < 0 {
                    continue;
                }
                let x1 = ((x0 as f64 + w).floor() as usize).min(nx - 1);
                for bx in (x0 as usize + 1) / SB..=x1 / SB {
                    used_w[bx + bd[0] * (y / SB + bd[1] * (z / SB))] = true;
                }
            }
        }
        let map_w = SparseMap::new(bd, &used_w);
        let capw = band + w + 3.0;
        let cells_w: Vec<(u32, f32)> = (0..map_w.slots())
            .into_par_iter()
            .flat_map_iter(|slot| {
                let c = map_w.coords[slot];
                let old = if c[0] == u32::MAX { 0 } else { map.table[c[0] as usize + bd[0] * (c[1] as usize + bd[1] * c[2] as usize)] as usize };
                let (map_w, cells, mesh, in_wake) = (&map_w, &cells, &mesh, &in_wake);
                (0..SB3).map(move |l| {
                    let i = slot * SB3 + l;
                    let Some(p) = map_w.pos(i, dims) else { return (ghost, 1e3) };
                    let boundary = p[0] == 0 || p[1] == 0 || p[2] == 0 || p[0] == nx - 1 || p[1] == ny - 1 || p[2] == nz - 1;
                    let pf = p.map(|v| v as f64);
                    let cur = if old > 0 {
                        cells[old * SB3 + l]
                    } else if boundary {
                        (ghost, 1e3)
                    } else if body.inside(pf[0], pf[1], pf[2]) {
                        (SOLID, 0.0)
                    } else {
                        (ghost, 1e3)
                    };
                    if cur.0 == ghost && !boundary && in_wake(p) {
                        (0, mesh.distance(pf, capw) as f32)
                    } else {
                        cur
                    }
                })
            })
            .collect();
        temps("couche : sillage", tb);
        (map_w, cells_w)
    } else {
        (map, cells)
    };
    let ns = map.slots() * SB3;
    let mut wall: Vec<u32> = cells.iter().map(|c| c.0).collect();
    let dist: Vec<f32> = cells.iter().map(|c| c.1).collect();
    let solid_count = wall.iter().filter(|&&w| w == SOLID).count();
    let at = |x: i64, y: i64, z: i64| -> Option<usize> {
        (x >= 0 && y >= 0 && z >= 0 && (x as usize) < nx && (y as usize) < ny && (z as usize) < nz).then(|| map.id(x as usize, y as usize, z as usize))
    };
    // Liens de paroi des mailles fluides actives (en parallèle, puis numérotées dans l'ordre).
    let links: Vec<Option<([f32; 18], [f32; 18], f32, bool)>> = {
        let wall_ro = &wall;
        let solid_at = |x: i64, y: i64, z: i64| at(x, y, z).is_some_and(|c| wall_ro[c] == SOLID);
        (0..ns)
            .into_par_iter()
            .map(|i| {
                if wall_ro[i] != 0 {
                    return None;
                }
                let p = map.pos(i, dims)?;
                // Plus loin que √2 de la surface, aucun lien (longueur ≤ √2)
                // ne la coupe : seuls des voisins solides en feraient une
                // maille pariétale (même résultat, sans les 18 recherches
                // de triangles).
                if dist[i] > 1.5
                    && !(1..19).any(|k| {
                        let c = C3[k];
                        solid_at(p[0] as i64 - c[0] as i64, p[1] as i64 - c[1] as i64, p[2] as i64 - c[2] as i64)
                    })
                {
                    return None;
                }
                wall_links(body, p, &solid_at)
            })
            .collect()
    };
    temps("couche : liens", tb);
    let mut qs: Vec<f32> = Vec::new();
    let mut qn: Vec<f32> = Vec::new();
    let mut edge_cos: Vec<f32> = Vec::new();
    let mut thin: Vec<u32> = Vec::new();
    let mut wall_cells: Vec<u32> = Vec::new();
    let mut wall_pos: Vec<[f64; 3]> = Vec::new();
    for (i, l) in links.into_iter().enumerate() {
        if let Some((row, nrow, cmin, beside)) = l {
            let p = map.pos(i, dims).unwrap();
            wall_cells.push(i as u32);
            wall_pos.push(p.map(|v| v as f64));
            wall[i] = wall_cells.len() as u32;
            qs.extend_from_slice(&row);
            qn.extend_from_slice(&nrow);
            edge_cos.push(cmin);
            if !beside {
                thin.push(i as u32);
            }
        }
    }
    let wall_count = wall_cells.len();
    qs.extend_from_slice(&qn);
    // Fantômes à remplir : voisins (26) d'une maille active.
    let active = |c: usize| wall[c] < ghost;
    // En parallèle, dans l'ordre des rangs (même liste qu'en séquentiel).
    let found: Vec<(u32, u32)> = (0..ns)
        .into_par_iter()
        .filter_map(|i| {
            if wall[i] != ghost {
                return None;
            }
            let p = map.pos(i, dims)?;
            for dz in -1i64..=1 {
                for dy in -1i64..=1 {
                    for dx in -1i64..=1 {
                        if let Some(c) = at(p[0] as i64 + dx, p[1] as i64 + dy, p[2] as i64 + dz) {
                            if active(c) {
                                return Some((i as u32, (p[0] + nx * (p[1] + ny * p[2])) as u32));
                            }
                        }
                    }
                }
            }
            None
        })
        .collect();
    let (fill, fill_pos): (Vec<u32>, Vec<u32>) = found.into_iter().unzip();
    // Mailles « profondes » (recopiées vers la grille parente) : actives,
    // et à plus de (2 MARGIN + 1) √3 mailles de toute fantôme, la distance à
    // la paroi variant d'au plus une maille par maille parcourue.
    temps("couche : fantômes", tb);
    let r = (2 * crate::refine3d::MARGIN3 + 1) as f64 * (3.0f64).sqrt();
    let deep: Vec<bool> = if extended {
        // Région quelconque : la distance à la paroi ne dit plus où sont les
        // fantômes. Érosion du masque « non fantôme » par un cube de
        // demi-côté ⌈r⌉ (trois passes d'une dimension) : plus strict que la
        // boule de rayon r, donc sûr.
        let rr = r.ceil() as i64;
        let mut m: Vec<bool> = wall.iter().map(|&w| w != ghost).collect();
        for ax in 0..3 {
            let prev = m;
            m = (0..ns)
                .into_par_iter()
                .map(|i| {
                    if !prev[i] {
                        return false;
                    }
                    let Some(p) = map.pos(i, dims) else { return false };
                    (1..=rr).all(|t| {
                        [-t, t].iter().all(|&dt| {
                            let mut q = p.map(|v| v as i64);
                            q[ax] += dt;
                            at(q[0], q[1], q[2]).is_some_and(|c| prev[c])
                        })
                    })
                })
                .collect();
        }
        (0..ns).map(|i| active(i) && m[i]).collect()
    } else {
        (0..ns).map(|i| active(i) && (dist[i] as f64) <= band - r).collect()
    };
    // Blocs gardés pour de bon : ceux qui contiennent une maille active ou
    // une fantôme à remplir (la sélection par le centre du bloc est large).
    let mut keep = vec![false; map.slots()];
    keep[0] = true;
    for i in 0..ns {
        if active(i) {
            keep[i / SB3] = true;
        }
    }
    for &i in &fill {
        keep[i as usize / SB3] = true;
    }
    let mut used2 = vec![false; bd[0] * bd[1] * bd[2]];
    for (slot, c) in map.coords.iter().enumerate().skip(1) {
        if keep[slot] {
            used2[c[0] as usize + bd[0] * (c[1] as usize + bd[1] * c[2] as usize)] = true;
        }
    }
    let map2 = SparseMap::new(bd, &used2);
    let ns2 = map2.slots() * SB3;
    // Nouveau rang de chaque ancien rang (blocs gardés, dans le même ordre).
    let mut new_slot = vec![0usize; map.slots()];
    let mut k = 0;
    for (slot, &kp) in keep.iter().enumerate() {
        if kp {
            new_slot[slot] = k;
            k += 1;
        }
    }
    let remap = |i: usize| new_slot[i / SB3] * SB3 + i % SB3;
    let mut wall2 = vec![ghost; ns2];
    let mut dist2 = vec![1e3f32; ns2];
    let mut deep2 = vec![false; ns2];
    for i in 0..ns {
        if keep[i / SB3] && i >= SB3 {
            let j = remap(i);
            wall2[j] = wall[i];
            dist2[j] = dist[i];
            deep2[j] = deep[i];
        }
    }
    let (wall, dist, deep, ns) = (wall2, dist2, deep2, ns2);
    if std::env::var("SOUFFLERIE_TEMPS").is_ok() {
        let act = wall.iter().filter(|&&w| w < ghost).count();
        eprintln!("  [mise en place] couche : {act} mailles actives, {ns} stockées (rapport {:.2})", ns as f64 / act.max(1) as f64);
    }
    let wall_cells: Vec<u32> = wall_cells.iter().map(|&i| remap(i as usize) as u32).collect();
    let thin: Vec<u32> = thin.iter().map(|&i| remap(i as usize) as u32).collect();
    let fill: Vec<u32> = fill.iter().map(|&i| remap(i as usize) as u32).collect();
    let map = map2;
    let slab = slab_of(
        ny,
        (0..ns).filter(|&i| wall[i] != 0 && wall[i] < ghost || wall[i] == SOLID).filter_map(|i| map.pos(i, dims)).map(|p| [p[0], p[1]]),
    );
    Ok(CellSetup {
        sparse: Some(map),
        ns,
        wall,
        qs,
        wall_count,
        wall_cells,
        wall_pos,
        edge_cos,
        thin,
        solid_count,
        band_fill: Some(fill),
        band_fill_pos: Some(fill_pos),
        band_deep: Some(deep),
        dist,
        slab,
        surf: None,
    })
}

/// Hauteur de raccord de la loi de paroi, en mailles de la grille.
pub const WM_H: f64 = 2.0;

/// u_τ de la loi de Spalding pour la vitesse tangentielle `ut` à la distance `y`.

pub fn spalding_utau(ut: f64, y: f64, nu: f64) -> f64 {
    if !(ut > 1e-12) {
        return 0.0;
    }
    let (k, eb) = (0.41, 0.128_734_90);
    let rey = y * ut / nu;
    let mut up = rey.sqrt();
    if rey >= 30.0 {
        up = rey.ln() / k;
    }
    for _ in 0..50 {
        let ku = k * up;
        let e = ku.exp();
        let yy = up + eb * (e - 1.0 - ku - ku * ku / 2.0 - ku * ku * ku / 6.0);
        let dy = 1.0 + eb * k * (e - 1.0 - ku - ku * ku / 2.0);
        let h = up * yy - rey;
        let next = up - h / (yy + up * dy);
        up = if next > 0.0 { next } else { 0.5 * up };
        if h.abs() < 1e-10 * rey {
            break;
        }
    }
    ut / up
}

/// Écart de la loi de Townsend–Mellor à la loi sans gradient en y (comme
/// `mellor_du` de `WALL3_WGSL`) : (2/κ)[√τ − u_τ + u_τ ln(2u_τ/(√τ + u_τ))],
/// τ = u_τ² + P y borné à 0, P = (1/ρ) dp/ds.
pub fn mellor_du(u: f64, pg: f64, y: f64) -> f64 {
    let k = 0.41;
    if !(u > 1e-12) {
        return 2.0 / k * (pg * y).max(0.0).sqrt();
    }
    let r = (pg * y / (u * u)).max(-1.0);
    let a = r / ((1.0 + r).sqrt() + 1.0);
    2.0 / k * u * (a - (0.5 * a).ln_1p())
}

/// u_τ de la loi de paroi avec gradient de pression : u_τ = u_τ,Spalding(u_t − Δ(u_τ)),
/// résolu par dichotomie dans l'encadrement de `mellor_utau` du shader.
pub fn mellor_utau(ut: f64, y: f64, nu: f64, pg: f64) -> f64 {
    let u0 = spalding_utau(ut, y, nu);
    if pg == 0.0 || !(u0 > 1e-12) {
        return u0;
    }
    let f = |u: f64| u - spalding_utau((ut - mellor_du(u, pg, y)).max(0.0), y, nu);
    let (mut lo, mut hi) = if pg > 0.0 { (0.0, u0) } else { (u0, 2.0 * u0) };
    if f(hi) <= 0.0 {
        return hi;
    }
    if f(lo) >= 0.0 {
        return lo;
    }
    for _ in 0..60 {
        let m = 0.5 * (lo + hi);
        if f(m) < 0.0 {
            lo = m;
        } else {
            hi = m;
        }
    }
    0.5 * (lo + hi)
}

pub struct WallForces {
    pub pos: Vec<[f64; 3]>,
    pub force: Vec<[f64; 3]>,
}

fn storage(binding: u32, read_only: bool) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Storage { read_only },
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    }
}

fn uniform(binding: u32, dynamic: bool, size: u64) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Uniform,
            has_dynamic_offset: dynamic,
            min_binding_size: wgpu::BufferSize::new(size),
        },
        count: None,
    }
}

fn buf_entry(binding: u32, buffer: &wgpu::Buffer, size: Option<u64>) -> wgpu::BindGroupEntry<'_> {
    wgpu::BindGroupEntry {
        binding,
        resource: wgpu::BindingResource::Buffer(wgpu::BufferBinding {
            buffer,
            offset: 0,
            size: size.and_then(wgpu::BufferSize::new),
        }),
    }
}

/// Normale unitaire codée en un flottant (projection octaédrique, 2 × 11 bits
/// en un entier exact) ; `oct_decode` dans les shaders.
fn oct_encode(n: [f64; 3]) -> f32 {
    let s = n[0].abs() + n[1].abs() + n[2].abs();
    let (mut u, mut v) = if s > 0.0 { (n[0] / s, n[1] / s) } else { (0.0, 0.0) };
    if n[2] < 0.0 {
        let (a, b) = (u, v);
        u = (1.0 - b.abs()) * a.signum();
        v = (1.0 - a.abs()) * b.signum();
    }
    let q = |x: f64| ((x * 0.5 + 0.5) * 2047.0).round().clamp(0.0, 2047.0);
    (q(u) * 2048.0 + q(v)) as f32
}

/// Fraction q ∈ ]0, 1] du lien p → p + d côté fluide (dichotomie).
fn link_q(body: &Body3D, p: [f64; 3], d: [f64; 3]) -> f32 {
    let (mut lo, mut hi) = (0.0, 1.0);
    for _ in 0..12 {
        let m = 0.5 * (lo + hi);
        if body.inside(p[0] + m * d[0], p[1] + m * d[1], p[2] + m * d[2]) {
            hi = m;
        } else {
            lo = m;
        }
    }
    (0.5 * (lo + hi)).clamp(0.02, 1.0) as f32
}

impl Gpu3D {
    /// Veine de la maquette `body`, viscosité `nu` (réseau), vitesse d'entrée `u_in`.
    pub fn new(body: &Body3D, nu: f64, u_in: f32, fp16: bool) -> Result<Gpu3D, String> {
        Self::with_turbulence(body, nu, u_in, fp16, false)
    }

    /// Idem, avec Spalart–Allmaras et loi de paroi si `sa`.
    pub fn with_turbulence(body: &Body3D, nu: f64, u_in: f32, fp16: bool, sa: impl Into<Turbulence>) -> Result<Gpu3D, String> {
        Self::with_refinement(body, nu, u_in, fp16, sa, 0)
    }

    /// Idem, avec `depth` grilles imbriquées autour de la maquette (0 à 2),
    /// chacune deux fois plus fine que la précédente.
    pub fn with_refinement(
        body: &Body3D,
        nu: f64,
        u_in: f32,
        fp16: bool,
        sa: impl Into<Turbulence>,
        depth: usize,
    ) -> Result<Gpu3D, String> {
        Self::with_refinement_le(body, nu, u_in, fp16, sa, depth, false)
    }

    /// Idem ; `le` ajoute une grille deux fois plus fine encore, limitée au
    /// bord d'attaque d'une maquette maillée (faces tournées vers l'amont,
    /// n_x < −0,5) : le rayon du bord d'attaque d'une aile (≈ 1 % de la
    /// corde) tient sinon dans moins d'une maille, et le pic de dépression,
    /// qui équilibre la traînée de pression de l'arrière de l'extrados, est
    /// écrêté. Les efforts de la surface viennent, triangle par triangle, de
    /// la grille la plus fine qui l'entoure.
    pub fn with_refinement_le(
        body: &Body3D,
        nu: f64,
        u_in: f32,
        fp16: bool,
        sa: impl Into<Turbulence>,
        depth: usize,
        le: bool,
    ) -> Result<Gpu3D, String> {
        Self::with_refinement_edges(body, nu, u_in, fp16, sa, depth, le, false, &[])
    }

    /// Idem, avec des grilles deux fois plus fines encore, limitées aux bords
    /// d'attaque (`le`) et de fuite (`te`) d'une maquette maillée : plusieurs
    /// boîtes qui suivent chaque bord en envergure (voir `edge_boxes`). Au
    /// bord de fuite, la portance dépend de la place du bord dans la maille
    /// (condition de Kutta discrète : ±8 % à 32 mailles par corde, ±3 % à 64,
    /// aile AR 4) ; l'erreur, d'ordre un en maille, se réduit en raffinant là.
    #[allow(clippy::too_many_arguments)]
    pub fn with_refinement_edges(
        body: &Body3D,
        nu: f64,
        u_in: f32,
        fp16: bool,
        sa: impl Into<Turbulence>,
        depth: usize,
        le: bool,
        te: bool,
        bands: &[f64],
    ) -> Result<Gpu3D, String> {
        let depth = depth.min(crate::refine3d::EXTENT3.len());
        let plan = RefinePlan {
            boxes: crate::refine3d::EXTENT3[..depth].to_vec(),
            layers: nested_bands(bands).into_iter().map(|band| LayerSpec { band, ..Default::default() }).collect(),
            le,
            te,
        };
        Self::with_plan(body, nu, u_in, fp16, sa, &plan)
    }

    /// Idem, d'après un plan de raffinement : grilles boîtes d'emprise
    /// donnée, puis couches pariétales (épaisseur, sillage, bords), ou
    /// grilles des bords. `with_refinement_edges` et le raffinement
    /// automatique (`raff_auto`) en sont deux cas.
    pub fn with_plan(body: &Body3D, nu: f64, u_in: f32, fp16: bool, sa: impl Into<Turbulence>, plan: &RefinePlan) -> Result<Gpu3D, String> {
        let sa: Turbulence = sa.into();
        let (le, te) = (plan.le, plan.te);
        let (device, queue, name) = Self::open_device()?;
        let mut top = Self::build(device.clone(), queue.clone(), name.clone(), body, nu, u_in, fp16, sa, false, &LayerSpec::default())?;
        let depth = plan.boxes.len();
        if depth == 0 {
            return Ok(top);
        }
        if body.ground {
            return Err("le raffinement local ne gère pas encore le sol mobile".into());
        }
        let mut levels: Vec<Gpu3D> = Vec::new();
        let mut links = Vec::new();
        let mut parent_body = body.clone();
        let l0 = body.l as f64;
        let mut origin: Vec<[f64; 3]> = Vec::new();
        let mut porg = [0.0f64; 3];
        for k in 1..=depth {
            let (up, down, side) = plan.boxes[k - 1];
            // Emprise en mailles de la grille parente (échelle 2^(k−1)).
            let sp = (1u64 << (k - 1)) as f64;
            let pb = parent_body.bounds;
            let (pnx, pny, pnz) = (parent_body.nx, parent_body.ny, parent_body.nz);
            let lo = [
                (pb.x0 - up * l0 * sp).floor().max(3.0) as usize,
                (pb.y0 - side * l0 * sp).floor().max(3.0) as usize,
                (pb.z0 - side * l0 * sp).floor().max(3.0) as usize,
            ];
            let hi = [
                ((pb.x1 + down * l0 * sp).ceil() as usize).min(pnx - 4),
                ((pb.y1 + side * l0 * sp).ceil() as usize).min(pny - 4),
                ((pb.z1 + side * l0 * sp).ceil() as usize).min(pnz - 4),
            ];
            let m = 2 * crate::refine3d::MARGIN3 + 2;
            if (0..3).any(|a| hi[a] <= lo[a] + m) {
                return Err("grille fine trop petite".into());
            }
            let dims = [0, 1, 2].map(|a| 2 * (hi[a] - lo[a]) + 1);
            let cb = parent_body.refined([lo[0] as f64, lo[1] as f64, lo[2] as f64], dims);
            let nu_k = nu * (1u64 << k) as f64;
            let mut child = Self::build(device.clone(), queue.clone(), name.clone(), &cb, nu_k, u_in, fp16, sa, true, &LayerSpec::default())?;
            child.les = top.les;
            {
                let parent = if k == 1 { &top } else { &levels[k - 2] };
                links.push({ let tl = std::time::Instant::now(); let l = crate::refine3d::Link3::new(parent, &child, lo); temps("raccord", tl); l });
            }
            levels.push(child);
            parent_body = cb;
            porg = [0, 1, 2].map(|a| porg[a] + lo[a] as f64 / sp);
            origin.push(porg);
        }
        let full = levels.len();
        let mut parent = (0..full).collect::<Vec<usize>>();
        let mut lev_depth: Vec<u32> = (1..=full as u32).collect();
        let mut core: Vec<Option<[[f64; 3]; 2]>> = vec![None; full];
        if !plan.layers.is_empty() {
            // Couches pariétales imbriquées : chacune deux fois plus fine que
            // la précédente, sur une épaisseur `bands[i]` (mailles de cette
            // couche) autour de toute la surface ; boîte qui entoure la
            // maquette, calcul et stockage limités à la couche.
            for (bi, spec) in plan.layers.iter().enumerate() {
                let band = spec.band;
                let k = depth + 1 + bi;
                let sp = (1u64 << (k - 1)) as f64;
                let pb = parent_body.bounds;
                let (pnx, pny, pnz) = (parent_body.nx, parent_body.ny, parent_body.nz);
                // Marge (mailles parentes) : demi-épaisseur de la couche, plus
                // la bande de raccord.
                let mg = 0.5 * band + (2 * crate::refine3d::MARGIN3 + 3) as f64;
                let lo = [
                    (pb.x0 - mg).floor().max(3.0) as usize,
                    (pb.y0 - mg).floor().max(3.0) as usize,
                    (pb.z0 - mg).floor().max(3.0) as usize,
                ];
                // Sillage : la couche s'étend vers l'aval de `wake` mailles.
                let hi = [
                    ((pb.x1 + mg + 0.5 * spec.wake).ceil() as usize).min(pnx - 4),
                    ((pb.y1 + mg).ceil() as usize).min(pny - 4),
                    ((pb.z1 + mg).ceil() as usize).min(pnz - 4),
                ];
                let dims = [0, 1, 2].map(|a| 2 * (hi[a] - lo[a]) + 1);
                let cb = parent_body.refined([lo[0] as f64, lo[1] as f64, lo[2] as f64], dims);
                let nu_k = nu * (1u64 << k) as f64;
                // Arêtes des bords : en mailles de cette couche.
                let o = [0, 1, 2].map(|a| porg[a] + lo[a] as f64 / sp);
                let sk = (1u64 << k) as f64;
                let spec_k = LayerSpec {
                    edges: spec.edges.as_ref().map(|e| {
                        std::sync::Arc::new(EdgeRegion {
                            set: e.set.transformed(|q| [0, 1, 2].map(|a| (q[a] - o[a]) * sk)),
                            radius: e.radius,
                        })
                    }),
                    ..spec.clone()
                };
                let mut child = Self::build(device.clone(), queue.clone(), name.clone(), &cb, nu_k, u_in, fp16, sa, true, &spec_k)?;
                child.les = top.les;
                // Mémoire du GPU : environ 100 octets par maille en demi-précision
                // (200 en simple) ; au-delà de 40 M mailles en tout, on s'arrête.
                let total: usize = std::iter::once(top.n).chain(levels.iter().map(|g| g.n)).chain(std::iter::once(child.n)).sum();
                if total > 40_000_000 / if fp16 { 1 } else { 2 } {
                    return Err(format!(
                        "couche pariétale {} (×{}) : {:.0} M mailles en tout, au-delà de la mémoire du GPU ; moins de couches ou plus minces",
                        bi + 1,
                        1u64 << k,
                        total as f64 / 1e6
                    ));
                }
                let pi = levels.len();
                {
                    let par = if pi == 0 { &top } else { &levels[pi - 1] };
                    links.push({ let tl = std::time::Instant::now(); let l = crate::refine3d::Link3::new(par, &child, lo); temps("raccord", tl); l });
                }
                levels.push(child);
                origin.push(o);
                parent.push(pi);
                lev_depth.push(k as u32);
                // Domaine propre : toute la boîte (les efforts viennent des
                // mailles pariétales, toutes dans la couche ; la couche la
                // plus fine l'emporte).
                core.push(Some([o, [0, 1, 2].map(|a| o[a] + (dims[a] - 1) as f64 / (2.0 * sp))]));
                parent_body = cb;
                porg = o;
            }
        } else if le || te {
            let k = depth + 1;
            let sp = (1u64 << (k - 1)) as f64;
            let mesh = parent_body.mesh.clone().ok_or("raffinement des bords : maquette maillée seulement")?;
            let (pnx, pny, pnz) = (parent_body.nx, parent_body.ny, parent_body.nz);
            let boxes = edge_boxes(&mesh.tris, parent_body.l as f64, [pnx, pny, pnz], le, te)?;
            let n_tot: usize = boxes.iter().map(|b| (0..3).map(|a| 2 * (b.1[a] - b.0[a]) + 1).product::<usize>()).sum();
            // Au-delà de quatre fois la grille parente, les boîtes coûtent
            // autant qu'un raffinement complet (environ cinq fois) et risquent
            // d'épuiser la mémoire du GPU.
            if n_tot > 4 * pnx * pny * pnz {
                return Err(format!(
                    "raffinement des bords : {:.0} M mailles, plus de quatre fois la grille parente : raffiner toute la maquette (--raff) plutôt",
                    n_tot as f64 / 1e6
                ));
            }
            for (lo, hi, cr) in boxes {
                let dims = [0, 1, 2].map(|a| 2 * (hi[a] - lo[a]) + 1);
                let cb = parent_body.refined([lo[0] as f64, lo[1] as f64, lo[2] as f64], dims);
                let nu_k = nu * (1u64 << k) as f64;
                let mut child = Self::build(device.clone(), queue.clone(), name.clone(), &cb, nu_k, u_in, fp16, sa, true, &LayerSpec::default())?;
                child.les = top.les;
                {
                    let par = if full == 0 { &top } else { &levels[full - 1] };
                    links.push({ let tl = std::time::Instant::now(); let l = crate::refine3d::Link3::new(par, &child, lo); temps("raccord", tl); l });
                }
                levels.push(child);
                origin.push([0, 1, 2].map(|a| porg[a] + lo[a] as f64 / sp));
                parent.push(full);
                lev_depth.push(k as u32);
                // Domaine propre en mailles de la veine.
                core.push(Some([[0, 1, 2].map(|a| porg[a] + cr[0][a] / sp), [0, 1, 2].map(|a| porg[a] + cr[1][a] / sp)]));
            }
        }
        // Mailles des grilles parentes recopiées des grilles fines et dont
        // tous les voisins le sont aussi : jamais lues avant d'être écrasées.
        if std::env::var("SOUFFLERIE_RECOUVERT").map_or(true, |v| v != "0") {
            for g in 0..=levels.len() {
                let lists: Vec<&[u32]> = links.iter().zip(&parent).filter(|(_, &p)| p == g).map(|(l, _)| l.restrict_list()).collect();
                if lists.is_empty() {
                    continue;
                }
                let grid = if g == 0 { &top } else { &levels[g - 1] };
                grid.mark_covered(&lists);
            }
        }
        top.refine = Some(Box::new(crate::refine3d::Refinement {
            levels,
            links,
            fresh: true,
            origin,
            full,
            parent,
            depth: lev_depth,
            core,
        }));
        Ok(top)
    }

    /// Ouvre le GPU : Vulkan ou Metal d'abord.
    fn open_device() -> Result<(wgpu::Device, wgpu::Queue, String), String> {
        // Vulkan ou Metal d'abord : sous DirectX 12, ce gros shader (tableaux
        // locaux indexés) est environ 100 fois plus lent. La variable
        // WGPU_BACKEND, que le lanceur de l'interface fixe pour la fenêtre, est
        // ignorée ; SOUFFLERIE_GPU3D=dx12|vulkan|gl en impose un.
        let forced = std::env::var("SOUFFLERIE_GPU3D").ok().map(|s| wgpu::Backends::from_comma_list(&s));
        let tries = match forced {
            Some(b) => vec![b],
            None => vec![wgpu::Backends::VULKAN | wgpu::Backends::METAL, wgpu::Backends::all()],
        };
        let mut found = None;
        let mut last_err = String::new();
        for backends in tries {
            let mut d = wgpu::InstanceDescriptor::new_without_display_handle();
            d.backends = backends;
            let instance = wgpu::Instance::new(d);
            match pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                ..Default::default()
            })) {
                Ok(a) => {
                    found = Some(a);
                    break;
                }
                Err(e) => last_err = e.to_string(),
            }
        }
        let adapter = found.ok_or_else(|| format!("aucun adaptateur GPU : {last_err}"))?;
        let info = adapter.get_info();
        let adapter_name = format!("{} ({:?})", info.name, info.backend);
        let al = adapter.limits();
        let mut limits = wgpu::Limits::default();
        limits.max_storage_buffer_binding_size = al.max_storage_buffer_binding_size;
        limits.max_buffer_size = al.max_buffer_size;
        limits.max_storage_buffers_per_shader_stage = limits.max_storage_buffers_per_shader_stage.max(16);
        let i16 = adapter.features().contains(wgpu::Features::SHADER_I16);
        POP16.store(i16 && std::env::var("SOUFFLERIE_U16").map_or(true, |v| v != "0"), std::sync::atomic::Ordering::Relaxed);
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            label: Some("soufflerie 3D"),
            required_limits: limits,
            required_features: if i16 { wgpu::Features::SHADER_I16 } else { wgpu::Features::empty() },
            ..Default::default()
        }))
        .map_err(|e| format!("GPU indisponible : {e}"))?;
        Ok((device, queue, adapter_name))
    }

    /// Une grille (veine ou grille fine imbriquée) sur un périphérique ouvert.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn build(
        device: wgpu::Device,
        queue: wgpu::Queue,
        adapter_name: String,
        body: &Body3D,
        nu: f64,
        u_in: f32,
        fp16: bool,
        turb: Turbulence,
        nested: bool,
        layer: &LayerSpec,
    ) -> Result<Gpu3D, String> {
        let sa = turb != Turbulence::Laminar;
        let vles = turb == Turbulence::Vles;
        let keps = turb == Turbulence::KEpsRng || vles;
        let band = layer.band;
        let (nx, ny, nz) = (body.nx, body.ny, body.nz);
        let al = device.limits();
        let words_a = if fp16 { 5 } else { 10 };
        let words_b = if fp16 { 5 } else { 9 };
        let max_bind = (al.max_storage_buffer_binding_size as u64).min(al.max_buffer_size);
        let st = wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC;
        let mk = |label: &str, size: u64, usage: wgpu::BufferUsages| {
            device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label),
                size: size.max(16),
                usage,
                mapped_at_creation: false,
            })
        };
        let params = mk(
            "params3d",
            SLOT * BATCH as u64,
            wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        );

        // Voxelisation, liens de paroi et couche pariétale : tableaux par
        // maille rangés comme sur le GPU (boîte, ou blocs près de la surface).
        let tt = std::time::Instant::now();
        let cs = if body.surfels {
            if nested && band > 0.0 {
                return Err("paroi volumique (surfels) : pas de couche pariétale creuse (grilles boîtes seulement)".into());
            }
            surfel_setup(body, nested)?
        } else if nested && band > 0.0 {
            band_setup(body, layer)?
        } else {
            dense_setup(body, nested)
        };
        temps("mailles", tt);
        let CellSetup {
            sparse,
            ns,
            wall: wall_gpu,
            mut qs,
            wall_count,
            wall_cells,
            wall_pos,
            edge_cos,
            thin,
            solid_count,
            band_fill,
            band_fill_pos,
            band_deep,
            dist: dist_gpu,
            slab,
            surf,
        } = cs;
        // Populations une par case de 16 bits sur les grilles boîtes ; les
        // couches creuses gardent deux populations par mot (12 % plus rapide
        // pour elles, voir `POP16`).
        let p16 = fp16 && pop16() && sparse.is_none();
        let biggest = (words_a * ns * 4).max(16 * ns) as u64;
        if biggest > max_bind {
            return Err(format!(
                "{:.1} M mailles : tampon de {} Mo, au-delà de la limite de l'adaptateur ({} Mo)",
                ns as f64 / 1e6,
                biggest >> 20,
                max_bind >> 20
            ));
        }
        let pops = [
            [mk("a0", (words_a * ns * 4) as u64, st), mk("b0", (words_b * ns * 4) as u64, st)],
            [mk("a1", (words_a * ns * 4) as u64, st), mk("b1", (words_b * ns * 4) as u64, st)],
        ];
        let macro_buf = mk("macro3d", 16 * ns as u64, st);
        // Table des blocs pour les shaders (vide : grille boîte).
        let blk_data: Vec<u32> = match sparse.as_ref() {
            Some(m) => m.gpu_table(),
            None => vec![0u32; 4],
        };
        let blk_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("blocs3d"),
            contents: bytemuck::cast_slice(&blk_data),
            usage: wgpu::BufferUsages::STORAGE,
        });
        // Maille de chaque lien pariétal (numéro → maille), pour les efforts par maille.
        let rest_force: Vec<[f32; 3]> = if let Some((st, _)) = &surf {
            st.rest_force.iter().map(|r| r.map(|v| v as f32)).collect()
        } else {
            qs
            .chunks_exact(18)
            .take(wall_count)
            .map(|row| {
                let mut r = [0f32; 3];
                for i in 1..19 {
                    if row[i - 1] > 0.0 {
                        for a in 0..3 {
                            r[a] -= (2.0 * w3(i)) as f32 * C3[i][a] as f32;
                        }
                    }
                }
                r
            })
            .collect()
        };
        if qs.is_empty() {
            qs.resize(18, 0.0);
        }
        let solid: Vec<bool> = wall_gpu.iter().map(|&w| w == SOLID).collect();
        let ghost: Vec<bool> = if sparse.is_some() { wall_gpu.iter().map(|&w| w == crate::refine3d::GHOST).collect() } else { Vec::new() };
        // Affichage des mailles du calcul : 2 bits par maille (1 solide, 2 paroi mince).
        let mut vox = vec![0u32; ns.div_ceil(16)];
        for (c, &s) in solid.iter().enumerate() {
            if s {
                vox[c >> 4] |= 1 << ((c & 15) * 2);
            }
        }
        for &c in &thin {
            let c = c as usize;
            vox[c >> 4] |= 2 << ((c & 15) * 2);
        }
        let vox_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("vox3d"),
            contents: bytemuck::cast_slice(&vox),
            usage: wgpu::BufferUsages::STORAGE,
        });
        // Affichage : distance signée tous les L/64 (≥ 0,2 maille, ≤ 0,5).
        let g = if nested {
            crate::bodies3d::SdfGrid {
                origin: [0.0; 3],
                h: 1.0,
                dims: [1, 1, 1],
                d: vec![1.0],
            }
        } else {
            body.sdf_grid((body.l as f64 / 64.0).clamp(0.2, 0.5))
        };
        let sdf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("sdf3d"),
            contents: bytemuck::cast_slice(&g.d),
            usage: wgpu::BufferUsages::STORAGE,
        });
        let sdf_geom = (
            [g.origin[0] as f32, g.origin[1] as f32, g.origin[2] as f32, g.h as f32],
            // Surface affichée : S < w ; les maillages sont épaissis de 0,6 pas de grille
            // pour que leurs parois plus minces que la grille restent visibles.
            [g.dims[0] as f32, g.dims[1] as f32, g.dims[2] as f32, if body.mesh.is_some() { 0.6 * g.h as f32 } else { 0.0 }],
        );
        let plain_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("plain3d"),
            contents: bytemuck::cast_slice(&plain_mask(&wall_gpu)),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        });
        let wall_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("wall3d"),
            contents: bytemuck::cast_slice(&wall_gpu),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST,
        });
        let qbuf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("q3d"),
            contents: bytemuck::cast_slice(&qs),
            usage: wgpu::BufferUsages::STORAGE,
        });
        // Spalart–Allmaras : normale et distance des mailles pariétales (dans
        // l'ordre de leur numéro).
        let geo: Vec<[f32; 4]> = if sa {
            // Distance négative : maille d'arête vive (faces touchées à plus
            // de 60° l'une de l'autre).
            wall_pos
                .par_iter()
                .enumerate()
                .map(|(w, &p)| {
                    let (nrm, d) = body.wall_normal(p);
                    let d = (d as f32).max(0.1);
                    // Avec les arêtes vives du maillage : seulement à portée
                    // de lien (√2) de l'une d'elles (pas autour d'un bord
                    // d'attaque arrondi plus fin qu'une maille).
                    let sharp = edge_cos[w] < 0.5
                        && body.features.as_ref().is_none_or(|f| f.distance(p, 2.0) <= std::f64::consts::SQRT_2);
                    let d = if sharp { -d } else { d };
                    // Surfels : réceptrice sans lien vers la paroi, hors loi de paroi.
                    let d = if surf.as_ref().is_some_and(|(_, law)| !law[w]) { 0.0 } else { d };
                    [nrm[0] as f32, nrm[1] as f32, nrm[2] as f32, d]
                })
                .collect()
        } else {
            Vec::new()
        };
        let slab_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("slab3d"),
            contents: bytemuck::cast_slice(&slab),
            usage: wgpu::BufferUsages::STORAGE,
        });
        let init_or = |label: &str, data: &[u8]| {
            if data.is_empty() {
                mk(label, 16, st)
            } else {
                device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some(label),
                    contents: data,
                    usage: st,
                })
            }
        };
        let wall_geo = init_or("wallgeo3d", bytemuck::cast_slice(&geo));
        let dist_buf = init_or("dist3d", bytemuck::cast_slice(&dist_gpu));
        // Par maille pariétale : u_τ et vitesse de glissement relaxée, puis
        // (à la suite, une valeur par maille) la vitesse de transpiration,
        // puis la viscosité de la couche pariétale (passe `wall3d`).
        let utau_buf = mk("utau3d", (24 * wall_count) as u64, st);
        let nsa = if sa { 4 * ns as u64 * (1 + keps as u64 + vles as u64) } else { 16 };
        let nut = [mk("nut0", nsa, st), mk("nut1", nsa, st)];
        let nut0: Vec<f32> = if keps {
            let (k, e) = crate::sa::ke_inf(u_in as f64, nu);
            let mut v: Vec<f32> = wall_gpu.iter().flat_map(|&w| if w == SOLID { [0.0, 0.0] } else { [k, e] }).collect();
            if vles {
                // ν_t amont : le filtre n'agit pas (Δ ε∞/k∞^{3/2} ≫ 1).
                v.extend(wall_gpu.iter().map(|&w| if w == SOLID { 0.0 } else { crate::sa::nu_ke(k, e) }));
            }
            v
        } else if sa {
            let inf = crate::sa::nu_tilde_inf(nu);
            wall_gpu.iter().map(|&w| if w == SOLID { 0.0 } else { inf }).collect()
        } else {
            Vec::new()
        };
        // Surfels : un effort par surfel, puis Γᵒᵘᵗ/P (12 cases par surfel).
        let n_surf = surf.as_ref().map_or(0, |(st, _)| st.area.len());
        let fsize = if surf.is_some() { (16 * 4 * n_surf).max(16) } else { (16 * wall_count).max(16) } as u64;
        let force = mk("force3d", fsize, st);
        let force_read = mk("force3d_read", if surf.is_some() { (16 * n_surf).max(16) as u64 } else { fsize }, wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST);

        let psize = std::mem::size_of::<Params>() as u64;
        let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("d3q19"),
            entries: &[
                uniform(0, true, psize),
                storage(1, true),
                storage(2, true),
                storage(3, false),
                storage(4, false),
                storage(5, true),
                storage(6, true),
                storage(7, false),
                storage(8, false),
                storage(9, true),
                storage(10, true),
                storage(11, true),
                storage(12, true),
                storage(13, true),
                storage(14, true),
                storage(15, true),
                storage(16, true),
            ],
        });
        let sa_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("sa3d"),
            entries: &[
                uniform(0, true, psize),
                storage(1, true),
                storage(2, true),
                storage(3, true),
                storage(4, true),
                storage(5, true),
                storage(6, true),
                storage(7, false),
                storage(8, true),
                storage(9, true),
                storage(10, true),
            ],
        });
        let wall_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("wall3d"),
            entries: &[
                uniform(0, true, psize),
                storage(1, true),
                storage(2, true),
                storage(3, true),
                storage(4, false),
                storage(5, true),
                storage(6, true),
                storage(7, true),
                storage(8, true),
            ],
        });
        let init_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("d3q19_init"),
            entries: &[uniform(0, true, psize), storage(1, false), storage(2, false)],
        });
        // Passes du calcul : sans contrôles d'exécution (`compute_module`) ;
        // affichage (coupe, rendu) : contrôles gardés.
        let pipe = |label: &str, layout: &wgpu::BindGroupLayout, code: String, trusted: bool| {
            let pl = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some(label),
                bind_group_layouts: &[Some(layout)],
                immediate_size: 0,
            });
            let module = if trusted {
                compute_module(&device, label, code)
            } else {
                device.create_shader_module(wgpu::ShaderModuleDescriptor {
                    label: Some(label),
                    source: wgpu::ShaderSource::Wgsl(code.into()),
                })
            };
            device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some(label),
                layout: Some(&pl),
                module: &module,
                entry_point: Some("main"),
                compilation_options: Default::default(),
                cache: None,
            })
        };
        temps("tampons", tt);
        let step_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("d3q19"),
            bind_group_layouts: &[Some(&layout)],
            immediate_size: 0,
        });
        // Collision spécialisée : paramètres fixes de la grille en constantes
        // (sous-maille selon la viscosité, comme `set_viscosity`).
        let spec = [((0.5 + 3.0 * nu) < 0.51) as u32, nested as u32 | (((body.ground && !nested) as u32) << 1), 0];
        let sparse_bnx = sparse.as_ref().map_or(0, |m| m.bd[0] as u32);
        // Grille creuse : le cœur de la collision ne lance que ses mailles à
        // calculer (fluides sans lien de paroi), rangées après les mailles
        // pariétales ; dans une couche, un tiers des mailles stockées sont
        // fantômes ou solides et occupaient des fils pour rien.
        let active: Vec<u32> = if sparse.is_some() {
            (0..ns as u32).filter(|&i| wall_gpu[i as usize] == 0).collect()
        } else {
            Vec::new()
        };
        let n_active = active.len() as u32;
        let step_list: Vec<u32> = wall_cells.iter().chain(&active).copied().chain(std::iter::once(0)).collect();
        let step_list_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("step_list3d"),
            contents: bytemuck::cast_slice(&step_list),
            usage: wgpu::BufferUsages::STORAGE,
        });
        let (step_pipe, step_wall_pipe) = step_pipelines(
            &device,
            &step_layout,
            fp16,
            p16,
            wall_count > 0,
            surf.is_some(),
            [sa as u32 + keps as u32 + vles as u32, spec[0], spec[1], spec[2], sparse_bnx, wall_count as u32, n_active],
        );
        let init_pipe = pipe("d3q19_init", &init_layout, wgsl_init(fp16, p16), true);
        // ν̃ (ou k, ε) par briques sur les grilles boîtes à grands plans.
        let bricks = sparse.is_none() && nx * ny > SA_PLAN;
        let sa_pipe = if keps {
            let ke = KE3_WGSL.replace("SPARSE_FNS", SPARSE_FNS).replace("KE_FNS", KE_FNS);
            let ke = if vles { ke.replace("const VLES = false;", "const VLES = true;") } else { ke };
            pipe("ke3d", &sa_layout, format!("{PARAMS_WGSL}{}", bricked(ke, bricks)), true)
        } else {
            pipe("sa3d", &sa_layout, format!("{PARAMS_WGSL}{}", bricked(SA3_WGSL.replace("SPARSE_FNS", SPARSE_FNS), bricks)), true)
        };
        // Rang des mailles pariétales (collision des mailles pariétales, loi de paroi).
        let wcell = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("wall_cells3d"),
            contents: bytemuck::cast_slice(if wall_cells.is_empty() { &[0u32][..] } else { &wall_cells[..] }),
            usage: wgpu::BufferUsages::STORAGE,
        });
        // Loi de paroi : passe sur les mailles pariétales (avec Spalart–Allmaras).
        let wall_pass = (sa && wall_count > 0).then(|| {
            let code = format!(
                "{PARAMS_WGSL}{}",
                WALL3_WGSL.replace("SPARSE_FNS", SPARSE_FNS).replace("LAM_FNS", LAM_FNS).replace("SA_FNS", SA_FNS)
            );
            let pipe_w = pipe("wall3d", &wall_layout, code, true);
            pipe_w
        });
        let slice_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("slice"),
            entries: &[
                uniform(0, false, std::mem::size_of::<SliceParams>() as u64),
                storage(1, true),
                storage(2, true),
                storage(3, false),
            ],
        });
        let slice_pipe = pipe("slice", &slice_layout, SLICE_WGSL.to_string(), false);
        let auto = |label: &str, code: &str| {
            let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some(label),
                source: wgpu::ShaderSource::Wgsl(code.into()),
            });
            device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some(label),
                layout: None,
                module: &module,
                entry_point: Some("main"),
                compilation_options: Default::default(),
                cache: None,
            })
        };
        temps("pipelines calcul", tt);
        let derive_pipe = auto("derive3d", DERIVE_WGSL);
        let render_pipe = auto("render3d", RENDER_WGSL);
        temps("pipelines affichage", tt);
        let pb = Some(psize);
        let bind = (0..2)
            .map(|s| {
                device.create_bind_group(&wgpu::BindGroupDescriptor {
                    label: Some("d3q19"),
                    layout: &layout,
                    entries: &[
                        buf_entry(0, &params, pb),
                        buf_entry(1, &pops[s][0], None),
                        buf_entry(2, &pops[s][1], None),
                        buf_entry(3, &pops[1 - s][0], None),
                        buf_entry(4, &pops[1 - s][1], None),
                        buf_entry(5, &wall_buf, None),
                        buf_entry(6, &qbuf, None),
                        buf_entry(7, &force, None),
                        buf_entry(8, &macro_buf, None),
                        buf_entry(9, &nut[s], None),
                        buf_entry(10, &wall_geo, None),
                        buf_entry(11, &utau_buf, None),
                        buf_entry(12, &dist_buf, None),
                        buf_entry(13, &slab_buf, None),
                        buf_entry(14, &blk_buf, None),
                        buf_entry(15, &step_list_buf, None),
                        buf_entry(16, &plain_buf, None),
                    ],
                })
            })
            .collect();
        let sa_bind = (0..2)
            .map(|s| {
                device.create_bind_group(&wgpu::BindGroupDescriptor {
                    label: Some("sa3d"),
                    layout: &sa_layout,
                    entries: &[
                        buf_entry(0, &params, pb),
                        buf_entry(1, &macro_buf, None),
                        buf_entry(2, &wall_buf, None),
                        buf_entry(3, &wall_geo, None),
                        buf_entry(4, &utau_buf, None),
                        buf_entry(5, &dist_buf, None),
                        buf_entry(6, &nut[s], None),
                        buf_entry(7, &nut[1 - s], None),
                        buf_entry(8, &slab_buf, None),
                        buf_entry(9, &blk_buf, None),
                        buf_entry(10, &plain_buf, None),
                    ],
                })
            })
            .collect();
        let surf_pass = surf.map(|(st, law)| {
            let ns_ = st.area.len();
            let home = st.home_cells(&law);
            let mut geo: Vec<[f32; 4]> = Vec::with_capacity(8 * ns_.max(1));
            for k in 0..ns_ {
                let (n, c) = (st.normal[k], st.center[k]);
                geo.push([n[0] as f32, n[1] as f32, n[2] as f32, st.area[k] as f32]);
                geo.push([c[0] as f32, c[1] as f32, c[2] as f32, f32::from_bits(home[k])]);
                let m = crate::surfel3d::wall_law_matrix(n);
                for col in m {
                    geo.push(col.map(|v| v as f32));
                }
                for r in 0..3 {
                    geo.push([0, 1, 2, 3].map(|c| if 4 * r + c < 9 { st.scat_inv[9 * k + 4 * r + c] as f32 } else { 0.0 }));
                }
            }
            if geo.is_empty() {
                geo.resize(8, [0.0; 4]);
            }
            let mk_init = |label: &str, data: &[u8]| {
                device.create_buffer_init(&wgpu::util::BufferInitDescriptor { label: Some(label), contents: data, usage: wgpu::BufferUsages::STORAGE })
            };
            let geo_buf = mk_init("surfel_geo", bytemuck::cast_slice(&geo));
            let gst_buf = mk_init("surfel_gst", bytemuck::cast_slice(&st.gath_start));
            let mut gat: Vec<u32> = st.gath.iter().flat_map(|&(c, w)| [c, w.to_bits()]).collect();
            if gat.is_empty() {
                gat.resize(2, 0);
            }
            let gat_buf = mk_init("surfel_gat", bytemuck::cast_slice(&gat));
            let pst_buf = mk_init("surfel_pst", bytemuck::cast_slice(&st.probe_start));
            let mut prb: Vec<u32> = st.probe.iter().flat_map(|&(c, a, b)| [c, a.to_bits(), b.to_bits()]).collect();
            if prb.is_empty() {
                prb.resize(3, 0);
            }
            let prb_buf = mk_init("surfel_prb", bytemuck::cast_slice(&prb));
            let lay = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("surfel3d"),
                entries: &[
                    uniform(0, true, psize),
                    storage(1, true),
                    storage(2, true),
                    storage(3, false),
                    storage(4, true),
                    storage(5, true),
                    storage(6, true),
                    storage(7, true),
                    storage(8, true),
                    storage(9, true),
                    storage(10, true),
                    storage(11, true),
                ],
            });
            let pipe_s = pipe("surfel3d", &lay, wgsl_surfel(fp16, p16, ns_ as u32, surfels_o2()), true);
            let bind = (0..2)
                .map(|s| {
                    device.create_bind_group(&wgpu::BindGroupDescriptor {
                        label: Some("surfel3d"),
                        layout: &lay,
                        entries: &[
                            buf_entry(0, &params, pb),
                            buf_entry(1, &pops[s][0], None),
                            buf_entry(2, &pops[s][1], None),
                            buf_entry(3, &force, None),
                            buf_entry(4, &geo_buf, None),
                            buf_entry(5, &gst_buf, None),
                            buf_entry(6, &gat_buf, None),
                            buf_entry(7, &utau_buf, None),
                            buf_entry(8, &wall_geo, None),
                            buf_entry(9, &macro_buf, None),
                            buf_entry(10, &pst_buf, None),
                            buf_entry(11, &prb_buf, None),
                        ],
                    })
                })
                .collect();
            let groups = (ns_ as u32).div_ceil(64).max(1);
            let dx = groups.min(65535);
            let mut vol = vec![0f32; st.wall_cells.len()];
            for (r, &c) in st.wall_cells.iter().enumerate() {
                vol[r] = st.vol[c as usize];
            }
            Box::new(SurfPass {
                pipe: pipe_s,
                bind,
                disp: (dx, groups.div_ceil(dx)),
                n: ns_,
                center: st.center.iter().map(|c| c.map(|v| v as f32)).collect(),
                vol,
                stats: st.stats,
                gpu_bytes: st.gpu_bytes(),
                _bufs: [geo_buf, gst_buf, gat_buf, pst_buf, prb_buf],
            })
        });
        let wall_disp = {
            let groups = (wall_count as u32).div_ceil(64).max(1);
            let dx = groups.min(65535);
            (dx, groups.div_ceil(dx))
        };
        let wall_pass = wall_pass.map(|pipe_w| {
            let b = device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("wall3d"),
                layout: &wall_layout,
                entries: &[
                    buf_entry(0, &params, pb),
                    buf_entry(1, &macro_buf, None),
                    buf_entry(2, &wall_buf, None),
                    buf_entry(3, &wall_geo, None),
                    buf_entry(4, &utau_buf, None),
                    buf_entry(5, &dist_buf, None),
                    buf_entry(6, &slab_buf, None),
                    buf_entry(7, &blk_buf, None),
                    buf_entry(8, &wcell, None),
                ],
            });
            WallPass { pipe: pipe_w, bind: b, disp: wall_disp }
        });
        let init_bind = (0..2)
            .map(|s| {
                device.create_bind_group(&wgpu::BindGroupDescriptor {
                    label: Some("d3q19_init"),
                    layout: &init_layout,
                    entries: &[
                        buf_entry(0, &params, pb),
                        buf_entry(1, &pops[s][0], None),
                        buf_entry(2, &pops[s][1], None),
                    ],
                })
            })
            .collect();
        let groups = (ns as u32).div_ceil(64);
        let dx = groups.min(65535);
        let disp = (dx, groups.div_ceil(dx));
        let sa_disp = if bricks {
            let (by, bz) = (SA_BY, SA_BZ);
            let groups = (nx as u32 * (ny as u32).div_ceil(by) * by * (nz as u32).div_ceil(bz) * bz).div_ceil(64);
            let dx = groups.min(65535);
            (dx, groups.div_ceil(dx))
        } else {
            disp
        };
        let mut s = Gpu3D {
            device,
            queue,
            adapter_name,
            nx,
            ny,
            nz,
            n: ns,
            fp16,
            p16,
            u_in,
            tau0: 0.5,
            tau_sponge: 0.8,
            les: false,
            sponge_start: if nested { nx as u32 } else { (nx as f64 * 0.85).round() as u32 },
            cs: 0.15,
            pops,
            wall: wall_buf,
            sdf,
            sdf_geom,
            vox: vox_buf,
            _qbuf: qbuf,
            sa,
            keps,
            vles,
            nu,
            l: body.l,
            les_a: 1.0,
            xtr: 0.0,
            _slab: slab_buf,
            nut,
            _dist: dist_buf,
            _wall_geo: wall_geo,
            _utau: utau_buf,
            sa_pipe,
            sa_bind,
            wall_pass,
            nut0,
            force,
            force_read,
            params,
            macro_buf,
            step_pipe,
            step_wall_pipe,
            step_layout,
            spec,
            wall_disp,
            n_active,
            _step_list: step_list_buf,
            plain: plain_buf,
            disp_step: {
                let groups = (if n_active > 0 { n_active } else { ns as u32 }).div_ceil(step_wg());
                let dx = groups.min(65535);
                (dx, groups.div_ceil(dx))
            },
            _wcell: wcell,
            init_pipe,
            slice_pipe,
            derive_pipe,
            render_pipe,
            deriv: None,
            render_bufs: None,
            bind,
            init_bind,
            disp,
            sa_disp,
            cur: 0,
            time: 0,
            force_steps: 0,
            wall_count,
            solid_count,
            solid,
            ghost,
            lost: false,
            nested,
            ground: body.ground && !nested,
            periodic_span: false,
            periodic_x: false,
            essai: 0,
            essai_nu: 0.0,
            refine: None,
            wall_cells,
            band_fill,
            band_deep,
            band_fill_pos,
            sparse: sparse.map(std::sync::Arc::new),
            blk_buf,
            edge_cos,
            rest_force,
            surface: body.mesh.clone(),
            transp: None,
            tble: None,
            surf: surf_pass,
            accel_x: 0.0,
        };
        s.set_viscosity(nu);
        s.reset(0.0);
        temps("fin", tt);
        Ok(s)
    }

    /// Marque « recouvertes » (0xfffffffd) les mailles fluides de cette
    /// grille que les grilles fines recopient (listes de positions `lists`)
    /// et dont les six voisins sont aussi recopiés : la collision et le
    /// transport de ν̃ les sautent. Leur nouvel état n'était lu par
    /// personne avant la recopie qui l'écrase (le transport de ν̃ et la loi
    /// de paroi ne lisent que les six voisins et le champ d'avant le pas ;
    /// le tirage, l'état d'avant le pas) : résultats identiques au bit.
    pub(crate) fn mark_covered(&self, lists: &[&[u32]]) {
        let (nx, ny, nz) = (self.nx, self.ny, self.nz);
        let rank = |x: usize, y: usize, z: usize| match self.sparse.as_ref() {
            Some(m) => m.id(x, y, z),
            None => x + nx * (y + ny * z),
        };
        let mut inr = vec![false; self.n];
        for l in lists {
            for &p in l.iter() {
                let p = p as usize;
                inr[rank(p % nx, (p / nx) % ny, p / (nx * ny))] = true;
            }
        }
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("wall3d_read"),
            size: 4 * self.n as u64,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let mut wall: Vec<u32> = self.read(&self.wall, &staging, 4 * self.n as u64).iter().map(|v| v.to_bits()).collect();
        if wall.len() != self.n {
            return;
        }
        let mut count = 0usize;
        for l in lists {
            for &p in l.iter() {
                let p = p as usize;
                let (x, y, z) = (p % nx, (p / nx) % ny, p / (nx * ny));
                let r = rank(x, y, z);
                if wall[r] != 0 || x == 0 || y == 0 || z == 0 || x + 1 >= nx || y + 1 >= ny || z + 1 >= nz {
                    continue;
                }
                let nb = [
                    rank(x - 1, y, z),
                    rank(x + 1, y, z),
                    rank(x, y - 1, z),
                    rank(x, y + 1, z),
                    rank(x, y, z - 1),
                    rank(x, y, z + 1),
                ];
                if nb.iter().all(|&c| inr[c]) {
                    wall[r] = COVERED;
                    count += 1;
                }
            }
        }
        if std::env::var("SOUFFLERIE_TEMPS").is_ok() {
            eprintln!("  [mise en place] {count} mailles recouvertes sur {}", self.n);
        }
        self.queue.write_buffer(&self.wall, 0, bytemuck::cast_slice(&wall));
        self.queue.write_buffer(&self.plain, 0, bytemuck::cast_slice(&plain_mask(&wall)));
    }

    /// Variantes de diagnostic (enquête sur l'épaisseur de la couche limite,
    /// `plaque3d`), sur toutes les grilles ; 0 : calcul normal. Bits :
    /// 1 ν_t de Spalart–Allmaras retiré de la collision (sauf premier nœud) ;
    /// 2 sous-maille amortie près de la paroi (Van Driest) avec loi de paroi ;
    /// 4 premier nœud : viscosité de Spalart–Allmaras au lieu de κ u_τ h ;
    /// 8 glissement non relaxé ; 16 glissement libre (témoin non visqueux) ;
    /// 32 sans sous-maille ; 64 et 128 constante de Smagorinsky divisée
    /// par 2 et par 4 ; 256 viscosité de la couche pariétale κ u_τ y_n ;
    /// 512 glissement libre là où la paroi n'est pas horizontale (|n_z| < 0,999) ;
    /// 1024 viscosité de la couche pariétale plafonnée à `nu_cap` (unités de
    /// la veine) ; 2048 ν̃ de Spalart–Allmaras plafonné à `nu_cap` partout ;
    /// 4096 sans la borne de monotonie du glissement ; 8192 ancienne
    /// borne (u_h − pente · h) ; 16384 contrainte pariétale imposée comme force
    /// au nœud, paroi en glissement libre ; 32768 (avec 16384) viscosité du
    /// premier nœud nulle ; 65536 (avec 16) glissement libre à la vitesse
    /// relevée à h ; 131072 profil de la couche pariétale ancré à h seulement ;
    /// 262144 loi de paroi aussi aux arêtes vives ; 524288 ν̃ du premier nœud
    /// κ u_τ y_n (distance du nœud) ; 1048576 advection de ν̃ de Lax–Wendroff ;
    /// 4194304 sans le plancher du glissement (profil ancré à h) ;
    /// 8388608 loi de paroi avec gradient de pression longitudinal
    /// (Townsend–Mellor, `WALL3_WGSL`) ; 16777216 (avec 8388608) viscosité
    /// de la couche pariétale κ u_τ h gardée ; 33554432 (avec 8388608) u_τ
    /// de Spalding, gradient dans le seul profil de la couche pariétale ;
    /// 67108864 glissement signé, porté par la vitesse relevée à h ;
    /// 134217728 (avec 67108864) sans plancher ; 2097152 ν_t = max(ν_SA,
    /// ν_Smagorinsky) au lieu de la somme ; 536870912 (avec 8388608) gradient
    /// de pression filtré sur ±2 mailles le long de la paroi ; 1073741824
    /// (avec les deux précédents) gradient de pression totale à h.
    /// `MACH3_BIT` (2²⁸) correction du défaut cubique du D3Q19 (domaine de
    /// Mach, README).
    pub fn set_essai(&mut self, bits: u32, nu_cap: f64) {
        self.essai = bits;
        self.essai_nu = nu_cap as f32;
        if let Some(r) = self.refine.as_mut() {
            for (k, g) in r.levels.iter_mut().enumerate() {
                g.essai = bits;
                // Viscosité en unités de cette grille (2^niveau fois plus grande).
                g.essai_nu = (nu_cap * (1u64 << r.depth[k]) as f64) as f32;
            }
        }
    }

    /// Loi de paroi hors équilibre à équations de couche limite 1D (TBLE,
    /// `crate::tble`) sur toutes les grilles, au lieu de u_τ de Spalding ;
    /// None : retour à la passe par défaut. L'état des profils repart du
    /// profil d'équilibre au premier pas qui suit.
    pub fn set_paroi_tble(&mut self, mode: Option<crate::tble::ParoiTble>) {
        self.set_paroi_tble_one(mode);
        if let Some(r) = self.refine.as_mut() {
            for g in r.levels.iter_mut() {
                g.set_paroi_tble_one(mode);
            }
        }
    }

    /// Loi de paroi choisie (None : Spalding).
    pub fn paroi_tble(&self) -> Option<crate::tble::ParoiTble> {
        self.tble.as_ref().map(|t| t.0)
    }

    fn set_paroi_tble_one(&mut self, mode: Option<crate::tble::ParoiTble>) {
        use crate::tble;
        if self.wall_pass.is_none() {
            return;
        }
        let psize = std::mem::size_of::<Params>() as u64;
        let mut entries = vec![uniform(0, true, psize)];
        for (b, ro) in [(1, true), (2, true), (3, true), (4, false), (5, true), (6, true), (7, true), (8, true)] {
            entries.push(storage(b, ro));
        }
        let base = WALL3_WGSL.replace("SPARSE_FNS", SPARSE_FNS).replace("LAM_FNS", LAM_FNS);
        let code = match mode {
            None => {
                self.tble = None;
                format!("{PARAMS_WGSL}{}", base.replace("SA_FNS", SA_FNS))
            }
            Some(m) => {
                entries.push(storage(9, false));
                let ys = tble::grille(WM_H, tble::premier_point(self.nu, self.u_in as f64));
                let rep = |s: String, a: &str, b: &str| {
                    assert_eq!(s.matches(a).count(), 1, "TBLE : ancre « {a} » introuvable ou multiple");
                    s.replace(a, &format!("{a}{b}"))
                };
                let c = base.replace("SA_FNS", &format!("{SA_FNS}{}", tble::wgsl_fns(m, &ys)));
                let c = rep(c, "    var pg = 0.0;\n", tble::WGSL_DECL);
                let c = rep(c, "        yw = WM_H;\n", tble::WGSL_SOURCE);
                let c = rep(c, "    var ut2 = mellor_utau(ut, yw, select(pg, 0.0, (P.essai & 33554432u) != 0u));\n", tble::WGSL_UTAU);
                let size = (4 * tble::TB_N * self.wall_count.max(1)) as u64;
                let buf = self.device.create_buffer(&wgpu::BufferDescriptor {
                    label: Some("tble3d"),
                    size,
                    usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC,
                    mapped_at_creation: false,
                });
                self.tble = Some((m, buf));
                format!("{PARAMS_WGSL}{c}")
            }
        };
        let layout = self.device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor { label: Some("wall3d"), entries: &entries });
        let pl = self.device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("wall3d"),
            bind_group_layouts: &[Some(&layout)],
            immediate_size: 0,
        });
        let module = compute_module(&self.device, "wall3d", code);
        let pipe = self.device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("wall3d"),
            layout: Some(&pl),
            module: &module,
            entry_point: Some("main"),
            compilation_options: Default::default(),
            cache: None,
        });
        let pb = Some(psize);
        let mut be = vec![
            buf_entry(0, &self.params, pb),
            buf_entry(1, &self.macro_buf, None),
            buf_entry(2, &self.wall, None),
            buf_entry(3, &self._wall_geo, None),
            buf_entry(4, &self._utau, None),
            buf_entry(5, &self._dist, None),
            buf_entry(6, &self._slab, None),
            buf_entry(7, &self.blk_buf, None),
            buf_entry(8, &self._wcell, None),
        ];
        if let Some((_, b)) = &self.tble {
            be.push(buf_entry(9, b, None));
        }
        let bind = self.device.create_bind_group(&wgpu::BindGroupDescriptor { label: Some("wall3d"), layout: &layout, entries: &be });
        let disp = self.wall_disp;
        self.wall_pass = Some(WallPass { pipe, bind, disp });
    }

    /// Accélération uniforme constante selon x (force volumique ρ a ;
    /// validation : écoulement de Poiseuille en veine périodique).
    pub fn set_accel(&mut self, a: f64) {
        self.accel_x = a as f32;
    }

    /// Veine périodique dans les trois directions (cas de validation sans
    /// maquette : onde de cisaillement, tourbillon convecté). Grille boîte
    /// sans raffinement ; ni entrée, ni sortie, ni éponge.
    pub fn set_periodic_x(&mut self) {
        self.periodic_x = true;
        self.sponge_start = self.nx as u32;
    }

    /// Champ initial donné par maille (u_x, u_y, u_z, ρ), rangé x + nx (y +
    /// ny z) : populations à l'équilibre de la collision (Hermite d'ordre 3
    /// partiel). Grille boîte en simple précision, sans raffinement.
    pub fn set_state(&mut self, m: &[[f32; 4]]) -> Result<(), String> {
        if self.fp16 || self.sparse.is_some() || self.refine.is_some() || m.len() != self.n {
            return Err("set_state : grille boîte simple précision sans raffinement, une valeur par maille".into());
        }
        let n = self.n;
        let mut e = vec![0f32; 19 * n];
        for (c, &[ux, uy, uz, rho]) in m.iter().enumerate() {
            let (ux, uy, uz, rho) = (ux as f64, uy as f64, uz as f64, rho as f64);
            let usq = 1.5 * (ux * ux + uy * uy + uz * uz);
            let a = [
                rho * ux * ux * uy,
                rho * ux * uy * uy,
                rho * ux * ux * uz,
                rho * ux * uz * uz,
                rho * uy * uy * uz,
                rho * uy * uz * uz,
            ];
            for i in 0..19 {
                let [cx, cy, cz] = C3[i].map(|v| v as f64);
                let w = w3(i);
                let cu = cx * ux + cy * uy + cz * uz;
                let h = |p: f64, q: f64| p * (q * q - 1.0 / 3.0);
                let h3 = h(cy, cx) * a[0] + h(cx, cy) * a[1] + h(cz, cx) * a[2] + h(cx, cz) * a[3] + h(cz, cy) * a[4] + h(cy, cz) * a[5];
                e[i * n + c] = (w * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - usq) + 13.5 * w * h3 - w) as f32;
            }
        }
        let a: Vec<u32> = e[..10 * n].iter().map(|v| v.to_bits()).collect();
        let b: Vec<u32> = e[10 * n..].iter().map(|v| v.to_bits()).collect();
        self.queue.write_buffer(&self.pops[self.cur][0], 0, bytemuck::cast_slice(&a));
        self.queue.write_buffer(&self.pops[self.cur][1], 0, bytemuck::cast_slice(&b));
        Ok(())
    }

    pub fn adapter_name(&self) -> &str {
        &self.adapter_name
    }

    pub fn set_viscosity(&mut self, nu: f64) {
        self.set_viscosity_one(nu);
        if let Some(r) = self.refine.as_mut() {
            for (k, g) in r.levels.iter_mut().enumerate() {
                g.set_viscosity_one(nu * (1u64 << r.depth[k]) as f64);
                g.les = self.les;
            }
        }
    }

    fn set_viscosity_one(&mut self, nu: f64) {
        self.nu = nu;
        // Longueur d'amortissement pariétal : 8 L/√Re mailles, bornée à [1, 4] (comme en 2D).
        let re = self.u_in as f64 * self.l as f64 / nu;
        self.les_a = (8.0 * self.l as f64 / re.sqrt()).clamp(1.0, 4.0);
        self.tau0 = 0.5 + 3.0 * nu;
        self.tau_sponge = self.tau0.max(0.8);
        // Sous-maille seulement aux Reynolds de maille élevés : sans amortissement
        // pariétal, elle épaissirait la couche limite et gonflerait la traînée.
        self.les = self.tau0 < 0.51;
    }

    pub fn is_solid(&self, x: usize, y: usize, z: usize) -> bool {
        self.solid[x + self.nx * (y + self.ny * z)]
    }

    pub(crate) fn params_for(&self, u: f32, accel: f32, write_macro: bool) -> Params {
        self.params(u, accel, write_macro)
    }

    fn params(&self, u: f32, accel: f32, write_macro: bool) -> Params {
        Params {
            nx: self.nx as u32,
            ny: self.ny as u32,
            nz: self.nz as u32,
            n: self.n as u32,
            tau0: self.tau0 as f32,
            tau_sponge: self.tau_sponge as f32,
            sponge_start: self.sponge_start,
            u_in: u,
            accel: accel + self.accel_x,
            c18: (18.0 * self.cs * self.cs) as f32,
            write_macro: write_macro as u32,
            les: self.les as u32,
            rows: self.nested as u32 | ((self.ground as u32) << 1) | ((self.periodic_span as u32) << 2) | ((self.periodic_x as u32) << 3),
            sa: self.sa as u32 + self.keps as u32 + self.vles as u32,
            nu: self.nu as f32,
            les_a: self.les_a as f32,
            xtr: self.xtr as f32,
            band: (self.l as f32 / 10.0).max(3.0),
            essai: self.essai,
            essai_nu: self.essai_nu,
            bnx: self.sparse.as_ref().map_or(0, |m| m.bd[0] as u32),
            bny: self.sparse.as_ref().map_or(0, |m| m.bd[1] as u32),
            bnz: self.sparse.as_ref().map_or(0, |m| m.bd[2] as u32),
            bcoord: self.sparse.as_ref().map_or(0, |m| m.table.len() as u32),
        }
    }

    /// Fluide au repos (ou à la vitesse u0).
    pub fn reset(&mut self, u0: f32) {
        self.reset_one(u0);
        if let Some(r) = self.refine.as_mut() {
            for g in r.levels.iter_mut() {
                g.reset_one(u0);
            }
            r.fresh = true;
        }
    }

    fn reset_one(&mut self, u0: f32) {
        self.queue.write_buffer(&self.params, 0, bytemuck::bytes_of(&self.params(u0, 0.0, false)));
        let mut enc = self.device.create_command_encoder(&Default::default());
        {
            let mut pass = enc.begin_compute_pass(&Default::default());
            pass.set_pipeline(&self.init_pipe);
            for b in &self.init_bind {
                pass.set_bind_group(0, b, &[0]);
                pass.dispatch_workgroups(self.disp.0, self.disp.1, 1);
            }
        }
        enc.clear_buffer(&self.force, 0, None);
        enc.clear_buffer(&self.macro_buf, 0, None);
        if let Some((_, b)) = &self.tble {
            enc.clear_buffer(b, 0, None);
        }
        self.queue.submit([enc.finish()]);
        if self.sa {
            for b in &self.nut {
                self.queue.write_buffer(b, 0, bytemuck::cast_slice(&self.nut0));
            }
        }
        self.cur = 0;
        self.time = 0;
        self.force_steps = 0;
    }

    /// `n` pas ; la vitesse d'entrée suit `ramp`. Le champ macroscopique est
    /// écrit au dernier pas (pour les coupes).
    pub fn run(&mut self, n: usize, ramp: &Ramp) {
        self.ensure_spec();
        if let Some(mut r) = self.refine.take() {
            crate::refine3d::run(self, &mut r, n, ramp);
            self.refine = Some(r);
            return;
        }
        // Soumissions courtes (≈ 3·10⁷ mises à jour de maille) : sous Windows, une
        // soumission de plus de 2 s fait réinitialiser le pilote, ce qui peut
        // arriver sur un GPU lent ou bridé en température.
        let per_submit = (30_000_000 / self.n.max(1)).clamp(1, BATCH);
        // Soumissions en vol bornées : sans attente, le processeur empile
        // toutes les soumissions (une par pas au-delà de 1,5·10⁷ mailles) et
        // la mémoire de l'hôte s'épuise (≈ 2 Mo par soumission en attente à
        // 2·10⁷ mailles : 7 Go après 3 t*, puis échec d'allocation).
        let mut en_vol = std::collections::VecDeque::new();
        let mut left = n;
        while left > 0 {
            let k = left.min(per_submit);
            let mut bytes = vec![0u8; SLOT as usize * k];
            for s in 0..k {
                let (u, a) = ramp.at(self.time + s as u64);
                let p = self.params(u, a, left == k && s == k - 1);
                let o = s * SLOT as usize;
                bytes[o..o + std::mem::size_of::<Params>()].copy_from_slice(bytemuck::bytes_of(&p));
            }
            self.queue.write_buffer(&self.params, 0, &bytes);
            let mut enc = self.device.create_command_encoder(&Default::default());
            {
                let mut pass = enc.begin_compute_pass(&Default::default());
                for s in 0..k {
                    let off = (s as u64 * SLOT) as u32;
                    if let Some(w) = &self.wall_pass {
                        w.encode(&mut pass, off);
                    }
                    self.encode_step(&mut pass, off);
                    if self.sa {
                        pass.set_pipeline(&self.sa_pipe);
                        pass.set_bind_group(0, &self.sa_bind[self.cur], &[off]);
                        pass.dispatch_workgroups(self.sa_disp.0, self.sa_disp.1, 1);
                    }
                    self.cur = 1 - self.cur;
                }
            }
            en_vol.push_back(self.queue.submit([enc.finish()]));
            if en_vol.len() > EN_VOL_MAX {
                let _ = self.device.poll(wgpu::PollType::Wait {
                    submission_index: en_vol.pop_front(),
                    timeout: None,
                });
            }
            self.time += k as u64;
            self.force_steps += k as u64;
            left -= k;
        }
    }

    fn read(&self, src: &wgpu::Buffer, staging: &wgpu::Buffer, size: u64) -> Vec<f32> {
        self.read_at(src, 0, staging, size)
    }

    fn read_at(&self, src: &wgpu::Buffer, offset: u64, staging: &wgpu::Buffer, size: u64) -> Vec<f32> {
        let mut enc = self.device.create_command_encoder(&Default::default());
        enc.copy_buffer_to_buffer(src, offset, staging, 0, size);
        self.queue.submit([enc.finish()]);
        let (tx, rx) = mpsc::channel();
        staging.map_async(wgpu::MapMode::Read, 0..size, move |r| {
            let _ = tx.send(r);
        });
        let _ = self.device.poll(wgpu::PollType::wait_indefinitely());
        let ok = rx.recv().map(|r| r.is_ok()).unwrap_or(false);
        let out = if ok {
            let view = staging.get_mapped_range(0..size).expect("lecture GPU");
            bytemuck::cast_slice::<u8, f32>(&view).to_vec()
        } else {
            Vec::new()
        };
        staging.unmap();
        out
    }

    /// Force moyenne par pas depuis l'appel précédent (lecture synchrone, attend
    /// la fin du calcul), puis remise à zéro de l'accumulateur.
    /// Efforts moyens par pas ; avec des grilles fines, ceux de la plus fine,
    /// ramenés aux unités de la veine (surfaces 4ᴷ fois plus grandes).
    pub fn take_forces(&mut self) -> Force3 {
        if let Some(mut r) = self.refine.take() {
            // Grille la plus fine qui entoure toute la maquette (celle du
            // bord d'attaque n'en couvre qu'une partie).
            let k = r.full as i32;
            let _ = self.take_forces_one();
            let s = 4f64.powi(k);
            if r.levels.len() == r.full {
                let f = r.levels[r.full - 1].take_forces_one();
                self.lost |= r.levels.iter().any(|g| g.lost);
                self.refine = Some(r);
                return Force3 {
                    fx: f.fx / s,
                    fy: f.fy / s,
                    fz: f.fz / s,
                    steps: f.steps >> k,
                };
            }
            // Grilles partielles (bords d'attaque et de fuite) : leurs mailles
            // pariétales dans leur domaine propre ; la grille entière ailleurs
            // (sous une grille fine, elle reçoit l'état recopié).
            let (fo, fs) = (r.origin[r.full - 1], r.scale(r.full - 1));
            let cores: Vec<[[f64; 3]; 2]> = r.core.iter().flatten().copied().collect();
            let in_core = |q: [f64; 3], c: &[[f64; 3]; 2]| (0..3).all(|a| q[a] >= c[0][a] && q[a] < c[1][a]);
            let a = r.levels[r.full - 1].take_forces_masked(|p| {
                let q = [0, 1, 2].map(|i| fo[i] + p[i] / fs);
                !cores.iter().any(|c| in_core(q, c))
            });
            let mut f = [a.fx / s, a.fy / s, a.fz / s];
            for k in r.full..r.levels.len() {
                let (lo, ls) = (r.origin[k], r.scale(k));
                let core = r.core[k];
                // Domaines des grilles plus fines (couches imbriquées) : comptés là-bas.
                let deeper: Vec<[[f64; 3]; 2]> =
                    (r.full..r.levels.len()).filter(|&j| r.depth[j] > r.depth[k]).filter_map(|j| r.core[j]).collect();
                let b = r.levels[k].take_forces_masked(|p| {
                    let q = [0, 1, 2].map(|i| lo[i] + p[i] / ls);
                    core.is_some_and(|c| in_core(q, &c)) && !deeper.iter().any(|c| in_core(q, c))
                });
                let sb = ls * ls;
                f[0] += b.fx / sb;
                f[1] += b.fy / sb;
                f[2] += b.fz / sb;
            }
            self.lost |= r.levels.iter().any(|g| g.lost);
            self.refine = Some(r);
            return Force3 {
                fx: f[0],
                fy: f[1],
                fz: f[2],
                steps: a.steps >> k,
            };
        }
        self.take_forces_one()
    }

    /// Couche limite laminaire imposée du bord d'attaque jusqu'à la fraction
    /// `xtr` de la corde locale (étendue en x du corps dans chaque tranche
    /// d'envergure) : frottement de Blasius à la paroi, pas de production de
    /// Spalart–Allmaras autour. 0 : turbulente dès l'avant.
    pub fn set_transition(&mut self, xtr: f64) {
        self.xtr = xtr.clamp(0.0, 1.0);
        if let Some(r) = self.refine.as_mut() {
            for g in r.levels.iter_mut() {
                g.xtr = self.xtr;
            }
        }
    }

    /// Grilles fines : (dimensions, mailles pariétales) de la plus fine, et facteur.
    /// Champs (u, ρ) de la grille la plus fine, dimensions, origine (mailles
    /// de la veine) et rapport de raffinement (diagnostic).
    pub fn finest_macros(&mut self) -> (Vec<[f32; 4]>, [usize; 3], [f64; 3], f64) {
        match self.refine.as_mut() {
            Some(r) => {
                let k = r.full;
                let org = r.origin[k - 1];
                let g = &mut r.levels[k - 1];
                let d = [g.nx, g.ny, g.nz];
                (g.macros(), d, org, (1u64 << k) as f64)
            }
            None => {
                let d = [self.nx, self.ny, self.nz];
                (self.macros(), d, [0.0; 3], 1.0)
            }
        }
    }

    /// ν̃ de Spalart–Allmaras par maille (vide sans SA) ; avec k-ε RNG,
    /// ν_t = C_μ k²/ε (que `nu_turb` laisse presque inchangé loin de la paroi).
    pub fn nu_tilde(&self) -> Vec<f32> {
        if !self.sa {
            return Vec::new();
        }
        if self.vles {
            let ns = self.nut[self.cur].size() / 12;
            let size = 4 * self.n as u64;
            let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("nut3d_read"),
                size,
                usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });
            return self.read_at(&self.nut[self.cur], 8 * ns, &staging, size);
        }
        if self.keps {
            return self.k_eps().iter().map(|&[k, e]| crate::sa::nu_ke(k, e)).collect();
        }
        let size = 4 * self.n as u64;
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("nut3d_read"),
            size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        self.read(&self.nut[self.cur], &staging, size)
    }

    /// (k, ε) de k-ε RNG par maille (vide sans ce modèle).
    pub fn k_eps(&self) -> Vec<[f32; 2]> {
        if !self.keps {
            return Vec::new();
        }
        let size = 8 * self.n as u64;
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("ke3d_read"),
            size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        self.read(&self.nut[self.cur], &staging, size).chunks_exact(2).map(|c| [c[0], c[1]]).collect()
    }

    /// Traînée en champ lointain (Betz, Maskell) : bilan de quantité de
    /// mouvement entre deux plans transversaux x = `x_in` (amont) et
    /// x = `x_out` (aval), en mailles de la veine, sur toute la section
    /// (bords latéraux périodiques). Chaque maille de la veine du plan est
    /// intégrée sur la grille la plus fine qui la couvre (trapèzes).
    /// Composantes de la traînée (Betz, Maskell ; le bilan complet entre les
    /// deux plans n'est pas exploitable ici : les bandes absorbantes des bords
    /// latéraux exercent leurs propres efforts, et les ondes de pression de la
    /// veine le font osciller) :
    /// - perte de pression totale ∫ (p₀,∞ − p₀) dS, p₀,∞ pris dans le même
    ///   plan hors du sillage, séparée entre le
    ///   sillage visqueux (ν̃ > 10 ν, fluide passé dans la couche limite) et le
    ///   reste de l'écoulement, qui devrait être sans perte : ce reste est la
    ///   traînée parasite de la discrétisation (dissipation numérique) ;
    /// - traînée induite : énergie de l'écoulement transversal ½ ∫ ρ(v² + w²) dS
    ///   (moins celle du plan amont).
    /// Unités du réseau de la veine (effort = coefficient × ½ ρ U² S).
    pub fn wake_drag(&mut self, x_in: usize, x_out: &[usize]) -> Vec<WakeDrag> {
        struct Lev {
            mac: Vec<[f32; 4]>,
            nut: Vec<f32>,
            d: [usize; 3],
            org: [f64; 3],
            s: f64,
            nu: f64,
            m: usize,
            map: Option<std::sync::Arc<SparseMap>>,
        }
        let mut levs = vec![Lev {
            mac: self.macros(),
            nut: self.nu_tilde(),
            d: [self.nx, self.ny, self.nz],
            org: [0.0; 3],
            s: 1.0,
            nu: self.nu,
            m: 0,
            map: None,
        }];
        if let Some(r) = self.refine.as_mut() {
            for (k, g) in r.levels.iter_mut().enumerate() {
                levs.push(Lev {
                    mac: g.macros(),
                    nut: g.nu_tilde(),
                    d: [g.nx, g.ny, g.nz],
                    org: r.origin[k],
                    s: (1u64 << r.depth[k]) as f64,
                    nu: g.nu,
                    // Enveloppe fantôme et zone de recouvrement exclues ; une
                    // couche pariétale ne couvre pas le plan (marge infinie).
                    m: if g.sparse.is_some() { usize::MAX / 4 } else { 2 * crate::refine3d::MARGIN3 + 2 },
                    map: g.sparse.clone(),
                });
            }
        }
        let (ny, nz) = (self.ny, self.nz);
        // Intégrales sur un plan : [∫(p + ρu²), ∫p₀, ∫ν̃-sillage p₀, aire sillage,
        // ∫½ρ(v² + w²), ∫ρu, aire, ∫ρuw, ∫ρuv].
        let plane = |x: usize| -> ([f64; 9], Vec<f64>) {
            let mut acc = [0.0f64; 9];
            let mut free = Vec::new();
            for z in 0..nz {
                for y in 0..ny {
                    // Grille la plus fine qui couvre la maille (y ± ½, z ± ½).
                    let mut pick = 0;
                    for (k, l) in levs.iter().enumerate().skip(1) {
                        let fx = (x as f64 - l.org[0]) * l.s;
                        let lo = [(y as f64 - 0.5 - l.org[1]) * l.s, (z as f64 - 0.5 - l.org[2]) * l.s];
                        let hi = [(y as f64 + 0.5 - l.org[1]) * l.s, (z as f64 + 0.5 - l.org[2]) * l.s];
                        let m = l.m as f64;
                        let inside = |v: f64, n: usize| v >= m && v <= (n - 1) as f64 - m;
                        if inside(fx, l.d[0]) && inside(lo[0], l.d[1]) && inside(hi[0], l.d[1]) && inside(lo[1], l.d[2]) && inside(hi[1], l.d[2]) {
                            pick = k;
                        }
                    }
                    let l = &levs[pick];
                    let half = (l.s / 2.0).round() as i64;
                    let (cx, cy, cz) = (
                        ((x as f64 - l.org[0]) * l.s).round() as i64,
                        ((y as f64 - l.org[1]) * l.s).round() as i64,
                        ((z as f64 - l.org[2]) * l.s).round() as i64,
                    );
                    let (j0, j1) = if pick == 0 { (0, 0) } else { (-half, half) };
                    for dj in j0..=j1 {
                        for dk in j0..=j1 {
                            let wy = if pick > 0 && dj.abs() == half { 0.5 } else { 1.0 };
                            let wz = if pick > 0 && dk.abs() == half { 0.5 } else { 1.0 };
                            let w = wy * wz / (l.s * l.s);
                            let (qx, qy, qz) = (cx as usize, (cy + dj) as usize, (cz + dk) as usize);
                            let c = match l.map.as_ref() {
                                Some(m) => m.id(qx, qy, qz),
                                None => qx + l.d[0] * (qy + l.d[1] * qz),
                            };
                            let v = l.mac[c];
                            let (u, vv, ww, rho) = (v[0] as f64, v[1] as f64, v[2] as f64, v[3] as f64);
                            let p = (rho - 1.0) / 3.0;
                            let p0 = p + 0.5 * rho * (u * u + vv * vv + ww * ww);
                            acc[0] += w * (p + rho * u * u);
                            acc[1] += w * p0;
                            let wake = l.nut.get(c).is_some_and(|&t| t as f64 > 10.0 * l.nu);
                            if wake {
                                acc[2] += w * p0;
                                acc[3] += w;
                            } else if dj == 0 && dk == 0 {
                                free.push(p0);
                            }
                            acc[4] += w * 0.5 * rho * (vv * vv + ww * ww);
                            acc[5] += w * rho * u;
                            acc[6] += w;
                            acc[7] += w * rho * u * ww;
                            acc[8] += w * rho * u * vv;
                        }
                    }
                }
            }
            (acc, free)
        };
        let (a, _) = plane(x_in);
        x_out
            .iter()
            .map(|&x| {
                let (b, mut free) = plane(x);
                // Pression totale de l'écoulement non perturbé dans ce plan
                // (médiane hors du sillage) : les ondes de pression de la
                // veine, uniformes sur la section, s'y retrouvent et s'annulent.
                free.sort_by(f64::total_cmp);
                let p0_ref = free.get(free.len() / 2).copied().unwrap_or(0.0);
                let loss = p0_ref * b[6] - b[1];
                let loss_wake = p0_ref * b[3] - b[2];
                WakeDrag {
                    loss_wake,
                    loss_outside: loss - loss_wake,
                    induced: b[4] - a[4],
                    mass_ratio: b[5] / a[5],
                    wake_area: b[3],
                }
            })
            .collect()
    }

    /// Efforts sur la maquette par bilan de quantité de mouvement sur une
    /// boîte qui l'entoure à `margin` longueur de référence, intégré sur la
    /// grille la plus fine qui entoure toute la maquette : −∮ (ρ u u·n + p n) dA.
    /// Indépendant du traitement des parois ; une onde de pression uniforme
    /// s'annule entre faces opposées. Unités du réseau de la veine.
    pub fn box_forces(&mut self, margin: f64) -> Option<[f64; 3]> {
        self.box_forces_vol(&[margin]).pop().flatten().map(|(f, _)| f)
    }

    /// Idem pour plusieurs marges sur un même champ : efforts et volume de
    /// la boîte (mailles de la veine).
    pub fn box_forces_vol(&mut self, margins: &[f64]) -> Vec<Option<([f64; 3], f64)>> {
        let (g, org, s) = match self.refine.as_mut() {
            Some(r) => {
                let k = r.full;
                (&mut r.levels[k - 1], r.origin[k - 1], (1u64 << k) as f64)
            }
            None => (&mut *self, [0.0; 3], 1.0),
        };
        let Some(mesh) = g.surface.clone() else { return vec![None; margins.len()] };
        let mac = g.macros();
        let d = [g.nx, g.ny, g.nz];
        let gl = g.l as f64;
        let _ = org;
        margins.iter().map(|&margin| {
        let mg = margin * gl;
        let mut lo = [0usize; 3];
        let mut hi = [0usize; 3];
        for a in 0..3 {
            let l = (mesh.lo[a] - mg).floor();
            let h = (mesh.hi[a] + mg).ceil();
            if l < 3.0 || h > (d[a] - 4) as f64 {
                return None;
            }
            lo[a] = l as usize;
            hi[a] = h as usize;
        }
        let at = |x: usize, y: usize, z: usize| mac[x + d[0] * (y + d[1] * z)];
        let mut flux = [0.0f64; 3];
        for (ax, sgn, pos) in [(0, -1.0, lo[0]), (0, 1.0, hi[0]), (1, -1.0, lo[1]), (1, 1.0, hi[1]), (2, -1.0, lo[2]), (2, 1.0, hi[2])] {
            let (a1, a2) = ((ax + 1) % 3, (ax + 2) % 3);
            for i in lo[a1]..=hi[a1] {
                for j in lo[a2]..=hi[a2] {
                    // Trapèzes : demi-poids sur les arêtes de la face.
                    let w = (if i == lo[a1] || i == hi[a1] { 0.5 } else { 1.0 }) * (if j == lo[a2] || j == hi[a2] { 0.5 } else { 1.0 });
                    let mut p3 = [0usize; 3];
                    p3[ax] = pos;
                    p3[a1] = i;
                    p3[a2] = j;
                    let c = at(p3[0], p3[1], p3[2]);
                    let (u, rho) = ([c[0] as f64, c[1] as f64, c[2] as f64], c[3] as f64);
                    let un = sgn * u[ax];
                    for k in 0..3 {
                        flux[k] += w * (rho * u[k] * un + if k == ax { sgn * (rho - 1.0) / 3.0 } else { 0.0 });
                    }
                }
            }
        }
        let vol: f64 = (0..3).map(|a| (hi[a] - lo[a]) as f64 / s).product();
        Some((flux.map(|f| -f / (s * s)), vol))
        }).collect()
    }

    /// Mailles pariétales de la grille la plus fine (diagnostic) : position
    /// (mailles de la veine) et plus petit cosinus entre les faces touchées.
    /// Grille d'indice `k` (None : la veine, Some(i) : niveau i du raffinement).
    fn grid_mut(&mut self, k: Option<usize>) -> &mut Gpu3D {
        match k {
            None => self,
            Some(i) => &mut self.refine.as_mut().expect("raffinement").levels[i],
        }
    }

    /// Grille entière la plus fine (indice pour `grid_mut`), son origine
    /// (mailles de la veine) et son facteur.
    fn full_grid(&self) -> (Option<usize>, [f64; 3], f64) {
        match self.refine.as_ref() {
            Some(r) => (Some(r.full - 1), r.origin[r.full - 1], (1u64 << r.full) as f64),
            None => (None, [0.0; 3], 1.0),
        }
    }

    /// Décambrage visqueux : prépare la couche limite intégrale sur le
    /// maillage de la maquette (grille entière la plus fine) et le triangle
    /// de chaque maille pariétale, sur cette grille et les grilles partielles
    /// plus fines (bord d'attaque). Les grilles plus grossières, recouvertes
    /// par la grille fine autour de la maquette, n'en reçoivent pas. `xtr` :
    /// transition imposée (fraction de corde ; 0 : turbulent dès le point
    /// d'arrêt). Erreur sans maillage ou sans loi de paroi. Rend le nombre
    /// de triangles.
    pub fn enable_transpiration(&mut self, xtr: f64, resolved: bool) -> Result<usize, String> {
        if !self.sa {
            return Err("la transpiration demande la loi de paroi (Spalart–Allmaras)".into());
        }
        let (kf, org_f, s_f) = self.full_grid();
        let mesh = self.grid_mut(kf).surface.clone().ok_or("maquette sans maillage")?;
        let topo = crate::ibl3d::Topology::new(&mesh.tris);
        let mut grids = vec![kf];
        if let Some(r) = self.refine.as_ref() {
            grids.extend((r.full..r.levels.len()).map(Some));
        }
        let mut maps = Vec::new();
        for k in grids {
            let (org, s) = match (k, self.refine.as_ref()) {
                (Some(i), Some(r)) => (r.origin[i], r.scale(i)),
                _ => ([0.0; 3], 1.0),
            };
            let g = self.grid_mut(k);
            // Position dans la grille entière la plus fine.
            let tri: Vec<u32> = g
                .wall_cells
                .iter()
                .map(|&c| {
                    let p = g.cell_xyz(c as usize).map(|v| v as f64);
                    let q = [0, 1, 2].map(|a| (org[a] + p[a] / s - org_f[a]) * s_f);
                    topo.nearest(q, 4.0 * s_f / s).map_or(u32::MAX, |t| t as u32)
                })
                .collect();
            let n = tri.len();
            maps.push((k, tri, vec![0f32; n]));
        }
        let nt = topo.c.len();
        self.transp = Some(Box::new(Transpiration { topo, maps, xtr, resolved }));
        Ok(nt)
    }

    /// Un couplage visqueux–non visqueux : vitesse extérieure relevée à la
    /// hauteur de raccord de la loi de paroi (WM_H) le long de la normale de
    /// chaque triangle, couche limite intégrale (`ibl3d::march`), vitesse de
    /// transpiration v_n = ∇·(u_e δ* ê) (Lighthill 1958) imposée aux mailles
    /// pariétales, sous-relaxée : v_n ← v_n + ω (v_n,calculée − v_n). None si
    /// `enable_transpiration` n'a pas été appelée.
    pub fn couple_transpiration(&mut self, relax: f64) -> Option<TranspirationReport> {
        let mut tr = self.transp.take()?;
        let (kf, _, _) = self.full_grid();
        let u_in = self.u_in as f64;
        let g = self.grid_mut(kf);
        let mac = g.macros();
        let (d, nu, lref) = ([g.nx, g.ny, g.nz], g.nu, g.l as f64);
        let solid = &g.solid;
        let sample = |p: [f64; 3]| -> Option<[f64; 3]> {
            let i = p.map(|v| v.floor() as i64);
            let w = [0, 1, 2].map(|a| p[a] - i[a] as f64);
            let mut out = [0.0; 3];
            for c in 0..8 {
                let o = [c & 1, (c >> 1) & 1, (c >> 2) & 1];
                let q = [0, 1, 2].map(|a| i[a] + o[a] as i64);
                if (0..3).any(|a| q[a] < 0 || q[a] >= d[a] as i64) {
                    return None;
                }
                let c = q[0] as usize + d[0] * (q[1] as usize + d[1] * q[2] as usize);
                if solid[c] {
                    return None;
                }
                let ww: f64 = (0..3).map(|a| if o[a] == 1 { w[a] } else { 1.0 - w[a] }).product();
                for k in 0..3 {
                    out[k] += ww * mac[c][k] as f64;
                }
            }
            Some(out)
        };
        let t = &tr.topo;
        let nt = t.c.len();
        // Profil le long de la normale, de 1 à 8 mailles par demi-maille :
        // vitesse à WM_H (plus loin si un nœud voisin est solide), bord de la
        // couche déficitaire résolue (premier maximum de u_t au-dessus de
        // WM_H) et son épaisseur de déplacement, la couche [0, WM_H] de la loi
        // de paroi prise linéaire (vitesse à la paroi extrapolée de 1 et 2).
        let prof: Vec<Option<([f64; 3], [f64; 3], f64, f64)>> = (0..nt)
            .map(|k| {
                if t.area[k] <= 0.0 {
                    return None;
                }
                let (c, n) = (t.c[k], t.n[k]);
                let pts: Vec<(f64, Option<[f64; 3]>)> = (2..=16)
                    .map(|i| {
                        let h = 0.5 * i as f64;
                        let v = sample([0, 1, 2].map(|a| c[a] + h * n[a])).map(|v| {
                            let vn = v[0] * n[0] + v[1] * n[1] + v[2] * n[2];
                            [0, 1, 2].map(|a| v[a] - vn * n[a])
                        });
                        (h, v)
                    })
                    .collect();
                let i0 = pts.iter().position(|(h, v)| *h >= WM_H && v.is_some())?;
                if pts[i0].0 > WM_H + 2.0 {
                    return None;
                }
                let uwm = pts[i0].1?;
                let mag = |v: [f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                let mut ie = i0;
                while ie + 1 < pts.len() {
                    match pts[ie + 1].1 {
                        Some(v) if mag(v) > mag(pts[ie].1.unwrap()) => ie += 1,
                        _ => break,
                    }
                }
                let (he, uev) = (pts[ie].0, pts[ie].1.unwrap());
                let ue = mag(uev).max(1e-12);
                // Points du profil jusqu'à h_e, vitesse à la paroi extrapolée.
                let mut pr: Vec<(f64, f64)> = pts[..=ie].iter().filter_map(|(h, v)| v.map(|v| (*h, mag(v)))).collect();
                let u0 = match (pr.iter().find(|p| p.0 == 1.0), pr.iter().find(|p| p.0 == 2.0)) {
                    (Some(a), Some(b)) => (2.0 * a.1 - b.1).max(0.0),
                    _ => pr[0].1,
                };
                pr.insert(0, (0.0, u0));
                let ds: f64 = pr.windows(2).map(|w| 0.5 * (w[1].0 - w[0].0) * (2.0 - (w[0].1 + w[1].1) / ue)).sum();
                Some((uwm, uev, he, ds.max(0.0)))
            })
            .collect();
        let resolved = tr.resolved;
        let mut ue: Vec<Option<[f64; 3]>> = prof.iter().map(|p| p.map(|p| if resolved { p.1 } else { p.0 })).collect();
        let dstar_res: Vec<f64> = prof.iter().map(|p| p.map_or(0.0, |p| p.3)).collect();
        let h_edge: Vec<f64> = prof.iter().map(|p| p.map_or(0.0, |p| p.2)).collect();
        let missed_area: f64 = (0..nt).filter(|&k| ue[k].is_none()).map(|k| t.area[k]).sum();
        // Relevés manquants : moyenne des voisins relevés, projetée dans le plan.
        for _ in 0..20 {
            let prev = ue.clone();
            let mut left = false;
            for k in 0..nt {
                if prev[k].is_some() || t.area[k] <= 0.0 {
                    continue;
                }
                let (mut s, mut m) = ([0.0; 3], 0.0);
                for e in &t.edges[k] {
                    if let Some(v) = prev.get(e.nb as usize).copied().flatten() {
                        s = [0, 1, 2].map(|a| s[a] + v[a]);
                        m += 1.0;
                    }
                }
                if m > 0.0 {
                    let n = t.n[k];
                    let vn = (s[0] * n[0] + s[1] * n[1] + s[2] * n[2]) / m;
                    ue[k] = Some([0, 1, 2].map(|a| s[a] / m - vn * n[a]));
                } else {
                    left = true;
                }
            }
            if !left {
                break;
            }
        }
        let ue: Vec<[f64; 3]> = ue.into_iter().map(|v| v.unwrap_or([0.0; 3])).collect();
        let layer = crate::ibl3d::march(t, &ue, &crate::ibl3d::Options { nu, xtr: tr.xtr });
        // Déplacement imposé : δ* entier, ou ce qui manque au champ résolu.
        let vn_t = if resolved {
            let m: Vec<f64> = (0..nt)
                .map(|k| {
                    let u = (ue[k][0] * ue[k][0] + ue[k][1] * ue[k][1] + ue[k][2] * ue[k][2]).sqrt();
                    u * (layer.dstar[k] - dstar_res[k])
                })
                .collect();
            crate::ibl3d::transpiration(t, &ue, &m).0
        } else {
            layer.vn.clone()
        };
        // Moyenne par l'aire sur le triangle et ses voisins : même débit
        // injecté, mais un triangle effilé (aire presque nulle) ne donne plus
        // une vitesse démesurée à la maille qui le touche.
        let vn_t: Vec<f64> = (0..nt)
            .map(|k| {
                let (mut q, mut ar) = (vn_t[k] * t.area[k], t.area[k]);
                for e in &t.edges[k] {
                    if let (Some(&v), Some(&a)) = (vn_t.get(e.nb as usize), t.area.get(e.nb as usize)) {
                        q += v * a;
                        ar += a;
                    }
                }
                if ar > 0.0 { q / ar } else { 0.0 }
            })
            .collect();
        // Mailles pariétales : v_n du triangle, sous-relaxée, envoyée au GPU.
        let maps = std::mem::take(&mut tr.maps);
        let mut out_maps = Vec::with_capacity(maps.len());
        for (k, tri, mut vn) in maps {
            for (v, &ti) in vn.iter_mut().zip(&tri) {
                let target = if ti == u32::MAX { 0.0 } else { vn_t[ti as usize] as f32 };
                *v += relax as f32 * (target - *v);
            }
            let g = self.grid_mut(k);
            if g.wall_count > 0 {
                g.queue.write_buffer(&g._utau, (16 * g.wall_count) as u64, bytemuck::cast_slice(&vn));
            }
            out_maps.push((k, tri, vn));
        }
        tr.maps = out_maps;
        let t = &tr.topo;
        let area: f64 = t.area.iter().sum();
        let (mut sb, mut ab, mut vmax, mut sep) = (0.0, 0.0, 0.0f64, 0.0);
        for k in 0..nt {
            if vn_t[k] > 0.0 {
                sb += vn_t[k] * t.area[k];
                ab += t.area[k];
            }
            vmax = vmax.max(vn_t[k].abs());
            if layer.separated[k] {
                sep += t.area[k];
            }
        }
        let rep = TranspirationReport {
            missed: missed_area / area.max(1e-12),
            vn_mean: sb / ab.max(1e-12) / u_in,
            vn_max: vmax / u_in,
            separated: sep / area.max(1e-12),
            outflow: layer.outflow / (u_in * lref * lref),
            l: lref,
            area: t.area.clone(),
            centers: t.c.clone(),
            normals: t.n.clone(),
            dstar_res,
            h_edge,
            vn: vn_t,
            ue,
            layer,
        };
        self.transp = Some(tr);
        Some(rep)
    }

    /// ν̃ de Spalart–Allmaras de la grille entière la plus fine (vide sans SA).
    pub fn finest_nu_tilde(&self) -> Vec<f32> {
        match self.refine.as_ref() {
            Some(r) => r.levels[r.full - 1].nu_tilde(),
            None => self.nu_tilde(),
        }
    }

    /// Diagnostic : par maille pariétale de la grille entière la plus fine
    /// (ordre de `finest_wall_cells`), u_τ et vitesse de glissement relaxée.
    pub fn finest_wall_state(&self) -> Vec<[f32; 4]> {
        let g = match self.refine.as_ref() {
            Some(r) => &r.levels[r.full - 1],
            None => self,
        };
        if g.wall_count == 0 || !g.sa {
            return Vec::new();
        }
        let size = 16 * g.wall_count as u64;
        let staging = g.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("utau_read"),
            size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        g.read(&g._utau, &staging, size).chunks(4).map(|c| [c[0], c[1], c[2], c[3]]).collect()
    }

    pub fn finest_wall_cells(&self) -> Vec<([f64; 3], f32)> {
        let (g, org, s) = match self.refine.as_ref() {
            Some(r) => {
                let k = r.full;
                (&r.levels[k - 1], r.origin[k - 1], (1u64 << k) as f64)
            }
            None => (self, [0.0; 3], 1.0),
        };
        g.wall_cells
            .iter()
            .zip(&g.edge_cos)
            .map(|(&c, &e)| {
                let c = c as usize;
                let p = [(c % g.nx) as f64, ((c / g.nx) % g.ny) as f64, (c / (g.nx * g.ny)) as f64];
                ([org[0] + p[0] / s, org[1] + p[1] / s, org[2] + p[2] / s], e)
            })
            .collect()
    }

    /// Position de la maille de rang `c` (origine de la grille si `c` est
    /// dans le bloc vide).
    pub(crate) fn cell_xyz(&self, c: usize) -> [usize; 3] {
        match self.sparse.as_ref() {
            Some(m) => m.pos(c, [self.nx, self.ny, self.nz]).unwrap_or([0; 3]),
            None => [c % self.nx, (c / self.nx) % self.ny, c / (self.nx * self.ny)],
        }
    }

    /// Grilles partielles des bords, une ligne par boîte (diagnostic) :
    /// origine, étendue et domaine propre en mailles de la veine.
    pub fn edge_box_report(&self) -> Vec<String> {
        let Some(r) = self.refine.as_ref() else { return Vec::new() };
        (r.full..r.levels.len())
            .map(|k| {
                let (g, o, s) = (&r.levels[k], r.origin[k], r.scale(k));
                let c = r.core[k].unwrap_or([[0.0; 3]; 2]);
                format!(
                    "x {:.1}..{:.1} y {:.1}..{:.1} z {:.1}..{:.1} · propre x {:.1}..{:.1} y {:.1}..{:.1} z {:.1}..{:.1} · {} mailles pariétales",
                    o[0], o[0] + (g.nx - 1) as f64 / s, o[1], o[1] + (g.ny - 1) as f64 / s, o[2], o[2] + (g.nz - 1) as f64 / s,
                    c[0][0], c[1][0], c[0][1], c[1][1], c[0][2], c[1][2], g.wall_count
                )
            })
            .collect()
    }

    /// Grilles partielles des bords : nombre de boîtes, mailles, mailles
    /// pariétales, facteur.
    pub fn le_grid(&self) -> Option<(usize, usize, usize, usize)> {
        let r = self.refine.as_ref()?;
        let parts = &r.levels[r.full..];
        if parts.is_empty() {
            return None;
        }
        let n = parts.iter().map(|g| g.n).sum();
        let w = parts.iter().map(|g| g.wall_count).sum();
        Some((parts.len(), n, w, r.scale(r.levels.len() - 1) as usize))
    }

    pub fn finest(&self) -> Option<([usize; 3], usize, usize)> {
        let r = self.refine.as_ref()?;
        let g = r.levels.get(r.full.checked_sub(1)?)?;
        Some(([g.nx, g.ny, g.nz], g.wall_count, 1 << r.full))
    }

    /// Mailles calculées de chaque grille (veine, puis grilles fines) et
    /// leur niveau (banc de mesure : mises à jour de maille par pas de la
    /// veine, Σ n_k 2^niveau).
    pub fn grid_cells(&self) -> Vec<(usize, u32)> {
        let mut v = vec![(self.n, 0)];
        if let Some(r) = self.refine.as_ref() {
            v.extend(r.levels.iter().zip(&r.depth).map(|(g, &d)| (g.n, d)));
        }
        v
    }

    /// Profilage (banc de mesure) : durée moyenne de chaque passe, mesurée
    /// isolément sur `reps` lancements (pas de collision, transport de ν̃,
    /// raccords des grilles). Modifie l'état : à appeler après les relevés.
    pub fn profile_passes(&mut self, reps: usize) -> Vec<(String, usize, f64)> {
        let time = |g: &Gpu3D, f: &dyn Fn(&mut wgpu::ComputePass)| {
            let _ = g.device.poll(wgpu::PollType::wait_indefinitely());
            let t0 = std::time::Instant::now();
            let mut enc = g.device.create_command_encoder(&Default::default());
            {
                let mut pass = enc.begin_compute_pass(&Default::default());
                for _ in 0..reps {
                    f(&mut pass);
                }
            }
            g.queue.submit([enc.finish()]);
            let _ = g.device.poll(wgpu::PollType::wait_indefinitely());
            t0.elapsed().as_secs_f64() * 1e3 / reps as f64
        };
        let one = |g: &Gpu3D, name: &str, out: &mut Vec<(String, usize, f64)>| {
            if let Some(w) = &g.wall_pass {
                out.push((format!("{name} loi de paroi"), g.wall_count, time(g, &|p| w.encode(p, 0))));
            }
            out.push((
                format!("{name} collision"),
                g.n,
                time(g, &|p| {
                    p.set_pipeline(&g.step_pipe);
                    p.set_bind_group(0, &g.bind[g.cur], &[0]);
                    p.dispatch_workgroups(g.disp_step.0, g.disp_step.1, 1);
                }),
            ));
            if let Some(w) = &g.step_wall_pipe {
                out.push((
                    format!("{name} collision pariétale"),
                    g.wall_count,
                    time(g, &|p| {
                        p.set_pipeline(w);
                        p.set_bind_group(0, &g.bind[g.cur], &[0]);
                        p.dispatch_workgroups(g.wall_disp.0, g.wall_disp.1, 1);
                    }),
                ));
            }
            if g.sa {
                out.push((
                    format!("{name} ν̃"),
                    g.n,
                    time(g, &|p| {
                        p.set_pipeline(&g.sa_pipe);
                        p.set_bind_group(0, &g.sa_bind[g.cur], &[0]);
                        p.dispatch_workgroups(g.sa_disp.0, g.sa_disp.1, 1);
                    }),
                ));
            }
        };
        self.ensure_spec();
        if let Some(r) = self.refine.as_mut() {
            for g in r.levels.iter_mut() {
                g.ensure_spec();
            }
        }
        let mut out = Vec::new();
        one(self, "veine", &mut out);
        if let Some(r) = self.refine.as_ref() {
            for (k, g) in r.levels.iter().enumerate() {
                one(g, &format!("grille {} (×{})", k + 1, 1u64 << r.depth[k]), &mut out);
            }
            for (k, link) in r.links.iter().enumerate() {
                let p = if r.parent[k] == 0 { &*self } else { &r.levels[r.parent[k] - 1] };
                let c = &r.levels[k];
                for (name, n, ms) in link.profile(p, c, &time) {
                    out.push((format!("raccord {} {name}", k + 1), n, ms));
                }
            }
        }
        out
    }

    /// Efforts sur une maquette maillée par intégration sur ses triangles
    /// (champ du dernier pas ; chaque triangle sur la grille la plus fine qui
    /// l'entoure) : −p_p n dA, p_p extrapolée linéairement jusqu'à la paroi
    /// depuis les relevés à ½ et 1 maille le long de la normale, interpolés
    /// sur les seuls nœuds fluides de ce côté de la paroi (relevée plus loin,
    /// là où tout le voisinage est fluide, la dépression du bord d'attaque
    /// était écrêtée et la portance sous-estimée de 10 % par rapport au bilan
    /// de quantité de mouvement), plus τ_p dA, τ_p = ρ u_τ² de la loi de
    /// Spalding dans la direction de la vitesse tangentielle relevée à WM_H.
    /// Plus juste que l'échange de quantité de mouvement pour la traînée : sur
    /// l'escalier des mailles, les liens presque parallèles à la paroi
    /// renvoient la pression et la quantité de mouvement longitudinale dans
    /// leur propre direction, d'où une traînée parasite sur les marches.
    /// None sans maillage.
    pub fn surface_forces(&mut self) -> Option<SurfaceForces> {
        // Grilles utilisées : la plus fine qui entoure toute la maquette, puis
        // les grilles partielles (bord d'attaque), de la moins à la plus fine.
        struct Lev {
            mac: Vec<[f32; 4]>,
            solid: Vec<bool>,
            d: [usize; 3],
            org: [f64; 3],
            s: f64,
            nu: f64,
            mesh: std::sync::Arc<crate::stl::MeshBody>,
            /// Bande exclue près des bords (fantômes, recouvrement) ; 0 : grille entière.
            m: f64,
            /// Grille partielle : domaine propre (mailles de la veine).
            core: Option<[[f64; 3]; 2]>,
            map: Option<std::sync::Arc<SparseMap>>,
        }
        let mut levs: Vec<Lev> = Vec::new();
        match self.refine.as_mut() {
            Some(r) => {
                for k in r.full..=r.levels.len() {
                    let sk = r.scale(k - 1);
                    let g = &mut r.levels[k - 1];
                    let mesh = g.surface.clone()?;
                    // Fantômes d'une couche creuse : ni calculées ni relevées
                    // (couche limitée à des zones ou des bords : la grille
                    // parente relève les triangles hors de la couche).
                    let solid = if g.ghost.is_empty() { g.solid.clone() } else { g.solid.iter().zip(&g.ghost).map(|(&a, &b)| a || b).collect() };
                    levs.push(Lev {
                        mac: g.macros(),
                        solid,
                        d: [g.nx, g.ny, g.nz],
                        org: r.origin[k - 1],
                        s: sk,
                        nu: g.nu,
                        mesh,
                        m: if k == r.full { 0.0 } else { (2 * crate::refine3d::MARGIN3 + 2) as f64 },
                        core: r.core[k - 1],
                        map: g.sparse.clone(),
                    });
                }
            }
            None => {
                let mesh = self.surface.clone()?;
                levs.push(Lev {
                    mac: self.macros(),
                    solid: self.solid.clone(),
                    d: [self.nx, self.ny, self.nz],
                    org: [0.0; 3],
                    s: 1.0,
                    nu: self.nu,
                    mesh,
                    m: 0.0,
                    core: None,
                    map: None,
                });
            }
        }
        let normal = |t: &crate::stl::Tri| {
            let e1 = [t[1][0] - t[0][0], t[1][1] - t[0][1], t[1][2] - t[0][2]];
            let e2 = [t[2][0] - t[0][0], t[2][1] - t[0][1], t[2][2] - t[0][2]];
            [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]]
        };
        // Normales sortantes si le volume signé est positif.
        let vol: f64 = levs[0]
            .mesh
            .tris
            .iter()
            .map(|t| {
                let n = normal(t);
                (0..3).map(|k| t[0][k] * n[k]).sum::<f64>()
            })
            .sum();
        let sgn = if vol >= 0.0 { 1.0 } else { -1.0 };
        let outside = |l: &Lev, q: [i64; 3]| (0..3).any(|a| (q[a] as f64) < l.m || q[a] as f64 > (l.d[a] - 1) as f64 - l.m);
        let sample = |l: &Lev, p: [f64; 3]| -> Option<[f64; 4]> {
            let i = p.map(|v| v.floor() as i64);
            let w = [0, 1, 2].map(|a| p[a] - i[a] as f64);
            let mut out = [0.0; 4];
            for c in 0..8 {
                let o = [c & 1, (c >> 1) & 1, (c >> 2) & 1];
                let q = [i[0] + o[0] as i64, i[1] + o[1] as i64, i[2] + o[2] as i64];
                if outside(l, q) {
                    return None;
                }
                let c = match l.map.as_ref() {
                    Some(m) => m.id(q[0] as usize, q[1] as usize, q[2] as usize),
                    None => q[0] as usize + l.d[0] * (q[1] as usize + l.d[1] * q[2] as usize),
                };
                if l.solid[c] {
                    return None;
                }
                let ww: f64 = (0..3).map(|a| if o[a] == 1 { w[a] } else { 1.0 - w[a] }).product();
                for k in 0..4 {
                    out[k] += ww * l.mac[c][k] as f64;
                }
            }
            Some(out)
        };
        // Pression au plus près de la paroi : interpolation trilinéaire sur
        // les seuls nœuds fluides du côté de la normale (pas ceux de l'autre
        // face d'une paroi mince), poids renormalisés.
        // Hauteurs des deux relevés ; SOUFFLERIE_PWH="h1,h2" les change
        // (diagnostic, exemple pwh3d ; h2 ≤ h1 : ancien relevé seul).
        let (h1, h2) = std::env::var("SOUFFLERIE_PWH")
            .ok()
            .and_then(|v| {
                let mut it = v.split(',').map(|x| x.parse::<f64>().ok());
                Some((it.next()??, it.next()??))
            })
            .unwrap_or((PW_H1, PW_H2));
        let near = |l: &Lev, c: [f64; 3], n: [f64; 3], h: f64| -> Option<[f64; 4]> {
            let p = [0, 1, 2].map(|k| c[k] + h * n[k]);
            let i = p.map(|v| v.floor() as i64);
            let w = [0, 1, 2].map(|a| p[a] - i[a] as f64);
            let (mut acc, mut ws) = ([0.0f64; 4], 0.0);
            for k in 0..8 {
                let o = [k & 1, (k >> 1) & 1, (k >> 2) & 1];
                let q = [i[0] + o[0] as i64, i[1] + o[1] as i64, i[2] + o[2] as i64];
                if outside(l, q) {
                    return None;
                }
                let cell = match l.map.as_ref() {
                    Some(m) => m.id(q[0] as usize, q[1] as usize, q[2] as usize),
                    None => q[0] as usize + l.d[0] * (q[1] as usize + l.d[1] * q[2] as usize),
                };
                let side: f64 = (0..3).map(|a| (q[a] as f64 - c[a]) * n[a]).sum();
                if l.solid[cell] || side <= 0.0 {
                    continue;
                }
                let ww: f64 = (0..3).map(|a| if o[a] == 1 { w[a] } else { 1.0 - w[a] }).product();
                let m = l.mac[cell];
                for a in 0..3 {
                    acc[a] += ww * m[a] as f64;
                }
                acc[3] += ww * (m[3] as f64 - 1.0) / 3.0;
                ws += ww;
            }
            (ws >= 0.25).then(|| acc.map(|v| v / ws))
        };
        // Relevés d'un triangle sur une grille : p à h, p à la paroi, vitesse
        // et hauteur du frottement, relevé à h ; None si la grille ne
        // l'entoure pas.
        let probe = |l: &Lev, c: [f64; 3], n: [f64; 3]| -> Option<(f64, f64, [f64; 4], f64, [f64; 4])> {
            // Une maille le long de la normale ; plus loin si un nœud voisin est solide.
            let mut h = 1.0;
            let mut v = None;
            while v.is_none() && h <= 3.0 {
                v = sample(l, [0, 1, 2].map(|k| c[k] + h * n[k]));
                if v.is_none() {
                    h += 0.5;
                }
            }
            let v = v?;
            let p = (v[3] - 1.0) / 3.0;
            // Pression extrapolée linéairement jusqu'à la paroi depuis les
            // relevés à ½ et 1 maille (nœuds fluides de ce côté seulement) ;
            // à défaut, depuis h et h + 1.
            let mut pw = match sample(l, [0, 1, 2].map(|k| c[k] + (h + 1.0) * n[k])) {
                Some(v2) => p + h * (p - (v2[3] - 1.0) / 3.0),
                None => p,
            };
            if let (true, Some(v1), Some(v2)) = (h2 > h1, near(l, c, n, h1), near(l, c, n, h2)) {
                pw = v1[3] + h1 * (v1[3] - v2[3]) / (h2 - h1);
            }
            // Frottement : loi de paroi alimentée à la hauteur de raccord du
            // calcul (WM_H mailles), comme le glissement.
            let (vf, hf) = match sample(l, [0, 1, 2].map(|k| c[k] + WM_H.max(h) * n[k])) {
                Some(vf) => (vf, WM_H.max(h)),
                None => (v, h),
            };
            Some((p, pw, vf, hf, v))
        };
        let mut out = SurfaceForces::default();
        let mut area = 0.0;
        let nl = levs.len();
        for ti in 0..levs[0].mesh.tris.len() {
            // Grille la plus fine qui entoure le triangle.
            let mut got = None;
            for li in (0..nl).rev() {
                let l = &levs[li];
                let t = &l.mesh.tris[ti];
                let n0 = normal(t);
                let a2 = (n0[0] * n0[0] + n0[1] * n0[1] + n0[2] * n0[2]).sqrt();
                if a2 < 1e-12 {
                    break;
                }
                let n = n0.map(|v| sgn * v / a2);
                let c = [0, 1, 2].map(|k| (t[0][k] + t[1][k] + t[2][k]) / 3.0);
                // Grille partielle : triangle entier loin de ses bords, centre
                // dans son domaine propre (une seule des grilles sœurs).
                if l.m > 0.0 && t.iter().any(|v| (0..3).any(|a| v[a] < l.m + 3.0 || v[a] > (l.d[a] - 1) as f64 - l.m - 3.0)) {
                    continue;
                }
                if let Some(cr) = l.core {
                    let q = [0, 1, 2].map(|a| l.org[a] + c[a] / l.s);
                    if !(0..3).all(|a| q[a] >= cr[0][a] && q[a] < cr[1][a]) {
                        continue;
                    }
                }
                let ar = 0.5 * a2;
                match probe(l, c, n) {
                    Some(r) => {
                        got = Some((li, n, c, ar, Some(r)));
                        break;
                    }
                    None if li == 0 => got = Some((li, n, c, ar, None)),
                    None => {}
                }
            }
            let Some((li, _, c0, ar0, r0)) = got else { continue };
            let l = &levs[li];
            // Surfaces ramenées aux mailles de la veine.
            let s2 = l.s * l.s;
            area += ar0 / s2;
            if r0.is_none() {
                out.missed += ar0 / s2;
                continue;
            }
            // Relevés : le triangle est découpé en m² sous-triangles de côté
            // au plus une maille de cette grille (m ≤ 8), chacun relevé en son
            // centre le long de la normale lisse (normales des coins
            // interpolées), au lieu d'un seul relevé au centre le long de la
            // normale plate de la facette.
            let t = l.mesh.tris[ti];
            let edge = |a: usize, b: usize| (0..3).map(|k| (t[a][k] - t[b][k]).powi(2)).sum::<f64>().sqrt();
            let m = edge(0, 1).max(edge(1, 2)).max(edge(2, 0)).ceil().clamp(1.0, 8.0) as usize;
            let sub_ar = ar0 / (m * m) as f64 / s2;
            let mf = m as f64;
            // Centres des sous-triangles en coordonnées barycentriques.
            let mut centers: Vec<[f64; 3]> = Vec::with_capacity(m * m);
            for i in 0..m {
                for j in 0..m - i {
                    let (fi, fj) = (i as f64, j as f64);
                    let up = [fi + 1.0 / 3.0, fj + 1.0 / 3.0];
                    centers.push([1.0 - (up[0] + up[1]) / mf, up[0] / mf, up[1] / mf]);
                    if i + j + 1 < m {
                        let dn = [fi + 2.0 / 3.0, fj + 2.0 / 3.0];
                        centers.push([1.0 - (dn[0] + dn[1]) / mf, dn[0] / mf, dn[1] / mf]);
                    }
                }
            }
            let mut fl = [0.0; 3];
            let mut sh = [0.0; 3];
            let (mut pw_acc, mut n_acc, mut got_any) = (0.0, [0.0; 3], false);
            for w in &centers {
                let c = [0, 1, 2].map(|k| w[0] * t[0][k] + w[1] * t[1][k] + w[2] * t[2][k]);
                let n = l.mesh.smooth_normal_bary(ti, *w).map(|v| sgn * v);
                let Some((p, pw, vf, hf, v)) = probe(l, c, n) else {
                    out.missed += sub_ar;
                    continue;
                };
                got_any = true;
                let ar = sub_ar;
                let un = vf[0] * n[0] + vf[1] * n[1] + vf[2] * n[2];
                let ut = [0, 1, 2].map(|k| vf[k] - un * n[k]);
                let utm = (ut[0] * ut[0] + ut[1] * ut[1] + ut[2] * ut[2]).sqrt();
                let mut utau = spalding_utau(utm, hf, l.nu);
                if self.essai & 8_388_608 != 0 && utm > 0.0 {
                    // Gradient de pression le long de la vitesse tangentielle
                    // relevée (différence centrée sur ±1 maille), comme la
                    // loi de paroi du calcul.
                    let at = |sg: f64| sample(l, [0, 1, 2].map(|k| c[k] + hf * n[k] + sg * ut[k] / utm));
                    if let (Some(a), Some(b)) = (at(1.0), at(-1.0)) {
                        let pg = (a[3] - b[3]) / (6.0 * v[3]);
                        utau = mellor_utau(utm, hf, l.nu, pg);
                    }
                }
                let tau = v[3] * utau * utau;
                for k in 0..3 {
                    out.pressure[k] -= p * n[k] * ar;
                    out.pressure_wall[k] -= pw * n[k] * ar;
                    fl[k] -= pw * n[k] * ar;
                    if utm > 0.0 {
                        out.friction[k] += tau * ut[k] / utm * ar;
                        fl[k] += tau * ut[k] / utm * ar;
                        sh[k] += tau * ut[k] / utm * ar;
                    }
                    n_acc[k] += n[k] * ar;
                }
                pw_acc += pw * ar;
            }
            if !got_any {
                continue;
            }
            let nl_ = (n_acc[0] * n_acc[0] + n_acc[1] * n_acc[1] + n_acc[2] * n_acc[2]).sqrt().max(1e-300);
            let ar_tot = sub_ar * centers.len() as f64;
            out.loads.push(([0, 1, 2].map(|k| l.org[k] + c0[k] / l.s), fl));
            out.walls.push((n_acc.map(|v| v / nl_), pw_acc / ar_tot));
            let mut ys = l.mesh.tris[ti].map(|v| l.org[1] + v[1] / l.s);
            ys.sort_by(f64::total_cmp);
            out.spans.push(ys);
            out.tri.push(ti);
            out.shear.push(sh.map(|v| v / ar_tot));
            out.wm_h.push(WM_H / l.s);
        }
        out.missed /= area.max(1e-12);
        Some(out)
    }

    /// Efforts par maille pariétale depuis le dernier relevé (grille la plus
    /// fine), moyennés par pas et ramenés aux unités de la veine.
    pub fn take_wall_forces(&mut self) -> WallForces {
        let (g, org, s) = match self.refine.as_mut() {
            Some(r) => {
                let k = r.full;
                let org = r.origin[k - 1];
                (&mut r.levels[k - 1], org, (1u64 << k) as f64)
            }
            None => (&mut *self, [0.0; 3], 1.0),
        };
        let steps = std::mem::take(&mut g.force_steps).max(1) as f64;
        let mut out = WallForces { pos: Vec::new(), force: Vec::new() };
        if g.n_force() == 0 {
            return out;
        }
        let v = g.read(&g.force, &g.force_read, 16 * g.n_force() as u64);
        let mut enc = g.device.create_command_encoder(&Default::default());
        enc.clear_buffer(&g.force, 0, None);
        g.queue.submit([enc.finish()]);
        for (k, c) in v.chunks_exact(4).enumerate() {
            let p = g.force_pos(k);
            out.pos.push([0, 1, 2].map(|a| org[a] + p[a] / s));
            let f = s * s * steps;
            let r = g.rest_force[k];
            out.force.push([0, 1, 2].map(|a| (c[a] as f64 - steps * r[a] as f64) / f));
        }
        if self.refine.is_some() {
            self.force_steps = 0;
        }
        out
    }

    /// Force moyenne par pas des seules mailles pariétales dont la position
    /// (mailles de cette grille) satisfait `keep`, puis remise à zéro.
    fn take_forces_masked(&mut self, keep: impl Fn([f64; 3]) -> bool) -> Force3 {
        let steps = std::mem::take(&mut self.force_steps);
        if steps == 0 || self.n_force() == 0 {
            return Force3 { steps, ..Default::default() };
        }
        let v = self.read(&self.force, &self.force_read, 16 * self.n_force() as u64);
        if v.is_empty() {
            self.lost = true;
            return Force3 { fx: f64::NAN, fy: f64::NAN, fz: f64::NAN, steps };
        }
        let mut enc = self.device.create_command_encoder(&Default::default());
        enc.clear_buffer(&self.force, 0, None);
        self.queue.submit([enc.finish()]);
        let mut f = [0.0f64; 3];
        for (k, c) in v.chunks_exact(4).enumerate() {
            if keep(self.force_pos(k)) {
                // Sans la part de la pression de référence (nulle seulement
                // en somme sur une surface fermée, pas sur une partie).
                let r = self.rest_force[k];
                for a in 0..3 {
                    f[a] += c[a] as f64 - steps as f64 * r[a] as f64;
                }
            }
        }
        let s = steps as f64;
        Force3 { fx: f[0] / s, fy: f[1] / s, fz: f[2] / s, steps }
    }

    /// Accumulateurs d'efforts : mailles pariétales, ou surfels.
    fn n_force(&self) -> usize {
        self.surf.as_ref().map_or(self.wall_count, |s| s.n)
    }

    /// Position (mailles de cette grille) de l'accumulateur `k`.
    fn force_pos(&self, k: usize) -> [f64; 3] {
        match &self.surf {
            Some(s) => s.center[k].map(|v| v as f64),
            None => self.cell_xyz(self.wall_cells[k] as usize).map(|v| v as f64),
        }
    }

    /// Paroi volumique : (surfels, mailles réceptrices, écarts géométriques,
    /// octets des tables sur le GPU) de la grille la plus fine.
    pub fn surfel_info(&self) -> Option<(crate::surfel3d::Stats, usize)> {
        let g = match &self.refine {
            Some(r) => r.levels.last()?,
            None => self,
        };
        g.surf.as_ref().map(|s| (s.stats, s.gpu_bytes))
    }

    /// Masse du fluide Σ V ρ (grille boîte sans raffinement ; volumes
    /// fluides des voxels coupés avec les surfels, 1 sinon, 0 au solide).
    pub fn fluid_mass(&mut self) -> f64 {
        let m = self.macros();
        let vol = self.cell_volumes();
        m.iter().zip(&vol).map(|(q, &v)| v as f64 * q[3] as f64).sum()
    }

    /// Volume fluide de chaque maille de cette grille : 0 solide, fraction
    /// des voxels coupés avec les surfels, 1 sinon.
    pub fn cell_volumes(&self) -> Vec<f32> {
        let mut vol = vec![1f32; self.n];
        for (c, s) in self.solid.iter().enumerate() {
            if *s {
                vol[c] = 0.0;
            }
        }
        if let Some(sp) = &self.surf {
            for (r, &c) in self.wall_cells.iter().enumerate() {
                vol[c as usize] = sp.vol[r];
            }
        }
        vol
    }

    fn take_forces_one(&mut self) -> Force3 {
        let steps = std::mem::take(&mut self.force_steps);
        if steps == 0 || self.n_force() == 0 {
            let _ = self.device.poll(wgpu::PollType::wait_indefinitely());
            return Force3 {
                steps,
                ..Default::default()
            };
        }
        let v = self.read(&self.force, &self.force_read, 16 * self.n_force() as u64);
        if v.is_empty() {
            // Lecture impossible : GPU perdu (pilote réinitialisé).
            self.lost = true;
            return Force3 {
                fx: f64::NAN,
                fy: f64::NAN,
                fz: f64::NAN,
                steps,
            };
        }
        let mut enc = self.device.create_command_encoder(&Default::default());
        enc.clear_buffer(&self.force, 0, None);
        self.queue.submit([enc.finish()]);
        let (mut fx, mut fy, mut fz) = (0.0f64, 0.0f64, 0.0f64);
        for c in v.chunks_exact(4) {
            fx += c[0] as f64;
            fy += c[1] as f64;
            fz += c[2] as f64;
        }
        let s = steps as f64;
        Force3 {
            fx: fx / s,
            fy: fy / s,
            fz: fz / s,
            steps,
        }
    }

    /// Coupe normale à l'axe `axis` (0 : x, 1 : y, 2 : z) en `pos` (maille).
    /// Plan x : (y, z) ; plan y : (x, z) ; plan z : (x, y).
    /// Champ macroscopique de la veine (u_x, u_y, u_z, ρ) par maille, tel
    /// qu'écrit au dernier pas (sous la zone raffinée : recopie de la grille fine).
    pub fn macros(&mut self) -> Vec<[f32; 4]> {
        let size = 16 * self.n as u64;
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("macro3d_read"),
            size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let v = self.read(&self.macro_buf, &staging, size);
        v.chunks_exact(4).map(|c| [c[0], c[1], c[2], c[3]]).collect()
    }

    pub fn slice(&mut self, axis: usize, pos: usize, field: SliceField, l: f64) -> Slice {
        let (w, h) = match axis {
            0 => (self.ny, self.nz),
            1 => (self.nx, self.nz),
            _ => (self.nx, self.ny),
        };
        let lim = [self.nx, self.ny, self.nz][axis.min(2)];
        let p = SliceParams {
            nx: self.nx as u32,
            ny: self.ny as u32,
            nz: self.nz as u32,
            axis: axis.min(2) as u32,
            pos: pos.min(lim - 1) as u32,
            field: match field {
                SliceField::Speed => 0,
                SliceField::Vorticity => 1,
                SliceField::Q => 2,
            },
            w: w as u32,
            h: h as u32,
            u: self.u_in,
            l: l as f32,
            pad: [0.0; 2],
        };
        let pbuf = self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("slice_params"),
            contents: bytemuck::bytes_of(&p),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        let size = (4 * w * h) as u64;
        let out = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("slice"),
            size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("slice_read"),
            size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let bg = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("slice"),
            layout: &self.slice_pipe.get_bind_group_layout(0),
            entries: &[
                buf_entry(0, &pbuf, None),
                buf_entry(1, &self.macro_buf, None),
                buf_entry(2, &self.wall, None),
                buf_entry(3, &out, None),
            ],
        });
        let mut enc = self.device.create_command_encoder(&Default::default());
        {
            let mut pass = enc.begin_compute_pass(&Default::default());
            pass.set_pipeline(&self.slice_pipe);
            pass.set_bind_group(0, &bg, &[]);
            pass.dispatch_workgroups((w as u32).div_ceil(8), (h as u32).div_ceil(8), 1);
        }
        self.queue.submit([enc.finish()]);
        let data = self.read(&out, &staging, size);
        Slice { w, h, data }
    }

    /// Image de la veine (RGBA, ligne du haut en premier) depuis le dernier
    /// champ écrit par `run`.
    pub fn render(&mut self, o: &RenderOptions) -> Vec<u8> {
        let (w, h) = (o.width.max(1), o.height.max(1));
        let n = self.n as u64;
        if self.deriv.is_none() {
            self.deriv = Some(self.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("deriv3d"),
                size: 16 * n,
                usage: wgpu::BufferUsages::STORAGE,
                mapped_at_creation: false,
            }));
        }
        let dims = [self.nx as u32, self.ny as u32, self.nz as u32, self.n as u32];
        let dbuf = self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("derive_params"),
            contents: bytemuck::cast_slice(&dims),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        // Mailles affichées : celles de la grille la plus fine qui entoure toute
        // la maquette (là où le calcul la voit), sinon celles de la veine.
        let vk = match self.refine.as_ref() {
            Some(r) if r.full > 0 && r.levels[r.full - 1].sparse.is_none() => Some(r.full - 1),
            _ => None,
        };
        let (vorg, vs) = match vk {
            Some(k) => {
                let r = self.refine.as_ref().unwrap();
                (r.origin[k], (1u64 << (k + 1)) as f64)
            }
            None => ([0.0; 3], 1.0),
        };
        let vdims = {
            let g = match vk {
                Some(k) => &self.refine.as_ref().unwrap().levels[k],
                None => &*self,
            };
            [g.nx as f32, g.ny as f32, g.nz as f32, 0.0]
        };
        let c = &o.camera;
        let v4 = |a: [f64; 3], w: f64| [a[0] as f32, a[1] as f32, a[2] as f32, w as f32];
        let ru = RenderUniform {
            eye: v4(c.eye, 1.0),
            fwd: v4(c.fwd, 0.0),
            right: v4(c.right, w as f64 / h as f64),
            up: v4(c.up, c.tan_half),
            dims: [self.nx as f32, self.ny as f32, self.nz as f32, o.l as f32],
            bg: [o.background[0], o.background[1], o.background[2], 1.0],
            body: [o.body[0], o.body[1], o.body[2], 1.0],
            plane: [o.plane_axis.min(2) as f32, o.plane_pos as f32, o.vorticity as u8 as f32, self.u_in],
            iso: [
                (o.iso_level * (self.u_in as f64 / o.l).powi(2)) as f32,
                o.iso as u8 as f32,
                0.8,
                4.0,
            ],
            size: [w as f32, h as f32, o.body_view as u8 as f32, 0.0],
            sdf0: self.sdf_geom.0,
            sdfn: self.sdf_geom.1,
            vorg: [vorg[0] as f32, vorg[1] as f32, vorg[2] as f32, vs as f32],
            vdims,
        };
        let rbuf = self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("render_params"),
            contents: bytemuck::bytes_of(&ru),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        let size = (4 * w * h) as u64;
        if !matches!(&self.render_bufs, Some((bw, bh, _, _)) if *bw == w && *bh == h) {
            let img = self.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("image3d"),
                size,
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
                mapped_at_creation: false,
            });
            let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("image3d_read"),
                size,
                usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });
            self.render_bufs = Some((w, h, img, staging));
        }
        let (_, _, img, staging) = self.render_bufs.as_ref().unwrap();
        let deriv = self.deriv.as_ref().unwrap();
        let vox = match vk {
            Some(k) => &self.refine.as_ref().unwrap().levels[k].vox,
            None => &self.vox,
        };
        let dbg = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("derive3d"),
            layout: &self.derive_pipe.get_bind_group_layout(0),
            entries: &[
                buf_entry(0, &dbuf, None),
                buf_entry(1, &self.macro_buf, None),
                buf_entry(2, &self.wall, None),
                buf_entry(3, deriv, None),
            ],
        });
        let rbg = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("render3d"),
            layout: &self.render_pipe.get_bind_group_layout(0),
            entries: &[
                buf_entry(0, &rbuf, None),
                buf_entry(1, &self.macro_buf, None),
                buf_entry(3, deriv, None),
                buf_entry(4, img, None),
                buf_entry(5, &self.sdf, None),
                buf_entry(6, vox, None),
            ],
        });
        let mut enc = self.device.create_command_encoder(&Default::default());
        {
            let mut pass = enc.begin_compute_pass(&Default::default());
            pass.set_pipeline(&self.derive_pipe);
            pass.set_bind_group(0, &dbg, &[]);
            pass.dispatch_workgroups(self.disp.0, self.disp.1, 1);
            pass.set_pipeline(&self.render_pipe);
            pass.set_bind_group(0, &rbg, &[]);
            pass.dispatch_workgroups((w as u32).div_ceil(8), (h as u32).div_ceil(8), 1);
        }
        self.queue.submit([enc.finish()]);
        let v = self.read(img, staging, size);
        if v.is_empty() {
            self.lost = true;
            return vec![0; 4 * w * h];
        }
        bytemuck::cast_slice::<f32, u8>(&v).to_vec()
    }

    /// Le calcul est-il resté borné ?
    pub fn is_finite(&self, f: &Force3) -> bool {
        f.fx.is_finite() && f.fy.is_finite() && f.fz.is_finite() && f.fx.abs() < 1e9
    }

    /// Le WGSL engendré (diagnostic).
    pub fn shader_source(fp16: bool) -> String {
        specialize(&wgsl_step(fp16, false, false, false).replace("WGu", "64u"), [1, 0, 0, 0, 0, 0, 0])
    }
}

#[cfg(test)]
mod tests {
    use super::oct_encode;

    /// Loi de Townsend–Mellor : nulle sans gradient, P y/(2κu_τ) aux faibles
    /// gradients, (2/κ)√(P y) à u_τ = 0 ; u_τ abaissé en gradient adverse,
    /// relevé en favorable, et racine de u = u_τ,S(u_t − Δ(u)).
    #[test]
    fn loi_mellor() {
        use super::{mellor_du, mellor_utau, spalding_utau};
        let k = 0.41;
        assert_eq!(mellor_du(0.002, 0.0, 2.0), 0.0);
        let (u, pg, y) = (0.002, 1e-9, 2.0);
        let lin = pg * y / (2.0 * k * u);
        assert!((mellor_du(u, pg, y) - lin).abs() < 1e-3 * lin);
        assert!((mellor_du(0.0, 4e-6, 2.0) - 2.0 / k * 8e-6f64.sqrt()).abs() < 1e-12);
        let (ut, nu) = (0.05, 1e-6);
        let u0 = spalding_utau(ut, y, nu);
        assert_eq!(mellor_utau(ut, y, nu, 0.0), u0);
        for pg in [2e-6, 2e-5, -2e-6] {
            let m = mellor_utau(ut, y, nu, pg);
            assert!(if pg > 0.0 { m < u0 } else { m > u0 }, "{pg} {m} {u0}");
            let r = m - spalding_utau(ut - mellor_du(m, pg, y), y, nu);
            assert!(r.abs() < 1e-9, "{pg} {r}");
        }
    }

    /// Décodage identique à `oct_decode` (WGSL).
    fn oct_decode(k: f32) -> [f64; 3] {
        let k = k as f64;
        let iu = (k / 2048.0).floor();
        let iv = k - iu * 2048.0;
        let (mut u, mut v) = (iu / 2047.0 * 2.0 - 1.0, iv / 2047.0 * 2.0 - 1.0);
        let z = 1.0 - u.abs() - v.abs();
        if z < 0.0 {
            let a = u;
            u = (1.0 - v.abs()) * a.signum();
            v = (1.0 - a.abs()) * v.signum();
        }
        let l = (u * u + v * v + z * z).sqrt();
        [u / l, v / l, z / l]
    }

    #[test]
    fn normale_octaedrique() {
        // Normales réparties sur la sphère : erreur angulaire < 0,2°.
        let mut worst: f64 = 1.0;
        for i in 0..400 {
            let t = (i as f64 + 0.5) / 400.0;
            let z = 1.0 - 2.0 * t;
            let r = (1.0 - z * z).sqrt();
            let ph = i as f64 * 2.399_963;
            let n = [r * ph.cos(), r * ph.sin(), z];
            let k = oct_encode(n);
            assert_eq!(k, k.round(), "entier exact en f32");
            let m = oct_decode(k);
            worst = worst.min(n[0] * m[0] + n[1] * m[1] + n[2] * m[2]);
        }
        assert!(worst > (0.2f64).to_radians().cos(), "cos = {worst}");
    }
}

#[cfg(test)]
mod tests_envergure {
    use super::span_distribution;

    /// Bande rectangulaire découpée en stations de 0,087 (deux triangles
    /// par station), chargement uniforme par unité d'aire : la répartition
    /// par tranches de 1/32 doit être plate (au centre de centre posé, elle
    /// oscillait du simple au double).
    #[test]
    fn repartition_plate() {
        let (h, n) = (0.087, 46);
        let mut spans = Vec::new();
        let mut vals = Vec::new();
        for i in 0..n {
            let (a, b) = (i as f64 * h, (i + 1) as f64 * h);
            // Deux triangles d'aire h/2 (corde 1) : sommets en y (a, a, b) et (a, b, b).
            spans.push([a, a, b]);
            spans.push([a, b, b]);
            vals.push([0.0, 0.0, 0.5 * h]);
            vals.push([0.0, 0.0, 0.5 * h]);
        }
        let dy = 1.0 / 32.0;
        let nb = ((n as f64 * h) / dy).floor() as usize;
        let d = span_distribution(&spans, &vals, 0.0, dy, nb + 1);
        for (k, v) in d.iter().take(nb).enumerate() {
            assert!((v[2] / dy - 1.0).abs() < 1e-9, "tranche {k} : {}", v[2] / dy);
        }
    }
}
