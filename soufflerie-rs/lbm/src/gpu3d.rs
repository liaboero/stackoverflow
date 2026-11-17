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

pub(crate) const SOLID: u32 = 0xffff_ffff;

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
  xtr: f32, band: f32, pad0: f32, pad1: f32,
};
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

fn cell(xi: i32, yi: i32, zi: i32) -> u32 {
  let ny = i32(P.ny);
  let nz = i32(P.nz);
  let yy = (yi + ny) % ny;
  var zz = (zi + nz) % nz;
  if ((P.rows & 2u) != 0u) { zz = clamp(zi, 0, nz - 1); }
  let xx = clamp(xi, 0, i32(P.nx) - 1);
  return u32(xx) + P.nx * (u32(yy) + P.ny * u32(zz));
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let idx = gid.x + gid.y * nwg.x * 64u;
  if (idx >= P.N) { return; }
  let x = i32(idx % P.nx);
  let y = i32((idx / P.nx) % P.ny);
  let z = i32(idx / (P.nx * P.ny));
  let wk0 = wall[idx];
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
  } else if (wk0 != 0u) {
    // Premier nœud : ν̃ = κ u_τ h, viscosité turbulente de la couche
    // pariétale [0, h] de la loi de paroi (hauteur de raccord fixe : la
    // valeur ne saute pas quand la paroi traverse une rangée de nœuds).
    let k = wk0 - 1u;
    out = min(NU_TILDE_MAX, KAPPA * utau[4u * k] * WM_H);
  } else {
    let ixm = cell(x - 1, y, z);
    let ixp = cell(x + 1, y, z);
    let iym = cell(x, y - 1, z);
    let iyp = cell(x, y + 1, z);
    let izm = cell(x, y, z - 1);
    let izp = cell(x, y, z + 1);
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
  }
  nutOut[idx] = out;
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

fn wgsl_step(fp16: bool) -> String {
    let words = if fp16 { 10 } else { 19 };
    let half = if fp16 { 5 } else { 10 };
    let load = if fp16 {
        format!(
            "fn pop(i: u32, c: u32) -> f32 {{
  let w = i >> 1u;
  var v: vec2<f32>;
  if (w < {half}u) {{ v = unpack2x16float(srcA[w * P.N + c]); }} else {{ v = unpack2x16float(srcB[(w - {half}u) * P.N + c]); }}
  return select(v.x, v.y, (i & 1u) == 1u);
}}"
        )
    } else {
        format!(
            "fn pop(i: u32, c: u32) -> f32 {{
  if (i < {half}u) {{ return bitcast<f32>(srcA[i * P.N + c]); }}
  return bitcast<f32>(srcB[(i - {half}u) * P.N + c]);
}}"
        )
    };
    let store = if fp16 {
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
        let _ = write!(
            pull,
            "  {{
    var v: f32;
    var q = 0.0;
    if (hasWall) {{ q = qbuf[wk * 18u + {im1}u]; }}
    if (q > 0.0) {{
      let fj = pop({j}u, idx);
      let qo = qbuf[wk * 18u + {jm1}u];
      // Paroi mobile (glissement de la loi de paroi) : 6 w_i ρ (c_i · u_w),
      // u_w ramenée dans le plan de la face que touche le lien (aux arêtes,
      // la normale de la maille n'est celle d'aucune face : sa composante
      // normale ferait souffler la paroi).
      var uwl = uw;
      if (P.sa == 1u) {{
        let nf = oct_decode(qbuf[arrayLength(&qbuf) / 2u + wk * 18u + {im1}u]);
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
    }}{ground} else if (x - {cx} < 0) {{
      v = feqIn({i}u);
    }} else if (x - {cx} >= i32(P.nx)) {{
      v = pop({i}u, idx);
    }} else {{
      v = pop({i}u, cell(x - {cx}, y - {cy}, z - {cz}));
    }}
    g[{i}] = v;
  }}
",
            im1 = i - 1,
            jm1 = j - 1,
            ground = ground_code(i),
            w2 = f(2.0 * w3(i)),
            w6 = f(6.0 * w3(i)),
            ccx = f(cx as f64),
            ccy = f(cy as f64),
            ccz = f(cz as f64),
            mcx = f(-cx as f64),
            mcy = f(-cy as f64),
            mcz = f(-cz as f64),
        );
    }
    let dotc = |i: usize, a: usize, b: usize| C3[i][a] * C3[i][b];
    let mut coll = String::new();
    for i in 0..19 {
        let [cx, cy, cz] = C3[i];
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
        let (cxf, cyf, czf) = (cx as f64, cy as f64, cz as f64);
        let qp = format!(
            "{} * pxx + {} * pyy + {} * pzz + {} * pxy + {} * pxz + {} * pyz",
            f(cxf * cxf - 1.0 / 3.0),
            f(cyf * cyf - 1.0 / 3.0),
            f(czf * czf - 1.0 / 3.0),
            f(2.0 * cxf * cyf),
            f(2.0 * cxf * czf),
            f(2.0 * cyf * czf)
        );
        let w = w3(i);
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
    let feq = {wf} * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - usq);
    let fnew = feq + k * {w45} * ({qp}) + {w135} * ({h3}) + {w3f} * rho * ({cxs} * P.accel);
    o[{i}] = fnew - {wf};
    if (sig > 0.0) {{
      let ci = {cxs} * P.uIn;
      o[{i}] -= sig * (o[{i}] - {wf} * (3.0 * ci + 4.5 * ci * ci - 1.5 * P.uIn * P.uIn));
    }}
  }}
",
            wf = f(w),
            w45 = f(4.5 * w),
            w135 = f(13.5 * w),
            w3f = f(3.0 * w),
            cxs = f(cxf),
        );
    }
    let out: String = if fp16 {
        (0..words)
            .map(|w| {
                let b = if 2 * w + 1 < 19 {
                    format!("o[{}]", 2 * w + 1)
                } else {
                    "0.0".into()
                };
                format!("  putWord({w}u, idx, o[{}], {b});\n", 2 * w)
            })
            .collect()
    } else {
        (0..19).map(|i| format!("  putPop({i}u, idx, o[{i}]);\n")).collect()
    };
    let feq_in: String = (1..19)
        .map(|i| format!("if (i == {i}u) {{ cx = {}; w = {}; }}\n  ", f(C3[i][0] as f64), f(w3(i))))
        .collect();
    let fv: String = (0..19).map(|i| format!("fv[{i}] = g[{i}] + {}; ", f(w3(i)))).collect();
    format!(
        "{PARAMS_WGSL}
@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read> srcA: array<u32>;
@group(0) @binding(2) var<storage, read> srcB: array<u32>;
@group(0) @binding(3) var<storage, read_write> dstA: array<u32>;
@group(0) @binding(4) var<storage, read_write> dstB: array<u32>;
@group(0) @binding(5) var<storage, read> wall: array<u32>;
@group(0) @binding(6) var<storage, read> qbuf: array<f32>;
@group(0) @binding(7) var<storage, read_write> force: array<vec4<f32>>;
@group(0) @binding(8) var<storage, read_write> macroBuf: array<vec4<f32>>;
// Spalart–Allmaras : ν̃ par maille ; par maille pariétale, normale et distance
// du nœud à la paroi (x, y, z, d), et u_τ écrit pour la passe de transport.
@group(0) @binding(9) var<storage, read> nut: array<f32>;
@group(0) @binding(10) var<storage, read> wallGeo: array<vec4<f32>>;
@group(0) @binding(11) var<storage, read_write> utau: array<f32>;
// Distance à la paroi de chaque maille (amortissement de la sous-maille).
@group(0) @binding(12) var<storage, read> wdist: array<f32>;
// Étendue en x du corps par tranche y (zone laminaire).
@group(0) @binding(13) var<storage, read> slab: array<vec2<f32>>;
{lam_fns}
{sa_fns}
{load}
{store}

fn cell(xi: i32, yi: i32, zi: i32) -> u32 {{
  let ny = i32(P.ny);
  let nz = i32(P.nz);
  let yy = (yi + ny) % ny;
  let zz = (zi + nz) % nz;
  let xx = clamp(xi, 0, i32(P.nx) - 1);
  return u32(xx) + P.nx * (u32(yy) + P.ny * u32(zz));
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
  var s = max(max(1.0 - dy / wy, sz), 1.0 - f32(x) / wx);
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

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {{
  let idx = gid.x + gid.y * nwg.x * 64u;
  if (idx >= P.N) {{ return; }}
  let wk0 = wall[idx];
  // Solide, ou maille fantôme d'une grille fine (remplie par refine3d).
  if (wk0 >= 0xfffffffeu) {{ return; }}
  let wk = select(wk0 - 1u, 0xffffffffu, wk0 == 0u);
  let hasWall = wk != 0xffffffffu;
  let x = i32(idx % P.nx);
  let y = i32((idx / P.nx) % P.ny);
  let z = i32(idx / (P.nx * P.ny));
  var g: array<f32, 19>;
  var fsum = vec3<f32>(0.0);
  // Loi de paroi : glissement tangentiel tiré de la vitesse du pas précédent.
  var uw = vec3<f32>(0.0);
  // Viscosité turbulente de la couche pariétale (κ u_τ h), < 0 hors paroi.
  var nuw = -1.0;
  if (P.sa == 1u && hasWall) {{
    let gw0 = wallGeo[wk];
    // Arête vive (distance négative) : pas de loi de paroi, paroi adhérente.
    let edge = gw0.w < 0.0;
    let gw = vec4<f32>(gw0.xyz, abs(gw0.w));
    let u0 = macroBuf[idx].xyz;
    // Vitesse tangentielle du nœud et sa distance à la paroi.
    let tv0 = u0 - dot(u0, gw.xyz) * gw.xyz;
    let ut0 = length(tv0);
    let y0 = max(0.1, gw.w);
    var ut = ut0;
    var yw = y0;
    {{
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
      var ok = true;
      for (var c = 0u; c < 8u; c++) {{
        let o = vec3<i32>(i32(c & 1u), i32((c >> 1u) & 1u), i32((c >> 2u) & 1u));
        let q3 = i0 + o;
        if (q3.x < 0 || q3.x >= i32(P.nx) || q3.z < 0 || q3.z >= i32(P.nz)) {{ ok = false; break; }}
        let ci = cell(q3.x, q3.y, q3.z);
        if (wall[ci] >= 0xfffffffeu) {{ ok = false; break; }}
        let wv = select(1.0 - fr, fr, vec3<bool>(o.x == 1, o.y == 1, o.z == 1));
        us3 += wv.x * wv.y * wv.z * macroBuf[ci].xyz;
        ds += wv.x * wv.y * wv.z * wdist[ci];
      }}
      if (ok && abs(ds - WM_H) < 0.5) {{
        ut = length(us3 - dot(us3, gw.xyz) * gw.xyz);
        yw = WM_H;
      }}
    }}
    var ut2 = spalding_utau(ut, yw);
    if (P.xtr > 0.0) {{
      // Zone laminaire : Cf de Blasius, 0,664/√Re_x, plafonné au point d'arrêt.
      let xl = lam_xl(x, y);
      if (xl.y > 0.0 && xl.x < P.xtr * xl.y) {{
        let rex = max(1.0, P.uIn * max(xl.x, 0.5) / P.nu);
        let cf = min(0.02, 0.664 / sqrt(rex));
        ut2 = ut * sqrt(0.5 * cf);
      }}
    }}
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
    var usv = vec3<f32>(0.0);
    if (edge) {{
      // Arête vive : l'écoulement y décolle et la loi de paroi (couche
      // limite attachée le long d'une face) ne s'applique pas ; sans elle,
      // la paroi presque glissante laisse croître la vitesse de contournement
      // de l'arête (singulière en non visqueux) jusqu'à la divergence.
      nuw = -1.0;
    }} else if (ut0 >= 1e-12) {{
      // Profil monotone : la paroi ne glisse pas plus vite que ne le permet
      // aussi la vitesse relevée à h (u_h − pente · h). Les deux estimations
      // coïncident en régime établi ; la borne arrête un nœud qui
      // s'emballerait (un nœud presque à la paroi glisse presque librement
      // et rien d'autre ne le retient).
      let us = max(0.0, min(ut0 - ut2 * ut2 * y0 / (P.nu + nuw), ut - ut2 * ut2 * yw / (P.nu + nuw)));
      usv = tv0 * (us / ut0);
    }}
    // Relaxation sur une dizaine de pas : une paroi qui suit à chaque pas la
    // vitesse de son propre nœud (glissement presque libre près de la paroi)
    // n'a pas de rappel et diverge ; en régime établi, sans effet (β = 0,1 ou
    // 0,02 : mêmes efforts à 0,3 % près).
    let prev = vec3<f32>(utau[4u * wk + 1u], utau[4u * wk + 2u], utau[4u * wk + 3u]);
    uw = prev + 0.1 * (usv - prev);
    utau[4u * wk + 1u] = uw.x;
    utau[4u * wk + 2u] = uw.y;
    utau[4u * wk + 3u] = uw.z;
  }}
  // Décambrage visqueux : vitesse de transpiration de la maille pariétale,
  // rangée après u_τ et le glissement (4 + 1 valeurs par maille).
  var vn = 0.0;
  if (P.sa == 1u && hasWall) {{
    vn = utau[arrayLength(&utau) / 5u * 4u + wk];
  }}
  g[0] = pop(0u, idx);
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
  if (P.sa == 1u) {{
    // Nœud pariétal : la viscosité de la couche [0, h] supposée par le
    // glissement, pour que la contrainte transmise soit τ_p.
    tau += 3.0 * select(nu_turb(nut[idx]), nuw, nuw >= 0.0);
  }}
  if (u32(x) > P.spongeStart) {{
    let s = f32(u32(x) - P.spongeStart) / f32(max(1u, P.nx - 1u - P.spongeStart));
    tau += (P.tauSponge - P.tau0) * s * s;
  }}
  var om = 1.0 / tau;
  if (P.les == 1u) {{
    let qn = sqrt(2.0 * (pxx * pxx + pyy * pyy + pzz * pzz + 2.0 * (pxy * pxy + pxz * pxz + pyz * pyz)));
    // Amortissement de Van Driest près de la paroi, (1 − e^(−d/A))², comme en
    // 2D, sans loi de paroi seulement : avec elle, la sous-maille stabilise les
    // premiers nœuds à τ ≈ 0,5 (aile à Re = 10⁶, L = 32 : divergence sinon).
    var damp = 1.0;
    if (P.sa == 0u) {{ damp = 1.0 - exp(-wdist[idx] / P.lesA); }}
    om = 2.0 / (tau + sqrt(tau * tau + P.c18 * damp * damp * qn * ir));
  }}
  let k = 1.0 - om;
  // Régularisation récursive (Malaspinas 2015) : moments d'ordre 3 tirés de
  // ρ, u et de la partie hors équilibre d'ordre 2 après collision.
  let bxx = k * pxx;
  let byy = k * pyy;
  let bzz = k * pzz;
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
{out}
  if (hasWall) {{
    force[wk] += vec4<f32>(fsum, 0.0);
  }}
  if (P.writeMacro == 1u || P.sa == 1u) {{
    macroBuf[idx] = vec4<f32>(u, rho);
  }}
}}
",
        sa_fns = SA_FNS,
        lam_fns = LAM_FNS,
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

fn wgsl_init(fp16: bool) -> String {
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
    let put: String = if fp16 {
        (0..words)
            .map(|w| {
                let b = if 2 * w + 1 < 19 {
                    format!("e{}", 2 * w + 1)
                } else {
                    "0.0".into()
                };
                let v = format!("pack2x16float(vec2<f32>(e{}, {b}))", 2 * w);
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
    format!(
        "{PARAMS_WGSL}
@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read_write> A: array<u32>;
@group(0) @binding(2) var<storage, read_write> B: array<u32>;
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
fn voxKind(c: vec3<i32>) -> u32 {
  let n = vec3<i32>(V.dims.xyz);
  if (any(c < vec3<i32>(0)) || any(c >= n)) { return 0u; }
  let i = u32(c.x) + u32(n.x) * (u32(c.y) + u32(n.y) * u32(c.z));
  return (vox[i >> 4u] >> ((i & 15u) * 2u)) & 3u;
}

// Première maille non fluide le long du rayon (traversée d'Amanatides et Woo,
// mailles centrées sur les entiers) : (t, axe de la face d'entrée, type).
fn voxHit(o: vec3<f32>, dir: vec3<f32>, t0: f32, t1: f32) -> vec3<f32> {
  let p = o + dir * t0;
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
  for (var k = 0; k < 6000; k++) {
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
      let fr = abs(fract(p + vec3<f32>(0.5)) - vec3<f32>(0.5));
      var e = 1.0;
      for (var b = 0; b < 3; b++) {
        if (b != a) { e = min(e, 0.5 - fr[b]); }
      }
      if (e < 0.05) { col = mix(col, V.bg.rgb, 0.55); }
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
      let inBody = select(S(p) < V.sdfn.w, voxKind(vec3<i32>(round(p))) != 0u, mode == 2u);
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
    pad: [f32; 2],
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

pub struct Gpu3D {
    pub(crate) device: wgpu::Device,
    pub(crate) queue: wgpu::Queue,
    adapter_name: String,
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub n: usize,
    pub fp16: bool,
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
    /// ν̃ initial (3ν dans le fluide, 0 dans le corps), pour `reset`.
    nut0: Vec<f32>,
    force: wgpu::Buffer,
    force_read: wgpu::Buffer,
    pub(crate) params: wgpu::Buffer,
    pub(crate) macro_buf: wgpu::Buffer,
    pub(crate) step_pipe: wgpu::ComputePipeline,
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
    pub(crate) cur: usize,
    pub time: u64,
    pub(crate) force_steps: u64,
    pub wall_count: usize,
    pub solid_count: usize,
    /// Mailles solides (pour l'affichage).
    solid: Vec<bool>,
    /// Le GPU ne répond plus (lecture impossible, pilote réinitialisé).
    pub lost: bool,
    /// Grille fine imbriquée : enveloppe de mailles fantômes, ni entrée, ni
    /// sortie, ni bandes absorbantes.
    pub(crate) nested: bool,
    /// Sol mobile en z = −½ et plafond glissant (au lieu de bords périodiques en z).
    ground: bool,
    /// Grilles fines autour de la maquette (raffinement local).
    pub(crate) refine: Option<Box<crate::refine3d::Refinement>>,
    /// Maille de chaque maille pariétale (dans l'ordre de leur numéro).
    wall_cells: Vec<u32>,
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
    pub fn with_turbulence(body: &Body3D, nu: f64, u_in: f32, fp16: bool, sa: bool) -> Result<Gpu3D, String> {
        Self::with_refinement(body, nu, u_in, fp16, sa, 0)
    }

    /// Idem, avec `depth` grilles imbriquées autour de la maquette (0 à 2),
    /// chacune deux fois plus fine que la précédente.
    pub fn with_refinement(
        body: &Body3D,
        nu: f64,
        u_in: f32,
        fp16: bool,
        sa: bool,
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
        sa: bool,
        depth: usize,
        le: bool,
    ) -> Result<Gpu3D, String> {
        let (device, queue, name) = Self::open_device()?;
        let mut top = Self::build(device.clone(), queue.clone(), name.clone(), body, nu, u_in, fp16, sa, false)?;
        let depth = depth.min(crate::refine3d::EXTENT3.len());
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
            let (up, down, side) = crate::refine3d::EXTENT3[k - 1];
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
            let mut child = Self::build(device.clone(), queue.clone(), name.clone(), &cb, nu_k, u_in, fp16, sa, true)?;
            child.les = top.les;
            {
                let parent = if k == 1 { &top } else { &levels[k - 2] };
                links.push(crate::refine3d::Link3::new(parent, &child, lo));
            }
            levels.push(child);
            parent_body = cb;
            porg = [0, 1, 2].map(|a| porg[a] + lo[a] as f64 / sp);
            origin.push(porg);
        }
        let full = levels.len();
        if le {
            let k = depth + 1;
            let sp = (1u64 << (k - 1)) as f64;
            let mesh = parent_body.mesh.clone().ok_or("raffinement du bord d'attaque : maquette maillée seulement")?;
            let nrm = |t: &crate::stl::Tri| {
                let e1 = [0, 1, 2].map(|a| t[1][a] - t[0][a]);
                let e2 = [0, 1, 2].map(|a| t[2][a] - t[0][a]);
                [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]]
            };
            let vol: f64 = mesh.tris.iter().map(|t| { let n = nrm(t); (0..3).map(|a| t[0][a] * n[a]).sum::<f64>() }).sum();
            let sgn = if vol >= 0.0 { 1.0 } else { -1.0 };
            let (mut a0, mut a1) = ([f64::MAX; 3], [f64::MIN; 3]);
            for t in &mesh.tris {
                let n = nrm(t);
                let m = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
                if m < 1e-12 || sgn * n[0] / m > -0.5 {
                    continue;
                }
                for v in t {
                    for a in 0..3 {
                        a0[a] = a0[a].min(v[a]);
                        a1[a] = a1[a].max(v[a]);
                    }
                }
            }
            if a0[0] > a1[0] {
                return Err("raffinement du bord d'attaque : aucune face tournée vers l'amont".into());
            }
            // Marges (mailles de la grille parente) : de quoi loger le pic de
            // dépression et le début de la recompression (10 % de la longueur
            // de référence vers l'aval), plus la zone de recouvrement.
            let e = (crate::refine3d::MARGIN3 + 4) as f64;
            let down = e.max(0.1 * parent_body.l as f64);
            let (pnx, pny, pnz) = (parent_body.nx, parent_body.ny, parent_body.nz);
            let lo = [
                (a0[0] - e).floor().max(3.0) as usize,
                (a0[1] - e).floor().max(3.0) as usize,
                (a0[2] - e).floor().max(3.0) as usize,
            ];
            let hi = [
                ((a1[0] + down).ceil() as usize).min(pnx - 4),
                ((a1[1] + e).ceil() as usize).min(pny - 4),
                ((a1[2] + e).ceil() as usize).min(pnz - 4),
            ];
            let dims = [0, 1, 2].map(|a| 2 * (hi[a] - lo[a]) + 1);
            // Une seule boîte alignée sur les axes : sur une aile en flèche,
            // elle couvre presque toute l'aile. Au-delà de quatre fois la
            // grille parente, elle coûte autant qu'un raffinement complet
            // (environ cinq fois) et risque d'épuiser la mémoire du GPU.
            let n_le = dims.iter().product::<usize>();
            if n_le > 4 * pnx * pny * pnz {
                return Err(format!(
                    "raffinement du bord d'attaque : {:.0} M mailles, plus de quatre fois la grille parente (bord d'attaque en flèche ?) : raffiner toute la maquette (--raff) plutôt",
                    n_le as f64 / 1e6
                ));
            }
            let cb = parent_body.refined([lo[0] as f64, lo[1] as f64, lo[2] as f64], dims);
            let nu_k = nu * (1u64 << k) as f64;
            let mut child = Self::build(device.clone(), queue.clone(), name.clone(), &cb, nu_k, u_in, fp16, sa, true)?;
            child.les = top.les;
            {
                let parent = if k == 1 { &top } else { &levels[k - 2] };
                links.push(crate::refine3d::Link3::new(parent, &child, lo));
            }
            levels.push(child);
            porg = [0, 1, 2].map(|a| porg[a] + lo[a] as f64 / sp);
            origin.push(porg);
        }
        top.refine = Some(Box::new(crate::refine3d::Refinement {
            levels,
            links,
            fresh: true,
            origin,
            full,
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
        limits.max_storage_buffers_per_shader_stage = limits.max_storage_buffers_per_shader_stage.max(13);
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            label: Some("soufflerie 3D"),
            required_limits: limits,
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
        sa: bool,
        nested: bool,
    ) -> Result<Gpu3D, String> {
        let (nx, ny, nz) = (body.nx, body.ny, body.nz);
        let n = nx * ny * nz;
        let al = device.limits();
        let words_a = if fp16 { 5 } else { 10 };
        let words_b = if fp16 { 5 } else { 9 };
        let biggest = (words_a * n * 4).max(16 * n) as u64;
        let max_bind = (al.max_storage_buffer_binding_size as u64).min(al.max_buffer_size);
        if biggest > max_bind {
            return Err(format!(
                "{:.1} M mailles : tampon de {} Mo, au-delà de la limite de l'adaptateur ({} Mo)",
                n as f64 / 1e6,
                biggest >> 20,
                max_bind >> 20
            ));
        }
        let st = wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC;
        let mk = |label: &str, size: u64, usage: wgpu::BufferUsages| {
            device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label),
                size: size.max(16),
                usage,
                mapped_at_creation: false,
            })
        };
        let pops = [
            [mk("a0", (words_a * n * 4) as u64, st), mk("b0", (words_b * n * 4) as u64, st)],
            [mk("a1", (words_a * n * 4) as u64, st), mk("b1", (words_b * n * 4) as u64, st)],
        ];
        let params = mk(
            "params3d",
            SLOT * BATCH as u64,
            wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        );
        let macro_buf = mk("macro3d", 16 * n as u64, st);

        // Voxelisation : mailles solides, puis mailles pariétales et fractions q.
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
        // Normale de la face touchée par chaque lien (codage octaédrique).
        let mut qn: Vec<f32> = Vec::new();
        // Plus petit cosinus entre les normales des faces que touchent les
        // liens de chaque maille pariétale (arêtes vives : < 0,5).
        let mut edge_cos: Vec<f32> = Vec::new();
        let mut k = 0u32;
        // Mailles fluides qui ne portent que des liens de paroi mince (affichage).
        let mut thin: Vec<usize> = Vec::new();
        for z in z0.saturating_sub(1)..=(z1 + 1).min(nz - 1) {
            for y in y0.saturating_sub(1)..=(y1 + 1).min(ny - 1) {
                for x in (x0 - 1).max(1)..=(x1 + 1).min(nx - 2) {
                    let c = idx(x, y, z);
                    if wall[c] == SOLID {
                        continue;
                    }
                    let mut row = [0f32; 18];
                    let mut nrow = [0f32; 18];
                    let mut nraw: Vec<[f64; 3]> = Vec::new();
                    let mut any = false;
                    let mut beside_solid = false;
                    for i in 1..19 {
                        let [cx, cy, cz] = C3[i];
                        let sx = (x as i32 - cx) as usize;
                        let sy = ((y as i32 - cy + ny as i32) % ny as i32) as usize;
                        let sz = ((z as i32 - cz + nz as i32) % nz as i32) as usize;
                        let p = [x as f64, y as f64, z as f64];
                        let d = [-cx as f64, -cy as f64, -cz as f64];
                        // Maillage : intersection exacte du lien avec les
                        // triangles, y compris entre deux mailles fluides
                        // (paroi plus mince qu'une maille : aile, dérive…).
                        let hit = body.mesh.as_ref().and_then(|m| m.segment_hit_normal(p, d));
                        let solid = wall[idx(sx, sy, sz)] == SOLID;
                        beside_solid |= solid;
                        let q = match (hit, solid) {
                            (Some((t, n)), _) => Some(((t as f32).clamp(0.02, 1.0), Some(n))),
                            (None, true) => Some((link_q(body, p, d), None)),
                            (None, false) => None,
                        };
                        if let Some((q, n)) = q {
                            row[i - 1] = q;
                            // Sans triangle : gradient de la distance juste devant la paroi.
                            let n = n.unwrap_or_else(|| {
                                let t = (q as f64 - 0.25).max(0.0);
                                body.wall_normal([p[0] + t * d[0], p[1] + t * d[1], p[2] + t * d[2]]).0
                            });
                            nrow[i - 1] = oct_encode(n);
                            nraw.push(n);
                            any = true;
                        }
                    }
                    if any {
                        k += 1;
                        wall[c] = k;
                        qs.extend_from_slice(&row);
                        qn.extend_from_slice(&nrow);
                        let mut cmin = 1.0f64;
                        for a in 0..nraw.len() {
                            for b in a + 1..nraw.len() {
                                cmin = cmin.min((0..3).map(|t| nraw[a][t] * nraw[b][t]).sum::<f64>());
                            }
                        }
                        edge_cos.push(cmin as f32);
                        if !beside_solid {
                            thin.push(c);
                        }
                    }
                }
            }
        }
        let wall_count = k as usize;
        // Normales des liens rangées après les fractions q (même tampon).
        qs.extend_from_slice(&qn);
        // Maille de chaque lien pariétal (numéro → maille), pour les efforts par maille.
        let rest_force: Vec<[f32; 3]> = qs
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
            .collect();
        let mut wall_cells = vec![0u32; wall_count];
        for (c, &w) in wall.iter().enumerate() {
            if w != 0 && w < SOLID - 1 {
                wall_cells[w as usize - 1] = c as u32;
            }
        }
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
        if qs.is_empty() {
            qs.resize(18, 0.0);
        }
        let solid: Vec<bool> = wall.iter().map(|&w| w == SOLID).collect();
        // Affichage des mailles du calcul : 2 bits par maille (1 solide, 2 paroi mince).
        let mut vox = vec![0u32; n.div_ceil(16)];
        for (c, &s) in solid.iter().enumerate() {
            if s {
                vox[c >> 4] |= 1 << ((c & 15) * 2);
            }
        }
        for &c in &thin {
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
        let wall_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("wall3d"),
            contents: bytemuck::cast_slice(&wall),
            usage: wgpu::BufferUsages::STORAGE,
        });
        let qbuf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("q3d"),
            contents: bytemuck::cast_slice(&qs),
            usage: wgpu::BufferUsages::STORAGE,
        });
        // Spalart–Allmaras : normale et distance des mailles pariétales (dans
        // l'ordre de leur numéro), distance à la paroi de toutes les mailles.
        let (geo,) = if sa {
            let mut cells = vec![0usize; wall_count];
            for (c, &w) in wall.iter().enumerate() {
                if w != 0 && w < crate::refine3d::GHOST {
                    cells[w as usize - 1] = c;
                }
            }
            // Distance négative : maille d'arête vive (faces touchées à plus
            // de 60° l'une de l'autre).
            let geo: Vec<[f32; 4]> = cells
                .par_iter()
                .enumerate()
                .map(|(w, &c)| {
                    let p = [(c % nx) as f64, ((c / nx) % ny) as f64, (c / (nx * ny)) as f64];
                    let (nrm, d) = body.wall_normal(p);
                    let d = (d as f32).max(0.1);
                    // Avec les arêtes vives du maillage : seulement à portée
                    // de lien (√2) de l'une d'elles (pas autour d'un bord
                    // d'attaque arrondi plus fin qu'une maille).
                    let sharp = edge_cos[w] < 0.5
                        && body.features.as_ref().is_none_or(|f| f.distance(p, 2.0) <= std::f64::consts::SQRT_2);
                    let d = if sharp { -d } else { d };
                    [nrm[0] as f32, nrm[1] as f32, nrm[2] as f32, d]
                })
                .collect();
            (geo,)
        } else {
            (Vec::new(),)
        };
        let solid_mask: Vec<bool> = wall.iter().map(|&w| w == SOLID).collect();
        let dist = body.wall_distance(&solid_mask);
        // Étendue en x du corps (mailles solides et pariétales) dans chaque tranche y.
        let mut slab = vec![[1.0f32, 0.0f32]; ny];
        for (c, &w) in wall.iter().enumerate() {
            if w != 0 && w != crate::refine3d::GHOST {
                let (x, y) = ((c % nx) as f32, (c / nx) % ny);
                if slab[y][0] > slab[y][1] {
                    slab[y] = [x, x];
                } else {
                    slab[y] = [slab[y][0].min(x), slab[y][1].max(x)];
                }
            }
        }
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
        let dist_buf = init_or("dist3d", bytemuck::cast_slice(&dist));
        // Par maille pariétale : u_τ et vitesse de glissement relaxée, puis
        // (à la suite, une valeur par maille) la vitesse de transpiration.
        let utau_buf = mk("utau3d", (20 * wall_count) as u64, st);
        let nsa = if sa { 4 * n as u64 } else { 16 };
        let nut = [mk("nut0", nsa, st), mk("nut1", nsa, st)];
        let nut0: Vec<f32> = if sa {
            let inf = crate::sa::nu_tilde_inf(nu);
            wall.iter().map(|&w| if w == SOLID { 0.0 } else { inf }).collect()
        } else {
            Vec::new()
        };
        let fsize = (16 * wall_count).max(16) as u64;
        let force = mk("force3d", fsize, st);
        let force_read = mk("force3d_read", fsize, wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST);

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
                storage(11, false),
                storage(12, true),
                storage(13, true),
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
            ],
        });
        let init_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("d3q19_init"),
            entries: &[uniform(0, true, psize), storage(1, false), storage(2, false)],
        });
        let pipe = |label: &str, layout: &wgpu::BindGroupLayout, code: String| {
            let pl = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some(label),
                bind_group_layouts: &[Some(layout)],
                immediate_size: 0,
            });
            let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some(label),
                source: wgpu::ShaderSource::Wgsl(code.into()),
            });
            device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some(label),
                layout: Some(&pl),
                module: &module,
                entry_point: Some("main"),
                compilation_options: Default::default(),
                cache: None,
            })
        };
        let step_pipe = pipe("d3q19", &layout, wgsl_step(fp16));
        let init_pipe = pipe("d3q19_init", &init_layout, wgsl_init(fp16));
        let sa_pipe = pipe("sa3d", &sa_layout, format!("{PARAMS_WGSL}{SA3_WGSL}"));
        let slice_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("slice"),
            entries: &[
                uniform(0, false, std::mem::size_of::<SliceParams>() as u64),
                storage(1, true),
                storage(2, true),
                storage(3, false),
            ],
        });
        let slice_pipe = pipe("slice", &slice_layout, SLICE_WGSL.to_string());
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
        let derive_pipe = auto("derive3d", DERIVE_WGSL);
        let render_pipe = auto("render3d", RENDER_WGSL);
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
                    ],
                })
            })
            .collect();
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
        let groups = (n as u32).div_ceil(64);
        let dx = groups.min(65535);
        let disp = (dx, groups.div_ceil(dx));
        let mut s = Gpu3D {
            device,
            queue,
            adapter_name,
            nx,
            ny,
            nz,
            n,
            fp16,
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
            nut0,
            force,
            force_read,
            params,
            macro_buf,
            step_pipe,
            init_pipe,
            slice_pipe,
            derive_pipe,
            render_pipe,
            deriv: None,
            render_bufs: None,
            bind,
            init_bind,
            disp,
            cur: 0,
            time: 0,
            force_steps: 0,
            wall_count,
            solid_count,
            solid,
            lost: false,
            nested,
            ground: body.ground && !nested,
            refine: None,
            wall_cells,
            edge_cos,
            rest_force,
            surface: body.mesh.clone(),
            transp: None,
        };
        s.set_viscosity(nu);
        s.reset(0.0);
        Ok(s)
    }

    pub fn adapter_name(&self) -> &str {
        &self.adapter_name
    }

    pub fn set_viscosity(&mut self, nu: f64) {
        self.set_viscosity_one(nu);
        if let Some(r) = self.refine.as_mut() {
            for (k, g) in r.levels.iter_mut().enumerate() {
                g.set_viscosity_one(nu * (2u64 << k) as f64);
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
            accel,
            c18: (18.0 * self.cs * self.cs) as f32,
            write_macro: write_macro as u32,
            les: self.les as u32,
            rows: self.nested as u32 | ((self.ground as u32) << 1),
            sa: self.sa as u32,
            nu: self.nu as f32,
            les_a: self.les_a as f32,
            xtr: self.xtr as f32,
            band: (self.l as f32 / 10.0).max(3.0),
            pad: [0.0; 2],
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
        if let Some(mut r) = self.refine.take() {
            crate::refine3d::run(self, &mut r, n, ramp);
            self.refine = Some(r);
            return;
        }
        // Soumissions courtes (≈ 3·10⁷ mises à jour de maille) : sous Windows, une
        // soumission de plus de 2 s fait réinitialiser le pilote, ce qui peut
        // arriver sur un GPU lent ou bridé en température.
        let per_submit = (30_000_000 / self.n.max(1)).clamp(1, BATCH);
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
                pass.set_pipeline(&self.step_pipe);
                for s in 0..k {
                    let off = (s as u64 * SLOT) as u32;
                    pass.set_bind_group(0, &self.bind[self.cur], &[off]);
                    pass.dispatch_workgroups(self.disp.0, self.disp.1, 1);
                    if self.sa {
                        pass.set_pipeline(&self.sa_pipe);
                        pass.set_bind_group(0, &self.sa_bind[self.cur], &[off]);
                        pass.dispatch_workgroups(self.disp.0, self.disp.1, 1);
                        pass.set_pipeline(&self.step_pipe);
                    }
                    self.cur = 1 - self.cur;
                }
            }
            self.queue.submit([enc.finish()]);
            self.time += k as u64;
            self.force_steps += k as u64;
            left -= k;
        }
    }

    fn read(&self, src: &wgpu::Buffer, staging: &wgpu::Buffer, size: u64) -> Vec<f32> {
        let mut enc = self.device.create_command_encoder(&Default::default());
        enc.copy_buffer_to_buffer(src, 0, staging, 0, size);
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
            // Grille du bord d'attaque : ses mailles pariétales hors de sa
            // bande de bord, et celles de la grille entière ailleurs (sous la
            // grille fine, la grille entière reçoit l'état recopié).
            let (fo, fs) = (r.origin[r.full - 1], (1u64 << r.full) as f64);
            let (lo, ls) = (r.origin[r.full], (1u64 << (r.full + 1)) as f64);
            let ld = {
                let g = &r.levels[r.full];
                [g.nx, g.ny, g.nz]
            };
            let m = (2 * crate::refine3d::MARGIN3 + 2) as f64;
            let inner = |q: [f64; 3]| (0..3).all(|a| q[a] >= m && q[a] <= (ld[a] - 1) as f64 - m);
            let a = r.levels[r.full - 1].take_forces_masked(|p| !inner([0, 1, 2].map(|i| (fo[i] + p[i] / fs - lo[i]) * ls)));
            let b = r.levels[r.full].take_forces_masked(|p| inner(p));
            for g in r.levels[r.full + 1..].iter_mut() {
                let _ = g.take_forces_one();
            }
            self.lost |= r.levels.iter().any(|g| g.lost);
            self.refine = Some(r);
            let sb = 4.0 * s;
            return Force3 {
                fx: a.fx / s + b.fx / sb,
                fy: a.fy / s + b.fy / sb,
                fz: a.fz / s + b.fz / sb,
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

    /// ν̃ de Spalart–Allmaras par maille (vide sans SA).
    pub fn nu_tilde(&self) -> Vec<f32> {
        if !self.sa {
            return Vec::new();
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
        }
        let mut levs = vec![Lev {
            mac: self.macros(),
            nut: self.nu_tilde(),
            d: [self.nx, self.ny, self.nz],
            org: [0.0; 3],
            s: 1.0,
            nu: self.nu,
            m: 0,
        }];
        if let Some(r) = self.refine.as_mut() {
            for (k, g) in r.levels.iter_mut().enumerate() {
                levs.push(Lev {
                    mac: g.macros(),
                    nut: g.nu_tilde(),
                    d: [g.nx, g.ny, g.nz],
                    org: r.origin[k],
                    s: (1u64 << (k + 1)) as f64,
                    nu: g.nu,
                    // Enveloppe fantôme et zone de recouvrement exclues.
                    m: 2 * crate::refine3d::MARGIN3 + 2,
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
                            let c = cx as usize + l.d[0] * ((cy + dj) as usize + l.d[1] * (cz + dk) as usize);
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
                (Some(i), Some(r)) => (r.origin[i], (1u64 << (i + 1)) as f64),
                _ => ([0.0; 3], 1.0),
            };
            let g = self.grid_mut(k);
            let (nx, ny) = (g.nx, g.ny);
            // Position dans la grille entière la plus fine.
            let tri: Vec<u32> = g
                .wall_cells
                .iter()
                .map(|&c| {
                    let c = c as usize;
                    let p = [(c % nx) as f64, ((c / nx) % ny) as f64, (c / (nx * ny)) as f64];
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

    /// Grille partielle du bord d'attaque : dimensions, mailles pariétales, facteur.
    pub fn le_grid(&self) -> Option<([usize; 3], usize, usize)> {
        let r = self.refine.as_ref()?;
        let g = r.levels.get(r.full)?;
        Some(([g.nx, g.ny, g.nz], g.wall_count, 1 << (r.full + 1)))
    }

    pub fn finest(&self) -> Option<([usize; 3], usize, usize)> {
        let r = self.refine.as_ref()?;
        let g = r.levels.get(r.full.checked_sub(1)?)?;
        Some(([g.nx, g.ny, g.nz], g.wall_count, 1 << r.full))
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
        }
        let mut levs: Vec<Lev> = Vec::new();
        match self.refine.as_mut() {
            Some(r) => {
                for k in r.full..=r.levels.len() {
                    let g = &mut r.levels[k - 1];
                    let mesh = g.surface.clone()?;
                    levs.push(Lev {
                        mac: g.macros(),
                        solid: g.solid.clone(),
                        d: [g.nx, g.ny, g.nz],
                        org: r.origin[k - 1],
                        s: (1u64 << k) as f64,
                        nu: g.nu,
                        mesh,
                        m: if k == r.full { 0.0 } else { (2 * crate::refine3d::MARGIN3 + 2) as f64 },
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
                let c = q[0] as usize + l.d[0] * (q[1] as usize + l.d[1] * q[2] as usize);
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
                let cell = q[0] as usize + l.d[0] * (q[1] as usize + l.d[1] * q[2] as usize);
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
                // Grille partielle : triangle entier loin de ses bords.
                if l.m > 0.0 && t.iter().any(|v| (0..3).any(|a| v[a] < l.m + 3.0 || v[a] > (l.d[a] - 1) as f64 - l.m - 3.0)) {
                    continue;
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
            let Some((li, n, c, ar, r)) = got else { continue };
            let l = &levs[li];
            // Surfaces ramenées aux mailles de la veine.
            let ar = ar / (l.s * l.s);
            area += ar;
            let Some((p, pw, vf, hf, v)) = r else {
                out.missed += ar;
                continue;
            };
            let un = vf[0] * n[0] + vf[1] * n[1] + vf[2] * n[2];
            let ut = [0, 1, 2].map(|k| vf[k] - un * n[k]);
            let utm = (ut[0] * ut[0] + ut[1] * ut[1] + ut[2] * ut[2]).sqrt();
            let utau = spalding_utau(utm, hf, l.nu);
            let tau = v[3] * utau * utau;
            let mut fl = [0.0; 3];
            for k in 0..3 {
                out.pressure[k] -= p * n[k] * ar;
                out.pressure_wall[k] -= pw * n[k] * ar;
                fl[k] = -pw * n[k] * ar;
                if utm > 0.0 {
                    out.friction[k] += tau * ut[k] / utm * ar;
                    fl[k] += tau * ut[k] / utm * ar;
                }
            }
            out.loads.push(([0, 1, 2].map(|k| l.org[k] + c[k] / l.s), fl));
            out.walls.push((n, pw));
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
        if g.wall_count == 0 {
            return out;
        }
        let v = g.read(&g.force, &g.force_read, 16 * g.wall_count as u64);
        let mut enc = g.device.create_command_encoder(&Default::default());
        enc.clear_buffer(&g.force, 0, None);
        g.queue.submit([enc.finish()]);
        let (nx, ny) = (g.nx, g.ny);
        for (k, c) in v.chunks_exact(4).enumerate() {
            let cell = g.wall_cells[k] as usize;
            let p = [(cell % nx) as f64, ((cell / nx) % ny) as f64, (cell / (nx * ny)) as f64];
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
        if steps == 0 || self.wall_count == 0 {
            return Force3 { steps, ..Default::default() };
        }
        let v = self.read(&self.force, &self.force_read, 16 * self.wall_count as u64);
        if v.is_empty() {
            self.lost = true;
            return Force3 { fx: f64::NAN, fy: f64::NAN, fz: f64::NAN, steps };
        }
        let mut enc = self.device.create_command_encoder(&Default::default());
        enc.clear_buffer(&self.force, 0, None);
        self.queue.submit([enc.finish()]);
        let (nx, ny) = (self.nx, self.ny);
        let mut f = [0.0f64; 3];
        for (k, c) in v.chunks_exact(4).enumerate() {
            let cell = self.wall_cells[k] as usize;
            if keep([(cell % nx) as f64, ((cell / nx) % ny) as f64, (cell / (nx * ny)) as f64]) {
                for a in 0..3 {
                    f[a] += c[a] as f64;
                }
            }
        }
        let s = steps as f64;
        Force3 { fx: f[0] / s, fy: f[1] / s, fz: f[2] / s, steps }
    }

    fn take_forces_one(&mut self) -> Force3 {
        let steps = std::mem::take(&mut self.force_steps);
        if steps == 0 || self.wall_count == 0 {
            let _ = self.device.poll(wgpu::PollType::wait_indefinitely());
            return Force3 {
                steps,
                ..Default::default()
            };
        }
        let v = self.read(&self.force, &self.force_read, 16 * self.wall_count as u64);
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
                buf_entry(6, &self.vox, None),
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
        wgsl_step(fp16)
    }
}

#[cfg(test)]
mod tests {
    use super::oct_encode;

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
