//! Raffinement local de la soufflerie 3D : grilles D3Q19 imbriquées autour de
//! la maquette, chacune deux fois plus fine en espace et en temps que sa
//! parente (même schéma qu'en 2D, voir refine.rs). Couplage par les moments
//! (ρ, u, partie hors équilibre de Π remise à l'échelle, ν̃) : la grille fine
//! a une enveloppe de mailles fantômes remplies par interpolation cubique
//! (dans le plan de chaque face) et linéaire en temps ; sous la zone
//! raffinée, la grille parente reçoit l'état de la grille fine. Les efforts
//! viennent de la grille la plus fine, l'affichage de la veine.

use std::fmt::Write as _;

use bytemuck::{Pod, Zeroable};

use crate::gpu3d::{f, w3, Gpu3D, C3};

/// Maille fantôme d'une grille fine (type de maille).
pub const GHOST: u32 = 0xffff_fffe;
/// Interpolation des fantômes : fantômes au plus par groupe de travail, et
/// nœuds parents au plus (moments en mémoire partagée, 48 o par nœud).
/// Essais (couche ×16 de `ar4c3`) : 64/192 et 128/288 à 3 %, 64/64 quatre
/// fois plus lent (un pochoir cubique remplit seul un groupe).
const GWG: usize = 128;
const GNODES: usize = 256;
/// Marge (mailles parentes) entre le bord de la grille fine et la zone recopiée.
const MARGIN: usize = 2;

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct LinkParams {
    pnx: u32,
    pny: u32,
    pnz: u32,
    pn: u32,
    cnx: u32,
    cny: u32,
    cnz: u32,
    cn: u32,
    ox: u32,
    oy: u32,
    oz: u32,
    ng: u32,
    rx: u32,
    ry: u32,
    rz: u32,
    margin: u32,
    sa: u32,
    alpha: f32,
    nu_p: f32,
    nu_c: f32,
    /// Grille fine en stockage creux : blocs selon x, y, z (0 : boîte) et
    /// début des coordonnées des blocs dans sa table.
    cbnx: u32,
    cbny: u32,
    cbnz: u32,
    cbcoord: u32,
    /// Idem pour la grille parente (couches pariétales imbriquées).
    pbnx: u32,
    pbny: u32,
    pbnz: u32,
    pbcoord: u32,
    /// Traitement pariétal du tirage (`pre_at`) : 0, rebond à mi-chemin ;
    /// 1, liens de Bouzidi et loi de paroi de la grille parente.
    pwm: u32,
    /// Propagation en place : parité du dernier pas de la grille parente,
    /// de la grille fine ; phase de la recopie (0 : mailles fluides, après
    /// le premier pas fin ; 1 : les autres, après le second).
    tp: u32,
    tc: u32,
    phase: u32,
}

const LINK_WGSL: &str = "
struct L {
  pnx: u32, pny: u32, pnz: u32, pn: u32,
  cnx: u32, cny: u32, cnz: u32, cn: u32,
  ox: u32, oy: u32, oz: u32, ng: u32,
  rx: u32, ry: u32, rz: u32, margin: u32,
  sa: u32, alpha: f32, nu_p: f32, nu_c: f32,
  cbnx: u32, cbny: u32, cbnz: u32, cbcoord: u32,
  pbnx: u32, pbny: u32, pbnz: u32, pbcoord: u32,
  pwm: u32, tp: u32, tc: u32, phase: u32,
};
@group(0) @binding(0) var<uniform> P: L;
const CV1_3 = 357.911;
fn nu_turb(v: f32, nu: f32) -> f32 {
  if (v <= 0.0) { return 0.0; }
  let chi = v / nu;
  let c3 = chi * chi * chi;
  return v * c3 / (c3 + CV1_3);
}
// k-ε RNG (P.sa = 2) : ν_t = C_μ k²/ε, comme `KE_FNS` de gpu3d.
const CMU = 0.0845;
fn nu_ke(k: f32, e: f32) -> f32 {
  if (k <= 0.0 || e <= 0.0) { return 0.0; }
  return min(0.08, CMU * k * k / e);
}
// Avec k-ε, b.w = k et c.w = ε (grille fine) ; τ s'en déduit.
fn tau_of(b: vec4<f32>, c: vec4<f32>) -> f32 {
  if (P.sa == 2u) { return 0.5 + 3.0 * (P.nu_c + nu_ke(b.w, c.w)); }
  return c.w;
}
// a = (ρ, u), b = (Π_xx, Π_yy, Π_zz hors équilibre, ν̃), c = (Π_xy, Π_xz, Π_yz, τ).
struct Mom {
  a: vec4<f32>,
  b: vec4<f32>,
  c: vec4<f32>,
};
fn flat(gid: vec3<u32>, nwg: vec3<u32>) -> u32 {
  return gid.x + gid.y * nwg.x * 64u;
}
";

/// Constantes D3Q19 et fonctions communes (moments, populations régularisées).
fn common() -> String {
    let list = |v: Vec<String>| v.join(", ");
    let cx = list((0..19).map(|i| C3[i][0].to_string()).collect());
    let cy = list((0..19).map(|i| C3[i][1].to_string()).collect());
    let cz = list((0..19).map(|i| C3[i][2].to_string()).collect());
    let w = list((0..19).map(|i| f(w3(i))).collect());
    let op = list(
        (0..19)
            .map(|i| {
                let c = C3[i];
                C3.iter().position(|d| d[0] == -c[0] && d[1] == -c[1] && d[2] == -c[2]).unwrap().to_string()
            })
            .collect(),
    );
    let hs = crate::gpu3d::hs_wgsl();
    let mut regul = String::new();
    for i in 0..19 {
        let [x, y, z] = C3[i].map(|v| v as f64);
        let cu: Vec<String> = [(x, "u.x"), (y, "u.y"), (z, "u.z")]
            .iter()
            .filter(|(c, _)| *c != 0.0)
            .map(|(c, s)| if *c > 0.0 { s.to_string() } else { format!("-{s}") })
            .collect();
        let cu = if cu.is_empty() { "0.0".into() } else { cu.join(" + ").replace("+ -", "- ") };
        let qp = format!(
            "{} * bxx + {} * byy + {} * bzz + {} * bxy + {} * bxz + {} * byz",
            f(x * x - 1.0 / 3.0),
            f(y * y - 1.0 / 3.0),
            f(z * z - 1.0 / 3.0),
            f(2.0 * x * y),
            f(2.0 * x * z),
            f(2.0 * y * z)
        );
        let h = |a: f64, b: f64| a * (b * b - 1.0 / 3.0);
        let h3: Vec<String> = [
            (h(y, x), "axxy"),
            (h(x, y), "axyy"),
            (h(z, x), "axxz"),
            (h(x, z), "axzz"),
            (h(z, y), "ayyz"),
            (h(y, z), "ayzz"),
        ]
        .iter()
        .filter(|(c, _)| c.abs() > 1e-12)
        .map(|(c, v)| format!("{} * {v}", f(*c)))
        .collect();
        let h3 = if h3.is_empty() { "0.0".into() } else { h3.join(" + ") };
        let wi = w3(i);
        let _ = writeln!(
            regul,
            "  {{ let cu = {cu}; p[{i}] = {} * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - usq) + {} * ({qp}) + {} * ({h3}); }}",
            f(wi),
            f(4.5 * wi),
            f(13.5 * wi)
        );
    }
    format!(
        "{LINK_WGSL}
const CX = array<i32, 19>({cx});
const CY = array<i32, 19>({cy});
const CZ = array<i32, 19>({cz});
const W = array<f32, 19>({w});
const OP = array<u32, 19>({op});
{hs}
fn moments(f: array<f32, 19>) -> Mom {{
  var rho = 0.0;
  var j = vec3<f32>(0.0);
  var pxx = 0.0; var pyy = 0.0; var pzz = 0.0; var pxy = 0.0; var pxz = 0.0; var pyz = 0.0;
  for (var i = 0u; i < 19u; i++) {{
    let c = vec3<f32>(f32(CX[i]), f32(CY[i]), f32(CZ[i]));
    let v = f[i];
    rho += v;
    j += c * v;
    pxx += c.x * c.x * v; pyy += c.y * c.y * v; pzz += c.z * c.z * v;
    pxy += c.x * c.y * v; pxz += c.x * c.z * v; pyz += c.y * c.z * v;
  }}
  let u = j / rho;
  return Mom(
    vec4<f32>(rho, u),
    vec4<f32>(pxx - rho * (1.0 / 3.0 + u.x * u.x), pyy - rho * (1.0 / 3.0 + u.y * u.y), pzz - rho * (1.0 / 3.0 + u.z * u.z), 0.0),
    vec4<f32>(pxy - rho * u.x * u.y, pxz - rho * u.x * u.z, pyz - rho * u.y * u.z, 0.0));
}}

// Populations régularisées récursives (comme la collision) pour ρ, u et la
// partie hors équilibre (b) d'ordre 2.
fn regul(rho: f32, u: vec3<f32>, bxx: f32, byy: f32, bzz: f32, bxy: f32, bxz: f32, byz: f32) -> array<f32, 19> {{
  let usq = 1.5 * dot(u, u);
  let axxy = rho * u.x * u.x * u.y + u.y * bxx + 2.0 * u.x * bxy;
  let axyy = rho * u.x * u.y * u.y + u.x * byy + 2.0 * u.y * bxy;
  let axxz = rho * u.x * u.x * u.z + u.z * bxx + 2.0 * u.x * bxz;
  let axzz = rho * u.x * u.z * u.z + u.x * bzz + 2.0 * u.z * bxz;
  let ayyz = rho * u.y * u.y * u.z + u.z * byy + 2.0 * u.y * byz;
  let ayzz = rho * u.y * u.z * u.z + u.y * bzz + 2.0 * u.z * byz;
  var p: array<f32, 19>;
{regul}  return p;
}}
"
    )
}

/// Accès aux populations `name` (deux tampons A et B, plans de `n` mailles),
/// stockées g = f − w en simple ou demi-précision ; `ld_` rend f.
fn access(name: &str, ba: u32, bb: u32, fp16: bool, p16: bool, write: bool, n: &str) -> String {
    let mode = if write { "read_write" } else { "read" };
    let pt = if p16 { "u16" } else { "u32" };
    let mut s = format!(
        "@group(0) @binding({ba}) var<storage, {mode}> {name}A: array<{pt}>;
@group(0) @binding({bb}) var<storage, {mode}> {name}B: array<{pt}>;
"
    );
    if p16 {
        s += &format!(
            "fn ld_{name}(i: u32, c: u32) -> f32 {{
  if (i < 10u) {{ return unpack2x16float(u32({name}A[i * {n} + c])).x + W[i]; }}
  return unpack2x16float(u32({name}B[(i - 10u) * {n} + c])).x + W[i];
}}
"
        );
        if write {
            s += &format!("fn st_{name}(c: u32, p: array<f32, 19>) {{
");
            for i in 0..19 {
                let (buf, k) = if i < 10 { ("A", i) } else { ("B", i - 10) };
                s += &format!("  {name}{buf}[{k}u * {n} + c] = u16(pack2x16float(vec2<f32>(p[{i}] - W[{i}], 0.0)) & 0xffffu);
");
            }
            s += "}
";
        }
    } else if fp16 {
        s += &format!(
            "fn ld_{name}(i: u32, c: u32) -> f32 {{
  let h = HS[i];
  let w = h >> 1u;
  var v: vec2<f32>;
  if (w < 5u) {{ v = unpack2x16float({name}A[w * {n} + c]); }} else {{ v = unpack2x16float({name}B[(w - 5u) * {n} + c]); }}
  return select(v.x, v.y, (h & 1u) == 1u) + W[i];
}}
"
        );
        if write {
            s += &format!("fn st_{name}(c: u32, p: array<f32, 19>) {{\n");
            let hd = crate::gpu3d::hdir();
            let pw = |h: usize| hd[h].map_or("0.0".to_string(), |i| format!("p[{i}] - W[{i}]"));
            for wd in 0..10 {
                let (buf, k) = if wd < 5 { ("A", wd) } else { ("B", wd - 5) };
                s += &format!("  {name}{buf}[{k}u * {n} + c] = pack2x16float(vec2<f32>({}, {}));\n", pw(2 * wd), pw(2 * wd + 1));
            }
            s += "}\n";
        }
    } else {
        s += &format!(
            "fn ld_{name}(i: u32, c: u32) -> f32 {{
  if (i < 10u) {{ return bitcast<f32>({name}A[i * {n} + c]) + W[i]; }}
  return bitcast<f32>({name}B[(i - 10u) * {n} + c]) + W[i];
}}
"
        );
        if write {
            s += &format!("fn st_{name}(c: u32, p: array<f32, 19>) {{\n");
            for i in 0..19 {
                let (buf, k) = if i < 10 { ("A", i) } else { ("B", i - 10) };
                s += &format!("  {name}{buf}[{k}u * {n} + c] = bitcast<u32>(p[{i}] - W[{i}]);\n");
            }
            s += "}\n";
        }
    }
    s
}

/// Accès par coordonnées : `LD{name}(i, x, y, z)` (population i sortante de la
/// maille (x, y, z) au dernier pas) et `ST{name}(x, y, z, rang, p)`. Deux jeux :
/// rang par `idx` ; en place (voir `gpu3d::pull_ep`) : population positive
/// rangée chez la voisine x + c_i, plan selon la parité `par`. `idx` prend des
/// coordonnées i32 (repliées dans la boîte s'il le faut).
fn access_xyz(name: &str, ep: bool, p16: bool, write: bool, n: &str, idx: &str, par: &str) -> String {
    if !ep {
        // Fantômes (écriture seule) : pas de lecture, ni de rang par coordonnées.
        let mut s = if write && idx == "cidw" {
            String::new()
        } else {
            format!(
                "fn LD{name}(i: u32, x: i32, y: i32, z: i32) -> f32 {{ return ld_{name}(i, {idx}(x, y, z)); }}
"
            )
        };
        if write {
            s += &format!(
                "fn ST{name}(x: i32, y: i32, z: i32, c: u32, p: array<f32, 19>) {{ st_{name}(c, p); }}
"
            );
        }
        return s;
    }
    let raw = if p16 {
        format!(
            "  if (pl < 10u) {{ return unpack2x16float(u32({name}A[pl * {n} + c])).x + W[i]; }}
  return unpack2x16float(u32({name}B[(pl - 10u) * {n} + c])).x + W[i];"
        )
    } else {
        format!(
            "  if (pl < 10u) {{ return bitcast<f32>({name}A[pl * {n} + c]) + W[i]; }}
  return bitcast<f32>({name}B[(pl - 10u) * {n} + c]) + W[i];"
        )
    };
    let mut s = format!(
        "fn LD{name}(i: u32, x: i32, y: i32, z: i32) -> f32 {{
  var c = {idx}(x, y, z);
  if ((i & 1u) == 1u) {{ c = {idx}(x + CX[i], y + CY[i], z + CZ[i]); }}
  let pl = select(i, OP[i], {par} == 1u);
{raw}
}}
"
    );
    if write {
        let put = if p16 {
            format!(
                "    let v = u16(pack2x16float(vec2<f32>(p[i] - W[i], 0.0)) & 0xffffu);
    if (pl < 10u) {{ {name}A[pl * {n} + c] = v; }} else {{ {name}B[(pl - 10u) * {n} + c] = v; }}"
            )
        } else {
            format!(
                "    let v = bitcast<u32>(p[i] - W[i]);
    if (pl < 10u) {{ {name}A[pl * {n} + c] = v; }} else {{ {name}B[(pl - 10u) * {n} + c] = v; }}"
            )
        };
        s += &format!(
            "fn ST{name}(x: i32, y: i32, z: i32, c0: u32, p: array<f32, 19>) {{
  for (var i = 0u; i < 19u; i++) {{
    var c = c0;
    if ((i & 1u) == 1u) {{ c = {idx}(x + CX[i], y + CY[i], z + CZ[i]); }}
    let pl = select(i, OP[i], {par} == 1u);
{put}
  }}
}}
"
        );
    }
    s
}

/// Tirage, première passe : moments avant collision de chaque nœud parent
/// des pochoirs (liste `pnode`, chaque nœud une seule fois ; un nœud sert
/// à jusqu'à 64 mailles fantômes).
const PRE: &str = "
@group(0) @binding(1) var<storage, read_write> pm: array<vec4<f32>>;
@group(0) @binding(2) var<storage, read> pnode: array<i32>;
@group(0) @binding(5) var<storage, read> pnut: array<f32>;
@group(0) @binding(6) var<storage, read> pwall: array<u32>;
@group(0) @binding(8) var<storage, read> pblk: array<u32>;
@group(0) @binding(9) var<storage, read> pq: array<f32>;
@group(0) @binding(10) var<storage, read> putau: array<f32>;
// Rang d'une maille de la grille parente (boîte, ou couche par blocs).
fn pcell(x: i32, y: i32, z: i32) -> u32 {
  let q = vec3<u32>(clamp(vec3<i32>(x, y, z), vec3<i32>(0), vec3<i32>(i32(P.pnx), i32(P.pny), i32(P.pnz)) - vec3<i32>(1)));
  if (P.pbnx == 0u) { return q.x + P.pnx * (q.y + P.pny * q.z); }
  let s = pblk[(q.x >> 3u) + P.pbnx * ((q.y >> 3u) + P.pbny * (q.z >> 3u))];
  return s * 512u + (q.x & 7u) + 8u * ((q.y & 7u) + 8u * (q.z & 7u));
}

// Moments avant collision au nœud parent, ramenés à la grille fine. Une
// population tirée d'un nœud solide (grille fine qui coupe la maquette)
// est renvoyée par la paroi : population opposée du nœud lui-même (rebond à
// mi-chemin), au lieu du contenu sans signification du nœud solide.
fn oct_dec(k: f32) -> vec3<f32> {
  let iu = floor(k / 2048.0);
  let iv = k - iu * 2048.0;
  var u = iu / 2047.0 * 2.0 - 1.0;
  var v = iv / 2047.0 * 2.0 - 1.0;
  let z = 1.0 - abs(u) - abs(v);
  if (z < 0.0) {
    let a = u;
    u = (1.0 - abs(v)) * sign(a);
    v = (1.0 - abs(a)) * sign(v);
  }
  return normalize(vec3<f32>(u, v, z));
}

// `solid` : bit i mis si le voisin d'où vient la population i est solide
// (précalculé, géométrie fixe).
fn pre_f(x: i32, y: i32, z: i32, solid: u32) -> Mom {
  var f: array<f32, 19>;
  let here = pcell(x, y, z);
  // Maille pariétale de la grille parente (avec `pwm`) : populations venues
  // de la paroi comme dans son propre pas (liens de Bouzidi, paroi qui
  // glisse à la vitesse de la loi de paroi), et viscosité de sa couche
  // pariétale ; sinon rebond à mi-chemin sur paroi fixe.
  let wid = pwall[here];
  let wm = P.pwm == 1u && wid != 0u && wid < 0xfffffffdu;
  let wk = select(0u, wid - 1u, wm);
  var uw = vec3<f32>(0.0);
  var vn = 0.0;
  if (wm) {
    uw = vec3<f32>(putau[4u * wk + 1u], putau[4u * wk + 2u], putau[4u * wk + 3u]);
    vn = putau[arrayLength(&putau) / 9u * 4u + wk];
  }
  for (var i = 0u; i < 19u; i++) {
    var q = 0.0;
    if (wm && i > 0u) { q = pq[wk * 18u + i - 1u]; }
    if (q > 0.0) {
      let j = OP[i];
      let fj = LDpf(j, x, y, z);
      let qo = pq[wk * 18u + j - 1u];
      let nf = oct_dec(pq[arrayLength(&pq) / 2u + wk * 18u + i - 1u]);
      let ci = vec3<f32>(f32(CX[i]), f32(CY[i]), f32(CZ[i]));
      var uwl = uw - dot(uw, nf) * nf + vn * nf;
      // Paroi virtuelle à q = 1 sous loi de paroi (`wall_q1` de gpu3d).
      let nuw = putau[arrayLength(&putau) / 9u * 5u + wk];
      let q1 = WALLQ1_C != 0u && nuw >= 0.0;
      if (q1) {
        let o9 = arrayLength(&putau) / 9u * 6u + 3u * wk;
        let g = vec3<f32>(putau[o9], putau[o9 + 1u], putau[o9 + 2u]);
        uwl -= (g - dot(g, nf) * nf) * (max(0.0, dot(ci, nf)) * (1.0 - q));
      }
      let corr = 6.0 * W[i] * dot(ci, uwl);
      if (q1) {
        f[i] = 0.5 * (fj + LDpf(i, x, y, z) + corr);
      } else if (q < 0.5 && qo == 0.0) {
        f[i] = 2.0 * q * fj + (1.0 - 2.0 * q) * LDpf(j, x + CX[i], y + CY[i], z + CZ[i]) + corr;
      } else if (q < 0.5) {
        f[i] = fj + corr;
      } else {
        f[i] = (fj + (2.0 * q - 1.0) * LDpf(i, x, y, z) + corr) / (2.0 * q);
      }
    } else if (((solid >> i) & 1u) != 0u) {
      f[i] = LDpf(OP[i], x, y, z);
    } else {
      f[i] = LDpf(i, x - CX[i], y - CY[i], z - CZ[i]);
    }
  }
  return moments(f);
}

// Mise à l'échelle de la grille fine (ν̃ ou k, ε de la grille parente après
// son pas, viscosité de la couche pariétale).
fn pre_scale(m0: Mom, x: i32, y: i32, z: i32) -> Mom {
  var m = m0;
  let here = pcell(x, y, z);
  let wid = pwall[here];
  let wm = P.pwm == 1u && wid != 0u && wid < 0xfffffffdu;
  let wk = select(0u, wid - 1u, wm);
  var v = 3.0 * P.nu_p;
  var vt = 0.0;
  var ef = 0.0;
  if (P.sa == 1u) {
    v = pnut[pcell(x, y, z)];
    vt = nu_turb(v, P.nu_p);
  }
  // Couche pariétale de la loi de paroi : viscosité ν_w de la collision.
  if (wm && P.sa == 1u) {
    let nuw = putau[arrayLength(&putau) / 9u * 5u + wk];
    if (nuw >= 0.0) { vt = nuw; }
  }
  if (P.sa == 2u) {
    // k inchangé d'une grille à l'autre (Δt ∝ Δx), ε réseau divisé par 2.
    let pc = pcell(x, y, z);
    v = 0.5 * pnut[2u * pc];
    ef = 0.5 * pnut[2u * pc + 1u];
    vt = nu_ke(2.0 * v, 2.0 * ef);
  }
  let tau_p = 0.5 + 3.0 * (P.nu_p + vt);
  let tau_c = 0.5 + 3.0 * (P.nu_c + 2.0 * vt);
  let k = tau_c / (2.0 * tau_p);
  m.b = vec4<f32>(m.b.xyz * k, 2.0 * v);
  m.c = vec4<f32>(m.c.xyz * k, select(tau_c, ef, P.sa == 2u));
  return m;
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let k = flat(gid, nwg);
  if (k >= arrayLength(&pnode) / 4u) { return; }
  let x = pnode[4u * k];
  let y = pnode[4u * k + 1u];
  let z = pnode[4u * k + 2u];
  // Deux jeux : populations d'avant le pas de la grille parente (l'autre
  // jeu) et ν̃ d'après. En place, l'état d'avant le pas est écrasé : moments
  // relevés avant le pas (phase 0), mis à l'échelle après (phase 1).
  var m: Mom;
  if (EP_PRE && P.phase == 0u) {
    m = pre_f(x, y, z, u32(pnode[4u * k + 3u]));
  } else if (EP_PRE) {
    m = pre_scale(Mom(pm[3u * k], pm[3u * k + 1u], pm[3u * k + 2u]), x, y, z);
  } else {
    m = pre_scale(pre_f(x, y, z, u32(pnode[4u * k + 3u])), x, y, z);
  }
  pm[3u * k] = m.a;
  pm[3u * k + 1u] = m.b;
  pm[3u * k + 2u] = m.c;
}
";

/// Tirage, seconde passe : interpolation des moments des nœuds parents
/// (`pm`) sur chaque maille fantôme, nœuds du pochoir dans l'ordre de
/// `gslot` (u32::MAX : nœud solide sauté), pochoir linéaire si le bit 31
/// de `gstart` est mis. Mêmes opérations, dans le même ordre, que le calcul
/// direct d'autrefois (moments de chaque nœud recalculés pour chaque maille
/// fantôme) : résultats identiques au bit près.
const GATHER: &str = "
// Moments interpolés : nouveaux (écrits ici), anciens (tirage précédent),
// deux tampons en alternance.
@group(0) @binding(1) var<storage, read_write> gm: array<vec4<f32>>;
@group(0) @binding(13) var<storage, read> gmo: array<vec4<f32>>;
@group(0) @binding(2) var<storage, read> pm: array<vec4<f32>>;
@group(0) @binding(7) var<storage, read> gpos: array<u32>;
@group(0) @binding(9) var<storage, read> gstart: array<u32>;
@group(0) @binding(10) var<storage, read> gslot: array<u32>;
// Remplissage du premier demi-pas de la grille fine (α = 0), fait ici.
@group(0) @binding(8) var<storage, read> gidx: array<u32>;
@group(0) @binding(5) var<storage, read_write> cnut: array<f32>;
@group(0) @binding(6) var<storage, read_write> cmacro: array<vec4<f32>>;

struct St { n: i32, o: vec4<i32>, w: vec4<f32>, };
// Pochoir le long d'un axe : nœud parent (coordonnée paire), ou interpolation
// cubique sur quatre nœuds (impaire).
fn stencil(x2: i32) -> St {
  if ((x2 & 1) == 0) {
    return St(1, vec4<i32>(x2 / 2, 0, 0, 0), vec4<f32>(1.0, 0.0, 0.0, 0.0));
  }
  let i0 = (x2 - 1) / 2;
  return St(4, vec4<i32>(i0 - 1, i0, i0 + 1, i0 + 2), vec4<f32>(-0.0625, 0.5625, 0.5625, -0.0625));
}

// Pochoir linéaire le long d'un axe (près d'une paroi).
fn lin(x2: i32) -> St {
  if ((x2 & 1) == 0) {
    return St(1, vec4<i32>(x2 / 2, 0, 0, 0), vec4<f32>(1.0, 0.0, 0.0, 0.0));
  }
  let i0 = (x2 - 1) / 2;
  return St(2, vec4<i32>(i0, i0 + 1, 0, 0), vec4<f32>(0.5, 0.5, 0.0, 0.0));
}

@group(0) @binding(11) var<storage, read> tile: array<u32>;
@group(0) @binding(12) var<storage, read> tnode: array<u32>;
// Moments des nœuds parents du groupe.
var<workgroup> sm: array<vec4<f32>, GNODES3>;

@compute @workgroup_size(GWG)
fn main(@builtin(workgroup_id) wid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>, @builtin(local_invocation_index) li: u32) {
  let t = wid.x + wid.y * nwg.x;
  if (t >= arrayLength(&tile) / 2u - 1u) { return; }
  let n0 = tile[2u * t + 1u];
  for (var e = n0 + li; e < tile[2u * t + 3u]; e += GWGu) {
    let l = 3u * (e - n0);
    let s = 3u * tnode[e];
    sm[l] = pm[s];
    sm[l + 1u] = pm[s + 1u];
    sm[l + 2u] = pm[s + 2u];
  }
  workgroupBarrier();
  for (var g = tile[2u * t] + li; g < tile[2u * t + 2u]; g += GWGu) {
    ghost(g);
  }
}

// Rang du nœud `e` dans la liste du groupe : un octet, quatre par mot
// (0xff : nœud solide sauté).
fn slot(e: u32) -> u32 {
  return (gslot[e >> 2u] >> (8u * (e & 3u))) & 0xffu;
}

fn ghost(g: u32) {
  let oa = gmo[3u * g];
  let ob = gmo[3u * g + 1u];
  let oc = gmo[3u * g + 2u];
  // Position de la maille fantôme dans la boîte de la grille fine.
  let c = gpos[g];
  let x = c % P.cnx;
  let y = (c / P.cnx) % P.cny;
  let z = c / (P.cnx * P.cny);
  let sx = stencil(i32(2u * P.ox + x));
  let sy = stencil(i32(2u * P.oy + y));
  let sz = stencil(i32(2u * P.oz + z));
  var m = Mom(vec4<f32>(0.0), vec4<f32>(0.0), vec4<f32>(0.0));
  // Pochoir cubique s'il ne touche aucun nœud solide ; sinon linéaire sur
  // les seuls nœuds fluides, poids renormalisés (grille fine qui coupe la
  // maquette : l'interpolation ne lit pas l'intérieur du corps).
  let st = gstart[g];
  var ptr = st & 0x7fffffffu;
  if ((st >> 31u) == 0u) {
    for (var k = 0; k < sz.n; k++) {
      for (var j = 0; j < sy.n; j++) {
        for (var i = 0; i < sx.n; i++) {
          let w = sx.w[i] * sy.w[j] * sz.w[k];
          let s = 3u * slot(ptr);
          ptr++;
          m.a += w * sm[s];
          m.b += w * sm[s + 1u];
          m.c += w * sm[s + 2u];
        }
      }
    }
  } else {
    let lx = lin(i32(2u * P.ox + x));
    let ly = lin(i32(2u * P.oy + y));
    let lz = lin(i32(2u * P.oz + z));
    var ws = 0.0;
    for (var k = 0; k < lz.n; k++) {
      for (var j = 0; j < ly.n; j++) {
        for (var i = 0; i < lx.n; i++) {
          let sl = slot(ptr);
          ptr++;
          if (sl == 0xffu) { continue; }
          let s = 3u * sl;
          let w = lx.w[i] * ly.w[j] * lz.w[k];
          m.a += w * sm[s];
          m.b += w * sm[s + 1u];
          m.c += w * sm[s + 2u];
          ws += w;
        }
      }
    }
    if (ws > 1e-6) {
      m.a /= ws;
      m.b /= ws;
      m.c /= ws;
    } else {
      // Maille fantôme dans le corps : état au repos (jamais tiré, les liens
      // qui y mènent sont des liens de paroi).
      m = Mom(vec4<f32>(1.0, 0.0, 0.0, 0.0), vec4<f32>(0.0), vec4<f32>(0.0, 0.0, 0.0, 0.5 + 3.0 * P.nu_c));
    }
  }
  gm[3u * g] = m.a;
  gm[3u * g + 1u] = m.b;
  gm[3u * g + 2u] = m.c;
  // Remplissage à α = P.alpha (0) : mêmes opérations que la passe FILL.
  let a = mix(oa, m.a, P.alpha);
  let b = mix(ob, m.b, P.alpha);
  let c3 = mix(oc, m.c, P.alpha);
  let kk = 1.0 - 1.0 / tau_of(b, c3);
  let p = regul(a.x, a.yzw, kk * b.x, kk * b.y, kk * b.z, kk * c3.x, kk * c3.y, kk * c3.z);
  let idx = gidx[g];
  STcf(i32(x), i32(y), i32(z), idx, p);
  if (P.sa == 1u) { cnut[idx] = b.w; }
  if (P.sa == 2u) {
    cnut[2u * idx] = max(b.w, 0.0);
    cnut[2u * idx + 1u] = max(c3.w, 0.0);
  }
  cmacro[idx] = vec4<f32>(a.yzw, a.x);
}
";

const FILL: &str = "
@group(0) @binding(1) var<storage, read> gmo: array<vec4<f32>>;
@group(0) @binding(7) var<storage, read> gm: array<vec4<f32>>;
@group(0) @binding(2) var<storage, read> gidx: array<u32>;
@group(0) @binding(5) var<storage, read_write> cnut: array<f32>;
@group(0) @binding(6) var<storage, read_write> cmacro: array<vec4<f32>>;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let g = flat(gid, nwg);
  if (g >= P.ng) { return; }
  let a = mix(gmo[3u * g], gm[3u * g], P.alpha);
  let b = mix(gmo[3u * g + 1u], gm[3u * g + 1u], P.alpha);
  let c = mix(gmo[3u * g + 2u], gm[3u * g + 2u], P.alpha);
  let k = 1.0 - 1.0 / tau_of(b, c);
  let p = regul(a.x, a.yzw, k * b.x, k * b.y, k * b.z, k * c.x, k * c.y, k * c.z);
  let idx = gidx[g];
  STcf(gx(g), gy(g), gz(g), idx, p);
  if (P.sa == 1u) { cnut[idx] = b.w; }
  if (P.sa == 2u) {
    cnut[2u * idx] = max(b.w, 0.0);
    cnut[2u * idx + 1u] = max(c.w, 0.0);
  }
  cmacro[idx] = vec4<f32>(a.yzw, a.x);
}
";

const RESTRICT: &str = "
@group(0) @binding(5) var<storage, read> cwall: array<u32>;
@group(0) @binding(6) var<storage, read> cnut: array<f32>;
@group(0) @binding(9) var<storage, read_write> pnut: array<f32>;
@group(0) @binding(10) var<storage, read_write> pmacro: array<vec4<f32>>;
@group(0) @binding(11) var<storage, read> pwall: array<u32>;
@group(0) @binding(12) var<storage, read> rlist: array<u32>;
@group(0) @binding(13) var<storage, read> cblk: array<u32>;
@group(0) @binding(14) var<storage, read> pblk: array<u32>;

fn pid(x: u32, y: u32, z: u32) -> u32 {
  if (P.pbnx == 0u) { return x + P.pnx * (y + P.pny * z); }
  let s = pblk[(x >> 3u) + P.pbnx * ((y >> 3u) + P.pbny * (z >> 3u))];
  return s * 512u + (x & 7u) + 8u * ((y & 7u) + 8u * (z & 7u));
}

// Rang d'une maille de la grille fine (boîte, ou stockage creux par blocs).
fn cid(x: u32, y: u32, z: u32) -> u32 {
  if (P.cbnx == 0u) { return x + P.cnx * (y + P.cny * z); }
  let s = cblk[(x >> 3u) + P.cbnx * ((y >> 3u) + P.cbny * (z >> 3u))];
  return s * 512u + (x & 7u) + 8u * ((y & 7u) + 8u * (z & 7u));
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let r = flat(gid, nwg);
  if (r >= P.rx) { return; }
  // Nœuds parents recopiés : liste (boîte moins sa marge, ou cœur de la
  // couche pariétale fine).
  // Liste : positions dans la boîte parente ; rang par `pid` (couche
  // parente par blocs). Pas de recopie dans un nœud solide ou fantôme.
  let ppos = rlist[r];
  let px = ppos % P.pnx;
  let py = (ppos / P.pnx) % P.pny;
  let pz = ppos / (P.pnx * P.pny);
  let pidx = pid(px, py, pz);
  if (pwall[pidx] >= 0xfffffffeu) { return; }
  let cx = 2u * (px - P.ox);
  let cy = 2u * (py - P.oy);
  let cz = 2u * (pz - P.oz);
  let cidx = cid(cx, cy, cz);
  EP_PHASE
  var v = 0.0;
  var vt = 0.0;
  var ev = 0.0;
  if (P.sa == 1u) {
    v = cnut[cidx];
    vt = nu_turb(v, P.nu_c);
  }
  if (P.sa == 2u) {
    v = cnut[2u * cidx];
    ev = cnut[2u * cidx + 1u];
    vt = nu_ke(v, ev);
  }
  var f: array<f32, 19>;
  var p: array<f32, 19>;
  var m: Mom;
  // Maille fluide, ou recouverte par une grille plus fine (même traitement).
  if (cwall[cidx] == 0u || cwall[cidx] == 0xfffffffdu) {
    // Populations avant collision au temps t + 1, tirées du demi-pas t + ½.
    EP_FLUID
    let tau_c = 0.5 + 3.0 * (P.nu_c + vt);
    let tau_p = 0.5 + 3.0 * (P.nu_p + 0.5 * vt);
    let k = 2.0 * (tau_p - 1.0) / tau_c;
    p = regul(m.a.x, m.a.yzw, k * m.b.x, k * m.b.y, k * m.b.z, k * m.c.x, k * m.c.y, k * m.c.z);
  } else if (cwall[cidx] < 0xfffffffeu) {
    // Maille pariétale de la grille fine (fluide, liens de paroi : pas de
    // tirage possible) : moments après collision, partie hors équilibre
    // ramenée avant collision (÷ (1 − 1/τ_c)), puis à la grille parente comme
    // pour les autres mailles. Une grille fine qui coupe la maquette rend
    // ainsi à la grille parente la contrainte de paroi, et pas un équilibre.
    for (var i = 0u; i < 19u; i++) {
      f[i] = LDcc(i, i32(cx), i32(cy), i32(cz));
    }
    m = moments(f);
    let tau_c = 0.5 + 3.0 * (P.nu_c + vt);
    let tau_p = 0.5 + 3.0 * (P.nu_p + 0.5 * vt);
    let post = 1.0 - 1.0 / tau_c;
    var k = 0.0;
    if (abs(post) > 0.1) { k = 2.0 * (tau_p - 1.0) / tau_c / post; }
    p = regul(m.a.x, m.a.yzw, k * m.b.x, k * m.b.y, k * m.b.z, k * m.c.x, k * m.c.y, k * m.c.z);
  } else {
    // Maille solide ou fantôme dans la grille fine : équilibre (affichage).
    for (var i = 0u; i < 19u; i++) {
      f[i] = LDcc(i, i32(cx), i32(cy), i32(cz));
    }
    m = moments(f);
    p = regul(m.a.x, m.a.yzw, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  }
  STpf(i32(px), i32(py), i32(pz), pidx, p);
  EP_SIDE
  if (P.sa == 1u) { pnut[pidx] = 0.5 * v; }
  if (P.sa == 2u) {
    pnut[2u * pidx] = v;
    pnut[2u * pidx + 1u] = 2.0 * ev;
  }
  pmacro[pidx] = vec4<f32>(m.a.yzw, m.a.x);
}
";

/// Couplage d'une grille (parente) et de la suivante (fine).
pub(crate) struct Link3 {
    ox: usize,
    oy: usize,
    oz: usize,
    ng: u32,
    /// Nœuds parents recopiés depuis la grille fine (et leurs positions).
    nr: u32,
    rl_pos: Vec<u32>,
    rlist: wgpu::Buffer,
    /// En place : moments des nœuds recopiés entre les deux pas fins.
    rmom: wgpu::Buffer,
    /// Paramètres, indice (pcur · 2 + ccur) · 2 + a : états des grilles
    /// parente et fine (parités en place), α = a / 2 (remplissage) ou phase
    /// a de la recopie.
    params: Vec<wgpu::Buffer>,
    /// Moments des fantômes aux deux derniers tirages, en alternance ;
    /// `gnew` : tampon du dernier.
    gm: [wgpu::Buffer; 2],
    gnew: std::cell::Cell<usize>,
    gidx: wgpu::Buffer,
    /// Position (indice de la boîte fine) de chaque maille fantôme.
    gpos: wgpu::Buffer,
    /// Tirage en deux passes : nœuds parents des pochoirs (x, y, z), leurs
    /// moments, et pour chaque maille fantôme le début de sa liste de nœuds
    /// (bit 31 : pochoir linéaire) dans `gslot`.
    nn: u32,
    pnode: wgpu::Buffer,
    pm: wgpu::Buffer,
    gstart: wgpu::Buffer,
    gslot: wgpu::Buffer,
    /// Tuiles de l'interpolation (voir `new`).
    nt: u32,
    tile: wgpu::Buffer,
    tnode: wgpu::Buffer,
    pre: wgpu::ComputePipeline,
    pre_bind: Vec<wgpu::BindGroup>,
    gather: wgpu::ComputePipeline,
    fill: wgpu::ComputePipeline,
    restrict: wgpu::ComputePipeline,
    gather_bind: Vec<wgpu::BindGroup>,
    fill_bind: Vec<wgpu::BindGroup>,
    restrict_bind: Vec<wgpu::BindGroup>,
}

/// Lancement de `n` groupes de travail (au plus 65535 par dimension).
fn disp_groups(n: u32) -> (u32, u32) {
    let groups = n.max(1);
    let dx = groups.min(65535);
    (dx, groups.div_ceil(dx))
}

fn disp(n: u32) -> (u32, u32) {
    let groups = n.div_ceil(64).max(1);
    let dx = groups.min(65535);
    (dx, groups.div_ceil(dx))
}

impl Link3 {
    /// Positions (boîte parente) des nœuds recopiés depuis la grille fine.
    pub(crate) fn restrict_list(&self) -> &[u32] {
        &self.rl_pos
    }

    pub(crate) fn new(parent: &Gpu3D, child: &Gpu3D, o: [usize; 3]) -> Link3 {
        let d = &parent.device;
        let fp16 = parent.fp16;
        let (cnx, cny, cnz) = (child.nx, child.ny, child.nz);
        let mut ghosts = Vec::new();
        let shell = child.band_fill.is_none();
        let mut gpos: Vec<u32> = Vec::new();
        if let (Some(f), Some(p)) = (child.band_fill.as_ref(), child.band_fill_pos.as_ref()) {
            ghosts = f.clone();
            gpos = p.clone();
        }
        for z in (0..cnz).filter(|_| shell) {
            for y in 0..cny {
                for x in 0..cnx {
                    if x == 0 || y == 0 || z == 0 || x == cnx - 1 || y == cny - 1 || z == cnz - 1 {
                        ghosts.push((x + cnx * (y + cny * z)) as u32);
                    }
                }
            }
        }
        if shell {
            gpos = ghosts.clone();
        }
        let ng = ghosts.len() as u32;
        // Pochoirs du tirage (comme `stencil` et `lin` du shader) : nœuds
        // parents distincts, et liste des nœuds de chaque maille fantôme.
        let pd = [parent.nx as i32, parent.ny as i32, parent.nz as i32];
        let prank = |p: [i32; 3]| -> usize {
            let q = [0, 1, 2].map(|a| p[a].clamp(0, pd[a] - 1) as usize);
            match parent.sparse.as_ref() {
                Some(m) => m.id(q[0], q[1], q[2]),
                None => q[0] + parent.nx * (q[1] + parent.ny * q[2]),
            }
        };
        let psolid = |p: [i32; 3]| parent.solid[prank(p)];
        let cubic = |x2: i32| -> Vec<i32> {
            if x2 & 1 == 0 {
                vec![x2 / 2]
            } else {
                let i0 = (x2 - 1) / 2;
                vec![i0 - 1, i0, i0 + 1, i0 + 2]
            }
        };
        let linear = |x2: i32| -> Vec<i32> {
            if x2 & 1 == 0 {
                vec![x2 / 2]
            } else {
                let i0 = (x2 - 1) / 2;
                vec![i0, i0 + 1]
            }
        };
        let mut node_of: std::collections::HashMap<[i32; 3], u32> = std::collections::HashMap::new();
        let mut pnode: Vec<i32> = Vec::new();
        let mut gstart: Vec<u32> = Vec::with_capacity(gpos.len());
        let mut gslot: Vec<u32> = Vec::new();
        for &c in &gpos {
            let c = c as usize;
            let x2 = [2 * o[0] + c % cnx, 2 * o[1] + (c / cnx) % cny, 2 * o[2] + c / (cnx * cny)].map(|v| v as i32);
            let s3 = x2.map(cubic);
            let mut solid = false;
            for &z in &s3[2] {
                for &y in &s3[1] {
                    for &x in &s3[0] {
                        solid |= psolid([x, y, z]);
                    }
                }
            }
            gstart.push(gslot.len() as u32 | if solid { 1 << 31 } else { 0 });
            let a3 = if solid { x2.map(linear) } else { s3 };
            for &z in &a3[2] {
                for &y in &a3[1] {
                    for &x in &a3[0] {
                        let p = [x, y, z];
                        if solid && psolid(p) {
                            gslot.push(u32::MAX);
                            continue;
                        }
                        let k = *node_of.entry(p).or_insert_with(|| {
                            // Voisins solides d'où viennent les populations
                            // (comme `pwall[src] == SOLID` du shader).
                            let mask = (0..19).filter(|&i| psolid([0, 1, 2].map(|a| p[a] - C3[i][a]))).fold(0u32, |m, i| m | 1 << i);
                            pnode.extend_from_slice(&p);
                            pnode.push(mask as i32);
                            (pnode.len() / 4 - 1) as u32
                        });
                        gslot.push(k);
                    }
                }
            }
        }
        let nn = (pnode.len() / 4) as u32;
        // Groupes de l'interpolation : au plus GWG fantômes consécutifs
        // (couche : rangés par blocs de 8³ ; boîte : par rangées) dont les
        // pochoirs touchent au plus GNODES nœuds parents. Fantômes [tile[2t], tile[2t + 2]), nœuds [tile[2t + 1],
        // tile[2t + 3]) de `tnode` (rangs dans `pm`) ; `gslot` passe en rangs
        // dans la liste du groupe, lue une fois en mémoire partagée.
        let mut tile: Vec<u32> = Vec::new();
        let mut tnode: Vec<u32> = Vec::new();
        let mut local: std::collections::HashMap<u32, u32> = std::collections::HashMap::new();
        let total = gslot.len();
        let ends = |g: usize| {
            let st = (gstart[g] & 0x7fff_ffff) as usize;
            let en = if g + 1 < gstart.len() { (gstart[g + 1] & 0x7fff_ffff) as usize } else { total };
            st..en
        };
        let mut count = 0usize;
        for g in 0..gpos.len() {
            let fresh = ends(g).filter(|&e| gslot[e] != u32::MAX && !local.contains_key(&gslot[e])).map(|e| gslot[e]).collect::<std::collections::BTreeSet<u32>>().len();
            // Rangs sur un octet, 0xff réservé : au plus 255 nœuds.
            if g == 0 || count == GWG || local.len() + fresh > GNODES.min(255) {
                tile.push(g as u32);
                tile.push(tnode.len() as u32);
                local.clear();
                count = 0;
            }
            count += 1;
            let n0 = tile[tile.len() - 1];
            for e in ends(g) {
                let k = gslot[e];
                if k == u32::MAX {
                    continue;
                }
                let l = *local.entry(k).or_insert_with(|| {
                    tnode.push(k);
                    tnode.len() as u32 - 1 - n0
                });
                gslot[e] = l;
            }
        }
        let nt = tile.len() as u32 / 2;
        tile.push(gpos.len() as u32);
        tile.push(tnode.len() as u32);
        if tnode.is_empty() {
            tnode.push(0);
        }
        if pnode.is_empty() {
            pnode.extend_from_slice(&[0, 0, 0, 0]);
        }
        if gstart.is_empty() {
            gstart.push(0);
        }
        // Rangs dans le groupe sur un octet (quatre par mot).
        let gslot: Vec<u32> = {
            let b: Vec<u8> = gslot.iter().map(|&k| if k == u32::MAX { 0xff } else { k as u8 }).collect();
            let mut w: Vec<u32> = b.chunks(4).map(|c| c.iter().enumerate().fold(0u32, |a, (i, &v)| a | (v as u32) << (8 * i))).collect();
            if w.is_empty() {
                w.push(0);
            }
            w
        };
        let base = common();
        let module = |label: &str, body: String| {
            crate::gpu3d::compute_module(d, label, format!("{}{base}{body}", crate::gpu3d::pop16_enable(parent.p16 || child.p16)))
        };
        let ep = parent.ep;
        let pre = module(
            "pre3d",
            format!(
                "const EP_PRE = {ep};\nconst WALLQ1_C = {}u;\n{}{}{PRE}",
                u32::from(crate::gpu3d::wall_q1()),
                access("pf", 3, 4, fp16, parent.p16, false, "P.pn"),
                access_xyz("pf", ep, parent.p16, false, "P.pn", "pcell", "P.tp")
            ),
        );
        // Tirage et remplissage du premier demi-pas de la grille fine en une
        // passe (un lancement de moins par pas de la grille parente).
        // En place : rang des voisines des fantômes dans la grille fine
        // (coordonnées repliées dans la boîte, stockage creux par blocs).
        let cidw = |blk: u32| {
            if ep {
                format!(
                    "@group(0) @binding({blk}) var<storage, read> cblk: array<u32>;
fn cidw(x0: i32, y0: i32, z0: i32) -> u32 {{
  let x = u32((x0 + i32(P.cnx)) % i32(P.cnx));
  let y = u32((y0 + i32(P.cny)) % i32(P.cny));
  let z = u32((z0 + i32(P.cnz)) % i32(P.cnz));
  if (P.cbnx == 0u) {{ return x + P.cnx * (y + P.cny * z); }}
  let s = cblk[(x >> 3u) + P.cbnx * ((y >> 3u) + P.cbny * (z >> 3u))];
  return s * 512u + (x & 7u) + 8u * ((y & 7u) + 8u * (z & 7u));
}}
"
                )
            } else {
                String::new()
            }
        };
        let gather = module(
            "gather3d",
            format!(
                "{}{}{}{}",
                access("cf", 3, 4, fp16, child.p16, true, "P.cn"),
                cidw(14),
                access_xyz("cf", ep, child.p16, true, "P.cn", "cidw", "P.tc"),
                GATHER.replace("GNODES3", &(3 * GNODES).to_string()).replace("GWG", &GWG.to_string())
            ),
        );
        let fill_pos = if ep {
            "@group(0) @binding(8) var<storage, read> gpos: array<u32>;
fn gx(g: u32) -> i32 { return i32(gpos[g] % P.cnx); }
fn gy(g: u32) -> i32 { return i32((gpos[g] / P.cnx) % P.cny); }
fn gz(g: u32) -> i32 { return i32(gpos[g] / (P.cnx * P.cny)); }
"
        } else {
            "fn gx(g: u32) -> i32 { return 0; }
fn gy(g: u32) -> i32 { return 0; }
fn gz(g: u32) -> i32 { return 0; }
"
        };
        let fill = module(
            "fill3d",
            format!(
                "{}{}{fill_pos}{}{FILL}",
                access("cf", 3, 4, fp16, child.p16, true, "P.cn"),
                cidw(9),
                access_xyz("cf", ep, child.p16, true, "P.cn", "cidw", "P.tc")
            ),
        );
        // En place : la recopie se fait en deux phases (voir `encode`) ; les
        // sortantes gardées des mailles pariétales de la grille parente
        // suivent les populations recopiées.
        let (phase, side) = if ep {
            (
                "if (P.phase == 0u && !(cwall[cidx] == 0u || cwall[cidx] == 0xfffffffdu)) { return; }",
                format!(
                    "{{
    let w = pwall[pidx];
    if (w != 0u && w < 0xfffffffdu) {{
      for (var k = 1u; k < 19u; k++) {{
        pside[(w - 1u) * 18u + k - 1u] = {};
      }}
    }}
  }}",
                    if parent.p16 { "unpack2x16float(pack2x16float(vec2<f32>(p[k] - W[k], 0.0)) & 0xffffu).x" } else { "p[k] - W[k]" }
                ),
            )
        } else {
            ("", String::new())
        };
        // Mailles fluides : moments avant collision au temps t + 1, tirés du
        // demi-pas ; en place, relevés après le premier pas fin (le second
        // les écrase) et gardés jusqu'à la recopie, après le second (ν̃ au
        // temps t + 1, comme à deux jeux).
        let pull = "for (var i = 0u; i < 19u; i++) {
      let s = vec3<i32>(i32(cx) - CX[i], i32(cy) - CY[i], i32(cz) - CZ[i]);
      f[i] = LDcs(i, s.x, s.y, s.z);
    }
    m = moments(f);";
        let fluid = if ep {
            format!(
                "if (P.phase == 0u) {{
      {pull}
      rmom[3u * r] = m.a;
      rmom[3u * r + 1u] = m.b;
      rmom[3u * r + 2u] = m.c;
      return;
    }}
    m = Mom(rmom[3u * r], rmom[3u * r + 1u], rmom[3u * r + 2u]);"
            )
        } else {
            pull.to_string()
        };
        let restrict_src = RESTRICT.replace("EP_PHASE", phase).replace("EP_SIDE", &side).replace("EP_FLUID", &fluid);
        let pidw = "fn pidw(x: i32, y: i32, z: i32) -> u32 { return pid(u32(x), u32(y), u32(z)); }
fn cidr(x: i32, y: i32, z: i32) -> u32 { return cid(u32(x), u32(y), u32(z)); }
";
        let restrict = module(
            "restrict3d",
            format!(
                "{}{}{}{}{}{}{}{}{restrict_src}",
                access("cs", 1, 2, fp16, child.p16, false, "P.cn"),
                access("cc", 3, 4, fp16, child.p16, false, "P.cn"),
                access("pf", 7, 8, fp16, parent.p16, true, "P.pn"),
                pidw,
                access_xyz("cs", ep, child.p16, false, "P.cn", "cidr", "P.tc"),
                access_xyz("cc", ep, child.p16, false, "P.cn", "cidr", "P.tc"),
                access_xyz("pf", ep, parent.p16, true, "P.pn", "pidw", "P.tp"),
                if ep {
                    "@group(0) @binding(15) var<storage, read_write> pside: array<f32>;
@group(0) @binding(16) var<storage, read_write> rmom: array<vec4<f32>>;
"
                } else {
                    ""
                },
            ),
        );
        let pipe = |m: &wgpu::ShaderModule, label: &str| {
            d.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some(label),
                layout: None,
                module: m,
                entry_point: Some("main"),
                compilation_options: Default::default(),
                cache: None,
            })
        };
        use wgpu::util::DeviceExt;
        let uni = wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST;
        let mk = |label: &str, size: u64, usage| {
            d.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label),
                size: size.max(16),
                usage,
                mapped_at_creation: false,
            })
        };
        let half = [(cnx - 1) / 2, (cny - 1) / 2, (cnz - 1) / 2];
        // Nœuds parents recopiés : la boîte moins sa marge ; avec une couche
        // pariétale, seulement ceux dont le nœud fin est loin de ses mailles
        // fantômes (`band_deep`).
        let (pnx, pny) = (parent.nx, parent.ny);
        let mut rl: Vec<u32> = Vec::new();
        for hz in MARGIN..=half[2] - MARGIN {
            for hy in MARGIN..=half[1] - MARGIN {
                for hx in MARGIN..=half[0] - MARGIN {
                    let keep = match child.band_deep.as_ref() {
                        Some(d) => d[child.sparse.as_ref().map_or(2 * hx + cnx * (2 * hy + cny * 2 * hz), |m| m.id(2 * hx, 2 * hy, 2 * hz))],
                        None => true,
                    };
                    if keep {
                        rl.push(((o[0] + hx) + pnx * ((o[1] + hy) + pny * (o[2] + hz))) as u32);
                    }
                }
            }
        }
        // Ordre de stockage de la grille parente (blocs de 8³ d'une couche) :
        // écritures des nœuds recopiés contiguës dans chaque bloc.
        if let Some(m) = parent.sparse.as_ref() {
            rl.sort_by_key(|&q| {
                let q = q as usize;
                m.id(q % pnx, (q / pnx) % pny, q / (pnx * pny))
            });
        }
        let nr = rl.len() as u32;
        let rl_pos = rl.clone();
        if rl.is_empty() {
            rl.push(0);
        }
        Link3 {
            rl_pos,
            ox: o[0],
            oy: o[1],
            oz: o[2],
            ng,
            nr,
            rmom: mk("restrict_mom3d", if ep { 48 * nr as u64 } else { 16 }, wgpu::BufferUsages::STORAGE),
            rlist: d.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("restrict_list3d"),
                contents: bytemuck::cast_slice(&rl),
                usage: wgpu::BufferUsages::STORAGE,
            }),
            params: (0..8).map(|_| mk("link3", std::mem::size_of::<LinkParams>() as u64, uni)).collect(),
            gm: [mk("ghost3d_0", 48 * ng as u64, wgpu::BufferUsages::STORAGE), mk("ghost3d_1", 48 * ng as u64, wgpu::BufferUsages::STORAGE)],
            gnew: std::cell::Cell::new(1),
            gidx: d.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("ghost_idx3d"),
                contents: bytemuck::cast_slice(&ghosts),
                usage: wgpu::BufferUsages::STORAGE,
            }),
            gpos: d.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("ghost_pos3d"),
                contents: bytemuck::cast_slice(&gpos),
                usage: wgpu::BufferUsages::STORAGE,
            }),
            nn,
            pm: mk("pre_mom3d", 48 * nn as u64, wgpu::BufferUsages::STORAGE),
            pnode: d.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("pre_node3d"),
                contents: bytemuck::cast_slice(&pnode),
                usage: wgpu::BufferUsages::STORAGE,
            }),
            gstart: d.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("ghost_start3d"),
                contents: bytemuck::cast_slice(&gstart),
                usage: wgpu::BufferUsages::STORAGE,
            }),
            gslot: d.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("ghost_slot3d"),
                contents: bytemuck::cast_slice(&gslot),
                usage: wgpu::BufferUsages::STORAGE,
            }),
            nt,
            tile: d.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("ghost_tile3d"),
                contents: bytemuck::cast_slice(&tile),
                usage: wgpu::BufferUsages::STORAGE,
            }),
            tnode: d.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("ghost_tile_node3d"),
                contents: bytemuck::cast_slice(&tnode),
                usage: wgpu::BufferUsages::STORAGE,
            }),
            pre: pipe(&pre, "pre3d"),
            pre_bind: Vec::new(),
            gather: pipe(&gather, "gather3d"),
            fill: pipe(&fill, "fill3d"),
            restrict: pipe(&restrict, "restrict3d"),
            gather_bind: Vec::new(),
            fill_bind: Vec::new(),
            restrict_bind: Vec::new(),
        }
    }

    /// Paramètres et groupes de liaison (à refaire si un tampon change).
    pub(crate) fn prepare(&mut self, p: &Gpu3D, c: &Gpu3D) {
        for (i, buf) in self.params.iter().enumerate() {
            let (pc, cc, a) = (i / 4, (i / 2) % 2, i % 2);
            let lp = LinkParams {
                pnx: p.nx as u32,
                pny: p.ny as u32,
                pnz: p.nz as u32,
                pn: p.n as u32,
                cnx: c.nx as u32,
                cny: c.ny as u32,
                cnz: c.nz as u32,
                cn: c.n as u32,
                ox: self.ox as u32,
                oy: self.oy as u32,
                oz: self.oz as u32,
                ng: self.ng,
                rx: self.nr,
                ry: 1,
                rz: 1,
                margin: MARGIN as u32,
                sa: if p.sa && c.sa { 1 + p.keps as u32 } else { 0 },
                alpha: 0.5 * a as f32,
                nu_p: p.nu as f32,
                nu_c: c.nu as f32,
                cbnx: c.sparse.as_ref().map_or(0, |m| m.bd()[0] as u32),
                cbny: c.sparse.as_ref().map_or(0, |m| m.bd()[1] as u32),
                cbnz: c.sparse.as_ref().map_or(0, |m| m.bd()[2] as u32),
                cbcoord: c.sparse.as_ref().map_or(0, |m| m.table_len() as u32),
                pbnx: p.sparse.as_ref().map_or(0, |m| m.bd()[0] as u32),
                pbny: p.sparse.as_ref().map_or(0, |m| m.bd()[1] as u32),
                pbnz: p.sparse.as_ref().map_or(0, |m| m.bd()[2] as u32),
                pbcoord: p.sparse.as_ref().map_or(0, |m| m.table_len() as u32),
                pwm: (crate::gpu3d::raccord_paroi() && p.sa && !p.keps && p.surf.is_none()) as u32,
                tp: 1 - pc as u32,
                tc: 1 - cc as u32,
                phase: a as u32,
            };
            p.queue.write_buffer(buf, 0, bytemuck::bytes_of(&lp));
        }
        let d = &p.device;
        let group = |pipe: &wgpu::ComputePipeline, res: Vec<(u32, &wgpu::Buffer)>| {
            let entries: Vec<wgpu::BindGroupEntry> = res
                .into_iter()
                .map(|(binding, b)| wgpu::BindGroupEntry {
                    binding,
                    resource: b.as_entire_binding(),
                })
                .collect();
            d.create_bind_group(&wgpu::BindGroupDescriptor {
                label: None,
                layout: &pipe.get_bind_group_layout(0),
                entries: &entries,
            })
        };
        // Après un pas, `cur` désigne le nouvel état ; l'ancien fournit les tirages.
        // Indice cur · 2 + phase (en place : avant, après le pas parent).
        self.pre_bind = (0..4)
            .map(|i| {
                let cur = i / 2;
                group(
                    &self.pre,
                    vec![
                        (0, &self.params[cur * 4 + i % 2]),
                        (1, &self.pm),
                        (2, &self.pnode),
                        (3, &p.pops[1 - cur][0]),
                        (4, &p.pops[1 - cur][1]),
                        (5, &p.nut[cur]),
                        (6, &p.wall),
                        (8, &p.blk_buf),
                        (9, &p._qbuf),
                        (10, &p._utau),
                    ],
                )
            })
            .collect();
        // Indice cur · 2 + n : état de la grille fine, tampon des nouveaux moments.
        self.gather_bind = (0..4)
            .map(|i| {
                let (cur, n) = (i / 2, i % 2);
                let mut e = vec![
                        (0, &self.params[cur * 2]),
                        (1, &self.gm[n]),
                        (13, &self.gm[1 - n]),
                        (2, &self.pm),
                        (3, &c.pops[cur][0]),
                        (4, &c.pops[cur][1]),
                        (5, &c.nut[cur]),
                        (6, &c.macro_buf),
                        (7, &self.gpos),
                        (8, &self.gidx),
                        (9, &self.gstart),
                        (10, &self.gslot),
                        (11, &self.tile),
                        (12, &self.tnode),
                    ];
                if p.ep {
                    e.push((14, &c.blk_buf));
                }
                group(&self.gather, e)
            })
            .collect();
        // Indice (cur · 2 + a) · 2 + n.
        self.fill_bind = (0..8)
            .map(|i| {
                let (cur, a, n) = (i / 4, (i / 2) % 2, i % 2);
                let mut e = vec![
                        (0, &self.params[cur * 2 + a]),
                        (1, &self.gm[1 - n]),
                        (7, &self.gm[n]),
                        (2, &self.gidx),
                        (3, &c.pops[cur][0]),
                        (4, &c.pops[cur][1]),
                        (5, &c.nut[cur]),
                        (6, &c.macro_buf),
                    ];
                if p.ep {
                    e.push((8, &self.gpos));
                    e.push((9, &c.blk_buf));
                }
                group(&self.fill, e)
            })
            .collect();
        // Indice (pc · 2 + cc) · 2 + phase.
        self.restrict_bind = (0..8)
            .map(|i| {
                let (pc, cc) = (i / 4, (i / 2) % 2);
                let mut e = vec![
                        (0, &self.params[i]),
                        (1, &c.pops[1 - cc][0]),
                        (2, &c.pops[1 - cc][1]),
                        (3, &c.pops[cc][0]),
                        (4, &c.pops[cc][1]),
                        (5, &c.wall),
                        (6, &c.nut[cc]),
                        (7, &p.pops[pc][0]),
                        (8, &p.pops[pc][1]),
                        (9, &p.nut[pc]),
                        (10, &p.macro_buf),
                        (11, &p.wall),
                        (12, &self.rlist),
                        (13, &c.blk_buf),
                        (14, &p.blk_buf),
                    ];
                if p.ep {
                    e.push((15, &p.side));
                    e.push((16, &self.rmom));
                }
                group(&self.restrict, e)
            })
            .collect();
    }

    /// Profilage (banc de mesure) : durée de chaque passe du raccord.
    pub(crate) fn profile(
        &self,
        p: &Gpu3D,
        c: &Gpu3D,
        time: &dyn Fn(&Gpu3D, &dyn Fn(&mut wgpu::ComputePass)) -> f64,
    ) -> Vec<(String, usize, f64)> {
        let (dx, dy) = disp(self.ng);
        let (rx, ry) = disp(self.nr);
        vec![
            (
                "tirage, moments des nœuds".into(),
                self.nn as usize,
                time(p, &|pass| {
                    let (dx, dy) = disp(self.nn);
                    pass.set_pipeline(&self.pre);
                    pass.set_bind_group(0, &self.pre_bind[p.cur * 2], &[]);
                    pass.dispatch_workgroups(dx, dy, 1);
                }),
            ),
            (
                format!("tirage, interpolation et remplissage ({:.0} fantômes par groupe)", self.ng as f64 / self.nt.max(1) as f64),
                self.ng as usize,
                time(p, &|pass| {
                    let (dx, dy) = disp_groups(self.nt);
                    pass.set_pipeline(&self.gather);
                    pass.set_bind_group(0, &self.gather_bind[c.cur * 2], &[]);
                    pass.dispatch_workgroups(dx, dy, 1);
                }),
            ),
            (
                "remplissage".into(),
                self.ng as usize,
                time(p, &|pass| {
                    pass.set_pipeline(&self.fill);
                    pass.set_bind_group(0, &self.fill_bind[(c.cur * 2 + 1) * 2], &[]);
                    pass.dispatch_workgroups(dx, dy, 1);
                }),
            ),
            (
                "recopie".into(),
                self.nr as usize,
                time(p, &|pass| {
                    pass.set_pipeline(&self.restrict);
                    pass.set_bind_group(0, &self.restrict_bind[(p.cur * 2 + c.cur) * 2], &[]);
                    pass.dispatch_workgroups(rx, ry, 1);
                }),
            ),
        ]
    }

    /// Tirage (et remplissage du premier demi-pas de la grille fine, état `ccur`).
    fn gather(&self, pass: &mut wgpu::ComputePass, pcur: usize, ccur: usize, ep: bool) {
        if self.nn > 0 {
            let (dx, dy) = disp(self.nn);
            pass.set_pipeline(&self.pre);
            pass.set_bind_group(0, &self.pre_bind[pcur * 2 + ep as usize], &[]);
            pass.dispatch_workgroups(dx, dy, 1);
        }
        let n = 1 - self.gnew.get();
        self.gnew.set(n);
        let (dx, dy) = disp_groups(self.nt);
        pass.set_pipeline(&self.gather);
        pass.set_bind_group(0, &self.gather_bind[ccur * 2 + n], &[]);
        pass.dispatch_workgroups(dx, dy, 1);
    }

    /// En place : moments des nœuds parents des pochoirs, sur l'état de la
    /// grille parente d'avant son pas (`pcur` : état avant le pas).
    fn pre_before(&self, pass: &mut wgpu::ComputePass, pcur: usize) {
        if self.nn > 0 {
            let (dx, dy) = disp(self.nn);
            pass.set_pipeline(&self.pre);
            pass.set_bind_group(0, &self.pre_bind[pcur * 2], &[]);
            pass.dispatch_workgroups(dx, dy, 1);
        }
    }

    /// Remplissage des fantômes à α = a / 2 entre les deux derniers tirages.
    fn fill_at(&self, pass: &mut wgpu::ComputePass, ccur: usize, a: usize) {
        let (dx, dy) = disp(self.ng);
        pass.set_pipeline(&self.fill);
        pass.set_bind_group(0, &self.fill_bind[(ccur * 2 + a) * 2 + self.gnew.get()], &[]);
        pass.dispatch_workgroups(dx, dy, 1);
    }
}

/// Grilles fines d'une veine 3D et leurs couplages.
pub(crate) struct Refinement {
    /// Grilles 1, 2… (2ᵏ fois plus fines que la veine).
    pub(crate) levels: Vec<Gpu3D>,
    pub(crate) links: Vec<Link3>,
    pub(crate) fresh: bool,
    /// Nœud 0 de chaque grille fine, en mailles de la veine.
    pub(crate) origin: Vec<[f64; 3]>,
    /// Grilles qui entourent toute la maquette (les suivantes n'en couvrent
    /// qu'une partie : bords d'attaque et de fuite).
    pub(crate) full: usize,
    /// Parent de chaque grille fine : 0 pour la veine, i + 1 pour `levels[i]`
    /// (les grilles partielles d'un même niveau sont sœurs).
    pub(crate) parent: Vec<usize>,
    /// Niveau de chaque grille fine : 2^niveau fois plus fine que la veine.
    pub(crate) depth: Vec<u32>,
    /// Grilles partielles : domaine propre (mailles de la veine, coins bas et
    /// haut), où leurs efforts comptent ; les domaines de sœurs ne se
    /// recouvrent pas. None pour les grilles entières.
    pub(crate) core: Vec<Option<[[f64; 3]; 2]>>,
}

impl Refinement {
    /// Rapport de raffinement de `levels[i]` par rapport à la veine.
    pub(crate) fn scale(&self, i: usize) -> f64 {
        (1u64 << self.depth[i]) as f64
    }
}

/// Emprise des grilles fines autour de la maquette, en longueurs de
/// référence de la veine : amont, aval, de part et d'autre.
pub const EXTENT3: [(f64, f64, f64); 3] = [(0.6, 1.5, 0.6), (0.35, 0.7, 0.4), (0.2, 0.4, 0.25)];
pub const MARGIN3: usize = MARGIN;

/// Un pas de la grille `k` (0 : veine, i + 1 : `levels[i]`), suivi de deux
/// pas de chacune de ses grilles filles ; `links[c − 1]` relie la grille c à
/// son parent.
fn encode(levels: &mut [&mut Gpu3D], links: &[Link3], children: &[Vec<usize>], pass: &mut wgpu::ComputePass, k: usize, slot: u32) {
    {
        let lev = &mut *levels[k];
        if let Some(w) = &lev.wall_pass {
            w.encode(pass, slot);
        }
        if lev.ep {
            for &c in &children[k] {
                links[c - 1].pre_before(pass, lev.cur);
            }
        }
        lev.encode_step(pass, slot);
        if lev.sa {
            pass.set_pipeline(&lev.sa_pipe);
            pass.set_bind_group(0, &lev.sa_bind[lev.cur], &[slot]);
            pass.dispatch_workgroups(lev.sa_disp.0, lev.sa_disp.1, 1);
        }
        lev.cur = 1 - lev.cur;
    }
    for &c in &children[k] {
        let link = &links[c - 1];
        // Tirage, avec le remplissage du premier demi-pas (α = 0).
        link.gather(pass, levels[k].cur, levels[c].cur, levels[k].ep);
        let (rx, ry) = disp(link.nr);
        for a in 0..2 {
            if a == 1 {
                link.fill_at(pass, levels[c].cur, a);
            }
            encode(levels, links, children, pass, c, slot);
            // Recopie : en place, les mailles fluides après le premier pas fin
            // (leurs populations avant collision au temps t + 1 sont alors en
            // mémoire, le second pas les écrase), les autres après le second ;
            // à deux jeux, tout après le second.
            if a == 0 && !levels[k].ep {
                continue;
            }
            pass.set_pipeline(&link.restrict);
            pass.set_bind_group(0, &link.restrict_bind[(levels[k].cur * 2 + levels[c].cur) * 2 + a], &[]);
            pass.dispatch_workgroups(rx, ry, 1);
        }
    }
}

/// `n` pas de la veine `top` et de ses grilles fines.
pub(crate) fn run(top: &mut Gpu3D, r: &mut Refinement, n: usize, ramp: &crate::cpu::Ramp) {
    for g in r.levels.iter_mut() {
        g.ensure_spec();
    }
    for k in 0..r.links.len() {
        let p = if r.parent[k] == 0 { &*top } else { &r.levels[r.parent[k] - 1] };
        let link = &mut r.links[k];
        link.prepare(p, &r.levels[k]);
    }
    // Filles de chaque grille (indices de `all`).
    let mut children = vec![Vec::new(); r.levels.len() + 1];
    for k in 0..r.levels.len() {
        children[r.parent[k]].push(k + 1);
    }
    // Niveau de chaque grille de `all`.
    let lev_depth: Vec<u32> = std::iter::once(0).chain(r.depth.iter().copied()).collect();
    let device = top.device.clone();
    let queue = top.queue.clone();
    let mut all: Vec<&mut Gpu3D> = std::iter::once(top).chain(r.levels.iter_mut()).collect();
    if r.fresh {
        let mut enc = device.create_command_encoder(&Default::default());
        {
            let mut pass = enc.begin_compute_pass(&Default::default());
            for (k, link) in r.links.iter().enumerate() {
                let p = &all[r.parent[k]];
                for _ in 0..2 {
                    if p.ep {
                        link.pre_before(&mut pass, p.cur);
                    }
                    link.gather(&mut pass, p.cur, all[k + 1].cur, p.ep);
                }
            }
        }
        queue.submit([enc.finish()]);
        r.fresh = false;
    }
    let cost: usize = all.iter().enumerate().map(|(k, g)| g.n << lev_depth[k]).sum();
    let per_submit = (30_000_000 / cost.max(1)).clamp(1, crate::gpu3d::BATCH);
    let mut left = n;
    while left > 0 {
        let kc = left.min(per_submit);
        let t0 = all[0].time;
        for (k, g) in all.iter().enumerate() {
            let mut bytes = vec![0u8; crate::gpu3d::SLOT as usize * kc];
            for s in 0..kc {
                let (u, a) = ramp.at(t0 + s as u64);
                let last = left == kc && s == kc - 1;
                let p = g.params_for(u, a / (1u64 << lev_depth[k]) as f32, last);
                let o = s * crate::gpu3d::SLOT as usize;
                bytes[o..o + std::mem::size_of_val(&p)].copy_from_slice(bytemuck::bytes_of(&p));
            }
            queue.write_buffer(&g.params, 0, &bytes);
        }
        let mut enc = device.create_command_encoder(&Default::default());
        {
            let mut pass = enc.begin_compute_pass(&Default::default());
            for s in 0..kc {
                encode(&mut all, &r.links, &children, &mut pass, 0, (s as u64 * crate::gpu3d::SLOT) as u32);
            }
        }
        queue.submit([enc.finish()]);
        for (k, g) in all.iter_mut().enumerate() {
            let steps = (kc as u64) << lev_depth[k];
            g.time += steps;
            g.force_steps += steps;
        }
        left -= kc;
    }
}
