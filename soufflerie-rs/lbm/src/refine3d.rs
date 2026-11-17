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
};
@group(0) @binding(0) var<uniform> P: L;
const CV1_3 = 357.911;
fn nu_turb(v: f32, nu: f32) -> f32 {
  if (v <= 0.0) { return 0.0; }
  let chi = v / nu;
  let c3 = chi * chi * chi;
  return v * c3 / (c3 + CV1_3);
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
fn access(name: &str, ba: u32, bb: u32, fp16: bool, write: bool, n: &str) -> String {
    let mode = if write { "read_write" } else { "read" };
    let mut s = format!(
        "@group(0) @binding({ba}) var<storage, {mode}> {name}A: array<u32>;
@group(0) @binding({bb}) var<storage, {mode}> {name}B: array<u32>;
"
    );
    if fp16 {
        s += &format!(
            "fn ld_{name}(i: u32, c: u32) -> f32 {{
  let w = i >> 1u;
  var v: vec2<f32>;
  if (w < 5u) {{ v = unpack2x16float({name}A[w * {n} + c]); }} else {{ v = unpack2x16float({name}B[(w - 5u) * {n} + c]); }}
  return select(v.x, v.y, (i & 1u) == 1u) + W[i];
}}
"
        );
        if write {
            s += &format!("fn st_{name}(c: u32, p: array<f32, 19>) {{\n");
            for wd in 0..10 {
                let (a, b) = (2 * wd, 2 * wd + 1);
                let second = if b < 19 { format!("p[{b}] - W[{b}]") } else { "0.0".into() };
                let (buf, k) = if wd < 5 { ("A", wd) } else { ("B", wd - 5) };
                s += &format!("  {name}{buf}[{k}u * {n} + c] = pack2x16float(vec2<f32>(p[{a}] - W[{a}], {second}));\n");
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

/// Tirage, première passe : moments avant collision de chaque nœud parent
/// des pochoirs (liste `pnode`, chaque nœud une seule fois ; un nœud sert
/// à jusqu'à 64 mailles fantômes).
const PRE: &str = "
@group(0) @binding(1) var<storage, read_write> pm: array<vec4<f32>>;
@group(0) @binding(2) var<storage, read> pnode: array<i32>;
@group(0) @binding(5) var<storage, read> pnut: array<f32>;
@group(0) @binding(6) var<storage, read> pwall: array<u32>;
@group(0) @binding(8) var<storage, read> pblk: array<u32>;
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
fn pre_at(x: i32, y: i32, z: i32) -> Mom {
  var f: array<f32, 19>;
  let here = pcell(x, y, z);
  for (var i = 0u; i < 19u; i++) {
    let src = pcell(x - CX[i], y - CY[i], z - CZ[i]);
    if (pwall[src] == 0xffffffffu) {
      f[i] = ld_pf(OP[i], here);
    } else {
      f[i] = ld_pf(i, src);
    }
  }
  var m = moments(f);
  var v = 3.0 * P.nu_p;
  var vt = 0.0;
  if (P.sa != 0u) {
    v = pnut[pcell(x, y, z)];
    vt = nu_turb(v, P.nu_p);
  }
  let tau_p = 0.5 + 3.0 * (P.nu_p + vt);
  let tau_c = 0.5 + 3.0 * (P.nu_c + 2.0 * vt);
  let k = tau_c / (2.0 * tau_p);
  m.b = vec4<f32>(m.b.xyz * k, 2.0 * v);
  m.c = vec4<f32>(m.c.xyz * k, tau_c);
  return m;
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let k = flat(gid, nwg);
  if (k >= arrayLength(&pnode) / 3u) { return; }
  let m = pre_at(pnode[3u * k], pnode[3u * k + 1u], pnode[3u * k + 2u]);
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
@group(0) @binding(1) var<storage, read_write> gm: array<vec4<f32>>;
@group(0) @binding(2) var<storage, read> pm: array<vec4<f32>>;
@group(0) @binding(7) var<storage, read> gpos: array<u32>;
@group(0) @binding(9) var<storage, read> gstart: array<u32>;
@group(0) @binding(10) var<storage, read> gslot: array<u32>;

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

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let g = flat(gid, nwg);
  if (g >= P.ng) { return; }
  gm[6u * g] = gm[6u * g + 3u];
  gm[6u * g + 1u] = gm[6u * g + 4u];
  gm[6u * g + 2u] = gm[6u * g + 5u];
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
          let s = 3u * gslot[ptr];
          ptr++;
          m.a += w * pm[s];
          m.b += w * pm[s + 1u];
          m.c += w * pm[s + 2u];
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
          let sl = gslot[ptr];
          ptr++;
          if (sl == 0xffffffffu) { continue; }
          let s = 3u * sl;
          let w = lx.w[i] * ly.w[j] * lz.w[k];
          m.a += w * pm[s];
          m.b += w * pm[s + 1u];
          m.c += w * pm[s + 2u];
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
  gm[6u * g + 3u] = m.a;
  gm[6u * g + 4u] = m.b;
  gm[6u * g + 5u] = m.c;
}
";

const FILL: &str = "
@group(0) @binding(1) var<storage, read> gm: array<vec4<f32>>;
@group(0) @binding(2) var<storage, read> gidx: array<u32>;
@group(0) @binding(5) var<storage, read_write> cnut: array<f32>;
@group(0) @binding(6) var<storage, read_write> cmacro: array<vec4<f32>>;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let g = flat(gid, nwg);
  if (g >= P.ng) { return; }
  let a = mix(gm[6u * g], gm[6u * g + 3u], P.alpha);
  let b = mix(gm[6u * g + 1u], gm[6u * g + 4u], P.alpha);
  let c = mix(gm[6u * g + 2u], gm[6u * g + 5u], P.alpha);
  let k = 1.0 - 1.0 / c.w;
  let p = regul(a.x, a.yzw, k * b.x, k * b.y, k * b.z, k * c.x, k * c.y, k * c.z);
  let idx = gidx[g];
  st_cf(idx, p);
  if (P.sa != 0u) { cnut[idx] = b.w; }
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
  var v = 0.0;
  var vt = 0.0;
  if (P.sa != 0u) {
    v = cnut[cidx];
    vt = nu_turb(v, P.nu_c);
  }
  var f: array<f32, 19>;
  var p: array<f32, 19>;
  var m: Mom;
  if (cwall[cidx] == 0u) {
    // Populations avant collision au temps t + 1, tirées du demi-pas t + ½.
    for (var i = 0u; i < 19u; i++) {
      let s = vec3<i32>(i32(cx) - CX[i], i32(cy) - CY[i], i32(cz) - CZ[i]);
      f[i] = ld_cs(i, cid(u32(s.x), u32(s.y), u32(s.z)));
    }
    m = moments(f);
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
      f[i] = ld_cc(i, cidx);
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
      f[i] = ld_cc(i, cidx);
    }
    m = moments(f);
    p = regul(m.a.x, m.a.yzw, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  }
  st_pf(pidx, p);
  if (P.sa != 0u) { pnut[pidx] = 0.5 * v; }
  pmacro[pidx] = vec4<f32>(m.a.yzw, m.a.x);
}
";

/// Couplage d'une grille (parente) et de la suivante (fine).
pub(crate) struct Link3 {
    ox: usize,
    oy: usize,
    oz: usize,
    ng: u32,
    /// Nœuds parents recopiés depuis la grille fine.
    nr: u32,
    rlist: wgpu::Buffer,
    params: [wgpu::Buffer; 2],
    gm: wgpu::Buffer,
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
    pre: wgpu::ComputePipeline,
    pre_bind: Vec<wgpu::BindGroup>,
    gather: wgpu::ComputePipeline,
    fill: wgpu::ComputePipeline,
    restrict: wgpu::ComputePipeline,
    gather_bind: Vec<wgpu::BindGroup>,
    fill_bind: Vec<wgpu::BindGroup>,
    restrict_bind: Vec<wgpu::BindGroup>,
}

fn disp(n: u32) -> (u32, u32) {
    let groups = n.div_ceil(64).max(1);
    let dx = groups.min(65535);
    (dx, groups.div_ceil(dx))
}

impl Link3 {
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
                            pnode.extend_from_slice(&p);
                            (pnode.len() / 3 - 1) as u32
                        });
                        gslot.push(k);
                    }
                }
            }
        }
        let nn = (pnode.len() / 3) as u32;
        if pnode.is_empty() {
            pnode.extend_from_slice(&[0, 0, 0]);
        }
        if gstart.is_empty() {
            gstart.push(0);
        }
        if gslot.is_empty() {
            gslot.push(0);
        }
        let base = common();
        let module = |label: &str, body: String| {
            crate::gpu3d::compute_module(d, label, format!("{base}{body}"))
        };
        let pre = module("pre3d", format!("{}{PRE}", access("pf", 3, 4, fp16, false, "P.pn")));
        let gather = module("gather3d", GATHER.to_string());
        let fill = module("fill3d", format!("{}{FILL}", access("cf", 3, 4, fp16, true, "P.cn")));
        let restrict = module(
            "restrict3d",
            format!(
                "{}{}{}{RESTRICT}",
                access("cs", 1, 2, fp16, false, "P.cn"),
                access("cc", 3, 4, fp16, false, "P.cn"),
                access("pf", 7, 8, fp16, true, "P.pn")
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
        let nr = rl.len() as u32;
        if rl.is_empty() {
            rl.push(0);
        }
        Link3 {
            ox: o[0],
            oy: o[1],
            oz: o[2],
            ng,
            nr,
            rlist: d.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("restrict_list3d"),
                contents: bytemuck::cast_slice(&rl),
                usage: wgpu::BufferUsages::STORAGE,
            }),
            params: [
                mk("link3_0", std::mem::size_of::<LinkParams>() as u64, uni),
                mk("link3_1", std::mem::size_of::<LinkParams>() as u64, uni),
            ],
            gm: mk("ghost3d", 96 * ng as u64, wgpu::BufferUsages::STORAGE),
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
        for (a, buf) in self.params.iter().enumerate() {
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
                sa: (p.sa && c.sa) as u32,
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
        self.pre_bind = (0..2)
            .map(|cur| {
                group(
                    &self.pre,
                    vec![
                        (0, &self.params[0]),
                        (1, &self.pm),
                        (2, &self.pnode),
                        (3, &p.pops[1 - cur][0]),
                        (4, &p.pops[1 - cur][1]),
                        (5, &p.nut[cur]),
                        (6, &p.wall),
                        (8, &p.blk_buf),
                    ],
                )
            })
            .collect();
        self.gather_bind = vec![group(
            &self.gather,
            vec![
                (0, &self.params[0]),
                (1, &self.gm),
                (2, &self.pm),
                (7, &self.gpos),
                (9, &self.gstart),
                (10, &self.gslot),
            ],
        )];
        self.fill_bind = (0..4)
            .map(|i| {
                let (cur, a) = (i / 2, i % 2);
                group(
                    &self.fill,
                    vec![
                        (0, &self.params[a]),
                        (1, &self.gm),
                        (2, &self.gidx),
                        (3, &c.pops[cur][0]),
                        (4, &c.pops[cur][1]),
                        (5, &c.nut[cur]),
                        (6, &c.macro_buf),
                    ],
                )
            })
            .collect();
        self.restrict_bind = (0..4)
            .map(|i| {
                let (pc, cc) = (i / 2, i % 2);
                group(
                    &self.restrict,
                    vec![
                        (0, &self.params[0]),
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
                    ],
                )
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
            ("tirage".into(), self.ng as usize, time(p, &|pass| self.gather(pass, p.cur))),
            (
                "remplissage".into(),
                self.ng as usize,
                time(p, &|pass| {
                    pass.set_pipeline(&self.fill);
                    pass.set_bind_group(0, &self.fill_bind[c.cur * 2], &[]);
                    pass.dispatch_workgroups(dx, dy, 1);
                }),
            ),
            (
                "recopie".into(),
                self.nr as usize,
                time(p, &|pass| {
                    pass.set_pipeline(&self.restrict);
                    pass.set_bind_group(0, &self.restrict_bind[p.cur * 2 + c.cur], &[]);
                    pass.dispatch_workgroups(rx, ry, 1);
                }),
            ),
        ]
    }

    fn gather(&self, pass: &mut wgpu::ComputePass, pcur: usize) {
        if self.nn > 0 {
            let (dx, dy) = disp(self.nn);
            pass.set_pipeline(&self.pre);
            pass.set_bind_group(0, &self.pre_bind[pcur], &[]);
            pass.dispatch_workgroups(dx, dy, 1);
        }
        let (dx, dy) = disp(self.ng);
        pass.set_pipeline(&self.gather);
        pass.set_bind_group(0, &self.gather_bind[0], &[]);
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
pub(crate) const EXTENT3: [(f64, f64, f64); 3] = [(0.6, 1.5, 0.6), (0.35, 0.7, 0.4), (0.2, 0.4, 0.25)];
pub(crate) const MARGIN3: usize = MARGIN;

/// Un pas de la grille `k` (0 : veine, i + 1 : `levels[i]`), suivi de deux
/// pas de chacune de ses grilles filles ; `links[c − 1]` relie la grille c à
/// son parent.
fn encode(levels: &mut [&mut Gpu3D], links: &[Link3], children: &[Vec<usize>], pass: &mut wgpu::ComputePass, k: usize, slot: u32) {
    {
        let lev = &mut *levels[k];
        if let Some(w) = &lev.wall_pass {
            w.encode(pass, slot);
        }
        lev.encode_step(pass, slot);
        if lev.sa {
            pass.set_pipeline(&lev.sa_pipe);
            pass.set_bind_group(0, &lev.sa_bind[lev.cur], &[slot]);
            pass.dispatch_workgroups(lev.disp.0, lev.disp.1, 1);
        }
        lev.cur = 1 - lev.cur;
    }
    for &c in &children[k] {
        let link = &links[c - 1];
        link.gather(pass, levels[k].cur);
        let (dx, dy) = disp(link.ng);
        for a in 0..2 {
            pass.set_pipeline(&link.fill);
            pass.set_bind_group(0, &link.fill_bind[levels[c].cur * 2 + a], &[]);
            pass.dispatch_workgroups(dx, dy, 1);
            encode(levels, links, children, pass, c, slot);
        }
        let nr = link.nr;
        let (rx, ry) = disp(nr);
        pass.set_pipeline(&link.restrict);
        pass.set_bind_group(0, &link.restrict_bind[levels[k].cur * 2 + levels[c].cur], &[]);
        pass.dispatch_workgroups(rx, ry, 1);
    }
}

/// `n` pas de la veine `top` et de ses grilles fines.
pub(crate) fn run(top: &mut Gpu3D, r: &mut Refinement, n: usize, ramp: &crate::cpu::Ramp) {
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
                for _ in 0..2 {
                    link.gather(&mut pass, all[r.parent[k]].cur);
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
