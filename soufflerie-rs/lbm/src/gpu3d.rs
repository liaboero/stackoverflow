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

fn w3(i: usize) -> f64 {
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

const SOLID: u32 = 0xffff_ffff;

/// Littéral flottant WGSL.
fn f(v: f64) -> String {
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
  rows: u32, pad0: u32, pad1: u32, pad2: u32,
};
";

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
      if (q < 0.5 && qo == 0.0) {{
        let nb = cell(x + {cx}, y + {cy}, z + {cz});
        v = 2.0 * q * fj + (1.0 - 2.0 * q) * pop({j}u, nb);
      }} else if (q < 0.5) {{
        v = fj;
      }} else {{
        v = (fj + (2.0 * q - 1.0) * pop({i}u, idx)) / (2.0 * q);
      }}
      let mom = fj + v + {w2};
      fsum += vec3<f32>({mcx}, {mcy}, {mcz}) * mom;
    }} else if (x - {cx} < 0) {{
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
            w2 = f(2.0 * w3(i)),
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
        let _ = write!(
            coll,
            "  {{
    let cu = {cu};
    let feq = {wf} * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - usq);
    let fnew = feq + k * {w45} * ({qp}) + {w3f} * rho * ({cxs} * P.accel);
    o[{i}] = fnew - {wf};
  }}
",
            wf = f(w),
            w45 = f(4.5 * w),
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
  if (wk0 == {SOLID}u) {{ return; }}
  let wk = select(wk0 - 1u, 0xffffffffu, wk0 == 0u);
  let hasWall = wk != 0xffffffffu;
  let x = i32(idx % P.nx);
  let y = i32((idx / P.nx) % P.ny);
  let z = i32(idx / (P.nx * P.ny));
  var g: array<f32, 19>;
  var fsum = vec3<f32>(0.0);
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
  if (u32(x) > P.spongeStart) {{
    let s = f32(u32(x) - P.spongeStart) / f32(max(1u, P.nx - 1u - P.spongeStart));
    tau += (P.tauSponge - P.tau0) * s * s;
  }}
  var om = 1.0 / tau;
  if (P.les == 1u) {{
    let qn = sqrt(2.0 * (pxx * pxx + pyy * pyy + pzz * pzz + 2.0 * (pxy * pxy + pxz * pxz + pyz * pyz)));
    om = 2.0 / (tau + sqrt(tau * tau + P.c18 * qn * ir));
  }}
  let k = 1.0 - om;
  var o: array<f32, 19>;
{coll}
{out}
  if (hasWall) {{
    force[wk] += vec4<f32>(fsum, 0.0);
  }}
  if (P.writeMacro == 1u) {{
    macroBuf[idx] = vec4<f32>(u, rho);
  }}
}}
",
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
};
@group(0) @binding(0) var<uniform> V: R;
@group(0) @binding(1) var<storage, read> macroBuf: array<vec4<f32>>;
@group(0) @binding(2) var<storage, read> wall: array<u32>;
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
fn solidI(p: vec3<i32>) -> f32 { return select(0.0, 1.0, wall[id(p)] == 0xffffffffu); }

// Interpolation trilinéaire (centres de maille aux coordonnées entières).
fn occ(p: vec3<f32>) -> f32 {
  let b = vec3<i32>(floor(p));
  let f = p - floor(p);
  let c00 = mix(solidI(b), solidI(b + vec3<i32>(1, 0, 0)), f.x);
  let c10 = mix(solidI(b + vec3<i32>(0, 1, 0)), solidI(b + vec3<i32>(1, 1, 0)), f.x);
  let c01 = mix(solidI(b + vec3<i32>(0, 0, 1)), solidI(b + vec3<i32>(1, 0, 1)), f.x);
  let c11 = mix(solidI(b + vec3<i32>(0, 1, 1)), solidI(b + vec3<i32>(1, 1, 1)), f.x);
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
    for (var k = 0; k < 4000; k++) {
      if (t > tend) { break; }
      let p = o + dir * t;
      if (isoOn && ti > 1e29 && Q(p) > level) { ti = t; tip = tp; }
      if (solidI(vec3<i32>(round(p))) > 0.5) {
        // Surface lissée : bisection sur l'occupation trilinéaire.
        var a = tp;
        var b = t;
        for (var j = 0; j < 6; j++) {
          let m = 0.5 * (a + b);
          if (occ(o + dir * m) > 0.5) { b = m; } else { a = m; }
        }
        tb = b;
        break;
      }
      tp = t;
      t += stp;
    }
    if (tb < 1e29) {
      let p = o + dir * tb;
      let e = 1.0;
      let g = vec3<f32>(
        occ(p + vec3<f32>(e, 0.0, 0.0)) - occ(p - vec3<f32>(e, 0.0, 0.0)),
        occ(p + vec3<f32>(0.0, e, 0.0)) - occ(p - vec3<f32>(0.0, e, 0.0)),
        occ(p + vec3<f32>(0.0, 0.0, e)) - occ(p - vec3<f32>(0.0, 0.0, e)));
      let n = -normalize(g + vec3<f32>(1e-9));
      let lam = abs(dot(n, -dir));
      let key = max(0.0, dot(n, normalize(vec3<f32>(-0.4, -0.3, 0.85))));
      col = V.body.rgb * (0.35 + 0.45 * lam + 0.25 * key);
    } else if (ts < 1e29) {
      let p = o + dir * ts;
      if (solidI(vec3<i32>(round(p))) > 0.5) {
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
}

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct Params {
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
    pad: [u32; 4],
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
const BATCH: usize = 64;
const SLOT: u64 = 256;

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
    device: wgpu::Device,
    queue: wgpu::Queue,
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
    _pops: [[wgpu::Buffer; 2]; 2],
    wall: wgpu::Buffer,
    _qbuf: wgpu::Buffer,
    force: wgpu::Buffer,
    force_read: wgpu::Buffer,
    params: wgpu::Buffer,
    macro_buf: wgpu::Buffer,
    step_pipe: wgpu::ComputePipeline,
    init_pipe: wgpu::ComputePipeline,
    slice_pipe: wgpu::ComputePipeline,
    derive_pipe: wgpu::ComputePipeline,
    render_pipe: wgpu::ComputePipeline,
    /// ω et Q par maille (alloué au premier rendu 3D).
    deriv: Option<wgpu::Buffer>,
    /// Tampons du rendu réutilisés d'une image à l'autre : taille, image, relecture.
    render_bufs: Option<(usize, usize, wgpu::Buffer, wgpu::Buffer)>,
    bind: Vec<wgpu::BindGroup>,
    init_bind: Vec<wgpu::BindGroup>,
    disp: (u32, u32),
    cur: usize,
    pub time: u64,
    force_steps: u64,
    pub wall_count: usize,
    pub solid_count: usize,
    /// Mailles solides (pour l'affichage).
    solid: Vec<bool>,
    /// Le GPU ne répond plus (lecture impossible, pilote réinitialisé).
    pub lost: bool,
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
        let (nx, ny, nz) = (body.nx, body.ny, body.nz);
        let n = nx * ny * nz;
        let al = adapter.limits();
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
        let mut limits = wgpu::Limits::default();
        limits.max_storage_buffer_binding_size = al.max_storage_buffer_binding_size;
        limits.max_buffer_size = al.max_buffer_size;
        limits.max_storage_buffers_per_shader_stage = limits.max_storage_buffers_per_shader_stage.max(8);
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            label: Some("soufflerie 3D"),
            required_limits: limits,
            ..Default::default()
        }))
        .map_err(|e| format!("GPU indisponible : {e}"))?;
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
        let mut k = 0u32;
        for z in z0.saturating_sub(1)..=(z1 + 1).min(nz - 1) {
            for y in y0.saturating_sub(1)..=(y1 + 1).min(ny - 1) {
                for x in (x0 - 1).max(1)..=(x1 + 1).min(nx - 2) {
                    let c = idx(x, y, z);
                    if wall[c] == SOLID {
                        continue;
                    }
                    let mut row = [0f32; 18];
                    let mut any = false;
                    for i in 1..19 {
                        let [cx, cy, cz] = C3[i];
                        let sx = (x as i32 - cx) as usize;
                        let sy = ((y as i32 - cy + ny as i32) % ny as i32) as usize;
                        let sz = ((z as i32 - cz + nz as i32) % nz as i32) as usize;
                        if wall[idx(sx, sy, sz)] == SOLID {
                            row[i - 1] = link_q(
                                body,
                                [x as f64, y as f64, z as f64],
                                [-cx as f64, -cy as f64, -cz as f64],
                            );
                            any = true;
                        }
                    }
                    if any {
                        k += 1;
                        wall[c] = k;
                        qs.extend_from_slice(&row);
                    }
                }
            }
        }
        let wall_count = k as usize;
        if qs.is_empty() {
            qs.resize(18, 0.0);
        }
        let solid: Vec<bool> = wall.iter().map(|&w| w == SOLID).collect();
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
            sponge_start: (nx as f64 * 0.85).round() as u32,
            cs: 0.15,
            _pops: pops,
            wall: wall_buf,
            _qbuf: qbuf,
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
        };
        s.set_viscosity(nu);
        s.reset(0.0);
        Ok(s)
    }

    pub fn adapter_name(&self) -> &str {
        &self.adapter_name
    }

    pub fn set_viscosity(&mut self, nu: f64) {
        self.tau0 = 0.5 + 3.0 * nu;
        self.tau_sponge = self.tau0.max(0.8);
        // Sous-maille seulement aux Reynolds de maille élevés : sans amortissement
        // pariétal, elle épaissirait la couche limite et gonflerait la traînée.
        self.les = self.tau0 < 0.51;
    }

    pub fn is_solid(&self, x: usize, y: usize, z: usize) -> bool {
        self.solid[x + self.nx * (y + self.ny * z)]
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
            pad: [0; 4],
        }
    }

    /// Fluide au repos (ou à la vitesse u0).
    pub fn reset(&mut self, u0: f32) {
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
        self.cur = 0;
        self.time = 0;
        self.force_steps = 0;
    }

    /// `n` pas ; la vitesse d'entrée suit `ramp`. Le champ macroscopique est
    /// écrit au dernier pas (pour les coupes).
    pub fn run(&mut self, n: usize, ramp: &Ramp) {
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
                    pass.set_bind_group(0, &self.bind[self.cur], &[(s as u64 * SLOT) as u32]);
                    pass.dispatch_workgroups(self.disp.0, self.disp.1, 1);
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
    pub fn take_forces(&mut self) -> Force3 {
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
            size: [w as f32, h as f32, 0.0, 0.0],
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
                buf_entry(2, &self.wall, None),
                buf_entry(3, deriv, None),
                buf_entry(4, img, None),
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
