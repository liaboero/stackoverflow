// Spalart–Allmaras : un pas explicite du transport de ν̃ (même schéma que
// sa.rs). Lit `sa_in` = (u_x, u_y, u_τ, ν̃) du pas LBM qui vient de finir et
// écrit `sa_out` (vitesses et u_τ recopiées, nouveau ν̃).

struct Params {
  nx: u32, ny: u32, n: u32, flags: u32,
  u_in: f32, accel: f32, c18: f32, spin: f32,
  rot_x: f32, rot_y: f32, ref_x: f32, ref_y: f32,
  nu: f32, pad0: f32, pad1: f32, pad2: f32,
};

@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read> kind: array<u32>;
@group(0) @binding(2) var<storage, read> aux: array<f32>;
@group(0) @binding(3) var<storage, read> cell: array<vec4<f32>>;
@group(0) @binding(4) var<storage, read> sa_in: array<vec4<f32>>;
@group(0) @binding(5) var<storage, read_write> sa_out: array<vec4<f32>>;

const CB1 = 0.1355;
const CB2 = 0.622;
const SIGMA = 0.6666667;
const KAPPA = 0.41;
const CW2 = 0.3;
const CW3_6 = 64.0;
const CV1_3 = 357.911;
const NU_TILDE_MAX = 0.15;

fn nt(i: u32) -> f32 {
  return sa_in[i].w;
}

@compute @workgroup_size(64, 1, 1)
fn sa_step(@builtin(global_invocation_id) gid: vec3<u32>) {
  let x = gid.x;
  let y = gid.y;
  if (x >= P.nx || y >= P.ny) {
    return;
  }
  let nx = P.nx;
  let idx = y * nx + x;
  let me = sa_in[idx];
  let k = kind[idx];
  let nu = P.nu;
  let inf = 3.0 * nu;
  var out: f32;
  if ((k & 3u) == 1u) {
    out = 0.0;
  } else if (x == 0u) {
    out = inf;
  } else if (x == nx - 1u) {
    out = nt(idx - 1u);
  } else if (y == 0u) {
    out = nt(idx + nx);
  } else if (y == P.ny - 1u) {
    out = nt(idx - nx);
  } else if ((k >> 2u) != 0u) {
    // Premier nœud : ν̃ ≈ ν_t = κ u_τ y (zone logarithmique).
    out = min(NU_TILDE_MAX, KAPPA * me.z * max(0.1, cell[idx].y));
  } else {
    let cg = cell[idx];
    let v = me.w;
    let vl = nt(idx - 1u);
    let vr = nt(idx + 1u);
    let vd = nt(idx - nx);
    let vu = nt(idx + nx);
    let a = me.x;
    let b = me.y;
    var adv = 0.0;
    if (a > 0.0) {
      adv = a * (v - vl);
    } else {
      adv = a * (vr - v);
    }
    if (b > 0.0) {
      adv += b * (v - vd);
    } else {
      adv += b * (vu - v);
    }
    var diff = (nu + 0.5 * (v + vr)) * (vr - v) - (nu + 0.5 * (v + vl)) * (v - vl)
      + (nu + 0.5 * (v + vu)) * (vu - v) - (nu + 0.5 * (v + vd)) * (v - vd);
    let gx = 0.5 * (vr - vl);
    let gy = 0.5 * (vu - vd);
    diff = (diff + CB2 * (gx * gx + gy * gy)) / SIGMA;
    let omega = abs(0.5 * (sa_in[idx + 1u].y - sa_in[idx - 1u].y) - 0.5 * (sa_in[idx + nx].x - sa_in[idx - nx].x));
    let d = max(0.5, cg.y);
    let kd2 = KAPPA * KAPPA * d * d;
    let chi = v / nu;
    let c3 = chi * chi * chi;
    let fv1 = c3 / (c3 + CV1_3);
    let fv2 = 1.0 - chi / (1.0 + chi * fv1);
    let s_tilde = max(omega + v / kd2 * fv2, 0.3 * omega);
    var prod = CB1 * s_tilde * v;
    if (cg.z > 0.0) {
      prod = 0.0;
    }
    var r = 10.0;
    if (s_tilde > 1e-12) {
      r = min(10.0, v / (s_tilde * kd2));
    }
    let r2 = r * r;
    let g = r + CW2 * (r2 * r2 * r2 - r);
    let g2 = g * g;
    let fw = g * pow((1.0 + CW3_6) / (g2 * g2 * g2 + CW3_6), 1.0 / 6.0);
    let cw1 = CB1 / (KAPPA * KAPPA) + (1.0 + CB2) / SIGMA;
    var nv = (v + prod + diff - adv) / (1.0 + cw1 * fw * v / (d * d));
    let sg = aux[nx + y];
    if (sg > 0.0) {
      nv -= sg * (nv - inf);
    }
    out = clamp(nv, 0.0, NU_TILDE_MAX);
  }
  sa_out[idx] = vec4<f32>(me.xyz, out);
}
