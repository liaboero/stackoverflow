// D2Q9 : propagation « pull » et collision fusionnées (même schéma que cpu.rs).
// Disposition : f[i * n + y * nx + x]. Type de maille dans les 2 bits de poids
// faible de `kind` (0 fluide, 1 solide, 2 voisin de paroi), indice de maille
// pariétale + 1 au-dessus (efforts).

struct Params {
  nx: u32, ny: u32, n: u32, flags: u32,
  u_in: f32, accel: f32, c18: f32, spin: f32,
  rot_x: f32, rot_y: f32, ref_x: f32, ref_y: f32,
};

@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read> src: array<f32>;
@group(0) @binding(2) var<storage, read_write> dst: array<f32>;
@group(0) @binding(3) var<storage, read> kind: array<u32>;
@group(0) @binding(4) var<storage, read> wq: array<f32>;
@group(0) @binding(5) var<storage, read> aux: array<f32>;
@group(0) @binding(6) var<storage, read> damp: array<f32>;
@group(0) @binding(7) var<storage, read_write> force: array<vec4<f32>>;

const CX = array<i32, 9>(0, 1, 0, -1, 0, 1, -1, -1, 1);
const CY = array<i32, 9>(0, 0, 1, 0, -1, 1, 1, -1, -1);
const W = array<f32, 9>(0.44444445, 0.11111111, 0.11111111, 0.11111111, 0.11111111, 0.027777778, 0.027777778, 0.027777778, 0.027777778);
const OPP = array<u32, 9>(0u, 3u, 4u, 1u, 2u, 7u, 8u, 5u, 6u);
const REFY = array<u32, 9>(0u, 1u, 4u, 3u, 2u, 8u, 7u, 6u, 5u);

fn ld(i: u32, x: u32, y: u32) -> f32 {
  return src[i * P.n + y * P.nx + x];
}

fn feq(rho: f32, ux: f32, uy: f32) -> array<f32, 9> {
  var e: array<f32, 9>;
  let usq = 1.5 * (ux * ux + uy * uy);
  for (var i = 0u; i < 9u; i++) {
    let cu = f32(CX[i]) * ux + f32(CY[i]) * uy;
    e[i] = W[i] * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - usq);
  }
  return e;
}

fn collide(f: array<f32, 9>, tau: f32, c18: f32, sigma: f32) -> array<f32, 9> {
  let rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
  let ir = 1.0 / rho;
  let ux = (f[1] - f[3] + f[5] - f[6] - f[7] + f[8]) * ir;
  let uy = (f[2] - f[4] + f[5] + f[6] - f[7] - f[8]) * ir;
  let les = (P.flags & 2u) != 0u;
  var p: array<f32, 9>;
  if ((P.flags & 1u) != 0u) {
    let d4 = f[5] + f[6] + f[7] + f[8];
    let rux = rho * ux;
    let ruy = rho * uy;
    let a1xx = f[1] + f[3] + d4 - rho / 3.0 - rux * ux;
    let a1yy = f[2] + f[4] + d4 - rho / 3.0 - ruy * uy;
    let a1xy = f[5] - f[6] + f[7] - f[8] - rux * uy;
    var om = 1.0 / tau;
    if (les) {
      let qn = sqrt(2.0 * (a1xx * a1xx + a1yy * a1yy + 2.0 * a1xy * a1xy));
      om = 2.0 / (tau + sqrt(tau * tau + c18 * qn * ir));
    }
    let k = 1.0 - om;
    let bxx = k * a1xx;
    let byy = k * a1yy;
    let bxy = k * a1xy;
    let pxx = rux * ux + bxx;
    let pyy = ruy * uy + byy;
    let pxy = rux * uy + bxy;
    let qxxy = rux * ux * uy + 2.0 * ux * bxy + uy * bxx;
    let qxyy = rux * uy * uy + 2.0 * uy * bxy + ux * byy;
    let sp = pxx + pyy;
    let ax = 3.0 * rux;
    let ay = 3.0 * ruy;
    let tx = rho + 3.0 * pxx - 1.5 * pyy;
    let ty = rho + 3.0 * pyy - 1.5 * pxx;
    let td = rho + 3.0 * sp;
    let qs = 9.0 * (qxxy + qxyy);
    let qd = 9.0 * (qxxy - qxyy);
    let pxy9 = 9.0 * pxy;
    p[0] = W[0] * (rho - 1.5 * sp);
    p[1] = W[1] * (tx + ax - 4.5 * qxyy);
    p[2] = W[2] * (ty + ay - 4.5 * qxxy);
    p[3] = W[3] * (tx - ax + 4.5 * qxyy);
    p[4] = W[4] * (ty - ay + 4.5 * qxxy);
    p[5] = W[5] * (td + ax + ay + pxy9 + qs);
    p[6] = W[6] * (td - ax + ay - pxy9 + qd);
    p[7] = W[7] * (td - ax - ay + pxy9 - qs);
    p[8] = W[8] * (td + ax - ay - pxy9 - qd);
  } else {
    let e = feq(rho, ux, uy);
    var om = 1.0 / tau;
    if (les) {
      let pd = (f[5] - e[5]) + (f[6] - e[6]) + (f[7] - e[7]) + (f[8] - e[8]);
      let pxx = (f[1] - e[1]) + (f[3] - e[3]) + pd;
      let pyy = (f[2] - e[2]) + (f[4] - e[4]) + pd;
      let pxy = (f[5] - e[5]) - (f[6] - e[6]) + (f[7] - e[7]) - (f[8] - e[8]);
      let qn = sqrt(2.0 * (pxx * pxx + pyy * pyy + 2.0 * pxy * pxy));
      om = 2.0 / (tau + sqrt(tau * tau + c18 * qn * ir));
    }
    for (var i = 0u; i < 9u; i++) {
      p[i] = f[i] + om * (e[i] - f[i]);
    }
  }
  // Forçage uniforme selon x (mise en vitesse).
  let g1 = P.accel / 3.0;
  let g2 = P.accel / 12.0;
  p[1] += g1;
  p[3] -= g1;
  p[5] += g2;
  p[8] += g2;
  p[6] -= g2;
  p[7] -= g2;
  // Bandes absorbantes : rappel vers l'équilibre amont.
  if (sigma > 0.0) {
    let q = feq(1.0, P.u_in, 0.0);
    for (var i = 0u; i < 9u; i++) {
      p[i] -= sigma * (p[i] - q[i]);
    }
  }
  return p;
}

@compute @workgroup_size(64, 1, 1)
fn step(@builtin(global_invocation_id) gid: vec3<u32>) {
  let x = gid.x;
  let y = gid.y;
  if (x < 1u || x + 1u >= P.nx || y >= P.ny) {
    return;
  }
  let idx = y * P.nx + x;
  let k = kind[idx];
  let t = k & 3u;
  if (t == 1u) {
    return;
  }
  var f: array<f32, 9>;
  var fsum = vec3<f32>(0.0);
  if (t == 0u) {
    for (var i = 0u; i < 9u; i++) {
      f[i] = ld(i, u32(i32(x) - CX[i]), u32(i32(y) - CY[i]));
    }
  } else {
    f[0] = ld(0u, x, y);
    for (var j = 1u; j < 9u; j++) {
      let cx = CX[j];
      let cy = CY[j];
      let sy = i32(y) - cy;
      let sx = u32(i32(x) - cx);
      if (sy < 0 || sy >= i32(P.ny)) {
        // Parois haute et basse glissantes : réflexion spéculaire.
        f[j] = ld(REFY[j], sx, y);
        continue;
      }
      let q = wq[j * P.n + idx];
      if (q == 0.0) {
        f[j] = ld(j, sx, u32(sy));
        continue;
      }
      // Lien vers la paroi : rebond interpolé de Bouzidi, paroi mobile si rotation.
      let i = OPP[j];
      let fi = ld(i, x, y);
      let wx = f32(x) - q * f32(cx);
      let wy = f32(y) - q * f32(cy);
      var uw = vec2<f32>(0.0);
      if (P.spin != 0.0) {
        uw = vec2<f32>(-P.spin * (wy - P.rot_y), P.spin * (wx - P.rot_x));
      }
      let corr = 6.0 * W[j] * (f32(cx) * uw.x + f32(cy) * uw.y);
      var v: f32;
      if ((P.flags & 4u) == 0u) {
        v = fi + corr;
      } else if (q < 0.5) {
        let bx = i32(x) + cx;
        let by = i32(y) + cy;
        if (by >= 0 && by < i32(P.ny) && bx >= 0 && bx < i32(P.nx) && (kind[u32(by) * P.nx + u32(bx)] & 3u) != 1u) {
          v = 2.0 * q * fi + (1.0 - 2.0 * q) * ld(i, u32(bx), u32(by)) + corr;
        } else {
          v = fi + corr;
        }
      } else {
        v = (fi + (2.0 * q - 1.0) * ld(j, x, y) + corr) / (2.0 * q);
      }
      f[j] = v;
      let s = fi + v;
      let d = fi - v;
      let fxl = -f32(cx) * s - uw.x * d;
      let fyl = -f32(cy) * s - uw.y * d;
      fsum += vec3<f32>(fxl, fyl, (wx - P.ref_x) * fyl - (wy - P.ref_y) * fxl);
    }
  }
  var c18 = 0.0;
  if ((P.flags & 2u) != 0u) {
    c18 = P.c18 * damp[idx];
  }
  let p = collide(f, aux[x], c18, aux[P.nx + y]);
  for (var i = 0u; i < 9u; i++) {
    dst[i * P.n + idx] = p[i];
  }
  let wid = k >> 2u;
  if (wid != 0u) {
    force[wid - 1u] += vec4<f32>(fsum, 0.0);
  }
}

// Entrée (vitesse imposée, densité extrapolée) et sortie (ρ = 1, vitesse
// extrapolée) : partie hors équilibre copiée du voisin (Guo 2002).
@compute @workgroup_size(64, 1, 1)
fn boundary(@builtin(global_invocation_id) gid: vec3<u32>) {
  let y = gid.x;
  if (y >= P.ny) {
    return;
  }
  for (var side = 0u; side < 2u; side++) {
    var b = 0u;
    var nb = 1u;
    if (side == 1u) {
      b = P.nx - 1u;
      nb = P.nx - 2u;
    }
    var fn_: array<f32, 9>;
    for (var i = 0u; i < 9u; i++) {
      fn_[i] = dst[i * P.n + y * P.nx + nb];
    }
    let rho = fn_[0] + fn_[1] + fn_[2] + fn_[3] + fn_[4] + fn_[5] + fn_[6] + fn_[7] + fn_[8];
    let ux = (fn_[1] - fn_[3] + fn_[5] - fn_[6] - fn_[7] + fn_[8]) / rho;
    let uy = (fn_[2] - fn_[4] + fn_[5] + fn_[6] - fn_[7] - fn_[8]) / rho;
    var eb: array<f32, 9>;
    if (side == 0u) {
      eb = feq(rho, P.u_in, 0.0);
    } else {
      eb = feq(1.0, ux, uy);
    }
    let en = feq(rho, ux, uy);
    for (var i = 0u; i < 9u; i++) {
      dst[i * P.n + y * P.nx + b] = eb[i] + fn_[i] - en[i];
    }
  }
}
