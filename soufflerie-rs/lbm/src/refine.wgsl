// Raffinement local : couplage d'une grille fine (pas d'espace et de temps
// moitiés) avec sa grille parente, en moments (ρ, u, partie hors équilibre
// du tenseur Π), populations reconstruites par régularisation comme dans la
// collision (Lagrava et al. 2012, Dupuis et Chopard 2003).
//
// Le nœud fin (x, y) est au point parent (ox + x/2, oy + y/2). Le pourtour
// de la grille fine est fait de mailles fantômes (type 3) : `gather` relève
// après chaque pas parent les moments de la grille parente à leurs positions
// (interpolation cubique le long du bord), `fill` reconstruit avant chaque
// pas fin leurs populations post-collision, au temps t (alpha = 0) puis
// t + ½ (alpha = ½, moyenne de t et t + 1). Après les deux pas fins,
// `restrict` recopie la grille fine dans la grille parente sous la zone
// raffinée (à `margin` mailles parentes du bord).
//
// Mise à l'échelle (unités de réseau de chaque grille) : u identique,
// ν × 2, gradients ÷ 2 ; partie hors équilibre avant collision
// Π_fin = Π_parent · τ_fin / (2 τ_parent).

struct L {
  pnx: u32, pny: u32, pn: u32, cnx: u32,
  cny: u32, cn: u32, ox: u32, oy: u32,
  ng: u32, margin: u32, sa_on: u32, nrx: u32,
  alpha: f32, nu_p: f32, nu_c: f32, nry: u32,
};

@group(0) @binding(0) var<uniform> P: L;

const CX = array<i32, 9>(0, 1, 0, -1, 0, 1, -1, -1, 1);
const CY = array<i32, 9>(0, 0, 1, 0, -1, 1, 1, -1, -1);
const W = array<f32, 9>(0.44444445, 0.11111111, 0.11111111, 0.11111111, 0.11111111, 0.027777778, 0.027777778, 0.027777778, 0.027777778);
const CV1_3 = 357.911;

fn nu_turb(v: f32, nu: f32) -> f32 {
  if (v <= 0.0) {
    return 0.0;
  }
  let chi = v / nu;
  let c3 = chi * chi * chi;
  return v * c3 / (c3 + CV1_3);
}

// a = (ρ, u_x, u_y, ν̃), b = (Π_xx, Π_yy, Π_xy hors équilibre, τ).
struct Mom {
  a: vec4<f32>,
  b: vec4<f32>,
};

fn moments(f: array<f32, 9>) -> Mom {
  let rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
  let ux = (f[1] - f[3] + f[5] - f[6] - f[7] + f[8]) / rho;
  let uy = (f[2] - f[4] + f[5] + f[6] - f[7] - f[8]) / rho;
  let d4 = f[5] + f[6] + f[7] + f[8];
  let axx = f[1] + f[3] + d4 - rho / 3.0 - rho * ux * ux;
  let ayy = f[2] + f[4] + d4 - rho / 3.0 - rho * uy * uy;
  let axy = f[5] - f[6] + f[7] - f[8] - rho * ux * uy;
  return Mom(vec4<f32>(rho, ux, uy, 0.0), vec4<f32>(axx, ayy, axy, 0.0));
}

// Populations régularisées récursives (même reconstruction que la collision)
// pour ρ, u et la partie hors équilibre (bxx, byy, bxy).
fn regularized(rho: f32, ux: f32, uy: f32, bxx: f32, byy: f32, bxy: f32) -> array<f32, 9> {
  let rux = rho * ux;
  let ruy = rho * uy;
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
  var p: array<f32, 9>;
  p[0] = W[0] * (rho - 1.5 * sp);
  p[1] = W[1] * (tx + ax - 4.5 * qxyy);
  p[2] = W[2] * (ty + ay - 4.5 * qxxy);
  p[3] = W[3] * (tx - ax + 4.5 * qxyy);
  p[4] = W[4] * (ty - ay + 4.5 * qxxy);
  p[5] = W[5] * (td + ax + ay + pxy9 + qs);
  p[6] = W[6] * (td - ax + ay - pxy9 + qd);
  p[7] = W[7] * (td - ax - ay + pxy9 - qs);
  p[8] = W[8] * (td + ax - ay - pxy9 - qd);
  return p;
}

// Maille fantôme numéro g : lignes basse et haute, puis colonnes gauche et droite.
fn ghost_xy(g: u32) -> vec2<u32> {
  if (g < P.cnx) {
    return vec2<u32>(g, 0u);
  }
  if (g < 2u * P.cnx) {
    return vec2<u32>(g - P.cnx, P.cny - 1u);
  }
  let r = g - 2u * P.cnx;
  if (r < P.cny - 2u) {
    return vec2<u32>(0u, r + 1u);
  }
  return vec2<u32>(P.cnx - 1u, r - (P.cny - 2u) + 1u);
}
