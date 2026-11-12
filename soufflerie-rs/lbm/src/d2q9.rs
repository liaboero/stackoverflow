//! Réseau D2Q9 : vitesses, poids, équilibre.

/// Directions : 0 repos, 1 (+x), 2 (+y), 3 (−x), 4 (−y), 5 (+x+y), 6 (−x+y), 7 (−x−y), 8 (+x−y).
pub const CX: [i32; 9] = [0, 1, 0, -1, 0, 1, -1, -1, 1];
pub const CY: [i32; 9] = [0, 0, 1, 0, -1, 1, 1, -1, -1];
pub const W: [f32; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];
/// Direction opposée.
pub const OPP: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];
/// Réflexion spéculaire selon y (parois haute et basse glissantes).
pub const REFY: [usize; 9] = [0, 1, 4, 3, 2, 8, 7, 6, 5];

/// Équilibre de Maxwell–Boltzmann tronqué à l'ordre 2.
#[inline(always)]
pub fn equilibrium(rho: f32, ux: f32, uy: f32) -> [f32; 9] {
    let usq = 1.5 * (ux * ux + uy * uy);
    let mut out = [0.0; 9];
    for i in 0..9 {
        let cu = CX[i] as f32 * ux + CY[i] as f32 * uy;
        out[i] = W[i] * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - usq);
    }
    out
}

/// Masse volumique et vitesse.
#[inline(always)]
pub fn moments(f: &[f32; 9]) -> (f32, f32, f32) {
    let rho: f32 = f.iter().sum();
    let ux = (f[1] - f[3] + f[5] - f[6] - f[7] + f[8]) / rho;
    let uy = (f[2] - f[4] + f[5] + f[6] - f[7] - f[8]) / rho;
    (rho, ux, uy)
}

/// Paramètres locaux de collision.
#[derive(Clone, Copy)]
pub struct Relax {
    /// Temps de relaxation de base (éponge comprise).
    pub tau: f32,
    /// 18·C_s²·amortissement pariétal (0 : pas de sous-maille).
    pub c18: f32,
    /// Accélération uniforme selon x (mise en vitesse).
    pub gx: f32,
    /// Coefficient de rappel des bandes absorbantes.
    pub sigma: f32,
}

/// Collision, régularisée récursive (`RR`) ou BGK, avec sous-maille (`LES`),
/// forçage et rappel vers l'équilibre amont `q` des bandes absorbantes.
#[inline(always)]
pub fn collide<const RR: bool, const LES: bool>(f: &[f32; 9], r: Relax, q: &[f32; 9]) -> [f32; 9] {
    let [f0, f1, f2, f3, f4, f5, f6, f7, f8] = *f;
    let rho = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
    let ir = 1.0 / rho;
    let ux = (f1 - f3 + f5 - f6 - f7 + f8) * ir;
    let uy = (f2 - f4 + f5 + f6 - f7 - f8) * ir;
    let (w0, w1, w2) = (4.0 / 9.0f32, 1.0 / 9.0f32, 1.0 / 36.0f32);
    let mut p: [f32; 9];
    if RR {
        // Malaspinas (2015) : seuls ρ, u et a1 = Π − ρc_s²I − ρuu sont conservés.
        let d4 = f5 + f6 + f7 + f8;
        let rux = rho * ux;
        let ruy = rho * uy;
        let a1xx = f1 + f3 + d4 - rho / 3.0 - rux * ux;
        let a1yy = f2 + f4 + d4 - rho / 3.0 - ruy * uy;
        let a1xy = f5 - f6 + f7 - f8 - rux * uy;
        let om = if LES {
            let qn = (2.0 * (a1xx * a1xx + a1yy * a1yy + 2.0 * a1xy * a1xy)).sqrt();
            2.0 / (r.tau + (r.tau * r.tau + r.c18 * qn * ir).sqrt())
        } else {
            1.0 / r.tau
        };
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
        p = [
            w0 * (rho - 1.5 * sp),
            w1 * (tx + ax - 4.5 * qxyy),
            w1 * (ty + ay - 4.5 * qxxy),
            w1 * (tx - ax + 4.5 * qxyy),
            w1 * (ty - ay + 4.5 * qxxy),
            w2 * (td + ax + ay + pxy9 + qs),
            w2 * (td - ax + ay - pxy9 + qd),
            w2 * (td - ax - ay + pxy9 - qs),
            w2 * (td + ax - ay - pxy9 - qd),
        ];
    } else {
        let usq = 1.5 * (ux * ux + uy * uy);
        let (r0, r1, r2) = (w0 * rho, w1 * rho, w2 * rho);
        let (ux3, uy3) = (3.0 * ux, 3.0 * uy);
        let uxx = 4.5 * ux * ux - usq;
        let uyy = 4.5 * uy * uy - usq;
        let s = ux + uy;
        let ss = 4.5 * s * s - usq;
        let d = uy - ux;
        let dd = 4.5 * d * d - usq;
        let e = [
            r0 * (1.0 - usq),
            r1 * (1.0 + ux3 + uxx),
            r1 * (1.0 + uy3 + uyy),
            r1 * (1.0 - ux3 + uxx),
            r1 * (1.0 - uy3 + uyy),
            r2 * (1.0 + 3.0 * s + ss),
            r2 * (1.0 + 3.0 * d + dd),
            r2 * (1.0 - 3.0 * s + ss),
            r2 * (1.0 - 3.0 * d + dd),
        ];
        let om = if LES {
            let n: [f32; 9] = std::array::from_fn(|i| f[i] - e[i]);
            let pd = n[5] + n[6] + n[7] + n[8];
            let pxx = n[1] + n[3] + pd;
            let pyy = n[2] + n[4] + pd;
            let pxy = n[5] - n[6] + n[7] - n[8];
            let qn = (2.0 * (pxx * pxx + pyy * pyy + 2.0 * pxy * pxy)).sqrt();
            2.0 / (r.tau + (r.tau * r.tau + r.c18 * qn * ir).sqrt())
        } else {
            1.0 / r.tau
        };
        p = std::array::from_fn(|i| f[i] + om * (e[i] - f[i]));
    }
    // Forçage uniforme selon x (3 w_i c_ix g).
    let g1 = r.gx / 3.0;
    let g2 = r.gx / 12.0;
    p[1] += g1;
    p[3] -= g1;
    p[5] += g2;
    p[8] += g2;
    p[6] -= g2;
    p[7] -= g2;
    // Bandes absorbantes.
    for i in 0..9 {
        p[i] -= r.sigma * (p[i] - q[i]);
    }
    p
}
