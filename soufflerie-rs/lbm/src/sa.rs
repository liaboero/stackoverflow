//! Turbulence de Spalart–Allmaras (version standard, sans terme ft2) et loi
//! de paroi de Spalding, communes aux deux moteurs.
//!
//! Sur le réseau (Δx = Δt = 1) : advection amont, diffusion conservative,
//! destruction semi-implicite. Au premier nœud fluide, ν̃ = κ u_τ y (zone
//! logarithmique) et la paroi glisse à la vitesse u_s qui fait transmettre au
//! réseau la contrainte u_τ² de la loi de Spalding. En amont de la transition
//! (zone laminaire donnée par la couche limite intégrale), pas de production
//! turbulente et frottement laminaire τ_w = ½ Cf u_t².

use rayon::prelude::*;

use crate::geometry::{Geometry, SOLID};

pub const CB1: f32 = 0.1355;
pub const CB2: f32 = 0.622;
pub const SIGMA: f32 = 2.0 / 3.0;
pub const KAPPA: f32 = 0.41;
pub const CW2: f32 = 0.3;
pub const CW3: f32 = 2.0;
pub const CV1: f32 = 7.1;
pub const CW1: f32 = CB1 / (KAPPA * KAPPA) + (1.0 + CB2) / SIGMA;
pub const CV1_3: f32 = CV1 * CV1 * CV1;
const CW3_6: f32 = CW3 * CW3 * CW3 * CW3 * CW3 * CW3;
/// Plafond de ν̃ (unités réseau) : stabilité de la diffusion explicite.
pub const NU_TILDE_MAX: f32 = 0.15;
/// Constante de la zone logarithmique.
const B: f64 = 5.0;

/// ν̃ au large : 3ν (valeur de la NASA TMR).
pub fn nu_tilde_inf(nu: f64) -> f32 {
    (3.0 * nu) as f32
}

/// k-ε RNG (Yakhot et al. 1992) : C_μ, plafond de ν_t (réseau) et état
/// amont (Tu = 0,1 %, ν_t/ν = 0,2), comme `KE_FNS` de gpu3d.
pub const CMU: f32 = 0.0845;
pub const NUT_KE_MAX: f32 = 0.08;
pub const TU_INF: f64 = 0.001;
pub const RNUT_INF: f64 = 0.2;

/// (k∞, ε∞) amont pour la vitesse `u` et la viscosité `nu` (réseau).
pub fn ke_inf(u: f64, nu: f64) -> (f32, f32) {
    let k = 1.5 * (TU_INF * u).powi(2);
    (k as f32, (CMU as f64 * k * k / (RNUT_INF * nu)) as f32)
}

/// ν_t = C_μ k²/ε du modèle k-ε RNG.
#[inline(always)]
pub fn nu_ke(k: f32, e: f32) -> f32 {
    if k > 0.0 && e > 0.0 {
        (CMU * k * k / e).min(NUT_KE_MAX)
    } else {
        0.0
    }
}

/// VLES (`Turbulence::Vles`) : ν_t de k-ε réduit par le filtre de maille de
/// Johansen et al. (2004), F = min(1, Δ ε/k^{3/2}) (Δ = 1, pas du réseau),
/// pondéré par la fonction f_d de la DDES (0 en couche limite attachée),
/// comme `nu_vles` de `KE_FNS` (gpu3d).
pub fn nu_vles(k: f32, e: f32, fd: f32) -> f32 {
    let nr = nu_ke(k, e);
    if nr <= 0.0 {
        return 0.0;
    }
    let f = (e / (k * k.sqrt())).min(1.0);
    nr * (1.0 - fd * (1.0 - f))
}

/// Viscosité turbulente ν_t = ν̃ f_v1(ν̃/ν).
#[inline(always)]
pub fn nu_turb(v: f32, inv_nu: f32) -> f32 {
    if v > 0.0 {
        let chi = v * inv_nu;
        let c3 = chi * chi * chi;
        v * c3 / (c3 + CV1_3)
    } else {
        0.0
    }
}

/// Vitesse de frottement u_τ telle que la loi de Spalding passe par (y, u_t) :
/// y⁺ = u⁺ + e^(−κB) [e^(κu⁺) − 1 − κu⁺ − (κu⁺)²/2 − (κu⁺)³/6].
/// Valable de la sous-couche visqueuse à la zone logarithmique.
pub fn spalding_utau(ut: f64, y: f64, nu: f64) -> f64 {
    if !(ut > 1e-12) {
        return 0.0;
    }
    let kappa = KAPPA as f64;
    let e_kb = (-kappa * B).exp();
    let rey = y * ut / nu; // = u⁺ y⁺
    let mut up = if rey < 30.0 { rey.sqrt() } else { rey.ln() / kappa };
    for _ in 0..20 {
        let ku = kappa * up;
        let e = ku.exp();
        let yy = up + e_kb * (e - 1.0 - ku - ku * ku / 2.0 - ku * ku * ku / 6.0);
        let dy = 1.0 + e_kb * kappa * (e - 1.0 - ku - ku * ku / 2.0);
        let h = up * yy - rey;
        let next = up - h / (yy + up * dy);
        up = if next > 0.0 { next } else { 0.5 * up };
        if h.abs() < 1e-7 * rey {
            break;
        }
    }
    ut / up
}

/// Vitesse de glissement de la loi de paroi en une maille pariétale.
///
/// `(ux, uy)` : vitesse du fluide relative à la paroi au pas précédent ;
/// `normal` : normale sortante ; `yw` : distance à la paroi ; `cf_lam` > 0 en
/// zone laminaire. Renvoie (u_sx, u_sy, u_τ). La contrainte transmise par le
/// réseau, (ν + κ u_τ y)(u_t − u_s)/y, égale alors u_τ².
pub fn wall_slip(ux: f64, uy: f64, normal: (f64, f64), yw: f64, cf_lam: f64, nu: f64) -> (f64, f64, f64) {
    let (nxw, nyw) = normal;
    let un = ux * nxw + uy * nyw;
    let tx = ux - un * nxw;
    let ty = uy - un * nyw;
    let ut = (tx * tx + ty * ty).sqrt();
    let utau = if cf_lam > 0.0 {
        ut * (0.5 * cf_lam).sqrt()
    } else {
        spalding_utau(ut, yw, nu)
    };
    if ut < 1e-12 {
        return (0.0, 0.0, utau);
    }
    // Viscosité de paroi κ u_τ y dans les deux cas : la couche limite laminaire,
    // plus mince qu'une maille, n'est pas résolue ; seule sa contrainte compte.
    let us = (ut - utau * utau * yw / (nu + KAPPA as f64 * utau * yw)).max(0.0);
    (us * tx / ut, us * ty / ut, utau)
}

/// Champs du modèle (lignes de nx mailles).
pub struct SaFields {
    /// ν̃ courant et tampon du pas suivant.
    pub nu_tilde: Vec<f32>,
    pub nu_tilde2: Vec<f32>,
    /// Vitesse avant collision du dernier pas.
    pub usx: Vec<f32>,
    pub usy: Vec<f32>,
    /// u_τ des mailles pariétales.
    pub utau: Vec<f32>,
}

impl SaFields {
    pub fn new(n: usize) -> SaFields {
        SaFields {
            nu_tilde: vec![0.0; n],
            nu_tilde2: vec![0.0; n],
            usx: vec![0.0; n],
            usy: vec![0.0; n],
            utau: vec![0.0; n],
        }
    }

    /// ν̃ = 3ν au large, nul dans le corps.
    pub fn init(&mut self, geo: &Geometry, nu: f64) {
        let v = nu_tilde_inf(nu);
        for (t, &k) in self.nu_tilde.iter_mut().zip(&geo.kind) {
            *t = if k == SOLID { 0.0 } else { v };
        }
        self.utau.iter_mut().for_each(|u| *u = 0.0);
    }

    /// Un pas explicite de Spalart–Allmaras. `sigma_row` : rappel des bandes
    /// absorbantes (vers ν̃∞). Les intervalles de mailles fluides ordinaires
    /// passent par une boucle sans branchement (vectorisée), les mailles
    /// pariétales et les bords à part.
    pub fn step(&mut self, geo: &Geometry, nu: f64, sigma_row: &[f32]) {
        let (nx, ny) = (geo.nx, geo.ny);
        let SaFields {
            nu_tilde: vt,
            nu_tilde2,
            usx: ux,
            usy: uy,
            utau,
        } = self;
        let (vt, ux, uy, utau) = (&*vt, &*ux, &*uy, &*utau);
        let nu = nu as f32;
        let inf = 3.0 * nu;
        nu_tilde2.par_chunks_mut(nx).enumerate().for_each(|(y, out)| {
            let row = y * nx;
            if y == 0 || y == ny - 1 {
                // Parois haute et basse glissantes : ν̃ de la ligne intérieure.
                let r = if y == 0 { row + nx } else { row - nx };
                out.copy_from_slice(&vt[r..r + nx]);
                return;
            }
            // Solides : 0 ; puis entrée, sortie, intervalles fluides, paroi.
            out.iter_mut().for_each(|v| *v = 0.0);
            out[0] = inf;
            out[nx - 1] = vt[row + nx - 2];
            let c = Rows {
                vt: &vt[row..row + nx],
                vd: &vt[row - nx..row],
                vu: &vt[row + nx..row + 2 * nx],
                ux: &ux[row..row + nx],
                uy: &uy[row..row + nx],
                uxd: &ux[row - nx..row],
                uxu: &ux[row + nx..row + 2 * nx],
                dist: &geo.dist[row..row + nx],
                lam: &geo.lam[row..row + nx],
            };
            let sg = sigma_row[y];
            for &(x0, x1) in &geo.runs[y] {
                c.run(&mut out[x0..x1], x0, sg, nu, inf);
            }
            for &x in &geo.near[y] {
                let idx = row + x;
                out[x] = if geo.is_wall[idx] {
                    // Premier nœud : ν̃ ≈ ν_t = κ u_τ y (zone logarithmique).
                    (KAPPA * utau[idx] * geo.dist[idx].max(0.1)).min(NU_TILDE_MAX)
                } else {
                    let mut one = [0.0];
                    c.run(&mut one, x, sg, nu, inf);
                    one[0]
                };
            }
        });
        std::mem::swap(&mut self.nu_tilde, &mut self.nu_tilde2);
    }
}

/// Lignes y − 1, y, y + 1 des champs lus par un pas.
struct Rows<'a> {
    vt: &'a [f32],
    vd: &'a [f32],
    vu: &'a [f32],
    ux: &'a [f32],
    uy: &'a [f32],
    uxd: &'a [f32],
    uxu: &'a [f32],
    dist: &'a [f32],
    lam: &'a [f32],
}

impl Rows<'_> {
    /// Mailles x0 .. x0 + out.len() (1 ≤ x0, x0 + len ≤ nx − 1).
    #[inline(always)]
    fn run(&self, out: &mut [f32], x0: usize, sg: f32, nu: f32, inf: f32) {
        let inv_nu = 1.0 / nu;
        let n = out.len();
        let x1 = x0 + n;
        let (v, vl, vr) = (&self.vt[x0..x1], &self.vt[x0 - 1..x1 - 1], &self.vt[x0 + 1..x1 + 1]);
        let (vd, vu) = (&self.vd[x0..x1], &self.vu[x0..x1]);
        let (a, b) = (&self.ux[x0..x1], &self.uy[x0..x1]);
        let (uyl, uyr) = (&self.uy[x0 - 1..x1 - 1], &self.uy[x0 + 1..x1 + 1]);
        let (uxd, uxu) = (&self.uxd[x0..x1], &self.uxu[x0..x1]);
        let (dist, lam) = (&self.dist[x0..x1], &self.lam[x0..x1]);
        for k in 0..n {
            let omega = (0.5 * (uyr[k] - uyl[k]) - 0.5 * (uxu[k] - uxd[k])).abs();
            out[k] = update(
                v[k], vl[k], vr[k], vd[k], vu[k], a[k], b[k], omega, dist[k], lam[k], sg, nu, inv_nu, inf,
            );
        }
    }
}

/// z^(−1/6) pour z ≥ ~1 : estimation par les bits de l'exposant (à ~10 %),
/// puis quatre itérations de Newton sans division, y ← y (7 − z y⁶)/6.
#[inline(always)]
fn inv_pow_1_6(z: f32) -> f32 {
    // Bits de y ≈ (7/6)·bits(1) − bits(z)/6, calculés en flottant (vectorisable).
    const ONE_7_6: f32 = (0x3f80_0000u64 * 7 / 6) as f32;
    let mut y = f32::from_bits((ONE_7_6 - (z.to_bits() as i32) as f32 * (1.0 / 6.0)) as i32 as u32);
    for _ in 0..4 {
        let y2 = y * y;
        y *= (7.0 - z * y2 * y2 * y2) * (1.0 / 6.0);
    }
    y
}

/// Nouveau ν̃ d'une maille courante : advection amont, diffusion
/// conservative, production (nulle en zone laminaire), destruction
/// semi-implicite, rappel des bandes absorbantes.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn update(
    v: f32,
    vl: f32,
    vr: f32,
    vd: f32,
    vu: f32,
    a: f32,
    b: f32,
    omega: f32,
    dist: f32,
    lam: f32,
    sg: f32,
    nu: f32,
    inv_nu: f32,
    inf: f32,
) -> f32 {
    let adv = a.max(0.0) * (v - vl) + a.min(0.0) * (vr - v) + b.max(0.0) * (v - vd) + b.min(0.0) * (vu - v);
    let mut diff = (nu + 0.5 * (v + vr)) * (vr - v) - (nu + 0.5 * (v + vl)) * (v - vl)
        + (nu + 0.5 * (v + vu)) * (vu - v)
        - (nu + 0.5 * (v + vd)) * (v - vd);
    let gx = 0.5 * (vr - vl);
    let gy = 0.5 * (vu - vd);
    diff = (diff + CB2 * (gx * gx + gy * gy)) * (1.0 / SIGMA);
    let d = dist.max(0.5);
    let id2 = 1.0 / (d * d);
    let ikd2 = id2 * (1.0 / (KAPPA * KAPPA));
    let chi = v * inv_nu;
    let c3 = chi * chi * chi;
    // f_v2 = 1 − χ/(1 + χ f_v1), f_v1 = χ³/(χ³ + c_v1³), en une division.
    let fv2 = 1.0 - chi * (c3 + CV1_3) / (c3 * (1.0 + chi) + CV1_3);
    let s_tilde = (omega + v * ikd2 * fv2).max(0.3 * omega);
    let prod = if lam > 0.0 { 0.0 } else { CB1 * s_tilde * v };
    let r = if s_tilde > 1e-12 {
        (v * ikd2 / s_tilde).min(10.0)
    } else {
        10.0
    };
    let r2 = r * r;
    let g = r + CW2 * (r2 * r2 * r2 - r);
    let g2 = g * g;
    // f_w = g [(1 + c_w3⁶)/(g⁶ + c_w3⁶)]^(1/6).
    let fw = g * inv_pow_1_6((g2 * g2 * g2 + CW3_6) * (1.0 / (1.0 + CW3_6)));
    let mut nv = (v + prod + diff - adv) / (1.0 + CW1 * fw * v * id2);
    nv -= sg * (nv - inf);
    nv.clamp(0.0, NU_TILDE_MAX)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn racine_sixieme() {
        for k in 0..=3000 {
            let z = 0.98f32 * 10f32.powf(k as f32 * 0.01);
            let y = inv_pow_1_6(z);
            let exact = (z as f64).powf(-1.0 / 6.0);
            let e = (y as f64 - exact).abs() / exact;
            assert!(e < 3e-7, "z = {z} : {y} (erreur {e})");
        }
    }

    #[test]
    fn k_eps_amont_et_paroi() {
        // Amont : ν_t/ν = 0,2.
        let (u, nu) = (0.1, 2e-6);
        let (k, e) = ke_inf(u, nu);
        assert!((nu_ke(k, e) as f64 / nu - 0.2).abs() < 1e-3);
        // Paroi : k = u_τ²/√C_μ, ε = u_τ³/(κ h) ⇒ ν_t = κ u_τ h.
        let (ut, h) = (0.004f32, 2.0f32);
        let (k, e) = (ut * ut / CMU.sqrt(), ut * ut * ut / (KAPPA * h));
        assert!((nu_ke(k, e) - KAPPA * ut * h).abs() < 1e-6);
    }

    #[test]
    fn vles_filtre() {
        // Zone logarithmique à y = 10 mailles : l = k^{3/2}/ε = κ y/C_μ^{3/4}
        // dépasse la maille ; protégée (f_d = 0) le RANS reste entier, sinon
        // ν_t → C_μ Δ √k = C_μ^{3/4} u_τ.
        let (ut, y) = (0.004f32, 10.0f32);
        let (k, e) = (ut * ut / CMU.sqrt(), ut * ut * ut / (KAPPA * y));
        assert!((nu_vles(k, e, 0.0) - nu_ke(k, e)).abs() < 1e-9);
        assert!((nu_vles(k, e, 1.0) - CMU.powf(0.75) * ut).abs() < 1e-7);
        // Amont : échelle k^{3/2}/ε ≪ Δ, aucun effet.
        let (k, e) = ke_inf(0.1, 2e-6);
        assert_eq!(nu_vles(k, e, 1.0), nu_ke(k, e));
    }

    #[test]
    fn spalding_inverse() {
        // Sous-couche visqueuse : u⁺ = y⁺ ⇒ u_τ = √(ν u / y).
        let (nu, y) = (1e-5, 0.5);
        let ut = 1e-4;
        let utau = spalding_utau(ut, y, nu);
        assert!((utau - (nu * ut / y).sqrt()).abs() / utau < 0.02);
        // Zone logarithmique : u⁺ = ln(y⁺)/κ + B.
        let nu = 1e-7;
        let ut = 0.1;
        let utau = spalding_utau(ut, y, nu);
        let yp = y * utau / nu;
        let up = ut / utau;
        assert!((up - (yp.ln() / KAPPA as f64 + B)).abs() < 0.05, "u+ {up} y+ {yp}");
    }
}
