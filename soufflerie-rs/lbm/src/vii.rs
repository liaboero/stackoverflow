//! Couplage visqueux–non visqueux simultané, résolu par Newton (méthode de
//! XFOIL : Drela & Giles 1987, Drela 1989).
//!
//! - Écoulement potentiel : les panneaux tourbillons de [`crate::panel`]
//!   (Kuethe & Chow), plus des sources uniformes sur chaque panneau du profil
//!   et du sillage. La densité de source est la transpiration d(ue δ*)/ds :
//!   ue = ue_inv + D·m, m = ue δ* (défaut de débit), D matrice d'influence
//!   (sources → tourbillons via la condition de glissement et de Kutta → vitesse).
//! - Sillage : ligne de courant non visqueuse issue du bord de fuite, une
//!   corde de long.
//! - Couche limite : équations intégrales de quantité de mouvement et
//!   d'énergie cinétique, différences logarithmiques de XFOIL (BLDIF) ;
//!   fermetures de Drela & Giles (laminaire) et de Drela (turbulent :
//!   frottement de Swafford, H* , dissipation, contrainte de cisaillement
//!   retardée de Green) ; amplification eᴺ de Drela & Giles ; point de
//!   transition interpolé dans l'intervalle (TRDIF).
//! - Newton sur toutes les inconnues (θ, m, N ou √Cτ à chaque station,
//!   profil et sillage) : ue n'est figé nulle part, le bord de fuite fermé
//!   est traité comme le reste.
//! - Traînée : Squire–Young au bout du sillage.

use crate::panel::{airfoil_points, leading_panel, Friction, PanelOptions, PanelResult, PanelShape};
use std::f64::consts::PI;

/// Panneaux du sillage, longueur du sillage (cordes).
const NW: usize = 40;
const WAKE_LEN: f64 = 1.0;
const MAX_IT: usize = 60;

// ------------------------------------------------------------- algèbre

/// Factorisation LU à pivot partiel, matrice carrée rangée par lignes.
struct Lu {
    n: usize,
    a: Vec<f64>,
    p: Vec<usize>,
}

impl Lu {
    fn new(mut a: Vec<f64>, n: usize) -> Lu {
        let mut p: Vec<usize> = (0..n).collect();
        for c in 0..n {
            let mut piv = c;
            for r in c + 1..n {
                if a[r * n + c].abs() > a[piv * n + c].abs() {
                    piv = r;
                }
            }
            if piv != c {
                for k in 0..n {
                    a.swap(c * n + k, piv * n + k);
                }
                p.swap(c, piv);
            }
            let d = a[c * n + c];
            for r in c + 1..n {
                let f = a[r * n + c] / d;
                a[r * n + c] = f;
                if f == 0.0 {
                    continue;
                }
                let (top, bot) = a.split_at_mut(r * n);
                let rc = &top[c * n..c * n + n];
                let rr = &mut bot[..n];
                for k in c + 1..n {
                    rr[k] -= f * rc[k];
                }
            }
        }
        Lu { n, a, p }
    }

    fn solve(&self, b: &[f64]) -> Vec<f64> {
        let n = self.n;
        let mut x: Vec<f64> = self.p.iter().map(|&i| b[i]).collect();
        for r in 0..n {
            let row = &self.a[r * n..r * n + r];
            let s: f64 = row.iter().zip(&x[..r]).map(|(a, b)| a * b).sum();
            x[r] -= s;
        }
        for r in (0..n).rev() {
            let row = &self.a[r * n + r + 1..r * n + n];
            let s: f64 = row.iter().zip(&x[r + 1..]).map(|(a, b)| a * b).sum();
            x[r] = (x[r] - s) / self.a[r * n + r];
        }
        x
    }
}

// ------------------------------------------------------------ singularités

/// Panneau tourbillon linéaire j (Kuethe & Chow) vu du point P, composantes
/// normale (gauche) et tangentielle selon la direction `tht` : [cn1, cn2, ct1, ct2].
#[allow(clippy::too_many_arguments)]
fn kc(px: f64, py: f64, tht: f64, x0: f64, y0: f64, thj: f64, sj: f64) -> [f64; 4] {
    let dx = px - x0;
    let dy = py - y0;
    let a = -dx * thj.cos() - dy * thj.sin();
    let b = dx * dx + dy * dy;
    let c = (tht - thj).sin();
    let d = (tht - thj).cos();
    let e = dx * thj.sin() - dy * thj.cos();
    let f = (1.0 + sj * (sj + 2.0 * a) / b).ln();
    let g = (e * sj).atan2(b + a * sj);
    let p = dx * (tht - 2.0 * thj).sin() + dy * (tht - 2.0 * thj).cos();
    let q = dx * (tht - 2.0 * thj).cos() - dy * (tht - 2.0 * thj).sin();
    let cn2 = d + 0.5 * q * f / sj - (a * c + d * e) * g / sj;
    let cn1 = 0.5 * d * f + c * g - cn2;
    let ct2 = c + 0.5 * p * f / sj + (a * d - c * e) * g / sj;
    let ct1 = 0.5 * c * f - d * g - ct2;
    [cn1, cn2, ct1, ct2]
}

/// Source à densité linéaire sur le segment A→B (σa en A, σb en B) : vitesses
/// (u, v) en P pour σa = 1 et pour σb = 1. `own` : P est le milieu du segment,
/// vu du côté droit (intérieur du profil).
fn source_lin(p: [f64; 2], a: [f64; 2], b: [f64; 2], own: bool) -> [[f64; 2]; 2] {
    let l = (b[0] - a[0]).hypot(b[1] - a[1]);
    let (tx, ty) = ((b[0] - a[0]) / l, (b[1] - a[1]) / l);
    let (nx, ny) = (-ty, tx);
    let (rx, ry) = (p[0] - a[0], p[1] - a[1]);
    let xi = rx * tx + ry * ty;
    let eta = rx * nx + ry * ny;
    let (l0, beta) = if own {
        (0.0, -PI)
    } else {
        let r1 = xi * xi + eta * eta;
        let r2 = (xi - l) * (xi - l) + eta * eta;
        (0.5 * (r1 / r2).ln(), (eta * l).atan2(xi * (xi - l) + eta * eta))
    };
    let eta = if own { 0.0 } else { eta };
    // ∫σ(ξ−s)/r² ds et ∫σ η/r² ds, σ = σa + (σb − σa) s/l.
    let l1 = xi * l0 - l + eta * beta;
    let b1 = xi * beta - eta * l0;
    let k = 1.0 / (2.0 * PI);
    let (ua, ub) = (k * (l0 - l1 / l), k * l1 / l);
    let (va, vb) = (k * (beta - b1 / l), k * b1 / l);
    [[ua * tx + va * nx, ua * ty + va * ny], [ub * tx + vb * nx, ub * ty + vb * ny]]
}

// ------------------------------------------------------------ modèle

/// Profil, sillage et matrices d'influence pour une incidence.
struct Model {
    n: usize,
    /// Nœuds du profil (utilisés par les essais).
    #[cfg_attr(not(test), allow(dead_code))]
    pts: Vec<[f64; 2]>,
    xc: Vec<f64>,
    yc: Vec<f64>,
    th: Vec<f64>,
    len: Vec<f64>,
    /// Abscisse curviligne des nœuds (n+1) et des points de contrôle (n).
    c_node: Vec<f64>,
    c_ctrl: Vec<f64>,
    /// Distance au bord de fuite du milieu de chaque panneau de sillage, longueurs.
    wd: Vec<f64>,
    wl: Vec<f64>,
    wx: Vec<f64>,
    /// Vitesse non visqueuse aux stations : sens du contour sur le profil,
    /// sens de l'écoulement dans le sillage.
    q_inv: Vec<f64>,
    /// dq_i/dσ_j (σ aux nœuds du profil puis du sillage), (n+nw) × (n+nw+2), par lignes.
    b: Vec<f64>,
}

impl Model {
    fn new(pts: Vec<[f64; 2]>, alpha: f64) -> Model {
        let n = pts.len() - 1;
        let (mut xc, mut yc, mut th, mut len) = (vec![0.0; n], vec![0.0; n], vec![0.0; n], vec![0.0; n]);
        for j in 0..n {
            let ([x0, y0], [x1, y1]) = (pts[j], pts[j + 1]);
            xc[j] = 0.5 * (x0 + x1);
            yc[j] = 0.5 * (y0 + y1);
            len[j] = (x1 - x0).hypot(y1 - y0);
            th[j] = (y1 - y0).atan2(x1 - x0);
        }
        let mut c_node = vec![0.0; n + 1];
        for j in 0..n {
            c_node[j + 1] = c_node[j] + len[j];
        }
        let c_ctrl: Vec<f64> = (0..n).map(|j| c_node[j] + 0.5 * len[j]).collect();
        // Matrices des tourbillons aux points de contrôle.
        let m1 = n + 1;
        let mut an = vec![0.0; m1 * m1];
        let mut at = vec![0.0; n * m1];
        for i in 0..n {
            for j in 0..n {
                let [cn1, cn2, ct1, ct2] = if i == j {
                    [-1.0, 1.0, PI / 2.0, PI / 2.0]
                } else {
                    kc(xc[i], yc[i], th[i], pts[j][0], pts[j][1], th[j], len[j])
                };
                an[i * m1 + j] += cn1;
                an[i * m1 + j + 1] += cn2;
                at[i * m1 + j] += ct1;
                at[i * m1 + j + 1] += ct2;
            }
        }
        an[n * m1] = 1.0;
        an[n * m1 + n] = 1.0;
        let lu = Lu::new(an, m1);
        let mut rhs: Vec<f64> = (0..n).map(|i| (th[i] - alpha).sin()).collect();
        rhs.push(0.0);
        let g0 = lu.solve(&rhs);
        let vel = |px: f64, py: f64| {
            let (mut u, mut v) = (alpha.cos(), alpha.sin());
            for j in 0..n {
                let [cn1, cn2, ct1, ct2] = kc(px, py, 0.0, pts[j][0], pts[j][1], th[j], len[j]);
                u += ct1 * g0[j] + ct2 * g0[j + 1];
                v += cn1 * g0[j] + cn2 * g0[j + 1];
            }
            (u, v)
        };
        // Sillage : bissectrice du bord de fuite, puis ligne de courant.
        let te = pts[0];
        let t0 = (th[0].cos(), th[0].sin());
        let t1 = (th[n - 1].cos(), th[n - 1].sin());
        let (bx, by) = (t1.0 - t0.0, t1.1 - t0.1);
        let bl = bx.hypot(by);
        let mut dir = (bx / bl, by / bl);
        let l0 = 0.5 * (len[0] + len[n - 1]);
        let total = |r: f64| l0 * (r.powi(NW as i32) - 1.0) / (r - 1.0);
        let (mut lo, mut hi) = (1.0 + 1e-9, 2.0);
        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            if total(mid) > WAKE_LEN {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let ratio = 0.5 * (lo + hi);
        let mut wk = vec![te];
        let mut wl = Vec::with_capacity(NW);
        let mut p = te;
        for k in 0..NW {
            let l = l0 * ratio.powi(k as i32);
            if (p[0] - te[0]).hypot(p[1] - te[1]) > 0.02 {
                let (u, v) = vel(p[0] + 0.5 * l * dir.0, p[1] + 0.5 * l * dir.1);
                let q = u.hypot(v).max(1e-12);
                dir = (u / q, v / q);
            }
            p = [p[0] + l * dir.0, p[1] + l * dir.1];
            wk.push(p);
            wl.push(l);
        }
        let mut wd = vec![0.0; NW];
        let mut acc = 0.0;
        for k in 0..NW {
            wd[k] = acc + 0.5 * wl[k];
            acc += wl[k];
        }
        let wmid: Vec<[f64; 2]> = (0..NW)
            .map(|k| [0.5 * (wk[k][0] + wk[k + 1][0]), 0.5 * (wk[k][1] + wk[k + 1][1])])
            .collect();
        let wth: Vec<f64> = (0..NW).map(|k| (wk[k + 1][1] - wk[k][1]).atan2(wk[k + 1][0] - wk[k][0])).collect();
        let mut cw = vec![0.0; NW * m1];
        for k in 0..NW {
            for j in 0..n {
                let [_, _, ct1, ct2] = kc(wmid[k][0], wmid[k][1], wth[k], pts[j][0], pts[j][1], th[j], len[j]);
                cw[k * m1 + j] += ct1;
                cw[k * m1 + j + 1] += ct2;
            }
        }
        let ns = n + NW;
        let mut q_inv = vec![0.0; ns];
        for i in 0..n {
            q_inv[i] = (th[i] - alpha).cos() + (0..m1).map(|j| at[i * m1 + j] * g0[j]).sum::<f64>();
        }
        for k in 0..NW {
            q_inv[n + k] = (wth[k] - alpha).cos() + (0..m1).map(|j| cw[k * m1 + j] * g0[j]).sum::<f64>();
        }
        // Influence des sources (densité linéaire continue, valeurs aux nœuds du
        // profil puis du sillage) : la condition de glissement redistribue les tourbillons.
        let nsrc = n + 1 + NW + 1;
        let mut b = vec![0.0; ns * nsrc];
        for jj in 0..nsrc {
            // Panneaux adjacents au nœud : (A, B, rang du nœud : 0 en A, 1 en B, indice de station du panneau).
            let mut adj: Vec<([f64; 2], [f64; 2], usize, usize)> = Vec::new();
            if jj <= n {
                if jj >= 1 {
                    adj.push((pts[jj - 1], pts[jj], 1, jj - 1));
                }
                if jj < n {
                    adj.push((pts[jj], pts[jj + 1], 0, jj));
                }
            } else {
                let k = jj - n - 1;
                if k >= 1 {
                    adj.push((wk[k - 1], wk[k], 1, n + k - 1));
                }
                if k < NW {
                    adj.push((wk[k], wk[k + 1], 0, n + k));
                }
            }
            let mut rhs = vec![0.0; m1];
            let mut tan = vec![0.0; ns];
            for &(a, bb, end, own_st) in &adj {
                for i in 0..n {
                    let [u, v] = source_lin([xc[i], yc[i]], a, bb, i == own_st)[end];
                    rhs[i] -= -th[i].sin() * u + th[i].cos() * v;
                    tan[i] += th[i].cos() * u + th[i].sin() * v;
                }
                for k in 0..NW {
                    let [u, v] = source_lin(wmid[k], a, bb, n + k == own_st)[end];
                    tan[n + k] += wth[k].cos() * u + wth[k].sin() * v;
                }
            }
            let g = lu.solve(&rhs);
            for i in 0..n {
                b[i * nsrc + jj] = tan[i] + (0..m1).map(|j| at[i * m1 + j] * g[j]).sum::<f64>();
            }
            for k in 0..NW {
                b[(n + k) * nsrc + jj] = tan[n + k] + (0..m1).map(|j| cw[k * m1 + j] * g[j]).sum::<f64>();
            }
        }
        let wx = wmid.iter().map(|p| p[0]).collect();
        Model {
            n,
            pts,
            xc,
            yc,
            th,
            len,
            c_node,
            c_ctrl,
            wd,
            wl,
            wx,
            q_inv,
            b,
        }
    }

    fn ns(&self) -> usize {
        self.n + NW
    }

    fn nsrc(&self) -> usize {
        self.n + NW + 2
    }

    /// Densités de source aux nœuds (profil puis sillage) d'une distribution de
    /// m aux stations : différences compactes entre stations voisines, σ = dm/dξ,
    /// le point d'arrêt (m = 0) étant à l'abscisse curviligne `cs`, entre les
    /// points de contrôle `ist` et `ist + 1`.
    fn sigma(&self, m: &[f64], ist: usize, cs: f64) -> Vec<f64> {
        let n = self.n;
        let mut s = vec![0.0; self.nsrc()];
        for k in 1..n {
            let (ca, cb) = (self.c_ctrl[k - 1], self.c_ctrl[k]);
            s[k] = if k - 1 == ist {
                if self.c_node[k] <= cs {
                    m[k - 1] / (cs - ca).max(1e-12)
                } else {
                    m[k] / (cb - cs).max(1e-12)
                }
            } else if k - 1 < ist {
                (m[k - 1] - m[k]) / (cb - ca)
            } else {
                (m[k] - m[k - 1]) / (cb - ca)
            };
        }
        // Bord de fuite : m constant sur le dernier demi-panneau.
        s[0] = 0.0;
        s[n] = 0.0;
        let w0 = n + 1;
        s[w0] = (m[n] - m[0] - m[n - 1]) / (0.5 * self.wl[0]);
        for k in 1..NW {
            s[w0 + k] = (m[n + k] - m[n + k - 1]) / (0.5 * (self.wl[k - 1] + self.wl[k]));
        }
        s[w0 + NW] = s[w0 + NW - 1];
        s
    }

    /// D = dq/dm (q : vitesses signées aux stations), par lignes.
    fn mass_matrix(&self, ist: usize, cs: f64) -> Vec<f64> {
        let ns = self.ns();
        let nsrc = self.nsrc();
        let mut d = vec![0.0; ns * ns];
        let mut unit = vec![0.0; ns];
        for k in 0..ns {
            unit[k] = 1.0;
            let s = self.sigma(&unit, ist, cs);
            unit[k] = 0.0;
            for (j, &sj) in s.iter().enumerate() {
                if sj == 0.0 {
                    continue;
                }
                for i in 0..ns {
                    d[i * ns + k] += self.b[i * nsrc + j] * sj;
                }
            }
        }
        d
    }
}

// ------------------------------------------------------------ fermetures

#[derive(Clone, Copy, PartialEq, Debug)]
enum Kind {
    Lam,
    Turb,
    Wake,
}

#[derive(Clone, Copy, Default, Debug)]
struct Cl {
    th: f64,
    ds: f64,
    ue: f64,
    c: f64,
    h: f64,
    hk: f64,
    hs: f64,
    us: f64,
    cf: f64,
    /// 2 CD / H*.
    di: f64,
    /// √Cτ d'équilibre.
    cq: f64,
    de: f64,
    rt: f64,
    /// dN/ds (laminaire).
    ax: f64,
}

/// Constantes de Drela (1989) : G–β (A = 6,7, B = 0,75), K_C = 5,6.
const GACON: f64 = 6.70;
const GBCON: f64 = 0.75;
const GCCON: f64 = 18.0;
const SCCON: f64 = 5.6;
const DLCON: f64 = 0.9;
const CTRCON: f64 = 1.8;
const CTRCEX: f64 = 3.3;

fn ctcon() -> f64 {
    0.5 / (GACON * GACON * GBCON)
}

fn hs_lam(hk: f64) -> f64 {
    if hk < 4.0 {
        1.515 + 0.076 * (4.0 - hk).powi(2) / hk
    } else {
        1.515 + 0.040 * (hk - 4.0).powi(2) / hk
    }
}

fn cf_lam(hk: f64, rt: f64) -> f64 {
    if hk < 5.5 {
        (0.0727 * (5.5 - hk).powi(3) / (hk + 1.0) - 0.07) / rt
    } else {
        (0.015 * (1.0 - 1.0 / (hk - 4.5)).powi(2) - 0.07) / rt
    }
}

fn di_lam(hk: f64, rt: f64) -> f64 {
    if hk < 4.0 {
        (0.207 + 0.00205 * (4.0 - hk).powf(5.5)) / rt
    } else {
        let d = (hk - 4.0).powi(2);
        (0.207 - 0.003 * d / (1.0 + 0.02 * d)) / rt
    }
}

fn hs_turb(hk: f64, rt: f64) -> f64 {
    let (hsmin, dhsinf) = (1.5, 0.015);
    let ho = if rt > 400.0 { 3.0 + 400.0 / rt } else { 4.0 };
    let rtz = rt.max(200.0);
    if hk < ho {
        let hr = (ho - hk) / (ho - 1.0);
        (2.0 - hsmin - 4.0 / rtz) * hr * hr * 1.5 / (hk + 0.5) + hsmin + 4.0 / rtz
    } else {
        let grt = rtz.ln();
        let hdif = hk - ho;
        let rtmp = hk - ho + 4.0 / grt;
        let htmp = 0.007 * grt / (rtmp * rtmp) + dhsinf / hk;
        hdif * hdif * htmp + hsmin + 4.0 / rtz
    }
}

fn cf_turb(hk: f64, rt: f64) -> f64 {
    let grt = rt.max(1e-9).ln().max(3.0);
    let gex = -1.74 - 0.31 * hk;
    let arg = (-1.33 * hk).max(-20.0);
    let thk = (4.0 - hk / 0.875).tanh();
    0.3 * arg.exp() * (grt / std::f64::consts::LN_10).powf(gex) + 1.1e-4 * (thk - 1.0)
}

/// Taux d'amplification de l'enveloppe eᴺ (Drela & Giles 1987), avec la
/// rampe de XFOIL autour du Re_θ critique (±0,08 en log₁₀).
fn amp_rate(hk: f64, th: f64, rt: f64) -> f64 {
    let hm = (hk - 1.0).max(0.05);
    let log10_re0 = (1.415 / hm - 0.489) * (20.0 / hm - 12.9).tanh() + 3.295 / hm + 0.44;
    let dgr = 0.08;
    let gr = rt.max(1.0).log10();
    if gr < log10_re0 - dgr {
        return 0.0;
    }
    let rn = (gr - (log10_re0 - dgr)) / (2.0 * dgr);
    let rfac = if rn >= 1.0 { 1.0 } else { 3.0 * rn * rn - 2.0 * rn * rn * rn };
    let dn_dre = 0.01 * ((2.4 * hk - 3.7 + 2.5 * (1.5 * hk - 4.65).tanh()).powi(2) + 0.25).sqrt();
    let ell = (6.54 * hk - 14.07) / (hk * hk);
    let m = (0.058 * (hk - 4.0).powi(2) / hm - 0.068) / ell;
    rfac * dn_dre * 0.5 * (m + 1.0) * ell / th
}

fn closure(kind: Kind, th: f64, ds: f64, c: f64, ue: f64, re: f64) -> Cl {
    let th = th.max(1e-12);
    let ue = ue.max(1e-5);
    let h = ds / th;
    let hk = if kind == Kind::Wake { h.max(1.00005) } else { h.max(1.05) };
    let rt = (re * ue * th).max(1e-6);
    let mut k = Cl {
        th,
        ds,
        ue,
        c,
        h,
        hk,
        rt,
        ..Cl::default()
    };
    match kind {
        Kind::Lam => {
            k.hs = hs_lam(hk);
            k.cf = cf_lam(hk, rt);
            k.di = di_lam(hk, rt);
            k.ax = amp_rate(hk, th, rt);
        }
        Kind::Turb | Kind::Wake => {
            let wake = kind == Kind::Wake;
            let hh = h.max(1.00005);
            k.hs = hs_turb(hk, rt);
            // Frottement turbulent borné par la valeur laminaire (faibles Re_θ, comme XFOIL).
            k.cf = if wake { 0.0 } else { cf_turb(hk, rt).max(cf_lam(hk, rt)) };
            k.us = (0.5 * k.hs * (1.0 - (hk - 1.0) / (GBCON * hh))).min(if wake { 0.99995 } else { 0.98 });
            let hkc = if wake { hk - 1.0 } else { (hk - 1.0 - GCCON / rt).max(0.01) };
            k.cq = (ctcon() * k.hs * (hk - 1.0) * hkc * hkc / ((1.0 - k.us) * hh * hk * hk)).max(0.0).sqrt();
            let s = c.max(1e-8);
            k.c = s;
            if wake {
                k.di = 2.0 * (s * s * (1.0 - k.us)) * 2.0 / k.hs;
            } else {
                k.di = ((0.5 * k.cf * k.us + s * s * (1.0 - k.us)) * 2.0 / k.hs).max(di_lam(hk, rt));
            }
            k.de = ((3.15 + 1.72 / (hk - 1.0)) * th + ds).min(12.0 * th);
        }
    }
    k
}

/// Point de transition dans l'intervalle 1 → 2 (TRCHEK2 de XFOIL) : N₂ est
/// défini implicitement, N₂ = N₁ + ⟨dN/ds⟩₁ₜ (x₂ − x₁), le taux étant moyenné
/// entre la station 1 (laminaire) et le point T où N atteint N_crit, état
/// interpolé entre 1 et 2. Renvoie la fraction w de l'intervalle (1 : pas de
/// transition libre) et si la transition a lieu ; `trip` (x/c) l'avance.
#[allow(clippy::too_many_arguments)]
fn transition_weight(a: &Cl, s2: [f64; 3], x1: f64, x2: f64, xs1: f64, xs2: f64, trip: f64, re: f64, ncrit: f64) -> (f64, bool) {
    let dx = x2 - x1;
    let res = |n2: f64| {
        let w = if n2 <= ncrit || n2 <= a.c { 1.0 } else { ((ncrit - a.c) / (n2 - a.c)).clamp(0.0, 1.0) };
        let t = closure(
            Kind::Lam,
            a.th + w * (s2[0] - a.th),
            a.ds + w * (s2[1] - a.ds),
            0.0,
            a.ue + w * (s2[2] - a.ue),
            re,
        );
        let axa = (0.5 * (a.ax * a.ax + t.ax * t.ax)).sqrt();
        (n2 - a.c - axa * dx, w)
    };
    // Départ : taux de la station 1 seule, puis Newton (dérivée numérique).
    let mut n2 = a.c + a.ax * dx;
    for _ in 0..40 {
        let (r, _) = res(n2);
        let h = 1e-7 * n2.abs().max(1.0);
        let dr = (res(n2 + h).0 - r) / h;
        let step = if dr.abs() > 1e-12 { -r / dr } else { -r };
        n2 += step.clamp(-5.0, 5.0);
        if step.abs() < 1e-13 * n2.abs().max(1.0) {
            break;
        }
    }
    let (_, mut w) = res(n2);
    let mut found = n2 > ncrit;
    if !found {
        w = 1.0;
    }
    if trip.is_finite() && trip < xs2 {
        found = true;
        let wt = if trip <= xs1 { 0.0 } else { (trip - xs1) / (xs2 - xs1) };
        w = w.min(wt);
    }
    (w.clamp(0.0, 1.0), found)
}

/// √Cτ initial au point de transition (XFOIL : 1,8 e^(−3,3/(Hk−1)) √Cτ_eq).
fn s_transition(th: f64, ds: f64, ue: f64, re: f64) -> f64 {
    let t = closure(Kind::Turb, th, ds, 1.0, ue, re);
    (CTRCON * (-CTRCEX / (t.hk - 1.0)).exp() * t.cq).max(1e-6)
}

/// Équations d'un intervalle 1 → 2 (BLDIF de XFOIL, sans compressibilité) :
/// quantité de mouvement, énergie cinétique, amplification ou retard de Cτ.
fn interval(kind: Kind, x1: f64, x2: f64, a: &Cl, b: &Cl) -> [f64; 3] {
    let xlog = (x2 / x1).ln();
    let ulog = (b.ue / a.ue).ln();
    let tlog = (b.th / a.th).ln();
    let hlog = (b.hs / a.hs).ln();
    let hl = ((b.hk - 1.0) / (a.hk - 1.0)).abs().max(1e-30).ln();
    let hdcon = if kind == Kind::Wake { 1.0 } else { 5.0 } / (b.hk * b.hk);
    let upw = 1.0 - 0.5 * (-(hl * hl).min(15.0) * hdcon).exp();
    let xa = 0.5 * (x1 + x2);
    let ta = 0.5 * (a.th + b.th);
    let ha = 0.5 * (a.h + b.h);
    let hka = 0.5 * (a.hk + b.hk);
    let rta = 0.5 * (a.rt + b.rt);
    let cfm = match kind {
        Kind::Lam => cf_lam(hka, rta),
        Kind::Turb => cf_turb(hka, rta).max(cf_lam(hka, rta)),
        Kind::Wake => 0.0,
    };
    let (xt1, xt2) = (x1 / a.th, x2 / b.th);
    let cfx = 0.5 * cfm * xa / ta + 0.25 * (a.cf * xt1 + b.cf * xt2);
    let r1 = tlog + (ha + 2.0) * ulog - xlog * 0.5 * cfx;
    let cfx2 = (1.0 - upw) * a.cf * xt1 + upw * b.cf * xt2;
    let dix = (1.0 - upw) * a.di * xt1 + upw * b.di * xt2;
    let r2 = hlog + (1.0 - ha) * ulog + xlog * (0.5 * cfx2 - dix);
    let dx = x2 - x1;
    let r3 = if kind == Kind::Lam {
        let axa = (0.5 * (a.ax * a.ax + b.ax * b.ax)).sqrt();
        b.c - a.c - axa * dx
    } else {
        let w = |p: f64, q: f64| (1.0 - upw) * p + upw * q;
        let sa = w(a.c, b.c);
        let cqa = w(a.cq, b.cq);
        let cfa = w(a.cf, b.cf);
        let hkw = w(a.hk, b.hk);
        let usa = 0.5 * (a.us + b.us);
        let dea = 0.5 * (a.de + b.de);
        let da = 0.5 * (a.ds + b.ds);
        let ald = if kind == Kind::Wake { DLCON } else { 1.0 };
        let hkc = if kind == Kind::Turb {
            (hkw - 1.0 - GCCON / rta).max(0.01)
        } else {
            hkw - 1.0
        };
        let hr = hkc / (GACON * ald * hkw);
        let uq = (0.5 * cfa - hr * hr) / (GBCON * da);
        let scc = SCCON * 1.333 / (1.0 + usa);
        let slog = (b.c / a.c).ln();
        scc * (cqa - sa * ald) * dx - dea * 2.0 * slog + dea * 2.0 * (uq * dx - ulog)
    };
    [r1, r2, r3]
}

// ------------------------------------------------------------ assemblage

#[derive(Clone, Copy, Debug)]
enum Group {
    /// Première station d'une face : similitude de Hiemenz.
    Simi(usize),
    Lam(usize, usize),
    Turb(usize, usize),
    /// Transition dans l'intervalle (station 1 laminaire, 2 turbulente).
    Trans(usize, usize, f64),
    /// Première station du sillage : jonction des deux faces au bord de fuite.
    Te {
        lo: usize,
        up: usize,
        w: usize,
        lo_turb: bool,
        up_turb: bool,
    },
    Wake(usize, usize),
}

impl Group {
    fn deps(&self) -> Vec<usize> {
        match *self {
            Group::Simi(p) => vec![p],
            Group::Lam(a, b) | Group::Turb(a, b) | Group::Trans(a, b, _) | Group::Wake(a, b) => vec![a, b],
            Group::Te { lo, up, w, .. } => vec![lo, up, w],
        }
    }
}

struct Ctx {
    re: f64,
    n_crit: f64,
    /// Abscisse curviligne depuis le point d'arrêt, et x/c, aux stations.
    xi: Vec<f64>,
    xs: Vec<f64>,
}

impl Ctx {
    /// v : [θ, m, c, ue, ξ] des stations de `g.deps()`.
    fn residual(&self, g: &Group, v: &[[f64; 5]]) -> [f64; 3] {
        let re = self.re;
        let cl = |kind, s: &[f64; 5]| closure(kind, s[0], s[1] / s[3].max(1e-5), s[2], s[3], re);
        match *g {
            Group::Simi(_) => {
                let k = cl(Kind::Lam, &v[0]);
                let x = v[0][4];
                let cfx = x * k.cf / k.th;
                [k.h + 2.0 - 0.5 * cfx, 1.0 - k.h + 0.5 * cfx - x * k.di / k.th, k.c]
            }
            Group::Lam(..) => interval(Kind::Lam, v[0][4], v[1][4], &cl(Kind::Lam, &v[0]), &cl(Kind::Lam, &v[1])),
            Group::Turb(..) => interval(Kind::Turb, v[0][4], v[1][4], &cl(Kind::Turb, &v[0]), &cl(Kind::Turb, &v[1])),
            Group::Wake(..) => interval(Kind::Wake, v[0][4], v[1][4], &cl(Kind::Wake, &v[0]), &cl(Kind::Wake, &v[1])),
            Group::Trans(p, q, trip) => {
                let (x1, x2) = (v[0][4], v[1][4]);
                let a = cl(Kind::Lam, &v[0]);
                let bl = cl(Kind::Lam, &v[1]);
                let (w, _) = transition_weight(
                    &a,
                    [bl.th, bl.ds, bl.ue],
                    x1,
                    x2,
                    self.xs[p],
                    self.xs[q],
                    trip,
                    re,
                    self.n_crit,
                );
                let xt = x1 + w * (x2 - x1);
                let lerp = |i: usize| v[0][i] + w * (v[1][i] - v[0][i]);
                let (tt, ut) = (lerp(0), lerp(3));
                let dt = a.ds + w * (bl.ds - a.ds);
                let tl = closure(Kind::Lam, tt, dt, self.n_crit, ut, re);
                let st = s_transition(tt, dt, ut, re);
                let tu = closure(Kind::Turb, tt, dt, st, ut, re);
                let b = cl(Kind::Turb, &v[1]);
                let l = interval(Kind::Lam, x1, xt, &a, &tl);
                let t = interval(Kind::Turb, xt, x2, &tu, &b);
                [l[0] + t[0], l[1] + t[1], t[2]]
            }
            Group::Te { lo_turb, up_turb, .. } => {
                let side = |s: &[f64; 5], turb: bool| {
                    let ds = s[1] / s[3].max(1e-5);
                    let c = if turb { s[2] } else { s_transition(s[0], ds, s[3], re) };
                    (s[0], ds, c, s[3].max(1e-5))
                };
                let (tl, dl, sl, ul) = side(&v[0], lo_turb);
                let (tu, du, su, uu) = side(&v[1], up_turb);
                let tv = tl + tu;
                let a = closure(Kind::Wake, tv, dl + du, (sl * tl + su * tu) / tv, 0.5 * (ul + uu), re);
                let b = cl(Kind::Wake, &v[2]);
                let xv = 0.5 * (v[0][4] + v[1][4]);
                interval(Kind::Wake, xv, v[2][4], &a, &b)
            }
        }
    }
}

/// Newton local à pas limité (marche initiale). `pos` : inconnues positives.
fn local_newton<F: Fn(&[f64]) -> Vec<f64>>(x: &mut [f64], floor: &[f64], pos: &[bool], f: F) -> bool {
    let n = x.len();
    for _ in 0..50 {
        let r = f(x);
        if r.iter().any(|v| !v.is_finite()) {
            return false;
        }
        let mut jac = vec![0.0; n * n];
        for j in 0..n {
            let h = 1e-6 * x[j].abs().max(floor[j]);
            let mut xp = x.to_vec();
            let mut xm = x.to_vec();
            xp[j] += h;
            xm[j] -= h;
            let (rp, rm) = (f(&xp), f(&xm));
            for i in 0..n {
                jac[i * n + j] = (rp[i] - rm[i]) / (2.0 * h);
            }
        }
        let neg: Vec<f64> = r.iter().map(|v| -v).collect();
        let dx = Lu::new(jac, n).solve(&neg);
        if dx.iter().any(|v| !v.is_finite()) {
            return false;
        }
        let mut rlx: f64 = 1.0;
        for j in 0..n {
            if pos[j] {
                let q = dx[j] / x[j];
                if q < -0.5 {
                    rlx = rlx.min(-0.5 / q);
                }
                if q > 1.0 {
                    rlx = rlx.min(1.0 / q);
                }
            }
        }
        let mut step: f64 = 0.0;
        for j in 0..n {
            x[j] += rlx * dx[j];
            step = step.max((rlx * dx[j]).abs() / x[j].abs().max(floor[j]));
        }
        if step < 1e-9 {
            return true;
        }
    }
    false
}

/// Résultat du couplage simultané.
pub struct Coupled {
    pub cl: f64,
    pub cm: f64,
    pub cd: f64,
    /// Traînée de frottement (intégrale de Cf sur le profil).
    pub cdf: f64,
    pub cl_inviscid: f64,
    pub transition_upper: f64,
    pub transition_lower: f64,
    pub separation_upper: f64,
    pub separation_lower: f64,
    pub friction_upper: Friction,
    pub friction_lower: Friction,
    pub cp_x: Vec<f64>,
    pub cp_y: Vec<f64>,
    pub cp: Vec<f64>,
    pub i_le: usize,
    pub converged: bool,
    pub iterations: usize,
}

struct Solver<'a> {
    md: &'a Model,
    o: PanelOptions,
    ist: usize,
    cs: f64,
    /// Nombre de stations laminaires de chaque face (indice de la station de transition).
    itr: [usize; 2],
    th: Vec<f64>,
    m: Vec<f64>,
    c: Vec<f64>,
    /// ue courante (XFOIL : l'écart à ue_inv + D·m est résorbé par le pas de Newton).
    ue: Vec<f64>,
}

impl<'a> Solver<'a> {
    fn sides(&self) -> [Vec<usize>; 2] {
        let n = self.md.n;
        [(0..=self.ist).rev().collect(), (self.ist + 1..n).collect()]
    }

    fn sign(&self, p: usize) -> f64 {
        if p <= self.ist {
            -1.0
        } else {
            1.0
        }
    }

    fn trip(&self, side: usize) -> f64 {
        let t = if side == 0 { self.o.x_trip_lower } else { self.o.x_trip_upper };
        t.unwrap_or(f64::NAN)
    }

    fn ctx(&self) -> Ctx {
        let md = self.md;
        let n = md.n;
        let mut xi = vec![0.0; md.ns()];
        let mut xs = vec![0.0; md.ns()];
        for p in 0..n {
            xi[p] = (md.c_ctrl[p] - self.cs).abs().max(1e-7);
            xs[p] = md.xc[p];
        }
        let xv = 0.5 * (xi[0] + xi[n - 1]);
        let half = 0.25 * (md.len[0] + md.len[n - 1]);
        for k in 0..NW {
            xi[n + k] = xv + half + md.wd[k];
            xs[n + k] = md.wx[k];
        }
        Ctx {
            re: self.o.re,
            n_crit: self.o.n_crit,
            xi,
            xs,
        }
    }

    fn kinds(&self) -> Vec<Kind> {
        let mut k = vec![Kind::Wake; self.md.ns()];
        for (s, side) in self.sides().iter().enumerate() {
            for (i, &p) in side.iter().enumerate() {
                k[p] = if i < self.itr[s] { Kind::Lam } else { Kind::Turb };
            }
        }
        k
    }

    fn groups(&self) -> Vec<Group> {
        let n = self.md.n;
        let mut g = vec![Group::Simi(0); self.md.ns()];
        for (s, side) in self.sides().iter().enumerate() {
            for (i, &p) in side.iter().enumerate() {
                g[p] = if i == 0 {
                    Group::Simi(p)
                } else if i < self.itr[s] {
                    Group::Lam(side[i - 1], p)
                } else if i == self.itr[s] {
                    Group::Trans(side[i - 1], p, self.trip(s))
                } else {
                    Group::Turb(side[i - 1], p)
                };
            }
        }
        let [lo_side, up_side] = self.sides();
        g[n] = Group::Te {
            lo: 0,
            up: n - 1,
            w: n,
            lo_turb: lo_side.len() > self.itr[0],
            up_turb: up_side.len() > self.itr[1],
        };
        for k in 1..NW {
            g[n + k] = Group::Wake(n + k - 1, n + k);
        }
        g
    }

    /// Point d'arrêt : changement de signe de q le plus proche du bord d'attaque.
    fn locate_stagnation(&self, q: &[f64]) -> Option<(usize, f64)> {
        let md = self.md;
        let ile = leading_panel(&md.xc);
        let mut best: Option<(usize, f64)> = None;
        for i in 0..md.n - 1 {
            if q[i] <= 0.0 && q[i + 1] > 0.0 && best.is_none_or(|(b, _)| i.abs_diff(ile) < b.abs_diff(ile)) {
                let f = q[i] / (q[i] - q[i + 1]);
                best = Some((i, md.c_ctrl[i] + f * (md.c_ctrl[i + 1] - md.c_ctrl[i])));
            }
        }
        best
    }

    /// Marche directe (inverse au décollement) le long des faces puis du sillage,
    /// à vitesse non visqueuse : point de départ du Newton.
    fn march(&mut self) -> bool {
        let md = self.md;
        let n = md.n;
        let ctx = self.ctx();
        let re = self.o.re;
        let ue_inv: Vec<f64> = (0..md.ns()).map(|p| (self.sign(p) * md.q_inv[p]).max(1e-4)).collect();
        let mut ue = ue_inv.clone();
        let mut ds = vec![0.0; md.ns()];
        let sides = self.sides();
        for (s, side) in sides.iter().enumerate() {
            let p0 = side[0];
            let u0 = ue[p0];
            let x0 = ctx.xi[p0];
            let mut x = [(0.085 * x0 / (re * u0)).sqrt(), 0.0];
            x[1] = 2.2 * x[0];
            let ok = local_newton(&mut x, &[1e-9, 1e-9], &[true, true], |x| {
                let r = ctx.residual(&Group::Simi(p0), &[[x[0], u0 * x[1], 0.0, u0, x0]]);
                vec![r[0], r[1]]
            });
            if !ok {
                return false;
            }
            self.th[p0] = x[0];
            ds[p0] = x[1];
            self.c[p0] = 0.0;
            self.itr[s] = side.len();
            let mut turb = false;
            for i in 1..side.len() {
                let (pp, p) = (side[i - 1], side[i]);
                let mut g = if turb { Group::Turb(pp, p) } else { Group::Lam(pp, p) };
                let mut c0 = self.c[pp];
                let solved = self.march_station(&ctx, g, pp, p, c0, &mut ue, &mut ds, &ue_inv);
                if !turb {
                    let forced = self.trip(s).is_finite() && md.xc[p] >= self.trip(s);
                    if !solved || self.c[p] >= self.o.n_crit || forced {
                        turb = true;
                        self.itr[s] = i;
                        g = Group::Trans(pp, p, self.trip(s));
                        c0 = s_transition(self.th[pp], ds[pp], ue[pp], re);
                        self.march_station(&ctx, g, pp, p, c0, &mut ue, &mut ds, &ue_inv);
                    }
                }
            }
        }
        // Sillage.
        let kinds = self.kinds();
        for k in 0..NW {
            let p = n + k;
            let g = if k == 0 {
                Group::Te {
                    lo: 0,
                    up: n - 1,
                    w: n,
                    lo_turb: kinds[0] == Kind::Turb,
                    up_turb: kinds[n - 1] == Kind::Turb,
                }
            } else {
                Group::Wake(p - 1, p)
            };
            if k == 0 {
                let side_s = |q: usize| {
                    if kinds[q] == Kind::Turb {
                        self.c[q]
                    } else {
                        s_transition(self.th[q], ds[q], ue[q], re)
                    }
                };
                let sv = (side_s(0) * self.th[0] + side_s(n - 1) * self.th[n - 1]) / (self.th[0] + self.th[n - 1]);
                self.th[p] = self.th[0] + self.th[n - 1];
                ds[p] = ds[0] + ds[n - 1];
                self.c[p] = sv;
            } else {
                self.th[p] = self.th[p - 1];
                ds[p] = ds[p - 1];
                self.c[p] = self.c[p - 1];
            }
            let c0 = self.c[p];
            self.march_station(&ctx, g, if k == 0 { 0 } else { p - 1 }, p, c0, &mut ue, &mut ds, &ue_inv);
        }
        for p in 0..md.ns() {
            self.m[p] = ue[p] * ds[p];
            self.ue[p] = ue[p];
        }
        true
    }

    /// Une station de la marche : direct, ou inverse (Hk imposé, ue inconnue)
    /// si Hk dépasse 3,8 (laminaire) ou 2,5 (turbulent).
    #[allow(clippy::too_many_arguments)]
    fn march_station(
        &mut self,
        ctx: &Ctx,
        g: Group,
        pp: usize,
        p: usize,
        c0: f64,
        ue: &mut [f64],
        ds: &mut [f64],
        ue_inv: &[f64],
    ) -> bool {
        let deps = g.deps();
        let kind = match g {
            Group::Lam(..) => Kind::Lam,
            Group::Wake(..) | Group::Te { .. } => Kind::Wake,
            _ => Kind::Turb,
        };
        let state = |q: usize, th: &[f64], c: &[f64], ue: &[f64], ds: &[f64]| [th[q], ue[q] * ds[q], c[q], ue[q]];
        let (th0, ds0) = if matches!(g, Group::Te { .. }) {
            (self.th[p], ds[p])
        } else {
            (self.th[pp], ds[pp])
        };
        let u = ue_inv[p];
        let cpos = kind != Kind::Lam;
        let cfloor = if cpos { 1e-4 } else { 1e-3 };
        let (thv, cv, uev, dsv) = (self.th.clone(), self.c.clone(), ue.to_vec(), ds.to_vec());
        let eval = |x: &[f64], ue_p: f64| {
            let v: Vec<[f64; 5]> = deps
                .iter()
                .map(|&q| {
                    if q == p {
                        [x[0], ue_p * x[1], x[2], ue_p, ctx.xi[q]]
                    } else {
                        let s = state(q, &thv, &cv, &uev, &dsv);
                        [s[0], s[1], s[2], s[3], ctx.xi[q]]
                    }
                })
                .collect();
            ctx.residual(&g, &v)
        };
        let mut x = vec![th0, ds0, c0];
        let ok = local_newton(&mut x, &[1e-9, 1e-9, cfloor], &[true, true, cpos], |x| eval(x, u).to_vec());
        let hmax = if kind == Kind::Lam { 3.8 } else { 2.5 };
        let mut res = (x[0], x[1], x[2], u);
        if !ok || x[1] / x[0] > hmax || !x.iter().all(|v| v.is_finite()) {
            // Inverse : Hk imposé.
            let hk1 = ds0 / th0;
            let dx = ctx.xi[p] - ctx.xi[pp].min(ctx.xi[p]);
            let htarg = match kind {
                Kind::Lam => (hk1 + 0.03 * dx / th0).max(hmax),
                Kind::Turb => (hk1 - 0.15 * dx / th0).max(hmax),
                Kind::Wake => (hk1 - 0.03 * dx / th0).max(1.01),
            };
            let up = if p > 0 && pp < ue.len() { ue[pp] } else { u };
            let mut y = vec![th0, htarg * th0, c0, up];
            let ok2 = local_newton(&mut y, &[1e-9, 1e-9, cfloor, 1e-3], &[true, true, cpos, true], |y| {
                let r = eval(y, y[3]);
                vec![r[0], r[1], r[2], y[1] / y[0] - htarg]
            });
            if ok2 && y.iter().all(|v| v.is_finite()) {
                res = (y[0], y[1], y[2], y[3]);
            } else if !ok {
                // Repli : on recopie la station précédente.
                res = (th0, ds0, c0, ue[pp.min(ue.len() - 1)]);
            }
        }
        self.th[p] = res.0;
        ds[p] = res.1;
        self.c[p] = res.2;
        ue[p] = res.3;
        ok
    }

    /// Vitesses signées aux stations et ue.
    fn speeds(&self, d: &[f64]) -> Vec<f64> {
        let ns = self.md.ns();
        (0..ns)
            .map(|i| self.md.q_inv[i] + (0..ns).map(|k| d[i * ns + k] * self.m[k]).sum::<f64>())
            .collect()
    }

    fn newton(&mut self) -> (bool, usize) {
        let md = self.md;
        let n = md.n;
        let ns = md.ns();
        let mut kinds_prev = self.kinds();
        for it in 0..MAX_IT {
            // Point d'arrêt.
            let d;
            let qs: Vec<f64> = (0..ns).map(|p| self.sign(p) * self.ue[p]).collect();
            if let Some((ist, cs)) = self.locate_stagnation(&qs) {
                if ist != self.ist {
                    let shift = ist as isize - self.ist as isize;
                    let (a, b) = (ist.min(self.ist), ist.max(self.ist));
                    for p in a..=b {
                        self.ue[p] = self.ue[p].abs();
                    }
                    self.itr[0] = (self.itr[0] as isize + shift).max(1) as usize;
                    self.itr[1] = (self.itr[1] as isize - shift).max(1) as usize;
                    self.ist = ist;
                }
                self.cs = cs;
            }
            d = md.mass_matrix(self.ist, self.cs);
            let ue = self.ue.clone();
            let ctx = self.ctx();
            // Transition : remonte librement, descend d'une station par itération.
            let sides = self.sides();
            let mut n_lam = vec![0.0; ns];
            for s in 0..2 {
                let side = &sides[s];
                let old = self.itr[s].min(side.len());
                let mut new = None;
                for k in 1..old.min(side.len()) {
                    let (pp, p) = (side[k - 1], side[k]);
                    let a = closure(Kind::Lam, self.th[pp], self.m[pp] / ue[pp].max(1e-5), self.c[pp], ue[pp], self.o.re);
                    let b = closure(Kind::Lam, self.th[p], self.m[p] / ue[p].max(1e-5), 0.0, ue[p], self.o.re);
                    let (_, found) = transition_weight(
                        &a,
                        [b.th, b.ds, b.ue],
                        ctx.xi[pp],
                        ctx.xi[p],
                        md.xc[pp],
                        md.xc[p],
                        self.trip(s),
                        self.o.re,
                        self.o.n_crit,
                    );
                    if found {
                        new = Some(k);
                        break;
                    }
                }
                self.itr[s] = match new {
                    // Une station par itération au plus, pour éviter les va-et-vient.
                    Some(k) => k.max(old.saturating_sub(1)).max(1),
                    None if old < side.len() && old > 0 => {
                        // Pas de transition jusqu'à `old` : marche laminaire à ue courante
                        // jusqu'à ce que N atteigne N_crit (TRCHEK de proche en proche).
                        let mut uea = self.ue.clone();
                        let mut dsa: Vec<f64> = (0..ns).map(|q| self.m[q] / self.ue[q].max(1e-5)).collect();
                        let mut k = old;
                        while k < side.len() {
                            let (pp, p) = (side[k - 1], side[k]);
                            let ue0 = uea.clone();
                            let saved = (self.th[p], self.c[p], uea[p], dsa[p]);
                            self.march_station(&ctx, Group::Lam(pp, p), pp, p, self.c[pp], &mut uea, &mut dsa, &ue0);
                            let a = closure(Kind::Lam, self.th[pp], dsa[pp], self.c[pp], uea[pp], self.o.re);
                            let (_, found) = transition_weight(
                                &a,
                                [self.th[p], dsa[p], uea[p]],
                                ctx.xi[pp],
                                ctx.xi[p],
                                md.xc[pp],
                                md.xc[p],
                                self.trip(s),
                                self.o.re,
                                self.o.n_crit,
                            );
                            if found || self.c[p] >= self.o.n_crit || !self.th[p].is_finite() {
                                if k == old || !self.th[p].is_finite() {
                                    // La transition reste là : état turbulent conservé.
                                    self.th[p] = saved.0;
                                    self.c[p] = saved.1;
                                } else {
                                    // Nouvelle station de transition : départ laminaire marché.
                                    self.m[p] = uea[p] * dsa[p];
                                    self.ue[p] = uea[p];
                                    self.c[p] = s_transition(self.th[p], dsa[p], uea[p], self.o.re);
                                }
                                break;
                            }
                            n_lam[p] = self.c[p];
                            self.m[p] = uea[p] * dsa[p];
                            self.ue[p] = uea[p];
                            k += 1;
                            if k > old {
                                break;
                            }
                        }
                        k
                    }
                    None => old,
                };
            }
            // Variables dont la nature a changé.
            let kinds = self.kinds();
            for p in 0..ns {
                if kinds[p] != kinds_prev[p] {
                    let ds = self.m[p] / ue[p].max(1e-5);
                    self.c[p] = if kinds[p] == Kind::Lam {
                        n_lam[p]
                    } else {
                        s_transition(self.th[p], ds, ue[p], self.o.re)
                    };
                }
            }
            kinds_prev = kinds.clone();
            // La marche de transition a pu changer m et ue.
            let q = self.speeds(&d);
            let gap: Vec<f64> = (0..ns).map(|p| self.sign(p) * q[p] - self.ue[p]).collect();
            let ue = self.ue.clone();
            let groups = self.groups();
            // Résidus et jacobienne.
            let nv = 3 * ns;
            let mut jac = vec![0.0; nv * nv];
            let mut r = vec![0.0; nv];
            let mut r_ue = vec![0.0; nv];
            // Sensibilité du point d'arrêt aux vitesses des deux stations qui l'encadrent.
            let (ia, ib) = (self.ist, self.ist + 1);
            let (ua, ub) = (ue[ia].max(1e-6), ue[ib].max(1e-6));
            let dc = md.c_ctrl[ib] - md.c_ctrl[ia];
            let dcs_a = dc * ub / ((ua + ub) * (ua + ub));
            let dcs_b = -dc * ua / ((ua + ub) * (ua + ub));
            for p in 0..ns {
                let g = &groups[p];
                let deps = g.deps();
                let v: Vec<[f64; 5]> =
                    deps.iter().map(|&q| [self.th[q], self.m[q], self.c[q], ue[q], ctx.xi[q]]).collect();
                let r0 = ctx.residual(g, &v);
                r[3 * p..3 * p + 3].copy_from_slice(&r0);
                for (di, &qd) in deps.iter().enumerate() {
                    for var in 0..5 {
                        // ξ ne dépend du point d'arrêt que sur le profil.
                        if var == 4 && qd >= n {
                            continue;
                        }
                        let floor = [1e-8, 1e-8, 1e-3, 1e-3, 1e-6][var];
                        let h = 1e-6 * v[di][var].abs().max(floor);
                        let mut vp = v.clone();
                        let mut vm = v.clone();
                        vp[di][var] += h;
                        vm[di][var] -= h;
                        let (rp, rm) = (ctx.residual(g, &vp), ctx.residual(g, &vm));
                        for e in 0..3 {
                            let dr = (rp[e] - rm[e]) / (2.0 * h);
                            if dr == 0.0 {
                                continue;
                            }
                            let row = (3 * p + e) * nv;
                            if var < 3 {
                                jac[row + 3 * qd + var] += dr;
                            } else if var == 4 {
                                // ξ = ±(c − c_s), c_s interpolé entre ue[ist] et ue[ist+1].
                                let sx = if qd <= self.ist { 1.0 } else { -1.0 };
                                for (&sta, &dcs) in [ia, ib].iter().zip(&[dcs_a, dcs_b]) {
                                    let f = sx * dr * dcs;
                                    r_ue[3 * p + e] += f * gap[sta];
                                    let sg = self.sign(sta) * f;
                                    for k in 0..ns {
                                        jac[row + 3 * k + 1] += sg * d[sta * ns + k];
                                    }
                                }
                            } else {
                                r_ue[3 * p + e] += dr * gap[qd];
                                let sg = self.sign(qd) * dr;
                                for k in 0..ns {
                                    jac[row + 3 * k + 1] += sg * d[qd * ns + k];
                                }
                            }
                        }
                    }
                }
            }
            if r.iter().any(|v| !v.is_finite()) {
                return (false, it);
            }
            let neg: Vec<f64> = r.iter().zip(&r_ue).map(|(v, w)| -v - w).collect();
            let dx = Lu::new(jac, nv).solve(&neg);
            if dx.iter().any(|v| !v.is_finite()) {
                return (false, it);
            }
            // Pas limité (UPDATE de XFOIL).
            let dm: Vec<f64> = (0..ns).map(|p| dx[3 * p + 1]).collect();
            let mut rlx: f64 = 1.0;
            let mut lim = |dn: f64| {
                if dn > 1.5 {
                    rlx = rlx.min(1.5 / dn);
                }
                if dn < -0.5 {
                    rlx = rlx.min(-0.5 / dn);
                }
            };
            let mut rms = 0.0;
            let mut dues = vec![0.0; ns];
            for p in 0..ns {
                let due = gap[p] + self.sign(p) * (0..ns).map(|k| d[p * ns + k] * dm[k]).sum::<f64>();
                dues[p] = due;
                let u = ue[p].max(1e-5);
                let ds = self.m[p] / u;
                let dds = (dm[p] - ds * due) / u;
                let dn1 = if kinds[p] == Kind::Lam {
                    dx[3 * p + 2] / 10.0
                } else {
                    dx[3 * p + 2] / self.c[p]
                };
                let dn2 = dx[3 * p] / self.th[p];
                let dn3 = dds / ds;
                let dn4 = due.abs() / 0.25;
                for dn in [dn1, dn2, dn3, dn4] {
                    lim(dn);
                    rms += dn * dn;
                }
            }
            let rms = (rms / (4 * ns) as f64).sqrt();
            for p in 0..ns {
                self.th[p] = (self.th[p] + rlx * dx[3 * p]).max(1e-10);
                self.m[p] += rlx * dx[3 * p + 1];
                self.c[p] += rlx * dx[3 * p + 2];
                self.ue[p] += rlx * dues[p];
                if kinds[p] != Kind::Lam {
                    self.c[p] = self.c[p].max(1e-7);
                }
            }
            // Hk minimal.
            for p in 0..ns {
                let u = self.ue[p].max(1e-5);
                let hmin = if p >= n { 1.00005 } else { 1.02 };
                let ds = (self.m[p] / u).max(hmin * self.th[p]);
                self.m[p] = u * ds;
            }
            if rlx * rms < 1e-6 {
                return (true, it + 1);
            }
        }
        (false, MAX_IT)
    }
}

/// État de la couche limite (transportable d'une incidence à l'autre : mêmes panneaux).
#[derive(Clone)]
struct BlState {
    ist: usize,
    cs: f64,
    itr: [usize; 2],
    th: Vec<f64>,
    m: Vec<f64>,
    c: Vec<f64>,
    ue: Vec<f64>,
}

impl<'a> Solver<'a> {
    fn from_state(md: &'a Model, o: &PanelOptions, s: BlState) -> Solver<'a> {
        Solver {
            md,
            o: *o,
            ist: s.ist,
            cs: s.cs,
            itr: s.itr,
            th: s.th,
            m: s.m,
            c: s.c,
            ue: s.ue,
        }
    }

    fn state(&self) -> BlState {
        BlState {
            ist: self.ist,
            cs: self.cs,
            itr: self.itr,
            th: self.th.clone(),
            m: self.m.clone(),
            c: self.c.clone(),
            ue: self.ue.clone(),
        }
    }
}

/// Une résolution : marche puis Newton, ou Newton depuis un état donné.
fn attempt(md: &Model, o: &PanelOptions, init: Option<&BlState>) -> (BlState, bool, usize) {
    let ns = md.ns();
    let mut sv = match init {
        Some(s) => Solver::from_state(md, o, s.clone()),
        None => {
            let mut sv = Solver {
                md,
                o: *o,
                ist: 0,
                cs: 0.0,
                itr: [0, 0],
                th: vec![0.0; ns],
                m: vec![0.0; ns],
                c: vec![0.0; ns],
                ue: vec![0.0; ns],
            };
            let ile = leading_panel(&md.xc);
            let (ist, cs) = sv.locate_stagnation(&md.q_inv).unwrap_or((ile, md.c_ctrl[ile]));
            sv.ist = ist;
            sv.cs = cs;
            if !sv.march() {
                return (sv.state(), false, 0);
            }
            sv
        }
    };
    let (c, i) = sv.newton();
    (sv.state(), c, i)
}

/// Couplage simultané pour une incidence (degrés).
pub fn solve_coupled(shape: &PanelShape, alpha_deg: f64, o: &PanelOptions) -> Coupled {
    let alpha = alpha_deg.to_radians();
    let pts = airfoil_points(shape, 80);
    let md = Model::new(pts.clone(), alpha);
    let n = md.n;
    let q0 = md.q_inv.clone();
    let (mut st, mut converged, mut iterations) = attempt(&md, o, None);
    if !converged && alpha_deg.abs() > 1.0 {
        // Continuation en incidence depuis 0°, par pas d'au plus 2°.
        let steps = (alpha_deg.abs() / 2.0).ceil() as usize;
        let mut prev: Option<BlState> = None;
        let mut ok = true;
        for k in 0..steps {
            let ak = alpha_deg * k as f64 / steps as f64;
            let mk = Model::new(pts.clone(), ak.to_radians());
            let (s, c, i) = attempt(&mk, o, prev.as_ref());
            iterations += i;
            if !c {
                ok = false;
                break;
            }
            prev = Some(s);
        }
        if ok {
            let (s, c, i) = attempt(&md, o, prev.as_ref());
            iterations += i;
            if c {
                st = s;
                converged = true;
            }
        }
    }
    let sv = Solver::from_state(&md, o, st);
    let ns = md.ns();
    // Efforts.
    let d = md.mass_matrix(sv.ist, sv.cs);
    let q = sv.speeds(&d);
    let integrate = |q: &[f64]| {
        let (mut cx, mut cy, mut cm) = (0.0, 0.0, 0.0);
        let mut cp = vec![0.0; n];
        for i in 0..n {
            cp[i] = 1.0 - q[i] * q[i];
            let (nx, ny) = (-md.th[i].sin(), md.th[i].cos());
            let fx = -cp[i] * md.len[i] * nx;
            let fy = -cp[i] * md.len[i] * ny;
            cx += fx;
            cy += fy;
            cm += -(md.xc[i] - 0.25) * fy + md.yc[i] * fx;
        }
        (cy * alpha.cos() - cx * alpha.sin(), cm, cp)
    };
    let (cl, cm, cp) = integrate(&q);
    let (cl_inviscid, _, _) = integrate(&q0);
    let ue: Vec<f64> = (0..ns).map(|p| (sv.sign(p) * q[p]).max(1e-5)).collect();
    let kinds = sv.kinds();
    let cl_at = |p: usize| closure(kinds[p], sv.th[p], sv.m[p] / ue[p], sv.c[p], ue[p], o.re);
    let last = n + NW - 1;
    let wl = cl_at(last);
    let cd = 2.0 * wl.th * wl.ue.powf(0.5 * (wl.h + 5.0));
    let sides = sv.sides();
    let ctx = sv.ctx();
    let mut cdf = 0.0;
    let mut out = Vec::new();
    for (s, side) in sides.iter().enumerate() {
        let mut fr = Friction::default();
        let mut sep = f64::NAN;
        for &p in side {
            let k = cl_at(p);
            fr.x.push(md.xc[p]);
            fr.cf.push(k.cf);
            fr.ue.push(k.ue);
            if k.cf < 0.0 && kinds[p] == Kind::Turb && !sep.is_finite() {
                sep = md.xc[p];
            }
            // Projection du frottement sur la direction de l'écoulement libre.
            let t = sv.sign(p) * (md.th[p] - alpha).cos();
            cdf += k.cf * k.ue * k.ue * md.len[p] * t;
        }
        // Point de transition : interpolé dans l'intervalle (comme dans Trans).
        let itr = sv.itr[s];
        let xtr = if itr < side.len() {
            let (pp, p) = (side[itr - 1], side[itr]);
            let a = cl_at(pp);
            let bl = closure(Kind::Lam, sv.th[p], sv.m[p] / ue[p], 0.0, ue[p], o.re);
            let (w, _) = transition_weight(
                &a,
                [bl.th, bl.ds, bl.ue],
                ctx.xi[pp],
                ctx.xi[p],
                md.xc[pp],
                md.xc[p],
                sv.trip(s),
                o.re,
                o.n_crit,
            );
            md.xc[pp] + w * (md.xc[p] - md.xc[pp])
        } else {
            1.0
        };
        out.push((fr, sep, xtr));
    }
    let (fu, su, xu) = out.pop().unwrap();
    let (fl, sl, xl) = out.pop().unwrap();
    Coupled {
        cl,
        cm,
        cd,
        cdf,
        cl_inviscid,
        transition_upper: xu,
        transition_lower: xl,
        separation_upper: su,
        separation_lower: sl,
        friction_upper: fu,
        friction_lower: fl,
        cp_x: md.xc.clone(),
        cp_y: md.yc.clone(),
        cp,
        i_le: leading_panel(&md.xc),
        converged,
        iterations,
    }
}


impl PanelResult {
    /// Résultat du couplage simultané au format commun.
    pub(crate) fn from_coupled(alpha_deg: f64, c: Coupled) -> PanelResult {
        PanelResult {
            alpha_deg,
            cl: if c.converged { c.cl } else { f64::NAN },
            cd: if c.converged { c.cd } else { f64::NAN },
            cm: c.cm,
            cl_inviscid: c.cl_inviscid,
            transition_upper: c.transition_upper,
            transition_lower: c.transition_lower,
            separation_upper: c.separation_upper,
            separation_lower: c.separation_lower,
            separated: false,
            friction_upper: c.friction_upper,
            friction_lower: c.friction_lower,
            cp_x: c.cp_x,
            cp_y: c.cp_y,
            cp: c.cp,
            i_le: c.i_le,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Une bosse fermée de δ* sur une face (sans sillage) : la transpiration
    /// doit rendre la portance du corps épaissi géométriquement.
    #[test]
    fn transpiration_contre_corps_epaissi() {
        let shape = PanelShape::Naca(crate::shapes::Naca::parse("0012").unwrap());
        let al = 4f64.to_radians();
        let md = Model::new(airfoil_points(&shape, 80), al);
        let n = md.n;
        let ns = md.ns();
        let ist = (0..n - 1).find(|&i| md.q_inv[i] <= 0.0 && md.q_inv[i + 1] > 0.0).unwrap();
        let f = md.q_inv[ist] / (md.q_inv[ist] - md.q_inv[ist + 1]);
        let cs = md.c_ctrl[ist] + f * (md.c_ctrl[ist + 1] - md.c_ctrl[ist]);
        let d = md.mass_matrix(ist, cs);
        let bump = |x: f64| if x > 0.6 { 0.01 * (PI * (x - 0.6) / 0.4).sin() } else { 0.0 };
        let cl = |q: &[f64]| {
            let (mut cx, mut cy) = (0.0, 0.0);
            for i in 0..n {
                let cp = 1.0 - q[i] * q[i];
                cx += cp * md.len[i] * md.th[i].sin();
                cy += -cp * md.len[i] * md.th[i].cos();
            }
            cy * al.cos() - cx * al.sin()
        };
        for lower in [true, false] {
            let m: Vec<f64> = (0..ns)
                .map(|p| if p < n && (p <= ist) == lower { bump(md.xc[p]) * md.q_inv[p].abs() } else { 0.0 })
                .collect();
            let q: Vec<f64> = (0..ns).map(|i| md.q_inv[i] + (0..ns).map(|k| d[i * ns + k] * m[k]).sum::<f64>()).collect();
            let pts: Vec<[f64; 2]> = md
                .pts
                .iter()
                .enumerate()
                .map(|(k, &[x, y])| {
                    let upper = k > n / 2;
                    let dd = if upper != lower { bump(x) } else { 0.0 };
                    [x, if upper { y + dd } else { y - dd }]
                })
                .collect();
            let g = crate::panel::potential_flow(&pts, al).cl;
            println!("bosse {} : transpiration {:.4}, corps épaissi {:.4}", if lower { "intrados" } else { "extrados" }, cl(&q), g);
            assert!((cl(&q) - g).abs() < 0.05 * (g - 0.4797).abs());
        }
    }
}
