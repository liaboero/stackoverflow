//! Méthode des panneaux et couche limite intégrale (profils minces, Re ≈ 10⁴
//! à 10⁷) : portage de `js/panel.js` de l'application web.
//!
//! - Écoulement potentiel : panneaux tourbillons à densité linéaire, condition
//!   de Kutta (Kuethe & Chow).
//! - Couche limite laminaire : Thwaites (corrélations de Cebeci–Bradshaw).
//! - Transition : méthode eᴺ enveloppe de Drela & Giles (1987), N_crit = 9 par
//!   défaut, ou au décollement laminaire, ou imposée.
//! - Couche limite turbulente : entraînement de Head, frottement de
//!   Ludwieg–Tillmann ; décollement si H > 2,4.
//! - Traînée : Squire–Young près du bord de fuite.
//! - Couplage faible : profil épaissi de δ*, une itération sous-relaxée.
//!
//! Coordonnées en cordes, x de 0 (bord d'attaque) à 1 (bord de fuite).

use crate::shapes::{Naca, ShapeKind, ShapeParams};
use std::f64::consts::PI;

/// Profil traité par la méthode des panneaux.
#[derive(Clone, Copy, Debug)]
pub enum PanelShape {
    Naca(Naca),
    /// Ellipse d'épaisseur relative donnée.
    Ellipse(f64),
    /// Plaque : ellipse de 2 % d'épaisseur.
    Plate,
}

impl PanelShape {
    /// Formes de la veine que la méthode sait traiter.
    pub fn from_shape(kind: ShapeKind, p: &ShapeParams) -> Option<PanelShape> {
        match kind {
            ShapeKind::Naca => Some(PanelShape::Naca(p.naca)),
            ShapeKind::Ellipse => Some(PanelShape::Ellipse(p.thickness)),
            ShapeKind::Plate => Some(PanelShape::Plate),
            _ => None,
        }
    }
}

/// Points du profil : bord de fuite → intrados → bord d'attaque → extrados → bord de fuite.
pub fn airfoil_points(shape: &PanelShape, n: usize) -> Vec<[f64; 2]> {
    // Demi-épaisseur et cambrure en x/c.
    let t = t_of(shape);
    let half = |x: f64| match *shape {
        PanelShape::Naca(nc) => nc.thickness(x),
        _ => t * (x * (1.0 - x)).max(0.0).sqrt(),
    };
    let camber = |x: f64| match *shape {
        PanelShape::Naca(nc) => nc.camber(x),
        _ => 0.0,
    };
    let mut upper = Vec::with_capacity(n + 1);
    let mut lower = Vec::with_capacity(n + 1);
    for k in 0..=n {
        let x = 0.5 * (1.0 - (PI * k as f64 / n as f64).cos());
        upper.push([x, camber(x) + half(x)]);
        lower.push([x, camber(x) - half(x)]);
    }
    let mut pts: Vec<[f64; 2]> = lower.iter().rev().cloned().collect();
    pts.extend_from_slice(&upper[1..]);
    pts
}

fn t_of(shape: &PanelShape) -> f64 {
    match *shape {
        PanelShape::Ellipse(t) => t.clamp(0.05, 1.0),
        _ => 0.02,
    }
}

/// Système linéaire dense (élimination de Gauss, pivot partiel).
fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for c in 0..n {
        let mut piv = c;
        for r in c + 1..n {
            if a[r][c].abs() > a[piv][c].abs() {
                piv = r;
            }
        }
        a.swap(c, piv);
        b.swap(c, piv);
        let d = a[c][c];
        for r in c + 1..n {
            let f = a[r][c] / d;
            if f == 0.0 {
                continue;
            }
            let (top, bottom) = a.split_at_mut(r);
            let (ac, ar) = (&top[c], &mut bottom[0]);
            for k in c..n {
                ar[k] -= f * ac[k];
            }
            b[r] -= f * b[c];
        }
    }
    let mut x = vec![0.0; n];
    for r in (0..n).rev() {
        let mut s = b[r];
        for k in r + 1..n {
            s -= a[r][k] * x[k];
        }
        x[r] = s / a[r][r];
    }
    x
}

/// Solution potentielle : vitesse tangentielle vt/V∞ et Cp aux points de contrôle.
#[derive(Clone, Debug)]
pub struct Flow {
    pub xc: Vec<f64>,
    pub yc: Vec<f64>,
    pub vt: Vec<f64>,
    pub cp: Vec<f64>,
    pub cl: f64,
    pub cm: f64,
}

/// Écoulement potentiel autour d'un contour fermé (incidence en radians).
pub fn potential_flow(pts: &[[f64; 2]], alpha: f64) -> Flow {
    let n = pts.len() - 1;
    let (mut xc, mut yc, mut th, mut s) = (vec![0.0; n], vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    for j in 0..n {
        let ([x0, y0], [x1, y1]) = (pts[j], pts[j + 1]);
        xc[j] = 0.5 * (x0 + x1);
        yc[j] = 0.5 * (y0 + y1);
        s[j] = (x1 - x0).hypot(y1 - y0);
        th[j] = (y1 - y0).atan2(x1 - x0);
    }
    let mut an = vec![vec![0.0; n + 1]; n];
    let mut at = vec![vec![0.0; n + 1]; n];
    let (mut cn1, mut cn2, mut ct1, mut ct2) = (vec![0.0; n], vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    for i in 0..n {
        for j in 0..n {
            if i == j {
                cn1[j] = -1.0;
                cn2[j] = 1.0;
                ct1[j] = PI / 2.0;
                ct2[j] = PI / 2.0;
                continue;
            }
            let dx = xc[i] - pts[j][0];
            let dy = yc[i] - pts[j][1];
            let a = -dx * th[j].cos() - dy * th[j].sin();
            let b = dx * dx + dy * dy;
            let c = (th[i] - th[j]).sin();
            let d = (th[i] - th[j]).cos();
            let e = dx * th[j].sin() - dy * th[j].cos();
            let f = (1.0 + s[j] * (s[j] + 2.0 * a) / b).ln();
            let g = (e * s[j]).atan2(b + a * s[j]);
            let p = dx * (th[i] - 2.0 * th[j]).sin() + dy * (th[i] - 2.0 * th[j]).cos();
            let q = dx * (th[i] - 2.0 * th[j]).cos() - dy * (th[i] - 2.0 * th[j]).sin();
            cn2[j] = d + 0.5 * q * f / s[j] - (a * c + d * e) * g / s[j];
            cn1[j] = 0.5 * d * f + c * g - cn2[j];
            ct2[j] = c + 0.5 * p * f / s[j] + (a * d - c * e) * g / s[j];
            ct1[j] = 0.5 * c * f - d * g - ct2[j];
        }
        an[i][0] = cn1[0];
        at[i][0] = ct1[0];
        for j in 1..n {
            an[i][j] = cn1[j] + cn2[j - 1];
            at[i][j] = ct1[j] + ct2[j - 1];
        }
        an[i][n] = cn2[n - 1];
        at[i][n] = ct2[n - 1];
    }
    let mut m = an.clone();
    let mut rhs: Vec<f64> = (0..n).map(|i| (th[i] - alpha).sin()).collect();
    // Condition de Kutta : γ₁ + γ_{N+1} = 0.
    let mut kutta = vec![0.0; n + 1];
    kutta[0] = 1.0;
    kutta[n] = 1.0;
    m.push(kutta);
    rhs.push(0.0);
    let g = solve(m, rhs);
    let (mut vt, mut cp) = (vec![0.0; n], vec![0.0; n]);
    let (mut cx, mut cy, mut cm) = (0.0, 0.0, 0.0);
    for i in 0..n {
        let mut v = (th[i] - alpha).cos();
        for j in 0..=n {
            v += at[i][j] * g[j];
        }
        vt[i] = v;
        cp[i] = 1.0 - v * v;
        // Contour horaire : normale extérieure (−sin θ, cos θ).
        let (nxo, nyo) = (-th[i].sin(), th[i].cos());
        let fx = -cp[i] * s[i] * nxo;
        let fy = -cp[i] * s[i] * nyo;
        cx += fx;
        cy += fy;
        cm += -(xc[i] - 0.25) * fy + yc[i] * fx;
    }
    let cl = cy * alpha.cos() - cx * alpha.sin();
    Flow { xc, yc, vt, cp, cl, cm }
}

// ------------------------------------------------------------ couche limite

fn thwaites_hl(lam: f64) -> (f64, f64) {
    if lam >= 0.0 {
        let l = lam.min(0.1);
        (2.61 - 3.75 * l + 5.24 * l * l, 0.22 + 1.57 * l - 1.8 * l * l)
    } else {
        let l = lam.max(-0.1);
        (2.088 + 0.0731 / (l + 0.14), 0.22 + 1.402 * l + 0.018 * l / (l + 0.107))
    }
}

/// Taux d'amplification de l'enveloppe eᴺ (Drela & Giles 1987).
fn dn_ds(h: f64, theta: f64, re_theta: f64) -> f64 {
    let hm = h - 1.0;
    let log10_re0 = (1.415 / hm - 0.489) * (20.0 / hm - 12.9).tanh() + 3.295 / hm + 0.44;
    if re_theta.max(1.0).log10() < log10_re0 {
        return 0.0;
    }
    let dn_dre = 0.01 * ((2.4 * h - 3.7 + 2.5 * (1.5 * h - 4.65).tanh()).powi(2) + 0.25).sqrt();
    let ell = (6.54 * h - 14.07) / (h * h);
    let m = (0.058 * (h - 4.0).powi(2) / hm - 0.068) / ell;
    dn_dre * 0.5 * (m + 1.0) * ell / theta
}

fn head_h1(h: f64) -> f64 {
    if h <= 1.6 {
        3.3 + 0.8234 * (h - 1.1).powf(-1.287)
    } else {
        3.3 + 1.5501 * (h - 0.6778).powf(-3.064)
    }
}

fn head_h(h1: f64) -> f64 {
    if h1 <= 3.32 {
        3.0
    } else if h1 >= 5.3 {
        0.86 * (h1 - 3.3).powf(-0.777) + 1.1
    } else {
        1.1538 * (h1 - 3.3).powf(-0.326) + 0.6778
    }
}

struct Layer {
    theta: Vec<f64>,
    h: Vec<f64>,
    cf: Vec<f64>,
    transition_x: f64,
    separation_x: f64,
}

/// Couche limite le long d'une face, depuis le point d'arrêt.
fn boundary_layer(s: &[f64], ue: &[f64], xs: &[f64], re: f64, n_crit: f64, x_trip: f64) -> Layer {
    let n = s.len();
    let (mut theta, mut hv, mut cf) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    let mut transition_x = f64::NAN;
    let mut separation_x = f64::NAN;
    let mut laminar = true;
    let mut amp = 0.0;
    let mut integral = 0.0;
    let mut h1 = 0.0;
    for k in 0..n {
        let u = ue[k].max(1e-4);
        if laminar {
            if k > 0 {
                let ds = s[k] - s[k - 1];
                integral += 0.5 * (ue[k - 1].max(0.0).powi(5) + u.powi(5)) * ds;
            }
            let th2 = if k == 0 || integral == 0.0 {
                // Point d'arrêt (Hiemenz) : θ² = 0,075 ν / (due/ds).
                let due = if k + 1 < n {
                    (ue[k + 1] - ue[k]) / (s[k + 1] - s[k]).max(1e-9)
                } else {
                    1.0
                };
                0.075 / (re * due.max(1e-3))
            } else {
                0.45 * integral / (re * u.powi(6))
            };
            theta[k] = th2.sqrt();
            let kk = k.max(1).min(n - 2);
            let due = (ue[kk + 1] - ue[kk - 1]) / (s[kk + 1] - s[kk - 1]).max(1e-9);
            let lam = th2 * due * re;
            let (h, ell) = thwaites_hl(lam);
            hv[k] = h;
            let re_theta = (u * theta[k] * re).max(1.0);
            cf[k] = 2.0 * ell / re_theta;
            if k > 0 {
                amp += dn_ds(h, theta[k], re_theta) * (s[k] - s[k - 1]);
            }
            let forced = x_trip.is_finite() && xs[k] >= x_trip;
            if k > 0 && (amp >= n_crit || lam < -0.09 || forced) {
                laminar = false;
                transition_x = xs[k];
                hv[k] = 1.4;
                h1 = head_h1(1.4);
            }
        } else {
            let ds = s[k] - s[k - 1];
            let u0 = ue[k - 1].max(1e-4);
            let due = (u - u0) / ds.max(1e-9);
            // Euler explicite en sous-pas.
            let (mut th, mut hh, mut h1v) = (theta[k - 1], hv[k - 1], h1);
            let sub = 4;
            for q in 0..sub {
                let uq = u0 + (u - u0) * (q as f64 + 0.5) / sub as f64;
                let re_t = (uq * th * re).max(10.0);
                let cfq = 0.246 * 10f64.powf(-0.678 * hh) * re_t.powf(-0.268);
                let dth = cfq / 2.0 - (hh + 2.0) * th * due / uq;
                let du_h1 = uq * 0.0306 * (h1v - 3.0).max(0.01).powf(-0.6169);
                let th_n = (th + dth * ds / sub as f64).max(1e-7);
                let uth_h1 = uq * th * h1v + du_h1 * ds / sub as f64;
                let u_n = u0 + (u - u0) * (q as f64 + 1.0) / sub as f64;
                h1v = uth_h1 / (u_n * th_n);
                th = th_n;
                hh = head_h(h1v);
            }
            theta[k] = th;
            hv[k] = hh;
            h1 = h1v;
            cf[k] = 0.246 * 10f64.powf(-0.678 * hh) * (u * th * re).max(10.0).powf(-0.268);
            if hh > 2.4 && !separation_x.is_finite() {
                separation_x = xs[k];
            }
        }
    }
    Layer {
        theta,
        h: hv,
        cf,
        transition_x,
        separation_x,
    }
}

struct Side {
    s: Vec<f64>,
    ue: Vec<f64>,
    xs: Vec<f64>,
    indices: Vec<usize>,
}

struct Surface {
    upper: Layer,
    lower: Layer,
    up: Side,
    lo: Side,
    cd: f64,
    dstar_node: Vec<f64>,
}

/// Couche limite sur les deux faces à partir de la solution potentielle.
fn surface_bl(flow: &Flow, re: f64, o: &PanelOptions) -> Surface {
    let n = flow.vt.len();
    // Point d'arrêt : changement de signe de vt le plus proche du bord d'attaque.
    let ile = leading_panel(&flow.xc);
    let mut ist = ile;
    let mut best = usize::MAX;
    for i in 0..n - 1 {
        if flow.vt[i] <= 0.0 && flow.vt[i + 1] > 0.0 && i.abs_diff(ile) < best {
            best = i.abs_diff(ile);
            ist = i;
        }
    }
    let side = |indices: Vec<usize>| {
        let m = indices.len() + 1;
        let (mut s, mut ue, mut xs) = (vec![0.0; m], vec![0.0; m], vec![0.0; m]);
        let (a, b) = (flow.vt[ist], flow.vt[(ist + 1).min(n - 1)]);
        let f = if a != b { a / (a - b) } else { 0.0 };
        let i1 = (ist + 1).min(n - 1);
        let xst = flow.xc[ist] + f * (flow.xc[i1] - flow.xc[ist]);
        let yst = flow.yc[ist] + f * (flow.yc[i1] - flow.yc[ist]);
        let (mut px, mut py) = (xst, yst);
        xs[0] = xst;
        for (k, &i) in indices.iter().enumerate() {
            s[k + 1] = s[k] + (flow.xc[i] - px).hypot(flow.yc[i] - py);
            ue[k + 1] = flow.vt[i].abs();
            xs[k + 1] = flow.xc[i];
            px = flow.xc[i];
            py = flow.yc[i];
        }
        Side { s, ue, xs, indices }
    };
    let up = side((ist + 1..n).collect());
    let lo = side((0..=ist).rev().collect());
    let upper = boundary_layer(&up.s, &up.ue, &up.xs, re, o.n_crit, o.x_trip_upper.unwrap_or(f64::NAN));
    let lower = boundary_layer(&lo.s, &lo.ue, &lo.xs, re, o.n_crit, o.x_trip_lower.unwrap_or(f64::NAN));
    let squire = |bl: &Layer, sd: &Side| {
        // Évalué à x/c ≈ 0,98 : au-delà, la décélération potentielle vers le
        // bord de fuite fermé n'existe pas dans l'écoulement réel.
        let mut k = sd.s.len() - 1;
        while k > 1 && sd.xs[k] > o.squire_x {
            k -= 1;
        }
        let th = bl.theta[k];
        let ue_te = sd.ue[k].max(0.05);
        let h = bl.h[k].min(3.0);
        2.0 * th * ue_te.powf((h + 5.0) / 2.0)
    };
    let cd = squire(&upper, &up) + squire(&lower, &lo);
    // δ* aux nœuds (moyenne des panneaux adjacents), figé au-delà de x/c = 0,95.
    let mut dp = vec![0.0; n];
    for (sd, bl) in [(&up, &upper), (&lo, &lower)] {
        for (k, &i) in sd.indices.iter().enumerate() {
            dp[i] = bl.theta[k + 1] * bl.h[k + 1];
        }
        let mut r = 0.0;
        for (k, &i) in sd.indices.iter().enumerate() {
            if sd.xs[k + 1] <= 0.95 {
                r = dp[i];
            } else {
                dp[i] = r;
            }
        }
    }
    let dstar_node = (0..=n)
        .map(|k| (0.5 * (dp[k.saturating_sub(1)] + dp[k.min(n - 1)])).min(0.05))
        .collect();
    Surface {
        upper,
        lower,
        up,
        lo,
        cd,
        dstar_node,
    }
}

fn leading_panel(xc: &[f64]) -> usize {
    let mut k = 0;
    for i in 1..xc.len() {
        if xc[i] < xc[k] {
            k = i;
        }
    }
    k
}

#[derive(Clone, Copy, Debug)]
pub struct PanelOptions {
    pub re: f64,
    pub n_crit: f64,
    /// Transition imposée (x/c) ; None : libre.
    pub x_trip_upper: Option<f64>,
    pub x_trip_lower: Option<f64>,
    pub iterations: usize,
    pub relax: f64,
    pub squire_x: f64,
}

impl PanelOptions {
    /// Transition libre (eᴺ, N = 9).
    pub fn free(re: f64) -> PanelOptions {
        PanelOptions {
            re,
            n_crit: 9.0,
            x_trip_upper: None,
            x_trip_lower: None,
            iterations: 1,
            relax: 0.5,
            squire_x: 0.98,
        }
    }

    /// Transition imposée sur les deux faces (0 : couche limite turbulente partout).
    pub fn tripped(re: f64, x: f64) -> PanelOptions {
        PanelOptions {
            x_trip_upper: Some(x),
            x_trip_lower: Some(x),
            ..PanelOptions::free(re)
        }
    }
}

/// Frottement local (rapporté à ue) le long d'une face ; cf·ue² le rapporte à V∞.
#[derive(Clone, Debug, Default)]
pub struct Friction {
    pub x: Vec<f64>,
    pub cf: Vec<f64>,
    pub ue: Vec<f64>,
}

#[derive(Clone, Debug)]
pub struct PanelResult {
    pub alpha_deg: f64,
    pub cl: f64,
    pub cd: f64,
    pub cm: f64,
    pub cl_inviscid: f64,
    pub transition_upper: f64,
    pub transition_lower: f64,
    pub separation_upper: f64,
    pub separation_lower: f64,
    /// Décollement turbulent avant 95 % de corde (hors du domaine de validité).
    pub separated: bool,
    pub friction_upper: Friction,
    pub friction_lower: Friction,
    /// Cp aux points de contrôle (bord de fuite → intrados → bord d'attaque → extrados).
    pub cp_x: Vec<f64>,
    pub cp_y: Vec<f64>,
    pub cp: Vec<f64>,
    pub i_le: usize,
}

/// Point de polaire : potentiel + couche limite + couplage faible.
pub fn panel_solve(shape: &PanelShape, alpha_deg: f64, o: &PanelOptions) -> PanelResult {
    let base = airfoil_points(shape, 80);
    let alpha = alpha_deg.to_radians();
    let inviscid = potential_flow(&base, alpha);
    let mut flow = inviscid.clone();
    let nb = base.len();
    let mut dstar = vec![0.0; nb];
    let mut it = 0;
    let bl = loop {
        let bl = surface_bl(&flow, o.re, o);
        if it == o.iterations {
            break bl;
        }
        it += 1;
        for k in 0..nb {
            dstar[k] = (1.0 - o.relax) * dstar[k] + o.relax * bl.dstar_node[k];
        }
        // Déplacement selon la normale extérieure (contour horaire : (−t_y, t_x)).
        let disp: Vec<[f64; 2]> = (0..nb)
            .map(|k| {
                let (kp, km) = ((k + 1).min(nb - 1), k.saturating_sub(1));
                let (tx, ty) = (base[kp][0] - base[km][0], base[kp][1] - base[km][1]);
                let tl = tx.hypot(ty).max(1e-300);
                [-dstar[k] * ty / tl, dstar[k] * tx / tl]
            })
            .collect();
        // Bord de fuite refermé sur la position moyenne des deux faces épaissies.
        let dte = [
            0.5 * (disp[0][0] + disp[nb - 1][0]),
            0.5 * (disp[0][1] + disp[nb - 1][1]),
        ];
        let pts: Vec<[f64; 2]> = base
            .iter()
            .zip(&disp)
            .map(|(&[x, y], &[mut dx, mut dy])| {
                if x > 0.9 {
                    let w = ((x - 0.9) / 0.1).powi(2);
                    dx += w * (dte[0] - dx);
                    dy += w * (dte[1] - dy);
                }
                [x + dx, y + dy]
            })
            .collect();
        flow = potential_flow(&pts, alpha);
    };
    let fr = |sd: &Side, l: &Layer| Friction {
        x: sd.xs.clone(),
        cf: l.cf.clone(),
        ue: sd.ue.clone(),
    };
    let sep = |x: f64| x.is_finite() && x < 0.95;
    PanelResult {
        alpha_deg,
        cl: flow.cl,
        cd: bl.cd,
        cm: flow.cm,
        cl_inviscid: inviscid.cl,
        transition_upper: bl.upper.transition_x,
        transition_lower: bl.lower.transition_x,
        separation_upper: bl.upper.separation_x,
        separation_lower: bl.lower.separation_x,
        separated: sep(bl.upper.separation_x) || sep(bl.lower.separation_x),
        friction_upper: fr(&bl.up, &bl.upper),
        friction_lower: fr(&bl.lo, &bl.lower),
        i_le: leading_panel(&flow.xc),
        cp_x: flow.xc,
        cp_y: flow.yc,
        cp: flow.cp,
    }
}

/// Coefficients de section corrigés du décollement.
#[derive(Clone, Copy, Debug)]
pub struct Section {
    pub cl: f64,
    pub cd: f64,
    pub cm: f64,
    /// Position relative du point de décollement (1 : attaché).
    pub f: f64,
    pub kirchhoff: bool,
}

/// Écoulement de Kirchhoff (forme de Beddoes et Leishman) : Cn = Cn_att ((1 + √f)/2)²,
/// aspiration Cc = Cc_att √f. Modèle approché : allure du décrochage, pas son incidence exacte.
pub fn section_coefficients(r: &PanelResult) -> Section {
    let f_of = |x: f64| if x.is_finite() { (x / 0.97).clamp(0.0, 1.0) } else { 1.0 };
    let f = if r.separated {
        f_of(r.separation_upper).min(f_of(r.separation_lower))
    } else {
        1.0
    };
    if f >= 1.0 {
        return Section {
            cl: r.cl,
            cd: r.cd,
            cm: r.cm,
            f: 1.0,
            kirchhoff: false,
        };
    }
    let a = r.alpha_deg.to_radians();
    let cn_att = r.cl * a.cos();
    let cc_att = r.cl * a.sin();
    let sf = f.sqrt();
    let cn = cn_att * ((1.0 + sf) / 2.0).powi(2);
    let cc = cc_att * sf;
    Section {
        cl: cn * a.cos() + cc * a.sin(),
        cd: r.cd + (cn * a.sin() - cc * a.cos()).max(0.0),
        cm: r.cm - cn * 0.15 * (1.0 - f),
        f,
        kirchhoff: true,
    }
}

/// Cp critique (écoulement localement sonique), γ = 1,4.
pub fn critical_cp(mach: f64) -> f64 {
    let g = 1.4;
    let m2 = mach * mach;
    2.0 / (g * m2) * (((2.0 + (g - 1.0) * m2) / (g + 1.0)).powf(g / (g - 1.0)) - 1.0)
}

/// Corrections de compressibilité subsonique (Prandtl–Glauert, Kármán–Tsien).
#[derive(Clone, Copy, Debug)]
pub struct Compressible {
    pub beta: f64,
    /// Facteur de Prandtl–Glauert sur Cz et Cm.
    pub factor: f64,
    /// Cp minimal corrigé (Kármán–Tsien) et Cp critique.
    pub cp_min: f64,
    pub cp_star: f64,
    /// Zone supersonique locale : correction hors domaine.
    pub critical: bool,
}

pub fn compressible_correction(mach: f64, cp_min: f64) -> Compressible {
    let m = mach.min(0.95);
    let beta = (1.0 - m * m).sqrt();
    let kt = |cp: f64| cp / (beta + m * m / (1.0 + beta) * (cp / 2.0));
    let cp_min_c = kt(cp_min);
    let cp_star = critical_cp(m);
    Compressible {
        beta,
        factor: 1.0 / beta,
        cp_min: cp_min_c,
        cp_star,
        critical: cp_min_c < cp_star,
    }
}
