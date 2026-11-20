//! Aéroacoustique 2D en post-traitement (processeur) : analogie de
//! Ffowcs Williams–Hawkings sur une surface d'intégration fixe dans un
//! écoulement uniforme U selon x, dans le domaine fréquentiel (Lockard 2000,
//! J. Sound Vib. 229) ; Curle compact à partir des efforts ; transformée de
//! Fourier discrète fenêtrée.
//!
//! Convention : x(t) = Re[x̂ e^{iωt}]. Repère de la maquette : la surface et
//! les observateurs sont fixes, le milieu s'écoule à U. Équation :
//! (1/c²)(∂t + U∂x)² p′ − ∇²p′ = ∂t[Q δ(f)] − ∂i[F_i δ(f)] (quadripôles
//! extérieurs à la surface négligés), avec n normale sortante (vers le fluide) :
//!   Q   = ρ u_n − ρ0 U_n,
//!   F_i = p′ n_i + ρ (u_i − 2 U_i) u_n + ρ0 U_i U_n,
//! d'où p̂(x) = ∮ (iω Q̂ G + F̂_i ∂G/∂y_i) dl, avec la fonction de Green
//! convectée G = −i/(4β) · e^{i M k d₁/β²} · H₀⁽²⁾(k r̄/β²),
//! d = x − y, r̄ = √(d₁² + β² d₂²), β = √(1 − M²).
//!
//! La formulation 1A de Farassat est temporelle et tridimensionnelle (temps
//! retardé unique) ; en 2D la fonction de Green garde une « queue » et la
//! forme fréquentielle ci-dessus en est l'équivalent exact.

use std::f64::consts::PI;
use std::ops::{Add, AddAssign, Mul, Sub};

/// Nombre complexe minimal.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Cx {
    pub re: f64,
    pub im: f64,
}

impl Cx {
    pub const fn new(re: f64, im: f64) -> Cx {
        Cx { re, im }
    }
    pub fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
    pub fn arg(self) -> f64 {
        self.im.atan2(self.re)
    }
    pub fn expi(phi: f64) -> Cx {
        Cx::new(phi.cos(), phi.sin())
    }
    pub fn scale(self, s: f64) -> Cx {
        Cx::new(self.re * s, self.im * s)
    }
    pub fn conj(self) -> Cx {
        Cx::new(self.re, -self.im)
    }
}
impl Add for Cx {
    type Output = Cx;
    fn add(self, o: Cx) -> Cx {
        Cx::new(self.re + o.re, self.im + o.im)
    }
}
impl Sub for Cx {
    type Output = Cx;
    fn sub(self, o: Cx) -> Cx {
        Cx::new(self.re - o.re, self.im - o.im)
    }
}
impl Mul for Cx {
    type Output = Cx;
    fn mul(self, o: Cx) -> Cx {
        Cx::new(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    }
}
impl AddAssign for Cx {
    fn add_assign(&mut self, o: Cx) {
        self.re += o.re;
        self.im += o.im;
    }
}
const I: Cx = Cx::new(0.0, 1.0);

// Fonctions de Bessel J0, J1, Y0, Y1 d'argument réel positif (approximations
// rationnelles de Numerical Recipes, erreur relative ≈ 1e-8).
pub fn bessel_j0(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 8.0 {
        let y = x * x;
        let a = 57568490574.0
            + y * (-13362590354.0 + y * (651619640.7 + y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
        let b = 57568490411.0 + y * (1029532985.0 + y * (9494680.718 + y * (59272.64853 + y * (267.8532712 + y))));
        a / b
    } else {
        let z = 8.0 / ax;
        let y = z * z;
        let xx = ax - 0.785398164;
        let a = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4 + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
        let b = -0.1562499995e-1 + y * (0.1430488765e-3 + y * (-0.6911147651e-5 + y * (0.7621095161e-6 - y * 0.934935152e-7)));
        (0.636619772 / ax).sqrt() * (xx.cos() * a - z * xx.sin() * b)
    }
}

pub fn bessel_j1(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 8.0 {
        let y = x * x;
        let a = x
            * (72362614232.0
                + y * (-7895059235.0 + y * (242396853.1 + y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))));
        let b = 144725228442.0 + y * (2300535178.0 + y * (18583304.74 + y * (99447.43394 + y * (376.9991397 + y))));
        a / b
    } else {
        let z = 8.0 / ax;
        let y = z * z;
        let xx = ax - 2.356194491;
        let a = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4 + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
        let b = 0.04687499995 + y * (-0.2002690873e-3 + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)));
        let r = (0.636619772 / ax).sqrt() * (xx.cos() * a - z * xx.sin() * b);
        if x < 0.0 {
            -r
        } else {
            r
        }
    }
}

pub fn bessel_y0(x: f64) -> f64 {
    if x < 8.0 {
        let y = x * x;
        let a = -2957821389.0
            + y * (7062834065.0 + y * (-512359803.6 + y * (10879881.29 + y * (-86327.92757 + y * 228.4622733))));
        let b = 40076544269.0 + y * (745249964.8 + y * (7189466.438 + y * (47447.26470 + y * (226.1030244 + y))));
        a / b + 0.636619772 * bessel_j0(x) * x.ln()
    } else {
        let z = 8.0 / x;
        let y = z * z;
        let xx = x - 0.785398164;
        let a = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4 + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
        let b = -0.1562499995e-1 + y * (0.1430488765e-3 + y * (-0.6911147651e-5 + y * (0.7621095161e-6 + y * (-0.934945152e-7))));
        (0.636619772 / x).sqrt() * (xx.sin() * a + z * xx.cos() * b)
    }
}

pub fn bessel_y1(x: f64) -> f64 {
    if x < 8.0 {
        let y = x * x;
        let a = x
            * (-0.4900604943e13
                + y * (0.1275274390e13 + y * (-0.5153438139e11 + y * (0.7349264551e9 + y * (-0.4237922726e7 + y * 0.8511937935e4)))));
        let b = 0.2499580570e14
            + y * (0.4244419664e12 + y * (0.3733650367e10 + y * (0.2245904002e8 + y * (0.1020426050e6 + y * (0.3549632885e3 + y)))));
        a / b + 0.636619772 * (bessel_j1(x) * x.ln() - 1.0 / x)
    } else {
        let z = 8.0 / x;
        let y = z * z;
        let xx = x - 2.356194491;
        let a = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4 + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
        let b = 0.04687499995 + y * (-0.2002690873e-3 + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)));
        (0.636619772 / x).sqrt() * (xx.sin() * a + z * xx.cos() * b)
    }
}

/// Milieu : masse volumique et célérité au repos, vitesse de l'écoulement selon x.
#[derive(Clone, Copy, Debug)]
pub struct Medium {
    pub rho0: f64,
    pub c: f64,
    pub u: f64,
}

impl Medium {
    /// Unités du réseau D2Q9 : ρ0 = 1, c = 1/√3.
    pub fn lattice(u: f64) -> Medium {
        Medium { rho0: 1.0, c: 1.0 / 3f64.sqrt(), u }
    }
    pub fn mach(&self) -> f64 {
        self.u / self.c
    }
}

/// Fonction de Green convectée G(x, y, ω) et son gradient par rapport à la
/// source y ; `d` = x − y.
pub fn green(m: &Medium, omega: f64, d: [f64; 2]) -> (Cx, [Cx; 2]) {
    let mach = m.mach();
    let b2 = 1.0 - mach * mach;
    let b = b2.sqrt();
    let k = omega / m.c;
    let rb = (d[0] * d[0] + b2 * d[1] * d[1]).sqrt();
    let z = k * rb / b2;
    let h0 = Cx::new(bessel_j0(z), -bessel_y0(z));
    let h1 = Cx::new(bessel_j1(z), -bessel_y1(z));
    let a = Cx::new(0.0, -1.0 / (4.0 * b)) * Cx::expi(mach * k * d[0] / b2);
    let g = a * h0;
    // ∂G/∂d₁ = A e^{iφ}[i M k/β² H₀ − H₁ (k/β²) d₁/r̄], ∂G/∂d₂ = −A e^{iφ} H₁ k d₂/r̄.
    let gd1 = a * ((I * h0).scale(mach * k / b2) - h1.scale(k / b2 * d[0] / rb));
    let gd2 = a * h1.scale(-k * d[1] / rb);
    // ∂/∂y = −∂/∂d.
    (g, [gd1.scale(-1.0), gd2.scale(-1.0)])
}

/// Amplitude complexe x̂ à la pulsation ω (x(t) ≈ Re[x̂ e^{iωt}]) d'une série
/// échantillonnée au pas `dt` à partir de t0, moyenne retirée, fenêtre de Hann.
pub fn dft_amp(x: &[f64], t0: f64, dt: f64, omega: f64) -> Cx {
    let n = x.len();
    if n < 2 {
        return Cx::default();
    }
    let m = x.iter().sum::<f64>() / n as f64;
    let mut s = Cx::default();
    let mut ws = 0.0;
    for (j, &v) in x.iter().enumerate() {
        let w = 0.5 - 0.5 * (2.0 * PI * j as f64 / (n - 1) as f64).cos();
        ws += w;
        s += Cx::expi(-omega * (t0 + j as f64 * dt)).scale(w * (v - m));
    }
    s.scale(2.0 / ws)
}

/// Fréquence (cycles par unité de temps) du pic d'une série, affinée par
/// recherche de l'amplitude maximale de `dft_amp` autour de `f_guess`.
pub fn refine_peak(x: &[f64], dt: f64, f_guess: f64) -> f64 {
    let t = x.len() as f64 * dt;
    let (mut lo, mut hi) = (f_guess - 1.0 / t, f_guess + 1.0 / t);
    let amp = |f: f64| dft_amp(x, 0.0, dt, 2.0 * PI * f).abs();
    // Section dorée.
    let g = 0.618_033_988_75;
    let (mut a, mut b) = (hi - g * (hi - lo), lo + g * (hi - lo));
    let (mut fa, mut fb) = (amp(a), amp(b));
    for _ in 0..60 {
        if fa > fb {
            hi = b;
            b = a;
            fb = fa;
            a = hi - g * (hi - lo);
            fa = amp(a);
        } else {
            lo = a;
            a = b;
            fa = fb;
            b = lo + g * (hi - lo);
            fb = amp(b);
        }
    }
    0.5 * (lo + hi)
}

/// Surface d'intégration : points, normales sortantes unitaires, longueurs.
#[derive(Clone, Debug)]
pub struct Surface {
    pub pts: Vec<[f64; 2]>,
    pub normals: Vec<[f64; 2]>,
    pub dl: Vec<f64>,
}

impl Surface {
    /// Cercle de centre `c`, rayon `r`, `n` points (règle des rectangles,
    /// spectrale pour une fonction périodique régulière).
    pub fn circle(c: [f64; 2], r: f64, n: usize) -> Surface {
        let mut s = Surface { pts: Vec::new(), normals: Vec::new(), dl: Vec::new() };
        for j in 0..n {
            let a = 2.0 * PI * (j as f64 + 0.5) / n as f64;
            let (sa, ca) = a.sin_cos();
            s.pts.push([c[0] + r * ca, c[1] + r * sa]);
            s.normals.push([ca, sa]);
            s.dl.push(2.0 * PI * r / n as f64);
        }
        s
    }
    pub fn len(&self) -> usize {
        self.pts.len()
    }
    pub fn is_empty(&self) -> bool {
        self.pts.is_empty()
    }
}

/// Termes sources Q̂, F̂₁, F̂₂ de chaque point de la surface à la pulsation ω.
/// `rho`, `ux`, `uy` : séries [échantillon][point], relevées aux instants t0 + j·dt.
#[allow(clippy::too_many_arguments)]
pub fn fwh_sources(
    m: &Medium,
    s: &Surface,
    rho: &[Vec<f64>],
    ux: &[Vec<f64>],
    uy: &[Vec<f64>],
    t0: f64,
    dt: f64,
    omega: f64,
) -> Vec<[Cx; 3]> {
    let ns = rho.len();
    let c2 = m.c * m.c;
    (0..s.len())
        .map(|j| {
            let [n1, n2] = s.normals[j];
            let un0 = m.u * n1;
            let mut q = vec![0.0; ns];
            let mut f1 = vec![0.0; ns];
            let mut f2 = vec![0.0; ns];
            for t in 0..ns {
                let (r, u, v) = (rho[t][j], ux[t][j], uy[t][j]);
                let un = u * n1 + v * n2;
                let p = c2 * (r - m.rho0);
                q[t] = r * un - m.rho0 * un0;
                f1[t] = p * n1 + r * (u - 2.0 * m.u) * un + m.rho0 * m.u * un0;
                f2[t] = p * n2 + r * v * un;
            }
            [dft_amp(&q, t0, dt, omega), dft_amp(&f1, t0, dt, omega), dft_amp(&f2, t0, dt, omega)]
        })
        .collect()
}

/// Pression acoustique p̂ à l'observateur `x` à partir des sources de `fwh_sources`.
pub fn fwh_pressure(m: &Medium, s: &Surface, src: &[[Cx; 3]], omega: f64, x: [f64; 2]) -> Cx {
    let mut p = Cx::default();
    for j in 0..s.len() {
        let y = s.pts[j];
        let (g, dg) = green(m, omega, [x[0] - y[0], x[1] - y[1]]);
        let [q, f1, f2] = src[j];
        p += ((I * q).scale(omega) * g + f1 * dg[0] + f2 * dg[1]).scale(s.dl[j]);
    }
    p
}

/// Curle compact : dipôle de la force F̂ exercée par le fluide sur la maquette
/// (par unité d'envergure), placé en `y`.
pub fn curle_compact(m: &Medium, force: [Cx; 2], omega: f64, y: [f64; 2], x: [f64; 2]) -> Cx {
    let (_, dg) = green(m, omega, [x[0] - y[0], x[1] - y[1]]);
    (force[0] * dg[0] + force[1] * dg[1]).scale(-1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bessel_values() {
        // Valeurs de référence (Abramowitz & Stegun).
        assert!((bessel_j0(1.0) - 0.765_197_686_6).abs() < 1e-7);
        assert!((bessel_y0(1.0) - 0.088_256_964_2).abs() < 1e-7);
        assert!((bessel_j1(2.5) - 0.497_094_102_5).abs() < 1e-7);
        assert!((bessel_y1(10.0) - 0.249_015_424_2).abs() < 1e-7);
        assert!((bessel_j0(12.0) - 0.047_689_310_8).abs() < 1e-7);
    }

    /// Résidu de l'opérateur convecté sur G loin de la source (différences finies).
    #[test]
    fn green_satisfies_convected_helmholtz() {
        let m = Medium::lattice(0.2 / 3f64.sqrt());
        let omega = 0.01;
        let k = omega / m.c;
        let h = 0.05;
        let gat = |x: f64, y: f64| green(&m, omega, [x, y]).0;
        for &(x, y) in &[(30.0, 10.0), (-25.0, 5.0), (3.0, -40.0)] {
            let g = gat(x, y);
            let gxx = (gat(x + h, y) + gat(x - h, y) - g.scale(2.0)).scale(1.0 / (h * h));
            let gyy = (gat(x, y + h) + gat(x, y - h) - g.scale(2.0)).scale(1.0 / (h * h));
            let gx = (gat(x + h, y) - gat(x - h, y)).scale(0.5 / h);
            // (1/c²)(iω + U∂x)² G − ∇²G = −k²G + 2ikM ∂xG + M² ∂xxG − ∇²G.
            let mach = m.mach();
            let r = g.scale(-k * k) + (I * gx).scale(2.0 * k * mach) + gxx.scale(mach * mach) - gxx - gyy;
            assert!(r.abs() < 1e-4 * (k * k * g.abs()), "résidu {r:?} en ({x}, {y})");
            // Gradient analytique contre différences finies (∂/∂y = −∂/∂x).
            let dg = green(&m, omega, [x, y]).1;
            assert!((dg[0] + gx).abs() < 1e-5 * gx.abs().max(1e-6));
        }
    }

    /// Source ponctuelle de potentiel φ = G(x, x_s) : p′ = −ρ0 (iω + U∂x) φ,
    /// u′ = ∇φ ; FW-H sur un cercle autour de la source doit rendre p̂ dehors
    /// et zéro dedans.
    #[test]
    fn fwh_reproduces_monopole() {
        let m = Medium::lattice(0.2 / 3f64.sqrt());
        let omega = 2.0 * PI * 0.002;
        let xs = [1.5, -0.7];
        let h = 1e-3;
        let phi = |x: f64, y: f64| green(&m, omega, [x - xs[0], y - xs[1]]).0;
        let fields = |x: f64, y: f64| {
            let ph = phi(x, y);
            let px = (phi(x + h, y) - phi(x - h, y)).scale(0.5 / h);
            let py = (phi(x, y + h) - phi(x, y - h)).scale(0.5 / h);
            let p = ((I * ph).scale(omega) + px.scale(m.u)).scale(-m.rho0);
            (p, px, py)
        };
        let s = Surface::circle([0.0, 0.0], 12.0, 128);
        // Séries temporelles de petite amplitude : 8 périodes, 64 échantillons par période.
        let eps = 1e-3;
        let per = 2.0 * PI / omega;
        let dt = per / 64.0;
        let ns = 8 * 64 + 1;
        let mut rho = vec![vec![0.0; s.len()]; ns];
        let mut ux = rho.clone();
        let mut uy = rho.clone();
        for (j, y) in s.pts.iter().enumerate() {
            let (p, u, v) = fields(y[0], y[1]);
            for t in 0..ns {
                let e = Cx::expi(omega * t as f64 * dt);
                rho[t][j] = m.rho0 + eps * (p * e).re / (m.c * m.c);
                ux[t][j] = m.u + eps * (u * e).re;
                uy[t][j] = eps * (v * e).re;
            }
        }
        let src = fwh_sources(&m, &s, &rho, &ux, &uy, 0.0, dt, omega);
        for &x in &[[60.0, 20.0], [-50.0, -30.0], [5.0, 80.0]] {
            let want = fields(x[0], x[1]).0.scale(eps);
            let got = fwh_pressure(&m, &s, &src, omega, x);
            let err = (got - want).abs() / want.abs();
            assert!(err < 2e-2, "observateur {x:?} : {got:?} contre {want:?} ({err:.2e})");
        }
        let inside = fwh_pressure(&m, &s, &src, omega, [-4.0, 3.0]);
        let ref_in = fields(-4.0, 3.0).0.scale(eps).abs();
        assert!(inside.abs() < 2e-2 * ref_in, "intérieur {inside:?}");
    }
}
