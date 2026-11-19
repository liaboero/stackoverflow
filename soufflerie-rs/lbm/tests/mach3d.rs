//! Domaine de Mach du D3Q19 (ignoré sans adaptateur GPU) : défaut cubique du
//! moment d'équilibre (ρ u_α³ manquant) et sa correction (`MACH3_BIT`).
//! Veine périodique 3D sans maquette (`set_periodic_x`, `set_state`).
//!
//! 1. Onde de cisaillement décroissante advectée selon x à Ma 0 / 0,1 / 0,2
//!    / 0,3 : nombre d'onde selon (1, 1, 0), vitesse selon (1, −1, 0) ; le
//!    réseau voit des déformations normales ∂_x u_x = −∂_y u_y, où le
//!    défaut réduit σ_xx du facteur 1 − 4,5 u_x². Viscosité apparente
//!    attendue sans correction : ν (1 − 2,25 U²) (README, « Domaine de Mach »).
//!    Onde alignée (k selon x, vitesse selon y) : sans défaut (moments
//!    xxy portés par l'équilibre), témoin.
//! 2. Tourbillon isentrope (isotherme) convecté d'une période de veine :
//!    écart au même tourbillon au repos après le même temps (invariance
//!    galiléenne), avec et sans correction.
//!
//! Par défaut Ma = 0,3 seul ; tableau complet (Ma 0,1 / 0,2 / 0,3, ≈ 9 min
//! sur GTX 1070 partagée) avec `MACH3_TABLEAU=1`.
#![cfg(feature = "gpu")]

use lbm::bodies3d::{Body3D, Bounds, Shape3D};
use lbm::cpu::Ramp;
use lbm::gpu3d::{Gpu3D, MACH3_BIT};
use std::f64::consts::PI;

/// Veine périodique nx × ny × nz sans maquette, viscosité ν.
/// Nombres de Mach essayés.
fn machs() -> Vec<f64> {
    if std::env::var("MACH3_TABLEAU").is_ok() {
        vec![0.1, 0.2, 0.3]
    } else {
        vec![0.3]
    }
}

fn veine(nx: usize, ny: usize, nz: usize, nu: f64) -> Option<Gpu3D> {
    let mut body = Body3D::new(Shape3D::Sphere, 4, 0.0);
    body.nx = nx;
    body.ny = ny;
    body.nz = nz;
    // Sphère hors de la veine : aucune maille solide.
    body.center = [1e6; 3];
    body.bounds = Bounds {
        x0: 2.0,
        x1: 2.0,
        y0: 0.0,
        y1: 0.0,
        z0: 0.0,
        z1: 0.0,
    };
    match Gpu3D::new(&body, nu, 0.0, false) {
        Ok(mut g) => {
            assert_eq!(g.solid_count, 0);
            g.set_periodic_x();
            Some(g)
        }
        Err(e) => {
            eprintln!("GPU ignoré : {e}");
            None
        }
    }
}

const RAMP0: Ramp = Ramp { steps: 1, u: 0.0 };

/// Onde de cisaillement : vitesse d'advection U selon x, nombre d'onde
/// (kx, ky) en périodes de veine, amplitude δ selon la direction normale à
/// k dans le plan (x, y). Taux de décroissance mesuré / ν k² (théorique).
fn onde(n: usize, nu: f64, u: f64, per: [f64; 2], corr: bool) -> Option<f64> {
    let mut g = veine(n, n, 2, nu)?;
    g.set_essai(if corr { MACH3_BIT } else { 0 }, 0.0);
    let kv = [2.0 * PI * per[0] / n as f64, 2.0 * PI * per[1] / n as f64];
    let k = (kv[0] * kv[0] + kv[1] * kv[1]).sqrt();
    let t = [-kv[1] / k, kv[0] / k];
    let d = 1e-3;
    let mut m = vec![[0f32; 4]; g.n];
    for z in 0..2 {
        for y in 0..n {
            for x in 0..n {
                let s = d * (kv[0] * x as f64 + kv[1] * y as f64).sin();
                m[x + n * (y + n * z)] = [(u + s * t[0]) as f32, (s * t[1]) as f32, 0.0, 1.0];
            }
        }
    }
    g.set_state(&m).unwrap();
    // Amplitude de la composante de Fourier (kx, ky) de u · t, plan z = 0.
    let amp = |g: &mut Gpu3D| {
        let m = g.macros();
        let (mut a, mut b) = (0.0, 0.0);
        for y in 0..n {
            for x in 0..n {
                let c = m[x + n * y];
                let v = (c[0] as f64 - u) * t[0] + c[1] as f64 * t[1];
                let ph = kv[0] * x as f64 + kv[1] * y as f64;
                a += v * ph.cos();
                b += v * ph.sin();
            }
        }
        (a * a + b * b).sqrt()
    };
    // Décroissance de e^{-1} environ ; ajustement entre 20 % et 100 %.
    let tt = (1.0 / (nu * k * k)).round() as usize;
    let t0 = tt / 5;
    g.run(t0, &RAMP0);
    let a0 = amp(&mut g);
    g.run(tt - t0, &RAMP0);
    let a1 = amp(&mut g);
    assert!(!g.lost && a1.is_finite());
    let rate = (a0 / a1).ln() / (tt - t0) as f64;
    Some(rate / (nu * k * k))
}

#[test]
fn onde_cisaillement_advectee() {
    let (n, tau) = (64usize, 0.8);
    let nu = (tau - 0.5) / 3.0;
    let cs = 1.0 / 3f64.sqrt();
    // Références au repos (erreur de discrétisation, retirée des rapports).
    let Some(r_diag) = onde(n, nu, 0.0, [1.0, 1.0], false) else { return };
    let r_al = onde(n, nu, 0.0, [1.0, 0.0], false).unwrap();
    println!("au repos : ν_app/ν = {r_diag:.5} (k selon (1,1)), {r_al:.5} (k selon x)");
    println!("  Ma  | diag. sans | diag. avec | attendu sans | aligné sans | aligné avec");
    for ma in machs() {
        let u = ma * cs;
        let ds = onde(n, nu, u, [1.0, 1.0], false).unwrap() / r_diag;
        let da = onde(n, nu, u, [1.0, 1.0], true).unwrap() / r_diag;
        let als = onde(n, nu, u, [1.0, 0.0], false).unwrap() / r_al;
        let ala = onde(n, nu, u, [1.0, 0.0], true).unwrap() / r_al;
        let th = 1.0 - 2.25 * u * u;
        println!(" {ma:.1}  |  {ds:.5}   |  {da:.5}   |   {th:.5}    |  {als:.5}    |  {ala:.5}");
        // Défaut mesuré conforme à l'analyse (à 15 % près de l'écart).
        assert!(((1.0 - ds) / (1.0 - th) - 1.0).abs() < 0.15, "défaut sans correction {ds} contre {th}");
        // Correction : écart résiduel inférieur au dixième du défaut.
        assert!((da - 1.0).abs() < 0.1 * (1.0 - th), "avec correction {da}");
        // Onde alignée : pas de défaut cubique, correction sans effet notable.
        assert!((als - 1.0).abs() < 2e-3 && (ala - 1.0).abs() < 2e-3, "onde alignée {als} {ala}");
    }
}

/// Tourbillon isotherme (Gaussien, u_θ = ε (r/R) e^{(1 − r²/R²)/2}, ρ en
/// équilibre radial) dans une veine périodique n × n, convecté à U pendant
/// `steps` pas. Champ de perturbation (u − U, ρ) final.
fn tourbillon(n: usize, nu: f64, u: f64, steps: usize, corr: bool) -> Option<Vec<[f64; 3]>> {
    let mut g = veine(n, n, 2, nu)?;
    g.set_essai(if corr { MACH3_BIT } else { 0 }, 0.0);
    let (eps, r) = (0.02, n as f64 / 10.0);
    let c = n as f64 / 2.0;
    let mut m = vec![[0f32; 4]; g.n];
    for z in 0..2 {
        for y in 0..n {
            for x in 0..n {
                let (dx, dy) = (x as f64 - c, y as f64 - c);
                let q = (dx * dx + dy * dy) / (r * r);
                let e = (0.5 * (1.0 - q)).exp();
                let rho = (-1.5 * eps * eps * (1.0 - q).exp()).exp();
                m[x + n * (y + n * z)] = [(u - eps * dy / r * e) as f32, (eps * dx / r * e) as f32, 0.0, rho as f32];
            }
        }
    }
    g.set_state(&m).unwrap();
    g.run(steps, &RAMP0);
    assert!(!g.lost);
    let m = g.macros();
    Some((0..n * n).map(|i| [m[i][0] as f64 - u, m[i][1] as f64, m[i][3] as f64]).collect())
}

#[test]
fn tourbillon_convecte() {
    let (n, tau) = (64usize, 0.56);
    let nu = (tau - 0.5) / 3.0;
    let cs = 1.0 / 3f64.sqrt();
    println!("  Ma  | écart sans | écart avec (norme L2 de u − U rapportée au tourbillon au repos)");
    let mut res = Vec::new();
    for ma in machs() {
        // Vitesse ajustée pour un nombre entier de pas par période de veine.
        let steps = (n as f64 / (ma * cs)).round() as usize;
        let u = n as f64 / steps as f64;
        let Some(rest) = tourbillon(n, nu, 0.0, steps, false) else { return };
        let norm: f64 = rest.iter().map(|p| p[0] * p[0] + p[1] * p[1]).sum::<f64>().sqrt();
        let err = |f: &[[f64; 3]]| f.iter().zip(&rest).map(|(a, b)| (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sum::<f64>().sqrt() / norm;
        let es = err(&tourbillon(n, nu, u, steps, false).unwrap());
        let ea = err(&tourbillon(n, nu, u, steps, true).unwrap());
        println!(" {ma:.1}  |  {es:.2e}  |  {ea:.2e}");
        res.push((es, ea));
    }
    // À Ma 0,3, la correction réduit l'écart galiléen.
    let (es, ea) = *res.last().unwrap();
    assert!(ea < es, "tourbillon : correction sans effet ({es} → {ea})");
}

#[test]
fn collision_par_defaut_sans_correction() {
    // Le code de la correction n'est compilé qu'avec le bit (défaut inchangé).
    let s = Gpu3D::shader_source(false);
    assert!(!s.contains("sxx") && !s.contains("M3["));
}
