//! Validation physique du moteur CPU : mêmes cas et mêmes bornes que la suite
//! de l'application web (tests/validate.mjs), valeurs JS rappelées.

use lbm::cpu::Macros;
use lbm::stats::{analyze, Analysis};
use lbm::wall::WallSampler;
use lbm::{EngineKind, Naca, ShapeKind, ShapeParams, Transition, TunnelConfig, Turbulence, WindTunnel};

/// Fait tourner une configuration et analyse les `avg_t` derniers temps convectifs.
pub fn run(cfg: TunnelConfig, t_end: f64, avg_t: f64) -> (Analysis, bool) {
    let mut wt = WindTunnel::new(cfg).expect("veine");
    let spt = wt.steps_per_tstar();
    let every = ((spt / 40.0).round() as usize).max(1);
    let total = (t_end * spt).round() as usize;
    let (mut t, mut cx, mut cz, mut cm) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let mut done = 0;
    while done < total {
        let n = every.min(total - done);
        wt.run(n);
        done += n;
        let sum = wt.engine.take_forces();
        let c = wt.coefficients(sum);
        let ts = wt.t_star();
        if ts >= t_end - avg_t {
            t.push(ts);
            cx.push(c.cx);
            cz.push(c.cz);
            cm.push(c.cm);
        }
    }
    let ok = wt.engine.is_finite();
    (analyze(&t, &cx, &cz, &cm), ok)
}

fn cyl(re: f64, l: usize) -> TunnelConfig {
    TunnelConfig {
        l,
        re,
        shape: ShapeKind::Cylinder,
        engine: EngineKind::Cpu { threads: 0 },
        ..Default::default()
    }
}

fn naca(re: f64, l: usize, alpha: f64) -> TunnelConfig {
    TunnelConfig {
        l,
        re,
        shape: ShapeKind::Naca,
        params: ShapeParams {
            alpha_deg: alpha,
            naca: Naca::parse("0012").unwrap(),
            ..Default::default()
        },
        engine: EngineKind::Cpu { threads: 0 },
        ..Default::default()
    }
}

#[test]
fn cylindre_re20_stationnaire() {
    // JS : Cx = 2,373 (2,05–2,09 en milieu infini, confinement de la veine).
    let (a, ok) = run(cyl(20.0, 16), 5.0 + 30.0, 10.0);
    println!("Cylindre Re = 20 : Cx = {:.3}, Cz = {:.4}", a.cx, a.cz);
    assert!(ok && a.cx > 2.0 && a.cx < 2.6 && a.cz.abs() < 0.02);
}

#[test]
fn cylindre_re100_karman() {
    // JS : Cx = 1,496, St = 0,176, Cz' = 0,375.
    let (a, ok) = run(cyl(100.0, 16), 90.0, 25.0);
    println!(
        "Cylindre Re = 100 : Cx = {:.3}, St = {:.3}, Cz' = {:.3}",
        a.cx, a.st, a.cz_amp
    );
    assert!(
        ok && a.periodic && a.cx > 1.3 && a.cx < 1.7 && a.st > 0.15 && a.st < 0.2 && a.cz_amp > 0.25 && a.cz_amp < 0.5
    );
}

#[test]
fn naca_antisymetrie() {
    // JS : Cz(±4°) = ±0,1972 à Re = 500, L = 12.
    let (p, ok1) = run(naca(500.0, 12, 4.0), 5.0 + 30.0, 10.0);
    let (m, ok2) = run(naca(500.0, 12, -4.0), 5.0 + 30.0, 10.0);
    println!(
        "NACA 0012 Re = 500 : Cz(+4°) = {:.4}, Cz(−4°) = {:.4}, Cx = {:.4}",
        p.cz, m.cz, p.cx
    );
    assert!(ok1 && ok2 && p.cz > 0.15 && p.cz < 0.26 && (p.cz + m.cz).abs() < 0.01);
}

#[test]
fn naca_re1000_traine() {
    // JS : Cx = 0,1454 (L = 16, sous-maille amortie) ; référence fine 0,119.
    let (a, ok) = run(naca(1000.0, 16, 0.0), 5.0 + 30.0, 10.0);
    println!("NACA 0012 Re = 1000 : Cx = {:.4}, Cz = {:.4}", a.cx, a.cz);
    assert!(ok && a.cx > 0.12 && a.cx < 0.17 && a.cz.abs() < 0.01);
}

#[test]
fn magnus() {
    // JS : Cz = 2,642 (≈ 2,5 d'après Mittal & Kumar 2003).
    let mut c = cyl(100.0, 12);
    c.params.spin_ratio = 1.0;
    let (a, ok) = run(c, 70.0, 25.0);
    println!("Cylindre tournant ΩR/U = 1 : Cz = {:.3}, Cx = {:.3}", a.cz, a.cx);
    assert!(ok && a.cz > 2.0 && a.cz < 3.2);
}

// Spalart–Allmaras + loi de paroi, transition imposée par la méthode des
// panneaux. Comme dans l'application web, seules la portance et la
// stabilité sont vérifiées : la traînée LBM d'un profil mince n'est pas
// convergée en maillage à ces résolutions.
#[test]
fn spalart_allmaras_transition_libre() {
    // JS : Cz = 0,491 (panneaux 0,471 ; mesures ≈ 0,44).
    let cfg = TunnelConfig {
        turbulence: Turbulence::Sa,
        transition: Transition::Free,
        ..naca(1e6, 16, 4.0)
    };
    let wt = WindTunnel::new(cfg.clone()).unwrap();
    let p = wt.panel.as_ref().expect("référence panneaux");
    assert!((p.cl - 0.471).abs() < 0.01, "panneaux Cz = {}", p.cl);
    drop(wt);
    let (a, ok) = run(cfg, 20.0, 5.0);
    println!(
        "SA, Re = 10⁶, α = 4°, transition libre : Cz = {:.4} (JS 0,491), Cx = {:.4}",
        a.cz, a.cx
    );
    assert!(ok && a.cz > 0.42 && a.cz < 0.56);
    assert!((a.cz - 0.491).abs() < 0.01);
}

#[test]
fn spalart_allmaras_re1e7() {
    // JS : Cz = 0,934 (stabilité à α = 8°).
    let cfg = TunnelConfig {
        turbulence: Turbulence::Sa,
        transition: Transition::Turbulent,
        ..naca(1e7, 16, 8.0)
    };
    let (a, ok) = run(cfg, 15.0, 5.0);
    println!("SA, Re = 10⁷, α = 8° : Cz = {:.4} (JS 0,934), Cx = {:.4}", a.cz, a.cx);
    assert!(ok && a.cz > 0.75 && a.cz < 1.1);
    assert!((a.cz - 0.934).abs() < 0.015);
}

/// Pression pariétale d'un cylindre stationnaire : point d'arrêt, symétrie, aspiration.
#[test]
fn paroi_cylindre_re40() {
    let mut wt = WindTunnel::new(cyl(40.0, 16)).unwrap();
    let spt = wt.steps_per_tstar();
    wt.run((25.0 * spt) as usize);
    let mut ws = WallSampler::new(&wt);
    let mut m = Macros::default();
    for _ in 0..20 {
        wt.run((0.25 * spt) as usize);
        wt.engine.macros(&mut m);
        ws.add(&wt, &m);
    }
    let r = ws.result();
    let near = |side: &[lbm::wall::WallPoint], s: f64| {
        side.iter()
            .min_by(|a, b| (a.s - s).abs().total_cmp(&(b.s - s).abs()))
            .unwrap()
            .cp
    };
    let cp0 = 0.5 * (near(&r.upper, 0.0) + near(&r.lower, 0.0));
    let cp_min = r.upper.iter().map(|p| p.cp).fold(f64::INFINITY, f64::min);
    let asym = (near(&r.upper, 90.0) - near(&r.lower, 90.0)).abs();
    let cf_max = r.upper.iter().map(|p| p.cf).fold(f64::NEG_INFINITY, f64::max);
    println!(
        "Cylindre Re = 40 : Cp(0°) = {cp0:.3}, Cp min = {cp_min:.3}, écart haut/bas à 90° = {asym:.4}, Cf max = {cf_max:.4}, {} mailles",
        r.upper.len() + r.lower.len()
    );
    assert!(cp0 > 0.9 && cp0 < 1.35);
    assert!(cp_min < -0.6 && cp_min > -2.0);
    assert!(asym < 0.02);
    assert!(cf_max > 0.0);
}
