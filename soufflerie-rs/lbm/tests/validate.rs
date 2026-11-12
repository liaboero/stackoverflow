//! Validation physique du moteur CPU : mêmes cas et mêmes bornes que la suite
//! de l'application web (tests/validate.mjs), valeurs JS rappelées.

use lbm::stats::{analyze, Analysis};
use lbm::{EngineKind, Naca, ShapeKind, ShapeParams, TunnelConfig, WindTunnel};

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
