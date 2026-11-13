//! Aile finie : mêmes cas et mêmes valeurs que l'application web (js/vlm.js).
use lbm::wing::{solve_wing, wing_polar, PolarTable, SectionCoef, Sections, Wing, WingParams};
use lbm::Naca;

fn lin(_: usize, a: f64) -> SectionCoef {
    SectionCoef {
        cl: 2.0 * std::f64::consts::PI * a.to_radians(),
        cd: 0.0,
        cm: 0.0,
    }
}

fn wing(ar: f64, taper: f64) -> Wing {
    Wing::new(WingParams {
        ar,
        taper,
        ..Default::default()
    })
}

#[test]
fn weissinger_lineaire() {
    let w = wing(6.0, 1.0);
    let r = solve_wing(&w, 4.0, &lin, None);
    let slope = (solve_wing(&w, 5.0, &lin, None).cl - solve_wing(&w, 3.0, &lin, None).cl) / 2.0 * 180.0 / std::f64::consts::PI;
    println!("A = 6 : CL = {:.5} (JS 0,29376), e = {:.4} (JS 0,9955), CLα = {slope:.3}/rad", r.cl, r.e);
    assert!((r.cl - 0.293_759).abs() < 1e-5 && (r.e - 0.995_485).abs() < 1e-4 && r.converged);
    assert!(slope > 4.0 && slope < 4.5);
    let t = solve_wing(&wing(8.0, 0.4), 4.0, &lin, None);
    println!("A = 8, λ = 0,4 : CL = {:.5} (JS 0,33134), e = {:.4} (JS 1,0067)", t.cl, t.e);
    assert!((t.cl - 0.331_336).abs() < 1e-5 && (t.e - 1.006_713).abs() < 1e-4);
}

#[test]
fn fleche_diedre_vrillage() {
    let w = Wing::new(WingParams {
        ar: 6.0,
        taper: 0.5,
        sweep: 30.0,
        dihedral: 5.0,
        twist: -3.0,
        n: 40,
    });
    let r = solve_wing(&w, 4.0, &lin, None);
    println!(
        "Flèche 30°, dièdre 5°, vrillage −3° : CL = {:.5} (JS 0,19494), CDi = {:.6} (JS 0,002062), Cm = {:.5} (JS 0,01408)",
        r.cl, r.cdi, r.cm
    );
    assert!((r.cl - 0.194_942).abs() < 1e-5 && (r.cdi - 0.002_0618).abs() < 2e-6 && (r.cm - 0.014_079).abs() < 1e-5);
}

#[test]
fn naca_2412_non_lineaire() {
    let tab = PolarTable::naca(Naca::parse("2412").unwrap(), 1e6, false);
    let w = wing(8.0, 0.5);
    // JS : une seule table (Re = 10⁶) pour toutes les bandes.
    let section = |_: usize, a: f64| tab.at(a);
    let js = [
        (4.0, 0.505_600, 0.010_112, 0.007_106),
        (12.0, 1.123_585, 0.050_786, 0.012_457),
        (16.0, 1.365_944, 0.076_099, 0.017_588),
    ];
    for (a, cl, cdi, cdv) in js {
        let r = solve_wing(&w, a, &section, None);
        println!(
            "α = {a}° : CL = {:.4} (JS {cl:.4}), CDi = {:.5} (JS {cdi:.5}), CDp = {:.5} (JS {cdv:.5}), {} itérations",
            r.cl, r.cdi, r.cdv, r.iterations
        );
        assert!(r.converged && (r.cl - cl).abs() < 2e-3 && (r.cdi - cdi).abs() < 2e-4 && (r.cdv - cdv).abs() < 2e-4);
    }
}

#[test]
fn polaire_effilee() {
    // Deux tables (emplanture, bout) et continuation en incidence.
    let w = wing(8.0, 0.5);
    let (rr, rt) = Sections::reynolds(&w, 1e6);
    let naca = Naca::parse("2412").unwrap();
    let sec = Sections::new(&w, 1e6, PolarTable::naca(naca, rr, false), rt.map(|r| PolarTable::naca(naca, r, false)));
    let alphas: Vec<f64> = (-6..=20).map(|a| a as f64).collect();
    let p = wing_polar(&w, &sec, &alphas);
    let best = p.iter().max_by(|a, b| a.cl.total_cmp(&b.cl)).unwrap();
    println!("CL max = {:.3} à {}°", best.cl, best.alpha);
    assert!(p.iter().filter(|r| r.alpha <= 12.0).all(|r| r.converged));
    assert!(best.cl > 1.2 && best.cl < 1.6);
}
