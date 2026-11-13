//! Méthode des panneaux : mêmes cas et mêmes valeurs que l'application web.
use lbm::panel::{airfoil_points, panel_solve, potential_flow, section_coefficients, PanelOptions, PanelShape};
use lbm::Naca;

fn naca(code: &str) -> PanelShape {
    PanelShape::Naca(Naca::parse(code).unwrap())
}

#[test]
fn non_visqueux() {
    let f = potential_flow(&airfoil_points(&naca("0012"), 80), 4f64.to_radians());
    println!("NACA 0012, α = 4°, non visqueux : Cz = {:.4} (JS 0,4797)", f.cl);
    assert!((f.cl - 0.4797).abs() < 5e-4);
}

#[test]
fn turbulent_6e6() {
    let o = PanelOptions::tripped(6e6, 0.0);
    let a = panel_solve(&naca("0012"), 0.0, &o);
    let b = panel_solve(&naca("0012"), 10.0, &o);
    println!(
        "Re = 6·10⁶ turbulent : Cx(0°) = {:.5} (JS 0,0082) ; Cz, Cx(10°) = {:.3}, {:.5} (JS 1,168, 0,0129)",
        a.cd, b.cl, b.cd
    );
    assert!((a.cd - 0.0082).abs() < 2e-4 && (b.cl - 1.168).abs() < 5e-3 && (b.cd - 0.0129).abs() < 3e-4);
}

#[test]
fn transition_libre() {
    let r = panel_solve(&naca("0012"), 0.0, &PanelOptions::free(1e6));
    println!(
        "Re = 10⁶ libre : x_tr = {:.3}, Cx = {:.5} (JS 0,607, 0,0062)",
        r.transition_upper, r.cd
    );
    assert!((r.transition_upper - 0.607).abs() < 0.01 && (r.cd - 0.0062).abs() < 2e-4);
    let c = panel_solve(&naca("2412"), 0.0, &PanelOptions::free(3e6));
    println!(
        "NACA 2412, Re = 3·10⁶ : Cz, Cx, Cm = {:.3}, {:.4}, {:.3} (JS 0,253, 0,0060, −0,055)",
        c.cl, c.cd, c.cm
    );
    assert!((c.cl - 0.253).abs() < 5e-3 && (c.cd - 0.006).abs() < 3e-4);
}

#[test]
fn kirchhoff() {
    // JS : NACA 0012, Re = 10⁶ : Cz = 1,40 à 16° (décollé), 1,13 à 20°.
    let o = PanelOptions::free(1e6);
    let a = section_coefficients(&panel_solve(&naca("0012"), 16.0, &o));
    let b = section_coefficients(&panel_solve(&naca("0012"), 20.0, &o));
    println!("Kirchhoff : Cz(16°) = {:.3}, Cz(20°) = {:.3}", a.cl, b.cl);
    assert!(a.kirchhoff && (a.cl - 1.40).abs() < 0.02 && (b.cl - 1.13).abs() < 0.02);
}
