//! Balayage d'incidence jusqu'au décrochage, 2D (panneaux + couche limite
//! intégrale, couplage simultané Newton), contre Ladson (NASA TM-4074,
//! transition déclenchée, grains 80 / 120 / 180) et TMR (Spalart–Allmaras).
//!
//! decrochage2d [NACA, défaut 0012] [Re, défaut 6e6] [incidences a:b:pas,
//! défaut 0:20:1] [x_tr, défaut 0 : déclenchée au bord d'attaque ; « libre »]
//!
//! Colonnes : α, Cz, Cx, Cm (quart de corde, cabreur positif), Cz/Cz_inv,
//! abscisse du décollement turbulent (premier Cf < 0, x/c) sur l'extrados et
//! l'intrados (1 : attaché), transition extrados ; `NaN` : Newton n'a pas
//! convergé. Puis Cz max et son incidence (sur les points convergés).
use lbm::panel::{panel_solve, Coupling, PanelOptions, PanelShape};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let naca = lbm::shapes::Naca::parse(a.get(1).map_or("0012", |s| s.as_str())).unwrap();
    let shape = PanelShape::Naca(naca);
    let re: f64 = a.get(2).map_or(6e6, |s| s.parse().unwrap());
    let r: Vec<f64> = a.get(3).map_or("0:20:1", |s| s.as_str()).split(':').map(|v| v.parse().unwrap()).collect();
    let base = match a.get(4).map(|s| s.as_str()) {
        Some("libre") => PanelOptions::free(re),
        Some(x) => PanelOptions::tripped(re, x.parse().unwrap()),
        None => PanelOptions::tripped(re, 0.0),
    };
    let o = PanelOptions { coupling: Coupling::Newton, ..base };
    println!("NACA {} · Re = {re:.0e} · couplage Newton", naca.code());
    println!("  alpha |   Cz    |    Cx    |    Cm    | Cz/Cz_inv | décoll. ext. / int. | x_tr ext.");
    let (mut best, mut n) = ((f64::NAN, f64::NEG_INFINITY), 0);
    loop {
        let al = r[0] + n as f64 * r[2];
        if al > r[1] + 1e-9 {
            break;
        }
        n += 1;
        let t = std::time::Instant::now();
        let s = panel_solve(&shape, al, &o);
        let sep = |x: f64| if x.is_finite() { x } else { 1.0 };
        println!(
            "{al:7.2} | {:7.4} | {:8.5} | {:8.5} | {:9.3} | {:8.3} / {:8.3} | {:6.3}   ({:.1} s)",
            s.cl,
            s.cd,
            s.cm,
            if s.cl_inviscid.abs() > 0.05 { s.cl / s.cl_inviscid } else { f64::NAN },
            sep(s.separation_upper),
            sep(s.separation_lower),
            s.transition_upper,
            t.elapsed().as_secs_f64()
        );
        if s.cl.is_finite() && s.cl > best.1 {
            best = (al, s.cl);
        }
    }
    println!("Cz max {:.4} à α = {}°", best.1, best.0);
}
