//! Décambrage visqueux 2D (panneaux + couche limite intégrale) : portance non
//! visqueuse, couplage faible (profil épaissi, 1 et 20 itérations) et couplage
//! simultané (Newton, à la XFOIL), traînée et transition.
//!
//! decamber2d [NACA, défaut 0012] [Re, défaut 6e6] [x_tr, défaut 0 : déclenchée
//! au bord d'attaque ; « libre » : eᴺ, N = 9] [incidences, défaut 0:15:1]
use lbm::panel::{panel_solve, Coupling, PanelOptions, PanelShape};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let naca = lbm::shapes::Naca::parse(a.get(1).map_or("0012", |s| s.as_str())).unwrap();
    let shape = PanelShape::Naca(naca);
    let re: f64 = a.get(2).map_or(6e6, |s| s.parse().unwrap());
    let base = match a.get(3).map(|s| s.as_str()) {
        Some("libre") => PanelOptions::free(re),
        Some(x) => PanelOptions::tripped(re, x.parse().unwrap()),
        None => PanelOptions::tripped(re, 0.0),
    };
    let r: Vec<f64> = a.get(4).map_or("0:15:1", |s| s.as_str()).split(':').map(|v| v.parse().unwrap()).collect();
    println!(
        "NACA {} · Re = {re:.0e} · {}",
        naca.code(),
        if base.x_trip_upper.is_some() { "transition déclenchée" } else { "transition libre (N = 9)" }
    );
    println!(" alpha | Cz inv. | faible 1 it. | faible 20 it. |  Newton  Cz/Cz_inv     Cx      Cm   x_tr ext./int.");
    let mut al = r[0];
    while al <= r[1] + 1e-9 {
        let w1 = panel_solve(&shape, al, &base);
        let w20 = panel_solve(&shape, al, &PanelOptions { iterations: 20, ..base });
        let o = PanelOptions { coupling: Coupling::Newton, ..base };
        let n = panel_solve(&shape, al, &o);
        println!(
            "{al:6.2} | {:7.4} | {:12.4} | {:13.4} | {:7.4}  {:9.3}  {:8.5} {:7.4}  {:5.3} / {:5.3}",
            w1.cl_inviscid,
            w1.cl,
            w20.cl,
            n.cl,
            if n.cl_inviscid.abs() > 0.05 { n.cl / n.cl_inviscid } else { f64::NAN },
            n.cd,
            n.cm,
            n.transition_upper,
            n.transition_lower
        );
        al += r[2];
    }
}
