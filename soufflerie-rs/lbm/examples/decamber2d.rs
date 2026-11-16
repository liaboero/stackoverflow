//! Décambrage visqueux 2D (panneaux + couche limite intégrale, couplage par
//! profil épaissi) : portance non visqueuse contre visqueuse, selon le nombre
//! d'itérations du couplage. decamber2d [NACA, défaut 0012]
use lbm::panel::{panel_solve, PanelOptions, PanelShape};

fn main() {
    let naca = lbm::shapes::Naca::parse(&std::env::args().nth(1).unwrap_or("0012".into())).unwrap();
    let shape = PanelShape::Naca(naca);
    for re in [1e6, 4e6, 6e6] {
        for a in [4.0, 8.0, 10.0] {
            let mut row = format!("Re {re:.0e} α {a:4.1} :");
            for it in [1usize, 5, 20] {
                let o = PanelOptions { iterations: it, ..PanelOptions::tripped(re, 0.0) };
                let r = panel_solve(&shape, a, &o);
                if it == 1 {
                    row += &format!(" non visqueux {:.4} ;", r.cl_inviscid);
                }
                row += &format!(" {it} it. {:.4} ({:+.1} %)", r.cl, 100.0 * (r.cl / r.cl_inviscid - 1.0));
            }
            println!("{row}");
        }
    }
}
