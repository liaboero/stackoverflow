//! Diagnostic : aile uniforme L = 32, centre décalé d'une demi-maille (comme une grille raffinée).
use lbm::bodies3d::{Body3D, Shape3D};
use lbm::gpu3d::Gpu3D;
fn main() {
    let re: f64 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(300.0);
    let c = Body3D::new(Shape3D::Wing, 16, 8.0);
    let b = c.refined([0.0; 3], [2 * c.nx, 2 * c.ny, 2 * c.nz]);
    let (u, l) = (0.1f32, 32usize);
    let spt = l as f64 / u as f64;
    let mut g = Gpu3D::with_turbulence(&b, u as f64 * l as f64 / re, u, true, re >= 5e4).unwrap();
    let ramp = lbm::cpu::Ramp { steps: (4.0 * spt) as u64, u };
    g.run((12.0 * spt) as usize, &ramp);
    let _ = g.take_forces();
    g.run((4.0 * spt) as usize, &ramp);
    let f = g.take_forces();
    let q = 0.5 * (u as f64).powi(2) * b.area;
    println!("uniforme L 32 décalé ½ Re {re} : Cx {:.4} Cz {:.4} ({} pariétales)", f.fx / q, f.fz / q, g.wall_count);
}
