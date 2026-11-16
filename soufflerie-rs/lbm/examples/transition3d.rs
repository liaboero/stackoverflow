//! Diagnostic : transition imposée en 3D (sphère, aile).
use lbm::bodies3d::{Body3D, Shape3D};
use lbm::gpu3d::Gpu3D;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let shape = Shape3D::from_id(a.get(1).map(String::as_str).unwrap_or("sphere")).unwrap();
    let l: usize = a.get(2).and_then(|s| s.parse().ok()).unwrap_or(16);
    let re: f64 = a.get(3).and_then(|s| s.parse().ok()).unwrap_or(1e5);
    let alpha: f64 = a.get(4).and_then(|s| s.parse().ok()).unwrap_or(0.0);
    let u = 0.1f32;
    let spt = l as f64 / u as f64;
    let body = Body3D::new(shape, l, alpha);
    for xtr in [0.0, 0.3, 0.6, 1.0] {
        let mut g = Gpu3D::with_turbulence(&body, u as f64 * l as f64 / re, u, true, true).unwrap();
        g.set_transition(xtr);
        let ramp = lbm::cpu::Ramp { steps: (4.0 * spt) as u64, u };
        g.run((22.0 * spt) as usize, &ramp);
        let _ = g.take_forces();
        g.run((8.0 * spt) as usize, &ramp);
        let f = g.take_forces();
        let q = 0.5 * (u as f64).powi(2) * body.area;
        println!("{shape:?} L {l} Re {re:.0e} xtr {xtr} : Cx {:.4} Cz {:.4}", f.fx / q, f.fz / q);
    }
}
