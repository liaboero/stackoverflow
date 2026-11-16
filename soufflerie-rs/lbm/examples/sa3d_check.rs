//! Diagnostic : Spalart–Allmaras 3D, traînée avec et sans loi de paroi.
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
    for sa in [false, true] {
        let t0 = std::time::Instant::now();
        let mut g = Gpu3D::with_turbulence(&body, u as f64 * l as f64 / re, u, true, sa).unwrap();
        let prep = t0.elapsed().as_secs_f64();
        let ramp = lbm::cpu::Ramp { steps: (4.0 * spt) as u64, u };
        let q = 0.5 * (u as f64).powi(2) * body.area;
        let mut hist = String::new();
        let (mut cx, mut cz, mut n) = (0.0, 0.0, 0);
        for k in 0..15 {
            g.run((2.0 * spt) as usize, &ramp);
            let f = g.take_forces();
            hist += &format!(" {:.3}", f.fx / q);
            if k >= 9 {
                cx += f.fx / q;
                cz += f.fz / q;
                n += 1;
            }
        }
        println!(
            "{:?} L {l} Re {re:.0e} SA {sa} : Cx {:.4} Cz {:.4} (t* 20–30) · préparation {prep:.1} s · histoire{hist}",
            shape,
            cx / n as f64,
            cz / n as f64
        );
    }
}
