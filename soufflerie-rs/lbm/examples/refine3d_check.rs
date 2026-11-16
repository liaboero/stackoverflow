//! Diagnostic : raffinement local 3D contre grille uniforme.
use lbm::bodies3d::{Body3D, Shape3D};
use lbm::gpu3d::Gpu3D;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let shape = Shape3D::from_id(a.get(1).map(String::as_str).unwrap_or("sphere")).unwrap();
    let l: usize = a.get(2).and_then(|s| s.parse().ok()).unwrap_or(8);
    let re: f64 = a.get(3).and_then(|s| s.parse().ok()).unwrap_or(100.0);
    let alpha: f64 = a.get(4).and_then(|s| s.parse().ok()).unwrap_or(0.0);
    let depth: usize = a.get(5).and_then(|s| s.parse().ok()).unwrap_or(1);
    let t_end: f64 = a.get(6).and_then(|s| s.parse().ok()).unwrap_or(30.0);
    let u = 0.1f32;
    let spt = l as f64 / u as f64;
    let body = Body3D::new(shape, l, alpha);
    let t0 = std::time::Instant::now();
    let mut g = Gpu3D::with_refinement(&body, u as f64 * l as f64 / re, u, true, re >= 5e4, depth).unwrap();
    let prep = t0.elapsed().as_secs_f64();
    let ramp = lbm::cpu::Ramp { steps: (4.0 * spt) as u64, u };
    let q = 0.5 * (u as f64).powi(2) * body.area;
    let t1 = std::time::Instant::now();
    let mut hist = String::new();
    let (mut cx, mut cz, mut n) = (0.0, 0.0, 0);
    let chunks = (t_end / 2.0) as usize;
    for k in 0..chunks {
        g.run((2.0 * spt) as usize, &ramp);
        let f = g.take_forces();
        hist += &format!(" {:.3}", f.fx / q);
        if k >= chunks * 2 / 3 {
            cx += f.fx / q;
            cz += f.fz / q;
            n += 1;
        }
    }
    let fine = g.finest().map(|(d, w, s)| format!("grille fine ×{s} {}×{}×{} ({w} pariétales)", d[0], d[1], d[2])).unwrap_or_default();
    println!(
        "{shape:?} L {l} raff {depth} Re {re:.0e} : Cx {:.4} Cz {:.4} · {:.0} s (préparation {prep:.1} s) · {fine} ·{hist}",
        cx / n as f64,
        cz / n as f64,
        t1.elapsed().as_secs_f64()
    );
}
