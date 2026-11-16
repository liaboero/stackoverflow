//! Diagnostic : bruit du calcul 3D (sphère) — fluctuation rapide de Cx et
//! vorticité résiduelle loin du corps.
use lbm::bodies3d::{Body3D, Shape3D};
use lbm::gpu3d::{Gpu3D, SliceField};
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let shape = Shape3D::from_id(a.get(1).map(String::as_str).unwrap_or("sphere")).unwrap();
    let l: usize = a.get(2).and_then(|s| s.parse().ok()).unwrap_or(16);
    let re: f64 = a.get(3).and_then(|s| s.parse().ok()).unwrap_or(1e6);
    let alpha: f64 = a.get(4).and_then(|s| s.parse().ok()).unwrap_or(0.0);
    let sa = re >= 5e4;
    let u = 0.1f32;
    let spt = l as f64 / u as f64;
    let body = Body3D::new(shape, l, alpha);
    let mut g = Gpu3D::with_turbulence(&body, u as f64 * l as f64 / re, u, true, sa).unwrap();
    let ramp = lbm::cpu::Ramp { steps: (4.0 * spt) as u64, u };
    let q = 0.5 * (u as f64).powi(2) * body.area;
    g.run((24.0 * spt) as usize, &ramp);
    let _ = g.take_forces();
    let mut v = Vec::new();
    for _ in 0..60 {
        g.run((0.1 * spt) as usize, &ramp);
        v.push(g.take_forces().fx / q);
    }
    let m = v.iter().sum::<f64>() / v.len() as f64;
    // Fluctuation rapide : écart entre relevés successifs.
    let hf = (v.windows(2).map(|w| (w[1] - w[0]).powi(2)).sum::<f64>() / (v.len() - 1) as f64).sqrt();
    let s = g.slice(1, body.center[1].round() as usize, SliceField::Vorticity, l as f64);
    // Vorticité (ωy·L/U) au loin : en amont du corps et sur les côtés.
    let x_body = (body.bounds.x0 - 1.5 * l as f64) as usize;
    let mut acc = (0.0, 0usize);
    for j in 0..s.h {
        for i in 2..x_body.min(s.w) {
            let w = s.data[i + s.w * j];
            if w.abs() < 1e29 {
                acc.0 += (w as f64).powi(2);
                acc.1 += 1;
            }
        }
    }
    println!(
        "{shape:?} L {l} Re {re:.0e} SA {sa} : Cx {m:.4} · fluctuation rapide {hf:.4} · vorticité amont rms {:.4}",
        (acc.0 / acc.1.max(1) as f64).sqrt()
    );
}
