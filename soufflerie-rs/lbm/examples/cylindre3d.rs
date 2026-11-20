//! Cylindre 3D à Re = 3900 (Parnaudeau et al. 2008) : Strouhal, Cd, Cl′,
//! résolution légère, sous-maille de Smagorinsky, envergure périodique 4D,
//! veine 16 × 4 × 8 D (blocage D/h = 1/8).
//!
//! cylindre3d [D=16] [Re=3900] [t=150] [moyenne=100] [fp16=0]
//!
//! `SOUFFLERIE_RANS=sa|keps|vles` : modèle RANS avec loi de paroi en plus de
//! la sous-maille (sans cette variable : sous-maille seule, comme avant).
use lbm::bodies3d::{Body3D, Shape3D};
use lbm::fwh::refine_peak;
use lbm::gpu3d::Gpu3D;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let g = |i: usize, d: f64| a.get(i).and_then(|s| s.parse().ok()).unwrap_or(d);
    let (l, re, t, avg, half) = (g(1, 16.0) as usize, g(2, 3900.0), g(3, 150.0), g(4, 100.0), g(5, 0.0) != 0.0);
    let u = 0.1f32;
    let spt = l as f64 / u as f64;
    let body = Body3D::new(Shape3D::Cylinder, l, 0.0);
    let turb = if std::env::var("SOUFFLERIE_RANS").is_ok() { lbm::gpu3d::Turbulence::rans_from_env() } else { lbm::gpu3d::Turbulence::Laminar };
    let mut sim = Gpu3D::with_turbulence(&body, u as f64 * l as f64 / re, u, half, turb).expect("GPU");
    let ramp = lbm::cpu::Ramp { steps: (4.0 * spt) as u64, u };
    let q = 0.5 * (u as f64).powi(2) * body.area;
    let clock = std::time::Instant::now();
    sim.run(((t - avg) * spt) as usize, &ramp);
    let _ = sim.take_forces();
    let dt = 0.05;
    let n = (avg / dt) as usize;
    let (mut cx, mut cz, mut ts) = (Vec::new(), Vec::new(), Vec::new());
    for k in 0..n {
        sim.run((dt * spt) as usize, &ramp);
        let f = sim.take_forces();
        cx.push(f.fx / q);
        cz.push(f.fz / q);
        ts.push((k as f64 + 0.5) * dt);
    }
    let el = clock.elapsed().as_secs_f64();
    let _ = ts;
    // Pic du spectre de Cl entre St = 0,1 et 0,4, puis affinage.
    let spec = |f: f64| lbm::fwh::dft_amp(&cz, 0.0, dt, 2.0 * std::f64::consts::PI * f).abs();
    let guess = (0..=300).map(|k| 0.1 + 0.001 * k as f64).fold((0.0, 0.0), |b, f| {
        let v = spec(f);
        if v > b.1 { (f, v) } else { b }
    });
    let st = refine_peak(&cz, dt, guess.0);
    let m = cx.iter().sum::<f64>() / n as f64;
    let mz = cz.iter().sum::<f64>() / n as f64;
    let rms = (cz.iter().map(|v| (v - mz).powi(2)).sum::<f64>() / n as f64).sqrt();
    let rmx = (cx.iter().map(|v| (v - m).powi(2)).sum::<f64>() / n as f64).sqrt();
    println!(
        "Cylindre 3D, Re {re}, D = {l} mailles, ν = {:.2e} : St = {st:.4} · Cd = {m:.4} (rms {rmx:.4}) · Cl rms = {rms:.4} · {:.0} s ({avg} D/U moyennés)",
        u as f64 * l as f64 / re,
        el
    );
}
