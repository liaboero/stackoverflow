//! Diagnostic : instant et grille où le calcul raffiné diverge.
use lbm::cpu::Macros;
use lbm::tunnel::Engine;
use lbm::{EngineKind, Naca, ShapeKind, ShapeParams, TunnelConfig, Turbulence, WindTunnel};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let l: usize = a.get(1).and_then(|s| s.parse().ok()).unwrap_or(64);
    let refine: usize = a.get(2).and_then(|s| s.parse().ok()).unwrap_or(3);
    let half = a.get(3).map(|s| s != "fp32").unwrap_or(true);
    let sa = a.get(4).map(|s| s != "nosa").unwrap_or(true);
    let mut cfg = TunnelConfig { l, re: 1e6, shape: ShapeKind::Naca, engine: EngineKind::Gpu, gpu_half: half, refine, ..Default::default() };
    if !sa { cfg.turbulence = Turbulence::Off; }
    cfg.params = ShapeParams { alpha_deg: 4.0, naca: Naca::parse("2412").unwrap(), ..cfg.params.clone() };
    let mut wt = WindTunnel::new(cfg).unwrap();
    let spt = wt.steps_per_tstar();
    let mut m = Macros::default();
    while wt.t_star() < 12.0 {
        let dt = if wt.t_star() > 2.8 { spt / 64.0 } else { spt / 4.0 };
        wt.run(dt.max(1.0) as usize);
        let mut line = format!("t* {:5.2}", wt.t_star());
        let Engine::Refined(r) = &mut wt.engine else { return };
        let mut bad = false;
        for k in 0..r.levels.len() {
            r.level_macros(k, &mut m);
            let (mut rmin, mut rmax, mut umax) = (f32::MAX, f32::MIN, 0f32);
            let mut nan = 0;
            let mut at = 0usize;
            for i in 0..m.rho.len() {
                if !m.rho[i].is_finite() { nan += 1; continue; }
                rmin = rmin.min(m.rho[i]); rmax = rmax.max(m.rho[i]);
                let u = m.ux[i].hypot(m.uy[i]);
                if u > umax { umax = u; at = i; }
            }
            line += &format!(" | k{k} ρ {rmin:.3}..{rmax:.3} |u| {umax:.3}@({},{}) nan {nan}", at % m.nx, at / m.nx);
            bad |= nan > 0;
        }
        println!("{line}");
        if bad {
            // Où, dans la grille la plus fine ?
            for k in 0..r.levels.len() {
                r.level_macros(k, &mut m);
                if let Some(i) = m.rho.iter().position(|v| !v.is_finite()) {
                    println!("  k{k} : premier NaN en ({}, {}) sur {}×{}", i % m.nx, i / m.nx, m.nx, m.ny);
                }
            }
            break;
        }
    }
}
