//! Diagnostic : traînée de frottement (∫Cf dx) et de pression, LBM et panneaux.
use lbm::cpu::Macros;
use lbm::panel::section_coefficients;
use lbm::wall::WallSampler;
use lbm::{EngineKind, Naca, ShapeKind, ShapeParams, Transition, TunnelConfig, WindTunnel};

fn integ(x: &[f64], y: &[f64]) -> f64 {
    (1..x.len()).map(|k| 0.5 * (y[k] + y[k - 1]) * (x[k] - x[k - 1])).sum()
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let l: usize = a.get(1).and_then(|s| s.parse().ok()).unwrap_or(128);
    let turb = a.get(2).map(|s| s == "turb").unwrap_or(false);
    let refine: usize = a.get(3).and_then(|s| s.parse().ok()).unwrap_or(0);
    let t_end: f64 = a.get(4).and_then(|s| s.parse().ok()).unwrap_or(26.0);
    let mut cfg = TunnelConfig { l, re: 1e6, shape: ShapeKind::Naca, engine: EngineKind::Gpu, gpu_half: true, refine, ..Default::default() };
    if turb { cfg.transition = Transition::Turbulent; }
    cfg.params = ShapeParams { alpha_deg: 4.0, naca: Naca::parse("2412").unwrap(), ..cfg.params.clone() };
    let mut wt = WindTunnel::new(cfg).unwrap();
    // Longueur d'amortissement pariétal de la sous-maille imposée (fraction de corde).
    if let Some(ac) = a.get(5).and_then(|s| s.parse::<f64>().ok()) {
        let lw = wt.wall_l() as f64;
        let l0 = wt.cfg.l as f64;
        match &mut wt.engine {
            lbm::tunnel::Engine::Refined(r) => r.set_les_wall_lengths(|k| (ac * l0 * (1u64 << k) as f64).max(1.0)),
            lbm::tunnel::Engine::Gpu(g) => g.set_les_wall_length((ac * lw).max(1.0)),
            _ => {}
        }
        eprintln!("amortissement A = {ac} c = {:.1} mailles près du corps", ac * lw);
    }
    let spt = wt.steps_per_tstar();
    eprintln!("{} · L près du corps {}", wt.engine.name(), wt.wall_l());
    let t0 = std::time::Instant::now();
    let mut hist = String::new();
    while wt.t_star() < t_end - 8.0 {
        wt.run((2.0 * spt) as usize);
        let f = wt.engine.take_forces();
        let c = wt.coefficients(f);
        hist += &format!(" {:.0}:{:.4}", wt.t_star(), c.cx);
    }
    eprintln!("Cx par 2 t* :{hist}");
    let mut ws = WallSampler::new(&wt);
    let mut m = Macros::default();
    let (mut fx, mut st) = (0.0, 0u64);
    while wt.t_star() < t_end {
        wt.run((spt / 4.0) as usize);
        let s = wt.engine.take_forces();
        fx += s.fx; st += s.steps;
        wt.engine.macros(&mut m);
        ws.add(&mut wt, &m);
    }
    eprintln!("{:.0} s", t0.elapsed().as_secs_f64());
    let cx = wt.coefficients(lbm::cpu::ForceSum { fx, fy: 0.0, mz: 0.0, steps: st }).cx;
    let w = ws.result();
    let f = |v: &Vec<lbm::wall::WallPoint>| integ(&v.iter().map(|p| p.s).collect::<Vec<_>>(), &v.iter().map(|p| p.cf).collect::<Vec<_>>());
    let (fu, fl) = (f(&w.upper), f(&w.lower));
    let ca = 4f64.to_radians().cos();
    let lbm_f = (fu + fl) * ca;
    let r = wt.panel.as_ref().unwrap();
    let finf = |f: &lbm::panel::Friction| -> f64 {
        let mut v: Vec<(f64, f64)> = f.x.iter().zip(f.cf.iter().zip(&f.ue)).map(|(&x, (&c, &u))| (x, c * u * u)).collect();
        v.sort_by(|a, b| a.0.total_cmp(&b.0));
        (1..v.len()).map(|k| 0.5 * (v[k].1 + v[k - 1].1) * (v[k].0 - v[k - 1].0)).sum()
    };
    let (pu, pl) = (finf(&r.friction_upper), finf(&r.friction_lower));
    let k = section_coefficients(r);
    let pan_f = (pu + pl) * ca;
    println!("L {l} {}: LBM Cx {cx:.5} = frottement {lbm_f:.5} (extrados {fu:.5}, intrados {fl:.5}) + pression {:.5}", if turb {"turb"} else {"libre"}, cx - lbm_f);
    println!("        panneaux Cx {:.5} = frottement {pan_f:.5} (extrados {pu:.5}, intrados {pl:.5}) + pression {:.5}", k.cd, k.cd - pan_f);
    // Cf le long de l'extrados (quelques points).
    for xs in [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9] {
        let near = |v: &Vec<lbm::wall::WallPoint>| v.iter().filter(|p| (p.s - xs).abs() < 0.02).map(|p| p.cf).sum::<f64>() / v.iter().filter(|p| (p.s - xs).abs() < 0.02).count().max(1) as f64;
        let pn = |x: &[f64], c: &[f64]| { let i = x.iter().enumerate().min_by(|a, b| (a.1 - xs).abs().total_cmp(&(b.1 - xs).abs())).unwrap().0; c[i] };
        let ue2 = |f: &lbm::panel::Friction| f.cf.iter().zip(&f.ue).map(|(c, u)| c * u * u).collect::<Vec<_>>();
        let (cu, cl) = (ue2(&r.friction_upper), ue2(&r.friction_lower));
        println!("   x/c {xs:.2} : Cf extrados LBM {:.5} panneaux {:.5} · intrados LBM {:.5} panneaux {:.5}", near(&w.upper), pn(&r.friction_upper.x, &cu), near(&w.lower), pn(&r.friction_lower.x, &cl));
    }
}
