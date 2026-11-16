//! Diagnostic 2D : historique de Cz et Cp/Cf pariétaux moyens contre les panneaux.
//! cargo run --release --example cp2d -- <L> <raff> <alpha> <re> <t> <avg> <csv>
use lbm::cpu::Macros;
use lbm::wall::WallSampler;
use lbm::{EngineKind, Naca, ShapeKind, ShapeParams, Transition, TunnelConfig, Turbulence, WindTunnel};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let g = |i: usize, d: f64| a.get(i).and_then(|v| v.parse().ok()).unwrap_or(d);
    let (l, raff, alpha, re, t_end, avg) = (g(1, 32.0) as usize, g(2, 2.0) as usize, g(3, 10.0), g(4, 6e6), g(5, 30.0), g(6, 10.0));
    let csv = a.get(7).cloned().unwrap_or_else(|| "cp2d.csv".into());
    let cfg = TunnelConfig {
        l,
        re,
        u_in: 0.1,
        shape: ShapeKind::Naca,
        params: ShapeParams {
            alpha_deg: alpha,
            naca: Naca::parse("0012").unwrap(),
            thickness: 0.4,
            spin_ratio: 0.0,
        },
        les: lbm::cpu::LesMode::Auto,
        collision: lbm::cpu::Collision::Rr,
        bouzidi: true,
        engine: EngineKind::Gpu,
        turbulence: if re < 5e4 { Turbulence::Off } else { Turbulence::Sa },
        transition: Transition::Turbulent,
        gpu_half: false,
        refine: raff,
    };
    let mut wt = WindTunnel::new(cfg).unwrap();
    let spt = wt.steps_per_tstar();
    let mut sampler = WallSampler::new(&wt);
    let mut m = Macros::default();
    let every = (spt / 4.0).round() as usize;
    let total = (t_end * spt) as usize;
    let mut done = 0;
    let (mut sx, mut sz, mut n) = (0.0, 0.0, 0.0);
    while done < total {
        wt.run(every);
        done += every;
        let s = wt.engine.take_forces();
        let c = wt.coefficients(s);
        let ts = done as f64 / spt;
        if (ts * 4.0).round() as usize % 8 == 0 {
            eprintln!("t* {ts:5.1} Cx {:.4} Cz {:.4}", c.cx, c.cz);
        }
        if ts > t_end - avg {
            wt.engine.macros(&mut m);
            sampler.add(&mut wt, &m);
            sx += c.cx;
            sz += c.cz;
            n += 1.0;
        }
    }
    println!("L {l} raff {raff} : Cx {:.4} Cz {:.4}", sx / n, sz / n);
    if let Some(p) = &wt.panel {
        println!("panneaux : Cl {:.4} (non visqueux {:.4}) Cd {:.4}", p.cl, p.cl_inviscid, p.cd);
    }
    let pr = sampler.result();
    let mut out = String::from("face,s,cp,cf\n");
    for (f, v) in [("ext", &pr.upper), ("int", &pr.lower)] {
        for p in v.iter() {
            out += &format!("{f},{:.5},{:.5},{:.6}\n", p.s, p.cp, p.cf);
        }
    }
    if let Some(p) = &wt.panel {
        for k in 0..p.cp.len() {
            let f = if k < p.i_le { "pan_int" } else { "pan_ext" };
            out += &format!("{f},{:.5},{:.5},0\n", p.cp_x[k], p.cp[k]);
        }
    }
    std::fs::write(&csv, out).unwrap();
    profiles(&mut wt);
}

/// Profils de couche limite (grille la plus fine) à l'extrados : u_t et ν_t/ν le long de la normale.
fn profiles(wt: &mut WindTunnel) {
    let mut m = Macros::default();
    wt.engine.wall_macros(&mut m);
    let st = match &mut wt.engine {
        lbm::tunnel::Engine::Refined(s) => s.levels.last_mut().unwrap().sa_state(),
        lbm::tunnel::Engine::Gpu(s) => s.sa_state(),
        _ => Vec::new(),
    };
    let geo = wt.engine.wall_geometry();
    let body = wt.wall_body().clone();
    let nu = wt.cfg.u_in as f64 * wt.wall_l() as f64 / wt.cfg.re;
    let loc = &body.outline_local;
    let nl = loc.len();
    let (mut i_le, mut xmin, mut xmax) = (0, f64::INFINITY, f64::NEG_INFINITY);
    for (k, p) in loc.iter().enumerate() {
        if p[0] < xmin {
            xmin = p[0];
            i_le = k;
        }
        xmax = xmax.max(p[0]);
    }
    let nx = m.nx;
    let at = |f: &dyn Fn(usize) -> f64, x: f64, y: f64| {
        let (i, j) = (x.floor() as usize, y.floor() as usize);
        let (fx, fy) = (x - i as f64, y - j as f64);
        let v = |a: usize, b: usize| f(b * nx + a);
        (1.0 - fx) * (1.0 - fy) * v(i, j) + fx * (1.0 - fy) * v(i + 1, j) + (1.0 - fx) * fy * v(i, j + 1) + fx * fy * v(i + 1, j + 1)
    };
    for target in [0.05, 0.3, 0.6, 0.9, -0.3] {
        let lower = target < 0.0;
        let target: f64 = f64::abs(target);
        let mut best = (f64::INFINITY, 0usize);
        for &idx in &geo.wall_cells {
            let pos = geo.near_pos[idx] as f64;
            if pos < 0.0 {
                continue;
            }
            let k = (pos.floor() as usize).min(nl - 1);
            if (k >= i_le) != lower {
                continue;
            }
            let t = pos - k as f64;
            let bx = loc[k][0] + t * (loc[(k + 1) % nl][0] - loc[k][0]);
            let sc = (bx - xmin) / (xmax - xmin);
            if (sc - target).abs() < best.0 {
                best = ((sc - target).abs(), idx);
            }
        }
        let idx = best.1;
        let (px, py) = ((idx % nx) as f64, (idx / nx) as f64);
        let (nnx, nny) = (geo.normal[2 * idx] as f64, geo.normal[2 * idx + 1] as f64);
        let d = geo.dist[idx] as f64;
        let (wx, wy) = (px - d * nnx, py - d * nny);
        let (tx, ty) = (-nny, nnx);
        println!("x/c = {target} {} (maille à {d:.2} de la paroi), L fin = {}", if lower { "intrados" } else { "extrados" }, wt.wall_l());
        let mut line = String::new();
        let mut uprev = 0.0;
        let mut ue = 0.0;
        let mut prof = Vec::new();
        for k in 1..=60 {
            let y = 0.5 * k as f64;
            let (x, yy) = (wx + y * nnx, wy + y * nny);
            let ux = at(&|i| m.ux[i] as f64, x, yy);
            let uy = at(&|i| m.uy[i] as f64, x, yy);
            let ut = (ux * tx + uy * ty).abs();
            let nt = if st.is_empty() { 0.0 } else { at(&|i| lbm::sa::nu_turb(st[i][3], (1.0 / nu) as f32) as f64, x, yy) / nu };
            prof.push((y, ut, nt));
            if k % 2 == 0 && k <= 40 {
                line += &format!(" {y:.0}:{:.3}/{:.0}", ut / wt.cfg.u_in as f64, nt);
            }
            if ut >= uprev {
                ue = ut;
            }
            uprev = ut;
        }
        let ue = prof.iter().map(|p| p.1).fold(0.0, f64::max).max(ue);
        // δ* et θ jusqu'au maximum de vitesse (u/ue).
        let kmax = prof.iter().position(|p| p.1 >= 0.995 * ue).unwrap_or(prof.len() - 1);
        let (mut ds, mut th) = (0.0, 0.0);
        let (mut yprev, mut fprev, mut gprev) = (0.0, 1.0, 0.0);
        for p in &prof[..=kmax] {
            let r = p.1 / ue;
            ds += 0.5 * (p.0 - yprev) * (fprev + 1.0 - r);
            th += 0.5 * (p.0 - yprev) * (gprev + r * (1.0 - r));
            yprev = p.0;
            fprev = 1.0 - r;
            gprev = r * (1.0 - r);
        }
        let l = wt.wall_l() as f64;
        println!("  ue/U {:.3}  δ*/c {:.4}  θ/c {:.4}  H {:.2}  δ99/c {:.4}", ue / wt.cfg.u_in as f64, ds / l, th / l, ds / th.max(1e-9), prof[kmax].0 / l);
        println!("  y:u/U/νt/ν{line}");
    }
}
