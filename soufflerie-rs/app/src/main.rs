//! Soufflerie de poche, version native.
//!
//! Sans argument : interface graphique. Sous-commandes :
//!   soufflerie run    [options]          un calcul, moyennes et Strouhal
//!   soufflerie polar  --alpha a:b:pas    polaire (incidences calculées en parallèle), CSV
//!   soufflerie bench  [--l 32]           vitesse des moteurs (millions de mailles par seconde)
//!
//! Options : --shape cylinder|square|ellipse|plate|naca|triangle  --naca 2412
//!           --alpha 4  --re 1000  --l 32  --t 40  --avg 15  --engine cpu|gpu
//!           --threads 0  --spin 1  --thickness 0.4  --collision rr|bgk
//!           --les auto|on|off  --u 0.1

use std::collections::HashMap;
use std::time::Instant;

use lbm::cpu::{Collision, LesMode};
use lbm::stats::{analyze, Analysis};
use lbm::{EngineKind, Naca, ShapeKind, ShapeParams, TunnelConfig, WindTunnel};
use rayon::prelude::*;

#[cfg(feature = "gui")]
mod gui;
#[cfg(feature = "gui")]
mod worker;

fn usage() -> ! {
    eprintln!(
        "Soufflerie de poche (native)\n\n\
         soufflerie                       interface graphique\n\
         soufflerie run   [options]       un calcul\n\
         soufflerie polar --alpha -4:12:2 [options]   polaire, CSV\n\
         soufflerie bench [--l 32] [--steps 400]       vitesse des moteurs\n\n\
         Options : --shape naca --naca 2412 --alpha 4 --re 1000 --l 32 --t 40 --avg 15\n\
         \x20         --engine cpu|gpu --threads 0 --spin 0 --thickness 0.4\n\
         \x20         --collision rr|bgk --les auto|on|off --u 0.1"
    );
    std::process::exit(2)
}

struct Args {
    cmd: String,
    opts: HashMap<String, String>,
}

impl Args {
    fn parse() -> Args {
        let mut it = std::env::args().skip(1).peekable();
        let cmd = match it.peek() {
            Some(s) if !s.starts_with("--") => it.next().unwrap(),
            _ => String::new(),
        };
        let mut opts = HashMap::new();
        while let Some(k) = it.next() {
            let Some(key) = k.strip_prefix("--") else { usage() };
            let v = it.next().unwrap_or_else(|| usage());
            opts.insert(key.to_string(), v);
        }
        Args { cmd, opts }
    }

    fn get<T: std::str::FromStr>(&self, k: &str, d: T) -> T {
        match self.opts.get(k) {
            Some(v) => v.parse().unwrap_or_else(|_| {
                eprintln!("valeur invalide pour --{k} : {v}");
                usage()
            }),
            None => d,
        }
    }

    fn config(&self) -> TunnelConfig {
        let shape = ShapeKind::from_id(&self.get("shape", "cylinder".to_string())).unwrap_or_else(|| usage());
        let naca = Naca::parse(&self.get("naca", "0012".to_string())).unwrap_or_else(|| usage());
        let engine = match self.get("engine", "cpu".to_string()).as_str() {
            "gpu" => EngineKind::Gpu,
            _ => EngineKind::Cpu {
                threads: self.get("threads", 0usize),
            },
        };
        TunnelConfig {
            l: self.get("l", 16usize),
            re: self.get("re", 100.0f64),
            u_in: self.get("u", 0.1f32),
            shape,
            params: ShapeParams {
                // Polaire : --alpha est une plage a:b:pas, lue à part.
                alpha_deg: self.opts.get("alpha").and_then(|v| v.parse().ok()).unwrap_or(0.0),
                naca,
                thickness: self.get("thickness", 0.4f64),
                spin_ratio: self.get("spin", 0.0f64),
            },
            les: match self.get("les", "auto".to_string()).as_str() {
                "on" => LesMode::On,
                "off" => LesMode::Off,
                _ => LesMode::Auto,
            },
            collision: if self.get("collision", "rr".to_string()) == "bgk" {
                Collision::Bgk
            } else {
                Collision::Rr
            },
            bouzidi: true,
            engine,
        }
    }
}

/// Calcul jusqu'à `t_end` temps convectifs ; analyse des `avg` derniers.
fn simulate(cfg: TunnelConfig, t_end: f64, avg: f64, quiet: bool) -> Result<(Analysis, f64, bool), String> {
    let mut wt = WindTunnel::new(cfg)?;
    let spt = wt.steps_per_tstar();
    let every = ((spt / 40.0).round() as usize).max(1);
    let total = (t_end * spt).round() as usize;
    let (mut t, mut cx, mut cz, mut cm) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let t0 = Instant::now();
    let mut done = 0;
    let mut next_print = 5.0;
    while done < total {
        let n = every.min(total - done);
        wt.run(n);
        done += n;
        let s = wt.engine.take_forces();
        let c = wt.coefficients(s);
        let ts = wt.t_star();
        if ts >= t_end - avg {
            t.push(ts);
            cx.push(c.cx);
            cz.push(c.cz);
            cm.push(c.cm);
        }
        if !quiet && ts >= next_print {
            eprintln!("  t* = {ts:5.1}  Cx = {:8.4}  Cz = {:8.4}", c.cx, c.cz);
            next_print += 5.0;
        }
        if !c.cx.is_finite() {
            break;
        }
    }
    let secs = t0.elapsed().as_secs_f64();
    let mlups = (wt.nx * wt.ny) as f64 * done as f64 / secs / 1e6;
    let ok = wt.engine.is_finite();
    Ok((analyze(&t, &cx, &cz, &cm), mlups, ok))
}

fn cmd_run(a: &Args) {
    let cfg = a.config();
    let wt = WindTunnel::new(cfg.clone()).unwrap_or_else(|e| {
        eprintln!("{e}");
        std::process::exit(1)
    });
    println!(
        "{} · Re = {} · L = {} mailles · grille {}×{} · {} · τ = {:.4} · sous-maille {}",
        wt.body.label(),
        cfg.re,
        cfg.l,
        wt.nx,
        wt.ny,
        wt.engine.name(),
        wt.engine.tau0(),
        if wt.engine.les_active() { "active" } else { "inactive" }
    );
    drop(wt);
    let t_end = a.get("t", 40.0f64);
    let avg = a.get("avg", 15.0f64);
    match simulate(cfg, t_end, avg, false) {
        Ok((r, mlups, ok)) => {
            println!(
                "Cx = {:.4}  Cz = {:.4}  Cm = {:.4}{}{}",
                r.cx,
                r.cz,
                r.cm,
                if r.periodic {
                    format!("  St = {:.4} ({} périodes, Cz' = ±{:.3})", r.st, r.periods, r.cz_amp)
                } else {
                    String::new()
                },
                if ok { String::new() } else { "  [DIVERGENCE]".into() }
            );
            println!("{mlups:.0} millions de mailles par seconde");
        }
        Err(e) => eprintln!("{e}"),
    }
}

fn cmd_polar(a: &Args) {
    let spec = a.get("alpha", "-4:12:2".to_string());
    let p: Vec<f64> = spec.split(':').filter_map(|s| s.parse().ok()).collect();
    if p.len() != 3 || p[2] <= 0.0 {
        usage()
    }
    let mut alphas = Vec::new();
    let mut x = p[0];
    while x <= p[1] + 1e-9 {
        alphas.push(x);
        x += p[2];
    }
    let base = a.config();
    let t_end = a.get("t", 40.0f64);
    let avg = a.get("avg", 15.0f64);
    let gpu = base.engine == EngineKind::Gpu;
    eprintln!(
        "{} incidences, {}",
        alphas.len(),
        if gpu {
            "GPU, l'une après l'autre"
        } else {
            "en parallèle sur les cœurs"
        }
    );
    let t0 = Instant::now();
    let run = |alpha: f64| {
        let mut cfg = base.clone();
        cfg.params.alpha_deg = alpha;
        if !gpu {
            // Un fil par incidence : les incidences se partagent les cœurs.
            cfg.engine = EngineKind::Cpu { threads: 1 };
        }
        (alpha, simulate(cfg, t_end, avg, true))
    };
    let res: Vec<_> = if gpu {
        alphas.iter().map(|&a| run(a)).collect()
    } else {
        alphas.par_iter().map(|&a| run(a)).collect()
    };
    println!("alpha_deg,Re,L,Cx,Cz,Cm,St,Cz_Cx");
    for (alpha, r) in res {
        match r {
            Ok((r, _, ok)) => println!(
                "{alpha},{},{},{:.5},{:.5},{:.5},{},{}",
                base.re,
                base.l,
                r.cx,
                r.cz,
                r.cm,
                if r.periodic {
                    format!("{:.4}", r.st)
                } else {
                    String::new()
                },
                if ok && r.cx.abs() > 1e-9 {
                    format!("{:.3}", r.cz / r.cx)
                } else {
                    "divergence".into()
                }
            ),
            Err(e) => println!("{alpha},,,,,,,{e}"),
        }
    }
    eprintln!("{:.1} s", t0.elapsed().as_secs_f64());
}

fn cmd_bench(a: &Args) {
    let l = a.get("l", 32usize);
    let steps = a.get("steps", 400usize);
    let base = TunnelConfig {
        l,
        re: 1000.0,
        shape: ShapeKind::Naca,
        params: ShapeParams {
            alpha_deg: 4.0,
            ..Default::default()
        },
        ..Default::default()
    };
    let cores = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    let mut cases = vec![
        ("CPU, 1 fil".to_string(), EngineKind::Cpu { threads: 1 }),
        (format!("CPU, {cores} fils"), EngineKind::Cpu { threads: cores }),
    ];
    if a.get("gpu", 1u8) == 1 {
        cases.push(("GPU".to_string(), EngineKind::Gpu));
    }
    println!(
        "NACA 0012, α = 4°, Re = 1000, L = {l} : grille {}×{} = {:.2} M mailles",
        24 * l,
        12 * l,
        (24 * l * 12 * l) as f64 / 1e6
    );
    for (name, engine) in cases {
        let mut wt = match WindTunnel::new(TunnelConfig { engine, ..base.clone() }) {
            Ok(w) => w,
            Err(e) => {
                println!("{name:12} : indisponible ({e})");
                continue;
            }
        };
        wt.run(20);
        let _ = wt.engine.take_forces();
        let t0 = Instant::now();
        wt.run(steps);
        let _ = wt.engine.take_forces(); // attend la fin du calcul GPU
        let secs = t0.elapsed().as_secs_f64();
        let mlups = (wt.nx * wt.ny) as f64 * steps as f64 / secs / 1e6;
        let tps = steps as f64 / secs / wt.steps_per_tstar();
        println!(
            "{name:12} : {mlups:8.1} MLUPS · {:7.1} pas/s · {tps:6.2} t*/s  [{}]",
            steps as f64 / secs,
            wt.engine.name()
        );
    }
}

fn main() {
    let a = Args::parse();
    match a.cmd.as_str() {
        "run" => cmd_run(&a),
        "polar" => cmd_polar(&a),
        "bench" => cmd_bench(&a),
        "" => {
            #[cfg(feature = "gui")]
            gui::main();
            #[cfg(not(feature = "gui"))]
            usage();
        }
        "help" | "-h" | "--help" => usage(),
        _ => usage(),
    }
}
