//! Soufflerie de poche, version native.
//!
//! Sans argument : interface graphique. Sous-commandes :
//!   soufflerie run    [options]          un calcul, moyennes et Strouhal
//!   soufflerie polar  --alpha a:b:pas    polaire (incidences calculées en parallèle), CSV
//!   soufflerie panel  --alpha a:b:pas    polaire instantanée : panneaux + couche limite, CSV
//!   soufflerie bench  [--l 32]           vitesse des moteurs (millions de mailles par seconde)
//!
//! Options : --shape cylinder|square|ellipse|plate|naca|triangle  --naca 2412
//!           --alpha 4  --re 1000  --l 32  --t 40  --avg 15  --engine cpu|gpu
//!           --threads 0  --spin 1  --thickness 0.4  --collision rr|bgk
//!           --les auto|on|off  --u 0.1  --turb auto|sa|off  --transition free|turb

use std::collections::HashMap;
use std::time::Instant;

use lbm::cpu::{Collision, LesMode};
use lbm::panel::{compressible_correction, panel_solve, section_coefficients, PanelOptions, PanelResult, PanelShape};
use lbm::stats::{analyze, Analysis};
use lbm::{EngineKind, Naca, ShapeKind, ShapeParams, Transition, TunnelConfig, Turbulence, WindTunnel};
use rayon::prelude::*;

#[cfg(feature = "gui")]
mod gui;
#[cfg(feature = "gui")]
mod launch;
#[cfg(feature = "gui")]
mod worker;

fn usage() -> ! {
    eprintln!(
        "Soufflerie de poche (native)\n\n\
         soufflerie [options]             interface graphique (--field eddy, --graph cp|cf)\n\
         soufflerie run   [options]       un calcul\n\
         soufflerie polar --alpha -4:12:2 [options]   polaire, CSV\n\
         soufflerie panel --alpha -4:16:1 [options]   polaire panneaux + couche limite, CSV\n\
         soufflerie bench [--l 32] [--steps 400]       vitesse des moteurs\n\n\
         Options : --shape naca --naca 2412 --alpha 4 --re 1000 --l 32 --t 40 --avg 15\n\
         \x20         --engine cpu|gpu --threads 0 --spin 0 --thickness 0.4\n\
         \x20         --collision rr|bgk --les auto|on|off --u 0.1\n\
         \x20         --turb auto|sa|off (Spalart–Allmaras dès Re = 5·10⁴ en auto)\n\
         \x20         --transition free|turb (libre : prédite par la couche limite intégrale)\n\
         \x20         --mach 0.3 (panel : correction de compressibilité)\n\
         \x20         --rendu auto|dx12|vulkan|glow (interface : rendu graphique)"
    );
    std::process::exit(2)
}

pub struct Args {
    cmd: String,
    pub opts: HashMap<String, String>,
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

    pub fn config(&self) -> TunnelConfig {
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
            turbulence: match self.get("turb", "auto".to_string()).as_str() {
                "sa" => Turbulence::Sa,
                "off" | "none" => Turbulence::Off,
                _ => Turbulence::Auto,
            },
            transition: if self.get("transition", "free".to_string()) == "turb" {
                Transition::Turbulent
            } else {
                Transition::Free
            },
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
    if wt.sa_active() {
        let profiled = wt.panel.is_some();
        println!(
            "Spalart–Allmaras + loi de paroi · y⁺ ≈ {:.0} au premier nœud · transition {}",
            wt.first_cell_yplus(),
            match (cfg.transition, profiled) {
                (Transition::Turbulent, _) => "au bord d'attaque".to_string(),
                (Transition::Free, true) => {
                    let r = wt.panel.as_ref().unwrap();
                    format!(
                        "libre (extrados {}, intrados {})",
                        pct(r.transition_upper),
                        pct(r.transition_lower)
                    )
                }
                (Transition::Free, false) => "non imposée (couche limite turbulente)".to_string(),
            }
        );
    }
    if let Some(r) = &wt.panel {
        let k = section_coefficients(r);
        println!(
            "Panneaux + couche limite : Cx = {:.4}  Cz = {:.4}  Cm = {:.4}{}",
            k.cd,
            k.cl,
            k.cm,
            if k.kirchhoff {
                format!(
                    "  (décollement à {} de corde : modèle de Kirchhoff)",
                    pct(r.separation_upper.min(r.separation_lower))
                )
            } else {
                String::new()
            }
        );
    }
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

/// Plage d'incidences a:b:pas (ou une seule valeur).
fn alpha_range(a: &Args, default: &str) -> Vec<f64> {
    let spec = a.get("alpha", default.to_string());
    let p: Vec<f64> = spec.split(':').filter_map(|s| s.parse().ok()).collect();
    match p.len() {
        1 => return p,
        3 if p[2] > 0.0 => {}
        _ => usage(),
    }
    let mut alphas = Vec::new();
    let mut k = 0;
    loop {
        let x = p[0] + k as f64 * p[2];
        if x > p[1] + 1e-9 {
            break;
        }
        alphas.push((x * 1e6).round() / 1e6);
        k += 1;
    }
    alphas
}

/// Méthode des panneaux + couche limite intégrale pour la maquette de `cfg`.
fn panel_of(cfg: &TunnelConfig, alpha: f64) -> Option<PanelResult> {
    let shape = PanelShape::from_shape(cfg.shape, &cfg.params)?;
    let opts = match cfg.transition {
        Transition::Free => PanelOptions::free(cfg.re),
        Transition::Turbulent => PanelOptions::tripped(cfg.re, 0.0),
    };
    let r = panel_solve(&shape, alpha, &opts);
    (r.cl.is_finite() && r.cd.is_finite()).then_some(r)
}

fn pct(v: f64) -> String {
    if v.is_finite() {
        format!("{:.0} %", v.clamp(0.0, 1.0) * 100.0)
    } else {
        "—".into()
    }
}

fn cmd_panel(a: &Args) {
    let mut cfg = a.config();
    if !a.opts.contains_key("shape") {
        cfg.shape = ShapeKind::Naca;
    }
    if PanelShape::from_shape(cfg.shape, &cfg.params).is_none() {
        eprintln!("La méthode des panneaux traite les profils NACA, les ellipses et les plaques.");
        std::process::exit(2);
    }
    if !a.opts.contains_key("re") {
        cfg.re = 1e6;
    }
    let mach = a.get("mach", 0.0f64);
    let alphas = alpha_range(a, "-4:16:1");
    eprintln!(
        "{} · Re = {} · transition {}{}",
        match cfg.shape {
            ShapeKind::Naca => format!("NACA {}", cfg.params.naca.code()),
            k => k.label().to_string(),
        },
        cfg.re,
        if cfg.transition == Transition::Free {
            "libre (eᴺ, N = 9)"
        } else {
            "au bord d'attaque"
        },
        if mach > 0.0 {
            format!(" · Mach {mach}")
        } else {
            String::new()
        }
    );
    println!("alpha_deg,Cz,Cx,Cm,Cz_attache,x_tr_extrados,x_tr_intrados,x_decollement,kirchhoff,Cp_min,critique");
    let rows: Vec<_> = alphas.par_iter().map(|&al| (al, panel_of(&cfg, al))).collect();
    for (al, r) in rows {
        let Some(r) = r else {
            println!("{al},,,,,,,,,,");
            continue;
        };
        let k = section_coefficients(&r);
        let cp_min = r.cp.iter().cloned().fold(f64::INFINITY, f64::min);
        let (mut cl, mut cm, mut crit) = (k.cl, k.cm, false);
        let mut cp_min_c = cp_min;
        if mach > 0.0 {
            let c = compressible_correction(mach, cp_min);
            cl *= c.factor;
            cm *= c.factor;
            cp_min_c = c.cp_min;
            crit = c.critical;
        }
        let f = |v: f64| {
            if v.is_finite() {
                format!("{v:.4}")
            } else {
                String::new()
            }
        };
        println!(
            "{al},{cl:.4},{:.5},{cm:.4},{:.4},{},{},{},{},{cp_min_c:.3},{}",
            k.cd,
            r.cl,
            f(r.transition_upper),
            f(r.transition_lower),
            f(r.separation_upper.min(r.separation_lower)),
            k.kirchhoff as u8,
            crit as u8
        );
    }
}

fn cmd_polar(a: &Args) {
    let alphas = alpha_range(a, "-4:12:2");
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
    println!("alpha_deg,Re,L,Cx,Cz,Cm,St,Cz_Cx,Cx_panneaux,Cz_panneaux");
    for (alpha, r) in res {
        let pan = match panel_of(&base, alpha).filter(|_| base.re >= lbm::tunnel::RE_PANEL_MIN) {
            Some(p) => {
                let k = section_coefficients(&p);
                format!(",{:.5},{:.4}", k.cd, k.cl)
            }
            None => ",,".into(),
        };
        match r {
            Ok((r, _, ok)) => println!(
                "{alpha},{},{},{:.5},{:.5},{:.5},{},{}{pan}",
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
            Err(e) => println!("{alpha},,,,,,,{e}{pan}"),
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
        "panel" => cmd_panel(&a),
        "bench" => cmd_bench(&a),
        "" => {
            #[cfg(feature = "gui")]
            launch::gui(&a);
            #[cfg(not(feature = "gui"))]
            usage();
        }
        "help" | "-h" | "--help" => usage(),
        _ => usage(),
    }
}
