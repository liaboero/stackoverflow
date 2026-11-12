//! Moteur GPU contre moteur processeur (ignoré sans adaptateur GPU).
#![cfg(feature = "gpu")]

use lbm::{EngineKind, Naca, ShapeKind, ShapeParams, TunnelConfig, WindTunnel};

fn coeffs(cfg: TunnelConfig, steps: usize, every: usize) -> Option<Vec<(f64, f64, f64)>> {
    let mut wt = match WindTunnel::new(cfg) {
        Ok(w) => w,
        Err(e) => {
            eprintln!("GPU ignoré : {e}");
            return None;
        }
    };
    println!("moteur : {}", wt.engine.name());
    let mut out = Vec::new();
    let mut done = 0;
    while done < steps {
        wt.run(every);
        done += every;
        let s = wt.engine.take_forces();
        let c = wt.coefficients(s);
        out.push((c.cx, c.cz, c.cm));
    }
    Some(out)
}

#[test]
fn gpu_egal_cpu() {
    for (shape, re, alpha, spin) in [
        (ShapeKind::Cylinder, 100.0, 0.0, 0.0),
        (ShapeKind::Naca, 2000.0, 6.0, 0.0),
        (ShapeKind::Cylinder, 100.0, 0.0, 1.0),
    ] {
        let base = TunnelConfig {
            l: 10,
            re,
            shape,
            params: ShapeParams {
                alpha_deg: alpha,
                naca: Naca::parse("2412").unwrap(),
                spin_ratio: spin,
                ..Default::default()
            },
            ..Default::default()
        };
        let cpu = coeffs(
            TunnelConfig {
                engine: EngineKind::Cpu { threads: 0 },
                ..base.clone()
            },
            2000,
            100,
        )
        .unwrap();
        let Some(gpu) = coeffs(
            TunnelConfig {
                engine: EngineKind::Gpu,
                ..base
            },
            2000,
            100,
        ) else {
            return;
        };
        let mut dmax: f64 = 0.0;
        for (a, b) in cpu.iter().zip(&gpu) {
            dmax = dmax
                .max((a.0 - b.0).abs())
                .max((a.1 - b.1).abs())
                .max((a.2 - b.2).abs());
        }
        let (a, b) = (cpu.last().unwrap(), gpu.last().unwrap());
        println!("{shape:?} Re = {re} : CPU Cx, Cz, Cm = {:.4}, {:.4}, {:.4} ; GPU = {:.4}, {:.4}, {:.4} ; écart max {dmax:.2e}", a.0, a.1, a.2, b.0, b.1, b.2);
        assert!(dmax < 2e-3, "écart GPU/CPU {dmax}");
    }
}
