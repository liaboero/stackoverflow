//! Moteur GPU contre moteur processeur (ignoré sans adaptateur GPU).
#![cfg(feature = "gpu")]

use lbm::{EngineKind, Naca, ShapeKind, ShapeParams, Transition, TunnelConfig, Turbulence, WindTunnel};

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
    // Spalart–Allmaras : transition libre (zone laminaire de la couche limite
    // intégrale), puis turbulent avec paroi mobile.
    for (shape, re, alpha, spin, sa) in [
        (ShapeKind::Cylinder, 100.0, 0.0, 0.0, false),
        (ShapeKind::Naca, 2000.0, 6.0, 0.0, false),
        (ShapeKind::Cylinder, 100.0, 0.0, 1.0, false),
        (ShapeKind::Naca, 1e6, 4.0, 0.0, true),
        (ShapeKind::Cylinder, 1e5, 0.0, 1.0, true),
    ] {
        let base = TunnelConfig {
            l: 10,
            re,
            shape,
            turbulence: if sa { Turbulence::Sa } else { Turbulence::Off },
            transition: if spin > 0.0 {
                Transition::Turbulent
            } else {
                Transition::Free
            },
            params: ShapeParams {
                alpha_deg: alpha,
                naca: Naca::parse(if sa { "0012" } else { "2412" }).unwrap(),
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
        println!("{shape:?} Re = {re} SA {sa} : CPU Cx, Cz, Cm = {:.4}, {:.4}, {:.4} ; GPU = {:.4}, {:.4}, {:.4} ; écart max {dmax:.2e}", a.0, a.1, a.2, b.0, b.1, b.2);
        assert!(dmax < 2e-3, "écart GPU/CPU {dmax}");
    }
}

#[test]
fn gpu_demi_precision() {
    // Populations en demi-précision : mêmes efforts à quelques 10⁻³ près.
    for (shape, re, alpha, sa) in [
        (ShapeKind::Cylinder, 100.0, 0.0, false),
        (ShapeKind::Naca, 1e6, 4.0, true),
    ] {
        let base = TunnelConfig {
            l: 10,
            re,
            shape,
            engine: EngineKind::Gpu,
            turbulence: if sa { Turbulence::Sa } else { Turbulence::Off },
            params: ShapeParams {
                alpha_deg: alpha,
                naca: Naca::parse("0012").unwrap(),
                ..Default::default()
            },
            ..Default::default()
        };
        let Some(single) = coeffs(base.clone(), 2000, 100) else {
            return;
        };
        let half = coeffs(
            TunnelConfig {
                gpu_half: true,
                ..base
            },
            2000,
            100,
        )
        .unwrap();
        let (a, b) = (single.last().unwrap(), half.last().unwrap());
        let d = (a.0 - b.0).abs().max((a.1 - b.1).abs());
        println!(
            "{shape:?} Re = {re} : simple Cx, Cz = {:.4}, {:.4} ; demi = {:.4}, {:.4} ; écart {d:.2e}",
            a.0, a.1, b.0, b.1
        );
        assert!(d < 5e-3 + 0.01 * a.0.abs(), "écart demi/simple {d}");
    }
}
