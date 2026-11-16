//! Diagnostic : efforts par bilan de quantité de mouvement sur une boîte,
//! contre l'échange de quantité de mouvement aux parois.
use lbm::bodies3d::{Body3D, Shape3D};
use lbm::gpu3d::Gpu3D;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let l: usize = a.get(1).and_then(|s| s.parse().ok()).unwrap_or(16);
    let depth: usize = a.get(2).and_then(|s| s.parse().ok()).unwrap_or(1);
    let re: f64 = a.get(3).and_then(|s| s.parse().ok()).unwrap_or(1e6);
    let alpha: f64 = a.get(4).and_then(|s| s.parse().ok()).unwrap_or(4.0);
    let u = 0.1f32;
    let spt = l as f64 / u as f64;
    let body = Body3D::new(Shape3D::Wing, l, alpha);
    let mut g = Gpu3D::with_refinement(&body, u as f64 * l as f64 / re, u, true, re >= 5e4, depth).unwrap();
    let ramp = lbm::cpu::Ramp { steps: (4.0 * spt) as u64, u };
    g.run((16.0 * spt) as usize, &ramp);
    let _ = g.take_forces();
    let q = 0.5 * (u as f64).powi(2) * body.area;
    let (nx, ny, nz) = (g.nx, g.ny, g.nz);
    let bb = body.bounds;
    for m in [0.4, 0.8] {
        let mut acc = [0.0f64; 3];
        let mut fex = [0.0f64; 3];
        let snaps = 8;
        for _ in 0..snaps {
            g.run((0.5 * spt) as usize, &ramp);
            let f = g.take_forces();
            fex[0] += f.fx / q / snaps as f64;
            fex[2] += f.fz / q / snaps as f64;
            let mac = g.macros();
            let lo = [bb.x0 - m * l as f64, bb.y0 - m * l as f64, bb.z0 - m * l as f64].map(|v| v.round() as usize);
            let hi = [bb.x1 + m * l as f64, bb.y1 + m * l as f64, bb.z1 + m * l as f64].map(|v| v.round() as usize);
            let hi = [hi[0].min(nx - 2), hi[1].min(ny - 2), hi[2].min(nz - 2)];
            let at = |x: usize, y: usize, z: usize| mac[x + nx * (y + ny * z)];
            // Flux sortant ∮ (ρ u u·n + p n) dA, p = ρ/3 ; force sur le corps = −flux.
            let mut flux = [0.0f64; 3];
            for (ax, sgn, pos) in [(0, -1.0, lo[0]), (0, 1.0, hi[0]), (1, -1.0, lo[1]), (1, 1.0, hi[1]), (2, -1.0, lo[2]), (2, 1.0, hi[2])] {
                let (a1, a2) = ((ax + 1) % 3, (ax + 2) % 3);
                for i in lo[a1]..=hi[a1] {
                    for j in lo[a2]..=hi[a2] {
                        let mut p3 = [0usize; 3];
                        p3[ax] = pos;
                        p3[a1] = i;
                        p3[a2] = j;
                        let c = at(p3[0], p3[1], p3[2]);
                        let (uu, rho) = ([c[0] as f64, c[1] as f64, c[2] as f64], c[3] as f64);
                        let un = sgn * uu[ax];
                        for k in 0..3 {
                            flux[k] += rho * uu[k] * un + if k == ax { sgn * rho / 3.0 } else { 0.0 };
                        }
                    }
                }
            }
            for k in 0..3 {
                acc[k] += -flux[k] / q / snaps as f64;
            }
        }
        println!(
            "L {l} raff {depth} Re {re:.0e} α {alpha} boîte ±{m} L : bilan Cx {:.4} Cz {:.4} · échange Cx {:.4} Cz {:.4}",
            acc[0], acc[2], fex[0], fex[2]
        );
    }
}
