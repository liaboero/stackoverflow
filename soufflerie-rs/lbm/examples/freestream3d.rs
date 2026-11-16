//! Écoulement amont réellement vu par la maquette : ρ et u le long de la
//! veine (ligne à mi-envergure, à mi-chemin entre la maquette et le bord
//! haut), et moyennes sur des sections transversales.
//! freestream3d <stl> <unite> <lref> <sref> <L> <raff> <re> <alpha> <t>
use lbm::bodies3d::Body3D;
use lbm::gpu3d::Gpu3D;
use lbm::stl::{MeshSetup, Unit};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let g = |i: usize, d: f64| a.get(i).and_then(|v| v.parse().ok()).unwrap_or(d);
    let tris = lbm::stl::read_stl(std::path::Path::new(&a[1])).unwrap();
    let setup = MeshSetup {
        unit: Unit::parse(&a[2]).unwrap(),
        ref_length: Some(g(3, 1.0)),
        ref_area: Some(g(4, 4.0)),
        margin: 0.4,
        alpha_deg: g(8, 4.0),
        ..Default::default()
    };
    let (l, raff, re, t) = (g(5, 16.0) as usize, g(6, 1.0) as usize, g(7, 1e6), g(9, 16.0));
    let (body, _) = Body3D::from_mesh(&tris, l, &setup).unwrap();
    let u = 0.1f32;
    let mut sim = Gpu3D::with_refinement(&body, u as f64 * l as f64 / re, u, true, true, raff).unwrap();
    let spt = l as f64 / u as f64;
    let ramp = lbm::cpu::Ramp { steps: (4.0 * spt) as u64, u };
    sim.run((t * spt) as usize, &ramp);
    let n = 8;
    let (nx, ny, nz) = (sim.nx, sim.ny, sim.nz);
    let b = body.bounds;
    let yc = (0.5 * (b.y0 + b.y1)) as usize;
    let zt = ((b.z1 + nz as f64 * 0.75) / 2.0).min(nz as f64 - 2.0) as usize;
    let xs: Vec<usize> = [-3.0, -2.0, -1.0, -0.5, 0.5, 1.5].iter().map(|f: &f64| ((if *f < 0.0 { b.x0 } else { b.x1 }) + f * l as f64).clamp(1.0, nx as f64 - 2.0) as usize).collect();
    let mut line = vec![[0.0f64; 2]; xs.len()];
    let mut sec = vec![[0.0f64; 3]; xs.len()];
    for _ in 0..n {
        sim.run((0.5 * spt) as usize, &ramp);
        let m = sim.macros();
        for (k, &x) in xs.iter().enumerate() {
            let c = m[x + nx * (yc + ny * zt)];
            line[k][0] += c[0] as f64 / n as f64;
            line[k][1] += c[3] as f64 / n as f64;
            // Section entière : ρ moyen, débit ρu moyen, u moyen.
            let (mut r, mut q, mut uu) = (0.0, 0.0, 0.0);
            for z in 0..nz {
                for y in 0..ny {
                    let c = m[x + nx * (y + ny * z)];
                    r += c[3] as f64;
                    q += (c[3] * c[0]) as f64;
                    uu += c[0] as f64;
                }
            }
            let s = (ny * nz) as f64 * n as f64;
            sec[k][0] += r / s;
            sec[k][1] += q / s;
            sec[k][2] += uu / s;
        }
    }
    println!("u_in = {u} ; ligne au-dessus de la maquette (z = {zt}) et sections :");
    for (k, &x) in xs.iter().enumerate() {
        println!(
            "x = {x:4} : ligne u/u_in {:.4} ρ {:.4} · section ρ {:.4} ρu/u_in {:.4} u/u_in {:.4}",
            line[k][0] / u as f64, line[k][1], sec[k][0], sec[k][1] / u as f64, sec[k][2] / u as f64
        );
    }
}
