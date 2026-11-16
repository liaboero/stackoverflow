//! Hauteurs de relevé de la pression pariétale contre le bilan sur une boîte.
//! (relit SOUFFLERIE_PWH="h1,h2" dans Gpu3D::surface_forces)
//! pwh3d <stl> <unite> <lref> <sref> <L> <raff> <re> <alpha> <t> <avg>
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
    let (l, raff, re, t, avg) = (g(5, 16.0) as usize, g(6, 2.0) as usize, g(7, 1e6), g(9, 20.0), g(10, 8.0));
    let (body, _) = Body3D::from_mesh(&tris, l, &setup).unwrap();
    let u = 0.1f32;
    let mut sim = Gpu3D::with_refinement(&body, u as f64 * l as f64 / re, u, true, true, raff).unwrap();
    let spt = l as f64 / u as f64;
    let ramp = lbm::cpu::Ramp { steps: (4.0 * spt) as u64, u };
    sim.run(((t - avg) * spt) as usize, &ramp);
    let q = 0.5 * (u as f64).powi(2) * body.area;
    let hs = ["0.5,1", "0.75,1.5", "1,2", "1.5,2.5", "9,9"];
    let n = (2.0 * avg) as usize;
    let mut acc = vec![[0.0f64; 2]; hs.len() + 1];
    for _ in 0..n {
        sim.run((0.5 * spt) as usize, &ramp);
        for (k, h) in hs.iter().enumerate() {
            std::env::set_var("SOUFFLERIE_PWH", h);
            let s = sim.surface_forces().unwrap();
            acc[k][0] += s.total()[0] / q / n as f64;
            acc[k][1] += s.total()[2] / q / n as f64;
        }
        let b = sim.box_forces(0.15).unwrap();
        acc[hs.len()][0] += b[0] / q / n as f64;
        acc[hs.len()][1] += b[2] / q / n as f64;
    }
    for (k, h) in hs.iter().enumerate() {
        println!("relevés {h} (9,9 : ancien) : Cx {:.4} Cz {:.4}", acc[k][0], acc[k][1]);
    }
    println!("boîte 0,15 : Cx {:.4} Cz {:.4}", acc[hs.len()][0], acc[hs.len()][1]);
}
