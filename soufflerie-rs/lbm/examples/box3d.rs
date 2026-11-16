//! Bilan sur des boîtes de tailles croissantes : extrapolation à volume nul.
//! box3d <stl> <unite> <lref> <sref> <L> <raff> <re> <alpha> <t> <avg>
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
    let ms = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35];
    let n = (2.0 * avg) as usize;
    let mut acc = vec![[0.0f64; 3]; ms.len()];
    let mut ok = vec![true; ms.len()];
    for _ in 0..n {
        sim.run((0.5 * spt) as usize, &ramp);
        for (k, r) in sim.box_forces_vol(&ms).into_iter().enumerate() {
            match r {
                Some((f, v)) => {
                    acc[k][0] += f[0] / q / n as f64;
                    acc[k][1] += f[2] / q / n as f64;
                    acc[k][2] = v / (body.area * l as f64);
                }
                None => ok[k] = false,
            }
        }
    }
    let pts: Vec<[f64; 3]> = acc.iter().zip(&ok).filter(|(_, &o)| o).map(|(a, _)| *a).collect();
    for (m, (a, o)) in ms.iter().zip(acc.iter().zip(&ok)) {
        if *o {
            println!("marge {m:.2} : V/(S L) {:.3}  Cx {:.4}  Cz {:.4}", a[2], a[0], a[1]);
        }
    }
    // Moindres carrés Cx = Cx0 + G V.
    let np = pts.len() as f64;
    let (sx, sy) = (pts.iter().map(|p| p[2]).sum::<f64>(), pts.iter().map(|p| p[0]).sum::<f64>());
    let sxx = pts.iter().map(|p| p[2] * p[2]).sum::<f64>();
    let sxy = pts.iter().map(|p| p[2] * p[0]).sum::<f64>();
    let gsl = (np * sxy - sx * sy) / (np * sxx - sx * sx);
    println!("extrapolation à volume nul : Cx0 = {:.4} (pente {:.5} par S L)", (sy - gsl * sx) / np, gsl);
}
