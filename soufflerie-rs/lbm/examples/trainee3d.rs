//! Diagnostic : traînée 3D d'une maquette STL — échange aux parois, gradient
//! de pression le long de la veine loin du corps, et poussée correspondante
//! sur le volume du corps (−V dp/dx).
use lbm::bodies3d::Body3D;
use lbm::gpu3d::Gpu3D;
use lbm::stl::{MeshSetup, Unit};
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let path = a.get(1).expect("fichier STL");
    let marge: f64 = a.get(2).and_then(|s| s.parse().ok()).unwrap_or(0.4);
    let alpha: f64 = a.get(3).and_then(|s| s.parse().ok()).unwrap_or(4.0);
    let l: usize = a.get(4).and_then(|s| s.parse().ok()).unwrap_or(16);
    let depth: usize = a.get(5).and_then(|s| s.parse().ok()).unwrap_or(1);
    let re = 1e6;
    let tris = lbm::stl::read_stl(std::path::Path::new(path)).unwrap();
    let setup = MeshSetup {
        unit: Unit::M,
        alpha_deg: alpha,
        ref_length: Some(1.0),
        ref_area: Some(4.0),
        margin: marge,
        ..Default::default()
    };
    let (body, _) = Body3D::from_mesh(&tris, l, &setup).unwrap();
    let u = 0.1f32;
    let spt = l as f64 / u as f64;
    let mut g = Gpu3D::with_refinement(&body, u as f64 * l as f64 / re, u, true, true, depth).unwrap();
    let ramp = lbm::cpu::Ramp { steps: (4.0 * spt) as u64, u };
    g.run((16.0 * spt) as usize, &ramp);
    let _ = g.take_forces();
    g.run((4.0 * spt) as usize, &ramp);
    let f = g.take_forces();
    let q = 0.5 * (u as f64).powi(2) * body.area;
    let mac = g.macros();
    let (nx, ny, nz) = (g.nx, g.ny, g.nz);
    // ρ moyen par plan x, loin du corps et hors des bandes absorbantes.
    let (y0, y1, z0, z1) = (ny / 8 + 2, ny - ny / 8 - 2, nz / 8 + 2, nz - nz / 8 - 2);
    let bb = body.bounds;
    let lf = l as f64;
    let prof: Vec<f64> = (0..nx)
        .map(|x| {
            let (mut s, mut n) = (0.0, 0);
            for z in z0..z1 {
                for y in y0..y1 {
                    let (yf, zf) = (y as f64, z as f64);
                    let far = yf < bb.y0 - 0.5 * lf || yf > bb.y1 + 0.5 * lf || zf < bb.z0 - 0.5 * lf || zf > bb.z1 + 0.5 * lf;
                    if far {
                        s += mac[x + nx * (y + ny * z)][0] as f64 / u as f64;
                        n += 1;
                    }
                }
            }
            s / n.max(1) as f64
        })
        .collect();
    let xa = (bb.x0.round() as usize).saturating_sub(l / 2);
    let xb = (bb.x1.round() as usize + l / 2).min(nx - 1);
    let dpdx = 0.0 * (prof[xb] - prof[xa]) / (xb - xa) as f64;
    let buoy = -body.mesh_volume * dpdx / q;
    let s: Vec<String> = (0..nx).step_by((nx / 12).max(1)).map(|x| format!("{x}:{:.4}", prof[x])).collect();
    println!(
        "marge {marge} grille {nx}×{ny}×{nz} : Cx échange {:.4} Cz {:.4} · dp/dx {:.3e} · poussée −V dp/dx = {:.4} en Cx · u/U loin du corps le long de x : {}",
        f.fx / q,
        f.fz / q,
        dpdx,
        buoy,
        s.join(" ")
    );
}
