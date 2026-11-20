//! Canal plan (Poiseuille) entre deux parois, droites ou inclinées sur la
//! grille, paroi volumique par surfels contre rebond de Bouzidi.
//!
//! canal_surfels <pente 0 | 0.5> <c₁> <largeur selon z> <ν> <pas>
//!
//! Fluide : c₁ < z − pente·y < c₁ + largeur (mod nz = 20) ; écoulement selon
//! x entretenu par une force volumique. Débit Σ V u_x contre a H³ L / (12 ν),
//! et écart L2 du profil aux nœuds des voxels entièrement fluides.
use std::sync::Arc;

use lbm::bodies3d::{Body3D, Bounds, Shape3D};
use lbm::cpu::Ramp;
use lbm::gpu3d::{Gpu3D, Turbulence};
use lbm::stl::{MeshBody, Tri};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let g = |i: usize, d: f64| a.get(i).and_then(|v| v.parse().ok()).unwrap_or(d);
    let slope = g(1, 0.5);
    let c1 = g(2, 3.3);
    let hz = g(3, 12.0);
    let nu = g(4, 0.1);
    let steps = g(5, 8000.0) as usize;
    let nz = 20usize;
    let ny = if slope == 0.0 { 4 } else { (nz as f64 / slope).round() as usize };
    let nx = 4usize;
    let c2 = c1 + hz;
    let cs = 1.0 / (1.0 + slope * slope).sqrt();
    let h = hz * cs;
    let mut tris: Vec<Tri> = Vec::new();
    let (x0, x1, y0, y1) = (-1.0, nx as f64 + 1.0, -3.0, ny as f64 + 3.0);
    let plane = |c: f64, up: bool, t: &mut Vec<Tri>| {
        let p = |x: f64, y: f64| [x, y, slope * y + c];
        let (a, b, cc, d) = (p(x0, y0), p(x1, y0), p(x1, y1), p(x0, y1));
        if up {
            t.push([a, b, cc]);
            t.push([a, cc, d]);
        } else {
            t.push([a, cc, b]);
            t.push([a, d, cc]);
        }
    };
    let kmin = -((slope * (ny as f64 + 3.0) / nz as f64).ceil() as i32) - 2;
    for k in kmin..=1 {
        let k = k as f64 * nz as f64;
        plane(c2 + k, false, &mut tris);
        plane(c1 + nz as f64 + k, true, &mut tris);
    }
    let acc = 1.2e-4 * (nu / 0.1);
    let lw = if slope == 0.0 { ny as f64 } else { (ny as f64).hypot(nz as f64) };
    let qa = acc * h.powi(3) / (12.0 * nu) * lw;
    let umax = acc * h * h / (8.0 * nu);
    for surf in [true, false] {
        let mut b = Body3D::new(Shape3D::Sphere, 8, 0.0);
        let mesh = MeshBody::new(tris.clone(), true);
        b.shape = Shape3D::Mesh;
        b.nx = nx;
        b.ny = ny;
        b.nz = nz;
        b.bounds = Bounds { x0: mesh.lo[0], x1: mesh.hi[0], y0: mesh.lo[1], y1: mesh.hi[1], z0: mesh.lo[2], z1: mesh.hi[2] };
        b.mesh = Some(Arc::new(mesh));
        b.features = None;
        b.surfels = surf;
        let mut sim = Gpu3D::with_turbulence(&b, nu, 0.0, false, Turbulence::Laminar).unwrap();
        sim.set_periodic_x();
        sim.set_accel(acc);
        sim.run(steps, &Ramp { steps: 0, u: 0.0 });
        let m = sim.macros();
        let vols = sim.cell_volumes();
        let (mut q, mut e2, mut n2) = (0.0, 0.0, 0.0);
        let mut prof = Vec::new();
        for z in 0..nz {
            for y in 0..ny {
                let c = nx * (y + ny * z);
                q += vols[c] as f64 * m[c][0] as f64;
                let eta = ((z as f64 - slope * y as f64 - c1).rem_euclid(nz as f64)) * cs;
                if vols[c] == 1.0 && eta < h {
                    let ua = acc * eta * (h - eta) / (2.0 * nu);
                    e2 += (m[c][0] as f64 - ua).powi(2);
                    n2 += umax * umax;
                }
                if y == 0 {
                    prof.push(format!("{z}:{:.3}/{:.3}", m[c][0] as f64 / umax, if eta < h { acc * eta * (h - eta) / (2.0 * nu) / umax } else { 0.0 }));
                }
            }
        }
        println!(
            "{} : débit {:+.2} %, profil L2/u_max {:.2} %, mailles pariétales {}, {:?}",
            if surf { "surfels" } else { "Bouzidi" },
            100.0 * (q / qa - 1.0),
            100.0 * (e2 / n2).sqrt(),
            sim.wall_count,
            sim.surfel_info().map(|s| s.0)
        );
        if std::env::var("CANAL_PROFIL").is_ok() {
            println!("  {}", prof.join(" "));
        }
    }
}
