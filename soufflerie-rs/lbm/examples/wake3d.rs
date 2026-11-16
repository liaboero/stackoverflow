//! Traînée en champ lointain (plans transversaux) contre l'intégration sur la surface.
//! wake3d <stl> <unite> <lref> <sref> <L> <raff> <re> <alpha> <t> <avg>
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
    let _ = sim.take_forces();
    let b = body.bounds;
    let lf = l as f64;
    let x_in = (b.x0 - 1.5 * lf).round() as usize;
    let fr = [0.1, 0.25, 0.5, 1.0, 2.0];
    let x_out: Vec<usize> = fr.iter().map(|f| (b.x1 + f * lf).round() as usize).collect();
    let n = avg as usize;
    let q = 0.5 * (u as f64).powi(2) * body.area;
    let mut acc = vec![[0.0f64; 7]; fr.len()];
    let mut surf = [0.0f64; 4];
    let mut bx = [[0.0f64; 2]; 2];
    for _ in 0..n {
        sim.run(spt as usize, &ramp);
        let s = sim.surface_forces().unwrap();
        surf[0] += s.pressure_wall[0] / q / n as f64;
        surf[1] += s.friction[0] / q / n as f64;
        surf[2] += (s.pressure_wall[2] + s.friction[2]) / q / n as f64;
        let wd = sim.wake_drag(x_in, &x_out);
        for mgn in [0.1, 0.2] {
            if let Some(f) = sim.box_forces(mgn) {
                bx[if mgn < 0.15 { 0 } else { 1 }] = [bx[if mgn < 0.15 { 0 } else { 1 }][0] + f[0] / q / n as f64, bx[if mgn < 0.15 { 0 } else { 1 }][1] + f[2] / q / n as f64];
                println!("  boîte {mgn} : Cx {:.4} Cz {:.4}", f[0] / q, f[2] / q);
            }
        }
        println!("  t {:.1} surface Cx {:.4} Cz {:.4} | plans {}", sim.time as f64 / spt, (s.pressure_wall[0] + s.friction[0]) / q, (s.pressure_wall[2] + s.friction[2]) / q, wd.iter().map(|w| format!("profil {:.4} hors {:+.4} induite {:.4}", w.loss_wake / q, w.loss_outside / q, w.induced / q)).collect::<Vec<_>>().join(" "));
        for (k, w) in wd.iter().enumerate() {
            let v = [w.loss_wake + w.loss_outside + w.induced, w.loss_wake, w.loss_outside, w.induced, 0.0, w.mass_ratio, w.wake_area];
            for i in 0..7 {
                acc[k][i] += v[i] / n as f64 / if i < 5 { q } else { 1.0 };
            }
        }
    }
    let f = sim.take_forces();
    let cz = surf[2];
    let span = b.y1 - b.y0;
    let ar = span * span / body.area;
    println!("surface : Cx pression {:.4} + frottement {:.4} = {:.4}, Cz {:.4} ; échange Cx {:.4}", surf[0], surf[1], surf[0] + surf[1], cz, f.fx / q);
    println!("boîte 0,1 L : Cx {:.4} Cz {:.4} ; boîte 0,2 L : Cx {:.4} Cz {:.4}", bx[0][0], bx[0][1], bx[1][0], bx[1][1]);
    println!("allongement {:.2}, CL²/(π A) = {:.4}", ar, cz * cz / (std::f64::consts::PI * ar));
    println!("plan (cordes après BF)  total  sillage  hors-sillage  induite  reste  débit  aire sillage");
    for (k, v) in acc.iter().enumerate() {
        println!("{:5.2}  {:.4}  {:.4}  {:+.4}  {:.4}  {:+.4}  {:.5}  {:.1}", fr[k], v[0], v[1], v[2], v[3], v[4], v[5], v[6]);
    }
}
