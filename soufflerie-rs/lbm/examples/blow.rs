//! Diagnostic : où naît la divergence 3D (maquette STL).
use lbm::bodies3d::Body3D;
use lbm::gpu3d::Gpu3D;
use lbm::stl::{MeshSetup, Unit};
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let g = |i: usize, d: f64| a.get(i).and_then(|v| v.parse().ok()).unwrap_or(d);
    let tris = lbm::stl::read_stl(std::path::Path::new(&a[1])).unwrap();
    let setup = MeshSetup { unit: Unit::parse(&a[2]).unwrap(), ref_length: Some(g(3, 1.0)), ref_area: Some(g(4, 4.0)), margin: 0.4, alpha_deg: g(8, 4.0), ..Default::default() };
    let (l, raff, re) = (g(5, 16.0) as usize, g(6, 2.0) as usize, g(7, 1e6));
    let (body, _) = Body3D::from_mesh(&tris, l, &setup).unwrap();
    let b = body.bounds;
    println!("corps (veine) x {:.1}..{:.1} y {:.1}..{:.1} z {:.1}..{:.1}", b.x0, b.x1, b.y0, b.y1, b.z0, b.z1);
    let u = 0.1f32;
    let mut sim = Gpu3D::with_refinement(&body, u as f64 * l as f64 / re, u, true, true, raff).unwrap();
    let spt = l as f64 / u as f64;
    let ramp = lbm::cpu::Ramp { steps: (4.0 * spt) as u64, u };
    let (t0, dt) = (g(9, 0.0), g(10, 0.05));
    sim.run((t0 * spt) as usize, &ramp);
    let mut t = t0;
    loop {
        sim.run((dt * spt) as usize, &ramp);
        t += dt;
        let (m, d, org, s) = sim.finest_macros();
        let mut v: Vec<(f32, usize)> = m.iter().enumerate().map(|(i, c)| ((c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt(), i)).collect();
        v.sort_by(|a, b| b.0.total_cmp(&a.0));
        let (_, i0) = v[0];
        println!("t* {t:.2} max |u|/U {:.2} en x {} y {} z {}", v[0].0 / u, i0 % d[0], (i0 / d[0]) % d[1], i0 / (d[0] * d[1]));
        if !(v[0].0 < 0.3) || t > 30.0 {
            println!("t* {t:.2} max |u|/U {:.2}", v[0].0 / u);
            for &(sp, i) in v.iter().take(12) {
                let p = [i % d[0], (i / d[0]) % d[1], i / (d[0] * d[1])];
                let w = [0, 1, 2].map(|k| org[k] + p[k] as f64 / s);
                println!("  |u|/U {:.2} veine ({:.2}, {:.2}, {:.2}) x/c {:.3} y/b {:.3} z-zb {:.2}", sp / u, w[0], w[1], w[2], (w[0] - b.x0) / (b.x1 - b.x0), (w[1] - b.y0) / (b.y1 - b.y0), w[2] - 0.5 * (b.z0 + b.z1));
            }
            break;
        }
    }
}
