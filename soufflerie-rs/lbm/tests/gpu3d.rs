//! Soufflerie 3D (ignorée sans adaptateur GPU) : sphère à Re = 100, mêmes
//! valeurs que l'application web (Cx = 1,26 à 8 mailles par diamètre).
#![cfg(feature = "gpu")]

use lbm::bodies3d::{Body3D, Shape3D};
use lbm::cpu::Ramp;
use lbm::gpu3d::{Gpu3D, SliceField};

#[test]
fn sphere_re100() {
    let (l, u) = (8usize, 0.1f32);
    let body = Body3D::new(Shape3D::Sphere, l, 0.0);
    let mut g = match Gpu3D::new(&body, u as f64 * l as f64 / 100.0, u, true) {
        Ok(g) => g,
        Err(e) => {
            eprintln!("GPU ignoré : {e}");
            return;
        }
    };
    println!("{} · {} mailles pariétales", g.adapter_name(), g.wall_count);
    let spt = l as f64 / u as f64;
    let ramp = Ramp {
        steps: (4.0 * spt) as u64,
        u,
    };
    let q = 0.5 * (u as f64).powi(2) * body.area;
    // 30 t*, moyenne sur les 10 derniers.
    g.run((20.0 * spt) as usize, &ramp);
    let _ = g.take_forces();
    g.run((10.0 * spt) as usize, &ramp);
    let f = g.take_forces();
    let (cx, cz, cy) = (f.fx / q, f.fz / q, f.fy / q);
    println!("Cx = {cx:.4} (web 1,26), Cz = {cz:.4}, Cy = {cy:.4}");
    assert!(!g.lost && (cx - 1.26).abs() < 0.03 && cz.abs() < 0.03 && cy.abs() < 0.03);
    // Coupe : corps au centre, sillage à vitesse réduite.
    let s = g.slice(1, body.center[1].round() as usize, SliceField::Speed, l as f64);
    let c = body.center;
    let at = |x: f64, z: f64| s.data[x.round() as usize + s.w * z.round() as usize];
    assert!(at(c[0], c[2]) > 1e29, "le centre de la sphère doit être solide");
    let wake = at(c[0] + 1.0 * l as f64, c[2]);
    println!("|u|/U à 1 D en aval : {wake:.3}");
    assert!(wake < 0.5);
}
