//! Soufflerie 3D (ignorée sans adaptateur GPU) : sphère à Re = 100, mêmes
//! valeurs que l'application web (Cx = 1,26 à 8 mailles par diamètre).
#![cfg(feature = "gpu")]

use lbm::bodies3d::{Body3D, Shape3D};
use lbm::cpu::Ramp;
use lbm::gpu3d::{Camera, Gpu3D, RenderOptions, SliceField};

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
    // Vue 3D vue de face (de l'amont) : sphère sombre au centre, fond clair au bord.
    let (w, h) = (96usize, 64usize);
    let img = g.render(&RenderOptions {
        camera: Camera::orbit(c, 0.0, 0.0, 4.0 * l as f64),
        width: w,
        height: h,
        plane_axis: 1,
        plane_pos: 0.0,
        vorticity: true,
        iso: false,
        iso_level: 0.4,
        l: l as f64,
        background: [1.0, 1.0, 1.0],
        body: [0.2, 0.25, 0.3],
    });
    let px = |x: usize, y: usize| &img[4 * (x + w * y)..4 * (x + w * y) + 3];
    println!("pixel central {:?}, coin {:?}", px(w / 2, h / 2), px(2, 2));
    assert!(px(w / 2, h / 2).iter().all(|&v| v < 140), "la sphère doit apparaître au centre");
    assert!(img.len() == 4 * w * h && !g.lost);
}

#[test]
fn sphere_stl() {
    // Sphère maillée (icosphère, 5 120 triangles) contre sphère analytique :
    // même voxelisation à une poignée de mailles près, même traînée.
    use lbm::bodies3d::{Body3D, Shape3D};
    use lbm::gpu3d::Gpu3D;
    let l = 8;
    let tris = lbm::stl::icosphere(1.0, [0.0; 3], 4);
    let (mesh, _) = Body3D::from_mesh(&tris, l, &lbm::stl::MeshSetup::default()).unwrap();
    let ana = Body3D::new(Shape3D::Sphere, l, 0.0);
    assert_eq!([mesh.nx, mesh.ny, mesh.nz], [ana.nx, ana.ny, ana.nz]);
    let s = std::f64::consts::PI * 16.0;
    assert!((mesh.area / s - 1.0).abs() < 0.02, "maître-couple {}", mesh.area);
    let u = 0.1f32;
    let nu = u as f64 * l as f64 / 100.0;
    let ramp = lbm::cpu::Ramp { steps: 4 * l as u64 * 10, u };
    let mut cx = Vec::new();
    for b in [&ana, &mesh] {
        let Ok(mut g) = Gpu3D::new(b, nu, u, true) else {
            eprintln!("GPU ignoré");
            return;
        };
        let steps = (30.0 * l as f64 / u as f64) as usize;
        g.run(steps - 800, &ramp);
        let _ = g.take_forces();
        g.run(800, &ramp);
        let f = g.take_forces();
        cx.push(f.fx / (0.5 * (u as f64).powi(2) * b.area));
        println!("{:?} : {} mailles solides, Cx = {:.4}", b.shape, g.solid_count, cx.last().unwrap());
    }
    assert!((cx[1] / cx[0] - 1.0).abs() < 0.03, "Cx maillé {} contre analytique {}", cx[1], cx[0]);
}
