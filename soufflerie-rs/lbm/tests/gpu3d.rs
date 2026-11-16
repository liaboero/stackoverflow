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
        body_view: Default::default(),
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

/// Pavé [x0, x1] × [y0, y1] × [z0, z1] en triangles.
fn pave(x: [f64; 2], y: [f64; 2], z: [f64; 2]) -> Vec<lbm::stl::Tri> {
    let p = |i: usize, j: usize, k: usize| [x[i], y[j], z[k]];
    let faces = [
        [p(0, 0, 0), p(0, 1, 0), p(0, 1, 1), p(0, 0, 1)],
        [p(1, 0, 0), p(1, 0, 1), p(1, 1, 1), p(1, 1, 0)],
        [p(0, 0, 0), p(0, 0, 1), p(1, 0, 1), p(1, 0, 0)],
        [p(0, 1, 0), p(1, 1, 0), p(1, 1, 1), p(0, 1, 1)],
        [p(0, 0, 0), p(1, 0, 0), p(1, 1, 0), p(0, 1, 0)],
        [p(0, 0, 1), p(0, 1, 1), p(1, 1, 1), p(1, 0, 1)],
    ];
    faces.iter().flat_map(|[a, b, c, d]| [[*a, *b, *c], [*a, *c, *d]]).collect()
}

#[test]
fn paroi_mince_stl() {
    // Plaque carrée face au vent, 0,08 maille d'épaisseur (aucune maille
    // solide) : liens de paroi mince, même traînée qu'une plaque épaisse.
    use lbm::bodies3d::Body3D;
    use lbm::gpu3d::Gpu3D;
    let l = 10;
    let setup = lbm::stl::MeshSetup {
        ref_length: Some(100.0),
        ..Default::default()
    };
    let u = 0.1f32;
    let nu = u as f64 * l as f64 / 50.0;
    let ramp = lbm::cpu::Ramp { steps: 4 * l as u64 * 10, u };
    let mut cx = Vec::new();
    for t in [0.8, 20.0] {
        let tris = pave([0.0, t], [-50.0, 50.0], [-50.0, 50.0]);
        let (b, _) = Body3D::from_mesh(&tris, l, &setup).unwrap();
        let Ok(mut g) = Gpu3D::new(&b, nu, u, true) else {
            eprintln!("GPU ignoré");
            return;
        };
        let steps = (25.0 * l as f64 / u as f64) as usize;
        g.run(steps - 500, &ramp);
        let _ = g.take_forces();
        g.run(500, &ramp);
        let f = g.take_forces();
        cx.push(f.fx / (0.5 * (u as f64).powi(2) * b.area));
        println!("épaisseur {t} mm : {} mailles solides, {} pariétales, Cx = {:.4}", g.solid_count, g.wall_count, cx.last().unwrap());
        if t < 1.0 {
            assert_eq!(g.solid_count, 0);
        }
    }
    assert!((cx[0] / cx[1] - 1.0).abs() < 0.06, "plaque mince {} contre épaisse {}", cx[0], cx[1]);
}

#[test]
fn loi_de_paroi_3d() {
    // Aile NACA 0012 d'allongement 4 à 8°, Re = 10⁶ : sans modèle de paroi, la
    // sous-maille seule écrase la portance ; Spalart–Allmaras et la loi de
    // Spalding la ramènent près de la ligne portante (CL ≈ 0,56).
    use lbm::bodies3d::{Body3D, Shape3D};
    use lbm::gpu3d::Gpu3D;
    let l = 12;
    let u = 0.1f32;
    let spt = l as f64 / u as f64;
    let body = Body3D::new(Shape3D::Wing, l, 8.0);
    let ramp = lbm::cpu::Ramp { steps: (4.0 * spt) as u64, u };
    let mut cz = Vec::new();
    for sa in [false, true] {
        let Ok(mut g) = Gpu3D::with_turbulence(&body, u as f64 * l as f64 / 1e6, u, true, sa) else {
            eprintln!("GPU ignoré");
            return;
        };
        g.run((20.0 * spt) as usize, &ramp);
        let _ = g.take_forces();
        g.run((6.0 * spt) as usize, &ramp);
        let f = g.take_forces();
        assert!(g.is_finite(&f));
        cz.push(f.fz / (0.5 * (u as f64).powi(2) * body.area));
        println!("SA {sa} : Cz = {:.4}", cz.last().unwrap());
    }
    assert!(cz[1] > 2.0 * cz[0] && (0.35..0.75).contains(&cz[1]), "Cz {cz:?}");
}

#[test]
fn raffinement_3d() {
    // Sphère à Re = 100 : veine D = 8 avec une grille fine ×2 autour de la
    // sphère, contre veine uniforme D = 16.
    use lbm::bodies3d::{Body3D, Shape3D};
    use lbm::gpu3d::Gpu3D;
    let u = 0.1f32;
    let mut cx = Vec::new();
    for (l, depth) in [(16usize, 0usize), (8, 1)] {
        let body = Body3D::new(Shape3D::Sphere, l, 0.0);
        let spt = l as f64 / u as f64;
        let Ok(mut g) = Gpu3D::with_refinement(&body, u as f64 * l as f64 / 100.0, u, true, false, depth) else {
            eprintln!("GPU ignoré");
            return;
        };
        let ramp = lbm::cpu::Ramp { steps: (4.0 * spt) as u64, u };
        g.run((20.0 * spt) as usize, &ramp);
        let _ = g.take_forces();
        g.run((4.0 * spt) as usize, &ramp);
        let f = g.take_forces();
        cx.push(f.fx / (0.5 * (u as f64).powi(2) * body.area));
        println!("D = {l}, raffinement {depth} : Cx = {:.4}", cx.last().unwrap());
    }
    assert!((cx[1] / cx[0] - 1.0).abs() < 0.015, "Cx {cx:?}");
}

#[test]
fn transition_3d() {
    // Aile à Re = 10⁶ avec loi de paroi (L = 16, α = 8°) : une couche limite
    // laminaire sur toute la corde retire environ l'écart de frottement
    // Blasius / turbulent (≈ 0,006). Sur une grille plus grossière, le
    // glissement laminaire change aussi la traînée de pression (≈ 0,013 avec
    // le bord de fuite calé à 0,7 maille après un nœud, ≈ 0,017 avec la loi de
    // paroi alimentée une maille plus loin).
    use lbm::bodies3d::{Body3D, Shape3D};
    use lbm::gpu3d::Gpu3D;
    let l = 16;
    let u = 0.1f32;
    let spt = l as f64 / u as f64;
    let body = Body3D::new(Shape3D::Wing, l, 8.0);
    let mut cx = Vec::new();
    for xtr in [0.0, 1.0] {
        let Ok(mut g) = Gpu3D::with_turbulence(&body, u as f64 * l as f64 / 1e6, u, true, true) else {
            eprintln!("GPU ignoré");
            return;
        };
        g.set_transition(xtr);
        let ramp = lbm::cpu::Ramp { steps: (4.0 * spt) as u64, u };
        g.run((20.0 * spt) as usize, &ramp);
        let _ = g.take_forces();
        g.run((6.0 * spt) as usize, &ramp);
        let f = g.take_forces();
        cx.push(f.fx / (0.5 * (u as f64).powi(2) * body.area));
        println!("xtr {xtr} : Cx = {:.4}", cx.last().unwrap());
    }
    let d = cx[0] - cx[1];
    assert!((0.002..0.02).contains(&d), "baisse de traînée {d}");
}
