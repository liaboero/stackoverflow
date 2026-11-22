//! Paroi volumique par surfels (`Body3D::surfels`, ignorée sans adaptateur
//! GPU) : conservation de la masse et écoulement de Poiseuille dans un canal
//! incliné sur la grille.
#![cfg(feature = "gpu")]

use std::sync::Arc;

use lbm::bodies3d::{Body3D, Bounds, Shape3D};
use lbm::cpu::Ramp;
use lbm::gpu3d::{Gpu3D, Turbulence};
use lbm::stl::{MeshBody, Tri};

/// Maquette maillée en coordonnées de la grille (pas de mise à l'échelle).
fn body_of(tris: Vec<Tri>, dims: [usize; 3], surfels: bool, robust: bool) -> Body3D {
    let mut b = Body3D::new(Shape3D::Sphere, 8, 0.0);
    let mesh = MeshBody::new(tris, robust);
    b.shape = Shape3D::Mesh;
    b.nx = dims[0];
    b.ny = dims[1];
    b.nz = dims[2];
    b.bounds = Bounds { x0: mesh.lo[0], x1: mesh.hi[0], y0: mesh.lo[1], y1: mesh.hi[1], z0: mesh.lo[2], z1: mesh.hi[2] };
    b.center = [0, 1, 2].map(|a| 0.5 * (mesh.lo[a] + mesh.hi[a]));
    b.mesh = Some(Arc::new(mesh));
    b.features = None;
    b.surfels = surfels;
    b
}

fn gpu(body: &Body3D, nu: f64) -> Option<Gpu3D> {
    match Gpu3D::with_turbulence(body, nu, 0.0, false, Turbulence::Laminar) {
        Ok(g) => Some(g),
        Err(e) => {
            eprintln!("GPU ignoré : {e}");
            None
        }
    }
}

/// Sphère maillée dans une veine périodique, écoulement entretenu par une
/// force volumique : dérive relative de la masse Σ V ρ en 4 000 pas, contre
/// celle de la même veine vide (arrondi de la simple précision, ≈ 1,5·10⁻⁸
/// par pas). Les surfels n'y ajoutent rien ; le rebond de Bouzidi, si.
#[test]
fn surfels_masse() {
    let dims = [40, 32, 32];
    let mut out = Vec::new();
    // Veine vide (maillage hors de la grille), surfels, Bouzidi.
    for (dx, surf) in [(1000.0, true), (0.0, true), (0.0, false)] {
        let tris = lbm::stl::icosphere(6.3, [17.3 + dx, 16.2, 15.9], 3);
        let body = body_of(tris, dims, surf, false);
        let Some(mut g) = gpu(&body, 0.02) else { return };
        g.set_periodic_x();
        g.set_accel(2e-5);
        let ramp = Ramp { steps: 0, u: 0.0 };
        g.run(1, &ramp);
        let m0 = g.fluid_mass();
        g.run(4000, &ramp);
        let m1 = g.fluid_mass();
        let f = g.take_forces();
        let drift = m1 / m0 - 1.0;
        println!("vide {} surfels {surf} : masse {m0:.3} → {m1:.3} (dérive {drift:.2e}), Fx = {:.3e}, {:?}", dx != 0.0, f.fx, g.surfel_info());
        assert!(!g.lost && f.fx.is_finite());
        out.push(drift);
    }
    println!("écart à la veine vide : surfels {:.2e}, Bouzidi {:.2e}", out[1] - out[0], out[2] - out[0]);
    assert!((out[1] - out[0]).abs() < 1e-5, "surfels : dérive de masse {:.2e} contre {:.2e} sans maquette", out[1], out[0]);
}

/// Canal plan de Poiseuille, fluide c₁ < z − pente·y < c₁ + largeur (mod
/// nz = 20), écoulement selon x entretenu par une force volumique : écarts
/// relatifs du débit Σ V u_x à a H³ L / (12 ν) et du profil (nœuds des voxels
/// entièrement fluides, L2 rapporté à u_max). None sans GPU.
fn poiseuille(slope: f64, c1: f64, hz: f64) -> Option<(f64, f64)> {
    let (nx, nz) = (4usize, 20usize);
    let ny = if slope == 0.0 { 4 } else { (nz as f64 / slope).round() as usize };
    let c2 = c1 + hz;
    let cs = 1.0 / (1.0 + slope * slope).sqrt();
    let h = hz * cs;
    let mut tris: Vec<Tri> = Vec::new();
    let (x0, x1, y0, y1) = (-1.0, nx as f64 + 1.0, -3.0, ny as f64 + 3.0);
    let plane = |c: f64, up: bool, t: &mut Vec<Tri>| {
        let p = |x: f64, y: f64| [x, y, slope * y + c];
        let (a, b, cc, d) = (p(x0, y0), p(x1, y0), p(x1, y1), p(x0, y1));
        // Normale sortante du solide : vers +z (face haute) ou −z (face basse).
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
        // Solide entre le haut du canal et le bas du canal suivant.
        plane(c2 + k, false, &mut tris);
        plane(c1 + nz as f64 + k, true, &mut tris);
    }
    let (nu, a) = (0.1, 1.2e-4);
    let lw = if slope == 0.0 { ny as f64 } else { (ny as f64).hypot(nz as f64) };
    let qa = a * h.powi(3) / (12.0 * nu) * lw;
    let umax = a * h * h / (8.0 * nu);
    let body = body_of(tris, [nx, ny, nz], true, true);
    let mut g = gpu(&body, nu)?;
    g.set_periodic_x();
    g.set_accel(a);
    g.run(8000, &Ramp { steps: 0, u: 0.0 });
    assert!(!g.lost);
    let m = g.macros();
    let vols = g.cell_volumes();
    let (mut q, mut e2, mut n2) = (0.0, 0.0, 0.0);
    for z in 0..nz {
        for y in 0..ny {
            let c = nx * (y + ny * z);
            q += vols[c] as f64 * m[c][0] as f64;
            let eta = ((z as f64 - slope * y as f64 - c1).rem_euclid(nz as f64)) * cs;
            if vols[c] == 1.0 && eta < h {
                let ua = a * eta * (h - eta) / (2.0 * nu);
                e2 += (m[c][0] as f64 - ua).powi(2);
                n2 += umax * umax;
            }
        }
    }
    let r = (q / qa - 1.0, (e2 / n2).sqrt());
    println!("Poiseuille pente {slope}, c₁ = {c1} : débit {:+.2} %, profil {:.2} %", 100.0 * r.0, 100.0 * r.1);
    Some(r)
}

/// Parois à mi-chemin entre deux nœuds : le rebond sur les surfels est le
/// rebond à mi-chemin, exact à l'ordre 2. Paroi inclinée (normale
/// (0, −1, 2)/√5) : la paroi volumique est d'ordre 1 près de la paroi (voxels
/// coupés mélangés, README) ; le débit perd environ 10 % à H ≈ 10,7 mailles.
#[test]
fn surfels_poiseuille_incline() {
    let Some((q, e)) = poiseuille(0.0, 3.5, 12.0) else { return };
    assert!(q.abs() < 0.01 && e < 0.01, "canal droit : débit {q}, profil {e}");
    let Some((q, e)) = poiseuille(0.5, 3.3, 12.0) else { return };
    assert!(q < 0.0 && q > -0.13 && e < 0.09, "canal incliné : débit {q}, profil {e}");
}

/// Sphère maillée à Re = 100, 8 mailles par diamètre (comme `sphere_stl`) :
/// avec les surfels, l'effort par échange de quantité de mouvement est celui
/// du bilan sur boîte (à 1 % près). Cx ≈ 1,47 contre 1,24 avec Bouzidi
/// (1,26 dans l'application web) : la paroi volumique se comporte comme une
/// paroi déplacée d'environ 0,2 maille vers le fluide (README).
#[test]
fn surfels_sphere_re100() {
    let l = 8;
    let tris = lbm::stl::icosphere(1.0, [0.0; 3], 4);
    let u = 0.1f32;
    let nu = u as f64 * l as f64 / 100.0;
    let ramp = Ramp { steps: 4 * l as u64 * 10, u };
    let mut cx = Vec::new();
    let mut bil = Vec::new();
    for surf in [true, false] {
        let (mut b, _) = Body3D::from_mesh(&tris, l, &lbm::stl::MeshSetup::default()).unwrap();
        b.surfels = surf;
        let Ok(mut g) = Gpu3D::new(&b, nu, u, true) else {
            eprintln!("GPU ignoré");
            return;
        };
        let steps = (30.0 * l as f64 / u as f64) as usize;
        g.run(steps - 800, &ramp);
        let _ = g.take_forces();
        g.run(800, &ramp);
        let f = g.take_forces();
        let q = 0.5 * (u as f64).powi(2) * b.area;
        let bx = g.box_forces(0.5).map_or(f64::NAN, |v| v[0] / q);
        cx.push(f.fx / q);
        bil.push(bx);
        println!("surfels {surf} : Cx échange {:.4}, bilan sur boîte {bx:.4}, {} mailles pariétales", f.fx / q, g.wall_count);
        assert!(!g.lost);
    }
    assert!((cx[0] / bil[0] - 1.0).abs() < 0.01, "surfels : échange {} contre bilan {}", cx[0], bil[0]);
    assert!(cx[0] > cx[1] && cx[0] < 1.6, "Cx surfels {} contre Bouzidi {}", cx[0], cx[1]);
}
