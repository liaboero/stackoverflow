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
    poiseuille_n(slope, c1, hz, 20, 8000, 1.2e-4)
}

/// Canal de Poiseuille de période `nz`, `steps` pas, accélération `a`.
fn poiseuille_n(slope: f64, c1: f64, hz: f64, nz: usize, steps: usize, a: f64) -> Option<(f64, f64)> {
    let nx = 4usize;
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
    let nu = 0.1;
    let lw = if slope == 0.0 { ny as f64 } else { (ny as f64).hypot(nz as f64) };
    let qa = a * h.powi(3) / (12.0 * nu) * lw;
    let umax = a * h * h / (8.0 * nu);
    let body = body_of(tris, [nx, ny, nz], true, true);
    let mut g = gpu(&body, nu)?;
    g.set_periodic_x();
    g.set_accel(a);
    g.run(steps, &Ramp { steps: 0, u: 0.0 });
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

/// Parois à mi-chemin entre deux nœuds, puis inclinées (normale
/// (0, −1, 2)/√5) : la contrainte tangentielle relevée au-dessus des
/// surfels place la paroi à sa position exacte (ordre 2, `surfels_convergence`) ;
/// le schéma d'origine (`SOUFFLERIE_SURFELS_O1=1`) perdait 10 % de débit.
#[test]
fn surfels_poiseuille_incline() {
    let o1 = !lbm::gpu3d::surfels_o2();
    let Some((q, e)) = poiseuille(0.0, 3.5, 12.0) else { return };
    assert!(q.abs() < 0.01 && e < 0.01, "canal droit : débit {q}, profil {e}");
    let Some((q, e)) = poiseuille(0.5, 3.3, 12.0) else { return };
    if o1 {
        assert!(q < 0.0 && q > -0.13 && e < 0.09, "canal incliné : débit {q}, profil {e}");
    } else {
        assert!(q.abs() < 0.015 && e < 0.005, "canal incliné : débit {q}, profil {e}");
    }
}

/// Sphère maillée à Re = 100, 8 mailles par diamètre (comme `sphere_stl`) :
/// avec les surfels, l'effort par échange de quantité de mouvement est celui
/// du bilan sur boîte (à 1 % près), et Cx celui du rebond de Bouzidi (1,236
/// contre 1,236 ; 1,471 avec le schéma d'origine, paroi déplacée vers le
/// fluide).
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
    if lbm::gpu3d::surfels_o2() {
        assert!((cx[0] / cx[1] - 1.0).abs() < 0.03, "Cx surfels {} contre Bouzidi {}", cx[0], cx[1]);
    } else {
        assert!(cx[0] > cx[1] && cx[0] < 1.6, "Cx surfels {} contre Bouzidi {}", cx[0], cx[1]);
    }
}

/// Conduite circulaire (Hagen–Poiseuille) d'axe x, rayon `r`, axe en
/// (yc, zc) = centre + `off` : écarts relatifs du débit Σ V u_x à
/// π a R⁴ / (8 ν) et du profil (nœuds entièrement fluides, L2 rapporté à
/// u_max). Solide : une boîte hors de la grille moins le cylindre (maillage
/// fermé, faces d'extrémité en x inutiles : colonnes selon y ou z).
fn conduite(r: f64, off: [f64; 2], steps: usize) -> Option<(f64, f64)> {
    let nx = 4usize;
    let n = (2.0 * r).ceil() as usize + 5;
    let (yc, zc) = (0.5 * (n - 1) as f64 + off[0], 0.5 * (n - 1) as f64 + off[1]);
    let (x0, x1) = (-1.0, nx as f64 + 1.0);
    let (lo, hi) = (-3.0, n as f64 + 2.0);
    let mut tris: Vec<Tri> = Vec::new();
    // Faces de la boîte, normales sortantes.
    let quad = |a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3], t: &mut Vec<Tri>| {
        t.push([a, b, c]);
        t.push([a, c, d]);
    };
    quad([x0, lo, lo], [x1, lo, lo], [x1, lo, hi], [x0, lo, hi], &mut tris); // y = lo, normale −y
    quad([x0, hi, lo], [x0, hi, hi], [x1, hi, hi], [x1, hi, lo], &mut tris); // y = hi, +y
    quad([x0, lo, lo], [x0, hi, lo], [x1, hi, lo], [x1, lo, lo], &mut tris); // z = lo, −z
    quad([x0, lo, hi], [x1, lo, hi], [x1, hi, hi], [x0, hi, hi], &mut tris); // z = hi, +z
    // Cylindre, normales vers l'axe (sortantes du solide).
    let m = (2.0 * std::f64::consts::PI * r * 2.0).ceil() as usize;
    for k in 0..m {
        let (t0, t1) = (2.0 * std::f64::consts::PI * k as f64 / m as f64, 2.0 * std::f64::consts::PI * (k + 1) as f64 / m as f64);
        let p = |x: f64, t: f64| [x, yc + r * t.cos(), zc + r * t.sin()];
        tris.push([p(x0, t0), p(x1, t0), p(x1, t1)]);
        tris.push([p(x0, t0), p(x1, t1), p(x0, t1)]);
    }
    let nu = 0.1;
    let a = 4.0 * nu * 0.03 / (r * r);
    // Comparaison : rebond de Bouzidi (SOUFFLERIE_TEST_BOUZIDI=1).
    let surf = std::env::var("SOUFFLERIE_TEST_BOUZIDI").is_err();
    let body = body_of(tris, [nx, n, n], surf, true);
    let mut g = gpu(&body, nu)?;
    g.set_periodic_x();
    g.set_accel(a);
    g.run(steps, &Ramp { steps: 0, u: 0.0 });
    assert!(!g.lost);
    let mm = g.macros();
    let vols = g.cell_volumes();
    // Rayon effectif du polygone inscrit : aire π R_e² = m/2 R² sin(2π/m).
    let re = r * ((m as f64 / (2.0 * std::f64::consts::PI)) * (2.0 * std::f64::consts::PI / m as f64).sin()).sqrt();
    let qa = std::f64::consts::PI * a * re.powi(4) / (8.0 * nu);
    let umax = a * re * re / (4.0 * nu);
    let (mut q, mut e2, mut n2) = (0.0, 0.0, 0.0);
    for z in 0..n {
        for y in 0..n {
            let c = nx * (y + n * z);
            q += vols[c] as f64 * mm[c][0] as f64;
            let rr = (y as f64 - yc).hypot(z as f64 - zc);
            if vols[c] == 1.0 && rr < re {
                let ua = a * (re * re - rr * rr) / (4.0 * nu);
                e2 += (mm[c][0] as f64 - ua).powi(2);
                n2 += umax * umax;
            }
        }
    }
    let res = (q / qa - 1.0, (e2 / n2).sqrt());
    println!("conduite R = {r}, décalage {off:?} : débit {:+.3} %, profil {:.3} %", 100.0 * res.0, 100.0 * res.1);
    Some(res)
}

/// Convergence en maillage (paroi volumique) : canal incliné (normale
/// (0, −1, 2)/√5) et conduite circulaire, trois résolutions, moyennes
/// quadratiques des écarts sur plusieurs places de la paroi dans la maille.
/// Ordre mesuré entre résolutions successives (README).
#[test]
fn surfels_convergence() {
    let nu = 0.1;
    let mut rows = Vec::new();
    for k in [1usize, 2, 4] {
        let kf = k as f64;
        let hz = 6.0 * kf;
        let h = hz * 2.0 / 5f64.sqrt();
        let steps = (8.0 * h * h / nu) as usize;
        let a = 1.2e-4 * (12.0 / hz).powi(2);
        let (mut eq, mut ep) = (0.0, 0.0);
        let offs = [0.0, 0.125, 0.25, 0.375];
        for o in offs {
            let Some((q, e)) = poiseuille_n(0.5, 1.6 * kf + o, hz, 10 * k, steps, a) else { return };
            eq += q * q;
            ep += e * e;
        }
        rows.push(("canal", h, (eq / offs.len() as f64).sqrt(), (ep / offs.len() as f64).sqrt()));
    }
    for r in [5.0, 10.0, 20.0] {
        let steps = (2.5 * r * r / nu) as usize;
        let (mut eq, mut ep) = (0.0, 0.0);
        let offs = [[0.0, 0.0], [0.23, 0.11], [0.41, 0.37], [0.13, 0.46]];
        for o in offs {
            let Some((q, e)) = conduite(r, o, steps) else { return };
            eq += q * q;
            ep += e * e;
        }
        rows.push(("conduite", r, (eq / offs.len() as f64).sqrt(), (ep / offs.len() as f64).sqrt()));
    }
    let mut orders = Vec::new();
    for w in rows.windows(2) {
        let (a, b) = (w[0], w[1]);
        if a.0 != b.0 {
            println!("{} H = {:.2} : débit {:.3e}, profil {:.3e}", a.0, a.1, a.2, a.3);
            continue;
        }
        let l = (b.1 / a.1).ln();
        let (oq, op) = ((a.2 / b.2).ln() / l, (a.3 / b.3).ln() / l);
        println!("{} H = {:.2} : débit {:.3e}, profil {:.3e} ; ordre jusqu'à H = {:.2} : débit {oq:.2}, profil {op:.2}", a.0, a.1, a.2, a.3, b.1);
        orders.push((a.0, oq, op));
    }
    let l = rows.last().unwrap();
    println!("{} H = {:.2} : débit {:.3e}, profil {:.3e}", l.0, l.1, l.2, l.3);
    if lbm::gpu3d::surfels_o2() && std::env::var("SOUFFLERIE_TEST_BOUZIDI").is_err() {
        for (cas, oq, _) in orders {
            assert!(oq > 1.7, "{cas} : ordre du débit {oq:.2}");
        }
    }
}
