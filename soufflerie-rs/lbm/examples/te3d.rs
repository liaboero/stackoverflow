//! Diagnostic 3D : écoulement autour du bord de fuite d'une aile STL, à mi-envergure.
//! te3d <stl> <unite> <lref> <sref> <L> <raff> <re> <alpha> <t> <essai> [décalage x] [décalage z]
//!
//! Affiche, dans la tranche y la plus proche du milieu de l'envergure de la
//! grille la plus fine : la place du bord de fuite dans la maille, la carte
//! des nœuds (solide, fluide) autour de lui, Cp, u/U et w/U, et la portance
//! du bilan sur une boîte.
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
        offset: [g(11, 0.0), 0.0, g(12, 0.0)],
        ..Default::default()
    };
    let (l, raff, re, t) = (g(5, 16.0) as usize, g(6, 1.0) as usize, g(7, 1e6), g(9, 20.0));
    let (body, _) = Body3D::from_mesh(&tris, l, &setup).unwrap();
    let u = 0.1f32;
    let mut sim = Gpu3D::with_refinement(&body, u as f64 * l as f64 / re, u, true, true, raff).unwrap();
    sim.set_essai(g(10, 0.0) as u32, 0.0);
    let spt = l as f64 / u as f64;
    let ramp = lbm::cpu::Ramp { steps: (4.0 * spt) as u64, u };
    sim.run((t * spt) as usize, &ramp);
    let q = 0.5 * (u as f64).powi(2);
    let bx = sim.box_forces_vol(&[0.1]);
    if let Some(Some((f, _))) = bx.first() {
        println!("Cz boîte {:.4}", f[2] / (q * body.area));
    }
    let (mac, d, org, sc) = sim.finest_macros();
    // Bord de fuite à mi-envergure : sommet le plus en aval des triangles
    // proches du milieu.
    let m = body.mesh.as_ref().unwrap();
    let ym = 0.5 * (m.lo[1] + m.hi[1]);
    let mut te = [f64::MIN, 0.0, 0.0];
    for t in &m.tris {
        for p in t.iter() {
            if (p[1] - ym).abs() < 0.3 * l as f64 && p[0] > te[0] {
                te = *p;
            }
        }
    }
    let tf = [0, 1, 2].map(|k| (te[k] - org[k]) * sc);
    println!(
        "bord de fuite (grille fine) x = {:.3} (fraction {:.2}), z = {:.3} (fraction {:.2})",
        tf[0],
        tf[0].rem_euclid(1.0),
        tf[2],
        tf[2].rem_euclid(1.0)
    );
    let yj = (ym - org[1]) * sc;
    let yj = yj.round() as usize;
    let xi = tf[0].floor() as i64;
    let zk = tf[2].round() as i64;
    let at = |x: i64, z: i64| mac[x as usize + d[0] * (yj + d[1] * z as usize)];
    // Fenêtre : TE_AMONT mailles en amont du bord de fuite (défaut 5).
    let up: i64 = std::env::var("TE_AMONT").ok().and_then(|v| v.parse().ok()).unwrap_or(5);
    for (name, f) in [
        ("Cp", 0usize),
        ("u/U", 1),
        ("w/U", 2),
        ("perte de pression totale (1 − Cp − V²/U²)", 3),
    ] {
        println!("{name} (x de {} à {}, z décroissant de {} à {}) :", xi - up, xi + 6, zk + 5, zk - 5);
        for z in (zk - 5..=zk + 5).rev() {
            let mut line = format!("z {:+3} |", z - zk);
            for x in xi - up..=xi + 6 {
                let v = at(x, z);
                if v[3] == 0.0 {
                    line += "     ##";
                    continue;
                }
                let val = match f {
                    0 => (v[3] as f64 - 1.0) / 3.0 / q,
                    1 => v[0] as f64 / u as f64,
                    2 => v[2] as f64 / u as f64,
                    _ => {
                        let cp = (v[3] as f64 - 1.0) / 3.0 / q;
                        let v2 = (v[0] as f64).powi(2) + (v[1] as f64).powi(2) + (v[2] as f64).powi(2);
                        1.0 - cp - v2 / (u as f64).powi(2)
                    }
                };
                line += &format!(" {val:+6.2}");
            }
            println!("{line}");
        }
    }
    // Circulation autour d'une tranche : ∮ u·dl sur un rectangle de la
    // grille fine entourant le profil (mailles fluides seulement).
    let xs = (body.bounds.x0 - org[0]) * sc;
    let x0 = xs.floor() as i64 - 3;
    let x1 = xi + 4;
    let (z0, z1) = (zk - (0.2 * l as f64 * sc) as i64, zk + (0.2 * l as f64 * sc) as i64);
    let mut gam = 0.0;
    for x in x0..x1 {
        gam += 0.5 * (at(x, z0)[0] + at(x + 1, z0)[0]) as f64;
        gam -= 0.5 * (at(x, z1)[0] + at(x + 1, z1)[0]) as f64;
    }
    for z in z0..z1 {
        gam += 0.5 * (at(x1, z)[2] + at(x1, z + 1)[2]) as f64;
        gam -= 0.5 * (at(x0, z)[2] + at(x0, z + 1)[2]) as f64;
    }
    // Cl local = 2 Γ / (U c), c = longueur de référence en mailles fines.
    println!("circulation mi-envergure : Cl local = {:.4}", 2.0 * gam / (u as f64 * l as f64 * sc));
}
