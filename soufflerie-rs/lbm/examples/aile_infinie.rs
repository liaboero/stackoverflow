//! Diagnostic 3D : aile NACA 0012 d'envergure infinie (veine périodique en y),
//! contre les références 2D (TMR, Ladson). Sans traînée induite : la traînée
//! est celle du profil, pression et frottement.
//!
//! aile_infinie <L> <raff> <Re> <alpha> <t> [bandes, ex. 16,16 ou -] [ny] [hauteur en cordes] [essai]
//!
//! L mailles par corde (veine), `raff` grilles fines entières, `bandes` :
//! couches pariétales imbriquées (comme `run3d --bande`). Efforts intégrés
//! sur la surface (triangles de la tranche centrale de 4 mailles : les bords
//! de la grille fine en envergure ne sont pas périodiques et faussent les
//! tranches voisines) et échange aux parois (toute la veine, indicatif),
//! moyennés sur la dernière corde parcourue (4 relevés).
use lbm::bodies3d::Body3D;
use lbm::gpu3d::Gpu3D;
use lbm::shapes::Naca;
use lbm::stl::{MeshBody, MeshSetup, Tri, Unit};
use std::sync::Arc;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let g = |i: usize, d: f64| a.get(i).and_then(|v| v.parse().ok()).unwrap_or(d);
    let l = g(1, 16.0) as usize;
    let raff = g(2, 2.0) as usize;
    let re = g(3, 6e6);
    let alpha = g(4, 0.0);
    let t_end = g(5, 12.0);
    let bands: Vec<f64> = a.get(6).map(|s| s.split(',').filter_map(|v| v.parse().ok()).collect()).unwrap_or_default();
    let ny = g(7, 8.0) as usize;
    let haut = g(8, 12.0);
    let essai = g(9, 0.0) as u32;
    let lf = l as f64;

    // Veine : 4 cordes en amont, 8 en aval, `haut` cordes de haut.
    let (nx, nz) = ((12.0 * lf).round() as usize, (haut * lf).round() as usize);
    let naca = Naca::parse("0012").unwrap();
    // Contour (x/c, z/c) en espacement cosinus, extrados du bord de fuite au
    // bord d'attaque puis intrados ; bord de fuite fermé.
    let n = 120;
    let mut prof: Vec<[f64; 2]> = Vec::new();
    for k in 0..n {
        let x = 0.5 * (1.0 + (std::f64::consts::PI * k as f64 / n as f64).cos());
        prof.push([x, naca.thickness(x)]);
    }
    for k in 0..n {
        let x = 0.5 * (1.0 - (std::f64::consts::PI * k as f64 / n as f64).cos());
        prof.push([x, -naca.thickness(x)]);
    }
    // Rotation autour du quart de corde (bord d'attaque relevé pour α > 0),
    // en mailles.
    let (sa, ca) = alpha.to_radians().sin_cos();
    let rot = |p: [f64; 2]| {
        let (x, z) = ((p[0] - 0.25) * lf, p[1] * lf);
        [x * ca + z * sa, -x * sa + z * ca]
    };
    let te = rot([1.0, 0.0]);
    // Bord de fuite à 0,7 maille après un nœud (plateau de portance), centre
    // en hauteur au milieu de la veine.
    let xq = (4.0 * lf + 0.75 * lf).floor() + 0.7 - te[0];
    let zq = 0.5 * nz as f64 + 0.37 + std::env::var("AILE_DZ").ok().and_then(|v| v.parse::<f64>().ok()).unwrap_or(0.0);
    let xq = xq + std::env::var("AILE_DX").ok().and_then(|v| v.parse::<f64>().ok()).unwrap_or(0.0);
    let pts: Vec<[f64; 2]> = prof.iter().map(|&p| {
        let r = rot(p);
        [xq + r[0], zq + r[1]]
    }).collect();
    // Extrusion en stations d'envergure de 2 mailles, de −4 à ny + 4 : les
    // triangles dont le centre tombe dans [0, ny) couvrent exactement la veine.
    let ys: Vec<f64> = (-2..=(ny as i64 / 2 + 2)).map(|k| 2.0 * k as f64).collect();
    let mut tris: Vec<Tri> = Vec::new();
    let m = pts.len();
    for w in ys.windows(2) {
        let (y0, y1) = (w[0], w[1]);
        for k in 0..m {
            let (p, q) = (pts[k], pts[(k + 1) % m]);
            let (a0, b0, a1, b1) = ([p[0], y0, p[1]], [q[0], y0, q[1]], [p[0], y1, p[1]], [q[0], y1, q[1]]);
            tris.push([a0, b0, b1]);
            tris.push([a0, b1, a1]);
        }
    }
    let setup = MeshSetup { unit: Unit::parse("m").unwrap(), margin: 0.2, ..Default::default() };
    let (mut body, _) = Body3D::from_mesh(&tris, l, &setup).unwrap();
    let mesh = MeshBody::new(tris, false);
    body.nx = nx;
    body.ny = ny;
    body.nz = nz;
    body.bounds = lbm::bodies3d::Bounds { x0: mesh.lo[0], x1: mesh.hi[0], y0: 0.0, y1: ny as f64 - 1.0, z0: mesh.lo[2], z1: mesh.hi[2] };
    body.center = [xq, 0.5 * ny as f64, zq];
    body.area = lf * ny as f64;
    body.mesh = Some(Arc::new(mesh));

    let u = 0.1f32;
    let nu = u as f64 * lf / re;
    let spt = lf / u as f64;
    let ramp = lbm::cpu::Ramp { steps: (2.0 * spt) as u64, u };
    let t0 = std::time::Instant::now();
    let mut sim = Gpu3D::with_refinement_edges(&body, nu, u, true, true, raff, false, false, &bands).unwrap();
    sim.periodic_span = true;
    if essai != 0 {
        sim.set_essai(essai, 0.0);
    }
    sim.run(((t_end - 1.0) * spt) as usize, &ramp);
    let _ = sim.take_forces();
    let q = 0.5 * (u as f64).powi(2);
    let nrel = 4;
    let (mut cx_p, mut cx_f, mut cz, mut cm_n) = (0.0, 0.0, 0.0, 0usize);
    for _ in 0..nrel {
        sim.run((spt / nrel as f64) as usize, &ramp);
        let Some(s) = sim.surface_forces() else { continue };
        // Pression et frottement séparés : `loads` les mêlent, on les
        // recompose par la part de chaque triangle dans la veine.
        // Seuls les triangles couverts par une grille fine sont relevés :
        // effort par unité d'envergure sur les bandes relevées (écoulement
        // plan), chaque bande de 2 mailles comptant 2 m triangles.
        let (mut fx, mut fz, mut cnt) = (0.0, 0.0, 0usize);
        // Tranche centrale de 4 mailles seulement : les bords de la grille
        // fine en envergure (non périodiques) faussent les tranches voisines.
        let (yc0, yc1) = (2.0 * (ny / 4) as f64 - 2.0, 2.0 * (ny / 4) as f64 + 2.0);
        for (p, f) in &s.loads {
            if p[1] >= yc0 && p[1] < yc1 {
                fx += f[0];
                fz += f[2];
                cnt += 1;
            }
        }
        let span = cnt as f64 / m as f64;
        let qa = q * lf * span.max(1e-9);
        // Part de frottement : celle du total, ramenée à la veine.
        let tot = s.pressure_wall[0] + s.friction[0];
        let fr = if tot.abs() > 0.0 { s.friction[0] / tot } else { 0.0 };
        cx_p += fx * (1.0 - fr) / qa;
        cx_f += fx * fr / qa;
        cz += fz / qa;
        cm_n += 1;
        if std::env::var("AILE_BANDES").is_ok() {
            // Cx et Cz par bande d'envergure de 2 mailles (veine).
            let nb = ny / 2;
            let mut b = vec![[0.0f64; 3]; nb];
            for (p, f) in &s.loads {
                if p[1] >= 0.0 && p[1] < ny as f64 {
                    let k = ((p[1] / 2.0) as usize).min(nb - 1);
                    b[k][0] += f[0];
                    b[k][1] += f[2];
                    b[k][2] += 1.0;
                }
            }
            let v: Vec<String> = b.iter().map(|x| if x[2] > 0.0 { format!("{:.4}/{:.3}", x[0] / (q * lf * 2.0), x[1] / (q * lf * 2.0)) } else { "-".into() }).collect();
            println!("  bandes : {}", v.join(" "));
        }
        if std::env::var("AILE_DEBUG").is_ok() {
            println!(
                "  totaux : pression_wall {:?} frottement {:?} manqué {:.3} · {} triangles relevés",
                s.pressure_wall.map(|v| v / qa),
                s.friction.map(|v| v / qa),
                s.missed,
                s.loads.len()
            );
        }
    }
    let f = sim.take_forces();
    let nn = cm_n.max(1) as f64;
    let (mac, ..) = sim.finest_macros();
    let ok = mac.iter().all(|m| m[0].is_finite());
    println!(
        "L = {l} raff {raff} bandes {:?} Re = {re:e} α = {alpha}° ny = {ny} haut {haut} : surface Cz {:.4} Cx {:.5} (pression {:.5} + frottement {:.5}) · échange Cz {:.4} Cx {:.5}{} · {:.0} s",
        bands,
        cz / nn,
        (cx_p + cx_f) / nn,
        cx_p / nn,
        cx_f / nn,
        f.fz / (q * body.area),
        f.fx / (q * body.area),
        if ok { "" } else { " · DIVERGÉ" },
        t0.elapsed().as_secs_f64()
    );
}
