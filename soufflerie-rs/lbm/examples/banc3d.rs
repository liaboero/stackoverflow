//! Banc de mesure de la soufflerie 3D : vitesse (mises à jour de maille par
//! seconde, temps de mise en place, de calcul et de relevé) et résultats à
//! temps fixé, à comparer d'une version à l'autre (aucune optimisation ne
//! doit changer les efforts au-delà du bruit d'arrondi).
//!
//!   cargo run --release --example banc3d -- ar4 [--t 6] [--avg 2]
//!       [--sortie ref.txt] [--ref ref.txt] [--champs f.bin] [--profil 1]
//!
//! Cas : `ar4` (aile AR 4 du STL, 4°, Re = 10⁶, L = 16, grille ×2),
//! `ar4b` (idem plus une couche pariétale ×4 de 16 mailles), `ar4r2` (grilles
//! ×2 et ×4), `ar4c3` (grille ×2 et trois couches ×4, ×8, ×16 : `--bande auto`), `ar4n` (sans
//! raffinement), `sphere` (Re = 1000, L = 16, sous-maille sans loi de paroi).
//! Cas proches des calculs de validation : `ahmed` (corps d'Ahmed 25°, sol
//! mobile, L = 64, sans raffinement), `ahmedr` (L = 32, sans sol, grille ×2
//! et couche ×4 de 16 mailles), `aile` (aile infinie NACA 0012, 4°, Re = 6·10⁶,
//! envergure 40 mailles, L = 48, sans raffinement), `ailer` (L = 16, grilles ×2 et ×4, couches ×8
//! et ×16 de 16 mailles : `aile_infinie 16 2 6e6 4 t 16,16`).
//! `--raff-auto N` remplace le raffinement du cas par le raffinement
//! automatique (`lbm::raff_auto`, N mailles par longueur à la paroi).
//! `SOUFFLERIE_RANS=keps|vles` : k-ε RNG ou VLES au lieu de Spalart–Allmaras.
//! La boucle reprend celle de `soufflerie run3d` : efforts relevés tous les
//! quarts de temps convectif, intégration sur la surface et bilan sur boîtes
//! tous les demi-temps convectifs pendant la fenêtre de moyenne.
use std::time::Instant;

use std::sync::Arc;

use lbm::bodies3d::{ahmed, ahmed_tris, Body3D, Shape3D};
use lbm::cpu::Ramp;
use lbm::gpu3d::Gpu3D;
use lbm::shapes::Naca;
use lbm::stl::{MeshBody, MeshSetup, Tri, Unit};

/// Corps d'Ahmed (lunette 25°, pieds, sol mobile), comme l'exemple `ahmed` ;
/// sans sol (veine libre) si `sol` est faux : le raffinement local ne gère
/// pas le sol mobile.
fn ahmed_body(l: usize, marge: f64, sol: bool) -> Body3D {
    let h = ahmed::LONGUEUR / l as f64;
    let tris = ahmed_tris(25.0, true, (2.0 * h).min(10.0), h);
    let setup = MeshSetup {
        unit: Unit::Mm,
        ref_length: Some(ahmed::LONGUEUR),
        ref_area: Some(ahmed::LARGEUR * ahmed::HAUTEUR),
        margin: marge,
        ground: sol.then_some(-h - 0.25 * h),
        ..Default::default()
    };
    Body3D::from_mesh(&tris, l, &setup).expect("maquette").0
}

/// Aile infinie NACA 0012 (veine périodique en envergure), comme l'exemple
/// `aile_infinie` (12 cordes de haut, envergure `ny` mailles).
fn aile_body(l: usize, ny: usize, alpha: f64) -> Body3D {
    let lf = l as f64;
    let (nx, nz) = ((12.0 * lf).round() as usize, (12.0 * lf).round() as usize);
    let naca = Naca::parse("0012").unwrap();
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
    let (sa, ca) = alpha.to_radians().sin_cos();
    let rot = |p: [f64; 2]| {
        let (x, z) = ((p[0] - 0.25) * lf, p[1] * lf);
        [x * ca + z * sa, -x * sa + z * ca]
    };
    let te = rot([1.0, 0.0]);
    let xq = (4.0 * lf + 0.75 * lf).floor() + 0.7 - te[0];
    let zq = 0.5 * nz as f64 + 0.37;
    let pts: Vec<[f64; 2]> = prof.iter().map(|&p| {
        let r = rot(p);
        [xq + r[0], zq + r[1]]
    }).collect();
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
    body.features = None;
    body
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let case = args.first().cloned().unwrap_or_else(|| "ar4".into());
    let opt = |k: &str| args.iter().position(|a| a == &format!("--{k}")).and_then(|i| args.get(i + 1)).cloned();
    let num = |k: &str, d: f64| opt(k).and_then(|v| v.parse().ok()).unwrap_or(d);
    let t_end = num("t", 6.0);
    let avg = num("avg", 2.0);
    const U: f32 = 0.1;
    let l = match case.as_str() {
        "ahmed" => 64usize,
        "aile" => 48,
        "ahmedr" => 32,
        _ => 16,
    };
    let (body, re, depth, bands): (Body3D, f64, usize, Vec<f64>) = match case.as_str() {
        "ahmed" => (ahmed_body(l, 0.4, true), 2.78e6, 0, vec![]),
        "ahmedr" => (ahmed_body(l, 0.4, false), 2.78e6, 1, vec![16.0]),
        "aile" => (aile_body(l, 40, 4.0), 6e6, 0, vec![]),
        "ailer" => (aile_body(l, 40, 4.0), 6e6, 2, vec![16.0, 16.0]),
        "ar4" | "ar4b" | "ar4n" | "ar4r2" | "ar4c3" => {
            let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../banc/wing_ar4.stl");
            let tris = lbm::stl::read_stl(&path).expect("banc/wing_ar4.stl");
            let setup = MeshSetup {
                unit: Unit::M,
                alpha_deg: 4.0,
                ref_length: Some(1.0),
                ref_area: Some(4.0),
                margin: 0.4,
                ..Default::default()
            };
            let (b, _) = Body3D::from_mesh(&tris, l, &setup).expect("maquette");
            let (depth, bands) = match case.as_str() {
                "ar4" => (1, vec![]),
                "ar4b" => (1, vec![16.0]),
                "ar4r2" => (2, vec![]),
                "ar4c3" => (1, vec![16.0, 16.0, 16.0]),
                _ => (0, vec![]),
            };
            (b, 1e6, depth, bands)
        }
        "sphere" => (Body3D::new(Shape3D::Sphere, l, 0.0), 1000.0, 0, vec![]),
        _ => {
            eprintln!("cas : ar4, ar4b, ar4n, ar4r2, ar4c3, sphere, ahmed, ahmedr, aile, ailer");
            std::process::exit(2)
        }
    };
    let mesh = body.mesh.is_some();
    let nu = U as f64 * l as f64 / re;
    // Modèle RANS : Spalart–Allmaras, ou SOUFFLERIE_RANS=keps|vles.
    let sa = if re >= lbm::tunnel::RE_SA_AUTO { lbm::gpu3d::Turbulence::rans_from_env() } else { lbm::gpu3d::Turbulence::Laminar };
    let t0 = Instant::now();
    let mut sim = match lbm::raff_auto::options_from_args(&args) {
        Some(o) => lbm::raff_auto::build(&body, nu, U, true, sa, &o).expect("GPU"),
        None => Gpu3D::with_refinement_edges(&body, nu, U, true, sa, depth, false, false, &bands).expect("GPU"),
    };
    let prep = t0.elapsed().as_secs_f64();
    eprintln!("veine {} × {} × {}", sim.nx, sim.ny, sim.nz);
    let cells = sim.grid_cells();
    let per_step: f64 = cells.iter().map(|&(n, d)| (n as f64) * (1u64 << d) as f64).sum();
    eprintln!(
        "{case} : {} grilles, {:.2} M mailles, {:.2} M mises à jour par pas de la veine · {}",
        cells.len(),
        cells.iter().map(|c| c.0).sum::<usize>() as f64 / 1e6,
        per_step / 1e6,
        sim.adapter_name()
    );
    let spt = l as f64 / U as f64;
    let ramp = Ramp { steps: (4.0 * spt).round() as u64, u: U };
    let q = 0.5 * (U as f64).powi(2) * body.area;
    let chunk = (0.25 * spt).round().max(1.0) as usize;
    let total = (t_end * spt).round() as u64;
    let mut next_surf = t_end - avg;
    let (mut t_run, mut t_read) = (0.0f64, 0.0f64);
    let (mut fx, mut fz, mut nf) = (0.0, 0.0, 0usize);
    let mut last = (0.0, 0.0, 0.0);
    let mut surf = Vec::new();
    let mut boxes: Vec<[f64; 6]> = Vec::new();
    let t_loop = Instant::now();
    while sim.time < total {
        let ts = Instant::now();
        sim.run(chunk, &ramp);
        let f = sim.take_forces();
        t_run += ts.elapsed().as_secs_f64();
        if !sim.is_finite(&f) {
            eprintln!("divergence");
            std::process::exit(1);
        }
        let t = sim.time as f64 / spt;
        last = (f.fx / q, f.fy / q, f.fz / q);
        if t >= t_end - avg {
            fx += f.fx / q;
            fz += f.fz / q;
            nf += 1;
        }
        if mesh && t >= next_surf {
            let tr = Instant::now();
            if let Some(s) = sim.surface_forces() {
                surf.push(s);
            }
            if let [Some((f1, v1)), Some((f2, v2))] = sim.box_forces_vol(&[0.1, 0.15])[..] {
                boxes.push([f1[0], f1[2], v1, f2[0], f2[2], v2]);
            }
            t_read += tr.elapsed().as_secs_f64();
            next_surf += 0.5;
        }
    }
    let loop_s = t_loop.elapsed().as_secs_f64();
    let steps = sim.time as f64;
    let mut res: Vec<(String, f64)> = vec![
        ("cx_echange".into(), fx / nf.max(1) as f64),
        ("cz_echange".into(), fz / nf.max(1) as f64),
        ("cx_dernier".into(), last.0),
        ("cy_dernier".into(), last.1),
        ("cz_dernier".into(), last.2),
    ];
    if !surf.is_empty() {
        let n = surf.len() as f64;
        let m = |f: &dyn Fn(&lbm::gpu3d::SurfaceForces) -> f64| surf.iter().map(f).sum::<f64>() / n / q;
        res.push(("cx_surface".into(), m(&|s| s.pressure_wall[0] + s.friction[0])));
        res.push(("cz_surface".into(), m(&|s| s.pressure_wall[2] + s.friction[2])));
        res.push(("cx_frottement".into(), m(&|s| s.friction[0])));
    }
    if !boxes.is_empty() {
        let nb = boxes.len() as f64;
        let b = |a: usize| boxes.iter().map(|f| f[a]).sum::<f64>() / nb;
        let (cx1, cx2, v1, v2) = (b(0) / q, b(3) / q, b(2), b(5));
        res.push(("cx_boite".into(), cx1 - v1 * (cx2 - cx1) / (v2 - v1)));
        res.push(("cz_boite".into(), 0.5 * (b(1) + b(4)) / q));
    }
    // Champ de la grille entière la plus fine : sommes de contrôle, et
    // fichier brut pour comparer maille par maille.
    let (mac, _, _, _) = sim.finest_macros();
    let mut sum = [0.0f64; 4];
    let mut sq = [0.0f64; 4];
    for c in &mac {
        for k in 0..4 {
            sum[k] += c[k] as f64;
            sq[k] += (c[k] as f64).powi(2);
        }
    }
    for (k, nm) in ["ux", "uy", "uz", "rho"].iter().enumerate() {
        res.push((format!("champ_somme_{nm}"), sum[k]));
        res.push((format!("champ_carres_{nm}"), sq[k]));
    }
    if let Some(path) = opt("champs") {
        let bytes: Vec<u8> = mac.iter().flat_map(|c| c.iter().flat_map(|v| v.to_le_bytes())).collect();
        if let Some(r) = opt("champs-ref") {
            let old = std::fs::read(&r).expect("champ de référence");
            let a: Vec<f32> = old.chunks_exact(4).map(|b| f32::from_le_bytes([b[0], b[1], b[2], b[3]])).collect();
            let bv: Vec<f32> = mac.iter().flatten().copied().collect();
            if a.len() == bv.len() {
                let mut dmax = [0.0f32; 4];
                let mut ndiff = 0usize;
                for (i, (x, y)) in a.iter().zip(&bv).enumerate() {
                    let d = (x - y).abs();
                    if d > 0.0 || x.to_bits() != y.to_bits() {
                        ndiff += 1;
                    }
                    dmax[i % 4] = dmax[i % 4].max(d);
                }
                eprintln!(
                    "champ : {ndiff} valeurs différentes sur {} · écart max u {:.2e} {:.2e} {:.2e} ρ {:.2e}",
                    a.len(),
                    dmax[0],
                    dmax[1],
                    dmax[2],
                    dmax[3]
                );
            } else {
                eprintln!("champ : tailles différentes ({} et {})", a.len(), bv.len());
            }
        }
        std::fs::write(&path, bytes).expect("écriture du champ");
    }
    let mlups = per_step * steps / t_run / 1e6;
    println!(
        "{case} · t* = {t_end} ({steps} pas) · préparation {prep:.2} s · calcul {t_run:.2} s ({mlups:.0} MLUPS) · relevés {t_read:.2} s · boucle {loop_s:.2} s"
    );
    let reference: Vec<(String, f64)> = opt("ref")
        .map(|p| {
            std::fs::read_to_string(&p)
                .expect("référence")
                .lines()
                .filter_map(|l| {
                    let mut it = l.split_whitespace();
                    Some((it.next()?.to_string(), it.next()?.parse().ok()?))
                })
                .collect()
        })
        .unwrap_or_default();
    let text: String = res.iter().map(|(k, v)| format!("{k} {v:.10e}\n")).collect();
    if let Some(p) = opt("sortie") {
        std::fs::write(p, &text).expect("écriture");
    }
    for (k, v) in &res {
        match reference.iter().find(|(r, _)| r == k) {
            Some((_, r)) => {
                let d = if *r != 0.0 { (v - r) / r.abs() } else { v - r };
                println!("  {k:18} {v:14.8} réf. {r:14.8} écart relatif {d:+.2e}");
            }
            None => println!("  {k:18} {v:14.8}"),
        }
    }
    if num("profil", 0.0) > 0.0 {
        let reps = num("reps", 50.0) as usize;
        let mut tot = 0.0;
        println!("Profil (ms par lancement, {reps} lancements) :");
        for (name, n, ms) in sim.profile_passes(reps) {
            println!("  {name:32} {:8.3} M  {ms:8.3} ms  {:7.0} M/s", n as f64 / 1e6, n as f64 / ms / 1e3);
            tot += ms;
        }
        println!("  somme des passes (un lancement chacune) : {tot:.3} ms");
    }
}
