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
//! ×2 et ×4), `ar4n` (sans
//! raffinement), `sphere` (Re = 1000, L = 16, sous-maille sans loi de paroi).
//! La boucle reprend celle de `soufflerie run3d` : efforts relevés tous les
//! quarts de temps convectif, intégration sur la surface et bilan sur boîtes
//! tous les demi-temps convectifs pendant la fenêtre de moyenne.
use std::time::Instant;

use lbm::bodies3d::{Body3D, Shape3D};
use lbm::cpu::Ramp;
use lbm::gpu3d::Gpu3D;
use lbm::stl::{MeshSetup, Unit};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let case = args.first().cloned().unwrap_or_else(|| "ar4".into());
    let opt = |k: &str| args.iter().position(|a| a == &format!("--{k}")).and_then(|i| args.get(i + 1)).cloned();
    let num = |k: &str, d: f64| opt(k).and_then(|v| v.parse().ok()).unwrap_or(d);
    let t_end = num("t", 6.0);
    let avg = num("avg", 2.0);
    const U: f32 = 0.1;
    let l = 16usize;
    let (body, re, depth, bands): (Body3D, f64, usize, Vec<f64>) = match case.as_str() {
        "ar4" | "ar4b" | "ar4n" | "ar4r2" => {
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
                _ => (0, vec![]),
            };
            (b, 1e6, depth, bands)
        }
        "sphere" => (Body3D::new(Shape3D::Sphere, l, 0.0), 1000.0, 0, vec![]),
        _ => {
            eprintln!("cas : ar4, ar4b, ar4n, ar4r2, sphere");
            std::process::exit(2)
        }
    };
    let mesh = body.mesh.is_some();
    let nu = U as f64 * l as f64 / re;
    let sa = re >= lbm::tunnel::RE_SA_AUTO;
    let t0 = Instant::now();
    let mut sim = Gpu3D::with_refinement_edges(&body, nu, U, true, sa, depth, false, false, &bands).expect("GPU");
    let prep = t0.elapsed().as_secs_f64();
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
