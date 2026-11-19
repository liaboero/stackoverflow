//! Validation automobile : corps d'Ahmed (lunette 25° ou 35°) sur sol mobile.
//!
//! Géométrie : `bodies3d::ahmed_tris` (1044 × 389 × 288 mm, avant en congés
//! de 100 mm, garde au sol 50 mm, quatre pieds de 30 mm). Référence de
//! longueur : la longueur de la maquette ; de surface : 389 × 288 mm.
//! Mesures :
//! - Cx par bilan de quantité de mouvement sur des boîtes autour de la
//!   maquette (moyenne des flux instantanés), face inférieure à mi-hauteur
//!   de la garde au sol : elle coupe les pieds, dont la partie haute reste
//!   dans la boîte (retirée ensuite, Ahmed corrige sa traînée de celle des
//!   pieds) ;
//! - Cx intégré sur la surface (`Gpu3D::surface_forces`), séparé comme chez
//!   Ahmed en pression sur l'avant (x < 100 mm, C_K), la lunette (C_S), le
//!   culot (C_B) et le reste (C_R = Cx − C_K − C_S − C_B, frottement surtout) ;
//! - profils de vitesse U, W dans le plan de symétrie, aux points mesurés
//!   par Lienhart et Schenk (fichiers `validation/ahmed/lienhart_*_y0.csv`).
//!
//! `cargo run --release --example ahmed -- --phi 25 --l 64` (léger) ;
//! options : --marge (veine, 0,5), --t (durée, t* = L/U), --moy (fenêtre de
//! moyenne), --re (2,78e6), --pieds 0|1, --memoire (Go, plafond), --sortie.
use std::collections::HashMap;
use std::io::Write;
use std::time::Instant;

use lbm::bodies3d::{ahmed, ahmed_tris, Body3D};
use lbm::cpu::Ramp;
use lbm::gpu3d::Gpu3D;
use lbm::stl::{MeshSetup, Unit};

struct Args(Vec<String>);

impl Args {
    fn get<T: std::str::FromStr>(&self, k: &str, d: T) -> T {
        let key = format!("--{k}");
        self.0.iter().position(|a| *a == key).and_then(|i| self.0.get(i + 1)).and_then(|v| v.parse().ok()).unwrap_or(d)
    }
}

/// Flux sortant ∮ (ρ u u·n + p n) dA sur la boîte [lo, hi] (nœuds), nœuds
/// solides exclus (la face coupe les pieds) ; force sur la maquette = −flux.
fn box_flux(mac: &[[f32; 4]], solid: &dyn Fn(usize, usize, usize) -> bool, d: [usize; 3], lo: [usize; 3], hi: [usize; 3]) -> [f64; 3] {
    let mut flux = [0.0f64; 3];
    for (ax, sgn, pos) in [(0, -1.0, lo[0]), (0, 1.0, hi[0]), (1, -1.0, lo[1]), (1, 1.0, hi[1]), (2, -1.0, lo[2]), (2, 1.0, hi[2])] {
        let (a1, a2) = ((ax + 1) % 3, (ax + 2) % 3);
        for i in lo[a1]..=hi[a1] {
            for j in lo[a2]..=hi[a2] {
                let w = (if i == lo[a1] || i == hi[a1] { 0.5 } else { 1.0 }) * (if j == lo[a2] || j == hi[a2] { 0.5 } else { 1.0 });
                let mut p = [0usize; 3];
                p[ax] = pos;
                p[a1] = i;
                p[a2] = j;
                if solid(p[0], p[1], p[2]) {
                    continue;
                }
                let c = mac[p[0] + d[0] * (p[1] + d[1] * p[2])];
                let (u, rho) = ([c[0] as f64, c[1] as f64, c[2] as f64], c[3] as f64);
                let un = sgn * u[ax];
                for k in 0..3 {
                    flux[k] += w * (rho * u[k] * un + if k == ax { sgn * (rho - 1.0) / 3.0 } else { 0.0 });
                }
            }
        }
    }
    flux
}

/// Parts de la traînée (unités du réseau), sommées sur un relevé de surface.
#[derive(Default, Clone, Copy)]
struct Parts {
    /// Pression : avant, lunette, culot, autres faces.
    k: f64,
    s: f64,
    b: f64,
    autres: f64,
    /// Frottement de la caisse.
    frott: f64,
    /// Pieds : total, et partie au-dessus de la face basse des boîtes.
    pieds: f64,
    pieds_haut: f64,
}

fn main() {
    let a = Args(std::env::args().collect());
    let phi: f64 = a.get("phi", 25.0);
    let l: usize = a.get("l", 64);
    let marge: f64 = a.get("marge", 0.5);
    let t_end: f64 = a.get("t", 12.0);
    let moy: f64 = a.get("moy", 4.0);
    let re: f64 = a.get("re", 2.78e6);
    let pieds = a.get("pieds", 1u32) == 1;
    let mem_max: f64 = a.get("memoire", 6.5);
    let sortie: String = a.get("sortie", "target/ahmed".to_string());
    let u = 0.1f32;
    std::fs::create_dir_all(&sortie).unwrap();

    // Maquette : mm par maille, facettes d'au plus deux mailles (la surface
    // est relevée maille par maille dans chaque facette).
    let h = ahmed::LONGUEUR / l as f64;
    let ds = (2.0 * h).min(10.0);
    let depth = h;
    let tris = ahmed_tris(phi, pieds, ds, depth);
    lbm::stl::write_stl(std::path::Path::new(&format!("{sortie}/ahmed_{phi}.stl")), &tris).unwrap();
    // Garde : le bas de la boîte englobante (bas des pieds, à −50 − depth mm)
    // est posé à `ground` + ¼ de maille au-dessus du plancher (z = −½) ; la
    // caisse est alors exactement à 50 mm du plancher.
    let setup = MeshSetup {
        unit: Unit::Mm,
        ref_length: Some(ahmed::LONGUEUR),
        ref_area: Some(ahmed::LARGEUR * ahmed::HAUTEUR),
        margin: marge,
        ground: Some(-depth - 0.25 * h),
        ..Default::default()
    };
    let (body, info) = Body3D::from_mesh(&tris, l, &setup).unwrap_or_else(|e| panic!("maquette : {e}"));
    let n = body.nx * body.ny * body.nz;
    let mem = n as f64 * (lbm::gpu3d::bytes_per_cell(true) + 32) as f64 / 1e9;
    println!(
        "Ahmed {phi}° · L = {l} mailles (h = {h:.2} mm) · veine {}×{}×{} = {:.1} M mailles ≈ {mem:.1} Go · blocage {:.1} % · {} triangles{}",
        body.nx,
        body.ny,
        body.nz,
        n as f64 / 1e6,
        100.0 * info.blockage,
        tris.len(),
        if info.diag.watertight() { "" } else { " (maillage ouvert !)" }
    );
    if mem > mem_max {
        eprintln!("Mémoire estimée {mem:.1} Go > --memoire {mem_max} : réduire L ou la marge.");
        std::process::exit(1);
    }
    let nu = u as f64 * l as f64 / re;
    let mut sim = Gpu3D::with_turbulence(&body, nu, u, true, true).unwrap_or_else(|e| panic!("GPU : {e}"));
    println!("{} · τ₀ = {:.6} · liens de paroi {}", sim.adapter_name(), sim.tau0, sim.wall_count);
    let spt = l as f64 / u as f64;
    let ramp = Ramp { steps: (4.0 * spt) as u64, u };
    let q = 0.5 * (u as f64).powi(2) * body.area;
    let d = [sim.nx, sim.ny, sim.nz];

    // Repère de la maquette (mm) depuis le réseau : x depuis l'avant, y depuis
    // le plan de symétrie, z depuis le dessous de la caisse.
    let bb = body.bounds;
    let ymid = 0.5 * (bb.y0 + bb.y1);
    let z_caisse = bb.z0 + (ahmed::GARDE + depth) / h;
    let to_mm = |p: [f64; 3]| [(p[0] - bb.x0) * h, (p[1] - ymid) * h, (p[2] - z_caisse) * h];
    let xs = ahmed::LONGUEUR - ahmed::LUNETTE * phi.to_radians().cos();
    // Aire de chaque triangle placé, retrouvé par son centre.
    let key = |c: [f64; 3]| c.map(|v| v.to_bits());
    let mut aire: HashMap<[u64; 3], f64> = HashMap::new();
    for t in &body.mesh.as_ref().unwrap().tris {
        let c = [0, 1, 2].map(|k| (t[0][k] + t[1][k] + t[2][k]) / 3.0);
        let e1 = [0, 1, 2].map(|k| t[1][k] - t[0][k]);
        let e2 = [0, 1, 2].map(|k| t[2][k] - t[0][k]);
        let cr = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]];
        aire.insert(key(c), 0.5 * (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt());
    }
    // Boîtes : marges en longueurs de maquette ; face basse à mi-garde.
    let marges = [0.1, 0.2, 0.3];
    let z_cut = (z_caisse - 0.5 * ahmed::GARDE / h).round() as usize;
    let z_cut_mm = (z_cut as f64 - z_caisse) * h;
    let lf = l as f64;
    let boites: Vec<Option<([usize; 3], [usize; 3])>> = marges
        .iter()
        .map(|&m| {
            let lo = [bb.x0 - m * lf, bb.y0 - m * lf, z_cut as f64];
            let hi = [bb.x1 + m * lf, bb.y1 + m * lf, bb.z1 + m * lf];
            let ok = (0..3).all(|k| lo[k] >= (if k == 2 { 1.0 } else { 3.0 }) && hi[k] <= (d[k] - 4) as f64);
            ok.then(|| (lo.map(|v| v.floor() as usize), hi.map(|v| v.ceil() as usize)))
        })
        .collect();

    let t0 = Instant::now();
    let chunk = (0.25 * spt).round() as usize;
    let total = (t_end * spt).round() as usize;
    let debut = ((t_end - moy) * spt).round() as usize;
    let mut step = 0usize;
    let mut somme: Vec<[f32; 4]> = Vec::new();
    let mut nsnap = 0usize;
    let mut boite_acc = vec![[0.0f64; 3]; marges.len()];
    let mut parts = Parts::default();
    let mut nsurf = 0usize;
    let mut manque = 0.0;
    let mut next_print = spt;
    while step < total {
        if step == debut {
            let _ = sim.take_forces();
        }
        sim.run(chunk, &ramp);
        step += chunk;
        if step as f64 >= next_print {
            next_print += spt;
            let el = t0.elapsed().as_secs_f64();
            println!(
                "t* = {:.1} · {:.0} s · {:.0} M mailles/s",
                step as f64 / spt,
                el,
                step as f64 * n as f64 / el / 1e6
            );
        }
        // Relevés tous les ½ t* pendant la fenêtre de moyenne.
        if step > debut && (step - debut) % (2 * chunk) == 0 {
            let mac = sim.macros();
            if !mac[n / 2][3].is_finite() {
                eprintln!("Calcul divergé à t* = {:.1}", step as f64 / spt);
                std::process::exit(2);
            }
            if somme.is_empty() {
                somme = vec![[0.0; 4]; n];
            }
            for (s, m) in somme.iter_mut().zip(&mac) {
                for k in 0..4 {
                    s[k] += m[k];
                }
            }
            nsnap += 1;
            let solid = |x: usize, y: usize, z: usize| sim.is_solid(x, y, z);
            for (i, b) in boites.iter().enumerate() {
                if let Some((lo, hi)) = b {
                    let f = box_flux(&mac, &solid, d, *lo, *hi);
                    for k in 0..3 {
                        boite_acc[i][k] -= f[k];
                    }
                }
            }
            drop(mac);
            if let Some(sf) = sim.surface_forces() {
                manque += sf.missed;
                nsurf += 1;
                for ((c, f), (nrm, pw)) in sf.loads.iter().zip(&sf.walls) {
                    let ar = aire.get(&key(*c)).copied().unwrap_or(0.0);
                    let fp = -pw * nrm[0] * ar;
                    let fr = f[0] - fp;
                    let p = to_mm(*c);
                    if p[2] < -0.5 {
                        parts.pieds += f[0];
                        if p[2] >= z_cut_mm {
                            parts.pieds_haut += f[0];
                        }
                        continue;
                    }
                    parts.frott += fr;
                    if p[0] < ahmed::RAYON {
                        parts.k += fp;
                    } else if nrm[0] > 0.9 && p[0] > ahmed::LONGUEUR - 1.0 {
                        parts.b += fp;
                    } else if nrm[0] > 0.05 && nrm[2] > 0.05 && p[0] > xs {
                        parts.s += fp;
                    } else {
                        parts.autres += fp;
                    }
                }
            }
        }
    }
    let fex = sim.take_forces();
    let duree = t0.elapsed().as_secs_f64();
    let ns = nsurf.max(1) as f64;
    let c = |v: f64| v / ns / q;
    let (ck, cs, cb, cautres, cf) = (c(parts.k), c(parts.s), c(parts.b), c(parts.autres), c(parts.frott));
    let cx_surf = ck + cs + cb + cautres + cf;
    let (cp_pieds, cp_haut) = (c(parts.pieds), c(parts.pieds_haut));
    println!("\n== Ahmed {phi}°, L = {l}, Re = {re:.3e}, {t_end} t* (moyenne sur {moy} t*, {nsnap} relevés), {duree:.0} s");
    println!(
        "Surface (caisse seule) : Cx = {cx_surf:.4} · C_K avant {ck:.4} · C_S lunette {cs:.4} · C_B culot {cb:.4} · pression ailleurs {cautres:.4} · frottement {cf:.4} · C_R = Cx − C_K − C_S − C_B = {:.4} (non relevé {:.2} %)",
        cx_surf - ck - cs - cb,
        100.0 * manque / ns
    );
    println!("Pieds (surface) : Cx {cp_pieds:.4}, dont {cp_haut:.4} au-dessus de la face basse des boîtes (z = {z_cut_mm:.0} mm)");
    println!("Échange de quantité de mouvement (caisse et pieds) : Cx = {:.4}", fex.fx / q);
    let mut bx = Vec::new();
    for (i, b) in boites.iter().enumerate() {
        if let Some((lo, hi)) = b {
            let cx = boite_acc[i][0] / nsnap.max(1) as f64 / q;
            let vol: f64 = (0..3).map(|k| (hi[k] - lo[k]) as f64).product();
            bx.push((vol, cx - cp_haut));
            println!("Boîte ±{} L : Cx = {:.4} (pieds compris), {:.4} sans le haut des pieds", marges[i], cx, cx - cp_haut);
        }
    }
    // Extrapolation linéaire à volume nul (moindres carrés).
    let cx_box0 = if bx.len() >= 2 {
        let nb = bx.len() as f64;
        let (mv, mc) = (bx.iter().map(|b| b.0).sum::<f64>() / nb, bx.iter().map(|b| b.1).sum::<f64>() / nb);
        let sl = bx.iter().map(|b| (b.0 - mv) * (b.1 - mc)).sum::<f64>() / bx.iter().map(|b| (b.0 - mv).powi(2)).sum::<f64>();
        let v = mc - sl * mv;
        println!("Boîte extrapolée à volume nul, sans les pieds : Cx = {v:.4}");
        v
    } else {
        f64::NAN
    };

    // Profils dans le plan de symétrie (champ moyen).
    let inv = 1.0 / nsnap.max(1) as f32;
    let mean = |x: usize, y: usize, z: usize| -> Option<[f64; 4]> {
        if sim.is_solid(x, y, z) {
            return None;
        }
        let s = somme[x + d[0] * (y + d[1] * z)];
        Some(s.map(|v| (v * inv) as f64))
    };
    let interp = |p: [f64; 3]| -> Option<[f64; 4]> {
        let i = p.map(|v| v.floor() as i64);
        if (0..3).any(|k| i[k] < 0 || i[k] + 1 >= d[k] as i64) {
            return None;
        }
        let w = [0, 1, 2].map(|k| p[k] - i[k] as f64);
        let mut out = [0.0; 4];
        for cc in 0..8 {
            let o = [cc & 1, (cc >> 1) & 1, (cc >> 2) & 1];
            let v = mean((i[0] + o[0] as i64) as usize, (i[1] + o[1] as i64) as usize, (i[2] + o[2] as i64) as usize)?;
            let ww: f64 = (0..3).map(|k| if o[k] == 1 { w[k] } else { 1.0 - w[k] }).product();
            for k in 0..4 {
                out[k] += ww * v[k];
            }
        }
        Some(out)
    };
    let ub = 40.0;
    let fic = format!("{}/../validation/ahmed/lienhart_{}_y0.csv", env!("CARGO_MANIFEST_DIR"), phi.round() as i64);
    match std::fs::read_to_string(&fic) {
        Err(_) => println!("(pas de mesures pour {phi}° : {fic})"),
        Ok(txt) => {
            let path = format!("{sortie}/profils_{phi}_L{l}.csv");
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "releve,x_mm,z_mm,U_mes,W_mes,U_calc,W_calc").unwrap();
            // Par station : écart moyen |ΔU|/Ub et nombre de points.
            let mut st: Vec<(String, f64, f64, f64, usize)> = Vec::new();
            for line in txt.lines().filter(|s| !s.starts_with('#') && !s.starts_with("releve")) {
                let v: Vec<&str> = line.split(',').collect();
                let (x, z, um, wm): (f64, f64, f64, f64) = (v[1].parse().unwrap(), v[2].parse().unwrap(), v[3].parse().unwrap(), v[4].parse().unwrap());
                // Repère de Lienhart : x = 0 au culot, z = 0 sur le plancher (z = −½ du réseau).
                let p = [bb.x0 + (x + ahmed::LONGUEUR) / h, ymid, -0.5 + z / h];
                let Some(m) = interp(p) else { continue };
                let (uc, wc) = (ub * m[0] / u as f64, ub * m[2] / u as f64);
                writeln!(f, "{},{x},{z},{um},{wm},{uc:.3},{wc:.3}", v[0]).unwrap();
                let du = (uc - um).abs() / ub;
                match st.iter_mut().find(|s| s.0 == v[0] && s.1 == x) {
                    Some(s) => {
                        s.2 += du;
                        s.3 += (um < 0.0) as u8 as f64 - (uc < 0.0) as u8 as f64;
                        s.4 += 1;
                    }
                    None => st.push((v[0].to_string(), x, du, (um < 0.0) as u8 as f64 - (uc < 0.0) as u8 as f64, 1)),
                }
            }
            st.sort_by(|p, q| p.1.partial_cmp(&q.1).unwrap());
            println!("Profils U (plan de symétrie) : écart moyen |ΔU|/Ub par station (x depuis le culot, mm) ; ± : points en retour mesurés − calculés");
            let mut tot = (0.0, 0usize);
            let mut lunette = (0.0, 0usize);
            let mut sillage = (0.0, 0usize);
            for s in &st {
                print!("  {} x={:.0}: {:.3} ({} pts, ±{:.0})", s.0, s.1, s.2 / s.4 as f64, s.4, s.3);
                println!();
                tot.0 += s.2;
                tot.1 += s.4;
                if s.1 < 0.0 && s.1 > -260.0 {
                    lunette.0 += s.2;
                    lunette.1 += s.4;
                } else if s.1 > 0.0 {
                    sillage.0 += s.2;
                    sillage.1 += s.4;
                }
            }
            let mean = |p: (f64, usize)| p.0 / p.1.max(1) as f64;
            println!(
                "Écart moyen |ΔU|/Ub : tous {:.3} ({} pts) · au-dessus de la lunette {:.3} · sillage {:.3} → {path}",
                mean(tot),
                tot.1,
                mean(lunette),
                mean(sillage)
            );
        }
    }
    // Longueur de la bulle de culot dans le plan de symétrie : dernier point
    // en retour (U < 0) à mi-hauteur du culot.
    let zb = ahmed::GARDE + 0.5 * (ahmed::HAUTEUR - ahmed::LUNETTE * phi.to_radians().sin());
    let mut lrec = 0.0;
    let mut x = 0.0;
    while x < 1000.0 {
        if let Some(m) = interp([bb.x0 + (x + ahmed::LONGUEUR) / h, ymid, -0.5 + zb / h]) {
            if m[0] < 0.0 {
                lrec = x;
            }
        }
        x += 0.5 * h;
    }
    println!("Bulle de culot (U < 0 à z = {zb:.0} mm, y = 0) jusqu'à x = {lrec:.0} mm derrière le culot");
    let res = format!("{sortie}/resultats.csv");
    let neuf = !std::path::Path::new(&res).exists();
    let mut f = std::fs::OpenOptions::new().create(true).append(true).open(&res).unwrap();
    if neuf {
        writeln!(f, "phi,L,marge,t,moy,mailles_M,blocage,cx_surface,c_k,c_s,c_b,c_r,cx_boite0,cx_echange,cx_pieds,bulle_mm,duree_s").unwrap();
    }
    writeln!(
        f,
        "{phi},{l},{marge},{t_end},{moy},{:.1},{:.4},{cx_surf:.4},{ck:.4},{cs:.4},{cb:.4},{:.4},{cx_box0:.4},{:.4},{cp_pieds:.4},{lrec:.0},{duree:.0}",
        n as f64 / 1e6,
        info.blockage,
        cx_surf - ck - cs - cb,
        fex.fx / q
    )
    .unwrap();
}
