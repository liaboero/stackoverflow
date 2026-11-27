//! Validation du décollement d'arête vive : marche descendante de Driver et
//! Seegmiller (1985).
//!
//! Canal plan (extrudé, envergure périodique) : hauteur 8 h en amont de la
//! marche, 9 h en aval (rapport d'expansion 1,125), Re_h = U_ref h / ν ≈ 36 000.
//! Les deux parois du canal sont visqueuses (loi de paroi) ; la couche limite
//! amont se développe depuis l'entrée (vitesse uniforme) sur `--amont` hauteurs
//! de marche, pour approcher δ ≈ 1,5 h à x = −4 h.
//!
//! Mesures (moyennes en temps et en envergure) :
//! - couche limite amont à x = −4 h : U_ref (cœur), δ₉₉, δ*, θ, Re_θ ;
//! - longueur de recollement x_r/h : changement de signe de u au premier rang
//!   de nœuds au-dessus de la paroi basse, et du frottement (échange de
//!   quantité de mouvement) ; mesure 6,26 ± 0,10 ;
//! - Cf(x) et Cp(x) sur la paroi basse en aval, Cp rapporté à x = −4 h ; Cf
//!   par l'échange de quantité de mouvement (moyenne) et par la loi de paroi
//!   (2 u_τ², état final, colonne `cf_loi`).
//!
//! `cargo run --release --example marche -- --n 16` ; options : --n (mailles
//! par hauteur de marche), --amont (longueur amont, h), --aval (h), --ny
//! (envergure, mailles), --t (durée, h/U), --moy (fenêtre de moyenne),
//! --re (Re_h), --essai (bits de `set_essai`), --q (position de la paroi basse
//! dans la maille, ½ par défaut), --fp16 (0 : simple précision), --sortie ;
//! `SOUFFLERIE_TBLE` (eq, gp, gpc) : loi de paroi hors équilibre.
//! Profils aval (u, ν_t, contrainte, −u'w' résolue) : `marche_profils_n*.csv`.
use std::io::Write;
use std::sync::Arc;

use lbm::bodies3d::{Body3D, Bounds};
use lbm::cpu::Ramp;
use lbm::gpu3d::Gpu3D;
use lbm::stl::{MeshBody, MeshSetup, Tri, Unit};

struct Args(Vec<String>);

impl Args {
    fn get<T: std::str::FromStr>(&self, k: &str, d: T) -> T {
        let key = format!("--{k}");
        self.0.iter().position(|a| *a == key).and_then(|i| self.0.get(i + 1)).and_then(|v| v.parse().ok()).unwrap_or(d)
    }
}

/// Prisme extrudé en y de y0 à y1, de section le polygone (x, z) convexe par
/// morceaux `prof` (sens trigonométrique vu de +y) ; `caps` : triangles (x, z)
/// des deux faces d'extrémité.
fn prism(prof: &[[f64; 2]], caps: &[[[f64; 2]; 3]], y0: f64, y1: f64, tris: &mut Vec<Tri>) {
    let n = prof.len();
    for k in 0..n {
        let (p, q) = (prof[k], prof[(k + 1) % n]);
        let (a0, b0, a1, b1) = ([p[0], y0, p[1]], [q[0], y0, q[1]], [p[0], y1, p[1]], [q[0], y1, q[1]]);
        tris.push([a0, a1, b1]);
        tris.push([a0, b1, b0]);
    }
    for t in caps {
        tris.push([[t[0][0], y0, t[0][1]], [t[1][0], y0, t[1][1]], [t[2][0], y0, t[2][1]]]);
        tris.push([[t[0][0], y1, t[0][1]], [t[2][0], y1, t[2][1]], [t[1][0], y1, t[1][1]]]);
    }
}

fn main() {
    let a = Args(std::env::args().collect());
    let n: usize = a.get("n", 16);
    let amont: f64 = a.get("amont", 80.0);
    let aval: f64 = a.get("aval", 40.0);
    let ny: usize = a.get("ny", 8);
    let t_end: f64 = a.get("t", 240.0);
    let moy: f64 = a.get("moy", 60.0);
    let re: f64 = a.get("re", 36000.0);
    let essai: u32 = a.get("essai", 0);
    let fp16 = a.get("fp16", 1u32) != 0;
    let sortie: String = a.get("sortie", "target/marche".to_string());
    std::fs::create_dir_all(&sortie).unwrap();
    let hf = n as f64;
    let u = 0.1f32;

    // Veine (mailles) : parois pleines en haut et en bas, plus épaisses que les
    // bandes absorbantes de la veine (nz/8) ; bords z périodiques, sans effet
    // dans les parois.
    let nx = ((amont + aval) * hf).round() as usize;
    let tw = (1.5 * hf + 3.0).ceil();
    let nz = (9.0 * hf + 2.0 * tw).round() as usize;
    assert!(tw > (nz as f64 / 8.0).max(2.0) + 1.0);
    // Parois à mi-chemin entre deux nœuds (q = ½) : paroi basse aval en z_lo,
    // dessus de la marche en z_lo + h, paroi haute en z_lo + 9 h ; arête de la
    // marche en x_s (entre deux nœuds aussi).
    let qz: f64 = a.get("q", 0.5);
    let z_lo = tw + qz;
    let z_st = z_lo + hf;
    let z_hi = z_lo + 9.0 * hf;
    let x_s = (amont * hf).round() + 0.5;
    let (xa, xb) = (-20.0, nx as f64 + 20.0);
    let (zb, zt) = (-5.0, nz as f64 + 5.0);
    let mut tris: Vec<Tri> = Vec::new();
    let (y0, y1) = (-20.0, ny as f64 + 20.0);
    // Paroi basse en L : P0..P6.
    let p = [[xa, zb], [x_s, zb], [xb, zb], [xb, z_lo], [x_s, z_lo], [x_s, z_st], [xa, z_st]];
    let caps = [[p[1], p[2], p[3]], [p[1], p[3], p[4]], [p[0], p[1], p[4]], [p[0], p[4], p[5]], [p[0], p[5], p[6]]];
    prism(&p, &caps, y0, y1, &mut tris);
    let q = [[xa, z_hi], [xb, z_hi], [xb, zt], [xa, zt]];
    prism(&q, &[[q[0], q[1], q[2]], [q[0], q[2], q[3]]], y0, y1, &mut tris);

    let setup = MeshSetup { unit: Unit::parse("m").unwrap(), margin: 0.2, ..Default::default() };
    let (mut body, _) = Body3D::from_mesh(&tris, 16, &setup).unwrap();
    let mesh = MeshBody::new(tris, false);
    body.nx = nx;
    body.ny = ny;
    body.nz = nz;
    body.l = n;
    body.bounds = Bounds { x0: 0.0, x1: nx as f64 - 1.0, y0: 0.0, y1: ny as f64 - 1.0, z0: 0.0, z1: nz as f64 - 1.0 };
    body.center = [x_s, 0.5 * ny as f64, 0.5 * (z_lo + z_hi)];
    body.area = hf * ny as f64;
    body.mesh = Some(Arc::new(mesh));
    body.features = None;

    let nu = u as f64 * hf / re;
    let cells = nx * ny * nz;
    println!(
        "Marche de Driver et Seegmiller · h = {n} mailles · veine {nx}×{ny}×{nz} = {:.2} M mailles · amont {amont} h, aval {aval} h · Re_h = {re:.0}",
        cells as f64 / 1e6
    );
    let mut sim = Gpu3D::with_turbulence(&body, nu, u, fp16, lbm::gpu3d::Turbulence::rans_from_env()).unwrap_or_else(|e| panic!("GPU : {e}"));
    sim.periodic_span = true;
    if essai != 0 {
        sim.set_essai(essai, 0.0);
    }
    // Loi de paroi hors équilibre (`SOUFFLERIE_TBLE` = eq, gp, gpc).
    if let Some(m) = lbm::tble::ParoiTble::from_env() {
        sim.set_paroi_tble(Some(m));
        println!("Loi de paroi TBLE : {m:?}");
    }
    println!("{} · τ₀ = {:.6} · liens de paroi {}", sim.adapter_name(), sim.tau0, sim.wall_count);
    let spt = hf / u as f64;
    let ramp = Ramp { steps: (10.0 * spt) as u64, u };
    let id = |x: usize, y: usize, z: usize| x + nx * (y + ny * z);

    let t0 = std::time::Instant::now();
    let chunk = (2.0 * spt).round() as usize;
    let total = (t_end * spt).round() as usize;
    let debut = ((t_end - moy) * spt).round() as usize;
    // Moyennes d'envergure du plan (x, z) : u, w, ρ.
    let mut somme = vec![[0.0f64; 4]; nx * nz];
    let mut nsnap = 0usize;
    // Frottement de la paroi basse par colonne x.
    let mut fric = vec![0.0f64; nx];
    let mut step = 0usize;
    let mut next_print = 20.0 * spt;
    while step < total {
        if step <= debut && step + chunk > debut {
            let _ = sim.take_wall_forces();
        }
        sim.run(chunk, &ramp);
        step += chunk;
        if step as f64 >= next_print {
            next_print += 20.0 * spt;
            let el = t0.elapsed().as_secs_f64();
            println!("t = {:.0} h/U · {:.0} s · {:.0} M mailles/s", step as f64 / spt, el, step as f64 * cells as f64 / el / 1e6);
        }
        if step > debut {
            let mac = sim.macros();
            if !mac[id(nx / 2, 0, nz / 2)][3].is_finite() {
                eprintln!("Calcul divergé à t = {:.0} h/U", step as f64 / spt);
                std::process::exit(2);
            }
            for z in 0..nz {
                for x in 0..nx {
                    let s = &mut somme[x + nx * z];
                    for y in 0..ny {
                        let m = mac[id(x, y, z)];
                        s[0] += m[0] as f64;
                        s[1] += m[2] as f64;
                        s[2] += m[3] as f64;
                        s[3] += m[0] as f64 * m[2] as f64;
                    }
                }
            }
            nsnap += 1;
        }
    }
    let wf = sim.take_wall_forces();
    // Mode pair-impair en x des premiers rangs (paroi à q ≈ ½ sous loi de
    // paroi) : |u − moyenne des deux voisines|/U, état final.
    {
        let mac = sim.macros();
        let z1 = z_lo.ceil() as usize;
        let zst = z_st.ceil() as usize;
        let oe = |x: usize, z: usize| {
            let g = |x: usize| mac[id(x, 0, z)][0] as f64;
            ((g(x) - 0.5 * (g(x - 1) + g(x + 1))) / u as f64).abs()
        };
        let xs_i = x_s.ceil() as usize;
        let pts: Vec<String> = [-10.0, 15.0]
            .iter()
            .map(|&xh| {
                let x = (x_s + xh * hf).round() as usize;
                let zb = if x < xs_i { zst } else { z1 };
                format!("{xh} h : {:.4} / {:.4}", oe(x, zb), oe(x, zb + 1))
            })
            .collect();
        println!("Mode pair-impair (rangs 1 / 2) : {}", pts.join(" · "));
    }
    // u_τ de la loi de paroi (état final), paroi basse aval.
    let mut cfl = vec![0.0f64; nx];
    for ((pp, _), st) in sim.finest_wall_cells().iter().zip(sim.finest_wall_state()) {
        if pp[2] < z_lo + 1.0 && pp[0] > x_s {
            let i = pp[0].round() as usize;
            if i < nx {
                cfl[i] += 2.0 * (st[0] as f64).powi(2) * st[1].signum() as f64;
            }
        }
    }
    for (pp, f) in wf.pos.iter().zip(&wf.force) {
        // Mailles de la paroi basse aval (premier rang au-dessus de z_lo).
        if pp[2] < z_lo + 1.0 && pp[0] > x_s {
            let i = pp[0].round() as usize;
            if i < nx {
                fric[i] += f[0];
            }
        }
    }
    let duree = t0.elapsed().as_secs_f64();
    let norm = 1.0 / (nsnap.max(1) * ny) as f64;
    let mean = |x: usize, z: usize| somme[x + nx * z].map(|v| v * norm);
    let solid = |x: usize, z: usize| sim.is_solid(x, 0, z);

    // Couche limite amont à x = −4 h (paroi basse : dessus de la marche).
    let xr = (x_s - 4.0 * hf).round() as usize;
    let zc = (0.5 * (z_st + z_hi)).round() as usize;
    let ucore = (z_st.ceil() as usize..zc).map(|z| mean(xr, z)[0]).fold(0.0f64, f64::max);
    let (mut d1, mut th, mut d99) = (0.0, 0.0, 0.0);
    {
        // u(0) extrapolé (paroi glissante de la loi de paroi) : nœuds seuls.
        let mut prev = (0.0f64, 0.0f64);
        for z in z_st.ceil() as usize..zc {
            let yv = z as f64 - z_st;
            let r = mean(xr, z)[0] / ucore;
            if d99 == 0.0 && r >= 0.99 {
                d99 = yv;
            }
            let dy = yv - prev.0;
            let r0 = if prev.0 == 0.0 { r } else { prev.1 };
            d1 += 0.5 * dy * ((1.0 - r0) + (1.0 - r));
            th += 0.5 * dy * (r0 * (1.0 - r0) + r * (1.0 - r));
            prev = (yv, r);
        }
    }
    let ret = ucore * th / nu;
    let re_h = ucore * hf / nu;
    let q = 0.5 * ucore * ucore;
    // Pression de référence : premier nœud au-dessus du dessus de la marche
    // à x = −4 h (station de référence de Driver et Seegmiller).
    let zw = z_st.ceil() as usize;
    let p_ref = (mean(xr, zw)[2] - 1.0) / 3.0;
    println!("\n== Marche, h = {n}, {t_end} h/U (moyenne {moy} h/U, {nsnap} relevés), {duree:.0} s");
    println!(
        "Amont x = −4 h : U_ref/U_entrée = {:.3}, Re_h = {re_h:.0}, δ₉₉ = {:.2} h, δ* = {:.3} h, θ = {:.3} h, H = {:.2}, Re_θ = {ret:.0} (Driver : δ ≈ 1,5 h, Re_θ ≈ 5000)",
        ucore / u as f64,
        d99 / hf,
        d1 / hf,
        th / hf,
        d1 / th.max(1e-12)
    );

    // Profil le long de la paroi basse aval : premier rang de nœuds.
    let z1 = z_lo.ceil() as usize;
    let xs_i = x_s.ceil() as usize;
    let xh = |x: usize| (x as f64 - x_s) / hf;
    let path = format!("{sortie}/marche_n{n}.csv");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "x_h,u1_uref,cf,cf_loi,cp").unwrap();
    let mut x_r_u = f64::NAN;
    let mut x_r_cf = f64::NAN;
    let mut u_prev = f64::NAN;
    let mut cf_prev = f64::NAN;
    let mut cf_min = (0.0f64, 0.0f64);
    let mut cp_min = (f64::MAX, 0.0f64);
    let mut cp_max = (f64::MIN, 0.0f64);
    let xe = nx - (0.15 * nx as f64) as usize - 1;
    for x in xs_i..xe {
        if solid(x, z1) {
            continue;
        }
        let m = mean(x, z1);
        let u1 = m[0] / ucore;
        let cf = fric[x] / ny as f64 / q;
        let cp = ((m[2] - 1.0) / 3.0 - p_ref) / q;
        let cfloi = cfl[x] / ny as f64 / (ucore * ucore);
        writeln!(f, "{:.4},{u1:.5},{cf:.6},{cfloi:.6},{cp:.5}", xh(x)).unwrap();
        // Recollement : dernier passage de u < 0 à u > 0 dans les 20 h aval.
        if xh(x) < 20.0 {
            if u_prev < 0.0 && u1 >= 0.0 {
                x_r_u = xh(x) - 1.0 / hf * u1 / (u1 - u_prev);
            }
            if cf_prev < 0.0 && cf >= 0.0 {
                x_r_cf = xh(x) - 1.0 / hf * cf / (cf - cf_prev);
            }
            if cf < cf_min.0 {
                cf_min = (cf, xh(x));
            }
            if cp < cp_min.0 {
                cp_min = (cp, xh(x));
            }
            if cp > cp_max.0 {
                cp_max = (cp, xh(x));
            }
        }
        u_prev = u1;
        cf_prev = cf;
    }
    let cf_at = |xv: f64| {
        let x = (x_s + xv * hf).round() as usize;
        fric[x] / ny as f64 / q
    };
    println!("Recollement x_r/h : {x_r_u:.2} (u au premier nœud), {x_r_cf:.2} (frottement) · mesure 6,26 ± 0,10");
    println!(
        "Cf min {:.5} à x/h = {:.2} · Cf à 10 h {:.5}, 20 h {:.5} · Cp min {:.3} à {:.2}, Cp max {:.3} à {:.2} → {path}",
        cf_min.0,
        cf_min.1,
        cf_at(10.0),
        cf_at(20.0),
        cp_min.0,
        cp_min.1,
        cp_max.0,
        cp_max.1
    );
    // Profils en aval (x = 10, 15, 20 h) : u moyen, ν_t de SA (état final,
    // ν̃ f_v1, sans Smagorinsky), contrainte (ν + ν_t) du/dy rapportée à q.
    let nt = sim.nu_tilde();
    if !nt.is_empty() {
        let pp = format!("{sortie}/marche_profils_n{n}.csv");
        let mut g = std::fs::File::create(&pp).unwrap();
        writeln!(g, "x_h,y_h,u_uref,nut_nu,tau_q,uw_res_q,cf_echange").unwrap();
        for xv in [10.0, 15.0, 20.0] {
            let x = (x_s + xv * hf).round() as usize;
            let nut = |z: usize| {
                let s: f64 = (0..ny).map(|y| nt[id(x, y, z)] as f64).sum::<f64>() / ny as f64;
                let chi = s / nu;
                s * chi.powi(3) / (chi.powi(3) + 7.1f64.powi(3))
            };
            for z in z1..(z_lo + 3.0 * hf) as usize {
                let yv = z as f64 - z_lo;
                let dudy = if z == z1 { mean(x, z + 1)[0] - mean(x, z)[0] } else { 0.5 * (mean(x, z + 1)[0] - mean(x, z - 1)[0]) };
                let tau = (nu + nut(z)) * dudy / q;
                let m = mean(x, z);
                let uw = -(m[3] - m[0] * m[1]) / q;
                writeln!(g, "{xv},{:.4},{:.5},{:.1},{:.6},{:.6},{:.6}", yv / hf, m[0] / ucore, nut(z) / nu, tau, uw, fric[x] / ny as f64 / q).unwrap();
            }
        }
    }
    let res = format!("{sortie}/resultats.csv");
    let neuf = !std::path::Path::new(&res).exists();
    let mut g = std::fs::OpenOptions::new().create(true).append(true).open(&res).unwrap();
    if neuf {
        writeln!(g, "n,amont,ny,t,moy,essai,rans,re_h,d99_h,theta_h,re_theta,xr_u,xr_cf,cf_min,cp_min,cp_max,duree_s").unwrap();
    }
    writeln!(
        g,
        "{n},{amont},{ny},{t_end},{moy},{essai},{},{re_h:.0},{:.3},{:.4},{ret:.0},{x_r_u:.3},{x_r_cf:.3},{:.5},{:.4},{:.4},{duree:.0}",
        std::env::var("SOUFFLERIE_RANS").unwrap_or_else(|_| "sa".into()),
        d99 / hf,
        th / hf,
        cf_min.0,
        cp_min.0,
        cp_max.0
    )
    .unwrap();
}
