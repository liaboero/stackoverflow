//! Diagnostic 3D : couche limite turbulente d'une plaque plane, contre la
//! théorie (loi en puissance 1/7, profil de Coles).
//!
//! plaque3d <L> <raff> <Re_L> <t> [essai] [demi-précision 1/0] [ny] [cale] [csv ou -] [nez]
//!
//! Plaque épaisse (épaisseur L/16) à bords d'attaque et de fuite
//! demi-elliptiques 4:1, alignée sur l'écoulement, d'envergure périodique (elle
//! traverse les bords y de la veine, sans bandes absorbantes latérales) :
//! (allongement `nez` des ellipses, 4 par défaut)
//! écoulement plan. L mailles par longueur de plaque, `raff` grilles fines,
//! Re_L = U L / ν, `t` longueurs de plaque parcourues. `essai` : variantes de
//! diagnostic de `Gpu3D::set_essai`. `cale` : distance (mailles de la veine)
//! du premier nœud au-dessus de la face supérieure. `csv` : profils u(y).
//!
//! Profils u(y) relevés sur l'extrados (moyenne sur l'envergure) à plusieurs
//! abscisses. Sous le premier nœud, profil de Spalding avec le u_τ de la loi
//! de paroi ; au-dessus, les nœuds. δ*, θ, H ; Cf par la loi de paroi (u_τ
//! tiré de la vitesse à WM_H) et par le bilan de quantité de mouvement
//! (dθ/dx + (2 + H) θ/U_e dU_e/dx = Cf/2).
use lbm::bodies3d::Body3D;
use lbm::gpu3d::{Gpu3D, WM_H};
use lbm::stl::{MeshBody, MeshSetup, Tri, Unit};
use std::fmt::Write as _;
use std::sync::Arc;

const KAPPA: f64 = 0.41;
const B_LOG: f64 = 5.0;
/// Paramètre de sillage de Coles aux Re_θ élevés.
const PI_COLES: f64 = 0.55;

/// u⁺ de la loi de Spalding en y⁺ (Newton).
fn spalding_uplus(yp: f64) -> f64 {
    let e_kb = (-KAPPA * B_LOG).exp();
    let mut up = if yp < 10.0 { yp } else { yp.ln() / KAPPA + B_LOG };
    for _ in 0..50 {
        let ku = KAPPA * up;
        let e = ku.exp();
        let f = up + e_kb * (e - 1.0 - ku - ku * ku / 2.0 - ku.powi(3) / 6.0) - yp;
        let df = 1.0 + e_kb * KAPPA * (e - 1.0 - ku - ku * ku / 2.0);
        up -= f / df;
        if f.abs() < 1e-10 * yp.max(1.0) {
            break;
        }
    }
    up.max(0.0)
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let g = |i: usize, d: f64| a.get(i).and_then(|v| v.parse().ok()).unwrap_or(d);
    let l = g(1, 64.0) as usize;
    let raff = g(2, 0.0) as usize;
    let re = g(3, 1e6);
    let t_end = g(4, 12.0);
    let essai = g(5, 0.0) as u32;
    let fp16 = g(6, 1.0) > 0.5;
    let ny = g(7, if raff > 0 { 20.0 } else { 8.0 }) as usize;
    let cale = g(8, 0.3);
    let csv = a.get(9).cloned().filter(|c| c != "-");
    // Allongement des bords d'attaque et de fuite demi-elliptiques.
    let nez = g(10, 4.0);
    let lf = l as f64;
    // Plafond de viscosité turbulente des essais 1024 et 2048 : valeur de
    // Clauser 0,0168 U δ* de la plaque théorique à mi-longueur (unités de la veine).
    let nu_cap = 0.0168 * 0.1 * 0.046 * 0.5 * (0.5 * re).powf(-0.2) * lf * std::env::var("PLAQUE_CAP").ok().and_then(|v| v.parse::<f64>().ok()).unwrap_or(1.0);

    // Veine : 1,5 L en amont, 2,5 L en aval de la plaque (zone éponge à
    // 0,85 nx, au-delà des grilles fines), 3 L de haut.
    let (nx, nz) = ((5.0 * lf).round() as usize, (3.0 * lf).round() as usize);
    let x_le = 1.5 * lf + 0.25;
    // Épaisseur ≈ L/16, ajustée pour que le premier nœud soit à `cale` de
    // la paroi sur les deux faces : 2b = n + 1 − 2 cale. Sinon les deux faces
    // n'ont pas la même couche limite numérique et la plaque porte.
    let b = 0.5 * (((lf / 16.0).round() - 1.0).max(1.0) + 1.0 - 2.0 * cale);
    let a_ell = nez * b;
    let z_top = (0.5 * nz as f64).floor() + 1.0 - cale;
    let zc = z_top - b;
    // Contour (x, z) : extrados du bord d'attaque au bord de fuite, puis intrados.
    let ne = 24;
    let mut prof: Vec<[f64; 2]> = Vec::new();
    for k in 0..=ne {
        let th = std::f64::consts::PI * (1.0 - k as f64 / ne as f64) / 2.0;
        prof.push([x_le + a_ell * (1.0 - th.sin()), zc + b * th.cos()]);
    }
    for k in 1..=ne {
        let ph = std::f64::consts::FRAC_PI_2 * k as f64 / ne as f64;
        prof.push([x_le + lf - a_ell + a_ell * ph.sin(), zc + b * ph.cos()]);
    }
    let top: Vec<[f64; 2]> = prof.clone();
    // Intrados : symétrique, du bord de fuite au bord d'attaque.
    for p in top.iter().rev().skip(1).take(top.len() - 2) {
        prof.push([p[0], 2.0 * zc - p[1]]);
    }
    let (y0, y1) = (-20.0, ny as f64 + 20.0);
    let mut tris: Vec<Tri> = Vec::new();
    let n = prof.len();
    for k in 0..n {
        let (p, q) = (prof[k], prof[(k + 1) % n]);
        let (a0, b0, a1, b1) = ([p[0], y0, p[1]], [q[0], y0, q[1]], [p[0], y1, p[1]], [q[0], y1, q[1]]);
        tris.push([a0, b0, b1]);
        tris.push([a0, b1, a1]);
    }
    let cx = prof.iter().map(|p| p[0]).sum::<f64>() / n as f64;
    for k in 0..n {
        let (p, q) = (prof[k], prof[(k + 1) % n]);
        tris.push([[cx, y0, zc], [q[0], y0, q[1]], [p[0], y0, p[1]]]);
        tris.push([[cx, y1, zc], [p[0], y1, p[1]], [q[0], y1, q[1]]]);
    }
    let setup = MeshSetup { unit: Unit::parse("m").unwrap(), margin: 0.2, ..Default::default() };
    let (mut body, _) = Body3D::from_mesh(&tris, l, &setup).unwrap();
    let mesh = MeshBody::new(tris, false);
    body.nx = nx;
    body.ny = ny;
    body.nz = nz;
    body.bounds = lbm::bodies3d::Bounds { x0: mesh.lo[0], x1: mesh.hi[0], y0: 0.0, y1: ny as f64 - 1.0, z0: mesh.lo[2], z1: mesh.hi[2] };
    body.center = [x_le + 0.5 * lf, 0.5 * ny as f64, zc];
    body.area = lf * ny as f64;
    body.mesh = Some(Arc::new(mesh));

    let u = 0.1f32;
    let nu = u as f64 * lf / re;
    let spt = lf / u as f64;
    let ramp = lbm::cpu::Ramp { steps: (2.0 * spt) as u64, u };
    // Profils u(y) de l'extrados (moyenne sur la moitié centrale de
    // l'envergure de la grille la plus fine) : (distance à la paroi en mailles
    // fines, u) par abscisse ; facteur d'échelle de la grille ; traînée
    // (échange de quantité de mouvement).
    let stations: Vec<f64> = (0..15).map(|k| 0.15 + 0.05 * k as f64).collect();
    let run = |essai: u32| -> Option<(Vec<Vec<(f64, f64)>>, f64, [usize; 3], f64, f64, f64)> {
        let t0 = std::time::Instant::now();
        let mut sim = Gpu3D::with_refinement(&body, nu, u, fp16, true, raff).unwrap();
        sim.periodic_span = true;
        sim.set_essai(essai, nu_cap);
        sim.run(((t_end - 2.0) * spt) as usize, &ramp);
        let _ = sim.take_forces();
        sim.run((2.0 * spt) as usize, &ramp);
        let f = sim.take_forces();
        let fx = f.fx;
        // Frottement transmis par la paroi (échange de quantité de mouvement)
        // sur la partie plane de l'extrados, 0,3 à 0,7 L, rapporté à ½ρU².
        sim.run((0.5 * spt) as usize, &ramp);
        let wf = sim.take_wall_forces();
        let fx_top: f64 = wf
            .pos
            .iter()
            .zip(&wf.force)
            .filter(|(p, _)| p[2] > zc && p[0] > x_le + 0.3 * lf && p[0] < x_le + 0.7 * lf)
            .map(|(_, f)| f[0])
            .sum();
        let cf_ex = fx_top / (0.5 * (u as f64).powi(2) * 0.4 * lf * ny as f64);
        if std::env::var("PLAQUE_NEZ").is_ok() {
            // Effort (x, z) des mailles pariétales de l'extrados par maille de
            // la veine le long de x, rapporté à ½ρU² × envergure (mailles).
            let mut bins = vec![[0.0f64; 2]; 16];
            for (p, f) in wf.pos.iter().zip(&wf.force) {
                let k = (p[0] - x_le + 1.0).floor();
                if p[2] > zc && k >= 0.0 && (k as usize) < bins.len() {
                    bins[k as usize][0] += f[0];
                    bins[k as usize][1] += f[2];
                }
            }
            let q = 0.5 * (u as f64).powi(2) * ny as f64;
            let v: Vec<String> = bins.iter().enumerate().map(|(k, b)| format!("{}:{:+.3}/{:+.3}", k as i64 - 1, b[0] / q, b[1] / q)).collect();
            println!("essai {essai} nez, Fx/Fz par maille : {}", v.join(" "));
        }
        if std::env::var("PLAQUE_COLONNES").is_ok() {
            let wc = sim.finest_wall_cells();
            let edges = wc.iter().filter(|c| c.1 < 0.5).count();
            let flat = wc.iter().filter(|c| c.1 < 0.5 && c.0[0] > x_le + 0.2 * lf && c.0[0] < x_le + 0.8 * lf).count();
            println!("essai {essai} : Cz {:.4}, {} mailles pariétales, {edges} d'arête dont {flat} sur la partie plane", f.fz / (0.5 * (u as f64).powi(2) * body.area), wc.len());
        }
        let (mac, d, org, s) = sim.finest_macros();
        if mac.iter().any(|m| !m[0].is_finite()) {
            return None;
        }
        if std::env::var("PLAQUE_COLONNES").is_ok() {
            // Colonnes u(z) entières (moyenne d'envergure) à quelques abscisses.
            for xl in [-0.3, 0.1, 0.5, 0.9, 1.5] {
                let i = ((x_le + xl * lf - org[0]) * s).round() as usize;
                let col: Vec<String> = (0..d[2])
                    .step_by((d[2] / 24).max(1))
                    .map(|k| {
                        let v: f64 = (0..d[1]).map(|j| mac[i + d[0] * (j + d[1] * k)][0] as f64).sum::<f64>() / d[1] as f64;
                        format!("{:.3}", v / u as f64)
                    })
                    .collect();
                println!("x/L {xl:+.1} : {}", col.join(" "));
            }
        }
        let (j0, j1) = (d[1] / 4, (3 * d[1] / 4).max(d[1] / 4 + 1));
        if std::env::var("PLAQUE_PERTE").is_ok() {
            // Épaisseur de perte de pression totale de la moitié supérieure,
            // ∫ (1 − p₀/p₀∞) dz jusqu'à 0,2 L au-dessus de l'axe (mailles de
            // la veine), le long de x : où la perte naît-elle ?
            let q0 = 0.5 * (u as f64).powi(2);
            let mut line = String::new();
            let zk = |k: usize| org[2] + k as f64 / s;
            for xi in 0..60 {
                let xl = -0.05 + 0.01 * xi as f64;
                let i = ((x_le + xl * lf - org[0]) * s).round();
                if i < 0.0 || i as usize >= d[0] {
                    continue;
                }
                let i = i as usize;
                let mut acc = 0.0;
                for k in 0..d[2] {
                    if zk(k) <= zc || zk(k) > zc + 0.2 * lf {
                        continue;
                    }
                    for j in j0..j1 {
                        let m = mac[i + d[0] * (j + d[1] * k)];
                        if m[3] == 0.0 {
                            continue;
                        }
                        let r = m[3] as f64;
                        let v2 = (m[0] as f64).powi(2) + (m[1] as f64).powi(2) + (m[2] as f64).powi(2);
                        let p0 = (r - 1.0) / 3.0 + 0.5 * r * v2;
                        acc += (1.0 - p0 / q0) / s;
                    }
                }
                let _ = write!(line, " {:+.2}:{:.3}", xl, acc / (j1 - j0) as f64);
            }
            println!("essai {essai} perte de p₀ (x/L : mailles) :{line}");
        }
        let ye = 0.2 * lf * s;
        let profs = stations
            .iter()
            .map(|&xl| {
                let i = (((x_le + xl * lf - org[0]) * s).round() as usize).min(d[0] - 1);
                (0..d[2])
                    .filter_map(|k| {
                        let dist = (org[2] + k as f64 / s - z_top) * s;
                        if dist <= 0.0 || dist > ye {
                            return None;
                        }
                        let acc: f64 = (j0..j1).map(|j| mac[i + d[0] * (j + d[1] * k)][0] as f64).sum();
                        Some((dist, acc / (j1 - j0) as f64))
                    })
                    .collect()
            })
            .collect();
        Some((profs, s, d, fx, cf_ex, t0.elapsed().as_secs_f64()))
    };
    let Some((profs, s, d, fx, cf_ex, el)) = run(essai) else {
        println!("divergence");
        return;
    };
    // Écoulement extérieur : la plaque épaisse a sa propre survitesse
    // potentielle, qui décroît en s'éloignant de la paroi. Droite ajustée sur
    // les relevés entre 0,1 et 0,2 L (hors de la couche lente) et prolongée
    // jusqu'à la paroi : u_r(y), et δ* = ∫ (1 − u/u_r) dy. (Le glissement
    // libre ne convient pas comme référence : sans condition de Kutta au bord
    // de fuite arrondi, sa circulation est arbitraire.)
    let refp: Vec<Vec<(f64, f64)>> = profs
        .iter()
        .map(|pr: &Vec<(f64, f64)>| {
            let pts: Vec<&(f64, f64)> = pr.iter().filter(|p| p.0 >= 0.1 * lf * s).collect();
            let n = pts.len().max(1) as f64;
            let (my, mu) = (pts.iter().map(|p| p.0).sum::<f64>() / n, pts.iter().map(|p| p.1).sum::<f64>() / n);
            let sxy: f64 = pts.iter().map(|p| (p.0 - my) * (p.1 - mu)).sum();
            let sxx: f64 = pts.iter().map(|p| (p.0 - my).powi(2)).sum::<f64>().max(1e-12);
            let b = sxy / sxx;
            pr.iter().map(|p| (p.0, mu + b * (p.0 - my))).collect()
        })
        .collect();
    let nuf = nu * s; // viscosité en unités de la grille la plus fine
    struct St {
        xl: f64,
        ue: f64,
        ds: f64,
        th: f64,
        ds_nodes: f64,
        utau: f64,
        prof: Vec<(f64, f64, f64)>,
    }
    let mut res: Vec<St> = Vec::new();
    for ((&xl, pr), rp) in stations.iter().zip(&profs).zip(&refp) {
        if pr.len() < 4 || rp.len() != pr.len() {
            continue;
        }
        let pr: Vec<(f64, f64, f64)> = pr.iter().zip(rp).map(|(a, b)| (a.0, a.1, b.1)).collect();
        // Vitesse extérieure : celle de la référence au premier nœud.
        let ue = pr[0].2;
        let uh = {
            let k = pr.iter().position(|p| p.0 >= WM_H).unwrap_or(1).max(1);
            let (p0, p1) = (pr[k - 1], pr[k]);
            p0.1 + (p1.1 - p0.1) * (WM_H - p0.0) / (p1.0 - p0.0)
        };
        let utau = lbm::gpu3d::spalding_utau(uh, WM_H, nuf);
        // Sous le premier nœud : profil de Spalding de ce u_τ ; au-dessus,
        // trapèzes sur les nœuds rapportés à la référence.
        let (mut ds, mut th) = (0.0, 0.0);
        let d1 = pr[0].0;
        let m = 200;
        for q in 0..m {
            let y = d1 * (q as f64 + 0.5) / m as f64;
            let r = (utau * spalding_uplus(y * utau / nuf) / ue).min(1.0);
            ds += (1.0 - r) * d1 / m as f64;
            th += r * (1.0 - r) * d1 / m as f64;
        }
        // « Nœuds seuls » : u(0) extrapolé linéairement des deux premiers nœuds.
        let u0 = (pr[0].1 - (pr[1].1 - pr[0].1) / (pr[1].0 - pr[0].0) * pr[0].0).max(0.0);
        let mut ds_nodes = 0.5 * d1 * (2.0 - (u0 + pr[0].1) / ue);
        for w in pr.windows(2) {
            let (r0, r1) = (w[0].1 / w[0].2, w[1].1 / w[1].2);
            let h = w[1].0 - w[0].0;
            ds += 0.5 * h * ((1.0 - r0) + (1.0 - r1));
            th += 0.5 * h * (r0 * (1.0 - r0) + r1 * (1.0 - r1));
            ds_nodes += 0.5 * h * ((1.0 - r0) + (1.0 - r1));
        }
        res.push(St { xl, ue, ds: ds / s, th: th / s, ds_nodes: ds_nodes / s, utau, prof: pr });
    }

    let ui = u as f64;
    println!(
        "Plaque plane L = {l} ×{} ({} grilles fines), Re_L = {re:.0e}, t = {t_end} L/U, essai {essai}, calage {cale}, {} · {:.0} s · grille fine {}×{}×{}",
        1 << raff,
        raff,
        if fp16 { "demi-précision" } else { "simple précision" },
        el,
        d[0],
        d[1],
        d[2]
    );
    println!("δ*, θ en mailles de la veine ; théorie : 0,046 et 0,036 x Re_x^−0,2, Cf = 0,026 Re_x^−1/7");
    println!(" x/L   Re_x    Ue/U    δ*    th.   δ*nœuds   θ    th.     H     Cf loi  Cf bilan  Cf th.");
    for (k, r) in res.iter().enumerate() {
        let rex = re * r.xl;
        let ds_th = 0.046 * r.xl * rex.powf(-0.2) * lf;
        let th_th = 0.036 * r.xl * rex.powf(-0.2) * lf;
        let cf_th = 0.026 * rex.powf(-1.0 / 7.0);
        let cf_wl = 2.0 * (r.utau / r.ue).powi(2);
        // Bilan : différences centrées (x en mailles de la veine).
        let cf_mom = if k > 0 && k + 1 < res.len() {
            let (p, q) = (&res[k - 1], &res[k + 1]);
            let dx = (q.xl - p.xl) * lf;
            let dth = (q.th - p.th) / dx;
            let due = (q.ue - p.ue) / dx;
            let h = r.ds / r.th;
            format!("{:8.5}", 2.0 * (dth + (2.0 + h) * r.th / r.ue * due))
        } else {
            "    —   ".into()
        };
        println!(
            " {:.2}  {:.2e} {:.4} {:6.3} {:5.3}  {:6.3}  {:6.3} {:5.3} {:5.2} {:.5}  {cf_mom}  {:.5}",
            r.xl,
            rex,
            r.ue / ui,
            r.ds,
            ds_th,
            r.ds_nodes,
            r.th,
            th_th,
            r.ds / r.th,
            cf_wl,
            cf_th
        );
    }
    // Traînée (échange de quantité de mouvement, deux faces) ; frottement
    // moyen de la loi 1/7 : 2 × (7/6) 0,026 Re_L^−1/7.
    let qa = 0.5 * ui * ui * body.area;
    let cf_wl_moy = {
        let v: Vec<f64> = res.iter().filter(|r| r.xl > 0.29 && r.xl < 0.71).map(|r| 2.0 * (r.utau / r.ue).powi(2) * (r.ue / ui).powi(2)).collect();
        v.iter().sum::<f64>() / v.len().max(1) as f64
    };
    println!(
        "Frottement de l'extrados entre 0,3 et 0,7 L (rapporté à ½ρU²) : échange {:.5}, loi de paroi ρu_τ² {:.5}",
        cf_ex, cf_wl_moy
    );
    println!(
        "Cx échange : {:.5} (frottement seul en théorie : {:.5})",
        fx / qa,
        2.0 * 7.0 / 6.0 * 0.026 * re.powf(-1.0 / 7.0)
    );
    // Profil à mi-plaque contre Coles (loi log + sillage) avec le Cf et le δ
    // théoriques : U_e⁺ = ln(δ⁺)/κ + B + 2Π/κ.
    if let Some(r) = res.iter().min_by(|p, q| (p.xl - 0.6).abs().partial_cmp(&(q.xl - 0.6).abs()).unwrap()) {
        let rex = re * r.xl;
        let cf = 0.026 * rex.powf(-1.0 / 7.0);
        let ut = (0.5 * cf).sqrt(); // u_τ / U_e
        let uep = 1.0 / ut;
        let ln_dp = KAPPA * (uep - B_LOG - 2.0 * PI_COLES / KAPPA);
        let delta_p = ln_dp.exp(); // δ⁺
        let delta = delta_p * nu / (ut * r.ue); // mailles de la veine
        println!(
            "Profil à x/L = {:.2} : u/U_e calculé et Coles (δ = {:.2} mailles, Π = {PI_COLES}, Cf = {:.5}) ; y en mailles de la veine",
            r.xl, delta, cf
        );
        let mut line = String::new();
        for &(y, uu, ur) in r.prof.iter().take(14) {
            let yb = y / s;
            let yp = yb * ut * r.ue / nu;
            let coles = if yb < delta {
                ut * (spalding_uplus(yp) + 2.0 * PI_COLES / KAPPA * (std::f64::consts::FRAC_PI_2 * yb / delta).sin().powi(2))
            } else {
                1.0
            };
            let _ = write!(line, "  y {:5.2} : {:.3} / {:.3}\n", yb, uu / ur, coles.min(1.0));
        }
        print!("{line}");
    }
    if let Some(path) = csv {
        let mut out = String::from("x_L,y_mailles_veine,u_Ue,u_ref_Ue\n");
        for r in &res {
            for &(y, uu, ur) in &r.prof {
                let _ = writeln!(out, "{:.3},{:.4},{:.5},{:.5}", r.xl, y / s, uu / r.ue, ur / r.ue);
            }
        }
        std::fs::write(&path, out).unwrap();
    }
}
