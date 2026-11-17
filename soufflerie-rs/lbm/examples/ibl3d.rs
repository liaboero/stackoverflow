//! Diagnostic 3D : couche limite intégrale sur le maillage (décambrage par
//! transpiration), relevée sur le champ d'un calcul établi, sans l'imposer.
//! ibl3d <stl> <unite> <lref> <sref> <L> <raff> <re> <alpha> <t> <csv> [xtr] [couplages] [relaxation] [mode 1|2]
//! CSV par triangle : centre, normale, u_e, θ, H, δ*, v_n (unités de la grille
//! entière la plus fine), laminaire, décollé.
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
        ..Default::default()
    };
    let (l, raff, re, t) = (g(5, 16.0) as usize, g(6, 1.0) as usize, g(7, 1e6), g(9, 10.0));
    let (body, _) = Body3D::from_mesh(&tris, l, &setup).unwrap();
    let u = 0.1f32;
    let mut sim = Gpu3D::with_refinement(&body, u as f64 * l as f64 / re, u, true, true, raff).unwrap();
    let spt = l as f64 / u as f64;
    let ramp = lbm::cpu::Ramp { steps: (4.0 * spt) as u64, u };
    sim.run((t * spt) as usize, &ramp);
    sim.enable_transpiration(g(11, 0.0), g(14, 1.0) > 1.5).unwrap();
    let r = sim.couple_transpiration(0.0).unwrap();
    eprintln!(
        "v_n/u∞ moyenne {:.5} max {:.4} · sillage {:.5} · décollé {:.1} % · non relevé {:.1} % · L = {}",
        r.vn_mean,
        r.vn_max,
        r.outflow,
        100.0 * r.separated,
        100.0 * r.missed,
        r.l
    );
    dump(&r, &a[10]);
    // Profils de vitesse tangentielle le long de la normale, extrados à
    // mi-envergure (grille entière la plus fine), à quelques x/c.
    {
        let (mac, d, _, _) = sim.finest_macros();
        let samp = |p: [f64; 3]| -> Option<[f64; 3]> {
            let i = p.map(|v| v.floor() as i64);
            let w = [0, 1, 2].map(|a| p[a] - i[a] as f64);
            let mut o = [0.0; 3];
            for c in 0..8 {
                let q = [i[0] + (c & 1) as i64, i[1] + ((c >> 1) & 1) as i64, i[2] + ((c >> 2) & 1) as i64];
                if (0..3).any(|a| q[a] < 0 || q[a] >= d[a] as i64) {
                    return None;
                }
                let m = mac[q[0] as usize + d[0] * (q[1] as usize + d[1] * q[2] as usize)];
                let ww: f64 = (0..3).map(|a| if (c >> a) & 1 == 1 { w[a] } else { 1.0 - w[a] }).product();
                for k in 0..3 {
                    o[k] += ww * m[k] as f64;
                }
            }
            Some(o)
        };
        let ys: Vec<f64> = r.centers.iter().map(|c| c[1]).collect();
        let yc = 0.5 * (ys.iter().cloned().fold(f64::MAX, f64::min) + ys.iter().cloned().fold(f64::MIN, f64::max));
        let sel: Vec<usize> = (0..r.centers.len()).filter(|&k| (r.centers[k][1] - yc).abs() < 2.0).collect();
        let x0 = sel.iter().map(|&k| r.centers[k][0]).fold(f64::MAX, f64::min);
        let x1 = sel.iter().map(|&k| r.centers[k][0]).fold(f64::MIN, f64::max);
        for side in [1.0, -1.0] {
            for xc in [0.3, 0.6, 0.8, 0.9, 0.95] {
                let Some(&k) = sel
                    .iter()
                    .filter(|&&k| r.normals[k][2] * side > 0.0)
                    .min_by(|&&a, &&b| {
                        let f = |k: usize| ((r.centers[k][0] - x0) / (x1 - x0) - xc).abs();
                        f(a).partial_cmp(&f(b)).unwrap()
                    })
                else {
                    continue;
                };
                let (c, n) = (r.centers[k], r.normals[k]);
                let mut line = format!("{} x/c {:.2} u_t/U :", if side > 0.0 { "ext" } else { "int" }, (c[0] - x0) / (x1 - x0));
                for h in [1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0] {
                    if let Some(v) = samp([0, 1, 2].map(|a| c[a] + h * n[a])) {
                        let vn = v[0] * n[0] + v[1] * n[1] + v[2] * n[2];
                        let ut = (0..3).map(|a| (v[a] - vn * n[a]).powi(2)).sum::<f64>().sqrt();
                        line += &format!(" {h}:{:.3}", ut / u as f64);
                    }
                }
                eprintln!("{line}");
            }
        }
    }
    // Couplages imposés ensuite : <n> tous les 0,5 temps convectif, relaxation <ω>.
    let (nc, relax) = (g(12, 0.0) as usize, g(13, 0.5));
    let q = 0.5 * (u as f64).powi(2) * body.area;
    let _ = sim.take_forces();
    for k in 0..nc {
        let r = sim.couple_transpiration(relax).unwrap();
        sim.run((0.5 * spt) as usize, &ramp);
        let f = sim.take_forces();
        let bx = sim.box_forces_vol(&[0.1]);
        let bz = bx[0].map_or(f64::NAN, |(f, _)| f[2] / q);
        println!(
            "{k:3} t* = {:5.1}  Cz échange {:.4}  boîte {:.4}  Cx {:.4} · v_n moy {:.5} max {:.4} · sillage {:.5} · décollé {:.1} % · limité {:.1} %",
            sim.time as f64 / spt,
            f.fz / q,
            bz,
            f.fx / q,
            r.vn_mean,
            r.vn_max,
            r.outflow,
            100.0 * r.separated,
            100.0 * r.layer.limited.iter().zip(&r.area).filter(|x| *x.0).map(|x| x.1).sum::<f64>() / r.area.iter().sum::<f64>()
        );
        if !sim.is_finite(&f) {
            break;
        }
        if k + 1 == nc {
            dump(&sim.couple_transpiration(relax).unwrap(), &format!("{}.fin", a[10]));
        }
    }
}

fn dump(r: &lbm::gpu3d::TranspirationReport, path: &str) {
    let ly = &r.layer;
    let mut csv = String::from("x,y,z,nx,ny,nz,ux,uy,uz,area,theta,h,dstar,vn,lam,sep,dres,he\n");
    for k in 0..r.centers.len() {
        let (c, n, v) = (r.centers[k], r.normals[k], r.ue[k]);
        csv += &format!(
            "{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.5},{:.5},{:.5},{:.5},{:.6e},{:.4},{:.6e},{:.6e},{},{},{:.4},{:.2}\n",
            c[0],
            c[1],
            c[2],
            n[0],
            n[1],
            n[2],
            v[0],
            v[1],
            v[2],
            r.area[k],
            ly.theta[k],
            ly.h[k],
            ly.dstar[k],
            r.vn[k],
            ly.laminar[k] as u8,
            ly.separated[k] as u8,
            r.dstar_res[k],
            r.h_edge[k]
        );
    }
    std::fs::write(path, csv).unwrap();
}
