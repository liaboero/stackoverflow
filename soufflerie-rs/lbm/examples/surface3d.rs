//! Diagnostic : efforts par intégration sur la surface du maillage STL —
//! pression relevée à h mailles de la paroi, frottement par la loi de Spalding
//! sur la vitesse tangentielle au même point — comparés à l'échange de
//! quantité de mouvement.
use lbm::bodies3d::Body3D;
use lbm::gpu3d::Gpu3D;
use lbm::stl::{MeshSetup, Unit};

fn spalding_utau(ut: f64, y: f64, nu: f64) -> f64 {
    if ut <= 1e-12 {
        return 0.0;
    }
    let (k, eb) = (0.41, 0.12873490);
    let rey = y * ut / nu;
    let mut up = rey.sqrt();
    if rey >= 30.0 {
        up = rey.ln() / k;
    }
    for _ in 0..50 {
        let ku = k * up;
        let e = ku.exp();
        let yy = up + eb * (e - 1.0 - ku - ku * ku / 2.0 - ku * ku * ku / 6.0);
        let dy = 1.0 + eb * k * (e - 1.0 - ku - ku * ku / 2.0);
        let h = up * yy - rey;
        let next = up - h / (yy + up * dy);
        up = if next > 0.0 { next } else { 0.5 * up };
        if h.abs() < 1e-10 * rey {
            break;
        }
    }
    ut / up
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let path = &a[1];
    let alpha: f64 = a.get(2).and_then(|s| s.parse().ok()).unwrap_or(4.0);
    let l: usize = a.get(3).and_then(|s| s.parse().ok()).unwrap_or(32);
    let re: f64 = a.get(4).and_then(|s| s.parse().ok()).unwrap_or(1e6);
    let t: f64 = a.get(5).and_then(|s| s.parse().ok()).unwrap_or(12.0);
    let tris = lbm::stl::read_stl(std::path::Path::new(path)).unwrap();
    let setup = MeshSetup { unit: Unit::M, alpha_deg: alpha, ref_length: Some(1.0), ref_area: Some(4.0), margin: 0.4, ..Default::default() };
    let (body, _) = Body3D::from_mesh(&tris, l, &setup).unwrap();
    let u = 0.1f32;
    let nu = u as f64 * l as f64 / re;
    let spt = l as f64 / u as f64;
    let mut g = Gpu3D::with_refinement(&body, nu, u, true, true, 0).unwrap();
    let ramp = lbm::cpu::Ramp { steps: (4.0 * spt) as u64, u };
    g.run((t * spt) as usize, &ramp);
    let _ = g.take_forces();
    g.run((2.0 * spt) as usize, &ramp);
    let f = g.take_forces();
    let q = 0.5 * (u as f64).powi(2) * body.area;
    let mac = g.macros();
    let (nx, ny, nz) = (g.nx, g.ny, g.nz);
    let m = body.mesh.as_ref().unwrap();
    // Orientation : normales sortantes si Σ c·n A > 0.
    let tri_n = |t: &[[f64; 3]; 3]| {
        let e1 = [t[1][0] - t[0][0], t[1][1] - t[0][1], t[1][2] - t[0][2]];
        let e2 = [t[2][0] - t[0][0], t[2][1] - t[0][1], t[2][2] - t[0][2]];
        [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]]
    };
    let vol: f64 = m.tris.iter().map(|t| { let n = tri_n(t); (0..3).map(|k| t[0][k] * n[k]).sum::<f64>() }).sum();
    let sgn = if vol > 0.0 { 1.0 } else { -1.0 };
    let sample = |p: [f64; 3]| -> Option<[f64; 4]> {
        let i = [p[0].floor() as i64, p[1].floor() as i64, p[2].floor() as i64];
        let w = [p[0] - i[0] as f64, p[1] - i[1] as f64, p[2] - i[2] as f64];
        let mut out = [0.0; 4];
        for c in 0..8 {
            let o = [c & 1, (c >> 1) & 1, (c >> 2) & 1];
            let (x, y, z) = (i[0] + o[0] as i64, i[1] + o[1] as i64, i[2] + o[2] as i64);
            if x < 0 || y < 0 || z < 0 || x >= nx as i64 || y >= ny as i64 || z >= nz as i64 { return None; }
            if body.inside(x as f64, y as f64, z as f64) { return None; }
            let ww: f64 = (0..3).map(|k| if o[k] == 1 { w[k] } else { 1.0 - w[k] }).product();
            let v = mac[x as usize + nx * (y as usize + ny * z as usize)];
            for k in 0..4 { out[k] += ww * v[k] as f64; }
        }
        Some(out)
    };
    for h in [1.0f64, 1.5, 2.0, 3.0] {
        let (mut fp, mut ff, mut area, mut miss) = ([0.0f64; 3], [0.0f64; 3], 0.0, 0.0);
        for t in &m.tris {
            let n0 = tri_n(t);
            let a2 = (n0[0] * n0[0] + n0[1] * n0[1] + n0[2] * n0[2]).sqrt();
            if a2 < 1e-12 { continue; }
            let n = n0.map(|v| sgn * v / a2);
            let ar = 0.5 * a2;
            area += ar;
            let c = [0, 1, 2].map(|k| (t[0][k] + t[1][k] + t[2][k]) / 3.0);
            let mut s = None;
            let mut hh = h;
            while s.is_none() && hh < h + 3.0 {
                s = sample([0, 1, 2].map(|k| c[k] + hh * n[k]));
                if s.is_none() { hh += 0.5; }
            }
            let Some(v) = s else { miss += ar; continue; };
            let p = (v[3] - 1.0) / 3.0;
            let un = v[0] * n[0] + v[1] * n[1] + v[2] * n[2];
            let ut = [v[0] - un * n[0], v[1] - un * n[1], v[2] - un * n[2]];
            let utm = (ut[0] * ut[0] + ut[1] * ut[1] + ut[2] * ut[2]).sqrt();
            let ut_tau = spalding_utau(utm, hh, nu);
            let tau = v[3] * ut_tau * ut_tau;
            for k in 0..3 {
                fp[k] += -p * n[k] * ar;
                if utm > 0.0 { ff[k] += tau * ut[k] / utm * ar; }
            }
        }
        println!(
            "h {h} : Cx pression {:.4} frottement {:.4} total {:.4} · Cz pression {:.4} frottement {:.4} total {:.4} · surface {:.0} mailles² (non relevée {:.1} %)",
            fp[0] / q, ff[0] / q, (fp[0] + ff[0]) / q, fp[2] / q, ff[2] / q, (fp[2] + ff[2]) / q, area, 100.0 * miss / area
        );
    }
    println!("échange : Cx {:.4} Cz {:.4}", f.fx / q, f.fz / q);
}
