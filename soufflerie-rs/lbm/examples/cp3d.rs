//! Diagnostic 3D : Cp pariétal d'une maquette STL dans une tranche d'envergure.
//! cp3d <stl> <unite> <lref> <sref> <L> <raff> <re> <alpha> <y/b centre> <t> <avg> [demi-précision 1/0]
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
    let (l, raff, re, eta, t, avg) = (g(5, 16.0) as usize, g(6, 2.0) as usize, g(7, 1e6), g(9, 0.0), g(10, 20.0), g(11, 6.0));
    let (body, _) = Body3D::from_mesh(&tris, l, &setup).unwrap();
    let u = 0.1f32;
    let edges = std::env::var("SOUFFLERIE_BORDS").unwrap_or_default();
    let mut sim = Gpu3D::with_refinement_edges(&body, u as f64 * l as f64 / re, u, g(12, 1.0) > 0.5, true, raff, edges.contains('a'), edges.contains('f'), &std::env::var("SOUFFLERIE_BANDE").unwrap_or_default().split(',').filter_map(|v| v.parse().ok()).collect::<Vec<f64>>()).unwrap();
    let spt = l as f64 / u as f64;
    let ramp = lbm::cpu::Ramp { steps: (4.0 * spt) as u64, u };
    sim.run(((t - avg) * spt) as usize, &ramp);
    let _ = sim.take_forces();
    let n = (2.0 * avg) as usize;
    let mut acc: Vec<f64> = Vec::new();
    let mut last = None;
    for _ in 0..n {
        sim.run((0.5 * spt) as usize, &ramp);
        let s = sim.surface_forces().unwrap();
        if acc.is_empty() {
            acc = vec![0.0; s.walls.len()];
        }
        for (k, w) in s.walls.iter().enumerate() {
            acc[k] += w.1 / n as f64;
        }
        last = Some(s);
    }
    // Traînée par tranche de corde locale : Fx de chaque triangle (pression +
    // frottement, dernier relevé), x/c dans sa bande d'envergure.
    {
        let s = last.as_ref().unwrap();
        let nb = 24;
        let ys: Vec<f64> = s.loads.iter().map(|l| l.0[1]).collect();
        let (y0, y1) = ys.iter().fold((f64::MAX, f64::MIN), |m, &y| (m.0.min(y), m.1.max(y)));
        let mut strip = vec![(f64::MAX, f64::MIN); nb];
        let sb = |y: f64| (((y - y0) / (y1 - y0 + 1e-9)) * nb as f64) as usize;
        for l in &s.loads {
            let k = sb(l.0[1]).min(nb - 1);
            strip[k] = (strip[k].0.min(l.0[0]), strip[k].1.max(l.0[0]));
        }
        let bins = [0.0, 0.02, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.97, 1.01];
        let mut fx = vec![[0.0f64; 2]; bins.len() - 1];
        for (l, w) in s.loads.iter().zip(&s.walls) {
            let k = sb(l.0[1]).min(nb - 1);
            let xc = (l.0[0] - strip[k].0) / (strip[k].1 - strip[k].0).max(1e-9);
            let b = bins.windows(2).position(|w| xc >= w[0] && xc < w[1]).unwrap_or(bins.len() - 2);
            fx[b][if w.0[2] > 0.0 { 0 } else { 1 }] += l.1[0];
        }
        let qa = 0.5 * (u as f64).powi(2) * body.area;
        println!("Cx par tranche de corde (extrados / intrados) :");
        for (w, f) in bins.windows(2).zip(&fx) {
            println!("  {:.2}-{:.2}  {:+.4} {:+.4}", w[0], w[1], f[0] / qa, f[1] / qa);
        }
    }

    // Profils de pression le long de la normale (extrados, mi-tranche).
    let q = 0.5 * (u as f64).powi(2);
    {
        let (mac, d, org, sc) = sim.finest_macros();
        let s0 = sim.surface_forces().unwrap();
        let samp = |p: [f64; 3]| -> Option<f64> {
            let i = p.map(|v| v.floor() as i64);
            let w = [0, 1, 2].map(|a| p[a] - i[a] as f64);
            let (mut acc, mut ws) = (0.0, 0.0);
            for c in 0..8 {
                let o = [c & 1, (c >> 1) & 1, (c >> 2) & 1];
                let q = [i[0] + o[0] as i64, i[1] + o[1] as i64, i[2] + o[2] as i64];
                if (0..3).any(|a| q[a] < 0 || q[a] >= d[a] as i64) { return None; }
                let m = mac[q[0] as usize + d[0] * (q[1] as usize + d[1] * q[2] as usize)];
                if m[3] == 0.0 || (m[0] == 0.0 && m[1] == 0.0 && m[2] == 0.0) { continue; }
                let ww: f64 = (0..3).map(|a| if o[a] == 1 { w[a] } else { 1.0 - w[a] }).product();
                acc += ww * (m[3] as f64 - 1.0) / 3.0;
                ws += ww;
            }
            if ws < 0.3 { None } else { Some(acc / ws) }
        };
        let ys: Vec<f64> = s0.loads.iter().map(|l| l.0[1]).collect();
        let (y0, y1) = ys.iter().fold((f64::MAX, f64::MIN), |m, &y| (m.0.min(y), m.1.max(y)));
        let yc = 0.5 * (y0 + y1) + eta * 0.5 * (y1 - y0);
        let band = 0.02 * (y1 - y0);
        let sel: Vec<usize> = (0..ys.len()).filter(|&k| (ys[k] - yc).abs() < band).collect();
        let xs: Vec<f64> = sel.iter().map(|&k| s0.loads[k].0[0]).collect();
        let (xa, xb) = xs.iter().fold((f64::MAX, f64::MIN), |m, &x| (m.0.min(x), m.1.max(x)));
        let hs = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0];
        println!("profils Cp(n), n en mailles fines : {:?}", hs);
        let mut done = vec![false; 8];
        let tg = [0.0, 0.003, 0.006, 0.01, 0.02, 0.04, 0.08, 0.2];
        for &k in &sel {
            let xc = (s0.loads[k].0[0] - xa) / (xb - xa);
            let n = s0.walls[k].0;
            let up = n[2] > 0.0;
            let Some(b) = tg.iter().position(|&t| (xc - t).abs() < 0.0025 + 0.05 * t) else { continue };
            if done[b] || !(up || b == 0) { continue; }
            done[b] = true;
            let c = [0, 1, 2].map(|a| (s0.loads[k].0[a] - org[a]) * sc);
            let v: Vec<String> = hs.iter().map(|&h| match samp([0, 1, 2].map(|a| c[a] + h * n[a])) { Some(p) => format!("{:6.2}", p / q), None => "   -  ".into() }).collect();
            println!("x/c {:.4} n=({:+.2},{:+.2}) pw {:6.2} | {}", xc, n[0], n[2], s0.walls[k].1 / q, v.join(" "));
        }
    }
    let f = sim.take_forces();
    let q = 0.5 * (u as f64).powi(2);
    println!("Cz échange {:.4}", f.fz / (q * body.area));
    let s = last.unwrap();
    // Tranche d'envergure : centres de triangles dans |y − y0| < 0,05 b.
    let ys: Vec<f64> = s.loads.iter().map(|l| l.0[1]).collect();
    let (y0, y1) = ys.iter().fold((f64::MAX, f64::MIN), |m, &y| (m.0.min(y), m.1.max(y)));
    let yc = 0.5 * (y0 + y1) + eta * 0.5 * (y1 - y0);
    let band = 0.05 * (y1 - y0);
    let sel: Vec<usize> = (0..ys.len()).filter(|&k| (ys[k] - yc).abs() < band).collect();
    let xs: Vec<f64> = sel.iter().map(|&k| s.loads[k].0[0]).collect();
    let (xa, xb) = xs.iter().fold((f64::MAX, f64::MIN), |m, &x| (m.0.min(x), m.1.max(x)));
    let mut pts: Vec<(f64, bool, f64, f64)> = sel
        .iter()
        .map(|&k| ((s.loads[k].0[0] - xa) / (xb - xa), s.walls[k].0[2] > 0.0, acc[k] / q, s.walls[k].0[0]))
        .collect();
    pts.sort_by(|a, b| a.0.total_cmp(&b.0));
    let bins = [0.0, 0.01, 0.03, 0.06, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 1.01];
    println!("tranche y/b = {eta} ({} triangles)\nx/c          Cp ext   Cp int", sel.len());
    for w in bins.windows(2) {
        let m = |up: bool| {
            let v: Vec<f64> = pts.iter().filter(|p| p.1 == up && p.0 >= w[0] && p.0 < w[1]).map(|p| p.2).collect();
            if v.is_empty() { f64::NAN } else { v.iter().sum::<f64>() / v.len() as f64 }
        };
        println!("{:.2}-{:.2}  {:8.3} {:8.3}", w[0], w[1], m(true), m(false));
    }
    for up in [true, false] {
        let line: Vec<String> = pts.iter().filter(|p| p.1 == up).map(|p| format!("{:.3}:{:.2}", p.0, p.2)).collect();
        println!("{} {}", if up { "ext" } else { "int" }, line.join(" "));
    }
}
