//! Validation aéroacoustique 2D : cylindre en écoulement uniforme (Inoue &
//! Hatakeyama 2002 : Re = 150, M = 0,2). Rayonnement direct relevé par des
//! sondes GPU, comparé à Ffowcs Williams–Hawkings (surfaces poreuses) et à
//! Curle compact (efforts), tous deux calculés en post-traitement (`lbm::fwh`).
//!
//! acoustique2d [--d 16] [--re 150] [--mach 0.2] [--rmax 75] [--couche 30]
//!              [--t 300] [--rec 120] [--every 8] [--bgk] [--sans-couche] [--sortie DIR]
//! acoustique2d --relire DIR      (analyse seule des relevés enregistrés)
//!
//! Domaine carré, cylindre au centre, couches absorbantes de `couche` D sur
//! les quatre côtés (rappel vers l'écoulement amont) ; longueurs en D, temps
//! en D/U. Relevés : cercles de sondes (r = 5 … rmax D, tous les 10°), ligne
//! θ = 90° (tous les D), surfaces FW-H (cercles de rayon 0,75 / 1,5 / 3 / 6 D,
//! 256 points). θ est compté depuis l'aval (+x), sens trigonométrique.

use lbm::cpu::{Collision, LesMode, Ramp, SolverConfig};
use lbm::fwh::{curle_compact, dft_amp, fwh_pressure, fwh_sources, refine_peak, Cx, Medium, Surface};
use lbm::gpu::GpuSolver;
use lbm::shapes::Body;
use lbm::{ShapeKind, ShapeParams};
use std::f64::consts::PI;
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};

const RINGS: [f64; 8] = [5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 75.0];
const SURF: [f64; 4] = [0.75, 1.5, 3.0, 6.0];
const NSURF: usize = 256;

#[derive(Clone, Copy)]
struct Probe {
    x: f64,
    y: f64,
    /// 0 : cercle d'observation, 1 : ligne θ = 90°, 2 : surface FW-H.
    kind: u8,
    /// r/D (0, 1) ou indice de surface (2).
    a: f64,
    /// θ en degrés (0, 1) ou indice du point (2).
    b: f64,
}

struct Run {
    d: f64,
    re: f64,
    u: f64,
    cx: f64,
    cy: f64,
    every: u64,
    probes: Vec<Probe>,
    /// Pas des relevés, valeurs [relevé][sonde][ρ, u_x, u_y].
    times: Vec<u64>,
    data: Vec<f32>,
    /// (pas de fin, pas cumulés, Fx, Fy) par relevé d'efforts.
    forces: Vec<(u64, u64, f64, f64)>,
}

fn arg<T: std::str::FromStr>(a: &[String], k: &str, d: T) -> T {
    a.iter()
        .position(|s| s == k)
        .and_then(|i| a.get(i + 1))
        .and_then(|s| s.parse().ok())
        .unwrap_or(d)
}

fn layout(d: f64, rmax: f64, cx: f64, cy: f64) -> Vec<Probe> {
    let mut p = Vec::new();
    for &r in RINGS.iter().filter(|&&r| r <= rmax + 1e-9) {
        for k in 0..36 {
            let th = 10.0 * k as f64;
            let t = th.to_radians();
            p.push(Probe { x: cx + r * d * t.cos(), y: cy + r * d * t.sin(), kind: 0, a: r, b: th });
        }
    }
    let mut r = 1.0;
    while r <= rmax + 1e-9 {
        p.push(Probe { x: cx, y: cy + r * d, kind: 1, a: r, b: 90.0 });
        r += 1.0;
    }
    for (si, &rs) in SURF.iter().enumerate() {
        let s = Surface::circle([cx, cy], rs * d, NSURF);
        for (j, q) in s.pts.iter().enumerate() {
            p.push(Probe { x: q[0], y: q[1], kind: 2, a: si as f64, b: j as f64 });
        }
    }
    p
}

fn simulate(a: &[String], out: &Path) -> Run {
    let d: usize = arg(a, "--d", 16);
    let re: f64 = arg(a, "--re", 150.0);
    let mach: f64 = arg(a, "--mach", 0.2);
    let rmax: f64 = arg(a, "--rmax", 75.0);
    let layer: f64 = arg(a, "--couche", 30.0);
    let t_end: f64 = arg(a, "--t", 300.0);
    let rec: f64 = arg(a, "--rec", 120.0);
    let every: u64 = arg(a, "--every", 8);
    let bgk = a.iter().any(|s| s == "--bgk");
    // Bords de la veine ordinaire : bandes haute et basse de ny/8, éponge
    // visqueuse de sortie sur 12 %, entrée et sortie de Guo sans rappel.
    let plain = a.iter().any(|s| s == "--sans-couche");
    let u = mach / 3f64.sqrt();
    let df = d as f64;
    let half = ((rmax + 5.0 + layer) * df).ceil() as usize;
    let (nx, ny) = (2 * half + 1, 2 * half + 1);
    let w = (layer * df).round() as usize;
    let cfg = SolverConfig {
        nx,
        ny,
        u_in: u as f32,
        nu: u * df / re,
        les: LesMode::Off,
        cs: 0.15,
        les_wall_a: 1.0,
        collision: if bgk { Collision::Bgk } else { Collision::Rr },
        bouzidi: true,
        // Viscosité accrue dans la couche de sortie seulement.
        sponge_fraction: if plain { 0.12 } else { w as f64 / nx as f64 },
        absorb: if plain { None } else { Some(w) },
        half: false,
    };
    let mut g = GpuSolver::new(cfg).expect("GPU");
    let body = Body::new(ShapeKind::Cylinder, &ShapeParams::default(), df, half as f64, half as f64, u);
    g.set_body(&body);
    g.reset(0.0);
    let (cx, cy) = (body.ref_x, body.ref_y);
    let spt = df / u;
    let ramp = Ramp { steps: (5.0 * spt) as u64, u: u as f32 };
    eprintln!(
        "{} · {nx}×{ny} ({:.1} M mailles), D = {d}, Re = {re}, M = {mach}, τ = {:.4}, {:.0} pas",
        g.adapter_name(),
        (nx * ny) as f64 / 1e6,
        g.tau0(),
        t_end * spt
    );
    let probes = layout(df, rmax, cx, cy);
    let total = (t_end * spt) as u64;
    let rec_start = total.saturating_sub((rec * spt) as u64) / every * every;
    let mut forces = Vec::new();
    let clock = std::time::Instant::now();
    let mut armed = false;
    while g.time < total {
        if !armed && g.time >= rec_start {
            let pos: Vec<[f32; 2]> = probes.iter().map(|p| [p.x as f32, p.y as f32]).collect();
            g.set_probes(&pos, every);
            armed = true;
        }
        g.run(every as usize, &ramp);
        for (t, f) in g.poll_forces() {
            forces.push((t, f.steps, f.fx, f.fy));
        }
        if g.time % (20 * spt as u64).max(every) < every {
            eprintln!("t = {:.0} D/U ({:.0} s)", g.time as f64 / spt, clock.elapsed().as_secs_f64());
            if !g.is_finite() {
                panic!("calcul divergent");
            }
        }
    }
    let _ = g.take_forces();
    let el = clock.elapsed().as_secs_f64();
    eprintln!("{:.0} s, {:.0} M mailles/s", el, (nx * ny) as f64 * total as f64 / el / 1e6);
    let (times, data) = g.take_probes();
    let r = Run { d: df, re, u, cx, cy, every, probes, times, data, forces };
    save(&r, out, el);
    r
}

fn save(r: &Run, out: &Path, elapsed: f64) {
    std::fs::create_dir_all(out).unwrap();
    let mut f = std::fs::File::create(out.join("parametres.txt")).unwrap();
    writeln!(f, "{} {} {} {} {} {} {}", r.d, r.re, r.u, r.cx, r.cy, r.every, elapsed).unwrap();
    let mut f = std::io::BufWriter::new(std::fs::File::create(out.join("sondes.csv")).unwrap());
    for p in &r.probes {
        writeln!(f, "{},{},{},{},{}", p.x, p.y, p.kind, p.a, p.b).unwrap();
    }
    let mut f = std::io::BufWriter::new(std::fs::File::create(out.join("efforts.csv")).unwrap());
    for e in &r.forces {
        writeln!(f, "{},{},{},{}", e.0, e.1, e.2, e.3).unwrap();
    }
    let mut f = std::io::BufWriter::new(std::fs::File::create(out.join("releves.bin")).unwrap());
    f.write_all(&(r.times.len() as u64).to_le_bytes()).unwrap();
    for t in &r.times {
        f.write_all(&t.to_le_bytes()).unwrap();
    }
    for v in &r.data {
        f.write_all(&v.to_le_bytes()).unwrap();
    }
}

fn load(dir: &Path) -> Run {
    let s = std::fs::read_to_string(dir.join("parametres.txt")).unwrap();
    let v: Vec<f64> = s.split_whitespace().map(|x| x.parse().unwrap()).collect();
    let rd = |n: &str| std::io::BufReader::new(std::fs::File::open(dir.join(n)).unwrap()).lines();
    let probes = rd("sondes.csv")
        .map(|l| {
            let c: Vec<f64> = l.unwrap().split(',').map(|x| x.parse().unwrap()).collect();
            Probe { x: c[0], y: c[1], kind: c[2] as u8, a: c[3], b: c[4] }
        })
        .collect();
    let forces = rd("efforts.csv")
        .map(|l| {
            let c: Vec<f64> = l.unwrap().split(',').map(|x| x.parse().unwrap()).collect();
            (c[0] as u64, c[1] as u64, c[2], c[3])
        })
        .collect();
    let b = std::fs::read(dir.join("releves.bin")).unwrap();
    let n = u64::from_le_bytes(b[0..8].try_into().unwrap()) as usize;
    let times = (0..n).map(|k| u64::from_le_bytes(b[8 + 8 * k..16 + 8 * k].try_into().unwrap())).collect();
    let data = b[8 + 8 * n..].chunks_exact(4).map(|c| f32::from_le_bytes(c.try_into().unwrap())).collect();
    Run { d: v[0], re: v[1], u: v[2], cx: v[3], cy: v[4], every: v[5] as u64, probes, times, data, forces }
}

fn analyze(r: &Run) {
    let m = Medium::lattice(r.u);
    let (d, u) = (r.d, r.u);
    let np = r.probes.len();
    let ns = r.times.len();
    assert!(ns > 16 && r.data.len() == ns * np * 3, "relevés incomplets");
    let t0 = r.times[0] as f64;
    let dt = r.every as f64;
    let q = 0.5 * u * u * d;
    // Efforts sur la fenêtre des relevés (instant milieu de chaque cumul).
    let fw: Vec<&(u64, u64, f64, f64)> = r.forces.iter().filter(|e| e.1 > 0 && e.0 > r.times[0]).collect();
    let tf: Vec<f64> = fw.iter().map(|e| e.0 as f64 - 0.5 * e.1 as f64).collect();
    let cd: Vec<f64> = fw.iter().map(|e| e.2 / e.1 as f64 / q).collect();
    let cl: Vec<f64> = fw.iter().map(|e| e.3 / e.1 as f64 / q).collect();
    let fdt = (tf[tf.len() - 1] - tf[0]) / (tf.len() - 1) as f64;
    let ts: Vec<f64> = tf.iter().map(|t| t * u / d).collect();
    let st0 = lbm::stats::analyze(&ts, &cd, &cl, &vec![0.0; cd.len()]);
    let f0 = refine_peak(&cl, fdt, st0.st * u / d);
    let st = f0 * d / u;
    let w0 = 2.0 * PI * f0;
    let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
    let amp = |v: &[f64]| 0.5 * (v.iter().cloned().fold(f64::MIN, f64::max) - v.iter().cloned().fold(f64::MAX, f64::min));
    let cl1 = dft_amp(&cl, tf[0], fdt, w0);
    let cd2 = dft_amp(&cd, tf[0], fdt, 2.0 * w0);
    println!("## Cylindre 2D, Re = {}, M = {:.3}, D = {} mailles", r.re, m.mach(), d);
    println!(
        "St = {st:.4} · Cd moyen = {:.4} · Cd' = {:.4} (½ crête à crête) · Cl' = {:.4} (½ crête à crête) · |Ĉl(f0)| = {:.4} · |Ĉd(2f0)| = {:.4} · {:.1} périodes relevées",
        mean(&cd),
        amp(&cd),
        amp(&cl),
        cl1.abs(),
        cd2.abs(),
        (tf[tf.len() - 1] - tf[0]) * f0
    );
    // Séries des sondes.
    let series = |k: usize, c: usize| -> Vec<f64> { (0..ns).map(|j| r.data[3 * (j * np + k) + c] as f64).collect() };
    // Pression directe p̂/(ρ0 c²) = ρ̂′.
    let direct = |k: usize, w: f64| dft_amp(&series(k, 0), t0, dt, w);
    // Sources FW-H par surface, à f0 et 2f0.
    let mut surf_names = Vec::new();
    let mut surf_src: Vec<(Surface, [Vec<[Cx; 3]>; 2])> = Vec::new();
    {
        for (si, &rs) in SURF.iter().enumerate() {
            let s = Surface::circle([r.cx, r.cy], rs * d, NSURF);
            let idx: Vec<usize> = (0..np).filter(|&k| r.probes[k].kind == 2 && r.probes[k].a as usize == si).collect();
            let get = |c: usize| -> Vec<Vec<f64>> {
                (0..ns).map(|j| idx.iter().map(|&k| r.data[3 * (j * np + k) + c] as f64).collect()).collect()
            };
            let (rho, ux, uy) = (get(0), get(1), get(2));
            let a = fwh_sources(&m, &s, &rho, &ux, &uy, t0, dt, w0);
            let b = fwh_sources(&m, &s, &rho, &ux, &uy, t0, dt, 2.0 * w0);
            surf_names.push(format!("FW-H {rs} D"));
            surf_src.push((s, [a, b]));
        }
    }
    // Efforts exercés par le fluide sur le cylindre (unités du réseau).
    let fx: Vec<f64> = cd.iter().map(|c| c * q).collect();
    let fy: Vec<f64> = cl.iter().map(|c| c * q).collect();
    let force = |w: f64| [dft_amp(&fx, tf[0], fdt, w), dft_amp(&fy, tf[0], fdt, w)];
    let (fo1, fo2) = (force(w0), force(2.0 * w0));
    let rc2 = m.rho0 * m.c * m.c;
    let pred = |x: [f64; 2], h: usize| -> (Vec<Cx>, Cx) {
        let w = if h == 0 { w0 } else { 2.0 * w0 };
        let fw = surf_src.iter().map(|(s, src)| fwh_pressure(&m, s, &src[h], w, x).scale(1.0 / rc2)).collect();
        let cu = curle_compact(&m, if h == 0 { fo1 } else { fo2 }, w, [r.cx, r.cy], x).scale(1.0 / rc2);
        (fw, cu)
    };
    for (h, name) in [(0usize, "f0 (portance)"), (1, "2f0 (traînée)")] {
        let w = if h == 0 { w0 } else { 2.0 * w0 };
        let rmax = r.probes.iter().filter(|p| p.kind == 0).map(|p| p.a).fold(0.0, f64::max);
        let mut radii = vec![rmax];
        if rmax > 40.0 {
            radii.insert(0, 40.0);
        }
        for rr in radii {
            println!("\n### Directivité à r = {rr} D, {name} : |p̂|/(ρ0c²) × 10⁵");
            println!("| θ | direct | {} | Curle compact |", surf_names.join(" | "));
            println!("|---|---|{}---|", "---|".repeat(surf_names.len()));
            // Énergie sur le cercle et profils (moyenne des deux demi-plans).
            let mut e = vec![0.0; surf_names.len() + 2];
            let mut prof = vec![[0.0f64; 2]; 19];
            for k in 0..np {
                let p = r.probes[k];
                if p.kind != 0 || (p.a - rr).abs() > 1e-9 {
                    continue;
                }
                let dv = direct(k, w).abs();
                let (fw, cu) = pred([p.x, p.y], h);
                e[0] += dv * dv;
                for (i, v) in fw.iter().chain(std::iter::once(&cu)).enumerate() {
                    e[i + 1] += v.abs().powi(2);
                }
                let b = if p.b > 180.0 { 360.0 - p.b } else { p.b };
                let j = (b / 10.0).round() as usize;
                let wgt = if j == 0 || j == 18 { 1.0 } else { 0.5 };
                prof[j][0] += wgt * dv;
                prof[j][1] += wgt * cu.abs();
                if (p.b as i64) % 30 == 0 {
                    let cells: Vec<String> = fw.iter().map(|v| format!("{:.3}", 1e5 * v.abs())).collect();
                    println!("| {:.0}° | {:.3} | {} | {:.3} |", p.b, 1e5 * dv, cells.join(" | "), 1e5 * cu.abs());
                }
            }
            // Angle du maximum (interpolation parabolique sur la grille de 10°).
            let peak = |c: usize| {
                let j = (1..18).max_by(|&a, &b| prof[a][c].total_cmp(&prof[b][c])).unwrap();
                let (a, b, cc) = (prof[j - 1][c], prof[j][c], prof[j + 1][c]);
                10.0 * (j as f64 + 0.5 * (a - cc) / (a - 2.0 * b + cc))
            };
            let rs: Vec<String> = e[1..].iter().map(|v| format!("{:.3}", (v / e[0]).sqrt())).collect();
            println!(
                "rapport des moyennes quadratiques sur le cercle, prédiction / direct ({}, Curle) : {} ; maximum (0–180°) direct {:.1}°, Curle {:.1}°",
                surf_names.join(", "),
                rs.join(", "),
                peak(0),
                peak(1)
            );
        }
        // Décroissance selon θ = 90°.
        let line: Vec<(f64, f64, Cx)> = (0..np)
            .filter(|&k| r.probes[k].kind == 1)
            .map(|k| {
                let p = r.probes[k];
                (p.a, direct(k, w).abs(), pred([p.x, p.y], h).0[1])
            })
            .collect();
        let fit: Vec<&(f64, f64, Cx)> = line.iter().filter(|l| l.0 >= 10.0 && l.0 <= 0.9 * rmax).collect();
        let (mut sx, mut sy, mut sxx, mut sxy) = (0.0, 0.0, 0.0, 0.0);
        for l in &fit {
            let (x, y) = (l.0.ln(), l.1.ln());
            sx += x;
            sy += y;
            sxx += x * x;
            sxy += x * y;
        }
        let n = fit.len() as f64;
        let slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
        let sr: Vec<f64> = fit.iter().map(|l| l.1 * l.0.sqrt()).collect();
        let (lo, hi) = sr.iter().fold((f64::MAX, f64::MIN), |a, &v| (a.0.min(v), a.1.max(v)));
        println!(
            "\n### Décroissance selon θ = 90°, {name} : |p̂|/(ρ0c²)·√(r/D) × 10⁵ ; pente log–log (10 à {:.0} D) = {slope:.3} ; ondulation de |p̂|√r ±{:.1} %",
            0.9 * rmax,
            100.0 * (hi - lo) / (hi + lo)
        );
        println!("| r/D | direct | FW-H 1,5 D |");
        println!("|---|---|---|");
        for l in line.iter().filter(|l| [2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 75.0].contains(&l.0)) {
            println!("| {:.0} | {:.3} | {:.3} |", l.0, 1e5 * l.1 * l.0.sqrt(), 1e5 * l.2.abs() * l.0.sqrt());
        }
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let r = if let Some(i) = a.iter().position(|s| s == "--relire") {
        load(Path::new(&a[i + 1]))
    } else {
        let out: PathBuf = arg(a.as_slice(), "--sortie", String::from("acoustique2d")).into();
        simulate(&a, &out)
    };
    analyze(&r);
}
