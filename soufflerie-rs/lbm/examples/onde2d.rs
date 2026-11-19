//! Atténuation et dispersion numériques d'une onde acoustique plane par le
//! D2Q9 (collisions du solveur, `d2q9::collide`), sur une ligne périodique.
//!
//! onde2d [nu] [mach]
//!
//! Onde progressive ρ′ = ε cos(kx), u′ = c ρ′ (vers +x), dans un écoulement
//! U = M·c de même sens. Théorie (Navier–Stokes isotherme du LBM, viscosité
//! de volume 2ν/D) : amplitude ∝ exp(−ν k² t), célérité c + U.
use lbm::d2q9::{collide, equilibrium, Relax, CX};
use std::f64::consts::PI;

fn run<const RR: bool>(n: usize, nu: f64, u0: f64, periods: f64) -> (f64, f64) {
    let c = 1.0 / 3f64.sqrt();
    let k = 2.0 * PI / n as f64;
    let eps = 1e-3;
    let tau = (3.0 * nu + 0.5) as f32;
    let mut f: Vec<[f32; 9]> = (0..n)
        .map(|x| {
            let r = eps * (k * x as f64).cos();
            equilibrium((1.0 + r) as f32, (u0 + c * r) as f32, 0.0)
        })
        .collect();
    let q = equilibrium(1.0, u0 as f32, 0.0);
    let relax = Relax { tau, c18: 0.0, gx: 0.0, sigma: 0.0 };
    let omega = k * (c + u0);
    // Assez de pas pour une décroissance théorique d'au moins e^{-1/2}.
    let steps = ((periods * 2.0 * PI / omega).max(0.5 / (nu * k * k))).round() as usize;
    // Composante progressive R̂ = ½ (ρ̂ + û/c), relevée à chaque pas.
    let mut series = Vec::with_capacity(steps);
    let mut g = f.clone();
    for t in 0..steps {
        for x in 0..n {
            f[x] = collide::<RR, false>(&f[x], relax, &q);
        }
        for x in 0..n {
            for i in 0..9 {
                let sx = (x as i64 - CX[i] as i64).rem_euclid(n as i64) as usize;
                g[x][i] = f[sx][i];
            }
        }
        std::mem::swap(&mut f, &mut g);
        let (mut re, mut im) = (0.0, 0.0);
        for x in 0..n {
            let (rho, ux, _) = lbm::d2q9::moments(&f[x]);
            let v = 0.5 * ((rho as f64 - 1.0) + (ux as f64 - u0) / c);
            re += v * (k * x as f64).cos();
            im -= v * (k * x as f64).sin();
        }
        series.push((t + 1, re.hypot(im), im.atan2(re)));
    }
    // Ajustements sur la seconde moitié : log |R̂| et phase déroulée.
    let half = &series[steps / 2..];
    let fit = |ys: &[f64]| {
        let m = ys.len() as f64;
        let xs: Vec<f64> = half.iter().map(|s| s.0 as f64).collect();
        let (sx, sy) = (xs.iter().sum::<f64>() / m, ys.iter().sum::<f64>() / m);
        let sxy: f64 = xs.iter().zip(ys).map(|(x, y)| (x - sx) * (y - sy)).sum();
        let sxx: f64 = xs.iter().map(|x| (x - sx).powi(2)).sum();
        sxy / sxx
    };
    let logs: Vec<f64> = half.iter().map(|s| s.1.ln()).collect();
    let mut ph = Vec::with_capacity(half.len());
    let mut acc = 0.0;
    let mut prev = half[0].2;
    for s in half {
        let mut d = s.2 - prev;
        while d > PI {
            d -= 2.0 * PI;
        }
        while d < -PI {
            d += 2.0 * PI;
        }
        acc += d;
        prev = s.2;
        ph.push(acc);
    }
    let gamma = -fit(&logs);
    let om = -fit(&ph);
    (gamma / (nu * k * k), om / omega)
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let nu: f64 = a.get(1).and_then(|s| s.parse().ok()).unwrap_or(0.0123);
    let mach: f64 = a.get(2).and_then(|s| s.parse().ok()).unwrap_or(0.2);
    let c = 1.0 / 3f64.sqrt();
    println!("ν = {nu} (τ = {:.4}), Mach de l'écoulement {mach}", 3.0 * nu + 0.5);
    println!("points/λ | BGK : γ/(νk²)  ω/(k(c+U)) | RR : γ/(νk²)  ω/(k(c+U)) | perte RR par λ parcourue (dB)");
    for &n in &[4usize, 6, 8, 12, 16, 24, 32, 64, 128, 256] {
        let periods = 30.0;
        let (gb, wb) = run::<false>(n, nu, mach * c, periods);
        let (gr, wr) = run::<true>(n, nu, mach * c, periods);
        // Perte par longueur d'onde parcourue : γ · λ/(c+U), en dB.
        let k = 2.0 * PI / n as f64;
        let db = 20.0 / 10f64.ln() * gr * nu * k * k * n as f64 / (c * (1.0 + mach));
        println!("{n:8} | {gb:12.3} {wb:12.5} | {gr:11.3} {wr:12.5} | {db:.3e}");
    }
}
