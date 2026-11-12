//! Moyennes sur un nombre entier de périodes, amplitude, nombre de Strouhal.

#[derive(Clone, Copy, Debug, Default)]
pub struct Analysis {
    pub n: usize,
    pub cx: f64,
    pub cz: f64,
    pub cm: f64,
    /// Demi-amplitudes crête à crête.
    pub cx_amp: f64,
    pub cz_amp: f64,
    /// Strouhal f·L/U (t en temps convectifs).
    pub st: f64,
    pub periodic: bool,
    pub periods: usize,
}

fn mean_std(v: &[f64]) -> (f64, f64) {
    let n = v.len().max(1) as f64;
    let m = v.iter().sum::<f64>() / n;
    let s = (v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / n).sqrt();
    (m, s)
}

/// `t` en temps convectifs t* ; séries de même longueur.
pub fn analyze(t: &[f64], cx: &[f64], cz: &[f64], cm: &[f64]) -> Analysis {
    let n = t.len();
    let mut a = Analysis {
        n,
        cx: f64::NAN,
        cz: f64::NAN,
        cm: f64::NAN,
        st: f64::NAN,
        ..Default::default()
    };
    if n < 8 {
        return a;
    }
    // Oscillation sur Cz, à défaut sur Cx, détectée avec hystérésis.
    let (mz, sz) = mean_std(cz);
    let (mx, sx) = mean_std(cx);
    let (sig, m, s) = if sz < 1e-4 && sx > 1e-4 {
        (cx, mx, sx)
    } else {
        (cz, mz, sz)
    };
    let mut ups = Vec::new();
    if s > 1e-4 {
        let h = 0.25 * s;
        let mut state = 0i8;
        for k in 0..n {
            let v = sig[k] - m;
            if v < -h {
                state = -1;
            } else if v > h && state == -1 {
                state = 1;
                // Passage montant interpolé entre k−1 et k.
                let (v0, v1) = (sig[k.max(1) - 1] - m, v);
                let f = if v1 != v0 {
                    (-v0 / (v1 - v0)).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                let tk = if k > 0 { t[k - 1] + f * (t[k] - t[k - 1]) } else { t[k] };
                ups.push((k, tk));
            }
        }
    }
    let (i0, i1) = if ups.len() >= 3 {
        a.periodic = true;
        a.periods = ups.len() - 1;
        let (k0, t0) = ups[0];
        let (k1, t1) = ups[ups.len() - 1];
        a.st = a.periods as f64 / (t1 - t0);
        (k0, k1)
    } else {
        (0, n)
    };
    let span = |v: &[f64]| -> (f64, f64) {
        let w = &v[i0..i1.max(i0 + 1).min(n)];
        let m = w.iter().sum::<f64>() / w.len() as f64;
        let lo = w.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = w.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        (m, 0.5 * (hi - lo))
    };
    let (cxm, cxa) = span(cx);
    let (czm, cza) = span(cz);
    let (cmm, _) = span(cm);
    a.cx = cxm;
    a.cz = czm;
    a.cm = cmm;
    a.cx_amp = cxa;
    a.cz_amp = cza;
    // Oscillation d'amplitude négligeable (bruit de bord de fuite, reste du
    // transitoire) : écoulement considéré comme établi, sans Strouhal.
    if a.periodic && cza < 0.01 && cxa < 0.005 {
        a.periodic = false;
        a.periods = 0;
        a.st = f64::NAN;
    }
    a
}
