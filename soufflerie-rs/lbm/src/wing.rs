//! Aile finie : surface portante de Weissinger (un tourbillon en fer à cheval
//! par bande d'envergure, segment lié au quart de corde, condition au point
//! des trois quarts de corde), rendue non linéaire par « décambrage » : la
//! circulation de chaque bande vérifie la loi de Kutta–Joukowski avec la
//! polaire 2D de sa section, à l'incidence effective α(3/4 c) + Cl/a0. En
//! régime linéaire, on retrouve exactement la méthode de Weissinger ; au-delà,
//! la polaire (décrochage compris) fixe la charge. Flèche, dièdre, effilement
//! et vrillage sont pris en compte. Système non linéaire résolu par Newton.
//!
//! Portage de `js/vlm.js` et du calcul de `js/wing-tab.js`.
//!
//! Repère : x vers l'aval, y vers la droite (envergure), z vers le haut.
//! Corde moyenne géométrique c̄ = 1, donc envergure b = allongement et S = b.

use std::f64::consts::PI;

use rayon::prelude::*;

use crate::panel::{panel_solve, section_coefficients, PanelOptions, PanelShape};
use crate::shapes::Naca;

type V3 = [f64; 3];

fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}
fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn unit(a: V3) -> V3 {
    let n = norm(a);
    scale(a, 1.0 / if n > 0.0 { n } else { 1.0 })
}
fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Géométrie en plan trapézoïdale.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WingParams {
    /// Allongement b²/S.
    pub ar: f64,
    /// Effilement c_bout / c_emplanture.
    pub taper: f64,
    /// Flèche du quart de corde (°).
    pub sweep: f64,
    /// Dièdre (°).
    pub dihedral: f64,
    /// Vrillage en bout (°, linéaire ; négatif = « wash-out »).
    pub twist: f64,
    /// Bandes par demi-aile.
    pub n: usize,
}

impl Default for WingParams {
    fn default() -> Self {
        WingParams {
            ar: 8.0,
            taper: 0.5,
            sweep: 0.0,
            dihedral: 0.0,
            twist: 0.0,
            n: 40,
        }
    }
}

/// Bande d'envergure : fer à cheval A → B au quart de corde.
#[derive(Clone, Debug)]
pub struct Strip {
    pub a: V3,
    pub b: V3,
    /// Milieu du segment lié.
    pub p: V3,
    /// Point de contrôle aux trois quarts de corde.
    pub p34: V3,
    pub dl: V3,
    /// Axe de corde (vers l'aval), normale, direction du segment.
    pub ua: V3,
    pub un: V3,
    pub us: V3,
    pub c: f64,
    /// Largeur projetée dans le plan (y, z).
    pub ds: f64,
    /// Position relative en envergure, de −1 à 1.
    pub eta: f64,
    /// Vrillage local (°).
    pub twist: f64,
}

#[derive(Clone, Debug)]
pub struct Wing {
    pub params: WingParams,
    pub b: f64,
    pub s: f64,
    pub cr: f64,
    pub ct: f64,
    /// Corde aérodynamique moyenne.
    pub mac: f64,
    /// Point de référence des moments (quart de la corde moyenne).
    pub reference: V3,
    pub strips: Vec<Strip>,
}

impl Wing {
    pub fn new(o: WingParams) -> Wing {
        let n = o.n.max(1);
        let b = o.ar;
        let half = b / 2.0;
        let cr = 2.0 / (1.0 + o.taper);
        let tan_s = o.sweep.to_radians().tan();
        let tan_d = o.dihedral.to_radians().tan();
        let twist_tip = o.twist.to_radians();
        let chord = |y: f64| cr * (1.0 - (1.0 - o.taper) * y.abs() / half);
        let c4 = |y: f64| [y.abs() * tan_s, y, y.abs() * tan_d];
        // Bandes uniformes : une bande très étroite devant la corde (la
        // répartition en cosinus en crée en bout d'aile) rend la condition aux
        // trois quarts de corde mal posée et fausse la traînée induite.
        let ys: Vec<f64> = (0..=2 * n).map(|k| half * (2.0 * k as f64 / (2 * n) as f64 - 1.0)).collect();
        let strips = (0..2 * n)
            .map(|i| {
                let (ya, yb) = (ys[i], ys[i + 1]);
                let ym = 0.5 * (ya + yb);
                let a = c4(ya);
                let bb = c4(yb);
                let p = scale(add(a, bb), 0.5);
                let dl = sub(bb, a);
                let us = unit(dl);
                let ua0 = [1.0, 0.0, 0.0];
                let un0 = unit(cross(ua0, us));
                let th = twist_tip * ym.abs() / half;
                let (sn, cs) = th.sin_cos();
                let ua = sub(scale(ua0, cs), scale(un0, sn));
                let un = add(scale(ua0, sn), scale(un0, cs));
                let c = chord(ym);
                let p34 = add(p, scale(ua, 0.5 * c));
                Strip {
                    a,
                    b: bb,
                    p,
                    p34,
                    dl,
                    ua,
                    un,
                    us,
                    c,
                    ds: dl[1].hypot(dl[2]),
                    eta: ym / half,
                    twist: th.to_degrees(),
                }
            })
            .collect();
        let t = o.taper;
        let mac = (2.0 / 3.0) * cr * (1.0 + t + t * t) / (1.0 + t);
        let y_mac = (b / 6.0) * (1.0 + 2.0 * t) / (1.0 + t);
        Wing {
            params: o,
            b,
            s: b,
            cr,
            ct: cr * t,
            mac,
            reference: [y_mac * tan_s, 0.0, 0.0],
            strips,
        }
    }

    /// Corde locale en y.
    pub fn chord_at(&self, y: f64) -> f64 {
        self.cr * (1.0 - (1.0 - self.params.taper) * y.abs() / (self.b / 2.0))
    }

    /// Quart de corde en y.
    pub fn c4_at(&self, y: f64) -> V3 {
        let p = &self.params;
        [y.abs() * p.sweep.to_radians().tan(), y, y.abs() * p.dihedral.to_radians().tan()]
    }
}

/// Coefficients d'une section 2D.
#[derive(Clone, Copy, Debug, Default)]
pub struct SectionCoef {
    pub cl: f64,
    pub cd: f64,
    pub cm: f64,
}

/// Vitesse induite en P par un segment tourbillonnaire A → B de circulation 1.
fn segment(p: V3, a: V3, b: V3) -> V3 {
    let r1 = sub(p, a);
    let r2 = sub(p, b);
    let c = cross(r1, r2);
    let c2 = dot(c, c);
    let n1 = norm(r1);
    let n2 = norm(r2);
    if c2 < 1e-10 || n1 < 1e-9 || n2 < 1e-9 {
        return [0.0; 3];
    }
    let r0 = sub(b, a);
    let k = (dot(r0, r1) / n1 - dot(r0, r2) / n2) / (4.0 * PI * c2);
    scale(c, k)
}

/// Tourbillon semi-infini partant de Q dans la direction d (unitaire), circulation 1.
fn semi_infinite(p: V3, q: V3, d: V3) -> V3 {
    let r = sub(p, q);
    let c = cross(d, r);
    let c2 = dot(c, c);
    let nr = norm(r);
    if c2 < 1e-10 || nr < 1e-9 {
        return [0.0; 3];
    }
    let k = (1.0 + dot(r, d) / nr) / (4.0 * PI * c2);
    scale(c, k)
}

/// Influence d'un fer à cheval (tourbillons libres selon d) en P.
fn horseshoe(p: V3, s: &Strip, d: V3) -> V3 {
    let v1 = segment(p, s.a, s.b);
    let v2 = semi_infinite(p, s.b, d);
    let v3 = semi_infinite(p, s.a, d);
    [v1[0] + v2[0] - v3[0], v1[1] + v2[1] - v3[1], v1[2] + v2[2] - v3[2]]
}

/// Résultat d'une bande.
#[derive(Clone, Copy, Debug)]
pub struct StripResult {
    pub eta: f64,
    pub c: f64,
    pub cl: f64,
    pub cd: f64,
    /// Incidence effective (°).
    pub alpha: f64,
    /// Charge locale c·Cl/c̄ (c̄ = 1), au sens de Kutta–Joukowski.
    pub ccl: f64,
}

#[derive(Clone, Debug)]
pub struct WingResult {
    pub alpha: f64,
    pub cl: f64,
    /// Traînée induite (plan de Trefftz).
    pub cdi: f64,
    /// Traînée de profil (polaires de section).
    pub cdv: f64,
    pub cd: f64,
    pub cm: f64,
    /// Coefficient d'Oswald CL²/(π A CDi).
    pub e: f64,
    pub strips: Vec<StripResult>,
    pub gamma: Vec<f64>,
    pub converged: bool,
    pub residual: f64,
    pub iterations: usize,
}

/// Résout l'aile à une incidence `alpha_deg` (emplanture). `section(i, α)`
/// donne la polaire 2D de la bande i ; `gamma0` amorce Newton (continuation).
pub fn solve_wing(
    wing: &Wing,
    alpha_deg: f64,
    section: &(dyn Fn(usize, f64) -> SectionCoef + Sync),
    gamma0: Option<&[f64]>,
) -> WingResult {
    let strips = &wing.strips;
    let n = strips.len();
    let a = alpha_deg.to_radians();
    let vinf = [a.cos(), 0.0, a.sin()];
    // Sillage plan, dans le prolongement de la corde (hypothèse classique de la
    // surface portante linéarisée).
    let d = [1.0, 0.0, 0.0];
    // Influences : au milieu du segment lié (Kutta–Joukowski) et aux trois
    // quarts de corde (incidence locale, comme Weissinger).
    let rows: Vec<(Vec<V3>, Vec<V3>)> = strips
        .par_iter()
        .map(|si| {
            let vb = strips.iter().map(|sj| horseshoe(si.p, sj, d)).collect();
            let v3 = strips.iter().map(|sj| horseshoe(si.p34, sj, d)).collect();
            (vb, v3)
        })
        .collect();
    let (vb_m, v3_m): (Vec<Vec<V3>>, Vec<Vec<V3>>) = rows.into_iter().unzip();
    // Partie linéaire de chaque polaire : pente a0 (par radian).
    let a0: Vec<f64> = (0..n)
        .map(|i| {
            let c1 = section(i, -2.0).cl;
            let c2 = section(i, 4.0).cl;
            ((c2 - c1) / 6.0).max(0.02).to_degrees()
        })
        .collect();
    let mut g: Vec<f64> = match gamma0 {
        Some(g) if g.len() == n => g.to_vec(),
        _ => vec![0.0; n],
    };
    // Résidu d'une bande : Cl de Kutta–Joukowski moins Cl de la polaire à
    // l'incidence effective α_e = α(3/4 c) + Cl/a0 (« décambrage »).
    let local = |i: usize, wbv: V3, w3v: V3, gi: f64| -> (f64, f64) {
        let s = &strips[i];
        let ck = 2.0 * gi * norm(cross(wbv, s.dl)) / (s.c * s.ds);
        let a34 = dot(w3v, s.un).atan2(dot(w3v, s.ua));
        let ae = a34 + ck / a0[i];
        (ae, ck - section(i, ae.to_degrees()).cl)
    };
    let acc = |m: &Vec<Vec<V3>>, i: usize, gv: &[f64]| -> V3 {
        let mut v = vinf;
        for (vj, &gj) in m[i].iter().zip(gv) {
            v[0] += vj[0] * gj;
            v[1] += vj[1] * gj;
            v[2] += vj[2] * gj;
        }
        v
    };
    let mut wb = vec![[0.0; 3]; n];
    let mut w3 = vec![[0.0; 3]; n];
    let mut aeff = vec![0.0; n];
    let mut res_v = vec![0.0; n];
    let evaluate = |gv: &[f64], wb: &mut [V3], w3: &mut [V3], aeff: &mut [f64], rv: &mut [f64]| -> f64 {
        let mut r2 = 0.0;
        for i in 0..n {
            wb[i] = acc(&vb_m, i, gv);
            w3[i] = acc(&v3_m, i, gv);
            let (ae, r) = local(i, wb[i], w3[i], gv[i]);
            aeff[i] = ae;
            rv[i] = r;
            r2 += r * r;
        }
        (r2 / n as f64).sqrt()
    };
    let mut res = evaluate(&g, &mut wb, &mut w3, &mut aeff, &mut res_v);
    let mut iter = 0;
    let max_iter = 40;
    let h = 1e-6;
    while res > 1e-9 && iter < max_iter {
        iter += 1;
        // Jacobien par différences finies (la ligne i ne dépend que des
        // vitesses induites en i, perturbées colonne par colonne).
        let jac: Vec<Vec<f64>> = (0..n)
            .into_par_iter()
            .map(|i| {
                (0..n)
                    .map(|j| {
                        let vb = vb_m[i][j];
                        let v3 = v3_m[i][j];
                        let wbp = add(wb[i], scale(vb, h));
                        let w3p = add(w3[i], scale(v3, h));
                        let gi = g[i] + if i == j { h } else { 0.0 };
                        (local(i, wbp, w3p, gi).1 - res_v[i]) / h
                    })
                    .collect()
            })
            .collect();
        let rhs: Vec<f64> = res_v.iter().map(|r| -r).collect();
        let Some(dg) = solve_linear(jac, rhs) else { break };
        // Newton amorti : on réduit le pas tant que le résidu ne baisse pas.
        let g0 = g.clone();
        let mut step = 1.0;
        for t in 0..8 {
            for i in 0..n {
                g[i] = g0[i] + step * dg[i];
            }
            let r = evaluate(&g, &mut wb, &mut w3, &mut aeff, &mut res_v);
            if r < res || t == 7 {
                res = r;
                break;
            }
            step *= 0.5;
        }
    }
    // Efforts : Kutta–Joukowski avec la vitesse locale, traînée de profil des
    // polaires de section.
    let el = [-a.sin(), 0.0, a.cos()];
    let (mut lift, mut dv, mut my) = (0.0, 0.0, 0.0);
    let mut per = Vec::with_capacity(n);
    for i in 0..n {
        let s = &strips[i];
        let f = scale(cross(wb[i], s.dl), g[i]);
        let sec = section(i, aeff[i].to_degrees());
        lift += dot(f, el);
        dv += 0.5 * sec.cd * s.c * s.ds;
        let r = sub(s.p, wing.reference);
        my += r[2] * f[0] - r[0] * f[2] + 0.5 * sec.cm * s.c * s.c * s.ds * s.us[1];
        per.push(StripResult {
            eta: s.eta,
            c: s.c,
            cl: sec.cl,
            cd: sec.cd,
            alpha: aeff[i].to_degrees(),
            ccl: 2.0 * dot(f, el) / s.ds,
        });
    }
    // Traînée induite dans le plan de Trefftz (plus robuste que
    // Kutta–Joukowski pour les ailes en flèche).
    let di = trefftz_drag(strips, &g);
    let q = 0.5;
    let cl = lift / (q * wing.s);
    let cdi = di / (q * wing.s);
    let cdv = dv / (q * wing.s);
    WingResult {
        alpha: alpha_deg,
        cl,
        cdi,
        cdv,
        cd: cdi + cdv,
        cm: my / (q * wing.s * wing.mac),
        e: if cdi > 1e-9 {
            cl * cl / (PI * wing.params.ar * cdi)
        } else {
            f64::NAN
        },
        strips: per,
        gamma: g,
        converged: res < 1e-6,
        residual: res,
        iterations: iter,
    }
}

/// D_i = −½ ρ Σ Γ_i w_n,i ds_i, w_n vitesse normale à la nappe (ρ = V∞ = 1).
fn trefftz_drag(strips: &[Strip], g: &[f64]) -> f64 {
    let n = strips.len();
    // Nœuds (y, z) et tourbillons libres (vorticité selon +x).
    let nodes: Vec<(f64, f64, f64)> = (0..=n)
        .map(|k| {
            let p = if k < n { strips[k].a } else { strips[n - 1].b };
            let gam = (if k > 0 { g[k - 1] } else { 0.0 }) - (if k < n { g[k] } else { 0.0 });
            (p[1], p[2], gam)
        })
        .collect();
    let mut d = 0.0;
    for (i, s) in strips.iter().enumerate() {
        let (ym, zm) = (s.p[1], s.p[2]);
        let (mut vy, mut vz) = (0.0, 0.0);
        for &(ny, nz, gam) in &nodes {
            let ry = ym - ny;
            let rz = zm - nz;
            let r2 = ry * ry + rz * rz;
            if r2 < 1e-14 {
                continue;
            }
            let f = gam / (2.0 * PI * r2);
            vy -= f * rz;
            vz += f * ry;
        }
        let sy = s.dl[1] / s.ds;
        let sz = s.dl[2] / s.ds;
        let wn = -sz * vy + sy * vz;
        d -= 0.5 * g[i] * wn * s.ds;
    }
    d
}

/// Élimination de Gauss avec pivot partiel ; None si singulier.
fn solve_linear(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for c in 0..n {
        let p = (c..n).max_by(|&x, &y| a[x][c].abs().total_cmp(&a[y][c].abs()))?;
        if a[p][c].abs() < 1e-14 {
            return None;
        }
        a.swap(c, p);
        b.swap(c, p);
        let (top, rest) = a.split_at_mut(c + 1);
        let pivot = &top[c];
        for (k, row) in rest.iter_mut().enumerate() {
            let f = row[c] / pivot[c];
            if f == 0.0 {
                continue;
            }
            for (x, &y) in row[c..].iter_mut().zip(&pivot[c..]) {
                *x -= f * y;
            }
            b[c + 1 + k] -= f * b[c];
        }
    }
    let mut x = vec![0.0; n];
    for r in (0..n).rev() {
        let s: f64 = b[r] - (r + 1..n).map(|k| a[r][k] * x[k]).sum::<f64>();
        x[r] = s / a[r][r];
    }
    Some(x)
}

/// Polaire 2D tabulée aux incidences entières.
#[derive(Clone, Debug)]
pub struct PolarTable {
    pub lo: i32,
    pub hi: i32,
    pub cl: Vec<f64>,
    pub cd: Vec<f64>,
    pub cm: Vec<f64>,
    pub clmax: f64,
}

/// Bornes des tables de section (°).
pub const TABLE_LO: i32 = -14;
pub const TABLE_HI: i32 = 24;

impl PolarTable {
    /// Table à partir d'une fonction ; une valeur absente (None) prend celle
    /// de l'incidence précédente.
    pub fn new(lo: i32, hi: i32, f: impl Fn(f64) -> Option<SectionCoef> + Sync) -> PolarTable {
        let rows: Vec<Option<SectionCoef>> = (lo..=hi).into_par_iter().map(|a| f(a as f64)).collect();
        let mut last = SectionCoef {
            cl: 0.0,
            cd: 0.02,
            cm: 0.0,
        };
        let (mut cl, mut cd, mut cm) = (Vec::new(), Vec::new(), Vec::new());
        for r in rows {
            let r = r.unwrap_or(last);
            last = r;
            cl.push(r.cl);
            cd.push(r.cd);
            cm.push(r.cm);
        }
        let clmax = cl.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        PolarTable {
            lo,
            hi,
            cl,
            cd,
            cm,
            clmax,
        }
    }

    /// Polaire d'un profil NACA : méthode des panneaux, couche limite intégrale
    /// et décrochage de Kirchhoff, de −14° à 24°.
    pub fn naca(naca: Naca, re: f64, tripped: bool) -> PolarTable {
        let shape = PanelShape::Naca(naca);
        let opts = if tripped {
            PanelOptions::tripped(re, 0.0)
        } else {
            PanelOptions::free(re)
        };
        PolarTable::new(TABLE_LO, TABLE_HI, |a| {
            let r = panel_solve(&shape, a, &opts);
            if !(r.cl.is_finite() && r.cd.is_finite()) {
                return None;
            }
            let k = section_coefficients(&r);
            Some(SectionCoef {
                cl: k.cl,
                cd: k.cd,
                cm: k.cm,
            })
        })
    }

    /// Valeur interpolée ; au-delà, raccord progressif vers la plaque plane
    /// décollée (Cn ≈ 2 sin α).
    pub fn at(&self, a_deg: f64) -> SectionCoef {
        let (lo, hi) = (self.lo as f64, self.hi as f64);
        if a_deg >= lo && a_deg <= hi {
            let k = ((a_deg - lo).floor() as usize).min(self.cl.len() - 2);
            let f = a_deg - lo - k as f64;
            let lerp = |v: &[f64]| v[k] + f * (v[k + 1] - v[k]);
            return SectionCoef {
                cl: lerp(&self.cl),
                cd: lerp(&self.cd),
                cm: lerp(&self.cm),
            };
        }
        let end = if a_deg > hi { self.cl.len() - 1 } else { 0 };
        let a0 = self.lo as f64 + end as f64;
        let r = a_deg.to_radians();
        let cn = 1.98 * r.sin();
        let blend = ((a_deg - a0).abs() / 15.0).min(1.0);
        SectionCoef {
            cl: (1.0 - blend) * self.cl[end] + blend * cn * r.cos(),
            cd: (1.0 - blend) * self.cd[end] + blend * (cn * r.sin() + 0.02),
            cm: (1.0 - blend) * self.cm[end] + blend * (-0.25 * cn),
        }
    }
}

/// Polaires de section le long de l'envergure : interpolation en log Re entre
/// l'emplanture et le bout d'aile.
pub struct Sections {
    pub root: PolarTable,
    pub tip: Option<PolarTable>,
    /// Poids de la table d'emplanture, par bande.
    pub weight: Vec<f64>,
    /// Cl max local, par bande.
    pub clmax: Vec<f64>,
}

impl Sections {
    /// Reynolds des polaires nécessaires (rapportés à la corde moyenne
    /// géométrique) : emplanture, et bout d'aile si l'aile est effilée.
    pub fn reynolds(wing: &Wing, re: f64) -> (f64, Option<f64>) {
        let tip = (wing.ct < wing.cr * 0.98).then_some(re * wing.ct);
        (re * wing.cr, tip)
    }

    pub fn new(wing: &Wing, re: f64, root: PolarTable, tip: Option<PolarTable>) -> Sections {
        let (rr, rt) = Sections::reynolds(wing, re);
        let weight: Vec<f64> = match rt {
            Some(rt) if tip.is_some() => {
                let (lr, lt) = (rr.ln(), rt.ln());
                wing.strips
                    .iter()
                    .map(|s| (((re * s.c).ln() - lt) / (lr - lt)).clamp(0.0, 1.0))
                    .collect()
            }
            _ => vec![1.0; wing.strips.len()],
        };
        let clmax = weight
            .iter()
            .map(|&w| w * root.clmax + (1.0 - w) * tip.as_ref().map_or(root.clmax, |t| t.clmax))
            .collect();
        Sections {
            root,
            tip,
            weight,
            clmax,
        }
    }

    pub fn at(&self, i: usize, a_deg: f64) -> SectionCoef {
        let p = self.root.at(a_deg);
        let w = self.weight[i];
        match &self.tip {
            Some(t) if w < 1.0 => {
                let q = t.at(a_deg);
                SectionCoef {
                    cl: w * p.cl + (1.0 - w) * q.cl,
                    cd: w * p.cd + (1.0 - w) * q.cd,
                    cm: w * p.cm + (1.0 - w) * q.cm,
                }
            }
            _ => p,
        }
    }
}

/// Polaire de l'aile par continuation en incidence (chaque point amorce le suivant).
pub fn wing_polar(wing: &Wing, sections: &Sections, alphas: &[f64]) -> Vec<WingResult> {
    let section = |i: usize, a: f64| sections.at(i, a);
    let mut g: Option<Vec<f64>> = None;
    alphas
        .iter()
        .map(|&a| {
            let r = solve_wing(wing, a, &section, g.as_deref());
            g = Some(r.gamma.clone());
            r
        })
        .collect()
}
