//! Polaire automatique de la soufflerie 2D (portage du balayage de js/app.js) :
//! pour chaque incidence (ou rotation du cylindre), stabilisation puis mesure,
//! avec des durées adaptatives : un point stationnaire convergé n'attend pas
//! la durée maximale, un point périodique attend trois périodes.

use lbm::stats::{analyze, Analysis};

use crate::worker::Frame;

/// Durées maximales et minimales (temps convectifs).
const SETTLE_T: f64 = 12.0;
const MEASURE_T: f64 = 15.0;
const SETTLE_MIN_T: f64 = 5.0;
const MEASURE_MIN_T: f64 = 6.0;
/// Nombre maximal de points.
const MAX_POINTS: usize = 60;

/// Variable balayée.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SweepVar {
    /// Incidence (°).
    Alpha,
    /// Rotation du cylindre ΩR/U.
    Spin,
}

impl SweepVar {
    pub fn name(self, v: f64) -> String {
        match self {
            SweepVar::Alpha => format!("α = {v:.1}°"),
            SweepVar::Spin => format!("ΩR/U = {v:.2}"),
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Phase {
    Settle,
    Measure,
}

/// Point de polaire enregistré.
#[derive(Clone, Debug)]
pub struct PolarPoint {
    pub var: SweepVar,
    pub value: f64,
    pub label: String,
    pub re: f64,
    pub l: usize,
    pub cx: f64,
    pub cz: f64,
    pub cm: f64,
    /// Strouhal, si l'écoulement est périodique.
    pub st: Option<f64>,
    /// Référence « panneaux + couche limite » (Cx, Cz), si disponible.
    pub panel: Option<(f64, f64)>,
    /// Point de balayage ou ajout manuel.
    pub auto: bool,
    /// Traînée LBM exploitable (faux pour un profil attaché sous loi de paroi).
    pub cx_ok: bool,
}

impl PolarPoint {
    /// Traînée retenue : celle des panneaux quand la traînée LBM n'est pas fiable.
    pub fn cx_retained(&self) -> f64 {
        match (self.cx_ok, self.panel) {
            (false, Some((cd, _))) => cd,
            _ => self.cx,
        }
    }
}

/// La traînée LBM d'un profil est-elle exploitable ? Non quand la loi de
/// paroi est active et l'écoulement attaché : le bord d'attaque n'a qu'une
/// maille ou moins, et une erreur de quelques dixièmes de degré sur la
/// direction de la résultante bascule de la portance dans la traînée (NACA
/// 2412, Re = 10⁶, 4° : Cx LBM de 0,016 à 0,022 de L = 32 à 96, sans
/// converger, pour 0,008 attendu). En écoulement décollé, la traînée de
/// pression domine et le calcul LBM reste la meilleure estimation.
pub fn cx_reliable(fr: &Frame) -> bool {
    match &fr.panel {
        Some(r) if fr.sa => lbm::panel::section_coefficients(r).kirchhoff,
        _ => true,
    }
}

struct Running {
    var: SweepVar,
    values: Vec<f64>,
    i: usize,
    phase: Phase,
    t0: f64,
}

/// État du balayage et points enregistrés.
pub struct Sweep {
    pub from: f64,
    pub to: f64,
    pub step: f64,
    run: Option<Running>,
    pub points: Vec<PolarPoint>,
    pub message: String,
}

/// Analyse de l'historique publié depuis `t0` (t*), et dérive : écart entre
/// les moyennes des deux moitiés de la fenêtre (Cx et Cz).
fn analyze_since(fr: &Frame, t0: f64) -> (Analysis, f64) {
    let w: Vec<&[f64; 4]> = fr.history.iter().filter(|h| h[0] >= t0).collect();
    let col = |k: usize| w.iter().map(|h| h[k]).collect::<Vec<_>>();
    let a = analyze(&col(0), &col(1), &col(2), &col(3));
    let n = w.len();
    let drift = if n >= 4 && !a.periodic {
        let mean = |k: usize, r: std::ops::Range<usize>| {
            let len = r.len().max(1) as f64;
            w[r].iter().map(|h| h[k]).sum::<f64>() / len
        };
        let h = n / 2;
        (mean(1, 0..h) - mean(1, h..n)).abs().max((mean(2, 0..h) - mean(2, h..n)).abs())
    } else {
        0.0
    };
    (a, drift)
}

fn converged(a: &Analysis, drift: f64) -> bool {
    !a.periodic && a.n >= 8 && drift < 0.003 + 0.003 * a.cx.abs()
}

impl Sweep {
    pub fn new() -> Sweep {
        Sweep {
            from: -4.0,
            to: 12.0,
            step: 2.0,
            run: None,
            points: Vec::new(),
            message: String::new(),
        }
    }

    pub fn active(&self) -> bool {
        self.run.is_some()
    }

    /// Valeurs à calculer, bornées à la plage de la maquette.
    fn values(&self, lim: (f64, f64)) -> Vec<f64> {
        let step = self.step.abs();
        if !(step > 0.0) || !self.from.is_finite() || !self.to.is_finite() {
            return Vec::new();
        }
        let dir = if self.to >= self.from { 1.0 } else { -1.0 };
        let mut out: Vec<f64> = Vec::new();
        let mut x = self.from;
        while dir * (self.to - x) >= -1e-9 && out.len() < MAX_POINTS {
            let v = (x.clamp(lim.0, lim.1) * 100.0).round() / 100.0;
            if !out.contains(&v) {
                out.push(v);
            }
            x += dir * step;
        }
        out
    }

    /// Démarre le balayage ; renvoie la première valeur à imposer.
    pub fn start(&mut self, var: SweepVar, lim: (f64, f64), t_star: f64) -> Option<f64> {
        let values = self.values(lim);
        let Some(&first) = values.first() else {
            self.message = "Bornes de balayage invalides.".into();
            return None;
        };
        self.run = Some(Running {
            var,
            values,
            i: 0,
            phase: Phase::Settle,
            t0: t_star,
        });
        self.message.clear();
        Some(first)
    }

    pub fn stop(&mut self, msg: &str) {
        if self.run.take().is_some() {
            self.message = msg.into();
        }
    }

    /// À chaque image : avance le balayage. Renvoie la prochaine valeur à
    /// imposer quand un point vient d'être enregistré.
    pub fn tick(&mut self, fr: &Frame, re: f64, l: usize) -> Option<f64> {
        let run = self.run.as_mut()?;
        if fr.diverged {
            self.stop("Polaire interrompue : le calcul a divergé.");
            return None;
        }
        let now = fr.t_star;
        // La maquette vient de changer : le fil de calcul n'a pas encore avancé.
        if now < run.t0 {
            run.t0 = now;
        }
        let el = now - run.t0;
        match run.phase {
            Phase::Settle => {
                let ok = el >= SETTLE_MIN_T && {
                    let (a, d) = analyze_since(fr, now - 4.0);
                    converged(&a, d)
                };
                if el >= SETTLE_T || ok {
                    run.phase = Phase::Measure;
                    run.t0 = now;
                }
                None
            }
            Phase::Measure => {
                let (a, d) = analyze_since(fr, run.t0);
                let enough = (el >= MEASURE_MIN_T && converged(&a, d))
                    || (el >= MEASURE_MIN_T && a.periodic && a.periods >= 3)
                    || (el >= MEASURE_T && (!a.periodic || a.periods >= 3))
                    || el >= 2.5 * MEASURE_T;
                if !(enough && a.cx.is_finite()) {
                    return None;
                }
                let (var, value) = (run.var, run.values[run.i]);
                self.points.push(point(fr, &a, var, value, re, l, true));
                run.i += 1;
                if run.i >= run.values.len() {
                    let n = run.values.len();
                    self.stop(&format!("Polaire terminée : {n} points enregistrés."));
                    return None;
                }
                run.phase = Phase::Settle;
                run.t0 = now;
                Some(run.values[run.i])
            }
        }
    }

    /// Avancement (0 à 1) et description du point en cours.
    pub fn progress(&self, fr: &Frame) -> Option<(f32, String)> {
        let run = self.run.as_ref()?;
        let n = run.values.len();
        let el = (fr.t_star - run.t0).max(0.0);
        let frac = match run.phase {
            Phase::Settle => 0.45 * (el / SETTLE_T).min(1.0),
            Phase::Measure => 0.45 + 0.55 * (el / MEASURE_T).min(1.0),
        };
        Some((
            ((run.i as f64 + frac) / n as f64) as f32,
            format!(
                "Point {}/{n} · {} · {}",
                run.i + 1,
                run.var.name(run.values[run.i]),
                if run.phase == Phase::Settle { "stabilisation" } else { "mesure" }
            ),
        ))
    }

    /// Ajoute le point courant (moyenne établie) à la main.
    pub fn add_current(&mut self, fr: &Frame, var: SweepVar, value: f64, re: f64, l: usize) {
        self.points.push(point(fr, &fr.analysis, var, value, re, l, false));
    }

    /// Points en CSV.
    pub fn csv(&self) -> String {
        let mut s = String::from("maquette,variable,valeur,Re,L,Cx,Cz,Cm,St,Cx_panneaux,Cz_panneaux,mode,Cx_retenu\n");
        for p in &self.points {
            let opt = |v: Option<f64>, d: usize| v.map(|x| format!("{x:.d$}")).unwrap_or_default();
            s += &format!(
                "{},{},{},{},{},{:.5},{:.5},{:.5},{},{},{},{},{:.5}\n",
                p.label,
                match p.var {
                    SweepVar::Alpha => "alpha_deg",
                    SweepVar::Spin => "rotation",
                },
                p.value,
                p.re,
                p.l,
                p.cx,
                p.cz,
                p.cm,
                opt(p.st, 4),
                opt(p.panel.map(|q| q.0), 5),
                opt(p.panel.map(|q| q.1), 4),
                if p.auto { "polaire" } else { "manuel" },
                p.cx_retained()
            );
        }
        s
    }
}

fn point(fr: &Frame, a: &Analysis, var: SweepVar, value: f64, re: f64, l: usize, auto: bool) -> PolarPoint {
    let panel = fr.panel.as_ref().map(|r| {
        let k = lbm::panel::section_coefficients(r);
        (k.cd, k.cl)
    });
    PolarPoint {
        var,
        value,
        label: fr.label.clone(),
        re,
        l,
        cx: a.cx,
        cz: a.cz,
        cm: a.cm,
        st: (a.periodic && a.st.is_finite()).then_some(a.st),
        panel,
        auto,
        cx_ok: cx_reliable(fr),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Historique synthétique : Cx et Cz constants, 16 échantillons par t*.
    fn frame(t_end: f64, cx: f64, cz: f64) -> Frame {
        let history = (0..=(t_end * 16.0) as usize)
            .map(|k| [k as f64 / 16.0, cx, cz, 0.0])
            .collect();
        Frame {
            history,
            t_star: t_end,
            label: "test".into(),
            ..Default::default()
        }
    }

    #[test]
    fn valeurs_bornees() {
        let mut s = Sweep::new();
        (s.from, s.to, s.step) = (-4.0, 30.0, 10.0);
        assert_eq!(s.values((-20.0, 25.0)), vec![-4.0, 6.0, 16.0, 25.0]);
        (s.from, s.to, s.step) = (2.0, -2.0, 1.0);
        assert_eq!(s.values((-20.0, 25.0)), vec![2.0, 1.0, 0.0, -1.0, -2.0]);
    }

    #[test]
    fn point_stationnaire_converge_vite() {
        let mut s = Sweep::new();
        (s.from, s.to, s.step) = (0.0, 2.0, 2.0);
        assert_eq!(s.start(SweepVar::Alpha, (-20.0, 25.0), 10.0), Some(0.0));
        // Stabilisation : écoulement établi dès 5 t*.
        assert_eq!(s.tick(&frame(15.5, 1.2, 0.3), 1000.0, 32), None);
        assert!(s.progress(&frame(15.5, 1.2, 0.3)).unwrap().1.contains("mesure"));
        // Mesure : 6 t* suffisent, point enregistré, valeur suivante demandée.
        assert_eq!(s.tick(&frame(22.0, 1.2, 0.3), 1000.0, 32), Some(2.0));
        assert_eq!(s.points.len(), 1);
        assert!((s.points[0].cx - 1.2).abs() < 1e-9 && s.points[0].value == 0.0);
        // Second et dernier point.
        s.tick(&frame(28.0, 1.1, 0.5), 1000.0, 32);
        s.tick(&frame(35.0, 1.1, 0.5), 1000.0, 32);
        assert!(!s.active() && s.points.len() == 2);
        assert!(s.csv().lines().count() == 3);
    }
}
