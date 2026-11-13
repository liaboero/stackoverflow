//! Onglet « Aile finie » : surface portante non linéaire alimentée par les
//! polaires 2D de la méthode des panneaux (lbm::wing, portage de
//! js/wing-tab.js). Les polaires de section et la polaire de l'aile sont
//! calculées dans un fil séparé ; l'incidence courante est résolue à chaque
//! changement (quelques millisecondes).

use std::collections::HashMap;
use std::sync::{Arc, Mutex};

use eframe::egui::{self, Color32, Pos2, Sense, Shape, Stroke, Vec2};
use egui_plot::{Legend, Line, LineStyle, Plot, PlotPoints, Points, VLine};
use lbm::wing::{solve_wing, wing_polar, PolarTable, Sections, Wing, WingParams, WingResult};
use lbm::Naca;

/// Reynolds proposés (rapportés à la corde moyenne géométrique).
const W_RE: [f64; 8] = [5e4, 1e5, 2e5, 5e5, 1e6, 2e6, 5e6, 1e7];
const POLAR_LO: i32 = -6;
const POLAR_HI: i32 = 20;
const N_HALF: usize = 40;
/// Viscosité cinématique de l'air à 15 °C (m²/s).
const NU_AIR: f64 = 1.46e-5;
const RHO_AIR: f64 = 1.225;

const ACCENT: Color32 = Color32::from_rgb(0xc8, 0x46, 0x1c);
const LIFT: Color32 = Color32::from_rgb(0x1f, 0x6f, 0x9c);
const WARN: Color32 = Color32::from_rgb(0x93, 0x65, 0x0f);
const INK: Color32 = Color32::from_rgb(0x13, 0x20, 0x2a);

/// Rampe de couleur Cl/Cl max : froid → chaud (rouge = décrochage local).
const RAMP: [(f32, [u8; 3]); 4] = [
    (0.0, [64, 112, 170]),
    (0.6, [96, 170, 160]),
    (0.85, [226, 190, 76]),
    (1.0, [206, 66, 46]),
];

fn ramp(t: f64) -> Color32 {
    let x = t.clamp(0.0, 1.0) as f32;
    let mut k = 0;
    while k < RAMP.len() - 2 && x > RAMP[k + 1].0 {
        k += 1;
    }
    let (x0, c0) = RAMP[k];
    let (x1, c1) = RAMP[k + 1];
    let f = (x - x0) / (x1 - x0);
    let m = |i: usize| (c0[i] as f32 + f * (c1[i] as f32 - c0[i] as f32)).round() as u8;
    Color32::from_rgb(m(0), m(1), m(2))
}

/// Ce qui détermine les polaires et la polaire de l'aile.
#[derive(Clone, PartialEq)]
struct Key {
    naca: String,
    re: f64,
    tripped: bool,
    params: WingParams,
}

/// Résultat du fil de calcul.
struct Computed {
    key: Key,
    wing: Wing,
    sections: Sections,
    polar: Vec<WingResult>,
}

#[derive(Default)]
struct Background {
    gen: u64,
    done: Option<Arc<Computed>>,
    /// Polaires de section déjà calculées : (NACA, Re arrondi, transition imposée).
    tables: HashMap<(String, u64, bool), Arc<PolarTable>>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum PolarMode {
    Alpha,
    Eiffel,
}

pub struct WingTab {
    params: WingParams,
    naca: Naca,
    naca_text: String,
    re_index: usize,
    tripped: bool,
    alpha: f64,
    /// Corde moyenne (m), pour les grandeurs physiques.
    chord: f64,
    polar_mode: PolarMode,
    /// Caméra orbitale : azimut et élévation (°).
    cam: (f32, f32),
    bg: Arc<Mutex<Background>>,
    requested: Option<Key>,
    /// Incidence courante résolue : (calcul source, incidence, résultat).
    current: Option<(Arc<Computed>, f64, WingResult)>,
}

impl WingTab {
    pub fn new() -> WingTab {
        WingTab {
            params: WingParams {
                n: N_HALF,
                ..Default::default()
            },
            naca: Naca::parse("2412").unwrap(),
            naca_text: "2412".into(),
            re_index: 4,
            tripped: false,
            alpha: 4.0,
            chord: 0.2,
            polar_mode: PolarMode::Alpha,
            cam: (32.0, 34.0),
            bg: Arc::new(Mutex::new(Background::default())),
            requested: None,
            current: None,
        }
    }

    fn key(&self) -> Key {
        Key {
            naca: self.naca.code(),
            re: W_RE[self.re_index],
            tripped: self.tripped,
            params: self.params,
        }
    }

    /// Lance le calcul des polaires (manquantes) et de la polaire de l'aile.
    fn request(&mut self, ctx: &egui::Context) {
        let key = self.key();
        if self.requested.as_ref() == Some(&key) {
            return;
        }
        self.requested = Some(key.clone());
        let gen = {
            let mut b = self.bg.lock().unwrap();
            b.gen += 1;
            b.gen
        };
        let bg = self.bg.clone();
        let naca = self.naca;
        let ctx = ctx.clone();
        std::thread::spawn(move || {
            let wing = Wing::new(key.params);
            let (rr, rt) = Sections::reynolds(&wing, key.re);
            let table = |re: f64| -> Option<PolarTable> {
                let k = (key.naca.clone(), re.round() as u64, key.tripped);
                if let Some(t) = bg.lock().unwrap().tables.get(&k) {
                    return Some((**t).clone());
                }
                let t = PolarTable::naca(naca, re, key.tripped);
                let mut b = bg.lock().unwrap();
                if b.gen != gen {
                    // Réglages changés entre-temps : la table reste utile.
                    b.tables.insert(k, Arc::new(t));
                    return None;
                }
                b.tables.insert(k, Arc::new(t.clone()));
                Some(t)
            };
            let Some(root) = table(rr) else { return };
            let tip = match rt {
                Some(r) => match table(r) {
                    Some(t) => Some(t),
                    None => return,
                },
                None => None,
            };
            if bg.lock().unwrap().gen != gen {
                return;
            }
            let sections = Sections::new(&wing, key.re, root, tip);
            let alphas: Vec<f64> = (POLAR_LO..=POLAR_HI).map(|a| a as f64).collect();
            let polar = wing_polar(&wing, &sections, &alphas);
            let mut b = bg.lock().unwrap();
            if b.gen == gen {
                b.done = Some(Arc::new(Computed {
                    key,
                    wing,
                    sections,
                    polar,
                }));
                ctx.request_repaint();
            }
        });
    }

    /// Résout l'incidence courante si le calcul ou l'incidence ont changé.
    fn update_current(&mut self, comp: &Arc<Computed>) {
        if let Some((c, a, _)) = &self.current {
            if Arc::ptr_eq(c, comp) && *a == self.alpha {
                return;
            }
        }
        let g0 = comp
            .polar
            .iter()
            .min_by(|p, q| (p.alpha - self.alpha).abs().total_cmp(&(q.alpha - self.alpha).abs()))
            .map(|p| p.gamma.clone());
        let section = |i: usize, a: f64| comp.sections.at(i, a);
        let r = solve_wing(&comp.wing, self.alpha, &section, g0.as_deref());
        self.current = Some((comp.clone(), self.alpha, r));
    }

    pub fn ui(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let key = self.key();
        let done = self.bg.lock().unwrap().done.clone();
        let up_to_date = done.as_ref().is_some_and(|d| d.key == key);
        if !up_to_date {
            self.request(&ctx);
        }
        if let Some(d) = &done {
            self.update_current(d);
        }
        egui::Panel::left("aile_reglages")
            .default_size(340.0)
            .resizable(true)
            .show(ui, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| self.controls(ui, done.as_deref(), up_to_date));
            });
        egui::Panel::bottom("aile_graphes")
            .default_size(300.0)
            .resizable(true)
            .show(ui, |ui| {
                ui.columns(2, |cols| {
                    self.span_chart(&mut cols[0], done.as_deref());
                    self.polar_chart(&mut cols[1], done.as_deref());
                });
            });
        egui::CentralPanel::default().show(ui, |ui| self.view(ui, done.as_deref(), up_to_date));
    }

    fn controls(&mut self, ui: &mut egui::Ui, done: Option<&Computed>, up_to_date: bool) {
        ui.heading("Profil");
        ui.horizontal(|ui| {
            ui.label("NACA");
            let r = ui.add(egui::TextEdit::singleline(&mut self.naca_text).desired_width(60.0));
            if r.changed() {
                if let Some(n) = Naca::parse(self.naca_text.trim()) {
                    self.naca = n;
                }
            }
            for p in ["0012", "2412", "4412", "0009"] {
                if ui.small_button(p).clicked() {
                    self.naca_text = p.into();
                    self.naca = Naca::parse(p).unwrap();
                }
            }
        });
        if Naca::parse(self.naca_text.trim()).is_none() {
            ui.colored_label(Color32::RED, "Code NACA à 4 chiffres invalide.");
        }
        ui.horizontal(|ui| {
            ui.label("Transition");
            ui.selectable_value(&mut self.tripped, false, "libre (eᴺ)");
            ui.selectable_value(&mut self.tripped, true, "bord d'attaque");
        });
        ui.horizontal(|ui| {
            ui.label("Re (corde moyenne)");
            let mut k = self.re_index as f64;
            ui.add(
                egui::Slider::new(&mut k, 0.0..=(W_RE.len() - 1) as f64)
                    .step_by(1.0)
                    .show_value(false),
            );
            self.re_index = k.round() as usize;
            ui.label(crate::gui::fmt_re(W_RE[self.re_index]));
        });
        ui.separator();
        ui.heading("Aile");
        let p = &mut self.params;
        ui.add(egui::Slider::new(&mut self.alpha, -8.0..=22.0).step_by(0.5).text("incidence α (°)"));
        ui.add(egui::Slider::new(&mut p.ar, 2.0..=20.0).step_by(0.5).text("allongement A"));
        ui.add(egui::Slider::new(&mut p.taper, 0.1..=1.0).step_by(0.05).text("effilement λ"));
        ui.add(egui::Slider::new(&mut p.sweep, -30.0..=45.0).step_by(1.0).text("flèche (°)"));
        ui.add(egui::Slider::new(&mut p.dihedral, -10.0..=15.0).step_by(0.5).text("dièdre (°)"));
        ui.add(egui::Slider::new(&mut p.twist, -8.0..=4.0).step_by(0.5).text("vrillage en bout (°)"));
        ui.separator();
        ui.heading("Résultats");
        let Some((_, _, r)) = &self.current else {
            ui.label("Calcul des polaires de section…");
            return;
        };
        if !up_to_date {
            ui.label(egui::RichText::new("Mise à jour…").weak());
        }
        egui::Grid::new("aile_coef")
            .num_columns(2)
            .spacing([16.0, 2.0])
            .show(ui, |ui| {
                let row = |ui: &mut egui::Ui, name: &str, c: Option<Color32>, v: String| {
                    match c {
                        Some(c) => ui.colored_label(c, name),
                        None => ui.label(name),
                    };
                    ui.label(v);
                    ui.end_row();
                };
                row(ui, "CL", Some(LIFT), format!("{:.4}", r.cl));
                row(ui, "CD", Some(ACCENT), format!("{:.5}", r.cd));
                row(ui, "  induite CDi", None, format!("{:.5}", r.cdi));
                row(ui, "  de profil CDp", None, format!("{:.5}", r.cdv));
                row(ui, "Cm (quart de c moyenne)", None, format!("{:.4}", r.cm));
                row(
                    ui,
                    "finesse L/D",
                    None,
                    if r.cd.abs() > 1e-6 {
                        format!("{:.1}", r.cl / r.cd)
                    } else {
                        "—".into()
                    },
                );
                row(
                    ui,
                    "Oswald e",
                    None,
                    if r.e.is_finite() && r.cl.abs() > 0.05 {
                        format!("{:.3}", r.e)
                    } else {
                        "—".into()
                    },
                );
            });
        if let Some(d) = done {
            if let Some(best) = d.polar.iter().max_by(|a, b| a.cl.total_cmp(&b.cl)) {
                ui.label(
                    egui::RichText::new(format!(
                        "CL max ≈ {:.2} à α = {}° (Cl max du profil : {:.2}).",
                        best.cl, best.alpha, d.sections.root.clmax
                    ))
                    .small(),
                );
            }
        }
        // Grandeurs physiques dans l'air.
        ui.add_space(6.0);
        ui.horizontal(|ui| {
            ui.label("Corde moyenne (m)");
            ui.add(egui::DragValue::new(&mut self.chord).range(0.01..=20.0).speed(0.01));
        });
        let re = W_RE[self.re_index];
        let v = re * NU_AIR / self.chord;
        let q = 0.5 * RHO_AIR * v * v;
        let s = self.params.ar * self.chord * self.chord;
        ui.label(
            egui::RichText::new(format!(
                "Air à 15 °C : V = {v:.1} m/s ({:.0} km/h), envergure {:.2} m, S = {s:.3} m², portance {:.1} N, traînée {:.2} N{}",
                v * 3.6,
                self.params.ar * self.chord,
                q * s * r.cl,
                q * s * r.cd,
                if v / 340.0 > 0.3 {
                    format!(" · Mach {:.2} : compressibilité non prise en compte", v / 340.0)
                } else {
                    String::new()
                }
            ))
            .small(),
        );
        // Avertissements.
        let mut msgs = Vec::new();
        if !r.converged {
            msgs.push("Solution non convergée : après le décrochage, la ligne portante admet plusieurs solutions ; résultat indicatif.".to_string());
        }
        if let Some((c, _, _)) = &self.current {
            let (worst, eta) = r
                .strips
                .iter()
                .zip(&c.sections.clmax)
                .map(|(s, m)| (s.cl / m, s.eta.abs()))
                .fold((0.0, 0.0), |acc, x| if x.0 > acc.0 { x } else { acc });
            if worst > 0.97 {
                msgs.push(format!(
                    "Décrochage local amorcé vers 2y/b ≈ {eta:.2} (Cl local au maximum de la polaire du profil)."
                ));
            }
        }
        if self.params.sweep.abs() > 30.0 {
            msgs.push("Flèche > 30° : la ligne portante représente mal le centre d'une aile en flèche ; ordre de grandeur seulement.".into());
        }
        if self.params.ar < 4.0 {
            msgs.push("Allongement < 4 : la ligne portante surestime la pente de portance.".into());
        }
        if re < 2e5 {
            msgs.push("Re < 2·10⁵ : bulbes de décollement laminaire non modélisés dans les polaires de section.".into());
        }
        for m in msgs {
            ui.colored_label(WARN, egui::RichText::new(m).small());
        }
    }

    /// Vue 3D orbitale : aile colorée par Cl/Cl max, tourbillons marginaux.
    fn view(&mut self, ui: &mut egui::Ui, done: Option<&Computed>, up_to_date: bool) {
        let (resp, painter) = ui.allocate_painter(ui.available_size(), Sense::drag());
        let rect = resp.rect;
        painter.rect_filled(rect, 0.0, Color32::from_rgb(0xe6, 0xeb, 0xee));
        if resp.dragged() {
            let d = resp.drag_delta();
            self.cam.0 -= d.x * 0.5;
            self.cam.1 = (self.cam.1 + d.y * 0.4).clamp(-30.0, 85.0);
        }
        let (Some(comp), Some((_, _, r))) = (done, &self.current) else {
            painter.text(
                rect.center(),
                egui::Align2::CENTER_CENTER,
                "Calcul des polaires de section…",
                egui::FontId::proportional(16.0),
                Color32::GRAY,
            );
            return;
        };
        // La géométrie affichée suit les réglages (immédiat), la couleur le dernier calcul.
        let wing = Wing::new(self.params);
        let colors: Vec<Color32> = if up_to_date {
            r.strips
                .iter()
                .zip(&comp.sections.clmax)
                .map(|(s, m)| ramp(s.cl / m))
                .collect()
        } else {
            vec![Color32::from_rgb(150, 160, 170); wing.strips.len()]
        };
        // Incidence : rotation autour de y (bord d'attaque vers le haut).
        let (sa, ca) = self.alpha.to_radians().sin_cos();
        let pitch = |p: [f64; 3]| [p[0] * ca + p[2] * sa, p[1], -p[0] * sa + p[2] * ca];
        let (saz, caz) = (self.cam.0 as f64).to_radians().sin_cos();
        let (sel, cel) = (self.cam.1 as f64).to_radians().sin_cos();
        // Projection orthographique : (écran horizontal, écran vertical, profondeur).
        let proj = |p: [f64; 3]| {
            let u = -p[0] * saz + p[1] * caz;
            let w = p[0] * caz + p[1] * saz;
            let v = p[2] * cel - w * sel;
            let depth = w * cel + p[2] * sel;
            (u, v, depth)
        };
        let scale = 0.8 * rect.width().min(rect.height() * 1.6) as f64 / (wing.b + 2.0);
        let center = rect.center();
        let cx = 0.5 * wing.cr;
        let to_screen = |p: [f64; 3]| {
            let q = pitch([p[0] - cx, p[1], p[2]]);
            let (u, v, d) = proj(q);
            (Pos2::new(center.x + (u * scale) as f32, center.y - (v * scale) as f32), d)
        };
        // Panneaux : bord d'attaque à c/4 en amont du segment lié.
        let mut quads: Vec<(f64, [Pos2; 4], Color32)> = wing
            .strips
            .iter()
            .zip(&colors)
            .map(|(s, &col)| {
                let (ca_, cb) = (wing.chord_at(s.a[1]), wing.chord_at(s.b[1]));
                let pts = [
                    [s.a[0] - 0.25 * ca_ * s.ua[0], s.a[1], s.a[2] - 0.25 * ca_ * s.ua[2]],
                    [s.b[0] - 0.25 * cb * s.ua[0], s.b[1], s.b[2] - 0.25 * cb * s.ua[2]],
                    [s.b[0] + 0.75 * cb * s.ua[0], s.b[1], s.b[2] + 0.75 * cb * s.ua[2]],
                    [s.a[0] + 0.75 * ca_ * s.ua[0], s.a[1], s.a[2] + 0.75 * ca_ * s.ua[2]],
                ];
                let sp: Vec<(Pos2, f64)> = pts.iter().map(|&p| to_screen(p)).collect();
                let depth = sp.iter().map(|x| x.1).sum::<f64>() / 4.0;
                (depth, [sp[0].0, sp[1].0, sp[2].0, sp[3].0], col)
            })
            .collect();
        quads.sort_by(|a, b| b.0.total_cmp(&a.0));
        let painter = painter.with_clip_rect(rect);
        // Tourbillons marginaux : dans le prolongement de la corde, sur 3 cordes.
        let tips = [wing.strips.first().unwrap().a, wing.strips.last().unwrap().b];
        for t in tips {
            let c = wing.ct;
            let p0 = to_screen([t[0] + 0.75 * c, t[1], t[2]]).0;
            let p1 = to_screen([t[0] + 0.75 * c + 3.0, t[1], t[2]]).0;
            painter.line_segment([p0, p1], Stroke::new(1.2, Color32::from_rgb(0x6b, 0x7a, 0x86)));
        }
        let edge = Stroke::new(0.5, Color32::from_rgba_premultiplied(13, 22, 28, 120));
        for (_, q, col) in &quads {
            painter.add(Shape::convex_polygon(q.to_vec(), *col, edge));
        }
        // Vent relatif.
        let a0 = Pos2::new(rect.left() + 24.0, rect.top() + 40.0);
        let (u, v, _) = proj([1.0, 0.0, 0.0]);
        let dir = Vec2::new(u as f32, -v as f32);
        let dir = if dir.length() > 0.1 { dir.normalized() } else { Vec2::X };
        painter.arrow(a0, dir * 40.0, Stroke::new(1.5, INK));
        painter.text(a0 + Vec2::new(0.0, -14.0), egui::Align2::LEFT_BOTTOM, "V∞", egui::FontId::proportional(12.0), INK);
        painter.text(
            rect.left_bottom() + Vec2::new(8.0, -8.0),
            egui::Align2::LEFT_BOTTOM,
            format!(
                "Couleur : Cl local / Cl max du profil (rouge = décrochage local)   ·   glisser pour tourner   ·   α = {}°",
                self.alpha
            ),
            egui::FontId::proportional(13.0),
            INK,
        );
    }

    /// Répartition en envergure : Cl local, Cl max, charge c·Cl/c̄.
    fn span_chart(&self, ui: &mut egui::Ui, done: Option<&Computed>) {
        ui.label(egui::RichText::new("Répartition en envergure (2y/b)").strong());
        let (Some(comp), Some((_, _, r))) = (done, &self.current) else { return };
        let cl: Vec<[f64; 2]> = r.strips.iter().map(|s| [s.eta, s.cl]).collect();
        let ccl: Vec<[f64; 2]> = r.strips.iter().map(|s| [s.eta, s.ccl]).collect();
        let cm: Vec<[f64; 2]> = r.strips.iter().zip(&comp.sections.clmax).map(|(s, &m)| [s.eta, m]).collect();
        Plot::new("aile_envergure")
            .height(ui.available_height().max(80.0))
            .legend(Legend::default())
            .allow_scroll(false)
            .x_axis_label("2y/b")
            .show(ui, |p| {
                p.line(Line::new("Cl local", PlotPoints::from(cl)).color(LIFT));
                p.line(Line::new("charge c·Cl/c̄", PlotPoints::from(ccl)).color(INK));
                p.line(
                    Line::new("Cl max du profil", PlotPoints::from(cm))
                        .color(ACCENT)
                        .style(LineStyle::dashed_dense()),
                );
            });
    }

    /// Polaire de l'aile : CL(α) et CD(α), ou CL(CD) ; profil 2D en tirets.
    fn polar_chart(&mut self, ui: &mut egui::Ui, done: Option<&Computed>) {
        ui.horizontal(|ui| {
            ui.label(egui::RichText::new("Polaire").strong());
            ui.selectable_value(&mut self.polar_mode, PolarMode::Alpha, "CL (α)");
            ui.selectable_value(&mut self.polar_mode, PolarMode::Eiffel, "CL (CD)");
        });
        let Some(comp) = done else { return };
        let root = &comp.sections.root;
        let eiffel = self.polar_mode == PolarMode::Eiffel;
        let cur = self.current.as_ref().map(|c| (c.2.cl, c.2.cd));
        Plot::new(if eiffel { "aile_eiffel" } else { "aile_polaire" })
            .height(ui.available_height().max(80.0))
            .legend(Legend::default())
            .allow_scroll(false)
            .x_axis_label(if eiffel { "CD" } else { "α (°)" })
            .show(ui, |p| {
                if eiffel {
                    let w: Vec<[f64; 2]> = comp.polar.iter().map(|r| [r.cd, r.cl]).collect();
                    let s: Vec<[f64; 2]> = (POLAR_LO..=POLAR_HI)
                        .map(|a| {
                            let k = root.at(a as f64);
                            [k.cd, k.cl]
                        })
                        .collect();
                    p.line(Line::new("aile", PlotPoints::from(w)).color(LIFT));
                    p.line(
                        Line::new("profil 2D (emplanture)", PlotPoints::from(s))
                            .color(Color32::GRAY)
                            .style(LineStyle::dashed_dense()),
                    );
                    if let Some((cl, cd)) = cur {
                        p.points(Points::new("α courant", PlotPoints::from(vec![[cd, cl]])).color(ACCENT).radius(4.0));
                    }
                } else {
                    let w: Vec<[f64; 2]> = comp.polar.iter().map(|r| [r.alpha, r.cl]).collect();
                    let d: Vec<[f64; 2]> = comp.polar.iter().map(|r| [r.alpha, 10.0 * r.cd]).collect();
                    let s: Vec<[f64; 2]> = (POLAR_LO..=POLAR_HI).map(|a| [a as f64, root.at(a as f64).cl]).collect();
                    p.line(Line::new("CL aile", PlotPoints::from(w)).color(LIFT));
                    p.line(Line::new("10 × CD aile", PlotPoints::from(d)).color(ACCENT));
                    p.line(
                        Line::new("Cl profil 2D (emplanture)", PlotPoints::from(s))
                            .color(Color32::GRAY)
                            .style(LineStyle::dashed_dense()),
                    );
                    p.vline(VLine::new("α courant", self.alpha).color(INK).style(LineStyle::dotted_dense()));
                }
            });
    }
}
