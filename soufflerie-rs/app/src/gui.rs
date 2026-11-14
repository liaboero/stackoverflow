//! Interface graphique (egui) : veine d'essai zoomable, traceurs, réglages,
//! coefficients, référence « panneaux + couche limite », historique des
//! efforts et distributions pariétales.

use std::sync::{Arc, Mutex};

use eframe::egui::{self, Color32, Pos2, Rect, Sense, Shape, Stroke, Vec2};
use egui_plot::{Legend, Line, LineStyle, Plot, PlotPoints, Points, VLine};
use lbm::cpu::{Collision, LesMode};
use lbm::panel::{section_coefficients, PanelShape};
use lbm::tunnel::RE_PANEL_MIN;
use lbm::wall::WallAxis;
use lbm::{EngineKind, Naca, ShapeKind, ShapeParams, Transition, TunnelConfig, Turbulence};

use crate::sweep::{Sweep, SweepVar};
use crate::tab3d::Tab3D;
use crate::wing_tab::WingTab;
use crate::theme;
use crate::worker::{self, Field, Frame, Requests, Settings, Shared, Tracers};

const RES: [usize; 9] = [12, 16, 24, 32, 48, 64, 96, 128, 192];

/// Réglages initiaux de l'interface (options de la ligne de commande).
#[derive(Clone)]
pub struct Start {
    pub cfg: Option<TunnelConfig>,
    /// vorticity | speed | pressure | eddy
    pub field: Option<String>,
    /// history | cp | cf
    pub graph: Option<String>,
    /// Onglet initial : tunnel | wing
    pub tab: Option<String>,
}

/// Ouvre l'interface ; false si aucun rendu n'a pu être créé.
/// Rendu wgpu par défaut, OpenGL (glow) en repli ; SOUFFLERIE_RENDU=glow|wgpu en impose un.
pub fn main(start: Start) -> bool {
    let renderers: &[eframe::Renderer] = match std::env::var("SOUFFLERIE_RENDU").as_deref() {
        Ok("glow") => &[eframe::Renderer::Glow],
        Ok("wgpu") => &[eframe::Renderer::Wgpu],
        _ => &[eframe::Renderer::Wgpu, eframe::Renderer::Glow],
    };
    for (k, &renderer) in renderers.iter().enumerate() {
        let opts = eframe::NativeOptions {
            viewport: egui::ViewportBuilder::default()
                .with_inner_size([1360.0, 860.0])
                .with_title("Soufflerie de poche"),
            renderer,
            ..Default::default()
        };
        // Chaque tentative de rendu reçoit sa copie des réglages.
        let st = start.clone();
        match eframe::run_native(
            "Soufflerie de poche",
            opts,
            Box::new(move |cc| Ok(Box::new(App::new(cc, st)))),
        ) {
            Ok(()) => return true,
            Err(e) if k + 1 < renderers.len() => eprintln!("Rendu {renderer:?} indisponible ({e}) : essai suivant."),
            Err(e) => eprintln!("Rendu {renderer:?} indisponible : {e}"),
        }
    }
    false
}

struct App {
    shared: Arc<Mutex<Shared>>,
    settings: Settings,
    naca_text: String,
    texture: Option<egui::TextureHandle>,
    image_id: u64,
    /// Vue : centre (réseau) et largeur affichée (mailles) ; None = cadrage maquette.
    view: Option<(Pos2, f32)>,
    show_force: bool,
    bottom: Bottom,
    tab: Tab,
    wing: WingTab,
    tab3d: Tab3D,
    /// Polaire automatique et points enregistrés.
    sweep: Sweep,
    /// Échelle physique : fluide et longueur de référence (m).
    fluid: Fluid,
    chord: f64,
    /// État de marche de la soufflerie 2D avant le passage à un autre onglet.
    resume_tunnel: bool,
}

/// Onglet affiché.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Tab {
    Tunnel,
    Wing,
    Tunnel3D,
}

/// Graphique du bas.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Bottom {
    History,
    Cp,
    Cf,
    Polar,
}

/// Fluide de l'échelle physique.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Fluid {
    Air,
    Water,
}

impl Fluid {
    fn label(self) -> &'static str {
        match self {
            Fluid::Air => "Air 15 °C",
            Fluid::Water => "Eau 20 °C",
        }
    }
    /// Viscosité cinématique (m²/s) et masse volumique (kg/m³).
    fn props(self) -> (f64, f64) {
        match self {
            Fluid::Air => (1.46e-5, 1.225),
            Fluid::Water => (1.0e-6, 998.0),
        }
    }
}

/// Nombre lisible à trois chiffres significatifs.
fn sig3(v: f64) -> String {
    if !v.is_finite() {
        return "—".into();
    }
    let a = v.abs();
    if a >= 1e4 || (a > 0.0 && a < 1e-3) {
        return format!("{v:.2e}").replace('e', "·10^");
    }
    let d = if a >= 100.0 {
        0
    } else if a >= 10.0 {
        1
    } else {
        2
    };
    format!("{v:.d$}").replace('.', ",")
}

/// Durée lisible.
fn duration(sec: f64) -> String {
    if !sec.is_finite() {
        "—".into()
    } else if sec < 60.0 {
        format!("{} s", sig3(sec))
    } else if sec < 3600.0 {
        format!("{} min {:02} s", (sec / 60.0).floor(), (sec % 60.0).round())
    } else {
        format!("{:.1} h", sec / 3600.0).replace('.', ",")
    }
}

impl App {
    fn new(cc: &eframe::CreationContext<'_>, start: Start) -> App {
        cc.egui_ctx.set_theme(egui::ThemePreference::System);
        match &cc.wgpu_render_state {
            Some(rs) => {
                let i = rs.adapter.get_info();
                eprintln!("Fenêtre ouverte : rendu wgpu, {} ({:?}).", i.name, i.backend);
            }
            None => eprintln!("Fenêtre ouverte : rendu OpenGL."),
        }
        let cfg = start.cfg.unwrap_or(TunnelConfig {
            l: 32,
            re: 1000.0,
            shape: ShapeKind::Naca,
            params: ShapeParams {
                alpha_deg: 6.0,
                naca: Naca::parse("2412").unwrap(),
                ..Default::default()
            },
            ..Default::default()
        });
        let field = match start.field.as_deref() {
            Some("speed") => Field::Speed,
            Some("pressure") => Field::Pressure,
            Some("eddy") => Field::Eddy,
            _ => Field::Vorticity,
        };
        let bottom = match start.graph.as_deref() {
            Some("cp") => Bottom::Cp,
            Some("cf") => Bottom::Cf,
            _ => Bottom::History,
        };
        let naca_text = cfg.params.naca.code();
        let tab = match start.tab.as_deref() {
            Some("wing") => Tab::Wing,
            Some("3d") => Tab::Tunnel3D,
            _ => Tab::Tunnel,
        };
        let settings = Settings {
            cfg,
            field,
            tracers: Tracers::Fronts,
            // Sur les autres onglets, la soufflerie 2D attend.
            running: tab == Tab::Tunnel,
            wall_corr: true,
        };
        let shared = Arc::new(Mutex::new(Shared {
            settings: settings.clone(),
            req: Requests::default(),
            frame: Arc::new(Frame::default()),
        }));
        let ctx = cc.egui_ctx.clone();
        worker::spawn(shared.clone(), move || ctx.request_repaint());
        let mut tab3d = Tab3D::new();
        if tab == Tab::Tunnel3D {
            tab3d.set_active(true, &cc.egui_ctx);
        }
        App {
            shared,
            settings,
            naca_text,
            texture: None,
            image_id: 0,
            view: None,
            show_force: true,
            bottom,
            tab,
            wing: WingTab::new(),
            tab3d,
            sweep: Sweep::new(),
            fluid: Fluid::Air,
            chord: 0.2,
            resume_tunnel: true,
        }
    }

    /// Change d'onglet ; seul le calcul de l'onglet affiché tourne.
    fn set_tab(&mut self, tab: Tab, ctx: &egui::Context) {
        if tab == self.tab {
            return;
        }
        if self.tab == Tab::Tunnel {
            self.resume_tunnel = self.settings.running;
            self.settings.running = false;
        } else if tab == Tab::Tunnel {
            self.settings.running = self.resume_tunnel;
        }
        self.tab3d.set_active(tab == Tab::Tunnel3D, ctx);
        self.tab = tab;
        self.push(|_| {});
    }

    /// Transmet les réglages et une demande au fil de calcul.
    fn push(&self, f: impl FnOnce(&mut Requests)) {
        let mut s = self.shared.lock().unwrap();
        s.settings = self.settings.clone();
        f(&mut s.req);
    }

    fn controls(&mut self, ui: &mut egui::Ui, fr: &Frame) {
        self.sweep_tick(fr);
        let sweeping = self.sweep.active();
        ui.add_enabled_ui(!sweeping, |ui| {
            ui.heading("Maquette");
            let mut shape = self.settings.cfg.shape;
            egui::ComboBox::from_id_salt("shape")
                .selected_text(shape.label())
                .show_ui(ui, |ui| {
                    for k in ShapeKind::ALL {
                        ui.selectable_value(&mut shape, k, k.label());
                    }
                });
            let mut body_changed = shape != self.settings.cfg.shape;
            self.settings.cfg.shape = shape;
            if shape == ShapeKind::Naca {
                ui.horizontal(|ui| {
                    ui.label("NACA");
                    let r = ui.add(egui::TextEdit::singleline(&mut self.naca_text).desired_width(60.0));
                    if r.changed() {
                        if let Some(n) = Naca::parse(&self.naca_text) {
                            self.settings.cfg.params.naca = n;
                            body_changed = true;
                        }
                    }
                    for code in ["0012", "2412", "4412"] {
                        if ui.small_button(code).clicked() {
                            self.naca_text = code.into();
                            self.settings.cfg.params.naca = Naca::parse(code).unwrap();
                            body_changed = true;
                        }
                    }
                });
            }
            if let Some((lo, hi)) = shape.alpha_range() {
                let r = ui.add(
                    egui::Slider::new(&mut self.settings.cfg.params.alpha_deg, lo..=hi)
                        .step_by(0.5)
                        .text("incidence α (°)"),
                );
                body_changed |= r.drag_stopped() || (r.changed() && !r.dragged());
            }
            if shape == ShapeKind::Ellipse {
                let r = ui.add(egui::Slider::new(&mut self.settings.cfg.params.thickness, 0.1..=1.0).text("épaisseur e/c"));
                body_changed |= r.drag_stopped() || (r.changed() && !r.dragged());
            }
            if shape == ShapeKind::Cylinder {
                let r = ui.add(
                    egui::Slider::new(&mut self.settings.cfg.params.spin_ratio, -3.0..=3.0)
                        .step_by(0.25)
                        .text("rotation ΩR/U"),
                );
                body_changed |= r.drag_stopped() || (r.changed() && !r.dragged());
            }
            if body_changed {
                self.push(|r| r.body = true);
            }

            ui.separator();
            ui.heading("Écoulement");
            let r = ui.add(
                egui::Slider::new(&mut self.settings.cfg.re, 10.0..=1e7)
                    .logarithmic(true)
                    .custom_formatter(|v, _| fmt_re(v))
                    .text("Reynolds"),
            );
            if r.drag_stopped() || (r.changed() && !r.dragged()) {
                self.push(|q| q.reynolds = true);
            }
            let t0 = (self.settings.cfg.turbulence, self.settings.cfg.transition);
            ui.horizontal(|ui| {
                ui.label("Turbulence");
                let t = &mut self.settings.cfg.turbulence;
                ui.selectable_value(t, Turbulence::Auto, "auto")
                    .on_hover_text("Spalart–Allmaras avec loi de paroi dès Re = 5·10⁴, aucun modèle en dessous");
                ui.selectable_value(t, Turbulence::Sa, "Spalart–Allmaras");
                ui.selectable_value(t, Turbulence::Off, "aucune");
            });
            let profiled = PanelShape::from_shape(shape, &self.settings.cfg.params).is_some();
            if fr.sa && profiled {
                ui.horizontal(|ui| {
                    ui.label("Transition");
                    let t = &mut self.settings.cfg.transition;
                    ui.selectable_value(t, Transition::Free, "libre").on_hover_text(
                        "Position prédite par la méthode des panneaux et la couche limite intégrale (critère e^N, N = 9) ; \
                         paroi laminaire en amont",
                    );
                    ui.selectable_value(t, Transition::Turbulent, "bord d'attaque")
                        .on_hover_text("Couche limite turbulente partout");
                });
            }
            if (self.settings.cfg.turbulence, self.settings.cfg.transition) != t0 {
                self.push(|q| q.turbulence = true);
            }
            let mut rebuild = false;
            let l0 = self.settings.cfg.l;
            egui::ComboBox::from_label("mailles par longueur L")
                .selected_text(format!("{l0}"))
                .show_ui(ui, |ui| {
                    for l in RES {
                        let cells = 24 * l * 12 * l;
                        ui.selectable_value(
                            &mut self.settings.cfg.l,
                            l,
                            format!("{l}  ({:.2} M mailles)", cells as f64 / 1e6),
                        );
                    }
                });
            rebuild |= self.settings.cfg.l != l0;
            let e0 = self.settings.cfg.engine;
            ui.horizontal(|ui| {
                ui.label("Moteur");
                let cpu = EngineKind::Cpu { threads: 0 };
                ui.selectable_value(&mut self.settings.cfg.engine, cpu, "CPU (tous les cœurs)");
                if cfg!(feature = "gpu") {
                    ui.selectable_value(&mut self.settings.cfg.engine, EngineKind::Gpu, "GPU");
                }
            });
            if self.settings.cfg.engine == EngineKind::Gpu {
                let h0 = self.settings.cfg.gpu_half;
                ui.horizontal(|ui| {
                    ui.label("Précision GPU");
                    ui.selectable_value(&mut self.settings.cfg.gpu_half, false, "simple");
                    ui.selectable_value(&mut self.settings.cfg.gpu_half, true, "demi")
                        .on_hover_text("Populations sur 2 octets : deux fois moins de mémoire, plus rapide ; écart ≈ 0,2 % sur les efforts moyens");
                });
                rebuild |= self.settings.cfg.gpu_half != h0;
                let r0 = self.settings.cfg.refine;
                ui.horizontal(|ui| {
                    ui.label("Près de la maquette");
                    let l = self.settings.cfg.l;
                    ui.selectable_value(&mut self.settings.cfg.refine, 0, "uniforme");
                    ui.selectable_value(&mut self.settings.cfg.refine, 1, format!("×2 ({})", 2 * l))
                        .on_hover_text("Grille deux fois plus fine autour de la maquette");
                    ui.selectable_value(&mut self.settings.cfg.refine, 2, format!("×4 ({})", 4 * l))
                        .on_hover_text(
                            "Deux grilles imbriquées, 2 et 4 fois plus fines : couche limite, bords d'attaque \
                             et de fuite résolus à 4 L pour un coût proche de la grille L ; \
                             traînée des profils à grand Reynolds nettement plus juste",
                        );
                });
                rebuild |= self.settings.cfg.refine != r0;
            }
            rebuild |= self.settings.cfg.engine != e0;
            let c0 = (self.settings.cfg.collision, self.settings.cfg.les);
            ui.horizontal(|ui| {
                ui.label("Collision");
                ui.selectable_value(&mut self.settings.cfg.collision, Collision::Rr, "régularisée");
                ui.selectable_value(&mut self.settings.cfg.collision, Collision::Bgk, "BGK");
            });
            ui.horizontal(|ui| {
                ui.label("Sous-maille");
                ui.selectable_value(&mut self.settings.cfg.les, LesMode::Auto, "auto");
                ui.selectable_value(&mut self.settings.cfg.les, LesMode::On, "oui");
                ui.selectable_value(&mut self.settings.cfg.les, LesMode::Off, "non");
            });
            rebuild |= (self.settings.cfg.collision, self.settings.cfg.les) != c0;
            ui.horizontal(|ui| {
                let label = if self.settings.running { "Pause" } else { "Reprendre" };
                if ui.button(label).clicked() {
                    self.settings.running = !self.settings.running;
                    self.push(|_| {});
                }
                if ui.button("Relancer").clicked() {
                    rebuild = true;
                }
            });
            if rebuild {
                self.push(|r| r.rebuild = true);
            }

        });
        self.sweep_panel(ui, fr);
        self.physical(ui, fr);
        ui.separator();
        ui.heading("Affichage");
        let f0 = (self.settings.field, self.settings.tracers);
        ui.horizontal_wrapped(|ui| {
            ui.selectable_value(&mut self.settings.field, Field::Vorticity, "Vorticité");
            ui.selectable_value(&mut self.settings.field, Field::Speed, "Vitesse");
            ui.selectable_value(&mut self.settings.field, Field::Pressure, "Pression");
            if fr.sa || self.settings.field == Field::Eddy {
                ui.selectable_value(&mut self.settings.field, Field::Eddy, "Viscosité turbulente");
            }
        });
        ui.horizontal(|ui| {
            ui.label("Traceurs");
            ui.selectable_value(&mut self.settings.tracers, Tracers::Off, "aucun");
            ui.selectable_value(&mut self.settings.tracers, Tracers::Streak, "filets");
            ui.selectable_value(&mut self.settings.tracers, Tracers::Fronts, "fronts");
        });
        if (self.settings.field, self.settings.tracers) != f0 {
            self.push(|_| {});
        }
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.show_force, "Force résultante");
            if ui.button("Recadrer").clicked() {
                self.view = None;
            }
        });
        ui.label(
            egui::RichText::new("Molette : zoom · glisser : déplacer · double-clic : recadrer")
                .small()
                .weak(),
        );

        ui.separator();
        ui.heading("Coefficients");
        let a = &fr.analysis;
        let status = if fr.diverged {
            "Divergence : réduisez Re ou affinez la grille"
        } else if !self.settings.running {
            "En pause"
        } else if fr.ramping {
            "Mise en vitesse"
        } else if !fr.ready {
            "Transitoire"
        } else if a.periodic && a.cz_amp > 0.01 {
            "Périodique"
        } else {
            "Établi"
        };
        ui.label(egui::RichText::new(status).strong());
        egui::Grid::new("coeffs")
            .num_columns(2)
            .spacing([16.0, 4.0])
            .show(ui, |ui| {
                let v = |x: f64| {
                    if fr.ready && x.is_finite() {
                        format!("{x:.4}")
                    } else {
                        "—".into()
                    }
                };
                let amp = |x: f64| {
                    if fr.ready && a.periodic && x > 5e-4 {
                        format!("± {x:.4}")
                    } else {
                        String::new()
                    }
                };
                let cx_ok = crate::sweep::cx_reliable(fr);
                if cx_ok {
                    ui.colored_label(theme::accent(), "Cx (traînée)");
                    ui.label(format!("{} {}", v(a.cx), amp(a.cx_amp)));
                    ui.end_row();
                } else {
                    let why = "Profil attaché sous loi de paroi : le bord d'attaque tient dans une maille, \
                               et quelques dixièmes de degré d'erreur sur la direction de la résultante \
                               basculent de la portance dans la traînée. La traînée LBM ne converge pas \
                               en affinant la grille (NACA 2412, Re = 10⁶, 4° : 0,016 à 0,022 de L = 32 à 96, \
                               pour 0,008 attendu). Traînée retenue : panneaux + couche limite.";
                    ui.colored_label(theme::warn(), "Cx LBM ⚠").on_hover_text(why);
                    ui.colored_label(theme::warn(), format!("{} (non fiable)", v(a.cx))).on_hover_text(why);
                    ui.end_row();
                    if let Some(r) = &fr.panel {
                        ui.colored_label(theme::accent(), "Cx retenu");
                        ui.label(format!("{:.4} (panneaux)", section_coefficients(r).cd));
                        ui.end_row();
                    }
                }
                ui.colored_label(theme::lift(), "Cz (portance)");
                ui.label(format!("{} {}", v(a.cz), amp(a.cz_amp)));
                ui.end_row();
                ui.label("Cm (tangage)");
                ui.label(v(a.cm));
                ui.end_row();
                ui.label("St (Strouhal)");
                ui.label(if fr.ready && a.periodic {
                    format!("{:.3} ({} périodes)", a.st, a.periods)
                } else {
                    "—".into()
                });
                ui.end_row();
                let cx_f = match (&fr.panel, crate::sweep::cx_reliable(fr)) {
                    (Some(r), false) => section_coefficients(r).cd,
                    _ => a.cx,
                };
                ui.label("Cz/Cx");
                ui.label(if fr.ready && cx_f.abs() > 1e-6 {
                    format!("{:.2}", a.cz / cx_f)
                } else {
                    "—".into()
                });
                ui.end_row();
            });
        let w0 = self.settings.wall_corr;
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.settings.wall_corr, "Corriger le blocage de la veine").on_hover_text(
                "La veine ne fait que 12 longueurs de haut : l'écoulement accélère autour de la maquette \
                 et de son sillage. Correction d'Allen et Vincenti (vitesse effective U(1 + ε)) ; \
                 sans elle, la traînée d'un cylindre est surestimée de 8 à 15 %.",
            );
            if self.settings.wall_corr && fr.ready && fr.blockage > 0.0 {
                ui.label(egui::RichText::new(format!("ε = {:.1} %", 100.0 * fr.blockage).replace('.', ",")).small().weak());
            }
        });
        if self.settings.wall_corr != w0 {
            self.push(|_| {});
        }
        self.panel_reference(ui, fr);
        if fr.ready && !self.sweep.active() && ui.button("Ajouter à la polaire").clicked() {
            let (var, value) = self.sweep_var();
            self.sweep.add_current(fr, var, value, self.settings.cfg.re, self.settings.cfg.l);
            self.bottom = Bottom::Polar;
        }
        if fr.ready && ui.button("Copier (CSV)").clicked() {
            let c = &self.settings.cfg;
            ui.ctx().copy_text(format!(
                "maquette,alpha_deg,Re,L,Cx,Cz,Cm,St\n{},{},{},{},{:.5},{:.5},{:.5},{}",
                fr.label,
                c.params.alpha_deg,
                c.re,
                c.l,
                a.cx,
                a.cz,
                a.cm,
                if a.periodic {
                    format!("{:.4}", a.st)
                } else {
                    String::new()
                }
            ));
        }

        ui.separator();
        ui.label(
            egui::RichText::new(format!(
                "{} · grille {}×{} · τ = {:.4} · sous-maille {} · {:.0} MLUPS · t* = {:.1}",
                fr.engine,
                fr.nx,
                fr.ny,
                fr.tau0,
                if fr.les { "active" } else { "inactive" },
                fr.mlups,
                fr.t_star
            ))
            .small()
            .monospace(),
        );
        if fr.sa {
            ui.label(
                egui::RichText::new(format!(
                    "Spalart–Allmaras + loi de paroi · y+ ≈ {} au premier nœud",
                    if fr.yplus < 100.0 {
                        format!("{:.3}", fr.yplus)
                    } else {
                        format!("{:.0}", (fr.yplus / 10.0).round() * 10.0)
                    }
                ))
                .small(),
            );
            ui.label(
                egui::RichText::new(
                    "La couche limite, plus mince qu'une maille, n'est pas résolue : un glissement de paroi \
                     transmet le frottement de la loi de Spalding. Portance et décollement qualitatifs. \
                     La traînée LBM d'un profil porteur dépend de la finesse du bord d'attaque (NACA 0012 \
                     à 4°, Re = 10⁶ : 0,006 à L = 32, 0,020 à L = 48, pour 0,008 attendu) : en écoulement \
                     attaché, prendre la traînée des panneaux.",
                )
                .small()
                .weak(),
            );
        } else if self.settings.cfg.re > fr.reliable_re && fr.reliable_re > 0.0 {
            ui.colored_label(
                theme::warn(),
                format!(
                    "Au-delà de Re ≈ {:.0} avec {} mailles par longueur, la couche limite est trop peu résolue : traînée surestimée, résultats qualitatifs.",
                    fr.reliable_re, self.settings.cfg.l
                ),
            );
        }
        if let Some(e) = &fr.error {
            ui.colored_label(Color32::RED, e);
        }
    }

    /// Variable de polaire de la maquette et sa valeur courante.
    fn sweep_var(&self) -> (SweepVar, f64) {
        let p = &self.settings.cfg.params;
        if self.settings.cfg.shape.alpha_range().is_some() {
            (SweepVar::Alpha, p.alpha_deg)
        } else {
            (SweepVar::Spin, p.spin_ratio)
        }
    }

    /// Impose une valeur de la variable balayée.
    fn set_sweep_value(&mut self, v: f64) {
        match self.sweep_var().0 {
            SweepVar::Alpha => self.settings.cfg.params.alpha_deg = v,
            SweepVar::Spin => self.settings.cfg.params.spin_ratio = v,
        }
        if !self.settings.running {
            self.settings.running = true;
        }
        self.push(|r| r.body = true);
    }

    fn sweep_tick(&mut self, fr: &Frame) {
        if let Some(v) = self.sweep.tick(fr, self.settings.cfg.re, self.settings.cfg.l) {
            self.set_sweep_value(v);
        }
    }

    /// Polaire automatique : bornes, lancement, avancement.
    fn sweep_panel(&mut self, ui: &mut egui::Ui, fr: &Frame) {
        ui.separator();
        let (var, _) = self.sweep_var();
        ui.heading(match var {
            SweepVar::Alpha => "Polaire automatique (incidence)",
            SweepVar::Spin => "Polaire automatique (rotation)",
        });
        let sweeping = self.sweep.active();
        ui.add_enabled_ui(!sweeping, |ui| {
            ui.horizontal(|ui| {
                let unit = if var == SweepVar::Alpha { "°" } else { "" };
                ui.label("de");
                ui.add(egui::DragValue::new(&mut self.sweep.from).speed(0.5).suffix(unit));
                ui.label("à");
                ui.add(egui::DragValue::new(&mut self.sweep.to).speed(0.5).suffix(unit));
                ui.label("pas");
                ui.add(egui::DragValue::new(&mut self.sweep.step).speed(0.1).range(0.05..=45.0).suffix(unit));
            });
        });
        ui.horizontal(|ui| {
            if ui.button(if sweeping { "Arrêter" } else { "Lancer" }).clicked() {
                if sweeping {
                    self.sweep.stop("Polaire interrompue.");
                } else {
                    let lim = self.settings.cfg.shape.alpha_range().unwrap_or((-3.0, 3.0));
                    if let Some(v) = self.sweep.start(var, lim, fr.t_star) {
                        self.set_sweep_value(v);
                        self.bottom = Bottom::Polar;
                    }
                }
            }
            if !self.sweep.points.is_empty() && !sweeping {
                ui.label(egui::RichText::new(format!("{} points", self.sweep.points.len())).weak());
            }
        });
        if let Some((frac, text)) = self.sweep.progress(fr) {
            ui.add(egui::ProgressBar::new(frac).text(text));
        } else if !self.sweep.message.is_empty() {
            ui.label(egui::RichText::new(&self.sweep.message).small());
        }
        ui.label(
            egui::RichText::new(
                "Chaque point : stabilisation (5 à 12 t*), puis mesure (6 à 15 t*, ou trois périodes). \
                 Les autres réglages sont verrouillés pendant la polaire.",
            )
            .small()
            .weak(),
        );
    }

    /// Échelle physique : vitesse, durées, efforts par mètre d'envergure, Mach réel.
    fn physical(&mut self, ui: &mut egui::Ui, fr: &Frame) {
        ui.separator();
        ui.heading("Échelle physique");
        ui.horizontal(|ui| {
            for f in [Fluid::Air, Fluid::Water] {
                ui.selectable_value(&mut self.fluid, f, f.label());
            }
        });
        ui.horizontal(|ui| {
            ui.label("Longueur de référence L (corde, diamètre)");
            ui.add(egui::DragValue::new(&mut self.chord).range(0.001..=50.0).speed(0.01).suffix(" m"));
        });
        let (nu, rho) = self.fluid.props();
        let re = self.settings.cfg.re;
        let v = re * nu / self.chord;
        let t_star = self.chord / v;
        // Pas de calcul par temps convectif : L / U (réseau).
        let steps = fr.l / self.settings.cfg.u_in as f64;
        let mut parts = vec![
            format!("V = {} m/s ({} km/h)", sig3(v), sig3(v * 3.6)),
            format!("1 t* = L/V = {}", duration(t_star)),
            format!("1 pas = {} s", sig3(t_star / steps)),
        ];
        if fr.ready {
            let q = 0.5 * rho * v * v * self.chord;
            parts.push(format!(
                "traînée {} N/m, portance {} N/m d'envergure",
                sig3(fr.analysis.cx * q),
                sig3(fr.analysis.cz * q)
            ));
            if fr.t_star > 0.0 {
                parts.push(format!("temps simulé {}", duration(fr.t_star * t_star)));
            }
        }
        ui.label(egui::RichText::new(parts.join(" · ")).small());
        let mach = if self.fluid == Fluid::Air { v / 340.0 } else { 0.0 };
        if mach > 0.3 {
            ui.colored_label(
                theme::warn(),
                egui::RichText::new(format!(
                    "Mach réel ≈ {} : compressibilité non représentée par le calcul (incompressible). \
                     La référence des panneaux peut être corrigée : soufflerie panel --mach {:.2}.",
                    sig3(mach),
                    mach.min(0.85)
                ))
                .small(),
            );
        }
    }

    /// Polaire enregistrée : Cz et Cx selon la variable, Cz(Cx), points des panneaux.
    fn polar_chart(&mut self, ui: &mut egui::Ui) {
        let pts = &self.sweep.points;
        let (var, _) = self.sweep_var();
        let mut copy = false;
        let mut clear = false;
        ui.horizontal(|ui| {
            copy = ui.button("Copier (CSV)").clicked();
            clear = ui.button("Effacer").clicked();
            ui.label(
                egui::RichText::new(format!(
                    "{} points · ronds : LBM ; tirets : panneaux + couche limite (profils, Re ≥ 2·10⁴) ; \
                     profil attaché sous loi de paroi : Cx retenu = panneaux, Cx LBM en gris",
                    pts.len()
                ))
                .small()
                .weak(),
            );
        });
        if copy {
            ui.ctx().copy_text(self.sweep.csv());
        }
        let mut sel: Vec<&crate::sweep::PolarPoint> = pts.iter().filter(|p| p.var == var).collect();
        sel.sort_by(|a, b| a.value.total_cmp(&b.value));
        let xs = |f: &dyn Fn(&crate::sweep::PolarPoint) -> Option<f64>| -> Vec<[f64; 2]> {
            sel.iter().filter_map(|p| f(p).map(|y| [p.value, y])).collect()
        };
        let cz = xs(&|p| Some(p.cz));
        let cx = xs(&|p| p.cx_ok.then_some(p.cx));
        let cx_bad = xs(&|p| (!p.cx_ok).then_some(p.cx));
        let pcz = xs(&|p| p.panel.map(|q| q.1));
        let pcx = xs(&|p| p.panel.map(|q| q.0));
        let eiffel: Vec<[f64; 2]> = sel.iter().map(|p| [p.cx_retained(), p.cz]).collect();
        let eiffel_bad: Vec<[f64; 2]> = sel.iter().filter(|p| !p.cx_ok).map(|p| [p.cx, p.cz]).collect();
        let grey = Color32::from_gray(140);
        let x_name = if var == SweepVar::Alpha { "α (°)" } else { "ΩR/U" };
        ui.columns(2, |cols| {
            Plot::new("polaire_var")
                .height(cols[0].available_height().max(80.0))
                .legend(Legend::default())
                .allow_scroll(false)
                .x_axis_label(x_name)
                .show(&mut cols[0], |p| {
                    let dash = LineStyle::dashed_dense();
                    p.line(Line::new("Cz", PlotPoints::from(cz.clone())).color(theme::lift()));
                    p.points(Points::new("Cz", PlotPoints::from(cz)).color(theme::lift()).radius(3.0));
                    if !cx.is_empty() {
                        p.line(Line::new("Cx", PlotPoints::from(cx.clone())).color(theme::accent()));
                        p.points(Points::new("Cx", PlotPoints::from(cx)).color(theme::accent()).radius(3.0));
                    }
                    if !cx_bad.is_empty() {
                        p.points(
                            Points::new("Cx LBM (non fiable)", PlotPoints::from(cx_bad))
                                .color(grey)
                                .filled(false)
                                .radius(3.0),
                        );
                    }
                    if !pcz.is_empty() {
                        p.line(Line::new("Cz panneaux", PlotPoints::from(pcz)).color(theme::lift()).style(dash));
                        p.line(Line::new("Cx panneaux", PlotPoints::from(pcx)).color(theme::accent()).style(dash));
                    }
                });
            Plot::new("polaire_eiffel")
                .height(cols[1].available_height().max(80.0))
                .allow_scroll(false)
                .legend(Legend::default())
                .x_axis_label("Cx")
                .y_axis_label("Cz")
                .show(&mut cols[1], |p| {
                    let name = if eiffel_bad.is_empty() { "Cz(Cx)" } else { "Cz(Cx retenu)" };
                    p.line(Line::new(name, PlotPoints::from(eiffel.clone())).color(theme::lift()));
                    p.points(Points::new(name, PlotPoints::from(eiffel)).color(theme::lift()).radius(3.0));
                    if !eiffel_bad.is_empty() {
                        p.points(
                            Points::new("Cx LBM (non fiable)", PlotPoints::from(eiffel_bad))
                                .color(grey)
                                .filled(false)
                                .radius(3.0),
                        );
                    }
                });
        });
        if clear {
            self.sweep.points.clear();
        }
    }

    /// Référence « panneaux + couche limite » (profils, Re ≥ 2·10⁴).
    fn panel_reference(&self, ui: &mut egui::Ui, fr: &Frame) {
        let profiled = PanelShape::from_shape(self.settings.cfg.shape, &self.settings.cfg.params).is_some();
        if !profiled {
            return;
        }
        ui.add_space(4.0);
        let Some(r) = &fr.panel else {
            if self.settings.cfg.re < RE_PANEL_MIN {
                ui.label(
                    egui::RichText::new("Référence « panneaux + couche limite » à partir de Re = 2·10⁴.")
                        .small()
                        .weak(),
                );
            }
            return;
        };
        let k = section_coefficients(r);
        ui.label(egui::RichText::new("Panneaux + couche limite").strong());
        let pct = |v: f64| {
            if v.is_finite() {
                format!("{:.0} %", v.clamp(0.0, 1.0) * 100.0)
            } else {
                "—".into()
            }
        };
        egui::Grid::new("panneaux")
            .num_columns(2)
            .spacing([16.0, 2.0])
            .show(ui, |ui| {
                ui.colored_label(theme::accent(), "Cx");
                ui.label(format!("{:.4}", k.cd));
                ui.end_row();
                ui.colored_label(theme::lift(), "Cz");
                ui.label(format!("{:.4}", k.cl));
                ui.end_row();
                ui.label("Cm");
                ui.label(format!("{:.4}", k.cm));
                ui.end_row();
                ui.label("Transition");
                ui.label(if self.settings.cfg.transition == Transition::Turbulent {
                    "bord d'attaque".to_string()
                } else {
                    format!(
                        "extrados {}, intrados {}",
                        pct(r.transition_upper),
                        pct(r.transition_lower)
                    )
                });
                ui.end_row();
            });
        let note = if k.kirchhoff {
            Some(format!(
                "Décollement à {} de corde : décrochage approché (modèle de Kirchhoff).",
                pct(r.separation_upper.min(r.separation_lower).min(1.0))
            ))
        } else if self.settings.cfg.re < 2e5 && self.settings.cfg.transition == Transition::Free {
            Some("Bulbes de décollement laminaire non modélisés.".to_string())
        } else {
            None
        };
        if let Some(n) = note {
            ui.colored_label(theme::warn(), egui::RichText::new(n).small());
        }
    }

    fn bottom(&mut self, ui: &mut egui::Ui, fr: &Frame) {
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.bottom, Bottom::History, "Historique");
            if fr.wall.is_some() {
                ui.selectable_value(&mut self.bottom, Bottom::Cp, "Pression pariétale Cp");
                ui.selectable_value(&mut self.bottom, Bottom::Cf, "Frottement Cf");
            }
            ui.selectable_value(&mut self.bottom, Bottom::Polar, format!("Polaire ({})", self.sweep.points.len()));
        });
        match self.bottom {
            Bottom::Polar => self.polar_chart(ui),
            Bottom::Cp | Bottom::Cf if fr.wall.is_some() => self.wall_chart(ui, fr),
            _ => self.history(ui, fr),
        }
    }

    fn history(&self, ui: &mut egui::Ui, fr: &Frame) {
        let cx: Vec<[f64; 2]> = fr.history.iter().map(|h| [h[0], h[1]]).collect();
        let cz: Vec<[f64; 2]> = fr.history.iter().map(|h| [h[0], h[2]]).collect();
        Plot::new("hist")
            .height(ui.available_height().max(80.0))
            .legend(Legend::default())
            .allow_scroll(false)
            .show(ui, |p| {
                p.line(Line::new("Cx", PlotPoints::from(cx)).color(theme::accent()));
                p.line(Line::new("Cz", PlotPoints::from(cz)).color(theme::lift()));
            });
    }

    /// −Cp (aspiration vers le haut) ou Cf×1000 le long de la paroi : points
    /// LBM, tirets de la méthode des panneaux (ou écoulement potentiel).
    fn wall_chart(&self, ui: &mut egui::Ui, fr: &Frame) {
        let Some(w) = &fr.wall else { return };
        let cp = self.bottom == Bottom::Cp;
        let cyl = w.axis == WallAxis::Angle;
        let (up_name, lo_name) = if cyl {
            ("dessus", "dessous")
        } else {
            ("extrados", "intrados")
        };
        let val = |p: &lbm::wall::WallPoint| if cp { -p.cp } else { 1000.0 * p.cf };
        let up: Vec<[f64; 2]> = w.upper.iter().map(|p| [p.s, val(p)]).collect();
        let lo: Vec<[f64; 2]> = w.lower.iter().map(|p| [p.s, val(p)]).collect();
        let dash = LineStyle::dashed_dense();
        let mut refs: Vec<(String, Vec<[f64; 2]>, Color32)> = Vec::new();
        if cyl {
            if cp {
                let pts = (0..=90)
                    .map(|k| {
                        let t = 2.0 * k as f64;
                        [t, -(1.0 - 4.0 * t.to_radians().sin().powi(2))]
                    })
                    .collect();
                refs.push(("potentiel".into(), pts, Color32::GRAY));
            }
        } else if let Some(r) = &fr.panel {
            if cp {
                let (mut pu, mut pl) = (Vec::new(), Vec::new());
                for i in 0..r.cp_x.len() {
                    if i >= r.i_le {
                        pu.push([r.cp_x[i], -r.cp[i]]);
                    }
                    if i <= r.i_le {
                        pl.push([r.cp_x[i], -r.cp[i]]);
                    }
                }
                refs.push((format!("panneaux, {up_name}"), pu, theme::lift()));
                refs.push((format!("panneaux, {lo_name}"), pl, theme::accent()));
            } else {
                let side = |f: &lbm::panel::Friction| {
                    (0..f.x.len())
                        .map(|i| [f.x[i], 1000.0 * f.cf[i] * f.ue[i] * f.ue[i]])
                        .collect::<Vec<_>>()
                };
                refs.push((format!("panneaux, {up_name}"), side(&r.friction_upper), theme::lift()));
                refs.push((format!("panneaux, {lo_name}"), side(&r.friction_lower), theme::accent()));
            }
        }
        let marks: Vec<(f64, Color32)> = match &fr.panel {
            Some(r) if !cyl && fr.sa && self.settings.cfg.transition == Transition::Free => {
                [(r.transition_upper, theme::lift()), (r.transition_lower, theme::accent())]
                    .into_iter()
                    .filter(|m| m.0.is_finite())
                    .collect()
            }
            _ => Vec::new(),
        };
        let x_name = if cyl { "angle depuis l'arrêt (°)" } else { "x/c" };
        let y_name = if cp { "−Cp" } else { "Cf × 1000" };
        ui.label(
            egui::RichText::new(if w.n > 0 {
                format!(
                    "{y_name} selon {x_name} · points : mailles voisines de la paroi, moyenne de {} relevés{}",
                    w.n,
                    if cp {
                        " ; pic d'aspiration du bord d'attaque écrêté"
                    } else if fr.sa {
                        " ; frottement de la loi de paroi (laminaire en amont de la transition)"
                    } else {
                        ""
                    }
                )
            } else {
                "Moyenne disponible après le transitoire".into()
            })
            .small()
            .weak(),
        );
        Plot::new(if cp { "cp" } else { "cf" })
            .height(ui.available_height().max(80.0))
            .legend(Legend::default())
            .allow_scroll(false)
            .x_axis_label(x_name)
            .y_axis_label(y_name)
            .show(ui, |p| {
                for (name, pts, color) in refs {
                    p.line(Line::new(name, PlotPoints::from(pts)).color(color).style(dash));
                }
                p.points(
                    Points::new(format!("LBM, {up_name}"), PlotPoints::from(up))
                        .color(theme::lift())
                        .radius(1.8),
                );
                p.points(
                    Points::new(format!("LBM, {lo_name}"), PlotPoints::from(lo))
                        .color(theme::accent())
                        .radius(1.8),
                );
                for (x, c) in marks {
                    p.vline(VLine::new("transition", x).color(c).style(LineStyle::dotted_dense()));
                }
            });
    }

    /// Veine d'essai : image du champ, traceurs, maquette, force ; zoom et déplacement.
    fn tunnel(&mut self, ui: &mut egui::Ui, fr: &Frame) {
        let avail = ui.available_size();
        let (resp, painter) = ui.allocate_painter(avail, Sense::click_and_drag());
        let rect = resp.rect;
        painter.rect_filled(rect, 0.0, theme::canvas());
        if fr.nx == 0 {
            painter.text(
                rect.center(),
                egui::Align2::CENTER_CENTER,
                "Calcul en cours…",
                egui::FontId::proportional(16.0),
                Color32::GRAY,
            );
            return;
        }
        if fr.image_id != self.image_id {
            self.image_id = fr.image_id;
            let img = egui::ColorImage::from_rgba_unmultiplied([fr.nx, fr.ny], &fr.rgba);
            match &mut self.texture {
                Some(t) => t.set(img, egui::TextureOptions::LINEAR),
                None => self.texture = Some(ui.ctx().load_texture("champ", img, egui::TextureOptions::LINEAR)),
            }
        }
        // Cadrage par défaut : 2 L en amont, 6 L en aval de la maquette.
        let l = fr.l as f32;
        let default_view = (Pos2::new(fr.refp[0] as f32 + 2.0 * l, fr.refp[1] as f32), 8.0 * l);
        let (mut c, mut w) = self.view.unwrap_or(default_view);
        // Zoom (molette, pincement) autour du pointeur ; déplacement au glisser.
        let s0 = rect.width() / w;
        let to_lat =
            |p: Pos2, c: Pos2, s: f32| Pos2::new(c.x + (p.x - rect.center().x) / s, c.y - (p.y - rect.center().y) / s);
        if resp.hovered() {
            let (scroll, zoom) = ui.input(|i| (i.smooth_scroll_delta.y, i.zoom_delta()));
            let factor = zoom * (scroll * 0.0015).exp();
            if (factor - 1.0).abs() > 1e-4 {
                if let Some(p) = resp.hover_pos() {
                    let a = to_lat(p, c, s0);
                    w = (w / factor).clamp(0.4 * l, fr.nx as f32);
                    let s1 = rect.width() / w;
                    c = Pos2::new(a.x - (p.x - rect.center().x) / s1, a.y + (p.y - rect.center().y) / s1);
                    self.view = Some((c, w));
                }
            }
        }
        if resp.dragged() {
            let d = resp.drag_delta();
            c = Pos2::new(c.x - d.x / s0, c.y + d.y / s0);
            self.view = Some((c, w));
        }
        if resp.double_clicked() {
            self.view = None;
            (c, w) = default_view;
        }
        let s = rect.width() / w;
        let to_screen = |x: f32, y: f32| {
            Pos2::new(
                rect.center().x + (x + 0.5 - c.x) * s,
                rect.center().y - (y + 0.5 - c.y) * s,
            )
        };
        let painter = painter.with_clip_rect(rect);
        if let Some(t) = &self.texture {
            // Image : maille (i, j) → texel (i, ny − 1 − j).
            let p0 = to_screen(-0.5, fr.ny as f32 - 0.5);
            let p1 = to_screen(fr.nx as f32 - 0.5, -0.5);
            painter.image(
                t.id(),
                Rect::from_two_pos(p0, p1),
                Rect::from_min_max(Pos2::ZERO, Pos2::new(1.0, 1.0)),
                Color32::WHITE,
            );
        }
        // Traceurs.
        let lw = 1.3;
        for (pl, &odd) in fr.lines.iter().zip(&fr.line_parity) {
            let pts: Vec<Pos2> = pl.iter().map(|p| to_screen(p[0], p[1])).collect();
            painter.add(Shape::line(pts, Stroke::new(lw, if odd { theme::accent() } else { theme::tracer() })));
        }
        for p in &fr.points {
            let q = to_screen(p[0], p[1]);
            if rect.contains(q) {
                painter.rect_filled(Rect::from_center_size(q, Vec2::splat(1.6)), 0.0, theme::tracer());
            }
        }
        // Maquette : éventail depuis le point de référence (formes étoilées).
        let body = theme::body();
        let r = to_screen(fr.refp[0] as f32, fr.refp[1] as f32);
        let ol: Vec<Pos2> = fr.outline.iter().map(|p| to_screen(p[0] as f32, p[1] as f32)).collect();
        let mut mesh = egui::Mesh::default();
        mesh.colored_vertex(r, body);
        for p in &ol {
            mesh.colored_vertex(*p, body);
        }
        let m = ol.len() as u32;
        for k in 0..m {
            mesh.add_triangle(0, 1 + k, 1 + (k + 1) % m);
        }
        painter.add(Shape::mesh(mesh));
        let mut closed = ol.clone();
        if let Some(f) = ol.first() {
            closed.push(*f);
        }
        painter.add(Shape::line(
            closed,
            Stroke::new(1.0, theme::body_edge()),
        ));
        painter.circle_filled(r, 2.5, theme::body_edge());
        // Force aérodynamique moyenne : 1 unité de coefficient ↔ 1,2 L.
        let a = &fr.analysis;
        if self.show_force && fr.ready && a.cx.is_finite() && a.cz.is_finite() {
            let k = 1.2 * l;
            let tip = to_screen(fr.refp[0] as f32 + a.cx as f32 * k, fr.refp[1] as f32 + a.cz as f32 * k);
            painter.line_segment([r, Pos2::new(tip.x, r.y)], Stroke::new(1.2, theme::accent()));
            painter.line_segment([r, Pos2::new(r.x, tip.y)], Stroke::new(1.2, theme::lift()));
            painter.arrow(r, tip - r, Stroke::new(2.0, theme::ink()));
        }
        // Légende du champ.
        let (lo, hi, _) = worker::field_range(self.settings.field);
        let name = match self.settings.field {
            Field::Vorticity => "Vorticité ω·L/U",
            Field::Speed => "Vitesse |u|/U",
            Field::Pressure => "Pression Cp",
            Field::Eddy => "Viscosité turbulente log₁₀(1 + νt/ν)",
        };
        painter.text(
            rect.left_bottom() + Vec2::new(8.0, -8.0),
            egui::Align2::LEFT_BOTTOM,
            format!("{name} : de {lo} à {hi}   ·   {}   ·   t* = {:.1}", fr.label, fr.t_star),
            egui::FontId::proportional(13.0),
            theme::ink(),
        );
    }
}

impl eframe::App for App {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        theme::sync(ui.ctx());
        let mut tab = self.tab;
        egui::Panel::top("onglets").show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.selectable_value(&mut tab, Tab::Tunnel, "Soufflerie 2D");
                ui.selectable_value(&mut tab, Tab::Wing, "Aile finie");
                ui.selectable_value(&mut tab, Tab::Tunnel3D, "Soufflerie 3D");
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), theme::switch);
            });
        });
        let ctx = ui.ctx().clone();
        self.set_tab(tab, &ctx);
        match self.tab {
            Tab::Wing => return self.wing.ui(ui),
            Tab::Tunnel3D => return self.tab3d.ui(ui),
            Tab::Tunnel => {}
        }
        let fr = self.shared.lock().unwrap().frame.clone();
        egui::Panel::left("reglages")
            .default_size(340.0)
            .resizable(true)
            .show(ui, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| self.controls(ui, &fr));
            });
        egui::Panel::bottom("historique")
            .default_size(190.0)
            .resizable(true)
            .show(ui, |ui| self.bottom(ui, &fr));
        egui::CentralPanel::default().show(ui, |ui| self.tunnel(ui, &fr));
    }
}

/// Reynolds lisible : 20 000, 2·10⁵, 10⁶…
pub(crate) fn fmt_re(v: f64) -> String {
    if v < 1e5 {
        let r = if v < 100.0 {
            v.round()
        } else {
            (v / 10.0).round() * 10.0
        };
        let s = format!("{r:.0}");
        let n = s.len();
        s.chars()
            .enumerate()
            .flat_map(|(k, c)| {
                let sep = k > 0 && (n - k) % 3 == 0;
                sep.then_some('\u{202f}').into_iter().chain(std::iter::once(c))
            })
            .collect()
    } else {
        let e = (v.log10() + 1e-9).floor() as i32;
        let m = (v / 10f64.powi(e) * 10.0).round() / 10.0;
        let sup: String = e
            .to_string()
            .chars()
            .map(|d| "⁰¹²³⁴⁵⁶⁷⁸⁹".chars().nth(d as usize - '0' as usize).unwrap())
            .collect();
        if (m - 1.0).abs() < 1e-9 {
            format!("10{sup}")
        } else {
            format!("{}·10{sup}", format!("{m}").replace('.', ","))
        }
    }
}
