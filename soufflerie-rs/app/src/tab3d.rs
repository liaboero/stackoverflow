//! Onglet « Soufflerie 3D » : calcul D3Q19 sur GPU dans un fil dédié, plan de
//! coupe (vitesse, vorticité normale au plan, critère Q), efforts et
//! historique. Portage de js/lab3d.js ; le rendu 3D par lancer de rayons de
//! l'application web est remplacé par une coupe plane.

use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use eframe::egui::{self, Color32, Pos2, Rect, Sense, Vec2};
use egui_plot::{Legend, Line, Plot, PlotPoints};

use crate::theme;

fn rgbf(c: Color32) -> [f32; 3] {
    [c.r() as f32 / 255.0, c.g() as f32 / 255.0, c.b() as f32 / 255.0]
}
use lbm::bodies3d::{Body3D, Shape3D};
use lbm::fluid::Fluid;
use lbm::stl::{Axis, MeshInfo, MeshSetup, Tri, Unit};
use lbm::cpu::Ramp;
use lbm::gpu3d::{bytes_per_cell, Camera, Gpu3D, RenderOptions, Slice, SliceField};
use lbm::stats::{analyze, Analysis};

/// Reynolds proposés ; Spalart–Allmaras et loi de paroi dès 5·10⁴.
const RE3: [f64; 14] = [
    10.0, 20.0, 50.0, 100.0, 200.0, 300.0, 500.0, 1000.0, 2000.0, 5000.0, 1e4, 1e5, 1e6, 1e7,
];
const RES3: [(usize, &str); 5] = [
    (8, "Aperçu"),
    (12, "Rapide"),
    (16, "Standard"),
    (24, "Fine"),
    (32, "Très fine"),
];
const U: f32 = 0.1;
const RAMP_T: f64 = 4.0;
const SETTLE_T: f64 = 6.0;
const WINDOW_T: f64 = 20.0;

/// Plan de coupe.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Plane {
    /// Longitudinal (x, z) à y donné.
    Y,
    /// Horizontal (x, y) à z donné.
    Z,
    /// Transversal (y, z) à x donné.
    X,
}

impl Plane {
    fn axis(self) -> usize {
        match self {
            Plane::X => 0,
            Plane::Y => 1,
            Plane::Z => 2,
        }
    }
}

/// Maquette STL chargée (comparée par numéro de chargement).
#[derive(Clone)]
struct StlFile {
    name: String,
    tris: Arc<Vec<Tri>>,
    id: u64,
}

impl PartialEq for StlFile {
    fn eq(&self, o: &StlFile) -> bool {
        self.id == o.id
    }
}

/// Résolution de l'aperçu de mise en place (les dimensions de la veine
/// sont proportionnelles à L).
const PREVIEW_L: usize = 32;

type Preview = Result<([usize; 3], MeshInfo), String>;

/// Réglages transmis au fil de calcul.
#[derive(Clone, PartialEq)]
struct Settings {
    shape: Shape3D,
    /// Maquette STL et sa mise en place (incidence : `alpha`).
    stl: Option<StlFile>,
    setup: MeshSetup,
    /// Échelle physique de la maquette STL : le Reynolds découle de la vitesse.
    phys: bool,
    speed: f64,
    fluid: Fluid,
    re: f64,
    l: usize,
    alpha: f64,
    fp16: bool,
    running: bool,
    field: SliceField,
    plane: Plane,
    /// Position relative du plan (0 à 1) ; None : maquette ou sillage.
    pos: Option<f64>,
    /// Incrémenté pour relancer le calcul.
    reset: u64,
    /// Vue 3D (maquette, plan de coupe, isosurfaces) ou coupe plane.
    view3d: bool,
    /// Caméra : azimut, élévation (°), distance (en L).
    az: f64,
    el: f64,
    dist: f64,
    iso: bool,
    /// Seuil des isosurfaces, en (U/L)².
    iso_level: f64,
    /// Taille de l'image 3D (pixels).
    img: (usize, usize),
}

/// Ce que le fil de calcul publie.
#[derive(Default)]
struct Frame {
    slice: Option<Arc<Slice>>,
    /// Vue 3D : largeur, hauteur, RGBA, et la caméra qui l'a produite.
    image: Option<Arc<(usize, usize, Vec<u8>)>>,
    camera: Option<Camera>,
    slice_id: u64,
    /// Plan effectivement affiché (axe, maille).
    plane_pos: (usize, usize),
    dims: [usize; 3],
    center: [f64; 3],
    l: usize,
    t_star: f64,
    /// t*, Cx, Cz, Cy.
    hist: Vec<[f64; 4]>,
    analysis: Analysis,
    ready: bool,
    status: String,
    info: String,
    error: Option<String>,
    /// Avertissement sur la maquette (parties perdues à la voxelisation).
    warning: Option<String>,
    mlups: f64,
}

struct Shared {
    settings: Settings,
    frame: Frame,
    /// Fil de calcul actif (onglet affiché).
    active: bool,
}

pub struct Tab3D {
    shared: Arc<Mutex<Shared>>,
    settings: Settings,
    /// Chargements STL successifs (identifiants).
    stl_count: u64,
    stl_error: Option<String>,
    /// Mise en place calculée pour l'aperçu : (fichier, réglage, incidence) →
    /// dimensions de la veine à L = 32 et grandeurs de référence.
    preview: Option<((u64, MeshSetup, f64), Preview)>,
    /// Réglage en attente de stabilisation (gros maillages) et depuis quand.
    pending: Option<((u64, MeshSetup, f64), Instant)>,
    texture: Option<egui::TextureHandle>,
    slice_id: u64,
    started: bool,
}

impl Tab3D {
    pub fn new() -> Tab3D {
        let settings = Settings {
            shape: Shape3D::Sphere,
            stl: None,
            setup: MeshSetup::default(),
            phys: false,
            speed: 20.0,
            fluid: Fluid::Air,
            re: 300.0,
            l: 16,
            alpha: 8.0,
            fp16: true,
            running: true,
            field: SliceField::Vorticity,
            plane: Plane::Y,
            pos: None,
            reset: 0,
            view3d: true,
            az: 62.0,
            el: 24.0,
            dist: 6.0,
            iso: true,
            iso_level: 0.4,
            img: (800, 450),
        };
        Tab3D {
            shared: Arc::new(Mutex::new(Shared {
                settings: settings.clone(),
                frame: Frame::default(),
                active: false,
            })),
            settings,
            stl_count: 0,
            stl_error: None,
            preview: None,
            pending: None,
            texture: None,
            slice_id: 0,
            started: false,
        }
    }

    /// Onglet affiché ou non (le calcul s'arrête quand il est masqué).
    pub fn set_active(&mut self, on: bool, ctx: &egui::Context) {
        self.shared.lock().unwrap().active = on;
        if on && !self.started {
            self.started = true;
            let shared = self.shared.clone();
            let ctx = ctx.clone();
            std::thread::spawn(move || worker(shared, move || ctx.request_repaint()));
        }
    }

    fn push(&self) {
        self.shared.lock().unwrap().settings = self.settings.clone();
    }

    /// Lit un fichier STL et en fait la maquette courante.
    pub fn load_stl(&mut self, path: &std::path::Path) {
        match lbm::stl::read_stl(path) {
            Ok(tris) => {
                self.stl_count += 1;
                self.settings.stl = Some(StlFile {
                    name: path
                        .file_name()
                        .map_or_else(|| path.display().to_string(), |n| n.to_string_lossy().into()),
                    tris: Arc::new(tris),
                    id: self.stl_count,
                });
                self.settings.shape = Shape3D::Mesh;
                self.settings.pos = None;
                self.stl_error = None;
                self.push();
            }
            Err(e) => self.stl_error = Some(e),
        }
    }

    /// Mise en place de la maquette STL, recalculée quand le réglage change ;
    /// pour un gros maillage, une fois le réglage stable depuis 0,25 s (le
    /// curseur reste fluide, l'aperçu précédent reste affiché entre-temps).
    fn mesh_preview(&mut self, ctx: &egui::Context) -> Option<Preview> {
        let s = &self.settings;
        let f = s.stl.as_ref()?;
        let key = (f.id, s.setup.clone(), s.alpha);
        let stale = self.preview.as_ref().is_none_or(|(k, _)| *k != key);
        let mut due = self.preview.is_none() || f.tris.len() < 20_000;
        if stale && !due {
            if self.pending.as_ref().is_none_or(|(k, _)| *k != key) {
                self.pending = Some((key.clone(), Instant::now()));
            }
            let wait = Duration::from_millis(250);
            let since = self.pending.as_ref().map_or(wait, |(_, t)| t.elapsed());
            due = since >= wait;
            if !due {
                ctx.request_repaint_after(wait - since);
            }
        }
        if stale && due {
            self.pending = None;
            let setup = MeshSetup {
                alpha_deg: s.alpha,
                ..s.setup.clone()
            };
            let r = Body3D::from_mesh(&f.tris, PREVIEW_L, &setup).map(|(b, i)| ([b.nx, b.ny, b.nz], i));
            self.preview = Some((key, r));
        }
        self.preview.as_ref().map(|(_, r)| r.clone())
    }

    pub fn ui(&mut self, ui: &mut egui::Ui) {
        let dropped: Vec<std::path::PathBuf> =
            ui.ctx().input(|i| i.raw.dropped_files.iter().map(|f| f.path().to_path_buf()).collect());
        if let Some(p) = dropped
            .iter()
            .find(|p| p.extension().is_some_and(|e| e.eq_ignore_ascii_case("stl")))
        {
            self.load_stl(p);
        }
        egui::Panel::left("s3_reglages")
            .default_size(340.0)
            .resizable(true)
            .show(ui, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| self.controls(ui));
            });
        egui::Panel::bottom("s3_historique")
            .default_size(180.0)
            .resizable(true)
            .show(ui, |ui| self.history(ui));
        egui::CentralPanel::default().show(ui, |ui| self.view(ui));
    }

    fn controls(&mut self, ui: &mut egui::Ui) {
        let before = self.settings.clone();
        ui.heading("Maquette");
        let mut pick = false;
        {
            let s = &mut self.settings;
            ui.horizontal_wrapped(|ui| {
                for k in Shape3D::ALL {
                    if ui.selectable_value(&mut s.shape, k, k.label()).changed() {
                        s.pos = None;
                    }
                }
                if let Some(f) = &s.stl {
                    if ui
                        .selectable_value(&mut s.shape, Shape3D::Mesh, format!("STL : {}", f.name))
                        .changed()
                    {
                        s.pos = None;
                    }
                }
                pick = ui
                    .button("Charger un STL…")
                    .on_hover_text(
                        "Maquette maillée (STL binaire ou ASCII, volume fermé) ; on peut aussi glisser le fichier sur la fenêtre",
                    )
                    .clicked();
            });
        }
        if pick {
            if let Some(p) = rfd::FileDialog::new()
                .add_filter("Maillage STL", &["stl", "STL"])
                .pick_file()
            {
                self.load_stl(&p);
            }
        }
        if let Some(e) = &self.stl_error {
            ui.colored_label(Color32::RED, e);
        }
        let mesh = self.settings.shape == Shape3D::Mesh && self.settings.stl.is_some();
        let preview = if mesh { self.mesh_preview(ui.ctx()) } else { None };
        let s = &mut self.settings;
        if s.shape.has_alpha() {
            ui.add(egui::Slider::new(&mut s.alpha, -10.0..=20.0).step_by(1.0).text("incidence (°)"));
        }
        if mesh {
            mesh_controls(ui, s, preview.as_ref());
        }
        if mesh && s.phys {
            if let Some(Ok((_, info))) = &preview {
                let (nu, _) = s.fluid.props();
                s.re = (s.speed * info.ref_length_m / nu).max(1.0);
            }
            ui.label(format!("Reynolds {:.0} (vitesse et longueur de référence)", s.re));
            if s.re > 5000.0 {
                ui.label(
                    egui::RichText::new(
                        "Au-delà de quelques milliers, la couche limite n'est pas résolue (sous-maille seule) : coefficients indicatifs.",
                    )
                    .small()
                    .color(theme::warn()),
                );
            }
        } else {
            ui.horizontal(|ui| {
                ui.label("Reynolds");
                let mut k = RE3.iter().position(|&r| r == s.re).unwrap_or(5) as f64;
                ui.add(egui::Slider::new(&mut k, 0.0..=(RE3.len() - 1) as f64).step_by(1.0).show_value(false));
                s.re = RE3[k.round() as usize];
                ui.label(if s.re >= 1e5 { format!("{:.0e}", s.re) } else { format!("{}", s.re) });
                if s.re >= lbm::tunnel::RE_SA_AUTO {
                    ui.label(egui::RichText::new("loi de paroi").small().weak())
                        .on_hover_text("Spalart–Allmaras et loi de Spalding : couche limite turbulente supposée dès l'avant");
                }
            });
        }
        ui.label("Résolution (mailles par longueur de référence)");
        ui.horizontal_wrapped(|ui| {
            for (l, name) in RES3 {
                let [nx, ny, nz] = match (&preview, mesh) {
                    (Some(Ok((d, _))), true) => d.map(|v| v * l / PREVIEW_L),
                    _ => s.shape.dims(l),
                };
                let n = nx * ny * nz;
                ui.selectable_value(&mut s.l, l, format!("{name} · {l}")).on_hover_text(format!(
                    "{:.2} M mailles, ≈ {} Mo de mémoire GPU",
                    n as f64 / 1e6,
                    n * bytes_per_cell(s.fp16) / 1_000_000
                ));
            }
        });
        ui.horizontal(|ui| {
            ui.label("Précision");
            ui.selectable_value(&mut s.fp16, true, "demi (2 octets)");
            ui.selectable_value(&mut s.fp16, false, "simple");
        });
        ui.horizontal(|ui| {
            let label = if s.running { "Pause" } else { "Reprendre" };
            if ui.button(label).clicked() {
                s.running = !s.running;
            }
            if ui.button("Recommencer").clicked() {
                s.reset += 1;
            }
        });
        ui.separator();
        ui.heading("Affichage");
        ui.horizontal(|ui| {
            ui.selectable_value(&mut s.view3d, true, "Vue 3D");
            ui.selectable_value(&mut s.view3d, false, "Coupe à plat");
        });
        if s.view3d {
            ui.horizontal(|ui| {
                ui.checkbox(&mut s.iso, "Isosurfaces du critère Q");
            });
            if s.iso {
                ui.add(
                    egui::Slider::new(&mut s.iso_level, 0.05..=5.0)
                        .logarithmic(true)
                        .text("seuil Q·(L/U)²"),
                );
            }
            ui.label(
                egui::RichText::new("Glisser : tourner · molette : rapprocher · double clic : vue par défaut")
                    .small()
                    .weak(),
            );
        }
        ui.label("Plan de coupe");
        ui.horizontal(|ui| {
            ui.selectable_value(&mut s.field, SliceField::Vorticity, "Vorticité");
            ui.selectable_value(&mut s.field, SliceField::Speed, "Vitesse");
            if !s.view3d {
                ui.selectable_value(&mut s.field, SliceField::Q, "Critère Q");
            }
        });
        ui.horizontal(|ui| {
            for (p, name) in [(Plane::Y, "Longitudinal"), (Plane::Z, "Horizontal"), (Plane::X, "Transversal")] {
                if ui.selectable_value(&mut s.plane, p, name).changed() {
                    s.pos = None;
                }
            }
        });
        let fr = self.shared.lock().unwrap();
        let f = &fr.frame;
        let (ax, cell) = f.plane_pos;
        let n = f.dims[ax].max(2);
        let mut rel = s.pos.unwrap_or(cell as f64 / (n - 1) as f64);
        let off = if f.l > 0 {
            (cell as f64 - f.center[ax]) / f.l as f64
        } else {
            0.0
        };
        if ui
            .add(egui::Slider::new(&mut rel, 0.0..=1.0).show_value(false).text(format!("position {off:+.1} L")))
            .changed()
        {
            s.pos = Some(rel);
        }
        ui.separator();
        ui.heading("Efforts");
        ui.horizontal(|ui| {
            ui.label(egui::RichText::new(&f.status).strong());
            ui.label(format!("t* = {:.1}", f.t_star));
        });
        let a = &f.analysis;
        let last = f.hist.last().copied().unwrap_or([0.0, f64::NAN, f64::NAN, f64::NAN]);
        let (cx, cz, cy) = if f.ready { (a.cx, a.cz, a.cm) } else { (last[1], last[2], last[3]) };
        egui::Grid::new("s3_coef").num_columns(2).spacing([16.0, 2.0]).show(ui, |ui| {
            let v = |x: f64| if x.is_finite() { format!("{x:.4}") } else { "—".into() };
            ui.colored_label(theme::accent(), "Cx (traînée)");
            ui.label(if f.ready && a.periodic && a.cx_amp > 0.001 {
                format!("{} ± {:.4}", v(cx), a.cx_amp)
            } else {
                v(cx)
            });
            ui.end_row();
            ui.colored_label(theme::lift(), "Cz (portance)");
            ui.label(if f.ready && a.periodic && a.cz_amp > 0.001 {
                format!("{} ± {:.4}", v(cz), a.cz_amp)
            } else {
                v(cz)
            });
            ui.end_row();
            ui.label("Cy (latéral)");
            ui.label(v(cy));
            ui.end_row();
            if let (true, true, Some(Ok((_, info)))) = (mesh, s.phys, &preview) {
                let (_, rho) = s.fluid.props();
                let qs = 0.5 * rho * s.speed * s.speed * info.ref_area_m2;
                let n = |c: f64| if c.is_finite() { format!("{} N", sig(c * qs)) } else { "—".into() };
                for (name, c) in [("Traînée", cx), ("Portance", cz), ("Effort latéral", cy)] {
                    ui.label(format!("{name} à {} m/s", s.speed));
                    ui.label(n(c));
                    ui.end_row();
                }
            }
            ui.label("Strouhal");
            ui.label(if f.ready && a.periodic && a.st.is_finite() {
                format!("{:.3} ({} périodes)", a.st, a.periods)
            } else {
                "—".into()
            });
            ui.end_row();
        });
        ui.label(
            egui::RichText::new(if f.ready { "moyenne sur la fenêtre établie" } else { "valeurs instantanées" })
                .small()
                .weak(),
        );
        ui.add_space(4.0);
        ui.label(egui::RichText::new(reference(s.shape, s.re)).small());
        ui.add_space(4.0);
        ui.label(egui::RichText::new(&f.info).small().weak());
        if let Some(w) = &f.warning {
            ui.colored_label(theme::warn(), w);
        }
        if let Some(e) = &f.error {
            ui.colored_label(Color32::RED, e);
        }
        drop(fr);
        if *s != before {
            self.push();
        }
    }

    fn view(&mut self, ui: &mut egui::Ui) {
        if self.settings.view3d {
            return self.view3d(ui);
        }
        let (resp, painter) = ui.allocate_painter(ui.available_size(), Sense::hover());
        let rect = resp.rect;
        painter.rect_filled(rect, 0.0, theme::canvas());
        let (slice, id, plane_pos, l) = {
            let s = self.shared.lock().unwrap();
            (s.frame.slice.clone(), s.frame.slice_id, s.frame.plane_pos, s.frame.l)
        };
        let Some(sl) = slice else {
            painter.text(
                rect.center(),
                egui::Align2::CENTER_CENTER,
                "Préparation du calcul 3D…",
                egui::FontId::proportional(16.0),
                Color32::GRAY,
            );
            return;
        };
        let field = self.settings.field;
        if id != self.slice_id || self.texture.is_none() {
            self.slice_id = id;
            let (lo, hi, div) = match field {
                SliceField::Speed => (0.0, 1.6, false),
                SliceField::Vorticity => (-4.0, 4.0, true),
                SliceField::Q => (-1.0, 1.0, true),
            };
            let mut px = vec![0u8; 4 * sl.w * sl.h];
            for j in 0..sl.h {
                let row = sl.h - 1 - j;
                for i in 0..sl.w {
                    let v = sl.data[i + sl.w * j];
                    let c = if v > 1e29 {
                        { let c = theme::body_rgb(); [c[0], c[1], c[2], 255] }
                    } else {
                        let t = (v - lo) / (hi - lo);
                        theme::color(if div { theme::diverging() } else { theme::sequential() }, t)
                    };
                    let o = 4 * (i + sl.w * row);
                    px[o..o + 4].copy_from_slice(&c);
                }
            }
            let img = egui::ColorImage::from_rgba_unmultiplied([sl.w, sl.h], &px);
            match &mut self.texture {
                Some(t) => t.set(img, egui::TextureOptions::LINEAR),
                None => self.texture = Some(ui.ctx().load_texture("coupe3d", img, egui::TextureOptions::LINEAR)),
            }
        }
        // Image à l'échelle, centrée (mailles carrées).
        let margin = 12.0;
        let avail = rect.shrink(margin);
        let sc = (avail.width() / sl.w as f32).min((avail.height() - 24.0) / sl.h as f32);
        let size = Vec2::new(sl.w as f32 * sc, sl.h as f32 * sc);
        let img_rect = Rect::from_center_size(avail.center() - Vec2::new(0.0, 10.0), size);
        if let Some(t) = &self.texture {
            painter.image(
                t.id(),
                img_rect,
                Rect::from_min_max(Pos2::ZERO, Pos2::new(1.0, 1.0)),
                Color32::WHITE,
            );
        }
        painter.rect_stroke(img_rect, 0.0, (1.0, theme::ink()), egui::StrokeKind::Outside);
        let (h_axis, v_axis) = match plane_pos.0 {
            0 => ("y", "z"),
            1 => ("x", "z"),
            _ => ("x", "y"),
        };
        let comp = ["ωx", "ωy", "ωz"][plane_pos.0];
        let legend = match field {
            SliceField::Speed => "|u|/U de 0 à 1,6".to_string(),
            SliceField::Vorticity => format!("{comp}·L/U de −4 à +4"),
            SliceField::Q => "critère Q·(L/U)² de −1 à +1 (positif : tourbillons)".to_string(),
        };
        painter.text(
            rect.left_bottom() + Vec2::new(8.0, -8.0),
            egui::Align2::LEFT_BOTTOM,
            format!(
                "Coupe ({h_axis}, {v_axis}) en {} = {} · {legend} · L = {l} mailles · écoulement de gauche à droite{}",
                ["x", "y", "z"][plane_pos.0],
                plane_pos.1,
                if plane_pos.0 == 0 { " (vu de l'amont)" } else { "" }
            ),
            egui::FontId::proportional(13.0),
            theme::ink(),
        );
    }

    /// Vue 3D : image calculée sur le GPU, cadre de la veine, caméra orbitale.
    fn view3d(&mut self, ui: &mut egui::Ui) {
        let (resp, painter) = ui.allocate_painter(ui.available_size(), Sense::click_and_drag());
        let rect = resp.rect;
        painter.rect_filled(rect, 0.0, theme::canvas());
        // Caméra.
        let mut changed = false;
        if resp.dragged() {
            let d = resp.drag_delta();
            self.settings.az -= d.x as f64 * 0.4;
            self.settings.el = (self.settings.el + d.y as f64 * 0.3).clamp(-85.0, 85.0);
            changed = true;
        }
        if resp.hovered() {
            let (scroll, zoom) = ui.input(|i| (i.smooth_scroll_delta.y, i.zoom_delta()));
            let f = (-(scroll as f64) * 0.0015).exp() / zoom as f64;
            if (f - 1.0).abs() > 1e-4 {
                self.settings.dist = (self.settings.dist * f).clamp(1.5, 20.0);
                changed = true;
            }
        }
        if resp.double_clicked() {
            (self.settings.az, self.settings.el, self.settings.dist) = (62.0, 24.0, 6.0);
            changed = true;
        }
        // Image à ~900 pixels de large au plus (le lancer de rayons coûte par
        // pixel), à demi-résolution pendant qu'on tourne ou rapproche la vue.
        let moving = resp.dragged() || changed;
        let k = (900.0 / rect.width()).min(1.0) * if moving { 0.5 } else { 1.0 };
        let want = (
            ((rect.width() * k) as usize).max(64),
            ((rect.height() * k) as usize).max(64),
        );
        if want != self.settings.img {
            self.settings.img = want;
            changed = true;
        }
        if changed {
            self.push();
        }
        let (image, id, cam, dims, l, plane_pos) = {
            let s = self.shared.lock().unwrap();
            let f = &s.frame;
            (f.image.clone(), f.slice_id, f.camera, f.dims, f.l, f.plane_pos)
        };
        let Some(img) = image else {
            painter.text(
                rect.center(),
                egui::Align2::CENTER_CENTER,
                "Préparation du calcul 3D…",
                egui::FontId::proportional(16.0),
                Color32::GRAY,
            );
            return;
        };
        let (w, h, px) = (&img.0, &img.1, &img.2);
        if id != self.slice_id || self.texture.is_none() {
            self.slice_id = id;
            let ci = egui::ColorImage::from_rgba_unmultiplied([*w, *h], px);
            match &mut self.texture {
                Some(t) => t.set(ci, egui::TextureOptions::LINEAR),
                None => self.texture = Some(ui.ctx().load_texture("vue3d", ci, egui::TextureOptions::LINEAR)),
            }
        }
        // L'image garde le rapport de l'instant où elle a été calculée.
        let aspect = *w as f32 / *h as f32;
        let size = if rect.width() / rect.height() > aspect {
            Vec2::new(rect.height() * aspect, rect.height())
        } else {
            Vec2::new(rect.width(), rect.width() / aspect)
        };
        let img_rect = Rect::from_center_size(rect.center(), size);
        if let Some(t) = &self.texture {
            painter.image(
                t.id(),
                img_rect,
                Rect::from_min_max(Pos2::ZERO, Pos2::new(1.0, 1.0)),
                Color32::WHITE,
            );
        }
        // Cadre de la veine.
        if let Some(c) = cam {
            let a = aspect as f64;
            let to = |p: [f64; 3]| {
                c.project(p, a).map(|q| {
                    Pos2::new(
                        img_rect.center().x + q[0] as f32 * img_rect.width() / 2.0,
                        img_rect.center().y - q[1] as f32 * img_rect.height() / 2.0,
                    )
                })
            };
            let [x, y, z] = [dims[0] as f64 - 1.0, dims[1] as f64 - 1.0, dims[2] as f64 - 1.0];
            let pts = [
                [0.0, 0.0, 0.0],
                [x, 0.0, 0.0],
                [x, y, 0.0],
                [0.0, y, 0.0],
                [0.0, 0.0, z],
                [x, 0.0, z],
                [x, y, z],
                [0.0, y, z],
            ];
            let edges = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)];
            let p = painter.with_clip_rect(img_rect);
            let stroke = (1.0, theme::edge());
            for (i, j) in edges {
                if let (Some(a), Some(b)) = (to(pts[i]), to(pts[j])) {
                    p.line_segment([a, b], stroke);
                }
            }
            // Flèche du vent à l'entrée.
            let mid = [0.0, y / 2.0, z / 2.0];
            if let (Some(a), Some(b)) = (to(mid), to([2.0 * l as f64, y / 2.0, z / 2.0])) {
                p.arrow(a, b - a, (2.0, theme::ink()));
                p.text(a, egui::Align2::RIGHT_BOTTOM, "V∞ ", egui::FontId::proportional(13.0), theme::ink());
            }
        }
        let st = &self.settings;
        let comp = ["ωx", "ωy", "ωz"][plane_pos.0];
        let plane = match st.field {
            SliceField::Speed => "coupe |u|/U de 0 à 1,6".to_string(),
            _ => format!("coupe {comp}·L/U de −4 à +4"),
        };
        painter.text(
            rect.left_bottom() + Vec2::new(8.0, -8.0),
            egui::Align2::LEFT_BOTTOM,
            format!(
                "{plane}{} · L = {l} mailles",
                if st.iso {
                    format!(" · isosurfaces Q = {:.2} (U/L)², teintées par la vitesse", st.iso_level)
                } else {
                    String::new()
                }
            ),
            egui::FontId::proportional(13.0),
            theme::ink(),
        );
    }

    fn history(&self, ui: &mut egui::Ui) {
        let s = self.shared.lock().unwrap();
        let h = &s.frame.hist;
        let t_max = h.last().map_or(10.0, |p| p[0]).max(10.0);
        let t_min = (t_max - 40.0).max(0.0);
        let pts = |k: usize| -> Vec<[f64; 2]> { h.iter().filter(|p| p[0] >= t_min).map(|p| [p[0], p[k]]).collect() };
        Plot::new("s3_hist")
            .height(ui.available_height().max(80.0))
            .legend(Legend::default())
            .allow_scroll(false)
            .x_axis_label("t*")
            .show(ui, |p| {
                p.line(Line::new("Cx", PlotPoints::from(pts(1))).color(theme::accent()));
                p.line(Line::new("Cz", PlotPoints::from(pts(2))).color(theme::lift()));
            });
    }
}

/// Nombre à trois chiffres significatifs.
fn sig(v: f64) -> String {
    let a = v.abs();
    if a == 0.0 {
        "0".into()
    } else if !(1e-3..1e5).contains(&a) {
        format!("{v:.2e}")
    } else {
        let d = (2 - a.log10().floor() as i32).max(0) as usize;
        format!("{v:.d$}")
    }
}

/// Mise en place de la maquette STL : dérapage, axes, unité, références,
/// échelle physique, et ce qu'elle donne (dimensions, blocage).
fn mesh_controls(ui: &mut egui::Ui, s: &mut Settings, preview: Option<&Preview>) {
    let info = preview.and_then(|p| p.as_ref().ok()).map(|(_, i)| i.clone());
    let st = &mut s.setup;
    ui.add(egui::Slider::new(&mut st.beta_deg, -30.0..=30.0).step_by(1.0).text("dérapage (°)"))
        .on_hover_text("Positif : vent venant de la droite de la maquette");
    let axis_box = |ui: &mut egui::Ui, id: &str, v: &mut Axis| {
        egui::ComboBox::from_id_salt(id)
            .width(52.0)
            .selected_text(v.label())
            .show_ui(ui, |ui| {
                for a in Axis::ALL {
                    ui.selectable_value(v, a, a.label());
                }
            });
    };
    ui.horizontal(|ui| {
        ui.label("Écoulement selon");
        axis_box(ui, "stl_flow", &mut st.flow);
        ui.label("haut selon");
        axis_box(ui, "stl_up", &mut st.up);
    })
    .response
    .on_hover_text("Axes du fichier : sens de l'écoulement (de l'avant vers l'arrière) et verticale");
    ui.horizontal(|ui| {
        ui.label("Unité du fichier");
        for u in Unit::ALL {
            ui.selectable_value(&mut st.unit, u, u.label());
        }
    });
    let u = st.unit.label();
    let m = st.unit.meters();
    ui.horizontal(|ui| {
        let mut on = st.ref_length.is_some();
        ui.checkbox(&mut on, "Longueur de réf.");
        if !on {
            st.ref_length = None;
            ui.label(egui::RichText::new("longueur hors tout").weak());
        } else if st.ref_length.is_none() {
            st.ref_length = info.as_ref().map(|i| i.ref_length_m / m);
        }
        if let Some(v) = &mut st.ref_length {
            ui.add(egui::DragValue::new(v).speed(1.0).range(1e-6..=1e9).suffix(format!(" {u}")));
        }
    });
    ui.horizontal(|ui| {
        let mut on = st.ref_area.is_some();
        ui.checkbox(&mut on, "Surface de réf.");
        if !on {
            st.ref_area = None;
            ui.label(egui::RichText::new("maître-couple").weak());
        } else if st.ref_area.is_none() {
            st.ref_area = info.as_ref().map(|i| i.frontal_m2 / (m * m));
        }
        if let Some(v) = &mut st.ref_area {
            ui.add(egui::DragValue::new(v).speed(10.0).range(1e-9..=1e15).suffix(format!(" {u}²")));
        }
    });
    match preview {
        Some(Ok((_, i))) => {
            let [a, b, c] = i.size_m.map(|v| sig(v / m));
            ui.label(
                egui::RichText::new(format!(
                    "{} triangles · {a} × {b} × {c} {u} (longueur, largeur, hauteur) · maître-couple {} {u}² · blocage {:.1} %",
                    i.triangles,
                    sig(i.frontal_m2 / (m * m)),
                    100.0 * i.blockage
                ))
                .small(),
            );
            if i.blockage > 0.05 {
                ui.colored_label(theme::warn(), "Blocage de la veine supérieur à 5 % : coefficients surestimés.");
            }
            let d = i.diag.summary();
            if !d.is_empty() {
                let more = if i.diag.watertight() {
                    ""
                } else {
                    " : maillage ouvert, appartenance votée sur trois rayons (x, y, z)"
                };
                ui.colored_label(theme::warn(), format!("{d}{more}"));
            }
        }
        Some(Err(e)) => {
            ui.colored_label(Color32::RED, e);
        }
        None => {}
    }
    ui.horizontal(|ui| {
        ui.checkbox(&mut s.phys, "Échelle physique");
        if s.phys {
            ui.add(egui::DragValue::new(&mut s.speed).speed(0.5).range(0.01..=1000.0).suffix(" m/s"));
            egui::ComboBox::from_id_salt("stl_fluid")
                .selected_text(s.fluid.label())
                .show_ui(ui, |ui| {
                    for f in [Fluid::Air, Fluid::Water] {
                        ui.selectable_value(&mut s.fluid, f, f.label());
                    }
                });
        }
    })
    .response
    .on_hover_text("Le Reynolds découle de la vitesse et de la longueur de référence ; efforts en newtons");
}

/// Repères de la littérature pour la configuration courante.
fn reference(shape: Shape3D, re: f64) -> String {
    match shape {
        Shape3D::Sphere => {
            let cd = crate::sphere_cd(re);
            let regime = if re < 24.0 {
                "écoulement attaché, sans recirculation"
            } else if re > 270.0 {
                "lâcher de tourbillons en épingle (St ≈ 0,13 à Re = 300)"
            } else if re > 210.0 {
                "sillage stationnaire oblique (portance latérale)"
            } else {
                "sillage stationnaire axisymétrique"
            };
            format!(
                "Référence sphère (Schiller–Naumann, Clift) : Cx ≈ {cd:.2} ; {regime}. Le confinement latéral (veine 6D périodique) augmente Cx de quelques %."
            )
        }
        Shape3D::Cylinder => {
            if re < 190.0 {
                "Cylindre : sillage 2D (allée de von Kármán dès Re ≈ 47) ; Cx ≈ 1,3 à Re = 100.".into()
            } else {
                "Cylindre : au-delà de Re ≈ 190, le sillage devient tridimensionnel (mode A, puis B vers 260) ; Cx ≈ 1,2 mesuré, contre ≈ 1,4 en 2D.".into()
            }
        }
        Shape3D::Cube => "Cube face au vent : décollement fixé par les arêtes ; Cx ≈ 1,05 aux grands Reynolds, plus élevé aux faibles.".into(),
        Shape3D::Wing => "Aile rectangulaire d'allongement 4 : tourbillons marginaux (plan transversal, ωx) ; à faible Reynolds, la couche limite épaisse réduit nettement la portance.".into(),
        Shape3D::Mesh => "Maquette maillée (STL) : coefficients rapportés au maître-couple.".into(),
    }
}

/// Fil de calcul : construit la veine, avance par lots d'environ 30 ms,
/// relit les efforts et publie une coupe à chaque lot.
fn worker(shared: Arc<Mutex<Shared>>, repaint: impl Fn() + Send + 'static) {
    let mut sim: Option<Gpu3D> = None;
    let mut built: Option<Settings> = None;
    let mut body: Option<Body3D> = None;
    let mut hist: Vec<[f64; 4]> = Vec::new();
    let mut t_change = 0.0f64;
    let mut steps = 8usize;
    let mut slice_id = 0u64;
    let mut perf = (Instant::now(), 0u64, 0.0f64);
    // Après une divergence, le calcul attend « Recommencer ».
    let mut halted = false;
    loop {
        let (st, active) = {
            let s = shared.lock().unwrap();
            (s.settings.clone(), s.active)
        };
        if !active {
            std::thread::sleep(Duration::from_millis(100));
            continue;
        }
        // Reconstruction si la maquette, la résolution ou la précision changent.
        let need = match &built {
            None => true,
            Some(b) => {
                b.shape != st.shape
                    || b.l != st.l
                    || b.alpha != st.alpha
                    || b.fp16 != st.fp16
                    || b.reset != st.reset
                    || (b.re >= lbm::tunnel::RE_SA_AUTO) != (st.re >= lbm::tunnel::RE_SA_AUTO)
                    || (st.shape == Shape3D::Mesh && (b.stl != st.stl || b.setup != st.setup))
            }
        };
        if need {
            sim = None;
            {
                let mut s = shared.lock().unwrap();
                s.frame.status = "Préparation".into();
                s.frame.error = None;
            }
            repaint();
            let mut mesh_info = None;
            let b = match (st.shape, &st.stl) {
                (Shape3D::Mesh, Some(f)) => {
                    let setup = MeshSetup {
                        alpha_deg: st.alpha,
                        ..st.setup.clone()
                    };
                    Body3D::from_mesh(&f.tris, st.l, &setup).map(|(b, i)| {
                        mesh_info = Some(i);
                        b
                    })
                }
                (Shape3D::Mesh, None) => Err("aucun fichier STL chargé".into()),
                (k, _) => Ok(Body3D::new(k, st.l, st.alpha)),
            };
            let nu = U as f64 * st.l as f64 / st.re;
            let sa = st.re >= lbm::tunnel::RE_SA_AUTO;
            match b.and_then(|b| Gpu3D::with_turbulence(&b, nu, U, st.fp16, sa).map(|g| (g, b))) {
                Ok((g, b)) => {
                    let lost = mesh_info.as_ref().and_then(|i: &lbm::stl::MeshInfo| i.volume_lost(g.solid_count, st.l));
                    shared.lock().unwrap().frame.warning = match lost {
                        _ if b.mesh.is_some() && g.wall_count == 0 => {
                            Some("Aucun lien de paroi : maquette invisible à cette résolution, l'augmenter.".into())
                        }
                        Some(v) if v > 0.2 => Some(format!(
                            "{:.0} % du volume tient dans des parties plus fines qu'une maille, traitées en parois minces (moins précises) : augmenter la résolution pour les résoudre.",
                            100.0 * v
                        )),
                        _ => None,
                    };
                    sim = Some(g);
                    body = Some(b);
                }
                Err(e) => {
                    let mut s = shared.lock().unwrap();
                    s.frame.error = Some(format!("Calcul 3D impossible : {e}"));
                    s.frame.status = "Indisponible".into();
                    drop(s);
                    repaint();
                }
            }
            built = Some(st.clone());
            hist.clear();
            halted = false;
            t_change = 0.0;
            steps = 8;
            continue;
        }
        let (Some(g), Some(b)) = (sim.as_mut(), body.as_ref()) else {
            std::thread::sleep(Duration::from_millis(100));
            continue;
        };
        let spt = st.l as f64 / U as f64;
        if built.as_ref().is_some_and(|x| x.re != st.re) {
            g.set_viscosity(U as f64 * st.l as f64 / st.re);
            t_change = g.time as f64 / spt;
            if let Some(x) = built.as_mut() {
                x.re = st.re;
            }
        }
        let ramp = Ramp {
            steps: (RAMP_T * spt).round() as u64,
            u: U,
        };
        let t0 = Instant::now();
        if st.running && !halted {
            g.run(steps, &ramp);
        } else {
            // Une passe sans avancer n'existe pas : en pause, on ne fait que la coupe.
            std::thread::sleep(Duration::from_millis(40));
        }
        let f = g.take_forces();
        let dt = t0.elapsed().as_secs_f64();
        if st.running && !halted {
            if dt < 0.02 {
                steps = (steps as f64 * 1.25).ceil().min(2000.0) as usize;
            } else if dt > 0.04 {
                steps = ((steps as f64 * 0.75).floor() as usize).max(1);
            }
            perf.1 += steps as u64;
        }
        let q = 0.5 * (U as f64).powi(2) * b.area;
        let t_star = g.time as f64 / spt;
        let diverged = halted || (f.steps > 0 && !g.is_finite(&f));
        halted = diverged;
        if f.steps > 0 && !diverged {
            hist.push([t_star, f.fx / q, f.fz / q, f.fy / q]);
            if hist.len() > 20000 {
                hist.drain(..5000);
            }
        }
        // Plan de coupe : position réglée, sinon centre de la maquette (sillage à 2,5 L en transversal).
        let ax = st.plane.axis();
        let n = [g.nx, g.ny, g.nz][ax];
        let cell = match st.pos {
            Some(p) => (p * (n - 1) as f64).round() as usize,
            None => {
                let c = b.center[ax];
                let c = if ax == 0 { (c + 2.5 * st.l as f64).min(n as f64 - 2.0) } else { c };
                c.round() as usize
            }
        }
        .min(n - 1);
        let (sl, image, camera) = if st.view3d {
            let cam = Camera::orbit(b.center, st.az, st.el, st.dist * st.l as f64);
            let rgba = g.render(&RenderOptions {
                camera: cam,
                width: st.img.0,
                height: st.img.1,
                plane_axis: ax,
                plane_pos: cell as f64,
                vorticity: st.field != SliceField::Speed,
                iso: st.iso,
                iso_level: st.iso_level,
                l: st.l as f64,
                background: rgbf(theme::canvas()),
                body: rgbf(theme::body()),
            });
            (None, Some(Arc::new((st.img.0, st.img.1, rgba))), Some(cam))
        } else {
            (Some(Arc::new(g.slice(ax, cell, st.field, st.l as f64))), None, None)
        };
        slice_id += 1;
        // Analyse sur la fenêtre établie.
        let settle = if t_change == 0.0 { RAMP_T + SETTLE_T } else { SETTLE_T };
        let start = (t_change + settle).max(t_star - WINDOW_T);
        let w: Vec<&[f64; 4]> = hist.iter().filter(|p| p[0] >= start).collect();
        let tt: Vec<f64> = w.iter().map(|p| p[0]).collect();
        let a = analyze(
            &tt,
            &w.iter().map(|p| p[1]).collect::<Vec<_>>(),
            &w.iter().map(|p| p[2]).collect::<Vec<_>>(),
            &w.iter().map(|p| p[3]).collect::<Vec<_>>(),
        );
        let ready = t_star - start >= 3.0 && a.n >= 8 && a.cx.is_finite();
        let status = if g.lost {
            "GPU perdu"
        } else if diverged {
            "Divergence"
        } else if !st.running {
            "En pause"
        } else if t_star < RAMP_T {
            "Mise en vitesse"
        } else if ready {
            if a.periodic && a.cz_amp.max(a.cx_amp) > 0.01 {
                "Instationnaire"
            } else {
                "Établi"
            }
        } else {
            "Transitoire"
        };
        let el = perf.0.elapsed().as_secs_f64();
        if el > 1.0 {
            perf.2 = perf.1 as f64 * g.n as f64 / el / 1e6;
            perf = (Instant::now(), 0, perf.2);
        }
        {
            let mut s = shared.lock().unwrap();
            let fr = &mut s.frame;
            fr.slice = sl;
            fr.image = image;
            fr.camera = camera;
            fr.slice_id = slice_id;
            fr.plane_pos = (ax, cell);
            fr.dims = [g.nx, g.ny, g.nz];
            fr.center = b.center;
            fr.l = st.l;
            fr.t_star = t_star;
            fr.hist = hist.clone();
            fr.analysis = a;
            fr.ready = ready;
            fr.status = status.into();
            fr.mlups = perf.2;
            fr.info = format!(
                "Grille {}×{}×{} ({:.2} M mailles) · {} mailles pariétales · τ = {:.4} · sous-maille {}{} · {} · {} · {:.0} M mailles/s",
                g.nx,
                g.ny,
                g.nz,
                g.n as f64 / 1e6,
                g.wall_count,
                g.tau0,
                if g.les { "active" } else { "inactive" },
                if g.sa { " · Spalart–Allmaras + loi de paroi" } else { "" },
                if g.fp16 { "demi-précision" } else { "simple précision" },
                g.adapter_name(),
                perf.2
            );
            if g.lost {
                fr.error = Some("Le GPU ne répond plus (pilote réinitialisé, souvent après une surchauffe) : Recommencer.".into());
            } else if diverged {
                fr.error = Some("Le calcul a divergé : baisser le Reynolds ou augmenter la résolution, puis Recommencer.".into());
            }
        }
        repaint();
        if diverged {
            std::thread::sleep(Duration::from_millis(200));
        }
    }
}
