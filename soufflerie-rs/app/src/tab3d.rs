//! Onglet « Soufflerie 3D » : maquette → réglages → calcul → résultats.
//!
//! - barre du haut : ouverture d'une maquette, préréglages, lancement, pause,
//!   progression (temps restant) ;
//! - à gauche, les réglages en trois groupes repliables (maquette,
//!   écoulement, maillage et durée) avec l'estimation du coût (mailles,
//!   mémoire du GPU, durée) avant le lancement ; les réglages physiques
//!   modifiés pendant un calcul attendent « Appliquer » ;
//! - au centre, la vue 3D (caméra orbitale, vues prédéfinies, échelle de
//!   couleur) ou la coupe à plat ;
//! - à droite, les résultats : coefficients selon chaque mesure (surface,
//!   échange aux parois, bilan sur boîtes), moment de tangage, répartition en
//!   envergure, références et exports ;
//! - en bas, la convergence des efforts et la fenêtre de moyenne.
//!
//! Le calcul tourne dans un fil dédié (`tab3d_calc`).

use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use eframe::egui::{self, Color32, Key, Modifiers, Pos2, Rect, RichText, Sense, Vec2};
use egui_plot::{Legend, Line, LineStyle, Plot, PlotPoints, Span, VLine};

use lbm::bodies3d::{Body3D, Bounds, Shape3D};
use lbm::fluid::Fluid;
use lbm::gpu3d::{bytes_per_cell, BodyView, SliceField};
use lbm::stl::{Axis, MeshInfo, MeshSetup, Unit};

use crate::tab3d_calc::{self as calc, Frame, Plane, Settings, Shared, StlFile, U};
use crate::theme;

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

/// Résolution de l'aperçu de mise en place (les dimensions de la veine
/// sont proportionnelles à L).
const PREVIEW_L: usize = 32;

/// Mises à jour de maille par seconde supposées tant qu'aucun calcul n'a
/// été mesuré dans l'interface (GTX 1070, image comprise).
const MLUPS_GUESS: [f64; 2] = [500.0, 300.0];

/// Au-delà (mailles), les réglages ne s'appliquent pas d'eux-mêmes : la
/// mise en place prend du temps et de la mémoire du GPU.
const HEAVY_CELLS: f64 = 8e6;

/// Caméra par défaut : azimut, élévation, distance.
const CAM_DEFAULT: (f64, f64, f64) = (62.0, 24.0, 6.0);
/// Vues prédéfinies : nom, touche, azimut, élévation, distance.
const VIEWS: [(&str, &str, f64, f64, f64); 5] = [
    ("3/4", "1", 62.0, 24.0, 6.0),
    ("Face", "2", 0.0, 4.0, 5.0),
    ("Côté", "3", 90.0, 0.0, 6.0),
    ("Dessus", "4", 90.0, 85.0, 6.0),
    ("Arrière", "5", 215.0, 20.0, 6.0),
];

/// Préréglage du maillage et de la durée.
struct Preset {
    name: &'static str,
    hint: &'static str,
    l: usize,
    refine: usize,
    bands: u8,
    chord: bool,
    t_end: f64,
}

const PRESETS: [Preset; 4] = [
    Preset {
        name: "Aperçu",
        hint: "L = 12, grille uniforme, 15 t* : mise en place et premier regard en quelques secondes ; coefficients grossiers.",
        l: 12,
        refine: 0,
        bands: 0,
        chord: false,
        t_end: 15.0,
    },
    Preset {
        name: "Rapide",
        hint: "L = 16 et grille ×2 autour de la maquette, 25 t* : tendances (portance à quelques % près).",
        l: 16,
        refine: 1,
        bands: 0,
        chord: false,
        t_end: 25.0,
    },
    Preset {
        name: "Précis",
        hint: "L = 16, grille ×2 et une couche pariétale ×4 (maquette STL), 30 t* : couche limite et bord d'attaque mieux résolus.",
        l: 16,
        refine: 1,
        bands: 1,
        chord: true,
        t_end: 30.0,
    },
    Preset {
        name: "Validation",
        hint: "Règle retenue pour les ailes : grille ×2 et couches pariétales automatiques d'au moins 1/8 de corde \
               (run3d --raff 1 --bande auto --bande-corde 0.125), 40 t*. Long : souvent une heure ou plus.",
        l: 16,
        refine: 1,
        bands: 4,
        chord: true,
        t_end: 40.0,
    },
];

/// Mise en place calculée pour l'aperçu (à L = PREVIEW_L).
#[derive(Clone)]
struct PreviewInfo {
    dims: [usize; 3],
    bounds: Bounds,
    /// Maquette STL : grandeurs de référence et diagnostic.
    info: Option<MeshInfo>,
    /// Plus petit rayon de bord d'attaque (mailles à PREVIEW_L), pour les
    /// couches automatiques.
    le_radius: f64,
}

type Preview = Result<PreviewInfo, String>;

/// Coût estimé d'un calcul.
struct Estimate {
    /// Mailles, niveau (×2^niveau), nature.
    grids: Vec<(f64, u32, &'static str)>,
    total: f64,
    /// Mises à jour de maille par pas de la veine.
    per_step: f64,
    mem_mb: f64,
    seconds: f64,
    bands: usize,
}

/// Demande d'export d'image en attente de la capture de la fenêtre.
struct ImageRequest {
    path: PathBuf,
    rect: Rect,
}

pub struct Tab3D {
    shared: Arc<Mutex<Shared>>,
    /// Réglages affichés (brouillon pour la partie physique).
    settings: Settings,
    /// Réglages physiques transmis au calcul.
    applied: Settings,
    /// Dernière modification des réglages physiques (application
    /// automatique quand le calcul ne tourne pas).
    changed_at: Option<Instant>,
    /// Chargements STL successifs (identifiants).
    stl_count: u64,
    stl_error: Option<String>,
    /// Mise en place de l'aperçu : (forme, fichier, réglage, incidence) → résultat.
    preview: Option<((Shape3D, u64, MeshSetup, f64), Preview)>,
    /// Réglage en attente de stabilisation (gros maillages) et depuis quand.
    pending: Option<((Shape3D, u64, MeshSetup, f64), Instant)>,
    texture: Option<egui::TextureHandle>,
    slice_id: u64,
    started: bool,
    /// Caméra visée (les vues prédéfinies y glissent en douceur).
    cam_target: (f64, f64, f64),
    /// Zone de la vue (export d'image).
    view_rect: Rect,
    /// Débit mesuré (M mises à jour par seconde) du dernier calcul.
    measured_mlups: Option<(bool, f64)>,
    /// Autoriser un calcul au-delà des garde-fous.
    force: bool,
    /// Réglages affichés au-delà des garde-fous.
    guarded: bool,
    /// Calcul lourd (mise en place longue) : pas d'application automatique.
    heavy: bool,
    /// Point de référence du moment (fraction de L depuis l'avant).
    xref: f64,
    /// Valeurs de référence saisies (mesure, autre code) : Cz, Cx.
    refs: (Option<f64>, Option<f64>),
    help: bool,
    message: Arc<Mutex<Option<(String, Instant)>>>,
    image_req: Arc<Mutex<Option<ImageRequest>>>,
    last_frame: Instant,
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
            reset: 0,
            refine: 0,
            bands: 0,
            band_chord: true,
            wall: 0,
            xtr: 0.0,
            // Le calcul attend « Lancer » : le GPU est souvent partagé.
            running: false,
            t_end: 40.0,
            field: SliceField::Vorticity,
            plane: Plane::Y,
            pos: None,
            view3d: true,
            az: CAM_DEFAULT.0,
            el: CAM_DEFAULT.1,
            dist: CAM_DEFAULT.2,
            iso: true,
            body_view: BodyView::Smooth,
            iso_level: 0.4,
            img: (800, 450),
        };
        Tab3D {
            shared: Arc::new(Mutex::new(Shared {
                settings: settings.clone(),
                frame: Frame::default(),
                active: false,
            })),
            applied: settings.clone(),
            settings,
            changed_at: None,
            stl_count: 0,
            stl_error: None,
            preview: None,
            pending: None,
            texture: None,
            slice_id: 0,
            started: false,
            cam_target: CAM_DEFAULT,
            view_rect: Rect::NOTHING,
            measured_mlups: None,
            force: false,
            guarded: false,
            heavy: false,
            xref: 0.25,
            refs: (None, None),
            help: false,
            message: Arc::new(Mutex::new(None)),
            image_req: Arc::new(Mutex::new(None)),
            last_frame: Instant::now(),
        }
    }

    /// Réglages de la ligne de commande : préréglage (`--prereglage
    /// apercu|rapide|precis|validation`), lancement immédiat (`--lancer 1`),
    /// écoulement (`--re`, `--alpha`) et mise en place de la maquette STL,
    /// avec les mêmes options que `run3d` (`--unite`, `--ecoulement`,
    /// `--haut`, `--beta`, `--lref`, `--sref`, `--marge`, `--lisser`).
    pub fn start_options(&mut self, opts: &std::collections::HashMap<String, String>) {
        let num = |k: &str| opts.get(k).and_then(|v| v.parse::<f64>().ok());
        let norm = |s: &str| s.to_lowercase().replace(['é', 'è'], "e").replace('ç', "c");
        if let Some(p) = opts.get("prereglage").and_then(|n| PRESETS.iter().find(|p| norm(p.name) == norm(n))) {
            self.apply_preset(p);
        }
        let s = &mut self.settings;
        if let Some(v) = num("re") {
            s.re = v;
        }
        if let Some(v) = num("alpha") {
            s.alpha = v;
        }
        let st = &mut s.setup;
        if let Some(u) = opts.get("unite").and_then(|v| Unit::parse(v)) {
            st.unit = u;
        }
        if let Some(a) = opts.get("ecoulement").and_then(|v| Axis::parse(v)) {
            st.flow = a;
        }
        if let Some(a) = opts.get("haut").and_then(|v| Axis::parse(v)) {
            st.up = a;
        }
        if let Some(v) = num("beta") {
            st.beta_deg = v;
        }
        if let Some(v) = num("lref") {
            st.ref_length = Some(v);
        }
        if let Some(v) = num("sref") {
            st.ref_area = Some(v);
        }
        if let Some(v) = num("marge") {
            st.margin = v;
        }
        if let Some(v) = num("lisser") {
            st.smooth = v > 0.0;
        }
        s.running = opts.get("lancer").is_some_and(|v| v != "0");
        self.apply();
    }

    /// Onglet affiché ou non (le calcul s'arrête quand il est masqué).
    pub fn set_active(&mut self, on: bool, ctx: &egui::Context) {
        self.shared.lock().unwrap().active = on;
        if on && !self.started {
            self.started = true;
            let shared = self.shared.clone();
            let ctx = ctx.clone();
            std::thread::spawn(move || calc::worker(shared, move || ctx.request_repaint()));
        }
    }

    /// Transmet au calcul les réglages appliqués et l'affichage courant.
    fn push(&self) {
        let mut s = self.applied.clone();
        s.take_view(&self.settings);
        self.shared.lock().unwrap().settings = s;
    }

    /// Applique les réglages physiques affichés.
    fn apply(&mut self) {
        self.applied = self.settings.clone();
        self.changed_at = None;
        self.push();
    }

    /// Rien n'a encore été calculé (t* = 0) : appliquer ne perd rien.
    fn untouched(&self) -> bool {
        let s = self.shared.lock().unwrap();
        s.frame.t_star == 0.0 || s.frame.error.is_some()
    }

    /// Réglages physiques modifiés et pas encore appliqués ?
    fn dirty(&self) -> bool {
        !self.settings.same_physics(&self.applied)
    }

    fn flash(&self, msg: String) {
        *self.message.lock().unwrap() = Some((msg, Instant::now()));
    }

    /// Lit un fichier STL, ou de géométrie (`.sgeo` : mise en place et
    /// arêtes vives comprises), et en fait la maquette courante.
    pub fn load_stl(&mut self, path: &Path) {
        let read = if lbm::sgeo::is_sgeo(path) {
            lbm::sgeo::GeoFile::read(path).map(|g| {
                let st = &mut self.settings.setup;
                st.unit = g.unit.unwrap_or(st.unit);
                st.flow = g.flow.unwrap_or(st.flow);
                st.up = g.up.unwrap_or(st.up);
                st.ref_length = g.ref_length;
                st.ref_area = g.ref_area;
                st.beta_deg = g.beta_deg.unwrap_or(0.0);
                st.ground = g.ground;
                st.edges = Some(Arc::new(g.feature_segments().iter().map(|e| [e.a, e.b]).collect()));
                if let Some(a) = g.alpha_deg {
                    self.settings.alpha = a;
                }
                g.tris()
            })
        } else {
            self.settings.setup.edges = None;
            lbm::stl::read_stl(path)
        };
        match read {
            Ok(tris) => {
                self.stl_count += 1;
                let wetted = tris
                    .iter()
                    .map(|t| {
                        let u = [0, 1, 2].map(|k| t[1][k] - t[0][k]);
                        let v = [0, 1, 2].map(|k| t[2][k] - t[0][k]);
                        let c = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
                        0.5 * (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt()
                    })
                    .sum();
                self.settings.stl = Some(StlFile {
                    name: path
                        .file_name()
                        .map_or_else(|| path.display().to_string(), |n| n.to_string_lossy().into()),
                    path: path.to_path_buf(),
                    tris: Arc::new(tris),
                    id: self.stl_count,
                    wetted,
                });
                self.settings.shape = Shape3D::Mesh;
                self.settings.pos = None;
                // Nouvelle maquette : on la met en place, le calcul attend.
                self.settings.running = false;
                self.stl_error = None;
                self.apply();
            }
            Err(e) => self.stl_error = Some(e),
        }
    }

    fn open_dialog(&mut self) {
        if let Some(p) = rfd::FileDialog::new()
            .add_filter("Maillage STL ou géométrie", &["stl", "STL", "sgeo"])
            .pick_file()
        {
            self.load_stl(&p);
        }
    }

    /// Mise en place de la maquette courante (aperçu à L = 32), recalculée
    /// quand le réglage change ; pour un gros maillage, une fois le réglage
    /// stable depuis 0,25 s (le curseur reste fluide, l'aperçu précédent
    /// reste affiché entre-temps).
    fn mesh_preview(&mut self, ctx: &egui::Context) -> Option<Preview> {
        let s = &self.settings;
        let mesh = s.shape == Shape3D::Mesh;
        let f = s.stl.as_ref();
        if mesh && f.is_none() {
            return None;
        }
        let key = (s.shape, f.map_or(0, |f| f.id), if mesh { s.setup.clone() } else { MeshSetup::default() }, s.alpha);
        let stale = self.preview.as_ref().is_none_or(|(k, _)| *k != key);
        let mut due = self.preview.is_none() || f.is_none_or(|f| f.tris.len() < 20_000);
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
            let r = match f.filter(|_| mesh) {
                Some(f) => {
                    let setup = MeshSetup {
                        alpha_deg: s.alpha,
                        ..s.setup.clone()
                    };
                    Body3D::from_mesh(&f.tris, PREVIEW_L, &setup).map(|(b, i)| PreviewInfo {
                        dims: [b.nx, b.ny, b.nz],
                        bounds: b.bounds,
                        le_radius: lbm::gpu3d::auto_band_count(&b, 0).1,
                        info: Some(i),
                    })
                }
                None => {
                    let b = Body3D::new(s.shape, PREVIEW_L, s.alpha);
                    Ok(PreviewInfo {
                        dims: [b.nx, b.ny, b.nz],
                        bounds: b.bounds,
                        info: None,
                        le_radius: f64::INFINITY,
                    })
                }
            };
            self.preview = Some((key, r));
        }
        self.preview.as_ref().map(|(_, r)| r.clone())
    }

    /// Débit supposé (M mises à jour par seconde) : mesuré au dernier calcul
    /// de même précision, sinon valeur type.
    fn mlups(&self, fp16: bool) -> f64 {
        match self.measured_mlups {
            Some((p, v)) if p == fp16 && v > 0.0 => v,
            _ => MLUPS_GUESS[if fp16 { 0 } else { 1 }],
        }
    }

    /// Coût estimé d'un calcul : mêmes grilles boîtes que
    /// `Gpu3D::with_refinement_edges`, couches pariétales estimées par la
    /// surface mouillée (stockage creux).
    fn estimate(&self, s: &Settings, p: &PreviewInfo) -> Estimate {
        let k = s.l as f64 / PREVIEW_L as f64;
        let mut dims = p.dims.map(|v| v * s.l / PREVIEW_L);
        let b = p.bounds;
        let mut bb = [b.x0 * k, b.y0 * k, b.z0 * k, b.x1 * k, b.y1 * k, b.z1 * k];
        let mut grids = vec![((dims[0] * dims[1] * dims[2]) as f64, 0u32, "veine")];
        let l0 = s.l as f64;
        let depth = s.refine.min(lbm::refine3d::EXTENT3.len());
        let mut ok = true;
        for lev in 1..=depth {
            let (up, down, side) = lbm::refine3d::EXTENT3[lev - 1];
            let sp = (1u64 << (lev - 1)) as f64;
            let lo = [
                (bb[0] - up * l0 * sp).floor().max(3.0),
                (bb[1] - side * l0 * sp).floor().max(3.0),
                (bb[2] - side * l0 * sp).floor().max(3.0),
            ];
            let hi = [
                (bb[3] + down * l0 * sp).ceil().min(dims[0] as f64 - 4.0),
                (bb[4] + side * l0 * sp).ceil().min(dims[1] as f64 - 4.0),
                (bb[5] + side * l0 * sp).ceil().min(dims[2] as f64 - 4.0),
            ];
            if (0..3).any(|a| hi[a] <= lo[a] + 6.0) {
                ok = false;
                break;
            }
            dims = [0, 1, 2].map(|a| (2.0 * (hi[a] - lo[a]) + 1.0) as usize);
            bb = [
                2.0 * (bb[0] - lo[0]),
                2.0 * (bb[1] - lo[1]),
                2.0 * (bb[2] - lo[2]),
                2.0 * (bb[3] - lo[0]),
                2.0 * (bb[4] - lo[1]),
                2.0 * (bb[5] - lo[2]),
            ];
            grids.push(((dims[0] * dims[1] * dims[2]) as f64, lev as u32, "boîte"));
        }
        // Couches pariétales : surface mouillée × (épaisseur + raccords et
        // fantômes ; recalé sur l'aile AR 4 du banc, ±25 %).
        let mut nb = 0;
        if let (true, Some(f), Some(info), true) = (s.shape == Shape3D::Mesh, &s.stl, &p.info, ok && depth >= 1) {
            nb = match s.bands {
                0 => 0,
                4 => {
                    let r = p.le_radius * k;
                    let mut n = 0;
                    while r.is_finite() && n < 4 && r * ((1u64 << (depth + n)) as f64) < 2.0 {
                        n += 1;
                    }
                    n
                }
                n => n as usize,
            };
            let lref = info.ref_length_m / s.setup.unit.meters();
            let wet = f.wetted / (lref * lref).max(1e-30);
            for (i, th) in calc::band_thickness(nb, s.l, depth, s.band_chord).iter().enumerate() {
                let lev = depth + 1 + i;
                let lk = l0 * (1u64 << lev) as f64;
                grids.push((wet * lk * lk * (th + 16.0), lev as u32, "couche"));
            }
        }
        let total: f64 = grids.iter().map(|g| g.0).sum();
        let per_step: f64 = grids.iter().map(|g| g.0 * (1u64 << g.1) as f64).sum();
        let steps = s.t_end.max(0.0) * l0 / U as f64;
        Estimate {
            mem_mb: total * bytes_per_cell(s.fp16) as f64 / 1e6,
            seconds: if s.t_end > 0.0 { steps * per_step / (self.mlups(s.fp16) * 1e6) } else { f64::NAN },
            grids,
            total,
            per_step,
            bands: nb,
        }
    }

    /// Le calcul dépasse-t-il les garde-fous ? (message)
    fn guard(e: &Estimate, s: &Settings) -> Option<String> {
        let limit = if s.fp16 { 40e6 } else { 20e6 };
        if e.total > limit {
            Some(format!(
                "≈ {:.0} M mailles : au-delà de la limite du calcul ({:.0} M en {} précision). Réduire la résolution, le raffinement ou les couches.",
                e.total / 1e6,
                limit / 1e6,
                if s.fp16 { "demi" } else { "simple" }
            ))
        } else if e.mem_mb > 6500.0 {
            Some(format!("≈ {:.1} Go de mémoire GPU : proche de la carte (8 Go).", e.mem_mb / 1e3))
        } else {
            None
        }
    }

    fn apply_preset(&mut self, p: &Preset) {
        let s = &mut self.settings;
        s.l = p.l;
        s.refine = p.refine;
        s.bands = p.bands;
        s.band_chord = p.chord;
        s.t_end = p.t_end;
    }

    fn preset_matches(&self, p: &Preset) -> bool {
        let s = &self.settings;
        s.l == p.l && s.refine == p.refine && s.t_end == p.t_end && (s.bands == p.bands || s.shape != Shape3D::Mesh)
            && (s.bands == 0 || s.band_chord == p.chord)
    }

    /// Bouton principal : lancer, reprendre, mettre en pause ou relancer.
    fn toggle_run(&mut self) {
        if self.dirty() && self.guarded && !self.force {
            return;
        }
        let (finished, t_star) = {
            let s = self.shared.lock().unwrap();
            (s.frame.finished, s.frame.t_star)
        };
        if self.dirty() {
            self.apply();
            self.settings.running = true;
        } else if finished {
            // Terminé : on prolonge de la moitié de la durée.
            self.settings.t_end = (t_star + (0.5 * self.settings.t_end).max(10.0)).round();
            self.settings.running = true;
        } else {
            self.settings.running = !self.settings.running;
        }
        self.push();
    }

    fn restart(&mut self) {
        if self.dirty() && self.guarded && !self.force {
            return;
        }
        if self.dirty() {
            self.applied = self.settings.clone();
        }
        self.settings.reset += 1;
        self.applied.reset = self.settings.reset;
        self.settings.running = true;
        self.changed_at = None;
        self.push();
    }

    fn set_view(&mut self, v: (f64, f64, f64)) {
        self.cam_target = v;
        self.settings.view3d = true;
    }

    /// Raccourcis clavier (hors saisie de texte).
    fn shortcuts(&mut self, ctx: &egui::Context) {
        if ctx.egui_wants_keyboard_input() {
            return;
        }
        let k = |m: Modifiers, key: Key| ctx.input_mut(|i| i.consume_key(m, key));
        if k(Modifiers::NONE, Key::Space) {
            self.toggle_run();
        }
        if k(Modifiers::COMMAND, Key::O) {
            self.open_dialog();
        }
        if k(Modifiers::COMMAND, Key::R) {
            self.restart();
        }
        if k(Modifiers::COMMAND, Key::Enter) && self.dirty() {
            self.apply();
        }
        if k(Modifiers::COMMAND, Key::E) {
            self.export_csv();
        }
        if k(Modifiers::COMMAND, Key::S) {
            self.export_image(ctx);
        }
        if k(Modifiers::NONE, Key::V) {
            self.settings.view3d = !self.settings.view3d;
        }
        if k(Modifiers::NONE, Key::Q) {
            self.settings.iso = !self.settings.iso;
        }
        if k(Modifiers::NONE, Key::F1) {
            self.help = !self.help;
        }
        for (i, key) in [Key::Num1, Key::Num2, Key::Num3, Key::Num4, Key::Num5].into_iter().enumerate() {
            if k(Modifiers::NONE, key) {
                let v = VIEWS[i];
                self.set_view((v.2, v.3, v.4));
            }
        }
        if k(Modifiers::NONE, Key::Num0) || k(Modifiers::NONE, Key::Home) {
            self.set_view(CAM_DEFAULT);
        }
    }

    pub fn ui(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();
        let dropped: Vec<PathBuf> = ctx.input(|i| i.raw.dropped_files.iter().map(|f| f.path().to_path_buf()).collect());
        if let Some(p) = dropped
            .iter()
            .find(|p| p.extension().is_some_and(|e| e.eq_ignore_ascii_case("stl") || e.eq_ignore_ascii_case("sgeo")))
        {
            self.load_stl(p);
        }
        let before = self.settings.clone();
        self.shortcuts(&ctx);
        let preview = self.mesh_preview(&ctx);
        let est = preview.as_ref().and_then(|p| p.as_ref().ok()).map(|p| self.estimate(&self.settings, p));
        self.guarded = est.as_ref().is_some_and(|e| Tab3D::guard(e, &self.settings).is_some());
        self.heavy = est.as_ref().is_some_and(|e| e.total > HEAVY_CELLS);
        if !self.guarded {
            self.force = false;
        }
        egui::Panel::top("s3_barre").show(ui, |ui| self.toolbar(ui, preview.as_ref()));
        egui::Panel::left("s3_reglages")
            .default_size(330.0)
            .resizable(true)
            .show(ui, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| self.controls(ui, preview.as_ref()));
            });
        egui::Panel::right("s3_resultats")
            .default_size(340.0)
            .resizable(true)
            .show(ui, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| self.results(ui, preview.as_ref()));
            });
        egui::Panel::bottom("s3_historique")
            .default_size(170.0)
            .resizable(true)
            .show(ui, |ui| self.history(ui));
        egui::CentralPanel::default().show(ui, |ui| self.view(ui));
        self.help_window(&ctx);
        // Réglages : l'affichage passe tout de suite ; la physique s'applique
        // d'elle-même quand le calcul ne tourne pas (après 0,4 s de calme,
        // bouton relâché, sous les garde-fous), sinon attend « Appliquer ».
        if self.settings.same_physics(&before) {
            if self.settings != before {
                self.push();
            }
        } else {
            self.changed_at = Some(Instant::now());
        }
        if let Some(t) = self.changed_at {
            let idle = !ctx.input(|i| i.pointer.any_down());
            let wait = Duration::from_millis(400);
            if !self.dirty() {
                self.changed_at = None;
            } else if !self.settings.running && !self.guarded && !self.heavy && self.untouched() {
                if idle && t.elapsed() >= wait {
                    self.apply();
                } else {
                    ctx.request_repaint_after(wait);
                }
            }
        }
        // Débit mesuré, pour les estimations suivantes.
        {
            let s = self.shared.lock().unwrap();
            if s.frame.advancing && s.frame.mlups > 0.0 {
                self.measured_mlups = Some((s.frame.fp16, s.frame.mlups));
            }
        }
        if let Some((_, t)) = &*self.message.lock().unwrap() {
            if t.elapsed() < Duration::from_secs(6) {
                ctx.request_repaint_after(Duration::from_millis(500));
            }
        }
    }

    // ------------------------------------------------------------------
    // Barre du haut
    // ------------------------------------------------------------------

    fn toolbar(&mut self, ui: &mut egui::Ui, preview: Option<&Preview>) {
        ui.add_space(3.0);
        ui.horizontal(|ui| {
            if ui
                .button("📂 Ouvrir…")
                .on_hover_text("Maquette STL ou géométrie .sgeo (Ctrl+O) ; on peut aussi glisser le fichier sur la fenêtre")
                .clicked()
            {
                self.open_dialog();
            }
            let model = match (&self.settings.stl, self.settings.shape) {
                (Some(f), Shape3D::Mesh) => f.name.clone(),
                (_, k) => k.label().to_string(),
            };
            ui.label(RichText::new(model).strong()).on_hover_text(match &self.settings.stl {
                Some(f) if self.settings.shape == Shape3D::Mesh => f.path.display().to_string(),
                _ => self.settings.shape.reference().to_string(),
            });
            ui.separator();
            ui.label(RichText::new("Préréglage").weak());
            let est = preview.and_then(|p| p.as_ref().ok());
            for p in &PRESETS {
                let on = self.preset_matches(p);
                let mut tip = p.hint.to_string();
                if let Some(pi) = est {
                    let mut s = self.settings.clone();
                    s.l = p.l;
                    s.refine = p.refine;
                    s.bands = p.bands;
                    s.band_chord = p.chord;
                    s.t_end = p.t_end;
                    let e = self.estimate(&s, pi);
                    tip += &format!(
                        "\n\nIci : ≈ {} mailles, {} de mémoire GPU, ≈ {}",
                        cells(e.total),
                        mem(e.mem_mb),
                        crate::gui::duration(e.seconds)
                    );
                }
                if ui.selectable_label(on, p.name).on_hover_text(tip).clicked() {
                    self.apply_preset(p);
                }
            }
            ui.separator();
            let (status, advancing, finished, building, t_star, t_rate, error) = {
                let s = self.shared.lock().unwrap();
                let f = &s.frame;
                (
                    f.status.clone(),
                    f.advancing,
                    f.finished,
                    f.building,
                    f.t_star,
                    f.t_rate,
                    f.error.is_some(),
                )
            };
            let dirty = self.dirty();
            let (label, tip) = if dirty {
                ("▶ Appliquer et lancer", "Applique les réglages modifiés (le calcul repart de zéro si la maquette ou le maillage changent) et lance (Espace)")
            } else if self.settings.running && !finished {
                ("⏸ Pause", "Met le calcul en pause ; la vue reste manipulable (Espace)")
            } else if finished {
                ("▶ Prolonger", "Poursuit le calcul de la moitié de sa durée (Espace)")
            } else if t_star > 0.0 {
                ("▶ Reprendre", "Reprend le calcul (Espace)")
            } else {
                ("▶ Lancer", "Lance le calcul (Espace)")
            };
            let primary = !self.settings.running || dirty || finished;
            let mut btn = egui::Button::new(RichText::new(label).strong());
            if primary {
                btn = btn.fill(theme::alpha(theme::ok(), 70));
            }
            let blocked = dirty && self.guarded && !self.force;
            if ui
                .add_enabled(!building && !blocked, btn)
                .on_hover_text(tip)
                .on_disabled_hover_text("Au-delà des garde-fous (voir « Coût estimé ») : cocher « Lancer quand même »")
                .clicked()
            {
                self.toggle_run();
            }
            if ui
                .button("⟲ Recommencer")
                .on_hover_text("Repart de t* = 0 avec les réglages affichés (Ctrl+R)")
                .clicked()
            {
                self.restart();
            }
            ui.separator();
            // Progression.
            let color = if error {
                theme::err()
            } else if advancing {
                theme::ok()
            } else if finished {
                theme::lift()
            } else {
                theme::warn()
            };
            let (r, _) = ui.allocate_exact_size(Vec2::splat(10.0), Sense::hover());
            ui.painter().circle_filled(r.center(), 5.0, color);
            ui.label(RichText::new(&status).strong());
            let t_end = self.settings.t_end;
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui.button("?").on_hover_text("Raccourcis et aide (F1)").clicked() {
                    self.help = !self.help;
                }
                let w = ui.available_width().max(120.0);
                if building {
                    ui.add(egui::Spinner::new());
                    ui.label("mise en place de la veine…");
                } else {
                    let frac = if t_end > 0.0 { (t_star / t_end).clamp(0.0, 1.0) as f32 } else { 0.0 };
                    let left = if t_end > 0.0 && t_rate > 0.0 && advancing {
                        format!(" · reste ≈ {}", crate::gui::duration((t_end - t_star).max(0.0) / t_rate))
                    } else {
                        String::new()
                    };
                    let text = if t_end > 0.0 {
                        format!("t* {:.1} / {:.0}{left}", t_star, t_end)
                    } else {
                        format!("t* {:.1} (sans fin)", t_star)
                    };
                    ui.add(egui::ProgressBar::new(frac).desired_width(w).text(text).animate(advancing && t_end <= 0.0))
                        .on_hover_text(format!(
                            "Temps convectifs t* = t·U/L ; {:.2} t* par seconde",
                            t_rate
                        ));
                }
            });
        });
        ui.add_space(2.0);
    }

    // ------------------------------------------------------------------
    // Réglages
    // ------------------------------------------------------------------

    fn controls(&mut self, ui: &mut egui::Ui, preview: Option<&Preview>) {
        if self.dirty() {
            egui::Frame::new()
                .fill(theme::alpha(theme::warn(), 40))
                .corner_radius(6)
                .inner_margin(8)
                .show(ui, |ui| {
                    ui.horizontal_wrapped(|ui| {
                        ui.label(RichText::new("Réglages modifiés").strong());
                        if self.settings.running || !self.untouched() {
                            ui.label(RichText::new("les appliquer relance le calcul (résultats actuels perdus si la maquette ou le maillage changent)").small());
                        } else if self.heavy || self.guarded {
                            ui.label(RichText::new("calcul lourd : appliqués au lancement").small());
                        } else {
                            ui.label(RichText::new("appliqués dès que le réglage est stable").small());
                        }
                    });
                    ui.horizontal(|ui| {
                        if ui.button("Appliquer").on_hover_text("Ctrl+Entrée").clicked() {
                            self.apply();
                        }
                        if ui.button("Annuler").clicked() {
                            let mut s = self.applied.clone();
                            s.take_view(&self.settings);
                            self.settings = s;
                            self.changed_at = None;
                        }
                    });
                });
            ui.add_space(4.0);
        }
        let mesh = self.settings.shape == Shape3D::Mesh && self.settings.stl.is_some();
        section(ui, "1 · Maquette", true, |ui| self.model_section(ui, preview, mesh));
        section(ui, "2 · Écoulement", true, |ui| self.flow_section(ui, preview, mesh));
        section(ui, "3 · Maillage et durée", true, |ui| self.mesh_section(ui, preview, mesh));
    }

    fn model_section(&mut self, ui: &mut egui::Ui, preview: Option<&Preview>, mesh: bool) {
        let s = &mut self.settings;
        ui.horizontal_wrapped(|ui| {
            for k in Shape3D::ALL {
                if ui.selectable_value(&mut s.shape, k, k.label()).on_hover_text(k.reference()).changed() {
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
        });
        if s.stl.is_none() {
            ui.label(
                RichText::new("Maquette maillée : « Ouvrir… », Ctrl+O, ou glisser un fichier .stl / .sgeo sur la fenêtre.")
                    .small()
                    .weak(),
            );
        }
        if let Some(e) = &self.stl_error {
            banner(ui, theme::err(), "⚠", e);
        }
        let s = &mut self.settings;
        if s.shape.has_alpha() {
            ui.add(egui::Slider::new(&mut s.alpha, -10.0..=20.0).step_by(0.5).text("incidence α (°)"))
                .on_hover_text("Rotation à cabrer autour du centre de la maquette");
        }
        if mesh {
            let info = preview.and_then(|p| p.as_ref().ok()).and_then(|p| p.info.clone());
            mesh_info(ui, s, info.as_ref(), preview);
            egui::CollapsingHeader::new("Mise en place du fichier")
                .default_open(false)
                .show(ui, |ui| mesh_controls(ui, s, info.as_ref()));
        } else {
            ui.label(RichText::new(format!("Longueur de référence : {}", s.shape.reference())).small().weak());
        }
    }

    fn flow_section(&mut self, ui: &mut egui::Ui, preview: Option<&Preview>, mesh: bool) {
        let s = &mut self.settings;
        if mesh {
            ui.horizontal(|ui| {
                ui.checkbox(&mut s.phys, "Échelle physique")
                    .on_hover_text("Le Reynolds découle de la vitesse et de la longueur de référence ; efforts en newtons");
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
            });
        }
        if mesh && s.phys {
            if let Some(Ok(PreviewInfo { info: Some(info), .. })) = preview {
                let (nu, _) = s.fluid.props();
                s.re = (s.speed * info.ref_length_m / nu).max(1.0);
            }
            ui.label(format!("Reynolds {} (vitesse × longueur de référence / ν)", crate::gui::fmt_re(s.re)));
        } else {
            ui.horizontal(|ui| {
                ui.label("Reynolds");
                let mut k = RE3.iter().position(|&r| r == s.re).unwrap_or(5) as f64;
                ui.add(egui::Slider::new(&mut k, 0.0..=(RE3.len() - 1) as f64).step_by(1.0).show_value(false))
                    .on_hover_text("Re = U·L/ν, L : longueur de référence");
                s.re = RE3[k.round() as usize];
                ui.label(RichText::new(crate::gui::fmt_re(s.re)).strong());
            });
        }
        ui.horizontal(|ui| {
            ui.label("Paroi");
            ui.selectable_value(&mut s.wall, 0, "auto")
                .on_hover_text("Loi de paroi dès Re = 5·10⁴, paroi sans glissement en dessous");
            ui.selectable_value(&mut s.wall, 1, "loi de paroi").on_hover_text(
                "Spalart–Allmaras et loi de Spalding : couche limite turbulente modélisée (régime d'après la crise de traînée pour les corps épais)",
            );
            ui.selectable_value(&mut s.wall, 2, "sans").on_hover_text(
                "Paroi sans glissement, sous-maille seule : couche limite épaisse, décollement précoce (proche du régime sous-critique des corps épais)",
            );
        });
        let wm = calc::wall_model(s);
        ui.label(
            RichText::new(if wm {
                "Couche limite turbulente modélisée (Spalart–Allmaras + loi de paroi) : la grille n'a pas à la résoudre."
            } else if s.re > 5000.0 {
                "Sans loi de paroi au-delà de quelques milliers : couche limite non résolue (sous-maille seule), coefficients indicatifs."
            } else {
                "Couche limite résolue par la grille (sous-maille si nécessaire)."
            })
            .small()
            .color(if !wm && s.re > 5000.0 { theme::warn() } else { ui.visuals().weak_text_color() }),
        );
        if wm {
            ui.add(egui::Slider::new(&mut s.xtr, 0.0..=1.0).step_by(0.05).text("transition (x/c)"))
                .on_hover_text(
                    "Couche limite laminaire (frottement de Blasius, sans production turbulente) du bord d'attaque jusqu'à cette fraction de la corde locale ; 0 : turbulente dès l'avant",
                );
        }
    }

    fn mesh_section(&mut self, ui: &mut egui::Ui, preview: Option<&Preview>, mesh: bool) {
        let p_ok = preview.and_then(|p| p.as_ref().ok()).cloned();
        let s = &mut self.settings;
        ui.label("Résolution : mailles par longueur de référence");
        ui.horizontal_wrapped(|ui| {
            for (l, name) in RES3 {
                ui.selectable_value(&mut s.l, l, format!("{l}")).on_hover_text(name);
            }
        });
        ui.label("Près de la maquette").on_hover_text(
            "Grilles boîtes imbriquées, chacune deux fois plus fine : couche limite et parois minces mieux résolues pour une fraction du coût d'une grille uniforme",
        );
        ui.horizontal_wrapped(|ui| {
            ui.selectable_value(&mut s.refine, 0, "uniforme");
            for d in 1..=3usize {
                ui.selectable_value(&mut s.refine, d, format!("×{} ({})", 1 << d, s.l << d))
                    .on_hover_text(match d {
                        1 => "Une grille deux fois plus fine autour de la maquette",
                        2 => "Deux grilles imbriquées, 2 et 4 fois plus fines",
                        _ => "Trois grilles imbriquées : calcul long, mémoire GPU à surveiller",
                    });
            }
        });
        // Couches pariétales : grilles deux fois plus fines chacune, limitées
        // à une couche autour de la paroi (stockage creux), sur une maquette STL.
        let can_band = mesh && s.refine >= 1;
        ui.add_enabled_ui(can_band, |ui| {
            ui.label("Couches pariétales").on_hover_text(
                "Couches minces autour de la paroi, chacune deux fois plus fine : la résolution d'une grille entière près de la surface pour une fraction des mailles",
            );
            ui.horizontal_wrapped(|ui| {
                ui.selectable_value(&mut s.bands, 0, "aucune");
                for k in 1..=3u8 {
                    let f = 1usize << (s.refine + k as usize);
                    ui.selectable_value(&mut s.bands, k, format!("{k} (×{f})"));
                }
                ui.selectable_value(&mut s.bands, 4, "auto").on_hover_text(
                    "Assez de couches pour mettre 2 mailles dans le plus petit rayon de bord d'attaque (au plus 4), comme --bande auto",
                );
            });
            if s.bands > 0 {
                ui.checkbox(&mut s.band_chord, "Épaisseur d'au moins 1/8 de la longueur de réf.").on_hover_text(
                    "Comme --bande-corde 0.125 (règle retenue) : une couche de 16 mailles très fines ne couvre que quelques % de corde, \
                     moins que la région où couche limite et bord de fuite fixent la portance",
                );
            }
        });
        if !can_band {
            ui.label(
                RichText::new(if mesh {
                    "Couches pariétales : choisir d'abord une grille ×2 au moins."
                } else {
                    "Couches pariétales : maquettes STL seulement."
                })
                .small()
                .weak(),
            );
        }
        ui.horizontal(|ui| {
            ui.label("Précision");
            ui.selectable_value(&mut s.fp16, true, "demi").on_hover_text("2 octets par population : deux fois moins de mémoire, écarts de l'ordre de 10⁻³");
            ui.selectable_value(&mut s.fp16, false, "simple").on_hover_text("4 octets par population");
        });
        ui.horizontal(|ui| {
            ui.label("Durée");
            let mut endless = s.t_end <= 0.0;
            if !endless {
                ui.add(egui::DragValue::new(&mut s.t_end).speed(1.0).range(5.0..=1000.0).suffix(" t*"))
                    .on_hover_text(format!(
                        "Temps convectifs (t·U/L) ; mise en vitesse {:.0} t*, moyenne sur les {:.0} derniers une fois établi",
                        calc::RAMP_T,
                        calc::WINDOW_T
                    ));
            }
            if ui.checkbox(&mut endless, "sans fin").changed() {
                s.t_end = if endless { 0.0 } else { 40.0 };
            }
        });
        // Coût estimé.
        if let Some(p) = p_ok {
            let s = self.settings.clone();
            let e = self.estimate(&s, &p);
            ui.add_space(4.0);
            egui::Frame::group(ui.style()).show(ui, |ui| {
                ui.set_width(ui.available_width());
                ui.label(RichText::new("Coût estimé").strong());
                egui::Grid::new("s3_cout").num_columns(2).spacing([12.0, 2.0]).show(ui, |ui| {
                    for (n, lev, kind) in &e.grids {
                        ui.label(RichText::new(format!("{kind} ×{}", 1u64 << lev)).small());
                        ui.label(RichText::new(cells(*n)).small());
                        ui.end_row();
                    }
                    ui.label("Mailles");
                    ui.label(RichText::new(cells(e.total)).strong());
                    ui.end_row();
                    ui.label("Mémoire GPU");
                    ui.label(RichText::new(mem(e.mem_mb)).strong());
                    ui.end_row();
                    ui.label("Durée");
                    ui.label(RichText::new(format!("≈ {}", crate::gui::duration(e.seconds))).strong())
                        .on_hover_text(format!(
                            "{} mises à jour de maille par pas de la veine, {:.0} M par seconde ({}) ; sans la mise en place",
                            cells(e.per_step),
                            self.mlups(s.fp16),
                            if self.measured_mlups.is_some_and(|m| m.0 == s.fp16) { "mesuré au dernier calcul" } else { "valeur type d'une GTX 1070" }
                        ));
                    ui.end_row();
                });
                if s.bands == 4 && mesh {
                    ui.label(RichText::new(format!("couches automatiques : {}", e.bands)).small().weak());
                }
                if let Some(g) = Tab3D::guard(&e, &s) {
                    banner(ui, theme::err(), "⚠", &g);
                    ui.checkbox(&mut self.force, "Lancer quand même");
                } else if mesh && s.setup.margin >= 0.99 && e.grids[0].0 > 0.6 * e.total && e.total > 8e6 {
                    ui.label(
                        RichText::new("La veine fait l'essentiel des mailles : pour une maquette allongée, réduire la « marge de la veine » (mise en place).")
                            .small()
                            .color(theme::warn()),
                    );
                } else if e.seconds > 1800.0 {
                    ui.label(RichText::new("Calcul long : le GPU sera occupé tout ce temps.").small().color(theme::warn()));
                }
            });
        }
    }

    // ------------------------------------------------------------------
    // Résultats
    // ------------------------------------------------------------------

    fn results(&mut self, ui: &mut egui::Ui, preview: Option<&Preview>) {
        ui.add_space(4.0);
        ui.heading("Résultats");
        let mesh = self.applied.shape == Shape3D::Mesh && self.applied.stl.is_some();
        let fr = self.shared.lock().unwrap();
        let f = &fr.frame;
        if let Some(e) = &f.error {
            banner(ui, theme::err(), "⚠", e);
        }
        if let Some(w) = &f.warning {
            banner(ui, theme::warn(), "⚠", w);
        }
        let start = f.window;
        let ready = f.ready;
        ui.label(
            RichText::new(if ready {
                format!("Moyennes sur t* = {:.0} à {:.0}", start, f.t_star)
            } else if f.hist.is_empty() {
                "Pas encore de mesure".to_string()
            } else {
                "Valeurs instantanées (régime pas encore établi)".to_string()
            })
            .small()
            .weak(),
        );
        // Moyennes de la fenêtre (ou dernière valeur).
        let avg = |v: &[[f64; 4]], k: usize| mean(v.iter().map(|p| (p[0], p[k])), start, ready);
        let main = [avg(&f.hist, 1), avg(&f.hist, 2), avg(&f.hist, 3)];
        let cm = mean(f.moment.iter().map(|p| (p[0], p[1] + self.xref * p[2])), start, ready);
        let v = |x: f64| if x.is_finite() { format!("{x:.4}") } else { "—".into() };
        let a = &f.analysis;
        egui::Frame::group(ui.style()).show(ui, |ui| {
            ui.set_width(ui.available_width());
            egui::Grid::new("s3_coef").num_columns(5).striped(true).spacing([10.0, 3.0]).show(ui, |ui| {
                ui.label("");
                ui.colored_label(theme::accent(), RichText::new("Cx").strong()).on_hover_text("Traînée");
                ui.colored_label(theme::lift(), RichText::new("Cz").strong()).on_hover_text("Portance");
                ui.label(RichText::new("Cy").strong()).on_hover_text("Effort latéral");
                ui.colored_label(theme::moment(), RichText::new("Cm").strong())
                    .on_hover_text("Moment de tangage, positif à cabrer (maquettes STL)");
                ui.end_row();
                let row = |ui: &mut egui::Ui, name: &str, tip: &str, c: [f64; 4], strong: bool| {
                    let t = RichText::new(name);
                    ui.label(if strong { t.strong() } else { t.small() }).on_hover_text(tip);
                    for x in c {
                        let t = RichText::new(v(x));
                        ui.label(if strong { t.strong() } else { t });
                    }
                    ui.end_row();
                };
                if mesh {
                    row(
                        ui,
                        "Surface",
                        "Intégration de la pression et du frottement sur les triangles de la maquette (mesure de référence de run3d)",
                        [main[0], main[1], main[2], cm],
                        true,
                    );
                    let sp = |k: usize| mean(f.split.iter().map(|p| (p[0], p[k])), start, ready);
                    row(ui, "  pression", "Part de pression de l'intégration sur la surface", [sp(1), sp(3), f64::NAN, f64::NAN], false);
                    row(ui, "  frottement", "Part du frottement pariétal", [sp(2), sp(4), f64::NAN, f64::NAN], false);
                    row(
                        ui,
                        "Échange",
                        "Échange de quantité de mouvement aux liens de paroi (à chaque pas)",
                        [avg(&f.exch, 1), avg(&f.exch, 2), avg(&f.exch, 3), f64::NAN],
                        false,
                    );
                    let (bx, bz) = box_balance(&f.boxes, start);
                    row(
                        ui,
                        "Boîtes",
                        "Bilan de quantité de mouvement sur deux boîtes (0,1 et 0,15 L autour de la maquette) : \
                         Cx extrapolé à volume nul, Cz moyen ; référence de portance. Relevé pendant la fenêtre de moyenne.",
                        [bx, bz, f64::NAN, f64::NAN],
                        false,
                    );
                } else {
                    row(
                        ui,
                        "Échange",
                        "Échange de quantité de mouvement aux liens de paroi",
                        [main[0], main[1], main[2], f64::NAN],
                        true,
                    );
                }
            });
            if ready && a.periodic {
                ui.label(
                    RichText::new(format!(
                        "Instationnaire : Cx ± {:.4}, Cz ± {:.4} · Strouhal {:.3} ({} périodes)",
                        a.cx_amp, a.cz_amp, a.st, a.periods
                    ))
                    .small(),
                );
            }
            if mesh {
                let (_, bz) = box_balance(&f.boxes, start);
                if bz.is_finite() && main[1].is_finite() && bz.abs() > 1e-3 {
                    ui.label(
                        RichText::new(format!("Écart surface / boîtes sur Cz : {:+.1} %", 100.0 * (main[1] / bz - 1.0)))
                            .small()
                            .weak(),
                    );
                }
                ui.horizontal(|ui| {
                    ui.label(RichText::new("Cm autour de").small());
                    ui.add(egui::DragValue::new(&mut self.xref).speed(0.01).range(-1.0..=2.0).max_decimals(2))
                        .on_hover_text("Point de référence en longueurs de référence depuis l'avant de la maquette (0,25 : quart de corde), au centre en z");
                    ui.label(RichText::new("L de l'avant").small());
                });
            }
        });
        // Efforts physiques.
        if let (true, true, Some(Ok(PreviewInfo { info: Some(info), .. }))) = (mesh, self.applied.phys, preview) {
            let s = &self.applied;
            let (_, rho) = s.fluid.props();
            let qs = 0.5 * rho * s.speed * s.speed * info.ref_area_m2;
            let n = |c: f64| if c.is_finite() { format!("{} N", sig(c * qs)) } else { "—".into() };
            ui.label(format!(
                "À {} m/s ({}) : traînée {}, portance {}, latéral {}",
                s.speed,
                s.fluid.label(),
                n(main[0]),
                n(main[1]),
                n(main[2])
            ));
        }
        // Références.
        ui.add_space(4.0);
        if mesh {
            ui.label(RichText::new("Comparaison").strong());
            ui.horizontal(|ui| {
                for (name, val, c) in [("Cz réf.", &mut self.refs.0, main[1]), ("Cx réf.", &mut self.refs.1, main[0])] {
                    let mut on = val.is_some();
                    if ui.checkbox(&mut on, name).on_hover_text("Valeur de référence (mesure, VSPAERO…) à comparer").changed() {
                        *val = on.then_some(if c.is_finite() { (c * 1000.0).round() / 1000.0 } else { 0.0 });
                    }
                    if let Some(r) = val {
                        ui.add(egui::DragValue::new(r).speed(0.001).max_decimals(4));
                        if c.is_finite() && r.abs() > 1e-9 {
                            ui.label(RichText::new(format!("{:+.1} %", 100.0 * (c / *r - 1.0))).small());
                        }
                    }
                }
            });
        }
        ui.label(RichText::new(reference(self.applied.shape, self.applied.re)).small());
        // Répartition en envergure.
        if mesh && !f.span.is_empty() {
            ui.add_space(4.0);
            ui.label(RichText::new("Répartition en envergure").strong())
                .on_hover_text("c·cl / c_ref par tranche (efforts par triangle de la surface), moyenne de la fenêtre");
            let n = f.span_eta.len();
            let mut acc = vec![0.0; n];
            let mut cnt = 0.0f64;
            for (t, b) in &f.span {
                if !ready || *t >= start {
                    for k in 0..n.min(b.len()) {
                        acc[k] += b[k];
                    }
                    cnt += 1.0;
                }
            }
            if !ready {
                // Dernier relevé seulement.
                if let Some((_, b)) = f.span.last() {
                    acc = b.clone();
                    cnt = 1.0;
                }
            }
            let pts: Vec<[f64; 2]> = (0..n).map(|k| [f.span_eta[k], acc[k] / cnt.max(1.0)]).filter(|p| p[1] != 0.0).collect();
            Plot::new("s3_envergure")
                .height(130.0)
                .allow_scroll(false)
                .allow_drag(false)
                .allow_zoom(false)
                .x_axis_label("η = 2y/b")
                .show(ui, |p| {
                    p.line(Line::new("c·cl/c_ref", PlotPoints::from(pts)).color(theme::lift()));
                });
        }
        // Exports.
        ui.add_space(6.0);
        drop(fr);
        ui.horizontal(|ui| {
            if ui.button("⤓ CSV").on_hover_text("Historique des efforts et résumé des mesures (Ctrl+E)").clicked() {
                self.export_csv();
            }
            if ui.button("🖼 Image").on_hover_text("Image PNG de la vue avec son échelle (Ctrl+S)").clicked() {
                self.export_image(ui.ctx());
            }
        });
        if let Some((m, t)) = &*self.message.lock().unwrap() {
            if t.elapsed() < Duration::from_secs(6) {
                ui.label(RichText::new(m).small().color(theme::ok()));
            }
        }
        ui.add_space(4.0);
        let fr = self.shared.lock().unwrap();
        let f = &fr.frame;
        egui::CollapsingHeader::new("Détails du calcul").default_open(false).show(ui, |ui| {
            if !f.grids.is_empty() {
                egui::Grid::new("s3_grilles").num_columns(3).spacing([12.0, 2.0]).show(ui, |ui| {
                    for (i, (n, d)) in f.grids.iter().enumerate() {
                        let kind = if i == 0 {
                            "veine".to_string()
                        } else if i <= f.boxes_n {
                            "boîte".to_string()
                        } else {
                            "couche".to_string()
                        };
                        ui.label(RichText::new(kind).small());
                        ui.label(RichText::new(format!("×{}", 1u64 << d)).small());
                        ui.label(RichText::new(cells(*n as f64)).small());
                        ui.end_row();
                    }
                });
                let total: usize = f.grids.iter().map(|g| g.0).sum();
                ui.label(
                    RichText::new(format!(
                        "{} mailles · {} · {} {}",
                        cells(total as f64),
                        mem((total * bytes_per_cell(f.fp16)) as f64 / 1e6),
                        if f.fp16 { "demi-précision" } else { "simple précision" },
                        format!("· grille {}×{}×{}", f.dims[0], f.dims[1], f.dims[2])
                    ))
                    .small(),
                );
            }
            ui.label(RichText::new(&f.info).small());
            ui.label(
                RichText::new(format!(
                    "{} · {:.0} M mises à jour de maille par seconde · mise en place {:.1} s · {:.1} t* en {}",
                    f.adapter,
                    f.mlups,
                    f.setup_s,
                    f.t_star,
                    crate::gui::duration(f.run_wall)
                ))
                .small()
                .weak(),
            );
            if let Some(n) = &f.note {
                ui.label(RichText::new(n).small());
            }
        });
    }

    // ------------------------------------------------------------------
    // Exports
    // ------------------------------------------------------------------

    fn export_csv(&mut self) {
        let stem = self
            .applied
            .stl
            .as_ref()
            .filter(|_| self.applied.shape == Shape3D::Mesh)
            .and_then(|f| f.path.file_stem().map(|s| s.to_string_lossy().to_string()))
            .unwrap_or_else(|| self.applied.shape.label().to_lowercase());
        let Some(path) = rfd::FileDialog::new()
            .set_file_name(format!("{stem}_a{:.0}_L{}.csv", self.applied.alpha, self.applied.l))
            .add_filter("CSV", &["csv"])
            .save_file()
        else {
            return;
        };
        let csv = {
            let fr = self.shared.lock().unwrap();
            let f = &fr.frame;
            let s = &self.applied;
            let mut out = String::new();
            out += &format!(
                "# Soufflerie de poche 3D : {} · alpha {} ° · Re {} · L {} · raffinement {} · couches {} · {}\n",
                s.stl.as_ref().filter(|_| s.shape == Shape3D::Mesh).map_or(s.shape.label().to_string(), |f| f.name.clone()),
                s.alpha,
                s.re,
                s.l,
                s.refine,
                match s.bands {
                    4 => "auto".to_string(),
                    n => n.to_string(),
                },
                if s.fp16 { "demi-précision" } else { "simple précision" }
            );
            let start = f.window;
            let ready = f.ready;
            let m = |v: &[[f64; 4]], k: usize| mean(v.iter().map(|p| (p[0], p[k])), start, ready);
            out += &format!("# moyenne sur t* {:.1} à {:.1} (établi : {})\n", start, f.t_star, if ready { "oui" } else { "non" });
            out += &format!("# reference: Cx {:.5} Cz {:.5} Cy {:.5}\n", m(&f.hist, 1), m(&f.hist, 2), m(&f.hist, 3));
            if !f.exch.is_empty() {
                out += &format!("# echange: Cx {:.5} Cz {:.5}\n", m(&f.exch, 1), m(&f.exch, 2));
                let (bx, bz) = box_balance(&f.boxes, start);
                out += &format!("# boites: Cx {bx:.5} Cz {bz:.5}\n");
                let cm = mean(f.moment.iter().map(|p| (p[0], p[1] + self.xref * p[2])), start, ready);
                out += &format!("# Cm {cm:.5} autour de {:.2} L de l'avant\n", self.xref);
            }
            out += "t_star,cx,cz,cy";
            let with_cm = f.moment.len() == f.hist.len() && !f.moment.is_empty();
            out += if with_cm { ",cm\n" } else { "\n" };
            for (i, p) in f.hist.iter().enumerate() {
                out += &format!("{:.4},{:.6},{:.6},{:.6}", p[0], p[1], p[2], p[3]);
                if with_cm {
                    let q = f.moment[i];
                    out += &format!(",{:.6}", q[1] + self.xref * q[2]);
                }
                out += "\n";
            }
            out
        };
        let msg = match std::fs::write(&path, csv) {
            Ok(()) => format!("Enregistré : {}", path.display()),
            Err(e) => format!("Écriture impossible : {e}"),
        };
        self.flash(msg);
    }

    /// Image PNG de la vue (avec son échelle) : capture de la fenêtre,
    /// recadrée sur la vue.
    fn export_image(&mut self, ctx: &egui::Context) {
        let Some(path) = rfd::FileDialog::new()
            .set_file_name("vue3d.png")
            .add_filter("Image PNG", &["png"])
            .save_file()
        else {
            return;
        };
        *self.image_req.lock().unwrap() = Some(ImageRequest {
            path,
            rect: self.view_rect,
        });
        let req = self.image_req.clone();
        let message = self.message.clone();
        let ppp = ctx.pixels_per_point();
        let c = ctx.clone();
        ctx.request_screenshot(move |img| {
            let Some(r) = req.lock().unwrap().take() else { return };
            let [w, h] = img.size;
            let x0 = ((r.rect.min.x * ppp).round().max(0.0) as usize).min(w);
            let y0 = ((r.rect.min.y * ppp).round().max(0.0) as usize).min(h);
            let x1 = ((r.rect.max.x * ppp).round() as usize).clamp(x0, w);
            let y1 = ((r.rect.max.y * ppp).round() as usize).clamp(y0, h);
            let mut rgba = Vec::with_capacity(4 * (x1 - x0) * (y1 - y0));
            for y in y0..y1 {
                for x in x0..x1 {
                    rgba.extend_from_slice(&img.pixels[x + w * y].to_array());
                }
            }
            let msg = match crate::capture::write_png(&r.path, x1 - x0, y1 - y0, &rgba) {
                Ok(()) => format!("Image enregistrée : {}", r.path.display()),
                Err(e) => format!("Image impossible : {e}"),
            };
            *message.lock().unwrap() = Some((msg, Instant::now()));
            c.request_repaint();
        });
    }

    fn help_window(&mut self, ctx: &egui::Context) {
        let mut open = self.help;
        egui::Window::new("Soufflerie 3D : raccourcis")
            .open(&mut open)
            .resizable(false)
            .collapsible(false)
            .show(ctx, |ui| {
                egui::Grid::new("s3_aide").num_columns(2).spacing([16.0, 3.0]).show(ui, |ui| {
                    for (k, d) in [
                        ("Espace", "lancer, mettre en pause, reprendre"),
                        ("Ctrl+R", "recommencer à t* = 0"),
                        ("Ctrl+Entrée", "appliquer les réglages modifiés"),
                        ("Ctrl+O", "ouvrir une maquette (ou glisser un .stl / .sgeo)"),
                        ("Ctrl+E", "exporter les efforts (CSV)"),
                        ("Ctrl+S", "exporter l'image de la vue (PNG)"),
                        ("F12", "capture de la fenêtre"),
                        ("V", "vue 3D ou coupe à plat"),
                        ("Q", "isosurfaces du critère Q"),
                        ("1 à 5", "vues 3/4, face, côté, dessus, arrière"),
                        ("0, Origine", "vue par défaut"),
                        ("glisser", "tourner la vue ; molette : rapprocher"),
                        ("double clic", "vue par défaut"),
                    ] {
                        ui.label(RichText::new(k).monospace().strong());
                        ui.label(d);
                        ui.end_row();
                    }
                });
                ui.add_space(6.0);
                ui.label(
                    RichText::new(
                        "Déroulement : maquette (forme ou STL) → écoulement (Reynolds, paroi) → maillage (préréglage, coût estimé) → Lancer. \
                         Les réglages de la physique modifiés pendant un calcul attendent « Appliquer » ; l'affichage change tout de suite.",
                    )
                    .small(),
                );
            });
        self.help = open;
    }

    // ------------------------------------------------------------------
    // Vue
    // ------------------------------------------------------------------

    fn view_bar(&mut self, ui: &mut egui::Ui) {
        let s = &mut self.settings;
        ui.horizontal(|ui| {
            ui.selectable_value(&mut s.view3d, true, "Vue 3D").on_hover_text("V");
            ui.selectable_value(&mut s.view3d, false, "Coupe à plat").on_hover_text("V");
            ui.separator();
            ui.selectable_value(&mut s.field, SliceField::Vorticity, "Vorticité")
                .on_hover_text("Composante de la vorticité normale au plan de coupe");
            ui.selectable_value(&mut s.field, SliceField::Speed, "Vitesse");
            if !s.view3d {
                ui.selectable_value(&mut s.field, SliceField::Q, "Critère Q");
            } else if s.field == SliceField::Q {
                s.field = SliceField::Vorticity;
            }
            ui.separator();
            for (p, name, tip) in [
                (Plane::Y, "Long.", "Plan longitudinal (x, z)"),
                (Plane::Z, "Horiz.", "Plan horizontal (x, y)"),
                (Plane::X, "Transv.", "Plan transversal (y, z), dans le sillage par défaut"),
            ] {
                if ui.selectable_value(&mut s.plane, p, name).on_hover_text(tip).changed() {
                    s.pos = None;
                }
            }
        });
        ui.horizontal(|ui| {
            let fr = self.shared.lock().unwrap();
            let f = &fr.frame;
            let (ax, cell) = f.plane_pos;
            let n = f.dims[ax].max(2);
            let mut rel = s.pos.unwrap_or(cell as f64 / (n - 1) as f64);
            let off = if f.l > 0 { (cell as f64 - f.center[ax]) / f.l as f64 } else { 0.0 };
            drop(fr);
            ui.label(RichText::new("Plan").weak());
            ui.spacing_mut().slider_width = 110.0;
            if ui
                .add(egui::Slider::new(&mut rel, 0.0..=1.0).show_value(false).text(format!("{off:+.1} L")))
                .on_hover_text("Position du plan, en longueurs de référence depuis le centre de la maquette")
                .changed()
            {
                s.pos = Some(rel);
            }
            if s.pos.is_some() && ui.small_button("⟲").on_hover_text("Plan au centre de la maquette (sillage en transversal)").clicked() {
                s.pos = None;
            }
            if s.view3d {
                ui.separator();
                egui::ComboBox::from_id_salt("s3_maquette")
                    .selected_text(match s.body_view {
                        BodyView::Smooth => "Maquette lisse",
                        BodyView::Grid => "Maquette + grille",
                        BodyView::Cells => "Mailles du calcul",
                    })
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut s.body_view, BodyView::Smooth, "Maquette lisse")
                            .on_hover_text("Surface exacte de la forme");
                        ui.selectable_value(&mut s.body_view, BodyView::Grid, "Maquette + grille")
                            .on_hover_text("Surface exacte, grille du calcul tracée dessus : taille des mailles");
                        ui.selectable_value(&mut s.body_view, BodyView::Cells, "Mailles du calcul").on_hover_text(
                            "Ce que voit le calcul : mailles solides (cubes) et, en orange, mailles qui ne portent que des liens de paroi mince",
                        );
                    });
                ui.checkbox(&mut s.iso, "Isosurfaces Q").on_hover_text("Tourbillons : critère Q > seuil (touche Q)");
                if s.iso {
                    ui.add(egui::Slider::new(&mut s.iso_level, 0.05..=5.0).logarithmic(true).max_decimals(2))
                        .on_hover_text("Seuil Q·(L/U)²");
                }
            }
        });
        if self.settings.view3d {
            ui.horizontal(|ui| {
                ui.label(RichText::new("Vue").weak());
                for (name, key, az, el, d) in VIEWS {
                    if ui.small_button(name).on_hover_text(format!("Touche {key}")).clicked() {
                        self.set_view((az, el, d));
                    }
                }
                if ui.small_button("⟲").on_hover_text("Vue par défaut (0, double clic)").clicked() {
                    self.set_view(CAM_DEFAULT);
                }
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(RichText::new("glisser : tourner · molette : rapprocher").small().weak());
                });
            });
        }
    }

    fn view(&mut self, ui: &mut egui::Ui) {
        self.view_bar(ui);
        if self.settings.view3d {
            self.view3d(ui)
        } else {
            self.view_slice(ui)
        }
    }

    fn view_slice(&mut self, ui: &mut egui::Ui) {
        let (resp, painter) = ui.allocate_painter(ui.available_size(), Sense::hover());
        let rect = resp.rect;
        self.view_rect = rect;
        painter.rect_filled(rect, 6.0, theme::canvas());
        let (slice, id, plane_pos, l) = {
            let s = self.shared.lock().unwrap();
            (s.frame.slice.clone(), s.frame.slice_id, s.frame.plane_pos, s.frame.l)
        };
        let Some(sl) = slice else {
            waiting(&painter, rect, "Préparation de la coupe…");
            return;
        };
        let field = self.settings.field;
        let (lo, hi, div) = match field {
            SliceField::Speed => (0.0, 1.6, false),
            SliceField::Vorticity => (-4.0, 4.0, true),
            SliceField::Q => (-1.0, 1.0, true),
        };
        if id != self.slice_id || self.texture.is_none() {
            self.slice_id = id;
            let mut px = vec![0u8; 4 * sl.w * sl.h];
            for j in 0..sl.h {
                let row = sl.h - 1 - j;
                for i in 0..sl.w {
                    let v = sl.data[i + sl.w * j];
                    let c = if v > 1e29 {
                        let c = theme::body_rgb();
                        [c[0], c[1], c[2], 255]
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
        // Image à l'échelle, centrée (mailles carrées), place pour l'échelle en bas.
        let avail = rect.shrink(12.0);
        let sc = (avail.width() / sl.w as f32).min((avail.height() - 56.0) / sl.h as f32);
        let size = Vec2::new(sl.w as f32 * sc, sl.h as f32 * sc);
        let img_rect = Rect::from_center_size(avail.center() - Vec2::new(0.0, 22.0), size);
        if let Some(t) = &self.texture {
            painter.image(t.id(), img_rect, Rect::from_min_max(Pos2::ZERO, Pos2::new(1.0, 1.0)), Color32::WHITE);
        }
        painter.rect_stroke(img_rect, 0.0, (1.0, theme::edge()), egui::StrokeKind::Outside);
        let (h_axis, v_axis) = match plane_pos.0 {
            0 => ("y", "z"),
            1 => ("x", "z"),
            _ => ("x", "y"),
        };
        let comp = ["ωx", "ωy", "ωz"][plane_pos.0];
        let (title, stops): (String, &theme::Stops) = match field {
            SliceField::Speed => ("|u| / U".into(), theme::sequential()),
            SliceField::Vorticity => (format!("{comp} · L / U"), theme::diverging()),
            SliceField::Q => ("Q · (L/U)² (positif : tourbillons)".into(), theme::diverging()),
        };
        let bar = Rect::from_min_size(
            Pos2::new(rect.right() - 232.0, rect.bottom() - 40.0),
            Vec2::new(200.0, 10.0),
        );
        colorbar(&painter, bar, &title, lo as f64, hi as f64, |t| {
            let c = theme::color(stops, t);
            Color32::from_rgb(c[0], c[1], c[2])
        });
        painter.text(
            rect.left_bottom() + Vec2::new(10.0, -10.0),
            egui::Align2::LEFT_BOTTOM,
            format!(
                "Coupe ({h_axis}, {v_axis}) en {} = {} · L = {l} mailles · écoulement de gauche à droite{}",
                ["x", "y", "z"][plane_pos.0],
                plane_pos.1,
                if plane_pos.0 == 0 { " (vu de l'amont)" } else { "" }
            ),
            egui::FontId::proportional(12.5),
            theme::ink(),
        );
    }

    /// Vue 3D : image calculée sur le GPU, cadre de la veine, caméra orbitale.
    fn view3d(&mut self, ui: &mut egui::Ui) {
        let (resp, painter) = ui.allocate_painter(ui.available_size(), Sense::click_and_drag());
        let rect = resp.rect;
        self.view_rect = rect;
        painter.rect_filled(rect, 6.0, theme::canvas());
        // Caméra : glisser et molette agissent directement, les vues
        // prédéfinies glissent en douceur vers leur cible.
        let now = Instant::now();
        let dt = now.duration_since(self.last_frame).as_secs_f64().min(0.1);
        self.last_frame = now;
        let mut moving = false;
        if resp.dragged() {
            let d = resp.drag_delta();
            self.cam_target.0 = self.settings.az - d.x as f64 * 0.4;
            self.cam_target.1 = (self.settings.el + d.y as f64 * 0.3).clamp(-85.0, 85.0);
            self.cam_target.2 = self.settings.dist;
            self.settings.az = self.cam_target.0;
            self.settings.el = self.cam_target.1;
            moving = true;
        }
        if resp.hovered() {
            let (scroll, zoom) = ui.input(|i| (i.smooth_scroll_delta.y, i.zoom_delta()));
            let f = (-(scroll as f64) * 0.0015).exp() / zoom as f64;
            if (f - 1.0).abs() > 1e-4 {
                self.cam_target.2 = (self.cam_target.2 * f).clamp(1.5, 20.0);
                self.settings.dist = self.cam_target.2;
                moving = true;
            }
        }
        if resp.double_clicked() {
            self.set_view(CAM_DEFAULT);
        }
        {
            let (az, el, d) = self.cam_target;
            let s = &mut self.settings;
            let mut daz = (az - s.az) % 360.0;
            if daz > 180.0 {
                daz -= 360.0;
            } else if daz < -180.0 {
                daz += 360.0;
            }
            let (del, dd) = (el - s.el, d - s.dist);
            if daz.abs() + del.abs() > 0.05 || dd.abs() > 0.005 {
                let k = 1.0 - (-dt * 9.0).exp();
                s.az += daz * k;
                s.el += del * k;
                s.dist += dd * k;
                moving = true;
                ui.ctx().request_repaint();
            } else {
                (s.az, s.el, s.dist) = (az, el, d);
            }
        }
        // Image à ~900 pixels de large au plus (le lancer de rayons coûte par
        // pixel), à demi-résolution pendant qu'on tourne ou rapproche la vue.
        let k = (900.0 / rect.width()).min(1.0) * if moving { 0.5 } else { 1.0 };
        let want = (((rect.width() * k) as usize).max(64), ((rect.height() * k) as usize).max(64));
        self.settings.img = want;
        let (image, id, cam, dims, l, plane_pos, center, status, t_star) = {
            let s = self.shared.lock().unwrap();
            let f = &s.frame;
            (f.image.clone(), f.slice_id, f.camera, f.dims, f.l, f.plane_pos, f.center, f.status.clone(), f.t_star)
        };
        let Some(img) = image else {
            waiting(&painter, rect, "Préparation du calcul 3D…");
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
            painter.image(t.id(), img_rect, Rect::from_min_max(Pos2::ZERO, Pos2::new(1.0, 1.0)), Color32::WHITE);
        }
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
            // Cadre de la veine (arêtes devant l'œil seulement).
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
            let stroke = (1.0, theme::alpha(theme::edge(), 90));
            for (i, j) in edges {
                if let (Some(a), Some(b)) = (to(pts[i]), to(pts[j])) {
                    p.line_segment([a, b], stroke);
                }
            }
            // Flèche du vent en amont de la maquette.
            let lf = l as f64;
            let a0 = [center[0] - 2.0 * lf, center[1], center[2] + 0.8 * lf];
            let a1 = [center[0] - 1.3 * lf, center[1], center[2] + 0.8 * lf];
            if let (Some(a), Some(b)) = (to(a0), to(a1)) {
                if (b - a).length() > 8.0 && (b - a).length() < 400.0 {
                    p.arrow(a, b - a, (2.0, theme::alpha(theme::ink(), 180)));
                    p.text(a, egui::Align2::RIGHT_CENTER, "V∞ ", egui::FontId::proportional(13.0), theme::ink());
                }
            }
            // Trièdre des axes de la veine (x : vers l'aval, z : vers le haut).
            let o = Pos2::new(rect.left() + 38.0, rect.bottom() - 38.0);
            let dot = |u: [f64; 3], v: [f64; 3]| u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
            for (k, (name, col)) in [("x", theme::accent()), ("y", theme::ok()), ("z", theme::lift())].into_iter().enumerate() {
                let mut e = [0.0; 3];
                e[k] = 1.0;
                let d = Vec2::new(dot(e, c.right) as f32, -dot(e, c.up) as f32) * 24.0;
                painter.arrow(o, d, (2.0, col));
                painter.text(o + d * 1.25, egui::Align2::CENTER_CENTER, name, egui::FontId::proportional(12.0), col);
            }
        }
        // Échelles de couleur.
        let st = &self.settings;
        let comp = ["ωx", "ωy", "ωz"][plane_pos.0];
        let bar = Rect::from_min_size(Pos2::new(rect.right() - 232.0, rect.bottom() - 34.0), Vec2::new(200.0, 10.0));
        match st.field {
            SliceField::Speed => colorbar(&painter, bar, "coupe : |u| / U", 0.0, 1.6, theme::render_seq),
            _ => colorbar(&painter, bar, &format!("coupe : {comp} · L / U"), -4.0, 4.0, |t| theme::render_div(2.0 * t - 1.0)),
        }
        if st.iso {
            // Teinte des isosurfaces : seqMap(0,15 + 0,75 |u| / 1,3 U).
            let bar2 = bar.translate(Vec2::new(0.0, -44.0));
            colorbar(
                &painter,
                bar2,
                &format!("isosurfaces Q = {:.2} (U/L)² : |u| / U", st.iso_level),
                0.0,
                1.3,
                |t| theme::render_seq(0.15 + 0.75 * t),
            );
        }
        painter.text(
            rect.left_top() + Vec2::new(10.0, 8.0),
            egui::Align2::LEFT_TOP,
            format!("L = {l} mailles"),
            egui::FontId::proportional(12.0),
            theme::alpha(theme::ink(), 160),
        );
        // Calcul pas encore lancé : invitation.
        if !self.settings.running && t_star == 0.0 && status == "Prêt" {
            let r = Rect::from_center_size(rect.center_top() + Vec2::new(0.0, 46.0), Vec2::new(330.0, 54.0));
            painter.rect_filled(r, 10.0, theme::alpha(ui.visuals().panel_fill, 230));
            painter.rect_stroke(r, 10.0, (1.0, theme::alpha(theme::ok(), 160)), egui::StrokeKind::Inside);
            painter.text(
                r.center() - Vec2::new(0.0, 9.0),
                egui::Align2::CENTER_CENTER,
                "Prêt : ▶ Lancer (Espace)",
                egui::FontId::proportional(16.0),
                theme::ink(),
            );
            painter.text(
                r.center() + Vec2::new(0.0, 13.0),
                egui::Align2::CENTER_CENTER,
                "vérifier la mise en place, puis lancer",
                egui::FontId::proportional(12.0),
                theme::alpha(theme::ink(), 160),
            );
        }
    }

    fn history(&self, ui: &mut egui::Ui) {
        let s = self.shared.lock().unwrap();
        let f = &s.frame;
        let mesh = !f.exch.is_empty();
        let t_end = self.settings.t_end;
        let t_now = f.hist.last().map_or(0.0, |p| p[0]);
        let t_max = if t_end > 0.0 { t_end.max(t_now) } else { t_now.max(10.0) };
        let t_min = if t_end > 0.0 { 0.0 } else { (t_max - 40.0).max(0.0) };
        let pts = |v: &[[f64; 4]], k: usize| -> Vec<[f64; 2]> { v.iter().filter(|p| p[0] >= t_min).map(|p| [p[0], p[k]]).collect() };
        let start = f.window;
        let ready = f.ready;
        let xref = self.xref;
        Plot::new("s3_hist")
            .height(ui.available_height().max(80.0))
            .legend(Legend::default().position(egui_plot::Corner::RightTop))
            .allow_scroll(false)
            .x_axis_label("t* (temps convectifs) · zone verte : fenêtre de moyenne, tirets : moyennes")
            .include_x(t_min)
            .include_x(t_max)
            .include_y(0.0)
            .include_y(0.5)
            .show(ui, |p| {
                if ready {
                    p.span(
                        Span::new("", start..=t_now)
                            .fill(theme::alpha(theme::ok(), 28))
                            .border_width(0.0),
                    );
                }
                if t_end > 0.0 {
                    p.vline(VLine::new("", t_end).color(theme::alpha(theme::ink(), 90)).style(LineStyle::dashed_dense()));
                }
                p.line(Line::new("Cx", PlotPoints::from(pts(&f.hist, 1))).color(theme::accent()).width(1.6));
                p.line(Line::new("Cz", PlotPoints::from(pts(&f.hist, 2))).color(theme::lift()).width(1.6));
                if mesh {
                    let cm: Vec<[f64; 2]> =
                        f.moment.iter().filter(|q| q[0] >= t_min).map(|q| [q[0], q[1] + xref * q[2]]).collect();
                    p.line(Line::new("Cm", PlotPoints::from(cm)).color(theme::moment()));
                    p.line(
                        Line::new("Cz échange", PlotPoints::from(pts(&f.exch, 2)))
                            .color(theme::alpha(theme::lift(), 110))
                            .style(LineStyle::dotted_loose()),
                    );
                }
                if ready {
                    for (k, c) in [(1, theme::accent()), (2, theme::lift())] {
                        let m = mean(f.hist.iter().map(|q| (q[0], q[k])), start, true);
                        if m.is_finite() {
                            p.line(
                                Line::new("", PlotPoints::from(vec![[start, m], [t_now, m]]))
                                    .color(c)
                                    .style(LineStyle::dashed_loose())
                                    .width(1.0),
                            );
                        }
                    }
                }
            });
    }
}

/// Groupe de réglages repliable.
fn section(ui: &mut egui::Ui, title: &str, open: bool, add: impl FnOnce(&mut egui::Ui)) {
    egui::CollapsingHeader::new(RichText::new(title).strong().size(15.0))
        .default_open(open)
        .show(ui, |ui| {
            add(ui);
            ui.add_space(4.0);
        });
    ui.separator();
}

/// Bandeau d'avertissement ou d'erreur.
fn banner(ui: &mut egui::Ui, color: Color32, icon: &str, text: &str) {
    egui::Frame::new()
        .fill(theme::alpha(color, 30))
        .stroke((1.0, theme::alpha(color, 120)))
        .corner_radius(6)
        .inner_margin(6)
        .show(ui, |ui| {
            ui.set_width(ui.available_width());
            ui.horizontal_wrapped(|ui| {
                ui.colored_label(color, RichText::new(icon).strong());
                ui.label(RichText::new(text).color(color));
            });
        });
}

/// Message d'attente au centre de la vue.
fn waiting(painter: &egui::Painter, rect: Rect, text: &str) {
    painter.text(rect.center(), egui::Align2::CENTER_CENTER, text, egui::FontId::proportional(16.0), Color32::GRAY);
}

/// Échelle de couleur horizontale avec titre et bornes.
fn colorbar(painter: &egui::Painter, r: Rect, title: &str, lo: f64, hi: f64, color: impl Fn(f32) -> Color32) {
    let bg = Rect::from_min_max(r.min - Vec2::new(8.0, 22.0), r.max + Vec2::new(8.0, 18.0));
    painter.rect_filled(bg, 6.0, theme::alpha(theme::canvas(), 210));
    let n = 48;
    for i in 0..n {
        let x0 = r.left() + r.width() * i as f32 / n as f32;
        let x1 = r.left() + r.width() * (i + 1) as f32 / n as f32 + 0.5;
        painter.rect_filled(
            Rect::from_min_max(Pos2::new(x0, r.top()), Pos2::new(x1, r.bottom())),
            0.0,
            color((i as f32 + 0.5) / n as f32),
        );
    }
    painter.rect_stroke(r, 0.0, (0.5, theme::alpha(theme::ink(), 120)), egui::StrokeKind::Outside);
    let font = egui::FontId::proportional(11.5);
    let ink = theme::ink();
    painter.text(r.left_top() - Vec2::new(0.0, 3.0), egui::Align2::LEFT_BOTTOM, title, font.clone(), ink);
    let fmt = |v: f64| format!("{}", (v * 10.0).round() / 10.0).replace('.', ",").replace('-', "−");
    painter.text(r.left_bottom() + Vec2::new(0.0, 2.0), egui::Align2::LEFT_TOP, fmt(lo), font.clone(), ink);
    painter.text(r.center_bottom() + Vec2::new(0.0, 2.0), egui::Align2::CENTER_TOP, fmt(0.5 * (lo + hi)), font.clone(), ink);
    painter.text(r.right_bottom() + Vec2::new(0.0, 2.0), egui::Align2::RIGHT_TOP, fmt(hi), font, ink);
}

/// Moyenne des valeurs relevées depuis `start` (régime établi), sinon la
/// dernière valeur.
fn mean(v: impl Iterator<Item = (f64, f64)>, start: f64, ready: bool) -> f64 {
    let (mut s, mut n, mut last) = (0.0, 0.0, f64::NAN);
    for (t, x) in v {
        last = x;
        if ready && t >= start && x.is_finite() {
            s += x;
            n += 1.0;
        }
    }
    if ready && n > 0.0 {
        s / n
    } else {
        last
    }
}

/// Bilan sur boîtes : Cx extrapolé à volume nul (le léger gradient de
/// pression de la veine pousse aussi le fluide de la boîte, ∝ volume), Cz
/// moyen des deux boîtes ; comme `run3d`.
fn box_balance(b: &[[f64; 7]], start: f64) -> (f64, f64) {
    let w: Vec<&[f64; 7]> = b.iter().filter(|p| p[0] >= start).collect();
    if w.is_empty() {
        return (f64::NAN, f64::NAN);
    }
    let n = w.len() as f64;
    let m = |k: usize| w.iter().map(|p| p[k]).sum::<f64>() / n;
    let (cx1, cz1, v1, cx2, cz2, v2) = (m(1), m(2), m(3), m(4), m(5), m(6));
    let cx0 = if (v2 - v1).abs() > 1e-9 { cx1 - v1 * (cx2 - cx1) / (v2 - v1) } else { cx1 };
    (cx0, 0.5 * (cz1 + cz2))
}

/// Nombre de mailles lisible.
fn cells(n: f64) -> String {
    if n >= 1e6 {
        format!("{} M", sig(n / 1e6)).replace('.', ",")
    } else if n >= 1e3 {
        format!("{:.0} k", n / 1e3)
    } else {
        format!("{n:.0}")
    }
}

/// Mémoire lisible (Mo).
fn mem(mb: f64) -> String {
    if mb >= 1000.0 {
        format!("{} Go", sig(mb / 1e3)).replace('.', ",")
    } else {
        format!("{mb:.0} Mo")
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

/// Résumé de la maquette STL : dimensions, blocage, diagnostic.
fn mesh_info(ui: &mut egui::Ui, s: &Settings, info: Option<&MeshInfo>, preview: Option<&Preview>) {
    let m = s.setup.unit.meters();
    let u = s.setup.unit.label();
    match (preview, info) {
        (Some(Ok(_)), Some(i)) => {
            let [a, b, c] = i.size_m.map(|v| sig(v / m));
            ui.label(
                RichText::new(format!(
                    "{} triangles · {a} × {b} × {c} {u} · réf. {} {u}, {} {u}² · blocage {:.1} %",
                    i.triangles,
                    sig(i.ref_length_m / m),
                    sig(i.ref_area_m2 / (m * m)),
                    100.0 * i.blockage
                ))
                .small(),
            )
            .on_hover_text("Longueur, largeur, hauteur ; longueur et surface de référence ; blocage : maître-couple / section de la veine");
            if i.blockage > 0.05 {
                banner(ui, theme::warn(), "⚠", "Blocage de la veine supérieur à 5 % : coefficients surestimés.");
            }
            let d = i.diag.summary();
            if !d.is_empty() {
                let more = if i.diag.watertight() {
                    ""
                } else {
                    " : maillage ouvert, appartenance votée sur trois rayons (x, y, z)"
                };
                banner(ui, theme::warn(), "⚠", &format!("{d}{more}"));
            }
        }
        (Some(Err(e)), _) => banner(ui, theme::err(), "⚠", e),
        _ => {}
    }
}

/// Mise en place de la maquette STL : dérapage, axes, unité, références.
fn mesh_controls(ui: &mut egui::Ui, s: &mut Settings, info: Option<&MeshInfo>) {
    let st = &mut s.setup;
    ui.add(egui::Slider::new(&mut st.beta_deg, -30.0..=30.0).step_by(1.0).text("dérapage β (°)"))
        .on_hover_text("Positif : vent venant de la droite de la maquette");
    let axis_box = |ui: &mut egui::Ui, id: &str, v: &mut Axis| {
        egui::ComboBox::from_id_salt(id).width(52.0).selected_text(v.label()).show_ui(ui, |ui| {
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
    ui.add(egui::Slider::new(&mut st.margin, 0.3..=1.0).step_by(0.05).text("marge de la veine")).on_hover_text(
        "Facteur sur les marges de la veine autour de la maquette (1 : 3,5 M en amont, 7,5 M en aval, 2,5 M de côté,          M : la plus grande de L et de la dimension de la maquette). Une aile très allongée tient dans une veine plus serrée          (0,4 pour l'aile AR 4 du banc) : surveiller le blocage affiché.",
    );
    ui.checkbox(&mut st.smooth, "Lisser les facettes").on_hover_text(
        "Facettes raffinées sur la surface lisse qui passe par les sommets du fichier (arêtes vives gardées) : un bord d'attaque décrit par quelques facettes retrouve son rayon",
    );
    ui.checkbox(&mut st.mesh_edges, "Arêtes vives du maillage").on_hover_text(
        "La loi de paroi n'est coupée qu'aux arêtes vives du maillage (bord de fuite, saumon plat), et non autour d'un bord d'attaque plus fin qu'une maille",
    );
    let u = st.unit.label();
    let m = st.unit.meters();
    ui.horizontal(|ui| {
        let mut on = st.ref_length.is_some();
        ui.checkbox(&mut on, "Longueur de réf.")
            .on_hover_text("Corde moyenne pour une aile ; sinon longueur hors tout");
        if !on {
            st.ref_length = None;
            ui.label(RichText::new("longueur hors tout").weak());
        } else if st.ref_length.is_none() {
            st.ref_length = info.map(|i| i.ref_length_m / m);
        }
        if let Some(v) = &mut st.ref_length {
            ui.add(egui::DragValue::new(v).speed(1.0).range(1e-6..=1e9).suffix(format!(" {u}")));
        }
    });
    ui.horizontal(|ui| {
        let mut on = st.ref_area.is_some();
        ui.checkbox(&mut on, "Surface de réf.")
            .on_hover_text("Surface alaire pour une aile ; sinon maître-couple");
        if !on {
            st.ref_area = None;
            ui.label(RichText::new("maître-couple").weak());
        } else if st.ref_area.is_none() {
            st.ref_area = info.map(|i| i.frontal_m2 / (m * m));
        }
        if let Some(v) = &mut st.ref_area {
            ui.add(egui::DragValue::new(v).speed(10.0).range(1e-9..=1e15).suffix(format!(" {u}²")));
        }
    });
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
                "Référence sphère (Schiller–Naumann, Clift) : Cx ≈ {cd:.3} ; {regime}. Le confinement latéral (veine 6D périodique) augmente Cx de quelques %."
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
        Shape3D::Mesh => "Maquette maillée : coefficients rapportés à la surface de référence (maître-couple par défaut) ; Cm à la longueur de référence.".into(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn moyenne_de_la_fenetre() {
        let v = [(0.0, 1.0), (5.0, 2.0), (10.0, 4.0), (15.0, 6.0)];
        assert_eq!(mean(v.into_iter(), 10.0, true), 5.0);
        // Régime pas encore établi : dernière valeur.
        assert_eq!(mean(v.into_iter(), 10.0, false), 6.0);
    }

    /// Bilan sur boîtes : Cx extrapolé à volume nul, comme run3d.
    #[test]
    fn bilan_sur_boites() {
        let b = [[12.0, 0.0723, 0.27, 1.0, 0.0719, 0.268, 2.0]];
        let (cx, cz) = box_balance(&b, 10.0);
        assert!((cx - 0.0727).abs() < 1e-9);
        assert!((cz - 0.269).abs() < 1e-9);
        assert!(box_balance(&b, 20.0).0.is_nan());
    }

    /// Estimation des grilles boîtes : même emprise que
    /// `Gpu3D::with_refinement_edges` (aile AR 4 du banc, L = 16, ×2 :
    /// 2,15 M + 0,78 M mailles, README « Vitesse de calcul 3D » : 2,91 M).
    #[test]
    fn estimation_des_grilles() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../banc/wing_ar4.stl");
        let Ok(tris) = lbm::stl::read_stl(&path) else { return };
        let mut t = Tab3D::new();
        t.settings.shape = Shape3D::Mesh;
        t.settings.setup = MeshSetup {
            unit: Unit::M,
            ref_length: Some(1.0),
            ref_area: Some(4.0),
            margin: 0.4,
            ..Default::default()
        };
        t.settings.alpha = 4.0;
        t.settings.stl = Some(StlFile {
            name: "wing_ar4.stl".into(),
            path: path.clone(),
            tris: Arc::new(tris),
            id: 1,
            wetted: 8.0,
        });
        t.settings.refine = 1;
        let ctx = egui::Context::default();
        let p = t.mesh_preview(&ctx).unwrap().unwrap();
        let e = t.estimate(&t.settings.clone(), &p);
        let real = 2.91e6;
        assert!((e.total - real).abs() / real < 0.03, "{} mailles", e.total);
    }
}
