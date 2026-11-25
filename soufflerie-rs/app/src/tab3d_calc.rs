//! Soufflerie 3D : réglages, état publié et fil de calcul (GPU).
//!
//! Le fil construit la veine, avance par lots d'environ 30 ms, relève les
//! efforts et publie l'image (vue 3D ou coupe) à chaque lot. Les mesures sont
//! celles de `soufflerie run3d` : intégration sur la surface (maquettes
//! STL), échange de quantité de mouvement aux parois, bilan sur deux boîtes
//! autour de la maquette (référence de portance), moment de tangage et
//! répartition en envergure. Le calcul lui-même n'est pas modifié.

use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use lbm::bodies3d::{Body3D, Shape3D};
use lbm::cpu::Ramp;
use lbm::fluid::Fluid;
use lbm::gpu3d::{BodyView, Camera, Gpu3D, RenderOptions, Slice, SliceField};
use lbm::stats::{analyze, Analysis};
use lbm::stl::{MeshSetup, Tri};

use crate::theme;

pub const U: f32 = 0.1;
pub const RAMP_T: f64 = 4.0;
pub const SETTLE_T: f64 = 6.0;
pub const WINDOW_T: f64 = 20.0;
/// Tranches de la répartition en envergure.
pub const SPAN_BINS: usize = 16;
/// Épaisseur des couches pariétales (mailles de la couche).
pub const BAND: f64 = 16.0;
/// Épaisseur minimale des couches en fraction de la longueur de référence
/// (`--bande-corde 0.125`, règle retenue).
pub const BAND_CHORD: f64 = 0.125;

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
    pub fn axis(self) -> usize {
        match self {
            Plane::X => 0,
            Plane::Y => 1,
            Plane::Z => 2,
        }
    }
}

/// Maquette STL chargée (comparée par numéro de chargement).
#[derive(Clone)]
pub struct StlFile {
    pub name: String,
    pub path: std::path::PathBuf,
    pub tris: Arc<Vec<Tri>>,
    pub id: u64,
    /// Surface mouillée (unités du fichier au carré).
    pub wetted: f64,
}

impl PartialEq for StlFile {
    fn eq(&self, o: &StlFile) -> bool {
        self.id == o.id
    }
}

/// Réglages transmis au fil de calcul. Les réglages « physiques » (maquette,
/// écoulement, maillage) ne passent au calcul qu'à l'application ; ceux de
/// l'affichage et de la marche passent tout de suite (`take_view`).
#[derive(Clone, PartialEq)]
pub struct Settings {
    pub shape: Shape3D,
    /// Maquette STL et sa mise en place (incidence : `alpha`).
    pub stl: Option<StlFile>,
    pub setup: MeshSetup,
    /// Échelle physique de la maquette STL : le Reynolds découle de la vitesse.
    pub phys: bool,
    pub speed: f64,
    pub fluid: Fluid,
    pub re: f64,
    pub l: usize,
    pub alpha: f64,
    pub fp16: bool,
    /// Incrémenté pour relancer le calcul.
    pub reset: u64,
    /// Grilles imbriquées autour de la maquette (0 à 3).
    pub refine: usize,
    /// Couches pariétales fines (maquette maillée, au-dessus d'au moins une
    /// grille boîte) : 0 aucune, 1 à 3 couches, 4 automatique (`--bande auto`).
    pub bands: u8,
    /// Épaisseur minimale des couches : 1/8 de la longueur de référence
    /// (`--bande-corde 0.125`).
    pub band_chord: bool,
    /// Paroi : 0 auto (loi de paroi dès Re = 5·10⁴), 1 loi de paroi, 2 sans.
    pub wall: u8,
    /// Transition imposée (fraction de la corde locale ; 0 : turbulent).
    pub xtr: f64,
    // --- marche et affichage (appliqués tout de suite) ---
    pub running: bool,
    /// Fin du calcul (temps convectifs) ; 0 : sans fin.
    pub t_end: f64,
    pub field: SliceField,
    pub plane: Plane,
    /// Position relative du plan (0 à 1) ; None : maquette ou sillage.
    pub pos: Option<f64>,
    /// Vue 3D (maquette, plan de coupe, isosurfaces) ou coupe plane.
    pub view3d: bool,
    /// Caméra : azimut, élévation (°), distance (en L).
    pub az: f64,
    pub el: f64,
    pub dist: f64,
    pub iso: bool,
    /// Maquette : lisse, avec la grille, ou mailles du calcul.
    pub body_view: BodyView,
    /// Seuil des isosurfaces, en (U/L)².
    pub iso_level: f64,
    /// Taille de l'image 3D (pixels).
    pub img: (usize, usize),
}

impl Settings {
    /// Recopie les réglages de marche et d'affichage de `v`.
    pub fn take_view(&mut self, v: &Settings) {
        self.running = v.running;
        self.t_end = v.t_end;
        self.field = v.field;
        self.plane = v.plane;
        self.pos = v.pos;
        self.view3d = v.view3d;
        self.az = v.az;
        self.el = v.el;
        self.dist = v.dist;
        self.iso = v.iso;
        self.body_view = v.body_view;
        self.iso_level = v.iso_level;
        self.img = v.img;
    }

    /// Mêmes réglages physiques (maquette, écoulement, maillage) ?
    pub fn same_physics(&self, o: &Settings) -> bool {
        let mut a = self.clone();
        a.take_view(o);
        a == *o
    }
}

/// Loi de paroi active pour ces réglages.
pub fn wall_model(s: &Settings) -> bool {
    match s.wall {
        1 => true,
        2 => false,
        _ => s.re >= lbm::tunnel::RE_SA_AUTO,
    }
}

/// Épaisseurs des couches pariétales (mailles de chaque couche), comme
/// `run3d --bande n|auto [--bande-corde 0.125]`.
pub fn band_thickness(n: usize, l: usize, depth: usize, chord: bool) -> Vec<f64> {
    let mut b = vec![BAND; n];
    if chord {
        for (i, v) in b.iter_mut().enumerate() {
            let lk = (l << (depth + 1 + i)) as f64;
            *v = v.max((BAND_CHORD * lk).round());
        }
    }
    lbm::gpu3d::nested_bands(&b)
}

/// Ce que le fil de calcul publie.
#[derive(Default)]
pub struct Frame {
    pub slice: Option<Arc<Slice>>,
    /// Vue 3D : largeur, hauteur, RGBA, et la caméra qui l'a produite.
    pub image: Option<Arc<(usize, usize, Vec<u8>)>>,
    pub camera: Option<Camera>,
    pub slice_id: u64,
    /// Plan effectivement affiché (axe, maille).
    pub plane_pos: (usize, usize),
    pub dims: [usize; 3],
    pub center: [f64; 3],
    pub l: usize,
    pub t_star: f64,
    /// Efforts de référence : t*, Cx, Cz, Cy (surface pour une maquette
    /// STL, échange de quantité de mouvement sinon).
    pub hist: Vec<[f64; 4]>,
    /// Maquette STL : échange de quantité de mouvement aux parois, à chaque lot.
    pub exch: Vec<[f64; 4]>,
    /// Maquette STL : t*, Cx pression, Cx frottement, Cz pression, Cz frottement.
    pub split: Vec<[f64; 5]>,
    /// Maquette STL : moment de tangage, t*, a, b ; Cm = a + xref·b
    /// (xref en longueurs de référence depuis l'avant de la maquette).
    pub moment: Vec<[f64; 3]>,
    /// Bilan sur boîtes (marges 0,1 et 0,15 longueur) : t*, Cx₁, Cz₁, V₁, Cx₂, Cz₂, V₂.
    pub boxes: Vec<[f64; 7]>,
    /// Répartition en envergure (c·cl/c_ref par tranche) relevée à t*.
    pub span: Vec<(f64, Vec<f64>)>,
    /// Centres des tranches en η = y / demi-envergure.
    pub span_eta: Vec<f64>,
    /// Début de la fenêtre de moyenne (t*).
    pub window: f64,
    pub analysis: Analysis,
    pub ready: bool,
    pub status: String,
    /// Le calcul avance (ni en pause, ni terminé, ni arrêté).
    pub advancing: bool,
    pub finished: bool,
    pub building: bool,
    /// Vitesse d'avancement : temps convectifs par seconde.
    pub t_rate: f64,
    pub info: String,
    pub error: Option<String>,
    /// Avertissement sur la maquette (parties perdues à la voxelisation).
    pub warning: Option<String>,
    /// Remarques de la mise en place (couches automatiques).
    pub note: Option<String>,
    pub mlups: f64,
    /// Mailles de chaque grille et son niveau.
    pub grids: Vec<(usize, u32)>,
    /// Grilles boîtes (les suivantes sont des couches pariétales).
    pub boxes_n: usize,
    pub fp16: bool,
    pub adapter: String,
    /// Durée de la mise en place (s).
    pub setup_s: f64,
    /// Mise en place en cours : étape, fraction faite, début.
    pub setup_stage: String,
    pub setup_frac: f64,
    pub setup_t0: Option<Instant>,
    /// Durée de calcul depuis t* = 0 (s).
    pub run_wall: f64,
}

/// Ce qui fixe l'image : refaite tout de suite si cela change.
#[derive(Clone, PartialEq)]
struct ViewKey {
    view3d: bool,
    cam: [f64; 3],
    img: (usize, usize),
    plane: Plane,
    pos: Option<f64>,
    field: SliceField,
    iso: Option<f64>,
    body_view: BodyView,
    dark: bool,
}

impl ViewKey {
    fn of(s: &Settings) -> ViewKey {
        ViewKey {
            view3d: s.view3d,
            cam: [s.az, s.el, s.dist],
            img: s.img,
            plane: s.plane,
            pos: s.pos,
            field: s.field,
            iso: s.iso.then_some(s.iso_level),
            body_view: s.body_view,
            dark: theme::dark(),
        }
    }
}

pub struct Shared {
    pub settings: Settings,
    pub frame: Frame,
    /// Fil de calcul actif (onglet affiché).
    pub active: bool,
}

fn rgbf(c: eframe::egui::Color32) -> [f32; 3] {
    [c.r() as f32 / 255.0, c.g() as f32 / 255.0, c.b() as f32 / 255.0]
}

/// Fil de calcul.
pub fn worker(shared: Arc<Mutex<Shared>>, repaint: impl Fn() + Send + 'static) {
    let rp_arc = Arc::new(repaint);
    let repaint = {
        let r = rp_arc.clone();
        move || (*r)()
    };
    // Étapes de la mise en place, publiées au fil de l'eau.
    {
        let (sh, rp) = (shared.clone(), rp_arc.clone());
        lbm::progress::set_sink(move |p| {
            if let Ok(mut s) = sh.lock() {
                s.frame.setup_stage = p.stage.clone();
                s.frame.setup_frac = p.frac;
            }
            (*rp)();
        });
    }
    let mut sim: Option<Gpu3D> = None;
    let mut built: Option<Settings> = None;
    let mut body: Option<Body3D> = None;
    let mut hist: Vec<[f64; 4]> = Vec::new();
    let mut exch: Vec<[f64; 4]> = Vec::new();
    let mut split: Vec<[f64; 5]> = Vec::new();
    let mut moment: Vec<[f64; 3]> = Vec::new();
    let mut boxes: Vec<[f64; 7]> = Vec::new();
    let mut span: Vec<(f64, Vec<f64>)> = Vec::new();
    let mut span_eta: Vec<f64> = Vec::new();
    let mut t_change = 0.0f64;
    let mut steps = 8usize;
    let mut slice_id = 0u64;
    let mut perf = (Instant::now(), 0u64, 0.0f64, 0.0f64);
    let mut t_rate = 0.0f64;
    // Après une divergence, le calcul attend « Recommencer ».
    let mut halted = false;
    // Maquette maillée : efforts intégrés sur la surface, relevés tous les
    // demi-temps convectifs (lecture du champ entier).
    let mut next_surf = 0.0f64;
    // Image : refaite à chaque lot si la vue change, sinon au plus toutes les
    // 120 ms pendant le calcul (le lancer de rayons et la relecture coûtent
    // autant que plusieurs lots).
    let mut last_render: Option<(Instant, ViewKey)> = None;
    // Durée de calcul (s) depuis t* = 0.
    let mut run_wall = 0.0f64;
    let mut announced = false;
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
                    || b.refine != st.refine
                    || b.bands != st.bands
                    || b.band_chord != st.band_chord
                    || wall_model(b) != wall_model(&st)
                    || (st.shape == Shape3D::Mesh && (b.stl != st.stl || b.setup != st.setup))
            }
        };
        if need {
            sim = None;
            {
                let mut s = shared.lock().unwrap();
                let fr = &mut s.frame;
                fr.status = "Préparation".into();
                fr.building = true;
                fr.setup_stage = "préparation de la maquette".into();
                fr.setup_frac = 0.0;
                fr.setup_t0 = Some(Instant::now());
                fr.advancing = false;
                fr.finished = false;
                fr.error = None;
                fr.note = None;
                fr.hist.clear();
                fr.exch.clear();
                fr.split.clear();
                fr.moment.clear();
                fr.boxes.clear();
                fr.span.clear();
                fr.t_star = 0.0;
                fr.ready = false;
                fr.analysis = Analysis::default();
            }
            repaint();
            let t_setup = Instant::now();
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
            let sa = wall_model(&st);
            // Couches pariétales (maquette maillée, au-dessus d'une grille
            // boîte), comme `run3d --bande n|auto [--bande-corde 0.125]`.
            let mut note = None;
            let bands_of = |b: &Body3D, note: &mut Option<String>| -> Vec<f64> {
                if b.mesh.is_none() || st.refine == 0 || st.bands == 0 {
                    return Vec::new();
                }
                let n = match st.bands {
                    4 => {
                        lbm::progress::prelude("couches automatiques : rayon du bord d'attaque", 0.8);
                        let (n, r, corner) = lbm::gpu3d::auto_band_count(b, st.refine);
                        *note = Some(format!(
                            "Couches automatiques : rayon de bord d'attaque minimal {:.2} mailles de la veine → {n} couche(s){}",
                            r,
                            if corner > 60.0 {
                                format!(" · bord d'attaque à facettes (coin de {corner:.0}°) : « Lisser les facettes » conseillé")
                            } else {
                                String::new()
                            }
                        ));
                        n
                    }
                    k => k as usize,
                };
                band_thickness(n, st.l, st.refine, st.band_chord)
            };
            match b.and_then(|b| {
                let bands = bands_of(&b, &mut note);
                Gpu3D::with_refinement_edges(&b, nu, U, st.fp16, sa, st.refine, false, false, &bands).map(|g| (g, b))
            }) {
                Ok((mut g, b)) => {
                    g.set_transition(st.xtr);
                    let lost = mesh_info.as_ref().and_then(|i: &lbm::stl::MeshInfo| i.volume_lost(g.solid_count, st.l));
                    let mut s = shared.lock().unwrap();
                    s.frame.warning = match lost {
                        _ if b.mesh.is_some() && g.wall_count == 0 => {
                            Some("Aucun lien de paroi : maquette invisible à cette résolution, l'augmenter.".into())
                        }
                        Some(v) if v > 0.2 => Some(format!(
                            "{:.0} % du volume tient dans des parties plus fines qu'une maille, traitées en parois minces (moins précises) : augmenter la résolution pour les résoudre.",
                            100.0 * v
                        )),
                        _ => None,
                    };
                    s.frame.note = note;
                    s.frame.grids = g.grid_cells();
                    s.frame.boxes_n = st.refine.min(s.frame.grids.len() - 1);
                    s.frame.fp16 = g.fp16;
                    s.frame.adapter = g.adapter_name().to_string();
                    s.frame.setup_s = t_setup.elapsed().as_secs_f64();
                    drop(s);
                    // Tranches de la répartition en envergure (maquette STL).
                    span_eta = (0..SPAN_BINS).map(|k| -1.0 + (k as f64 + 0.5) * 2.0 / SPAN_BINS as f64).collect();
                    sim = Some(g);
                    body = Some(b);
                }
                Err(e) => {
                    let mut s = shared.lock().unwrap();
                    s.frame.error = Some(format!("Calcul 3D impossible : {e}"));
                    s.frame.status = "Indisponible".into();
                    s.frame.building = false;
                    drop(s);
                    repaint();
                }
            }
            built = Some(st.clone());
            hist.clear();
            exch.clear();
            split.clear();
            moment.clear();
            boxes.clear();
            span.clear();
            halted = false;
            t_change = 0.0;
            next_surf = 0.0;
            last_render = None;
            run_wall = 0.0;
            announced = false;
            steps = 8;
            t_rate = 0.0;
            perf = (Instant::now(), 0, perf.2, 0.0);
            continue;
        }
        let (Some(g), Some(b)) = (sim.as_mut(), body.as_ref()) else {
            std::thread::sleep(Duration::from_millis(100));
            continue;
        };
        let spt = st.l as f64 / U as f64;
        if built.as_ref().is_some_and(|x| x.xtr != st.xtr) {
            g.set_transition(st.xtr);
            t_change = g.time as f64 / spt;
            if let Some(x) = built.as_mut() {
                x.xtr = st.xtr;
            }
        }
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
        let t_before = g.time as f64 / spt;
        let finished = st.t_end > 0.0 && t_before >= st.t_end;
        let advance = st.running && !halted && !finished;
        let t0 = Instant::now();
        if advance {
            // Pas plus loin que la fin demandée.
            let left = if st.t_end > 0.0 {
                ((st.t_end * spt).ceil() as u64).saturating_sub(g.time).max(1) as usize
            } else {
                usize::MAX
            };
            g.run(steps.min(left), &ramp);
        } else {
            // Une passe sans avancer n'existe pas : en pause, on ne fait que l'image.
            std::thread::sleep(Duration::from_millis(40));
        }
        let f = g.take_forces();
        let dt = t0.elapsed().as_secs_f64();
        if advance {
            run_wall += dt;
            if dt < 0.04 {
                steps = (steps as f64 * 1.25).ceil().min(4000.0) as usize;
            } else if dt > 0.07 {
                steps = ((steps as f64 * 0.75).floor() as usize).max(1);
            }
            perf.1 += f.steps;
        }
        let q = 0.5 * (U as f64).powi(2) * b.area;
        let t_star = g.time as f64 / spt;
        let diverged = halted || (f.steps > 0 && !g.is_finite(&f));
        halted = diverged;
        if t_star + 1.0 < next_surf {
            next_surf = t_star;
        }
        let settle = if t_change == 0.0 { RAMP_T + SETTLE_T } else { SETTLE_T };
        let start = (t_change + settle).max(t_star - WINDOW_T);
        if b.mesh.is_some() {
            if f.steps > 0 && !diverged {
                exch.push([t_star, f.fx / q, f.fz / q, f.fy / q]);
            }
            if f.steps > 0 && !diverged && t_star >= next_surf {
                next_surf = t_star + 0.5;
                if let Some(s) = g.surface_forces() {
                    let t = s.total();
                    hist.push([t_star, t[0] / q, t[2] / q, t[1] / q]);
                    split.push([
                        t_star,
                        s.pressure_wall[0] / q,
                        s.friction[0] / q,
                        s.pressure_wall[2] / q,
                        s.friction[2] / q,
                    ]);
                    // Moment de tangage autour de (x0 + xref·L, centre en z),
                    // positif à cabrer : Cm = a + xref·b.
                    let bb = b.bounds;
                    let lf = st.l as f64;
                    let zr = 0.5 * (bb.z0 + bb.z1);
                    let (mut a, mut fz) = (0.0, 0.0);
                    for (p, f) in &s.loads {
                        a += (p[2] - zr) * f[0] - (p[0] - bb.x0) * f[2];
                        fz += f[2];
                    }
                    // Effort de chaque triangle réparti sur son étendue en
                    // envergure (densité d'aire le long de y), comme
                    // `run3d --envergure`.
                    let dy = (bb.y1 - bb.y0) / SPAN_BINS as f64;
                    let forces: Vec<[f64; 3]> = s.loads.iter().map(|l| l.1).collect();
                    let bins = lbm::gpu3d::span_distribution(&s.spans, &forces, bb.y0, dy, SPAN_BINS);
                    moment.push([t_star, a / (q * lf), fz / q]);
                    // c·cl/c_ref = (dL/dy) / (q∞ c_ref), c_ref = L mailles.
                    let qc = 0.5 * (U as f64).powi(2) * lf;
                    span.push((t_star, bins.iter().map(|v| v[2] / dy / qc).collect()));
                    if span.len() > 200 {
                        span.drain(..50);
                    }
                }
                // Bilan sur boîtes pendant la fenêtre de moyenne.
                if t_star >= t_change + settle {
                    if let [Some((f1, v1)), Some((f2, v2))] = g.box_forces_vol(&[0.1, 0.15])[..] {
                        boxes.push([t_star, f1[0] / q, f1[2] / q, v1, f2[0] / q, f2[2] / q, v2]);
                    }
                }
            }
        } else if f.steps > 0 && !diverged {
            hist.push([t_star, f.fx / q, f.fz / q, f.fy / q]);
        }
        for v in [&mut hist, &mut exch] {
            if v.len() > 20000 {
                v.drain(..5000);
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
        let key = ViewKey::of(&st);
        let fresh = match &last_render {
            Some((t, k)) => advance && *k == key && t.elapsed() < Duration::from_millis(120),
            None => false,
        };
        let (sl, image, camera) = if fresh {
            (None, None, None)
        } else if st.view3d {
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
                body_view: st.body_view,
            });
            (None, Some(Arc::new((st.img.0, st.img.1, rgba))), Some(cam))
        } else {
            (Some(Arc::new(g.slice(ax, cell, st.field, st.l as f64))), None, None)
        };
        if !fresh {
            last_render = Some((Instant::now(), key));
            slice_id += 1;
        }
        // Analyse sur la fenêtre établie.
        let w: Vec<&[f64; 4]> = hist.iter().filter(|p| p[0] >= start).collect();
        let tt: Vec<f64> = w.iter().map(|p| p[0]).collect();
        let a = analyze(
            &tt,
            &w.iter().map(|p| p[1]).collect::<Vec<_>>(),
            &w.iter().map(|p| p[2]).collect::<Vec<_>>(),
            &w.iter().map(|p| p[3]).collect::<Vec<_>>(),
        );
        let ready = t_star - start >= 3.0 && a.n >= 8 && a.cx.is_finite();
        let finished = st.t_end > 0.0 && t_star >= st.t_end - 1e-9;
        if finished != announced {
            announced = finished;
            if finished {
                eprintln!("Soufflerie 3D : {t_star:.1} t* calculés en {run_wall:.1} s ({:.0} M mailles/s)", perf.2);
            }
        }
        let status = if g.lost {
            "GPU perdu"
        } else if diverged {
            "Divergence"
        } else if finished {
            "Terminé"
        } else if !st.running {
            if t_star == 0.0 {
                "Prêt"
            } else {
                "En pause"
            }
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
            // Toutes les grilles : n_k mailles, 2^niveau pas par pas de la veine.
            let per_step: f64 = g.grid_cells().iter().map(|&(n, d)| n as f64 * (1u64 << d) as f64).sum();
            if perf.1 > 0 {
                perf.2 = perf.1 as f64 * per_step / el / 1e6;
                let rate = perf.1 as f64 / spt / el;
                t_rate = if t_rate > 0.0 { 0.5 * (t_rate + rate) } else { rate };
            } else if !advance {
                t_rate = 0.0;
            }
            perf = (Instant::now(), 0, perf.2, 0.0);
        }
        {
            let mut s = shared.lock().unwrap();
            let fr = &mut s.frame;
            if !fresh {
                fr.slice = sl;
                fr.image = image;
                fr.camera = camera;
                fr.slice_id = slice_id;
                fr.plane_pos = (ax, cell);
            }
            fr.run_wall = run_wall;
            fr.dims = [g.nx, g.ny, g.nz];
            fr.center = b.center;
            fr.l = st.l;
            fr.t_star = t_star;
            fr.hist = hist.clone();
            fr.exch = exch.clone();
            fr.split = split.clone();
            fr.moment = moment.clone();
            fr.boxes = boxes.clone();
            fr.span = span.clone();
            fr.span_eta = span_eta.clone();
            fr.window = start;
            fr.analysis = a;
            fr.ready = ready;
            fr.status = status.into();
            fr.advancing = advance && !diverged;
            fr.finished = finished;
            fr.building = false;
            fr.t_rate = t_rate;
            fr.mlups = perf.2;
            fr.info = format!(
                "τ = {:.4} · sous-maille {}{}",
                g.tau0,
                if g.les { "active" } else { "inactive" },
                if g.sa { " · Spalart–Allmaras + loi de paroi" } else { "" },
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

#[cfg(test)]
mod tests {
    use super::*;

    fn base() -> Settings {
        Settings {
            shape: Shape3D::Sphere,
            stl: None,
            setup: MeshSetup::default(),
            phys: false,
            speed: 20.0,
            fluid: Fluid::Air,
            re: 300.0,
            l: 16,
            alpha: 0.0,
            fp16: true,
            reset: 0,
            refine: 0,
            bands: 0,
            band_chord: true,
            wall: 0,
            xtr: 0.0,
            running: false,
            t_end: 40.0,
            field: SliceField::Vorticity,
            plane: Plane::Y,
            pos: None,
            view3d: true,
            az: 0.0,
            el: 0.0,
            dist: 6.0,
            iso: true,
            body_view: BodyView::Smooth,
            iso_level: 0.4,
            img: (800, 450),
        }
    }

    /// L'affichage et la marche ne comptent pas comme un changement de la
    /// physique ; la résolution, le Reynolds ou la maquette, si.
    #[test]
    fn reglages_physiques() {
        let a = base();
        let mut b = a.clone();
        b.running = true;
        b.az = 120.0;
        b.t_end = 10.0;
        b.iso = false;
        assert!(a.same_physics(&b));
        let mut c = a.clone();
        c.re = 1000.0;
        assert!(!a.same_physics(&c));
        let mut d = a.clone();
        d.setup.margin = 0.4;
        assert!(!a.same_physics(&d));
        let mut e = c.clone();
        e.take_view(&b);
        assert_eq!(e.az, 120.0);
        assert_eq!(e.re, 1000.0);
    }

    /// Couches pariétales comme `run3d --bande auto --bande-corde 0.125` :
    /// L = 16 et grille ×2, trois couches → 19, 23, 32 mailles (README).
    #[test]
    fn epaisseur_des_couches() {
        assert_eq!(band_thickness(3, 16, 1, true), vec![19.0, 23.0, 32.0]);
        assert_eq!(band_thickness(1, 16, 1, false), vec![16.0]);
        assert!(band_thickness(0, 16, 1, true).is_empty());
    }
}
