//! Fil de calcul de l'interface : fait avancer la veine, échantillonne les
//! efforts, colore le champ, transporte les traceurs, et publie une image
//! prête à afficher. L'interface ne fait que lire la dernière image.

use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use lbm::cpu::Macros;
use lbm::geometry::SOLID;
use lbm::panel::PanelResult;
use lbm::stats::{analyze, Analysis};
use lbm::wall::{wall_supported, WallProfile, WallSampler};
use lbm::{TunnelConfig, WindTunnel};
use rayon::prelude::*;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Field {
    Vorticity,
    Speed,
    Pressure,
    /// Viscosité turbulente de Spalart–Allmaras, log₁₀(1 + ν_t/ν).
    Eddy,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Tracers {
    Off,
    Streak,
    Fronts,
}

/// Réglages écrits par l'interface.
#[derive(Clone)]
pub struct Settings {
    pub cfg: TunnelConfig,
    pub field: Field,
    pub tracers: Tracers,
    pub running: bool,
    /// Coefficients corrigés du blocage de la veine.
    pub wall_corr: bool,
}

/// Demandes ponctuelles de l'interface.
#[derive(Default)]
pub struct Requests {
    /// Recréer la veine (résolution, moteur, schéma).
    pub rebuild: bool,
    /// Changer de maquette sans relancer l'écoulement.
    pub body: bool,
    /// Nouveau Reynolds.
    pub reynolds: bool,
    /// Modèle de turbulence ou transition modifiés.
    pub turbulence: bool,
}

/// Image publiée par le fil de calcul.
#[derive(Default)]
pub struct Frame {
    pub nx: usize,
    pub ny: usize,
    pub l: f64,
    /// Champ coloré (RGBA, ligne du haut en premier).
    pub rgba: Vec<u8>,
    pub image_id: u64,
    /// Blocage de la veine ε (0 sans correction).
    pub blockage: f64,
    /// Contour de la maquette et point de référence (réseau).
    pub outline: Vec<[f64; 2]>,
    pub refp: [f64; 2],
    /// Traceurs : points (filets) et polylignes (fronts), réseau.
    pub points: Vec<[f32; 2]>,
    pub lines: Vec<Vec<[f32; 2]>>,
    pub line_parity: Vec<bool>,
    /// Historique (t*, Cx, Cz, Cm).
    pub history: Vec<[f64; 4]>,
    pub analysis: Analysis,
    pub ready: bool,
    pub t_star: f64,
    pub ramping: bool,
    pub mlups: f64,
    pub engine: String,
    pub tau0: f64,
    pub les: bool,
    pub reliable_re: f64,
    pub label: String,
    pub diverged: bool,
    pub error: Option<String>,
    /// Spalart–Allmaras actif, y⁺ attendu au premier nœud.
    pub sa: bool,
    pub yplus: f64,
    /// Référence « panneaux + couche limite » (profils, Re ≥ 2·10⁴).
    pub panel: Option<PanelResult>,
    /// Distributions pariétales moyennées (après le transitoire).
    pub wall: Option<WallProfile>,
}

pub struct Shared {
    pub settings: Settings,
    pub req: Requests,
    pub frame: Arc<Frame>,
}

const WINDOW_T: f64 = 20.0;
const SETTLE_T: f64 = 6.0;
const FRONT_N: usize = 480;
const FRONT_SPAN: f64 = 4.0;
const FRONT_DX: f64 = 0.5;
const MAX_FRONTS: usize = 40;
const MAX_STREAK: usize = 12000;

/// Bornes et palette d'une grandeur affichée.
pub fn field_range(f: Field) -> (f32, f32, bool) {
    match f {
        Field::Vorticity => (-5.0, 5.0, true),
        Field::Speed => (0.0, 1.6, false),
        Field::Pressure => (-1.5, 1.5, true),
        Field::Eddy => (0.0, 4.0, false),
    }
}

/// Points, polylignes et parité (couleur) des polylignes.
type TracerShapes = (Vec<[f32; 2]>, Vec<Vec<[f32; 2]>>, Vec<bool>);

struct Particles {
    x: Vec<f32>,
    y: Vec<f32>,
    head: usize,
    emit: f64,
    slot: usize,
    mode: Tracers,
}

impl Particles {
    fn new() -> Particles {
        Particles {
            x: Vec::new(),
            y: Vec::new(),
            head: 0,
            emit: 0.0,
            slot: 0,
            mode: Tracers::Off,
        }
    }

    fn reset(&mut self, mode: Tracers) {
        self.mode = mode;
        self.head = 0;
        self.emit = 0.0;
        self.slot = 0;
        let n = match mode {
            Tracers::Fronts => FRONT_N * MAX_FRONTS,
            Tracers::Streak => MAX_STREAK,
            Tracers::Off => 0,
        };
        self.x = vec![-1.0; n];
        self.y = vec![0.0; n];
    }

    /// Émission puis transport sur `steps` pas (vitesse figée, sous-pas ≤ ½ maille).
    fn advect(&mut self, wt: &WindTunnel, m: &Macros, kind: &[u8], steps: u64) {
        if self.mode == Tracers::Off || steps == 0 || m.ux.is_empty() {
            return;
        }
        let l = wt.cfg.l as f64;
        let u = wt.cfg.u_in as f64;
        self.emit += steps as f64 * u;
        match self.mode {
            Tracers::Fronts => {
                let spacing = FRONT_DX * l;
                let x0 = (wt.x0 - 3.0 * l).max(2.0);
                while self.emit >= spacing {
                    self.emit -= spacing;
                    let base = self.slot * FRONT_N;
                    for k in 0..FRONT_N {
                        self.x[base + k] = (x0 + self.emit) as f32;
                        self.y[base + k] =
                            (wt.y0 - FRONT_SPAN * l / 2.0 + FRONT_SPAN * l * (k as f64 + 0.5) / FRONT_N as f64) as f32;
                    }
                    self.slot = (self.slot + 1) % MAX_FRONTS;
                }
            }
            Tracers::Streak => {
                let seeds = 17;
                let span = 5.0 * l;
                let spacing = (l / 6.4).max(2.5);
                while self.emit >= spacing {
                    self.emit -= spacing;
                    for s in 0..seeds {
                        let y = wt.y0 - span / 2.0 + span * (s as f64 + 0.5) / seeds as f64;
                        self.x[self.head] = (1.0 + self.emit) as f32;
                        self.y[self.head] = y as f32;
                        self.head = (self.head + 1) % MAX_STREAK;
                    }
                }
            }
            Tracers::Off => {}
        }
        let (nx, ny) = (m.nx, m.ny);
        let nsub = ((steps as f64 * u * 2.0 / 0.5).ceil() as usize).clamp(1, 2000);
        let dt = steps as f32 / nsub as f32;
        let (ux, uy) = (&m.ux, &m.uy);
        self.x.par_iter_mut().zip(self.y.par_iter_mut()).for_each(|(px, py)| {
            let (mut x, mut y) = (*px, *py);
            if x < 0.0 {
                return;
            }
            for _ in 0..nsub {
                let (i, j) = (x.floor() as isize, y.floor() as isize);
                if i < 0 || j < 0 || i as usize >= nx - 1 || j as usize >= ny - 1 {
                    x = -1.0;
                    break;
                }
                let (i, j) = (i as usize, j as usize);
                let (fx, fy) = (x - i as f32, y - j as f32);
                let a = j * nx + i;
                let (b, c, d) = (a + 1, a + nx, a + nx + 1);
                let vx = (ux[a] * (1.0 - fx) + ux[b] * fx) * (1.0 - fy) + (ux[c] * (1.0 - fx) + ux[d] * fx) * fy;
                let vy = (uy[a] * (1.0 - fx) + uy[b] * fx) * (1.0 - fy) + (uy[c] * (1.0 - fx) + uy[d] * fx) * fy;
                x += vx * dt;
                y += vy * dt;
                let (ci, cj) = (x.round() as isize, y.round() as isize);
                if ci >= 0
                    && cj >= 0
                    && (ci as usize) < nx
                    && (cj as usize) < ny
                    && kind[cj as usize * nx + ci as usize] == SOLID
                {
                    x = -1.0;
                    break;
                }
            }
            *px = x;
            *py = y;
        });
    }

    /// Points (filets) ou polylignes (fronts : segments rompus quand ils
    /// enjambent le corps ou s'étirent dans la couche limite).
    fn shapes(&self, nx: usize, ny: usize, kind: &[u8], l: f64) -> TracerShapes {
        let mut points = Vec::new();
        let mut lines = Vec::new();
        let mut parity = Vec::new();
        match self.mode {
            Tracers::Streak => {
                for (&x, &y) in self.x.iter().zip(&self.y) {
                    if x >= 0.0 {
                        points.push([x, y]);
                    }
                }
            }
            Tracers::Fronts => {
                let max_seg = (0.25 * l).max(2.0) as f32;
                let solid = |x: f32, y: f32| {
                    let (i, j) = (x.round() as isize, y.round() as isize);
                    i >= 0
                        && j >= 0
                        && (i as usize) < nx
                        && (j as usize) < ny
                        && kind[j as usize * nx + i as usize] == SOLID
                };
                let ok = |a: usize, b: usize| {
                    let (x0, y0, x1, y1) = (self.x[a], self.y[a], self.x[b], self.y[b]);
                    if x0 < 0.0 || x1 < 0.0 {
                        return false;
                    }
                    let (dx, dy) = (x1 - x0, y1 - y0);
                    if dx * dx + dy * dy > max_seg * max_seg {
                        return false;
                    }
                    (1..5).all(|t| {
                        let t = t as f32 * 0.2;
                        !solid(x0 + t * dx, y0 + t * dy)
                    })
                };
                for s in 0..MAX_FRONTS {
                    let base = s * FRONT_N;
                    let mut cur: Vec<[f32; 2]> = Vec::new();
                    for k in 0..FRONT_N - 1 {
                        let (a, b) = (base + k, base + k + 1);
                        if ok(a, b) {
                            if cur.is_empty() {
                                cur.push([self.x[a], self.y[a]]);
                            }
                            cur.push([self.x[b], self.y[b]]);
                        } else {
                            if cur.len() > 1 {
                                lines.push(std::mem::take(&mut cur));
                                parity.push(s % 2 == 1);
                            }
                            cur.clear();
                            if self.x[a] >= 0.0 && (k == 0 || !ok(a - 1, a)) {
                                points.push([self.x[a], self.y[a]]);
                            }
                        }
                    }
                    if cur.len() > 1 {
                        lines.push(cur);
                        parity.push(s % 2 == 1);
                    }
                }
            }
            Tracers::Off => {}
        }
        (points, lines, parity)
    }
}

/// Lance le fil de calcul ; `repaint` réveille l'interface.
pub fn spawn(shared: Arc<Mutex<Shared>>, repaint: impl Fn() + Send + 'static) {
    std::thread::spawn(move || run(shared, repaint));
}

fn run(shared: Arc<Mutex<Shared>>, repaint: impl Fn()) {
    let mut luts = (crate::theme::lut(crate::theme::diverging()), crate::theme::lut(crate::theme::sequential()));
    let mut dark = crate::theme::dark();
    let mut wt: Option<WindTunnel> = None;
    let mut err: Option<String> = None;
    let mut history: Vec<[f64; 4]> = Vec::new();
    let mut t_change = 0.0f64;
    let mut parts = Particles::new();
    let mut macros = Macros::default();
    let mut image_id = 0u64;
    let mut perf = (Instant::now(), 0u64, 0.0f64);
    let mut batch = 8usize;
    let mut last_frame = Instant::now() - Duration::from_secs(1);
    let mut steps_since_frame = 0u64;
    let mut diverged = false;
    let mut wall: Option<WallSampler> = None;
    let mut eddy: Vec<f32> = Vec::new();
    loop {
        let (settings, req) = {
            let mut s = shared.lock().unwrap();
            (s.settings.clone(), std::mem::take(&mut s.req))
        };
        if req.rebuild || (wt.is_none() && err.is_none()) {
            wt = None;
            match WindTunnel::new(settings.cfg.clone()) {
                Ok(w) => {
                    wt = Some(w);
                    err = None;
                }
                Err(e) => err = Some(e),
            }
            history.clear();
            t_change = 0.0;
            diverged = false;
            wall = None;
            parts.reset(settings.tracers);
            batch = 8;
        }
        if parts.mode != settings.tracers {
            parts.reset(settings.tracers);
        }
        let Some(w) = wt.as_mut() else {
            publish_error(&shared, err.clone());
            repaint();
            std::thread::sleep(Duration::from_millis(100));
            continue;
        };
        if req.body {
            w.set_body(settings.cfg.shape, &settings.cfg.params);
            t_change = w.t_star();
            parts.reset(settings.tracers);
            wall = None;
        }
        if req.reynolds {
            w.set_reynolds(settings.cfg.re);
            t_change = w.t_star();
        }
        if req.turbulence {
            w.set_turbulence(settings.cfg.turbulence);
            w.set_transition_mode(settings.cfg.transition);
            t_change = w.t_star();
        }
        if req.reynolds || req.turbulence {
            if let Some(ws) = wall.as_mut() {
                ws.reset();
            }
        }
        // Calcul : lots de pas, ≈ 16 échantillons d'efforts par temps convectif.
        let spt = w.steps_per_tstar();
        let max_batch = ((spt / 16.0).floor() as usize).max(1);
        if settings.running && !diverged {
            let t0 = Instant::now();
            while t0.elapsed() < Duration::from_millis(28) {
                let n = batch.min(max_batch);
                w.run(n);
                let spt = w.steps_per_tstar();
                for (step, s) in w.engine.poll_forces() {
                    let c = if settings.wall_corr {
                        w.corrected(w.coefficients(s))
                    } else {
                        w.coefficients(s)
                    };
                    history.push([step as f64 / spt, c.cx, c.cz, c.cm]);
                    if !c.cx.is_finite() || !c.cz.is_finite() {
                        diverged = true;
                    }
                }
                steps_since_frame += n as u64;
                perf.1 += n as u64;
                if diverged {
                    break;
                }
                batch = n;
            }
            batch = (batch * 2).min(max_batch);
            if history.len() > 60000 {
                history.drain(..20000);
            }
        } else {
            std::thread::sleep(Duration::from_millis(30));
        }
        // Image (≈ 30 par seconde).
        if last_frame.elapsed() >= Duration::from_millis(33) || !settings.running {
            last_frame = Instant::now();
            w.engine.macros(&mut macros);
            let has_eddy = settings.field == Field::Eddy && w.engine.eddy_viscosity_ratio(&mut eddy);
            if !w.engine.is_finite() {
                diverged = true;
            }
            let kind = match &w.engine {
                lbm::tunnel::Engine::Cpu(s) => &s.geo.kind,
                #[cfg(feature = "gpu")]
                lbm::tunnel::Engine::Gpu(s) => &s.geo.kind,
            };
            parts.advect(w, &macros, kind, steps_since_frame);
            steps_since_frame = 0;
            if crate::theme::dark() != dark {
                dark = !dark;
                luts = (crate::theme::lut(crate::theme::diverging()), crate::theme::lut(crate::theme::sequential()));
            }
            let rgba = colorize(
                &macros,
                kind,
                settings.field,
                w.cfg.l as f32,
                w.cfg.u_in,
                if has_eddy { &eddy } else { &[] },
                &luts,
            );
            let (points, lines, line_parity) = parts.shapes(w.nx, w.ny, kind, w.cfg.l as f64);
            let el = perf.0.elapsed().as_secs_f64();
            if el > 1.0 {
                perf.2 = (w.nx * w.ny) as f64 * perf.1 as f64 / el / 1e6;
                perf = (Instant::now(), 0, perf.2);
            }
            // Moyennes : après la mise en vitesse et l'établissement, sur 20 t* au plus.
            let now = w.t_star();
            let settle = if t_change == 0.0 {
                lbm::tunnel::RAMP_T + SETTLE_T
            } else {
                SETTLE_T
            };
            // Distributions pariétales : moyenne depuis la fin du transitoire.
            if wall_supported(w.cfg.shape) && !diverged {
                let ws = wall.get_or_insert_with(|| WallSampler::new(w));
                if now < t_change + settle {
                    ws.reset();
                } else if settings.running && !ws.is_empty() {
                    ws.add(w, &macros);
                }
            } else {
                wall = None;
            }
            let start = (t_change + settle).max(now - WINDOW_T);
            let win: Vec<&[f64; 4]> = history.iter().filter(|h| h[0] >= start).collect();
            let col = |k: usize| win.iter().map(|h| h[k]).collect::<Vec<_>>();
            let mut analysis = analyze(&col(0), &col(1), &col(2), &col(3));
            // Blocage : vitesse effective U(1 + ε), Strouhal ramené à cette vitesse.
            let blockage = if settings.wall_corr && analysis.cx.is_finite() {
                let (es, ew) = w.blockage(analysis.cx);
                es + ew
            } else {
                0.0
            };
            analysis.st /= 1.0 + blockage;
            let ready = now - start >= 3.0 && analysis.n >= 8 && analysis.cx.is_finite();
            image_id += 1;
            let frame = Frame {
                nx: w.nx,
                ny: w.ny,
                l: w.cfg.l as f64,
                rgba,
                image_id,
                blockage,
                outline: w.body.outline.clone(),
                refp: [w.body.ref_x, w.body.ref_y],
                points,
                lines,
                line_parity,
                history: history.iter().rev().take(4000).rev().cloned().collect(),
                analysis,
                ready,
                t_star: now,
                ramping: w.ramping(),
                mlups: perf.2,
                engine: w.engine.name(),
                tau0: w.engine.tau0(),
                les: w.engine.les_active(),
                reliable_re: w.reliable_reynolds(),
                label: w.body.label(),
                diverged,
                error: None,
                sa: w.sa_active(),
                yplus: w.first_cell_yplus(),
                panel: w.panel.clone(),
                wall: wall.as_ref().map(|ws| ws.result()),
            };
            shared.lock().unwrap().frame = Arc::new(frame);
            repaint();
        }
    }
}

fn publish_error(shared: &Arc<Mutex<Shared>>, e: Option<String>) {
    let mut s = shared.lock().unwrap();
    if s.frame.error != e {
        s.frame = Arc::new(Frame {
            error: e,
            ..Default::default()
        });
    }
}

/// Champ coloré, ligne du haut en premier (y vers le haut à l'écran).
fn colorize(
    m: &Macros,
    kind: &[u8],
    field: Field,
    l: f32,
    u: f32,
    eddy: &[f32],
    luts: &(Vec<[u8; 4]>, Vec<[u8; 4]>),
) -> Vec<u8> {
    let (nx, ny) = (m.nx, m.ny);
    let mut out = vec![0u8; 4 * nx * ny];
    let (lo, hi, div) = field_range(field);
    let lut = if div { &luts.0 } else { &luts.1 };
    let body = crate::theme::body_rgb();
    let k = 255.0 / (hi - lo);
    out.par_chunks_mut(4 * nx).enumerate().for_each(|(row, px)| {
        let y = ny - 1 - row;
        for x in 0..nx {
            let idx = y * nx + x;
            let c = if kind[idx] == SOLID {
                [body[0], body[1], body[2], 255]
            } else {
                let v = match field {
                    Field::Speed => (m.ux[idx] * m.ux[idx] + m.uy[idx] * m.uy[idx]).sqrt() / u,
                    Field::Pressure => (m.rho[idx] - 1.0) / 3.0 / (0.5 * u * u),
                    Field::Eddy => eddy.get(idx).map(|&r| (1.0 + r).log10()).unwrap_or(0.0),
                    Field::Vorticity => {
                        let xl = if x > 0 { idx - 1 } else { idx };
                        let xr = if x + 1 < nx { idx + 1 } else { idx };
                        let yd = if y > 0 { idx - nx } else { idx };
                        let yu = if y + 1 < ny { idx + nx } else { idx };
                        let dvdx = (m.uy[xr] - m.uy[xl]) / ((xr - xl).max(1) as f32);
                        let dudy = (m.ux[yu] - m.ux[yd]) / (((yu - yd) / nx).max(1) as f32);
                        (dvdx - dudy) * l / u
                    }
                };
                lut[((v - lo) * k).clamp(0.0, 255.0) as usize]
            };
            px[4 * x..4 * x + 4].copy_from_slice(&c);
        }
    });
    out
}
