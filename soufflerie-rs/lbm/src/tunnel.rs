//! Veine d'essai : domaine, maquette, mise en vitesse, coefficients.

use crate::cpu::{Collision, CpuSolver, ForceSum, LesMode, Macros, Ramp, SolverConfig};
use crate::geometry::Geometry;
#[cfg(feature = "gpu")]
use crate::gpu::GpuSolver;
#[cfg(feature = "gpu")]
use crate::refine::RefinedGpu;
use crate::panel::{panel_solve, PanelOptions, PanelResult, PanelShape};
use crate::shapes::{Body, ShapeKind, ShapeParams};

/// Domaine en multiples de L ; maquette à 8 L de l'entrée, à mi-hauteur.
pub const DOMAIN: (usize, usize, usize) = (24, 12, 8);
/// Durée de la mise en vitesse (temps convectifs).
pub const RAMP_T: f64 = 5.0;

/// Reynolds à partir duquel Spalart–Allmaras s'active en mode automatique.
pub const RE_SA_AUTO: f64 = 5e4;
/// Reynolds minimal de la référence « panneaux + couche limite ».
pub const RE_PANEL_MIN: f64 = 2e4;

/// Modèle de turbulence RANS.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Turbulence {
    /// Spalart–Allmaras à partir de Re = 5·10⁴.
    Auto,
    /// Spalart–Allmaras et loi de paroi.
    Sa,
    Off,
}

impl Turbulence {
    pub fn resolve(self, re: f64) -> bool {
        match self {
            Turbulence::Auto => re >= RE_SA_AUTO,
            Turbulence::Sa => true,
            Turbulence::Off => false,
        }
    }
}

/// Transition laminaire–turbulent de la couche limite (avec Spalart–Allmaras).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Transition {
    /// Libre : position donnée par la couche limite intégrale (eᴺ), zone
    /// laminaire imposée en amont.
    Free,
    /// Couche limite turbulente dès le bord d'attaque.
    Turbulent,
}

/// Longueur d'amortissement pariétal de la sous-maille (mailles), bornée à [1, 4].
pub fn les_wall_length(l: usize, re: f64) -> f64 {
    (8.0 * l as f64 / re.sqrt()).clamp(1.0, 4.0)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EngineKind {
    /// Processeur ; `threads` = 0 : tous les cœurs.
    Cpu {
        threads: usize,
    },
    Gpu,
}

#[derive(Clone, Debug)]
pub struct TunnelConfig {
    /// Mailles par longueur de référence.
    pub l: usize,
    pub re: f64,
    pub u_in: f32,
    pub shape: ShapeKind,
    pub params: ShapeParams,
    pub les: LesMode,
    pub collision: Collision,
    pub bouzidi: bool,
    pub engine: EngineKind,
    pub turbulence: Turbulence,
    pub transition: Transition,
    /// Moteur GPU : populations en demi-précision.
    pub gpu_half: bool,
    /// Moteur GPU : nombre de grilles imbriquées autour de la maquette, chacune
    /// deux fois plus fine que la précédente (0 : grille uniforme).
    pub refine: usize,
}

impl Default for TunnelConfig {
    fn default() -> Self {
        TunnelConfig {
            l: 16,
            re: 100.0,
            u_in: 0.1,
            shape: ShapeKind::Cylinder,
            params: ShapeParams::default(),
            les: LesMode::Auto,
            collision: Collision::Rr,
            bouzidi: true,
            engine: EngineKind::Cpu { threads: 0 },
            turbulence: Turbulence::Auto,
            transition: Transition::Free,
            gpu_half: false,
            refine: 0,
        }
    }
}

/// Coefficients aérodynamiques (Cm positif cabreur).
#[derive(Clone, Copy, Debug, Default)]
pub struct Coefficients {
    pub cx: f64,
    pub cz: f64,
    pub cm: f64,
}

pub enum Engine {
    Cpu(Box<CpuSolver>),
    #[cfg(feature = "gpu")]
    Gpu(Box<GpuSolver>),
    #[cfg(feature = "gpu")]
    Refined(Box<RefinedGpu>),
}

macro_rules! dispatch {
    ($self:expr, $s:ident => $e:expr) => {
        match $self {
            Engine::Cpu($s) => $e,
            #[cfg(feature = "gpu")]
            Engine::Gpu($s) => $e,
            #[cfg(feature = "gpu")]
            Engine::Refined($s) => $e,
        }
    };
}

impl Engine {
    pub fn run(&mut self, n: usize, ramp: &Ramp) {
        dispatch!(self, s => s.run(n, ramp))
    }
    pub fn take_forces(&mut self) -> ForceSum {
        dispatch!(self, s => s.take_forces())
    }
    /// Efforts relevés sans bloquer le GPU : (pas atteint, somme) des relevés
    /// disponibles, éventuellement avec un ou deux relevés de retard.
    pub fn poll_forces(&mut self) -> Vec<(u64, ForceSum)> {
        match self {
            Engine::Cpu(s) => {
                let t = s.time;
                let f = s.take_forces();
                if f.steps > 0 {
                    vec![(t, f)]
                } else {
                    Vec::new()
                }
            }
            #[cfg(feature = "gpu")]
            Engine::Gpu(s) => s.poll_forces(),
            #[cfg(feature = "gpu")]
            Engine::Refined(s) => s.poll_forces(),
        }
    }
    pub fn macros(&mut self, out: &mut Macros) {
        dispatch!(self, s => s.macros(out))
    }
    pub fn time(&self) -> u64 {
        match self {
            Engine::Cpu(s) => s.time,
            #[cfg(feature = "gpu")]
            Engine::Gpu(s) => s.time,
            #[cfg(feature = "gpu")]
            Engine::Refined(s) => s.time(),
        }
    }
    pub fn tau0(&self) -> f64 {
        dispatch!(self, s => s.tau0())
    }
    pub fn les_active(&self) -> bool {
        dispatch!(self, s => s.les_active())
    }
    pub fn is_finite(&mut self) -> bool {
        dispatch!(self, s => s.is_finite())
    }
    pub fn set_viscosity(&mut self, nu: f64) {
        dispatch!(self, s => s.set_viscosity(nu))
    }
    /// Amortissement pariétal de la sous-maille pour L mailles par longueur
    /// de référence (chaque grille fine avec sa propre résolution).
    pub fn set_les_wall_length(&mut self, l: usize, re: f64) {
        match self {
            #[cfg(feature = "gpu")]
            Engine::Refined(s) => s.set_les_wall_lengths(|k| les_wall_length(l << k, re)),
            Engine::Cpu(s) => s.set_les_wall_length(les_wall_length(l, re)),
            #[cfg(feature = "gpu")]
            Engine::Gpu(s) => s.set_les_wall_length(les_wall_length(l, re)),
        }
    }
    pub fn set_turbulence(&mut self, on: bool) {
        dispatch!(self, s => s.set_turbulence(on))
    }
    pub fn sa_active(&self) -> bool {
        dispatch!(self, s => s.sa_active())
    }
    pub fn geometry(&self) -> &Geometry {
        match self {
            Engine::Cpu(s) => &s.geo,
            #[cfg(feature = "gpu")]
            Engine::Gpu(s) => &s.geo,
            #[cfg(feature = "gpu")]
            Engine::Refined(s) => s.geometry(),
        }
    }
    pub fn geometry_mut(&mut self) -> &mut Geometry {
        match self {
            Engine::Cpu(s) => &mut s.geo,
            #[cfg(feature = "gpu")]
            Engine::Gpu(s) => &mut s.geo,
            #[cfg(feature = "gpu")]
            Engine::Refined(s) => s.geometry_mut(),
        }
    }
    pub fn update_laminar_zone(&mut self) {
        dispatch!(self, s => s.update_laminar_zone())
    }
    /// Rapport ν_t/ν par maille (Spalart–Allmaras) ; false si le modèle est inactif.
    pub fn eddy_viscosity_ratio(&mut self, out: &mut Vec<f32>) -> bool {
        match self {
            Engine::Cpu(s) => {
                let inv_nu = (1.0 / s.cfg.nu) as f32;
                let Some(st) = s.sa_fields() else { return false };
                out.clear();
                out.extend(st.nu_tilde.iter().map(|&v| crate::sa::nu_turb(v, inv_nu) * inv_nu));
                true
            }
            #[cfg(feature = "gpu")]
            Engine::Gpu(_) | Engine::Refined(_) => {
                let inv_nu = (1.0 / self.nu()) as f32;
                let st = match self {
                    Engine::Gpu(s) => s.sa_state(),
                    Engine::Refined(s) => s.sa_state(),
                    _ => Vec::new(),
                };
                if st.is_empty() {
                    return false;
                }
                out.clear();
                out.extend(st.iter().map(|c| crate::sa::nu_turb(c[3], inv_nu) * inv_nu));
                true
            }
        }
    }
    /// Géométrie de la grille la plus fine (la veine sans raffinement).
    pub fn wall_geometry(&self) -> &Geometry {
        match self {
            #[cfg(feature = "gpu")]
            Engine::Refined(s) => &s.levels.last().unwrap().geo,
            e => e.geometry(),
        }
    }
    /// ρ, u de la grille la plus fine.
    pub fn wall_macros(&mut self, out: &mut Macros) {
        match self {
            #[cfg(feature = "gpu")]
            Engine::Refined(s) => {
                let k = s.levels.len() - 1;
                s.level_macros(k, out)
            }
            e => e.macros(out),
        }
    }
    /// Viscosité (réseau de la veine).
    pub fn nu(&self) -> f64 {
        match self {
            Engine::Cpu(s) => s.cfg.nu,
            #[cfg(feature = "gpu")]
            Engine::Gpu(s) => s.cfg.nu,
            #[cfg(feature = "gpu")]
            Engine::Refined(s) => s.levels[0].cfg.nu,
        }
    }
    pub fn name(&self) -> String {
        match self {
            Engine::Cpu(s) => format!("CPU ({} fils)", s.threads()),
            #[cfg(feature = "gpu")]
            Engine::Gpu(s) => format!("GPU ({})", s.adapter_name()),
            #[cfg(feature = "gpu")]
            Engine::Refined(s) => format!("GPU ({})", s.adapter_name()),
        }
    }
}

pub struct WindTunnel {
    pub cfg: TunnelConfig,
    pub nx: usize,
    pub ny: usize,
    pub x0: f64,
    pub y0: f64,
    pub body: Body,
    pub engine: Engine,
    pub ramp: Ramp,
    /// Référence « panneaux + couche limite » à l'incidence courante (profils, Re ≥ 2·10⁴).
    pub panel: Option<PanelResult>,
}

impl WindTunnel {
    pub fn new(cfg: TunnelConfig) -> Result<WindTunnel, String> {
        let l = cfg.l;
        let (nx, ny) = (DOMAIN.0 * l, DOMAIN.1 * l);
        let x0 = (DOMAIN.2 * l) as f64;
        let y0 = (ny as f64 - 1.0) / 2.0;
        let scfg = SolverConfig {
            nx,
            ny,
            u_in: cfg.u_in,
            nu: cfg.u_in as f64 * l as f64 / cfg.re,
            les: cfg.les,
            cs: 0.15,
            les_wall_a: les_wall_length(l, cfg.re),
            collision: cfg.collision,
            bouzidi: cfg.bouzidi,
            sponge_fraction: 0.12,
            half: cfg.gpu_half,
        };
        let body = Body::new(cfg.shape, &cfg.params, l as f64, x0, y0, cfg.u_in as f64);
        let mut engine = match cfg.engine {
            EngineKind::Cpu { threads } => Engine::Cpu(Box::new(CpuSolver::new(scfg, threads))),
            #[cfg(feature = "gpu")]
            EngineKind::Gpu if cfg.refine > 0 => {
                Engine::Refined(Box::new(RefinedGpu::new(scfg, cfg.refine, l, x0, y0)?))
            }
            #[cfg(feature = "gpu")]
            EngineKind::Gpu => Engine::Gpu(Box::new(GpuSolver::new(scfg)?)),
            #[cfg(not(feature = "gpu"))]
            EngineKind::Gpu => return Err("moteur GPU non compilé (option « gpu »)".into()),
        };
        let sa = cfg.turbulence.resolve(cfg.re);
        dispatch!(&mut engine, s => {
            s.set_body(&body);
            s.set_turbulence(sa);
            // Démarrage progressif depuis le repos.
            s.reset(0.0);
        });
        engine.set_les_wall_length(l, cfg.re);
        let ramp = Ramp {
            steps: (RAMP_T * l as f64 / cfg.u_in as f64).round() as u64,
            u: cfg.u_in,
        };
        let mut wt = WindTunnel {
            cfg,
            nx,
            ny,
            x0,
            y0,
            body,
            engine,
            ramp,
            panel: None,
        };
        wt.refresh_panel();
        Ok(wt)
    }

    pub fn steps_per_tstar(&self) -> f64 {
        self.cfg.l as f64 / self.cfg.u_in as f64
    }

    pub fn t_star(&self) -> f64 {
        self.engine.time() as f64 / self.steps_per_tstar()
    }

    pub fn ramping(&self) -> bool {
        self.engine.time() < self.ramp.steps
    }

    pub fn run(&mut self, n: usize) {
        self.engine.run(n, &self.ramp);
    }

    /// Change de maquette (ou d'incidence) sans relancer l'écoulement.
    pub fn set_body(&mut self, shape: ShapeKind, params: &ShapeParams) {
        self.cfg.shape = shape;
        self.cfg.params = params.clone();
        self.body = Body::new(shape, params, self.cfg.l as f64, self.x0, self.y0, self.cfg.u_in as f64);
        let body = &self.body;
        dispatch!(&mut self.engine, s => s.set_body(body));
        self.refresh_panel();
    }

    /// Reynolds au-delà duquel la traînée est surestimée de plus de ~10 %
    /// faute de résolution (mesuré sur le NACA 0012, voir l'application web).
    pub fn reliable_reynolds(&self) -> f64 {
        match self.cfg.l {
            12 => 350.0,
            16 => 650.0,
            24 => 1450.0,
            32 => 5000.0,
            48 | 64 => 10000.0,
            l => ((2.5 * (l * l) as f64 / 50.0).round() * 50.0).max(350.0),
        }
    }

    /// Change le Reynolds sans relancer le calcul.
    pub fn set_reynolds(&mut self, re: f64) {
        self.cfg.re = re;
        self.engine.set_viscosity(self.cfg.u_in as f64 * self.cfg.l as f64 / re);
        self.engine.set_les_wall_length(self.cfg.l, re);
        self.engine.set_turbulence(self.cfg.turbulence.resolve(re));
        self.refresh_panel();
    }

    /// Modèle de turbulence, sans relancer l'écoulement.
    pub fn set_turbulence(&mut self, t: Turbulence) {
        self.cfg.turbulence = t;
        self.engine.set_turbulence(t.resolve(self.cfg.re));
        self.refresh_panel();
    }

    /// Transition libre ou imposée au bord d'attaque.
    pub fn set_transition_mode(&mut self, t: Transition) {
        self.cfg.transition = t;
        self.refresh_panel();
    }

    pub fn sa_active(&self) -> bool {
        self.engine.sa_active()
    }

    /// Méthode des panneaux et couche limite intégrale pour la maquette
    /// courante à l'incidence `alpha_deg` (None : forme non profilée ou Re trop bas).
    pub fn panel_at(&self, alpha_deg: f64) -> Option<PanelResult> {
        let shape = PanelShape::from_shape(self.cfg.shape, &self.cfg.params)?;
        let re = self.cfg.re;
        if re < RE_PANEL_MIN {
            return None;
        }
        let opts = match self.cfg.transition {
            Transition::Free => PanelOptions::free(re),
            Transition::Turbulent => PanelOptions::tripped(re, 0.0),
        };
        let r = panel_solve(&shape, alpha_deg, &opts);
        (r.cl.is_finite() && r.cd.is_finite()).then_some(r)
    }

    /// Recalcule la référence et la zone laminaire imposée au calcul LBM.
    pub fn refresh_panel(&mut self) {
        self.panel = self.panel_at(self.body.alpha_deg);
        let couple = self.sa_active() && self.cfg.transition == Transition::Free;
        let bl = if couple { self.panel.clone() } else { None };
        self.set_transition(bl.as_ref());
    }

    /// Transition imposée à Spalart–Allmaras à partir de la couche limite
    /// intégrale : en amont de la transition, pas de production turbulente
    /// et frottement pariétal laminaire (Cf de la couche limite), dans une
    /// bande de max(3, L/10) mailles autour du corps. None : turbulent partout.
    pub fn set_transition(&mut self, bl: Option<&PanelResult>) {
        let l = self.cfg.l as f64;
        let by_index = self.cfg.shape == ShapeKind::Naca;
        let loc = &self.body.outline_local;
        match &mut self.engine {
            #[cfg(feature = "gpu")]
            Engine::Refined(r) => {
                for (k, lev) in r.levels.iter_mut().enumerate() {
                    laminar_zone(&mut lev.geo, loc, l * (1u64 << k) as f64, by_index, bl);
                }
            }
            e => laminar_zone(e.geometry_mut(), loc, l, by_index, bl),
        }
        self.engine.update_laminar_zone();
    }

    /// Maquette dans la grille la plus fine.
    pub fn wall_body(&self) -> &Body {
        match &self.engine {
            #[cfg(feature = "gpu")]
            Engine::Refined(r) => r.bodies.last().unwrap_or(&self.body),
            _ => &self.body,
        }
    }

    /// Mailles par longueur de référence au plus près du corps.
    pub fn wall_l(&self) -> usize {
        match &self.engine {
            #[cfg(feature = "gpu")]
            Engine::Refined(r) => self.cfg.l << (r.levels.len() - 1),
            _ => self.cfg.l,
        }
    }

    /// y⁺ attendu au premier nœud (plaque plane turbulente,
    /// Cf ≈ 0,026 Re^(−1/7), nœud à une demi-maille de la paroi).
    pub fn first_cell_yplus(&self) -> f64 {
        let re = self.cfg.re;
        let cf = 0.026 * re.powf(-1.0 / 7.0);
        0.5 / self.wall_l() as f64 * re * (cf / 2.0).sqrt()
    }

    /// Coefficients moyens d'une somme d'efforts.
    pub fn coefficients(&self, s: ForceSum) -> Coefficients {
        let steps = s.steps.max(1) as f64;
        let u = self.cfg.u_in as f64;
        let q = 0.5 * u * u;
        let l = self.cfg.l as f64;
        Coefficients {
            cx: s.fx / steps / (q * l),
            cz: s.fy / steps / (q * l),
            cm: -s.mz / steps / (q * l * l),
        }
    }

    /// Blocage de la veine (Allen et Vincenti, 1944) : (ε_s, ε_w), survitesse
    /// relative au droit de la maquette due au volume de la maquette et à son
    /// sillage, pour une veine de hauteur ny et un Cx brut donné.
    ///
    /// ε_s = (π/6)(1 + e/c)·A/h² (A : section, e/c : épaisseur sur longueur ;
    /// π²/12·(d/h)² pour un cylindre) ; ε_w = L·Cx/(4h).
    pub fn blockage(&self, cx: f64) -> (f64, f64) {
        let ol = &self.body.outline;
        let h = self.ny as f64;
        if ol.len() < 3 {
            return (0.0, 0.0);
        }
        let mut area = 0.0;
        let (mut x0, mut x1, mut y0, mut y1) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
        for (k, p) in ol.iter().enumerate() {
            let q = ol[(k + 1) % ol.len()];
            area += p[0] * q[1] - q[0] * p[1];
            (x0, x1, y0, y1) = (x0.min(p[0]), x1.max(p[0]), y0.min(p[1]), y1.max(p[1]));
        }
        let area = 0.5 * area.abs();
        let ratio = ((y1 - y0) / (x1 - x0).max(1e-9)).min(4.0);
        let eps_s = std::f64::consts::PI / 6.0 * (1.0 + ratio) * area / (h * h);
        let eps_w = (self.cfg.l as f64 * cx / (4.0 * h)).max(0.0);
        (eps_s, eps_w)
    }

    /// Coefficients corrigés du blocage (Allen et Vincenti) : ramenés à la
    /// vitesse effective U(1 + ε) ; la traînée perd en outre la poussée
    /// horizontale due au gradient de pression du sillage (ε_s).
    pub fn corrected(&self, c: Coefficients) -> Coefficients {
        let (es, ew) = self.blockage(c.cx);
        let k = 1.0 - 2.0 * (es + ew);
        Coefficients {
            cx: c.cx * (1.0 - 3.0 * es - 2.0 * ew),
            cz: c.cz * k,
            cm: c.cm * k,
        }
    }

    /// Mach réseau.
    pub fn mach(&self) -> f64 {
        self.cfg.u_in as f64 * 3f64.sqrt()
    }
}

/// Zone laminaire d'une grille : en amont de la transition, Cf de la couche
/// limite intégrale dans une bande de max(3, L/10) mailles autour du corps.
fn laminar_zone(geo: &mut Geometry, loc: &[[f64; 2]], l: f64, by_index: bool, bl: Option<&PanelResult>) {
    geo.lam.iter_mut().for_each(|v| *v = 0.0);
    if let Some(bl) = bl {
        let m = loc.len();
        let (mut i_le, mut xmin, mut xmax) = (0usize, f64::INFINITY, f64::NEG_INFINITY);
        for (k, p) in loc.iter().enumerate() {
            if p[0] < xmin {
                xmin = p[0];
                i_le = k;
            }
            xmax = xmax.max(p[0]);
        }
        let span = if xmax > xmin { xmax - xmin } else { 1.0 };
        let band = (0.1 * l).max(3.0) as f32;
        let cf_at = |x: &[f64], cf: &[f64], xc: f64| -> f32 {
            if x.len() < 2 || cf.len() < x.len() {
                return 0.0;
            }
            let mut k = 1;
            while k < x.len() - 1 && x[k] < xc {
                k += 1;
            }
            // Plafond : pic de Cf au point d'arrêt (abscisse curviligne nulle).
            cf[k].clamp(1e-5, 0.02) as f32
        };
        for i in 0..geo.n {
            let pos = geo.near_pos[i];
            if pos < 0.0 || geo.dist[i] > band {
                continue;
            }
            let k = (pos.round() as usize).min(m - 1);
            let [bx, by] = loc[k];
            // Faces : pour le NACA, extrados de l'indice 0 (bord de fuite)
            // au bord d'attaque ; ailleurs, signe de y dans le repère propre.
            let upper = if by_index { k <= i_le } else { by >= 0.0 };
            let xc = (bx - xmin) / span;
            let xtr = if upper {
                bl.transition_upper
            } else {
                bl.transition_lower
            };
            if xtr.is_finite() && xc < xtr {
                let side = if upper { &bl.friction_upper } else { &bl.friction_lower };
                geo.lam[i] = cf_at(&side.x, &side.cf, xc);
            }
        }
    }
}
