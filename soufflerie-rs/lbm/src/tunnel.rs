//! Veine d'essai : domaine, maquette, mise en vitesse, coefficients.

use crate::cpu::{Collision, CpuSolver, ForceSum, LesMode, Macros, Ramp, SolverConfig};
#[cfg(feature = "gpu")]
use crate::gpu::GpuSolver;
use crate::shapes::{Body, ShapeKind, ShapeParams};

/// Domaine en multiples de L ; maquette à 8 L de l'entrée, à mi-hauteur.
pub const DOMAIN: (usize, usize, usize) = (24, 12, 8);
/// Durée de la mise en vitesse (temps convectifs).
pub const RAMP_T: f64 = 5.0;

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
}

macro_rules! dispatch {
    ($self:expr, $s:ident => $e:expr) => {
        match $self {
            Engine::Cpu($s) => $e,
            #[cfg(feature = "gpu")]
            Engine::Gpu($s) => $e,
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
    pub fn macros(&mut self, out: &mut Macros) {
        dispatch!(self, s => s.macros(out))
    }
    pub fn time(&self) -> u64 {
        dispatch!(self, s => s.time)
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
    pub fn set_les_wall_length(&mut self, a: f64) {
        dispatch!(self, s => s.set_les_wall_length(a))
    }
    pub fn name(&self) -> String {
        match self {
            Engine::Cpu(s) => format!("CPU ({} fils)", s.threads()),
            #[cfg(feature = "gpu")]
            Engine::Gpu(s) => format!("GPU ({})", s.adapter_name()),
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
        };
        let body = Body::new(cfg.shape, &cfg.params, l as f64, x0, y0, cfg.u_in as f64);
        let mut engine = match cfg.engine {
            EngineKind::Cpu { threads } => Engine::Cpu(Box::new(CpuSolver::new(scfg, threads))),
            #[cfg(feature = "gpu")]
            EngineKind::Gpu => Engine::Gpu(Box::new(GpuSolver::new(scfg)?)),
            #[cfg(not(feature = "gpu"))]
            EngineKind::Gpu => return Err("moteur GPU non compilé (option « gpu »)".into()),
        };
        dispatch!(&mut engine, s => {
            s.set_body(&body);
            // Démarrage progressif depuis le repos.
            s.reset(0.0);
        });
        let ramp = Ramp {
            steps: (RAMP_T * l as f64 / cfg.u_in as f64).round() as u64,
            u: cfg.u_in,
        };
        Ok(WindTunnel {
            cfg,
            nx,
            ny,
            x0,
            y0,
            body,
            engine,
            ramp,
        })
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
        self.engine.set_les_wall_length(les_wall_length(self.cfg.l, re));
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

    /// Mach réseau.
    pub fn mach(&self) -> f64 {
        self.cfg.u_in as f64 * 3f64.sqrt()
    }
}
