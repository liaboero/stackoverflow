//! Moteur processeur : propagation « pull » et collision fusionnées, ligne par
//! ligne sur tous les cœurs (rayon). Les intervalles de mailles fluides
//! ordinaires passent par une boucle sans branchement, que le compilateur
//! vectorise ; les mailles voisines de la paroi (rebond de Bouzidi, parois
//! haute et basse) passent par le chemin lent.
//!
//! Disposition mémoire : `f[(y * 9 + i) * nx + x]` (neuf plans par ligne).

use rayon::prelude::*;

use crate::d2q9::{collide, equilibrium, moments, Relax, CX, CY, OPP, REFY, W};
use crate::geometry::{Geometry, NEAR, SOLID};
use crate::sa::{self, SaFields};
use crate::shapes::Body;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Collision {
    /// Régularisée récursive (Malaspinas 2015), par défaut.
    Rr,
    Bgk,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LesMode {
    Auto,
    On,
    Off,
}

/// Paramètres communs aux moteurs.
#[derive(Clone, Debug)]
pub struct SolverConfig {
    pub nx: usize,
    pub ny: usize,
    pub u_in: f32,
    pub nu: f64,
    pub les: LesMode,
    pub cs: f32,
    pub les_wall_a: f64,
    pub collision: Collision,
    pub bouzidi: bool,
    pub sponge_fraction: f64,
    /// Couches absorbantes : `None`, bandes haute et basse de ny/8 lignes ;
    /// `Some(w)` (moteur GPU seulement) : couches de w mailles sur les quatre
    /// côtés, entrée et sortie comprises (calculs aéroacoustiques).
    pub absorb: Option<usize>,
    /// Moteur GPU : populations en demi-précision (deux fois moins de mémoire).
    pub half: bool,
}

/// Somme des efforts sur un nombre de pas.
#[derive(Clone, Copy, Debug, Default)]
pub struct ForceSum {
    pub fx: f64,
    pub fy: f64,
    pub mz: f64,
    pub steps: u64,
}

/// Champs macroscopiques (lignes de nx mailles).
#[derive(Clone, Default)]
pub struct Macros {
    pub nx: usize,
    pub ny: usize,
    pub rho: Vec<f32>,
    pub ux: Vec<f32>,
    pub uy: Vec<f32>,
}

/// Loi de mise en vitesse : (1 − cos) sur `steps` pas jusqu'à `u`.
#[derive(Clone, Copy, Debug)]
pub struct Ramp {
    pub steps: u64,
    pub u: f32,
}

impl Ramp {
    /// Vitesse d'entrée et accélération uniforme au pas `t`.
    pub fn at(&self, t: u64) -> (f32, f32) {
        if t >= self.steps || self.steps == 0 {
            return (self.u, 0.0);
        }
        let tt = self.steps as f64;
        let ph = std::f64::consts::PI * (t as f64 + 0.5) / tt;
        let u = 0.5 * self.u as f64 * (1.0 - ph.cos());
        let a = 0.5 * self.u as f64 * std::f64::consts::PI / tt * ph.sin();
        (u as f32, a as f32)
    }
}

pub struct CpuSolver {
    pub cfg: SolverConfig,
    pub geo: Geometry,
    f: Vec<f32>,
    g: Vec<f32>,
    tau_col: Vec<f32>,
    sigma_row: Vec<f32>,
    tau0: f64,
    les_active: bool,
    pub time: u64,
    u_in: f32,
    accel: f32,
    ref_x: f64,
    ref_y: f64,
    spin: f64,
    rot_x: f64,
    rot_y: f64,
    forces: ForceSum,
    last: (f64, f64, f64),
    pool: Option<rayon::ThreadPool>,
    /// Spalart–Allmaras actif ; champs alloués à la première activation.
    sa_active: bool,
    sa: Option<SaFields>,
}

impl CpuSolver {
    /// `threads` : nombre de fils (0 : tous les cœurs).
    pub fn new(cfg: SolverConfig, threads: usize) -> CpuSolver {
        let (nx, ny) = (cfg.nx, cfg.ny);
        let n = nx * ny;
        assert!(cfg.absorb.is_none(), "couches absorbantes sur quatre côtés : moteur GPU seulement");
        let absorb = ((ny as f64) / 8.0).round();
        let sigma_row = (0..ny)
            .map(|y| {
                let d = y.min(ny - 1 - y) as f64;
                let s = if absorb > 0.0 { 1.0 - d / absorb } else { 0.0 };
                if s > 0.0 {
                    (0.05 * s * s) as f32
                } else {
                    0.0
                }
            })
            .collect();
        let pool = if threads > 0 {
            rayon::ThreadPoolBuilder::new().num_threads(threads).build().ok()
        } else {
            None
        };
        let mut s = CpuSolver {
            geo: Geometry::new(nx, ny, cfg.les_wall_a),
            f: vec![0.0; 9 * n],
            g: vec![0.0; 9 * n],
            tau_col: vec![0.0; nx],
            sigma_row,
            tau0: 0.5,
            les_active: false,
            time: 0,
            u_in: cfg.u_in,
            accel: 0.0,
            ref_x: nx as f64 / 4.0,
            ref_y: (ny as f64 - 1.0) / 2.0,
            spin: 0.0,
            rot_x: nx as f64 / 4.0,
            rot_y: (ny as f64 - 1.0) / 2.0,
            forces: ForceSum::default(),
            last: (0.0, 0.0, 0.0),
            pool,
            sa_active: false,
            sa: None,
            cfg,
        };
        s.set_viscosity(s.cfg.nu);
        s.reset(s.cfg.u_in);
        s
    }

    pub fn tau0(&self) -> f64 {
        self.tau0
    }

    /// Nombre de fils de calcul.
    pub fn threads(&self) -> usize {
        self.pool
            .as_ref()
            .map(|p| p.current_num_threads())
            .unwrap_or_else(rayon::current_num_threads)
    }

    pub fn les_active(&self) -> bool {
        self.les_active
    }

    /// Viscosité moléculaire ; τ par colonne (zone éponge comprise).
    pub fn set_viscosity(&mut self, nu: f64) {
        self.cfg.nu = nu;
        self.tau0 = 3.0 * nu + 0.5;
        let nx = self.cfg.nx;
        let s0 = (nx as f64 * (1.0 - self.cfg.sponge_fraction)).round() as usize;
        let tau_sponge = self.tau0.max(0.9);
        let len = (nx - 1).saturating_sub(s0).max(1) as f64;
        for x in 0..nx {
            let mut tau = self.tau0;
            if x > s0 {
                let s = (x - s0) as f64 / len;
                tau += (tau_sponge - self.tau0) * s * s;
            }
            self.tau_col[x] = tau as f32;
        }
        self.les_active = match self.cfg.les {
            LesMode::On => true,
            LesMode::Off => false,
            LesMode::Auto => self.tau0 < 0.51,
        };
    }

    pub fn set_les_wall_length(&mut self, a: f64) {
        self.geo.set_les_wall_length(a);
    }

    /// Active ou coupe Spalart–Allmaras avec loi de paroi.
    pub fn set_turbulence(&mut self, on: bool) {
        let was = self.sa_active;
        self.sa_active = on;
        if on && self.sa.is_none() {
            self.sa = Some(SaFields::new(self.cfg.nx * self.cfg.ny));
        }
        if on && !was {
            self.init_sa();
        }
    }

    pub fn sa_active(&self) -> bool {
        self.sa_active
    }

    /// Champs de Spalart–Allmaras (ν̃, u_τ…), s'il a été activé.
    pub fn sa_fields(&self) -> Option<&SaFields> {
        self.sa.as_ref().filter(|_| self.sa_active)
    }

    /// Zone laminaire modifiée (`geo.lam`) : rien à faire sur processeur.
    pub fn update_laminar_zone(&mut self) {}

    fn init_sa(&mut self) {
        if let Some(st) = self.sa.as_mut() {
            st.init(&self.geo, self.cfg.nu);
        }
    }

    /// Définit le corps.
    pub fn set_body(&mut self, body: &Body) {
        self.ref_x = body.ref_x;
        self.ref_y = body.ref_y;
        self.spin = body.spin;
        self.rot_x = body.ref_x;
        self.rot_y = body.ref_y;
        let revived = self.geo.build(body);
        let nx = self.cfg.nx;
        let eq = equilibrium(1.0, 0.0, 0.0);
        for idx in revived {
            let (x, y) = (idx % nx, idx / nx);
            for i in 0..9 {
                self.f[(y * 9 + i) * nx + x] = eq[i];
            }
        }
    }

    /// Fluide à la vitesse u0, avec une petite perturbation transverse dans
    /// le sillage pour déclencher rapidement l'éventuel lâcher tourbillonnaire.
    pub fn reset(&mut self, u0: f32) {
        let (nx, ny) = (self.cfg.nx, self.cfg.ny);
        let kick = 0.15 * u0.max(self.cfg.u_in).max(0.05);
        let r = ny as f64 / 8.0;
        let px = self.ref_x + 1.5 * r;
        let py = self.ref_y + 0.5 * r;
        for y in 0..ny {
            for x in 0..nx {
                let idx = y * nx + x;
                let eq = if self.geo.kind[idx] == SOLID {
                    equilibrium(1.0, 0.0, 0.0)
                } else {
                    let d2 = (x as f64 - px).powi(2) + (y as f64 - py).powi(2);
                    let vy = kick as f64 * (-d2 / (r * r)).exp();
                    equilibrium(1.0, u0, vy as f32)
                };
                for i in 0..9 {
                    self.f[(y * 9 + i) * nx + x] = eq[i];
                }
            }
        }
        self.g.copy_from_slice(&self.f);
        self.time = 0;
        self.forces = ForceSum::default();
        self.init_sa();
    }

    /// `n` pas ; la vitesse d'entrée suit `ramp`.
    pub fn run(&mut self, n: usize, ramp: &Ramp) {
        for _ in 0..n {
            let (u, a) = ramp.at(self.time);
            self.u_in = u;
            self.accel = a;
            match self.pool.take() {
                Some(p) => {
                    p.install(|| self.step());
                    self.pool = Some(p);
                }
                None => self.step(),
            }
        }
    }

    /// Efforts cumulés depuis le dernier appel, puis remise à zéro.
    pub fn take_forces(&mut self) -> ForceSum {
        std::mem::take(&mut self.forces)
    }

    /// Efforts du dernier pas.
    pub fn last_forces(&self) -> (f64, f64, f64) {
        self.last
    }

    fn step(&mut self) {
        let (nx, ny) = (self.cfg.nx, self.cfg.ny);
        let rr = self.cfg.collision == Collision::Rr;
        let les = self.les_active;
        let ctx = StepCtx {
            nx,
            ny,
            n: nx * ny,
            src: &self.f,
            geo: &self.geo,
            tau_col: &self.tau_col,
            sigma_row: &self.sigma_row,
            c18: 18.0 * self.cfg.cs * self.cfg.cs,
            gx: self.accel,
            q_eq: equilibrium(1.0, self.u_in, 0.0),
            bouzidi: self.cfg.bouzidi,
            spin: self.spin,
            rot: (self.rot_x, self.rot_y),
            refp: (self.ref_x, self.ref_y),
            nu: self.cfg.nu,
            inv_nu: (1.0 / self.cfg.nu) as f32,
            nu_tilde: &[],
        };
        let sum = |a: (f64, f64, f64), b: (f64, f64, f64)| (a.0 + b.0, a.1 + b.1, a.2 + b.2);
        let (fx, fy, mz) = match self.sa.as_mut().filter(|_| self.sa_active) {
            Some(st) => {
                let ctx = StepCtx {
                    nu_tilde: &st.nu_tilde,
                    ..ctx
                };
                self.g
                    .par_chunks_mut(9 * nx)
                    .zip(st.usx.par_chunks_mut(nx))
                    .zip(st.usy.par_chunks_mut(nx))
                    .zip(st.utau.par_chunks_mut(nx))
                    .enumerate()
                    .map(|(y, (((row, usx), usy), utau))| {
                        let mut w = SaRow { usx, usy, utau };
                        match (rr, les) {
                            (true, true) => ctx.row::<true, true, true>(y, row, &mut w),
                            (true, false) => ctx.row::<true, false, true>(y, row, &mut w),
                            (false, true) => ctx.row::<false, true, true>(y, row, &mut w),
                            (false, false) => ctx.row::<false, false, true>(y, row, &mut w),
                        }
                    })
                    .reduce(|| (0.0, 0.0, 0.0), sum)
            }
            None => self
                .g
                .par_chunks_mut(9 * nx)
                .enumerate()
                .map(|(y, row)| {
                    let mut w = SaRow::default();
                    match (rr, les) {
                        (true, true) => ctx.row::<true, true, false>(y, row, &mut w),
                        (true, false) => ctx.row::<true, false, false>(y, row, &mut w),
                        (false, true) => ctx.row::<false, true, false>(y, row, &mut w),
                        (false, false) => ctx.row::<false, false, false>(y, row, &mut w),
                    }
                })
                .reduce(|| (0.0, 0.0, 0.0), sum),
        };
        self.open_boundaries();
        if self.sa_active {
            if let Some(st) = self.sa.as_mut() {
                st.step(&self.geo, self.cfg.nu, &self.sigma_row);
            }
        }
        std::mem::swap(&mut self.f, &mut self.g);
        self.last = (fx, fy, mz);
        self.forces.fx += fx;
        self.forces.fy += fy;
        self.forces.mz += mz;
        self.forces.steps += 1;
        self.time += 1;
    }

    /// Entrée (x = 0) : vitesse imposée, densité extrapolée. Sortie : ρ = 1,
    /// vitesse extrapolée. Partie hors équilibre copiée du voisin (Guo 2002).
    fn open_boundaries(&mut self) {
        let (nx, ny) = (self.cfg.nx, self.cfg.ny);
        let u = self.u_in;
        let dst = &mut self.g;
        let mut sa = self.sa.as_mut().filter(|_| self.sa_active);
        for y in 0..ny {
            for (b, nb, inlet) in [(0usize, 1usize, true), (nx - 1, nx - 2, false)] {
                let fnb: [f32; 9] = std::array::from_fn(|i| dst[(y * 9 + i) * nx + nb]);
                let (rho, ux, uy) = moments(&fnb);
                let eb = if inlet {
                    equilibrium(rho, u, 0.0)
                } else {
                    equilibrium(1.0, ux, uy)
                };
                let en = equilibrium(rho, ux, uy);
                for i in 0..9 {
                    dst[(y * 9 + i) * nx + b] = eb[i] + fnb[i] - en[i];
                }
                if let Some(st) = sa.as_deref_mut() {
                    let (vx, vy) = if inlet { (u, 0.0) } else { (ux, uy) };
                    st.usx[y * nx + b] = vx;
                    st.usy[y * nx + b] = vy;
                }
            }
        }
    }

    /// ρ, u (mailles solides : repos).
    pub fn macros(&self, out: &mut Macros) {
        let (nx, ny) = (self.cfg.nx, self.cfg.ny);
        let n = nx * ny;
        out.nx = nx;
        out.ny = ny;
        out.rho.resize(n, 0.0);
        out.ux.resize(n, 0.0);
        out.uy.resize(n, 0.0);
        let f = &self.f;
        let kind = &self.geo.kind;
        out.rho
            .par_chunks_mut(nx)
            .zip(out.ux.par_chunks_mut(nx))
            .zip(out.uy.par_chunks_mut(nx))
            .enumerate()
            .for_each(|(y, ((r, ux), uy))| {
                for x in 0..nx {
                    if kind[y * nx + x] == SOLID {
                        r[x] = 1.0;
                        ux[x] = 0.0;
                        uy[x] = 0.0;
                        continue;
                    }
                    let fl: [f32; 9] = std::array::from_fn(|i| f[(y * 9 + i) * nx + x]);
                    let (a, b, c) = moments(&fl);
                    r[x] = a;
                    ux[x] = b;
                    uy[x] = c;
                }
            });
    }

    /// Divergence ? (efforts et quelques densités)
    pub fn is_finite(&self) -> bool {
        let (fx, fy, _) = self.last;
        if !fx.is_finite() || !fy.is_finite() {
            return false;
        }
        let (nx, ny) = (self.cfg.nx, self.cfg.ny);
        let step = (nx * ny / 97).max(1);
        let mut k = 0;
        while k < nx * ny {
            let (x, y) = (k % nx, k / nx);
            let r: f32 = (0..9).map(|i| self.f[(y * 9 + i) * nx + x]).sum();
            if !(r > 0.3 && r < 3.0) {
                return false;
            }
            k += step;
        }
        true
    }
}

/// Données en lecture seule d'un pas (partagées entre les fils).
struct StepCtx<'a> {
    nx: usize,
    ny: usize,
    n: usize,
    src: &'a [f32],
    geo: &'a Geometry,
    tau_col: &'a [f32],
    sigma_row: &'a [f32],
    c18: f32,
    gx: f32,
    q_eq: [f32; 9],
    bouzidi: bool,
    spin: f64,
    rot: (f64, f64),
    refp: (f64, f64),
    /// Viscosité moléculaire, son inverse et ν̃ (Spalart–Allmaras).
    nu: f64,
    inv_nu: f32,
    nu_tilde: &'a [f32],
}

/// Champs de Spalart–Allmaras d'une ligne, écrits pendant le pas.
#[derive(Default)]
struct SaRow<'a> {
    usx: &'a mut [f32],
    usy: &'a mut [f32],
    utau: &'a mut [f32],
}

impl StepCtx<'_> {
    #[inline(always)]
    fn at(&self, x: usize, y: usize, i: usize) -> f32 {
        self.src[(y * 9 + i) * self.nx + x]
    }

    /// Une ligne ; renvoie les efforts (Fx, Fy, Mz) de ses mailles pariétales.
    fn row<const RR: bool, const LES: bool, const SA: bool>(
        &self,
        y: usize,
        row: &mut [f32],
        w: &mut SaRow,
    ) -> (f64, f64, f64) {
        let nx = self.nx;
        let sigma = self.sigma_row[y];
        let damp = &self.geo.les_damp[y * nx..(y + 1) * nx];
        for &(x0, x1) in &self.geo.runs[y] {
            self.fluid_run::<RR, LES, SA>(y, row, x0, x1, damp, sigma, w);
        }
        let mut force = (0.0, 0.0, 0.0);
        for &x in &self.geo.near[y] {
            let idx = y * nx + x;
            debug_assert_eq!(self.geo.kind[idx], NEAR);
            let f = self.pull_near_wall::<SA>(x, y, &mut force, w);
            let mut tau = self.tau_col[x];
            if SA {
                let (_, ux, uy) = moments(&f);
                w.usx[x] = ux;
                w.usy[x] = uy;
                tau += 3.0 * sa::nu_turb(self.nu_tilde[idx], self.inv_nu);
            }
            let r = Relax {
                tau,
                c18: if LES { self.c18 * damp[x] } else { 0.0 },
                gx: self.gx,
                sigma,
            };
            let p = collide::<RR, LES>(&f, r, &self.q_eq);
            for i in 0..9 {
                row[i * nx + x] = p[i];
            }
        }
        force
    }

    /// Mailles fluides ordinaires [x0, x1) : boucle sans branchement.
    #[inline(always)]
    #[allow(clippy::too_many_arguments)]
    fn fluid_run<const RR: bool, const LES: bool, const SA: bool>(
        &self,
        y: usize,
        row: &mut [f32],
        x0: usize,
        x1: usize,
        damp: &[f32],
        sigma: f32,
        w: &mut SaRow,
    ) {
        let nx = self.nx;
        let n = x1 - x0;
        let src = self.src;
        // Plan i de la ligne y − c_iy, décalé de −c_ix.
        let s = |i: usize| -> &[f32] {
            let yy = (y as i64 - CY[i] as i64) as usize;
            let a = ((yy * 9 + i) * nx) as i64 + x0 as i64 - CX[i] as i64;
            &src[a as usize..a as usize + n]
        };
        let (s0, s1, s2, s3, s4, s5, s6, s7, s8) = (s(0), s(1), s(2), s(3), s(4), s(5), s(6), s(7), s(8));
        let mut planes = row.chunks_exact_mut(nx);
        let mut d = || &mut planes.next().unwrap()[x0..x1];
        let (d0, d1, d2, d3, d4, d5, d6, d7, d8) = (d(), d(), d(), d(), d(), d(), d(), d(), d());
        let tau = &self.tau_col[x0..x1];
        let dm = &damp[x0..x1];
        let q = &self.q_eq;
        let (c18, gx) = (self.c18, self.gx);
        let (usx, usy, nt): (&mut [f32], &mut [f32], &[f32]) = if SA {
            (
                &mut w.usx[x0..x1],
                &mut w.usy[x0..x1],
                &self.nu_tilde[y * nx + x0..y * nx + x1],
            )
        } else {
            (&mut [], &mut [], &[])
        };
        let inv_nu = self.inv_nu;
        for k in 0..n {
            let f = [s0[k], s1[k], s2[k], s3[k], s4[k], s5[k], s6[k], s7[k], s8[k]];
            let mut tk = tau[k];
            if SA {
                // Vitesse avant collision (loi de paroi, transport de ν̃) et
                // viscosité turbulente ajoutée au temps de relaxation.
                let (_, ux, uy) = moments(&f);
                usx[k] = ux;
                usy[k] = uy;
                tk += 3.0 * sa::nu_turb(nt[k], inv_nu);
            }
            let r = Relax {
                tau: tk,
                c18: if LES { c18 * dm[k] } else { 0.0 },
                gx,
                sigma,
            };
            let p = collide::<RR, LES>(&f, r, q);
            d0[k] = p[0];
            d1[k] = p[1];
            d2[k] = p[2];
            d3[k] = p[3];
            d4[k] = p[4];
            d5[k] = p[5];
            d6[k] = p[6];
            d7[k] = p[7];
            d8[k] = p[8];
        }
    }

    /// Propagation d'une maille voisine de la paroi : rebond interpolé de
    /// Bouzidi (paroi mobile si le cylindre tourne), réflexion spéculaire en
    /// haut et en bas, efforts par échange de quantité de mouvement.
    fn pull_near_wall<const SA: bool>(
        &self,
        x: usize,
        y: usize,
        force: &mut (f64, f64, f64),
        w: &mut SaRow,
    ) -> [f32; 9] {
        let (nx, ny, n) = (self.nx, self.ny, self.n);
        let idx = y * nx + x;
        let mut slip: Option<(f64, f64)> = None;
        let mut out = [0.0f32; 9];
        out[0] = self.at(x, y, 0);
        for j in 1..9 {
            let (cx, cy) = (CX[j], CY[j]);
            let sy = y as i64 - cy as i64;
            let sx = (x as i64 - cx as i64) as usize;
            if sy < 0 || sy >= ny as i64 {
                out[j] = self.at(sx, y, REFY[j]);
                continue;
            }
            let q = self.geo.wall_q[j * n + idx];
            if q == 0.0 {
                out[j] = self.at(sx, sy as usize, j);
                continue;
            }
            // Lien vers la paroi : i = opp(j), c_i = −c_j.
            let i = OPP[j];
            let fi = self.at(x, y, i);
            let wx = x as f64 - q as f64 * cx as f64;
            let wy = y as f64 - q as f64 * cy as f64;
            let (mut uwx, mut uwy) = (0.0f64, 0.0f64);
            if self.spin != 0.0 {
                uwx = -self.spin * (wy - self.rot.1);
                uwy = self.spin * (wx - self.rot.0);
            }
            if SA {
                // Loi de paroi : glissement tangentiel (vitesse du pas précédent).
                let (sx, sy) = *slip.get_or_insert_with(|| self.wall_slip(idx, x, y, w));
                uwx += sx;
                uwy += sy;
            }
            let corr = 6.0 * W[j] * (cx as f32 * uwx as f32 + cy as f32 * uwy as f32);
            let val = if !self.bouzidi {
                fi + corr
            } else if q < 0.5 {
                let bx = x as i64 + cx as i64;
                let by = y as i64 + cy as i64;
                if by >= 0
                    && by < ny as i64
                    && bx >= 0
                    && bx < nx as i64
                    && self.geo.kind[by as usize * nx + bx as usize] != SOLID
                {
                    2.0 * q * fi + (1.0 - 2.0 * q) * self.at(bx as usize, by as usize, i) + corr
                } else {
                    fi + corr
                }
            } else {
                (fi + (2.0 * q - 1.0) * self.at(x, y, j) + corr) / (2.0 * q)
            };
            out[j] = val;
            let sum = (fi + val) as f64;
            let dif = (fi - val) as f64;
            let fxl = -cx as f64 * sum - uwx * dif;
            let fyl = -cy as f64 * sum - uwy * dif;
            force.0 += fxl;
            force.1 += fyl;
            force.2 += (wx - self.refp.0) * fyl - (wy - self.refp.1) * fxl;
        }
        out
    }

    /// Glissement de la loi de paroi d'une maille pariétale ; met à jour u_τ.
    fn wall_slip(&self, idx: usize, x: usize, y: usize, w: &mut SaRow) -> (f64, f64) {
        let geo = self.geo;
        let mut ux = w.usx[x] as f64;
        let mut uy = w.usy[x] as f64;
        if self.spin != 0.0 {
            // Vitesse relative à la paroi tournante.
            ux += self.spin * (y as f64 - self.rot.1);
            uy -= self.spin * (x as f64 - self.rot.0);
        }
        let normal = (geo.normal[2 * idx] as f64, geo.normal[2 * idx + 1] as f64);
        let yw = (geo.dist[idx] as f64).max(0.1);
        let (sx, sy, utau) = sa::wall_slip(ux, uy, normal, yw, geo.lam[idx] as f64, self.nu);
        w.utau[x] = utau as f32;
        (sx, sy)
    }
}
