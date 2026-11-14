//! Distributions pariétales du calcul LBM : coefficients de pression Cp et
//! de frottement Cf le long du corps, moyennés dans le temps, à partir de ρ et
//! u des mailles voisines de la paroi (les deux moteurs).

use crate::cpu::Macros;
use crate::sa::spalding_utau;
use crate::shapes::ShapeKind;
use crate::tunnel::WindTunnel;

/// Abscisse le long de la paroi.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WallAxis {
    /// Angle depuis le point d'arrêt avant (cylindre), en degrés.
    Angle,
    /// x/c (profils, ellipses, plaques).
    Chord,
}

/// Formes pour lesquelles une abscisse le long de la paroi a un sens simple.
pub fn wall_supported(kind: ShapeKind) -> bool {
    matches!(
        kind,
        ShapeKind::Naca | ShapeKind::Ellipse | ShapeKind::Plate | ShapeKind::Cylinder
    )
}

#[derive(Clone, Copy, Debug)]
pub struct WallPoint {
    pub s: f64,
    pub cp: f64,
    pub cf: f64,
}

/// Moyennes, triées par abscisse, par face.
#[derive(Clone, Debug)]
pub struct WallProfile {
    pub axis: WallAxis,
    /// Nombre de relevés moyennés.
    pub n: usize,
    pub upper: Vec<WallPoint>,
    pub lower: Vec<WallPoint>,
}

struct Cell {
    idx: usize,
    s: f64,
    upper: bool,
    d: f64,
    nx: f64,
    ny: f64,
}

pub struct WallSampler {
    pub axis: WallAxis,
    /// Axe x du corps dans le repère du réseau (sens du frottement positif).
    ex: f64,
    ey: f64,
    cells: Vec<Cell>,
    sum_cp: Vec<f64>,
    sum_cf: Vec<f64>,
    n: usize,
}

impl WallSampler {
    pub fn new(wt: &WindTunnel) -> WallSampler {
        let geo = wt.engine.wall_geometry();
        let body = wt.wall_body();
        let axis = if body.kind == ShapeKind::Cylinder {
            WallAxis::Angle
        } else {
            WallAxis::Chord
        };
        let loc = &body.outline_local;
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
        let by_index = body.kind == ShapeKind::Naca;
        let a = body.alpha_deg.to_radians();
        let mut cells = Vec::new();
        for &idx in &geo.wall_cells {
            let pos = geo.near_pos[idx] as f64;
            if pos < 0.0 || m == 0 {
                continue;
            }
            let k = (pos.floor() as usize).min(m - 1);
            let t = pos - k as f64;
            let [ax, ay] = loc[k];
            let [bx2, by2] = loc[(k + 1) % m];
            let bx = ax + t * (bx2 - ax);
            let by = ay + t * (by2 - ay);
            let upper = if by_index {
                k < i_le || (k == i_le && t == 0.0)
            } else {
                by >= 0.0
            };
            let s = match axis {
                WallAxis::Angle => 180.0 - by.atan2(bx).abs().to_degrees(),
                WallAxis::Chord => (bx - xmin) / span,
            };
            cells.push(Cell {
                idx,
                s,
                upper,
                d: geo.dist[idx] as f64,
                nx: geo.normal[2 * idx] as f64,
                ny: geo.normal[2 * idx + 1] as f64,
            });
        }
        let n = cells.len();
        WallSampler {
            axis,
            ex: a.cos(),
            ey: -a.sin(),
            cells,
            sum_cp: vec![0.0; n],
            sum_cf: vec![0.0; n],
            n: 0,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.cells.is_empty()
    }

    pub fn reset(&mut self) {
        self.sum_cp.iter_mut().for_each(|v| *v = 0.0);
        self.sum_cf.iter_mut().for_each(|v| *v = 0.0);
        self.n = 0;
    }

    /// Ajoute un relevé : `m`, champs macroscopiques complets de la veine
    /// (pression de référence) ; avec des grilles fines, la paroi est relevée
    /// dans la plus fine.
    pub fn add(&mut self, wt: &mut WindTunnel, m: &Macros) {
        if m.rho.len() != m.nx * m.ny || m.nx == 0 {
            return;
        }
        let rho_ref = upstream_density(m);
        let fine;
        let m = if wt.wall_l() != wt.cfg.l {
            let mut f = Macros::default();
            wt.engine.wall_macros(&mut f);
            fine = f;
            &fine
        } else {
            m
        };
        if m.rho.len() != m.nx * m.ny {
            return;
        }
        let geo = wt.engine.wall_geometry();
        let u = wt.cfg.u_in as f64;
        let nu = u * wt.wall_l() as f64 / wt.cfg.re;
        let sa = wt.sa_active();
        let q = 0.5 * rho_ref * u * u;
        for (k, c) in self.cells.iter().enumerate() {
            let i = c.idx;
            let (rho, ux, uy) = (m.rho[i] as f64, m.ux[i] as f64, m.uy[i] as f64);
            // p = ρ c_s², c_s² = 1/3.
            self.sum_cp[k] += (rho - rho_ref) / (3.0 * q);
            let un = ux * c.nx + uy * c.ny;
            let (tx, ty) = (ux - un * c.nx, uy - un * c.ny);
            let ut = tx.hypot(ty);
            let sign = if tx * self.ex + ty * self.ey >= 0.0 { 1.0 } else { -1.0 };
            let tau = if sa {
                let cf_lam = geo.lam[i] as f64;
                if cf_lam > 0.0 {
                    0.5 * cf_lam * ut * ut
                } else {
                    let ut2 = spalding_utau(ut, c.d.max(0.1), nu);
                    ut2 * ut2
                }
            } else {
                // Profil linéaire entre la paroi (u = 0) et la maille.
                nu * ut / c.d.max(0.3)
            };
            self.sum_cf[k] += sign * rho * tau / q;
        }
        self.n += 1;
    }

    pub fn result(&self) -> WallProfile {
        let mut out = WallProfile {
            axis: self.axis,
            n: self.n,
            upper: Vec::new(),
            lower: Vec::new(),
        };
        if self.n == 0 {
            return out;
        }
        let n = self.n as f64;
        for (k, c) in self.cells.iter().enumerate() {
            let p = WallPoint {
                s: c.s,
                cp: self.sum_cp[k] / n,
                cf: self.sum_cf[k] / n,
            };
            if c.upper {
                out.upper.push(p)
            } else {
                out.lower.push(p)
            }
        }
        out.upper.sort_by(|a, b| a.s.total_cmp(&b.s));
        out.lower.sort_by(|a, b| a.s.total_cmp(&b.s));
        out
    }
}

/// Masse volumique moyenne de la moitié centrale de la colonne x = 3 (référence de pression).
pub fn upstream_density(m: &Macros) -> f64 {
    let (nx, ny) = (m.nx, m.ny);
    if nx < 4 || m.rho.len() < nx * ny {
        return 1.0;
    }
    let (j0, j1) = (ny / 4, (3 * ny).div_ceil(4));
    let k = j1.saturating_sub(j0);
    if k == 0 {
        return 1.0;
    }
    (j0..j1).map(|j| m.rho[j * nx + 3] as f64).sum::<f64>() / k as f64
}
