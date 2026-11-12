//! Géométrie du corps sur le réseau, commune aux deux moteurs : mailles
//! solides, liens de paroi interpolés (fraction q de Bouzidi, y compris à
//! travers les parois plus minces qu'une maille), distance à la paroi et
//! amortissement pariétal de la sous-maille.

use crate::d2q9::{CX, CY};
use crate::shapes::Body;

pub const FLUID: u8 = 0;
pub const SOLID: u8 = 1;
/// Fluide traité par le chemin lent : voisin d'un solide, ou ligne haute/basse.
pub const NEAR: u8 = 2;

pub struct Geometry {
    pub nx: usize,
    pub ny: usize,
    pub n: usize,
    pub kind: Vec<u8>,
    /// q du lien tiré dans la direction j : `wall_q[j * n + idx]` ; 0 si pas de paroi.
    pub wall_q: Vec<f32>,
    /// Distance (mailles) à la paroi du corps.
    pub dist: Vec<f32>,
    /// Amortissement de Van Driest de la sous-maille, (1 − e^(−d/A))².
    pub les_damp: Vec<f32>,
    pub les_wall_a: f64,
    /// Par ligne : intervalles [x0, x1) de mailles fluides ordinaires.
    pub runs: Vec<Vec<(usize, usize)>>,
    /// Par ligne : mailles du chemin lent (x).
    pub near: Vec<Vec<usize>>,
    /// Mailles portant au moins un lien de paroi avec le corps.
    pub wall_cells: Vec<usize>,
    pub solid_count: usize,
    pub link_count: usize,
}

impl Geometry {
    pub fn new(nx: usize, ny: usize, les_wall_a: f64) -> Geometry {
        let n = nx * ny;
        let mut g = Geometry {
            nx,
            ny,
            n,
            kind: vec![FLUID; n],
            wall_q: vec![0.0; 9 * n],
            dist: vec![1e6; n],
            les_damp: vec![1.0; n],
            les_wall_a,
            runs: Vec::new(),
            near: Vec::new(),
            wall_cells: Vec::new(),
            solid_count: 0,
            link_count: 0,
        };
        for x in 0..nx {
            g.kind[x] = NEAR;
            g.kind[(ny - 1) * nx + x] = NEAR;
        }
        g.index_rows();
        g
    }

    /// Construit la géométrie ; renvoie les mailles solides redevenues fluides.
    pub fn build(&mut self, body: &Body) -> Vec<usize> {
        let (nx, ny, n) = (self.nx, self.ny, self.n);
        let mut revived = Vec::new();
        let mut solids = 0;
        let inside = |x: f64, y: f64| body.inside(x, y);
        for y in 0..ny {
            for x in 0..nx {
                let idx = y * nx + x;
                let solid = x >= 2 && x + 3 <= nx && y >= 2 && y + 3 <= ny && inside(x as f64, y as f64);
                if solid {
                    self.kind[idx] = SOLID;
                    solids += 1;
                } else {
                    if self.kind[idx] == SOLID {
                        revived.push(idx);
                    }
                    self.kind[idx] = FLUID;
                }
            }
        }
        self.solid_count = solids;
        let [bx0, bx1, by0, by1] = body.bbox();
        let (bx0, bx1) = (bx0.floor() - 2.0, bx1.ceil() + 2.0);
        let (by0, by1) = (by0.floor() - 2.0, by1.ceil() + 2.0);
        self.wall_q.iter_mut().for_each(|q| *q = 0.0);
        let mut links = 0;
        for y in 0..ny {
            for x in 1..nx - 1 {
                let idx = y * nx + x;
                if self.kind[idx] == SOLID {
                    continue;
                }
                let mut near = y == 0 || y == ny - 1;
                let (xf, yf) = (x as f64, y as f64);
                let in_box = xf >= bx0 && xf <= bx1 && yf >= by0 && yf <= by1;
                for j in 1..9 {
                    let sx = x as i64 - CX[j] as i64;
                    let sy = y as i64 - CY[j] as i64;
                    if sy < 0 || sy >= ny as i64 {
                        continue;
                    }
                    let (dx, dy) = (-CX[j] as f64, -CY[j] as f64);
                    let q = if self.kind[sy as usize * nx + sx as usize] == SOLID {
                        link_fraction(&inside, xf, yf, dx, dy)
                    } else if in_box {
                        // Lien fluide–fluide qui traverse une paroi mince.
                        crossing_fraction(&inside, xf, yf, dx, dy)
                    } else {
                        0.0
                    };
                    if q > 0.0 {
                        near = true;
                        self.wall_q[j * n + idx] = q as f32;
                        links += 1;
                    }
                }
                if near {
                    self.kind[idx] = NEAR;
                }
            }
        }
        // Colonnes d'entrée et de sortie : traitées à part.
        for y in 0..ny {
            self.kind[y * nx] = FLUID;
            self.kind[y * nx + nx - 1] = FLUID;
        }
        self.link_count = links;
        self.wall_distance(body);
        self.update_les_damping();
        self.index_rows();
        revived
    }

    /// Distance exacte au contour près du corps, distance de chanfrein ailleurs.
    fn wall_distance(&mut self, body: &Body) {
        let (nx, ny, n) = (self.nx, self.ny, self.n);
        // Chanfrein (1, √2) aux mailles solides, deux passes.
        let big = 1e6f32;
        let d2 = std::f32::consts::SQRT_2;
        let mut d: Vec<f32> = (0..n).map(|i| if self.kind[i] == SOLID { 0.0 } else { big }).collect();
        for y in 0..ny {
            for x in 0..nx {
                let i = y * nx + x;
                let mut v = d[i];
                if x > 0 {
                    v = v.min(d[i - 1] + 1.0);
                }
                if y > 0 {
                    v = v.min(d[i - nx] + 1.0);
                    if x > 0 {
                        v = v.min(d[i - nx - 1] + d2);
                    }
                    if x + 1 < nx {
                        v = v.min(d[i - nx + 1] + d2);
                    }
                }
                d[i] = v;
            }
        }
        for y in (0..ny).rev() {
            for x in (0..nx).rev() {
                let i = y * nx + x;
                let mut v = d[i];
                if x + 1 < nx {
                    v = v.min(d[i + 1] + 1.0);
                }
                if y + 1 < ny {
                    v = v.min(d[i + nx] + 1.0);
                    if x + 1 < nx {
                        v = v.min(d[i + nx + 1] + d2);
                    }
                    if x > 0 {
                        v = v.min(d[i + nx - 1] + d2);
                    }
                }
                d[i] = v;
            }
        }
        let has_solid = self.solid_count > 0;
        for i in 0..n {
            self.dist[i] = if self.kind[i] == SOLID {
                0.0
            } else if has_solid {
                (d[i] - 0.5).max(0.5)
            } else {
                1e6
            };
        }
        // Distance exacte au polygone dans une boîte autour du corps.
        let o = &body.outline;
        let m = o.len();
        if m < 2 {
            return;
        }
        let [bx0, bx1, by0, by1] = body.bbox();
        let margin = 8f64.max(0.6 * (bx1 - bx0).max(by1 - by0));
        let x0 = (bx0 - margin).floor().max(0.0) as usize;
        let x1 = ((bx1 + margin).ceil() as usize).min(nx - 1);
        let y0 = (by0 - margin).floor().max(0.0) as usize;
        let y1 = ((by1 + margin).ceil() as usize).min(ny - 1);
        for y in y0..=y1 {
            for x in x0..=x1 {
                let i = y * nx + x;
                if self.kind[i] == SOLID {
                    continue;
                }
                let (px, py) = (x as f64, y as f64);
                let mut best = f64::INFINITY;
                for k in 0..m {
                    let [ax, ay] = o[k];
                    let [cx, cy] = o[(k + 1) % m];
                    let (ex, ey) = (cx - ax, cy - ay);
                    let l2 = (ex * ex + ey * ey).max(1e-12);
                    let t = (((px - ax) * ex + (py - ay) * ey) / l2).clamp(0.0, 1.0);
                    let dx = px - ax - t * ex;
                    let dy = py - ay - t * ey;
                    best = best.min(dx * dx + dy * dy);
                }
                self.dist[i] = (best.sqrt() as f32).max(1e-3);
            }
        }
    }

    pub fn set_les_wall_length(&mut self, a: f64) {
        self.les_wall_a = a;
        self.update_les_damping();
    }

    fn update_les_damping(&mut self) {
        let a = self.les_wall_a as f32;
        if !(a > 0.0) || self.solid_count == 0 {
            self.les_damp.iter_mut().for_each(|v| *v = 1.0);
            return;
        }
        for (v, &d) in self.les_damp.iter_mut().zip(&self.dist) {
            let f = 1.0 - (-d / a).exp();
            *v = f * f;
        }
    }

    /// Index des lignes : intervalles fluides (chemin rapide) et mailles du chemin lent.
    fn index_rows(&mut self) {
        let (nx, ny) = (self.nx, self.ny);
        self.runs = vec![Vec::new(); ny];
        self.near = vec![Vec::new(); ny];
        self.wall_cells.clear();
        for y in 0..ny {
            let mut start: Option<usize> = None;
            for x in 1..nx - 1 {
                let idx = y * nx + x;
                let k = self.kind[idx];
                if k == FLUID {
                    if start.is_none() {
                        start = Some(x);
                    }
                } else {
                    if let Some(s) = start.take() {
                        self.runs[y].push((s, x));
                    }
                    if k == NEAR {
                        self.near[y].push(x);
                        if (1..9).any(|j| self.wall_q[j * self.n + idx] > 0.0) {
                            self.wall_cells.push(idx);
                        }
                    }
                }
            }
            if let Some(s) = start {
                self.runs[y].push((s, nx - 1));
            }
        }
    }
}

/// Fraction q ∈ ]0, 1] du lien (x, y) → (x + dx, y + dy) située côté fluide.
fn link_fraction(inside: &impl Fn(f64, f64) -> bool, x: f64, y: f64, dx: f64, dy: f64) -> f64 {
    let (mut lo, mut hi) = (0.0, 1.0);
    for _ in 0..14 {
        let mid = 0.5 * (lo + hi);
        if inside(x + mid * dx, y + mid * dy) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    (0.5 * (lo + hi)).clamp(0.02, 1.0)
}

/// Lien entre deux nœuds fluides : fraction de la première intersection avec le corps, ou 0.
fn crossing_fraction(inside: &impl Fn(f64, f64) -> bool, x: f64, y: f64, dx: f64, dy: f64) -> f64 {
    let n = 16;
    for k in 1..n {
        let t = k as f64 / n as f64;
        if inside(x + t * dx, y + t * dy) {
            let (mut lo, mut hi) = ((k - 1) as f64 / n as f64, t);
            for _ in 0..10 {
                let mid = 0.5 * (lo + hi);
                if inside(x + mid * dx, y + mid * dy) {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            return (0.5 * (lo + hi)).clamp(0.02, 1.0);
        }
    }
    0.0
}
