//! Maquettes 2D : repère propre adimensionné par la longueur de référence L,
//! origine au point de réduction des moments ; placement dans le réseau avec
//! une incidence α (positive = cabrée, bord d'attaque vers le haut).

use std::f64::consts::PI;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ShapeKind {
    Cylinder,
    Square,
    Ellipse,
    Plate,
    Naca,
    Triangle,
}

impl ShapeKind {
    pub const ALL: [ShapeKind; 6] = [
        ShapeKind::Cylinder,
        ShapeKind::Square,
        ShapeKind::Ellipse,
        ShapeKind::Plate,
        ShapeKind::Naca,
        ShapeKind::Triangle,
    ];

    pub fn id(self) -> &'static str {
        match self {
            ShapeKind::Cylinder => "cylinder",
            ShapeKind::Square => "square",
            ShapeKind::Ellipse => "ellipse",
            ShapeKind::Plate => "plate",
            ShapeKind::Naca => "naca",
            ShapeKind::Triangle => "triangle",
        }
    }

    pub fn from_id(s: &str) -> Option<ShapeKind> {
        ShapeKind::ALL.into_iter().find(|k| k.id() == s)
    }

    pub fn label(self) -> &'static str {
        match self {
            ShapeKind::Cylinder => "Cylindre",
            ShapeKind::Square => "Carré",
            ShapeKind::Ellipse => "Ellipse",
            ShapeKind::Plate => "Plaque plane",
            ShapeKind::Naca => "Profil NACA",
            ShapeKind::Triangle => "Triangle",
        }
    }

    /// Plage d'incidence (degrés) [min, max], ou None.
    pub fn alpha_range(self) -> Option<(f64, f64)> {
        match self {
            ShapeKind::Cylinder => None,
            ShapeKind::Square => Some((-45.0, 45.0)),
            ShapeKind::Ellipse | ShapeKind::Plate => Some((-90.0, 90.0)),
            ShapeKind::Naca => Some((-20.0, 25.0)),
            ShapeKind::Triangle => Some((-180.0, 180.0)),
        }
    }

    /// Point local du bord de fuite (recalé sur le réseau).
    fn trailing_edge(self) -> Option<[f64; 2]> {
        match self {
            ShapeKind::Naca => Some([0.75, 0.0]),
            ShapeKind::Plate => Some([0.5, 0.0]),
            _ => None,
        }
    }
}

/// Profil NACA à 4 chiffres.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Naca {
    pub m: f64,
    pub p: f64,
    pub t: f64,
}

impl Naca {
    pub fn parse(code: &str) -> Option<Naca> {
        let s = code.trim();
        if s.len() != 4 || !s.bytes().all(|b| b.is_ascii_digit()) {
            return None;
        }
        let d: Vec<f64> = s.bytes().map(|b| (b - b'0') as f64).collect();
        let m = d[0] / 100.0;
        let p = d[1] / 10.0;
        let t = (d[2] * 10.0 + d[3]) / 100.0;
        if !(0.04..=0.4).contains(&t) || (m > 0.0 && (p <= 0.0 || p >= 1.0)) {
            return None;
        }
        Some(Naca { m, p, t })
    }

    pub fn thickness(&self, x: f64) -> f64 {
        // Bord de fuite fermé (coefficient −0,1036).
        5.0 * self.t * (0.2969 * x.sqrt() - 0.126 * x - 0.3516 * x * x + 0.2843 * x.powi(3) - 0.1036 * x.powi(4))
    }

    pub fn camber(&self, x: f64) -> f64 {
        let (m, p) = (self.m, self.p);
        if m == 0.0 {
            0.0
        } else if x < p {
            m / (p * p) * (2.0 * p * x - x * x)
        } else {
            m / ((1.0 - p) * (1.0 - p)) * (1.0 - 2.0 * p + 2.0 * p * x - x * x)
        }
    }

    pub fn code(&self) -> String {
        format!(
            "{}{}{:02}",
            (self.m * 100.0).round(),
            (self.p * 10.0).round(),
            (self.t * 100.0).round()
        )
    }
}

/// Paramètres de forme.
#[derive(Clone, Debug)]
pub struct ShapeParams {
    pub alpha_deg: f64,
    pub naca: Naca,
    /// Épaisseur relative de l'ellipse.
    pub thickness: f64,
    /// Cylindre : vitesse périphérique / U (> 0 : portance positive).
    pub spin_ratio: f64,
}

impl Default for ShapeParams {
    fn default() -> Self {
        ShapeParams {
            alpha_deg: 0.0,
            naca: Naca::parse("0012").unwrap(),
            thickness: 0.4,
            spin_ratio: 0.0,
        }
    }
}

/// Position de la pointe du bord de fuite par rapport au réseau (x entier,
/// y demi-entier), comme dans l'application web : la portance varie alors
/// régulièrement avec α.
const TE_ANCHOR: [f64; 2] = [0.0, 0.5];

const TRI_LEN: f64 = 0.866_025_403_784_438_6; // √3/2

fn triangle() -> [[f64; 2]; 3] {
    let v = [[-TRI_LEN / 2.0, 0.0], [TRI_LEN / 2.0, 0.5], [TRI_LEN / 2.0, -0.5]];
    let gx = (v[0][0] + v[1][0] + v[2][0]) / 3.0;
    v.map(|[x, y]| [x - gx, y])
}

fn circle(rx: f64, ry: f64, n: usize) -> Vec<[f64; 2]> {
    (0..n)
        .map(|k| {
            let a = k as f64 / n as f64 * 2.0 * PI;
            [rx * a.cos(), ry * a.sin()]
        })
        .collect()
}

/// Maquette placée dans le réseau.
#[derive(Clone, Debug)]
pub struct Body {
    pub kind: ShapeKind,
    pub l: f64,
    /// Point de référence (réseau).
    pub ref_x: f64,
    pub ref_y: f64,
    pub alpha_deg: f64,
    ca: f64,
    sa: f64,
    params: ShapeParams,
    plate_half: f64,
    /// Contour (réseau), sens horaire pour les profils.
    pub outline: Vec<[f64; 2]>,
    /// Contour dans le repère propre.
    pub outline_local: Vec<[f64; 2]>,
    /// Rotation du cylindre (rad/pas, > 0 sens trigonométrique).
    pub spin: f64,
}

impl Body {
    /// Place la forme : point de référence voisin de (x0, y0), longueur L (mailles).
    pub fn new(kind: ShapeKind, params: &ShapeParams, l: f64, x0: f64, y0: f64, u_in: f64) -> Body {
        let alpha_deg = if kind.alpha_range().is_some() {
            params.alpha_deg
        } else {
            0.0
        };
        let a = alpha_deg.to_radians();
        let (sa, ca) = a.sin_cos();
        let mut x0 = x0;
        let mut y0 = y0;
        if let Some([bx, by]) = kind.trailing_edge() {
            let [ax, ay] = TE_ANCHOR;
            let tx = x0 + l * (ca * bx + sa * by);
            let ty = y0 + l * (-sa * bx + ca * by);
            x0 += (tx - ax + 0.5).floor() + ax - tx;
            y0 += (ty - ay + 0.5).floor() + ay - ty;
        }
        let params = ShapeParams {
            thickness: params.thickness.clamp(0.1, 1.0),
            ..params.clone()
        };
        let plate_half = (1.25 / l).max(0.03);
        let local = match kind {
            ShapeKind::Cylinder => circle(0.5, 0.5, 96),
            ShapeKind::Square => vec![[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]],
            ShapeKind::Ellipse => circle(0.5, 0.5 * params.thickness, 96),
            ShapeKind::Plate => vec![
                [-0.5, -plate_half],
                [0.5, -plate_half],
                [0.5, plate_half],
                [-0.5, plate_half],
            ],
            ShapeKind::Triangle => triangle().to_vec(),
            ShapeKind::Naca => {
                let n = 80;
                let nc = params.naca;
                let mut upper = Vec::new();
                let mut lower = Vec::new();
                for k in 0..=n {
                    let xc = 0.5 * (1.0 - (k as f64 / n as f64 * PI).cos());
                    let yc = nc.camber(xc);
                    let yt = nc.thickness(xc);
                    upper.push([xc - 0.25, yc + yt]);
                    lower.push([xc - 0.25, yc - yt]);
                }
                upper.reverse();
                upper.extend_from_slice(&lower[1..]);
                upper
            }
        };
        let outline = local
            .iter()
            .map(|&[bx, by]| [x0 + l * (ca * bx + sa * by), y0 + l * (-sa * bx + ca * by)])
            .collect();
        let spin = if kind == ShapeKind::Cylinder {
            -params.spin_ratio * u_in / (l / 2.0)
        } else {
            0.0
        };
        Body {
            kind,
            l,
            ref_x: x0,
            ref_y: y0,
            alpha_deg,
            ca,
            sa,
            params,
            plate_half,
            outline,
            outline_local: local,
            spin,
        }
    }

    /// Test d'appartenance (coordonnées réseau, nœuds aux entiers).
    pub fn inside(&self, x: f64, y: f64) -> bool {
        let dx = (x - self.ref_x) / self.l;
        let dy = (y - self.ref_y) / self.l;
        let bx = self.ca * dx - self.sa * dy;
        let by = self.sa * dx + self.ca * dy;
        match self.kind {
            ShapeKind::Cylinder => bx * bx + by * by <= 0.25,
            ShapeKind::Square => bx.abs() <= 0.5 && by.abs() <= 0.5,
            ShapeKind::Ellipse => {
                let b = 0.5 * self.params.thickness;
                (bx / 0.5).powi(2) + (by / b).powi(2) <= 1.0
            }
            ShapeKind::Plate => bx.abs() <= 0.5 && by.abs() <= self.plate_half,
            ShapeKind::Naca => {
                let xc = bx + 0.25;
                if !(0.0..=1.0).contains(&xc) {
                    return false;
                }
                let nc = self.params.naca;
                (by - nc.camber(xc)).abs() <= nc.thickness(xc)
            }
            ShapeKind::Triangle => {
                let t = triangle();
                let mut sign = 0.0f64;
                for k in 0..3 {
                    let [ax, ay] = t[k];
                    let [cx, cy] = t[(k + 1) % 3];
                    let c = (cx - ax) * (by - ay) - (cy - ay) * (bx - ax);
                    if c != 0.0 {
                        let s = c.signum();
                        if sign == 0.0 {
                            sign = s;
                        } else if s != sign {
                            return false;
                        }
                    }
                }
                true
            }
        }
    }

    /// Paramètres de la forme.
    pub fn params(&self) -> &ShapeParams {
        &self.params
    }

    /// Boîte englobante (réseau).
    pub fn bbox(&self) -> [f64; 4] {
        let mut b = [f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY];
        for &[x, y] in &self.outline {
            b[0] = b[0].min(x);
            b[1] = b[1].max(x);
            b[2] = b[2].min(y);
            b[3] = b[3].max(y);
        }
        b
    }

    /// Désignation lisible.
    pub fn label(&self) -> String {
        match self.kind {
            ShapeKind::Naca => format!("NACA {}", self.params.naca.code()),
            ShapeKind::Ellipse => format!("Ellipse e/c = {:.2}", self.params.thickness),
            k => k.label().to_string(),
        }
    }
}
