//! Maquettes de la soufflerie 3D : domaine, test d'appartenance (voxelisation)
//! et surface de référence. Coordonnées réseau : x vers l'aval, y latéral
//! (envergure), z vers le haut. Portage de `js/bodies3d.js`.

/// Maquettes disponibles.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Shape3D {
    Sphere,
    Cube,
    Cylinder,
    Wing,
}

impl Shape3D {
    pub const ALL: [Shape3D; 4] = [Shape3D::Sphere, Shape3D::Cube, Shape3D::Cylinder, Shape3D::Wing];

    pub fn from_id(s: &str) -> Option<Shape3D> {
        match s {
            "sphere" => Some(Shape3D::Sphere),
            "cube" => Some(Shape3D::Cube),
            "cylinder" => Some(Shape3D::Cylinder),
            "wing" => Some(Shape3D::Wing),
            _ => None,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Shape3D::Sphere => "Sphère",
            Shape3D::Cube => "Cube",
            Shape3D::Cylinder => "Cylindre",
            Shape3D::Wing => "Aile",
        }
    }

    /// Longueur de référence.
    pub fn reference(self) -> &'static str {
        match self {
            Shape3D::Sphere => "diamètre D",
            Shape3D::Cube => "arête a",
            Shape3D::Cylinder => "diamètre D (envergure périodique 4D)",
            Shape3D::Wing => "corde c (NACA 0012, allongement 4)",
        }
    }

    /// Dimensions de la veine en longueurs de référence.
    pub fn domain(self) -> [f64; 3] {
        match self {
            Shape3D::Sphere | Shape3D::Cube => [12.0, 6.0, 6.0],
            Shape3D::Cylinder => [16.0, 4.0, 8.0],
            Shape3D::Wing => [9.0, 8.0, 4.0],
        }
    }

    /// Position de la maquette en longueurs de référence.
    fn at(self) -> [f64; 3] {
        match self {
            Shape3D::Sphere | Shape3D::Cube => [4.0, 3.0, 3.0],
            Shape3D::Cylinder => [5.0, 2.0, 4.0],
            Shape3D::Wing => [2.5, 4.0, 2.0],
        }
    }

    /// L'incidence est-elle réglable ?
    pub fn has_alpha(self) -> bool {
        matches!(self, Shape3D::Cube | Shape3D::Wing)
    }

    /// Nombre de mailles pour la longueur de référence `l`.
    pub fn dims(self, l: usize) -> [usize; 3] {
        let d = self.domain();
        [
            (d[0] * l as f64).round() as usize,
            (d[1] * l as f64).round() as usize,
            (d[2] * l as f64).round() as usize,
        ]
    }
}

/// Boîte englobante (coordonnées réseau).
#[derive(Clone, Copy, Debug)]
pub struct Bounds {
    pub x0: f64,
    pub x1: f64,
    pub y0: f64,
    pub y1: f64,
    pub z0: f64,
    pub z1: f64,
}

#[derive(Clone, Debug)]
pub struct Body3D {
    pub shape: Shape3D,
    pub l: usize,
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub center: [f64; 3],
    pub bounds: Bounds,
    /// Surface de référence (mailles²).
    pub area: f64,
    alpha: f64,
    /// Demi-envergure de l'aile (mailles).
    half: f64,
}

fn naca00(t: f64, x: f64) -> f64 {
    5.0 * t * (0.2969 * x.sqrt() - 0.126 * x - 0.3516 * x * x + 0.2843 * x.powi(3) - 0.1036 * x.powi(4))
}

impl Body3D {
    /// `alpha_deg` : incidence (cube, aile), cabrage bord de fuite vers le bas.
    pub fn new(shape: Shape3D, l: usize, alpha_deg: f64) -> Body3D {
        let [nx, ny, nz] = shape.dims(l);
        let at = shape.at();
        let lf = l as f64;
        // Centre décalé d'un quart de maille : évite les symétries d'arrondi exactes.
        let c = [at[0] * lf + 0.25, at[1] * lf + 0.25, at[2] * lf + 0.25];
        let alpha = if shape.has_alpha() { alpha_deg.to_radians() } else { 0.0 };
        let r = lf / 2.0;
        let half = 2.0 * lf; // allongement 4
        let (bounds, area) = match shape {
            Shape3D::Sphere => (
                Bounds {
                    x0: c[0] - r,
                    x1: c[0] + r,
                    y0: c[1] - r,
                    y1: c[1] + r,
                    z0: c[2] - r,
                    z1: c[2] + r,
                },
                std::f64::consts::PI * lf * lf / 4.0,
            ),
            Shape3D::Cube => {
                let e = r * std::f64::consts::SQRT_2;
                (
                    Bounds {
                        x0: c[0] - e,
                        x1: c[0] + e,
                        y0: c[1] - r,
                        y1: c[1] + r,
                        z0: c[2] - e,
                        z1: c[2] + e,
                    },
                    lf * lf,
                )
            }
            Shape3D::Cylinder => (
                Bounds {
                    x0: c[0] - r,
                    x1: c[0] + r,
                    y0: 0.0,
                    y1: ny as f64 - 1.0,
                    z0: c[2] - r,
                    z1: c[2] + r,
                },
                lf * ny as f64,
            ),
            Shape3D::Wing => {
                let rr = 0.8 * lf;
                (
                    Bounds {
                        x0: c[0] - rr,
                        x1: c[0] + rr,
                        y0: c[1] - half,
                        y1: c[1] + half,
                        z0: c[2] - rr,
                        z1: c[2] + rr,
                    },
                    lf * 2.0 * half,
                )
            }
        };
        Body3D {
            shape,
            l,
            nx,
            ny,
            nz,
            center: c,
            bounds,
            area,
            alpha,
            half,
        }
    }

    /// Le point (réseau) est-il dans la maquette ?
    pub fn inside(&self, x: f64, y: f64, z: f64) -> bool {
        let [cx, cy, cz] = self.center;
        let lf = self.l as f64;
        let r = lf / 2.0;
        let (sa, ca) = self.alpha.sin_cos();
        match self.shape {
            Shape3D::Sphere => (x - cx).powi(2) + (y - cy).powi(2) + (z - cz).powi(2) <= r * r,
            Shape3D::Cube => {
                let (xx, zz) = (x - cx, z - cz);
                let xb = xx * ca - zz * sa;
                let zb = xx * sa + zz * ca;
                xb.abs() <= r && (y - cy).abs() <= r && zb.abs() <= r
            }
            Shape3D::Cylinder => (x - cx).powi(2) + (z - cz).powi(2) <= r * r,
            Shape3D::Wing => {
                if (y - cy).abs() > self.half {
                    return false;
                }
                let (xx, zz) = (x - cx, z - cz);
                let xb = xx * ca - zz * sa;
                let zb = xx * sa + zz * ca;
                // Origine au quart de corde.
                let xc = xb / lf + 0.25;
                if !(0.0..=1.0).contains(&xc) {
                    return false;
                }
                (zb / lf).abs() <= naca00(0.12, xc)
            }
        }
    }
}
