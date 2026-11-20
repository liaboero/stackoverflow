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
    /// Maillage STL.
    Mesh,
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
            Shape3D::Mesh => "Maillage STL",
        }
    }

    /// Longueur de référence.
    pub fn reference(self) -> &'static str {
        match self {
            Shape3D::Sphere => "diamètre D",
            Shape3D::Cube => "arête a",
            Shape3D::Cylinder => "diamètre D (envergure périodique 4D)",
            Shape3D::Wing => "corde c (NACA 0012, allongement 4)",
            Shape3D::Mesh => "références de la mise en place",
        }
    }

    /// Dimensions de la veine en longueurs de référence.
    pub fn domain(self) -> [f64; 3] {
        match self {
            Shape3D::Sphere | Shape3D::Cube | Shape3D::Mesh => [12.0, 6.0, 6.0],
            Shape3D::Cylinder => [16.0, 4.0, 8.0],
            Shape3D::Wing => [9.0, 8.0, 4.0],
        }
    }

    /// Position de la maquette en longueurs de référence.
    fn at(self) -> [f64; 3] {
        match self {
            Shape3D::Sphere | Shape3D::Cube | Shape3D::Mesh => [4.0, 3.0, 3.0],
            Shape3D::Cylinder => [5.0, 2.0, 4.0],
            Shape3D::Wing => [2.5, 4.0, 2.0],
        }
    }

    /// L'incidence est-elle réglable ?
    pub fn has_alpha(self) -> bool {
        matches!(self, Shape3D::Cube | Shape3D::Wing | Shape3D::Mesh)
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
    /// Maillage placé dans la veine (forme `Mesh`).
    pub mesh: Option<Arc<MeshBody>>,
    /// Volume du maillage (mailles³) ; 0 pour les formes analytiques.
    pub mesh_volume: f64,
    /// Sol mobile sous la maquette (plancher de la veine en z = −½).
    pub ground: bool,
    /// Arêtes vives du maillage (coordonnées réseau), si la mise en place les
    /// demande (`MeshSetup::mesh_edges`).
    pub features: Option<Arc<crate::geom3d::FeatureSet>>,
    /// Paroi volumique par surfels (Chen et al. 1998, `surfel3d`) au lieu
    /// du rebond interpolé de Bouzidi ; maquette maillée seulement. Gardée
    /// par les grilles fines (`refined`).
    pub surfels: bool,
}

use std::sync::Arc;

use rayon::prelude::*;

use crate::stl::{MeshBody, MeshInfo, MeshSetup, Tri};

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn naca00(t: f64, x: f64) -> f64 {
    5.0 * t * (0.2969 * x.sqrt() - 0.126 * x - 0.3516 * x * x + 0.2843 * x.powi(3) - 0.1036 * x.powi(4))
}

/// Décalage en x (entre −0,5 et 0,5 maille) qui place le point le plus en
/// aval `x_te` à 0,7 maille après un nœud.
pub fn te_snap(x_te: f64) -> f64 {
    let d = 0.7 - x_te.rem_euclid(1.0);
    d - d.round()
}

impl Body3D {
    /// `alpha_deg` : incidence (cube, aile), cabrage bord de fuite vers le bas.
    pub fn new(shape: Shape3D, l: usize, alpha_deg: f64) -> Body3D {
        let [nx, ny, nz] = shape.dims(l);
        let at = shape.at();
        let lf = l as f64;
        // Centre décalé d'un quart de maille : évite les symétries d'arrondi exactes.
        let mut c = [at[0] * lf + 0.25, at[1] * lf + 0.25, at[2] * lf + 0.25];
        if shape == Shape3D::Wing {
            // Bord de fuite calé comme pour les maillages (voir `te_snap`).
            c[0] += te_snap(c[0] + 0.75 * lf * alpha_deg.to_radians().cos());
        }
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
            Shape3D::Mesh => panic!("maillage : Body3D::from_mesh"),
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
            mesh: None,
            mesh_volume: 0.0,
            ground: false,
            features: None,
            surfels: false,
        }
    }

    /// Maquette maillée placée selon `setup` : axes du fichier ramenés au
    /// repère de la veine (x vers l'aval, z vers le haut), longueur de
    /// référence ramenée à `l` mailles, incidence puis dérapage autour du
    /// centre de la boîte englobante. Veine : 3,5 M en amont, 7,5 M en aval,
    /// 2,5 M de part et d'autre (M : la plus grande de L et des dimensions de
    /// la maquette dans la direction concernée), au moins 6 L de large et de haut.
    pub fn from_mesh(tris: &[Tri], l: usize, setup: &MeshSetup) -> Result<(Body3D, MeshInfo), String> {
        let (tris, diag) = crate::stl::clean(tris.to_vec());
        if tris.is_empty() {
            return Err("aucun triangle valable dans le maillage".into());
        }
        // Arêtes vives (coordonnées du fichier), avant le raffinement lisse
        // qui les garde droites.
        let edges: Option<Vec<[[f64; 3]; 2]>> = setup.mesh_edges.then(|| match &setup.edges {
            Some(e) => e.as_ref().clone(),
            None => crate::geom3d::feature_edges(&crate::geom3d::Indexed::from_tris(&tris), crate::geom3d::FEATURE_DEG)
                .iter()
                .map(|e| [e.a, e.b])
                .collect(),
        });
        let tris = if setup.smooth { crate::geom3d::smooth_tris(&tris).0 } else { tris };
        // Volume enclos, après le raffinement lisse s'il y a lieu.
        let volume = if setup.smooth { crate::geom3d::Indexed::from_tris(&tris).signed_volume().abs() } else { diag.volume };
        let ef = setup.flow.vec();
        let ez = setup.up.vec();
        if dot3(ef, ez) != 0.0 {
            return Err("l'axe vers le haut doit être perpendiculaire à l'écoulement".into());
        }
        let ey = [
            ez[1] * ef[2] - ez[2] * ef[1],
            ez[2] * ef[0] - ez[0] * ef[2],
            ez[0] * ef[1] - ez[1] * ef[0],
        ];
        let body: Vec<Tri> = tris
            .iter()
            .map(|t| t.map(|p| [dot3(p, ef), dot3(p, ey), dot3(p, ez)]))
            .collect();
        let (lo, hi) = crate::stl::bbox(&body);
        let length = hi[0] - lo[0];
        let lref = setup.ref_length.unwrap_or(length);
        if !(length > 0.0) || !(lref > 0.0) {
            return Err("maillage sans longueur selon l'écoulement".into());
        }
        let lf = l as f64;
        let s = lf / lref;
        let c = [0, 1, 2].map(|a| 0.5 * (lo[a] + hi[a]));
        let alpha = setup.alpha_deg.to_radians();
        let (sa, ca) = alpha.sin_cos();
        let (sb, cb) = setup.beta_deg.to_radians().sin_cos();
        let place = |p: [f64; 3]| {
            let (xb, yb, zb) = (s * (p[0] - c[0]), s * (p[1] - c[1]), s * (p[2] - c[2]));
            // Incidence (avant vers le haut), puis dérapage (nez vers −y).
            let (x1, z1) = (xb * ca + zb * sa, -xb * sa + zb * ca);
            [x1 * cb + yb * sb, -x1 * sb + yb * cb, z1]
        };
        let mut placed: Vec<Tri> = body.iter().map(|t| t.map(place)).collect();
        let (rl, rh) = crate::stl::bbox(&placed);
        let e = [rh[0] - rl[0], rh[1] - rl[1], rh[2] - rl[2]];
        let k = setup.margin.clamp(0.2, 2.0);
        let mx = k * e[0].max(lf);
        let ms = k * e[1].max(e[2]).max(lf);
        let nx = (3.5 * mx + e[0] + 7.5 * mx).round() as usize;
        let ny = (6.0 * k * lf).max(e[1] + 5.0 * ms).round() as usize;
        // Avec sol : maquette posée à `ground` (unités du fichier) au-dessus du
        // plancher, 2,5 M d'air au-dessus.
        let clear = setup.ground.map(|g| g * s);
        let nz = match clear {
            Some(c) => (3.0 * k * lf).max(c + e[2] + 2.5 * ms).round() as usize,
            None => (6.0 * k * lf).max(e[2] + 5.0 * ms).round() as usize,
        };
        let cells = nx as f64 * ny as f64 * nz as f64;
        if cells > 4e8 {
            return Err(format!(
                "veine de {:.0} M mailles : réduire la résolution (L = {l})",
                cells / 1e6
            ));
        }
        // Translation : avant de la maquette à 3,5 M, centrée en y et z
        // (décalage d'un quart de maille contre les symétries d'arrondi).
        // En x, l'arrière de la maquette est calé à 0,7 maille après un nœud :
        // la portance d'une aile à bord de fuite droit dépend de la place du
        // bord de fuite dans la maille (marche d'une maille ≈ volet Gurney),
        // et 0,7 reste loin du cas défavorable (bord de fuite juste avant un
        // nœud) sur la grille de base comme sur les grilles fines (×2 : 0,4 ;
        // ×4 : 0,8 ; ×8 : 0,6).
        let sx = 3.5 * mx - rl[0];
        let shift = [
            sx + te_snap(rh[0] + sx) + setup.offset[0],
            0.5 * ny as f64 - 0.5 * (rl[1] + rh[1]) + 0.25 + setup.offset[1],
            match clear {
                Some(c) => c - 0.5 - rl[2] + 0.25,
                None => 0.5 * nz as f64 - 0.5 * (rl[2] + rh[2]) + 0.25 + setup.offset[2],
            },
        ];
        for t in &mut placed {
            for p in t.iter_mut() {
                for a in 0..3 {
                    p[a] += shift[a];
                }
            }
        }
        let features = edges.map(|es| {
            let to = |p: [f64; 3]| {
                let q = place([dot3(p, ef), dot3(p, ey), dot3(p, ez)]);
                [q[0] + shift[0], q[1] + shift[1], q[2] + shift[2]]
            };
            let fe = es
                .iter()
                .map(|[a, b]| crate::geom3d::FeatureEdge { a: to(*a), b: to(*b), angle_deg: 180.0, convex: true })
                .collect();
            Arc::new(crate::geom3d::FeatureSet::new(fe))
        });
        let mesh = MeshBody::new(placed, !diag.watertight());
        let bounds = Bounds {
            x0: mesh.lo[0],
            x1: mesh.hi[0],
            y0: mesh.lo[1],
            y1: mesh.hi[1],
            z0: mesh.lo[2],
            z1: mesh.hi[2],
        };
        let center = [0, 1, 2].map(|a| 0.5 * (mesh.lo[a] + mesh.hi[a]));
        let frontal = mesh.frontal_area();
        if frontal <= 0.0 {
            return Err("maître-couple nul".into());
        }
        let area = setup.ref_area.map(|a| a * s * s).unwrap_or(frontal);
        // Mètres par maille.
        let h = setup.unit.meters() / s;
        let info = MeshInfo {
            triangles: tris.len() + diag.degenerate + diag.duplicates,
            size_m: e.map(|v| v * h),
            ref_length_m: lref * setup.unit.meters(),
            ref_area_m2: area * h * h,
            frontal_m2: frontal * h * h,
            blockage: frontal / (ny as f64 * nz as f64),
            volume_m3: volume * setup.unit.meters().powi(3),
            diag: diag.clone(),
        };
        let b = Body3D {
            shape: Shape3D::Mesh,
            l,
            nx,
            ny,
            nz,
            center,
            bounds,
            area,
            alpha,
            half: 0.0,
            mesh: Some(Arc::new(mesh)),
            mesh_volume: volume * s * s * s,
            ground: setup.ground.is_some(),
            features,
            surfels: false,
        };
        Ok((b, info))
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
            Shape3D::Mesh => self.mesh.as_ref().is_some_and(|m| m.inside(x, y, z)),
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

/// Distance signée à la maquette (mailles, négative dedans), échantillonnée
/// tous les `h` mailles : sert seulement à l'affichage d'une surface lisse.
pub struct SdfGrid {
    pub origin: [f64; 3],
    pub h: f64,
    pub dims: [usize; 3],
    pub d: Vec<f32>,
}

/// Distance d'un point à une ligne polygonale fermée.
fn poly_dist(p: [f64; 2], poly: &[[f64; 2]]) -> f64 {
    let mut best = f64::INFINITY;
    for k in 0..poly.len() {
        let a = poly[k];
        let b = poly[(k + 1) % poly.len()];
        let (ex, ey) = (b[0] - a[0], b[1] - a[1]);
        let (wx, wy) = (p[0] - a[0], p[1] - a[1]);
        let t = ((wx * ex + wy * ey) / (ex * ex + ey * ey).max(1e-30)).clamp(0.0, 1.0);
        best = best.min((wx - t * ex).hypot(wy - t * ey));
    }
    best
}

impl Body3D {
    /// La même maquette dans une grille deux fois plus fine dont le nœud 0 est
    /// au nœud `origin` de celle-ci, de dimensions `dims`.
    pub fn refined(&self, origin: [f64; 3], dims: [usize; 3]) -> Body3D {
        let t = |p: [f64; 3]| [0, 1, 2].map(|a| 2.0 * (p[a] - origin[a]));
        let b = self.bounds;
        let [x0, y0, z0] = t([b.x0, b.y0, b.z0]);
        let [x1, y1, z1] = t([b.x1, b.y1, b.z1]);
        Body3D {
            shape: self.shape,
            l: 2 * self.l,
            nx: dims[0],
            ny: dims[1],
            nz: dims[2],
            center: t(self.center),
            bounds: Bounds { x0, x1, y0, y1, z0, z1 },
            area: 4.0 * self.area,
            alpha: self.alpha,
            half: 2.0 * self.half,
            mesh: self.mesh.as_ref().map(|m| Arc::new(m.transformed(|p| t(p)))),
            mesh_volume: 8.0 * self.mesh_volume,
            ground: false,
            features: self.features.as_ref().map(|f| Arc::new(f.transformed(|p| t(p)))),
            surfels: self.surfels,
        }
    }

    /// Profil NACA 0012 de l'aile dans son repère (mailles, origine au quart de corde).
    fn foil(&self) -> Vec<[f64; 2]> {
        let lf = self.l as f64;
        let n = 160;
        let xs: Vec<f64> = (0..=n)
            .map(|k| 0.5 * (1.0 - (k as f64 / n as f64 * std::f64::consts::PI).cos()))
            .collect();
        let mut v: Vec<[f64; 2]> = xs.iter().map(|&x| [(x - 0.25) * lf, naca00(0.12, x) * lf]).collect();
        v.extend(xs.iter().rev().skip(1).map(|&x| [(x - 0.25) * lf, -naca00(0.12, x) * lf]));
        v
    }

    /// Distance à la surface (mailles, sans signe), plafonnée à `cap`.
    /// `foil` : profil de l'aile (`Body3D::foil`), ignoré pour les autres formes.
    fn surface_distance(&self, p: [f64; 3], cap: f64, foil: &[[f64; 2]]) -> f64 {
        if let Some(m) = &self.mesh {
            return m.distance(p, cap);
        }
        let [cx, cy, cz] = self.center;
        let (x, y, z) = (p[0], p[1], p[2]);
        let lf = self.l as f64;
        let r = lf / 2.0;
        let (sa, ca) = self.alpha.sin_cos();
        let d = match self.shape {
            Shape3D::Sphere => ((x - cx).powi(2) + (y - cy).powi(2) + (z - cz).powi(2)).sqrt() - r,
            Shape3D::Cylinder => (x - cx).hypot(z - cz) - r,
            Shape3D::Cube => {
                let (xx, zz) = (x - cx, z - cz);
                let q = [(xx * ca - zz * sa).abs() - r, (y - cy).abs() - r, (xx * sa + zz * ca).abs() - r];
                let out = q.iter().map(|v| v.max(0.0).powi(2)).sum::<f64>().sqrt();
                out + q[0].max(q[1]).max(q[2]).min(0.0)
            }
            Shape3D::Wing => {
                let (xx, zz) = (x - cx, z - cz);
                let xb = xx * ca - zz * sa;
                let zb = xx * sa + zz * ca;
                let sec = poly_dist([xb, zb], foil);
                let w = [sec, (y - cy).abs() - self.half];
                w[0].max(w[1]).min(0.0) + w[0].max(0.0).hypot(w[1].max(0.0))
            }
            Shape3D::Mesh => cap,
        };
        d.abs().min(cap)
    }

    /// Distance à la paroi de chaque maille (mailles) pour Spalart–Allmaras :
    /// exacte à moins de 3 mailles du corps, propagée au-delà par une
    /// transformée de chanfrein (1, √2, √3) dans une boîte de 1,5 L autour de
    /// la maquette ; 10⁴ plus loin, 0 dans les mailles solides.
    pub fn wall_distance(&self, solid: &[bool]) -> Vec<f32> {
        let (nx, ny, nz) = (self.nx, self.ny, self.nz);
        let foil = if self.shape == Shape3D::Wing { self.foil() } else { Vec::new() };
        let b = self.bounds;
        let r = 1.5 * self.l as f64;
        let lo = [
            (b.x0 - r).floor().max(0.0) as usize,
            (b.y0 - r).floor().max(0.0) as usize,
            (b.z0 - r).floor().max(0.0) as usize,
        ];
        let hi = [
            ((b.x1 + r).ceil() as usize).min(nx - 1),
            ((b.y1 + r).ceil() as usize).min(ny - 1),
            ((b.z1 + r).ceil() as usize).min(nz - 1),
        ];
        let far = 1e4f32;
        let idx = |x: usize, y: usize, z: usize| x + nx * (y + ny * z);
        let mut d = vec![far; nx * ny * nz];
        // Graines exactes près du corps (en parallèle, tranche par tranche).
        let near = [
            (b.x0 - 3.0).floor().max(0.0) as usize,
            (b.y0 - 3.0).floor().max(0.0) as usize,
            (b.z0 - 3.0).floor().max(0.0) as usize,
        ];
        let near_hi = [
            ((b.x1 + 3.0).ceil() as usize).min(nx - 1),
            ((b.y1 + 3.0).ceil() as usize).min(ny - 1),
            ((b.z1 + 3.0).ceil() as usize).min(nz - 1),
        ];
        let seeds: Vec<(usize, f32)> = (near[2]..=near_hi[2])
            .into_par_iter()
            .flat_map_iter(|z| {
                let foil = &foil;
                (near[1]..=near_hi[1]).flat_map(move |y| {
                    (near[0]..=near_hi[0]).filter_map(move |x| {
                        let c = idx(x, y, z);
                        if solid[c] {
                            return Some((c, 0.0));
                        }
                        let v = self.surface_distance([x as f64, y as f64, z as f64], 3.0, foil);
                        (v < 3.0).then_some((c, v as f32))
                    })
                })
            })
            .collect();
        for (c, v) in seeds {
            d[c] = v;
        }
        // Chanfrein : passe avant puis arrière dans la boîte.
        let (s2, s3) = (2f32.sqrt(), 3f32.sqrt());
        let nb: Vec<([i64; 3], f32)> = (-1..=1i64)
            .flat_map(|k| (-1..=1i64).flat_map(move |j| (-1..=1i64).map(move |i| [i, j, k])))
            .filter(|o| *o != [0, 0, 0])
            .map(|o| {
                let m = o.iter().map(|v| v.abs()).sum::<i64>();
                (o, [0.0, 1.0, s2, s3][m as usize])
            })
            .collect();
        let before = |o: &[i64; 3]| o[2] < 0 || (o[2] == 0 && (o[1] < 0 || (o[1] == 0 && o[0] < 0)));
        let fwd: Vec<_> = nb.iter().filter(|(o, _)| before(o)).cloned().collect();
        let bwd: Vec<_> = nb.iter().filter(|(o, _)| !before(o)).cloned().collect();
        let pass = |d: &mut Vec<f32>, offs: &[([i64; 3], f32)], rev: bool| {
            let zs: Vec<usize> = if rev { (lo[2]..=hi[2]).rev().collect() } else { (lo[2]..=hi[2]).collect() };
            for z in zs {
                let ys: Vec<usize> = if rev { (lo[1]..=hi[1]).rev().collect() } else { (lo[1]..=hi[1]).collect() };
                for y in ys {
                    let xs: Vec<usize> = if rev { (lo[0]..=hi[0]).rev().collect() } else { (lo[0]..=hi[0]).collect() };
                    for x in xs {
                        let c = idx(x, y, z);
                        let mut best = d[c];
                        for (o, w) in offs {
                            let (px, py, pz) = (x as i64 + o[0], y as i64 + o[1], z as i64 + o[2]);
                            if px < lo[0] as i64 || py < lo[1] as i64 || pz < lo[2] as i64 {
                                continue;
                            }
                            if px > hi[0] as i64 || py > hi[1] as i64 || pz > hi[2] as i64 {
                                continue;
                            }
                            best = best.min(d[idx(px as usize, py as usize, pz as usize)] + w);
                        }
                        d[c] = best;
                    }
                }
            }
        };
        pass(&mut d, &fwd, false);
        pass(&mut d, &bwd, true);
        d
    }

    /// Normale sortante de la paroi et distance à celle-ci au nœud p (mailles) :
    /// gradient de la distance exacte, pour la loi de paroi.
    pub fn wall_normal(&self, p: [f64; 3]) -> ([f64; 3], f64) {
        let foil = if self.shape == Shape3D::Wing { self.foil() } else { Vec::new() };
        let dist = |q: [f64; 3]| self.surface_distance(q, 4.0, &foil);
        // Gradient de la distance signée (négative dans le corps) : à moins
        // de e de la paroi, les différences centrées de la distance non
        // signée enjambent la surface et la composante normale du gradient
        // s'annule (normale fausse, relevé de la loi de paroi dans le corps).
        let sdist = |q: [f64; 3]| if self.inside(q[0], q[1], q[2]) { -dist(q) } else { dist(q) };
        let e = 0.2;
        let mut g = [0.0; 3];
        for a in 0..3 {
            let (mut pp, mut pm) = (p, p);
            pp[a] += e;
            pm[a] -= e;
            g[a] = sdist(pp) - sdist(pm);
        }
        let n = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
        let g = if n > 1e-9 { g.map(|v| v / n) } else { [0.0, 0.0, 1.0] };
        (g, dist(p))
    }

    /// Distance signée, `h` mailles entre échantillons, autour de la maquette.
    pub fn sdf_grid(&self, h: f64) -> SdfGrid {
        let b = self.bounds;
        let m = 2.0;
        let lo = [(b.x0 - m).max(0.0), (b.y0 - m).max(0.0), (b.z0 - m).max(0.0)];
        let hi = [
            (b.x1 + m).min(self.nx as f64 - 1.0),
            (b.y1 + m).min(self.ny as f64 - 1.0),
            (b.z1 + m).min(self.nz as f64 - 1.0),
        ];
        let dims = [0, 1, 2].map(|a| ((hi[a] - lo[a]) / h).ceil() as usize + 1);
        let [cx, cy, cz] = self.center;
        let lf = self.l as f64;
        let r = lf / 2.0;
        let (sa, ca) = self.alpha.sin_cos();
        // Profil de l'aile dans son repère (mailles, origine au quart de corde).
        let foil: Vec<[f64; 2]> = if self.shape == Shape3D::Wing {
            let n = 160;
            let xs: Vec<f64> = (0..=n)
                .map(|k| 0.5 * (1.0 - (k as f64 / n as f64 * std::f64::consts::PI).cos()))
                .collect();
            let mut v: Vec<[f64; 2]> = xs.iter().map(|&x| [(x - 0.25) * lf, naca00(0.12, x) * lf]).collect();
            v.extend(xs.iter().rev().skip(1).map(|&x| [(x - 0.25) * lf, -naca00(0.12, x) * lf]));
            v
        } else {
            Vec::new()
        };
        // Section (x, z) de l'aile, indépendante de y : calculée une fois.
        let section: Vec<f64> = if self.shape == Shape3D::Wing {
            let mut s = vec![0.0; dims[0] * dims[2]];
            for k in 0..dims[2] {
                for i in 0..dims[0] {
                    let (xx, zz) = (lo[0] + i as f64 * h - cx, lo[2] + k as f64 * h - cz);
                    let xb = xx * ca - zz * sa;
                    let zb = xx * sa + zz * ca;
                    let xc = xb / lf + 0.25;
                    let inside = (0.0..=1.0).contains(&xc) && (zb / lf).abs() <= naca00(0.12, xc);
                    let d = poly_dist([xb, zb], &foil);
                    s[k * dims[0] + i] = if inside { -d } else { d };
                }
            }
            s
        } else {
            Vec::new()
        };
        let mut d = vec![0f32; dims[0] * dims[1] * dims[2]];
        if let Some(m) = &self.mesh {
            // Distance aux triangles (plafonnée à une maille : seule la surface compte pour l'affichage), signe par appartenance.
            d.par_iter_mut().enumerate().for_each(|(n, v)| {
                let i = n % dims[0];
                let j = (n / dims[0]) % dims[1];
                let k = n / (dims[0] * dims[1]);
                let p = [lo[0] + i as f64 * h, lo[1] + j as f64 * h, lo[2] + k as f64 * h];
                let dist = m.distance(p, 1.0);
                *v = if m.inside(p[0], p[1], p[2]) { -dist } else { dist } as f32;
            });
            return SdfGrid { origin: lo, h, dims, d };
        }
        for k in 0..dims[2] {
            for j in 0..dims[1] {
                for i in 0..dims[0] {
                    let (x, y, z) = (lo[0] + i as f64 * h, lo[1] + j as f64 * h, lo[2] + k as f64 * h);
                    let v = match self.shape {
                        Shape3D::Sphere => ((x - cx).powi(2) + (y - cy).powi(2) + (z - cz).powi(2)).sqrt() - r,
                        Shape3D::Cylinder => (x - cx).hypot(z - cz) - r,
                        Shape3D::Cube => {
                            let (xx, zz) = (x - cx, z - cz);
                            let q = [(xx * ca - zz * sa).abs() - r, (y - cy).abs() - r, (xx * sa + zz * ca).abs() - r];
                            let out = q.iter().map(|v| v.max(0.0).powi(2)).sum::<f64>().sqrt();
                            out + q[0].max(q[1]).max(q[2]).min(0.0)
                        }
                        Shape3D::Mesh => 0.0,
                        Shape3D::Wing => {
                            // Extrusion de la section sur l'envergure.
                            let w = [section[k * dims[0] + i], (y - cy).abs() - self.half];
                            w[0].max(w[1]).min(0.0) + w[0].max(0.0).hypot(w[1].max(0.0))
                        }
                    };
                    d[i + dims[0] * (j + dims[1] * k)] = v as f32;
                }
            }
        }
        SdfGrid { origin: lo, h, dims, d }
    }
}

/// Corps d'Ahmed (Ahmed, Ramm et Faltin 1984, SAE 840300) : dimensions (mm).
pub mod ahmed {
    pub const LONGUEUR: f64 = 1044.0;
    pub const LARGEUR: f64 = 389.0;
    pub const HAUTEUR: f64 = 288.0;
    /// Rayon des arrondis de l'avant.
    pub const RAYON: f64 = 100.0;
    /// Longueur de la lunette, mesurée le long de la pente.
    pub const LUNETTE: f64 = 222.0;
    /// Garde au sol.
    pub const GARDE: f64 = 50.0;
    pub const PIED_DIAMETRE: f64 = 30.0;
    /// Centres des quatre pieds (x depuis l'avant, y depuis le plan de
    /// symétrie). En x : 202 et 672 mm (entraxe 470 mm, dessin d'Ahmed) ; en
    /// y : ±163,5 mm, valeur reprise de la littérature et non vérifiée à la
    /// source (Ahmed corrige sa traînée de celle des pieds).
    pub const PIEDS: [[f64; 2]; 4] = [[202.0, -163.5], [202.0, 163.5], [672.0, -163.5], [672.0, 163.5]];
    /// Exposant de la super-ellipse des coins avant à mi-congé (x = 50 mm),
    /// ajusté sur la géométrie d'origine (fichier `Ahmed_Front_Geo.dat`
    /// d'ERCOFTAC, données de S. R. Ahmed) : les coins s'arrondissent à
    /// l'avant et deviennent vifs à x = 100 mm.
    pub const EXPOSANT_COIN: f64 = 8.0;
}

/// Nœuds de [a, b] espacés d'au plus `ds`, avec les valeurs imposées
/// `forced` ; un nœud régulier trop proche d'un nœud imposé est retiré.
fn ahmed_nodes(a: f64, b: f64, forced: &[f64], ds: f64) -> Vec<f64> {
    let n = ((b - a) / ds).ceil().max(1.0) as usize;
    let mut f: Vec<f64> = forced.iter().copied().filter(|&v| v > a && v < b).collect();
    f.push(a);
    f.push(b);
    let mut out: Vec<f64> = (1..n)
        .map(|k| a + (b - a) * k as f64 / n as f64)
        .filter(|&v| f.iter().all(|&g| (v - g).abs() > 0.3 * ds))
        .collect();
    out.extend(f);
    out.sort_by(|p, q| p.partial_cmp(q).unwrap());
    out.dedup_by(|p, q| (*p - *q).abs() < 1e-9);
    out
}

/// Corps d'Ahmed maillé (mm) : x de l'avant (0) vers l'arrière, y latéral,
/// z vers le haut depuis le dessous de la caisse (le sol est à z = −50).
/// Lunette de `slant_deg` degrés ; avec `stilts`, quatre pieds de 30 mm qui
/// descendent jusqu'à z = −50 − `depth` (ils traversent le plancher de
/// `depth` pour ne pas laisser de fente sous eux). Facettes d'au plus `ds`
/// environ : l'intégration des efforts sur la surface relève la pression au
/// centre de chaque triangle. Maillage fermé, normales sortantes.
///
/// Avant : congés de 100 mm sur les quatre arêtes de la face avant ; la
/// section à l'abscisse x < 100 est le rectangle intérieur (|y| ≤ 94,5,
/// 100 ≤ z ≤ 188) élargi de r = √(100² − (100 − x)²) avec des coins en
/// super-ellipse d'exposant 4 / w, w = (100 − x)/100 (8 à mi-congé, vifs à
/// x = 100), ce qui reproduit la géométrie d'origine à quelques millimètres
/// près. Arêtes longitudinales, de lunette et de culot vives.
pub fn ahmed_tris(slant_deg: f64, stilts: bool, ds: f64, depth: f64) -> Vec<Tri> {
    use ahmed::*;
    let (l, hw, h, r0) = (LONGUEUR, LARGEUR / 2.0, HAUTEUR, RAYON);
    let (cy, zlo, zhi) = (hw - r0, r0, h - r0);
    let phi = slant_deg.to_radians();
    let xs = l - LUNETTE * phi.cos();
    let ztop = |x: f64| if x <= xs { h } else { h - (x - xs) * phi.tan() };
    let a = 25.0; // demi-côté de la fenêtre carrée autour de chaque pied
    let rp = PIED_DIAMETRE / 2.0;
    let pieds: Vec<[f64; 2]> = if stilts { PIEDS.to_vec() } else { Vec::new() };
    // Pourtour de la section (rectangle hw × h) : dessous de −hw à hw, côté
    // droit, dessus, côté gauche.
    let mut fy = vec![-cy, cy];
    for p in &pieds {
        fy.extend([p[1] - a, p[1] + a]);
    }
    let ys = ahmed_nodes(-hw, hw, &fy, ds);
    let zs = ahmed_nodes(0.0, h, &[zlo, zhi], ds);
    let mut perim: Vec<[f64; 2]> = ys.iter().map(|&y| [y, 0.0]).collect();
    let nbot = perim.len();
    perim.extend(zs[1..].iter().map(|&z| [hw, z]));
    perim.extend(ys.iter().rev().skip(1).map(|&y| [y, h]));
    perim.extend(zs.iter().rev().skip(1).take(zs.len() - 2).map(|&z| [-hw, z]));
    let np = perim.len();
    // Sections : congé avant (pas angulaire), puis x ≥ 100.
    let nf = (std::f64::consts::FRAC_PI_2 * r0 / ds).ceil() as usize;
    let mut stations: Vec<f64> = (0..nf).map(|j| r0 * (1.0 - (std::f64::consts::FRAC_PI_2 * j as f64 / nf as f64).cos())).collect();
    let mut fx = vec![xs];
    for p in &pieds {
        fx.extend([p[0] - a, p[0] + a]);
    }
    stations.extend(ahmed_nodes(r0, l, &fx, ds));
    let section = |x: f64| -> Vec<[f64; 3]> {
        if x >= r0 {
            let k = ztop(x) / h;
            return perim.iter().map(|&[y, z]| [x, y, z * k]).collect();
        }
        let w = (r0 - x) / r0;
        let r = r0 * (1.0 - w * w).max(0.0).sqrt();
        let p = 0.5 * EXPOSANT_COIN / w.max(1e-9);
        perim
            .iter()
            .map(|&[y, z]| {
                let dy = (y.abs() - cy).max(0.0);
                let dz = if z < zlo { zlo - z } else if z > zhi { z - zhi } else { 0.0 };
                let m = dy.max(dz);
                let norm = if p > 500.0 { m } else { m * ((dy / m).powf(p) + (dz / m).powf(p)).powf(1.0 / p) };
                let s = if norm > 0.0 { r / norm } else { 0.0 };
                let yy = if dy > 0.0 { y.signum() * (cy + dy * s) } else { y };
                let zz = if z < zlo { zlo - dz * s } else if z > zhi { zhi + dz * s } else { z };
                [x, yy, zz]
            })
            .collect()
    };
    let rings: Vec<Vec<[f64; 3]>> = stations.iter().map(|&x| section(x)).collect();
    let mut out: Vec<Tri> = Vec::new();
    // Triangle orienté selon `dir` (normale sortante approchée).
    let mut push = |t: Tri, dir: [f64; 3]| {
        let e1 = [t[1][0] - t[0][0], t[1][1] - t[0][1], t[1][2] - t[0][2]];
        let e2 = [t[2][0] - t[0][0], t[2][1] - t[0][1], t[2][2] - t[0][2]];
        let n = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]];
        out.push(if dot3(n, dir) >= 0.0 { t } else { [t[0], t[2], t[1]] });
    };
    let quad = |a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]| [[a, b, c], [a, c, d]];
    let hole = |x0: f64, x1: f64, y0: f64, y1: f64| {
        pieds.iter().any(|p| x0 >= p[0] - a - 1e-9 && x1 <= p[0] + a + 1e-9 && y0 >= p[1] - a - 1e-9 && y1 <= p[1] + a + 1e-9)
    };
    // Bandes entre sections.
    for i in 0..rings.len() - 1 {
        let (ra, rb) = (&rings[i], &rings[i + 1]);
        let (xa, xb) = (stations[i], stations[i + 1]);
        let xm = 0.5 * (xa + xb);
        let zc = if xm < r0 { 0.5 * h } else { 0.5 * ztop(xm) };
        for j in 0..np {
            let k = (j + 1) % np;
            if j + 1 < nbot && xa >= r0 && hole(xa, xb, ys[j], ys[j + 1]) {
                continue;
            }
            for t in quad(ra[j], ra[k], rb[k], rb[j]) {
                let c = [0, 1, 2].map(|q| (t[0][q] + t[1][q] + t[2][q]) / 3.0);
                push(t, [-(r0 - xm).max(0.0), c[1], c[2] - zc]);
            }
        }
    }
    // Bouchons plans par anneaux homothétiques (face avant, culot).
    let cap = |ring: &[[f64; 3]], dir: [f64; 3], push: &mut dyn FnMut(Tri, [f64; 3])| {
        let n = ring.len() as f64;
        let c = [0, 1, 2].map(|q| ring.iter().map(|p| p[q]).sum::<f64>() / n);
        let rmax = ring.iter().map(|p| ((p[1] - c[1]).powi(2) + (p[2] - c[2]).powi(2)).sqrt()).fold(0.0, f64::max);
        let nk = (rmax / ds).ceil().max(1.0) as usize;
        let at = |p: [f64; 3], s: f64| [0, 1, 2].map(|q| c[q] + s * (p[q] - c[q]));
        for k in 0..nk {
            let (s0, s1) = (k as f64 / nk as f64, (k + 1) as f64 / nk as f64);
            for j in 0..ring.len() {
                let jn = (j + 1) % ring.len();
                for t in quad(at(ring[j], s0), at(ring[jn], s0), at(ring[jn], s1), at(ring[j], s1)) {
                    push(t, dir);
                }
            }
        }
    };
    cap(&rings[0], [-1.0, 0.0, 0.0], &mut push);
    cap(&rings[rings.len() - 1], [1.0, 0.0, 0.0], &mut push);
    // Pieds : couronne carré → cercle dans le dessous, fût, fond.
    for p in &pieds {
        let ia: Vec<usize> = (0..stations.len()).filter(|&i| stations[i] >= p[0] - a - 1e-9 && stations[i] <= p[0] + a + 1e-9).collect();
        let ja: Vec<usize> = (0..nbot).filter(|&j| ys[j] >= p[1] - a - 1e-9 && ys[j] <= p[1] + a + 1e-9).collect();
        let (y0, y1) = (ys[ja[0]], ys[*ja.last().unwrap()]);
        let (x0, x1) = (stations[ia[0]], stations[*ia.last().unwrap()]);
        let mut sq: Vec<[f64; 3]> = Vec::new();
        sq.extend(ia.iter().map(|&i| [stations[i], y0, 0.0]));
        sq.extend(ja.iter().skip(1).map(|&j| [x1, ys[j], 0.0]));
        sq.extend(ia.iter().rev().skip(1).map(|&i| [stations[i], y1, 0.0]));
        sq.extend(ja.iter().rev().skip(1).take(ja.len() - 2).map(|&j| [x0, ys[j], 0.0]));
        let circ = |q: [f64; 3], z: f64| {
            let (dx, dy) = (q[0] - p[0], q[1] - p[1]);
            let d = dx.hypot(dy);
            [p[0] + rp * dx / d, p[1] + rp * dy / d, z]
        };
        let zb = -GARDE - depth;
        let nzs = ((GARDE + depth) / ds).ceil().max(1.0) as usize;
        let n = sq.len();
        for j in 0..n {
            let k = (j + 1) % n;
            for t in quad(sq[j], sq[k], circ(sq[k], 0.0), circ(sq[j], 0.0)) {
                push(t, [0.0, 0.0, -1.0]);
            }
            for m in 0..nzs {
                let (z0, z1) = (-(GARDE + depth) * m as f64 / nzs as f64, -(GARDE + depth) * (m + 1) as f64 / nzs as f64);
                let dir = [sq[j][0] - p[0], sq[j][1] - p[1], 0.0];
                for t in quad(circ(sq[j], z0), circ(sq[k], z0), circ(sq[k], z1), circ(sq[j], z1)) {
                    push(t, dir);
                }
            }
        }
        let bottom: Vec<[f64; 3]> = sq.iter().map(|&q| circ(q, zb)).collect();
        cap(&bottom, [0.0, 0.0, -1.0], &mut push);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stl::{icosphere, Axis, MeshSetup, Unit};

    /// Ellipsoïde 100 × 25 × 18 mm (fichier : x vers l'aval, z vers le haut).
    fn fuseau() -> Vec<Tri> {
        icosphere(1.0, [0.0; 3], 3)
            .into_iter()
            .map(|t| t.map(|p| [50.0 * p[0], 12.5 * p[1], 9.0 * p[2]]))
            .collect()
    }

    #[test]
    fn ahmed_ferme() {
        use super::ahmed::*;
        for (phi, pieds) in [(25.0, true), (35.0, false)] {
            let (t, d) = crate::stl::clean(ahmed_tris(phi, pieds, 8.0, 5.0));
            assert!(d.watertight(), "{phi}° : {}", d.summary());
            assert_eq!(d.flipped_edges, 0);
            let (lo, hi) = crate::stl::bbox(&t);
            assert!((hi[0] - lo[0] - LONGUEUR).abs() < 1e-9 && (hi[1] - lo[1] - LARGEUR).abs() < 1e-9);
            assert!((hi[2] - HAUTEUR).abs() < 1e-9);
            let p = phi.to_radians();
            let v0 = LONGUEUR * LARGEUR * HAUTEUR - 0.5 * LUNETTE * LUNETTE * p.sin() * p.cos() * LARGEUR;
            // Congés de l'avant : de l'ordre de 2 % du volume.
            assert!(d.volume > 0.95 * v0 && d.volume < 0.995 * v0, "{phi}° : {} contre {v0}", d.volume);
        }
    }

    fn voxels(b: &Body3D) -> usize {
        let bb = b.bounds;
        let mut n = 0;
        for z in bb.z0.floor() as i64..=bb.z1.ceil() as i64 {
            for y in bb.y0.floor() as i64..=bb.y1.ceil() as i64 {
                for x in bb.x0.floor() as i64..=bb.x1.ceil() as i64 {
                    n += b.inside(x as f64, y as f64, z as f64) as usize;
                }
            }
        }
        n
    }

    #[test]
    fn placement_stl() {
        let l = 24;
        let base = MeshSetup {
            alpha_deg: 6.0,
            ..Default::default()
        };
        let (a, ia) = Body3D::from_mesh(&fuseau(), l, &base).unwrap();
        // Même maquette décrite avec d'autres axes : écoulement selon −y, haut selon +x.
        let turned: Vec<Tri> = fuseau().into_iter().map(|t| t.map(|p| [p[2], -p[0], -p[1]])).collect();
        let setup = MeshSetup {
            flow: Axis::Ym,
            up: Axis::Xp,
            ..base.clone()
        };
        let (b, ib) = Body3D::from_mesh(&turned, l, &setup).unwrap();
        assert_eq!([a.nx, a.ny, a.nz], [b.nx, b.ny, b.nz]);
        assert_eq!(voxels(&a), voxels(&b));
        assert!((ia.ref_length_m - 0.1).abs() < 1e-3, "longueur {}", ia.ref_length_m);
        assert!((ia.frontal_m2 - ib.frontal_m2).abs() < 1e-9);
        // Surface imposée, en cm.
        let (c, ic) = Body3D::from_mesh(
            &fuseau(),
            l,
            &MeshSetup {
                unit: Unit::Cm,
                ref_area: Some(400.0),
                ..base.clone()
            },
        )
        .unwrap();
        assert!((ic.ref_area_m2 - 0.04).abs() < 1e-9);
        assert!((c.area - 400.0 * (l as f64 / 100.0).powi(2)).abs() < 1e-6);
        // Dérapage : la maquette s'élargit dans la veine.
        let (_, id) = Body3D::from_mesh(
            &fuseau(),
            l,
            &MeshSetup {
                beta_deg: 30.0,
                ..base.clone()
            },
        )
        .unwrap();
        assert!(id.size_m[1] > 2.0 * ia.size_m[1]);
        // Axes non perpendiculaires.
        let bad = MeshSetup {
            up: Axis::Xm,
            ..base
        };
        assert!(Body3D::from_mesh(&fuseau(), l, &bad).is_err());
    }

    /// Aile NACA 0012 grossière (corde 100 mm, envergure 200 mm), arêtes vives
    /// du maillage, raffinement lisse et lecture `.sgeo`.
    #[test]
    fn pretraitement_stl() {
        let tris = crate::geom3d::tests::naca_wing(0.12, 100.0, 200.0, 16);
        let l = 32;
        let base = MeshSetup {
            alpha_deg: 4.0,
            margin: 0.4,
            ..Default::default()
        };
        let (a, _) = Body3D::from_mesh(&tris, l, &base).unwrap();
        assert!(a.features.is_none());
        // Arêtes vives : bord de fuite et saumons, placés avec la maquette.
        let (b, _) = Body3D::from_mesh(&tris, l, &MeshSetup { mesh_edges: true, ..base.clone() }).unwrap();
        let f = b.features.clone().unwrap();
        let te = [b.bounds.x1, 0.5 * (b.bounds.y0 + b.bounds.y1), 0.0];
        // Point le plus en aval (bord de fuite, sous l'incidence) : sur une arête vive.
        let mesh = b.mesh.as_ref().unwrap();
        let p = mesh.tris.iter().flatten().copied().filter(|p| (p[1] - te[1]).abs() < 4.0).max_by(|u, v| u[0].total_cmp(&v[0])).unwrap();
        assert!(f.distance(p, 2.0) < 1e-6, "bord de fuite {}", f.distance(p, 2.0));
        // Bord d'attaque : loin de toute arête vive.
        let q = mesh.tris.iter().flatten().copied().filter(|p| (p[1] - te[1]).abs() < 4.0).min_by(|u, v| u[0].total_cmp(&v[0])).unwrap();
        assert!(f.distance(q, 5.0) >= 5.0);
        // Grille fine : arêtes transformées avec la maquette.
        let r = b.refined([10.0, 10.0, 10.0], [b.nx, b.ny, b.nz]);
        let p2 = [0, 1, 2].map(|k| 2.0 * (p[k] - 10.0));
        assert!(r.features.as_ref().unwrap().distance(p2, 2.0) < 1e-6);
        // Raffinement lisse : volume plus proche du profil exact (le polygone
        // facetté est inscrit dans le profil), même mise en place.
        let (c, _) = Body3D::from_mesh(&tris, l, &MeshSetup { smooth: true, ..base.clone() }).unwrap();
        let exact = crate::geom3d::Indexed::from_tris(&crate::geom3d::tests::naca_wing(0.12, 100.0, 200.0, 400)).signed_volume();
        let s3 = (l as f64 / 100.0).powi(3);
        assert!(a.mesh_volume < c.mesh_volume, "{} {}", a.mesh_volume, c.mesh_volume);
        assert!((c.mesh_volume - exact * s3).abs() < 0.5 * (a.mesh_volume - exact * s3).abs());
        assert_eq!([a.nx, a.ny, a.nz], [c.nx, c.ny, c.nz]);
        // Fichier .sgeo : mêmes triangles et même voxelisation que le STL.
        let g = crate::sgeo::GeoFile::from_tris(&tris);
        let h = crate::sgeo::GeoFile::parse(&g.to_text()).unwrap();
        let (d, _) = Body3D::from_mesh(&h.tris(), l, &base).unwrap();
        assert_eq!(voxels(&a), voxels(&d));
        let segs = std::sync::Arc::new(h.feature_segments().iter().map(|e| [e.a, e.b]).collect());
        let (e, _) = Body3D::from_mesh(&h.tris(), l, &MeshSetup { mesh_edges: true, edges: Some(segs), ..base }).unwrap();
        assert!(e.features.unwrap().distance(p, 2.0) < 1e-6);
    }
}
