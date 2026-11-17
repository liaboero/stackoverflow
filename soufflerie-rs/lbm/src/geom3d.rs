//! Prétraitement géométrique des maquettes maillées : maillage indexé
//! (sommets soudés, adjacence des faces), arêtes vives tirées de l'angle
//! dièdre entre faces voisines, et contrôle de résolution d'une aile par
//! tranches d'envergure (rayon du bord d'attaque, épaisseur du bord de fuite,
//! longueur de paroi plus mince qu'une maille, taille des facettes).
//!
//! Toutes les grandeurs sont dans le repère où les triangles sont donnés :
//! sur une maquette placée (`MeshBody`), en mailles de la veine, x vers
//! l'aval, y en envergure, z vers le haut.

use std::collections::{HashMap, HashSet};

use crate::stl::{bbox, Tri};

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Maillage indexé : sommets soudés à 10⁻⁷ de la diagonale près (comme
/// `stl::clean`), faces par indices de sommets.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Indexed {
    pub verts: Vec<[f64; 3]>,
    pub faces: Vec<[u32; 3]>,
}

impl Indexed {
    pub fn from_tris(tris: &[Tri]) -> Indexed {
        let (lo, hi) = bbox(tris);
        let eps = 1e-7 * norm(sub(hi, lo)).max(1e-30);
        let mut ids: HashMap<[i64; 3], u32> = HashMap::new();
        let mut verts = Vec::new();
        let mut faces = Vec::with_capacity(tris.len());
        for t in tris {
            let f = t.map(|p| {
                let k = p.map(|c| (c / eps).round() as i64);
                *ids.entry(k).or_insert_with(|| {
                    verts.push(p);
                    (verts.len() - 1) as u32
                })
            });
            if f[0] != f[1] && f[1] != f[2] && f[0] != f[2] {
                faces.push(f);
            }
        }
        Indexed { verts, faces }
    }

    pub fn tris(&self) -> Vec<Tri> {
        self.faces.iter().map(|f| f.map(|i| self.verts[i as usize])).collect()
    }

    /// Normale non unitaire (deux fois l'aire) de la face `f`.
    pub fn face_normal(&self, f: usize) -> [f64; 3] {
        let [a, b, c] = self.faces[f].map(|i| self.verts[i as usize]);
        cross(sub(b, a), sub(c, a))
    }

    /// Volume signé (positif : normales vers l'extérieur).
    pub fn signed_volume(&self) -> f64 {
        self.faces
            .iter()
            .map(|f| {
                let [a, b, c] = f.map(|i| self.verts[i as usize]);
                dot(a, cross(b, c)) / 6.0
            })
            .sum()
    }

    /// Faces de chaque arête (sommets triés) : une (bord ouvert), deux, ou plus.
    pub fn edge_faces(&self) -> HashMap<(u32, u32), Vec<u32>> {
        let mut m: HashMap<(u32, u32), Vec<u32>> = HashMap::new();
        for (k, f) in self.faces.iter().enumerate() {
            for (a, b) in [(f[0], f[1]), (f[1], f[2]), (f[2], f[0])] {
                m.entry((a.min(b), a.max(b))).or_default().push(k as u32);
            }
        }
        m
    }
}

/// Arête vive : segment du maillage où la surface tourne brusquement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FeatureEdge {
    pub a: [f64; 3],
    pub b: [f64; 3],
    /// Angle entre les normales des deux faces (°) ; 180 pour un bord ouvert.
    pub angle_deg: f64,
    /// Arête saillante (bord de fuite, saumon plat) ; sinon rentrante
    /// (raccord aile–fuselage) ou bord ouvert.
    pub convex: bool,
}

/// Seuil d'arête vive (°) : celui du critère par maille du solveur (faces
/// touchées par les liens d'une maille à plus de 60° l'une de l'autre). Les
/// facettes d'une surface lisse font quelques degrés entre elles, une arête
/// vive 90° (saumon plat) à 165° (bord de fuite) : le choix entre 30° et 80°
/// ne change pas le classement (voir `feature_edges` et ses tests).
pub const FEATURE_DEG: f64 = 60.0;

/// Arêtes vives (`crease_set`) avec leur angle et leur sens. Les normales
/// sont orientées vers l'extérieur d'après le signe du volume.
pub fn feature_edges(m: &Indexed, min_deg: f64) -> Vec<FeatureEdge> {
    let sgn = if m.signed_volume() >= 0.0 { 1.0 } else { -1.0 };
    let ef = m.edge_faces();
    let mut keys: Vec<_> = crease_set(m, min_deg).into_iter().collect();
    keys.sort_unstable();
    keys.into_iter()
        .map(|(i, j)| {
            let (a, b) = (m.verts[i as usize], m.verts[j as usize]);
            let fs = &ef[&(i, j)];
            if fs.len() != 2 {
                return FeatureEdge { a, b, angle_deg: 180.0, convex: false };
            }
            let n0 = m.face_normal(fs[0] as usize);
            let n1 = m.face_normal(fs[1] as usize);
            let c = dot(n0, n1) / (norm(n0) * norm(n1)).max(1e-300);
            // Saillante si le sommet opposé de la face 1 est derrière le plan
            // de la face 0 (normales extérieures).
            let f1 = m.faces[fs[1] as usize];
            let o = f1.iter().find(|&&v| v != i && v != j).map(|&v| m.verts[v as usize]).unwrap_or(a);
            let convex = sgn * dot(n0, sub(o, a)) < 0.0;
            FeatureEdge { a, b, angle_deg: c.clamp(-1.0, 1.0).acos().to_degrees(), convex }
        })
        .collect()
}

/// Distance d'un point à un segment.
pub fn seg_dist(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let e = sub(b, a);
    let t = (dot(sub(p, a), e) / dot(e, e).max(1e-300)).clamp(0.0, 1.0);
    norm(sub(p, [a[0] + t * e[0], a[1] + t * e[1], a[2] + t * e[2]]))
}

/// Arêtes vives rangées dans une grille de cases unité, pour la distance.
#[derive(Clone, Debug, Default)]
pub struct FeatureSet {
    pub edges: Vec<FeatureEdge>,
    lo: [f64; 3],
    dim: [usize; 3],
    bins: Vec<Vec<u32>>,
}

impl FeatureSet {
    pub fn new(edges: Vec<FeatureEdge>) -> FeatureSet {
        if edges.is_empty() {
            return FeatureSet::default();
        }
        let mut lo = [f64::MAX; 3];
        let mut hi = [f64::MIN; 3];
        for e in &edges {
            for p in [e.a, e.b] {
                for k in 0..3 {
                    lo[k] = lo[k].min(p[k]);
                    hi[k] = hi[k].max(p[k]);
                }
            }
        }
        let dim = [0, 1, 2].map(|k| ((hi[k] - lo[k]).floor() as usize + 1).min(4096));
        let mut bins = vec![Vec::new(); dim[0] * dim[1] * dim[2]];
        let cellof = |v: f64, k: usize| (((v - lo[k]).floor().max(0.0)) as usize).min(dim[k] - 1);
        for (n, e) in edges.iter().enumerate() {
            let g0 = [0, 1, 2].map(|k| cellof(e.a[k].min(e.b[k]), k));
            let g1 = [0, 1, 2].map(|k| cellof(e.a[k].max(e.b[k]), k));
            for z in g0[2]..=g1[2] {
                for y in g0[1]..=g1[1] {
                    for x in g0[0]..=g1[0] {
                        bins[x + dim[0] * (y + dim[1] * z)].push(n as u32);
                    }
                }
            }
        }
        FeatureSet { edges, lo, dim, bins }
    }

    /// Distance à l'arête vive la plus proche, plafonnée à `cap`.
    pub fn distance(&self, p: [f64; 3], cap: f64) -> f64 {
        if self.edges.is_empty() {
            return cap;
        }
        let mut best = cap;
        let r = |v: f64, k: usize| (v - self.lo[k]).floor();
        let g0 = [0, 1, 2].map(|k| r(p[k] - cap, k));
        let g1 = [0, 1, 2].map(|k| r(p[k] + cap, k));
        if (0..3).any(|k| g1[k] < 0.0 || g0[k] > (self.dim[k] - 1) as f64) {
            return cap;
        }
        let c = |v: f64, k: usize| v.clamp(0.0, (self.dim[k] - 1) as f64) as usize;
        for z in c(g0[2], 2)..=c(g1[2], 2) {
            for y in c(g0[1], 1)..=c(g1[1], 1) {
                for x in c(g0[0], 0)..=c(g1[0], 0) {
                    for &n in &self.bins[x + self.dim[0] * (y + self.dim[1] * z)] {
                        let e = &self.edges[n as usize];
                        best = best.min(seg_dist(p, e.a, e.b));
                    }
                }
            }
        }
        best
    }

    /// Le même ensemble transformé point par point.
    pub fn transformed(&self, t: impl Fn([f64; 3]) -> [f64; 3]) -> FeatureSet {
        FeatureSet::new(self.edges.iter().map(|e| FeatureEdge { a: t(e.a), b: t(e.b), ..*e }).collect())
    }
}

/// Coupe d'une tranche d'envergure y = cte : segments (x, z).
pub fn slice_y(tris: &[Tri], y: f64) -> Vec<[[f64; 2]; 2]> {
    let mut out = Vec::new();
    for t in tris {
        let mut pts: Vec<[f64; 2]> = Vec::with_capacity(2);
        for k in 0..3 {
            let (p, q) = (t[k], t[(k + 1) % 3]);
            let (dp, dq) = (p[1] - y, q[1] - y);
            if (dp < 0.0) != (dq < 0.0) {
                let s = dp / (dp - dq);
                pts.push([p[0] + s * (q[0] - p[0]), p[2] + s * (q[2] - p[2])]);
            }
        }
        if pts.len() == 2 {
            out.push([pts[0], pts[1]]);
        }
    }
    out
}

/// Profil d'une tranche d'envergure, longueurs dans l'unité des triangles.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Section {
    pub y: f64,
    /// Corde (distance du bord de fuite au point le plus éloigné).
    pub chord: f64,
    /// Épaisseur relative maximale.
    pub thickness: f64,
    /// Rayon du bord d'attaque, ajusté sur y_demi² = 2 R s entre 0,2 % et
    /// 2 % de corde (s : distance au bord d'attaque le long de la corde).
    pub le_radius: f64,
    /// Épaisseur du culot du bord de fuite (0 : bord de fuite fermé).
    pub te_thickness: f64,
    /// Longueur (le long de la corde, depuis le bord de fuite) où l'épaisseur
    /// locale est inférieure à 1 et à 2 unités.
    pub thin1: f64,
    pub thin2: f64,
    /// Longueur moyenne des facettes au bord d'attaque (segments de la coupe
    /// à moins de 2 % de corde du bord d'attaque).
    pub le_facet: f64,
    /// Plus grand angle (°) entre deux facettes successives de la coupe à
    /// moins de 2 % de corde du bord d'attaque : quelques degrés pour une
    /// surface lisse finement facettée ; au-delà de 60°, le critère d'arête
    /// vive par maille du solveur peut y couper la loi de paroi.
    pub le_corner_deg: f64,
}

/// Analyse d'une tranche : None si la coupe est vide ou dégénérée.
pub fn section(tris: &[Tri], y: f64) -> Option<Section> {
    let segs = slice_y(tris, y);
    if segs.len() < 6 {
        return None;
    }
    let pts: Vec<[f64; 2]> = segs.iter().flatten().copied().collect();
    // Bord de fuite : point le plus en aval ; bord d'attaque : le plus éloigné de lui.
    let te = *pts.iter().max_by(|a, b| a[0].total_cmp(&b[0]))?;
    let d2 = |p: &[f64; 2], q: &[f64; 2]| (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2);
    let le = *pts.iter().max_by(|a, b| d2(a, &te).total_cmp(&d2(b, &te)))?;
    let chord = d2(&le, &te).sqrt();
    if chord <= 0.0 {
        return None;
    }
    // Repère de corde : s le long de la corde depuis le bord d'attaque, n normal.
    let ex = [(te[0] - le[0]) / chord, (te[1] - le[1]) / chord];
    let to_sn = |p: [f64; 2]| {
        let (dx, dz) = (p[0] - le[0], p[1] - le[1]);
        [dx * ex[0] + dz * ex[1], -dx * ex[1] + dz * ex[0]]
    };
    let sn: Vec<[[f64; 2]; 2]> = segs.iter().map(|s| s.map(to_sn)).collect();
    // Épaisseur locale à l'abscisse s : écart entre les coupes extrêmes.
    let width = |s: f64| {
        let (mut a, mut b) = (f64::MAX, f64::MIN);
        for g in &sn {
            let (p, q) = (g[0], g[1]);
            if (p[0] - s) * (q[0] - s) <= 0.0 && (p[0] - q[0]).abs() > 1e-300 {
                let t = (s - p[0]) / (q[0] - p[0]);
                let n = p[1] + t * (q[1] - p[1]);
                a = a.min(n);
                b = b.max(n);
            }
        }
        (b > a).then_some(b - a)
    };
    let mut tmax = 0.0f64;
    let mut thin1 = 0.0;
    let mut thin2 = 0.0;
    let ns = 400;
    for k in 1..ns {
        let s = chord * k as f64 / ns as f64;
        let w = width(s).unwrap_or(0.0);
        tmax = tmax.max(w);
        let ds = chord / ns as f64;
        if w < 1.0 && s > 0.5 * chord {
            thin1 += ds;
        }
        if w < 2.0 && s > 0.5 * chord {
            thin2 += ds;
        }
    }
    // Culot : épaisseur juste avant le bord de fuite.
    let te_thickness = width(chord * (1.0 - 1e-6)).unwrap_or(0.0);
    // Rayon du bord d'attaque : développement de la demi-épaisseur près du
    // bord d'attaque, y = a √s + b s + …, d'où y² = a² s + 2ab s^1,5 + … et
    // R = a²/2 ; moindres carrés à deux termes (sans le second, le rayon
    // d'un NACA 0012 sort 10 % trop petit).
    let mut m = [[0.0f64; 2]; 2];
    let mut r = [0.0f64; 2];
    // Abscisses réparties en progression géométrique de 0,2 % à 2 % de corde :
    // l'épaisseur y est celle du polygone des facettes, que voit le solveur.
    for k in 0..24 {
        let s = chord * 0.002 * 10f64.powf(k as f64 / 23.0);
        if let Some(w) = width(s) {
            let yh = 0.5 * w;
            let f = [s, s * s.sqrt()];
            for i in 0..2 {
                for j in 0..2 {
                    m[i][j] += f[i] * f[j];
                }
                r[i] += f[i] * yh * yh;
            }
        }
    }
    let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    let le_radius = if det.abs() > 1e-300 { 0.5 * (r[0] * m[1][1] - r[1] * m[0][1]) / det } else { 0.0 };
    // Facettes du bord d'attaque.
    let (mut fl, mut fn_) = (0.0, 0);
    for g in &sn {
        if g[0][0] < 0.02 * chord && g[1][0] < 0.02 * chord {
            fl += d2(&g[0], &g[1]).sqrt();
            fn_ += 1;
        }
    }
    let le_facet = if fn_ > 0 { fl / fn_ as f64 } else { 0.0 };
    // Coins de la coupe : directions des deux segments issus de chaque sommet.
    let key = |p: [f64; 2]| [(p[0] * 1e6).round() as i64, (p[1] * 1e6).round() as i64];
    let mut at: HashMap<[i64; 2], Vec<[f64; 2]>> = HashMap::new();
    for g in &sn {
        for (a, b) in [(g[0], g[1]), (g[1], g[0])] {
            let l = d2(&a, &b).sqrt();
            if l > 1e-9 * chord {
                at.entry(key(a)).or_default().push([(b[0] - a[0]) / l, (b[1] - a[1]) / l]);
            }
        }
    }
    let mut le_corner_deg = 0.0f64;
    for g in &sn {
        for p in g {
            if p[0] < 0.02 * chord {
                if let Some(d) = at.get(&key(*p)).filter(|d| d.len() == 2) {
                    let c = (d[0][0] * d[1][0] + d[0][1] * d[1][1]).clamp(-1.0, 1.0);
                    le_corner_deg = le_corner_deg.max(180.0 - c.acos().to_degrees());
                }
            }
        }
    }
    Some(Section {
        y,
        chord,
        thickness: tmax / chord,
        le_radius,
        te_thickness,
        thin1,
        thin2,
        le_facet,
        le_corner_deg,
    })
}

/// Facteur de la tessellation de Phong (Boubekeur et Alexa, *ACM Trans.
/// Graph.* 27, 2008) : 1/2 rend un arc de cercle entre deux sommets dont les
/// normales sont radiales, au second ordre près (la flèche d'une corde de
/// demi-angle θ/2 est R θ²/8 ; la projection sur les plans tangents des deux
/// extrémités pousse le milieu de R θ²/4, dont la moitié). La valeur 3/4 de
/// l'article, choisie pour l'aspect, gonflerait les arêtes de 50 %.
pub const PHONG_ALPHA: f64 = 0.5;

/// Angle maximal (°) entre sous-facettes voisines après raffinement : la
/// flèche résiduelle d'un arc de rayon R est R (10°)²/8 ≈ 0,4 % de R.
pub const SMOOTH_DEG: f64 = 10.0;

/// Angle (°) entre les normales des deux faces de chaque arête ; 180 pour un
/// bord ouvert ou une arête partagée par plus de deux faces.
pub fn dihedral_angles(m: &Indexed) -> HashMap<(u32, u32), f64> {
    m.edge_faces()
        .into_iter()
        .map(|(e, fs)| {
            let ang = if fs.len() == 2 {
                let (a, b) = (m.face_normal(fs[0] as usize), m.face_normal(fs[1] as usize));
                (dot(a, b) / (norm(a) * norm(b)).max(1e-300)).clamp(-1.0, 1.0).acos().to_degrees()
            } else {
                180.0
            };
            (e, ang)
        })
        .collect()
}

/// Arêtes vives (paires de sommets triées) : bords ouverts, et arêtes dont
/// l'angle dièdre dépasse `min_deg`.
///
/// L'angle seul ne distingue pas une arête vive d'un bord d'attaque arrondi
/// décrit par trop peu de facettes. Sur les maillages exportés d'OpenVSP, le
/// sommet du bord d'attaque fait 48° (aile AR 4) et 59° (TR-1208) d'une
/// facette à l'autre, juste sous le seuil, les arêtes vives 67° à 165°. Un
/// critère de concentration du virage (rangée voisine deux fois moins
/// courbée) a été essayé : sur TR-1208, le rapport vaut exactement 2 au bord
/// d'attaque, il ne sépare rien de façon sûre. Pour un maillage plus
/// grossier, les arêtes vives se corrigent à la main dans un fichier `.sgeo`.
pub fn crease_set(m: &Indexed, min_deg: f64) -> HashSet<(u32, u32)> {
    let ef = m.edge_faces();
    dihedral_angles(m)
        .into_iter()
        .filter(|(e, th)| ef[e].len() != 2 || *th > min_deg)
        .map(|(e, _)| e)
        .collect()
}

/// Plus grand angle (°) entre faces voisines de part et d'autre d'une arête
/// qui n'est pas vive.
pub fn max_smooth_angle(m: &Indexed, creases: &HashSet<(u32, u32)>) -> f64 {
    m.edge_faces()
        .into_iter()
        .filter(|(e, fs)| fs.len() == 2 && !creases.contains(e))
        .map(|(_, fs)| {
            let (a, b) = (m.face_normal(fs[0] as usize), m.face_normal(fs[1] as usize));
            (dot(a, b) / (norm(a) * norm(b)).max(1e-300)).clamp(-1.0, 1.0).acos().to_degrees()
        })
        .fold(0.0, f64::max)
}

fn find(p: &mut [usize], mut x: usize) -> usize {
    while p[x] != x {
        p[x] = p[p[x]];
        x = p[x];
    }
    x
}

/// Normales aux coins des faces (indice 3 f + k) : moyenne, pondérée par
/// l'angle au sommet, des normales des faces qui entourent le sommet sans
/// traverser d'arête vive (groupe lisse du coin).
pub fn corner_normals(m: &Indexed, creases: &HashSet<(u32, u32)>) -> Vec<[f64; 3]> {
    let nc = 3 * m.faces.len();
    let mut parent: Vec<usize> = (0..nc).collect();
    let corner = |f: u32, v: u32| 3 * f as usize + m.faces[f as usize].iter().position(|&w| w == v).unwrap_or(0);
    for (e, fs) in m.edge_faces() {
        if fs.len() == 2 && !creases.contains(&e) {
            for v in [e.0, e.1] {
                let (a, b) = (corner(fs[0], v), corner(fs[1], v));
                let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
                parent[ra] = rb;
            }
        }
    }
    let mut acc = vec![[0.0f64; 3]; nc];
    for (f, face) in m.faces.iter().enumerate() {
        let n = m.face_normal(f);
        let ln = norm(n).max(1e-300);
        for k in 0..3 {
            let p = m.verts[face[k] as usize];
            let a = sub(m.verts[face[(k + 1) % 3] as usize], p);
            let b = sub(m.verts[face[(k + 2) % 3] as usize], p);
            let ang = (dot(a, b) / (norm(a) * norm(b)).max(1e-300)).clamp(-1.0, 1.0).acos();
            let r = find(&mut parent, 3 * f + k);
            for c in 0..3 {
                acc[r][c] += ang * n[c] / ln;
            }
        }
    }
    (0..nc)
        .map(|c| {
            let v = acc[find(&mut parent, c)];
            let l = norm(v).max(1e-300);
            v.map(|x| x / l)
        })
        .collect()
}

/// Point de Phong : barycentre `w` des sommets `p`, projeté sur les plans
/// tangents des trois sommets (normales `n`), mêlé au barycentre plat.
fn phong(p: [[f64; 3]; 3], n: [[f64; 3]; 3], w: [f64; 3]) -> [f64; 3] {
    let x = [0, 1, 2].map(|c| w[0] * p[0][c] + w[1] * p[1][c] + w[2] * p[2][c]);
    let mut q = [0.0; 3];
    for i in 0..3 {
        let d = dot(sub(x, p[i]), n[i]);
        for c in 0..3 {
            q[c] += w[i] * (x[c] - d * n[i][c]);
        }
    }
    [0, 1, 2].map(|c| (1.0 - PHONG_ALPHA) * x[c] + PHONG_ALPHA * q[c])
}

/// Raffinement lisse d'un maillage facetté : chaque face est découpée en n²
/// sous-faces dont les sommets sont placés sur la surface de Phong qui
/// interpole les sommets d'origine et leurs normales. Les arêtes vives
/// (`creases`) restent droites, et le maillage reste fermé : les points d'une
/// arête qui n'est pas vive ne dépendent que de ses deux bouts et de leurs
/// normales, communes aux deux faces (même groupe lisse). Rend une surface
/// lisse là où le fichier ne donne qu'un polygone (bord d'attaque facetté en
/// quelques cordes), sans rien changer aux arêtes vives (bord de fuite,
/// saumons) ni aux faces planes.
pub fn smooth_refine(m: &Indexed, creases: &HashSet<(u32, u32)>, n: usize) -> Indexed {
    let n = n.max(1);
    if n == 1 {
        return m.clone();
    }
    let cn = corner_normals(m, creases);
    let mut verts = m.verts.clone();
    let mut faces = Vec::with_capacity(m.faces.len() * n * n);
    // Points intérieurs des arêtes, rangés depuis le plus petit indice.
    let mut edge_pts: HashMap<(u32, u32), Vec<u32>> = HashMap::new();
    let nf = n as f64;
    // Rang de (i, j), i + j ≤ n, rangées de j constant.
    let idx = |i: usize, j: usize| j * (n + 1) - j * j.saturating_sub(1) / 2 + i;
    for (f, face) in m.faces.iter().enumerate() {
        let p = face.map(|v| m.verts[v as usize]);
        let nn = [cn[3 * f], cn[3 * f + 1], cn[3 * f + 2]];
        let mut ids = vec![0u32; (n + 1) * (n + 2) / 2];
        for j in 0..=n {
            for i in 0..=(n - j) {
                let id = if (i, j) == (0, 0) {
                    face[0]
                } else if (i, j) == (n, 0) {
                    face[1]
                } else if (i, j) == (0, n) {
                    face[2]
                } else {
                    // Arête (k0 → k1, s-ième point sur n), ou intérieur.
                    let edge = if j == 0 {
                        Some((0, 1, i))
                    } else if i + j == n {
                        Some((1, 2, j))
                    } else if i == 0 {
                        Some((2, 0, n - j))
                    } else {
                        None
                    };
                    match edge {
                        None => {
                            let w = [(n - i - j) as f64 / nf, i as f64 / nf, j as f64 / nf];
                            verts.push(phong(p, nn, w));
                            (verts.len() - 1) as u32
                        }
                        Some((k0, k1, s)) => {
                            let (a, b) = (face[k0], face[k1]);
                            let key = (a.min(b), a.max(b));
                            let list = edge_pts.entry(key).or_insert_with(|| {
                                let (ka, kb) = if a == key.0 { (k0, k1) } else { (k1, k0) };
                                let (pa, pb) = (p[ka], p[kb]);
                                let crease = creases.contains(&key);
                                (1..n)
                                    .map(|t| {
                                        let u = t as f64 / nf;
                                        let q = if crease {
                                            [0, 1, 2].map(|c| (1.0 - u) * pa[c] + u * pb[c])
                                        } else {
                                            phong([pa, pb, pa], [nn[ka], nn[kb], nn[ka]], [1.0 - u, u, 0.0])
                                        };
                                        verts.push(q);
                                        (verts.len() - 1) as u32
                                    })
                                    .collect()
                            });
                            list[if a == key.0 { s } else { n - s } - 1]
                        }
                    }
                };
                ids[idx(i, j)] = id;
            }
        }
        for j in 0..n {
            for i in 0..(n - j) {
                faces.push([ids[idx(i, j)], ids[idx(i + 1, j)], ids[idx(i, j + 1)]]);
                if i + j + 1 < n {
                    faces.push([ids[idx(i + 1, j)], ids[idx(i + 1, j + 1)], ids[idx(i, j + 1)]]);
                }
            }
        }
    }
    Indexed { verts, faces }
}

/// Nombre de découpes par arête pour que les sous-facettes ne fassent pas
/// plus de `SMOOTH_DEG` entre elles (1 : rien à faire), plafonné à 16.
pub fn smooth_steps(m: &Indexed, creases: &HashSet<(u32, u32)>) -> usize {
    ((max_smooth_angle(m, creases) / SMOOTH_DEG).ceil() as usize).clamp(1, 16)
}

/// Raffinement lisse de triangles (STL) : arêtes vives au-delà de
/// `FEATURE_DEG`, découpes selon `smooth_steps`. Rend aussi le nombre de
/// découpes par arête.
pub fn smooth_tris(tris: &[Tri]) -> (Vec<Tri>, usize) {
    let m = Indexed::from_tris(tris);
    let cr = crease_set(&m, FEATURE_DEG);
    let n = smooth_steps(&m, &cr);
    (smooth_refine(&m, &cr, n).tris(), n)
}

/// Volume (unités³) de boîtes alignées sur les axes couvrant les faces
/// tournées vers l'amont (n_x < −0,5, comme le raffinement du bord
/// d'attaque), coupées en `n` tronçons d'envergure égaux : marge `e` de
/// chaque côté, `down` vers l'aval. n = 1 : la boîte unique du solveur.
pub fn le_cover(tris: &[Tri], n: usize, e: f64, down: f64) -> f64 {
    let m = Indexed::from_tris(tris);
    let sgn = if m.signed_volume() >= 0.0 { 1.0 } else { -1.0 };
    let up: Vec<Tri> = tris
        .iter()
        .filter(|t| {
            let nv = cross(sub(t[1], t[0]), sub(t[2], t[0]));
            let l = norm(nv);
            l > 1e-12 && sgn * nv[0] / l < -0.5
        })
        .copied()
        .collect();
    if up.is_empty() {
        return 0.0;
    }
    let (lo, hi) = bbox(&up);
    let mut v = 0.0;
    for k in 0..n {
        let (y0, y1) = (lo[1] + (hi[1] - lo[1]) * k as f64 / n as f64, lo[1] + (hi[1] - lo[1]) * (k + 1) as f64 / n as f64);
        let (mut a, mut b) = ([f64::MAX; 3], [f64::MIN; 3]);
        for t in &up {
            for p in t {
                if p[1] >= y0 - 1e-9 && p[1] <= y1 + 1e-9 {
                    for c in 0..3 {
                        a[c] = a[c].min(p[c]);
                        b[c] = b[c].max(p[c]);
                    }
                }
            }
        }
        if a[0] > b[0] {
            continue;
        }
        // Tronçons voisins : la marge latérale n'est comptée qu'aux bouts.
        let ey = if n == 1 { 2.0 * e } else if k == 0 || k == n - 1 { e } else { 0.0 };
        v += (b[0] - a[0] + e + down) * (y1 - y0 + ey) * (b[2] - a[2] + 2.0 * e);
    }
    v
}

/// Tranches régulières sur l'envergure (bouts exclus : 2 % de chaque côté).
pub fn sections(tris: &[Tri], n: usize) -> Vec<Section> {
    let (lo, hi) = bbox(tris);
    (0..n)
        .filter_map(|k| {
            let f = 0.02 + 0.96 * (k as f64 + 0.5) / n as f64;
            section(tris, lo[1] + f * (hi[1] - lo[1]))
        })
        .collect()
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// Aile rectangulaire NACA 00tt (bord de fuite fermé) à saumons plats,
    /// `nc` points par face, corde c, envergure b.
    pub(crate) fn naca_wing(t: f64, c: f64, b: f64, nc: usize) -> Vec<Tri> {
        let yt = |x: f64| {
            5.0 * t * c * (0.2969 * x.sqrt() - 0.126 * x - 0.3516 * x * x + 0.2843 * x.powi(3) - 0.1036 * x.powi(4))
        };
        // Cosinus : points serrés aux bords d'attaque et de fuite ; bord de fuite fermé.
        let xs: Vec<f64> = (0..=nc).map(|k| 0.5 * (1.0 - (std::f64::consts::PI * k as f64 / nc as f64).cos())).collect();
        let mut prof: Vec<[f64; 2]> = Vec::new();
        for &x in xs.iter().rev() {
            prof.push([x * c, if x >= 1.0 { 0.0 } else { yt(x) }]);
        }
        for &x in xs.iter().skip(1).take(nc - 1) {
            prof.push([x * c, -yt(x)]);
        }
        let n = prof.len();
        let p = |k: usize, y: f64| [prof[k % n][0], y, prof[k % n][1]];
        let mut tris = Vec::new();
        // Quatre tronçons d'envergure.
        for s in 0..4 {
            let (y0, y1) = (b * s as f64 / 4.0, b * (s + 1) as f64 / 4.0);
            for k in 0..n {
                tris.push([p(k, y0), p(k + 1, y1), p(k + 1, y0)]);
                tris.push([p(k, y0), p(k, y1), p(k + 1, y1)]);
            }
        }
        // Saumons : éventail depuis le bord d'attaque.
        for k in 1..n - 1 {
            tris.push([p(nc, 0.0), p(nc + k, 0.0), p(nc + k + 1, 0.0)]);
            tris.push([p(nc, b), p(nc + k + 1, b), p(nc + k, b)]);
        }
        tris
    }

    #[test]
    fn aretes_vives_aile() {
        let tris = naca_wing(0.12, 10.0, 40.0, 40);
        let m = Indexed::from_tris(&tris);
        let (_, d) = crate::stl::clean(tris.clone());
        assert!(d.watertight() && d.flipped_edges == 0, "{d:?}");
        assert!(m.signed_volume() > 0.0);
        let f = feature_edges(&m, FEATURE_DEG);
        // Bord de fuite (une arête, 165°) et pourtours des deux saumons (90°).
        let te: Vec<_> = f.iter().filter(|e| e.a[0] > 9.99 && e.b[0] > 9.99 && (e.a[1] - e.b[1]).abs() > 1.0).collect();
        assert_eq!(te.len(), 4, "{te:?}");
        assert!(te[0].angle_deg > 150.0 && te[0].convex);
        let tips = f.iter().filter(|e| e.a[1] == e.b[1]).count();
        assert_eq!(tips, 2 * (2 * 40), "arêtes de saumon");
        assert!(f.iter().all(|e| e.convex), "toutes saillantes");
        // Le seuil ne change rien entre 30° et 80°.
        assert_eq!(feature_edges(&m, 30.0).len(), f.len());
        assert_eq!(feature_edges(&m, 80.0).len(), f.len());
        // Distance au bord de fuite.
        let fs = FeatureSet::new(f);
        assert!((fs.distance([11.0, 20.0, 0.0], 5.0) - 1.0).abs() < 1e-9);
        assert_eq!(fs.distance([5.0, 20.0, 3.0], 2.0), 2.0);
    }

    #[test]
    fn coupe_naca0012() {
        // Corde 100 : rayon du bord d'attaque 1,1019 t² c = 1,587.
        let tris = naca_wing(0.12, 100.0, 40.0, 200);
        let s = section(&tris, 20.0).unwrap();
        assert!((s.chord - 100.0).abs() < 1e-6, "{s:?}");
        assert!((s.thickness - 0.12).abs() < 0.002, "{s:?}");
        assert!((s.le_radius / 1.587 - 1.0).abs() < 0.1, "rayon {}", s.le_radius);
        assert!(s.te_thickness < 1e-3, "{s:?}");
        // Épaisseur < 1 sur les derniers ≈ 3,6 % de corde (yt(0,964) ≈ 0,5).
        assert!(s.thin1 > 2.5 && s.thin1 < 5.0, "{s:?}");
    }

    #[test]
    fn raffinement_lisse_sphere() {
        // Icosaèdre subdivisé une fois : 80 faces, ≈ 32° entre faces voisines.
        let tris = crate::stl::icosphere(10.0, [0.0; 3], 1);
        let m = Indexed::from_tris(&tris);
        let cr = crease_set(&m, FEATURE_DEG);
        assert!(cr.is_empty());
        let n = smooth_steps(&m, &cr);
        assert!((3..=6).contains(&n), "{n}");
        let r = smooth_refine(&m, &cr, n);
        assert_eq!(r.faces.len(), 80 * n * n);
        let (_, d) = crate::stl::clean(r.tris());
        assert!(d.watertight() && d.flipped_edges == 0, "{d:?}");
        let exact = 4.0 / 3.0 * std::f64::consts::PI * 1000.0;
        let (v0, v1) = (m.signed_volume(), r.signed_volume());
        // Polygone : 14 % de volume en moins ; surface de Phong : moins de 2 %.
        assert!(v0 < 0.9 * exact && (v1 / exact - 1.0).abs() < 0.03, "{v0} {v1} {exact}");
        let dev = r.verts.iter().map(|p| (norm(*p) - 10.0).abs()).fold(0.0, f64::max);
        assert!(dev < 0.25, "écart au rayon {dev}");
    }

    #[test]
    fn raffinement_lisse_aile() {
        // NACA 0012 grossier (24 points par face) : polygone au bord
        // d'attaque, 53° d'une facette à l'autre.
        let tris = naca_wing(0.12, 100.0, 40.0, 24);
        let m = Indexed::from_tris(&tris);
        let cr = crease_set(&m, FEATURE_DEG);
        let n = smooth_steps(&m, &cr);
        let r = smooth_refine(&m, &cr, n);
        let (_, d) = crate::stl::clean(r.tris());
        assert!(d.watertight() && d.flipped_edges == 0, "{d:?}");
        let le = dihedral_angles(&m).values().filter(|&&a| a > 40.0 && a < FEATURE_DEG).count();
        assert!(le > 0, "bord d'attaque trop finement facetté pour l'essai");
        // Arêtes vives conservées : bord de fuite et pourtours des saumons.
        let fe0 = feature_edges(&m, FEATURE_DEG).len();
        let fe1 = feature_edges(&r, FEATURE_DEG).len();
        assert_eq!(fe1, fe0 * n, "arêtes vives {fe0} → {fe1}");
        let s0 = section(&m.tris(), 20.3).unwrap();
        let s1 = section(&r.tris(), 20.3).unwrap();
        let fine = section(&naca_wing(0.12, 100.0, 40.0, 400), 20.3).unwrap();
        // Bord d'attaque : coin réduit, bord de fuite toujours fermé.
        assert!(s1.le_corner_deg < 0.5 * s0.le_corner_deg, "{s0:?} {s1:?}");
        assert!(s1.te_thickness < 1e-3, "{s1:?}");
        // Plus près du profil exact : épaisseur, et volume.
        let exact = naca_wing(0.12, 100.0, 40.0, 400);
        let ve = Indexed::from_tris(&exact).signed_volume();
        let (v0, v1) = (m.signed_volume(), r.signed_volume());
        assert!((v1 - ve).abs() < 0.5 * (v0 - ve).abs(), "volumes {v0} {v1} {ve}");
        assert!((s1.thickness - fine.thickness).abs() <= (s0.thickness - fine.thickness).abs() + 1e-4, "{s0:?} {s1:?}");
    }
}
