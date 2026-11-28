//! Éléments de surface (« surfels ») de la paroi volumique de Chen et al.
//! (1998), celle de PowerFLOW : géométrie seulement, côté processeur.
//!
//! Chaque triangle du maillage est découpé par les plans de la grille
//! (Sutherland–Hodgman) ; chaque morceau, contenu dans un seul voxel, est un
//! surfel : voxel, aire, normale (sortante du corps, donc vers le fluide),
//! centre. Le voxel du nœud (i, j, k) est le cube [i − ½, i + ½] × … centré
//! sur le nœud.
//!
//! Pour chaque surfel et chaque direction D3Q19 cᵢ (cᵢ·n ≠ 0), le
//! parallélépipède balayé pendant un pas de temps, côté fluide, a pour volume
//! Pᵢ = |cᵢ·n| A (Δt = 1). Il sert au rassemblement des populations qui
//! arrivent sur le surfel (cᵢ·n < 0) et à la dispersion de celles qui en
//! repartent (cᵢ·n > 0) : c'est le même volume pour cᵢ et −cᵢ (le surfel
//! poussé de 0 à c dans le fluide). Les poids sont les volumes exacts de son
//! intersection avec chaque voxel (découpe d'un polyèdre convexe).
//!
//! La fraction solide de chaque voxel coupé est exacte (théorème de la
//! divergence le long des colonnes en x, à partir des seuls morceaux de
//! surface) ; elle suppose un maillage fermé.
//!
//! Raffinement imbriqué : chaque niveau (grille 2ᵏ fois plus fine) a un
//! domaine propre ; un morceau de surface appartient au niveau le plus fin
//! dont le domaine propre le contient (découpe exacte par les plans de ces
//! domaines), et ses surfels sont dans la grille de ce niveau.

use std::collections::HashMap;

use rayon::prelude::*;

use crate::stl::Tri;

type P = [f64; 3];

/// Directions D3Q19, même ordre que `gpu3d::C3` : la direction opposée de
/// 2p + 1 est 2p + 2.
pub const C19: [[i32; 3]; 19] = [
    [0, 0, 0],
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1],
    [1, 1, 0],
    [-1, -1, 0],
    [1, 0, 1],
    [-1, 0, -1],
    [0, 1, 1],
    [0, -1, -1],
    [1, -1, 0],
    [-1, 1, 0],
    [1, 0, -1],
    [-1, 0, 1],
    [0, 1, -1],
    [0, -1, 1],
];

fn sub(a: P, b: P) -> P {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn add(a: P, b: P) -> P {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn dot(a: P, b: P) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: P, b: P) -> P {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn scale(a: P, s: f64) -> P {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn norm(a: P) -> f64 {
    dot(a, a).sqrt()
}
fn cvec(i: usize) -> P {
    C19[i].map(|v| v as f64)
}

/// Vecteur aire (½ Σ vₖ × vₖ₊₁) et centre d'un polygone plan.
fn area_centroid(poly: &[P]) -> (P, P) {
    let o = poly[0];
    let mut a = [0.0; 3];
    let mut m = [0.0; 3];
    for k in 1..poly.len() - 1 {
        let t = cross(sub(poly[k], o), sub(poly[k + 1], o));
        a = add(a, t);
        let g = scale(add(add(poly[k], poly[k + 1]), o), 1.0 / 3.0);
        m = add(m, scale(sub(g, o), norm(t)));
    }
    let s: f64 = (1..poly.len() - 1).map(|k| norm(cross(sub(poly[k], o), sub(poly[k + 1], o)))).sum();
    let c = if s > 0.0 { add(o, scale(m, 1.0 / s)) } else { o };
    (scale(a, 0.5), c)
}

/// Coupe un polygone convexe par le plan x_a = x : (côté bas, côté haut).
/// Les points d'intersection sont posés exactement sur le plan.
fn split(poly: &[P], a: usize, x: f64) -> (Vec<P>, Vec<P>) {
    let (mut lo, mut hi) = (Vec::with_capacity(poly.len() + 1), Vec::with_capacity(poly.len() + 1));
    for k in 0..poly.len() {
        let (p, q) = (poly[k], poly[(k + 1) % poly.len()]);
        let (dp, dq) = (p[a] - x, q[a] - x);
        if dp <= 0.0 {
            lo.push(p);
        }
        if dp >= 0.0 {
            hi.push(p);
        }
        if (dp < 0.0 && dq > 0.0) || (dp > 0.0 && dq < 0.0) {
            let t = dp / (dp - dq);
            let mut r = add(p, scale(sub(q, p), t));
            r[a] = x;
            lo.push(r);
            hi.push(r);
        }
    }
    (lo, hi)
}

/// Morceaux d'un polygone convexe, un par voxel traversé (plans i + ½).
fn voxel_pieces(poly: Vec<P>, out: &mut Vec<Vec<P>>) {
    let mut cur = vec![poly];
    for a in 0..3 {
        let mut next = Vec::new();
        for p in cur {
            let mn = p.iter().map(|v| v[a]).fold(f64::INFINITY, f64::min);
            let mx = p.iter().map(|v| v[a]).fold(f64::NEG_INFINITY, f64::max);
            let (i0, i1) = ((mn + 0.5).floor() as i64, (mx + 0.5).floor() as i64);
            let mut rest = p;
            for i in i0..i1 {
                let (lo, hi) = split(&rest, a, i as f64 + 0.5);
                if lo.len() >= 3 {
                    next.push(lo);
                }
                rest = hi;
                if rest.len() < 3 {
                    break;
                }
            }
            if rest.len() >= 3 {
                next.push(rest);
            }
        }
        cur = next;
    }
    out.extend(cur);
}

/// Niveau de grille : 2ᵏ fois plus fine que la veine (`scale`), nœud 0 au
/// point `origin` et domaine propre `core` (coins bas et haut), tous deux en
/// mailles de la veine ; mêmes conventions que `refine3d::Refinement`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Level {
    pub scale: f64,
    pub origin: P,
    pub core: [P; 2],
}

impl Level {
    /// La veine : tout l'espace.
    pub fn base() -> Level {
        Level { scale: 1.0, origin: [0.0; 3], core: [[f64::NEG_INFINITY; 3], [f64::INFINITY; 3]] }
    }
    /// Point en mailles de la veine → coordonnées de ce niveau.
    pub fn to_level(&self, p: P) -> P {
        [0, 1, 2].map(|a| self.scale * (p[a] - self.origin[a]))
    }
    fn contains(&self, p: P) -> bool {
        (0..3).all(|a| p[a] >= self.core[0][a] && p[a] <= self.core[1][a])
    }
}

/// Niveau propriétaire d'un point : le plus fin dont le domaine propre le
/// contient (le premier à finesse égale).
fn owner(levels: &[Level], p: P) -> usize {
    let mut best = 0;
    for (k, l) in levels.iter().enumerate() {
        if l.contains(p) && l.scale > levels[best].scale {
            best = k;
        }
    }
    best
}

/// Morceau de surface dans un voxel.
#[derive(Clone, Copy, Debug)]
pub struct Surfel {
    pub voxel: [i32; 3],
    /// Triangle d'origine.
    pub tri: u32,
    /// Aire (mailles du niveau au carré).
    pub area: f64,
    /// Normale unitaire, sortante du corps (vers le fluide).
    pub normal: P,
    /// Centre de gravité (coordonnées du niveau).
    pub center: P,
    /// Sommets : `verts[v0..v1]` de son niveau.
    pub v0: u32,
    pub v1: u32,
}

/// Part d'un parallélépipède balayé dans un voxel.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Overlap {
    pub voxel: [i32; 3],
    pub vol: f64,
}

/// Fraction solide des voxels coupés, rangés par colonne en x.
#[derive(Clone, Debug, Default)]
pub struct Fractions {
    /// Axe des colonnes (0 : x par défaut) ; (j, k) sont les deux autres
    /// coordonnées dans l'ordre circulaire.
    axis: usize,
    /// (j, k) → voxels coupés de la colonne, triés par i : (i, part solide
    /// de la face x = i + ½, volume solide du voxel).
    cols: HashMap<(i32, i32), Vec<(i32, f64, f64)>>,
}

impl Fractions {
    /// Construite à partir des morceaux de surface (voxel, vecteur aire,
    /// centre) d'un maillage fermé : dans une colonne, la part solide de la
    /// face aval vaut a(i + ½) = a(i − ½) − Σ Aₓ, et le volume solide
    /// V = Σ (cₓ − (i − ½)) Aₓ + a(i + ½) (divergence de (x − xᵢ₋½) eₓ).
    /// Colonnes selon l'axe `ax` (même calcul, coordonnées permutées) : un
    /// maillage périodique parallèle à x (canal incliné) n'a aucune
    /// composante Aₓ et demande des colonnes selon y ou z.
    fn build(pieces: &[([i32; 3], P, P)], ax: usize) -> Fractions {
        let (a1, a2) = ((ax + 1) % 3, (ax + 2) % 3);
        let mut by: HashMap<(i32, i32), Vec<(i32, f64, f64)>> = HashMap::new();
        for (v, a, c) in pieces {
            by.entry((v[a1], v[a2])).or_default().push((v[ax], a[ax], (c[ax] - (v[ax] as f64 - 0.5)) * a[ax]));
        }
        let cols = by
            .into_par_iter()
            .map(|(jk, mut list)| {
                list.sort_by_key(|e| e.0);
                let mut out: Vec<(i32, f64, f64)> = Vec::new();
                let mut a = 0.0;
                let mut k = 0;
                while k < list.len() {
                    let i = list[k].0;
                    let (mut ax, mut m) = (0.0, 0.0);
                    while k < list.len() && list[k].0 == i {
                        ax += list[k].1;
                        m += list[k].2;
                        k += 1;
                    }
                    a -= ax;
                    out.push((i, a, m + a));
                }
                (jk, out)
            })
            .collect();
        Fractions { axis: ax, cols }
    }

    /// Axe des colonnes.
    pub fn axis(&self) -> usize {
        self.axis
    }

    /// Volume solide du voxel (0 à 1) ; hors des voxels coupés, celui du
    /// dernier voxel coupé en amont dans la colonne (0 ou 1).
    pub fn solid(&self, v: [i32; 3]) -> f64 {
        let ax = self.axis;
        let Some(col) = self.cols.get(&(v[(ax + 1) % 3], v[(ax + 2) % 3])) else { return 0.0 };
        match col.binary_search_by_key(&v[ax], |e| e.0) {
            Ok(k) => col[k].2,
            Err(0) => 0.0,
            Err(k) => col[k - 1].1.round(),
        }
    }
    /// Volume fluide du voxel.
    pub fn fluid(&self, v: [i32; 3]) -> f64 {
        1.0 - self.solid(v)
    }
    /// Voxels coupés (voxel, volume solide).
    pub fn cut(&self) -> impl Iterator<Item = ([i32; 3], f64)> + '_ {
        let ax = self.axis;
        self.cols.iter().flat_map(move |(&(j, k), c)| {
            c.iter().map(move |e| {
                let mut v = [0; 3];
                v[ax] = e.0;
                v[(ax + 1) % 3] = j;
                v[(ax + 2) % 3] = k;
                (v, e.2)
            })
        })
    }
    /// Volume solide total : voxels coupés et voxels pleins entre eux.
    pub fn total_solid(&self) -> f64 {
        self.cols
            .values()
            .map(|c| {
                let mut v = c[0].2;
                for w in c.windows(2) {
                    v += w[1].2 + (w[1].0 - w[0].0 - 1) as f64 * w[0].1.round();
                }
                v
            })
            .sum()
    }
    /// Part solide de la face aval du dernier voxel de chaque colonne : nulle
    /// pour un maillage fermé (somme des écarts absolus).
    pub fn leak(&self) -> f64 {
        self.cols.values().filter_map(|c| c.last()).map(|e| e.1.abs()).sum()
    }
}

/// Surfels d'un niveau de grille.
#[derive(Clone, Debug)]
pub struct LevelSurfels {
    pub level: Level,
    pub surfels: Vec<Surfel>,
    pub verts: Vec<P>,
    /// Poids : `overlaps[start[9 s + p] .. start[9 s + p + 1]]` pour le
    /// surfel s et la paire de directions (2p + 1, 2p + 2).
    pub overlaps: Vec<Overlap>,
    pub start: Vec<u32>,
    /// Fractions solides des voxels coupés de ce niveau (toute la surface,
    /// pas seulement celle du domaine propre).
    pub fractions: Fractions,
}

impl LevelSurfels {
    /// Polygone du surfel `s`.
    pub fn poly(&self, s: usize) -> &[P] {
        let f = &self.surfels[s];
        &self.verts[f.v0 as usize..f.v1 as usize]
    }
    /// Volume balayé Pᵢ = |cᵢ·n| A (Δt = 1).
    pub fn swept(&self, s: usize, i: usize) -> f64 {
        let f = &self.surfels[s];
        dot(cvec(i), f.normal).abs() * f.area
    }
    /// Poids du surfel `s` pour la direction `i` (1 à 18) : parts de son
    /// parallélépipède dans chaque voxel (mêmes pour i et son opposée).
    pub fn weights(&self, s: usize, i: usize) -> &[Overlap] {
        assert!(i > 0);
        let k = 9 * s + (i - 1) / 2;
        &self.overlaps[self.start[k] as usize..self.start[k + 1] as usize]
    }
    /// Voxels touchés par les parallélépipèdes du surfel `s`, avec leur
    /// volume fluide (0 : voxel entièrement solide).
    /// Dispersion « miroir » du surfel `s`, paire `p` (`mirror_overlaps`) :
    /// (voxel rassemblé y, voxel d'arrivée x, volume).
    pub fn mirror(&self, s: usize, p: usize) -> Vec<([i32; 3], [i32; 3], f64)> {
        let f = &self.surfels[s];
        let (c1, c2) = (cvec(2 * p + 1), cvec(2 * p + 2));
        let c = if dot(c1, f.normal) > 0.0 { c1 } else { c2 };
        mirror_overlaps(self.poly(s), f.normal, c, self.weights(s, 2 * p + 1))
    }
    pub fn touched(&self, s: usize) -> Vec<([i32; 3], f64)> {
        let mut v: Vec<[i32; 3]> = (0..9).flat_map(|p| self.weights(s, 2 * p + 1).iter().map(|o| o.voxel)).collect();
        v.sort();
        v.dedup();
        v.into_iter().map(|x| (x, self.fractions.fluid(x))).collect()
    }
    /// Σ A n (nul pour un corps fermé dont tous les surfels sont à ce niveau).
    pub fn area_vector(&self) -> P {
        self.surfels.iter().fold([0.0; 3], |s, f| add(s, scale(f.normal, f.area)))
    }
    pub fn area(&self) -> f64 {
        self.surfels.iter().map(|f| f.area).sum()
    }
}

/// Polyèdre convexe donné par ses faces (polygones, sens indifférent).
type Poly3 = Vec<Vec<P>>;

const EPS: f64 = 1e-12;

fn unit_box(c: P) -> Poly3 {
    let v = |dx: f64, dy: f64, dz: f64| [c[0] + dx, c[1] + dy, c[2] + dz];
    let h = 0.5;
    vec![
        vec![v(-h, -h, -h), v(-h, h, -h), v(-h, h, h), v(-h, -h, h)],
        vec![v(h, -h, -h), v(h, h, -h), v(h, h, h), v(h, -h, h)],
        vec![v(-h, -h, -h), v(h, -h, -h), v(h, -h, h), v(-h, -h, h)],
        vec![v(-h, h, -h), v(h, h, -h), v(h, h, h), v(-h, h, h)],
        vec![v(-h, -h, -h), v(h, -h, -h), v(h, h, -h), v(-h, h, -h)],
        vec![v(-h, -h, h), v(h, -h, h), v(h, h, h), v(-h, h, h)],
    ]
}

/// Garde la partie n·x ≤ d (n unitaire) ; la coupe ferme le polyèdre.
fn clip(faces: Poly3, n: P, d: f64) -> Poly3 {
    let mut out: Poly3 = Vec::with_capacity(faces.len() + 1);
    let mut cap: Vec<P> = Vec::new();
    for f in faces {
        let mut g = Vec::with_capacity(f.len() + 1);
        for k in 0..f.len() {
            let (p, q) = (f[k], f[(k + 1) % f.len()]);
            let (dp, dq) = (dot(n, p) - d, dot(n, q) - d);
            if dp <= EPS {
                g.push(p);
                if dp >= -EPS {
                    cap.push(p);
                }
            }
            if (dp < -EPS && dq > EPS) || (dp > EPS && dq < -EPS) {
                let r = add(p, scale(sub(q, p), dp / (dp - dq)));
                g.push(r);
                cap.push(r);
            }
        }
        if g.len() >= 3 {
            out.push(g);
        }
    }
    if cap.len() >= 3 {
        // Points de la coupe rangés par angle autour de leur centre.
        let m = scale(cap.iter().fold([0.0; 3], |s, &p| add(s, p)), 1.0 / cap.len() as f64);
        let t = if n[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
        let u = cross(n, t);
        let w = cross(n, u);
        let mut a: Vec<(f64, P)> = cap.iter().map(|&p| (dot(sub(p, m), w).atan2(dot(sub(p, m), u)), p)).collect();
        a.sort_by(|x, y| x.0.total_cmp(&y.0));
        let mut c: Vec<P> = Vec::with_capacity(a.len());
        for (_, p) in a {
            if c.last().is_none_or(|&q| norm(sub(p, q)) > 1e-10) {
                c.push(p);
            }
        }
        while c.len() > 1 && norm(sub(c[0], *c.last().unwrap())) <= 1e-10 {
            c.pop();
        }
        if c.len() >= 3 {
            out.push(c);
        }
    }
    out
}

/// Volume d'un polyèdre convexe : pyramides depuis un point intérieur.
fn volume(faces: &Poly3) -> f64 {
    let (mut g, mut n) = ([0.0; 3], 0.0);
    for f in faces {
        for &p in f {
            g = add(g, p);
            n += 1.0;
        }
    }
    if n == 0.0 {
        return 0.0;
    }
    let g = scale(g, 1.0 / n);
    faces
        .iter()
        .filter(|f| f.len() >= 3)
        .map(|f| {
            let (a, _) = area_centroid(f);
            dot(sub(f[0], g), a).abs() / 3.0
        })
        .sum()
}

/// Parts du parallélépipède {x + t c : x ∈ polygone, t ∈ [0, 1]} dans les
/// voxels (`n` : normale unitaire du polygone, c·n > 0).
fn sweep_overlaps(poly: &[P], n: P, c: P) -> Vec<Overlap> {
    let cn = dot(c, n);
    if cn <= 0.0 {
        return Vec::new();
    }
    let planes = sweep_planes(poly, n, c);
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for &p in poly {
        for q in [p, add(p, c)] {
            for a in 0..3 {
                lo[a] = lo[a].min(q[a]);
                hi[a] = hi[a].max(q[a]);
            }
        }
    }
    let r0 = lo.map(|v| (v + 0.5).floor() as i32);
    let r1 = hi.map(|v| (v + 0.5).floor() as i32);
    let mut out = Vec::new();
    for k in r0[2]..=r1[2] {
        for j in r0[1]..=r1[1] {
            for i in r0[0]..=r1[0] {
                let ctr = [i as f64, j as f64, k as f64];
                // Plans qui coupent vraiment le voxel ; voxel exclu si l'un le laisse dehors.
                let mut cut: Vec<(P, f64)> = Vec::new();
                let mut empty = false;
                for &(m, d) in &planes {
                    let r = 0.5 * (m[0].abs() + m[1].abs() + m[2].abs());
                    let s = dot(m, ctr) - d;
                    if s - r >= -EPS {
                        empty = true;
                        break;
                    }
                    if s + r > EPS {
                        cut.push((m, d));
                    }
                }
                if empty {
                    continue;
                }
                let vol = if cut.is_empty() {
                    1.0
                } else {
                    let mut f = unit_box(ctr);
                    for (m, d) in cut {
                        f = clip(f, m, d);
                        if f.len() < 4 {
                            break;
                        }
                    }
                    if f.len() < 4 { 0.0 } else { volume(&f) }
                };
                if vol > 1e-13 {
                    out.push(Overlap { voxel: [i, j, k], vol });
                }
            }
        }
    }
    out
}

/// Demi-espaces m·x ≤ d (m unitaire) du parallélépipède
/// {x + t c : x ∈ polygone, t ∈ [0, 1]} (c·n > 0).
fn sweep_planes(poly: &[P], n: P, c: P) -> Vec<(P, f64)> {
    let cn = dot(c, n);
    let g = scale(poly.iter().fold([0.0; 3], |s, &p| add(s, p)), 1.0 / poly.len() as f64);
    let mut planes: Vec<(P, f64)> = vec![(scale(n, -1.0), -dot(n, poly[0])), (n, dot(n, poly[0]) + cn)];
    for k in 0..poly.len() {
        let (a, b) = (poly[k], poly[(k + 1) % poly.len()]);
        let m = cross(sub(b, a), c);
        let l = norm(m);
        if l <= 1e-14 * (1.0 + norm(sub(b, a))) {
            continue;
        }
        let mut m = scale(m, 1.0 / l);
        if dot(m, sub(g, a)) > 0.0 {
            m = scale(m, -1.0);
        }
        planes.push((m, dot(m, a)));
    }
    planes
}

/// Rebond à sa place (dispersion « miroir ») : une particule entrante
/// (direction −c) partie de x à la hauteur η = n·(x − x₀) au-dessus du
/// surfel le frappe au temps η/h (h = c·n) et revient en sens inverse : au
/// pas suivant elle est en φ(x) = x + (1 − 2η/h) c, à la hauteur h − η.
/// φ est une involution affine du parallélépipède. Rend, pour chaque couple
/// de voxels (y, x) de ses parts, le volume des points de la part de y que
/// φ envoie dans x : la population rassemblée dans y ressort dans x avec ce
/// poids, au lieu d'être répartie uniformément sur tout le parallélépipède
/// (mélange sur une maille à chaque pas : diffusion numérique d'ordre 1).
/// Pour chaque y, la somme sur x redonne sa part ; pour chaque x aussi
/// (φ conserve le volume).
/// Comme `clip`, sans rien faire si le plan laisse tout le polyèdre du bon
/// côté (un plan confondu avec une face la dédoublerait), et vide s'il le
/// laisse tout dehors.
fn clip_safe(faces: Poly3, n: P, d: f64) -> Poly3 {
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for f in &faces {
        for &p in f {
            let v = dot(n, p) - d;
            lo = lo.min(v);
            hi = hi.max(v);
        }
    }
    if hi <= 1e-10 {
        faces
    } else if lo >= -1e-10 {
        Vec::new()
    } else {
        clip(faces, n, d)
    }
}

pub fn mirror_overlaps(poly: &[P], n: P, c: P, parts: &[Overlap]) -> Vec<([i32; 3], [i32; 3], f64)> {
    let h = dot(c, n);
    if h <= 0.0 || parts.is_empty() {
        return Vec::new();
    }
    let planes = sweep_planes(poly, n, c);
    let d0 = dot(n, poly[0]);
    // φ(x)_a = g_a·x + k_a, g_a = e_a − (2 c_a / h) n, k_a = c_a (1 + 2 d₀ / h).
    let g: [P; 3] = [0, 1, 2].map(|a| {
        let mut e = [0.0; 3];
        e[a] = 1.0;
        sub(e, scale(n, 2.0 * c[a] / h))
    });
    let k: [f64; 3] = [0, 1, 2].map(|a| c[a] * (1.0 + 2.0 * d0 / h));
    let mut out = Vec::new();
    for py in parts {
        let ctr = py.voxel.map(|v| v as f64);
        let mut fy = unit_box(ctr);
        for &(m, d) in &planes {
            fy = clip_safe(fy, m, d);
            if fy.len() < 4 {
                break;
            }
        }
        if fy.len() < 4 {
            continue;
        }
        for px in parts {
            let mut f = fy.clone();
            for a in 0..3 {
                let l = norm(g[a]);
                let (m, xa) = (scale(g[a], 1.0 / l), px.voxel[a] as f64);
                // g_a·x + k_a ≤ x_a + ½ et ≥ x_a − ½.
                f = clip_safe(f, m, (xa + 0.5 - k[a]) / l);
                if f.len() < 4 {
                    break;
                }
                f = clip_safe(f, scale(m, -1.0), -(xa - 0.5 - k[a]) / l);
                if f.len() < 4 {
                    break;
                }
            }
            let v = if f.len() < 4 { 0.0 } else { volume(&f) };
            if v > 1e-13 {
                out.push((py.voxel, px.voxel, v));
            }
        }
    }
    out
}


/// Surfels du maillage `tris` (coordonnées de la veine, normales sortantes)
/// à chaque niveau de `levels` (le premier est en général `Level::base()`).
pub fn build(tris: &[Tri], levels: &[Level]) -> Vec<LevelSurfels> {
    build_axis(tris, levels, 0)
}

/// Comme `build`, fractions solides intégrées en colonnes selon l'axe `ax`.
pub fn build_axis(tris: &[Tri], levels: &[Level], ax: usize) -> Vec<LevelSurfels> {
    // Plans des domaines propres finis (coordonnées de la veine).
    let mut planes: Vec<(usize, f64)> = Vec::new();
    for l in levels {
        for s in 0..2 {
            for a in 0..3 {
                if l.core[s][a].is_finite() {
                    planes.push((a, l.core[s][a]));
                }
            }
        }
    }
    // Morceaux de chaque triangle par domaine propre : (triangle, niveau, polygone).
    let parts: Vec<(u32, usize, Vec<P>)> = tris
        .par_iter()
        .enumerate()
        .flat_map_iter(|(t, tri)| {
            let mut cur = vec![tri.to_vec()];
            for &(a, x) in &planes {
                let mut next = Vec::with_capacity(cur.len() + 1);
                for p in cur {
                    let mn = p.iter().map(|v| v[a]).fold(f64::INFINITY, f64::min);
                    let mx = p.iter().map(|v| v[a]).fold(f64::NEG_INFINITY, f64::max);
                    if mn < x && mx > x {
                        let (lo, hi) = split(&p, a, x);
                        next.extend([lo, hi].into_iter().filter(|q| q.len() >= 3));
                    } else {
                        next.push(p);
                    }
                }
                cur = next;
            }
            cur.into_iter()
                .filter_map(|p| {
                    let (a, c) = area_centroid(&p);
                    (norm(a) > 0.0).then(|| (t as u32, owner(levels, c), p))
                })
                .collect::<Vec<_>>()
        })
        .collect();
    levels.iter().enumerate().map(|(li, l)| build_level(tris, &parts, li, *l, ax)).collect()
}

fn build_level(tris: &[Tri], parts: &[(u32, usize, Vec<P>)], li: usize, l: Level, ax: usize) -> LevelSurfels {
    // Fenêtre (y, z) des colonnes utiles aux fractions : domaine propre et
    // 3 voxels de marge (portée des parallélépipèdes : √2 maille).
    let win: [[f64; 2]; 3] = [0, 1, 2].map(|a| {
        let (u, v) = (l.scale * (l.core[0][a] - l.origin[a]), l.scale * (l.core[1][a] - l.origin[a]));
        [u - 3.0, v + 3.0]
    });
    let normals: Vec<P> = tris
        .iter()
        .map(|t| {
            let a = cross(sub(t[1], t[0]), sub(t[2], t[0]));
            let n = norm(a);
            if n > 0.0 { scale(a, 1.0 / n) } else { [0.0; 3] }
        })
        .collect();
    // Morceaux par voxel : (voxel, triangle, à ce niveau ?, polygone).
    let pieces: Vec<([i32; 3], u32, bool, Vec<P>)> = parts
        .par_iter()
        .filter(|(_, _, p)| {
            (1..3).all(|a| {
                let mn = p.iter().map(|v| l.scale * (v[a] - l.origin[a])).fold(f64::INFINITY, f64::min);
                let mx = p.iter().map(|v| l.scale * (v[a] - l.origin[a])).fold(f64::NEG_INFINITY, f64::max);
                mx >= win[a][0] && mn <= win[a][1]
            })
        })
        .flat_map_iter(|(t, o, p)| {
            let mut out = Vec::new();
            voxel_pieces(p.iter().map(|&v| l.to_level(v)).collect(), &mut out);
            out.into_iter().filter_map(move |q| {
                let (_, c) = area_centroid(&q);
                let v = c.map(|x| (x + 0.5).floor() as i32);
                Some((v, *t, *o == li, q))
            })
        })
        .collect();
    let frac: Vec<([i32; 3], P, P)> = pieces
        .par_iter()
        .filter(|(v, ..)| (1..3).all(|a| (v[a] as f64) >= win[a][0] && (v[a] as f64) <= win[a][1]))
        .map(|(v, _, _, q)| {
            let (a, c) = area_centroid(q);
            (*v, a, c)
        })
        .collect();
    let fractions = Fractions::build(&frac, ax);
    let mut surfels = Vec::new();
    let mut verts = Vec::new();
    for (v, t, own, q) in pieces {
        if !own {
            continue;
        }
        let (a, c) = area_centroid(&q);
        let n = normals[t as usize];
        let area = dot(a, n);
        if !(area > 1e-14) {
            continue;
        }
        let v0 = verts.len() as u32;
        verts.extend_from_slice(&q);
        surfels.push(Surfel { voxel: v, tri: t, area, normal: n, center: c, v0, v1: verts.len() as u32 });
    }
    let lists: Vec<Vec<Vec<Overlap>>> = surfels
        .par_iter()
        .map(|f| {
            let poly = &verts[f.v0 as usize..f.v1 as usize];
            (0..9)
                .map(|p| {
                    let (c1, c2) = (cvec(2 * p + 1), cvec(2 * p + 2));
                    let c = if dot(c1, f.normal) > 0.0 { c1 } else { c2 };
                    sweep_overlaps(poly, f.normal, c)
                })
                .collect()
        })
        .collect();
    let mut overlaps = Vec::new();
    let mut start = vec![0u32];
    for s in lists {
        for p in s {
            overlaps.extend(p);
            start.push(overlaps.len() as u32);
        }
    }
    LevelSurfels { level: l, surfels, verts, overlaps, start, fractions }
}

/// Bilans d'un niveau : sommes des volumes balayés, écart à la partition, et
/// volume balayé en trop par rapport au volume fluide des voxels.
#[derive(Clone, Copy, Debug, Default)]
pub struct Balance {
    /// Σ Pᵢ sur les surfels et les 18 directions.
    pub swept: f64,
    /// max |Σ poids − Pᵢ| sur les surfels et directions.
    pub partition_err: f64,
    /// Σ sur voxels et directions de max(0, Σ poids − volume fluide) : volume
    /// balayé qui tombe dans le solide ou se recouvre.
    pub excess: f64,
    /// Volume balayé dans des voxels entièrement solides.
    pub in_solid: f64,
}

impl LevelSurfels {
    pub fn balance(&self) -> Balance {
        let mut b = Balance::default();
        for i in 1..19 {
            let mut per: HashMap<[i32; 3], f64> = HashMap::new();
            for s in 0..self.surfels.len() {
                let p = self.swept(s, i);
                b.swept += p;
                let out = dot(cvec(i), self.surfels[s].normal) > 0.0;
                let w = self.weights(s, i);
                let tot: f64 = w.iter().map(|o| o.vol).sum();
                b.partition_err = b.partition_err.max((tot - p).abs());
                if out {
                    for o in w {
                        *per.entry(o.voxel).or_default() += o.vol;
                    }
                }
            }
            for (v, w) in per {
                let f = self.fractions.fluid(v);
                b.excess += (w - f).max(0.0);
                if f <= 1e-12 {
                    b.in_solid += w;
                }
            }
        }
        b
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Rebond miroir : paroi plane z = z_w, voxel coupé k0 de volume fluide
    /// V = ½ (paroi au nœud) ; direction normale : la moitié basse de k0 reçoit
    /// la moitié haute du parallélépipède (voxel k0 + 1) et inversement ;
    /// direction oblique : sommes par y et par x égales aux parts.
    #[test]
    fn miroir_paroi_plane() {
        for &zw in &[0.0, 0.2, -0.35] {
            let poly = vec![[-0.5, -0.5, zw], [0.5, -0.5, zw], [0.5, 0.5, zw], [-0.5, 0.5, zw]];
            let n = [0.0, 0.0, 1.0];
            for c in [[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, -1.0, 1.0]] {
                let parts = sweep_overlaps(&poly, n, c);
                let m = mirror_overlaps(&poly, n, c, &parts);
                for o in &parts {
                    let sy: f64 = m.iter().filter(|e| e.0 == o.voxel).map(|e| e.2).sum();
                    let sx: f64 = m.iter().filter(|e| e.1 == o.voxel).map(|e| e.2).sum();
                    assert!((sy - o.vol).abs() < 1e-10 && (sx - o.vol).abs() < 1e-10, "{zw} {c:?} {o:?} {sy} {sx}");
                }
                if zw == 0.0 && c == [0.0, 0.0, 1.0] {
                    // y = k0 + 1 (z = 1) → x = k0 (z = 0) : volume ½.
                    let v: f64 = m.iter().filter(|e| e.0[2] == 1 && e.1[2] == 0).map(|e| e.2).sum();
                    assert!((v - 0.5).abs() < 1e-10, "{m:?}");
                }
            }
        }
    }

    fn tri_area(t: &Tri) -> f64 {
        0.5 * norm(cross(sub(t[1], t[0]), sub(t[2], t[0])))
    }
    fn mesh_volume(tris: &[Tri]) -> f64 {
        tris.iter().map(|t| dot(t[0], cross(t[1], t[2])) / 6.0).sum()
    }
    /// Cube d'arête `a`, centre `c`, tourné de `th` autour de (1, 2, 3).
    fn cube(a: f64, c: P, th: f64) -> Vec<Tri> {
        let k = scale([1.0, 2.0, 3.0], 1.0 / 14f64.sqrt());
        let rot = |v: P| {
            let (s, co) = th.sin_cos();
            add(add(scale(v, co), scale(cross(k, v), s)), scale(k, dot(k, v) * (1.0 - co)))
        };
        let h = a / 2.0;
        let p = |x: f64, y: f64, z: f64| add(rot([x * h, y * h, z * h]), c);
        let mut t = Vec::new();
        // Faces : sommets dans le sens direct vu de l'extérieur.
        let quads = [
            [p(-1., -1., -1.), p(-1., -1., 1.), p(-1., 1., 1.), p(-1., 1., -1.)],
            [p(1., -1., -1.), p(1., 1., -1.), p(1., 1., 1.), p(1., -1., 1.)],
            [p(-1., -1., -1.), p(1., -1., -1.), p(1., -1., 1.), p(-1., -1., 1.)],
            [p(-1., 1., -1.), p(-1., 1., 1.), p(1., 1., 1.), p(1., 1., -1.)],
            [p(-1., -1., -1.), p(-1., 1., -1.), p(1., 1., -1.), p(1., -1., -1.)],
            [p(-1., -1., 1.), p(1., -1., 1.), p(1., 1., 1.), p(-1., 1., 1.)],
        ];
        for q in quads {
            t.push([q[0], q[1], q[2]]);
            t.push([q[0], q[2], q[3]]);
        }
        t
    }
    /// Volume solide des voxels de [0, n)³.
    fn solid_sum(l: &LevelSurfels, n: i32) -> f64 {
        let mut v = 0.0;
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    v += l.fractions.solid([i, j, k]);
                }
            }
        }
        v
    }
    fn check_closed(tris: &[Tri], levels: &[Level]) -> Vec<LevelSurfels> {
        let ls = build(tris, levels);
        let a0: f64 = tris.iter().map(tri_area).sum();
        let a: f64 = ls.iter().map(|l| l.area() / (l.level.scale * l.level.scale)).sum();
        assert!((a - a0).abs() < 1e-10 * a0, "aire {a} contre {a0}");
        let an = ls.iter().fold([0.0; 3], |s, l| add(s, scale(l.area_vector(), 1.0 / (l.level.scale * l.level.scale))));
        assert!(norm(an) < 1e-10 * a0, "Σ A n = {an:?}");
        ls
    }

    #[test]
    fn cube_area_volume_fractions() {
        // Cube aligné : volume solide de chaque voxel = produit des recouvrements 1D.
        let (a, c) = (5.3, [10.37, 11.71, 9.13]);
        let tris = cube(a, c, 0.0);
        let ls = check_closed(&tris, &[Level::base()]);
        let l = &ls[0];
        assert!((l.area() - 6.0 * a * a).abs() < 1e-10);
        let ov = |i: i32, ax: usize| {
            let (u, v) = (c[ax] - a / 2.0, c[ax] + a / 2.0);
            ((i as f64 + 0.5).min(v) - (i as f64 - 0.5).max(u)).max(0.0)
        };
        for (v, s) in l.fractions.cut() {
            let e = ov(v[0], 0) * ov(v[1], 1) * ov(v[2], 2);
            assert!((s - e).abs() < 1e-12, "voxel {v:?} : {s} contre {e}");
        }
        for i in 0..25 {
            for j in 0..25 {
                for k in 0..25 {
                    let e = ov(i, 0) * ov(j, 1) * ov(k, 2);
                    assert!((l.fractions.solid([i, j, k]) - e).abs() < 1e-12);
                }
            }
        }
        let vs = solid_sum(l, 25);
        assert!((vs - a * a * a).abs() < 1e-10, "volume {vs}");
        assert!(l.fractions.leak() < 1e-12);
        // Direction +x : le parallélépipède de la face +x est la tranche
        // [x₁, x₁ + 1] × face ; poids par voxel exacts.
        let mut per: HashMap<[i32; 3], f64> = HashMap::new();
        for s in 0..l.surfels.len() {
            if l.surfels[s].normal[0] > 0.5 {
                for o in l.weights(s, 1) {
                    *per.entry(o.voxel).or_default() += o.vol;
                }
            }
        }
        let x1 = c[0] + a / 2.0;
        for (v, w) in &per {
            let ex = ((v[0] as f64 + 0.5).min(x1 + 1.0) - (v[0] as f64 - 0.5).max(x1)).max(0.0);
            let e = ex * ov(v[1], 1) * ov(v[2], 2);
            assert!((w - e).abs() < 1e-11, "voxel {v:?} : {w} contre {e}");
        }
        let b = l.balance();
        assert!(b.partition_err < 1e-11, "{b:?}");
        assert!(b.excess < 1e-10 && b.in_solid == 0.0, "{b:?}");
    }

    #[test]
    fn rotated_cube() {
        let a = 4.7;
        let tris = cube(a, [8.2, 7.9, 8.6], 0.61);
        let ls = check_closed(&tris, &[Level::base()]);
        let l = &ls[0];
        let v = solid_sum(l, 20);
        assert!((v - a * a * a).abs() < 1e-9, "volume {v}");
        let b = l.balance();
        // Corps convexe : les parallélépipèdes d'une direction sortante sont
        // disjoints et dans le fluide.
        assert!(b.partition_err < 1e-11, "{b:?}");
        assert!(b.excess < 1e-9 && b.in_solid == 0.0, "{b:?}");
        // Σ Pᵢ = Σ |cᵢ·n| A, sommé sur les directions.
        let e: f64 = (0..6).map(|f| 2.0 * a * a * {
            let n = cross(sub(tris[2 * f][1], tris[2 * f][0]), sub(tris[2 * f][2], tris[2 * f][0]));
            let n = scale(n, 1.0 / norm(n));
            (1..19).map(|i| dot(cvec(i), n).abs()).sum::<f64>() / 2.0
        }).sum();
        assert!((b.swept - e).abs() < 1e-9 * e, "{} contre {e}", b.swept);
    }

    #[test]
    fn sphere_area_volume() {
        let r = 6.3;
        let tris = crate::stl::icosphere(r, [11.3, 12.1, 10.7], 3);
        let a0: f64 = tris.iter().map(tri_area).sum();
        let v0 = mesh_volume(&tris);
        // Polyèdre proche de la sphère exacte.
        assert!((a0 / (4.0 * std::f64::consts::PI * r * r) - 1.0).abs() < 0.02);
        let ls = check_closed(&tris, &[Level::base()]);
        let l = &ls[0];
        let v = solid_sum(l, 25);
        assert!((v - v0).abs() < 1e-9 * v0, "volume {v} contre {v0}");
        assert!(l.fractions.leak() < 1e-10);
        for (_, s) in l.fractions.cut() {
            assert!((-1e-12..=1.0 + 1e-12).contains(&s));
        }
        let b = l.balance();
        assert!(b.partition_err < 1e-11, "{b:?}");
        assert!(b.excess < 1e-9 && b.in_solid == 0.0, "{b:?}");
    }

    #[test]
    fn nested_levels() {
        // Sphère coupée par le domaine propre d'un niveau ×2 puis ×4 : l'aire
        // et Σ A n se conservent entre niveaux.
        let tris = crate::stl::icosphere(5.0, [10.2, 10.4, 10.1], 3);
        let l1 = Level { scale: 2.0, origin: [3.0, 2.0, 4.0], core: [[9.1, 4.0, 4.0], [20.0, 16.0, 16.0]] };
        let l2 = Level { scale: 4.0, origin: [8.0, 6.0, 7.0], core: [[12.3, 6.0, 9.6], [20.0, 13.0, 16.0]] };
        let ls = check_closed(&tris, &[Level::base(), l1, l2]);
        assert!(ls.iter().all(|l| !l.surfels.is_empty()));
        // Les fractions solides de chaque niveau donnent le volume de la
        // partie de sphère dans une boîte de ce niveau.
        let v0 = mesh_volume(&tris);
        let l = &ls[1];
        let v = solid_sum(l, 50);
        assert!((v / 8.0 - v0).abs() < 1e-9 * v0, "volume ×2 {} contre {v0}", v / 8.0);
        for l in &ls {
            let b = l.balance();
            assert!(b.partition_err < 1e-11, "{b:?}");
        }
    }

    #[test]
    fn plane_through_nodes() {
        // Faces posées exactement sur des plans de voxels et des nœuds.
        let tris = cube(4.0, [10.5, 10.5, 10.0], 0.0);
        let ls = check_closed(&tris, &[Level::base()]);
        let b = ls[0].balance();
        assert!(b.partition_err < 1e-11 && b.excess < 1e-10, "{b:?}");
        let v = solid_sum(&ls[0], 20);
        assert!((v - 64.0).abs() < 1e-10, "volume {v}");
    }
}
