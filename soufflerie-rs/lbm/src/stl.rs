//! Maquettes 3D décrites par un maillage STL (binaire ou ASCII) : lecture,
//! mise en place dans la veine et requêtes géométriques accélérées par des
//! grilles de rangement (appartenance par lancer de rayon selon x, distance
//! à la surface, maître-couple).

use std::path::Path;

pub type Tri = [[f64; 3]; 3];

/// Axe orienté du repère du fichier.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Axis {
    Xp,
    Xm,
    Yp,
    Ym,
    Zp,
    Zm,
}

impl Axis {
    pub const ALL: [Axis; 6] = [Axis::Xp, Axis::Xm, Axis::Yp, Axis::Ym, Axis::Zp, Axis::Zm];

    /// « x », « +x », « -x »…
    pub fn parse(s: &str) -> Option<Axis> {
        Some(match s.trim().to_ascii_lowercase().as_str() {
            "x" | "+x" => Axis::Xp,
            "-x" => Axis::Xm,
            "y" | "+y" => Axis::Yp,
            "-y" => Axis::Ym,
            "z" | "+z" => Axis::Zp,
            "-z" => Axis::Zm,
            _ => return None,
        })
    }
    pub fn label(self) -> &'static str {
        match self {
            Axis::Xp => "+x",
            Axis::Xm => "−x",
            Axis::Yp => "+y",
            Axis::Ym => "−y",
            Axis::Zp => "+z",
            Axis::Zm => "−z",
        }
    }
    pub fn vec(self) -> [f64; 3] {
        match self {
            Axis::Xp => [1.0, 0.0, 0.0],
            Axis::Xm => [-1.0, 0.0, 0.0],
            Axis::Yp => [0.0, 1.0, 0.0],
            Axis::Ym => [0.0, -1.0, 0.0],
            Axis::Zp => [0.0, 0.0, 1.0],
            Axis::Zm => [0.0, 0.0, -1.0],
        }
    }
}

/// Unité des coordonnées du fichier.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Unit {
    Mm,
    Cm,
    M,
    Inch,
}

impl Unit {
    pub const ALL: [Unit; 4] = [Unit::Mm, Unit::Cm, Unit::M, Unit::Inch];

    pub fn parse(s: &str) -> Option<Unit> {
        Some(match s.trim().to_ascii_lowercase().as_str() {
            "mm" => Unit::Mm,
            "cm" => Unit::Cm,
            "m" => Unit::M,
            "in" | "inch" | "po" => Unit::Inch,
            _ => return None,
        })
    }
    pub fn label(self) -> &'static str {
        match self {
            Unit::Mm => "mm",
            Unit::Cm => "cm",
            Unit::M => "m",
            Unit::Inch => "in",
        }
    }
    /// Mètres par unité.
    pub fn meters(self) -> f64 {
        match self {
            Unit::Mm => 1e-3,
            Unit::Cm => 1e-2,
            Unit::M => 1.0,
            Unit::Inch => 0.0254,
        }
    }
}

/// Mise en place d'une maquette maillée dans la veine.
#[derive(Clone, Debug, PartialEq)]
pub struct MeshSetup {
    /// Axe du fichier dans le sens de l'écoulement (de l'avant vers l'arrière).
    pub flow: Axis,
    /// Axe du fichier vers le haut.
    pub up: Axis,
    pub unit: Unit,
    /// Incidence (°) : cabrage, avant vers le haut.
    pub alpha_deg: f64,
    /// Dérapage (°) : vent venant de la droite de la maquette.
    pub beta_deg: f64,
    /// Longueur de référence (unités du fichier) ; None : longueur hors tout
    /// selon l'écoulement.
    pub ref_length: Option<f64>,
    /// Surface de référence (unités du fichier au carré) ; None : maître-couple.
    pub ref_area: Option<f64>,
    /// Facteur sur les marges de la veine autour de la maquette (1 : 3,5 M en
    /// amont, 7,5 M en aval, 2,5 M de côté). Une aile très allongée tient dans
    /// une veine plus serrée : le blocage affiché dit si c'est raisonnable.
    pub margin: f64,
    /// Sol mobile : garde au sol (unités du fichier) ; None : veine libre.
    pub ground: Option<f64>,
    /// Décalage de la maquette dans la maille (mailles, x y z), ajouté au
    /// calage par défaut : sert à mesurer la sensibilité au placement.
    pub offset: [f64; 3],
}

impl Default for MeshSetup {
    fn default() -> Self {
        MeshSetup {
            flow: Axis::Xp,
            up: Axis::Zp,
            unit: Unit::Mm,
            alpha_deg: 0.0,
            beta_deg: 0.0,
            ref_length: None,
            ref_area: None,
            margin: 1.0,
            ground: None,
            offset: [0.0; 3],
        }
    }
}

/// Grandeurs de référence d'une maquette maillée placée.
#[derive(Clone, Debug)]
pub struct MeshInfo {
    pub triangles: usize,
    /// Longueur hors tout selon l'écoulement, envergure, hauteur (m), maquette placée.
    pub size_m: [f64; 3],
    /// Longueur (m) et surface (m²) de référence.
    pub ref_length_m: f64,
    pub ref_area_m2: f64,
    /// Maître-couple (m²).
    pub frontal_m2: f64,
    /// Rapport du maître-couple à la section de la veine.
    pub blockage: f64,
    /// Volume enclos (m³).
    pub volume_m3: f64,
    /// Défauts du maillage.
    pub diag: MeshDiag,
}

impl MeshInfo {
    /// Part du volume perdue à la voxelisation (`solid` mailles solides, `l`
    /// mailles par longueur de référence) : parties plus fines qu'une maille.
    /// None si le volume est inconnu.
    pub fn volume_lost(&self, solid: usize, l: usize) -> Option<f64> {
        let h = self.ref_length_m / l as f64;
        let v = self.volume_m3 / (h * h * h);
        (v > 0.0).then(|| 1.0 - solid as f64 / v)
    }
}

/// Triangles lus dans un fichier STL (coordonnées du fichier).
pub fn read_stl(path: &Path) -> Result<Vec<Tri>, String> {
    let data = std::fs::read(path).map_err(|e| format!("lecture de {} impossible : {e}", path.display()))?;
    parse_stl(&data)
}

/// STL binaire (en-tête de 80 octets, nombre de triangles, 50 octets par
/// triangle) ou ASCII (`facet … vertex x y z …`).
pub fn parse_stl(data: &[u8]) -> Result<Vec<Tri>, String> {
    if data.len() >= 84 {
        let n = u32::from_le_bytes([data[80], data[81], data[82], data[83]]) as usize;
        // Taille exacte ; ou octets en trop après les triangles, sauf texte ASCII.
        let ascii = data.starts_with(b"solid") && std::str::from_utf8(data).is_ok();
        if data.len() == 84 + 50 * n || (n > 0 && data.len() > 84 + 50 * n && !ascii) {
            let mut tris = Vec::with_capacity(n);
            for k in 0..n {
                let o = 84 + 50 * k + 12;
                let f = |i: usize| {
                    let b = &data[o + 4 * i..o + 4 * i + 4];
                    f32::from_le_bytes([b[0], b[1], b[2], b[3]]) as f64
                };
                tris.push([[f(0), f(1), f(2)], [f(3), f(4), f(5)], [f(6), f(7), f(8)]]);
            }
            return check(tris);
        }
    }
    let text = std::str::from_utf8(data).map_err(|_| "fichier STL illisible (ni binaire ni ASCII)".to_string())?;
    let mut tris = Vec::new();
    let mut v: Vec<[f64; 3]> = Vec::new();
    let mut words = text.split_ascii_whitespace();
    while let Some(w) = words.next() {
        if w.eq_ignore_ascii_case("vertex") {
            let mut p = [0.0; 3];
            for c in &mut p {
                *c = words
                    .next()
                    .and_then(|s| s.parse().ok())
                    .ok_or_else(|| "sommet STL mal formé".to_string())?;
            }
            v.push(p);
        } else if w.eq_ignore_ascii_case("endfacet") {
            if v.len() == 3 {
                tris.push([v[0], v[1], v[2]]);
            }
            v.clear();
        }
    }
    check(tris)
}

fn check(tris: Vec<Tri>) -> Result<Vec<Tri>, String> {
    let tris: Vec<Tri> = tris
        .into_iter()
        .filter(|t| t.iter().flatten().all(|c| c.is_finite()))
        .collect();
    if tris.is_empty() {
        return Err("aucun triangle dans le fichier STL".into());
    }
    Ok(tris)
}

pub fn bbox(tris: &[Tri]) -> ([f64; 3], [f64; 3]) {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in tris.iter().flatten() {
        for a in 0..3 {
            lo[a] = lo[a].min(p[a]);
            hi[a] = hi[a].max(p[a]);
        }
    }
    (lo, hi)
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Distance au carré d'un point à un triangle (Ericson, Real-Time Collision Detection, 5.1.5).
fn tri_dist2(p: [f64; 3], t: &Tri) -> f64 {
    let [a, b, c] = *t;
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = sub(p, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    let q = if d1 <= 0.0 && d2 <= 0.0 {
        a
    } else {
        let bp = sub(p, b);
        let d3 = dot(ab, bp);
        let d4 = dot(ac, bp);
        if d3 >= 0.0 && d4 <= d3 {
            b
        } else {
            let vc = d1 * d4 - d3 * d2;
            if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
                let v = d1 / (d1 - d3);
                [a[0] + v * ab[0], a[1] + v * ab[1], a[2] + v * ab[2]]
            } else {
                let cp = sub(p, c);
                let d5 = dot(ab, cp);
                let d6 = dot(ac, cp);
                if d6 >= 0.0 && d5 <= d6 {
                    c
                } else {
                    let vb = d5 * d2 - d1 * d6;
                    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
                        let w = d2 / (d2 - d6);
                        [a[0] + w * ac[0], a[1] + w * ac[1], a[2] + w * ac[2]]
                    } else {
                        let va = d3 * d6 - d5 * d4;
                        if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
                            let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
                            let bc = sub(c, b);
                            [b[0] + w * bc[0], b[1] + w * bc[1], b[2] + w * bc[2]]
                        } else {
                            let denom = 1.0 / (va + vb + vc);
                            let v = vb * denom;
                            let w = vc * denom;
                            [
                                a[0] + ab[0] * v + ac[0] * w,
                                a[1] + ab[1] * v + ac[1] * w,
                                a[2] + ab[2] * v + ac[2] * w,
                            ]
                        }
                    }
                }
            }
        }
    };
    let d = sub(p, q);
    dot(d, d)
}

/// Décalage des rayons d'appartenance : évite de passer exactement par une
/// arête ou un sommet (coordonnées entières des mailles).
const JU: f64 = 1.37e-5;
const JV: f64 = 2.71e-5;

/// Défauts d'un maillage, après nettoyage.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MeshDiag {
    /// Triangles retirés : dégénérés (aire nulle) et en double.
    pub degenerate: usize,
    pub duplicates: usize,
    /// Arêtes bordant un seul triangle (trous), plus de deux (non variété).
    pub open_edges: usize,
    pub nonmanifold_edges: usize,
    /// Arêtes dont les deux triangles sont orientés en sens contraires.
    pub flipped_edges: usize,
    /// Volume enclos (unités du fichier au cube), par le théorème de la divergence.
    pub volume: f64,
}

impl MeshDiag {
    /// Volume fermé et bien orienté.
    pub fn watertight(&self) -> bool {
        self.open_edges == 0 && self.nonmanifold_edges == 0
    }

    /// Résumé des défauts (vide si aucun).
    pub fn summary(&self) -> String {
        let mut v = Vec::new();
        if self.open_edges > 0 {
            v.push(format!("{} arêtes ouvertes (trous)", self.open_edges));
        }
        if self.nonmanifold_edges > 0 {
            v.push(format!("{} arêtes partagées par plus de deux faces", self.nonmanifold_edges));
        }
        if self.flipped_edges > 0 {
            v.push(format!("{} arêtes entre faces d'orientations opposées", self.flipped_edges));
        }
        if self.degenerate + self.duplicates > 0 {
            v.push(format!(
                "{} triangles dégénérés et {} en double retirés",
                self.degenerate, self.duplicates
            ));
        }
        v.join(" · ")
    }
}

/// Retire les triangles dégénérés et en double (sommets soudés à 10⁻⁷ de la
/// diagonale près) et dresse le bilan des arêtes.
pub fn clean(tris: Vec<Tri>) -> (Vec<Tri>, MeshDiag) {
    use std::collections::HashMap;
    let (lo, hi) = bbox(&tris);
    let diag = dot(sub(hi, lo), sub(hi, lo)).sqrt().max(1e-30);
    let eps = 1e-7 * diag;
    let mut ids: HashMap<[i64; 3], u32> = HashMap::new();
    let mut id = |p: [f64; 3]| -> u32 {
        let k = p.map(|c| (c / eps).round() as i64);
        let n = ids.len() as u32;
        *ids.entry(k).or_insert(n)
    };
    let mut d = MeshDiag::default();
    let mut seen = std::collections::HashSet::new();
    let mut edges: HashMap<(u32, u32), (u32, i32)> = HashMap::new();
    let mut out = Vec::with_capacity(tris.len());
    for t in tris {
        let v = t.map(&mut id);
        let ab = sub(t[1], t[0]);
        let ac = sub(t[2], t[0]);
        let cr = [
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        ];
        if v[0] == v[1] || v[1] == v[2] || v[0] == v[2] || dot(cr, cr).sqrt() < 1e-12 * diag * diag {
            d.degenerate += 1;
            continue;
        }
        let mut key = v;
        key.sort_unstable();
        if !seen.insert(key) {
            d.duplicates += 1;
            continue;
        }
        for (a, b) in [(v[0], v[1]), (v[1], v[2]), (v[2], v[0])] {
            let e = edges.entry((a.min(b), a.max(b))).or_insert((0, 0));
            e.0 += 1;
            e.1 += if a < b { 1 } else { -1 };
        }
        d.volume += dot(t[0], cr_of(t[1], t[2])) / 6.0;
        out.push(t);
    }
    for &(n, s) in edges.values() {
        match n {
            1 => d.open_edges += 1,
            2 if s != 0 => d.flipped_edges += 1,
            2 => {}
            _ => d.nonmanifold_edges += 1,
        }
    }
    d.volume = d.volume.abs();
    (out, d)
}

fn cr_of(b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
    [b[1] * c[2] - b[2] * c[1], b[2] * c[0] - b[0] * c[2], b[0] * c[1] - b[1] * c[0]]
}

/// Rangement des triangles pour des rayons parallèles à l'axe `axis` : cases
/// unité dans le plan des deux autres axes (u, v).
#[derive(Debug)]
struct RayBins {
    axis: usize,
    u: usize,
    v: usize,
    lo: [f64; 2],
    n: [usize; 2],
    bins: Vec<Vec<u32>>,
}

impl RayBins {
    fn new(tris: &[Tri], lo: [f64; 3], hi: [f64; 3], axis: usize) -> RayBins {
        let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
        let n = [
            (hi[u] - lo[u]).floor() as usize + 1,
            (hi[v] - lo[v]).floor() as usize + 1,
        ];
        let mut bins = vec![Vec::new(); n[0] * n[1]];
        for (k, t) in tris.iter().enumerate() {
            let (tl, th) = bbox(std::slice::from_ref(t));
            let i0 = ((tl[u] - lo[u]).floor() as usize).min(n[0] - 1);
            let i1 = ((th[u] - lo[u] + JU).floor() as usize).min(n[0] - 1);
            let j0 = ((tl[v] - lo[v]).floor() as usize).min(n[1] - 1);
            let j1 = ((th[v] - lo[v] + JV).floor() as usize).min(n[1] - 1);
            for j in j0..=j1 {
                for i in i0..=i1 {
                    bins[j * n[0] + i].push(k as u32);
                }
            }
        }
        RayBins {
            axis,
            u,
            v,
            lo: [lo[u], lo[v]],
            n,
            bins,
        }
    }

    fn bin(&self, pu: f64, pv: f64) -> Option<&[u32]> {
        let (fu, fv) = (pu - self.lo[0], pv - self.lo[1]);
        if fu < 0.0 || fv < 0.0 {
            return None;
        }
        let (i, j) = (fu.floor() as usize, fv.floor() as usize);
        (i < self.n[0] && j < self.n[1]).then(|| self.bins[j * self.n[0] + i].as_slice())
    }

    /// Coordonnée selon l'axe où la droite (pu, pv) traverse le triangle, si elle le traverse.
    fn cross(&self, t: &Tri, pu: f64, pv: f64) -> Option<f64> {
        let [a, b, c] = *t;
        let (u, v, w) = (self.u, self.v, self.axis);
        let det = (b[u] - a[u]) * (c[v] - a[v]) - (c[u] - a[u]) * (b[v] - a[v]);
        if det.abs() < 1e-14 {
            return None;
        }
        let s = ((pu - a[u]) * (c[v] - a[v]) - (c[u] - a[u]) * (pv - a[v])) / det;
        let r = ((b[u] - a[u]) * (pv - a[v]) - (pu - a[u]) * (b[v] - a[v])) / det;
        if s < 0.0 || r < 0.0 || s + r > 1.0 {
            return None;
        }
        Some(a[w] + s * (b[w] - a[w]) + r * (c[w] - a[w]))
    }

    /// Nombre de traversées du rayon partant de p vers les coordonnées croissantes.
    fn crossings(&self, tris: &[Tri], p: [f64; 3]) -> usize {
        let (pu, pv) = (p[self.u] + JU, p[self.v] + JV);
        let Some(bin) = self.bin(pu, pv) else { return 0 };
        bin.iter()
            .filter(|&&k| self.cross(&tris[k as usize], pu, pv).is_some_and(|c| c > p[self.axis]))
            .count()
    }
}

/// Maillage placé dans la veine (coordonnées réseau), avec ses grilles de rangement.
#[derive(Debug)]
pub struct MeshBody {
    pub tris: Vec<Tri>,
    pub lo: [f64; 3],
    pub hi: [f64; 3],
    /// Rayons selon x ; avec un maillage ouvert, aussi selon y et z (vote).
    rays: Vec<RayBins>,
    /// Rangement 3D, cases de `cell` mailles, pour la distance.
    grid: Vec<Vec<u32>>,
    gdim: [usize; 3],
    cell: f64,
}

impl MeshBody {
    /// `robust` : maillage ouvert ou mal formé, appartenance votée sur trois
    /// rayons (x, y, z) au lieu de la parité d'un seul.
    pub fn new(tris: Vec<Tri>, robust: bool) -> MeshBody {
        let (lo, hi) = bbox(&tris);
        let cell = 1.0;
        let gdim = [0, 1, 2].map(|a| ((hi[a] - lo[a]) / cell).floor() as usize + 1);
        let mut grid = vec![Vec::new(); gdim[0] * gdim[1] * gdim[2]];
        for (k, t) in tris.iter().enumerate() {
            let (tl, th) = bbox(std::slice::from_ref(t));
            let g0 = [0, 1, 2].map(|a| (((tl[a] - lo[a]) / cell).floor() as usize).min(gdim[a] - 1));
            let g1 = [0, 1, 2].map(|a| (((th[a] - lo[a]) / cell).floor() as usize).min(gdim[a] - 1));
            for gz in g0[2]..=g1[2] {
                for gy in g0[1]..=g1[1] {
                    for gx in g0[0]..=g1[0] {
                        grid[gx + gdim[0] * (gy + gdim[1] * gz)].push(k as u32);
                    }
                }
            }
        }
        let axes: &[usize] = if robust { &[0, 1, 2] } else { &[0] };
        let rays = axes.iter().map(|&a| RayBins::new(&tris, lo, hi, a)).collect();
        MeshBody {
            tris,
            lo,
            hi,
            rays,
            grid,
            gdim,
            cell,
        }
    }

    /// Le même maillage transformé point par point (grilles de rangement refaites).
    pub fn transformed(&self, t: impl Fn([f64; 3]) -> [f64; 3]) -> MeshBody {
        MeshBody::new(self.tris.iter().map(|tr| tr.map(&t)).collect(), self.rays.len() == 3)
    }

    /// Point dans le volume : parité des traversées d'un rayon selon x, ou
    /// vote à la majorité sur trois rayons pour un maillage ouvert.
    pub fn inside(&self, x: f64, y: f64, z: f64) -> bool {
        let p = [x, y, z];
        if (0..3).any(|a| p[a] < self.lo[a] || p[a] > self.hi[a]) {
            return false;
        }
        let votes = self.rays.iter().filter(|r| r.crossings(&self.tris, p) % 2 == 1).count();
        2 * votes > self.rays.len()
    }

    /// Première traversée de la surface par le segment p → p + d : fraction
    /// t ∈ ]0, 1] (Möller et Trumbore), ou None.
    pub fn segment_hit(&self, p: [f64; 3], d: [f64; 3]) -> Option<f64> {
        self.segment_hit_normal(p, d).map(|h| h.0)
    }

    /// Idem, avec la normale unitaire du triangle touché, orientée vers p.
    pub fn segment_hit_normal(&self, p: [f64; 3], d: [f64; 3]) -> Option<(f64, [f64; 3])> {
        let q = [p[0] + d[0], p[1] + d[1], p[2] + d[2]];
        let lo = [0, 1, 2].map(|a| p[a].min(q[a]));
        let hi = [0, 1, 2].map(|a| p[a].max(q[a]));
        if (0..3).any(|a| hi[a] < self.lo[a] || lo[a] > self.hi[a]) {
            return None;
        }
        let r = |a: usize, v: f64| (((v - self.lo[a]) / self.cell).floor().max(0.0) as usize).min(self.gdim[a] - 1);
        let mut best: Option<(f64, [f64; 3])> = None;
        for gz in r(2, lo[2])..=r(2, hi[2]) {
            for gy in r(1, lo[1])..=r(1, hi[1]) {
                for gx in r(0, lo[0])..=r(0, hi[0]) {
                    for &k in &self.grid[gx + self.gdim[0] * (gy + self.gdim[1] * gz)] {
                        let [a, b, c] = self.tris[k as usize];
                        let e1 = sub(b, a);
                        let e2 = sub(c, a);
                        let h = cr_of(d, e2);
                        let det = dot(e1, h);
                        if det.abs() < 1e-14 {
                            continue;
                        }
                        let inv = 1.0 / det;
                        let s = sub(p, a);
                        let u = inv * dot(s, h);
                        if !(0.0..=1.0).contains(&u) {
                            continue;
                        }
                        let qv = cr_of(s, e1);
                        let v = inv * dot(d, qv);
                        if v < 0.0 || u + v > 1.0 {
                            continue;
                        }
                        let t = inv * dot(e2, qv);
                        if t > 1e-9 && t <= 1.0 && best.is_none_or(|b| t < b.0) {
                            let nv = cr_of(e1, e2);
                            let ln = dot(nv, nv).sqrt().max(1e-300);
                            let sg = if dot(nv, d) > 0.0 { -1.0 } else { 1.0 };
                            best = Some((t, nv.map(|v| sg * v / ln)));
                        }
                    }
                }
            }
        }
        best
    }

    /// Distance à la surface, plafonnée à `cap` mailles.
    pub fn distance(&self, p: [f64; 3], cap: f64) -> f64 {
        let mut best = cap * cap;
        let g0 = [0, 1, 2].map(|a| ((p[a] - cap - self.lo[a]) / self.cell).floor());
        let g1 = [0, 1, 2].map(|a| ((p[a] + cap - self.lo[a]) / self.cell).floor());
        let r = |a: usize, v: f64| v.clamp(0.0, (self.gdim[a] - 1) as f64) as usize;
        if (0..3).any(|a| g1[a] < 0.0 || g0[a] > (self.gdim[a] - 1) as f64) {
            return cap;
        }
        for gz in r(2, g0[2])..=r(2, g1[2]) {
            for gy in r(1, g0[1])..=r(1, g1[1]) {
                for gx in r(0, g0[0])..=r(0, g1[0]) {
                    for &k in &self.grid[gx + self.gdim[0] * (gy + self.gdim[1] * gz)] {
                        best = best.min(tri_dist2(p, &self.tris[k as usize]));
                    }
                }
            }
        }
        best.sqrt()
    }

    /// Maître-couple : aire (mailles²) de la projection sur le plan (y, z),
    /// échantillonnée au quart de maille.
    pub fn frontal_area(&self) -> f64 {
        let r = &self.rays[0];
        let s = 0.25;
        let (my, mz) = (
            ((self.hi[1] - self.lo[1]) / s).ceil() as usize,
            ((self.hi[2] - self.lo[2]) / s).ceil() as usize,
        );
        let mut count = 0usize;
        for k in 0..mz {
            for j in 0..my {
                let y = self.lo[1] + (j as f64 + 0.5) * s + JU;
                let z = self.lo[2] + (k as f64 + 0.5) * s + JV;
                if let Some(bin) = r.bin(y, z) {
                    if bin.iter().any(|&t| r.cross(&self.tris[t as usize], y, z).is_some()) {
                        count += 1;
                    }
                }
            }
        }
        count as f64 * s * s
    }
}

/// Icosphère de rayon r centrée en c.
pub fn icosphere(r: f64, c: [f64; 3], sub: usize) -> Vec<Tri> {
    let t = (1.0 + 5f64.sqrt()) / 2.0;
    let mut v: Vec<[f64; 3]> = vec![
        [-1.0, t, 0.0],
        [1.0, t, 0.0],
        [-1.0, -t, 0.0],
        [1.0, -t, 0.0],
        [0.0, -1.0, t],
        [0.0, 1.0, t],
        [0.0, -1.0, -t],
        [0.0, 1.0, -t],
        [t, 0.0, -1.0],
        [t, 0.0, 1.0],
        [-t, 0.0, -1.0],
        [-t, 0.0, 1.0],
    ];
    let mut f: Vec<[usize; 3]> = vec![
        [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9], [5, 11, 4], [11, 10, 2],
        [10, 7, 6], [7, 1, 8], [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9], [4, 9, 5],
        [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],
    ];
    for _ in 0..sub {
        let mut nf = Vec::new();
        let mut mid = std::collections::HashMap::new();
        let mut m = |a: usize, b: usize, v: &mut Vec<[f64; 3]>| -> usize {
            *mid.entry((a.min(b), a.max(b))).or_insert_with(|| {
                v.push([(v[a][0] + v[b][0]) / 2.0, (v[a][1] + v[b][1]) / 2.0, (v[a][2] + v[b][2]) / 2.0]);
                v.len() - 1
            })
        };
        for [a, b, c] in f {
            let (ab, bc, ca) = (m(a, b, &mut v), m(b, c, &mut v), m(c, a, &mut v));
            nf.extend([[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]);
        }
        f = nf;
    }
    let p = |i: usize| {
        let n = dot(v[i], v[i]).sqrt();
        [c[0] + r * v[i][0] / n, c[1] + r * v[i][1] / n, c[2] + r * v[i][2] / n]
    };
    f.iter().map(|&[a, b, cc]| [p(a), p(b), p(cc)]).collect()
}

/// Écrit des triangles en STL binaire.
pub fn write_stl(path: &Path, tris: &[Tri]) -> std::io::Result<()> {
    let mut b = vec![0u8; 80];
    b.extend((tris.len() as u32).to_le_bytes());
    for t in tris {
        b.extend([0u8; 12]);
        for c in t.iter().flatten() {
            b.extend((*c as f32).to_le_bytes());
        }
        b.extend([0u8; 2]);
    }
    std::fs::write(path, b)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sphere_maillee() {
        let m = MeshBody::new(icosphere(10.0, [20.0, 20.0, 20.0], 4), false);
        assert!(m.inside(20.0, 20.0, 20.0));
        assert!(m.inside(29.0, 20.0, 20.0));
        assert!(!m.inside(31.0, 20.0, 20.0));
        assert!(!m.inside(28.0, 28.0, 20.0));
        let d = m.distance([20.0, 20.0, 33.0], 5.0);
        assert!((d - 3.0).abs() < 0.05, "distance {d}");
        let a = m.frontal_area();
        let exact = std::f64::consts::PI * 100.0;
        assert!((a / exact - 1.0).abs() < 0.02, "maître-couple {a} contre {exact}");
        // Volume par comptage de mailles.
        let mut n = 0;
        for z in 5..36 {
            for y in 5..36 {
                for x in 5..36 {
                    n += m.inside(x as f64, y as f64, z as f64) as usize;
                }
            }
        }
        let v = 4.0 / 3.0 * std::f64::consts::PI * 1000.0;
        assert!((n as f64 / v - 1.0).abs() < 0.03, "volume {n} contre {v}");
    }

    #[test]
    fn maillage_troue() {
        let mut tris = icosphere(10.0, [20.0, 20.0, 20.0], 4);
        let full = tris.len();
        // Défauts : un triangle en double, un dégénéré, puis des trous.
        tris.push(tris[7]);
        tris.push([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]]);
        let (ok, d) = clean(tris.clone());
        assert_eq!((ok.len(), d.duplicates, d.degenerate), (full, 1, 1));
        assert!(d.watertight() && d.flipped_edges == 0, "{d:?}");
        let v = 4.0 / 3.0 * std::f64::consts::PI * 1000.0;
        assert!((d.volume / v - 1.0).abs() < 0.01, "volume {}", d.volume);
        // Une triangle sur 40 retiré : maillage ouvert.
        let holed: Vec<Tri> = ok.iter().enumerate().filter(|(k, _)| k % 40 != 0).map(|(_, t)| *t).collect();
        let (holed, d) = clean(holed);
        assert!(d.open_edges > 0);
        let count = |m: &MeshBody| {
            let mut n = 0;
            for z in 5..36 {
                for y in 5..36 {
                    for x in 5..36 {
                        n += m.inside(x as f64, y as f64, z as f64) as usize;
                    }
                }
            }
            n as f64
        };
        let single = count(&MeshBody::new(holed.clone(), false));
        let voted = count(&MeshBody::new(holed, true));
        println!("volume exact {v:.0}, un rayon {single}, vote {voted}");
        assert!((voted / v - 1.0).abs() < 0.03, "vote {voted} contre {v}");
        assert!((voted / v - 1.0).abs() < (single / v - 1.0).abs(), "un rayon {single}, vote {voted}");
    }

    #[test]
    fn lecture_ascii_et_binaire() {
        let tris = icosphere(1.0, [0.0; 3], 1);
        let mut ascii = String::from("solid s\n");
        for t in &tris {
            ascii += "facet normal 0 0 0\nouter loop\n";
            for p in t {
                ascii += &format!("vertex {} {} {}\n", p[0], p[1], p[2]);
            }
            ascii += "endloop\nendfacet\n";
        }
        ascii += "endsolid s\n";
        assert_eq!(parse_stl(ascii.as_bytes()).unwrap().len(), tris.len());
        let mut bin = vec![0u8; 80];
        bin.extend((tris.len() as u32).to_le_bytes());
        for t in &tris {
            bin.extend([0u8; 12]);
            for p in t {
                for c in p {
                    bin.extend((*c as f32).to_le_bytes());
                }
            }
            bin.extend([0u8; 2]);
        }
        let b = parse_stl(&bin).unwrap();
        assert_eq!(b.len(), tris.len());
        assert!((b[3][1][2] - tris[3][1][2]).abs() < 1e-6);
    }
}
