//! Maquettes 3D décrites par un maillage STL (binaire ou ASCII) : lecture,
//! mise en place dans la veine et requêtes géométriques accélérées par des
//! grilles de rangement (appartenance par lancer de rayon selon x, distance
//! à la surface, maître-couple).

use std::path::Path;

pub type Tri = [[f64; 3]; 3];

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
        if data.len() == 84 + 50 * n {
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

/// Décalage du rayon d'appartenance : évite de passer exactement par une
/// arête ou un sommet (coordonnées entières des mailles).
const JY: f64 = 1.37e-5;
const JZ: f64 = 2.71e-5;

/// Maillage placé dans la veine (coordonnées réseau), avec ses grilles de rangement.
#[derive(Debug)]
pub struct MeshBody {
    pub tris: Vec<Tri>,
    pub lo: [f64; 3],
    pub hi: [f64; 3],
    /// Rangement en (y, z), mailles unité : triangles dont la projection touche la case.
    yz: Vec<Vec<u32>>,
    ny: usize,
    nz: usize,
    /// Rangement 3D, cases de `cell` mailles, pour la distance.
    grid: Vec<Vec<u32>>,
    gdim: [usize; 3],
    cell: f64,
}

impl MeshBody {
    pub fn new(tris: Vec<Tri>) -> MeshBody {
        let (lo, hi) = bbox(&tris);
        let ny = (hi[1] - lo[1]).floor() as usize + 1;
        let nz = (hi[2] - lo[2]).floor() as usize + 1;
        let mut yz = vec![Vec::new(); ny * nz];
        let cell = 2.0;
        let gdim = [0, 1, 2].map(|a| ((hi[a] - lo[a]) / cell).floor() as usize + 1);
        let mut grid = vec![Vec::new(); gdim[0] * gdim[1] * gdim[2]];
        for (k, t) in tris.iter().enumerate() {
            let (tl, th) = bbox(std::slice::from_ref(t));
            let j0 = ((tl[1] - lo[1]).floor() as usize).min(ny - 1);
            let j1 = ((th[1] - lo[1] + JY).floor() as usize).min(ny - 1);
            let k0 = ((tl[2] - lo[2]).floor() as usize).min(nz - 1);
            let k1 = ((th[2] - lo[2] + JZ).floor() as usize).min(nz - 1);
            for kz in k0..=k1 {
                for jy in j0..=j1 {
                    yz[kz * ny + jy].push(k as u32);
                }
            }
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
        MeshBody {
            tris,
            lo,
            hi,
            yz,
            ny,
            nz,
            grid,
            gdim,
            cell,
        }
    }

    /// Abscisse où la droite (y, z) parallèle à x traverse le triangle k, si elle le traverse.
    fn cross_x(&self, k: u32, y: f64, z: f64) -> Option<f64> {
        let [a, b, c] = self.tris[k as usize];
        // Coordonnées barycentriques dans la projection (y, z).
        let det = (b[1] - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (b[2] - a[2]);
        if det.abs() < 1e-14 {
            return None;
        }
        let u = ((y - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (z - a[2])) / det;
        let v = ((b[1] - a[1]) * (z - a[2]) - (y - a[1]) * (b[2] - a[2])) / det;
        if u < 0.0 || v < 0.0 || u + v > 1.0 {
            return None;
        }
        Some(a[0] + u * (b[0] - a[0]) + v * (c[0] - a[0]))
    }

    fn yz_bin(&self, y: f64, z: f64) -> Option<&[u32]> {
        let (fy, fz) = (y - self.lo[1], z - self.lo[2]);
        if fy < 0.0 || fz < 0.0 {
            return None;
        }
        let (j, k) = (fy.floor() as usize, fz.floor() as usize);
        (j < self.ny && k < self.nz).then(|| self.yz[k * self.ny + j].as_slice())
    }

    /// Point dans le volume fermé : parité des traversées d'un rayon vers +x.
    pub fn inside(&self, x: f64, y: f64, z: f64) -> bool {
        if x < self.lo[0] || x > self.hi[0] {
            return false;
        }
        let (y, z) = (y + JY, z + JZ);
        let Some(bin) = self.yz_bin(y, z) else { return false };
        let mut n = 0;
        for &k in bin {
            if let Some(xc) = self.cross_x(k, y, z) {
                if xc > x {
                    n += 1;
                }
            }
        }
        n % 2 == 1
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
        let s = 0.25;
        let (my, mz) = (
            ((self.hi[1] - self.lo[1]) / s).ceil() as usize,
            ((self.hi[2] - self.lo[2]) / s).ceil() as usize,
        );
        let mut count = 0usize;
        for k in 0..mz {
            for j in 0..my {
                let y = self.lo[1] + (j as f64 + 0.5) * s + JY;
                let z = self.lo[2] + (k as f64 + 0.5) * s + JZ;
                if let Some(bin) = self.yz_bin(y, z) {
                    if bin.iter().any(|&t| self.cross_x(t, y, z).is_some()) {
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
        let m = MeshBody::new(icosphere(10.0, [20.0, 20.0, 20.0], 4));
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
