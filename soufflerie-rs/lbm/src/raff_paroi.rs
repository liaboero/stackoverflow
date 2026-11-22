//! Raffinement de paroi par zones de gradient de pression adverse
//! (`aile_infinie --raff-paroi n`).
//!
//! La loi de paroi d'équilibre (Spalding) suppose une contrainte constante
//! entre la paroi et le point de raccord : τ(h) = τ_p. Avec un gradient de
//! pression, la couche interne porte τ(y) = τ_p + y dp/ds (Townsend, Mellor) :
//! la loi n'est plus valide au raccord dès que h dp/ds dépasse τ_p, soit
//!
//!   β_h = h · (dp/ds) / τ_p > 1,
//!
//! ce qui arrive en gradient adverse quand h n'est plus dans la couche
//! interne (README, « Pourquoi la soufflerie ne décolle pas » : h ≥ δ,
//! P h / u_τ² ≈ 10). Ce critère, sans constante de réglage, choisit les
//! zones où la couche limite doit être résolue au-dessus de h : un premier
//! calcul (la grille de base) donne, par triangle de la surface, la pression
//! pariétale et la contrainte de la loi de paroi moyennées dans le temps, et
//! β à la hauteur de raccord de cette grille. Chaque niveau ajouté divise h
//! par deux : le niveau j (j = 1…n, 2ʲ fois plus fin que la grille de base)
//! couvre les triangles où β à la hauteur du niveau j − 1 dépasse 1, soit
//! β₀ > 2ʲ⁻¹. On raffine ainsi tant que la loi d'équilibre est fausse au
//! raccord, dans la limite de n niveaux (coût).
//!
//! Gradient le long de la surface : dp/ds = ∇p · ŝ, ŝ direction de la
//! contrainte pariétale (vitesse relevée à h), par moindres carrés à une
//! dimension sur les triangles voisins (centres à moins de 2h, de part et
//! d'autre de h le long de ŝ, faces de même côté d'une arête vive :
//! normales à moins de `geom3d::FEATURE_DEG`).
//!
//! Les niveaux sont des couches creuses (`LayerSpec::zone`) imbriquées sous
//! les grilles boîtes de base : épaisseur selon la règle de `raff_auto`
//! (max(16 mailles, 0,125 L) au niveau le plus fin, puis `nested_bands`),
//! sillage en fraction de l'épaisseur (0 conseillé : sur l'aile infinie, le
//! sillage balayé des couches creuses ajoute au bilan une traînée qui croît
//! avec la marge, README « Raffinement de paroi par zones »).
//!
//! Limite connue : là où la zone s'arrête sur la paroi, le raccord coupe la
//! surface (comme les grilles des bords) et ajoute de la traînée.

use crate::gpu3d::{LayerSpec, RefinePlan, SurfaceForces, ZoneRegion};
use crate::stl::Tri;

/// Moyennes temporelles par triangle de la pression pariétale, de la
/// contrainte de la loi de paroi et de la hauteur de raccord.
#[derive(Clone, Debug)]
pub struct ZoneStats {
    p: Vec<f64>,
    tau: Vec<[f64; 3]>,
    h: Vec<f64>,
    n: Vec<u32>,
}

impl ZoneStats {
    pub fn new(ntri: usize) -> ZoneStats {
        ZoneStats { p: vec![0.0; ntri], tau: vec![[0.0; 3]; ntri], h: vec![0.0; ntri], n: vec![0; ntri] }
    }

    /// Ajoute un relevé de `Gpu3D::surface_forces`.
    pub fn add(&mut self, s: &SurfaceForces) {
        for (k, &t) in s.tri.iter().enumerate() {
            if t >= self.n.len() {
                continue;
            }
            self.p[t] += s.walls[k].1;
            for a in 0..3 {
                self.tau[t][a] += s.shear[k][a];
            }
            self.h[t] += s.wm_h[k];
            self.n[t] += 1;
        }
    }

    /// β_h = h (dp/ds)⁺ / |τ_p| par triangle (NaN sans relevé), au raccord
    /// de la grille qui a relevé le triangle. `tris` : maillage de la
    /// maquette (mailles de la veine, ordre des relevés).
    pub fn beta(&self, tris: &[Tri]) -> Vec<f64> {
        use rayon::prelude::*;
        let nt = tris.len().min(self.n.len());
        let center = |t: &Tri| [0, 1, 2].map(|a| (t[0][a] + t[1][a] + t[2][a]) / 3.0);
        let normal = |t: &Tri| {
            let e1 = [t[1][0] - t[0][0], t[1][1] - t[0][1], t[1][2] - t[0][2]];
            let e2 = [t[2][0] - t[0][0], t[2][1] - t[0][1], t[2][2] - t[0][2]];
            let n = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]];
            let m = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt().max(1e-300);
            n.map(|v| v / m)
        };
        let c: Vec<[f64; 3]> = tris[..nt].iter().map(center).collect();
        let nn: Vec<[f64; 3]> = tris[..nt].iter().map(normal).collect();
        let ok: Vec<bool> = (0..nt).map(|i| self.n[i] > 0).collect();
        let pm: Vec<f64> = (0..nt).map(|i| if ok[i] { self.p[i] / self.n[i] as f64 } else { 0.0 }).collect();
        let hm: Vec<f64> = (0..nt).map(|i| if ok[i] { self.h[i] / self.n[i] as f64 } else { 0.0 }).collect();
        let tm: Vec<[f64; 3]> = (0..nt).map(|i| if ok[i] { self.tau[i].map(|v| v / self.n[i] as f64) } else { [0.0; 3] }).collect();
        // Rangement des centres par cases de côté 2 h max.
        let r_max = 2.0 * hm.iter().cloned().fold(0.0, f64::max).max(1e-9);
        let key = |p: [f64; 3]| p.map(|v| (v / r_max).floor() as i64);
        let mut grid: std::collections::HashMap<[i64; 3], Vec<usize>> = std::collections::HashMap::new();
        for i in 0..nt {
            if ok[i] {
                grid.entry(key(c[i])).or_default().push(i);
            }
        }
        let cos_f = crate::geom3d::FEATURE_DEG.to_radians().cos();
        (0..nt)
            .into_par_iter()
            .map(|i| {
                if !ok[i] {
                    return f64::NAN;
                }
                let t = tm[i];
                let tn = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt();
                if tn <= 0.0 {
                    return f64::NAN;
                }
                let s = t.map(|v| v / tn);
                let r = 2.0 * hm[i];
                let k = key(c[i]);
                let (mut num, mut den) = (0.0, 0.0);
                for dz in -1..=1 {
                    for dy in -1..=1 {
                        for dx in -1..=1 {
                            let Some(v) = grid.get(&[k[0] + dx, k[1] + dy, k[2] + dz]) else { continue };
                            for &j in v {
                                if j == i {
                                    continue;
                                }
                                let d = [0, 1, 2].map(|a| c[j][a] - c[i][a]);
                                if d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > r * r {
                                    continue;
                                }
                                if nn[i][0] * nn[j][0] + nn[i][1] * nn[j][1] + nn[i][2] * nn[j][2] < cos_f {
                                    continue;
                                }
                                let ds = d[0] * s[0] + d[1] * s[1] + d[2] * s[2];
                                num += (pm[j] - pm[i]) * ds;
                                den += ds * ds;
                            }
                        }
                    }
                }
                if den <= 0.0 {
                    return f64::NAN;
                }
                hm[i] * (num / den).max(0.0) / tn
            })
            .collect()
    }
}

/// Masques de triangles des niveaux 1…n : β₀ > 2ʲ⁻¹ (et `keep`).
pub fn zones(beta: &[f64], n: usize, keep: impl Fn(usize) -> bool) -> Vec<Vec<bool>> {
    (1..=n).map(|j| beta.iter().enumerate().map(|(i, &b)| b > (1u64 << (j - 1)) as f64 && keep(i)).collect()).collect()
}

/// Épaisseurs (mailles de chaque niveau) des couches par zones sous `depth`
/// grilles boîtes d'une veine de `l` mailles par longueur de référence : au
/// niveau le plus fin, max(16 mailles, 0,125 L) (`raff_auto`), ou `finest`
/// s'il est donné ; puis emboîtement de `nested_bands`.
pub fn bands(l: f64, depth: usize, n: usize, finest: Option<f64>) -> Vec<f64> {
    let s = (1u64 << (depth + n)) as f64;
    let d = finest.unwrap_or_else(|| (crate::raff_auto::OFFSET_CELLS).max((crate::raff_auto::OFFSET_CHORD * l * s).ceil()));
    let mut b = vec![0.0; n];
    b[n - 1] = d;
    crate::gpu3d::nested_bands(&b)
}

/// Plan : `depth` boîtes de base (`refine3d::EXTENT3`), puis une couche par
/// zone (masques `zones`, épaisseurs `bands`), sillage de `wake` fois
/// l'épaisseur.
pub fn plan(depth: usize, zones: &[Vec<bool>], bands: &[f64], wake: f64) -> RefinePlan {
    let depth = depth.min(crate::refine3d::EXTENT3.len());
    RefinePlan {
        boxes: crate::refine3d::EXTENT3[..depth].to_vec(),
        layers: zones
            .iter()
            .zip(bands)
            .map(|(z, &b)| LayerSpec { band: b, wake: (wake * b).ceil(), edges: None, zone: Some(std::sync::Arc::new(ZoneRegion { tris: z.clone() })) })
            .collect(),
        le: false,
        te: false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn epaisseurs_emboitees() {
        // L = 16, trois boîtes (×8), deux niveaux : ×32 de 0,125 L = 64
        // mailles, ×16 d'au moins 64/2 + 7 = 39.
        let b = bands(16.0, 3, 2, None);
        assert_eq!(b, vec![39.0, 64.0]);
        assert_eq!(bands(16.0, 3, 1, Some(24.0)), vec![24.0]);
    }

    #[test]
    fn beta_rampe_de_pression() {
        // Plaque plane z = 0 en bandes : p croît linéairement le long de x
        // (gradient adverse dans le sens de τ), τ_p = 0,01 · x̂, h = 1 :
        // β = h dp/dx / τ = 0,5 pour dp/dx = 0,005.
        let mut tris = Vec::new();
        for i in 0..40 {
            for j in 0..4 {
                let (x, y) = (0.25 * i as f64, j as f64);
                tris.push([[x, y, 0.0], [x + 0.25, y, 0.0], [x + 0.25, y + 1.0, 0.0]]);
                tris.push([[x, y, 0.0], [x + 0.25, y + 1.0, 0.0], [x, y + 1.0, 0.0]]);
            }
        }
        let mut s = SurfaceForces::default();
        for (k, t) in tris.iter().enumerate() {
            let cx = (t[0][0] + t[1][0] + t[2][0]) / 3.0;
            s.tri.push(k);
            s.walls.push(([0.0, 0.0, 1.0], 0.005 * cx));
            s.shear.push([0.01, 0.0, 0.0]);
            s.wm_h.push(1.0);
        }
        let mut z = ZoneStats::new(tris.len());
        z.add(&s);
        let b = z.beta(&tris);
        for v in &b {
            assert!((v - 0.5).abs() < 1e-9, "{v}");
        }
        // Gradient favorable : β = 0.
        let mut s2 = s.clone();
        for w in &mut s2.walls {
            w.1 = -w.1;
        }
        let mut z2 = ZoneStats::new(tris.len());
        z2.add(&s2);
        assert!(z2.beta(&tris).iter().all(|&v| v == 0.0));
        let zs = zones(&b.iter().map(|v| v * 5.0).collect::<Vec<_>>(), 3, |_| true);
        assert!(zs[0].iter().all(|&k| k) && zs[1].iter().all(|&k| k) && zs[2].iter().all(|&k| !k));
    }
}
