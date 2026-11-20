//! Raffinement automatique d'après la géométrie (`run3d --raff-auto N`), à
//! la manière des régions décalées de la paroi des codes Boltzmann
//! commerciaux : un seul paramètre, N mailles par longueur de référence à
//! la paroi, et des règles de décalage tirées des essais déjà faits.
//!
//! Niveau k (2ᵏ fois plus fin que la veine, k = 1…n, 2ⁿ L = N) : points à
//! moins de d_k de la surface, plus un sillage balayé vers l'aval sur
//! W_k = s · d_k (s = 1 par défaut, `--raff-auto-sillage`). Règles :
//! - couche la plus fine : d_n = max(16 mailles ×2ⁿ, 0,125 L_ref), épaisseur
//!   sous laquelle le raccord tombe dans la couche limite et la zone de la loi
//!   de paroi (README, « Couche pariétale fine » et `--bande-corde`) ;
//! - niveaux plus grossiers : d_k = d_{k+1} + m_k mailles ×2ᵏ, marge du
//!   raccord : pochoir cubique (±2 mailles parentes) et zone non recopiée
//!   (MARGIN), soit m = 2 MARGIN = 4 pour une boîte ; pour une couche, la
//!   règle de `nested_bands`, m = 2 MARGIN + 3 = 7 ;
//! - d_k croît donc avec la taille de maille : chaque raccord est à nombre
//!   de mailles à peu près constant de la paroi.
//!
//! Les grilles imbriquées sont des boîtes ou des couches creuses. Une boîte
//! est l'enveloppe exacte de sa région : maquette ± d_k, plus W_k vers
//! l'aval. Par défaut, les niveaux 1 à n − 1 sont des boîtes et seul le
//! niveau de la paroi est une couche creuse (`LayerSpec`, sillage compris) :
//! le bilan sur boîte, mesure de référence de la traînée, s'intègre sur la
//! grille entière la plus fine, qui doit rester à un niveau de la paroi
//! (sur la grille ×2 de TR-1208, il donne un Cx négatif). Avec
//! `--raff-auto-creux 1`, un niveau n'est une boîte que si elle ne fait pas
//! plus de deux fois le volume de sa région (`BOX_RATIO` : le stockage par
//! blocs de 8³ garde 1,6 fois les mailles actives, blocs entamés et
//! fantômes) ; le premier niveau reste une boîte.
//!
//! Bords (`--raff-auto-bords 1`, essai) : un niveau de plus, limité aux
//! arêtes vives (angle dièdre ≥ 60°, `geom3d::FEATURE_DEG` : bord de fuite,
//! saumon plat) et aux arêtes dont le rayon de courbure local est sous deux
//! mailles de la couche la plus fine (bord d'attaque), à moins de 16 mailles
//! de ce niveau. Cette région coupe la surface, comme les grilles des bords
//! (`--ba`, `--bf`) : voir le README avant de s'y fier.

use crate::bodies3d::Body3D;
use crate::geom3d::{FeatureEdge, FeatureSet, Indexed, FEATURE_DEG};
use crate::gpu3d::{EdgeRegion, LayerSpec, RefinePlan};

/// Épaisseur minimale de la couche la plus fine, en ses mailles.
pub const OFFSET_CELLS: f64 = 16.0;
/// Idem en longueur de référence.
pub const OFFSET_CHORD: f64 = 0.125;
/// Marge de raccord (mailles parentes) sous une boîte et sous une couche.
pub const MARGIN_BOX: f64 = (2 * crate::refine3d::MARGIN3) as f64;
pub const MARGIN_LAYER: f64 = (2 * crate::refine3d::MARGIN3 + 3) as f64;
/// Une boîte reste préférée à une couche creuse jusqu'à ce rapport des volumes.
pub const BOX_RATIO: f64 = 2.0;

#[derive(Clone, Debug)]
pub struct AutoOptions {
    /// Mailles par longueur de référence à la paroi (L · 2ⁿ).
    pub wall_cells: usize,
    /// Longueur du sillage en décalages : W_k = wake · d_k.
    pub wake: f64,
    /// Niveau de plus près des bords.
    pub edges: bool,
    /// Niveaux intermédiaires en couches creuses dès que la boîte dépasse
    /// `BOX_RATIO` fois la région (sinon boîtes jusqu'au niveau n − 1).
    pub sparse: bool,
}

impl Default for AutoOptions {
    fn default() -> Self {
        AutoOptions { wall_cells: 0, wake: 1.0, edges: false, sparse: false }
    }
}

/// Un niveau du plan (longueurs en longueurs de référence).
#[derive(Clone, Debug)]
pub struct AutoLevel {
    pub scale: u64,
    pub offset: f64,
    pub wake: f64,
    pub boxed: bool,
    /// Volumes de la boîte et de la région décalée (L_ref³).
    pub box_volume: f64,
    pub region_volume: f64,
}

#[derive(Clone, Debug)]
pub struct AutoPlan {
    pub levels: Vec<AutoLevel>,
    /// Niveau des bords : nombre d'arêtes retenues et rayon (L_ref).
    pub edges: Option<(usize, f64)>,
    pub plan: RefinePlan,
}

impl AutoPlan {
    /// Résumé lisible (une ligne par niveau).
    pub fn describe(&self) -> String {
        let mut s = String::new();
        for l in &self.levels {
            s += &format!(
                "  ×{} : {} · décalage {:.3} L · sillage {:.3} L · boîte {:.2} L³, région {:.2} L³ (rapport {:.1})\n",
                l.scale,
                if l.boxed { "boîte" } else { "couche" },
                l.offset,
                l.wake,
                l.box_volume,
                l.region_volume,
                l.box_volume / l.region_volume.max(1e-12)
            );
        }
        if let Some((n, r)) = self.edges {
            s += &format!("  bords : {n} arêtes, rayon {r:.3} L\n");
        }
        s
    }
}

/// Décalages d_k (L_ref) des niveaux 1…n (indice k − 1), pour des niveaux
/// boîtes jusqu'à `nb` compris, couches au-delà ; `extra` : décalage d'un
/// niveau n + 1 (bords) à loger dans le niveau n.
pub fn offsets(l: f64, n: usize, nb: usize, extra: Option<f64>) -> Vec<f64> {
    let cell = |k: usize| 1.0 / (l * (1u64 << k) as f64);
    let mut d = vec![0.0; n];
    d[n - 1] = (OFFSET_CELLS * cell(n)).max(OFFSET_CHORD);
    if let Some(e) = extra {
        let m = if n <= nb { MARGIN_BOX } else { MARGIN_LAYER };
        d[n - 1] = d[n - 1].max(e + m * cell(n));
    }
    for k in (1..n).rev() {
        let m = if k <= nb { MARGIN_BOX } else { MARGIN_LAYER };
        d[k - 1] = d[k] + m * cell(k);
    }
    d
}

/// Volume (L_ref³) des points hors de la maquette à moins de `d` (L_ref) de
/// sa surface ou dans le sillage balayé sur `w`, par échantillonnage à
/// d/3 (au plus une maille de la veine).
fn region_volume(body: &Body3D, d: f64, w: f64) -> f64 {
    use rayon::prelude::*;
    let Some(mesh) = body.mesh.as_ref() else {
        return f64::NAN;
    };
    let l = body.l as f64;
    let (dc, wc) = (d * l, w * l);
    let h = (dc / 3.0).min(1.0);
    let b = body.bounds;
    let lo = [b.x0 - dc, b.y0 - dc, b.z0 - dc];
    let hi = [b.x1 + dc + wc, b.y1 + dc, b.z1 + dc];
    let n = [0, 1, 2].map(|a| ((hi[a] - lo[a]) / h).ceil() as usize + 1);
    let count: usize = (0..n[1] * n[2])
        .into_par_iter()
        .map(|r| {
            let (y, z) = (lo[1] + (r % n[1]) as f64 * h, lo[2] + (r / n[1]) as f64 * h);
            let mut last = f64::NEG_INFINITY;
            let mut c = 0;
            for i in 0..n[0] {
                let x = lo[0] + i as f64 * h;
                if body.inside(x, y, z) {
                    continue;
                }
                if x <= b.x1 + dc && mesh.distance([x, y, z], dc + 1.0) <= dc {
                    last = x;
                    c += 1;
                } else if x - last <= wc {
                    c += 1;
                }
            }
            c
        })
        .sum();
    count as f64 * (h / l).powi(3)
}

/// Arêtes à raffiner (coordonnées de la veine) : arêtes vives, et arêtes dont
/// le rayon de courbure local (hauteur moyenne des deux faces sur l'angle
/// entre leurs normales) est sous `r_max` mailles de la veine.
pub fn edge_set(body: &Body3D, r_max: f64) -> FeatureSet {
    let Some(mesh) = body.mesh.as_ref() else {
        return FeatureSet::default();
    };
    let m = Indexed::from_tris(&mesh.tris);
    let sgn = if m.signed_volume() >= 0.0 { 1.0 } else { -1.0 };
    let mut keys: Vec<_> = m.edge_faces().into_iter().collect();
    keys.sort_unstable_by_key(|(k, _)| *k);
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let norm = |a: [f64; 3]| (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let mut out = Vec::new();
    for ((i, j), fs) in keys {
        let (a, b) = (m.verts[i as usize], m.verts[j as usize]);
        let le = norm(sub(b, a));
        if fs.len() != 2 || le <= 0.0 {
            continue;
        }
        let (n0, n1) = (m.face_normal(fs[0] as usize), m.face_normal(fs[1] as usize));
        let (m0, m1) = (norm(n0), norm(n1));
        if m0 <= 0.0 || m1 <= 0.0 {
            continue;
        }
        let ang = (dot(n0, n1) / (m0 * m1)).clamp(-1.0, 1.0).acos();
        // Hauteur de chaque face au-dessus de l'arête : 2 aire / longueur.
        let h = 0.5 * (m0 + m1) / le;
        let f1 = m.faces[fs[1] as usize];
        let o = f1.iter().find(|&&v| v != i && v != j).map(|&v| m.verts[v as usize]).unwrap_or(a);
        let convex = sgn * dot(n0, sub(o, a)) < 0.0;
        let sharp = ang.to_degrees() >= FEATURE_DEG;
        let curved = convex && ang > 1e-9 && h / ang < r_max;
        if sharp || curved {
            out.push(FeatureEdge { a, b, angle_deg: ang.to_degrees(), convex });
        }
    }
    FeatureSet::new(out)
}

/// Plan automatique pour la maquette `body` (veine de L = body.l mailles par
/// longueur de référence).
pub fn auto_plan(body: &Body3D, o: &AutoOptions) -> Result<AutoPlan, String> {
    let l = body.l as f64;
    let ratio = o.wall_cells as f64 / l;
    let n = ratio.log2().round();
    if !(n >= 1.0) || ((1u64 << n as u32) as f64 * l - o.wall_cells as f64).abs() > 0.5 {
        return Err(format!("--raff-auto {} : il faut L · 2ⁿ mailles par longueur (L = {}, n ≥ 1)", o.wall_cells, body.l));
    }
    let n = n as usize;
    if n > 4 {
        return Err("--raff-auto : au plus quatre niveaux (×16)".into());
    }
    let extra = o.edges.then(|| OFFSET_CELLS / (l * (1u64 << (n + 1)) as f64));
    let b = body.bounds;
    let ext = [(b.x1 - b.x0) / l, (b.y1 - b.y0) / l, (b.z1 - b.z0) / l];
    let box_vol = |d: f64, w: f64| (ext[0] + 2.0 * d + w) * (ext[1] + 2.0 * d) * (ext[2] + 2.0 * d);
    // Niveaux boîtes : préfixe 1..nb ; décidé sur les décalages des couches
    // (les plus grands), puis décalages recalculés.
    let d_layer = offsets(l, n, 0, extra);
    let mut nb = 1;
    let mesh = body.mesh.is_some();
    let need = if o.sparse { 1 } else { n.saturating_sub(1).max(1) };
    while nb < n {
        let d = d_layer[nb];
        if nb < need {
            nb += 1;
            continue;
        }
        let w = o.wake * d;
        if mesh && box_vol(d, w) > BOX_RATIO * region_volume(body, d, w) {
            break;
        }
        nb += 1;
    }
    if !mesh && (nb < n || o.edges) {
        return Err("--raff-auto : couches et bords pour une maquette maillée seulement".into());
    }
    let d = offsets(l, n, nb, extra);
    let mut levels = Vec::new();
    let mut plan = RefinePlan::default();
    for k in 1..=n {
        let (dk, wk) = (d[k - 1], o.wake * d[k - 1]);
        let s = (1u64 << k) as f64;
        levels.push(AutoLevel {
            scale: 1 << k,
            offset: dk,
            wake: wk,
            boxed: k <= nb,
            box_volume: box_vol(dk, wk),
            region_volume: if mesh { region_volume(body, dk, wk) } else { box_vol(dk, wk) },
        });
        if k <= nb {
            plan.boxes.push((dk, dk + wk, dk));
        } else {
            plan.layers.push(LayerSpec { band: (dk * l * s).ceil(), wake: (wk * l * s).ceil(), edges: None });
        }
    }
    let mut edges = None;
    if let Some(e) = extra {
        let s = (1u64 << (n + 1)) as f64;
        // Rayon de courbure non résolu : moins de deux mailles du niveau n.
        let set = edge_set(body, 2.0 / (1u64 << n) as f64);
        if set.edges.is_empty() {
            return Err("--raff-auto-bords : aucune arête vive ni courbure sous deux mailles".into());
        }
        edges = Some((set.edges.len(), e));
        let band = (e * l * s).ceil();
        plan.layers.push(LayerSpec {
            band,
            wake: 0.0,
            edges: Some(std::sync::Arc::new(EdgeRegion { set, radius: band })),
        });
    }
    Ok(AutoPlan { levels, edges, plan })
}

/// `--raff-auto N` dans une liste d'arguments (N, ou None), avec
/// `--raff-auto-sillage s`, `--raff-auto-bords 1` et `--raff-auto-creux 1` :
/// options des exemples.
pub fn options_from_args(args: &[String]) -> Option<AutoOptions> {
    let val = |k: &str| args.iter().position(|a| a == k).and_then(|i| args.get(i + 1)).and_then(|v| v.parse::<f64>().ok());
    let n = val("--raff-auto")?;
    Some(AutoOptions {
        wall_cells: n as usize,
        wake: val("--raff-auto-sillage").unwrap_or(1.0),
        edges: val("--raff-auto-bords").unwrap_or(0.0) > 0.0,
        sparse: val("--raff-auto-creux").unwrap_or(0.0) > 0.0,
    })
}

/// Veine raffinée d'après `o` (plan affiché sur la sortie d'erreur).
pub fn build(body: &Body3D, nu: f64, u: f32, fp16: bool, sa: bool, o: &AutoOptions) -> Result<crate::gpu3d::Gpu3D, String> {
    let ap = auto_plan(body, o)?;
    eprint!("Raffinement automatique ({} mailles par longueur à la paroi) :\n{}", o.wall_cells, ap.describe());
    crate::gpu3d::Gpu3D::with_plan(body, nu, u, fp16, sa, &ap.plan)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decalages_emboites() {
        // L = 8, paroi ×8 : couche ×8 de 0,25 L (16 mailles), boîtes ×4 et
        // ×2 à 4 mailles parentes de plus chacune.
        let d = offsets(8.0, 3, 2, None);
        assert!((d[2] - 0.25).abs() < 1e-12);
        assert!((d[1] - (0.25 + 4.0 / 32.0)).abs() < 1e-12);
        assert!((d[0] - (0.375 + 4.0 / 16.0)).abs() < 1e-12);
        // Couches : règle de `nested_bands` (D_parent ≥ D/2 + 7) en mailles.
        let d = offsets(8.0, 3, 1, None);
        let (b2, b3) = (d[1] * 32.0, d[2] * 64.0);
        assert!(b2 >= 0.5 * b3 + 7.0 - 1e-9);
        // L = 16, ×4 : 16 mailles ×4 = 0,25 L, au-dessus de 0,125 L.
        assert!((offsets(16.0, 2, 1, None)[1] - 0.25).abs() < 1e-12);
        // Croissance avec la taille de maille.
        for w in offsets(8.0, 4, 1, Some(0.0625)).windows(2) {
            assert!(w[0] > w[1]);
        }
    }

    #[test]
    fn plan_aile() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../banc/wing_ar4.stl");
        let tris = crate::stl::read_stl(&path).unwrap();
        let setup = crate::stl::MeshSetup {
            unit: crate::stl::Unit::M,
            alpha_deg: 4.0,
            ref_length: Some(1.0),
            ref_area: Some(4.0),
            margin: 0.4,
            ..Default::default()
        };
        let (body, _) = Body3D::from_mesh(&tris, 16, &setup).unwrap();
        let p = auto_plan(&body, &AutoOptions { wall_cells: 64, ..Default::default() }).unwrap();
        assert_eq!(p.levels.len(), 2);
        assert!(p.levels[0].boxed);
        // Région décalée plus petite que sa boîte, et emboîtement des volumes.
        for lv in &p.levels {
            assert!(lv.region_volume <= lv.box_volume * 1.01, "{lv:?}");
        }
        assert!(p.levels[1].region_volume < p.levels[0].region_volume);
        assert_eq!(p.plan.boxes.len() + p.plan.layers.len(), 2);
        // Bords : bord de fuite et saumons (arêtes vives) au moins.
        let e = edge_set(&body, 2.0 / 4.0);
        assert!(e.edges.len() > 10);
        assert!(auto_plan(&body, &AutoOptions { wall_cells: 48, ..Default::default() }).is_err());
    }
}
