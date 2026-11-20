//! Paroi volumique par surfels (Chen, Teixeira et Molvig 1998, « Digital
//! physics approach to computational fluid dynamics » ; celle de PowerFLOW) :
//! préparation sur le processeur des tables du GPU (`gpu3d`, passe `surfel3d`).
//!
//! Géométrie : `surfel::build_axis` (morceaux du maillage par voxel,
//! parallélépipèdes balayés Pᵢ = |cᵢ·n| A et leurs parts exactes dans chaque
//! voxel, volumes fluides V des voxels coupés).
//!
//! Un pas, en populations f (nombres de particules V f) :
//!
//! - rassemblement : le surfel s reçoit, pour chaque direction entrante i
//!   (cᵢ·n < 0), Γᵢⁿ = Σ_y w_{s,i}(y) f*ᵢ(y), parts de son parallélépipède
//!   dans les voxels y, populations après collision ;
//! - réflexion : Γᵒᵘᵗ de la direction sortante opposée = Γⁱⁿ (rebond, masse
//!   et quantité de mouvement normale renvoyées exactement), plus, avec la loi
//!   de paroi, une perturbation de forme d'équilibre (moments 0 et 1) qui
//!   rend la part spéculaire de la quantité de mouvement tangentielle : la
//!   contrainte transmise au corps est ρ u_τ² A (voir `wall_law_matrix`) ;
//! - dispersion et propagation : la maille x reçoit, dans la direction i,
//!   la part transmise de la maille amont y = x − cᵢ, Tᵢ(x) f*ᵢ(y) avec
//!   Tᵢ(x) = V_y − Σ_s w_{s,i}(y) (ce qui n'a touché aucun surfel), plus
//!   Σ_s w_{s,i}(x) Γᵒᵘᵗ_{s,i} / P_{s,i}, le tout divisé par
//!   Dᵢ(x) = Tᵢ(x) + Σ_s w_{s,i}(x).
//!
//! Géométrie exacte : Dᵢ(x) = V_x (ce qui entre dans le volume fluide de x
//! est ce qui sort de celui de y, moins ce qui frappe la paroi, plus ce
//! qu'elle émet), et la masse Σ V f est conservée exactement. Le diviseur
//! Dᵢ(x) plutôt que V_x garde une combinaison convexe (état uniforme au
//! repos conservé, pas d'amplification dans les petits voxels) là où la
//! géométrie discrète s'écarte de l'identité (parallélépipèdes qui se
//! recouvrent dans un coin rentrant, bord de fuite plus mince qu'un voxel) ;
//! l'écart est relevé (`Stats`).

use std::collections::HashMap;

use rayon::prelude::*;

use crate::stl::Tri;
use crate::surfel::{self, Level, C19};

/// Type d'une maille solide (comme `gpu3d::SOLID`).
pub const SOLID: u32 = 0xffff_ffff;
/// Volume fluide en dessous duquel un voxel coupé est solide (zéro numérique
/// de la découpe, pas un réglage).
pub const VMIN: f64 = 1e-9;

pub fn w19(i: usize) -> f64 {
    let c = C19[i];
    match c[0].abs() + c[1].abs() + c[2].abs() {
        0 => 1.0 / 3.0,
        1 => 1.0 / 18.0,
        _ => 1.0 / 36.0,
    }
}

/// Direction opposée (ordre de `C19` : 2p + 1 ↔ 2p + 2).
pub fn opp19(i: usize) -> usize {
    if i == 0 {
        0
    } else if i % 2 == 1 {
        i + 1
    } else {
        i - 1
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cvec(i: usize) -> [f64; 3] {
    C19[i].map(|v| v as f64)
}

/// Valeurs et vecteurs propres d'une matrice symétrique 4 × 4 (Jacobi).
fn jacobi4(mut a: [[f64; 4]; 4]) -> ([f64; 4], [[f64; 4]; 4]) {
    let mut v = [[0.0; 4]; 4];
    for (k, r) in v.iter_mut().enumerate() {
        r[k] = 1.0;
    }
    for _ in 0..60 {
        let off: f64 = (0..4).flat_map(|p| (p + 1..4).map(move |q| (p, q))).map(|(p, q)| a[p][q] * a[p][q]).sum();
        if off < 1e-30 {
            break;
        }
        for p in 0..4 {
            for q in p + 1..4 {
                if a[p][q].abs() < 1e-300 {
                    continue;
                }
                let th = 0.5 * (a[q][q] - a[p][p]) / a[p][q];
                let t = th.signum() / (th.abs() + (th * th + 1.0).sqrt());
                let t = if th == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..4 {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..4 {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for r in v.iter_mut() {
                    let (vp, vq) = (r[p], r[q]);
                    r[p] = c * vp - s * vq;
                    r[q] = s * vp + c * vq;
                }
            }
        }
    }
    ([a[0][0], a[1][1], a[2][2], a[3][3]], v)
}

/// Loi de paroi : perturbation des populations sortantes d'un surfel de
/// normale `n`, δⱼ = Pⱼ wⱼ (λ₀ + Λ·cⱼ) (Pⱼ = |cⱼ·n| A, directions sortantes),
/// qui ne change ni la masse ni la quantité de mouvement normale renvoyées
/// et ajoute ΔM (tangentiel) à la quantité de mouvement sortante :
/// (λ₀, Λ) = M ΔM / A. Rend M (colonnes : ΔM selon x, y, z ; composante
/// normale de ΔM ignorée). G = Σ |c·n| w φ φᵀ, φ = (1, c), est singulière
/// pour une normale d'axe (c·n = 1 pour toutes les sortantes : masse et
/// quantité de mouvement normale ne font qu'une contrainte) : pseudo-inverse.
pub fn wall_law_matrix(n: [f64; 3]) -> [[f64; 4]; 3] {
    let mut g = [[0.0; 4]; 4];
    for i in 1..19 {
        let c = cvec(i);
        let cn = dot(c, n);
        if cn <= 0.0 {
            continue;
        }
        let phi = [1.0, c[0], c[1], c[2]];
        for a in 0..4 {
            for b in 0..4 {
                g[a][b] += cn * w19(i) * phi[a] * phi[b];
            }
        }
    }
    let (ev, v) = jacobi4(g);
    let emax = ev.iter().fold(0.0f64, |m, &e| m.max(e.abs()));
    let mut gp = [[0.0; 4]; 4];
    for k in 0..4 {
        if ev[k].abs() <= 1e-10 * emax {
            continue;
        }
        for a in 0..4 {
            for b in 0..4 {
                gp[a][b] += v[a][k] * v[b][k] / ev[k];
            }
        }
    }
    // Colonne a : G⁺ (0, (I − n nᵀ) e_a).
    let mut m = [[0.0; 4]; 3];
    for (a, col) in m.iter_mut().enumerate() {
        for b in 0..3 {
            let pba = (if a == b { 1.0 } else { 0.0 }) - n[a] * n[b];
            for k in 0..4 {
                col[k] += gp[k][1 + b] * pba;
            }
        }
    }
    m
}

/// Écarts de la géométrie discrète, relevés à la préparation.
#[derive(Clone, Copy, Debug, Default)]
pub struct Stats {
    pub surfels: usize,
    pub receivers: usize,
    /// Volume balayé tombé hors de la grille ou dans un voxel solide.
    pub dropped: f64,
    /// Volume rassemblé en trop (parallélépipèdes qui se recouvrent dans un
    /// voxel au-delà de son volume fluide), retiré par mise à l'échelle.
    pub over_gathered: f64,
    /// max et somme de |V_x − Dᵢ(x)| sur les mailles réceptrices.
    pub max_mismatch: f64,
    pub sum_mismatch: f64,
    /// Directions sans aucune source (Dᵢ = 0) : rebond simple du nœud.
    pub fallback: usize,
    /// Volume transmis vers une maille solide (somme sur les directions de
    /// V_y − Σ w, y fluide, y + cᵢ solide) : perdu (nul pour une géométrie
    /// exacte ; les bords d'une grille fine relèvent des raccords).
    pub leak: f64,
}

/// Tables de la paroi volumique d'une grille boîte.
pub struct SurfelSetup {
    pub dims: [usize; 3],
    /// Volume fluide de chaque maille (rang x + nx (y + ny z)) : 1 fluide, 0 solide.
    pub vol: Vec<f32>,
    /// Mailles réceptrices (rang croissant) : leur propagation passe par les surfels.
    pub wall_cells: Vec<u32>,
    /// Par maille réceptrice et direction 1..18 : poids de la population
    /// propagée T/D, ou −1 (aucune source : rebond simple).
    pub recv_a: Vec<f32>,
    /// Entrées de chaque (maille réceptrice, direction) : `recv[recv_start[18 r + i − 1] ..
    /// recv_start[18 r + i]]`, (case 9 s + p de Γᵒᵘᵗ/P, poids w/D).
    pub recv_start: Vec<u32>,
    pub recv: Vec<(u32, f32)>,
    /// Surfels gardés : normale, aire, centre (coordonnées de la grille).
    pub normal: Vec<[f64; 3]>,
    pub area: Vec<f64>,
    pub center: Vec<[f64; 3]>,
    /// Rassemblement de la paire p du surfel s : `gath[gath_start[9 s + p] .. gath_start[9 s + p + 1]]`,
    /// (maille, poids w).
    pub gath_start: Vec<u32>,
    pub gath: Vec<(u32, f32)>,
    /// Effort de la pression de référence (ρ = 1, repos) par surfel : nul en
    /// somme sur une surface fermée entière.
    pub rest_force: Vec<[f64; 3]>,
    /// Par surfel et paire : 1 / Σ des poids de dispersion gardés (mailles
    /// réceptrices) ; Γᵒᵘᵗ fois cette valeur est la densité dispersée, et la
    /// masse émise est rendue exactement (0 : rien de gardé).
    pub scat_inv: Vec<f64>,
    pub stats: Stats,
}

/// Tables d'une grille de dimensions `dims` pour le maillage `tris`
/// (coordonnées de la grille, normales sortantes). `wrap` : axes
/// périodiques (sinon, hors de la grille : perdu) ; `ghost` : couronne de
/// mailles fantômes (grille fine imbriquée), ni réceptrices ni calculées.
/// Colonnes des fractions solides selon l'axe de plus grande aire projetée.
pub fn prepare(tris: &[Tri], dims: [usize; 3], wrap: [bool; 3], ghost: bool) -> SurfelSetup {
    let [nx, ny, nz] = dims;
    let ns = nx * ny * nz;
    // Axe des colonnes : plus grande aire projetée Σ |A_a|.
    let mut proj = [0.0f64; 3];
    for t in tris {
        let (u, v) = ([0, 1, 2].map(|a| t[1][a] - t[0][a]), [0, 1, 2].map(|a| t[2][a] - t[0][a]));
        let cr = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
        for a in 0..3 {
            proj[a] += cr[a].abs();
        }
    }
    let ax = (0..3).fold(0, |b, a| if proj[a] > proj[b] { a } else { b });
    // Orientation : volume Σ A_ax c_ax (divergence, colonnes selon ax) ;
    // négatif, les normales rentrent (profil extrudé parcouru dans l'autre
    // sens) : triangles retournés.
    let vsign: f64 = tris
        .iter()
        .map(|t| {
            let (u, v) = ([0, 1, 2].map(|a| t[1][a] - t[0][a]), [0, 1, 2].map(|a| t[2][a] - t[0][a]));
            let cr = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
            0.5 * cr[ax] * (t[0][ax] + t[1][ax] + t[2][ax]) / 3.0
        })
        .sum();
    let flipped: Vec<Tri>;
    let tris = if vsign < 0.0 {
        flipped = tris.iter().map(|t| [t[0], t[2], t[1]]).collect();
        &flipped[..]
    } else {
        tris
    };
    let ls = surfel::build_axis(tris, &[Level::base()], ax).remove(0);
    let rank = |x: usize, y: usize, z: usize| x + nx * (y + ny * z);
    let cell_of = |v: [i32; 3]| -> Option<usize> {
        let mut p = [0usize; 3];
        for a in 0..3 {
            let n = dims[a] as i32;
            let q = if wrap[a] { v[a].rem_euclid(n) } else { v[a] };
            if q < 0 || q >= n {
                return None;
            }
            p[a] = q as usize;
        }
        Some(rank(p[0], p[1], p[2]))
    };
    let is_ghost = |c: usize| {
        ghost && {
            let (x, y, z) = (c % nx, (c / nx) % ny, c / (nx * ny));
            x == 0 || y == 0 || z == 0 || x == nx - 1 || y == ny - 1 || z == nz - 1
        }
    };
    // Volumes fluides : dans la boîte englobante du maillage (+ 2), le reste fluide.
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for t in tris {
        for p in t {
            for a in 0..3 {
                lo[a] = lo[a].min(p[a]);
                hi[a] = hi[a].max(p[a]);
            }
        }
    }
    let b0: [usize; 3] = [0, 1, 2].map(|a| ((lo[a] - 2.5).floor().max(0.0) as usize).min(dims[a] - 1));
    let b1: [usize; 3] = [0, 1, 2].map(|a| ((hi[a] + 2.5).ceil().max(0.0) as usize).min(dims[a] - 1));
    let mut vol = vec![1.0f64; ns];
    if tris.iter().len() > 0 {
        vol.par_chunks_mut(nx * ny).enumerate().for_each(|(z, sl)| {
            if z < b0[2] || z > b1[2] {
                return;
            }
            for y in b0[1]..=b1[1] {
                for x in b0[0]..=b1[0] {
                    let f = ls.fractions.fluid([x as i32, y as i32, z as i32]);
                    sl[x + nx * y] = if f < VMIN { 0.0 } else { f.min(1.0) };
                }
            }
        });
    }
    let vol64 = |c: usize| vol[c];
    let mut stats = Stats::default();
    // Surfels gardés : voxel dans la grille (sans repli périodique, une seule
    // image de la surface).
    let keep: Vec<usize> = (0..ls.surfels.len())
        .filter(|&s| (0..3).all(|a| ls.surfels[s].voxel[a] >= 0 && (ls.surfels[s].voxel[a] as usize) < dims[a]))
        .collect();
    // Parts par paire : (maille, w) dans le fluide de la grille.
    let mut pairs: Vec<[Vec<(usize, f64)>; 9]> = Vec::with_capacity(keep.len());
    let mut dropped = 0.0;
    for &s in &keep {
        let mut pr: [Vec<(usize, f64)>; 9] = Default::default();
        for (p, slot) in pr.iter_mut().enumerate() {
            for o in ls.weights(s, 2 * p + 1) {
                match cell_of(o.voxel) {
                    Some(c) if vol[c] > 0.0 => slot.push((c, o.vol)),
                    _ => dropped += o.vol,
                }
            }
        }
        pairs.push(pr);
    }
    stats.dropped = dropped;
    // Directions entrante et sortante de chaque paire.
    let dirs = |s: usize, p: usize| -> Option<(usize, usize)> {
        let cn = dot(cvec(2 * p + 1), ls.surfels[s].normal);
        if cn > 0.0 {
            Some((2 * p + 1, 2 * p + 2))
        } else if cn < 0.0 {
            Some((2 * p + 2, 2 * p + 1))
        } else {
            None
        }
    };
    // Volumes rassemblés (direction entrante) et émis (sortante) par maille.
    let mut win: HashMap<usize, [f64; 19]> = HashMap::new();
    let mut wout: HashMap<usize, [f64; 19]> = HashMap::new();
    for (k, &s) in keep.iter().enumerate() {
        for p in 0..9 {
            let Some((o, i)) = dirs(s, p) else { continue };
            for &(c, w) in &pairs[k][p] {
                win.entry(c).or_insert([0.0; 19])[i] += w;
                if !is_ghost(c) {
                    wout.entry(c).or_insert([0.0; 19])[o] += w;
                }
            }
        }
    }
    // Rassemblement plus grand que le volume fluide : parts réduites
    // d'autant (la maille ne cède pas plus qu'elle n'a).
    let mut gscale: HashMap<(usize, usize), f64> = HashMap::new();
    for (&c, w) in &win {
        for i in 1..19 {
            if w[i] > vol64(c) + 1e-12 {
                stats.over_gathered += w[i] - vol64(c);
                gscale.insert((c, i), vol64(c) / w[i]);
            }
        }
    }
    let win_at = |c: usize, i: usize| -> f64 { win.get(&c).map_or(0.0, |w| w[i]).min(vol64(c)) };
    // Tables du rassemblement.
    let mut gath_start = vec![0u32];
    let mut gath: Vec<(u32, f32)> = Vec::new();
    let mut rest_force = Vec::with_capacity(keep.len());
    for (k, &s) in keep.iter().enumerate() {
        let mut fr = [0.0f64; 3];
        for p in 0..9 {
            if let Some((_, i)) = dirs(s, p) {
                let ci = cvec(i);
                for &(c, w) in &pairs[k][p] {
                    let w = w * gscale.get(&(c, i)).copied().unwrap_or(1.0);
                    gath.push((c as u32, w as f32));
                    // Repos : Γⁱⁿ = Γᵒᵘᵗ = wᵢ Σ w, effort 2 cᵢ Γⁱⁿ.
                    for a in 0..3 {
                        fr[a] += 2.0 * ci[a] * w19(i) * w;
                    }
                }
            }
            gath_start.push(gath.len() as u32);
        }
        rest_force.push(fr);
    }
    // Entrées de dispersion par (maille, direction sortante) : (case 9 k + p, w).
    let mut scat: HashMap<usize, Vec<(usize, u32, f64)>> = HashMap::new();
    for (k, &s) in keep.iter().enumerate() {
        for p in 0..9 {
            let Some((o, _)) = dirs(s, p) else { continue };
            for &(c, w) in &pairs[k][p] {
                if !is_ghost(c) {
                    scat.entry(c).or_default().push((o, (9 * k + p) as u32, w));
                }
            }
        }
    }
    // Mailles candidates : fluides, calculées, dans la boîte (+ 2), dont une
    // source amont est solide, coupée ou rassemblée, ou qui reçoivent des surfels.
    let candidates: Vec<usize> = (b0[2]..=b1[2])
        .into_par_iter()
        .flat_map_iter(|z| {
            let mut out = Vec::new();
            for y in b0[1]..=b1[1] {
                for x in b0[0]..=b1[0] {
                    let c = rank(x, y, z);
                    if vol[c] == 0.0 || is_ghost(c) {
                        continue;
                    }
                    let mut need = wout.contains_key(&c) || vol[c] < 1.0;
                    if !need {
                        for i in 1..19 {
                            let ci = C19[i];
                            let Some(yc) = cell_of([x as i32 - ci[0], y as i32 - ci[1], z as i32 - ci[2]]) else { continue };
                            if vol[yc] < 1.0 || win.get(&yc).is_some_and(|w| w[i] > 0.0) {
                                need = true;
                                break;
                            }
                        }
                    }
                    if need {
                        out.push(c);
                    }
                }
            }
            out
        })
        .collect();
    // Somme des poids de dispersion gardés par (surfel, paire).
    let mut scat_sum = vec![0.0f64; 9 * keep.len()];
    for v in scat.values() {
        for &(_, k, w) in v {
            scat_sum[k as usize] += w;
        }
    }
    let scat_inv: Vec<f64> = scat_sum.iter().map(|&s| if s > 0.0 { 1.0 / s } else { 0.0 }).collect();
    // Volume transmis vers le solide.
    stats.leak = (b0[2]..=b1[2])
        .into_par_iter()
        .map(|z| {
            let mut l = 0.0;
            for y in b0[1]..=b1[1] {
                for x in b0[0]..=b1[0] {
                    let c = rank(x, y, z);
                    if vol[c] == 0.0 {
                        continue;
                    }
                    for i in 1..19 {
                        let ci = C19[i];
                        let dst = cell_of([x as i32 + ci[0], y as i32 + ci[1], z as i32 + ci[2]]);
                        if dst.is_some_and(|d| vol[d] == 0.0) {
                            l += (vol64(c) - win_at(c, i)).max(0.0);
                        }
                    }
                }
            }
            l
        })
        .sum();
    let mut wall_cells = Vec::new();
    let mut recv_a = Vec::new();
    let mut recv_start = vec![0u32];
    let mut recv = Vec::new();
    for c in candidates {
        let (x, y, z) = ((c % nx) as i32, ((c / nx) % ny) as i32, (c / (nx * ny)) as i32);
        let mut a_row = [1.0f32; 18];
        let mut rows: Vec<Vec<(u32, f32)>> = vec![Vec::new(); 18];
        let mut plain = true;
        for i in 1..19 {
            let ci = C19[i];
            let t = match cell_of([x - ci[0], y - ci[1], z - ci[2]]) {
                Some(yc) if vol[yc] > 0.0 => (vol64(yc) - win_at(yc, i)).max(0.0),
                _ => 0.0,
            };
            let ents: Vec<(u32, f64)> = scat.get(&c).map_or(Vec::new(), |v| v.iter().filter(|e| e.0 == i).map(|e| (e.1, e.2)).collect());
            let wo: f64 = ents.iter().map(|e| e.1).sum();
            let d = t + wo;
            if (t - 1.0).abs() < 1e-12 && ents.is_empty() {
                continue;
            }
            plain = false;
            let mm = (vol64(c) - d).abs();
            stats.max_mismatch = stats.max_mismatch.max(mm);
            stats.sum_mismatch += mm;
            if d < 1e-12 {
                a_row[i - 1] = -1.0;
                stats.fallback += 1;
                continue;
            }
            a_row[i - 1] = (t / d) as f32;
            rows[i - 1] = ents.iter().map(|e| (e.0, (e.1 / d) as f32)).collect();
        }
        if plain {
            continue;
        }
        wall_cells.push(c as u32);
        recv_a.extend_from_slice(&a_row);
        for r in rows {
            recv.extend(r);
            recv_start.push(recv.len() as u32);
        }
    }
    stats.surfels = keep.len();
    stats.receivers = wall_cells.len();
    SurfelSetup {
        dims,
        vol: vol.iter().map(|&v| v as f32).collect(),
        wall_cells,
        recv_a,
        recv_start,
        recv,
        normal: keep.iter().map(|&s| ls.surfels[s].normal).collect(),
        area: keep.iter().map(|&s| ls.surfels[s].area).collect(),
        center: keep.iter().map(|&s| ls.surfels[s].center).collect(),
        gath_start,
        gath,
        rest_force,
        scat_inv,
        stats,
    }
}

impl SurfelSetup {
    /// Maille pariétale (rang parmi `wall_cells`) dont chaque surfel prend
    /// u_τ : celle qui reçoit la plus grande part de ses parallélépipèdes
    /// parmi celles où `law` est vrai (loi de paroi), sinon parmi toutes ;
    /// u32::MAX si aucune.
    pub fn home_cells(&self, law: &[bool]) -> Vec<u32> {
        let ns = self.area.len();
        let mut best = vec![(u32::MAX, -1.0f64, false); ns];
        for r in 0..self.wall_cells.len() {
            for k in 0..18 {
                let (e0, e1) = (self.recv_start[18 * r + k] as usize, self.recv_start[18 * r + k + 1] as usize);
                for &(slot, w) in &self.recv[e0..e1] {
                    let s = slot as usize / 9;
                    let (b, bw, bl) = best[s];
                    let l = law[r];
                    // Une maille de la loi de paroi l'emporte toujours.
                    if (l && !bl) || (l == bl && (w as f64) > bw) || b == u32::MAX {
                        best[s] = (r as u32, w as f64, l);
                    }
                }
            }
        }
        best.into_iter().map(|b| b.0).collect()
    }

    /// Octets de ces tables sur le GPU.
    pub fn gpu_bytes(&self) -> usize {
        4 * (self.recv_a.len() + self.recv_start.len() + 2 * self.recv.len())
            + 4 * (self.gath_start.len() + 2 * self.gath.len())
            + self.area.len() * (8 * 16 + 16 + 12 * 4)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matrice_loi_de_paroi() {
        // Pour des normales quelconques et d'axe : la perturbation n'a ni
        // masse ni quantité de mouvement normale, et donne ΔM tangentiel.
        let ns = [[0.0, 0.0, 1.0], [0.3, -0.5, 0.81], [1.0, 1.0, 0.0], [-0.2, 0.9, -0.1]];
        for n in ns {
            let l = dot(n, n).sqrt();
            let n = n.map(|v| v / l);
            let m = wall_law_matrix(n);
            for dm in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.3, -0.7, 0.2]] {
                let dn = dot(dm, n);
                let dt = [0, 1, 2].map(|a| dm[a] - dn * n[a]);
                let lam = [0, 1, 2, 3].map(|k| (0..3).map(|a| m[a][k] * dm[a]).sum::<f64>());
                let (mut mass, mut mom) = (0.0, [0.0; 3]);
                for i in 1..19 {
                    let c = cvec(i);
                    let cn = dot(c, n);
                    if cn <= 0.0 {
                        continue;
                    }
                    let d = cn * w19(i) * (lam[0] + lam[1] * c[0] + lam[2] * c[1] + lam[3] * c[2]);
                    mass += d;
                    for a in 0..3 {
                        mom[a] += c[a] * d;
                    }
                }
                assert!(mass.abs() < 1e-12, "masse {mass} pour n = {n:?}");
                for a in 0..3 {
                    assert!((mom[a] - dt[a]).abs() < 1e-10, "n = {n:?}, ΔM = {dm:?} : {mom:?}");
                }
            }
        }
    }

    /// Sphère maillée, veine périodique : identité géométrique Dᵢ(x) = V_x et
    /// conservation de la masse par les tables seules (populations au repos
    /// et quelconques).
    #[test]
    fn tables_conservatives() {
        let tris = crate::stl::icosphere(5.3, [12.3, 11.6, 12.9], 3);
        let dims = [26, 25, 26];
        let st = prepare(&tris, dims, [true; 3], false);
        assert!(st.stats.surfels > 100 && st.stats.fallback == 0, "{:?}", st.stats);
        assert!(st.stats.max_mismatch < 1e-9 && st.stats.leak < 1e-9 && st.stats.dropped < 1e-12 && st.stats.over_gathered < 1e-5, "{:?}", st.stats);
        // Normales rentrantes : même découpe après retournement.
        let inv: Vec<Tri> = tris.iter().map(|t| [t[0], t[2], t[1]]).collect();
        let st2 = prepare(&inv, dims, [true; 3], false);
        assert!(st2.wall_cells == st.wall_cells && st2.stats.max_mismatch < 1e-9, "{:?}", st2.stats);
        // Un pas de propagation à populations pseudo-aléatoires : masse V f
        // conservée exactement par rassemblement, rebond et dispersion.
        let ns = dims[0] * dims[1] * dims[2];
        let f = |c: usize, i: usize| 0.05 + 0.01 * (((c * 31 + i * 17) % 23) as f64);
        let mass0: f64 = (0..ns).map(|c| st.vol[c] as f64 * (1..19).map(|i| f(c, i)).sum::<f64>()).sum();
        let nk = st.area.len();
        let mut gout = vec![0.0f64; 9 * nk];
        for k in 0..nk {
            for p in 0..9 {
                let (e0, e1) = (st.gath_start[9 * k + p] as usize, st.gath_start[9 * k + p + 1] as usize);
                let cn = dot(cvec(2 * p + 1), st.normal[k]);
                if cn == 0.0 {
                    continue;
                }
                let i = if cn > 0.0 { 2 * p + 2 } else { 2 * p + 1 };
                let g: f64 = st.gath[e0..e1].iter().map(|&(c, w)| w as f64 * f(c as usize, i)).sum();
                gout[9 * k + p] = g * st.scat_inv[9 * k + p];
            }
        }
        let [nx, ny, _] = dims;
        let mut recv_of = HashMap::new();
        for (r, &c) in st.wall_cells.iter().enumerate() {
            recv_of.insert(c as usize, r);
        }
        let mut mass1 = 0.0;
        for c in 0..ns {
            if st.vol[c] == 0.0 {
                continue;
            }
            let (x, y, z) = ((c % nx) as i32, ((c / nx) % ny) as i32, (c / (nx * ny)) as i32);
            for i in 1..19 {
                let ci = C19[i];
                let (sx, sy, sz) = (x - ci[0], y - ci[1], z - ci[2]);
                let (sx, sy, sz) = (sx.rem_euclid(nx as i32) as usize, sy.rem_euclid(ny as i32) as usize, sz.rem_euclid(dims[2] as i32) as usize);
                let sc = sx + nx * (sy + ny * sz);
                let src = if st.vol[sc] == 0.0 { 0.0 } else { f(sc, i) };
                let v = match recv_of.get(&c) {
                    Some(&r) => {
                        let a = st.recv_a[18 * r + i - 1] as f64;
                        let (e0, e1) = (st.recv_start[18 * r + i - 1] as usize, st.recv_start[18 * r + i] as usize);
                        a * src + st.recv[e0..e1].iter().map(|&(k, w)| w as f64 * gout[k as usize]).sum::<f64>()
                    }
                    None => src,
                };
                mass1 += st.vol[c] as f64 * v;
            }
        }
        assert!((mass1 / mass0 - 1.0).abs() < 2e-6, "masse {mass0} → {mass1}");
    }
}
