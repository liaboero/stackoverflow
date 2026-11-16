//! Couche limite intégrale sur le maillage triangulé d'une maquette, pour le
//! décambrage visqueux du calcul 3D (transpiration, Lighthill 1958).
//!
//! La vitesse extérieure u_e est donnée par triangle (relevée dans le champ
//! du calcul). La couche limite est marchée d'amont en aval sur le graphe des
//! triangles voisins, en volumes finis décentrés amont sur les équations
//! intégrales en lignes de courant, sans écoulement transversal dans la
//! couche (Cousteix, Myring) :
//!
//! - quantité de mouvement : ∇·(u_e² θ ê) + δ* u_e (ê·∇u_e) = τ_p/ρ ;
//! - entraînement (Head) : ∇·(u_e θ H₁ ê) = u_e C_E ;
//!
//! ê = u_e/|u_e|. La divergence des lignes de courant (expansion au bord
//! d'attaque d'une aile en flèche, convergence vers le saumon) y est donc
//! comprise. Fermetures : Thwaites en laminaire, Head (H₁, C_E) et
//! Ludwieg–Tillmann (C_f) en turbulent, comme la méthode des panneaux 2D
//! (`panel.rs`), dont elles sont copiées pour que les deux modules évoluent
//! séparément.
//!
//! Résultat : δ* = Hθ et m = u_e δ* par triangle, et la vitesse normale de
//! transpiration v_n = ∇·(m ê) (la même divergence décentrée), qui impose aux
//! parois du calcul le déplacement de la couche limite.

use crate::stl::Tri;

/// Pas de voisin (bord ouvert, arête non manifold).
const NONE: u32 = u32::MAX;

/// Deux faces plus inclinées que 60° l'une sur l'autre : arête vive (bord de
/// fuite, saumon), comme pour la loi de paroi du calcul.
const CREASE_COS: f64 = 0.5;

/// H turbulent au-delà duquel la couche limite est tenue pour décollée
/// (Head : décollement entre 1,8 et 2,4) ; H y est plafonné.
pub const H_SEP: f64 = 2.4;

/// H de départ de la couche limite turbulente à la transition (panel.rs).
const H_TR: f64 = 1.4;

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

// ------------------------------------------------------------ fermetures

/// Thwaites : H(λ) et ℓ(λ) = C_f Re_θ / 2 (copie de panel.rs).
pub fn thwaites_hl(lam: f64) -> (f64, f64) {
    if lam >= 0.0 {
        let l = lam.min(0.1);
        (2.61 - 3.75 * l + 5.24 * l * l, 0.22 + 1.57 * l - 1.8 * l * l)
    } else {
        let l = lam.max(-0.1);
        (2.088 + 0.0731 / (l + 0.14), 0.22 + 1.402 * l + 0.018 * l / (l + 0.107))
    }
}

/// Head : H₁ = (δ − δ*)/θ en fonction de H (copie de panel.rs).
pub fn head_h1(h: f64) -> f64 {
    if h <= 1.6 {
        3.3 + 0.8234 * (h - 1.1).powf(-1.287)
    } else {
        3.3 + 1.5501 * (h - 0.6778).powf(-3.064)
    }
}

/// Head : H en fonction de H₁ (copie de panel.rs).
pub fn head_h(h1: f64) -> f64 {
    if h1 <= 3.32 {
        3.0
    } else if h1 >= 5.3 {
        0.86 * (h1 - 3.3).powf(-0.777) + 1.1
    } else {
        1.1538 * (h1 - 3.3).powf(-0.326) + 0.6778
    }
}

/// Coefficient d'entraînement de Head.
pub fn head_ce(h1: f64) -> f64 {
    0.0306 * (h1 - 3.0).max(1e-9).powf(-0.6169)
}

/// Ludwieg–Tillmann : C_f(H, Re_θ), Re_θ ≥ 10 comme dans panel.rs.
pub fn ludwieg_tillmann(h: f64, re_theta: f64) -> f64 {
    0.246 * 10f64.powf(-0.678 * h) * re_theta.max(10.0).powf(-0.268)
}

// ------------------------------------------------------------ topologie

/// Arête d'un triangle : voisin, normale à l'arête dans le plan du triangle,
/// sortante, multipliée par la longueur de l'arête, et arête vive.
#[derive(Clone, Copy, Debug)]
pub struct Edge {
    pub nb: u32,
    pub nl: [f64; 3],
    pub crease: bool,
}

/// Triangles, voisinage par les arêtes et rangement pour la recherche du
/// triangle le plus proche.
pub struct Topology {
    pub tris: Vec<Tri>,
    pub c: Vec<[f64; 3]>,
    /// Normale sortante (volume signé positif).
    pub n: Vec<[f64; 3]>,
    pub area: Vec<f64>,
    pub edges: Vec<[Edge; 3]>,
    /// Fraction de la corde locale (tranche en y) au centre de chaque triangle.
    pub chord: Vec<f64>,
    lo: [f64; 3],
    gdim: [usize; 3],
    cell: f64,
    grid: Vec<Vec<u32>>,
}

impl Topology {
    pub fn new(tris: &[Tri]) -> Topology {
        let nt = tris.len();
        let raw: Vec<[f64; 3]> = tris.iter().map(|t| cross(sub(t[1], t[0]), sub(t[2], t[0]))).collect();
        let vol: f64 = tris.iter().zip(&raw).map(|(t, n)| dot(t[0], *n)).sum();
        let sgn = if vol >= 0.0 { 1.0 } else { -1.0 };
        let area: Vec<f64> = raw.iter().map(|n| 0.5 * norm(*n)).collect();
        let n: Vec<[f64; 3]> = raw
            .iter()
            .map(|v| {
                let l = norm(*v);
                if l > 0.0 {
                    v.map(|x| sgn * x / l)
                } else {
                    [0.0; 3]
                }
            })
            .collect();
        let c: Vec<[f64; 3]> = tris.iter().map(|t| [0, 1, 2].map(|k| (t[0][k] + t[1][k] + t[2][k]) / 3.0)).collect();
        // Sommets soudés : mêmes coordonnées (à 10⁻⁶ maille près).
        let key = |p: [f64; 3]| p.map(|v| (v * 1e6).round() as i64);
        let mut ids: std::collections::HashMap<[i64; 3], u32> = std::collections::HashMap::new();
        let vid: Vec<[u32; 3]> = tris
            .iter()
            .map(|t| {
                t.map(|p| {
                    let k = ids.len() as u32;
                    *ids.entry(key(p)).or_insert(k)
                })
            })
            .collect();
        let mut emap: std::collections::HashMap<(u32, u32), Vec<(u32, u8)>> = std::collections::HashMap::new();
        for (ti, v) in vid.iter().enumerate() {
            for e in 0..3 {
                let (a, b) = (v[e], v[(e + 1) % 3]);
                if a != b {
                    emap.entry((a.min(b), a.max(b))).or_default().push((ti as u32, e as u8));
                }
            }
        }
        let mut edges = vec![
            [Edge {
                nb: NONE,
                nl: [0.0; 3],
                crease: false,
            }; 3];
            nt
        ];
        for ti in 0..nt {
            for e in 0..3 {
                let (a, b) = (tris[ti][e], tris[ti][(e + 1) % 3]);
                let ab = sub(b, a);
                let mut nl = cross(ab, n[ti]);
                let mid = [0, 1, 2].map(|k| 0.5 * (a[k] + b[k]));
                if dot(nl, sub(mid, c[ti])) < 0.0 {
                    nl = nl.map(|v| -v);
                }
                // |ab × n| = |ab| : normale déjà multipliée par la longueur.
                edges[ti][e].nl = nl;
            }
        }
        for list in emap.values() {
            if let [(a, ea), (b, eb)] = list[..] {
                let crease = dot(n[a as usize], n[b as usize]) < CREASE_COS;
                edges[a as usize][ea as usize].nb = b;
                edges[a as usize][ea as usize].crease = crease;
                edges[b as usize][eb as usize].nb = a;
                edges[b as usize][eb as usize].crease = crease;
            }
        }
        // Fraction de corde : étendue en x des triangles par tranche d'une maille en y.
        let (lo, hi) = crate::stl::bbox(tris);
        let nb = ((hi[1] - lo[1]).ceil() as usize).max(1) + 1;
        let mut span = vec![(f64::MAX, f64::MIN); nb];
        let bin = |y: f64| (((y - lo[1]).floor().max(0.0)) as usize).min(nb - 1);
        for t in tris {
            for p in t {
                let s = &mut span[bin(p[1])];
                s.0 = s.0.min(p[0]);
                s.1 = s.1.max(p[0]);
            }
        }
        let chord = c
            .iter()
            .map(|p| {
                let s = span[bin(p[1])];
                if s.1 > s.0 {
                    ((p[0] - s.0) / (s.1 - s.0)).clamp(0.0, 1.0)
                } else {
                    0.0
                }
            })
            .collect();
        // Rangement par cases de 2 mailles (triangle dans toutes les cases de sa boîte).
        let cell = 2.0;
        let gdim = [0, 1, 2].map(|a| ((hi[a] - lo[a]) / cell).floor() as usize + 1);
        let mut grid = vec![Vec::new(); gdim[0] * gdim[1] * gdim[2]];
        for (k, t) in tris.iter().enumerate() {
            let (tl, th) = crate::stl::bbox(std::slice::from_ref(t));
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
        Topology {
            tris: tris.to_vec(),
            c,
            n,
            area,
            edges,
            chord,
            lo,
            gdim,
            cell,
            grid,
        }
    }

    /// Triangle le plus proche de `p` (à moins de `cap` mailles) parmi ceux
    /// dont la face extérieure regarde `p` (une paroi mince, au bord de
    /// fuite, ne renvoie pas à la face opposée) ; à défaut, le plus proche.
    pub fn nearest(&self, p: [f64; 3], cap: f64) -> Option<usize> {
        let g0 = [0, 1, 2].map(|a| ((p[a] - cap - self.lo[a]) / self.cell).floor());
        let g1 = [0, 1, 2].map(|a| ((p[a] + cap - self.lo[a]) / self.cell).floor());
        if (0..3).any(|a| g1[a] < 0.0 || g0[a] > (self.gdim[a] - 1) as f64) {
            return None;
        }
        let r = |a: usize, v: f64| v.clamp(0.0, (self.gdim[a] - 1) as f64) as usize;
        let (mut best, mut any) = ((cap * cap, NONE), (cap * cap, NONE));
        for gz in r(2, g0[2])..=r(2, g1[2]) {
            for gy in r(1, g0[1])..=r(1, g1[1]) {
                for gx in r(0, g0[0])..=r(0, g1[0]) {
                    for &k in &self.grid[gx + self.gdim[0] * (gy + self.gdim[1] * gz)] {
                        let d2 = crate::stl::tri_dist2(p, &self.tris[k as usize]);
                        if d2 < any.0 {
                            any = (d2, k);
                        }
                        if d2 < best.0 && dot(sub(p, self.c[k as usize]), self.n[k as usize]) >= 0.0 {
                            best = (d2, k);
                        }
                    }
                }
            }
        }
        let k = if best.1 != NONE { best.1 } else { any.1 };
        (k != NONE).then_some(k as usize)
    }
}

// ------------------------------------------------------------ marche

/// Couche limite par triangle.
#[derive(Clone, Debug, Default)]
pub struct Layer3 {
    pub theta: Vec<f64>,
    pub h: Vec<f64>,
    pub dstar: Vec<f64>,
    /// Défaut de débit m = u_e δ*.
    pub m: Vec<f64>,
    /// Vitesse de transpiration v_n = ∇·(m ê), positive vers le fluide.
    pub vn: Vec<f64>,
    pub laminar: Vec<bool>,
    /// H plafonné à H_SEP (décollement turbulent), et en aval.
    pub separated: Vec<bool>,
    /// Bilan de quantité de mouvement limité par la garde numérique
    /// (décélération plus forte que ne le permet la maille).
    pub limited: Vec<bool>,
    /// Débit de déplacement qui quitte la surface (bords de fuite) : Σ v_n dA.
    pub outflow: f64,
}

/// Échanges de flux entre triangles voisins pour un champ de directions ê.
struct Fluxes {
    /// Par triangle : (amont, flux géométrique ê·n L > 0) entrants.
    inflow: Vec<Vec<(u32, f64)>>,
    /// Par triangle : flux sortant total (vers un voisin ou hors de la surface).
    out: Vec<f64>,
    /// Par triangle : part du flux sortant qui quitte la surface (arête vive
    /// convergente, bord ouvert).
    leave: Vec<f64>,
}

fn fluxes(t: &Topology, e: &[[f64; 3]]) -> Fluxes {
    let nt = t.c.len();
    let mut inflow = vec![Vec::new(); nt];
    let mut out = vec![0.0; nt];
    let mut leave = vec![0.0; nt];
    for a in 0..nt {
        for ed in &t.edges[a] {
            let pa = dot(e[a], ed.nl);
            if ed.nb == NONE {
                if pa > 0.0 {
                    out[a] += pa;
                    leave[a] += pa;
                }
                continue;
            }
            let b = ed.nb as usize;
            // Flux vu de b à travers la même arête.
            let eb = t.edges[b].iter().find(|x| x.nb == a as u32).map_or([0.0; 3], |x| x.nl);
            let pb = dot(e[b], eb);
            if ed.crease && pa > 0.0 && pb > 0.0 {
                // Bord de fuite : les deux faces convergent, le déplacement
                // part dans le sillage.
                out[a] += pa;
                leave[a] += pa;
                continue;
            }
            if ed.crease && pa <= 0.0 && pb <= 0.0 {
                continue;
            }
            let f = 0.5 * (pa - pb);
            if f > 0.0 {
                out[a] += f;
                inflow[b].push((a as u32, f));
            }
        }
    }
    Fluxes { inflow, out, leave }
}

/// Options de la marche : viscosité (unités du maillage) et transition
/// imposée (fraction de corde ; 0 : turbulent dès le point d'arrêt).
#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub nu: f64,
    pub xtr: f64,
}

/// Couche limite intégrale pour la vitesse extérieure `ue` (par triangle,
/// tangente à la surface).
pub fn march(t: &Topology, ue: &[[f64; 3]], o: &Options) -> Layer3 {
    let nt = t.c.len();
    let u: Vec<f64> = ue.iter().map(|v| norm(*v)).collect();
    let e: Vec<[f64; 3]> = ue.iter().zip(&u).map(|(v, &m)| if m > 1e-12 { v.map(|x| x / m) } else { [0.0; 3] }).collect();
    let fx = fluxes(t, &e);
    // Gradient de |u_e| (Green–Gauss, moyenne aux arêtes) : ê·∇u et |∇u|.
    let mut g = vec![0.0; nt];
    let mut gabs = vec![0.0; nt];
    for a in 0..nt {
        if t.area[a] <= 0.0 {
            continue;
        }
        let mut gr = [0.0; 3];
        for ed in &t.edges[a] {
            let uf = if ed.nb == NONE || ed.crease { u[a] } else { 0.5 * (u[a] + u[ed.nb as usize]) };
            for k in 0..3 {
                gr[k] += uf * ed.nl[k];
            }
        }
        let gr = gr.map(|v| v / t.area[a]);
        g[a] = dot(e[a], gr);
        gabs[a] = norm(gr);
    }
    // Ordre de la marche : temps de parcours depuis les sources (triangles
    // sans amont : lignes de partage au bord d'attaque, bords d'entrée),
    // Dijkstra sur les arêtes où le flux passe d'un triangle à l'autre ; pas
    // de cycle même si des boucles de flux existent (décollement, bruit).
    let mut down: Vec<Vec<u32>> = vec![Vec::new(); nt];
    for b in 0..nt {
        for &(a, _) in &fx.inflow[b] {
            down[a as usize].push(b as u32);
        }
    }
    let mut time = vec![f64::INFINITY; nt];
    let mut heap = std::collections::BinaryHeap::new();
    for a in 0..nt {
        if fx.inflow[a].is_empty() {
            time[a] = 0.0;
            heap.push(std::cmp::Reverse((0u64, a as u32)));
        }
    }
    let mut seen = vec![false; nt];
    while let Some(std::cmp::Reverse((_, a))) = heap.pop() {
        let a = a as usize;
        if seen[a] {
            continue;
        }
        seen[a] = true;
        for &b in &down[a] {
            let b = b as usize;
            let dt = norm(sub(t.c[b], t.c[a])) / (0.5 * (u[a] + u[b])).max(1e-9);
            let tb = time[a] + dt;
            if tb < time[b] {
                time[b] = tb;
                heap.push(std::cmp::Reverse((tb.to_bits(), b as u32)));
            }
        }
    }
    // Tri topologique (Kahn) : un triangle après tous ses amonts ; parmi les
    // prêts, le plus tôt atteint. Boucle de flux (aucun prêt) : rompue au
    // triangle restant atteint le plus tôt, ses amonts non traités ignorés.
    let key = |a: usize| if time[a].is_finite() { time[a].to_bits() } else { u64::MAX };
    let mut indeg: Vec<usize> = fx.inflow.iter().map(|v| v.len()).collect();
    let mut ready = std::collections::BinaryHeap::new();
    let mut all = std::collections::BinaryHeap::new();
    for a in 0..nt {
        all.push(std::cmp::Reverse((key(a), a as u32)));
        if indeg[a] == 0 {
            ready.push(std::cmp::Reverse((key(a), a as u32)));
        }
    }
    let mut order = Vec::with_capacity(nt);
    let mut placed = vec![false; nt];
    while order.len() < nt {
        let a = match ready.pop() {
            Some(std::cmp::Reverse((_, a))) => a as usize,
            None => loop {
                let std::cmp::Reverse((_, a)) = all.pop().expect("triangle restant");
                if !placed[a as usize] {
                    break a as usize;
                }
            },
        };
        if placed[a] {
            continue;
        }
        placed[a] = true;
        order.push(a);
        for &b in &down[a] {
            let b = b as usize;
            indeg[b] -= 1;
            if indeg[b] == 0 && !placed[b] {
                ready.push(std::cmp::Reverse((key(b), b as u32)));
            }
        }
    }

    let mut ly = Layer3 {
        theta: vec![0.0; nt],
        h: vec![0.0; nt],
        dstar: vec![0.0; nt],
        m: vec![0.0; nt],
        vn: vec![0.0; nt],
        laminar: vec![false; nt],
        separated: vec![false; nt],
        limited: vec![false; nt],
        outflow: 0.0,
    };
    let mut h1 = vec![0.0; nt];
    let mut done = vec![false; nt];
    let nu = o.nu;
    for &b in &order {
        done[b] = true;
        let (ub, ar, fo) = (u[b].max(1e-9), t.area[b], fx.out[b]);
        if ar <= 0.0 || u[b] <= 1e-9 {
            continue;
        }
        // Apports des triangles amont déjà traités.
        let (mut min, mut ein, mut fin, mut lam_in, mut sep_in) = (0.0, 0.0, 0.0, 0.0, 0.0);
        let (mut th_in, mut h_in) = (0.0, 0.0);
        for &(a, f) in &fx.inflow[b] {
            let a = a as usize;
            if !done[a] || !(ly.theta[a] > 0.0) {
                continue;
            }
            let ua = u[a];
            min += f * ua * ua * ly.theta[a];
            // Transition : la couche turbulente part de H = 1,4 (panel.rs).
            let h1a = if ly.laminar[a] { head_h1(H_TR) } else { h1[a] };
            ein += f * ua * ly.theta[a] * h1a;
            fin += f;
            th_in += f * ly.theta[a];
            h_in += f * ly.h[a];
            if ly.laminar[a] {
                lam_in += f;
            }
            if ly.separated[a] {
                sep_in += f;
            }
        }
        let source = fin <= 0.0;
        let laminar_up = source || lam_in > 0.5 * fin;
        let want_lam = source || (laminar_up && o.xtr > 0.0 && t.chord[b] < o.xtr);
        // Sortie presque nulle (point de convergence) : état amont recopié.
        if !source && fo < 1e-6 * ar.sqrt() {
            let th = th_in / fin;
            let h = h_in / fin;
            ly.theta[b] = th;
            ly.h[b] = h;
            h1[b] = head_h1(h.min(H_SEP));
            ly.laminar[b] = laminar_up;
            ly.separated[b] = sep_in > 0.5 * fin;
            continue;
        }
        let mut lam_ok = want_lam;
        if want_lam {
            // Thwaites en volumes finis : (u²F_s + A H u g) θ² − M_e θ − A u ℓ ν = 0.
            // Sans amont (M_e = 0) au point d'arrêt (u → 0, g = a) :
            // θ² = ℓν/(Ha) ≈ 0,09 ν/a, la solution de Hiemenz ; au bord
            // d'une plaque, θ² ≈ 2ℓ ν Δx / u, celle de Blasius.
            let mut th = if source { (nu / gabs[b].max(1e-12)).sqrt().min(ar.sqrt()) } else { th_in / fin };
            let mut hl = (2.61, 0.22);
            for _ in 0..30 {
                let lam = (th * th * g[b] / nu).clamp(-0.1, 0.1);
                hl = thwaites_hl(lam);
                let a2 = ub * ub * fo + ar * hl.0 * ub * g[b];
                // λ < −0,09 : décollement laminaire (ℓ s'annule vers −0,09).
                if a2 <= 0.0 || hl.1 <= 0.0 || lam < -0.09 {
                    lam_ok = false;
                    break;
                }
                let tn = (min + (min * min + 4.0 * a2 * ar * ub * hl.1 * nu).sqrt()) / (2.0 * a2);
                th = 0.5 * th + 0.5 * tn;
            }
            let lam = th * th * g[b] / nu;
            if lam_ok && (source || lam >= -0.09) {
                ly.theta[b] = th;
                ly.h[b] = hl.0;
                h1[b] = head_h1(H_TR);
                ly.laminar[b] = true;
            } else if source {
                // Décollement laminaire dès la source : état de Hiemenz.
                ly.theta[b] = (0.075 * nu / gabs[b].max(1e-12)).sqrt();
                ly.h[b] = 2.35;
                h1[b] = head_h1(H_TR);
                ly.laminar[b] = true;
            } else {
                lam_ok = false;
            }
        }
        if !lam_ok && !source {
            // Turbulent : quantité de mouvement pour θ, entraînement pour H₁.
            let mut th = th_in / fin;
            let mut hh = if laminar_up { H_TR } else { (h_in / fin).min(H_SEP) };
            let mut h1v = head_h1(hh);
            let mut limited = false;
            let mut sep = false;
            for _ in 0..30 {
                let cf = ludwieg_tillmann(hh, ub * th / nu);
                let mut den = ub * ub * fo + ar * hh * ub * g[b];
                // Garde numérique : une décélération qui annulerait le
                // dénominateur (décollement dans la maille) double au plus
                // la quantité de mouvement d'une maille à l'autre.
                if den < 0.5 * ub * ub * fo {
                    den = 0.5 * ub * ub * fo;
                    limited = true;
                }
                let tn = (min + ar * ub * ub * 0.5 * cf) / den;
                th = 0.5 * th + 0.5 * tn.max(1e-12);
                // u θ F_s H₁ − A u C_E(H₁) = E_e, croissant en H₁ : bissection.
                let res = |x: f64| ub * th * fo * x - ar * ub * head_ce(x) - ein;
                let (mut lo, mut hi) = (3.0 + 1e-9, 3.3 + 1.5501 * (1.1f64 - 0.6778).powf(-3.064));
                // H₁ maximal : celui de H = 1,1 ; au-delà, la couche est
                // plus pleine qu'aucun profil turbulent.
                if res(hi) < 0.0 {
                    lo = hi;
                } else {
                    for _ in 0..60 {
                        let mid = 0.5 * (lo + hi);
                        if res(mid) > 0.0 {
                            hi = mid;
                        } else {
                            lo = mid;
                        }
                    }
                }
                h1v = 0.5 * (lo + hi);
                hh = head_h(h1v);
            }
            ly.theta[b] = th;
            if hh > H_SEP {
                hh = H_SEP;
                h1v = head_h1(H_SEP);
                sep = true;
            }
            ly.h[b] = hh;
            h1[b] = h1v;
            ly.laminar[b] = false;
            ly.separated[b] = sep || sep_in > 0.5 * fin;
            ly.limited[b] = limited;
        }
    }
    for a in 0..nt {
        if !(ly.theta[a].is_finite() && ly.h[a].is_finite()) {
            ly.theta[a] = 0.0;
            ly.h[a] = 0.0;
            ly.limited[a] = true;
        }
        ly.dstar[a] = ly.h[a] * ly.theta[a];
        ly.m[a] = u[a] * ly.dstar[a];
    }
    let (vn, outflow) = transpiration_fx(t, &fx, &ly.m);
    ly.vn = vn;
    ly.outflow = outflow;
    ly
}

/// Vitesse de transpiration v_n = ∇·(m ê) par triangle pour le défaut de
/// débit `m` (décentré amont, mêmes flux que la marche), et le débit qui
/// quitte la surface aux bords de fuite.
pub fn transpiration(t: &Topology, ue: &[[f64; 3]], m: &[f64]) -> (Vec<f64>, f64) {
    let e: Vec<[f64; 3]> = ue
        .iter()
        .map(|v| {
            let n = norm(*v);
            if n > 1e-12 {
                v.map(|x| x / n)
            } else {
                [0.0; 3]
            }
        })
        .collect();
    transpiration_fx(t, &fluxes(t, &e), m)
}

fn transpiration_fx(t: &Topology, fx: &Fluxes, m: &[f64]) -> (Vec<f64>, f64) {
    let nt = t.c.len();
    let mut vn = vec![0.0; nt];
    let mut outflow = 0.0;
    for b in 0..nt {
        if t.area[b] <= 0.0 {
            continue;
        }
        let mut s = fx.out[b] * m[b];
        for &(a, f) in &fx.inflow[b] {
            s -= f * m[a as usize];
        }
        vn[b] = s / t.area[b];
        outflow += fx.leave[b] * m[b];
    }
    (vn, outflow)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Plaque plane x ∈ [0, lx], y ∈ [0, ly], normale +z, carrés de côté `d`.
    fn plate(lx: f64, ly: f64, d: f64) -> Vec<Tri> {
        let (nx, ny) = ((lx / d).round() as usize, (ly / d).round() as usize);
        let mut t = Vec::new();
        for j in 0..ny {
            for i in 0..nx {
                let p = |a: usize, b: usize| [a as f64 * d, b as f64 * d, 0.0];
                t.push([p(i, j), p(i + 1, j), p(i + 1, j + 1)]);
                t.push([p(i, j), p(i + 1, j + 1), p(i, j + 1)]);
            }
        }
        t
    }

    fn at_x(t: &Topology, v: &[f64], x: f64) -> f64 {
        let (mut s, mut n) = (0.0, 0.0);
        for (k, c) in t.c.iter().enumerate() {
            if (c[0] - x).abs() < 0.5 && c[1] > 2.0 && c[1] < 8.0 {
                s += v[k];
                n += 1.0;
            }
        }
        s / n
    }

    #[test]
    fn plaque_laminaire_blasius() {
        let tris = plate(400.0, 10.0, 1.0);
        let t = Topology::new(&tris);
        let (u, nu) = (0.1, 0.1 * 400.0 / 2e5);
        let ue = vec![[u, 0.0, 0.0]; tris.len()];
        let l = march(&t, &ue, &Options { nu, xtr: 1.0 });
        for x in [100.0, 300.0] {
            let th = at_x(&t, &l.theta, x);
            let blasius = 0.664 * (nu * x / u).sqrt();
            assert!((th / blasius - 1.0).abs() < 0.04, "x = {x} : θ = {th}, Blasius {blasius}");
            assert!((at_x(&t, &l.h, x) - 2.61).abs() < 1e-6);
        }
    }

    #[test]
    fn plaque_turbulente_et_bilan_de_debit() {
        let tris = plate(400.0, 10.0, 1.0);
        let t = Topology::new(&tris);
        let (u, nu) = (0.1, 0.1 * 400.0 / 4e6);
        let ue = vec![[u, 0.0, 0.0]; tris.len()];
        let l = march(&t, &ue, &Options { nu, xtr: 0.0 });
        for x in [100.0, 300.0] {
            // Loi en puissance 1/7 : θ = 0,036 x Re_x^(−1/5).
            let th = at_x(&t, &l.theta, x);
            let p7 = 0.036 * x * (u * x / nu).powf(-0.2);
            assert!((th / p7 - 1.0).abs() < 0.15, "x = {x} : θ = {th}, 1/7 : {p7}");
            let h = at_x(&t, &l.h, x);
            assert!(h > 1.25 && h < 1.5, "H = {h}");
        }
        // Le débit injecté est celui qui quitte la surface au bord aval.
        let inj: f64 = l.vn.iter().zip(&t.area).map(|(v, a)| v * a).sum();
        assert!((inj / l.outflow - 1.0).abs() < 1e-9, "{inj} {}", l.outflow);
        assert!(l.separated.iter().all(|s| !s));
    }

    #[test]
    fn point_d_arret_hiemenz() {
        // Écoulement de point d'arrêt u = a x sur une plaque, de part et
        // d'autre de x = 0 : θ = 0,29 √(ν/a) (Hiemenz, θ²a/ν = 0,085).
        let tris = plate(40.0, 10.0, 0.5);
        let t = Topology::new(&tris);
        let (a, nu) = (1e-3, 1e-4);
        let ue: Vec<[f64; 3]> = t.c.iter().map(|c| [a * (c[0] - 20.0), 0.0, 0.0]).collect();
        let l = march(&t, &ue, &Options { nu, xtr: 1.0 });
        let th = at_x(&t, &l.theta, 25.0);
        let hiemenz = (0.085 * nu / a).sqrt();
        assert!((th / hiemenz - 1.0).abs() < 0.15, "θ = {th}, Hiemenz {hiemenz}");
    }
}
