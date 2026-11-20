//! Bilans des surfels (`lbm::surfel`) sur l'aile TR-1208 et sur une aile
//! NACA 0012 extrudée fermée : aires, Σ A n, volume par les fractions
//! solides, volumes balayés et poids, temps de prétraitement comparé à la
//! voxelisation actuelle (nœuds solides et coupes de Bouzidi), à la grille
//! la plus fine.
//! surfels_check [stl TR-1208] [L 8] [raff 2] [alpha 4]
use std::time::Instant;

use lbm::bodies3d::Body3D;
use lbm::shapes::Naca;
use lbm::stl::{MeshBody, MeshSetup, Tri, Unit};
use lbm::surfel::{build, Level, C19};
use rayon::prelude::*;

fn area(t: &Tri) -> f64 {
    let (u, v) = ([0, 1, 2].map(|a| t[1][a] - t[0][a]), [0, 1, 2].map(|a| t[2][a] - t[0][a]));
    let c = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    0.5 * (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt()
}

fn volume(tris: &[Tri]) -> f64 {
    tris.iter()
        .map(|t| {
            let [a, b, c] = *t;
            (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6.0
        })
        .sum()
}

/// Aile NACA 0012 extrudée fermée (corde 1, envergure `b`, unités du fichier).
fn naca_wing(b: f64) -> Vec<Tri> {
    let naca = Naca::parse("0012").unwrap();
    let n = 80;
    // Stations en cosinus ; bord de fuite fermé (épaisseur nulle en x = 1).
    let xs: Vec<f64> = (0..=n).map(|k| 0.5 * (1.0 - (std::f64::consts::PI * k as f64 / n as f64).cos())).collect();
    let up = |x: f64, y: f64| [x, y, naca.thickness(x)];
    let lo = |x: f64, y: f64| [x, y, -naca.thickness(x)];
    let mut t = Vec::new();
    for w in xs.windows(2) {
        let (x0, x1) = (w[0], w[1]);
        // Extrados, intrados.
        t.push([up(x0, 0.0), up(x1, b), up(x1, 0.0)]);
        t.push([up(x0, 0.0), up(x0, b), up(x1, b)]);
        t.push([lo(x0, 0.0), lo(x1, 0.0), lo(x1, b)]);
        t.push([lo(x0, 0.0), lo(x1, b), lo(x0, b)]);
        // Saumons y = 0 et y = b.
        t.push([up(x0, 0.0), up(x1, 0.0), lo(x1, 0.0)]);
        t.push([up(x0, 0.0), lo(x1, 0.0), lo(x0, 0.0)]);
        t.push([up(x0, b), lo(x1, b), up(x1, b)]);
        t.push([up(x0, b), lo(x0, b), lo(x1, b)]);
    }
    // Triangles dégénérés aux bords d'attaque et de fuite retirés ; sommets
    // remis dans le sens direct vu de l'extérieur.
    t.retain(|x| area(x) > 0.0);
    t.into_iter().map(|[a, b, c]| [a, c, b]).collect()
}

/// Voxelisation actuelle à la grille fine : nœuds solides (parité de rayon)
/// et coupes de Bouzidi des 18 liens de chaque nœud fluide voisin d'un
/// nœud solide. Renvoie (nœuds solides, liens coupés).
fn voxelize(mesh: &MeshBody) -> (usize, usize) {
    let lo = mesh.lo.map(|v| v.floor() as i64 - 1);
    let hi = mesh.hi.map(|v| v.ceil() as i64 + 1);
    let d = [0, 1, 2].map(|a| (hi[a] - lo[a] + 1) as usize);
    let solid: Vec<bool> = (0..d[0] * d[1] * d[2])
        .into_par_iter()
        .map(|c| {
            let (i, j, k) = (c % d[0], (c / d[0]) % d[1], c / (d[0] * d[1]));
            mesh.inside((lo[0] + i as i64) as f64, (lo[1] + j as i64) as f64, (lo[2] + k as i64) as f64)
        })
        .collect();
    let at = |i: i64, j: i64, k: i64| {
        (0..d[0] as i64).contains(&i) && (0..d[1] as i64).contains(&j) && (0..d[2] as i64).contains(&k) && solid[i as usize + d[0] * (j as usize + d[1] * k as usize)]
    };
    let links: usize = (0..d[0] * d[1] * d[2])
        .into_par_iter()
        .filter(|&c| !solid[c])
        .map(|c| {
            let (i, j, k) = ((c % d[0]) as i64, ((c / d[0]) % d[1]) as i64, (c / (d[0] * d[1])) as i64);
            let p = [(lo[0] + i) as f64, (lo[1] + j) as f64, (lo[2] + k) as f64];
            C19[1..]
                .iter()
                .filter(|e| at(i + e[0] as i64, j + e[1] as i64, k + e[2] as i64))
                .filter(|e| mesh.segment_hit(p, e.map(|v| v as f64)).is_some())
                .count()
        })
        .sum();
    (solid.iter().filter(|&&s| s).count(), links)
}

fn report(name: &str, tris: &[Tri], raff: u32) {
    let s = (1u64 << raff) as f64;
    let (lo, hi) = lbm::stl::bbox(tris);
    // Grille la plus fine autour de toute la maquette (comme les grilles
    // entières de refine3d) : tous les surfels y sont.
    let mut levels = vec![Level::base()];
    if raff > 0 {
        levels.push(Level {
            scale: s,
            origin: lo.map(|v| v.floor() - 3.0),
            core: [lo.map(|v| v - 2.0), hi.map(|v| v + 2.0)],
        });
    }
    let fine = levels.last().unwrap();
    let ftris: Vec<Tri> = tris.iter().map(|t| t.map(|p| fine.to_level(p))).collect();
    let a0: f64 = ftris.iter().map(area).sum();
    let v0 = volume(&ftris);
    let t0 = Instant::now();
    let ls = build(tris, &levels);
    let t_surf = t0.elapsed().as_secs_f64();
    let l = ls.last().unwrap();
    let t1 = Instant::now();
    let mesh = MeshBody::new(ftris.clone(), false);
    let (nsolid, nlinks) = voxelize(&mesh);
    let t_vox = t1.elapsed().as_secs_f64();
    let an = l.area_vector();
    let anorm = (an[0] * an[0] + an[1] * an[1] + an[2] * an[2]).sqrt();
    let vs = l.fractions.total_solid();
    let ncut = l.fractions.cut().count();
    println!("== {name} : {} triangles, grille ×{s} ({} niveaux)", tris.len(), levels.len());
    println!("aire STL {a0:.4}, Σ surfels {:.4} (écart relatif {:.1e})", l.area(), l.area() / a0 - 1.0);
    println!("|Σ A n| / A = {:.1e}  (Σ A n = [{:.2e}, {:.2e}, {:.2e}])", anorm / a0, an[0], an[1], an[2]);
    println!(
        "volume STL {v0:.4}, Σ fractions solides {vs:.4} (écart relatif {:.1e}) ; nœuds solides {nsolid} ; fuite des colonnes {:.1e}",
        vs / v0 - 1.0,
        l.fractions.leak()
    );
    let nov = l.overlaps.len();
    println!(
        "{} surfels ({:.2} par voxel coupé, {ncut} voxels coupés), {nov} poids ({:.1} par surfel)",
        l.surfels.len(),
        l.surfels.len() as f64 / ncut as f64,
        nov as f64 / l.surfels.len() as f64
    );
    let b = l.balance();
    println!(
        "Σ Pᵢ = {:.4} ; partition : max |Σ poids − Pᵢ| = {:.1e} ; volume balayé hors fluide (recouvrements, solide) {:.4} = {:.2e} de Σ Pᵢ ; dans des voxels pleins {:.2e}",
        b.swept,
        b.partition_err,
        b.excess,
        b.excess / b.swept,
        b.in_solid
    );
    // Part du volume balayé tombant dans des voxels dont le nœud est solide
    // (là où le solveur actuel ne garde pas de populations).
    let mut in_node = 0.0;
    let mut tot = 0.0;
    for k in 0..l.surfels.len() {
        for p in 0..9 {
            for o in l.weights(k, 2 * p + 1) {
                tot += o.vol;
                if mesh.inside(o.voxel[0] as f64, o.voxel[1] as f64, o.voxel[2] as f64) {
                    in_node += o.vol;
                }
            }
        }
    }
    println!("volume balayé dans des voxels à nœud solide : {:.2} %", 100.0 * in_node / tot);
    println!("temps : surfels {t_surf:.2} s ; voxelisation actuelle (nœuds + {nlinks} coupes de liens) {t_vox:.2} s");
    let mem = l.surfels.len() * std::mem::size_of::<lbm::surfel::Surfel>()
        + l.verts.len() * 24
        + nov * std::mem::size_of::<lbm::surfel::Overlap>()
        + l.start.len() * 4;
    println!("mémoire des surfels du niveau fin : {:.1} Mo", mem as f64 / 1e6);
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let g = |i: usize, d: f64| a.get(i).and_then(|v| v.parse().ok()).unwrap_or(d);
    let l = g(2, 8.0) as usize;
    let raff = g(3, 2.0) as u32;
    let alpha = g(4, 4.0);
    if let Some(path) = a.get(1).filter(|p| p.as_str() != "-") {
        let tris = lbm::stl::read_stl(std::path::Path::new(path)).unwrap();
        let setup = MeshSetup {
            unit: Unit::parse("in").unwrap(),
            ref_length: Some(16.672),
            ref_area: Some(2024.35),
            margin: 0.4,
            alpha_deg: alpha,
            ..Default::default()
        };
        let (body, info) = Body3D::from_mesh(&tris, l, &setup).unwrap();
        println!("TR-1208 : maillage {}", if info.diag.watertight() { "fermé".to_string() } else { info.diag.summary() });
        report("TR-1208", &body.mesh.as_ref().unwrap().tris, raff);
    }
    let setup = MeshSetup { unit: Unit::parse("m").unwrap(), alpha_deg: alpha, ..Default::default() };
    let (body, _) = Body3D::from_mesh(&naca_wing(2.0), 2 * l, &setup).unwrap();
    report("NACA 0012 extrudée (envergure 2 c)", &body.mesh.as_ref().unwrap().tris, raff);
}
