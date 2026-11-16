//! Contrôle géométrique d'une aile maillée : résolution par tranche
//! d'envergure (grille la plus fine), arêtes vives du maillage, et mailles
//! que le critère par maille du solveur classe en arête vive loin de toute
//! arête du maillage.
//! geom3d_check <stl> <unite> <lref> <sref> <L> <raff> <alpha> [gpu 0|1] [lissage 0|1]
use lbm::bodies3d::Body3D;
use lbm::geom3d::{feature_edges, sections, FeatureSet, Indexed, FEATURE_DEG};
use lbm::stl::{MeshSetup, Unit};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let g = |i: usize, d: f64| a.get(i).and_then(|v| v.parse().ok()).unwrap_or(d);
    let mut tris = lbm::stl::read_stl(std::path::Path::new(&a[1])).unwrap();
    if g(9, 0.0) > 0.0 {
        let (t, n) = lbm::geom3d::smooth_tris(&tris);
        println!("raffinement lisse : {n} découpes par arête, {} → {} triangles", tris.len(), t.len());
        tris = t;
    }
    let setup = MeshSetup {
        unit: Unit::parse(&a[2]).unwrap(),
        ref_length: Some(g(3, 1.0)),
        ref_area: Some(g(4, 4.0)),
        margin: 0.4,
        alpha_deg: g(7, 4.0),
        ..Default::default()
    };
    let (l, raff) = (g(5, 16.0) as usize, g(6, 2.0) as usize);
    let (body, info) = Body3D::from_mesh(&tris, l, &setup).unwrap();
    let mesh = body.mesh.as_ref().unwrap();
    let s = (1u64 << raff) as f64;
    // Maquette dans la grille la plus fine (mailles fines).
    let fine: Vec<_> = mesh.tris.iter().map(|t| t.map(|p| p.map(|c| c * s))).collect();
    println!("{} triangles, L = {l}, ×{} : {} mailles par longueur de référence", info.triangles, s, l as f64 * s);
    let ix = Indexed::from_tris(&mesh.tris);
    let fe = feature_edges(&ix, FEATURE_DEG);
    let (cv, open) = (fe.iter().filter(|e| e.convex).count(), fe.iter().filter(|e| e.angle_deg >= 180.0).count());
    println!(
        "arêtes vives (> {FEATURE_DEG}°) : {} ({} saillantes, {} rentrantes, {} bords ouverts)",
        fe.len(),
        cv,
        fe.len() - cv - open,
        open
    );
    // Histogramme des angles dièdres de toutes les arêtes.
    let mut hist = [0usize; 9];
    let ef = ix.edge_faces();
    for fs in ef.values() {
        if fs.len() == 2 {
            let n0 = ix.face_normal(fs[0] as usize);
            let n1 = ix.face_normal(fs[1] as usize);
            let c = (n0[0] * n1[0] + n0[1] * n1[1] + n0[2] * n1[2])
                / ((n0.iter().map(|v| v * v).sum::<f64>() * n1.iter().map(|v| v * v).sum::<f64>()).sqrt());
            let ang = c.clamp(-1.0, 1.0).acos().to_degrees();
            hist[((ang / 20.0) as usize).min(8)] += 1;
        }
    }
    println!("angles dièdres par tranches de 20° : {hist:?}");
    println!("y/b    corde   e/c    R_ba   facette  coin°  culot   e<1    e<2    (mailles fines ; e<k : longueur depuis le bord de fuite)");
    let secs = sections(&fine, 12);
    let (y0, y1) = (fine.iter().flatten().map(|p| p[1]).fold(f64::MAX, f64::min), fine.iter().flatten().map(|p| p[1]).fold(f64::MIN, f64::max));
    for c in &secs {
        println!(
            "{:.3}  {:6.1}  {:.3}  {:5.2}  {:6.2}  {:5.0}  {:5.2}  {:5.1}  {:5.1}",
            (c.y - y0) / (y1 - y0),
            c.chord,
            c.thickness,
            c.le_radius,
            c.le_facet,
            c.le_corner_deg,
            c.te_thickness,
            c.thin1,
            c.thin2
        );
    }
    // Raffinement du bord d'attaque (`--ba 1`) : boîte unique contre boîtes
    // par tronçons d'envergure, en mailles de la grille du bord d'attaque.
    let e = 6.0f64;
    let down = e.max(0.1 * l as f64 * s);
    println!("raffinement du bord d'attaque, M mailles de la grille ×{} :", 2.0 * s);
    for n in [1, 2, 4, 8, 16] {
        println!("  {n:2} boîte(s) : {:.1}", 8.0 * lbm::geom3d::le_cover(&fine, n, e, down) / 1e6);
    }
    if g(8, 1.0) == 0.0 {
        return;
    }
    let fset = FeatureSet::new(fe).transformed(|p| p.map(|c| c * s));
    let sim = lbm::gpu3d::Gpu3D::with_refinement(&body, 1e-5, 0.1, true, true, raff).unwrap();
    if let Some((d, w, f)) = sim.finest() {
        println!("grille ×{f} : {}×{}×{} = {:.1} M mailles, {w} pariétales", d[0], d[1], d[2], (d[0] * d[1] * d[2]) as f64 / 1e6);
    }
    let cells = sim.finest_wall_cells();
    let edge: Vec<_> = cells.iter().filter(|c| c.1 < 0.5).collect();
    let mut near = [0usize; 5];
    for c in &edge {
        let d = fset.distance(c.0.map(|v| v * s), 10.0);
        near[match d {
            d if d < 1.0 => 0,
            d if d < 1.5 => 1,
            d if d < 2.0 => 2,
            d if d < 4.0 => 3,
            _ => 4,
        }] += 1;
    }
    println!(
        "{} mailles pariétales, {} d'arête vive ; distance à une arête du maillage < 1, 1–1,5, 1,5–2, 2–4, > 4 mailles : {:?}",
        cells.len(),
        edge.len(),
        near
    );
}
