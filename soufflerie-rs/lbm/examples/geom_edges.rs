//! Arêtes à plus de 60° d'un STL : angle, angle des arêtes parallèles issues
//! des sommets opposés, position. geom_edges <stl>
use lbm::geom3d::{dihedral_angles, Indexed};
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let tris = lbm::stl::read_stl(std::path::Path::new(&a[1])).unwrap();
    let m = Indexed::from_tris(&tris);
    let ang = dihedral_angles(&m);
    let ef = m.edge_faces();
    let mut around: std::collections::HashMap<u32, Vec<(u32, u32)>> = Default::default();
    for &e in ang.keys() {
        around.entry(e.0).or_default().push(e);
        around.entry(e.1).or_default().push(e);
    }
    let unit = |e: (u32, u32)| {
        let p = m.verts[e.0 as usize];
        let q = m.verts[e.1 as usize];
        let d = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
        let l = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        d.map(|x| x / l)
    };
    let mut v: Vec<_> = ang.iter().filter(|(_, &t)| t > 40.0 && t < 140.0).collect();
    v.sort_by(|a, b| a.1.total_cmp(b.1));
    for (e, th) in v {
        let d = unit(*e);
        let mut sides = Vec::new();
        for &f in &ef[e] {
            let o = m.faces[f as usize].iter().copied().find(|&x| x != e.0 && x != e.1).unwrap();
            let par: Vec<f64> = around[&o]
                .iter()
                .filter(|&&g| (unit(g).iter().zip(&d).map(|(a, b)| a * b).sum::<f64>()).abs() > 0.8)
                .map(|g| (ang[g] * 10.0).round() / 10.0)
                .collect();
            sides.push(par);
        }
        let p = m.verts[e.0 as usize];
        println!("{th:6.1}  ({:.3} {:.3} {:.3})  {:?}", p[0], p[1], p[2], sides);
    }
}
