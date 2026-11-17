//! Boîtes des grilles des bords (diagnostic) : edge_boxes <stl> <unite> <lref> <sref> <L> <raff> <alpha> [ba] [bf]
use lbm::bodies3d::Body3D;
use lbm::gpu3d::Gpu3D;
use lbm::stl::{MeshSetup, Unit};
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let g = |i: usize, d: f64| a.get(i).and_then(|v| v.parse().ok()).unwrap_or(d);
    let tris = lbm::stl::read_stl(std::path::Path::new(&a[1])).unwrap();
    let setup = MeshSetup { unit: Unit::parse(&a[2]).unwrap(), ref_length: Some(g(3, 1.0)), ref_area: Some(g(4, 4.0)), margin: 0.4, alpha_deg: g(7, 4.0), ..Default::default() };
    let (body, _) = Body3D::from_mesh(&tris, g(5, 16.0) as usize, &setup).unwrap();
    println!("maquette x {:.1}..{:.1} y {:.1}..{:.1} z {:.1}..{:.1}", body.bounds.x0, body.bounds.x1, body.bounds.y0, body.bounds.y1, body.bounds.z0, body.bounds.z1);
    let sim = Gpu3D::with_refinement_edges(&body, 1e-5, 0.1, true, true, g(6, 1.0) as usize, g(8, 0.0) > 0.0, g(9, 1.0) > 0.0, &[g(10, 0.0)].into_iter().filter(|&b| b > 0.0).collect::<Vec<f64>>()).unwrap();
    for (k, b) in sim.edge_box_report().iter().enumerate() {
        println!("boîte {k} : {b}");
    }
}
