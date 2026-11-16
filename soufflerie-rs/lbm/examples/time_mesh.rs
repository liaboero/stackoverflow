//! Diagnostic : coût de la mise en place d'un gros maillage.
fn main() {
    let p = std::env::args().nth(1).unwrap();
    let t = std::time::Instant::now();
    let tris = lbm::stl::read_stl(std::path::Path::new(&p)).unwrap();
    println!("lecture {:.2} s", t.elapsed().as_secs_f64());
    let t = std::time::Instant::now();
    let (_, d) = lbm::stl::clean(tris.clone());
    println!("nettoyage {:.2} s ({d:?})", t.elapsed().as_secs_f64());
    let t = std::time::Instant::now();
    let (b, _) = lbm::bodies3d::Body3D::from_mesh(&tris, 32, &Default::default()).unwrap();
    println!("mise en place L = 32 {:.2} s", t.elapsed().as_secs_f64());
    let t = std::time::Instant::now();
    let g = b.sdf_grid(0.5);
    println!("distance signée {:.2} s ({} points)", t.elapsed().as_secs_f64(), g.d.len());
}
