//! Diagnostic : ce que la grille conserve d'une maquette STL (volume, vue de dessus).
use lbm::bodies3d::Body3D;
use lbm::stl::{MeshSetup, Unit};
fn main() {
    let p = std::env::args().nth(1).unwrap();
    let tris = lbm::stl::read_stl(std::path::Path::new(&p)).unwrap();
    let setup = MeshSetup { unit: Unit::M, alpha_deg: 8.0, ..Default::default() };
    for l in [16usize, 32, 64, 128] {
        let (b, info) = Body3D::from_mesh(&tris, l, &setup).unwrap();
        let bb = b.bounds;
        let (mut solid, mut top) = (0usize, 0usize);
        for y in bb.y0.floor() as i64..=bb.y1.ceil() as i64 {
            for x in bb.x0.floor() as i64..=bb.x1.ceil() as i64 {
                let mut col = false;
                for z in bb.z0.floor() as i64..=bb.z1.ceil() as i64 {
                    if b.inside(x as f64, y as f64, z as f64) {
                        solid += 1;
                        col = true;
                    }
                }
                top += col as usize;
            }
        }
        let h = info.ref_length_m / l as f64;
        let cells = b.nx as f64 * b.ny as f64 * b.nz as f64;
        println!(
            "L {l:3} : maille {:.2} m · veine {}×{}×{} = {:.0} M mailles · volume gardé {:.0} % · vue de dessus {:.0} m² · défauts : {}",
            h, b.nx, b.ny, b.nz, cells / 1e6,
            100.0 * solid as f64 / (info.volume_m3 / (h * h * h)),
            top as f64 * h * h, info.diag.summary()
        );
    }
}
