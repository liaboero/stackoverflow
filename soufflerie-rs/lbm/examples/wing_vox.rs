//! Diagnostic : aile L = 32 directe contre aile raffinée depuis L = 16.
use lbm::bodies3d::{Body3D, Shape3D};
fn count(b: &Body3D) -> (usize, [f64; 3]) {
    let bb = b.bounds;
    let mut n = 0;
    for z in bb.z0.floor() as i64 - 1..=bb.z1.ceil() as i64 + 1 {
        for y in bb.y0.floor() as i64 - 1..=bb.y1.ceil() as i64 + 1 {
            for x in bb.x0.floor() as i64 - 1..=bb.x1.ceil() as i64 + 1 {
                n += b.inside(x as f64, y as f64, z as f64) as usize;
            }
        }
    }
    (n, [bb.x1 - bb.x0, bb.y1 - bb.y0, bb.z1 - bb.z0])
}
fn main() {
    let a = Body3D::new(Shape3D::Wing, 32, 8.0);
    let c = Body3D::new(Shape3D::Wing, 16, 8.0);
    let r = c.refined([20.0, 20.0, 10.0], [200, 200, 100]);
    println!("direct : {:?} centre {:?} aire {}", count(&a), a.center, a.area);
    println!("raffiné : {:?} centre {:?} aire {}", count(&r), r.center, r.area);
}
