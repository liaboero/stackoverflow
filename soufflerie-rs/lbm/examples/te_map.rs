//! Diagnostic : classement des mailles autour du bord de fuite (mi-envergure).
use lbm::bodies3d::Body3D;
use lbm::stl::{MeshSetup, Unit};
const C: [[i32; 3]; 18] = [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1],[1,1,0],[-1,-1,0],[1,-1,0],[-1,1,0],[1,0,1],[-1,0,-1],[1,0,-1],[-1,0,1],[0,1,1],[0,-1,-1],[0,1,-1],[0,-1,1]];
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let tris = lbm::stl::read_stl(std::path::Path::new(&a[1])).unwrap();
    let off: f64 = a.get(2).and_then(|v| v.parse().ok()).unwrap_or(0.0);
    let setup = MeshSetup { unit: Unit::M, alpha_deg: 4.0, ref_length: Some(1.0), ref_area: Some(4.0), margin: 0.4, offset: [off, 0.0, 0.0], ..Default::default() };
    let (b, _) = Body3D::from_mesh(&tris, 32, &setup).unwrap();
    let m = b.mesh.as_ref().unwrap();
    let bb = b.bounds;
    let y = (0.5 * (bb.y0 + bb.y1)).round() as i32;
    println!("bounds x {:.2}..{:.2} z {:.2}..{:.2}", bb.x0, bb.x1, bb.z0, bb.z1);
    let (xa, xb) = (bb.x1.floor() as i32 - 8, bb.x1.ceil() as i32 + 2);
    for z in ((bb.z0.floor() as i32 - 1)..=(bb.z1.ceil() as i32 + 1)).rev() {
        let mut s = format!("{z:4} ");
        for x in xa..=xb {
            let p = [x as f64, y as f64, z as f64];
            let ch = if b.inside(p[0], p[1], p[2]) { '#' } else {
                let mut hit = false; let mut sol = false;
                for c in C { let d = [-c[0] as f64, -c[1] as f64, -c[2] as f64];
                    if m.segment_hit(p, d).is_some() { hit = true; }
                    if b.inside(p[0]+d[0], p[1]+d[1], p[2]+d[2]) { sol = true; } }
                match (hit, sol) { (_, true) => 'w', (true, false) => 't', _ => '.' }
            };
            s.push(ch);
        }
        println!("{s}");
    }
}
