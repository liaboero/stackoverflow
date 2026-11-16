//! Diagnostic : Σ 2 w_i c_i sur les liens pariétaux (poussée de la pression de référence).
use lbm::bodies3d::Body3D;
use lbm::stl::{MeshSetup, Unit};
const C: [([i32; 3], f64); 18] = [([1,0,0],1./18.),([-1,0,0],1./18.),([0,1,0],1./18.),([0,-1,0],1./18.),([0,0,1],1./18.),([0,0,-1],1./18.),([1,1,0],1./36.),([-1,-1,0],1./36.),([1,-1,0],1./36.),([-1,1,0],1./36.),([1,0,1],1./36.),([-1,0,-1],1./36.),([1,0,-1],1./36.),([-1,0,1],1./36.),([0,1,1],1./36.),([0,-1,-1],1./36.),([0,1,-1],1./36.),([0,-1,1],1./36.)];
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let tris = lbm::stl::read_stl(std::path::Path::new(&a[1])).unwrap();
    let setup = MeshSetup { unit: Unit::M, alpha_deg: 4.0, ref_length: Some(1.0), ref_area: Some(4.0), margin: 0.4, ..Default::default() };
    let (b, _) = Body3D::from_mesh(&tris, 32, &setup).unwrap();
    let m = b.mesh.as_ref().unwrap();
    let bb = b.bounds;
    let (mut s, mut n, mut asym) = ([0.0f64; 3], 0usize, 0usize);
    for z in (bb.z0.floor() as i32 - 2)..=(bb.z1.ceil() as i32 + 2) {
        for y in (bb.y0.floor() as i32 - 2)..=(bb.y1.ceil() as i32 + 2) {
            for x in (bb.x0.floor() as i32 - 2)..=(bb.x1.ceil() as i32 + 2) {
                let p = [x as f64, y as f64, z as f64];
                if b.inside(p[0], p[1], p[2]) { continue; }
                for (c, w) in C {
                    let d = [-c[0] as f64, -c[1] as f64, -c[2] as f64];
                    let q = [p[0] + d[0], p[1] + d[1], p[2] + d[2]];
                    let hit = m.segment_hit(p, d).is_some();
                    let sol = b.inside(q[0], q[1], q[2]);
                    if hit || sol {
                        n += 1;
                        // Quantité de mouvement cédée à la paroi : population sortante selon −c... signe sans importance ici.
                        for k in 0..3 { s[k] += 2.0 * w * d[k]; }
                        if !sol && m.segment_hit(q, c.map(|v| v as f64)).is_none() { asym += 1; }
                    }
                }
            }
        }
    }
    let qd = 0.5 * 0.01 * b.area;
    println!("{n} liens · Σ2wc = ({:.3}, {:.3}, {:.3}) · en coefficients ({:.3}, {:.3}, {:.3}) · liens minces non réciproques {asym}", s[0], s[1], s[2], s[0]/qd, s[1]/qd, s[2]/qd);
}
