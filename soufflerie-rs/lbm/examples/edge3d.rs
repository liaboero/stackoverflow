//! Mailles d'arête vive par tranche de corde.
//! edge3d <stl> <unite> <lref> <sref> <L> <raff> <alpha>
use lbm::bodies3d::Body3D;
use lbm::gpu3d::Gpu3D;
use lbm::stl::{MeshSetup, Unit};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let g = |i: usize, d: f64| a.get(i).and_then(|v| v.parse().ok()).unwrap_or(d);
    let tris = lbm::stl::read_stl(std::path::Path::new(&a[1])).unwrap();
    let setup = MeshSetup {
        unit: Unit::parse(&a[2]).unwrap(),
        ref_length: Some(g(3, 1.0)),
        ref_area: Some(g(4, 4.0)),
        margin: 0.4,
        alpha_deg: g(7, 4.0),
        ..Default::default()
    };
    let (l, raff) = (g(5, 16.0) as usize, g(6, 2.0) as usize);
    let (body, _) = Body3D::from_mesh(&tris, l, &setup).unwrap();
    let sim = Gpu3D::with_refinement(&body, 1e-5, 0.1, true, true, raff).unwrap();
    let cells = sim.finest_wall_cells();
    let (y0, y1) = cells.iter().fold((f64::MAX, f64::MIN), |m, c| (m.0.min(c.0[1]), m.1.max(c.0[1])));
    let nb = 40;
    let sb = |y: f64| ((((y - y0) / (y1 - y0 + 1e-9)) * nb as f64) as usize).min(nb - 1);
    let mut strip = vec![(f64::MAX, f64::MIN); nb];
    for c in &cells {
        let k = sb(c.0[1]);
        strip[k] = (strip[k].0.min(c.0[0]), strip[k].1.max(c.0[0]));
    }
    let bins = [0.0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.3, 0.7, 0.9, 0.97, 1.01];
    let mut tot = vec![0usize; bins.len() - 1];
    let mut edge = vec![0usize; bins.len() - 1];
    let mut hist = [0usize; 6];
    for c in &cells {
        let k = sb(c.0[1]);
        if k < 2 || k > nb - 3 {
            continue; // saumons
        }
        let xc = (c.0[0] - strip[k].0) / (strip[k].1 - strip[k].0).max(1e-9);
        let b = bins.windows(2).position(|w| xc >= w[0] && xc < w[1]).unwrap_or(bins.len() - 2);
        tot[b] += 1;
        if c.1 < 0.5 {
            edge[b] += 1;
        }
        if xc < 0.02 {
            let h = ((1.0 - c.1 as f64) / 0.4).floor().clamp(0.0, 5.0) as usize;
            hist[h] += 1;
        }
    }
    println!("{} mailles pariétales, {} d'arête", cells.len(), cells.iter().filter(|c| c.1 < 0.5).count());
    println!("x/c          parois  arêtes");
    for (w, (t, e)) in bins.windows(2).zip(tot.iter().zip(&edge)) {
        println!("{:.3}-{:.3}  {:6}  {:6}", w[0], w[1], t, e);
    }
    println!("bord d'attaque (x/c<0,02), cos min par tranches de 0,4 depuis 1 : {:?}", hist);
}
