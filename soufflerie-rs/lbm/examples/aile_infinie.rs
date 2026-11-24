//! Diagnostic 3D : aile NACA 0012 d'envergure infinie (veine périodique en y),
//! contre les références 2D (TMR, Ladson). Sans traînée induite : la traînée
//! est celle du profil, pression et frottement.
//!
//! aile_infinie <L> <raff> <Re> <alpha> <t> [bandes, ex. 16,16 ou -] [ny] [hauteur en cordes] [essai] [--raff-auto N]
//!
//! L mailles par corde (veine), `raff` grilles fines entières, `bandes` :
//! couches pariétales imbriquées (comme `run3d --bande`). Efforts intégrés
//! sur la surface (triangles de la tranche centrale de 4 mailles : les bords
//! de la grille fine en envergure ne sont pas périodiques et faussent les
//! tranches voisines), échange aux parois de la même tranche, et bilan de
//! quantité de mouvement sur des boîtes 3D autour de la tranche (mesure de
//! référence : l'intégration sur la surface sous-estime la traînée aux
//! faibles incidences), moyennés sur la dernière corde parcourue (4 relevés,
//! `AILE_NREL`). La grille fine ne couvre pas les bords de l'envergure
//! périodique (raccords) : l'aile n'y existe que sur la grille de la veine,
//! dont la circulation contamine la tranche centrale ; prendre ny ≥ 40.
//!
//! Décrochage : Cm au quart de corde (surface, cabreur positif), durée de
//! moyenne `AILE_MOY` cordes (défaut 1) pour les écoulements instationnaires,
//! avec l'écart type des relevés (`AILE_SERIE=1` : chaque relevé), et abscisse
//! du décollement : vitesse tangentielle moyennée dans le temps des nœuds
//! pariétaux de la rangée centrale, par 40e de corde ; le décollement commence
//! à la première tranche de la zone de retour qui touche le bord de fuite
//! (`AILE_DECOL=1` : détail par tranche). La loi de paroi impose un
//! glissement dans le sens de u_τ : la couche pariétale peut rester attachée
//! sous un écoulement de retour ; on relève donc aussi u·t le long de la
//! normale, jusqu'à `HAUT_DECOL` corde, et le décollement « au-dessus de la
//! paroi » (minimum sur la normale négatif).
use lbm::bodies3d::Body3D;
use lbm::gpu3d::Gpu3D;
use lbm::shapes::Naca;
use lbm::stl::{MeshBody, MeshSetup, Tri, Unit};
use std::sync::Arc;

/// Marges des contours du bilan, en cordes : les 4 premières dans la grille
/// entière la plus fine, les suivantes dans la veine.
const MARGES: [f64; 4] = [0.05, 0.1, 0.15, 0.2];

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let g = |i: usize, d: f64| a.get(i).and_then(|v| v.parse().ok()).unwrap_or(d);
    let l = g(1, 16.0) as usize;
    let raff = g(2, 2.0) as usize;
    let re = g(3, 6e6);
    let alpha = g(4, 0.0);
    let t_end = g(5, 12.0);
    let bands: Vec<f64> = a.get(6).map(|s| s.split(',').filter_map(|v| v.parse().ok()).collect()).unwrap_or_default();
    let ny = g(7, 8.0) as usize;
    let haut = g(8, 12.0);
    let essai = g(9, 0.0) as u32;
    let lf = l as f64;

    // Veine : 4 cordes en amont, 8 en aval, `haut` cordes de haut.
    let (nx, nz) = ((12.0 * lf).round() as usize, (haut * lf).round() as usize);
    let naca = Naca::parse("0012").unwrap();
    // Contour (x/c, z/c) en espacement cosinus, extrados du bord de fuite au
    // bord d'attaque puis intrados ; bord de fuite fermé.
    let n = 120;
    let mut prof: Vec<[f64; 2]> = Vec::new();
    for k in 0..n {
        let x = 0.5 * (1.0 + (std::f64::consts::PI * k as f64 / n as f64).cos());
        prof.push([x, naca.thickness(x)]);
    }
    for k in 0..n {
        let x = 0.5 * (1.0 - (std::f64::consts::PI * k as f64 / n as f64).cos());
        prof.push([x, -naca.thickness(x)]);
    }
    // Rotation autour du quart de corde (bord d'attaque relevé pour α > 0),
    // en mailles.
    let (sa, ca) = alpha.to_radians().sin_cos();
    let rot = |p: [f64; 2]| {
        let (x, z) = ((p[0] - 0.25) * lf, p[1] * lf);
        [x * ca + z * sa, -x * sa + z * ca]
    };
    let te = rot([1.0, 0.0]);
    // Bord de fuite à 0,7 maille après un nœud (plateau de portance), centre
    // en hauteur au milieu de la veine.
    let xq = (4.0 * lf + 0.75 * lf).floor() + 0.7 - te[0];
    let zq = 0.5 * nz as f64 + 0.37 + std::env::var("AILE_DZ").ok().and_then(|v| v.parse::<f64>().ok()).unwrap_or(0.0);
    let xq = xq + std::env::var("AILE_DX").ok().and_then(|v| v.parse::<f64>().ok()).unwrap_or(0.0);
    let pts: Vec<[f64; 2]> = prof.iter().map(|&p| {
        let r = rot(p);
        [xq + r[0], zq + r[1]]
    }).collect();
    // Extrusion en stations d'envergure de 2 mailles, de −4 à ny + 4 : les
    // triangles dont le centre tombe dans [0, ny) couvrent exactement la veine.
    let ys: Vec<f64> = (-2..=(ny as i64 / 2 + 2)).map(|k| 2.0 * k as f64).collect();
    let mut tris: Vec<Tri> = Vec::new();
    let m = pts.len();
    for w in ys.windows(2) {
        let (y0, y1) = (w[0], w[1]);
        for k in 0..m {
            let (p, q) = (pts[k], pts[(k + 1) % m]);
            let (a0, b0, a1, b1) = ([p[0], y0, p[1]], [q[0], y0, q[1]], [p[0], y1, p[1]], [q[0], y1, q[1]]);
            tris.push([a0, b0, b1]);
            tris.push([a0, b1, a1]);
        }
    }
    let setup = MeshSetup { unit: Unit::parse("m").unwrap(), margin: 0.2, ..Default::default() };
    let (mut body, _) = Body3D::from_mesh(&tris, l, &setup).unwrap();
    let mesh = MeshBody::new(tris, false);
    body.nx = nx;
    body.ny = ny;
    body.nz = nz;
    body.bounds = lbm::bodies3d::Bounds { x0: mesh.lo[0], x1: mesh.hi[0], y0: 0.0, y1: ny as f64 - 1.0, z0: mesh.lo[2], z1: mesh.hi[2] };
    body.center = [xq, 0.5 * ny as f64, zq];
    body.area = lf * ny as f64;
    body.mesh = Some(Arc::new(mesh));
    // Arête vive du bord de fuite (AILE_ARETE=1), comme `--aretes maillage`.
    if std::env::var("AILE_ARETE").is_ok() {
        let fe = lbm::geom3d::FeatureEdge { a: [te[0] + xq, -8.0, te[1] + zq], b: [te[0] + xq, ny as f64 + 8.0, te[1] + zq], angle_deg: 180.0, convex: true };
        body.features = Some(Arc::new(lbm::geom3d::FeatureSet::new(vec![fe])));
    } else {
        body.features = None;
    }
    // Paroi volumique par surfels (SOUFFLERIE_SURFELS=1) au lieu du rebond de Bouzidi.
    body.surfels = std::env::var("SOUFFLERIE_SURFELS").is_ok_and(|v| v == "1");

    let u = 0.1f32;
    let nu = u as f64 * lf / re;
    let spt = lf / u as f64;
    let ramp = lbm::cpu::Ramp { steps: (2.0 * spt * std::env::var("AILE_RAMPE").ok().and_then(|v| v.parse::<f64>().ok()).unwrap_or(1.0)) as u64, u };
    let t0 = std::time::Instant::now();
    // `--raff-auto N` (en fin de ligne) : raffinement automatique au lieu
    // de `raff` et `bandes`.
    // `--raff-paroi n` : raffinement de paroi par zones de gradient adverse
    // (`lbm::raff_paroi`), n niveaux de plus que la grille de base, choisis
    // sur un premier calcul de la grille de base (`--raff-paroi-t` cordes,
    // défaut min(t, 6)) ; `--raff-paroi-bande b` : épaisseur du niveau le
    // plus fin (mailles ; défaut max(16, 0,125 L)) ; `--raff-paroi-sillage s` :
    // sillage des couches de s fois leur épaisseur (défaut 0) ; AILE_ZONE_Y=w : zones
    // limitées à une largeur w (mailles de la veine) autour de la tranche
    // centrale (coût).
    let opt = |k: &str| a.iter().position(|v| v == k).and_then(|i| a.get(i + 1)).and_then(|v| v.parse::<f64>().ok());
    let paroi_plan = opt("--raff-paroi").map(|nz| {
        let nz = nz as usize;
        assert!(bands.is_empty(), "--raff-paroi : sans couches pariétales de base");
        let t1 = opt("--raff-paroi-t").unwrap_or(t_end.min(6.0));
        let tp = std::time::Instant::now();
        let tris = &body.mesh.as_ref().unwrap().tris;
        // AILE_BETA=fichier : β du premier calcul relu s'il existe (même
        // maquette, incidence et grille de base), sinon écrit.
        let cache = std::env::var("AILE_BETA").ok();
        let lu: Option<Vec<f64>> = cache.as_ref().and_then(|f| std::fs::read_to_string(f).ok()).map(|t| t.split_whitespace().map(|v| v.parse().unwrap_or(f64::NAN)).collect());
        let beta = match lu.filter(|b| b.len() == tris.len()) {
            Some(b) => b,
            None => {
        let mut s0 = Gpu3D::with_refinement_edges(&body, nu, u, std::env::var("AILE_F32").is_err(), lbm::gpu3d::Turbulence::rans_from_env(), raff, false, false, &[]).unwrap();
        s0.periodic_span = true;
        if essai != 0 {
            let cap = std::env::var("AILE_NUCAP").ok().and_then(|v| v.parse::<f64>().ok()).unwrap_or(0.0);
            s0.set_essai(essai, cap * u as f64 * lf / re);
        }
        s0.run(((t1 - 1.0).max(0.0) * spt) as usize, &ramp);
        let mut zs = lbm::raff_paroi::ZoneStats::new(tris.len());
        for _ in 0..4 {
            s0.run((0.25 * spt) as usize, &ramp);
            if let Some(sf) = s0.surface_forces() {
                zs.add(&sf);
            }
        }
        drop(s0);
        let b = zs.beta(tris);
        if let Some(f) = cache.as_ref() {
            let _ = std::fs::write(f, b.iter().map(|v| v.to_string()).collect::<Vec<_>>().join("
"));
        }
        b
            }
        };
        let yc = 2.0 * (ny / 4) as f64;
        let wy: Option<f64> = std::env::var("AILE_ZONE_Y").ok().and_then(|v| v.parse().ok());
        let cen = |t: &lbm::stl::Tri| [0, 1, 2].map(|k| (t[0][k] + t[1][k] + t[2][k]) / 3.0);
        let keep = |i: usize| wy.is_none_or(|w| (cen(&tris[i])[1] - yc).abs() <= 0.5 * w);
        // AILE_ZONE_TOUT=1 : toute la surface (contrôle : couche sans raccord
        // coupant la paroi, comme `--bande`).
        let tout = std::env::var("AILE_ZONE_TOUT").is_ok();
        let mut beta_z: Vec<f64> = if tout { vec![f64::INFINITY; beta.len()] } else { beta.clone() };
        // AILE_ZONE_XC=a,b (diagnostic des raccords qui coupent la paroi) :
        // zone = triangles de x/c dans [a, b], extrados et intrados.
        // AILE_ZONE_PLUS=a,b : idem, ajoutés à la zone du critère.
        let plus = std::env::var("AILE_ZONE_PLUS").is_ok();
        if let Some(r) = std::env::var(if plus { "AILE_ZONE_PLUS" } else { "AILE_ZONE_XC" }).ok().map(|v| v.split(',').filter_map(|t| t.parse::<f64>().ok()).collect::<Vec<_>>()) {
            let (le, tev) = (pts[n], pts[0]);
            let ch2 = (tev[0] - le[0]).powi(2) + (tev[1] - le[1]).powi(2);
            for (i, t) in tris.iter().enumerate() {
                let c = cen(t);
                let x = ((c[0] - le[0]) * (tev[0] - le[0]) + (c[2] - le[1]) * (tev[1] - le[1])) / ch2;
                beta_z[i] = if x >= r[0] && x <= r[1] { f64::INFINITY } else if plus { beta_z[i] } else { 0.0 };
            }
        }
        let zones = lbm::raff_paroi::zones(&beta_z, nz, keep);
        let bands_z = lbm::raff_paroi::bands(lf, raff, nz, opt("--raff-paroi-bande"));
        // Résumé sur la tranche centrale : part de corde couverte par
        // niveau, extrados et intrados, et β le long de l'extrados.
        let (le, tev) = (pts[n], pts[0]);
        let ch2 = (tev[0] - le[0]).powi(2) + (tev[1] - le[1]).powi(2);
        let xc_of = |c: [f64; 3]| ((c[0] - le[0]) * (tev[0] - le[0]) + (c[2] - le[1]) * (tev[1] - le[1])) / ch2;
        let upper = |c: [f64; 3]| (tev[0] - le[0]) * (c[2] - le[1]) - (tev[1] - le[1]) * (c[0] - le[0]) > 0.0;
        let row = |i: usize| (cen(&tris[i])[1] - yc).abs() < 1.0;
        println!("  raffinement de paroi : premier calcul {t1} cordes, {:.0} s ; épaisseurs {:?} mailles", tp.elapsed().as_secs_f64(), bands_z);
        for (j, z) in zones.iter().enumerate() {
            let (mut se, mut si, mut ne, mut ni) = (0.0, 0.0, 0.0, 0.0);
            let (mut xe, mut xi) = ([f64::INFINITY, f64::NEG_INFINITY], [f64::INFINITY, f64::NEG_INFINITY]);
            for (i, t) in tris.iter().enumerate() {
                if !row(i) {
                    continue;
                }
                let c = cen(t);
                let ar = 0.5 * {
                    let e1 = [t[1][0] - t[0][0], t[1][1] - t[0][1], t[1][2] - t[0][2]];
                    let e2 = [t[2][0] - t[0][0], t[2][1] - t[0][1], t[2][2] - t[0][2]];
                    ((e1[1] * e2[2] - e1[2] * e2[1]).powi(2) + (e1[2] * e2[0] - e1[0] * e2[2]).powi(2) + (e1[0] * e2[1] - e1[1] * e2[0]).powi(2)).sqrt()
                };
                let x = xc_of(c);
                let (s_, n_, xr) = if upper(c) { (&mut se, &mut ne, &mut xe) } else { (&mut si, &mut ni, &mut xi) };
                *n_ += ar;
                if z[i] {
                    *s_ += ar;
                    xr[0] = xr[0].min(x);
                    xr[1] = xr[1].max(x);
                }
            }
            println!(
                "  niveau +{} (×{}) : {:.1} % des triangles · tranche centrale : extrados {:.0} % (x/c {:.2}–{:.2}), intrados {:.0} % (x/c {:.2}–{:.2})",
                j + 1,
                1u64 << (raff + j + 1),
                100.0 * z.iter().filter(|&&k| k).count() as f64 / z.len() as f64,
                100.0 * se / ne.max(1e-12),
                xe[0],
                xe[1],
                100.0 * si / ni.max(1e-12),
                xi[0],
                xi[1]
            );
        }
        let mut bins = [[0.0f64; 2]; 20];
        for (i, t) in tris.iter().enumerate() {
            let c = cen(t);
            if row(i) && upper(c) && beta[i].is_finite() {
                let b = ((xc_of(c) * 20.0).floor().max(0.0) as usize).min(19);
                bins[b][0] += beta[i];
                bins[b][1] += 1.0;
            }
        }
        let v: Vec<String> = bins.iter().enumerate().map(|(b, x)| format!("{:.2}:{:.2}", (b as f64 + 0.5) / 20.0, x[0] / x[1].max(1.0))).collect();
        println!("  β extrados (x/c:β moyen) {}", v.join(" "));
        lbm::raff_paroi::plan(raff, &zones, &bands_z, opt("--raff-paroi-sillage").unwrap_or(0.0))
    });
    let mut sim = match (paroi_plan, lbm::raff_auto::options_from_args(&a)) {
        (Some(p), _) => Gpu3D::with_plan(&body, nu, u, std::env::var("AILE_F32").is_err(), lbm::gpu3d::Turbulence::rans_from_env(), &p).unwrap(),
        (None, Some(o)) => lbm::raff_auto::build(&body, nu, u, std::env::var("AILE_F32").is_err(), lbm::gpu3d::Turbulence::rans_from_env(), &o).unwrap(),
        (None, None) => Gpu3D::with_refinement_edges(&body, nu, u, std::env::var("AILE_F32").is_err(), lbm::gpu3d::Turbulence::rans_from_env(), raff, false, false, &bands).unwrap(),
    };
    {
        let gc = sim.grid_cells();
        let cells: usize = gc.iter().map(|c| c.0).sum();
        let upd: f64 = gc.iter().map(|c| c.0 as f64 * (1u64 << c.1) as f64).sum();
        let per: Vec<String> = gc.iter().map(|c| format!("×{} {:.2} M", 1u64 << c.1, c.0 as f64 / 1e6)).collect();
        println!("  mailles (toutes grilles) : {:.2} M ({}), {:.1} M mises à jour par pas de la veine, modèle {}", cells as f64 / 1e6, per.join(", "), upd / 1e6, if sim.keps { "k-ε RNG" } else { "Spalart–Allmaras" });
    }
    if let Ok(f) = std::env::var("AILE_DUMP") {
        // Champ fin (u, w, ρ de la rangée centrale) après AILE_DUMP_T cordes.
        let tt: f64 = std::env::var("AILE_DUMP_T").ok().and_then(|v| v.parse().ok()).unwrap_or(0.5);
        let ramp0 = lbm::cpu::Ramp { steps: (2.0 * spt) as u64, u };
        sim.run((tt * spt) as usize, &ramp0);
        let (mac, d, org, sc) = if std::env::var("AILE_DUMP_VEINE").is_ok() { (sim.macros(), [nx, ny, nz], [0.0; 3], 1.0) } else { sim.finest_macros() };
        let j = ((2.0 * (ny / 4) as f64 - org[1]) * sc).round() as usize;
        let mut out = format!("{} {} {} {} {}
", d[0], d[2], org[0], org[2], sc);
        for k in 0..d[2] {
            for i in 0..d[0] {
                let c = mac[i + d[0] * (j + d[1] * k)];
                out += &format!("{} {} {}
", c[0], c[2], c[3]);
            }
        }
        std::fs::write(&f, out).unwrap();
        return;
    }
    if std::env::var("AILE_VOX").is_ok() {
        // Empreinte des mailles pariétales fines relative au profil (rangée
        // centrale) : nombre et somme des positions, pour comparer deux
        // placements de même géométrie fine.
        let (_, _, org, sc) = sim.finest_macros();
        let yc = 2.0 * (ny / 4) as f64;
        let mut v: Vec<(i64, i64)> = sim.finest_wall_cells().iter().filter(|(p, _)| (p[1] - yc).abs() < 0.5 / sc).map(|(p, _)| (((p[0] - xq) * sc).round() as i64, ((p[2] - zq) * sc).round() as i64)).collect();
        v.sort();
        let h = v.iter().fold(0u64, |h, &(a, b)| h.wrapping_mul(1_000_003).wrapping_add((a * 7919 + b) as u64));
        println!("  mailles pariétales fines (rangée y = {yc}) : {} empreinte {h:x} ; origine z {:.3}", v.len(), org[2]);
    }
    {
        // La grille fine ne couvre pas toute l'envergure périodique (marges
        // de raccord) : l'aile n'existe ailleurs que sur la grille de la
        // veine, avec sa propre circulation, qui fausse la tranche centrale
        // (traînée +15 %, portance dépendante du placement : ny ≥ 40).
        let (_, d, org, sc) = sim.finest_macros();
        let (y0, y1) = (org[1], org[1] + (d[1] - 1) as f64 / sc);
        println!("  grille fine en envergure : y de {y0:.1} à {y1:.1} sur [0, {ny}) ({:.0} % de l'envergure)", 100.0 * (y1 - y0) / ny as f64);
    }
    if std::env::var("AILE_GRILLE").is_ok() {
        let (_, d, org, sc) = sim.finest_macros();
        let zl = pts.iter().map(|p| p[1]).fold(f64::INFINITY, f64::min);
        let zh = pts.iter().map(|p| p[1]).fold(f64::NEG_INFINITY, f64::max);
        println!("  grille fine : z de {:.3} à {:.3} ; profil de {:.3} à {:.3} : marges {:.3} dessous, {:.3} dessus", org[2], org[2] + (d[2] - 1) as f64 / sc, zl, zh, zl - org[2], org[2] + (d[2] - 1) as f64 / sc - zh);
    }
    sim.periodic_span = true;
    // Loi de paroi hors équilibre (TBLE) : SOUFFLERIE_TBLE=eq|gp|gpc.
    if let Some(m) = lbm::tble::ParoiTble::from_env() {
        sim.set_paroi_tble(Some(m));
        println!("  loi de paroi : {}", m.nom());
    }
    if essai != 0 {
        // AILE_NUCAP : plafond de ν (multiple de ν) pour les essais 1024 / 2048.
        let cap = std::env::var("AILE_NUCAP").ok().and_then(|v| v.parse::<f64>().ok()).unwrap_or(0.0);
        sim.set_essai(essai, cap * u as f64 * lf / re);
    }
    let moy: f64 = std::env::var("AILE_MOY").ok().and_then(|v| v.parse().ok()).unwrap_or(1.0);
    sim.run(((t_end - moy) * spt) as usize, &ramp);
    let _ = sim.take_forces();
    let _ = sim.take_wall_forces();
    let mut ex = [0.0f64; 2];
    let q = 0.5 * (u as f64).powi(2);
    let nrel: usize = std::env::var("AILE_NREL").ok().and_then(|v| v.parse().ok()).unwrap_or(4);
    let (mut cx_p, mut cx_f, mut cz, mut cm_n) = (0.0, 0.0, 0.0, 0usize);
    let mut cmq = 0.0;
    let mut bx = [[0.0f64; 3]; MARGES.len()];
    // Relevés successifs (surface Cz, Cm ; bilan Cz, Cx à la plus grande marge
    // valide) pour juger de l'instationnarité.
    let mut serie: Vec<[f64; 4]> = Vec::new();
    let mut decol = Decol::new(&mut sim, &pts, n, ny);
    for _ in 0..nrel {
        sim.run((moy * spt / nrel as f64) as usize, &ramp);
        let mut rel = [f64::NAN; 4];
        // Échange aux parois des seules mailles de la tranche centrale.
        {
            let (yc0, yc1) = (2.0 * (ny / 4) as f64 - 2.0, 2.0 * (ny / 4) as f64 + 2.0);
            let w = sim.take_wall_forces();
            let mut f = [0.0f64; 2];
            for (p, v) in w.pos.iter().zip(&w.force) {
                if p[1] >= yc0 && p[1] < yc1 {
                    f[0] += v[0];
                    f[1] += v[2];
                }
            }
            if std::env::var("AILE_XBINS").is_ok() {
                // Diagnostic : Cx d'échange par dixième de corde (tranche centrale).
                let mut b = [0.0f64; 12];
                for (p, v) in w.pos.iter().zip(&w.force) {
                    if p[1] >= yc0 && p[1] < yc1 {
                        let k = (((p[0] - (xq - 0.25 * lf)) / lf * 10.0).floor() + 1.0).clamp(0.0, 11.0) as usize;
                        b[k] += v[0] / (q * lf * (yc1 - yc0));
                    }
                }
                let v: Vec<String> = b.iter().map(|x| format!("{x:+.4}")).collect();
                println!("  Cx par dixième de corde : {}", v.join(" "));
            }
            ex[0] += f[0] / (q * lf * (yc1 - yc0));
            ex[1] += f[1] / (q * lf * (yc1 - yc0));
        }
        let Some(s) = sim.surface_forces() else { continue };
        // Pression et frottement séparés : `loads` les mêlent, on les
        // recompose par la part de chaque triangle dans la veine.
        // Seuls les triangles couverts par une grille fine sont relevés :
        // effort par unité d'envergure sur les bandes relevées (écoulement
        // plan), chaque bande de 2 mailles comptant 2 m triangles.
        let (mut fx, mut fz, mut my, mut cnt) = (0.0, 0.0, 0.0, 0usize);
        // Tranche centrale de 4 mailles seulement : les bords de la grille
        // fine en envergure (non périodiques) faussent les tranches voisines.
        let (yc0, yc1) = (2.0 * (ny / 4) as f64 - 2.0, 2.0 * (ny / 4) as f64 + 2.0);
        for (p, f) in &s.loads {
            if p[1] >= yc0 && p[1] < yc1 {
                fx += f[0];
                fz += f[2];
                // Moment au quart de corde (xq, zq), cabreur positif.
                my += (p[2] - zq) * f[0] - (p[0] - xq) * f[2];
                cnt += 1;
            }
        }
        let span = cnt as f64 / m as f64;
        let qa = q * lf * span.max(1e-9);
        // Part de frottement : celle du total, ramenée à la veine.
        let tot = s.pressure_wall[0] + s.friction[0];
        let fr = if tot.abs() > 0.0 { s.friction[0] / tot } else { 0.0 };
        cx_p += fx * (1.0 - fr) / qa;
        cx_f += fx * fr / qa;
        cz += fz / qa;
        cmq += my / (qa * lf);
        rel[0] = fz / qa;
        rel[1] = my / (qa * lf);
        cm_n += 1;
        // Bilan de quantité de mouvement sur des contours rectangulaires du
        // plan (x, z), moyenné sur les rangées de la tranche centrale :
        // −∮ (ρ u u·n + p n) dl, dans la grille entière la plus fine, plus le
        // flux à travers les faces y de la tranche (les bords non périodiques
        // de la grille fine y entretiennent un écoulement transversal qui
        // emporte de la quantité de mouvement x : sans lui, le bilan plan
        // trouve deux fois moins de traînée). Indépendant de la paroi.
        {
            let (xl, xh) = (pts.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min), pts.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max));
            let (zl, zh) = (pts.iter().map(|p| p[1]).fold(f64::INFINITY, f64::min), pts.iter().map(|p| p[1]).fold(f64::NEG_INFINITY, f64::max));
            let fine = sim.finest_macros();
            for (k, mgc) in MARGES.iter().enumerate() {
                let (mac, d, org, sc) = (&fine.0, fine.1, fine.2, fine.3);
                let (j0, j1) = (((yc0 - org[1]) * sc).ceil() as usize, ((yc1 - org[1]) * sc).floor() as usize);
                let mg = mgc * lf;
                let fx = |v: f64, a: usize| ((v - org[a]) * sc).round() as i64;
                let (x0, x1, z0, z1) = (fx(xl - mg, 0), fx(xh + mg, 0), fx(zl - mg, 2), fx(zh + mg, 2));
                if x0 < 2 || z0 < 2 || x1 > d[0] as i64 - 3 || z1 > d[2] as i64 - 3 || j1 <= j0 {
                    continue;
                }
                // Boîte 3D x0..x1, j0..j1, z0..z1 (trapèzes), flux de
                // quantité de mouvement ρ (u − U e_x) u·n + p n : le débit
                // net à travers la boîte fermée est nul, et la forme en
                // u − U le rend insensible à l'écart de fermeture.
                let uin = u as f64;
                let tw = |i: i64, a: i64, b: i64| if i == a || i == b { 0.5 } else { 1.0 };
                let (mut fl, mut mass) = ([0.0f64; 2], 0.0f64);
                // Détail (AILE_BILAN=2) : part x de chaque famille de faces
                // (x, z, y), convection et pression.
                let mut part = [[0.0f64; 2]; 3];
                let mut add = |x: i64, j: i64, z: i64, n: [f64; 3], w: f64| {
                    let c = mac[x as usize + d[0] * (j as usize + d[1] * z as usize)];
                    let (ux, uy, uz, r) = (c[0] as f64, c[1] as f64, c[2] as f64, c[3] as f64);
                    let un = ux * n[0] + uy * n[1] + uz * n[2];
                    let p = (r - 1.0) / 3.0;
                    fl[0] += w * (r * (ux - uin) * un + p * n[0]);
                    fl[1] += w * (r * uz * un + p * n[2]);
                    mass += w * r * un;
                    let fam = if n[0] != 0.0 { 0 } else if n[2] != 0.0 { 1 } else { 2 };
                    part[fam][0] += w * r * (ux - uin) * un;
                    part[fam][1] += w * p * n[0];
                };
                let (ja, jb) = (j0 as i64, j1 as i64);
                for j in ja..=jb {
                    for z in z0..=z1 {
                        let w = tw(j, ja, jb) * tw(z, z0, z1);
                        add(x0, j, z, [-1.0, 0.0, 0.0], w);
                        add(x1, j, z, [1.0, 0.0, 0.0], w);
                    }
                    for x in x0..=x1 {
                        let w = tw(j, ja, jb) * tw(x, x0, x1);
                        add(x, j, z0, [0.0, 0.0, -1.0], w);
                        add(x, j, z1, [0.0, 0.0, 1.0], w);
                    }
                }
                for z in z0..=z1 {
                    for x in x0..=x1 {
                        let w = tw(z, z0, z1) * tw(x, x0, x1);
                        add(x, ja, z, [0.0, -1.0, 0.0], w);
                        add(x, jb, z, [0.0, 1.0, 0.0], w);
                    }
                }
                let span = (jb - ja) as f64;
                if std::env::var("AILE_BILAN").is_ok_and(|v| v == "2") {
                    let k = -1.0 / span / (q * lf * sc);
                    print!(" [x {:.4}+{:.4} z {:.4}+{:.4} y {:.4}+{:.4}]", k * part[0][0], k * part[0][1], k * part[1][0], k * part[1][1], k * part[2][0], k * part[2][1]);
                }
                if std::env::var("AILE_BILAN").is_ok() {
                    print!(" {:.4}/{:.4} (débit {:.2e})", -fl[1] / span / (q * lf * sc), -fl[0] / span / (q * lf * sc), mass / span / (uin * lf * sc));
                    if k + 1 == MARGES.len() { println!(); }
                }
                bx[k][0] += -fl[0] / span / (q * lf * sc);
                bx[k][1] += -fl[1] / span / (q * lf * sc);
                bx[k][2] += 1.0;
                rel[2] = -fl[1] / span / (q * lf * sc);
                rel[3] = -fl[0] / span / (q * lf * sc);
            }
            decol.relever(&sim, &fine);
        }
        serie.push(rel);
        if std::env::var("AILE_SERIE").is_ok() {
            println!("  relevé {} : surface Cz {:.4} Cm {:.5} · bilan Cz {:.4} Cx {:.5}", serie.len(), rel[0], rel[1], rel[2], rel[3]);
        }
        if std::env::var("AILE_CORDE").is_ok() {
            // Cz par dixième de corde (tranche centrale), extrados/intrados.
            let xl = pts.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min);
            let xh = pts.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max);
            let zc = 0.5 * (pts[0][1] + pts[pts.len() / 2][1]);
            let mut b = [[0.0f64; 2]; 10];
            for (p, f) in &s.loads {
                if p[1] >= yc0 && p[1] < yc1 {
                    let k = (((p[0] - xl) / (xh - xl) * 10.0) as usize).min(9);
                    b[k][if p[2] > zc { 0 } else { 1 }] += f[2] / qa;
                }
            }
            let v: Vec<String> = b.iter().map(|x| format!("{:+.4}/{:+.4}", x[0], x[1])).collect();
            println!("  Cz par dixième de corde (extrados/intrados) {}", v.join(" "));
        }
        if std::env::var("AILE_BANDES").is_ok() {
            // Cx et Cz par bande d'envergure de 2 mailles (veine).
            let nb = ny / 2;
            let mut b = vec![[0.0f64; 3]; nb];
            for (p, f) in &s.loads {
                if p[1] >= 0.0 && p[1] < ny as f64 {
                    let k = ((p[1] / 2.0) as usize).min(nb - 1);
                    b[k][0] += f[0];
                    b[k][1] += f[2];
                    b[k][2] += 1.0;
                }
            }
            let v: Vec<String> = b.iter().map(|x| if x[2] > 0.0 { format!("{:.4}/{:.3}", x[0] / (q * lf * 2.0), x[1] / (q * lf * 2.0)) } else { "-".into() }).collect();
            println!("  bandes : {}", v.join(" "));
        }
        if std::env::var("AILE_DEBUG").is_ok() {
            println!(
                "  totaux : pression_wall {:?} frottement {:?} manqué {:.3} · {} triangles relevés",
                s.pressure_wall.map(|v| v / qa),
                s.friction.map(|v| v / qa),
                s.missed,
                s.loads.len()
            );
        }
    }
    // Durée de chaque passe (AILE_PROFIL=lancements), comme `banc3d --profil` ;
    // modifie le champ : relevés suivants (AILE_P0, AILE_W) sans valeur.
    if let Some(reps) = std::env::var("AILE_PROFIL").ok().and_then(|v| v.parse::<usize>().ok()) {
        let mut tot = 0.0;
        for (name, n, ms) in sim.profile_passes(reps) {
            println!("  passe {name:32} {:8.3} M  {ms:8.3} ms", n as f64 / 1e6);
            tot += ms;
        }
        println!("  somme des passes : {tot:.3} ms");
    }
    if std::env::var("AILE_P0").is_ok() {
        perte_p0(&mut sim, &pts, lf, ny, u as f64, nu);
    }
    if std::env::var("AILE_W").is_ok() {
        // Vitesse verticale moyenne de la veine sur des sections x fixées
        // (en cordes depuis l'entrée), toute la hauteur et l'envergure.
        let mac = sim.macros();
        let mut o = String::from("  w moyen (x/c depuis l'entrée : w/U)");
        for xc in [0.5, 1.0, 2.0, 3.0, 6.0, 8.0, 10.0, 11.5] {
            let i = ((xc * lf) as usize).min(nx - 1);
            let (mut sw, mut n) = (0.0, 0.0);
            for k in 1..nz - 1 {
                for j in 0..ny {
                    sw += mac[i + nx * (j + ny * k)][2] as f64;
                    n += 1.0;
                }
            }
            o += &format!(" {xc}: {:.5}", sw / n / u as f64);
        }
        println!("{o}");
    }
    {
        // Moyenne et écart type des relevés.
        let st = |k: usize| {
            let v: Vec<f64> = serie.iter().map(|r| r[k]).filter(|x| x.is_finite()).collect();
            let m = v.iter().sum::<f64>() / v.len().max(1) as f64;
            (m, (v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / v.len().max(1) as f64).sqrt())
        };
        let (a, b, c, d) = (st(0), st(1), st(2), st(3));
        println!(
            "  {} relevés sur {moy} corde(s) : surface Cz {:.4} ± {:.4}, Cm {:.5} ± {:.5} · bilan (marge {}) Cz {:.4} ± {:.4}, Cx {:.5} ± {:.5}",
            serie.len(),
            a.0,
            a.1,
            b.0,
            b.1,
            MARGES.iter().zip(&bx).filter(|(_, b)| b[2] > 0.0).map(|(m, _)| *m).fold(f64::NAN, f64::max),
            c.0,
            c.1,
            d.0,
            d.1
        );
        decol.rapport();
    }
    let bl: Vec<String> = bx.iter().zip(MARGES).filter(|(b, _)| b[2] > 0.0).map(|(b, m)| format!("{m}: {:.4}/{:.5}", b[1] / b[2], b[0] / b[2])).collect();
    println!("  échange (tranche centrale) Cz/Cx {:.4}/{:.5}", ex[1] / nrel as f64, ex[0] / nrel as f64);
    println!("  bilan (marge en cordes : Cz/Cx) {}", bl.join(" · "));
    let f = sim.take_forces();
    let nn = cm_n.max(1) as f64;
    let (mac, ..) = sim.finest_macros();
    let ok = mac.iter().all(|m| m[0].is_finite());
    println!(
        "L = {l} raff {raff} bandes {:?} Re = {re:e} α = {alpha}° ny = {ny} haut {haut} : surface Cz {:.4} Cx {:.5} (pression {:.5} + frottement {:.5}) Cm {:.5} · échange Cz {:.4} Cx {:.5}{} · {:.0} s",
        bands,
        cz / nn,
        (cx_p + cx_f) / nn,
        cx_p / nn,
        cx_f / nn,
        cmq / nn,
        f.fz / (q * body.area),
        f.fx / (q * body.area),
        if ok { "" } else { " · DIVERGÉ" },
        t0.elapsed().as_secs_f64()
    );
}

/// Diagnostic : où se perd la pression totale. Invariant de Bernoulli
/// isotherme H = c_s² ln ρ + u²/2 (exact pour l'écoulement non visqueux
/// stationnaire du schéma), déficit ρ∞ (H∞ − H) / q. Intégré sur des
/// verticales (x fixé, toute la hauteur de la grille fine, mailles fluides,
/// rangées de la tranche centrale) : sa croissance d'une station à l'autre
/// localise la perte, et sa valeur loin en aval approche la traînée de
/// profil. Puis le déficit moyen au premier rang de mailles fluides (à moins
/// de 1,5 maille fine de la paroi), par tranche de corde, extrados/intrados.
fn perte_p0(sim: &mut Gpu3D, pts: &[[f64; 2]], lf: f64, ny: usize, u: f64, nu: f64) {
    let (_, _, _, sc0) = sim.finest_macros();
    let nu_f = (nu * sc0) as f32;
    let (mac, d, org, sc) = sim.finest_macros();
    let q = 0.5 * u * u;
    let h = |c: &[f32; 4]| {
        let (ux, uy, uz, r) = (c[0] as f64, c[1] as f64, c[2] as f64, c[3] as f64);
        (r.ln() / 3.0 + 0.5 * (ux * ux + uy * uy + uz * uz), r)
    };
    let inside = |x: f64, z: f64| {
        let mut c = false;
        let n = pts.len();
        for k in 0..n {
            let (a, b) = (pts[k], pts[(k + n - 1) % n]);
            if (a[1] > z) != (b[1] > z) && x < (b[0] - a[0]) * (z - a[1]) / (b[1] - a[1]) + a[0] {
                c = !c;
            }
        }
        c
    };
    let dist = |x: f64, z: f64| {
        let mut m = f64::INFINITY;
        let n = pts.len();
        for k in 0..n {
            let (a, b) = (pts[k], pts[(k + 1) % n]);
            let (ex, ez) = (b[0] - a[0], b[1] - a[1]);
            let t = (((x - a[0]) * ex + (z - a[1]) * ez) / (ex * ex + ez * ez)).clamp(0.0, 1.0);
            m = m.min(((x - a[0] - t * ex).powi(2) + (z - a[1] - t * ez).powi(2)).sqrt());
        }
        m
    };
    let (yc0, yc1) = (2.0 * (ny / 4) as f64 - 2.0, 2.0 * (ny / 4) as f64 + 2.0);
    let (j0, j1) = (((yc0 - org[1]) * sc).ceil() as usize, ((yc1 - org[1]) * sc).floor() as usize);
    let xle = pts.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min);
    let xte = pts.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max);
    let ch = xte - xle;
    let at = |i: usize, j: usize, k: usize| &mac[i + d[0] * (j + d[1] * k)];
    // Référence : colonne amont (x = 3) de la grille fine.
    let (mut hinf, mut rinf, mut nn) = (0.0, 0.0, 0.0);
    for j in j0..=j1 {
        for k in 3..d[2] - 3 {
            let (hh, r) = h(at(3, j, k));
            hinf += hh;
            rinf += r;
            nn += 1.0;
        }
    }
    hinf /= nn;
    rinf /= nn;
    let xs = [-0.2, -0.05, 0.0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.98, 1.0, 1.02, 1.05, 1.1, 1.25, 1.5, 2.0];
    let mut line = Vec::new();
    for xc in xs {
        let x = xle + xc * ch;
        let i = ((x - org[0]) * sc).round() as i64;
        if i < 3 || i > d[0] as i64 - 4 {
            continue;
        }
        let i = i as usize;
        let xr = org[0] + i as f64 / sc;
        let mut s = 0.0;
        for j in j0..=j1 {
            for k in 3..d[2] - 3 {
                let zr = org[2] + k as f64 / sc;
                if inside(xr, zr) {
                    continue;
                }
                s += rinf * (hinf - h(at(i, j, k)).0);
            }
        }
        line.push(format!("{xc}: {:.4}", s / ((j1 - j0 + 1) as f64) / sc / (q * lf)));
    }
    println!("  perte de p₀ intégrée (x/c : ∫/(q c)) {}", line.join(" · "));
    // Premier rang de mailles fluides.
    let nb = 20;
    let mut bins = vec![[0.0f64; 4]; nb];
    for j in j0..=j1 {
        for k in 3..d[2] - 3 {
            for i in 3..d[0] - 3 {
                let (xr, zr) = (org[0] + i as f64 / sc, org[2] + k as f64 / sc);
                if xr < xle - 1.0 || xr > xte + 1.0 || inside(xr, zr) {
                    continue;
                }
                let dw = dist(xr, zr) * sc;
                if dw > 1.5 {
                    continue;
                }
                let b = (((xr - xle) / ch * nb as f64).floor().max(0.0) as usize).min(nb - 1);
                // Côté : au-dessus ou au-dessous de la corde.
                let zc = pts[0][1] + (pts[pts.len() / 2][1] - pts[0][1]) * (xr - pts[0][0]) / (pts[pts.len() / 2][0] - pts[0][0]);
                let o = if zr > zc { 0 } else { 2 };
                bins[b][o] += rinf * (hinf - h(at(i, j, k)).0) / q;
                bins[b][o + 1] += 1.0;
            }
        }
    }
    let v: Vec<String> = bins.iter().enumerate().map(|(b, x)| format!("{:.2}:{:.3}/{:.3}", (b as f64 + 0.5) / nb as f64, x[0] / x[1].max(1.0), x[2] / x[3].max(1.0))).collect();
    if let Ok(v) = std::env::var("AILE_COL") {
        // Colonnes x/c ∈ [a, b] : 6 premiers nœuds fluides au-dessus de
        // l'extrados, distance à la paroi (mailles fines) et déficit.
        let ab: Vec<f64> = v.split(',').filter_map(|t| t.parse().ok()).collect();
        let j = (j0 + j1) / 2;
        let ws = sim.finest_wall_state();
        let nt = sim.finest_nu_tilde();
        let wc: std::collections::HashMap<(i64, i64, i64), String> = sim.finest_wall_cells().into_iter().enumerate().map(|(k, (p, e))| {
            let st = ws.get(k).map(|v| format!(" uτ {:.4} uw {:.4}", v[0], v[1])).unwrap_or_default();
            (((p[0] * sc).round() as i64, (p[1] * sc).round() as i64, (p[2] * sc).round() as i64), format!(" e{e:.2}{st}"))
        }).collect();
        let (ia, ib) = (((xle + ab[0] * ch - org[0]) * sc) as usize, ((xle + ab[1] * ch - org[0]) * sc) as usize);
        for i in ia..=ib {
            let xr = org[0] + i as f64 / sc;
            let mut k = d[2] - 4;
            while k > 3 && !inside(xr, org[2] + (k - 1) as f64 / sc) && (org[2] + k as f64 / sc) > zq_mid(pts) {
                k -= 1;
            }
            let mut o = format!("    x/c {:.3} :", (xr - xle) / ch);
            for kk in k..(k + 10).min(d[2] - 3) {
                let zr = org[2] + kk as f64 / sc;
                let c = at(i, j, kk);
                let yr = org[1] + j as f64 / sc;
                let e = wc.get(&((xr * sc).round() as i64, (yr * sc).round() as i64, (zr * sc).round() as i64)).cloned().unwrap_or_default();
                let nv = nt.get(i + d[0] * (j + d[1] * kk)).map(|v| format!(" ν̃/ν {:.0}", v / nu_f)).unwrap_or_default();
                o += &format!(" [{:.2} {:.3} u {:.4}{nv}{e}]", dist(xr, zr) * sc, rinf * (hinf - h(c).0) / q, c[0]);
            }
            println!("{o}");
        }
    }
    if let Ok(v) = std::env::var("AILE_PAROI") {
        // Mailles pariétales de la rangée centrale, intrados (i) ou extrados
        // (e) selon AILE_PAROI : x/c, distance (mailles fines), u_τ, glissement,
        // vitesse du nœud et perte de p₀ (÷ q).
        let j = (j0 + j1) / 2;
        let ws = sim.finest_wall_state();
        let mut rows: Vec<(f64, String)> = Vec::new();
        for (k, (p, _)) in sim.finest_wall_cells().into_iter().enumerate() {
            if ((p[1] - org[1]) * sc).round() as usize != j {
                continue;
            }
            let (i, kk) = (((p[0] - org[0]) * sc).round() as usize, ((p[2] - org[2]) * sc).round() as usize);
            let zc = pts[0][1] + (pts[pts.len() / 2][1] - pts[0][1]) * (p[0] - pts[0][0]) / (pts[pts.len() / 2][0] - pts[0][0]);
            if (p[2] > zc) != (v == "e") {
                continue;
            }
            let c = at(i, j, kk);
            let w = ws.get(k).copied().unwrap_or([0.0; 4]);
            let xc = (p[0] - xle) / ch;
            rows.push((xc, format!("{:.4} {:.3} {:.5} {:.4} {:.4} {:.3}", xc, dist(p[0], p[2]) * sc, w[0], w[1], (c[0] * c[0] + c[2] * c[2]).sqrt(), rinf * (hinf - h(c).0) / q)));
        }
        rows.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        for r in rows {
            println!("  paroi {}", r.1);
        }
    }
    if std::env::var("AILE_TE").is_ok() {
        // Carte du bord de fuite (grille fine, rangée centrale) : '#' dans
        // le profil, maille pariétale : u en millièmes (et '!' si arête vive,
        // cos < 0,5), '.' fluide.
        let j = (j0 + j1) / 2;
        let wc: std::collections::HashMap<(i64, i64), f32> = sim.finest_wall_cells().into_iter().filter(|(p, _)| ((p[1] - org[1]) * sc).round() as usize == j).map(|(p, e)| ((((p[0] - org[0]) * sc).round() as i64, ((p[2] - org[2]) * sc).round() as i64), e)).collect();
        // AILE_CARTE=x/c : carte centrée sur l'intrados à cette abscisse.
        let pte = match std::env::var("AILE_CARTE").ok().and_then(|v| v.parse::<f64>().ok()) {
            Some(xc) => {
                let xr = xle + xc * ch;
                let mut zr = pts[0][1] + (pts[pts.len() / 2][1] - pts[0][1]) * (xr - pts[0][0]) / (pts[pts.len() / 2][0] - pts[0][0]);
                while inside(xr, zr) {
                    zr -= 0.01 / sc;
                }
                [xr, zr]
            }
            None => pts[0],
        };
        let (it, kt) = (((pte[0] - org[0]) * sc).round() as i64, ((pte[1] - org[2]) * sc).round() as i64);
        println!("  bord de fuite en ({:.2}, {:.2}) mailles fines", (pte[0] - org[0]) * sc, (pte[1] - org[2]) * sc);
        for k in (kt - 4..=kt + 4).rev() {
            let mut o = format!("  {:>4.2} ", org[2] + k as f64 / sc - pte[1]);
            for i in it - 10..=it + 6 {
                let (xr, zr) = (org[0] + i as f64 / sc, org[2] + k as f64 / sc);
                let c = at(i as usize, j, k as usize);
                let t = if let Some(e) = wc.get(&(i, k)) {
                    format!("{:>4}{}", (c[0] * 1000.0).round(), if *e < 0.5 { '!' } else { ' ' })
                } else if inside(xr, zr) {
                    "   # ".into()
                } else {
                    format!("{:>4}.", (c[0] * 1000.0).round())
                };
                o += &t;
            }
            println!("{o}");
        }
    }
    println!("  perte de p₀ au premier rang (x/c:extrados/intrados) {}", v.join(" "));
}

fn zq_mid(pts: &[[f64; 2]]) -> f64 {
    0.5 * (pts[0][1] + pts[pts.len() / 2][1])
}

/// Tranches de corde du relevé de décollement.
const NB_DECOL: usize = 40;
/// Hauteur du relevé le long de la normale, en cordes.
const HAUT_DECOL: f64 = 0.08;

/// Mailles pariétales de la rangée centrale de la grille fine, avec la
/// tangente au profil orientée vers le bord de fuite, pour situer le
/// décollement (vitesse tangentielle moyenne négative).
struct Decol {
    /// (indice dans `finest_wall_cells`, indice de maille fine, face : 0
    /// extrados, 1 intrados ; tranche de corde ; tangente (x, z)).
    cells: Vec<(usize, usize, usize, usize, [f64; 2])>,
    /// Par maille : indices des nœuds fluides le long de la normale (une
    /// maille fine par pas).
    cols: Vec<Vec<usize>>,
    /// Sommes par face et tranche : u·t du nœud, du glissement, nombre.
    acc: Vec<[f64; 3]>,
    /// Sommes par face, tranche et hauteur : u·t, nombre.
    hacc: Vec<[f64; 2]>,
    nh: usize,
}

impl Decol {
    fn new(sim: &mut Gpu3D, pts: &[[f64; 2]], n: usize, ny: usize) -> Decol {
        let (_, d, org, sc) = sim.finest_macros();
        let (yc0, yc1) = (2.0 * (ny / 4) as f64 - 2.0, 2.0 * (ny / 4) as f64 + 2.0);
        let (j0, j1) = (((yc0 - org[1]) * sc).ceil() as usize, ((yc1 - org[1]) * sc).floor() as usize);
        let j = (j0 + j1) / 2;
        let xle = pts.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min);
        let xte = pts.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max);
        let m = pts.len();
        let mut cells = Vec::new();
        let mut cols = Vec::new();
        let nh = (HAUT_DECOL * (xte - xle) * sc).round().max(1.0) as usize;
        let inside = |x: f64, z: f64| {
            let mut c = false;
            for k in 0..pts.len() {
                let (a, b) = (pts[k], pts[(k + pts.len() - 1) % pts.len()]);
                if (a[1] > z) != (b[1] > z) && x < (b[0] - a[0]) * (z - a[1]) / (b[1] - a[1]) + a[0] {
                    c = !c;
                }
            }
            c
        };
        for (w, (p, _)) in sim.finest_wall_cells().into_iter().enumerate() {
            if ((p[1] - org[1]) * sc).round() as usize != j {
                continue;
            }
            // Segment le plus proche : extrados (du bord de fuite au bord
            // d'attaque, indices < n) ou intrados.
            let (mut best, mut kb) = (f64::INFINITY, 0);
            for k in 0..m {
                let (a, b) = (pts[k], pts[(k + 1) % m]);
                let (ex, ez) = (b[0] - a[0], b[1] - a[1]);
                let t = (((p[0] - a[0]) * ex + (p[2] - a[1]) * ez) / (ex * ex + ez * ez)).clamp(0.0, 1.0);
                let dd = (p[0] - a[0] - t * ex).powi(2) + (p[2] - a[1] - t * ez).powi(2);
                if dd < best {
                    best = dd;
                    kb = k;
                }
            }
            let (a, b) = (pts[kb], pts[(kb + 1) % m]);
            let face = if kb < n { 0 } else { 1 };
            let s = if face == 0 { -1.0 } else { 1.0 };
            let l = ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2)).sqrt();
            let t = [s * (b[0] - a[0]) / l, s * (b[1] - a[1]) / l];
            let xc = ((p[0] - xle) / (xte - xle)).clamp(0.0, 1.0);
            let bin = ((xc * NB_DECOL as f64) as usize).min(NB_DECOL - 1);
            let (i, kk) = (((p[0] - org[0]) * sc).round() as usize, ((p[2] - org[2]) * sc).round() as usize);
            cells.push((w, i + d[0] * (j + d[1] * kk), face, bin, t));
            // Normale extérieure : extrados (−t_z, t_x), intrados (t_z, −t_x).
            let nn = if face == 0 { [-t[1], t[0]] } else { [t[1], -t[0]] };
            let mut col = Vec::new();
            for h in 1..=nh {
                let (x, z) = (p[0] + nn[0] * h as f64 / sc, p[2] + nn[1] * h as f64 / sc);
                let (ii, kz) = (((x - org[0]) * sc).round(), ((z - org[2]) * sc).round());
                if ii < 0.0 || kz < 0.0 || ii >= d[0] as f64 || kz >= d[2] as f64 {
                    break;
                }
                let (xr, zr) = (org[0] + ii / sc, org[2] + kz / sc);
                col.push(if inside(xr, zr) { usize::MAX } else { ii as usize + d[0] * (j + d[1] * kz as usize) });
            }
            cols.push(col);
        }
        Decol { cells, cols, acc: vec![[0.0; 3]; NB_DECOL * 2], hacc: vec![[0.0; 2]; NB_DECOL * 2 * nh], nh }
    }

    /// Ajoute un relevé : par tranche et par face, somme de u·t du nœud
    /// pariétal et du glissement de la loi de paroi, et nombre de mailles.
    fn relever(&mut self, sim: &Gpu3D, fine: &(Vec<[f32; 4]>, [usize; 3], [f64; 3], f64)) {
        let ws = sim.finest_wall_state();
        for (&(w, idx, face, bin, t), col) in self.cells.iter().zip(&self.cols) {
            for (h, &k) in col.iter().enumerate() {
                if k != usize::MAX {
                    let c = fine.0[k];
                    let a = &mut self.hacc[(face * NB_DECOL + bin) * self.nh + h];
                    a[0] += c[0] as f64 * t[0] + c[2] as f64 * t[1];
                    a[1] += 1.0;
                }
            }
            let c = fine.0[idx];
            let a = &mut self.acc[face * NB_DECOL + bin];
            a[0] += c[0] as f64 * t[0] + c[2] as f64 * t[1];
            if let Some(v) = ws.get(w) {
                a[1] += v[1] as f64 * t[0] + v[3] as f64 * t[1];
            }
            a[2] += 1.0;
        }
    }

    fn rapport(&self) {
        let acc = &self.acc;
        // Minimum sur la normale de u·t moyen, par face et tranche (NaN :
        // aucun relevé).
        let hmin = |face: usize, b: usize| {
            let mut m = f64::NAN;
            for h in 0..self.nh {
                let a = self.hacc[(face * NB_DECOL + b) * self.nh + h];
                if a[1] > 0.0 {
                    let v = a[0] / a[1];
                    if !(v >= m) {
                        m = v;
                    }
                }
            }
            m
        };
        let debut_h = |face: usize| {
            let mut x = 1.0;
            for b in (0..NB_DECOL).rev() {
                let v = hmin(face, b);
                if v.is_nan() {
                    continue;
                }
                if v >= 0.0 {
                    break;
                }
                x = b as f64 / NB_DECOL as f64;
            }
            x
        };
        println!(
            "  décollement x/c au-dessus de la paroi (u·t < 0 jusqu'à {HAUT_DECOL} c) : extrados {:.3} · intrados {:.3}",
            debut_h(0),
            debut_h(1)
        );
        // Début de la zone de retour (u·t < 0) qui touche le bord de fuite
        // (1 : attaché).
        let debut = |face: usize, k: usize| {
            let mut x = 1.0;
            for b in (0..NB_DECOL).rev() {
                let a = acc[face * NB_DECOL + b];
                if a[2] == 0.0 {
                    continue;
                }
                if a[k] / a[2] >= 0.0 {
                    break;
                }
                x = b as f64 / NB_DECOL as f64;
            }
            x
        };
        println!(
            "  décollement x/c (vitesse du nœud pariétal ; glissement) : extrados {:.3} ; {:.3} · intrados {:.3} ; {:.3}",
            debut(0, 0),
            debut(0, 1),
            debut(1, 0),
            debut(1, 1)
        );
        if std::env::var("AILE_DECOL").is_ok() {
            for face in 0..2 {
                let v: Vec<String> = (0..NB_DECOL)
                    .map(|b| {
                        let a = acc[face * NB_DECOL + b];
                        if a[2] > 0.0 {
                            format!("{:.3}:{:+.4}/{:+.4}/{:+.4}", (b as f64 + 0.5) / NB_DECOL as f64, a[0] / a[2], a[1] / a[2], hmin(face, b))
                        } else {
                            "-".into()
                        }
                    })
                    .collect();
                println!("  u·t {} (x/c:nœud/glissement/min. sur la normale) {}", if face == 0 { "extrados" } else { "intrados" }, v.join(" "));
            }
        }
    }
}
