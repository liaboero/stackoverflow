//! Soufflerie de poche, version native.
//!
//! Sans argument : interface graphique. Sous-commandes :
//!   soufflerie run    [options]          un calcul, moyennes et Strouhal
//!   soufflerie polar  --alpha a:b:pas    polaire (incidences calculées en parallèle), CSV
//!   soufflerie panel  --alpha a:b:pas    polaire instantanée : panneaux + couche limite, CSV
//!   soufflerie bench  [--l 32]           vitesse des moteurs (millions de mailles par seconde)
//!   soufflerie wing   --alpha a:b:pas    aile finie : surface portante + polaires 2D, CSV
//!   soufflerie run3d  [options]          soufflerie 3D D3Q19 sur GPU : efforts moyens
//!
//! Options : --shape cylinder|square|ellipse|plate|naca|triangle  --naca 2412
//!           --alpha 4  --re 1000  --l 32  --t 40  --avg 15  --engine cpu|gpu
//!           --threads 0  --spin 1  --thickness 0.4  --collision rr|bgk
//!           --les auto|on|off  --u 0.1  --turb auto|sa|off  --transition free|turb

use std::collections::HashMap;
use std::time::Instant;

use lbm::cpu::{Collision, LesMode};
use lbm::panel::{compressible_correction, panel_solve, section_coefficients, PanelOptions, PanelResult, PanelShape};
use lbm::stats::{analyze, Analysis};
use lbm::wing::{wing_polar, PolarTable, Sections, Wing, WingParams};
use lbm::{EngineKind, Naca, ShapeKind, ShapeParams, Transition, TunnelConfig, Turbulence, WindTunnel};
use rayon::prelude::*;

#[cfg(feature = "gui")]
mod gui;
#[cfg(feature = "gui")]
mod launch;
#[cfg(feature = "gui")]
mod wing_tab;
#[cfg(feature = "gui")]
mod sweep;
#[cfg(feature = "gui")]
mod tab3d;
#[cfg(feature = "gui")]
mod worker;
#[cfg(feature = "gui")]
mod theme;

fn usage() -> ! {
    eprintln!(
        "Soufflerie de poche (native)\n\n\
         soufflerie [options]             interface graphique (--field eddy, --graph cp|cf)\n\
         soufflerie run   [options]       un calcul\n\
         soufflerie polar --alpha -4:12:2 [options]   polaire, CSV\n\
         soufflerie panel --alpha -4:16:1 [options]   polaire panneaux + couche limite, CSV\n\
         soufflerie bench [--l 32] [--steps 400]       vitesse des moteurs\n\
         soufflerie wing  --alpha -6:20:1 [--span 1]   aile finie (surface portante), CSV\n\
         soufflerie run3d --shape sphere|cube|cylinder|wing | --stl fichier.stl --re 100 --l 12 [--alpha 8] [--t 40] [--avg 15] [--fp32 1]\n\
         \x20                                  soufflerie 3D (D3Q19, GPU)\n\
         \x20         --naca 2412 --ar 8 --taper 0.5 --sweep 0 --dihedral 0 --twist 0 --re 1e6 --strips 40\n\n\
         Options : --shape naca --naca 2412 --alpha 4 --re 1000 --l 32 --t 40 --avg 15\n\
         \x20         --engine cpu|gpu --threads 0 --spin 0 --thickness 0.4\n\
         \x20         --collision rr|bgk --les auto|on|off --u 0.1\n\
         \x20         --turb auto|sa|off (Spalart–Allmaras et loi de paroi dès Re = 5·10⁴ en auto, 2D et 3D)\n\
         \x20         --transition free|turb (libre : prédite par la couche limite intégrale)\n\
         \x20         --mach 0.3 (panel : correction de compressibilité)\n\
         \x20         --xtr 0.3 (run3d, loi de paroi : couche limite laminaire jusqu'à 30 % de la corde locale)\n\
         \x20         --stl : --unite mm|cm|m|in --ecoulement +x --haut +z --alpha 0 --beta 0 (dérapage)\n\
         \x20                 --lref (longueur de réf., unités du fichier) --sref (surface) --vitesse 30 (m/s, fixe Re) --fluide air|eau\n\
         \x20         --raff 1|2 (GPU, 2D et 3D : grilles 2 ou 4 fois plus fines autour de la maquette)\n\
         \x20         --paroi 0 (run, polar : coefficients bruts, sans correction du blocage de la veine)\n\
         \x20         --rendu auto|dx12|vulkan|glow (interface : rendu graphique)"
    );
    std::process::exit(2)
}

pub struct Args {
    cmd: String,
    pub opts: HashMap<String, String>,
}

impl Args {
    fn parse() -> Args {
        let mut it = std::env::args().skip(1).peekable();
        let cmd = match it.peek() {
            Some(s) if !s.starts_with("--") => it.next().unwrap(),
            _ => String::new(),
        };
        let mut opts = HashMap::new();
        while let Some(k) = it.next() {
            let Some(key) = k.strip_prefix("--") else { usage() };
            let v = it.next().unwrap_or_else(|| usage());
            opts.insert(key.to_string(), v);
        }
        Args { cmd, opts }
    }

    fn get<T: std::str::FromStr>(&self, k: &str, d: T) -> T {
        match self.opts.get(k) {
            Some(v) => v.parse().unwrap_or_else(|_| {
                eprintln!("valeur invalide pour --{k} : {v}");
                usage()
            }),
            None => d,
        }
    }

    pub fn config(&self) -> TunnelConfig {
        let shape = ShapeKind::from_id(&self.get("shape", "cylinder".to_string())).unwrap_or_else(|| usage());
        let naca = Naca::parse(&self.get("naca", "0012".to_string())).unwrap_or_else(|| usage());
        let engine = match self.get("engine", "cpu".to_string()).as_str() {
            "gpu" => EngineKind::Gpu,
            _ => EngineKind::Cpu {
                threads: self.get("threads", 0usize),
            },
        };
        TunnelConfig {
            l: self.get("l", 16usize),
            re: self.get("re", 100.0f64),
            u_in: self.get("u", 0.1f32),
            shape,
            params: ShapeParams {
                // Polaire : --alpha est une plage a:b:pas, lue à part.
                alpha_deg: self.opts.get("alpha").and_then(|v| v.parse().ok()).unwrap_or(0.0),
                naca,
                thickness: self.get("thickness", 0.4f64),
                spin_ratio: self.get("spin", 0.0f64),
            },
            les: match self.get("les", "auto".to_string()).as_str() {
                "on" => LesMode::On,
                "off" => LesMode::Off,
                _ => LesMode::Auto,
            },
            collision: if self.get("collision", "rr".to_string()) == "bgk" {
                Collision::Bgk
            } else {
                Collision::Rr
            },
            bouzidi: true,
            engine,
            turbulence: match self.get("turb", "auto".to_string()).as_str() {
                "sa" => Turbulence::Sa,
                "off" | "none" => Turbulence::Off,
                _ => Turbulence::Auto,
            },
            transition: if self.get("transition", "free".to_string()) == "turb" {
                Transition::Turbulent
            } else {
                Transition::Free
            },
            gpu_half: self.get("half", 0u8) == 1,
            refine: self.get("raff", 0usize),
        }
    }
}

/// Calcul jusqu'à `t_end` temps convectifs ; analyse des `avg` derniers.
/// `corr` : coefficients corrigés du blocage de la veine (Allen et Vincenti).
fn simulate(cfg: TunnelConfig, t_end: f64, avg: f64, quiet: bool, corr: bool) -> Result<(Analysis, f64, bool), String> {
    let mut wt = WindTunnel::new(cfg)?;
    let spt = wt.steps_per_tstar();
    let every = ((spt / 40.0).round() as usize).max(1);
    let total = (t_end * spt).round() as usize;
    let (mut t, mut cx, mut cz, mut cm) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let t0 = Instant::now();
    let mut done = 0;
    let mut next_print = 5.0;
    while done < total {
        let n = every.min(total - done);
        wt.run(n);
        done += n;
        let mut polled = wt.engine.poll_forces();
        if done >= total {
            // Derniers relevés encore en vol.
            let s = wt.engine.take_forces();
            if s.steps > 0 {
                polled.push((wt.engine.time(), s));
            }
        }
        let mut bad = false;
        for (step, s) in polled {
            let c = if corr { wt.corrected(wt.coefficients(s)) } else { wt.coefficients(s) };
            let ts = step as f64 / spt;
            if ts >= t_end - avg {
                t.push(ts);
                cx.push(c.cx);
                cz.push(c.cz);
                cm.push(c.cm);
            }
            if !quiet && ts >= next_print {
                eprintln!("  t* = {ts:5.1}  Cx = {:8.4}  Cz = {:8.4}", c.cx, c.cz);
                next_print += 5.0;
            }
            bad |= !c.cx.is_finite();
        }
        if bad {
            break;
        }
    }
    let secs = t0.elapsed().as_secs_f64();
    let mlups = (wt.nx * wt.ny) as f64 * done as f64 / secs / 1e6;
    let ok = wt.engine.is_finite();
    let mut r = analyze(&t, &cx, &cz, &cm);
    if corr && r.cx.is_finite() {
        let (es, ew) = wt.blockage(r.cx);
        r.st /= 1.0 + es + ew;
        if !quiet {
            eprintln!("Blocage de la veine : ε = {:.2} % (volume {:.2} %, sillage {:.2} %), coefficients corrigés (--paroi 0 : bruts)", 100.0 * (es + ew), 100.0 * es, 100.0 * ew);
        }
    }
    Ok((r, mlups, ok))
}

fn cmd_run(a: &Args) {
    let cfg = a.config();
    let wt = WindTunnel::new(cfg.clone()).unwrap_or_else(|e| {
        eprintln!("{e}");
        std::process::exit(1)
    });
    println!(
        "{} · Re = {} · L = {} mailles · grille {}×{} · {} · τ = {:.4} · sous-maille {}",
        wt.body.label(),
        cfg.re,
        cfg.l,
        wt.nx,
        wt.ny,
        wt.engine.name(),
        wt.engine.tau0(),
        if wt.engine.les_active() { "active" } else { "inactive" }
    );
    if wt.sa_active() {
        let profiled = wt.panel.is_some();
        println!(
            "Spalart–Allmaras + loi de paroi · y⁺ ≈ {:.0} au premier nœud · transition {}",
            wt.first_cell_yplus(),
            match (cfg.transition, profiled) {
                (Transition::Turbulent, _) => "au bord d'attaque".to_string(),
                (Transition::Free, true) => {
                    let r = wt.panel.as_ref().unwrap();
                    format!(
                        "libre (extrados {}, intrados {})",
                        pct(r.transition_upper),
                        pct(r.transition_lower)
                    )
                }
                (Transition::Free, false) => "non imposée (couche limite turbulente)".to_string(),
            }
        );
    }
    if let Some(r) = &wt.panel {
        let k = section_coefficients(r);
        println!(
            "Panneaux + couche limite : Cx = {:.4}  Cz = {:.4}  Cm = {:.4}{}",
            k.cd,
            k.cl,
            k.cm,
            if k.kirchhoff {
                format!(
                    "  (décollement à {} de corde : modèle de Kirchhoff)",
                    pct(r.separation_upper.min(r.separation_lower))
                )
            } else {
                String::new()
            }
        );
    }
    drop(wt);
    let t_end = a.get("t", 40.0f64);
    let avg = a.get("avg", 15.0f64);
    match simulate(cfg, t_end, avg, false, a.get("paroi", 1usize) == 1) {
        Ok((r, mlups, ok)) => {
            println!(
                "Cx = {:.4}  Cz = {:.4}  Cm = {:.4}{}{}",
                r.cx,
                r.cz,
                r.cm,
                if r.periodic {
                    format!("  St = {:.4} ({} périodes, Cz' = ±{:.3})", r.st, r.periods, r.cz_amp)
                } else {
                    String::new()
                },
                if ok { String::new() } else { "  [DIVERGENCE]".into() }
            );
            println!("{mlups:.0} millions de mailles par seconde");
        }
        Err(e) => eprintln!("{e}"),
    }
}

/// Plage d'incidences a:b:pas (ou une seule valeur).
fn alpha_range(a: &Args, default: &str) -> Vec<f64> {
    let spec = a.get("alpha", default.to_string());
    let p: Vec<f64> = spec.split(':').filter_map(|s| s.parse().ok()).collect();
    match p.len() {
        1 => return p,
        3 if p[2] > 0.0 => {}
        _ => usage(),
    }
    let mut alphas = Vec::new();
    let mut k = 0;
    loop {
        let x = p[0] + k as f64 * p[2];
        if x > p[1] + 1e-9 {
            break;
        }
        alphas.push((x * 1e6).round() / 1e6);
        k += 1;
    }
    alphas
}

/// Méthode des panneaux + couche limite intégrale pour la maquette de `cfg`.
fn panel_of(cfg: &TunnelConfig, alpha: f64) -> Option<PanelResult> {
    let shape = PanelShape::from_shape(cfg.shape, &cfg.params)?;
    let opts = match cfg.transition {
        Transition::Free => PanelOptions::free(cfg.re),
        Transition::Turbulent => PanelOptions::tripped(cfg.re, 0.0),
    };
    let r = panel_solve(&shape, alpha, &opts);
    (r.cl.is_finite() && r.cd.is_finite()).then_some(r)
}

fn pct(v: f64) -> String {
    if v.is_finite() {
        format!("{:.0} %", v.clamp(0.0, 1.0) * 100.0)
    } else {
        "—".into()
    }
}

fn cmd_panel(a: &Args) {
    let mut cfg = a.config();
    if !a.opts.contains_key("shape") {
        cfg.shape = ShapeKind::Naca;
    }
    if PanelShape::from_shape(cfg.shape, &cfg.params).is_none() {
        eprintln!("La méthode des panneaux traite les profils NACA, les ellipses et les plaques.");
        std::process::exit(2);
    }
    if !a.opts.contains_key("re") {
        cfg.re = 1e6;
    }
    let mach = a.get("mach", 0.0f64);
    let alphas = alpha_range(a, "-4:16:1");
    eprintln!(
        "{} · Re = {} · transition {}{}",
        match cfg.shape {
            ShapeKind::Naca => format!("NACA {}", cfg.params.naca.code()),
            k => k.label().to_string(),
        },
        cfg.re,
        if cfg.transition == Transition::Free {
            "libre (eᴺ, N = 9)"
        } else {
            "au bord d'attaque"
        },
        if mach > 0.0 {
            format!(" · Mach {mach}")
        } else {
            String::new()
        }
    );
    println!("alpha_deg,Cz,Cx,Cm,Cz_attache,x_tr_extrados,x_tr_intrados,x_decollement,kirchhoff,Cp_min,critique");
    let rows: Vec<_> = alphas.par_iter().map(|&al| (al, panel_of(&cfg, al))).collect();
    for (al, r) in rows {
        let Some(r) = r else {
            println!("{al},,,,,,,,,,");
            continue;
        };
        let k = section_coefficients(&r);
        let cp_min = r.cp.iter().cloned().fold(f64::INFINITY, f64::min);
        let (mut cl, mut cm, mut crit) = (k.cl, k.cm, false);
        let mut cp_min_c = cp_min;
        if mach > 0.0 {
            let c = compressible_correction(mach, cp_min);
            cl *= c.factor;
            cm *= c.factor;
            cp_min_c = c.cp_min;
            crit = c.critical;
        }
        let f = |v: f64| {
            if v.is_finite() {
                format!("{v:.4}")
            } else {
                String::new()
            }
        };
        println!(
            "{al},{cl:.4},{:.5},{cm:.4},{:.4},{},{},{},{},{cp_min_c:.3},{}",
            k.cd,
            r.cl,
            f(r.transition_upper),
            f(r.transition_lower),
            f(r.separation_upper.min(r.separation_lower)),
            k.kirchhoff as u8,
            crit as u8
        );
    }
}

/// Géométrie de l'aile finie lue sur la ligne de commande.
fn wing_params(a: &Args) -> WingParams {
    let d = WingParams::default();
    WingParams {
        ar: a.get("ar", d.ar),
        taper: a.get("taper", d.taper),
        sweep: a.get("sweep", d.sweep),
        dihedral: a.get("dihedral", d.dihedral),
        twist: a.get("twist", d.twist),
        n: a.get("strips", d.n),
    }
}

fn cmd_wing(a: &Args) {
    let p = wing_params(a);
    if !(p.ar >= 1.0 && p.ar <= 40.0 && p.taper > 0.0 && p.taper <= 1.0 && p.n >= 2) {
        eprintln!("Géométrie invalide : 1 ≤ --ar ≤ 40, 0 < --taper ≤ 1, --strips ≥ 2.");
        std::process::exit(2);
    }
    let naca = Naca::parse(&a.get("naca", "2412".to_string())).unwrap_or_else(|| usage());
    let re = a.get("re", 1e6f64);
    let tripped = a.get("transition", "free".to_string()) == "turb";
    let alphas = alpha_range(a, "-6:20:1");
    let t0 = Instant::now();
    let w = Wing::new(p);
    let (rr, rt) = Sections::reynolds(&w, re);
    let sec = Sections::new(
        &w,
        re,
        PolarTable::naca(naca, rr, tripped),
        rt.map(|r| PolarTable::naca(naca, r, tripped)),
    );
    eprintln!(
        "Aile NACA {} · A = {} · λ = {} · flèche {}° · dièdre {}° · vrillage {}° · Re = {re} (corde moyenne) · transition {}",
        naca.code(),
        p.ar,
        p.taper,
        p.sweep,
        p.dihedral,
        p.twist,
        if tripped { "au bord d'attaque" } else { "libre" }
    );
    let polar = wing_polar(&w, &sec, &alphas);
    if a.opts.contains_key("span") {
        // Répartition en envergure à une incidence.
        let r = polar.last().unwrap();
        println!("eta,corde,Cl,Cl_max,cCl,alpha_eff_deg,Cd");
        for (s, cm) in r.strips.iter().zip(&sec.clmax) {
            println!(
                "{:.4},{:.4},{:.4},{:.4},{:.4},{:.3},{:.5}",
                s.eta, s.c, s.cl, cm, s.ccl, s.alpha, s.cd
            );
        }
        eprintln!(
            "α = {}° : CL = {:.4}, CD = {:.5}, e = {:.3}",
            r.alpha, r.cl, r.cd, r.e
        );
    } else {
        println!("alpha_deg,CL,CD,CDi,CDp,Cm,e,L_D,converge");
        for r in &polar {
            println!(
                "{},{:.4},{:.5},{:.5},{:.5},{:.4},{},{:.2},{}",
                r.alpha,
                r.cl,
                r.cd,
                r.cdi,
                r.cdv,
                r.cm,
                if r.cl.abs() > 0.05 && r.e.is_finite() {
                    format!("{:.3}", r.e)
                } else {
                    String::new()
                },
                r.cl / r.cd,
                r.converged as u8
            );
        }
    }
    eprintln!("{:.2} s", t0.elapsed().as_secs_f64());
}

/// Coefficients de référence 3D : Schiller–Naumann et Clift pour la sphère.
pub fn sphere_cd(re: f64) -> f64 {
    24.0 / re * (1.0 + 0.15 * re.powf(0.687)) + 0.42 / (1.0 + 42500.0 * re.powf(-1.16))
}

#[cfg(feature = "gpu")]
fn cmd_run3d(a: &Args) {
    use lbm::bodies3d::{Body3D, Shape3D};
    use lbm::cpu::Ramp;
    use lbm::gpu3d::Gpu3D;
    const U: f32 = 0.1;
    const RAMP_T: f64 = 4.0;
    let l = a.get("l", 12usize);
    let stl = a.opts.get("stl").cloned();
    let shape = if stl.is_some() {
        Shape3D::Mesh
    } else {
        Shape3D::from_id(&a.get("shape", "sphere".to_string())).unwrap_or_else(|| {
            eprintln!("--shape sphere|cube|cylinder|wing, ou --stl fichier.stl");
            std::process::exit(2)
        })
    };
    let mut re = a.get("re", 100.0f64);
    let t_end = a.get("t", 40.0f64);
    let avg = a.get("avg", 15.0f64);
    let fp16 = a.get("fp32", 0u8) == 0;
    let fluid = lbm::fluid::Fluid::from_id(&a.get("fluide", "air".to_string())).unwrap_or_else(|| {
        eprintln!("--fluide air|eau");
        std::process::exit(2)
    });
    // Échelle physique (maquettes STL) : vitesse (m/s), longueur et surface de référence (m, m²).
    let mut phys: Option<(f64, f64)> = None;
    let mut mesh_info: Option<lbm::stl::MeshInfo> = None;
    let body = match &stl {
        Some(path) => {
            use lbm::stl::{Axis, MeshSetup, Unit};
            let axis = |k: &str, d: Axis| {
                a.opts.get(k).map_or(d, |v| {
                    Axis::parse(v).unwrap_or_else(|| {
                        eprintln!("--{k} : +x, -x, +y, -y, +z ou -z");
                        std::process::exit(2)
                    })
                })
            };
            let setup = MeshSetup {
                flow: axis("ecoulement", Axis::Xp),
                up: axis("haut", Axis::Zp),
                unit: Unit::parse(&a.get("unite", "mm".to_string())).unwrap_or_else(|| {
                    eprintln!("--unite mm|cm|m|in");
                    std::process::exit(2)
                }),
                alpha_deg: a.get("alpha", 0.0f64),
                beta_deg: a.get("beta", 0.0f64),
                ref_length: a.opts.get("lref").map(|_| a.get("lref", 0.0f64)),
                ref_area: a.opts.get("sref").map(|_| a.get("sref", 0.0f64)),
                margin: a.get("marge", 1.0f64),
                ground: a.opts.get("sol").map(|_| a.get("sol", 0.0f64)),
            };
            let tris = lbm::stl::read_stl(std::path::Path::new(path)).unwrap_or_else(|e| {
                eprintln!("{e}");
                std::process::exit(1)
            });
            let (b, info) = Body3D::from_mesh(&tris, l, &setup).unwrap_or_else(|e| {
                eprintln!("{path} : {e}");
                std::process::exit(1)
            });
            let u = setup.unit.label();
            let m = setup.unit.meters();
            let [ex, ey, ez] = info.size_m.map(|v| v / m);
            eprintln!(
                "{path} : {} triangles · écoulement selon {} du fichier, haut selon {} · maquette {ex:.4} × {ey:.4} × {ez:.4} {u} (long., larg., haut.) · référence : longueur {:.4} {u} = {l} mailles, surface {:.4} {u}² ({}) · maître-couple {:.4} {u}² · blocage {:.1} %",
                info.triangles,
                setup.flow.label(),
                setup.up.label(),
                info.ref_length_m / m,
                info.ref_area_m2 / (m * m),
                if setup.ref_area.is_some() { "imposée" } else { "maître-couple" },
                info.frontal_m2 / (m * m),
                100.0 * info.blockage
            );
            if info.blockage > 0.05 {
                eprintln!("Attention : blocage de la veine supérieur à 5 %, coefficients surestimés.");
            }
            let d = info.diag.summary();
            if !d.is_empty() {
                eprintln!(
                    "Maillage : {d}{}",
                    if info.diag.watertight() { "" } else { " ; appartenance votée sur trois rayons (x, y, z)" }
                );
            }
            if let Some(v) = a.opts.get("vitesse").map(|_| a.get("vitesse", 0.0f64)) {
                let (nu_f, _) = fluid.props();
                re = v * info.ref_length_m / nu_f;
                phys = Some((v, info.ref_area_m2));
                let mach = v / fluid.sound_speed();
                eprintln!(
                    "{} à {v} m/s : Re = {re:.4e} sur la longueur de référence{}",
                    fluid.label(),
                    if mach > 0.3 { format!(" · Mach {mach:.2} : compressibilité ignorée") } else { String::new() }
                );
            }
            mesh_info = Some(info);
            b
        }
        None => Body3D::new(shape, l, a.get("alpha", 0.0f64)),
    };
    let nu = U as f64 * l as f64 / re;
    let t0 = Instant::now();
    // Spalart–Allmaras et loi de paroi : dès Re = 5·10⁴ en automatique, comme en 2D.
    let sa = match a.get("turb", "auto".to_string()).as_str() {
        "sa" => true,
        "off" | "none" => false,
        _ => re >= lbm::tunnel::RE_SA_AUTO,
    };
    let depth = a.get("raff", 0usize).min(3);
    let mut sim = match Gpu3D::with_refinement(&body, nu, U, fp16, sa, depth) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Calcul 3D impossible : {e}");
            std::process::exit(1);
        }
    };
    let xtr = a.get("xtr", 0.0f64);
    if sa && xtr > 0.0 {
        sim.set_transition(xtr);
        eprintln!("Couche limite laminaire (Blasius) jusqu'à {:.0} % de la corde locale", 100.0 * xtr);
    }
    eprintln!(
        "{} · Re = {} · L = {l} mailles ({}) · grille {}×{}×{} = {:.2} M mailles · {} mailles pariétales · τ = {:.4} · sous-maille {} · {} · {}",
        shape.label(),
        if re < 1e5 { format!("{re:.0}") } else { format!("{re:.3e}") },
        shape.reference(),
        sim.nx,
        sim.ny,
        sim.nz,
        sim.n as f64 / 1e6,
        sim.wall_count,
        sim.tau0,
        if sim.les { "active" } else { "inactive" },
        if fp16 { "demi-précision" } else { "simple précision" },
        sim.adapter_name()
    );
    if let Some((d, w, s)) = sim.finest() {
        eprintln!(
            "Grille fine ×{s} autour de la maquette : {}×{}×{} = {:.2} M mailles, {w} mailles pariétales, {} mailles par longueur",
            d[0],
            d[1],
            d[2],
            (d[0] * d[1] * d[2]) as f64 / 1e6,
            l * s
        );
    }
    if sa {
        // y⁺ attendu au premier nœud (plaque plane turbulente, nœud à une demi-maille).
        let cf = 0.026 * re.powf(-1.0 / 7.0);
        eprintln!(
            "Spalart–Allmaras + loi de paroi (Spalding) · y⁺ ≈ {:.0} au premier nœud",
            0.5 / (l << depth) as f64 * re * (cf / 2.0).sqrt()
        );
    }
    if let Some(lost) = mesh_info.as_ref().and_then(|i| i.volume_lost(sim.solid_count, l)) {
        if lost > 0.2 {
            eprintln!(
                "Note : {:.0} % du volume de la maquette tient dans des parties plus fines qu'une maille, traitées en parois minces (moins précises) : augmenter --l pour les résoudre.",
                100.0 * lost
            );
        }
    }
    if sim.wall_count == 0 {
        eprintln!("Aucun lien de paroi : maquette invisible à cette résolution.");
        std::process::exit(1);
    }
    let spt = l as f64 / U as f64;
    let ramp = Ramp {
        steps: (RAMP_T * spt).round() as u64,
        u: U,
    };
    let q = 0.5 * (U as f64).powi(2) * body.area;
    let chunk = (0.25 * spt).round().max(1.0) as usize;
    let total = (t_end * spt).round() as u64;
    let (mut ts, mut cx, mut cz, mut cy) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let t_run = Instant::now();
    let mut next_print = 5.0;
    while sim.time < total {
        sim.run(chunk, &ramp);
        let f = sim.take_forces();
        if sim.lost {
            eprintln!("Le GPU ne répond plus à t* = {:.1} (pilote réinitialisé ?)", sim.time as f64 / spt);
            std::process::exit(1);
        }
        if !sim.is_finite(&f) {
            eprintln!("Divergence à t* = {:.1}", sim.time as f64 / spt);
            std::process::exit(1);
        }
        let t = sim.time as f64 / spt;
        ts.push(t);
        cx.push(f.fx / q);
        cz.push(f.fz / q);
        cy.push(f.fy / q);
        if t >= next_print {
            eprintln!("  t* = {t:5.1}  Cx = {:8.4}  Cz = {:8.4}  Cy = {:8.4}", f.fx / q, f.fz / q, f.fy / q);
            next_print += 5.0;
        }
    }
    // Répartition de la portance en envergure : efforts par maille pariétale
    // moyennés sur deux temps convectifs de plus, regroupés par demi-maille en y.
    // Moment de tangage (efforts par maille pariétale sur deux temps convectifs
    // de plus), autour du point situé à xref longueur de référence de l'avant
    // de la maquette (quart de corde par défaut), au centre en y et z ;
    // positif à cabrer.
    let wf = if a.opts.contains_key("envergure") || a.opts.contains_key("xref") || shape == Shape3D::Mesh {
        sim.run((2.0 * spt) as usize, &ramp);
        Some(sim.take_wall_forces())
    } else {
        None
    };
    if let Some(w) = &wf {
        let bb = body.bounds;
        let xr = bb.x0 + a.get("xref", 0.25f64) * l as f64;
        let zr = 0.5 * (bb.z0 + bb.z1);
        let my: f64 = w
            .pos
            .iter()
            .zip(&w.force)
            .map(|(p, f)| (p[2] - zr) * f[0] - (p[0] - xr) * f[2])
            .sum();
        eprintln!(
            "Cm = {:.4} (autour de {:.2} L de l'avant, rapporté à la surface et à la longueur de référence)",
            my / (q * l as f64),
            a.get("xref", 0.25f64)
        );
    }
    if let (Some(path), Some(w)) = (a.opts.get("envergure"), &wf) {
        let bb = body.bounds;
        let (yc, half) = (0.5 * (bb.y0 + bb.y1), 0.5 * (bb.y1 - bb.y0));
        let dy = 0.5;
        let nb = ((bb.y1 - bb.y0) / dy).ceil() as usize + 2;
        let mut fz = vec![0.0f64; nb];
        let mut fx = vec![0.0f64; nb];
        for (p, f) in w.pos.iter().zip(&w.force) {
            let k = (((p[1] - bb.y0) / dy).floor().max(0.0) as usize).min(nb - 1);
            fz[k] += f[2];
            fx[k] += f[0];
        }
        // c·cl/c_ref = (dL/dy) / (q c_ref), c_ref = L mailles.
        let qc = 0.5 * (U as f64).powi(2) * l as f64;
        let mut csv = String::from("eta,ccl_cref,ccd_cref\n");
        for k in 0..nb {
            let y = bb.y0 + (k as f64 + 0.5) * dy;
            if fz[k] != 0.0 || fx[k] != 0.0 {
                csv += &format!("{:.4},{:.5},{:.5}\n", (y - yc) / half, fz[k] / dy / qc, fx[k] / dy / qc);
            }
        }
        match std::fs::write(path, csv) {
            Ok(()) => eprintln!("Répartition en envergure écrite dans {path}"),
            Err(e) => eprintln!("Écriture de {path} impossible : {e}"),
        }
    }
    let secs = t_run.elapsed().as_secs_f64();
    let from = ts.iter().position(|&t| t >= t_end - avg).unwrap_or(0);
    let r = analyze(&ts[from..], &cx[from..], &cz[from..], &cy[from..]);
    println!(
        "Cx = {:.4}  Cz = {:.4}  Cy = {:.4}{}{}",
        r.cx,
        r.cz,
        r.cm,
        if r.periodic {
            format!("  St = {:.4}  Cz' = ±{:.4}", r.st, r.cz_amp)
        } else {
            String::new()
        },
        match phys {
            Some((v, area)) => {
                let (_, rho) = fluid.props();
                let qs = 0.5 * rho * v * v * area;
                format!(
                    "\nEfforts à {v} m/s : traînée {:.3e} N, portance {:.3e} N, latéral {:.3e} N",
                    r.cx * qs,
                    r.cz * qs,
                    r.cm * qs
                )
            }
            None => String::new(),
        }
    );
    if shape == Shape3D::Sphere {
        println!("Référence sphère (Schiller–Naumann, Clift) : Cx ≈ {:.3}", sphere_cd(re));
    }
    // Coupe longitudinale de la vorticité ωy au centre, en image PPM.
    if let Some(path) = a.opts.get("coupe") {
        use lbm::gpu3d::SliceField;
        let s = sim.slice(1, body.center[1].round() as usize, SliceField::Vorticity, l as f64);
        let mut img = format!("P6 {} {} 255\n", s.w, s.h).into_bytes();
        for j in (0..s.h).rev() {
            for i in 0..s.w {
                let v = s.data[i + s.w * j];
                let c = if v > 1e29 {
                    [0x24, 0x31, 0x3a]
                } else {
                    let t = (v / 4.0).clamp(-1.0, 1.0);
                    let k = |a: f32, b: f32| (a + (b - a) * t.abs()) as u8;
                    if t >= 0.0 {
                        [k(243.0, 140.0), k(246.0, 36.0), k(247.0, 16.0)]
                    } else {
                        [k(243.0, 18.0), k(246.0, 63.0), k(247.0, 102.0)]
                    }
                };
                img.extend_from_slice(&c);
            }
        }
        match std::fs::write(path, img) {
            Ok(()) => eprintln!("Coupe ωy écrite dans {path}"),
            Err(e) => eprintln!("Écriture de {path} impossible : {e}"),
        }
    }
    // Vue 3D (maquette, coupe de vorticité, isosurfaces Q) en image PPM.
    if let Some(path) = a.opts.get("image") {
        use lbm::gpu3d::{Camera, RenderOptions};
        let (w, h) = (960usize, 540usize);
        let cam = Camera::orbit(
            body.center,
            a.get("az", 62.0f64),
            a.get("el", 24.0f64),
            a.get("dist", 6.0f64) * l as f64,
        );
        let rgba = sim.render(&RenderOptions {
            camera: cam,
            width: w,
            height: h,
            plane_axis: 1,
            plane_pos: body.center[1],
            vorticity: true,
            iso: true,
            iso_level: a.get("q", 0.4f64),
            l: l as f64,
            background: [0.9, 0.92, 0.93],
            body: [0.2, 0.25, 0.3],
            body_view: match a.get("maquette", "lisse".to_string()).as_str() {
                "grille" => lbm::gpu3d::BodyView::Grid,
                "mailles" => lbm::gpu3d::BodyView::Cells,
                _ => lbm::gpu3d::BodyView::Smooth,
            },
        });
        let mut img = format!("P6 {w} {h} 255\n").into_bytes();
        for px in rgba.chunks_exact(4) {
            img.extend_from_slice(&px[..3]);
        }
        match std::fs::write(path, img) {
            Ok(()) => eprintln!("Vue 3D écrite dans {path}"),
            Err(e) => eprintln!("Écriture de {path} impossible : {e}"),
        }
    }
    println!(
        "{:.0} millions de mailles par seconde · {:.1} s (dont préparation {:.1} s)",
        sim.n as f64 * sim.time as f64 / secs / 1e6,
        t0.elapsed().as_secs_f64(),
        t0.elapsed().as_secs_f64() - secs
    );
}

fn cmd_polar(a: &Args) {
    let alphas = alpha_range(a, "-4:12:2");
    let base = a.config();
    let t_end = a.get("t", 40.0f64);
    let avg = a.get("avg", 15.0f64);
    let gpu = base.engine == EngineKind::Gpu;
    eprintln!(
        "{} incidences, {}",
        alphas.len(),
        if gpu {
            "GPU, l'une après l'autre"
        } else {
            "en parallèle sur les cœurs"
        }
    );
    let t0 = Instant::now();
    let corr = a.get("paroi", 1usize) == 1;
    let run = |alpha: f64| {
        let mut cfg = base.clone();
        cfg.params.alpha_deg = alpha;
        if !gpu {
            // Un fil par incidence : les incidences se partagent les cœurs.
            cfg.engine = EngineKind::Cpu { threads: 1 };
        }
        (alpha, simulate(cfg, t_end, avg, true, corr))
    };
    let res: Vec<_> = if gpu {
        alphas.iter().map(|&a| run(a)).collect()
    } else {
        alphas.par_iter().map(|&a| run(a)).collect()
    };
    println!("alpha_deg,Re,L,Cx,Cz,Cm,St,Cz_Cx,Cx_panneaux,Cz_panneaux");
    for (alpha, r) in res {
        let pan = match panel_of(&base, alpha).filter(|_| base.re >= lbm::tunnel::RE_PANEL_MIN) {
            Some(p) => {
                let k = section_coefficients(&p);
                format!(",{:.5},{:.4}", k.cd, k.cl)
            }
            None => ",,".into(),
        };
        match r {
            Ok((r, _, ok)) => println!(
                "{alpha},{},{},{:.5},{:.5},{:.5},{},{}{pan}",
                base.re,
                base.l,
                r.cx,
                r.cz,
                r.cm,
                if r.periodic {
                    format!("{:.4}", r.st)
                } else {
                    String::new()
                },
                if ok && r.cx.abs() > 1e-9 {
                    format!("{:.3}", r.cz / r.cx)
                } else {
                    "divergence".into()
                }
            ),
            Err(e) => println!("{alpha},,,,,,,{e}{pan}"),
        }
    }
    eprintln!("{:.1} s", t0.elapsed().as_secs_f64());
}

fn cmd_bench(a: &Args) {
    let l = a.get("l", 32usize);
    let steps = a.get("steps", 400usize);
    let base = TunnelConfig {
        l,
        re: 1000.0,
        shape: ShapeKind::Naca,
        params: ShapeParams {
            alpha_deg: 4.0,
            ..Default::default()
        },
        ..Default::default()
    };
    let cores = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    let mut cases = vec![
        ("CPU, 1 fil".to_string(), EngineKind::Cpu { threads: 1 }, false),
        (format!("CPU, {cores} fils"), EngineKind::Cpu { threads: cores }, false),
    ];
    if a.get("gpu", 1u8) == 1 {
        cases.push(("GPU".to_string(), EngineKind::Gpu, false));
        cases.push(("GPU demi".to_string(), EngineKind::Gpu, true));
    }
    println!(
        "NACA 0012, α = 4°, Re = 1000, L = {l} : grille {}×{} = {:.2} M mailles",
        24 * l,
        12 * l,
        (24 * l * 12 * l) as f64 / 1e6
    );
    for (name, engine, gpu_half) in cases {
        let mut wt = match WindTunnel::new(TunnelConfig {
            engine,
            gpu_half,
            ..base.clone()
        }) {
            Ok(w) => w,
            Err(e) => {
                println!("{name:12} : indisponible ({e})");
                continue;
            }
        };
        wt.run(20);
        let _ = wt.engine.take_forces();
        let t0 = Instant::now();
        wt.run(steps);
        let _ = wt.engine.take_forces(); // attend la fin du calcul GPU
        let secs = t0.elapsed().as_secs_f64();
        let mlups = (wt.nx * wt.ny) as f64 * steps as f64 / secs / 1e6;
        let tps = steps as f64 / secs / wt.steps_per_tstar();
        println!(
            "{name:12} : {mlups:8.1} MLUPS · {:7.1} pas/s · {tps:6.2} t*/s  [{}]",
            steps as f64 / secs,
            wt.engine.name()
        );
    }
}

fn main() {
    let a = Args::parse();
    match a.cmd.as_str() {
        "run" => cmd_run(&a),
        "polar" => cmd_polar(&a),
        "panel" => cmd_panel(&a),
        "bench" => cmd_bench(&a),
        "wing" => cmd_wing(&a),
        #[cfg(feature = "gpu")]
        "run3d" => cmd_run3d(&a),
        "" => {
            #[cfg(feature = "gui")]
            launch::gui(&a);
            #[cfg(not(feature = "gui"))]
            usage();
        }
        "help" | "-h" | "--help" => usage(),
        _ => usage(),
    }
}
