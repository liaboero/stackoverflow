//! Lancement surveillé de l'interface graphique.
//!
//! Un pilote graphique peut faire planter le processus sans aucun message
//! (accès mémoire invalide à la création de la fenêtre, extension Vulkan
//! défaillante…). L'interface tourne donc dans un processus enfant : s'il
//! s'arrête anormalement, le code d'erreur est expliqué et un autre rendu est
//! essayé. Sous Windows : wgpu sur DirectX 12, puis OpenGL, puis wgpu avec
//! tous les pilotes (Vulkan compris).

use std::process::Command;

use crate::{gui, Args};

/// Code de sortie de l'enfant quand aucun rendu n'a pu être créé.
const NO_RENDERER: i32 = 3;
const CHILD: &str = "SOUFFLERIE_ENFANT";

/// Rendu essayé : nom affiché, valeur de SOUFFLERIE_RENDU, pilotes wgpu (WGPU_BACKEND).
type Attempt = (&'static str, &'static str, Option<&'static str>);

pub fn gui(a: &Args) {
    if std::env::var_os(CHILD).is_some() {
        child(a);
        return;
    }
    let attempts = attempts(a.opts.get("rendu").map(String::as_str));
    let (exe, args) = match std::env::current_exe() {
        Ok(e) => (e, std::env::args().skip(1).collect::<Vec<_>>()),
        Err(_) => {
            child(a);
            return;
        }
    };
    let user_backend = std::env::var_os("WGPU_BACKEND").is_some();
    eprintln!("Soufflerie de poche : ouverture de l'interface…");
    let mut failures = Vec::new();
    for &(name, renderer, backend) in &attempts {
        eprintln!("  essai du rendu {name}");
        let mut cmd = Command::new(&exe);
        cmd.args(&args).env(CHILD, "1").env("SOUFFLERIE_RENDU", renderer);
        if let (Some(b), false) = (backend, user_backend) {
            cmd.env("WGPU_BACKEND", b);
        }
        match cmd.status() {
            Ok(s) if s.success() => return,
            Ok(s) => {
                let why = match s.code() {
                    Some(NO_RENDERER) => "indisponible".to_string(),
                    Some(c) => explain(c),
                    None => "arrêt par un signal".to_string(),
                };
                eprintln!("  rendu {name} : {why}");
                failures.push(format!("{name} : {why}"));
            }
            Err(e) => {
                eprintln!("  impossible de relancer {} : {e}", exe.display());
                child(a);
                return;
            }
        }
    }
    eprintln!(
        "\nL'interface graphique n'a pas pu s'ouvrir ({}).\n\
         Le calcul reste disponible en ligne de commande, par exemple :\n\
         \x20 soufflerie run --shape naca --naca 2412 --alpha 6 --re 1000 --l 16 --t 10\n\
         \x20 soufflerie panel --naca 2412 --re 3e6 --alpha -4:16:2\n\
         Si ces commandes fonctionnent, le problème vient du pilote graphique :\n\
         le mettre à jour, ou essayer --rendu glow, --rendu dx12, --rendu vulkan.",
        failures.join(" ; ")
    );
    // Lancée par double-clic, la console se fermerait aussitôt : on attend.
    eprintln!("\nAppuyez sur Entrée pour fermer.");
    let _ = std::io::stdin().read_line(&mut String::new());
    std::process::exit(1);
}

/// Rendus à essayer, dans l'ordre.
fn attempts(forced: Option<&str>) -> Vec<Attempt> {
    match forced {
        Some("glow") | Some("opengl") => vec![("OpenGL", "glow", None)],
        Some("dx12") => vec![("wgpu / DirectX 12", "wgpu", Some("dx12"))],
        Some("vulkan") => vec![("wgpu / Vulkan", "wgpu", Some("vulkan"))],
        Some("wgpu") => vec![("wgpu", "wgpu", None)],
        _ if cfg!(windows) => vec![
            ("wgpu / DirectX 12", "wgpu", Some("dx12")),
            ("OpenGL", "glow", None),
            ("wgpu (tous pilotes)", "wgpu", None),
        ],
        _ => vec![("wgpu", "wgpu", None), ("OpenGL", "glow", None)],
    }
}

/// Processus enfant : l'interface elle-même.
fn child(a: &Args) {
    let start = gui::Start {
        cfg: a
            .opts
            .keys()
            .any(|k| !matches!(k.as_str(), "field" | "graph" | "rendu" | "tab" | "stl"))
            .then(|| {
                let mut c = a.config();
                if !a.opts.contains_key("l") {
                    c.l = 32;
                }
                c
            }),
        field: a.opts.get("field").cloned(),
        graph: a.opts.get("graph").cloned(),
        tab: a.opts.get("tab").cloned(),
        stl: a.opts.get("stl").map(std::path::PathBuf::from),
    };
    if !gui::main(start) {
        std::process::exit(NO_RENDERER);
    }
}

/// Code de sortie anormal, en clair (codes NTSTATUS de Windows).
fn explain(code: i32) -> String {
    let u = code as u32;
    let what = match u {
        101 => "erreur interne (message ci-dessus)",
        0xC000_0005 => "accès mémoire invalide, le plus souvent dans le pilote graphique",
        0xC000_001D => "instruction non prise en charge par ce processeur",
        0xC000_00FD => "dépassement de pile",
        0xC000_0409 => "arrêt d'urgence du programme",
        0xC000_0135 | 0xC000_0139 => "bibliothèque (DLL) introuvable",
        0xC000_0142 => "échec d'initialisation d'une bibliothèque (DLL)",
        0xC000_0022 => "accès refusé (stratégie de sécurité ?)",
        _ => "arrêt anormal",
    };
    if u > 0xFFFF {
        format!("{what} (code 0x{u:08X})")
    } else {
        format!("{what} (code {code})")
    }
}
