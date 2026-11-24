//! Captures d'écran et images PNG.
//!
//! - `write_png` : image RGBA vers un fichier PNG (vue 3D exportée, capture) ;
//! - capture automatique de la fenêtre : `SOUFFLERIE_CAPTURE=fichier.png`
//!   (et `SOUFFLERIE_CAPTURE_T=secondes`, 8 par défaut) enregistre la fenêtre
//!   au bout de ce délai puis la ferme ; sert aux captures de la
//!   documentation et aux vérifications visuelles ;
//! - F12 : capture de la fenêtre dans le dossier courant ;
//! - essais d'interface : `SOUFFLERIE_TOUCHES="3:espace;8:v;9:2"` presse ces
//!   touches à ces instants (secondes), comme au clavier.

use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use eframe::egui;

/// Écrit une image RGBA (8 bits par canal) en PNG.
pub fn write_png(path: &Path, w: usize, h: usize, rgba: &[u8]) -> Result<(), String> {
    let file = std::fs::File::create(path).map_err(|e| format!("{} : {e}", path.display()))?;
    let mut enc = png::Encoder::new(std::io::BufWriter::new(file), w as u32, h as u32);
    enc.set_color(png::ColorType::Rgba);
    enc.set_depth(png::BitDepth::Eight);
    let mut wr = enc.write_header().map_err(|e| e.to_string())?;
    wr.write_image_data(rgba).map_err(|e| e.to_string())?;
    wr.finish().map_err(|e| e.to_string())
}

/// Capture de la fenêtre : demandée au clavier (F12) ou par l'environnement.
pub struct Capture {
    auto: Option<(PathBuf, f64)>,
    /// Touches à presser : (instant, modificateurs, touche).
    keys: Vec<(f64, egui::Modifiers, egui::Key)>,
    start: Instant,
    /// Capture demandée, en attente de l'image (chemin, fermer ensuite).
    pending: Arc<Mutex<Option<(PathBuf, bool)>>>,
    requested: bool,
    /// Dernier message (fichier écrit ou erreur) et son heure.
    pub message: Arc<Mutex<Option<(String, Instant)>>>,
}

impl Capture {
    pub fn from_env() -> Capture {
        let auto = std::env::var_os("SOUFFLERIE_CAPTURE").map(|p| {
            let t = std::env::var("SOUFFLERIE_CAPTURE_T").ok().and_then(|v| v.parse().ok()).unwrap_or(8.0);
            (PathBuf::from(p), t)
        });
        let keys = std::env::var("SOUFFLERIE_TOUCHES")
            .map(|v| {
                v.split(';')
                    .filter_map(|e| {
                        let (t, k) = e.split_once(':')?;
                        let (m, k) = match k.strip_prefix("ctrl+") {
                            Some(k) => (egui::Modifiers::COMMAND, k),
                            None => (egui::Modifiers::NONE, k),
                        };
                        let key = match k {
                            "espace" => egui::Key::Space,
                            "entree" => egui::Key::Enter,
                            k => egui::Key::from_name(&k.to_uppercase()).or_else(|| egui::Key::from_name(k))?,
                        };
                        Some((t.trim().parse().ok()?, m, key))
                    })
                    .collect()
            })
            .unwrap_or_default();
        Capture {
            auto,
            keys,
            start: Instant::now(),
            pending: Arc::new(Mutex::new(None)),
            requested: false,
            message: Arc::new(Mutex::new(None)),
        }
    }

    /// À appeler à chaque image.
    pub fn tick(&mut self, ctx: &egui::Context) {
        let el = self.start.elapsed().as_secs_f64();
        while let Some(&(t, m, key)) = self.keys.first() {
            if t > el {
                ctx.request_repaint_after(std::time::Duration::from_secs_f64((t - el).min(0.5)));
                break;
            }
            self.keys.remove(0);
            eprintln!("essai : touche {key:?} à {el:.1} s");
            ctx.input_mut(|i| {
                i.events.push(egui::Event::Key {
                    key,
                    physical_key: None,
                    pressed: true,
                    repeat: false,
                    modifiers: m,
                })
            });
        }
        if let Some((path, t)) = &self.auto {
            if !self.requested {
                let el = self.start.elapsed().as_secs_f64();
                if el >= *t {
                    self.request(ctx, path.clone(), true);
                } else {
                    ctx.request_repaint_after(std::time::Duration::from_secs_f64((*t - el).min(0.5)));
                }
            }
        }
        if ctx.input(|i| i.key_pressed(egui::Key::F12)) {
            let name = format!(
                "soufflerie-{}.png",
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .map_or(0, |d| d.as_secs())
            );
            self.request(ctx, PathBuf::from(name), false);
        }
    }

    fn request(&mut self, ctx: &egui::Context, path: PathBuf, close: bool) {
        self.requested = true;
        *self.pending.lock().unwrap() = Some((path, close));
        let pending = self.pending.clone();
        let message = self.message.clone();
        let c = ctx.clone();
        ctx.request_screenshot(move |img| {
            let Some((path, close)) = pending.lock().unwrap().take() else {
                return;
            };
            let [w, h] = img.size;
            let rgba: Vec<u8> = img.pixels.iter().flat_map(|p| p.to_array()).collect();
            let msg = match write_png(&path, w, h, &rgba) {
                Ok(()) => format!("Capture enregistrée : {}", path.display()),
                Err(e) => format!("Capture impossible : {e}"),
            };
            eprintln!("{msg}");
            *message.lock().unwrap() = Some((msg, Instant::now()));
            if close {
                c.send_viewport_cmd(egui::ViewportCommand::Close);
            }
            c.request_repaint();
        });
        if !close {
            self.requested = false;
        }
    }
}
