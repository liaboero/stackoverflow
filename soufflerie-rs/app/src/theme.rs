//! Thème clair ou sombre : couleurs de l'interface et palettes des champs,
//! reprises de la version web (css/app.css, js/render.js).
//!
//! Le thème courant est mémorisé dans un drapeau global, mis à jour à chaque
//! image par l'interface : les fils de calcul (champ coloré, rendu 3D) le
//! lisent sans passer par l'état de l'interface.

use std::sync::atomic::{AtomicBool, Ordering};

use eframe::egui::{self, Color32};

static DARK: AtomicBool = AtomicBool::new(false);

/// Palette de couleurs : (position, RVB).
pub type Stops = [(f32, [u8; 3]); 5];

/// Relève le thème de l'interface (à appeler à chaque image).
pub fn sync(ctx: &egui::Context) {
    DARK.store(ctx.global_style().visuals.dark_mode, Ordering::Relaxed);
}

pub fn dark() -> bool {
    DARK.load(Ordering::Relaxed)
}

const fn rgb(c: u32) -> Color32 {
    Color32::from_rgb((c >> 16) as u8, (c >> 8) as u8, c as u8)
}

fn pick(light: Color32, dark_: Color32) -> Color32 {
    if dark() {
        dark_
    } else {
        light
    }
}

/// Traînée, couleur d'accent.
pub fn accent() -> Color32 {
    pick(rgb(0xc8461c), rgb(0xf07a43))
}
/// Portance.
pub fn lift() -> Color32 {
    pick(rgb(0x1f6f9c), rgb(0x5eb1e4))
}
/// Avertissement.
pub fn warn() -> Color32 {
    pick(rgb(0x93650f), rgb(0xe2b04a))
}
/// Erreur.
pub fn err() -> Color32 {
    pick(rgb(0xb3261e), rgb(0xff8a80))
}
/// Calcul en marche, résultat établi.
pub fn ok() -> Color32 {
    pick(rgb(0x2e7d4f), rgb(0x6cc795))
}
/// Moment de tangage.
pub fn moment() -> Color32 {
    pick(rgb(0x6b4fa0), rgb(0xb59cf0))
}
/// Couleur translucide.
pub fn alpha(c: Color32, a: u8) -> Color32 {
    Color32::from_rgba_unmultiplied(c.r(), c.g(), c.b(), a)
}

/// Style commun : espacements un peu plus aérés, coins arrondis, boutons
/// plus hauts (cibles plus faciles à atteindre).
pub fn setup_style(ctx: &egui::Context) {
    ctx.all_styles_mut(|s| {
        s.spacing.item_spacing = egui::vec2(8.0, 5.0);
        s.spacing.button_padding = egui::vec2(8.0, 3.0);
        s.spacing.interact_size.y = 22.0;
        s.spacing.slider_width = 130.0;
        let r = egui::CornerRadius::same(5);
        s.visuals.widgets.inactive.corner_radius = r;
        s.visuals.widgets.hovered.corner_radius = r;
        s.visuals.widgets.active.corner_radius = r;
        s.visuals.widgets.open.corner_radius = r;
        s.visuals.widgets.noninteractive.corner_radius = r;
    });
}

/// Palettes du rendu 3D (lbm/src/gpu3d.rs, seqMap et divMap), pour les
/// échelles de couleur de la vue 3D.
pub fn render_seq(t: f32) -> Color32 {
    let t = t.clamp(0.0, 1.0);
    let a = [0.05, 0.11, 0.27];
    let b = [0.10, 0.45, 0.62];
    let c = [0.42, 0.78, 0.62];
    let d = [0.99, 0.88, 0.35];
    let mix = |p: [f32; 3], q: [f32; 3], u: f32| [0, 1, 2].map(|k| p[k] + (q[k] - p[k]) * u);
    let v = if t < 0.33 {
        mix(a, b, t / 0.33)
    } else if t < 0.66 {
        mix(b, c, (t - 0.33) / 0.33)
    } else {
        mix(c, d, (t - 0.66) / 0.34)
    };
    Color32::from_rgb((v[0] * 255.0) as u8, (v[1] * 255.0) as u8, (v[2] * 255.0) as u8)
}

/// Palette divergente du rendu 3D, t ∈ [−1, 1] ; valeur nulle : teinte du fond.
pub fn render_div(t: f32) -> Color32 {
    let t = t.clamp(-1.0, 1.0);
    let bg = canvas();
    let mid = [bg.r(), bg.g(), bg.b()].map(|c| (c as f32 / 255.0 * 1.04 + 0.01).clamp(0.0, 1.0));
    let (end, u) = if t < 0.0 { ([0.13, 0.36, 0.70], -t) } else { ([0.78, 0.20, 0.15], t) };
    let v = [0, 1, 2].map(|k| mid[k] + (end[k] - mid[k]) * u);
    Color32::from_rgb((v[0] * 255.0) as u8, (v[1] * 255.0) as u8, (v[2] * 255.0) as u8)
}
/// Texte et traits sur le fond de la veine.
pub fn ink() -> Color32 {
    pick(rgb(0x13202a), rgb(0xe1e8ed))
}
/// Fond de la veine et des vues.
pub fn canvas() -> Color32 {
    pick(rgb(0xe6ebee), rgb(0x0c1216))
}
/// Remplissage et contour des maquettes.
pub fn body() -> Color32 {
    pick(rgb(0x24313a), rgb(0xc2cdd5))
}
pub fn body_edge() -> Color32 {
    pick(rgb(0x0d161c), rgb(0xeef3f6))
}
pub fn body_rgb() -> [u8; 3] {
    let c = body();
    [c.r(), c.g(), c.b()]
}
/// Traceurs (particules, lignes d'émission).
pub fn tracer() -> Color32 {
    pick(
        Color32::from_rgba_premultiplied(16, 28, 38, 170),
        Color32::from_rgba_unmultiplied(236, 243, 247, 150),
    )
}
/// Arêtes discrètes (maillage de l'aile, cadre 3D).
pub fn edge() -> Color32 {
    pick(
        Color32::from_rgba_premultiplied(13, 22, 28, 120),
        Color32::from_rgba_unmultiplied(225, 232, 237, 110),
    )
}

const DIV_LIGHT: Stops = [
    (0.0, [0x12, 0x3f, 0x66]),
    (0.22, [0x4f, 0x86, 0xb8]),
    (0.5, [0xf3, 0xf6, 0xf7]),
    (0.78, [0xde, 0x84, 0x56]),
    (1.0, [0x8c, 0x24, 0x10]),
];
const SEQ_LIGHT: Stops = [
    (0.0, [0xf3, 0xf6, 0xf7]),
    (0.3, [0xbc, 0xd8, 0xdf]),
    (0.6, [0x5f, 0x9f, 0xba]),
    (0.85, [0x24, 0x55, 0x7d]),
    (1.0, [0x10, 0x23, 0x3f]),
];
const DIV_DARK: Stops = [
    (0.0, [0x9f, 0xd3, 0xfb]),
    (0.22, [0x2f, 0x74, 0xa8]),
    (0.5, [0x0e, 0x16, 0x1b]),
    (0.78, [0xb4, 0x50, 0x2a]),
    (1.0, [0xff, 0xc0, 0x8f]),
];
const SEQ_DARK: Stops = [
    (0.0, [0x0e, 0x16, 0x1b]),
    (0.3, [0x17, 0x3f, 0x52]),
    (0.6, [0x2f, 0x82, 0x99]),
    (0.85, [0x7c, 0xc6, 0xc6]),
    (1.0, [0xe9, 0xf6, 0xd8]),
];

/// Palette divergente (vorticité, pression) du thème courant.
pub fn diverging() -> &'static Stops {
    if dark() {
        &DIV_DARK
    } else {
        &DIV_LIGHT
    }
}

/// Palette séquentielle (vitesse, viscosité turbulente) du thème courant.
pub fn sequential() -> &'static Stops {
    if dark() {
        &SEQ_DARK
    } else {
        &SEQ_LIGHT
    }
}

/// Couleur d'une palette en t ∈ [0, 1].
pub fn color(stops: &Stops, t: f32) -> [u8; 4] {
    let t = t.clamp(0.0, 1.0);
    let mut s = 0;
    while s < stops.len() - 2 && t > stops[s + 1].0 {
        s += 1;
    }
    let (p0, c0) = stops[s];
    let (p1, c1) = stops[s + 1];
    let u = ((t - p0) / (p1 - p0)).clamp(0.0, 1.0);
    let c = |i: usize| (c0[i] as f32 + (c1[i] as f32 - c0[i] as f32) * u).round() as u8;
    [c(0), c(1), c(2), 255]
}

/// Table de 256 couleurs.
pub fn lut(stops: &Stops) -> Vec<[u8; 4]> {
    (0..256).map(|k| color(stops, k as f32 / 255.0)).collect()
}

/// Bouton de choix du thème (clair, sombre, système).
pub fn switch(ui: &mut egui::Ui) {
    egui::widgets::global_theme_preference_buttons(ui);
}
