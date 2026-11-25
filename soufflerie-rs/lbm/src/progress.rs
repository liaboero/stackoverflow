//! Avancement de la mise en place (voxelisation, grilles, raccords), pour
//! l'interface : le fil qui construit la veine installe un récepteur, les
//! étapes lui sont signalées au fil de l'eau. Sans récepteur, rien.

use std::cell::{Cell, RefCell};

/// Où en est la mise en place.
#[derive(Clone, Debug, Default)]
pub struct Progress {
    /// Étape en cours (« voxelisation », « grille 2/4 : tampons »…).
    pub stage: String,
    /// Fraction faite, entre 0 et 1 (au jugé : grilles de poids égaux).
    pub frac: f64,
}

type Sink = Box<dyn Fn(&Progress)>;

thread_local! {
    static SINK: RefCell<Option<Sink>> = const { RefCell::new(None) };
    /// Grille en construction (indice, nombre) ; (0, 0) avant les grilles.
    static GRID: Cell<(usize, usize)> = const { Cell::new((0, 0)) };
}

/// Installe le récepteur du fil courant (remplace le précédent).
pub fn set_sink(f: impl Fn(&Progress) + 'static) {
    SINK.with(|s| *s.borrow_mut() = Some(Box::new(f)));
    GRID.with(|g| g.set((0, 0)));
}

/// Retire le récepteur du fil courant.
pub fn clear_sink() {
    SINK.with(|s| *s.borrow_mut() = None);
}

fn send(stage: String, frac: f64) {
    SINK.with(|s| {
        if let Some(f) = s.borrow().as_ref() {
            f(&Progress { stage, frac: frac.clamp(0.0, 1.0) });
        }
    });
}

/// Part de la voxelisation et des préparatifs avant la première grille.
const PRE: f64 = 0.1;

/// Étape hors des grilles (maquette) : fraction de la part préparatoire.
pub fn prelude(stage: &str, f: f64) {
    send(stage.to_string(), PRE * f);
}

/// Début de la grille `k` (0 = veine) sur `n`.
pub fn grid(k: usize, n: usize) {
    GRID.with(|g| g.set((k, n.max(1))));
    stage("voxelisation et liens", 0.0);
}

/// Étape dans la grille en cours, `f` : part de cette grille déjà faite.
pub fn stage(label: &str, f: f64) {
    let (k, n) = GRID.with(|g| g.get());
    if n == 0 {
        send(label.to_string(), PRE);
        return;
    }
    let name = if k == 0 { "veine".to_string() } else { format!("grille {}/{}", k + 1, n) };
    send(format!("{name} : {label}"), PRE + (1.0 - PRE) * (k as f64 + f) / n as f64);
}
