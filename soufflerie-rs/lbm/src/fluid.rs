//! Fluides de l'échelle physique.

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Fluid {
    Air,
    Water,
}

impl Fluid {
    pub fn from_id(s: &str) -> Option<Fluid> {
        match s {
            "air" => Some(Fluid::Air),
            "eau" | "water" => Some(Fluid::Water),
            _ => None,
        }
    }
    pub fn label(self) -> &'static str {
        match self {
            Fluid::Air => "Air 15 °C",
            Fluid::Water => "Eau 20 °C",
        }
    }
    /// Viscosité cinématique (m²/s) et masse volumique (kg/m³).
    pub fn props(self) -> (f64, f64) {
        match self {
            Fluid::Air => (1.46e-5, 1.225),
            Fluid::Water => (1.0e-6, 998.0),
        }
    }
    /// Vitesse du son (m/s) ; l'eau est traitée comme incompressible.
    pub fn sound_speed(self) -> f64 {
        match self {
            Fluid::Air => 340.0,
            Fluid::Water => 1480.0,
        }
    }
}
