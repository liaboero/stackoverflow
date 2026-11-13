//! Soufflerie numérique 2D par la méthode de Boltzmann sur réseau (D2Q9).
//!
//! Portage du moteur de l'application web « Soufflerie de poche » : collision
//! régularisée récursive (ou BGK), modèle de sous-maille de Smagorinsky amorti
//! près de la paroi, paroi interpolée de Bouzidi (y compris les parois plus
//! minces qu'une maille), entrée et sortie de Guo, bandes absorbantes, zone
//! éponge, mise en vitesse progressive, efforts par échange de quantité de
//! mouvement.
//!
//! Deux moteurs partagent la géométrie et la veine d'essai :
//! - [`cpu::CpuSolver`] : processeur, tous les cœurs (rayon), boucles
//!   vectorisables sur les mailles fluides ;
//! - `gpu::GpuSolver` (option `gpu`) : wgpu (Vulkan, Metal, DirectX 12).

// Boucles indicées sur les 9 directions : plus lisibles ici que des itérateurs.
#![allow(clippy::needless_range_loop)]
// `!(a > 0.0)` : traite aussi NaN comme invalide.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

pub mod bodies3d;
pub mod cpu;
pub mod d2q9;
pub mod geometry;
#[cfg(feature = "gpu")]
pub mod gpu;
#[cfg(feature = "gpu")]
pub mod gpu3d;
pub mod panel;
pub mod sa;
pub mod shapes;
pub mod stats;
pub mod tunnel;
pub mod wall;
pub mod wing;

pub use shapes::{Naca, ShapeKind, ShapeParams};
pub use tunnel::{Coefficients, EngineKind, Transition, TunnelConfig, Turbulence, WindTunnel};
