//! Format de géométrie de la soufflerie (`.sgeo`) : maillage indexé en double
//! précision, avec sa mise en place (unité, axes, références, attitude) et
//! ses arêtes vives, pour que la ligne de commande et l'interface n'aient ni
//! à les redemander ni à les deviner.
//!
//! Fichier texte UTF-8, un mot-clé par ligne, `#` pour un commentaire :
//!
//! ```text
//! soufflerie-geometrie 1
//! nom TR-1208
//! unite in
//! ecoulement +x
//! haut +z
//! lref 16.672
//! sref 2024.35
//! alpha 4
//! beta 0
//! sol 1.5
//! sommets 3
//! 0 0 0
//! …
//! faces 1
//! 0 1 2
//! aretes 1
//! 0 1
//! fin
//! ```
//!
//! Les faces sont orientées vers l'extérieur ; `aretes` liste les arêtes vives
//! (paires de sommets), calculées à l'export d'après l'angle dièdre et
//! modifiables à la main (ajouter un bord de fuite émoussé, retirer une arête).
//! Toutes les clés de mise en place sont facultatives ; la ligne de commande
//! les remplace quand elle les donne.

use std::fmt::Write as _;
use std::path::Path;

use crate::geom3d::{feature_edges, FeatureEdge, Indexed, FEATURE_DEG};
use crate::stl::{Axis, Tri, Unit};

pub const MAGIC: &str = "soufflerie-geometrie";
pub const VERSION: u32 = 1;

/// Contenu d'un fichier `.sgeo`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GeoFile {
    pub name: Option<String>,
    pub unit: Option<Unit>,
    pub flow: Option<Axis>,
    pub up: Option<Axis>,
    pub ref_length: Option<f64>,
    pub ref_area: Option<f64>,
    pub alpha_deg: Option<f64>,
    pub beta_deg: Option<f64>,
    pub ground: Option<f64>,
    pub mesh: Indexed,
    /// Arêtes vives (indices de sommets).
    pub edges: Vec<[u32; 2]>,
}

impl GeoFile {
    /// Depuis des triangles (STL) : sommets soudés, arêtes vives d'après
    /// l'angle dièdre (`FEATURE_DEG`).
    pub fn from_tris(tris: &[Tri]) -> GeoFile {
        let mesh = Indexed::from_tris(tris);
        let edges = feature_indices(&mesh, &feature_edges(&mesh, FEATURE_DEG));
        GeoFile {
            mesh,
            edges,
            ..Default::default()
        }
    }

    pub fn tris(&self) -> Vec<Tri> {
        self.mesh.tris()
    }

    /// Arêtes vives en segments (coordonnées du fichier), sans angle ni sens
    /// connus : marquées saillantes à 180°.
    pub fn feature_segments(&self) -> Vec<FeatureEdge> {
        self.edges
            .iter()
            .map(|&[a, b]| FeatureEdge {
                a: self.mesh.verts[a as usize],
                b: self.mesh.verts[b as usize],
                angle_deg: 180.0,
                convex: true,
            })
            .collect()
    }

    pub fn to_text(&self) -> String {
        let mut s = format!("{MAGIC} {VERSION}\n");
        s += "# Soufflerie de poche : maillage indexé, faces orientées vers l'extérieur.\n";
        if let Some(n) = &self.name {
            let _ = writeln!(s, "nom {}", n.replace('\n', " "));
        }
        if let Some(u) = self.unit {
            let _ = writeln!(s, "unite {}", u.label());
        }
        let ax = |a: Axis| a.label().replace('−', "-");
        if let Some(a) = self.flow {
            let _ = writeln!(s, "ecoulement {}", ax(a));
        }
        if let Some(a) = self.up {
            let _ = writeln!(s, "haut {}", ax(a));
        }
        for (k, v) in [
            ("lref", self.ref_length),
            ("sref", self.ref_area),
            ("alpha", self.alpha_deg),
            ("beta", self.beta_deg),
            ("sol", self.ground),
        ] {
            if let Some(v) = v {
                let _ = writeln!(s, "{k} {v:?}");
            }
        }
        let _ = writeln!(s, "sommets {}", self.mesh.verts.len());
        for p in &self.mesh.verts {
            // `{:?}` : écriture la plus courte qui relit le même f64.
            let _ = writeln!(s, "{:?} {:?} {:?}", p[0], p[1], p[2]);
        }
        let _ = writeln!(s, "faces {}", self.mesh.faces.len());
        for f in &self.mesh.faces {
            let _ = writeln!(s, "{} {} {}", f[0], f[1], f[2]);
        }
        let _ = writeln!(s, "aretes {}", self.edges.len());
        for e in &self.edges {
            let _ = writeln!(s, "{} {}", e[0], e[1]);
        }
        s += "fin\n";
        s
    }

    pub fn parse(text: &str) -> Result<GeoFile, String> {
        let mut lines = text
            .lines()
            .enumerate()
            .map(|(k, l)| (k + 1, l.split('#').next().unwrap_or("").trim()))
            .filter(|(_, l)| !l.is_empty());
        let err = |k: usize, m: &str| format!("ligne {k} : {m}");
        let (k, head) = lines.next().ok_or("fichier vide")?;
        let mut w = head.split_ascii_whitespace();
        if w.next() != Some(MAGIC) {
            return Err(err(k, &format!("« {MAGIC} » attendu en tête")));
        }
        let v: u32 = w.next().and_then(|v| v.parse().ok()).ok_or_else(|| err(k, "version illisible"))?;
        if v > VERSION {
            return Err(err(k, &format!("version {v} plus récente que celle du programme ({VERSION})")));
        }
        let mut g = GeoFile::default();
        let mut done = false;
        while let Some((k, l)) = lines.next() {
            let (key, rest) = l.split_once(char::is_whitespace).unwrap_or((l, ""));
            let rest = rest.trim();
            let num = |r: &str| r.parse::<f64>().ok().filter(|v| v.is_finite()).ok_or_else(|| err(k, &format!("nombre illisible : {r}")));
            let count = |r: &str| r.parse::<usize>().map_err(|_| err(k, &format!("nombre d'éléments illisible : {r}")));
            match key {
                "nom" => g.name = Some(rest.to_string()),
                "unite" => g.unit = Some(Unit::parse(rest).ok_or_else(|| err(k, "unité : mm, cm, m ou in"))?),
                "ecoulement" => g.flow = Some(Axis::parse(rest).ok_or_else(|| err(k, "axe : +x, -x, +y, -y, +z ou -z"))?),
                "haut" => g.up = Some(Axis::parse(rest).ok_or_else(|| err(k, "axe : +x, -x, +y, -y, +z ou -z"))?),
                "lref" => g.ref_length = Some(num(rest)?),
                "sref" => g.ref_area = Some(num(rest)?),
                "alpha" => g.alpha_deg = Some(num(rest)?),
                "beta" => g.beta_deg = Some(num(rest)?),
                "sol" => g.ground = Some(num(rest)?),
                "sommets" => {
                    let n = count(rest)?;
                    g.mesh.verts.reserve(n);
                    for _ in 0..n {
                        let (k, l) = lines.next().ok_or("fin de fichier dans les sommets")?;
                        let v: Vec<f64> = l.split_ascii_whitespace().filter_map(|t| t.parse().ok()).collect();
                        if v.len() != 3 || v.iter().any(|c| !c.is_finite()) {
                            return Err(err(k, "sommet : trois nombres attendus"));
                        }
                        g.mesh.verts.push([v[0], v[1], v[2]]);
                    }
                }
                "faces" | "aretes" => {
                    let n = count(rest)?;
                    let m = if key == "faces" { 3 } else { 2 };
                    for _ in 0..n {
                        let (k, l) = lines.next().ok_or("fin de fichier dans les faces ou les arêtes")?;
                        let v: Vec<u32> = l.split_ascii_whitespace().filter_map(|t| t.parse().ok()).collect();
                        if v.len() != m || v.iter().any(|&i| i as usize >= g.mesh.verts.len()) {
                            return Err(err(k, &format!("{m} numéros de sommets (déclarés avant) attendus")));
                        }
                        if m == 3 {
                            g.mesh.faces.push([v[0], v[1], v[2]]);
                        } else {
                            g.edges.push([v[0], v[1]]);
                        }
                    }
                }
                "fin" => {
                    done = true;
                    break;
                }
                _ => return Err(err(k, &format!("mot-clé inconnu : {key}"))),
            }
        }
        if !done {
            return Err("« fin » manquant : fichier tronqué ?".into());
        }
        if g.mesh.faces.is_empty() {
            return Err("aucune face".into());
        }
        Ok(g)
    }

    pub fn read(path: &Path) -> Result<GeoFile, String> {
        let t = std::fs::read_to_string(path).map_err(|e| format!("lecture de {} impossible : {e}", path.display()))?;
        GeoFile::parse(&t).map_err(|e| format!("{} : {e}", path.display()))
    }

    pub fn write(&self, path: &Path) -> std::io::Result<()> {
        std::fs::write(path, self.to_text())
    }
}

/// Indices des arêtes vives (sommets soudés de `m`).
fn feature_indices(m: &Indexed, fe: &[FeatureEdge]) -> Vec<[u32; 2]> {
    let mut at = std::collections::HashMap::new();
    for (k, p) in m.verts.iter().enumerate() {
        at.insert(p.map(f64::to_bits), k as u32);
    }
    fe.iter()
        .filter_map(|e| Some([*at.get(&e.a.map(f64::to_bits))?, *at.get(&e.b.map(f64::to_bits))?]))
        .collect()
}

/// Fichier de géométrie d'après l'extension (`.sgeo`), sinon STL.
pub fn is_sgeo(path: &Path) -> bool {
    path.extension().is_some_and(|e| e.eq_ignore_ascii_case("sgeo"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aller_retour() {
        let tris = crate::geom3d::tests::naca_wing(0.12, 10.0, 40.0, 20);
        let mut g = GeoFile::from_tris(&tris);
        g.name = Some("aile d'essai".into());
        g.unit = Some(Unit::Inch);
        g.flow = Some(Axis::Xp);
        g.up = Some(Axis::Zm);
        g.ref_length = Some(16.672);
        g.ref_area = Some(2024.35);
        g.alpha_deg = Some(4.0);
        assert!(!g.edges.is_empty());
        let t = g.to_text();
        let h = GeoFile::parse(&t).unwrap();
        assert_eq!(h, g);
        // Mêmes triangles, au bit près (double précision).
        assert_eq!(h.tris(), Indexed::from_tris(&tris).tris());
        assert!(t.contains("haut -z"));
        // Défauts signalés avec leur ligne.
        assert!(GeoFile::parse("solid x").unwrap_err().contains("ligne 1"));
        let cut = &t[..t.len() - 5];
        assert!(GeoFile::parse(cut).is_err());
        let bad = t.replacen("
faces ", "
fâces ", 1);
        assert!(GeoFile::parse(&bad).unwrap_err().contains("inconnu"));
    }
}
