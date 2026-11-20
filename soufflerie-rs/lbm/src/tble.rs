//! Loi de paroi hors équilibre : équations de couche limite à une dimension
//! (TBLE, Balaras, Cabot et Moin 1996 ; Kawai et Larsson 2012 ; Park et Moin
//! 2014), résolues sous chaque nœud pariétal sur une grille 1D étirée entre
//! la paroi et la hauteur de raccord h (`WM_H`) :
//!
//!   ∂u/∂t + F = ∂/∂y [(ν + ν_t) ∂u/∂y],   u(0) = 0,   u(h) = u_h,
//!
//! u le long de la vitesse tangentielle relevée à h, ν_t = κ y u_τ (1 −
//! e^{−y⁺/A})² (κ = 0,41, A = 17 : Kawai et Larsson), u_τ = √|τ_p| tiré du
//! profil courant. F, constant le long de la normale, vient du champ LBM au
//! point relevé : 0 (équilibre), (1/ρ) ∂p/∂s (gradient de pression), ou
//! (1/ρ) ∂p/∂s + u_h ∂u/∂s + v_n ∂u/∂n (gradient et convection). Pas de temps :
//! celui de la grille (implicite, algorithme de Thomas). Le frottement
//! τ_p = ν ∂u/∂y|₀ remplace u_τ de Spalding dans le glissement de la couche
//! pariétale (`WALL3_WGSL`) ; le reste de la passe est inchangé.
//!
//! Ce module porte la version processeur (tests) et le code WGSL que la
//! passe de paroi reçoit quand le modèle est choisi (`Gpu3D::set_paroi_tble`).

/// Intervalles de la grille 1D (N + 1 points, dont la paroi et h).
pub const TB_N: usize = 32;
/// Constante de von Kármán et longueur d'amortissement (Kawai et Larsson 2012).
pub const KAPPA: f64 = 0.41;
pub const A_PLUS: f64 = 17.0;

/// Termes sources pris au point relevé.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ParoiTble {
    /// F = 0 : équilibre (doit redonner une loi logarithmique proche de Spalding).
    Equilibre,
    /// F = (1/ρ) ∂p/∂s.
    Gradient,
    /// F = (1/ρ) ∂p/∂s + u_h ∂u/∂s + v_n ∂u/∂n.
    GradientConvection,
}

impl ParoiTble {
    /// `SOUFFLERIE_TBLE` = `eq`, `gp` ou `gpc` ; absent : loi d'équilibre de
    /// Spalding (défaut).
    pub fn from_env() -> Option<ParoiTble> {
        match std::env::var("SOUFFLERIE_TBLE").as_deref() {
            Ok("eq") => Some(ParoiTble::Equilibre),
            Ok("gp") => Some(ParoiTble::Gradient),
            Ok("gpc") => Some(ParoiTble::GradientConvection),
            _ => None,
        }
    }

    pub fn code(self) -> u32 {
        match self {
            ParoiTble::Equilibre => 1,
            ParoiTble::Gradient => 2,
            ParoiTble::GradientConvection => 3,
        }
    }

    pub fn nom(self) -> &'static str {
        match self {
            ParoiTble::Equilibre => "TBLE équilibre",
            ParoiTble::Gradient => "TBLE gradient de pression",
            ParoiTble::GradientConvection => "TBLE gradient et convection",
        }
    }
}

/// Grille géométrique de N + 1 points de 0 à h, premier point à y₁ (uniforme
/// si h/y₁ ≤ N).
pub fn grille(h: f64, y1: f64) -> Vec<f64> {
    let n = TB_N as f64;
    let mut ys = vec![0.0; TB_N + 1];
    if h / y1 <= n * 1.0001 {
        for (j, y) in ys.iter_mut().enumerate() {
            *y = h * j as f64 / n;
        }
        return ys;
    }
    // y₁ (r^N − 1)/(r − 1) = h.
    let s = |r: f64| y1 * (r.powf(n) - 1.0) / (r - 1.0) - h;
    let (mut lo, mut hi) = (1.0 + 1e-9, 2.0);
    while s(hi) < 0.0 {
        hi *= 1.5;
    }
    for _ in 0..200 {
        let m = 0.5 * (lo + hi);
        if s(m) < 0.0 {
            lo = m;
        } else {
            hi = m;
        }
    }
    let r = 0.5 * (lo + hi);
    for j in 1..=TB_N {
        ys[j] = y1 * (r.powi(j as i32) - 1.0) / (r - 1.0);
    }
    ys[TB_N] = h;
    ys
}

/// Premier point de la grille : y⁺ = 1 pour u_τ = 0,1 U (frottement bien
/// au-delà de celui des couches limites étudiées).
pub fn premier_point(nu: f64, u_in: f64) -> f64 {
    nu / (0.1 * u_in.abs().max(1e-6))
}

/// Viscosité turbulente de la grille 1D.
pub fn nu_t(y: f64, ut: f64, nu: f64) -> f64 {
    let d = 1.0 - (-y * ut / (nu * A_PLUS)).exp();
    KAPPA * y * ut * d * d
}

/// Profil d'équilibre (contrainte constante u_τ²) mis à l'échelle de u_h :
/// état initial d'un nœud.
pub fn profil_equilibre(ys: &[f64], nu: f64, ut: f64, uh: f64) -> Vec<f64> {
    let mut u = vec![0.0; ys.len()];
    let mut g0 = 1.0 / nu;
    let mut acc = 0.0;
    for j in 1..ys.len() {
        let g = 1.0 / (nu + nu_t(ys[j], ut, nu));
        acc += 0.5 * (g0 + g) * (ys[j] - ys[j - 1]);
        u[j] = acc;
        g0 = g;
    }
    let s = if acc > 0.0 { uh / acc } else { 0.0 };
    u.iter_mut().for_each(|v| *v *= s);
    u
}

/// Un pas (Δt = 1, implicite) ; u[0] = 0, u[N] imposé à u_h. Rend τ_p signé.
pub fn pas(u: &mut [f64], ys: &[f64], nu: f64, uh: f64, f: f64) -> f64 {
    let n = TB_N;
    u[0] = 0.0;
    u[n] = uh;
    let us = (nu * u[1] / ys[1]).abs().sqrt();
    let nf: Vec<f64> = (0..n).map(|j| nu + nu_t(0.5 * (ys[j] + ys[j + 1]), us, nu)).collect();
    let (mut cp, mut dp) = (vec![0.0; n + 1], vec![0.0; n + 1]);
    for j in 1..n {
        let dyc = 0.5 * (ys[j + 1] - ys[j - 1]);
        let lw = nf[j - 1] / ((ys[j] - ys[j - 1]) * dyc);
        let up = nf[j] / ((ys[j + 1] - ys[j]) * dyc);
        let mut d = u[j] - f;
        let mut c = -up;
        if j == n - 1 {
            d += up * uh;
            c = 0.0;
        }
        let m = 1.0 + lw + up + lw * cp[j - 1];
        cp[j] = c / m;
        dp[j] = (d + lw * dp[j - 1]) / m;
    }
    u[n - 1] = dp[n - 1];
    for j in (1..n - 1).rev() {
        u[j] = dp[j] - cp[j] * u[j + 1];
    }
    nu * u[1] / ys[1]
}

/// Code WGSL ajouté à la passe de paroi : tampon d'état (N valeurs par nœud
/// pariétal ; la première sert de marque d'initialisation), grille, et un pas
/// de la TBLE (même calcul que `pas`). `ys` : grille de cette grille de
/// calcul (mailles).
pub fn wgsl_fns(mode: ParoiTble, ys: &[f64]) -> String {
    let yl: Vec<String> = ys.iter().map(|y| format!("{:e}", *y as f32)).collect();
    format!(
        "
@group(0) @binding(9) var<storage, read_write> tble: array<f32>;
const TB_N = {n}u;
const TB_MODE = {mode}u;
const TB_A = {a:?};
const TB_YS = array<f32, {n1}>({ys});
fn tb_nut(y: f32, us: f32) -> f32 {{
  let d = 1.0 - exp(-y * us / (P.nu * TB_A));
  return KAPPA * y * us * d * d;
}}
// Un pas implicite de la TBLE du nœud pariétal wk ; rend τ_p signé.
fn tble_step(wk: u32, uh: f32, f: f32, ut0: f32) -> f32 {{
  var yv = TB_YS;
  var u: array<f32, {n1}>;
  let base = wk * TB_N;
  u[0] = 0.0;
  u[TB_N] = uh;
  if (tble[base] == 0.0) {{
    // Première fois : profil d'équilibre de u_τ de Spalding, à l'échelle de u_h.
    var acc = 0.0;
    var g0 = 1.0 / P.nu;
    for (var j = 1u; j <= TB_N; j++) {{
      let g = 1.0 / (P.nu + tb_nut(yv[j], ut0));
      acc += 0.5 * (g0 + g) * (yv[j] - yv[j - 1u]);
      u[j] = acc;
      g0 = g;
    }}
    let s = select(0.0, uh / acc, acc > 0.0);
    for (var j = 1u; j < TB_N; j++) {{ u[j] *= s; }}
    tble[base] = 1.0;
  }} else {{
    for (var j = 1u; j < TB_N; j++) {{ u[j] = tble[base + j]; }}
  }}
  let us = sqrt(abs(P.nu * u[1] / yv[1]));
  var nf: array<f32, {n}>;
  for (var j = 0u; j < TB_N; j++) {{ nf[j] = P.nu + tb_nut(0.5 * (yv[j] + yv[j + 1u]), us); }}
  var cp: array<f32, {n1}>;
  var dp: array<f32, {n1}>;
  cp[0] = 0.0;
  dp[0] = 0.0;
  for (var j = 1u; j < TB_N; j++) {{
    let dyc = 0.5 * (yv[j + 1u] - yv[j - 1u]);
    let lw = nf[j - 1u] / ((yv[j] - yv[j - 1u]) * dyc);
    let up = nf[j] / ((yv[j + 1u] - yv[j]) * dyc);
    var d = u[j] - f;
    var c = -up;
    if (j == TB_N - 1u) {{ d += up * uh; c = 0.0; }}
    let m = 1.0 + lw + up + lw * cp[j - 1u];
    cp[j] = c / m;
    dp[j] = (d + lw * dp[j - 1u]) / m;
  }}
  u[TB_N - 1u] = dp[TB_N - 1u];
  for (var j = TB_N - 2u; j >= 1u; j--) {{ u[j] = dp[j] - cp[j] * u[j + 1u]; }}
  for (var j = 1u; j < TB_N; j++) {{ tble[base + j] = u[j]; }}
  return P.nu * u[1] / yv[1];
}}
",
        n = TB_N,
        n1 = TB_N + 1,
        mode = mode.code(),
        a = A_PLUS as f32,
        ys = yl.join(", ")
    )
}

/// Déclarations ajoutées au début de la loi de paroi (après `var pg = 0.0;`).
pub const WGSL_DECL: &str = "
    // TBLE : terme source au point relevé, et relevé valide.
    var tbF = 0.0;
    var tbOk = false;
";

/// Termes sources au point relevé (après `yw = WM_H;`) : gradient de pression
/// filtré (différence centrée de ρ à ±2 mailles le long de u_h, comme l'essai
/// 536870912), convection u_h ∂u/∂s (mêmes relevés) et v_n ∂u/∂n (±1 maille
/// le long de la normale ; décentrée vers le haut si le point bas touche un
/// nœud solide).
pub const WGSL_SOURCE: &str = "
        tbOk = true;
        if (TB_MODE >= 2u && ut > 1e-12) {
          let te = th / ut;
          let tp = mac_at(ps + 2.0 * te);
          let tm = mac_at(ps - 2.0 * te);
          if (tp.w > 0.0 && tm.w > 0.0) {
            tbF = (tp.w - tm.w) / (12.0 * rh);
            if (TB_MODE >= 3u) {
              let dus = (dot(tp.xyz, te) - dot(tm.xyz, te)) / 4.0;
              let vn = dot(us3, gw.xyz);
              let a1 = mac_at(ps + gw.xyz);
              let a0 = mac_at(ps - gw.xyz);
              var dun = 0.0;
              if (a1.w > 0.0 && a0.w > 0.0) {
                dun = (dot(a1.xyz, te) - dot(a0.xyz, te)) / 2.0;
              } else if (a1.w > 0.0) {
                dun = dot(a1.xyz, te) - ut;
              }
              tbF += ut * dus + vn * dun;
            }
          }
        }
";

/// u_τ de la TBLE à la place de celui de Spalding (après le calcul de `ut2`).
pub const WGSL_UTAU: &str = "
    if (tbOk) {
      let tw = tble_step(wk, ut, tbF, ut2);
      ut2 = sqrt(max(tw, 0.0));
    }
";

#[cfg(test)]
mod tests {
    use super::*;

    fn spalding(ut: f64, y: f64, nu: f64) -> f64 {
        // Loi de Spalding par dichotomie sur u_τ.
        let (k, eb) = (0.41f64, 0.128_734_90f64);
        let mut lo = 1e-12;
        let mut hi = ut;
        for _ in 0..200 {
            let m = 0.5 * (lo + hi);
            let up = ut / m;
            let ku = k * up;
            let yp = up + eb * (ku.exp() - 1.0 - ku - ku * ku / 2.0 - ku * ku * ku / 6.0);
            if yp > y * m / nu {
                lo = m;
            } else {
                hi = m;
            }
        }
        0.5 * (lo + hi)
    }

    /// Régime établi de la TBLE d'équilibre contre Spalding, de y⁺(h) ≈ 30 à
    /// 4·10⁴ (aile infinie : 2000 à 5000). Écart borné dans la couche
    /// logarithmique (y⁺(h) > 300) ; plus bas, la zone tampon des deux lois
    /// diffère (A = 17 contre la forme de Spalding).
    #[test]
    fn equilibre_contre_spalding() {
        let (h, uin) = (2.0, 0.1);
        for re_h in [5e2, 5e3, 5e4, 2e5, 1e6] {
            let nu = uin * h / re_h;
            let uh = 0.9 * uin;
            let ys = grille(h, premier_point(nu, uin));
            let us = spalding(uh, h, nu);
            let mut u = profil_equilibre(&ys, nu, us, uh);
            let mut tw = 0.0;
            for _ in 0..200_000 {
                tw = pas(&mut u, &ys, nu, uh, 0.0);
            }
            let ut = tw.sqrt();
            let e = ut / us - 1.0;
            println!("y⁺(h) {:.0} : u_τ TBLE / Spalding = {:.4} (y₁⁺ {:.2})", h * us / nu, 1.0 + e, ys[1] * us / nu);
            if h * us / nu > 300.0 {
                assert!(e.abs() < 0.03, "écart {e}");
            }
        }
    }

    /// Gradient adverse : frottement réduit, puis négatif (retour près de la
    /// paroi) quand F h dépasse l'ordre de u_τ².
    #[test]
    fn gradient_adverse() {
        let (h, uin) = (2.0, 0.1);
        let nu = uin * 128.0 / 6e6;
        let uh = 0.9 * uin;
        let ys = grille(h, premier_point(nu, uin));
        let us = spalding(uh, h, nu);
        let run = |f: f64| {
            let mut u = profil_equilibre(&ys, nu, us, uh);
            let mut tw = 0.0;
            for _ in 0..100_000 {
                tw = pas(&mut u, &ys, nu, uh, f);
            }
            tw
        };
        let t0 = run(0.0);
        let t1 = run(0.5 * us * us / h);
        let t2 = run(20.0 * us * us / h);
        println!("τ_p / u_τ² : {:.3} {:.3} {:.3}", t0 / (us * us), t1 / (us * us), t2 / (us * us));
        assert!(t1 < t0 && t1 > 0.0);
        assert!(t2 < 0.0);
    }
}
