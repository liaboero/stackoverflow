/**
 * Méthode des panneaux et couche limite intégrale (profils minces, Re ≈ 10⁴
 * à 10⁷), complément instantané du calcul Lattice Boltzmann.
 *
 *  - Écoulement potentiel : panneaux tourbillons à densité linéaire, condition
 *    de Kutta (Kuethe & Chow, « Foundations of Aerodynamics »).
 *  - Couche limite laminaire : Thwaites (corrélations de Cebeci–Bradshaw).
 *  - Transition : méthode e^N enveloppe de Drela & Giles (1987), N_crit = 9
 *    par défaut, ou au décollement laminaire, ou imposée.
 *  - Couche limite turbulente : méthode d'entraînement de Head, frottement
 *    de Ludwieg–Tillmann ; décollement si H > 2,4.
 *  - Traînée : Squire–Young au bord de fuite.
 *  - Couplage faible : profil épaissi de δ* (figé à 95 % de corde, bord de
 *    fuite refermé sur la ligne moyenne), une itération sous-relaxée.
 *    NACA 0012, Re = 6·10⁶, turbulent : Cz = 1,17 / Cx = 0,0129 à 10°
 *    (NASA TMR : 1,09 / 0,0123), Cx = 0,0082 à 0° (0,0081).
 *
 * Coordonnées en cordes, x de 0 (bord d'attaque) à 1 (bord de fuite). Sans
 * interaction forte, le décrochage n'est pas prédit : au-delà du décollement
 * turbulent, les résultats sont signalés comme hors du domaine de validité.
 */

/** Points du profil, bord de fuite → intrados → bord d'attaque → extrados → bord de fuite. */
export function airfoilPoints(shape, params = {}, n = 80) {
  const lower = [];
  const upper = [];
  let camber = () => 0;
  let half;
  if (shape === 'naca') {
    const s = String(params.naca ?? '0012');
    const m = Number(s[0]) / 100;
    const p = Number(s[1]) / 10;
    const t = Number(s.slice(2)) / 100;
    half = (x) => 5 * t * (0.2969 * Math.sqrt(x) - 0.126 * x - 0.3516 * x * x + 0.2843 * x ** 3 - 0.1036 * x ** 4);
    if (m > 0) camber = (x) => (x < p ? (m / (p * p)) * (2 * p * x - x * x) : (m / (1 - p) ** 2) * (1 - 2 * p + 2 * p * x - x * x));
  } else if (shape === 'ellipse' || shape === 'plate') {
    // Plaque : ellipse de 2 % d'épaisseur (bords arrondis, plaque mince).
    const t = shape === 'plate' ? 0.02 : Math.min(1, Math.max(0.05, params.thickness ?? 0.4));
    half = (x) => t * Math.sqrt(Math.max(0, x * (1 - x)));
  } else {
    return null;
  }
  for (let k = 0; k <= n; k++) {
    const x = 0.5 * (1 - Math.cos((Math.PI * k) / n));
    upper.push([x, camber(x) + half(x)]);
    lower.push([x, camber(x) - half(x)]);
  }
  // Bord de fuite fermé : un seul point de départ et d'arrivée.
  const pts = [];
  for (let k = n; k >= 0; k--) pts.push(lower[k]);
  for (let k = 1; k <= n; k++) pts.push(upper[k]);
  return pts;
}

/** Résolution d'un système linéaire dense par élimination de Gauss (pivot partiel). */
function solve(A, b) {
  const n = b.length;
  for (let c = 0; c < n; c++) {
    let piv = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(A[r][c]) > Math.abs(A[piv][c])) piv = r;
    [A[c], A[piv]] = [A[piv], A[c]];
    [b[c], b[piv]] = [b[piv], b[c]];
    const d = A[c][c];
    for (let r = c + 1; r < n; r++) {
      const f = A[r][c] / d;
      if (f === 0) continue;
      const Ar = A[r];
      const Ac = A[c];
      for (let k = c; k < n; k++) Ar[k] -= f * Ac[k];
      b[r] -= f * b[c];
    }
  }
  const x = new Float64Array(n);
  for (let r = n - 1; r >= 0; r--) {
    let s = b[r];
    for (let k = r + 1; k < n; k++) s -= A[r][k] * x[k];
    x[r] = s / A[r][r];
  }
  return x;
}

/**
 * Écoulement potentiel autour d'un profil (panneaux tourbillons linéaires).
 * @param {Array<[number, number]>} pts contour fermé (voir airfoilPoints)
 * @param {number} alpha incidence (rad)
 * @returns {{xc, yc, theta, len, vt, cp, cl, cm}} vitesse tangentielle vt/V∞ aux points de contrôle
 */
export function potentialFlow(pts, alpha) {
  const N = pts.length - 1;
  const X = pts.map((p) => p[0]);
  const Y = pts.map((p) => p[1]);
  const xc = new Float64Array(N);
  const yc = new Float64Array(N);
  const th = new Float64Array(N);
  const S = new Float64Array(N);
  for (let j = 0; j < N; j++) {
    xc[j] = 0.5 * (X[j] + X[j + 1]);
    yc[j] = 0.5 * (Y[j] + Y[j + 1]);
    S[j] = Math.hypot(X[j + 1] - X[j], Y[j + 1] - Y[j]);
    th[j] = Math.atan2(Y[j + 1] - Y[j], X[j + 1] - X[j]);
  }
  const An = [];
  const At = [];
  for (let i = 0; i < N; i++) {
    const cn1 = new Float64Array(N);
    const cn2 = new Float64Array(N);
    const ct1 = new Float64Array(N);
    const ct2 = new Float64Array(N);
    for (let j = 0; j < N; j++) {
      if (i === j) {
        cn1[j] = -1;
        cn2[j] = 1;
        ct1[j] = Math.PI / 2;
        ct2[j] = Math.PI / 2;
        continue;
      }
      const dx = xc[i] - X[j];
      const dy = yc[i] - Y[j];
      const A = -dx * Math.cos(th[j]) - dy * Math.sin(th[j]);
      const B = dx * dx + dy * dy;
      const C = Math.sin(th[i] - th[j]);
      const D = Math.cos(th[i] - th[j]);
      const E = dx * Math.sin(th[j]) - dy * Math.cos(th[j]);
      const F = Math.log(1 + (S[j] * (S[j] + 2 * A)) / B);
      const G = Math.atan2(E * S[j], B + A * S[j]);
      const P = dx * Math.sin(th[i] - 2 * th[j]) + dy * Math.cos(th[i] - 2 * th[j]);
      const Q = dx * Math.cos(th[i] - 2 * th[j]) - dy * Math.sin(th[i] - 2 * th[j]);
      cn2[j] = D + (0.5 * Q * F) / S[j] - ((A * C + D * E) * G) / S[j];
      cn1[j] = 0.5 * D * F + C * G - cn2[j];
      ct2[j] = C + (0.5 * P * F) / S[j] + ((A * D - C * E) * G) / S[j];
      ct1[j] = 0.5 * C * F - D * G - ct2[j];
    }
    const an = new Float64Array(N + 1);
    const at = new Float64Array(N + 1);
    an[0] = cn1[0];
    at[0] = ct1[0];
    for (let j = 1; j < N; j++) {
      an[j] = cn1[j] + cn2[j - 1];
      at[j] = ct1[j] + ct2[j - 1];
    }
    an[N] = cn2[N - 1];
    at[N] = ct2[N - 1];
    An.push(an);
    At.push(at);
  }
  const M = An.map((r) => Float64Array.from(r));
  const rhs = new Float64Array(N + 1);
  for (let i = 0; i < N; i++) rhs[i] = Math.sin(th[i] - alpha);
  // Condition de Kutta : γ au bord de fuite nul (γ₁ + γ_{N+1} = 0).
  const kutta = new Float64Array(N + 1);
  kutta[0] = 1;
  kutta[N] = 1;
  M.push(kutta);
  rhs[N] = 0;
  const g = solve(M, rhs);
  const vt = new Float64Array(N);
  const cp = new Float64Array(N);
  let cx = 0;
  let cy = 0;
  let cm = 0;
  for (let i = 0; i < N; i++) {
    let v = Math.cos(th[i] - alpha);
    for (let j = 0; j <= N; j++) v += At[i][j] * g[j];
    vt[i] = v;
    cp[i] = 1 - v * v;
    // Force de pression sur le panneau (contour horaire : normale extérieure à gauche, (−sin θ, cos θ)).
    const nxo = -Math.sin(th[i]);
    const nyo = Math.cos(th[i]);
    const fx = -cp[i] * S[i] * nxo;
    const fy = -cp[i] * S[i] * nyo;
    cx += fx;
    cy += fy;
    cm += -(xc[i] - 0.25) * fy + yc[i] * fx;
  }
  const cl = cy * Math.cos(alpha) - cx * Math.sin(alpha);
  return { xc, yc, theta: th, len: S, vt, cp, cl, cm };
}

// ------------------------------------------------------------ couche limite

function thwaitesHL(lam) {
  if (lam >= 0) {
    const l = Math.min(lam, 0.1);
    return { H: 2.61 - 3.75 * l + 5.24 * l * l, ell: 0.22 + 1.57 * l - 1.8 * l * l };
  }
  const l = Math.max(lam, -0.1);
  return { H: 2.088 + 0.0731 / (l + 0.14), ell: 0.22 + 1.402 * l + (0.018 * l) / (l + 0.107) };
}

/** Taux d'amplification de l'enveloppe e^N (Drela & Giles 1987). */
function dNdS(H, theta, reTheta) {
  const hm = H - 1;
  const log10Re0 = (1.415 / hm - 0.489) * Math.tanh(20 / hm - 12.9) + 3.295 / hm + 0.44;
  if (Math.log10(Math.max(reTheta, 1)) < log10Re0) return 0;
  const dNdRe = 0.01 * Math.sqrt((2.4 * H - 3.7 + 2.5 * Math.tanh(1.5 * H - 4.65)) ** 2 + 0.25);
  const ell = (6.54 * H - 14.07) / (H * H);
  const m = (0.058 * ((H - 4) ** 2) / hm - 0.068) / ell;
  return (dNdRe * 0.5 * (m + 1) * ell) / theta;
}

function headH1(H) {
  return H <= 1.6 ? 3.3 + 0.8234 * (H - 1.1) ** -1.287 : 3.3 + 1.5501 * (H - 0.6778) ** -3.064;
}

function headH(H1) {
  if (H1 <= 3.32) return 3;
  return H1 >= 5.3 ? 0.86 * (H1 - 3.3) ** -0.777 + 1.1 : 1.1538 * (H1 - 3.3) ** -0.326 + 0.6778;
}

/**
 * Couche limite le long d'une face, depuis le point d'arrêt.
 * @param {Float64Array} s abscisse curviligne (cordes), croissante
 * @param {Float64Array} ue vitesse extérieure / V∞ (≥ 0)
 * @param {Float64Array} xs abscisse x/c correspondante
 * @param {object} o { re, nCrit, xTrip }
 */
function boundaryLayer(s, ue, xs, o) {
  const n = s.length;
  const theta = new Float64Array(n);
  const H = new Float64Array(n);
  const cf = new Float64Array(n);
  const re = o.re;
  let transitionX = NaN;
  let separationX = NaN;
  let laminar = true;
  let N = 0;
  let integral = 0;
  let h1 = 0;
  for (let k = 0; k < n; k++) {
    const u = Math.max(ue[k], 1e-4);
    if (laminar) {
      if (k > 0) {
        const ds = s[k] - s[k - 1];
        integral += 0.5 * (Math.max(ue[k - 1], 0) ** 5 + u ** 5) * ds;
      }
      let th2;
      if (k === 0 || integral === 0) {
        // Point d'arrêt (Hiemenz) : θ² = 0,075 ν / (due/ds).
        const due = k + 1 < n ? (ue[k + 1] - ue[k]) / Math.max(s[k + 1] - s[k], 1e-9) : 1;
        th2 = 0.075 / (re * Math.max(due, 1e-3));
      } else {
        th2 = (0.45 * integral) / (re * u ** 6);
      }
      theta[k] = Math.sqrt(th2);
      const kk = Math.min(Math.max(k, 1), n - 2);
      const due = (ue[kk + 1] - ue[kk - 1]) / Math.max(s[kk + 1] - s[kk - 1], 1e-9);
      const lam = th2 * due * re;
      const { H: h, ell } = thwaitesHL(lam);
      H[k] = h;
      const reTheta = Math.max(u * theta[k] * re, 1);
      cf[k] = (2 * ell) / reTheta;
      if (k > 0) N += dNdS(h, theta[k], reTheta) * (s[k] - s[k - 1]);
      const forced = Number.isFinite(o.xTrip) && xs[k] >= o.xTrip;
      if (k > 0 && (N >= o.nCrit || lam < -0.09 || forced)) {
        laminar = false;
        transitionX = xs[k];
        H[k] = 1.4;
        h1 = headH1(1.4);
      }
    } else {
      const ds = s[k] - s[k - 1];
      const u0 = Math.max(ue[k - 1], 1e-4);
      const due = (u - u0) / Math.max(ds, 1e-9);
      // Euler explicite en sous-pas.
      let th = theta[k - 1];
      let hh = H[k - 1];
      let H1v = h1;
      const sub = 4;
      for (let q = 0; q < sub; q++) {
        const uq = u0 + ((u - u0) * (q + 0.5)) / sub;
        const reT = Math.max(uq * th * re, 10);
        const cfq = 0.246 * 10 ** (-0.678 * hh) * reT ** -0.268;
        const dth = cfq / 2 - ((hh + 2) * th * due) / uq;
        const dUH1 = uq * 0.0306 * Math.max(H1v - 3, 0.01) ** -0.6169;
        const thN = Math.max(th + (dth * ds) / sub, 1e-7);
        // d(ue θ H1)/ds
        const uthH1 = uq * th * H1v + (dUH1 * ds) / sub;
        const uN = u0 + ((u - u0) * (q + 1)) / sub;
        H1v = uthH1 / (uN * thN);
        th = thN;
        hh = headH(H1v);
      }
      theta[k] = th;
      H[k] = hh;
      h1 = H1v;
      cf[k] = 0.246 * 10 ** (-0.678 * hh) * Math.max(u * th * re, 10) ** -0.268;
      if (hh > 2.4 && !Number.isFinite(separationX)) separationX = xs[k];
    }
  }
  return { theta, H, cf, transitionX, separationX };
}

/**
 * Polaire d'un point : potentiel + couche limite + couplage faible.
 * @param {string} shape 'naca' | 'ellipse' | 'plate'
 * @param {object} o { params, alphaDeg, re, nCrit = 9, xTrip (x/c, NaN = libre), iterations = 3 }
 */
export function panelSolve(shape, o) {
  const base = airfoilPoints(shape, o.params);
  if (!base) return null;
  const alpha = (o.alphaDeg * Math.PI) / 180;
  const nCrit = o.nCrit ?? 9;
  const xTripU = o.xTripUpper ?? o.xTrip ?? NaN;
  const xTripL = o.xTripLower ?? o.xTrip ?? NaN;
  const re = o.re;
  const inviscid = potentialFlow(base, alpha);
  let pts = base;
  let flow = inviscid;
  let bl = null;
  // Une itération de couplage sous-relaxée : au-delà, sans interaction forte,
  // la correction de portance n'est pas plus juste (essais contre NASA TMR).
  const iterations = o.iterations ?? 1;
  let dstar = new Float64Array(base.length);
  for (let it = 0; it <= iterations; it++) {
    bl = surfaceBL(base, flow, re, nCrit, xTripU, xTripL, o);
    if (it === iterations) break;
    // Profil épaissi de δ* (relaxation 0,5) ; bord de fuite ouvert.
    const nb = base.length;
    const next = new Float64Array(nb);
    const rlx = o.relax ?? 0.5;
    for (let k = 0; k < nb; k++) next[k] = (1 - rlx) * dstar[k] + rlx * bl.dstarNode[k];
    dstar = next;
    // Déplacement de chaque nœud selon la normale extérieure (contour horaire :
    // à gauche du sens de parcours, (−t_y, t_x)).
    const disp = base.map(([x, y], k) => {
      const kp = Math.min(nb - 1, k + 1);
      const km = Math.max(0, k - 1);
      const tx = base[kp][0] - base[km][0];
      const ty = base[kp][1] - base[km][1];
      const tl = Math.hypot(tx, ty) || 1;
      return [(-dstar[k] * ty) / tl, (dstar[k] * tx) / tl];
    });
    // Bord de fuite refermé sur la position moyenne des deux faces épaissies :
    // la condition de Kutta reste posée sur un bord de fuite aigu, déplacé
    // comme la ligne moyenne (décambrage visqueux).
    const dTE = [0.5 * (disp[0][0] + disp[nb - 1][0]), 0.5 * (disp[0][1] + disp[nb - 1][1])];
    pts = base.map(([x, y], k) => {
      let [dx, dy] = disp[k];
      if (x > 0.9) {
        const w = ((x - 0.9) / 0.1) ** 2;
        dx += w * (dTE[0] - dx);
        dy += w * (dTE[1] - dy);
      }
      return [x + dx, y + dy];
    });
    flow = potentialFlow(pts, alpha);
  }
  const cd = bl.cd;
  return {
    cl: flow.cl,
    cd,
    cm: flow.cm,
    clInviscid: inviscid.cl,
    transitionUpper: bl.upper.transitionX,
    transitionLower: bl.lower.transitionX,
    separationUpper: bl.upper.separationX,
    separationLower: bl.lower.separationX,
    separated:
      (Number.isFinite(bl.upper.separationX) && bl.upper.separationX < 0.95) ||
      (Number.isFinite(bl.lower.separationX) && bl.lower.separationX < 0.95),
    // Frottement local (rapporté à la vitesse extérieure locale ue) le long de
    // chaque face ; cf·ue² le rapporte à la vitesse amont.
    friction: {
      upper: { x: Array.from(bl.upperX), cf: Array.from(bl.upper.cf), ue: Array.from(bl.upperUe) },
      lower: { x: Array.from(bl.lowerX), cf: Array.from(bl.lower.cf), ue: Array.from(bl.lowerUe) },
    },
    // Cp aux points de contrôle, du bord de fuite par l'intrados au bord
    // d'attaque puis par l'extrados ; iLE : panneau le plus en avant.
    cp: { x: Array.from(flow.xc), y: Array.from(flow.yc), cp: Array.from(flow.cp), iLE: leadingPanel(flow.xc) },
    alphaDeg: o.alphaDeg,
  };
}

function leadingPanel(xc) {
  let k = 0;
  for (let i = 1; i < xc.length; i++) if (xc[i] < xc[k]) k = i;
  return k;
}

/**
 * Coefficients de section corrigés du décollement (écoulement de Kirchhoff,
 * forme de Beddoes et Leishman) : f, position du point de décollement de la
 * couche limite intégrale, réduit l'effort normal, Cn = Cn_att ((1 + √f)/2)²,
 * et l'aspiration de bord d'attaque, Cc = Cc_att √f. Sans décollement (f = 1)
 * on retrouve le résultat attaché. Modèle approché : il donne l'allure du
 * décrochage, pas son incidence exacte.
 * @param {object} r résultat de panelSolve
 * @returns {{cl:number, cd:number, cm:number, f:number, kirchhoff:boolean}}
 */
export function sectionCoefficients(r) {
  const fOf = (x) => (Number.isFinite(x) ? Math.min(1, Math.max(0, x / 0.97)) : 1);
  const f = r.separated ? Math.min(fOf(r.separationUpper), fOf(r.separationLower)) : 1;
  if (f >= 1) return { cl: r.cl, cd: r.cd, cm: r.cm, f: 1, kirchhoff: false };
  const a = (r.alphaDeg * Math.PI) / 180;
  const cnAtt = r.cl * Math.cos(a);
  const ccAtt = r.cl * Math.sin(a);
  const sf = Math.sqrt(f);
  const cn = cnAtt * ((1 + sf) / 2) ** 2;
  const cc = ccAtt * sf;
  const cl = cn * Math.cos(a) + cc * Math.sin(a);
  const cd = r.cd + Math.max(0, cn * Math.sin(a) - cc * Math.cos(a));
  // Le centre de poussée recule vers ≈ 0,4 c quand le décollement remonte.
  const cm = r.cm - cn * 0.15 * (1 - f);
  return { cl, cd, cm, f, kirchhoff: true };
}

/** Cp critique (écoulement localement sonique), γ = 1,4. */
export function criticalCp(mach) {
  const g = 1.4;
  const m2 = mach * mach;
  return (2 / (g * m2)) * (((2 + (g - 1) * m2) / (g + 1)) ** (g / (g - 1)) - 1);
}

/**
 * Corrections de compressibilité subsonique : Prandtl–Glauert pour Cz et Cm,
 * Kármán–Tsien pour Cp. `critical` : Cp minimal corrigé sous le Cp critique
 * (zone supersonique locale, hors du domaine de ces corrections).
 */
export function compressibleCorrection(mach, cpMin) {
  const m = Math.min(mach, 0.95);
  const beta = Math.sqrt(1 - m * m);
  const kt = (cp) => cp / (beta + ((m * m) / (1 + beta)) * (cp / 2));
  const cpMinC = Number.isFinite(cpMin) ? kt(cpMin) : NaN;
  return { beta, factor: 1 / beta, cp: kt, cpMin: cpMinC, cpStar: criticalCp(m), critical: cpMinC < criticalCp(m) };
}

/** Couche limite sur les deux faces à partir de la solution potentielle. */
function surfaceBL(pts, flow, re, nCrit, xTripU, xTripL, o2 = {}) {
  const N = flow.vt.length;
  // Point d'arrêt : changement de signe de la vitesse tangentielle le plus
  // proche du bord d'attaque (un profil épaissi de δ* a un culot au bord de
  // fuite, où la vitesse peut aussi changer de signe).
  let ile = 0;
  for (let i = 1; i < N; i++) if (flow.xc[i] < flow.xc[ile]) ile = i;
  let ist = ile;
  let bestGap = Infinity;
  for (let i = 0; i < N - 1; i++) {
    if (flow.vt[i] <= 0 && flow.vt[i + 1] > 0 && Math.abs(i - ile) < bestGap) {
      bestGap = Math.abs(i - ile);
      ist = i;
    }
  }
  const side = (indices) => {
    const n = indices.length + 1;
    const s = new Float64Array(n);
    const ue = new Float64Array(n);
    const xs = new Float64Array(n);
    // Premier point : le point d'arrêt (entre ist et ist + 1).
    const a = flow.vt[ist];
    const b = flow.vt[ist + 1];
    const f = a / (a - b);
    const xst = flow.xc[ist] + f * (flow.xc[ist + 1] - flow.xc[ist]);
    const yst = flow.yc[ist] + f * (flow.yc[ist + 1] - flow.yc[ist]);
    let px = xst;
    let py = yst;
    xs[0] = xst;
    for (let k = 0; k < indices.length; k++) {
      const i = indices[k];
      s[k + 1] = s[k] + Math.hypot(flow.xc[i] - px, flow.yc[i] - py);
      ue[k + 1] = Math.abs(flow.vt[i]);
      xs[k + 1] = flow.xc[i];
      px = flow.xc[i];
      py = flow.yc[i];
    }
    return { s, ue, xs, indices };
  };
  const upIdx = [];
  for (let i = ist + 1; i < N; i++) upIdx.push(i);
  const loIdx = [];
  for (let i = ist; i >= 0; i--) loIdx.push(i);
  const U = side(upIdx);
  const Lw = side(loIdx);
  const blU = boundaryLayer(U.s, U.ue, U.xs, { re, nCrit, xTrip: xTripU });
  const blL = boundaryLayer(Lw.s, Lw.ue, Lw.xs, { re, nCrit, xTrip: xTripL });
  const squire = (bl, side) => {
    // Évalué à x/c ≈ 0,98 : au-delà, la décélération potentielle vers le bord
    // de fuite fermé n'existe pas dans l'écoulement visqueux réel.
    let k = side.s.length - 1;
    while (k > 1 && side.xs[k] > (o2.squireX ?? 0.98)) k--;
    const th = bl.theta[k];
    const ueTE = Math.max(side.ue[k], 0.05);
    const H = Math.min(bl.H[k], 3);
    return 2 * th * ueTE ** ((H + 5) / 2);
  };
  // δ* aux nœuds du contour (moyenne des panneaux adjacents).
  const dPanel = new Float64Array(N);
  U.indices.forEach((i, k) => (dPanel[i] = blU.theta[k + 1] * blU.H[k + 1]));
  Lw.indices.forEach((i, k) => (dPanel[i] = blL.theta[k + 1] * blL.H[k + 1]));
  // Près du bord de fuite, l'écoulement potentiel d'un profil fermé décélère
  // jusqu'à l'arrêt et gonfle artificiellement δ* : on le fige à sa valeur à
  // x/c = 0,95 sur chaque face.
  const freeze = (side, bl) => {
    let ref = 0;
    side.indices.forEach((i, k) => {
      if (side.xs[k + 1] <= 0.95) ref = dPanel[i];
      else dPanel[i] = ref;
    });
    return bl;
  };
  freeze(U, blU);
  freeze(Lw, blL);
  const dstarNode = new Float64Array(N + 1);
  for (let k = 0; k <= N; k++) {
    const a = dPanel[Math.max(0, k - 1)];
    const b = dPanel[Math.min(N - 1, k)];
    dstarNode[k] = Math.min(0.05, 0.5 * (a + b));
  }
  return {
    upper: blU,
    lower: blL,
    upperX: U.xs,
    lowerX: Lw.xs,
    upperUe: U.ue,
    lowerUe: Lw.ue,
    cd: squire(blU, U) + squire(blL, Lw),
    dstarNode,
  };
}
