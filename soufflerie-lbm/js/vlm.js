/**
 * Aile finie : surface portante de Weissinger (un tourbillon en fer à cheval
 * par bande d'envergure, segment lié au quart de corde, condition au point
 * des trois quarts de corde), rendue non linéaire par « décambrage » : la
 * circulation de chaque bande vérifie la loi de Kutta–Joukowski avec la
 * polaire 2D de sa section, à l'incidence effective α(3/4 c) + Cl/a0. En
 * régime linéaire, on retrouve exactement la méthode de Weissinger ; au-delà,
 * la polaire (décrochage compris) fixe la charge. Flèche, dièdre, effilement
 * et vrillage sont pris en compte. Système non linéaire résolu par Newton.
 *
 * Repère : x vers l'aval, y vers la droite (envergure), z vers le haut.
 * Corde moyenne géométrique c̄ = 1, donc envergure b = allongement et S = b.
 */

const FOUR_PI = 4 * Math.PI;
const DEG = Math.PI / 180;

const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const norm = (a) => Math.hypot(a[0], a[1], a[2]);
const scale = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
const unit = (a) => scale(a, 1 / (norm(a) || 1));

/**
 * Géométrie en plan trapézoïdale.
 * @param {object} o { ar, taper, sweep (°, quart de corde), dihedral (°),
 *        twist (° en bout, linéaire, négatif = vrillage « wash-out »), n (bandes par demi-aile) }
 */
export function buildWing(o) {
  const ar = o.ar;
  const taper = o.taper;
  const n = o.n ?? 40;
  const b = ar;
  const half = b / 2;
  const cr = 2 / (1 + taper);
  const tanS = Math.tan((o.sweep ?? 0) * DEG);
  const tanD = Math.tan((o.dihedral ?? 0) * DEG);
  const twistTip = (o.twist ?? 0) * DEG;
  const chord = (y) => cr * (1 - ((1 - taper) * Math.abs(y)) / half);
  const c4 = (y) => [Math.abs(y) * tanS, y, Math.abs(y) * tanD];
  const ys = [];
  // Bandes uniformes par défaut : une bande très étroite devant la corde (la
  // répartition en cosinus en crée en bout d'aile) rend la condition aux
  // trois quarts de corde mal posée et fausse la traînée induite.
  const blend = o.blend ?? 0;
  for (let k = 0; k <= 2 * n; k++) {
    const u = k / (2 * n);
    ys.push(half * (blend * -Math.cos(Math.PI * u) + (1 - blend) * (2 * u - 1)));
  }
  const strips = [];
  for (let i = 0; i < 2 * n; i++) {
    const ya = ys[i];
    const yb = ys[i + 1];
    const ym = 0.5 * (ya + yb);
    const A = c4(ya);
    const B = c4(yb);
    const P = scale([A[0] + B[0], A[1] + B[1], A[2] + B[2]], 0.5);
    const dl = sub(B, A);
    const us = unit(dl);
    const ua0 = [1, 0, 0];
    const un0 = unit(cross(ua0, us));
    const th = (twistTip * Math.abs(ym)) / half;
    const ua = [
      Math.cos(th) * ua0[0] - Math.sin(th) * un0[0],
      Math.cos(th) * ua0[1] - Math.sin(th) * un0[1],
      Math.cos(th) * ua0[2] - Math.sin(th) * un0[2],
    ];
    const un = [
      Math.sin(th) * ua0[0] + Math.cos(th) * un0[0],
      Math.sin(th) * ua0[1] + Math.cos(th) * un0[1],
      Math.sin(th) * ua0[2] + Math.cos(th) * un0[2],
    ];
    const c34 = 0.5 * chord(ym);
    const P34 = [P[0] + c34 * ua[0], P[1] + c34 * ua[1], P[2] + c34 * ua[2]];
    strips.push({ A, B, P, P34, dl, ua, un, us, c: chord(ym), ds: Math.hypot(dl[1], dl[2]), eta: ym / half, twist: th / DEG });
  }
  const mac = ((2 / 3) * cr * (1 + taper + taper * taper)) / (1 + taper);
  const yMac = ((b / 6) * (1 + 2 * taper)) / (1 + taper);
  return {
    ...o,
    b,
    S: b,
    cr,
    ct: cr * taper,
    mac,
    ref: [yMac * tanS, 0, 0],
    strips,
    chordAt: chord,
    c4At: c4,
  };
}

/** Vitesse induite en P par un segment tourbillonnaire A → B de circulation 1. */
function segment(P, A, B) {
  const r1 = sub(P, A);
  const r2 = sub(P, B);
  const c = cross(r1, r2);
  const c2 = dot(c, c);
  const n1 = norm(r1);
  const n2 = norm(r2);
  if (c2 < 1e-10 || n1 < 1e-9 || n2 < 1e-9) return [0, 0, 0];
  const r0 = sub(B, A);
  const k = (dot(r0, r1) / n1 - dot(r0, r2) / n2) / (FOUR_PI * c2);
  return scale(c, k);
}

/** Tourbillon semi-infini partant de Q dans la direction d (unitaire), circulation 1. */
function semiInfinite(P, Q, d) {
  const r = sub(P, Q);
  const c = cross(d, r);
  const c2 = dot(c, c);
  const nr = norm(r);
  if (c2 < 1e-10 || nr < 1e-9) return [0, 0, 0];
  const k = (1 + dot(r, d) / nr) / (FOUR_PI * c2);
  return scale(c, k);
}

/** Influence d'un fer à cheval (tourbillons libres selon d) en P. */
function horseshoe(P, s, d) {
  const v1 = segment(P, s.A, s.B);
  const v2 = semiInfinite(P, s.B, d);
  const v3 = semiInfinite(P, s.A, d);
  return [v1[0] + v2[0] - v3[0], v1[1] + v2[1] - v3[1], v1[2] + v2[2] - v3[2]];
}

/**
 * Résout l'aile à une incidence.
 * @param {object} wing résultat de buildWing
 * @param {number} alphaDeg incidence de l'emplanture (°)
 * @param {(i:number, alphaDeg:number) => {cl:number, cd:number, cm:number}} section
 *        polaire 2D de la bande i
 * @param {object} [o] { gamma0 (solution initiale), maxIter }
 */
export function solveWing(wing, alphaDeg, section, o = {}) {
  const { strips } = wing;
  const N = strips.length;
  const a = alphaDeg * DEG;
  const Vinf = [Math.cos(a), 0, Math.sin(a)];
  // Sillage plan, dans le prolongement de la corde (hypothèse classique de la
  // surface portante linéarisée ; un sillage aligné sur l'écoulement amont
  // ajouterait une non-linéarité de quelques pour cent).
  const d = [1, 0, 0];
  // Influences des fers à cheval : au milieu du segment lié (force de
  // Kutta–Joukowski) et au point aux trois quarts de corde (incidence locale,
  // comme dans la méthode de Weissinger : pas de singularité au raccord d'une
  // aile en flèche).
  const Vb = new Array(N);
  const V3 = new Array(N);
  for (let i = 0; i < N; i++) {
    Vb[i] = new Array(N);
    V3[i] = new Array(N);
    for (let j = 0; j < N; j++) {
      Vb[i][j] = horseshoe(strips[i].P, strips[j], d);
      V3[i][j] = horseshoe(strips[i].P34, strips[j], d);
    }
  }
  // Partie linéaire de chaque polaire : pente a0 (par radian) et incidence de
  // portance nulle α0.
  const a0 = new Float64Array(N);
  const al0 = new Float64Array(N);
  for (let i = 0; i < N; i++) {
    const c1 = section(i, -2).cl;
    const c2 = section(i, 4).cl;
    const sl = Math.max((c2 - c1) / 6, 0.02);
    a0[i] = sl / DEG;
    al0[i] = (-2 - c1 / sl) * DEG;
  }
  const G = o.gamma0 && o.gamma0.length === N ? Float64Array.from(o.gamma0) : new Float64Array(N);
  const wb = new Array(N);
  const w3 = new Array(N);
  const aeff = new Float64Array(N);
  const R = new Float64Array(N);
  // Résidu d'une bande : Cl de Kutta–Joukowski moins Cl de la polaire à
  // l'incidence effective α_e = α(3/4 c) + Cl/a0 (méthode de « décambrage » :
  // en régime linéaire, elle se réduit à la condition de Weissinger
  // α(3/4 c) = α0 ; au-delà, elle suit la polaire non linéaire).
  const local = (i, wbv, w3v, g) => {
    const s = strips[i];
    const ck = (2 * g * norm(cross(wbv, s.dl))) / (s.c * s.ds);
    const a34 = Math.atan2(dot(w3v, s.un), dot(w3v, s.ua));
    const ae = a34 + ck / a0[i];
    return [ae, ck - section(i, ae / DEG).cl];
  };
  const acc = (M, i, Gv) => {
    let x = Vinf[0];
    let y = Vinf[1];
    let z = Vinf[2];
    for (let j = 0; j < N; j++) {
      const v = M[i][j];
      x += v[0] * Gv[j];
      y += v[1] * Gv[j];
      z += v[2] * Gv[j];
    }
    return [x, y, z];
  };
  const evaluate = (Gv) => {
    let r2 = 0;
    for (let i = 0; i < N; i++) {
      wb[i] = acc(Vb, i, Gv);
      w3[i] = acc(V3, i, Gv);
      const [ae, r] = local(i, wb[i], w3[i], Gv[i]);
      aeff[i] = ae;
      R[i] = r;
      r2 += r * r;
    }
    return Math.sqrt(r2 / N);
  };
  let res = evaluate(G);
  let iter = 0;
  const maxIter = o.maxIter ?? 40;
  const J = Array.from({ length: N }, () => new Float64Array(N));
  const h = 1e-6;
  while (res > 1e-9 && iter < maxIter) {
    iter++;
    // Jacobien par différences finies, colonne par colonne (O(N²)).
    for (let j = 0; j < N; j++) {
      for (let i = 0; i < N; i++) {
        const vb = Vb[i][j];
        const v3 = V3[i][j];
        const wbp = [wb[i][0] + h * vb[0], wb[i][1] + h * vb[1], wb[i][2] + h * vb[2]];
        const w3p = [w3[i][0] + h * v3[0], w3[i][1] + h * v3[1], w3[i][2] + h * v3[2]];
        J[i][j] = (local(i, wbp, w3p, G[i] + (i === j ? h : 0))[1] - R[i]) / h;
      }
    }
    const dG = solveLinear(J.map((r) => Float64Array.from(r)), Float64Array.from(R, (r) => -r));
    if (!dG) break;
    // Newton amorti : on réduit le pas tant que le résidu ne baisse pas.
    let step = 1;
    const G0 = Float64Array.from(G);
    for (let t = 0; t < 8; t++) {
      for (let i = 0; i < N; i++) G[i] = G0[i] + step * dG[i];
      const r = evaluate(G);
      if (r < res || t === 7) {
        res = r;
        break;
      }
      step *= 0.5;
    }
  }
  const w = wb;
  const alpha = aeff;
  // Efforts : Kutta–Joukowski avec la vitesse locale (portance + traînée induite),
  // traînée de profil des polaires de section.
  const eL = [-Math.sin(a), 0, Math.cos(a)];
  let L = 0;
  let Dv = 0;
  let My = 0;
  const per = [];
  for (let i = 0; i < N; i++) {
    const s = strips[i];
    const F = scale(cross(w[i], s.dl), G[i]);
    const sec = section(i, alpha[i] / DEG);
    L += dot(F, eL);
    const dv = 0.5 * sec.cd * s.c * s.ds;
    Dv += dv;
    const r = sub(s.P, wing.ref);
    My += r[2] * F[0] - r[0] * F[2] + 0.5 * sec.cm * s.c * s.c * s.ds * s.us[1];
    per.push({
      eta: s.eta,
      c: s.c,
      cl: sec.cl,
      cd: sec.cd,
      alpha: alpha[i] / DEG,
      // Charge locale c·Cl/c̄ (c̄ = 1), au sens de Kutta–Joukowski.
      ccl: (2 * dot(F, eL)) / s.ds,
    });
  }
  // Traînée induite dans le plan de Trefftz (plus robuste que Kutta–Joukowski
  // pour les ailes en flèche) : tourbillons libres aux nœuds d'envergure,
  // d'intensité égale au saut de circulation.
  const DiT = trefftzDrag(strips, G);
  const q = 0.5;
  const S = wing.S;
  const CL = L / (q * S);
  const CDi = DiT / (q * S);
  const CDv = Dv / (q * S);
  return {
    alpha: alphaDeg,
    CL,
    CDi,
    CDv,
    CD: CDi + CDv,
    Cm: My / (q * S * wing.mac),
    e: CDi > 1e-9 ? (CL * CL) / (Math.PI * wing.ar * CDi) : NaN,
    strips: per,
    gamma: Array.from(G),
    converged: res < 1e-6,
    residual: res,
    iterations: iter,
  };
}

/** D_i = −½ ρ Σ Γ_i w_n,i ds_i, w_n vitesse normale à la nappe (ρ = V∞ = 1). */
function trefftzDrag(strips, G) {
  const N = strips.length;
  // Nœuds (y, z) et tourbillons libres (vorticité selon +x).
  const ny = [];
  const nz = [];
  const gam = [];
  for (let k = 0; k <= N; k++) {
    const p = k < N ? strips[k].A : strips[N - 1].B;
    ny.push(p[1]);
    nz.push(p[2]);
    gam.push((k > 0 ? G[k - 1] : 0) - (k < N ? G[k] : 0));
  }
  let D = 0;
  for (let i = 0; i < N; i++) {
    const s = strips[i];
    const ym = s.P[1];
    const zm = s.P[2];
    let vy = 0;
    let vz = 0;
    for (let k = 0; k <= N; k++) {
      const ry = ym - ny[k];
      const rz = zm - nz[k];
      const r2 = ry * ry + rz * rz;
      if (r2 < 1e-14) continue;
      const f = gam[k] / (2 * Math.PI * r2);
      vy -= f * rz;
      vz += f * ry;
    }
    const sy = s.dl[1] / s.ds;
    const sz = s.dl[2] / s.ds;
    const wn = -sz * vy + sy * vz;
    D -= 0.5 * G[i] * wn * s.ds;
  }
  return D;
}

/** Élimination de Gauss avec pivot partiel ; null si singulier. */
function solveLinear(A, b) {
  const n = b.length;
  for (let c = 0; c < n; c++) {
    let p = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(A[r][c]) > Math.abs(A[p][c])) p = r;
    if (Math.abs(A[p][c]) < 1e-14) return null;
    [A[c], A[p]] = [A[p], A[c]];
    [b[c], b[p]] = [b[p], b[c]];
    for (let r = c + 1; r < n; r++) {
      const f = A[r][c] / A[c][c];
      if (f === 0) continue;
      for (let k = c; k < n; k++) A[r][k] -= f * A[c][k];
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
 * Table de polaire 2D (incidences entières) et interpolation linéaire.
 * @param {(alphaDeg:number) => {cl, cd, cm}} fn
 */
export function polarTable(fn, lo = -14, hi = 24) {
  const al = [];
  const cl = [];
  const cd = [];
  const cm = [];
  for (let a = lo; a <= hi; a++) {
    const r = fn(a);
    al.push(a);
    cl.push(r.cl);
    cd.push(r.cd);
    cm.push(r.cm);
  }
  return { lo, hi, al, cl, cd, cm };
}

/** Valeur interpolée d'une table ; au-delà, prolongement de plaque plane décollée. */
export function polarAt(t, aDeg) {
  if (aDeg >= t.lo && aDeg <= t.hi) {
    const k = Math.min(t.al.length - 2, Math.floor(aDeg - t.lo));
    const f = aDeg - t.lo - k;
    return {
      cl: t.cl[k] + f * (t.cl[k + 1] - t.cl[k]),
      cd: t.cd[k] + f * (t.cd[k + 1] - t.cd[k]),
      cm: t.cm[k] + f * (t.cm[k + 1] - t.cm[k]),
    };
  }
  // Raccord progressif vers la plaque plane décollée (Cn ≈ 2 sin α).
  const end = aDeg > t.hi ? t.al.length - 1 : 0;
  const a0 = t.al[end];
  const r = aDeg * DEG;
  const cn = 1.98 * Math.sin(r);
  const blend = Math.min(1, Math.abs(aDeg - a0) / 15);
  return {
    cl: (1 - blend) * t.cl[end] + blend * cn * Math.cos(r),
    cd: (1 - blend) * t.cd[end] + blend * (cn * Math.sin(r) + 0.02),
    cm: (1 - blend) * t.cm[end] + blend * (-0.25 * cn),
  };
}
