/**
 * Validation du solveur (Node ≥ 18, sans dépendance) : `npm test`.
 *
 * Les bornes tiennent compte du confinement du domaine (24 L × 12 L) : les
 * coefficients sont un peu supérieurs aux valeurs en milieu infini, rappelées
 * en commentaire.
 */
import { LBM } from '../js/lbm.js';
import { WindTunnel, RAMP_T } from '../js/sim.js';
import { ForceHistory, analyze } from '../js/stats.js';
import { panelSolve, potentialFlow, airfoilPoints, sectionCoefficients } from '../js/panel.js';
import { buildWing, solveWing, polarTable, polarAt } from '../js/vlm.js';

let failures = 0;
const t0 = performance.now();

function check(name, ok, detail) {
  console.log(`${ok ? '✔' : '✘'} ${name} — ${detail}`);
  if (!ok) failures++;
}

/** Fait tourner une configuration et renvoie l'analyse des `avgT` derniers t*. */
function run({ shape = 'cylinder', re, L = 16, tEnd, avgT = 20, params = {}, turbulence, transition = false }) {
  const wt = new WindTunnel({ L, re, shape, params, turbulence });
  if (transition) wt.setTransition(panelSolve(shape, { params, alphaDeg: params.alpha ?? 0, re }));
  const hist = new ForceHistory(32768);
  const spt = wt.stepsPerTStar;
  const every = Math.max(1, Math.round(spt / 40));
  const total = Math.round(tEnd * spt);
  for (let s = 1; s <= total; s++) {
    wt.step();
    if (s % every === 0) {
      const c = wt.averagedCoefficients();
      hist.push(wt.lbm.time, c.cx, c.cz, c.cm);
    }
  }
  const a = analyze(hist.since(wt.lbm.time - avgT * spt), wt.uIn, wt.L);
  a.finite = wt.lbm.isFinite();
  return a;
}

// 0. Méthode des panneaux + couche limite intégrale (instantané).
{
  const inv = potentialFlow(airfoilPoints('naca', { naca: '0012' }), (4 * Math.PI) / 180);
  check('Panneaux : NACA 0012 non visqueux, α = 4°', Math.abs(inv.cl - 0.48) < 0.01, `Cz = ${inv.cl.toFixed(4)} (≈ 0,48)`);
  const t0 = panelSolve('naca', { params: { naca: '0012' }, alphaDeg: 0, re: 6e6, xTrip: 0 });
  const t10 = panelSolve('naca', { params: { naca: '0012' }, alphaDeg: 10, re: 6e6, xTrip: 0 });
  check(
    'Panneaux + CL : NACA 0012, Re = 6·10⁶ turbulent',
    Math.abs(t0.cd - 0.0081) < 0.0008 && Math.abs(t10.cl - 1.09) < 0.12 && Math.abs(t10.cd - 0.0123) < 0.002,
    `Cx(0°) = ${t0.cd.toFixed(4)} (NASA TMR 0,0081) ; Cz, Cx(10°) = ${t10.cl.toFixed(3)}, ${t10.cd.toFixed(4)} (1,09, 0,0123)`,
  );
  const free = panelSolve('naca', { params: { naca: '0012' }, alphaDeg: 0, re: 1e6 });
  check(
    'Panneaux + CL : transition naturelle, Re = 10⁶',
    free.transitionUpper > 0.4 && free.transitionUpper < 0.8 && free.cd > 0.0045 && free.cd < 0.0075,
    `x_tr = ${free.transitionUpper.toFixed(2)}, Cx = ${free.cd.toFixed(4)} (XFOIL ≈ 0,0054)`,
  );
}

// 0 bis. Aile finie (surface portante de Weissinger non linéaire).
{
  const lin = (i, a) => ({ cl: (2 * Math.PI * a * Math.PI) / 180, cd: 0, cm: 0 });
  const slope = (w) => ((solveWing(w, 5, lin).CL - solveWing(w, 3, lin).CL) / 2) * (180 / Math.PI);
  const rect = buildWing({ ar: 6, taper: 1 });
  const r6 = solveWing(rect, 4, lin);
  const a6 = slope(rect);
  check('Aile rectangulaire A = 6 (profil 2π)', a6 > 4.0 && a6 < 4.5 && r6.e > 0.95 && r6.e < 1.03, `CLα = ${a6.toFixed(2)}/rad (surface portante ≈ 4,2 ; ligne portante 4,5), e = ${r6.e.toFixed(3)}`);
  const tap = solveWing(buildWing({ ar: 8, taper: 0.4 }), 4, lin);
  check('Aile effilée λ = 0,4 : quasi elliptique', tap.e > 0.98 && tap.e < 1.03, `e = ${tap.e.toFixed(3)}`);
  const sw = solveWing(buildWing({ ar: 6, taper: 0.5, sweep: 30 }), 4, lin);
  const root = sw.strips[sw.strips.length / 2].ccl;
  const peak = Math.max(...sw.strips.map((s) => s.ccl));
  check('Aile en flèche 30° : charge régulière à l’emplanture', sw.e > 0.9 && sw.e < 1.03 && root > 0.9 * peak, `e = ${sw.e.toFixed(3)}, charge emplanture/max = ${(root / peak).toFixed(2)}`);
  const tab = polarTable((a) => sectionCoefficients(panelSolve('naca', { params: { naca: '2412' }, alphaDeg: a, re: 1e6 })), -6, 22);
  const nl = solveWing(buildWing({ ar: 8, taper: 0.5 }), 4, (i, a) => polarAt(tab, a));
  check('Aile NACA 2412, A = 8, Re = 10⁶, α = 4°', nl.converged && nl.CL > 0.4 && nl.CL < 0.65 && nl.CDv > 0.004 && nl.CDv < 0.012, `CL = ${nl.CL.toFixed(3)}, CDi = ${nl.CDi.toFixed(4)}, CDp = ${nl.CDv.toFixed(4)}`);
}

// 1. Écoulement uniforme sans obstacle : conditions d'entrée, de sortie et
//    bandes absorbantes ne doivent pas le perturber.
{
  const lbm = new LBM({ nx: 120, ny: 60, uIn: 0.1, nu: 0.02 });
  lbm.setGeometry(() => false);
  lbm.reset();
  for (let k = 0; k < 2000; k++) lbm.step();
  lbm.computeMacros();
  let err = 0;
  for (let i = 0; i < lbm.N; i++) err = Math.max(err, Math.abs(lbm.ux[i] - 0.1), Math.abs(lbm.uy[i]), Math.abs(lbm.rho[i] - 1));
  check('Écoulement uniforme', err < 1e-4, `écart max ${err.toExponential(2)}`);
}

// 2. Cylindre, Re = 20 : écoulement stationnaire symétrique (Cx ≈ 2,05 en milieu infini).
{
  const a = run({ re: 20, tEnd: RAMP_T + 30, avgT: 10 });
  check(
    'Cylindre Re = 20',
    a.finite && !a.periodic && a.cx > 2.1 && a.cx < 2.6 && Math.abs(a.cz) < 0.01,
    `Cx = ${a.cx.toFixed(3)}, Cz = ${a.cz.toFixed(4)}`,
  );
}

// 3. Cylindre, Re = 100 : allée de von Kármán (Cx ≈ 1,33, St ≈ 0,165, Cz' ≈ 0,33 en milieu infini).
{
  const a = run({ re: 100, tEnd: 90, avgT: 25 });
  check(
    'Cylindre Re = 100',
    a.finite && a.periodic && a.cx > 1.35 && a.cx < 1.6 && a.st > 0.155 && a.st < 0.19 && a.czAmp > 0.25 && a.czAmp < 0.45,
    `Cx = ${a.cx.toFixed(3)}, St = ${a.st.toFixed(3)}, Cz' = ${a.czAmp.toFixed(3)}`,
  );
}

// 4. NACA 0012 : portance nulle à incidence nulle, antisymétrie en ±α.
{
  const p = (alpha) => run({ shape: 'naca', re: 500, L: 12, tEnd: RAMP_T + 30, avgT: 10, params: { naca: '0012', alpha } });
  const a0 = p(0);
  const ap = p(4);
  const am = p(-4);
  check('NACA 0012, α = 0°', a0.finite && Math.abs(a0.cz) < 0.005, `Cz = ${a0.cz.toFixed(4)}, Cx = ${a0.cx.toFixed(3)}`);
  check(
    'NACA 0012, α = ±4°',
    ap.cz > 0.08 && Math.abs(ap.cz + am.cz) < 0.02 * ap.cz + 0.003 && Math.abs(ap.cx - am.cx) < 0.01,
    `Cz(+4°) = ${ap.cz.toFixed(4)}, Cz(−4°) = ${am.cz.toFixed(4)}`,
  );
}

// 4 ter. Linéarité de la portance aux petites incidences : la pointe du bord
//    de fuite est recalée sur le réseau, sinon Cz(α) varie en dents de scie.
{
  const cz = [1, 2, 3].map((alpha) => run({ shape: 'naca', re: 1000, L: 16, tEnd: RAMP_T + 25, avgT: 8, params: { naca: '0012', alpha } }).cz);
  const slopes = cz.map((c, k) => c / (k + 1));
  check(
    'NACA 0012, Cz(α) régulier',
    cz[0] < cz[1] && cz[1] < cz[2] && slopes.every((v) => v > 0.025 && v < 0.065),
    `Cz(1°, 2°, 3°) = ${cz.map((c) => c.toFixed(3)).join(', ')}`,
  );
}

// 4 bis. Traînée de frottement : le modèle de sous-maille, amorti près de la
//    paroi, ne doit pas épaissir la couche limite (Cx ≈ 0,144 sans LES ;
//    référence fine 0,119, Di Ilio et al. 2018 et XFOIL).
{
  const a = run({ shape: 'naca', re: 1000, L: 16, tEnd: RAMP_T + 30, avgT: 10, params: { naca: '0012', alpha: 0 } });
  check('NACA 0012, Re = 1000, LES amortie', a.finite && a.cx > 0.119 && a.cx < 0.155, `Cx = ${a.cx.toFixed(4)}`);
}

// 5. Cylindre tournant ΩR/U = 1, Re = 100 : effet Magnus (Cz ≈ 2,5 d'après Mittal & Kumar 2003).
{
  const a = run({ re: 100, L: 12, tEnd: 70, avgT: 25, params: { spinRatio: 1 } });
  check('Effet Magnus ΩR/U = 1', a.finite && a.cz > 2 && a.cz < 3.2, `Cz = ${a.cz.toFixed(3)}, Cx = ${a.cx.toFixed(3)}`);
}

// 6. Stabilité à Reynolds élevé avec le modèle de sous-maille.
{
  const a = run({ re: 5000, L: 12, tEnd: 30, avgT: 10 });
  check('Re = 5000 (LES)', a.finite && Number.isFinite(a.cx), `Cx = ${a.cx.toFixed(3)}`);
}

// 7. Spalart–Allmaras + loi de paroi, transition imposée par la méthode des
//    panneaux. Seules la portance et la stabilité sont vérifiées : la traînée
//    LBM d'un profil mince n'est pas convergée en maillage à ces résolutions
//    (voir README, « Étape A »).
{
  const a = run({ shape: 'naca', re: 1e6, L: 16, tEnd: 20, avgT: 5, turbulence: 'sa', transition: true, params: { naca: '0012', alpha: 4 } });
  check('SA, NACA 0012, Re = 10⁶, α = 4°, transition libre', a.finite && a.cz > 0.42 && a.cz < 0.56, `Cz = ${a.cz.toFixed(3)} (panneaux 0,471 ; mesures ≈ 0,44), Cx = ${a.cx.toFixed(4)}`);
  const b = run({ shape: 'naca', re: 1e7, L: 16, tEnd: 15, avgT: 5, turbulence: 'sa', params: { naca: '0012', alpha: 8 } });
  check('SA, NACA 0012, Re = 10⁷, α = 8° (stabilité)', b.finite && b.cz > 0.75 && b.cz < 1.1, `Cz = ${b.cz.toFixed(3)} (≈ 0,85 à 0,9)`);
}

console.log(`\n${failures ? `${failures} échec(s)` : 'Tous les tests passent'} en ${((performance.now() - t0) / 1000).toFixed(0)} s`);
process.exit(failures ? 1 : 0);
