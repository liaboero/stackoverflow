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

let failures = 0;
const t0 = performance.now();

function check(name, ok, detail) {
  console.log(`${ok ? '✔' : '✘'} ${name} — ${detail}`);
  if (!ok) failures++;
}

/** Fait tourner une configuration et renvoie l'analyse des `avgT` derniers t*. */
function run({ shape = 'cylinder', re, L = 16, tEnd, avgT = 20, params = {} }) {
  const wt = new WindTunnel({ L, re, shape, params });
  const hist = new ForceHistory(32768);
  const spt = wt.stepsPerTStar;
  const every = Math.max(1, Math.round(spt / 40));
  const total = Math.round(tEnd * spt);
  for (let s = 1; s <= total; s++) {
    wt.step();
    if (s % every === 0) {
      const c = wt.coefficients();
      hist.push(wt.lbm.time, c.cx, c.cz, c.cm);
    }
  }
  const a = analyze(hist.since(wt.lbm.time - avgT * spt), wt.uIn, wt.L);
  a.finite = wt.lbm.isFinite();
  return a;
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

console.log(`\n${failures ? `${failures} échec(s)` : 'Tous les tests passent'} en ${((performance.now() - t0) / 1000).toFixed(0)} s`);
process.exit(failures ? 1 : 0);
