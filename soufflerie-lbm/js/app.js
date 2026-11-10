/**
 * Soufflerie de poche — interface mobile.
 */
import { WindTunnel, RESOLUTIONS, RAMP_T, RE_WAKE_3D } from './sim.js';
import { SHAPES, SHAPE_BY_ID, parseNaca, shapeIconPath } from './shapes.js';
import { ForceHistory, analyze } from './stats.js';
import { FlowView, FIELDS, VIEWS, gradientCss } from './render.js';
import { drawHistory, drawPolar, drawWall } from './chart.js';
import { GpuLBM, gpuSupported } from './lbm-gpu.js';
import { GpuFlowView } from './render-gpu.js';
import { panelSolve, sectionCoefficients, compressibleCorrection } from './panel.js';
import { WallSampler, WALL_SHAPES, upstreamDensity } from './wall.js';
import { WingTab } from './wing-tab.js';
import { Lab3D } from './lab3d.js';

const RE_VALUES = [
  10, 20, 40, 60, 100, 150, 200, 300, 500, 1000, 1500, 2000, 3000, 5000, 10000, 15000, 20000, 50000, 1e5, 2e5,
  5e5, 1e6, 2e6, 5e6, 1e7,
];
/** Reynolds à partir duquel « Auto » active Spalart–Allmaras et la loi de paroi. */
const RE_SA_AUTO = 50000;
/** Formes traitées par la méthode des panneaux et la couche limite intégrale. */
const PANEL_SHAPES = ['naca', 'ellipse', 'plate'];
/** Reynolds minimal de la référence « panneaux + couche limite ». */
const RE_PANEL_MIN = 20000;
const U_VALUES = [0.05, 0.08, 0.1, 0.12];
/** Facteurs de vitesse : temps physique simulé par seconde écoulée. */
const SPEEDS = ['max', 10, 2, 1, 0.5, 0.1];
/** Fluides : viscosité cinématique (m²/s) et masse volumique (kg/m³). */
const FLUIDS = {
  air: { label: 'Air 15 °C', nu: 1.46e-5, rho: 1.225 },
  water: { label: 'Eau 20 °C', nu: 1.0e-6, rho: 998 },
};
const WINDOW_T = 20; // fenêtre glissante de moyenne (t*)
const SETTLE_T = 6; // transitoire ignoré après un changement (t*)
const HISTORY_SPAN_T = 40; // largeur du graphique (t*)
const SWEEP_SETTLE_T = 12; // durées maximales (t*)
const SWEEP_MEASURE_T = 15;
const SWEEP_SETTLE_MIN_T = 5; // durées minimales si l'écoulement a convergé
const SWEEP_MEASURE_MIN_T = 6;
const FRAME_BUDGET_MS = 11;
const PREFS_KEY = 'soufflerie.prefs.v1';
const POINTS_KEY = 'soufflerie.points.v1';

const $ = (id) => document.getElementById(id);

// ---------------------------------------------------------------- état

const DEFAULTS = {
  shape: 'cylinder',
  alphaByShape: {},
  naca: '0012',
  thickness: 0.4,
  spin: 0,
  re: 100,
  res: 'standard',
  field: 'vort',
  tracers: true,
  arrow: true,
  wall: 'bouzidi',
  les: 'auto',
  u: 0.1,
  polarMode: 'var',
  view: 'near',
  engine: 'auto',
  collision: 'rr',
  speed: 'max',
  chord: 0.1,
  fluid: 'air',
  turb: 'auto',
  transition: 'free',
  wallQ: 'cp',
  tab: '2d',
};

function storageGet(key) {
  try {
    const v = localStorage.getItem(key);
    return v ? JSON.parse(v) : null;
  } catch {
    return null;
  }
}

function storageSet(key, value) {
  try {
    localStorage.setItem(key, JSON.stringify(value));
  } catch {
    /* stockage indisponible : sans conséquence */
  }
}

const state = { ...DEFAULTS, ...(storageGet(PREFS_KEY) || {}) };
if (!SHAPE_BY_ID[state.shape]) state.shape = DEFAULTS.shape;
if (!RESOLUTIONS[state.res]) state.res = DEFAULTS.res;
if (!RE_VALUES.includes(state.re)) state.re = DEFAULTS.re;
if (!U_VALUES.includes(state.u)) state.u = DEFAULTS.u;
if (!FIELDS[state.field]) state.field = DEFAULTS.field;
if (!VIEWS[state.view]) state.view = DEFAULTS.view;
if (!SPEEDS.includes(state.speed)) state.speed = DEFAULTS.speed;
if (!FLUIDS[state.fluid]) state.fluid = DEFAULTS.fluid;
if (!(state.chord > 0 && state.chord <= 100)) state.chord = DEFAULTS.chord;
if (!parseNaca(state.naca)) state.naca = DEFAULTS.naca;
if (!['auto', 'none', 'sa'].includes(state.turb)) state.turb = DEFAULTS.turb;
if (!['free', 'turb'].includes(state.transition)) state.transition = DEFAULTS.transition;
if (!['cp', 'cf'].includes(state.wallQ)) state.wallQ = DEFAULTS.wallQ;
if (!['2d', 'wing', '3d'].includes(state.tab)) state.tab = DEFAULTS.tab;

let points = storageGet(POINTS_KEY) || [];

let savePending = 0;
function savePrefs() {
  clearTimeout(savePending);
  savePending = setTimeout(() => storageSet(PREFS_KEY, state), 400);
}

const shapeDef = () => SHAPE_BY_ID[state.shape];

function alphaOf(id = state.shape) {
  const s = SHAPE_BY_ID[id];
  if (!s.alpha) return 0;
  const v = state.alphaByShape[id];
  return Number.isFinite(v) ? Math.min(s.alpha[1], Math.max(s.alpha[0], v)) : s.alpha[2];
}

function bodyParams() {
  return { alpha: alphaOf(), naca: state.naca, thickness: state.thickness, spinRatio: state.spin };
}

/** Modèle de turbulence effectif : Spalart–Allmaras ou aucun (LES seule). */
function turbMode() {
  if (state.turb === 'auto') return state.re >= RE_SA_AUTO ? 'sa' : 'none';
  return state.turb;
}

/** Reynolds lisible : 20 000, 2·10⁵, 10⁶… */
function fmtRe(v) {
  if (v < 1e5) return v.toLocaleString('fr-FR');
  const e = Math.floor(Math.log10(v) + 1e-9);
  const m = Math.round(v / 10 ** e);
  const sup = String(e).replace(/\d/g, (d) => '⁰¹²³⁴⁵⁶⁷⁸⁹'[d]);
  return m === 1 ? `10${sup}` : `${m}·10${sup}`;
}

/** Variable de la polaire : incidence, ou rotation pour le cylindre. */
function sweepVar() {
  return shapeDef().alpha ? 'alpha' : 'spin';
}

// ---------------------------------------------------------------- simulation

let wt;
let hist;
let view;
let running = true;
let diverged = false;
let tChange = 0; // pas de temps du dernier changement de configuration
let markers = []; // t* des changements (graphique)
let sampleEvery = 4;
let lastCoeffs = { cx: NaN, cz: NaN, cm: NaN };
let analysis = null;
let pendingBody = false;
let stepsSinceRender = 0;
let sweep = null;
let wakeLock = null;

// Moteur GPU (WebGL2) si disponible, sinon CPU.
let gl = null;
let gpuOk = false;
let cpuView = null;
let gpuView = null;
let gpuSteps = 16;

function initGpu() {
  try {
    gl = $('flow-gl').getContext('webgl2', {
      alpha: false,
      antialias: false,
      depth: false,
      stencil: false,
      preserveDrawingBuffer: false,
      powerPreference: 'high-performance',
    });
    gpuOk = gpuSupported(gl);
    if (gpuOk) gpuView = new GpuFlowView($('flow-gl'), $('flow'), gl);
  } catch {
    gpuOk = false;
    gpuView = null;
  }
}

function engine() {
  return gpuOk && state.engine !== 'cpu' ? 'gpu' : 'cpu';
}

function resAllowed(key) {
  return !!RESOLUTIONS[key] && (engine() === 'gpu' || !RESOLUTIONS[key].gpuOnly);
}

function buildSim() {
  if (wt) wt.dispose();
  if (!resAllowed(state.res)) state.res = 'tresfine';
  const eng = engine();
  const L = RESOLUTIONS[state.res].L;
  const cfg = {
    L,
    re: state.re,
    uIn: state.u,
    shape: state.shape,
    params: bodyParams(),
    les: state.les,
    bouzidi: state.wall === 'bouzidi',
    collision: state.collision,
    turbulence: turbMode(),
  };
  try {
    wt = new WindTunnel(eng === 'gpu' ? { ...cfg, gpu: { Solver: GpuLBM, gl } } : cfg);
  } catch (e) {
    if (eng !== 'gpu') throw e;
    // GPU refusé (mémoire, taille de texture…) : repli sur le CPU.
    gpuOk = false;
    wt = null;
    buildResSeg();
    buildSim();
    showToast(`Calcul GPU impossible (${e.message}) : passage au CPU.`, 5000);
    return;
  }
  view = wt.engine === 'gpu' ? gpuView : cpuView;
  $('flow-gl').hidden = wt.engine !== 'gpu';
  gpuSteps = 16;
  hist = new ForceHistory(16384);
  // 40 échantillons par t* sur CPU (moyenne de tous les pas, gratuite) ; 16 sur
  // GPU, chacun moyenné sur deux pas (quatre tracés) : assez pour un Strouhal
  // jusqu'à ~3.
  sampleEvery = Math.max(2, Math.round(wt.stepsPerTStar / (wt.engine === 'gpu' ? 16 : 40)));
  halfSample = null;
  tChange = 0;
  markers = [];
  analysis = null;
  diverged = false;
  lastCoeffs = { cx: NaN, cz: NaN, cm: NaN };
  view.attach(wt);
  hideToast();
  wallSampler = null;
  refreshPanel();
  updateInfo();
}

function markChange() {
  wallSampler?.reset();
  tChange = wt.lbm.time;
  markers.push(wt.tStar);
  if (markers.length > 40) markers.shift();
}

function applyBody() {
  wt.setBody(state.shape, bodyParams());
  wallSampler = null;
  refreshPanel();
  markChange();
  updateInfo();
}

// ---------------------------------------------------------------- méthode des panneaux

/**
 * Méthode des panneaux + couche limite intégrale (js/panel.js) : référence
 * en écoulement attaché, et position de la transition imposée au calcul LBM
 * quand Spalart–Allmaras est actif (zone laminaire en amont).
 */
let panelRef = null;
const panelCache = new Map();

function panelKey(alpha = alphaOf()) {
  const s = shapeDef();
  const extra = s.naca ? state.naca : s.thickness ? state.thickness.toFixed(2) : '';
  return [state.shape, extra, state.re, state.transition, alpha].join('|');
}

function panelAt(alpha) {
  if (!PANEL_SHAPES.includes(state.shape) || state.re < RE_PANEL_MIN) return null;
  const key = panelKey(alpha);
  if (panelCache.has(key)) return panelCache.get(key);
  let r = null;
  try {
    r = panelSolve(state.shape, {
      params: { naca: state.naca, thickness: state.thickness },
      alphaDeg: alpha,
      re: state.re,
      // Transition « turbulente » : déclenchement au bord d'attaque.
      xTrip: state.transition === 'turb' ? 0 : undefined,
    });
    if (r && !(Number.isFinite(r.cl) && Number.isFinite(r.cd))) r = null;
  } catch {
    r = null;
  }
  if (panelCache.size > 400) panelCache.clear();
  panelCache.set(key, r);
  return r;
}

/** Recalcule la référence et la zone laminaire imposée au calcul LBM. */
function refreshPanel() {
  panelRef = panelAt(alphaOf());
  if (!wt) return;
  const couple = wt.lbm.saActive && state.transition === 'free' && panelRef;
  wt.setTransition(couple ? panelRef : null);
  syncPanelRef();
}

function pct(v) {
  return Number.isFinite(v) ? `${Math.round(Math.min(1, Math.max(0, v)) * 100)} %` : '—';
}

/** Mach réel de l'essai (air seulement ; l'eau reste incompressible). */
function realMach() {
  return state.fluid === 'air' ? physical().V / 340 : 0;
}

function syncPanelRef() {
  const el = $('panel-ref');
  if (!panelRef) {
    el.hidden = true;
    return;
  }
  const r = panelRef;
  const k = sectionCoefficients(r);
  const mach = realMach();
  let cl = k.cl;
  let cm = k.cm;
  let comp = '';
  if (mach > 0.7) {
    comp = ` · Mach ${fmt(mach, 2)} : transsonique, hors du domaine de la méthode`;
  } else if (mach >= 0.1) {
    const c = compressibleCorrection(mach, Math.min(...r.cp.cp));
    cl *= c.factor;
    cm *= c.factor;
    comp = c.critical
      ? ` · Mach ${fmt(mach, 2)} : zone supersonique locale (Cp min < Cp* = ${fmt(c.cpStar, 2)}), correction de compressibilité hors domaine`
      : ` · corrigé de la compressibilité (Prandtl–Glauert, Mach ${fmt(mach, 2)})`;
  }
  const tr = state.transition === 'turb' ? 'transition au bord d’attaque' : `transition extrados ${pct(r.transitionUpper)}, intrados ${pct(r.transitionLower)}`;
  let txt = `Panneaux + couche limite : Cz ${fmt(cl)} · Cx ${fmt(k.cd, 4)} · Cm ${fmt(cm)} · ${tr}${comp}`;
  if (k.kirchhoff) txt += ` · décollement à ${pct(Math.min(r.separationUpper ?? 1, r.separationLower ?? 1))} de corde : décrochage approché (modèle de Kirchhoff)`;
  else if (state.re < 2e5 && state.transition === 'free') txt += ' · bulbes laminaires non modélisés';
  el.textContent = txt;
  el.classList.toggle('caution', k.kirchhoff || mach > 0.7);
  el.hidden = false;
}

// ---------------------------------------------------------------- distributions pariétales

let wallSampler = null;
let wallVisible = false;
let wallTick = 0;

/** Échantillonne Cp et Cf à la paroi (moyenne depuis la fin du transitoire). */
function sampleWall() {
  if (!wallVisible || diverged || !WALL_SHAPES.includes(state.shape)) return;
  const spt = wt.stepsPerTStar;
  const settle = tChange === 0 ? RAMP_T + SETTLE_T : SETTLE_T;
  if (wt.lbm.time < tChange + settle * spt) {
    wallSampler?.reset();
    return;
  }
  if (!wallSampler) wallSampler = new WallSampler(wt);
  if (!wallSampler.region) return;
  const lbm = wt.lbm;
  if (wt.engine === 'gpu') {
    // Lecture synchrone de la seule zone pariétale, une fois par seconde.
    if (++wallTick % 4 !== 0) return;
    const { x0, y0, w, h } = wallSampler.region;
    const f = lbm.readRegion(x0, y0, w, h);
    const col = lbm.readRegion(3, 0, 1, lbm.ny).rho;
    wallSampler.add(f, (x, y) => (y - y0) * w + (x - x0), upstreamDensity(col));
  } else {
    const nx = lbm.nx;
    const col = new Float32Array(lbm.ny);
    for (let y = 0; y < lbm.ny; y++) col[y] = lbm.rho[y * nx + 3];
    wallSampler.add(lbm, (x, y) => y * nx + x, upstreamDensity(col));
  }
}

function drawWallChart() {
  const card = $('card-wall');
  const ok = WALL_SHAPES.includes(state.shape);
  card.hidden = !ok;
  if (!ok || !wallVisible || !wt) return;
  const cyl = state.shape === 'cylinder';
  const cp = state.wallQ === 'cp';
  const res = wallSampler ? wallSampler.result() : { upper: [], lower: [], n: 0 };
  const key = (p) => (cp ? p.cp : 1000 * p.cf);
  const series = [];
  if (cyl) {
    if (cp) {
      const pts = [];
      for (let t = 0; t <= 180; t += 2) pts.push([t, 1 - 4 * Math.sin((t * Math.PI) / 180) ** 2]);
      series.push({ pts, color: themeColors.muted, style: 'dash' });
    }
  } else if (panelRef) {
    const r = panelRef;
    if (cp) {
      const { x, cp: c, iLE } = r.cp;
      const up = [];
      const lo = [];
      for (let i = 0; i < x.length; i++) {
        if (i >= iLE) up.push([x[i], c[i]]);
        if (i <= iLE) lo.unshift([x[i], c[i]]);
      }
      series.push({ pts: up, color: themeColors.lift, style: 'dash' }, { pts: lo, color: themeColors.accent, style: 'dash' });
    } else {
      const side = (f) => f.x.map((xv, i) => [xv, 1000 * f.cf[i] * f.ue[i] * f.ue[i]]);
      series.push(
        { pts: side(r.friction.upper), color: themeColors.lift, style: 'dash' },
        { pts: side(r.friction.lower), color: themeColors.accent, style: 'dash' },
      );
    }
  }
  series.push(
    { pts: res.upper.map((p) => [p.s, key(p)]), color: themeColors.lift, style: 'dots' },
    { pts: res.lower.map((p) => [p.s, key(p)]), color: themeColors.accent, style: 'dots' },
  );
  const marks = [];
  if (!cyl && panelRef && wt.lbm.saActive && state.transition === 'free') {
    if (Number.isFinite(panelRef.transitionUpper)) marks.push({ x: panelRef.transitionUpper, color: themeColors.lift });
    if (Number.isFinite(panelRef.transitionLower)) marks.push({ x: panelRef.transitionLower, color: themeColors.accent });
  }
  drawWall(
    $('wallchart'),
    {
      series,
      xLo: 0,
      xHi: cyl ? 180 : 1,
      xUnit: cyl ? '°' : '',
      yLabel: cp ? 'Cp' : 'Cf × 1000',
      invert: cp,
      marks,
      empty: 'Moyenne disponible après le transitoire',
    },
    themeColors,
  );
  $('key-up').textContent = cyl ? 'dessus' : 'extrados';
  $('key-lo').textContent = cyl ? 'dessous' : 'intrados';
  $('key-dash').textContent = cyl ? 'potentiel' : 'panneaux';
  const what = cyl ? 'angle depuis le point d’arrêt avant' : 'x/c';
  const lbmNote = res.n
    ? `Points : mailles voisines de la paroi, moyenne de ${res.n} relevés.`
    : 'Points : mailles voisines de la paroi (après le transitoire).';
  const refNote = cyl ? ' Tirets : écoulement potentiel (1 − 4 sin²θ).' : panelRef ? ' Tirets : méthode des panneaux + couche limite.' : '';
  const qNote = cp
    ? cyl
      ? ' Cp lu à moins d’une maille de la paroi ; l’écart au potentiel au-delà de ≈ 80° vient du décollement.'
      : ' Cp lu à moins d’une maille de la paroi : le pic d’aspiration du bord d’attaque est écrêté.'
    : ' Cf rapporté à ½ρV² ; négatif = écoulement de retour (décollement).' +
      (wt.lbm.saActive ? ' Frottement donné par la loi de paroi (laminaire en amont de la transition).' : '');
  $('wall-hint').textContent = `Abscisse : ${what}. ${lbmNote}${refNote}${qNote}`;
}

/** Polaire de référence, calculée par petits morceaux entre deux images. */
let refPolar = null;
function refPolarTick() {
  const s = shapeDef();
  if (!s.alpha || !PANEL_SHAPES.includes(state.shape) || state.re < RE_PANEL_MIN) {
    refPolar = null;
    return;
  }
  const key = panelKey('*');
  if (!refPolar || refPolar.key !== key) {
    const lo = Math.max(s.alpha[0], -12);
    const hi = Math.min(s.alpha[1], 20);
    const alphas = [];
    for (let a = Math.ceil(lo); a <= hi; a++) alphas.push(a);
    refPolar = { key, alphas, pts: [], next: 0 };
  }
  const t0 = performance.now();
  while (refPolar.next < refPolar.alphas.length && performance.now() - t0 < 20) {
    const a = refPolar.alphas[refPolar.next++];
    const r = panelAt(a);
    if (r) {
      const k = sectionCoefficients(r);
      refPolar.pts.push({ x: a, cx: k.cd, cz: k.cl, sep: false });
    }
  }
}

function sample() {
  const c = wt.averagedCoefficients();
  lastCoeffs = c;
  hist.push(wt.lbm.time, c.cx, c.cz, c.cm);
}

let halfSample = null;
/** Moyenne des efforts GPU relus aux pas n et n + 1. */
function pushGpuSample(smp) {
  if (halfSample && smp.time === halfSample.time + 1) {
    const c = wt.coefficientsFrom(
      0.5 * (smp.fx + halfSample.fx),
      0.5 * (smp.fy + halfSample.fy),
      0.5 * (smp.mz + halfSample.mz),
    );
    halfSample = null;
    lastCoeffs = c;
    hist.push(smp.time, c.cx, c.cz, c.cm);
  } else {
    halfSample = smp;
  }
}

function currentAnalysis() {
  const spt = wt.stepsPerTStar;
  const now = wt.lbm.time;
  // Au démarrage : mise en vitesse puis établissement du sillage.
  const settle = tChange === 0 ? RAMP_T + SETTLE_T : SETTLE_T;
  const start = Math.max(tChange + settle * spt, now - WINDOW_T * spt);
  const a = analyze(hist.since(start), wt.uIn, wt.L);
  a.tStart = start / spt;
  a.tEnd = now / spt;
  a.ready = now - start >= 3 * spt && a.n >= 20 && Number.isFinite(a.cx);
  return a;
}

function onDiverged() {
  diverged = true;
  running = false;
  stopSweep('Polaire interrompue : le calcul a divergé.');
  showToast(
    'Le calcul a divergé : Reynolds trop élevé pour cette résolution. Activez le modèle de sous-maille ou Spalart–Allmaras, réduisez Re ou affinez la grille, puis relancez.',
  );
  syncPlay();
}

// ---------------------------------------------------------------- thème

function themeName() {
  const t = document.documentElement.getAttribute('data-theme');
  if (t === 'dark' || t === 'light') return t;
  return window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
}

let themeColors = {};
function applyTheme() {
  const cs = getComputedStyle(document.documentElement);
  const v = (n) => cs.getPropertyValue(n).trim();
  themeColors = {
    ink: v('--ink'),
    muted: v('--muted'),
    grid: v('--grid'),
    accent: v('--accent'),
    lift: v('--lift'),
    body: v('--body-fill'),
    bodyEdge: v('--body-edge'),
    tracer: v('--tracer'),
    stage: v('--surface-2'),
    font: v('--font-num') || 'monospace',
  };
  cpuView.setTheme(themeName(), themeColors);
  if (gpuView) gpuView.setTheme(themeName(), themeColors);
  updateColorbar();
  drawCharts();
  wingTab?.draw();
}

// ---------------------------------------------------------------- onglets

let wingTab = null;
let lab3d = null;
const TAB_SUB = {
  '2d': 'CFD 2D Lattice Boltzmann · D2Q9',
  wing: 'Aile finie · surface portante + polaires 2D',
  '3d': 'CFD 3D Lattice Boltzmann · D3Q19 · WebGPU',
};

function setTab(tab) {
  state.tab = tab;
  for (const t of ['2d', 'wing', '3d']) {
    $(`tab-${t}`).hidden = t !== tab;
    $(`tabbtn-${t}`).setAttribute('aria-selected', String(t === tab));
  }
  $('app-sub').textContent = TAB_SUB[tab];
  if (tab === 'wing') wingTab.activate();
  else wingTab.deactivate();
  if (tab === '3d') lab3d.activate();
  else lab3d.deactivate();
  if (tab === '2d') {
    view.resize();
    drawCharts();
  }
  savePrefs();
}

// ---------------------------------------------------------------- formatage

function fmt(v, d = 3) {
  if (!Number.isFinite(v)) return '—';
  const s = Math.abs(v).toFixed(d);
  return (v < 0 && Number(s) !== 0 ? '−' : '') + s;
}

const DEC = (1.5).toLocaleString(navigator.language || 'fr-FR').includes(',') ? ',' : '.';
function fmtCsv(v, d = 4) {
  return Number.isFinite(v) ? v.toFixed(d).replace('.', DEC) : '';
}

// ---------------------------------------------------------------- interface : construction

function segmented(container, options, current, onPick) {
  container.textContent = '';
  for (const o of options) {
    const b = document.createElement('button');
    b.type = 'button';
    b.setAttribute('role', 'radio');
    b.dataset.value = String(o.value);
    b.textContent = o.label;
    b.setAttribute('aria-checked', String(o.value === current));
    b.addEventListener('click', () => {
      if (b.disabled) return;
      for (const c of container.children) c.setAttribute('aria-checked', String(c === b));
      onPick(o.value);
    });
    container.appendChild(b);
  }
}

function setSegValue(container, value) {
  for (const c of container.children) c.setAttribute('aria-checked', String(c.dataset.value === String(value)));
}

function buildShapeList() {
  const list = $('shape-list');
  list.textContent = '';
  for (const s of SHAPES) {
    const b = document.createElement('button');
    b.type = 'button';
    b.className = 'shape-btn';
    b.setAttribute('role', 'radio');
    b.dataset.shape = s.id;
    b.setAttribute('aria-checked', String(s.id === state.shape));
    b.setAttribute('aria-label', s.name);
    b.innerHTML = `<svg viewBox="0 0 24 24" aria-hidden="true"><path d="${shapeIconPath(s.id)}"/></svg><span>${s.short ?? s.name}</span>`;
    b.addEventListener('click', () => {
      if (sweep) return;
      state.shape = s.id;
      for (const c of list.children) c.setAttribute('aria-checked', String(c === b));
      syncShapeControls();
      pendingBody = true;
      savePrefs();
    });
    list.appendChild(b);
  }
}

function syncShapeControls() {
  const s = shapeDef();
  $('row-alpha').hidden = !s.alpha;
  $('row-naca').hidden = !s.naca;
  $('row-thick').hidden = !s.thickness;
  $('row-spin').hidden = !s.spin;
  if (s.alpha) {
    const a = $('alpha');
    a.min = s.alpha[0];
    a.max = s.alpha[1];
    a.value = alphaOf();
    $('alpha-out').textContent = `${fmt(alphaOf(), 1)}°`;
  }
  $('naca').value = state.naca;
  syncNacaHint();
  $('thick').value = state.thickness;
  $('thick-out').textContent = state.thickness.toFixed(2);
  $('spin').value = state.spin;
  $('spin-out').textContent = fmt(state.spin, 2);
  const L = RESOLUTIONS[state.res].L;
  $('ref-hint').textContent = `Longueur de référence : ${s.ref} (${L} mailles). Moments ${s.momentRef}.`;
  syncSweepLabels();
  drawCharts();
}

function syncNacaHint() {
  const p = parseNaca($('naca').value);
  const hint = $('naca-hint');
  const input = $('naca');
  if (!p) {
    input.setAttribute('aria-invalid', 'true');
    hint.textContent = 'Saisissez 4 chiffres, épaisseur 04 à 40 % (ex. 0012, 2412).';
    return;
  }
  input.removeAttribute('aria-invalid');
  const t = `épaisseur ${Math.round(p.t * 100)} %`;
  hint.textContent =
    p.m === 0
      ? `Symétrique, ${t}.`
      : `Cambrure ${Math.round(p.m * 100)} % à ${Math.round(p.p * 100)} % de corde, ${t}.`;
}

function syncSweepLabels() {
  const v = sweepVar();
  const unit = v === 'alpha' ? 'α' : 'ΩR/U';
  const labels = document.querySelectorAll('#sweep .sweep-row label');
  const words = [`${unit} de`, 'à', 'pas'];
  labels.forEach((l, i) => {
    l.firstChild.textContent = `${words[i]} `;
  });
  $('sweep').querySelector('h3').textContent =
    v === 'alpha' ? 'Polaire automatique en incidence' : 'Balayage automatique en rotation (Magnus)';
  if (!sweep) {
    if (v === 'alpha') {
      if ($('sw-from').dataset.var !== 'alpha') {
        $('sw-from').value = -4;
        $('sw-to').value = 12;
        $('sw-step').value = 2;
      }
    } else if ($('sw-from').dataset.var !== 'spin') {
      $('sw-from').value = 0;
      $('sw-to').value = 2;
      $('sw-step').value = 0.5;
    }
    $('sw-from').dataset.var = v;
  }
  setSegValue($('polar-seg'), state.polarMode);
  const first = $('polar-seg').children[0];
  if (first) first.textContent = v === 'alpha' ? 'Cx, Cz (α)' : 'Cx, Cz (ΩR/U)';
}

function updateColorbar() {
  const f = FIELDS[state.field];
  $('cb-label').textContent = f.label;
  $('cb-min').textContent = fmt(f.min, f.min % 1 ? 1 : 0);
  $('cb-max').textContent = fmt(f.max, f.max % 1 ? 1 : 0);
  $('cb-grad').style.background = gradientCss(themeName(), f.map);
}

function updateInfo() {
  if (!wt) return;
  const i = wt.info();
  const sa = !!wt.lbm.saActive;
  const s = shapeDef();
  const L = i.L;
  // Couche limite turbulente de plaque plane (ordre de grandeur) : épaisseur
  // δ ≈ 0,37 L Re^(−1/5), frottement Cf ≈ 0,026 Re^(−1/7), y⁺ du premier nœud
  // (½ maille de la paroi en moyenne).
  const re = state.re;
  const cf = 0.026 * re ** (-1 / 7);
  const yPlus = (0.5 / L) * re * Math.sqrt(cf / 2);
  const deltaCells = 0.37 * L * re ** -0.2;
  $('row-trans').hidden = !(sa && PANEL_SHAPES.includes(state.shape));
  $('numinfo').textContent =
    `${wt.engine === 'gpu' ? 'GPU' : 'CPU'} · collision ${state.collision === 'rr' ? 'régularisée' : 'BGK'} · ` +
    `Grille ${i.nx}×${i.ny} · L = ${L} mailles · τ = ${i.tau < 0.5001 ? `½ + ${sig3(i.tau - 0.5)}` : i.tau.toFixed(4)} · Ma = ${i.mach.toFixed(2)} · ` +
    `confinement ${(i.blockage * 100).toFixed(1)} % · sous-maille ${i.les ? 'active' : 'inactive'} · ` +
    (sa
      ? `Spalart–Allmaras + loi de paroi · y⁺ ≈ ${yPlus < 100 ? sig3(yPlus) : (Math.round(yPlus / 10) * 10).toLocaleString('fr-FR')} · ` +
        `couche limite ≈ ${sig3(deltaCells)} maille${deltaCells >= 2 ? 's' : ''}`
      : `Re fiable ≲ ${fmtRe(i.reLimit)}`);
  const warn = $('warn');
  const msgs = [];
  const profile = !['cylinder', 'square', 'triangle'].includes(s.id);
  if (sa) {
    msgs.push(
      `Turbulence modélisée (RANS Spalart–Allmaras, couche limite plus mince qu’une maille : frottement donné par la loi de paroi). ` +
        `Portance et décrochage sont qualitatifs ; la traînée LBM d’un profil mince reste approximative (de −50 % à +300 % à L = 32, ` +
        `davantage sur une grille plus grossière).` +
        (PANEL_SHAPES.includes(s.id) && panelRef ? ' Pour la traînée en écoulement attaché, fiez-vous à la référence « panneaux + couche limite ».' : ''),
    );
    if (!profile || (['plate', 'ellipse'].includes(s.id) && Math.abs(alphaOf()) >= 20)) {
      msgs.push(
        'Corps non profilé : la couche limite est supposée turbulente dès l’amont, comme au-delà de la crise de traînée ' +
          '(cylindre lisse : Re ≳ 3·10⁵, Cx ≈ 0,3 à 0,6). En dessous, préférez « Aucun » : le décollement laminaire domine.',
      );
    }
    if (L < 24) msgs.push(`Avec ${L} mailles par longueur, le bord d’attaque n’est pas résolu : passez en résolution Très fine ou au-delà.`);
  } else if (state.re > i.reLimit) {
    const next =
      state.res === 'rapide'
        ? 'passez en résolution Standard ou Fine'
        : state.res === 'standard'
          ? 'passez en résolution Fine'
          : state.res === 'fine'
            ? 'passez en résolution Très fine (sur ordinateur) ou réduisez Re'
            : 'réduisez Re';
    msgs.push(
      `Au-delà de Re ≈ ${fmtRe(i.reLimit)} avec ${L} mailles par longueur, la couche limite est trop peu résolue : ` +
        `la traînée est surestimée (de +20 % à +100 %) et les résultats sont qualitatifs ; ${next}` +
        (state.re >= 20000 ? ', ou activez Spalart–Allmaras (Modèle de turbulence).' : '.'),
    );
  }
  const a = Math.abs(alphaOf());
  if (!sa && state.re >= 10000 && profile) {
    msgs.push(
      'Vers Re ≈ 10⁴ apparaissent sur un profil réel des bulbes de décollement laminaire et la transition, ' +
        'phénomènes tridimensionnels que ce calcul 2D ne représente pas : résultats indicatifs, surtout en incidence.',
    );
  }
  const bluff = !profile || (['plate', 'ellipse'].includes(s.id) && a >= 20);
  if (!sa && bluff && state.re > RE_WAKE_3D) {
    msgs.push(
      `Au-delà de Re ≈ ${RE_WAKE_3D}, le sillage réel d’un corps non profilé devient tridimensionnel : ` +
        'le calcul 2D surestime Cx, St et l’amplitude de Cz (cylindre vers Re = 1000 : Cx ≈ 1,4 en 2D, ≈ 1 mesuré).',
    );
  }
  if (!i.les && i.tau < 0.515) msgs.push('τ proche de 0,5 sans modèle de sous-maille : risque de divergence.');
  warn.hidden = msgs.length === 0;
  warn.textContent = msgs.join(' ');
  $('ref-hint').textContent = `Longueur de référence : ${shapeDef().ref} (${i.L} mailles). Moments ${shapeDef().momentRef}.`;
}

// ---------------------------------------------------------------- interface : affichage

let lastUiUpdate = 0;

function statusOf(a) {
  if (diverged) return ['diverged', 'Divergence'];
  if (!running) return ['paused', 'En pause'];
  if (wt.ramping) return ['transient', 'Mise en vitesse'];
  if (!a || !a.ready) return ['transient', 'Transitoire'];
  if (a.periodic) {
    // Oscillation d'amplitude négligeable (bruit de bord de fuite, reste du transitoire).
    if (a.czAmp < 0.01 && a.cxAmp < 0.005) return ['steady', 'Quasi stationnaire'];
    return ['periodic', 'Périodique'];
  }
  const tol = 0.004 + 0.004 * Math.abs(a.cx);
  if (a.drift < tol) return ['steady', 'Stationnaire'];
  return ['transient', 'En convergence'];
}

function updateReadouts() {
  analysis = currentAnalysis();
  const a = analysis;
  const [st, label] = statusOf(a);
  const status = $('status');
  status.dataset.state = st;
  status.textContent = label;
  $('hud-time').textContent = `t* ${wt.tStar.toFixed(1)} · ${fmtDuration(wt.tStar * physical().secondsPerTStar)}`;

  const useMean = a.ready;
  const cx = useMean ? a.cx : lastCoeffs.cx;
  const cz = useMean ? a.cz : lastCoeffs.cz;
  const cm = useMean ? a.cm : lastCoeffs.cm;
  $('ro-cx').textContent = fmt(cx);
  $('ro-cz').textContent = fmt(cz);
  $('ro-cm').textContent = fmt(cm);
  if (useMean) {
    $('ro-cx-amp').textContent = a.periodic && a.cxAmp > 0.0005 ? `± ${fmt(a.cxAmp)}` : 'moyenne';
    $('ro-cz-amp').textContent = a.periodic && a.czAmp > 0.0005 ? `± ${fmt(a.czAmp)}` : 'moyenne';
  } else {
    $('ro-cx-amp').textContent = 'instantané';
    $('ro-cz-amp').textContent = 'instantané';
  }
  $('ro-cm-sub').textContent = shapeDef().momentRef.replace('au ', '').replace('quart de corde', 'c/4');
  if (useMean && a.periodic && Number.isFinite(a.st)) {
    $('ro-st').textContent = fmt(a.st, 3);
    $('ro-st-sub').textContent = `${a.periods} périodes`;
  } else {
    $('ro-st').textContent = '—';
    $('ro-st-sub').textContent = useMean ? 'permanent' : '';
  }
}

let perfSteps = 0;
let perfT0 = performance.now();
let perfText = '';
let achievedSpeed = NaN;
function perfRate() {
  const now = performance.now();
  if (now - perfT0 > 1000) {
    const sps = (perfSteps * 1000) / (now - perfT0);
    const mlups = (sps * wt.nx * wt.ny) / 1e6;
    const tps = (sps * wt.uIn) / wt.L;
    achievedSpeed = running && !diverged ? sps * physical().secondsPerStep : NaN;
    perfText = `${Math.round(sps)} pas/s · ${mlups.toFixed(1)} MLUPS · ${tps.toFixed(2)} t*/s`;
    perfSteps = 0;
    perfT0 = now;
  }
  return perfText;
}

// ---------------------------------------------------------------- échelle physique

/** Grandeurs dimensionnelles déduites de Re, de L (m) et du fluide. */
function physical() {
  const f = FLUIDS[state.fluid];
  const V = (state.re * f.nu) / state.chord; // m/s
  const secondsPerTStar = state.chord / V;
  return {
    V,
    rho: f.rho,
    secondsPerTStar,
    secondsPerStep: wt ? secondsPerTStar / wt.stepsPerTStar : NaN,
  };
}

/** Nombre avec 3 chiffres significatifs, sans notation exponentielle si possible. */
function sig3(v) {
  if (!Number.isFinite(v)) return '—';
  const a = Math.abs(v);
  if (a >= 1000 || (a > 0 && a < 1e-3)) return v.toExponential(2).replace('e', '·10^').replace('+', '');
  return Number(v.toPrecision(3)).toString();
}

function fmtDuration(sec) {
  if (!Number.isFinite(sec)) return '—';
  if (sec < 60) return `${sig3(sec)} s`;
  if (sec < 3600) return `${Math.floor(sec / 60)} min ${String(Math.round(sec % 60)).padStart(2, '0')} s`;
  return `${(sec / 3600).toFixed(1)} h`;
}

function speedLabel(v) {
  if (v === 'max') return 'Max';
  if (v === 1) return 'Temps réel';
  return `×${v}`;
}

/** Infos physiques (vitesse, durées, efforts par mètre d'envergure). */
function updatePhysical() {
  const p = physical();
  const parts = [
    `V = ${sig3(p.V)} m/s`,
    `1 t* = L/V = ${fmtDuration(p.secondsPerTStar)}`,
    `1 pas = ${sig3(p.secondsPerStep)} s`,
  ];
  const a = analysis;
  if (a && a.ready) {
    const q = 0.5 * p.rho * p.V * p.V * state.chord; // N/m par unité de coefficient
    parts.push(`traînée ${sig3(a.cx * q)} N/m, portance ${sig3(a.cz * q)} N/m d\u2019envergure`);
  }
  // Écoulement réel compressible (air, Mach > 0,3) : le calcul reste incompressible.
  const mach = state.fluid === 'air' ? p.V / 340 : 0;
  if (mach > 0.3) parts.push(`Mach réel ≈ ${mach >= 10 ? Math.round(mach) : sig3(mach)} : compressibilité non représentée`);
  $('phys-info').textContent = parts.join(' · ');
  const hint = $('speed-hint');
  if (sweep) {
    hint.textContent = 'Polaire en cours : calcul à vitesse maximale.';
  } else if (state.speed === 'max') {
    hint.textContent = Number.isFinite(achievedSpeed)
      ? `Calcul aussi rapide que possible : ×${sig3(achievedSpeed)} par rapport au temps réel.`
      : 'Calcul aussi rapide que possible.';
  } else {
    const ok = Number.isFinite(achievedSpeed) && achievedSpeed >= 0.95 * state.speed;
    hint.textContent = Number.isFinite(achievedSpeed)
      ? `Cible ${speedLabel(state.speed).toLowerCase()}, atteint ×${sig3(achievedSpeed)}${ok ? '.' : ' : le calcul ne suit pas, réduisez la résolution ou Re, ou augmentez L.'}`
      : `Cible ${speedLabel(state.speed).toLowerCase()}.`;
  }
}

function drawCharts() {
  if (!wt || !hist) return;
  drawWallChart();
  const spt = wt.stepsPerTStar;
  const tMax = Math.max(wt.tStar, 10);
  const tMin = Math.max(0, tMax - HISTORY_SPAN_T);
  const w = hist.since(tMin * spt);
  const t = new Float64Array(w.t.length);
  for (let k = 0; k < t.length; k++) t[k] = w.t[k] / spt;
  drawHistory(
    $('chart'),
    {
      t,
      series: [
        { values: w.cx, color: themeColors.accent },
        { values: w.cz, color: themeColors.lift },
      ],
      tMin,
      tMax,
      markers,
      meanWindow: analysis && analysis.ready ? [analysis.tStart, analysis.tEnd] : null,
    },
    themeColors,
  );
  const key = configKey();
  const v = sweepVar();
  const pts = points
    .filter((p) => p.key === key)
    .map((p) => ({ x: v === 'alpha' ? p.alpha : p.spin, cx: p.cx, cz: p.cz }));
  const ref = v === 'alpha' && refPolar && refPolar.pts.length > 1 ? refPolar.pts : null;
  drawPolar(
    $('polar'),
    pts,
    state.polarMode,
    themeColors,
    { lift: themeColors.lift, drag: themeColors.accent },
    v === 'alpha' ? '°' : '',
    ref,
  );
  const refTxt = ref ? ' Tirets : méthode des panneaux + couche limite (écoulement attaché).' : '';
  $('polar-hint').textContent = pts.length
    ? `${pts.length} point${pts.length > 1 ? 's' : ''} pour ${wt.body.label} à Re ${fmtRe(state.re)} (résolution ${RESOLUTIONS[state.res].label.toLowerCase()}).${refTxt}`
    : `La polaire trace les points enregistrés pour la configuration courante (forme, Re, résolution).${refTxt}`;
}

function renderTable() {
  const tbody = $('results').querySelector('tbody');
  tbody.textContent = '';
  if (!points.length) {
    const tr = document.createElement('tr');
    tr.innerHTML = '<td class="empty" colspan="8">Aucun point enregistré pour l’instant.</td>';
    tbody.appendChild(tr);
    return;
  }
  for (const p of [...points].reverse()) {
    const tr = document.createElement('tr');
    const cells = [
      p.label,
      p.alpha === null ? '—' : `${fmt(p.alpha, 1)}°`,
      fmtRe(p.re) + (p.turb && p.turb !== 'aucun' ? ' SA' : ''),
      fmt(p.cx),
      fmt(p.cz),
      fmt(p.cm),
      Number.isFinite(p.st) ? fmt(p.st) : '—',
      Math.abs(p.cx) > 1e-6 ? fmt(p.cz / p.cx, 2) : '—',
    ];
    for (const c of cells) {
      const td = document.createElement('td');
      td.textContent = c;
      tr.appendChild(td);
    }
    tbody.appendChild(tr);
  }
}

function configKey() {
  const s = shapeDef();
  const extra = s.naca ? state.naca : s.thickness ? state.thickness.toFixed(2) : '';
  const spin = s.spin && sweepVar() !== 'spin' ? state.spin : '';
  const turb = turbMode() === 'sa' ? `sa-${state.transition}` : '';
  return [state.shape, extra, spin, state.re, RESOLUTIONS[state.res].L, state.u, state.wall, turb].join('|');
}

let toastTimer = 0;
function showToast(msg, ms = 0) {
  const t = $('toast');
  t.textContent = msg;
  t.hidden = false;
  clearTimeout(toastTimer);
  if (ms) toastTimer = setTimeout(hideToast, ms);
}

function hideToast() {
  $('toast').hidden = true;
}

function syncPlay() {
  const btn = $('btn-play');
  btn.setAttribute('aria-label', running ? 'Mettre en pause' : 'Reprendre le calcul');
  btn.querySelector('.ic-pause').toggleAttribute('hidden', !running);
  btn.querySelector('.ic-play').toggleAttribute('hidden', running);
}

// ---------------------------------------------------------------- mesures

function recordPoint(a, source = 'manuel') {
  const s = shapeDef();
  const p = {
    key: configKey(),
    label: wt.body.label,
    shape: state.shape,
    alpha: s.alpha ? alphaOf() : null,
    spin: s.spin ? state.spin : null,
    re: state.re,
    L: wt.L,
    cx: a.cx,
    cz: a.cz,
    cm: a.cm,
    st: a.periodic ? a.st : NaN,
    turb: turbMode() === 'sa' ? (state.transition === 'free' ? 'SA, transition libre' : 'SA, turbulent') : 'aucun',
    source,
  };
  points.push(p);
  if (points.length > 500) points.shift();
  storageSet(POINTS_KEY, points);
  renderTable();
  drawCharts();
}

function copyTable() {
  const head = ['maquette', 'alpha_deg', 'rotation_OmegaR_U', 'Re', 'L_mailles', 'turbulence', 'Cx', 'Cz', 'Cm', 'St', 'Cz_Cx'];
  const rows = points.map((p) =>
    [
      p.label,
      p.alpha === null ? '' : fmtCsv(p.alpha, 1),
      p.spin === null ? '' : fmtCsv(p.spin, 2),
      p.re,
      p.L,
      p.turb ?? 'aucun',
      fmtCsv(p.cx),
      fmtCsv(p.cz),
      fmtCsv(p.cm),
      fmtCsv(p.st),
      Math.abs(p.cx) > 1e-6 ? fmtCsv(p.cz / p.cx, 3) : '',
    ].join('\t'),
  );
  const text = [head.join('\t'), ...rows].join('\n');
  const fallback = () => {
    const ta = $('csv-fallback');
    ta.hidden = false;
    ta.value = text;
    ta.focus();
    ta.select();
    showToast('Le tableau est sélectionné : utilisez Copier.', 2500);
  };
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard.writeText(text).then(
      () => {
        $('btn-copy').textContent = 'Copié';
        setTimeout(() => ($('btn-copy').textContent = 'Copier le tableau'), 1600);
      },
      fallback,
    );
  } else {
    fallback();
  }
}

let clearArmed = 0;
function clearPoints() {
  const btn = $('btn-clear');
  if (!clearArmed) {
    btn.textContent = 'Confirmer l’effacement';
    btn.classList.add('danger');
    clearArmed = setTimeout(() => {
      clearArmed = 0;
      btn.textContent = 'Effacer';
      btn.classList.remove('danger');
    }, 3000);
    return;
  }
  clearTimeout(clearArmed);
  clearArmed = 0;
  btn.textContent = 'Effacer';
  btn.classList.remove('danger');
  points = [];
  storageSet(POINTS_KEY, points);
  renderTable();
  drawCharts();
}

// ---------------------------------------------------------------- balayage automatique

function setVar(v) {
  if (sweep.var === 'alpha') {
    state.alphaByShape[state.shape] = v;
    $('alpha').value = v;
    $('alpha-out').textContent = `${fmt(v, 1)}°`;
  } else {
    state.spin = v;
    $('spin').value = v;
    $('spin-out').textContent = fmt(v, 2);
  }
  applyBody();
}

function lockControls(locked) {
  const ids = ['alpha', 'naca', 'thick', 'spin', 're'];
  for (const id of ids) $(id).disabled = locked;
  document
    .querySelectorAll(
      '.shape-btn, .step, .chip-btn, #res-seg button, #u-seg button, #wall-seg button, #les-seg button, #engine-seg button, #coll-seg button, #btn-reset',
    )
    .forEach((b) => (b.disabled = locked));
}

function startSweep() {
  if (sweep) {
    stopSweep('Polaire interrompue.');
    return;
  }
  const from = Number($('sw-from').value);
  const to = Number($('sw-to').value);
  const step = Math.abs(Number($('sw-step').value));
  if (![from, to, step].every(Number.isFinite) || step <= 0) {
    showToast('Bornes de balayage invalides.', 2500);
    return;
  }
  const v = sweepVar();
  const lim = v === 'alpha' ? shapeDef().alpha : [-3, 3];
  const values = [];
  const dir = to >= from ? 1 : -1;
  for (let x = from; dir * (to - x) >= -1e-9 && values.length < 60; x += dir * step) {
    values.push(Math.round(Math.min(lim[1], Math.max(lim[0], x)) * 100) / 100);
  }
  const uniq = values.filter((x, i) => values.indexOf(x) === i);
  if (!uniq.length) return;
  if (!running && !diverged) {
    running = true;
    syncPlay();
  }
  sweep = { var: v, values: uniq, i: 0, phase: 'settle', t0: wt.lbm.time };
  setVar(uniq[0]);
  sweep.t0 = wt.lbm.time;
  lockControls(true);
  $('btn-sweep').textContent = 'Arrêter';
  $('sweep-progress').hidden = false;
  if (navigator.wakeLock && navigator.wakeLock.request) {
    navigator.wakeLock.request('screen').then(
      (l) => (wakeLock = l),
      () => {},
    );
  }
  updateSweepUi();
}

function stopSweep(msg) {
  if (!sweep) return;
  sweep = null;
  lockControls(false);
  $('btn-sweep').textContent = 'Lancer';
  $('sweep-text').textContent = msg || 'Polaire terminée.';
  $('sweep-bar').style.width = '100%';
  if (wakeLock) {
    wakeLock.release().catch(() => {});
    wakeLock = null;
  }
}

function sweepTick() {
  if (!sweep || !running) return;
  const spt = wt.stepsPerTStar;
  const el = (wt.lbm.time - sweep.t0) / spt;
  // Critère adaptatif : un point stationnaire convergé n'attend pas la durée maximale.
  const converged = (a) => !a.periodic && a.drift < 0.003 + 0.003 * Math.abs(a.cx);
  if (sweep.phase === 'settle') {
    const recent = el >= SWEEP_SETTLE_MIN_T ? analyze(hist.since(wt.lbm.time - 4 * spt), wt.uIn, wt.L) : null;
    if (el >= SWEEP_SETTLE_T || (recent && converged(recent))) {
      sweep.phase = 'measure';
      sweep.t0 = wt.lbm.time;
    }
  } else {
    const a = analyze(hist.since(sweep.t0), wt.uIn, wt.L);
    const enough =
      (el >= SWEEP_MEASURE_MIN_T && converged(a)) ||
      (el >= SWEEP_MEASURE_MIN_T && a.periodic && a.periods >= 3) ||
      (el >= SWEEP_MEASURE_T && (!a.periodic || a.periods >= 3)) ||
      el >= 2.5 * SWEEP_MEASURE_T;
    if (enough && Number.isFinite(a.cx)) {
      recordPoint(a, 'polaire');
      sweep.i++;
      if (sweep.i >= sweep.values.length) {
        stopSweep(`Polaire terminée : ${sweep.values.length} points enregistrés.`);
        return;
      }
      setVar(sweep.values[sweep.i]);
      sweep.phase = 'settle';
      sweep.t0 = wt.lbm.time;
    }
  }
  updateSweepUi();
}

function updateSweepUi() {
  if (!sweep) return;
  const n = sweep.values.length;
  const spt = wt.stepsPerTStar;
  const el = (wt.lbm.time - sweep.t0) / spt;
  const frac = sweep.phase === 'settle' ? 0.45 * Math.min(1, el / SWEEP_SETTLE_T) : 0.45 + 0.55 * Math.min(1, el / SWEEP_MEASURE_T);
  $('sweep-bar').style.width = `${(((sweep.i + frac) / n) * 100).toFixed(1)}%`;
  const v = sweep.values[sweep.i];
  const name = sweep.var === 'alpha' ? `α = ${fmt(v, 1)}°` : `ΩR/U = ${fmt(v, 2)}`;
  $('sweep-text').textContent = `Point ${sweep.i + 1}/${n} · ${name} · ${sweep.phase === 'settle' ? 'stabilisation' : 'mesure'}`;
}

// ---------------------------------------------------------------- boucle principale

let lastFrameTs = 0;
let frameMs = 16;

let stepBudget = 0;

/** Pas autorisés pour cette image quand une vitesse de simulation est imposée. */
function allowedSteps(dtMs) {
  if (state.speed === 'max' || sweep) {
    stepBudget = 0;
    return Infinity;
  }
  const perSecond = state.speed / physical().secondsPerStep;
  // Au plus 100 ms de retard accumulé : on ne rattrape pas une longue pause.
  stepBudget = Math.min(stepBudget + (perSecond * Math.min(dtMs, 100)) / 1000, perSecond * 0.1 + 1);
  return Math.floor(stepBudget);
}

function frame(ts) {
  if (state.tab !== '2d') {
    // Soufflerie 2D en sommeil pendant les autres onglets.
    lastFrameTs = 0;
    requestAnimationFrame(frame);
    return;
  }
  const dtFrame = lastFrameTs ? ts - lastFrameTs : 16;
  if (lastFrameTs) frameMs = 0.8 * frameMs + 0.2 * Math.min(200, ts - lastFrameTs);
  lastFrameTs = ts;
  if (pendingBody && !diverged) {
    pendingBody = false;
    applyBody();
  }
  if (running && !diverged) {
    const t0 = performance.now();
    let n = 0;
    const lbm = wt.lbm;
    const allowed = allowedSteps(dtFrame);
    if (wt.engine === 'gpu') {
      lbm.beginFrame();
      const todo = Math.min(gpuSteps, allowed);
      for (; n < todo; n++) {
        wt.step();
        // Deux pas consécutifs par échantillon : la moyenne annule une
        // éventuelle alternance pair/impair des efforts.
        const r = lbm.time % sampleEvery;
        if (r === 0 || r === 1) lbm.sampleForces();
      }
      for (const smp of lbm.readSamples()) pushGpuSample(smp);
      // Réglage du nombre de pas par image : le GPU doit rester occupé (une
      // relecture en attente au plus) sans faire chuter l'affichage sous
      // ~30 images/s. Relecture asynchrone : gpuLag = relectures en vol.
      // Une image lente alors que le GPU a fini (lag ≤ 1) vient du rendu ou
      // de la page, pas du calcul : on n'y réduit pas le nombre de pas.
      const lag = lbm.gpuLag;
      if (lag >= 3) gpuSteps = Math.max(1, Math.floor(gpuSteps * 0.75));
      else if (lag === 2 && frameMs > 40) gpuSteps = Math.max(1, Math.floor(gpuSteps * 0.85));
      else if (lag <= 1 && frameMs < 60 && n === gpuSteps) gpuSteps = Math.min(4000, Math.ceil(gpuSteps * 1.15) + 1);
    } else {
      while (n < allowed && (n === 0 || performance.now() - t0 < FRAME_BUDGET_MS)) {
        wt.step();
        n++;
        if (lbm.time % sampleEvery === 0) sample();
      }
    }
    if (allowed !== Infinity) stepBudget -= n;
    stepsSinceRender += n;
    perfSteps += n;
    if (!lbm.isFinite()) onDiverged();
  }

  const ok = diverged || wt.lbm.computeMacros();
  if (!ok) onDiverged();
  if (!diverged) {
    if (state.tracers) view.advect(stepsSinceRender);
    stepsSinceRender = 0;
  }
  const a = analysis;
  const force = state.arrow && a && a.ready ? { cx: a.cx, cz: a.cz } : null;
  view.draw({ field: state.field, tracers: state.tracers, force });

  const now = performance.now();
  if (now - lastUiUpdate > 250) {
    lastUiUpdate = now;
    updateReadouts();
    drawCharts();
    sweepTick();
    refPolarTick();
    sampleWall();
    $('perf').textContent = perfRate() || '—';
    updatePhysical();
  }
  requestAnimationFrame(frame);
}

// ---------------------------------------------------------------- évènements

function bind() {
  $('btn-play').addEventListener('click', () => {
    if (diverged) {
      buildSim();
      running = true;
    } else {
      running = !running;
    }
    syncPlay();
  });
  $('btn-reset').addEventListener('click', () => {
    buildSim();
    running = true;
    syncPlay();
  });

  const alpha = $('alpha');
  alpha.addEventListener('input', () => {
    const v = Number(alpha.value);
    state.alphaByShape[state.shape] = v;
    $('alpha-out').textContent = `${fmt(v, 1)}°`;
    pendingBody = true;
    savePrefs();
  });
  document.querySelectorAll('.step').forEach((b) =>
    b.addEventListener('click', () => {
      const input = $(b.dataset.target);
      const v = Math.min(Number(input.max), Math.max(Number(input.min), Number(input.value) + Number(b.dataset.delta)));
      input.value = v;
      input.dispatchEvent(new Event('input'));
    }),
  );

  const naca = $('naca');
  naca.addEventListener('input', () => {
    naca.value = naca.value.replace(/\D/g, '').slice(0, 4);
    syncNacaHint();
    if (parseNaca(naca.value)) {
      state.naca = naca.value;
      pendingBody = true;
      savePrefs();
    }
  });
  document.querySelectorAll('[data-naca]').forEach((b) =>
    b.addEventListener('click', () => {
      naca.value = b.dataset.naca;
      naca.dispatchEvent(new Event('input'));
    }),
  );

  $('thick').addEventListener('input', () => {
    state.thickness = Number($('thick').value);
    $('thick-out').textContent = state.thickness.toFixed(2);
    pendingBody = true;
    savePrefs();
  });
  $('spin').addEventListener('input', () => {
    state.spin = Number($('spin').value);
    $('spin-out').textContent = fmt(state.spin, 2);
    pendingBody = true;
    savePrefs();
  });

  const re = $('re');
  re.max = RE_VALUES.length - 1;
  re.value = RE_VALUES.indexOf(state.re);
  $('re-out').textContent = fmtRe(state.re);
  re.addEventListener('input', () => {
    state.re = RE_VALUES[Number(re.value)];
    $('re-out').textContent = fmtRe(state.re);
    wt.setReynolds(state.re);
    syncTurbulence();
    refreshPanel();
    markChange();
    updateInfo();
    drawCharts();
    savePrefs();
  });
  segmented(
    $('turb-seg'),
    [
      { value: 'auto', label: 'Auto' },
      { value: 'none', label: 'Aucun' },
      { value: 'sa', label: 'Spalart–Allmaras' },
    ],
    state.turb,
    (v) => {
      state.turb = v;
      syncTurbulence();
      refreshPanel();
      markChange();
      updateInfo();
      drawCharts();
      savePrefs();
    },
  );
  segmented(
    $('trans-seg'),
    [
      { value: 'free', label: 'Libre (eᴺ)' },
      { value: 'turb', label: 'Bord d’attaque' },
    ],
    state.transition,
    (v) => {
      state.transition = v;
      refreshPanel();
      markChange();
      updateInfo();
      drawCharts();
      savePrefs();
    },
  );

  segmented(
    $('speed-seg'),
    SPEEDS.map((v) => ({ value: v, label: speedLabel(v) })),
    state.speed,
    (v) => {
      state.speed = v;
      stepBudget = 0;
      updatePhysical();
      savePrefs();
    },
  );
  segmented(
    $('fluid-seg'),
    Object.entries(FLUIDS).map(([k, f]) => ({ value: k, label: f.label })),
    state.fluid,
    (v) => {
      state.fluid = v;
      updatePhysical();
      syncPanelRef();
      savePrefs();
    },
  );
  const chord = $('chord');
  chord.value = state.chord;
  chord.addEventListener('input', () => {
    const v = Number(chord.value);
    if (v > 0 && v <= 100) {
      state.chord = v;
      chord.removeAttribute('aria-invalid');
      updatePhysical();
      syncPanelRef();
      savePrefs();
    } else {
      chord.setAttribute('aria-invalid', 'true');
    }
  });

  buildResSeg();
  segmented(
    $('engine-seg'),
    [
      { value: 'gpu', label: 'GPU (WebGL2)' },
      { value: 'cpu', label: 'CPU' },
    ],
    engine(),
    (v) => {
      state.engine = v;
      buildResSeg();
      buildSim();
      syncShapeControls();
      savePrefs();
    },
  );
  if (!gpuOk) {
    const b = $('engine-seg').querySelector('[data-value="gpu"]');
    b.disabled = true;
    b.title = 'WebGL2 avec rendu flottant indisponible sur cet appareil';
    $('engine-hint').textContent = 'Ce navigateur n\u2019offre pas WebGL2 avec rendu flottant : calcul sur le processeur.';
  }
  segmented(
    $('coll-seg'),
    [
      { value: 'rr', label: 'Régularisée' },
      { value: 'bgk', label: 'BGK' },
    ],
    state.collision,
    (v) => {
      state.collision = v;
      wt.lbm.collision = v;
      markChange();
      updateInfo();
      savePrefs();
    },
  );
  segmented(
    $('field-seg'),
    [
      { value: 'vort', label: 'Vorticité' },
      { value: 'speed', label: 'Vitesse' },
      { value: 'cp', label: 'Pression' },
    ],
    state.field,
    (v) => {
      state.field = v;
      updateColorbar();
      savePrefs();
    },
  );
  segmented(
    $('view-seg'),
    Object.entries(VIEWS).map(([k, v]) => ({ value: k, label: v.label })),
    state.view,
    (v) => {
      state.view = v;
      cpuView.setView(v);
      if (gpuView) gpuView.setView(v);
      savePrefs();
    },
  );
  segmented(
    $('wall-seg'),
    [
      { value: 'bouzidi', label: 'Interpolée (Bouzidi)' },
      { value: 'bb', label: 'Marches' },
    ],
    state.wall,
    (v) => {
      state.wall = v;
      wt.lbm.bouzidi = v === 'bouzidi';
      markChange();
      drawCharts();
      savePrefs();
    },
  );
  segmented(
    $('les-seg'),
    [
      { value: 'auto', label: 'Auto' },
      { value: 'on', label: 'Activé' },
      { value: 'off', label: 'Désactivé' },
    ],
    state.les,
    (v) => {
      state.les = v;
      wt.lbm.setLes(v);
      markChange();
      updateInfo();
      savePrefs();
    },
  );
  segmented(
    $('u-seg'),
    U_VALUES.map((u) => ({ value: u, label: `${u.toFixed(2)} · Ma ${(u * Math.sqrt(3)).toFixed(2)}` })),
    state.u,
    (v) => {
      state.u = v;
      buildSim();
      savePrefs();
    },
  );
  segmented(
    $('wallq-seg'),
    [
      { value: 'cp', label: 'Cp' },
      { value: 'cf', label: 'Cf' },
    ],
    state.wallQ,
    (v) => {
      state.wallQ = v;
      drawWallChart();
      savePrefs();
    },
  );
  if ('IntersectionObserver' in window) {
    new IntersectionObserver((entries) => {
      wallVisible = entries.some((e) => e.isIntersecting);
      if (wallVisible) drawWallChart();
    }).observe($('card-wall'));
  } else {
    wallVisible = true;
  }
  segmented(
    $('polar-seg'),
    [
      { value: 'var', label: 'Cx, Cz (α)' },
      { value: 'eiffel', label: 'Polaire Cz(Cx)' },
    ],
    state.polarMode,
    (v) => {
      state.polarMode = v;
      drawCharts();
      savePrefs();
    },
  );

  $('tracers').checked = state.tracers;
  $('tracers').addEventListener('change', () => {
    state.tracers = $('tracers').checked;
    if (!state.tracers) view.clearParticles();
    savePrefs();
  });
  $('arrow').checked = state.arrow;
  $('arrow').addEventListener('change', () => {
    state.arrow = $('arrow').checked;
    savePrefs();
  });

  $('btn-record').addEventListener('click', () => {
    const a = currentAnalysis();
    if (!a.ready) {
      showToast('Pas encore de moyenne : laissez l’écoulement s’établir quelques t*.', 2500);
      return;
    }
    recordPoint(a);
    $('btn-record').textContent = 'Point enregistré';
    setTimeout(() => ($('btn-record').textContent = 'Enregistrer le point'), 1400);
  });
  $('btn-copy').addEventListener('click', copyTable);
  $('btn-clear').addEventListener('click', clearPoints);
  $('btn-sweep').addEventListener('click', startSweep);
  $('toast').addEventListener('click', hideToast);

  const ro = new ResizeObserver(() => {
    view.resize();
    drawCharts();
  });
  ro.observe($('flow'));
  ro.observe($('chart'));
  ro.observe($('wallchart'));

  if (window.matchMedia) {
    const mq = window.matchMedia('(prefers-color-scheme: dark)');
    mq.addEventListener?.('change', applyTheme);
  }
  new MutationObserver(applyTheme).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
}

/** Applique le modèle de turbulence effectif (Auto suit le Reynolds). */
function syncTurbulence() {
  const mode = turbMode();
  if ((mode === 'sa') !== !!wt.lbm.saActive) wt.lbm.setTurbulence(mode);
}

// ---------------------------------------------------------------- démarrage

/** Résolutions proposées selon le moteur (48 et 64 mailles : GPU seulement). */
function buildResSeg() {
  const opts = Object.entries(RESOLUTIONS)
    .filter(([k]) => resAllowed(k))
    .map(([k, r]) => ({ value: k, label: `${r.label} · ${r.L}` }));
  if (!resAllowed(state.res)) state.res = 'tresfine';
  segmented($('res-seg'), opts, state.res, (v) => {
    state.res = v;
    buildSim();
    syncShapeControls();
    savePrefs();
  });
  setSegValue($('engine-seg'), engine());
}

function start() {
  cpuView = new FlowView($('flow'));
  cpuView.setView(state.view);
  initGpu();
  if (gpuView) gpuView.setView(state.view);
  // Première visite avec GPU : résolution plus fine par défaut.
  if (!storageGet(PREFS_KEY) && gpuOk) state.res = 'tresfine';
  view = cpuView;
  $('flow-gl').addEventListener('webglcontextlost', (e) => {
    e.preventDefault();
    gpuOk = false;
    wt = null;
    buildResSeg();
    buildSim();
    syncShapeControls();
    showToast('Le GPU a été réinitialisé par le navigateur : le calcul continue sur le processeur.', 5000);
  });
  // Ligne de performance sous les coefficients (créée ici pour rester hors du HTML statique).
  const perf = document.createElement('p');
  perf.className = 'numinfo num';
  perf.id = 'perf';
  perf.textContent = '—';
  $('numinfo').after(perf);

  buildShapeList();
  bind();
  wingTab = new WingTab({
    segmented,
    theme: () => themeColors,
    fmt,
    fmtRe,
    sig3,
    fluids: FLUIDS,
    storageGet,
    storageSet,
  });
  wingTab.init();
  lab3d = new Lab3D({ segmented, theme: () => themeColors, fmt, sig3, storageGet, storageSet });
  lab3d.init();
  for (const t of ['2d', 'wing', '3d']) $(`tabbtn-${t}`).addEventListener('click', () => setTab(t));
  buildSim();
  syncShapeControls();
  applyTheme();
  renderTable();
  syncPlay();
  setTab(state.tab);
  requestAnimationFrame(frame);

  if ('serviceWorker' in navigator && !window.__SOUFFLERIE_INLINE__) {
    const local = location.hostname === 'localhost' || location.hostname === '127.0.0.1';
    if (location.protocol === 'https:' || local) navigator.serviceWorker.register('./sw.js').catch(() => {});
  }
}

start();
