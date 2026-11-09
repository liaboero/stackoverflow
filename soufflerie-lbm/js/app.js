/**
 * Soufflerie de poche — interface mobile.
 */
import { WindTunnel, RESOLUTIONS, RAMP_T, RE_WAKE_3D } from './sim.js';
import { SHAPES, SHAPE_BY_ID, parseNaca, shapeIconPath } from './shapes.js';
import { ForceHistory, analyze } from './stats.js';
import { FlowView, FIELDS, VIEWS, gradientCss } from './render.js';
import { drawHistory, drawPolar } from './chart.js';
import { GpuLBM, gpuSupported } from './lbm-gpu.js';
import { GpuFlowView } from './render-gpu.js';

const RE_VALUES = [10, 20, 40, 60, 100, 150, 200, 300, 500, 1000, 1500, 2000, 3000, 5000, 10000, 15000, 20000];
const U_VALUES = [0.05, 0.08, 0.1, 0.12];
const WINDOW_T = 20; // fenêtre glissante de moyenne (t*)
const SETTLE_T = 6; // transitoire ignoré après un changement (t*)
const HISTORY_SPAN_T = 40; // largeur du graphique (t*)
const SWEEP_SETTLE_T = 12;
const SWEEP_MEASURE_T = 15;
const FRAME_BUDGET_MS = 11;
const GPU_FRAME_MS = 14; // temps de calcul GPU visé par image
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
if (!parseNaca(state.naca)) state.naca = DEFAULTS.naca;

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
  sampleEvery = Math.max(1, Math.round(wt.stepsPerTStar / 40));
  tChange = 0;
  markers = [];
  analysis = null;
  diverged = false;
  lastCoeffs = { cx: NaN, cz: NaN, cm: NaN };
  view.attach(wt);
  hideToast();
  updateInfo();
}

function markChange() {
  tChange = wt.lbm.time;
  markers.push(wt.tStar);
  if (markers.length > 40) markers.shift();
}

function applyBody() {
  wt.setBody(state.shape, bodyParams());
  markChange();
  updateInfo();
}

function sample() {
  const c = wt.coefficients();
  lastCoeffs = c;
  hist.push(wt.lbm.time, c.cx, c.cz, c.cm);
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
    'Le calcul a divergé : Reynolds trop élevé pour cette résolution. Activez le modèle de sous-maille, réduisez Re ou affinez la grille, puis relancez.',
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
    font: v('--font-num') || 'monospace',
  };
  cpuView.setTheme(themeName(), themeColors);
  if (gpuView) gpuView.setTheme(themeName(), themeColors);
  updateColorbar();
  drawCharts();
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
  $('numinfo').textContent =
    `${wt.engine === 'gpu' ? 'GPU' : 'CPU'} · collision ${state.collision === 'rr' ? 'régularisée' : 'BGK'} · ` +
    `Grille ${i.nx}×${i.ny} · L = ${i.L} mailles · τ = ${i.tau.toFixed(4)} · Ma = ${i.mach.toFixed(2)} · ` +
    `confinement ${(i.blockage * 100).toFixed(1)} % · sous-maille ${i.les ? 'active' : 'inactive'} · ` +
    `Re fiable ≲ ${i.reLimit}`;
  const warn = $('warn');
  const msgs = [];
  if (state.re > i.reLimit) {
    const next =
      state.res === 'rapide'
        ? 'passez en résolution Standard ou Fine'
        : state.res === 'standard'
          ? 'passez en résolution Fine'
          : state.res === 'fine'
            ? 'passez en résolution Très fine (sur ordinateur) ou réduisez Re'
            : 'réduisez Re';
    msgs.push(
      `Au-delà de Re ≈ ${i.reLimit} avec ${i.L} mailles par longueur, la couche limite est trop peu résolue : ` +
        `la traînée est surestimée (de +20 % à +100 %) et les résultats sont qualitatifs ; ${next}.`,
    );
  }
  const s = shapeDef();
  const a = Math.abs(alphaOf());
  if (state.re >= 10000 && !['cylinder', 'square', 'triangle'].includes(s.id)) {
    msgs.push(
      'Vers Re ≈ 10⁴ apparaissent sur un profil réel des bulbes de décollement laminaire et la transition, ' +
        'phénomènes tridimensionnels que ce calcul 2D ne représente pas : résultats indicatifs, surtout en incidence.',
    );
  }
  const bluff = ['cylinder', 'square', 'triangle'].includes(s.id) || (['plate', 'ellipse'].includes(s.id) && a >= 20);
  if (bluff && state.re > RE_WAKE_3D) {
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
  if (a.periodic) return ['periodic', 'Périodique'];
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
  $('hud-time').textContent = `t* ${wt.tStar.toFixed(1)}`;

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
function perfRate() {
  const now = performance.now();
  if (now - perfT0 > 1000) {
    const sps = (perfSteps * 1000) / (now - perfT0);
    const mlups = (sps * wt.nx * wt.ny) / 1e6;
    const tps = (sps * wt.uIn) / wt.L;
    perfText = `${Math.round(sps)} pas/s · ${mlups.toFixed(1)} MLUPS · ${tps.toFixed(2)} t*/s`;
    perfSteps = 0;
    perfT0 = now;
  }
  return perfText;
}

function drawCharts() {
  if (!wt || !hist) return;
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
  drawPolar($('polar'), pts, state.polarMode, themeColors, { lift: themeColors.lift, drag: themeColors.accent }, v === 'alpha' ? '°' : '');
  $('polar-hint').textContent = pts.length
    ? `${pts.length} point${pts.length > 1 ? 's' : ''} pour ${wt.body.label} à Re ${state.re} (résolution ${RESOLUTIONS[state.res].label.toLowerCase()}).`
    : `La polaire trace les points enregistrés pour la configuration courante (forme, Re, résolution).`;
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
      String(p.re),
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
  return [state.shape, extra, spin, state.re, RESOLUTIONS[state.res].L, state.u, state.wall].join('|');
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
    source,
  };
  points.push(p);
  if (points.length > 500) points.shift();
  storageSet(POINTS_KEY, points);
  renderTable();
  drawCharts();
}

function copyTable() {
  const head = ['maquette', 'alpha_deg', 'rotation_OmegaR_U', 'Re', 'L_mailles', 'Cx', 'Cz', 'Cm', 'St', 'Cz_Cx'];
  const rows = points.map((p) =>
    [
      p.label,
      p.alpha === null ? '' : fmtCsv(p.alpha, 1),
      p.spin === null ? '' : fmtCsv(p.spin, 2),
      p.re,
      p.L,
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
  if (sweep.phase === 'settle') {
    if (el >= SWEEP_SETTLE_T) {
      sweep.phase = 'measure';
      sweep.t0 = wt.lbm.time;
    }
  } else {
    const a = analyze(hist.since(sweep.t0), wt.uIn, wt.L);
    const enough = (el >= SWEEP_MEASURE_T && (!a.periodic || a.periods >= 3)) || el >= 2.5 * SWEEP_MEASURE_T;
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

function frame() {
  if (pendingBody && !diverged) {
    pendingBody = false;
    applyBody();
  }
  if (running && !diverged) {
    const t0 = performance.now();
    let n = 0;
    const lbm = wt.lbm;
    if (wt.engine === 'gpu') {
      for (; n < gpuSteps; n++) {
        wt.step();
        if (lbm.time % sampleEvery === 0) lbm.sampleForces();
      }
      // Relecture des efforts : synchronise avec le GPU une fois par image.
      for (const smp of lbm.readSamples()) {
        const c = wt.coefficientsFrom(smp.fx, smp.fy, smp.mz);
        lastCoeffs = c;
        hist.push(smp.time, c.cx, c.cz, c.cm);
      }
      const dt = Math.max(1, performance.now() - t0);
      const k = Math.min(1.25, Math.max(0.6, GPU_FRAME_MS / dt));
      gpuSteps = Math.max(1, Math.min(4000, Math.round(gpuSteps * k + (k > 1 ? 1 : 0))));
    } else {
      do {
        wt.step();
        n++;
        if (lbm.time % sampleEvery === 0) sample();
      } while (performance.now() - t0 < FRAME_BUDGET_MS);
    }
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
    $('perf').textContent = perfRate() || '—';
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
  $('re-out').textContent = String(state.re);
  re.addEventListener('input', () => {
    state.re = RE_VALUES[Number(re.value)];
    $('re-out').textContent = String(state.re);
    wt.setReynolds(state.re);
    markChange();
    updateInfo();
    drawCharts();
    savePrefs();
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

  if (window.matchMedia) {
    const mq = window.matchMedia('(prefers-color-scheme: dark)');
    mq.addEventListener?.('change', applyTheme);
  }
  new MutationObserver(applyTheme).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
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
  buildSim();
  syncShapeControls();
  applyTheme();
  renderTable();
  syncPlay();
  requestAnimationFrame(frame);

  if ('serviceWorker' in navigator && !window.__SOUFFLERIE_INLINE__) {
    const local = location.hostname === 'localhost' || location.hostname === '127.0.0.1';
    if (location.protocol === 'https:' || local) navigator.serviceWorker.register('./sw.js').catch(() => {});
  }
}

start();
