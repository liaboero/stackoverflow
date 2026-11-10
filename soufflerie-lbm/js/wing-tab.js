/**
 * Onglet « Aile finie » : ligne portante non linéaire alimentée par les
 * polaires 2D de la méthode des panneaux (js/vlm.js, js/panel.js).
 */
import { buildWing, solveWing, polarAt } from './vlm.js';
import { panelSolve, sectionCoefficients } from './panel.js';
import { drawWall } from './chart.js';

const W_RE = [5e4, 1e5, 2e5, 5e5, 1e6, 2e6, 5e6, 1e7];
const KEY = 'soufflerie.wing.v1';
const DEF = {
  alpha: 4,
  naca: '2412',
  ar: 8,
  taper: 0.5,
  sweep: 0,
  dihedral: 0,
  twist: 0,
  re: 1e6,
  transition: 'free',
  chord: 0.2,
  fluid: 'air',
  polarMode: 'var',
};
const A_LO = -14;
const A_HI = 24;
const POLAR_LO = -6;
const POLAR_HI = 20;
const N_HALF = 40;
const $ = (id) => document.getElementById(id);
const DEG = Math.PI / 180;

/** Rampe de couleur Cl/Cl max : froid → chaud (rouge = décrochage local). */
const RAMP = [
  [0, [64, 112, 170]],
  [0.6, [96, 170, 160]],
  [0.85, [226, 190, 76]],
  [1, [206, 66, 46]],
];
function rampColor(t, alpha = 1) {
  const x = Math.min(1, Math.max(0, t));
  let k = 0;
  while (k < RAMP.length - 2 && x > RAMP[k + 1][0]) k++;
  const [x0, c0] = RAMP[k];
  const [x1, c1] = RAMP[k + 1];
  const f = (x - x0) / (x1 - x0 || 1);
  const c = c0.map((v, i) => Math.round(v + f * (c1[i] - v)));
  return `rgba(${c[0]}, ${c[1]}, ${c[2]}, ${alpha})`;
}

function validNaca(s) {
  if (!/^\d{4}$/.test(s)) return false;
  const t = Number(s.slice(2));
  return t >= 4 && t <= 30 && (s[0] === '0' || Number(s[1]) > 0);
}

export class WingTab {
  /**
   * @param {object} ctx { segmented, theme: () => couleurs, fmt, fmtRe, sig3, fluids, storageGet, storageSet }
   */
  constructor(ctx) {
    this.ctx = ctx;
    const saved = ctx.storageGet(KEY) || {};
    this.st = { ...DEF, ...saved };
    const st = this.st;
    if (!W_RE.includes(st.re)) st.re = DEF.re;
    if (!validNaca(st.naca)) st.naca = DEF.naca;
    if (!ctx.fluids[st.fluid]) st.fluid = DEF.fluid;
    this.tables = new Map();
    this.pending = null;
    this.result = null;
    this.polar = null;
    this.active = false;
    this.cam = { az: 32, el: 34 };
    this.gen = 0;
  }

  save() {
    this.ctx.storageSet(KEY, this.st);
  }

  // ------------------------------------------------------------ interface

  init() {
    const st = this.st;
    const fmt = this.ctx.fmt;
    const slider = (id, key, out, onChange = () => this.update()) => {
      const el = $(id);
      el.value = st[key];
      $(`${id}-out`).textContent = out(st[key]);
      el.addEventListener('input', () => {
        st[key] = Number(el.value);
        $(`${id}-out`).textContent = out(st[key]);
        this.save();
        onChange();
      });
    };
    slider('w-alpha', 'alpha', (v) => `${fmt(v, 1)}°`, () => this.update(false));
    slider('w-ar', 'ar', (v) => fmt(v, v % 1 ? 1 : 0));
    slider('w-taper', 'taper', (v) => fmt(v, 2));
    slider('w-sweep', 'sweep', (v) => `${fmt(v, 0)}°`);
    slider('w-dihedral', 'dihedral', (v) => `${fmt(v, 1)}°`);
    slider('w-twist', 'twist', (v) => `${fmt(v, 1)}°`);
    const re = $('w-re');
    re.max = W_RE.length - 1;
    re.value = W_RE.indexOf(st.re);
    $('w-re-out').textContent = this.ctx.fmtRe(st.re);
    re.addEventListener('input', () => {
      st.re = W_RE[Number(re.value)];
      $('w-re-out').textContent = this.ctx.fmtRe(st.re);
      this.save();
      this.update();
    });
    const naca = $('w-naca');
    naca.value = st.naca;
    const setNaca = (v) => {
      if (!validNaca(v)) {
        naca.setAttribute('aria-invalid', 'true');
        return;
      }
      naca.removeAttribute('aria-invalid');
      st.naca = v;
      naca.value = v;
      this.save();
      this.update();
    };
    naca.addEventListener('change', () => setNaca(naca.value.trim()));
    naca.addEventListener('input', () => {
      if (naca.value.length === 4) setNaca(naca.value.trim());
    });
    for (const b of $('w-naca-presets').querySelectorAll('[data-naca]')) b.addEventListener('click', () => setNaca(b.dataset.naca));
    this.ctx.segmented(
      $('w-trans-seg'),
      [
        { value: 'free', label: 'Libre (eᴺ)' },
        { value: 'turb', label: 'Bord d’attaque' },
      ],
      st.transition,
      (v) => {
        st.transition = v;
        this.save();
        this.update();
      },
    );
    this.ctx.segmented(
      $('w-fluid-seg'),
      Object.entries(this.ctx.fluids).map(([k, f]) => ({ value: k, label: f.label })),
      st.fluid,
      (v) => {
        st.fluid = v;
        this.save();
        this.showPhysical();
      },
    );
    const chord = $('w-chord');
    chord.value = st.chord;
    chord.addEventListener('input', () => {
      const v = Number(chord.value);
      if (v > 0 && v <= 20) {
        st.chord = v;
        chord.removeAttribute('aria-invalid');
        this.save();
        this.showPhysical();
      } else {
        chord.setAttribute('aria-invalid', 'true');
      }
    });
    this.ctx.segmented(
      $('w-polar-seg'),
      [
        { value: 'var', label: 'CL (α)' },
        { value: 'eiffel', label: 'CL (CD)' },
      ],
      st.polarMode,
      (v) => {
        st.polarMode = v;
        this.save();
        this.drawPolar();
      },
    );
    this.bindOrbit($('wing-view'));
    const ro = new ResizeObserver(() => this.draw());
    ['wing-view', 'w-span', 'w-polar'].forEach((id) => ro.observe($(id)));
    $('wing-grad').style.background = `linear-gradient(90deg, ${RAMP.map(([x, c]) => `rgb(${c.join(',')}) ${x * 100}%`).join(', ')})`;
  }

  bindOrbit(canvas) {
    let drag = null;
    canvas.addEventListener('pointerdown', (e) => {
      drag = { x: e.clientX, y: e.clientY, az: this.cam.az, el: this.cam.el };
      canvas.setPointerCapture(e.pointerId);
    });
    canvas.addEventListener('pointermove', (e) => {
      if (!drag) return;
      this.cam.az = drag.az - (e.clientX - drag.x) * 0.5;
      this.cam.el = Math.min(85, Math.max(-30, drag.el + (e.clientY - drag.y) * 0.4));
      this.drawView();
    });
    const end = () => (drag = null);
    canvas.addEventListener('pointerup', end);
    canvas.addEventListener('pointercancel', end);
  }

  activate() {
    this.active = true;
    this.update();
  }

  deactivate() {
    this.active = false;
  }

  // ------------------------------------------------------------ calcul

  wing() {
    const st = this.st;
    return buildWing({ ar: st.ar, taper: st.taper, sweep: st.sweep, dihedral: st.dihedral, twist: st.twist, n: N_HALF });
  }

  tableKey(re) {
    return `${this.st.naca}|${Math.round(re)}|${this.st.transition}`;
  }

  /** Reynolds des polaires nécessaires : emplanture et bout d'aile. */
  tableRes(wing) {
    const re = this.st.re;
    const r = [re * wing.cr];
    if (wing.ct < wing.cr * 0.98) r.push(re * wing.ct);
    return r;
  }

  /**
   * Relance le calcul. Les polaires de section manquantes sont calculées par
   * petits morceaux (≈ 20 ms) pour ne pas figer l'interface.
   * @param {boolean} [geometryChanged=true] si faux, seule l'incidence a changé
   */
  update(geometryChanged = true) {
    if (!this.active) return;
    const gen = ++this.gen;
    const wing = this.wing();
    this.wingGeo = wing;
    if (geometryChanged) this.polar = null;
    const missing = this.tableRes(wing).filter((re) => !this.tables.has(this.tableKey(re)));
    if (missing.length) {
      $('wing-status').textContent = 'Calcul des polaires de section…';
      this.computeTables(missing, gen, () => this.update(geometryChanged));
      this.draw();
      return;
    }
    this.section = this.makeSection(wing);
    const g0 = this.polar ? this.nearestGamma(this.st.alpha) : null;
    this.result = solveWing(wing, this.st.alpha, this.section, { gamma0: g0 });
    $('wing-status').textContent = this.result.converged ? 'Surface portante · Weissinger' : 'Non convergé (après décrochage)';
    this.showResult();
    this.draw();
    if (!this.polar) this.computePolar(wing, gen);
  }

  computeTables(list, gen, done) {
    const jobs = list.map((re) => ({ re, a: A_LO, rows: [] }));
    const xTrip = this.st.transition === 'turb' ? 0 : undefined;
    const naca = this.st.naca;
    const tick = () => {
      if (gen !== this.gen) return;
      const t0 = performance.now();
      while (jobs.length && performance.now() - t0 < 20) {
        const j = jobs[0];
        let row = null;
        try {
          const r = panelSolve('naca', { params: { naca }, alphaDeg: j.a, re: j.re, xTrip });
          if (r && Number.isFinite(r.cl) && Number.isFinite(r.cd)) row = sectionCoefficients(r);
        } catch {
          row = null;
        }
        j.rows.push(row);
        j.a++;
        if (j.a > A_HI) {
          this.tables.set(this.tableKey(j.re), finishTable(j.rows));
          jobs.shift();
        }
      }
      if (jobs.length) setTimeout(tick, 0);
      else done();
    };
    setTimeout(tick, 0);
  }

  /** Polaire de chaque bande : interpolation en log Re entre emplanture et bout. */
  makeSection(wing) {
    const res = this.tableRes(wing);
    const tR = this.tables.get(this.tableKey(res[0]));
    const tT = res.length > 1 ? this.tables.get(this.tableKey(res[1])) : tR;
    const lr = Math.log(res[0]);
    const lt = Math.log(res[res.length - 1]);
    const wgt = wing.strips.map((s) => {
      if (tR === tT) return 1;
      const li = Math.log(this.st.re * s.c);
      return Math.min(1, Math.max(0, (li - lt) / (lr - lt)));
    });
    const clmax = wing.strips.map((_, i) => wgt[i] * tR.clmax + (1 - wgt[i]) * tT.clmax);
    const fn = (i, a) => {
      const p = polarAt(tR, a);
      if (wgt[i] === 1) return p;
      const q = polarAt(tT, a);
      const w = wgt[i];
      return { cl: w * p.cl + (1 - w) * q.cl, cd: w * p.cd + (1 - w) * q.cd, cm: w * p.cm + (1 - w) * q.cm };
    };
    fn.clmax = clmax;
    fn.root = tR;
    return fn;
  }

  nearestGamma(a) {
    let best = null;
    for (const p of this.polar.pts) if (!best || Math.abs(p.alpha - a) < Math.abs(best.alpha - a)) best = p;
    return best ? best.gamma : null;
  }

  /** Polaire de l'aile, par continuation en incidence (un point par tranche). */
  computePolar(wing, gen) {
    const pts = [];
    let a = POLAR_LO;
    let g = null;
    const section = this.section;
    const tick = () => {
      if (gen !== this.gen) return;
      const t0 = performance.now();
      while (a <= POLAR_HI && performance.now() - t0 < 25) {
        const r = solveWing(wing, a, section, { gamma0: g });
        g = r.gamma;
        pts.push({ alpha: a, CL: r.CL, CD: r.CD, gamma: r.gamma, converged: r.converged });
        a++;
      }
      if (a <= POLAR_HI) setTimeout(tick, 0);
      else {
        this.polar = { pts };
        this.drawPolar();
      }
      this.drawPolar(pts);
    };
    setTimeout(tick, 0);
  }

  // ------------------------------------------------------------ affichage

  showResult() {
    const r = this.result;
    const { fmt } = this.ctx;
    $('w-cl').textContent = fmt(r.CL);
    $('w-cl-sub').textContent = `Cm ${fmt(r.Cm)}`;
    $('w-cd').textContent = fmt(r.CD, 4);
    $('w-cd-sub').textContent = `CDi ${fmt(r.CDi, 4)}`;
    $('w-ld').textContent = Math.abs(r.CD) > 1e-6 ? fmt(r.CL / r.CD, 1) : '—';
    $('w-ld-sub').textContent = `CDp ${fmt(r.CDv, 4)}`;
    $('w-e').textContent = Number.isFinite(r.e) && Math.abs(r.CL) > 0.05 ? fmt(r.e, 3) : '—';
    $('w-e-sub').textContent = `A = ${fmt(this.st.ar, this.st.ar % 1 ? 1 : 0)}`;
    this.showPhysical();
    this.showWarnings();
  }

  showPhysical() {
    const st = this.st;
    const f = this.ctx.fluids[st.fluid];
    const V = (st.re * f.nu) / st.chord;
    const b = st.ar * st.chord;
    const S = b * st.chord;
    const q = 0.5 * f.rho * V * V;
    const { sig3 } = this.ctx;
    const parts = [`envergure ${sig3(b)} m`, `S = ${sig3(S)} m²`, `V = ${sig3(V)} m/s`];
    if (this.result) parts.push(`portance ${sig3(q * S * this.result.CL)} N`, `traînée ${sig3(q * S * this.result.CD)} N`);
    if (st.fluid === 'air' && V / 340 > 0.3) parts.push(`Mach ${sig3(V / 340)} : compressibilité non prise en compte`);
    $('w-phys').textContent = parts.join(' · ');
  }

  showWarnings() {
    const st = this.st;
    const r = this.result;
    const msgs = [];
    if (r && !r.converged) msgs.push('Solution non convergée : après le décrochage, la ligne portante admet plusieurs solutions ; résultat indicatif.');
    if (r && this.section) {
      let worst = 0;
      let eta = 0;
      r.strips.forEach((s, i) => {
        const k = s.cl / this.section.clmax[i];
        if (k > worst) {
          worst = k;
          eta = Math.abs(s.eta);
        }
      });
      if (worst > 0.97) msgs.push(`Décrochage local amorcé vers 2y/b ≈ ${this.ctx.fmt(eta, 2)} (Cl local au maximum de la polaire du profil).`);
    }
    if (Math.abs(st.sweep) > 30) msgs.push('Flèche > 30° : la ligne portante représente mal le centre d’une aile en flèche ; ordre de grandeur seulement.');
    if (st.ar < 4) msgs.push('Allongement < 4 : la ligne portante surestime la pente de portance (une méthode de surface portante serait plus juste).');
    if (st.re < 2e5) msgs.push('Re < 2·10⁵ : bulbes de décollement laminaire non modélisés dans les polaires de section.');
    const w = $('w-warn');
    w.hidden = !msgs.length;
    w.textContent = msgs.join(' ');
  }

  draw() {
    if (!this.active) return;
    this.drawView();
    this.drawSpan();
    this.drawPolar(this.polar?.pts);
  }

  drawSpan() {
    const r = this.result;
    const th = this.ctx.theme();
    if (!r || !this.section) {
      drawWall($('w-span'), { series: [], xLo: -1, xHi: 1, empty: 'Calcul en cours…' }, th);
      return;
    }
    const ell = [];
    for (let k = 0; k <= 80; k++) {
      const e = -1 + k / 40;
      ell.push([e, ((4 * r.CL) / Math.PI) * Math.sqrt(Math.max(0, 1 - e * e))]);
    }
    const S = r.strips;
    drawWall(
      $('w-span'),
      {
        series: [
          { pts: ell, color: th.muted, style: 'dash' },
          { pts: S.map((s, i) => [s.eta, this.section.clmax[i]]), color: th.accent, style: 'dash' },
          { pts: S.map((s) => [s.eta, s.cl]), color: th.accent, style: 'line' },
          { pts: S.map((s) => [s.eta, s.ccl]), color: th.lift, style: 'line' },
        ],
        xLo: -1,
        xHi: 1,
        yLabel: 'abscisse 2y/b',
      },
      th,
    );
    $('w-span-hint').textContent =
      `Charge c·Cl/c̄ (proportionnelle à la circulation) comparée à la répartition elliptique de même CL, qui minimise la traînée induite. ` +
      `Le décrochage commence là où le Cl local atteint le Cl max du profil.`;
  }

  drawPolar(pts) {
    if (!this.active) return;
    const th = this.ctx.theme();
    const eiffel = this.st.polarMode === 'eiffel';
    const series = [];
    const marks = [];
    const sec = this.section?.root;
    if (sec) {
      const sp = [];
      for (let k = 0; k < sec.al.length; k++) {
        const a = sec.al[k];
        if (a < POLAR_LO || a > POLAR_HI) continue;
        sp.push(eiffel ? [sec.cd[k], sec.cl[k]] : [a, sec.cl[k]]);
      }
      series.push({ pts: sp, color: th.muted, style: 'dash' });
    }
    if (pts && pts.length) {
      series.push({ pts: pts.map((p) => (eiffel ? [p.CD, p.CL] : [p.alpha, p.CL])), color: th.lift, style: 'line' });
      if (!eiffel) series.push({ pts: pts.map((p) => [p.alpha, p.CD]), color: th.accent, style: 'line' });
    }
    if (this.result) {
      if (eiffel) series.push({ pts: [[this.result.CD, this.result.CL]], color: th.ink, style: 'dots' });
      else marks.push({ x: this.st.alpha, color: th.ink });
    }
    let xLo = POLAR_LO;
    let xHi = POLAR_HI;
    if (eiffel) {
      xLo = 0;
      xHi = 0.01;
      for (const s of series) for (const [x] of s.pts) if (Number.isFinite(x) && x < 0.3) xHi = Math.max(xHi, x);
      xHi *= 1.1;
    }
    drawWall($('w-polar'), { series, xLo, xHi, xUnit: eiffel ? '' : '°', marks, empty: 'Calcul en cours…' }, th);
    const best = pts && pts.length ? pts.reduce((m, p) => (p.CD > 0 && p.CL / p.CD > m.ld ? { ld: p.CL / p.CD, a: p.alpha, cl: p.CL } : m), { ld: -Infinity }) : null;
    const clmax = pts && pts.length ? pts.reduce((m, p) => (p.CL > m.CL ? p : m), pts[0]) : null;
    const { fmt } = this.ctx;
    $('w-polar-hint').textContent =
      (eiffel ? 'Polaire d’Eiffel CL(CD) de l’aile ; tirets : profil 2D au Reynolds d’emplanture. ' : 'CL et CD de l’aile en fonction de α ; tirets : Cl du profil 2D. ') +
      (best && Number.isFinite(best.ld) ? `Finesse max ≈ ${fmt(best.ld, 1)} à α = ${best.a}° (CL = ${fmt(best.cl, 2)}). ` : '') +
      (clmax && this.polar ? `CL max ≈ ${fmt(clmax.CL, 2)} à α = ${clmax.alpha}° (Cl max du profil : ${fmt(sec?.clmax ?? NaN, 2)}).` : '');
  }

  /** Vue 3D : surface colorée par Cl/Cl max, charge, nappe tourbillonnaire. */
  drawView() {
    const canvas = $('wing-view');
    const dpr = Math.min(window.devicePixelRatio || 1, 2.5);
    const W = Math.max(1, Math.round(canvas.clientWidth * dpr));
    const H = Math.max(1, Math.round(canvas.clientHeight * dpr));
    if (canvas.width !== W || canvas.height !== H) {
      canvas.width = W;
      canvas.height = H;
    }
    const g = canvas.getContext('2d');
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    const w = W / dpr;
    const h = H / dpr;
    g.clearRect(0, 0, w, h);
    const th = this.ctx.theme();
    const wing = this.wingGeo || this.wing();
    const r = this.result;
    const a = (this.st.alpha || 0) * DEG;
    const az = this.cam.az * DEG;
    const el = this.cam.el * DEG;
    // Incidence : rotation cabreuse autour de y (bord de fuite vers le bas).
    const pitch = ([x, y, z]) => [x * Math.cos(a) + z * Math.sin(a), y, -x * Math.sin(a) + z * Math.cos(a)];
    // Projection d'un point du repère « soufflerie » (écoulement amont selon x).
    const viewW = ([x, y, z]) => {
      const x1 = x * Math.cos(az) - y * Math.sin(az);
      const y1 = x * Math.sin(az) + y * Math.cos(az);
      return { u: y1, v: z * Math.cos(el) + x1 * Math.sin(el), depth: x1 * Math.cos(el) - z * Math.sin(el) };
    };
    const half = wing.b / 2;
    const span = wing.b;
    // Cadrage automatique sur l'aile, sa charge et le début du sillage.
    let u0 = Infinity;
    let u1 = -Infinity;
    let v0 = Infinity;
    let v1 = -Infinity;
    const fit = (p) => {
      const q = viewW(p);
      u0 = Math.min(u0, q.u);
      u1 = Math.max(u1, q.u);
      v0 = Math.min(v0, q.v);
      v1 = Math.max(v1, q.v);
    };
    for (const y of [-half, 0, half]) {
      const c4 = wing.c4At(y);
      const c = wing.chordAt(y);
      fit(pitch([c4[0] - 0.25 * c, y, c4[2]]));
      fit(pitch([c4[0] + 0.75 * c, y, c4[2]]));
      fit(pitch([c4[0] + 0.75 * c, y, c4[2] + 0.16 * span]));
      fit([c4[0] + 0.5 * span, y, c4[2]]);
    }
    fit([-0.45 * span, -half * 0.9, 0]);
    const sc = Math.min((0.9 * w) / (u1 - u0 || 1), (0.78 * h) / (v1 - v0 || 1));
    const cx = w / 2 - ((u0 + u1) / 2) * sc;
    const cy = h * 0.52 + ((v0 + v1) / 2) * sc;
    const PW = (p) => {
      const q = viewW(p);
      return [cx + q.u * sc, cy - q.v * sc, q.depth];
    };
    // Point du repère de l'aile (tangué de α).
    const P = (p) => PW(pitch(p));
    // Nœuds d'envergure : bord d'attaque et bord de fuite (vrillage compris).
    const ys = [];
    const n = N_HALF;
    for (let k = 0; k <= 2 * n; k++) ys.push(half * ((2 * k) / (2 * n) - 1));
    const edge = (y) => {
      const c = wing.chordAt(y);
      const p = wing.c4At(y);
      const t = (this.st.twist * DEG * Math.abs(y)) / half;
      const ux = Math.cos(t);
      const uz = -Math.sin(t);
      return {
        le: [p[0] - 0.25 * c * ux, p[1], p[2] - 0.25 * c * uz],
        te: [p[0] + 0.75 * c * ux, p[1], p[2] + 0.75 * c * uz],
        c4: p,
      };
    };
    const nodes = ys.map(edge);
    // Sillage : tourbillons libres, intensité ∝ saut de circulation.
    if (r) {
      const G = r.gamma;
      const jumps = [];
      for (let k = 0; k <= 2 * n; k++) jumps.push(Math.abs((k < 2 * n ? G[k] : 0) - (k > 0 ? G[k - 1] : 0)));
      const jmax = Math.max(...jumps, 1e-9);
      // Les tourbillons libres suivent l'écoulement amont (horizontal).
      for (let k = 0; k <= 2 * n; k++) {
        const s = jumps[k] / jmax;
        if (s < 0.04) continue;
        const te = pitch(nodes[k].te);
        const start = PW(te);
        const end = PW([te[0] + 0.5 * span, te[1], te[2]]);
        g.strokeStyle = th.lift;
        g.globalAlpha = 0.15 + 0.6 * s;
        g.lineWidth = 0.6 + 1.6 * s;
        g.beginPath();
        g.moveTo(start[0], start[1]);
        g.lineTo(end[0], end[1]);
        g.stroke();
      }
      g.globalAlpha = 1;
    }
    // Surface : quadrilatères triés du plus loin au plus proche.
    const quads = [];
    for (let i = 0; i < 2 * n; i++) {
      const A = nodes[i];
      const B = nodes[i + 1];
      const pts = [P(A.le), P(B.le), P(B.te), P(A.te)];
      const depth = (pts[0][2] + pts[1][2] + pts[2][2] + pts[3][2]) / 4;
      const ratio = r && this.section ? r.strips[i].cl / this.section.clmax[i] : 0.3;
      quads.push({ pts, depth, ratio });
    }
    quads.sort((p, q) => q.depth - p.depth);
    for (const qd of quads) {
      g.fillStyle = rampColor(qd.ratio, 0.92);
      g.strokeStyle = th.bodyEdge || th.ink;
      g.lineWidth = 0.4;
      g.beginPath();
      qd.pts.forEach(([x, y], k) => (k ? g.lineTo(x, y) : g.moveTo(x, y)));
      g.closePath();
      g.fill();
      g.stroke();
    }
    // Contour.
    g.strokeStyle = th.ink;
    g.lineWidth = 1.1;
    g.beginPath();
    nodes.forEach((nd, k) => {
      const [x, y] = P(nd.le);
      if (k) g.lineTo(x, y);
      else g.moveTo(x, y);
    });
    for (let k = nodes.length - 1; k >= 0; k--) {
      const [x, y] = P(nodes[k].te);
      g.lineTo(x, y);
    }
    g.closePath();
    g.stroke();
    // Charge : barres verticales (repère aile) proportionnelles à c·Cl.
    if (r) {
      const cmax = Math.max(...r.strips.map((s) => Math.abs(s.ccl)), 1e-6);
      const hScale = (0.15 * span) / cmax;
      g.strokeStyle = th.lift;
      g.globalAlpha = 0.45;
      g.lineWidth = 1;
      g.beginPath();
      r.strips.forEach((s, i) => {
        if (i % 2) return;
        const y = 0.5 * (ys[i] + ys[i + 1]);
        const p = wing.c4At(y);
        const top = [p[0], p[1], p[2] + s.ccl * hScale];
        const [x0, y0] = P(p);
        const [x1, y1] = P(top);
        g.moveTo(x0, y0);
        g.lineTo(x1, y1);
      });
      g.stroke();
      g.globalAlpha = 1;
      g.lineWidth = 1.8;
      g.beginPath();
      r.strips.forEach((s, i) => {
        const y = 0.5 * (ys[i] + ys[i + 1]);
        const p = wing.c4At(y);
        const [x1, y1] = P([p[0], p[1], p[2] + s.ccl * hScale]);
        if (i) g.lineTo(x1, y1);
        else g.moveTo(x1, y1);
      });
      g.stroke();
    }
    // Écoulement amont.
    const p0 = PW([-0.45 * span, -half * 0.9, 0]);
    const p1 = PW([-0.1 * span, -half * 0.9, 0]);
    g.strokeStyle = th.muted;
    g.fillStyle = th.muted;
    g.lineWidth = 1.5;
    g.beginPath();
    g.moveTo(p0[0], p0[1]);
    g.lineTo(p1[0], p1[1]);
    g.stroke();
    const ang = Math.atan2(p1[1] - p0[1], p1[0] - p0[0]);
    g.beginPath();
    g.moveTo(p1[0], p1[1]);
    g.lineTo(p1[0] - 8 * Math.cos(ang - 0.4), p1[1] - 8 * Math.sin(ang - 0.4));
    g.lineTo(p1[0] - 8 * Math.cos(ang + 0.4), p1[1] - 8 * Math.sin(ang + 0.4));
    g.closePath();
    g.fill();
    g.font = `11px ${th.font}`;
    g.fillText('V∞', p0[0] - 4, p0[1] - 8);
  }
}

/** Table complète : trous éventuels comblés par le voisin, Cl max relevé. */
function finishTable(rows) {
  const al = [];
  const cl = [];
  const cd = [];
  const cm = [];
  let last = null;
  rows.forEach((r, k) => {
    const row = r || last || { cl: 0, cd: 0.02, cm: 0 };
    last = row;
    al.push(A_LO + k);
    cl.push(row.cl);
    cd.push(row.cd);
    cm.push(row.cm);
  });
  return { lo: A_LO, hi: A_HI, al, cl, cd, cm, clmax: Math.max(...cl) };
}
