/**
 * Onglet « Soufflerie 3D » : pilotage du calcul D3Q19 WebGPU et de son rendu.
 */
import { Lbm3D, webgpuAdapter, bytesPerCell } from './lbm3d.js';
import { makeBody3D, BODIES3D } from './bodies3d.js';
import { Render3D, cssRgb } from './render3d.js';
import { ForceHistory, analyze } from './stats.js';
import { drawHistory } from './chart.js';

const KEY = 'soufflerie.3d.v1';
const RE3 = [10, 20, 50, 100, 200, 300, 500, 1000, 2000];
const RES3 = [
  { L: 8, label: 'Aperçu' },
  { L: 12, label: 'Rapide' },
  { L: 16, label: 'Standard' },
  { L: 24, label: 'Fine' },
  { L: 32, label: 'Très fine' },
];
const U = 0.1;
const RAMP_T = 4;
const SETTLE_T = 6;
const WINDOW_T = 20;
const DEF = { shape: 'sphere', re: 300, L: 16, alpha: 8, field: 'vort', plane: 'y', pos: null, iso: true, isoLevel: 0.4, fp16: true };
const $ = (id) => document.getElementById(id);

export class Lab3D {
  /** @param {object} ctx { segmented, theme, fmt, sig3, storageGet, storageSet } */
  constructor(ctx) {
    this.ctx = ctx;
    this.st = { ...DEF, ...(ctx.storageGet(KEY) || {}) };
    if (!BODIES3D[this.st.shape]) this.st.shape = DEF.shape;
    if (!RE3.includes(this.st.re)) this.st.re = DEF.re;
    if (!RES3.some((r) => r.L === this.st.L)) this.st.L = DEF.L;
    this.cam = { az: 62, el: 24, dist: 6 };
    this.active = false;
    this.running = true;
    this.ready = false;
    this.sim = null;
    this.inflight = false;
    this.steps = 4;
    this.perf = { steps: 0, t0: performance.now(), mlups: 0 };
    this.raf = 0;
  }

  save() {
    this.ctx.storageSet(KEY, this.st);
  }

  msg(text) {
    const m = $('s3-msg');
    m.hidden = !text;
    m.textContent = text || '';
  }

  async init() {
    const st = this.st;
    const fmt = this.ctx.fmt;
    const seg = this.ctx.segmented;
    seg(
      $('s3-shape-seg'),
      Object.entries(BODIES3D).map(([k, b]) => ({ value: k, label: b.label })),
      st.shape,
      (v) => {
        st.shape = v;
        st.pos = null;
        this.save();
        this.rebuild();
      },
    );
    const alpha = $('s3-alpha');
    alpha.value = st.alpha;
    $('s3-alpha-out').textContent = `${st.alpha}°`;
    alpha.addEventListener('input', () => {
      st.alpha = Number(alpha.value);
      $('s3-alpha-out').textContent = `${st.alpha}°`;
      this.save();
      clearTimeout(this.alphaTimer);
      this.alphaTimer = setTimeout(() => this.rebuild(), 300);
    });
    const re = $('s3-re');
    re.max = RE3.length - 1;
    re.value = RE3.indexOf(st.re);
    $('s3-re-out').textContent = String(st.re);
    re.addEventListener('input', () => {
      st.re = RE3[Number(re.value)];
      $('s3-re-out').textContent = String(st.re);
      this.save();
      if (this.sim) {
        this.sim.setViscosity((U * st.L) / st.re);
        this.markChange();
        this.showInfo();
      }
    });
    seg(
      $('s3-field-seg'),
      [
        { value: 'speed', label: 'Vitesse' },
        { value: 'vort', label: 'Vorticité' },
      ],
      st.field,
      (v) => {
        st.field = v;
        this.save();
        this.legend();
      },
    );
    seg(
      $('s3-plane-seg'),
      [
        { value: 'y', label: 'Longitudinal' },
        { value: 'z', label: 'Horizontal' },
        { value: 'x', label: 'Transversal' },
      ],
      st.plane,
      (v) => {
        st.plane = v;
        st.pos = null;
        this.save();
        this.syncPos();
        this.legend();
      },
    );
    $('s3-pos').addEventListener('input', () => {
      st.pos = Number($('s3-pos').value);
      this.save();
      this.syncPos();
    });
    $('s3-iso').checked = st.iso;
    $('s3-iso').addEventListener('change', () => {
      st.iso = $('s3-iso').checked;
      this.save();
    });
    const lvl = $('s3-iso-level');
    lvl.value = st.isoLevel;
    $('s3-iso-out').textContent = `Q = ${fmt(st.isoLevel, 2)} (U/L)²`;
    lvl.addEventListener('input', () => {
      st.isoLevel = Number(lvl.value);
      $('s3-iso-out').textContent = `Q = ${fmt(st.isoLevel, 2)} (U/L)²`;
      this.save();
    });
    $('s3-play').addEventListener('click', () => {
      this.running = !this.running;
      this.syncPlay();
    });
    $('s3-reset').addEventListener('click', () => this.rebuild());
    this.bindCamera($('view3d'));
  }

  syncPlay() {
    const b = $('s3-play');
    b.setAttribute('aria-label', this.running ? 'Mettre en pause' : 'Reprendre le calcul');
    b.querySelector('.ic-pause').toggleAttribute('hidden', !this.running);
    b.querySelector('.ic-play').toggleAttribute('hidden', this.running);
  }

  bindCamera(canvas) {
    const ptrs = new Map();
    let base = null;
    canvas.addEventListener('pointerdown', (e) => {
      ptrs.set(e.pointerId, { x: e.clientX, y: e.clientY });
      canvas.setPointerCapture(e.pointerId);
      base = { ...this.cam, pts: [...ptrs.values()].map((p) => ({ ...p })) };
    });
    canvas.addEventListener('pointermove', (e) => {
      if (!ptrs.has(e.pointerId) || !base) return;
      ptrs.set(e.pointerId, { x: e.clientX, y: e.clientY });
      const pts = [...ptrs.values()];
      if (pts.length === 1 && base.pts.length === 1) {
        this.cam.az = base.az - (pts[0].x - base.pts[0].x) * 0.4;
        this.cam.el = Math.min(85, Math.max(-85, base.el + (pts[0].y - base.pts[0].y) * 0.3));
      } else if (pts.length >= 2 && base.pts.length >= 2) {
        const d0 = Math.hypot(base.pts[0].x - base.pts[1].x, base.pts[0].y - base.pts[1].y);
        const d1 = Math.hypot(pts[0].x - pts[1].x, pts[0].y - pts[1].y);
        if (d0 > 10) this.cam.dist = Math.min(20, Math.max(1.5, (base.dist * d0) / d1));
      }
    });
    const end = (e) => {
      ptrs.delete(e.pointerId);
      base = { ...this.cam, pts: [...ptrs.values()].map((p) => ({ ...p })) };
    };
    canvas.addEventListener('pointerup', end);
    canvas.addEventListener('pointercancel', end);
    canvas.addEventListener(
      'wheel',
      (e) => {
        e.preventDefault();
        this.cam.dist = Math.min(20, Math.max(1.5, this.cam.dist * Math.exp(e.deltaY * 0.001)));
      },
      { passive: false },
    );
  }

  async activate() {
    this.active = true;
    if (!this.device && !this.starting) {
      this.starting = true;
      await this.start();
      this.starting = false;
    }
    if (this.device && !this.raf) this.raf = requestAnimationFrame((t) => this.frame(t));
  }

  deactivate() {
    this.active = false;
  }

  async start() {
    const adapter = await webgpuAdapter();
    if (!adapter) {
      this.msg(
        'WebGPU n’est pas disponible dans ce navigateur : la soufflerie 3D demande Chrome ou Edge récents (ordinateur, Android 12+), Safari 26 (iPhone, iPad, Mac) ou Firefox récent (Windows, Mac).',
      );
      $('s3-status').textContent = 'WebGPU indisponible';
      return;
    }
    this.limits = {
      buffer: Math.min(adapter.limits.maxBufferSize, adapter.limits.maxStorageBufferBindingSize),
      tex: adapter.limits.maxTextureDimension3D,
    };
    try {
      this.device = await adapter.requestDevice({
        requiredLimits: {
          maxStorageBufferBindingSize: adapter.limits.maxStorageBufferBindingSize,
          maxBufferSize: adapter.limits.maxBufferSize,
        },
      });
    } catch (e) {
      this.msg(`Impossible d’ouvrir le GPU : ${e.message}`);
      return;
    }
    this.device.addEventListener?.('uncapturederror', (e) => {
      console.error(`WebGPU : ${e.error.message}`);
      this.msg(`Erreur WebGPU : ${e.error.message.split('\n')[0]}`);
    });
    this.device.lost?.then((info) => {
      this.msg(`Le GPU a été réinitialisé (${info.message || info.reason}). Rechargez la page.`);
      this.device = null;
    });
    // Rendu hors écran : réservé aux tests automatiques (navigateur sans affichage).
    this.renderer = new Render3D(this.device, $('view3d'), { offscreen: !!this.st.offscreen });
    this.buildResSeg();
    this.rebuild();
  }

  /** Mailles et mémoire d'une résolution pour la maquette courante. */
  sizeOf(L) {
    const b = BODIES3D[this.st.shape];
    const N = Math.round(b.domain[0] * L) * Math.round(b.domain[1] * L) * Math.round(b.domain[2] * L);
    const bytes = N * (bytesPerCell(this.st.fp16) + 16);
    const perBuffer = N * (this.st.fp16 ? 20 : 40);
    const dims = b.domain.map((d) => Math.round(d * L));
    const ok = !this.limits || (perBuffer <= this.limits.buffer && Math.max(...dims) <= this.limits.tex);
    return { N, bytes, ok };
  }

  buildResSeg() {
    const opts = RES3.map((r) => ({ value: r.L, label: `${r.label} · ${r.L}` }));
    this.ctx.segmented($('s3-res-seg'), opts, this.st.L, (v) => {
      this.st.L = v;
      this.save();
      this.rebuild();
    });
    this.syncResSeg();
  }

  syncResSeg() {
    for (const r of RES3) {
      const btn = $('s3-res-seg').querySelector(`[data-value="${r.L}"]`);
      if (!btn) continue;
      const s = this.sizeOf(r.L);
      btn.disabled = !s.ok;
      btn.title = `${(s.N / 1e6).toFixed(2)} M mailles, ≈ ${Math.round(s.bytes / 1e6)} Mo${s.ok ? '' : ' : dépasse la mémoire GPU disponible'}`;
    }
    const s = this.sizeOf(this.st.L);
    $('s3-res-hint').textContent = `${(s.N / 1e6).toFixed(2)} million de mailles, ≈ ${Math.round(s.bytes / 1e6)} Mo de mémoire GPU. Chaque palier coûte ≈ 3 à 5 fois plus.`;
  }

  rebuild() {
    if (!this.device) return;
    const st = this.st;
    if (!this.sizeOf(st.L).ok) {
      const fit = [...RES3].reverse().find((r) => this.sizeOf(r.L).ok);
      st.L = fit ? fit.L : 8;
    }
    this.syncResSeg();
    this.sim?.dispose();
    this.sim = null;
    const L = st.L;
    const body = makeBody3D(st.shape, L, { alpha: st.alpha });
    $('s3-row-alpha').hidden = !BODIES3D[st.shape].alpha;
    $('s3-shape-hint').textContent = `Longueur de référence : ${body.ref} = ${L} mailles. Veine ${BODIES3D[st.shape].domain.join(' × ')} L (x, y, z), bords latéraux périodiques.`;
    try {
      this.sim = new Lbm3D(this.device, { nx: body.nx, ny: body.ny, nz: body.nz, nu: (U * L) / st.re, uIn: U, fp16: st.fp16 });
      this.sim.setGeometry(body.inside, body.box);
      this.sim.reset(0);
    } catch (e) {
      this.msg(`Calcul 3D impossible : ${e.message}`);
      return;
    }
    this.body = body;
    this.renderer.attach(this.sim, body);
    this.hist = new ForceHistory(8192);
    this.tChange = 0;
    this.steps = 4;
    this.inflight = false;
    this.msg('');
    this.syncPos();
    this.legend();
    this.showInfo();
    this.ready = true;
  }

  markChange() {
    this.tChange = this.sim.time;
  }

  uAt(t) {
    const ramp = (RAMP_T * this.st.L) / U;
    return t >= ramp ? U : 0.5 * U * (1 - Math.cos((Math.PI * t) / ramp));
  }

  /** Position du plan de coupe (mailles), par défaut au centre de la maquette ou dans le sillage. */
  planePos() {
    const st = this.st;
    const s = this.sim;
    const ax = { x: 0, y: 1, z: 2 }[st.plane];
    const n = [s.nx, s.ny, s.nz][ax];
    if (st.pos === null || st.pos === undefined) {
      const c = this.body.center[ax];
      return st.plane === 'x' ? Math.min(n - 2, c + 2.5 * this.body.L) : c;
    }
    return st.pos * (n - 1);
  }

  syncPos() {
    if (!this.sim) return;
    const ax = { x: 0, y: 1, z: 2 }[this.st.plane];
    const n = [this.sim.nx, this.sim.ny, this.sim.nz][ax];
    const p = this.planePos();
    $('s3-pos').value = p / (n - 1);
    const rel = (p - this.body.center[ax]) / this.body.L;
    $('s3-pos-out').textContent = `${rel >= 0 ? '+' : '−'}${this.ctx.fmt(Math.abs(rel), 1)} L`;
  }

  legend() {
    const st = this.st;
    const comp = { x: 'ωx', y: 'ωy', z: 'ωz' }[st.plane];
    $('s3-legend').textContent = st.field === 'speed' ? 'Coupe : |u|/U de 0 à 1,6' : `Coupe : ${comp}·L/U de −4 à +4`;
  }

  showInfo() {
    const s = this.sim;
    if (!s) return;
    $('s3-info').textContent =
      `Grille ${s.nx}×${s.ny}×${s.nz} (${(s.N / 1e6).toFixed(2)} M mailles) · ${s.wallCount} mailles pariétales · ` +
      `τ = ${s.tau0.toFixed(4)} · sous-maille ${s.les ? 'active' : 'inactive'} · ${s.fp16 ? 'demi-précision' : 'simple précision'}` +
      (this.perf.mlups ? ` · ${this.ctx.sig3(this.perf.mlups)} MLUPS` : '');
  }

  async frame(ts) {
    this.raf = 0;
    if (!this.active || !this.device) return;
    this.raf = requestAnimationFrame((t) => this.frame(t));
    if (!this.sim || this.inflight) return;
    const sim = this.sim;
    let n = 0;
    if (this.running) {
      n = this.steps;
      sim.step(n, (t) => this.uAt(t));
    }
    const th = this.ctx.theme();
    this.renderer.render({
      az: this.cam.az,
      el: this.cam.el,
      dist: this.cam.dist,
      plane: { axis: this.st.plane, pos: this.planePos() },
      field: this.st.field,
      iso: { on: this.st.iso, level: this.st.isoLevel, step: 0.8 },
      theme: { bg: cssRgb(th.stage, [0.9, 0.92, 0.93]), body: cssRgb(th.body, [0.2, 0.25, 0.3]), line: cssRgb(th.muted) },
      U,
    });
    this.inflight = true;
    const t0 = performance.now();
    await this.device.queue.onSubmittedWorkDone();
    this.inflight = false;
    const dt = performance.now() - t0;
    // Pas par image : ≈ 25 ms de calcul GPU, sans laisser l'affichage décrocher.
    if (this.running) {
      if (dt < 20) this.steps = Math.min(400, Math.ceil(this.steps * 1.25));
      else if (dt > 40) this.steps = Math.max(1, Math.floor(this.steps * 0.75));
      this.perf.steps += n;
    }
    const now = performance.now();
    if (now - this.perf.t0 > 1000) {
      this.perf.mlups = (this.perf.steps * sim.N) / ((now - this.perf.t0) * 1e3);
      this.perf.steps = 0;
      this.perf.t0 = now;
      this.showInfo();
    }
    if (!this.lastRead || now - this.lastRead > 250) {
      this.lastRead = now;
      const f = await sim.readForce();
      if (f && sim === this.sim) this.addForce(f);
    }
  }

  addForce(f) {
    const q = 0.5 * U * U * this.body.area;
    this.hist.push(f.time, f.fx / q, f.fz / q, f.fy / q);
    this.updateReadouts();
  }

  updateReadouts() {
    const sim = this.sim;
    const spt = this.st.L / U;
    const tStar = sim.time / spt;
    const { fmt, sig3 } = this.ctx;
    $('s3-time').textContent = `t* ${tStar.toFixed(1)}`;
    const settle = this.tChange === 0 ? RAMP_T + SETTLE_T : SETTLE_T;
    const start = Math.max(this.tChange + settle * spt, sim.time - WINDOW_T * spt);
    const a = analyze(this.hist.since(start), U, this.st.L);
    const ready = sim.time - start >= 3 * spt && a.n >= 8 && Number.isFinite(a.cx);
    const status = $('s3-status');
    let label = 'Transitoire';
    let state = 'transient';
    if (!this.running) [state, label] = ['paused', 'En pause'];
    else if (sim.time < (RAMP_T * this.st.L) / U) label = 'Mise en vitesse';
    else if (ready) [state, label] = a.periodic && a.czAmp > 0.01 ? ['periodic', 'Instationnaire'] : ['steady', 'Établi'];
    status.dataset.state = state;
    status.textContent = label;
    const last = this.hist.since(sim.time - spt);
    const lastV = (k) => (last[k].length ? last[k][last[k].length - 1] : NaN);
    $('s3-cx').textContent = fmt(ready ? a.cx : lastV('cx'));
    $('s3-cz').textContent = fmt(ready ? a.cz : lastV('cz'));
    $('s3-cy').textContent = fmt(ready ? a.cm : lastV('cm'));
    $('s3-cx-sub').textContent = ready ? (a.periodic && a.cxAmp > 0.001 ? `± ${fmt(a.cxAmp)}` : 'moyenne') : 'instantané';
    $('s3-cz-sub').textContent = ready ? (a.periodic && a.czAmp > 0.001 ? `± ${fmt(a.czAmp)}` : 'moyenne') : 'instantané';
    $('s3-cy-sub').textContent = ready ? 'moyenne' : 'instantané';
    $('s3-st').textContent = ready && a.periodic && Number.isFinite(a.st) ? fmt(a.st, 3) : '—';
    $('s3-st-sub').textContent = ready && a.periodic ? `${a.periods} périodes` : '';
    $('s3-ref').textContent = this.reference();
    // Historique.
    const tMax = Math.max(tStar, 10);
    const tMin = Math.max(0, tMax - 40);
    const w = this.hist.since(tMin * spt);
    const t = new Float64Array(w.t.length);
    for (let k = 0; k < t.length; k++) t[k] = w.t[k] / spt;
    const th = this.ctx.theme();
    drawHistory(
      $('s3-chart'),
      { t, series: [{ values: w.cx, color: th.accent }, { values: w.cz, color: th.lift }], tMin, tMax, markers: [], meanWindow: ready ? [start / spt, sim.time / spt] : null },
      th,
    );
    void sig3;
  }

  /** Repères de la littérature pour la configuration courante. */
  reference() {
    const { shape, re } = this.st;
    const fmt = this.ctx.fmt;
    if (shape === 'sphere') {
      const cd = (24 / re) * (1 + 0.15 * re ** 0.687) + 0.42 / (1 + 42500 * re ** -1.16);
      let regime = 'sillage stationnaire axisymétrique';
      if (re < 24) regime = 'écoulement attaché, sans recirculation';
      else if (re > 270) regime = 'lâcher de tourbillons en épingle (St ≈ 0,13 à Re = 300)';
      else if (re > 210) regime = 'sillage stationnaire oblique (portance latérale)';
      return `Référence sphère (Schiller–Naumann, Clift) : Cx ≈ ${fmt(cd, 2)} ; ${regime}. Le confinement latéral (veine 6D périodique) augmente Cx de quelques %.`;
    }
    if (shape === 'cylinder') {
      return re < 190
        ? 'Cylindre : sillage 2D (allée de von Kármán dès Re ≈ 47) ; Cx ≈ 1,3 à Re = 100.'
        : `Cylindre : au-delà de Re ≈ 190, le sillage devient tridimensionnel (mode A, puis B vers 260) ; Cx ≈ 1,2 mesuré, contre ≈ 1,4 en 2D.`;
    }
    if (shape === 'cube') return 'Cube face au vent : décollement fixé par les arêtes ; Cx ≈ 1,05 aux grands Reynolds, plus élevé aux faibles.';
    return 'Aile rectangulaire d’allongement 4 : tourbillons marginaux (plan transversal, ωx) ; à faible Reynolds, la couche limite épaisse réduit nettement la portance.';
  }

  dispose() {
    this.sim?.dispose();
    this.renderer?.dispose();
  }
}
