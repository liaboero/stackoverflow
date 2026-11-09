/**
 * Rendu du champ (vorticité, vitesse, pression), de la maquette, des traceurs
 * et de la force résultante sur un <canvas>.
 */
import { SOLID } from './lbm.js';

/** Palettes : arrêts [position, '#rrggbb']. Diverging centrée sur 0.5. */
export const PALETTES = {
  light: {
    diverging: [
      [0, '#123f66'],
      [0.22, '#4f86b8'],
      [0.5, '#f3f6f7'],
      [0.78, '#de8456'],
      [1, '#8c2410'],
    ],
    sequential: [
      [0, '#f3f6f7'],
      [0.3, '#bcd8df'],
      [0.6, '#5f9fba'],
      [0.85, '#24557d'],
      [1, '#10233f'],
    ],
  },
  dark: {
    diverging: [
      [0, '#9fd3fb'],
      [0.22, '#2f74a8'],
      [0.5, '#0e161b'],
      [0.78, '#b4502a'],
      [1, '#ffc08f'],
    ],
    sequential: [
      [0, '#0e161b'],
      [0.3, '#173f52'],
      [0.6, '#2f8299'],
      [0.85, '#7cc6c6'],
      [1, '#e9f6d8'],
    ],
  },
};

export const FIELDS = {
  vort: { label: 'Vorticité ω·L/U', min: -5, max: 5, map: 'diverging' },
  speed: { label: 'Vitesse |u|/U', min: 0, max: 1.6, map: 'sequential' },
  cp: { label: 'Pression Cp', min: -1.5, max: 1.5, map: 'diverging' },
};

export function hexToRgb(h) {
  const n = parseInt(h.slice(1), 16);
  return [(n >> 16) & 255, (n >> 8) & 255, n & 255];
}

export function buildLut(stops, n = 256) {
  const lut = new Uint32Array(n);
  const rgb = stops.map(([p, c]) => [p, hexToRgb(c)]);
  for (let k = 0; k < n; k++) {
    const t = k / (n - 1);
    let s = 0;
    while (s < rgb.length - 2 && t > rgb[s + 1][0]) s++;
    const [p0, c0] = rgb[s];
    const [p1, c1] = rgb[s + 1];
    const u = Math.min(1, Math.max(0, (t - p0) / (p1 - p0)));
    const r = Math.round(c0[0] + (c1[0] - c0[0]) * u);
    const g = Math.round(c0[1] + (c1[1] - c0[1]) * u);
    const b = Math.round(c0[2] + (c1[2] - c0[2]) * u);
    lut[k] = (255 << 24) | (b << 16) | (g << 8) | r; // RGBA petit-boutiste
  }
  return lut;
}

export function gradientCss(theme, map) {
  const stops = PALETTES[theme][map];
  return `linear-gradient(90deg, ${stops.map(([p, c]) => `${c} ${(p * 100).toFixed(0)}%`).join(', ')})`;
}

const MAX_PARTICLES = 3200;

/** Cadrages : [amont, aval, demi-hauteur] en multiples de L autour de la maquette. */
export const VIEWS = {
  near: { label: 'Maquette', box: [2, 6, 2] },
  wake: { label: 'Sillage', box: [3, 13, 4] },
  full: { label: 'Veine', box: null },
};

export class FlowView {
  constructor(canvas) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.off = document.createElement('canvas');
    this.offCtx = this.off.getContext('2d');
    this.px = new Float32Array(MAX_PARTICLES);
    this.py = new Float32Array(MAX_PARTICLES);
    this.pCount = 0;
    this.pHead = 0;
    this.emitAcc = 0;
    this.theme = 'light';
    this.colors = { body: '#24313a', bodyEdge: '#0d161c', tracer: 'rgba(16,28,38,0.6)', accent: '#c8461c', lift: '#1f6f9c', ink: '#13202a' };
    this.setTheme('light', this.colors);
  }

  setTheme(theme, colors) {
    this.theme = theme;
    this.colors = { ...this.colors, ...colors };
    this.luts = {
      diverging: buildLut(PALETTES[theme].diverging),
      sequential: buildLut(PALETTES[theme].sequential),
    };
    const b = hexToRgb(this.colors.body.startsWith('#') ? this.colors.body : '#24313a');
    this.bodyPixel = (255 << 24) | (b[2] << 16) | (b[1] << 8) | b[0];
  }

  /** Fenêtre affichée (mailles) : x ∈ [wx, wx + ww[, y ∈ [wy, wy + wh[. */
  setView(mode) {
    this.viewMode = VIEWS[mode] ? mode : 'near';
    if (this.wt) this.win = computeWindow(this.wt, this.viewMode);
  }

  attach(wt) {
    this.wt = wt;
    this.off.width = wt.nx;
    this.off.height = wt.ny;
    this.image = this.offCtx.createImageData(wt.nx, wt.ny);
    this.pix = new Uint32Array(this.image.data.buffer);
    this.setView(this.viewMode ?? 'near');
    this.clearParticles();
    this.resize();
  }

  clearParticles() {
    this.pCount = 0;
    this.pHead = 0;
    this.emitAcc = 0;
  }

  resize() {
    const dpr = Math.min(window.devicePixelRatio || 1, 2.5);
    const w = Math.max(1, Math.round(this.canvas.clientWidth * dpr));
    const h = Math.max(1, Math.round(this.canvas.clientHeight * dpr));
    if (this.canvas.width !== w || this.canvas.height !== h) {
      this.canvas.width = w;
      this.canvas.height = h;
    }
    this.dpr = dpr;
  }

  /** Avance les traceurs de `steps` pas de temps (vitesses déjà calculées). */
  advect(steps) {
    const wt = this.wt;
    if (!wt || steps <= 0) return;
    const { nx, ny, lbm } = wt;
    const { ux, uy, type } = lbm;
    const U = wt.uIn;
    // Émission : une particule tous les ~2.5 mailles parcourues, sur un peigne.
    const seeds = 17;
    const span = 5 * wt.L;
    this.emitAcc += steps * U;
    const spacing = Math.max(2.5, wt.L / 6.4);
    while (this.emitAcc >= spacing) {
      this.emitAcc -= spacing;
      for (let s = 0; s < seeds; s++) {
        const y = wt.y0 - span / 2 + (span * (s + 0.5)) / seeds;
        this._add(1 + this.emitAcc, y);
      }
    }
    const nsub = Math.max(1, Math.ceil((steps * U * 2) / 0.8));
    const dt = steps / nsub;
    const px = this.px;
    const py = this.py;
    for (let k = 0; k < this.pCount; k++) {
      let x = px[k];
      let y = py[k];
      if (x < 0) continue;
      for (let s = 0; s < nsub; s++) {
        const i = Math.floor(x);
        const j = Math.floor(y);
        if (i < 0 || j < 0 || i >= nx - 1 || j >= ny - 1) {
          x = -1;
          break;
        }
        const fx = x - i;
        const fy = y - j;
        const a = j * nx + i;
        const b = a + 1;
        const c = a + nx;
        const d = c + 1;
        const vx = (ux[a] * (1 - fx) + ux[b] * fx) * (1 - fy) + (ux[c] * (1 - fx) + ux[d] * fx) * fy;
        const vy = (uy[a] * (1 - fx) + uy[b] * fx) * (1 - fy) + (uy[c] * (1 - fx) + uy[d] * fx) * fy;
        x += vx * dt;
        y += vy * dt;
        const ci = Math.round(x);
        const cj = Math.round(y);
        if (ci >= 0 && cj >= 0 && ci < nx && cj < ny && type[cj * nx + ci] === SOLID) {
          x = -1;
          break;
        }
      }
      px[k] = x;
      py[k] = y;
    }
  }

  _add(x, y) {
    const k = this.pHead;
    this.px[k] = x;
    this.py[k] = y;
    this.pHead = (k + 1) % MAX_PARTICLES;
    if (this.pCount < MAX_PARTICLES) this.pCount++;
  }

  /**
   * @param {object} o
   * @param {'vort'|'speed'|'cp'} o.field
   * @param {boolean} o.tracers
   * @param {{cx:number, cz:number}|null} o.force coefficient moyen à représenter
   */
  draw(o) {
    const wt = this.wt;
    if (!wt) return;
    const { nx, ny, lbm, L } = wt;
    const { ux, uy, rho, type } = lbm;
    const U = wt.uIn;
    const field = FIELDS[o.field] ?? FIELDS.vort;
    const lut = this.luts[field.map];
    const pix = this.pix;
    const n = lut.length - 1;
    const lo = field.min;
    const k = n / (field.max - field.min);
    const body = this.bodyPixel;
    const { wx, wy, ww, wh } = this.win;

    for (let y = wy; y < wy + wh; y++) {
      const row = y * nx;
      const prow = (ny - 1 - y) * nx;
      for (let x = wx; x < wx + ww; x++) {
        const idx = row + x;
        if (type[idx] === SOLID) {
          pix[prow + x] = body;
          continue;
        }
        let v;
        if (o.field === 'speed') {
          v = Math.sqrt(ux[idx] * ux[idx] + uy[idx] * uy[idx]) / U;
        } else if (o.field === 'cp') {
          v = (rho[idx] - 1) / 3 / (0.5 * U * U);
        } else {
          const xl = x > 0 ? idx - 1 : idx;
          const xr = x < nx - 1 ? idx + 1 : idx;
          const yd = y > 0 ? idx - nx : idx;
          const yu = y < ny - 1 ? idx + nx : idx;
          const dvdx = (uy[xr] - uy[xl]) / (xr - xl || 1);
          const dudy = (ux[yu] - ux[yd]) / ((yu - yd) / nx || 1);
          v = ((dvdx - dudy) * L) / U;
        }
        let c = ((v - lo) * k) | 0;
        if (c < 0) c = 0;
        else if (c > n) c = n;
        pix[prow + x] = lut[c];
      }
    }
    const srcY = ny - wy - wh; // ligne image du haut de la fenêtre
    this.offCtx.putImageData(this.image, 0, 0, wx, srcY, ww, wh);

    const ctx = this.ctx;
    const W = this.canvas.width;
    const H = this.canvas.height;
    const sx = W / ww;
    const sy = H / wh;
    ctx.imageSmoothingEnabled = true;
    ctx.imageSmoothingQuality = 'high';
    ctx.drawImage(this.off, wx, srcY, ww, wh, 0, 0, W, H);

    const toX = (x) => (x + 0.5 - wx) * sx;
    const toY = (y) => (ny - 0.5 - y - srcY) * sy;

    // Traceurs
    if (o.tracers && this.pCount) {
      ctx.fillStyle = this.colors.tracer;
      const r = Math.max(1, 1.1 * this.dpr);
      ctx.beginPath();
      for (let i = 0; i < this.pCount; i++) {
        const x = this.px[i];
        if (x < 0) continue;
        ctx.rect(toX(x) - r / 2, toY(this.py[i]) - r / 2, r, r);
      }
      ctx.fill();
    }

    drawOverlay(ctx, wt, this.win, W, H, this.dpr, this.colors, o.force);
  }
}

/** Fenêtre affichée (mailles) pour un cadrage : x ∈ [wx, wx + ww[, y ∈ [wy, wy + wh[. */
export function computeWindow(wt, mode) {
  const box = (VIEWS[mode] ?? VIEWS.near).box;
  if (!box) return { wx: 0, wy: 0, ww: wt.nx, wh: wt.ny };
  const [up, down, half] = box;
  const L = wt.L;
  const ww = Math.min(wt.nx, Math.round((up + down) * L));
  const wh = Math.min(wt.ny, Math.round(2 * half * L));
  const wx = Math.max(0, Math.min(wt.nx - ww, Math.round(wt.x0 - up * L)));
  const wy = Math.max(0, Math.min(wt.ny - wh, Math.round(wt.y0 + 0.5 - half * L)));
  return { wx, wy, ww, wh };
}

/**
 * Maquette (contour lisse), point de réduction des moments et force moyenne,
 * dessinés en 2D par-dessus le champ (rendu CPU ou GPU).
 */
export function drawOverlay(ctx, wt, win, W, H, dpr, colors, force) {
  const { wx, wy, ww, wh } = win;
  const { ny, L } = wt;
  const sx = W / ww;
  const sy = H / wh;
  const srcY = ny - wy - wh;
  const toX = (x) => (x + 0.5 - wx) * sx;
  const toY = (y) => (ny - 0.5 - y - srcY) * sy;
  {
    const o = { force };
    const body2 = wt.body;
    ctx.beginPath();
    body2.outline.forEach(([x, y], i) => (i ? ctx.lineTo(toX(x), toY(y)) : ctx.moveTo(toX(x), toY(y))));
    ctx.closePath();
    ctx.fillStyle = colors.body;
    ctx.fill();
    ctx.lineWidth = Math.max(1, dpr);
    ctx.strokeStyle = colors.bodyEdge;
    ctx.stroke();

    // Point de réduction des moments
    const rx = toX(body2.refX);
    const ry = toY(body2.refY);
    ctx.fillStyle = colors.bodyEdge;
    ctx.beginPath();
    ctx.arc(rx, ry, 2 * dpr, 0, Math.PI * 2);
    ctx.fill();

    // Force aérodynamique moyenne : 1 unité de coefficient ↔ 1,2 L
    if (o.force && Number.isFinite(o.force.cx) && Number.isFinite(o.force.cz)) {
      const scale = 1.2 * L;
      const ex = toX(body2.refX + o.force.cx * scale);
      const ey = toY(body2.refY + o.force.cz * scale);
      const lw = Math.max(1.5, 2 * dpr);
      // Composantes
      ctx.setLineDash([4 * dpr, 3 * dpr]);
      ctx.lineWidth = Math.max(1, dpr);
      ctx.strokeStyle = colors.accent;
      ctx.beginPath();
      ctx.moveTo(rx, ry);
      ctx.lineTo(ex, ry);
      ctx.stroke();
      ctx.strokeStyle = colors.lift;
      ctx.beginPath();
      ctx.moveTo(rx, ry);
      ctx.lineTo(rx, ey);
      ctx.stroke();
      ctx.setLineDash([]);
      arrow(ctx, rx, ry, ex, ey, lw, colors.ink, 7 * dpr);
    }
  }
}

function arrow(ctx, x0, y0, x1, y1, lw, color, head) {
  const dx = x1 - x0;
  const dy = y1 - y0;
  const len = Math.hypot(dx, dy);
  if (len < 2) return;
  const ux = dx / len;
  const uy = dy / len;
  const h = Math.min(head, len * 0.5);
  ctx.strokeStyle = color;
  ctx.fillStyle = color;
  ctx.lineWidth = lw;
  ctx.lineCap = 'round';
  ctx.beginPath();
  ctx.moveTo(x0, y0);
  ctx.lineTo(x1 - ux * h * 0.8, y1 - uy * h * 0.8);
  ctx.stroke();
  ctx.beginPath();
  ctx.moveTo(x1, y1);
  ctx.lineTo(x1 - ux * h - uy * h * 0.5, y1 - uy * h + ux * h * 0.5);
  ctx.lineTo(x1 - ux * h + uy * h * 0.5, y1 - uy * h - ux * h * 0.5);
  ctx.closePath();
  ctx.fill();
}
