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

const MAX_PARTICLES = 9600;

/** Cadrages : [amont, aval, demi-hauteur] en multiples de L autour de la maquette. */
export const VIEWS = {
  near: { label: 'Maquette', box: [2, 6, 2] },
  wake: { label: 'Sillage', box: [3, 13, 4] },
  full: { label: 'Veine', box: null },
};

/**
 * Fronts de traceurs (lignes de temps) : FRONT_N particules réparties sur
 * FRONT_SPAN L de hauteur, lâchées ensemble à FRONT_X L en amont de la
 * maquette, un front tous les FRONT_DX L parcourus par l'écoulement amont.
 */
export const FRONT_N = 480;
export const FRONT_SPAN = 4;
export const FRONT_X = 3;
export const FRONT_DX = 0.5;

/**
 * Longueur (mailles) au-delà de laquelle un segment de front est rompu : il
 * enjamberait le corps ou la couche limite (particules voisines au départ
 * mais séparées par la paroi ou par le fort cisaillement pariétal).
 */
export function frontMaxSeg(L) {
  return Math.max(2, 0.25 * L);
}

/** Largeur minimale de la fenêtre (en L) : zoom maximal. */
export const ZOOM_MIN_W = 0.4;

/**
 * Caméra 2D, relative à la maquette et en unités de L : centre (dx, dy)
 * par rapport au point de référence, largeur w de la fenêtre.
 */
export function presetCamera(wt, mode) {
  const box = (VIEWS[mode] ?? VIEWS.near).box;
  const L = wt.L;
  if (!box) {
    return { dx: (wt.nx / 2 - (wt.x0 + 0.5)) / L, dy: (wt.ny / 2 - (wt.y0 + 0.5)) / L, w: wt.nx / L };
  }
  const [up, down] = box;
  return { dx: (down - up) / 2, dy: 0, w: up + down };
}

/**
 * Fenêtre affichée (mailles, réels) pour une caméra et un rapport hauteur /
 * largeur du canevas ; la fenêtre reste dans la veine. Renvoie aussi la
 * caméra effective (après butées).
 */
export function cameraWindow(wt, cam, aspect) {
  const L = wt.L;
  let ww = Math.min(wt.nx, Math.max(ZOOM_MIN_W * L, cam.w * L));
  let wh = ww * aspect;
  if (wh > wt.ny) {
    wh = wt.ny;
    ww = wh / aspect;
  }
  const cx = wt.x0 + 0.5 + cam.dx * L;
  const cy = wt.y0 + 0.5 + cam.dy * L;
  const wx = Math.max(0, Math.min(wt.nx - ww, cx - ww / 2));
  const wy = Math.max(0, Math.min(wt.ny - wh, cy - wh / 2));
  return {
    wx,
    wy,
    ww,
    wh,
    cam: { dx: (wx + ww / 2 - wt.x0 - 0.5) / L, dy: (wy + wh / 2 - wt.y0 - 0.5) / L, w: ww / L },
  };
}

/**
 * Caméra de largeur w (en L) qui place le point (px, py) du réseau à la
 * fraction (fx, fy) du canevas (depuis le coin haut gauche).
 */
export function anchoredCamera(wt, w, aspect, px, py, fx, fy) {
  const L = wt.L;
  const ww = w * L;
  const wh = ww * aspect;
  const cx = px - (fx - 0.5) * ww;
  const cy = py + (fy - 0.5) * wh;
  return cameraWindow(wt, { dx: (cx - wt.x0 - 0.5) / L, dy: (cy - wt.y0 - 0.5) / L, w }, aspect).cam;
}

/** Point du réseau sous la fraction (fx, fy) du canevas. */
export function windowPoint(win, fx, fy) {
  return [win.wx + fx * win.ww, win.wy + (1 - fy) * win.wh];
}

/** Ordonnées initiales des particules d'un front. */
function frontY(wt, k) {
  return wt.y0 - (FRONT_SPAN * wt.L) / 2 + (FRONT_SPAN * wt.L * (k + 0.5)) / FRONT_N;
}

/** Abscisse de lâcher des fronts. */
function frontX(wt) {
  return Math.max(2, wt.x0 - FRONT_X * wt.L);
}

/**
 * Segment de front à tracer ? Deux particules vivantes, pas trop éloignées,
 * et un segment qui ne traverse pas le corps (quatre points testés : un
 * segment qui enjambe un profil mince est rompu).
 */
function frontSegment(wt, type, x0, y0, x1, y1, maxSeg) {
  if (x0 < 0 || x1 < 0) return false;
  const dx = x1 - x0;
  const dy = y1 - y0;
  if (dx * dx + dy * dy > maxSeg * maxSeg) return false;
  for (let t = 0.2; t < 0.9; t += 0.2) {
    const mx = Math.round(x0 + t * dx);
    const my = Math.round(y0 + t * dy);
    if (mx >= 0 && my >= 0 && mx < wt.nx && my < wt.ny && type[my * wt.nx + mx] === SOLID) return false;
  }
  return true;
}

export { frontY, frontX };

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
    this.tracerMode = 'streak';
    this.frontSlot = 0;
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

  /** Cadrage prédéfini (remplace le zoom manuel). */
  setView(mode) {
    this.viewMode = VIEWS[mode] ? mode : 'near';
    this.cam = null;
  }

  /** Caméra manuelle (dx, dy, w en L), ou null pour le cadrage prédéfini. */
  setCamera(cam) {
    this.cam = cam;
  }

  /** Fenêtre courante (mailles) d'après la caméra et la forme du canevas. */
  window() {
    const wt = this.wt;
    const aspect = this.canvas.clientHeight / Math.max(1, this.canvas.clientWidth) || 0.5;
    const cam = this.cam ?? presetCamera(wt, this.viewMode);
    return cameraWindow(wt, cam, aspect);
  }

  /** 'streak' : peigne continu (lignes d'émission) ; 'fronts' : lignes de temps. */
  setTracerMode(mode) {
    if (mode === this.tracerMode) return;
    this.tracerMode = mode === 'fronts' ? 'fronts' : 'streak';
    this.clearParticles();
  }

  attach(wt) {
    this.wt = wt;
    this.off.width = wt.nx;
    this.off.height = wt.ny;
    this.image = this.offCtx.createImageData(wt.nx, wt.ny);
    this.pix = new Uint32Array(this.image.data.buffer);
    this.viewMode = this.viewMode ?? 'near';
    this.clearParticles();
    this.resize();
  }

  clearParticles() {
    this.pHead = 0;
    this.emitAcc = 0;
    this.frontSlot = 0;
    if (this.tracerMode === 'fronts') {
      // Emplacements fixes de FRONT_N particules, morts au départ.
      this.px.fill(-1);
      this.pCount = Math.floor(MAX_PARTICLES / FRONT_N) * FRONT_N;
    } else {
      this.pCount = 0;
    }
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
    this.emitAcc += steps * U;
    if (this.tracerMode === 'fronts') {
      // Un front complet tous les FRONT_DX L parcourus.
      const spacing = FRONT_DX * wt.L;
      const slots = Math.floor(MAX_PARTICLES / FRONT_N);
      const x0 = frontX(wt);
      while (this.emitAcc >= spacing) {
        this.emitAcc -= spacing;
        const base = this.frontSlot * FRONT_N;
        for (let k = 0; k < FRONT_N; k++) {
          this.px[base + k] = x0 + this.emitAcc;
          this.py[base + k] = frontY(wt, k);
        }
        this.frontSlot = (this.frontSlot + 1) % slots;
      }
    } else {
      // Émission : une particule tous les ~2.5 mailles parcourues, sur un peigne.
      const seeds = 17;
      const span = 5 * wt.L;
      const spacing = Math.max(2.5, wt.L / 6.4);
      while (this.emitAcc >= spacing) {
        this.emitAcc -= spacing;
        for (let s = 0; s < seeds; s++) {
          const y = wt.y0 - span / 2 + (span * (s + 0.5)) / seeds;
          this._add(1 + this.emitAcc, y);
        }
      }
    }
    const nsub = Math.max(1, Math.ceil((steps * U * 2) / 0.5));
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
    const win = this.window();
    this.win = win;
    const { wx, wy, ww, wh } = win;
    // Mailles couvrant la fenêtre (réelle) : bornes entières, avec une maille
    // de marge pour l'interpolation de l'image.
    const x0 = Math.max(0, Math.floor(wx) - 1);
    const x1 = Math.min(nx, Math.ceil(wx + ww) + 1);
    const y0 = Math.max(0, Math.floor(wy) - 1);
    const y1 = Math.min(ny, Math.ceil(wy + wh) + 1);

    // Mailles solides en bordure du corps : couleur du fluide voisin. Le
    // contour lisse de la maquette, tracé ensuite, recouvre l'intérieur ; ainsi
    // l'escalier des mailles ne dépasse pas du contour quand on zoome.
    const vals = this.vals && this.vals.length === lbm.N ? this.vals : (this.vals = new Float32Array(lbm.N));
    const solids = [];
    for (let y = y0; y < y1; y++) {
      const row = y * nx;
      const prow = (ny - 1 - y) * nx;
      for (let x = x0; x < x1; x++) {
        const idx = row + x;
        if (type[idx] === SOLID) {
          pix[prow + x] = body;
          vals[idx] = NaN;
          solids.push(idx);
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
        vals[idx] = v;
        let c = ((v - lo) * k) | 0;
        if (c < 0) c = 0;
        else if (c > n) c = n;
        pix[prow + x] = lut[c];
      }
    }
    for (const idx of solids) {
      const x = idx % nx;
      const y = (idx - x) / nx;
      let sum = 0;
      let cnt = 0;
      for (const j of [idx - 1, idx + 1, idx - nx, idx + nx]) {
        if (j < 0 || j >= lbm.N) continue;
        const v = vals[j];
        if (v === v) {
          sum += v;
          cnt++;
        }
      }
      if (!cnt) continue;
      let c = ((sum / cnt - lo) * k) | 0;
      if (c < 0) c = 0;
      else if (c > n) c = n;
      pix[(ny - 1 - y) * nx + x] = lut[c];
    }
    const srcY = ny - wy - wh; // ligne image du haut de la fenêtre
    this.offCtx.putImageData(this.image, 0, 0, x0, ny - y1, x1 - x0, y1 - y0);

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
    if (o.tracers && this.pCount && this.tracerMode === 'fronts') {
      const maxSeg = frontMaxSeg(L);
      const px = this.px;
      const py = this.py;
      const r = Math.max(1, this.dpr);
      ctx.lineWidth = Math.max(1, 1.3 * this.dpr);
      ctx.lineJoin = 'round';
      // Couleurs alternées d'un front à l'autre (même front suivi de part et
      // d'autre du profil).
      for (const parity of [0, 1]) {
        const color = parity ? this.colors.accent : this.colors.tracer;
        ctx.strokeStyle = color;
        ctx.fillStyle = color;
        const dots = [];
        ctx.beginPath();
        for (let base = parity * FRONT_N; base < this.pCount; base += 2 * FRONT_N) {
          let pen = false;
          for (let k = 0; k < FRONT_N - 1; k++) {
            const a = base + k;
            if (!frontSegment(wt, type, px[a], py[a], px[a + 1], py[a + 1], maxSeg)) {
              // Particule isolée : point.
              if (!pen && px[a] >= 0) dots.push(a);
              pen = false;
              continue;
            }
            if (!pen) ctx.moveTo(toX(px[a]), toY(py[a]));
            ctx.lineTo(toX(px[a + 1]), toY(py[a + 1]));
            pen = true;
          }
        }
        ctx.stroke();
        ctx.beginPath();
        for (const a of dots) ctx.rect(toX(px[a]) - r / 2, toY(py[a]) - r / 2, r, r);
        ctx.fill();
      }
    } else if (o.tracers && this.pCount) {
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

    drawOverlay(ctx, wt, win, W, H, this.dpr, this.colors, o.force);
  }
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
