/**
 * Rendu sur GPU (WebGL2) pour le moteur GpuLBM : champ (vorticité, vitesse,
 * pression) calculé et coloré par shaders, traceurs advectés et tracés sur
 * GPU. La maquette, le point de réduction et la force sont dessinés en 2D sur
 * un calque superposé, avec le même code que le rendu CPU.
 */
import { createProgram, createTexture, createFramebuffer } from './lbm-gpu.js';
import { FIELDS, PALETTES, VIEWS, buildLut, computeWindow, drawOverlay } from './render.js';

const PW = 128;
const PH = 128;
const NP = PW * PH;

const VS_QUAD = `#version 300 es
void main() {
  vec2 p = vec2(float((gl_VertexID << 1) & 2), float(gl_VertexID & 2));
  gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0);
}`;

const HEAD = `#version 300 es
precision highp float;
precision highp int;
precision highp sampler2D;
uniform ivec2 uSize;
uniform sampler2D uC;
uniform sampler2D uGeo;
int cellType(ivec2 q) {
  return int(texelFetch(uGeo, q, 0).r * 255.0 + 0.5);
}
`;

// Grandeurs macroscopiques (0, ρ, ux, uy) à partir des populations, une fois par image.
const FS_MACRO = `#version 300 es
precision highp float;
precision highp int;
precision highp sampler2D;
uniform sampler2D uA;
uniform sampler2D uB;
uniform sampler2D uF8;
out vec4 o;
void main() {
  ivec2 p = ivec2(gl_FragCoord.xy);
  vec4 a = texelFetch(uA, p, 0);
  vec4 b = texelFetch(uB, p, 0);
  float f8 = texelFetch(uF8, p, 0).r;
  float rho = a.x + a.y + a.z + a.w + b.x + b.y + b.z + b.w + f8;
  float ux = (a.y - a.w + b.y - b.z - b.w + f8) / rho;
  float uy = (a.z - b.x + b.y + b.z - b.w - f8) / rho;
  o = vec4(0.0, rho, ux, uy);
}
`;

// Grandeur affichée, à la résolution du réseau (sentinelle −1e30 : solide).
const FS_SCALAR = `${HEAD}
uniform int uField;
uniform float uU;
uniform float uL;
out vec4 o;
void main() {
  ivec2 p = ivec2(gl_FragCoord.xy);
  if (cellType(p) == 1) {
    o = vec4(-1e30);
    return;
  }
  vec4 c = texelFetch(uC, p, 0);
  float v;
  if (uField == 1) {
    v = length(c.zw) / uU;
  } else if (uField == 2) {
    v = (c.y - 1.0) / 3.0 / (0.5 * uU * uU);
  } else {
    ivec2 xl = ivec2(max(p.x - 1, 0), p.y);
    ivec2 xr = ivec2(min(p.x + 1, uSize.x - 1), p.y);
    ivec2 yd = ivec2(p.x, max(p.y - 1, 0));
    ivec2 yu = ivec2(p.x, min(p.y + 1, uSize.y - 1));
    float dvdx = (texelFetch(uC, xr, 0).w - texelFetch(uC, xl, 0).w) / float(max(1, xr.x - xl.x));
    float dudy = (texelFetch(uC, yu, 0).z - texelFetch(uC, yd, 0).z) / float(max(1, yu.y - yd.y));
    v = (dvdx - dudy) * uL / uU;
  }
  o = vec4(v, 0.0, 0.0, 1.0);
}
`;

// Affichage : interpolation bilinéaire (sans les mailles solides) + palette.
const FS_DISPLAY = `#version 300 es
precision highp float;
precision highp int;
precision highp sampler2D;
uniform sampler2D uS;
uniform sampler2D uLut;
uniform ivec2 uSize;
uniform vec4 uWin;
uniform vec2 uCanvas;
uniform float uMin;
uniform float uMax;
uniform vec4 uBody;
out vec4 o;
void main() {
  vec2 cell = uWin.xy + gl_FragCoord.xy / uCanvas * uWin.zw - 0.5;
  ivec2 i0 = ivec2(floor(cell));
  vec2 fr = cell - vec2(i0);
  float sum = 0.0;
  float wsum = 0.0;
  for (int k = 0; k < 4; k++) {
    ivec2 d = ivec2(k & 1, k >> 1);
    ivec2 q = clamp(i0 + d, ivec2(0), uSize - 1);
    float w = (d.x == 1 ? fr.x : 1.0 - fr.x) * (d.y == 1 ? fr.y : 1.0 - fr.y);
    float v = texelFetch(uS, q, 0).r;
    if (v > -1e29) {
      sum += w * v;
      wsum += w;
    }
  }
  if (wsum < 0.5) {
    o = uBody;
    return;
  }
  float t = clamp((sum / wsum - uMin) / (uMax - uMin), 0.0, 1.0);
  o = vec4(texture(uLut, vec2(t * (255.0 / 256.0) + 0.5 / 256.0, 0.5)).rgb, 1.0);
}
`;

// Advection des traceurs (sous-pas, vitesse interpolée bilinéairement).
const FS_PUPDATE = `${HEAD}
uniform sampler2D uP;
uniform float uDt;
uniform int uSub;
out vec4 o;
vec2 vel(vec2 pos) {
  ivec2 i = ivec2(floor(pos));
  vec2 f = pos - vec2(i);
  vec2 a = texelFetch(uC, i, 0).zw;
  vec2 b = texelFetch(uC, i + ivec2(1, 0), 0).zw;
  vec2 c = texelFetch(uC, i + ivec2(0, 1), 0).zw;
  vec2 d = texelFetch(uC, i + ivec2(1, 1), 0).zw;
  return mix(mix(a, b, f.x), mix(c, d, f.x), f.y);
}
void main() {
  vec4 s = texelFetch(uP, ivec2(gl_FragCoord.xy), 0);
  if (s.z < 0.5) {
    o = s;
    return;
  }
  vec2 pos = s.xy;
  for (int k = 0; k < uSub; k++) {
    if (pos.x < 0.0 || pos.y < 0.0 || pos.x >= float(uSize.x - 1) || pos.y >= float(uSize.y - 1)) {
      o = vec4(-1.0, -1.0, 0.0, 0.0);
      return;
    }
    pos += vel(pos) * uDt;
    ivec2 r = ivec2(floor(pos + 0.5));
    if (r.x >= 0 && r.y >= 0 && r.x < uSize.x && r.y < uSize.y && cellType(r) == 1) {
      o = vec4(-1.0, -1.0, 0.0, 0.0);
      return;
    }
  }
  o = vec4(pos, 1.0, 0.0);
}
`;

const VS_POINTS = `#version 300 es
precision highp float;
precision highp int;
precision highp sampler2D;
uniform sampler2D uP;
uniform vec4 uWin;
uniform float uPointSize;
void main() {
  vec4 s = texelFetch(uP, ivec2(gl_VertexID % ${PW}, gl_VertexID / ${PW}), 0);
  if (s.z < 0.5) {
    gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
    gl_PointSize = 0.0;
    return;
  }
  vec2 c = (s.xy + 0.5 - uWin.xy) / uWin.zw * 2.0 - 1.0;
  gl_Position = vec4(c, 0.0, 1.0);
  gl_PointSize = uPointSize;
}`;

const FS_POINTS = `#version 300 es
precision highp float;
uniform vec4 uColor;
out vec4 o;
void main() {
  o = uColor;
}`;

const FIELD_INDEX = { vort: 0, speed: 1, cp: 2 };

let colorCtx = null;
/** Couleur CSS → [r, g, b, a] dans [0, 1]. */
function parseColor(str) {
  if (!colorCtx) {
    const c = document.createElement('canvas');
    c.width = c.height = 1;
    colorCtx = c.getContext('2d', { willReadFrequently: true });
  }
  colorCtx.clearRect(0, 0, 1, 1);
  colorCtx.fillStyle = '#000';
  colorCtx.fillStyle = str;
  colorCtx.fillRect(0, 0, 1, 1);
  const d = colorCtx.getImageData(0, 0, 1, 1).data;
  return [d[0] / 255, d[1] / 255, d[2] / 255, d[3] / 255];
}

export class GpuFlowView {
  /**
   * @param {HTMLCanvasElement} glCanvas canevas WebGL2 (celui du solveur)
   * @param {HTMLCanvasElement} overlay canevas 2D superposé
   * @param {WebGL2RenderingContext} gl
   */
  constructor(glCanvas, overlay, gl) {
    this.glCanvas = glCanvas;
    this.canvas = overlay;
    this.ctx = overlay.getContext('2d');
    this.gl = gl;
    this.vao = gl.createVertexArray();
    this.progMacro = createProgram(gl, VS_QUAD, FS_MACRO);
    this.progScalar = createProgram(gl, VS_QUAD, FS_SCALAR);
    this.progDisplay = createProgram(gl, VS_QUAD, FS_DISPLAY);
    this.progUpdate = createProgram(gl, VS_QUAD, FS_PUPDATE);
    this.progPoints = createProgram(gl, VS_POINTS, FS_POINTS);
    this.lut = {
      diverging: createTexture(gl, 256, 1, gl.RGBA8, gl.LINEAR),
      sequential: createTexture(gl, 256, 1, gl.RGBA8, gl.LINEAR),
    };
    this.pTex = [createTexture(gl, PW, PH, gl.RGBA32F), createTexture(gl, PW, PH, gl.RGBA32F)];
    this.pFbo = this.pTex.map((t) => createFramebuffer(gl, [t]));
    this.pCur = 0;
    this.pHead = 0;
    this.emitAcc = 0;
    this.scalarTex = null;
    this.scalarFbo = null;
    this.macroTex = null;
    this.macroFbo = null;
    this._macroTime = -1;
    this.size = [0, 0];
    this.viewMode = 'near';
    this.dpr = 1;
    this.colors = { body: '#24313a', bodyEdge: '#0d161c', tracer: 'rgba(16,28,38,0.6)', accent: '#c8461c', lift: '#1f6f9c', ink: '#13202a' };
    this.setTheme('light', this.colors);
    this.clearParticles();
  }

  setTheme(theme, colors) {
    const gl = this.gl;
    this.theme = theme;
    this.colors = { ...this.colors, ...colors };
    gl.activeTexture(gl.TEXTURE0);
    for (const map of ['diverging', 'sequential']) {
      const lut = buildLut(PALETTES[theme][map]);
      gl.bindTexture(gl.TEXTURE_2D, this.lut[map]);
      gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
      gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, 256, 1, gl.RGBA, gl.UNSIGNED_BYTE, new Uint8Array(lut.buffer));
    }
    this.bodyRgba = parseColor(this.colors.body);
    this.tracerRgba = parseColor(this.colors.tracer);
  }

  setView(mode) {
    this.viewMode = VIEWS[mode] ? mode : 'near';
    if (this.wt) this.win = computeWindow(this.wt, this.viewMode);
  }

  attach(wt) {
    const gl = this.gl;
    this.wt = wt;
    if (this.size[0] !== wt.nx || this.size[1] !== wt.ny) {
      if (this.scalarTex) {
        gl.deleteTexture(this.scalarTex);
        gl.deleteFramebuffer(this.scalarFbo);
        gl.deleteTexture(this.macroTex);
        gl.deleteFramebuffer(this.macroFbo);
      }
      this.scalarTex = createTexture(gl, wt.nx, wt.ny, gl.R32F);
      this.scalarFbo = createFramebuffer(gl, [this.scalarTex]);
      this.macroTex = createTexture(gl, wt.nx, wt.ny, gl.RGBA32F);
      this.macroFbo = createFramebuffer(gl, [this.macroTex]);
      this.size = [wt.nx, wt.ny];
    }
    this._macroTime = -1;
    this.setView(this.viewMode);
    this.clearParticles();
    this.resize();
  }

  /** ρ, u sur tout le domaine, recalculés si le calcul a avancé. */
  _macros() {
    const sim = this.wt.lbm;
    if (this._macroTime === sim.time) return;
    this._macroTime = sim.time;
    const gl = this.gl;
    const P = this.progMacro;
    const st = sim.state;
    gl.bindVertexArray(this.vao);
    gl.useProgram(P.p);
    ['uA', 'uB', 'uF8'].forEach((name, k) => {
      gl.activeTexture(gl.TEXTURE0 + k);
      gl.bindTexture(gl.TEXTURE_2D, st[k]);
      gl.uniform1i(P.loc[name], k);
    });
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.macroFbo);
    gl.viewport(0, 0, this.wt.nx, this.wt.ny);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
  }

  clearParticles() {
    const gl = this.gl;
    const zero = new Float32Array(4 * NP);
    gl.activeTexture(gl.TEXTURE0);
    for (const t of this.pTex) {
      gl.bindTexture(gl.TEXTURE_2D, t);
      gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, PW, PH, gl.RGBA, gl.FLOAT, zero);
    }
    this.pHead = 0;
    this.emitAcc = 0;
  }

  resize() {
    const dpr = Math.min(window.devicePixelRatio || 1, 2.5);
    for (const c of [this.glCanvas, this.canvas]) {
      const w = Math.max(1, Math.round(c.clientWidth * dpr));
      const h = Math.max(1, Math.round(c.clientHeight * dpr));
      if (c.width !== w || c.height !== h) {
        c.width = w;
        c.height = h;
      }
    }
    this.dpr = dpr;
  }

  /** Émet de nouveaux traceurs puis les advecte de `steps` pas de temps. */
  advect(steps) {
    const wt = this.wt;
    if (!wt || steps <= 0) return;
    const gl = this.gl;
    const U = wt.uIn;
    const seeds = 17;
    const span = 5 * wt.L;
    const spacing = Math.max(2.5, wt.L / 6.4);
    this.emitAcc += steps * U;
    const fresh = [];
    while (this.emitAcc >= spacing && fresh.length < 2 * NP) {
      this.emitAcc -= spacing;
      for (let s = 0; s < seeds; s++) fresh.push(1 + this.emitAcc, wt.y0 - span / 2 + (span * (s + 0.5)) / seeds);
    }
    this.emitAcc = Math.min(this.emitAcc, spacing);
    const n = fresh.length / 2;
    if (n) {
      // Écriture dans l'anneau, par segments de lignes contigus.
      const data = new Float32Array(4 * n);
      for (let k = 0; k < n; k++) {
        data[4 * k] = fresh[2 * k];
        data[4 * k + 1] = fresh[2 * k + 1];
        data[4 * k + 2] = 1;
      }
      gl.activeTexture(gl.TEXTURE0);
      gl.bindTexture(gl.TEXTURE_2D, this.pTex[this.pCur]);
      let k = 0;
      while (k < n) {
        const idx = (this.pHead + k) % NP;
        const x = idx % PW;
        const y = Math.floor(idx / PW);
        const len = Math.min(PW - x, n - k);
        gl.texSubImage2D(gl.TEXTURE_2D, 0, x, y, len, 1, gl.RGBA, gl.FLOAT, data.subarray(4 * k, 4 * (k + len)));
        k += len;
      }
      this.pHead = (this.pHead + n) % NP;
    }

    this._macros();
    const nsub = Math.min(64, Math.max(1, Math.ceil((steps * U * 2) / 0.8)));
    const P = this.progUpdate;
    gl.bindVertexArray(this.vao);
    gl.useProgram(P.p);
    this._bindSim(P);
    gl.activeTexture(gl.TEXTURE2);
    gl.bindTexture(gl.TEXTURE_2D, this.pTex[this.pCur]);
    gl.uniform1i(P.loc.uP, 2);
    gl.uniform1f(P.loc.uDt, steps / nsub);
    gl.uniform1i(P.loc.uSub, nsub);
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.pFbo[1 - this.pCur]);
    gl.viewport(0, 0, PW, PH);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    this.pCur = 1 - this.pCur;
  }

  _bindSim(P) {
    const gl = this.gl;
    const sim = this.wt.lbm;
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, this.macroTex);
    gl.uniform1i(P.loc.uC, 0);
    gl.activeTexture(gl.TEXTURE1);
    gl.bindTexture(gl.TEXTURE_2D, sim.geoTex);
    gl.uniform1i(P.loc.uGeo, 1);
    gl.uniform2i(P.loc.uSize, this.wt.nx, this.wt.ny);
  }

  /**
   * @param {object} o
   * @param {'vort'|'speed'|'cp'} o.field
   * @param {boolean} o.tracers
   * @param {{cx:number, cz:number}|null} o.force
   */
  draw(o) {
    const wt = this.wt;
    if (!wt) return;
    const gl = this.gl;
    const field = FIELDS[o.field] ?? FIELDS.vort;
    const { wx, wy, ww, wh } = this.win;
    this._macros();
    gl.bindVertexArray(this.vao);

    // 1. Grandeur à la résolution du réseau
    let P = this.progScalar;
    gl.useProgram(P.p);
    this._bindSim(P);
    gl.uniform1i(P.loc.uField, FIELD_INDEX[o.field] ?? 0);
    gl.uniform1f(P.loc.uU, wt.uIn);
    gl.uniform1f(P.loc.uL, wt.L);
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.scalarFbo);
    gl.viewport(0, 0, wt.nx, wt.ny);
    gl.drawArrays(gl.TRIANGLES, 0, 3);

    // 2. Affichage dans la fenêtre
    const W = this.glCanvas.width;
    const H = this.glCanvas.height;
    P = this.progDisplay;
    gl.useProgram(P.p);
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, this.scalarTex);
    gl.uniform1i(P.loc.uS, 0);
    gl.activeTexture(gl.TEXTURE1);
    gl.bindTexture(gl.TEXTURE_2D, this.lut[field.map]);
    gl.uniform1i(P.loc.uLut, 1);
    gl.uniform2i(P.loc.uSize, wt.nx, wt.ny);
    gl.uniform4f(P.loc.uWin, wx, wy, ww, wh);
    gl.uniform2f(P.loc.uCanvas, W, H);
    gl.uniform1f(P.loc.uMin, field.min);
    gl.uniform1f(P.loc.uMax, field.max);
    gl.uniform4f(P.loc.uBody, ...this.bodyRgba);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.viewport(0, 0, W, H);
    gl.drawArrays(gl.TRIANGLES, 0, 3);

    // 3. Traceurs
    if (o.tracers) {
      P = this.progPoints;
      gl.useProgram(P.p);
      gl.activeTexture(gl.TEXTURE0);
      gl.bindTexture(gl.TEXTURE_2D, this.pTex[this.pCur]);
      gl.uniform1i(P.loc.uP, 0);
      gl.uniform4f(P.loc.uWin, wx, wy, ww, wh);
      gl.uniform1f(P.loc.uPointSize, Math.max(1, 1.2 * this.dpr));
      const [r, g, b, a] = this.tracerRgba;
      gl.uniform4f(P.loc.uColor, r * a, g * a, b * a, a);
      gl.enable(gl.BLEND);
      gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
      gl.drawArrays(gl.POINTS, 0, NP);
      gl.disable(gl.BLEND);
    }

    // 4. Calque 2D : maquette, point de réduction, force
    const ctx = this.ctx;
    ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
    drawOverlay(ctx, wt, this.win, this.canvas.width, this.canvas.height, this.dpr, this.colors, o.force);
  }

  dispose() {
    const gl = this.gl;
    [this.scalarTex, this.macroTex, this.lut.diverging, this.lut.sequential, ...this.pTex].forEach((t) => t && gl.deleteTexture(t));
    [this.scalarFbo, this.macroFbo, ...this.pFbo].forEach((f) => f && gl.deleteFramebuffer(f));
    [this.progMacro, this.progScalar, this.progDisplay, this.progUpdate, this.progPoints].forEach((p) => gl.deleteProgram(p.p));
    gl.deleteVertexArray(this.vao);
    this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
  }
}
