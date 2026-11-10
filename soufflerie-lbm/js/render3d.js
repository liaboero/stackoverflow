/**
 * Rendu WebGPU de la soufflerie 3D : maquette éclairée, plan de coupe
 * (vitesse ou vorticité), isosurfaces du critère Q par lancer de rayons à
 * travers le champ, cadre de la veine. Coordonnées monde = mailles du réseau.
 */

const VIEW_BYTES = 256;

const COMMON = /* wgsl */ `
struct View {
  viewProj: mat4x4<f32>,
  invViewProj: mat4x4<f32>,
  eye: vec4<f32>,
  dims: vec4<f32>,
  bg: vec4<f32>,
  bodyCol: vec4<f32>,
  lineCol: vec4<f32>,
  plane: vec4<f32>,
  iso: vec4<f32>,
  viewport: vec4<f32>,
};
@group(0) @binding(0) var<uniform> V: View;

fn seqMap(t0: f32) -> vec3<f32> {
  let t = clamp(t0, 0.0, 1.0);
  let a = vec3<f32>(0.05, 0.11, 0.27);
  let b = vec3<f32>(0.10, 0.45, 0.62);
  let c = vec3<f32>(0.42, 0.78, 0.62);
  let d = vec3<f32>(0.99, 0.88, 0.35);
  if (t < 0.33) { return mix(a, b, t / 0.33); }
  if (t < 0.66) { return mix(b, c, (t - 0.33) / 0.33); }
  return mix(c, d, (t - 0.66) / 0.34);
}

fn divMap(t0: f32) -> vec3<f32> {
  let t = clamp(t0, -1.0, 1.0);
  let mid = vec3<f32>(0.95, 0.95, 0.93);
  if (t < 0.0) { return mix(mid, vec3<f32>(0.13, 0.36, 0.70), -t); }
  return mix(mid, vec3<f32>(0.78, 0.20, 0.15), t);
}
`;

const DERIVE = /* wgsl */ `
@group(0) @binding(0) var macroTex: texture_3d<f32>;
@group(0) @binding(1) var derivTex: texture_storage_3d<rgba16float, write>;
fn U(p: vec3<i32>, dm: vec3<i32>) -> vec3<f32> {
  return textureLoad(macroTex, clamp(p, vec3<i32>(0), dm), 0).xyz;
}
@compute @workgroup_size(4, 4, 4)
fn main(@builtin(global_invocation_id) g: vec3<u32>) {
  let d = textureDimensions(macroTex);
  if (g.x >= d.x || g.y >= d.y || g.z >= d.z) { return; }
  let p = vec3<i32>(g);
  let dm = vec3<i32>(d) - vec3<i32>(1);
  let dx = 0.5 * (U(p + vec3<i32>(1, 0, 0), dm) - U(p - vec3<i32>(1, 0, 0), dm));
  let dy = 0.5 * (U(p + vec3<i32>(0, 1, 0), dm) - U(p - vec3<i32>(0, 1, 0), dm));
  let dz = 0.5 * (U(p + vec3<i32>(0, 0, 1), dm) - U(p - vec3<i32>(0, 0, 1), dm));
  // dx = ∂u/∂x (vecteur), etc.
  let w = vec3<f32>(dy.z - dz.y, dz.x - dx.z, dx.y - dy.x);
  let s2 = dx.x * dx.x + dy.y * dy.y + dz.z * dz.z
    + 0.5 * ((dx.y + dy.x) * (dx.y + dy.x) + (dx.z + dz.x) * (dx.z + dz.x) + (dy.z + dz.y) * (dy.z + dz.y));
  var q = 0.5 * (0.5 * dot(w, w) - s2);
  // Au contact du corps, les différences finies traversent la paroi (u = 0
  // dans le solide) : Q y est artificiel, on l'annule.
  let r = textureLoad(macroTex, p, 0).w;
  var nearSolid = r < 0.5;
  for (var k = 0; k < 3; k++) {
    var e = vec3<i32>(0);
    e[k] = 1;
    if (textureLoad(macroTex, clamp(p + e, vec3<i32>(0), dm), 0).w < 0.5 || textureLoad(macroTex, clamp(p - e, vec3<i32>(0), dm), 0).w < 0.5) {
      nearSolid = true;
    }
  }
  if (nearSolid) { q = 0.0; }
  textureStore(derivTex, g, vec4<f32>(w, q));
}
`;

const BODY = /* wgsl */ `
${COMMON}
struct VOut { @builtin(position) pos: vec4<f32>, @location(0) n: vec3<f32>, @location(1) w: vec3<f32> };
@vertex fn vs(@location(0) p: vec3<f32>, @location(1) n: vec3<f32>) -> VOut {
  var o: VOut;
  o.pos = V.viewProj * vec4<f32>(p, 1.0);
  o.n = n;
  o.w = p;
  return o;
}
@fragment fn fs(i: VOut) -> @location(0) vec4<f32> {
  let n = normalize(i.n);
  let l = normalize(V.eye.xyz - i.w);
  let lam = abs(dot(n, l));
  let key = max(0.0, dot(n, normalize(vec3<f32>(-0.4, -0.3, 0.85))));
  let c = V.bodyCol.rgb * (0.35 + 0.45 * lam + 0.25 * key);
  return vec4<f32>(c, 1.0);
}
`;

const LINES = /* wgsl */ `
${COMMON}
@vertex fn vs(@location(0) p: vec3<f32>) -> @builtin(position) vec4<f32> {
  return V.viewProj * vec4<f32>(p, 1.0);
}
@fragment fn fs() -> @location(0) vec4<f32> {
  return vec4<f32>(V.lineCol.rgb, 1.0);
}
`;

const SLICE = /* wgsl */ `
${COMMON}
@group(0) @binding(1) var macroTex: texture_3d<f32>;
@group(0) @binding(2) var derivTex: texture_3d<f32>;
@group(0) @binding(3) var samp: sampler;
struct VOut { @builtin(position) pos: vec4<f32>, @location(0) w: vec3<f32> };
@vertex fn vs(@builtin(vertex_index) vi: u32) -> VOut {
  let corners = array<vec2<f32>, 6>(vec2<f32>(0.0, 0.0), vec2<f32>(1.0, 0.0), vec2<f32>(1.0, 1.0), vec2<f32>(0.0, 0.0), vec2<f32>(1.0, 1.0), vec2<f32>(0.0, 1.0));
  let c = corners[vi];
  let d = V.dims.xyz - vec3<f32>(1.0);
  var p: vec3<f32>;
  let axis = i32(V.plane.x);
  if (axis == 0) { p = vec3<f32>(V.plane.y, c.x * d.y, c.y * d.z); }
  else if (axis == 1) { p = vec3<f32>(c.x * d.x, V.plane.y, c.y * d.z); }
  else { p = vec3<f32>(c.x * d.x, c.y * d.y, V.plane.y); }
  var o: VOut;
  o.pos = V.viewProj * vec4<f32>(p, 1.0);
  o.w = p;
  return o;
}
@fragment fn fs(i: VOut) -> @location(0) vec4<f32> {
  let uvw = (i.w + vec3<f32>(0.5)) / V.dims.xyz;
  let m = textureSampleLevel(macroTex, samp, uvw, 0.0);
  let dv = textureSampleLevel(derivTex, samp, uvw, 0.0);
  if (m.w < 0.5) { return vec4<f32>(V.bodyCol.rgb * 0.6, 1.0); }
  let U0 = V.plane.w;
  let L = V.dims.w;
  var col: vec3<f32>;
  if (V.plane.z < 0.5) {
    col = seqMap(length(m.xyz) / (1.6 * U0));
  } else {
    let axis = i32(V.plane.x);
    var w = dv.y;
    if (axis == 0) { w = dv.x; } else if (axis == 2) { w = dv.z; }
    col = divMap(w * L / U0 / V.iso.w);
  }
  return vec4<f32>(col, 1.0);
}
`;

const ISO = /* wgsl */ `
${COMMON}
@group(0) @binding(1) var derivTex: texture_3d<f32>;
@group(0) @binding(2) var samp: sampler;
@group(0) @binding(3) var depthTex: texture_depth_2d;
@group(0) @binding(4) var macroTex: texture_3d<f32>;
@vertex fn vs(@builtin(vertex_index) vi: u32) -> @builtin(position) vec4<f32> {
  let p = array<vec2<f32>, 3>(vec2<f32>(-1.0, -1.0), vec2<f32>(3.0, -1.0), vec2<f32>(-1.0, 3.0));
  return vec4<f32>(p[vi], 0.0, 1.0);
}
fn Q(p: vec3<f32>) -> f32 {
  return textureSampleLevel(derivTex, samp, (p + vec3<f32>(0.5)) / V.dims.xyz, 0.0).w;
}
fn unproj(ndc: vec3<f32>) -> vec3<f32> {
  let h = V.invViewProj * vec4<f32>(ndc, 1.0);
  return h.xyz / h.w;
}
@fragment fn fs(@builtin(position) fc: vec4<f32>) -> @location(0) vec4<f32> {
  let ndc = vec2<f32>(fc.x / V.viewport.x * 2.0 - 1.0, 1.0 - fc.y / V.viewport.y * 2.0);
  let pn = unproj(vec3<f32>(ndc, 0.0));
  let pf = unproj(vec3<f32>(ndc, 1.0));
  let o = V.eye.xyz;
  let dir = normalize(pf - pn);
  let lo = vec3<f32>(0.0);
  let hi = V.dims.xyz - vec3<f32>(1.0);
  let inv = 1.0 / dir;
  let t0v = (lo - o) * inv;
  let t1v = (hi - o) * inv;
  let tmin = max(max(min(t0v.x, t1v.x), min(t0v.y, t1v.y)), min(t0v.z, t1v.z));
  let tmax0 = min(min(max(t0v.x, t1v.x), max(t0v.y, t1v.y)), max(t0v.z, t1v.z));
  let dpt = textureLoad(depthTex, vec2<i32>(fc.xy), 0);
  var tmax = tmax0;
  if (dpt < 1.0) { tmax = min(tmax, length(unproj(vec3<f32>(ndc, dpt)) - o)); }
  var t = max(tmin, 0.0);
  if (t >= tmax) { discard; }
  let level = V.iso.x;
  let stp = V.iso.z;
  var hit = false;
  var tp = t;
  for (var k = 0; k < 1200; k++) {
    if (t > tmax) { break; }
    if (Q(o + dir * t) > level) { hit = true; break; }
    tp = t;
    t += stp;
  }
  if (!hit) { discard; }
  var a = tp;
  var b = t;
  for (var k = 0; k < 5; k++) {
    let m = 0.5 * (a + b);
    if (Q(o + dir * m) > level) { b = m; } else { a = m; }
  }
  let p = o + dir * b;
  let e = 1.0;
  let gq = vec3<f32>(
    Q(p + vec3<f32>(e, 0.0, 0.0)) - Q(p - vec3<f32>(e, 0.0, 0.0)),
    Q(p + vec3<f32>(0.0, e, 0.0)) - Q(p - vec3<f32>(0.0, e, 0.0)),
    Q(p + vec3<f32>(0.0, 0.0, e)) - Q(p - vec3<f32>(0.0, 0.0, e)));
  let n = -normalize(gq + vec3<f32>(1e-9));
  let lam = abs(dot(n, -dir));
  // Teinte : vitesse locale (lent = bleu, rapide = jaune).
  let u = textureSampleLevel(macroTex, samp, (p + vec3<f32>(0.5)) / V.dims.xyz, 0.0).xyz;
  let base = seqMap(0.15 + 0.75 * length(u) / (1.3 * V.plane.w));
  let c = base * (0.35 + 0.75 * lam) + vec3<f32>(0.25) * pow(lam, 24.0);
  return vec4<f32>(c, 0.88);
}
`;

// ------------------------------------------------------------ matrices (colonnes)

function perspective(fovy, aspect, near, far) {
  const f = 1 / Math.tan(fovy / 2);
  return new Float32Array([f / aspect, 0, 0, 0, 0, f, 0, 0, 0, 0, far / (near - far), -1, 0, 0, (far * near) / (near - far), 0]);
}

function lookAt(eye, at, up) {
  const z = norm3(sub3(eye, at));
  const x = norm3(cross3(up, z));
  const y = cross3(z, x);
  return new Float32Array([x[0], y[0], z[0], 0, x[1], y[1], z[1], 0, x[2], y[2], z[2], 0, -dot3(x, eye), -dot3(y, eye), -dot3(z, eye), 1]);
}

const sub3 = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const dot3 = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross3 = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const norm3 = (a) => {
  const l = Math.hypot(a[0], a[1], a[2]) || 1;
  return [a[0] / l, a[1] / l, a[2] / l];
};

function mul4(a, b) {
  const o = new Float32Array(16);
  for (let c = 0; c < 4; c++)
    for (let r = 0; r < 4; r++) {
      let s = 0;
      for (let k = 0; k < 4; k++) s += a[k * 4 + r] * b[c * 4 + k];
      o[c * 4 + r] = s;
    }
  return o;
}

function invert4(m) {
  const inv = new Float32Array(16);
  const a = m;
  inv[0] = a[5] * a[10] * a[15] - a[5] * a[11] * a[14] - a[9] * a[6] * a[15] + a[9] * a[7] * a[14] + a[13] * a[6] * a[11] - a[13] * a[7] * a[10];
  inv[4] = -a[4] * a[10] * a[15] + a[4] * a[11] * a[14] + a[8] * a[6] * a[15] - a[8] * a[7] * a[14] - a[12] * a[6] * a[11] + a[12] * a[7] * a[10];
  inv[8] = a[4] * a[9] * a[15] - a[4] * a[11] * a[13] - a[8] * a[5] * a[15] + a[8] * a[7] * a[13] + a[12] * a[5] * a[11] - a[12] * a[7] * a[9];
  inv[12] = -a[4] * a[9] * a[14] + a[4] * a[10] * a[13] + a[8] * a[5] * a[14] - a[8] * a[6] * a[13] - a[12] * a[5] * a[10] + a[12] * a[6] * a[9];
  inv[1] = -a[1] * a[10] * a[15] + a[1] * a[11] * a[14] + a[9] * a[2] * a[15] - a[9] * a[3] * a[14] - a[13] * a[2] * a[11] + a[13] * a[3] * a[10];
  inv[5] = a[0] * a[10] * a[15] - a[0] * a[11] * a[14] - a[8] * a[2] * a[15] + a[8] * a[3] * a[14] + a[12] * a[2] * a[11] - a[12] * a[3] * a[10];
  inv[9] = -a[0] * a[9] * a[15] + a[0] * a[11] * a[13] + a[8] * a[1] * a[15] - a[8] * a[3] * a[13] - a[12] * a[1] * a[11] + a[12] * a[3] * a[9];
  inv[13] = a[0] * a[9] * a[14] - a[0] * a[10] * a[13] - a[8] * a[1] * a[14] + a[8] * a[2] * a[13] + a[12] * a[1] * a[10] - a[12] * a[2] * a[9];
  inv[2] = a[1] * a[6] * a[15] - a[1] * a[7] * a[14] - a[5] * a[2] * a[15] + a[5] * a[3] * a[14] + a[13] * a[2] * a[7] - a[13] * a[3] * a[6];
  inv[6] = -a[0] * a[6] * a[15] + a[0] * a[7] * a[14] + a[4] * a[2] * a[15] - a[4] * a[3] * a[14] - a[12] * a[2] * a[7] + a[12] * a[3] * a[6];
  inv[10] = a[0] * a[5] * a[15] - a[0] * a[7] * a[13] - a[4] * a[1] * a[15] + a[4] * a[3] * a[13] + a[12] * a[1] * a[7] - a[12] * a[3] * a[5];
  inv[14] = -a[0] * a[5] * a[14] + a[0] * a[6] * a[13] + a[4] * a[1] * a[14] - a[4] * a[2] * a[13] - a[12] * a[1] * a[6] + a[12] * a[2] * a[5];
  inv[3] = -a[1] * a[6] * a[11] + a[1] * a[7] * a[10] + a[5] * a[2] * a[11] - a[5] * a[3] * a[10] - a[9] * a[2] * a[7] + a[9] * a[3] * a[6];
  inv[7] = a[0] * a[6] * a[11] - a[0] * a[7] * a[10] - a[4] * a[2] * a[11] + a[4] * a[3] * a[10] + a[8] * a[2] * a[7] - a[8] * a[3] * a[6];
  inv[11] = -a[0] * a[5] * a[11] + a[0] * a[7] * a[9] + a[4] * a[1] * a[11] - a[4] * a[3] * a[9] - a[8] * a[1] * a[7] + a[8] * a[3] * a[5];
  inv[15] = a[0] * a[5] * a[10] - a[0] * a[6] * a[9] - a[4] * a[1] * a[10] + a[4] * a[2] * a[9] + a[8] * a[1] * a[6] - a[8] * a[2] * a[5];
  let det = a[0] * inv[0] + a[1] * inv[4] + a[2] * inv[8] + a[3] * inv[12];
  det = det ? 1 / det : 0;
  for (let i = 0; i < 16; i++) inv[i] *= det;
  return inv;
}

/** Couleur CSS (#rgb, #rrggbb, rgb()) → [r, g, b] dans [0, 1]. */
export function cssRgb(s, fallback = [0.5, 0.5, 0.5]) {
  if (!s) return fallback;
  s = s.trim();
  let m = s.match(/^#([0-9a-f]{3})$/i);
  if (m) return [...m[1]].map((h) => parseInt(h + h, 16) / 255);
  m = s.match(/^#([0-9a-f]{6})/i);
  if (m) return [0, 2, 4].map((k) => parseInt(m[1].slice(k, k + 2), 16) / 255);
  m = s.match(/rgba?\(([^)]+)\)/);
  if (m) return m[1].split(/[ ,/]+/).slice(0, 3).map((v) => Number(v) / 255);
  return fallback;
}

export class Render3D {
  /**
   * @param {GPUDevice} device
   * @param {HTMLCanvasElement} canvas
   */
  constructor(device, canvas, o = {}) {
    this.device = device;
    this.canvas = canvas;
    // Rendu hors écran (tests) : image relue par readPixels().
    this.offscreen = !!o.offscreen;
    if (this.offscreen) {
      this.ctx = null;
      this.format = 'rgba8unorm';
    } else {
      this.ctx = canvas.getContext('webgpu');
      this.format = navigator.gpu.getPreferredCanvasFormat();
      this.ctx.configure({ device, format: this.format, alphaMode: 'opaque' });
    }
    this.view = device.createBuffer({ size: VIEW_BYTES, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
    this.sampler = device.createSampler({ magFilter: 'linear', minFilter: 'linear' });
    const mod = (code) => device.createShaderModule({ code });
    const depth = { format: 'depth32float', depthWriteEnabled: true, depthCompare: 'less' };
    const bodyMod = mod(BODY);
    this.bodyPipe = device.createRenderPipeline({
      layout: 'auto',
      vertex: {
        module: bodyMod,
        entryPoint: 'vs',
        buffers: [{ arrayStride: 24, attributes: [{ shaderLocation: 0, offset: 0, format: 'float32x3' }, { shaderLocation: 1, offset: 12, format: 'float32x3' }] }],
      },
      fragment: { module: bodyMod, entryPoint: 'fs', targets: [{ format: this.format }] },
      primitive: { topology: 'triangle-list', cullMode: 'none' },
      depthStencil: depth,
    });
    const lineMod = mod(LINES);
    this.linePipe = device.createRenderPipeline({
      layout: 'auto',
      vertex: { module: lineMod, entryPoint: 'vs', buffers: [{ arrayStride: 12, attributes: [{ shaderLocation: 0, offset: 0, format: 'float32x3' }] }] },
      fragment: { module: lineMod, entryPoint: 'fs', targets: [{ format: this.format }] },
      primitive: { topology: 'line-list' },
      depthStencil: depth,
    });
    const sliceMod = mod(SLICE);
    this.slicePipe = device.createRenderPipeline({
      layout: 'auto',
      vertex: { module: sliceMod, entryPoint: 'vs' },
      fragment: { module: sliceMod, entryPoint: 'fs', targets: [{ format: this.format }] },
      primitive: { topology: 'triangle-list', cullMode: 'none' },
      depthStencil: depth,
    });
    const isoMod = mod(ISO);
    this.isoPipe = device.createRenderPipeline({
      layout: 'auto',
      vertex: { module: isoMod, entryPoint: 'vs' },
      fragment: {
        module: isoMod,
        entryPoint: 'fs',
        targets: [
          {
            format: this.format,
            blend: {
              color: { srcFactor: 'src-alpha', dstFactor: 'one-minus-src-alpha', operation: 'add' },
              alpha: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha', operation: 'add' },
            },
          },
        ],
      },
      primitive: { topology: 'triangle-list' },
    });
    this.derivePipe = device.createComputePipeline({ layout: 'auto', compute: { module: mod(DERIVE), entryPoint: 'main' } });
    this.depth = null;
  }

  /** Associe un calcul et sa maquette. */
  attach(sim, body) {
    const d = this.device;
    this.sim = sim;
    this.body = body;
    this.derivTex?.destroy();
    this.derivTex = d.createTexture({
      size: [sim.nx, sim.ny, sim.nz],
      dimension: '3d',
      format: 'rgba16float',
      usage: GPUTextureUsage.STORAGE_BINDING | GPUTextureUsage.TEXTURE_BINDING,
    });
    this.mesh?.destroy();
    this.mesh = d.createBuffer({ size: Math.max(24, body.mesh.byteLength), usage: GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST });
    d.queue.writeBuffer(this.mesh, 0, body.mesh);
    this.meshCount = body.mesh.length / 6;
    const X = sim.nx - 1;
    const Y = sim.ny - 1;
    const Z = sim.nz - 1;
    const c = [
      [0, 0, 0], [X, 0, 0], [X, Y, 0], [0, Y, 0],
      [0, 0, Z], [X, 0, Z], [X, Y, Z], [0, Y, Z],
    ];
    const e = [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4], [0, 4], [1, 5], [2, 6], [3, 7]];
    const lines = new Float32Array(e.flatMap(([a, b]) => [...c[a], ...c[b]]));
    this.lines?.destroy();
    this.lines = d.createBuffer({ size: lines.byteLength, usage: GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST });
    d.queue.writeBuffer(this.lines, 0, lines);
    this.deriveBG = d.createBindGroup({
      layout: this.derivePipe.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: sim.macroTex.createView() },
        { binding: 1, resource: this.derivTex.createView() },
      ],
    });
    this.sliceBG = d.createBindGroup({
      layout: this.slicePipe.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: this.view } },
        { binding: 1, resource: sim.macroTex.createView() },
        { binding: 2, resource: this.derivTex.createView() },
        { binding: 3, resource: this.sampler },
      ],
    });
    this.bodyBG = d.createBindGroup({ layout: this.bodyPipe.getBindGroupLayout(0), entries: [{ binding: 0, resource: { buffer: this.view } }] });
    this.lineBG = d.createBindGroup({ layout: this.linePipe.getBindGroupLayout(0), entries: [{ binding: 0, resource: { buffer: this.view } }] });
    this.isoBG = null;
  }

  resize() {
    const dpr = Math.min(window.devicePixelRatio || 1, 1.5);
    const w = Math.max(1, Math.round(this.canvas.clientWidth * dpr));
    const h = Math.max(1, Math.round(this.canvas.clientHeight * dpr));
    if (this.canvas.width !== w || this.canvas.height !== h || !this.depth) {
      this.canvas.width = w;
      this.canvas.height = h;
      this.depth?.destroy();
      this.depth = this.device.createTexture({ size: [w, h], format: 'depth32float', usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING });
      if (this.offscreen) {
        this.colorTex?.destroy();
        this.colorTex = this.device.createTexture({ size: [w, h], format: this.format, usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_SRC });
      }
      this.isoBG = null;
    }
  }

  /**
   * @param {object} o { az, el, dist (en L), plane: {axis, pos (mailles)}, field: 'speed'|'vort',
   *        iso: {on, level}, theme: {bg, body, line}, U }
   */
  render(o) {
    if (!this.sim) return;
    const d = this.device;
    this.resize();
    const { sim, body } = this;
    const w = this.canvas.width;
    const h = this.canvas.height;
    const L = body.L;
    const target = body.center;
    const r = o.dist * L;
    const az = (o.az * Math.PI) / 180;
    const el = (o.el * Math.PI) / 180;
    const eye = [target[0] - r * Math.cos(el) * Math.cos(az), target[1] - r * Math.cos(el) * Math.sin(az), target[2] + r * Math.sin(el)];
    const proj = perspective((40 * Math.PI) / 180, w / h, Math.max(0.5, r * 0.02), r * 4 + 4 * Math.max(sim.nx, sim.ny, sim.nz));
    const viewM = lookAt(eye, target, [0, 0, 1]);
    const vp = mul4(proj, viewM);
    const ivp = invert4(vp);
    const b = new Float32Array(VIEW_BYTES / 4);
    b.set(vp, 0);
    b.set(ivp, 16);
    b.set([...eye, 1], 32);
    b.set([sim.nx, sim.ny, sim.nz, L], 36);
    b.set([...o.theme.bg, 1], 40);
    b.set([...o.theme.body, 1], 44);
    b.set([...o.theme.line, 1], 48);
    const axis = { x: 0, y: 1, z: 2 }[o.plane.axis];
    b.set([axis, o.plane.pos, o.field === 'vort' ? 1 : 0, o.U], 52);
    b.set([o.iso.level * (o.U / L) ** 2, o.iso.on ? 1 : 0, o.iso.step ?? 0.75, 4], 56);
    b.set([w, h, 0, 0], 60);
    d.queue.writeBuffer(this.view, 0, b);
    const enc = d.createCommandEncoder();
    const cp = enc.beginComputePass();
    cp.setPipeline(this.derivePipe);
    cp.setBindGroup(0, this.deriveBG);
    cp.dispatchWorkgroups(Math.ceil(sim.nx / 4), Math.ceil(sim.ny / 4), Math.ceil(sim.nz / 4));
    cp.end();
    const colorView = this.offscreen ? this.colorTex.createView() : this.ctx.getCurrentTexture().createView();
    const [br, bg, bb] = o.theme.bg;
    const pass = enc.beginRenderPass({
      colorAttachments: [{ view: colorView, clearValue: { r: br, g: bg, b: bb, a: 1 }, loadOp: 'clear', storeOp: 'store' }],
      depthStencilAttachment: { view: this.depth.createView(), depthClearValue: 1, depthLoadOp: 'clear', depthStoreOp: 'store' },
    });
    pass.setPipeline(this.bodyPipe);
    pass.setBindGroup(0, this.bodyBG);
    pass.setVertexBuffer(0, this.mesh);
    pass.draw(this.meshCount);
    pass.setPipeline(this.slicePipe);
    pass.setBindGroup(0, this.sliceBG);
    pass.draw(6);
    pass.setPipeline(this.linePipe);
    pass.setBindGroup(0, this.lineBG);
    pass.setVertexBuffer(0, this.lines);
    pass.draw(24);
    pass.end();
    if (o.iso.on) {
      if (!this.isoBG) {
        this.isoBG = d.createBindGroup({
          layout: this.isoPipe.getBindGroupLayout(0),
          entries: [
            { binding: 0, resource: { buffer: this.view } },
            { binding: 1, resource: this.derivTex.createView() },
            { binding: 2, resource: this.sampler },
            { binding: 3, resource: this.depth.createView() },
            { binding: 4, resource: sim.macroTex.createView() },
          ],
        });
      }
      const p2 = enc.beginRenderPass({ colorAttachments: [{ view: colorView, loadOp: 'load', storeOp: 'store' }] });
      p2.setPipeline(this.isoPipe);
      p2.setBindGroup(0, this.isoBG);
      p2.draw(3);
      p2.end();
    }
    d.queue.submit([enc.finish()]);
  }

  /** Image du rendu hors écran : { w, h, data (RGBA 8 bits) }. */
  async readPixels() {
    const d = this.device;
    const w = this.canvas.width;
    const h = this.canvas.height;
    const bpr = Math.ceil((w * 4) / 256) * 256;
    const buf = d.createBuffer({ size: bpr * h, usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ });
    const enc = d.createCommandEncoder();
    enc.copyTextureToBuffer({ texture: this.colorTex }, { buffer: buf, bytesPerRow: bpr }, [w, h]);
    d.queue.submit([enc.finish()]);
    await buf.mapAsync(GPUMapMode.READ);
    const src = new Uint8Array(buf.getMappedRange());
    const data = new Uint8ClampedArray(w * h * 4);
    for (let y = 0; y < h; y++) data.set(src.subarray(y * bpr, y * bpr + w * 4), y * w * 4);
    buf.unmap();
    buf.destroy();
    return { w, h, data };
  }

  dispose() {
    this.colorTex?.destroy();
    this.derivTex?.destroy();
    this.mesh?.destroy();
    this.lines?.destroy();
    this.depth?.destroy();
    this.view.destroy();
  }
}
