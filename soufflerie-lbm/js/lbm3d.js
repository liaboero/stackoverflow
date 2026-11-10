/**
 * Soufflerie 3D : Lattice Boltzmann D3Q19 sur WebGPU (calcul).
 *
 * Un pas = une passe de calcul : propagation « pull » et collision fusionnées.
 * Populations stockées décalées de leur poids, g = f − w, en demi-précision
 * (deux par mot de 32 bits, pack2x16float), ou en simple précision au choix.
 * Collision régularisée (projection d'ordre 2) avec Smagorinsky aux Reynolds
 * de maille élevés. Parois : rebond interpolé de Bouzidi (fractions q des
 * liens calculées sur le CPU), efforts par échange de quantité de mouvement
 * accumulés par maille pariétale. Entrée à vitesse imposée (équilibre),
 * sortie à gradient nul avec zone éponge, bords latéraux périodiques.
 */

export const C3 = [
  [0, 0, 0],
  [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1],
  [1, 1, 0], [-1, -1, 0], [1, 0, 1], [-1, 0, -1], [0, 1, 1], [0, -1, -1],
  [1, -1, 0], [-1, 1, 0], [1, 0, -1], [-1, 0, 1], [0, 1, -1], [0, -1, 1],
];
export const W3 = C3.map((c, i) => (i === 0 ? 1 / 3 : Math.abs(c[0]) + Math.abs(c[1]) + Math.abs(c[2]) === 1 ? 1 / 18 : 1 / 36));
const OPP3 = C3.map((c) => C3.findIndex((d) => d[0] === -c[0] && d[1] === -c[1] && d[2] === -c[2]));
const SOLID = 0xffffffff;

/** WebGPU disponible ? Renvoie l'adaptateur ou null. */
export async function webgpuAdapter() {
  if (typeof navigator === 'undefined' || !navigator.gpu) return null;
  try {
    return await navigator.gpu.requestAdapter({ powerPreference: 'high-performance' });
  } catch {
    return null;
  }
}

/** Octets par maille (populations ping-pong + table pariétale). */
export function bytesPerCell(fp16) {
  return (fp16 ? 2 * 10 * 4 : 2 * 19 * 4) + 4;
}

function wgslStep(fp16) {
  const words = fp16 ? 10 : 19;
  const half = fp16 ? 5 : 10; // mots dans le premier tampon
  const f = (v) => (Number.isInteger(v) ? `${v}.0` : `${v}`);
  // Lecture d'une population i en maille c.
  const load = fp16
    ? `fn pop(i: u32, c: u32) -> f32 {
  let w = i >> 1u;
  var v: vec2<f32>;
  if (w < ${half}u) { v = unpack2x16float(srcA[w * P.N + c]); } else { v = unpack2x16float(srcB[(w - ${half}u) * P.N + c]); }
  return select(v.x, v.y, (i & 1u) == 1u);
}`
    : `fn pop(i: u32, c: u32) -> f32 {
  if (i < ${half}u) { return bitcast<f32>(srcA[i * P.N + c]); }
  return bitcast<f32>(srcB[(i - ${half}u) * P.N + c]);
}`;
  const store = fp16
    ? `fn putWord(w: u32, c: u32, a: f32, b: f32) {
  let v = pack2x16float(vec2<f32>(a, b));
  if (w < ${half}u) { dstA[w * P.N + c] = v; } else { dstB[(w - ${half}u) * P.N + c] = v; }
}`
    : `fn putPop(i: u32, c: u32, a: f32) {
  if (i < ${half}u) { dstA[i * P.N + c] = bitcast<u32>(a); } else { dstB[(i - ${half}u) * P.N + c] = bitcast<u32>(a); }
}`;
  const C = C3;
  const lines = [];
  // Propagation : pour chaque direction i, population venant de x − c_i.
  for (let i = 1; i < 19; i++) {
    const [cx, cy, cz] = C[i];
    const j = OPP3[i];
    lines.push(`  {
    // direction ${i} (${cx}, ${cy}, ${cz})
    var v: f32;
    var q = 0.0;
    if (hasWall) { q = qbuf[wk * 18u + ${i - 1}u]; }
    if (q > 0.0) {
      // Lien vers la paroi : rebond interpolé (Bouzidi), g = f − w affine.
      let fj = pop(${j}u, idx);
      let qo = qbuf[wk * 18u + ${j - 1}u];
      if (q < 0.5 && qo == 0.0) {
        let nb = cell(x + ${cx}, y + ${cy}, z + ${cz});
        v = 2.0 * q * fj + (1.0 - 2.0 * q) * pop(${j}u, nb);
      } else if (q < 0.5) {
        v = fj;
      } else {
        v = (fj + (2.0 * q - 1.0) * pop(${i}u, idx)) / (2.0 * q);
      }
      let mom = fj + v + ${f(2 * W3[i])};
      fsum += vec3<f32>(${f(-cx)}, ${f(-cy)}, ${f(-cz)}) * mom;
    } else if (x - ${cx} < 0) {
      v = feqIn(${i}u);
    } else if (x - ${cx} >= i32(P.nx)) {
      v = pop(${i}u, idx);
    } else {
      v = pop(${i}u, cell(x - ${cx}, y - ${cy}, z - ${cz}));
    }
    g[${i}] = v;
  }`);
  }
  const pull = lines.join('\n');
  const dot = (i, a, b) => C[i][a] * C[i][b];
  const sumExpr = (fn) => {
    const t = [];
    for (let i = 0; i < 19; i++) {
      const k = fn(i);
      if (k === 0) continue;
      t.push(`${k === 1 ? '' : k === -1 ? '-' : `${k} * `}fv[${i}]`.replace('- ', '-'));
    }
    return t.join(' + ').replace(/\+ -/g, '- ') || '0.0';
  };
  const coll = [];
  for (let i = 0; i < 19; i++) {
    const [cx, cy, cz] = C[i];
    const cu = [cx && `${cx === 1 ? '' : '-'}u.x`, cy && `${cy === 1 ? '' : '-'}u.y`, cz && `${cz === 1 ? '' : '-'}u.z`].filter(Boolean).join(' + ').replace(/\+ -/g, '- ') || '0.0';
    const qp = `${f(cx * cx - 1 / 3)} * pxx + ${f(cy * cy - 1 / 3)} * pyy + ${f(cz * cz - 1 / 3)} * pzz + ${f(2 * cx * cy)} * pxy + ${f(2 * cx * cz)} * pxz + ${f(2 * cy * cz)} * pyz`;
    coll.push(`  {
    let cu = ${cu};
    let feq = ${f(W3[i])} * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - usq);
    let fnew = feq + k * ${f(4.5 * W3[i])} * (${qp}) + ${f(3 * W3[i])} * rho * (${f(cx)} * P.accel);
    o[${i}] = fnew - ${f(W3[i])};
  }`);
  }
  const out = fp16
    ? Array.from({ length: words }, (_, w) => `  putWord(${w}u, idx, o[${2 * w}], ${2 * w + 1 < 19 ? `o[${2 * w + 1}]` : '0.0'});`).join('\n')
    : Array.from({ length: 19 }, (_, i) => `  putPop(${i}u, idx, o[${i}]);`).join('\n');
  return /* wgsl */ `
struct Params {
  nx: u32, ny: u32, nz: u32, N: u32,
  tau0: f32, tauSponge: f32, spongeStart: u32, uIn: f32,
  accel: f32, c18: f32, writeMacro: u32, les: u32,
  rows: u32, pad0: u32, pad1: u32, pad2: u32,
};
@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read> srcA: array<u32>;
@group(0) @binding(2) var<storage, read> srcB: array<u32>;
@group(0) @binding(3) var<storage, read_write> dstA: array<u32>;
@group(0) @binding(4) var<storage, read_write> dstB: array<u32>;
@group(0) @binding(5) var<storage, read> wall: array<u32>;
@group(0) @binding(6) var<storage, read> qbuf: array<f32>;
@group(0) @binding(7) var<storage, read_write> force: array<vec4<f32>>;
@group(0) @binding(8) var macroTex: texture_storage_3d<rgba16float, write>;

${load}
${store}

fn cell(xi: i32, yi: i32, zi: i32) -> u32 {
  let ny = i32(P.ny);
  let nz = i32(P.nz);
  let yy = (yi + ny) % ny;
  let zz = (zi + nz) % nz;
  let xx = clamp(xi, 0, i32(P.nx) - 1);
  return u32(xx) + P.nx * (u32(yy) + P.ny * u32(zz));
}

fn feqIn(i: u32) -> f32 {
  // Équilibre amont (ρ = 1, u = (U, 0, 0)), décalé de w_i.
  var cx = 0.0;
  var w = 1.0 / 36.0;
  ${C.map((c, i) => (i === 0 ? '' : `if (i == ${i}u) { cx = ${f(c[0])}; w = ${f(W3[i])}; }`)).filter(Boolean).join('\n  ')}
  let cu = cx * P.uIn;
  return w * (3.0 * cu + 4.5 * cu * cu - 1.5 * P.uIn * P.uIn);
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let idx = gid.x + gid.y * nwg.x * 64u;
  if (idx >= P.N) { return; }
  let wk0 = wall[idx];
  if (wk0 == ${SOLID}u) { return; }
  // wk : indice de maille pariétale, ou 0xffffffff si aucun lien de paroi.
  let wk = select(wk0 - 1u, 0xffffffffu, wk0 == 0u);
  let hasWall = wk != 0xffffffffu;
  let x = i32(idx % P.nx);
  let y = i32((idx / P.nx) % P.ny);
  let z = i32(idx / (P.nx * P.ny));
  var g: array<f32, 19>;
  var fsum = vec3<f32>(0.0);
  g[0] = pop(0u, idx);
${pull}
  var fv: array<f32, 19>;
  ${Array.from({ length: 19 }, (_, i) => `fv[${i}] = g[${i}] + ${f(W3[i])};`).join(' ')}
  let rho = ${sumExpr(() => 1)};
  let ir = 1.0 / rho;
  let u = vec3<f32>(${sumExpr((i) => C[i][0])}, ${sumExpr((i) => C[i][1])}, ${sumExpr((i) => C[i][2])}) * ir;
  let usq = 1.5 * dot(u, u);
  // Tenseur hors équilibre Π − Π_eq.
  let pxx = ${sumExpr((i) => dot(i, 0, 0))} - rho * (1.0 / 3.0 + u.x * u.x);
  let pyy = ${sumExpr((i) => dot(i, 1, 1))} - rho * (1.0 / 3.0 + u.y * u.y);
  let pzz = ${sumExpr((i) => dot(i, 2, 2))} - rho * (1.0 / 3.0 + u.z * u.z);
  let pxy = ${sumExpr((i) => dot(i, 0, 1))} - rho * u.x * u.y;
  let pxz = ${sumExpr((i) => dot(i, 0, 2))} - rho * u.x * u.z;
  let pyz = ${sumExpr((i) => dot(i, 1, 2))} - rho * u.y * u.z;
  var tau = P.tau0;
  if (u32(x) > P.spongeStart) {
    let s = f32(u32(x) - P.spongeStart) / f32(max(1u, P.nx - 1u - P.spongeStart));
    tau += (P.tauSponge - P.tau0) * s * s;
  }
  var om = 1.0 / tau;
  if (P.les == 1u) {
    let qn = sqrt(2.0 * (pxx * pxx + pyy * pyy + pzz * pzz + 2.0 * (pxy * pxy + pxz * pxz + pyz * pyz)));
    om = 2.0 / (tau + sqrt(tau * tau + P.c18 * qn * ir));
  }
  let k = 1.0 - om;
  var o: array<f32, 19>;
${coll.join('\n')}
${out}
  if (hasWall) {
    force[wk] += vec4<f32>(fsum, 0.0);
  }
  if (P.writeMacro == 1u) {
    textureStore(macroTex, vec3<u32>(u32(x), u32(y), u32(z)), vec4<f32>(u, rho));
  }
}
`;
}

function wgslInit(fp16) {
  const words = fp16 ? 10 : 19;
  const half = fp16 ? 5 : 10;
  const f = (v) => (Number.isInteger(v) ? `${v}.0` : `${v}`);
  const eq = W3.map((w, i) => `let e${i} = ${f(w)} * (3.0 * (${f(C3[i][0])} * u0) + 4.5 * (${f(C3[i][0])} * u0) * (${f(C3[i][0])} * u0) - 1.5 * u0 * u0);`).join('\n  ');
  const put = fp16
    ? Array.from({ length: words }, (_, w) => {
        const v = `pack2x16float(vec2<f32>(e${2 * w}, ${2 * w + 1 < 19 ? `e${2 * w + 1}` : '0.0'}))`;
        return w < half ? `A[${w}u * P.N + idx] = ${v};` : `B[${w - half}u * P.N + idx] = ${v};`;
      }).join('\n  ')
    : Array.from({ length: 19 }, (_, i) => (i < half ? `A[${i}u * P.N + idx] = bitcast<u32>(e${i});` : `B[${i - half}u * P.N + idx] = bitcast<u32>(e${i});`)).join('\n  ');
  return /* wgsl */ `
struct Params {
  nx: u32, ny: u32, nz: u32, N: u32,
  tau0: f32, tauSponge: f32, spongeStart: u32, uIn: f32,
  accel: f32, c18: f32, writeMacro: u32, les: u32,
  rows: u32, pad0: u32, pad1: u32, pad2: u32,
};
@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read_write> A: array<u32>;
@group(0) @binding(2) var<storage, read_write> B: array<u32>;
@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let idx = gid.x + gid.y * nwg.x * 64u;
  if (idx >= P.N) { return; }
  let u0 = P.uIn;
  ${eq}
  ${put}
}
`;
}

/** Fraction q ∈ ]0, 1] du lien p → p + d côté fluide (dichotomie). */
function linkQ(inside, x, y, z, dx, dy, dz) {
  let lo = 0;
  let hi = 1;
  for (let k = 0; k < 12; k++) {
    const m = 0.5 * (lo + hi);
    if (inside(x + m * dx, y + m * dy, z + m * dz)) hi = m;
    else lo = m;
  }
  return Math.min(1, Math.max(0.02, 0.5 * (lo + hi)));
}

export class Lbm3D {
  /**
   * @param {GPUDevice} device
   * @param {object} o { nx, ny, nz, nu, uIn, fp16 }
   */
  constructor(device, o) {
    this.device = device;
    this.nx = o.nx;
    this.ny = o.ny;
    this.nz = o.nz;
    this.N = o.nx * o.ny * o.nz;
    this.fp16 = o.fp16 !== false;
    this.uIn = o.uIn ?? 0.1;
    this.uNow = 0;
    this.accel = 0;
    this.time = 0;
    this.cs = 0.15;
    this.spongeStart = Math.round(this.nx * 0.85);
    this.setViscosity(o.nu);
    const N = this.N;
    const wordsA = this.fp16 ? 5 : 10;
    const wordsB = this.fp16 ? 5 : 9;
    const S = GPUBufferUsage.STORAGE;
    const mk = (words) => device.createBuffer({ size: Math.max(16, words * N * 4), usage: S | GPUBufferUsage.COPY_DST | GPUBufferUsage.COPY_SRC });
    this.pops = [
      [mk(wordsA), mk(wordsB)],
      [mk(wordsA), mk(wordsB)],
    ];
    this.wall = device.createBuffer({ size: N * 4, usage: S | GPUBufferUsage.COPY_DST });
    this.params = device.createBuffer({ size: 64, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
    this.macroTex = device.createTexture({
      size: [this.nx, this.ny, this.nz],
      dimension: '3d',
      format: 'rgba16float',
      usage: GPUTextureUsage.STORAGE_BINDING | GPUTextureUsage.TEXTURE_BINDING,
    });
    this.stepModule = device.createShaderModule({ code: wgslStep(this.fp16) });
    this.stepPipe = device.createComputePipeline({ layout: 'auto', compute: { module: this.stepModule, entryPoint: 'main' } });
    this.initPipe = device.createComputePipeline({
      layout: 'auto',
      compute: { module: device.createShaderModule({ code: wgslInit(this.fp16) }), entryPoint: 'main' },
    });
    this.cur = 0;
    this.wallCount = 0;
    this.qbuf = null;
    this.force = null;
    this.forceSteps = 0;
    const groups = Math.ceil(N / 64);
    this.dispX = Math.min(groups, 65535);
    this.dispY = Math.ceil(groups / this.dispX);
    this.pendingRead = null;
  }

  setViscosity(nu) {
    this.nu = nu;
    this.tau0 = 0.5 + 3 * nu;
    this.tauSponge = Math.max(this.tau0, 0.8);
    // Sous-maille seulement aux Reynolds de maille élevés : sans amortissement
    // pariétal, elle épaissirait la couche limite et gonflerait la traînée.
    this.les = this.tau0 < 0.51;
  }

  /**
   * Voxelise le corps : mailles solides, fractions q des liens pariétaux.
   * @param {(x:number, y:number, z:number) => boolean} inside (coordonnées réseau)
   * @param {{x0,x1,y0,y1,z0,z1}} box boîte englobante (réseau)
   */
  setGeometry(inside, box) {
    const { nx, ny, nz, N, device } = this;
    const wall = new Uint32Array(N);
    const X0 = Math.max(1, Math.floor(box.x0) - 1);
    const X1 = Math.min(nx - 2, Math.ceil(box.x1) + 1);
    const Y0 = Math.max(0, Math.floor(box.y0) - 1);
    const Y1 = Math.min(ny - 1, Math.ceil(box.y1) + 1);
    const Z0 = Math.max(0, Math.floor(box.z0) - 1);
    const Z1 = Math.min(nz - 1, Math.ceil(box.z1) + 1);
    const id = (x, y, z) => x + nx * (y + ny * z);
    let solids = 0;
    for (let z = Z0; z <= Z1; z++)
      for (let y = Y0; y <= Y1; y++)
        for (let x = X0; x <= X1; x++)
          if (inside(x, y, z)) {
            wall[id(x, y, z)] = SOLID;
            solids++;
          }
    const qs = [];
    let k = 0;
    for (let z = Math.max(0, Z0 - 1); z <= Math.min(nz - 1, Z1 + 1); z++)
      for (let y = Math.max(0, Y0 - 1); y <= Math.min(ny - 1, Y1 + 1); y++)
        for (let x = Math.max(1, X0 - 1); x <= Math.min(nx - 2, X1 + 1); x++) {
          const c = id(x, y, z);
          if (wall[c] === SOLID) continue;
          let any = false;
          const row = new Float32Array(18);
          for (let i = 1; i < 19; i++) {
            const [cx, cy, cz] = C3[i];
            const sx = x - cx;
            const sy = (y - cy + ny) % ny;
            const sz = (z - cz + nz) % nz;
            if (wall[id(sx, sy, sz)] === SOLID) {
              row[i - 1] = linkQ(inside, x, y, z, -cx, -cy, -cz);
              any = true;
            }
          }
          if (any) {
            wall[c] = ++k;
            qs.push(row);
          }
        }
    this.wallCount = k;
    this.solidCount = solids;
    const q = new Float32Array(Math.max(18, 18 * k));
    qs.forEach((row, i) => q.set(row, 18 * i));
    this.qbuf?.destroy();
    this.force?.destroy();
    this.qbuf = device.createBuffer({ size: q.byteLength, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST });
    device.queue.writeBuffer(this.qbuf, 0, q);
    const fsize = Math.max(16, 16 * k);
    this.force = device.createBuffer({ size: fsize, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC | GPUBufferUsage.COPY_DST });
    this.forceRead = [0, 1].map(() => device.createBuffer({ size: fsize, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST }));
    this.readBusy = [false, false];
    this.readSlot = 0;
    device.queue.writeBuffer(this.wall, 0, wall);
    this.wallHost = wall;
    this.bindGroups = [0, 1].map((s) =>
      device.createBindGroup({
        layout: this.stepPipe.getBindGroupLayout(0),
        entries: [
          { binding: 0, resource: { buffer: this.params } },
          { binding: 1, resource: { buffer: this.pops[s][0] } },
          { binding: 2, resource: { buffer: this.pops[s][1] } },
          { binding: 3, resource: { buffer: this.pops[1 - s][0] } },
          { binding: 4, resource: { buffer: this.pops[1 - s][1] } },
          { binding: 5, resource: { buffer: this.wall } },
          { binding: 6, resource: { buffer: this.qbuf } },
          { binding: 7, resource: { buffer: this.force } },
          { binding: 8, resource: this.macroTex.createView() },
        ],
      }),
    );
  }

  /** Fluide au repos (ou à la vitesse u0). */
  reset(u0 = 0) {
    const d = this.device;
    this.uNow = u0;
    this.writeParams(0);
    for (const s of [0, 1]) {
      const bg = d.createBindGroup({
        layout: this.initPipe.getBindGroupLayout(0),
        entries: [
          { binding: 0, resource: { buffer: this.params } },
          { binding: 1, resource: { buffer: this.pops[s][0] } },
          { binding: 2, resource: { buffer: this.pops[s][1] } },
        ],
      });
      const enc = d.createCommandEncoder();
      const pass = enc.beginComputePass();
      pass.setPipeline(this.initPipe);
      pass.setBindGroup(0, bg);
      pass.dispatchWorkgroups(this.dispX, this.dispY);
      pass.end();
      d.queue.submit([enc.finish()]);
    }
    this.cur = 0;
    this.time = 0;
    this.clearForces();
  }

  clearForces() {
    const enc = this.device.createCommandEncoder();
    enc.clearBuffer(this.force);
    this.device.queue.submit([enc.finish()]);
    this.forceSteps = 0;
  }

  writeParams(writeMacro) {
    const b = new ArrayBuffer(64);
    const u = new Uint32Array(b);
    const f = new Float32Array(b);
    u[0] = this.nx;
    u[1] = this.ny;
    u[2] = this.nz;
    u[3] = this.N;
    f[4] = this.tau0;
    f[5] = this.tauSponge;
    u[6] = this.spongeStart;
    f[7] = this.uNow;
    f[8] = this.accel;
    f[9] = 18 * this.cs * this.cs;
    u[10] = writeMacro ? 1 : 0;
    u[11] = this.les ? 1 : 0;
    this.device.queue.writeBuffer(this.params, 0, b);
  }

  /**
   * n pas ; l'entrée suit la loi de vitesse uAt(time) (mise en vitesse).
   * Le champ macroscopique (rendu) est écrit au dernier pas.
   */
  step(n, uAt) {
    const d = this.device;
    for (let k = 0; k < n; k++) {
      const u = uAt ? uAt(this.time) : this.uIn;
      const u1 = uAt ? uAt(this.time + 1) : this.uIn;
      this.uNow = u;
      this.accel = u1 - u;
      this.writeParams(k === n - 1);
      const enc = d.createCommandEncoder();
      const pass = enc.beginComputePass();
      pass.setPipeline(this.stepPipe);
      pass.setBindGroup(0, this.bindGroups[this.cur]);
      pass.dispatchWorkgroups(this.dispX, this.dispY);
      pass.end();
      d.queue.submit([enc.finish()]);
      this.cur = 1 - this.cur;
      this.time++;
      this.forceSteps++;
    }
  }

  /**
   * Relit (asynchrone) la force moyenne depuis la dernière relecture, puis
   * remet l'accumulateur à zéro. Renvoie null si une relecture est en cours.
   */
  async readForce() {
    if (this.pendingRead || !this.wallCount || !this.forceSteps) return null;
    const d = this.device;
    const buf = this.forceRead[this.readSlot];
    this.readSlot = 1 - this.readSlot;
    const steps = this.forceSteps;
    const time = this.time;
    const enc = d.createCommandEncoder();
    enc.copyBufferToBuffer(this.force, 0, buf, 0, 16 * this.wallCount);
    enc.clearBuffer(this.force);
    d.queue.submit([enc.finish()]);
    this.forceSteps = 0;
    this.pendingRead = buf.mapAsync(GPUMapMode.READ);
    try {
      await this.pendingRead;
    } catch {
      this.pendingRead = null;
      return null;
    }
    const a = new Float32Array(buf.getMappedRange().slice(0));
    buf.unmap();
    this.pendingRead = null;
    let fx = 0;
    let fy = 0;
    let fz = 0;
    for (let k = 0; k < this.wallCount; k++) {
      fx += a[4 * k];
      fy += a[4 * k + 1];
      fz += a[4 * k + 2];
    }
    return { fx: fx / steps, fy: fy / steps, fz: fz / steps, steps, time };
  }

  dispose() {
    for (const p of this.pops) for (const b of p) b.destroy();
    this.wall.destroy();
    this.params.destroy();
    this.qbuf?.destroy();
    this.force?.destroy();
    this.forceRead?.forEach((b) => b.destroy());
    this.macroTex.destroy();
  }
}
