/**
 * Maquettes de la soufflerie 3D : domaine, test d'appartenance (voxelisation),
 * surface de référence et maillage d'affichage. Coordonnées réseau : x vers
 * l'aval, y latéral (envergure), z vers le haut.
 */

const DEG = Math.PI / 180;

/** Catalogue. domain : dimensions de la veine en longueurs de référence L. */
export const BODIES3D = {
  sphere: { label: 'Sphère', ref: 'diamètre D', domain: [12, 6, 6], at: [4, 3, 3], alpha: false },
  cube: { label: 'Cube', ref: 'arête a', domain: [12, 6, 6], at: [4, 3, 3], alpha: true },
  cylinder: { label: 'Cylindre', ref: 'diamètre D (envergure périodique 4D)', domain: [16, 4, 8], at: [5, 2, 4], alpha: false },
  wing: { label: 'Aile', ref: 'corde c (NACA 0012, allongement 4)', domain: [9, 8, 4], at: [2.5, 4, 2], alpha: true },
};

function naca00(t, x) {
  return 5 * t * (0.2969 * Math.sqrt(x) - 0.126 * x - 0.3516 * x * x + 0.2843 * x ** 3 - 0.1036 * x ** 4);
}

/**
 * @param {string} id
 * @param {number} L longueur de référence (mailles)
 * @param {object} p { alpha (°), ar (allongement de l'aile) }
 */
export function makeBody3D(id, L, p = {}) {
  const b = BODIES3D[id];
  const [dx, dy, dz] = b.domain;
  const nx = Math.round(dx * L);
  const ny = Math.round(dy * L);
  const nz = Math.round(dz * L);
  // Centre décalé d'un quart de maille : évite les symétries d'arrondi exactes.
  const cx = b.at[0] * L + 0.25;
  const cy = b.at[1] * L + 0.25;
  const cz = b.at[2] * L + 0.25;
  const a = (b.alpha ? p.alpha ?? 0 : 0) * DEG;
  const ca = Math.cos(a);
  const sa = Math.sin(a);
  const R = L / 2;
  let inside;
  let box;
  let area;
  let mesh;
  if (id === 'sphere') {
    inside = (x, y, z) => (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2 <= R * R;
    box = { x0: cx - R, x1: cx + R, y0: cy - R, y1: cy + R, z0: cz - R, z1: cz + R };
    area = (Math.PI * L * L) / 4;
    mesh = sphereMesh(cx, cy, cz, R);
  } else if (id === 'cube') {
    inside = (x, y, z) => {
      const X = x - cx;
      const Z = z - cz;
      const xb = X * ca - Z * sa;
      const zb = X * sa + Z * ca;
      return Math.abs(xb) <= R && Math.abs(y - cy) <= R && Math.abs(zb) <= R;
    };
    const e = R * Math.SQRT2;
    box = { x0: cx - e, x1: cx + e, y0: cy - R, y1: cy + R, z0: cz - e, z1: cz + e };
    area = L * L;
    mesh = boxMesh(cx, cy, cz, R, ca, sa);
  } else if (id === 'cylinder') {
    inside = (x, y, z) => (x - cx) ** 2 + (z - cz) ** 2 <= R * R;
    box = { x0: cx - R, x1: cx + R, y0: 0, y1: ny - 1, z0: cz - R, z1: cz + R };
    area = L * ny;
    mesh = cylinderMesh(cx, cz, R, ny);
  } else {
    const ar = p.ar ?? 4;
    const half = (ar * L) / 2;
    const t = 0.12;
    // Origine au quart de corde ; cabrage : bord de fuite vers le bas.
    inside = (x, y, z) => {
      if (Math.abs(y - cy) > half) return false;
      const X = x - cx;
      const Z = z - cz;
      const xb = X * ca - Z * sa;
      const zb = X * sa + Z * ca;
      const xc = xb / L + 0.25;
      if (xc < 0 || xc > 1) return false;
      return Math.abs(zb / L) <= naca00(t, xc);
    };
    const r = 0.8 * L;
    box = { x0: cx - r, x1: cx + r, y0: cy - half, y1: cy + half, z0: cz - r, z1: cz + r };
    area = L * 2 * half;
    mesh = wingMesh(cx, cy, cz, L, half, t, ca, sa);
  }
  return { id, L, nx, ny, nz, inside, box, area, mesh, center: [cx, cy, cz], label: b.label, ref: b.ref };
}

// ------------------------------------------------------------ maillages (triangles, normales)

function pushTri(out, p, n) {
  for (let k = 0; k < 3; k++) out.push(p[k][0], p[k][1], p[k][2], n[k][0], n[k][1], n[k][2]);
}

function sphereMesh(cx, cy, cz, R) {
  const out = [];
  const nt = 18;
  const np = 36;
  const P = (i, j) => {
    const th = (Math.PI * i) / nt;
    const ph = (2 * Math.PI * j) / np;
    const n = [Math.sin(th) * Math.cos(ph), Math.sin(th) * Math.sin(ph), Math.cos(th)];
    return [[cx + R * n[0], cy + R * n[1], cz + R * n[2]], n];
  };
  for (let i = 0; i < nt; i++)
    for (let j = 0; j < np; j++) {
      const a = P(i, j);
      const b = P(i + 1, j);
      const c = P(i + 1, j + 1);
      const d = P(i, j + 1);
      pushTri(out, [a[0], b[0], c[0]], [a[1], b[1], c[1]]);
      pushTri(out, [a[0], c[0], d[0]], [a[1], c[1], d[1]]);
    }
  return new Float32Array(out);
}

function boxMesh(cx, cy, cz, R, ca, sa) {
  const out = [];
  const tr = ([x, y, z]) => [cx + x * ca + z * sa, cy + y, cz - x * sa + z * ca];
  const faces = [
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[-1, 0, 0], [0, 0, 1], [0, 1, 0]],
    [[0, 1, 0], [0, 0, 1], [1, 0, 0]],
    [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
    [[0, 0, 1], [1, 0, 0], [0, 1, 0]],
    [[0, 0, -1], [0, 1, 0], [1, 0, 0]],
  ];
  for (const [n, u, v] of faces) {
    const c = (s, t) => tr([R * (n[0] + s * u[0] + t * v[0]), R * (n[1] + s * u[1] + t * v[1]), R * (n[2] + s * u[2] + t * v[2])]);
    const nn = [n[0] * ca + n[2] * sa, n[1], -n[0] * sa + n[2] * ca];
    const N3 = [nn, nn, nn];
    pushTri(out, [c(-1, -1), c(1, -1), c(1, 1)], N3);
    pushTri(out, [c(-1, -1), c(1, 1), c(-1, 1)], N3);
  }
  return new Float32Array(out);
}

function cylinderMesh(cx, cz, R, ny) {
  const out = [];
  const n = 48;
  for (let k = 0; k < n; k++) {
    const a0 = (2 * Math.PI * k) / n;
    const a1 = (2 * Math.PI * (k + 1)) / n;
    const n0 = [Math.cos(a0), 0, Math.sin(a0)];
    const n1 = [Math.cos(a1), 0, Math.sin(a1)];
    const p = (nn, y) => [cx + R * nn[0], y, cz + R * nn[2]];
    pushTri(out, [p(n0, 0), p(n1, 0), p(n1, ny - 1)], [n0, n1, n1]);
    pushTri(out, [p(n0, 0), p(n1, ny - 1), p(n0, ny - 1)], [n0, n1, n0]);
  }
  return new Float32Array(out);
}

function wingMesh(cx, cy, cz, L, half, t, ca, sa) {
  const out = [];
  const m = 40;
  const sec = [];
  for (let k = 0; k <= m; k++) {
    const x = 0.5 * (1 - Math.cos((Math.PI * k) / m));
    sec.push([x, naca00(t, x)]);
  }
  // Contour fermé : extrados du bord de fuite au bord d'attaque, puis intrados.
  const loop = [];
  for (let k = m; k >= 0; k--) loop.push([sec[k][0], sec[k][1]]);
  for (let k = 1; k <= m; k++) loop.push([sec[k][0], -sec[k][1]]);
  const tr = (xc, zc, y) => {
    const xb = (xc - 0.25) * L;
    const zb = zc * L;
    return [cx + xb * ca + zb * sa, y, cz - xb * sa + zb * ca];
  };
  const rot = (nx, nz) => [nx * ca + nz * sa, 0, -nx * sa + nz * ca];
  for (let k = 0; k < loop.length - 1; k++) {
    const [x0, z0] = loop[k];
    const [x1, z1] = loop[k + 1];
    const tx = x1 - x0;
    const tz = z1 - z0;
    const l = Math.hypot(tx, tz) || 1;
    // Contour parcouru dans le sens horaire (vu de +y) : normale extérieure (−t_z, t_x).
    const nrm = rot(-tz / l, tx / l);
    const a = tr(x0, z0, cy - half);
    const b = tr(x1, z1, cy - half);
    const c = tr(x1, z1, cy + half);
    const d = tr(x0, z0, cy + half);
    pushTri(out, [a, b, c], [nrm, nrm, nrm]);
    pushTri(out, [a, c, d], [nrm, nrm, nrm]);
  }
  // Saumons.
  for (const s of [-1, 1]) {
    const y = cy + s * half;
    const nrm = [0, s, 0];
    const c0 = tr(0.35, 0, y);
    for (let k = 0; k < loop.length - 1; k++) {
      pushTri(out, [c0, tr(loop[k][0], loop[k][1], y), tr(loop[k + 1][0], loop[k + 1][1], y)], [nrm, nrm, nrm]);
    }
  }
  return new Float32Array(out);
}
