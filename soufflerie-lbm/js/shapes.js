/**
 * Formes de maquettes. Chaque forme est décrite dans son repère propre,
 * adimensionnée par sa longueur de référence L, origine au point de réduction
 * des moments. `makeBody` la place dans le domaine (coordonnées réseau),
 * avec une incidence α (positive = cabrée, bord d'attaque vers le haut).
 */

const TAU = Math.PI * 2;

/** Profil NACA 4 chiffres : épaisseur et cambrure en fonction de x ∈ [0, 1]. */
export function parseNaca(code) {
  const s = String(code ?? '').trim();
  if (!/^\d{4}$/.test(s)) return null;
  const m = Number(s[0]) / 100;
  const p = Number(s[1]) / 10;
  const t = Number(s.slice(2)) / 100;
  if (t < 0.04 || t > 0.4) return null;
  if (m > 0 && (p <= 0 || p >= 1)) return null;
  return { m, p, t, code: s };
}

function nacaThickness(t, x) {
  // Bord de fuite fermé (coefficient −0.1036).
  return 5 * t * (0.2969 * Math.sqrt(x) - 0.126 * x - 0.3516 * x * x + 0.2843 * x ** 3 - 0.1036 * x ** 4);
}

function nacaCamber(m, p, x) {
  if (m === 0) return 0;
  return x < p ? (m / (p * p)) * (2 * p * x - x * x) : (m / ((1 - p) ** 2)) * (1 - 2 * p + 2 * p * x - x * x);
}

/**
 * Position de la pointe du bord de fuite par rapport au réseau (parties
 * fractionnaires x, y) : abscisse d'un nœud, à mi-hauteur entre deux nœuds,
 * comme à α = 0 (symétrie haut/bas). Parmi les positions essayées, c'est
 * celle qui donne la portance la plus régulière en fonction de α. Elle peut
 * exciter près de la pointe un mode alterné pair/impair : les efforts sont
 * donc toujours moyennés sur plusieurs pas (voir WindTunnel).
 */
const TE_ANCHOR = [0, 0.5];

/** Triangle isocèle, pointe vers l'amont, hauteur frontale 1, centré sur son centre de gravité. */
const TRI_LEN = Math.sqrt(3) / 2;
const TRI = (() => {
  const v = [
    [-TRI_LEN / 2, 0],
    [TRI_LEN / 2, 0.5],
    [TRI_LEN / 2, -0.5],
  ];
  const gx = (v[0][0] + v[1][0] + v[2][0]) / 3;
  return v.map(([x, y]) => [x - gx, y]);
})();

function insideTriangle(x, y) {
  let sign = 0;
  for (let k = 0; k < 3; k++) {
    const [ax, ay] = TRI[k];
    const [bx, by] = TRI[(k + 1) % 3];
    const c = (bx - ax) * (y - ay) - (by - ay) * (x - ax);
    if (c !== 0) {
      const s = c > 0 ? 1 : -1;
      if (sign === 0) sign = s;
      else if (s !== sign) return false;
    }
  }
  return true;
}

function circle(rx, ry, n = 96) {
  const pts = [];
  for (let k = 0; k < n; k++) {
    const a = (k / n) * TAU;
    pts.push([rx * Math.cos(a), ry * Math.sin(a)]);
  }
  return pts;
}

/**
 * Catalogue. `alpha` : [min, max, défaut] en degrés, ou null si sans objet.
 * `ref` : nom de la longueur de référence ; `momentRef` : point de réduction.
 */
export const SHAPES = [
  {
    id: 'cylinder',
    name: 'Cylindre',
    ref: 'diamètre D',
    momentRef: 'au centre',
    alpha: null,
    spin: true,
    inside: (x, y) => x * x + y * y <= 0.25,
    outline: () => circle(0.5, 0.5),
  },
  {
    id: 'square',
    name: 'Carré',
    ref: 'côté a',
    momentRef: 'au centre',
    alpha: [-45, 45, 0],
    inside: (x, y) => Math.abs(x) <= 0.5 && Math.abs(y) <= 0.5,
    outline: () => [
      [-0.5, -0.5],
      [0.5, -0.5],
      [0.5, 0.5],
      [-0.5, 0.5],
    ],
  },
  {
    id: 'ellipse',
    name: 'Ellipse',
    ref: 'grand axe c',
    momentRef: 'au centre',
    alpha: [-90, 90, 0],
    thickness: true,
    inside: (x, y, p) => {
      const b = 0.5 * p.thickness;
      return (x / 0.5) ** 2 + (y / b) ** 2 <= 1;
    },
    outline: (p) => circle(0.5, 0.5 * p.thickness),
  },
  {
    id: 'plate',
    name: 'Plaque plane',
    trailingEdge: [0.5, 0],
    short: 'Plaque',
    ref: 'corde c',
    momentRef: 'au centre',
    alpha: [-90, 90, 10],
    inside: (x, y, p) => Math.abs(x) <= 0.5 && Math.abs(y) <= p.plateHalf,
    outline: (p) => [
      [-0.5, -p.plateHalf],
      [0.5, -p.plateHalf],
      [0.5, p.plateHalf],
      [-0.5, p.plateHalf],
    ],
  },
  {
    id: 'naca',
    name: 'Profil NACA',
    short: 'NACA',
    ref: 'corde c',
    momentRef: 'au quart de corde',
    alpha: [-20, 25, 4],
    naca: true,
    trailingEdge: [0.75, 0],
    inside: (x, y, p) => {
      const xc = x + 0.25;
      if (xc < 0 || xc > 1) return false;
      const { m, p: pp, t } = p.nacaParams;
      return Math.abs(y - nacaCamber(m, pp, xc)) <= nacaThickness(t, xc);
    },
    outline: (p) => {
      const { m, p: pp, t } = p.nacaParams;
      const n = 80;
      const upper = [];
      const lower = [];
      for (let k = 0; k <= n; k++) {
        const xc = 0.5 * (1 - Math.cos((k / n) * Math.PI));
        const yc = nacaCamber(m, pp, xc);
        const yt = nacaThickness(t, xc);
        upper.push([xc - 0.25, yc + yt]);
        lower.push([xc - 0.25, yc - yt]);
      }
      return upper.reverse().concat(lower.slice(1));
    },
  },
  {
    id: 'triangle',
    name: 'Triangle',
    ref: 'hauteur frontale h',
    momentRef: 'au centre de gravité',
    alpha: [-180, 180, 0],
    inside: insideTriangle,
    outline: () => TRI.map((v) => v.slice()),
  },
];

export const SHAPE_BY_ID = Object.fromEntries(SHAPES.map((s) => [s.id, s]));

/**
 * Place une forme dans le domaine.
 * @param {string} id
 * @param {object} o
 * @param {number} o.L longueur de référence (mailles)
 * @param {number} o.x0 position du point de référence
 * @param {number} o.y0
 * @param {number} [o.alpha=0] incidence (degrés)
 * @param {string} [o.naca='0012']
 * @param {number} [o.thickness=0.4] épaisseur relative (ellipse)
 * @param {number} [o.spinRatio=0] cylindre : vitesse périphérique / U (>0 : portance positive)
 * @param {number} [o.uIn=0.1] vitesse amont (pour la rotation)
 */
export function makeBody(id, o) {
  const shape = SHAPE_BY_ID[id] ?? SHAPES[0];
  const L = o.L;
  const alphaDeg = shape.alpha ? o.alpha ?? 0 : 0;
  const a = (alphaDeg * Math.PI) / 180;
  const ca = Math.cos(a);
  const sa = Math.sin(a);
  const p = {
    thickness: Math.min(1, Math.max(0.1, o.thickness ?? 0.4)),
    plateHalf: Math.max(0.03, 1.25 / L),
    nacaParams: parseNaca(o.naca) ?? parseNaca('0012'),
  };
  let x0 = o.x0;
  let y0 = o.y0;
  // Bord de fuite aigu : décalage de moins d'une demi-maille pour que sa pointe
  // occupe toujours la même position par rapport au réseau, quelle que soit
  // l'incidence (x entier, y demi-entier, comme à α = 0). Sans cela, la
  // position effective du bord de fuite saute d'une maille à l'autre et la
  // portance varie de façon erratique avec α (≈ ±1° d'incidence effective à
  // 16 mailles). Les efforts ne dépendent pas de cette translation.
  if (shape.trailingEdge && o.anchorTrailingEdge !== false) {
    const [bx, by] = shape.trailingEdge;
    const [ax, ay] = o.teAnchor ?? TE_ANCHOR;
    const tx = x0 + L * (ca * bx + sa * by);
    const ty = y0 + L * (-sa * bx + ca * by);
    x0 += Math.floor(tx - ax + 0.5) + ax - tx;
    y0 += Math.floor(ty - ay + 0.5) + ay - ty;
  }
  const test = shape.inside;
  const inside = (x, y) => {
    const dx = (x - x0) / L;
    const dy = (y - y0) / L;
    return test(ca * dx - sa * dy, sa * dx + ca * dy, p);
  };
  const outline = shape.outline(p).map(([bx, by]) => [x0 + L * (ca * bx + sa * by), y0 + L * (-sa * bx + ca * by)]);
  let ymin = Infinity;
  let ymax = -Infinity;
  for (const [, y] of outline) {
    ymin = Math.min(ymin, y);
    ymax = Math.max(ymax, y);
  }
  // Rotation horaire pour spinRatio > 0 : portance positive (effet Magnus).
  const spin = shape.spin && o.spinRatio ? (-(o.spinRatio ?? 0) * (o.uIn ?? 0.1)) / (L / 2) : 0;
  return {
    shape,
    inside,
    outline,
    L,
    refX: x0,
    refY: y0,
    alpha: alphaDeg,
    spin,
    frontal: ymax - ymin,
    label: describe(shape, alphaDeg, p, o),
  };
}

function describe(shape, alpha, p, o) {
  let s = shape.name;
  if (shape.naca) s = `NACA ${p.nacaParams.code}`;
  if (shape.thickness) s += ` e/c = ${p.thickness.toFixed(2)}`;
  if (shape.spin && o.spinRatio) s += ` ΩR/U = ${Number(o.spinRatio).toFixed(1)}`;
  return s;
}

/** Chemin SVG (repère propre, y vers le haut) pour les vignettes de l'interface. */
export function shapeIconPath(id, size = 24) {
  const shape = SHAPE_BY_ID[id];
  const p = { thickness: 0.4, plateHalf: 0.06, nacaParams: parseNaca('2412') };
  const pts = shape.outline(p);
  const k = size * 0.8;
  const c = size / 2;
  let cx = 0;
  if (id === 'naca') cx = 0.25; // milieu de corde
  return (
    pts
      .map(([x, y], i) => `${i ? 'L' : 'M'}${(c + (x - cx) * k).toFixed(2)} ${(c - y * k).toFixed(2)}`)
      .join(' ') + ' Z'
  );
}
