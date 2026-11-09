/**
 * Historique des coefficients aérodynamiques et analyse : moyennes sur un
 * nombre entier de périodes quand l'écoulement est périodique, amplitude,
 * nombre de Strouhal.
 */

export class ForceHistory {
  constructor(capacity = 8192) {
    this.cap = capacity;
    this.t = new Float64Array(capacity);
    this.cx = new Float32Array(capacity);
    this.cz = new Float32Array(capacity);
    this.cm = new Float32Array(capacity);
    this.clear();
  }

  clear() {
    this.head = 0;
    this.count = 0;
  }

  push(t, cx, cz, cm) {
    const k = this.head;
    this.t[k] = t;
    this.cx[k] = cx;
    this.cz[k] = cz;
    this.cm[k] = cm;
    this.head = (k + 1) % this.cap;
    if (this.count < this.cap) this.count++;
  }

  get lastT() {
    return this.count ? this.t[(this.head - 1 + this.cap) % this.cap] : -Infinity;
  }

  /** Échantillons (ordre chronologique) avec t ≥ tMin. */
  since(tMin = -Infinity) {
    const n = this.count;
    const start = (this.head - n + this.cap) % this.cap;
    // Recherche du premier indice (les temps sont croissants).
    let lo = 0;
    let hi = n;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (this.t[(start + mid) % this.cap] < tMin) lo = mid + 1;
      else hi = mid;
    }
    const m = n - lo;
    const out = {
      t: new Float64Array(m),
      cx: new Float32Array(m),
      cz: new Float32Array(m),
      cm: new Float32Array(m),
    };
    for (let k = 0; k < m; k++) {
      const s = (start + lo + k) % this.cap;
      out.t[k] = this.t[s];
      out.cx[k] = this.cx[s];
      out.cz[k] = this.cz[s];
      out.cm[k] = this.cm[s];
    }
    return out;
  }
}

/** Écart-type minimal pour parler d'oscillation (en unités de coefficient). */
const MIN_STD = 0.006;

function meanStd(a, i0, i1) {
  let s = 0;
  let s2 = 0;
  const n = i1 - i0;
  if (n <= 0) return { mean: NaN, std: NaN, min: NaN, max: NaN };
  let min = Infinity;
  let max = -Infinity;
  for (let k = i0; k < i1; k++) {
    const v = a[k];
    s += v;
    s2 += v * v;
    if (v < min) min = v;
    if (v > max) max = v;
  }
  const mean = s / n;
  return { mean, std: Math.sqrt(Math.max(0, s2 / n - mean * mean)), min, max };
}

/**
 * Analyse d'une fenêtre d'échantillons.
 * @param {{t:Float64Array, cx:Float32Array, cz:Float32Array, cm:Float32Array}} w
 * @param {number} U vitesse amont (réseau)
 * @param {number} L longueur de référence (mailles)
 */
export function analyze(w, U, L) {
  const n = w.t.length;
  const res = {
    n,
    duration: n > 1 ? ((w.t[n - 1] - w.t[0]) * U) / L : 0,
    cx: NaN,
    cz: NaN,
    cm: NaN,
    cxAmp: 0,
    czAmp: 0,
    st: NaN,
    periodic: false,
    periods: 0,
    drift: NaN,
  };
  if (n < 8) return res;

  // Détection d'oscillation sur Cz (à défaut sur Cx), avec hystérésis.
  const z = meanStd(w.cz, 0, n);
  const x = meanStd(w.cx, 0, n);
  let sig = w.cz;
  let ref = z;
  if (z.std < MIN_STD && x.std > MIN_STD) {
    sig = w.cx;
    ref = x;
  }
  const crossings = [];
  if (ref.std > MIN_STD) {
    const h = 0.25 * ref.std;
    let state = 0;
    for (let k = 1; k < n; k++) {
      const v = sig[k] - ref.mean;
      if (v > h && state <= 0) {
        if (state === -1) {
          // Interpolation linéaire du passage par la moyenne.
          let j = k;
          while (j > 0 && sig[j - 1] - ref.mean > 0) j--;
          const a = sig[j - 1] - ref.mean;
          const b = sig[j] - ref.mean;
          const tt = j > 0 && b !== a ? w.t[j - 1] + ((w.t[j] - w.t[j - 1]) * -a) / (b - a) : w.t[k];
          crossings.push({ t: tt, k: j });
        }
        state = 1;
      } else if (v < -h && state >= 0) {
        state = -1;
      }
    }
  }

  // Régularité des périodes : rejette le bruit et les dérives lentes.
  let regular = false;
  if (crossings.length >= 3) {
    const dts = [];
    for (let k = 1; k < crossings.length; k++) dts.push(crossings[k].t - crossings[k - 1].t);
    const m = dts.reduce((a, b) => a + b, 0) / dts.length;
    const sd = Math.sqrt(dts.reduce((a, b) => a + (b - m) ** 2, 0) / dts.length);
    regular = sd < 0.2 * m;
  }

  let i0 = 0;
  let i1 = n;
  if (regular) {
    const first = crossings[0];
    const last = crossings[crossings.length - 1];
    const periods = crossings.length - 1;
    const T = (last.t - first.t) / periods; // en pas de temps
    res.periodic = true;
    res.periods = periods;
    const f = 1 / T;
    // Si le signal analysé est Cx (détachement symétrique), sa fréquence est 2 f.
    res.st = ((sig === w.cx ? f / 2 : f) * L) / U;
    i0 = first.k;
    i1 = last.k;
  }

  const sx = meanStd(w.cx, i0, i1);
  const sz = meanStd(w.cz, i0, i1);
  const sm = meanStd(w.cm, i0, i1);
  res.cx = sx.mean;
  res.cz = sz.mean;
  res.cm = sm.mean;
  res.cxAmp = 0.5 * (sx.max - sx.min);
  res.czAmp = 0.5 * (sz.max - sz.min);

  // Dérive : écart entre les moyennes des deux moitiés de la fenêtre.
  if (!res.periodic) {
    const h = (n / 2) | 0;
    const a = meanStd(w.cx, 0, h).mean;
    const b = meanStd(w.cx, h, n).mean;
    const c = meanStd(w.cz, 0, h).mean;
    const d = meanStd(w.cz, h, n).mean;
    res.drift = Math.max(Math.abs(a - b), Math.abs(c - d));
  } else {
    res.drift = 0;
  }
  return res;
}
