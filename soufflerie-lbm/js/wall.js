/**
 * Distributions pariétales du calcul LBM : coefficient de pression Cp et de
 * frottement Cf le long du corps, moyennés dans le temps. Fonctionne avec les
 * deux moteurs : seuls ρ et u des mailles voisines de la paroi sont lus.
 */
import { spaldingUtau } from './lbm.js';

/** Formes pour lesquelles une abscisse le long de la paroi a un sens simple. */
export const WALL_SHAPES = ['naca', 'ellipse', 'plate', 'cylinder'];

export class WallSampler {
  /** @param {import('./sim.js').WindTunnel} wt */
  constructor(wt) {
    this.wt = wt;
    const lbm = wt.lbm;
    const geo = lbm.geo;
    const body = wt.body;
    const nx = lbm.nx;
    this.kind = body.shape.id === 'cylinder' ? 'angle' : 'chord';
    const loc = body.outlineLocal;
    const m = loc.length;
    let iLE = 0;
    let xmin = Infinity;
    let xmax = -Infinity;
    for (let k = 0; k < m; k++) {
      if (loc[k][0] < xmin) {
        xmin = loc[k][0];
        iLE = k;
      }
      xmax = Math.max(xmax, loc[k][0]);
    }
    const span = xmax - xmin || 1;
    const byIndex = body.shape.id === 'naca';
    const a = (body.alpha * Math.PI) / 180;
    // Axe x du corps dans le repère du réseau (voir makeBody).
    this.ex = Math.cos(a);
    this.ey = -Math.sin(a);
    const cells = [];
    let x0 = Infinity;
    let x1 = -Infinity;
    let y0 = Infinity;
    let y1 = -Infinity;
    for (const idx of geo.wallCells()) {
      const pos = geo.nearPos[idx];
      if (pos < 0) continue;
      const k = Math.min(m - 1, Math.floor(pos));
      const t = pos - k;
      const [ax, ay] = loc[k];
      const [bx2, by2] = loc[(k + 1) % m];
      const bx = ax + t * (bx2 - ax);
      const by = ay + t * (by2 - ay);
      const upper = byIndex ? k < iLE || (k === iLE && t === 0) : by >= 0;
      const s =
        this.kind === 'angle' ? 180 - (Math.abs(Math.atan2(by, bx)) * 180) / Math.PI : (bx - xmin) / span;
      const x = idx % nx;
      const y = (idx - x) / nx;
      cells.push({ idx, x, y, s, upper, d: geo.dist[idx], nx: geo.normal[2 * idx], ny: geo.normal[2 * idx + 1] });
      x0 = Math.min(x0, x);
      x1 = Math.max(x1, x);
      y0 = Math.min(y0, y);
      y1 = Math.max(y1, y);
    }
    this.cells = cells;
    this.region = cells.length ? { x0, y0, w: x1 - x0 + 1, h: y1 - y0 + 1 } : null;
    this.sumCp = new Float64Array(cells.length);
    this.sumCf = new Float64Array(cells.length);
    this.n = 0;
  }

  reset() {
    this.sumCp.fill(0);
    this.sumCf.fill(0);
    this.n = 0;
  }

  /**
   * Ajoute un échantillon.
   * @param {{rho, ux, uy}} f champs ; at(x, y) donne l'indice d'une maille
   * @param {(x:number, y:number) => number} at
   * @param {number} rhoRef masse volumique de référence (amont)
   */
  add(f, at, rhoRef) {
    const wt = this.wt;
    const lbm = wt.lbm;
    const geo = lbm.geo;
    const U = wt.uIn;
    const nu = wt.viscosityFor(wt.re);
    const sa = !!lbm.saActive;
    const q = 0.5 * rhoRef * U * U;
    const { ex, ey } = this;
    for (let k = 0; k < this.cells.length; k++) {
      const c = this.cells[k];
      const i = at(c.x, c.y);
      const rho = f.rho[i];
      const ux = f.ux[i];
      const uy = f.uy[i];
      // p = ρ c_s², c_s² = 1/3.
      this.sumCp[k] += (rho - rhoRef) / (3 * q);
      const un = ux * c.nx + uy * c.ny;
      const tx = ux - un * c.nx;
      const ty = uy - un * c.ny;
      const ut = Math.hypot(tx, ty);
      const sign = tx * ex + ty * ey >= 0 ? 1 : -1;
      let tau;
      if (sa) {
        const cfLam = geo.lam[c.idx];
        if (cfLam > 0) tau = 0.5 * cfLam * ut * ut;
        else {
          const ut2 = spaldingUtau(ut, Math.max(0.1, c.d), nu);
          tau = ut2 * ut2;
        }
      } else {
        // Profil linéaire entre la paroi (u = 0) et la maille.
        tau = (nu * ut) / Math.max(0.3, c.d);
      }
      this.sumCf[k] += (sign * rho * tau) / q;
    }
    this.n++;
  }

  /** Moyennes, triées par abscisse, séparées en faces. */
  result() {
    const out = { kind: this.kind, n: this.n, upper: [], lower: [] };
    if (!this.n) return out;
    for (let k = 0; k < this.cells.length; k++) {
      const c = this.cells[k];
      const p = { s: c.s, cp: this.sumCp[k] / this.n, cf: this.sumCf[k] / this.n, d: c.d };
      (c.upper ? out.upper : out.lower).push(p);
    }
    out.upper.sort((a, b) => a.s - b.s);
    out.lower.sort((a, b) => a.s - b.s);
    return out;
  }
}

/** Masse volumique moyenne d'une colonne amont (référence de pression). */
export function upstreamDensity(rhoColumn) {
  const n = rhoColumn.length;
  let s = 0;
  let k = 0;
  for (let j = Math.floor(n / 4); j < Math.ceil((3 * n) / 4); j++) {
    s += rhoColumn[j];
    k++;
  }
  return k ? s / k : 1;
}
