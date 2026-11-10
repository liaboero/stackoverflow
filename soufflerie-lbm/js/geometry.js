/**
 * Géométrie du corps sur le réseau, commune aux moteurs CPU et GPU :
 * mailles solides, liens de paroi interpolés (fraction q de Bouzidi, y compris
 * à travers les parois plus minces qu'une maille), distance à la paroi et
 * amortissement pariétal de la viscosité de sous-maille.
 */

export const CX = [0, 1, 0, -1, 0, 1, -1, -1, 1];
export const CY = [0, 0, 1, 0, -1, 1, 1, -1, -1];

export const FLUID = 0;
export const SOLID = 1;
/** Fluide traité par le chemin lent : voisin d'un solide ou ligne haute/basse. */
export const NEAR_WALL = 2;

export class Geometry {
  constructor(nx, ny, lesWallA = 4) {
    this.nx = nx;
    this.ny = ny;
    const N = (this.N = nx * ny);
    this.type = new Uint8Array(N);
    /** q du lien tiré dans la direction j : wallQ[j·N + idx] ; 0 si pas de paroi. */
    this.wallQ = new Float32Array(9 * N);
    this.wallDist = new Float32Array(N);
    /** Distance (mailles) du centre de la maille à la paroi : exacte près du corps. */
    this.dist = new Float32Array(N).fill(1e6);
    /** Normale unitaire paroi → maille (près du corps), pour la loi de paroi. */
    this.normal = new Float32Array(2 * N);
    /** Position k + t du point de paroi le plus proche sur le contour (−1 : loin). */
    this.nearPos = new Float32Array(N).fill(-1);
    /**
     * Zone laminaire imposée (transition) : > 0 dans la zone, où la valeur est
     * le frottement laminaire Cf fourni par la couche limite intégrale ; 0
     * ailleurs (loi de paroi turbulente).
     */
    this.lam = new Float32Array(N);
    /** Facteur d'amortissement de la LES (0 à la paroi → 1 au large). */
    this.lesDamp = new Float32Array(N).fill(1);
    this.lesWallA = lesWallA;
    this.solidCount = 0;
    this.linkCount = 0;
    this.thinLinkCount = 0;
    // Lignes haute et basse : chemin lent (réflexion spéculaire).
    for (let x = 0; x < nx; x++) {
      this.type[x] = NEAR_WALL;
      this.type[(ny - 1) * nx + x] = NEAR_WALL;
    }
  }

  /**
   * Construit la géométrie.
   * @param {(x:number, y:number) => boolean} inside test d'appartenance
   * @param {{x0:number,x1:number,y0:number,y1:number}} [bbox] boîte englobante
   *        du corps, où chercher les parois minces
   * @returns {number[]} indices des mailles solides redevenues fluides
   */
  /**
   * @param {Array<[number, number]>} [outline] contour du corps (coordonnées
   *        réseau), pour les distances exactes et les normales de paroi
   */
  build(inside, bbox, outline) {
    const { nx, ny, N, type, wallQ } = this;
    const revived = [];
    let solids = 0;
    // Zone laminaire liée à l'ancienne géométrie : à redéfinir (setTransition).
    this.lam.fill(0);
    for (let y = 0; y < ny; y++) {
      for (let x = 0; x < nx; x++) {
        const idx = y * nx + x;
        const solid = x >= 2 && x <= nx - 3 && y >= 2 && y <= ny - 3 && inside(x, y);
        if (solid) {
          type[idx] = SOLID;
          solids++;
        } else {
          if (type[idx] === SOLID) revived.push(idx);
          type[idx] = FLUID;
        }
      }
    }
    this.solidCount = solids;

    const bx0 = bbox ? Math.floor(bbox.x0) - 2 : 0;
    const bx1 = bbox ? Math.ceil(bbox.x1) + 2 : -1;
    const by0 = bbox ? Math.floor(bbox.y0) - 2 : 0;
    const by1 = bbox ? Math.ceil(bbox.y1) + 2 : -1;

    wallQ.fill(0);
    let links = 0;
    let thin = 0;
    for (let y = 0; y < ny; y++) {
      for (let x = 1; x < nx - 1; x++) {
        const idx = y * nx + x;
        if (type[idx] === SOLID) continue;
        let near = y === 0 || y === ny - 1;
        const inBox = x >= bx0 && x <= bx1 && y >= by0 && y <= by1;
        for (let j = 1; j < 9; j++) {
          const sx = x - CX[j];
          const sy = y - CY[j];
          if (sy < 0 || sy >= ny) continue;
          // Lien du nœud fluide vers la paroi : direction −c_j.
          let q = 0;
          if (type[sy * nx + sx] === SOLID) {
            q = linkFraction(inside, x, y, -CX[j], -CY[j]);
          } else if (inBox) {
            // Lien fluide–fluide qui traverse une paroi mince.
            q = crossingFraction(inside, x, y, -CX[j], -CY[j]);
            if (q > 0) thin++;
          }
          if (q > 0) {
            near = true;
            wallQ[j * N + idx] = q;
            links++;
          }
        }
        if (near) type[idx] = NEAR_WALL;
      }
    }
    // Colonnes d'entrée et de sortie : traitées à part par les solveurs.
    for (let y = 0; y < ny; y++) {
      type[y * nx] = FLUID;
      type[y * nx + nx - 1] = FLUID;
    }
    this.linkCount = links;
    this.thinLinkCount = thin;
    this._computeWallDistance();
    this._exactWallDistance(outline, bbox);
    this.updateLesDamping();
    return revived;
  }

  /** Mailles portant au moins un lien de paroi avec le corps (calcul des efforts). */
  wallCells() {
    const { N, wallQ, type } = this;
    const out = [];
    for (let idx = 0; idx < N; idx++) {
      if (type[idx] !== NEAR_WALL) continue;
      for (let j = 1; j < 9; j++) {
        if (wallQ[j * N + idx] > 0) {
          out.push(idx);
          break;
        }
      }
    }
    return out;
  }

  setLesWallLength(A) {
    this.lesWallA = A;
    this.updateLesDamping();
  }

  /**
   * Facteur d'amortissement de type Van Driest, (1 − e^(−d/A))², où d est la
   * distance (en mailles) à la paroi du corps. Le modèle de Smagorinsky
   * ajoute sinon de la viscosité dans la couche limite laminaire, où le
   * cisaillement est fort mais l'écoulement résolu, et surestime la traînée
   * de frottement.
   */
  updateLesDamping() {
    const { N, lesDamp } = this;
    const A = this.lesWallA;
    if (!(A > 0) || this.solidCount === 0) {
      lesDamp.fill(1);
      return;
    }
    const d = this.dist;
    for (let i = 0; i < N; i++) {
      const f = 1 - Math.exp(-d[i] / A);
      lesDamp[i] = f * f;
    }
  }

  /**
   * Distance exacte au contour polygonal dans une boîte autour du corps (les
   * modèles de paroi et de turbulence y sont sensibles), distance de chanfrein
   * corrigée ailleurs. Normale et point de paroi le plus proche dans la boîte.
   */
  _exactWallDistance(outline, bbox) {
    const { nx, ny, N, type, wallDist, dist, normal, nearPos } = this;
    for (let i = 0; i < N; i++) dist[i] = type[i] === SOLID ? 0 : this.solidCount ? Math.max(0.5, wallDist[i] - 0.5) : 1e6;
    normal.fill(0);
    nearPos.fill(-1);
    if (!outline || outline.length < 2 || !bbox) return;
    const m = outline.length;
    const ax = new Float64Array(m);
    const ay = new Float64Array(m);
    const ex = new Float64Array(m);
    const ey = new Float64Array(m);
    const l2 = new Float64Array(m);
    for (let k = 0; k < m; k++) {
      const [x0, y0] = outline[k];
      const [x1, y1] = outline[(k + 1) % m];
      ax[k] = x0;
      ay[k] = y0;
      ex[k] = x1 - x0;
      ey[k] = y1 - y0;
      l2[k] = Math.max(1e-12, ex[k] * ex[k] + ey[k] * ey[k]);
    }
    const margin = Math.max(8, 0.6 * Math.max(bbox.x1 - bbox.x0, bbox.y1 - bbox.y0));
    const X0 = Math.max(0, Math.floor(bbox.x0 - margin));
    const X1 = Math.min(nx - 1, Math.ceil(bbox.x1 + margin));
    const Y0 = Math.max(0, Math.floor(bbox.y0 - margin));
    const Y1 = Math.min(ny - 1, Math.ceil(bbox.y1 + margin));
    for (let y = Y0; y <= Y1; y++) {
      for (let x = X0; x <= X1; x++) {
        const idx = y * nx + x;
        if (type[idx] === SOLID) continue;
        let best = Infinity;
        let bk = 0;
        let bt = 0;
        for (let k = 0; k < m; k++) {
          let t = ((x - ax[k]) * ex[k] + (y - ay[k]) * ey[k]) / l2[k];
          t = t < 0 ? 0 : t > 1 ? 1 : t;
          const dx = x - ax[k] - t * ex[k];
          const dy = y - ay[k] - t * ey[k];
          const d2 = dx * dx + dy * dy;
          if (d2 < best) {
            best = d2;
            bk = k;
            bt = t;
          }
        }
        const d = Math.sqrt(best);
        dist[idx] = Math.max(d, 1e-3);
        const qx = ax[bk] + bt * ex[bk];
        const qy = ay[bk] + bt * ey[bk];
        if (d > 1e-6) {
          normal[2 * idx] = (x - qx) / d;
          normal[2 * idx + 1] = (y - qy) / d;
        } else {
          // Point sur la paroi : normale du segment.
          const len = Math.sqrt(l2[bk]);
          normal[2 * idx] = ey[bk] / len;
          normal[2 * idx + 1] = -ex[bk] / len;
        }
        nearPos[idx] = bk + bt;
      }
    }
  }

  /** Distance de chanfrein (1, √2) aux mailles solides, en deux passes. */
  _computeWallDistance() {
    const { nx, ny, N, type } = this;
    const d = this.wallDist;
    const BIG = 1e6;
    for (let i = 0; i < N; i++) d[i] = type[i] === SOLID ? 0 : BIG;
    const D = Math.SQRT2;
    for (let y = 0; y < ny; y++) {
      for (let x = 0; x < nx; x++) {
        const i = y * nx + x;
        let v = d[i];
        if (x > 0) v = Math.min(v, d[i - 1] + 1);
        if (y > 0) {
          v = Math.min(v, d[i - nx] + 1);
          if (x > 0) v = Math.min(v, d[i - nx - 1] + D);
          if (x < nx - 1) v = Math.min(v, d[i - nx + 1] + D);
        }
        d[i] = v;
      }
    }
    for (let y = ny - 1; y >= 0; y--) {
      for (let x = nx - 1; x >= 0; x--) {
        const i = y * nx + x;
        let v = d[i];
        if (x < nx - 1) v = Math.min(v, d[i + 1] + 1);
        if (y < ny - 1) {
          v = Math.min(v, d[i + nx] + 1);
          if (x < nx - 1) v = Math.min(v, d[i + nx + 1] + D);
          if (x > 0) v = Math.min(v, d[i + nx - 1] + D);
        }
        d[i] = v;
      }
    }
  }
}

/**
 * Pour un lien entre deux nœuds fluides : fraction q de la première
 * intersection avec le corps, ou 0 si le lien ne le traverse pas.
 */
function crossingFraction(inside, x, y, dx, dy) {
  const n = 16;
  for (let k = 1; k < n; k++) {
    const t = k / n;
    if (inside(x + t * dx, y + t * dy)) {
      let lo = (k - 1) / n;
      let hi = t;
      for (let it = 0; it < 10; it++) {
        const mid = 0.5 * (lo + hi);
        if (inside(x + mid * dx, y + mid * dy)) hi = mid;
        else lo = mid;
      }
      return Math.min(1, Math.max(0.02, 0.5 * (lo + hi)));
    }
  }
  return 0;
}

/**
 * Fraction q ∈ ]0, 1] du lien (x, y) → (x + dx, y + dy) située côté fluide,
 * obtenue par dichotomie sur le test d'appartenance.
 */
function linkFraction(inside, x, y, dx, dy) {
  let lo = 0;
  let hi = 1;
  for (let k = 0; k < 14; k++) {
    const mid = 0.5 * (lo + hi);
    if (inside(x + mid * dx, y + mid * dy)) hi = mid;
    else lo = mid;
  }
  return Math.min(1, Math.max(0.02, 0.5 * (lo + hi)));
}
