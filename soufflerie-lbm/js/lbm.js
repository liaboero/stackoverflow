/**
 * Solveur Lattice Boltzmann D2Q9 pour l'écoulement 2D autour d'un corps fixe
 * (ou d'un cylindre en rotation).
 *
 *  - Collision BGK ; viscosité turbulente de Smagorinsky optionnelle (LES).
 *  - Parois : rebond à mi-chemin, interpolé (Bouzidi, Firdaouss & Lallemand 2001)
 *    ou simple (« escalier »), avec correction de paroi mobile.
 *  - Efforts : échange de quantité de mouvement (Ladd 1994), forme
 *    Galiléenne-invariante de Wen et al. (2014) pour les parois mobiles.
 *  - Entrée : vitesse imposée ; sortie : pression imposée (extrapolation hors
 *    équilibre de Guo et al. 2002) précédée d'une zone éponge.
 *  - Parois haute et basse : glissement (réflexion spéculaire).
 *
 * Unités réseau : Δx = Δt = 1, ρ∞ = 1, c_s² = 1/3. L'axe y est orienté vers le
 * haut (ligne y = 0 en bas du domaine).
 *
 * Aucune dépendance au DOM : le module tourne aussi sous Node (tests).
 */

export const CX = [0, 1, 0, -1, 0, 1, -1, -1, 1];
export const CY = [0, 0, 1, 0, -1, 1, 1, -1, -1];
export const W = [4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36];
export const OPP = [0, 3, 4, 1, 2, 7, 8, 5, 6];
/** Direction réfléchie par un plan horizontal (cy → −cy). */
const REFY = [0, 1, 4, 3, 2, 8, 7, 6, 5];

export const FLUID = 0;
export const SOLID = 1;
/** Fluide traité par le chemin lent : voisin d'un solide ou ligne haute/basse. */
export const NEAR_WALL = 2;

/** Équilibre D2Q9 écrit dans `out` (longueur 9). */
export function equilibrium(rho, ux, uy, out) {
  const usq = 1.5 * (ux * ux + uy * uy);
  for (let i = 0; i < 9; i++) {
    const cu = CX[i] * ux + CY[i] * uy;
    out[i] = W[i] * rho * (1 + 3 * cu + 4.5 * cu * cu - usq);
  }
  return out;
}

export class LBM {
  /**
   * @param {object} o
   * @param {number} o.nx  nombre de mailles en x (sens de l'écoulement)
   * @param {number} o.ny  nombre de mailles en y
   * @param {number} [o.uIn=0.1] vitesse d'entrée (unités réseau)
   * @param {number} [o.nu=0.02] viscosité cinématique (unités réseau)
   * @param {'auto'|'on'|'off'} [o.les='auto'] modèle de Smagorinsky
   * @param {number} [o.cs=0.15] constante de Smagorinsky
   * @param {boolean} [o.bouzidi=true] paroi interpolée
   * @param {number} [o.spongeFraction=0.12] longueur de la zone éponge / nx
   */
  constructor({
    nx,
    ny,
    uIn = 0.1,
    nu = 0.02,
    les = 'auto',
    cs = 0.15,
    lesWallA = 4,
    bouzidi = true,
    spongeFraction = 0.12,
    absorbWidth = Math.round(ny / 8),
    absorbStrength = 0.05,
  }) {
    this.nx = nx;
    this.ny = ny;
    const N = (this.N = nx * ny);
    this.f = new Float32Array(9 * N);
    this.g = new Float32Array(9 * N);
    this.type = new Uint8Array(N);
    this.wallQ = new Float32Array(9 * N); // 0 : pas de lien paroi
    this.tauCol = new Float64Array(nx);
    this.omegaCol = new Float64Array(nx);
    // Bandes absorbantes haute et basse : relaxation vers l'écoulement amont,
    // pour évacuer les ondes acoustiques piégées entre les parois.
    this.sigmaRow = new Float64Array(ny);
    // Amortissement pariétal de la viscosité de sous-maille (0 à la paroi → 1 au large).
    this.lesDamp = new Float32Array(N).fill(1);
    this.lesWallA = lesWallA;
    this.rho = new Float32Array(N);
    this.ux = new Float32Array(N);
    this.uy = new Float32Array(N);
    this._tmp = new Float64Array(9);
    this._eqA = new Float64Array(9);
    this._eqB = new Float64Array(9);

    this.uIn = uIn;
    this.cs = cs;
    this.lesMode = les;
    this.bouzidi = bouzidi;
    this.spongeStart = Math.round(nx * (1 - spongeFraction));
    for (let y = 0; y < ny; y++) {
      const d = Math.min(y, ny - 1 - y); // distance à la paroi
      const s = absorbWidth > 0 ? 1 - d / absorbWidth : 0;
      this.sigmaRow[y] = s > 0 ? absorbStrength * s * s : 0;
    }

    // Corps : point de référence des moments, rotation éventuelle (rad/pas).
    this.refX = nx / 4;
    this.refY = (ny - 1) / 2;
    this.spin = 0;
    this.rotX = this.refX;
    this.rotY = this.refY;
    this.solidCount = 0;
    this.linkCount = 0;

    this.fx = 0;
    this.fy = 0;
    this.mz = 0;
    this.time = 0;
    /** Accélération uniforme selon x (démarrage progressif), en unités réseau. */
    this.accelX = 0;

    this.setViscosity(nu);
    this.reset();
  }

  /** Viscosité moléculaire ; met à jour τ par colonne (zone éponge comprise). */
  setViscosity(nu) {
    this.nu = nu;
    this.tau0 = 3 * nu + 0.5;
    const tauSponge = Math.max(this.tau0, 0.9);
    const s0 = this.spongeStart;
    const len = Math.max(1, this.nx - 1 - s0);
    for (let x = 0; x < this.nx; x++) {
      let tau = this.tau0;
      if (x > s0) {
        const s = (x - s0) / len;
        tau += (tauSponge - this.tau0) * s * s;
      }
      this.tauCol[x] = tau;
      this.omegaCol[x] = 1 / tau;
    }
    this._updateLes();
  }

  setLes(mode) {
    this.lesMode = mode;
    this._updateLes();
  }

  _updateLes() {
    this.lesActive = this.lesMode === 'on' || (this.lesMode === 'auto' && this.tau0 < 0.51);
  }

  setInletVelocity(u) {
    this.uIn = u;
  }

  /**
   * Initialise un écoulement uniforme de vitesse u0 (par défaut la vitesse
   * d'entrée), avec une petite perturbation transverse dans le sillage pour
   * déclencher rapidement l'éventuel lâcher tourbillonnaire.
   */
  reset(u0 = this.uIn) {
    const { nx, ny, N, f, type } = this;
    const eq = this._eqA;
    const U = u0;
    const kick = 0.15 * Math.max(u0, this.uIn, 0.05);
    const px = this.refX + 1.5 * (ny / 8);
    const py = this.refY + 0.5 * (ny / 8);
    const r2 = (ny / 8) ** 2;
    for (let y = 0; y < ny; y++) {
      for (let x = 0; x < nx; x++) {
        const idx = y * nx + x;
        if (type[idx] === SOLID) {
          equilibrium(1, 0, 0, eq);
        } else {
          const d2 = (x - px) ** 2 + (y - py) ** 2;
          const vy = kick * Math.exp(-d2 / r2);
          equilibrium(1, U, vy, eq);
        }
        for (let i = 0; i < 9; i++) f[i * N + idx] = eq[i];
      }
    }
    this.g.set(f);
    this.time = 0;
    this.fx = this.fy = this.mz = 0;
  }

  /**
   * Définit le corps solide.
   * @param {(x:number, y:number) => boolean} inside  test d'appartenance (coordonnées réseau,
   *        centres de mailles aux entiers)
   * @param {object} [o]
   * @param {number} [o.refX] point de réduction des moments
   * @param {number} [o.refY]
   * @param {number} [o.spin=0] vitesse de rotation Ωz (rad/pas, >0 sens trigonométrique)
   * @param {number} [o.rotX] centre de rotation
   * @param {number} [o.rotY]
   */
  setGeometry(inside, o = {}) {
    const { nx, ny, N, type, wallQ, f } = this;
    this.refX = o.refX ?? this.refX;
    this.refY = o.refY ?? this.refY;
    this.spin = o.spin ?? 0;
    this.rotX = o.rotX ?? this.refX;
    this.rotY = o.rotY ?? this.refY;

    const eq = equilibrium(1, 0, 0, this._eqA);
    let solids = 0;
    for (let y = 0; y < ny; y++) {
      for (let x = 0; x < nx; x++) {
        const idx = y * nx + x;
        const solid = x >= 2 && x <= nx - 3 && y >= 2 && y <= ny - 3 && inside(x, y);
        if (solid) {
          type[idx] = SOLID;
          solids++;
        } else {
          // Maille qui redevient fluide : on la remet au repos.
          if (type[idx] === SOLID) for (let i = 0; i < 9; i++) f[i * N + idx] = eq[i];
          type[idx] = FLUID;
        }
      }
    }
    this.solidCount = solids;

    // Zone où chercher les parois plus minces qu'une maille (bord de fuite,
    // plaque) : boîte englobante du corps, élargie de 2 mailles.
    const bb = o.bbox;
    const bx0 = bb ? Math.floor(bb.x0) - 2 : 0;
    const bx1 = bb ? Math.ceil(bb.x1) + 2 : -1;
    const by0 = bb ? Math.floor(bb.y0) - 2 : 0;
    const by1 = bb ? Math.ceil(bb.y1) + 2 : -1;

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
          // Lien du nœud fluide vers la paroi : direction i = OPP[j].
          let q = 0;
          if (type[sy * nx + sx] === SOLID) {
            q = linkFraction(inside, x, y, -CX[j], -CY[j]);
          } else if (inBox) {
            // Lien fluide–fluide qui traverse une paroi mince : sans ce
            // traitement, l'écoulement passerait à travers.
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
    this.linkCount = links;
    this.thinLinkCount = thin;
    this._computeWallDistance();
    this._updateLesDamping();
  }

  /** Longueur d'amortissement pariétal de la LES (mailles). */
  setLesWallLength(A) {
    this.lesWallA = A;
    this._updateLesDamping();
  }

  /**
   * Facteur d'amortissement de type Van Driest, (1 − e^(−d/A))², où d est la
   * distance (en mailles) à la paroi du corps. Le modèle de Smagorinsky
   * ajoute sinon de la viscosité dans la couche limite laminaire, où le
   * cisaillement est fort mais l'écoulement résolu, et surestime la traînée
   * de frottement.
   */
  _updateLesDamping() {
    const { N, lesDamp } = this;
    const d = this.wallDist;
    const A = this.lesWallA;
    if (!d || !(A > 0) || this.solidCount === 0) {
      lesDamp.fill(1);
      return;
    }
    for (let i = 0; i < N; i++) {
      // Distance à la paroi ≈ distance au centre de la maille solide − ½.
      const w = Math.max(0, d[i] - 0.5);
      const f = 1 - Math.exp(-w / A);
      lesDamp[i] = f * f;
    }
  }

  /** Distance de chanfrein (1, √2) aux mailles solides, en deux passes. */
  _computeWallDistance() {
    const { nx, ny, N, type } = this;
    const d = (this.wallDist ??= new Float32Array(N));
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

  /** Chemin lent : propagation (pull) pour une maille proche d'une paroi. */
  _pullNearWall(idx, x, y, src, out) {
    const { N, nx, ny, wallQ, type } = this;
    const bouzidi = this.bouzidi;
    const spin = this.spin;
    out[0] = src[idx];
    for (let j = 1; j < 9; j++) {
      const cx = CX[j];
      const cy = CY[j];
      const sy = y - cy;
      if (sy < 0 || sy >= ny) {
        // Paroi de glissement : réflexion spéculaire.
        out[j] = src[REFY[j] * N + y * nx + (x - cx)];
        continue;
      }
      const q = wallQ[j * N + idx];
      if (q === 0) {
        out[j] = src[j * N + sy * nx + (x - cx)];
        continue;
      }
      // Lien vers la paroi : i = OPP[j], c_i = −c_j.
      const i = OPP[j];
      const fi = src[i * N + idx];
      // Point d'impact sur la paroi et vitesse de paroi (rotation).
      const wx = x - q * cx;
      const wy = y - q * cy;
      let uwx = 0;
      let uwy = 0;
      if (spin !== 0) {
        uwx = -spin * (wy - this.rotY);
        uwy = spin * (wx - this.rotX);
      }
      const corr = 6 * W[j] * (cx * uwx + cy * uwy);
      let val;
      if (!bouzidi) {
        val = fi + corr;
      } else if (q < 0.5) {
        const bx = x + cx;
        const by = y + cy;
        if (by >= 0 && by < ny && bx >= 0 && bx < nx && type[by * nx + bx] !== SOLID) {
          val = 2 * q * fi + (1 - 2 * q) * src[i * N + by * nx + bx] + corr;
        } else {
          val = fi + corr;
        }
      } else {
        val = (fi + (2 * q - 1) * src[j * N + idx] + corr) / (2 * q);
      }
      out[j] = val;
      // Échange de quantité de mouvement : F = (c_i − u_w) f_i* − (c_j − u_w) f_j.
      const sum = fi + val;
      const dif = fi - val;
      const fxl = -cx * sum - uwx * dif;
      const fyl = -cy * sum - uwy * dif;
      this._fx += fxl;
      this._fy += fyl;
      this._mz += (wx - this.refX) * fyl - (wy - this.refY) * fxl;
    }
  }

  /** Un pas de temps : propagation + collision fusionnées, puis conditions entrée/sortie. */
  step() {
    const nx = this.nx;
    const ny = this.ny;
    const N = this.N;
    const src = this.f;
    const dst = this.g;
    const type = this.type;
    const tauCol = this.tauCol;
    const omegaCol = this.omegaCol;
    const tmp = this._tmp;
    const les = this.lesActive;
    const c18 = 18 * this.cs * this.cs;
    const lesDamp = this.lesDamp;
    const N2 = 2 * N, N3 = 3 * N, N4 = 4 * N, N5 = 5 * N, N6 = 6 * N, N7 = 7 * N, N8 = 8 * N;
    const w0 = 4 / 9, w1 = 1 / 9, w2 = 1 / 36;

    this._fx = 0;
    this._fy = 0;
    this._mz = 0;
    const ff = equilibrium(1, this.uIn, 0, this._eqB);
    const gx = this.accelX;
    const g1 = gx / 3; // 3 w_i c_ix g : axes
    const g2 = gx / 12; // diagonales
    const q0 = ff[0], q1 = ff[1], q2 = ff[2], q3 = ff[3], q4 = ff[4];
    const q5 = ff[5], q6 = ff[6], q7 = ff[7], q8 = ff[8];

    for (let y = 0; y < ny; y++) {
      const row = y * nx;
      const sg = this.sigmaRow[y];
      for (let x = 1; x < nx - 1; x++) {
        const idx = row + x;
        const t = type[idx];
        if (t === SOLID) continue;
        let f0, f1, f2, f3, f4, f5, f6, f7, f8;
        if (t === FLUID) {
          f0 = src[idx];
          f1 = src[N + idx - 1];
          f2 = src[N2 + idx - nx];
          f3 = src[N3 + idx + 1];
          f4 = src[N4 + idx + nx];
          f5 = src[N5 + idx - nx - 1];
          f6 = src[N6 + idx - nx + 1];
          f7 = src[N7 + idx + nx + 1];
          f8 = src[N8 + idx + nx - 1];
        } else {
          this._pullNearWall(idx, x, y, src, tmp);
          f0 = tmp[0]; f1 = tmp[1]; f2 = tmp[2]; f3 = tmp[3]; f4 = tmp[4];
          f5 = tmp[5]; f6 = tmp[6]; f7 = tmp[7]; f8 = tmp[8];
        }

        const rho = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
        const ir = 1 / rho;
        const ux = (f1 - f3 + f5 - f6 - f7 + f8) * ir;
        const uy = (f2 - f4 + f5 + f6 - f7 - f8) * ir;
        const usq = 1.5 * (ux * ux + uy * uy);
        const r0 = w0 * rho, r1 = w1 * rho, r2 = w2 * rho;
        const ux3 = 3 * ux, uy3 = 3 * uy;
        const uxx = 4.5 * ux * ux - usq;
        const uyy = 4.5 * uy * uy - usq;
        const s = ux + uy, s3 = 3 * s, ss = 4.5 * s * s - usq;
        const d = uy - ux, d3 = 3 * d, dd = 4.5 * d * d - usq;
        const e0 = r0 * (1 - usq);
        const e1 = r1 * (1 + ux3 + uxx);
        const e3 = r1 * (1 - ux3 + uxx);
        const e2 = r1 * (1 + uy3 + uyy);
        const e4 = r1 * (1 - uy3 + uyy);
        const e5 = r2 * (1 + s3 + ss);
        const e7 = r2 * (1 - s3 + ss);
        const e6 = r2 * (1 + d3 + dd);
        const e8 = r2 * (1 - d3 + dd);

        let om;
        if (les) {
          const n1 = f1 - e1, n2 = f2 - e2, n3 = f3 - e3, n4 = f4 - e4;
          const n5 = f5 - e5, n6 = f6 - e6, n7 = f7 - e7, n8 = f8 - e8;
          const pd = n5 + n6 + n7 + n8;
          const pxx = n1 + n3 + pd;
          const pyy = n2 + n4 + pd;
          const pxy = n5 - n6 + n7 - n8;
          const qn = Math.sqrt(2 * (pxx * pxx + pyy * pyy + 2 * pxy * pxy));
          const t0 = tauCol[x];
          om = 2 / (t0 + Math.sqrt(t0 * t0 + c18 * lesDamp[idx] * qn * ir));
        } else {
          om = omegaCol[x];
        }

        let p0 = f0 + om * (e0 - f0);
        let p1 = f1 + om * (e1 - f1);
        let p2 = f2 + om * (e2 - f2);
        let p3 = f3 + om * (e3 - f3);
        let p4 = f4 + om * (e4 - f4);
        let p5 = f5 + om * (e5 - f5);
        let p6 = f6 + om * (e6 - f6);
        let p7 = f7 + om * (e7 - f7);
        let p8 = f8 + om * (e8 - f8);
        if (gx !== 0) {
          p1 += g1;
          p3 -= g1;
          p5 += g2;
          p8 += g2;
          p6 -= g2;
          p7 -= g2;
        }
        if (sg > 0) {
          p0 -= sg * (p0 - q0);
          p1 -= sg * (p1 - q1);
          p2 -= sg * (p2 - q2);
          p3 -= sg * (p3 - q3);
          p4 -= sg * (p4 - q4);
          p5 -= sg * (p5 - q5);
          p6 -= sg * (p6 - q6);
          p7 -= sg * (p7 - q7);
          p8 -= sg * (p8 - q8);
        }
        dst[idx] = p0;
        dst[N + idx] = p1;
        dst[N2 + idx] = p2;
        dst[N3 + idx] = p3;
        dst[N4 + idx] = p4;
        dst[N5 + idx] = p5;
        dst[N6 + idx] = p6;
        dst[N7 + idx] = p7;
        dst[N8 + idx] = p8;
      }
    }

    this._openBoundaries(dst);

    this.f = dst;
    this.g = src;
    this.fx = this._fx;
    this.fy = this._fy;
    this.mz = this._mz;
    this.time++;
  }

  /**
   * Entrée (x = 0) : vitesse imposée, densité extrapolée.
   * Sortie (x = nx−1) : densité ρ = 1 imposée, vitesse extrapolée.
   * Partie hors équilibre copiée du voisin (Guo, Zheng & Shi 2002).
   */
  _openBoundaries(dst) {
    const { nx, ny, N } = this;
    const U = this.uIn;
    const eqB = this._eqA;
    const eqN = this._eqB;
    for (let y = 0; y < ny; y++) {
      // Entrée
      let b = y * nx;
      let n = b + 1;
      let m = moments(dst, n, N);
      equilibrium(m.rho, U, 0, eqB);
      equilibrium(m.rho, m.ux, m.uy, eqN);
      for (let i = 0; i < 9; i++) dst[i * N + b] = eqB[i] + dst[i * N + n] - eqN[i];
      // Sortie
      b = y * nx + nx - 1;
      n = b - 1;
      m = moments(dst, n, N);
      equilibrium(1, m.ux, m.uy, eqB);
      equilibrium(m.rho, m.ux, m.uy, eqN);
      for (let i = 0; i < 9; i++) dst[i * N + b] = eqB[i] + dst[i * N + n] - eqN[i];
    }
  }

  /**
   * Calcule ρ, ux, uy pour l'affichage. Renvoie false si le calcul a divergé.
   */
  computeMacros() {
    const { N, f, type, rho, ux, uy } = this;
    let ok = true;
    for (let idx = 0; idx < N; idx++) {
      if (type[idx] === SOLID) {
        rho[idx] = 1;
        ux[idx] = 0;
        uy[idx] = 0;
        continue;
      }
      const f1 = f[N + idx], f2 = f[2 * N + idx], f3 = f[3 * N + idx], f4 = f[4 * N + idx];
      const f5 = f[5 * N + idx], f6 = f[6 * N + idx], f7 = f[7 * N + idx], f8 = f[8 * N + idx];
      const r = f[idx] + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
      if (!(r > 0.3 && r < 3)) ok = false;
      rho[idx] = r;
      ux[idx] = (f1 - f3 + f5 - f6 - f7 + f8) / r;
      uy[idx] = (f2 - f4 + f5 + f6 - f7 - f8) / r;
    }
    return ok;
  }

  /** Vérification rapide de stabilité (échantillon de mailles + efforts). */
  isFinite() {
    if (!Number.isFinite(this.fx) || !Number.isFinite(this.fy)) return false;
    const { N, f } = this;
    const stride = 97;
    for (let idx = 0; idx < N; idx += stride) {
      const v = f[idx] + f[N + idx] + f[3 * N + idx];
      if (!(v > 0.1 && v < 3)) {
        if (this.type[idx] !== SOLID) return false;
      }
    }
    return true;
  }
}

function moments(f, idx, N) {
  const f0 = f[idx], f1 = f[N + idx], f2 = f[2 * N + idx], f3 = f[3 * N + idx], f4 = f[4 * N + idx];
  const f5 = f[5 * N + idx], f6 = f[6 * N + idx], f7 = f[7 * N + idx], f8 = f[8 * N + idx];
  const rho = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
  return {
    rho,
    ux: (f1 - f3 + f5 - f6 - f7 + f8) / rho,
    uy: (f2 - f4 + f5 + f6 - f7 - f8) / rho,
  };
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
