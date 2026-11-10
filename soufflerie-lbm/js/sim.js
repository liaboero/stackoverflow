/**
 * Mise en place d'une « soufflerie » : domaine, maquette, conversion en
 * coefficients aérodynamiques. Partagé par l'application et les tests.
 */
import { LBM } from './lbm.js';
import { makeBody } from './shapes.js';

/** Résolutions : nombre de mailles sur la longueur de référence L. */
export const RESOLUTIONS = {
  rapide: { L: 12, label: 'Rapide' },
  standard: { L: 16, label: 'Standard' },
  fine: { L: 24, label: 'Fine' },
  // ≈ 4 fois plus lente que Fine : destinée à un ordinateur.
  tresfine: { L: 32, label: 'Très fine' },
  // Moteur GPU seulement.
  ultra: { L: 48, label: 'Ultra', gpuOnly: true },
  max: { L: 64, label: 'Max', gpuOnly: true },
};

/**
 * Domaine en multiples de L ; maquette à 8 L de l'entrée, à mi-hauteur.
 * (Cylindre à Re = 100 : Cx ≈ 1,48 et St ≈ 0,176 dans ce domaine, contre
 * 1,43 / 0,167 dans un domaine 40 L × 24 L : l'écart est dû au confinement.)
 */
export const DOMAIN = { length: 24, height: 12, x0: 8 };

/** Durée de la mise en vitesse, en temps convectifs L/U. */
export const RAMP_T = 5;

/**
 * Longueur d'amortissement pariétal de la LES, en mailles : ~1,6 fois
 * l'épaisseur de couche limite laminaire 5 L/√Re, bornée à [1, 4]. Plus
 * longue, la LES serait trop faible pour stabiliser les couches cisaillées
 * décollées aux Reynolds élevés sur grille grossière.
 */
export function lesWallLength(L, re) {
  return Math.min(4, Math.max(1, (8 * L) / Math.sqrt(re)));
}

/**
 * Reynolds au-delà duquel la traînée d'un profil est surestimée de plus de
 * ~10 % faute de résolution. Mesuré sur le NACA 0012 à α = 0 (références
 * Di Ilio et al. 2018 / XFOIL ; voir README) : ≈ 2,5 L² jusqu'à 24 mailles,
 * 5000 à 32 mailles (−4 % à Re = 2000, +1 % à 5000), 10 000 à 48 et 64
 * mailles (collision régularisée : Cx = 0,042 et 0,035 contre 0,037–0,040).
 * Faute de référence, rien n'est garanti au-delà.
 */
const RELIABLE_RE = { 12: 350, 16: 650, 24: 1450, 32: 5000, 48: 10000, 64: 10000 };

export function reliableReynolds(L) {
  return RELIABLE_RE[L] ?? Math.round((2.5 * L * L) / 50) * 50;
}

/** Reynolds d'apparition du sillage tridimensionnel d'un cylindre (mode A). */
export const RE_WAKE_3D = 190;

export class WindTunnel {
  /**
   * @param {object} cfg
   * @param {number} cfg.L mailles par longueur de référence
   * @param {number} cfg.re nombre de Reynolds U L / ν
   * @param {number} [cfg.uIn=0.1] vitesse réseau
   * @param {string} cfg.shape
   * @param {object} [cfg.params] alpha, naca, thickness, spinRatio
   * @param {'auto'|'on'|'off'} [cfg.les]
   * @param {boolean} [cfg.bouzidi]
   */
  constructor(cfg) {
    this.L = cfg.L;
    this.uIn = cfg.uIn ?? 0.1;
    this.nx = DOMAIN.length * this.L;
    this.ny = DOMAIN.height * this.L;
    this.x0 = DOMAIN.x0 * this.L;
    this.y0 = (this.ny - 1) / 2;
    this.re = cfg.re;
    this._acc = { fx: 0, fy: 0, mz: 0, n: 0 };
    const opts = {
      nx: this.nx,
      ny: this.ny,
      uIn: this.uIn,
      nu: this.viscosityFor(cfg.re),
      les: cfg.les ?? 'auto',
      lesWallA: lesWallLength(this.L, cfg.re),
      bouzidi: cfg.bouzidi ?? true,
      collision: cfg.collision ?? 'rr',
    };
    // Moteur : solveur CPU, ou classe GPU fournie avec son contexte WebGL2.
    this.engine = cfg.gpu ? 'gpu' : 'cpu';
    this.lbm = cfg.gpu ? new cfg.gpu.Solver(cfg.gpu.gl, opts) : new LBM(opts);
    this.lbm.setTurbulence?.(cfg.turbulence ?? 'none');
    this.setBody(cfg.shape, cfg.params ?? {});
    // Démarrage progressif : fluide au repos, vitesse d'entrée et accélération
    // uniforme du fluide en rampe (1 − cos) sur RAMP_T temps convectifs. Évite
    // le choc du démarrage impulsif et les ondes acoustiques qui en résultent.
    this.rampSteps = cfg.ramp === false ? 0 : Math.round(RAMP_T * this.stepsPerTStar);
    if (this.rampSteps > 0) {
      this.lbm.setInletVelocity(0);
      this.lbm.reset(0);
    } else {
      this.lbm.reset();
    }
  }

  /** Libère les ressources GPU éventuelles. */
  dispose() {
    this.lbm.dispose?.();
  }

  viscosityFor(re) {
    return (this.uIn * this.L) / re;
  }

  setReynolds(re) {
    this.re = re;
    this.lbm.setViscosity(this.viscosityFor(re));
    this.lbm.setLesWallLength(lesWallLength(this.L, re));
  }

  setInletVelocity(u) {
    this.uIn = u;
    if (!this.ramping) this.lbm.setInletVelocity(u);
    this.lbm.setViscosity(this.viscosityFor(this.re));
    if (this.body) this.setBody(this.shapeId, this.params);
  }

  setBody(shape, params) {
    this.shapeId = shape;
    this.params = { ...params };
    this.body = makeBody(shape, { ...params, L: this.L, x0: this.x0, y0: this.y0, uIn: this.uIn });
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    for (const [x, y] of this.body.outline) {
      x0 = Math.min(x0, x);
      x1 = Math.max(x1, x);
      y0 = Math.min(y0, y);
      y1 = Math.max(y1, y);
    }
    this.lbm.setGeometry(this.body.inside, {
      refX: this.body.refX,
      refY: this.body.refY,
      spin: this.body.spin,
      bbox: { x0, x1, y0, y1 },
      outline: this.body.outline,
    });
  }

  /**
   * Transition laminaire–turbulent imposée au modèle de Spalart–Allmaras à
   * partir de la couche limite intégrale (méthode des panneaux) : en amont de
   * la transition, pas de production turbulente et frottement pariétal pris
   * égal au Cf laminaire calculé.
   * @param {object|null} bl résultat de panelSolve (null : entièrement turbulent)
   */
  setTransition(bl) {
    const geo = this.lbm.geo;
    geo.lam.fill(0);
    if (bl && bl.friction) {
      const body = this.body;
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
      const band = Math.max(3, 0.1 * this.L);
      const cfAt = (side, xc) => {
        const { x, cf } = side;
        let k = 1;
        while (k < x.length - 1 && x[k] < xc) k++;
        // Plafond : pic de Cf au point d'arrêt (abscisse curviligne nulle).
        return Math.min(0.02, Math.max(1e-5, cf[k]));
      };
      const { N } = this.lbm;
      for (let i = 0; i < N; i++) {
        const pos = geo.nearPos[i];
        if (pos < 0 || geo.dist[i] > band) continue;
        const k = Math.min(m - 1, Math.round(pos));
        const [bx, by] = loc[k];
        // Faces : pour le NACA, extrados de l'indice 0 (bord de fuite) au bord
        // d'attaque ; ailleurs (formes symétriques), signe de y dans le repère propre.
        const upper = byIndex ? k <= iLE : by >= 0;
        const xc = (bx - xmin) / span;
        const xtr = upper ? bl.transitionUpper : bl.transitionLower;
        if (Number.isFinite(xtr) && xc < xtr) {
          geo.lam[i] = cfAt(upper ? bl.friction.upper : bl.friction.lower, xc);
        }
      }
    }
    this.lbm.updateLaminarZone?.();
  }

  get ramping() {
    return this.lbm.time < this.rampSteps;
  }

  step(n = 1) {
    const lbm = this.lbm;
    for (let k = 0; k < n; k++) {
      if (lbm.time < this.rampSteps) {
        const T = this.rampSteps;
        const ph = (Math.PI * (lbm.time + 0.5)) / T;
        lbm.uIn = 0.5 * this.uIn * (1 - Math.cos(ph));
        lbm.accelX = ((0.5 * this.uIn * Math.PI) / T) * Math.sin(ph);
      } else if (lbm.accelX !== 0 || lbm.uIn !== this.uIn) {
        lbm.uIn = this.uIn;
        lbm.accelX = 0;
      }
      lbm.step();
      if (!lbm.isGpu) {
        const a = this._acc;
        a.fx += lbm.fx;
        a.fy += lbm.fy;
        a.mz += lbm.mz;
        a.n++;
      }
    }
  }

  /**
   * Coefficients moyennés sur les pas écoulés depuis l'appel précédent (moteur
   * CPU). Un effort lu sur un seul pas peut porter une alternance d'un pas à
   * l'autre (mode pair/impair près d'une paroi mince) qu'un échantillonnage
   * régulier replierait en un biais ; la moyenne l'annule.
   */
  averagedCoefficients() {
    const a = this._acc;
    if (!a.n) return this.coefficients();
    const c = this.coefficientsFrom(a.fx / a.n, a.fy / a.n, a.mz / a.n);
    a.fx = a.fy = a.mz = 0;
    a.n = 0;
    return c;
  }

  /** Temps convectif t* = t U / L. */
  get tStar() {
    return (this.lbm.time * this.uIn) / this.L;
  }

  /** Pas de temps par unité de temps convectif. */
  get stepsPerTStar() {
    return this.L / this.uIn;
  }

  /** Coefficients instantanés (ρ∞ = 1). Cm positif cabreur. */
  coefficients() {
    const { fx, fy, mz } = this.lbm;
    return this.coefficientsFrom(fx, fy, mz);
  }

  coefficientsFrom(fx, fy, mz) {
    const q = 0.5 * this.uIn * this.uIn;
    return {
      cx: fx / (q * this.L),
      cz: fy / (q * this.L),
      cm: -mz / (q * this.L * this.L),
    };
  }

  /** Infos numériques pour l'affichage. */
  info() {
    const lbm = this.lbm;
    return {
      nx: this.nx,
      ny: this.ny,
      L: this.L,
      tau: lbm.tau0,
      mach: this.uIn * Math.sqrt(3),
      blockage: this.body.frontal / this.ny,
      les: lbm.lesActive,
      reLimit: reliableReynolds(this.L),
    };
  }
}
