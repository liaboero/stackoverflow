/**
 * Solveur Lattice Boltzmann D2Q9 sur GPU (WebGL2), même modèle que lbm.js :
 * collision régularisée récursive ou BGK, Smagorinsky amorti à la paroi,
 * parois de Bouzidi (y compris minces) avec rotation, entrée/sortie de Guo,
 * bandes absorbantes, zone éponge, accélération de démarrage.
 *
 * Stockage par tampon (ping-pong) : A = (f0, f1, f2, f3), B = (f4, f5, f6, f7)
 * en RGBA32F et C = f8 en R32F, soit 36 octets par maille et par tampon.
 * Un pas = un tracé plein écran avec trois cibles de rendu. Les textures sont
 * liées une fois pour toutes à des unités réservées et chaque sens du
 * ping-pong a sa copie du programme : un pas coûte 3 à 4 appels WebGL.
 *
 * Efforts : calculés à la demande (pas d'échantillonnage) par une passe sur
 * les seules mailles de paroi, réduite dans un accumulateur, puis relus de
 * façon asynchrone (tampon de pixels + barrière) une à deux images plus tard.
 * Différence avec le CPU : les colonnes d'entrée et de sortie sont
 * extrapolées depuis l'état du pas précédent (décalage d'un pas).
 */
import { Geometry } from './geometry.js';

const LIST_W = 256;
const MAX_SAMPLES = 512;
/** Unités de texture réservées au solveur (le rendu utilise 0 à 3). */
const U_SET = [8, 11]; // A, B, C du tampon 0 : 8, 9, 10 ; du tampon 1 : 11, 12, 13
const U_GEO = 14;
const U_Q0 = 15;
const U_Q1 = 16;
const U_LIST = 17;
const U_PART = 18;
const U_NT = [19, 20]; // ν̃ (Spalart–Allmaras), ping-pong synchrone avec les populations
const U_WALL = 21; // distance, normale, zone laminaire
const PBO_COUNT = 3;

const VS = `#version 300 es
void main() {
  vec2 p = vec2(float((gl_VertexID << 1) & 2), float(gl_VertexID & 2));
  gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0);
}`;

const COMMON = `#version 300 es
precision highp float;
precision highp int;
precision highp sampler2D;
uniform sampler2D uA;
uniform sampler2D uB;
uniform sampler2D uC;
uniform sampler2D uGeo;
uniform sampler2D uQ0;
uniform sampler2D uQ1;
uniform ivec2 uSize;
uniform int uBouzidi;
uniform float uSpin;
uniform vec2 uRot;
uniform sampler2D uWall;
uniform int uSA;
uniform float uNu;

const ivec2 CV[9] = ivec2[9](ivec2(0, 0), ivec2(1, 0), ivec2(0, 1), ivec2(-1, 0), ivec2(0, -1),
                             ivec2(1, 1), ivec2(-1, 1), ivec2(-1, -1), ivec2(1, -1));
const float WT[9] = float[9](4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
                             1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0);
const int OPP[9] = int[9](0, 3, 4, 1, 2, 7, 8, 5, 6);
const int REFY[9] = int[9](0, 1, 4, 3, 2, 8, 7, 6, 5);

float getF(ivec2 q, int i) {
  if (i < 4) return texelFetch(uA, q, 0)[i];
  if (i < 8) return texelFetch(uB, q, 0)[i - 4];
  return texelFetch(uC, q, 0).x;
}

int cellType(ivec2 q) {
  return int(texelFetch(uGeo, q, 0).r * 255.0 + 0.5);
}

float linkQ(ivec2 p, int j) {
  if (j < 5) return texelFetch(uQ0, p, 0)[j - 1];
  return texelFetch(uQ1, p, 0)[j - 5];
}

const float KAPPA = 0.41;
const float E_KB = 0.12873490358780418; // e^(−κB), B = 5

// Vitesse de frottement par la loi de Spalding (voir lbm.js).
float spaldingUtau(float ut, float y, float nu) {
  if (ut < 1e-12) return 0.0;
  float rey = y * ut / nu;
  float up = rey < 30.0 ? sqrt(rey) : log(rey) / KAPPA;
  for (int k = 0; k < 20; k++) {
    float ku = KAPPA * up;
    float e = exp(ku);
    float Y = up + E_KB * (e - 1.0 - ku - ku * ku / 2.0 - ku * ku * ku / 6.0);
    float dY = 1.0 + E_KB * KAPPA * (e - 1.0 - ku - ku * ku / 2.0);
    float h = up * Y - rey;
    float nxt = up - h / (Y + up * dY);
    up = nxt > 0.0 ? nxt : 0.5 * up;
    if (abs(h) < 1e-6 * rey) break;
  }
  return ut / up;
}

// Vitesse (populations de l'état lu) d'une maille.
vec2 cellVelocity(ivec2 q) {
  vec4 a = texelFetch(uA, q, 0);
  vec4 b = texelFetch(uB, q, 0);
  float f8 = texelFetch(uC, q, 0).r;
  float rho = a.x + a.y + a.z + a.w + b.x + b.y + b.z + b.w + f8;
  return vec2(a.y - a.w + b.y - b.z - b.w + f8, a.z - b.x + b.y + b.z - b.w - f8) / rho;
}

// Loi de paroi : glissement tangentiel tel que (ν + ν_t)(u_t − u_s)/y = u_τ².
// Renvoie (glissement x, glissement y) ; u_τ dans utau.
vec2 wallSlip(ivec2 p, out float utau) {
  utau = 0.0;
  vec4 w = texelFetch(uWall, p, 0);
  if (uSA == 0) return vec2(0.0);
  vec2 n = w.yz;
  vec2 u = cellVelocity(p);
  if (uSpin != 0.0) u -= uSpin * vec2(-(float(p.y) - uRot.y), float(p.x) - uRot.x);
  vec2 tv = u - dot(u, n) * n;
  float ut = length(tv);
  float yw = max(0.1, w.x);
  // Zone laminaire (w.w = Cf de la couche limite intégrale) : τ_w = ½ Cf u_t².
  bool lam = w.w > 0.0;
  utau = lam ? ut * sqrt(0.5 * w.w) : spaldingUtau(ut, yw, uNu);
  if (ut < 1e-12) return vec2(0.0);
  // Viscosité de paroi κ u_τ y dans les deux cas : la couche limite laminaire,
  // plus mince qu'une maille, n'est pas résolue ; seule sa contrainte compte.
  float us = max(0.0, ut - utau * utau * yw / (uNu + KAPPA * utau * yw));
  return us * tv / ut;
}

// Population j tirée en p à travers une paroi située à la fraction q du lien
// (Bouzidi linéaire, correction de paroi mobile). fi : population incidente.
float bounce(ivec2 p, int j, float q, vec2 slip, out float fi, out vec2 uw, out vec2 wp) {
  ivec2 c = CV[j];
  int i = OPP[j];
  fi = getF(p, i);
  wp = vec2(p) - q * vec2(c);
  uw = uSpin * vec2(-(wp.y - uRot.y), wp.x - uRot.x) + slip;
  float corr = 6.0 * WT[j] * dot(vec2(c), uw);
  if (uBouzidi == 0) return fi + corr;
  if (q < 0.5) {
    ivec2 b = p + c;
    if (b.y >= 0 && b.y < uSize.y && b.x >= 0 && b.x < uSize.x && cellType(b) != 1) {
      return 2.0 * q * fi + (1.0 - 2.0 * q) * getF(b, i) + corr;
    }
    return fi + corr;
  }
  return (fi + (2.0 * q - 1.0) * getF(p, j) + corr) / (2.0 * q);
}

void equilibrium(float rho, vec2 u, out float e[9]) {
  float usq = 1.5 * dot(u, u);
  for (int i = 0; i < 9; i++) {
    float cu = dot(vec2(CV[i]), u);
    e[i] = WT[i] * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - usq);
  }
}
`;

const FS_STEP = `${COMMON}
uniform float uU;
uniform float uTau0;
uniform float uTauSponge;
uniform int uSpongeStart;
uniform float uAbsW;
uniform float uAbsS;
uniform int uLes;
uniform float uC18;
uniform float uAccel;
uniform int uRR;
uniform sampler2D uNt;
layout(location = 0) out vec4 o0;
layout(location = 1) out vec4 o1;
layout(location = 2) out float o2;

void main() {
  ivec2 p = ivec2(gl_FragCoord.xy);
  int x = p.x;
  int y = p.y;
  vec4 geo = texelFetch(uGeo, p, 0);
  int t = int(geo.r * 255.0 + 0.5);
  if (t == 1) {
    // Solide : repos (prêt si la maille redevient fluide).
    o0 = vec4(4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0);
    o1 = vec4(1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0);
    o2 = 1.0 / 36.0;
    return;
  }
  float f[9];

  if (x == 0 || x == uSize.x - 1) {
    // Entrée (vitesse imposée) / sortie (ρ = 1) : extrapolation hors équilibre.
    ivec2 n = ivec2(x == 0 ? 1 : uSize.x - 2, y);
    vec4 a = texelFetch(uA, n, 0);
    vec4 b = texelFetch(uB, n, 0);
    float c = texelFetch(uC, n, 0).r;
    f = float[9](a.x, a.y, a.z, a.w, b.x, b.y, b.z, b.w, c);
    float rho = 0.0;
    for (int i = 0; i < 9; i++) rho += f[i];
    vec2 u = vec2(f[1] - f[3] + f[5] - f[6] - f[7] + f[8], f[2] - f[4] + f[5] + f[6] - f[7] - f[8]) / rho;
    float en[9];
    float eb[9];
    equilibrium(rho, u, en);
    float rb = x == 0 ? rho : 1.0;
    vec2 ub = x == 0 ? vec2(uU, 0.0) : u;
    equilibrium(rb, ub, eb);
    for (int i = 0; i < 9; i++) f[i] = eb[i] + f[i] - en[i];
    o0 = vec4(f[0], f[1], f[2], f[3]);
    o1 = vec4(f[4], f[5], f[6], f[7]);
    o2 = f[8];
    return;
  }

  if (t == 0) {
    f[0] = texelFetch(uA, p, 0).x;
    f[1] = texelFetch(uA, p + ivec2(-1, 0), 0).y;
    f[2] = texelFetch(uA, p + ivec2(0, -1), 0).z;
    f[3] = texelFetch(uA, p + ivec2(1, 0), 0).w;
    f[4] = texelFetch(uB, p + ivec2(0, 1), 0).x;
    f[5] = texelFetch(uB, p + ivec2(-1, -1), 0).y;
    f[6] = texelFetch(uB, p + ivec2(1, -1), 0).z;
    f[7] = texelFetch(uB, p + ivec2(1, 1), 0).w;
    f[8] = texelFetch(uC, p + ivec2(-1, 1), 0).x;
  } else {
    f[0] = texelFetch(uA, p, 0).x;
    float utau;
    vec2 slip = vec2(0.0);
    bool slipDone = false;
    for (int j = 1; j < 9; j++) {
      ivec2 c = CV[j];
      ivec2 s = p - c;
      if (s.y < 0 || s.y >= uSize.y) {
        // Parois haute et basse : glissement (réflexion spéculaire).
        f[j] = getF(ivec2(x - c.x, y), REFY[j]);
        continue;
      }
      float q = linkQ(p, j);
      if (q <= 0.0) {
        f[j] = getF(s, j);
        continue;
      }
      if (!slipDone) {
        slip = wallSlip(p, utau);
        slipDone = true;
      }
      float fi;
      vec2 uw;
      vec2 wp;
      f[j] = bounce(p, j, q, slip, fi, uw, wp);
    }
  }

  float rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
  float ir = 1.0 / rho;
  float ux = (f[1] - f[3] + f[5] - f[6] - f[7] + f[8]) * ir;
  float uy = (f[2] - f[4] + f[5] + f[6] - f[7] - f[8]) * ir;

  // τ local : zone éponge près de la sortie.
  float tau = uTau0;
  if (x > uSpongeStart) {
    float sp = float(x - uSpongeStart) / float(max(1, uSize.x - 1 - uSpongeStart));
    tau += (uTauSponge - uTau0) * sp * sp;
  }
  if (uSA == 1) {
    // Viscosité turbulente ν_t = ν̃ f_v1 de Spalart–Allmaras.
    float v = texelFetch(uNt, p, 0).r;
    float chi = v / uNu;
    float c3 = chi * chi * chi;
    tau += 3.0 * v * c3 / (c3 + 357.911);
  }

  float d4 = f[5] + f[6] + f[7] + f[8];
  float rux = rho * ux;
  float ruy = rho * uy;
  float a1xx = f[1] + f[3] + d4 - rho / 3.0 - rux * ux;
  float a1yy = f[2] + f[4] + d4 - rho / 3.0 - ruy * uy;
  float a1xy = f[5] - f[6] + f[7] - f[8] - rux * uy;
  float om = 1.0 / tau;
  if (uLes == 1) {
    float qn = sqrt(2.0 * (a1xx * a1xx + a1yy * a1yy + 2.0 * a1xy * a1xy));
    om = 2.0 / (tau + sqrt(tau * tau + uC18 * geo.g * qn * ir));
  }

  float pp[9];
  if (uRR == 1) {
    // Collision régularisée récursive, Hermite d'ordre 3 (Malaspinas 2015).
    float k = 1.0 - om;
    float bxx = k * a1xx;
    float byy = k * a1yy;
    float bxy = k * a1xy;
    float Pxx = rux * ux + bxx;
    float Pyy = ruy * uy + byy;
    float Pxy = rux * uy + bxy;
    float Qxxy = rux * ux * uy + 2.0 * ux * bxy + uy * bxx;
    float Qxyy = rux * uy * uy + 2.0 * uy * bxy + ux * byy;
    float sP = Pxx + Pyy;
    float ax = 3.0 * rux;
    float ay = 3.0 * ruy;
    float tx = rho + 3.0 * Pxx - 1.5 * Pyy;
    float ty = rho + 3.0 * Pyy - 1.5 * Pxx;
    float td = rho + 3.0 * sP;
    float qs = 9.0 * (Qxxy + Qxyy);
    float qd = 9.0 * (Qxxy - Qxyy);
    float pxy9 = 9.0 * Pxy;
    pp[0] = (4.0 / 9.0) * (rho - 1.5 * sP);
    pp[1] = (1.0 / 9.0) * (tx + ax - 4.5 * Qxyy);
    pp[3] = (1.0 / 9.0) * (tx - ax + 4.5 * Qxyy);
    pp[2] = (1.0 / 9.0) * (ty + ay - 4.5 * Qxxy);
    pp[4] = (1.0 / 9.0) * (ty - ay + 4.5 * Qxxy);
    pp[5] = (1.0 / 36.0) * (td + ax + ay + pxy9 + qs);
    pp[7] = (1.0 / 36.0) * (td - ax - ay + pxy9 - qs);
    pp[6] = (1.0 / 36.0) * (td - ax + ay - pxy9 + qd);
    pp[8] = (1.0 / 36.0) * (td + ax - ay - pxy9 - qd);
  } else {
    float e[9];
    equilibrium(rho, vec2(ux, uy), e);
    for (int i = 0; i < 9; i++) pp[i] = f[i] + om * (e[i] - f[i]);
  }

  if (uAccel != 0.0) {
    float g1 = uAccel / 3.0;
    float g2 = uAccel / 12.0;
    pp[1] += g1;
    pp[3] -= g1;
    pp[5] += g2;
    pp[8] += g2;
    pp[6] -= g2;
    pp[7] -= g2;
  }

  // Bandes absorbantes haute et basse : relaxation vers l'écoulement amont.
  float dw = float(min(y, uSize.y - 1 - y));
  float sb = uAbsW > 0.0 ? 1.0 - dw / uAbsW : 0.0;
  if (sb > 0.0) {
    float sg = uAbsS * sb * sb;
    float ff[9];
    equilibrium(1.0, vec2(uU, 0.0), ff);
    for (int i = 0; i < 9; i++) pp[i] -= sg * (pp[i] - ff[i]);
  }

  o0 = vec4(pp[0], pp[1], pp[2], pp[3]);
  o1 = vec4(pp[4], pp[5], pp[6], pp[7]);
  o2 = pp[8];
}
`;

const FS_INIT = `${COMMON}
uniform float uU0;
uniform float uKick;
uniform vec2 uKickPos;
uniform float uKickR2;
layout(location = 0) out vec4 o0;
layout(location = 1) out vec4 o1;
layout(location = 2) out float o2;
void main() {
  ivec2 p = ivec2(gl_FragCoord.xy);
  float e[9];
  float rho = 1.0;
  vec2 u = vec2(0.0);
  if (cellType(p) != 1) {
    vec2 d = vec2(p) - uKickPos;
    u = vec2(uU0, uKick * exp(-dot(d, d) / uKickR2));
  }
  equilibrium(rho, u, e);
  o0 = vec4(e[0], e[1], e[2], e[3]);
  o1 = vec4(e[4], e[5], e[6], e[7]);
  o2 = e[8];
}
`;

// Contribution de chaque maille de paroi à la force et au moment.
const FS_FORCE = `${COMMON}
uniform sampler2D uList;
uniform int uCount;
uniform vec2 uRef;
out vec4 o;
void main() {
  ivec2 fp = ivec2(gl_FragCoord.xy);
  int k = fp.y * ${LIST_W} + fp.x;
  if (k >= uCount) {
    o = vec4(0.0);
    return;
  }
  ivec2 p = ivec2(texelFetch(uList, fp, 0).xy + 0.5);
  vec3 acc = vec3(0.0);
  float utau;
  vec2 slip = vec2(0.0);
  bool slipDone = false;
  for (int j = 1; j < 9; j++) {
    ivec2 s = p - CV[j];
    if (s.y < 0 || s.y >= uSize.y) continue;
    float q = linkQ(p, j);
    if (q <= 0.0) continue;
    if (!slipDone) {
      slip = wallSlip(p, utau);
      slipDone = true;
    }
    float fi;
    vec2 uw;
    vec2 wp;
    float val = bounce(p, j, q, slip, fi, uw, wp);
    vec2 cc = vec2(CV[j]);
    float sum = fi + val;
    float dif = fi - val;
    vec2 fl = -cc * sum - uw * dif;
    acc.xy += fl;
    acc.z += (wp.x - uRef.x) * fl.y - (wp.y - uRef.y) * fl.x;
  }
  o = vec4(acc, 0.0);
}
`;

// Spalart–Allmaras : un pas explicite de ν̃ (même schéma que lbm.js).
const FS_SA = `${COMMON}
uniform sampler2D uNt;
uniform float uAbsW;
uniform float uAbsS;
out float o;
const float CB1 = 0.1355;
const float CB2 = 0.622;
const float ISIG = 1.5;
const float CW1 = 0.1355 / (0.41 * 0.41) + (1.0 + 0.622) * 1.5;
const float CW2 = 0.3;
const float CW3_6 = 64.0;
const float CV1_3 = 357.911;
const float NT_MAX = 0.15;
float nt(ivec2 q) {
  return texelFetch(uNt, q, 0).r;
}
void main() {
  ivec2 p = ivec2(gl_FragCoord.xy);
  int x = p.x;
  int y = p.y;
  float inf = 3.0 * uNu;
  int t = cellType(p);
  if (t == 1) {
    o = 0.0;
    return;
  }
  if (x == 0) {
    o = inf;
    return;
  }
  if (x == uSize.x - 1) {
    o = nt(p + ivec2(-1, 0));
    return;
  }
  if (y == 0 || y == uSize.y - 1) {
    o = nt(p + ivec2(0, y == 0 ? 1 : -1));
    return;
  }
  vec4 w = texelFetch(uWall, p, 0);
  if (t == 2) {
    bool wall = false;
    for (int j = 1; j < 9; j++) if (linkQ(p, j) > 0.0) wall = true;
    if (wall) {
      float utau;
      wallSlip(p, utau);
      o = min(NT_MAX, KAPPA * utau * max(0.1, w.x));
      return;
    }
  }
  float v = nt(p);
  float vl = nt(p + ivec2(-1, 0));
  float vr = nt(p + ivec2(1, 0));
  float vd = nt(p + ivec2(0, -1));
  float vu = nt(p + ivec2(0, 1));
  vec2 u = cellVelocity(p);
  vec2 ul = cellVelocity(p + ivec2(-1, 0));
  vec2 ur = cellVelocity(p + ivec2(1, 0));
  vec2 ud = cellVelocity(p + ivec2(0, -1));
  vec2 uu = cellVelocity(p + ivec2(0, 1));
  float adv = (u.x > 0.0 ? u.x * (v - vl) : u.x * (vr - v)) + (u.y > 0.0 ? u.y * (v - vd) : u.y * (vu - v));
  float diff = (uNu + 0.5 * (v + vr)) * (vr - v) - (uNu + 0.5 * (v + vl)) * (v - vl)
             + (uNu + 0.5 * (v + vu)) * (vu - v) - (uNu + 0.5 * (v + vd)) * (v - vd);
  float gx = 0.5 * (vr - vl);
  float gy = 0.5 * (vu - vd);
  diff = (diff + CB2 * (gx * gx + gy * gy)) * ISIG;
  float omega = abs(0.5 * (ur.y - ul.y) - 0.5 * (uu.x - ud.x));
  float d = max(0.5, w.x);
  float kd2 = KAPPA * KAPPA * d * d;
  float chi = v / uNu;
  float c3 = chi * chi * chi;
  float fv1 = c3 / (c3 + CV1_3);
  float fv2 = 1.0 - chi / (1.0 + chi * fv1);
  float sT = max(omega + v / kd2 * fv2, 0.3 * omega);
  float prod = w.w > 0.0 ? 0.0 : CB1 * sT * v;
  float r = sT > 1e-12 ? min(v / (sT * kd2), 10.0) : 10.0;
  float g = r + CW2 * (pow(r, 6.0) - r);
  float fw = g * pow((1.0 + CW3_6) / (pow(g, 6.0) + CW3_6), 1.0 / 6.0);
  float nv = (v + prod + diff - adv) / (1.0 + CW1 * fw * v / (d * d));
  float dw = float(min(y, uSize.y - 1 - y));
  float sb = uAbsW > 0.0 ? 1.0 - dw / uAbsW : 0.0;
  if (sb > 0.0) nv -= uAbsS * sb * sb * (nv - inf);
  o = clamp(nv, 0.0, NT_MAX);
}
`;

// Somme des contributions + indicateur de santé (somme de ρ sur 16 mailles).
const FS_REDUCE = `${COMMON}
uniform sampler2D uPart;
uniform int uCount;
out vec4 o;
void main() {
  vec3 s = vec3(0.0);
  for (int k = 0; k < uCount; k++) s += texelFetch(uPart, ivec2(k % ${LIST_W}, k / ${LIST_W}), 0).xyz;
  float h = 0.0;
  for (int k = 0; k < 16; k++) {
    ivec2 q = ivec2((2 * k + 1) * uSize.x / 32, ((k * 7) % 16 * 2 + 1) * uSize.y / 32);
    vec4 a = texelFetch(uA, q, 0);
    vec4 b = texelFetch(uB, q, 0);
    h += a.x + a.y + a.z + a.w + b.x + b.y + b.z + b.w + texelFetch(uC, q, 0).r;
  }
  o = vec4(s, h);
}
`;

function compile(gl, type, src) {
  const sh = gl.createShader(type);
  gl.shaderSource(sh, src);
  gl.compileShader(sh);
  if (!gl.getShaderParameter(sh, gl.COMPILE_STATUS)) {
    const log = gl.getShaderInfoLog(sh);
    gl.deleteShader(sh);
    throw new Error(`Shader : ${log}`);
  }
  return sh;
}

export function createProgram(gl, vsSrc, fsSrc) {
  const p = gl.createProgram();
  const vs = compile(gl, gl.VERTEX_SHADER, vsSrc);
  const fs = compile(gl, gl.FRAGMENT_SHADER, fsSrc);
  gl.attachShader(p, vs);
  gl.attachShader(p, fs);
  gl.linkProgram(p);
  gl.deleteShader(vs);
  gl.deleteShader(fs);
  if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error(`Programme : ${gl.getProgramInfoLog(p)}`);
  const n = gl.getProgramParameter(p, gl.ACTIVE_UNIFORMS);
  const loc = {};
  for (let i = 0; i < n; i++) {
    const info = gl.getActiveUniform(p, i);
    loc[info.name.replace(/\[0\]$/, '')] = gl.getUniformLocation(p, info.name);
  }
  return { p, loc };
}

export function createTexture(gl, w, h, internal, filter = gl.NEAREST) {
  const t = gl.createTexture();
  // Unité 0 : ne jamais écraser les liaisons des unités réservées au solveur.
  gl.activeTexture(gl.TEXTURE0);
  gl.bindTexture(gl.TEXTURE_2D, t);
  gl.texStorage2D(gl.TEXTURE_2D, 1, internal, w, h);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, filter);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, filter);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  return t;
}

export function createFramebuffer(gl, textures) {
  const fb = gl.createFramebuffer();
  gl.bindFramebuffer(gl.FRAMEBUFFER, fb);
  textures.forEach((t, i) => gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0 + i, gl.TEXTURE_2D, t, 0));
  gl.drawBuffers(textures.map((_, i) => gl.COLOR_ATTACHMENT0 + i));
  const status = gl.checkFramebufferStatus(gl.FRAMEBUFFER);
  gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  if (status !== gl.FRAMEBUFFER_COMPLETE) throw new Error(`Framebuffer incomplet (${status})`);
  return fb;
}

/** WebGL2 avec rendu dans des textures flottantes 32 bits ? */
export function gpuSupported(gl) {
  return !!(gl && gl.getExtension('EXT_color_buffer_float'));
}

export class GpuLBM {
  /**
   * @param {WebGL2RenderingContext} gl
   * @param {object} o mêmes options que LBM
   */
  constructor(
    gl,
    {
      nx,
      ny,
      uIn = 0.1,
      nu = 0.02,
      les = 'auto',
      cs = 0.15,
      lesWallA = 4,
      collision = 'rr',
      bouzidi = true,
      spongeFraction = 0.12,
      absorbWidth = Math.round(ny / 8),
      absorbStrength = 0.05,
    },
  ) {
    if (!gpuSupported(gl)) throw new Error('WebGL2 sans rendu flottant');
    this.gl = gl;
    this.isGpu = true;
    this.nx = nx;
    this.ny = ny;
    this.N = nx * ny;
    const maxTex = gl.getParameter(gl.MAX_TEXTURE_SIZE);
    if (nx > maxTex || ny > maxTex) throw new Error(`Grille ${nx}×${ny} trop grande pour ce GPU`);

    this.geo = new Geometry(nx, ny, lesWallA);
    this.type = this.geo.type;
    this.uIn = uIn;
    this.cs = cs;
    this.lesMode = les;
    this.collision = collision;
    this.bouzidi = bouzidi;
    this.spongeStart = Math.round(nx * (1 - spongeFraction));
    this.absorbWidth = absorbWidth;
    this.absorbStrength = absorbStrength;
    this.refX = nx / 4;
    this.refY = (ny - 1) / 2;
    this.spin = 0;
    this.rotX = this.refX;
    this.rotY = this.refY;
    this.fx = 0;
    this.fy = 0;
    this.mz = 0;
    this.health = 16;
    this.time = 0;
    this.accelX = 0;
    this.wallCount = 0;
    /** Nombre de relectures d'efforts en attente côté GPU (charge du GPU). */
    this.gpuLag = 0;

    this.vao = gl.createVertexArray();
    // Une copie du programme de pas par sens du ping-pong.
    this.progStep = [0, 1].map(() => createProgram(gl, VS, FS_STEP));
    this.progInit = createProgram(gl, VS, FS_INIT);
    this.progForce = createProgram(gl, VS, FS_FORCE);
    this.progReduce = createProgram(gl, VS, FS_REDUCE);
    this.progSA = [0, 1].map(() => createProgram(gl, VS, FS_SA));
    // Spalart–Allmaras : ν̃ en ping-pong, géométrie de paroi (d, n_x, n_y, laminaire).
    this.ntTex = [createTexture(gl, nx, ny, gl.R32F), createTexture(gl, nx, ny, gl.R32F)];
    this.ntFbo = this.ntTex.map((t) => createFramebuffer(gl, [t]));
    this.wallTex = createTexture(gl, nx, ny, gl.RGBA32F);
    this.turbulence = 'none';
    this.saActive = false;

    this.tex = [0, 1].map(() => [
      createTexture(gl, nx, ny, gl.RGBA32F),
      createTexture(gl, nx, ny, gl.RGBA32F),
      createTexture(gl, nx, ny, gl.R32F),
    ]);
    this.fbo = this.tex.map((set) => createFramebuffer(gl, set));
    this.cur = 0;
    this.geoTex = createTexture(gl, nx, ny, gl.RGBA8);
    this.qTex = [createTexture(gl, nx, ny, gl.RGBA8), createTexture(gl, nx, ny, gl.RGBA8)];
    this.accTex = createTexture(gl, MAX_SAMPLES, 1, gl.RGBA32F);
    this.accFbo = createFramebuffer(gl, [this.accTex]);
    this.sampleTimes = [];
    this.sampleBuf = new Float32Array(4 * MAX_SAMPLES);
    this.listTex = null;
    this.partTex = null;
    this.partFbo = null;
    this.listRows = 0;
    this.pboFree = [];
    for (let k = 0; k < PBO_COUNT; k++) {
      const b = gl.createBuffer();
      gl.bindBuffer(gl.PIXEL_PACK_BUFFER, b);
      gl.bufferData(gl.PIXEL_PACK_BUFFER, 16 * MAX_SAMPLES, gl.STREAM_READ);
      this.pboFree.push(b);
    }
    gl.bindBuffer(gl.PIXEL_PACK_BUFFER, null);
    this.pending = [];
    this._dirty = true;

    this.setViscosity(nu);
    this._uploadGeometry();
    this._initPrograms();
    this.reset();
  }

  // ---------------------------------------------------------------- paramètres

  setViscosity(nu) {
    this.nu = nu;
    this.tau0 = 3 * nu + 0.5;
    this.tauSponge = Math.max(this.tau0, 0.9);
    this._updateLes();
  }

  setLes(mode) {
    this.lesMode = mode;
    this._updateLes();
  }

  _updateLes() {
    // Avec Spalart–Allmaras, la LES reste active loin de la paroi (voir lbm.js).
    this.lesActive = this.lesMode === 'on' || (this.lesMode === 'auto' && this.tau0 < 0.51);
  }

  setInletVelocity(u) {
    this.uIn = u;
  }

  get solidCount() {
    return this.geo.solidCount;
  }

  get linkCount() {
    return this.geo.linkCount;
  }

  get lesWallA() {
    return this.geo.lesWallA;
  }

  setLesWallLength(A) {
    this.geo.setLesWallLength(A);
    this._uploadGeometry(false);
  }

  /** Textures A, B, C de l'état courant (lecture par le rendu). */
  get state() {
    return this.tex[this.cur];
  }

  // ---------------------------------------------------------------- géométrie

  setGeometry(inside, o = {}) {
    this.refX = o.refX ?? this.refX;
    this.refY = o.refY ?? this.refY;
    this.spin = o.spin ?? 0;
    this.rotX = o.rotX ?? this.refX;
    this.rotY = o.rotY ?? this.refY;
    // Les mailles solides portent déjà l'équilibre au repos (écrit à chaque pas).
    this.geo.build(inside, o.bbox, o.outline);
    this._uploadGeometry(true);
  }

  _uploadGeometry(withLinks = true) {
    const gl = this.gl;
    const { nx, ny, N } = this;
    const { type, lesDamp, wallQ } = this.geo;
    const g = new Uint8Array(4 * N);
    for (let i = 0; i < N; i++) {
      g[4 * i] = type[i];
      g[4 * i + 1] = Math.round(lesDamp[i] * 255);
      g[4 * i + 3] = 255;
    }
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, this.geoTex);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, nx, ny, gl.RGBA, gl.UNSIGNED_BYTE, g);
    if (withLinks) {
      for (let t = 0; t < 2; t++) {
        const q = new Uint8Array(4 * N);
        for (let c = 0; c < 4; c++) {
          const base = (1 + 4 * t + c) * N;
          for (let i = 0; i < N; i++) {
            const v = wallQ[base + i];
            if (v > 0) q[4 * i + c] = Math.max(1, Math.round(v * 255));
          }
        }
        gl.bindTexture(gl.TEXTURE_2D, this.qTex[t]);
        gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, nx, ny, gl.RGBA, gl.UNSIGNED_BYTE, q);
      }

      // Liste des mailles de paroi pour le calcul des efforts.
      const cells = this.geo.wallCells();
      this.wallCount = cells.length;
      const rows = Math.max(1, Math.ceil(cells.length / LIST_W));
      if (rows !== this.listRows) {
        if (this.listTex) {
          gl.deleteTexture(this.listTex);
          gl.deleteTexture(this.partTex);
          gl.deleteFramebuffer(this.partFbo);
        }
        this.listTex = createTexture(gl, LIST_W, rows, gl.RGBA32F);
        this.partTex = createTexture(gl, LIST_W, rows, gl.RGBA32F);
        this.partFbo = createFramebuffer(gl, [this.partTex]);
        this.listRows = rows;
      }
      const list = new Float32Array(4 * LIST_W * rows);
      cells.forEach((idx, k) => {
        list[4 * k] = idx % nx;
        list[4 * k + 1] = Math.floor(idx / nx);
      });
      gl.activeTexture(gl.TEXTURE0);
      gl.bindTexture(gl.TEXTURE_2D, this.listTex);
      gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, LIST_W, rows, gl.RGBA, gl.FLOAT, list);
      this._uploadWall();
    }
    this._bindUnits();
  }

  /** Distance à la paroi, normale et zone laminaire (modèle de turbulence). */
  _uploadWall() {
    const gl = this.gl;
    const { nx, ny, N } = this;
    const { dist, normal, lam } = this.geo;
    const w = new Float32Array(4 * N);
    for (let i = 0; i < N; i++) {
      w[4 * i] = Math.min(dist[i], 1e5);
      w[4 * i + 1] = normal[2 * i];
      w[4 * i + 2] = normal[2 * i + 1];
      w[4 * i + 3] = lam[i];
    }
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, this.wallTex);
    gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, nx, ny, gl.RGBA, gl.FLOAT, w);
  }

  /** Zone laminaire imposée (transition) : à appeler après modification de geo.lam. */
  updateLaminarZone() {
    this._uploadWall();
    this._bindUnits();
  }

  /** Active ou coupe le modèle de Spalart–Allmaras avec loi de paroi. */
  setTurbulence(mode) {
    this.turbulence = mode === 'sa' ? 'sa' : 'none';
    this.saActive = this.turbulence === 'sa';
    this._updateLes();
    if (this.saActive) this._initSA();
  }

  /** ν̃ = 3ν partout (les solides sont remis à zéro au premier pas). */
  _initSA() {
    const gl = this.gl;
    for (const fb of this.ntFbo) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, fb);
      gl.clearBufferfv(gl.COLOR, 0, new Float32Array([3 * this.nu, 0, 0, 0]));
    }
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    this._dirty = true;
  }

  /** Lie toutes les textures du solveur à leurs unités réservées. */
  _bindUnits() {
    const gl = this.gl;
    const bind = (unit, tex) => {
      gl.activeTexture(gl.TEXTURE0 + unit);
      gl.bindTexture(gl.TEXTURE_2D, tex);
    };
    for (let s = 0; s < 2; s++) for (let k = 0; k < 3; k++) bind(U_SET[s] + k, this.tex[s][k]);
    bind(U_GEO, this.geoTex);
    bind(U_Q0, this.qTex[0]);
    bind(U_Q1, this.qTex[1]);
    if (this.ntTex) {
      bind(U_NT[0], this.ntTex[0]);
      bind(U_NT[1], this.ntTex[1]);
      bind(U_WALL, this.wallTex);
    }
    if (this.listTex) {
      bind(U_LIST, this.listTex);
      bind(U_PART, this.partTex);
    }
    gl.activeTexture(gl.TEXTURE0);
  }

  /** Échantillonneurs et uniformes constants, fixés une fois pour toutes. */
  _initPrograms() {
    const gl = this.gl;
    const setSamplers = (P, set) => {
      gl.useProgram(P.p);
      const L = P.loc;
      if (L.uA) gl.uniform1i(L.uA, U_SET[set]);
      if (L.uB) gl.uniform1i(L.uB, U_SET[set] + 1);
      if (L.uC) gl.uniform1i(L.uC, U_SET[set] + 2);
      if (L.uGeo) gl.uniform1i(L.uGeo, U_GEO);
      if (L.uQ0) gl.uniform1i(L.uQ0, U_Q0);
      if (L.uQ1) gl.uniform1i(L.uQ1, U_Q1);
      if (L.uList) gl.uniform1i(L.uList, U_LIST);
      if (L.uPart) gl.uniform1i(L.uPart, U_PART);
      if (L.uNt) gl.uniform1i(L.uNt, U_NT[set]);
      if (L.uWall) gl.uniform1i(L.uWall, U_WALL);
      if (L.uSize) gl.uniform2i(L.uSize, this.nx, this.ny);
      P.cache = {};
    };
    this.progStep.forEach((P, set) => {
      setSamplers(P, set);
      gl.uniform1i(P.loc.uSpongeStart, this.spongeStart);
      gl.uniform1f(P.loc.uAbsW, this.absorbWidth);
      gl.uniform1f(P.loc.uAbsS, this.absorbStrength);
    });
    this.progSA.forEach((P, set) => {
      setSamplers(P, set);
      // Après un pas, les populations sont dans le tampon `set` et ν̃ courant
      // dans l'autre : la passe lit ν̃[1 − set] et écrit ν̃[set].
      gl.uniform1i(P.loc.uNt, U_NT[1 - set]);
      gl.uniform1f(P.loc.uAbsW, this.absorbWidth);
      gl.uniform1f(P.loc.uAbsS, this.absorbStrength);
    });
    setSamplers(this.progInit, 0);
    setSamplers(this.progForce, 0);
    setSamplers(this.progReduce, 0);
  }

  /** Met à jour un uniforme seulement s'il a changé. */
  _uf(P, name, v) {
    if (P.cache[name] !== v) {
      P.cache[name] = v;
      this.gl.uniform1f(P.loc[name], v);
    }
  }

  _ui(P, name, v) {
    if (P.cache[name] !== v) {
      P.cache[name] = v;
      this.gl.uniform1i(P.loc[name], v);
    }
  }

  _wallUniforms(P) {
    this._ui(P, 'uBouzidi', this.bouzidi ? 1 : 0);
    if (P.loc.uSA) this._ui(P, 'uSA', this.saActive ? 1 : 0);
    if (P.loc.uNu) this._uf(P, 'uNu', this.nu);
    this._uf(P, 'uSpin', this.spin);
    if (P.cache.rotX !== this.rotX || P.cache.rotY !== this.rotY) {
      P.cache.rotX = this.rotX;
      P.cache.rotY = this.rotY;
      this.gl.uniform2f(P.loc.uRot, this.rotX, this.rotY);
    }
  }

  /** Les échantillonneurs des passes d'efforts suivent l'état courant. */
  _followState(P) {
    const set = this.cur;
    if (P.cache.set === set) return;
    P.cache.set = set;
    const gl = this.gl;
    gl.uniform1i(P.loc.uA, U_SET[set]);
    gl.uniform1i(P.loc.uB, U_SET[set] + 1);
    gl.uniform1i(P.loc.uC, U_SET[set] + 2);
  }

  // ---------------------------------------------------------------- calcul

  /** À appeler en début d'image : d'autres tracés ont pu changer l'état WebGL. */
  beginFrame() {
    this._dirty = true;
  }

  /** Initialise l'écoulement (même état initial que le solveur CPU). */
  reset(u0 = this.uIn) {
    const gl = this.gl;
    const P = this.progInit;
    gl.bindVertexArray(this.vao);
    gl.useProgram(P.p);
    const s = this.ny / 8;
    gl.uniform1f(P.loc.uU0, u0);
    gl.uniform1f(P.loc.uKick, 0.15 * Math.max(u0, this.uIn, 0.05));
    gl.uniform2f(P.loc.uKickPos, this.refX + 1.5 * s, this.refY + 0.5 * s);
    gl.uniform1f(P.loc.uKickR2, s * s);
    gl.viewport(0, 0, this.nx, this.ny);
    for (const fb of this.fbo) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, fb);
      gl.drawArrays(gl.TRIANGLES, 0, 3);
    }
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    this.time = 0;
    this.fx = this.fy = this.mz = 0;
    this.health = 16;
    if (this.saActive) this._initSA();
    this._drainPending(false);
    this.sampleTimes = [];
    this._dirty = true;
  }

  step() {
    const gl = this.gl;
    const P = this.progStep[this.cur];
    if (this._dirty) {
      gl.bindVertexArray(this.vao);
      gl.viewport(0, 0, this.nx, this.ny);
      this._dirty = false;
    }
    gl.useProgram(P.p);
    this._uf(P, 'uU', this.uIn);
    this._uf(P, 'uTau0', this.tau0);
    this._uf(P, 'uTauSponge', this.tauSponge);
    this._ui(P, 'uLes', this.lesActive ? 1 : 0);
    this._uf(P, 'uC18', 18 * this.cs * this.cs);
    this._uf(P, 'uAccel', this.accelX);
    this._ui(P, 'uRR', this.collision === 'rr' ? 1 : 0);
    this._wallUniforms(P);
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.fbo[1 - this.cur]);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
    this.cur = 1 - this.cur;
    this.time++;
    if (this.saActive) {
      // ν̃ : lit les populations du nouvel état et l'ancien ν̃, écrit ν̃[cur].
      const S = this.progSA[this.cur];
      gl.useProgram(S.p);
      this._wallUniforms(S);
      gl.bindFramebuffer(gl.FRAMEBUFFER, this.ntFbo[this.cur]);
      gl.drawArrays(gl.TRIANGLES, 0, 3);
    }
  }

  /** Programme le calcul des efforts sur l'état courant (relu par readSamples). */
  sampleForces() {
    if (this.sampleTimes.length >= MAX_SAMPLES) return;
    const gl = this.gl;
    gl.bindVertexArray(this.vao);
    // Contributions par maille de paroi
    let P = this.progForce;
    gl.useProgram(P.p);
    this._followState(P);
    this._wallUniforms(P);
    this._ui(P, 'uCount', this.wallCount);
    if (P.cache.refX !== this.refX || P.cache.refY !== this.refY) {
      P.cache.refX = this.refX;
      P.cache.refY = this.refY;
      gl.uniform2f(P.loc.uRef, this.refX, this.refY);
    }
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.partFbo);
    gl.viewport(0, 0, LIST_W, this.listRows);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
    // Réduction dans l'accumulateur
    P = this.progReduce;
    gl.useProgram(P.p);
    this._followState(P);
    this._ui(P, 'uCount', this.wallCount);
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.accFbo);
    gl.viewport(this.sampleTimes.length, 0, 1, 1);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
    this.sampleTimes.push(this.time);
    this._dirty = true;
  }

  /**
   * Renvoie, dans l'ordre chronologique, les efforts dont la relecture est
   * terminée, puis lance la relecture asynchrone des échantillons programmés
   * depuis l'appel précédent. Ne bloque pas, sauf si le GPU a plus de
   * PBO_COUNT images de retard.
   * @returns {Array<{time:number, fx:number, fy:number, mz:number}>}
   */
  readSamples() {
    const gl = this.gl;
    const out = [];
    this._collect(out, false);
    const n = this.sampleTimes.length;
    if (n) {
      if (!this.pboFree.length) this._collect(out, true);
      const pbo = this.pboFree.pop();
      gl.bindFramebuffer(gl.FRAMEBUFFER, this.accFbo);
      gl.bindBuffer(gl.PIXEL_PACK_BUFFER, pbo);
      gl.readPixels(0, 0, n, 1, gl.RGBA, gl.FLOAT, 0);
      gl.bindBuffer(gl.PIXEL_PACK_BUFFER, null);
      gl.bindFramebuffer(gl.FRAMEBUFFER, null);
      const sync = gl.fenceSync(gl.SYNC_GPU_COMMANDS_COMPLETE, 0);
      gl.flush();
      this.pending.push({ pbo, sync, times: this.sampleTimes });
      this.sampleTimes = [];
    }
    this.gpuLag = this.pending.length;
    return out;
  }

  /** Récupère les relectures terminées (toutes si `wait`). */
  _collect(out, wait) {
    const gl = this.gl;
    while (this.pending.length) {
      const r = this.pending[0];
      if (!wait) {
        const st = gl.clientWaitSync(r.sync, 0, 0);
        if (st === gl.TIMEOUT_EXPIRED) break;
      }
      const n = r.times.length;
      gl.bindBuffer(gl.PIXEL_PACK_BUFFER, r.pbo);
      gl.getBufferSubData(gl.PIXEL_PACK_BUFFER, 0, this.sampleBuf, 0, 4 * n);
      gl.bindBuffer(gl.PIXEL_PACK_BUFFER, null);
      gl.deleteSync(r.sync);
      this.pboFree.push(r.pbo);
      this.pending.shift();
      const b = this.sampleBuf;
      for (let k = 0; k < n; k++) {
        out.push({ time: r.times[k], fx: b[4 * k], fy: b[4 * k + 1], mz: b[4 * k + 2] });
      }
      this.fx = b[4 * (n - 1)];
      this.fy = b[4 * (n - 1) + 1];
      this.mz = b[4 * (n - 1) + 2];
      this.health = b[4 * (n - 1) + 3];
    }
  }

  /** Abandonne ou attend les relectures en cours. */
  _drainPending(keep) {
    const out = [];
    if (keep) this._collect(out, true);
    const gl = this.gl;
    for (const r of this.pending) {
      gl.deleteSync(r.sync);
      this.pboFree.push(r.pbo);
    }
    this.pending = [];
    return out;
  }

  /** Attend la fin des calculs en cours et renvoie les derniers efforts (tests). */
  finish() {
    const out = this.readSamples();
    out.push(...this._drainPending(true));
    this.gpuLag = 0;
    return out;
  }

  /** Relit ρ, ux, uy de tout le domaine (tests ; lent). */
  readMacros() {
    const gl = this.gl;
    const { nx, ny, N } = this;
    const bufs = [0, 1, 2].map(() => new Float32Array(4 * N));
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.fbo[this.cur]);
    bufs.forEach((buf, k) => {
      gl.readBuffer(gl.COLOR_ATTACHMENT0 + k);
      gl.readPixels(0, 0, nx, ny, gl.RGBA, gl.FLOAT, buf);
    });
    gl.readBuffer(gl.COLOR_ATTACHMENT0);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    const rho = new Float32Array(N);
    const ux = new Float32Array(N);
    const uy = new Float32Array(N);
    const [A, B, C] = bufs;
    for (let i = 0; i < N; i++) {
      const f = [A[4 * i], A[4 * i + 1], A[4 * i + 2], A[4 * i + 3], B[4 * i], B[4 * i + 1], B[4 * i + 2], B[4 * i + 3], C[4 * i]];
      const r = f.reduce((a, b) => a + b, 0);
      rho[i] = r;
      ux[i] = (f[1] - f[3] + f[5] - f[6] - f[7] + f[8]) / r;
      uy[i] = (f[2] - f[4] + f[5] + f[6] - f[7] - f[8]) / r;
    }
    this._dirty = true;
    return { rho, ux, uy };
  }

  /**
   * ρ, ux, uy d'une zone rectangulaire (lecture synchrone, à réserver aux
   * petites zones : mailles pariétales, colonne amont).
   */
  readRegion(x0, y0, w, h) {
    const gl = this.gl;
    const n = w * h;
    const bufs = [0, 1, 2].map(() => new Float32Array(4 * n));
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.fbo[this.cur]);
    bufs.forEach((buf, k) => {
      gl.readBuffer(gl.COLOR_ATTACHMENT0 + k);
      gl.readPixels(x0, y0, w, h, gl.RGBA, gl.FLOAT, buf);
    });
    gl.readBuffer(gl.COLOR_ATTACHMENT0);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    const rho = new Float32Array(n);
    const ux = new Float32Array(n);
    const uy = new Float32Array(n);
    const [A, B, C] = bufs;
    for (let i = 0; i < n; i++) {
      const a = 4 * i;
      const f1 = A[a + 1], f2 = A[a + 2], f3 = A[a + 3], f4 = B[a];
      const f5 = B[a + 1], f6 = B[a + 2], f7 = B[a + 3], f8 = C[a];
      const r = A[a] + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
      rho[i] = r;
      ux[i] = (f1 - f3 + f5 - f6 - f7 + f8) / r;
      uy[i] = (f2 - f4 + f5 + f6 - f7 - f8) / r;
    }
    this._dirty = true;
    return { rho, ux, uy };
  }

  isFinite() {
    return Number.isFinite(this.fx) && Number.isFinite(this.fy) && this.health > 8 && this.health < 32;
  }

  /** Pas de lecture des champs côté CPU : le rendu GPU lit les textures. */
  computeMacros() {
    return this.isFinite();
  }

  dispose() {
    const gl = this.gl;
    this._drainPending(false);
    for (const b of this.pboFree) gl.deleteBuffer(b);
    for (const set of this.tex) for (const t of set) gl.deleteTexture(t);
    for (const fb of this.fbo) gl.deleteFramebuffer(fb);
    [this.geoTex, ...this.qTex, this.accTex, this.listTex, this.partTex, ...this.ntTex, this.wallTex].forEach((t) => t && gl.deleteTexture(t));
    [this.accFbo, this.partFbo, ...this.ntFbo].forEach((f) => f && gl.deleteFramebuffer(f));
    [...this.progStep, ...this.progSA, this.progInit, this.progForce, this.progReduce].forEach((p) => gl.deleteProgram(p.p));
    gl.deleteVertexArray(this.vao);
  }
}
