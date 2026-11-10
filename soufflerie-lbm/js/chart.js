/**
 * Petits graphiques canvas : historique Cx/Cz et polaires.
 */

function setup(canvas) {
  const dpr = Math.min(window.devicePixelRatio || 1, 2.5);
  const w = Math.max(1, Math.round(canvas.clientWidth * dpr));
  const h = Math.max(1, Math.round(canvas.clientHeight * dpr));
  if (canvas.width !== w || canvas.height !== h) {
    canvas.width = w;
    canvas.height = h;
  }
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { ctx, w: w / dpr, h: h / dpr };
}

/** Graduations « rondes » couvrant [lo, hi]. */
export function niceTicks(lo, hi, count = 4) {
  if (!(hi > lo)) {
    hi = lo + 1;
  }
  const raw = (hi - lo) / count;
  const mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= raw) ?? 10 * mag;
  const start = Math.floor(lo / step) * step;
  const end = Math.ceil(hi / step) * step;
  const ticks = [];
  for (let v = start; v <= end + step * 1e-6; v += step) ticks.push(Math.abs(v) < step * 1e-9 ? 0 : v);
  return { ticks, step, lo: start, hi: end };
}

function fmt(v, step) {
  const d = Math.max(0, -Math.floor(Math.log10(step) + 1e-9));
  return v.toFixed(Math.min(3, d)).replace('-', '−');
}

/**
 * Historique temporel.
 * @param {HTMLCanvasElement} canvas
 * @param {object} d  { t: t* (Float64Array), series: [{values, color, label}], tMin, tMax, markers: [t*], meanWindow: [t0, t1] }
 * @param {object} theme { ink, muted, grid, font }
 */
export function drawHistory(canvas, d, theme) {
  const { ctx, w, h } = setup(canvas);
  ctx.clearRect(0, 0, w, h);
  const padL = 40;
  const padR = 10;
  const padT = 8;
  const padB = 22;
  const pw = w - padL - padR;
  const ph = h - padT - padB;
  const { t } = d;

  // Échelle robuste : quantiles 1 % – 99 % pour ne pas être écrasée par un pic
  // isolé (démarrage, changement brutal de géométrie).
  const all = [];
  for (const s of d.series) {
    for (let k = 0; k < s.values.length; k++) if (Number.isFinite(s.values[k])) all.push(s.values[k]);
  }
  all.sort((a, b) => a - b);
  let lo = all.length ? all[Math.floor(0.01 * (all.length - 1))] : NaN;
  let hi = all.length ? all[Math.ceil(0.99 * (all.length - 1))] : NaN;
  if (Number.isFinite(lo)) {
    const pad = 0.12 * (hi - lo || 1);
    lo -= pad;
    hi += pad;
  }
  if (!Number.isFinite(lo)) {
    lo = -0.5;
    hi = 1.5;
  }
  if (hi - lo < 0.2) {
    const m = 0.5 * (hi + lo);
    lo = m - 0.1;
    hi = m + 0.1;
  }
  const yt = niceTicks(lo, hi, 4);
  const xt = niceTicks(d.tMin, d.tMax, 5);
  const X = (v) => padL + ((v - d.tMin) / (d.tMax - d.tMin || 1)) * pw;
  const Y = (v) => padT + (1 - (v - yt.lo) / (yt.hi - yt.lo)) * ph;

  ctx.font = `11px ${theme.font}`;
  ctx.textBaseline = 'middle';
  ctx.lineWidth = 1;
  // Grille horizontale
  for (const v of yt.ticks) {
    const y = Math.round(Y(v)) + 0.5;
    ctx.strokeStyle = v === 0 ? theme.muted : theme.grid;
    ctx.globalAlpha = v === 0 ? 0.6 : 1;
    ctx.beginPath();
    ctx.moveTo(padL, y);
    ctx.lineTo(w - padR, y);
    ctx.stroke();
    ctx.globalAlpha = 1;
    ctx.fillStyle = theme.muted;
    ctx.textAlign = 'right';
    ctx.fillText(fmt(v, yt.step), padL - 6, y);
  }
  // Axe du temps
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  for (const v of xt.ticks) {
    if (v < d.tMin - 1e-9 || v > d.tMax + 1e-9) continue;
    const x = X(v);
    if (x > w - padR - 24) continue; // place pour l'unité
    ctx.fillStyle = theme.muted;
    ctx.fillText(fmt(v, xt.step), x, h - padB + 6);
  }
  ctx.textAlign = 'right';
  ctx.fillText('t*', w - padR, h - padB + 6);

  // Fenêtre de moyenne
  if (d.meanWindow && d.meanWindow[1] > d.meanWindow[0]) {
    const x0 = Math.max(padL, X(d.meanWindow[0]));
    const x1 = Math.min(w - padR, X(d.meanWindow[1]));
    ctx.fillStyle = theme.grid;
    ctx.globalAlpha = 0.45;
    ctx.fillRect(x0, padT, x1 - x0, ph);
    ctx.globalAlpha = 1;
  }
  // Changements de configuration
  ctx.strokeStyle = theme.muted;
  ctx.setLineDash([2, 3]);
  for (const m of d.markers || []) {
    if (m < d.tMin || m > d.tMax) continue;
    const x = Math.round(X(m)) + 0.5;
    ctx.beginPath();
    ctx.moveTo(x, padT);
    ctx.lineTo(x, padT + ph);
    ctx.stroke();
  }
  ctx.setLineDash([]);

  ctx.save();
  ctx.beginPath();
  ctx.rect(padL, padT, pw, ph);
  ctx.clip();
  for (const s of d.series) {
    const v = s.values;
    if (!v.length) continue;
    ctx.strokeStyle = s.color;
    ctx.lineWidth = 1.6;
    ctx.lineJoin = 'round';
    ctx.beginPath();
    // Décimation : au plus ~2 points par pixel.
    const stride = Math.max(1, Math.floor(v.length / (pw * 2)));
    let started = false;
    for (let k = 0; k < v.length; k += stride) {
      if (!Number.isFinite(v[k])) continue;
      const x = X(t[k]);
      const y = Y(v[k]);
      if (!started) {
        ctx.moveTo(x, y);
        started = true;
      } else ctx.lineTo(x, y);
    }
    ctx.stroke();
    // Point final
    const k = v.length - 1;
    if (Number.isFinite(v[k])) {
      ctx.fillStyle = s.color;
      ctx.beginPath();
      ctx.arc(X(t[k]), Y(v[k]), 3, 0, Math.PI * 2);
      ctx.fill();
    }
  }
  ctx.restore();
}

/**
 * Polaire : soit Cx et Cz en fonction de la variable balayée (α ou rotation),
 * soit Cz en fonction de Cx (polaire d'Eiffel).
 * @param {Array<{x:number, cx:number, cz:number}>} pts
 * @param {'var'|'eiffel'} mode
 * @param {string} unit suffixe des graduations de la variable (« ° » pour α)
 */
/**
 * @param {Array<{x:number,cx:number,cz:number,sep?:boolean}>|null} [ref] courbe
 *        de référence (méthode des panneaux) tracée en tirets ; les points
 *        décollés (sep) sont omis.
 */
export function drawPolar(canvas, pts, mode, theme, colors, unit = '', ref = null) {
  const { ctx, w, h } = setup(canvas);
  ctx.clearRect(0, 0, w, h);
  const padL = 40;
  const padR = 12;
  const padT = 10;
  const padB = 26;
  const pw = w - padL - padR;
  const ph = h - padT - padB;
  ctx.font = `11px ${theme.font}`;

  const refPts = ref ? ref.filter((p) => !p.sep) : [];
  if (pts.length < 1 && refPts.length < 2) {
    ctx.fillStyle = theme.muted;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillText('Aucun point pour cette configuration', w / 2, h / 2);
    return;
  }
  const sorted = [...pts].sort((a, b) => a.x - b.x);
  let xs;
  let ys;
  let rxs = [];
  let rys = [];
  if (mode === 'eiffel') {
    xs = sorted.map((p) => p.cx);
    ys = [sorted.map((p) => p.cz)];
    rxs = refPts.map((p) => p.cx);
    rys = [refPts.map((p) => p.cz)];
  } else {
    xs = sorted.map((p) => p.x);
    ys = [sorted.map((p) => p.cz), sorted.map((p) => p.cx)];
    rxs = refPts.map((p) => p.x);
    rys = [refPts.map((p) => p.cz), refPts.map((p) => p.cx)];
  }
  // Axes ajustés aux points mesurés ; la référence n'élargit l'échelle que
  // s'il n'y a pas encore de point (elle est sinon tronquée au cadre).
  const flat = pts.length ? ys.flat() : rys.flat();
  const allX = pts.length ? xs : rxs;
  let ylo = Math.min(0, ...flat);
  let yhi = Math.max(...flat);
  let xlo = Math.min(...allX);
  let xhi = Math.max(...allX);
  if (mode === 'eiffel') xlo = Math.min(0, xlo);
  if (xhi - xlo < 1e-6) {
    xlo -= mode === 'eiffel' ? 0.1 : 2;
    xhi += mode === 'eiffel' ? 0.1 : 2;
  }
  if (yhi - ylo < 0.2) yhi = ylo + 0.2;
  const xt = niceTicks(xlo, xhi, 5);
  const yt = niceTicks(ylo, yhi, 4);
  const X = (v) => padL + ((v - xt.lo) / (xt.hi - xt.lo)) * pw;
  const Y = (v) => padT + (1 - (v - yt.lo) / (yt.hi - yt.lo)) * ph;

  ctx.lineWidth = 1;
  ctx.textBaseline = 'middle';
  for (const v of yt.ticks) {
    const y = Math.round(Y(v)) + 0.5;
    ctx.strokeStyle = v === 0 ? theme.muted : theme.grid;
    ctx.beginPath();
    ctx.moveTo(padL, y);
    ctx.lineTo(w - padR, y);
    ctx.stroke();
    ctx.fillStyle = theme.muted;
    ctx.textAlign = 'right';
    ctx.fillText(fmt(v, yt.step), padL - 6, y);
  }
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  for (const v of xt.ticks) {
    const x = Math.round(X(v)) + 0.5;
    ctx.strokeStyle = v === 0 ? theme.muted : theme.grid;
    ctx.beginPath();
    ctx.moveTo(x, padT);
    ctx.lineTo(x, padT + ph);
    ctx.stroke();
    ctx.fillStyle = theme.muted;
    // Dernière graduation au bord droit : alignée à droite pour rester visible.
    const label = fmt(v, xt.step) + (mode === 'var' ? unit : '');
    const lw = ctx.measureText(label).width;
    ctx.fillText(label, Math.min(x, w - 2 - lw / 2), h - padB + 6);
  }

  const seriesColors = mode === 'var' ? [colors.lift, colors.drag] : [colors.lift];
  if (rxs.length > 1) {
    ctx.save();
    ctx.beginPath();
    ctx.rect(padL, padT, pw, ph);
    ctx.clip();
    ctx.setLineDash([5, 4]);
    ctx.lineWidth = 1.3;
    ctx.globalAlpha = 0.75;
    rys.forEach((vals, s) => {
      ctx.strokeStyle = seriesColors[s];
      ctx.beginPath();
      vals.forEach((v, k) => (k ? ctx.lineTo(X(rxs[k]), Y(v)) : ctx.moveTo(X(rxs[k]), Y(v))));
      ctx.stroke();
    });
    ctx.restore();
  }
  ys.forEach((vals, s) => {
    ctx.strokeStyle = seriesColors[s];
    ctx.fillStyle = seriesColors[s];
    ctx.lineWidth = 1.6;
    ctx.beginPath();
    vals.forEach((v, k) => (k ? ctx.lineTo(X(xs[k]), Y(v)) : ctx.moveTo(X(xs[k]), Y(v))));
    ctx.stroke();
    vals.forEach((v, k) => {
      ctx.beginPath();
      ctx.arc(X(xs[k]), Y(v), 3, 0, Math.PI * 2);
      ctx.fill();
    });
  });
}
