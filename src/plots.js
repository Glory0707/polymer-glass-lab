/**
 * plots.js — 两块 Canvas 2D 图表
 *  1) MSD–τ 双对数曲线（当前温度实时 + 历史温度幽灵曲线 + 扩散参考线）
 *  2) 热历史图：固定滞后窗口 MSD 与每珠势能 vs 温度，两段式拟合标注 Tg
 */
import { tColorCss, thermal } from './analysis.js?v=37';

/** 曲线用：热成像提亮，保证深底可读 */
function curveColor(T, alpha = 1) {
  const x = Math.min(1, Math.max(0, (T - 0.05) / 1.45));
  let [r, g, b] = thermal(x);
  r = r * 0.6 + 0.4; g = g * 0.6 + 0.4; b = b * 0.6 + 0.4;
  const R = Math.round(r * 255), G = Math.round(g * 255), B = Math.round(b * 255);
  return alpha >= 1 ? `rgb(${R},${G},${B})` : `rgba(${R},${G},${B},${alpha})`;
}

function prep(canvas) {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const w = canvas.clientWidth;
  const h = canvas.clientHeight;
  if (!w || !h) return null;
  const W = Math.round(w * dpr), H = Math.round(h * dpr);
  if (canvas.width !== W || canvas.height !== H) {
    canvas.width = W; canvas.height = H;
  }
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  /* 透明底：玻璃面板透出 */
  return { ctx, w, h };
}

function grid(ctx, m, w, h, xTicks, yTicks, X, Y) {
  ctx.strokeStyle = 'rgba(255,255,255,0.06)';
  ctx.lineWidth = 1;
  ctx.fillStyle = 'rgba(233,235,242,0.42)';
  ctx.font = '10px "IBM Plex Mono", ui-monospace, Consolas, monospace';
  ctx.textAlign = 'center';
  for (const [v, label] of xTicks) {
    const x = X(v);
    ctx.beginPath();
    ctx.moveTo(x, m.t); ctx.lineTo(x, m.t + (h - m.t - m.b));
    ctx.stroke();
    ctx.fillText(label, x, h - m.b + 14);
  }
  ctx.textAlign = 'right';
  for (const [v, label] of yTicks) {
    const y = Y(v);
    ctx.beginPath();
    ctx.moveTo(m.l, y); ctx.lineTo(w - m.r, y);
    ctx.stroke();
    ctx.fillText(label, m.l - 6, y + 3);
  }
}

function poly(ctx, pts, X, Y) {
  ctx.beginPath();
  let started = false;
  for (const [a, b] of pts) {
    if (!isFinite(a) || !isFinite(b)) continue;
    const x = X(a), y = Y(b);
    if (!started) { ctx.moveTo(x, y); started = true; }
    else ctx.lineTo(x, y);
  }
  ctx.stroke();
}

/** MSD–τ 双对数图。ghosts: [{T, pts}], active: {T, pts} | null；pts = [[τ, msd], ...] */
export function drawMSDPlot(canvas, ghosts, active) {
  const g = prep(canvas);
  if (!g) return;
  const { ctx, w, h } = g;
  const m = { l: 48, r: 16, t: 14, b: 32 };
  const pw = w - m.l - m.r, ph = h - m.t - m.b;

  let maxTau = 10, maxMsd = 10;
  for (const c of [...ghosts, active].filter(Boolean)) {
    for (const [tau, msd] of c.pts) {
      if (tau > maxTau) maxTau = tau;
      if (msd > maxMsd) maxMsd = msd;
    }
  }
  const xlo = -1, xhi = Math.log10(maxTau * 1.35);
  const ylo = -2.5, yhi = Math.log10(maxMsd * 1.6);
  const X = (tau) => m.l + (Math.log10(Math.max(tau, 1e-9)) - xlo) / (xhi - xlo) * pw;
  const Y = (msd) => m.t + (1 - (Math.log10(Math.max(msd, 1e-9)) - ylo) / (yhi - ylo)) * ph;

  const xTicks = [], yTicks = [];
  for (let e = Math.ceil(xlo); e <= Math.floor(xhi); e++) {
    if ((xhi - xlo) > 6 && e % 2) continue;
    xTicks.push([Math.pow(10, e), fmtPow(e)]);
  }
  for (let e = Math.ceil(ylo); e <= Math.floor(yhi); e++) {
    if ((yhi - ylo) > 6 && e % 2) continue;
    yTicks.push([Math.pow(10, e), fmtPow(e)]);
  }
  grid(ctx, m, w, h, xTicks, yTicks, X, Y);

  ctx.save();
  ctx.beginPath();
  ctx.rect(m.l, m.t, pw, ph);
  ctx.clip();

  // 扩散参考线：MSD ∝ τ（经活性曲线 τ≈1 处的值）
  if (active && active.pts.length > 2) {
    let anchor = null;
    for (const [tau, msd] of active.pts) {
      if (tau >= 0.8) { anchor = [tau, msd]; break; }
    }
    if (anchor) {
      const A = anchor[1] / anchor[0];
      ctx.strokeStyle = 'rgba(255,255,255,0.22)';
      ctx.lineWidth = 1;
      ctx.setLineDash([4, 4]);
      poly(ctx, [[0.15, A * 0.15], [maxTau * 1.3, A * maxTau * 1.3]], X, Y);
      ctx.setLineDash([]);
      ctx.fillStyle = 'rgba(255,255,255,0.45)';
      ctx.font = '9.5px "IBM Plex Mono", ui-monospace, Consolas, monospace';
      ctx.textAlign = 'left';
      ctx.fillText('斜率 1', X(0.5), Y(A * 0.5) - 7);
    }
  }

  // 幽灵曲线（历史温度），末端标 T 值
  for (const c of ghosts) {
    ctx.strokeStyle = curveColor(c.T, 0.3);
    ctx.lineWidth = 1.2;
    poly(ctx, c.pts, X, Y);
    if (c.pts.length) {
      const last = c.pts[c.pts.length - 1];
      ctx.fillStyle = curveColor(c.T, 0.55);
      ctx.font = '9px "IBM Plex Mono", ui-monospace, Consolas, monospace';
      ctx.textAlign = 'left';
      ctx.fillText(`T=${c.T.toFixed(2)}`, X(last[0]) + 4, Y(last[1]) + 3);
    }
  }

  // 当前曲线
  if (active && active.pts.length > 1) {
    ctx.strokeStyle = curveColor(active.T, 1);
    ctx.shadowColor = curveColor(active.T, 0.8);
    ctx.shadowBlur = 10;
    ctx.lineWidth = 2.2;
    poly(ctx, active.pts, X, Y);
    const last = active.pts[active.pts.length - 1];
    ctx.fillStyle = curveColor(active.T, 1);
    ctx.beginPath();
    ctx.arc(X(last[0]), Y(last[1]), 2.4 + Math.sin(performance.now() / 280) * 0.9 + 0.9, 0, Math.PI * 2);
    ctx.fill();
    ctx.shadowBlur = 0;
  }
  ctx.restore();

  ctx.fillStyle = 'rgba(233,235,242,0.55)';
  ctx.font = '10px "IBM Plex Mono", ui-monospace, Consolas, monospace';
  ctx.textAlign = 'left';
  ctx.fillText('MSD / σ²', m.l + 6, m.t + 2);
  ctx.textAlign = 'right';
  ctx.fillText('τ (LJ 时间)', w - m.r, h - m.b + 14);
  if (active) {
    ctx.fillStyle = curveColor(active.T, 1);
    ctx.beginPath();
    ctx.arc(w - m.r - 54, m.t + 8, 2.5, 0, Math.PI * 2);
    ctx.fill();
    ctx.font = '600 10.5px "IBM Plex Mono", ui-monospace, Consolas, monospace';
    ctx.textAlign = 'right';
    ctx.fillText(`T = ${active.T.toFixed(2)}`, w - m.r - 4, m.t + 12);
  }
}

function fmtPow(e) {
  if (e === 0) return '1';
  const sup = '⁰¹²³⁴⁵⁶⁷⁸⁹';
  const es = String(Math.abs(e)).split('').map((d) => sup[+d]).join('');
  return e > 0 ? `10${es}` : `10⁻${es}`;
}

/**
 * 非高斯参数 α₂–τ：α₂ = 3⟨Δr⁴⟩ / 5⟨Δr²⟩² − 1，高斯动力学为 0
 */
export function drawA2Plot(canvas, pts) {
  const g = prep(canvas);
  if (!g) return;
  const { ctx, w, h } = g;
  const m = { l: 36, r: 10, t: 12, b: 24 };
  const pw = w - m.l - m.r, ph = h - m.t - m.b;
  const font = '10px "IBM Plex Mono", ui-monospace, Consolas, monospace';

  if (!pts || pts.length < 2) {
    ctx.fillStyle = 'rgba(233,235,242,0.35)';
    ctx.font = '11px "IBM Plex Mono", ui-monospace, Consolas, monospace';
    ctx.textAlign = 'center';
    ctx.fillText('采样中…', m.l + pw / 2, m.t + ph / 2);
    return;
  }
  let maxTau = 1, hi = 0.2, lo = 0;
  for (const [t, a] of pts) {
    if (t > maxTau) maxTau = t;
    if (a > hi) hi = a;
    if (a < lo) lo = a;
  }
  hi *= 1.12;
  const lx = Math.log10(maxTau * 1.25);
  const X = (tau) => m.l + (Math.log10(Math.max(tau, 1)) / lx) * pw;
  const Y = (a) => m.t + (1 - (a - lo) / (hi - lo || 1)) * ph;

  ctx.strokeStyle = 'rgba(255,255,255,0.06)';
  ctx.fillStyle = 'rgba(233,235,242,0.42)';
  ctx.font = font;
  ctx.textAlign = 'center';
  for (let d = -1; d <= Math.ceil(lx); d++) {
    const tau = Math.pow(10, d);
    if (tau > maxTau * 1.25) break;
    const x = X(tau);
    ctx.beginPath(); ctx.moveTo(x, m.t); ctx.lineTo(x, m.t + ph); ctx.stroke();
    ctx.fillText(tau >= 1 ? tau.toFixed(0) : tau.toFixed(1), x, h - m.b + 14);
  }
  // 零线（高斯基准）
  ctx.strokeStyle = 'rgba(255,255,255,0.18)';
  ctx.beginPath(); ctx.moveTo(m.l, Y(0)); ctx.lineTo(m.l + pw, Y(0)); ctx.stroke();
  ctx.fillStyle = 'rgba(233,235,242,0.35)';
  ctx.textAlign = 'left';
  ctx.fillText('α₂ = 0', m.l + 6, Y(0) - 5);

  // 原始曲线（弱化）+ 平滑主线（窗口 5），短 τ 涨落大，平滑后趋势可读
  ctx.strokeStyle = 'rgba(242,244,250,0.28)';
  ctx.lineWidth = 1;
  ctx.beginPath();
  pts.forEach(([t, a], idx) => {
    const x = X(t), y = Y(a);
    idx === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
  });
  ctx.stroke();
  const win = 2;
  ctx.strokeStyle = '#f2f4fa';
  ctx.lineWidth = 2;
  ctx.beginPath();
  pts.forEach(([t, a], idx) => {
    let s = 0, n = 0;
    for (let k = Math.max(0, idx - win); k <= Math.min(pts.length - 1, idx + win); k++) { s += pts[k][1]; n++; }
    const x = X(t), y = Y(s / n);
    idx === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
  });
  ctx.stroke();

  ctx.fillStyle = 'rgba(233,235,242,0.55)';
  ctx.textAlign = 'left';
  ctx.fillText('α₂', m.l + 6, m.t + 2);
  ctx.textAlign = 'right';
  ctx.fillText('τ (LJ 时间)', m.l + pw - 4, h - m.b + 14);
}

/**
 * 热历史图。bins: [{T, msd, pe, n}], fit: {Tg} | null,
 * opts: { peLo, peHi, hasData }
 */
export function drawHistoryPlot(canvas, bins, fit, opts = {}) {
  const g = prep(canvas);
  if (!g) return;
  const { ctx, w, h } = g;
  const m = { l: 48, r: 46, t: 16, b: 32 };
  const pw = w - m.l - m.r, ph = h - m.t - m.b;

  const TMIN = 0, TMAX = 1.6;
  const X = (T) => m.l + (T - TMIN) / (TMAX - TMIN) * pw;
  const ylo = -2, yhi = 2; // MSD: 1e-2 .. 1e2
  const Y = (msd) => m.t + (1 - (Math.log10(Math.max(msd, 1e-9)) - ylo) / (yhi - ylo)) * ph;
  const peLo = opts.peLo ?? -3, peHi = opts.peHi ?? 0;
  const Ype = (pe) => m.t + (1 - (pe - peLo) / (peHi - peLo || 1)) * ph;

  // 网格
  ctx.strokeStyle = 'rgba(255,255,255,0.06)';
  ctx.lineWidth = 1;
  ctx.fillStyle = 'rgba(233,235,242,0.42)';
  ctx.font = '10px "IBM Plex Mono", ui-monospace, Consolas, monospace';
  ctx.textAlign = 'center';
  for (let T = 0.4; T <= 1.41; T += 0.4) {
    const x = X(T);
    ctx.beginPath();
    ctx.moveTo(x, m.t); ctx.lineTo(x, m.t + ph);
    ctx.stroke();
    ctx.fillText(T.toFixed(1), x, h - m.b + 14);
  }
  ctx.textAlign = 'right';
  for (let e = ylo; e <= yhi; e++) {
    const y = Y(Math.pow(10, e));
    ctx.beginPath();
    ctx.moveTo(m.l, y); ctx.lineTo(m.l + pw, y);
    ctx.stroke();
    ctx.fillText(fmtPow(e), m.l - 6, y + 3);
  }

  // 势能（右轴，DSC 类比）
  if (opts.hasData) {
    ctx.strokeStyle = 'rgba(232,163,61,0.85)';
    ctx.lineWidth = 1.2;
    ctx.beginPath();
    let started = false;
    for (const b of bins) {
      const x = X(b.T), y = Ype(b.pe);
      if (!started) { ctx.moveTo(x, y); started = true; }
      else ctx.lineTo(x, y);
    }
    ctx.stroke();
  }

  // 原始样本散点
  if (opts.raw) {
    ctx.fillStyle = 'rgba(255,255,255,0.16)';
    for (const s of opts.raw) {
      ctx.fillRect(X(s.T) - 1, Y(s.msd) - 1, 2, 2);
    }
  }

  // 分箱均值线
  if (bins.length > 1) {
    ctx.strokeStyle = '#f2f4fa';
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(X(bins[0].T), Y(bins[0].msd));
    for (const b of bins) ctx.lineTo(X(b.T), Y(b.msd));
    ctx.stroke();
  }

  // 两段式拟合 + Tg 标注
  if (fit) {
    ctx.strokeStyle = 'rgba(255,255,255,0.28)';
    ctx.lineWidth = 1;
    ctx.setLineDash([5, 4]);
    for (const seg of [fit.seg1, fit.seg2]) {
      ctx.beginPath();
      ctx.moveTo(X(seg[0]), Y(seg[1]));
      ctx.lineTo(X(seg[2]), Y(seg[3]));
      ctx.stroke();
    }
    ctx.setLineDash([]);
    ctx.strokeStyle = 'rgba(255,107,94,0.9)';
    ctx.beginPath();
    ctx.moveTo(X(fit.Tg), m.t);
    ctx.lineTo(X(fit.Tg), m.t + ph);
    ctx.stroke();
    ctx.fillStyle = '#ff8577';
    ctx.font = 'bold 12px "IBM Plex Mono", ui-monospace, Consolas, monospace';
    ctx.textAlign = fit.Tg > 1.1 ? 'right' : 'left';
    ctx.fillText(`Tg ≈ ${fit.Tg.toFixed(2)}`, X(fit.Tg) + (fit.Tg > 1.1 ? -6 : 6), m.t + 14);
  } else if (!opts.hasData) {
    ctx.fillStyle = 'rgba(233,235,242,0.42)';
    ctx.font = '12px "IBM Plex Mono", ui-monospace, Consolas, monospace';
    ctx.textAlign = 'center';
    ctx.fillText('按「降温」扫一遍温度，拐点即 Tg', m.l + pw / 2, m.t + ph / 2);
  }

  // 轴标 + 图例
  // 图例：白=MSD@20τ 橙=势能/珠
  ctx.font = '9.5px "IBM Plex Mono", ui-monospace, Consolas, monospace';
  ctx.textAlign = 'left';
  ctx.strokeStyle = '#f2f4fa';
  ctx.beginPath(); ctx.moveTo(m.l + 6, m.t + 6); ctx.lineTo(m.l + 22, m.t + 6); ctx.stroke();
  ctx.fillStyle = 'rgba(233,235,242,0.6)';
  ctx.fillText('MSD@20τ', m.l + 26, m.t + 9);
  ctx.strokeStyle = 'rgba(232,163,61,0.9)';
  ctx.beginPath(); ctx.moveTo(m.l + 86, m.t + 6); ctx.lineTo(m.l + 102, m.t + 6); ctx.stroke();
  ctx.fillStyle = 'rgba(232,163,61,0.9)';
  ctx.fillText('势能/珠', m.l + 106, m.t + 9);
  if (opts.peTicks) {
    ctx.fillStyle = 'rgba(232,163,61,0.65)';
    for (const [v, label] of opts.peTicks) {
      ctx.fillText(label, w - m.r + 8, Ype(v) + 3);
    }
  }

  ctx.fillStyle = 'rgba(233,235,242,0.55)';
  ctx.font = '10px "IBM Plex Mono", ui-monospace, Consolas, monospace';
  ctx.textAlign = 'right';
  ctx.fillText('T', w - m.r, h - m.b + 14);
}

/**
 * VFT 图：log10(τα) 对 1/T。super-Arrhenius 上弯 = 玻璃化本质；
 * 白虚线 VFT 拟合 log τ = A + B/(T−T0)，灰虚线高温段 Arrhenius 对照。
 */
export function drawVFTPlot(canvas, points, fit, arrhenius) {
  const g = prep(canvas);
  if (!g) return;
  const { ctx, w, h } = g;
  const m = { l: 44, r: 12, t: 12, b: 26 };
  const pw = w - m.l - m.r, ph = h - m.t - m.b;
  const font = '10px "IBM Plex Mono", ui-monospace, Consolas, monospace';

  if (!points || points.length < 2) {
    ctx.fillStyle = 'rgba(233,235,242,0.35)';
    ctx.font = '11px "IBM Plex Mono", ui-monospace, Consolas, monospace';
    ctx.textAlign = 'center';
    ctx.fillText('降温扫描中自动积累 τp(T) 数据点', m.l + pw / 2, m.t + ph / 2);
    return;
  }
  const xs = points.map((p) => 1 / p.T);
  const ys = points.map((p) => p.y);
  let xlo = Math.min(...xs), xhi = Math.max(...xs);
  let ylo = Math.min(...ys), yhi = Math.max(...ys);
  const padX = (xhi - xlo) * 0.08 + 0.02, padY = (yhi - ylo) * 0.12 + 0.05;
  xlo -= padX; xhi += padX; ylo -= padY; yhi += padY;
  const X = (x) => m.l + (x - xlo) / (xhi - xlo) * pw;
  const Y = (y) => m.t + (1 - (y - ylo) / (yhi - ylo)) * ph;

  ctx.strokeStyle = 'rgba(255,255,255,0.06)';
  ctx.fillStyle = 'rgba(233,235,242,0.42)';
  ctx.font = font;
  ctx.textAlign = 'center';
  for (const xv of niceTicks(xlo, xhi, 4)) {
    const x = X(xv);
    ctx.beginPath(); ctx.moveTo(x, m.t); ctx.lineTo(x, m.t + ph); ctx.stroke();
    ctx.fillText(xv.toFixed(1), x, h - m.b + 14);
  }
  ctx.textAlign = 'right';
  for (const yv of niceTicks(ylo, yhi, 3)) {
    const y = Y(yv);
    ctx.beginPath(); ctx.moveTo(m.l, y); ctx.lineTo(m.l + pw, y); ctx.stroke();
    ctx.fillText(yv.toFixed(1), m.l - 5, y + 3);
  }

  // Arrhenius 对照线（高温段直线外推）
  if (arrhenius) {
    ctx.strokeStyle = 'rgba(255,255,255,0.22)';
    ctx.lineWidth = 1;
    ctx.setLineDash([3, 4]);
    poly(ctx, [[xlo, arrhenius.a * xlo + arrhenius.b], [xhi, arrhenius.a * xhi + arrhenius.b]], X, Y);
    ctx.setLineDash([]);
  }

  // VFT 拟合曲线（全温区）
  if (fit) {
    ctx.strokeStyle = 'rgba(255,255,255,0.55)';
    ctx.lineWidth = 1.2;
    ctx.setLineDash([6, 4]);
    const pts = [];
    for (let k = 0; k <= 50; k++) {
      const x = xlo + (xhi - xlo) * k / 50;
      const T = 1 / x;
      if (T > fit.T0 + 0.004) pts.push([x, fit.A + fit.B / (T - fit.T0)]);
    }
    poly(ctx, pts, X, Y);
    ctx.setLineDash([]);
  }

  // 数据点（按温度热成像着色）
  for (const p of points) {
    ctx.fillStyle = curveColor(p.T, 0.95);
    ctx.beginPath();
    ctx.arc(X(1 / p.T), Y(p.y), 3, 0, Math.PI * 2);
    ctx.fill();
  }

  ctx.fillStyle = 'rgba(233,235,242,0.55)';
  ctx.font = font;
  ctx.textAlign = 'left';
  ctx.fillText('log₁₀ τp', m.l + 4, m.t + 2);
  ctx.textAlign = 'right';
  ctx.fillText('1/T', m.l + pw - 2, h - m.b + 14);
  if (fit) {
    ctx.fillStyle = 'rgba(255,255,255,0.6)';
    ctx.fillText(`VFT: T₀ = ${fit.T0.toFixed(2)}`, m.l + pw - 4, m.t + 12);
  }
}

function niceTicks(lo, hi, n) {
  const span = hi - lo;
  if (span <= 0) return [lo];
  const step = Math.pow(10, Math.floor(Math.log10(span / n)));
  const err = span / n / step;
  const mult = err >= 7.5 ? 10 : err >= 3.5 ? 5 : err >= 1.5 ? 2 : 1;
  const s = step * mult;
  const out = [];
  for (let v = Math.ceil(lo / s) * s; v <= hi + 1e-9; v += s) out.push(+v.toFixed(6));
  return out;
}

/**
 * 自中间散射函数 Fs(q*, τ)：半对数横轴，两步弛豫（β 平台 → α 衰减）。
 * 1/e 参考线给出 τα 的图解定义。
 */
export function drawFsqPlot(canvas, pts) {
  const g = prep(canvas);
  if (!g) return;
  const { ctx, w, h } = g;
  const m = { l: 40, r: 12, t: 12, b: 26 };
  const pw = w - m.l - m.r, ph = h - m.t - m.b;
  const font = '10px "IBM Plex Mono", ui-monospace, Consolas, monospace';

  if (!pts || pts.length < 2) {
    ctx.fillStyle = 'rgba(233,235,242,0.35)';
    ctx.font = '11px "IBM Plex Mono", ui-monospace, Consolas, monospace';
    ctx.textAlign = 'center';
    ctx.fillText('采样中…', m.l + pw / 2, m.t + ph / 2);
    return;
  }
  const maxTau = Math.max(...pts.map((p) => p[0])) * 1.25;
  const lx = Math.log10(Math.max(maxTau, 10));
  const X = (t) => m.l + Math.log10(Math.max(t, 1)) / lx * pw;
  const Y = (v) => m.t + (1 - Math.min(1, Math.max(-0.05, v))) * ph;

  // 1/e 参考线
  ctx.strokeStyle = 'rgba(255,107,94,0.45)';
  ctx.setLineDash([4, 4]);
  ctx.beginPath(); ctx.moveTo(m.l, Y(1 / Math.E)); ctx.lineTo(m.l + pw, Y(1 / Math.E)); ctx.stroke();
  ctx.setLineDash([]);
  ctx.fillStyle = 'rgba(255,107,94,0.75)';
  ctx.font = font;
  ctx.textAlign = 'left';
  ctx.fillText('1/e（τα 定义）', m.l + 4, Y(1 / Math.E) - 4);

  ctx.strokeStyle = '#f2f4fa';
  ctx.lineWidth = 1.8;
  ctx.beginPath();
  pts.forEach(([t, v], i) => (i === 0 ? ctx.moveTo(X(t), Y(v)) : ctx.lineTo(X(t), Y(v))));
  ctx.stroke();

  ctx.fillStyle = 'rgba(233,235,242,0.42)';
  ctx.textAlign = 'center';
  for (let d = 0; d <= Math.ceil(lx); d++) {
    const x = X(Math.pow(10, d));
    if (x > m.l + pw + 2) break;
    ctx.fillText(d === 0 ? '1' : '10' + (d === 1 ? '' : d), x, h - m.b + 14);
  }
  ctx.fillStyle = 'rgba(233,235,242,0.55)';
  ctx.textAlign = 'left';
  ctx.fillText('Fs', m.l + 4, m.t + 2);
  ctx.textAlign = 'right';
  ctx.fillText('τ (LJ 时间)', m.l + pw - 4, h - m.b + 14);
}

/**
 * 薄膜迁移率剖面：每层平均位移² 沿 z。自由表面层活动性高、中部低——
 * 薄膜 Tg 下降的直接证据。两端色偏琥珀 = 自由表面。
 */
export function drawProfilePlot(canvas, prof) {
  const g = prep(canvas);
  if (!g) return;
  const { ctx, w, h } = g;
  const m = { l: 40, r: 12, t: 12, b: 26 };
  const pw = w - m.l - m.r, ph = h - m.t - m.b;
  const font = '10px "IBM Plex Mono", ui-monospace, Consolas, monospace';

  if (!prof || !prof.length) {
    ctx.fillStyle = 'rgba(233,235,242,0.35)';
    ctx.font = '11px "IBM Plex Mono", ui-monospace, Consolas, monospace';
    ctx.textAlign = 'center';
    ctx.fillText('样品面板切换到「薄膜」后出现', m.l + pw / 2, m.t + ph / 2);
    return;
  }
  const nb = prof.length;
  const maxV = Math.max(...prof) * 1.15 || 1;
  const bw = pw / nb;
  for (let b = 0; b < nb; b++) {
    const v = Math.max(0, prof[b]);
    const bh = Math.sqrt(v / maxV) * ph; // √ 标度
    const surf = Math.abs(b + 0.5 - nb / 2) / (nb / 2);
    const rr = 0.35 + 0.55 * surf * surf, gg = 0.6, bb = 0.95 - 0.5 * surf * surf;
    ctx.fillStyle = `rgba(${(rr * 255) | 0},${(gg * 255) | 0},${(bb * 255) | 0},0.75)`;
    ctx.fillRect(m.l + b * bw + 1, m.t + ph - bh, bw - 2, bh);
  }
  ctx.fillStyle = 'rgba(233,235,242,0.45)';
  ctx.font = font;
  ctx.textAlign = 'left';
  ctx.fillText('⟨Δr²⟩层 · √标度', m.l + 4, m.t + 2);
  ctx.textAlign = 'center';
  ctx.fillText('中面', m.l + pw / 2, h - m.b + 14);
  ctx.fillStyle = 'rgba(232,163,61,0.8)';
  ctx.textAlign = 'left';
  ctx.fillText('自由面', m.l + 2, h - m.b + 14);
  ctx.textAlign = 'right';
  ctx.fillText('自由面', m.l + pw - 2, h - m.b + 14);
}

/**
 * 力学响应图：单轴/循环为 σ–ε；orientPts 存在时叠加键取向 P2（右轴）
 */
export function drawStressPlot(canvas, pts, orientPts) {
  const g = prep(canvas);
  if (!g) return;
  const { ctx, w, h } = g;
  const m = { l: 44, r: 12, t: 12, b: 26 };
  const pw = w - m.l - m.r, ph = h - m.t - m.b;
  const font = '10px "IBM Plex Mono", ui-monospace, Consolas, monospace';

  // 应变从未扫开（无力学模式）→ 引导文案，避免画出退化坐标轴与纯噪声线
  let spanX = 0;
  if (pts && pts.length) {
    spanX = Math.max(...pts.map((p) => Math.abs(p[0] - pts[0][0])));
  }
  if (!pts || pts.length < 2 || spanX < 1e-10) {
    ctx.fillStyle = 'rgba(233,235,242,0.35)';
    ctx.font = '11px "IBM Plex Mono", ui-monospace, Consolas, monospace';
    ctx.textAlign = 'center';
    ctx.fillText('在「样品」面板选力学模式后实时采样', m.l + pw / 2, m.t + ph / 2);
    return;
  }
  let xlo = Infinity, xhi = -Infinity, ylo = Infinity, yhi = -Infinity;
  for (const [x, y] of pts) {
    if (x < xlo) xlo = x; if (x > xhi) xhi = x;
    if (y < ylo) ylo = y; if (y > yhi) yhi = y;
  }
  const X = (x) => m.l + (x - xlo) / (xhi - xlo || 1) * pw;
  const Y = (y) => m.t + (1 - (y - ylo) / (yhi - ylo || 1)) * ph;

  // 纵轴范围刻度（弱化）
  ctx.fillStyle = 'rgba(233,235,242,0.35)';
  ctx.font = font;
  ctx.textAlign = 'right';
  ctx.fillText(yhi.toFixed(2), m.l - 5, m.t + 8);
  ctx.fillText(ylo.toFixed(2), m.l - 5, m.t + ph);

  if (ylo < 0 && yhi > 0) {
    ctx.strokeStyle = 'rgba(255,255,255,0.08)';
    ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(m.l, Y(0)); ctx.lineTo(m.l + pw, Y(0)); ctx.stroke();
  }

  ctx.strokeStyle = '#f2f4fa';
  ctx.lineWidth = 1.8;
  ctx.beginPath();
  pts.forEach(([x, y], i) => (i === 0 ? ctx.moveTo(X(x), Y(y)) : ctx.lineTo(X(x), Y(y))));
  ctx.stroke();

  // 键取向 P2 叠加（琥珀，右轴 0..0.5）——应力光学的双折射对应量
  const hasP2 = orientPts && orientPts.length > 2;
  if (hasP2) {
    const Y2 = (v) => m.t + (1 - v / 0.5) * ph;
    ctx.strokeStyle = 'rgba(232,163,61,0.9)';
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    let started = false;
    for (const [x, v] of orientPts) {
      if (!isFinite(v)) continue;
      const px = X(x), py = Y2(Math.max(0, Math.min(0.5, v)));
      started ? ctx.lineTo(px, py) : ctx.moveTo(px, py);
      started = true;
    }
    ctx.stroke();
    ctx.fillStyle = 'rgba(232,163,61,0.75)';
    ctx.textAlign = 'right';
    ctx.fillText('0.5', m.l + pw + 9, m.t + 8);
    ctx.fillText('0', m.l + pw + 9, m.t + ph);
  }

  ctx.fillStyle = 'rgba(233,235,242,0.55)';
  ctx.font = '10px "IBM Plex Mono", ui-monospace, Consolas, monospace';
  ctx.textAlign = 'left';
  ctx.fillText('偏应力', m.l + 4, m.t + 2);
  if (hasP2) {
    ctx.fillStyle = 'rgba(232,163,61,0.9)';
    ctx.fillText('取向 P2', m.l + 48, m.t + 2);
  }
  ctx.textAlign = 'right';
  ctx.fillStyle = 'rgba(233,235,242,0.55)';
  ctx.fillText('应变 ε', m.l + pw - 2, h - m.b + 16);
  ctx.textAlign = 'left';
  ctx.fillStyle = 'rgba(233,235,242,0.45)';
  ctx.fillText(xlo.toFixed(2), m.l + 2, h - m.b + 16);
}

/**
 * 记忆实验：PE(t) 白线（左轴）+ T(t) 琥珀线（独立标度）
 */
export function drawProtoPlot(canvas, pts) {
  const g = prep(canvas);
  if (!g) return;
  const { ctx, w, h } = g;
  const m = { l: 44, r: 12, t: 12, b: 26 };
  const pw = w - m.l - m.r, ph = h - m.t - m.b;
  const font = '10px "IBM Plex Mono", ui-monospace, Consolas, monospace';

  if (!pts || pts.length < 2) {
    ctx.fillStyle = 'rgba(233,235,242,0.35)';
    ctx.font = '11px "IBM Plex Mono", ui-monospace, Consolas, monospace';
    ctx.textAlign = 'center';
    ctx.fillText('下方输入温度阶梯，点「跑协议」', m.l + pw / 2, m.t + ph / 2);
    return;
  }
  const t0 = pts[0][0], t1 = pts[pts.length - 1][0];
  const X = (t) => m.l + (t - t0) / (t1 - t0 || 1) * pw;
  let peLo = Infinity, peHi = -Infinity, tLo = Infinity, tHi = -Infinity;
  for (const [, pe, T] of pts) {
    if (pe < peLo) peLo = pe; if (pe > peHi) peHi = pe;
    if (T != null) { if (T < tLo) tLo = T; if (T > tHi) tHi = T; }
  }
  const pad = (peHi - peLo) * 0.15 + 0.01;
  const Ype = (pe) => m.t + (1 - (pe - (peLo - pad)) / (peHi - peLo + 2 * pad)) * ph;
  const hasT = tHi > tLo;
  const Yt = (T) => m.t + (1 - (T - tLo) / (tHi - tLo || 1)) * ph;

  ctx.strokeStyle = 'rgba(255,255,255,0.08)';
  ctx.beginPath(); ctx.moveTo(m.l, Ype(0)); ctx.lineTo(m.l + pw, Ype(0)); ctx.stroke();

  // T(t) 琥珀线
  if (hasT) {
    ctx.strokeStyle = 'rgba(232,163,61,0.85)';
    ctx.lineWidth = 1.4;
    ctx.beginPath();
    pts.forEach(([t, , T], i) => (i === 0 ? ctx.moveTo(X(t), Yt(T)) : ctx.lineTo(X(t), Yt(T))));
    ctx.stroke();
  }

  ctx.strokeStyle = '#f2f4fa';
  ctx.lineWidth = 1.6;
  ctx.beginPath();
  pts.forEach(([t, pe], i) => (i === 0 ? ctx.moveTo(X(t), Ype(pe)) : ctx.lineTo(X(t), Ype(pe))));
  ctx.stroke();

  ctx.font = font;
  ctx.textAlign = 'left';
  ctx.fillStyle = 'rgba(233,235,242,0.6)';
  ctx.fillText('势能/珠', m.l + 4, m.t + 10);
  if (hasT) {
    ctx.fillStyle = 'rgba(232,163,61,0.85)';
    ctx.fillText('T', m.l + 48, m.t + 10);
  }
  ctx.textAlign = 'right';
  ctx.fillStyle = 'rgba(233,235,242,0.45)';
  ctx.fillText('τ=' + t0.toFixed(0) + '→' + t1.toFixed(0), m.l + pw - 4, m.t + 10);
}

/**
 * van Hove 自相关函数 G_s(r, t) 直方图 + 高斯参考线
 */
export function drawVHPlot(canvas, vhBins, vhMax, msdEst) {
  const g = prep(canvas);
  if (!g) return;
  const { ctx, w, h } = g;
  const m = { l: 40, r: 12, t: 12, b: 26 };
  const pw = w - m.l - m.r, ph = h - m.t - m.b;
  const font = '10px "IBM Plex Mono", ui-monospace, Consolas, monospace';

  if (!vhBins || !vhBins.length) {
    ctx.fillStyle = 'rgba(233,235,242,0.35)';
    ctx.font = '11px "IBM Plex Mono", ui-monospace, Consolas, monospace';
    ctx.textAlign = 'center';
    ctx.fillText('采样中…', m.l + pw / 2, m.t + ph / 2);
    return;
  }
  const nb = vhBins.length;
  const dr = vhMax / nb;
  const maxCount = Math.max(...vhBins) * 1.15 || 1;
  // √ 纵轴：冷态主峰极高，线性标度下重尾完全不可见——开方后尾部结构可读
  const Y = (v) => m.t + (1 - Math.sqrt(Math.max(v, 0) / maxCount)) * ph;
  const X = (r) => m.l + (r / vhMax) * pw;

  const bw = pw / nb;
  for (let b = 0; b < nb; b++) {
    const v = vhBins[b];
    if (v <= 0) continue;
    const bh = (v / maxCount) * ph;
    ctx.fillStyle = 'rgba(120, 170, 255, 0.4)';
    ctx.fillRect(X(b * dr) + 1, m.t + ph - bh, bw - 1, bh);
  }
  const msdEst2 = msdEst || 1;
  ctx.strokeStyle = 'rgba(255,107,94,0.8)';
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  for (let k = 0; k <= 60; k++) {
    const r = k / 60 * vhMax;
    const exp3 = 3 / (2 * Math.PI * msdEst2);
    const gv = 4 * Math.PI * r * r * Math.pow(exp3, 1.5) * Math.exp(-3 * r * r / (2 * msdEst2));
    const y = Y(Math.min(gv, maxCount * 1.05));
    const x = X(r);
    k === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
  }
  ctx.stroke();

  ctx.fillStyle = 'rgba(233,235,242,0.45)';
  ctx.font = font;
  ctx.textAlign = 'left';
  ctx.fillText('G_s(r) · √标度', m.l + 4, m.t + 2);
  ctx.fillText('r/σ', m.l + pw - 24, h - m.b + 14);
  ctx.textAlign = 'right';
  ctx.fillStyle = 'rgba(255,107,94,0.7)';
  ctx.fillText('高斯', m.l + pw - 4, m.t + 14);
}
