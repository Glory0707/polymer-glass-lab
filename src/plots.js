/**
 * plots.js — 两块 Canvas 2D 图表
 *  1) MSD–τ 双对数曲线（当前温度实时 + 历史温度幽灵曲线 + 扩散参考线）
 *  2) 热历史图：固定滞后窗口 MSD 与每珠势能 vs 温度，两段式拟合标注 Tg
 */
import { tColorCss, thermal } from './analysis.js?v=19';

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

  // 幽灵曲线（历史温度）
  for (const c of ghosts) {
    ctx.strokeStyle = curveColor(c.T, 0.3);
    ctx.lineWidth = 1.2;
    poly(ctx, c.pts, X, Y);
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
  ctx.fillText('α₂ = 0', m.l + pw - 8, Y(0) - 5);

  // 曲线
  ctx.strokeStyle = '#f2f4fa';
  ctx.lineWidth = 2;
  ctx.beginPath();
  pts.forEach(([t, a], idx) => {
    const x = X(t), y = Y(a);
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
}
