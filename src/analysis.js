/**
 * analysis.js — 纯函数分析工具：温度分箱、Tg 两段式拟合、温度→颜色映射
 * 不依赖 DOM / Canvas，浏览器与 Node 复用。
 */

/** 把 {T, msd, pe} 样本按温度分箱求均值 */
export function binByT(history, width = 0.04, tLo = 0, tHi = 1.6, minCount = 3) {
  const nb = Math.max(1, Math.ceil((tHi - tLo) / width));
  const acc = Array.from({ length: nb }, () => ({ msd: 0, pe: 0, n: 0 }));
  for (const s of history) {
    if (!isFinite(s.msd) || !isFinite(s.pe)) continue;
    let k = Math.floor((s.T - tLo) / width);
    if (k < 0 || k >= nb) continue;
    acc[k].msd += s.msd;
    acc[k].pe += s.pe;
    acc[k].n++;
  }
  const bins = [];
  for (let k = 0; k < nb; k++) {
    if (acc[k].n < minCount) continue;
    bins.push({
      T: tLo + (k + 0.5) * width,
      msd: acc[k].msd / acc[k].n,
      pe: acc[k].pe / acc[k].n,
      n: acc[k].n,
    });
  }
  return bins;
}

/** 最小二乘直线拟合 y = a x + b，返回 {a, b, sse} */
export function linFit(xs, ys) {
  const n = xs.length;
  let sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (let i = 0; i < n; i++) {
    sx += xs[i]; sy += ys[i];
    sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i];
  }
  const den = n * sxx - sx * sx;
  const a = den === 0 ? 0 : (n * sxy - sx * sy) / den;
  const b = (sy - a * sx) / n;
  let sse = 0;
  for (let i = 0; i < n; i++) {
    const e = ys[i] - (a * xs[i] + b);
    sse += e * e;
  }
  return { a, b, sse };
}

/**
 * 两段式折线拟合找拐点（Tg 的玩具级估计）：
 * 在 log10(MSD)–T 上枚举分割点做两段最小二乘——MSD 跨两个数量级，
 * 对数空间才能把“玻璃侧陡、液体侧缓”的弯折变成干净的两条直线，
 * Tg 取两直线交点。拐点不显著时返回 null。
 */
export function twoSegmentFit(bins) {
  const v = bins.filter((b) => isFinite(b.msd) && b.msd > 0);
  if (v.length < 8) return null;
  const span = v[v.length - 1].T - v[0].T;
  if (span < 0.4) return null;

  const xsAll = v.map((b) => b.T);
  const ysAll = v.map((b) => Math.log10(b.msd));
  const single = linFit(xsAll, ysAll);

  let best = null;
  for (let k = 3; k <= v.length - 4; k++) {
    const Lx = xsAll.slice(0, k + 1), Ly = ysAll.slice(0, k + 1);
    const Rx = xsAll.slice(k + 1), Ry = ysAll.slice(k + 1);
    const f1 = linFit(Lx, Ly);
    const f2 = linFit(Rx, Ry);
    const sse = f1.sse + f2.sse;
    if (!best || sse < best.sse) best = { k, f1, f2, sse };
  }
  if (!best) return null;
  if (single.sse <= 0 || (single.sse - best.sse) / single.sse < 0.5) return null;

  const a1 = best.f1.a, b1 = best.f1.b; // 低温段（log 斜率更陡）
  const a2 = best.f2.a, b2 = best.f2.b; // 高温段
  if (a1 <= 0 || a2 <= 0 || a1 < a2 * 1.3) return null;
  let Tg = (b2 - b1) / (a1 - a2);
  const lo = v[2].T, hi = v[v.length - 3].T;
  if (Tg < lo) Tg = lo;
  if (Tg > hi) Tg = hi;
  return {
    Tg, splitT: v[best.k].T,
    slopeGlass: a1, slopeLiquid: a2, // 注意：log10(MSD) 空间的斜率
    // 画拟合线用：MSD = 10^(a·T + b)
    glass: { a: a1, b: b1, toMsd: (T) => Math.pow(10, a1 * T + b1) },
    liquid: { a: a2, b: b2, toMsd: (T) => Math.pow(10, a2 * T + b2) },
  };
}

// 热成像色标：深钢蓝(冻结) → 冰白 → 琥珀(活跃)。整页只讲冷热一件事
const THERMAL = [
  [0.055, 0.16, 0.30],   // 深钢蓝
  [0.80, 0.87, 0.93],    // 冰白
  [1.00, 0.55, 0.30],    // 琥珀
];

/** 归一化 t ∈ [0,1] → 热成像 RGB（0..1 浮点） */
export function thermal(t) {
  t = Math.min(1, Math.max(0, t));
  const x = t * 2;
  const i = Math.min(1, Math.floor(x));
  const f = x - i;
  const a = THERMAL[i], b = THERMAL[i + 1];
  return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f];
}

/** 256 级查找表（珠子逐帧着色用） */
export const THERMAL_LUT = (() => {
  const lut = new Float32Array(256 * 3);
  for (let i = 0; i < 256; i++) {
    const [r, g, b] = thermal(i / 255);
    lut[i * 3] = r; lut[i * 3 + 1] = g; lut[i * 3 + 2] = b;
  }
  return lut;
})();

/**
 * 温度 → CSS 颜色，全站统一 viridis 色标。
 * lift: 向白色提亮 0..1（曲线上暗紫色在深底不可读，画线时提 0.45）
 */
export function tColorCss(T, alpha = 1, lift = 0) {
  const x = Math.min(1, Math.max(0, (T - 0.05) / 1.45));
  let [r, g, b] = thermal(x);
  if (lift > 0) {
    r = r * (1 - lift) + lift;
    g = g * (1 - lift) + lift;
    b = b * (1 - lift) + lift;
  }
  const R = Math.round(r * 255), G = Math.round(g * 255), B = Math.round(b * 255);
  return alpha >= 1 ? `rgb(${R},${G},${B})` : `rgba(${R},${G},${B},${alpha})`;
}

/** hsl → rgb 写入 Float32Array（sRGB，0..1）；着色用 */
export function hsl2rgb(h, s, l, out) {
  h = ((h % 360) + 360) % 360 / 360;
  if (s === 0) {
    out[0] = out[1] = out[2] = l;
    return out;
  }
  const q = l < 0.5 ? l * (1 + s) : l + s - l * s;
  const p = 2 * l - q;
  out[0] = hue2rgb(p, q, h + 1 / 3);
  out[1] = hue2rgb(p, q, h);
  out[2] = hue2rgb(p, q, h - 1 / 3);
  return out;
}

function hue2rgb(p, q, t) {
  if (t < 0) t += 1;
  if (t > 1) t -= 1;
  if (t < 1 / 6) return p + (q - p) * 6 * t;
  if (t < 1 / 2) return q;
  if (t < 2 / 3) return p + (q - p) * (2 / 3 - t) * 6;
  return p;
}
