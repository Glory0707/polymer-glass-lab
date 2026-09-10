/**
 * main.js — 应用主控：模拟循环、UI 接线、MSD 采样、热历史记录与 Tg 拟合
 */
import { KGSim } from './md.js?v=8';
import { GlassRenderer } from './renderer.js?v=8';
import { drawMSDPlot, drawHistoryPlot } from './plots.js?v=8';
import { binByT, twoSegmentFit, linFit, tColorCss, hsl2rgb, VIRIDIS_LUT } from './analysis.js?v=8';

const $ = (id) => document.getElementById(id);
const T_MIN = 0.05, T_MAX = 1.5;

const state = {
  sim: null,
  renderer: null,
  numChains: 60,
  seed: 20260910,
  paused: false,
  mode: 'free',        // free | cool | heat
  rate: 1e-3,          // ε/τ
  speed: 16,           // MD 步/帧
  colorMode: 'mobility',
  showBonds: true,
  // MSD 曲线
  msdPts: [],
  ghosts: [],
  nextSampleStep: 8,
  // 热历史
  history: [],
  lastRecStep: 0,
  bins: [],
  fit: null,
  fitDirty: false,
  peLo: -3, peHi: 0,
  chainColors: null,
  perf: { tau: 0, fps: 0 },
  annealLeft: 0, // >0 时正在分帧做初始化退火
};

/* ---------------- 构建 / 重建模拟 ---------------- */

function rebuild({ newSeed = false, keepT = true } = {}) {
  if (newSeed) state.seed = (Math.random() * 0x7fffffff) | 0;
  const prevT = state.sim ? state.sim.T : 1.0;
  if (state.renderer) {
    state.renderer.dispose();
    state.renderer = null;
  }
  state.sim = new KGSim({
    numChains: state.numChains,
    chainLen: 40,
    seed: state.seed,
    temperature: keepT ? Math.min(T_MAX, Math.max(T_MIN, prevT)) : 1.0,
    annealSteps: 0, // 浏览器端分帧退火，避免卡死页面
  });
  state.annealLeft = 4000;
  $('anneal').hidden = false;
  state.renderer = new GlassRenderer($('viewport'), state.sim);
  const sampleEl = $('sampleInfo');
  if (sampleEl) {
    sampleEl.textContent =
      `${state.sim.numChains} 链 × ${state.sim.chainLen} 珠 · ρ = 1.0`;
  }
  syncSliderToSim();

  // 链分色查找表（金角分布色相，预先线性化）
  const sim = state.sim;
  state.chainColors = new Float32Array(sim.N * 3);
  const rgb = [0, 0, 0];
  for (let c = 0; c < sim.numChains; c++) {
    hsl2rgb((c * 137.508) % 360, 0.62, 0.58, rgb);
    const r = Math.pow(rgb[0], 2.2), g = Math.pow(rgb[1], 2.2), b = Math.pow(rgb[2], 2.2);
    for (let k = 0; k < 40 && (c * 40 + k) < sim.N; k++) {
      const i3 = (c * 40 + k) * 3;
      state.chainColors[i3] = r; state.chainColors[i3 + 1] = g; state.chainColors[i3 + 2] = b;
    }
  }

  state.msdPts = [];
  state.ghosts = [];
  state.nextSampleStep = 8;
  state.history = [];
  state.lastRecStep = 0;
  state.bins = [];
  state.fit = null;
  state.fitDirty = false;
  $('seedVal').textContent = String(state.seed);
}

function syncSliderToSim() {
  const T = state.sim.T;
  $('tempSlider').value = String(T);
  const v = $('vT');
  v.textContent = T.toFixed(2);
  v.style.color = tColorCss(T, 1, 0.35); // 大数字在深底上需要提亮
  const b = $('railBubble');
  if (b) {
    b.style.left = ((T - 0.05) / 1.45 * 100).toFixed(2) + '%';
    b.textContent = T.toFixed(2);
  }
}

/* ---------------- MSD 曲线采样与幽灵存档 ---------------- */

function sampleMSD() {
  const sim = state.sim;
  const st = sim.stepCount - sim.refStep;
  if (st < state.nextSampleStep) return;
  const tau = st * sim.dt;
  state.msdPts.push([tau, sim.msdRef()]);
  state.nextSampleStep = Math.max(state.nextSampleStep + 8, Math.ceil(state.nextSampleStep * 1.12));
  if (state.msdPts.length > 500) state.msdPts.shift();
}

/** 把当前 MSD 曲线存为幽灵线并重置参考点（换温度时） */
function archiveGhost() {
  if (state.msdPts.length > 5) {
    state.ghosts.push({ T: state.sim.refT, pts: state.msdPts });
    while (state.ghosts.length > 3) state.ghosts.shift();
  }
  state.sim.resetRef();
  state.msdPts = [];
  state.nextSampleStep = 8;
}

/** 温度改变（滑块拖动/急冷/熔化/自动升降温）统一入口 */
function setTemperature(T, { archive = true } = {}) {
  const sim = state.sim;
  T = Math.min(T_MAX, Math.max(T_MIN, T));
  // 拖动中累计漂移过大也先归档一次，避免曲线无法分辨新旧温度
  if (archive && Math.abs(sim.T - sim.refT) > 0.15) archiveGhost();
  sim.T = T;
  syncSliderToSim();
}

/* ---------------- 热历史记录与 Tg 拟合 ---------------- */

function recordHistory() {
  const sim = state.sim;
  if (sim.stepCount - state.lastRecStep < 150) return;
  state.lastRecStep = sim.stepCount;
  // 窗口基本满（≥80%，滞后 16–20τ）即可记录
  if (sim.lagAge() < sim.lagSteps * sim.dt * 0.8) return;
  state.history.push({ T: sim.T, msd: sim.msdLag(), pe: sim.pePerBead });
  if (state.history.length > 6000) {
    // 减半防炸：保留偶数位样本
    state.history = state.history.filter((_, i) => i % 2 === 0);
    state.lastRecStep = sim.stepCount;
  }
  state.fitDirty = true;
}

function refreshFit() {
  state.bins = binByT(state.history, 0.05, 0, 1.6, 1);
  state.fit = twoSegmentFit(state.bins);
  // 快速降温等稀疏场景：分箱加宽再试一次
  if (!state.fit && state.history.length >= 30) {
    const wideBins = binByT(state.history, 0.12, 0, 1.6, 1);
    state.fit = twoSegmentFit(wideBins);
  }
  let peLo = Infinity, peHi = -Infinity;
  for (const b of state.bins) {
    if (b.pe < peLo) peLo = b.pe;
    if (b.pe > peHi) peHi = b.pe;
  }
  if (peLo < peHi) {
    const pad = Math.max(0.05, (peHi - peLo) * 0.15);
    state.peLo = peLo - pad;
    state.peHi = peHi + pad;
  }
  const vTg = $('vTg');
  vTg.textContent = state.fit ? '≈ ' + state.fit.Tg.toFixed(2) : '—';
  vTg.style.color = state.fit ? 'var(--amber)' : '';
}

/** 画热历史图的拟合线段端点（拟合在对数空间，端点换算回 MSD 值） */
function fitSegments(fit, bins) {
  if (!fit || bins.length < 2) return null;
  const t0 = bins[0].T, tSplit = fit.splitT, t1 = bins[bins.length - 1].T;
  return {
    seg1: [t0, fit.glass.toMsd(t0), tSplit, fit.glass.toMsd(tSplit)],
    seg2: [tSplit, fit.liquid.toMsd(tSplit), t1, fit.liquid.toMsd(t1)],
    Tg: fit.Tg,
  };
}

/* ---------------- 珠子着色 ---------------- */

const _rgb = [0, 0, 0];
function updateColors() {
  const sim = state.sim;
  const ct = state.renderer.colorTarget;
  if (state.colorMode === 'chain') {
    ct.set(state.chainColors);
    return;
  }
  if (sim.mobAge() < sim.dt) return; // 窗口尚未建立
  const u = sim.upos, s = sim.snapMob, lut = VIRIDIS_LUT;
  for (let i3 = 0; i3 < u.length; i3 += 3) {
    const dx = u[i3] - s[i3], dy = u[i3 + 1] - s[i3 + 1], dz = u[i3 + 2] - s[i3 + 2];
    const m2 = dx * dx + dy * dy + dz * dz;
    // 10^-3 .. 10^0.25 对数映射：viridis 深紫(冻结) → 亮黄(活跃)
    const x = Math.min(1, Math.max(0, (Math.log10(m2 + 1e-9) + 3) / 3.1));
    const idx = (x * 255) | 0;
    ct[i3] = Math.pow(lut[idx * 3], 2.2);
    ct[i3 + 1] = Math.pow(lut[idx * 3 + 1], 2.2);
    ct[i3 + 2] = Math.pow(lut[idx * 3 + 2], 2.2);
  }
}

/* ---------------- 主循环 ---------------- */

let lastFrame = performance.now();
let frameNo = 0;

function frame(now) {
  requestAnimationFrame(frame);
  const dtms = now - lastFrame;
  lastFrame = now;
  if (dtms > 0 && dtms < 500) {
    state.perf.fps = state.perf.fps * 0.92 + (1000 / dtms) * 0.08;
  }
  frameNo++;

  try {
    // 初始化退火（分帧）：dt/4 小步长化解随机游走初始构象的重叠
    if (state.annealLeft > 0) {
      const sim = state.sim;
      const baseDt = sim.dt;
      sim.dt = baseDt / 4;
      const chunk = Math.min(state.annealLeft, 100);
      for (let s = 0; s < chunk; s++) sim.step();
      state.annealLeft -= chunk;
      sim.dt = baseDt;
      $('annealPct').textContent = Math.round((1 - state.annealLeft / 4000) * 100) + '%';
      state.renderer.update({ showBonds: state.showBonds, showBox: $('boxChk').checked });
      if (state.annealLeft === 0) {
        $('anneal').hidden = true;
        sim.resetRef();
        state.nextSampleStep = 8;
      }
      return;
    }

    if (!state.paused) {
      const sim = state.sim;
      const n = state.speed;
      for (let s = 0; s < n; s++) {
        if (state.mode === 'cool') {
          const next = sim.T - state.rate * sim.dt;
          if (next <= T_MIN) { sim.T = T_MIN; setMode('free'); }
          else sim.T = next;
        } else if (state.mode === 'heat') {
          const next = sim.T + state.rate * sim.dt;
          if (next >= T_MAX) { sim.T = T_MAX; setMode('free'); }
          else sim.T = next;
        }
        // 升降温时参考点跟随温度分段归档，MSD 曲线颜色/标注才与当前温度一致
        if (state.mode !== 'free' && Math.abs(sim.T - sim.refT) > 0.15) archiveGhost();
        sim.step();
        sampleMSD();
        recordHistory();
      }
      // 滑块跟随自动降温/升温
      if (state.mode !== 'free') syncSliderToSim();
    }

    updateColors();
    state.renderer.update({ showBonds: state.showBonds, showBox: $('boxChk').checked });

    if (frameNo % 2 === 0) {
      drawMSDPlot($('msdPlot'), state.ghosts, { T: state.sim.refT, pts: state.msdPts });
    }
    if (state.fitDirty && frameNo % 45 === 0) {
      refreshFit();
      state.fitDirty = false;
    }
    if (frameNo % 3 === 0) {
      const segs = fitSegments(state.fit, state.bins);
      drawHistoryPlot($('histPlot'), state.bins, segs, {
        peLo: state.peLo,
        peHi: state.peHi,
        peTicks: peTicks(),
        raw: state.history,
        hasData: state.history.length > 0,
      });
    }
    if (frameNo % 15 === 0) updateStats();
  } catch (err) {
    showFatal(err);
  }
}

function peTicks() {
  if (!state.history.length) return [];
  const { peLo, peHi } = state;
  if (!isFinite(peLo) || !isFinite(peHi) || peHi <= peLo) return [];
  const mid = (peLo + peHi) / 2;
  return [[peLo, peLo.toFixed(1)], [mid, mid.toFixed(1)], [peHi, peHi.toFixed(1)]];
}

function updateStats() {
  const sim = state.sim;
  const vT = $('vT');
  vT.textContent = sim.T.toFixed(2);
  vT.style.color = tColorCss(sim.T);
  $('vTmeas').textContent = `(${sim.keTemp.toFixed(2)})`;
  $('vPE').textContent = sim.pePerBead.toFixed(2);
  const lastMsd = state.msdPts.length ? state.msdPts[state.msdPts.length - 1][1] : NaN;
  $('vMSD').textContent = isFinite(lastMsd) ? lastMsd.toFixed(2) : '—';
  $('vTau').textContent = sim.time.toFixed(0);
  $('vFps').textContent = Math.round(state.perf.fps) + 'fps';
}

function showFatal(err) {
  console.error(err);
  const box = $('fatal');
  $('fatalMsg').textContent = String(err && err.message ? err.message : err);
  box.hidden = false;
}

/* ---------------- UI 接线 ---------------- */

function setMode(mode) {
  state.mode = mode;
  $('btnCool').classList.toggle('active', mode === 'cool');
  $('btnHeat').classList.toggle('active', mode === 'heat');
}

function bindUI() {
  const tempSlider = $('tempSlider');
  tempSlider.addEventListener('input', () => {
    setMode('free');
    setTemperature(parseFloat(tempSlider.value));
  });
  tempSlider.addEventListener('change', () => archiveGhost());

  $('btnCool').addEventListener('click', () => setMode('cool'));
  $('btnHeat').addEventListener('click', () => setMode('heat'));

  const rateSlider = $('rateSlider');
  const rateLabel = $('rateVal');
  const showRate = () => {
    rateLabel.textContent = state.rate >= 1e-2
      ? state.rate.toFixed(1) + '×10⁻²'
      : (state.rate * 1e3).toFixed(1) + '×10⁻³';
  };
  rateSlider.addEventListener('input', () => {
    state.rate = Math.pow(10, -4 + (parseInt(rateSlider.value, 10) / 100) * 2.3);
    showRate();
  });
  state.rate = Math.pow(10, -4 + (43 / 100) * 2.3); // ≈1.0×10⁻³
  showRate();

  $('btnQuench').addEventListener('click', () => { setMode('free'); setTemperature(T_MIN); archiveGhost(); });
  $('btnMelt').addEventListener('click', () => { setMode('free'); setTemperature(T_MAX); archiveGhost(); });

  const pauseBtn = $('btnPause');
  pauseBtn.addEventListener('click', () => {
    state.paused = !state.paused;
    pauseBtn.textContent = state.paused ? '继续' : '暂停';
  });

  $('speedSlider').addEventListener('input', (e) => {
    state.speed = parseInt(e.target.value, 10);
    $('speedVal').textContent = String(state.speed);
  });

  $('sizeSel').addEventListener('change', (e) => {
    state.numChains = parseInt(e.target.value, 10);
    rebuild();
  });
  $('btnNew').addEventListener('click', () => rebuild({ newSeed: true }));
  $('btnResetRef').addEventListener('click', () => archiveGhost());

  $('colorSel').addEventListener('change', (e) => { state.colorMode = e.target.value; });
  $('bondsChk').addEventListener('change', (e) => { state.showBonds = e.target.checked; });

  $('figsToggle').addEventListener('click', () => $('figs').classList.toggle('open'));
  $('figsClose').addEventListener('click', () => $('figs').classList.remove('open'));
  $('paramsToggle').addEventListener('click', () => $('paramsPop').classList.toggle('hidden'));
  $('intro').addEventListener('click', () => $('intro').classList.add('gone'));
  setTimeout(() => $('intro')?.classList.add('gone'), 14000);

  window.addEventListener('keydown', (e) => {
    if (e.target.tagName === 'INPUT' || e.target.tagName === 'SELECT') return;
    if (e.code === 'Space') { e.preventDefault(); pauseBtn.click(); }
    else if (e.key === 'c') $('btnQuench').click();
    else if (e.key === 'h') $('btnMelt').click();
    else if (e.key === 'r') $('btnResetRef').click();
  });

  const repo = $('repoLink');
  if (repo) repo.href = 'https://github.com/Glory0707/polymer-glass-lab';
}

/* ---------------- 启动 ---------------- */

function boot() {
  window.addEventListener('error', (e) => {
    const msg = e.error?.message || e.message || '';
    // 忽略无消息的资源加载错误与 ResizeObserver 的良性循环警告
    if (!msg || /ResizeObserver loop/.test(msg)) return;
    showFatal(msg);
  });
  try {
    bindUI();
    rebuild();
    window.__lab = { state, setTemperature, archiveGhost, refreshFit, setMode }; // 调试/无头测试句柄
    requestAnimationFrame(frame);
  } catch (err) {
    showFatal(err);
  }
}

boot();
