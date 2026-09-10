/**
 * main.js — 渲染、HUD 与 UI 接线
 * MD 内核运行在 Web Worker（sim.worker.js），本线程只做渲染与交互。
 */
import { GlassRenderer } from './renderer.js?v=30';
import { drawMSDPlot, drawHistoryPlot, drawA2Plot, drawVHPlot, drawStressPlot, drawProtoPlot } from './plots.js?v=30';
import { THERMAL_LUT } from './analysis.js?v=30';

const $ = (id) => document.getElementById(id);
const T_MIN = 0.05, T_MAX = 1.5;

const worker = new Worker(new URL('./sim.worker.js?v=30', import.meta.url), { type: 'module' });

/* 渲染所需的场景镜像（由 worker 消息填充） */
const view = {
  pos: new Float32Array(0),
  sigma: null,
  N: 0,
  Lx: 1, Ly: 1, Lz: 1,
  bondPairs: new Int32Array(0),
  msdPts: [], a2Pts: [], ghosts: [], history: [],
  refT: 1.0, fit: null, bins: [],
  density: 1.0,
};

const state = {
  numChains: 60,
  seed: 20260910,
  temperature: 1.0,
  stiffness: 0,
  smallFrac: 0,
  npt: false,
  targetP: 10,
  rate: 1e-3,
  speed: 16,
  paused: false,
  mode: 'free',
  showBonds: true,
  colorMode: 'mobility',
  perf: { fps: 0 },
  defRate: 0.02,
  defAmp: 0.12,
  defFreq: 0.5,
};

let renderer = null;
let frameNo = 0;
let lastFrame = performance.now();
let chainColors = null;

/* ---------------- 小工具 ---------------- */
const _rgb = [0, 0, 0];
function thermalLUTColor(x, out) {
  const idx = (Math.min(1, Math.max(0, x)) * 255) | 0;
  out[0] = Math.pow(THERMAL_LUT[idx * 3] * 0.62 + 0.38, 2.2);
  out[1] = Math.pow(THERMAL_LUT[idx * 3 + 1] * 0.62 + 0.38, 2.2);
  out[2] = Math.pow(THERMAL_LUT[idx * 3 + 2] * 0.62 + 0.38, 2.2);
}
function tColorCss(T, alpha = 1, lift = 0.35) {
  const x = Math.min(1, Math.max(0, (T - T_MIN) / (T_MAX - T_MIN)));
  const h = 220 - 210 * x;
  return alpha >= 1 ? `hsl(${h.toFixed(0)},85%,58%)` : `hsla(${h.toFixed(0)},85%,58%,${alpha})`;
}
const wsend = (msg, transfer) => worker.postMessage(msg, transfer || []);

function download(name, text, mime) {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([text], { type: mime }));
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}
const stamp = () => new Date().toISOString().slice(0, 19).replace(/[:T]/g, '-');

function showFatal(err) {
  console.error(err);
  const box = $('fatal');
  $('fatalMsg').textContent = String(err && err.message ? err.message : err);
  box.hidden = false;
}

/* ---------------- worker 消息 ---------------- */
worker.onmessage = (e) => {
  const m = e.data;
  switch (m.type) {
    case 'ready': {
      view.N = m.N;
      view.Lx = m.Lx; view.Ly = m.Ly; view.Lz = m.Lz;
      view.sigma = new Float32Array(m.sigma);
      view.bondPairs = new Int32Array(m.bondPairs);
      view.pos = new Float32Array(m.N * 3);
      if (renderer) { renderer.dispose(); renderer = null; }
      renderer = new GlassRenderer($('viewport'), view);
      applyChainColors();
      hideAnneal();
      $('seedVal').textContent = String(state.seed);
      break;
    }
    case 'anneal': {
      const box = $('anneal');
      box.hidden = false;
      $('annealPct').textContent = Math.round(m.pct * 100) + '%';
      break;
    }
    case 'anneal-done': {
      hideAnneal();
      break;
    }
    case 'frame': {
      view.pos = new Float32Array(m.pos);
      // NaN 防护
      for (let a = 0; a < view.pos.length; a++) {
        if (!Number.isFinite(view.pos[a])) view.pos[a] = 0;
      }
      view.chi = m.chi ? new Float32Array(m.chi) : view.chi;
      view.mob = m.mob ? new Float32Array(m.mob) : view.mob;
      if (m.vhBins) { view.vhBins = new Float32Array(m.vhBins); view.vhMax = m.vhMax; view.vhN = m.vhN; }
      view.density = m.stats.density;
      // 盒子尺寸变化（密度调整/NPT）：同步括号框与键镜像基准
      if (renderer && (Math.abs(m.stats.Lx - view.Lx) > 1e-6 || Math.abs(m.stats.Ly - view.Ly) > 1e-6 || Math.abs(m.stats.Lz - view.Lz) > 1e-6)) {
        view.Lx = m.stats.Lx; view.Ly = m.stats.Ly; view.Lz = m.stats.Lz;
        renderer.setBoxDims(view.Lx, view.Ly, view.Lz);
      }
      updateStats(m.stats);
      if (renderer) {
        updateColors();
        renderer.update({ showBonds: state.showBonds, showBox: $('boxChk').checked });
      }
      if (frameNo % 2 === 0) {
        drawMSDPlot($('msdPlot'), view.ghosts, { T: view.refT, pts: view.msdPts });
      }
      if (frameNo % 3 === 0) {
        drawHistoryPlot($('histPlot'), view.bins, fitSegments(view.fit, view.bins), {
          peLo: view.peLo ?? -3, peHi: view.peHi ?? 0,
          peTicks: peTicks(),
          raw: view.history,
          hasData: view.history.length > 0,
        });
      }
      if (frameNo % 4 === 0) drawA2Plot($('a2Plot'), view.a2Pts);
      if (frameNo % 4 === 2) drawStressPlot($('stressPlot'), view.stressPts);
      if (frameNo % 4 === 3) drawVHPlot($('vhPlot'), view.vhBins, view.vhMax, view.msdPts);
      if (frameNo % 4 === 1) drawProtoPlot($('protoPlot'), view.protoPts);
      frameNo++;
      break;
    }
    case 'samples': {
      view.stressPts = m.stressPts || [];
      view.protoPts = m.protoPts || [];
      view.msdPts = m.msdPts;
      view.a2Pts = m.a2Pts;
      view.ghosts = m.ghosts;
      view.history = m.history;
      view.refT = m.refT;
      view.fit = m.fit;
      view.bins = binByTLocal(m.history);
      const vTg = $('vTg');
      vTg.textContent = m.fit ? '≈ ' + m.fit.Tg.toFixed(2) : '—';
      vTg.style.color = m.fit ? 'var(--accent)' : '';
      break;
    }
    case 'proto-done': {
      const el = $('protoState');
      if (el) el.textContent = '完成';
      break;
    }
    case 'fatal':
      showFatal(m.msg);
      break;
  }
};
worker.onerror = (e) => showFatal(e.message || 'worker 错误');

function binByTLocal(history, width = 0.05) {
  const bins = [];
  const nb = Math.ceil(1.6 / width);
  for (let k = 0; k < nb; k++) {
    const lo = k * width, hi = lo + width;
    let msd = 0, n = 0;
    for (const s of history) {
      if (s.T >= lo && s.T < hi) { msd += s.msd; n++; }
    }
    if (n) bins.push({ T: lo + width / 2, msd: msd / n, n });
  }
  return bins;
}

function applyChainColors() {
  if (!renderer) return;
  const per = Math.floor(view.N / state.numChains);
  if (per < 1) return;
  chainColors = new Float32Array(view.N * 3);
  for (let c = 0; c < state.numChains; c++) {
    const h = ((c * 137.508 / 360) % 1) * 6;
    const seg = Math.floor(h), f = h - seg;
    const q = 1 - f;
    let r = 0.35, g = 0.35, b = 0.35;
    if (seg === 0) { r += 0.45 * q; g += 0.15 * q; }
    else if (seg === 1) { g += 0.35 * q; b += 0.1 * q; }
    else if (seg === 2) { g += 0.4 * q; }
    else if (seg === 3) { b += 0.4 * q; r += 0.1 * q; }
    else if (seg === 4) { b += 0.45 * q; r += 0.15 * q; }
    else { r += 0.4 * q; b += 0.2 * q; }
    const R = Math.pow(r, 2.2), G = Math.pow(g, 2.2), B = Math.pow(b, 2.2);
    for (let k = c * per; k < Math.min((c + 1) * per, view.N); k++) {
      chainColors[k * 3] = R; chainColors[k * 3 + 1] = G; chainColors[k * 3 + 2] = B;
    }
  }
  const ct2 = renderer.colorTarget;
  if (ct2) ct2.set(chainColors);
}

/* ---------------- 珠子着色 ---------------- */
function updateColors() {
  if (!renderer) return;
  const ct = renderer.colorTarget;
  if (!ct) return;
  if (state.colorMode === 'chain') {
    if (chainColors) ct.set(chainColors);
    return;
  }
  const src = state.colorMode === 'chi' ? view.chi : view.mob;
  if (!src) return;
  for (let i = 0; i < view.N; i++) {
    thermalLUTColor((Math.log10(src[i] + 1e-9) + 3) / 3.1, _rgb);
    ct[i * 3] = _rgb[0]; ct[i * 3 + 1] = _rgb[1]; ct[i * 3 + 2] = _rgb[2];
  }
}

/* ---------------- HUD ---------------- */
function updateStats(st) {
  const vT = $('vT');
  vT.textContent = st.T.toFixed(2);
  vT.style.color = tColorCss(st.T);
  $('vTmeas').textContent = '(' + st.Tmeas.toFixed(2) + ')';
  $('vPE').textContent = st.pe.toFixed(2);
  const lastMsd = view.msdPts.length ? view.msdPts[view.msdPts.length - 1][1] : st.msd;
  $('vMSD').textContent = lastMsd.toFixed(2);
  $('vTau').textContent = st.tau.toFixed(0);
  $('vFps').textContent = Math.round(state.perf.fps) + 'fps';
  $('densityVal').textContent = st.density.toFixed(2);
}

function peTicks() {
  if (!view.history.length) return [];
  const peLo = view.peLo ?? -3, peHi = view.peHi ?? 0;
  if (!isFinite(peLo) || !isFinite(peHi) || peHi <= peLo) return [];
  const mid = (peLo + peHi) / 2;
  return [[peLo, peLo.toFixed(1)], [mid, mid.toFixed(1)], [peHi, peHi.toFixed(1)]];
}

function fitSegments(fit, bins) {
  if (!fit || bins.length < 2) return null;
  const t0 = bins[0].T, tSplit = fit.splitT, t1 = bins[bins.length - 1].T;
  const at = (line, T) => Math.pow(10, line.a * T + line.b);
  return {
    seg1: [t0, at(fit.glass, t0), tSplit, at(fit.glass, tSplit)],
    seg2: [tSplit, at(fit.liquid, tSplit), t1, at(fit.liquid, t1)],
    Tg: fit.Tg,
  };
}

/* ---------------- UI 接线 ---------------- */
function setMode(mode) {
  state.mode = mode;
  $('btnCool').classList.toggle('active', mode === 'cool');
  $('btnHeat').classList.toggle('active', mode === 'heat');
  wsend({ cmd: 'mode', v: mode });
}

function bindUI() {
  const tempSlider = $('tempSlider');
  tempSlider.addEventListener('input', () => {
    setMode('free');
    const T = parseFloat(tempSlider.value);
    state.temperature = T;
    wsend({ cmd: 'temp', T });
    const v = $('vT');
    v.textContent = T.toFixed(2);
    v.style.color = tColorCss(T);
  });
  tempSlider.addEventListener('change', () => wsend({ cmd: 'archive' }));

  $('btnCool').addEventListener('click', () => setMode('cool'));
  $('btnHeat').addEventListener('click', () => setMode('heat'));

  const rateSlider = $('rateSlider');
  const showRate = () => {
    $('rateVal').textContent = state.rate >= 1e-2
      ? state.rate.toFixed(1) + '×10⁻²'
      : (state.rate * 1e3).toFixed(1) + '×10⁻³';
  };
  rateSlider.addEventListener('input', () => {
    state.rate = Math.pow(10, -4 + (parseInt(rateSlider.value, 10) / 100) * 2.3);
    showRate();
    wsend({ cmd: 'rate', v: state.rate });
  });
  state.rate = Math.pow(10, -4 + (43 / 100) * 2.3);
  showRate();

  $('btnQuench').addEventListener('click', () => {
    setMode('free'); state.temperature = T_MIN;
    tempSlider.value = String(T_MIN);
    wsend({ cmd: 'temp', T: T_MIN });
    wsend({ cmd: 'archive' });
  });
  $('btnMelt').addEventListener('click', () => {
    setMode('free'); state.temperature = T_MAX;
    tempSlider.value = String(T_MAX);
    wsend({ cmd: 'temp', T: T_MAX });
    wsend({ cmd: 'archive' });
  });

  const pauseBtn = $('btnPause');
  pauseBtn.addEventListener('click', () => {
    state.paused = !state.paused;
    pauseBtn.textContent = state.paused ? '继续' : '暂停';
    wsend({ cmd: 'pause', v: state.paused });
  });

  $('speedSlider').addEventListener('input', (e) => {
    state.speed = parseInt(e.target.value, 10);
    $('speedVal').textContent = String(state.speed);
    wsend({ cmd: 'speed', v: state.speed });
  });

  $('sizeSel').addEventListener('change', (e) => {
    state.numChains = parseInt(e.target.value, 10);
    sendRebuild();
  });

  $('compSel').addEventListener('change', (e) => {
    state.smallFrac = parseFloat(e.target.value);
    sendRebuild();
  });
  $('densitySlider').addEventListener('input', (e) => {
    wsend({ cmd: 'density-target', v: parseFloat(e.target.value) });
  });

  $('nptChk').addEventListener('change', (e) => {
    state.npt = e.target.checked;
    wsend({ cmd: 'npt', v: state.npt, p0: state.targetP });
  });
  $('p0Slider').addEventListener('input', (e) => {
    state.targetP = parseFloat(e.target.value);
    $('p0Val').textContent = String(state.targetP);
    wsend({ cmd: 'npt', v: state.npt, p0: state.targetP });
  });

  const stiffnessSlider = $('stiffnessSlider');
  stiffnessSlider.addEventListener('input', () => {
    state.stiffness = parseFloat(stiffnessSlider.value);
    $('stiffnessVal').textContent = String(parseFloat(state.stiffness.toFixed(2)));
    wsend({ cmd: 'stiffness', v: state.stiffness });
  });

  $('btnNew').addEventListener('click', () => {
    state.seed = (Math.random() * 0x7fffffff) | 0;
    sendRebuild();
  });
  $('btnResetRef').addEventListener('click', () => wsend({ cmd: 'archive' }));

  $('colorSel').addEventListener('change', (e) => { state.colorMode = e.target.value; });
  $('bondsChk').addEventListener('change', (e) => { state.showBonds = e.target.checked; });
  $('boxChk').addEventListener('change', (e) => { state.showBox = e.target.checked; });

  $('figsToggle').addEventListener('click', () => $('figs').classList.toggle('open'));
  $('figsClose').addEventListener('click', () => $('figs').classList.remove('open'));

  const pop = $('paramsPop');
  const ptoggle = $('paramsToggle');
  ptoggle.addEventListener('click', (e) => {
    e.stopPropagation();
    const open = pop.classList.toggle('hidden') === false;
    ptoggle.classList.toggle('active', open);
  });
  document.addEventListener('click', (e) => {
    if (pop.classList.contains('hidden')) return;
    if (pop.contains(e.target) || ptoggle.contains(e.target)) return;
    pop.classList.add('hidden');
    ptoggle.classList.remove('active');
  });
  window.addEventListener('keydown', (e) => {
    if (e.key === 'Escape') {
      pop.classList.add('hidden');
      ptoggle.classList.remove('active');
      $('figs').classList.remove('open');
    }
  });

  window.addEventListener('keydown', (e) => {
    if (e.target.tagName === 'INPUT' || e.target.tagName === 'SELECT') return;
    if (e.code === 'Space') { e.preventDefault(); pauseBtn.click(); }
    else if (e.key === 'c') $('btnQuench').click();
    else if (e.key === 'h') $('btnMelt').click();
  });

  const defModeSel = $('defModeSel');
  defModeSel.addEventListener('change', (e) => {
    const mode = e.target.value;
    wsend({ cmd: 'deform', mode, rate: state.defRate, amp: state.defAmp, freq: state.defFreq });
    const grp = document.getElementById('defRateGroup');
    if (grp) grp.hidden = mode === 'none';
    const grp2 = document.getElementById('defCycGroup');
    if (grp2) grp2.hidden = mode !== 'cyclic';
  });
  $('defRateSlider').addEventListener('input', (e) => {
    state.defRate = parseFloat(e.target.value);
    $('defRateVal').textContent = state.defRate.toFixed(3);
    wsend({ cmd: 'deform', mode: defModeSel.value, rate: state.defRate, amp: state.defAmp, freq: state.defFreq });
  });
  $('defAmpSlider').addEventListener('input', (e) => {
    state.defAmp = parseFloat(e.target.value);
    $('defAmpVal').textContent = state.defAmp.toFixed(2);
    wsend({ cmd: 'deform', mode: defModeSel.value, rate: state.defRate, amp: state.defAmp, freq: state.defFreq });
  });
  $('defFreqSlider').addEventListener('input', (e) => {
    state.defFreq = parseFloat(e.target.value);
    $('defFreqVal').textContent = state.defFreq.toFixed(2);
    wsend({ cmd: 'deform', mode: defModeSel.value, rate: state.defRate, amp: state.defAmp, freq: state.defFreq });
  });
  $('btnRelease').addEventListener('click', () => wsend({ cmd: 'deform-release' }));
  $('btnProto').addEventListener('click', () => {
    const seq = $('protoSeq').value.split(',').map((p2) => p2.trim().split(/\s+/)).filter((p2) => p2.length === 2)
      .map((p2) => ({ T: Math.min(1.5, Math.max(0.05, parseFloat(p2[0]))), dur: Math.round(parseFloat(p2[1]) * 125) }));
    if (seq.length) wsend({ cmd: 'protocol', seq });
    const psEl = $('protoState');
    if (psEl) psEl.textContent = '运行中';
  });
  $('btnProtoStop').addEventListener('click', () => {
    wsend({ cmd: 'protocol-stop' });
    const psEl = $('protoState');
    if (psEl) psEl.textContent = '';
  });
  $('btnCsv').addEventListener('click', exportCsv);
  $('btnJson').addEventListener('click', exportJson);
}

function sendRebuild() {
  showAnneal();
  wsend({
    cmd: 'rebuild',
    numChains: state.numChains, seed: state.seed,
    temperature: state.temperature, smallFrac: state.smallFrac,
    stiffness: state.stiffness, npt: state.npt, targetP: state.targetP,
    annealSteps: 4000,
  });
}

function showAnneal() {
  $('anneal').hidden = false;
  $('annealPct').textContent = '0%';
}
function hideAnneal() {
  $('anneal').hidden = true;
  $('seedVal').textContent = String(state.seed);
}

/* ---------------- 导出 ---------------- */
function exportCsv() {
  const L = [];
  L.push('# polymer-glass-lab export');
  L.push('# N=' + view.N + ' chains=' + state.numChains + ' chainLen=40 density=' + view.density.toFixed(3) + ' kappa=' + state.stiffness + ' seed=' + state.seed);
  L.push('');
  L.push('# msd-tau (refT=' + view.refT.toFixed(3) + ')');
  L.push('tau,msd');
  for (const [t, v] of view.msdPts) L.push(t.toFixed(4) + ',' + v.toPrecision(6));
  L.push('');
  L.push('# thermal-history');
  L.push('T,msd20,pe');
  for (const s of view.history) L.push(s.T.toFixed(3) + ',' + s.msd.toPrecision(6) + ',' + s.pe.toPrecision(6));
  L.push('');
  L.push('# nongaussian-a2');
  L.push('tau,a2');
  for (const [t, a] of view.a2Pts) L.push(t.toFixed(4) + ',' + a.toPrecision(6));
  download('polymer-glass-' + stamp() + '.csv', L.join(String.fromCharCode(10)), 'text/csv');
}

function exportJson() {
  const data = {
    meta: {
      N: view.N, chains: state.numChains, chainLen: 40,
      density: view.density, stiffness: state.stiffness, seed: state.seed,
      npt: state.npt, targetP: state.targetP,
      exportedAt: new Date().toISOString(),
    },
    msdTau: { refT: view.refT, pts: view.msdPts },
    nongaussianA2: view.a2Pts,
    stressResponse: view.stressPts,
    memoryProtocol: view.protoPts,
    thermalHistory: view.history,
  };
  download('polymer-glass-' + stamp() + '.json', JSON.stringify(data, null, 2), 'application/json');
}

/* ---------------- 帧循环（渲染 + HUD） ---------------- */
function frame(now) {
  requestAnimationFrame(frame);
  const dtms = now - lastFrame;
  lastFrame = now;
  if (dtms > 0 && dtms < 500) state.perf.fps = state.perf.fps * 0.92 + (1000 / dtms) * 0.08;
  frameNo++;
  if (renderer) renderer.update({ showBonds: state.showBonds, showBox: $('boxChk').checked });
}

/* ---------------- 启动 ---------------- */
function boot() {
  window.addEventListener('error', (e) => {
    const msg = (e.error && e.error.message) || e.message || '';
    if (!msg || /ResizeObserver loop/.test(msg)) return;
    showFatal(msg);
  });
  try {
    bindUI();
    sendRebuild();
    requestAnimationFrame(frame);
    $('repoLink').href = 'https://github.com/Glory0707/polymer-glass-lab';
  } catch (err) {
    showFatal(err);
  }
}

boot();
window.__lab = { view, state, worker }; // 调试句柄
