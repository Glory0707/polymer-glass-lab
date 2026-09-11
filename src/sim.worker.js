/**
 * sim.worker.js — MD 内核宿主：模拟循环、采样、热历史与 Tg 拟合全部在此线程运行，
 * 主线程只负责渲染与 UI。与主线程通过消息通信。
 *
 * 协议：
 *   main → worker: {cmd:'init'|'rebuild'|'temp'|'archive'|'mode'|'rate'|'speed'|'pause'
 *                   |'stiffness'|'npt'|'density-target'|'deform'|'deform-release'
 *                   |'protocol'|'protocol-stop'|'heat-brush'|'grab'|'grab-move'|'grab-release'}
 *   worker → main: {type:'ready'|'anneal'|'anneal-done'|'frame'|'samples'|'proto-done'|'fatal'}
 */
import { KGSim } from './md.js?v=37';
import { binByT, twoSegmentFit, tauFromMsd, vftFit } from './analysis.js?v=37';

let sim = null;
const cfg = {
  speed: 16,
  paused: false,
  mode: 'free',          // 温度扫描：free | cool | heat
  rate: 1e-3,
  T_MIN: 0.05,
  T_MAX: 1.5,
  densityTarget: null,
  deform: { mode: 'none', rate: 0.02, amp: 0.12, freq: 0.5 },
  protocol: null,        // 记忆实验温度阶梯 [{T, dur}]
  film: false,
};
let annealLeft = 0;
let refT = 1.0;
let msdPts = [], a2Pts = [], fsqPts = [], ghosts = [], history = [], stressPts = [], orientPts = [], protoPts = [];
let vftPts = [];         // [{T, tau}] α 弛豫时间随温度（VFT 图）
let nextSampleStep = 8, lastRecStep = 0, lastStressStep = 0;
let protocol = null, protoIdx = -1, protoStepLeft = 0, protoT0 = 0;
let emaDev = 0, emaXY = 0;
let tickNo = 0;

const post = (msg, transfer) => self.postMessage(msg, transfer || []);

function archive() {
  // VFT：换参考前若当前窗口已测到 τα（MSD 穿越 1σ²），记入 {T, tau}
  const tau = tauFromMsd(msdPts);
  if (tau != null) {
    vftPts = vftPts.filter((p) => Math.abs(p.T - refT) > 0.02);
    vftPts.push({ T: refT, tau });
    if (vftPts.length > 90) vftPts.shift();
  }
  if (msdPts.length > 5) {
    ghosts.push({ T: refT, pts: msdPts });
    while (ghosts.length > 3) ghosts.shift();
  }
  sim.resetRef();
  refT = sim.T;
  msdPts = [];
  a2Pts = [];   // α₂ 与 MSD 同窗口同参考，必须一起清，否则多温度段拼成乱线
  fsqPts = [];  // Fs(q,t) 同理
  nextSampleStep = 8;
}

function sample() {
  const st = sim.stepCount - sim.refStep;
  if (st < nextSampleStep) return;
  const tau = st * sim.dt;
  const msd = sim.msdRef();
  const m4 = sim.msd4Ref();
  const a2 = msd > 1e-9 ? (3 * m4) / (5 * msd * msd) - 1 : 0;
  msdPts.push([tau, msd]);
  a2Pts.push([tau, a2]);
  fsqPts.push([tau, sim.fsqRef()]);
  nextSampleStep = Math.max(nextSampleStep + 8, Math.ceil(nextSampleStep * 1.12));
  if (msdPts.length > 500) { msdPts.shift(); a2Pts.shift(); fsqPts.shift(); }
}

function record() {
  if (sim.stepCount - lastRecStep < 150) return;
  lastRecStep = sim.stepCount;
  if (sim.lagAge() < sim.lagSteps * sim.dt * 0.8) return;
  history.push({ T: sim.T, msd: sim.msdLag(), pe: sim.pePerBead });
  if (history.length > 6000) history = history.filter((_, i) => i % 2 === 0);
}

function sampleStress() {
  if (sim.stepCount - lastStressStep < 100) return;
  lastStressStep = sim.stepCount;
  const ema = 0.08;
  if (cfg.deform.mode === 'shear') {
    emaXY += (sim.stressXY - emaXY) * ema;
    stressPts.push([sim.deform.gamma, emaXY]);
  } else {
    emaDev += (sim.stressDev - emaDev) * ema;
    // 循环模式的应变是 epsCur（strain 只累计单轴拉伸）
    const eps = cfg.deform.mode === 'cyclic' ? sim.deform.epsCur : sim.deform.strain;
    stressPts.push([eps, emaDev]);
    // 键取向 P2（应力光学对应）：与应力同步采样
    if (cfg.deform.mode !== 'none') {
      orientPts.push([eps, sim.bondP2()]);
      if (orientPts.length > 500) orientPts.shift();
    }
  }
  if (stressPts.length > 500) stressPts.shift();
}

function tickProtocol() {
  if (!protocol) return;
  const seg = protocol[Math.min(protoIdx, protocol.length - 1)];
  if (protoStepLeft <= 0) {
    protoIdx++;
    if (protoIdx >= protocol.length) {
      protocol = null;
      post({ type: 'proto-done' });
      return;
    }
    sim.T = protocol[protoIdx].T;
    refT = sim.T;
    archive();
    protoStepLeft = protocol[protoIdx].dur;
    protoT0 = sim.time;
  }
  const chunk = Math.min(protoStepLeft, 400);
  for (let i = 0; i < chunk; i++) {
    sim.step();
    sample();
    record();
    if (i % 50 === 0) protoPts.push([sim.time - protoT0, sim.pePerBead, sim.T]);
  }
  protoStepLeft -= chunk;
  if (protoPts.length > 3000) protoPts = protoPts.filter((_, i) => i % 2 === 0);
}

function tick() {
  schedule();
  if (!sim) return;
  tickNo++;

  // 制备中的退火优先于协议；暂停冻结协议推进
  if (annealLeft > 0) {
    const chunk = Math.min(annealLeft, 600);
    sim.anneal(chunk);
    annealLeft -= chunk;
    post({ type: 'anneal', pct: Math.min(1, 1 - annealLeft / 4000) });
    if (annealLeft <= 0) { archive(); post({ type: 'anneal-done' }); }
    return;
  }

  if (cfg.paused) { pushSamples(); return; }

  if (protocol) { tickProtocol(); pushFrame(); pushSamples(); return; }

  if (cfg.densityTarget != null) {
    const remain = Math.log(cfg.densityTarget / sim.density);
    if (Math.abs(remain) < 0.004) {
      sim.setDensity(cfg.densityTarget);
      cfg.densityTarget = null;
      archive();
    } else {
      sim.setDensity(sim.density * Math.exp(remain * 0.1));
    }
  }

  const n = cfg.speed;
  // 卸载完成同步：md 侧 stretch 回到 ε=0 后自动转 none
  if (cfg.deform.target === 0 && sim.deform.mode === 'none' && cfg.deform.mode === 'stretch') {
    cfg.deform.mode = 'none'; cfg.deform.target = null;
  }
  for (let s = 0; s < n; s++) {
    if (cfg.mode === 'cool') {
      const next = sim.T - cfg.rate * sim.dt;
      sim.T = next <= cfg.T_MIN ? cfg.T_MIN : next;
      if (sim.T <= cfg.T_MIN) cfg.mode = 'free';
    } else if (cfg.mode === 'heat') {
      const next = sim.T + cfg.rate * sim.dt;
      sim.T = next >= cfg.T_MAX ? cfg.T_MAX : next;
      if (sim.T >= cfg.T_MAX) cfg.mode = 'free';
    }
    if (cfg.mode !== 'free' && Math.abs(sim.T - refT) > 0.15) archive();
    sim.step();
    sample();
    record();
    sampleStress();
  }

  pushFrame();
  if (tickNo % 30 === 0) pushSamples();
}

function pushFrame() {
  const posCopy = sim.pos.slice();
  // NaN 防护：不发送含 NaN 的位置
  for (let a = 0; a < posCopy.length; a++) {
    if (!Number.isFinite(posCopy[a])) { console.error('NaN in pos[' + a + ']'); posCopy[a] = 0; }
  }
  const u = sim.upos, sm = sim.snapMob;
  const mob = new Float32Array(sim.N);
  let sumD2 = 0;
  for (let i = 0; i < sim.N; i++) {
    const dx = u[i * 3] - sm[i * 3], dy = u[i * 3 + 1] - sm[i * 3 + 1], dz = u[i * 3 + 2] - sm[i * 3 + 2];
    const d2 = dx * dx + dy * dy + dz * dz;
    mob[i] = d2;
    sumD2 += d2;
  }
  const meanD2 = sumD2 / sim.N;
  let chi = null;
  if (tickNo % 20 === 0) chi = sim.smoothMobility();
  // van Hove 直方图（每 10 帧）
  let vhBins = null, vhMax = 3.5;
  if (tickNo % 10 === 0) {
    vhBins = new Float64Array(28);
    for (let i = 0; i < sim.N; i++) {
      const dx = u[i * 3] - sm[i * 3], dy = u[i * 3 + 1] - sm[i * 3 + 1], dz = u[i * 3 + 2] - sm[i * 3 + 2];
      const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
      const b = Math.min(27, (r / vhMax * 28) | 0);
      vhBins[b]++;
    }
    for (let b = 0; b < 28; b++) vhBins[b] /= sim.N * (vhMax / 28);
  }
  // 薄膜迁移率剖面：z 分 24 层（限膜内范围），层内平均位移²（每 15 帧）
  let mobProf = null;
  if (cfg.film && tickNo % 15 === 0) {
    const NB = 24, Lz = sim.Lz, m = sim.wallMargin;
    const zLo = m - 0.5, zHi = Lz - m + 0.5;
    const acc = new Float64Array(NB), cnt = new Float64Array(NB);
    for (let i = 0; i < sim.N; i++) {
      const i3 = i * 3;
      const b = Math.min(NB - 1, Math.max(0, ((sim.pos[i3 + 2] - zLo) / (zHi - zLo) * NB) | 0));
      const dx = u[i3] - sm[i3], dy = u[i3 + 1] - sm[i3 + 1], dz = u[i3 + 2] - sm[i3 + 2];
      acc[b] += dx * dx + dy * dy + dz * dz;
      cnt[b]++;
    }
    const prof = new Float64Array(NB);
    for (let b = 0; b < NB; b++) prof[b] = cnt[b] ? acc[b] / cnt[b] : 0;
    mobProf = prof;
  }
  const transfers = [posCopy.buffer, mob.buffer];
  if (chi) transfers.push(chi.buffer);
  if (vhBins) transfers.push(vhBins.buffer);
  if (mobProf) transfers.push(mobProf.buffer);
  post({
    type: 'frame',
    pos: posCopy.buffer,
    mob: mob.buffer,
    chi: chi ? chi.buffer : null,
    vhBins: vhBins ? vhBins.buffer : null,
    vhMax,
    vhN: sim.N,
    mobProf: mobProf ? mobProf.buffer : null,
    stats: {
      tau: sim.time, T: sim.T, Tmeas: sim.keTemp,
      pe: sim.pePerBead, msd: meanD2,
      density: sim.density, strain: sim.deform.strain,
      Lx: sim.Lx, Ly: sim.Ly, Lz: sim.Lz,
      stressDev: emaDev, stressXY: emaXY,
      deformMode: cfg.deform.mode, deformGamma: sim.deform.gamma,
    },
  }, transfers);
}

function pushSamples() {
  const vf = vftFit(vftPts);
  post({
    type: 'samples',
    msdPts, a2Pts, ghosts, history,
    stressPts, orientPts, protoPts,
    fsqPts,
    vftPts: vf.points, vftFit: vf.fit, vftArr: vf.arrhenius,
    refT, mode: cfg.mode,
    protocolActive: !!protocol,
    fit: fitFor(history),
    q6Mean: sim.q6Mean ?? null,
  });
}

function fitFor(hist) {
  const bins = binByT(hist, 0.05, 0, 1.6, 1);
  let fit = twoSegmentFit(bins, true);
  if (!fit && hist.length >= 30) fit = twoSegmentFit(binByT(hist, 0.12, 0, 1.6, 1), true);
  if (!fit) return null;
  return {
    Tg: fit.Tg, splitT: fit.splitT,
    glass: { a: fit.glass.a, b: fit.glass.b },
    liquid: { a: fit.liquid.a, b: fit.liquid.b },
  };
}

function schedule() { setTimeout(tick, 1000 / 60); }

function emitFatal(msg) { post({ type: 'fatal', msg }); }

function emitReady() {
  const sigmaCopy = sim.sigma.slice();
  const bondCopy = sim.bondPairs.slice();
  post({
    type: 'ready',
    N: sim.N, Lx: sim.Lx, Ly: sim.Ly, Lz: sim.Lz,
    film: sim.film,
    sigma: sigmaCopy.buffer,
    bondPairs: bondCopy.buffer,
  }, [sigmaCopy.buffer, bondCopy.buffer]);
}

self.onmessage = (e) => {
  const m = e.data;
  try {
    switch (m.cmd) {
      case 'init':
      case 'rebuild': {
        sim = new KGSim({
          numChains: m.numChains, chainLen: 40, seed: m.seed,
          temperature: m.temperature, smallFrac: m.smallFrac,
          stiffness: m.stiffness, npt: m.npt, targetP: m.targetP,
          film: m.film ?? false,
          annealSteps: 0,
        });
        cfg.film = sim.film;
        annealLeft = m.annealSteps;
        msdPts = []; a2Pts = []; fsqPts = []; ghosts = []; history = [];
        stressPts = []; orientPts = []; protoPts = []; vftPts = [];
        nextSampleStep = 8; lastRecStep = 0; lastStressStep = 0;
        refT = sim.T;
        protocol = null;
        emitReady();
        break;
      }
      case 'anneal-chunk': sim.anneal(m.steps || 400); break;
      case 'speed': cfg.speed = m.v; break;
      case 'pause': cfg.paused = m.v; break;
      case 'mode': cfg.mode = m.v; break;
      case 'rate': cfg.rate = m.v; break;
      case 'temp': sim.T = m.T; break;
      case 'archive': archive(); break;
      case 'stiffness': sim.stiffness = m.v; break;
      case 'npt': sim.npt = m.v; sim.targetP = m.p0; break;
      case 'density-target': cfg.densityTarget = m.v; break;
      case 'deform':
        cfg.deform.mode = m.mode; cfg.deform.rate = m.rate; cfg.deform.amp = m.amp; cfg.deform.freq = m.freq;
        // _deformStep 读的是 sim.deform，必须同步，否则形变永远不会发生
        if (sim) { sim.deform.mode = m.mode; sim.deform.rate = m.rate; sim.deform.amp = m.amp; sim.deform.freq = m.freq; }
        break;
      case 'deform-release':
        // 释放 = 平滑卸载：反向拉伸回 ε=0（应力-应变图画出卸载曲线），回到 0 后自动停
        cfg.deform.target = 0;
        if (sim) {
          if (sim.deform.mode === 'stretch') {
            sim.deform.target = 0;
            sim.deform.rate = -Math.abs(sim.deform.rate) * 5;
          } else {
            sim.deform.mode = 'none';
          }
        }
        break;
      case 'protocol': protocol = m.seq; protoIdx = -1; protoStepLeft = 0; break;
      case 'protocol-stop': protocol = null; break;
      case 'heat-brush': if (sim) sim.heatBrush(m.x, m.y, m.z, m.r ?? 2.5, m.dT ?? 0.15); break;
      case 'grab': if (sim) sim.grabPick(m.x, m.y, m.z, m.r ?? 2.2); break;
      case 'grab-move': if (sim) sim.grabMove(m.x, m.y, m.z); break;
      case 'grab-release': if (sim) sim.grabRelease(); break;
      case 'reset-ref': sim.resetRef(); msdPts = []; a2Pts = []; fsqPts = []; nextSampleStep = 8; break;
    }
  } catch (err) {
    emitFatal(String(err && err.message ? err.message : err));
  }
};

schedule();
