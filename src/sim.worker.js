/**
 * sim.worker.js — MD 内核宿主：模拟循环、采样、热历史与 Tg 拟合全部在此线程运行，
 * 主线程只负责渲染与 UI。与主线程通过消息通信。
 *
 * 协议：
 *   main → worker: {cmd:'init'|'rebuild'|'temp'|'archive'|'mode'|'rate'|'speed'|'pause'
 *                   |'stiffness'|'npt'|'density-target'|'deform'|'deform-release'
 *                   |'protocol'|'protocol-stop'}
 *   worker → main: {type:'ready'|'anneal'|'frame'|'samples'|'proto-done'|'fatal'}
 */
import { KGSim } from './md.js?v=22';
import { binByT, twoSegmentFit } from './analysis.js?v=22';

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
};
let annealLeft = 0;
let refT = 1.0;
let msdPts = [], a2Pts = [], ghosts = [], history = [], stressPts = [], protoPts = [];
let nextSampleStep = 8, lastRecStep = 0, lastStressStep = 0;
let protocol = null, protoIdx = -1, protoStepLeft = 0, protoT0 = 0;
let emaDev = 0, emaXY = 0;
let tickNo = 0;

const post = (msg, transfer) => self.postMessage(msg, transfer || []);

function archive() {
  if (msdPts.length > 5) {
    ghosts.push({ T: refT, pts: msdPts });
    while (ghosts.length > 3) ghosts.shift();
  }
  sim.resetRef();
  refT = sim.T;
  msdPts = [];
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
  nextSampleStep = Math.max(nextSampleStep + 8, Math.ceil(nextSampleStep * 1.12));
  if (msdPts.length > 500) { msdPts.shift(); a2Pts.shift(); }
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
    stressPts.push([sim.deform.strain, emaDev]);
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

  if (protocol) { tickProtocol(); pushFrame(); pushSamples(); return; }

  if (annealLeft > 0) {
    const chunk = Math.min(annealLeft, 600);
    sim.anneal(chunk);
    annealLeft -= chunk;
    post({ type: 'anneal', pct: Math.min(1, 1 - annealLeft / 4000) });
    if (annealLeft <= 0) { archive(); post({ type: 'ready' }); }
    return;
  }

  if (cfg.paused) { pushSamples(); return; }

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
  // van Hove 直方图（每 5 帧）
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
  const transfers = [posCopy.buffer, mob.buffer];
  if (chi) transfers.push(chi.buffer);
  if (vhBins) transfers.push(vhBins.buffer);
  post({
    type: 'frame',
    pos: posCopy.buffer,
    mob: mob.buffer,
    chi: chi ? chi.buffer : null,
    vhBins: vhBins ? vhBins.buffer : null,
    vhMax,
    vhN: sim.N,
    stats: {
      tau: sim.time, T: sim.T, Tmeas: sim.keTemp,
      pe: sim.pePerBead, msd: meanD2,
      density: sim.density, strain: sim.deform.strain,
      stressDev: emaDev, stressXY: emaXY,
      deformMode: cfg.deform.mode, deformGamma: sim.deform.gamma,
    },
  }, transfers);
}

function pushSamples() {
  post({
    type: 'samples',
    msdPts, a2Pts, ghosts, history,
    stressPts, protoPts,
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
          annealSteps: 0,
        });
        annealLeft = m.annealSteps;
        msdPts = []; a2Pts = []; ghosts = []; history = []; stressPts = []; protoPts = [];
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
      case 'reset-ref': sim.resetRef(); msdPts = []; nextSampleStep = 8; break;
      case 'stiffness': sim.stiffness = m.v; break;
      case 'npt': sim.npt = m.v; sim.targetP = m.p0; break;
      case 'density-target': cfg.densityTarget = m.v; break;
      case 'deform': cfg.deform.mode = m.mode; cfg.deform.rate = m.rate; cfg.deform.amp = m.amp; cfg.deform.freq = m.freq; break;
      case 'deform-release': cfg.deform.mode = 'none'; cfg.deform.target = 0; break;
      case 'protocol': protocol = m.seq; protoIdx = -1; protoStepLeft = 0; break;
      case 'protocol-stop': protocol = null; break;
    }
  } catch (err) {
    emitFatal(String(err && err.message ? err.message : err));
  }
};

schedule();
