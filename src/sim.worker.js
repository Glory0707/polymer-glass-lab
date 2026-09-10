/**
 * sim.worker.js — MD 内核宿主：模拟循环、采样、热历史与 Tg 拟合全部在此线程运行，
 * 主线程只负责渲染与 UI。与主线程通过消息通信。
 */
import { KGSim } from './md.js?v=19';
import { binByT, twoSegmentFit } from './analysis.js?v=19';

let sim = null;
const cfg = {
  speed: 16,
  paused: false,
  mode: 'free',        // free | cool | heat
  rate: 1e-3,          // ε/τ
  T_MIN: 0.05,
  T_MAX: 1.5,
  densityTarget: null,
};
let annealLeft = 0;
let refT = 1.0;
let msdPts = [], a2Pts = [], ghosts = [], history = [];
let nextSampleStep = 8, lastRecStep = 0;
let tickNo = 0;

const post = (msg, transfer) => self.postMessage(msg, transfer || []);

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

function tick() {
  schedule();
  if (!sim) return;
  tickNo++;

  // 退火分块
  if (annealLeft > 0) {
    const chunk = Math.min(annealLeft, 600);
    sim.anneal(chunk);
    annealLeft -= chunk;
    post({ type: 'anneal', pct: 1 - annealLeft / 4000 });
    if (annealLeft <= 0) { archive(); post({ type: 'ready' }); }
    return;
  }

  if (cfg.paused) return;

  // 密度渐变逼近目标
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
  }

  if (tickNo % 2 === 0) {
    const posCopy = sim.pos.slice();
    // 珠色数据：迁移率 x（每帧）与 χ₄ 平滑（每 10 帧）
    let mob = null, chi = null;
    const u = sim.upos, sm = sim.snapMob;
    mob = new Float32Array(sim.N);
    for (let i = 0; i < sim.N; i++) {
      const dx = u[i * 3] - sm[i * 3], dy = u[i * 3 + 1] - sm[i * 3 + 1], dz = u[i * 3 + 2] - sm[i * 3 + 2];
      mob[i] = (Math.log10(dx * dx + dy * dy + dz * dz + 1e-9) + 3) / 3.1;
    }
    if (tickNo % 20 === 0) chi = sim.smoothMobility();
    const transfers = [posCopy.buffer];
    if (mob) transfers.push(mob.buffer);
    if (chi) transfers.push(chi.buffer);
    post({
      type: 'frame',
      pos: posCopy.buffer,
      mob: mob ? mob.buffer : null,
      chi: chi ? chi.buffer : null,
      stats: {
        tau: sim.time, T: sim.T, Tmeas: sim.keTemp,
        pe: sim.pePerBead, msd: sim.msdRef(),
        density: sim.density, npt: sim.npt,
        a2: a2Pts.length ? a2Pts[a2Pts.length - 1][1] : 0,
      },
    }, transfers);
  }
  if (tickNo % 30 === 0) {
    post({
      type: 'samples',
      msdPts, a2Pts, ghosts, history,
      refT,
      fit: fitFor(history),
      mode: cfg.mode,
    });
  }
}

function schedule() { setTimeout(tick, 1000 / 60); }

function fitFor(hist) {
  const bins = binByT(hist, 0.05, 0, 1.6, 1);
  let fit = twoSegmentFit(bins);
  if (!fit && hist.length >= 30) fit = twoSegmentFit(binByT(hist, 0.12, 0, 1.6, 1));
  if (!fit) return null;
  return {
    Tg: fit.Tg, splitT: fit.splitT,
    glass: { a: fit.glass.a, b: fit.glass.b },
    liquid: { a: fit.liquid.a, b: fit.liquid.b },
  };
}

self.onmessage = (e) => {
  const m = e.data;
  switch (m.cmd) {
    case 'init': {
      sim = new KGSim({
        numChains: m.numChains, chainLen: 40, seed: m.seed,
        temperature: m.temperature, smallFrac: m.smallFrac,
        stiffness: m.stiffness, npt: m.npt, targetP: m.targetP,
        annealSteps: 0,
      });
      annealLeft = m.annealSteps;
      msdPts = []; a2Pts = []; ghosts = []; history = [];
      nextSampleStep = 8; lastRecStep = 0; refT = sim.T;
      const sigmaCopy = sim.sigma.slice();
      const bondCopy = sim.bondPairs.slice();
      post({
        type: 'ready',
        N: sim.N, Lx: sim.Lx, Ly: sim.Ly, Lz: sim.Lz,
        sigma: sigmaCopy.buffer,
        bondPairs: bondCopy.buffer,
      }, [sigmaCopy.buffer, bondCopy.buffer]);
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
    case 'rebuild': {
      const keepT = sim ? sim.T : 1.0;
      sim = new KGSim({
        numChains: m.numChains, chainLen: 40, seed: m.seed,
        temperature: m.temperature, smallFrac: m.smallFrac,
        stiffness: m.stiffness, npt: m.npt, targetP: m.targetP,
        annealSteps: 0,
      });
      annealLeft = m.annealSteps;
      msdPts = []; a2Pts = []; ghosts = []; history = [];
      nextSampleStep = 8; lastRecStep = 0; refT = sim.T;
      const sigmaCopy = sim.sigma.slice();
      const bondCopy = sim.bondPairs.slice();
      post({
        type: 'ready',
        N: sim.N, Lx: sim.Lx, Ly: sim.Ly, Lz: sim.Lz,
        sigma: sigmaCopy.buffer,
        bondPairs: bondCopy.buffer,
      }, [sigmaCopy.buffer, bondCopy.buffer]);
      break;
    }
  }
};

schedule();
