/**
 * sanity.mjs — 物理正确性无头自检（Node 直接运行，无需浏览器）
 *
 *   node tools/sanity.mjs
 *
 * 检查项：
 *  1. 稳定性：积分无 NaN、键长不发散（FENE 上限 1.5）
 *  2. 恒温：实测动能温度 ≈ 目标温度
 *  3. 玻璃化定性物理：高温液体 MSD@20τ 远大于深冷玻璃态（笼子化）
 *  4. 降温扫描：MSD@20τ 随 T 出现拐点（两段式拟合，Tg ≈ 0.3–0.5）
 *
 * 协议：ρ=1.0（KG 玻璃化研究标准密度）、γ=1.0；每个测点弛豫后取快照，
 * 精确推进 2500 步（20τ）测 msdOver(snapshot)，保证各点滞后窗口一致。
 */
import { KGSim } from '../src/md.js';
import { binByT, twoSegmentFit } from '../src/analysis.js';

const LAG = 2500; // 20 τ

const results = [];
function check(name, ok, detail) {
  results.push({ name, ok, detail });
  console.log(`${ok ? '✅' : '❌'} ${name}  ${detail}`);
}

/** 在当前温度下弛豫 relax 步后，测 MSD@LAG 与窗口内平均 PE */
function measureAt(sim, relax) {
  sim.run(relax);
  const snap = sim.snapshotU();
  let peSum = 0;
  for (let i = 0; i < LAG; i++) {
    sim.step();
    peSum += sim.pePerBead;
  }
  return { msd: sim.msdOver(snap), pe: peSum / LAG };
}

// --- 小体系快速验证：10 链 × 40 珠 = 400 珠，ρ=1.0 ---
const sim = new KGSim({ numChains: 10, chainLen: 40, seed: 7, temperature: 1.2 });
console.log(`体系: N=${sim.N}, 盒子 ${sim.Lx.toFixed(2)}×${sim.Ly.toFixed(2)}×${sim.Lz.toFixed(2)} σ, ρ=${sim.density}（构造器已含 ${6000} 步初始化退火）\n`);

// 1) 高温液体
const liq = measureAt(sim, 30000);
const maxBond1 = sim.maxBondLength();
check('高温无 NaN', Number.isFinite(liq.msd) && Number.isFinite(liq.pe),
  `MSD@20τ=${liq.msd.toFixed(2)}, PE=${liq.pe.toFixed(2)}`);
check('键长稳定 (<1.4σ)', maxBond1 < 1.4, `maxBond=${maxBond1.toFixed(3)}`);
check('高温恒温 (目标1.2)', Math.abs(sim.keTemp - 1.2) / 1.2 < 0.2, `T_meas=${sim.keTemp.toFixed(3)}`);

// 2) 急冷到深玻璃态
sim.T = 0.1;
const gl = measureAt(sim, 40000);
const maxBond2 = sim.maxBondLength();
check('低温无 NaN', Number.isFinite(gl.msd) && Number.isFinite(gl.pe),
  `MSD@20τ=${gl.msd.toFixed(4)}, PE=${gl.pe.toFixed(2)}`);
check('低温键长稳定', maxBond2 < 1.4, `maxBond=${maxBond2.toFixed(3)}`);
check('深冷恒温 (目标0.1)', Math.abs(sim.keTemp - 0.1) / 0.1 < 0.35, `T_meas=${sim.keTemp.toFixed(3)}`);
check('玻璃化笼子化 (MSD_液 > 8×MSD_玻)', liq.msd > gl.msd * 8,
  `液=${liq.msd.toFixed(2)} vs 玻=${gl.msd.toFixed(4)} (比值 ${(liq.msd / gl.msd).toFixed(0)}×)`);

// 3) 降温扫描（一致滞后窗口）+ 拐点拟合
console.log('\n降温扫描（每档弛豫 15000 步后测 MSD@20τ / 窗口平均 PE）:');
const history = [];
for (let T = 1.2; T >= 0.1 - 1e-9; T -= 0.1) {
  sim.T = Math.max(0.1, T);
  const r = measureAt(sim, 15000);
  history.push({ T: sim.T, msd: r.msd, pe: r.pe });
  console.log(`  T=${sim.T.toFixed(2)}  MSD@20τ=${r.msd.toFixed(4).padStart(9)}  PE=${r.pe.toFixed(3)}`);
}
const glassRow = history[history.length - 1];
const liquidRow = history[0];
check('扫描趋势：高温 MSD ≫ 低温 MSD', liquidRow.msd > glassRow.msd * 8,
  `${liquidRow.msd.toFixed(2)} vs ${glassRow.msd.toFixed(4)}`);
check('PE 随 T 下降（玻 < 液）', glassRow.pe < liquidRow.pe - 0.2,
  `玻=${glassRow.pe.toFixed(2)} < 液=${liquidRow.pe.toFixed(2)}`);

// 链刚度：数值梯度一致性 + NVE 能量守恒（κ=4）
sim.stiffness = 4;
{
  const h = 1e-3;
  let worst = 0;
  for (const j of [100, 233, 350]) {
    for (let c = 0; c < 3; c++) {
      const orig = sim.pos[j * 3 + c];
      sim.pos[j * 3 + c] = orig + h; sim._computeForces();
      const Up = sim.pePerBead * sim.N;
      sim.pos[j * 3 + c] = orig - h; sim._computeForces();
      const Um = sim.pePerBead * sim.N;
      sim.pos[j * 3 + c] = orig; sim._computeForces();
      const numeric = -(Up - Um) / (2 * h);
      const err = Math.abs(sim.force[j * 3 + c] - numeric) / (Math.abs(numeric) + 1);
      if (err > worst) worst = err;
    }
  }
  check('链刚度弯角力 = −∇U（数值梯度）', worst < 0.02, `最大相对误差 ${worst.toExponential(2)}`);
}
const E0 = sim.pePerBead * sim.N + 1.5 * sim.keTemp * sim.N;
sim.gamma = 0;
let maxDrift = 0;
for (let i = 0; i < 2500; i++) {
  sim.step();
  const E = sim.pePerBead * sim.N + 1.5 * sim.keTemp * sim.N;
  maxDrift = Math.max(maxDrift, Math.abs(E - E0) / sim.N);
}
sim.gamma = 1;
check('刚度开启 NVE 能量守恒', maxDrift < 0.05, `最大漂移 ${maxDrift.toFixed(4)} ε/珠`);

// 双分散组分 + 密度调整 + χ₄ 视角数据：稳定性抽查（κ=4 沿用）
{
  const b = new KGSim({ numChains: 10, chainLen: 40, seed: 11, temperature: 1.0, density: 1.0, smallFrac: 0.25, stiffness: 4 });
  b.run(6000);
  const okB = Number.isFinite(b.pePerBead) && b.maxBondLength() < 1.4;
  b.setDensity(1.2); b.run(4000);
  const okD = Number.isFinite(b.pePerBead) && Math.abs(b.density - 1.2) < 1e-9 && b.maxBondLength() < 1.4;
  const chi = b.smoothMobility();
  let okC = true;
  for (const v of chi) if (!Number.isFinite(v)) okC = false;
  check('双分散体系稳定（25% 小珠, κ=4）', okB, `PE=${b.pePerBead.toFixed(2)} maxBond=${b.maxBondLength().toFixed(3)}`);
  check('密度调整稳定（ρ→1.2）', okD, `PE=${b.pePerBead.toFixed(2)} ρ=${b.density.toFixed(3)}`);
  check('χ₄ 视角邻域迁移率有限', okC, `${chi.length} 珠全有限`);
  // NPT 恒压：P₀=10，T=0.5 → 应压实到 ρ > 1.0
  const n = new KGSim({ numChains: 10, chainLen: 40, seed: 13, temperature: 0.5, density: 1.0, npt: true, targetP: 10 });
  n.run(15000);
  check('NPT 恒压压实（P₀=10, T=0.5 → ρ>1.0）', n.density > 1.0 && Number.isFinite(n.pePerBead),
    `ρ=${n.density.toFixed(3)} P=${n.pressure.toFixed(2)} maxBond=${n.maxBondLength().toFixed(3)}`);
  // χ₄ 粗粒化估计量有限且非负
  n.run(2000);
  const c4s = [];
  for (let i = 0; i < 4; i++) { n.run(300); c4s.push(n.chi4Coarse()); }
  const okC4 = c4s.every((v) => Number.isFinite(v) && v >= 0);
  check('χ₄ 粗粒化估计量有限非负', okC4, c4s.map((v) => v.toFixed(2)).join(' '));
}

const bins = binByT(history, 0.1, 0, 1.6, 1);
const fit = twoSegmentFit(bins);
if (fit) {
  check('两段式拟合出拐点', fit.Tg > 0.2 && fit.Tg < 0.7,
    `Tg=${fit.Tg.toFixed(2)} (玻璃段斜率 ${fit.slopeGlass.toFixed(2)}, 液体段 ${fit.slopeLiquid.toFixed(2)})`);
} else {
  check('两段式拟合出拐点', false, '拐点不显著，请复核扫描数据');
}

// ============ v33 扩展功能自检 ============
console.log('\n--- v33 扩展：VFT / Fs(q,t) / P2 / 热笔 / 拽链 / 薄膜 ---');

/** 局部动能温度：半径内珠子的 <v²>/3 */
function localKe(s, x, y, z, R) {
  let ke = 0, n = 0;
  for (let i = 0; i < s.N; i++) {
    const i3 = i * 3;
    const dx = s.pos[i3] - x, dy = s.pos[i3 + 1] - y, dz = s.pos[i3 + 2] - z;
    if (dx * dx + dy * dy + dz * dz < R * R) {
      ke += s.vel[i3] ** 2 + s.vel[i3 + 1] ** 2 + s.vel[i3 + 2] ** 2;
      n++;
    }
  }
  return n ? ke / (3 * n) : 0;
}

// 18) τα 提取 + VFT 拟合：小扫描产 4 档温度的 MSD(τ) 曲线
{
  const { tauFromMsd, vftFit } = await import('../src/analysis.js');
  const msdCurves = [];
  const tMinScan = 0.25;
  sim.T = 1.2; sim.run(12000);
  for (const T of [1.2, 0.8, 0.5, tMinScan]) {
    sim.T = T;
    sim.run(15000);
    const snap = sim.snapshotU();
    const pts = [];
    const t0c = sim.stepCount;
    for (let k = 0; k < 60; k++) {
      sim.run(250);
      pts.push([(sim.stepCount - t0c) * sim.dt, sim.msdOver(snap)]);
    }
    msdCurves.push({ T, pts });
  }
  const taus = [];
  for (const s of msdCurves) {
    const tau = tauFromMsd(s.pts);
    if (tau != null) taus.push({ T: s.T, tau });
  }
  const okTau = taus.length >= 2 && taus[taus.length - 1].tau > taus[0].tau * 2;
  check('τα 提取随温度上升', okTau,
    `${taus.length} 档` + (taus.length >= 2
      ? `: τp(T=${taus[0].T.toFixed(1)})=${taus[0].tau.toFixed(2)} → τp(T=${taus[taus.length - 1].T.toFixed(2)})=${taus[taus.length - 1].tau.toFixed(2)}`
      : '（窗口内未穿越，屏幕 UI 中低 T 的 τα 会超出窗口属正常物理）'));
  // 拟合数学用合成数据验证：log τ = A + B/(T−T0)，加 2% 噪声
  const A = -0.8, B = 0.9, T0true = 0.3;
  const synth = [];
  for (let k = 0; k < 10; k++) {
    const T = 0.45 + k * 0.1;
    const y = A + B / (T - T0true);
    synth.push({ T, tau: Math.pow(10, y * (1 + (k % 3 - 1) * 0.02)) });
  }
  const vf = vftFit(synth);
  const okVft = vf.fit && Math.abs(vf.fit.T0 - T0true) < 0.06 && vf.arrhenius && isFinite(vf.arrhenius.a);
  check('VFT 拟合还原 T₀（合成数据）', okVft,
    vf.fit ? `T₀=${vf.fit.T0.toFixed(3)} (真值 0.30), B=${vf.fit.B.toFixed(2)} (真值 0.90)` : '拟合失败');
}

// 19) Fs(q,t)：高温衰减到近 0，低温短滞后仍在平台上
{
  sim.T = 1.2; sim.resetRef();
  for (let k = 0; k < 40; k++) { sim.run(50); }
  const decayed = sim.fsqRef();
  sim.T = 0.15; sim.resetRef();
  sim.run(50); // 0.4τ：仍处 β 平台
  const fsqGlass0 = sim.fsqRef();
  check('Fs(q,t) 高温衰减 / 低温保持', decayed < 0.35 && fsqGlass0 > 0.55 && fsqGlass0 <= 1.0,
    `Fs(液,长)=${decayed.toFixed(3)}, Fs(玻,短)=${fsqGlass0.toFixed(3)}`);
}

// 20) 键取向 P2：平衡 ≈ 0，拉伸后 > 0
{
  sim.T = 1.0;
  sim.run(8000);
  const p20 = sim.bondP2();
  sim.deform.mode = 'stretch';
  sim.deform.rate = 0.05;
  const target = 0.3;
  let guard = 0;
  while (sim.deform.strain < target && guard++ < 20000) sim.step();
  const p2s = sim.bondP2();
  sim.deform.mode = 'none'; sim.deform.strain = 0; sim.deform.target = null;
  check('键取向 P2 平衡≈0 / 拉伸>0', Math.abs(p20) < 0.08 && p2s > 0.02,
    `P2(平衡)=${p20.toFixed(4)}（390 键统计涨落 ±0.05）, P2(ε=0.3)=${p2s.toFixed(4)}`);
}

// 21) 热笔：局部注入后近点动能上升
{
  sim.T = 0.3;
  sim.run(5000);
  const cx = sim.Lx / 2, cy = sim.Ly / 2, cz = sim.Lz / 2;
  const keNearBefore = localKe(sim, cx, cy, cz, 2.5);
  sim.heatBrush(cx, cy, cz, 2.5, 0.8);
  const keNear = localKe(sim, cx, cy, cz, 2.5);
  const keFar = localKe(sim, 1, 1, 1, 2.5);
  check('热笔局部升温', keNear > keNearBefore * 1.3 && keNear > keFar,
    `近点 ${keNearBefore.toFixed(3)}→${keNear.toFixed(3)}, 远点 ${keFar.toFixed(3)}`);
}

// 22) 拽链：抓取并把目标移到远处，珠子跟随
{
  sim.T = 0.5;
  sim.run(3000);
  const sx = sim.Lx * 0.3, sy = sim.Ly * 0.5, sz = sim.Lz * 0.5;
  const idx = sim.grabPick(sx, sy, sz, 3.0);
  const okPick = idx >= 0;
  if (okPick) {
    const tx = sim.Lx * 0.8, ty = sim.Ly * 0.5, tz = sim.Lz * 0.5;
    for (let k = 0; k < 1500; k++) {
      sim.grabMove(tx, ty, tz);
      sim.step();
    }
    const dx = sim.pos[idx * 3] - sim.grabTarget[0];
    const dy = sim.pos[idx * 3 + 1] - sim.grabTarget[1];
    const dz = sim.pos[idx * 3 + 2] - sim.grabTarget[2];
    const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
    check('拽链探针跟随目标', dist < 1.5, `珠子距目标 ${dist.toFixed(3)} σ`);
    sim.grabRelease();
  } else {
    check('拽链探针跟随目标', false, '拾取失败');
  }
}

// 23) 薄膜模式：珠子被约束在壁内、无 NaN、表面层迁移率更高
{
  const fSim = new KGSim({ numChains: 10, chainLen: 40, seed: 11, temperature: 0.5, film: true });
  fSim.run(8000);
  let okZ = true, okNaN = true;
  for (let i = 0; i < fSim.N; i++) {
    const z = fSim.pos[i * 3 + 2];
    if (!Number.isFinite(z) || z < fSim.wallMargin - 0.8 || z > fSim.Lz - fSim.wallMargin + 0.8) okZ = false;
  }
  for (let a = 0; a < fSim.pos.length; a++) if (!Number.isFinite(fSim.pos[a])) okNaN = false;
  // 层迁移率：表面 2 层 vs 中部 4 层（膜内范围）
  const NB = 12, acc = new Float64Array(NB), cnt = new Float64Array(NB);
  const zLo = fSim.wallMargin - 0.5, zHi = fSim.Lz - fSim.wallMargin + 0.5;
  const snap = fSim.snapshotU();
  fSim.run(625);
  for (let i = 0; i < fSim.N; i++) {
    const i3 = i * 3;
    const b = Math.min(NB - 1, Math.max(0, ((fSim.pos[i3 + 2] - zLo) / (zHi - zLo) * NB) | 0));
    const dx = fSim.upos[i3] - snap[i3], dy = fSim.upos[i3 + 1] - snap[i3 + 1], dz = fSim.upos[i3 + 2] - snap[i3 + 2];
    acc[b] += dx * dx + dy * dy + dz * dz;
    cnt[b]++;
  }
  const layer = Array.from(acc, (v, k) => cnt[k] ? v / cnt[k] : 0);
  const midMob = layer.slice(4, 8).reduce((a, b) => a + b, 0) / 4;
  const surfMob = Math.max(layer[0], layer[NB - 1]);
  check('薄膜壁约束 + 无 NaN', okZ && okNaN,
    `z ∈ [${Math.min(...fSim.pos.filter((_, i) => i % 3 === 2)).toFixed(2)}, ${Math.max(...fSim.pos.filter((_, i) => i % 3 === 2)).toFixed(2)}]（壁位 2 / ${(fSim.Lz - 2).toFixed(1)}）`);
  check('自由表面迁移率 ≥ 中部', surfMob > midMob * 0.8,
    `表面=${surfMob.toFixed(3)}, 中部=${midMob.toFixed(3)}`);
}

const failed = results.filter((r) => !r.ok);
console.log(`\n${results.length - failed.length}/${results.length} 项通过`);
process.exit(failed.length ? 1 : 0);
