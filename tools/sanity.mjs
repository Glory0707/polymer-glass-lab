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
}

const bins = binByT(history, 0.1, 0, 1.6, 1);
const fit = twoSegmentFit(bins);
if (fit) {
  check('两段式拟合出拐点', fit.Tg > 0.2 && fit.Tg < 0.7,
    `Tg=${fit.Tg.toFixed(2)} (玻璃段斜率 ${fit.slopeGlass.toFixed(2)}, 液体段 ${fit.slopeLiquid.toFixed(2)})`);
} else {
  check('两段式拟合出拐点', false, '拐点不显著，请复核扫描数据');
}

const failed = results.filter((r) => !r.ok);
console.log(`\n${results.length - failed.length}/${results.length} 项通过`);
process.exit(failed.length ? 1 : 0);
