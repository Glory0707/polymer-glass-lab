/**
 * md.js — Kremer–Grest 珠簧模型分子动力学核心
 * 纯计算模块，不依赖 DOM，可在浏览器与 Node 中复用（见 tools/sanity.mjs）。
 *
 * 模型（LJ 约化单位 σ = ε = m = k_B = 1）：
 *   - 非键合珠子：WCA 排斥势（LJ 在 r_c = 2^(1/6) 截断并平移至零，纯排斥）
 *   - 链上相邻珠子：FENE 键（k = 30, R0 = 1.5）
 *     Grest & Kremer, J. Chem. Phys. 86, 4068 (1986)
 *     Kremer & Grest, J. Chem. Phys. 92, 5057 (1990)
 * 积分：velocity-Verlet + Langevin 恒温器（每整步一次 Ornstein–Uhlenbeck 速度更新）
 * 边界：三方向周期性边界 + 最小镜像；元胞列表邻居搜索
 */

export const RC = Math.pow(2, 1 / 6); // WCA 截断距离 ≈ 1.12246
const RC2 = RC * RC;
const K_FENE = 30;
const R0 = 1.5;
const R02 = R0 * R0;

function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/** 把珠子数 n 分解为最接近立方体的 nx*ny*nz */
export function latticeDims(n) {
  let best = null, bestScore = Infinity;
  const lim = Math.ceil(Math.cbrt(n)) + 1;
  for (let nx = 1; nx <= lim; nx++) {
    if (n % nx) continue;
    const m = n / nx;
    for (let ny = nx; ny * ny <= m; ny++) {
      if (m % ny) continue;
      const nz = m / ny;
      const s = [nx, ny, nz].sort((p, q) => p - q);
      const score = s[2] / s[0];
      if (score < bestScore) { bestScore = score; best = s; }
    }
  }
  return { nx: best[0], ny: best[1], nz: best[2] };
}

export class KGSim {
  constructor(opts = {}) {
    this.chainLen = opts.chainLen ?? 40;
    this.numChains = opts.numChains ?? 60;
    // ρ=1.0 是 KG 玻璃化研究的标准密度（有效堆积分数使 Tg 落在 ~0.4）；
    // ρ=0.85 的熔体要降到 ~0.1 才玻璃化
    this.density = opts.density ?? 1.0;
    this.dt = opts.dt ?? 0.008;
    this.gamma = opts.gamma ?? 1.0;
    // 随机游走初始化会有同格点堆叠，需要在 dt/4 下小步退火化解，
    // 否则大力 + 大步长会让珠子隧穿势墙正反馈爆掉
    this.annealSteps = opts.annealSteps ?? 6000;
    // MSD@τ_lag 观测窗口（默认 2500 步 ≈ 20 τ，画“热历史 vs T”用）
    this.lagSteps = opts.lagSteps ?? 2500;
    // 链段迁移率着色窗口（默认 625 步 ≈ 5 τ）
    this.mobSteps = opts.mobSteps ?? 625;
    this.T = opts.temperature ?? 1.0;
    // 链刚度 κ：U_bend = κ(1 - cosφ)，φ 为相邻两键夹角（0 = 柔性，越大越挺直）
    this.stiffness = opts.stiffness ?? 0;
    // 双分散：小珠（σ = 0.7σ₀）占比，增塑效应使 Tg 下移
    this.smallFrac = opts.smallFrac ?? 0;
    // 恒压（NPT）：Berendsen 各向同性弱耦合，P₀ = 0 时密度自调
    this.npt = opts.npt ?? false;
    this.targetP = opts.targetP ?? 0;
    this.tauP = 1.0;
    this.virial = 0;
    // 维里张量分量（应力用）与形变状态（力学轴）
    this.wxx = 0; this.wyy = 0; this.wzz = 0; this.wxy = 0;
    this.deform = { mode: 'none', rate: 0.02, amp: 0.12, freq: 0.5, strain: 0, epsCur: 0, phase: 0, gamma: 0, target: null };
    this.seed = (opts.seed ?? 20260910) >>> 0;

    this.rng = mulberry32(this.seed);
    this._gaussSaved = null;
    this.stepCount = 0;
    this.lagAnchor = 0;
    this.mobAnchor = 0;
    this.refStep = 0;
    this.refT = this.T;
    this.pePerBead = 0;
    this.keTemp = this.T;

    this._initChains();
    this._allocCells();

    const N3 = this.N * 3;
    this.refPos = new Float64Array(N3); // MSD 参考构象（非折叠坐标）
    this.snapLag = new Float64Array(N3);
    this.snapMob = new Float64Array(N3);

    this._buildCells();
    this._computeForces();
    if (this.annealSteps > 0) this.anneal(this.annealSteps);
    this.resetRef();
    this.snapLag.set(this.upos);
    this.snapMob.set(this.upos);
  }

  get N() { return this.numChains * this.chainLen; }
  get time() { return this.stepCount * this.dt; }

  /**
   * 初始构象：晶格上的非回溯随机游走链（优先空位），加小扰动。
   * 注意不能用拉直的"蛇形棒"初始化——巨末端距的键张力会在低温下持续驱动
   * 链段重排，制造不冻结的 MSD 底噪（标准 KG 协议用无规线团起始）。
   */
  _initChains() {
    const N = this.N;
    const L = this.chainLen;
    const { nx, ny, nz } = latticeDims(N);
    const a = Math.pow(1 / this.density, 1 / 3);
    this.Lx = nx * a; this.Ly = ny * a; this.Lz = nz * a;

    this.pos = new Float32Array(N * 3);   // 折叠坐标（用于力计算/渲染）
    this.upos = new Float64Array(N * 3);  // 非折叠坐标（用于 MSD）
    this.vel = new Float32Array(N * 3);
    this.force = new Float32Array(N * 3);

    const siteIdx = (x, y, z) =>
      (((x % nx) + nx) % nx) + nx * ((((y % ny) + ny) % ny) + ny * (((z % nz) + nz) % nz));
    const DX = [1, -1, 0, 0, 0, 0];
    const DY = [0, 0, 1, -1, 0, 0];
    const DZ = [0, 0, 0, 0, 1, -1];
    const occ = new Uint8Array(N);

    for (let c = 0; c < this.numChains; c++) {
      // 随机挑一个空位点作起点
      let s = (this.rng() * N) | 0;
      for (let guard = 0; occ[s] && guard < 20 * N; guard++) s = (this.rng() * N) | 0;
      let x = s % nx;
      let y = ((s / nx) | 0) % ny;
      let z = (s / (nx * ny)) | 0;
      occ[s] = 1;
      let last = -1;

      for (let k = 0; k < L; k++) {
        if (k > 0) {
          // 非回溯随机游走，优先走向未占用的格点
          let chosen = -1;
          for (let t = 0; t < 12 && chosen < 0; t++) {
            const d = (this.rng() * 6) | 0;
            if (d === (last ^ 1)) continue;
            if (!occ[siteIdx(x + DX[d], y + DY[d], z + DZ[d])]) chosen = d;
          }
          if (chosen < 0) {
            do { chosen = (this.rng() * 6) | 0; } while (chosen === (last ^ 1));
          }
          x += DX[chosen]; y += DY[chosen]; z += DZ[chosen];
          last = chosen;
          occ[siteIdx(x, y, z)] = 1;
        }
        const i3 = (c * L + k) * 3;
        this.pos[i3] = (x + 0.5) * a + (this.rng() - 0.5) * 0.2;
        this.pos[i3 + 1] = (y + 0.5) * a + (this.rng() - 0.5) * 0.2;
        this.pos[i3 + 2] = (z + 0.5) * a + (this.rng() - 0.5) * 0.2;
        this.upos[i3] = this.pos[i3];
        this.upos[i3 + 1] = this.pos[i3 + 1];
        this.upos[i3 + 2] = this.pos[i3 + 2];
      }
    }

    // 逐珠尺寸参数 σ（双分散时按 smallFrac 随机指派小珠）
    this.sigma = new Float32Array(N);
    for (let i = 0; i < N; i++) {
      this.sigma[i] = this.rng() < this.smallFrac ? 0.7 : 1.0;
    }

    // 键表：链内相邻珠子
    const nb = N - this.numChains;
    this.bondPairs = new Int32Array(nb * 2);
    let b = 0;
    for (let i = 0; i < N; i++) {
      if (i % L === L - 1) continue;
      this.bondPairs[b * 2] = i;
      this.bondPairs[b * 2 + 1] = i + 1;
      b++;
    }
  }

  _allocCells() {
    // 每维元胞数取 floor(L/rc)，保证元胞边长 ≥ r_c，27 元胞搜索即完备
    this.ncx = Math.max(1, Math.floor(this.Lx / RC));
    this.ncy = Math.max(1, Math.floor(this.Ly / RC));
    this.ncz = Math.max(1, Math.floor(this.Lz / RC));
    const total = this.ncx * this.ncy * this.ncz;
    if (!this.cellHead || this.cellHead.length !== total) {
      this.cellHead = new Int32Array(total);
    }
    this.cellNext = new Int32Array(this.N);
  }

  _buildCells() {
    const head = this.cellHead;
    head.fill(-1);
    const { Lx, Ly, Lz, ncx, ncy, ncz } = this;
    const p = this.pos;
    for (let i = 0; i < this.N; i++) {
      const i3 = i * 3;
      let cx = (p[i3] / Lx * ncx) | 0;
      let cy = (p[i3 + 1] / Ly * ncy) | 0;
      let cz = (p[i3 + 2] / Lz * ncz) | 0;
      if (cx >= ncx) cx = ncx - 1; if (cx < 0) cx = 0;
      if (cy >= ncy) cy = ncy - 1; if (cy < 0) cy = 0;
      if (cz >= ncz) cz = ncz - 1; if (cz < 0) cz = 0;
      const c = cx + ncx * (cy + ncy * cz);
      this.cellNext[i] = head[c];
      head[c] = i;
    }
  }

  _computeForces() {
    const f = this.force;
    f.fill(0);
    this.wxx = 0; this.wyy = 0; this.wzz = 0; this.wxy = 0;
    const shearOff = this.deform.gamma * this.Ly; // y 镜像行的 x 偏移 = γ·Ly
    const p = this.pos;
    const { Lx, Ly, Lz, ncx, ncy, ncz, N } = this;
    let pe = 0;

    for (let i = 0; i < N; i++) {
      const i3 = i * 3;
      const xi = p[i3], yi = p[i3 + 1], zi = p[i3 + 2];
      let cx = (xi / Lx * ncx) | 0;
      let cy = (yi / Ly * ncy) | 0;
      let cz = (zi / Lz * ncz) | 0;
      if (cx >= ncx) cx = ncx - 1; if (cx < 0) cx = 0;
      if (cy >= ncy) cy = ncy - 1; if (cy < 0) cy = 0;
      if (cz >= ncz) cz = ncz - 1; if (cz < 0) cz = 0;

      for (let oz = -1; oz <= 1; oz++) {
        let z2 = cz + oz; if (z2 < 0) z2 += ncz; else if (z2 >= ncz) z2 -= ncz;
        for (let oy = -1; oy <= 1; oy++) {
          let y2 = cy + oy; if (y2 < 0) y2 += ncy; else if (y2 >= ncy) y2 -= ncy;
          for (let ox = -1; ox <= 1; ox++) {
            let x2 = cx + ox; if (x2 < 0) x2 += ncx; else if (x2 >= ncx) x2 -= ncx;
            let j = this.cellHead[x2 + ncx * (y2 + ncy * z2)];
            while (j !== -1) {
              if (j > i) {
                const j3 = j * 3;
                let dx = xi - p[j3];
                let dy = yi - p[j3 + 1];
                let dz = zi - p[j3 + 2];
                if (shearOff !== 0) dx -= shearOff * Math.round(dy / Ly);
                dx -= Lx * Math.round(dx / Lx);
                dy -= Ly * Math.round(dy / Ly);
                dz -= Lz * Math.round(dz / Lz);
                const r2 = dx * dx + dy * dy + dz * dz;
                const s2ij = (this.sigma[i] + this.sigma[j]) * 0.5;
                const rc2ij = 1.2599210498948732 * s2ij * s2ij; // (2^(1/6)·σij)²
                if (r2 < rc2ij) {
                  const inv2 = 1 / r2;
                  const inv6 = inv2 * inv2 * inv2;
                  const s6 = s2ij * s2ij * s2ij;
                  const i6 = s6 * inv6;
                  const i12 = i6 * i6;
                  let fc = 24 * (2 * i12 - i6) * inv2; // WCA
                  // 数值护栏：仅在初始化瞬间的深重叠时触发（正常动力学 fc ≲ 100）
                  if (fc > 2e3) fc = 2e3;
                  const fxc = fc * dx, fyc = fc * dy, fzc = fc * dz;
                  f[i3] += fxc; f[i3 + 1] += fyc; f[i3 + 2] += fzc;
                  f[j3] -= fxc; f[j3 + 1] -= fyc; f[j3 + 2] -= fzc;
                  let uw = 4 * (i12 - i6) + 1; // 已平移至截断处为零
                  if (uw > 2e3) uw = 2e3;
                  pe += uw;
                  this.wxx += fxc * dx; this.wxy += fxc * dy;
                  this.wyy += fyc * dy; this.wzz += fzc * dz;
                }
              }
              j = this.cellNext[j];
            }
          }
        }
      }
    }

    // FENE 键
    const L = this.chainLen;
    for (let i = 0; i < N; i++) {
      if (i % L === L - 1) continue;
      const i3 = i * 3, j3 = i3 + 3;
      let dx = p[i3] - p[j3];
      let dy = p[i3 + 1] - p[j3 + 1];
      let dz = p[i3 + 2] - p[j3 + 2];
      if (shearOff !== 0) dx -= shearOff * Math.round(dy / Ly);
      dx -= Lx * Math.round(dx / Lx);
      dy -= Ly * Math.round(dy / Ly);
      dz -= Lz * Math.round(dz / Lz);
      let r2 = dx * dx + dy * dy + dz * dz;
      if (r2 > 0.98 * R02) r2 = 0.98 * R02; // 数值护栏：log 参数恒为正
      const denom = 1 - r2 / R02;
      const fc = -K_FENE / denom;
      const fx = fc * dx, fy = fc * dy, fz = fc * dz;
      f[i3] += fx; f[i3 + 1] += fy; f[i3 + 2] += fz;
      f[j3] -= fx; f[j3 + 1] -= fy; f[j3 + 2] -= fz;
      pe += -0.5 * K_FENE * R02 * Math.log(denom);
      this.wxx += fx * dx; this.wxy += fx * dy;
      this.wyy += fy * dy; this.wzz += fz * dz;
    }

    // 弯角势：相邻两键夹角的 cos 型弯曲能（半柔性链）
    if (this.stiffness > 0) {
      const k = this.stiffness;
      for (let i = 0; i < N; i++) {
        if (i % L === 0 || i % L === L - 1) continue; // 链端无弯角
        const a3 = (i - 1) * 3, i3 = i * 3, b3 = (i + 1) * 3;
        let b1x = p[i3] - p[a3];     b1x -= Lx * Math.round(b1x / Lx);
        let b1y = p[i3 + 1] - p[a3 + 1]; b1y -= Ly * Math.round(b1y / Ly);
        let b1z = p[i3 + 2] - p[a3 + 2]; b1z -= Lz * Math.round(b1z / Lz);
        let b2x = p[b3] - p[i3];
        let b2y = p[b3 + 1] - p[i3 + 1];
        let b2z = p[b3 + 2] - p[i3 + 2];
        if (shearOff !== 0) b2x -= shearOff * Math.round(b2y / Ly);
        b2x -= Lx * Math.round(b2x / Lx);
        b2y -= Ly * Math.round(b2y / Ly);
        b2z -= Lz * Math.round(b2z / Lz);
        const inv1 = 1 / Math.sqrt(b1x * b1x + b1y * b1y + b1z * b1z);
        const inv2 = 1 / Math.sqrt(b2x * b2x + b2y * b2y + b2z * b2z);
        const n1x = b1x * inv1, n1y = b1y * inv1, n1z = b1z * inv1;
        const n2x = b2x * inv2, n2y = b2y * inv2, n2z = b2z * inv2;
        let c = n1x * n2x + n1y * n2y + n1z * n2z;
        if (c > 1) c = 1; else if (c < -1) c = -1;
        // dU/db = -κ · ∂cosφ/∂b，∇cosφ 见下
        const g1x = (n2x - c * n1x) * inv1, g1y = (n2y - c * n1y) * inv1, g1z = (n2z - c * n1z) * inv1;
        const g2x = (n1x - c * n2x) * inv2, g2y = (n1y - c * n2y) * inv2, g2z = (n1z - c * n2z) * inv2;
        f[a3] += -k * g1x; f[a3 + 1] += -k * g1y; f[a3 + 2] += -k * g1z;
        f[i3] += k * (g1x - g2x); f[i3 + 1] += k * (g1y - g2y); f[i3 + 2] += k * (g1z - g2z);
        f[b3] += k * g2x; f[b3 + 1] += k * g2y; f[b3 + 2] += k * g2z;
        pe += k * (1 - c);
      }
    }

    this.pePerBead = pe / N;
  }

  _gauss() {
    if (this._gaussSaved !== null) {
      const g = this._gaussSaved;
      this._gaussSaved = null;
      return g;
    }
    let u1 = 0;
    do { u1 = this.rng(); } while (u1 === 0);
    const u2 = this.rng();
    const r = Math.sqrt(-2 * Math.log(u1));
    const th = 2 * Math.PI * u2;
    this._gaussSaved = r * Math.sin(th);
    return r * Math.cos(th);
  }

  /** 推进一个时间步：半踢 → Langevin OU → 漂移 → 重算力 → 半踢 */
  step() {
    const dt = this.dt;
    const half = dt * 0.5;
    const damp = Math.max(0, 1 - this.gamma * dt);
    const noise = Math.sqrt(2 * this.gamma * this.T * dt);
    const p = this.pos, v = this.vel, f = this.force, u = this.upos;

    let ke2 = 0;
    for (let a = 0; a < v.length; a++) {
      const nv = v[a] * damp + noise * this._gauss() + f[a] * half;
      v[a] = nv;
      ke2 += nv * nv;
    }
    this.keTemp = ke2 / (3 * this.N); // k_B=1, m=1 → T = <v²>/3

    for (let i = 0; i < this.N; i++) {
      const i3 = i * 3;
      u[i3] += v[i3] * dt;
      u[i3 + 1] += v[i3 + 1] * dt;
      u[i3 + 2] += v[i3 + 2] * dt;
      let np = p[i3] + v[i3] * dt;
      if (np < 0) np += this.Lx; else if (np >= this.Lx) np -= this.Lx;
      p[i3] = np;
      np = p[i3 + 1] + v[i3 + 1] * dt;
      if (np < 0) np += this.Ly; else if (np >= this.Ly) np -= this.Ly;
      p[i3 + 1] = np;
      np = p[i3 + 2] + v[i3 + 2] * dt;
      if (np < 0) np += this.Lz; else if (np >= this.Lz) np -= this.Lz;
      p[i3 + 2] = np;
    }

    this._deformStep();
    this._buildCells();
    this._computeForces();

    for (let a = 0; a < v.length; a++) v[a] += f[a] * half;

    this.stepCount++;
    if (this.stepCount - this.lagAnchor >= this.lagSteps) {
      this.snapLag.set(this.upos);
      this.lagAnchor = this.stepCount;
    }
    if (this.stepCount - this.mobAnchor >= this.mobSteps) {
      this.snapMob.set(this.upos);
      this.mobAnchor = this.stepCount;
    }
    // NPT：Berendsen 各向同性弱耦合（每个整步一次）
    if (this.npt) {
      const beta = 0.5, lambda = Math.max(0.98, Math.min(1.02,
        1 - beta * (this.dt / this.tauP) * (this.targetP - this.pressure)));
      this._scaleBox(lambda, lambda, lambda);
    }
  }

  run(nSteps) { for (let i = 0; i < nSteps; i++) this.step(); }

  /** 初始化退火：dt/4 小步长化解随机游走初始化的深重叠（浏览器可分帧调用） */
  anneal(nSteps) {
    const dt0 = this.dt;
    this.dt = dt0 / 4;
    this.run(nSteps);
    this.dt = dt0;
  }

  /** 重置 MSD 参考点（换温度时调用） */
  resetRef() {
    this.refPos.set(this.upos);
    this.refStep = this.stepCount;
    this.refT = this.T;
  }

  /** MSD 相对参考点（当前温度曲线） */
  msdRef() {
    let s = 0;
    const u = this.upos, r = this.refPos;
    for (let a = 0; a < u.length; a++) {
      const d = u[a] - r[a];
      s += d * d;
    }
    return s / this.N;
  }

  /** 固定滞后窗口 MSD（热历史图用），lagAge() 达到窗口时长后有效 */
  msdLag() {
    let s = 0;
    const u = this.upos, r = this.snapLag;
    for (let a = 0; a < u.length; a++) {
      const d = u[a] - r[a];
      s += d * d;
    }
    return s / this.N;
  }

  lagAge() { return (this.stepCount - this.lagAnchor) * this.dt; }
  mobAge() { return (this.stepCount - this.mobAnchor) * this.dt; }

  /**
   * 调整体系密度：按 f = (ρ_old/ρ_new)^(1/3) 等比缩放盒子与全部坐标。
   * 由调用方渐进调用（每帧一小步），避免密度突变冲击体系。
   */
  setDensity(rhoNew) {
    const f = Math.pow(this.density / rhoNew, 1 / 3);
    this._scaleBox(f, f, f);
  }

  /** 各向异性盒缩放：单轴/循环形变、密度调整与 NPT 共用 */
  _scaleBox(fx, fy, fz) {
    this.Lx *= fx; this.Ly *= fy; this.Lz *= fz;
    this.density = this.N / (this.Lx * this.Ly * this.Lz);
    for (let a = 0; a < this.pos.length; a += 3) {
      this.pos[a] *= fx; this.pos[a + 1] *= fy; this.pos[a + 2] *= fz;
      this.upos[a] *= fx; this.upos[a + 1] *= fy; this.upos[a + 2] *= fz;
    }
    this._allocCells();
  }

  /** 每步形变：单轴拉伸 / 循环 / LE 剪切（在漂移后、元胞重建前调用） */
  _deformStep() {
    const d = this.deform;
    if (d.mode === 'stretch') {
      if (d.target != null && d.strain >= d.target) {
        const back = d.target - d.strain;
        this._scaleBox(Math.exp(back), Math.exp(-back / 2), Math.exp(-back / 2));
        d.strain = d.target;
        d.mode = 'none';
        return;
      }
      const f = Math.exp(d.rate * this.dt);
      this._scaleBox(f, 1 / Math.sqrt(f), 1 / Math.sqrt(f));
      d.strain += d.rate * this.dt;
    } else if (d.mode === 'cyclic') {
      d.phase += 2 * Math.PI * d.freq * this.dt;
      const eps = d.amp * Math.sin(d.phase);
      this._scaleBox((1 + eps) / (1 + d.epsCur), Math.sqrt((1 + d.epsCur) / (1 + eps)), Math.sqrt((1 + d.epsCur) / (1 + eps)));
      d.epsCur = eps;
    }
  }

  /**
   * Steinhardt |q₆|（逐珠，1.4σ 邻域）：晶体序指标
   * 返回 { q6: Float32Array(N), mean: 全局平均 }
   */
  q6PerBead() {
    const N = this.N, p = this.pos;
    const rc = 1.4, rc2 = rc * rc;
    const C = new Float64Array(7 * N), S = new Float64Array(7 * N);
    const cnt = new Int32Array(N);
    this._buildCells();
    const Lx = this.Lx, Ly = this.Ly, Lz = this.Lz;
    // l=6 复球谐归一化系数（m ≥ 0）
    const NRM = new Float64Array(7);
    for (let m = 0; m <= 6; m++) {
      NRM[m] = Math.sqrt((13 / (4 * Math.PI)) * NFAC[6 - m] / NFAC[6 + m]);
    }
    for (let i = 0; i < N; i++) {
      const i3 = i * 3;
      const xi = p[i3], yi = p[i3 + 1], zi = p[i3 + 2];
      let cx = (xi / Lx * this.ncx) | 0, cy = (yi / this.Ly * this.ncy) | 0, cz = (zi / this.Lz * this.ncz) | 0;
      if (cx >= this.ncx) cx = this.ncx - 1;
      if (cy >= this.ncy) cy = this.ncy - 1;
      if (cz >= this.ncz) cz = this.ncz - 1;
      for (let oz = -1; oz <= 1; oz++) {
        let z2 = cz + oz; if (z2 < 0) z2 += this.ncz; else if (z2 >= this.ncz) z2 -= this.ncz;
        for (let oy = -1; oy <= 1; oy++) {
          let y2 = cy + oy; if (y2 < 0) y2 += this.ncy; else if (y2 >= this.ncy) y2 -= this.ncy;
          for (let ox = -1; ox <= 1; ox++) {
            let x2 = cx + ox; if (x2 < 0) x2 += this.ncx; else if (x2 >= this.ncx) x2 -= this.ncx;
            let j = this.cellHead[x2 + this.ncx * (y2 + this.ncy * z2)];
            while (j !== -1) {
              if (j > i) {
                const j3 = j * 3;
                let dx = xi - p[j3];     dx -= Lx * Math.round(dx / Lx);
                let dy = yi - p[j3 + 1]; dy -= Ly * Math.round(dy / Ly);
                let dz = zi - p[j3 + 2]; dz -= Lz * Math.round(dz / Lz);
                const r2 = dx * dx + dy * dy + dz * dz;
                if (r2 < rc2 && r2 > 1e-12) {
                  const inv = 1 / Math.sqrt(r2);
                  const ct = Math.min(1, Math.max(-1, dz * inv));
                  const st = Math.sqrt(Math.max(0, 1 - ct * ct));
                  const phi = Math.atan2(dy, dx);
                  for (let m = 0; m <= 6; m++) {
                    const y6c = NRM[m] * pl6m(m, ct) * Math.cos(m * phi);
                    const y6s = NRM[m] * pl6m(m, ct) * Math.sin(m * phi);
                    C[m * N + i] += y6c; S[m * N + i] += y6s;
                    C[m * N + j] += y6c; S[m * N + j] += y6s;
                    cnt[i]++; cnt[j]++;
                  }
                }
              }
              j = this.cellNext[j];
            }
          }
        }
      }
    }
    const q6 = new Float32Array(N);
    let mean = 0;
    const norm = Math.sqrt(4 * Math.PI / 13);
    for (let i = 0; i < N; i++) {
      const n_i = Math.max(1, cnt[i]);
      let acc = 0;
      for (let m = 0; m <= 6; m++) {
        const cr = C[m * N + i] / n_i, ci = S[m * N + i] / n_i;
        acc += (m === 0 ? 1 : 2) * (cr * cr + ci * ci);
      }
      const w = m === 0 ? 1 : 2;
      q6[i] = Math.sqrt(Math.max(0, norm * acc));
      mean += q6[i];
    }
    mean /= N;
    return { q6, mean };
  }

  /** 邻域平滑迁移率：每珠取 2σ 邻域内位移平方的均值（χ₄ 视角逐珠可视化） */  /** 邻域平滑迁移率：每珠取 2σ 邻域内位移平方的均值（χ₄ 视角逐珠可视化） */  /** 邻域平滑迁移率：每珠取 2σ 邻域内位移平方的均值（χ₄ 视角逐珠可视化） */
  smoothMobility() {
    const N = this.N;
    const p = this.pos, u = this.upos, sm = this.snapMob;
    const { Lx, Ly, Lz, ncx, ncy, ncz } = this;
    const m2 = new Float32Array(N);
    for (let i = 0; i < N; i++) {
      const dx = u[i * 3] - sm[i * 3], dy = u[i * 3 + 1] - sm[i * 3 + 1], dz = u[i * 3 + 2] - sm[i * 3 + 2];
      m2[i] = dx * dx + dy * dy + dz * dz;
    }
    const out = new Float32Array(N);
    this._buildCells();
    const count = new Int32Array(N);
    for (let i = 0; i < N; i++) {
      const i3 = i * 3;
      const xi = p[i3], yi = p[i3 + 1], zi = p[i3 + 2];
      let cx = (xi / Lx * ncx) | 0, cy = (yi / Ly * ncy) | 0, cz = (zi / Lz * ncz) | 0;
      if (cx >= ncx) cx = ncx - 1; if (cy >= ncy) cy = ncy - 1; if (cz >= ncz) cz = ncz - 1;
      let sum = 0, cnt = 0;
      for (let oz = -1; oz <= 1; oz++) {
        let z2 = cz + oz; if (z2 < 0) z2 += ncz; else if (z2 >= ncz) z2 -= ncz;
        for (let oy = -1; oy <= 1; oy++) {
          let y2 = cy + oy; if (y2 < 0) y2 += ncy; else if (y2 >= ncy) y2 -= ncy;
          for (let ox = -1; ox <= 1; ox++) {
            let x2 = cx + ox; if (x2 < 0) x2 += ncx; else if (x2 >= ncx) x2 -= ncx;
            let j = this.cellHead[x2 + ncx * (y2 + ncy * z2)];
            while (j !== -1) {
              if (j !== i) {
                const j3 = j * 3;
                let dx = xi - p[j3];     dx -= Lx * Math.round(dx / Lx);
                let dy = yi - p[j3 + 1]; dy -= Ly * Math.round(dy / Ly);
                let dz = zi - p[j3 + 2]; dz -= Lz * Math.round(dz / Lz);
                if (dx * dx + dy * dy + dz * dz < 4.0) { // 2σ 邻域
                  sum += m2[j];
                  cnt++;
                }
              }
              j = this.cellNext[j];
            }
          }
        }
      }
      out[i] = cnt ? sum / cnt : m2[i];
    }
    return out;
  }

  /** 瞬时压强（理想项 + 维里项），LJ 约化单位 */
  get pressure() {
    const V = this.Lx * this.Ly * this.Lz;
    return (this.N * this.keTemp + (this.wxx + this.wyy + this.wzz) / 3) / (3 * V);
  }

  /** 偏应力（单轴用）与剪切应力（LE 用） */
  get stressDev() {
    const V = this.Lx * this.Ly * this.Lz;
    return (this.N * this.keTemp + this.wxx - 0.5 * (this.wyy + this.wzz)) / V;
  }
  get stressXY() {
    return this.wxy / (this.Lx * this.Ly * this.Lz);
  }

  /**
   * 粗粒化四点易感性 χ₄(s=5τ)：以 5τ 迁移率窗口为滞后，位移阈值 0.3σ，
   * 3σ 粗粒化格子上的占有分数涨落 × N（实用估计量）
   */
  chi4Coarse() {
    const cs = 3, a2 = 0.09;
    const sm = this.snapMob;
    const m3x = Math.max(1, Math.floor(this.Lx / cs));
    const m3y = Math.max(1, Math.floor(this.Ly / cs));
    const m3z = Math.max(1, Math.floor(this.Lz / cs));
    const M = m3x * m3y * m3z;
    const qsum = new Float64Array(M), cnt = new Float64Array(M);
    const p = this.pos, u = this.upos;
    for (let i = 0; i < this.N; i++) {
      const i3 = i * 3;
      const cx = Math.min(m3x - 1, (p[i3] / this.Lx * m3x) | 0);
      const cy = Math.min(m3y - 1, (p[i3 + 1] / this.Ly * m3y) | 0);
      const cz = Math.min(m3z - 1, (p[i3 + 2] / this.Lz * m3z) | 0);
      const ci = cx + m3x * (cy + m3y * cz);
      const dx = u[i3] - sm[i3], dy = u[i3 + 1] - sm[i3 + 1], dz = u[i3 + 2] - sm[i3 + 2];
      const occ = dx * dx + dy * dy + dz * dz < a2 ? 1 : 0;
      qsum[ci] += occ;
      cnt[ci] += 1;
    }
    let wq = 0, wqq = 0, wt = 0;
    for (let k = 0; k < M; k++) {
      if (cnt[k] === 0) continue;
      const q = qsum[k] / cnt[k];
      wq += cnt[k] * q;
      wqq += cnt[k] * q * q;
      wt += cnt[k];
    }
    if (wt === 0) return 0;
    const mean = wq / wt;
    const varW = wqq / wt - mean * mean;
    return this.N * Math.max(0, varW) / (wt / Math.max(1, m3x * m3y * m3z * 0.999));
  }

  /** 非折叠坐标快照（配合 msdOver 做自定义窗口的 MSD 测量） */
  snapshotU() { return this.upos.slice(); }

  /** 相对参考点的四阶位移矩 ⟨|Δr|⁴⟩（非高斯参数 α₂ 用） */
  msd4Ref() {
    let s = 0;
    const u = this.upos, r = this.refPos;
    for (let i = 0; i < this.N; i++) {
      const i3 = i * 3;
      const dx = u[i3] - r[i3], dy = u[i3 + 1] - r[i3 + 1], dz = u[i3 + 2] - r[i3 + 2];
      const d2 = dx * dx + dy * dy + dz * dz;
      s += d2 * d2;
    }
    return s / this.N;
  }

  /** 相对给定快照的均方位移 */
  msdOver(snap) {
    let s = 0;
    const u = this.upos;
    for (let a = 0; a < u.length; a++) {
      const d = u[a] - snap[a];
      s += d * d;
    }
    return s / this.N;
  }

  /** 最长键长（数值稳定性体检指标） */
  maxBondLength() {
    const p = this.pos;
    const { Lx, Ly, Lz } = this;
    let mx = 0;
    for (let b = 0; b < this.bondPairs.length; b += 2) {
      const i3 = this.bondPairs[b] * 3, j3 = this.bondPairs[b + 1] * 3;
      let dx = p[i3] - p[j3];
      let dy = p[i3 + 1] - p[j3 + 1];
      let dz = p[i3 + 2] - p[j3 + 2];
      dx -= Lx * Math.round(dx / Lx);
      dy -= Ly * Math.round(dy / Ly);
      dz -= Lz * Math.round(dz / Lz);
      const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
      if (r > mx) mx = r;
    }
    return mx;
  }
}
