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
    this.cellHead = new Int32Array(this.ncx * this.ncy * this.ncz);
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
                dx -= Lx * Math.round(dx / Lx);
                dy -= Ly * Math.round(dy / Ly);
                dz -= Lz * Math.round(dz / Lz);
                const r2 = dx * dx + dy * dy + dz * dz;
                if (r2 < RC2) {
                  const inv2 = 1 / r2;
                  const inv6 = inv2 * inv2 * inv2;
                  let fc = 24 * (2 * inv6 * inv6 - inv6) * inv2; // WCA
                  // 数值护栏：仅在初始化瞬间的深重叠时触发（正常动力学 fc ≲ 100）
                  if (fc > 2e3) fc = 2e3;
                  const fxc = fc * dx, fyc = fc * dy, fzc = fc * dz;
                  f[i3] += fxc; f[i3 + 1] += fyc; f[i3 + 2] += fzc;
                  f[j3] -= fxc; f[j3 + 1] -= fyc; f[j3 + 2] -= fzc;
                  let uw = 4 * (inv6 * inv6 - inv6) + 1; // 已平移至截断处为零
                  if (uw > 2e3) uw = 2e3;
                  pe += uw;
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

  /** 非折叠坐标快照（配合 msdOver 做自定义窗口的 MSD 测量） */
  snapshotU() { return this.upos.slice(); }

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
