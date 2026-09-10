/**
 * renderer.js — Three.js 实时渲染
 * 珠子用 InstancedMesh（一次 draw call），键用顶点着色 LineSegments，
 * 颜色缓冲区直接暴露给外部按迁移率/链写入。
 */
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

const BEAD_R = 0.5;

export class GlassRenderer {
  constructor(container, sim) {
    this.container = container;
    this.sim = sim;

    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color(0x0c0d10);
    this.scene.fog = new THREE.Fog(0x0c0d10, 0, 0); // far 在 _fitCamera 里设

    this.camera = new THREE.PerspectiveCamera(50, 1, 0.1, 500);

    const hemi = new THREE.HemisphereLight(0x99b7ff, 0x1a2030, 1.25);
    this.scene.add(hemi);
    const key = new THREE.DirectionalLight(0xffffff, 2.6);
    const fill = new THREE.DirectionalLight(0x88aaff, 0.8);
    this.scene.add(key, fill);
    this._key = key; this._fill = fill;

    this._m = new THREE.Matrix4();

    this._buildSimObjects(sim);
    this._fitCamera(sim);

    this.renderer = new THREE.WebGLRenderer({ antialias: true });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    this.renderer.toneMapping = THREE.ACESFilmicToneMapping;
    this.renderer.toneMappingExposure = 1.1;
    container.appendChild(this.renderer.domElement);
    this.renderer.domElement.style.display = 'block';

    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.12;
    this.controls.enablePan = false;
    this.controls.target.copy(this._center);
    // 空闲自转：页面永远是活的；用户一上手就停，放开 6 秒后恢复
    this.controls.autoRotate = true;
    this.controls.autoRotateSpeed = 0.45;
    this._idleTimer = null;
    this.controls.addEventListener('start', () => {
      this.controls.autoRotate = false;
      clearTimeout(this._idleTimer);
    });
    this.controls.addEventListener('end', () => {
      clearTimeout(this._idleTimer);
      this._idleTimer = setTimeout(() => { this.controls.autoRotate = true; }, 6000);
    });

    this._resize();
    this._ro = new ResizeObserver(() => this._resize());
    this._ro.observe(container);
  }

  /** 颜色缓冲区（sRGB 数值线性化前的 pow2.2 由写入方处理） */
  get colorTarget() { return this.beadMesh.instanceColor.array; }

  _buildSimObjects(sim) {
    const geo = new THREE.SphereGeometry(BEAD_R, 14, 10);
    const mat = new THREE.MeshStandardMaterial({ roughness: 0.38, metalness: 0.06 });
    const mesh = new THREE.InstancedMesh(geo, mat, sim.N);
    mesh.frustumCulled = false;
    const white = new THREE.Color(0xffffff);
    for (let i = 0; i < sim.N; i++) mesh.setColorAt(i, white);
    this.beadMesh = mesh;
    this.scene.add(mesh);

    // 键线（含跨周期边界的镜像补画）：中性灰弱化存在感，让珠色独占数据表达
    const nb = sim.bondPairs.length / 2;
    this.bondPos = new Float32Array(nb * 6);
    const bgeo = new THREE.BufferGeometry();
    bgeo.setAttribute('position', new THREE.BufferAttribute(this.bondPos, 3));
    this.bondLines = new THREE.LineSegments(
      bgeo,
      new THREE.LineBasicMaterial({ color: 0x52545a, transparent: true, opacity: 0.28 })
    );
    this.bondLines.frustumCulled = false;
    this.scene.add(this.bondLines);

    // 盒子参考框：只画 8 个角的取景括号——整根棱边在近距透视下会把远端
    // 投影得很大，看起来像从珠子块里辐射出去的长线
    const hw = sim.Lx / 2, hh = sim.Ly / 2, hd = sim.Lz / 2;
    const k = 0.9; // 括号臂长 (σ)
    const pts = [];
    for (const sx of [-1, 1]) {
      for (const sy of [-1, 1]) {
        for (const sz of [-1, 1]) {
          const bx = sx * hw, by = sy * hh, bz = sz * hd;
          pts.push(bx, by, bz, bx - sx * k, by, bz);
          pts.push(bx, by, bz, bx, by - sy * k, bz);
          pts.push(bx, by, bz, bx, by, bz - sz * k);
        }
      }
    }
    const bracketGeo = new THREE.BufferGeometry();
    bracketGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(pts), 3));
    this.boxHelper = new THREE.LineSegments(
      bracketGeo,
      new THREE.LineBasicMaterial({ color: 0x3a3d45, transparent: true, opacity: 0.9 })
    );
    this.boxHelper.frustumCulled = false;
    this.scene.add(this.boxHelper);
  }

  _fitCamera(sim) {
    const c = new THREE.Vector3(sim.Lx / 2, sim.Ly / 2, sim.Lz / 2);
    this._center = c;
    this.boxHelper.position.copy(c);
    const Lmax = Math.max(sim.Lx, sim.Ly, sim.Lz);
    this.camera.position.set(c.x + Lmax * 0.85, c.y + Lmax * 0.6, c.z + Lmax * 1.25);
    this.camera.lookAt(c);
    this.camera.near = Lmax * 0.01;
    this.camera.far = Lmax * 20;
    this.camera.updateProjectionMatrix();
    if (this.controls) this.controls.target.copy(c);
    this.scene.fog.near = Lmax * 1.2;
    this.scene.fog.far = Lmax * 4;
    this._key.position.set(c.x + Lmax, c.y + 2 * Lmax, c.z + 1.2 * Lmax);
    this._fill.position.set(c.x - Lmax, c.y - Lmax, c.z - 1.5 * Lmax);
  }

  /** 每帧同步：位置矩阵、实例颜色、键线几何 */
  update(opts = {}) {
    const sim = this.sim;
    const p = sim.pos;
    const m = this._m;
    const mesh = this.beadMesh;
    const sig = sim.sigma;
    for (let i = 0; i < sim.N; i++) {
      const i3 = i * 3;
      const s = sig ? sig[i] : 1;
      m.makeScale(s, s, s);
      m.elements[12] = p[i3];
      m.elements[13] = p[i3 + 1];
      m.elements[14] = p[i3 + 2];
      mesh.setMatrixAt(i, m);
    }
    mesh.instanceMatrix.needsUpdate = true;
    if (mesh.instanceColor) mesh.instanceColor.needsUpdate = true;

    if (opts.showBonds) {
      const arr = this.bondPos;
      const bp = sim.bondPairs;
      const { Lx, Ly, Lz } = sim;
      for (let b = 0, w = 0; b < bp.length; b += 2, w += 6) {
        const i3 = bp[b] * 3, j3 = bp[b + 1] * 3;
        arr[w] = p[i3]; arr[w + 1] = p[i3 + 1]; arr[w + 2] = p[i3 + 2];
        // 跨周期边界的键：端点 j 取其最近镜像（与力计算的最小镜像约定一致），
        // 否则跨界键会被画成横跨整个盒子的长线
        arr[w + 3] = p[j3]     - Lx * Math.round((p[j3] - p[i3]) / Lx);
        arr[w + 4] = p[j3 + 1] - Ly * Math.round((p[j3 + 1] - p[i3 + 1]) / Ly);
        arr[w + 5] = p[j3 + 2] - Lz * Math.round((p[j3 + 2] - p[i3 + 2]) / Lz);
      }
      this.bondLines.geometry.attributes.position.needsUpdate = true;
      this.bondLines.visible = true;
    } else {
      this.bondLines.visible = false;
    }
    this.boxHelper.visible = opts.showBox !== false;

    this.controls.update();
    this.renderer.render(this.scene, this.camera);
  }

  resetView() {
    this._fitCamera(this.sim);
    this.controls.target.copy(this._center);
  }

  _resize() {
    const w = this.container.clientWidth || 1;
    const h = this.container.clientHeight || 1;
    this.renderer.setSize(w, h); // 同步更新 canvas CSS 尺寸，高分缩放屏不再溢出
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
    this._refitDistance();
  }

  /** 按盒子包围球与当前视场角重算相机距离：任意窗口尺寸（含全屏）都完整取景，视角方向不变 */
  _refitDistance() {
    const sim = this.sim;
    const radius = 0.5 * Math.hypot(sim.Lx, sim.Ly, sim.Lz);
    const fovY = (this.camera.fov * Math.PI) / 180;
    const fovX = 2 * Math.atan(Math.tan(fovY / 2) * this.camera.aspect);
    const dist = (radius / Math.sin(Math.min(fovY, fovX) / 2)) * 1.04;
    const dir = this.camera.position.clone().sub(this.controls.target);
    if (dir.lengthSq() < 1e-6) dir.set(0.85, 0.6, 1.25);
    dir.normalize();
    this.controls.target.set(sim.Lx / 2, sim.Ly / 2, sim.Lz / 2);
    this.camera.position.copy(this.controls.target).add(dir.multiplyScalar(dist));
  }

  dispose() {
    this._ro.disconnect();
    clearTimeout(this._idleTimer);
    this.controls.dispose();
    this.scene.traverse((o) => {
      if (o.geometry) o.geometry.dispose();
      if (o.material) o.material.dispose();
    });
    this.renderer.dispose();
    this.renderer.domElement.remove();
  }
}
