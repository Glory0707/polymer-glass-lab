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
    this.scene.background = new THREE.Color(0x0b0e14);
    this.scene.fog = new THREE.Fog(0x0b0e14, 0, 0); // far 在 _fitCamera 里设

    this.camera = new THREE.PerspectiveCamera(50, 1, 0.1, 500);

    const hemi = new THREE.HemisphereLight(0x99b7ff, 0x1a2030, 1.1);
    this.scene.add(hemi);
    const key = new THREE.DirectionalLight(0xffffff, 2.2);
    const fill = new THREE.DirectionalLight(0x88aaff, 0.7);
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
    this.controls.target.copy(this._center);

    this._resize();
    this._ro = new ResizeObserver(() => this._resize());
    this._ro.observe(container);
  }

  /** 颜色缓冲区（sRGB 数值线性化前的 pow2.2 由写入方处理） */
  get colorTarget() { return this.beadMesh.instanceColor.array; }

  _buildSimObjects(sim) {
    const geo = new THREE.SphereGeometry(BEAD_R, 14, 10);
    const mat = new THREE.MeshLambertMaterial();
    const mesh = new THREE.InstancedMesh(geo, mat, sim.N);
    mesh.frustumCulled = false;
    const white = new THREE.Color(0xffffff);
    for (let i = 0; i < sim.N; i++) mesh.setColorAt(i, white);
    this.beadMesh = mesh;
    this.scene.add(mesh);

    // 键线（含跨周期边界的镜像补画）
    const nb = sim.bondPairs.length / 2;
    this.bondPos = new Float32Array(nb * 6);
    this.bondCol = new Float32Array(nb * 6);
    const bgeo = new THREE.BufferGeometry();
    bgeo.setAttribute('position', new THREE.BufferAttribute(this.bondPos, 3));
    bgeo.setAttribute('color', new THREE.BufferAttribute(this.bondCol, 3));
    this.bondLines = new THREE.LineSegments(
      bgeo,
      new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.45 })
    );
    this.bondLines.frustumCulled = false;
    this.scene.add(this.bondLines);

    // 盒子线框
    const boxGeo = new THREE.BoxGeometry(sim.Lx, sim.Ly, sim.Lz);
    const edges = new THREE.EdgesGeometry(boxGeo);
    this.boxHelper = new THREE.LineSegments(
      edges,
      new THREE.LineBasicMaterial({ color: 0x39435a, transparent: true, opacity: 0.8 })
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
    for (let i = 0; i < sim.N; i++) {
      const i3 = i * 3;
      m.makeTranslation(p[i3], p[i3 + 1], p[i3 + 2]);
      mesh.setMatrixAt(i, m);
    }
    mesh.instanceMatrix.needsUpdate = true;
    if (mesh.instanceColor) mesh.instanceColor.needsUpdate = true;

    if (opts.showBonds) {
      const arr = this.bondPos, col = this.bondCol;
      const bp = sim.bondPairs;
      const ct = mesh.instanceColor ? mesh.instanceColor.array : null;
      const { Lx, Ly, Lz } = sim;
      for (let b = 0, w = 0; b < bp.length; b += 2, w += 6) {
        const i3 = bp[b] * 3, j3 = bp[b + 1] * 3;
        arr[w] = p[i3]; arr[w + 1] = p[i3 + 1]; arr[w + 2] = p[i3 + 2];
        // 跨周期边界的键用最小镜像补画，避免横穿盒子的长线
        let sx = 0, sy = 0, sz = 0;
        let dx = p[j3] - p[i3];
        if (dx > Lx / 2) sx = Lx; else if (dx < -Lx / 2) sx = -Lx;
        let dy = p[j3 + 1] - p[i3 + 1];
        if (dy > Ly / 2) sy = Ly; else if (dy < -Ly / 2) sy = -Ly;
        let dz = p[j3 + 2] - p[i3 + 2];
        if (dz > Lz / 2) sz = Lz; else if (dz < -Lz / 2) sz = -Lz;
        arr[w + 3] = p[j3] + sx;
        arr[w + 4] = p[j3 + 1] + sy;
        arr[w + 5] = p[j3 + 2] + sz;
        if (ct) {
          const r = (ct[i3] + ct[j3]) * 0.5;
          const g = (ct[i3 + 1] + ct[j3 + 1]) * 0.5;
          const b2 = (ct[i3 + 2] + ct[j3 + 2]) * 0.5;
          col[w] = r; col[w + 1] = g; col[w + 2] = b2;
          col[w + 3] = r; col[w + 4] = g; col[w + 5] = b2;
        }
      }
      this.bondLines.geometry.attributes.position.needsUpdate = true;
      this.bondLines.geometry.attributes.color.needsUpdate = true;
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
    this.renderer.setSize(w, h, false);
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
  }

  dispose() {
    this._ro.disconnect();
    this.controls.dispose();
    this.scene.traverse((o) => {
      if (o.geometry) o.geometry.dispose();
      if (o.material) o.material.dispose();
    });
    this.renderer.dispose();
    this.renderer.domElement.remove();
  }
}
