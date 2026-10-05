// The computational task network.
//
// What is drawn here is the pipeline's real structure with its real state:
// one node per stage, coloured by that stage's aggregate execution state,
// connected by the dependencies the pipeline itself declares. Edges carry
// small travelling marks only when a task on them has actually just
// finished — if nothing is running, the network is still, because inventing
// motion would be the single most dishonest thing this interface could do.
//
// Rendering is Canvas 2D with additive glow. That is a deliberate choice
// over WebGL: the scene is 16 nodes and ~32 edges, so the bottleneck would
// be the layout maths, not fill rate, and Canvas 2D keeps the whole thing
// dependency-free and debuggable.

import { STAGE_LIST } from './data.js';
import { nodeColor, stateColor, stateStyle, isLive } from './palette.js';
import { clamp } from './util.js';

const TAU = Math.PI * 2;

export class NetworkView {
  constructor(host, canvas, options = {}) {
    this.host = host;
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.onSelect = options.onSelect || (() => {});
    this.onHover = options.onHover || (() => {});

    this.nodes = [];
    this.edges = [];
    this.subjectStates = {};   // "stage|subject" -> state, for genome view
    this.subjects = [];

    this.scale = 1;
    this.offsetX = 0;
    this.offsetY = 0;
    this.fitted = false;

    this.selected = null;     // stage key
    this.hovered = null;
    this.edgeFilter = 'both';
    this.view = 'stage';

    this.pulses = [];         // {edge, t, life}
    this.particles = [];      // ambient field, deterministic
    this.time = 0;
    this.lastFrame = 0;
    this.frozen = false;
    this.idle = true;

    this._bindInput();
    this._seedParticles();
    this.resize();
  }

  // ── data ────────────────────────────────────────────────────────

  setGraph(payload) {
    const nodes = payload.nodes || [];
    const edges = payload.edges || [];
    this.nodes = nodes.map((n) => ({ ...n, counts: n.counts || {}, r: 15 }));
    this.edges = edges;
    const index = new Map(this.nodes.map((n) => [n.stage, n]));
    this.byStage = index;
    // Cache the endpoints so the draw loop never does a lookup.
    this.edges = edges.map((e) => ({ ...e, a: index.get(e.source), b: index.get(e.target) }))
      .filter((e) => e.a && e.b);
    this.stages = payload.stages || STAGE_LIST;
    if (!this.fitted && this.nodes.length) this.fit();
  }

  setState(snapshot) {
    const byStage = new Map((snapshot.nodes || []).map((n) => [n.stage, n]));
    for (const node of this.nodes) {
      const live = byStage.get(node.stage);
      if (live) {
        node.state = live.state;
        node.counts = live.counts || {};
        node.total = live.total;
      }
    }
    this.subjects = snapshot.subjects || [];
    this.subjectStates = {};
    for (const task of snapshot.tasks || []) {
      this.subjectStates[`${task.stage}|${task.subject}`] = task.state;
    }
    this.runKey = snapshot.run_key;
    this.idle = !!snapshot.idle;
    this.finishedTotal = (snapshot.counts || {}).SUCCEEDED || 0;
  }

  /**
   * Launch a travelling mark along the edge that a finished task used.
   *
   * Called from a real TASK_COMPLETED / TASK_VALIDATED event, so the motion
   * is a consequence of observed work rather than a loop.
   */
  markActivity(stage, targetStage) {
    if (targetStage) {
      const edge = this.edges.find(
        (e) => e.source === stage && e.target === targetStage,
      ) || this.edges.find((e) => e.target === targetStage);
      if (edge) {
        this.pulses.push({ edge, t: 0, life: 1.5 });
        return;
      }
    }
    // No declared downstream edge: pulse the node's own ring instead of
    // inventing a connection.
    const node = this.byStage && this.byStage.get(stage);
    if (node) node.flash = 1.6;
  }

  setEdgeFilter(mode) { this.edgeFilter = mode; }
  setView(mode) { this.view = mode; }
  setFrozen(on) { this.frozen = on; }
  select(stage) { this.selected = stage; }

  // ── camera ──────────────────────────────────────────────────────

  resize() {
    const ratio = window.devicePixelRatio || 1;
    const rect = this.host.getBoundingClientRect();
    this.width = Math.max(1, rect.width);
    this.height = Math.max(1, rect.height);
    this.canvas.width = Math.round(this.width * ratio);
    this.canvas.height = Math.round(this.height * ratio);
    this.canvas.style.width = `${this.width}px`;
    this.canvas.style.height = `${this.height}px`;
    this.ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
  }

  fit() {
    if (!this.nodes.length) return;
    let minX = Infinity; let maxX = -Infinity;
    let minY = Infinity; let maxY = -Infinity;
    for (const n of this.nodes) {
      minX = Math.min(minX, n.x); maxX = Math.max(maxX, n.x);
      minY = Math.min(minY, n.y); maxY = Math.max(maxY, n.y);
    }
    const pad = 130;
    const w = (maxX - minX) || 1;
    const h = (maxY - minY) || 1;
    this.scale = clamp(
      Math.min((this.width - pad * 2) / w, (this.height - pad * 2) / h), 0.12, 2.2,
    );
    this.offsetX = this.width / 2 - ((minX + maxX) / 2) * this.scale;
    this.offsetY = this.height / 2 - ((minY + maxY) / 2) * this.scale;
    this.fitted = true;
  }

  zoomBy(factor, cx, cy) {
    const before = this.scale;
    this.scale = clamp(this.scale * factor, 0.1, 4);
    if (cx === undefined) { cx = this.width / 2; cy = this.height / 2; }
    // Keep the point under the cursor fixed.
    this.offsetX = cx - (cx - this.offsetX) * (this.scale / before);
    this.offsetY = cy - (cy - this.offsetY) * (this.scale / before);
  }

  focus(stage) {
    const node = this.byStage && this.byStage.get(stage);
    if (!node) return;
    this.scale = Math.max(this.scale, 0.85);
    this.offsetX = this.width / 2 - node.x * this.scale;
    this.offsetY = this.height / 2 - node.y * this.scale;
  }

  toScreen(n) {
    return { x: n.x * this.scale + this.offsetX, y: n.y * this.scale + this.offsetY };
  }

  pick(px, py) {
    let best = null;
    let bestD = 30 * 30;
    for (const n of this.nodes) {
      const p = this.toScreen(n);
      const d = (p.x - px) ** 2 + (p.y - py) ** 2;
      if (d < bestD) { bestD = d; best = n; }
    }
    return best;
  }

  // ── input ───────────────────────────────────────────────────────

  _bindInput() {
    const c = this.canvas;
    let dragging = false;
    let lastX = 0;
    let lastY = 0;
    let moved = 0;

    c.addEventListener('pointerdown', (e) => {
      dragging = true; moved = 0;
      lastX = e.clientX; lastY = e.clientY;
      c.setPointerCapture(e.pointerId);
    });

    c.addEventListener('pointermove', (e) => {
      const rect = c.getBoundingClientRect();
      const px = e.clientX - rect.left;
      const py = e.clientY - rect.top;
      if (dragging) {
        const dx = e.clientX - lastX;
        const dy = e.clientY - lastY;
        moved += Math.abs(dx) + Math.abs(dy);
        this.offsetX += dx;
        this.offsetY += dy;
        lastX = e.clientX; lastY = e.clientY;
        return;
      }
      const node = this.pick(px, py);
      const stage = node ? node.stage : null;
      if (stage !== this.hovered) {
        this.hovered = stage;
        this.onHover(stage, node ? this.toScreen(node) : null);
        c.style.cursor = node ? 'pointer' : 'default';
      }
    });

    const endDrag = (e) => {
      if (!dragging) return;
      dragging = false;
      try { c.releasePointerCapture(e.pointerId); } catch { /* already released */ }
      if (moved < 5) {
        const rect = c.getBoundingClientRect();
        const node = this.pick(e.clientX - rect.left, e.clientY - rect.top);
        this.onSelect(node ? node.stage : null);
      }
    };
    c.addEventListener('pointerup', endDrag);
    c.addEventListener('pointercancel', endDrag);
    c.addEventListener('pointerleave', () => {
      if (this.hovered) { this.hovered = null; this.onHover(null, null); }
    });

    c.addEventListener('wheel', (e) => {
      e.preventDefault();
      const rect = c.getBoundingClientRect();
      this.zoomBy(e.deltaY < 0 ? 1.12 : 1 / 1.12, e.clientX - rect.left, e.clientY - rect.top);
    }, { passive: false });
  }

  // ── ambient field ───────────────────────────────────────────────

  /**
   * A faint static starfield, generated once from a fixed seed.
   *
   * Deterministic on purpose: a field that reshuffles on every reload looks
   * like activity, and this network must not look busy when it is idle.
   */
  _seedParticles() {
    let seed = 20240601;
    const rand = () => {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      return seed / 0x7fffffff;
    };
    this.particles = Array.from({ length: 130 }, () => ({
      u: rand(), v: rand(), s: 0.4 + rand() * 1.1, a: 0.05 + rand() * 0.16,
    }));
  }

  // ── draw ────────────────────────────────────────────────────────

  start() {
    const loop = (ts) => {
      const dt = this.lastFrame ? Math.min(0.05, (ts - this.lastFrame) / 1000) : 0.016;
      this.lastFrame = ts;
      if (!this.frozen) {
        this.time += dt;
        this._advance(dt);
      }
      this.draw();
      this.raf = requestAnimationFrame(loop);
    };
    this.raf = requestAnimationFrame(loop);
  }

  stop() { if (this.raf) cancelAnimationFrame(this.raf); }

  _advance(dt) {
    for (let i = this.pulses.length - 1; i >= 0; i -= 1) {
      this.pulses[i].t += dt;
      if (this.pulses[i].t >= this.pulses[i].life) this.pulses.splice(i, 1);
    }
    for (const n of this.nodes) {
      if (n.flash > 0) n.flash = Math.max(0, n.flash - dt);
    }
  }

  _visibleEdges() {
    if (this.edgeFilter === 'both') return this.edges;
    return this.edges.filter((e) => e.kind === this.edgeFilter);
  }

  draw() {
    const { ctx } = this;
    ctx.clearRect(0, 0, this.width, this.height);

    this._drawField();
    this._drawEdges();
    this._drawPulses();
    this._drawNodes();
  }

  _drawField() {
    const { ctx } = this;
    // A faint grid gives depth without becoming a chart.
    const step = 64 * this.scale;
    if (step > 14) {
      ctx.strokeStyle = 'rgba(59,130,246,0.045)';
      ctx.lineWidth = 1;
      ctx.beginPath();
      const ox = this.offsetX % step;
      const oy = this.offsetY % step;
      for (let x = ox; x < this.width; x += step) {
        ctx.moveTo(x, 0); ctx.lineTo(x, this.height);
      }
      for (let y = oy; y < this.height; y += step) {
        ctx.moveTo(0, y); ctx.lineTo(this.width, y);
      }
      ctx.stroke();
    }
    for (const p of this.particles) {
      const x = p.u * this.width;
      const y = p.v * this.height;
      ctx.fillStyle = `rgba(125,211,252,${p.a})`;
      ctx.beginPath();
      ctx.arc(x, y, p.s, 0, TAU);
      ctx.fill();
    }
  }

  /** A quadratic curve, so the network reads as flow rather than a tree. */
  _curve(a, b, bend) {
    const mx = (a.x + b.x) / 2;
    const my = (a.y + b.y) / 2;
    const dx = b.x - a.x;
    const dy = b.y - a.y;
    const len = Math.hypot(dx, dy) || 1;
    return {
      cx: mx - (dy / len) * bend,
      cy: my + (dx / len) * bend,
    };
  }

  _pointOnCurve(a, c, b, t) {
    const u = 1 - t;
    return {
      x: u * u * a.x + 2 * u * t * c.x + t * t * b.x,
      y: u * u * a.y + 2 * u * t * c.y + t * t * b.y,
    };
  }

  _drawEdges() {
    const { ctx } = this;
    for (const e of this._visibleEdges()) {
      const a = this.toScreen(e.a);
      const b = this.toScreen(e.b);
      const seq = e.kind === 'sequence';
      const dist = Math.hypot(b.x - a.x, b.y - a.y);
      const { cx, cy } = this._curve(a, b, dist * 0.12);

      const fromLive = isLive(e.a.state);
      const toLive = isLive(e.b.state);
      let alpha = seq ? 0.055 : 0.13;
      if (!seq && (fromLive || toLive)) alpha = 0.3;
      if (e.source === this.selected || e.target === this.selected) alpha += 0.3;

      ctx.save();
      ctx.setLineDash(seq ? [2, 6] : []);
      ctx.lineWidth = seq ? 1 : 1.3;
      ctx.strokeStyle = seq
        ? `rgba(100,116,139,${alpha})`
        : `rgba(34,211,238,${alpha})`;
      ctx.shadowBlur = 0;
      ctx.beginPath();
      ctx.moveTo(a.x, a.y);
      ctx.quadraticCurveTo(cx, cy, b.x, b.y);
      ctx.stroke();
      ctx.restore();
    }
  }

  _drawPulses() {
    const { ctx } = this;
    for (const p of this.pulses) {
      const a = this.toScreen(p.edge.a);
      const b = this.toScreen(p.edge.b);
      const dist = Math.hypot(b.x - a.x, b.y - a.y);
      const { cx, cy } = this._curve(a, b, dist * 0.12);
      const t = p.t / p.life;
      const pt = this._pointOnCurve(a, { x: cx, y: cy }, b, t);
      const fade = 1 - t;
      ctx.save();
      ctx.globalCompositeOperation = 'lighter';
      ctx.fillStyle = `rgba(74,222,128,${0.85 * fade})`;
      ctx.shadowColor = '#4ade80';
      ctx.shadowBlur = 14;
      ctx.beginPath();
      ctx.arc(pt.x, pt.y, 2.6, 0, TAU);
      ctx.fill();
      ctx.restore();
    }
  }

  _drawNodes() {
    const { ctx } = this;
    const showGenome = this.view === 'genome' && this.subjects.length > 0;

    for (const node of this.nodes) {
      const p = this.toScreen(node);
      if (p.x < -160 || p.x > this.width + 160 || p.y < -160 || p.y > this.height + 160) {
        continue;
      }
      const color = nodeColor(node.stage, node.state);
      const style = stateStyle(node.state);
      const selected = node.stage === this.selected;
      const hovered = node.stage === this.hovered;
      const live = isLive(node.state);
      const base = 13 + Math.min(9, (node.total || 0) * 0.12);
      const radius = base * (selected ? 1.22 : hovered ? 1.1 : 1);

      ctx.save();
      ctx.globalCompositeOperation = 'lighter';

      // Halo. Its strength is the state's glow, so a running stage is
      // visibly hotter than a completed one at a glance.
      if (style.glow > 0) {
        const pulse = live ? 0.82 + 0.18 * Math.sin(this.time * 2.6) : 1;
        const grad = ctx.createRadialGradient(p.x, p.y, 0, p.x, p.y, radius * 5.4);
        grad.addColorStop(0, this._alpha(color, 0.5 * style.glow * pulse));
        grad.addColorStop(0.42, this._alpha(color, 0.14 * style.glow * pulse));
        grad.addColorStop(1, this._alpha(color, 0));
        ctx.fillStyle = grad;
        ctx.beginPath();
        ctx.arc(p.x, p.y, radius * 5.4, 0, TAU);
        ctx.fill();
      }

      // Filaments: fixed spokes, denser when work is in flight. These are
      // static structure, not motion, so an idle node stays quiet.
      const spokes = live ? 13 : 8;
      ctx.strokeStyle = this._alpha(color, live ? 0.42 : 0.2);
      ctx.lineWidth = 1;
      for (let i = 0; i < spokes; i += 1) {
        const ang = (i / spokes) * TAU + node.index * 0.37;
        const inner = radius * 1.25;
        const outer = radius * (live ? 4.3 : 3.1);
        ctx.beginPath();
        ctx.moveTo(p.x + Math.cos(ang) * inner, p.y + Math.sin(ang) * inner);
        ctx.lineTo(p.x + Math.cos(ang) * outer, p.y + Math.sin(ang) * outer);
        ctx.stroke();
      }
      ctx.restore();

      // Core.
      ctx.save();
      ctx.globalCompositeOperation = 'lighter';
      ctx.fillStyle = this._alpha(color, selected || live ? 0.95 : 0.7);
      ctx.shadowColor = color;
      ctx.shadowBlur = selected ? 26 : live ? 20 : 11;
      ctx.beginPath();
      ctx.arc(p.x, p.y, radius, 0, TAU);
      ctx.fill();
      ctx.restore();

      // Ring: selection and hover are geometry, not colour, so they stay
      // legible for a stage whose own hue is already busy.
      if (selected || hovered) {
        ctx.strokeStyle = selected ? '#e2e8f0' : 'rgba(226,232,240,0.45)';
        ctx.lineWidth = selected ? 1.8 : 1;
        ctx.setLineDash(selected ? [] : [3, 4]);
        ctx.beginPath();
        ctx.arc(p.x, p.y, radius + 8, 0, TAU);
        ctx.stroke();
        ctx.setLineDash([]);
      }

      if (node.flash > 0) {
        const t = node.flash / 1.6;
        ctx.strokeStyle = `rgba(74,222,128,${t * 0.8})`;
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.arc(p.x, p.y, radius + 8 + (1 - t) * 34, 0, TAU);
        ctx.stroke();
      }

      this._drawLabel(node, p, color, radius, showGenome);
    }
  }

  _drawLabel(node, p, color, radius, showGenome) {
    if (this.scale < 0.3) return;
    const { ctx } = this;
    const count = showGenome ? (node.counts && Object.values(node.counts).reduce((a, b) => a + b, 0)) : node.total;
    const main = node.label;
    const sub = count ? `(${count})` : '';
    const padX = 7;
    const lineH = 13;
    const wName = ctx.measureText(main).width;
    const wSub = sub ? ctx.measureText(sub).width : 0;
    const boxW = Math.max(wName + wSub + padX * 2, 58);
    const boxH = lineH + 8;
    const bx = p.x - boxW / 2;
    const by = p.y + radius + 13;

    ctx.save();
    ctx.fillStyle = 'rgba(8,13,22,0.86)';
    ctx.strokeStyle = this._alpha(color, 0.5);
    ctx.lineWidth = 1;
    const r = 4;
    ctx.beginPath();
    ctx.moveTo(bx + r, by);
    ctx.arcTo(bx + boxW, by, bx + boxW, by + boxH, r);
    ctx.arcTo(bx + boxW, by + boxH, bx, by + boxH, r);
    ctx.arcTo(bx, by + boxH, bx, by, r);
    ctx.arcTo(bx, by, bx + boxW, by, r);
    ctx.closePath();
    ctx.fill();
    ctx.stroke();

    ctx.font = '600 10.5px -apple-system, "Segoe UI", sans-serif';
    ctx.textBaseline = 'middle';
    ctx.fillStyle = '#dbe4f0';
    ctx.fillText(main, bx + padX, by + boxH / 2 + 0.5);
    if (sub) {
      ctx.fillStyle = this._alpha(color, 0.95);
      ctx.fillText(sub, bx + padX + wName + 3, by + boxH / 2 + 0.5);
    }
    ctx.restore();
  }

  _alpha(hex, a) {
    const h = hex.replace('#', '');
    const n = parseInt(h.length === 3 ? h.split('').map((c) => c + c).join('') : h, 16);
    return `rgba(${(n >> 16) & 255},${(n >> 8) & 255},${n & 255},${a})`;
  }
}
