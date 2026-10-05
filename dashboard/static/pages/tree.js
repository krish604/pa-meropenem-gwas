/* pages/tree.js — V1, the stage-9 phylogeny viewer.
 *
 * The trap this file is written around (UI-D6):
 *
 *   Node identity is the node's POSITION in the tree — the server assigns
 *   `node_id` as a postorder index — and it is NEVER the internal label.
 *   `test_data/similarity/real_tree.nwk` repeats `100/100` across three
 *   internal nodes and leaves one unlabelled, so a viewer that keys a Map by
 *   `node.label` collapses or mislabels. Every Map and Set below is keyed on
 *   `node_id`. The only place a label is read is:
 *     - to DRAW it, and
 *     - to JOIN a tip to its isolate row by `sample_id` (a join, not identity).
 *   `assertNodeIdentity` checks the three structural invariants the whole
 *   module depends on and refuses to draw if any fails.
 *
 * Distances come from branch lengths (`length`), never from labels. When the
 * tree carries no branch lengths the horizontal axis falls back to node depth
 * and the page says so rather than drawing a phylogram from nothing.
 *
 * Rendering is canvas. Zoom/pan is hand-written: `d3-zoom` cannot be used here
 * because the vendored set omits `d3-dispatch`, `d3-drag` and `d3-transition`,
 * so `d3.zoom` throws `e.dispatch is not a function` at construction (measured).
 * d3-scale / d3-array / d3-format / d3-interpolate are used where they help.
 */

import {
  ABSENT,
  absentSpan,
  clear,
  h,
  kv,
  message,
  note,
  panel,
  probeList,
  serverMessage,
} from '../app.js';

const MARGIN = 18;
const ROW_H = 16;
const X_SPAN = 820;
const LABEL_MIN_PX = 8;

const SIR_TOKENS = {
  R: 'tone-bad-fg',
  I: 'tone-warn-fg',
  S: 'tone-ok-fg',
  SDD: 'tone-busy-fg',
  ND: 'tone-unknown-fg',
};

const OPRD_TOKENS = {
  intact: 'tone-ok-fg',
  disrupted: 'tone-bad-fg',
  absent: 'tone-warn-fg',
  not_assessed: 'tone-unknown-fg',
};

function debounce(fn, ms) {
  let timer = null;
  return (...args) => {
    if (timer !== null) clearTimeout(timer);
    timer = setTimeout(() => {
      timer = null;
      fn(...args);
    }, ms);
  };
}

/**
 * The D6 guard. Checks, on the server's own structural JSON:
 *   1. every `node_id` is unique;
 *   2. every child and parent id resolves to a node;
 *   3. a repeated label maps to DISTINCT ids (so a label-keyed viewer would
 *      collide, which is exactly why this one does not).
 * Returns a report; the caller refuses to draw when `ok` is false.
 */
export function assertNodeIdentity(nodes) {
  const problems = [];
  const byId = new Map();
  for (const node of nodes) {
    if (byId.has(node.node_id)) problems.push(`node_id ${node.node_id} appears twice`);
    byId.set(node.node_id, node);
  }
  for (const node of nodes) {
    for (const child of node.children || []) {
      if (!byId.has(child)) problems.push(`node ${node.node_id} names child ${child}, which does not exist`);
    }
    if (node.parent_id !== null && node.parent_id !== undefined && !byId.has(node.parent_id)) {
      problems.push(`node ${node.node_id} names parent ${node.parent_id}, which does not exist`);
    }
  }
  const byLabel = new Map();
  for (const node of nodes) {
    if (node.label === null || node.label === undefined) continue;
    if (!byLabel.has(node.label)) byLabel.set(node.label, []);
    byLabel.get(node.label).push(node.node_id);
  }
  const duplicateLabels = [];
  for (const [label, ids] of byLabel) {
    if (ids.length > 1) duplicateLabels.push({ label, ids });
  }
  const ids = [...byId.keys()];
  return {
    ok: problems.length === 0,
    problems,
    n_nodes: nodes.length,
    n_unique_ids: new Set(ids).size,
    duplicateLabels,
  };
}

export function mount(container, ctx) {
  const abort = new AbortController();
  let disposed = false;
  const disposers = [];
  const controlDisposers = [];
  const observers = [];

  const state = {
    layout: 'rect',
    colorBy: 'sir',
    collapsed: new Set(), // node_id
    selected: null, // node_id
    searchMatch: null, // node_id
    searchTotal: null,
    k: 1,
    tx: 0,
    ty: 0,
    dragging: false,
    moved: false,
    lastX: 0,
    lastY: 0,
    hoverId: null,
  };

  /* ------------------------------------------------------------ chrome */

  const head = h(
    'div',
    { class: 'page-head' },
    h(
      'div',
      { class: 'page-head-text' },
      h('h2', null, 'Phylogeny'),
      h(
        'p',
        { class: 'page-note' },
        'The stage-9 tree as structural JSON. Node identity is the server’s postorder ' +
          'index (UI-D6): selection, collapse and zoom key on node_id, never on the ' +
          'Newick label. Branch lengths are substitutions per site, not SNP counts.'
      )
    )
  );
  container.appendChild(head);

  const alertHost = h('div');
  const controls = h('div', { class: 'panel' });
  const searchStatus = h('span', { class: 'faint' });
  const canvasWrap = h('div', {
    class: 'panel',
    style: { padding: '0', position: 'relative', overflow: 'hidden', height: '620px' },
  });
  const canvas = h('canvas', { style: { display: 'block', width: '100%', height: '100%' } });
  const tooltip = h('div', {
    style: {
      position: 'absolute',
      'pointer-events': 'none',
      display: 'none',
      'max-width': '26rem',
      'z-index': '5',
      background: 'var(--bg-raised)',
      border: '1px solid var(--border-strong)',
      'border-radius': 'var(--radius)',
      'box-shadow': 'var(--shadow)',
      padding: 'var(--space-2)',
      'font-size': 'var(--fs-xs)',
    },
  });
  canvasWrap.append(canvas, tooltip);
  const legendHost = h('div', { class: 'panel' });
  const metaHost = h('div', { class: 'panel' });
  container.append(alertHost, controls, canvasWrap, legendHost, metaHost);

  /* ------------------------------------------------------------- data */

  let tree = null;
  let nodes = [];
  let byId = new Map(); // node_id -> node   (NEVER keyed by label)
  let rootId = null;
  let isolates = new Map(); // sample_id -> row (a join, not identity)
  let layout = null; // {mode, positions:Map, terminals:[], modelW, modelH, ...}
  let renderMs = null;

  const ctx2d = canvas.getContext('2d');

  function tokens() {
    try {
      return ctx.theme.tokens || {};
    } catch {
      return {};
    }
  }

  function tipMeta(node) {
    return node && node.tip_metadata ? node.tip_metadata : null;
  }

  function rowFor(node) {
    if (!node || !node.label) return null;
    return isolates.get(node.label) || null;
  }

  function tipField(node, key) {
    const row = rowFor(node);
    if (row && row[key] !== undefined) return row[key];
    const meta = tipMeta(node);
    if (meta && meta[key] !== undefined) return meta[key];
    return null;
  }

  /* ------------------------------------------------------- colour map */

  function hashHue(text) {
    let hash = 0;
    const value = String(text);
    for (let i = 0; i < value.length; i += 1) {
      hash = (hash * 31 + value.charCodeAt(i)) % 360;
    }
    return hash;
  }

  function categorical(text) {
    if (text === null || text === undefined || text === '') return tokens()['tone-unknown-fg'] || '#888';
    return `hsl(${hashHue(text)}, 55%, 45%)`;
  }

  function colourFor(node) {
    const t = tokens();
    if (node.is_tip) {
      if (state.colorBy === 'sir') {
        const sir = tipField(node, 'phenotype_sir');
        return t[SIR_TOKENS[sir]] || t['tone-unknown-fg'] || '#888';
      }
      if (state.colorBy === 'oprd') {
        const oprd = tipField(node, 'oprd_state');
        return t[OPRD_TOKENS[oprd]] || t['tone-unknown-fg'] || '#888';
      }
      if (state.colorBy === 'st') return categorical(tipField(node, 'st'));
      if (state.colorBy === 'lineage') return categorical(tipField(node, 'lineage'));
      return t.fg || '#333';
    }
    return t['fg-faint'] || '#999';
  }

  /* ----------------------------------------------------------- layout */

  function buildFullLeafOrder() {
    const result = new Map(); // node_id -> {min,max,mean} leaf index
    const visitOrder = [];
    const walk = [rootId];
    while (walk.length) {
      const id = walk.pop();
      visitOrder.push(id);
      const node = byId.get(id);
      if (node && node.children) {
        for (let i = node.children.length - 1; i >= 0; i -= 1) walk.push(node.children[i]);
      }
    }
    let leafIndex = 0;
    for (const id of visitOrder.reverse()) {
      const node = byId.get(id);
      if (!node) continue;
      if (node.is_tip) {
        result.set(id, { min: leafIndex, max: leafIndex, mean: leafIndex });
        leafIndex += 1;
      } else {
        let min = Infinity;
        let max = -Infinity;
        let sum = 0;
        let count = 0;
        for (const child of node.children) {
          const info = result.get(child);
          if (!info) continue;
          min = Math.min(min, info.min);
          max = Math.max(max, info.max);
          sum += info.mean;
          count += 1;
        }
        result.set(id, {
          min: min === Infinity ? 0 : min,
          max: max === -Infinity ? 0 : max,
          mean: count ? sum / count : 0,
        });
      }
    }
    return { info: result, nLeaves: leafIndex };
  }

  function computeDistances() {
    const dist = new Map();
    const depth = new Map();
    const stack = [{ id: rootId, d: 0, depth: 0 }];
    let maxDist = 0;
    let maxDepth = 0;
    while (stack.length) {
      const { id, d, depth: dep } = stack.pop();
      const node = byId.get(id);
      if (!node) continue;
      const length = typeof node.length === 'number' && isFinite(node.length) ? node.length : 0;
      const total = d + length;
      dist.set(id, total);
      depth.set(id, dep);
      if (total > maxDist) maxDist = total;
      if (dep > maxDepth) maxDepth = dep;
      for (const child of node.children) stack.push({ id: child, d: total, depth: dep + 1 });
    }
    return { dist, depth, maxDist, maxDepth };
  }

  function collectVisibleTerminals() {
    const out = [];
    const stack = [rootId];
    while (stack.length) {
      const id = stack.pop();
      const node = byId.get(id);
      if (!node) continue;
      if (node.is_tip || state.collapsed.has(id)) {
        out.push(id);
        continue;
      }
      for (let i = node.children.length - 1; i >= 0; i -= 1) stack.push(node.children[i]);
    }
    return out;
  }

  function computeLayout() {
    if (!tree || !rootId) return null;
    const { info, nLeaves } = buildFullLeafOrder();
    const { dist, depth, maxDist, maxDepth } = computeDistances();
    const useBranchLengths = tree.has_branch_lengths && maxDist > 0;
    const terminals = collectVisibleTerminals();
    const slotOf = new Map();
    terminals.forEach((id, index) => slotOf.set(id, index));

    const positions = new Map();
    const stack = [{ id: rootId, visited: false }];
    const postorder = [];
    while (stack.length) {
      const frame = stack.pop();
      if (!frame.visited) {
        stack.push({ id: frame.id, visited: true });
        const node = byId.get(frame.id);
        if (node && !node.is_tip && !state.collapsed.has(frame.id)) {
          for (const child of node.children) stack.push({ id: child, visited: false });
        }
      } else {
        postorder.push(frame.id);
      }
    }
    for (const id of postorder) {
      const node = byId.get(id);
      if (!node) continue;
      const terminal = node.is_tip || state.collapsed.has(id);
      const info1 = info.get(id) || { mean: 0 };
      if (terminal) {
        positions.set(id, { y: slotOf.get(id) !== undefined ? slotOf.get(id) : info1.mean, x: useBranchLengths ? dist.get(id) : depth.get(id) });
      } else {
        let sum = 0;
        let count = 0;
        for (const child of node.children) {
          const p = positions.get(child);
          if (p) {
            sum += p.y;
            count += 1;
          }
        }
        positions.set(id, { y: count ? sum / count : info1.mean, x: useBranchLengths ? dist.get(id) : depth.get(id) });
      }
    }

    const xMax = useBranchLengths ? maxDist : maxDepth || 1;
    const nTerminals = Math.max(1, terminals.length);
    const modelW = X_SPAN + MARGIN * 2;
    const modelH = nTerminals * ROW_H + MARGIN * 2;
    const model = new Map();
    for (const [id, p] of positions) {
      model.set(id, {
        x: MARGIN + (p.x / xMax) * X_SPAN,
        y: MARGIN + (p.y + 0.5) * ROW_H,
      });
    }
    return {
      mode: state.layout,
      positions: model,
      raw: positions,
      info,
      terminals,
      nLeaves,
      nTerminals,
      modelW,
      modelH,
      useBranchLengths,
      maxDist,
      maxDepth,
      radius: 0,
    };
  }

  function computeCircular() {
    const base = computeLayout();
    if (!base) return null;
    const nLeaves = Math.max(1, base.nLeaves);
    const radius = 360;
    const model = new Map();
    for (const [id, p] of base.raw) {
      const angle = (p.y + 0.5) / nLeaves * Math.PI * 2 - Math.PI / 2;
      const r = MARGIN + (base.useBranchLengths
        ? (base.maxDist > 0 ? p.x / base.maxDist : 0)
        : (base.maxDepth > 0 ? p.x / base.maxDepth : 0)) * radius;
      model.set(id, { x: 400 + r * Math.cos(angle), y: 400 + r * Math.sin(angle), angle, r });
    }
    return Object.assign({}, base, {
      mode: 'circ',
      positions: model,
      modelW: 800,
      modelH: 800,
      radius,
    });
  }

  function relayout() {
    layout = state.layout === 'circ' ? computeCircular() : computeLayout();
    fit();
    draw();
    updateMeta();
  }

  /* ------------------------------------------------------- transforms */

  function viewSize() {
    const rect = canvasWrap.getBoundingClientRect
      ? canvasWrap.getBoundingClientRect()
      : { width: 900, height: 640 };
    return {
      width: Math.max(320, Math.floor(rect.width || 900)),
      height: Math.max(360, Math.floor(rect.height || 640)),
    };
  }

  function fit() {
    if (!layout) return;
    const { width, height } = viewSize();
    const k = Math.min(width / layout.modelW, height / layout.modelH);
    state.k = k > 0 && isFinite(k) ? k : 1;
    state.tx = (width - layout.modelW * state.k) / 2;
    state.ty = (height - layout.modelH * state.k) / 2;
  }

  function screenOf(model) {
    return { x: model.x * state.k + state.tx, y: model.y * state.k + state.ty };
  }

  /* ------------------------------------------------------------ draw */

  function resizeCanvas() {
    const { width, height } = viewSize();
    const dpr = (globalThis.devicePixelRatio || 1);
    canvas.width = Math.floor(width * dpr);
    canvas.height = Math.floor(height * dpr);
    canvas.style.width = `${width}px`;
    canvas.style.height = `${height}px`;
    if (ctx2d.setTransform) ctx2d.setTransform(dpr, 0, 0, dpr, 0, 0);
  }

  function draw() {
    if (!layout || !ctx2d) return;
    const started = typeof performance !== 'undefined' ? performance.now() : Date.now();
    resizeCanvas();
    const { width, height } = viewSize();
    const t = tokens();
    ctx2d.clearRect(0, 0, width, height);
    ctx2d.save();
    ctx2d.translate(state.tx, state.ty);
    ctx2d.scale(state.k, state.k);
    ctx2d.lineWidth = Math.max(0.4, 1 / state.k);
    ctx2d.strokeStyle = t.border || '#ccc';
    ctx2d.lineCap = 'round';

    const visible = new Set();
    for (const id of layout.terminals) visible.add(id);
    // every ancestor of a visible terminal is drawn
    for (const id of layout.terminals) {
      let cursor = byId.get(id);
      while (cursor && cursor.parent_id !== null && cursor.parent_id !== undefined) {
        visible.add(cursor.parent_id);
        cursor = byId.get(cursor.parent_id);
      }
    }

    // edges
    ctx2d.beginPath();
    for (const id of visible) {
      const node = byId.get(id);
      if (!node || node.parent_id === null || node.parent_id === undefined) continue;
      if (!visible.has(node.parent_id)) continue;
      const from = layout.positions.get(node.parent_id);
      const to = layout.positions.get(id);
      if (!from || !to) continue;
      if (layout.mode === 'circ') {
        ctx2d.moveTo(from.x, from.y);
        ctx2d.lineTo(to.x, to.y);
      } else {
        ctx2d.moveTo(from.x, from.y);
        ctx2d.lineTo(to.x, from.y);
        ctx2d.lineTo(to.x, to.y);
      }
    }
    ctx2d.stroke();

    // nodes
    const radius = Math.max(0.6, 2.4 / state.k);
    for (const id of visible) {
      const node = byId.get(id);
      const p = layout.positions.get(id);
      if (!node || !p) continue;
      ctx2d.beginPath();
      ctx2d.fillStyle = colourFor(node);
      if (node.is_tip) {
        ctx2d.arc(p.x, p.y, radius, 0, Math.PI * 2);
      } else {
        ctx2d.rect(p.x - radius, p.y - radius, radius * 2, radius * 2);
      }
      ctx2d.fill();
      if (state.selected === id) {
        ctx2d.beginPath();
        ctx2d.strokeStyle = t.accent || '#06c';
        ctx2d.lineWidth = Math.max(1, 2 / state.k);
        ctx2d.arc(p.x, p.y, radius + 3 / state.k, 0, Math.PI * 2);
        ctx2d.stroke();
      }
      if (state.searchMatch === id) {
        ctx2d.beginPath();
        ctx2d.strokeStyle = t['tone-warn-fg'] || '#a80';
        ctx2d.lineWidth = Math.max(1, 2 / state.k);
        ctx2d.arc(p.x, p.y, radius + 6 / state.k, 0, Math.PI * 2);
        ctx2d.stroke();
      }
    }

    // labels: only when there is room, or for the selected/searched tip
    const showAllLabels = ROW_H * state.k >= LABEL_MIN_PX;
    if (showAllLabels || state.selected !== null || state.searchMatch !== null) {
      ctx2d.font = `${Math.max(7, 10 / state.k)}px ${t.font || 'sans-serif'}`;
      ctx2d.textBaseline = 'middle';
      ctx2d.fillStyle = t.fg || '#111';
      for (const id of visible) {
        const node = byId.get(id);
        const p = layout.positions.get(id);
        if (!node || !p || !node.is_tip) continue;
        const wanted = showAllLabels || state.selected === id || state.searchMatch === id;
        if (!wanted) continue;
        const text = node.label || '';
        if (!text) continue;
        if (layout.mode === 'circ') {
          const angle = p.angle || 0;
          ctx2d.textAlign = Math.cos(angle) >= 0 ? 'left' : 'right';
          ctx2d.fillText(text, p.x + (Math.cos(angle) >= 0 ? 4 / state.k : -4 / state.k), p.y);
        } else {
          ctx2d.textAlign = 'left';
          ctx2d.fillText(text, p.x + 4 / state.k, p.y);
        }
      }
    }
    ctx2d.restore();
    renderMs = (typeof performance !== 'undefined' ? performance.now() : Date.now()) - started;
    try {
      container.dataset.treeRenderMs = String(Math.round(renderMs));
    } catch {
      /* dataset is best-effort */
    }
  }

  /* -------------------------------------------------------- hit test */

  function hitTest(sx, sy) {
    if (!layout) return null;
    let best = null;
    let bestDist = 12;
    for (const [id, p] of layout.positions) {
      const s = screenOf(p);
      const dx = s.x - sx;
      const dy = s.y - sy;
      const distance = Math.sqrt(dx * dx + dy * dy);
      if (distance < bestDist) {
        bestDist = distance;
        best = id;
      }
    }
    return best;
  }

  /* --------------------------------------------------------- tooltip */

  function showTooltip(node, sx, sy) {
    if (!node || !node.is_tip) {
      tooltip.style.display = 'none';
      return;
    }
    clear(tooltip);
    tooltip.appendChild(h('strong', null, node.label || '(unlabelled tip)'));
    const rows = [
      ['lineage', tipField(node, 'lineage')],
      ['ST', tipField(node, 'st')],
      ['imipenem SIR', tipField(node, 'phenotype_sir')],
      ['oprD state', tipField(node, 'oprd_state')],
      ['oprD verdict', tipField(node, 'oprd_verdict')],
      ['virulence factors', tipField(node, 'n_virulence')],
    ];
    tooltip.appendChild(kv(rows.map(([k, v]) => [
      k,
      v === null || v === undefined || v === ''
        ? absentSpan('the joined isolate table records no value for this field', ABSENT.NOT_ASSESSED)
        : document.createTextNode(String(v)),
    ])));
    const s = screenOf(layout.positions.get(node.node_id) || { x: 0, y: 0 });
    tooltip.style.display = 'block';
    tooltip.style.left = `${Math.round(s.x + 12)}px`;
    tooltip.style.top = `${Math.round(s.y + 12)}px`;
    void sx;
    void sy;
  }

  function hideTooltip() {
    tooltip.style.display = 'none';
  }

  /* ---------------------------------------------------------- legend */

  function renderLegend() {
    clear(legendHost);
    legendHost.appendChild(h('h3', null, `Colour by ${state.colorBy}`));
    const row = h('div', { class: 'chip-row' });
    const t = tokens();
    const add = (colour, label) => {
      const swatch = h('span', {
        class: 'chip',
        style: { 'border-left': `0.9rem solid ${colour}`, 'border-left-width': '0.9rem' },
      }, label);
      row.appendChild(swatch);
    };
    if (state.colorBy === 'sir') {
      add(t['tone-bad-fg'] || '#a00', 'R — resistant');
      add(t['tone-warn-fg'] || '#a80', 'I — intermediate');
      add(t['tone-ok-fg'] || '#080', 'S — susceptible');
      add(t['tone-busy-fg'] || '#088', 'SDD');
      add(t['tone-unknown-fg'] || '#888', 'ND / not assessed');
    } else if (state.colorBy === 'oprd') {
      add(t['tone-ok-fg'] || '#080', 'intact');
      add(t['tone-bad-fg'] || '#a00', 'disrupted');
      add(t['tone-warn-fg'] || '#a80', 'absent');
      add(t['tone-unknown-fg'] || '#888', 'not_assessed (a refusal maps here, never absent)');
    } else {
      const values = new Set();
      for (const id of layout ? layout.terminals : []) {
        const node = byId.get(id);
        if (!node || !node.is_tip) continue;
        const v = state.colorBy === 'st' ? tipField(node, 'st') : tipField(node, 'lineage');
        if (v !== null && v !== undefined && v !== '') values.add(String(v));
      }
      if (values.size === 0) {
        row.appendChild(absentSpan('no value for this field is recorded for any tip', ABSENT.NOT_ASSESSED));
      } else {
        for (const value of [...values].sort().slice(0, 24)) add(categorical(value), value);
        if (values.size > 24) row.appendChild(h('span', { class: 'faint' }, `and ${values.size - 24} more`));
      }
    }
    legendHost.appendChild(row);
    legendHost.appendChild(
      note(
        'Status is carried by colour AND a word, never by hue alone. The categorical ' +
          'palette is derived in the browser because the token palette names the six ' +
          'status tones and not an arbitrary number of STs.',
        'faint'
      )
    );
  }

  /* ------------------------------------------------------------ meta */

  function updateMeta() {
    clear(metaHost);
    metaHost.appendChild(h('h3', null, 'What this view is'));
    const identity = tree && tree.identity ? tree.identity : '(the server did not send an identity field)';
    const pairs = [
      ['identity', h('code', null, identity)],
      ['tips', h('b', null, String(tree ? tree.n_tips : 0))],
      ['internal nodes', h('b', null, String(tree ? tree.n_internal : 0))],
      ['branch lengths', tree && tree.has_branch_lengths ? 'present — distances are summed substitutions per site' : 'absent — the horizontal axis is node depth, not a distance'],
      ['layout', state.layout === 'circ' ? 'circular' : 'rectangular'],
      ['collapsed clades', h('b', null, String(state.collapsed.size))],
      ['render time', renderMs === null ? absentSpan('not rendered yet') : h('b', null, `${renderMs.toFixed(1)} ms for ${layout ? layout.positions.size : 0} nodes`)],
    ];
    metaHost.appendChild(kv(pairs));
    if (tree && tree.tip_check) {
      const tc = tree.tip_check;
      const box = h('div');
      box.appendChild(h('h4', null, 'Tip join check (phylogeny.validate_tree_samples)'));
      if (tc.matched) {
        box.appendChild(note('the tip set equals the cohort set; every tip joins to an isolate by sample_id'));
      } else {
        box.appendChild(
          message(
            'warn',
            'This tree’s tips cannot be joined to the isolates',
            serverMessage(tc.reason || 'the server reported a set difference'),
            note(`missing from tree: ${tc.n_missing_in_tree}; extra in tree: ${tc.n_extra_in_tree}`)
          )
        );
      }
      metaHost.appendChild(box);
    }
    if (tree && tree.metadata_reason) {
      metaHost.appendChild(note(tree.metadata_reason, 'faint'));
    }
    const ident = assertNodeIdentity(nodes);
    metaHost.appendChild(h('h4', null, 'Node identity assertion (UI-D6)'));
    if (ident.ok) {
      metaHost.appendChild(
        note(
          `${ident.n_nodes} nodes, ${ident.n_unique_ids} unique node_ids. ` +
            `${ident.duplicateLabels.length} internal label(s) repeat ` +
            `(${ident.duplicateLabels.slice(0, 4).map((d) => `${JSON.stringify(d.label)} ×${d.ids.length}`).join(', ') || 'none'}); ` +
            'each repeated label resolves to distinct node_ids, and this view keys on node_id only.'
        )
      );
    } else {
      metaHost.appendChild(
        message('error', 'The server’s structural JSON failed the node-identity assertion', note(ident.problems.join('; ')))
      );
    }
  }

  /* --------------------------------------------------------- controls */

  function renderControls() {
    // The controls are rebuilt on every selection/collapse; release the previous
    // build's listeners first so a long session does not accumulate closures on
    // detached inputs.
    for (const fn of controlDisposers.splice(0)) {
      try {
        fn();
      } catch {
        /* already gone */
      }
    }
    clear(controls);
    controls.appendChild(h('h3', null, 'View'));
    const row = h('div', { class: 'toolbar' });

    const mkToggle = (label, active, onClick) => {
      const button = h('button', { type: 'button', class: active ? 'primary' : null }, label);
      button.addEventListener('click', onClick);
      controlDisposers.push(() => button.removeEventListener('click', onClick));
      return button;
    };
    row.appendChild(mkToggle('rectangular', state.layout === 'rect', () => {
      state.layout = 'rect';
      relayout();
      renderControls();
    }));
    row.appendChild(mkToggle('circular', state.layout === 'circ', () => {
      state.layout = 'circ';
      relayout();
      renderControls();
    }));

    const colour = h('select', { 'aria-label': 'colour tips by' });
    for (const [value, label] of [
      ['sir', 'imipenem SIR'],
      ['st', 'sequence type (ST)'],
      ['lineage', 'lineage'],
      ['oprd', 'oprD verdict'],
      ['none', 'no colour'],
    ]) {
      colour.appendChild(h('option', { value, selected: state.colorBy === value ? true : null }, label));
    }
    const onColour = () => {
      state.colorBy = colour.value;
      draw();
      renderLegend();
    };
    colour.addEventListener('change', onColour);
    controlDisposers.push(() => colour.removeEventListener('change', onColour));
    row.appendChild(h('label', null, 'colour', colour));

    const search = h('input', {
      type: 'search',
      placeholder: 'search a tip (server-side)',
      'aria-label': 'search a tip',
      spellcheck: false,
      autocomplete: 'off',
    });
    const onSearch = debounce(async () => {
      const q = search.value.trim();
      if (!q) {
        state.searchMatch = null;
        state.searchTotal = null;
        searchStatus.textContent = '';
        draw();
        return;
      }
      searchStatus.textContent = 'searching…';
      try {
        const body = await ctx.api.get('/tree/tips', { q, limit: 1 }, abort.signal);
        if (disposed) return;
        const first = body.items && body.items[0];
        state.searchMatch = first ? first.node_id : null;
        state.searchTotal = body.meta ? body.meta.total : null;
        if (state.searchMatch !== null) {
          state.selected = state.searchMatch;
          centreOn(state.searchMatch);
          searchStatus.textContent = `${state.searchTotal} tip(s) match (the server’s count over the whole tree); the first is selected`;
        } else {
          searchStatus.textContent = `0 tips match “${q}” (the server’s count over the whole tree)`;
        }
        draw();
      } catch (error) {
        if (error && error.name === 'AbortError') return;
        state.searchTotal = null;
        searchStatus.textContent = `the search could not be run: ${error && error.message ? error.message : error}`;
      }
    }, 250);
    search.addEventListener('input', onSearch);
    controlDisposers.push(() => search.removeEventListener('input', onSearch));
    row.appendChild(h('label', null, 'find tip', search));
    row.appendChild(searchStatus);

    const mkButton = (label, onClick) => {
      const button = h('button', { type: 'button' }, label);
      button.addEventListener('click', onClick);
      controlDisposers.push(() => button.removeEventListener('click', onClick));
      return button;
    };
    row.appendChild(mkButton('fit', () => {
      fit();
      draw();
    }));
    row.appendChild(mkButton('expand every clade', () => {
      state.collapsed.clear();
      relayout();
    }));
    if (state.selected !== null) {
      const node = byId.get(state.selected);
      if (node && !node.is_tip) {
        row.appendChild(mkButton(
          state.collapsed.has(state.selected) ? 'expand selected clade' : 'collapse selected clade',
          () => toggleCollapse(state.selected)
        ));
      }
    }
    controls.appendChild(row);
    controls.appendChild(
      note(
        'Click a tip to select it and show its isolate summary; click an internal node to ' +
          'collapse or expand its clade. Drag to pan, scroll to zoom. The zoom is hand-written ' +
          'because the vendored d3 set omits d3-dispatch/d3-drag/d3-transition, which d3-zoom needs.'
      )
    );
  }

  function centreOn(id) {
    const p = layout && layout.positions.get(id);
    if (!p) return;
    const { width, height } = viewSize();
    state.tx = width / 2 - p.x * state.k;
    state.ty = height / 2 - p.y * state.k;
  }

  function toggleCollapse(id) {
    if (state.collapsed.has(id)) state.collapsed.delete(id);
    else state.collapsed.add(id);
    relayout();
    renderControls();
  }

  /* --------------------------------------------------------- events */

  function onWheel(event) {
    event.preventDefault();
    const rect = canvas.getBoundingClientRect();
    const sx = event.clientX - rect.left;
    const sy = event.clientY - rect.top;
    const factor = event.deltaY < 0 ? 1.15 : 1 / 1.15;
    const nextK = Math.min(4000, Math.max(0.02, state.k * factor));
    const mx = (sx - state.tx) / state.k;
    const my = (sy - state.ty) / state.k;
    state.k = nextK;
    state.tx = sx - mx * nextK;
    state.ty = sy - my * nextK;
    draw();
  }

  function onPointerDown(event) {
    state.dragging = true;
    state.moved = false;
    state.lastX = event.clientX;
    state.lastY = event.clientY;
    if (canvas.setPointerCapture && event.pointerId !== undefined) {
      try {
        canvas.setPointerCapture(event.pointerId);
      } catch {
        /* not fatal */
      }
    }
  }

  function onPointerMove(event) {
    const rect = canvas.getBoundingClientRect();
    const sx = event.clientX - rect.left;
    const sy = event.clientY - rect.top;
    if (state.dragging) {
      const dx = event.clientX - state.lastX;
      const dy = event.clientY - state.lastY;
      if (Math.abs(dx) + Math.abs(dy) > 3) state.moved = true;
      state.tx += dx;
      state.ty += dy;
      state.lastX = event.clientX;
      state.lastY = event.clientY;
      draw();
      return;
    }
    const id = hitTest(sx, sy);
    if (id !== state.hoverId) {
      state.hoverId = id;
      const node = id === null ? null : byId.get(id);
      if (node && node.is_tip) showTooltip(node, sx, sy);
      else hideTooltip();
    } else if (id !== null) {
      const node = byId.get(id);
      if (node && node.is_tip) showTooltip(node, sx, sy);
    }
  }

  function onPointerUp() {
    state.dragging = false;
  }

  function onClick(event) {
    if (state.moved) return;
    const rect = canvas.getBoundingClientRect();
    const id = hitTest(event.clientX - rect.left, event.clientY - rect.top);
    if (id === null) {
      state.selected = null;
      hideTooltip();
      draw();
      return;
    }
    const node = byId.get(id);
    state.selected = id;
    if (node && !node.is_tip) toggleCollapse(id);
    else draw();
    renderControls();
  }

  function onPointerLeave() {
    state.dragging = false;
    hideTooltip();
  }

  function onResize() {
    draw();
  }

  canvas.addEventListener('wheel', onWheel, { passive: false });
  canvas.addEventListener('pointerdown', onPointerDown);
  canvas.addEventListener('pointermove', onPointerMove);
  canvas.addEventListener('pointerup', onPointerUp);
  canvas.addEventListener('click', onClick);
  canvas.addEventListener('pointerleave', onPointerLeave);
  disposers.push(() => {
    canvas.removeEventListener('wheel', onWheel);
    canvas.removeEventListener('pointerdown', onPointerDown);
    canvas.removeEventListener('pointermove', onPointerMove);
    canvas.removeEventListener('pointerup', onPointerUp);
    canvas.removeEventListener('click', onClick);
    canvas.removeEventListener('pointerleave', onPointerLeave);
  });

  if (typeof ResizeObserver === 'function') {
    const observer = new ResizeObserver(() => onResize());
    observer.observe(canvasWrap);
    observers.push(observer);
  } else if (typeof window !== 'undefined' && window.addEventListener) {
    window.addEventListener('resize', onResize);
    disposers.push(() => window.removeEventListener('resize', onResize));
  }

  if (typeof MutationObserver === 'function' && typeof document !== 'undefined') {
    const observer = new MutationObserver(() => draw());
    observer.observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
    observers.push(observer);
  }
  if (typeof window !== 'undefined' && window.matchMedia) {
    const query = window.matchMedia('(prefers-color-scheme: dark)');
    if (query.addEventListener) {
      const handler = () => draw();
      query.addEventListener('change', handler);
      disposers.push(() => query.removeEventListener('change', handler));
    }
  }

  /* ------------------------------------------------------------ load */

  (async () => {
    let treeBody;
    let isolateBody = null;
    try {
      const [treeResult, isolateResult] = await Promise.all([
        ctx.api.get('/tree', { source: 'stage9' }, abort.signal),
        ctx.api.get(
          '/isolates',
          { limit: 1000, columns: 'sample_id,phenotype_sir,lineage,st,oprd_state,oprd_verdict,n_virulence' },
          abort.signal
        ).catch((error) => (error && error.name === 'AbortError' ? null : null)),
      ]);
      treeBody = treeResult;
      isolateBody = isolateResult;
    } catch (error) {
      if (error && error.name === 'AbortError') return;
      alertHost.appendChild(
        message('error', 'The tree could not be read', note(String((error && error.message) || error)))
      );
      return;
    }
    if (disposed) return;

    if (!treeBody || treeBody.present === false) {
      alertHost.appendChild(
        message(
          'info',
          'This tree is not produced',
          serverMessage(treeBody && treeBody.reason ? treeBody.reason : 'the server gave no reason'),
          treeBody && treeBody.probes ? probeList(treeBody.probes) : null
        )
      );
      return;
    }

    tree = treeBody;
    nodes = Array.isArray(tree.nodes) ? tree.nodes : [];
    byId = new Map();
    for (const node of nodes) byId.set(node.node_id, node);
    const roots = nodes.filter((n) => n.parent_id === null || n.parent_id === undefined);
    rootId = roots.length ? roots[0].node_id : (nodes.length ? nodes[nodes.length - 1].node_id : null);

    if (isolateBody && Array.isArray(isolateBody.items)) {
      for (const row of isolateBody.items) {
        if (row && row.sample_id !== undefined) isolates.set(row.sample_id, row);
      }
    }

    const identity = assertNodeIdentity(nodes);
    if (!identity.ok) {
      alertHost.appendChild(
        message(
          'error',
          'The server’s structural JSON failed the node-identity assertion; nothing is drawn',
          note(identity.problems.join('; '))
        )
      );
      updateMeta();
      return;
    }

    layout = state.layout === 'circ' ? computeCircular() : computeLayout();
    fit();
    renderControls();
    renderLegend();
    draw();
    updateMeta();
  })();

  /* -------------------------------------------------------- teardown */

  return () => {
    disposed = true;
    abort.abort();
    hideTooltip();
    for (const observer of observers.splice(0)) {
      try {
        observer.disconnect();
      } catch {
        /* already gone */
      }
    }
    for (const fn of controlDisposers.splice(0)) {
      try {
        fn();
      } catch (error) {
        console.error('[tree] a control disposer threw', error);
      }
    }
    for (const fn of disposers.splice(0)) {
      try {
        fn();
      } catch (error) {
        console.error('[tree] a disposer threw', error);
      }
    }
  };
}

export default { mount };
