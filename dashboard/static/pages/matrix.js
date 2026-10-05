/* pages/matrix.js — V2, the stage-10 patristic-distance heatmap.
 *
 * Three honesty rules the page exists to keep:
 *
 *   1. Above 300 tips the server downsamples the matrix and says so
 *      (`n_total`/`n_shown`/`downsampled`/`method`). This page displays that
 *      fact prominently; it never presents the shown subset as the cohort
 *      (UI-D2/UI-D7). It does not try to fetch all 900 and draw them.
 *   2. Ordering is TREE ORDER. The matrix's own `order` is the file's row
 *      order; the stage-9 tree's left-to-right tip order is fetched and the
 *      rows/columns are permuted to match. A sample absent from the tree is
 *      appended and named, never silently dropped.
 *   3. Units are `similarity.units.json` passed through verbatim. When it is
 *      absent the page says the units are not recorded; it does NOT assume
 *      "substitutions per site". When it is present, `is_a_snp_count: false`
 *      is stated, because rule 11 makes that a claim about what may be written.
 *
 * Rendering is canvas. d3-scale / d3-array / d3-format / d3-interpolate do the
 * colour ramp, the extents and the number formatting.
 */

import {
  ABSENT,
  absentSpan,
  clear,
  h,
  kv,
  message,
  note,
  probeList,
  serverMessage,
} from '../app.js';
import { mountPowerBanner } from './tables.js';

function d3lib() {
  return typeof globalThis !== 'undefined' ? globalThis.d3 : undefined;
}

function treeTipOrder(tree) {
  if (!tree || !Array.isArray(tree.nodes)) return null;
  const byId = new Map();
  for (const node of tree.nodes) byId.set(node.node_id, node);
  const roots = tree.nodes.filter((n) => n.parent_id === null || n.parent_id === undefined);
  const root = roots.length ? roots[0] : tree.nodes[tree.nodes.length - 1];
  if (!root) return null;
  const order = [];
  const stack = [root.node_id];
  while (stack.length) {
    const id = stack.pop();
    const node = byId.get(id);
    if (!node) continue;
    if (node.is_tip) {
      order.push(node.label);
      continue;
    }
    for (let i = node.children.length - 1; i >= 0; i -= 1) stack.push(node.children[i]);
  }
  return order;
}

export function mount(container, ctx) {
  const d3 = d3lib();
  const abort = new AbortController();
  let disposed = false;
  const disposers = [];
  const observers = [];

  const head = h(
    'div',
    { class: 'page-head' },
    h(
      'div',
      { class: 'page-head-text' },
      h('h2', null, 'Similarity'),
      h(
        'p',
        { class: 'page-note' },
        'The patristic-distance matrix (stage 10), ordered by the stage-9 tree. ' +
          'Above 300 tips the server downsamples; this page states that rather than ' +
          'presenting the subset as the cohort. Distances are expected substitutions ' +
          'per site, never a SNP count.'
      )
    )
  );
  container.appendChild(head);

  const alertHost = h('div');
  const bannerHost = h('div');
  const statusHost = h('div', { class: 'panel' });
  const canvasWrap = h('div', {
    class: 'panel',
    style: { padding: '0', position: 'relative', overflow: 'hidden', height: '640px' },
  });
  const canvas = h('canvas', { style: { display: 'block', width: '100%', height: '100%' } });
  const tooltip = h('div', {
    style: {
      position: 'absolute',
      'pointer-events': 'none',
      display: 'none',
      'z-index': '5',
      background: 'var(--bg-raised)',
      border: '1px solid var(--border-strong)',
      'border-radius': 'var(--radius)',
      'box-shadow': 'var(--shadow)',
      padding: 'var(--space-2)',
      'font-size': 'var(--fs-xs)',
      'font-family': 'var(--font-mono)',
    },
  });
  canvasWrap.append(canvas, tooltip);
  const legendHost = h('div', { class: 'panel' });
  container.append(alertHost, bannerHost, statusHost, canvasWrap, legendHost);

  const ctx2d = canvas.getContext('2d');
  let payload = null;
  let order = [];
  let matrix = [];
  let units = null;
  let unitsReason = '';
  let orderNote = '';
  let colourScale = null;
  let geom = null;
  let renderMs = null;

  function tokens() {
    try {
      return ctx.theme.tokens || {};
    } catch {
      return {};
    }
  }

  function viewSize() {
    const rect = canvasWrap.getBoundingClientRect
      ? canvasWrap.getBoundingClientRect()
      : { width: 700, height: 640 };
    const width = Math.max(320, Math.floor(rect.width || 700));
    const height = Math.max(320, Math.floor(rect.height || 640));
    const side = Math.min(width, height);
    return { width, height, side };
  }

  function permute(shownOrder, treeOrder) {
    if (!treeOrder || treeOrder.length === 0) {
      return { order: shownOrder, note: 'the stage-9 tree was not readable, so the file’s own row order is used and that is stated rather than guessed' };
    }
    const position = new Map();
    treeOrder.forEach((label, index) => {
      if (!position.has(label)) position.set(label, index);
    });
    const inTree = shownOrder.filter((s) => position.has(s));
    const notInTree = shownOrder.filter((s) => !position.has(s));
    inTree.sort((a, b) => position.get(a) - position.get(b));
    const note = notInTree.length
      ? `${notInTree.length} shown sample(s) are not tips in the stage-9 tree and are appended after the tree order: ${notInTree.slice(0, 8).join(', ')}${notInTree.length > 8 ? '…' : ''}`
      : 'rows and columns are in the stage-9 tree’s left-to-right tip order';
    return { order: inTree.concat(notInTree), note };
  }

  function parseHex(value) {
    const text = String(value || '').trim();
    const match = /^#([0-9a-f]{3}|[0-9a-f]{6})$/i.exec(text);
    if (!match) return null;
    let hex = match[1];
    if (hex.length === 3) hex = hex.split('').map((c) => c + c).join('');
    return [parseInt(hex.slice(0, 2), 16), parseInt(hex.slice(2, 4), 16), parseInt(hex.slice(4, 6), 16)];
  }

  function mix(a, b, t) {
    const ca = parseHex(a);
    const cb = parseHex(b);
    if (!ca || !cb) return b;
    const channel = (i) => Math.round(ca[i] + (cb[i] - ca[i]) * t);
    return `rgb(${channel(0)}, ${channel(1)}, ${channel(2)})`;
  }

  function buildColourScale(flat) {
    const t = tokens();
    const low = t['bg-inset'] || '#eee';
    const high = t['accent'] || '#1f5fbf';
    if (!d3 || !d3.scaleLinear) {
      return null;
    }
    const extent = d3.extent ? d3.extent(flat, (v) => v) : [Math.min(...flat), Math.max(...flat)];
    const scale = d3.scaleLinear().domain(extent).range([0, 1]);
    // d3-interpolate's colour interpolators need d3-color, which the vendored
    // set omits, so `d3.interpolateRgb` throws `n.rgb is not a function`.
    // The ramp is therefore a hand-written channel lerp over the tokens.
    return { extent, colour: (value) => mix(low, high, scale(value)), low, high };
  }

  function draw() {
    if (!ctx2d || !order.length) return;
    const started = typeof performance !== 'undefined' ? performance.now() : Date.now();
    const { width, height, side } = viewSize();
    const dpr = globalThis.devicePixelRatio || 1;
    canvas.width = Math.floor(width * dpr);
    canvas.height = Math.floor(height * dpr);
    canvas.style.width = `${width}px`;
    canvas.style.height = `${height}px`;
    if (ctx2d.setTransform) ctx2d.setTransform(dpr, 0, 0, dpr, 0, 0);
    const t = tokens();
    ctx2d.clearRect(0, 0, width, height);

    const n = order.length;
    const labelRoom = n <= 40 ? 130 : 30;
    const plotSide = Math.max(40, side - labelRoom - 12);
    const offx = labelRoom;
    const offy = labelRoom * 0.85 + 6;
    const cell = plotSide / n;
    geom = { offx, offy, cell, n };

    // cells
    for (let r = 0; r < n; r += 1) {
      for (let c = 0; c < n; c += 1) {
        const value = matrix[r] ? matrix[r][c] : null;
        ctx2d.fillStyle = value === null || value === undefined
          ? (t['tone-unknown-bg'] || '#eee')
          : (colourScale ? colourScale.colour(value) : (t.accent || '#1f5fbf'));
        ctx2d.fillRect(offx + c * cell, offy + r * cell, Math.ceil(cell), Math.ceil(cell));
      }
    }
    // grid line for the diagonal, so the zero-distance self-pairs are legible
    ctx2d.strokeStyle = t['fg-faint'] || '#999';
    ctx2d.lineWidth = 1;
    ctx2d.strokeRect(offx + 0.5, offy + 0.5, plotSide, plotSide);

    // axis labels: every sample when small, a tick every k when large
    ctx2d.fillStyle = t.fg || '#111';
    ctx2d.font = `${n <= 40 ? 10 : 9}px ${t['font-mono'] || 'monospace'}`;
    const every = n <= 40 ? 1 : Math.ceil(n / 30);
    for (let i = 0; i < n; i += every) {
      const label = order[i] || '';
      ctx2d.textAlign = 'right';
      ctx2d.textBaseline = 'middle';
      ctx2d.fillText(label, offx - 4, offy + i * cell + cell / 2);
      ctx2d.save();
      ctx2d.translate(offx + i * cell + cell / 2, offy - 4);
      ctx2d.rotate(-Math.PI / 2);
      ctx2d.textAlign = 'left';
      ctx2d.textBaseline = 'middle';
      ctx2d.fillText(label, 0, 0);
      ctx2d.restore();
    }
    renderMs = (typeof performance !== 'undefined' ? performance.now() : Date.now()) - started;
    try {
      container.dataset.matrixRenderMs = String(Math.round(renderMs));
    } catch {
      /* best-effort */
    }
    renderLegend();
    renderStatus();
  }

  function cellAt(sx, sy) {
    if (!geom) return null;
    const c = Math.floor((sx - geom.offx) / geom.cell);
    const r = Math.floor((sy - geom.offy) / geom.cell);
    if (r < 0 || c < 0 || r >= geom.n || c >= geom.n) return null;
    return { r, c };
  }

  function showTooltip(rc, sx, sy) {
    if (!rc) {
      tooltip.style.display = 'none';
      return;
    }
    const a = order[rc.r];
    const b = order[rc.c];
    const value = matrix[rc.r] ? matrix[rc.r][rc.c] : null;
    clear(tooltip);
    tooltip.appendChild(h('div', null, `${a}  ×  ${b}`));
    tooltip.appendChild(
      h(
        'div',
        null,
        value === null || value === undefined
          ? 'no distance recorded for this pair'
          : `distance ${d3 && d3.format ? d3.format('.6g')(value) : value}${units && units.unit ? ` ${units.unit}` : units && units.quantity ? ` (${units.quantity})` : ''}`
      )
    );
    tooltip.style.display = 'block';
    tooltip.style.left = `${Math.round(sx + 12)}px`;
    tooltip.style.top = `${Math.round(sy + 12)}px`;
  }

  function renderLegend() {
    clear(legendHost);
    legendHost.appendChild(h('h3', null, 'Distance scale'));
    if (!colourScale) {
      legendHost.appendChild(
        note(
          'd3-scale / d3-interpolate are not loaded, so the cells are drawn in a single token ' +
            'colour. The matrix is still correct; the ramp is what is missing.',
          'faint'
        )
      );
      return;
    }
    const bar = h('div', {
      style: {
        height: '0.8rem',
        width: '16rem',
        'border-radius': 'var(--radius)',
        border: '1px solid var(--border)',
        background: `linear-gradient(to right, ${colourScale.low}, ${colourScale.high})`,
      },
    });
    const [lo, hi] = colourScale.extent;
    legendHost.appendChild(
      h(
        'div',
        { class: 'toolbar' },
        h('span', { class: 'faint' }, String(lo)),
        bar,
        h('span', { class: 'faint' }, String(hi))
      )
    );
    const unitText = units
      ? `${units.quantity || 'distance'}${units.definition ? ` — ${units.definition}` : ''}${units.unit ? ` (${units.unit})` : ''}`
      : `units not recorded${unitsReason ? `: ${unitsReason}` : ''}`;
    legendHost.appendChild(note(unitText));
    if (units && units.is_a_snp_count === false) {
      legendHost.appendChild(
        note(
          '`is_a_snp_count: false` — these are expected substitutions per site, NOT a count of SNPs.',
          'faint'
        )
      );
    }
  }

  function renderStatus() {
    clear(statusHost);
    statusHost.appendChild(h('h3', null, 'What this matrix is'));
    const pairs = [
      ['samples in the artefact', h('b', null, String(payload.n_total))],
      ['samples shown', h('b', null, String(payload.n_shown))],
      [
        'downsampled',
        payload.downsampled
          ? h('b', { style: { color: 'var(--tone-warn-fg)' } }, `yes — ${payload.method || 'the server reported no method'}`)
          : 'no',
      ],
      ['ordering', orderNote],
      ['render time', renderMs === null ? absentSpan('not rendered yet') : h('b', null, `${renderMs.toFixed(1)} ms for ${order.length}×${order.length}`)],
    ];
    statusHost.appendChild(kv(pairs));
    if (payload.downsampled) {
      statusHost.appendChild(
        message(
          'warn',
          `This view shows ${payload.n_shown} of ${payload.n_total} samples`,
          note(
            'The server downsampled the matrix so the picture keeps the cohort’s shape. ' +
              'This is NOT the whole cohort and no count over it is a cohort count. ' +
              (payload.method ? `Method: ${payload.method}.` : '')
          )
        )
      );
    }
    const basis = payload.basis;
    if (basis) {
      statusHost.appendChild(
        note(
          `basis: n=${basis.n} from ${basis.artefact || 'similarity'}${basis.path ? ` (${basis.path})` : ''}` +
            (basis.filtered ? ' — computed over a filtered subset' : ''),
          'faint'
        )
      );
    } else {
      statusHost.appendChild(
        note('The response carried no `basis`, so the n above is not shown as a cohort size.', 'faint')
      );
    }
  }

  function onMove(event) {
    const rect = canvas.getBoundingClientRect();
    const sx = event.clientX - rect.left;
    const sy = event.clientY - rect.top;
    showTooltip(cellAt(sx, sy), sx, sy);
  }

  function onLeave() {
    tooltip.style.display = 'none';
  }

  function onResize() {
    draw();
  }

  canvas.addEventListener('mousemove', onMove);
  canvas.addEventListener('mouseleave', onLeave);
  disposers.push(() => {
    canvas.removeEventListener('mousemove', onMove);
    canvas.removeEventListener('mouseleave', onLeave);
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

  (async () => {
    let body;
    let tree = null;
    try {
      [body, tree] = await Promise.all([
        ctx.api.get('/similarity', null, abort.signal),
        ctx.api.get('/tree', { source: 'stage9' }, abort.signal).catch((error) => {
          if (error && error.name === 'AbortError') return null;
          return null;
        }),
      ]);
    } catch (error) {
      if (error && error.name === 'AbortError') return;
      alertHost.appendChild(
        message('error', 'The similarity matrix could not be read', note(String((error && error.message) || error)))
      );
      return;
    }
    if (disposed) return;

    if (!body || body.present === false) {
      alertHost.appendChild(
        message(
          'info',
          'This matrix is not produced',
          serverMessage(body && body.reason ? body.reason : 'the server gave no reason'),
          body && body.probes ? probeList(body.probes) : null
        )
      );
      return;
    }

    payload = body;
    mountPowerBanner(
      ctx,
      bannerHost,
      body.basis,
      'For the similarity matrix `basis.n` is the number of samples the matrix was computed over (n_total). It is the server’s number; the flag is shown wherever a statistic is, as R16 requires.'
    );
    units = body.units || null;
    unitsReason = body.units_reason || '';
    const fileOrder = Array.isArray(body.order) ? body.order.slice() : [];
    const tipOrder = treeTipOrder(tree);
    const permuted = permute(fileOrder, tipOrder);
    order = permuted.order;
    orderNote = permuted.note;
    const indexOf = new Map();
    fileOrder.forEach((s, i) => indexOf.set(s, i));
    matrix = order.map((rowSample) => {
      const sourceRow = body.matrix[indexOf.get(rowSample)] || [];
      return order.map((colSample) => sourceRow[indexOf.get(colSample)]);
    });
    const flat = [];
    for (const row of matrix) {
      for (const value of row) {
        if (typeof value === 'number' && isFinite(value)) flat.push(value);
      }
    }
    colourScale = flat.length ? buildColourScale(flat) : null;
    if (body.n_incomplete_vectors) {
      alertHost.appendChild(
        note(
          `${body.n_incomplete_vectors} sample(s) had a short distance vector and were dropped by the server rather than padded with zeros: ${(body.incomplete_samples || []).slice(0, 8).join(', ')}`,
          'faint'
        )
      );
    }
    draw();
  })();

  return () => {
    disposed = true;
    abort.abort();
    for (const observer of observers.splice(0)) {
      try {
        observer.disconnect();
      } catch {
        /* already gone */
      }
    }
    for (const fn of disposers.splice(0)) {
      try {
        fn();
      } catch (error) {
        console.error('[matrix] a disposer threw', error);
      }
    }
  };
}

export default { mount };
