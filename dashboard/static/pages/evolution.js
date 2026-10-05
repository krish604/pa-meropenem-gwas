/* pages/evolution.js — V4, convergence and co-occurrence.
 *
 * Two server-paged tables (reusing DATA's `createPagedGrid`, so paging,
 * sorting, filtering and the per-response basis are identical to the rest of
 * the dashboard) plus a circular co-occurrence network drawn on canvas.
 *
 * Every statistic here mounts the one D3 banner, with `N` from `basis.n`.
 * Convergence and co-occurrence are descriptions of THIS run's isolates, not
 * findings, and the pages say so in the server's own words
 * (`limit_note` / `interpretation_limit`).
 */

import {
  clear,
  h,
  message,
  note,
  probeList,
  serverMessage,
} from '../app.js';
import { createPagedGrid, mountPowerBanner } from './tables.js';

const CONVERGENCE_CATEGORIES = [
  'widespread_background',
  'rare_isolated',
  'lineage_associated',
  'recurrent_convergent',
  'unknown',
];

export function mount(container, ctx) {
  const abort = new AbortController();
  let disposed = false;
  const disposers = [];
  const observers = [];
  let grids = [];

  const head = h(
    'div',
    { class: 'page-head' },
    h(
      'div',
      { class: 'page-head-text' },
      h('h2', null, 'Convergence & co-occurrence'),
      h(
        'p',
        { class: 'page-note' },
        'Stage 13 convergence classifications and stage 14 co-occurring feature pairs, ' +
          'both server-paged. A pattern in the data is not a demonstration of a shared ' +
          'selective cause, and co-occurrence is association only — never an interaction.'
      )
    )
  );
  container.appendChild(head);

  const alertHost = h('div');
  const convergenceHost = h('div', { class: 'panel' });
  const coocHost = h('div', { class: 'panel' });
  const networkHost = h('div', { class: 'panel' });
  container.append(alertHost, convergenceHost, coocHost, networkHost);

  function tokens() {
    try {
      return ctx.theme.tokens || {};
    } catch {
      return {};
    }
  }

  function buildConvergenceGrid() {
    convergenceHost.appendChild(h('h3', null, 'Convergence (stage 13)'));
    convergenceHost.appendChild(
      note(
        'The five categories are mutually exclusive and evaluated in a fixed order, which is ' +
          'why the order is not re-sorted for display. Filtering is server-side.'
      )
    );
    const categoryNote = h('div');
    const grid = createPagedGrid({
      ctx,
      path: '/convergence',
      columns: [
        { key: 'determinant', label: 'determinant', sortKey: 'determinant' },
        { key: 'independent_lineages', label: 'independent lineages', sortKey: 'independent_lineages' },
        { key: 'branch_count', label: 'branch count', sortKey: 'branch_count' },
        { key: 'distribution', label: 'distribution', sortKey: 'distribution' },
        { key: 'convergence_category', label: 'category', sortKey: 'convergence_category' },
      ],
      project: false,
      state: { offset: 0, limit: 100, sort: null, q: null, filters: {}, hidden: new Set() },
      search: { label: 'search (server-side substring)', placeholder: 'determinant, category…' },
      filters: [
        {
          param: 'convergence_category',
          label: 'category',
          kind: 'select',
          options: CONVERGENCE_CATEGORIES,
          hint: 'the five mutually-exclusive categories, in their evaluated order',
        },
      ],
      exportName: 'convergence',
      exportLabel: 'determinants',
      onLoad: (payload) => {
        if (payload && payload.limit_note) {
          clear(categoryNote);
          categoryNote.appendChild(note(payload.limit_note));
          if (Array.isArray(payload.category_order)) {
            categoryNote.appendChild(
              note(`evaluated order: ${payload.category_order.join(' → ')}`, 'faint')
            );
          }
        }
      },
      bannerNote:
        'Convergence is a pattern in this run’s data. `basis.n` is the number of determinant rows the statistic covers, as the server computed it.',
    });
    convergenceHost.append(grid.el, categoryNote);
    grids.push(grid);
    void grid.load();
  }

  function buildCooccurrenceGrid() {
    coocHost.appendChild(h('h3', null, 'Co-occurrence (stage 14)'));
    coocHost.appendChild(
      note(
        '`interpretation_limit` is rendered beside every pair. Co-occurrence is association ' +
          'only and is never presented as an interaction.'
      )
    );
    const grid = createPagedGrid({
      ctx,
      path: '/cooccurrence',
      columns: [
        { key: 'feature_a', label: 'feature A', sortKey: 'feature_a' },
        { key: 'feature_b', label: 'feature B', sortKey: 'feature_b' },
        { key: 'feature_type', label: 'type', sortKey: 'feature_type' },
        { key: 'n_a', label: 'n A', sortKey: 'n_a' },
        { key: 'n_b', label: 'n B', sortKey: 'n_b' },
        { key: 'n_both', label: 'n both', sortKey: 'n_both' },
        { key: 'statistic', label: 'statistic', sortKey: 'statistic' },
        { key: 'statistic_value', label: 'value', sortKey: 'statistic_value' },
        { key: 'adjusted_p_value', label: 'adj p', sortKey: 'adjusted_p_value' },
        {
          key: 'interpretation_limit',
          label: 'interpretation limit',
          sortKey: 'interpretation_limit',
          render: (row) => h('span', { class: 'faint' }, String(row.interpretation_limit || 'not reported')),
        },
      ],
      project: false,
      state: { offset: 0, limit: 100, sort: 'adjusted_p_value', q: null, filters: {}, hidden: new Set() },
      search: { label: 'search (server-side substring)', placeholder: 'feature, type…' },
      exportName: 'cooccurrence',
      exportLabel: 'pairs',
      bannerNote:
        '`basis.n` is the number of feature pairs the statistic covers, as the server computed it — not a cohort size.',
    });
    coocHost.appendChild(grid.el);
    grids.push(grid);
    void grid.load();
  }

  /* --------------------------------------------------------- network */

  let networkPayload = null;
  const networkCanvasWrap = h('div', {
    style: { position: 'relative', height: '560px', overflow: 'hidden' },
  });
  const networkCanvas = h('canvas', { style: { display: 'block', width: '100%', height: '100%' } });
  const networkTip = h('div', {
    style: {
      position: 'absolute',
      'pointer-events': 'none',
      display: 'none',
      'z-index': '5',
      background: 'var(--bg-raised)',
      border: '1px solid var(--border-strong)',
      'border-radius': 'var(--radius)',
      'box-shadow': 'var(--shadow)',
      padding: 'var(--space-1) var(--space-2)',
      'font-size': 'var(--fs-xs)',
      'font-family': 'var(--font-mono)',
    },
  });
  const networkBanner = h('div');
  const networkLegend = h('div');
  const networkNote = h('div');
  networkCanvasWrap.append(networkCanvas, networkTip);
  networkHost.append(
    h('h3', null, 'Co-occurrence network'),
    note(
      'A circular layout, not a force layout: every pair is association, and a force layout’s ' +
        'distances would imply a structure the test did not measure. Node size is the pair count; ' +
        'edge width is −log10(adjusted p), capped. Only the most significant pairs are drawn.'
    ),
    networkBanner,
    networkCanvasWrap,
    networkLegend,
    networkNote
  );

  const networkCtx = networkCanvas.getContext('2d');
  let nodePositions = new Map();
  let networkEdges = [];

  function buildNetwork(payload) {
    const t = tokens();
    const items = (payload.items || []).filter(
      (row) => row.feature_a && row.feature_b && typeof row.adjusted_p_value === 'number'
    );
    // cap to the most significant pairs so the picture stays readable
    const sorted = items.slice().sort((a, b) => a.adjusted_p_value - b.adjusted_p_value);
    const shown = sorted.slice(0, 120);
    networkEdges = shown;
    const nodes = new Map();
    for (const row of shown) {
      for (const feature of [row.feature_a, row.feature_b]) {
        if (!nodes.has(feature)) nodes.set(feature, { id: feature, count: 0 });
        nodes.get(feature).count += 1;
      }
    }
    const ids = [...nodes.keys()].sort();
    const rect = networkCanvasWrap.getBoundingClientRect
      ? networkCanvasWrap.getBoundingClientRect()
      : { width: 760, height: 560 };
    const width = Math.max(320, Math.floor(rect.width || 760));
    const height = Math.max(320, Math.floor(rect.height || 560));
    const dpr = globalThis.devicePixelRatio || 1;
    networkCanvas.width = Math.floor(width * dpr);
    networkCanvas.height = Math.floor(height * dpr);
    networkCanvas.style.width = `${width}px`;
    networkCanvas.style.height = `${height}px`;
    if (networkCtx.setTransform) networkCtx.setTransform(dpr, 0, 0, dpr, 0, 0);
    networkCtx.clearRect(0, 0, width, height);

    const cx = width / 2;
    const cy = height / 2;
    const radius = Math.max(60, Math.min(width, height) / 2 - 40);
    nodePositions = new Map();
    ids.forEach((id, i) => {
      const angle = (i / Math.max(1, ids.length)) * Math.PI * 2 - Math.PI / 2;
      nodePositions.set(id, {
        x: cx + radius * Math.cos(angle),
        y: cy + radius * Math.sin(angle),
        angle,
        count: nodes.get(id).count,
      });
    });

    const maxWeight = Math.max(1, ...shown.map((row) => -Math.log10(Math.max(row.adjusted_p_value, 1e-12))));
    networkCtx.lineCap = 'round';
    for (const row of shown) {
      const a = nodePositions.get(row.feature_a);
      const b = nodePositions.get(row.feature_b);
      if (!a || !b) continue;
      const weight = -Math.log10(Math.max(row.adjusted_p_value, 1e-12));
      const alpha = 0.15 + 0.6 * (weight / maxWeight);
      networkCtx.strokeStyle = t['tone-warn-fg'] || '#a80';
      networkCtx.globalAlpha = alpha;
      networkCtx.lineWidth = 0.5 + 3 * (weight / maxWeight);
      networkCtx.beginPath();
      networkCtx.moveTo(a.x, a.y);
      networkCtx.lineTo(b.x, b.y);
      networkCtx.stroke();
    }
    networkCtx.globalAlpha = 1;
    for (const [id, position] of nodePositions) {
      networkCtx.beginPath();
      networkCtx.fillStyle = t.accent || '#1f5fbf';
      const nodeRadius = 2.5 + Math.min(6, position.count);
      networkCtx.arc(position.x, position.y, nodeRadius, 0, Math.PI * 2);
      networkCtx.fill();
      if (ids.length <= 60) {
        networkCtx.fillStyle = t.fg || '#111';
        networkCtx.font = `9px ${t['font-mono'] || 'monospace'}`;
        networkCtx.textAlign = Math.cos(position.angle) >= 0 ? 'left' : 'right';
        networkCtx.textBaseline = 'middle';
        networkCtx.fillText(id, position.x + (Math.cos(position.angle) >= 0 ? 5 : -5), position.y);
      }
    }

    clear(networkLegend);
    networkLegend.appendChild(
      note(
        `drawn: ${shown.length} most-significant pair(s) of ${items.length} with a numeric adjusted p; ${ids.length} feature(s). Edges fade and thin as adjusted p grows.`,
        'faint'
      )
    );
    clear(networkNote);
    networkNote.appendChild(
      note(
        'Co-occurrence is association only and is never presented as an interaction. A pair here ' +
          'is not evidence that one feature acts on the other.',
        'faint'
      )
    );
    mountPowerBanner(
      ctx,
      networkBanner,
      payload.basis || (payload.meta && payload.meta.basis),
      '`basis.n` is the number of feature pairs the statistic covers, as the server computed it — not a cohort size.'
    );
  }

  function networkHit(mx, my) {
    for (const [id, position] of nodePositions) {
      const dx = position.x - mx;
      const dy = position.y - my;
      if (dx * dx + dy * dy <= 100) return id;
    }
    return null;
  }

  function onNetworkMove(event) {
    const rect = networkCanvas.getBoundingClientRect();
    const mx = event.clientX - rect.left;
    const my = event.clientY - rect.top;
    const id = networkHit(mx, my);
    if (!id) {
      networkTip.style.display = 'none';
      return;
    }
    const related = networkEdges.filter(
      (row) => row.feature_a === id || row.feature_b === id
    );
    clear(networkTip);
    networkTip.appendChild(h('div', null, id));
    networkTip.appendChild(
      document.createTextNode(`${related.length} drawn pair(s) in this view`)
    );
    networkTip.style.display = 'block';
    networkTip.style.left = `${Math.round(mx + 12)}px`;
    networkTip.style.top = `${Math.round(my + 12)}px`;
  }

  function onNetworkLeave() {
    networkTip.style.display = 'none';
  }

  networkCanvas.addEventListener('mousemove', onNetworkMove);
  networkCanvas.addEventListener('mouseleave', onNetworkLeave);
  disposers.push(() => {
    networkCanvas.removeEventListener('mousemove', onNetworkMove);
    networkCanvas.removeEventListener('mouseleave', onNetworkLeave);
  });

  if (typeof ResizeObserver === 'function') {
    const observer = new ResizeObserver(() => {
      if (networkPayload) buildNetwork(networkPayload);
    });
    observer.observe(networkCanvasWrap);
    observers.push(observer);
  }
  if (typeof MutationObserver === 'function' && typeof document !== 'undefined') {
    const observer = new MutationObserver(() => {
      if (networkPayload) buildNetwork(networkPayload);
    });
    observer.observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
    observers.push(observer);
  }

  buildConvergenceGrid();
  buildCooccurrenceGrid();

  (async () => {
    try {
      const body = await ctx.api.get('/cooccurrence', { limit: 1000, sort: 'adjusted_p_value' }, abort.signal);
      if (disposed) return;
      if (!body || body.present === false) {
        clear(networkCanvasWrap);
        networkBanner.appendChild(
          message(
            'info',
            'The co-occurrence table is not produced',
            serverMessage(body && body.reason ? body.reason : 'the server gave no reason'),
            body && body.probes ? probeList(body.probes) : null
          )
        );
        return;
      }
      networkPayload = body;
      if (body.limit_note) {
        networkNote.appendChild(note(body.limit_note));
      }
      buildNetwork(body);
    } catch (error) {
      if (error && error.name === 'AbortError') return;
      clear(networkCanvasWrap);
      alertHost.appendChild(
        message('error', 'The co-occurrence network could not be read', note(String((error && error.message) || error)))
      );
    }
  })();

  return () => {
    disposed = true;
    abort.abort();
    for (const grid of grids) {
      try {
        grid.destroy();
      } catch (error) {
        console.error('[evolution] a grid disposer threw', error);
      }
    }
    grids = [];
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
        console.error('[evolution] a disposer threw', error);
      }
    }
  };
}

export default { mount };
