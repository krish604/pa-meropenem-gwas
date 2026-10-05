/* pages/gwas.js — V3, the association results.
 *
 * Four views over `/api/gwas` and `/api/gwas/top`:
 *   - a Manhattan-style plot,
 *   - a QQ plot,
 *   - an effect-versus-significance plot,
 *   - a top-hits table plus the full server-paged association table.
 *
 * The D3 rule (UI-D3) is enforced by mounting `mountPowerBanner` on every
 * statistic, and `N` is `basis.n` from the response — never the page's own
 * count. The GWAS table carries no genomic coordinates, so the x-axis of the
 * "Manhattan" is the feature's index in the server's own order; calling it a
 * genomic Manhattan would be a claim the data does not support, and the page
 * says so.
 *
 * `model` is shown verbatim. Every reference-engine row is stamped
 * `reference_fisher:NO_KINSHIP_CORRECTION`; a row whose model does not name its
 * correction is flagged, never hidden, and never presented as a finding.
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
import { createPagedGrid, mountPowerBanner } from './tables.js';

function d3lib() {
  return typeof globalThis !== 'undefined' ? globalThis.d3 : undefined;
}

function hashHue(text) {
  let hash = 0;
  const value = String(text);
  for (let i = 0; i < value.length; i += 1) hash = (hash * 31 + value.charCodeAt(i)) % 360;
  return hash;
}

function categorical(text, fallback) {
  if (text === null || text === undefined || text === '') return fallback || '#888';
  return `hsl(${hashHue(text)}, 55%, 45%)`;
}

function numberOrNull(value) {
  if (value === null || value === undefined || value === '') return null;
  const n = Number(value);
  return isFinite(n) ? n : null;
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
      h('h2', null, 'GWAS'),
      h(
        'p',
        { class: 'page-note' },
        'Association results (stage 12) from the reference engine. Every plot carries the ' +
          'R16 flag with n taken from the response’s basis, never from what is on the page. ' +
          'No row is a finding: the flag is shown wherever a statistic is.'
      )
    )
  );
  container.appendChild(head);

  const alertHost = h('div');
  const bannerHost = h('div');
  const plotGrid = h('div', { class: 'panel-grid' });
  const topHost = h('div', { class: 'panel' });
  const gridHost = h('div', { class: 'panel' });
  container.append(alertHost, bannerHost, plotGrid, topHost, gridHost);

  const plots = [];
  let gwasPayload = null;
  let grid = null;

  function tokens() {
    try {
      return ctx.theme.tokens || {};
    } catch {
      return {};
    }
  }

  function makePlot(title, subtitle) {
    const wrap = h('div', { class: 'panel' });
    wrap.appendChild(h('h3', null, title));
    if (subtitle) wrap.appendChild(note(subtitle));
    const canvasWrap = h('div', {
      style: { position: 'relative', height: '340px', overflow: 'hidden' },
    });
    const canvas = h('canvas', { style: { display: 'block', width: '100%', height: '100%' } });
    const tip = h('div', {
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
    canvasWrap.append(canvas, tip);
    const noteHost = h('div');
    wrap.append(canvasWrap, noteHost);
    plotGrid.appendChild(wrap);
    const plot = { wrap, canvasWrap, canvas, tip, noteHost, describe: null, hit: null };
    // One hover listener per plot, attached once. `describe` is swapped by
    // `renderPlots`; re-attaching a listener on every redraw would leak.
    const onMove = (event) => {
      if (!plot.hit || !plot.describe) return;
      const rect = canvas.getBoundingClientRect();
      const mx = event.clientX - rect.left;
      const my = event.clientY - rect.top;
      let best = null;
      let bestD = 12;
      for (const point of plot.hit.points) {
        const dx = plot.hit.sx(point.x) - mx;
        const dy = plot.hit.sy(point.y) - my;
        const distance = Math.sqrt(dx * dx + dy * dy);
        if (distance < bestD) {
          bestD = distance;
          best = point;
        }
      }
      if (!best) {
        tip.style.display = 'none';
        return;
      }
      clear(tip);
      tip.appendChild(document.createTextNode(plot.describe(best)));
      tip.style.display = 'block';
      tip.style.left = `${Math.round(mx + 12)}px`;
      tip.style.top = `${Math.round(my + 12)}px`;
    };
    const onLeave = () => {
      tip.style.display = 'none';
    };
    canvas.addEventListener('mousemove', onMove);
    canvas.addEventListener('mouseleave', onLeave);
    disposers.push(() => {
      canvas.removeEventListener('mousemove', onMove);
      canvas.removeEventListener('mouseleave', onLeave);
    });
    plots.push(plot);
    return plot;
  }

  /**
   * Draw a scatter with hand-rolled canvas axes. d3-scale supplies the tick
   * values, d3-format the labels. No SVG, no d3-axis: rendering is canvas.
   */
  function drawScatter(plot, points, opts) {
    const canvas = plot.canvas;
    const ctx2d = canvas.getContext('2d');
    const t = tokens();
    const rect = plot.canvasWrap.getBoundingClientRect
      ? plot.canvasWrap.getBoundingClientRect()
      : { width: 460, height: 340 };
    const width = Math.max(280, Math.floor(rect.width || 460));
    const height = Math.max(240, Math.floor(rect.height || 340));
    const dpr = globalThis.devicePixelRatio || 1;
    canvas.width = Math.floor(width * dpr);
    canvas.height = Math.floor(height * dpr);
    canvas.style.width = `${width}px`;
    canvas.style.height = `${height}px`;
    if (ctx2d.setTransform) ctx2d.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx2d.clearRect(0, 0, width, height);
    const margin = { left: 62, right: 14, top: 14, bottom: 46 };
    const plotW = width - margin.left - margin.right;
    const plotH = height - margin.top - margin.bottom;

    const xs = points.map((p) => p.x);
    const ys = points.map((p) => p.y);
    const xDomain = opts.xDomain || [Math.min(...xs, 0), Math.max(...xs, 1)];
    const yDomain = opts.yDomain || [0, Math.max(...ys, 1)];
    const sx = (v) => margin.left + ((v - xDomain[0]) / (xDomain[1] - xDomain[0] || 1)) * plotW;
    const sy = (v) => margin.top + plotH - ((v - yDomain[0]) / (yDomain[1] - yDomain[0] || 1)) * plotH;

    ctx2d.fillStyle = t.fg || '#111';
    ctx2d.strokeStyle = t.border || '#ccc';
    ctx2d.lineWidth = 1;
    ctx2d.font = `10px ${t['font-mono'] || 'monospace'}`;

    // axes frame
    ctx2d.strokeRect(margin.left + 0.5, margin.top + 0.5, plotW, plotH);

    // ticks
    const xTicks = d3 && d3.scaleLinear
      ? d3.scaleLinear().domain(xDomain).ticks(Math.min(8, Math.max(2, Math.floor(plotW / 70))))
      : [xDomain[0], (xDomain[0] + xDomain[1]) / 2, xDomain[1]];
    const yTicks = d3 && d3.scaleLinear
      ? d3.scaleLinear().domain(yDomain).ticks(5)
      : [yDomain[0], (yDomain[0] + yDomain[1]) / 2, yDomain[1]];
    ctx2d.textAlign = 'center';
    ctx2d.textBaseline = 'top';
    for (const tick of xTicks) {
      if (tick < xDomain[0] || tick > xDomain[1]) continue;
      const px = sx(tick);
      ctx2d.beginPath();
      ctx2d.moveTo(px, margin.top + plotH);
      ctx2d.lineTo(px, margin.top + plotH + 4);
      ctx2d.stroke();
      ctx2d.fillText(d3 && d3.format ? d3.format(opts.xFormat || '~s')(tick) : String(tick), px, margin.top + plotH + 6);
    }
    ctx2d.textAlign = 'right';
    ctx2d.textBaseline = 'middle';
    for (const tick of yTicks) {
      if (tick < yDomain[0] || tick > yDomain[1]) continue;
      const py = sy(tick);
      ctx2d.beginPath();
      ctx2d.moveTo(margin.left - 4, py);
      ctx2d.lineTo(margin.left, py);
      ctx2d.stroke();
      ctx2d.fillText(d3 && d3.format ? d3.format(opts.yFormat || '.1f')(tick) : String(tick), margin.left - 6, py);
    }

    // identity diagonal (QQ)
    if (opts.diagonal) {
      ctx2d.save();
      ctx2d.strokeStyle = t['fg-faint'] || '#999';
      ctx2d.setLineDash([3, 3]);
      ctx2d.beginPath();
      ctx2d.moveTo(sx(xDomain[0]), sy(yDomain[0]));
      ctx2d.lineTo(sx(xDomain[1]), sy(yDomain[1]));
      ctx2d.stroke();
      ctx2d.restore();
    }

    // threshold line
    if (opts.threshold !== undefined && opts.threshold !== null && opts.threshold >= yDomain[0] && opts.threshold <= yDomain[1]) {
      ctx2d.save();
      ctx2d.strokeStyle = t['tone-warn-fg'] || '#a80';
      ctx2d.setLineDash([4, 3]);
      ctx2d.beginPath();
      ctx2d.moveTo(margin.left, sy(opts.threshold));
      ctx2d.lineTo(margin.left + plotW, sy(opts.threshold));
      ctx2d.stroke();
      ctx2d.restore();
      ctx2d.fillStyle = t['tone-warn-fg'] || '#a80';
      ctx2d.textAlign = 'left';
      ctx2d.textBaseline = 'bottom';
      ctx2d.fillText(opts.thresholdLabel || '', margin.left + 4, sy(opts.threshold) - 2);
    }

    // points
    for (const point of points) {
      ctx2d.beginPath();
      ctx2d.fillStyle = point.colour || (t.accent || '#1f5fbf');
      ctx2d.globalAlpha = point.alpha === undefined ? 0.75 : point.alpha;
      const r = point.r || 2.2;
      ctx2d.arc(sx(point.x), sy(point.y), r, 0, Math.PI * 2);
      ctx2d.fill();
    }
    ctx2d.globalAlpha = 1;

    // axis titles
    ctx2d.fillStyle = t['fg-muted'] || '#555';
    ctx2d.textAlign = 'center';
    ctx2d.textBaseline = 'top';
    ctx2d.fillText(opts.xLabel || '', margin.left + plotW / 2, height - 14);
    ctx2d.save();
    ctx2d.translate(12, margin.top + plotH / 2);
    ctx2d.rotate(-Math.PI / 2);
    ctx2d.textBaseline = 'middle';
    ctx2d.fillText(opts.yLabel || '', 0, 0);
    ctx2d.restore();

    plot.hit = { points, sx, sy, margin, plotW, plotH, xDomain, yDomain };
  }

  function renderPlots() {
    const t = tokens();
    const items = (gwasPayload.items || []).map((row, index) => ({ row, index }));
    const withP = items
      .map((entry) => ({ ...entry, p: numberOrNull(entry.row.p_value) }))
      .filter((entry) => entry.p !== null && entry.p > 0);
    const missingP = items.length - withP.length;

    // Manhattan
    const manhattan = plots[0];
    const mPoints = withP.map((entry) => ({
      x: entry.index,
      y: -Math.log10(entry.p),
      colour: categorical(entry.row.feature_type, t['fg-faint']),
      row: entry.row,
      index: entry.index,
    }));
    drawScatter(manhattan, mPoints, {
      xLabel: 'feature index (the server’s own order — this table carries no genomic coordinate)',
      yLabel: '−log10(p_value)',
      yDomain: [0, Math.max(1, ...mPoints.map((p) => p.y)) * 1.05],
      threshold: -Math.log10(0.05),
      thresholdLabel: 'nominal p = 0.05',
      xFormat: 'd',
    });
    manhattan.describe = (p) => `${p.row.feature}  p=${p.row.p_value}  −log10=${p.y.toFixed(2)}`;
    clear(manhattan.noteHost);
    if (missingP) {
      manhattan.noteHost.appendChild(
        note(
          `${missingP} row(s) have no numeric p_value, or a p_value of exactly 0 (which a log axis cannot place). They are counted here, never drawn at zero.`,
          'faint'
        )
      );
    }

    // QQ
    const qq = plots[1];
    const rawP = Array.isArray(gwasPayload.p_values) ? gwasPayload.p_values.slice().sort((a, b) => a - b) : [];
    const pValues = rawP.filter((p) => typeof p === 'number' && p > 0);
    const nTests = typeof gwasPayload.n_tests === 'number' ? gwasPayload.n_tests : pValues.length;
    const qPoints = pValues.map((p, i) => {
      const expected = (i + 0.5) / nTests;
      return { x: -Math.log10(expected), y: -Math.log10(p), p, i };
    });
    const maxAxis = Math.max(1, ...qPoints.map((p) => Math.max(p.x, p.y))) * 1.05;
    drawScatter(qq, qPoints, {
      xLabel: 'expected −log10(p) under the null (uniform, from n_tests)',
      yLabel: 'observed −log10(p_value)',
      xDomain: [0, maxAxis],
      yDomain: [0, maxAxis],
      xFormat: '.1f',
      yFormat: '.1f',
      diagonal: true,
    });
    qq.describe = (p) => `observed p=${p.p}  expected rank ${p.i + 1} of ${nTests}`;
    clear(qq.noteHost);
    qq.noteHost.appendChild(
      note(
        `Expected quantiles use the server’s n_tests = ${nTests} (the row count over the whole artefact), not the ${pValues.length} values on the page.`,
        'faint'
      )
    );
    if (gwasPayload.p_values_truncated) {
      qq.noteHost.appendChild(
        message(
          'warn',
          'The p-value array is truncated',
          note(
            'The server caps the QQ array, and it caps it after sorting, so it holds the most significant values. ' +
              'The QQ over them is therefore not the whole test set; the n_tests denominator above is the whole set.'
          )
        )
      );
    }

    // Effect vs significance
    const effect = plots[2];
    const ePoints = [];
    let missingEffect = 0;
    for (const entry of withP) {
      const effectSize = numberOrNull(entry.row.effect_size);
      if (effectSize === null) {
        missingEffect += 1;
        continue;
      }
      const corrected = entry.row.model_states_kinship_correction === true;
      ePoints.push({
        x: effectSize,
        y: -Math.log10(entry.p),
        colour: corrected ? (t['tone-ok-fg'] || '#080') : (t['tone-bad-fg'] || '#a00'),
        row: entry.row,
      });
    }
    const xExtent = ePoints.length
      ? [Math.min(...ePoints.map((p) => p.x)), Math.max(...ePoints.map((p) => p.x))]
      : [-1, 1];
    drawScatter(effect, ePoints, {
      xLabel: 'effect_size',
      yLabel: '−log10(p_value)',
      xDomain: xExtent[0] === xExtent[1] ? [xExtent[0] - 1, xExtent[1] + 1] : xExtent,
      yDomain: [0, Math.max(1, ...ePoints.map((p) => p.y)) * 1.05],
      xFormat: '.2g',
    });
    effect.describe = (p) => `${p.row.feature}  effect=${p.x}  −log10(p)=${p.y.toFixed(2)}`;
    clear(effect.noteHost);
    effect.noteHost.appendChild(
      note(
        `green = model names a kinship correction; red = it does not (e.g. reference_fisher). ${missingEffect} row(s) have no numeric effect_size and are not drawn at zero.`,
        'faint'
      )
    );
  }

  function renderTop(topPayload) {
    clear(topHost);
    topHost.appendChild(h('h3', null, 'Top hits by adjusted p-value'));
    const banner = h('div');
    topHost.appendChild(banner);
    mountPowerBanner(
      ctx,
      banner,
      topPayload.basis,
      'For the top-hits endpoint `basis.n` is the number of hits RETURNED by this request, not the number of tests and not the cohort. It is shown because it is the server’s number; the plot banner above carries the test count.'
    );
    const items = topPayload.items || [];
    if (!items.length) {
      topHost.appendChild(absentSpan(topPayload.reason || 'the endpoint returned no hit', ABSENT.NOT_PRODUCED));
      return;
    }
    const table = h('table', { class: 'data-table' });
    const thead = h('thead');
    const hr = h('tr');
    for (const label of ['feature', 'adjusted p', 'p', 'dominant lineage share', 'lineage-linked']) {
      hr.appendChild(h('th', null, label));
    }
    thead.appendChild(hr);
    const tbody = h('tbody');
    for (const row of items) {
      const tr = h('tr');
      tr.appendChild(h('td', { class: 'mono' }, row.feature || ''));
      tr.appendChild(h('td', { class: 'mono' }, row.adjusted_p_value === null || row.adjusted_p_value === undefined ? absentSpan(ABSENT.NOT_REPORTED) : String(row.adjusted_p_value)));
      tr.appendChild(h('td', { class: 'mono' }, row.p_value === null || row.p_value === undefined ? absentSpan(ABSENT.NOT_REPORTED) : String(row.p_value)));
      tr.appendChild(h('td', { class: 'mono' }, row.dominant_lineage_share === null || row.dominant_lineage_share === undefined ? absentSpan(ABSENT.NOT_ASSESSED) : String(row.dominant_lineage_share)));
      const linkedCell = h('td');
      if (row.lineage_linked) {
        linkedCell.appendChild(
          h('span', { class: 'chip', dataset: { flag: 'true' } }, `lineage-linked — relabelled, never a determinant (threshold ${topPayload.lineage_confound_threshold})`)
        );
      } else {
        linkedCell.appendChild(document.createTextNode('not lineage-linked'));
      }
      tr.appendChild(linkedCell);
      tbody.appendChild(tr);
    }
    table.append(thead, tbody);
    topHost.appendChild(h('div', { class: 'table-wrap' }, table));
    topHost.appendChild(
      note(
        'A feature carried only by one lineage cannot be separated from that lineage by any ' +
          'association test. Such a feature is RELABELLED, never deleted, and never presented ' +
          'as a resistance determinant.'
      )
    );
  }

  function renderGrid() {
    const columns = [
      { key: 'feature', label: 'feature', sortKey: 'feature' },
      { key: 'feature_type', label: 'type', sortKey: 'feature_type' },
      { key: 'effect', label: 'effect', sortKey: 'effect' },
      { key: 'p_value', label: 'p', sortKey: 'p_value' },
      { key: 'adjusted_p_value', label: 'adj p', sortKey: 'adjusted_p_value' },
      { key: 'effect_size', label: 'effect size', sortKey: 'effect_size' },
      { key: 'frequency', label: 'frequency', sortKey: 'frequency' },
      { key: 'lineage_distribution', label: 'lineage distribution', sortKey: 'lineage_distribution' },
      {
        key: 'model',
        label: 'model (verbatim)',
        sortKey: 'model',
        render: (row) => {
          const span = h('span', { class: 'mono' }, String(row.model || ''));
          if (row.model_states_kinship_correction !== true) {
            span.appendChild(h('span', { class: 'chip', dataset: { flag: 'true' } }, ' no kinship correction named'));
          }
          return span;
        },
      },
    ];
    gridHost.appendChild(
      h('h3', null, 'Every association, server-paged')
    );
    gridHost.appendChild(
      note(
        'Sorting, searching and paging are the server’s. The count beside the page is the ' +
          'server’s count over the whole table, not the rows on the page (UI-D7).'
      )
    );
    grid = createPagedGrid({
      ctx,
      path: '/gwas',
      columns,
      project: false,
      state: { offset: 0, limit: 100, sort: 'p_value', q: null, filters: {}, hidden: new Set() },
      search: { label: 'search (server-side substring)', placeholder: 'feature, type, model…' },
      exportName: 'gwas',
      exportLabel: 'associations',
      bannerNote:
        'For the GWAS table `basis.n` is the number of TESTS in the artefact (the server’s n_tests), not the number of isolates. It is shown because it is what the server reported; it is not a cohort size.',
    });
    gridHost.appendChild(grid.el);
    void grid.load();
  }

  makePlot(
    'Manhattan-style',
    '−log10(p) by feature index. A true genomic Manhattan needs chromosome and position, which the stage-12 table does not carry, so the axis is the feature’s index and is labelled as such.'
  );
  makePlot('QQ plot', 'Observed versus expected p-value quantiles. The expected line is y = x.');
  makePlot(
    'Effect versus significance',
    'effect_size against −log10(p). Colour marks whether the model names a kinship correction.'
  );

  if (typeof ResizeObserver === 'function') {
    const observer = new ResizeObserver(() => {
      if (gwasPayload) renderPlots();
    });
    observer.observe(plotGrid);
    observers.push(observer);
  }
  if (typeof MutationObserver === 'function' && typeof document !== 'undefined') {
    const observer = new MutationObserver(() => {
      if (gwasPayload) renderPlots();
    });
    observer.observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
    observers.push(observer);
  }

  (async () => {
    let body;
    let top;
    try {
      [body, top] = await Promise.all([
        ctx.api.get('/gwas', { limit: 1000 }, abort.signal),
        ctx.api.get('/gwas/top', { limit: 25 }, abort.signal).catch((error) => (error && error.name === 'AbortError' ? null : null)),
      ]);
    } catch (error) {
      if (error && error.name === 'AbortError') return;
      alertHost.appendChild(
        message('error', 'The GWAS results could not be read', note(String((error && error.message) || error)))
      );
      return;
    }
    if (disposed) return;

    if (!body || body.present === false) {
      alertHost.appendChild(
        message(
          'info',
          'The GWAS table is not produced',
          serverMessage(body && body.reason ? body.reason : 'the server gave no reason'),
          body && body.probes ? probeList(body.probes) : null
        )
      );
      return;
    }

    gwasPayload = body;
    mountPowerBanner(
      ctx,
      bannerHost,
      body.meta && body.meta.basis ? body.meta.basis : body.basis,
      'For the GWAS table `basis.n` is the number of TESTS (n_tests), which is what the server computed the flag on. The cohort size is a separate number and is not substituted here (UI-D3).'
    );
    renderPlots();
    if (top && top.present !== false) renderTop(top);
    else if (top) {
      clear(topHost);
      topHost.appendChild(message('info', 'Top hits are not produced', serverMessage(top.reason || 'the server gave no reason')));
    }
    renderGrid();
  })();

  return () => {
    disposed = true;
    abort.abort();
    if (grid) grid.destroy();
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
        console.error('[gwas] a disposer threw', error);
      }
    }
  };
}

export default { mount };
