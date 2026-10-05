/* pages/pipeline.js — the 16-stage DAG, and the live run monitor.
 *
 * Two exports, both on the §7 contract:
 *
 *   export function mount(container, ctx)          -> the DAG page
 *   export function mountMonitor(container, ctx)   -> the run monitor
 *
 * Three rules this file holds itself to:
 *
 *  1. **The stage list is the server's.** `STAGE_ORDER` arrives as
 *     `/api/stages`'s `stage_order` and as each row's `order`. Nothing here
 *     carries a stage name, a stage number or a filename.
 *  2. **Every badge goes through `ctx.badges`.** No status sentence is
 *     assembled in this file (UI-D2). Where the server already classified a
 *     stage, that classification is used verbatim, reason and all.
 *  3. **The D3 banner is mounted wherever a statistic renders**, with `N`
 *     taken from `basis.n` — the cohort the statistic was actually computed
 *     on, never the page's total (§7.1).
 */

import {
  ABSENT,
  absentSpan,
  basisOf,
  clear,
  dataTable,
  h,
  kv,
  message,
  note,
  probeList,
  serverMessage,
  svg,
} from '../app.js';

const DAG_NODE_W = 168;
const DAG_NODE_H = 42;
const DAG_GAP_X = 74;
const DAG_GAP_Y = 30;

/* ------------------------------------------------------------------ DAG */

/**
 * Observed transitions, derived from the event log's own timestamps.
 *
 * The API exposes **no `prerequisites` field**: `/api/stages` and `/api/run`
 * carry order, state, reason, presence and row counts, and nothing that
 * declares which stage needs which. So the edges drawn here are *observed*
 * transitions — a stage whose last event is timestamped after another
 * stage's — and the legend says which is which. Drawing a declared edge would
 * mean inventing a list the server never sent, and drawing nothing would hide
 * a real DAG.
 */
export function deriveEdges(stages, snapshot) {
  const state = (snapshot && snapshot.stage_event_state) || {};
  const timed = [];
  for (const stage of stages) {
    const entry = state[stage.name];
    const last = entry && entry.last_event ? entry.last_event : null;
    if (last && last.t) timed.push({ name: stage.name, t: String(last.t) });
  }
  timed.sort((a, b) => (a.t < b.t ? -1 : a.t > b.t ? 1 : 0));
  const observed = [];
  for (let i = 1; i < timed.length; i += 1) {
    observed.push({ from: timed[i - 1].name, to: timed[i].name, evidence: 'observed' });
  }
  const hasLog = Boolean(snapshot && (snapshot.cursor !== undefined || snapshot.n_events !== undefined));
  const edges = observed.slice();
  if (observed.length === 0 && stages.length > 1) {
    // No transition evidence. The order edges below are the server's declared
    // order, drawn dashed, and the legend says they are order and not need.
    for (let i = 1; i < stages.length; i += 1) {
      edges.push({ from: stages[i - 1].name, to: stages[i].name, evidence: 'order' });
    }
  }
  return { edges, evidence: observed.length > 0 ? 'observed' : hasLog ? 'order' : 'unavailable' };
}

function layout(stages) {
  const byName = new Map();
  const positions = new Map();
  stages.forEach((stage, index) => {
    byName.set(stage.name, stage);
    const column = index;
    const row = column % 2;
    positions.set(stage.name, {
      x: column * (DAG_NODE_W + DAG_GAP_X),
      y: row * (DAG_NODE_H + DAG_GAP_Y),
      stage,
    });
  });
  return { byName, positions };
}

function drawDag(svgEl, stages, edges, selected, onSelect) {
  const { positions } = layout(stages);
  const columns = stages.length;
  const rows = stages.length > 1 ? 2 : 1;
  svgEl.setAttribute(
    'viewBox',
    `0 0 ${columns * (DAG_NODE_W + DAG_GAP_X)} ${rows * (DAG_NODE_H + DAG_GAP_Y)}`
  );
  svgEl.setAttribute('width', String(columns * (DAG_NODE_W + DAG_GAP_X)));
  svgEl.setAttribute('height', String(rows * (DAG_NODE_H + DAG_GAP_Y)));
  svgEl.setAttribute('role', 'group');
  svgEl.setAttribute('aria-label', 'The pipeline DAG, in the order the server declared');

  const edgeLayer = svg('g', { class: 'dag-edges' });
  const nodeLayer = svg('g', { class: 'dag-nodes' });

  for (const edge of edges) {
    const from = positions.get(edge.from);
    const to = positions.get(edge.to);
    if (!from || !to) continue;
    const x1 = from.x + DAG_NODE_W;
    const y1 = from.y + DAG_NODE_H / 2;
    const x2 = to.x;
    const y2 = to.y + DAG_NODE_H / 2;
    const mid = (x1 + x2) / 2;
    edgeLayer.appendChild(
      svg('path', {
        class: 'dag-edge',
        'data-evidence': edge.evidence,
        d: `M${x1},${y1} C${mid},${y1} ${mid},${y2} ${x2},${y2}`,
      })
    );
  }

  for (const stage of stages) {
    const box = positions.get(stage.name);
    if (!box) continue;
    const group = svg('g', {
      class: 'dag-node',
      transform: `translate(${box.x},${box.y})`,
      'data-tone': badgeTone(stage),
      'data-stage': stage.name,
      tabindex: '0',
      role: 'button',
      'aria-label': `${stage.name}, ${badgeLabel(stage)}`,
    });
    if (stage.name === selected) group.setAttribute('aria-selected', 'true');
    group.appendChild(svg('rect', { width: DAG_NODE_W, height: DAG_NODE_H, rx: 4 }));
    // The glyph and the word are both on the node: a status carried by colour
    // alone disappears for a reader who cannot see the hue.
    group.appendChild(svg('text', { class: 'dag-glyph', x: 8, y: 16 }, glyphOf(stage)));
    group.appendChild(svg('text', { x: 24, y: 16 }, truncate(stage.name, 20)));
    group.appendChild(
      svg('text', { class: 'dag-spec', x: 8, y: 31 }, truncate(`#${stage.order} ${badgeLabel(stage)}`, 26))
    );
    const activate = () => onSelect(stage.name);
    group.addEventListener('click', activate);
    group.addEventListener('keydown', (event) => {
      if (event.key === 'Enter' || event.key === ' ') {
        event.preventDefault();
        activate();
      }
    });
    nodeLayer.appendChild(group);
  }

  clear(svgEl);
  svgEl.appendChild(edgeLayer);
  svgEl.appendChild(nodeLayer);
}

function truncate(text_, max) {
  const s = String(text_ === null || text_ === undefined ? '' : text_);
  return s.length > max ? `${s.slice(0, max - 1)}…` : s;
}

function badgeOf(stage) {
  return stage.badge && stage.badge.key ? stage.badge : stage.state;
}
function badgeTone(stage) {
  const badge = badgeOf(stage);
  const tone = typeof badge === 'string' ? badge : badge.tone;
  return ['ok', 'busy', 'bad', 'warn', 'unknown'].includes(tone) ? tone : 'unknown';
}
function badgeLabel(stage) {
  const badge = badgeOf(stage);
  if (typeof badge === 'string') return badge.replace(/_/g, ' ');
  return badge.label || badge.key || ABSENT.NOT_REPORTED;
}
function glyphOf(stage) {
  const badge = badgeOf(stage);
  return typeof badge === 'string' ? '?' : badge.glyph || '?';
}

/* ------------------------------------------------------- the DAG page */

export function mount(container, ctx) {
  const abort = new AbortController();
  const disposers = [];
  let disposed = false;
  let stages = [];
  let selected = null;
  let stageDetail = null;
  let snapshot = null;
  let run = null;

  container.appendChild(
    h(
      'div',
      { class: 'page-head' },
      h(
        'div',
        { class: 'page-head-text' },
        h('h2', null, 'Pipeline'),
        h(
          'p',
          { class: 'page-note' },
          'The ' +
            'stage list, the order and every state on it come from the server. ' +
            'Click a stage for its contract header, the header its file actually ' +
            'carries, its events, its folded-step tables and its outputs.'
        )
      )
    )
  );

  const bannerHost = h('div');
  const countsHost = h('div', { class: 'panel' });
  const dagHost = h('div', { class: 'panel' });
  const listHost = h('div', { class: 'panel' });
  const detailHost = h('div');
  container.append(bannerHost, countsHost, dagHost, listHost, detailHost);

  async function loadAll() {
    const [stagesResponse, runResponse] = await Promise.all([
      ctx.api.get('/stages', null, abort.signal),
      ctx.api.get('/run', null, abort.signal),
    ]);
    if (disposed) return;
    stages = Array.isArray(stagesResponse.items) ? stagesResponse.items : [];
    run = runResponse;
    renderCounts();
    renderDag();
    renderList();
  }

  function renderCounts() {
    clear(countsHost);
    countsHost.appendChild(
      h(
        'div',
        { class: 'panel-head' },
        h('h3', null, 'Stages'),
        h('span', { class: 'faint' }, `${stages.length} in the order the server declared`)
      )
    );
    if (stages.length === 0) {
      countsHost.appendChild(absentSpan('the server returned no stage rows', ABSENT.NOT_PRODUCED));
      return;
    }
    // The tally is a statistic, so it carries a basis and the banner fires.
    const row = h('div', { class: 'stat-row' });
    const basis = { n: stages.length, artefact: 'papipeline.run.STAGE_ORDER', min_samples: ctx.banner.minSamples };
    for (const stage of stages) {
      const stat = h('div', { class: 'stat', dataset: { tone: badgeTone(stage) } });
      stat.appendChild(h('span', { class: 'stat-value' }, stage.order === undefined ? '—' : String(stage.order)));
      stat.appendChild(h('span', { class: 'stat-label' }, stage.name));
      stat.appendChild(h('div', null, ctx.badges.render(stage)));
      row.appendChild(stat);
    }
    countsHost.appendChild(row);
    clear(bannerHost);
    ctx.banner.mount(bannerHost, null, basis);
    if (ctx.banner.mode) {
      countsHost.appendChild(
        note(
          `This run is ${ctx.banner.mode}. ${basis.n} stages is a count of stages, not of isolates; ` +
            'the isolate count and the R16 minimum are in the header.'
        )
      );
    }
  }

  function renderDag() {
    clear(dagHost);
    dagHost.appendChild(
      h(
        'div',
        { class: 'panel-head' },
        h('h3', null, 'DAG'),
        h('span', { class: 'faint' }, 'edges are observed transitions, drawn from event-log timestamps')
      )
    );
    if (stages.length === 0) {
      dagHost.appendChild(absentSpan('no stage rows, so there is no graph to draw', ABSENT.NOT_PRODUCED));
      return;
    }
    const { edges, evidence } = deriveEdges(stages, snapshot);
    const canvas = svg('svg', { class: 'dag' });
    drawDag(canvas, stages, edges, selected, (name) => select(name, true));
    dagHost.appendChild(h('div', { class: 'dag-canvas-wrap' }, canvas));

    const legend = h('div', { class: 'dag-legend' });
    const observed = h('span');
    observed.appendChild(h('span', { class: 'dag-swatch', dataset: { evidence: 'observed' } }));
    observed.appendChild(
      document.createTextNode(
        'observed transition — one stage’s last event is timestamped after another’s'
      )
    );
    const order = h('span');
    order.appendChild(h('span', { class: 'dag-swatch', dataset: { evidence: 'order' } }));
    order.appendChild(document.createTextNode('declared order only — no transition evidence in the log'));
    legend.append(observed, order);
    dagHost.appendChild(legend);

    if (evidence === 'unavailable' || evidence === 'order') {
      dagHost.appendChild(
        note(
          'The API exposes no `prerequisites` field on /api/stages or /api/run, so a ' +
            'declared edge cannot be drawn from the server. These connectors are the ' +
            "server's stage order, drawn dashed, and they are not claims about which " +
            'stage needs which. Proposal for BACKEND: add `prerequisites: [stage]` per ' +
            'stage so the DAG can show needs rather than order.',
          'panel-note'
        )
      );
    }
  }

  function renderList() {
    clear(listHost);
    listHost.appendChild(
      h(
        'div',
        { class: 'panel-head' },
        h('h3', null, 'Stage by stage'),
        h('span', { class: 'faint' }, 'state, rows, output path and the reason when there is one')
      )
    );
    if (stages.length === 0) {
      listHost.appendChild(absentSpan('the server returned no stage rows', ABSENT.NOT_PRODUCED));
      return;
    }
    const list = h('ul', { class: 'stage-list' });
    for (const stage of stages) {
      const row = h('li', { class: 'stage-row', dataset: { stage: stage.name } });
      const detailId = `stage-detail-${stage.name}`;
      const button = h(
        'button',
        {
          type: 'button',
          class: 'stage-button',
          'aria-expanded': stage.name === selected ? 'true' : 'false',
          'aria-controls': detailId,
        },
        h('span', 'stage-order', stage.order === undefined ? '—' : String(stage.order)),
        h(
          'span',
          null,
          h('span', { class: 'stage-name' }, stage.name),
          h('span', { class: 'stage-spec' }, specNameOf(stage))
        ),
        h('span', { class: 'stage-rows' }, rowsLabel(stage))
      );
      button.addEventListener('click', () => select(stage.name));
      row.appendChild(button);
      row.appendChild(
        h('div', { class: 'stage-row-state' }, ctx.badges.render(stage), reasonNode(stage))
      );
      if (stage.name === selected) {
        row.appendChild(h('div', { class: 'stage-detail', id: detailId }, detailNode()));
      }
      list.appendChild(row);
    }
    listHost.appendChild(list);
  }

  function specNameOf(stage) {
    if (!stage.spec_name || stage.spec_name === stage.name) {
      return h('span', { class: 'faint' }, 'no spec name is recorded for this stage');
    }
    // `spec_number` arrives as the string "1" .. "15" plus "6a": it is spec.md's
    // row number, so it is labelled as a row and not inflated into "D6a".
    const row_ = stage.spec_number ? `spec.md row ${stage.spec_number}` : 'no spec.md row number';
    return `${stage.spec_name} — ${row_}. The code name is the one shown because the manifest, the event log and every table key use it.`;
  }

  function rowsLabel(stage) {
    if (stage.present) return `${stage.n_rows} rows`;
    return ABSENT.NOT_PRODUCED;
  }

  function reasonNode(stage) {
    const badge = stage.badge && stage.badge.key ? stage.badge : null;
    if (badge && badge.key === 'completed') {
      return h('span', { class: 'faint' }, 'ran and its declared output is present and readable');
    }
    const reason = (badge && badge.reason) || stage.reason;
    if (!reason) {
      return absentSpan(
        'the server classified this stage without giving a reason, which is itself a defect in the record'
      );
    }
    return h('div', { class: 'panel-note' }, reason);
  }

  function detailNode() {
    if (!stageDetail) return h('p', { class: 'panel-note' }, 'Loading this stage…');
    const d = stageDetail;
    const box = h('div');

    box.appendChild(kv([
      ['state', ctx.badges.render(d)],
      ['contract header', d.declared_header && d.declared_header.length
        ? h('span', { class: 'mono' }, `${d.declared_header.length} columns: ${d.declared_header.join(', ')}`)
        : absentSpan('contracts.py declares no header for this stage')],
      ['header the file carries', d.actual_header && d.actual_header.length
        ? h('span', { class: 'mono' }, `${d.actual_header.length} columns: ${d.actual_header.join(', ')}`)
        : absentSpan(d.actual_header_reason || 'the file carries no readable header')],
      ['output', d.path ? h('code', null, d.path) : absentSpan('no path was reported for this stage')],
      ['rows', typeof d.n_rows === 'number' ? String(d.n_rows) : absentSpan('the table is absent, so there is no row count')],
      [
        'folded-step tables',
        d.internal_tables && d.internal_tables.length
          ? h(
              'ul',
              { class: 'gap-list' },
              d.internal_tables.map((t) =>
                h(
                  'li',
                  null,
                  h('code', null, t.key),
                  ' — ',
                  t.present ? `${t.n_rows} rows` : ABSENT.NOT_PRODUCED,
                  t.present ? null : ` (${t.reason})`,
                  ' ',
                  foldedButton(t.key)
                )
              )
            )
          : absentSpan('this stage owns no folded-step table'),
      ],
      [
        'spec.md name',
        d.spec_name
          ? h('code', null, `${d.spec_name}${d.spec_number ? ` — spec.md row ${d.spec_number}` : ''}`)
          : absentSpan('no spec name is recorded for this stage'),
      ],
    ]));

    const drift = Array.isArray(d.header_drift) ? d.header_drift : [];
    box.appendChild(h('h4', null, 'Header drift'));
    if (drift.length === 0) {
      box.appendChild(
        note(
          'None: the file carries exactly the columns the contract declares. Drift is ' +
            'reported as information, never used as a filter — a file whose columns ' +
            'have drifted is still the run’s record.',
          'panel-note'
        )
      );
    } else {
      const ul = h('ul', { class: 'gap-list' });
      for (const item of drift) {
        ul.appendChild(h('li', null, h('code', null, item.column), ` — ${item.kind.replace(/_/g, ' ')}`));
      }
      box.appendChild(ul);
    }

    box.appendChild(h('h4', null, 'Events for this stage'));
    const events = Array.isArray(d.events) ? d.events : [];
    if (events.length === 0) {
      box.appendChild(
        absentSpan(
          'the event log holds no events naming this stage. That is not "the stage did nothing": ' +
            'a stage that ran before the log existed, or in another entry point, leaves no event.',
          ABSENT.NOT_PRODUCED
        )
      );
    } else {
      box.appendChild(eventTable(events));
      box.appendChild(runtimeNode(events));
    }

    box.appendChild(h('h4', null, 'Logs'));
    box.appendChild(
      note(
        'No per-stage log file is contracted anywhere in contracts.py or the Snakefile, so ' +
          'there is no per-stage log to open. The only log files the API serves are the ' +
          'tripwire and watcher tails on the Provenance page, which name every path they ' +
          'probed.',
        'panel-note'
      )
    );

    box.appendChild(h('h4', null, 'Outputs'));
    const outputs = run && run.outputs ? run.outputs : {};
    const mine = Object.entries(outputs).filter(([key]) => key === d.name || key.startsWith(`${d.name}_`));
    if (mine.length === 0) {
      box.appendChild(
        absentSpan(
          run && run.manifest_writer === 'provenance'
            ? 'run_manifest.json was written by write_provenance.py before the run, and that writer records no outputs'
            : 'the manifest carries no output entry for this stage',
          ABSENT.NOT_REPORTED
        )
      );
    } else {
      box.appendChild(kv(mine.map(([key, path]) => [key, h('code', null, path)])));
    }
    return box;
  }

  function foldedButton(key) {
    const button = h('button', { type: 'button', class: 'ghost' }, 'rows');
    button.addEventListener('click', (event) => {
      event.stopPropagation();
      void showTable(key);
    });
    return button;
  }

  async function showTable(key) {
    clear(detailHost);
    const box = h('div', { class: 'panel' });
    box.appendChild(h('h3', null, `Folded-step table: ${key}`));
    const body = h('p', { class: 'panel-note' }, 'Loading…');
    box.appendChild(body);
    detailHost.appendChild(box);
    try {
      const payload = await ctx.api.get(`/tables/${encodeURIComponent(key)}`, { limit: 50 }, abort.signal);
      if (disposed) return;
      clear(body);
      const meta = payload.meta || {};
      const basis = basisOf(payload);
      if (basis) ctx.banner.mount(body, null, basis);
      body.appendChild(
        note(
          `${meta.total ?? ABSENT.NOT_REPORTED} rows in the whole artefact; ${(payload.items || []).length} shown. ` +
            `artefact: ${basis && basis.artefact ? basis.artefact : key}.`,
          'table-meta'
        )
      );
      if (!payload.present) {
        clear(body);
        if (basis) ctx.banner.mount(body, null, basis);
        body.appendChild(absentSpan(payload.reason || 'the table is not present', ABSENT.NOT_PRODUCED));
        body.appendChild(h('pre', null, payload.path || ''));
        return;
      }
      const columns = Array.isArray(payload.header) ? payload.header : Object.keys((payload.items || [])[0] || {});
      body.appendChild(dataTable(columns, payload.items, 'the response carries no rows'));
    } catch (error) {
      if (error && error.name === 'AbortError') return;
      clear(body);
      body.appendChild(message('error', 'This table could not be read', note(String(error.message || error))));
    }
  }

  /**
   * The stage's observed runtime, from the event log's own timestamps.
   *
   * The API exposes no runtime field, so this is computed from the two
   * timestamps the log already carries and is labelled as exactly that. When
   * the log holds fewer than two timestamps it reads `runtime not recorded`
   * with the reason — an estimate derived from a per-stage average would be a
   * measurement nobody made (UI-D2).
   */
  function runtimeNode(events) {
    const starts = events.filter((e) => e.event === 'start').map((e) => Date.parse(e.t)).filter(Number.isFinite);
    const ends = events
      .filter((e) => e.event === 'done' || e.event === 'fail')
      .map((e) => Date.parse(e.t))
      .filter(Number.isFinite);
    if (starts.length === 0 || ends.length === 0) {
      return h(
        'p',
        { class: 'panel-note' },
        absentSpan(
          starts.length === 0
            ? 'the log holds no start event for this stage, so no elapsed time can be computed'
            : 'the log holds start events but no terminal event, so the span is still open',
          'runtime not recorded'
        )
      );
    }
    const first = Math.min(...starts);
    const last = Math.max(...ends);
    const seconds = (last - first) / 1000;
    const fmt =
      seconds < 90
        ? `${seconds.toFixed(seconds < 10 ? 1 : 0)} s`
        : `${(seconds / 60).toFixed(1)} min`;
    return h(
      'p',
      { class: 'panel-note' },
      `Observed span: ${fmt} between ${new Date(first).toISOString()} and ${new Date(last).toISOString()}, ` +
        `over ${starts.length} start and ${ends.length} terminal events for this stage. This is the log's ` +
        'own timestamps, not a server-computed duration: the API exposes no per-stage runtime field.'
    );
  }

  function eventTable(events) {
    const rows = events
      .slice()
      .reverse()
      .slice(0, 200)
      .map((e) => ({
        t: e.t,
        event: e.event,
        rule: e.stage,
        code_stage: e.mapped_stage || `${e.stage} maps to no code stage name`,
        sample: e.sample,
        line: e.line,
      }));
    return dataTable(['t', 'event', 'rule', 'code_stage', 'sample', 'line'], rows, 'no event rows');
  }

  async function select(name, scroll) {
    selected = name;
    stageDetail = null;
    renderDag();
    renderList();
    clear(detailHost);
    if (scroll) {
      const row = container.querySelector(`.stage-row[data-stage="${CSS.escape(name)}"]`);
      if (row) row.scrollIntoView({ block: 'nearest' });
    }
    try {
      stageDetail = await ctx.api.get(`/stages/${encodeURIComponent(name)}`, null, abort.signal);
      if (disposed) return;
      renderList();
    } catch (error) {
      if (error && error.name === 'AbortError') return;
      stageDetail = null;
      renderList();
      detailHost.appendChild(
        message('error', `Stage "${name}" could not be read`, note(String(error.message || error)))
      );
    }
  }

  // Live state: a running stage's badge must move without a reload, and the
  // one subscription is the shell's (§7 `ctx.on`).
  disposers.push(
    ctx.on.runState((snapshotOrNull) => {
      if (disposed) return;
      snapshot = snapshotOrNull;
      if (!snapshot || !Array.isArray(snapshot.stages) || snapshot.stages.length === 0) return;
      // The snapshot's stage rows are the same payload /api/stages returns, so
      // they replace the fetched copy rather than being merged into it.
      stages = snapshot.stages;
      renderCounts();
      renderDag();
      renderList();
    })
  );

  (async () => {
    try {
      await loadAll();
    } catch (error) {
      if (error && error.name === 'AbortError') return;
      clear(listHost);
      listHost.appendChild(
        message(
          'error',
          'The stage list could not be read',
          note(String(error.message || error)),
          error.probes ? probeList(error.probes) : null,
          error.serverMessage ? serverMessage(error.serverMessage) : null
        )
      );
    }
  })();

  return () => {
    disposed = true;
    abort.abort();
    for (const fn of disposers.splice(0)) {
      try {
        fn();
      } catch (error) {
        console.error('[pipeline] a subscription disposer threw', error);
      }
    }
  };
}

/* ------------------------------------------------------- the run monitor */

const REAL_MODE_NOTE =
  'A REAL run needs all four of: the phrase typed exactly, the bigmachine ' +
  'overlay, the reuse mode shown, and PIPELINE_ALLOW_REAL_MODE set on the ' +
  'child process only. This build has no endpoint that starts one: a launch ' +
  'is a human command run in the repository.';

export function mountMonitor(container, ctx) {
  const abort = new AbortController();
  const disposers = [];
  let disposed = false;

  let snapshot = null;
  let logLines = [];
  let counts = null;
  let provenance = null;
  let capabilities = null;
  let launcherStatus = null;
  let preflight = null;
  let transportNote = 'the transport indicator in the header says which mechanism this page is reading';

  container.appendChild(
    h(
      'div',
      { class: 'page-head' },
      h(
        'div',
        { class: 'page-head-text' },
        h('h2', null, 'Run monitor'),
        h(
          'p',
          { class: 'page-note' },
          'Live status from the event log, over server-sent events with a polling ' +
            'fallback. The transport in the header says which one is in use, so a ' +
            'screenshot of this page says so too.'
        )
      )
    )
  );

  const bannerHost = h('div');
  const statsHost = h('div', { class: 'panel' });
  const transportHost = h('div', { class: 'panel' });
  const logHost = h('div', { class: 'panel' });
  const hostHost = h('div', { class: 'panel' });
  const launcherHost = h('div', { class: 'panel' });
  container.append(bannerHost, statsHost, transportHost, logHost, hostHost, launcherHost);

  /* -- rendering ---------------------------------------------------- */

  function renderStats() {
    clear(statsHost);
    statsHost.appendChild(
      h('div', { class: 'panel-head' }, h('h3', null, 'Stage states'), h('span', { class: 'faint' }, 'from the server, over all stages'))
    );
    const stages = snapshot && Array.isArray(snapshot.stages) ? snapshot.stages : [];
    if (stages.length === 0) {
      statsHost.appendChild(
        absentSpan(
          snapshot ? 'no stage rows arrived with the event state' : 'the event state has not arrived yet',
          ABSENT.NOT_REPORTED
        )
      );
      return;
    }
    // Counts here are the server's tally (`/api/run/counts`), never a count of
    // what is on screen (UI-D7). If the tally has not arrived, nothing is
    // counted here either.
    const row = h('div', { class: 'stat-row' });
    if (counts) {
      for (const key of [
        'n_stages_completed',
        'n_stages_running',
        'n_stages_failed',
        'n_stages_refused',
        'n_stages_not_assessed',
        'n_stages_not_run',
      ]) {
        if (typeof counts[key] !== 'number') continue;
        const stat = h('div', { class: 'stat' });
        stat.appendChild(h('span', { class: 'stat-value' }, String(counts[key])));
        stat.appendChild(h('span', { class: 'stat-label' }, key.replace(/^n_stages_/, '').replace(/_/g, ' ')));
        row.appendChild(stat);
      }
    } else {
      row.appendChild(absentSpan('the server tally has not arrived; nothing is counted in the browser', ABSENT.NOT_REPORTED));
    }
    statsHost.appendChild(row);

    const stages2 = h('div', { class: 'chip-row' });
    for (const stage of stages) {
      const chip = h('span', { class: 'chip' });
      chip.appendChild(h('code', null, stage.name));
      chip.appendChild(ctx.badges.render(stage));
      stages2.appendChild(chip);
    }
    statsHost.appendChild(stages2);

    const basis = counts ? basisOf(counts, 'stages') : { n: stages.length, artefact: 'stage_event_state' };
    clear(bannerHost);
    ctx.banner.mount(bannerHost, null, basis);

    const totals = snapshot && snapshot.totals ? snapshot.totals : null;
    if (totals) {
      statsHost.appendChild(
        note(
          `Event log totals: ${totals.start || 0} start, ${totals.done || 0} done, ${totals.fail || 0} fail. ` +
            'These are the log’s own counts, as replayed by the server.',
          'table-meta'
        )
      );
    }
  }

  function renderTransportPanel() {
    clear(transportHost);
    transportHost.appendChild(h('h3', null, 'Transport'));
    transportHost.appendChild(
      kv([
        ['in use', h('strong', null, transportNote)],
        [
          'cursor',
          snapshot && typeof snapshot.cursor === 'number'
            ? h('code', null, String(snapshot.cursor))
            : absentSpan('no cursor has arrived', ABSENT.NOT_REPORTED),
        ],
        [
          'events replayed',
          snapshot && typeof snapshot.n_events === 'number'
            ? h('b', null, String(snapshot.n_events))
            : absentSpan('no count has arrived', ABSENT.NOT_REPORTED),
        ],
        [
          'lines skipped',
          snapshot && typeof snapshot.n_lines_skipped === 'number'
            ? h('b', null, String(snapshot.n_lines_skipped))
            : absentSpan('the server has not reported a skip count', ABSENT.NOT_REPORTED),
        ],
      ])
    );

    const warnings = snapshot && Array.isArray(snapshot.skipped_warnings) ? snapshot.skipped_warnings : [];
    if (warnings.length > 0) {
      transportHost.appendChild(h('h4', null, 'Lines the reader skipped, verbatim'));
      transportHost.appendChild(
        note(
          'The log carries exactly four fields and three event values. A line that ' +
            'breaks that is skipped with a warning naming its line number — never ' +
            'coerced into the shape.'
        )
      );
      const ul = h('ul', { class: 'gap-list' });
      for (const warning of warnings.slice(0, 50)) ul.appendChild(h('li', { class: 'mono' }, warning));
      transportHost.appendChild(ul);
    }

    const unmapped = snapshot && Array.isArray(snapshot.unmapped_stages) ? snapshot.unmapped_stages : [];
    if (unmapped.length > 0) {
      transportHost.appendChild(h('h4', null, 'Rule names that map to no stage'));
      const ul = h('ul', { class: 'gap-list' });
      for (const item of unmapped.slice(0, 50)) ul.appendChild(h('li', { class: 'mono' }, item));
      transportHost.appendChild(ul);
    }

    const probes = snapshot && Array.isArray(snapshot.log_probes) ? snapshot.log_probes : [];
    if (probes.length > 0) {
      transportHost.appendChild(h('h4', null, 'Where the log was looked for'));
      transportHost.appendChild(probeList(probes));
    }
  }

  function renderLog() {
    clear(logHost);
    logHost.appendChild(
      h('div', { class: 'panel-head' }, h('h3', null, 'Log tail'), h('span', { class: 'faint' }, 'newest last'))
    );
    if (snapshot && snapshot.transport === 'none') {
      logHost.appendChild(
        message(
          'error',
          'Event stream unavailable',
          h('p', null, `The server's reason: ${snapshot.reason || 'it did not say'}.`),
          note(
            'The run summary above is still served from the manifest. This panel ' +
              'deliberately shows no log rather than an empty one: an empty log ' +
              'reads as "nothing has happened", which is a different claim.'
          )
        )
      );
      return;
    }
    if (logLines.length === 0) {
      logHost.appendChild(
        absentSpan(
          snapshot
            ? 'no event line has arrived on this transport. The log may be empty, may not exist, or may carry only lines the reader skipped — the panel above says which.'
            : 'the event state has not arrived yet',
          ABSENT.NOT_PRODUCED
        )
      );
      return;
    }
    const pre = h('div', { class: 'log-tail' });
    for (const line of logLines) {
      const row = h('span', { class: 'log-line', dataset: { event: line.event } });
      row.appendChild(h('span', 'log-time', `${line.t} `));
      row.appendChild(h('span', 'log-event', `${line.event.toUpperCase()} `));
      row.appendChild(
        h('span', null, `rule=${line.stage} sample=${line.sample}${line.mapped_stage ? '' : ' (maps to no stage)'}`)
      );
      pre.appendChild(row);
    }
    logHost.appendChild(pre);
    logHost.appendChild(
      note(
        `${logLines.length} line(s) held in this tab. A sample id is shown verbatim as the ` +
          'log carries it; a synthetic run uses fixture ids and a REAL run uses real accessions.',
        'table-meta'
      )
    );
  }

  function renderHost() {
    clear(hostHost);
    hostHost.appendChild(h('h3', null, 'Host'));
    if (!provenance) {
      hostHost.appendChild(absentSpan('/api/provenance has not answered', ABSENT.NOT_REPORTED));
      return;
    }
    const host = provenance.host || {};
    hostHost.appendChild(
      kv([
        [
          'free disk on the results volume',
          typeof host.disk_free_bytes === 'number'
            ? h('b', null, formatBytes(host.disk_free_bytes))
            : absentSpan(
                host.disk_free_bytes === undefined
                  ? 'the field is absent from the response'
                  : 'the volume could not be stat’d; 0 bytes free would be a claim about a filesystem nobody measured'
              ),
        ],
        [
          'volume total',
          typeof host.disk_total_bytes === 'number' ? h('b', null, formatBytes(host.disk_total_bytes)) : absentSpan('not reported'),
        ],
        [
          'load average',
          Array.isArray(host.load_average) && host.load_average.length
            ? h('code', null, host.load_average.map((v) => v.toFixed(2)).join(' · '))
            : absentSpan('this platform reports no load average; that is not a load of zero'),
        ],
      ])
    );
  }

  function renderLauncher() {
    clear(launcherHost);
    launcherHost.appendChild(
      h('div', { class: 'panel-head' }, h('h3', null, 'Launcher'))
    );
    if (capabilities === null) {
      launcherHost.appendChild(absentSpan('/api/launcher/capabilities has not answered', ABSENT.NOT_REPORTED));
      return;
    }

    // The panel is VISIBLE ONLY when the server says launching is enabled.
    // Otherwise it says why it is hidden — the boundary, not a blank.
    if (!capabilities.allow_launch_enabled) {
      const box = h('div', { class: 'launcher-locked' });
      box.appendChild(h('h4', null, 'Launcher disabled'));
      box.appendChild(h('p', { class: 'panel-note' }, capabilities.reason));
      box.appendChild(h('h4', null, 'What the launcher may not do here'));
      const ul = h('ul', { class: 'gap-list' });
      for (const item of capabilities.forbidden_actions || []) ul.appendChild(h('li', null, item));
      box.appendChild(ul);
      box.appendChild(h('h4', null, 'What it may do here'));
      const ul2 = h('ul', { class: 'gap-list' });
      for (const item of capabilities.allowed_actions || []) ul2.appendChild(h('li', null, item));
      box.appendChild(ul2);
      box.appendChild(h('h4', null, 'The REAL gates, stated in advance'));
      box.appendChild(phraseNote(capabilities));
      launcherHost.appendChild(box);
      renderRunState();
      return;
    }

    launcherHost.appendChild(note(capabilities.reason));

    const chips = h('div', { class: 'chip-row' });
    for (const tool of capabilities.expensive_tools || []) {
      chips.appendChild(h('span', { class: 'chip', dataset: { flag: 'true' } }, `⚠ ${tool}`));
    }
    launcherHost.appendChild(
      note(
        'Flagged before anything starts, not after: these dominate a run’s wall clock.',
        'panel-note'
      )
    );
    launcherHost.appendChild(chips);

    const controls = h('div', { class: 'toolbar' });
    const dryButton = h('button', { type: 'button', class: 'primary' }, 'Show the dry-run plan');
    dryButton.addEventListener('click', () => void runPreflight('dry_run'));
    controls.appendChild(dryButton);
    launcherHost.appendChild(controls);

    if (preflight) renderPreflight();

    launcherHost.appendChild(h('h4', null, 'Start a REAL run'));
    launcherHost.appendChild(realForm(capabilities));
    launcherHost.appendChild(note(REAL_MODE_NOTE));
    renderRunState();
  }

  /**
   * Whether a run is in progress, and the per-isolate timings.
   *
   * Shown whether or not the launcher is enabled: both are facts about the run
   * rather than about permission to start one, and the timings degrade to
   * `timings not recorded` either way. The timings are never derived from the
   * cohort size — an invented per-isolate estimate is a measurement nobody made.
   */
  function renderRunState() {
    if (!launcherStatus) {
      launcherHost.appendChild(
        absentSpan('/api/launcher/status has not answered', ABSENT.NOT_REPORTED)
      );
      return;
    }
    launcherHost.appendChild(h('h4', null, 'Is a run in progress'));
    const timings = launcherStatus.per_isolate_timings;
    launcherHost.appendChild(
      kv([
        [
          'running',
          typeof launcherStatus.running === 'boolean'
            ? h('strong', null, launcherStatus.running ? 'yes' : 'no')
            : absentSpan(
                'no event log was found, so this is unknown — it is not "nothing is running"'
              ),
        ],
        [
          'per-isolate timings',
          timings && Object.keys(timings).length
            ? h('b', null, `${Object.keys(timings).length} isolates`)
            : ctx.badges.render(ctx.badges.timings(launcherStatus.per_isolate_timings_label)),
        ],
        ['timings source', h('code', null, launcherStatus.timings_source || 'none')],
        [
          'last event',
          launcherStatus.last_event
            ? h('code', null, `${launcherStatus.last_event.t} ${launcherStatus.last_event.event} ${launcherStatus.last_event.stage}`)
            : absentSpan('no event was replayed from the log'),
        ],
      ])
    );
    const probes = launcherStatus.runbook_probes || [];
    if (probes.length) {
      launcherHost.appendChild(h('h4', null, 'Where a runbook was looked for'));
      launcherHost.appendChild(probeList(probes));
    }
  }

  function phraseNote(caps) {
    const box = h('div');
    box.appendChild(
      kv([
        ['typed phrase', h('code', null, caps.real_phrase || 'the server did not name one')],
        ['machine overlay', h('code', null, caps.required_overlay || 'the server did not name one')],
        [
          'child environment',
          h('code', null, 'PIPELINE_ALLOW_REAL_MODE'),
        ],
        ['reuse mode', 'shown from configuration before anything would run'],
      ])
    );
    return box;
  }

  function realForm(caps) {
    const box = h('div');
    const phrase = caps.real_phrase;
    const input = h('input', {
      type: 'text',
      class: 'phrase-input',
      id: 'real-phrase',
      placeholder: phrase || 'run real samples',
      autocomplete: 'off',
      autocapitalize: 'off',
      spellcheck: false,
      'aria-label': 'Type the confirmation phrase exactly',
    });
    const overlayInput = h('input', {
      type: 'text',
      class: 'phrase-input',
      autocomplete: 'off',
      spellcheck: false,
      'aria-label': 'Machine overlay',
    });
    // Assigned, not set as an attribute: an attribute would only set the input's
    // *default* value, and the overlay that is sent is the one in the field.
    overlayInput.value = caps.required_overlay || '';
    const button = h('button', { type: 'button', class: 'primary', disabled: true }, 'Check the REAL gates');
    const status = h('p', { class: 'phrase-exact' });
    let gateHost = h('div');

    // Exact equality. No trimming, no case folding, no prefix match: the value
    // that would be sent is the value that was typed, character for character,
    // so " run real samples", "Run Real Samples" and "run real" are all refused
    // by this button rather than being silently normalised into consent.
    const exact = () => input.value === phrase;
    const sync = () => {
      button.disabled = !exact();
      clear(status);
      status.appendChild(
        exact()
          ? h('span', { class: 'probe-ok' }, 'exact match — the button is enabled')
          : h(
              'span',
              { class: 'probe-miss' },
              `not an exact match. The phrase is typed character for character; ` +
                'surrounding spaces, capitals and a shorter phrase are all refused.'
            )
      );
    };
    input.addEventListener('input', sync);
    input.addEventListener('keydown', (event) => {
      if (event.key === 'Enter' && exact()) {
        event.preventDefault();
        button.click();
      }
    });
    button.addEventListener('click', () => void runPreflight('real', input.value, overlayInput.value));
    sync();

    box.appendChild(
      kv([
        ['type this exactly', input],
        ['machine overlay', overlayInput],
        ['', button],
      ])
    );
    box.appendChild(status);
    box.appendChild(gateHost);
    box.appendChild(gateHost);
    return box;
  }

  function renderPreflight() {
    const existing = launcherHost.querySelector('.preflight-block');
    if (existing) existing.remove();
    const box = h('div', { class: 'preflight-block' });
    box.appendChild(h('h4', null, 'Preflight — starts nothing'));

    const checks = Array.isArray(preflight.checks) ? preflight.checks : [];
    if (checks.length === 0) {
      box.appendChild(absentSpan('the server returned no checks', ABSENT.NOT_REPORTED));
    } else {
      const list = h('ul', { class: 'gate-list' });
      for (const check of checks) list.appendChild(gateRow(check));
      box.appendChild(list);
    }

    box.appendChild(h('h4', null, 'Plan'));
    if (!preflight.plan) {
      box.appendChild(
        absentSpan(
          'no plan was returned: the launcher is disabled, or the dry-run was not requested. Nothing would have been executed either way.',
          ABSENT.NOT_PRODUCED
        )
      );
    } else {
      box.appendChild(
        kv([
          ['kind', h('code', null, preflight.plan.kind)],
          ['command that would be run', h('code', null, (preflight.plan.argv || []).join(' '))],
          ['executed?', h('strong', null, 'no — the server ran nothing')],
          ['reuse mode', preflight.reuse_mode ? h('code', null, preflight.reuse_mode) : absentSpan('reuse_tool_output could not be read from config/science.yaml')],
          [
            'per-isolate timings',
            preflight.plan.per_isolate_timings_label
              ? ctx.badges.render(
                  preflight.plan.per_isolate_timings && Object.keys(preflight.plan.per_isolate_timings).length
                    ? 'completed'
                    : ctx.badges.timings(preflight.plan.per_isolate_timings_label)
                )
              : absentSpan('the plan carries no timings field', ABSENT.NOT_REPORTED),
          ],
        ])
      );
      box.appendChild(note(preflight.plan.note || ''));
    }

    const realGates = Array.isArray(preflight.real_gates) ? preflight.real_gates : [];
    if (realGates.length > 0) {
      box.appendChild(h('h4', null, 'The four gates a REAL run must pass'));
      const list = h('ul', { class: 'gate-list' });
      for (const gate of realGates) list.appendChild(gateRow(gate));
      box.appendChild(list);
      box.appendChild(
        note(
          preflight.real_gates_all_passed
            ? 'Every gate the server can check reports passed. It still has no endpoint that starts a run.'
            : 'At least one gate is open, so no REAL run could begin from here even if one existed.'
        )
      );
    }
    if (preflight.reason) box.appendChild(message('warn', 'Why no run starts', h('p', null, preflight.reason)));
    launcherHost.appendChild(box);
  }

  function gateRow(gate) {
    const row = h('li', { class: 'gate', dataset: { passed: String(Boolean(gate.passed)) } });
    row.appendChild(h('span', 'gate-glyph', gate.passed ? '✓' : '✗'));
    const body = h('div');
    body.appendChild(h('span', 'gate-name', gate.name));
    body.appendChild(h('div', 'gate-detail', gate.detail || ''));
    row.appendChild(body);
    return row;
  }

  async function runPreflight(mode, phrase, overlay) {
    const body = mode === 'real' ? { mode, phrase, overlay } : { mode };
    try {
      const payload = await ctx.api.post('/launcher/preflight', body, abort.signal);
      if (disposed) return;
      preflight = payload;
      renderLauncher();
    } catch (error) {
      if (error && error.name === 'AbortError') return;
      preflight = null;
      launcherHost.appendChild(
        message('error', 'The preflight request failed', note(String(error.message || error)))
      );
    }
  }

  /* -- subscriptions and loads --------------------------------------- */

  disposers.push(
    ctx.on.runState((state) => {
      if (disposed || !state) return;
      snapshot = state;
      if (state.transport) {
        transportNote =
          state.transport === 'sse'
            ? 'server-sent events (SSE): a frame arrived, so the 3 s fallback timer was stopped'
            : state.transport === 'poll'
              ? 'polling /api/events/state every 2 s, carrying the cursor'
              : `unavailable — ${state.reason || 'the server did not say why'}`;
      }
      renderStats();
      renderTransportPanel();
      renderLog();
    })
  );

  disposers.push(
    ctx.on.delta((delta, meta) => {
      if (disposed || !Array.isArray(delta)) return;
      for (const event of delta) logLines.push(event);
      // A bounded tail: a 900-isolate run emits tens of thousands of lines and
      // an unbounded array is a leak with a render loop attached.
      if (logLines.length > 500) logLines = logLines.slice(-500);
      renderLog();
      void meta;
    })
  );

  (async () => {
    try {
      const [countsResponse, provenanceResponse, capabilitiesResponse, statusResponse] =
        await Promise.all([
          ctx.api.get('/run/counts', null, abort.signal),
          ctx.api.get('/provenance', null, abort.signal),
          ctx.api.get('/launcher/capabilities', null, abort.signal),
          ctx.api.get('/launcher/status', null, abort.signal),
        ]);
      if (disposed) return;
      counts = countsResponse;
      provenance = provenanceResponse;
      capabilities = capabilitiesResponse;
      launcherStatus = statusResponse;
      renderStats();
      renderHost();
      renderLauncher();
    } catch (error) {
      if (error && error.name === 'AbortError') return;
      launcherHost.appendChild(
        message(
          'error',
          'The monitor could not read its endpoints',
          note(String(error.message || error))
        )
      );
    }
  })();

  return () => {
    disposed = true;
    abort.abort();
    for (const fn of disposers.splice(0)) {
      try {
        fn();
      } catch (error) {
        console.error('[monitor] a subscription disposer threw', error);
      }
    }
  };
}

function formatBytes(bytes) {
  const units = ['B', 'KiB', 'MiB', 'GiB', 'TiB'];
  let value = Number(bytes);
  let unit = 0;
  while (value >= 1024 && unit < units.length - 1) {
    value /= 1024;
    unit += 1;
  }
  return `${value.toFixed(unit === 0 ? 0 : 1)} ${units[unit]}`;
}
