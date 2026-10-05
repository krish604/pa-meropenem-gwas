// Wiring. Owns the connection, the selection, and the view state.

import { Stream, getJSON, startClock } from './data.js';
import { NetworkView } from './network.js';
import { STATE_LIST, STATES, isLive } from './palette.js';
import * as detail from './detail.js';
import * as metrics from './metrics.js';
import { byId, on, clamp } from './util.js';
import { count, EM_DASH } from './format.js';

const app = {
  stream: null,
  network: null,
  snapshot: null,
  selected: null,
  frozen: false,
};

// ── boot ──────────────────────────────────────────────────────────

async function boot() {
  startClock();
  detail.initDetail();

  app.network = new NetworkView(byId('canvas-host'), byId('network-canvas'), {
    onSelect: (stage) => selectStage(stage),
    onHover: showTip,
  });
  app.network.start();

  bindControls();
  renderLegend();
  renderIdleLegend();

  // Structure first, so the canvas has something real to draw immediately.
  try {
    const structure = await getJSON('/api/structure');
    app.network.setGraph(structure);
  } catch (err) {
    setIdleMessage(`Could not read the pipeline structure: ${err}`);
  }

  app.stream = new Stream({
    runKey: new URLSearchParams(location.search).get('run') || '',
    onSnapshot: applySnapshot,
    onEvent: onEvent,
    onStatus: setStatus,
    onHello: () => {},
  });
  await app.stream.refresh();
  app.stream.connect();

  setInterval(() => {
    detail.tickElapsed();
    metrics.metricsSnapshot();
  }, 1000);

  window.addEventListener('resize', () => {
    app.network.resize();
    if (app.network.fitted) app.network.fit();
  });
}

// ── snapshots ─────────────────────────────────────────────────────

function applySnapshot(snapshot) {
  app.snapshot = snapshot;
  app.network.setState(snapshot);
  metrics.update(snapshot);

  const store = snapshot.store || {};
  setText('ctx-kind', `Genomic Analysis Pipeline · run ${snapshot.run_key || '—'}`);
  setText('network-sub', describeNetwork(snapshot));

  const idleOverlay = byId('idle-overlay');
  if (idleOverlay) idleOverlay.hidden = !snapshot.idle;
  if (snapshot.idle) {
    setIdleMessage(snapshot.exists && snapshot.tasks && snapshot.tasks.length
      ? 'Nothing is running. Recorded state is shown above.'
      : (store.present
        ? 'The execution store has no rows for this run yet.'
        : 'No execution store found. Host metrics below are live.'));
  }

  // Follow the running work unless the operator has picked something.
  if (!app.selected) {
    const running = (snapshot.tasks || []).find((t) => isLive(t.state));
    if (running) showTask(running);
  } else {
    const fresh = (snapshot.tasks || []).find(
      (t) => t.stage === app.selected.stage && t.subject === app.selected.subject,
    );
    if (fresh) showTask(fresh);
  }
}

function describeNetwork(snapshot) {
  const counts = snapshot.counts || {};
  const total = Object.values(counts).reduce((a, b) => a + b, 0);
  if (!total) return 'Awaiting execution state';
  const bits = [];
  if (counts.RUNNING) bits.push(`${counts.RUNNING} running`);
  if (counts.SUCCEEDED) bits.push(`${counts.SUCCEEDED} validated`);
  if (counts.FAILED) bits.push(`${counts.FAILED} failed`);
  if (counts.INVALID) bits.push(`${counts.INVALID} invalid`);
  if (counts.INCOMPLETE) bits.push(`${counts.INCOMPLETE} incomplete`);
  const edges = (snapshot.edges || []).length;
  return `${snapshot.stages ? snapshot.stages.length : 0} stages · ${edges} declared connections · ${bits.join(' · ')}`;
}

function onEvent(ev) {
  detail.recordEvent(ev);
  // A finished task sends a travelling mark along its real downstream edge.
  if (ev.name === 'TASK_COMPLETED' || ev.name === 'TASK_VALIDATED') {
    app.network.markActivity(ev.stage, downstreamOf(ev.stage));
  }
  if (ev.name === 'TASK_STARTED' || ev.name === 'TASK_RETRY') {
    const node = app.network.byStage && app.network.byStage.get(ev.stage);
    if (node) node.flash = 1.0;
  }
}

function downstreamOf(stage) {
  const edge = (app.snapshot && app.snapshot.edges || [])
    .find((e) => e.source === stage && e.kind === 'dependency');
  return edge ? edge.target : null;
}

// ── selection ─────────────────────────────────────────────────────

function selectStage(stage) {
  if (!stage) { app.selected = null; return; }
  app.network.select(stage);
  const task = pickTask(stage);
  if (task) { app.selected = task; showTask(task); }
}

function pickTask(stage) {
  const tasks = (app.snapshot && app.snapshot.tasks) || [];
  const forStage = tasks.filter((t) => t.stage === stage);
  if (!forStage.length) return null;
  // Prefer something an operator would want to look at.
  const order = ['RUNNING', 'RETRYING', 'FAILED', 'INVALID', 'INCOMPLETE', 'SUCCEEDED'];
  for (const state of order) {
    const hit = forStage.find((t) => t.state === state);
    if (hit) return hit;
  }
  return forStage[0];
}

async function showTask(task) {
  let attempts = [];
  try {
    const full = await getJSON(
      `/api/tasks/${encodeURIComponent(task.stage)}/${encodeURIComponent(task.subject || '')}`,
    );
    attempts = full.attempts_history || [];
    task = full;
  } catch (err) {
    // The row from the snapshot is still useful; only history is missing.
  }
  app.selected = task;
  detail.showTask(task, attempts);
}

function showTip(stage, point) {
  const tip = byId('node-tip');
  if (!tip) return;
  if (!stage || !point) { tip.hidden = true; return; }
  const node = app.network.byStage.get(stage);
  if (!node) { tip.hidden = true; return; }
  tip.textContent = '';
  const title = document.createElement('b');
  title.textContent = node.label;
  const rows = [
    ['state', node.state],
    ['tasks', count(node.total)],
    ...Object.entries(node.counts || {}).map(([k, v]) => [k.toLowerCase(), count(v)]),
  ];
  for (const [k, v] of rows) {
    const row = document.createElement('div');
    row.className = 'row';
    row.textContent = `${k}: ${v}`;
    tip.appendChild(row);
  }
  tip.hidden = false;
  const host = byId('canvas-host').getBoundingClientRect();
  const x = clamp(point.x + 14, 8, host.width - tip.offsetWidth - 8);
  const y = clamp(point.y + 14, 8, host.height - tip.offsetHeight - 8);
  tip.style.left = `${x}px`;
  tip.style.top = `${y}px`;
}

// ── controls ──────────────────────────────────────────────────────

function bindControls() {
  on(byId('btn-fit'), 'click', () => app.network.fit());
  on(byId('btn-zoom-in'), 'click', () => app.network.zoomBy(1.2));
  on(byId('btn-zoom-out'), 'click', () => app.network.zoomBy(1 / 1.2));
  on(byId('btn-reconnect'), 'click', async () => {
    app.stream.close();
    await app.stream.refresh();
    app.stream.connect();
  });
  on(byId('btn-freeze'), 'click', (e) => {
    app.frozen = !app.frozen;
    app.network.setFrozen(app.frozen);
    e.currentTarget.classList.toggle('is-on', app.frozen);
  });

  document.querySelectorAll('[data-edges]').forEach((btn) => {
    on(btn, 'click', () => {
      document.querySelectorAll('[data-edges]').forEach((b) => b.classList.toggle('is-active', b === btn));
      app.network.setEdgeFilter(btn.dataset.edges);
    });
  });
  document.querySelectorAll('[data-view2]').forEach((btn) => {
    on(btn, 'click', () => {
      document.querySelectorAll('[data-view2]').forEach((b) => b.classList.toggle('is-active', b === btn));
      app.network.setView(btn.dataset.view2);
    });
  });
  document.querySelectorAll('.tabs .tab[data-view]').forEach((btn) => {
    on(btn, 'click', () => {
      document.querySelectorAll('.tabs .tab[data-view]').forEach((b) => b.classList.toggle('is-active', b === btn));
      // The network is the only implemented view; the others are stated as
      // unavailable rather than silently showing the same thing.
      if (btn.dataset.view !== 'overview') {
        setIdleMessage(`"${btn.textContent.trim()}" is not a separate view in this build — `
          + 'the overview already shows every task, its validation and its log.');
        const overlay = byId('idle-overlay');
        if (overlay) overlay.hidden = false;
      }
    });
  });

  on(window, 'keydown', (e) => {
    if (e.key === 'f') app.network.fit();
    if (e.key === 'Escape') { app.network.select(null); app.selected = null; detail.showTask(null, []); }
  });
}

function renderLegend() {
  const node = byId('legend');
  if (!node) return;
  node.textContent = '';
  for (const key of STATE_LIST) {
    const item = document.createElement('span');
    item.className = 'item';
    const sw = document.createElement('i');
    sw.className = 'swatch';
    sw.style.background = STATES[key].color;
    item.appendChild(sw);
    item.appendChild(document.createTextNode(STATES[key].label));
    node.appendChild(item);
  }
  const dep = document.createElement('span');
  dep.className = 'item';
  const dl = document.createElement('i');
  dl.className = 'line-swatch solid';
  dep.appendChild(dl);
  dep.appendChild(document.createTextNode('dependency'));
  node.appendChild(dep);

  const seq = document.createElement('span');
  seq.className = 'item';
  const sl = document.createElement('i');
  sl.className = 'line-swatch';
  seq.appendChild(sl);
  seq.appendChild(document.createTextNode('execution order only'));
  node.appendChild(seq);
}

function renderIdleLegend() { /* legend is static; kept for symmetry */ }

function setStatus(status, message) {
  const badge = byId('live-badge');
  const text = byId('live-text');
  if (!badge) return;
  badge.classList.toggle('is-down', status === 'down');
  badge.classList.toggle('is-stale', status === 'stale');
  if (text) {
    text.textContent = status === 'down' ? 'RECONNECTING' : status === 'stale' ? 'STALE' : 'LIVE';
  }
  badge.title = message || '';
}

function setIdleMessage(text) {
  const node = byId('idle-msg');
  if (node) node.textContent = text;
}

function setText(id, value) {
  const node = byId(id);
  if (node) node.textContent = value === null || value === undefined ? EM_DASH : String(value);
}

boot().catch((err) => {
  document.body.insertAdjacentHTML(
    'afterbegin',
    `<pre style="color:#ef4444;padding:16px;font:12px monospace">observatory failed to start: ${err}</pre>`,
  );
});
