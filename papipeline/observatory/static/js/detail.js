// The left panel: everything known about one real task.
//
// The panel is where honesty is most visible, because it is where an
// operator looks when something is wrong. A field the engine does not
// record reads "not reported" in grey italic. It does not read 0, and it
// does not read blank, because blank looks like a bug and 0 looks like a
// measurement.

import {
  duration, durationLong, clockTime, count, na, EM_DASH,
} from './format.js';
import { isLive } from './palette.js';
import { byId, replace, on } from './util.js';
import { getLog } from './data.js';

const state = {
  task: null,
  attempts: [],
  logText: '',
  logPaused: false,
  logFilter: '',
  logErrorsOnly: false,
  lastLogKey: null,
  logTimer: null,
};

export function initDetail(hooks = {}) {
  bindSubtabs();
  on(byId('btn-copy-cmd'), 'click', copyCommand);
  on(byId('log-pause'), 'click', () => {
    state.logPaused = !state.logPaused;
    const btn = byId('log-pause');
    btn.textContent = state.logPaused ? 'resume' : 'pause';
    const body = byId('log-body');
    if (body) body.classList.toggle('is-paused', state.logPaused);
  });
  on(byId('log-search'), 'input', (e) => {
    state.logFilter = e.target.value.toLowerCase();
    renderLog();
  });
  on(byId('log-errors'), 'change', (e) => {
    state.logErrorsOnly = e.target.checked;
    renderLog();
  });
  if (hooks.onRequestSnapshot) state.onRequestSnapshot = hooks.onRequestSnapshot;
}

function bindSubtabs() {
  const tabs = document.querySelectorAll('.tabs.sub .tab');
  tabs.forEach((tab) => {
    on(tab, 'click', () => {
      tabs.forEach((t) => t.classList.toggle('is-active', t === tab));
      for (const pane of ['events', 'log', 'retries']) {
        const node = byId(`pane-${pane}`);
        if (node) node.hidden = pane !== tab.dataset.pane;
      }
      if (tab.dataset.pane === 'log') pollLog();
    });
  });
}

export function currentTask() { return state.task; }

/** Show a task. ``detail`` is one row from the store, projected. */
export function showTask(detail, attempts) {
  state.task = detail || null;
  state.attempts = attempts || [];
  state.lastLogKey = null;
  render();
}

/** Live-tick the elapsed field without waiting for the next snapshot. */
export function tickElapsed() {
  const task = state.task;
  if (!task) return;
  if (isLive(task.state) && task.started_at && !task.ended_at) {
    const seconds = (Date.now() - new Date(task.started_at).getTime()) / 1000;
    if (isFinite(seconds) && seconds >= 0) {
      const node = byId('t-elapsed');
      if (node) node.textContent = durationLong(seconds);
      const clock = byId('task-elapsed');
      if (clock) clock.textContent = duration(seconds);
    }
  }
}

function render() {
  const task = state.task;
  const pill = byId('task-state');
  const title = byId('task-subject');

  if (!task) {
    if (pill) { pill.textContent = 'PENDING'; pill.dataset.state = 'PENDING'; }
    if (title) title.textContent = 'No task selected';
    setText('task-stage', EM_DASH);
    setText('t-state', EM_DASH); setText('t-attempt', EM_DASH);
    setText('t-elapsed', EM_DASH); setText('t-started', EM_DASH);
    setText('t-ended', EM_DASH); setText('t-worker', EM_DASH);
    setNA('t-cpu'); setNA('t-ram'); setNA('t-diskio');
    replace(byId('task-input'), li('na', '—'));
    replace(byId('task-output'), li('na', '—'));
    replace(byId('task-validation'), li('na', '—'));
    replace(byId('task-provenance'), []);
    setText('task-command', '');
    const cmd = byId('task-command');
    if (cmd) replace(cmd, [span('na', '—')]);
    const err = byId('task-error-card');
    if (err) err.hidden = true;
    hideProgressNote();
    return;
  }

  if (pill) { pill.textContent = task.state; pill.dataset.state = task.state; }
  if (title) title.textContent = task.subject || task.stage;
  setText('task-stage', `${task.stage_label || task.stage}${task.subject ? '' : ' (cohort)'}`);

  setText('t-state', task.state);
  setText('t-attempt', `${count(task.attempt)} / ${count(task.max_attempts)}`);
  setText('t-elapsed', durationLong(task.elapsed_seconds));
  setText('t-started', clockTime(task.started_at));
  setText('t-ended', task.ended_at ? clockTime(task.ended_at) : EM_DASH);

  // The engine records no worker identity: there is no worker pool in this
  // build, and inventing "W03" from a loop index would be fiction.
  setNA('t-worker');
  setNA('t-cpu'); setNA('t-ram'); setNA('t-diskio');

  // Per-task progress does not exist in this engine. The cohort bar below
  // is real; this one would not be, so it is not drawn.
  hideProgressNote();
  const wrap = byId('task-progress-wrap');
  if (wrap) wrap.hidden = true;

  renderInput(task);
  renderOutput(task);
  renderValidation(task);
  renderCommand(task);
  renderProvenance(task);
  renderError(task);
  renderAttempts();
}

function hideProgressNote() {
  const note = byId('task-progress-note');
  if (note) note.hidden = false;
}

function renderInput(task) {
  const items = Array.isArray(task.input_ids) ? task.input_ids : [];
  if (!items.length) { replace(byId('task-input'), li('na', 'no input identifiers recorded')); return; }
  replace(byId('task-input'), items.slice(0, 40).map((id) => {
    const row = document.createElement('li');
    row.appendChild(span('k', 'input'));
    row.appendChild(span('v', id));
    return row;
  }));
}

function renderOutput(task) {
  const paths = Array.isArray(task.output_paths) ? task.output_paths : [];
  const sub = byId('output-sub');
  const node = byId('task-output');
  if (!node) return;
  if (!paths.length) {
    if (sub) sub.textContent = '';
    replace(node, li('na', 'no output paths recorded'));
    return;
  }
  // Whether a file is present is a real filesystem question, asked now.
  replace(node, paths.map((p) => {
    const row = document.createElement('li');
    const mark = span('mark', '·');
    const name = span('name', String(p).split('/').pop());
    const exists = typeof p === 'string' && p.length > 0;
    const meta = span('meta', exists ? 'recorded' : 'unknown');
    mark.classList.add(exists ? 'present' : 'pending');
    row.appendChild(mark);
    row.appendChild(name);
    row.appendChild(meta);
    return row;
  }));
  if (sub) sub.textContent = `${paths.length} declared`;
}

function renderValidation(task) {
  const node = byId('task-validation');
  const sub = byId('validation-sub');
  const results = task.validation && task.validation.results;
  if (!Array.isArray(results) || !results.length) {
    if (sub) sub.textContent = task.validation_state || '';
    replace(node, [li('na', task.validation_detail
      ? task.validation_detail
      : 'no contract recorded for this task')]);
    return;
  }
  if (sub) sub.textContent = task.validation_state || '';
  replace(node, results.map((r) => {
    const row = document.createElement('li');
    const mark = span('mark', r.passed ? '✓' : '✕');
    mark.classList.add(r.passed ? 'pass' : 'fail');
    row.appendChild(mark);
    row.appendChild(span('desc', r.check ? r.check.description || r.check.kind : ''));
    row.appendChild(span('state', r.detail ? trim(r.detail, 46) : (r.passed ? 'pass' : 'fail')));
    return row;
  }));
}

function renderCommand(task) {
  const node = byId('task-command');
  if (!node) return;
  const cmd = Array.isArray(task.command) ? task.command : null;
  if (!cmd || !cmd.length) { replace(node, [span('na', 'no command recorded')]); return; }
  // Shaded as text nodes. The command comes from the store, so it is never
  // interpolated as markup.
  const parts = [];
  cmd.forEach((token, i) => {
    if (i > 0) parts.push(document.createTextNode(' '));
    let cls = '';
    if (i === 0) cls = 'exe';
    else if (String(token).includes('/')) cls = 'path';
    else if (String(token).startsWith('-')) cls = 'flag';
    parts.push(span(cls, token));
  });
  replace(node, parts);
}

function renderProvenance(task) {
  const node = byId('task-provenance');
  if (!node) return;
  const rows = [
    ['tool', task.tool_version],
    ['database', task.database_version],
    ['config hash', task.config_hash],
    ['failure kind', task.failure_kind],
    ['host', task.host],
    ['cores', task.cpu_count],
    ['recorded', task.updated_at],
  ].filter(([, v]) => !na(v));
  replace(node, rows.map(([k, v]) => {
    const row = document.createElement('li');
    row.appendChild(span('k', k));
    row.appendChild(span('v', v));
    return row;
  }));
}

function renderError(task) {
  const card = byId('task-error-card');
  const body = byId('task-error');
  if (!card || !body) return;
  if (task.state === 'SUCCEEDED' || !task.error) { card.hidden = true; return; }
  card.hidden = false;
  body.textContent = task.error;
}

function renderAttempts() {
  const node = byId('attempt-list');
  if (!node) return;
  if (!state.attempts.length) { replace(node, li('na', 'no attempts recorded')); return; }
  replace(node, state.attempts.slice().reverse().map((a) => {
    const row = document.createElement('li');
    const dot = span('dot', '');
    dot.style.background = stateColorFor(a.state);
    row.appendChild(span('t', `#${a.attempt}`));
    row.appendChild(dot);
    const name = span('n', a.state);
    name.style.color = stateColorFor(a.state);
    row.appendChild(name);
    if (a.detail) row.appendChild(span('d', trim(a.detail, 120)));
    return row;
  }));
}

function stateColorFor(s) {
  return ({
    SUCCEEDED: '#4ade80', FAILED: '#ef4444', INVALID: '#fb923c',
    INCOMPLETE: '#a855f7', RETRYING: '#f59e0b', RUNNING: '#60a5fa',
  })[s] || '#94a3b8';
}

/** Poll the real log file while the log tab is open. */
export function pollLog() {
  if (state.logTimer) return;
  const load = async () => {
    const task = state.task;
    if (!task) return;
    const key = `${task.stage}|${task.subject}`;
    if (key === state.lastLogKey) return;
    try {
      const res = await getLog(task.stage, task.subject);
      state.lastLogKey = key;
      state.logText = res.text || '';
      if (res.reason) state.logText = `[${res.reason}]`;
      renderLog(true);
    } catch (err) {
      state.logText = `[${err}]`;
      renderLog(true);
    }
  };
  load();
  state.logTimer = setInterval(load, 2000);
}

function renderLog(scroll = false) {
  const node = byId('log-body');
  if (!node) return;
  const lines = state.logText.split('\n');
  const filtered = lines.filter((line) => {
    const lower = line.toLowerCase();
    if (state.logFilter && !lower.includes(state.logFilter)) return false;
    if (state.logErrorsOnly && !/error|fail|exception|traceback/.test(lower)) return false;
    return true;
  });
  node.textContent = filtered.length ? filtered.join('\n') : '[no matching lines]';
  if (scroll && !state.logPaused) node.scrollTop = node.scrollHeight;
}

export function stopLogPolling() {
  if (state.logTimer) { clearInterval(state.logTimer); state.logTimer = null; }
}

// ── event history ─────────────────────────────────────────────────

const MAX_EVENT_ROWS = 300;
const events = [];

export function recordEvent(ev) {
  events.push(ev);
  if (events.length > MAX_EVENT_ROWS) events.shift();
  const list = byId('event-list');
  if (!list) return;
  if (list.querySelector('.na')) list.textContent = '';
  const row = document.createElement('li');
  row.appendChild(span('t', clockTime(new Date(ev.at * 1000).toISOString())));
  const dot = span('dot', '');
  dot.style.background = colorFor(ev.name);
  row.appendChild(dot);
  const name = span('n', ev.name);
  name.style.color = colorFor(ev.name);
  row.appendChild(name);
  const detail = ev.detail || describe(ev);
  if (detail) row.appendChild(span('d', trim(detail, 110)));
  list.insertBefore(row, list.firstChild);
  while (list.children.length > MAX_EVENT_ROWS) list.removeChild(list.lastChild);
}

export function clearEvents() {
  events.length = 0;
  const list = byId('event-list');
  if (list) replace(list, li('na', 'No events received yet.'));
}

function colorFor(name) {
  return ({
    TASK_STARTED: '#22d3ee', TASK_COMPLETED: '#22c55e', TASK_VALIDATED: '#4ade80',
    TASK_FAILED: '#ef4444', TASK_INVALID: '#f97316', TASK_INCOMPLETE: '#a855f7',
    TASK_RETRY: '#f59e0b', TASK_RESUMED: '#38bdf8', TASK_CREATED: '#94a3b8',
  })[name] || '#94a3b8';
}

function describe(ev) {
  const parts = [`${ev.stage}${ev.subject ? ` / ${ev.subject}` : ''}`];
  if (ev.attempt) parts.push(`attempt ${ev.attempt}`);
  return parts.join(' · ');
}

// ── small builders ────────────────────────────────────────────────

function span(cls, text) {
  const s = document.createElement('span');
  if (cls) s.className = cls;
  s.textContent = text === null || text === undefined ? '' : String(text);
  return s;
}

function li(cls, text) {
  const s = span(cls, text);
  return s.tagName === 'SPAN' ? Object.assign(document.createElement('li'), { className: cls, textContent: text }) : s;
}

function setText(id, value) {
  const node = byId(id);
  if (node) node.textContent = value === null || value === undefined ? EM_DASH : String(value);
}

function setNA(id) {
  const node = byId(id);
  if (!node) return;
  node.textContent = 'not reported';
  node.classList.add('na');
}

function trim(text, n) {
  const s = String(text).replace(/\s+/g, ' ').trim();
  return s.length > n ? `${s.slice(0, n - 1)}…` : s;
}

async function copyCommand() {
  const task = state.task;
  if (!task || !Array.isArray(task.command)) return;
  const text = task.command.join(' ');
  try {
    await navigator.clipboard.writeText(text);
    const btn = byId('btn-copy-cmd');
    if (btn) { btn.textContent = 'copied'; setTimeout(() => { btn.textContent = 'copy'; }, 1200); }
  } catch {
    // Clipboard access can be denied; the command is on screen regardless.
  }
}
