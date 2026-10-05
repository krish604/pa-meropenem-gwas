/* app.js — the shell.
 *
 * Builds the DESIGN §7 context and hands it to the router:
 *
 *     ctx = { api, router, theme, banner, badges, on }
 *
 * Exactly those six keys, and the object is frozen, so a page cannot replace
 * the badge vocabulary or the banner with its own and thereby disagree with
 * another page (UI-D2, UI-D3).
 *
 * Everything a page might want that is *not* one of the six lives in the
 * named exports below, so `ctx` stays the contract and the helpers stay
 * library functions.
 *
 * Boot is guarded by a flag `index.html` sets immediately before loading this
 * file. A page module may `import` from here for the helpers without starting
 * a second shell.
 */

/* ------------------------------------------------------------------ DOM */

/**
 * `h(tag, props, ...children)`.
 *
 * There is deliberately no `html` prop and no `innerHTML` anywhere in this
 * file: report text, probe messages and log lines are data, and data does not
 * become DOM (DESIGN §11, "a report containing <script>").
 */
export function h(tag, props, ...children) {
  const el = document.createElement(tag);
  applyProps(el, props);
  append(el, children);
  return el;
}

const SVG_NS = 'http://www.w3.org/2000/svg';

export function svg(tag, props, ...children) {
  const el = document.createElementNS(SVG_NS, tag);
  if (props) {
    for (const [key, value] of Object.entries(props)) {
      if (value === null || value === undefined || value === false) continue;
      if (key === 'class') el.setAttribute('class', String(value));
      else if (key.startsWith('on') && typeof value === 'function') {
        el.addEventListener(key.slice(2).toLowerCase(), value);
      } else if (key === 'text') el.textContent = String(value);
      else el.setAttribute(key, String(value));
    }
  }
  append(el, children);
  return el;
}

function applyProps(el, props) {
  if (!props) return;
  for (const [key, value] of Object.entries(props)) {
    if (value === null || value === undefined || value === false) continue;
    if (key === 'class' || key === 'className') el.className = String(value);
    else if (key === 'text') el.textContent = String(value);
    else if (key === 'dataset') {
      for (const [dk, dv] of Object.entries(value)) {
        if (dv !== null && dv !== undefined) el.dataset[dk] = String(dv);
      }
    } else if (key === 'style' && typeof value === 'object') {
      for (const [sk, sv] of Object.entries(value)) el.style.setProperty(sk, String(sv));
    } else if (key.startsWith('on') && typeof value === 'function') {
      el.addEventListener(key.slice(2).toLowerCase(), value);
    } else if (key in el && key !== 'list' && typeof value !== 'object') {
      try {
        el[key] = value;
      } catch {
        el.setAttribute(key, String(value));
      }
    } else {
      el.setAttribute(key, value === true ? '' : String(value));
    }
  }
}

function append(el, children) {
  for (const child of children.flat(4)) {
    if (child === null || child === undefined || child === false) continue;
    el.appendChild(child instanceof Node ? child : document.createTextNode(String(child)));
  }
}

export function clear(el) {
  while (el.firstChild) el.removeChild(el.firstChild);
  return el;
}

/* -------------------------------------------------------------- absence */

/** The words UI-D2 allows for an absent fact. Never `0`, `""`, `[]`, blank. */
export const ABSENT = Object.freeze({
  NOT_PRODUCED: 'not produced',
  NOT_RUN: 'not run',
  NOT_ASSESSED: 'not assessed',
  NOT_REPORTED: 'not reported',
  TIMINGS: 'timings not recorded',
});

/**
 * Render a value, or the sentence for its absence.
 *
 * `reason` is appended so the reader learns *why* the value is absent — the
 * whole point of D2. `emptyIsAbsent` is for arrays: an empty array is a
 * finding in this project (`[]` means "recorded nothing"), so callers pass
 * `false` for a list whose emptiness is a result.
 */
export function value(valueOrNull, reason, { emptyIsAbsent = true } = {}) {
  if (valueOrNull === null || valueOrNull === undefined) {
    return h('span', { class: 'absent', title: reason || null }, reason || ABSENT.NOT_REPORTED);
  }
  if (Array.isArray(valueOrNull) && emptyIsAbsent && valueOrNull.length === 0) {
    return h('span', { class: 'absent', title: reason || null }, reason || ABSENT.NOT_PRODUCED);
  }
  return document.createTextNode(String(valueOrNull));
}

export function absentSpan(reason, word = ABSENT.NOT_PRODUCED) {
  return h('span', { class: 'absent' }, reason ? `${word}: ${reason}` : word);
}

export function note(text_, className = 'panel-note') {
  return h('p', { class: className }, text_);
}

export function panel(title, ...children) {
  const box = h('div', { class: 'panel' });
  if (title) box.appendChild(h('h3', null, title));
  append(box, children);
  return box;
}

export function kv(pairs) {
  const dl = h('dl', { class: 'kv' });
  for (const [term, node] of pairs) {
    dl.appendChild(h('dt', null, term));
    const dd = h('dd');
    append(dd, [node]);
    dl.appendChild(dd);
  }
  return dl;
}

export function message(severity, title, ...children) {
  const box = h('div', { class: 'message', dataset: { severity } });
  if (title) box.appendChild(h('h3', null, title));
  append(box, children);
  return box;
}

/** The server's own sentence, verbatim, as text. */
export function serverMessage(text_) {
  return h('pre', { class: 'server-message' }, text_);
}

/** A probe/attempt list, in order, each with the reason it did or did not match. */
export function probeList(probes, { pathKey = 'path' } = {}) {
  const list = h('ol', { class: 'probe-list' });
  for (const probe of probes || []) {
    const item = h('li');
    item.appendChild(h('span', 'probe-n', probe.n === undefined ? '·' : String(probe.n)));
    const body = h('div');
    const verdict = probe.ok ? 'matched' : probe.reason || 'no match';
    body.appendChild(
      h('span', probe.ok ? 'probe-ok' : 'probe-miss', `${probe.probe || ''} — ${verdict}`)
    );
    const where = probe[pathKey];
    if (where) body.appendChild(h('div', 'probe-path', where));
    if (probe.note) body.appendChild(note(probe.note, 'faint'));
    item.appendChild(body);
    list.appendChild(item);
  }
  return list;
}

/** A table that refuses to invent rows. `rows` absent -> the reason. */
export function dataTable(columns, rows, absentReason) {
  if (!Array.isArray(rows)) {
    return absentSpan(absentReason || 'this endpoint returned no rows field');
  }
  if (rows.length === 0) {
    return absentSpan(
      absentReason || 'the response carries zero rows; an empty page is not a measurement of zero',
      ABSENT.NOT_PRODUCED
    );
  }
  const table = h('table');
  const thead = h('thead');
  const hr = h('tr');
  for (const column of columns) hr.appendChild(h('th', null, column));
  thead.appendChild(hr);
  const tbody = h('tbody');
  for (const row of rows) {
    const tr = h('tr');
    for (const column of columns) {
      const cell = h('td');
      const cellValue = row[column];
      if (cellValue === null || cellValue === undefined || cellValue === '') {
        cell.appendChild(absentSpan('this cell holds no value'));
      } else {
        cell.textContent = String(cellValue);
      }
      tr.appendChild(cell);
    }
    tbody.appendChild(tr);
  }
  table.appendChild(thead);
  table.appendChild(tbody);
  return h('div', { class: 'table-wrap' }, table);
}

/* ------------------------------------------------------------------ api */

export class ApiError extends Error {
  constructor(message, { status, payload, url, cause } = {}) {
    super(message);
    this.name = 'ApiError';
    this.status = status === undefined ? null : status;
    this.payload = payload === undefined ? null : payload;
    this.url = url;
    if (cause) this.cause = cause;
  }

  /** The §2.3 probe list, when this is a "no results source" failure. */
  get probes() {
    return (this.payload && Array.isArray(this.payload.probes) && this.payload.probes) || null;
  }

  get serverMessage() {
    return (this.payload && this.payload.message) || (this.payload && this.payload.error) || null;
  }
}

const TOKEN_HEADER = 'X-PA-Dash-Token';

/**
 * `api.get(path, params, signal)` / `api.post(path, body, signal)`.
 *
 * Every call takes an AbortSignal and every call is abortable: a page that
 * navigates away mid-flight aborts, and the rejection is an AbortError the
 * page must not render as a failure (DESIGN §7).
 *
 * `params` follows the §5.1 schema: `offset`, `limit`, `sort`, `q`,
 * `include_total`, `columns` (csv), `filter.<column>` (repeated, with
 * `__ne` / `__in` / `__null` suffixes), and `include_total=0`.
 */
export function createApi({ base = '/api' } = {}) {
  const token = { value: null };

  function buildQuery(params) {
    if (!params) return '';
    const search = new URLSearchParams();
    for (const [key, raw] of Object.entries(params)) {
      if (raw === null || raw === undefined || raw === '') continue;
      if (Array.isArray(raw)) {
        for (const item of raw) search.append(key, String(item));
      } else {
        search.append(key, String(raw));
      }
    }
    const qs = search.toString();
    return qs ? `?${qs}` : '';
  }

  async function request(path, { method = 'GET', params, body, signal } = {}) {
    const url = `${base}${path}${buildQuery(params)}`;
    const headers = { Accept: 'application/json' };
    if (token.value) headers[TOKEN_HEADER] = token.value;
    if (body !== undefined) headers['Content-Type'] = 'application/json';
    let response;
    try {
      response = await fetch(url, {
        method,
        headers,
        signal,
        credentials: 'same-origin',
        cache: 'no-store',
        body: body === undefined ? undefined : JSON.stringify(body),
      });
    } catch (error) {
      if (error && error.name === 'AbortError') throw error;
      throw new ApiError(
        `${method} ${url} could not be reached: ${error && error.message ? error.message : error}`,
        { url, cause: error }
      );
    }
    const text_ = await response.text();
    let payload = null;
    if (text_) {
      try {
        payload = JSON.parse(text_);
      } catch {
        payload = { error: text_.slice(0, 2000) };
      }
    }
    if (!response.ok) {
      const detail =
        (payload && (payload.error || payload.message || payload.detail)) ||
        `${method} ${url} answered ${response.status}`;
      throw new ApiError(String(detail), { status: response.status, payload, url });
    }
    return payload;
  }

  return {
    get: (path, params, signal) => request(path, { method: 'GET', params, signal }),
    post: (path, body, signal) => request(path, { method: 'POST', body, signal }),
    /**
     * The token, in memory only.
     *
     * It is never a query parameter (UI-D4.3) and never written to storage
     * (UI-D1): a token in localStorage is a token every script on the origin
     * can read, and this origin may serve report text.
     */
    setToken(value) {
      token.value = value || null;
    },
    hasToken: () => Boolean(token.value),
  };
}

/* --------------------------------------------------------------- badges */

/**
 * The six-badge vocabulary of `dashboard/server/badges.py` §4.2.
 *
 * `labels`, `glyphs` and `tones` are copied from that module and are used ONLY
 * when `/api/health` could not be reached — i.e. when there is no server
 * vocabulary to prefer. As soon as `/api/health` answers, the server's own
 * `badges` object replaces this wholesale, so the wording a reader sees comes
 * from the backend and two pages cannot disagree about what `not_assessed`
 * means (§7).
 */
const FALLBACK_VOCABULARY = Object.freeze({
  completed: { label: 'completed', glyph: '●', tone: 'ok' },
  running: { label: 'running', glyph: '◐', tone: 'busy' },
  failed: { label: 'failed', glyph: '✕', tone: 'bad' },
  refused: { label: 'refused', glyph: '⊘', tone: 'warn' },
  not_assessed: { label: 'not assessed', glyph: '○', tone: 'unknown' },
  not_run: { label: 'not run', glyph: '–', tone: 'unknown' },
});

const TONE_ORDER = ['ok', 'busy', 'bad', 'warn', 'unknown'];

/**
 * `ctx.badges` — the only place a state becomes a label, a colour token and a
 * reason. No page assembles a status sentence itself (UI-D2).
 */
export function createBadges(vocabulary) {
  let table = { ...FALLBACK_VOCABULARY };
  let source = 'offline-fallback';

  function adopt(serverBadges) {
    if (!serverBadges || typeof serverBadges !== 'object') return;
    const next = {};
    for (const [key, entry] of Object.entries(serverBadges)) {
      if (!entry || typeof entry !== 'object') continue;
      next[key] = {
        label: entry.label ?? key,
        glyph: entry.glyph ?? '?',
        tone: TONE_ORDER.includes(entry.tone) ? entry.tone : 'unknown',
        reason: typeof entry.reason === 'string' ? entry.reason : '',
      };
    }
    if (Object.keys(next).length > 0) {
      table = next;
      source = 'server';
    }
  }

  function unknown(reason, word = ABSENT.NOT_REPORTED, glyph = '?') {
    return { key: word, label: word, glyph, tone: 'unknown', reason: reason || 'no reason recorded' };
  }

  /**
   * `state(x)` accepts anything a page might have:
   *
   *  - a stage payload (`{state, reason, badge}`) — the server's own
   *    classification wins, because `badges.classify` also read the manifest
   *    and the table;
   *  - a badge object (`{key, label, glyph, tone, reason}`);
   *  - one of the six keys;
   *  - `null`/`undefined` — a missing value, which reads `not reported` with a
   *    reason, never a blank and never `0`.
   *
   * Anything else is *not* coerced into a badge: it is reported as
   * unrecognised, with the token that was not recognised.
   */
  function state(input, hint = '') {
    if (input === null || input === undefined || input === '') {
      return unknown(
        hint || 'no state was supplied by the server for this field, so none is claimed'
      );
    }
    if (typeof input === 'object') {
      if (input.badge && typeof input.badge === 'object' && input.badge.key) {
        return { ...input.badge, reason: input.badge.reason || input.reason || '' };
      }
      if (input.key && input.tone) {
        return { ...input, reason: input.reason || '' };
      }
      if (typeof input.state === 'string') {
        return state(input.state, input.reason || hint);
      }
      return unknown('the server sent an object with neither a badge nor a state key');
    }
    const key = String(input);
    const entry = table[key];
    if (!entry) {
      return unknown(
        `${JSON.stringify(key)} is not one of the six badges; it is shown as ` +
          `unrecognised rather than mapped onto one of them`,
        'state not recognised',
        '?'
      );
    }
    return { key, label: entry.label, glyph: entry.glyph, tone: entry.tone, reason: '' };
  }

  /** `tone(key)` -> the tone token, for CSS. Never a colour value. */
  function tone(key) {
    if (key && typeof key === 'object') return key.tone || 'unknown';
    const entry = table[String(key)];
    return entry ? entry.tone : 'unknown';
  }

  /**
   * Render a badge: colour AND glyph AND word, always together.
   *
   * The reason follows the badge for every badge that is not `completed`,
   * because a non-completed badge with an empty reason is indistinguishable
   * from a bug (`badges._with_reason`).
   */
  function render(input, hint) {
    const badge = state(input, hint);
    const el = h('span', { class: 'badge', dataset: { tone: badge.tone, badge: badge.key } });
    el.appendChild(h('span', { class: 'badge-glyph', 'aria-hidden': 'true' }, badge.glyph));
    el.appendChild(h('span', { class: 'badge-label' }, badge.label));
    el.title = badge.reason || '';
    if (badge.key !== 'completed' && badge.reason) {
      el.appendChild(h('span', { class: 'badge-reason' }, `— ${badge.reason}`));
    }
    return el;
  }

  /** The D2 wording for a table that exists but holds nothing / does not exist. */
  function notProduced(reason) {
    return unknown(reason || 'the server did not say why', ABSENT.NOT_PRODUCED, '∅');
  }

  function notAssessed(reason) {
    return unknown(reason || 'the server did not say why', ABSENT.NOT_ASSESSED, '○');
  }

  function notReported(reason) {
    return unknown(reason || 'the server did not say why', ABSENT.NOT_REPORTED, '∅');
  }

  function timings(reason) {
    return unknown(
      reason || 'no runbook supplied per-isolate timings; an estimate derived from the cohort size would be a measurement nobody made',
      ABSENT.TIMINGS,
      '∅'
    );
  }

  return Object.freeze({
    state,
    tone,
    render,
    adopt,
    notProduced,
    notAssessed,
    notReported,
    timings,
    vocabulary: () => ({ ...table }),
    source: () => source,
  });
}

/* --------------------------------------------------------------- banner */

/**
 * Fallback copies of `papipeline/stages/reporting.py:123-131`.
 *
 * Used only until the server sends a `power` block carrying its own `template`
 * and `meaning` — at which point the server's strings win. The wording is
 * DESIGN §7.1's, verbatim.
 */
const POWER_TEMPLATE_FALLBACK = 'n={n}, underpowered, not a finding';
const POWER_MEANING_FALLBACK =
  'Ruling R16. This is a description of the isolates in this run and nothing more. ' +
  'It is not evidence about Pseudomonas aeruginosa imipenem susceptibility in general, ' +
  'and no association, rate or difference below may be cited as a finding.';

/**
 * `ctx.banner` — the D3 component's parameters, and its renderer.
 *
 * `N` is always `basis.n`: the cohort the statistic was ACTUALLY computed on,
 * never the cohort the page was told about (§7.1). A page that renders a
 * filtered count beside the full-cohort flag is exactly the error R16 exists to
 * prevent, so `power()` takes an `n` and refuses to guess one.
 */
export function createBanner() {
  let template = POWER_TEMPLATE_FALLBACK;
  let meaning = POWER_MEANING_FALLBACK;
  let minSamples = null;
  let mode = null;
  let modeText = null;

  function adopt(payload) {
    if (!payload || typeof payload !== 'object') return;
    if (typeof payload.template === 'string' && payload.template.includes('{n}')) {
      template = payload.template;
    }
    if (typeof payload.meaning === 'string' && payload.meaning.trim()) {
      meaning = payload.meaning;
    }
    if (payload.min_samples !== undefined && payload.min_samples !== null) {
      minSamples = Number(payload.min_samples);
    }
  }

  /** `{n, flag, template, meaning, underpowered, min_samples}`. */
  function power(n, basis) {
    const resolved =
      n === null || n === undefined ? basis && basis.n !== undefined ? basis.n : null : Number(n);
    const floor =
      basis && basis.min_samples !== undefined && basis.min_samples !== null
        ? Number(basis.min_samples)
        : minSamples;
    const underpowered =
      resolved === null || floor === null || floor === undefined ? null : resolved < floor;
    // The pipeline's template reads "n={n}, underpowered, not a finding" for
    // every n. On a cohort at or above the configured minimum that word is
    // false, so it is dropped from the rendered flag; the R16 meaning is
    // unchanged. `underpowered` (from `basis.min_samples`) is the only thing
    // that decides it.
    let flag = resolved === null ? null : template.replace('{n}', String(resolved));
    if (flag !== null && underpowered === false && flag.includes('underpowered')) {
      flag = flag
        .replace(/,\s*underpowered\s*,/, ',')
        .replace(/underpowered\s*,/, '')
        .replace(/\s*underpowered/, '');
    }
    return {
      n: resolved,
      flag,
      template,
      meaning,
      underpowered,
      min_samples: floor === undefined ? null : floor,
    };
  }

  function underpowered(n, basis) {
    return power(n, basis).underpowered;
  }

  /** The two-line flag, as plain text, with the server's markdown emphasis removed. */
  function text(n, basis) {
    const p = power(n, basis);
    if (p.flag === null) {
      return `${template.replace('{n}', '(not reported)')} — the cohort this statistic was computed on is not recorded, so no flag is claimed.`;
    }
    return `${p.flag}\n${meaning.replace(/\*/g, '')}`;
  }

  /**
   * The one component. Mount it wherever a statistic renders; a page cannot
   * opt out because the call site is the statistic, not the page (UI-D3).
   */
  function mount(container, n, basis) {
    const p = power(n, basis);
    const el = h('div', {
      class: 'power-banner',
      dataset: { underpowered: p.underpowered === null ? 'unknown' : String(p.underpowered) },
    });
    el.appendChild(h('strong', { class: 'power-flag' }, p.flag ?? `${template.replace('{n}', '(not reported)')}`));
    el.appendChild(h('span', { class: 'power-meaning' }, meaning.replace(/\*/g, '')));
    if (p.min_samples !== null && p.min_samples !== undefined) {
      el.appendChild(
        h(
          'small',
          { class: 'faint' },
          ` n = ${p.n === null ? ABSENT.NOT_REPORTED : p.n}; the configured minimum is ${p.min_samples} (config/science.yaml).`
        )
      );
    }
    if (p.underpowered === false) {
      el.appendChild(
        h('small', { class: 'faint' }, ' n is at or above the configured minimum; the flag is still shown, as R16 requires.')
      );
    }
    if (p.n === null) {
      el.appendChild(
        h('small', { class: 'faint' }, ' The cohort this statistic was computed on was not reported.')
      );
    }
    if (container) container.appendChild(el);
    return el;
  }

  return Object.freeze({
    get mode() {
      return mode;
    },
    get text_() {
      return modeText;
    },
    get minSamples() {
      return minSamples;
    },
    setRun({ mode: m, banner } = {}) {
      if (typeof m === 'string') mode = m;
      if (typeof banner === 'string') modeText = banner;
    },
    adopt,
    power,
    underpowered,
    text,
    mount,
  });
}

/**
 * The `basis` a statistic was computed on, from whichever of the server's
 * three shapes carries it: `meta.basis` (every paged list), `power.n` (same
 * responses), or `bases.<key>` (`/api/run/counts`). Returns null rather than a
 * guess.
 */
export function basisOf(response, key) {
  if (!response || typeof response !== 'object') return null;
  if (response.meta && response.meta.basis) return response.meta.basis;
  if (key && response.bases && response.bases[key]) return response.bases[key];
  if (response.bases && !key) {
    const first = Object.values(response.bases)[0];
    if (first) return first;
  }
  if (response.power && typeof response.power.n === 'number') {
    return {
      n: response.power.n,
      artefact: null,
      path: null,
      rows_total: response.meta && response.meta.total !== undefined ? response.meta.total : null,
      filtered: Boolean(response.meta && response.meta.filters && Object.keys(response.meta.filters).length),
      min_samples: response.power.min_samples,
    };
  }
  if (response.basis) return response.basis;
  return null;
}

/* ---------------------------------------------------------------- theme */

/**
 * `ctx.theme` — `{tokens}` plus the session-only setter.
 *
 * `tokens` is read from the CSS custom properties on every access and is
 * frozen: it is a reading of the stylesheet, not a second palette that can
 * drift from it (§7, §13).
 *
 * Nothing is persisted. UI-D1 makes the state directory the only write
 * location, and it is not reachable from the browser, so the theme follows the
 * OS until the reader chooses otherwise for this session.
 */
export function createTheme({ root = document.documentElement } = {}) {
  const NAMES = [
    'bg', 'bg-raised', 'bg-sunken', 'bg-inset', 'border', 'border-strong',
    'fg', 'fg-muted', 'fg-faint', 'accent', 'accent-fg', 'accent-soft', 'link', 'focus',
    'tone-ok-fg', 'tone-ok-bg', 'tone-ok-line',
    'tone-busy-fg', 'tone-busy-bg', 'tone-busy-line',
    'tone-bad-fg', 'tone-bad-bg', 'tone-bad-line',
    'tone-warn-fg', 'tone-warn-bg', 'tone-warn-line',
    'tone-unknown-fg', 'tone-unknown-bg', 'tone-unknown-line',
    'power-fg', 'power-bg', 'power-line', 'power-glyph',
    'mode-fg', 'mode-bg', 'mode-line',
  ];

  function readTokens() {
    const style = getComputedStyle(root);
    const tokens = {};
    for (const name of NAMES) {
      tokens[name] = style.getPropertyValue(`--${name}`).trim();
    }
    return Object.freeze(tokens);
  }

  let cached = null;

  function tokens() {
    // Re-read on demand: a stylesheet edit or a scheme change must be visible
    // without a reload, and a cache would make the value a second palette.
    return readTokens();
  }

  function current() {
    const explicit = root.getAttribute('data-theme');
    if (explicit === 'light' || explicit === 'dark') return explicit;
    const prefersDark =
      typeof window.matchMedia === 'function' &&
      window.matchMedia('(prefers-color-scheme: dark)').matches;
    return prefersDark ? 'dark' : 'light';
  }

  function set(scheme) {
    if (scheme === 'light' || scheme === 'dark') root.setAttribute('data-theme', scheme);
    else root.removeAttribute('data-theme');
    cached = null;
    return current();
  }

  /** system -> light -> dark -> system */
  function cycle() {
    const now = root.getAttribute('data-theme');
    if (!now) return set('light');
    if (now === 'light') return set('dark');
    return set(null);
  }

  function describe() {
    const explicit = root.getAttribute('data-theme');
    if (explicit === 'light' || explicit === 'dark') return `theme: ${explicit}`;
    return 'theme: follow system';
  }

  return Object.freeze({
    get tokens() {
      if (!cached) cached = tokens();
      return cached;
    },
    read: tokens,
    current,
    set,
    cycle,
    describe,
    invalidate: () => {
      cached = null;
    },
  });
}

/* ------------------------------------------------------------ transport */

const SSE_FRAME_TIMEOUT_MS = 3000;
const POLL_INTERVAL_MS = 2000;

/**
 * The live transport, with the polling fallback, exactly as §6.3 describes it.
 *
 * Three outcomes, all explicit and all visible in the header:
 *
 *  1. `open` + any frame -> `sse`. The 3 s timer is stopped.
 *  2. timeout, then `/api/events/state` succeeds -> `poll`, every 2 s,
 *     carrying `cursor`.
 *  3. timeout AND `/api/events/state` fails -> `none`, with the reason. The
 *     run summary is still served; the event panel says the stream is
 *     unavailable rather than showing an empty log, because an empty log reads
 *     as "nothing has happened".
 */
export function createTransport({ api, onState, onDelta, onTransport } = {}) {
  const stateSubs = new Set();
  const deltaSubs = new Set();

  let source = null;
  let timer = null;
  let pollTimer = null;
  let abort = null;
  let cursor = 0;
  let transport = 'connecting';
  let reason = '';
  let closed = false;
  let lastSnapshot = null;

  function announce(next, why) {
    transport = next;
    reason = why || '';
    if (onTransport) {
      try {
        onTransport({ transport, reason, cursor });
      } catch (error) {
        console.error('[transport] the transport indicator threw', error);
      }
    }
  }

  function emitState(snapshot) {
    lastSnapshot = snapshot;
    for (const fn of Array.from(stateSubs)) {
      try {
        fn(snapshot);
      } catch (error) {
        console.error('[transport] a runState subscriber threw', error);
      }
    }
  }

  function emitDelta(delta, meta) {
    if (!delta || delta.length === 0) return;
    for (const fn of Array.from(deltaSubs)) {
      try {
        fn(delta, meta);
      } catch (error) {
        console.error('[transport] a delta subscriber threw', error);
      }
    }
  }

  function stopTimer() {
    if (timer !== null) {
      clearTimeout(timer);
      timer = null;
    }
  }

  function stopPolling() {
    if (pollTimer !== null) {
      clearInterval(pollTimer);
      pollTimer = null;
    }
    if (abort) {
      abort.abort();
      abort = null;
    }
  }

  function closeSource() {
    if (source) {
      source.close();
      source = null;
    }
  }

  function gotFrame() {
    stopTimer();
    if (transport !== 'sse') announce('sse', 'the event stream answered with a frame');
  }

  async function pollOnce() {
    stopPolling();
    abort = new AbortController();
    try {
      const payload = await api.get('/events/state', { cursor }, abort.signal);
      transport = 'poll';
      announce('poll', 'polling /api/events/state every 2 s with the cursor');
      if (typeof payload.cursor === 'number') cursor = payload.cursor;
      emitState(payload);
      emitDelta(Array.isArray(payload.delta) ? payload.delta : [], { source: 'poll' });
      if (closed) return;
      pollTimer = setInterval(pollOnce, POLL_INTERVAL_MS);
    } catch (error) {
      if (error && error.name === 'AbortError') return;
      stopPolling();
      announce('none', error && error.message ? error.message : String(error));
      emitState({
        transport: 'none',
        reason: error && error.message ? error.message : String(error),
        delta: [],
        stages: [],
        totals: null,
      });
    }
  }

  function start() {
    closed = false;
    announce('connecting', 'opening /api/events and arming a 3 s fallback timer');
    closeSource();
    stopPolling();
    cursor = 0;

    if (typeof EventSource === 'undefined') {
      stopTimer();
      void pollOnce();
      return;
    }

    try {
      source = new EventSource('/api/events');
    } catch (error) {
      stopTimer();
      announce('none', `EventSource could not be constructed: ${error && error.message}`);
      void pollOnce();
      return;
    }

    source.addEventListener('snapshot', (event) => {
      gotFrame();
      const snapshot = safeParse(event.data);
      if (!snapshot) return;
      if (typeof snapshot.cursor === 'number') cursor = snapshot.cursor;
      emitState(snapshot);
    });

    source.addEventListener('delta', (event) => {
      gotFrame();
      const payload = safeParse(event.data);
      if (!payload) return;
      if (Array.isArray(payload)) {
        emitDelta(payload, { source: 'sse' });
        return;
      }
      if (Array.isArray(payload.delta)) emitDelta(payload.delta, { source: 'sse' });
      if (typeof payload.cursor === 'number') cursor = payload.cursor;
      // A delta after a reconnect carries no full state, so ask for one: the
      // stream's snapshot is the only thing that carries every stage badge.
      void refreshSnapshot();
    });

    source.onopen = () => {
      /* `open` alone is not enough (§6.3 outcome 1 needs a frame). */
    };

    source.onerror = () => {
      if (closed) return;
      if (source && source.readyState === 2 /* CLOSED */) {
        stopTimer();
        closeSource();
        void pollOnce();
      }
    };

    // Armed simultaneously with the EventSource, per §6.3.
    timer = setTimeout(() => {
      timer = null;
      closeSource();
      void pollOnce();
    }, SSE_FRAME_TIMEOUT_MS);
  }

  async function refreshSnapshot() {
    const controller = new AbortController();
    try {
      const payload = await api.get('/events/state', { cursor: 0 }, controller.signal);
      emitState(payload);
    } catch (error) {
      if (error && error.name === 'AbortError') return;
      if (onState) onState(null, error);
    }
  }

  function safeParse(data) {
    try {
      return JSON.parse(data);
    } catch (error) {
      console.error('[transport] a frame was not JSON', error);
      return null;
    }
  }

  function close() {
    closed = true;
    stopTimer();
    stopPolling();
    closeSource();
    stateSubs.clear();
    deltaSubs.clear();
    transport = 'idle';
  }

  return {
    start,
    close,
    state: () => ({ transport, reason, cursor, snapshot: lastSnapshot }),
    subscribeState(fn) {
      stateSubs.add(fn);
      return () => stateSubs.delete(fn);
    },
    subscribeDelta(fn) {
      deltaSubs.add(fn);
      return () => deltaSubs.delete(fn);
    },
  };
}

/* --------------------------------------------------------------- routes */

/**
 * The route table.
 *
 * A page module another owner is writing RIGHT NOW is declared here but loaded
 * lazily, so a missing script is a named absence in one panel rather than a
 * blank page (D2). `owner` is only used for the nav's honest state.
 */
export const ROUTES = [
  { id: 'source', path: '/source', label: 'Results source', group: 'run', owner: 'SHELL', specifier: 'app.js#mountSourcePicker', load: () => Promise.resolve(import('./app.js')), exportName: 'mountSourcePicker' },
  { id: 'pipeline', path: '/pipeline', label: 'Pipeline', group: 'run', owner: 'SHELL', specifier: 'pages/pipeline.js', load: () => import('./pages/pipeline.js'), exportName: 'mount' },
  { id: 'monitor', path: '/monitor', label: 'Run monitor', group: 'run', owner: 'SHELL', specifier: 'pages/pipeline.js', load: () => import('./pages/pipeline.js'), exportName: 'mountMonitor' },
  { id: 'provenance', path: '/provenance', label: 'Provenance', group: 'run', owner: 'SHELL', specifier: 'pages/provenance.js', load: () => import('./pages/provenance.js'), exportName: 'mount' },
  { id: 'isolates', path: '/isolates', label: 'Isolates', group: 'results', owner: 'DATA', specifier: 'pages/isolates.js', load: () => import('./pages/isolates.js'), exportName: 'mount' },
  { id: 'isolate', path: '/isolates/:sample_id', label: 'Isolate detail', group: 'results', owner: 'DATA', hidden: true, specifier: 'pages/isolate.js', load: () => import('./pages/isolate.js'), exportName: 'mount' },
  { id: 'tables', path: '/tables', label: 'Tables', group: 'results', owner: 'DATA', specifier: 'pages/tables.js', load: () => import('./pages/tables.js'), exportName: 'mount' },
  { id: 'summaries', path: '/summaries', label: 'Summaries', group: 'results', owner: 'DATA', specifier: 'pages/summaries.js', load: () => import('./pages/summaries.js'), exportName: 'mount' },
  { id: 'oprd', path: '/oprd', label: 'oprD', group: 'results', owner: 'DATA', specifier: 'pages/oprd.js', load: () => import('./pages/oprd.js'), exportName: 'mount' },
  { id: 'tree', path: '/tree', label: 'Tree', group: 'analysis', owner: 'VIZ', specifier: 'pages/tree.js', load: () => import('./pages/tree.js'), exportName: 'mount' },
  { id: 'matrix', path: '/matrix', label: 'Similarity', group: 'analysis', owner: 'VIZ', specifier: 'pages/matrix.js', load: () => import('./pages/matrix.js'), exportName: 'mount' },
  { id: 'gwas', path: '/gwas', label: 'GWAS', group: 'analysis', owner: 'VIZ', specifier: 'pages/gwas.js', load: () => import('./pages/gwas.js'), exportName: 'mount' },
  { id: 'evolution', path: '/evolution', label: 'Convergence &amp; combination', group: 'analysis', owner: 'VIZ', specifier: 'pages/evolution.js', load: () => import('./pages/evolution.js'), exportName: 'mount' },
];

export const GROUPS = [
  { id: 'run', label: 'run' },
  { id: 'results', label: 'results' },
  { id: 'analysis', label: 'analysis' },
];

/* -------------------------------------------------- the source picker */

/**
 * S5 — the results-source panel.
 *
 * Three honest facts, in this order:
 *
 *  1. WHICH layout was detected — `live` or `bundle`, from
 *     `/api/health`'s `source_kind`, with the root and the manifest writer.
 *  2. The probe search, in order, each with the reason it did or did not
 *     match, whenever the search failed (§2.3).
 *  3. The recent roots — and, when the backend does not carry them, that they
 *     are **not persisted**, because the browser writes nothing locally
 *     (UI-D1: the state directory is not reachable from here).
 *
 * The path input does not, and cannot, re-point the running server: the
 * results root is read once at start-up and there is no endpoint that changes
 * it. The panel therefore shows the command a human runs instead of
 * pretending the field is wired.
 */
export function mountSourcePicker(container, ctx) {
  const abort = new AbortController();
  const state = { health: null, error: null, recents: null, recentsReason: null };
  const disposers = [];
  let disposed = false;

  const head = h('div', { class: 'page-head' });
  head.appendChild(
    h(
      'div',
      { class: 'page-head-text' },
      h('h2', null, 'Results source'),
      h(
        'p',
        { class: 'page-note' },
        'The dashboard opened one results root at start-up and reads every file ' +
          'under it read-only. This panel reports which layout was detected and ' +
          'names every path the search probed, in order.'
      )
    )
  );
  container.appendChild(head);

  const detected = h('div', { class: 'panel' });
  const search = h('div', { class: 'panel' });
  const recent = h('div', { class: 'panel' });
  container.append(detected, search, recent);

  function renderDetected() {
    clear(detected);
    detected.appendChild(h('h3', null, 'Detected layout'));
    const health = state.health;
    if (!health) {
      detected.appendChild(
        absentSpan(
          state.error
            ? `${state.error}`
            : 'the server has not reported which layout it opened',
          ABSENT.NOT_REPORTED
        )
      );
      return;
    }
    detected.appendChild(
      kv([
        ['layout', h('strong', null, health.source_kind === 'bundle' ? 'delivery bundle' : 'live results root')],
        ['root', h('code', null, health.root || ABSENT.NOT_REPORTED)],
        ['manifest writer', h('code', null, health.manifest_writer || 'none')],
        ['run mode', health.mode ? h('code', null, health.mode) : h('span', { class: 'absent' }, 'no manifest records a run mode')],
        ['pipeline version', health.pipeline_version ? h('code', null, health.pipeline_version) : h('span', { class: 'absent' }, ABSENT.NOT_REPORTED)],
        ['state directory', h('code', null, health.state_dir || ABSENT.NOT_REPORTED)],
        ['bind', h('code', null, health.loopback ? 'loopback (no token required)' : 'non-loopback (a token is required on every /api request)')],
      ])
    );
    detected.appendChild(
      note(
        'A bundle is a directory carrying 02_stage_outputs/ and the rest of the ' +
          'delivery layout; a live root is results/<mode>/ from the machine ' +
          "overlay. The two are told apart by where run_manifest.json landed, not " +
          'by which probe matched.'
      )
    );
  }

  function renderSearch() {
    clear(search);
    search.appendChild(h('h3', null, 'Probe search'));
    const probes = state.error ? state.error.probes : state.health && state.health.probes;
    if (!probes || probes.length === 0) {
      search.appendChild(
        note(
          'No probe list was returned: the server opened a root on the first ' +
            'probe, or the endpoint that carries the search was not reached.',
          'faint'
        )
      );
      return;
    }
    search.appendChild(
      note(
        `The ${probes.length} probes below are the ones the server tried, in ` +
          'order. Probes 1–3 report an absence of configuration; probes 4–6 ' +
          'report a directory that exists and holds nothing. Those are ' +
          'different facts and are not collapsed into "not found".'
      )
    );
    search.appendChild(probeList(probes));
    if (state.error && state.error.serverMessage) {
      search.appendChild(h('h4', null, "The server's message, verbatim"));
      search.appendChild(serverMessage(state.error.serverMessage));
    }
  }

  function renderRecent() {
    clear(recent);
    recent.appendChild(h('h3', null, 'Recent roots'));
    const pathInput = h('input', {
      type: 'text',
      id: 'results-root-input',
      class: 'phrase-input',
      placeholder: '/absolute/path/to/results',
      spellcheck: false,
      autocapitalize: 'off',
      autocomplete: 'off',
    });
    const status = h('p', { class: 'status-line' });
    const show = h('button', { type: 'button', class: 'primary' }, 'How do I open this path?');

    const onShow = () => {
      const typed = pathInput.value;
      clear(status);
      if (!typed.trim()) {
        status.appendChild(absentSpan('no path was typed', 'nothing to show'));
        return;
      }
      if (!typed.startsWith('/')) {
        status.appendChild(
          absentSpan(`"${typed}" is not an absolute path; the server resolves relative to its own working directory`, 'not an absolute path')
        );
        return;
      }
      status.appendChild(
        h(
          'span',
          null,
          'Point the server at it by restarting it: ',
          h('code', null, `python -m dashboard --results ${typed}`),
          ' or ',
          h('code', null, `PA_DASH_RESULTS_ROOT=${typed}`),
          '. Nothing in this page changed the running server: the results root is read once at start-up, and no endpoint changes it.'
        )
      );
    };
    const onKey = (event) => {
      if (event.key === 'Enter') {
        event.preventDefault();
        onShow();
      }
    };
    show.addEventListener('click', onShow);
    pathInput.addEventListener('keydown', onKey);
    disposers.push(() => {
      show.removeEventListener('click', onShow);
      pathInput.removeEventListener('keydown', onKey);
    });

    recent.appendChild(
      kv([
        ['path', h('div', null, pathInput)],
        ['action', show],
      ])
    );
    recent.appendChild(status);
    recent.appendChild(h('h4', null, 'Previously opened'));
    if (Array.isArray(state.recents) && state.recents.length > 0) {
      recent.appendChild(dataTable(['root', 'opened_at_utc', 'kind'], state.recents, 'no recends row'));
    } else {
      recent.appendChild(absentSpan(state.recentsReason, ABSENT.NOT_REPORTED));
      recent.appendChild(
        note(
          'The browser writes nothing locally — no cookies, no localStorage, no ' +
            'sessionStorage — because UI-D1 makes the state directory the only ' +
            'write location and it is not reachable from here. A recent-roots ' +
            'list therefore needs a backend endpoint to read and write ' +
            '<state>/session.json; none exists, so the list is reported absent ' +
            'rather than kept in this tab.'
        )
      );
    }
  }

  function renderAll() {
    renderDetected();
    renderSearch();
    renderRecent();
  }

  (async () => {
    try {
      state.health = await ctx.api.get('/health', null, abort.signal);
      state.error = null;
    } catch (error) {
      if (error && error.name === 'AbortError') return;
      state.health = null;
      state.error = error;
    }
    // `recents` would be carried on /api/health if the backend exposed it. It
    // does not today; the panel says so rather than probing an endpoint that
    // does not exist (which would be a 404 in every screenshot).
    state.recents = state.health && Array.isArray(state.health.recents) ? state.health.recents : null;
    state.recentsReason =
      'the backend carries no recents field: /api/health answers without one, ' +
      'and the only write location (the state directory) is not reachable from ' +
      'the browser';
    if (!disposed) renderAll();
  })();

  return () => {
    disposed = true;
    abort.abort();
    for (const fn of disposers) {
      try {
        fn();
      } catch (error) {
        console.error('[source-picker] a disposer threw', error);
      }
    }
  };
}

/* ----------------------------------------------------------------- boot */

function renderTransport(el, { transport, reason }) {
  el.dataset.transport = transport;
  const glyphs = { sse: '▶', poll: '⟳', none: '✕', connecting: '…', idle: '·' };
  const labels = {
    sse: 'transport: SSE (server-sent events)',
    poll: 'transport: polling every 2 s with the cursor',
    none: 'transport: unavailable',
    connecting: 'transport: connecting',
    idle: 'transport: idle',
  };
  el.querySelector('.transport-glyph').textContent = glyphs[transport] || '·';
  el.querySelector('.transport-label').textContent = labels[transport] || `transport: ${transport}`;
  el.querySelector('.transport-note').textContent = reason ? `— ${reason}` : '';
  el.title =
    'Which mechanism is feeding this page (§6.3). The client opens the stream ' +
    'and arms a 3 s timer; a frame means SSE, no frame plus a working ' +
    '/api/events/state means polling, and neither means this label.';
}

function fact(term, node) {
  return [term, node];
}

function renderFacts(el, health, run, minSamples, healthError) {
  clear(el);
  const pairs = [];
  if (health) {
    pairs.push(fact('root', health.root ? h('code', null, health.root) : h('span', { class: 'fact-unknown' }, 'no results source open')));
    pairs.push(fact('layout', h('code', null, health.source_kind || 'none')));
    pairs.push(fact('manifest', h('code', null, health.manifest_writer || 'none')));
  } else {
    pairs.push(
      fact(
        'root',
        h(
          'span',
          { class: 'fact-unknown' },
          healthError ? `no results source open — ${healthError.message}` : 'the server did not answer /api/health'
        )
      )
    );
  }
  if (run) {
    pairs.push(fact('mode', run.run_mode ? h('code', null, run.run_mode) : h('span', { class: 'fact-unknown' }, ABSENT.NOT_REPORTED)));
    pairs.push(
      fact(
        'antibiotic',
        run.antibiotic ? h('code', null, run.antibiotic) : h('span', { class: 'fact-unknown' }, ABSENT.NOT_REPORTED)
      )
    );
    pairs.push(
      fact(
        'cohort (n)',
        typeof run.n_samples === 'number'
          ? h('b', null, String(run.n_samples))
          : h('span', { class: 'fact-unknown' }, `${ABSENT.NOT_REPORTED} — no manifest key records it`)
      )
    );
    pairs.push(
      fact(
        'stages completed',
        Array.isArray(run.stages)
          ? h('b', null, `${run.stages.filter((s) => s.state === 'completed').length} of ${run.stages.length}`)
          : h('span', { class: 'fact-unknown' }, ABSENT.NOT_REPORTED)
      )
    );
  }
  if (minSamples !== null && minSamples !== undefined) {
    pairs.push(fact('min n (R16)', h('b', null, String(minSamples))));
  }
  for (const [term, node] of pairs) {
    el.appendChild(h('dt', null, term));
    const dd = h('dd');
    dd.appendChild(node);
    el.appendChild(dd);
  }
}

function renderNav(list, availability) {
  clear(list);
  for (const group of GROUPS) {
    list.appendChild(h('li', { class: 'nav-group-label' }, group.label));
    for (const route of ROUTES) {
      if (route.group !== group.id || route.hidden) continue;
      const item = h('li', { class: 'nav-item' });
      const link = h('a', { class: 'nav-link', href: `#${route.path}` }, route.label);
      if (availability && availability[route.id] === 'missing') {
        link.dataset.availability = 'pending';
        link.title = `${route.specifier} is declared but is not on disk yet`;
      }
      item.appendChild(link);
      list.appendChild(item);
    }
  }
}

async function boot() {
  const theme = createTheme();
  const badges = createBadges();
  const banner = createBanner();
  const api = createApi();

  const outlet = document.getElementById('page');
  const navList = document.getElementById('nav-list');
  const factsEl = document.getElementById('run-facts');
  const transportEl = document.getElementById('transport-indicator');
  const themeButton = document.getElementById('theme-toggle');

  let run = null;
  let health = null;

  function syncThemeButton() {
    themeButton.textContent = theme.describe();
    themeButton.setAttribute('aria-label', 'Switch colour scheme');
  }

  themeButton.addEventListener('click', () => {
    theme.cycle();
    syncThemeButton();
  });
  if (typeof window.matchMedia === 'function') {
    const query = window.matchMedia('(prefers-color-scheme: dark)');
    const onSchemeChange = () => {
      theme.invalidate();
      syncThemeButton();
    };
    if (typeof query.addEventListener === 'function') query.addEventListener('change', onSchemeChange);
  }
  syncThemeButton();

  // -- one abortable request per boot; the shell lives for the whole session.
  const bootAbort = new AbortController();

  let healthError = null;
  try {
    health = await api.get('/health', null, bootAbort.signal);
    badges.adopt(health.badges);
  } catch (error) {
    if (error && error.name === 'AbortError') return;
    health = null;
    healthError = error;
  }

  if (health && health.ok) {
    try {
      run = await api.get('/run', null, bootAbort.signal);
      banner.setRun({ mode: run.run_mode, banner: run.banner });
    } catch (error) {
      if (error && error.name === 'AbortError') return;
      run = null;
    }
    try {
      const counts = await api.get('/run/counts', null, bootAbort.signal);
      const basis = basisOf(counts, 'manifest');
      if (basis && basis.min_samples !== undefined) banner.adopt({ min_samples: basis.min_samples });
    } catch {
      /* The minimum is shown when it arrives; not having it is not fatal. */
    }
  }

  const minSamples = banner.minSamples;
  renderFacts(factsEl, health, run, minSamples, healthError);
  renderTransport(transportEl, { transport: 'connecting', reason: '' });

  // -- the transport, shared by every page through `ctx.on` (§7).
  const transport = createTransport({
    api,
    onTransport: (state) => renderTransport(transportEl, state),
  });
  transport.start();

  const on = Object.freeze({
    runState(fn) {
      return transport.subscribeState(fn);
    },
    delta(fn) {
      return transport.subscribeDelta(fn);
    },
  });

  const ctx = {
    api,
    // Filled in below and frozen before the router mounts anything, so a page
    // can never observe a half-built context.
    router: null,
    theme,
    banner,
    badges,
    on,
  };

  const { createRouter } = await import('./router.js');
  const router = createRouter({ routes: ROUTES, outlet, ctx });

  ctx.router = Object.freeze({
    go: (path, options) => router.go(path, options),
    current: () => router.current(),
    on: (event, fn) => router.on(event, fn),
    off: (event, fn) => router.off(event, fn),
  });
  Object.freeze(ctx);

  // Availability is learned from what the router actually loads, not by
  // pre-importing every module: a pre-import would run five other owners'
  // modules at start-up and would turn a missing script into an eager failure
  // nobody asked for.
  const availability = {};
  renderNav(navList, availability);
  router.on('moduleMissing', ({ route }) => {
    availability[route.id] = 'missing';
    renderNav(navList, availability);
  });
  router.on('modulePresent', ({ route }) => {
    if (availability[route.id] !== 'missing') return;
    availability[route.id] = 'present';
    renderNav(navList, availability);
  });

  router.on('navigated', () => {
    const currentPath = router.current() ? router.current().path : '';
    for (const link of navList.querySelectorAll('.nav-link')) {
      const href = link.getAttribute('href') || '';
      const path = href.replace(/^#/, '');
      if (path === currentPath) link.setAttribute('aria-current', 'page');
      else link.removeAttribute('aria-current');
    }
  });

  // -- the mode banner, on every view: a STUB run must not read as an
  //    analysis. It is the run's own sentence, not one composed here.
  if (run && run.banner) {
    const bannerEl = h('div', { class: 'mode-banner', role: 'note' });
    bannerEl.textContent = run.banner;
    outlet.insertBefore(bannerEl, outlet.firstChild);
  }
  if (healthError) {
    const box = message('warn', 'The server did not answer /api/health');
    box.appendChild(
      note(
        healthError.serverMessage ||
          'The API may not be running. Every page will report its own absence; ' +
            'nothing is shown as a zero.'
      )
    );
    if (healthError.probes) {
      box.appendChild(h('h4', null, 'The search the server reported'));
      box.appendChild(probeList(healthError.probes));
    }
    outlet.insertBefore(box, outlet.firstChild);
  }

  router.start();

  // The shell's own teardown, for a page that wants to stop the transport.
  globalThis.__paDashboard = Object.assign(globalThis.__paDashboard || {}, {
    ctx,
    router,
    transport,
    api,
    health: () => health,
    run: () => run,
  });
}

if (globalThis.__PA_DASHBOARD_BOOT__) {
  globalThis.__PA_DASHBOARD_BOOT__ = false;
  boot().catch((error) => {
    console.error('[shell] boot failed', error);
    const outlet = document.getElementById('page');
    if (!outlet) return;
    clear(outlet);
    outlet.appendChild(
      message(
        'error',
        'The dashboard shell could not start',
        note(String((error && error.message) || error))
      )
    );
  });
}
