/* pages/tables.js — browse ANY stage output, plus the shared paged-grid
 * primitive the other DATA pages are built on.
 *
 * Two things live here, deliberately:
 *
 *  1. `mount(container, ctx)` — the stage-tables browser (§5.1, §5.2).
 *  2. `createPagedGrid(...)` — one server-side paging/sorting/filtering
 *     controller, exported so `isolates.js`, `summaries.js` and `oprd.js` do
 *     not each grow their own. Two grids that page differently are two
 *     answers to "how many rows are there", and one of them will be wrong.
 *
 * The rules this file is built around:
 *
 * - **The browser never counts, sorts or filters a cohort.** Every sort,
 *   filter, search and page goes to the server as a query parameter, and every
 *   number on screen is `meta.total` / `meta.basis` as the server computed it.
 *   A client that sorts what is on the page reports the page, not the run
 *   (UI-D7). `createPagedGrid` has no sort/filter/count code path at all: the
 *   only thing it does to a response is render it.
 * - **`limit_applied`, `total_unfiltered`, `basis` and `meta.read` are shown,
 *   not swallowed.** A clamped page that says so, and a filter the server did
 *   not honour, is the difference between a page and a lie.
 * - **A filter the backend cannot answer is refused in the UI.** The stage-row
 *   endpoints page over a byte-offset index keyed on `sample_id`, so a filter
 *   on any other column is accepted and then silently NOT applied
 *   (`app._rows_via_index` drops the predicate and reports the column in
 *   `meta.read.unsupported_filter_columns`). This page reads that field and
 *   renders the refusal; it never presents the unfiltered count as a filtered
 *   one.
 */

import {
  ABSENT,
  absentSpan,
  basisOf,
  clear,
  h,
  kv,
  message,
  note,
  probeList,
  serverMessage,
} from '../app.js';

/* ------------------------------------------------------------------ limits */

/** `query.py` clamps `limit` to `MAX_LIMIT = 1000` and reports the clamp. */
export const MAX_LIMIT = 1000;
export const DEFAULT_PAGE_SIZES = Object.freeze([25, 50, 100, 200, 500, 1000]);

/**
 * Columns per page in a projected table.
 *
 * 40 names is ~800 bytes of query string, comfortably inside every server's
 * header limit, and about what a reader can scan. The feature matrix is the
 * reason this is a constant and not "however many there are".
 */
export const COLUMN_PAGE = 40;

/** Above this, a full export is refused rather than left to hang a tab. */
const EXPORT_ROW_CEILING = 100000;

/* ------------------------------------------------------------- absence ---- */

/**
 * One cell, with the D2 wording.
 *
 * The distinction this whole project turns on: a table that EXISTS and holds
 * no row for this sample is a measurement of zero; a table that does NOT exist
 * is not a measurement. The server already encodes that — `0`, `[]`, or the
 * string `not assessed` — so the only job here is to never turn an absent
 * value into a blank, a dash or a zero.
 *
 * @param {*} cell        the server's value
 * @param {string} reason  the server's sentence for its absence
 * @param {object} opts    `{emptyIsZeroText, title}`
 */
export function cell(cellValue, reason, opts = {}) {
  if (cellValue === null || cellValue === undefined) {
    return absentSpan(reason || 'the server sent no value and no reason', ABSENT.NOT_REPORTED);
  }
  if (Array.isArray(cellValue)) {
    if (cellValue.length === 0) {
      // `[]` is a FINDING in this project: the present table holds no row for
      // this sample. It must never read as an absent value, and it must never
      // read as an error either.
      return h(
        'span',
        { class: 'measured-zero', title: opts.title || reason || '' },
        opts.emptyIsZeroText || 'none recorded'
      );
    }
    return h(
      'span',
      { class: 'chip-row' },
      cellValue.map((entry) => h('span', { class: 'chip' }, String(entry)))
    );
  }
  if (cellValue === '') {
    return absentSpan(reason || 'the cell holds an empty string', ABSENT.NOT_REPORTED);
  }
  if (typeof cellValue === 'number' || typeof cellValue === 'boolean') {
    return h('span', { class: 'mono' }, String(cellValue));
  }
  return document.createTextNode(String(cellValue));
}

/** True when the server used one of its absence words rather than a value. */
export function isAbsenceWord(text) {
  if (typeof text !== 'string') return false;
  const lowered = text.trim().toLowerCase();
  return (
    lowered === 'not assessed' ||
    lowered === 'not produced' ||
    lowered === 'not reported' ||
    lowered === 'not run'
  );
}

/** Render the server's absence word in the absent style, keeping the word. */
export function absence(text, reason) {
  return absentSpan(reason || 'the server recorded this cell as absent', String(text));
}

/* ------------------------------------------------------------ D3 banner --- */

/**
 * Mount the one D3 banner, with `n` taken from `basis.n` and never from the
 * page's own count (§7.1). `noteText` is for the case where `basis.n` is not a
 * cohort size — the banner still shows the server's number, and the note says
 * what that number actually is, because silently reading "n=5, underpowered"
 * next to a pangenome computed over 900 isolates is its own kind of lie.
 */
export function mountPowerBanner(ctx, host, basis, noteText) {
  clear(host);
  if (!basis) {
    host.appendChild(
      message(
        'warn',
        'The cohort this statistic was computed on is not reported',
        note(
          'The response carried no `basis`, so the R16 flag cannot name an n. ' +
            'It is not shown with the cohort size the page was told about: ' +
            'that is the error the flag exists to prevent (UI-D3).'
        )
      )
    );
    return null;
  }
  const n = basis.n === null || basis.n === undefined ? null : Number(basis.n);
  ctx.banner.mount(host, n, basis);
  const pairs = [
    ['n (basis.n)', n === null ? absentSpan('the response recorded no n') : h('b', null, String(n))],
    ['artefact', basis.artefact ? h('code', null, basis.artefact) : absentSpan(ABSENT.NOT_REPORTED)],
    ['rows in the artefact', basis.rows_total === null || basis.rows_total === undefined
      ? absentSpan(ABSENT.NOT_REPORTED)
      : String(basis.rows_total)],
    ['computed over a filtered subset', basis.filtered ? 'yes — the count below is the filtered one' : 'no'],
  ];
  host.appendChild(kv(pairs));
  if (basis.path) {
    host.appendChild(h('div', { class: 'table-meta' }, h('span', { class: 'faint' }, 'computed from '), h('code', null, basis.path)));
  }
  if (noteText) host.appendChild(note(noteText));
  return host;
}

/* ------------------------------------------------------------------ CSV --- */

export function csvCell(value) {
  if (value === null || value === undefined) return '';
  const text = Array.isArray(value) ? value.join('; ') : String(value);
  return /[",\n]/.test(text) ? `"${text.replace(/"/g, '""')}"` : text;
}

/**
 * Offer a CSV for `rows`.
 *
 * The preamble is the point. A CSV that leaves the browser looks like the
 * cohort; a CSV whose first lines state the exact server query, the
 * server-computed `total`, and each row's `basis` does not. So the preamble
 * carries them, as `#` comment lines, which every spreadsheet treats as text.
 */
export function offerCsv(host, { filename, columns, rows, preamble }) {
  clear(host);
  const lines = [];
  for (const line of preamble || []) lines.push(`# ${line}`);
  lines.push(columns.map(csvCell).join(','));
  for (const row of rows) lines.push(columns.map((c) => csvCell(row[c])).join(','));
  const text = lines.join('\n');
  const status = h('p', { class: 'status-line' });
  host.appendChild(status);

  let url = null;
  const download = h('button', { type: 'button', class: 'primary' }, `Download ${filename}`);
  const onDownload = () => {
    if (url) return;
    try {
      const blob = new Blob([text], { type: 'text/csv;charset=utf-8' });
      url = URL.createObjectURL(blob);
    } catch (error) {
      status.textContent =
        `This browser could not produce a download for the CSV (${error && error.message}). ` +
        'The text is below instead — it is the same bytes.';
      host.appendChild(h('pre', { class: 'log-tail' }, text.slice(0, 200000)));
      return;
    }
    const link = h('a', { href: url, download: filename });
    host.appendChild(link);
    link.click();
    status.textContent = `Downloading ${filename} — ${rows.length} data row(s).`;
  };
  download.addEventListener('click', onDownload);
  host.appendChild(download);
  const PREVIEW_CHARS = 40000;
  const preview = text.length > PREVIEW_CHARS
    ? `${text.slice(0, PREVIEW_CHARS)}\n\n… the preview stops here. The download above is the whole file (${text.length} characters, ${rows.length} data row(s) plus the preamble and header).`
    : text;
  host.appendChild(
    h(
      'details',
      null,
      h('summary', null, `Show the CSV text — ${rows.length} data row(s), ${text.length} characters`),
      h('pre', { class: 'log-tail' }, preview)
    )
  );
  return {
    // Both halves of DESIGN §7: the listener goes with the button, and the
    // object URL does not outlive the page.
    dispose() {
      download.removeEventListener('click', onDownload);
      this.release();
    },
    release() {
      if (url) {
        try {
          URL.revokeObjectURL(url);
        } catch {
          /* already gone */
        }
        url = null;
      }
    },
  };
}

/* ------------------------------------------------------- the paged grid --- */

const SORT_NONE = 'none';
const SORT_ASC = 'asc';
const SORT_DESC = 'desc';

/**
 * `createPagedGrid`.
 *
 * @param {object} o
 * @param {object} o.ctx          the DESIGN §7 context
 * @param {string} o.path         endpoint, e.g. `/isolates`
 * @param {Array}  o.columns      `[{key, label, render?, sort?, title?, align?}]`
 * @param {Array}  o.filters      `[{param, label, kind, options?, hint?}]`
 * @param {object} o.state        `{offset, limit, sort, q, filters, hidden}`
 * @param {boolean} o.project     send the visible columns as `columns=`
 * @param {Function} o.onLoad     `(payload, meta) => void`, for page extras
 * @param {Function} o.paramsFor  extra static params merged into every query
 * @param {string} o.exportLabel  noun used in the CSV filename
 */
export function createPagedGrid(o) {
  const ctx = o.ctx;
  const columns = o.columns || [];
  const state = Object.assign(
    { offset: 0, limit: 100, sort: o.state && o.state.sort ? o.state.sort : null, q: null, filters: {}, hidden: new Set(), columnOffset: 0 },
    o.state || {}
  );
  if (!(state.hidden instanceof Set)) state.hidden = new Set(state.hidden || []);

  const listeners = [];
  let controller = null;
  let disposed = false;
  let token = 0;
  let payload = null;
  let csvHandle = null;
  let searchTimer = null;
  const timings = [];

  const el = h('div', { class: 'paged-grid' });
  const headHost = h('div');
  const toolbarHost = h('div', { class: 'toolbar' });
  const bannerHost = h('div');
  const alertHost = h('div');
  const metaHost = h('div', { class: 'table-meta' });
  const tableHost = h('div');
  const pagerHost = h('div', { class: 'toolbar' });
  const exportHost = h('div', { class: 'toolbar' });
  el.append(headHost, toolbarHost, bannerHost, alertHost, metaHost, tableHost, pagerHost, exportHost);
  if (o.head) headHost.appendChild(o.head);

  /* -- query construction ------------------------------------------- */

  /**
   * The columns on the current COLUMN PAGE, minus the hidden ones.
   *
   * **Columns page server-side, and this is the mechanism.** A feature matrix
   * with tens of thousands of columns cannot be requested whole: 4,001 column
   * names in a query string is a ~91 KB URL, and the framework answers that
   * with a `400` before a handler sees it (measured against the real server —
   * see the DATA report). So the projection this grid sends is bounded to
   * `COLUMN_PAGE` names, and the rest are reached with the column pager. A
   * table with more columns than that is normal; assuming otherwise is the bug
   * this guards.
   */
  function windowColumns() {
    const from = Math.max(0, state.columnOffset);
    return columns.slice(from, from + COLUMN_PAGE).filter((c) => !state.hidden.has(c.key));
  }

  function columnPages() {
    return Math.max(1, Math.ceil(columns.length / COLUMN_PAGE));
  }

  function columnPage() {
    return Math.floor(Math.max(0, state.columnOffset) / COLUMN_PAGE) + 1;
  }

  function queryParams() {
    const params = Object.assign({}, o.paramsFor ? o.paramsFor() : {});
    params.offset = state.offset;
    params.limit = state.limit;
    if (state.sort) params.sort = state.sort;
    if (state.q) params.q = state.q;
    for (const [key, value_] of Object.entries(state.filters)) {
      if (value_ === null || value_ === undefined || value_ === '') continue;
      params[`filter.${key}`] = value_;
    }
    if (o.columnFilter && state.filterColumn && state.filterValue) {
      params[`filter.${state.filterColumn}`] = state.filterValue;
    }
    if (o.project) {
      const visible = windowColumns().map((c) => c.key);
      // `query.py` refuses an unknown column, and every key here came from the
      // header the server itself reported.
      if (visible.length) params.columns = visible.join(',');
    }
    return params;
  }

  function queryDescription() {
    const bits = [`offset=${state.offset}`, `limit=${state.limit}`];
    if (state.sort) bits.push(`sort=${state.sort}`);
    if (state.q) bits.push(`q=${JSON.stringify(state.q)}`);
    if (o.columnFilter && state.filterColumn && state.filterValue) {
      bits.push(`filter.${state.filterColumn}=${JSON.stringify(state.filterValue)}`);
    }
    if (o.project) bits.push(`columns=${windowColumns().map((c) => c.key).join(',')}`);
    for (const [key, value_] of Object.entries(state.filters)) {
      if (value_ === null || value_ === undefined || value_ === '') continue;
      bits.push(`filter.${key}=${JSON.stringify(value_)}`);
    }
    return bits.join('  ');
  }

  /* -- fetch -------------------------------------------------------- */

  async function load() {
    if (disposed) return null;
    if (controller) controller.abort();
    controller = new AbortController();
    const mine = (token += 1);
    clear(alertHost);
    const started = performance.now();
    try {
      const body = await ctx.api.get(o.path, queryParams(), controller.signal);
      if (disposed || mine !== token) return null;
      payload = body;
      timings.push({ ms: performance.now() - started, rows: (body.items || body.rows || []).length });
      if (timings.length > 20) timings.shift();
      render();
      if (o.onLoad) o.onLoad(body, body.meta || null);
      return body;
    } catch (error) {
      if (error && error.name === 'AbortError') return null;
      if (disposed || mine !== token) return null;
      payload = null;
      renderError(error);
      if (o.onLoad) o.onLoad(null, null);
      return null;
    }
  }

  function renderError(error) {
    clear(tableHost);
    clear(metaHost);
    const box = message('error', `The server did not answer ${o.path}`);
    if (error && error.serverMessage) {
      box.appendChild(serverMessage(error.serverMessage));
    } else {
      box.appendChild(note(String((error && error.message) || error)));
    }
    if (error && error.status) {
      box.appendChild(note(`HTTP ${error.status}. The query was: ${queryDescription()}`));
    }
    if (error && error.probes) {
      box.appendChild(h('h4', null, 'The search the server reported'));
      box.appendChild(probeList(error.probes));
    }
    alertHost.appendChild(box);
  }

  /* -- toolbar ------------------------------------------------------ */

  function renderToolbar() {
    clear(toolbarHost);
    if (o.columnFilter) {
      // Filter on ANY column of the table, by name, with the name coming from
      // the header the server reported. A stage table's index is keyed on one
      // column, so a filter on another is accepted and then dropped by the
      // server — `renderAlerts` catches that and says so rather than showing
      // the unfiltered count as a filtered one.
      const select = h('select', { 'aria-label': 'filter on column' });
      for (const column of columns) {
        select.appendChild(h('option', { value: column.key }, column.key));
      }
      if (state.filterColumn) select.value = state.filterColumn;
      const value_ = h('input', {
        type: 'text',
        placeholder: 'exact value',
        value: state.filterValue || '',
        'aria-label': 'filter value',
        spellcheck: false,
        autocomplete: 'off',
      });
      const apply = () => {
        state.filterColumn = select.value || null;
        state.filterValue = value_.value === '' ? null : value_.value;
        state.offset = 0;
        void load();
      };
      select.addEventListener('change', apply);
      value_.addEventListener('change', apply);
      listeners.push(() => {
        select.removeEventListener('change', apply);
        value_.removeEventListener('change', apply);
      });
      toolbarHost.appendChild(h('label', { title: 'exact match; the server evaluates it, or refuses it and says which' }, 'filter column', select));
      toolbarHost.appendChild(h('label', null, 'equals', value_));
    }
    if (o.search) {
      const input = h('input', {
        type: 'search',
        placeholder: o.search.placeholder || 'search…',
        value: state.q || '',
        'aria-label': o.search.label || 'search',
        spellcheck: false,
        autocomplete: 'off',
      });
      const onInput = () => {
        if (searchTimer !== null) clearTimeout(searchTimer);
        searchTimer = setTimeout(() => {
          searchTimer = null;
          state.q = input.value.trim() || null;
          state.offset = 0;
          void load();
        }, 250);
      };
      input.addEventListener('input', onInput);
      listeners.push(() => {
        input.removeEventListener('input', onInput);
        if (searchTimer !== null) clearTimeout(searchTimer);
      });
      toolbarHost.appendChild(h('label', null, o.search.label || 'search', input));
    }

    for (const f of o.filters || []) {
      let input;
      if (f.kind === 'select') {
        input = h('select', { 'aria-label': f.label });
        input.appendChild(h('option', { value: '' }, f.anyLabel || 'any'));
        for (const option of f.options || []) {
          input.appendChild(h('option', { value: String(option), selected: state.filters[f.param] === String(option) ? true : null }, String(option)));
        }
        input.value = state.filters[f.param] || '';
      } else {
        input = h('input', {
          type: 'text',
          placeholder: f.placeholder || 'exact value',
          value: state.filters[f.param] || '',
          'aria-label': f.label,
          spellcheck: false,
          autocomplete: 'off',
        });
      }
      const onChange = () => {
        state.filters[f.param] = input.value === '' ? null : input.value;
        state.offset = 0;
        void load();
      };
      input.addEventListener('change', onChange);
      listeners.push(() => input.removeEventListener('change', onChange));
      toolbarHost.appendChild(
        h(
          'label',
          { title: f.hint || '' },
          f.label,
          input,
          f.hint ? h('span', { class: 'faint' }, ' ⓘ') : null
        )
      );
    }

    const size = h('select', { 'aria-label': 'rows per page' });
    const sizes = DEFAULT_PAGE_SIZES.includes(state.limit) ? DEFAULT_PAGE_SIZES : [...DEFAULT_PAGE_SIZES, state.limit];
    for (const value_ of sizes) {
      size.appendChild(h('option', { value: String(value_) }, `${value_} / page`));
    }
    size.value = String(state.limit);
    const onSize = () => {
      state.limit = Number(size.value);
      state.offset = 0;
      void load();
    };
    size.addEventListener('change', onSize);
    listeners.push(() => size.removeEventListener('change', onSize));
    toolbarHost.appendChild(h('label', null, 'page size', size));
    toolbarHost.appendChild(h('span', { class: 'spacer' }));

    // The column chooser. Hiding a column is a display decision and needs no
    // permission from the server; when `project` is on, the hidden columns are
    // also left out of the request, which is what keeps a 40,000-column table
    // off the wire.
    const details = h('details');
    details.appendChild(
      h(
        'summary',
        null,
        `columns (${windowColumns().length} of ${columns.length} shown, page ${columnPage()} of ${columnPages()})`
      )
    );
    const list = h('div', { class: 'toolbar' });
    for (const column of columns.slice(state.columnOffset, state.columnOffset + COLUMN_PAGE)) {
      const box = h('input', { type: 'checkbox', checked: state.hidden.has(column.key) ? null : true });
      const onToggle = () => {
        if (box.checked) state.hidden.delete(column.key);
        else state.hidden.add(column.key);
        void load();
      };
      box.addEventListener('change', onToggle);
      listeners.push(() => box.removeEventListener('change', onToggle));
      list.appendChild(h('label', { title: column.title || '' }, box, column.label));
    }
    if (columns.length > COLUMN_PAGE) {
      list.appendChild(
        h(
          'span',
          { class: 'faint' },
          `the chooser lists the ${COLUMN_PAGE} columns on the current page; use the column pager to reach the other ${columns.length - COLUMN_PAGE}`
        )
      );
    }
    const showAll = h('button', { type: 'button', class: 'ghost' }, 'show every column on this page');
    const onShowAll = () => {
      state.hidden.clear();
      void load();
    };
    showAll.addEventListener('click', onShowAll);
    listeners.push(() => showAll.removeEventListener('click', onShowAll));
    details.append(list, showAll);
    toolbarHost.appendChild(details);

    if (columns.length > COLUMN_PAGE) {
      const pages = columnPages();
      const goCols = (page) => {
        state.columnOffset = Math.max(0, Math.min(page - 1, pages - 1)) * COLUMN_PAGE;
        void load();
      };
      const mkCol = (label, page, disabled) => {
        const b = h('button', { type: 'button', disabled: disabled ? true : null }, label);
        if (!disabled) {
          const fn = () => goCols(page);
          b.addEventListener('click', fn);
          listeners.push(() => b.removeEventListener('click', fn));
        }
        return b;
      };
      toolbarHost.appendChild(
        h(
          'span',
          null,
          mkCol('‹ columns', columnPage() - 1, columnPage() === 1),
          h(
            'span',
            { class: 'faint' },
            ` columns ${state.columnOffset + 1}–${Math.min(state.columnOffset + COLUMN_PAGE, columns.length)} of ${columns.length} `,
          ),
          mkCol('columns ›', columnPage() + 1, columnPage() === pages)
        )
      );
      toolbarHost.appendChild(
        h(
          'span',
          { class: 'faint' },
          'columns page server-side: only these are requested, because a projection of ' +
            'every column in a wide feature matrix is a request the server refuses outright'
        )
      );
    }
  }

  /* -- meta line ---------------------------------------------------- */

  function renderMeta(body) {
    clear(metaHost);
    const meta = body.meta || {};
    const items = body.items || body.rows || [];
    const total = typeof meta.total === 'number' ? meta.total : null;
    const first = total === 0 ? 0 : state.offset + 1;
    const last = state.offset + items.length;
    const scope = h('span', null, `rows ${first}–${last}`);
    if (total !== null) {
      scope.appendChild(document.createTextNode(' of '));
      const strong = h('b', null, String(total));
      scope.appendChild(strong);
      scope.appendChild(
        document.createTextNode(
          ` (the server's count over the whole artefact${meta.filtered || state.q ? ', after the filter and search above' : ''}, not the ${items.length} row(s) on this page)`
        )
      );
      if (typeof meta.total_unfiltered === 'number') {
        scope.appendChild(document.createTextNode(`. Unfiltered total: ${meta.total_unfiltered}.`));
      }
    }
    metaHost.appendChild(h('span', null, scope));

    if (meta.limit_applied) {
      metaHost.appendChild(
        h(
          'span',
          { dataset: { flag: 'true' } },
          `you asked for ${state.limit} rows; the server capped the page at ${meta.limit_applied} (meta.limit_applied)`
        )
      );
    }
    const basis = meta.basis;
    if (basis) {
      // Every number reaches the reader with its denominator (UI-D7), even
      // where the panel-level banner already carries one: `total` is a
      // statement about `basis.n` rows of `basis.artefact`, and saying so here
      // costs one span.
      metaHost.appendChild(
        h(
          'span',
          null,
          `basis: n=${basis.n === null || basis.n === undefined ? ABSENT.NOT_REPORTED : basis.n}`,
          basis.artefact ? ` from ${basis.artefact}` : '',
          basis.filtered ? ' (computed over the filtered subset, not the whole artefact)' : ''
        )
      );
    }
    if (body.path) {
      metaHost.appendChild(h('span', { class: 'faint' }, 'file: '), h('code', null, body.path));
    }
    if (body.reason) {
      metaHost.appendChild(h('span', null, body.reason));
    }
    const read = meta.read;
    if (read) {
      const size = Number(read.file_size) || 0;
      const readBytes = Number(read.bytes_read) || 0;
      metaHost.appendChild(
        h(
          'span',
          { class: read.fully_loaded ? '' : 'faint' },
          read.fully_loaded
            ? `read ${readBytes} of ${size} bytes — the WHOLE file (${read.mechanism})`
            : `read ${readBytes} of ${size} bytes — not the whole file (${read.mechanism})`
        )
      );
    }
    const last3 = timings.slice(-3).map((t) => `${Math.round(t.ms)} ms`);
    if (last3.length) metaHost.appendChild(h('span', { class: 'faint' }, `page render: ${last3.join(', ')}`));
  }

  /* -- refusals ----------------------------------------------------- */

  function renderAlerts(body) {
    clear(alertHost);
    const meta = body.meta || {};
    const read = meta.read || {};
    const unsupported = Array.isArray(read.unsupported_filter_columns) ? read.unsupported_filter_columns : [];
    if (unsupported.length) {
      const box = message(
        'warn',
        'The server did not apply that filter',
        note(
          `This endpoint pages over a byte-offset index keyed on ` +
            `${read.key_column ? `\`${read.key_column}\`` : 'one column'}, so a filter on ` +
            `${unsupported.map((c) => `\`${c}\``).join(', ')} cannot be evaluated and was ` +
            'dropped. The counts above are therefore the UNFILTERED counts, not ' +
            'zero matches and not filtered ones. This page refuses to present ' +
            'them as a filtered result.'
        )
      );
      alertHost.appendChild(box);
    }
    if (body.present === false || (Array.isArray(body.items) && body.items.length === 0 && (meta.total === 0 || meta.total === undefined))) {
      if (body.reason) {
        const box = message('info', 'This table is not produced', serverMessage(body.reason));
        alertHost.appendChild(box);
      }
    }
    if (body.oprd_reason && body.oprd_available === false) {
      alertHost.appendChild(
        message('info', 'The oprD column is `not_assessed` for every isolate', serverMessage(body.oprd_reason))
      );
    }
  }

  /* -- the table ---------------------------------------------------- */

  function renderTable(body) {
    clear(tableHost);
    const items = body.items || body.rows || [];
    if (!Array.isArray(items)) {
      tableHost.appendChild(
        message('warn', 'The response carried no rows field', note('Nothing is drawn rather than an empty grid that reads as a measurement of zero.'))
      );
      return;
    }
    if (items.length === 0) {
      const total = body.meta && body.meta.total;
      const box = message('info', total === 0 ? 'No row matches' : 'This page is empty');
      box.appendChild(
        note(
          total === 0
            ? 'The server counted the whole artefact and found no row matching this query. That is a measurement — a filter that matched nothing — and it is reported as such rather than as an absent table.'
            : `The offset (${state.offset}) is past the end of the ${total} row(s) the server counted.`
        )
      );
      tableHost.appendChild(box);
      return;
    }
    const visible = windowColumns();
    const table = h('table', { class: 'data-table' });
    const thead = h('thead');
    const hr = h('tr');
    for (const column of visible) {
      const th = h('th', { scope: 'col', title: column.title || '' });
      if (column.sort === false || !column.sortKey) {
        th.appendChild(h('span', null, column.label));
      } else {
        const active = state.sort === column.sortKey || state.sort === `-${column.sortKey}`;
        const button = h(
          'button',
          { type: 'button', class: 'ghost', 'aria-label': `sort by ${column.label}` },
          h('span', null, column.label),
          ' ',
          h(
            'span',
            { class: 'mono', title: 'sort direction: glyph and word, never colour alone' },
            active ? (state.sort === `-${column.sortKey}` ? '▼ desc' : '▲ asc') : '↕ unsorted'
          )
        );
        const onSort = () => {
          const key = column.sortKey;
          if (state.sort === key) state.sort = `-${key}`;
          else if (state.sort === `-${key}`) state.sort = null;
          else state.sort = key;
          state.offset = 0;
          void load();
        };
        button.addEventListener('click', onSort);
        listeners.push(() => button.removeEventListener('click', onSort));
        th.appendChild(button);
      }
      hr.appendChild(th);
    }
    thead.appendChild(hr);
    const tbody = h('tbody');
    items.forEach((row, index) => {
      const tr = h('tr');
      for (const column of visible) {
        const td = h('td', { dataset: { column: column.key } });
        if (column.render) {
          const rendered = column.render(row, index);
          if (rendered === null || rendered === undefined) td.appendChild(absentSpan('this cell holds no value'));
          else if (rendered instanceof Node) td.appendChild(rendered);
          else td.appendChild(document.createTextNode(String(rendered)));
        } else {
          td.appendChild(cell(row[column.key], column.absentReason, { title: column.title }));
        }
        tr.appendChild(td);
      }
      tbody.appendChild(tr);
    });
    table.append(thead, tbody);
    // `th` is already `position: sticky` in base.css and `.table-wrap` is the
    // scroll container, so the header stays put with no extra stylesheet.
    tableHost.appendChild(h('div', { class: 'table-wrap' }, table));
  }

  /* -- pager -------------------------------------------------------- */

  function renderPager(body) {
    clear(pagerHost);
    const meta = body.meta || {};
    const total = typeof meta.total === 'number' ? meta.total : 0;
    const pages = Math.max(1, Math.ceil(total / state.limit));
    const page = Math.floor(state.offset / state.limit) + 1;
    const go = (offset) => {
      state.offset = Math.max(0, Math.min(offset, Math.max(0, total - 1)));
      void load();
    };
    const mk = (label, offset, disabled, title) => {
      const b = h('button', { type: 'button', disabled: disabled ? true : null, title: title || null }, label);
      if (!disabled) {
        const fn = () => go(offset);
        b.addEventListener('click', fn);
        listeners.push(() => b.removeEventListener('click', fn));
      }
      return b;
    };
    pagerHost.appendChild(mk('« first', 0, state.offset === 0));
    pagerHost.appendChild(mk('‹ prev', state.offset - state.limit, state.offset === 0));
    pagerHost.appendChild(h('span', null, `page ${page} of ${pages} · offset ${state.offset}`));
    pagerHost.appendChild(mk('next ›', state.offset + state.limit, state.offset + state.limit >= total));
    pagerHost.appendChild(mk('last »', (pages - 1) * state.limit, state.offset + state.limit >= total));
    const reset = h('button', { type: 'button', class: 'ghost' }, 'reset every filter, search and sort');
    const onReset = () => {
      state.filters = {};
      state.filterColumn = null;
      state.filterValue = null;
      state.q = null;
      state.sort = o.state && o.state.sort ? o.state.sort : null;
      state.offset = 0;
      void load();
    };
    reset.addEventListener('click', onReset);
    listeners.push(() => reset.removeEventListener('click', onReset));
    pagerHost.appendChild(reset);
  }

  /* -- export ------------------------------------------------------- */

  function renderExport(body) {
    clear(exportHost);
    if (csvHandle) {
      csvHandle.dispose();
      csvHandle = null;
    }
    const items = body.items || body.rows || [];
    if (items.length === 0) {
      exportHost.appendChild(h('span', { class: 'faint' }, 'nothing to export: this page holds no row'));
      return;
    }
    const meta = body.meta || {};
    const total = typeof meta.total === 'number' ? meta.total : items.length;
    const complete = total <= items.length;
    const noun = o.exportLabel || 'rows';
    const status = h('span', { class: 'status-line' });

    const pageOnly = h('button', { type: 'button' }, `export these ${items.length} ${noun}`);
    const onPageOnly = () => {
      csvHandle = offerCsv(exportHost, {
        filename: `${(o.exportName || 'export')}-page-only.csv`,
        columns: windowColumns().map((c) => c.key),
        rows: items,
        preamble: [
          'EXPORT SCOPE: this page only.',
          `The server counted ${total} row(s) for this query; this file holds ${items.length}.`,
          'This is NOT the whole table.',
          `Query: ${queryDescription()}`,
        ],
      });
      status.textContent = '';
    };
    pageOnly.addEventListener('click', onPageOnly);
    listeners.push(() => pageOnly.removeEventListener('click', onPageOnly));
    exportHost.appendChild(pageOnly);

    if (!complete) {
      const all = h('button', { type: 'button', class: 'primary' }, `export all ${total} ${noun} (server-paged)`);
      const onAll = () => {
        all.disabled = true;
        status.textContent = `fetching all ${total} row(s) from the server, ${MAX_LIMIT} at a time…`;
        void exportAll(total, status, all);
      };
      all.addEventListener('click', onAll);
      listeners.push(() => all.removeEventListener('click', onAll));
      exportHost.appendChild(all);
      exportHost.appendChild(
        h(
          'span',
          { class: 'faint' },
          'the export pages through the same server-side query and assembles the file here; nothing is counted in the browser'
        )
      );
    }
    exportHost.appendChild(status);
  }

  async function exportAll(total, status, button) {
    if (total > EXPORT_ROW_CEILING) {
      status.textContent =
        `refused: ${total} rows is above the ${EXPORT_ROW_CEILING}-row export ceiling. A file that large ` +
        'would have to be streamed by the server instead, and no endpoint streams one. The current page is ' +
        'still exportable and is labelled as a page.';
      button.disabled = false;
      return;
    }
    const collected = [];
    try {
      for (let offset = 0; offset < total; offset += MAX_LIMIT) {
        const params = queryParams();
        params.offset = offset;
        params.limit = MAX_LIMIT;
        params.include_total = 1;
        const body = await ctx.api.get(o.path, params, controller ? controller.signal : undefined);
        const rows = body.items || body.rows || [];
        if (rows.length === 0) break;
        collected.push(...rows);
        status.textContent = `fetched ${collected.length} of ${total}…`;
      }
    } catch (error) {
      if (error && error.name === 'AbortError') return;
      status.textContent = `the export stopped: ${(error && error.message) || error}. Nothing partial was offered as complete.`;
      button.disabled = false;
      return;
    }
    csvHandle = offerCsv(exportHost, {
      filename: `${(o.exportName || 'export')}-all.csv`,
      columns: windowColumns().map((c) => c.key),
      rows: collected,
      preamble: [
        'EXPORT SCOPE: the complete server-side result set for the query below.',
        `The server counted ${total} row(s); this file holds ${collected.length}.`,
        `Query: ${queryDescription()}`,
        `Counted by the server: ${JSON.stringify(meta_basis_of_export())}`,
      ],
    });
    status.textContent = `assembled ${collected.length} of ${total} row(s).`;
    button.disabled = false;
  }

  function meta_basis_of_export() {
    const meta = payload && payload.meta;
    return meta ? meta.basis : null;
  }

  /* -- render ------------------------------------------------------- */

  function render() {
    if (!payload) return;
    renderToolbar();
    renderAlerts(payload);
    renderMeta(payload);
    renderTable(payload);
    renderPager(payload);
    renderExport(payload);
    const basis = basisOf(payload);
    if (o.banner !== false) mountPowerBanner(ctx, bannerHost, basis, o.bannerNote);
  }

  renderToolbar();

  return {
    el,
    load,
    get payload() {
      return payload;
    },
    get state() {
      return state;
    },
    setFilter(param, value_) {
      state.filters[param] = value_;
      state.offset = 0;
      return load();
    },
    setHidden(key, hidden) {
      if (hidden) state.hidden.add(key);
      else state.hidden.delete(key);
      return load();
    },
    timings: () => timings.slice(),
    destroy() {
      disposed = true;
      if (searchTimer !== null) clearTimeout(searchTimer);
      if (controller) controller.abort();
      if (csvHandle) csvHandle.dispose();
      for (const fn of listeners.splice(0)) {
        try {
          fn();
        } catch (error) {
          console.error('[grid] a disposer threw', error);
        }
      }
    },
  };
}

/* ====================================================================== */
/* The stage-tables browser — `mount(container, ctx)`                      */
/* ====================================================================== */

/**
 * The 16 stage keys, in `STAGE_ORDER`, plus the four folded steps and the four
 * Snakefile-declared paths.
 *
 * The stage list is not hard-coded here: it comes from `/api/stages`, which
 * reads `STAGE_ORDER`, and the folded steps ride along on each stage's
 * `internal_tables`. The Snakefile-declared keys are discovered from
 * `/api/pangenome`'s `gene_tables`, which is the only endpoint that lists
 * them. All three are server-declared, so this page cannot fall out of step
 * with the pipeline.
 */
async function discoverTables(ctx, signal) {
  const [stages, pangenome] = await Promise.all([
    ctx.api.get('/stages', null, signal).catch((error) => (error && error.name === 'AbortError' ? null : { __error: error })),
    ctx.api.get('/pangenome', null, signal).catch((error) => (error && error.name === 'AbortError' ? null : { __error: error })),
  ]);
  const keys = [];
  const seen = new Set();
  const add = (key, group, owner) => {
    if (!key || seen.has(key)) return;
    seen.add(key);
    keys.push({ key, group, owner });
  };
  const stageItems = stages && Array.isArray(stages.items) ? stages.items : [];
  for (const stage of stageItems) {
    add(stage.name, 'stage', `${stage.spec_name} (spec D1 #${stage.spec_number})`);
    for (const folded of stage.internal_tables || []) add(folded.key, 'folded step', `folded into ${folded.owner}`);
  }
  const geneTables = pangenome && Array.isArray(pangenome.gene_tables) ? pangenome.gene_tables : [];
  for (const table of geneTables) add(table.key, 'Snakefile-declared', table.declared_by);
  return { keys, stages: stages, pangenomeError: pangenome && pangenome.__error ? pangenome.__error : null, stageError: stages && stages.__error ? stages.__error : null };
}

export function mount(container, ctx) {
  const abort = new AbortController();
  let disposed = false;
  const disposers = [];

  const head = h(
    'div',
    { class: 'page-head' },
    h(
      'div',
      { class: 'page-head-text' },
      h('h2', null, 'Stage tables'),
      h(
        'p',
        { class: 'page-note' },
        'Every stage output, folded step and Snakefile-declared table the pipeline ' +
          'writes — paged, sorted, searched and filtered by the server. Columns ' +
          'page too: a feature matrix with tens of thousands of columns is ' +
          'requested one page of columns at a time, never in full.'
      )
    )
  );
  container.appendChild(head);

  const alertHost = h('div');
  const pickerHost = h('div', { class: 'panel' });
  const gridHost = h('div');
  container.append(alertHost, pickerHost, gridHost);

  let keys = [];
  let stageItems = [];
  let grid = null;
  let current = null;

  function renderPicker() {
    clear(pickerHost);
    pickerHost.appendChild(h('h3', null, 'Which table'));
    if (keys.length === 0) {
      pickerHost.appendChild(
        message(
          'warn',
          'The server declared no table key',
          note(
            'Neither /api/stages nor /api/pangenome answered with a table list, so this page ' +
              'has nothing to browse. That is reported rather than replaced with an empty grid.'
          )
        )
      );
      return;
    }
    const groups = [
      ['stage', 'Stage outputs (papipeline/run.py STAGE_ORDER)'],
      ['folded step', 'Folded steps (contracts.INTERNAL_TABLES) — written, consumed, not stages'],
      ['Snakefile-declared', 'Snakefile-declared paths (report_tables.SNAKEFILE_ONLY_TABLES)'],
    ];
    for (const [group, label] of groups) {
      const inGroup = keys.filter((entry) => entry.group === group);
      if (inGroup.length === 0) continue;
      pickerHost.appendChild(h('h4', null, `${label} — ${inGroup.length}`));
      const row = h('div', { class: 'toolbar' });
      for (const entry of inGroup) {
        const stage = stageItems.find((s) => s.name === entry.key);
        const button = h(
          'button',
          { type: 'button', dataset: { key: entry.key, selected: current === entry.key ? 'true' : null }, title: entry.owner },
          entry.key,
          stage && stage.state ? h('span', { class: 'badge', dataset: { tone: ctx.badges.tone(stage.state), badge: stage.state } }, ` ${ctx.badges.state(stage.state).glyph} ${ctx.badges.state(stage.state).label}`) : null
        );
        const onPick = () => openTable(entry.key);
        button.addEventListener('click', onPick);
        disposers.push(() => button.removeEventListener('click', onPick));
        row.appendChild(button);
      }
      pickerHost.appendChild(row);
    }
  }

  function openTable(key) {
    current = key;
    clear(gridHost);
    clear(alertHost);
    if (grid) {
      grid.destroy();
      grid = null;
    }
    // A stage's own header comes from `/api/stages/{name}` — free, because the
    // endpoint reports `actual_header` without reading a single row. Anything
    // else is discovered by one unprojected probe row, which is the only
    // request on this page that can materialise a whole row.
    const stage = stageItems.find((s) => s.name === key);
    const preset = stage && Array.isArray(stage.declared_columns) && stage.declared_columns.length
      ? stage.declared_columns
      : null;
    const probe = h('div', { class: 'panel' }, h('p', { class: 'panel-note' }, `reading the header of ${key}…`));
    gridHost.appendChild(probe);

    (async () => {
      let header = preset;
      if (!header) {
        try {
          const body = await ctx.api.get(`/tables/${encodeURIComponent(key)}`, { limit: 1 }, abort.signal);
          header = Array.isArray(body.header) ? body.header : [];
          if (body.present === false) {
            clear(gridHost);
            const box = message('warn', `${key} is not produced`, serverMessage(body.reason || 'the server gave no reason'));
            alertHost.appendChild(box);
            return;
          }
        } catch (error) {
          if (error && error.name === 'AbortError') return;
          clear(gridHost);
          alertHost.appendChild(message('error', `the header of ${key} could not be read`, note(String(error.message || error))));
          return;
        }
      }
      if (disposed) return;
      clear(gridHost);
      const columns = header.map((name) => ({
        key: name,
        label: name,
        sortKey: name,
        title: `column \`${name}\` as the file itself carries it (DESIGN §3.6.3: the file's header wins over the contract, as information)`,
      }));
      grid = createPagedGrid({
        ctx,
        path: `/tables/${encodeURIComponent(key)}`,
        paramsFor: () => ({}),
        columns,
        project: true,
        state: { offset: 0, limit: 50, sort: null, q: null, filters: {}, hidden: new Set() },
        search: {
          label: 'search (server-side substring)',
          placeholder: 'substring over every column',
        },
        columnFilter: true,
        filters: [],
        exportName: key,
        exportLabel: 'rows',
        head: h(
          'div',
          null,
          h('h3', null, `${key} — server-side paging`),
          h(
            'p',
            { class: 'panel-note' },
            `${header.length} column(s) in the file. Only the ${Math.min(COLUMN_PAGE, header.length)} ` +
              'columns on the current column page are requested, so a wide feature matrix costs one ' +
              'page of columns rather than all of them — a projection naming every column is a ' +
              'request the server refuses outright. Sorting, searching and paging are the ' +
              "server's; this page counts nothing."
          )
        ),
        bannerNote:
          'For a stage table the reported n is the number of ROWS in that table, not the cohort ' +
          'the study is about. It is shown because it is what the server reported; it is not a ' +
          'statement about how many isolates were studied.',
      });
      gridHost.appendChild(grid.el);
      renderPicker();
      void grid.load();
    })().catch((error) => {
      if (error && error.name === 'AbortError') return;
      console.error('[tables] failed', error);
    });
  }

  (async () => {
    try {
      const discovered = await discoverTables(ctx, abort.signal);
      if (disposed) return;
      keys = discovered.keys;
      stageItems = discovered.stages && Array.isArray(discovered.stages.items) ? discovered.stages.items : [];
      if (discovered.stageError) {
        alertHost.appendChild(
          message('warn', '/api/stages did not answer', note(`Only the tables it declares are listed. ${discovered.stageError.message || discovered.stageError}`))
        );
      }
      renderPicker();
      if (keys.length) openTable(keys[0].key);
    } catch (error) {
      if (error && error.name === 'AbortError') return;
      clear(pickerHost);
      pickerHost.appendChild(message('error', 'The table list could not be read', note(String((error && error.message) || error))));
    }
  })();

  return () => {
    disposed = true;
    abort.abort();
    if (grid) grid.destroy();
    for (const fn of disposers.splice(0)) {
      try {
        fn();
      } catch (error) {
        console.error('[tables] a disposer threw', error);
      }
    }
  };
}

export default { mount };
