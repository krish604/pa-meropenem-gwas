/* pages/summaries.js — the five summaries, each with its basis and its flag.
 *
 * AMR (genes vs point mutations), virulence, variants (SNV vs non-SNV, per
 * isolate), pangenome (core / shell / cloud) and phenotype provenance.
 *
 * **No panel computes a count, a rate, a total or a share.** (UI-D7.) The
 * browser is not allowed to know how many isolates carry a gene, so it does not
 * count them: it asks the server, with the same filters, and gets a number with
 * a `basis` (UI-D7). Where no endpoint computes an aggregate the panel says so
 * and shows the underlying values instead of a client-side sum. That is the
 * difference between "4 of 20 isolates carry blaKPC-2" (measured by the server,
 * `basis.n` attached) and "about a fifth of the isolates I could see carry
 * blaKPC-2" (a count of a page, which is the failure D7 exists to stop).
 *
 * Every panel that displays a statistic mounts the one D3 banner, with `n` taken
 * from `basis.n` — never from a count on the page.
 *
 * The phenotype panel carries one specific honesty requirement: the phenotype
 * table has no column recording the **testing standard** (CLSI/EUCAST and the
 * breakpoint version behind an `R`/`I`/`S` call). When no such column exists in
 * the header the file actually carries, the panel says `testing standard not
 * recorded` and names the columns it looked at. It does not print `CLSI`, and it
 * does not print `none`.
 */

import { ABSENT, absentSpan, basisOf, clear, h, kv, message, note } from '../app.js';
import { cell, createPagedGrid, isAbsenceWord, mountPowerBanner } from './tables.js';

/**
 * The susceptibility categories DESIGN §5.2 declares for the run's antibiotic.
 *
 * This is a *declared vocabulary*, not a measurement, and it is labelled as one
 * in the UI. The counts that go with it come from the server, one request per
 * category, so no category's count is derived from another.
 */
const SIR_CATEGORIES = ['R', 'I', 'S', 'SDD', 'ND'];

/** Columns that would record a susceptibility testing standard, if one were written. */
const STANDARD_COLUMN_HINTS = [
  'standard',
  'standard_used',
  'standard_version',
  'breakpoint',
  'breakpoints',
  'guideline',
  'guidelines',
  'method',
  'method_standard',
  'clsi',
  'eucast',
  'ast_standard',
  'interpretation',
  'interpretive_criteria',
  'susceptibility_standard',
];

function hasStandardColumn(header) {
  const lowered = (header || []).map((c) => String(c).toLowerCase());
  return STANDARD_COLUMN_HINTS.filter((hint) => lowered.some((c) => c === hint || c.includes(hint)));
}

export function mount(container, ctx) {
  const abort = new AbortController();
  let disposed = false;
  const grids = [];
  const disposers = [];

  container.appendChild(
    h(
      'div',
      { class: 'page-head' },
      h(
        'div',
        { class: 'page-head-text' },
        h('h2', null, 'Summaries'),
        h(
          'p',
          { class: 'page-note' },
          'Five views over the run\'s results. Every number here was computed by ' +
            'the server over the whole artefact and arrives with the artefact it was ' +
            'computed from; where the pipeline records no such number, this page says ' +
            'so rather than counting the rows it can see.'
        )
      )
    )
  );

  const alertHost = h('div');
  const amrHost = h('div', { class: 'panel' });
  const virHost = h('div', { class: 'panel' });
  const varHost = h('div', { class: 'panel' });
  const panHost = h('div', { class: 'panel' });
  const pheHost = h('div', { class: 'panel' });
  container.append(alertHost, amrHost, virHost, varHost, panHost, pheHost);

  const banner = (host, basis, text_) => {
    const box = h('div');
    host.insertBefore(box, host.firstChild);
    mountPowerBanner(ctx, box, basis, text_);
  };

  const gridHost = (host) => {
    const box = h('div');
    host.appendChild(box);
    return box;
  };

  const addGrid = (host, options) => {
    const grid = createPagedGrid(Object.assign({ banner: false }, options));
    const host_ = gridHost(host);
    host_.appendChild(grid.el);
    grids.push(grid);
    void grid.load();
    return grid;
  };

  /* -- shared: the isolate table, restricted to a few columns ---------- */

  const isolateColumns = (columns) =>
    columns.map((c) =>
      Object.assign(
        {
          render: (row) => {
            const value_ = row[c.key];
            if (value_ === null || value_ === undefined) return absentSpan(`no value for \`${c.key}\``, ABSENT.NOT_ASSESSED);
            if (isAbsenceWord(value_)) return h('span', { class: 'absent', title: 'the server marked this cell not assessed' }, value_);
            return cell(value_, '', { title: c.title || '' });
          },
          title: c.title || '',
        },
        c
      )
    );

  const SAMPLE_COLUMN = {
    key: 'sample_id',
    label: 'isolate',
    sortKey: 'sample_id',
    render: (row) => h('a', { href: `#/isolates/${encodeURIComponent(row.sample_id)}` }, row.sample_id),
    title: 'open the full per-isolate detail',
  };

  /* ==================================================== 1. AMR ========= */

  async function renderAmr() {
    clear(amrHost);
    amrHost.appendChild(h('h3', null, 'AMR — genes vs point mutations'));
    amrHost.appendChild(
      note(
        '`04_amr.tsv` holds one row per (isolate, determinant) and the `variant` column is ' +
          'what separates an element the tool reported as mutated from one it did not. ' +
          'A non-empty `variant` means MUTATED — it does not mean disruptive, and ' +
          '`determinant_type` is `AMR` in both cases.'
      )
    );

    let stageBody = null;
    let isolates = null;
    try {
      [stageBody, isolates] = await Promise.all([
        ctx.api.get('/stages/amr/rows', { limit: 5 }, abort.signal),
        ctx.api.get('/isolates', { limit: 1 }, abort.signal),
      ]);
    } catch (error) {
      if (error && error.name === 'AbortError') return;
      amrHost.appendChild(message('error', 'The AMR summary could not be read', note(String((error && error.message) || error))));
      return;
    }
    if (disposed) return;

    banner(
      amrHost,
      basisOf(isolates),
      'n is the isolate count the isolate table was built over. The row count of ' +
        '`04_amr.tsv` is a different number and is shown separately, because a table with ' +
        'three rows per isolate has three times the rows of the cohort.'
    );

    const tableTotal = stageBody && stageBody.meta ? stageBody.meta.total : null;
    amrHost.appendChild(
      kv([
        [
          'rows in 04_amr.tsv (server count)',
          typeof tableTotal === 'number' ? h('b', null, String(tableTotal)) : absentSpan(ABSENT.NOT_REPORTED),
        ],
        ['path', stageBody && stageBody.path ? h('code', null, stageBody.path) : absentSpan(ABSENT.NOT_REPORTED)],
        [
          'the table itself',
          stageBody && stageBody.present === false
            ? absentSpan(stageBody.reason, ABSENT.NOT_PRODUCED)
            : note('the first rows are browsable in the tables browser; the sample below is the head of the table.'),
        ],
      ])
    );
    if (stageBody && stageBody.present !== false && Array.isArray(stageBody.items) && stageBody.items.length) {
      const table = h('table', { class: 'data-table' });
      const header = stageBody.header || [];
      const hr = h('tr');
      for (const column of header) hr.appendChild(h('th', { scope: 'col' }, column));
      table.appendChild(h('thead', null, hr));
      const tbody = h('tbody');
      for (const row of stageBody.items) {
        const tr = h('tr');
        for (const column of header) {
          const td = h('td');
          const value_ = row[column];
          if (value_ === null || value_ === undefined || value_ === '') td.appendChild(absentSpan(`the \`${column}\` cell holds no value`));
          else td.textContent = Array.isArray(value_) ? value_.join(', ') : String(value_);
          tr.appendChild(td);
        }
        tbody.appendChild(tr);
      }
      table.appendChild(tbody);
      amrHost.appendChild(h('div', { class: 'table-wrap' }, table));
    }

    // -- the aggregate the server will compute, and only for one gene ----
    const geneHost = h('div', { class: 'panel' });
    geneHost.appendChild(h('h4', null, 'Isolates carrying one named gene — a server count'));
    const input = h('input', { type: 'text', placeholder: 'oprD', spellcheck: false, autocomplete: 'off', 'aria-label': 'gene name' });
    const button = h('button', { type: 'button', class: 'primary' }, 'count isolates carrying it');
    const status = h('p', { class: 'status-line' });
    const onCount = async () => {
      const gene = input.value.trim();
      clear(status);
      if (!gene) {
        status.appendChild(absentSpan('no gene was typed, so nothing was counted', 'nothing to count'));
        return;
      }
      button.disabled = true;
      status.textContent = `asking the server for the isolates whose amr_genes cell contains "${gene}"…`;
      try {
        const body = await ctx.api.get('/isolates', { 'filter.amr_genes__in': gene, limit: 1 }, abort.signal);
        clear(status);
        const total = body.meta && body.meta.total;
        const b = body.meta && body.meta.basis;
        status.appendChild(
          h(
            'span',
            null,
            h('b', null, typeof total === 'number' ? String(total) : ABSENT.NOT_REPORTED),
            ` isolate(s) carry `,
            h('code', null, gene),
            ' — counted by the server over ',
            b && typeof b.n === 'number' ? `${b.n} row(s)` : 'the isolate table',
            '. `filter.<column>__in` is a list membership test, so this means "this gene is ' +
              'in the cell", not "the cell equals this gene".'
          )
        );
        const bannerBox = h('div');
        status.appendChild(bannerBox);
        mountPowerBanner(
          ctx,
          bannerBox,
          b,
          'n is the filtered isolate count, which is the denominator this number was ' +
            'computed over. A rate would need the unfiltered cohort as well and is not ' +
            'computed on this page.'
        );
        const link = h('a', { href: `#/isolates?q=${encodeURIComponent(gene)}` }, 'see those isolates in the isolate table');
        status.appendChild(h('div', null, link));
      } catch (error) {
        if (error && error.name === 'AbortError') return;
        status.textContent = '';
        status.appendChild(message('error', 'the server refused that query', note(String((error && error.message) || error))));
      } finally {
        button.disabled = false;
      }
    };
    const onKey = (event) => {
      if (event.key === 'Enter') {
        event.preventDefault();
        void onCount();
      }
    };
    button.addEventListener('click', onCount);
    input.addEventListener('keydown', onKey);
    disposers.push(() => {
      button.removeEventListener('click', onCount);
      input.removeEventListener('keydown', onKey);
    });
    geneHost.appendChild(h('div', { class: 'toolbar' }, h('label', null, 'gene', input), button));
    geneHost.appendChild(status);
    amrHost.appendChild(geneHost);

    amrHost.appendChild(
      message(
        'info',
        'There is no server-side genes-vs-point-mutations aggregate',
        note(
          'No endpoint in the contract counts "isolates with at least one point mutation" or ' +
            '"isolates with at least one gene". DESIGN §5.2 promises three filters for exactly ' +
            'that — `has_amr_gene`, `has_point_mutation`, `virulence_min` — and the backend ' +
            'accepts them, but no isolate row carries those keys, so each one matches nothing ' +
            'and returns a server count of 0. This page therefore shows the two columns per ' +
            'isolate, one named gene at a time, and does not sum them. See the DATA report.'
        )
      )
    );

    addGrid(amrHost, {
      ctx,
      path: '/isolates',
      project: false,
      columns: isolateColumns([
        SAMPLE_COLUMN,
        { key: 'amr_genes', label: 'AMR genes', sortKey: 'amr_genes', title: '`[]` (none recorded) means the present table holds no row for this isolate' },
        { key: 'amr_point_mutations', label: 'point mutations (reported as mutated)', sortKey: 'amr_point_mutations' },
      ]),
      state: { offset: 0, limit: 25, sort: 'sample_id', q: null, filters: {}, hidden: new Set() },
      search: { label: 'search isolates', placeholder: 'gene or mutation' },
      exportName: 'amr-genes-and-mutations',
      exportLabel: 'isolates',
      head: h('h4', null, 'Genes and point mutations, per isolate (paged by the server)'),
    });
  }

  /* ============================================== 2. Virulence ========= */

  async function renderVirulence() {
    clear(virHost);
    virHost.appendChild(h('h3', null, 'Virulence'));
    virHost.appendChild(
      note(
        '`08_virulence.tsv` is one row per (isolate, virulence factor). A factor carried ' +
          'by no isolate in this run produces no row, so the absence of a row is a finding and ' +
          'not a gap in the record.'
      )
    );
    let stageBody = null;
    let isolates = null;
    try {
      [stageBody, isolates] = await Promise.all([
        ctx.api.get('/stages/virulence/rows', { limit: 1 }, abort.signal),
        ctx.api.get('/isolates', { limit: 1 }, abort.signal),
      ]);
    } catch (error) {
      if (error && error.name === 'AbortError') return;
      virHost.appendChild(message('error', 'The virulence summary could not be read', note(String((error && error.message) || error))));
      return;
    }
    if (disposed) return;
    banner(virHost, basisOf(isolates), 'n is the isolate count. The row count of the virulence table is a different quantity and is shown as such.');
    const sources = (isolates && isolates.contributing_sources) || {};
    virHost.appendChild(
      kv([
        [
          'rows in 08_virulence.tsv (server count)',
          stageBody && typeof stageBody.meta && typeof stageBody.meta.total === 'number'
            ? h('b', null, String(stageBody.meta.total))
            : absentSpan(ABSENT.NOT_REPORTED),
        ],
        [
          'is the virulence table there at all?',
          sources.virulence === true
            ? h('span', { class: 'chip' }, 'present — so 0 factors for an isolate is a measurement')
            : absentSpan('absent — so a per-isolate count would be `not assessed`, not 0'),
        ],
        ['path', stageBody && stageBody.path ? h('code', null, stageBody.path) : absentSpan(ABSENT.NOT_REPORTED)],
      ])
    );
    virHost.appendChild(
      note(
        'No panel here adds up the per-isolate counts: the server reports no total, and a ' +
          'client-side sum over the visible page would report the page rather than the run ' +
          '(UI-D7). The per-isolate values are below, each with its own basis.'
      )
    );
    addGrid(virHost, {
      ctx,
      path: '/isolates',
      project: false,
      columns: isolateColumns([
        SAMPLE_COLUMN,
        { key: 'n_virulence', label: 'virulence factors', sortKey: 'n_virulence', title: '`0` is a real count when 08_virulence.tsv is present' },
      ]),
      state: { offset: 0, limit: 25, sort: '-n_virulence', q: null, filters: {}, hidden: new Set() },
      search: { label: 'search isolates', placeholder: 'sample id' },
      exportName: 'virulence',
      exportLabel: 'isolates',
      head: h('h4', null, 'Virulence factors per isolate (paged by the server)'),
    });
  }

  /* ================================================= 3. Variants ======= */

  async function renderVariants() {
    clear(varHost);
    varHost.appendChild(h('h3', null, 'Variants — SNV and non-SNV, per isolate'));
    varHost.appendChild(
      note(
        'The per-isolate SNV and non-SNV counts come from `variants_provenance.json`, whose ' +
          'location is not contracted (DESIGN A8). When it exists it writes every key, so `0` ' +
          'is a real measurement; when it does not, both columns read `not assessed` — because ' +
          '"called nothing" and "never screened" must not look alike.'
      )
    );
    let stageBody = null;
    let isolates = null;
    try {
      [stageBody, isolates] = await Promise.all([
        ctx.api.get('/stages/variants/rows', { limit: 1 }, abort.signal),
        ctx.api.get('/isolates', { limit: 1 }, abort.signal),
      ]);
    } catch (error) {
      if (error && error.name === 'AbortError') return;
      varHost.appendChild(message('error', 'The variant summary could not be read', note(String((error && error.message) || error))));
      return;
    }
    if (disposed) return;
    banner(varHost, basisOf(isolates), 'n is the isolate count the per-isolate variant columns were joined over.');
    varHost.appendChild(
      kv([
        [
          'rows in variants.tsv (server count)',
          stageBody && typeof stageBody.meta && typeof stageBody.meta.total === 'number'
            ? h('b', null, String(stageBody.meta.total))
            : absentSpan(ABSENT.NOT_REPORTED),
        ],
        [
          'variant accounting file',
          isolates && isolates.variants_reason
            ? absentSpan(isolates.variants_reason, ABSENT.NOT_ASSESSED)
            : h('span', { class: 'chip' }, 'variants_provenance.json read — every key is present at zero'),
        ],
        ['path', stageBody && stageBody.path ? h('code', null, stageBody.path) : absentSpan(ABSENT.NOT_REPORTED)],
      ])
    );
    varHost.appendChild(
      note(
        'No total SNV count and no SNV/non-SNV ratio is shown: no endpoint computes either, ' +
          'and a ratio assembled from the visible page would be a ratio of the page (UI-D7).'
      )
    );
    addGrid(varHost, {
      ctx,
      path: '/isolates',
      project: false,
      columns: isolateColumns([
        SAMPLE_COLUMN,
        { key: 'n_snv', label: 'SNVs', sortKey: 'n_snv' },
        { key: 'n_indel', label: 'non-SNV variants', sortKey: 'n_indel' },
      ]),
      state: { offset: 0, limit: 25, sort: '-n_snv', q: null, filters: {}, hidden: new Set() },
      search: { label: 'search isolates', placeholder: 'sample id' },
      exportName: 'variants',
      exportLabel: 'isolates',
      head: h('h4', null, 'SNV and non-SNV counts per isolate (paged by the server)'),
    });
  }

  /* ================================================ 4. Pangenome ======= */

  async function renderPangenome() {
    clear(panHost);
    panHost.appendChild(h('h3', null, 'Pangenome — core, shell and cloud'));
    let body = null;
    try {
      body = await ctx.api.get('/pangenome', null, abort.signal);
    } catch (error) {
      if (error && error.name === 'AbortError') return;
      panHost.appendChild(message('error', 'The pangenome summary could not be read', note(String((error && error.message) || error))));
      return;
    }
    if (disposed) return;
    banner(
      panHost,
      body.basis,
      'CAREFUL: `/api/pangenome` sets `basis.n` to the number of ROWS in ' +
        '`pangenome_summary.tsv` — the number of metrics — not the cohort the pangenome was ' +
        'computed over. The flag above therefore reads n=<metric rows>. That is the server\'s ' +
        'number and this page does not substitute a different one, because a client that ' +
        'replaces basis.n with a cohort size it was told about is the exact error R16 exists to ' +
        'prevent. Where the summary itself records a cohort (`n_samples`), it is shown below as ' +
        'a VALUE and clearly marked as the file\'s own metric, not as the flag\'s denominator.'
    );

    if (!body.present) {
      panHost.appendChild(absentSpan(body.reason, ABSENT.NOT_PRODUCED));
      return;
    }
    const summary = body.summary || {};
    const entries = Object.entries(summary);
    const classed = entries.filter(([metric]) => /core|shell|cloud|pangenome|total|sample/i.test(metric));
    panHost.appendChild(
      kv(
        (classed.length ? classed : entries).map(([metric, value_]) => [
          h('code', null, metric),
          value_ === null || value_ === undefined ? absentSpan(ABSENT.NOT_REPORTED) : h('span', { class: 'mono' }, String(value_)),
        ])
      )
    );
    panHost.appendChild(
      note(
        'Core / shell / cloud are shown as the file\'s own counts. Their SHARES are not ' +
          'computed here: a share is a rate, UI-D7 forbids computing one in the browser, and no ' +
          'endpoint reports one. A reader who wants a share can divide — and should carry the ' +
          'denominator with it.'
      )
    );

    panHost.appendChild(h('h4', null, 'Gene tables'));
    const geneTables = Array.isArray(body.gene_tables) ? body.gene_tables : [];
    panHost.appendChild(
      kv(
        geneTables.map((table) => [
          h('code', null, table.key),
          table.present
            ? h('span', null, h('b', null, String(table.n_rows)), ' row(s) — counted by the server, not by this page')
            : absentSpan(table.reason, ABSENT.NOT_PRODUCED),
        ])
      )
    );
    panHost.appendChild(
      note(
        'The presence/absence table is the one that grows quadratically with the cohort, which ' +
          'is why the endpoint counts it instead of paging it. The counts here are ' +
          '`n_rows` over the whole file.'
      )
    );
    const rest = entries.filter((entry) => !classed.includes(entry));
    if (rest.length) {
      panHost.appendChild(
        note(
          `Other metrics the file records, verbatim: ${rest.map(([metric, value_]) => `${metric}=${value_}`).join(', ')}.`
        )
      );
    }
  }

  /* ============================================ 5. Phenotype =========== */

  async function renderPhenotype() {
    clear(pheHost);
    pheHost.appendChild(h('h3', null, 'Phenotype provenance'));
    let stageBody = null;
    let isolates = null;
    try {
      [stageBody, isolates] = await Promise.all([
        ctx.api.get('/stages/phenotype/rows', { limit: 5 }, abort.signal),
        ctx.api.get('/isolates', { limit: 1 }, abort.signal),
      ]);
    } catch (error) {
      if (error && error.name === 'AbortError') return;
      pheHost.appendChild(message('error', 'The phenotype summary could not be read', note(String((error && error.message) || error))));
      return;
    }
    if (disposed) return;
    banner(
      pheHost,
      basisOf(isolates),
      'n is the isolate count the susceptibility column was joined over. Each category ' +
        'count below carries its own n, because a category count is computed over the ' +
        'isolates in that category.'
    );

    if (!stageBody || stageBody.present === false) {
      pheHost.appendChild(absentSpan((stageBody && stageBody.reason) || 'the phenotype table could not be read', ABSENT.NOT_PRODUCED));
    } else {
      const header = stageBody.header || [];
      const standardColumns = hasStandardColumn(header);
      pheHost.appendChild(
        kv([
          [
            'testing standard',
            standardColumns.length
              ? h('span', { class: 'chip' }, standardColumns.join(', '))
              : h('span', { class: 'chip', dataset: { flag: 'true' } }, 'testing standard not recorded'),
          ],
          [
            'why',
            standardColumns.length
              ? note('the table carries a column that records it; its values are browsable below')
              : note(
                  'No column in `11_phenotype.tsv` records the susceptibility testing standard ' +
                    '(CLSI/EUCAST and the breakpoint version) behind an R/I/S call. The columns ' +
                    'this page looked at are below. This is stated rather than filled in with a ' +
                    'guess: an `R` without its standard and breakpoints is not a finding, and ' +
                    'the pipeline never recorded which standard the cohort was read against.'
                ),
          ],
          ['columns in the table', h('code', null, header.join(', '))],
          [
            'rows (server count)',
            stageBody.meta && typeof stageBody.meta.total === 'number' ? h('b', null, String(stageBody.meta.total)) : absentSpan(ABSENT.NOT_REPORTED),
          ],
          ['path', stageBody.path ? h('code', null, stageBody.path) : absentSpan(ABSENT.NOT_REPORTED)],
        ])
      );
      if (!standardColumns.length) {
        pheHost.appendChild(
          message(
            'warn',
            'testing standard not recorded',
            note(
              'The `source` column records where a row came from; it is not the testing standard, ' +
                'and this page does not read it as one. The panel also never converts an MIC into ' +
                'a category and never converts a category into an MIC.'
            )
          )
        );
      }

      // -- per-category counts, one server request each ------------------
      const catHost = h('div');
      catHost.appendChild(h('h4', null, 'Susceptibility categories — one server count each'));
      catHost.appendChild(
        note(
          'Each count below is its own request with `filter.phenotype_sir`, so no category\'s ' +
            'count is derived from another and none is derived from the page. `ND` is a recorded ' +
            'value (the isolate was not determined), not an absent row.'
        )
      );
      const catRows = await Promise.all(
        SIR_CATEGORIES.map(async (category) => {
          try {
            const body = await ctx.api.get('/isolates', { 'filter.phenotype_sir': category, limit: 1 }, abort.signal);
            return { category, total: body.meta ? body.meta.total : null, basis: body.meta ? body.meta.basis : null, unfiltered: body.meta ? body.meta.total_unfiltered : null };
          } catch (error) {
            if (error && error.name === 'AbortError') return null;
            return { category, error: String((error && error.message) || error) };
          }
        })
      );
      if (disposed) return;
      const dl = h('dl', { class: 'kv' });
      for (const entry of catRows) {
        if (!entry) continue;
        dl.appendChild(h('dt', null, entry.category));
        const dd = h('dd');
        if (entry.error) {
          dd.appendChild(absentSpan(`the count could not be fetched: ${entry.error}`, ABSENT.NOT_REPORTED));
        } else {
          dd.appendChild(h('b', null, typeof entry.total === 'number' ? String(entry.total) : ABSENT.NOT_REPORTED));
          dd.appendChild(
            document.createTextNode(
              ` isolate(s) — server count, basis.n = ${entry.basis && entry.basis.n !== undefined ? entry.basis.n : ABSENT.NOT_REPORTED}`
            )
          );
          if (typeof entry.total === 'number') {
            const flag = ctx.banner.power(entry.basis ? entry.basis.n : null, entry.basis);
            if (flag && flag.flag) {
              dd.appendChild(h('span', { class: 'faint' }, ` — ${flag.flag}`));
            }
          }
        }
        dl.appendChild(dd);
      }
      catHost.appendChild(dl);
      catHost.appendChild(
        note(
          'The categories are the vocabulary DESIGN §5.2 declares for the run\'s antibiotic; they ' +
            'are a declared list, not something measured. A category that is absent from the table ' +
            'shows a server count of 0, which is a measurement of the file — not of the cohort\'s ' +
            'biology.'
        )
      );
      pheHost.appendChild(catHost);

      addGrid(pheHost, {
        ctx,
        path: '/stages/phenotype/rows',
        project: false,
        columns: (stageBody.header || []).map((column) => ({ key: column, label: column, sortKey: column })),
        state: { offset: 0, limit: 25, sort: null, q: null, filters: {}, hidden: new Set() },
        search: { label: 'search the phenotype table', placeholder: 'sample id, value, source' },
        exportName: 'phenotype',
        exportLabel: 'rows',
        head: h('h4', null, 'The phenotype table itself (paged by the server)'),
      });
    }
  }

  (async () => {
    // The five panels load independently: one unreadable artefact must not blank
    // the other four, and each reports its own absence.
    await Promise.all([
      renderAmr().catch((error) => error && error.name === 'AbortError' ? null : console.error('[summaries] AMR panel failed', error)),
      renderVirulence().catch((error) => error && error.name === 'AbortError' ? null : console.error('[summaries] virulence panel failed', error)),
      renderVariants().catch((error) => error && error.name === 'AbortError' ? null : console.error('[summaries] variants panel failed', error)),
      renderPangenome().catch((error) => error && error.name === 'AbortError' ? null : console.error('[summaries] pangenome panel failed', error)),
      renderPhenotype().catch((error) => error && error.name === 'AbortError' ? null : console.error('[summaries] phenotype panel failed', error)),
    ]);
    if (disposed) return;
  })().catch((error) => {
    if (error && error.name === 'AbortError') return;
    console.error('[summaries] failed', error);
    alertHost.appendChild(message('error', 'The summaries could not be assembled', note(String((error && error.message) || error))));
  });

  return () => {
    disposed = true;
    abort.abort();
    for (const grid of grids.splice(0)) {
      try {
        grid.destroy();
      } catch (error) {
        console.error('[summaries] a grid teardown threw', error);
      }
    }
    for (const fn of disposers.splice(0)) {
      try {
        fn();
      } catch (error) {
        console.error('[summaries] a disposer threw', error);
      }
    }
  };
}

export default { mount };
