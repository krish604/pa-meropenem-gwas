/* pages/isolates.js — the per-isolate table, 900 rows or 20.
 *
 * Every column of this grid comes from `/api/isolates`, which folds six
 * per-sample tables into one row per `sample_id` and then pages the result
 * server-side. **This file contains no sort, no filter and no count.** It builds
 * a query, hands it to `createPagedGrid`, and renders what comes back. That is
 * the whole trick behind "900 rows instantly": at 900 rows the join is the
 * server's job and the browser holds one page of it.
 *
 * The three column rules that are easy to get wrong, and what this file does:
 *
 * - `[]` in `amr_genes` is a **finding** (the present `04_amr.tsv` holds no row
 *   for this sample), and `not assessed` means the table is not there. The
 *   grid renders the first as a measured zero and the second in the absent
 *   style, and the header panel carries the server's own sentence for both.
 * - `amr_point_mutations` holds rows whose `variant` column is non-empty: the
 *   tool reported the element as **mutated**. It does not mean disruptive, and
 *   the column header says so, because that is the difference between a
 *   determinant call and a claim.
 * - `oprd_state` is `not_assessed` for every isolate unless a verdict source
 *   exists on disk (DESIGN §4.3). This page never fills it in from the master
 *   table's `chromosomal_mutation`, and never from `viz.oprd_status_per_sample`
 *   — on the ten-isolate cohort the latter reported `absent` for 8 isolates
 *   that all carry the gene, which is the study's central negative waiting to
 *   be manufactured. The cell links to the oprD panel, which shows the probes.
 */

import { ABSENT, absentSpan, basisOf, clear, h, kv, message, note, serverMessage } from '../app.js';
import {
  absence,
  cell,
  createPagedGrid,
  isAbsenceWord,
  mountPowerBanner,
} from './tables.js';

/** The four oprD display states (`viz.OPRD_STATE_ORDER`), for the filter select. */
const OPRD_DISPLAY_STATES = ['intact', 'disrupted', 'absent', 'not_assessed'];

/**
 * `oprd_state` is NOT one of the six stage badges, so `ctx.badges.state` cannot
 * be used for `intact`/`disrupted`/`absent`. `not_assessed` IS one of the six
 * and means exactly what the badge says, so it goes through `ctx.badges` like
 * every other state; the other three get a glyph + word chip of their own,
 * because inventing three more badge keys would give two vocabularies the same
 * authority (UI-D2). Colour never carries it alone.
 */
const OPRD_STATE_GLYPH = { intact: '●', disrupted: '✂', absent: '∅', not_assessed: '○' };

function oprdStateCell(row, ctx, oprdReason) {
  const state = row.oprd_state;
  const verdict = row.oprd_verdict;
  const box = h('span', { class: 'chip-row' });
  if (state === 'not_assessed' || state === null || state === undefined) {
    box.appendChild(
      ctx.badges.render(
        'not_assessed',
        oprdReason ||
          'no contracted path holds per-isolate oprD verdicts; the state is not assessed, and it is never `absent` (DESIGN §4.3)'
      )
    );
  } else if (OPRD_DISPLAY_STATES.includes(state)) {
    box.appendChild(
      h(
        'span',
        { class: 'badge', dataset: { tone: state === 'intact' ? 'ok' : state === 'disrupted' ? 'warn' : 'bad' }, title: 'an oprD display state, not one of the six stage badges' },
        h('span', { class: 'badge-glyph', 'aria-hidden': 'true' }, OPRD_STATE_GLYPH[state] || '?'),
        h('span', { class: 'badge-label' }, state.replace('_', ' '))
      )
    );
  } else {
    box.appendChild(absentSpan(`the server sent an oprD display state outside viz.OPRD_STATE_ORDER: ${JSON.stringify(state)}`, 'state not recognised'));
  }
  if (verdict === null || verdict === undefined) {
    box.appendChild(
      h('span', { class: 'faint' }, verdict === null ? 'no verdict recorded' : '')
    );
  } else {
    box.appendChild(h('code', { title: 'the verdict string from adapters/oprd_locus.Verdict, verbatim' }, String(verdict)));
  }
  return box;
}

function sirCell(row, absentBehaviour) {
  const value_ = row.phenotype_sir;
  if (value_ === null || value_ === undefined) {
    return absentSpan(
      'no row for this isolate in 11_phenotype.tsv, or the table is absent',
      ABSENT.NOT_ASSESSED
    );
  }
  if (isAbsenceWord(value_)) return absence(value_, absentBehaviour.phenotype_sir);
  // The code is shown verbatim and NOT translated. A susceptibility category
  // is a clinical interpretation; this page does not invent one, and `ND` is a
  // recorded value (the isolate was not determined), not an absent row.
  return h(
    'span',
    { class: 'chip', dataset: { sir: String(value_) }, title: absentBehaviour.phenotype_sir },
    String(value_)
  );
}

function link(sampleId) {
  const a = h('a', { href: `#/isolates/${encodeURIComponent(sampleId)}` }, sampleId);
  return a;
}

export function mount(container, ctx) {
  const abort = new AbortController();
  let disposed = false;

  container.appendChild(
    h(
      'div',
      { class: 'page-head' },
      h(
        'div',
        { class: 'page-head-text' },
        h('h2', null, 'Isolates'),
        h(
          'p',
          { class: 'page-note' },
          'One row per isolate, joined by the server from the MLST, AMR, ' +
            'virulence, phenotype, master-table and variant-accounting outputs. ' +
            'Paging, sorting, searching and filtering are the server\'s; this page ' +
            'never counts a cohort, so "of N" is always the server\'s N.'
        )
      )
    )
  );

  const alertHost = h('div');
  const summaryHost = h('div', { class: 'panel' });
  const gridHost = h('div');
  container.append(alertHost, summaryHost, gridHost);

  let grid = null;
  let body = null;

  const COLUMNS = [
    {
      key: 'sample_id',
      label: 'isolate',
      sortKey: 'sample_id',
      render: (row) => link(row.sample_id),
      title: 'the one identifier every stage uses; the row exists even when every other cell is absent',
    },
    {
      key: 'st',
      label: 'ST',
      sortKey: 'st',
      title: 'MLST sequence type. `not assessed` when MLST_status is no_call or the row is missing — the stage looked and could not call, which is not the same as never having run.',
    },
    {
      key: 'amr_genes',
      label: 'AMR genes',
      sortKey: 'amr_genes',
      render: (row) => cell(row.amr_genes, 'the AMR table is not there', { title: 'none recorded = the present table holds no row for this isolate' }),
      title: 'determinants from 04_amr.tsv. `none recorded` is a measurement of zero; `not assessed` would mean the table is absent.',
    },
    {
      key: 'amr_point_mutations',
      label: 'AMR point mutations (reported as mutated)',
      sortKey: 'amr_point_mutations',
      render: (row) => cell(row.amr_point_mutations, 'the AMR table is not there', { title: 'rows whose `variant` column is non-empty' }),
      title: 'rows whose `variant` column is non-empty, i.e. the tool reported the element as MUTATED. This is not a statement that the mutation is disruptive — `determinant_type` is `AMR` for a point mutation and for an intact gene alike.',
    },
    {
      key: 'n_virulence',
      label: 'virulence factors',
      sortKey: 'n_virulence',
      render: (row) =>
        isAbsenceWord(row.n_virulence)
          ? absence(row.n_virulence, '08_virulence.tsv is not there')
          : cell(row.n_virulence, 'no value reported'),
      title: 'rows in 08_virulence.tsv for this isolate. `0` is a real measurement when the table is present; `not assessed` means the table is absent.',
    },
    {
      key: 'oprd_state',
      label: 'oprD state / verdict',
      sortKey: 'oprd_state',
      render: (row) => oprdStateCell(row, ctx, body && body.oprd_reason),
      title: 'DESIGN §4.3. Only a verdict read from a real contracted file may say `intact`, `disrupted` or `absent`; every refusal is `not_assessed`.',
    },
    {
      key: 'phenotype_sir',
      label: 'imipenem SIR',
      sortKey: 'phenotype_sir',
      render: (row) => sirCell(row, (body && body.absent_behaviour) || { phenotype_sir: '' }),
      title: 'the run\'s own category for the run\'s antibiotic. Shown verbatim; `ND` is a recorded value and is never turned into an MIC.',
    },
    {
      key: 'lineage',
      label: 'lineage',
      sortKey: 'lineage',
      render: (row) =>
        isAbsenceWord(row.lineage) ? absence(row.lineage, 'neither 15_master_table.tsv nor tree_metadata.tsv gave a lineage for this isolate') : cell(row.lineage),
      title: 'from 15_master_table.tsv, else tree_metadata.tsv. The master table refuses a sentinel lineage, so a sentinel here is surfaced rather than rendered.',
    },
    {
      key: 'n_snv',
      label: 'SNVs',
      sortKey: 'n_snv',
      render: (row) => (isAbsenceWord(row.n_snv) ? absence(row.n_snv, 'variants_provenance.json is not readable') : cell(row.n_snv)),
      title: 'from variants_provenance.json. `0` is a real measurement when that file exists — it writes every key at zero. Absent file reads `not assessed`, because "called nothing" and "never screened" must not look alike.',
    },
    {
      key: 'n_indel',
      label: 'non-SNV variants',
      sortKey: 'n_indel',
      render: (row) => (isAbsenceWord(row.n_indel) ? absence(row.n_indel, 'variants_provenance.json is not readable') : cell(row.n_indel)),
      title: 'the non-SNV classes (indels and anything else the accounting file groups there), from the same source as the SNV column.',
    },
  ];

  function renderSummary(payload) {
    clear(summaryHost);
    summaryHost.appendChild(h('h3', null, 'What this table is made of'));
    const meta = payload.meta || {};
    const basis = basisOf(payload) || {};
    summaryHost.appendChild(
      kv([
        [
          'cohort',
          typeof meta.total === 'number'
            ? h('b', null, `${meta.total} isolate row(s)`)
            : absentSpan('the server reported no total', ABSENT.NOT_REPORTED),
        ],
        [
          'membership',
          payload.membership_source
            ? h('code', { title: 'which artefact defined the cohort: the run manifest, or the union of the per-sample tables' }, payload.membership_source)
            : absentSpan(ABSENT.NOT_REPORTED),
        ],
        ['n used by the flag', basis.n === null || basis.n === undefined ? absentSpan(ABSENT.NOT_REPORTED) : h('b', null, String(basis.n))],
      ])
    );

    const sources = payload.contributing_sources || {};
    const sourceRows = Object.entries(sources).map(([key, present]) => [
      key,
      present
        ? h('span', { class: 'chip', dataset: { flag: 'false' } }, 'present')
        : h('span', { class: 'chip', dataset: { flag: 'true' } }, 'not there'),
    ]);
    if (sourceRows.length) {
      summaryHost.appendChild(h('h4', null, 'Contributing tables'));
      summaryHost.appendChild(
        note(
          'A table marked "not there" makes every cell it feeds read `not assessed`. A table ' +
            'marked "present" that holds no row for an isolate makes that cell a measurement of ' +
            'zero. The two are not the same fact and this grid does not merge them.'
        )
      );
      summaryHost.appendChild(kv(sourceRows));
    }

    const behaviour = payload.absent_behaviour || {};
    const behaviourKeys = Object.keys(behaviour);
    if (behaviourKeys.length) {
      summaryHost.appendChild(h('h4', null, 'What each absent cell means — the server\'s own sentences'));
      summaryHost.appendChild(kv(behaviourKeys.map((key) => [key, h('span', null, behaviour[key])])));
    }

    if (payload.variants_reason) {
      summaryHost.appendChild(h('h4', null, 'Variant accounting'));
      summaryHost.appendChild(serverMessage(payload.variants_reason));
    }
    if (payload.oprd_reason) {
      const box = message('warn', 'The oprD columns are not assessed for any isolate');
      box.appendChild(serverMessage(payload.oprd_reason));
      if (payload.oprd_note) box.appendChild(note(payload.oprd_note));
      box.appendChild(
        h('p', null, h('a', { href: '#/oprd' }, 'open the oprD panel — it lists every probe path that was tried'))
      );
      summaryHost.appendChild(box);
    }
  }

  (async () => {
    // The display-state vocabulary for the oprD filter is the server's, from
    // /api/oprd. If that endpoint cannot be reached the select falls back to
    // the four literals `viz.OPRD_STATE_ORDER` declares, and says so.
    let displayStates = OPRD_DISPLAY_STATES;
    let vocabularyFrom = "viz.OPRD_STATE_ORDER, as DESIGN §4.3 quotes it";
    try {
      const oprd = await ctx.api.get('/oprd', { limit: 1 }, abort.signal);
      if (disposed) return;
      if (Array.isArray(oprd.display_states) && oprd.display_states.length) {
        displayStates = oprd.display_states;
        vocabularyFrom = 'the server (/api/oprd display_states)';
      }
    } catch {
      /* the fallback vocabulary stands, and the label says which was used */
    }
    if (disposed) return;

    clear(gridHost);
    grid = createPagedGrid({
      ctx,
      path: '/isolates',
      columns: COLUMNS,
      project: false,
      state: { offset: 0, limit: 100, sort: 'sample_id', q: null, filters: {}, hidden: new Set() },
      search: {
        label: 'search (server-side substring)',
        placeholder: 'sample id, ST, gene, mutation, SIR, lineage',
      },
      filters: [
        {
          param: 'phenotype_sir',
          label: 'SIR =',
          kind: 'text',
          placeholder: 'R',
          hint: 'exact match on 11_phenotype.tsv for the run\'s antibiotic; `ND` is a recorded value',
        },
        {
          param: 'lineage',
          label: 'lineage =',
          kind: 'text',
          placeholder: 'lineage label',
          hint: 'exact match on the lineage the master table (else the tree metadata) records',
        },
        {
          param: 'oprd_state',
          label: 'oprD state =',
          kind: 'select',
          options: displayStates,
          anyLabel: 'any',
          hint: `the oprD display states, from ${vocabularyFrom}`,
        },
        {
          param: 'amr_genes__in',
          label: 'carries gene =',
          kind: 'text',
          placeholder: 'oprD',
          hint: 'list membership: matches when ANY gene in the cell equals this value, so it means "carries this gene"',
        },
      ],
      exportName: 'isolates',
      exportLabel: 'isolates',
      head: h(
        'div',
        null,
        h('h3', null, 'Isolate table — server-side paging, sorting, filtering'),
        h(
          'p',
          { class: 'panel-note' },
          'The join of six per-sample tables happens on the server, once per request, and ' +
            'only the requested page is sent. The column chooser hides columns here; it ' +
            'does not ask the server for a projection, because /api/isolates builds its ' +
            'projection from the raw `columns` string rather than the parsed column list, ' +
            'and a projected request comes back with empty cells (see the DATA report).'
        )
      ),
      bannerNote:
        'This flag is about the cohort the table was computed on, which for this endpoint ' +
        'is the isolate count. When a filter is applied the n falls to the filtered count, ' +
        'which is what the statistic was actually computed over — that is the point of R16.',
      onLoad: (payload) => {
        body = payload;
        if (payload) renderSummary(payload);
      },
    });
    gridHost.appendChild(grid.el);
    void grid.load();
  })().catch((error) => {
    if (error && error.name === 'AbortError') return;
    console.error('[isolates] failed', error);
    alertHost.appendChild(message('error', 'The isolate table could not be loaded', note(String((error && error.message) || error))));
  });

  return () => {
    disposed = true;
    abort.abort();
    if (grid) grid.destroy();
  };
}

export default { mount };
