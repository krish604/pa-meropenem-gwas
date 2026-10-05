/* pages/isolate.js — every stage's result for ONE isolate, in one place.
 *
 * Data sources, all server-side, all abortable:
 *
 * - `/api/isolates/{sample_id}` — the joined row, the ten contributing tables
 *   (present or not, with their reasons), the variant accounting and the oprD
 *   block.
 * - `/api/stages` — the 16 stage badges, so a stage that never ran reads
 *   `not run` with the run's own reason and is visibly different from a stage
 *   that ran and whose output cannot answer.
 * - `/api/stages/{stage}/rows?filter.sample_id=…` — this isolate's own rows.
 *   The filter is only *applied* when `sample_id` is the table's index key
 *   column; otherwise the server reports it in `meta.read.unsupported_filter_columns`
 *   and the block says the table is cohort-level instead of pretending it holds
 *   a row for this isolate.
 * - `/api/tree/tips?q=…` — the tip's `node_id`, so the reader can go to the
 *   tree page and find the same node (UI-D6: the node id is a postorder index
 *   in one file, not a biological identity, and the block says so).
 *
 * Wording rules, which are the whole point of this page:
 *
 * - a table that is **not there** → `not produced: <the server's reason>`;
 * - a table that is **there and holds no row for this isolate** → "the present
 *   table records nothing for this isolate", which is a measurement of zero;
 * - a stage that **never ran** → the six-badge `not run` with its reason;
 * - a stage that **ran and cannot answer** → `not assessed` with its reason.
 *
 * None of those is a blank, a `0`, a dash, or each other.
 */

import { ABSENT, absentSpan, basisOf, clear, h, kv, message, note, serverMessage } from '../app.js';
import { absence, cell, isAbsenceWord, mountPowerBanner } from './tables.js';
import { requestFocus } from './oprd.js';

const ROWS_PER_STAGE = 25;
const COHORT_HEAD_ROWS = 5;

/** Which of the joined columns each contributing table feeds. */
const COLUMN_SOURCE = {
  st: 'mlst',
  amr_genes: 'amr',
  amr_point_mutations: 'amr',
  n_virulence: 'virulence',
  phenotype_sir: 'phenotype',
  lineage: 'master_table',
};

function tableStatus(tables, key) {
  const found = (tables || []).find((entry) => entry.key === key);
  return found || { key, present: false, reason: `the server returned no record for ${key}`, rows: [] };
}

/**
 * Why one cell of the joined row reads the way it does.
 *
 * The point is that `0` and `not assessed` are told apart by the *table's*
 * state, not by the cell: the same value means different things depending on
 * whether the artefact exists.
 */
function cellReason(row, column, tables, extra) {
  const key = COLUMN_SOURCE[column];
  const status = key ? tableStatus(tables, key) : null;
  if (status && !status.present) return status.reason;
  switch (column) {
    case 'st':
      return 'MLST_status is `no_call` for this isolate, or 03_mlst.tsv holds no row for it. The stage looked and could not call — which is not `not run` and not `absent`.';
    case 'amr_genes':
    case 'amr_point_mutations':
      return '04_amr.tsv is present and holds no row for this isolate. An empty list here is a measurement of zero, not a missing value.';
    case 'n_virulence':
      return '08_virulence.tsv is present and holds no row for this isolate. `0` is a real count.';
    case 'phenotype_sir':
      return '11_phenotype.tsv is present and holds no row for this isolate, so no susceptibility call was recorded. `ND` would be a recorded value and is shown as one.';
    case 'lineage':
      return 'neither 15_master_table.tsv nor tree_metadata.tsv records a lineage for this isolate';
    case 'n_snv':
    case 'n_indel':
      return (extra && extra.variantsReason) || 'variants_provenance.json is not readable, so the variant accounting is not assessed. `not assessed` here is deliberately not `0`: "called nothing" and "never screened" are different facts.';
    default:
      return 'the server sent no value and no reason for this cell';
  }
}

function joinedCell(row, column, tables, extra) {
  const value_ = row[column];
  const reason = cellReason(row, column, tables, extra);
  if (value_ === null || value_ === undefined) {
    // The artefact behind this column decides the word: a table that is not
    // there is `not produced`; a table that is there and holds no row is
    // `not assessed` for this cell with the zero already shown elsewhere.
    const key = COLUMN_SOURCE[column];
    const status = key ? tableStatus(tables, key) : null;
    return absentSpan(reason, status && !status.present ? ABSENT.NOT_PRODUCED : ABSENT.NOT_ASSESSED);
  }
  if (isAbsenceWord(value_)) return absence(value_, reason);
  return cell(value_, reason);
}

function miniTable(rows, columns) {
  if (!Array.isArray(rows) || rows.length === 0) return null;
  const table = h('table', { class: 'data-table' });
  const thead = h('thead');
  const hr = h('tr');
  for (const column of columns) hr.appendChild(h('th', { scope: 'col' }, column));
  thead.appendChild(hr);
  const tbody = h('tbody');
  for (const row of rows) {
    const tr = h('tr');
    for (const column of columns) {
      const td = h('td');
      const value_ = row[column];
      if (value_ === null || value_ === undefined || value_ === '') {
        td.appendChild(absentSpan(`the \`${column}\` cell of this row holds no value`));
      } else if (Array.isArray(value_)) {
        td.appendChild(cell(value_));
      } else {
        td.textContent = String(value_);
      }
      tr.appendChild(td);
    }
    tbody.appendChild(tr);
  }
  table.append(thead, tbody);
  return h('div', { class: 'table-wrap' }, table);
}

export function mount(container, ctx) {
  const abort = new AbortController();
  let disposed = false;
  const disposers = [];

  const current = ctx.router.current();
  const sampleId = current && current.params ? decodeURIComponent(current.params.sample_id || '') : '';

  const headHost = h('div', { class: 'page-head' });
  const alertHost = h('div');
  const summaryHost = h('div', { class: 'panel' });
  const oprdHost = h('div', { class: 'panel' });
  const treeHost = h('div', { class: 'panel' });
  const stagesHost = h('div');
  container.append(headHost, alertHost, summaryHost, oprdHost, treeHost, stagesHost);

  if (!sampleId) {
    headHost.appendChild(
      h(
        'div',
        { class: 'page-head-text' },
        h('h2', null, 'Isolate detail'),
        message(
          'warn',
          'No isolate was named in the address',
          note('This page is reached as `#/isolates/<sample_id>`. The isolate table links every row here.')
        )
      )
    );
    return () => {
      abort.abort();
      for (const fn of disposers.splice(0)) fn();
    };
  }

  headHost.appendChild(
    h(
      'div',
      { class: 'page-head-text' },
      h('h2', null, `Isolate ${sampleId}`),
      h(
        'p',
        { class: 'page-note' },
        'Every stage that has something to say about this isolate, and a named ' +
          'reason for every stage that does not.'
      )
    )
  );
  headHost.appendChild(
    h('div', { class: 'page-actions' }, h('a', { href: '#/isolates' }, '← back to the isolate table'))
  );

  /* -- the joined row ------------------------------------------------- */

  function renderSummary(detail) {
    clear(summaryHost);
    summaryHost.appendChild(h('h3', null, 'The joined row'));
    const tables = detail.tables || [];
    const extra = { variantsReason: detail.variants_reason };
    const present = tables.filter((t) => t.present).map((t) => t.key);
    const absent = tables.filter((t) => !t.present);
    summaryHost.appendChild(
      note(
        `${present.length} contributing table(s) answered (${present.join(', ')})` +
          (absent.length ? `; ${absent.length} did not: ${absent.map((t) => t.key).join(', ')}` : '') +
          '. A table that did not answer makes its cells `not produced`; a table that ' +
          'answered and holds no row for this isolate makes them a measurement of zero. ' +
          'The two are never merged.'
      )
    );
    const pairs = [
      ['isolate', h('code', null, detail.sample_id)],
      ['cohort membership', detail.membership_source ? h('code', null, detail.membership_source) : absentSpan(ABSENT.NOT_REPORTED)],
    ];
    for (const column of ['st', 'amr_genes', 'amr_point_mutations', 'n_virulence', 'phenotype_sir', 'lineage', 'n_snv', 'n_indel']) {
      const label =
        column === 'amr_point_mutations'
          ? 'AMR point mutations (reported as mutated)'
          : column === 'n_indel'
            ? 'non-SNV variants'
            : column === 'phenotype_sir'
              ? 'imipenem SIR'
              : column.replace(/^n_/, '').replace(/_/g, ' ');
      pairs.push([h('span', { title: cellReason(detail.row, column, tables, extra) }, label), joinedCell(detail.row, column, tables, extra)]);
    }
    summaryHost.appendChild(kv(pairs));
    if (absent.length) {
      summaryHost.appendChild(h('h4', null, 'Contributing tables that are not there, with the server\'s reason'));
      summaryHost.appendChild(
        kv(absent.map((t) => [h('code', null, t.key), h('span', null, t.reason || ABSENT.NOT_PRODUCED)]))
      );
    }
    if (detail.variants_reason) {
      summaryHost.appendChild(h('h4', null, 'Variant accounting'));
      summaryHost.appendChild(serverMessage(detail.variants_reason));
    }
    const provenance = detail.variants_provenance;
    summaryHost.appendChild(
      h('h4', null, 'variants_provenance.json entry for this isolate')
    );
    if (!provenance) {
      summaryHost.appendChild(
        absentSpan(
          (detail.variants_reason || ABSENT.NOT_PRODUCED) +
            '. Not `0` SNVs: an absent accounting file means the screening is not assessed.',
          ABSENT.NOT_ASSESSED
        )
      );
    } else {
      summaryHost.appendChild(
        kv(Object.entries(provenance).map(([key, value_]) => [h('code', null, key), cell(value_)]))
      );
    }
  }

  /* -- oprD ----------------------------------------------------------- */

  function renderOprd(detail) {
    clear(oprdHost);
    const block = detail.oprd || {};
    oprdHost.appendChild(h('h3', null, 'oprD locus'));
    if (block.available && (block.items || []).length) {
      oprdHost.appendChild(kv([['verdict', miniTable(block.items, ['verdict', 'lesion_type', 'position', 'truncation_aa', 'display_state'])]]));
    } else {
      oprdHost.appendChild(
        absentSpan(
          block.reason ||
            'no contracted path holds per-isolate oprD verdicts, so this isolate is not assessed',
          ABSENT.NOT_PRODUCED
        )
      );
    }
    // Stated in both states, because it is a standing refusal rather than a
    // statement about this run: the two fallbacks below are never used, whatever
    // this isolate's verdict happens to be.
    oprdHost.appendChild(
      note(
        'The master table\'s `chromosomal_mutation` column and `viz.oprd_status_per_sample` are ' +
          'deliberately NOT used to fill this in, in either state. The latter reported `absent` ' +
          'for 8 isolates that all carry the gene, because 34 of the 36 `gene=oprD` CDS features ' +
          'are OprD/OprP/OprQ paralogs — rendering that here would manufacture the study\'s ' +
          'central negative (DESIGN §4.3).'
      )
    );
    if (Array.isArray(block.probes) && block.probes.length) {
      const details = h('details');
      details.appendChild(h('summary', null, `probe paths tried (${block.probes.length})`));
      const list = h('ol', { class: 'probe-list' });
      for (const probe of block.probes) {
        const item = h('li');
        item.appendChild(h('span', { class: 'probe-n' }, String(probe.n)));
        const body = h('div');
        body.appendChild(h('span', { class: probe.ok ? 'probe-ok' : 'probe-miss' }, probe.ok ? 'matched' : probe.reason || 'no match'));
        body.appendChild(h('div', { class: 'probe-path' }, probe.path));
        item.appendChild(body);
        list.appendChild(item);
      }
      details.appendChild(list);
      oprdHost.appendChild(details);
    }
    if (block.available) {
      const button = h('button', { type: 'button' }, 'open this isolate in the oprD panel');
      const onClick = () => {
        requestFocus(detail.sample_id);
        void ctx.router.go('/oprd');
      };
      button.addEventListener('click', onClick);
      disposers.push(() => button.removeEventListener('click', onClick));
      oprdHost.appendChild(button);
    }
  }

  /* -- tree position --------------------------------------------------- */

  function renderTree(tips, sample) {
    clear(treeHost);
    treeHost.appendChild(h('h3', null, 'Position in the stage-9 tree'));
    if (!tips || tips.present === false) {
      treeHost.appendChild(absentSpan((tips && tips.reason) || 'no tree is available to place this isolate in', ABSENT.NOT_PRODUCED));
      return;
    }
    const tip = (tips.items || []).find((entry) => entry.sample_id === sample) || null;
    if (!tip) {
      treeHost.appendChild(
        absentSpan(
          'this isolate is not a tip of the stage-9 tree. Either it is absent from the tree, or the tip labels do not join to the sample ids — which is a data problem worth reporting, not an absence to render as one.',
          ABSENT.NOT_ASSESSED
        )
      );
      return;
    }
    treeHost.appendChild(
      kv([
        ['node_id', h('code', null, String(tip.node_id))],
        ['tip label', h('code', null, tip.label || ABSENT.NOT_REPORTED)],
        ['lineage label (tree metadata)', tip.lineage_label ? cell(tip.lineage_label) : absentSpan('tree_metadata.tsv records no lineage for this tip')],
        ['ST (tree metadata)', tip.st ? cell(tip.st) : absentSpan('tree_metadata.tsv records no ST for this tip')],
        ['metadata source', tip.source ? h('code', null, tip.source) : absentSpan(ABSENT.NOT_REPORTED)],
        ['identity', h('span', null, tips.identity || 'postorder index from the root, assigned server-side')],
      ])
    );
    treeHost.appendChild(
      note(
        `A node_id is this file's postorder index, not a biological identity and not a clade ` +
          'number (UI-D6). The tree page owns the drawing; this page links to it and names the ' +
          'node, because the router\'s hash has no query string to carry a focus through.'
      )
    );
    // A machine-readable breadcrumb for the tree page, without depending on a
    // router feature that does not exist yet.
    globalThis.__paDashboard = Object.assign(globalThis.__paDashboard || {}, {
      links: Object.assign((globalThis.__paDashboard || {}).links || {}, {
        tree_node: tip.node_id,
        tree_sample: sample,
      }),
    });
    const link = h('a', { href: '#/tree' }, 'open the tree');
    treeHost.appendChild(link);
    treeHost.appendChild(
      h('span', { class: 'faint' }, ` and look for node ${tip.node_id} / label ${tip.label || sample}`)
    );
  }

  /* -- per-stage blocks ------------------------------------------------ */

  function renderStageBlock(stage, body) {
    const box = h('section', { class: 'panel', dataset: { stage: stage.name } });
    const head = h('div', { class: 'panel-head' });
    head.appendChild(h('h3', null, `${stage.order}. ${stage.name}`));
    head.appendChild(h('span', { class: 'stage-spec' }, `${stage.spec_name} — spec.md D1 #${stage.spec_number}`));
    head.appendChild(ctx.badges.render(stage, stage.reason));
    box.appendChild(head);

    if (!body) {
      box.appendChild(note('this stage\'s table could not be read: the endpoint did not answer'));
      return box;
    }
    if (body.present === false) {
      box.appendChild(absentSpan(body.reason || 'the server gave no reason', ABSENT.NOT_PRODUCED));
      box.appendChild(
        note(
          'The stage is therefore `not assessed` rather than `not run`: the run recorded it, and ' +
            'its declared output is not there to answer from. The badge above carries the run\'s own reason.'
        )
      );
      return box;
    }

    const meta = body.meta || {};
    const unsupported = (meta.read && meta.read.unsupported_filter_columns) || [];
    const header = Array.isArray(body.header) ? body.header : [];
    const items = Array.isArray(body.items) ? body.items : [];
    const keyed = !unsupported.length && header.includes('sample_id');

    if (!keyed) {
      box.appendChild(
        h(
          'p',
          null,
          h('strong', null, 'cohort-level output: this table has no per-isolate row.'),
          ' ',
          unsupported.length
            ? absentSpan(`the server refused a sample_id filter on this table (${unsupported.join(', ')}) — its byte-offset index is keyed on another column`, ABSENT.NOT_ASSESSED)
            : absentSpan(`the header carries no sample_id column (${header.join(', ') || 'no header'})`, ABSENT.NOT_ASSESSED)
        )
      );
      box.appendChild(
        kv([
          ['rows in the table (server count)', typeof meta.total === 'number' ? h('b', null, String(meta.total)) : absentSpan(ABSENT.NOT_REPORTED)],
          ['path', body.path ? h('code', null, body.path) : absentSpan(ABSENT.NOT_REPORTED)],
        ])
      );
      box.appendChild(note(`The first ${Math.min(COHORT_HEAD_ROWS, items.length)} row(s) of the table, for context — not this isolate's row:`));
      const shown = items.slice(0, COHORT_HEAD_ROWS);
      const table = miniTable(shown, header.slice(0, 10));
      if (table) box.appendChild(table);
      else box.appendChild(absentSpan('the table carries no row to show', ABSENT.NOT_PRODUCED));
      if (header.length > 10) {
        box.appendChild(note(`${header.length - 10} further column(s) not shown here; the tables browser pages the whole header.`));
      }
      return box;
    }

    box.appendChild(
      kv([
        [
          `rows for ${sampleId} (server count)`,
          typeof meta.total === 'number'
            ? meta.total === 0
              ? h('b', null, '0')
              : h('b', null, String(meta.total))
            : absentSpan(ABSENT.NOT_REPORTED),
        ],
        ['of rows in the whole table', typeof meta.total_unfiltered === 'number' ? String(meta.total_unfiltered) : absentSpan('no filter was applied, so there is no separate unfiltered total')],
        ['path', body.path ? h('code', null, body.path) : absentSpan(ABSENT.NOT_REPORTED)],
      ])
    );
    if (meta.total === 0) {
      box.appendChild(
        note(
          'The table is present and holds no row for this isolate. That is a measurement of zero, ' +
            'and it is not the same as the table being absent (which would read `not produced`).'
        )
      );
      return box;
    }
    const table = miniTable(items, header);
    if (table) box.appendChild(table);
    if (meta.total > items.length) {
      box.appendChild(
        note(
          `Showing the first ${items.length} of ${meta.total} row(s) for this isolate, capped by this page ` +
            `at ${ROWS_PER_STAGE}; the stage-tables browser pages the rest.`
        )
      );
    }
    return box;
  }

  (async () => {
    let detail = null;
    try {
      detail = await ctx.api.get(`/isolates/${encodeURIComponent(sampleId)}`, null, abort.signal);
    } catch (error) {
      if (error && error.name === 'AbortError') return;
      clear(alertHost);
      const box = message('error', `The isolate detail for ${sampleId} could not be read`);
      if (error && error.serverMessage) box.appendChild(serverMessage(error.serverMessage));
      else box.appendChild(note(String((error && error.message) || error)));
      box.appendChild(h('p', null, h('a', { href: '#/isolates' }, 'back to the isolate table')));
      alertHost.appendChild(box);
      return;
    }
    if (disposed) return;

    renderSummary(detail);
    renderOprd(detail);

    const basis = basisOf(detail);
    if (basis) {
      const bannerHost = h('div');
      container.insertBefore(bannerHost, stagesHost);
      mountPowerBanner(
        ctx,
        bannerHost,
        basis,
        'This page shows one isolate. Its basis n is 1 because the statistic — if any — would be ' +
          'about one isolate; the cohort figures live on the isolate table, where n is the cohort.'
      );
    }

    // The tree position and the stage list, in parallel; both are abortable and
    // neither is allowed to fail the page.
    const [stagesBody, tipsBody] = await Promise.all([
      ctx.api.get('/stages', null, abort.signal).catch((error) => (error && error.name === 'AbortError' ? null : error)),
      ctx.api.get('/tree/tips', { q: sampleId, limit: 5 }, abort.signal).catch((error) => (error && error.name === 'AbortError' ? null : error)),
    ]);
    if (disposed) return;

    renderTree(tipsBody, sampleId);

    if (!stagesBody || !Array.isArray(stagesBody.items)) {
      stagesHost.appendChild(
        message(
          'warn',
          'The stage list could not be read',
          note(
            `Without /api/stages the badges below cannot be shown, and this page will not invent ` +
              `them. ${stagesBody && stagesBody.message ? stagesBody.message : ''}`
          )
        )
      );
      return;
    }

    const stageItems = stagesBody.items;
    const results = await Promise.all(
      stageItems.map((stage) =>
        ctx.api
          .get(`/stages/${encodeURIComponent(stage.name)}/rows`, { 'filter.sample_id': sampleId, limit: ROWS_PER_STAGE }, abort.signal)
          .catch((error) => (error && error.name === 'AbortError' ? null : { __error: error, present: false, reason: `the endpoint did not answer: ${error && error.message}` }))
      )
    );
    if (disposed) return;

    stagesHost.appendChild(
      h(
        'div',
        { class: 'page-head' },
        h(
          'div',
          { class: 'page-head-text' },
          h('h3', null, `Stage by stage — ${stageItems.length} stages, ${sampleId}`),
          h(
            'p',
            { class: 'page-note' },
            'Each block is this isolate\'s own rows from that stage\'s table, fetched with a ' +
              'server-side `filter.sample_id`. Where the table has no per-isolate row the block ' +
              'says so instead of showing a cohort-level row as if it were this isolate\'s.'
          )
        )
      )
    );
    stageItems.forEach((stage, index) => {
      const body = results[index];
      stagesHost.appendChild(
        body && body.__error
          ? renderStageBlock(stage, { present: false, reason: body.reason })
          : renderStageBlock(stage, body)
      );
    });
  })().catch((error) => {
    if (error && error.name === 'AbortError') return;
    console.error('[isolate] failed', error);
  });

  return () => {
    disposed = true;
    abort.abort();
    for (const fn of disposers.splice(0)) {
      try {
        fn();
      } catch (error) {
        console.error('[isolate] a disposer threw', error);
      }
    }
  };
}

export default { mount };
