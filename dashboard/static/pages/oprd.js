/* pages/oprd.js — the oprD locus, and the hole in it.
 *
 * DESIGN §4.3, and the reason this page has exactly two states.
 *
 * **There is no contracted file holding per-isolate oprD verdicts.**
 * `contracts.py` declares no oprD table in `STAGE_TABLES` or `INTERNAL_TABLES`,
 * `workflow/Snakefile` declares no oprD output, and `run.py`'s in-memory
 * `oprd_locus_resolution` is the only copy the pipeline ever keeps. So
 * `/api/oprd` returns either verdicts read from a real path, or `not produced`
 * with the pipeline's own sentence and every probe path that was tried.
 *
 * **Two fallbacks are refused here, permanently:**
 *
 *  1. the master table's `chromosomal_mutation` column, and
 *  2. `viz.oprd_status_per_sample`.
 *
 * The second is the function the `oprd_locus` adapter's docstring was written to
 * correct: on the ten-isolate cohort it reported `absent` for 8 isolates that
 * all carry the gene, because 34 of the 36 `gene=oprD` CDS features are
 * OprD/OprP/OprQ paralogs. Rendering that as a verdict would manufacture the
 * study's central negative out of a paralog mislabelling. This page has no code
 * path that could do it.
 *
 * The vocabulary is rendered verbatim and never coerced
 * (`adapters/oprd_locus.py:93`):
 *
 *     resolved                          <- the only one asserting the locus is present
 *     refused:no_orthologous_hit
 *     refused:insufficient_coverage
 *     refused:ambiguous_locus
 *     refused:ambiguous_second_locus
 *
 * and the display states are `intact, disrupted, absent, not_assessed`
 * (`viz.OPRD_STATE_ORDER`). Every refusal maps to `not_assessed`; `absent` and
 * `disrupted` come only from a confirmed lesion.
 */

import { ABSENT, absentSpan, clear, h, kv, message, note, probeList, serverMessage } from '../app.js';
import { absence, cell, createPagedGrid, isAbsenceWord, mountPowerBanner } from './tables.js';

/**
 * The focus channel.
 *
 * The router's `parseHash` has no query string, so `#/oprd?sample=X` is "no such
 * page" — a cross-page link cannot carry a parameter through the hash. Both
 * pages are this owner's, so the focus travels through here instead: the
 * isolate page calls `requestFocus`, the oprD page consumes it on mount. A
 * reload loses it, which is why the focus is also shown in the URL-independent
 * state line rather than only used to scroll.
 */
let pendingFocus = null;

export function requestFocus(sampleId) {
  pendingFocus = sampleId ? { sampleId: String(sampleId) } : null;
}

export function takeFocus() {
  const focus = pendingFocus;
  pendingFocus = null;
  return focus;
}

/** Evidence keys this page labels by name rather than dumping as a pair. */
const EVIDENCE_GROUPS = [
  {
    id: 'repeat',
    title: 'repeat / homopolymer evidence',
    test: /repeat|homopolymer|poly[acgt]/i,
    note: 'A repeat or homopolymer at the locus is a candidate mechanism for a slipped-strand misalignment, which is exactly the situation in which a short read cannot resolve the locus.',
  },
  {
    id: 'compensating',
    title: 'compensating-indel evidence',
    test: /compensat/i,
    note: 'A compensating indel elsewhere in the porin network is a separate claim from the oprD locus itself, and is not evidence about this locus.',
  },
  {
    id: 'tblastn',
    title: 'tblastn alignment summary',
    test: /tblastn|alignment|aln|hits|evalue/i,
    note: 'The alignment is what a `resolved` verdict rests on. An alignment against a paralog is not a resolution of the locus; the verdict string already records which happened.',
  },
];

function groupEvidence(evidence) {
  const groups = EVIDENCE_GROUPS.map((g) => ({ ...g, entries: [] }));
  const rest = [];
  for (const [key, value_] of Object.entries(evidence || {})) {
    const home = groups.find((g) => g.test.test(key));
    if (home) home.entries.push([key, value_]);
    else rest.push([key, value_]);
  }
  return { groups: groups.filter((g) => g.entries.length > 0), rest };
}

const VERDICT_NOTES = {
  resolved: 'the only verdict that asserts the locus is present',
  'refused:no_orthologous_hit': 'no orthologous hit: not a resolution, and NOT `absent`',
  'refused:insufficient_coverage': 'insufficient coverage: not a resolution, and NOT `absent`',
  'refused:ambiguous_locus': 'the locus itself is ambiguous: not a resolution, and NOT `absent`',
  'refused:ambiguous_second_locus': 'a second locus is ambiguous: not a resolution, and NOT `absent`',
};

const DISPLAY_STATE_NOTES = {
  intact: 'from `resolved` with no lesion recorded',
  disrupted: 'from `resolved` WITH a confirmed lesion',
  absent: 'only from a confirmed lesion; never from a refusal',
  not_assessed: 'every refusal maps here',
};

/** The reason as the server wrote it, prefixed only when it is not already. */
function reasonText(reason) {
  if (!reason) return `${ABSENT.NOT_PRODUCED}: the server gave no reason`;
  return String(reason).startsWith(ABSENT.NOT_PRODUCED) ? String(reason) : `${ABSENT.NOT_PRODUCED}: ${reason}`;
}

export function mount(container, ctx) {
  const abort = new AbortController();
  let disposed = false;
  const disposers = [];

  let body = null;
  let focus = takeFocus();
  if (!focus && typeof location !== 'undefined' && /sample=/.test(location.hash || '')) {
    // Mounted directly (a harness, or a router that grew a query string).
    const match = /[?&]sample=([^&]+)/.exec(location.hash);
    if (match) focus = { sampleId: decodeURIComponent(match[1]) };
  }

  container.appendChild(
    h(
      'div',
      { class: 'page-head' },
      h(
        'div',
        { class: 'page-head-text' },
        h('h2', null, 'oprD locus'),
        h(
          'p',
          { class: 'page-note' },
          'The oprD locus is the study\'s central question and the pipeline keeps its ' +
            'verdicts in memory only. This page reports what is on disk, which is ' +
            'either a contracted table of verdicts or nothing at all — and says which, ' +
            'with every probe path it tried.'
        )
      ),
      focus
        ? h(
            'div',
            { class: 'page-actions' },
            h('span', { class: 'chip' }, `focused on ${focus.sampleId}`),
            h('a', { href: `#/isolates/${encodeURIComponent(focus.sampleId)}` }, 'open the isolate detail')
          )
        : null
    )
  );

  const alertHost = h('div');
  const statusHost = h('div');
  const vocabHost = h('div', { class: 'panel' });
  const bannerHost = h('div');
  const evidenceHost = h('div');
  const gridHost = h('div');
  container.append(alertHost, statusHost, vocabHost, bannerHost, evidenceHost, gridHost);

  function renderVocabulary(body) {
    clear(vocabHost);
    vocabHost.appendChild(h('h3', null, 'The vocabulary, verbatim'));
    vocabHost.appendChild(
      note(
        'These strings are the pipeline\'s own (`adapters/oprd_locus.Verdict` and ' +
          '`viz.OPRD_STATE_ORDER`). They are reproduced, not paraphrased, and no value ' +
          'outside them is coerced into one.'
      )
    );
    const verdicts = Array.isArray(body.verdict_vocabulary) ? body.verdict_vocabulary : [];
    const states = Array.isArray(body.display_states) ? body.display_states : [];
    vocabHost.appendChild(
      kv(
        verdicts.map((verdict) => [
          h('code', null, verdict),
          h('span', { class: verdict === 'resolved' ? '' : 'muted' }, VERDICT_NOTES[verdict] || 'a value this page has no gloss for'),
        ])
      )
    );
    vocabHost.appendChild(h('h4', null, 'Display states'));
    vocabHost.appendChild(
      kv(states.map((state_) => [h('code', null, state_), h('span', { class: 'muted' }, DISPLAY_STATE_NOTES[state_] || '')]))
    );
  }

  function renderNotProduced(body) {
    clear(statusHost);
    clear(gridHost);
    clear(evidenceHost);
    const box = message('warn', 'oprD verdicts: not produced');
    box.appendChild(serverMessage(reasonText(body.reason)));
    statusHost.appendChild(box);

    statusHost.appendChild(
      message(
        'error',
        'No isolate on this dashboard has an oprD verdict',
        note(
          'Every isolate therefore reads `not_assessed`, which is one of the six badges ' +
            'and is not a claim that the gene is missing. `absent` is not among the ' +
            'available answers here, and this page has no path that could produce it.'
        )
      )
    );

    if (body.not_produced_note) {
      const refused = h('div', { class: 'panel' });
      refused.appendChild(h('h3', null, 'The fallbacks that are refused, and why'));
      refused.appendChild(serverMessage(body.not_produced_note));
      statusHost.appendChild(refused);
    }

    const probes = h('div', { class: 'panel' });
    probes.appendChild(h('h3', null, `Probe paths tried, in order (${(body.probes || []).length})`));
    probes.appendChild(
      note(
        'The endpoint probes for a contracted path first and a bundle-supplied table ' +
          'second. Every path it looked at is listed with the reason it did not answer, ' +
          'because "we looked here and found nothing" and "we never looked" are ' +
          'different facts.'
      )
    );
    probes.appendChild(probeList(body.probes || []));
    statusHost.appendChild(probes);
  }

  function renderFocusRow(body) {
    if (!focus) return;
    const item = (body.items || []).find((entry) => entry.sample_id === focus.sampleId);
    const host = h('div', { class: 'panel' });
    host.appendChild(h('h3', null, `${focus.sampleId} — verdict`));
    if (!body.available) {
      host.appendChild(
        absence(
          'not assessed',
          `no verdict source exists, so ${focus.sampleId} is not assessed — and the master table's chromosomal_mutation column is deliberately not used as a stand-in`
        )
      );
      host.appendChild(note(`The isolate detail page shows the same absence: #/isolates/${encodeURIComponent(focus.sampleId)}`));
      evidenceHost.appendChild(host);
      return;
    }
    if (!item) {
      host.appendChild(
        absentSpan(
          `the verdict source exists but holds no row for ${focus.sampleId}; that is not the ` +
            'same as a refusal, and it is not `absent`',
          ABSENT.NOT_ASSESSED
        )
      );
      evidenceHost.appendChild(host);
      return;
    }
    host.appendChild(verdictBlock(item));
    evidenceHost.appendChild(host);
  }

  function verdictBlock(item) {
    const box = h('div');
    box.appendChild(
      kv([
        ['sample', h('code', null, item.sample_id)],
        [
          'verdict',
          h('code', { title: VERDICT_NOTES[item.verdict] || '' }, String(item.verdict)),
        ],
        [
          'display state',
          isAbsenceWord(item.display_state)
            ? absence(item.display_state, 'the server recorded no display state for this verdict')
            : h('strong', null, String(item.display_state)),
        ],
        ['lesion type', item.lesion_type ? cell(item.lesion_type) : absentSpan('no lesion is recorded; the state above is therefore `intact` only because the verdict is `resolved`')],
        ['position', item.position === null || item.position === undefined ? absentSpan('no coordinate is recorded in the verdict table') : h('span', { class: 'mono' }, String(item.position))],
        ['truncation (aa)', item.truncation_aa === null || item.truncation_aa === undefined ? absentSpan('the verdict table records no truncation length') : h('span', { class: 'mono' }, String(item.truncation_aa))],
        ['identity %', typeof item.identity_pct === 'number' ? h('span', { class: 'mono' }, item.identity_pct.toFixed(2)) : absentSpan('not in the verdict table')],
        ['coverage %', typeof item.coverage_pct === 'number' ? h('span', { class: 'mono' }, item.coverage_pct.toFixed(2)) : absentSpan('not in the verdict table')],
      ])
    );
    const { groups, rest } = groupEvidence(item.evidence);
    for (const group of EVIDENCE_GROUPS) {
      if (groups.includes(group)) continue;
      // Named, because "the endpoint did not carry it" and "there is no such
      // evidence" are different facts. `oprd.py` filters the `evidence` block to
      // its own `EXPECTED_COLUMNS`, so a future verdict table that DOES carry a
      // repeat flag loses it on the way out; that is a gap to report, not
      // something to paper over with a dash.
      box.appendChild(h('h4', null, group.title));
      box.appendChild(
        absentSpan(
          `/api/oprd's evidence block carries no ${group.id} attribute for ` +
            `${item.sample_id}. The block is filtered to \`oprd.EXPECTED_COLUMNS\`, which ` +
            'lists neither a repeat/homopolymer flag nor a compensating-indel flag nor a ' +
            'tblastn summary, so even a verdict table that holds them loses them on the way out.',
          ABSENT.NOT_REPORTED
        )
      );
    }
    if (groups.length === 0 && rest.length === 0) {
      box.appendChild(
        note(
          'The verdict carries no evidence attributes beyond the columns above. A `resolved` ' +
            'verdict with no alignment evidence is a claim the table did not support, and ' +
            'this page shows that absence rather than implying evidence exists.'
        )
      );
      return box;
    }
    for (const group of groups) {
      box.appendChild(h('h4', null, group.title));
      box.appendChild(kv(group.entries.map(([key, value_]) => [h('code', null, key), cell(value_)])));
      box.appendChild(note(group.note));
    }
    if (rest.length) {
      box.appendChild(h('h4', null, 'Other evidence attributes, verbatim'));
      box.appendChild(kv(rest.map(([key, value_]) => [h('code', null, key), cell(value_)])));
    }
    return box;
  }

  function renderAvailable(body) {
    clear(statusHost);
    const intro = h('div', { class: 'panel' });
    intro.appendChild(h('h3', null, 'Verdicts are available'));
    intro.appendChild(
      note(
        `A contracted path answered, so ${(body.meta && body.meta.total) || 0} verdict row(s) ` +
          'are on disk. `absent` and `disrupted` below come only from a confirmed lesion; ' +
          'every `refused:*` row reads `not_assessed`.'
      )
    );
    statusHost.appendChild(intro);
    clear(evidenceHost);
    renderFocusRow(body);
  }

  (async () => {
    try {
      body = await ctx.api.get('/oprd', { limit: 200 }, abort.signal);
    } catch (error) {
      if (error && error.name === 'AbortError') return;
      clear(vocabHost);
      vocabHost.appendChild(
        message('error', '/api/oprd did not answer', note(
          `The oprD verdict page needs that endpoint. ${(error && error.message) || error}`
        ))
      );
      return;
    }
    if (disposed) return;
    clear(alertHost);
    renderVocabulary(body);
    if (!body.available) {
      renderNotProduced(body);
      // No banner: no statistic is displayed. The verdict count is zero
      // because nothing was produced, and "0 verdicts" is exactly the plausible
      // zero this project refuses (D2).
      return;
    }
    renderAvailable(body);

    mountPowerBanner(
      ctx,
      bannerHost,
      body.meta && body.meta.basis,
      'n here is the number of verdict rows the server counted, which is a count of verdicts ' +
        'rather than of isolates. No rate is computed from it on this page.'
    );

    clear(gridHost);
    const grid = createPagedGrid({
      ctx,
      path: '/oprd',
      columns: [
        { key: 'sample_id', label: 'isolate', sortKey: 'sample_id', render: (row) => h('a', { href: `#/isolates/${encodeURIComponent(row.sample_id)}` }, row.sample_id) },
        { key: 'verdict', label: 'verdict (verbatim)', sortKey: 'verdict', title: VERDICT_NOTES.resolved },
        {
          key: 'display_state',
          label: 'display state',
          sortKey: 'display_state',
          render: (row) =>
            isAbsenceWord(row.display_state) ? absence(row.display_state, 'no display state recorded') : h('strong', null, row.display_state),
          title: 'a refusal is `not_assessed`; `absent` needs a confirmed lesion',
        },
        { key: 'lesion_type', label: 'lesion type', sortKey: 'lesion_type' },
        { key: 'position', label: 'position', sortKey: 'position' },
        { key: 'truncation_aa', label: 'truncation (aa)', sortKey: 'truncation_aa' },
        { key: 'identity_pct', label: 'identity %', sortKey: 'identity_pct' },
        { key: 'coverage_pct', label: 'coverage %', sortKey: 'coverage_pct' },
        {
          key: 'evidence',
          label: 'evidence attributes',
          sort: false,
          render: (row) => {
            const { groups, rest } = groupEvidence(row.evidence);
            if (!groups.length && !rest.length) return absentSpan('the verdict row carries no evidence attribute beyond its own columns');
            return h(
              'div',
              { class: 'chip-row' },
              groups.map((g) =>
                h('span', { class: 'chip', title: g.note }, `${g.title}: ${g.entries.map(([k, v]) => `${k}=${v}`).join(', ')}`)
              )
            );
          },
        },
      ],
      project: false,
      state: { offset: 0, limit: 100, sort: 'sample_id', q: null, filters: {}, hidden: new Set() },
      search: { label: 'search (verdict table)', placeholder: 'sample id or verdict' },
      filters: [
        { param: 'verdict', label: 'verdict =', kind: 'text', placeholder: 'resolved', hint: 'exact match on the verbatim verdict string' },
        { param: 'display_state', label: 'display state =', kind: 'select', options: body.display_states || [], anyLabel: 'any' },
      ],
      exportName: 'oprd-verdicts',
      exportLabel: 'verdicts',
      head: h(
        'div',
        null,
        h('h3', null, 'Verdict table'),
        h('p', { class: 'panel-note' }, 'Paged and filtered by the server. A `refused:*` verdict is shown as its own string; it is never folded into `absent`.')
      ),
    });
    gridHost.appendChild(grid.el);
    void grid.load();
    disposers.push(() => grid.destroy());
  })().catch((error) => {
    if (error && error.name === 'AbortError') return;
    console.error('[oprd] failed', error);
  });

  return () => {
    disposed = true;
    abort.abort();
    for (const fn of disposers.splice(0)) {
      try {
        fn();
      } catch (error) {
        console.error('[oprd] a disposer threw', error);
      }
    }
  };
}

export default { mount };
