/* pages/provenance.js — what ran, from what, and what is not recorded.
 *
 * The honesty rule this file is built around:
 *
 * > **`bakta_executions` is in NEITHER manifest writer** (DESIGN assumption A6,
 * > `papipeline/run.py:2464` and `scripts/common/write_provenance.py:101`).
 * > When the field is ABSENT this page renders `not reported`. It never renders
 * > `0`. `0` says "Bakta was never invoked", which is a claim of proof, and
 * > `provenance.bakta_executions` returns `None` — not `0` — precisely because
 * > that claim is not available.
 *
 * The rule is not a special case in one widget: it is why every field on this
 * page goes through `ctx.badges.notReported` / `absentSpan`, so a missing value
 * cannot be rendered as a blank, a zero or an empty string (UI-D2).
 */

import {
  ABSENT,
  absentSpan,
  clear,
  dataTable,
  h,
  kv,
  message,
  note,
  probeList,
} from '../app.js';

export function mount(container, ctx) {
  const abort = new AbortController();
  let disposed = false;

  let provenance = null;
  let bakta = null;
  let digests = null;
  let logs = null;

  container.appendChild(
    h(
      'div',
      { class: 'page-head' },
      h(
        'div',
        { class: 'page-head-text' },
        h('h2', null, 'Provenance'),
        h(
          'p',
          { class: 'page-note' },
          'Tool versions, the annotation-reuse record, the "Bakta executions" ' +
            'proof, a SHA-256 per stage artefact, the git state of the opened ' +
            'root and the configuration that was loaded — and, for each of ' +
            'them, whether it was recorded at all.'
        )
      )
    )
  );

  const errorHost = h('div');
  const proofHost = h('div');
  const toolsHost = h('div', { class: 'panel' });
  const referencesHost = h('div', { class: 'panel' });
  const reuseHost = h('div', { class: 'panel' });
  const logsHost = h('div', { class: 'panel' });
  const digestsHost = h('div', { class: 'panel' });
  const gitHost = h('div', { class: 'panel' });
  const configHost = h('div', { class: 'panel' });
  container.append(errorHost, proofHost, toolsHost, referencesHost, reuseHost, logsHost, digestsHost, gitHost, configHost);

  /* -- the proof panel ------------------------------------------------ */

  /**
   * "Bakta executions: N", or `not reported`.
   *
   * The wording is the server's (`bakta_executions_label`), because the server
   * owns the template; this decides how it is *presented*: a reported count is
   * bold with its derivation, an unreported one is set in the absent style with
   * the reason attached, and there is no code path here that can produce a
   * bare `0`.
   */
  function renderProof() {
    clear(proofHost);
    if (bakta === null) {
      proofHost.appendChild(
        message('info', 'Proof', absentSpan('/api/provenance/bakta has not answered', ABSENT.NOT_REPORTED))
      );
      return;
    }
    const box = h('div', { class: 'proof-panel' });
    box.appendChild(h('h3', null, 'Annotation reuse — the proof panel'));

    const reported = typeof bakta.bakta_executions === 'number';
    const label = bakta.bakta_executions_label || (reported ? `Bakta executions: ${bakta.bakta_executions}` : 'not reported');

    box.appendChild(
      h('div', { class: 'proof-value', dataset: { reported: String(reported) } }, label)
    );
    box.appendChild(
      note(
        reported
          ? `The server derived this from the annotation-reuse record: ` +
              `${rowsOrUnknown(bakta.rows)} row(s), ` +
              `${countOrUnknown(bakta.n_reused)} reused, ` +
              `${countOrUnknown(bakta.n_not_reused)} not reused. Every non-reused sample ran the tool once.`
          : '**Neither manifest writer records this field.** ' +
              '`bakta_executions` appears in neither `papipeline/run.py:_write_run_manifest` ' +
              'nor `scripts/common/write_provenance.py`, and the reuse record on this page ' +
              'is not complete enough to mean a count (it needs one row per sample). ' +
              'So the count is reported absent. It is deliberately NOT shown as 0: ' +
              'a zero would assert that Bakta was never invoked, and nothing here can ' +
              'support that claim.'
      )
    );

    box.appendChild(
      kv([
        [
          'reuse record on disk',
          bakta.available
            ? h('strong', null, 'present')
            : absentSpan(bakta.reason || 'no reuse record was found'),
        ],
        ['reused from run', bakta.reused_from ? h('code', null, bakta.reused_from) : absentSpan('the record names no single source directory')],
        [
          'n reused',
          typeof bakta.n_reused === 'number' ? String(bakta.n_reused) : absentSpan('the record holds no rows, so no count of reused samples is a measurement'),
        ],
        [
          'n not reused',
          typeof bakta.n_not_reused === 'number' ? String(bakta.n_not_reused) : absentSpan('as above'),
        ],
        [
          'cohort',
          typeof bakta.n_samples === 'number' ? String(bakta.n_samples) : absentSpan('no manifest key records a cohort size'),
        ],
      ])
    );

    const rows = Array.isArray(bakta.rows) ? bakta.rows : [];
    if (rows.length > 0) {
      const columns = Object.keys(rows[0]);
      box.appendChild(h('h4', null, 'The reuse record, row by row'));
      box.appendChild(dataTable(columns, rows.slice(0, 200), 'the record holds no rows'));
      if (rows.length > 200) {
        box.appendChild(note(`Showing the first 200 of ${rows.length} rows the server returned.`, 'table-meta'));
      }
    }

    proofHost.appendChild(box);
  }

  /* -- tools ---------------------------------------------------------- */

  function renderTools() {
    clear(toolsHost);
    toolsHost.appendChild(
      h('div', { class: 'panel-head' }, h('h3', null, 'Tool versions'))
    );
    if (!provenance) {
      toolsHost.appendChild(absentSpan('/api/provenance has not answered', ABSENT.NOT_REPORTED));
      return;
    }
    const source = provenance.tool_source;
    toolsHost.appendChild(
      note(
        source === 'tools_detected'
          ? 'From `run_manifest.tools_detected`, written after the run by `papipeline/run.py`.'
          : source === 'stage_tool_requirements'
            ? 'From `run_manifest.stage_tool_requirements`, written before the run by `scripts/common/write_provenance.py`. The per-stage requirements are flattened to one row per tool, with the stages that declared it.'
            : 'No manifest key carries a tool matrix: the writer produced neither `tools_detected` nor `stage_tool_requirements`, or no manifest is readable. This is `not reported`, not an empty matrix meaning "no tools were needed".'
      )
    );
    const tools = provenance.tools && typeof provenance.tools === 'object' ? provenance.tools : {};
    const names = Object.keys(tools);
    if (names.length === 0) {
      toolsHost.appendChild(
        absentSpan(
          source === 'none'
            ? 'no manifest records any tool version or availability'
            : 'the manifest key is present but carries no entries'
        )
      );
      return;
    }
    const rows = names.sort().map((name) => {
      const info = tools[name] || {};
      return {
        tool: name,
        version: info.version === null || info.version === undefined ? '' : String(info.version),
        available:
          info.available === null || info.available === undefined ? '' : String(info.available),
        executable: info.executable === null || info.executable === undefined ? '' : String(info.executable),
        stages: Array.isArray(info.stages) ? info.stages.join(', ') : '',
      };
    });
    const table = h('div', { class: 'table-wrap' });
    const el = h('table', { class: 'tool-matrix' });
    el.appendChild(headRow(['tool', 'version', 'available', 'executable', 'declared by']));
    const tbody = h('tbody');
    for (const row of rows) {
      const tr = h('tr');
      tr.appendChild(h('td', { class: 'tool-name' }, row.tool));
      for (const key of ['version', 'available', 'executable', 'stages']) {
        const cell = h('td');
        if (row[key] === '') {
          cell.appendChild(
            absentSpan(
              key === 'version'
                ? 'the manifest records no version for this tool'
                : key === 'stages'
                  ? 'the run writer does not record which stage declared a tool'
                  : `the ${key === 'available' ? 'availability' : key} was not recorded`
            )
          );
        } else {
          cell.textContent = row[key];
        }
        tr.appendChild(cell);
      }
      tbody.appendChild(tr);
    }
    el.appendChild(tbody);
    table.appendChild(el);
    toolsHost.appendChild(table);
  }

  /* -- references ----------------------------------------------------- */

  function renderReferences() {
    clear(referencesHost);
    referencesHost.appendChild(h('h3', null, 'References'));
    if (!provenance) {
      referencesHost.appendChild(absentSpan('/api/provenance has not answered', ABSENT.NOT_REPORTED));
      return;
    }
    const references = Array.isArray(provenance.references) ? provenance.references : null;
    if (references === null) {
      referencesHost.appendChild(absentSpan('the manifest carries no `references` key'));
    } else if (references.length === 0) {
      referencesHost.appendChild(
        absentSpan(
          'the manifest carries `references` with nothing in it: the run recorded no reference row, which is a fact about the run and not an error'
        )
      );
    } else {
      const columns = ['reference_id', 'tool', 'tool_version', 'database', 'database_version', 'version_status'];
      const table = h('div', { class: 'table-wrap' });
      const el = h('table');
      el.appendChild(headRow(columns));
      const tbody = h('tbody');
      for (const row of references) {
        const tr = h('tr');
        for (const column of columns) {
          const cell = h('td');
          const v = row[column];
          if (v === null || v === undefined || v === '') {
            cell.appendChild(absentSpan(`this reference row records no ${column.replace(/_/g, ' ')}`));
          } else {
            cell.textContent = String(v);
          }
          tr.appendChild(cell);
        }
        tbody.appendChild(tr);
      }
      el.appendChild(tbody);
      table.appendChild(el);
      referencesHost.appendChild(table);
    }

    const unpinned = provenance.unpinned_references;
    referencesHost.appendChild(h('h4', null, 'Unpinned references'));
    if (unpinned === null || unpinned === undefined) {
      referencesHost.appendChild(
        absentSpan(
          'the manifest writer produced no `unpinned_references` key, so nothing here ' +
            'says whether any reference was unpinned. That is `not reported`, not an ' +
            'empty list meaning "none was unpinned".'
        )
      );
    } else if (unpinned.length === 0) {
      referencesHost.appendChild(
        absentSpan('the run recorded the key and it is empty: no reference was unpinned', 'none recorded')
      );
    } else {
      const ul = h('ul', { class: 'gap-list' });
      for (const id of unpinned) ul.appendChild(h('li', { class: 'mono' }, id));
      referencesHost.appendChild(ul);
      referencesHost.appendChild(
        message(
          'warn',
          'An unpinned reference is reported as a failure, not a silent pass',
          h('p', null, 'A database or tool version that is not pinned is not reproducible.')
        )
      );
    }
  }

  /* -- logs ----------------------------------------------------------- */

  function renderLogs() {
    clear(logsHost);
    logsHost.appendChild(
      h('div', { class: 'panel-head' }, h('h3', null, 'Tripwire and watcher logs'))
    );
    logsHost.appendChild(
      note(
        'No log path for either of these is contracted anywhere in contracts.py or the ' +
          'Snakefile (assumption A7), so an absent log is `not produced` with every probed ' +
          'path named. An empty panel would read as "nothing was reported", which is a ' +
          'different claim from "there is no log".'
      )
    );
    if (!logs) {
      logsHost.appendChild(absentSpan('/api/provenance/logs has not answered', ABSENT.NOT_REPORTED));
      return;
    }
    const items = Array.isArray(logs.items) ? logs.items : [];
    for (const item of items) {
      const box = h('div');
      box.appendChild(h('h4', null, item.name));
      if (!item.present) {
        box.appendChild(absentSpan(item.reason || 'no log was found at any probed path', ABSENT.NOT_PRODUCED));
      } else {
        box.appendChild(
          kv([
            ['path', h('code', null, item.path || 'the server reported no path')],
            ['lines tailed', h('b', null, String(item.n_lines))],
          ])
        );
        const pre = h('div', { class: 'log-tail' });
        for (const line of item.lines || []) pre.appendChild(h('span', { class: 'log-line' }, line));
        box.appendChild(pre);
      }
      const probes = Array.isArray(item.probes) ? item.probes : [];
      if (probes.length > 0) {
        box.appendChild(h('h5', null, 'Paths probed, in order'));
        box.appendChild(probeList(probes));
      }
      logsHost.appendChild(box);
    }
    if (items.length === 0) {
      logsHost.appendChild(absentSpan('the response carries no log entries', ABSENT.NOT_PRODUCED));
    }
  }

  /* -- digests -------------------------------------------------------- */

  function renderDigests() {
    clear(digestsHost);
    digestsHost.appendChild(
      h(
        'div',
        { class: 'panel-head' },
        h('h3', null, 'Artefact digests'),
        h('span', { class: 'faint' }, 'SHA-256, taken at read time')
      )
    );
    digestsHost.appendChild(
      note(
        'A digest describes the file as it is when it is read, not as it was when a run ' +
          'wrote it, so the size and mtime travel with it.'
      )
    );
    const items = digests && Array.isArray(digests.items) ? digests.items : null;
    if (items === null) {
      digestsHost.appendChild(absentSpan('/api/provenance/digests has not answered', ABSENT.NOT_REPORTED));
      return;
    }
    if (items.length === 0) {
      digestsHost.appendChild(
        absentSpan(
          'no stage artefact is readable, so there is nothing to hash. That is `not produced`, ' +
            'not "the run produced files with no content".',
          ABSENT.NOT_PRODUCED
        )
      );
      return;
    }
    const el = h('table', { class: 'digest-table' });
    el.appendChild(headRow(['key', 'sha256', 'size', 'path']));
    const tbody = h('tbody');
    for (const row of items) {
      const tr = h('tr');
      tr.appendChild(h('td', null, row.key));
      const hash = h('td', { class: 'digest-hash' });
      if (row.sha256) hash.textContent = row.sha256;
      else hash.appendChild(absentSpan('the file could not be hashed'));
      tr.appendChild(hash);
      tr.appendChild(h('td', null, typeof row.size === 'number' ? String(row.size) : ''));
      const path = h('td', { class: 'digest-hash' });
      if (row.path) path.textContent = row.path;
      else path.appendChild(absentSpan('no path was reported'));
      tr.appendChild(path);
      tbody.appendChild(tr);
    }
    el.appendChild(tbody);
    digestsHost.appendChild(h('div', { class: 'table-wrap' }, el));
  }

  /* -- git and config ------------------------------------------------- */

  function renderGit() {
    clear(gitHost);
    gitHost.appendChild(h('h3', null, 'Git state of the opened root'));
    if (!provenance || !provenance.git) {
      gitHost.appendChild(absentSpan('the response carries no git block', ABSENT.NOT_REPORTED));
      return;
    }
    const git = provenance.git;
    gitHost.appendChild(
      kv([
        ['HEAD', git.head ? h('code', null, git.head) : absentSpan(git.reason || 'no HEAD was recorded')],
        ['branch', git.branch ? h('code', null, git.branch) : absentSpan('no branch was recorded')],
        [
          'dirty',
          typeof git.dirty === 'boolean'
            ? h('strong', null, git.dirty ? 'yes — the working tree has uncommitted changes' : 'no')
            : absentSpan('git status could not be read'),
        ],
        [
          'git status --porcelain',
          git.status_lines && git.status_lines.length
            ? h('pre', null, git.status_lines.join('\n'))
            : absentSpan(
                typeof git.dirty === 'boolean' && !git.dirty
                  ? 'the working tree is clean, so status printed nothing'
                  : 'no status output was recorded'
              ),
        ],
      ])
    );
  }

  function renderConfig() {
    clear(configHost);
    configHost.appendChild(
      h(
        'div',
        { class: 'panel-head' },
        h('h3', null, 'Configuration snapshot'),
        h('span', { class: 'faint' }, 'verbatim file text')
      )
    );
    configHost.appendChild(
      note(
        'Verbatim text rather than a parsed re-serialisation: a snapshot whose ' +
          'formatting differs from the file it claims to snapshot is not a snapshot.'
      )
    );
    const snapshot = provenance && provenance.config_snapshot ? provenance.config_snapshot : null;
    if (snapshot === null) {
      configHost.appendChild(absentSpan('the response carries no configuration snapshot', ABSENT.NOT_REPORTED));
      return;
    }
    const names = Object.keys(snapshot);
    if (names.length === 0) {
      configHost.appendChild(
        absentSpan(
          'no configuration file was readable at the opened root, so the snapshot is empty. ' +
            'An empty snapshot is not a default configuration.'
        )
      );
      return;
    }
    for (const name of names) {
      configHost.appendChild(h('h4', null, name));
      configHost.appendChild(h('pre', { class: 'snapshot-pre' }, snapshot[name]));
    }
  }

  function renderManifestShape() {
    const box = h('div', { class: 'panel' });
    box.appendChild(h('h3', null, 'Manifest shape'));
    if (!provenance) {
      box.appendChild(absentSpan('/api/provenance has not answered', ABSENT.NOT_REPORTED));
      return box;
    }
    const writer = provenance.manifest_writer || 'none';
    box.appendChild(
      kv([
        ['writer', h('code', null, writer)],
        [
          'what that means',
          writer === 'run'
            ? h('span', null, 'written after the run by papipeline/run.py; it carries `stages`, `stages_skipped`, `outputs` and `tools_detected`')
            : writer === 'provenance'
              ? h('span', null, 'written BEFORE the run by scripts/common/write_provenance.py; it records no stage outcomes, so every stage reads `not run` with the pre-run reason — which is a different sentence from "the run did not record this stage"')
              : absentSpan('no manifest is readable, so there is no writer to name'),
        ],
      ])
    );
    return box;
  }

  function renderAll() {
    clear(errorHost);
    if (!provenance) {
      errorHost.appendChild(
        message(
          'error',
          'Provenance could not be read',
          note('Nothing on this page is shown from a cache or a default; each panel states its own absence.')
        )
      );
    }
    renderProof();
    renderTools();
    renderReferences();
    renderLogs();
    renderDigests();
    renderGit();
    renderConfig();
    const existing = container.querySelector('.manifest-shape-slot');
    if (existing) existing.remove();
    const slot = h('div', { class: 'manifest-shape-slot' });
    slot.appendChild(renderManifestShape());
    container.insertBefore(slot, toolsHost);
  }

  (async () => {
    const results = await Promise.allSettled([
      ctx.api.get('/provenance', null, abort.signal),
      ctx.api.get('/provenance/bakta', null, abort.signal),
      ctx.api.get('/provenance/digests', null, abort.signal),
      ctx.api.get('/provenance/logs', { name: 'any', limit: 50, tail_lines: 200 }, abort.signal),
    ]);
    if (disposed) return;
    const failures = [];
    results.forEach((result, index) => {
      const names = ['/api/provenance', '/api/provenance/bakta', '/api/provenance/digests', '/api/provenance/logs'];
      if (result.status === 'fulfilled') {
        if (index === 0) provenance = result.value;
        if (index === 1) bakta = result.value;
        if (index === 2) digests = result.value;
        if (index === 3) logs = result.value;
      } else if (result.reason && result.reason.name === 'AbortError') {
        return;
      } else {
        failures.push(`${names[index]}: ${result.reason && result.reason.message ? result.reason.message : result.reason}`);
      }
    });
    if (failures.length > 0) {
      clear(errorHost);
      errorHost.appendChild(
        message(
          'warn',
          'Some provenance endpoints did not answer',
          h('ul', { class: 'gap-list' }, failures.map((f) => h('li', { class: 'mono' }, f)))
        )
      );
    }
    renderAll();
  })();

  return () => {
    disposed = true;
    abort.abort();
  };
}

function headRow(columns) {
  const thead = h('thead');
  const tr = h('tr');
  for (const column of columns) tr.appendChild(h('th', null, column));
  thead.appendChild(tr);
  return thead;
}

function countOrUnknown(value) {
  return typeof value === 'number' ? String(value) : 'an unreported number';
}

function rowsOrUnknown(rows) {
  if (Array.isArray(rows)) return String(rows.length);
  return 'an unreported number of';
}
