# The results dashboard

A read-only browser over a run's artefacts: the 16 pipeline stages, the joined
per-isolate table, the trees, the similarity matrix, GWAS / convergence /
co-occurrence, and the provenance that makes all of it checkable.

It **never** writes to the results root. It **never** starts a run by default.
It **never** fetches anything from another origin. The contract it implements is
`dashboard/DESIGN.md`; the HTTP surface is `dashboard/openapi.yaml`.

---

## Start it

`python` is not on PATH in this project; the environment is a micromamba one.
Use the wrapper, which delegates every decision to `python -m dashboard`:

```sh
scripts/run_dashboard.sh --results <path-to-results-root>
```

or the module directly from an activated environment:

```sh
eval "$(micromamba shell hook -s bash)" && micromamba activate pa-amr
python -m dashboard --results <path> [--port 8765] [--allow-launch]
```

Then open <http://127.0.0.1:8765/>.

### Flags

| flag | default | what it does |
|---|---|---|
| `--results PATH` / `--results-root PATH` | auto-detect | the results root to read |
| `--bundle PATH` | auto-detect | a delivery-bundle root |
| `--state-dir PATH` | `~/.pa_dashboard` | index + log location; must be **outside** the results root |
| `--host HOST` | `127.0.0.1` | bind address; a non-loopback host requires `PA_DASH_TOKEN`, and `0.0.0.0` is always refused |
| `--port PORT` | `8765` | bind port |
| `--machine NAME` | config default | machine overlay for the config (`laptop`, `bigmachine`) |
| `--allow-launch` | off | enable the launcher panel (a dry run only; it starts nothing) |
| `--no-serve` | off | build the app, report what opened, exit |

### Environment variables

| variable | meaning |
|---|---|
| `PA_DASH_RESULTS_ROOT` | probe 2 of the auto-detect order |
| `PA_DASH_BUNDLE` | probe 3, a delivery bundle |
| `PA_DASH_STATE_DIR` | overrides `~/.pa_dashboard` |
| `PA_DASH_TOKEN` | required for a non-loopback bind; sent as the `X-PA-Dash-Token` header |
| `PIPELINE_RESULTS_ROOT` | honoured by the pipeline's own config loader |

### Auto-detect

With no `--results`, six probes are tried in order and **all** are recorded, so a
failure names every place that was looked at: `--results` → `PA_DASH_RESULTS_ROOT`
→ a bundle (`PA_DASH_BUNDLE` / `--bundle` / the flag's sibling) → the live REAL,
TEST and STUB roots from the machine overlay's `paths.results_root`. The failure
is a `503` whose message lists each probe and why it did not match; it keeps
*not configured* apart from *configured and holds nothing*.

---

## The two layouts

**Live results root** — what a run writes:

```
<root>/
  run_manifest.json                 # written by run.py:_write_run_manifest
  intermediate/stages/*.tsv         # the 16 stage tables + folded steps
  intermediate/phylogeny/*.nwk      # the stage-9 tree
  reports/                          # the run's own .md/.html/figures
  status/events.jsonl               # the event log the monitor tails
```

**Delivery bundle** — the A1–A5 hand-over layout:

```
<root>/
  01_bakta_input/
  02_stage_outputs/                 # the bundle's copy of intermediate/stages/
  03_report/                        # the .md/.html/figures
  04_run_info/                      # run_manifest.json, status/events.jsonl
  05_validation/
  06_for_900_isolates/RUNBOOK_900.md
```

The `kind` (`live` or `bundle`) is decided by **where `run_manifest.json`
landed**, not by which probe matched, so a bundle missing `02_stage_outputs/` but
carrying the manifest at its root is still served correctly. Every bundle
directory is optional and probed.

A committed 10-isolate pair of both layouts lives under
`dashboard/fixtures/small/` for development and the test suite.

---

## Viewing it on the big machine

The dashboard binds `127.0.0.1` by default and is meant to be reached over an
SSH tunnel. On the analysis machine:

```sh
scripts/run_dashboard.sh --results /path/to/results
```

From your workstation:

```sh
ssh -L 8765:127.0.0.1:8765 you@bigmachine
```

then open <http://127.0.0.1:8765/>. The token is **not** needed for a loopback
bind.

If you must bind a non-loopback address, set a token first — the server refuses
to start without one, and refuses `0.0.0.0` even with one:

```sh
PA_DASH_TOKEN=$(openssl rand -hex 24) \
  python -m dashboard --results /path --host 10.0.0.5
```

Every `/api/**` request then needs `X-PA-Dash-Token: <token>`. The token is
**never** accepted as a query parameter: a query string lands in shell history,
in the address bar and in every proxy log.

---

## Security model

- **UI-D1 — read-only.** The results root is opened once, `realpath()`d, and
  every file under it is opened `'r'`. `ResultsSource` exposes no write method,
  and a middleware refuses every write verb except `POST /api/launcher/preflight`
  (which starts nothing) with a `405`. The **only** write location is the state
  directory (`~/.pa_dashboard/`), which is outside the results root by
  construction; `python -m dashboard` refuses a `--state-dir` inside the root at
  startup.
- **UI-D4 — containment, allowlist, bind.** A served path is resolved (following
  symlinks) and required to be under the opened root; an absolute path or a `..`
  component is refused rather than normalised; only a fixed suffix allowlist is
  served; `text/html` is served as an attachment with `X-Content-Type-Options:
  nosniff`. A non-loopback bind requires `PA_DASH_TOKEN`, compared with
  `secrets.compare_digest`, header only.
- **UI-D5 — huge tables are never fully loaded.** A page is `seek` +
  `readline` over a byte-offset index cached in the state directory, keyed on
  `sha256(realpath + size + mtime_ns)`. `meta.read` on every paged response
  reports `bytes_read`, `file_size`, `fully_loaded` and the mechanism, so the
  claim is checkable per request.
- **UI-D6 — tree identity is a postorder index**, never the Newick label. The
  fixtures repeat `100/100` and omit labels precisely so a label-keyed client
  breaks loudly in the test suite.

---

## Data privacy (D7)

**Real accessions belong in the browser, not in this repository.** The committed
fixtures, the tests, the docs and every `pa-artifacts/ui/` report use synthetic
`TEST_PA_###` ids only. `dashboard/tests/test_q1_fixtures.py` asserts that no
`PDT…` / `GCF_…` / `SAMN…` / `SAMEA…` accession appears under
`dashboard/fixtures/`.

When you point the dashboard at a real results root, the accessions are served
from that root at request time and are never copied into the repository, the
state directory's index (which stores byte offsets and the key column's values
only), or any log.

---

## What is not built (honestly)

- **The launcher cannot start a run.** It is disabled unless `--allow-launch`,
  and even then the only action is a dry run that is *displayed*, never
  executed. A REAL run would need the phrase `run real samples` typed exactly,
  the `bigmachine` overlay, and `PIPELINE_ALLOW_REAL_MODE` on the child process
  — and no endpoint executes one. This is DESIGN §13's position.
- **No recents persistence.** The source picker renders the honest absence: the
  browser writes nothing locally, and there is no endpoint that re-points the
  running server at a different root.
- **No `prerequisites` endpoint.** The pipeline DAG's edges are derived from the
  event log's own timestamps and labelled "observed"; where the log carries no
  transition evidence, the server's declared stage order is drawn dashed and
  labelled "declared order only".
- **SSE is hand-rolled** on Starlette's `StreamingResponse` (no `sse_starlette`
  in this environment). The polling fallback (`/api/events/state`) is the
  client's fallback when the stream does not open.
- **No origin/CSRF check.** `POST /api/launcher/preflight` starts nothing and
  takes no state-changing action; a future launcher POST must add one before it
  may.
- **d3 zoom and colour interpolation are hand-written.** `d3-zoom` cannot
  construct without `d3-dispatch`/`d3-drag`/`d3-transition`, and
  `d3.interpolateRgb` needs `d3-color`; neither is vendored, so pan/zoom and the
  heatmap ramp are written directly.
- **A 4,001-column projection is a `400`** (the query string is ~91 KB before a
  handler runs). The tables page pages the projection 40 columns at a time.
- **The observatory is not read.** The dashboard reads `status/events.jsonl`;
  `papipeline/observatory/`'s in-process `EventBus` is a second mechanism that
  this build does not depend on.

---

## Tests

```sh
eval "$(micromamba shell hook -s bash)" && micromamba activate pa-amr
PA_FIXTURES_LARGE=1 python -m pytest dashboard/tests -q
node --test dashboard/tests/js/
```

`PA_FIXTURES_LARGE=1` generates the 900-isolate root in a temp dir (never
committed) and runs the three performance gates. The JS suite needs only
`node`, no npm install.
