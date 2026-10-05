# The Computational Observatory

A read-only interface onto PAPipeline's execution state. It shows what the
pipeline is doing: which stage, for which genome, in which state, having
run which command, with which result — and it shows nothing it cannot
substantiate.

It is not a simulator, not a dashboard with plausible numbers, and not a
neural animation. Every value on screen is read from the execution store or
from the host, at the moment it is displayed.

```
┌───────────────────────────────────────────────────────────────────────┐
│ PAPipeline · Computational Observatory      LIVE   2026-09-27 14:32:18│
├──────────────────────┬────────────────────────────────────────────────┤
│ CURRENT TASK         │  Computational Task Network                    │
│  genome, stage,      │                                                │
│  state, attempt,     │   stage nodes coloured by real state,          │
│  command, inputs,    │   connected by the pipeline's own declared      │
│  outputs, validation │   dependencies, with a mark travelling an       │
│  log, retry history  │   edge when a task on it really completes       │
│                      │                                                │
│                      ├────────────────────────────────────────────────┤
│                      │ cohort progress · counters · real host metrics  │
└──────────────────────┴────────────────────────────────────────────────┘
```

## Running it

```bash
# Watch a real cohort execute, in a temporary store
python -m papipeline.observatory --demo

# Read a real store from a real run
PAPIPELINE_OBSERVATORY_DB=/path/to/execution.db \
PAPIPELINE_OBSERVATORY_RUN=<run-key> \
python -m papipeline.observatory
```

Then open <http://127.0.0.1:8765>.

`--demo` is not a mock. It runs the real `run_task` against real
subprocesses, writing to a real SQLite store, and the interface reads that
store. To see the failure states, watch a subject that is retried, one
whose outputs fail their contract, and one killed part-way.

Requires `fastapi`, `uvicorn` and `psutil` (added to
`environment/environment.yml` as optional). Without `psutil` the interface
still runs and reports CPU, memory and disk I/O as unavailable.

## What is real

| Shown | Source |
| --- | --- |
| Task state, attempt, timings, command, versions | `execution` table, written by `run_task` |
| Retry history | `attempt_log`, one row per attempt |
| Input / output identifiers | `input_ids`, `output_paths` columns |
| Validation detail | `validation_json`, the real check results |
| Log text | the file the runner actually wrote |
| Node colour and counts | aggregated from the store, per stage |
| Stage success/failure | the contract in `papipeline/execution/contracts.py` |
| Edges | `papipeline.run.PREREQUISITES` and `EXECUTION_ORDER` |
| CPU, RAM, disk I/O, load, uptime | `psutil` / `getloadavg`, sampled live |
| Live events | transitions `run_task` emitted through `event_sink` |

## What is deliberately absent

These are not oversights. The engine does not produce them, so the
interface shows them as unavailable rather than inventing them.

- **Per-task progress percentage.** There is no progress reporting in
  `run_task` and no column for it. The cohort progress bar is real — it is
  counted from finished tasks — but there is no per-genome "72%".
- **Per-task CPU and memory.** The store keeps a static host capacity
  snapshot from write time, and the runner does not surface a child process
  id. So the bottom bar shows *host* metrics, and the per-task tiles read
  "not reported".
- **Worker identity.** There is no worker pool in this build. Labelling a
  task "W03" would be inventing a scheduling fact.
- **ETA below three completed tasks.** A rate from one or two samples is
  noise. `throughput.py` returns `None` and states why.

Each of these is reported in `/api/snapshot` under `capabilities`, so a
client can tell "unavailable" from "not yet wired".

## Watching a real pipeline run

`run.py` now runs every stage inside `run_task`, so a normal run populates
this interface. See `docs/execution_wiring.md` for the detail.

```bash
# terminal 1 — the run, observed
python3 scripts/common/run_pipeline.py --mode TEST --observatory \
    --observatory-db /tmp/observed.db --observatory-run demo

# terminal 2 — this interface, reading the same store
python3 -m papipeline.observatory --db /tmp/observed.db --run demo
```

The pipeline runs with the observatory off unless asked, and needs neither
the store nor the UI. `--observatory` is opt-in.

## Design decisions

**Zero-build frontend.** The repository had no JavaScript, no `package.json`
and no bundler. Adding React Three Fiber and Vite would have introduced a
second toolchain and a build step for a scene of 16 nodes and 32 edges. The
network is Canvas 2D with additive compositing; the layout is a
deterministic barycentric sweep with no random seeding, so the same
structure always draws the same picture — which is both calmer to look at
and testable.

**Two edge classes, never merged.** `dependency` edges come from
`PREREQUISITES` and mean the target reads the source's output.
`sequence` edges come from `EXECUTION_ORDER` and mean only "scheduled
first". A pair can be both — `phenotype → gwas` is adjacent *and* declared —
so both are drawn, styled differently, and the legend says which is which.
Collapsing them would overstate the pipeline's dependencies.

**Worst member wins.** A stage with 64 successes and 1 failure renders as
`FAILED`, not `SUCCEEDED`. The tooltip breaks the count down by state.

**Idleness is honest.** The network is still when nothing is running. An
ambient starfield is generated from a fixed seed and does not move, because
motion that is not caused by work is a lie about causation. Travelling
marks appear only on a real `TASK_COMPLETED`.

**The store is authoritative.** Events are notifications, not state. Every
SSE heartbeat re-sends a full snapshot, so a dropped or coalesced event
costs a frame of freshness and nothing else. A client that reconnects
re-reads rather than resuming on a stale picture.

## Layout

| Width | Behaviour |
| --- | --- |
| ≥ 1180px | left 33% detail, network and metrics stacked right |
| 860–1180px | left 40%, counters and sparklines reflow to two columns |
| < 860px | stacked, network above details |

Desktop-first, as a scientific instrument on a desk should be.

## Tests

`tests/unit/test_observatory.py`, 56 tests:

- structure comes only from `PREREQUISITES` / `EXECUTION_ORDER`; nothing invented
- state aggregation, including that a stage is as bad as its worst member
- unmeasurable metrics read `None`, never `0`
- throughput refuses a rate below three samples
- events are the runner's, and a broken subscriber cannot break a run
- every API value matches the store row it came from
- a real cohort's real states appear through the real API
- SSE frames delivered over real HTTP to a real uvicorn server
- the shipped assets contain none of the reference dashboard's numbers
