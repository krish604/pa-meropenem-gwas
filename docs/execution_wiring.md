# Production execution wiring

How `run.py` feeds the Computational Observatory, and how to run a real
observed pipeline.

## What changed

Before this milestone the execution substrate was real but orphaned:
`papipeline.execution` was imported by its own tests and nothing else, so a
production run wrote no execution state and the observatory had nothing to
show. `run.py` executed its sixteen stages through an `if`/`elif` chain and
returned.

Now:

```
run.py
  └─ _run_one_stage
       └─ run_task(...)          ← execution state, validation, retry, events
            └─ stage.run(...)    ← the pipeline's own stage body, unchanged
       └─ validate the stage's declared outputs
       └─ persist state, publish the transition
            ↓
       ExecutionStore (SQLite) → EventBus → SSE → Observatory
```

The stage body is not reimplemented, reordered or forked. It was moved into
a closure and is invoked **once**, through `run_task`. There is one task
engine.

### One engine, two invokers

`run_task` originally only knew how to spawn a subprocess, which suits Bakta
and the other external tools but not the fifteen stages that are ordinary
Python functions. Rather than add a second execution path, `run_task` now
accepts either:

- an argument vector, invoked as a subprocess; or
- a zero-argument callable, invoked in-process.

Both travel the identical path: same `TaskContext`, same `OutputSpec`, same
validation, same retry policy, same store rows, same events. A callable is
recorded as `<in-process> <qualname>` rather than being stringified into
something that looks like a shell command.

A caveat, stated rather than hidden: a callable cannot be killed the way a
process can, so `ctx.timeout` is recorded but not enforced for in-process
stages. That was already true before observability existed.

## Enabling it

Off by default, in three ways. The pipeline never requires the observatory,
and a normal run needs neither a store nor a UI.

```bash
# 1. no flag, no config: observability off, behaviour identical
python3 scripts/common/run_pipeline.py --mode TEST

# 2. config
#    config/config.yaml
#    runtime:
#      observatory:
#        enabled: true
#        db: ~/.local/share/papipeline/observatory.db
#        run_key: pipeline

# 3. command line
python3 scripts/common/run_pipeline.py --mode TEST --observatory
python3 scripts/common/run_pipeline.py --mode TEST \
    --observatory-db /tmp/observed.db --observatory-run my-run
```

An explicit `--observatory-db` or `--observatory-run` implies the flag, so
you do not need both.

## Watching a real run

```bash
# terminal 1 — the run, observed
python3 scripts/common/run_pipeline.py --mode TEST --observatory \
    --observatory-db /tmp/observed.db --observatory-run demo

# terminal 2 — the interface, reading that same store
python3 -m papipeline.observatory \
    --db /tmp/observed.db --run demo
```

The store is the only thing they share. The pipeline writes it; the
interface reads it and subscribes to the bus. Neither needs the other
running, which is why an interrupted run can be inspected afterwards.

## How events reach SSE

`run_task` announces each state transition it has already made and recorded.
There is no separate notifier and no second state machine.

```
run_task → TaskEvent → observatory.events.bus_sink → EventBus
        → ExecutionEvent → SSE /api/stream → browser
```

`SNAPSHOT` remains authoritative. An event says something changed; the
snapshot says what is now true. The stream re-sends a full snapshot on every
heartbeat tick, so a dropped or coalesced event costs a frame of freshness
and nothing else. A client that reconnects re-reads rather than resuming on
a stale picture.

Event names in use: `TASK_CREATED`, `TASK_STARTED`, `TASK_COMPLETED`,
`TASK_VALIDATED`, `TASK_FAILED`, `TASK_INVALID`, `TASK_INCOMPLETE`,
`TASK_RETRY`, `TASK_RESUMED`, plus `STAGE_UNBLOCKED` and `STAGE_FAILED` for
the pipeline's own transitions.

A bus belongs to one observed run. A client is only ever shown its own run's
events, so two observed pipelines can coexist without contaminating each
other.

## Task identity

`(run_key, stage, subject)` — the same triple the store is keyed on.

- `subject=""` for a cohort-level stage.
- `subject=<sample_id>` for annotation, because annotation is per genome:
  Bakta runs once per isolate, so that is the granularity at which a
  failure is attributable and the granularity at which it is observed.

No random component is introduced, deliberately: a random id would make
resume and provenance impossible. Attempt numbers are allocated by the store
and are monotonic for the life of the task, so a re-run continues the
history (1, 2, 3) rather than erasing it.

## What each stage is held to

`papipeline/execution/contracts.py` declares, per stage, the exact filename
and header of the table the pipeline writes. `run.py` writes through
`table_path(stage_dir, stage)` and the contract validates the same path, so
the file a stage writes and the file a contract checks cannot drift apart.
`tests/integration/test_run_observability.py::TestContractsMatchThePipeline`
fails if they ever do.

The default checks are: exists, non-empty, parses as TSV, carries the
declared header, holds at least one record.

Cohort coverage is **opt-in** and off by default. Whether every isolate has
a determinant call is a statement about the data — absence is a finding, not
a missing record — so asserting it as a default execution contract would
fail a correct run whose correct answer is "this isolate has no AMR gene".

## Bakta

`papipeline/adapters/bakta.py` builds the real invocation:

```
bakta --db <pinned db> --in <assembly.fna> --out <dir> --threads <config>
```

`--threads` is passed through from configuration unchanged. This milestone
does not manage Bakta's concurrency.

Bakta's output is converted into the pipeline's own schema by the same
`parse_bakta_tsv` / `standardise` functions TEST mode uses, and written where
`load_from_intermediate` reads it, so there is one downstream code path and
no way for downstream stages to tell how the records were obtained.

REAL mode is still gated by `runtime.allow_real_mode: false` in
`config.yaml`. Implementing the adapter did not unlock it.

## Execution semantics: NEW RUN, RE-RUN, RESUME, RETRY

These four were previously conflated. The discriminator is the **run key**.

| Situation | Definition | What happens | Run event |
| --- | --- | --- | --- |
| **NEW RUN** | the run key has no recorded state | everything executes | `RUN_STARTED` |
| **RE-RUN** | the run key has state, caller did *not* ask to resume | everything executes again | `RUN_STARTED` |
| **RESUME** | the run key has state *and* `--resume` was passed | validated stages are skipped, the rest execute | `RUN_RESUMED` |
| **RETRY** | one task, second attempt, inside one invocation | the task re-executes; history grows | `TASK_RETRY` |

The rule that was previously wrong, and is now enforced by test:

> **`TASK_RESUMED` is emitted only by a genuine resume.**

`announce_resume()` used to fire whenever an invocation found a prior
`SUCCEEDED` row for the same key — so an ordinary second run announced
`RUN_RESUMED`-equivalent events for work it had in fact redone, while its
own detail text said "re-executing". It now emits only from
`announce_skipped()`, which fires once per task a resume genuinely did not
re-execute. The count of `TASK_RESUMED` events is therefore the count of
work a resume actually saved.

```bash
# a re-run: executes everything, no TASK_RESUMED
python3 scripts/common/run_pipeline.py --mode TEST --observatory

# a resume: skips what still validates
python3 scripts/common/run_pipeline.py --mode TEST --observatory --resume
```

`papipeline/execution/semantics.py` owns the classification and the
skip decision, including the rule that a recorded `SUCCEEDED` whose files
were deleted is re-executed rather than skipped.

## Failure behaviour

| What happened | State | Retried? |
| --- | --- | --- |
| Stage completed, outputs satisfy the contract | `SUCCEEDED` | — |
| Stage completed, outputs absent or empty | `INCOMPLETE` | no¹ |
| Stage completed, outputs present but wrong | `INVALID` | no |
| Stage raised / tool exited non-zero, wrote nothing | `FAILED` | no¹ |
| Stage raised after writing part of its output | `INCOMPLETE` | no¹ |
| Validated success, re-run | `TASK_RESUMED` | — |

¹ The pipeline has always stopped at the first stage failure, so
`RetryPolicy(max_attempts=1)` is deliberate. Retrying a scientific stage
would be a behaviour change dressed up as resilience. Retry policy is
exercised on the external-tool path, where a transient tool failure is a
real possibility.

One correctness detail worth recording. A failed stage whose output table
already exists **from a previous run** used to be reported `INCOMPLETE`,
because the stale file was indistinguishable from one this attempt had
half-written. Outputs are now only counted as partial if they were modified
during the attempt, so an immediate failure is `FAILED` and a genuine
partial write is still `INCOMPLETE`.

## Overhead

Measured on the TEST-mode pipeline, 20 samples, all 16 stages, five repeats:

| | wall (median) | CPU (median) | peak RSS |
| --- | --- | --- | --- |
| observatory off | 5.380 s | 0.106 s | 81.4 MB |
| observatory on | 5.274 s | 0.124 s | 81.4 MB |

Wall-time difference (−106 ms) is inside the run-to-run standard deviation
(0.31–0.83 s), i.e. not distinguishable from noise. The measurable cost is
about 19 ms of CPU across 16 stages, roughly 1.2 ms per stage, and no
memory growth. Not yet optimised, and not yet worth optimising.

## Limitations

- **No per-task progress.** The engine reports none, so none is shown. The
  cohort completion bar is counted from real finished tasks.
- **No per-task CPU or memory.** The runner does not surface a child pid, so
  the bottom bar reports *host* metrics and per-task tiles read
  "not reported". There is no worker pool, so no worker identity is shown
  either.
- **Cohort stages are one task, not one per genome.** Stages 3–16 each write
  a single cohort-level table, so that is the granularity at which they are
  observed. Only annotation is per genome, because only annotation runs per
  genome. Finer granularity would mean changing what those stages write.
- **Coverage is not asserted.** See above.
- **A second run is a re-run, not a resume.** `TASK_RESUMED` is emitted when
  a prior `SUCCEEDED` row exists, and the attempt history continues. Stages
  are not skipped, because skipping them would change what the pipeline
  computes. `papipeline.execution.resume.decide` already implements
  validation-based skipping; wiring it to skip work is a separate decision
  and was not in scope here.
- **REAL mode still cannot run end to end** without the pinned toolchain
  (Bakta, MLST, AMRFinderPlus, minimap2, …). The Bakta path is implemented
  and tested against a real executable; the others are observed as in-process
  stages, not as external tool invocations.
