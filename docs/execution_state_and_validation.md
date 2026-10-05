# Execution state and output validation

This document describes the execution substrate added to PAPipeline: how a
stage's state is decided, how its outputs are declared and checked, how
failures are retried, and how a resumed run decides what to skip.

It is infrastructure only. No scientific algorithm, threshold, filter or
interpretation was changed to accommodate it.

## The problem this solves

The pipeline shipped a set of failures that all looked like successes:

| Symptom | Cause | What a file-existence check would say |
| --- | --- | --- |
| Every stage-6 locus scored 0/65 | `bakta_outputs` returned the *inference* table, so gene names were empty | "annotation present" |
| Per-sample VCFs byte-identical, one distinct sequence in the alignment | A leaked loop variable | "VCF present, 65 files" |
| Every kinship entry exactly 1.0 | Matrix built from a called/not-called mask | "kinship present" |
| Runs died at 90% and resumed as complete | Partial output left on disk | "files present" |

In each case the artefacts existed, were non-empty, and parsed. Existence
is not evidence of completion.

## The invariant

> `SUCCEEDED` means the process exited cleanly **and** every declared output
> check passed. No other combination produces `SUCCEEDED`.

The invariant is enforced structurally, not by convention:

- `TaskContext` requires an `OutputSpec` and rejects a spec with no checks.
  An empty contract would pass unconditionally, which is the same as having
  no contract.
- `decide()` requires a spec. Skipping a recorded `SUCCEEDED` without
  re-validating the files is exactly the failure this module exists to
  prevent, so it raises rather than skipping.
- `build_plan()` treats a task it has no contract for as `RUN`, not `SKIP`.

A stage therefore cannot report success it has not earned.

## States

| State | Meaning | Retryable |
| --- | --- | --- |
| `PENDING` | Never attempted in this run | yes |
| `RUNNING` | A process is executing | — |
| `SUCCEEDED` | Clean exit **and** contract passed | — |
| `INCOMPLETE` | Missing, empty or truncated output | yes |
| `INVALID` | Output is present but wrong | **no** |
| `FAILED` | The process itself failed | depends on cause |
| `RETRYING` | A retryable failure; another attempt is scheduled | — |

`INVALID` is never retried. Re-running a deterministic producer reproduces
the same wrong output; retrying only delays the report and hides the cause.

`INCOMPLETE` is retryable because a truncated file is usually the residue
of an interrupted run rather than a wrong answer.

### Two distinctions the return code cannot make

**Died part-way vs. never started.** A process that exits non-zero *and left
output on disk* has produced bytes that downstream stages would otherwise
read as authoritative. That is `INCOMPLETE` and it is retryable. A process
that failed having written nothing is `FAILED` — it did no work, and the
cause is usually deterministic. `_partial_output()` makes this explicit
rather than inferring it from the exit status.

**A clean exit is not success.** Exit 0 with missing output is
`INCOMPLETE`; exit 0 with malformed output is `INVALID`. Both are decided
after the process ends, never from its exit code.

## Declaring a contract

A contract is a set of checks over a stage's outputs:

```python
from papipeline.execution import Check, CheckKind, OutputSpec

spec = OutputSpec(
    stage="standardised_annotation",
    checks=(
        Check(CheckKind.NON_EMPTY, path=path),
        Check(CheckKind.PARSEABLE, path=path, parser="tsv"),
        Check(CheckKind.COLUMNS, path=path,
              columns=ANNOTATION_REQUIRED_COLUMNS, dialect="pipeline"),
        Check(CheckKind.MIN_ROWS, path=path, minimum=1, dialect="pipeline"),
        Check(CheckKind.NO_SENTINEL_COLUMN, path=path, key="gene_name",
              dialect="pipeline"),
    ),
    expectations={"sample_id": "S1"},
)
```

Checks are declarative data, so a contract can be stored, diffed and
compared across runs.

### Check kinds

| Kind | Fails when | Failure state |
| --- | --- | --- |
| `EXISTS` | the path is absent | `INCOMPLETE` |
| `NON_EMPTY` | the path is absent or zero bytes | `INCOMPLETE` |
| `MIN_SIZE` | the file is smaller than `minimum` | `INCOMPLETE` |
| `MIN_ROWS` | fewer than `minimum` records are readable | `INCOMPLETE` |
| `PARSEABLE` | the file cannot be read by the named `parser` | `INVALID` |
| `COLUMNS` | required `columns` are absent from the header | `INVALID` |
| `NO_SENTINEL_COLUMN` | every value in `key` is empty or a sentinel | `INVALID` |
| `SAMPLE_COUNT` | the file accounts for fewer than `expected` samples | `INVALID` |
| `IDENTIFIERS` | an identifier in `key` is absent from the manifest | `INVALID` |
| `MATRIX_NOT_CONSTANT` | every cell of a matrix holds the same value | `INVALID` |
| `IDENTIFIES_FILE` | the file is not named after the submitted genome | `INVALID` |

`MIN_ROWS` is what makes "killed at 90%" detectable: a truncated file has
fewer records than the contract requires.

`MATRIX_NOT_CONSTANT` guards the identity-by-state bug: a kinship matrix
whose every entry is 1.0 is a mask, not a relationship.

### Sibling degeneracy

Per-sample files are compared with each other through `SiblingSpec`
(`root`, `pattern`, `min_distinct`), attached to an `OutputSpec`. Two
isolates' VCFs that are byte-for-byte identical are not two isolates' VCFs,
and this is the direct guard against the leaked-loop-variable bug.

### Dialects

Checks read files through a named dialect rather than assuming TSV, because
the pipeline's formats disagree about columns:

| Dialect | Reads | Keyed by |
| --- | --- | --- |
| `pipeline` | TSV via the pipeline reader | pipeline column names |
| `bakta` | Bakta TSV, normalised | `seqid`, `gene_id`, `gene_name`, … |
| `bakta_raw` | Bakta TSV, unnormalised | the file's own names (`Gene`, `Locus Tag`) |
| `vcf` | VCF, one dict per record | the file's `#CHROM` column names |

`bakta` and `bakta_raw` are both required. The normalised view is what
stages 6 and 9 consume, but a feature table and an inference table
normalise to the *same* key set, so only the raw view distinguishes them.
`Gene` and `Score` are the discriminator.

Reading a VCF with the `pipeline` dialect promotes a data row to be the
header. The `vcf` dialect exists so that never happens.

## Running a task

```python
result = run_task(ctx, ["bakta", "--in", faa, "--out", out_dir], store=store)
```

`run_task` executes, validates, retries per policy, and records every
attempt. The verdict is reached only after the process has exited.

`env=` is layered over the parent environment rather than replacing it, so
a stage can be given a database directory or a tool path without silently
losing `PATH`.

## Retry

`RetryPolicy` retries `INCOMPLETE` and the retryable `FAILED` causes —
`RESOURCE`, `TIMEOUT` — with exponential backoff capped at `max_delay`.
Deterministic failures (`EXECUTION`, `CONFIG`, `TOOL_MISSING`) and `INVALID`
are not retried.

Every attempt is written to the store, including the failed ones, so a
late success still has complete provenance.

### Attempt numbering

Attempt numbers are allocated by the store and are monotonic for the life
of the `(run, stage, subject)` triple. A task that is interrupted and later
resumed continues its history (1, 2, 3, …) rather than restarting and
erasing the earlier evidence. A late success must be able to show what went
wrong before it.

## Provenance

The store is a SQLite database with WAL enabled and foreign keys enforced.
One row per task holds state, attempt count, command, exit code, tool and
database versions, configuration hash, input identifiers, output paths,
validation verdict and detail, error text, host, CPU count and start/end
times. `attempt_log` holds one row per attempt.

The configuration hash is what makes a cached success honest: if the
configuration changed, the recorded success no longer describes the current
run and the task is re-run.

## Resume

`decide()` returns one of `SKIP`, `RUN`, `RERUN`, `BLOCKED`.

| Recorded state | Action | Why |
| --- | --- | --- |
| nothing recorded | `RUN` | cold start |
| `SUCCEEDED`, outputs re-validate | `SKIP` | the claim holds up |
| `SUCCEEDED`, outputs do **not** re-validate | `RERUN` | deleted, truncated or corrupted |
| `SUCCEEDED`, configuration hash differs | `RERUN` | the result describes a different run |
| `FAILED` / `INCOMPLETE` / `INVALID` | `RERUN` | never completed |
| prerequisite did not succeed | `BLOCKED` | never execute a stage on a failed one |

A recorded `SUCCEEDED` is never believed without looking at the files. This
is the single most important line of behaviour in the module: the store
says what happened, the filesystem says what is true, and only the
filesystem is allowed to end a run.

## Verification

```
781 passed
```

- 641 pre-existing tests, unchanged.
- 140 new tests: state and retry, the validator, and integration tests that
  execute real subprocesses.

The integration tests run **real processes** and inspect the bytes they
leave on disk. The audit found 99 green tests over helper functions
coexisting with an engine at 0% coverage; the new tests close that gap for
the part of the system built so far.

`tests/integration/test_real_stage_execution.py` runs the actual
`papipeline.stages.annotation` code in a subprocess over real Bakta-format
files, and holds the result to its contract — including the killed-at-90%
case, where a plausible partial table is on disk and the process is dead.

## Limits

**This module is now wired into the pipeline.** `run.py` calls `run_task()`
for every stage; see `docs/execution_wiring.md`. Two things remain open,
both honest gaps rather than silent ones:

1. **Stages 3–16 are observed as in-process stages, not tool invocations.**
   `run_task` executes them as callables, so their state, validation and
   events are real, but their external tools are called by each stage's own
   code rather than through `run_task`'s subprocess path. Only annotation
   drives an external tool through `run_task` directly.
2. **No progress, per-task CPU, per-task memory or worker identity.** None
   of these is recorded by the engine, so the interface reports them as
   unavailable instead of estimating them.

The contracts in `specs.py` encode what the pipeline already requires. Where
the current pipeline does not consume an artefact — FAA and FFN are not
read by the standardisation path — the check is optional and defaults off
rather than being asserted.
