# Status

Last real-data run: round 12, `integration-r3` at `3f1c633`.

## Summary

| | |
|---|---|
| Stages implemented | 16 of 16 |
| Stages verified on **real** data | 1, 2, 3, 4, 4a, 5, 6, 6a (8 of 16) |
| Stages never run on real data | 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 |
| Tests | 3502 collected — **3434 passed, 68 skipped, 0 failed** |
| Cohort used | 10 isolates (`local/smoke_isolates.txt`) |
| `bakta` / `panaroo` executions | **0** — Bakta is an input, never re-run |

**This is not a complete 16-stage result and nothing here is a biological
finding.** n = 10 is grossly underpowered.

## Per-stage, on the 10-isolate smoke subset

| # | Stage | Real-data outcome | Evidence |
|---|---|---|---|
| 1 | `validation` | ran | 10/10 assemblies pass |
| 2 | `annotation` | ran | 10/10 `reused=true reason=verified`; no tool executed |
| 3 | `mlst` | ran | 10/10 typed, scheme *paeruginosa* |
| 4 | `amr` | ran | 199 rows; AMRFinderPlus 4.2.7, DB 2026-08-07.1 |
| 4a | structural variants | ran | 409 candidate calls, 132 unassessable regions |
| 5 | `virulence` | ran | VFDB 4573 factors |
| 6 | `variants` | ran | 509,612 alleles over 9 isolates; 100 coverage values (99 at 1.000) |
| 6a | `cohort_variants` | ran | 73,752 sites |
| 7 | `pangenome` | **refused** | `ToolNotAvailableError`: `panaroo` not importable, `cd-hit` not on `PATH` |
| 8 | `recombination` | **not run** | `gubbins` cannot run on macOS arm64 |
| 9–16 | phylogeny → report | **not run** | downstream of 7 and 8 |

The run exits non-zero and `run_manifest.json` is not written, because the run
failed before the manifest step. Per-stage status is in the per-stage log only.

## Why stages 7–16 are blocked

Platform, not code. Full evidence in
[`environment-arm64.md`](environment-arm64.md) §3.

- **`gubbins` (stage 8)** — every conda `osx-arm64` build is Python-3.10-only and
  this codebase is 3.11.16. Installed in a separate 3.10 env it still crashes:
  `run_gubbins.py` exits **139 (SIGSEGV)**, and `nm -u` reports `_gzopen`,
  `_gzread`, `_gzclose` undefined in `libgubbins.0.dylib` with no `gz*` exported
  by `libSystem`. Input-independent, reproduces in a pristine environment. No
  debugger backtrace is stored, so no claim is made about the faulting
  instruction.
- **`panaroo` (stage 7)** — pip-installable on arm64 but not via conda, and it
  was absent from the run environment. `cd-hit`, a separate dependency, was also
  missing.
- **Stages 9–16 consume stage 8**, so the laptop cannot complete a REAL run.

## Open defects and pending decisions

1. **n = 9, not 10, in the variant tables.** `bcftools mpileup` was SIGKILLed
   (`returncode=-9`) on one isolate; it stayed a cohort member, which is correct
   and disclosed, but any variant statistic is n = 9. **Cause unproven** — not
   reproduced at realistic geometry (94 MB RSS).
2. **`max_depth` = 250, inherited not chosen.** It is the only mpileup RSS lever
   and it subsamples reads, so changing it moves `DP`/`AD`. **Pending science
   decision.** Output size was bounded instead with `-O z`, verified to leave the
   callable variant set byte-identical.
3. **Four reporting fixes specified, not implemented** — chiefly carrying
   `call_status` and the exit signal beside the cohort denominator.
4. **`--max-BP`/`-b` does not exist** in the pinned bcftools 1.23.1. Verified
   absent from the usage block and rejected on invocation. Not used.
5. **`tools_detected` is mostly `UNKNOWN`** by design — probing is on-demand so no
   unneeded tool is executed. `config/references.tsv` is the only version record
   and is largely unpinned.

## Bakta / Panaroo guards

A REAL run must show zero tool executions. Three independent guards, all required:

1. **Config** — `annotation.reuse_tool_output: "require"`; stage 2 imports
   curated output and refuses when it is absent.
2. **Tripwire** — executable `bakta` / `panaroo` shims on `PATH` that log their
   argv and exit non-zero. Non-vacuity proven by firing them before the run.
   **Final `calls.log` must be 0 bytes.** It was 64 bytes / 2 lines before the
   on-demand-probe fix, and is 0 bytes after.
3. **Watcher** — background `ps` sampler. Honest limitation: a sampler **cannot
   prove a negative** for sub-second processes. It caught none because the two
   real invocations each finished inside one sampling interval. It widens
   coverage (an absolute-path invocation that bypasses the shim); it does not
   decide the question.

## Dashboard (read-only results UI)

A read-only monitoring dashboard lives in [`dashboard/`](../dashboard/). It reads
a results root or a delivery bundle, writes only its own state directory
(`~/.pa_dashboard/`), makes no network calls, and binds `127.0.0.1` by default.
Start it with:

```bash
scripts/run_dashboard.sh --results <results-dir> --port 8765
```

- **Acceptance: 13 of 13** items (UI-1 … UI-13). The stage view shows true
  per-stage states and named reasons; the launcher is disabled unless the server
  is started with `--allow-launch`, defaults to a dry run, and requires the exact
  typed phrase `run real samples` for a REAL run.
- **Real data limited to stages 1–6a.** Against the real 10-isolate bundle the
  dashboard shows stages 1–6/6a `completed`, stage 7 `failed`
  (`ToolNotAvailableError`), and stages 8–16 `not_run`. Tree, imipenem SIR and
  lineage read **`not produced`** with the reason, never a zero or a blank.
- **Never browser-painted.** The pages were exercised through the HTTP API and a
  `node --test` route-mount suite; no browser or screenshot pass was run, so CSS
  layout and paint are asserted only structurally, not visually.

## Test suite

| Invocation | Result |
|---|---|
| `pytest -q` | 3343 passed, 65 skipped, 0 failed |
| `pytest tests dashboard/tests -q` | 3434 passed, 68 skipped, 0 failed |
| `pytest dashboard/tests -q` | 112 passed, 3 skipped, 0 failed |
| `PA_FIXTURES_LARGE=1 pytest dashboard/tests -q` | 114 passed, 1 skipped, 0 failed |
| `node --test dashboard/tests/js/` | 32 passed, 0 failed |

The 3 dashboard skips in default mode are the UI-D5/perf gates, which need the
900-isolate set (`PA_FIXTURES_LARGE=1`); in large mode those run and the one skip
is the guard that asserts the large set is *not* built in a normal run.

Skip census (68): 21 real clinical file absent by design, 10 `gubbins`, 20
opt-in `pyseer` (`PAPIPELINE_TEST_PYSEER=1`), 12 empty-parametrize tests that
assert nothing in a default run, remainder `snakemake`/STUB-mode and absent
fixtures.