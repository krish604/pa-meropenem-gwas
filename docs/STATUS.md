# Status

Last real-data run: round 12, `integration-r3` at `3f1c633`.

## Summary

| | |
|---|---|
| Stages implemented | 16 of 16 |
| Stages verified on **real** data | 1, 2, 3, 4, 4a, 5, 6, 6a (8 of 16) |
| Stages never run on real data | 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 |
| Tests | 3764 collected — **33 failed, 3519 passed, 151 skipped, 61 errors** (all pre-existing; see [Test suite](#test-suite)) |
| Cohort used | 10 isolates (`local/smoke_isolates.txt`) |
| `bakta` / `panaroo` executions in a pipeline run | **0** — Bakta is an input, never re-run. `panaroo` has since been executed **once, outside any pipeline run**, against real Bakta GFF3 for this cohort (exit 0, 10,019 families); no REAL run has repeated stage 7 since |

**This is not a complete 16-stage result and nothing here is a biological
finding.** n = 10 is grossly underpowered.

## 2026-10-08 — meropenem GWAS build

Six build packages, seven commits: `1a50b60` + `552e56c` (meropenem config +
cohort gate), `098d3e2` (layer encoder), `248b4f3` (kinship + pyseer adapter),
`552ac7a` (downstream stats), `5ffcdda` (linux runner), `a379cf6` (serial
integration).

**Nothing ran on real data in this build.** No `data/`, no `db/`, no
`PDC_essential.tsv` was read, and no REAL run was started. The round-12 sections
below remain the only real-data evidence this repository has. Also unchanged:
`project.primary_antibiotic` is still `imipenem`, and the ~500-isolate
meropenem cohort has not been assembled or analysed — so no cohort-gate
threshold has ever been exercised against a real cohort.

### Weakest true state, per package

`built` = implemented and unit-tested, nothing calls it. `wired` = reachable
from the path a run takes. `verified` = wired and exercised end to end on
realistic input. **Nothing below is verified on real data.**

| Package | Weakest true state | Note |
|---|---|---|
| Meropenem enablement | **wired and verified (config only)** | `science.yaml` `antibiotics:` = `imipenem, meropenem`; `antibiotics.tsv` row enabled; all 22 `mechanisms.tsv` rows read `imipenem,meropenem` (21 rewritten, `oprD` already had both). **No analysis code changed.** |
| Cohort gate (`papipeline/cohort_gate.py`) | **wired, not verified on a real cohort** | Evaluated before stage 1 in every mode, recorded on `RunResult` and in `run_manifest.json`. `cohort_gate.intermediate_policy` (`exclude`), `cohort_gate.min_resistant` (`100`) — **enforced in REAL only** (the 20-sample TEST fixture is 7 R of 20, below it by construction); counts recorded in every mode so the gate is never silent. |
| Layer encoding (`papipeline/layers/`) | **wired** | Called at the pangenome stage in REAL and TEST, never STUB; a TEST run writes the 8 layer files. Keys `layers.min_carriers: 5`, `layers.max_prevalence: 0.98`. The **unitig adapter is built and deliberately unwired**: `unitig-caller` is absent and its flags were never verified with `--help`, so REAL raises `UnverifiedFlagsError` rather than guess. |
| Kinship (`papipeline/gwas_real/kinship.py`) | **built and unit-verified** | 13 tests, including the pre-existing spec test `tests/unit/test_gwas_real_kinship.py`, which passes with zero edits. Reuses `stages.similarity.patristic_distances` and Gower-centres the squared distances. `PREREQUISITES["gwas"]` now includes `"similarity"`. |
| pyseer two-pass adapter (`papipeline/gwas_real/adapter.py`) | **built, NOT wired** | Three input resolutions undecided; wiring it would risk a wrong `variants_path` silently narrowing which variant families are tested. Separately, `REAL_REFUSING_STAGES["gwas"]` says "no REAL engine" while `_build_gwas_engine` injects `PyseerEngine` in REAL when pyseer resolves — conditional on pyseer being absent, and **not re-verified by this build**. |
| Downstream (`papipeline/downstream/`) | **built; runner wired for REAL only** | Seven steps in the pinned order, positive-control gate first, writing `reports/downstream_evidence.tsv`. **Not exercised in TEST**: `run_control_gate` raises on the committed fixtures (`any_MBL` and `oprD_burden` both absent from stage 12's TEST results) and `build_report` re-runs the gate with no disable flag. Steps 2–7 execute correctly on TEST fixtures once the gate is fed a recovering baseline. Two tests pin the mode gate from both ends. |
| Linux runner (`scripts/linux/`, `config/machines/linux.yaml`, `docs/LINUX_RUN.md`, `docs/PIPELINE_FLOWCHART.md`) | **built, never executed on a Linux host** | The environment **fails to solve as pinned**: `micromamba create --dry-run --platform linux-64 -f environment/environment-linux.yml` → exit 1, `gubbins =3.4.1 * does not exist` (published linux-64 builds are py310); `samtools=0.1.19` needs `openssl <=1.1.1`. Core minus those pins solves (357 packages). **No pin was changed** — a re-solve is open work. `shellcheck` is not installed; substitute evidence is `bash -n`. The overlay sets `allow_real_mode: false`. |

### Verification run, 2026-10-08

| Check | Result |
|---|---|
| `python3 -m pytest tests --collect-only -q` | **3764 collected, 0 errors** |
| `python3 -m pytest tests -q --tb=no` | **33 failed, 3519 passed, 151 skipped, 61 errors** |
| non-passing node IDs vs `.build/BASELINE_PREEXISTING_FAILURES.txt` (106 ids, recorded **before** this build) | 94 ids, **0 new failures**, **12 previously-failing STUB tests now pass** |
| `python3 -m pytest dashboard/tests -q` | 112 passed, 3 skipped (unchanged) |
| `node --test dashboard/tests/js/` | 32 passed, 0 failed (unchanged) |
| `python3 scripts/common/run_pipeline.py --mode TEST` / `--mode STUB` | exit 0 both |
| `snakemake -n` with `machine=config/machines/laptop.yaml` and `…/linux.yaml`, `mode=STUB` | exit 0 both |

Open after this build: the `project.primary_antibiotic` flip **plus** its
fixture; the three `gwas_real` input resolutions; `unitig-caller` flags; the
Linux environment re-solve; and downstream's REAL-only path. Details and the
full open list are in `TASKS.md` Phase 8.

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
| 7 | `pangenome` | **refused** | `ToolNotAvailableError`: `panaroo` not importable, `cd-hit` not on `PATH`. **Record of round 12 only** — `panaroo` has since been installed from source (pinned) with `cd-hit` from conda and verified end-to-end on real Bakta GFF3 for this same cohort (exit 0, 10,019 families, 4,828 core); **no REAL run has been repeated since**, so nothing has been re-measured |
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
  missing. **Superseded for the tool, not for the result:** panaroo has since
  been installed from source (pinned, `cd-hit` from conda) and verified
  end-to-end on real Bakta GFF3 for this cohort, but **no REAL run has been
  repeated since**, so stage 7's round-12 refusal above still stands as the
  last thing a pipeline run measured.
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

Corrected 2026-10-08. **`pytest tests -q` (and a bare `pytest`, which is the
same invocation — `pytest.ini` sets `testpaths = tests`) currently reports
33 failed, 3519 passed, 151 skipped, 61 errors.** The suite is **not** green
and is not reported as green.

| Invocation | Result (2026-10-08) |
|---|---|
| `pytest tests --collect-only -q` | 3764 collected, **0 errors** |
| `pytest -q` / `pytest tests -q` | **33 failed**, 3519 passed, 151 skipped, **61 errors** |
| `pytest tests dashboard/tests -q` | 3879 collected — **33 failed**, 3631 passed, 154 skipped, **61 errors** |
| `pytest dashboard/tests -q` | 112 passed, 3 skipped, 0 failed (**unchanged**) |
| `PA_FIXTURES_LARGE=1 pytest dashboard/tests -q` | 114 passed, 1 skipped, 0 failed (**unchanged**) |
| `node --test dashboard/tests/js/` | 32 passed, 0 failed (**unchanged**) |

The dashboard suites are unchanged by this build. The `tests/` numbers above
are **all pre-existing**: every one of the 94 non-passing node IDs is listed in
`.build/BASELINE_PREEXISTING_FAILURES.txt`, which was recorded from the
untouched repository *before* the meropenem build — **0 new failures, and 12
previously-failing STUB tests now pass** (the STUB data-root fix). They fail
because they need `db/` (`db/reference/…`, `db/smoke_genomes`, the Bakta and
AMRFinderPlus databases) or `data/` — inputs that are **absent by design**,
since databases and genomes are never committed — so they fail on any clean
checkout. They are not evidence of a regression and they are not evidence of
health.

Skip census (151): 117 need uncommitted inputs (`db/`, `PDC_essential.tsv`, the
source-built `tools/gubbins`), 21 are the opt-in `pyseer` binary suite
(`PAPIPELINE_TEST_PYSEER=1`), 12 are empty-parametrize tests that assert
nothing in a default run, 1 needs `python-pptx`.

The 3 dashboard skips in default mode are the UI-D5/perf gates, which need the
900-isolate set (`PA_FIXTURES_LARGE=1`); in large mode those run and the one skip
is the guard that asserts the large set is *not* built in a normal run.