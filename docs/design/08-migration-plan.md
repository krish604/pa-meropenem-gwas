# Design 8 — Migration plan

Companion to [1 architecture](01-architecture.md).

Rule for every step: **the full existing test suite runs after each one, and
must stay green.** Baseline today: **641 passing**.

---

## Step 0 — Documentation (this milestone)

Nine design documents plus the reusable-module inventory. No code touched.
`git status` clean except `docs/design/`.

---

## Step 1 — Packaging

**Problem.** There is no `pyproject.toml`, `setup.py` or `setup.cfg` anywhere.
The package is not installable; every script and `tests/conftest.py` repeats
a `sys.path.insert(PIPELINE_ROOT)`. There is no console-script entry point.

**Change.** Add `pyproject.toml`:

- package `papipeline` and `papipeline_app`
- dependencies already in `environment/environment.yml`
- `[project.scripts] papipeline = "papipeline_app.cli.main:main"`
- keep `scripts/` importable so the one test that does
  `import scripts.run_annotation_sweep` keeps working

**Also:** update the `sys.path` hacks in `scripts/common/*.py` and
`tests/conftest.py` to fall back gracefully when the package is installed,
rather than removing them outright — a source checkout must still work.

**Gate.** 641 pass. `pip install -e .` then `python -c "import papipeline"`
from outside the repo.

**Risk.** Low. Nothing scientific changes.

---

## Step 2 — Persisted stage state *(the important one)*

**Problem.** Six of sixteen stages take their inputs as in-memory
arguments, so a process can only run one stage by re-running its
prerequisites. `Snakefile:14-29` documents this as the reason the DAG was
never exercised. It is also the reason resume is impossible.

**The good news:** the other ten stages already read from disk, and
`run.py` already *writes* `01_validation.tsv` … `15_master_table.tsv`. The
artefacts exist. Only the reading side is missing.

**Change.** For each of the six, add a loader beside the existing `run()`
and have `run()` use it when the argument is not supplied:

| Stage | Add | Reads |
|---|---|---|
| `mechanisms` | `load_amr_from_disk`, `load_regulators_from_disk`, `load_sv_from_disk`, `load_gene_names_from_disk` | `04_amr.tsv`, `06_regulators.tsv`, `07_structural_variants.tsv`, `02_annotation_summary.tsv` |
| `gwas` | reuse `load_features`; add `load_phenotype_calls_from_disk` | `12_gwas_features.tsv`, `11_phenotype.tsv` |
| `convergence` | `load_names_from_disk(amr, regulators)`, `load_lineages_from_disk` | `04_amr.tsv`, `06_regulators.tsv`, `10_phylogeny.tsv` |
| `cooccurrence` | `load_names_from_disk(amr, regulators, mechanisms)` | `04_amr.tsv`, `06_regulators.tsv`, `05_mechanisms.tsv` |
| `integration` | `load_all_from_disk` | the eleven stage TSVs |
| `reporting` | `build_context_from_disk` | the stage TSVs + `run_manifest.json` |

**Signature rule:** every argument becomes `Optional`, defaulting to the
disk loader. **No existing call site breaks, and no existing test changes.**
`run.py` keeps passing in-memory objects and is therefore untouched.

**Gate.** 641 pass. Plus new tests: for each stage, run it *alone* from
disk-only inputs and assert its output equals the output from the
in-process path. This is the test that proves the DAG can be incremental.

**Risk.** Medium — this touches the scientific engine's signatures. Mitigated
by the Optional-default rule and by the equality tests.

---

## Step 3 — One per-stage entry point

**Problem.** 18 wrappers in `scripts/<stage>/run_*.py` in two
near-identical templates: 10 call `stage.run()` directly (no prerequisite
expansion, so they are only correct for disk-reading stages) and 5 are
near-verbatim copies of `run_stage.py` with the stage name hard-coded.
Plus 6 ad-hoc top-level drivers totalling 2,032 lines of pilot-specific
code with hard-coded paths.

**Change.**

- Promote `scripts/common/run_stage.py` to be *the* per-stage entry point:
  add `--project`, `--out-root`, `--resume`, and machine-readable `--json`.
- Convert `scripts/<stage>/run_*.py` into three-line shims that call it, so
  existing invocations keep working. **Do not delete them yet** — the
  Snakefile and any muscle memory depend on the paths.
- Move the genuinely reusable logic out of the ad-hoc drivers:
  - `scripts/run_annotation_sweep.py` → `papipeline_app/resources/bakta.py`
  - `scripts/preflight.py` → `papipeline_app/tools/preflight.py`
  - `scripts/consolidate_results.py` → `papipeline_app/results/`
  - `scripts/pyseer_compat.py` → stays, but referenced by
    `papipeline_app/engine/` rather than by a hard-coded
    `parents[2]/"scripts"/...` path

**Gate.** 641 pass. `test_run_stage_cli.py` covers `main(argv=[...])` for
every stage, both `--mode` values, and the exit-code contract from design 5.

**Risk.** Low-medium. Shims preserve every existing invocation.

---

## Step 4 — Store, roots, and the project CLI

**Change.** `papipeline_app/store/*` per design 3; `paths.py`;
`backend.py` verbs `project.*`, `input.*`, `config.effective`,
`workflow.stages`; the CLI in `papipeline_app/cli/`.

**Gate.** 641 pass + new `tests/app/test_store.py`,
`test_projects.py`, `test_validation.py`, `test_backend_parity.py`.

---

## Step 5 — Snakemake becomes authoritative

**Blocked on step 2.** Once stages are disk-loadable:

- rewrite `Snakefile` as a thin include-only file;
- move the 20 rules into `rules/*.smk` with real `resources:`;
- generate prerequisite edges from `papipeline.run.PREREQUISITES`;
- add `profiles/{macbook,workstation,hpc}.yaml`;
- **add the missing `gene_families.tsv` to `CONFIG_FILES`** — a real
  dependency-tracking bug today;
- switch to `configfile:` so config changes invalidate correctly.

**Gate.** 641 pass + `test_snakemake_integration.py`:
DAG builds in TEST; `--dry-run` refuses REAL; a full `-j1` TEST run matches
`run_pipeline.py --mode TEST` byte-for-byte; touching
`config/gene_families.tsv` invalidates the AMR and mechanism rules; killing
mid-annotation and re-running schedules only pending genomes.

**This is the first time the DAG is ever executed.** Until those tests are
green, the Snakefile is not authoritative and the app says so in the UI.

---

## Step 6 — Resource-aware Bakta executor

`papipeline_app/resources/{detect,profiles,recommend}.py` and `bakta.py`,
lifting the logic out of `scripts/`. Shared DB, checkpointing, per-genome
logs and timing, retries with partial-output clearing, pruning, and
`benchmark` recommendations.

**Gate.** 641 pass + `test_resources.py`, `test_bakta_scheduler.py`,
`test_resume.py`, `test_checkpointing.py`, `test_failed_jobs.py`.

---

## Step 7 — Results registry and reports

`results/registry.py` discovers every declared artefact type by stage and
records type, path, format, row/column counts and a header fingerprint.
`reports/builder.py` composes the four report kinds **from the registry**,
so a report cannot silently omit a stage that produced output.

Prefer the engine's `stages/reporting.py` for the full report and build the
other three as views over the same registry. **Do not fork
`stages/reporting.py`** — it already carries the TEST/REAL banner logic and
the claim vocabulary, which are scientific rules.

**Gate.** 641 pass + `test_results_registry.py`, `test_reports.py`
(including: a report for a run with a FAILED stage says so and does not
present the partial result as complete).

---

## Step 8 — Provenance, and pinning `references.tsv`

`provenance/capture.py` records inputs + checksums, tool versions, database
versions, effective config, resolved resources, pipeline version and git
commit, timestamps, and every executed command.

Then: **fill in the real versions in `config/references.tsv`.** Today all 11
rows are `UNPINNED` or `na`, which is a standing violation of scientific
rule 8 ("the database version must always be recorded"). This is a
scientific-completeness task, not plumbing, and it is the last thing that
should be marked done.

Note `environment/environment.yml` disagrees with reality (design 1, D3):
it pins `bakta=1.9.3`, `amrfinder=3.12.0`, `pyseer=0.4.5`, and
`panaroo=4.0.0`/`snp-sites`/`iqtree2` which have no `osx-arm64` build.
Add `envs/constraints-osx-arm64.yml` and `envs/constraints-linux.yml`, and
make provenance record what was *actually used*, not what a manifest
requested.

**Gate.** 641 pass + `test_provenance.py`: a run's provenance alone is
sufficient to explain every number in its report.

---

## Step 9 — HTTP API

`papipeline_app/web/`, per design 4. Backend already exists, so this is
thin. Framework is open decision **D4**.

**Gate.** 641 pass + `test_web_api.py`, and a parity test asserting the HTTP
route set and the CLI verb set are the same registry.

---

## Step 10 — Web UI, one screen at a time

Order: Dashboard → New Project → Input Selection → Metadata/Phenotype →
Pipeline Designer → Resource Configuration → Run Monitor → Logs → Results
Explorer → Report Generator → Tool/Database Manager → Project Settings.

Each screen lands against a backend verb that already exists and is already
tested. **No screen gets its own logic.** If a screen needs something the
backend cannot do, the backend grows a verb and the CLI grows it too.

**Gate.** 641 pass + parity test green at every step.

---

## Retiring `papipeline/pilot/`

Not deleted — migrated, then deprecated, one concern at a time. It is the
only implementation of several things the app needs.

| Pilot asset | Destination | Note |
|---|---|---|
| `pilot/amr_detect.py` | `papipeline_app/engine/amrfinder.py` | The v4 parser. `stages/amr.py::parse_amrfinder_stdout` is v3-era and its `AmrFinderPlusAdapter` passes `--input`, which 4.x renamed to `--nucleotide`. The pilot module documents the divergence and is the correct implementation. |
| `pilot/cohort.py` | `papipeline_app/projects/validation.py` | Assembly-integrity assessment. Useful: 835 assemblies yielded 65 analysable, and the reason per file is worth keeping. |
| `pilot/pdc_fields.py` | `papipeline_app/inputs/ast.py` | PDC free-text AST/genotype parsing. |
| `run_annotation_sweep.py` | `resources/bakta.py` | Step 6. |
| `pilot/figures.py`, `pptx.py` | `papipeline_app/reports/` | The engine prepares figure *data* as JSON in `viz.py` and explicitly has no drawing layer; the pilot's matplotlib and python-pptx code fills that gap. |
| `pilot/analysis.py`, `mechanisms.py`, `report.py`, `orchestrator.py` | — | **Superseded.** `stages/{mechanisms,cooccurrence,reporting}.py` are the tested 16-stage versions; the pilot's co-occurrence in particular has no statistics. |

**Rule:** the pilot is never deleted while a test depends on it.
`tests/unit/test_pilot.py` has 92 tests. Those are migrated to the new
homes, and only then does the module go. Deleting tests is out of the
question; migrating them is the work.

**Recommendation on decision D5:** migrate and retire. The pilot cannot
produce stages 6, 7 or 8 at all, so keeping it guarantees two divergent
answers to "what does this pipeline do".

---

## Risk register

| Risk | Likelihood | Mitigation |
|---|---|---|
| Step 2 changes scientific behaviour | Medium | Optional-argument defaulting only; equality tests old-vs-new path; no scientific definition touched |
| Snakefile stays untested and becomes a second, wrong scheduler | Medium | Step 5 is gated on byte-identical TEST output; app reports which path ran |
| `pilot/` and `stages/` diverge further during migration | Medium | Migrate one asset per commit; never modify a stage to suit the app |
| REAL mode produces unpinned-version results | High today | Step 8 pins `references.tsv`; preflight warns; decision D2 decides refuse-vs-warn |
| The UI is built before the backend can tell the truth | Medium | Milestones order UI last; Run Monitor depends on step 2's persisted state |
| `results/` git-ignored hides real outputs | Certain today | Project directories live outside the repo; decision D6 places the DB under `PAPIPELINE_HOME` |
