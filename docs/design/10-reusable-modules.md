# Design 10 — Reusable module inventory

Companion to [1 architecture](01-architecture.md). Produced by reading the
repository, not by inferring from the README.

The brief asks which existing modules can be reused. This is the answer,
with the four categories stated explicitly so nothing is ambiguous.

| Category | Count | Meaning |
|---|---|---|
| **REUSE AS-IS** | 24 | Use directly. Do not edit. |
| **REUSE, EXTEND** | 6 | Correct shape, missing a capability the app needs. |
| **REUSE, MIGRATE** | 5 | Right behaviour, wrong location — move it in. |
| **SUPERSEDED** | 7 | A better implementation already exists. Retire, do not delete. |

---

## 1. REUSE AS-IS (24)

These are the assets. Most of the value of this repository is here, and the
brief's "do not rewrite existing modules" is easy to honour because there
is very little worth rewriting.

### Contracts — the highest-value assets

| Module | Why it is reused unchanged |
|---|---|
| `papipeline/models.py` (531 lines) | 9 enums + 14 frozen dataclasses, every one with `to_row()`. `ClaimStatus` having **no `CAUSAL` member** is a scientific guarantee encoded in the type system. The app's store must serialise these, never redefine them. |
| `papipeline/manifest.py` | `SampleManifest` frozen dataclass with `require()`, `index()`, `with_source()`; `validate_sample_id` encodes the real identifier rules. Used verbatim by project creation. |
| `papipeline/io/tsv.py` | The strictness rules *are* the data contract: missing sentinels → `None`, more-fields-than-header is an error, short rows pad, header-only is an error. Reused by every reader in the app. |
| `papipeline/io/fasta.py` | Streaming reader, `assembly_stats` (N50/N90/GC in one pass), `file_integrity` (returns a dict, never raises — exactly what input validation needs). |
| `papipeline/errors.py` | 14-class hierarchy rooted at `PipelineError(message, **context)`, context rendered into the message. The app's `AppError` mirrors it deliberately, so `str(exc)` is user-facing on both sides. |
| `papipeline/logging_utils.py` | Idempotent `configure_logging`, single `papipeline` namespace. |

### Scientific engine — all 16 stages

`validation`, `annotation` (parsers), `mlst`, `amr` (table path),
`mechanisms`, `regulators`, `sv`, `virulence`, `pangenome`, `phylogeny`,
`phenotype`, `gwas` (engines + parsers), `convergence`, `cooccurrence`,
`integration`, `reporting`.

Ten already read from disk and are independently runnable **today**. The
six that take in-memory arguments gain a loader in migration step 2 —
their scientific logic is untouched.

Two are worth singling out:

- **`stages/mechanisms.py`** is the scientific core. It enforces the
  claim-ceiling rule (`knowledge.ceiling_for`), keeps `oprD` presence at
  `locus_intact`/`DETECTED` rather than implying susceptibility, and refuses
  to emit a call for a `not_assessable` SV. Touch nothing.
- **`stages/phenotype.py`** has no `run()`; the entry point is
  `load_phenotype(config, phenotype_dir, antibiotic, ...)`. It already
  hard-errors on MIC + `ND` and on an unconfigured antibiotic, and
  `tests/unit/test_phenotype.py::TestNoFabricatedQuantitation` guards
  R/I/S → MIC.

### Configuration and knowledge

| Module | Reused for |
|---|---|
| `config/loader.py` (592 lines) | `PipelineConfig` + 4 sub-configs + 4 knowledge-table specs. `load_config()` takes an optional path; `run_pipeline(config=...)` accepts a preloaded config, so the app parses once. `tool_output_root(mode)` is what keeps a run from reading its own outputs as inputs. |
| `config/config.yaml` | Unchanged. `analysis` flags drive which stages run; `qc` thresholds drive validation. |
| `config/{mechanisms,regulators,antibiotics,references}.tsv` | Unchanged. Editable knowledge, per the brief. |
| `config/gene_families.tsv` | Unchanged, and it must be added to the Snakefile's `CONFIG_FILES` — currently missing. |
| `papipeline/knowledge/` | `ceiling_for()` is the enforcement point for claim strength. The app must not bypass it. |

### Orchestration and support

| Module | Reused for |
|---|---|
| `papipeline/run.py` | `STAGE_ORDER`, `EXECUTION_ORDER`, `PREREQUISITES`, `resolve_prerequisites`, `resolve_mode`, `run_pipeline`, `build_provenance`, `_write_run_manifest`. The stage graph is **data the app reads**, not logic the app reimplements. |
| `papipeline/viz.py` | 13 figure-data builders, pure data in / JSON out, no plotting. The GUI renders; this prepares. |
| `papipeline/testing/synthetic.py` (1,003 lines) | Seeded, byte-reproducible fixture generator. Every app test uses it instead of hand-written fixtures. |
| `papipeline/cli.py` | `describe`/`summarise`/`write_rows` helpers for thin wrappers. |
| `scripts/common/run_stage.py` | Promoted to the single per-stage entry point (step 3). |
| `scripts/common/run_pipeline.py` | The full-run entry point; `print_plan()` is a good model for the preflight surface. |
| `scripts/common/write_provenance.py` | `STAGE_TOOL_REQUIREMENTS` is the natural source for the tool registry. |
| `environment/versions.sh` | Already a provenance recorder that "only reads versions, never installs". Reuse, don't rewrite. |

---

## 2. REUSE, EXTEND (6)

Correct shape, one capability missing.

| Module | Gap | Change |
|---|---|---|
| `papipeline/adapters/external.py` (265 lines) | `ToolStatus`, `detect_all()`, `VERSION_FLAGS` for 15 tools — the right home for the registry. Missing: **version-floor comparison** and **resolution that ignores PATH order**. | Add `resolve(name, min_version) -> ToolStatus`, preferring any candidate that meets the floor. Concrete reason: this machine's env ships samtools/bcftools 0.1.19 (2014) that shadow working Homebrew 1.24 and **fail silently** on modern SAM. |
| `papipeline/stages/annotation.py` | `run()` raises `NotImplementedError` when `mode is RunMode.REAL` — *"REAL-mode annotation requires the Bakta adapter and a pinned database; wire it in before the analysis phase"*. Bakta is invoked from `scripts/run_annotation_sweep.py` instead, so the engine cannot annotate. | Implement the REAL path by delegating to the Bakta executor (design 7) and persisting to the paths `load_from_intermediate` already expects. No change to the parsers. |
| `papipeline/stages/amr.py` | `AmrFinderPlusAdapter` is v3-era: passes `--input`, which 4.x renamed to `--nucleotide`. `parse_amrfinder_stdout` expects the old column order. | Replace the adapter with the v4 implementation from `pilot/amr_detect.py` (see §4). The table-reading path (`load_amr_table`) is correct and stays. |
| `papipeline/run.py` | `resolve_mode()` hard-blocks REAL via `runtime.allow_real_mode: false`. Deliberate and correct for an unattended run; wrong for an interactive user who has explicitly created a project. | Open decision **D1**. Recommendation: add an app-level acknowledgement rather than removing the gate. |
| `papipeline/stages/gwas.py` | `PyseerEngine` shells out to `sys.executable scripts/pyseer_compat.py` via a hard-coded `parents[2]/"scripts"/...` path. Also `cooccurrence.benjamini_hochberg` is a duplicated copy of the one in `gwas.py`. | Move the compat wrapper's location into the package; de-duplicate the BH implementation. Both are mechanical. |
| `workflow/Snakefile` | 20 rules, never executed, `rules/` empty, no profiles, no `resources:`, no `configfile:`, and `CONFIG_FILES` omits `gene_families.tsv`. | Rewrite per design 6 — **gated on migration step 2**. |

---

## 3. REUSE, MIGRATE (5)

Right behaviour, currently in the wrong place. Move in; do not rewrite.

| From | To | Why |
|---|---|---|
| `scripts/run_annotation_sweep.py` (498 lines) | `papipeline_app/resources/bakta.py` | Already has: shared DB resolution, `database_ready()`, checkpointing via `is_complete()`, `prune_outputs()`, a `run_one` documented "never raises", per-genome logs and timing, state file, `--prune/--no-prune`. This **is** the resource-aware Bakta executor the brief asks for; it is in the wrong place and has no retry. |
| `scripts/preflight.py` (307 lines) | `papipeline_app/tools/preflight.py` | Already does tool version floors, database readiness, config load, cohort checks, disk headroom, and a GO/NO-GO. Hard-codes a `353 s/genome, 2 jobs × 4 threads` ETA that must come from `benchmark`. |
| `scripts/consolidate_results.py` (320 lines) | `papipeline_app/results/` | Correct discipline already: reads only files the pipeline wrote, leaves `"."` where a stage produced nothing, never imputes. That is the results registry's contract. |
| `scripts/pyseer_compat.py` (67 lines) | `papipeline_app/engine/` | Namespace shims so pyseer 1.1.2 runs under scipy 1.17 / numpy 2.4 / pandas 3.0 / statsmodels 0.15. Its own docstring asserts no fitted value is altered — worth preserving verbatim. |
| `papipeline/pilot/figures.py`, `pilot/pptx.py` | `papipeline_app/reports/` | The engine has `viz.py` (13 figure-*data* builders) and README says "the drawing layer is not implemented". The pilot's matplotlib and python-pptx code is that missing layer. |

---

## 4. SUPERSEDED (7)

A tested 16-stage implementation already exists. Retire the duplicate; do
not delete the file until its tests move.

| Module | Superseded by | Evidence |
|---|---|---|
| `pilot/analysis.py` (492) | `stages/{cooccurrence,phenotype}.py` | The pilot's gene-pair table is raw counts with **no statistics and no multiple-testing correction**; `stages/cooccurrence.py` does Fisher + BH and stamps `interpretation_limit="association_only_not_causal"`. Strictly worse. |
| `pilot/mechanisms.py` (332) | `stages/mechanisms.py` + `config/gene_families.tsv` | Two competing mechanism vocabularies. The stage version is config-driven and has the claim-ceiling enforcement. |
| `pilot/report.py` (254) | `stages/reporting.py` | The stage version has the TEST/REAL banner logic and reuses `io.tsv.write_tsv`. |
| `pilot/orchestrator.py` (478) | `papipeline/run.py` | A second orchestrator for a 7-step subset. Its only stage import is `stages/validation.py`. |
| `pilot/amr_detect.py` — *partly* | `stages/amr.py` | See below: the **parser** is better in the pilot; the orchestration is not. |
| `scripts/<stage>/run_*.py` (18 files) | `scripts/common/run_stage.py` | 10 in one template, 5 in another, 5 near-verbatim copies of `run_stage.py`. Converting to three-line shims removes the duplication without breaking any invocation. |
| `scripts/run_downstream_stages.py` (903) | `run_pipeline.py` + the stage disk loaders | An ad-hoc driver that re-implements stage ordering, has no `--config`/`--mode`/TEST-REAL gating, and hard-codes `results/pilot100/`. It is also where several real bugs were found and fixed — so its *lessons* are captured in the engine, and the file itself is temporary. |

### The `pilot/amr_detect.py` nuance

This one is genuinely two things:

- **Reuse:** the v4 report parser. Its module docstring states the
  divergence outright — *"AMRFinderPlus 4.x renamed `--input` to
  `--nucleotide` and re-ordered the report columns, so the v3-era parser in
  `papipeline.stages.amr` is not reused verbatim. This module produces the
  same `AmrDeterminant` records, so everything downstream of stage 4 is
  unchanged."* That is the correct implementation and it belongs in the
  engine as the stage-4 adapter.
- **Supersede:** the cohort sweep orchestration, which duplicates what
  Snakemake's per-genome rules do properly.

---

## 5. The single most important structural fact

Ten of sixteen stages **already read their input from a path** under the
intermediate root, and `run.py` **already writes** `01_validation.tsv` …
`15_master_table.tsv`:

| File-based today (10) | In-memory (6) |
|---|---|
| `validation` · `annotation` · `mlst` · `amr` · `regulators` | `mechanisms` |
| `sv` · `virulence` · `pangenome` · `phylogeny` · `phenotype` | `gwas` · `convergence` |
| | `cooccurrence` · `integration` · `reporting` |

So the artefacts for a resumable, incremental pipeline are **already being
written**. What is missing is only that the six in-memory stages do not read
their predecessors' artefacts back.

That makes the hardest requirement in the brief — resume, and a Snakemake
scheduler that can see between-stage edges — a **loader-addition**, not a
rewrite. It is migration step 2, it is the highest-value change in the
plan, and it is why this design does not propose touching a scientific
definition anywhere.

---

## 6. Not reusable, and why

| Thing | Why not |
|---|---|
| Genome assembly | Input is assembled FASTA. The brief forbids making assembly a stage. |
| `run_pdc.sh` | NCBI Datasets downloader hard-coded to `/home/thor/Downloads/...`. Data acquisition, not application. |
| `results/pilot100/` | Ad-hoc outputs from a bespoke driver. Useful as a fixture for L6 tests, not as an input to anything. |
| `environment/environment.yml` as an authority | It pins versions the machine does not have, three of which cannot be installed on `osx-arm64`. See decision D3. |
| The 92 tests in `test_pilot.py` | Not reusable *as they are* — but every one is **migrated**, not deleted, when `pilot/` retires. |
