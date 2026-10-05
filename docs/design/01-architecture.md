# Design 1 — Architecture

Status: **proposal**. Nothing here is implemented. First milestone of the
reusable-application work on `papipeline`.

Companion documents: [2 directory structure](02-directory-structure.md) ·
[3 data model](03-data-model.md) · [4 API design](04-api-design.md) ·
[5 CLI design](05-cli-design.md) · [6 Snakemake integration](06-snakemake-integration.md) ·
[7 Bakta optimisation](07-bakta-optimization.md) · [8 migration plan](08-migration-plan.md) ·
[9 test strategy](09-test-strategy.md) · [reusable module inventory](10-reusable-modules.md)

---

## 1. What exists today

Verified by reading the repository, not inferred from the README.

| Layer | Location | Size | State |
|---|---|---|---|
| Scientific engine | `papipeline/stages/` (16 modules) | ~4,200 lines | Implemented, **TEST mode only** |
| Orchestrator | `papipeline/run.py` | 1,136 lines | Implemented; in-memory state passing |
| Contracts | `papipeline/models.py` | 531 lines | 9 enums, 14 dataclasses |
| Pilot path | `papipeline/pilot/` (11 modules) | ~3,700 lines | **Parallel duplicate** of 5 stages |
| Workflow | `workflow/Snakefile` | 500 lines, 20 rules | **Never executed** |
| CLI wrappers | `scripts/common/` (4) + `scripts/<stage>/` (18) | ~2,000 lines | Working but 18 near-identical files |
| Ad-hoc drivers | `scripts/` top level (6) | 2,032 lines | Pilot-100 specific, hard-coded paths |
| Tests | `tests/` | **641 passing** | No Snakefile or CLI coverage |
| Knowledge tables | `config/*.tsv` | 6 files | `references.tsv` **entirely unpinned** |
| Docs | `docs/` | 5 files, 1,152 lines | Includes the 10 scientific rules |
| Packaging | — | — | **None. No `pyproject.toml` anywhere.** |

### Three findings that shape the whole design

**F1 — There are two pipelines, not one.**
`papipeline/stages/` is the real 16-stage engine. `papipeline/pilot/` is a
later parallel implementation that imports exactly **one** stage module
(`stages/validation.py`, for QC) and duplicates the rest:

| Stage | `stages/` | `pilot/` | Relationship |
|---|---|---|---|
| 1 validation | yes | reused | shared |
| 2 annotation | parsers only | external sweep script | partial |
| 4 amr | `amr.py` | `pilot/amr_detect.py` | **duplicate** (v3 vs v4 CLI) |
| 5 mechanisms | `mechanisms.py` | `pilot/mechanisms.py` | **duplicate** (different vocabulary) |
| 7 sv | `sv.py` | — | **missing from pilot** |
| 8 virulence | `virulence.py` | — | **missing from pilot** |
| 6 regulators | `regulators.py` | — | **missing from pilot** |
| 11 phenotype | `phenotype.py` | `pilot/analysis.py` | **duplicate** |
| 14 cooccurrence | `cooccurrence.py` | `pilot/analysis.py` | **duplicate** (pilot version has no statistics) |
| 15 integration | `integration.py` | bypassed by `scripts/consolidate_results.py` | **bypassed** |
| 16 reporting | `reporting.py` | `pilot/report.py` + `pptx.py` | **duplicate** |

The pilot is a dead end for a reusable application: it cannot produce
regulator, structural-variant or virulence output at all. **Decision: the
engine is `papipeline/stages/`. `papipeline/pilot/` is retired by migration,
not deleted** (see [8](08-migration-plan.md)).

**F2 — The engine is already ~60% file-based, which makes resume cheap.**
This is the most important structural fact. Stages that read their input
from a path under the intermediate root can each be run independently
today:

| File-based (10) | In-memory (6) |
|---|---|
| `validation` (FASTA paths) | `mechanisms` (amr + regulators + sv + genes) |
| `annotation` (`annotation/<id>.tsv`) | `gwas` (phenotype calls, lineages) |
| `mlst` (`mlst/mlst_results.tsv`) | `convergence` (amr, regulators, lineages) |
| `amr` (`amr/amr_determinants.tsv`) | `cooccurrence` (genes, variants, mechanisms, lineages) |
| `regulators` (`regulators/regulator_variants.tsv`) | `integration` (12 arguments) |
| `sv` (`structural_variants/structural_variants.tsv`) | `reporting` (`ReportContext`) |
| `virulence` (`virulence/virulence_factors.tsv`) | |
| `pangenome` (annotation dir) | |
| `phylogeny` (tree dir) | |
| `phenotype` (`<Antibiotic>_phenotype.tsv`) | |

`run.py` **already writes** `01_validation.tsv` … `15_master_table.tsv`.
So the artefacts exist. What is missing is only that the six in-memory
stages do not read their own predecessors' artefacts back. Adding a
`load_*_from_disk` sibling to each of those six converts the whole engine
into independently-resumable stages — **with no change to any scientific
definition.**

**F3 — The Snakefile's own documented limitation is the resume blocker.**
`workflow/Snakefile:14-29`, verbatim:

> The stage functions exchange in-memory structures, not files, so a
> process that runs only one stage cannot reconstruct the state of the
> stages before it. […] a full run through this Snakefile recomputes shared
> prerequisites more than once. […] it was not built here because snakemake
> is not installed in the build environment, so the DAG could not be
> exercised and an untested caching layer would be a liability rather than
> an asset.

That is F2 seen from the workflow side. Fixing F2 fixes this, and the DAG
becomes genuinely incremental. The comment also tells us the author already
identified per-stage state caching as the right fix.

---

## 2. Target architecture

```
┌──────────────────────────────────────────────────────────────────────┐
│  PRESENTATION          papipeline_app/web  ·  papipeline_app/cli      │
│  (two front doors, one backend — identical capabilities)             │
└───────────────────────────────┬──────────────────────────────────────┘
                                │  HTTP  /  in-process calls
┌───────────────────────────────▼──────────────────────────────────────┐
│  APPLICATION LAYER            papipeline_app/                        │
│                                                                       │
│   store/      SQLite: projects, samples, stages, jobs, runs, results │
│   engine/     adapter that drives the scientific engine               │
│   resources/  machine detection, resource profiles, Bakta executor    │
│   tools/      tool + database registry (discovery, versions)         │
│   results/    results registry (what outputs exist, typed)           │
│   provenance/ run manifests, checksums, config snapshots             │
│   reports/    report composition over the registry                    │
└───────────────────────────────┬──────────────────────────────────────┘
                                │  one interface:  run_stage(config, stage, roots) -> StageResult
┌───────────────────────────────▼──────────────────────────────────────┐
│  SCHEDULER                   Snakemake — authoritative for           │
│                              resources, dependencies, retries        │
│                              profiles: macbook · workstation · hpc    │
└───────────────────────────────┬──────────────────────────────────────┘
                                │  subprocess: scripts/common/run_stage.py --stage X
┌───────────────────────────────▼──────────────────────────────────────┐
│  SCIENTIFIC ENGINE            papipeline/stages/*.py  (UNCHANGED)     │
│  contracts: papipeline/models.py, config/, io/  (UNCHANGED)          │
│  orchestrator: papipeline/run.py                (UNCHANGED)          │
└──────────────────────────────────────────────────────────────────────┘
```

### The one seam that matters

Everything above the engine talks to it through a **single interface**:

```python
StageResult = engine.run_stage(
    stage: str,                 # one of STAGE_ORDER
    config: PipelineConfig,     # existing, already accepts a preloaded config
    manifest: SampleManifest,   # existing contract
    mode: RunMode,
    roots: ProjectRoots,        # NEW: input/config/results/logs/reports/provenance
) -> StageResult                # status, outputs, timing, tool+db versions, error
```

Why this is the right seam, in the vocabulary of `codebase-design`:

- **Depth.** Today the engine is reachable two ways — `run_pipeline()` for
  everything, or 18 hand-written per-stage wrappers. Both are shallow at
  the point of use: the caller must know the prerequisite graph, choose a
  mode, wire roots, and read stdout. `run_stage` puts the prerequisite
  expansion, mode resolution, logging and provenance behind one call.
- **Leverage.** One interface serves the CLI, the HTTP API, and Snakemake.
  Writing a new front end costs nothing; a new stage costs one function.
- **Locality.** Stage selection, prerequisites, resource declaration and
  provenance are decided in one place instead of in 18 wrappers plus a
  Snakefile plus a pilot orchestrator.
- **Deletion test.** Delete the 18 `scripts/<stage>/run_*.py` wrappers and
  the complexity reappears in every caller. They earn nothing today; they
  are 5 near-verbatim copies of `run_stage.py` and 10 near-verbatim copies
  of each other.

### Two adapters means a real seam

`codebase-design` says one adapter means a hypothetical seam. The engine
seam is real: there are already **three** callers of the same stages —
`run.py`, the Snakefile, and `scripts/run_downstream_stages.py`. The app
layer is the fourth, and it is the first to need per-stage results,
resumption and provenance as *data* rather than as files on disk.

### The GUI/CLI symmetry requirement

`codebase-design`'s "return results, don't produce side effects" applies
directly. The rule:

> Every CLI verb is a thin argument parser over one backend function.
> No scientific logic, no scheduling, and no formatting lives in a CLI verb.

`papipeline_app/backend.py` holds the verbs. `cli/` parses argv into a verb
call; `web/` calls the identical function. A test asserts the two surfaces
expose the same verb set, so the GUI can never drift from the CLI.

---

## 3. Layer contracts

| Layer | May import | Must not |
|---|---|---|
| `papipeline_app.web` | `papipeline_app.backend` | import `papipeline.stages` directly |
| `papipeline_app.cli` | `papipeline_app.backend` | import `papipeline.stages` directly |
| `papipeline_app.backend` | `store`, `engine`, `resources`, `tools`, `results`, `provenance` | shell out to tools directly |
| `papipeline_app.engine` | `papipeline.run`, `papipeline.stages`, `papipeline.models` | know about HTTP, SQLite, or the filesystem layout of a project |
| `papipeline_app.store` | stdlib `sqlite3` | know about stages or tools |
| `papipeline/stages/*` | unchanged | import anything from `papipeline_app` |

The last row is the load-bearing one: **the engine never imports the
application.** Dependency arrows point one way only, so the engine remains
usable with no app layer installed — which is what the existing 641 tests
already assume.

---

## 4. What must not change

Non-negotiable, carried forward from `docs/scientific_rules.md`:

1. **No scientific redefinition.** `ClaimStatus` keeps having no `CAUSAL`
   member, so causality stays inexpressible in the schema.
2. **No data-contract change.** The 14 dataclasses in `papipeline/models.py`,
   the TSV column sets, and the sentinel rules in `papipeline/io/tsv.py`
   (`""`, `.`, `na`, `nan`, `none`, `null`, `-` → `None`) are fixed.
3. **R/I/S never becomes an MIC.** `tests/unit/test_phenotype.py::
   TestNoFabricatedQuantitation` guards this.
4. **Missing stays missing.** A stage that cannot run is `FAILED` or
   `SKIPPED`; it is never silently completed.
5. **Database versions are always recorded**, and never updated implicitly
   mid-run (`annotation.allow_database_update: false` and the same for amr
   and virulence).
6. **Selection precedes metadata.** Cohort choice must not be informed by
   phenotype. `tests/unit/test_pilot.py` already asserts this for the pilot
   path; the same invariant is required of the new project-creation flow.

---

## 5. Open decisions (need your answer before implementation)

These are hard to reverse and I have not guessed them.

| # | Question | Why it blocks |
|---|---|---|
| D1 | **REAL mode.** `config/config.yaml:186` has `runtime.allow_real_mode: false`, and `resolve_mode()` raises before reading any data. A GUI that cannot run real data is not a usable application. Should the app layer (a) leave the gate and require an explicit config edit, (b) add a separate `allow_real_mode: app` acknowledgement, or (c) remove the gate? **Recommendation: (b)** — the gate exists to stop an unattended run reading real data, and an interactive user creating a project is exactly the case it should permit, while still requiring a deliberate act. |
| D2 | **`references.tsv` is entirely `UNPINNED`** (9 unpinned, 2 `na`). `docs/reproducibility.md` and scientific rule 8 require a recorded version. Scientific rule 8 says a version is *always recorded*; it does not say what to do when unpinned. **Should the app refuse to start a REAL run with unpinned references, or run and stamp the report?** Recommendation: warn loudly and stamp; refusing makes the app unusable on day one, and the rule's intent (never silently claim a version) is met by stamping. |
| D3 | **`environment/environment.yml` does not match the machine.** It pins `bakta=1.9.3`, `amrfinder=3.12.0`, `pyseer=0.4.5`, `panaroo=4.0.0`, `snp-sites=2.6.3`, `iqtree=2.2.6`. Installed and measured: `bakta 1.12.1`, `amrfinder 4.2.7`, `pyseer 1.1.2`; `panaroo`, `snp-sites` and `iqtree2` **have no installable `osx-arm64` build**. Should the app treat `environment.yml` as authoritative (and fail), or detect what is present and record reality in provenance? **Recommendation: detect and record reality**, and add an `environment/constraints-osx-arm64.yml` documenting what is actually obtainable. A single manifest cannot describe two platforms. |
| D4 | **Web framework.** Not chosen. Options: FastAPI+uvicorn (async, OpenAPI generated free, adds 2 deps) vs Flask (1 dep, sync, simpler) vs stdlib `http.server` (0 deps, not production-viable). **Recommendation: FastAPI** — the OpenAPI schema gives the GUI a generated client and makes "GUI and CLI expose the same verbs" mechanically checkable. |
| D5 | **`papipeline/pilot/`.** Migrate its unique value (cohort integrity assessment, PDC field parsing, the Bakta sweep executor, figure/PPTX output) into the app layer and retire it, or keep it as a supported second front end? **Recommendation: migrate and retire** — it cannot produce stages 6/7/8, so keeping it guarantees two divergent answers to "what does this pipeline do". |
| D6 | **Where does SQLite live?** `results/` and `reports/` are git-ignored (`.gitignore:20,22`) and a SQLite file is currently *not* ignored, so a DB at the repo root would be tracked by accident. **Recommendation: `~/.local/share/papipeline/app.db` by default, overridable by `PAPIPELINE_HOME`.** Per-project state is portable via the project directory, which is the thing that should be archived. |

---

## 6. Milestones

| # | Milestone | Gate |
|---|---|---|
| M1 | These nine documents | agreed |
| M2 | `pyproject.toml`, package installs, `papipeline` importable without `sys.path` hacking | 641 tests still pass |
| M3 | Persisted stage state: the six in-memory stages gain `load_*_from_disk`; `run_stage.py` becomes the single per-stage entry point | 641 tests pass; each stage demonstrably runnable alone |
| M4 | SQLite store + `ProjectRoots` on disk; CLI (`project create/list`, `run`, `status`, `resume`, `tools`, `report`) | new tests green, 641 pass |
| M5 | Snakemake profiles (`macbook`/`workstation`/`hpc`) with real `threads:`/`resources:`; DAG executed for the first time | DAG dry-run and real run both green |
| M6 | Resource-aware Bakta executor: shared DB, checkpointing, per-genome logs, measured-throughput recommendations | resume proven by killing mid-sweep |
| M7 | Results registry + report generation (full / AMR / GWAS / antibiotic-comparison) | registry discovers every declared artefact type |
| M8 | Provenance capture; `references.tsv` real versions filled in | a run is reproducible from its own provenance record |
| M9 | HTTP API | API contract test green |
| M10 | Web UI, one screen at a time, dashboard first | CLI/GUI verb-parity test green |

M2 and M3 come before any UI work deliberately: until a stage can be run
alone and resumed, a progress monitor has nothing truthful to show.
