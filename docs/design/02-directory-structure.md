# Design 2 — Proposed directory structure

Companion to [1 architecture](01-architecture.md). Nothing here exists yet.

## Principle

Three trees, three owners, three lifetimes:

| Tree | Owner | Lifetime | Git |
|---|---|---|---|
| `papipeline/` | scientific engine | permanent | tracked |
| `papipeline_app/` | application layer | permanent | tracked |
| project directory | one analysis | disposable, archivable | **not** tracked |

A project directory is self-contained: copy it and the run is portable.
That is what makes `papipeline project export` meaningful later.

---

## 1. Repository layout

```
papipeline/                          UNCHANGED — the scientific engine
├── stages/                          16 stage modules, do not rewrite
├── config/                          knowledge tables, do not rewrite
├── io/                              strict TSV/FASTA, do not rewrite
├── models.py                        the data contracts, do not rewrite
├── run.py                           native orchestrator, do not rewrite
├── pilot/                           RETIRED by migration (see design 8)
├── adapters/                        tool detection — extended, not replaced
├── knowledge/                       claim-ceiling enforcement — do not touch
├── testing/                         synthetic fixtures — reuse for app tests
├── manifest.py, errors.py, logging_utils.py, viz.py, cli.py

papipeline_app/                      NEW — the application layer
├── __init__.py                      __version__
├── backend.py                       every verb, once. CLI and web both call this.
├── errors.py                        AppError hierarchy, mapped to HTTP codes
├── paths.py                         ProjectRoots — the only place a layout is named
│
├── store/
│   ├── __init__.py
│   ├── schema.sql                   the DDL, versioned
│   ├── db.py                        connect(), migrate(), transaction helper
│   ├── projects.py                  project CRUD
│   ├── samples.py                   sample registration + validation results
│   ├── runs.py                      run lifecycle
│   ├── jobs.py                      per-stage job rows, status transitions
│   ├── results.py                   result registry rows
│   └── provenance.py                provenance records
│
├── engine/
│   ├── __init__.py
│   ├── adapter.py                   run_stage() — THE seam
│   ├── stage_state.py               load_*_from_disk for the 6 in-memory stages
│   └── stage_map.py                 stage name -> module, prerequisites, tool needs
│
├── resources/
│   ├── __init__.py
│   ├── detect.py                    CPU/RAM/disk/tmpdir detection
│   ├── profiles.py                  macbook / workstation / hpc
│   ├── recommend.py                 throughput-based recommendations
│   └── bakta.py                     the resource-aware Bakta executor
│
├── tools/
│   ├── __init__.py
│   ├── registry.py                  installed tools + versions
│   ├── databases.py                 database discovery + versions
│   └── preflight.py                 GO / NO-GO, generalised from scripts/preflight.py
│
├── results/
│   ├── __init__.py
│   ├── registry.py                  typed result discovery
│   └── types.py                     result type enum + schema per type
│
├── provenance/
│   ├── __init__.py
│   ├── capture.py                   versions, checksums, config snapshot
│   └── manifest.py                  run_manifest.json writer
│
├── reports/
│   ├── __init__.py
│   ├── builder.py                   full / AMR / GWAS / antibiotic-comparison
│   └── templates/
│
├── cli/
│   ├── __init__.py
│   ├── main.py                      argv -> backend verb
│   └── commands/                    one module per verb group
│
├── web/
│   ├── __init__.py
│   ├── app.py                       FastAPI app factory
│   ├── routes/                      project, run, status, logs, results, tools
│   ├── schemas.py                   request/response models
│   └── static/                      the UI
│
└── testing/
    ├── __init__.py
    └── factories.py                 project/run fixtures for app tests

workflow/                             Snakemake stays authoritative
├── Snakefile                        REWRITTEN to be incremental (see design 6)
├── rules/
│   ├── annotation.smk               per-genome Bakta, checkpointed
│   ├── stages.smk                   one rule per scientific stage
│   └── report.smk
├── profiles/
│   ├── macbook.yaml
│   ├── workstation.yaml
│   └── hpc.yaml
└── envs/                            per-platform tool constraints

docs/design/                         these nine documents
docs/                                existing docs, unchanged
config/                              unchanged
environment/
├── environment.yml                  UNCHANGED, but see decision D3
├── constraints-osx-arm64.yml        NEW — what is actually installable on arm64 macOS
└── versions.sh                      existing provenance recorder — reuse, don't rewrite

tests/                               existing 641 tests, never deleted
├── unit/                            existing 13 files
├── integration/                     existing 1 file
└── app/                             NEW — tests for papipeline_app
    ├── test_store.py
    ├── test_projects.py
    ├── test_validation.py
    ├── test_tools.py
    ├── test_resources.py
    ├── test_bakta_scheduler.py
    ├── test_resume.py
    ├── test_results_registry.py
    ├── test_provenance.py
    ├── test_backend_parity.py
    └── test_web_api.py
```

---

## 2. Project directory layout

Created by `papipeline project create`. Fixed, documented, and the only
layout the application ever assumes.

```
<project>/
├── project.yaml                    name, organism, antibiotics, created_at, schema_version
├── input/
│   ├── genomes/                    the assembled FASTA files (never modified)
│   ├── metadata/
│   │   └── sample_metadata.tsv     sample_id, assembly_path, lineage, data_class
│   ├── phenotype/
│   │   └── <Antibiotic>_phenotype.tsv
│   └── amr/                        OPTIONAL pre-existing AMR data
│       └── amr_determinants.tsv
├── config/
│   ├── config.yaml                 effective config, resolved
│   ├── antibiotics.tsv
│   ├── mechanisms.tsv
│   ├── regulators.tsv
│   ├── references.tsv
│   └── gene_families.tsv
├── results/
│   ├── intermediate/               the engine's contract — UNCHANGED layout
│   │   ├── annotation/<sample_id>.annotation.tsv
│   │   ├── mlst/mlst_results.tsv
│   │   ├── amr/amr_determinants.tsv
│   │   ├── regulators/regulator_variants.tsv
│   │   ├── structural_variants/structural_variants.tsv
│   │   ├── virulence/virulence_factors.tsv
│   │   ├── pangenome/
│   │   ├── phylogeny/
│   │   ├── gwas/
│   │   └── stages/01_validation.tsv … 15_master_table.tsv
│   └── figures/
├── logs/
│   ├── run-<run_id>.log
│   ├── stages/<stage>.log
│   └── bakta/<sample_id>.log       per-genome, as required
├── reports/
│   ├── full_analysis.md / .html
│   ├── amr.md
│   ├── gwas.md
│   └── antibiotic_comparison.md
├── provenance/
│   ├── run_manifest.json
│   ├── checksums.sha256
│   ├── config_snapshot/
│   └── versions.json
└── .papipeline/
    ├── state.json                  resumable cursor
    └── resources.json              the resource profile actually used
```

### Why `intermediate/` keeps the engine's exact layout

`papipeline/stages/*.py` read from hard-coded relative paths under the
intermediate root. If the project directory invents a new layout, every
stage needs changing — which the brief forbids and which would put the
scientific engine at risk. The project directory therefore *adopts* the
engine's layout rather than the reverse.

### Why genomes live in `input/genomes/`, not in place

Two reasons. First, the pipeline must never write to its inputs, and a
project that references files scattered across a user's disk cannot
guarantee that. Second, a project directory that does not contain its
inputs is not archivable. `papipeline project create` **copies or symlinks**
(caller's choice, recorded in `project.yaml`) and records a SHA-256 of
each file so provenance survives the copy.

---

## 3. What is deliberately absent

| Not added | Why |
|---|---|
| genome assembly | Input is assembled FASTA. The brief forbids making assembly a stage. |
| a second results tree | The engine's `intermediate/` layout is the contract. |
| a plugin system | No second implementation exists to justify a seam. `codebase-design`: one adapter means a hypothetical seam. |
| ORM (SQLAlchemy etc.) | `sqlite3` is stdlib and the schema is 8 tables. An ORM would be a dependency with nothing to abstract. |
| a job queue (Celery etc.) | Snakemake is the scheduler. A second scheduler would be a second source of truth. |
| `results/` inside the repo | `.gitignore:20` already ignores it; a project belongs outside the source tree. |
