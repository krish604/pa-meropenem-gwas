# Reproducibility

A result that cannot be reproduced is an anecdote. This document describes what
the pipeline records, how to reproduce a run, and — importantly — what is
**not yet** pinned.

## Current status

**All 10 references in `config/references.tsv` are unpinned.** Every run reports
this as a warning, and `run_manifest.json` records
`version_status: unpinned` for each. This is expected during the build phase
and is a blocker for any real analysis (scientific rule 8).

Reproducing a TEST-mode run today requires only Python and the pinned packages
in `environment/environment.yml`. Reproducing a REAL-mode run is not yet
possible, because no reference database is pinned and no analysis has been
authorised.

## Environment

```bash
micromamba env create -f environment/environment.yml
micromamba activate pa-amr
# or: conda env create -f environment/environment.yml
```

Package versions are pinned with `=<version>`. The YAML file is the
authoritative record of the software environment; `config/references.tsv` is
the authoritative record of the *databases*. Keep them reconciled — a tool
version without its database version does not describe a reproducible run.

## Recording versions

```bash
bash environment/versions.sh                # human-readable
bash environment/versions.sh --json         # machine-readable
bash environment/versions.sh --out versions_2026-06-01.txt
```

The script only *reads* versions. It never installs, updates or downloads
anything. It reports:

- the status and version of every tool the stages depend on;
- every reference in `config/references.tsv` with its `version_status`;
- a warning listing any unpinned reference.

Run it before and after any analysis and keep both outputs.

## What every run records

`results/<mode>/run_manifest.json`:

| Field | Meaning |
|---|---|
| `pipeline_version` | `papipeline.__version__` |
| `run_mode` | `TEST` or `REAL` |
| `generated_at_utc` | ISO-8601 UTC timestamp |
| `antibiotic` | the antibiotic analysed |
| `python`, `platform` | interpreter and OS |
| `stages` | per-stage status |
| `outputs` | every artefact path |
| `tools_detected` | executable path, detected version, availability |
| `references` | the pinned reference table merged with detected versions |

Where a reference is unpinned, the **detected** tool version is recorded but
`version_status` stays `unpinned`, so the gap is visible rather than papered
over.

Stage tables additionally carry `evidence_source`, `database` and
`database_version` per row (stages 4 and 8).

## Database discipline

1. Databases are **provisioned out of band**. Nothing downloads during a run.
2. All `allow_database_update` flags in `config.yaml` are `false`.
3. `adapters.ToolAdapter.ensure_database()` raises `ToolExecutionError` if a
   different version is requested mid-run.
4. `amr.allow_database_update: false` and `amr.allow_database_update: false`
   guard the AMR stage specifically.
5. An unpinned reference is surfaced as a run warning, never treated as pinned.

To pin a database: provision it, set `tool_version`, `database_version` and
`version_status=pinned` in `config/references.tsv`, and confirm with
`versions.sh`.

## Determinism

| Artefact | Deterministic | Notes |
|---|---|---|
| Synthetic fixtures | yes | single integer seed; byte-identical across runs |
| Figure data | yes | pure functions of stage outputs |
| Master table | yes | fixed column order, sorted multi-values |
| Stage 1 QC | yes | integer and float arithmetic only |
| Stage 13 convergence | yes | dict ordering is insertion-ordered |
| Stage 12 GWAS | yes for the reference engine | Fisher's exact is deterministic; `pyseer` is not |
| Stages using external tools | tool-dependent | depends on the pinned tool and database versions |

Verified by `tests/integration/test_pipeline.py::TestDeterminism`, which
regenerates the fixtures and re-runs the pipeline and compares byte-for-byte.

## Reproducing a TEST run

```bash
# 1. environment
micromamba env create -f environment/environment.yml && micromamba activate pa-amr

# 2. record the starting state
bash environment/versions.sh --out versions_before.txt

# 3. regenerate the fixtures (optional; they are checked in and deterministic)
python3 scripts/common/make_synthetic_data.py

# 4. verify
python3 -m pytest tests -q

# 5. run
python3 scripts/common/run_pipeline.py --mode TEST

# 6. record the ending state and the manifest
bash environment/versions.sh --out versions_after.txt
cat results/test/run_manifest.json
```

The same run via Snakemake:

```bash
snakemake --snakefile workflow/Snakefile --config mode=TEST --cores 4
```

## Mode separation

`TEST` reads `test_data/`, writes `results/test/`.
`REAL` reads `data/`, writes `results/real/`.

The roots are separate directories resolved by `PipelineConfig.data_root()`,
so synthetic and real data cannot be mixed. `REAL` additionally requires
`runtime.allow_real_mode: true`, currently `false`; attempting it raises
`ModeNotAllowedError` before any data is read, including during
`snakemake --dry-run`.

`synthetic.py.generate()` refuses to write into a directory named `data`.

## Connecting the real dataset (future phase)

Not done, and deliberately so. When authorised:

1. Build `data/metadata/sample_metadata.tsv` with one row per assembly.
2. Build `data/phenotype/imipenem_phenotype.tsv` to the contract in
   `docs/data_contract.md`. Do not derive categories from MIC or vice versa.
3. Provide `data/phylogeny/` once alignments and a tree exist.
4. Pin every reference in `config/references.tsv`.
5. Set `runtime.allow_real_mode: true`.
6. Dry-run first: `--mode REAL --dry-run` prints the plan and the pinned-state
   warnings without reading a genome.

The `data/` directory already holds 835 `GCA_*` assemblies downloaded by the
pre-existing `run_pdc.sh`. This build phase has not read, indexed or analysed
any of them. `data/` is not a pipeline input until step 1 above is done —
stage 1 reads assemblies through the manifest, not by scanning the directory.

## Known reproducibility gaps

1. **No database is pinned.** Blocks any real analysis.
2. **No reference genome is pinned.** `ref_reference_genome` is `UNPINNED`.
   Variant calling in stages 6 and 7 requires it.
3. **The Snakemake DAG is unexecuted.** Snakemake is not installed in the
   build environment, so `workflow/Snakefile` has not been run. The native
   orchestrator, which shares every stage function with it, is tested
   (498 tests). Treat the Snakefile as reviewed-but-unverified and run
   `--dry-run` first when Snakemake becomes available.
4. **The GWAS reference engine is not pyseer.** It has no kinship correction.
   It is a mechanics-testing stand-in, not a production GWAS engine. See
   `docs/architecture.md` limitation 2.
5. **No plotting backend.** `viz.py` prepares figure data and writes JSON; the
   drawing layer is not implemented.
