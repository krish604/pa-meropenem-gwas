# Architecture

## Purpose

An assembly-based *Pseudomonas aeruginosa* antimicrobial resistance pipeline
focused on **imipenem**, built so that a second antibiotic can be added through
configuration rather than code.

The pipeline answers a layered question:

```
imipenem phenotype
    ↓
AMR genes  →  acquired determinants
    ↓
chromosomal mutations  →  oprD, regulator loci
    ↓
resistance mechanisms  →  permeability, efflux, AmpC, acquired
    ↓
structural variants / MGE
    ↓
phylogenetic background
    ↓
GWAS associations
    ↓
convergent resistance
    ↓
co-occurring mechanisms
```

It never asserts causation. See `docs/scientific_rules.md`.

## Layering

The design rule is that **science lives in the library and orchestration lives
outside it**. Concretely:

```
papipeline/                 the importable library: all scientific logic
├── models.py               typed records + controlled vocabularies
├── errors.py               typed exception hierarchy
├── logging_utils.py        one logging configuration point
├── io/                     strict TSV and FASTA readers/writers
│   ├── tsv.py              the data contract, enforced in one place
│   └── fasta.py            streaming FASTA + assembly statistics
├── config/                 YAML + knowledge-table loading
│   └── loader.py           PipelineConfig, frozen and validated
├── knowledge/              knowledge-table queries
├── manifest.py             the master sample manifest (the join key)
├── adapters/               external tool detection + adapter base
├── testing/                the synthetic fixture generator
├── stages/                 the sixteen stages, one module each
├── viz.py                  figure data preparation (no plotting)
├── cli.py                  shared CLI output helpers
└── run.py                  the native orchestrator

workflow/Snakefile          thin Snakemake layer -> calls stages/
scripts/<stage>/            thin CLI wrappers  -> calls stages/
```

Three consequences of that rule:

1. **The stage functions are the single implementation.** The Snakefile and
   the CLI wrappers both call `papipeline.stages.*`. A rule that re-implemented
   a parser would be a second parser, and the two would drift.
2. **Parsers are pure functions over text.** They are tested against literal
   fixture strings, so the suite runs with no bioinformatics tool installed.
3. **A run needs only the pinned Python dependencies.** Snakemake, Bakta,
   AMRFinderPlus, Panaroo, IQ-TREE and pyseer are optional and detected at
   runtime.

## Execution paths

| Path | Command | Use when |
|---|---|---|
| Native orchestrator | `python3 scripts/common/run_pipeline.py --mode TEST` | Everyday use. Runs all sixteen stages once. |
| Snakemake | `snakemake --snakefile workflow/Snakefile --config mode=TEST` | Cluster submission, checkpointing, per-stage reruns. |
| Single stage | `python3 scripts/common/run_stage.py --mode TEST --stage amr` | Debugging one stage. |
| Dry run | `python3 scripts/common/run_pipeline.py --mode TEST --dry-run` | Check config and see the plan without running. |

Both paths execute the same functions in the same order.

## Stage ordering

Stages are **numbered** as in the specification but **executed** in dependency
order, because stage 5 consumes the outputs of stages 6 and 7.

| # | Stage | Module | Reads | Key output |
|---|---|---|---|---|
| 1 | Genome validation | `stages/validation.py` | FASTA | `01_validation.tsv` |
| 2 | Annotation | `stages/annotation.py` | Bakta output | `02_annotation_summary.tsv` |
| 3 | MLST | `stages/mlst.py` | mlst output | `03_mlst.tsv` |
| 4 | AMR detection | `stages/amr.py` | AMRFinderPlus | `04_amr.tsv` |
| 6 | Regulator screen | `stages/regulators.py` | variant calls | `06_regulators.tsv` |
| 7 | SV / MGE | `stages/sv.py` | SV calls | `07_structural_variants.tsv` |
| 5 | Mechanism interpretation | `stages/mechanisms.py` | 4, 6, 7, 2 | `05_mechanisms.tsv` |
| 8 | Virulence | `stages/virulence.py` | VFDB search | `08_virulence.tsv` |
| 9 | Pan-genome | `stages/pangenome.py` | annotations | `gene_presence_absence.tsv` + 3 |
| 10 | Phylogenomics | `stages/phylogeny.py` | alignments, tree | `10_phylogeny.tsv` |
| 11 | Phenotype | `stages/phenotype.py` | phenotype TSV | `11_phenotype.tsv` |
| 12 | GWAS | `stages/gwas.py` | 11 + features | `12_gwas.tsv` |
| 13 | Convergence | `stages/convergence.py` | 4, 6, 10 | `13_convergence.tsv` |
| 14 | Co-occurrence | `stages/cooccurrence.py` | 4, 5, 6 | `14_cooccurrence.tsv` |
| 15 | Integration | `stages/integration.py` | all of the above | `15_master_table.tsv` |
| 16 | Figures + report | `viz.py`, `stages/reporting.py` | 15 and others | `test_pipeline_report.md` |

`papipeline.run.EXECUTION_ORDER` is the authoritative list; the Snakefile
mirrors it. `docs/architecture.md` is documentation, not configuration.

## The claim vocabulary

One enum, `models.ClaimStatus`, governs every claim the pipeline makes:

| Status | Meaning | Reached by |
|---|---|---|
| `DETECTED` | directly observed in this sample | stages 4, 6, 7 |
| `PREDICTED` | inferred without direct observation | candidate SVs (stage 7) |
| `ASSOCIATED` | statistically associated in this cohort | stage 12, or stage 13 |
| `SUPPORTED` | associated **and** robust across lineages | stages 12 **and** 13 |
| `UNKNOWN` | nothing to report | default |

There is deliberately **no `CAUSAL` member**. `config/mechanisms.tsv` gives
every gene a `claim_ceiling` of `DETECTED`, and
`knowledge.ceiling_for()` clamps any attempt to report a raw detection as
something stronger. Only stages 12 and 13, which have phenotype and phylogeny
in hand, may raise a claim above `DETECTED`, and `SUPPORTED` requires both.

## Adding an antibiotic

No code change is required. Three configuration edits:

1. `config/config.yaml` — add the name to `antibiotics:`.
2. `config/antibiotics.tsv` — add a row with its class and phenotype standard.
3. `config/mechanisms.tsv` — extend the `antibiotic` column of the loci that
   apply. The column accepts a comma-separated list or the keyword `all`.

The phenotype filename follows `<antibiotic>_phenotype.tsv`. Tests
`test_enabling_meropenem_requires_no_code_change` and
`test_all_keyword_covers_any_antibiotic` cover this path.

## Mode separation

`RunMode.TEST` reads `test_data/`, `RunMode.REAL` reads `data/`. They are
different directories by construction, resolved by
`PipelineConfig.data_root()`. Results go to `results/test/` and `results/real/`.

`REAL` mode additionally requires `runtime.allow_real_mode: true` in
`config.yaml`, which is `false`. Attempting it raises `ModeNotAllowedError`
before any data is read — including during `snakemake --dry-run`, because the
Snakefile re-checks at DAG-build time.

`synthetic.py.generate()` refuses to write into a directory named `data`.

## Extension points

| To add | Do this | Not this |
|---|---|---|
| An antibiotic | edit 3 config tables | add a stage |
| A resistance locus | add a row to `mechanisms.tsv` + `regulators.tsv` | edit a stage |
| An AMR database | subclass `stages/amr.AmrAdapter` | edit the AMR stage |
| A QC tool | fill the hook in `stages/validation.py` | change thresholds |
| A figure | add a builder in `viz.py` | add plotting to a stage |

## Known limitations

These are real and deliberate, not oversights:

1. **Per-stage Snakemake rules recompute prerequisites.** The stage functions
   exchange in-memory structures, so a process running one stage cannot
   reconstruct its predecessors' state. Each rule calls `--only <stage>`, and
   the orchestrator expands that to include prerequisites. Results are correct;
   a full Snakemake run does redundant work. The native orchestrator is the
   efficient path. Per-stage state caching is the fix and is future work — it
   was not built here because Snakemake is not installed, so the DAG could not
   be exercised and an untested caching layer would be a liability.
2. **The GWAS reference engine has no kinship correction.** `gwas.ReferenceEngine`
   (Fisher's exact + Benjamini-Hochberg) exists so the pipeline mechanics are
   testable without pyseer. Its p-values are anti-conservative in a structured
   cohort. Every result it emits is stamped
   `reference_fisher:NO_KINSHIP_CORRECTION`. Use `PyseerEngine` for real work.
3. **Promoter-region assessment is disabled.** It needs a pinned reference
   genome and an aligner; neither is pinned. Promoter calls are dropped with a
   warning rather than guessed.
4. **Promoter calls and SV calls are ingested, not called.** Stages 6 and 7
   consume variant calls from an upstream tool. No variant caller is wired in;
   `config/references.tsv` has an `UNPINNED` `PAO1_reference` row to record
   that gap.
5. **All 10 references are unpinned.** This is reported as a run warning. It
   must be resolved before any real analysis (scientific rule 8).
6. **No plotting is implemented.** `viz.py` prepares the data for all thirteen
   figures and writes it to JSON; the drawing layer is future work. This is
   deliberate — it keeps the data inspectable without a plotting backend and
   keeps biology out of the plotting code.
