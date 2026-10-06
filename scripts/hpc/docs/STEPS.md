# Pipeline steps

What each stage does, what it consumes, and what it writes. Phases are split
at the GWAS gate enforced by `scripts/hpc/run_all.sh`.

## Phase 1 — safe to run unattended

| # | Snakemake rule | Reads | Writes |
|---|---|---|---|
| 1 | `validation` | `assemblies/*.fna` | `stage_tables/01_validation.tsv` |
| 2 | `annotation` | curated Bakta GFF/TSV | `stage_tables/02_annotation.tsv` |
| 3 | `mlst` | assemblies | `stage_tables/03_mlst.tsv` |
| 4 | `amr` | assemblies | `stage_tables/04_amr.tsv` |
| 5 | `virulence` | assemblies | `stage_tables/08_virulence.tsv` |
| 6 | `variants` | assemblies, PAO1 reference | `variants.tsv`, `06_regulators.tsv` |
| 6a | `cohort_variants` | `variants.tsv` | `cohort_variants.tsv` |
| 7 | `pangenome` | assemblies | panaroo graph + `07_gene_families.tsv` |
| 8 | `recombination` | pangenome alignment | gubbins-masked alignment |
| 9 | `phylogeny` | masked alignment | tree (Newick) |
| 10 | `similarity` | tree, alignment | distance matrix |
| 11 | `phenotype` | `config/phenotype.tsv` | `phenotype.tsv` |

Stage 2 **reuses** existing Bakta output and never executes Bakta. If you do
not have Bakta output yet, run `bakta` yourself over `assemblies/` first.

## Gate

`scripts/hpc/run_all.sh` stops here. Nothing beyond this point runs without
`results/.gwas-consent` containing exactly `APPROVED`.

## Phase 2 — requires approval

| # | Snakemake rule | Reads | Writes |
|---|---|---|---|
| 12 | `gwas` | gene families, mutations, phenotype | association tables |
| 13 | `convergence` | association output | convergence table |
| 14 | `cooccurrence` | feature matrix | pairwise association table |
| 15 | `reporting` | all tables | `REPORT.md` |
| 16 | `combine` | all tables | consolidated results |

Stage 14 refuses rather than runs if the feature universe exceeds
`max_features` (2,500 by default). This is intentional: the uncapped cost at
n=900 projects to ~101 million pairs.

## Resume semantics

`results/.phase1-complete` and `results/.phase2-complete` mark finished
phases. Re-running `scripts/hpc/run_all.sh` skips whatever has already finished,
so an interrupted multi-day run resumes rather than restarting.

## Running one step alone

```bash
source scripts/hpc/lib.sh
run_stage amr          # any rule name from the tables above
```

The GWAS gate is still enforced: `run_stage gwas` fails until consent exists.
