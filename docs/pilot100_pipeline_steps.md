# Pipeline steps and tools

> **Archived.** This is the stage-by-stage record of the *pilot-100* run, kept
> for provenance. It is **not** the current pipeline definition. The canonical
> stage list is the 15-stage imipenem GWAS pipeline specified in
> `.scratch/imipenem-gwas-dashboard/spec.md`. Where the two disagree, the spec
> wins.
>
> Moved here from `results/pilot100/PIPELINE_STEPS.md`, which sat under the
> gitignored `results/` directory. Content is otherwise unchanged.

Environment: `micromamba_pilot/pilot100`, Python 3.11.16, macOS arm64 (Apple Silicon).

```
                      data/  (835 assemblies)
                            |
  v  cohort selection: 65 analysable  [ papipeline.pilot.cohort ]
                            |
  ================================================================
  1  ASSEMBLY QC
     python (stdlib)  -  size / contig / N50 / decode integrity
     -> assembly_QC.tsv, assembly_integrity.tsv
  ================================================================
                            |
  2  ANNOTATION
     bakta 1.12.1  +  BaktaDB v6.0 "light" (2025-02-24)
     -> annotation/<sample>/*.gff3 + *.tsv
     -> parser: papipeline.stages.annotation
                            |
  ================================================================
  3  AMR GENE DETECTION
     AMRFinderPlus 4.2.7  +  DB 2026-08-07.1
     -> amrfinder/<sample>.amrfinder.tsv
                            |
  4  RESISTANCE MECHANISM CALLING
     python  -  config/mechanisms.tsv, config/gene_families.tsv
     -> Imipenem_mechanism_summary.tsv, Meropenem_mechanism_summary.tsv
                            |
  5  MLST
     mlst 2.33.1  +  PubMLST scheme "paeruginosa"
     -> MLST_calls.tsv
                            |
  6  PHENOTYPE
     python  -  PDC AST field, R/I/S kept distinct
     -> Imipenem_phenotype.tsv, Meropenem_phenotype.tsv
                            |
  7  GENE / MECHANISM TABLES + FIGURES
     python (matplotlib, python-pptx)
     -> *_gene_*.tsv, *_mechanism_*.tsv, PPTX, QC report
  ================================================================
                            |
  ---------------------- downstream driver -----------------------
  scripts/run_downstream_stages.py
  calls the existing papipeline.stages entry points in order
                            |
  8  PAN-GENOME
     python  -  papipeline.stages.pangenome (from stage 2 output)
     [ panaroo NOT USED - not installable on arm64 ]
     -> pangenome/{gene_presence_absence,core_genes,accessory_genes}.tsv
                            |
  9  CORE-SNP PHYLOGENY
     minimap2 2.31-r1302   assembly -> reference alignment
     samtools 1.24         sort, index
     bcftools 1.24         mpileup -> call -> consensus
     FastTree 2.2.0        tree from core SNPs
     [ snippy NOT USED - not installable on arm64 ]
     -> phylogeny/{tree.nwk,core_snp_alignment.fasta}
                            |
 10  GWAS
     pyseer 1.1.2          mixed model (LMM) + kinship
     run via scripts/pyseer_compat.py  (library-name shims only)
     -> GWAS_results.tsv, gwas/pyseer_results.tsv
                            |
 11  CONVERGENCE
     python  -  papipeline.stages.convergence
     -> convergence_calls.tsv
  ================================================================
                            |
  CONSOLIDATION
     scripts/consolidate_results.py
     -> PIPELINE_RESULTS.md, master_sample_table.tsv
```

## Tools actually invoked

| Step | Tool | Version | Database |
|---|---|---|---|
| 2 Annotation | Bakta | 1.12.1 | BaktaDB v6.0 light, 2025-02-24 |
| 3 AMR detection | AMRFinderPlus | 4.2.7 | 2026-08-07.1 |
| 5 MLST | mlst | 2.33.1 | PubMLST `paeruginosa` |
| 9 Alignment | minimap2 | 2.31-r1302 | - |
| 9 Sort/index | samtools | 1.24 | - |
| 9 Variant calling | bcftools | 1.24 | - |
| 9 Tree | FastTree | 2.2.0 | - |
| 10 GWAS | pyseer | 1.1.2 | - |
| 1, 4, 6, 7, 8, 11 | Python | 3.11.16 | - |

## Tools NOT used, and why

| Tool | Status | Reason |
|---|---|---|
| panaroo | not installed | No installable `osx-arm64` build: `prokka` pins `perl 5.26.2` exactly, older builds need Intel `mkl`. No PyPI release. Stage 8 builds from annotation instead. |
| snippy | not installed | No installable `osx-arm64` build (`perl-bioperl`, `tabixpp`, `vcflib` unavailable). Stage 9 uses minimap2 + bcftools + FastTree. |
| prokka | not installed | Transitive dependency of panaroo only. Bakta supersedes it. |
| samtools / bcftools 0.1.19 | present but rejected | The `pilot100` env ships 2014 builds that cannot read modern SAM and fail silently. Version-checked at run time; Homebrew 1.24 used instead. |

## Compatibility shims

`scripts/pyseer_compat.py` exists only because pyseer 1.1.2 (2021) predates this
environment's libraries. It binds three names to the identical functions they
used to refer to, and computes the kinship matrix pyseer requires as input:

- `smf.Logit(p, v)` / `smf.OLS(p, v)` - pyseer calls the *formula* API with raw
  arrays; the formula wrappers parse argument 1 as a formula string. Bound to
  `statsmodels.api.Logit` / `.OLS`, which take arrays.
- `scipy.arange` and friends - numpy names removed from the scipy namespace.
- Kinship - pyseer 1.1.2 takes the structure matrix as input rather than
  deriving it. Supplied as a normalised Gram matrix (PSD by construction).

No fitted value, test statistic or p-value is altered by these.

---

## Before a large cohort: `scripts/preflight.py`

Run this first. It fails in seconds on anything that would otherwise surface
hours into a sweep, and prints a single GO / NO-GO.

```
micromamba run -n pilot100 python scripts/preflight.py --expect 200
```

It checks: tool presence and version floors, the Bakta and AMRFinderPlus
databases, `config/config.yaml`, the manifest and the genome files it points
at, how many samples are already annotated, and free disk against the
measured per-genome cost.

## Disk: prune unused Bakta output

Bakta writes ~84 MB per genome. The pipeline reads 8.6 MB of that: the
annotation GFF, the main TSV, and (for provenance) the inference TSV. The
rest — `json` 20 MB, `embl` 15 MB, `gbff` 14 MB, `svg` 8 MB, `fna` 6 MB,
`ffn` 6 MB, `faa` 2 MB, `png` 2 MB — is unused.

```
# prune everything already annotated
python scripts/run_annotation_sweep.py --prune

# pruning is automatic after each genome; disable with --no-prune
```

Measured on the 65-genome cohort: **5.2 GB -> 695 MB, 4.85 GB reclaimed.**
For 200 genomes that is the difference between ~17.6 GB and ~1.9 GB.
Pruning keeps the GFF, the annotation TSV and the inference TSV, so it never
invalidates a sample and never triggers a re-annotation.

## Failure isolation

`run_one` is documented "never raises" and the whole body is now guarded.
A missing genome, an empty genome, an unwritable output directory, a missing
binary, a timeout or an unexpected exception all become a state record and
the sweep continues. Previously `bakta_executable()` was called outside the
`try`, so a missing binary would have abandoned every genome not yet
started.
