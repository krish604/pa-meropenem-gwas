#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_phase1.sh -- stages 1-11: per-isolate processing and cohort assembly.
#
# Everything here is PRE-GWAS. It reads assemblies and writes derived tables
# and alignments; it performs no statistical association, so it is safe to run
# without human approval.
#
#   1  validation      assemblies pass structural checks
#   2  annotation      reuse curated Bakta output; Bakta is never re-run
#   3  mlst            sequence typing (scheme paeruginosa)
#   4  amr             AMRFinderPlus determinants
#   5  virulence       VFDB screen
#   6  variants        minimap2/bcftools vs PAO1, per isolate
#   6a cohort_variants polymorphic sites across the cohort
#   7  pangenome       panaroo gene families
#   8  recombination   gubbins masking
#   9  phylogeny       snp-sites + iqtree
#  10  similarity      pairwise distances
#  11  phenotype       phenotype table joined to the sample list
#
# Usage:  scripts/hpc/run_phase1.sh
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

log "PHASE 1/2: isolates -> cohort tables (stages 1-11)"

run_stage \
    validation \
    annotation \
    mlst \
    amr \
    virulence \
    variants \
    cohort_variants \
    pangenome \
    recombination \
    phylogeny \
    similarity \
    phenotype

log "PHASE 1 COMPLETE"
log "outputs under: $RESULTS_ROOT"
