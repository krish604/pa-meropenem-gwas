#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_gwas.sh -- stages 12-16: the statistical association half of the run.
#
# THIS SCRIPT WILL NOT RUN WITHOUT EXPLICIT APPROVAL.
#
# GWAS is the point of no return for interpretation: it produces p-values,
# effect sizes and ranked candidate determinants that get read as findings.
# Running it silently on a cohort nobody has sanity-checked is how a cohort
# with a swapped case/control column turns into a published false positive.
# So the approval is a real, recorded event rather than a `read` prompt that
# a batch script can quietly answer for you.
#
#   12  gwas           pyseer unitig / gene / mutation association
#   13  convergence    convergent-evolution tests
#   14  cooccurrence   pairwise feature association (capped universe)
#   15  report         final tables
#   16  combine        consolidated results
#
# Approval is a file written by run_all.sh when a human answers yes. It is
# checked here too, so re-invoking this step directly cannot bypass the gate.
#
# Usage:  scripts/hpc/run_gwas.sh
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

CONSENT_PATH="$REPO_ROOT/$CONSENT_FILE"

if [[ ! -f "$CONSENT_PATH" ]]; then
    cat >&2 <<EOF

GWAS HAS NOT BEEN APPROVED -- refusing to run.

No consent file at:
    $CONSENT_PATH

Phases 1-2 have already written every cohort table, so you can inspect them
before deciding:

    less $RESULTS_ROOT/artifacts/stage_tables/04_amr.tsv
    less $RESULTS_ROOT/artifacts/stage_tables/phenotype.tsv

To proceed, approve explicitly with:

    scripts/hpc/run_all.sh --approve-gwas

or, if phase 1 finished in an earlier session, simply run:

    scripts/hpc/run_all.sh

which will prompt you. This script will keep refusing until that file exists.

EOF
    exit 3
fi

# The file must carry the exact approval token, not merely exist: an empty or
# truncated file (a crashed shell, a stray touch) is not consent.
if [[ "$(cat "$CONSENT_PATH")" != "APPROVED" ]]; then
    printf '\nCONSENT FILE MALFORMED -- refusing to run.\n' >&2
    printf 'Expected exactly: APPROVED\n' >&2
    printf 'Found:            %s\n' "$(cat "$CONSENT_PATH")" >&2
    printf 'Delete %s and re-approve.\n' "$CONSENT_FILE" >&2
    exit 4
fi

log "GWAS consent on file ($(cat "$CONSENT_PATH")): $(stat -f '%Sm' -t '%Y-%m-%d %H:%M:%S' "$CONSENT_PATH" 2>/dev/null || stat -c '%y' "$CONSENT_PATH")"
log "PHASE 2/2: cohort association (stages 12-16)"

run_stage \
    gwas \
    convergence \
    cooccurrence \
    reporting \
    combine

log "PHASE 2 COMPLETE"
log "final results under: $RESULTS_ROOT"
