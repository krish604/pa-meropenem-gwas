#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_all.sh -- the one command you need.
#
#   scripts/hpc/run_all.sh
#
# What it does, in order:
#
#   1. validates the cohort (no isolate listed twice, right size, all
#      assemblies present) and aborts before compute if anything is wrong
#   2. runs PHASE 1 -- stages 1-11 -- unless those outputs already exist
#   3. STOPS AND ASKS YOUR PERMISSION before GWAS
#   4. only if you approve, runs PHASE 2 -- stages 12-16
#   5. leaves every output under results/
#
# The stop is deliberate and cannot be skipped: answering "no" (or just
# pressing Enter) ends the run with all phase-1 results intact, and you can
# come back days later, look at the tables, and resume with the same command.
# A previously granted approval is remembered in results/.gwas-consent.
#
# Options:
#   --approve-gwas   record approval without prompting, then run phase 2
#                    (for resuming on a machine with no terminal attached)
#   --phase1-only    stop after stage 11 regardless of prior approval
#   --restart        delete consent and re-run from stage 1 (destructive)
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

PHASE1_MARKER="$RESULTS_ROOT/.phase1-complete"
PHASE2_MARKER="$RESULTS_ROOT/.phase2-complete"
FORCE_GWAS=0
PHASE1_ONLY=0
RESTART=0

for arg in "$@"; do
    case "$arg" in
        --approve-gwas) FORCE_GWAS=1 ;;
        --phase1-only)  PHASE1_ONLY=1 ;;
        --restart)      RESTART=1 ;;
        -h|--help)      sed -n '2,24p' "$0"; exit 0 ;;
        *) die "unknown option: $arg (try --help)" ;;
    esac
done

if (( RESTART )); then
    printf '\n--restart removes consent and completion markers, and deletes\n'
    printf 'everything under %s/.\n' "$RESULTS_DIR"
    read -r -p "Type 'yes' to confirm full wipe: " CONFIRM
    [[ "$CONFIRM" == "yes" ]] || { echo "aborted."; exit 1; }
    find "$RESULTS_ROOT" -mindepth 1 -maxdepth 1 ! -name '.gitkeep' -exec rm -rf {} +
    log "results cleared"
fi

mkdir -p "$RESULTS_ROOT"

# --- 0. preflight ------------------------------------------------------------
# Fail here, not four stages in: pipeline clone present, snakemake importable,
# reference databases provisioned at their pinned releases.
log "STEP 0  preflight (pipeline, environment, databases)"
require_pipeline
log "pipeline OK   $PIPELINE_DIR"
log "databases OK  $PIPELINE_DB_DIR"

# --- 1. inputs ---------------------------------------------------------------
log "STEP 1  validating cohort"
"$REPO_ROOT/scripts/hpc/check_inputs.sh"

# --- 2. phase 1 --------------------------------------------------------------
if [[ -f "$PHASE1_MARKER" ]]; then
    log "STEP 2  phase 1 already complete -- skipping to the GWAS gate"
else
    log "STEP 2  phase 1 (stages 1-11). This is the long part."
    "$REPO_ROOT/scripts/hpc/run_phase1.sh"
    date > "$PHASE1_MARKER"
fi

# --- 3. THE GATE -------------------------------------------------------------
if (( PHASE1_ONLY )); then
    log "STEP 3  --phase1-only: stopping here. Phase 1 results are in $RESULTS_DIR/"
    exit 0
fi

if [[ -f "$REPO_ROOT/$CONSENT_FILE" ]] && \
   [[ "$(cat "$REPO_ROOT/$CONSENT_FILE")" == "APPROVED" ]]; then
    log "STEP 3  GWAS already approved -- proceeding to phase 2"
elif (( FORCE_GWAS )); then
    printf 'APPROVED\n' > "$REPO_ROOT/$CONSENT_FILE"
    log "STEP 3  approval recorded via --approve-gwas"
else
    cat <<EOF

============================================================================
 PHASE 1 COMPLETE.  GWAS HAS NOT BEEN APPROVED.
============================================================================

Stages 1-11 are finished. Every cohort table is on disk under:

    $RESULTS_ROOT/artifacts/stage_tables/

The next step runs the genome-wide association (pyseer), convergence tests
and pairwise co-occurrence. Those produce p-values and ranked candidate
determinants -- the part of the output people quote as a finding -- so it
waits for you.

You can inspect the phase-1 results now, in another terminal, before deciding.

  Cohort size ......... $(sample_list | grep -vc '^$' 2>/dev/null || echo '?') isolates
  AMR determinants .... less $RESULTS_ROOT/artifacts/stage_tables/04_amr.tsv
  Phenotype table ..... less $RESULTS_ROOT/artifacts/stage_tables/phenotype.tsv
  Variants ............ less $RESULTS_ROOT/artifacts/stage_tables/variants.tsv

  RESULT HEAD ........ less $RESULTS_ROOT/artifacts/stage_tables/06_regulators.tsv

 NOTE: this cohort is sized for discovery, not confirmation. Interpret
 associations as candidates and correct for multiple testing before quoting
 any single locus.

============================================================================
EOF

    read -r -p "Proceed to GWAS? [y/N] " ANSWER
    case "${ANSWER:-}" in
        y|Y|yes|YES)
            printf 'APPROVED\n' > "$REPO_ROOT/$CONSENT_FILE"
            log "approval recorded at $CONSENT_FILE"
            ;;
        *)
            log "GWAS not approved -- stopping. Phase 1 results kept in $RESULTS_DIR/"
            log "Resume any time with the same command: scripts/hpc/run_all.sh"
            exit 0
            ;;
    esac
fi

# --- 4. phase 2 --------------------------------------------------------------
if [[ -f "$PHASE2_MARKER" ]]; then
    log "STEP 4  phase 2 already complete"
else
    log "STEP 4  phase 2 (stages 12-16)"
    "$REPO_ROOT/scripts/hpc/run_gwas.sh"
    date > "$PHASE2_MARKER"
fi

cat <<EOF

============================================================================
 RUN COMPLETE
============================================================================
All outputs are under:  $RESULTS_ROOT/

    stage tables ...... artifacts/stage_tables/
    GWAS results ...... artifacts/   (see REPORT.md if written)
    logs .............. $RESULTS_DIR/logs/

Re-running the same command is safe: completed phases are skipped.

============================================================================
EOF
