#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# lib.sh -- shared helpers. Sourced by every script in this directory.
#
# Responsibilities:
#   * load config/run.conf exactly once per invocation
#   * resolve PIPELINE_DIR / Snakefile and fail loudly if missing
#   * expose run_stage(), the single place a pipeline step is invoked
# ---------------------------------------------------------------------------
set -euo pipefail

# This directory IS the runner root: the scripts sit beside config/ and
# results/, so the root is the script's own directory, not its parent.
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config/run.conf
source "$REPO_ROOT/config/run.conf"

RESULTS_ROOT="$REPO_ROOT/$RESULTS_DIR"
SNAKEFILE="$PIPELINE_DIR/workflow/Snakefile"

log()  { printf '\n[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
die()  { printf '\nERROR: %s\n' "$*" >&2; exit 1; }

# The cohort's unique isolate IDs, one per line, in file order.
sample_list() { printf '%s\n' "$REPO_ROOT/config/samples.txt"; }

require_pipeline() {
    [[ -f "$SNAKEFILE" ]] || die "Snakefile not found at $SNAKEFILE
Set PIPELINE_DIR in config/run.conf to your clone of the pipeline."
    command -v snakemake >/dev/null 2>&1 || die "snakemake is not on PATH.
Create the environment from the pipeline's pinned file:

    conda env create -f $PIPELINE_DIR/environment/environment-linux.yml
    conda activate pa-amr"
    require_databases
}

# REAL mode is shut in every committed machine overlay (laptop, bigmachine,
# smoke), by design and by two tests that guard it. The pipeline offers a
# session-scoped escape hatch -- PIPELINE_ALLOW_REAL_MODE -- precisely so that
# a run can be authorised without editing a tracked file and dirtying the
# working tree for days.
#
# This runner opens that hatch ONLY when the caller left the variable unset.
# A value the caller wrote is passed through untouched, including a falsey one:
# overriding an explicit PIPELINE_ALLOW_REAL_MODE=0 with 1 would enable exactly
# what the caller asked to keep disabled. Whatever arrives is validated by the
# pipeline itself, which rejects anything that is neither truthy nor falsey
# rather than guessing (loader.py: _TRUTHY / _FALSEY).
open_real_mode() {
    if [[ -z "${PIPELINE_ALLOW_REAL_MODE+x}" ]]; then
        export PIPELINE_ALLOW_REAL_MODE=1
    fi
}

# Reference databases are never downloaded by the pipeline: configuration sets
# analysis.*.allow_database_update=false everywhere, and each stage runs a
# preflight that refuses a database whose release is not the pinned one in
# config/references.tsv. Failing here, with the provisioning command in hand,
# beats failing four stages in.
require_databases() {
    local missing=()

    # The pipeline reads its own paths from the machine overlay; we only check
    # the two databases whose absence is discovered farthest into the run.
    if [[ -d "$PIPELINE_DB_DIR" ]]; then
        [[ -f "$PIPELINE_DB_DIR/bakta_db/version.json" ]] || missing+=("Bakta DB  (light, 2025-02-24)")
        [[ -f "$PIPELINE_DB_DIR/amrfinderplus_db/version.json" ]] || missing+=("AMRFinderPlus DB (2026-08-07.1)")
    else
        missing+=("reference databases (no $PIPELINE_DB_DIR directory)")
    fi

    (( ${#missing[@]} )) || return 0

    printf '\nREFERENCE DATABASES NOT FOUND\n\n' >&2
    local m; for m in "${missing[@]}"; do printf '    %s\n' "$m" >&2; done
    cat >&2 <<EOF

The pipeline provisions nothing: databases are pinned in
config/references.tsv and a stage preflight refuses any other release.
Provision once, next to the pipeline:

    cd $PIPELINE_DIR
    bakta_db download --type light --output db/bakta_db
    amrfinder --update

EOF
    return 1
}

# run_stage <stage> [stage...]
#
# Invokes Snakemake for the named rules with this repo's result root. The
# result root is passed as an ENVIRONMENT variable rather than `--config`
# because the stage subprocesses recompute their own paths from config: a
# --config value would reach the DAG but not the child process that actually
# writes the file.
#
# Stages run with `--rerun-incomplete` so an interrupted 900-isolate run can be
# resumed by re-issuing the same command without redoing finished isolates.
run_stage() {
    (( $# )) || die "run_stage requires at least one stage name"
    require_pipeline
    open_real_mode

    log "running: $*  (machine=$MACHINE cores=$CORES)"
    PIPELINE_RESULTS_ROOT="$RESULTS_ROOT" \
    PIPELINE_ALLOW_REAL_MODE="$PIPELINE_ALLOW_REAL_MODE" \
    snakemake \
        --snakefile "$SNAKEFILE" \
        --config mode=REAL machine="$MACHINE" \
        --cores "$CORES" \
        --rerun-incomplete \
        --keep-going \
        "$@"
}

# Path to a stage's primary table, for existence assertions.
stage_table() { printf '%s/artifacts/stage_tables/%s\n' "$RESULTS_ROOT" "$1"; }
