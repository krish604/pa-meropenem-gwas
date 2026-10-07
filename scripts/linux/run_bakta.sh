#!/usr/bin/env bash
# =============================================================================
# run_bakta.sh -- annotate a manifest's assemblies with Bakta (db-light), in batch.
#
# Stage 2 is the one stage whose output is also an INPUT: the annotation stage
# reuses what is on disk rather than re-running Bakta, so this script is how a
# fresh cohort gets annotated before the pipeline is pointed at it. It annotates
# and nothing else - it never starts the pipeline.
#
# Layout, exactly as README.md and docs/INSTALL.md describe it:
#
#     results/<mode>/intermediate/bakta/<sample_id>/<assembly stem>.gff3
#                                                            .tsv
#                                                            .inference.tsv
#
# i.e. `results/real/intermediate/bakta/<sample_id>/`. That root is NOT written
# here: it comes from config/machines/linux.yaml's `paths.results_root` plus the
# mode, read through the pipeline's own configuration loader, so this script and
# the stage that later reads the same directory cannot disagree about where it
# is. `--out-root` overrides it when you need to annotate somewhere else.
#
# Sample IDs and assemblies come from the manifest, resolved by the pipeline's
# own code - papipeline.manifest for the cohort, papipeline.assemblies for the
# file - never by scanning a directory. A sample whose assembly cannot be
# located is reported by name and skipped; it is never silently dropped, and
# one missing genome does not end the batch.
#
# Reuse: a sample whose existing output passes the stage's own acceptance check
# (papipeline.stages.annotation.decide_reuse) is skipped. That check verifies
# the database named in the GFF header, so "already annotated" means verified,
# not merely present.
#
# Usage:
#   scripts/linux/run_bakta.sh --manifest data/metadata/sample_metadata.tsv
#   scripts/linux/run_bakta.sh --manifest data/metadata --dry-run   # print plan
#   scripts/linux/run_bakta.sh --manifest data/metadata --slurm \
#       --slurm-time 12:00:00 --slurm-mem-gb 32
#
#   --manifest PATH      REQUIRED. `sample_metadata.tsv`, or the directory
#                        containing it.
#   --machine PATH       machine overlay       (default: config/machines/linux.yaml)
#   --science PATH       science config        (default: config/science.yaml)
#   --db PATH            Bakta DB directory    (default: annotation.bakta_db,
#                        i.e. the db-light database, with the pinned-release
#                        preflight the stage itself runs)
#   --out-root PATH      where per-sample directories go
#   --threads N          threads per Bakta process (default: the overlay's
#                        runtime.threads)
#   --dry-run            resolve everything, print the commands, run nothing
#   --slurm              submit one batch job per pending sample with sbatch
#   --slurm-time T       #SBATCH --time for those jobs (omit for the cluster default)
#   --slurm-mem-gb N     #SBATCH --mem for those jobs (omit for the cluster default)
#   --slurm-log-dir PATH where job stdout/stderr go
#   --help
#
# SLURM: the directives used are the ones scripts/hpc/submit_slurm.sh already
# uses in this repository. `sbatch` was NOT available on the machine this
# script was written on, so the --slurm path is unverified here - it fails with
# a clear message if sbatch is not on PATH rather than pretending to submit.
#
# Nothing environmental is hard-coded: paths come from configuration, threads
# from the overlay, resources from arguments.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
REPO_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"  # scripts/linux -> repo root

MANIFEST=""
SCIENCE="$REPO_ROOT/config/science.yaml"
MACHINE="$REPO_ROOT/config/machines/linux.yaml"
DB=""
OUT_ROOT=""
THREADS=""
DRY_RUN=0
SLURM=0
SLURM_TIME=""
SLURM_MEM_GB=""
SLURM_LOG_DIR=""
PYTHON="${PYTHON:-python3}"

log()  { printf '[%s] %s\n' "$(date -u '+%H:%M:%S')" "$*" >&2; }
die()  { printf 'ERROR: %s\n' "$*" >&2; exit 2; }

usage() {
    sed -n '2,63p' "$0" | sed 's/^# \{0,1\}//'
}

while [ $# -gt 0 ]; do
    case "$1" in
        -h|--help) usage; exit 0 ;;
        --manifest)      [ $# -ge 2 ] || die "--manifest needs a path"; MANIFEST="$2"; shift 2 ;;
        --science)       [ $# -ge 2 ] || die "--science needs a path"; SCIENCE="$2"; shift 2 ;;
        --machine)       [ $# -ge 2 ] || die "--machine needs a path"; MACHINE="$2"; shift 2 ;;
        --db)            [ $# -ge 2 ] || die "--db needs a path"; DB="$2"; shift 2 ;;
        --out-root)      [ $# -ge 2 ] || die "--out-root needs a path"; OUT_ROOT="$2"; shift 2 ;;
        --threads)       [ $# -ge 2 ] || die "--threads needs a number"; THREADS="$2"; shift 2 ;;
        --slurm-time)    [ $# -ge 2 ] || die "--slurm-time needs a value"; SLURM_TIME="$2"; shift 2 ;;
        --slurm-mem-gb)  [ $# -ge 2 ] || die "--slurm-mem-gb needs a number"; SLURM_MEM_GB="$2"; shift 2 ;;
        --slurm-log-dir) [ $# -ge 2 ] || die "--slurm-log-dir needs a path"; SLURM_LOG_DIR="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        --slurm)   SLURM=1; shift ;;
        *) printf 'unknown argument: %s\n\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
done

# --- the contract: a manifest, named ----------------------------------------
if [ -z "$MANIFEST" ]; then
    printf 'ERROR: --manifest is required.\n' >&2
    printf '       Point it at sample_metadata.tsv (or the directory holding it),\n' >&2
    printf '       e.g. --manifest data/metadata/sample_metadata.tsv\n' >&2
    printf '       The contract is docs/data_contract.md; only `sample_id` is required.\n' >&2
    exit 2
fi
[ -e "$MANIFEST" ] || die "manifest not found: $MANIFEST"
[ -f "$SCIENCE" ]  || die "science config not found: $SCIENCE"
[ -f "$MACHINE" ]  || die "machine overlay not found: $MACHINE"
command -v "$PYTHON" >/dev/null 2>&1 \
    || die "'$PYTHON' not found. Activate the pipeline environment first (see docs/LINUX_RUN.md)."

# One helper, two modes: `layout` reports where to write and with what, `plan`
# emits one tab-separated row per sample (id, action, output dir,
# command-or-reason). Both go through papipeline for the cohort, the assembly
# lookup, the file stem, the reuse verdict and the command line, so nothing
# here can drift from what stage 2 will later look for.
#
# Overrides travel as environment variables rather than as more command-line
# flags, because `BAKTA_DB` is already this pipeline's documented way of
# pointing at a database (papipeline.stages.annotation._bakta_database reads it
# first), and BAKTA_OUT_ROOT / BAKTA_THREADS are this script's own.
PLAN_HELPER="$SCRIPT_DIR/bakta_plan.py"
[ -f "$PLAN_HELPER" ] || die "helper not found: $PLAN_HELPER"
[ -n "$DB" ] && export BAKTA_DB="$DB"
export BAKTA_OUT_ROOT="$OUT_ROOT" BAKTA_THREADS="$THREADS"

cd "$REPO_ROOT"
LAYOUT="$( "$PYTHON" "$PLAN_HELPER" layout --science "$SCIENCE" --machine "$MACHINE" )" \
    || die "could not resolve layout (Bakta database, output root or threads) from configuration"

OUT_ROOT="$(printf '%s\n' "$LAYOUT" | sed -n '1p')"
DB="$(printf '%s\n' "$LAYOUT" | sed -n '2p')"
THREADS="$(printf '%s\n' "$LAYOUT" | sed -n '3p')"
[ -n "$OUT_ROOT" ] && [ -n "$DB" ] && [ -n "$THREADS" ] \
    || die "internal: could not parse the layout output"

log "out root : $OUT_ROOT"
log "database : $DB"
log "threads  : $THREADS"

if [ "$SLURM" -eq 1 ]; then
    command -v sbatch >/dev/null 2>&1 \
        || die "sbatch not found on PATH, so nothing can be submitted.
Run without --slurm to annotate in this shell instead.
(The --slurm path was not executed on the machine this script was written on;
it is documented as unverified in docs/LINUX_RUN.md.)"
    [ -n "$SLURM_LOG_DIR" ] || SLURM_LOG_DIR="$OUT_ROOT/slurm-logs"
fi

PLAN="$(mktemp)"
trap 'rm -f "$PLAN"' EXIT

"$PYTHON" "$PLAN_HELPER" plan --science "$SCIENCE" --machine "$MACHINE" \
    --manifest "$MANIFEST" > "$PLAN" \
    || die "could not build the annotation plan from $MANIFEST"

run=0 reuse=0 missing=0 failed=0

while IFS=$'\t' read -r sample_id action outdir payload; do
    [ -n "${sample_id:-}" ] || continue
    case "$action" in
        reuse)
            reuse=$((reuse + 1))
            log "reuse   $sample_id"
            ;;
        missing)
            missing=$((missing + 1))
            printf 'MISSING %s: %s\n' "$sample_id" "$payload" >&2
            ;;
        mismatch)
            failed=$((failed + 1))
            printf 'MISMATCH %s: %s\n' "$sample_id" "$payload" >&2
            ;;
        run)
            run=$((run + 1))
            if [ "$DRY_RUN" -eq 1 ]; then
                printf '%s\n' "$payload"
                continue
            fi
            if [ "$SLURM" -eq 1 ]; then
                job="$OUT_ROOT/slurm/$sample_id.sh"
                mkdir -p "$OUT_ROOT/slurm" "$SLURM_LOG_DIR"
                {
                    printf '#!/bin/bash\n'
                    printf '#SBATCH --job-name=bakta_%s\n' "$sample_id"
                    printf '#SBATCH --output=%s/%s.out\n' "$SLURM_LOG_DIR" "$sample_id"
                    printf '#SBATCH --error=%s/%s.err\n' "$SLURM_LOG_DIR" "$sample_id"
                    printf '#SBATCH --cpus-per-task=%s\n' "$THREADS"
                    [ -n "$SLURM_MEM_GB" ] && printf '#SBATCH --mem=%sG\n' "$SLURM_MEM_GB"
                    [ -n "$SLURM_TIME" ] && printf '#SBATCH --time=%s\n' "$SLURM_TIME"
                    printf 'set -euo pipefail\n'
                    printf 'cd "%s"\n' "$REPO_ROOT"
                    printf '%s\n' "$payload"
                } > "$job"
                if sbatch "$job" >&2; then
                    log "submit  $sample_id ($job)"
                else
                    failed=$((failed + 1))
                    printf 'FAILED   %s: sbatch rejected %s\n' "$sample_id" "$job" >&2
                fi
                continue
            fi
            mkdir -p "$OUT_ROOT"
            log "bakta   $sample_id"
            # Shell-quoted by the planner, so a path with a space is one word.
            if bash -c "$payload"; then
                printf 'ok       %s\n' "$sample_id" >&2
            else
                failed=$((failed + 1))
                printf 'FAILED   %s (exit %s)\n' "$sample_id" "$?" >&2
            fi
            ;;
        *)
            die "internal: unknown plan action '$action' for sample '$sample_id'"
            ;;
    esac
done < "$PLAN"

printf '\nannotation plan: %s to run, %s reused, %s missing assembly, %s refused or failed\n' \
    "$run" "$reuse" "$missing" "$failed" >&2

# The verdict does not depend on --dry-run: a plan that cannot be executed
# cleanly should not read as success just because nothing was executed.
status=0
[ "$failed" -eq 0 ] || status=1
[ "$missing" -eq 0 ] || status=1

if [ "$DRY_RUN" -eq 1 ]; then
    log "--dry-run: nothing was executed"
fi
exit "$status"
