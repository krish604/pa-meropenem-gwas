#!/usr/bin/env bash
# =============================================================================
# run_meropenem_real.sh -- start the REAL analysis run on the Linux machine.
#
# The single authorised entry point for a real-genome run on this machine. It
# does three things, in this order, and the order is the safety mechanism:
#
#   1. refuses unless CONFIRM_REAL=yes is in the environment. Anything else -
#      unset, "1", "true", "YES" - stops here, with a message saying what is
#      needed. This is AGENTS.md rule 6: no real-genome run until the operator
#      has explicitly said `run real samples`. The confirmation lives in the
#      environment rather than on the command line so it cannot be forgotten by
#      a shell that reuses history, and so it never appears in a tracked file.
#   2. only then exports the session override that opens REAL mode for this one
#      process. config/machines/linux.yaml keeps runtime.allow_real_mode
#      closed - permanently, in a tracked file - and the override is how a run
#      is authorised without editing it. Both gates have to open: the operator's
#      word above, and the pipeline's own check at DAG-build time.
#   3. runs the analysis with the repository's workflow, targeting `full_run`
#      (one process, one pass, every stage exactly once) under
#      machine=config/machines/linux.yaml.
#
# Usage:
#   CONFIRM_REAL=yes scripts/linux/run_meropenem_real.sh
#   CONFIRM_REAL=yes scripts/linux/run_meropenem_real.sh --dry-run
#   CONFIRM_REAL=yes scripts/linux/run_meropenem_real.sh --cores 8
#
#   --cores N   cores for the workflow (default: the overlay's runtime.threads)
#   --dry-run   print the plan and the pinned-state warnings; execute nothing
#   --help
#
# Before the first run of a session, in order (docs/LINUX_RUN.md is the list):
#
#   * the environment exists and `bash environment/versions.sh` has recorded
#     what is in it;
#   * databases are provisioned and pinned - the pipeline never downloads one;
#   * `data/metadata/sample_metadata.tsv` and the phenotype table exist to the
#     contract in docs/data_contract.md; phenotypes are S/I/R, never derived
#     from an MIC, and an absent phenotype is never defaulted to S;
#   * stage 2 output exists, i.e. scripts/linux/run_bakta.sh has been run, if
#     you want it reused (annotation.reuse_tool_output must then say `prefer`
#     or `require`);
#   * `--dry-run` first: it prints the plan without reading a genome.
#
# Nothing environmental is hard-coded: the machine overlay supplies threads and
# paths, the cores count comes from configuration or from --cores, and the
# machine path is a repository-relative one.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
REPO_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"  # scripts/linux -> repo root

MACHINE="config/machines/linux.yaml"
SNAKEFILE="workflow/Snakefile"
TARGET="full_run"
CORES=""
DRY_RUN=0
PYTHON="${PYTHON:-python3}"

usage() {
    sed -n '2,48p' "$0" | sed 's/^# \{0,1\}//'
}

die() { printf 'ERROR: %s\n' "$*" >&2; exit 2; }
log() { printf '[%s] %s\n' "$(date -u '+%H:%M:%S')" "$*" >&2; }

while [ $# -gt 0 ]; do
    case "$1" in
        -h|--help) usage; exit 0 ;;
        --cores)   [ $# -ge 2 ] || die "--cores needs a number"; CORES="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        *) printf 'unknown argument: %s\n\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
done

# --- gate 1: the operator's word --------------------------------------------
if [ "${CONFIRM_REAL:-}" != "yes" ]; then
    cat >&2 <<'EOF'
ERROR: refused: this is a REAL-genome run and it has not been confirmed.

       AGENTS.md rule 6: do not start a real-genome run until the user has
       explicitly said `run real samples`. Nothing has been read, executed or
       written.

       To authorise ONE session, re-run with:

           CONFIRM_REAL=yes scripts/linux/run_meropenem_real.sh

       To inspect the plan without running anything:

           CONFIRM_REAL=yes scripts/linux/run_meropenem_real.sh --dry-run

       The machine overlay (config/machines/linux.yaml) keeps
       runtime.allow_real_mode: false permanently; this script opens it for one
       session only, and never by editing a tracked file.
EOF
    exit 3
fi

# --- gate 2: open REAL mode for this session --------------------------------
export PIPELINE_ALLOW_REAL_MODE=1

cd "$REPO_ROOT"
[ -f "$MACHINE" ]  || die "machine overlay not found: $MACHINE (run from the repository, or check the checkout)"
[ -f "$SNAKEFILE" ] || die "Snakefile not found: $SNAKEFILE"
command -v "$PYTHON" >/dev/null 2>&1 \
    || die "'$PYTHON' not found. Activate the pipeline environment first (docs/LINUX_RUN.md step 2)."

# Threads and paths come from the overlay, read through the pipeline's own
# loader so an invalid overlay fails here with its own message rather than
# halfway through the DAG.
resolve_cores() {
    "$PYTHON" - "$MACHINE" "${CORES:-}" <<'PY'
import sys
from pathlib import Path
from papipeline.config.loader import load_machine_config

machine = Path(sys.argv[1])
override = sys.argv[2]
if override:
    if not override.isdigit() or int(override) < 1:
        raise SystemExit(f"--cores must be a positive integer, got {override!r}")
    print(int(override))
else:
    print(int(load_machine_config(machine).threads or 0))
PY
}

CORES="$(resolve_cores)" || die "could not read runtime.threads from $MACHINE; pass --cores N"
[ "$CORES" -gt 0 ] 2>/dev/null || die "thread count from $MACHINE is not usable; pass --cores N"

command -v snakemake >/dev/null 2>&1 \
    || die "snakemake not found on PATH. Activate the environment created by
scripts/linux/bootstrap_linux.sh (it records what is in it with
environment/versions.sh)."

log "machine : $MACHINE"
log "cores   : $CORES"
log "mode    : REAL (authorised for this session only: CONFIRM_REAL=yes)"
if [ "$DRY_RUN" -eq 1 ]; then
    log "dry run : the plan and any pinned-state warnings, no genome read"
fi

set -- snakemake \
    --snakefile "$SNAKEFILE" \
    --config mode=REAL machine="$MACHINE" \
    --cores "$CORES" \
    --rerun-incomplete \
    --keep-going
[ "$DRY_RUN" -eq 1 ] && set -- "$@" --dry-run
set -- "$@" "$TARGET"

log "running: $*"
exec "$@"
