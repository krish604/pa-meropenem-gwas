#!/usr/bin/env bash
# =============================================================================
# bootstrap_linux.sh -- create the Linux analysis environment, record versions.
#
# What this does, in order:
#
#   1. reads the environment name out of environment/environment-linux.yml
#      (it is not written here, so the two cannot drift);
#   2. creates it with micromamba, skipping creation if it already exists;
#   3. runs environment/versions.sh INSIDE that environment, so every tool and
#      database version this machine is about to use is written down before
#      anything is analysed. A run whose provenance cannot be reproduced is not
#      a result, it is an anecdote.
#
# What it deliberately does NOT do: this script never starts an analysis, never
# touches data/, db/ or a phenotype table, and never downloads a database.
# Databases are provisioned out of band (see docs/LINUX_RUN.md) because the
# pipeline refuses a release that is not the one pinned in config/references.tsv.
#
# Usage:
#   scripts/linux/bootstrap_linux.sh                  # create, then record
#   scripts/linux/bootstrap_linux.sh --dry-run        # solve only, install nothing
#   scripts/linux/bootstrap_linux.sh --versions-out environment/versions.linux.txt
#
#   --file PATH       environment file        (default: environment/environment-linux.yml)
#   --name NAME       environment name        (default: read from the file)
#   --yes             non-interactive create
#   --dry-run         micromamba solves and reports, linking nothing
#   --versions-out F  write the version record to F as well as stdout
#
# Nothing environmental is hard-coded: the file, the name and the micromamba
# binary all come from configuration or from the caller's PATH.
#
# KNOWN STATE, MEASURED 2026-10-07: environment/environment-linux.yml does NOT
# solve for linux-64 as pinned (two blockers: gubbins has no linux-64 build for
# Python 3.11.16, and samtools=0.1.19 needs openssl <=1.1.1). See the header of
# that file and .build/linux-ready.solve.linux-64.fresh-full.txt. This script
# will therefore fail at step 2 until that is deliberately re-solved - which is
# the truthful outcome, not a bug in this script.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
REPO_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"  # scripts/linux -> repo root

ENV_FILE="$REPO_ROOT/environment/environment-linux.yml"
ENV_NAME=""
MICROMAMBA="${MICROMAMBA:-micromamba}"
DRY_RUN=0
ASSUME_YES=0
VERSIONS_OUT=""

log() { printf '[%s] %s\n' "$(date -u '+%H:%M:%S')" "$*" >&2; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

usage() {
    # The header block, lines 2-40. No `--help` short form exists for the
    # details: the state of environment-linux.yml is part of what the operator
    # needs before typing anything.
    sed -n '2,40p' "$0" | sed 's/^# \{0,1\}//'
}

while [ $# -gt 0 ]; do
    case "$1" in
        -h|--help) usage; exit 0 ;;
        --file) [ $# -ge 2 ] || die "--file needs a path"; ENV_FILE="$2"; shift 2 ;;
        --name) [ $# -ge 2 ] || die "--name needs a value"; ENV_NAME="$2"; shift 2 ;;
        --versions-out) [ $# -ge 2 ] || die "--versions-out needs a path"; VERSIONS_OUT="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        --yes|-y) ASSUME_YES=1; shift ;;
        *) printf 'unknown argument: %s\n\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
done

# --- 0. prerequisites --------------------------------------------------------
[ -f "$ENV_FILE" ] || die "environment file not found: $ENV_FILE"
command -v "$MICROMAMBA" >/dev/null 2>&1 \
    || die "micromamba not found (looked for '$MICROMAMBA').
Install it, or point MICROMAMBA at the binary:
    MICROMAMBA=/path/to/micromamba scripts/linux/bootstrap_linux.sh"
command -v bash >/dev/null 2>&1 || die "bash not found on PATH"

# The name comes from the file, not from this script. Reading it with PyYAML
# when available and a comment/`name:`-line parse otherwise keeps the two in
# step without requiring the environment we are about to build.
name_from_file() {
    if command -v python3 >/dev/null 2>&1 \
        && python3 -c 'import yaml' >/dev/null 2>&1; then
        python3 - "$1" <<'PY'
import sys, yaml
name = (yaml.safe_load(open(sys.argv[1], encoding="utf-8")) or {}).get("name")
if not name:
    raise SystemExit("no `name:` in " + sys.argv[1])
print(name)
PY
    else
        awk '$1 == "name:" { gsub(/["'\'']/, "", $2); print $2; exit }' "$1"
    fi
}

[ -n "$ENV_NAME" ] || ENV_NAME="$(name_from_file "$ENV_FILE")" \
    || die "could not read the environment name from $ENV_FILE; pass --name"
[ -n "$ENV_NAME" ] || die "the environment name in $ENV_FILE is empty"

log "environment file: $ENV_FILE"
log "environment name: $ENV_NAME"

# --- 1. create (or confirm) the environment ---------------------------------
exists=0
if "$MICROMAMBA" env list 2>/dev/null | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    exists=1
fi

if [ "$exists" -eq 1 ]; then
    log "environment '$ENV_NAME' already exists -- creation skipped"
else
    cmd=("$MICROMAMBA" env create -f "$ENV_FILE" -n "$ENV_NAME")
    [ "$ASSUME_YES" -eq 1 ] && cmd+=(--yes)
    [ "$DRY_RUN" -eq 1 ] && cmd+=(--dry-run)
    log "creating: ${cmd[*]}"
    "${cmd[@]}" || die "micromamba could not build '$ENV_NAME' from $ENV_FILE.
The solve has known blockers (see the header of that file). Reproduce the
failure on its own with:
    micromamba create --dry-run --platform linux-64 -f $ENV_FILE"
fi

if [ "$DRY_RUN" -eq 1 ]; then
    log "--dry-run: nothing was installed; stopping before the version record"
    exit 0
fi

# --- 2. record every tool and database version -------------------------------
# Run inside the environment, not with whatever python3 happens to be first on
# PATH: the point of the record is what THIS machine will actually execute.
VERSIONS="$REPO_ROOT/environment/versions.sh"
[ -f "$VERSIONS" ] || die "versions script not found: $VERSIONS"

log "recording versions via environment/versions.sh"
versions_cmd=("$MICROMAMBA" run -n "$ENV_NAME" bash "$VERSIONS")
if [ -n "$VERSIONS_OUT" ]; then
    versions_cmd+=(--out "$VERSIONS_OUT")
fi
(cd "$REPO_ROOT" && "${versions_cmd[@]}") \
    || die "version record failed; a run without recorded provenance is not a result"

log "bootstrap complete. Next: provision the databases, then see docs/LINUX_RUN.md"
