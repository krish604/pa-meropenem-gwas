#!/usr/bin/env bash
# Start the read-only results dashboard.
#
#   scripts/run_dashboard.sh [--results <path>] [--port 8765] [--allow-launch]
#
# The ONE file outside dashboard/ that BACKEND owns. Everything it does is
# delegation: the interpreter, the environment and the arguments go straight
# to `python -m dashboard`, which owns the bind check, the launcher gates and
# the auto-detect probe order. Nothing is decided here that is not decided
# there, because a shell wrapper that could disagree with the application is a
# second answer to the same question.
#
# Why a wrapper at all: `python` is not on PATH in this project. The pa-amr
# micromamba environment is, and this makes the invocation one line for a human
# instead of a sentence they have to remember.
#
# ---------------------------------------------------------------------------
# What this script deliberately does NOT do
# ---------------------------------------------------------------------------
# It never sets PIPELINE_ALLOW_REAL_MODE. That variable belongs on a REAL run's
# child process and nowhere else; setting it here would open the gate for the
# dashboard itself, and UI-D1's sibling rule says the dashboard must not be the
# thing that starts a real-genome run.
#
# It never runs snakemake, bakta, panaroo, gubbins, iqtree or pyseer. The
# launcher's default action is `snakemake --dry-run full_run`, and even that is
# only *displayed* by the API - this script does not execute it. See
# dashboard/DESIGN.md §13: a launcher that inferred a real run's command from
# the pipeline's CLI would be a second, unowned copy of that decision.
#
# `PA_DASH_STATE_DIR` is honoured if set, and deliberately defaults to
# ~/.pa_dashboard, which is outside every results root. The state directory is
# the only place the dashboard writes (UI-D1); `python -m dashboard` refuses a
# --state-dir that falls inside the opened results root, and this wrapper does
# not second-guess it.

set -euo pipefail

# -- repository root: this script's own grandparent directory ---------------
# Derived from the script's location rather than from $PWD, so the dashboard is
# the same dashboard whichever directory it is invoked from.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd -P)"

# -- the environment --------------------------------------------------------
# `MAMBA_ROOT_PREFIX` is unset in this project, so environments live under the
# micromamba binary's own prefix and `micromamba run -n pa-amr` works with no
# setup. `env -u PYTHONPATH` keeps an inherited PYTHONPATH from shadowing the
# environment's own packages.
MAMBA_BIN="${MAMBA_BIN:-micromamba}"
ENV_NAME="${PA_AMR_ENV:-pa-amr}"

if ! command -v "${MAMBA_BIN}" >/dev/null 2>&1; then
    echo "run_dashboard.sh: ${MAMBA_BIN} is not on PATH." >&2
    echo "  This project's environments are micromamba environments; python is" >&2
    echo "  deliberately not on PATH, so the dashboard cannot be started with a" >&2
    echo "  bare 'python -m dashboard'. Set MAMBA_BIN if micromamba lives" >&2
    echo "  somewhere unusual." >&2
    exit 127
fi

# -- go ---------------------------------------------------------------------
# exec, so the server owns the terminal's signals: a Ctrl-C reaches uvicorn
# directly rather than being relayed through a shell that may swallow it.
cd "${REPO_ROOT}"
exec env -u PYTHONPATH "${MAMBA_BIN}" run -n "${ENV_NAME}" \
    python -m dashboard "$@"