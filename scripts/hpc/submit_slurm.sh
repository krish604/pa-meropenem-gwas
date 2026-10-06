#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# submit_slurm.sh -- run the pipeline as batch jobs on an HPC cluster.
#
# The interactive script (run_all.sh) is fine on a workstation but not on a
# cluster: compute nodes have no terminal to prompt you in. This script splits
# the run into the two jobs the consent gate already implies.
#
#   scripts/hpc/submit_slurm.sh phase1     submit stages 1-11
#   scripts/hpc/submit_slurm.sh gwas       submit stages 12-16 (needs approval)
#   scripts/hpc/submit_slurm.sh status     show submitted jobs
#   scripts/hpc/submit_slurm.sh consent    record approval from a login node
#
# Typical session on the cluster:
#
#   1.  scripts/hpc/submit_slurm.sh phase1        -> job submitted
#   2.  ... days later ...
#   3.  inspect results/ by eye
#   4.  scripts/hpc/submit_slurm.sh consent       -> writes the approval file
#   5.  scripts/hpc/submit_slurm.sh gwas          -> job submitted
#
# Resources come from config/run.conf (CORES, MEM_GB, TIME_*). No scheduler
# numbers are written here.
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

MODE="${1:-}"
[[ -n "$MODE" ]] || { sed -n '2,25p' "$0"; exit 1; }

# Write a job file and submit it. $1 = label, $2 = script path.
submit() {
    local label="$1" runner="$2" wall cores mem
    case "$label" in
        phase1) wall="$TIME_PHASE1" ;;
        gwas)   wall="$TIME_PHASE2" ;;
    esac
    cores="$CORES"
    mem="${MEM_GB}G"

    local job="$RESULTS_ROOT/slurm_$label.sh"
    cat > "$job" <<EOF
#!/bin/bash
#SBATCH --job-name=pa900_$label
#SBATCH --output=$RESULTS_ROOT/logs/slurm_$label-%j.out
#SBATCH --error=$RESULTS_ROOT/logs/slurm_$label-%j.err
#SBATCH --cpus-per-task=$cores
#SBATCH --mem=$mem
#SBATCH --time=$wall
set -euo pipefail
cd "$REPO_ROOT"
exec "$runner"
EOF

    mkdir -p "$RESULTS_ROOT/logs"
    log "submitting $label  (cores=$cores mem=$mem time=$wall)"
    sbatch "$job"
}

case "$MODE" in
    phase1)
        "$REPO_ROOT/scripts/hpc/check_inputs.sh"
        submit phase1 "$REPO_ROOT/scripts/hpc/run_phase1.sh"
        ;;
    gwas)
        CONSENT_PATH="$REPO_ROOT/$CONSENT_FILE"
        if [[ ! -f "$CONSENT_PATH" || "$(cat "$CONSENT_PATH")" != "APPROVED" ]]; then
            printf '\nGWAS NOT APPROVED -- nothing submitted.\n\n' >&2
            printf 'The permission gate is checked here, not just in the script,\n' >&2
            printf 'because a compute node cannot prompt you.\n\n' >&2
            printf 'Approve from a login node first:\n' >&2
            printf '    scripts/hpc/submit_slurm.sh consent\n\n' >&2
            exit 3
        fi
        submit gwas "$REPO_ROOT/scripts/hpc/run_gwas.sh"
        ;;
    consent)
        cat <<EOF
This records your approval to run the GWAS stages (12-16) on:

    $REPO_ROOT/$CONSENT_FILE

Only do this after you have looked at the phase-1 results.

EOF
        read -r -p "Approve GWAS on this cohort? [y/N] " ANSWER
        if [[ "${ANSWER:-}" =~ ^([yY]|yes|YES)$ ]]; then
            printf 'APPROVED\n' > "$REPO_ROOT/$CONSENT_FILE"
            log "approval recorded -- you can now submit: scripts/hpc/submit_slurm.sh gwas"
        else
            log "not approved; nothing written"
            exit 0
        fi
        ;;
    status)
        command -v squeue >/dev/null 2>&1 || die "squeue not found"
        squeue -u "${USER:-$(id -un)}" -o '%.10i %.14j %.8T %.10M %.6D %R'
        ;;
    *) die "unknown mode: $MODE (use: phase1 | gwas | consent | status)" ;;
esac
