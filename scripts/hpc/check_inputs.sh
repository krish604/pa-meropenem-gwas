#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# check_inputs.sh -- input gate. Run before any compute.
#
# Verifies three things and exits non-zero on any failure, so a bad cohort is
# caught before a single CPU-minute is spent:
#
#   1. NO ISOLATE IS REUSED. The sample list may not contain the same
#      identifier twice. Duplicate IDs would silently pool the same assembly
#      into the cohort twice, which corrupts every downstream count --
#      allele frequencies, case/control splits and the GWAS denominator all
#      assume each row is a distinct organism.
#
#   2. The cohort size matches EXPECTED_ISOLATES (default 900), so a truncated
#      download or a botched concatenation is caught here rather than
#      discovered halfway through the run.
#
#   3. Every listed isolate has an assembly FASTA on disk.
#
# Usage:  scripts/hpc/check_inputs.sh [--size-only]
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

LIST="$(sample_list)"
[[ -f "$LIST" ]] || die "no sample list at $LIST
Create it with one isolate ID per line, e.g.  > config/samples.txt"

# --- 1. no isolate may appear twice ----------------------------------------
# Trailing whitespace and CR (a common result of editing on Windows) are
# normalised first: 'PDT001 ' and 'PDT001' are the same isolate and would
# otherwise slip past a plain duplicate check.
#
# Read line-by-line rather than via `mapfile`, which is a bash 4 builtin and
# is absent from the bash 3.2 that macOS still ships. A researcher validating
# a roster on a laptop before uploading it must not need bash 4 for that.
IDS=()
while IFS= read -r RAW || [[ -n "$RAW" ]]; do
    RAW="${RAW%$'\r'}"
    RAW="${RAW%"${RAW##*[![:space:]]}"}"   # strip trailing whitespace
    [[ -n "$RAW" ]] && IDS+=("$RAW")
done < "$LIST"
(( ${#IDS[@]} )) || die "sample list is empty: $LIST"

TOTAL=${#IDS[@]}
UNIQUE=$(printf '%s\n' "${IDS[@]}" | sort -u | wc -l | tr -d ' ')

if (( TOTAL != UNIQUE )); then
    printf '\nDUPLICATE ISOLATES DETECTED -- refusing to run.\n' >&2
    printf '%d lines but %d unique IDs. Offenders:\n\n' "$TOTAL" "$UNIQUE" >&2
    printf '%s\n' "${IDS[@]}" | sort | uniq -d | sed 's/^/    /' >&2
    printf '\nEach isolate must appear exactly once. Remove the repeats\n' >&2
    printf 'from config/samples.txt and re-run this check.\n' >&2
    exit 1
fi
log "no duplicate isolates: $TOTAL unique IDs"

# --- 2. expected cohort size ------------------------------------------------
if (( TOTAL != EXPECTED_ISOLATES )); then
    printf '\nCOHORT SIZE MISMATCH -- refusing to run.\n' >&2
    printf '  samples.txt : %d unique isolates\n' "$TOTAL" >&2
    printf '  expected    : %d (EXPECTED_ISOLATES in config/run.conf)\n' \
           "$EXPECTED_ISOLATES" >&2
    printf '\nIf %d is correct for this cohort, edit config/run.conf.\n' \
           "$EXPECTED_ISOLATES" >&2
    printf 'If not, your sample list is incomplete.\n' >&2
    exit 1
fi
log "cohort size confirmed: $TOTAL isolates"

[[ "${1:-}" == "--size-only" ]] && { log "check_inputs: OK (--size-only)"; exit 0; }

# --- 3. every isolate has an assembly ---------------------------------------
ASSEMBLIES="$REPO_ROOT/assemblies"
[[ -d "$ASSEMBLIES" ]] || die "no assemblies/ directory at $ASSEMBLIES
Expected FASTA files at assemblies/<isolate_id>.fna (or .fasta / .fa)"

MISSING=0
for id in "${IDS[@]}"; do
    found=0
    for ext in fna fasta fa fna.gz fasta.gz; do
        [[ -f "$ASSEMBLIES/$id.$ext" ]] && { found=1; break; }
    done
    if (( ! found )); then
        printf '    missing: %s\n' "$id" >&2
        MISSING=$((MISSING + 1))
        (( MISSING >= 20 )) && { printf '    ... (further omissions suppressed)\n' >&2; break; }
    fi
done

if (( MISSING )); then
    printf '\n%d isolate(s) listed in samples.txt have no assembly FASTA.\n' \
           "$MISSING" >&2
    printf 'Assemblies were expected under %s\n' "$ASSEMBLIES" >&2
    exit 1
fi
log "all $TOTAL assemblies present"
log "check_inputs: OK"
