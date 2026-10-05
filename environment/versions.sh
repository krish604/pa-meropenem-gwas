#!/usr/bin/env bash
# =============================================================================
# Record the versions of every tool and database this pipeline depends on.
#
# Run this BEFORE and AFTER any analysis and keep both outputs. A run whose
# provenance cannot be reproduced is not a result, it is an anecdote.
#
#   bash environment/versions.sh                 # human-readable
#   bash environment/versions.sh --json          # machine-readable
#   bash environment/versions.sh --out FILE      # write to FILE
#
# This script only *reads* versions. It never installs, updates or downloads
# anything. Provision databases out of band, then record them here.
# =============================================================================

set -uo pipefail

JSON=0
OUT=""
for arg in "$@"; do
    case "$arg" in
        --json) JSON=1 ;;
        --out) shift; OUT="${1:-}" ;;
        -h|--help)
            sed -n '2,14p' "$0" | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        *) echo "unknown argument: $arg" >&2; exit 2 ;;
    esac
done

# --- tools the pipeline invokes, in stage order -------------------------------
TOOLS=(
    # core / workflow
    python3 snakemake pytest
    # stage 2
    bakta
    # stage 3
    mlst
    # stage 4
    amrfinder rgi
    # stage 7 / 8
    minimap2 samtools makeblastdb blastn
    # stage 9
    panaroo
    # stage 10
    snp-sites iqtree2 iqtree
    # stage 12
    pyseer
)

# Flag to try per tool. Several tools reject --version and require -version.
version_of() {
    local tool="$1"
    local path
    path="$(command -v "$tool" 2>/dev/null)" || return 1
    [ -n "$path" ] || return 1
    local flag
    for flag in --version -version; do
        local out
        out="$("$tool" "$flag" 2>&1 | head -n 1)"
        if [ -n "$out" ] && ! printf '%s' "$out" | grep -qiE 'usage|unrecognized|invalid|unknown option|error:'; then
            printf '%s' "$out"
            return 0
        fi
    done
    # Some tools print their version with no flag at all.
    printf '%s' "$("$tool" 2>&1 | head -n 1)"
}

emit_text() {
    echo "================================================================"
    echo " Pipeline environment versions"
    echo " recorded: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    echo " host     : $(uname -srm)"
    echo "================================================================"
    echo
    printf '%-16s %-10s %s\n' "TOOL" "STATUS" "VERSION"
    printf '%-16s %-10s %s\n' "----------------" "----------" "-------"
    for tool in "${TOOLS[@]}"; do
        local_version="$(version_of "$tool" 2>/dev/null || true)"
        if [ -z "$local_version" ]; then
            printf '%-16s %-10s %s\n' "$tool" "MISSING" "-"
        else
            printf '%-16s %-10s %s\n' "$tool" "present" "$local_version"
        fi
    done
    echo
    echo "DATABASES"
    echo "----------"
    echo "Database versions are NOT auto-detected. They are pinned in"
    echo "config/references.tsv and must be provisioned out of band."
    echo
    if [ -f config/references.tsv ]; then
        # Print reference_id, database, database_version, version_status.
        awk -F'\t' '
            /^#/ || NF < 6 { next }
            !seen { for (i=1;i<=NF;i++) col[$i]=i; seen=1 }
            { printf "  %-24s %-28s %-14s %s\n", $col["reference_id"], $col["database"], $col["database_version"], $col["version_status"] }
        ' config/references.tsv
        echo
        unpinned="$(awk -F'\t' '!/^#/ && $6=="unpinned"' config/references.tsv | wc -l | tr -d ' ')"
        if [ "$unpinned" != "0" ]; then
            echo "WARNING: $unpinned reference(s) are UNPINNED."
            echo "         Pin every version before any real analysis"
            echo "         (docs/scientific_rules.md, rule 8)."
        else
            echo "All references are pinned."
        fi
    else
        echo "  config/references.tsv not found (run from the pipeline root)."
    fi
    echo
    echo "NOTE: nothing was installed, updated or downloaded by this script."
}

emit_json() {
    local inner
    inner="$(mktemp)"
    emit_text > "$inner"
    python3 - "$inner" <<'PY'
import json, re, subprocess, sys, datetime, platform

text = open(sys.argv[1]).read()
tools, dbs = {}, {}
section = None
for line in text.splitlines():
    if line.startswith("TOOL"):
        section = "tools"; continue
    if line.startswith("DATABASES"):
        section = "dbs"; continue
    if section == "tools":
        m = re.match(r"^(\S+)\s+(present|MISSING)\s+(\S.*)?$", line)
        if m:
            tools[m.group(1)] = {
                "status": m.group(2),
                "version": (m.group(3) or "").strip() or None,
            }
    elif section == "dbs":
        m = re.match(r"^\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$", line)
        if m:
            dbs[m.group(1)] = {
                "database": m.group(2),
                "database_version": m.group(3),
                "version_status": m.group(4),
            }

print(json.dumps({
    "recorded_utc": datetime.datetime.now(datetime.timezone.utc)
        .strftime("%Y-%m-%dT%H:%M:%SZ"),
    "host": platform.platform(),
    "python": sys.version.split()[0],
    "tools": tools,
    "references": dbs,
    "note": "Versions only. Nothing installed, updated or downloaded.",
}, indent=2))
PY
    rm -f "$inner"
}

if [ "$JSON" -eq 1 ]; then
    OUTPUT="$(emit_json)"
    if [ -n "$OUT" ]; then printf '%s\n' "$OUTPUT" > "$OUT"; echo "wrote $OUT";
    else printf '%s\n' "$OUTPUT"; fi
else
    if [ -n "$OUT" ]; then emit_text > "$OUT"; echo "wrote $OUT";
    else emit_text; fi
fi
