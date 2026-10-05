#!/usr/bin/env python3
"""Pre-flight check: everything that could fail after hours of runtime.

Run this before launching a large cohort. It verifies the tools, the
databases, the free disk, the configuration and the cohort itself, and
prints a single GO / NO-GO verdict with the reason.

The point is to fail in seconds here rather than part-way through a
multi-hour annotation sweep. Nothing is computed and no output is written;
this only inspects.

Usage:
    python scripts/preflight.py
    python scripts/preflight.py --expect 200
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Tuple

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

RESULTS = ROOT / "results" / "pilot100"
ANNOTATION = RESULTS / "annotation"

#: tool -> (minimum version, why it matters)
TOOLS = [
    ("bakta", "1.12", "stage 2 annotation"),
    ("amrfinder", "4.0", "stage 3 AMR detection"),
    ("mlst", "2.19", "stage 5 MLST"),
    ("pyseer", "1.1", "stage 10 GWAS"),
]

#: Tools that may live outside the env, with the version floor that matters.
#: samtools/bcftools 0.1.19 (2014) exist in pilot100 and fail silently on
#: modern SAM, so the floor is 1.0 and PATH order is not trusted.
HTS_TOOLS = [
    ("minimap2", "2.17", "stage 9 alignment"),
    ("samtools", "1.10", "stage 9 sort/index"),
    ("bcftools", "1.10", "stage 9 variant calling"),
    ("FastTree", "2.1", "stage 9 tree"),
]

#: Bytes of annotation output kept per genome after pruning, and the raw
#: figure before it, measured on this data (84 MB raw, 8.6 MB pruned).
RAW_BYTES_PER_GENOME = 84 * 1024 * 1024
PRUNED_BYTES_PER_GENOME = 9 * 1024 * 1024
PHYLOGENY_BYTES_PER_GENOME = 9 * 1024 * 1024


class Report:
    def __init__(self) -> None:
        self.lines: List[str] = []
        self.failures: List[str] = []
        self.warnings: List[str] = []

    def ok(self, label: str, detail: str = "") -> None:
        line = f"  [ ok ] {label:<28} {detail}"
        self.lines.append(line)
        print(line)

    def warn(self, label: str, detail: str) -> None:
        self.warnings.append(f"{label}: {detail}")
        line = f"  [warn] {label:<28} {detail}"
        self.lines.append(line)
        print(line)

    def fail(self, label: str, detail: str) -> None:
        self.failures.append(f"{label}: {detail}")
        line = f"  [FAIL] {label:<28} {detail}"
        self.lines.append(line)
        print(line)


def version_of(name: str, path: str) -> Optional[str]:
    """Best-effort version string.

    Not every tool implements ``--version``: FastTree prints usage and exits
    non-zero, and the 2014 samtools stub rejects the flag as an unknown
    command. Several probes are tried and the first version-looking token
    wins, so a tool is not reported as missing merely because it spelled
    its version flag differently.
    """
    probes = [[path, "--version"], [path, "-version"], [path, "version"], [path]]
    for probe in probes:
        try:
            proc = subprocess.run(probe, capture_output=True, text=True, timeout=60)
        except Exception:  # noqa: BLE001
            continue
        text = ((proc.stdout or "") + " " + (proc.stderr or "")).strip()
        for token in text.replace("(", " ").replace(")", " ").split():
            cleaned = token.strip(",;")
            if "." not in cleaned:
                continue
            head = cleaned.split(".")[0]
            if head.isdigit():
                # minimap2 reports "2.31-r1302"; the numeric prefix is enough.
                return cleaned
    return None


def _load_config(machine: str = "laptop"):
    """The overlay this pre-flight is checking. Tools resolve through it."""
    from papipeline.config.loader import load_config

    return load_config(ROOT / "config" / "science.yaml", machine=machine)


def candidates(name: str, config=None) -> List[str]:
    """Every location worth trying, best first.

    PATH is authoritative and comes first. After that, only the search
    directories this machine's overlay declares.

    The PATH-first ordering is not cosmetic: the pilot100 environment ships
    samtools/bcftools 0.1.19 from 2014, and those can shadow a working build,
    so the version floor below is what actually decides. The ordering only
    decides which candidate is tried first.

    An earlier revision of this function also tried a literal macOS Homebrew
    path and /usr/local/bin. Neither exists on the analysis machine, so a
    genuinely missing tool resolved to nothing and was reported as present -
    a pre-flight check that cannot fail is not a pre-flight check.
    """
    if config is None or config.machine is None:
        resolved = shutil.which(name)
        return [resolved] if resolved else []
    return list(config.machine.tool_candidates(name))


def check_tools(report: Report, config=None) -> None:
    def tup(v: str) -> Tuple[int, ...]:
        parts = []
        for chunk in v.split(".")[:2]:
            digits = "".join(ch for ch in chunk if ch.isdigit())
            if not digits:
                break
            parts.append(int(digits))
        return tuple(parts) or (0,)

    for name, floor, why in TOOLS + HTS_TOOLS:
        paths = candidates(name)
        if not paths:
            report.fail(name, f"not found ({why})")
            continue
        best: Optional[Tuple[Tuple[int, ...], str, str]] = None
        for path in paths:
            found = version_of(name, path) or ""
            if not found:
                continue
            score = tup(found)
            if best is None or score > best[0]:
                best = (score, found, path)
        if best is None:
            report.warn(name, f"found at {paths[0]} but no version could be read")
            continue
        score, found, path = best
        if score < tup(floor):
            report.fail(name, f"newest found is {found} < {floor} required ({why})")
        else:
            skipped = [q for q in paths if q != path and (version_of(name, q) or "")]
            note = ""
            if skipped:
                note = "  (ignored older: " + ", ".join(q.split("/")[-2] + "/" + q.split("/")[-1] for q in skipped) + ")"
            report.ok(name, f"{found}  [{path}]{note}")


def check_databases(report: Report) -> None:
    try:
        from scripts.run_annotation_sweep import default_db, database_ready
        db = default_db()
        ready, why = database_ready(db)
        if ready:
            report.ok("bakta database", f"ready ({db.name})")
        else:
            report.fail("bakta database", why)
    except Exception as exc:  # noqa: BLE001
        report.fail("bakta database", f"{type(exc).__name__}: {exc}")

    amr_db = Path.home() / "micromamba_pilot/envs/pilot100/share/amrfinderplus/latest"
    version_file = amr_db / "version.txt"
    if not amr_db.exists():
        report.fail("amrfinderplus db", f"missing {amr_db}")
    elif not version_file.exists():
        report.fail("amrfinderplus db", "version.txt absent; database looks partial")
    else:
        report.ok("amrfinderplus db", f"{version_file.read_text().strip()}")


def check_config(report: Report) -> None:
    try:
        from papipeline.config.loader import load_config
        from papipeline.models import RunMode
        config = load_config(ROOT / "config" / "science.yaml", machine=args.machine)
        report.ok("config", f"primary={config.antibiotics[0]} thresholds={config.qc.min_assembly_size}-{config.qc.max_assembly_size}")
    except Exception as exc:  # noqa: BLE001
        report.fail("config", f"{type(exc).__name__}: {exc}")


def check_cohort(report: Report, expect: Optional[int]) -> Tuple[int, int]:
    manifest = RESULTS / "pilot100_manifest.tsv"
    if not manifest.exists():
        report.fail("cohort manifest", f"missing {manifest}; run run_pilot100.py --prepare-only")
        return 0, 0
    lines = [l for l in manifest.read_text().splitlines()
             if l.strip() and not l.startswith("#")]
    header = lines[0].split("\t")
    rows = [dict(zip(header, l.split("\t"))) for l in lines[1:]]
    included = [r for r in rows
                if str(r.get("analysis_included", "")).upper() == "TRUE"]
    report.ok("cohort manifest", f"{len(included)} included of {len(rows)} rows")
    if expect is not None and len(included) != expect:
        report.warn("cohort size",
                    f"expected {expect} but manifest has {len(included)}; "
                    f"re-run prepare if this is not intended")
    missing = [r for r in included if not Path(r["Genome_path"]).exists()]
    if missing:
        report.fail("genome files", f"{len(missing)} listed in the manifest are absent")
    else:
        report.ok("genome files", f"all {len(included)} present")
    empty = [r for r in included
             if Path(r["Genome_path"]).exists()
             and Path(r["Genome_path"]).stat().st_size == 0]
    if empty:
        report.fail("genome files", f"{len(empty)} are zero bytes")
    return len(included), len(missing)


def check_disk(report: Report, n_genomes: int, done: int) -> None:
    pending = max(0, n_genomes - done)
    usage = shutil.disk_usage(ROOT)
    free = usage.free
    need = pending * (PRUNED_BYTES_PER_GENOME + PHYLOGENY_BYTES_PER_GENOME)
    # 20% headroom, because a full disk mid-sweep is unrecoverable.
    want = int(need * 1.2)
    print(f"  [disk] free={free/1e9:.1f} GB  pending={pending}  "
          f"need~{need/1e9:.1f} GB (with headroom {want/1e9:.1f} GB)")
    if pending == 0:
        report.ok("disk space", f"{free/1e9:.1f} GB free, nothing pending")
    elif free < need:
        report.fail("disk space",
                    f"only {free/1e9:.1f} GB free, need ~{need/1e9:.1f} GB for {pending} genomes")
    elif free < want:
        report.warn("disk space",
                    f"{free/1e9:.1f} GB free vs ~{need/1e9:.1f} GB needed; "
                    f"headroom is thin")
    else:
        report.ok("disk space", f"{free/1e9:.1f} GB free, ~{need/1e9:.1f} GB needed")
    print(f"  [note] unpruned Bakta output would be "
          f"{pending * RAW_BYTES_PER_GENOME/1e9:.1f} GB; pruning keeps it to "
          f"{pending * PRUNED_BYTES_PER_GENOME/1e9:.1f} GB")


def check_annotated(report: Report) -> int:
    if not ANNOTATION.exists():
        report.warn("annotation", "no annotation directory yet")
        return 0
    try:
        from scripts.run_annotation_sweep import is_complete
        dirs = [p for p in ANNOTATION.iterdir() if p.is_dir()]
        done = sum(1 for p in dirs if is_complete(p))
        report.ok("annotation", f"{done}/{len(dirs)} sample directories complete")
        return done
    except Exception as exc:  # noqa: BLE001
        report.warn("annotation", f"could not verify: {exc}")
        return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expect", type=int, default=None,
                        help="Warn if the prepared cohort is not this size.")
    parser.add_argument("--machine", default="laptop",
                        help="Machine overlay: laptop or bigmachine.")
    args = parser.parse_args(argv)

    try:
        config = _load_config(args.machine)
    except Exception as exc:  # noqa: BLE001 - reported, not raised: a pre-flight
        # that crashes tells you less than one that reports why it cannot run.
        print(f"Pre-flight\n\nERROR: cannot load the configuration for machine "
              f"{args.machine!r}: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1

    report = Report()
    print("Pre-flight\n")
    print("Tools")
    check_tools(report, config)
    print("\nDatabases")
    check_databases(report)
    print("\nConfiguration")
    check_config(report)
    print("\nCohort")
    n, _missing = check_cohort(report, args.expect)
    done = check_annotated(report)
    print("\nCapacity")
    check_disk(report, n, done)

    print("\n" + "=" * 66)
    if report.failures:
        print(f"NO-GO - {len(report.failures)} blocking problem(s):")
        for f in report.failures:
            print(f"  - {f}")
    else:
        print("GO - all blocking checks passed.")
    if report.warnings:
        print(f"\n{len(report.warnings)} warning(s):")
        for w in report.warnings:
            print(f"  - {w}")
    if n and done < n:
        pending = n - done
        per = 353.0
        print(f"\nAnnotation remaining: {pending} genome(s).")
        print(f"  At the measured 353 s/genome, 2 jobs x 4 threads: "
              f"~{pending * per / 2 / 3600:.1f} h wall.")
    print("=" * 66)
    return 1 if report.failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
