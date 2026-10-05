#!/usr/bin/env python3
"""Extract imipenem AST calls for the smoke cohort into a smoke-scoped table.

**Why this exists.** `data/phenotype/` is the real cohort's directory and is
empty. Stage 11 refuses without phenotype data, and stage 12 (GWAS) refuses
below `gwas.min_samples_per_group` per outcome group - so a smoke cohort of ten
isolates that were all imipenem-R produced *correct* data that still failed
stage 12. A single-group cohort cannot be tested. The composition was fixed
instead; this script writes the phenotype table for whatever cohort is in
`smoke_genome_dir`.

**The calls are measured, not assigned.** They come from the `AST phenotypes`
column of `PDC_essential.tsv`, which holds laboratory AST results per isolate.
Nothing here defaults, infers, or falls back: an isolate with no imipenem entry
is an error, because silently omitting it would shrink the group and make the
GWAS floor fail for a reason that looks like a pipeline bug.

The output is written to `db/smoke_phenotypes/` - a directory of its own, not
`data/phenotype/`. A bounded run's table must not sit in the full cohort's
place.

Re-runnable and deterministic: same inputs, same bytes out.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional

REPO = Path(__file__).resolve().parents[2]

#: The PDC isolate table. Read directly; `PDC_essential.tsv` is the source of
#: truth for both cohort membership and AST calls, and copying it into a
#: metadata file would create a second one.
PDC_TABLE = REPO / "PDC_essential.tsv"

#: The full cohort's phenotype directory. Read-only to this script, always.
REAL_PHENOTYPE_DIR = REPO / "data" / "phenotype"

#: The column holding per-drug AST results, e.g. ``amikacin=R,imipenem=S``.
AST_COLUMN = "AST phenotypes"

#: Entries this script will refuse to write rather than guess at. `I` and `ND`
#: are real categories and the pipeline accepts them, but the GWAS model is
#: binary R-vs-S, so writing an intermediate as either side would fabricate a
#: call the laboratory did not make.
SUPPORTED_CALLS = {"R", "S"}

HEADER = """\
# SMOKE RUN PHENOTYPE - measured values, not synthesised
#
# Extracted from `{pdc}` column `{column}` by
# `scripts/smoke/build_smoke_phenotype.py`. Every call below is a laboratory AST
# result for that isolate; none is assigned, defaulted or inferred.
#
# Scope: the isolates in `db/smoke_genomes/` only. This is NOT the full cohort's
# phenotype table - that is `data/phenotype/`, which this script does not write
# and must not be confused with. A report built from this data carries the
# SMOKE TEST marker.
#
# `MIC` is empty throughout: the AST gives a categorical call only, and
# converting R/S into a number would fabricate a measurement."""


def imipenem_call(ast_field: str) -> Optional[str]:
    """The isolate's imipenem call from its AST field, or None if absent.

    The exported field wraps its value in quotes and lists drugs comma-separated,
    so a plain split is not enough and the drug names contain hyphens.
    """
    text = (ast_field or "").strip().strip('"')
    match = re.search(r"(?:^|,)imipenem=([^,]+)", text)
    return match.group(1).strip() if match else None


def refuses_to_write_into(out: Path, real_dir: Path) -> Optional[str]:
    """Why `out` is the forbidden directory, or None if it is fine.

    Split out from `main` and given the forbidden directory as an argument so it
    can be tested without pointing the test at the real one. Tested by
    *calling* it with the real path instead: were this guard ever dropped, that
    test would write a ten-isolate table into `data/phenotype/` - the exact
    corruption the guard exists to prevent, caused by the test that checks it.

    Resolved before comparing, so `db/./smoke_phenotypes` and a trailing slash
    cannot slip past.
    """
    if out.resolve() == real_dir.resolve():
        return (
            f"refusing to write into {real_dir.name}/: that is the full cohort's "
            "phenotype directory, and a 10-isolate smoke table written there "
            "would be indistinguishable from the real one"
        )
    return None


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--genomes", type=Path, default=REPO / "db" / "smoke_genomes",
        help="Directory of smoke-cohort assemblies (stem = isolate).",
    )
    parser.add_argument(
        "--out", type=Path, default=REPO / "db" / "smoke_phenotypes",
        help="Smoke-scoped phenotype directory. Never data/phenotype/.",
    )
    parser.add_argument("--antibiotic", default="imipenem")
    parser.add_argument(
        "--pdc", type=Path, default=PDC_TABLE,
        help="PDC isolate table to read calls from. Injectable so the refusal "
             "paths can be tested without editing the real one.",
    )
    args = parser.parse_args(argv)
    pdc = args.pdc

    refusal = refuses_to_write_into(args.out, REAL_PHENOTYPE_DIR)
    if refusal:
        print(refusal, file=sys.stderr)
        return 2

    cohort = sorted(p.stem for p in args.genomes.glob("*.fna"))
    if not cohort:
        print(f"no assemblies in {args.genomes}", file=sys.stderr)
        return 2
    wanted = set(cohort)

    calls: Dict[str, str] = {}
    with pdc.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            isolate = row.get("Isolate")
            if isolate in wanted:
                calls[isolate] = imipenem_call(row.get(AST_COLUMN, ""))

    missing = sorted(wanted - set(calls))
    if missing:
        print(
            f"{len(missing)} cohort isolate(s) have no row in {pdc.name}: "
            f"{', '.join(missing)}",
            file=sys.stderr,
        )
        return 2

    unmeasured = sorted(i for i in cohort if calls.get(i) is None)
    if unmeasured:
        print(
            f"{len(unmeasured)} isolate(s) have no {args.antibiotic} AST entry: "
            f"{', '.join(unmeasured)}. Refusing rather than defaulting them - an "
            "omitted isolate would shrink an outcome group and make the GWAS "
            "floor fail for a reason that looks like a pipeline fault.",
            file=sys.stderr,
        )
        return 2

    unsupported = sorted(
        i for i in cohort if calls[i] not in SUPPORTED_CALLS
    )
    if unsupported:
        detail = ", ".join(f"{i}={calls[i]}" for i in unsupported)
        print(
            f"{len(unsupported)} isolate(s) carry a call this binary model does "
            f"not use: {detail}. Refusing rather than assigning them to R or S.",
            file=sys.stderr,
        )
        return 2

    args.out.mkdir(parents=True, exist_ok=True)
    path = args.out / f"{args.antibiotic}_phenotype.tsv"
    lines = [
        HEADER.format(pdc=pdc.name, column=AST_COLUMN),
        "sample_id\tantibiotic\tphenotype\tMIC\tMIC_unit\tsource",
    ]
    for isolate in cohort:
        lines.append(
            f"{isolate}\t{args.antibiotic}\t{calls[isolate]}\t.\t.\t"
            f"PDC_{AST_COLUMN.replace(' ', '_')}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    counts: Dict[str, int] = {}
    for isolate in cohort:
        counts[calls[isolate]] = counts.get(calls[isolate], 0) + 1
    print(f"wrote {path}")
    print(f"  cohort      : {len(cohort)} isolates")
    print(f"  composition : {dict(sorted(counts.items()))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
