"""§4.3 — the oprD endpoint: two honest states, and never a third.

`papipeline/stages/reporting.py` `OPRD_NOT_ON_DISK` states it in the pipeline's
own words: `contracts.py` declares no oprD table in `STAGE_TABLES` or
`INTERNAL_TABLES`, `workflow/Snakefile` declares no oprD output, and `run.py`'s
in-memory `oprd_locus_resolution` is the only copy of those verdicts.

So there is **no file** holding per-isolate oprD verdicts and their evidence,
and this endpoint has exactly two states:

- **verdicts available** — only if a future run writes them. The probe order is
  a contracted path first, then a bundle-supplied table.
- **not produced** — with the reason quoted verbatim, plus every probe path.

**What it must not do.** It must not fall back to the master table's
`chromosomal_mutation` column or to `viz.oprd_status_per_sample`. On the
ten-isolate cohort that function reported `absent` for 8 isolates that all
carry the gene: 34 of the 36 `gene=oprD` CDS features are OprD/OprP/OprQ
paralogs, and 8 of 10 carry the true locus unlabelled. Surfacing that as a
verdict would manufacture the study's central negative.

Verdicts are refused **by name** and never collapse to `absent`
(`adapters/oprd_locus.py:93`). Only `resolved` asserts the locus is present;
every `refused:*` maps to the display state `not_assessed`.
"""

from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from papipeline.adapters.oprd_locus import Verdict
from papipeline.stages.report_tables import NOT_PRODUCED, TableRead
from papipeline.viz import OPRD_STATE_ORDER

#: The reason, quoted. Not paraphrased: this sentence is the pipeline's own and
#: it names three things (the contracts, the Snakefile, the in-memory copy) that
#: a reader needs in order to trust the absence.
REASON = (
    f"{NOT_PRODUCED}: no contracted path exists for the oprD locus verdict. "
    f"`papipeline/execution/contracts.py` declares no oprD table in "
    f"STAGE_TABLES or INTERNAL_TABLES and `workflow/Snakefile` declares no "
    f"oprD output, so `run.py`'s in-memory `oprd_locus_resolution` is the only "
    f"copy of these verdicts and nothing on disk holds them."
)

#: Rendered verbatim. The order is the adapter's.
VERDICTS: Tuple[str, ...] = (
    Verdict.RESOLVED,
    Verdict.NO_HIT,
    Verdict.INSUFFICIENT,
    Verdict.AMBIGUOUS,
    Verdict.AMBIGUOUS_SECOND,
)

#: Display states, from `viz.OPRD_STATE_ORDER`.
DISPLAY_STATES: Tuple[str, ...] = tuple(OPRD_STATE_ORDER)

#: The columns a contracted oprD table would carry, if one existed. Used to
#: recognise a future table rather than to invent a schema for today's.
EXPECTED_COLUMNS: Tuple[str, ...] = (
    "sample_id",
    "verdict",
    "lesion_type",
    "position",
    "truncation_aa",
    "identity_pct",
    "coverage_pct",
)

#: The evidence families DESIGN §4.3 requires. A future contracted table that
#: carries them must be able to reach the UI, so they are named here and the
#: `evidence` block carries every additional column the row holds.
EVIDENCE_FAMILIES: Tuple[str, ...] = (
    "repeat_flag",
    "homopolymer_flag",
    "compensating_indel",
    "tblastn_summary",
)


def probe_paths(root) -> List[Any]:
    """Every path probed, in order.

    Contracted first, then bundle-supplied. A bundle may have collected the
    table under a name the pipeline does not declare, and the endpoint probes
    for that rather than assuming it does not exist.
    """
    from pathlib import Path

    root = Path(root)
    return [
        root / "intermediate" / "stages" / "oprd.tsv",
        root / "intermediate" / "stages" / "oprD.tsv",
        root / "intermediate" / "stages" / "oprd_locus.tsv",
        root / "intermediate" / "stages" / "17_oprd.tsv",
        root / "02_stage_outputs" / "oprd.tsv",
        root / "02_stage_outputs" / "oprd_locus.tsv",
        root / "03_report" / "oprd.tsv",
        root / "04_run_info" / "oprd.tsv",
        root / "05_validation" / "oprd.tsv",
    ]


def _read_candidate(path) -> Tuple[Optional[TableRead], Dict[str, Any]]:
    """Try one candidate path, and say why it did or did not yield verdicts."""
    from pathlib import Path

    from papipeline.stages.report_tables import read_table

    path = Path(path)
    if not path.exists():
        return None, {
            "path": str(path),
            "ok": False,
            "reason": "no file at this path",
        }
    read = read_table("oprd", path, required_columns=("sample_id", "verdict"))
    if not read.present:
        return None, {"path": str(path), "ok": False, "reason": read.reason}
    unknown = sorted(
        {
            str(row.get("verdict"))
            for row in read.rows
            if str(row.get("verdict")) not in VERDICTS
        }
    )
    if unknown:
        # A verdict outside the vocabulary is refused rather than coerced. A
        # file that speaks a different language here is not a source of
        # verdicts, and rendering its values would invent a claim.
        return None, {
            "path": str(path),
            "ok": False,
            "reason": (
                f"{path} carries verdict value(s) {unknown}, which are not in "
                f"the {len(VERDICTS)} of `adapters.oprd_locus.Verdict`. A "
                f"verdict outside that vocabulary is not coerced into one."
            ),
        }
    return read, {"path": str(path), "ok": True, "reason": f"read {read.n_rows} verdicts"}


def display_state_for(verdict: Optional[str], lesion_type: Optional[str] = None) -> str:
    """The display state a verdict maps to.

    Only `resolved` asserts presence. Every refusal is `not_assessed` — never
    `absent`, which is reserved for a confirmed structural lesion. Collapsing
    the two is what turned 8 gene-carrying isolates into `absent` on the
    ten-isolate cohort.
    """
    if verdict == Verdict.RESOLVED:
        if lesion_type:
            return "disrupted"
        return "intact"
    return "not_assessed"


def collect(root) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], Optional[str]]:
    """Probe for verdicts.

    Returns `(items, probes, reason)`. `items` is empty and `reason` is the
    quoted `OPRD_NOT_ON_DISK` sentence unless a source was found.
    """
    probes: List[Dict[str, Any]] = []
    for index, path in enumerate(probe_paths(root), start=1):
        read, record = _read_candidate(path)
        record["n"] = index
        probes.append({"n": index, **record})
        if read is None:
            continue
        items: List[Dict[str, Any]] = []
        for row in read.rows:
            verdict = str(row.get("verdict"))
            lesion = row.get("lesion_type")
            items.append(
                {
                    "sample_id": str(row.get("sample_id")),
                    "verdict": verdict,
                    "lesion_type": lesion,
                    "position": _int_or_none(row.get("position")),
                    "truncation_aa": _int_or_none(row.get("truncation_aa")),
                    "identity_pct": _float_or_none(row.get("identity_pct")),
                    "coverage_pct": _float_or_none(row.get("coverage_pct")),
                    "display_state": display_state_for(verdict, lesion),
                    # Every additional column the verdict row carries, so the
                    # repeat / compensating-indel / tblastn evidence DESIGN §4.3
                    # requires reaches the UI. A contracted table that names
                    # them differently still arrives: the whole row is carried,
                    # not just the names this build happens to know.
                    "evidence": {
                        k: v
                        for k, v in row.items()
                        if k not in ("sample_id", "verdict") and v is not None
                    }
                    or None,
                }
            )
        return items, probes, None
    return [], probes, REASON


def _int_or_none(value: Any) -> Optional[int]:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _float_or_none(value: Any) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def by_sample(items: Sequence[Mapping[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """The verdicts keyed on sample, for the isolate view's join.

    Empty when nothing was produced, so every isolate reads `not_assessed`
    rather than a fabricated state.
    """
    return {str(item["sample_id"]): dict(item) for item in items if item.get("sample_id")}


__all__ = [
    "DISPLAY_STATES",
    "EXPECTED_COLUMNS",
    "REASON",
    "VERDICTS",
    "by_sample",
    "collect",
    "display_state_for",
    "probe_paths",
]