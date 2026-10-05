"""Stage 1 - genome validation.

Computes assembly metrics in pure Python (no external tool required) and
applies the configured thresholds. A sample that fails a threshold is
*flagged*, never dropped: exclusion is an analysis decision, not a side
effect of QC.

The species / contamination / completeness / duplicate hooks are declared
here as extension points. They are disabled by default and each reports
``None`` (not ``pass``) until a tool is plugged in, so an unrun hook is
never mistaken for a clean hook.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from ..config.loader import PipelineConfig
from ..errors import DataContractError
from ..io.fasta import assembly_stats
from ..logging_utils import get_logger
from ..manifest import SampleManifest
from ..models import AssemblyQC, RunMode

LOGGER = get_logger("stages.validation")

QC_COLUMNS: Tuple[str, ...] = (
    "sample_id",
    "assembly_size",
    "contig_count",
    "n50",
    "n90",
    "largest_contig",
    "gc_content",
    "ambiguous_bases",
    "basic_quality_status",
    "species_confirmation",
    "contamination_status",
    "completeness_status",
    "duplicate_status",
    "flag_reasons",
)

#: Statuses that mean "usable".
PASS_STATUSES = frozenset({"pass", "not_run", "not_assessable"})


@dataclass(frozen=True)
class QcVerdict:
    """Outcome of threshold evaluation for one sample."""

    status: str
    reasons: Tuple[str, ...]

    @property
    def passed(self) -> bool:
        return not self.reasons


def evaluate_thresholds(
    size: Optional[int],
    n50: Optional[int],
    ambiguous: Optional[int],
    config: PipelineConfig,
) -> QcVerdict:
    """Apply the configured QC thresholds.

    A metric that could not be measured produces a ``not_measurable:<metric>``
    reason rather than being treated as a pass or a fail.
    """
    qc = config.qc
    reasons: List[str] = []

    if size is None:
        reasons.append("not_measurable:assembly_size")
    else:
        if qc.min_assembly_size is not None and size < qc.min_assembly_size:
            reasons.append(f"below_min_assembly_size:{size}<{qc.min_assembly_size}")
        if qc.max_assembly_size is not None and size > qc.max_assembly_size:
            reasons.append(f"above_max_assembly_size:{size}>{qc.max_assembly_size}")

    if n50 is None:
        reasons.append("not_measurable:n50")
    elif qc.min_n50 is not None and n50 < qc.min_n50:
        reasons.append(f"below_min_n50:{n50}<{qc.min_n50}")

    if ambiguous is None:
        reasons.append("not_measurable:ambiguous_bases")
    elif qc.max_ambiguous_bases is not None and ambiguous > qc.max_ambiguous_bases:
        reasons.append(
            f"above_max_ambiguous_bases:{ambiguous}>{qc.max_ambiguous_bases}"
        )

    status = "pass" if not reasons else "flagged"
    return QcVerdict(status=status, reasons=tuple(reasons))


def _hook_status(config: PipelineConfig, hook: str) -> Optional[str]:
    """Return the recorded status of a validation hook.

    ``None`` means the hook has not been run. This is deliberately distinct
    from ``pass``: an unrun hook is an absence of information.
    """
    if not config.qc.hooks.get(hook, False):
        return None
    return "not_run"


def validate_sample(
    config: PipelineConfig, sample_id: str, fasta_path: Optional[Path]
) -> AssemblyQC:
    """Validate one assembly.

    A missing, empty or corrupt assembly file yields a QC record with status
    ``unreadable_input`` and a reason; it does **not** raise. One bad file
    must not abort a cohort-wide run, and it must not be silently counted as
    a pass either. Such samples are excluded from any downstream analysis
    that requires sequence.
    """
    if fasta_path is None:
        LOGGER.warning("No assembly path for %s", sample_id)
        return AssemblyQC(
            sample_id=sample_id,
            basic_quality_status="missing_input",
            flag_reasons="no_assembly_path",
        )

    path = Path(fasta_path)
    if not path.exists():
        LOGGER.warning("Assembly file does not exist for %s: %s", sample_id, path)
        return AssemblyQC(
            sample_id=sample_id,
            basic_quality_status="missing_input",
            flag_reasons="assembly_file_not_found",
        )

    try:
        summary = assembly_stats(path)
    except (DataContractError, OSError) as exc:
        # A truncated or garbled assembly: record why, and keep the cohort
        # running. This is a real failure mode of bulk-downloaded genome sets.
        reason = getattr(exc, "context", {}) or {}
        detail = str(exc)
        LOGGER.error("Cannot parse assembly for %s: %s", sample_id, detail)
        return AssemblyQC(
            sample_id=sample_id,
            basic_quality_status="unreadable_input",
            assembly_size=path.stat().st_size or None,
            flag_reasons=f"unreadable:{detail}",
        )

    verdict = evaluate_thresholds(
        summary.assembly_size, summary.n50, summary.ambiguous_bases, config
    )
    if not verdict.passed:
        LOGGER.debug("%s QC: %s", sample_id, "; ".join(verdict.reasons))

    return AssemblyQC(
        sample_id=sample_id,
        assembly_size=summary.assembly_size,
        contig_count=summary.contig_count,
        n50=summary.n50,
        n90=summary.n90,
        largest_contig=summary.largest_contig,
        gc_content=summary.gc_content,
        ambiguous_bases=summary.ambiguous_bases,
        basic_quality_status=verdict.status,
        species_confirmation=_hook_status(config, "species_confirmation"),
        contamination_status=_hook_status(config, "contamination"),
        completeness_status=_hook_status(config, "completeness"),
        duplicate_status=_hook_status(config, "duplicate_detection"),
        flag_reasons=";".join(verdict.reasons),
    )


def run(
    config: PipelineConfig, manifest: SampleManifest, mode: RunMode
) -> List[AssemblyQC]:
    """Run stage 1 across the manifest.

    Args:
        config: Loaded configuration.
        manifest: The sample manifest.
        mode: Run mode. Recorded for context only; stage 1 does not branch
            on mode because it reads whatever manifest it is given, and the
            manifest is already mode-scoped.
    """
    results: List[AssemblyQC] = []
    for sample in manifest:
        results.append(validate_sample(config, sample.sample_id, sample.assembly_path))

    by_status: Dict[str, int] = {}
    for record in results:
        by_status[record.basic_quality_status] = (
            by_status.get(record.basic_quality_status, 0) + 1
        )
    LOGGER.info(
        "Stage 1 [%s]: %d assemblies | status %s",
        mode.value,
        len(results),
        ", ".join(f"{k}={v}" for k, v in sorted(by_status.items())),
    )
    n_unreadable = sum(1 for r in results if r.basic_quality_status == "unreadable_input")
    if n_unreadable:
        LOGGER.error(
            "%d of %d assemblies are unreadable (empty or corrupt). They are "
            "excluded from any analysis requiring sequence and must be "
            "re-downloaded.",
            n_unreadable,
            len(results),
        )
    return results


def summarise(records: Sequence[AssemblyQC]) -> Dict[str, object]:
    """Cohort-level QC summary, used by the report."""
    sizes = [r.assembly_size for r in records if r.assembly_size is not None]
    n50s = [r.n50 for r in records if r.n50 is not None]
    gcs = [r.gc_content for r in records if r.gc_content is not None]
    unreadable = [
        r.sample_id for r in records if r.basic_quality_status == "unreadable_input"
    ]
    missing = [
        r.sample_id for r in records if r.basic_quality_status == "missing_input"
    ]
    flagged = [
        r.sample_id
        for r in records
        if r.basic_quality_status not in ("pass", "unreadable_input", "missing_input")
    ]
    return {
        "n_samples": len(records),
        "n_flagged": len(flagged),
        "flagged_sample_ids": flagged,
        "n_unreadable": len(unreadable),
        "unreadable_sample_ids": unreadable,
        "n_missing_input": len(missing),
        "missing_input_sample_ids": missing,
        "n_usable": len(records) - len(unreadable) - len(missing),
        "min_assembly_size": min(sizes) if sizes else None,
        "max_assembly_size": max(sizes) if sizes else None,
        "median_assembly_size": sorted(sizes)[len(sizes) // 2] if sizes else None,
        "min_n50": min(n50s) if n50s else None,
        "max_n50": max(n50s) if n50s else None,
        "min_gc_content": round(min(gcs), 3) if gcs else None,
        "max_gc_content": round(max(gcs), 3) if gcs else None,
    }


def rows(records: Sequence[AssemblyQC]) -> List[Dict[str, object]]:
    """Convert QC records to output rows in the declared column order."""
    return [{column: record.to_row().get(column) for column in QC_COLUMNS} for record in records]
