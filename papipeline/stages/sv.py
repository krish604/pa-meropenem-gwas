"""Stage 7 - structural variation and mobile genetic elements.

Accepts assembly-based SV results and normalises them, with the three-state
call status preserved end to end:

``confirmed``
    supported by the configured evidence requirements;
``candidate``
    suggestive only. Reported as ``candidate`` and mapped to
    :attr:`ClaimStatus.PREDICTED` downstream. Never promoted, regardless of
    configuration, unless a human explicitly edits the input table.
``not_assessable``
    the region could not be evaluated (contig gap, no coverage). Produces no
    mechanism call and is *not* counted as "no variant".

``config.structural_variants.promote_candidate_calls`` is deliberately
ignored. Promotion is a scientific decision requiring orthogonal evidence;
the pipeline does not make it.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from ..adapters.sv import SvPreflightError, screen_assembly
from ..config.loader import PipelineConfig
from ..io.tsv import read_tsv
from ..logging_utils import get_logger
from ..manifest import SampleManifest
from ..models import RunMode, StructuralCallStatus, StructuralVariant, SvType

LOGGER = get_logger("stages.sv")

SV_COLUMNS: Tuple[str, ...] = (
    "sample_id",
    "variant_id",
    "variant_type",
    "position",
    "affected_gene",
    "size",
    "evidence",
    "confidence",
    "call_status",
    "mge",
)

REQUIRED = ("sample_id", "variant_id", "variant_type", "call_status")


def parse_call_status(value: Optional[str]) -> StructuralCallStatus:
    """Parse a call status.

    An unrecognised status is *not* silently upgraded to ``confirmed``; it
    becomes ``candidate``, the conservative choice.
    """
    if value is None:
        return StructuralCallStatus.CANDIDATE
    token = str(value).strip().lower()
    try:
        return StructuralCallStatus(token)
    except ValueError:
        LOGGER.warning(
            "Unrecognised SV call_status %r; downgrading to 'candidate'", value
        )
        return StructuralCallStatus.CANDIDATE


def parse_sv_type(value: Optional[str]) -> SvType:
    if value is None:
        return SvType.UNKNOWN
    try:
        return SvType(str(value).strip().lower())
    except ValueError:
        LOGGER.warning("Unrecognised SV type %r; recording as unknown", value)
        return SvType.UNKNOWN


def _to_int(value: Optional[str]) -> Optional[int]:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def load_structural_variants(config: PipelineConfig, path: Path) -> List[StructuralVariant]:
    """Load and normalise the stage 7 table."""
    rows = read_tsv(path, required_columns=REQUIRED)
    min_size = config.raw.get("structural_variants", {}).get("minimum_deletion_size")

    results: List[StructuralVariant] = []
    for row in rows:
        status = parse_call_status(row.get("call_status"))
        sv_type = parse_sv_type(row.get("variant_type"))
        size = _to_int(row.get("size"))

        if (
            min_size is not None
            and sv_type is SvType.DELETION
            and size is not None
            and size < int(min_size)
        ):
            LOGGER.warning(
                "Deletion %s is smaller than the configured minimum (%s < %s); "
                "keeping the call but flagging it",
                row.get("variant_id"),
                size,
                min_size,
            )

        results.append(
            StructuralVariant(
                sample_id=str(row["sample_id"]),
                variant_id=str(row["variant_id"]),
                variant_type=sv_type.value,
                position=_to_int(row.get("position")),
                affected_gene=row.get("affected_gene"),
                size=size,
                evidence=str(row.get("evidence") or "unknown"),
                confidence=row.get("confidence"),
                call_status=status,
                mge=row.get("mge"),
            )
        )

    by_status: Dict[str, int] = {}
    for record in results:
        key = record.call_status.value
        by_status[key] = by_status.get(key, 0) + 1
    LOGGER.info(
        "Stage 7: %d SV records | by call status: %s",
        len(results),
        ", ".join(f"{k}={v}" for k, v in sorted(by_status.items())) or "none",
    )
    return results


def group_by_sample(
    records: Sequence[StructuralVariant], manifest: SampleManifest
) -> Dict[str, List[StructuralVariant]]:
    """Group SVs by sample, restricted to the manifest."""
    grouped: Dict[str, List[StructuralVariant]] = {sid: [] for sid in manifest.sample_ids}
    orphans: List[str] = []
    for record in records:
        if record.sample_id in grouped:
            grouped[record.sample_id].append(record)
        else:
            orphans.append(record.sample_id)
    if orphans:
        LOGGER.warning(
            "%d SV rows reference samples not in the manifest", len(orphans)
        )
    return grouped


def confirmed_only(
    records: Sequence[StructuralVariant],
) -> List[StructuralVariant]:
    """Filter to confirmed calls only.

    Used wherever a downstream stage needs a variant it can rely on. An
    empty result from this function does not mean "no SV"; it means "no
    confirmed SV".
    """
    return [r for r in records if r.call_status is StructuralCallStatus.CONFIRMED]


def mge_summary(records: Sequence[StructuralVariant]) -> Dict[str, int]:
    """Count MGEs associated with confirmed SVs."""
    counts: Dict[str, int] = {}
    for record in confirmed_only(records):
        if record.mge:
            counts[record.mge] = counts.get(record.mge, 0) + 1
    return counts


def run(
    config: PipelineConfig,
    manifest: SampleManifest,
    mode: RunMode,
    intermediate_root: Path,
) -> Dict[str, List[StructuralVariant]]:
    """Stage 7 entry point.

    REAL screens each prepared assembly with `nucmer` against the pinned
    reference and writes ``structural_variants/structural_variants.tsv``. TEST
    keeps reading the committed fixture.
    """
    if mode is not RunMode.TEST:
        return _run_real(config, manifest, intermediate_root)

    path = (
        Path(intermediate_root) / "structural_variants" / "structural_variants.tsv"
    )
    return group_by_sample(load_structural_variants(config, path), manifest)


def _run_real(
    config: PipelineConfig,
    manifest: SampleManifest,
    intermediate_root: Path,
) -> Dict[str, List[StructuralVariant]]:
    """Assembly-based calling, candidate-only.

    The refusal this replaces (a50bc0f) was correct while there was no caller.
    There is now one, and it is deliberately narrow: this pipeline's only input
    is assembled genomes, so there is no read depth, no coverage and no
    split-read or discordant-pair evidence. Every dedicated SV caller is built
    on that evidence and none of them applies - `sniffles`, `delly` and `gridss`
    all install on osx-arm64 and all require reads, `sista` (the assembly-based
    caller) is in no channel, and `mummer4` does not solve. What remains is an
    assembly-to-reference alignment, which supports a breakpoint with aligned
    flanks and nothing more.

    So every call is a **candidate**. `confirmed_only()` is empty for every REAL
    cohort, by construction. That is a property of the input, not a tuning
    problem, and it is why this does not quietly promote anything.
    """
    section = dict(config.raw.get("structural_variants") or {})
    if not section.get("caller"):
        raise SvPreflightError(
            "structural_variants.caller is not configured, so there is no "
            "program to run. Set it in config/science.yaml rather than "
            "defaulting here.",
        )

    reference = config.reference_fasta()
    if not reference.is_file():
        raise SvPreflightError(
            f"Pinned reference not found at {reference}. Stage 7 compares each "
            "assembly against it; without it there is nothing to compare to.",
            path=str(reference),
        )

    nucmer = _first_candidate(config, "nucmer")
    show_coords = _first_candidate(config, "show-coords")
    if not nucmer or not show_coords:
        raise SvPreflightError(
            "nucmer and show-coords are both required and at least one was not "
            "found. On osx-arm64 mummer comes from homebrew, not conda "
            "(`mummer4` does not solve: libcxx >=18 is absent), so the path is "
            "declared in the machine overlay's tool_search_dirs.",
            searched=list(config.tool_candidates("nucmer")) or "PATH only",
        )

    threads = config.threads or 1
    genome_root = config.assembly_root(RunMode.REAL)
    workdir = Path(intermediate_root) / "structural_variants" / "nucmer"
    min_sv_size = int(section.get("min_sv_size") or 1000)

    records: List[StructuralVariant] = []
    unassessable = 0
    screened = 0
    for sample in manifest:
        # An isolate with no assembly is recorded as screened-with-nothing
        # rather than dropped, or the denominator for an SV rate shrinks.
        if not sample.assembly_path:
            LOGGER.warning(
                "Stage 7: no assembly for %s, recorded as screened with no "
                "structural variant found", sample.sample_id,
            )
            continue
        screened += 1
        candidates, regions = screen_assembly(
            sample_id=sample.sample_id,
            genome=Path(sample.assembly_path),
            reference=reference,
            workdir=workdir,
            nucmer=nucmer,
            show_coords=show_coords,
            minmatch=int(section.get("nucmer_minmatch") or 20),
            mincluster=int(section.get("nucmer_mincluster") or 65),
            maxgap=int(section.get("nucmer_maxgap") or 90),
            maxmatch=bool(section.get("nucmer_maxmatch", True)),
            threads=threads,
            min_sv_size=min_sv_size,
            min_flanking_identity_pct=float(
                section.get("min_flanking_identity_pct") or 90.0
            ),
        )
        # Surfaced as a run warning rather than folded into a count of "no SV":
        # rule 7 forbids conflating an unevaluated region with an absent one.
        unassessable += len(regions)
        records.extend(_to_variant(c) for c in candidates)

    target = Path(intermediate_root) / "structural_variants"
    target.mkdir(parents=True, exist_ok=True)
    _write_sv_tsv(target / "structural_variants.tsv", records)
    LOGGER.info(
        "Stage 7: screened %d assemblies against %s | %d candidate calls | "
        "%d unassessable contig-gap regions (assembly-only input: no call can "
        "reach `confirmed`)",
        screened, reference.name, len(records), unassessable,
    )
    if unassessable:
        LOGGER.warning(
            "Stage 7: %d region(s) could not be evaluated - contig gaps or "
            "unaligned ends. These are NOT 'no variant present' and are not "
            "counted as such.", unassessable,
        )
    return group_by_sample(records, manifest)


def _first_candidate(config: PipelineConfig, tool: str) -> Optional[str]:
    for candidate in config.tool_candidates(tool):
        if Path(candidate).exists():
            return candidate
    return None


def _to_variant(candidate) -> StructuralVariant:
    return StructuralVariant(
        sample_id=candidate.sample_id,
        variant_id=candidate.variant_id,
        variant_type=candidate.variant_type,
        position=candidate.position,
        affected_gene=candidate.affected_gene,
        size=candidate.size,
        evidence=candidate.evidence,
        confidence=candidate.confidence,
        call_status=StructuralCallStatus(candidate.call_status),
        mge=candidate.mge,
    )


def _write_sv_tsv(path: Path, records: Sequence[StructuralVariant]) -> None:
    import csv

    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(SV_COLUMNS), delimiter="\t"
        )
        writer.writeheader()
        for record in records:
            writer.writerow(record.to_row())
