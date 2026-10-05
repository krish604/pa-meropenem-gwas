"""Stage 3 - MLST.

Parses MLST output into :class:`~papipeline.models.MlstCall`, keeping the
per-locus allele profile so that a sequence type can always be traced back
to the alleles that produced it.

A partial profile is reported as ``partial`` with ``ST=None``, never as a
guessed sequence type. ``no_call`` is reserved for "the tool ran and
returned nothing for this sample", which is a different situation from
"the profile is incomplete".
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from ..config.loader import PipelineConfig
from ..errors import PipelineError
from ..io.tsv import read_tsv
from ..logging_utils import get_logger
from ..manifest import SampleManifest
from ..models import MlstCall, MlstStatus, RunMode

LOGGER = get_logger("stages.mlst")

MLST_COLUMNS: Tuple[str, ...] = (
    "sample_id",
    "ST",
    "alleles",
    "MLST_status",
    "mlst_scheme",
    "allele_database",
)

REQUIRED = ("sample_id", "ST", "alleles", "MLST_status")

#: Sentinel values that mean "no sequence type", not a number.
_NULL_ST_VALUES = frozenset({"", ".", "NA", "N/A", "none", "None", "-", "0"})


def parse_alleles(field: Optional[str]) -> Dict[str, str]:
    """Parse an ``locus:allele`` allele string.

    Accepts ``;`` or ``,`` separated pairs. Malformed tokens are skipped
    with a warning rather than aborting the sample.
    """
    if not field:
        return {}
    alleles: Dict[str, str] = {}
    for token in field.replace(",", ";").split(";"):
        token = token.strip()
        if not token:
            continue
        if ":" not in token:
            LOGGER.warning("Skipping malformed allele token: %r", token)
            continue
        locus, allele = token.split(":", 1)
        locus, allele = locus.strip(), allele.strip()
        if not locus or not allele:
            LOGGER.warning("Skipping incomplete allele token: %r", token)
            continue
        if allele in _NULL_ST_VALUES:
            LOGGER.warning(
                "Allele for locus %s is a null sentinel; locus treated as unmatched",
                locus,
            )
            continue
        alleles[locus] = allele
    return alleles


def normalise_scheme(config: PipelineConfig, scheme: Optional[str]) -> str:
    """Resolve a scheme alias to its canonical name."""
    if not scheme:
        return str(config.raw.get("mlst", {}).get("fallback_scheme", "unknown"))
    aliases = dict(config.raw.get("mlst", {}).get("scheme_aliases") or {})
    return aliases.get(scheme, scheme)


def parse_row(
    row: Dict[str, Optional[str]], config: PipelineConfig
) -> MlstCall:
    """Convert one MLST table row into a typed call."""
    sample_id = str(row.get("sample_id") or "").strip()
    if not sample_id:
        raise ValueError("MLST row has no sample_id")

    raw_status = (row.get("MLST_status") or "").strip().lower()
    try:
        status = MlstStatus(raw_status) if raw_status else MlstStatus.NO_CALL
    except ValueError:
        LOGGER.warning(
            "Unrecognised MLST_status %r for %s; treating as no_call",
            raw_status,
            sample_id,
        )
        status = MlstStatus.NO_CALL

    alleles = parse_alleles(row.get("alleles"))
    raw_st = (row.get("ST") or "").strip()

    if raw_st in _NULL_ST_VALUES:
        sequence_type: Optional[str] = None
        if status is MlstStatus.TYPED and alleles:
            # A typed row with no usable ST is internally inconsistent.
            LOGGER.warning(
                "MLST row for %s claims 'typed' but has no usable ST; demoting to partial",
                sample_id,
            )
            status = MlstStatus.PARTIAL
    else:
        sequence_type = raw_st

    # Cross-check the reported status against the number of matched loci.
    min_loci = int(config.raw.get("mlst", {}).get("min_typing_locus_match", 5))
    if status is MlstStatus.TYPED and len(alleles) < min_loci:
        LOGGER.warning(
            "MLST for %s is typed on only %d loci (minimum %d); demoting to partial",
            sample_id,
            len(alleles),
            min_loci,
        )
        status = MlstStatus.PARTIAL
        sequence_type = None

    if not alleles and status is not MlstStatus.NO_CALL:
        status = MlstStatus.NO_CALL
        sequence_type = None

    return MlstCall(
        sample_id=sample_id,
        sequence_type=sequence_type,
        alleles=alleles,
        mlst_status=status.value,
        mlst_scheme=normalise_scheme(config, row.get("mlst_scheme")),
        allele_database=row.get("allele_database"),
    )


def load_mlst(
    config: PipelineConfig, path: Path
) -> Dict[str, MlstCall]:
    """Load an MLST table keyed by sample ID.

    A sample appearing twice is a hard error: MLST calls are not
    mergeable.
    """
    rows = read_tsv(path, required_columns=REQUIRED, unique_columns=("sample_id",))
    calls: Dict[str, MlstCall] = {}
    for row in rows:
        call = parse_row(row, config)
        calls[call.sample_id] = call

    by_status: Dict[str, int] = {}
    for call in calls.values():
        by_status[call.mlst_status] = by_status.get(call.mlst_status, 0) + 1
    LOGGER.info(
        "Stage 3: %d MLST calls | %s",
        len(calls),
        ", ".join(f"{k}={v}" for k, v in sorted(by_status.items())),
    )
    return calls


def align_to_manifest(
    calls: Dict[str, MlstCall], manifest: SampleManifest
) -> Dict[str, Optional[MlstCall]]:
    """Restrict MLST calls to manifest samples.

    A manifest sample with no MLST row maps to ``None`` so that "not typed"
    is distinguishable from a missing key.
    """
    return {sample_id: calls.get(sample_id) for sample_id in manifest.sample_ids}


def _allele_database(config: PipelineConfig) -> Optional[str]:
    """Name the scheme that produced an ST, or the stage cannot be audited.

    The config's `paeruginosa` is a legacy PubMLST spelling that normalises to
    `pseudomonas_aeruginosa`; both resolve to one scheme, and the value is
    recorded so a profile says which spelling was in force when it was called.
    """
    scheme = config.raw.get("mlst", {}).get("scheme")
    if not scheme:
        return None
    return f"pubmlst:{normalise_scheme(config, scheme)}"


def _run_by_calling_tool(
    config: PipelineConfig,
    manifest: SampleManifest,
    data_root: Path,
    workdir: Path,
) -> Dict[str, Optional[MlstCall]]:
    """REAL path: call `mlst` per isolate and apply the same contract.

    Every value the tool needs comes from ``config.raw['mlst']``; nothing is
    restated here. Rows are handed to :func:`parse_row`, so a REAL profile and a
    TEST fixture with the same content produce the same :class:`MlstCall` -
    otherwise the fixture would stop describing the behaviour.

    **A failed isolate is recorded, not raised.** 82% of the cohort's
    assemblies are corrupt (see :func:`papipeline.stages.variants`), so a stage
    that stopped at the first failure could never complete and one that
    ignored them would report a smaller cohort without saying so. Each failure
    is logged with its isolate and reason, and the isolate stays a cohort member
    carrying nothing, which :func:`align_to_manifest` represents as ``None``.
    """
    from ..adapters import mlst as mlst_tool
    from ..adapters.external import detect_tools, require_tool
    from ..assemblies import locate_assembly

    settings = dict(config.raw.get("mlst") or {})
    scheme = settings.get("scheme")
    if not scheme:
        raise PipelineError(
            "config.mlst.scheme is not set, so there is no scheme to call; the "
            "tool cannot auto-detect one for a cohort-wide comparison without "
            "the same scheme applied to every isolate"
        )

    # On demand, for the one tool this stage is about to invoke - not a sweep
    # over everything the project has heard of. Stage 3 runs mlst for every
    # isolate below, so its version is asked for by a real execution.
    statuses = detect_tools(["mlst"], execute=True)
    tool = require_tool(statuses, "mlst", "stage 3 (MLST)")
    tools = {"mlst": tool.executable}
    threads = int(config.runtime.get("threads", 1) or 1)

    calls: Dict[str, Optional[MlstCall]] = {}
    failures: List[str] = []

    for sample_id in manifest.sample_ids:
        try:
            assembly = locate_assembly(manifest.require(sample_id), data_root)
            result = mlst_tool.call_isolate(
                sample_id,
                assembly=assembly,
                scheme=str(scheme),
                workdir=Path(workdir) / sample_id,
                tools=tools,
                threads=threads,
            )
        except PipelineError as exc:
            LOGGER.warning(
                "Stage 3: mlst produced no call for %s, recorded as a cohort "
                "member with no sequence type (%s)", sample_id, exc,
            )
            calls[sample_id] = None
            failures.append(sample_id)
            continue

        calls[sample_id] = parse_row(
            result.as_stage_row(_allele_database(config)), config
        )

    # The aggregate the DAG declares as this stage's input.
    #
    # `load_mlst` reads `intermediate/mlst/mlst_results.tsv` and the Snakemake
    # rule declares it as MLST_IN, but the REAL branch returned records without
    # ever writing it - so the file had producers in TEST only (the synthetic
    # fixture generator and a legacy script) and none in REAL. The same gap
    # stage 4 had for `amr_determinants.tsv`, fixed in 69eb707; this is the
    # sixth instance of that class, and the second found in one day.
    _write_mlst_results(workdir / "mlst_results.tsv", calls)

    typed = sum(
        1 for c in calls.values() if c is not None and c.sequence_type
    )
    LOGGER.info(
        "Stage 3: %d/%d isolates typed by mlst (scheme %s)%s",
        typed, len(manifest.sample_ids), scheme,
        f"; {len(failures)} failed" if failures else "",
    )
    return calls


def _write_mlst_results(path: Path, calls: Dict[str, Any]) -> None:
    """Write the table `load_mlst` and the DAG both read.

    Header written even when every call is empty, for the same reason stage 4
    does it: "typed nothing" and "stage 3 did not run" must not be the same file
    state, because the second is a failure the DAG should surface.
    """
    from ..io.tsv import write_tsv

    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        call.to_row()
        for call in calls.values()
        if call is not None
    ]
    write_tsv(path, rows, list(MLST_COLUMNS))


def run(
    config: PipelineConfig,
    manifest: SampleManifest,
    mode: RunMode,
    intermediate_root: Path,
    data_root: Optional[Path] = None,
) -> Dict[str, Optional[MlstCall]]:
    """Stage 3 entry point.

    ``TEST`` reads the committed fixture table; ``REAL`` calls the tool. The
    two paths converge on :func:`parse_row`, so the status contract is applied
    identically either way.
    """
    if mode is RunMode.REAL:
        if data_root is None:
            raise PipelineError(
                "REAL mode needs a data_root to locate assemblies; refusing to "
                "fall back to anything, because substituting another sample's "
                "sequence would shift every profile silently"
            )
        return _run_by_calling_tool(
            config,
            manifest,
            data_root,
            Path(intermediate_root) / "mlst",
        )

    path = Path(intermediate_root) / "mlst" / "mlst_results.tsv"
    calls = load_mlst(config, path)
    aligned = align_to_manifest(calls, manifest)
    untyped = [sid for sid, call in aligned.items() if call is None]
    if untyped:
        LOGGER.warning(
            "%d manifest samples have no MLST record: %s",
            len(untyped),
            ",".join(untyped[:10]) + ("..." if len(untyped) > 10 else ""),
        )
    return aligned
