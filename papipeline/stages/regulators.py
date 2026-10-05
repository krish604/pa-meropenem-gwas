"""Stage 6 - regulator / chromosomal mutation screen.

Consumes variant calls for the loci in ``config/regulators.tsv`` and
normalises them into :class:`~papipeline.models.RegulatorVariant` records.

**Who writes the table.** Until this round nothing did: ``run()`` read
``<intermediate_root>/regulators/regulator_variants.tsv`` and returned, and the
only code in the repository that ever wrote that filename was the TEST fixture
emulator, ``papipeline/testing/synthetic.py``. In REAL the DAG therefore had an
undeclared input, which is what ``workflow/Snakefile`` and
``papipeline/stages/cooccurrence.py`` were both documenting rather than fixing.

The producer is :func:`_run_real`, and it is deliberately thin. The coordinate
intersection and the variant classification are
:mod:`papipeline.adapters.gff`'s - ``load_gene_intervals``, ``assign``,
``codon_at``, ``classify`` - which were written, documented and unit-tested for
exactly this purpose and had **no caller anywhere in the repository**. This
module supplies the two things that were missing, a caller and a writer; it does
not reimplement the screen.

Its shape follows ``stages.sv._run_real``, the same pattern for the same
reason: refuse by name when an input is missing, derive the table in REAL, and
keep a fixture read for TEST only.

**Why the calls arrive as an argument.** Stage 6's calls live in the caller's
scope (``papipeline/run.py``, the ``variants`` branch) and are re-read from disk
by nobody, precisely so the cohort merge cannot disagree with them. So the
screen takes them as ``calls_by_isolate`` and **refuses by name** when they are
absent, rather than writing an empty table that would report the whole cohort as
carrying no regulator variant.

**Where the screen's own input comes from.** ``docs/data_contract.md:114``
declares ``regulators/regulator_variants.tsv`` a *stage input table*, and its
section says that in REAL the tools write those into the run's intermediate
directory. That is the table :func:`run` reads, in every mode: the producer
(:func:`produce_regulator_variants`, reached through ``papipeline/run.py``)
writes it, and the stage reads it back off disk. What arrives as an *argument*
is the calls the producer derives that table from, and passing
``calls_by_isolate=`` to :func:`run` is now a deprecated door onto the producer
rather than this stage's input.

**Empty is not the same as not-run.** A screen that runs and finds nothing
writes a zero-row table *and* a ``regulator_screen.json`` sidecar recording how
many isolates and loci it examined. A table with no sidecar is the ambiguous
case, and :func:`load_regulator_variants` refuses it rather than reading it as a
negative. See :func:`read_screen_report`.

Three safeguards, all pre-existing:

* a variant in a gene that is not on the regulator list is still reported,
  but tagged so the report can distinguish screened loci from incidental
  findings;
* promoter-region calls are refused unless
  ``config.regulators.promoter_assessment_enabled`` is true, because a
  reliable promoter call needs a pinned reference and an aligner, neither
  of which exists yet (scientific rule 7).

**Consequences come from coordinates, not from ``bcftools csq``.** ``csq`` is
installed (1.23.1) and is nevertheless unusable on this project's pinned
reference: the annotation is a RefSeq *genomic* GFF3 whose 5,573 CDS features
are parented to ``gene`` nodes with no ``mRNA`` transcript anywhere, so
``bcftools csq`` reports ``Indexed 13 transcripts, 13 exons, 0 CDSs, 0 UTRs`` and
annotates nothing while exiting 0. Consequences are therefore derived from the
CDS coordinates and the reference bases directly. The vocabulary is
:class:`~papipeline.models.VariantType`, which is what the reader parses; the
contract declares the ``variant_type`` and ``effect`` columns but no controlled
vocabulary for either, and this module invents none.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import (
    Any,
    Dict,
    Iterable,
    List,
    Mapping,
    Optional,
    Sequence,
    Set,
    Tuple,
)

from ..adapters import gff
from ..config.loader import PipelineConfig, RegulatorSpec
from ..errors import DataContractError
from ..io.tsv import read_tsv, write_tsv
from ..logging_utils import get_logger
from ..manifest import SampleManifest
from ..models import (
    DISRUPTIVE_VARIANT_TYPES,
    ClaimStatus,
    RegulatorVariant,
    RunMode,
    VariantType,
)

LOGGER = get_logger("stages.regulators")

REGULATOR_COLUMNS: Tuple[str, ...] = (
    "sample_id",
    "gene",
    "variant",
    "variant_type",
    "position",
    "reference",
    "alternate",
    "effect",
    "mechanism",
    "confidence",
    "evidence_source",
    "call_status",
)

REQUIRED = ("sample_id", "gene", "variant", "variant_type")

#: Variant classes that are treated as disruptive to a gene's function.
#:
#: Aliased from :data:`papipeline.models.DISRUPTIVE_VARIANT_TYPES` so this stays
#: the one definition consulted by :func:`is_disruptive`, by :func:`gene_status`
#: and by the synthetic emulator alike. It used to be spelled out here while
#: :func:`gene_status` ignored it and matched two type names by literal string,
#: which is how a real FRAMESHIFT came to be reported as an ordinary "variant".
DISRUPTIVE_TYPES = DISRUPTIVE_VARIANT_TYPES

# ---------------------------------------------------------------------------
# Per-locus coverage.
#
# The distinction this exists to protect: "we could not see it" must never be
# reported as "there is none". Stage 6 screens a locus by intersecting calls
# with it, and an intersection is silent about loci the aligner never reached -
# an unaligned locus and a clean locus produce identical output, because in both
# cases no call is emitted. On this project's own real data that is not
# hypothetical: 8 of the 10 smoke isolates have NO primary alignment over any
# screened locus and reach all ten on SUPPLEMENTARY records alone (measured:
# 80 of 100 (isolate, locus) pairs), and a supplementary alignment can start
# or stop in the middle of a CDS. Measured on the 10-isolate smoke set, one pair
# already sits at 0.9641 (PDT000034122.1 x ampR: the alignment begins at
# 4593022, 32 nt into a CDS that starts at 4592990), and the missing 32 nt are
# silent.
#
# The coverage figure counts ANY alignment over the locus - primary, secondary
# and supplementary alike - because they are all evidence about the locus's
# sequence. Excluding supplementary records, as some tools default to, would
# report most of this cohort's loci as unseen.
# ---------------------------------------------------------------------------

#: Minimum fraction of a locus's reference positions that must be covered by an
#: alignment before this stage will say anything about that locus.
#:
#: The default is 0.90 and the reason it is not higher is that the threshold
#: must not fire on data the screen can actually read. Every one of the 100
#: (isolate, locus) pairs measured on the 10-isolate smoke set is at 1.0000
#: except one at 0.9641, so 0.90 leaves the whole cohort assessed and the
#: threshold only fires on a locus that is genuinely unaligned.
#:
#: Its measured sensitivity, and it is not free: with the test written as
#: ``coverage >= MIN_LOCUS_COVERAGE_FRACTION``, any threshold in (0, 0.9641]
#: changes zero classifications on this cohort; any threshold in (0.9641, 1.00]
#: changes exactly one - PDT000034122.1 x ampR, and only because 32 nt at the
#: start of that CDS are unaligned. One observation therefore carries the whole
#: threshold's behaviour, and it is named here so that tightening the constant
#: is a decision about that observation rather than an unexplained change in
#: output. A threshold of 1.00 would demand that no aligner boundary ever fall
#: inside a CDS, which on draft assemblies is not a property of the biology.
MIN_LOCUS_COVERAGE_FRACTION: float = 0.90

#: The one reason this module gives for not assessing a locus.
#:
#: A NAMED reason, not a sentence and not a boolean. ``not_assessed`` with this
#: reason is a positive statement - "the locus was not sufficiently covered" -
#: and it is deliberately distinct from any statement about the sequence at the
#: locus. Nothing else in this module produces this value, so a reader can
#: distinguish "we did not look" from every other outcome by string equality.
COVERAGE_NOT_ASSESSED_REASON: str = "locus_coverage_below_minimum"

#: The status label for a locus the screen declines to speak about. Paired with
#: :data:`COVERAGE_NOT_ASSESSED_REASON` and never used on its own to mean
#: "clean": :func:`gene_status` reaches it only through
#: :func:`not_assessed_reason`.
NOT_ASSESSED: str = "not_assessed"


def covered_reference_positions(
    alignments: Iterable[Any],
    start: int,
    end: int,
    *,
    include_supplementary: bool = True,
    include_secondary: bool = True,
) -> Set[int]:
    """Reference positions of ``[start, end]`` that any alignment covers.

    ``start``/``end`` are 1-based inclusive, the frame every interval in this
    repository uses. Positions come back in the same frame, so a caller can
    compare the count against ``end - start + 1`` without an off-by-one.

    **Supplementary and secondary alignments count.** They are evidence about
    the locus's sequence: minimap2 splits a draft assembly into one primary
    alignment plus dozens of supplementary segments, so on this project's own
    real data 8 of 10 isolates reach oprD only on a supplementary record.
    Excluding them - which is what several tools do by default - would report
    most of this cohort's loci as unseen and turn the coverage check into a
    constant false alarm. The flags are parameters rather than a fixed choice
    only so the exclusion can be demonstrated rather than asserted; the default
    is the one this stage uses.

    Args:
        alignments: Anything iterable yielding objects exposing
            ``reference_start`` (0-based), ``reference_end`` (exclusive),
            ``is_supplementary`` and ``is_secondary``. A ``pysam``
            ``AlignmentFile`` over one contig satisfies this, and so does a list
            of four-field stand-ins, which is what the tests use.
        start: Locus start, 1-based inclusive.
        end: Locus end, 1-based inclusive.
        include_supplementary: Count records flagged supplementary.
        include_secondary: Count records flagged secondary.

    Returns:
        The set of covered positions. Empty when nothing overlaps.
    """
    covered: Set[int] = set()
    for record in alignments:
        if not include_supplementary and getattr(
            record, "is_supplementary", False
        ):
            continue
        if not include_secondary and getattr(record, "is_secondary", False):
            continue
        lo = max(int(record.reference_start) + 1, start)
        hi = min(int(record.reference_end), end)
        if hi >= lo:
            covered.update(range(lo, hi + 1))
    return covered


def locus_coverage_fraction(
    alignments: Iterable[Any],
    start: int,
    end: int,
    *,
    include_supplementary: bool = True,
    include_secondary: bool = True,
) -> float:
    """:func:`covered_reference_positions` as a fraction of the locus length.

    ``end < start`` - an interval given the wrong way round - yields ``0.0``
    rather than raising: a malformed interval means nothing was covered, and a
    zero is the honest reading of that. It is not the same as a locus that was
    measured and found empty, and the threshold treats them the same, so the
    distinction is recorded here rather than left to the reader.
    """
    if end < start:
        return 0.0
    return len(covered_reference_positions(
        alignments, start, end,
        include_supplementary=include_supplementary,
        include_secondary=include_secondary,
    )) / float(end - start + 1)


def not_assessed_reason(
    coverage: Optional[float],
    min_locus_coverage: float = MIN_LOCUS_COVERAGE_FRACTION,
) -> Optional[str]:
    """Why this locus is not assessed, or ``None`` when it is.

    Returns :data:`COVERAGE_NOT_ASSESSED_REASON` when the locus was covered
    below ``min_locus_coverage``, and ``None`` otherwise - including for
    ``coverage=None``, which means the caller supplied no coverage evidence at
    all. That is deliberate and it is the safe direction: with no coverage
    figure there is nothing to compare against, and inventing a passing value
    would assert an observation nobody made. Absence of coverage evidence is
    reported by :func:`locus_coverage_rows` as ``coverage=None``, not by
    silently passing.

    The comparison is ``coverage >= min_locus_coverage``, so a locus exactly at
    the threshold is assessed. A threshold is a floor on what may be claimed,
    and a value that meets it meets it; making the boundary exclusive would mean
    the constant's own documented default had to be written as 0.9000001 to
    include a fully-covered locus.

    Args:
        coverage: Fraction of the locus's reference positions covered by an
            alignment, or ``None`` when unknown.
        min_locus_coverage: The threshold. Defaults to
            :data:`MIN_LOCUS_COVERAGE_FRACTION`.

    Returns:
        The named reason string, or ``None``.
    """
    if coverage is None:
        return None
    return (
        None
        if coverage >= min_locus_coverage
        else COVERAGE_NOT_ASSESSED_REASON
    )


def locus_assessment(
    gene: str,
    records: Sequence[RegulatorVariant],
    *,
    coverage: Optional[float] = None,
    min_locus_coverage: float = MIN_LOCUS_COVERAGE_FRACTION,
) -> Dict[str, Any]:
    """One locus's verdict with the evidence that produced it.

    Combines what the calls say with whether the calls could have said anything.
    A locus whose coverage is below the threshold reports status
    :data:`NOT_ASSESSED` and reason
    :data:`COVERAGE_NOT_ASSESSED_REASON` **regardless of the records**, because
    a call in a locus the aligner did not fully read is not evidence that the
    locus is clean - it is evidence that part of the locus was never examined.

    Args:
        gene: The locus, as named in ``config/regulators.tsv``.
        records: This isolate's records at this locus.
        coverage: Covered fraction, or ``None`` when no coverage evidence was
            supplied.
        min_locus_coverage: Threshold; defaults to
            :data:`MIN_LOCUS_COVERAGE_FRACTION`.

    Returns:
        ``{"gene", "status", "reason", "coverage", "min_locus_coverage",
        "records"}``. ``reason`` is ``None`` whenever ``status`` is not
        :data:`NOT_ASSESSED` for want of coverage.
    """
    reason = not_assessed_reason(coverage, min_locus_coverage)
    status = gene_status(records).get(gene, NOT_ASSESSED)
    if reason is not None:
        status = NOT_ASSESSED
    return {
        "gene": gene,
        "status": status,
        "reason": reason,
        "coverage": coverage,
        "min_locus_coverage": min_locus_coverage,
        "records": len(records),
    }


def locus_coverage_rows(
    sample_ids: Sequence[str],
    genes: Sequence[str],
    coverage_by_isolate_gene: Optional[Mapping[Tuple[str, str], float]],
    *,
    min_locus_coverage: float = MIN_LOCUS_COVERAGE_FRACTION,
) -> List[Dict[str, Any]]:
    """The isolate x locus coverage table, one row per pair.

    This is the artefact that makes a below-threshold verdict checkable: a
    reader can see the covered fraction, the threshold it was compared against,
    and the verdict, for every pair the screen looked at - including the pairs
    it declined to assess, which is precisely the set a bare result table would
    have made invisible.

    Args:
        sample_ids: Cohort members, in the order the caller wants them.
        genes: Loci, likewise.
        coverage_by_isolate_gene: ``(sample_id, gene) -> covered fraction``, or
            ``None`` when no coverage evidence exists. ``None`` yields rows with
            ``coverage=None`` and ``reason=None``; it never yields a pass.
        min_locus_coverage: Threshold; defaults to
            :data:`MIN_LOCUS_COVERAGE_FRACTION`.

    Returns:
        One dict per ``(sample_id, gene)``, sample-major then gene, with keys
        ``sample_id``, ``gene``, ``coverage``, ``min_locus_coverage``,
        ``not_assessed`` (bool) and ``reason``.
    """
    rows: List[Dict[str, Any]] = []
    for sample_id in sample_ids:
        for gene in genes:
            coverage = (
                None
                if coverage_by_isolate_gene is None
                else coverage_by_isolate_gene.get((sample_id, gene))
            )
            reason = not_assessed_reason(coverage, min_locus_coverage)
            rows.append({
                "sample_id": sample_id,
                "gene": gene,
                "coverage": coverage,
                "min_locus_coverage": min_locus_coverage,
                "not_assessed": reason is not None,
                "reason": reason,
            })
    return rows


def parse_variant_type(value: Optional[str]) -> VariantType:
    """Parse a variant class, defaulting to ``UNKNOWN`` rather than raising."""
    if value is None:
        return VariantType.UNKNOWN
    try:
        return VariantType(str(value).strip().upper())
    except ValueError:
        LOGGER.warning("Unrecognised variant_type %r; recording as UNKNOWN", value)
        return VariantType.UNKNOWN


def _to_int(value: Optional[str]) -> Optional[int]:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def parse_claim(value: Optional[str]) -> ClaimStatus:
    if value is None:
        return ClaimStatus.DETECTED
    try:
        return ClaimStatus(str(value).strip().upper())
    except ValueError:
        LOGGER.warning("Unrecognised call_status %r; recording as DETECTED", value)
        return ClaimStatus.DETECTED


def load_regulator_variants(
    config: PipelineConfig, path: Path, intermediate_root: Optional[Path] = None
) -> List[RegulatorVariant]:
    """Load and normalise the stage 6 variant table.

    Args:
        path: The table itself.
        intermediate_root: The run's intermediate directory, used only to find
            the provenance sidecar. Optional, so existing callers are unaffected;
            when given, a zero-row table is reported as *"the screen ran and found
            nothing"* rather than as the generic parse failure below.

    Raises:
        DataContractError: The table is missing, malformed, or empty. The empty
            case is split in two, because those are different findings and
            collapsing them is how a screen that never ran comes to read as a
            cohort carrying no regulator variant:

            * sidecar present -> the screen **ran**, over the recorded number of
              isolates and loci, and found nothing. A result.
            * sidecar absent -> the file was **not produced by this screen**, and
              its emptiness is evidence of nothing at all.
    """
    try:
        rows = read_tsv(path, required_columns=REQUIRED)
    except DataContractError as exc:
        if "no data rows" not in str(exc):
            raise
        report = read_screen_report(intermediate_root) if intermediate_root else None
        if report is None:
            raise DataContractError(
                f"{path} holds a header and no rows, and carries no "
                f"{SCREEN_REPORT_NAME} beside it. An empty table with no "
                "provenance is ambiguous between 'the screen ran and found "
                "nothing' and 'the screen never ran', and the second must not be "
                "read as a cohort carrying no regulator variant. In REAL the "
                "table is written by this stage's own producer, which always "
                "writes the sidecar; if you are reading one produced elsewhere, "
                "say which producer wrote it.",
                path=str(path),
                reason="empty_table_without_provenance",
            ) from exc
        raise DataContractError(
            f"{path} holds a header and no rows. That is a result, not a "
            f"failure: the screen ran over {report.get('isolates_screened')} of "
            f"{report.get('isolates_in_manifest')} isolates against "
            f"{report.get('loci_screened')} loci and found no variant in any of "
            f"them. See {SCREEN_REPORT_NAME} for the full accounting, including "
            f"{report.get('calls_outside_screened_loci')} call(s) that fell "
            "outside the screened loci.",
            path=str(path),
            reason="screen_ran_and_found_nothing",
            isolates_screened=report.get("isolates_screened"),
            loci_screened=report.get("loci_screened"),
        ) from exc

    promoter_enabled = bool(
        config.raw.get("regulators", {}).get("promoter_assessment_enabled", False)
    )

    results: List[RegulatorVariant] = []
    for row in rows:
        sample_id = str(row["sample_id"])
        gene = str(row["gene"])
        variant_type = parse_variant_type(row.get("variant_type"))

        if variant_type is VariantType.PROMOTER_ALTERATION and not promoter_enabled:
            LOGGER.warning(
                "Promoter alteration for %s/%s dropped: promoter assessment is "
                "disabled in configuration",
                sample_id,
                gene,
            )
            continue

        spec = config.regulator(gene)
        if spec is None:
            LOGGER.warning(
                "Gene %s is not on the regulator screen list; recording the "
                "call but flagging it as incidental",
                gene,
            )
        elif variant_type.value not in spec.variant_classes:
            LOGGER.warning(
                "Variant type %s is outside the screened classes for %s",
                variant_type.value,
                gene,
            )

        mechanism = row.get("mechanism")
        if mechanism is None and spec is not None:
            mechanism = spec.mechanism

        results.append(
            RegulatorVariant(
                sample_id=sample_id,
                gene=gene,
                variant=str(row["variant"]),
                variant_type=variant_type.value,
                position=_to_int(row.get("position")),
                reference=row.get("reference"),
                alternate=row.get("alternate"),
                effect=row.get("effect"),
                mechanism=mechanism,
                confidence=row.get("confidence"),
                evidence_source=str(row.get("evidence_source") or "unknown"),
                call_status=parse_claim(row.get("call_status")),
            )
        )

    LOGGER.info("Stage 6: %d variant calls loaded", len(results))
    return results


def group_by_sample(
    records: Sequence[RegulatorVariant], manifest: SampleManifest
) -> Dict[str, List[RegulatorVariant]]:
    """Group variants by sample, restricted to the manifest."""
    grouped: Dict[str, List[RegulatorVariant]] = {sid: [] for sid in manifest.sample_ids}
    orphans: List[str] = []
    for record in records:
        if record.sample_id in grouped:
            grouped[record.sample_id].append(record)
        else:
            orphans.append(record.sample_id)
    if orphans:
        LOGGER.warning(
            "%d variant rows reference samples not in the manifest",
            len(orphans),
        )
    return grouped


def gene_status(
    records: Sequence[RegulatorVariant],
    *,
    coverage: Optional[Mapping[str, float]] = None,
    min_locus_coverage: float = MIN_LOCUS_COVERAGE_FRACTION,
) -> Dict[str, str]:
    """Per-gene status for one sample: intact / variant / absent / disrupted.

    ``absent`` and ``disrupted`` are kept separate from ``variant`` because
    only the former two indicate loss of the locus.

    What counts as ``disrupted`` is decided by :func:`is_disruptive`, not by a
    list of type names written out here. That is deliberate: the branch used to
    special-case GENE_ABSENCE and GENE_DISRUPTION only, so a FRAMESHIFT or a
    PREMATURE_STOP - both of which truncate the coding sequence, and both of
    which are in :data:`DISRUPTIVE_TYPES` - fell through to ``variant``, and
    ``oprD_LoF`` read 0 on real data for a genuinely lost oprD. GENE_ABSENCE is
    still tested first, because absence and disruption are different
    observations and must not collapse into one label.

    ``coverage`` is optional and additive. When supplied - ``gene -> covered
    fraction`` for THIS sample - any gene below ``min_locus_coverage`` is
    reported as :data:`NOT_ASSESSED`, overriding whatever its records imply.
    The override is what keeps "we could not see it" from being reported as
    "there is none"; without it a partially-aligned locus and a clean locus
    are indistinguishable in this return value. Genes absent from the mapping
    are assessed normally, because a locus nobody supplied coverage for is not
    a locus that was measured and found unseen.
    """
    status: Dict[str, str] = {}
    for record in records:
        current = status.get(record.gene)
        rank = {"intact": 0, "variant": 1, "disrupted": 2, "absent": 3}
        if record.variant_type == VariantType.GENE_ABSENCE.value:
            label = "absent"
        elif is_disruptive(record):
            label = "disrupted"
        elif current in ("absent", "disrupted"):
            continue
        else:
            label = "variant"
        if current is None or rank[label] > rank.get(current, 0):
            status[record.gene] = label
    if coverage:
        for gene, fraction in coverage.items():
            if not_assessed_reason(fraction, min_locus_coverage) is not None:
                status[gene] = NOT_ASSESSED
    return status


def is_disruptive(record: RegulatorVariant) -> bool:
    """Whether a call is expected to impair the gene's function.

    The single definition of "disruptive" in the pipeline. :func:`gene_status`
    routes its ``disrupted`` verdict through here so the two cannot disagree.
    """
    return record.variant_type in DISRUPTIVE_TYPES


def oprd_status(
    records: Sequence[RegulatorVariant],
    *,
    coverage: Optional[float] = None,
    min_locus_coverage: float = MIN_LOCUS_COVERAGE_FRACTION,
) -> str:
    """OprD status for one sample.

    Returns ``absent``, ``disrupted``, ``variant`` or ``not_assessed``. The
    last value is distinct from ``intact``: absence of a call is not
    evidence of an intact gene until stage 2 confirms the locus is present.

    ``coverage`` is the covered fraction of the oprD CDS. When it is supplied
    and falls below ``min_locus_coverage`` the answer is ``not_assessed``
    whatever the records say, and :func:`locus_assessment` is what pairs that
    verdict with its named reason.
    """
    status = gene_status(
        records,
        coverage=None if coverage is None else {"oprD": coverage},
        min_locus_coverage=min_locus_coverage,
    ).get("oprD")
    return status if status in ("absent", "disrupted", "variant") else NOT_ASSESSED


#: Feature labels for the two OprD observations.
#:
#: They are kept distinct end to end and neither is derived from the other: a
#: locus that is wholly absent and a locus that is present but carries a
#: truncating variant are different observations with different mechanisms, and
#: collapsing them destroys exactly the distinction this pair exists to carry.
OPRD_ABSENT_FEATURE = "oprD_absent"
OPRD_LOF_FEATURE = "oprD_LoF"


def oprd_feature_labels() -> Tuple[str, ...]:
    """The two OprD feature labels, absent first, then loss-of-function."""
    return (OPRD_ABSENT_FEATURE, OPRD_LOF_FEATURE)


def oprd_features(
    records: Sequence[RegulatorVariant],
    *,
    coverage: Optional[float] = None,
    min_locus_coverage: float = MIN_LOCUS_COVERAGE_FRACTION,
) -> Dict[str, int]:
    """The two OprD features for one sample, as binary values.

    Both features are read off the single :func:`oprd_status` verdict rather
    than being computed independently, so the two can never disagree and
    neither can be derived from the other.

    ``not_assessed`` and ``variant`` both yield ``0`` for each feature. That is
    the deliberate direction: this stage can only report a locus it actually
    observed a call in, so silence is recorded as "not observed" rather than
    guessed either way. Absence of a call is not evidence of an intact gene
    until stage 2 confirms the locus is present.

    A below-threshold ``coverage`` yields ``not_assessed`` for the same reason,
    and additionally suppresses any record-derived verdict, so a feature can
    never read 1 on the strength of evidence from a locus the aligner did not
    read. That case is not detectable from the features alone; the caller must
    read :func:`locus_coverage_rows` for the covered fraction and its named
    reason.
    """
    status = oprd_status(
        records, coverage=coverage, min_locus_coverage=min_locus_coverage
    )
    return {
        OPRD_ABSENT_FEATURE: int(status == "absent"),
        OPRD_LOF_FEATURE: int(status == "disrupted"),
    }


def oprd_feature_columns(feature_type: str) -> Tuple[str, ...]:
    """Column names for the two OprD features under ``feature_type``.

    The stage 12 contract names a feature ``<type>__<label>``, so the declared
    prefix has to be supplied by the caller rather than hard-coded here: which
    prefix an absent/LoF pair should carry is a schema decision, not a fact
    this stage can settle on its own.
    """
    return tuple(f"{feature_type}__{label}" for label in oprd_feature_labels())


def oprd_feature_rows(
    grouped: Mapping[str, Sequence[RegulatorVariant]],
    feature_type: str,
    *,
    coverage: Optional[Mapping[Tuple[str, str], float]] = None,
    min_locus_coverage: float = MIN_LOCUS_COVERAGE_FRACTION,
) -> List[Dict[str, object]]:
    """One feature row per sample, columns in a stable order.

    Mirrors the shape of the stage 12 feature table: a ``sample_id`` column
    plus one binary column per feature, ordered sample-major so the output does
    not depend on dict iteration order.

    ``coverage`` is the optional ``(sample_id, "oprD") -> covered fraction``
    input. Supplied, a sample whose oprD coverage is below
    ``min_locus_coverage`` gets both features 0 whatever its records say; the
    reason is not in this table, because the table's columns are fixed by the
    stage-12 contract, and it is in :func:`locus_coverage_rows`.
    """
    labels = oprd_feature_labels()
    columns = oprd_feature_columns(feature_type)
    rows: List[Dict[str, object]] = []
    for sample_id in sorted(grouped):
        features = oprd_features(
            grouped[sample_id],
            coverage=None if coverage is None else coverage.get((sample_id, "oprD")),
            min_locus_coverage=min_locus_coverage,
        )
        row: Dict[str, object] = {"sample_id": sample_id}
        for label, column in zip(labels, columns):
            row[column] = features[label]
        rows.append(row)
    return rows


#: The stage's REAL input table, relative to the root ``run`` is handed.
#:
#: Named in ``docs/data_contract.md:114`` under "Stage input tables (under
#: ``intermediate/``)". One source, so the reader's path and the producer's path
#: cannot be spelled two ways and disagree.
REGULATOR_TABLE_RELPATH = ("regulators", "regulator_variants.tsv")


def regulator_table_path(root: Path) -> Path:
    """Where this stage's table lives under ``root``.

    Args:
        root: The stage-input root (``config.tool_output_root(mode)``).

    Returns:
        ``<root>/regulators/regulator_variants.tsv``.
    """
    return Path(root).joinpath(*REGULATOR_TABLE_RELPATH)


def produce_regulator_variants(
    config: PipelineConfig,
    manifest: SampleManifest,
    intermediate_root: Path,
    *,
    calls_by_isolate: Mapping[str, Sequence[Mapping[str, str]]],
    reference_fasta: Optional[Path] = None,
    locus_coverage: Optional[Mapping[Tuple[str, str], float]] = None,
    min_locus_coverage: float = MIN_LOCUS_COVERAGE_FRACTION,
) -> Dict[str, List[RegulatorVariant]]:
    """Derive ``regulators/regulator_variants.tsv`` from this run's calls.

    The producer half of stage 6, split out from :func:`run` so that the stage's
    entry point can be a plain **reader** of the table while the wiring that
    feeds it lives at the call site (:mod:`papipeline.run`). Splitting them is
    what lets ``run.py`` own the on-disk read of ``variants.tsv`` without
    ``run()`` having to know where the calls came from.

    ``locus_coverage`` and ``min_locus_coverage`` carry the same meaning as on
    :func:`run` and are applied here, at the point where the calls and the
    coverage evidence meet. This is the only half that can suppress records
    before they reach the table.
    """
    return _run_real(
        config,
        manifest,
        intermediate_root,
        calls_by_isolate=calls_by_isolate,
        reference_fasta=reference_fasta,
        locus_coverage=locus_coverage,
        min_locus_coverage=min_locus_coverage,
    )


def run(
    config: PipelineConfig,
    manifest: SampleManifest,
    mode: RunMode,
    intermediate_root: Path,
    *,
    calls_by_isolate: Optional[Mapping[str, Sequence[Mapping[str, str]]]] = None,
    reference_fasta: Optional[Path] = None,
    locus_coverage: Optional[Mapping[Tuple[str, str], float]] = None,
    min_locus_coverage: float = MIN_LOCUS_COVERAGE_FRACTION,
) -> Dict[str, List[RegulatorVariant]]:
    """Stage 6 entry point. A **reader** of the stage's table, in every mode.

    ``TEST`` reads the committed fixture. Every other mode reads
    ``regulators/regulator_variants.tsv`` — the stage input table named at
    ``docs/data_contract.md:114`` — and refuses by name when it is absent.
    ``run.py`` produces that table from stage 6's own per-isolate calls and then
    calls this, so the value that reaches the rest of the pipeline is the one on
    disk rather than a value still held in a caller's scope.

    Args:
        calls_by_isolate: **Deprecated.** This run's stage-6 calls,
            ``sample_id -> rows`` in
            :data:`~papipeline.stages.variants.PER_ISOLATE_COLUMNS`. Retained
            only so a caller still holding the calls is served by the producer
            rather than silently ignored; the parameter is no longer this
            stage's input. ``papipeline/run.py`` does not pass it.
        reference_fasta: Override for the pinned reference. Defaults to
            ``config.reference_fasta()``; read from configuration rather than
            taken as a path, per the repo's no-hard-coded-paths rule. Used by
            the producer only.
        locus_coverage: Optional ``(sample_id, gene) -> covered fraction``.
            This stage cannot derive it - it sees calls, not alignments - so it
            is an input rather than a computation. Absent, the screen behaves
            exactly as before and the sidecar records
            ``coverage_evidence: absent``; supplied, a pair below
            ``min_locus_coverage`` is reported ``not_assessed`` with reason
            :data:`COVERAGE_NOT_ASSESSED_REASON` and its records are suppressed.
            Suppression is applied by whichever half is handed the calls: the
            ``calls_by_isolate`` door hands it to the producer, and the reader
            applies it to the records it has just read. Doing it on both paths
            is what keeps this parameter from being accepted and ignored.
        min_locus_coverage: Threshold, defaulting to
            :data:`MIN_LOCUS_COVERAGE_FRACTION`.

    Returns:
        ``sample_id -> records``, with a key for every manifest sample.

    Raises:
        DataContractError: Any named refusal in :func:`_run_real` when
            ``calls_by_isolate`` is supplied, or in :func:`load_regulator_variants`
            when the on-disk table is absent or unusable.
    """
    if calls_by_isolate is not None:
        # The compatibility door. Not the input: a caller that still has the
        # calls in scope gets the producer, so nothing that worked before this
        # split silently starts reading a table instead. The coverage arguments
        # are forwarded because the producer is the half that applies them.
        return _run_real(
            config,
            manifest,
            intermediate_root,
            calls_by_isolate=calls_by_isolate,
            reference_fasta=reference_fasta,
            locus_coverage=locus_coverage,
            min_locus_coverage=min_locus_coverage,
        )

    path = regulator_table_path(intermediate_root)
    if mode is RunMode.TEST:
        # Unchanged: the committed fixture, read exactly as before.
        records = load_regulator_variants(config, path)
        return group_by_sample(records, manifest)

    if not path.is_file():
        raise DataContractError(
            f"The regulator screen's stage input is absent: no {path.name} at "
            f"{path}. Stage 6's REAL input is that table - "
            f"{'/'.join(REGULATOR_TABLE_RELPATH)}, the stage input table named at "
            "docs/data_contract.md:114 - and nothing else. It is written by this "
            "stage's own producer, "
            "papipeline.stages.regulators.produce_regulator_variants, which "
            "papipeline/run.py calls after stage 6 writes its per-isolate calls "
            "to <intermediate>/stages/variants.tsv; run the `variants` stage of "
            "this pipeline rather than calling this stage on its own. Refusing "
            "here rather than returning an empty mapping, which would report the "
            "whole cohort as carrying no regulator variant. A caller that still "
            "holds the calls in memory may pass "
            "calls_by_isolate=variant_calls to invoke the producer directly; "
            "without either this stage has nothing to screen.",
            path=str(path),
            stage="regulators",
            missing="regulator_variants.tsv",
        )

    records = load_regulator_variants(
        config, path, intermediate_root=Path(intermediate_root)
    )
    if locus_coverage is not None:
        # The producer normally applied this threshold when it wrote the table,
        # so filtering again is idempotent. It is done anyway, because a caller
        # that supplies coverage evidence and reads a table written WITHOUT it
        # would otherwise be handed a locus reported as clean that was in fact
        # never read - the exact conflation
        # :data:`COVERAGE_NOT_ASSESSED_REASON` exists to prevent.
        records = [
            record for record in records
            if not_assessed_reason(
                locus_coverage.get((record.sample_id, record.gene)),
                min_locus_coverage,
            ) is None
        ]
    return group_by_sample(records, manifest)


#: Provenance sidecar written beside the table, so "the screen ran and found
#: nothing" is a recorded positive finding rather than an inference from an
#: absent row. See :func:`_write_screen_report`.
SCREEN_REPORT_NAME = "regulator_screen.json"

#: ``evidence_source`` for every call the screen emits. Names the mechanism, not
#: the run: a reader must be able to tell a reference-coordinate intersection
#: from an assembly-based call without consulting the manifest.
SCREEN_EVIDENCE_SOURCE = "stage6_call_x_reference_gff_cds"


def load_contig_sequences(fasta: Path) -> Dict[str, str]:
    """Every contig of the pinned reference, keyed by sequence name.

    Read once and passed down rather than re-read per isolate: the reference is
    6.26 Mb and the screen runs once per cohort member.

    Raises:
        DataContractError: The FASTA is absent, or holds no sequence. An empty
            parse would make every call fall outside every locus and report a
            cohort-wide clean screen - a fabricated negative on the study's
            central mechanisms.
    """
    path = Path(fasta)
    if not path.is_file():
        raise DataContractError(
            f"Pinned reference not found at {path}. The regulator screen "
            "intersects stage-6 calls with the reference's CDS features, so it "
            "needs both the reference sequence and its annotation. There is no "
            "fallback to the cohort's own assemblies: those sit on draft "
            "scaffolds in a different coordinate frame, and intersecting them "
            "would attribute variants to unrelated genes.",
            path=str(path),
        )

    sequences: Dict[str, str] = {}
    name: Optional[str] = None
    chunks: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks).upper()
                header = line[1:].split()
                name = header[0] if header else ""
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        sequences[name] = "".join(chunks).upper()

    if not sequences:
        raise DataContractError(
            f"No sequences parsed from {path}. Screening against an empty "
            "reference would place every call outside every screened locus and "
            "report the whole cohort as carrying no regulator variant, which is "
            "a fabricated negative rather than a result.",
            path=str(path),
        )
    return sequences


def describe_effect(
    variant_type: str, call: gff.VariantCall, interval: gff.GeneInterval
) -> str:
    """A precise, deterministic description of what the call does.

    ``effect`` has no declared vocabulary in the contract -
    ``docs/data_contract.md:114`` names the column and stops - so nothing here
    is validated against an enum and nothing downstream checks it. That is
    exactly why the string is built from measured quantities (codon number,
    signed length change, the observed alleles) rather than from prose: a claim
    nobody validates should at least be checkable by eye against its own row.

    Codon numbers come from :func:`gff.codon_number`, which is strand-aware.
    That is not pedantry: on the pinned reference ``PA0958`` (oprD) is on the
    minus strand, and counting from ``interval.start`` there labels its first
    codon 444 and its 368th codon 77 - the exact residue a truncating variant
    has to name correctly.
    """
    codon = gff.codon_number(interval, call.position)
    ref = (call.reference or "").upper()
    alt = (call.alternate or "").upper()
    label = gff.variant_label(call)

    if variant_type == VariantType.FRAMESHIFT.value:
        return f"frameshift_{len(alt) - len(ref):+d}nt_at_codon_{codon}"
    if variant_type == VariantType.INDEL.value:
        return f"inframe_indel_{len(alt) - len(ref):+d}nt_at_codon_{codon}"
    if variant_type == VariantType.PREMATURE_STOP.value:
        return f"stop_gained_at_codon_{codon}_{label}"
    if variant_type == VariantType.SNV.value:
        return f"substitution_at_codon_{codon}_{label}"
    return f"unclassified_{label}_at_codon_{codon}"


def variant_identifier(gene: str, call: gff.VariantCall) -> str:
    """A stable id for one call in one locus.

    Built from coordinates and alleles rather than a counter, so the same call
    yields the same id on a re-run and two runs of one cohort are diffable.
    ``variant`` is in the contract's REQUIRED tuple, so it must never be empty.
    """
    return "{gene}_{contig}_{pos}_{ref}_{alt}".format(
        gene=gene,
        contig=call.contig,
        pos=call.position,
        ref=(call.reference or "N").upper(),
        alt=(call.alternate or "N").upper(),
    )


def screen_calls(
    calls_by_isolate: Mapping[str, Sequence[Mapping[str, str]]],
    manifest: SampleManifest,
    *,
    intervals: Sequence[gff.GeneInterval],
    gene_by_tag: Mapping[str, str],
    specs: Mapping[str, "RegulatorSpec"],
    contigs: Mapping[str, str],
    locus_coverage: Optional[Mapping[Tuple[str, str], float]] = None,
    min_locus_coverage: float = MIN_LOCUS_COVERAGE_FRACTION,
) -> Tuple[List[RegulatorVariant], Dict[str, Any]]:
    """Intersect the cohort's calls with the screened loci.

    Returns the records plus an accounting dict. The accounting is returned
    rather than only logged because it is what makes an empty result
    interpretable: a caller receiving zero records can still say how many
    isolates were screened, how many loci were resolvable, and how many calls
    fell outside them.

    ``locus_coverage`` is optional and additive: ``(sample_id, gene) -> covered
    fraction``. When supplied, a pair whose coverage is below
    ``min_locus_coverage`` has its records **suppressed** and is counted under
    ``loci_not_assessed`` with the named reason. Suppression is the load-bearing
    part. Emitting a call from a locus the aligner only partly read would assert
    the sequence there, which is exactly the claim that is unsupported; and
    emitting nothing while saying nothing would report the locus as clean.
    """
    records: List[RegulatorVariant] = []
    unscreened: List[str] = []
    called_with_nothing: List[str] = []
    unassignable_contigs: Dict[str, int] = {}
    unparseable = 0
    outside_loci = 0
    out_of_scope: Dict[str, int] = {}
    not_assessed: Dict[str, Dict[str, Any]] = {}
    suppressed_by_coverage = 0

    for sample_id in manifest.sample_ids:
        rows = calls_by_isolate.get(sample_id)
        if rows is None:
            # Not "screened and clean": a cohort member the caller never passed
            # is an isolate whose calls do not exist. Refused in `_run_real`, so
            # this is belt and braces rather than the primary guard.
            unscreened.append(sample_id)
            continue
        if not rows:
            called_with_nothing.append(sample_id)

        for row in rows:
            contig = str(row.get("chrom") or "")
            try:
                position = int(str(row.get("pos") or ""))
            except ValueError:
                unparseable += 1
                continue

            if contig not in contigs:
                # The coordinate-frame violation this whole design exists to
                # catch: a draft-scaffold contig reaching a PAO1-coordinate
                # screen. Counted and surfaced, never dropped silently.
                key = contig or "(blank)"
                unassignable_contigs[key] = unassignable_contigs.get(key, 0) + 1
                continue

            call = gff.VariantCall(
                sample_id=sample_id,
                contig=contig,
                position=position,
                reference=str(row.get("ref") or ""),
                alternate=str(row.get("alt") or ""),
            )
            hits = gff.assign(call, intervals)
            if not hits:
                outside_loci += 1
                continue

            sequence = contigs[contig]
            for _name, interval in hits:
                # The gene symbol comes from config/regulators.tsv, never from
                # the GFF's Name=. For mexZ (PA2020) and dacB (PA3047) the GFF
                # carries the locus tag in Name= and no gene symbol at all, and
                # `gff.assign` returns `interval.name or interval.locus_tag`, so
                # trusting it would put "PA2020" in the gene column - which then
                # fails to join config.regulator("mexZ") and gets logged as "not
                # on the regulator screen list" for a locus that is on it.
                gene = gene_by_tag.get(interval.locus_tag)
                spec = specs.get(gene) if gene else None
                if spec is None:
                    continue

                coverage = (
                    None
                    if locus_coverage is None
                    else locus_coverage.get((sample_id, gene))
                )
                reason = not_assessed_reason(coverage, min_locus_coverage)
                if reason is not None:
                    key = f"{sample_id}:{gene}"
                    entry = not_assessed.setdefault(
                        key,
                        {
                            "sample_id": sample_id,
                            "gene": gene,
                            "coverage": coverage,
                            "min_locus_coverage": min_locus_coverage,
                            "reason": reason,
                            "suppressed_records": 0,
                        },
                    )
                    entry["suppressed_records"] += 1
                    suppressed_by_coverage += 1
                    continue

                variant_type = gff.classify(call, sequence, interval)
                if variant_type not in spec.variant_classes:
                    # Not emitted, but not invisible either. Every locus lists
                    # GENE_ABSENCE in screen_for and nothing can emit it from
                    # per-site calls, so a silently-dropped class would read as
                    # a screened-and-clean locus.
                    scope_key = f"{gene}:{variant_type}"
                    out_of_scope[scope_key] = out_of_scope.get(scope_key, 0) + 1
                    continue

                records.append(
                    RegulatorVariant(
                        sample_id=sample_id,
                        gene=gene,
                        variant=variant_identifier(gene, call),
                        variant_type=variant_type,
                        position=position,
                        reference=call.reference or None,
                        alternate=call.alternate or None,
                        effect=describe_effect(variant_type, call, interval),
                        mechanism=spec.mechanism,
                        # No calibrated confidence exists for a call derived this
                        # way, and inventing one would be a number nothing
                        # downstream can defend. Absent, not guessed.
                        confidence=None,
                        evidence_source=SCREEN_EVIDENCE_SOURCE,
                        # The locus's declared ceiling, never upgraded: a variant
                        # call in a locus cannot support a stronger claim than
                        # config/regulators.tsv allows for that locus.
                        call_status=spec.evidence_level,
                    )
                )

    records.sort(
        key=lambda r: (
            r.sample_id,
            r.gene,
            r.position if r.position is not None else 0,
            r.variant,
        )
    )
    accounting: Dict[str, Any] = {
        "isolates_in_manifest": len(manifest.sample_ids),
        "isolates_screened": len(manifest.sample_ids) - len(unscreened),
        "isolates_called_with_nothing": sorted(called_with_nothing),
        "isolates_not_screened": sorted(unscreened),
        "loci_screened": len(intervals),
        "calls_examined": sum(len(r) for r in calls_by_isolate.values()),
        "calls_outside_screened_loci": outside_loci,
        "calls_out_of_declared_scope": dict(sorted(out_of_scope.items())),
        "calls_unparseable_position": unparseable,
        "calls_on_unknown_contig": dict(sorted(unassignable_contigs.items())),
        "records_emitted": len(records),
        # Absent entirely when no coverage evidence was supplied, so a reader
        # can tell "nothing was below the threshold" from "nothing was
        # measured". The two are the whole point of this key.
        **({
            "min_locus_coverage": min_locus_coverage,
            "coverage_evidence": "supplied" if locus_coverage is not None
            else "absent",
            "loci_not_assessed": dict(sorted(not_assessed.items())),
            "records_suppressed_below_coverage": suppressed_by_coverage,
        } if locus_coverage is not None else {}),
    }
    return records, accounting


def _write_screen_report(path: Path, payload: Mapping[str, Any]) -> Path:
    """Write the provenance sidecar. Best effort, and never fatal.

    The table is the result; the sidecar is what makes a zero-row table
    interpretable. Losing the sidecar must not lose the run, but losing it must
    be logged, because a table with no sidecar is the ambiguous case this whole
    mechanism exists to remove.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        path.write_text(
            json.dumps(dict(payload), indent=2, sort_keys=True, default=str) + "\n",
            encoding="utf-8",
        )
    except OSError as exc:  # pragma: no cover - filesystem dependent
        LOGGER.warning(
            "Stage 6: could not write the screen provenance sidecar at %s (%s). "
            "The table is intact, but a zero-row table is now ambiguous between "
            "'screened and clean' and 'not run'.",
            path,
            exc,
        )
    return path


def read_screen_report(intermediate_root: Path) -> Optional[Dict[str, Any]]:
    """The sidecar's contents, or ``None`` when there is no sidecar.

    ``None`` is the load-bearing value: it is what distinguishes "this screen
    never ran" from "this screen ran and found nothing". A missing file and a
    malformed one are both ``None``, because neither carries evidence that the
    screen executed.
    """
    path = Path(intermediate_root) / "regulators" / SCREEN_REPORT_NAME
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        LOGGER.warning(
            "Stage 6: screen provenance sidecar at %s is unreadable (%s). Treating "
            "it as absent, so a zero-row table cannot be read as evidence the "
            "screen ran.", path, exc,
        )
        return None
    return payload if isinstance(payload, dict) else None


def _run_real(
    config: PipelineConfig,
    manifest: SampleManifest,
    intermediate_root: Path,
    *,
    calls_by_isolate: Optional[Mapping[str, Sequence[Mapping[str, str]]]] = None,
    reference_fasta: Optional[Path] = None,
    locus_coverage: Optional[Mapping[Tuple[str, str], float]] = None,
    min_locus_coverage: float = MIN_LOCUS_COVERAGE_FRACTION,
) -> Dict[str, List[RegulatorVariant]]:
    """REAL: derive the regulator table from this run's stage-6 calls.

    Modelled on ``stages.sv._run_real``, which has the same shape for the same
    reason: a stage that owns a scientific screen refuses by name when an input
    is missing, derives its own table in REAL, and keeps a fixture read for TEST
    only.

    The intersection and the classification are
    :mod:`papipeline.adapters.gff`'s - already written, already unit-tested, and
    until now never called by anything. This function supplies the two things it
    lacked: a caller, and a writer.

    **No empty table is ever produced by accident.** A missing input raises
    instead; a screen that genuinely finds nothing writes a zero-row table *and*
    a sidecar recording how many isolates and loci it examined, so the two cases
    cannot be confused by a downstream reader.
    """
    if calls_by_isolate is None:
        raise DataContractError(
            "The regulator screen was asked to run without stage-6 calls. It "
            "filters this run's variant calls down to the loci in "
            "config/regulators.tsv and has no other input - there is no file to "
            "fall back to, because the calls are held in the caller's scope and "
            "re-reading them would be a second read free to disagree with the "
            "cohort merge. Call papipeline.stages.regulators."
            "produce_regulator_variants(calls_by_isolate=...) from the wiring in "
            "papipeline/run.py rather than this function: `run` is the reader of "
            "the on-disk stage input table named at docs/data_contract.md:114. "
            "Refusing here rather than writing an empty table, which would "
            "report the whole cohort as carrying no regulator variant.",
            stage="regulators", missing="calls_by_isolate",
        )

    missing_samples = [s for s in manifest.sample_ids if s not in calls_by_isolate]
    if missing_samples:
        # Rule 5: sample IDs must match across stages, and a mismatch is a hard
        # failure with a clear message, never a silent drop. Seeding every
        # manifest sample is the caller's job - see the stage-6 docstring on why
        # the empty entries are load-bearing for the merge's denominator - so an
        # absent one means the two disagree about who is in the study.
        raise DataContractError(
            f"{len(missing_samples)} manifest sample(s) have no entry in the "
            f"stage-6 calls: {missing_samples[:10]}. The screen cannot tell an "
            "isolate that was screened and found nothing from one that was never "
            "screened, and recording the latter as clean would report an "
            "unanalysable genome as a negative. Stage 6 seeds every manifest "
            "sample for exactly this reason.",
            missing=",".join(missing_samples[:10]), n_missing=len(missing_samples),
        )

    foreign = sorted(set(calls_by_isolate) - set(manifest.sample_ids))
    if foreign:
        raise DataContractError(
            f"Stage-6 calls reference {len(foreign)} sample(s) the manifest "
            f"never names: {foreign[:10]}. The manifest decides cohort "
            "membership; a call from an isolate it does not name cannot be "
            "screened, because it would add a member the study does not contain.",
            foreign=",".join(foreign[:10]), n_foreign=len(foreign),
        )

    reference = Path(reference_fasta) if reference_fasta else config.reference_fasta()
    gff_path = config.reference_gff()

    specs = dict(config.regulators)
    screenable = {g: s for g, s in specs.items() if s.screenable}
    unresolvable = sorted(g for g, s in specs.items() if not s.screenable)

    if not screenable:
        # Would write a zero-row table and report a cohort-wide clean screen.
        raise DataContractError(
            "No locus in config/regulators.tsv carries a PAO1 locus_tag, so none "
            "can be intersected with reference-coordinate calls. An empty result "
            "here is a fabricated negative on the study's central mechanisms, not "
            "a finding. Resolve the tags with scripts/resolve_locus_tags.py.",
            n_loci=len(specs),
        )

    # The two refusals below are raised by the adapter and are already covered by
    # tests/unit/test_gff_interval_lookup.py: a missing GFF, and a configured tag
    # the pinned annotation does not contain.
    intervals_by_tag = gff.load_gene_intervals(
        gff_path, [s.locus_tag for s in screenable.values()]
    )
    intervals = [intervals_by_tag[s.locus_tag] for s in screenable.values()]
    gene_by_tag = {s.locus_tag: g for g, s in screenable.items() if s.locus_tag}

    contigs = load_contig_sequences(reference)

    records, accounting = screen_calls(
        calls_by_isolate, manifest,
        intervals=intervals, gene_by_tag=gene_by_tag, specs=specs,
        contigs=contigs,
        locus_coverage=locus_coverage, min_locus_coverage=min_locus_coverage,
    )

    target = Path(intermediate_root) / "regulators"
    target.mkdir(parents=True, exist_ok=True)
    table = target / "regulator_variants.tsv"

    write_tsv(
        table,
        [r.to_row() for r in records],
        list(REGULATOR_COLUMNS),
        [
            "Stage 6 regulator screen: this run's stage-6 calls intersected with "
            f"the CDS features of {gff_path.name}.",
            f"mode=REAL screened_isolates={accounting['isolates_screened']}"
            f"/{accounting['isolates_in_manifest']}"
            f" loci_screened={accounting['loci_screened']}"
            f" records={accounting['records_emitted']}",
            "A zero-row body means the screen ran and found nothing. The "
            f"sidecar {SCREEN_REPORT_NAME} records what it examined; no sidecar "
            "means the screen did not run.",
        ],
    )
    _write_screen_report(
        target / SCREEN_REPORT_NAME,
        {
            "mode": RunMode.REAL.value,
            "reference_fasta": str(reference),
            "reference_gff": str(gff_path),
            "evidence_source": SCREEN_EVIDENCE_SOURCE,
            # Named `loci_names`, not `loci_screened`: the accounting below
            # already carries `loci_screened` as a *count*, and the two collided -
            # the count overwrote the list, and a sidecar claiming to list the
            # screened loci held an integer.
            "loci_names": sorted(gene_by_tag.values()),
            "loci_without_pao1_counterpart": unresolvable,
            # Every (isolate, locus) pair, with the fraction covered and the
            # verdict. Present whatever the coverage input was, because the
            # pairs with no coverage evidence are exactly the ones a reader
            # would otherwise have to assume were measured and passed.
            "locus_coverage": locus_coverage_rows(
                list(manifest.sample_ids), sorted(gene_by_tag.values()),
                locus_coverage, min_locus_coverage=min_locus_coverage,
            ),
            **accounting,
        },
    )

    if unresolvable:
        # Reported, never dropped silently. `mexS` is the expected case and
        # config/regulators.tsv documents why: it has no PAO1 counterpart, so
        # there is no reference interval to intersect against. An absent tag is
        # not a screened-and-clean locus.
        LOGGER.warning(
            "Stage 6: %d configured locus/loci carry no PAO1 locus_tag and were "
            "NOT screened: %s. An empty locus_tag means the locus is absent from "
            "PAO1; it does not mean the locus was screened and found clean.",
            len(unresolvable), ",".join(unresolvable),
        )
    if accounting.get("coverage_evidence") == "supplied":
        LOGGER.info(
            "Stage 6: locus coverage threshold %s; %d (isolate, locus) pair(s) "
            "below it and reported not_assessed with reason %s",
            accounting["min_locus_coverage"],
            len(accounting["loci_not_assessed"]),
            COVERAGE_NOT_ASSESSED_REASON,
        )
    if accounting.get("loci_not_assessed"):
        LOGGER.warning(
            "Stage 6: %d (isolate, locus) pair(s) were NOT assessed because the "
            "alignment covered less than %.0f%% of them: %s. Their records were "
            "suppressed and they are reported `not_assessed`, which is a "
            "statement that the locus was not read, NOT that it is clean.",
            len(accounting["loci_not_assessed"]),
            100.0 * accounting["min_locus_coverage"],
            ",".join(sorted(accounting["loci_not_assessed"])[:10]),
        )
    LOGGER.info(
        "Stage 6: screened %d/%d isolates against %d loci from %s | %d regulator "
        "variant call(s) | %d call(s) outside the screened loci",
        accounting["isolates_screened"], accounting["isolates_in_manifest"],
        accounting["loci_screened"], gff_path.name,
        accounting["records_emitted"], accounting["calls_outside_screened_loci"],
    )
    if accounting["calls_on_unknown_contig"]:
        LOGGER.warning(
            "Stage 6: %d call(s) sat on contigs absent from the pinned reference "
            "and could not be assigned to any locus: %s. Stage-6 calls are PAO1 "
            "coordinates, so a draft-scaffold contig here means the coordinate "
            "frames disagree.",
            sum(accounting["calls_on_unknown_contig"].values()),
            ",".join(sorted(accounting["calls_on_unknown_contig"])),
        )

    return group_by_sample(records, manifest)
