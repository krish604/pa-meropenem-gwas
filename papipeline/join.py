"""Directional join: manifest (what exists) to phenotype (evidence about it).

The asymmetry is the whole point, and it is easy to get wrong because a
symmetric 1:1 assertion is the intuitive thing to write.

* **The manifest is the authority on what exists.** A genome with no phenotype
  row, or two, is a hard failure: the pipeline cannot know which cohort was
  meant, and guessing is how a GWAS ends up on a different set than the report
  describes.
* **The phenotype table is external evidence and may be larger than reality.**
  A row naming an isolate that was never downloaded is excluded, not a failure.
  That is a normal situation - you have a phenotype table for 966 isolates and
  835 assemblies in the dataset this was written against. The join does not
  read that number: it reads whatever the manifest it is given contains, and
  the figure is here only to say which cohort the notes describe.
* **An unusable category is excluded, not failed.** A row excluded this way
  keeps its reason.
* **A measured MIC makes ``I`` and ``SDD`` rows analysable.** Both are real
  determinations, not missing results. ``I`` is intermediate; ``SDD`` is
  susceptible-dose dependent, a separate CLSI category whose EUCAST
  equivalent is ``I`` again (Susceptible, increased exposure). An SDD call is
  derived from an MIC inside the susceptible range, so MIC+SDD is the normal
  combination. Under the continuous trait they contribute log2(MIC) and are
  kept; discarding them would throw away the only quantitative thing in the
  row. ``ND`` is the only category that can never enter: not determined carries
  no measurement by definition.

Nothing here is a fuzzy join. Sample identifiers match exactly or the run stops.

This module governs the manifest edge only. A sample-ID mismatch discovered
*inside* the pipeline, between assemblies, annotations, variant tables and the
final phenotype table, remains a hard failure in the stage that finds it.
"""

from __future__ import annotations

import csv
import enum
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from .errors import DataContractError
from .models import Phenotype, PhenotypeCall
from .stages.phenotype import interpret_mic

#: Categories that can never enter the analysis.
#:
#: ``ND`` means not determined. It is the only category that cannot carry a
#: measurement: the parser refuses an MIC alongside it, so this exclusion is a
#: defence rather than the primary control.
UNUSABLE_ALWAYS = (Phenotype.ND,)

#: Categories that are real determinations, and so keep their MIC when they
#: have one.
#:
#: ``I`` is intermediate in CLSI: a real determination that often carries a
#: measured MIC, and the obvious case here.
#:
#: ``SDD`` is *susceptible-dose dependent*, a distinct **CLSI** category, and it
#: is equally a determination rather than a missing result. The EUCAST
#: equivalent of SDD is ``I`` again - "Susceptible, increased exposure" - so
#: under EUCAST an SDD-like call arrives as I with the same meaning. Under
#: CLSI an SDD call is derived from an MIC falling inside the susceptible
#: range, so MIC+SDD is the normal combination and the parser deliberately
#: permits it. Treating SDD as not-determined would discard real measurements on
#: a false premise.
UNUSABLE_UNLESS_MEASURED = (Phenotype.I, Phenotype.SDD)


class ExclusionReason(str, enum.Enum):
    """Why a row did not enter the analysis. Never a free-text excuse."""

    NO_ASSEMBLY = "no_assembly"
    UNUSABLE_CATEGORY = "unusable_category"
    NO_TRAIT = "no_continuous_trait"


@dataclass(frozen=True)
class Exclusion:
    """One row that did not enter the analysis, and why."""

    sample_id: str
    reason: ExclusionReason
    detail: str

    def to_row(self) -> Dict[str, str]:
        return {
            "sample_id": self.sample_id,
            "reason": self.reason.value,
            "detail": self.detail,
        }


@dataclass(frozen=True)
class JoinedSample:
    """A manifest genome with its phenotype, ready for the analysis."""

    sample_id: str
    call: PhenotypeCall
    log2_mic: Optional[float]


@dataclass
class JoinCounts:
    """What the run did with what it was given."""

    total_in_manifest: int = 0
    total_joined: int = 0
    total_excluded: int = 0
    by_reason: Dict[ExclusionReason, int] = field(default_factory=dict)

    def record(self, reason: ExclusionReason) -> None:
        self.by_reason[reason] = self.by_reason.get(reason, 0) + 1
        self.total_excluded += 1

    def summary(self) -> str:
        parts = [
            f"Manifest {self.total_in_manifest} samples; "
            f"{self.total_joined} entered the analysis; "
            f"{self.total_excluded} excluded."
        ]
        for reason in sorted(self.by_reason, key=lambda r: r.value):
            count = self.by_reason[reason]
            parts.append(f"{reason.value}: {count}")
        return " ".join(parts)


@dataclass(frozen=True)
class JoinedCohort:
    """The analysis cohort, plus everything that was left out and why."""

    joined: Tuple[JoinedSample, ...]
    excluded: Mapping[str, Exclusion]
    counts: JoinCounts

    def sample_ids(self) -> List[str]:
        return [s.sample_id for s in self.joined]

    def write_exclusions(self, path: Path) -> Path:
        """Write the exclusion table. One row per excluded sample."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=("sample_id", "reason", "detail"),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for sample_id in sorted(self.excluded):
                writer.writerow(self.excluded[sample_id].to_row())
        return path


def _index_calls(
    calls: Sequence[PhenotypeCall], manifest_ids: Sequence[str]
) -> Dict[str, PhenotypeCall]:
    """Index phenotype rows by sample, failing loudly on ambiguity.

    Two rows for one sample is a hard failure here, before any categorisation
    happens, so the ambiguity is reported even when other problems exist too.
    """
    index: Dict[str, PhenotypeCall] = {}
    for call in calls:
        previous = index.get(call.sample_id)
        if previous is not None:
            raise DataContractError(
                "More than one phenotype row for a sample in the manifest; the "
                "cohort is ambiguous and the run cannot continue",
                sample_id=call.sample_id,
                first_source=previous.source,
                second_source=call.source,
            )
        index[call.sample_id] = call
    return index


def join_manifest_to_phenotype(
    manifest_ids: Sequence[str],
    calls: Iterable[PhenotypeCall],
    *,
    log_base: int = 2,
    require_trait: bool = True,
) -> JoinedCohort:
    """Join a manifest to its phenotype table, directionally.

    Args:
        manifest_ids: Sample identifiers, in manifest order. The result follows
            this order, not the order of the phenotype file.
        calls: Phenotype rows for one antibiotic.
        log_base: Base of the continuous trait.
        require_trait: When True, a row must carry a measured MIC to enter the
            analysis. Set False where the categorical call is wanted instead.

    Raises:
        DataContractError: A manifest genome has no phenotype row, or more than
            one. These are failures, not exclusions.
    """
    manifest_ids = list(manifest_ids)
    manifest_set = set(manifest_ids)
    calls = list(calls)
    index = _index_calls(calls, manifest_ids)

    missing = [s for s in manifest_ids if s not in index]
    if missing:
        raise DataContractError(
            "Every genome in the manifest must have exactly one phenotype row; "
            f"{len(missing)} have none: {', '.join(missing[:10])}"
            + (" ..." if len(missing) > 10 else ""),
            missing=len(missing),
            first_missing=missing[0],
        )

    joined: List[JoinedSample] = []
    excluded: Dict[str, Exclusion] = {}
    counts = JoinCounts(total_in_manifest=len(manifest_ids))

    for call in calls:
        if call.sample_id not in manifest_set:
            excluded[call.sample_id] = Exclusion(
                sample_id=call.sample_id,
                reason=ExclusionReason.NO_ASSEMBLY,
                detail=(
                    "phenotype row names an isolate with no assembly in the "
                    "manifest; excluded, not an error, because the phenotype "
                    "table legitimately describes more isolates than were "
                    "downloaded"
                ),
            )
            counts.record(ExclusionReason.NO_ASSEMBLY)

    for sample_id in manifest_ids:
        call = index[sample_id]
        if call.phenotype in UNUSABLE_ALWAYS or (
            call.phenotype in UNUSABLE_UNLESS_MEASURED and call.mic is None
        ):
            excluded[sample_id] = Exclusion(
                sample_id=sample_id,
                reason=ExclusionReason.UNUSABLE_CATEGORY,
                detail=(
                    f"category {call.phenotype.value} is not a susceptibility "
                    f"determination"
                    + (
                        " and carries no measured MIC"
                        if call.phenotype in UNUSABLE_UNLESS_MEASURED
                        else " (an MIC alongside it would be meaningless)"
                    )
                ),
            )
            counts.record(ExclusionReason.UNUSABLE_CATEGORY)
            continue

        trait = interpret_mic(call.mic, base=log_base)
        if require_trait and trait is None:
            excluded[sample_id] = Exclusion(
                sample_id=sample_id,
                reason=ExclusionReason.NO_TRAIT,
                detail=(
                    "no measured MIC, so there is no continuous trait to model"
                ),
            )
            counts.record(ExclusionReason.NO_TRAIT)
            continue

        joined.append(JoinedSample(sample_id=sample_id, call=call, log2_mic=trait))

    counts.total_joined = len(joined)
    return JoinedCohort(joined=tuple(joined), excluded=excluded, counts=counts)
