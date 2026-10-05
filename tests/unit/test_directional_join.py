"""The manifest-to-phenotype join is directional, and deliberately asymmetric.

Ticket 06. This is the rule most easily got wrong, because the intuitive thing
is a symmetric 1:1 assertion and the intuitive thing is *wrong here*:

* the manifest is the authority on what exists;
* the phenotype table is external evidence about it, and may legitimately
  describe isolates that were never downloaded.

So a phenotype row with no assembly is **excluded**, not a failure. A manifest
genome with no phenotype row, or two, **is** a failure - the pipeline cannot
know which cohort was meant. And a genome whose phenotype is unusable is
excluded with a reason, because "I" and "not determined" are not resistance and
must not be pooled into it.

Inside the pipeline any sample-ID mismatch between assemblies, annotations,
variant tables and the final phenotype table remains a hard failure. The
directional rule governs the manifest edge only.

Written before the implementation.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.errors import DataContractError
from papipeline.join import (
    ExclusionReason,
    JoinedCohort,
    join_manifest_to_phenotype,
)

G1, G2, G3 = "GCA_000000001.1", "GCA_000000002.1", "GCA_000000003.1"


def _phenotype(sample_id, value="R", mic=8.0):
    from papipeline.models import Phenotype, PhenotypeCall

    return PhenotypeCall(
        sample_id=sample_id,
        antibiotic="imipenem",
        phenotype=Phenotype(value),
        mic=mic,
    )


# ---------------------------------------------------------------------------
# the hard failures: the manifest edge
# ---------------------------------------------------------------------------


def test_a_genome_with_no_phenotype_row_is_a_hard_failure():
    with pytest.raises(DataContractError) as excinfo:
        join_manifest_to_phenotype([G1, G2], [_phenotype(G1)])
    assert G2 in str(excinfo.value), "the missing genome must be named"


def test_a_genome_with_two_phenotype_rows_is_a_hard_failure():
    calls = [_phenotype(G1), _phenotype(G2)]
    duplicate = _phenotype(G1)
    calls.append(duplicate)
    with pytest.raises(DataContractError) as excinfo:
        join_manifest_to_phenotype([G1, G2], calls)
    assert G1 in str(excinfo.value)


def test_duplicate_rows_are_detected_before_anything_else():
    """The ambiguity must be reported even when other problems exist."""
    with pytest.raises(DataContractError) as excinfo:
        join_manifest_to_phenotype([G1, G2], [_phenotype(G1), _phenotype(G1)])
    assert "duplicate" in str(excinfo.value).lower() or "more than one" in str(excinfo.value).lower()


# ---------------------------------------------------------------------------
# the exclusions: the phenotype edge, and unusable categories
# ---------------------------------------------------------------------------


def test_a_phenotype_row_with_no_assembly_is_excluded_not_failed():
    calls = [_phenotype(G1), _phenotype("GCA_000000099.1")]
    cohort = join_manifest_to_phenotype([G1], calls)  # G1 has a measured MIC
    assert [c.call.sample_id for c in cohort.joined] == [G1]
    assert cohort.excluded["GCA_000000099.1"].reason is ExclusionReason.NO_ASSEMBLY


@pytest.mark.parametrize("value", ["I", "SDD", "ND"])
def test_an_unusable_category_is_excluded_with_its_reason(value: str):
    calls = [_phenotype(G1), _phenotype(G2, value=value, mic=None)]
    cohort = join_manifest_to_phenotype([G1, G2], calls)
    assert [c.call.sample_id for c in cohort.joined] == [G1]
    assert cohort.excluded[G2].reason is ExclusionReason.UNUSABLE_CATEGORY
    assert value in cohort.excluded[G2].detail or "categor" in cohort.excluded[G2].detail.lower()


def test_a_genome_whose_row_is_missing_from_the_phenotype_table_fails_not_excludes():
    """It is already covered above; this pins the asymmetry explicitly.

    NO_ASSEMBLY is an exclusion, MISSING_PHENOTYPE is not reachable, because a
    manifest genome with no row is a hard failure. The distinction is the point
    of the whole ticket.
    """
    with pytest.raises(DataContractError):
        join_manifest_to_phenotype([G1, G2], [_phenotype(G1)])


# ---------------------------------------------------------------------------
# the continuous trait decides inclusion
# ---------------------------------------------------------------------------


def test_not_determined_is_excluded_even_if_an_mic_survived_somewhere():
    """SDD/ND mean not determined. The parser refuses MIC alongside them; the
    join refuses it too rather than trusting that the parser ran."""
    calls = [_phenotype(G1), _phenotype(G2, value="ND", mic=8.0)]
    cohort = join_manifest_to_phenotype([G1, G2], calls)
    assert [c.call.sample_id for c in cohort.joined] == [G1]
    assert cohort.excluded[G2].reason is ExclusionReason.UNUSABLE_CATEGORY


@pytest.mark.parametrize("value", ["I", "SDD"])
def test_intermediate_and_susceptible_dose_dependent_rows_keep_their_mic(value: str):
    """SDD is NOT "not determined".

    Susceptible-dose dependent is a *determination*, and under EUCAST an SDD
    call is derived from an MIC that falls inside the susceptible range. An MIC
    alongside SDD is therefore the normal combination, not a contradiction -
    and the parser deliberately permits it. Excluding SDD would discard real
    measurements, and would do so on a false premise.
    """
    calls = [_phenotype(G1, mic=8.0), _phenotype(G2, value=value, mic=1.0)]
    cohort = join_manifest_to_phenotype([G1, G2], calls)
    assert [c.call.sample_id for c in cohort.joined] == [G1, G2]
    assert cohort.joined[1].log2_mic == pytest.approx(0.0)


def test_not_determined_is_excluded_because_it_is_not_a_determination():
    """The only category that can never enter: ND, which by definition has
    no measurement. The parser refuses an MIC alongside it, so the join's
    exclusion is a defence, not the primary control."""
    calls = [_phenotype(G1, mic=8.0), _phenotype(G2, value="ND", mic=None)]
    cohort = join_manifest_to_phenotype([G1, G2], calls)
    assert [c.call.sample_id for c in cohort.joined] == [G1]
    assert cohort.excluded[G2].reason is ExclusionReason.UNUSABLE_CATEGORY


def test_a_row_with_a_measured_mic_is_usable_whatever_its_category():
    """`I` with a measurement is analysable as a continuous trait.

    This is why the trait change matters: under a binary R/S outcome an
    intermediate with an MIC would be discarded, and the MIC - the actual
    measurement - thrown away with it.
    """
    calls = [_phenotype(G1), _phenotype(G2, value="I", mic=2.4)]
    cohort = join_manifest_to_phenotype([G1, G2], calls, log_base=2)
    assert [c.call.sample_id for c in cohort.joined] == [G1, G2]
    assert cohort.joined[1].log2_mic == pytest.approx(1.2630, abs=1e-4)


def test_a_row_with_a_category_but_no_mic_is_excluded_when_a_trait_is_required():
    calls = [_phenotype(G1, mic=8.0), _phenotype(G2, value="S", mic=None)]
    cohort = join_manifest_to_phenotype([G1, G2], calls, require_trait=True)
    assert [c.call.sample_id for c in cohort.joined] == [G1]
    assert cohort.excluded[G2].reason is ExclusionReason.NO_TRAIT


def test_a_category_only_row_is_kept_when_no_trait_is_required():
    calls = [_phenotype(G1, mic=8.0), _phenotype(G2, value="S")]
    cohort = join_manifest_to_phenotype([G1, G2], calls, require_trait=False)
    assert [c.call.sample_id for c in cohort.joined] == [G1, G2]


# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------


def test_exclusion_counts_are_reported():
    calls = [
        _phenotype(G1, mic=8.0),
        _phenotype(G2, value="ND"),
        _phenotype("GCA_000000099.1", mic=1.0),
    ]
    cohort = join_manifest_to_phenotype([G1, G2], calls)
    assert cohort.counts.by_reason[ExclusionReason.NO_ASSEMBLY] == 1
    assert cohort.counts.by_reason[ExclusionReason.UNUSABLE_CATEGORY] == 1
    assert cohort.counts.total_excluded == 2
    assert cohort.counts.total_joined == 1
    assert cohort.counts.total_in_manifest == 2
    # the summary promises counts; the per-sample detail lives in the table
    summary = cohort.counts.summary()
    assert "2 excluded" in summary
    assert "no_assembly: 1" in summary
    assert "unusable_category: 1" in summary


def test_each_exclusion_keeps_its_own_detail():
    calls = [_phenotype(G1), _phenotype(G2, value="ND", mic=None),
             _phenotype("GCA_000000099.1", mic=1.0)]
    cohort = join_manifest_to_phenotype([G1, G2], calls)
    assert "ND" in cohort.excluded[G2].detail
    assert "no assembly" in cohort.excluded["GCA_000000099.1"].detail


def test_the_exclusion_table_is_writable_and_round_trips(tmp_path: Path):
    import csv

    calls = [_phenotype(G1, mic=8.0), _phenotype(G2, value="SDD", mic=None)]
    cohort = join_manifest_to_phenotype([G1, G2], calls)
    path = tmp_path / "exclusions.tsv"
    cohort.write_exclusions(path)
    rows = list(csv.DictReader(path.open(encoding="utf-8"), delimiter="\t"))
    assert [r["sample_id"] for r in rows] == [G2]
    assert rows[0]["reason"]


def test_an_empty_cohort_is_reported_as_empty_not_as_an_error():
    cohort = join_manifest_to_phenotype([], [])
    assert not cohort.joined
    assert cohort.counts.total_joined == 0
    assert cohort.counts.total_excluded == 0


def test_the_join_is_ordered_by_the_manifest_not_by_the_file():
    """A phenotype file sorted differently must not reorder the cohort."""
    calls = [_phenotype(G3, mic=1.0), _phenotype(G1, mic=8.0), _phenotype(G2, mic=2.0)]
    cohort = join_manifest_to_phenotype([G1, G2, G3], calls)
    assert [c.call.sample_id for c in cohort.joined] == [G1, G2, G3]
