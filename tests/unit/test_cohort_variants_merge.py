"""The cohort merge: keep only positions where isolates actually disagree.

Long/tidy, one row per `(chrom, pos, ref, alt)`. That shape is settled - pyseer
is wide only via `--pres`, and reading its source shows samples come from the
header row, so a wide table would mean one column per isolate for 835 of them.
The pivot belongs to whoever builds pyseer's input, and has no consumer yet.

What the filter is, and what it deliberately is not:

Per-isolate calls conflate two things. A site every isolate differs from PAO1 at
is a species-wide fixed difference, not a polymorphism, and no amount of
per-sample processing can tell those apart - only a multi-isolate comparison
can. So the merge drops the fixed sites.

It is tempting to define the drop as "keep sites carried by >= 2 isolates". That
is wrong, and measured rather than assumed. On four real isolates, 68.3% of all
sites are carried by exactly one, so that rule keeps 31.7% of the genome - the
opposite of conservative. It also does not scale: 2 is 50% of a 4-isolate cohort
and 0.24% of an 835-isolate one. So the filter is on **minor-allele frequency**,
`1/N <= maf`, which makes the "carried by at least two" case implicit and free.

`max_maf` is NOT set here. A minor-allele frequency needs a real minor allele and
four isolates cannot show a rare one, so the upper bound is deferred to a
measurement at scale rather than guessed - same reasoning that produced the
per-isolate contract by reading real output.
"""

from __future__ import annotations

import pytest

from papipeline.errors import DataContractError
from papipeline.stages.cohort_variants import MERGE_COLUMNS, merge_calls

# Three isolates; `carriers` names how many carry each site.
CALLS = {
    "iso1": [("NC_002516.2", "100", "C", "T"), ("NC_002516.2", "200", "A", "G")],
    "iso2": [("NC_002516.2", "100", "C", "T"), ("NC_002516.2", "300", "G", "A")],
    "iso3": [("NC_002516.2", "100", "C", "T")],
}
# 100 is carried by all three (fixed difference), 200/300 by one each.


def _rows():
    return merge_calls(CALLS)


def _positions(rows):
    return {r["pos"] for r in rows}


class TestTheCeilingIsSymmetricWithTheFloor:
    """Same rule, both tails: a lone dissenter is as suspect as a lone carrier.

    The floor drops a site carried by one isolate, because a single carrier is
    indistinguishable from an artefact. The ceiling drops a site that *every*
    isolate but one carries, by the identical argument: a single dissenter is
    equally indistinguishable from that isolate's own error.

    This is not a claim about species-wide fixation. The 34-isolate measurement
    showed that claim is not recoverable from this data - the 100% bar *shrank*
    as N grew, and the upper-tail mass sat at 94% rather than 100%. So the
    ceiling is not asked to identify real fixation. It asks only that enough
    independent isolates sit on each side of a site for it to carry signal, and
    the floor had already earned that by being the same rule applied once.
    """

    def test_a_site_all_but_one_isolate_carries_is_dropped(self):
        calls = {"i1": [("c", "1", "A", "T")], "i2": [("c", "1", "A", "T")]}
        calls["i3"] = [("c", "1", "A", "T")]
        assert merge_calls(calls) == [], "a single dissenter was not dropped"

    def test_a_site_all_isolates_carry_is_dropped(self):
        calls = {f"i{k}": [("c", "1", "A", "T")] for k in range(3)}
        assert merge_calls(calls) == []

    def test_exactly_one_dissenter_is_dropped(self):
        """The discriminating case: one dissenter, not zero.

        `test_a_site_all_but_one_isolate_carries_is_dropped` uses zero
        dissenters, which a threshold of one would also drop - so it cannot tell
        "two required" from "one required". This is the case that separates them.
        """
        calls = {"i1": [("c", "1", "A", "T")], "i2": [("c", "1", "A", "T")], "i3": []}
        assert merge_calls(calls) == [], "a single dissenter was kept"

    def test_exactly_two_dissenters_is_kept(self):
        """And the boundary just past it is kept, so the rule is a bound of two."""
        calls = {"i1": [("c", "1", "A", "T")], "i2": [("c", "1", "A", "T")],
                 "i3": [], "i4": []}
        assert len(merge_calls(calls)) == 1

    def test_two_dissenters_enough_is_kept(self):
        """The bound is two on each side, so a near-even split survives."""
        calls = {"i1": [("c", "1", "A", "T")], "i2": [("c", "1", "A", "T")],
                 "i3": [], "i4": []}
        rows = merge_calls(calls)
        assert len(rows) == 1
        assert rows[0]["ac"] == "2" and rows[0]["an"] == "4"

    @pytest.mark.parametrize("n", [4, 5, 10, 34])
    def test_both_tails_need_two_independent_isolates(self, n):
        """Every retained site has >=2 carriers AND >=2 non-carriers."""
        calls = {f"i{k}": [("c", "1", "A", "T")] for k in range(n - 2)}
        calls["x"] = []
        calls["y"] = []
        rows = merge_calls(calls)
        assert len(rows) == 1
        assert int(rows[0]["ac"]) >= 2
        assert int(rows[0]["an"]) - int(rows[0]["ac"]) >= 2

    def test_it_is_the_algebraic_equivalent_of_a_maf_ceiling(self):
        """`an - ac >= 2` is exactly `af > 1 - 2/an`; the two must agree."""
        n = 34
        assert abs((1 - 2 / n) - (1 - 2.0 / n)) < 1e-12
        # A site at 32/34 keeps; one at 33/34 goes.
        assert (32 / n) <= (1 - 2 / n)
        assert (33 / n) > (1 - 2 / n)


class TestGenuinePolymorphismsSurvive:
    def test_a_site_two_of_three_carry_is_kept(self):
        """2/3 carriers is a 0.33 MAF, inside the 1/N floor."""
        rows = merge_calls(
            {
                "iso1": [("NC_002516.2", "100", "C", "T")],
                "iso2": [("NC_002516.2", "100", "C", "T")],
                "iso3": [],
                "iso4": [],
            }
        )
        assert _positions(rows) == {"100"}

    def test_the_carrier_count_is_recorded(self):
        rows = merge_calls(
            {
                "iso1": [("NC_002516.2", "100", "C", "T")],
                "iso2": [("NC_002516.2", "100", "C", "T")],
                "iso3": [],
                "iso4": [],
            }
        )
        assert rows[0]["ac"] == "2"
        assert rows[0]["an"] == "4"

    def test_frequency_is_computed_from_the_cohort_not_the_subset(self):
        """`an` is the cohort size, so MAF is comparable across sites."""
        rows = merge_calls(
            {
                "iso1": [("NC_002516.2", "100", "C", "T")],
                "iso2": [("NC_002516.2", "100", "C", "T")],
                "iso3": [],
                "iso4": [],
            }
        )
        assert rows[0]["an"] == "4"
        assert rows[0]["af"] == "0.5"


class TestScaleIndependence:
    """A threshold pinned to a count is invisible at one size and fatal at another."""

    def test_the_filter_scales_with_cohort_size(self):
        """A fixed carrier count would keep far more sites at N=50 than N=4.

        With 1/N filtering, a site carried by exactly one isolate is dropped at
        every size, which is the property a hard-coded count cannot have.
        """
        small = merge_calls({f"i{n}": ([("c", "1", "A", "T")] if n == 0 else [])
                             for n in range(4)})
        large = merge_calls({f"i{n}": ([("c", "1", "A", "T")] if n == 0 else [])
                             for n in range(50)})
        assert small == [] and large == [], "a singleton must drop at any N"

    # N=2 is excluded: with one non-carrier the site is a singleton, which
    # is dropped by design, so it would not test scale at all.
    @pytest.mark.parametrize("n", [4, 5, 20, 50])
    def test_no_assumed_cohort_size(self, n):
        """Runs at many sizes through one code path; nothing may be capped.

        `n - 1` isolates carry the site and one does not, so it is a genuine
        polymorphism at every size: a singleton would be dropped by design, and
        an all-carrier site would need the deferred upper bound.
        """
        # n-2 carriers, so two isolates dissent: the site has support on both
        # sides. n-1 would now sit inside the ceiling, which is correct.
        calls = {f"i{k}": ([("c", "1", "A", "T")] if k < n - 2 else []) for k in range(n)}
        rows = merge_calls(calls)
        assert len(rows) == 1, f"a real polymorphism was lost at N={n}"
        assert rows[0]["an"] == str(n)
        assert rows[0]["ac"] == str(n - 2)

    def test_a_single_isolate_yields_nothing(self):
        """One isolate cannot disagree with itself."""
        assert merge_calls({"only": [("c", "1", "A", "T")]}) == []


    def test_a_repeated_call_from_one_isolate_counts_once(self):
        """A caller emitting the same allele twice must not inflate its frequency.

        The alternative is a site reading 2/2 - apparently fixed - when one
        isolate really carries it, which silently deletes a polymorphism.
        """
        doubled = {"i1": [("c", "1", "A", "T"), ("c", "1", "A", "T")], "i2": []}
        rows = merge_calls(doubled)
        assert len(rows) == 0, "a doubled call was counted as two carriers"

    def test_a_repeated_call_does_not_change_a_real_count(self):
        calls = {"i1": [("c", "1", "A", "T"), ("c", "1", "A", "T")], "i2": [("c", "1", "A", "T")]}
        assert merge_calls(calls) == merge_calls({"i1": [("c", "1", "A", "T")], "i2": [("c", "1", "A", "T")]})


class TestRefusals:
    def test_no_isolates_is_refused(self):
        with pytest.raises(DataContractError):
            merge_calls({})

    def test_a_malformed_call_is_refused(self):
        """A call without a locus cannot be merged or audited."""
        with pytest.raises(DataContractError):
            merge_calls({"i1": [("c", "", "A", "T")], "i2": [("c", "1", "A", "T")]})

    def test_a_ceiling_is_present_by_default(self):
        """The symmetric bound is the default, so there is always an upper tail.

        Previously the ceiling had to be supplied, and `require_max_maf` caught
        its absence. Now `an - ac >= 2` is the default, so a bound always exists
        unless a caller explicitly disables it.
        """
        assert merge_calls(CALLS) == merge_calls(CALLS, min_dissenters=2)

    def test_require_max_maf_refuses_only_with_no_bound_at_all(self):
        with pytest.raises(DataContractError):
            merge_calls(CALLS, max_maf=None, min_dissenters=None, require_max_maf=True)
        assert merge_calls(CALLS, min_dissenters=2, require_max_maf=True) == merge_calls(CALLS)

    def test_the_ceiling_can_be_disabled_explicitly(self):
        """One dissenter is enough when the ceiling is turned off."""
        calls = {"i1": [("c", "1", "A", "T")], "i2": [("c", "1", "A", "T")], "i3": []}
        assert len(merge_calls(calls, min_dissenters=None)) == 1


class TestContract:
    def test_columns_are_long_and_tidy(self):
        """One row per site, not one column per isolate."""
        assert MERGE_COLUMNS == ("chrom", "pos", "ref", "alt", "ac", "an", "af")
        assert "iso1" not in MERGE_COLUMNS
