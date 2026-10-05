"""TEST mode must exercise the real variant code, not a stand-in.

`variants` and `cohort_variants` have real implementations, and for a while the
taxonomy still called them unbuilt - so they refused outside STUB and TEST read
a labelled stand-in. That is the stand-in refusal doing its job, but it was
based on a fact the code no longer reflected: the parser and the merge exist.

So TEST gets real fixtures. This is not routing around the stand-in mechanism -
it is using it correctly. Stand-ins remain right for `recombination` and
`similarity`, which genuinely have no code, and `papipeline.standins` continues
to refuse anything in `test_data/standins/unbuilt_stages/`. The check below
asserts both halves, because "the stand-in is no longer used here" and "the
stand-in guard still works elsewhere" are independent properties and only the
second one is load-bearing for the unbuilt stages.

The fixtures are tiny and shaped like real tool output - a VCF with a proper
header, multiallelic ALT, `DP4` as a four-tuple - because the point is to drive
the real parser. A fixture that avoided the awkward shapes would test a parser
that never sees them in production.

The cohort fixtures exercise both tails of the symmetric filter: `ac >= 2` and
`an - ac >= 2`. A test that only used comfortably-variable sites would pass
against a filter that ignored the ceiling entirely, which is a mistake this repo
has made before.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.stages.cohort_variants import merge_calls
from papipeline.stages.variants import PER_ISOLATE_COLUMNS, parse_vcf

FIXTURES = Path("test_data/variants")


def _read(name: str) -> str:
    return (FIXTURES / name).read_text(encoding="utf-8")


class TestFixturesExistAndAreReal:
    def test_the_directory_is_committed(self):
        assert FIXTURES.is_dir(), "test_data/variants/ is missing"
        assert any(FIXTURES.glob("*.vcf")), "no per-isolate VCF fixtures"

    def test_they_are_not_stand_ins(self):
        """A stand-in here would make the whole exercise circular."""
        for path in FIXTURES.glob("*"):
            assert "standins" not in path.parts, f"{path} lives under standins/"
            assert "STAND-IN" not in path.read_text(encoding="utf-8")[:400], (
                f"{path.name} carries a stand-in banner"
            )

    def test_the_standin_guard_still_protects_the_unbuilt_stages(self):
        """The mechanism is untouched, just no longer used by these two."""
        from papipeline.standins import (
            STANDIN_BANNER,
            is_standin_table,
            load_stage_table,
        )

        root = Path("test_data/standins/unbuilt_stages")
        # `similarity` left this set in 2223332 and `recombination` left it when
        # `dag-resolve` dispatched stage 8, so the tuple is empty: every stage
        # in `STAGE_ORDER` computes its own output. Asserted rather than assumed,
        # because a tuple that quietly grew again is how a fabricated table gets
        # back into the DAG.
        from papipeline.standins import UNBUILT_WITH_STANDIN

        assert UNBUILT_WITH_STANDIN == (), (
            "no stage is unbuilt, so none may be stood in for. The files stay on "
            "disk and are unreferenced; see papipeline/standins.py."
        )
        # The mechanism still works, on the fixture files that remain: the banner
        # is still detected and the loader still refuses to hand one back as a
        # real result.
        assert is_standin_table(root / "recombination.tsv")
        with pytest.raises(Exception):
            load_stage_table("recombination", root / "recombination.tsv")
        assert STANDIN_BANNER in (root / "recombination.tsv").read_text(
            encoding="utf-8"
        )


class TestRealParserOnRealShapedVcf:
    def test_it_parses_the_committed_fixture(self):
        rows = parse_vcf(_read("TEST_PA_001.vcf"), sample_id="TEST_PA_001")
        assert rows, "the fixture produced no rows"
        for row in rows:
            assert set(row) == set(PER_ISOLATE_COLUMNS)

    def test_a_multiallelic_site_becomes_one_row_per_alt(self):
        """Real VCFs carry these; a fixture that avoided them proves nothing."""
        rows = parse_vcf(_read("TEST_PA_001.vcf"), sample_id="TEST_PA_001")
        by_pos = {}
        for r in rows:
            by_pos.setdefault(r["pos"], []).append(r["alt"])
        assert "." not in [a for v in by_pos.values() for a in v]
        # The fixture has a `T,G` site; it must yield exactly two rows, and the
        # second ALT must still be present. A parser that keeps only the first
        # allele passes any looser check, which is how one survived mutation.
        multi = by_pos["602"]
        assert multi == ["T", "G"], f"multiallelic site collapsed to {multi}"

    def test_dp4_stays_a_four_tuple(self):
        rows = parse_vcf(_read("TEST_PA_001.vcf"), sample_id="TEST_PA_001")
        with_dp4 = [r for r in rows if r["DP4"]]
        assert with_dp4, "fixture carries no DP4, so the tuple shape is untested"
        for row in with_dp4:
            assert row["DP4"].count(",") == 3


class TestMergeOnCommittedFixtures:
    def _calls(self):
        calls = {}
        for path in sorted(FIXTURES.glob("TEST_PA_*.vcf")):
            calls[path.stem] = [
                (r["chrom"], r["pos"], r["ref"], r["alt"])
                for r in parse_vcf(path.read_text(encoding="utf-8"), sample_id=path.stem)
            ]
        return calls

    def test_the_fixtures_merge_to_a_non_empty_table(self):
        """A merge that returns nothing would satisfy 'no crash' and be useless."""
        rows = merge_calls(self._calls())
        assert rows, "the committed cohort fixtures merged to nothing"

    def test_the_fixtures_contain_a_singleton_so_the_floor_is_exercised(self):
        """Without this the whole file passes against a filter with no floor."""
        calls = self._calls()
        n = len(calls)
        counts = {
            c: sum(1 for per in calls.values() if c in per)
            for iso, per in calls.items() for c in set(per)
        }
        assert any(k == 1 for k in counts.values()), (
            f"no singleton site among the fixtures: {sorted(set(counts.values()))}"
        )
        assert any(k == n for k in counts.values()), "no fixed difference to test the ceiling"

    def test_every_retained_site_has_support_on_both_tails(self):
        """The symmetric filter, on real fixtures rather than inline dicts."""
        for row in merge_calls(self._calls()):
            ac, an = int(row["ac"]), int(row["an"])
            assert ac >= 2, f"site kept with {ac} carrier(s)"
            assert an - ac >= 2, f"site kept with {an - ac} dissenter(s)"

    def test_the_cohort_size_is_the_number_of_fixtures(self):
        """`an` must be the cohort, not the number of surviving sites."""
        calls = self._calls()
        for row in merge_calls(calls):
            assert row["an"] == str(len(calls))

    def test_a_fixed_difference_in_the_fixtures_is_dropped(self):
        """The fixture must actually contain a site to exclude.

        Otherwise this test would pass against a filter with no ceiling at all,
        which is the exact failure the symmetric rule was added to prevent.
        """
        calls = self._calls()
        n = len(calls)
        every = [c for iso, per in calls.items() for c in per]
        shared = {c for c in every if sum(1 for per in calls.values() if c in per) == n}
        assert shared, "no fixed difference in the fixtures, so nothing tests the ceiling"
        kept = {(r["chrom"], r["pos"], r["ref"], r["alt"]) for r in merge_calls(calls)}
        assert not (shared & kept), "a fixed difference survived the merge"
