"""The `cohort_variants` stage entry point.

`stages/cohort_variants.py` had `merge_calls` — a set operation across isolates,
validated at N=4, 12 and 34 — and no caller. This tests the seam that makes it
reachable: `cohort_variants.run()`.

What is pinned here is the wiring, not the merge. The filter itself is proven in
`test_cohort_variants_merge.py` and, against committed fixtures, in
`test_variants_test_fixtures.py`. What was missing was the question this file
answers: **what does the merge receive, and what does the stage hand back?**

**The cohort denominator is the manifest, not the callers.** This is the one
property here that could silently inflate a result. `merge_calls` computes
`an = len(calls_by_isolate)`. If a sample with no calls were absent from the
mapping, a 20-isolate cohort in which 5 produced calls would report `an=5`, and
every allele frequency would be four times too large — with no error, and a
table that looks perfectly well formed. `run()` therefore seeds every manifest
sample, and this is asserted at the seam where it is preventable.

**The unmeasurable ceiling stays unset.** `max_maf` is deliberately not
defaulted. It has now been measured at three cohort sizes and could not be
resolved: the fixed-difference mode *shrank* as N grew, which a genuinely
species-fixed set cannot do, and the upper tail piles up at 94% rather than
100% (docs/design/cohort-variant-merge.md 5c). Fitting a ceiling to that bump
would fit an assembly artefact. So the symmetric bound applies and no ceiling
is passed — and if a caller ever sets `require_max_maf`, the merge refuses
rather than proceeding on a bound nobody measured.

**Every isolate is a member, so a sample that cannot be aligned is not a
finding.** An isolate with no assembly is an ordinary cohort member carrying no
variant. That is the same distinction the per-isolate stage draws, and it is
what keeps the two stages from disagreeing about who is in the cohort.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.manifest import SampleManifest
from papipeline.models import RunMode, Sample
from papipeline.stages import cohort_variants as stage_cohort

COHORT = tuple(f"TEST_PA_00{i}" for i in range(1, 6))


def _manifest(*sample_ids: str) -> SampleManifest:
    return SampleManifest(samples=[Sample(sample_id=s) for s in sample_ids])


def _calls(config, tmp_path, sample_ids=COHORT):
    """Per-isolate calls, produced by the stage this one consumes."""
    from papipeline.stages import variants as stage_variants

    return stage_variants.run(
        config, _manifest(*sample_ids), RunMode.TEST,
        data_root=Path("test_data"), workdir=tmp_path,
    )


class TestTheStageHasAnEntryPoint:
    def test_run_exists(self):
        assert callable(stage_cohort.run)


class TestItMergesWhatThePerIsolateStageProduced:
    def test_it_returns_rows(self, config, tmp_path):
        rows = stage_cohort.run(
            config, _manifest(*COHORT), _calls(config, tmp_path)
        )
        assert rows, "the committed cohort fixtures merged to nothing"

    def test_the_rows_carry_the_declared_merge_contract(self, config, tmp_path):
        rows = stage_cohort.run(config, _manifest(*COHORT), _calls(config, tmp_path))
        for row in rows:
            assert set(row) == set(stage_cohort.MERGE_COLUMNS)

    def test_the_merge_output_is_not_a_per_sample_table(self, config, tmp_path):
        """No `sample_id` column: this is one row per site, cohort-wide.

        Asserted because the temptation to make the merge's output look like
        every other stage table is exactly what would erase the distinction
        between "one isolate's calls" and "the cohort's polymorphisms".
        """
        rows = stage_cohort.run(config, _manifest(*COHORT), _calls(config, tmp_path))
        assert rows
        assert "sample_id" not in stage_cohort.MERGE_COLUMNS


class TestTheDenominatorIsTheCohort:
    def test_an_isolate_with_no_calls_is_still_counted(self, config, tmp_path):
        """`an` must equal the manifest, not the number of calling isolates.

        The fixtures cover five isolates. Adding one that calls nothing must
        leave every allele frequency *half* what it was — which is the whole
        point of including it, and the property most easily lost.
        """
        # The sixth isolate is in the manifest: it is a cohort member that
        # happens to carry no call, not an outsider.
        five = _manifest(*COHORT)
        six = _manifest(*(COHORT + ("TEST_PA_006",)))
        before = stage_cohort.run(
            config, five, _calls(config, tmp_path)
        )
        after = stage_cohort.run(
            config, six,
            _calls(config, tmp_path, COHORT + ("TEST_PA_006",)),
        )
        assert before and after
        assert {r["an"] for r in before} == {"5"}
        assert {r["an"] for r in after} == {"6"}, (
            "an isolate that produced no calls was dropped from the denominator"
        )
        # And the frequencies really did fall, rather than only the label moving.
        # Contig-agnostic: the committed TEST fixtures use their own contig name,
        # and a test that hard-coded the PAO1 accession would pass for the wrong
        # reason if the fixtures were ever regenerated against another reference.
        first = before[0]
        same_site = next(
            r for r in after
            if (r["chrom"], r["pos"], r["ref"], r["alt"])
            == (first["chrom"], first["pos"], first["ref"], first["alt"])
        )
        assert float(same_site["af"]) == pytest.approx(
            float(first["af"]) * 5 / 6, rel=1e-4
        )


class TestNoUnmeasuredCeiling:
    def test_it_does_not_invent_a_max_maf(self, config, tmp_path, monkeypatch):
        """A ceiling fitted to the 94% bump would be fitting an artefact.

        Measured at N=4, 12 and 34 and still unresolvable, so the merge runs
        with the symmetric bound only. Asserted by watching what the stage
        passes: if it ever supplies a ceiling, this fails and the decision has
        to be re-made deliberately rather than inherited from a default.
        """
        seen = {}
        real_merge = stage_cohort.merge_calls

        def spy(calls_by_isolate, **kwargs):
            seen.update(kwargs)
            return real_merge(calls_by_isolate, **kwargs)

        monkeypatch.setattr(stage_cohort, "merge_calls", spy)
        stage_cohort.run(config, _manifest(*COHORT), _calls(config, tmp_path))
        assert seen.get("max_maf") is None
        assert not seen.get("require_max_maf")


class TestRefusals:
    def test_a_cohort_of_zero_is_refused(self, config):
        """Nothing to compare with itself, so there is no merge to perform."""
        with pytest.raises(Exception) as excinfo:
            stage_cohort.run(config, _manifest(), {})
        assert "cohort" in str(excinfo.value).lower()

    def test_a_manifest_of_one_who_called_nothing_merges_to_nothing(
        self, config
    ):
        """Zero rows is the *right* answer here, and must not be an error.

        A cohort that exists but produced no calls is a finding - every isolate
        was screened and nothing was found. Raising here would make the stage
        refuse the one result it is most likely to produce on a real cohort,
        where 82% of assemblies are unusable.
        """
        assert stage_cohort.run(config, _manifest("TEST_PA_001"), {}) == []

    def test_calls_from_an_isolate_outside_the_manifest_are_refused(
        self, config, tmp_path
    ):
        """Hard rule 5: a membership mismatch fails loudly, never silently.

        The manifest decides who is in the study. A call from an isolate it does
        not name would otherwise be counted into `an`, inflating the denominator
        with a sample the cohort does not contain — and every allele frequency
        with it.
        """
        calls = _calls(config, tmp_path)
        calls["TEST_PA_999"] = [
            {"chrom": "TESTCHR", "pos": "500", "ref": "A", "alt": "T"}
        ]
        with pytest.raises(Exception) as excinfo:
            stage_cohort.run(config, _manifest(*COHORT), calls)
        message = str(excinfo.value)
        assert "TEST_PA_999" in message
        assert "manifest" in message.lower()
