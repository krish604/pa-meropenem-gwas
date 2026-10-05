"""Stage 7's mode branch: the guarantees that survive having a REAL caller.

This file used to assert that stage 7 *refused* REAL mode. That was correct
while there was no caller (a50bc0f) and the stage read a table nothing produced.
The assembly-based caller now exists, so the refusal is gone and three
assertions here were retired. Their reasoning is kept in `git log`; the coverage
they provided is now:

* `test_sv_real_execution.py::TestRealNoLongerRefuses` - REAL reaches the
  reference and tool lookup and fails with a provisioning error naming what is
  missing, rather than a missing-implementation error.
* `test_sv_real_execution.py::TestRealIsCandidateOnly` - the property that
  matters more than the refusal ever did: no REAL call can be `confirmed`.
* `test_sv_adapter.py` - the caller itself.

The two guarantees below are unaffected and are kept: TEST still reads the
committed fixture, and the stage still branches on mode. The second matters
because the branch is what routes REAL to the caller; losing it would put every
mode on the fixture path, which looks like a successful TEST run and a
fabricated REAL one.
"""

from __future__ import annotations

from papipeline.models import RunMode
from papipeline.stages import sv as stage_sv


class TestTestModeIsUnchanged:
    def test_test_mode_still_reads_the_fixture(self, config, tmp_path):
        """The REAL branch is additive; the fixture path is untouched.

        The header is the stage's real `REQUIRED` columns. An earlier version of
        this test invented its own and then "passed" a fixture the loader could
        never have read - which would have hidden a genuine regression behind a
        plausible-looking failure.
        """
        from papipeline.manifest import SampleManifest
        from papipeline.models import Sample

        target = tmp_path / "structural_variants"
        target.mkdir(parents=True)
        (target / "structural_variants.tsv").write_text(
            "sample_id\tvariant_id\tvariant_type\tcall_status\n"
            "PDT_A\tv1\tdeletion\tconfirmed\n",
            encoding="utf-8",
        )
        manifest = SampleManifest(samples=[Sample(sample_id="PDT_A")])
        result = stage_sv.run(config, manifest, RunMode.TEST, tmp_path)
        assert result["PDT_A"], "the fixture row was not returned"
        assert result["PDT_A"][0].call_status.value == "confirmed", (
            "TEST fixtures carry confirmed calls; the candidate-only rule is a "
            "property of the REAL caller, not of the stage's parser"
        )


class TestTheModeBranchSurvives:
    def test_sv_branches_on_mode(self):
        from pathlib import Path

        source = (
            Path(__file__).resolve().parents[2]
            / "papipeline" / "stages" / "sv.py"
        ).read_text(encoding="utf-8")
        assert "if mode is not RunMode.TEST:" in source, (
            "sv.py no longer branches on mode, so every mode would take the "
            "fixture path - a passing TEST run and a fabricated REAL one"
        )

    def test_the_real_branch_does_not_read_the_fixture_table(self):
        """REAL must reach nucmer, not the fixture reader.

        Worth keeping now that both paths exist: a stale
        `structural_variants.tsv` from an earlier run must not be able to
        satisfy a REAL run, which would report a real screen's results from
        data that was never screened.
        """
        from pathlib import Path

        source = (
            Path(__file__).resolve().parents[2]
            / "papipeline" / "stages" / "sv.py"
        ).read_text(encoding="utf-8")
        real_body = source.split("def _run_real(")[1]
        assert "load_structural_variants" not in real_body, (
            "the REAL branch reads the fixture table; it must call nucmer"
        )