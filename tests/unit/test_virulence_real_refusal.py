"""Stage 5's mode branch: the guarantees that survive having a REAL caller.

This file used to assert that stage 5 *refused* REAL mode. That was correct
while there was no caller (7e28840) and the stage read a table nothing produced.
The raw-blastn caller now exists, so the refusal is gone and four assertions
here were retired. Their reasoning is kept in `git log`; the coverage they
provided is now:

* `test_virulence_real_execution.py::TestRealNoLongerRefuses` - REAL reaches the
  database lookup and fails with a provisioning error naming the path, rather
  than a missing-table or missing-implementation error.
* `test_virulence_blast_adapter.py` - the caller itself.

The two guarantees below are unaffected by that change and are kept: TEST still
reads the committed fixture, and the stage still branches on mode. The second
matters because the branch is what routes REAL to the caller; losing it would
put every mode back on the fixture path, which looks like a successful TEST run
and a fabricated REAL one.
"""

from __future__ import annotations

import pytest

from papipeline.models import RunMode
from papipeline.stages import virulence as stage_virulence


def _manifest():
    from papipeline.manifest import SampleManifest
    from papipeline.models import Sample

    return SampleManifest(samples=[Sample(sample_id="PDT_A")])


class TestTestModeIsUnchanged:
    def test_test_mode_still_reads_the_fixture(self, config, tmp_path):
        """The REAL branch is additive; the fixture path is untouched.

        The header is the stage's real `REQUIRED` columns. An earlier version of
        this test invented its own (`gene, category, coverage_pct`) and then
        "passed" a fixture the loader could never have read - which would have
        hidden a genuine regression behind a plausible-looking failure.
        """
        target = tmp_path / "virulence"
        target.mkdir(parents=True)
        (target / "virulence_factors.tsv").write_text(
            "sample_id\tvirulence_factor\tgene\tdatabase\tdatabase_version\n"
            "PDT_A\tvfh\tvfh\tSYNTHETIC_VFDB\tv0-synthetic\n",
            encoding="utf-8",
        )
        result = stage_virulence.run(config, _manifest(), RunMode.TEST, tmp_path)
        assert set(result) == {"PDT_A"}
        assert result["PDT_A"], "the fixture row was not returned"


class TestTheModeBranchSurvives:
    def test_virulence_branches_on_mode(self):
        from pathlib import Path

        source = (
            Path(__file__).resolve().parents[2]
            / "papipeline" / "stages" / "virulence.py"
        ).read_text(encoding="utf-8")
        assert "if mode is not RunMode.TEST:" in source, (
            "virulence.py no longer branches on mode, so every mode would take "
            "the fixture path - a passing TEST run and a fabricated REAL one"
        )

    def test_the_real_branch_does_not_read_a_table(self):
        """REAL must reach the caller, not the fixture reader.

        Worth keeping now that both paths exist: a stale
        `virulence_factors.tsv` from an earlier run must not be able to satisfy
        a REAL run, which would report a real screen's results from data that
        was never screened.
        """
        from pathlib import Path

        source = (
            Path(__file__).resolve().parents[2]
            / "papipeline" / "stages" / "virulence.py"
        ).read_text(encoding="utf-8")
        real_body = source.split("def _run_real(")[1]
        assert "load_virulence_factors" not in real_body, (
            "the REAL branch reads the fixture table; it must screen with blastn"
        )