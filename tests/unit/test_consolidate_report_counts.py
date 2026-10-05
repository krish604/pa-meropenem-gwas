"""The cohort report must not state a count it never measured.

`scripts/consolidate_results.py` used to write a literal census into its output:

    - Assemblies present in `data/`: **835**
    - Not analysable: **835 - n** - a census of all 835 found 233 zero-byte, 51
      truncated/partial and 426 otherwise corrupt files, leaving 93 intact ...

That prose was a *pilot* result - `run_pilot100.py:216` computes a real census
through `cohort.assess_cohort` - quoted as though it were a census of `data/`.
No pipeline run produces one: the validation stage records QC for manifest
samples only, so there was no run-specific figure behind those numbers. The
report now says so rather than stating a figure it cannot source.

Scope note: this is a **source guard**, not a behavioural test. Driving
`_write_report` end to end needs the real shapes of seven nested inputs
(`downstream['annotation']['total_records']`, `downstream['mlst']['n_typed']`,
and more, several interpolated with `{:,}`), and reconstructing them is a
larger job than it looks - see the handoff. This guard pins the thing that
actually regressed: someone re-adding a hardcoded census.
"""

from __future__ import annotations

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "consolidate_results.py"


def _report_source() -> str:
    return SCRIPT.read_text(encoding="utf-8")


class TestNoUnmeasuredCensusInTheReport:
    def test_it_states_no_literal_assembly_total(self):
        """The cohort size must never be a constant in the report."""
        source = _report_source()
        body = source[source.index("def _write_report"):]
        assert not re.search(r"\b8\d\d\b", body), (
            "the report body contains a literal cohort size; it must report the "
            "cohort it was given"
        )

    def test_it_carries_no_census_figures(self):
        """Guard the *emitted* lines, so the comment explaining the history passes.

        The rationale above the block names the figures it removed; a comment is
        not report text, so the check has to look at the A("...") calls.
        """
        body = _report_source()
        body = body[body.index("def _write_report"):]
        emitted = "\n".join(
            match.group(1)
            for match in re.finditer(r'A\((f?)"([^"]*)"', body)
        )
        for figure in ("zero-byte", "truncated/partial", "otherwise corrupt"):
            assert figure not in emitted, (
                f"the report still states a census figure ({figure!r}) that no "
                "run measured"
            )

    def test_it_says_plainly_that_no_census_was_performed(self):
        body = _report_source()
        body = body[body.index("def _write_report"):]
        assert re.search(r"No census", body), (
            "the report must tell the reader the census was not performed, "
            "rather than dropping the section and leaving them to assume the "
            "directory was examined"
        )

    def test_it_names_where_a_real_census_can_be_obtained(self):
        body = _report_source()
        body = body[body.index("def _write_report"):]
        assert "assess_cohort" in body, (
            "the report should point at the code that does compute a census, "
            "rather than leaving the reader with nowhere to go"
        )

    def test_the_cohort_it_does_report_is_computed(self):
        body = _report_source()
        body = body[body.index("def _write_report"):]
        assert "n = len(rows)" in body, (
            "the reported cohort must come from the rows it was given"
        )
