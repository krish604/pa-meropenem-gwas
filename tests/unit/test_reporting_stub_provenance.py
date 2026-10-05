"""A STUB report must never claim REAL provenance.

`build_markdown` chose its banner with a **two-way** test - `mode is RunMode.TEST`
and everything else took the REAL branch - and `RunMode` has three members. So
STUB, whose entire purpose is to fabricate stage outputs with no genome, tool,
parser or fixture involved (`papipeline/stub.py`), rendered a report titled
`Pseudomonas aeruginosa Imipenem AMR Report` and bannered `REAL DATASET
REPORT`. The STUB report was indistinguishable from a real one in the two fields
a reader looks at first, and in the filename the operator finds on disk.

That matters more for STUB than for any other mode, and not only because the
content is fake. TEST's numbers come from committed fixtures and are therefore
*stable and reproducible* - "wrong, but reproducibly wrong", which a reader can
eventually discover. STUB's are not derived from anything at all: they are
whatever the shape of an empty table looks like. There is no input to re-run and
compare against, so a STUB number can never be shown to be fabricated after the
fact. The banner is the only thing standing between it and a reader.

STUB already had the right words. `papipeline.stub.STUB_BANNER` is stamped into
the report `stub.fabricate_report` writes - the report the *workflow* produces.
`stages.reporting.write_report` is a second, independent path to the same
artefact, and it is the one that reached for `REAL_BANNER`. These tests pin that
both paths say the same thing.

What this file pins:

* STUB's banner and title are its own, and are not REAL's;
* a STUB report on disk claims REAL nowhere - not in the Markdown, not in the
  HTML, not in the filename;
* TEST and REAL keep exactly the labels they had. A fix that silences the
  fall-through by making every non-REAL mode look synthetic would pass the first
  group and be wrong;
* the three modes are pairwise distinguishable, so "which run produced this" is
  answerable from the report alone;
* the banner choice is **fail-closed**: a mode the table does not know inherits
  no provenance at all rather than defaulting to REAL. The defect was a silent
  fall-through onto the most dangerous label, so the guard against its return is
  that the unreachable branch claims nothing.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.models import RunMode
from papipeline.stub import STUB_BANNER
from papipeline.stages.reporting import REAL_BANNER, TEST_BANNER, write_report

REAL_TITLE = "# Pseudomonas aeruginosa Imipenem AMR Report"


def _context(config, mode: RunMode):
    """A STUB-legal context: nothing declared, which is what STUB really is."""
    from papipeline.stages.reporting import ReportContext

    return ReportContext(
        mode=mode,
        config=config,
        generated_at="2026-01-01T00:00:00Z",
        antibiotic="imipenem",
        n_samples=20,
    )


def _markdown(config, mode: RunMode) -> str:
    from papipeline.stages.reporting import build_markdown

    return build_markdown(_context(config, mode))


# --------------------------------------------------------------------------
# The defect
# --------------------------------------------------------------------------


class TestAStubReportClaimsNoRealProvenance:
    def test_its_banner_is_not_the_real_banner(self, config):
        """The exact fall-through: `is RunMode.TEST` sent STUB to REAL_BANNER."""
        text = _markdown(config, RunMode.STUB)
        assert REAL_BANNER not in text, (
            "a STUB report is bannered REAL DATASET REPORT; STUB fabricates "
            "every output and has no provenance to report"
        )
        assert "REAL DATASET REPORT" not in text

    def test_its_title_is_not_the_real_title(self, config):
        assert not _markdown(config, RunMode.STUB).startswith(REAL_TITLE), (
            "a STUB report is titled as the real report, so the filename and "
            "the first line of the document agree about nothing"
        )

    def test_it_carries_the_stub_banner_the_workflow_already_uses(self, config):
        """One wording for both paths to the artefact.

        `stub.fabricate_report` stamps this into the report the Snakemake rule
        produces. Reporting is the other producer, and it must not invent a
        second, weaker phrasing of the same warning.
        """
        assert STUB_BANNER in _markdown(config, RunMode.STUB)

    def test_it_names_its_own_mode_in_the_run_metadata(self, config):
        """`run_mode` was already honest. Pin it, so a fix cannot regress it."""
        assert "| run_mode | STUB |" in _markdown(config, RunMode.STUB)


class TestAStubReportOnDiskClaimsNoRealProvenance:
    def test_neither_the_markdown_nor_the_html_claims_real(
        self, config, tmp_path: Path
    ):
        """The HTML is the half that is read first, and it was just as wrong.

        `markdown_to_html` copies the banner into a red-bordered `div`, so the
        mislabelling survived the render. Checking only the Markdown would have
        called this fixed while the browser still said REAL DATASET REPORT.
        """
        written = write_report(
            _context(config, RunMode.STUB), tmp_path, write_html=True
        )
        assert written["markdown"].exists()
        assert written["html"].exists()
        for artefact in ("markdown", "html"):
            text = written[artefact].read_text(encoding="utf-8")
            assert "REAL DATASET REPORT" not in text, (
                f"the {artefact} artefact claims REAL provenance"
            )
            assert STUB_BANNER in text, (
                f"the {artefact} artefact does not say it is fabricated"
            )

    def test_it_is_not_written_under_the_real_filename(self, config, tmp_path: Path):
        """A STUB report must not sit on disk named as the real analysis.

        `write_report` picked the real name with `is RunMode.REAL`, which STUB
        does not satisfy - so the name was already right. Pinned because the
        filename is how an operator tells the two apart once the file has been
        mailed or archived, and because it is the field most likely to be
        "tidied" into a mode test that reintroduces the fall-through.
        """
        written = write_report(
            _context(config, RunMode.STUB), tmp_path, write_html=False
        )
        assert written["markdown"].name != (
            "Pseudomonas_aeruginosa_Imipenem_AMR_Report.md"
        )
        assert written["markdown"].name == "test_pipeline_report.md"


# --------------------------------------------------------------------------
# The other two modes must be untouched
# --------------------------------------------------------------------------


class TestTestAndRealKeepTheirOwnLabels:
    def test_test_keeps_its_synthetic_banner(self, config):
        text = _markdown(config, RunMode.TEST)
        assert TEST_BANNER in text
        assert REAL_BANNER not in text

    def test_real_keeps_its_real_banner_and_title(self, config):
        text = _markdown(config, RunMode.REAL)
        assert REAL_BANNER in text
        assert text.startswith(REAL_TITLE)

    def test_the_three_modes_are_pairwise_distinguishable(self, config):
        """"Which run produced this" must be answerable from the report alone.

        Two modes sharing a title and a banner is precisely what let a STUB
        report pass for a real one, so this is the invariant that actually
        matters rather than any individual string.
        """
        rendered = {m: _markdown(config, m) for m in RunMode}
        for left, right in (("STUB", "TEST"), ("STUB", "REAL"), ("TEST", "REAL")):
            assert rendered[left] != rendered[right], (
                f"{left} and {right} render identically, so provenance cannot be "
                f"read off the report"
            )


# --------------------------------------------------------------------------
# The banner choice must fail closed
# --------------------------------------------------------------------------


class TestTheBannerChoiceFailsClosed:
    def test_a_mode_the_table_does_not_know_claims_no_provenance(self):
        """The defect was a silent fall-through onto the most dangerous label.

        `RunMode` is a closed enum, so no such mode exists today and this test
        cannot fail on the shipped code. It is here because the fix replaces a
        two-way `if` with a table lookup, and the failure mode that must not
        return is a lookup whose miss lands on REAL_BANNER. Under-claiming
        provenance is recoverable; over-claiming it is the whole bug.

        Asserted on `banner_for` rather than through `build_markdown`: the
        renderer requires a real `RunMode` (it reads `mode.value` for the run
        metadata), so passing it a stranger proves only that it validates its
        input - not how the table resolves a miss.
        """
        from papipeline.stages.reporting import banner_for

        assert banner_for(RunMode.STUB) == STUB_BANNER
        assert banner_for("BENCHMARK") == STUB_BANNER  # type: ignore[arg-type]
        assert banner_for(None) == STUB_BANNER  # type: ignore[arg-type]

    def test_an_unknown_mode_gets_no_title_claiming_real_provenance(self):
        from papipeline.stages.reporting import title_for

        assert title_for("BENCHMARK") == title_for(RunMode.STUB)  # type: ignore[arg-type]
        assert title_for("BENCHMARK") != title_for(RunMode.REAL)  # type: ignore[arg-type]

    @pytest.mark.parametrize("mode", list(RunMode))
    def test_every_mode_in_the_enum_has_a_banner(self, mode: RunMode):
        """The table covers the enum exhaustively, so no member can miss.

        A `RunMode` member with no entry is how the next mode silently acquires
        the wrong provenance; this asserts the table and the enum agree.
        """
        from papipeline.stages.reporting import _BANNERS, _TITLES

        assert mode in _BANNERS
        assert mode in _TITLES