"""A run-written report states stage 11's phenotype provenance.

**What this is.** `papipeline/stages/reporting.py:172` `phenotype_report` is the
PRODUCTION caller of `papipeline.stages/phenotype.provenance_report`, and it is
reached from `build_markdown` (`reporting.py:364`) on the reporting path of every
run. It reads exactly one thing: `ReportContext.phenotype_calls`
(`reporting.py:138`). That field defaulted to `None`, and `run.py` did not pass
it - so `phenotype_report` returned `None` on every run, and the provenance
section was never rendered.

The module documents the consequence of `None` in both directions, and both were
being taken. In REAL, `build_markdown`'s `elif context.mode is RunMode.REAL`
branch (`reporting.py:368-379`) then prints that "No stage 11 phenotype calls
were handed to this report, so it states nothing about the cohort's
susceptibility calls". **That was false of the run.** Stage 11 had loaded the
calls (`run.py:1289`), the join had kept them (`run.py:1309`), and
`_execute_stage` was still holding them in `phenotype_calls` - it even handed
the same list to stage 12a and to `stage_gwas.run`. The report said the run held
nothing it was in fact holding. In TEST and STUB nothing was printed at all,
because the absence branch is REAL-only, so the omission was invisible there
rather than absent.

**Why this is driven through `run_pipeline` and not through the context.** The
defect is a wiring defect. `_build_report_context` could assemble a
`ReportContext` with `phenotype_calls` populated perfectly and the run would
still write a report without the section, and conversely the field could be
dropped from the context construction while every test of
`reporting.phenotype_report` in isolation kept passing - those tests hand-build
a context, so they cannot see the orchestrator forget the argument. Only a
report that a run actually wrote distinguishes the two.

**What is asserted.** The committed TEST phenotype fixture
(`test_data/phenotype/imipenem_phenotype.tsv`, 20 rows) carries no
`ast_standard` column, so every row lacks one and
`reporting.NO_STANDARD_RECORDED` must appear - that is the strongest form of the
provenance wording, because it is the part a reader relies on and the part that
cannot be produced by a section that merely exists. The absence wording is
asserted *absent*, because that string is the specific false claim this fixes
and asserting only its opposite would let a future edit reintroduce it.

**Anti-vacuity.** `test_the_fixture_really_lacks_a_standard_column` reads the
fixture and fails if it ever gains an `ast_standard` column, which would make the
`NO_STANDARD_RECORDED` assertion vacuous. The run is also required to have
written the report file this file then reads; a report that is absent fails
rather than skips.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from papipeline.config.loader import RESULTS_ROOT_ENV, load_config
from papipeline.errors import PipelineError
from papipeline.run import run_pipeline
from papipeline.stages import reporting

REPO = Path(__file__).resolve().parents[2]
PHENOTYPE_FIXTURE = REPO / "test_data" / "phenotype" / "imipenem_phenotype.tsv"

#: The section heading `phenotype_provenance_section` builds. Asserted on its own
#: so a report carrying the section for an unrelated reason cannot pass.
PROVENANCE_SENTENCE = "enter the"


@pytest.fixture()
def project(tmp_path, monkeypatch):
    """A pipeline whose every root is a distinct path under ``tmp_path``.

    The same shape as `tests/unit/test_regulators_input_output_paths.py`: config
    is copied so the pipeline root moves and every accessor becomes a temporary
    path, ``test_data`` is symlinked because `tool_output_root(TEST)` has no
    redirect hook, and `PIPELINE_RESULTS_ROOT` moves the run's own output tree.
    Reports follow the results redirect (`config/loader.py:857`), so the report
    this file reads is under ``tmp_path`` and not in the repository.
    """
    root = tmp_path / "project"
    root.mkdir()
    shutil.copytree(REPO / "config", root / "config")
    (root / "test_data").symlink_to(REPO / "test_data", target_is_directory=True)
    monkeypatch.setenv(RESULTS_ROOT_ENV, str(tmp_path / "written"))
    return load_config(root / "config" / "science.yaml")


@pytest.fixture(scope="module")
def report_text(tmp_path_factory):
    """The Markdown a whole TEST run wrote. Module-scoped: one run, many reads."""
    import os

    written = tmp_path_factory.mktemp("written")
    root = tmp_path_factory.mktemp("project")
    shutil.copytree(REPO / "config", root / "config")
    (root / "test_data").symlink_to(REPO / "test_data", target_is_directory=True)
    previous = os.environ.get(RESULTS_ROOT_ENV)
    os.environ[RESULTS_ROOT_ENV] = str(written)
    try:
        config = load_config(root / "config" / "science.yaml")
        failure = ""
        try:
            run_pipeline(config=config, mode="TEST", write_html=False)
        except PipelineError as exc:
            failure = f"{type(exc).__name__}: {exc}"
        # Reported through an assertion below rather than raised here, so the
        # failure names the run's own error instead of a missing file.
        assert not failure, f"the TEST run did not complete: {failure}"
        report = written / "test" / "reports" / "test_pipeline_report.md"
        assert report.is_file(), (
            f"no report at {report}. `reports_root` follows the results redirect "
            "(config/loader.py:857), so this is where a TEST run writes it."
        )
        return report.read_text(encoding="utf-8")
    finally:
        if previous is None:
            os.environ.pop(RESULTS_ROOT_ENV, None)
        else:
            os.environ[RESULTS_ROOT_ENV] = previous


class TestTheFixtureMakesTheAssertionNonVacuous:
    """The premise. A fixture that gained a standard column voids the test."""

    def test_the_fixture_really_lacks_a_standard_column(self):
        header = [
            line
            for line in PHENOTYPE_FIXTURE.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")
        ][0].split("\t")
        assert "ast_standard" not in header, (
            "the TEST phenotype fixture now carries an ast_standard column, so "
            "NO_STANDARD_RECORDED is no longer guaranteed and the assertion "
            "below has become vacuous. Re-point it at a fixture that lacks one."
        )

    def test_the_provenance_wording_is_the_modules_own_constant(self):
        """The string is asserted against the module, not retyped.

        `NO_STANDARD_RECORDED` is documented as an assertion a reader may rely
        on, declared once for that reason. Hard-coding it here would make a
        deliberate rewording look like a regression and an accidental one look
        like a pass.
        """
        assert reporting.NO_STANDARD_RECORDED == (
            "testing standard not recorded in source data"
        )


class TestTheRunWrittenReportCarriesTheProvenance:
    """The assertion, on the artefact a run produced."""

    def test_the_provenance_section_is_present(self, report_text):
        assert "Stage 11 loaded" in report_text, (
            "the run wrote no phenotype provenance section. "
            "`reporting.phenotype_report` (reporting.py:172) returns None unless "
            "`ReportContext.phenotype_calls` is not None, so this means "
            "run.py's ReportContext construction is not passing it."
        )

    def test_the_section_reports_the_cohort_it_was_given(self, report_text):
        assert PROVENANCE_SENTENCE in report_text, (
            "the provenance section is present but states no R-vs-S contrast, so "
            "it is not the section `phenotype_provenance_section` builds."
        )

    def test_the_unrecorded_standard_is_stated_in_those_words(self, report_text):
        assert reporting.NO_STANDARD_RECORDED in report_text, (
            "the report does not state that the testing standard was not "
            "recorded. The fixture carries no ast_standard column, so every row "
            "lacks one and the report is required to say so rather than leave a "
            "reader to infer the gap."
        )

    def test_the_absence_wording_is_gone(self, report_text):
        assert "No stage 11 phenotype calls were handed to this report" not in (
            report_text
        ), (
            "the report still claims no phenotype calls were handed to it. That "
            "branch (reporting.py:368) fires only when phenotype_calls is None, "
            "so run.py is still not passing the calls it holds."
        )