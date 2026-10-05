"""REAL-mode reporting must refuse on INCOMPLETE input.

`convergence` and `cooccurrence` were given this guard by `46f5005` and its
follow-ups. Reporting had none, and the reason is structural rather than an
oversight: **reporting reads nothing from disk.** `write_report` takes a
`ReportContext` the orchestrator already assembled and renders it. So there was
no `load_*` call to guard, and nothing prompted the question.

The question still has an answer, and `_build_report_context` supplies it. Each
finding-bearing aggregate it renders is emitted **conditionally**
(`run.py:1615` gwas rows, `run.py:1653` convergence, `run.py:1670` significant
co-occurrence pairs). An empty one is dropped without a word, so the report that
reaches disk is titled `Pseudomonas aeruginosa Imipenem AMR Report`, banners
itself `REAL DATASET REPORT`, and is silently missing whole stages. That is the
artefact class `stages/real_inputs.py` exists to prevent.

What this file pins:

* it refuses **by input name**, with the path the aggregate is read from and
  what is holding it back - the shape `refuse_incomplete` already builds;
* an aggregate that is **present but degenerate** counts as absent. Two real
  cases, not invented ones: every stage-13 call classified `UNKNOWN`, and every
  stage-12 association with no adjusted p-value. Both render as a table that
  reads as "we looked and found nothing", which is indistinguishable from "we
  had no data" - the same failure `lineages={"S1": "unknown"}` already has a
  test for;
* it **stands down** when the declared inputs carry real content and the report
  is then written and computes. A guard never observed to stand down is not
  known to be a guard rather than a permanent refusal.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Sequence

import pytest

from papipeline.errors import StageError
from papipeline.models import (
    Cooccurrence,
    ConvergenceCall,
    ConvergenceCategory,
    RunMode,
)


# --------------------------------------------------------------------------
# Complete, real-shaped aggregates. These are what "everything present" looks
# like. Every value is real, not a sentinel - a fixture built from sentinels
# would be the degenerate case, not a passing one.
# --------------------------------------------------------------------------

SAMPLES = ("S1", "S2", "S3", "S4", "S5", "S6")


def _convergence_call(
    determinant: str = "blaOXA-1",
    category: ConvergenceCategory = ConvergenceCategory.RECURRENT_CONVERGENT,
) -> ConvergenceCall:
    return ConvergenceCall(
        determinant=determinant,
        independent_lineages=2,
        branch_count=2,
        distribution={"L1": 3, "L2": 3},
        convergence_category=category,
    )


def _cooccurrence() -> Cooccurrence:
    return Cooccurrence(
        feature_a="blaOXA-1",
        feature_b="fosA",
        feature_type="gene_gene",
        n_a=4,
        n_b=3,
        n_both=2,
        statistic="odds_ratio",
        statistic_value=4.0,
        adjusted_p_value=0.01,
    )


def _gwas_association() -> Dict[str, Any]:
    """A stage-12 association row, in the shape `figure_data` actually holds.

    Not a dataclass: `_prepare_figure_data` renders these as plain dicts
    (`run.py:1614` reads `figure_data["11_gwas_associations"]["rows"]`), and a
    fixture using a different shape would test the wrong object.
    """
    return {
        "feature": "gene__blaOXA-1",
        "feature_type": "acquired_gene",
        "adjusted_p_value": 0.001,
        "effect": 3.2,
        "frequency": 0.5,
        "lineage_linked": False,
        "passes_threshold": True,
        "dominant_lineage_share": 0.5,
    }


def complete_inputs() -> Dict[str, Any]:
    """Every declared REAL input, with real content."""
    return {
        "convergence_calls": [_convergence_call()],
        "cooccurrence": [_cooccurrence()],
        "gwas_associations": [_gwas_association()],
    }


def _context(config, mode: RunMode, declared: Dict[str, Any], **kw):
    from papipeline.stages.reporting import ReportContext

    return ReportContext(
        mode=mode,
        config=config,
        generated_at="2026-09-30T00:00:00Z",
        antibiotic="imipenem",
        n_samples=len(SAMPLES),
        declared_inputs=dict(declared),
        **kw,
    )


# --------------------------------------------------------------------------
# Refusal
# --------------------------------------------------------------------------


class TestReportingRefusesByInputName:
    def test_missing_convergence_calls_is_named_with_path_and_blocker(
        self, config, tmp_path: Path
    ):
        declared = complete_inputs()
        del declared["convergence_calls"]
        with pytest.raises(StageError) as excinfo:
            from papipeline.stages.reporting import write_report

            write_report(_context(config, RunMode.REAL, declared), tmp_path)
        message = str(excinfo.value)
        assert "convergence_calls" in message
        # The path it was expected at, and what is holding it back. A refusal
        # naming only the input is not actionable, and this one fires at the end
        # of a long DAG where the reader has no other signal.
        assert "13_convergence.tsv" in message
        assert "blocked by" in message

    def test_missing_cooccurrence_is_named_with_its_path(
        self, config, tmp_path: Path
    ):
        declared = complete_inputs()
        del declared["cooccurrence"]
        with pytest.raises(StageError) as excinfo:
            from papipeline.stages.reporting import write_report

            write_report(_context(config, RunMode.REAL, declared), tmp_path)
        message = str(excinfo.value)
        assert "cooccurrence" in message
        assert "14_cooccurrence.tsv" in message

    def test_missing_gwas_associations_is_named_with_its_path(
        self, config, tmp_path: Path
    ):
        declared = complete_inputs()
        del declared["gwas_associations"]
        with pytest.raises(StageError) as excinfo:
            from papipeline.stages.reporting import write_report

            write_report(_context(config, RunMode.REAL, declared), tmp_path)
        message = str(excinfo.value)
        assert "gwas_associations" in message
        assert "12_gwas.tsv" in message

    def test_it_names_every_absent_input_not_just_the_first(
        self, config, tmp_path: Path
    ):
        with pytest.raises(StageError) as excinfo:
            from papipeline.stages.reporting import write_report

            write_report(_context(config, RunMode.REAL, {}), tmp_path)
        assert excinfo.value.context["missing"].split(",") == [
            "convergence_calls",
            "cooccurrence",
            "gwas_associations",
        ]


# --------------------------------------------------------------------------
# The degenerate case: present, but carrying nothing usable
# --------------------------------------------------------------------------


class TestPresentButDegenerateCountsAsAbsent:
    def test_an_all_unknown_convergence_table_counts_as_absent(
        self, config, tmp_path: Path
    ):
        """The case a naive emptiness check would miss.

        Every call classified ``UNKNOWN`` is a stage *claiming* to have looked
        and reporting nothing classifiable. Rendered, that is a "Stage 13
        convergence classification" table full of ``unknown`` - which reads as
        a result. It is the same failure as ``lineages={"S1": "unknown"}``.
        """
        declared = complete_inputs()
        declared["convergence_calls"] = [
            _convergence_call(f"gene{i}", ConvergenceCategory.UNKNOWN)
            for i in range(3)
        ]
        with pytest.raises(StageError) as excinfo:
            from papipeline.stages.reporting import write_report

            write_report(_context(config, RunMode.REAL, declared), tmp_path)
        message = str(excinfo.value)
        assert "convergence_calls" in message
        assert "unknown" in message

    def test_a_single_real_classification_stands_the_sentinel_down(
        self, config, tmp_path: Path
    ):
        """One real row among UNKNOWNs is enough.

        ``is_effectively_empty`` is an ``all``, not an ``any``: a partially
        degenerate table still carries information, and refusing it would be a
        false alarm.
        """
        declared = complete_inputs()
        declared["convergence_calls"] = [
            _convergence_call("gene0", ConvergenceCategory.UNKNOWN),
            _convergence_call("blaOXA-1", ConvergenceCategory.RECURRENT_CONVERGENT),
        ]
        from papipeline.stages.reporting import write_report

        written = write_report(
            _context(config, RunMode.REAL, declared), tmp_path, write_html=False
        )
        assert written["markdown"].exists()

    def test_gwas_rows_with_no_adjusted_p_value_count_as_absent(
        self, config, tmp_path: Path
    ):
        """Associations that were never testable are not associations.

        Every row present but with ``adjusted_p_value=None`` means the stage
        computed and could test nothing - a table that renders as a result.
        """
        declared = complete_inputs()
        declared["gwas_associations"] = [dict(_gwas_association())]
        declared["gwas_associations"][0]["adjusted_p_value"] = None
        with pytest.raises(StageError) as excinfo:
            from papipeline.stages.reporting import write_report

            write_report(_context(config, RunMode.REAL, declared), tmp_path)
        assert "gwas_associations" in str(excinfo.value)

    def test_a_zero_row_aggregate_counts_as_absent(
        self, config, tmp_path: Path
    ):
        """The plain empty case: a table that exists and carries no rows."""
        declared = complete_inputs()
        declared["cooccurrence"] = []
        with pytest.raises(StageError) as excinfo:
            from papipeline.stages.reporting import write_report

            write_report(_context(config, RunMode.REAL, declared), tmp_path)
        assert "cooccurrence" in str(excinfo.value)


# --------------------------------------------------------------------------
# The guard must stand down
# --------------------------------------------------------------------------


class TestReportingRunsWhenInputsAreComplete:
    def test_it_writes_the_real_report_and_does_not_refuse(
        self, config, tmp_path: Path
    ):
        from papipeline.stages.reporting import write_report

        written = write_report(
            _context(config, RunMode.REAL, complete_inputs()),
            tmp_path,
            write_html=True,
        )
        assert written["markdown"].exists()
        assert written["html"].exists()
        text = written["markdown"].read_text(encoding="utf-8")
        # The real filename, i.e. the REAL branch was actually taken.
        assert written["markdown"].name == (
            "Pseudomonas_aeruginosa_Imipenem_AMR_Report.md"
        )
        assert "REAL DATASET REPORT" in text

    def test_the_aggregates_are_actually_rendered_when_present(
        self, config, tmp_path: Path
    ):
        """The guard standing down must not mean the report lost its content.

        Guards against a fix that satisfies the guard by dropping the tables.
        """
        from papipeline.stages.reporting import write_report

        context = _context(config, RunMode.REAL, complete_inputs())
        # As `_build_report_context` assembles them.
        context.tables = [
            (
                "Stage 13 convergence classification",
                ["determinant"],
                [c.to_row() for c in complete_inputs()["convergence_calls"]],
            )
        ]
        written = write_report(context, tmp_path, write_html=False)
        text = written["markdown"].read_text(encoding="utf-8")
        assert "Stage 13 convergence classification" in text
        assert "blaOXA-1" in text


# --------------------------------------------------------------------------
# The guard is REAL-only
# --------------------------------------------------------------------------


class TestOtherModesAreUnaffected:
    @pytest.mark.parametrize("mode", [RunMode.TEST, RunMode.STUB])
    def test_a_test_or_stub_report_is_never_refused(
        self, config, tmp_path: Path, mode: RunMode
    ):
        """TEST and STUB legitimately run on incomplete input.

        TEST is the 20 committed fixtures; STUB fabricates zero-row outputs by
        design. Refusing either would make the whole development and CI path
        unreachable, and would be exactly the "permanent refusal" this guard
        must not become.
        """
        from papipeline.stages.reporting import write_report

        written = write_report(
            _context(config, mode, {}), tmp_path, write_html=False
        )
        assert written["markdown"].exists()

    def test_a_test_report_keeps_its_synthetic_banner(self, config, tmp_path: Path):
        from papipeline.stages.reporting import write_report

        written = write_report(
            _context(config, RunMode.TEST, {}), tmp_path, write_html=False
        )
        assert "SYNTHETIC TEST DATA" in written["markdown"].read_text(
            encoding="utf-8"
        )


# --------------------------------------------------------------------------
# A refusal must not leave a partial artefact
# --------------------------------------------------------------------------


class TestRefusalLeavesNothingBehind:
    def test_no_markdown_and_no_html_reach_disk(
        self, config, tmp_path: Path
    ):
        """The whole point of refusing is that the artefact is absent.

        A half-written report stamped REAL is the dangerous artefact, not the
        exception - so the check is before the write, not after it.
        """
        from papipeline.stages.reporting import write_report

        out = tmp_path / "reports"
        with pytest.raises(StageError):
            write_report(_context(config, RunMode.REAL, {}), out)
        assert not list(out.glob("*.md")) if out.exists() else True
        assert not list(out.glob("*.html")) if out.exists() else True