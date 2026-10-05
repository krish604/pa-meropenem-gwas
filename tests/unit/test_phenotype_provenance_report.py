"""R7's accounting, and a production report that states it.

Two acceptance items live here because they are the same sentence in two
places.

**A1 - the primary phenotype is R vs S.** ``I`` and ``SDD`` are real, measured
categories. Ruling R7 excludes them from the contrast and **counts** them;
neither is ever merged into R or S. `binarise` already refused to guess, but a
refusal that leaves no trace is indistinguishable from a cohort with no
intermediates, so `provenance_report` now carries the counts and `load_phenotype`
logs them. A code that is not an AST category at all is a different mistake -
`RESISTANT` is not a spelling of `R`, and mapping it would be a guess about a
clinical record - so it gets its own refusal, and both name the code.

**A2 - the report says so.** `provenance_report` had no production caller: it
was reached from tests only, so the function that decides whether this cohort's
calls are comparable was never on the path that produces a report. It is now
called from `build_markdown`, and the REAL report states the AST provenance gap
in the words `NO_STANDARD_RECORDED`.

The real-data class reads the actual smoke AST table through the actual loader
and renders a REAL report from the result. It is not a reimplementation: the
counts come from `load_phenotype`, exactly as `run.py` calls it.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.config.loader import load_config
from papipeline.errors import PhenotypeError
from papipeline.models import (
    Cooccurrence,
    ConvergenceCall,
    ConvergenceCategory,
    RunMode,
)
from papipeline.stages import phenotype as stage
from papipeline.stages.reporting import NO_STANDARD_RECORDED, write_report

REPO = Path(__file__).resolve().parents[2]
SCIENCE = REPO / "config" / "science.yaml"
SMOKE_OVERLAY = REPO / "config" / "machines" / "smoke.yaml"
LAPTOP = REPO / "config" / "machines" / "laptop.yaml"

#: The measured composition of the smoke cohort. Fixed by the extractor's
#: refusal to fold an intermediate, so neither number may drift.
EXPECTED_R = 7
EXPECTED_S = 3
COHORT_N = EXPECTED_R + EXPECTED_S


def _call(sample_id, phenotype, **kw):
    return stage.PhenotypeCall(
        sample_id=sample_id,
        antibiotic="imipenem",
        phenotype=stage.Phenotype(phenotype),
        **kw,
    )


def _write_table(directory: Path, rows, header=None) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "imipenem_phenotype.tsv"
    lines = [header or "sample_id\tantibiotic\tphenotype\tMIC\tMIC_unit\tsource"]
    lines.extend(rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# --------------------------------------------------------------------------
# A1: the primary contrast, and what is held out of it
# --------------------------------------------------------------------------
class TestExcludedCategoriesAreCountedNotMerged:
    def test_r7_names_the_primary_contrast_and_the_exclusions(self):
        """The ruling, read from the place it is written down."""
        assert (stage.PRIMARY_POSITIVE, stage.PRIMARY_NEGATIVE) == ("R", "S")
        assert "I" in stage.EXCLUDED_CATEGORIES
        assert "SDD" in stage.EXCLUDED_CATEGORIES

    def test_the_module_constant_agrees_with_the_config_decision(self):
        """`EXCLUDED_CATEGORIES` restates `gwas.outcome.excluded`, so the two
        could drift - and a drift would mean the report counts a category the
        model actually folds. Asserted against the loaded config, so config
        stays the single decision and this module stays its restatement."""
        config = load_config(SCIENCE)
        assert set(stage.EXCLUDED_CATEGORIES) <= set(config.gwas.excluded)
        assert set(stage.EXCLUDED_CATEGORIES) <= set(
            p.value for p in config.allowed_phenotypes
        )
        # And nothing the model *does* compare may appear in the held-out set.
        assert not set(stage.EXCLUDED_CATEGORIES) & {
            stage.PRIMARY_POSITIVE,
            stage.PRIMARY_NEGATIVE,
        }

    def test_the_counts_separate_the_two_arms_from_the_held_out_rows(self):
        calls = (
            [_call(f"S{i}", "R") for i in range(7)]
            + [_call(f"T{i}", "S") for i in range(3)]
            + [_call("U1", "I"), _call("U2", "I"), _call("U3", "SDD")]
        )
        report = stage.provenance_report(calls)
        assert (report.positive, report.negative) == (7, 3)
        assert report.analysed == 10
        assert report.excluded == 3
        assert report.excluded_counts == {"I": 2, "SDD": 1}

    def test_the_excluded_rows_never_reach_either_arm(self):
        """Merging would show up as a larger cohort. Asserted on the count, not
        on the prose."""
        mixed = stage.provenance_report(
            [_call("A", "R"), _call("B", "S"), _call("C", "I")]
        )
        assert mixed.analysed == 2, "an intermediate must not enlarge the contrast"
        binarised = [stage.binarise(c, ["R"], ["S"]) for c in (
            _call("A", "R"), _call("B", "S"), _call("C", "I")
        )]
        assert binarised == [1, 0, None]

    def test_a_cohort_with_no_intermediates_says_none(self):
        report = stage.provenance_report([_call("A", "R"), _call("B", "S")])
        assert report.excluded == 0
        assert report.excluded_counts == {}
        assert "none" in report.summary()

    def test_the_summary_counts_both_arms_and_the_exclusions(self):
        report = stage.provenance_report(
            [_call("A", "R"), _call("B", "S"), _call("C", "SDD")]
        )
        text = report.summary()
        assert "R=1" in text and "S=1" in text
        assert "1 excluded and counted" in text
        assert "SDD=1" in text
        # The pre-existing provenance sentence must survive the addition.
        assert "may not be pooled confidently" in text

    def test_load_phenotype_logs_the_exclusions(self, caplog, tmp_path):
        """Counted in the log too: a run log listing I=6 beside R=7 and S=3,
        without saying so, reads as thirteen usable samples."""
        import logging

        config = load_config(SCIENCE)
        directory = tmp_path / "pheno"
        _write_table(directory, [
            "S1\timipenem\tR\t.\t.\treport",
            "S2\timipenem\tS\t.\t.\treport",
            "S3\timipenem\tI\t.\t.\treport",
        ])
        with caplog.at_level(logging.INFO, logger="papipeline.stages.phenotype"):
            stage.load_phenotype(config, directory, "imipenem")
        assert "excluded and counted, not merged into either arm" in caplog.text
        assert "I=1" in caplog.text

    def test_a_cohort_with_no_exclusions_logs_no_exclusion_line(self, caplog, tmp_path):
        """Absent, not zero. A line saying `excluded and counted: 0` on every
        ordinary run is noise that trains a reader to skip it."""
        import logging

        config = load_config(SCIENCE)
        directory = tmp_path / "pheno"
        _write_table(directory, [
            "S1\timipenem\tR\t.\t.\treport",
            "S2\timipenem\tS\t.\t.\treport",
        ])
        with caplog.at_level(logging.INFO, logger="papipeline.stages.phenotype"):
            stage.load_phenotype(config, directory, "imipenem")
        assert "excluded and counted" not in caplog.text


class TestUnknownCodesRefuseAndAreNamed:
    @pytest.mark.parametrize(
        "value", ["RESISTANT", "resistant", "RR", "1", "SENSITIVE?", "R/S"]
    )
    def test_a_code_that_is_not_an_ast_category_is_refused_and_named(self, value):
        config = load_config(SCIENCE)
        row = {"sample_id": "S1", "antibiotic": "imipenem", "phenotype": value}
        with pytest.raises(PhenotypeError) as excinfo:
            stage.parse_phenotype_row(
                row, config.allowed_phenotypes, config.antibiotics
            )
        message = str(excinfo.value)
        assert value.upper() in message, "the offending code must be named"
        assert "not a recognised AST category" in message
        # And the codes it would have been, so the reader can see the choice
        # being refused rather than infer it.
        for known in ("R", "S", "I", "SDD", "ND"):
            assert known in message

    def test_a_real_category_outside_the_allowed_set_is_refused_differently(self):
        """`I` is not a typo. It is a category the binary model refuses, so the
        message says that rather than calling it unrecognised - which would be
        a false statement about a perfectly valid AST code."""
        allowed = [stage.Phenotype.R, stage.Phenotype.S]
        row = {"sample_id": "S1", "antibiotic": "imipenem", "phenotype": "I"}
        with pytest.raises(PhenotypeError) as excinfo:
            stage.parse_phenotype_row(row, allowed, ["imipenem"])
        message = str(excinfo.value)
        assert "not in the allowed set" in message
        assert "real AST category" in message
        assert "not a recognised AST category" not in message
        assert "R7" in message, "the ruling that excludes it should be cited"

    @pytest.mark.parametrize("value", ["R", "I", "S", "SDD", "ND"])
    def test_every_ast_category_is_still_accepted_by_default(self, value):
        config = load_config(SCIENCE)
        row = {"sample_id": "S1", "antibiotic": "imipenem", "phenotype": value}
        call = stage.parse_phenotype_row(
            row, config.allowed_phenotypes, config.antibiotics
        )
        assert call.phenotype.value == value


# --------------------------------------------------------------------------
# A2: a PRODUCTION caller, and a REAL report that states the provenance
# --------------------------------------------------------------------------
def _complete_declared():
    """The three aggregates `refuse_if_incomplete` demands.

    Real-shaped and non-degenerate. They are declared here only to get past the
    guard so the provenance section can be reached: the assertions below are
    about what the report says about phenotype, not about convergence.
    """
    return {
        "convergence_calls": [
            ConvergenceCall(
                determinant="blaOXA-1",
                independent_lineages=2,
                branch_count=2,
                distribution={"L1": 3},
                convergence_category=ConvergenceCategory.RECURRENT_CONVERGENT,
            )
        ],
        "cooccurrence": [
            Cooccurrence(
                feature_a="blaOXA-1",
                feature_b="fosA",
                feature_type="gene_gene",
                n_a=4, n_b=3, n_both=2,
                statistic="odds_ratio",
                statistic_value=4.0,
                adjusted_p_value=0.01,
            )
        ],
        "gwas_associations": [
            {
                "feature": "gene__blaOXA-1",
                "feature_type": "acquired_gene",
                "adjusted_p_value": 0.001,
                "effect": 3.2,
                "frequency": 0.5,
                "lineage_linked": False,
                "passes_threshold": True,
                "dominant_lineage_share": 0.5,
            }
        ],
    }


def _real_context(config, calls, n_samples):
    from papipeline.stages.reporting import ReportContext

    return ReportContext(
        mode=RunMode.REAL,
        config=config,
        generated_at="2026-10-04T00:00:00Z",
        antibiotic="imipenem",
        n_samples=n_samples,
        declared_inputs=_complete_declared(),
        phenotype_calls=calls,
    )


class TestProvenanceReportHasAProductionCaller:
    def test_build_markdown_is_the_caller_not_a_test(self):
        """Named explicitly, because "a caller exists" was the finding.

        `provenance_report` was reachable only from tests. `build_markdown` is
        on the path of every report the pipeline writes, so this is the layer
        that makes the function production code.
        """
        import inspect

        from papipeline.stages import reporting

        source = inspect.getsource(reporting.build_markdown)
        assert "phenotype_report(context)" in source
        assert reporting.phenotype_report.__doc__ is not None
        # And it is the same function, not a reimplementation of the counting.
        assert (
            reporting.provenance_report is stage.provenance_report
        ), "reporting must call stage 11's function, not its own copy of it"

    def test_it_returns_none_when_no_phenotype_was_loaded(self):
        from papipeline.stages.reporting import phenotype_report

        config = load_config(SCIENCE)
        context = _real_context(config, [], COHORT_N)
        assert phenotype_report(context) is not None  # empty list, still a count
        context.phenotype_calls = None
        assert phenotype_report(context) is None

    def test_a_real_report_with_no_phenotype_states_the_absence(self, tmp_path):
        """An absence in REAL is a gap in the run. Silence would let a report
        assert an association without ever mentioning the phenotype it tested."""
        config = load_config(SCIENCE, machine=LAPTOP)
        context = _real_context(config, None, 100)
        written = write_report(context, tmp_path, write_html=False)
        text = written["markdown"].read_text(encoding="utf-8")
        assert "## Phenotype provenance" in text
        assert "No stage 11 phenotype calls were handed to this report" in text
        assert "This is an absence in the run, not a finding" in text

    def test_a_non_real_report_is_not_given_the_absence_wording(self):
        """The sentence is a REAL-mode claim about a real run. On a STUB report
        it would be noise, and it would be a false alarm."""
        from papipeline.stages.reporting import ReportContext, build_markdown

        config = load_config(SCIENCE, machine=LAPTOP)
        stub = ReportContext(
            mode=RunMode.STUB,
            config=config,
            generated_at="2026-10-04T00:00:00Z",
            antibiotic="imipenem",
            n_samples=0,
        )
        text = build_markdown(stub)
        assert "No stage 11 phenotype calls" not in text


# --------------------------------------------------------------------------
# The real cohort, end to end through the real loader
# --------------------------------------------------------------------------
needs_smoke_data = pytest.mark.skipif(
    not (REPO / "db" / "smoke_phenotypes" / "imipenem_phenotype.tsv").is_file(),
    reason="db/smoke_phenotypes/ is uncommitted by AGENTS.md rule 4 - absent "
           "in a fresh clone",
)


@needs_smoke_data
class TestTheRealCohortReport:
    """The acceptance criterion, as a reader would meet it."""

    @pytest.fixture(scope="class")
    def real_config(self):
        return load_config(SCIENCE, machine=SMOKE_OVERLAY)

    @pytest.fixture(scope="class")
    def calls(self, real_config):
        return stage.load_phenotype(
            real_config,
            real_config.phenotype_dir(RunMode.REAL),
            "imipenem",
        )

    def test_the_split_is_unchanged(self, calls):
        """R=7 / S=3, unchanged. The extractor's refusal to fold `I` is why."""
        report = stage.provenance_report(calls)
        assert (report.positive, report.negative) == (EXPECTED_R, EXPECTED_S)
        assert report.analysed == COHORT_N
        assert report.excluded == 0

    def test_every_row_lacks_an_ast_standard(self, calls):
        report = stage.provenance_report(calls)
        assert report.lacking_standard == report.total == COHORT_N
        assert not report.standard_recorded_everywhere

    def test_the_real_report_states_the_standard_was_not_recorded(self, tmp_path):
        """The literal words. Asserted as a phrase, not paraphrased, because the
        point is that the reader is told in English rather than left to infer a
        gap from an absent column."""
        config = load_config(SCIENCE, machine=LAPTOP)
        calls = stage.load_phenotype(
            load_config(SCIENCE, machine=SMOKE_OVERLAY),
            load_config(SCIENCE, machine=SMOKE_OVERLAY).phenotype_dir(RunMode.REAL),
            "imipenem",
        )
        context = _real_context(config, calls, COHORT_N)
        written = write_report(context, tmp_path, write_html=False)
        text = written["markdown"].read_text(encoding="utf-8")

        assert "## Phenotype provenance" in text
        assert NO_STANDARD_RECORDED in text
        assert "CLSI M100" in text, "the intended standard is named too"

    def test_the_report_labels_the_intended_standard_as_an_intention(
        self, tmp_path
    ):
        """`config/science.yaml` says `CLSI M100`. That is what this repository
        intended, not what the source laboratories did, and a reader who cannot
        tell those apart has been told a provenance fact that does not exist."""
        config = load_config(SCIENCE, machine=LAPTOP)
        calls = stage.load_phenotype(
            load_config(SCIENCE, machine=SMOKE_OVERLAY),
            load_config(SCIENCE, machine=SMOKE_OVERLAY).phenotype_dir(RunMode.REAL),
            "imipenem",
        )
        text = write_report(
            _real_context(config, calls, COHORT_N), tmp_path, write_html=False
        )["markdown"].read_text(encoding="utf-8")
        assert "intention held by this repository" in text
        assert "not a provenance fact about these isolates" in text

    def test_the_report_counts_the_arms_and_says_none_were_excluded(
        self, tmp_path
    ):
        config = load_config(SCIENCE, machine=LAPTOP)
        smoke = load_config(SCIENCE, machine=SMOKE_OVERLAY)
        calls = stage.load_phenotype(
            smoke, smoke.phenotype_dir(RunMode.REAL), "imipenem"
        )
        text = write_report(
            _real_context(config, calls, COHORT_N), tmp_path, write_html=False
        )["markdown"].read_text(encoding="utf-8")
        assert f"R (primary, positive arm) | {EXPECTED_R} |" in text
        assert f"S (primary, negative arm) | {EXPECTED_S} |" in text
        assert "excluded and counted, not merged" in text
        # And the absence is stated, so a reader knows this cohort has no
        # intermediates rather than wondering whether they went missing.
        assert "No record is excluded from that contrast" in text

    def test_a_smoke_report_carries_the_marker_and_the_provenance(
        self, tmp_path
    ):
        """Both together. The marker says how many isolates; the provenance
        says how comparable their calls are. A report with only the first is
        the failure the smoke overlay exists to prevent."""
        config = load_config(SCIENCE, machine=SMOKE_OVERLAY)
        calls = stage.load_phenotype(
            config, config.phenotype_dir(RunMode.REAL), "imipenem"
        )
        written = write_report(
            _real_context(config, calls, COHORT_N), tmp_path, write_html=False
        )
        text = written["markdown"].read_text(encoding="utf-8")
        assert "SMOKE TEST" in text
        assert NO_STANDARD_RECORDED in text

    def test_excluded_rows_reach_the_report_when_present(self, tmp_path):
        """Rendered, not just counted. Built by hand because the smoke cohort
        has no intermediates - which is exactly why this needs its own case."""
        config = load_config(SCIENCE, machine=LAPTOP)
        calls = (
            [_call(f"R{i}", "R") for i in range(EXPECTED_R)]
            + [_call(f"S{i}", "S") for i in range(EXPECTED_S)]
            + [_call("I1", "I"), _call("I2", "I"), _call("D1", "SDD")]
        )
        text = write_report(
            _real_context(config, calls, COHORT_N), tmp_path, write_html=False
        )["markdown"].read_text(encoding="utf-8")
        assert "I=2" in text
        assert "SDD=1" in text
        assert "3 record(s) are excluded from that contrast" in text
        assert "ruling R7 forbids merging" in text
        # And the arms are unaffected by the extras.
        assert f"R (primary, positive arm) | {EXPECTED_R} |" in text
        assert f"S (primary, negative arm) | {EXPECTED_S} |" in text

    def test_a_partially_sourced_cohort_is_not_called_wholly_unsourced(
        self, tmp_path
    ):
        """Two different sentences for two different facts: no row has a
        standard, or some rows do not. One sentence for both would over-claim
        on a partially sourced cohort."""
        config = load_config(SCIENCE, machine=LAPTOP)
        calls = [
            _call("A", "R", ast_standard="CLSI M100"),
            _call("B", "S"),
        ]
        report = stage.provenance_report(calls)
        assert report.lacking_standard == 1
        text = write_report(
            _real_context(config, calls, 2), tmp_path, write_html=False
        )["markdown"].read_text(encoding="utf-8")
        assert "1 of 2 records carry no AST standard" in text
        assert "**The testing standard not recorded in source data.**" not in text