"""Phenotype: source calls, measured MICs, and the continuous trait.

Ticket 07. Three things are being asserted here, and one of them is a change of
scientific approach the user asked for:

* The source laboratory's R/S/I call is authoritative and is never overwritten.
* MIC-to-R/S re-interpretation runs **only** on rows carrying a genuinely
  measured MIC, and only when the configuration supplies breakpoints. The
  thresholds are CLSI M100, which is paywalled, so no numbers are present in
  this repository. Breakpoints are therefore *optional*: with none configured
  the pipeline says so and proceeds, rather than inventing a threshold.
* The association model consumes a **continuous** trait, log2(MIC), not the
  binary R/S label. More power from the same data, and it is what pyseer's
  continuous mode expects. This is recorded in the spec as a change from the
  earlier binary-outcome decision.

AST provenance travels with each row, and the run reports how many rows lack it -
an MIC from one laboratory's method is not automatically comparable with
another's.

Written before the implementation.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from papipeline.config.loader import load_config
from papipeline.errors import ConfigError, PhenotypeError
from papipeline.stages.phenotype import (
    interpret_mic,
    load_phenotype,
    phenotype_matrix,
)

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
SCIENCE = PIPELINE_ROOT / "config" / "science.yaml"

HEADER = "sample_id\tantibiotic\tphenotype\tMIC\tMIC_unit\tsource\n"


def _write(tmp_path: Path, rows, header: str = HEADER, name="imipenem_phenotype.tsv"):
    directory = tmp_path / "phenotype"
    directory.mkdir(exist_ok=True)
    path = directory / name
    path.write_text(header + "".join(rows), encoding="utf-8")
    return directory


# ---------------------------------------------------------------------------
# log2 is the continuous trait
# ---------------------------------------------------------------------------


def test_log2_of_a_measured_mic_is_the_continuous_trait():
    """The model's phenotype is log2(MIC), not the R/S label."""
    assert interpret_mic(2.0, base=2) == pytest.approx(1.0)
    assert interpret_mic(8.0, base=2) == pytest.approx(3.0)
    assert interpret_mic(0.5, base=2) == pytest.approx(-1.0)


def test_log2_is_computed_numerically_not_by_a_lookup_table():
    # 0.0625 -> -4 exactly; 0.2 -> log2(0.2), not rounded to a breakpoint
    assert interpret_mic(0.0625, base=2) == pytest.approx(-4.0)
    assert interpret_mic(0.2, base=2) == pytest.approx(math.log2(0.2))


def test_an_unmeasurable_mic_has_no_continuous_trait():
    assert interpret_mic(0, base=2) is None, "log2(0) is undefined"
    assert interpret_mic(-1, base=2) is None, "a negative MIC is not a measurement"
    assert interpret_mic(None, base=2) is None


def test_the_natural_log_base_is_available_for_comparison():
    assert interpret_mic(8.0, base=10) == pytest.approx(math.log10(8.0))


# ---------------------------------------------------------------------------
# breakpoints are optional and config-driven
# ---------------------------------------------------------------------------


def test_the_science_config_declares_the_breakpoint_section():
    section = load_config(SCIENCE, machine="laptop").raw.get("breakpoints") or {}
    assert section, "science.yaml must declare a breakpoints section"
    assert "standard" in section
    assert "log_base" in section


def test_no_threshold_values_are_shipped_in_this_repository():
    """CLSI M100 is paywalled. Inventing a number here would be the worst
    possible failure: an unsourced threshold that reaches a result looks
    exactly like a sourced one."""
    import yaml

    section = yaml.safe_load(SCIENCE.read_text(encoding="utf-8"))["breakpoints"]
    thresholds = section.get("thresholds")
    assert not thresholds, (
        f"breakpoint thresholds must stay empty until supplied from a licensed "
        f"copy; found {thresholds!r}"
    )


def test_a_row_with_no_thresholds_configured_is_not_reinterpreted(tmp_path: Path):
    """With no breakpoints, the source call stands and the pipeline says so."""
    directory = _write(tmp_path, ["GCA_000000001.1\timipenem\tS\t8\tmg/L\treport\n"])
    config = load_config(SCIENCE, machine="laptop")
    calls = load_phenotype(config, directory, "imipenem")
    call = calls[0]
    assert call.phenotype.value == "S", "the source call must survive untouched"
    assert call.mic == 8.0
    assert call.reinterpreted_from_mic is False
    assert call.log2_mic == pytest.approx(3.0)


def test_the_run_reports_that_no_standard_was_applied(tmp_path: Path):
    from papipeline.stages.phenotype import BreakpointStatus

    directory = _write(tmp_path, ["GCA_000000001.1\timipenem\tS\t8\tmg/L\treport\n"])
    config = load_config(SCIENCE, machine="laptop")
    status = BreakpointStatus.from_config(config)
    assert status.is_configured is False
    assert status.standard == "CLSI M100"
    assert status.edition, "the intended edition is recorded even when unapplied"
    reason = status.reason.lower()
    assert "thresholds" in reason and "configured" in reason, (
        "the reason must say that no thresholds are configured"
    )
    assert "clsi m100" in reason, "and name the standard it intends to use"
    assert "source calls stand" in reason, "and say the source calls are untouched"


# ---------------------------------------------------------------------------
# the source call is authoritative
# ---------------------------------------------------------------------------


def test_a_categorical_call_never_gains_a_mic(tmp_path: Path):
    directory = _write(tmp_path, ["GCA_000000001.1\timipenem\tR\t.\t.\treport\n"])
    config = load_config(SCIENCE, machine="laptop")
    call = load_phenotype(config, directory, "imipenem")[0]
    assert call.mic is None
    assert call.log2_mic is None
    assert call.phenotype.value == "R"


def test_a_measured_mic_is_preserved_exactly(tmp_path: Path):
    directory = _write(tmp_path, ["GCA_000000001.1\timipenem\tR\t11.1\tmg/L\treport\n"])
    config = load_config(SCIENCE, machine="laptop")
    call = load_phenotype(config, directory, "imipenem")[0]
    assert call.mic == 11.1
    assert call.mic_unit == "mg/L"


# ---------------------------------------------------------------------------
# AST provenance
# ---------------------------------------------------------------------------


PROVENANCE_HEADER = (
    "sample_id\tantibiotic\tphenotype\tMIC\tMIC_unit\tast_method\tast_standard\tast_edition\n"
)


def test_ast_provenance_is_carried_per_row(tmp_path: Path):
    directory = _write(
        tmp_path,
        ["GCA_000000001.1\timipenem\tR\t8\tmg/L\tbroth microdilution\tCLSI M100\t35th ed.\n"],
        header=PROVENANCE_HEADER,
    )
    config = load_config(SCIENCE, machine="laptop")
    call = load_phenotype(config, directory, "imipenem")[0]
    assert call.ast_method == "broth microdilution"
    assert call.ast_standard == "CLSI M100"
    assert call.ast_edition == "35th ed."


def test_the_run_reports_how_many_rows_lack_provenance(tmp_path: Path):
    directory = _write(
        tmp_path,
        [
            "GCA_000000001.1\timipenem\tR\t8\tmg/L\tbroth microdilution\tCLSI M100\t35th ed.\n",
            "GCA_000000002.1\timipenem\tS\t.\t.\t.\t.\t.\n",
            "GCA_000000003.1\timipenem\tR\t.\t.\tdisk diffusion\tEUCAST\t12.0\n",
        ],
        header=PROVENANCE_HEADER,
    )
    config = load_config(SCIENCE, machine="laptop")
    calls = load_phenotype(config, directory, "imipenem")
    from papipeline.stages.phenotype import provenance_report

    report = provenance_report(calls)
    assert report.total == 3
    assert report.with_mic == 1
    assert report.lacking_method == 1
    assert report.lacking_standard == 1
    assert report.lacking_edition == 1


def test_absent_provenance_columns_are_reported_not_fabricated(tmp_path: Path):
    directory = _write(tmp_path, ["GCA_000000001.1\timipenem\tR\t8\tmg/L\treport\n"])
    config = load_config(SCIENCE, machine="laptop")
    calls = load_phenotype(config, directory, "imipenem")
    from papipeline.stages.phenotype import provenance_report

    report = provenance_report(calls)
    assert report.lacking_method == 1, "an absent column means every row lacks it"
    assert report.lacking_standard == 1
    assert report.lacking_edition == 1


# ---------------------------------------------------------------------------
# the continuous matrix pyseer will consume
# ---------------------------------------------------------------------------


def test_the_matrix_carries_log2_mic_where_measured(tmp_path: Path):
    directory = _write(
        tmp_path,
        [
            "GCA_000000001.1\timipenem\tR\t8\tmg/L\treport\n",
            "GCA_000000002.1\timipenem\tS\t0.5\tmg/L\treport\n",
        ],
    )
    config = load_config(SCIENCE, machine="laptop")
    calls = load_phenotype(config, directory, "imipenem")
    traits = phenotype_matrix(calls)
    assert traits["GCA_000000001.1"] == pytest.approx(3.0)
    assert traits["GCA_000000002.1"] == pytest.approx(-1.0)


def test_a_row_without_a_measured_mic_is_absent_from_the_matrix_not_zero(tmp_path: Path):
    """A categorical-only row has no continuous trait. Zero would be a
    measurement, and it would silently pull the model toward it."""
    directory = _write(
        tmp_path,
        [
            "GCA_000000001.1\timipenem\tR\t8\tmg/L\treport\n",
            "GCA_000000002.1\timipenem\tS\t.\t.\treport\n",
        ],
    )
    config = load_config(SCIENCE, machine="laptop")
    calls = load_phenotype(config, directory, "imipenem")
    traits = phenotype_matrix(calls)
    assert "GCA_000000002.1" not in traits
    assert 0.0 not in traits.values()


# ---------------------------------------------------------------------------
# unchanged strictness
# ---------------------------------------------------------------------------


def test_a_non_numeric_mic_is_still_rejected(tmp_path: Path):
    directory = _write(tmp_path, ["GCA_000000001.1\timipenem\tR\thigh\tmg/L\treport\n"])
    config = load_config(SCIENCE, machine="laptop")
    with pytest.raises(PhenotypeError):
        load_phenotype(config, directory, "imipenem")


def test_an_unknown_category_is_still_rejected(tmp_path: Path):
    directory = _write(tmp_path, ["GCA_000000001.1\timipenem\tVERY-R\t.\t.\treport\n"])
    config = load_config(SCIENCE, machine="laptop")
    with pytest.raises(PhenotypeError):
        load_phenotype(config, directory, "imipenem")
