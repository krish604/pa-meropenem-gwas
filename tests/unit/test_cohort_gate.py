"""The cohort join and size gate: stage 1 of the meropenem build.

The gate does one thing before anything else runs: it joins the assemblies
(manifest, what exists) to the meropenem S/I/R phenotype table (external
evidence) with the existing :func:`papipeline.join.join_manifest_to_phenotype`,
reports five counts - n tested, n with assembly, n R, n I, n S - and refuses a
cohort with too few resistant isolates to support an association scan.

Decisions pinned here:

* **Directional join, reused not reimplemented.** A manifest genome with no
  phenotype row (or two) is a hard failure; a phenotype row with no assembly
  is excluded, not failed - the phenotype table legitimately describes more
  isolates than were downloaded (spec.md D3).
* **Phenotype is binary S/I/R. There are no MICs.** Nothing here converts a
  category to a number, and ``require_trait=False`` is what the categorical
  join takes (spec.md D4).
* **The intermediate-isolate policy is configuration, not code**: the
  ``cohort_gate.intermediate_policy`` key, default ``exclude``. ``I`` isolates
  are always *counted*; the key decides whether they *enter* the analysis
  cohort. They are never folded into R or S either way.
* **The size threshold is configuration**: ``cohort_gate.min_resistant``,
  default 100. Falling below it is a refusal whose message names that key, so
  the fix (change the key) travels with the error.

Synthetic fixtures only: every manifest here is a hand-built tuple of
``TEST_PA_*`` identifiers and every phenotype table is written into
``tmp_path``. Nothing under ``data/``, ``db/`` or ``PDC_essential.tsv`` is
read, and no real-mode run happens.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Iterable, Optional, Sequence, Tuple

import pytest

from papipeline.cohort_gate import CohortGateError, evaluate_cohort_gate
from papipeline.config.loader import PipelineConfig, load_config
from papipeline.errors import (
    AntibioticError,
    ConfigError,
    DataContractError,
)

PIPELINE_ROOT = Path(__file__).resolve().parents[2]

#: The config keys the gate reads, as literals on purpose: the refusal message
#: must name a key a human can edit, so the test pins the name itself rather
#: than importing the constant the implementation also defines.
MIN_RESISTANT_KEY = "cohort_gate.min_resistant"
INTERMEDIATE_POLICY_KEY = "cohort_gate.intermediate_policy"

#: A four-isolate synthetic cohort: two resistant, one susceptible, one
#: intermediate - in manifest order.
MANIFEST_IDS: Tuple[str, ...] = (
    "TEST_PA_001",
    "TEST_PA_002",
    "TEST_PA_003",
    "TEST_PA_004",
)
MANIFEST_CALLS: Tuple[str, ...] = ("R", "R", "S", "I")

#: The phenotype rows the counts are asserted against: every manifest genome,
#: plus one row naming an isolate that was never assembled.
STANDARD_ROWS: Tuple[Tuple[str, str], ...] = tuple(
    zip(MANIFEST_IDS, MANIFEST_CALLS)
) + (("TEST_PA_099", "R"),)


def _gate_config(
    root: Path,
    *,
    min_resistant: Optional[int] = None,
    policy: Optional[str] = None,
    drop_policy: bool = False,
) -> PipelineConfig:
    """A self-contained copy of the repository configuration under ``root``.

    The cohort-gate keys are adjusted through the same YAML the run reads;
    nothing about the thresholds is hard-coded here beyond the text
    substitution, which fails loudly if the shipped defaults change.
    """
    root = Path(root)
    config_dir = root / "config"
    shutil.copytree(PIPELINE_ROOT / "config", config_dir)
    science = config_dir / "science.yaml"
    text = science.read_text(encoding="utf-8")
    if drop_policy:
        assert "  intermediate_policy: exclude\n" in text, (
            "the shipped science.yaml no longer carries the default policy line "
            "this test removes to prove the code default"
        )
        text = text.replace("  intermediate_policy: exclude\n", "")
    if policy is not None:
        assert "intermediate_policy: exclude" in text, (
            "the shipped science.yaml no longer carries `intermediate_policy: "
            "exclude`, which this test rewrites"
        )
        text = text.replace(
            "intermediate_policy: exclude", f"intermediate_policy: {policy}"
        )
    if min_resistant is not None:
        assert "min_resistant: 100" in text, (
            "the shipped science.yaml no longer carries `min_resistant: 100`, "
            "which this test rewrites"
        )
        text = text.replace("min_resistant: 100", f"min_resistant: {min_resistant}")
    science.write_text(text, encoding="utf-8")
    return load_config(science, machine=None)


def _write_meropenem_phenotype(
    phenotype_dir: Path, rows: Iterable[Tuple[str, str]]
) -> Path:
    """Write a synthetic meropenem phenotype TSV: sample_id, antibiotic, call."""
    phenotype_dir = Path(phenotype_dir)
    phenotype_dir.mkdir(parents=True, exist_ok=True)
    lines = ["sample_id\tantibiotic\tphenotype"]
    lines.extend(f"{sample_id}\tmeropenem\t{call}" for sample_id, call in rows)
    path = phenotype_dir / "meropenem_phenotype.tsv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _report(
    root: Path,
    rows: Sequence[Tuple[str, str]],
    *,
    run: str = "run",
    min_resistant: Optional[int] = None,
    policy: Optional[str] = None,
    drop_policy: bool = False,
    antibiotic: str = "meropenem",
):
    """Build a config copy and a phenotype table under ``root``, then gate."""
    _write_meropenem_phenotype(Path(root) / "phenotype", rows)
    config = _gate_config(
        Path(root) / run,
        min_resistant=min_resistant,
        policy=policy,
        drop_policy=drop_policy,
    )
    return evaluate_cohort_gate(
        config,
        manifest_ids=MANIFEST_IDS,
        phenotype_dir=Path(root) / "phenotype",
        antibiotic=antibiotic,
    )


class TestTheFiveCounts:
    """n tested, n with assembly, n R, n I, n S - reported, not guessed."""

    def test_it_reports_tested_assembly_and_category_counts(self, tmp_path):
        # Rows: 4 manifest genomes (R, R, S, I) + 1 row with no assembly (R).
        report = _report(tmp_path, STANDARD_ROWS, min_resistant=2)

        assert report.n_tested == 5
        assert report.n_with_assembly == 4
        assert report.n_no_assembly == 1
        assert report.n_resistant == 2
        assert report.n_intermediate == 1
        assert report.n_susceptible == 1
        assert report.n_in_cohort == 3
        # Default policy: the intermediate isolate is counted above and held
        # out of the cohort here.
        assert report.cohort_sample_ids == MANIFEST_IDS[:3]
        assert report.intermediate_policy == "exclude"
        assert report.min_resistant == 2

    def test_the_summary_line_carries_all_five_counts(self, tmp_path):
        report = _report(tmp_path, STANDARD_ROWS, min_resistant=2)

        summary = report.summary()
        for fragment in (
            "meropenem",
            "tested=5",
            "with_assembly=4",
            "R=2",
            "I=1",
            "S=1",
        ):
            assert fragment in summary, fragment


class TestTheIntermediatePolicy:
    """`cohort_gate.intermediate_policy` decides membership; never folding."""

    def test_intermediates_are_excluded_when_the_key_is_absent(self, tmp_path):
        # The key itself is removed from the config copy: what is asserted is
        # the CODE default, exclude, not the file's value.
        report = _report(tmp_path, STANDARD_ROWS, min_resistant=2, drop_policy=True)

        assert report.intermediate_policy == "exclude"
        assert MANIFEST_IDS[3] not in report.cohort_sample_ids
        assert report.n_intermediate == 1
        assert report.n_in_cohort == 3

    def test_intermediates_enter_the_cohort_when_the_key_says_include(self, tmp_path):
        report = _report(tmp_path, STANDARD_ROWS, min_resistant=2, policy="include")

        assert report.intermediate_policy == "include"
        assert MANIFEST_IDS[3] in report.cohort_sample_ids
        assert report.n_in_cohort == 4
        # Still counted as I, never re-labelled R or S.
        assert report.n_intermediate == 1
        assert report.n_resistant == 2
        assert report.n_susceptible == 1

    def test_an_unknown_policy_is_refused_naming_the_config_key(self, tmp_path):
        # `fold_into_r` is exactly the value that must never be accepted: it
        # would relabel a laboratory's I as a resistance call.
        with pytest.raises(ConfigError) as excinfo:
            _report(tmp_path, STANDARD_ROWS, min_resistant=2, policy="fold_into_r")
        assert INTERMEDIATE_POLICY_KEY in str(excinfo.value)


class TestTheSizeGate:
    """Below `cohort_gate.min_resistant`, the run stops and names the key."""

    def test_the_default_minimum_refuses_a_small_cohort_and_names_the_key(
        self, tmp_path
    ):
        # No overrides: the shipped configuration, minimum 100.
        with pytest.raises(CohortGateError) as excinfo:
            _report(tmp_path, STANDARD_ROWS)

        message = str(excinfo.value)
        assert MIN_RESISTANT_KEY in message
        assert "100" in message
        assert "meropenem" in message
        # The refusal carries the numbers, not just prose.
        assert excinfo.value.context["n_resistant"] == 2
        assert excinfo.value.context["min_resistant"] == 100
        assert excinfo.value.context["config_key"] == MIN_RESISTANT_KEY

    def test_the_minimum_is_configurable_through_the_key(self, tmp_path):
        with pytest.raises(CohortGateError) as excinfo:
            _report(tmp_path, STANDARD_ROWS, run="strict", min_resistant=3)
        assert MIN_RESISTANT_KEY in str(excinfo.value)
        assert excinfo.value.context["min_resistant"] == 3

        report = _report(tmp_path, STANDARD_ROWS, run="lenient", min_resistant=2)
        assert report.n_resistant == 2
        assert report.min_resistant == 2

    def test_the_gate_refuses_an_antibiotic_that_is_not_configured(self, tmp_path):
        with pytest.raises(AntibioticError):
            _report(
                tmp_path,
                STANDARD_ROWS,
                min_resistant=2,
                antibiotic="ceftazidime",
            )


class TestTheJoinIsTheExistingDirectionalOne:
    """The gate reuses papipeline.join; its hard failures stay hard."""

    def test_a_manifest_genome_with_no_phenotype_row_is_a_hard_failure(
        self, tmp_path
    ):
        # The fourth genome has no phenotype row: hard failure naming it.
        rows = tuple(zip(MANIFEST_IDS[:3], MANIFEST_CALLS[:3]))
        with pytest.raises(DataContractError) as excinfo:
            _report(tmp_path, rows, min_resistant=2)
        assert MANIFEST_IDS[3] in str(excinfo.value)

    def test_a_duplicate_phenotype_row_is_a_hard_failure(self, tmp_path):
        rows = STANDARD_ROWS + ((MANIFEST_IDS[0], "S"),)
        with pytest.raises(DataContractError):
            _report(tmp_path, rows, min_resistant=2)

    def test_a_row_with_no_assembly_is_excluded_not_failed(self, tmp_path):
        rows = STANDARD_ROWS + (("TEST_PA_098", "S"),)
        report = _report(tmp_path, rows, min_resistant=2)

        assert report.n_tested == 6
        assert report.n_with_assembly == 4
        assert report.n_no_assembly == 2
