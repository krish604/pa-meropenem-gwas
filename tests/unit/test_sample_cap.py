"""Tests for the laptop sample cap.

Ticket 04. A cap is only useful if it stops the run, and only tolerable if the
refusal tells the user what to do. Both halves are asserted here, because a cap
that fires with an unactionable message is an obstacle rather than a safeguard.

Written against the ticket's acceptance criteria.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.config.loader import load_config, load_machine_config
from papipeline.errors import SampleCapExceeded

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
SCIENCE = PIPELINE_ROOT / "config" / "science.yaml"
LAPTOP = PIPELINE_ROOT / "config" / "machines" / "laptop.yaml"
BIGMACHINE = PIPELINE_ROOT / "config" / "machines" / "bigmachine.yaml"


# ---------------------------------------------------------------------------
# the boundary
# ---------------------------------------------------------------------------


def test_exactly_twenty_samples_is_allowed():
    """The cap is inclusive; 20 is the documented limit, not one below it."""
    load_machine_config(LAPTOP).enforce_sample_cap(20)


def test_twenty_one_samples_is_refused():
    with pytest.raises(SampleCapExceeded):
        load_machine_config(LAPTOP).enforce_sample_cap(21)


def test_zero_and_one_are_allowed():
    laptop = load_machine_config(LAPTOP)
    laptop.enforce_sample_cap(0)
    laptop.enforce_sample_cap(1)


def test_the_full_analysis_cohort_is_refused_on_the_laptop():
    with pytest.raises(SampleCapExceeded):
        load_machine_config(LAPTOP).enforce_sample_cap(900)


# ---------------------------------------------------------------------------
# the refusal must be actionable
# ---------------------------------------------------------------------------


def test_the_refusal_names_the_limit():
    with pytest.raises(SampleCapExceeded) as excinfo:
        load_machine_config(LAPTOP).enforce_sample_cap(21)
    assert "20" in str(excinfo.value)


def test_the_refusal_names_the_observed_count():
    with pytest.raises(SampleCapExceeded) as excinfo:
        load_machine_config(LAPTOP).enforce_sample_cap(835)
    assert "835" in str(excinfo.value)


def test_the_refusal_names_the_key_and_the_file_to_change():
    with pytest.raises(SampleCapExceeded) as excinfo:
        load_machine_config(LAPTOP).enforce_sample_cap(21)
    message = str(excinfo.value)
    assert "max_samples" in message, "must name the config key"
    assert "laptop.yaml" in message, "must name the file to change"


def test_the_refusal_points_at_the_machine_that_can_do_it():
    """A refusal that does not say what to do instead is half a refusal."""
    with pytest.raises(SampleCapExceeded) as excinfo:
        load_machine_config(LAPTOP).enforce_sample_cap(900)
    assert "bigmachine" in str(excinfo.value)


def test_the_refusal_is_a_pipeline_error_so_callers_catch_one_type():
    from papipeline.errors import PipelineError

    with pytest.raises(PipelineError):
        load_machine_config(LAPTOP).enforce_sample_cap(21)


# ---------------------------------------------------------------------------
# the analysis machine is deliberately uncapped
# ---------------------------------------------------------------------------


def test_the_analysis_machine_accepts_the_full_cohort():
    load_machine_config(BIGMACHINE).enforce_sample_cap(900)


def test_uncapped_is_distinguishable_from_unset():
    """`None` must mean "checked, and there is no cap", not "not looked at"."""
    assert load_machine_config(BIGMACHINE).max_samples is None
    assert load_machine_config(LAPTOP).max_samples == 20


# ---------------------------------------------------------------------------
# reachable through the composed config
# ---------------------------------------------------------------------------


def test_the_cap_is_reachable_from_a_loaded_config():
    config = load_config(SCIENCE, machine="laptop")
    with pytest.raises(SampleCapExceeded):
        config.enforce_sample_cap(21)


def test_a_self_contained_config_has_no_cap_to_enforce():
    """machine=None means the file is everything; there is no overlay to cap it."""
    config = load_config(SCIENCE, machine=None)
    assert config.max_samples is None
    config.enforce_sample_cap(900)  # must not raise
