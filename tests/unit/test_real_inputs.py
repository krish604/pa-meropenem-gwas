"""The shared REAL-input refusal predicate.

Extracted from `test_convergence_real_inputs.py` so the predicate can be
committed - and reviewed - on its own, before either stage depends on it.

What is being pinned, and why it is a module rather than a copy in each stage:

* An **empty** mapping and a **degenerate** one are different failures with the
  same symptom, and only the second needs a sentinel. `{"S1": []}` is data
  saying "these samples carry nothing"; `{"S1": "unknown"}` is data *claiming*
  to carry lineage information and having none. Checking the first as empty
  would refuse a merely uninformative run; checking the second as present
  would classify every determinant `UNKNOWN` and read as "no convergence found".
* A refusal that does not **name** the missing input is not actionable, and
  these fire deep in a DAG where the reader has no other signal. So the message
  carries the input name, the path it was expected at, and what holds it back -
  the same shape `adapters.external.require_tool` uses.
"""

from __future__ import annotations

import pytest

from papipeline.errors import StageError
from papipeline.stages import real_inputs


class TestIsEffectivelyEmpty:
    def test_absent_mapping_is_empty(self):
        assert real_inputs.is_effectively_empty(None)
        assert real_inputs.is_effectively_empty({})

    def test_a_mapping_of_empty_lists_is_not_empty_without_a_sentinel(self):
        """The asymmetry is deliberate.

        ``{"S1": [], "S2": []}`` is *present* data saying these samples carry
        nothing. Checking it for emptiness would refuse a run that is merely
        uninformative, which is a different thing from an absent input - and it
        is why cooccurrence checks the built namespaces rather than the raw
        mappings it was handed.
        """
        assert not real_inputs.is_effectively_empty({"S1": [], "S2": []})

    def test_all_sentinel_is_empty(self):
        assert real_inputs.is_effectively_empty(
            {"S1": "unknown", "S2": "unknown"}, sentinel="unknown"
        )

    def test_one_real_value_defeats_the_sentinel(self):
        assert not real_inputs.is_effectively_empty(
            {"S1": "unknown", "S2": "L1"}, sentinel="unknown"
        )

    def test_a_mixed_sentinel_is_not_empty(self):
        assert not real_inputs.is_effectively_empty(
            {"S1": "unknown", "S2": "unknown", "S3": "L1"}, sentinel="unknown"
        )


class TestMissingRequiredInputs:
    def test_nothing_missing_is_empty(self):
        assert real_inputs.missing_required_inputs(
            {"amr_calls": {"S1": ["a"]}, "lineages": {"S1": "L1"}}
        ) == {}

    def test_it_reports_the_empty_one_by_name(self):
        missing = real_inputs.missing_required_inputs(
            {"amr_calls": {"S1": ["a"]}, "regulator_variants": {}}
        )
        assert list(missing) == ["regulator_variants"]
        assert missing["regulator_variants"]

    def test_an_all_sentinel_input_reports_the_sentinel_reason(self):
        missing = real_inputs.missing_required_inputs(
            {"lineages": {"S1": "unknown"}}, sentinels={"lineages": "unknown"}
        )
        assert "unknown" in missing["lineages"]
        assert "sentinel" in missing["lineages"]

    def test_a_configured_sentinel_is_not_applied_to_other_inputs(self):
        """`lineages: unknown` must not make an absent `amr_calls` pass."""
        missing = real_inputs.missing_required_inputs(
            {"amr_calls": {}, "lineages": {"S1": "L1"}},
            sentinels={"lineages": "unknown"},
        )
        assert list(missing) == ["amr_calls"]


class TestRefuseIncompleteNamesEverything:
    def test_it_names_each_missing_input(self):
        err = real_inputs.refuse_incomplete(
            "convergence",
            {"lineages": "nothing carries a lineage"},
            paths={"lineages": "/tmp/tree_metadata.tsv"},
            blocked_by={"lineages": "stage 9"},
        )
        message = str(err)
        assert "convergence" in message
        assert "lineages" in message
        assert "/tmp/tree_metadata.tsv" in message
        assert "stage 9" in message

    def test_it_reports_the_count_of_missing_inputs(self):
        err = real_inputs.refuse_incomplete(
            "cooccurrence", {"gene": "empty", "mutation": "empty"}
        )
        assert "2 of its declared inputs" in str(err)

    def test_context_is_machine_readable(self):
        err = real_inputs.refuse_incomplete("cooccurrence", {"gene": "empty"})
        assert isinstance(err, StageError)
        assert err.context["stage"] == "cooccurrence"
        assert err.context["missing"] == "gene"

    def test_it_raises_nothing_itself(self):
        """It *builds* the error so a caller can attach context and re-raise.

        Building rather than raising is what lets the two stages name their own
        paths and blockers in one place; the test is that a bare call does not
        raise on its own.
        """
        err = real_inputs.refuse_incomplete("convergence", {"lineages": "x"})
        assert isinstance(err, Exception)

    @pytest.mark.parametrize("name", ["lineages", "gene", "mutation"])
    def test_every_input_name_appears_verbatim(self, name):
        """A refusal that paraphrases the input name is not actionable.

        The name is what a reader greps for, so it has to survive verbatim into
        both the message and the context.
        """
        err = real_inputs.refuse_incomplete("stage", {name: "absent"})
        assert name in str(err)
        assert name in err.context["missing"]