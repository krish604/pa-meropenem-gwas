"""Stage 6 must say which isolates it could not read.

The four fixes specified in `pa-artifacts/round12/FIX-MPILEUP.md`. The trigger
is recorded here because the tests below only make sense against it: the
round-12 run shipped `variants.tsv` covering 9 of 10 genomes, the provenance
sidecar carried `no_assembly` for the tenth in a JSON file that exists only on
disk, the cohort table's `an=10` reported itself as complete, and a reader of
the output had no way to learn that one genome had never been read.

Everything asserted here is about **stating** a shortfall, never about changing
one. In particular `af` keeps its denominator: an isolate that could not be
read stays in `an`, because counting an unread isolate as a non-carrier is what
`docs/scientific_rules.md` rule 9 already requires of missing data.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict

import pytest

from papipeline.errors import PipelineError, ToolExecutionError
from papipeline.stages import cohort_variants, variants
from papipeline.stages.cohort_variants import merge_calls
from papipeline.stages.variants import (
    NOT_CALLED_STATUSES,
    STATUS_CALLED,
    STATUS_CALL_FAILED,
    STATUS_CALL_FAILED_SIGNAL,
    STATUS_NO_ASSEMBLY,
    call_provenance,
)


def _row(ref: str, alt: str, sample_id: str = "S1", pos: str = "1") -> Dict[str, str]:
    return {
        "sample_id": sample_id, "chrom": "NC_002516.2", "pos": pos,
        "ref": ref, "alt": alt, "qual": "30.4", "filter": "PASS", "GT": "1",
        "AC": "1", "AN": "1", "DP4": "0,0,1,0", "MQ": "60", "MQ0F": "0",
    }


# ---------------------------------------------------------------------------
# Fix 1: `call_status` in the per-isolate provenance record.
# ---------------------------------------------------------------------------


class TestCallStatusIsRecordedPerIsolate:
    def test_the_verdict_is_carried_by_the_isolate_that_had_it(self):
        calls = {"S1": [_row("A", "C")], "S2": []}
        statuses = {"S1": STATUS_CALLED, "S2": STATUS_NO_ASSEMBLY}
        payload = call_provenance(calls, statuses=statuses)
        by_id = {e["sample_id"]: e for e in payload["per_isolate"]}
        assert by_id["S1"]["call_status"] == "called"
        assert by_id["S2"]["call_status"] == "no_assembly"

    def test_an_empty_call_list_no_longer_proves_the_isolate_was_read(self):
        """The ambiguity this exists to break.

        Both an isolate that read cleanly and carried nothing, and one that was
        never read, have an empty row list. Before `call_status` the record was
        identical for both.
        """
        calls = {"S1": [], "S2": []}
        statuses = {"S1": STATUS_CALLED, "S2": STATUS_NO_ASSEMBLY}
        by_id = {
            e["sample_id"]: e
            for e in call_provenance(calls, statuses=statuses)["per_isolate"]
        }
        assert by_id["S1"]["n_alleles"] == by_id["S2"]["n_alleles"] == 0
        assert by_id["S1"]["call_status"] != by_id["S2"]["call_status"]

    def test_the_status_is_absent_not_invented_when_nobody_tracked_it(self):
        """Rule 9: unknown must not be dressed up as a result.

        `call_status` is omitted entirely rather than emitted as `"unknown"`,
        so a reader cannot mistake a placeholder for an observation.
        """
        payload = call_provenance({"S1": []})
        assert "call_status" not in payload["per_isolate"][0]

    def test_the_counts_of_called_and_not_called_appear_with_the_status(self):
        statuses = {"S1": STATUS_CALLED, "S2": STATUS_CALL_FAILED}
        payload = call_provenance({"S1": [_row("A", "C")], "S2": []}, statuses=statuses)
        assert payload["totals"]["n_called"] == 1
        assert payload["totals"]["n_not_called"] == 1

    def test_those_counts_are_absent_without_statuses_to_derive_them_from(self):
        assert "n_called" not in call_provenance({"S1": []})["totals"]


# ---------------------------------------------------------------------------
# Fix 1 (cont.): a tool killed by the OS is not a tool that exited.
# ---------------------------------------------------------------------------


class TestAKilledToolIsClassifiedDifferently:
    """A negative return code is a signal, and signals are not exit codes."""

    def test_a_signal_death_is_named_as_such(self):
        killed = ToolExecutionError("mpileup failed", returncode=-9)
        assert variants._failure_status(killed) == STATUS_CALL_FAILED_SIGNAL

    def test_an_ordinary_failure_is_not_reported_as_a_signal(self):
        exited = ToolExecutionError("mpileup failed", returncode=1)
        assert variants._failure_status(exited) == STATUS_CALL_FAILED

    def test_no_return_code_means_plain_failure(self):
        assert (
            variants._failure_status(ToolExecutionError("mpileup failed"))
            == STATUS_CALL_FAILED
        )

    def test_a_pre_tool_failure_is_not_attributed_to_a_tool(self):
        """A missing assembly never launched mpileup.

        Reporting it as `call_failed` would send a reader to a log line that
        cannot exist, so anything without a return code is a plain failure and
        the `no_assembly` branch upstream keeps handling its own case.
        """
        assert (
            variants._failure_status(PipelineError("no assembly for S1"))
            == STATUS_CALL_FAILED
        )

    def test_only_the_four_documented_statuses_exist(self):
        """A fifth would need a row in the contract, not a literal."""
        assert {STATUS_CALLED, STATUS_NO_ASSEMBLY, STATUS_CALL_FAILED,
                STATUS_CALL_FAILED_SIGNAL} == {
            "called", "no_assembly", "call_failed", "call_failed_signal",
        }
        # Every failure mode minus `called` is "not called", and that set is
        # what `an_calls` counts against.
        assert NOT_CALLED_STATUSES == {
            STATUS_NO_ASSEMBLY, STATUS_CALL_FAILED, STATUS_CALL_FAILED_SIGNAL,
        }


# ---------------------------------------------------------------------------
# Fix 3: `an_calls` beside `an`, with `af` left alone.
# ---------------------------------------------------------------------------

#: Round 12 as a fixture, because that is the run these fixes describe: ten
#: isolates, three carrying a site, `iso10` never read at all. The cohort is
#: sized so the site still clears `min_dissenters=2` - seven isolates do not
#: carry it - which is what makes it worth asserting anything about.
SITE = ("NC_002516.2", "100", "C", "T")

CARRIERS = ("iso1", "iso2", "iso3")
UNREAD = "iso10"


def _round12(all_called: bool = False):
    """The cohort, and stage 6's verdict on each isolate."""
    calls = {
        f"iso{i}": [SITE] if f"iso{i}" in CARRIERS else []
        for i in range(1, 11)
    }
    statuses = {f"iso{i}": STATUS_CALLED for i in range(1, 11)}
    if not all_called:
        statuses[UNREAD] = STATUS_NO_ASSEMBLY
    return calls, statuses


class TestAnCallsIsReportedNextToAn:
    def test_an_stays_the_cohort_size_and_an_calls_is_the_readable_part(self):
        calls, statuses = _round12()
        rows = merge_calls(calls, statuses=statuses)
        assert rows, "the site must survive the filter for this to test anything"
        for row in rows:
            assert row["an"] == "10"
            assert row["an_calls"] == "9"

    def test_the_unread_isolate_stays_in_the_af_denominator(self):
        """Fix 3's whole point: say it without changing it.

        `af` is still ``ac / an``. The unread isolate counts as a non-carrier
        because it was never observed to carry the site, and dropping it would
        renumber every frequency the pipeline has already reported.
        """
        calls, statuses = _round12()
        for row in merge_calls(calls, statuses=statuses):
            assert float(row["af"]) == pytest.approx(int(row["ac"]) / int(row["an"]))
            assert float(row["af"]) == pytest.approx(3 / 10)

    def test_nothing_is_reported_when_status_was_not_tracked(self):
        """Omitting `statuses` must not fabricate a shortfall."""
        calls, _ = _round12()
        for row in merge_calls(calls):
            assert row["an_calls"] == row["an"] == "10"

    def test_a_cohort_with_no_failures_reports_an_equal_to_an(self):
        calls, statuses = _round12(all_called=True)
        for row in merge_calls(calls, statuses=statuses):
            assert row["an"] == row["an_calls"] == "10"

    def test_the_column_is_in_the_documented_order(self):
        assert cohort_variants.MERGE_COLUMNS == (
            "chrom", "pos", "ref", "alt", "ac", "an", "an_calls", "af",
        )


# ---------------------------------------------------------------------------
# Fix 2: the denominator is stated in the log, and failures are NAMED.
# ---------------------------------------------------------------------------


class TestTheDenominatorIsStatedAloud:
    def test_the_warning_names_the_isolates_that_produced_none(self, caplog):
        """A count cannot be acted on; an isolate ID can.

        This is the test that would have caught round 12: `an=10` was printed,
        the fact that `iso10` had never been read was not.
        """
        calls, statuses = _round12()
        with caplog.at_level("WARNING", logger="stages.cohort_variants"):
            cohort_variants.run(None, _manifest(list(calls)), _calls(calls), statuses=statuses)
        warning = [r for r in caplog.records if r.levelname == "WARNING"]
        assert warning, "a shortfall must be warned about, not only counted"
        message = warning[0].getMessage()
        assert UNREAD in message
        assert "an=10" in message
        assert "9 with calls" in message

    def test_it_is_quiet_when_every_isolate_was_read(self, caplog):
        calls, statuses = _round12(all_called=True)
        with caplog.at_level("WARNING", logger="stages.cohort_variants"):
            cohort_variants.run(None, _manifest(list(calls)), _calls(calls), statuses=statuses)
        assert not [r for r in caplog.records if r.levelname == "WARNING"]


# ---------------------------------------------------------------------------
# Fix 4 / plumbing: `run` carries the statuses from stage 6 into stage 6a.
# ---------------------------------------------------------------------------


class _Manifest:
    def __init__(self, sample_ids):
        self.sample_ids = sample_ids


def _manifest(sample_ids):
    return _Manifest(sample_ids)


def _calls(calls):
    """Stage 6's shape: per-isolate rows as dicts, not tuples."""
    return {
        sample_id: [
            {
                "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                "qual": "30.4", "filter": "PASS", "GT": "1", "AC": "1",
                "AN": "1", "DP4": "0,0,1,0", "MQ": "60", "MQ0F": "0",
            }
            for (chrom, pos, ref, alt) in rows
        ]
        for sample_id, rows in calls.items()
    }


class TestRunHandsStatusesThrough:
    def test_run_reaches_merge_calls_with_the_statuses(self, caplog):
        """The in-memory path, not a re-read of the sidecar."""
        calls, statuses = _round12()
        with caplog.at_level("WARNING", logger="stages.cohort_variants"):
            rows = cohort_variants.run(
                None, _manifest(list(calls)), _calls(calls), statuses=statuses
            )
        assert rows
        for row in rows:
            assert row["an"] == "10"
            assert row["an_calls"] == "9"

    def test_run_still_works_for_callers_that_pass_no_statuses(self):
        """The signature keeps its old call shape working."""
        calls, _ = _round12()
        rows = cohort_variants.run(None, _manifest(list(calls)), _calls(calls))
        for row in rows:
            assert row["an"] == row["an_calls"]


# ---------------------------------------------------------------------------
# The REAL-mode loop itself, where the verdicts are actually assigned.
#
# The doubles in `tests/integration/` replace `stage_variants.run` wholesale, so
# none of them reach these branches. Everything above tests the helpers in
# isolation; this drives `_run_by_aligning` end to end with the aligner and the
# assembler removed, which is the only way to assert that a killed tool and an
# absent assembly are recorded where they claim to be.
# ---------------------------------------------------------------------------


class _LoopConfig:
    """Just enough config for `_run_by_aligning` to reach the loop body."""

    def __init__(self):
        self.runtime = {"threads": 1}

    def reference_fasta(self):
        return Path("pinned_reference.fa")


class _LoopManifest:
    sample_ids = ("ok", "noasm", "killed")

    def require(self, sample_id):
        return sample_id


# Verdicts the loop is expected to produce, one isolate per branch.
EXPECTED = {
    "ok": STATUS_CALLED,
    "noasm": STATUS_NO_ASSEMBLY,
    "killed": STATUS_CALL_FAILED_SIGNAL,
}


@pytest.fixture()
def real_loop(monkeypatch, tmp_path):
    """`_run_by_aligning` with only the two external tools removed."""
    import papipeline.adapters.minimap2 as aligner

    def locate(sample, data_root):
        if sample == "noasm":
            raise PipelineError(f"no assembly for {sample}")
        return tmp_path / f"{sample}.fna"

    def call(sample_id, **kwargs):
        if sample_id == "killed":
            # Exactly what the round-12 mpileup produced.
            raise ToolExecutionError("mpileup was killed", returncode=-9)
        # `_run_by_aligning` reads the file BEFORE handing its text to the
        # parser, so a path alone is not enough - it must exist.
        vcf = tmp_path / f"{sample_id}.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
        return type("Result", (), {"vcf": str(vcf)})()

    parsed = []

    def parse(text, sample_id=None):
        parsed.append(sample_id)
        return [{"sample_id": sample_id, "chrom": "NC_002516.2", "pos": "1",
                 "ref": "A", "alt": "G"}]

    monkeypatch.setattr(variants, "locate_assembly", locate)
    monkeypatch.setattr(variants, "call_isolate", call)
    monkeypatch.setattr(variants, "parse_vcf", parse)
    # `call_settings` has its own tests; this one is about verdicts, so it
    # takes a settings dict rather than a whole config's `raw` block.
    monkeypatch.setattr(variants, "call_settings", lambda config: {})
    monkeypatch.setattr(aligner, "require_callers", lambda tools: {})

    return tmp_path, parsed


class TestTheRealLoopRecordsEachVerdict:
    def test_each_isolate_gets_the_verdict_its_branch_reached(self, real_loop):
        workdir, _ = real_loop
        out: Dict[str, str] = {}
        variants._run_by_aligning(
            _LoopConfig(), _LoopManifest(), data_root=workdir,
            workdir=workdir / "work", statuses_out=out,
        )
        assert out == EXPECTED

    def test_an_isolate_that_died_still_becomes_a_cohort_member_with_no_calls(
        self, real_loop
    ):
        """82% of the round-12 assemblies were corrupt.

        The stage must complete anyway, and record the isolate as carrying
        nothing rather than dropping it - dropping it would inflate every
        `af` in the merge.
        """
        workdir, _ = real_loop
        calls = variants._run_by_aligning(
            _LoopConfig(), _LoopManifest(), data_root=workdir,
            workdir=workdir / "work",
        )
        assert set(calls) == {"ok", "noasm", "killed"}
        assert calls["killed"] == [] and calls["noasm"] == []
        assert calls["ok"], "the isolate that read cleanly still produces rows"

    def test_only_the_readable_isolate_reaches_the_vcf_parser(self, real_loop):
        workdir, parsed = real_loop
        variants._run_by_aligning(
            _LoopConfig(), _LoopManifest(), data_root=workdir,
            workdir=workdir / "work",
        )
        assert parsed == ["ok"]

    def test_the_sidecar_records_the_same_verdicts_in_memory(self, real_loop):
        """Log, sidecar and out-param must not be three accounts of one run."""
        workdir, _ = real_loop
        work = workdir / "work"
        variants._run_by_aligning(
            _LoopConfig(), _LoopManifest(), data_root=workdir,
            workdir=work, statuses_out={},
        )
        payload = json.loads(
            (work / variants.PROVENANCE_NAME).read_text(encoding="utf-8")
        )
        by_id = {e["sample_id"]: e for e in payload["per_isolate"]}
        assert {sid: by_id[sid]["call_status"] for sid in EXPECTED} == EXPECTED
        assert payload["totals"]["isolates"] == 3
        assert payload["totals"]["n_called"] == 1
        assert payload["totals"]["n_not_called"] == 2

    def test_the_out_param_is_replaced_not_appended_to(self, real_loop):
        """A rerun must not inherit the previous run's verdicts."""
        workdir, _ = real_loop
        stale = {"gone": STATUS_CALLED}
        variants._run_by_aligning(
            _LoopConfig(), _LoopManifest(), data_root=workdir,
            workdir=workdir / "work", statuses_out=stale,
        )
        assert set(stale) == set(EXPECTED)

    def test_no_out_param_is_fine(self, real_loop):
        workdir, _ = real_loop
        assert variants._run_by_aligning(
            _LoopConfig(), _LoopManifest(), data_root=workdir,
            workdir=workdir / "work",
        )
