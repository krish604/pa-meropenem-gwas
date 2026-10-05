"""The real ``mlst`` invocation, and the v4-style output it produces.

`science.yaml` already settles the science: ``mlst.tool: mlst``,
``mlst.scheme: "paeruginosa"``, ``min_typing_locus_match: 5``. What was missing
was the code that runs the tool. So what is tested here is that the code
*reads* that configuration rather than restating it, and that it reports what
the tool actually said.

Three properties are load-bearing.

**Every flag is one this ``mlst`` actually has.** Verified against the installed
binary's ``--help`` rather than assumed. In particular ``--legacy`` is *not* a
cosmetic choice: the help records that it "requires --scheme", and without it
the allele header row is absent, so the loci cannot be named. An ST is
therefore useless without it - the pipeline's contract keeps the per-locus
profile precisely so an ST can be traced back to the alleles behind it
(`MlstCall.alleles`), and the header is the only place the locus names come
from. A wrong flag here is worse than no code.

**A missing locus is not a zero allele.** ``mlst`` reports an unmatched locus
as ``-`` and a partial one as ``?``. Neither is an allele, and both must not
reach ``MlstCall.alleles``: ``parse_alleles`` treats a null sentinel as an
unmatched locus, so passing them through would be caught, but counting them
would inflate the matched-loci total and could promote a 2-of-7 profile to
``typed``.

**Fewer than `min_typing_locus_match` loci is ``partial``, with no ST.** Not
``typed``, and not a guessed ST. That rule already lives in
`stages.mlst.parse_row`, so the caller hands back raw tool output and lets the
stage apply it - one implementation of the status contract, not two.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.adapters import mlst as adapter
from papipeline.errors import PipelineError
from papipeline.config.loader import PipelineConfig

# Verbatim `mlst --legacy` output, captured from the installed 2.23 tool against
# db/smoke_genomes/PDT000167136.1.fna. Real, not invented: the ST and the seven
# allele calls are that isolate's.
LEGACY_TYPED = (
    "FILE\tSCHEME\tST\tacsA\taroE\tguaA\tmutL\tnuoD\tppsA\ttrpE\n"
    "PDT000167136.1\tpaeruginosa\t235\t38\t11\t3\t13\t1\t2\t4\n"
)


class TestTheCommandIsBuiltFromConfig:
    def test_the_scheme_comes_from_config(self, config: PipelineConfig):
        command = adapter.query_command(
            "mlst",
            scheme=config.raw["mlst"]["scheme"],
            assembly="a.fna",
            threads=4,
        )
        assert "--scheme" in command
        assert "paeruginosa" in command

    def test_legacy_is_requested_because_the_loci_come_from_the_header(
        self, config: PipelineConfig
    ):
        """Without the header row the loci are unnameable, so the profile is."""
        command = adapter.query_command(
            "mlst",
            scheme=config.raw["mlst"]["scheme"],
            assembly="a.fna",
            threads=1,
        )
        assert "--legacy" in command

    def test_threads_come_from_the_overlay(self, config: PipelineConfig):
        command = adapter.query_command(
            "mlst",
            scheme=config.raw["mlst"]["scheme"],
            assembly="a.fna",
            threads=int(config.runtime["threads"]),
        )
        assert "--threads" in command
        assert str(int(config.runtime["threads"])) in command

    def test_quiet_is_requested_so_stderr_carries_no_prose(
        self, config: PipelineConfig
    ):
        """mlst writes progress to stderr. Parsing stdout is the contract, and
        progress text there would be read as data."""
        command = adapter.query_command(
            "mlst",
            scheme=config.raw["mlst"]["scheme"],
            assembly="a.fna",
            threads=1,
        )
        assert "--quiet" in command

    def test_nothing_is_hard_coded(self, config: PipelineConfig):
        """No scheme literal in the adapter: it must not restate the config.

        A hard-coded scheme would pass every functional test while silently
        disagreeing with the file that records the science.
        """
        command = adapter.query_command(
            "mlst",
            scheme=config.raw["mlst"]["scheme"],
            assembly="a.fna",
            threads=1,
        )
        assert "paeruginosa" not in command[: command.index("--scheme")]

    def test_the_datadir_is_passed_when_given(self, tmp_path):
        """The scheme data must be the one the config provisioned, not whatever
        the tool happens to default to on this machine."""
        command = adapter.query_command(
            "mlst", scheme="paeruginosa", assembly="a.fna", threads=1,
            datadir=str(tmp_path / "pubmlst"),
        )
        assert "--datadir" in command
        assert str(tmp_path / "pubmlst") in command


class TestParsingRealOutput:
    def test_a_typed_profile_yields_its_alleles(self):
        parsed = adapter.parse_legacy_output(LEGACY_TYPED, "PDT000167136.1")
        assert parsed.st == "235"
        assert parsed.alleles == {
            "acsA": "38", "aroE": "11", "guaA": "3", "mutL": "13",
            "nuoD": "1", "ppsA": "2", "trpE": "4",
        }
        assert len(parsed.alleles) == 7

    def test_the_loci_come_from_the_header_not_from_a_hard_coded_list(self):
        """The scheme defines the loci. Hard-coding seven names would silently
        mis-parse any other scheme the config ever names."""
        header = "FILE\tSCHEME\tST\tygiZ\nS1\tscheme\t9\t7\n"
        parsed = adapter.parse_legacy_output(header, "S1")
        assert parsed.alleles == {"ygiZ": "7"}

    def test_an_unmatched_locus_is_not_an_allele(self):
        """mlst writes `-` for a locus it did not match. Counting it would
        inflate the matched-loci total and could promote a partial profile to
        typed."""
        row = "FILE\tSCHEME\tST\tacsA\taroE\tguaA\nS1\ts\t1\t-\t11\t-\n"
        parsed = adapter.parse_legacy_output(row, "S1")
        assert parsed.alleles == {"aroE": "11"}

    def test_a_partial_allele_is_kept_and_flagged(self):
        """`?` means a partial match, not a clean call.

        Dropping it would understate how much of the profile was actually
        determined, so it is reported as partial rather than discarded.
        """
        row = "FILE\tSCHEME\tST\tacsA\taroE\tguaA\nS1\ts\t-\t38\t?\t3\n"
        parsed = adapter.parse_legacy_output(row, "S1")
        assert parsed.partial_loci == ["aroE"]
        assert parsed.alleles.get("aroE") == "?"

    def test_no_st_when_the_profile_is_incomplete(self):
        """mlst writes `-` for an unresolved ST. It must not become a number."""
        row = "FILE\tSCHEME\tST\tacsA\taroE\nS1\ts\t-\t1\t2\n"
        assert adapter.parse_legacy_output(row, "S1").st is None

    def test_an_empty_result_is_not_an_error(self):
        """mlst exits 0 and prints a header with no rows for an assembly it
        cannot type. That is `no_call`, not a tool failure - conflating them
        would fail a whole cohort over one untypable genome."""
        parsed = adapter.parse_legacy_output(
            "FILE\tSCHEME\tST\tacsA\n", "S1"
        )
        assert parsed.st is None
        assert parsed.alleles == {}
        assert not parsed.found

    def test_a_row_is_matched_to_its_sample_not_to_line_order(self):
        """The FILE column is authoritative. A single-file call is the common
        case, but matching on position would silently mis-assign a result."""
        out = (
            "FILE\tSCHEME\tST\tacsA\n"
            "other.fna\ts\t9\t1\n"
            "wanted.fna\ts\t8\t2\n"
        )
        parsed = adapter.parse_legacy_output(out, "other.fna")
        assert parsed.st == "9"


class TestTheResultBecomesAStageRow:
    """The caller hands raw tool output to the stage, which owns the schema.

    One implementation of the status contract: a second copy of the
    `min_typing_locus_match` rule here would be free to drift from
    `stages.mlst.parse_row`.
    """

    def test_it_produces_a_row_parse_row_understands(self, config: PipelineConfig):
        from papipeline.stages.mlst import parse_row

        parsed = adapter.parse_legacy_output(LEGACY_TYPED, "PDT000167136.1")
        row = parsed.as_stage_row("pubmlst:paeruginosa:2026-03-11")
        call = parse_row(row, config)
        assert call.sample_id == "PDT000167136.1"
        assert call.sequence_type == "235"
        assert call.mlst_status == "typed"
        assert call.mlst_scheme == "pseudomonas_aeruginosa"

    def test_a_three_locus_profile_is_demoted_by_the_stage(
        self, config: PipelineConfig
    ):
        """Below min_typing_locus_match=5, so `partial` with no ST. Asserted
        through the real `parse_row`, so the rule is proven to be shared rather
        than duplicated."""
        from papipeline.stages.mlst import parse_row

        row = "FILE\tSCHEME\tST\tacsA\taroE\tguaA\tmutL\nS1\ts\t99\t1\t2\t3\t4\n"
        parsed = adapter.parse_legacy_output(row, "S1")
        call = parse_row(parsed.as_stage_row("db"), config)
        assert call.mlst_status == "partial"
        assert call.sequence_type is None, (
            "a 4-of-7 profile must not carry the ST mlst printed; the ST is "
            "only defined for a full profile and keeping it would be a guess"
        )


class TestTheStageRunsTheToolOnARealRun:
    """`stages.mlst.run` reads a table; on REAL there is no table to read.

    So the stage gains a REAL branch that calls the tool and applies the same
    `parse_row` contract the TEST branch uses. Reusing `parse_row` is the point:
    a REAL run and a TEST run with identical profiles must produce identical
    `MlstCall`s, or the fixture stops describing the behaviour.
    """

    def _manifest(self, *sample_ids):
        from papipeline.manifest import SampleManifest
        from papipeline.models import Sample

        return SampleManifest(samples=[Sample(sample_id=s) for s in sample_ids])

    def _assemblies(self, root, *sample_ids):
        """Real files under `<data_root>/<sample_id>/`, so locate_assembly works.

        The stage refuses to substitute one isolate's sequence for another's,
        which is right; a test that skipped this would be testing a refusal.
        """
        for sample_id in sample_ids:
            directory = root / sample_id
            directory.mkdir(parents=True, exist_ok=True)
            (directory / f"{sample_id}.fna").write_text(
                ">x\nACGT\n", encoding="utf-8"
            )
        return root

    def test_a_real_run_calls_the_tool_per_isolate(
        self, config: PipelineConfig, tmp_path, monkeypatch
    ):
        from papipeline.models import RunMode
        from papipeline.stages import mlst as stage

        seen = []

        def fake_call(sample_id, **kwargs):
            seen.append(sample_id)
            return adapter.parse_legacy_output(
                LEGACY_TYPED.replace("PDT000167136.1", sample_id), sample_id
            )

        monkeypatch.setattr(adapter, "call_isolate", fake_call)
        data_root = self._assemblies(tmp_path / "data", "PDT_A", "PDT_B")
        calls = stage.run(
            config, self._manifest("PDT_A", "PDT_B"), RunMode.REAL,
            intermediate_root=tmp_path, data_root=data_root,
        )
        assert sorted(seen) == ["PDT_A", "PDT_B"]
        assert calls["PDT_A"].sequence_type == "235"
        assert calls["PDT_B"].sequence_type == "235"

    def test_a_test_run_still_reads_the_fixture(
        self, config: PipelineConfig, tmp_path
    ):
        """The TEST branch is untouched. If REAL had replaced it, TEST would
        need a tool installed - which is exactly what it exists to avoid."""
        from papipeline.models import RunMode
        from papipeline.stages import mlst as stage

        target = tmp_path / "mlst"
        target.mkdir()
        (target / "mlst_results.tsv").write_text(
            "sample_id\tST\talleles\tMLST_status\tmlst_scheme\tallele_database\n"
            "PDT_A\t7\tacsA:1;aroE:2;guaA:3;mutL:4;nuoD:5\ttyped\t"
            "pseudomonas_aeruginosa\tsynthetic\n",
            encoding="utf-8",
        )
        calls = stage.run(
            config, self._manifest("PDT_A"), RunMode.TEST,
            intermediate_root=tmp_path, data_root=tmp_path,
        )
        assert calls["PDT_A"].sequence_type == "7"

    def test_an_isolate_with_no_assembly_is_recorded_not_dropped(
        self, config: PipelineConfig, tmp_path, monkeypatch
    ):
        """A key for every manifest sample, always.

        `align_to_manifest` already guarantees the key exists, so a failed call
        must not raise past it: one unanalysable genome should cost that genome
        its ST, not the cohort's report. The same reasoning the aligner
        documents for 82% of assemblies being corrupt.
        """
        from papipeline.models import RunMode
        from papipeline.stages import mlst as stage

        def fake_call(sample_id, **kwargs):
            if sample_id == "PDT_MISSING":
                raise PipelineError("no assembly for PDT_MISSING")
            return adapter.parse_legacy_output(
                LEGACY_TYPED.replace("PDT000167136.1", sample_id), sample_id
            )

        monkeypatch.setattr(adapter, "call_isolate", fake_call)
        data_root = self._assemblies(tmp_path / "data", "PDT_A")
        calls = stage.run(
            config, self._manifest("PDT_A", "PDT_MISSING"), RunMode.REAL,
            intermediate_root=tmp_path, data_root=data_root,
        )
        assert set(calls) == {"PDT_A", "PDT_MISSING"}
        assert calls["PDT_MISSING"] is None
        assert calls["PDT_A"] is not None

    def test_the_provenance_names_the_scheme_actually_used(
        self, config: PipelineConfig, tmp_path, monkeypatch
    ):
        """An ST is only reproducible if the scheme that produced it is named.

        The config's `paeruginosa` is a legacy spelling that normalises to
        `pseudomonas_aeruginosa`; the recorded value has to say which was used,
        or re-running against a different scheme is undetectable.
        """
        from papipeline.models import RunMode
        from papipeline.stages import mlst as stage

        monkeypatch.setattr(
            adapter, "call_isolate",
            lambda sample_id, **kw: adapter.parse_legacy_output(
                LEGACY_TYPED.replace("PDT000167136.1", sample_id), sample_id
            ),
        )
        data_root = self._assemblies(tmp_path / "data", "PDT_A")
        calls = stage.run(
            config, self._manifest("PDT_A"), RunMode.REAL,
            intermediate_root=tmp_path, data_root=data_root,
        )
        assert calls["PDT_A"].allele_database, (
            "the allele profile was recorded with no statement of which scheme "
            "produced it, so the ST cannot be reproduced or audited"
        )
