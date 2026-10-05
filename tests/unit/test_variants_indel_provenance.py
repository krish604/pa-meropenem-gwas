"""Ruling R6 end to end: `min_ireads` reaches the bcftools argv, and a REAL run
that calls no indels says so.

Two findings are being pinned here, and they are the same finding seen from
opposite ends.

**The pinned bcftools default for `--min-ireads` is 2, and that is what
suppressed every assembly indel.** Verified against the installed 1.23.1, whose
own usage line reads::

    -m, --min-ireads INT    Minimum number gapped reads for indel candidates [2]

(`bcftools mpileup --help` does not exist on 1.23.1 - it exits 1 with
``unrecognized option `--help'`` - so the usage block printed by a bare
`bcftools mpileup` is the authoritative text. That is the same observation
`tests/unit/test_minimap2_adapter.py::TestTheFlagsAreReal` records.)

The decisive measurement is a paired control on **identical alignments**: the
same BAM, the same `-q/-Q/-a/-d`, once with `-m 1` and once without. The SNV sets
come out byte-identical, so base quality, BAQ and `-Q 20` cannot be what
suppressed the indels - the only thing that changed is `-m`. Without it: 4
indels per isolate and **zero** inside oprD. With it: 284 and 277, including the
2 nt deletion and the 5 nt insertion at the two oprD positions the SAM CIGAR
walk found. So a test asserting the number appears in a command vector is
worth having - but only because the number is load-bearing, and this file also
asserts the whole chain that carries it.

**The stage cannot tell from its own output that it called nothing.** A REAL run
over the 10-isolate smoke subset produced 539,527 calls of which 28 were
non-SNV, none inside oprD - and the table was well formed, because the contract
carries no column saying what kind of allele a row is. That is the failure
`summarise_call_provenance` exists to make impossible to repeat quietly: the
per-isolate SNV/indel counts are recorded, and a REAL run with zero indels
cohort-wide logs a WARNING saying why that is not a result.

The warning is asserted by its substring, `NO_INDEL_WARNING`, and the code emits
that same constant - so the test cannot pass against a message that stopped
mentioning the reason.
"""

from __future__ import annotations

import json
import logging
import shutil
from pathlib import Path
from typing import Any, Dict, List, Mapping

import pytest

from papipeline.adapters import minimap2 as adapter
from papipeline.config.loader import DEFAULT_MPILEUP_MIN_IREADS, PipelineConfig
from papipeline.models import RunMode
from papipeline.stages import variants as stage

#: The installed binary's own usage text, verbatim, for the line R6 rests on.
#: Kept as a literal so the assertion below is against the tool's words and not
#: against this file's memory of them.
MPILEUP_MIN_IREADS_USAGE = (
    "  -m, --min-ireads INT    Minimum number gapped reads for indel "
    "candidates [2]"
)


def _pileup_argv(config: PipelineConfig) -> List[str]:
    """The mpileup argv the STAGE would build, end to end.

    Deliberately not a hand-assembled call. The chain under test is
    ``science.yaml -> PipelineConfig -> call_settings -> call_isolate ->
    mpileup_command``, and a test that builds the argv itself proves only that
    the last link works. This walks the first three and asserts on the vector,
    which is the thing the tool is actually handed.
    """
    settings = stage.call_settings(config)
    pileup_cfg = settings["bcftools"]["mpileup"]
    return adapter.mpileup_command(
        "bcftools",
        bam=Path("in.bam"),
        reference=Path("ref.fa"),
        min_mapq=int(pileup_cfg["min_mapq"]),
        min_bq=int(pileup_cfg["min_bq"]),
        annotate=str(pileup_cfg["annotate"]),
        max_depth=int(pileup_cfg["max_depth"]),
        min_ireads=int(pileup_cfg["min_ireads"]),
    )


def _row(ref: str, alt: str, sample_id: str = "S1", pos: str = "1") -> Dict[str, str]:
    return {
        "sample_id": sample_id, "chrom": "NC_002516.2", "pos": pos,
        "ref": ref, "alt": alt, "qual": "30.4", "filter": "PASS", "GT": "1",
        "AC": "1", "AN": "1", "DP4": "0,0,1,0", "MQ": "60", "MQ0F": "0",
    }


# ---------------------------------------------------------------------------
# A1: the tool's own default, and the configured value in the argv.
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not shutil.which("bcftools"),
    reason="bcftools is not installed on this machine, so its own default "
           "cannot be read; the pinned value is asserted from science.yaml "
           "by tests/unit/test_round11_config_keys.py instead",
)
class TestTheToolSaysItsDefaultIsTwo:
    """The installed binary is the authority on its own default."""

    def test_the_usage_line_says_two(self):
        assert MPILEUP_MIN_IREADS_USAGE in _installed_mpileup_usage()

    def test_help_is_not_an_option_so_the_usage_block_is_the_text(self):
        """`bcftools mpileup --help` exits 1 on 1.23.1; record that.

        Otherwise the obvious way to check a tool's default silently produces
        an error message and no usage at all.
        """
        import subprocess

        out = subprocess.run(
            ["bcftools", "mpileup", "--help"], capture_output=True, text=True,
            check=False,
        )
        assert out.returncode != 0
        assert "unrecognized option" in (out.stderr or "").lower()

    def test_the_default_bcftools_supplies_is_not_the_one_we_configure(self):
        """If these ever agreed, R6 would have nothing to enforce."""
        assert "[2]" in MPILEUP_MIN_IREADS_USAGE
        assert DEFAULT_MPILEUP_MIN_IREADS == 1


class TestTheConfiguredValueReachesTheArgv:
    """The whole chain, not just the last link."""

    def test_the_stage_builds_an_argv_carrying_the_configured_value(
        self, config: PipelineConfig
    ):
        command = _pileup_argv(config)
        assert "-m" in command, (
            "the stage's mpileup argv has no --min-ireads at all, so the "
            "pinned default of 2 applies whatever science.yaml says"
        )
        assert command[command.index("-m") + 1] == "1", command

    def test_it_is_not_the_tools_own_default(self, config: PipelineConfig):
        """Guards against the flag being present and set to the wrong number."""
        command = _pileup_argv(config)
        assert command[command.index("-m") + 1] != "2", command

    def test_the_value_survives_a_config_that_omits_the_key(
        self, config: PipelineConfig
    ):
        """The raw block without the key still yields the loader's answer.

        `science.yaml` states it, but a key that only works when present is a
        key whose absence silently reinstates the unusable default.
        """
        stripped = _config_without(config, "min_ireads")
        assert "min_ireads" not in (stripped.raw.get("variants") or {}).get(
            "bcftools", {}
        ).get("mpileup", {})
        command = _pileup_argv(stripped)
        assert command[command.index("-m") + 1] == str(
            DEFAULT_MPILEUP_MIN_IREADS
        )

    def test_call_settings_does_not_edit_the_shared_raw_config(
        self, config: PipelineConfig
    ):
        """`config.raw` is shared by every later stage in a run.

        Filling the default in place would be invisible here and would change
        the configuration for stages 6a onwards.
        """
        before = json.dumps(config.raw.get("variants"), sort_keys=True)
        stage.call_settings(config)
        assert json.dumps(config.raw.get("variants"), sort_keys=True) == before

    def test_the_adapter_reads_the_key_the_stage_wrote(self, config: PipelineConfig):
        """`call_isolate` and `call_settings` must name the same key.

        Two spellings of one setting is two settings; the stage writes
        `bcftools.mpileup.min_ireads` and the adapter reads it from there, so
        this asserts the path, not just the number.
        """
        settings = stage.call_settings(config)
        assert settings["bcftools"]["mpileup"]["min_ireads"] == 1

    def test_every_other_mpileup_key_still_comes_from_the_file(
        self, config: PipelineConfig
    ):
        settings = stage.call_settings(config)
        configured = config.raw["variants"]["bcftools"]["mpileup"]
        for key in ("min_mapq", "min_bq", "annotate", "max_depth"):
            assert settings["bcftools"]["mpileup"][key] == configured[key]


# ---------------------------------------------------------------------------
# A4: provenance, and the zero-indel warning.
# ---------------------------------------------------------------------------


class TestTheStageCountsWhatKindOfAlleleItCalled:
    def test_a_substitution_is_not_an_indel(self):
        assert stage.count_alleles([_row("A", "C")]) == {
            "n_snv": 1, "n_indel": 0, "n_alleles": 1
        }

    def test_a_deletion_is_an_indel(self):
        counts = stage.count_alleles([_row("CAT", "C")])
        assert (counts["n_snv"], counts["n_indel"]) == (0, 1)

    def test_an_insertion_is_an_indel(self):
        counts = stage.count_alleles([_row("C", "CT")])
        assert (counts["n_snv"], counts["n_indel"]) == (0, 1)

    def test_a_multi_base_substitution_is_not_reported_as_a_substitution(self):
        """`AC>GT` changes no length but is not one base for one base.

        Bucketing it with substitutions would understate the non-SNV content,
        which is the quantity this whole file is about.
        """
        counts = stage.count_alleles([_row("AC", "GT")])
        assert (counts["n_snv"], counts["n_indel"]) == (0, 1)

    def test_every_key_is_present_at_zero(self):
        """An absent count and a count of zero must not look alike."""
        assert stage.count_alleles([]) == {
            "n_snv": 0, "n_indel": 0, "n_alleles": 0
        }

    def test_the_counts_are_per_isolate(self):
        calls = {
            "S1": [_row("A", "C", "S1"), _row("AT", "A", "S1")],
            "S2": [_row("A", "C", "S2")],
            "S3": [],
        }
        payload = stage.call_provenance(calls)
        by_id = {e["sample_id"]: e for e in payload["per_isolate"]}
        assert by_id["S1"]["n_snv"] == 1 and by_id["S1"]["n_indel"] == 1
        assert by_id["S2"]["n_indel"] == 0
        assert by_id["S3"]["n_snv"] == 0 and by_id["S3"]["n_indel"] == 0
        assert payload["totals"] == {
            "isolates": 3, "n_calls": 3, "n_snv": 2, "n_indel": 1,
            "n_alleles": 3,
        }

    def test_the_totals_are_recomputed_not_accumulated(self):
        """Two counters free to drift are two numbers to distrust."""
        calls = {"S%d" % i: [_row("A", "C"), _row("A", "AT")] for i in range(4)}
        payload = stage.call_provenance(calls)
        assert payload["totals"]["n_alleles"] == sum(
            e["n_alleles"] for e in payload["per_isolate"]
        )

    def test_the_record_is_byte_stable_across_runs(self):
        calls = {"S2": [_row("A", "C", "S2")], "S1": [_row("A", "AT", "S1")]}
        assert stage.call_provenance(calls) == stage.call_provenance(calls)
        assert [e["sample_id"] for e in stage.call_provenance(calls)["per_isolate"]] == [
            "S1", "S2"
        ]


class TestTheSidecarIsWrittenBesideTheCalls:
    def test_it_lands_on_disk_and_reads_back(self, tmp_path: Path):
        calls = {"S1": [_row("A", "C"), _row("A", "AT")]}
        target = tmp_path / "nested" / stage.PROVENANCE_NAME
        payload = stage.write_provenance(calls, target)
        assert target.is_file()
        assert json.loads(target.read_text(encoding="utf-8")) == payload

    def test_an_unwritable_path_is_a_warning_not_a_stage_failure(
        self, tmp_path: Path, caplog
    ):
        """Losing the receipt must not lose the result."""
        blocker = tmp_path / "blocker"
        blocker.write_text("not a directory", encoding="utf-8")
        with caplog.at_level(logging.WARNING, logger="stages.variants"):
            payload = stage.write_provenance(
                {"S1": [_row("A", "C")]}, blocker / stage.PROVENANCE_NAME
            )
        assert payload["totals"]["n_snv"] == 1
        assert any("could not write" in r.message for r in caplog.records)


class TestARealRunWithNoIndelsSaysSo:
    """The warning, asserted on the message a reader would actually receive."""

    def _calls(self, *, n_snv: int, n_indel: int) -> Dict[str, List[Dict[str, str]]]:
        rows = [_row("A", "C") for _ in range(n_snv)]
        rows += [_row("AT", "A") for _ in range(n_indel)]
        return {"S1": rows, "S2": list(rows)}

    def test_it_warns_and_names_the_reason(
        self, tmp_path: Path, caplog
    ):
        with caplog.at_level(logging.WARNING, logger="stages.variants"):
            stage.summarise_call_provenance(
                self._calls(n_snv=500, n_indel=0),
                mode=RunMode.REAL,
                path=tmp_path / stage.PROVENANCE_NAME,
            )
        warnings = [r for r in caplog.records if r.levelno >= logging.WARNING]
        assert warnings, "a REAL run with substitutions and no indels was silent"
        assert any(
            stage.NO_INDEL_WARNING in r.getMessage() for r in warnings
        ), [r.getMessage() for r in warnings]

    def test_the_warning_names_min_ireads_because_that_is_the_cause(
        self, tmp_path: Path, caplog
    ):
        """A warning that says only "suspicious" makes the reader guess."""
        with caplog.at_level(logging.WARNING, logger="stages.variants"):
            stage.summarise_call_provenance(
                self._calls(n_snv=3, n_indel=0),
                mode=RunMode.REAL,
                path=tmp_path / stage.PROVENANCE_NAME,
            )
        text = " ".join(r.getMessage() for r in caplog.records)
        assert "min-ireads" in text

    def test_one_indel_anywhere_silences_it(self, tmp_path: Path, caplog):
        """The condition is the cohort total, not per isolate."""
        calls = self._calls(n_snv=10, n_indel=0)
        calls["S2"] = [_row("AT", "A", "S2")]
        with caplog.at_level(logging.WARNING, logger="stages.variants"):
            stage.summarise_call_provenance(
                calls, mode=RunMode.REAL, path=tmp_path / stage.PROVENANCE_NAME
            )
        assert not any(
            stage.NO_INDEL_WARNING in r.getMessage() for r in caplog.records
        )

    def test_test_mode_does_not_warn(self, tmp_path: Path, caplog):
        """A committed fixture's allele mix is a property of the fixture.

        The fixture generator is not obliged to emit indels, and warning about
        every TEST run would train a reader to ignore the warning.
        """
        with caplog.at_level(logging.WARNING, logger="stages.variants"):
            stage.summarise_call_provenance(
                self._calls(n_snv=5, n_indel=0), mode=RunMode.TEST
            )
        assert not any(
            stage.NO_INDEL_WARNING in r.getMessage() for r in caplog.records
        )

    def test_the_sidecar_is_written_on_the_real_path(self, tmp_path: Path):
        target = tmp_path / stage.PROVENANCE_NAME
        payload = stage.summarise_call_provenance(
            self._calls(n_snv=1, n_indel=1),
            mode=RunMode.REAL, path=target,
        )
        assert target.is_file()
        assert payload["totals"]["n_indel"] == 2

    def test_test_mode_writes_no_sidecar_when_no_path_is_given(
        self, tmp_path: Path
    ):
        """The TEST branch has no scratch workdir; it must not invent one."""
        before = set(tmp_path.iterdir())
        stage.summarise_call_provenance(
            self._calls(n_snv=1, n_indel=0), mode=RunMode.TEST
        )
        assert set(tmp_path.iterdir()) == before


# ---------------------------------------------------------------------------
# A5: the rows the new indels produce are accepted downstream.
# ---------------------------------------------------------------------------


class TestNonSnvRowsSurviveToTheMerge:
    """The real oprD rows from the A2 rerun, as input.

    Shapes taken from the run, not invented: a 2 nt net deletion written by
    bcftools left-aligned as `CAT>C`, and a 5 nt insertion written as a
    `GCCGGACCGGAC>GCCGGACCGGACCGGAC` repeat expansion. Both are the kind of row
    that a length-1 assumption would drop.
    """

    OPRD_ROWS = (
        ("1044210", "CAT", "C"),
        ("1044098", "GCCGGACCGGAC", "GCCGGACCGGACCGGAC"),
    )

    def _isolate(self) -> Dict[str, List[Dict[str, str]]]:
        rows = [_row("A", "C", "S1", "100")]
        rows += [_row(ref, alt, "S1", pos) for pos, ref, alt in self.OPRD_ROWS]
        return {"S1": rows}

    def test_the_parser_keeps_them(self):
        payload = stage.call_provenance(self._isolate())
        assert payload["totals"]["n_indel"] == 2

    def test_the_merge_counts_them_as_carriers(self, manifest):
        """Carried by a minority, so they read as polymorphic.

        `merge_calls` drops a site nearly every isolate carries: that is the
        fixed-difference rule, and giving all twenty isolates the same indel
        would exercise that rule instead of the indel path.
        """
        from papipeline.stages import cohort_variants

        carriers = set(manifest.sample_ids[:5])
        dissenters = set(manifest.sample_ids[5:8])
        calls: Dict[str, List[Any]] = {}
        for sample_id in manifest.sample_ids:
            tuples: List[Any] = []
            if sample_id in carriers:
                tuples += [
                    (r["chrom"], r["pos"], r["ref"], r["alt"])
                    for r in self._isolate()["S1"]
                ]
            if sample_id in dissenters:
                tuples.append(("NC_002516.2", "200", "A", "C"))
            calls[sample_id] = tuples
        merged = cohort_variants.merge_calls(calls)
        loci = {(r["pos"], r["ref"], r["alt"]) for r in merged}
        assert loci == {
            ("100", "A", "C"),
            ("200", "A", "C"),
            ("1044210", "CAT", "C"),
            ("1044098", "GCCGGACCGGAC", "GCCGGACCGGACCGGAC"),
        }
        for row in merged:
            if row["pos"] in {"1044210", "1044098"}:
                assert row["ac"] == "5" and row["an"] == str(len(manifest.sample_ids))

    def test_the_mutation_features_classify_them(self, config: PipelineConfig):
        """`gff.classify` is what the regulator screen reads each call through."""
        from papipeline.adapters import gff

        interval = gff.GeneInterval(
            locus_tag="PA0958", name="oprD", contig="NC_002516.2",
            start=1043983, end=1045314, strand="-",
        )
        deletion = gff.classify(
            gff.VariantCall("S1", "NC_002516.2", 1044210, "CAT", "C"), "", interval
        )
        insertion = gff.classify(
            gff.VariantCall("S1", "NC_002516.2", 1044098, "GCCGGACCGGAC",
                            "GCCGGACCGGACCGGAC"),
            "", interval,
        )
        # -2 nt and +5 nt are both frameshifts, and neither may be reported as
        # a substitution.
        assert deletion == "FRAMESHIFT"
        assert insertion == "FRAMESHIFT"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _installed_mpileup_usage() -> str:
    """The installed bcftools' own mpileup usage text.

    `bcftools mpileup --help` is not an option on 1.23.1 - it exits 1 with
    ``unrecognized option `--help'`` - so a bare `bcftools mpileup` is what
    prints the usage. Returns empty when bcftools is absent; the class is
    skipped in that case.
    """
    import subprocess

    if not shutil.which("bcftools"):
        return ""
    out = subprocess.run(
        ["bcftools", "mpileup"], capture_output=True, text=True, check=False
    )
    return (out.stdout or "") + (out.stderr or "")


def _config_without(config: PipelineConfig, key: str) -> PipelineConfig:
    """A copy of ``config`` with one mpileup key removed from ``raw``.

    Built by replacing ``raw`` wholesale rather than by mutating the dict in
    place, because mutating a frozen config's ``raw`` would leak into every
    other test in the session.
    """
    import copy

    raw = copy.deepcopy(dict(config.raw))
    raw.get("variants", {}).get("bcftools", {}).get("mpileup", {}).pop(key, None)
    return PipelineConfig(
        root=config.root,
        raw=raw,
        antibiotics=config.antibiotics,
        allowed_phenotypes=config.allowed_phenotypes,
        analysis=dict(config.analysis),
        mechanisms=dict(config.mechanisms),
        regulators=dict(config.regulators),
        references=dict(config.references),
        antibiotic_specs=dict(config.antibiotic_specs),
        organism=dict(config.organism),
        qc=config.qc, gwas=config.gwas, phylogeny=config.phylogeny,
        convergence=config.convergence,
        paths=dict(config.paths), runtime=dict(config.runtime),
    )