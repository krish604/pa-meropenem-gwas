"""Stage 7's REAL caller: panaroo's own table in, our partition out.

**What was missing.** `stages.pangenome.run()` had exactly one path. It built the
partition from the standardised annotation table, which is a TEST/STUB shape: in
that shape a "gene" is a literal annotation gene name. panaroo clusters orthologs
and names them `group_<n>`, so a REAL run that went through that path would have
reported gene *names* as if they were ortholog groups - every downstream feature
label resting on a correspondence the data does not contain. The adapter that
reads panaroo existed (`adapters.panaroo`) and nothing called it.

**What is wired here, and in what order.** The REAL branch is:

    gate (`runtime.allow_real_mode`) -> preflight() -> parse -> partition

`run_pipeline` then calls `write_outputs` on the result, so the chain is complete
without `run()` writing anything itself.

The order is load-bearing and each boundary is tested separately:

- The **gate comes first**. A shut gate must raise about `allow_real_mode`, not
  about a missing tool. If the preflight ran first, an operator whose gate is
  shut would be told to install panaroo - which is not the problem, and fixing it
  would not have helped.
- The **preflight comes before the parse**. panaroo shells out to `cd-hit`, so
  `panaroo` can be importable while the cohort cannot be processed at all. A
  valid presence table sitting on disk must not read as "this stage works".

**The gate is a gate, not a refusal.** `test_stage_taxonomy_is_self_verifying.py`
distinguishes the two from each stage module's AST: a raise guarded by
`allow_real_mode` is state 1 (runnable, gate shut) and must *not* appear in
`REAL_REFUSING_STAGES`. So this test also pins that `pangenome` is absent from
that set, which is what stops the wiring from quietly reclassifying a working
stage as a refusing one.

**No panaroo counts are asserted here.** The verified figures for the 10-isolate
smoke cohort are 10,019 families / 4,828 core / 5,191 accessory / 0
genome-specific, and they are real - but they belong to one dataset. A fixture
asserting them would fail the moment the cohort changed while teaching nothing
about the code; `test_pangenome_real_identifiers.py` says the same thing about
identifiers. What is asserted is that the metrics the summary reports are the
ones the partition computes.
"""

from __future__ import annotations

import copy
from pathlib import Path

import pytest

from papipeline.adapters import panaroo as adapter
from papipeline.config.loader import load_config
from papipeline.errors import (
    DataContractError,
    ModeNotAllowedError,
    ToolNotAvailableError,
)
from papipeline.manifest import SampleManifest
from papipeline.models import AnnotationRecord, RunMode, Sample
from papipeline.run import REAL_REFUSING_STAGES
from papipeline.stages import pangenome as stage

REPO = Path(__file__).resolve().parents[2]

#: The ten isolate ids of the smoke cohort, in the form panaroo names its
#: columns: the assembly file stem. Real ids, because the cohort-mismatch check
#: is a string comparison and made-up ids would not exercise it.
COHORT = tuple(
    f"PDT{n:09d}.1"
    for n in (34122, 167133, 167135, 641196, 730804, 876011, 902875, 907375, 949412, 988243)
)

#: Every isolate. Named so a "core" row reads as what it is.
ALL = frozenset(range(len(COHORT)))


def _manifest(sample_ids=COHORT) -> SampleManifest:
    return SampleManifest(
        samples=[Sample(sample_id=s, assembly_path=f"/nowhere/{s}.fna") for s in sample_ids]
    )


def _annotation(sample_id: str, gene_name: str) -> AnnotationRecord:
    """A minimal stage-2 record. Only `gene_name` matters to this stage."""
    return AnnotationRecord(
        sample_id=sample_id,
        contig_id="contig_1",
        gene_id=f"{gene_name}_{sample_id}",
        gene_name=gene_name,
        product=None,
        gene_type="CDS",
        start=1,
        end=300,
        strand="+",
        annotation_source="test",
    )


def _header(sample_ids=COHORT) -> str:
    return "Gene,Non-unique Gene name,Annotation," + ",".join(sample_ids)


def _presence_csv(rows, sample_ids=COHORT) -> str:
    return "\n".join([_header(sample_ids)] + rows) + "\n"


def _row(gene: str, carriers, annotation: str = "") -> str:
    """One panaroo presence row: gene, non-unique name, annotation, then genomes.

    `carriers` is a set of positions into ``COHORT``. Written from the cohort's
    width rather than by hand so a row can never be a comma short - a short row
    is read as a shorter cohort rather than as a malformed one, which is the same
    silent misalignment the adapter's `csv.reader` choice exists to prevent.
    """
    cells = ["1" if i in carriers else "0" for i in range(len(COHORT))]
    return ",".join([gene, "", annotation] + cells)


@pytest.fixture()
def shut():
    """The committed laptop overlay: `allow_real_mode` false."""
    config = load_config(REPO / "config" / "science.yaml", machine="laptop")
    assert config.runtime.get("allow_real_mode") is False
    return config


@pytest.fixture()
def open_gate(monkeypatch):
    """Open REAL through the real mechanism, not a patched config value.

    `PIPELINE_ALLOW_REAL_MODE` is the session override added for exactly this, so
    a test that bypassed it would not prove the stage reads what an operator
    sets. monkeypatch undoes it, and the other overlays stay shut.
    """
    monkeypatch.setenv("PIPELINE_ALLOW_REAL_MODE", "1")
    config = load_config(REPO / "config" / "science.yaml", machine="laptop")
    assert config.runtime.get("allow_real_mode") is True
    return config


@pytest.fixture()
def tools_present(monkeypatch):
    """Both tools present, so the preflight is not what is under test."""
    monkeypatch.setattr(adapter, "preflight", lambda: None)


@pytest.fixture()
def cd_hit_missing(monkeypatch):
    """panaroo importable, cd-hit absent - the case a single check would miss."""
    monkeypatch.setitem(__import__("sys").modules, "panaroo", object())
    monkeypatch.setattr(adapter.shutil, "which", lambda name: None)


def _write_panaroo(intermediate_root: Path, text: str) -> Path:
    """Put a panaroo-format table where the stage is expected to read it."""
    panaroo_dir = Path(intermediate_root) / adapter.OUTPUT_DIRNAME
    panaroo_dir.mkdir(parents=True, exist_ok=True)
    (panaroo_dir / adapter.PRESENCE_CSV).write_text(text, encoding="utf-8")
    return panaroo_dir


class TestTheGateComesFirst:
    def test_it_refuses_when_allow_real_mode_is_shut(self, shut, tmp_path):
        with pytest.raises(ModeNotAllowedError) as excinfo:
            stage.run(shut, _manifest(), RunMode.REAL, tmp_path)
        message = str(excinfo.value)
        assert "allow_real_mode" in message
        assert str(tmp_path) not in message, (
            "the message is about the gate, so it must not describe a path it "
            "never tried to read"
        )

    def test_a_shut_gate_says_so_rather_than_blaming_a_missing_tool(self, shut, tmp_path):
        """No panaroo table, no panaroo, no cd-hit - and it must still be the gate.

        If the preflight ran first the reader would be told to install panaroo.
        That is not the problem, and fixing it would not have unblocked anything.
        """
        monkeypatch = pytest.MonkeyPatch()
        try:
            monkeypatch.setattr(adapter.shutil, "which", lambda name: None)
            monkeypatch.delitem(__import__("sys").modules, "panaroo", raising=False)
            with pytest.raises(ModeNotAllowedError) as excinfo:
                stage.run(shut, _manifest(), RunMode.REAL, tmp_path)
        finally:
            monkeypatch.undo()
        assert "cd-hit" not in str(excinfo.value)

    def test_the_refusal_names_what_to_change(self, shut, tmp_path):
        """Both remedies, not just the name of the flag.

        The sibling assertion above proves the message mentions
        `allow_real_mode`; this proves it says how to act on that. AGENTS.md
        rule 7 requires exactly this of the sample-cap refusal - "a message that
        explains the limit and what to change" - and a refusal that only names the
        flag leaves an operator to grep the tree for the mechanism. Two remedies
        exist and both are real: the tracked overlay, and the session override the
        open-gate fixture above already uses.
        """
        with pytest.raises(ModeNotAllowedError) as excinfo:
            stage.run(shut, _manifest(), RunMode.REAL, tmp_path)
        message = str(excinfo.value)
        assert "allow_real_mode" in message
        assert "PIPELINE_ALLOW_REAL_MODE" in message, (
            "the session override is a supported way to open this gate, so a "
            "refusal that omits it sends the reader to edit a tracked file for "
            "something they can do with one environment variable"
        )
        assert shut.machine_name in message, (
            "the refusal must name the overlay that is shut, or the reader has "
            "three candidate files and no way to tell which one refused"
        )

    def test_it_is_a_gate_and_not_a_stage_refusal(self):
        """`pangenome` must stay out of `REAL_REFUSING_STAGES`.

        That set means "runnable in TEST, refuses REAL by design". A stage whose
        `run()` refuses only while `allow_real_mode` is false is runnable, and
        listing it there would be a false positive that invites the next reader
        to "fix" a stage that works.
        """
        assert "pangenome" not in REAL_REFUSING_STAGES


class TestThePreflightComesBeforeTheParse:
    def test_a_missing_cd_hit_stops_a_run_with_a_valid_table_on_disk(
        self, open_gate, cd_hit_missing, tmp_path
    ):
        """The table parses. The stage still refuses, because cd-hit is absent.

        This is the ordering that matters: a well-formed presence table on disk
        must not read as "this stage works" when the tool that produces it cannot
        run on this machine.
        """
        _write_panaroo(tmp_path, _presence_csv([_row("group_1375", ALL)]))
        with pytest.raises(ToolNotAvailableError) as excinfo:
            stage.run(open_gate, _manifest(), RunMode.REAL, tmp_path)
        assert "cd-hit" in str(excinfo.value)

    def test_the_preflight_runs_before_the_table_is_read(
        self, open_gate, cd_hit_missing, tmp_path
    ):
        """No table at all, and still a tool refusal rather than a data refusal.

        A `DataContractError` here would mean the parse ran first, which is the
        ordering bug this class exists to catch.
        """
        with pytest.raises(ToolNotAvailableError):
            stage.run(open_gate, _manifest(), RunMode.REAL, tmp_path)


class TestTheRealPathParsesAndPartitions:
    def test_it_partitions_panaroos_own_table(
        self, open_gate, tools_present, tmp_path
    ):
        rows = [
            # core: in all ten isolates
            _row("group_1375", ALL),
            # accessory: in two. The annotation carries a comma, quoted - the
            # case `csv.reader` exists for.
            _row("group_5706_ABC transporter", {0, 1}, "protein binding, ATP binding"),
            # genome-specific: in one
            _row("group_5289", {2}),
            # absent from every isolate, and still a row
            _row("group_9001", set()),
        ]
        _write_panaroo(tmp_path, _presence_csv(rows))

        result = stage.run(open_gate, _manifest(), RunMode.REAL, tmp_path)

        assert result.sample_ids == COHORT
        assert result.n_samples == len(COHORT)
        assert result.core == ("group_1375",)
        assert set(result.accessory) == {
            "group_5706_ABC transporter",
            "group_5289",
            "group_9001",
        }

    def test_real_identifiers_reach_the_matrix_unchanged(
        self, open_gate, tools_present, tmp_path
    ):
        """`gene` means an ortholog group id in REAL, not an annotation name.

        A normalisation on the way in would feed pyseer a feature matrix keyed on
        labels that are not orthologous to each other, and every count downstream
        would still look reasonable.
        """
        _write_panaroo(tmp_path, _presence_csv([_row("group_5706_ABC transporter", ALL)]))
        result = stage.run(open_gate, _manifest(), RunMode.REAL, tmp_path)
        assert result.genes == ("group_5706_ABC transporter",)

    def test_the_summary_metrics_are_the_ones_the_partition_computes(
        self, open_gate, tools_present, tmp_path
    ):
        """The three metrics a REAL run reports, read off the written summary.

        `n_genes_total`, `n_core_genes` and `n_accessory_genes` are what the
        verified 10-isolate figures (10,019 / 4,828 / 5,191) are read from. The
        counts are cohort facts and are not asserted here; what is asserted is
        that the file the run publishes carries exactly those three metrics and
        that they agree with the partition that produced them.
        """
        _write_panaroo(
            tmp_path,
            _presence_csv(
                [
                    _row("group_1", ALL),
                    _row("group_2", set(range(7))),
                    _row("group_3", {0}),
                ]
            ),
        )
        result = stage.run(open_gate, _manifest(), RunMode.REAL, tmp_path)
        written = stage.write_outputs(result, tmp_path / "out")

        summary = dict(
            line.split("\t")
            for line in Path(written["pangenome_summary"])
            .read_text(encoding="utf-8")
            .splitlines()[1:]
            if line.strip()
        )
        assert summary["n_samples"] == str(len(COHORT))
        assert summary["n_genes_total"] == "3"
        assert summary["n_core_genes"] == "1"
        assert summary["n_accessory_genes"] == "2"
        # The partition is exhaustive: core and accessory partition the total, so
        # a summary whose three numbers did not add up would be describing a
        # cohort other than this one.
        assert int(summary["n_core_genes"]) + int(summary["n_accessory_genes"]) == int(
            summary["n_genes_total"]
        )


class TestTheCohortIsStillEnforced:
    def test_a_missing_isolate_is_a_hard_failure(
        self, open_gate, tools_present, tmp_path
    ):
        """A pangenome over 9 of 10 isolates reports different core/accessory
        boundaries while every number still looks plausible."""
        short = COHORT[:-1]
        _write_panaroo(
            tmp_path,
            _presence_csv([_row("group_1", ALL)], sample_ids=short),
        )
        with pytest.raises(DataContractError) as excinfo:
            stage.run(open_gate, _manifest(), RunMode.REAL, tmp_path)
        assert COHORT[-1] in str(excinfo.value), "the missing isolate must be named"

    def test_an_unprepared_isolate_is_named_too(
        self, open_gate, tools_present, tmp_path
    ):
        _write_panaroo(tmp_path, _presence_csv([_row("group_1", ALL)]))
        short_manifest = _manifest(COHORT[:-1])
        with pytest.raises(DataContractError) as excinfo:
            stage.run(open_gate, short_manifest, RunMode.REAL, tmp_path)
        assert COHORT[-1] in str(excinfo.value)

    def test_a_cohort_larger_than_its_assemblies_refuses(
        self, open_gate, tools_present, tmp_path
    ):
        """The REAL path partitions over the *manifest*, not over what is on disk.

        A smoke overlay's cohort is deliberately the whole roster (967 isolates,
        every one a member) while only ten carry an assembly, so panaroo's table
        has ten columns and the manifest names 967. That combination refuses,
        naming the isolates that were never clustered.

        It has to refuse. A core/accessory boundary computed over ten of a
        cohort's members is a real partition of a *different* cohort, and
        reporting it as this cohort's pangenome would be a false claim that every
        number in it supports. Narrowing the cohort here instead would be a
        scientific decision - it changes the denominator the partition is
        reported over - and it belongs to the overlay, not to this stage.
        """
        roster = _manifest(COHORT + ("PDT999999999.1",))
        _write_panaroo(tmp_path, _presence_csv([_row("group_1", ALL)]))
        with pytest.raises(DataContractError) as excinfo:
            stage.run(open_gate, roster, RunMode.REAL, tmp_path)
        assert "PDT999999999.1" in str(excinfo.value)


class TestBothModesPublishOneFormat:
    """TEST and REAL must publish the *same* format, and it must survive a re-read.

    **The gap this closes.** `test_pangenome_gpa_roundtrip.py` round-trips a
    `Pangenome` built by calling `partition()` directly. It never calls `run()`,
    so it says nothing about either mode's actual output - and the REAL branch
    added a second way for a `Pangenome` to be born. Two producers, one declared
    format, no test spanning both, is precisely the shape that leaves a consumer
    reading TEST fixtures while REAL writes something else.

    **What "round-trips" has to mean here.** Not that the object equals itself:
    `partition()` is pure and would pass that trivially. The chain is
    `run()` -> `write_outputs` -> `read_gene_presence_absence` -> `partition`,
    and the claim is that the *file* still carries enough to rebuild the same
    partition. That is the only version of the claim a downstream consumer
    depends on, and the only one that notices a writer and reader disagreeing.

    **Why format-identity across modes is the load-bearing assertion.** In TEST a
    `gene` is an annotation gene name; in REAL it is `group_<n>`. Same column,
    different meaning. A consumer branching on mode is the bug this forbids, and
    it is invisible in a per-mode test because each mode is individually fine.
    """

    def _round_trip(self, result, out_dir):
        written = stage.write_outputs(result, out_dir)
        samples, presence = stage.read_gene_presence_absence(
            written["gene_presence_absence"]
        )
        return written, stage.partition(samples, presence)

    def _test_mode_pangenome(self, config, tmp_path):
        annotations = {
            "S1": [_annotation("S1", "OprD"), _annotation("S1", "abcA")],
            "S2": [_annotation("S2", "OprD")],
        }
        return stage.run(
            config,
            _manifest(("S1", "S2")),
            RunMode.TEST,
            tmp_path,
            annotations,
        )

    def test_test_mode_output_round_trips(self, shut, tmp_path):
        result = self._test_mode_pangenome(shut, tmp_path)
        # A gene no isolate carries. It must survive as a present-but-empty entry.
        #
        # This row is what gives the round trip teeth. A writer that emits only
        # `present == 1` rows loses exactly these genes and nothing else, because
        # absence is implicit in a sparse matrix: every surviving row still reads
        # back correctly and every carrier-bearing gene still round-trips. Asserted
        # over carrier-bearing genes alone, this suite passes on a writer that
        # silently drops the genes nothing was found for - checked by mutation,
        # not assumed.
        result = stage.partition(
            result.sample_ids, {**result.presence, "gene_absent_from_all": set()}
        )
        written, rebuilt = self._round_trip(result, tmp_path / "out")
        assert rebuilt.sample_ids == result.sample_ids
        assert rebuilt.genes == result.genes
        assert rebuilt.core == result.core
        assert rebuilt.accessory == result.accessory
        assert rebuilt.core == ("OprD",)
        assert rebuilt.presence["gene_absent_from_all"] == set(), (
            "a gene no isolate carries vanished across the round trip; it must "
            "survive as present-but-empty, or 'absent' and 'never looked for' "
            "become indistinguishable"
        )

    def test_real_mode_output_round_trips(self, open_gate, tools_present, tmp_path):
        _write_panaroo(
            tmp_path,
            _presence_csv(
                [_row("group_1", ALL), _row("group_2", {0, 1}), _row("group_3", set())]
            ),
        )
        result = stage.run(open_gate, _manifest(), RunMode.REAL, tmp_path)
        written, rebuilt = self._round_trip(result, tmp_path / "out")
        assert rebuilt.sample_ids == result.sample_ids
        assert rebuilt.genes == result.genes
        assert rebuilt.core == result.core
        assert rebuilt.accessory == result.accessory
        assert rebuilt.presence["group_3"] == set(), (
            "panaroo emitted a family row no isolate carries; it must survive "
            "the round trip as present-but-empty"
        )

    def test_a_real_gene_name_with_spaces_survives(self, open_gate, tools_present, tmp_path):
        """The realistic panaroo name, which is not a bare token.

        panaroo labels a family `group_<n>_<product>`, so spaces are the normal
        case rather than an exotic one. A reader that splits or strips on
        whitespace would silently rename every accessory gene, and the pyseer
        feature matrix D6 builds would be keyed on labels that no longer match
        the table it came from.
        """
        name = "group_5706_ABC transporter"
        _write_panaroo(tmp_path, _presence_csv([_row(name, ALL)]))
        result = stage.run(open_gate, _manifest(), RunMode.REAL, tmp_path)
        _, rebuilt = self._round_trip(result, tmp_path / "out")
        assert rebuilt.genes == (name,), "the family name was altered by the format"
        assert name in rebuilt.presence

    def test_both_modes_write_the_same_header(self, shut, open_gate, tools_present, tmp_path):
        """One format, byte-identical, so no consumer has to know the mode.

        The two modes disagree about what a `gene` *means* and must not disagree
        about the file. A REAL-only column, or a reordered one, would force every
        reader to branch on mode - and a reader that branches on mode is how a
        TEST-validated consumer meets a REAL table for the first time in
        production.
        """
        _write_panaroo(tmp_path, _presence_csv([_row("group_1", ALL)]))
        real = stage.run(open_gate, _manifest(), RunMode.REAL, tmp_path)
        test = self._test_mode_pangenome(shut, tmp_path)

        real_written = stage.write_outputs(real, tmp_path / "real")
        test_written = stage.write_outputs(test, tmp_path / "test")

        def header(path):
            return [
                line
                for line in path.read_text(encoding="utf-8").splitlines()
                if line and not line.startswith("#")
            ][0]

        assert header(real_written["gene_presence_absence"]) == header(
            test_written["gene_presence_absence"]
        )
        assert header(real_written["gene_presence_absence"]) == "\t".join(stage.GPA_COLUMNS)

    def test_every_declared_output_is_written_in_both_modes(
        self, shut, open_gate, tools_present, tmp_path
    ):
        """All four declared files, from both paths.

        A mode that skipped `write_outputs` would leave a downstream stage reading
        a previous run's table - stale but well-formed, which is the worst
        outcome available: no refusal, just the wrong cohort's numbers.
        """
        _write_panaroo(tmp_path, _presence_csv([_row("group_1", ALL)]))
        for result, sub in (
            (self._test_mode_pangenome(shut, tmp_path), "test"),
            (stage.run(open_gate, _manifest(), RunMode.REAL, tmp_path), "real"),
        ):
            written = stage.write_outputs(result, tmp_path / sub)
            assert set(written) == {
                "gene_presence_absence",
                "core_genes",
                "accessory_genes",
                "pangenome_summary",
            }, f"{sub} mode published a different set of outputs"
            for name, path in written.items():
                assert path.exists(), f"{sub} mode did not write {name}"


class TestWhereTheTableIsReadFrom:
    def test_the_intermediate_root_decides_which_table_is_read(
        self, open_gate, tools_present, tmp_path
    ):
        """The location comes from the run's intermediate root, not a constant.

        Two roots, one carrying the table. Reading the wrong one either fails or
        picks up a previous run's cohort - and both are silent if the path
        happens to resolve.
        """
        populated = tmp_path / "run-a"
        empty = tmp_path / "run-b"
        _write_panaroo(populated, _presence_csv([_row("group_1", ALL)]))
        empty.mkdir()

        assert stage.run(open_gate, _manifest(), RunMode.REAL, populated).n_samples == 10
        with pytest.raises(DataContractError):
            stage.run(open_gate, _manifest(), RunMode.REAL, empty)

    def test_the_layout_constant_is_the_adapters(self):
        """One declaration, so a future runner and this reader cannot disagree.

        A runner that wrote to `intermediate/panaroo_out` and a reader looking in
        `intermediate/panaroo` would both be reasonable and jointly broken, and
        the only symptom would be a stage that refuses for want of a file that
        was written moments earlier.
        """
        assert adapter.OUTPUT_DIRNAME == "panaroo"
        assert adapter.output_dir(Path("/tmp/x")) == Path("/tmp/x/panaroo")


class TestTestModeIsUntouched:
    def test_test_mode_still_builds_from_annotations(self, shut, tmp_path):
        """The gate is on the REAL path only.

        `allow_real_mode` is false in every committed overlay, so a gate that
        leaked into TEST would make the whole test suite unrunnable.
        """
        annotations = {
            "S1": [
                _annotation("S1", "abcA"),
                _annotation("S1", "OprD"),
            ],
            "S2": [_annotation("S2", "OprD")],
        }
        result = stage.run(
            shut, _manifest(("S1", "S2")), RunMode.TEST, tmp_path, annotations
        )
        assert result.core == ("OprD",), "abcA is in 1 of 2, so it is not core"
        assert set(result.accessory) == {"abcA"}

    def test_test_mode_never_calls_the_preflight(self, shut, tmp_path, monkeypatch):
        """TEST invokes no tools at all, so a preflight here would refuse a run
        that needs no tool - on a machine where no tool is installed."""

        def explode() -> None:
            raise AssertionError("TEST mode must not call the panaroo preflight")

        monkeypatch.setattr(adapter, "preflight", explode)
        result = stage.run(
            shut,
            _manifest(("S1",)),
            RunMode.TEST,
            tmp_path,
            {"S1": [_annotation("S1", "abcA")]},
        )
        assert result.genes == ("abcA",)

    def test_stub_fabricates_the_stage_itself(self, tmp_path):
        """STUB never reaches `stages.pangenome.run` - `run_pipeline` fabricates.

        Asserted here only because the gate above must not change that: if STUB
        ever fell through to `run()`, it would be reading a cohort in a mode
        whose whole purpose is to read none.
        """
        from papipeline import stub

        fabricated = stub.fabricate("pangenome", tmp_path)
        assert fabricated, "STUB must fabricate stage 7 without calling the stage"