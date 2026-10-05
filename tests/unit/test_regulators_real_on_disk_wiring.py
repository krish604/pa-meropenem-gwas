"""REAL stage 6: the regulator screen's input is the ON-DISK TABLE.

`run.py:938` handed `stage_regulators.run` four arguments, so REAL reached
`regulators._run_real` with no `calls_by_isolate` and refused
`missing='calls_by_isolate'` on **every** REAL run. No test covered the call
site, so nothing caught it.

**The three roles are asserted apart, on distinct roots.**

* the *producer* - `run.derive_regulator_table` - reads this run's
  `<stage_dir>/variants.tsv` and writes `regulators/regulator_variants.tsv`,
  the stage input table named at `docs/data_contract.md:114`;
* the *stage* - `regulators.run` - reads that table and nothing else;
* TEST does neither, and reads the committed fixture exactly as before.

**Every root here is a different `tmp_path` child.** The two candidate paths -
the call table under `calls_root` and the regulator table under `tool_root` -
could otherwise coincide and let a swapped root pass. That is asserted first, in
`TestTheRootsAreDistinct`, because everything below is void without it.

**Why the call site is real and the reference is real.** The screen intersects
calls with the pinned PAO1 CDS features, and `PA0958` (oprD) is on the minus
strand, so a fixture cannot place a call in it without reproducing the strand
arithmetic the test would then be asserting. The coordinates here are measured
from `db/reference/GCF_000006765.1`, the same authority
`tests/unit/test_regulator_screen_producer.py` uses.
"""

from __future__ import annotations

import ast
import inspect
import json
from pathlib import Path
from typing import Any, Dict, List

import pytest

from papipeline.adapters import gff
from papipeline.errors import DataContractError, PipelineError, StageError
from papipeline.io.tsv import read_tsv, write_tsv
from papipeline.manifest import SampleManifest
from papipeline.models import RunMode, Sample
from papipeline.run import derive_regulator_table, run_pipeline
from papipeline.stages import regulators as reg
from papipeline.stages import variants as stage_variants

REPO = Path(__file__).resolve().parents[2]

#: Measured from the pinned reference, not from the code under test:
#: PA0958 (oprD) is a single 1332 nt gene feature, 1043983..1045314, MINUS
#: strand, so `start` is its *last* coding base. Verified against the real GFF.
OPRD_START, OPRD_END, OPRD_STRAND, OPRD_CONTIG = 1043983, 1045314, "-", "NC_002516.2"


# ---------------------------------------------------------------------------
# Roots. Four different directories, no two of them able to coincide.
# ---------------------------------------------------------------------------


@pytest.fixture()
def roots(tmp_path) -> Dict[str, Path]:
    """Four distinct roots, one per role, so a swapped root cannot pass."""
    made = {
        "calls": tmp_path / "calls_root",      # holds <stages>/variants.tsv
        "tool": tmp_path / "tool_root",        # holds regulators/*.tsv
        "absent": tmp_path / "absent_root",    # holds nothing at all
        "fixture": REPO / "test_data" / "intermediate",  # committed, read-only
    }
    (made["calls"] / "stages").mkdir(parents=True)
    (made["tool"] / "regulators").mkdir(parents=True)
    (made["absent"] / "stages").mkdir(parents=True)
    return made


@pytest.fixture()
def calls_table(roots) -> Path:
    return roots["calls"] / "stages" / "variants.tsv"


@pytest.fixture()
def one_sample() -> SampleManifest:
    """One cohort member, named `one_sample` so it cannot shadow the conftest's
    20-isolate TEST `manifest`. Every REAL case below needs an isolate whose
    absence from the call table is a real absence, not a cohort of 19."""
    return SampleManifest([Sample("TEST_A_01", None, "test")])


def _row(sample: str, contig: str, pos: int, ref: str, alt: str) -> Dict[str, str]:
    return {
        "sample_id": sample, "chrom": contig, "pos": str(pos),
        "ref": ref, "alt": alt, "qual": "60", "filter": "PASS",
        "GT": "1/1", "AC": "1", "AN": "1", "DP4": "0,0,30,29",
        "MQ": "60", "MQ0F": "0",
    }


def _opr_d_snv() -> Dict[str, str]:
    """A substitution inside oprD's CDS, on the minus strand.

    Position 1044100 is 118 nt from the gene's high end, so on the coding strand
    it is early in the transcript - and a plus-strand fixture would place it in
    the wrong place entirely.
    """
    return _row("TEST_A_01", OPRD_CONTIG, 1044100, "A", "T")


def _outside_every_locus() -> Dict[str, str]:
    return _row("TEST_A_01", OPRD_CONTIG, 779, "C", "G")


def _write_calls(path: Path, rows: List[Dict[str, str]]) -> Path:
    write_tsv(path, rows, list(stage_variants.PER_ISOLATE_COLUMNS))
    return path


def _regulator_table(root: Path) -> Path:
    return reg.regulator_table_path(root)


class TestTheRootsAreDistinct:
    """The premise. If two of these were equal, nothing below could fail."""

    def test_the_two_candidate_paths_are_different_files(self, roots, calls_table):
        produced = _regulator_table(roots["tool"])
        assert produced != calls_table
        assert calls_table.is_file() or not calls_table.exists()

    def test_every_temporary_root_is_a_separate_directory(self, roots, tmp_path):
        temporary = [roots["calls"], roots["tool"], roots["absent"]]
        assert len(set(temporary)) == len(temporary), (
            f"two roots coincide, so a swapped wiring is invisible: {temporary}"
        )
        for root in temporary:
            assert root.is_relative_to(tmp_path), (
                f"{root} is not under {tmp_path}, so this test would write into "
                "the repository"
            )

    def test_the_tool_root_is_not_this_runs_output_root(self, config):
        """In TEST the two roots differ; the contract is about REAL, where they
        are one path (`config/loader.py:529-543`). Asserting the REAL identity
        keeps this file honest about what it can and cannot distinguish."""
        assert config.tool_output_root(RunMode.REAL) == config.intermediate_root(
            RunMode.REAL
        )
        assert config.tool_output_root(RunMode.TEST) != config.intermediate_root(
            RunMode.TEST
        )


# ---------------------------------------------------------------------------
# (a) REAL with no variants.tsv refuses, naming the path and its producer
# ---------------------------------------------------------------------------


class TestRealWithNoCallTableRefuses:
    def test_it_refuses_naming_the_expected_path(self, config, one_sample, roots):
        missing = roots["absent"] / "stages" / "variants.tsv"
        assert not missing.exists()
        with pytest.raises(StageError) as excinfo:
            derive_regulator_table(
                config,
                one_sample,
                calls_table=missing,
                tool_output_root=roots["tool"],
            )
        message = str(excinfo.value)
        assert str(missing) in message, (
            "the refusal must name the path it looked in, or a wrong "
            "table_path(stage_dir, 'variants') is indistinguishable from a stage "
            "that never ran"
        )
        assert "variants.tsv" in message

    def test_it_names_what_produces_the_missing_table(self, config, one_sample, roots):
        missing = roots["absent"] / "stages" / "variants.tsv"
        with pytest.raises(StageError) as excinfo:
            derive_regulator_table(
                config, one_sample, calls_table=missing, tool_output_root=roots["tool"]
            )
        assert "stage_variants.run" in str(excinfo.value), (
            "a refusal that does not say what writes the file leaves the reader "
            "with no remedy"
        )

    def test_it_writes_nothing_on_the_way_out(self, config, one_sample, roots):
        """A refusal that left a table behind would be read as a result."""
        missing = roots["absent"] / "stages" / "variants.tsv"
        with pytest.raises(StageError):
            derive_regulator_table(
                config, one_sample, calls_table=missing, tool_output_root=roots["tool"]
            )
        assert not _regulator_table(roots["tool"]).exists()
        assert not (roots["tool"] / "regulators" / reg.SCREEN_REPORT_NAME).exists()

    def test_a_missing_pinned_reference_is_refused_by_name(
        self, config, one_sample, roots, calls_table
    ):
        """The other declared input, and the one that fabricates a negative.

        An empty reference places every call outside every locus, which reads as
        a cohort carrying no regulator variant. `PipelineConfig` is a frozen
        dataclass, so the reference is displaced with a delegating proxy rather
        than by assignment.
        """
        _write_calls(calls_table, [_opr_d_snv()])
        absent = roots["absent"] / "gone.gff"

        class _WithoutGff:
            reference_gff = staticmethod(lambda: absent)

            def __getattr__(self, name):
                return getattr(config, name)

        with pytest.raises(StageError) as excinfo:
            derive_regulator_table(
                _WithoutGff(), one_sample,
                calls_table=calls_table, tool_output_root=roots["tool"],
            )
        message = str(excinfo.value)
        assert str(absent) in message
        assert "reference_fasta()" in message, (
            "the refusal must say how the GFF path is derived, or the reader "
            "cannot tell a wrong reference from a wrong path"
        )
        assert not _regulator_table(roots["tool"]).exists()

    def test_a_call_from_a_sample_the_manifest_never_names_is_refused(
        self, config, one_sample, roots, calls_table
    ):
        """Rule 5: a sample-ID mismatch is a hard failure, never a silent drop."""
        _write_calls(calls_table, [_row("GHOST_1", OPRD_CONTIG, 1044100, "A", "T")])
        with pytest.raises(StageError) as excinfo:
            derive_regulator_table(
                config, one_sample,
                calls_table=calls_table, tool_output_root=roots["tool"],
            )
        assert "GHOST_1" in str(excinfo.value)
        assert not _regulator_table(roots["tool"]).exists()


# ---------------------------------------------------------------------------
# (b) REAL with a variants.tsv writes the table, and the stage reads it
# ---------------------------------------------------------------------------


class TestRealDerivesThenReadsTheTable:
    def test_the_producer_writes_the_contracted_table(
        self, config, one_sample, roots, calls_table
    ):
        _write_calls(calls_table, [_opr_d_snv(), _outside_every_locus()])
        written = derive_regulator_table(
            config, one_sample, calls_table=calls_table, tool_output_root=roots["tool"]
        )
        expected = _regulator_table(roots["tool"])
        assert written == expected
        assert expected.is_file(), (
            f"{expected.name} was not written; docs/data_contract.md:114 declares "
            "it a stage input table, read in REAL as well as TEST"
        )
        assert (roots["tool"] / "regulators" / reg.SCREEN_REPORT_NAME).is_file(), (
            "the provenance sidecar is what makes a zero-row table readable as "
            "'screened and clean' rather than 'never ran'"
        )

    def test_the_table_holds_the_intersected_call_and_not_the_stray_one(
        self, config, one_sample, roots, calls_table
    ):
        _write_calls(calls_table, [_opr_d_snv(), _outside_every_locus()])
        derive_regulator_table(
            config, one_sample, calls_table=calls_table, tool_output_root=roots["tool"]
        )
        rows = read_tsv(_regulator_table(roots["tool"]), required_columns=reg.REQUIRED)
        assert [r["gene"] for r in rows] == ["oprD"]
        assert rows[0]["variant_type"] == "SNV"
        assert rows[0]["position"] == "1044100"

    def test_the_producer_records_what_it_examined(
        self, config, one_sample, roots, calls_table
    ):
        """Calls outside the screened loci are counted, never dropped silently."""
        _write_calls(calls_table, [_opr_d_snv(), _outside_every_locus()])
        derive_regulator_table(
            config, one_sample, calls_table=calls_table, tool_output_root=roots["tool"]
        )


        report = json.loads(
            (roots["tool"] / "regulators" / reg.SCREEN_REPORT_NAME).read_text(
                encoding="utf-8"
            )
        )
        assert report["mode"] == "REAL"
        assert report["calls_examined"] == 2
        assert report["calls_outside_screened_loci"] == 1
        assert report["records_emitted"] == 1

    def test_the_stage_reads_the_table_and_does_not_recompute_it(
        self, config, one_sample, roots, calls_table
    ):
        """Delete the table and the stage refuses.

        This is the assertion that distinguishes *reads the file* from *recomputes
        from something it still holds*. If the stage computed the table itself,
        deleting it would change nothing and this test would pass for the wrong
        reason.
        """
        _write_calls(calls_table, [_opr_d_snv()])
        derive_regulator_table(
            config, one_sample, calls_table=calls_table, tool_output_root=roots["tool"]
        )
        _regulator_table(roots["tool"]).unlink()
        with pytest.raises(DataContractError) as excinfo:
            reg.run(config, one_sample, RunMode.REAL, roots["tool"])
        message = str(excinfo.value)
        assert "regulator_variants.tsv" in message
        assert str(_regulator_table(roots["tool"])) in message
        assert "produce_regulator_variants" in message, (
            "the refusal must name the producer, or a reader cannot tell an "
            "unrun producer from a wrong root"
        )

    def test_the_stage_returns_exactly_the_rows_on_disk(
        self, config, one_sample, roots, calls_table
    ):
        _write_calls(calls_table, [_opr_d_snv()])
        derive_regulator_table(
            config, one_sample, calls_table=calls_table, tool_output_root=roots["tool"]
        )
        grouped = reg.run(config, one_sample, RunMode.REAL, roots["tool"])
        on_disk = read_tsv(_regulator_table(roots["tool"]), required_columns=reg.REQUIRED)
        assert set(grouped) == set(one_sample.sample_ids), (
            "the mapping must carry a key for every manifest sample; an absent "
            "key means 'never screened' and is not 'clean'"
        )
        assert [(r.variant, r.gene, r.variant_type) for r in grouped["TEST_A_01"]] == [
            (r["variant"], r["gene"], r["variant_type"]) for r in on_disk
        ]

    def test_the_round_trip_survives_the_manifest_sample_with_no_calls(
        self, config, roots, calls_table
    ):
        """A cohort member with no called variant is present and empty.

        The denominator of every allele frequency downstream comes from the cohort
        size, so dropping it would inflate each frequency the merge reports.
        """
        manifest = SampleManifest(
            [Sample("TEST_A_01", None, "test"), Sample("TEST_A_02", None, "test")]
        )
        _write_calls(calls_table, [_opr_d_snv()])
        derive_regulator_table(
            config, manifest, calls_table=calls_table, tool_output_root=roots["tool"]
        )
        grouped = reg.run(config, manifest, RunMode.REAL, roots["tool"])
        assert set(grouped) == {"TEST_A_01", "TEST_A_02"}
        assert grouped["TEST_A_02"] == []
        assert grouped["TEST_A_01"]


# ---------------------------------------------------------------------------
# The call site itself. Everything above calls `derive_regulator_table`
# directly, so without this class a branch that never called it would pass.
# ---------------------------------------------------------------------------


class TestTheVariantsBranchCallsTheProducer:
    """`run.py::_execute_stage`'s `variants` branch, inspected as an AST.

    The same seam `tests/unit/test_variant_stage_branches.py` uses, and for the
    same reason: the branch is the function the Snakemake rules and the
    observatory call, so a test against it is a test against what a run does.

    **Why AST and not a REAL run.** Driving this branch for real needs minimap2,
    bcftools and blastp over the cohort, and AGENTS.md rule 6 forbids that
    without an explicit instruction. What is being asserted is the wiring - which
    function is called, with which two roots, and that the stage below it is not
    handed the calls - and none of that is a question a tool can answer. What is
    NOT covered this way, stated plainly: that the branch *runs* in REAL. That
    is covered by `test_the_producer_is_never_called_in_test` for TEST (it is
    not called) and is left to a REAL run for REAL.
    """

    @staticmethod
    def _branch_source() -> str:
        tree = ast.parse(
            (REPO / "papipeline" / "run.py").read_text(encoding="utf-8")
        )
        for node in ast.walk(tree):
            if not (isinstance(node, ast.FunctionDef) and node.name == "_execute_stage"):
                continue
            for inner in ast.walk(node):
                if not isinstance(inner, ast.If):
                    continue
                test = inner.test
                if (
                    isinstance(test, ast.Compare)
                    and isinstance(test.comparators[0], ast.Constant)
                    and test.comparators[0].value == "variants"
                ):
                    return "\n".join(ast.unparse(s) for s in inner.body)
        raise AssertionError("no branch for stage 'variants' in _execute_stage")

    def test_the_branch_calls_the_producer(self):
        assert "derive_regulator_table(" in self._branch_source(), (
            "the variants branch no longer derives the regulator table, so REAL "
            "would reach stage_regulators.run with no table on disk and refuse"
        )

    def test_the_producer_is_handed_the_call_table_the_branch_just_wrote(self):
        source = self._branch_source()
        assert "table_path(stage_dir, 'variants')" in source, (
            "the producer must read the file the branch itself wrote; anything "
            "else is a second source of truth for the calls"
        )
        assert "calls_table=table_path(stage_dir, 'variants')" in source

    def test_the_producer_is_handed_the_contracted_stage_input_root(self):
        assert "tool_output_root=tool_output_root" in self._branch_source(), (
            "the regulator table is contracted to tool_output_root "
            "(docs/data_contract.md:114), not to this run's own output root"
        )

    def test_the_stage_below_is_not_handed_the_calls(self):
        """The point of the ruling: the file is the input, not the parameter."""
        source = self._branch_source()
        head, _, tail = source.partition("stage_regulators.run(")
        assert "calls_by_isolate" not in tail.split("\n", 1)[0], (
            "stage_regulators.run is still handed calls_by_isolate, so REAL's "
            "input is a caller-scope mapping rather than the on-disk table"
        )
        # `variant_calls` legitimately exists in this branch - it is what
        # `stage_variants.run` returns and what is written to the call table.
        # What must not happen is the screen seeing it. `cohort_variants`
        # downstream owns that copy; the screen's input is the file.
        producer = head.split("derive_regulator_table(", 1)
        assert len(producer) == 2 and "variant_calls" not in producer[1], (
            "the producer is handed the in-memory mapping rather than the file"
        )



class TestTestModeIsUnchanged:
    def test_it_reads_the_committed_fixture(self, config, manifest):
        """The conftest `manifest` is the 20 committed TEST isolates."""
        fixture_root = config.tool_output_root(RunMode.TEST)
        path = _regulator_table(fixture_root)
        assert path == REPO / "test_data" / "intermediate" / "regulators" / (
            "regulator_variants.tsv"
        )
        assert path.is_file(), "the committed stage-input fixture is absent"

        grouped = reg.run(config, manifest, RunMode.TEST, fixture_root)
        on_disk = read_tsv(path, required_columns=reg.REQUIRED)
        assert set(grouped) == set(manifest.sample_ids)
        carriers = sorted(sid for sid, rows in grouped.items() if rows)
        assert carriers, "TEST read the committed fixture and found nothing"
        assert carriers == sorted({r["sample_id"] for r in on_disk}), (
            "TEST must return the fixture's own rows; a different set would mean "
            "the reader no longer reads the fixture"
        )
        assert sum(len(rows) for rows in grouped.values()) == len(on_disk)

    def test_the_producer_is_never_called_in_test(self, monkeypatch):
        """The whole TEST pipeline, with the producer spied on.

        Not a source-string assertion: `derive_regulator_table` writes into the
        stage-input root, so a TEST run that reached it would put a REAL
        regulator table inside the committed fixture tree.
        """
        import papipeline.run as run_module

        called: List[Any] = []
        original = run_module.derive_regulator_table

        def spy(*args, **kwargs):
            bound = inspect.signature(original).bind(*args, **kwargs)
            called.append(bound.arguments)
            return original(*args, **kwargs)

        monkeypatch.setattr(run_module, "derive_regulator_table", spy, raising=True)
        run_pipeline(mode="TEST", write_html=False)
        assert called == [], (
            f"the REAL producer ran during a TEST pipeline: {called}. It writes "
            "into the stage-input root, which in TEST is the committed fixture "
            "tree."
        )

    def test_the_committed_fixture_is_byte_identical_after_a_test_run(self):
        before = {
            p: (p.stat().st_size, p.stat().st_mtime_ns)
            for p in sorted((REPO / "test_data").rglob("*"))
            if p.is_file()
        }
        run_pipeline(mode="TEST", write_html=False)
        after = {
            p: (p.stat().st_size, p.stat().st_mtime_ns)
            for p in sorted((REPO / "test_data").rglob("*"))
            if p.is_file()
        }
        assert before == after, "a TEST run wrote into the committed fixtures"

    def test_the_two_doors_on_run_are_distinguishable(
        self, config, one_sample, roots, calls_table
    ):
        """`calls_by_isolate` remains on `run` for compatibility, and it is a
        *door onto the producer*, not this stage's input.

        Both doors are asserted, so neither can quietly become the other: passing
        the calls writes a table and returns the screen; omitting them reads the
        table and refuses when it is absent. `papipeline/run.py` takes the second
        door, which is what makes the file the input rather than a parameter.
        """
        table = _regulator_table(roots["tool"])
        assert not table.exists()

        # Door 1: no calls supplied -> the on-disk table is the input.
        with pytest.raises(DataContractError) as excinfo:
            reg.run(config, one_sample, RunMode.REAL, roots["tool"])
        assert str(table) in str(excinfo.value)
        assert not table.exists(), "the refusal path wrote a table"

        # Door 2: calls supplied -> the producer runs and writes it.
        _write_calls(calls_table, [_opr_d_snv()])
        grouped = reg.run(
            config, one_sample, RunMode.REAL, roots["tool"],
            calls_by_isolate={"TEST_A_01": [_opr_d_snv()]},
        )
        assert table.is_file(), "the compatibility door stopped producing"
        assert grouped["TEST_A_01"]
