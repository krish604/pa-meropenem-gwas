"""A stand-in fixture must never be read as stage output.

The stand-ins under `test_data/standins/unbuilt_stages/` exist so the DAG can
declare the edges spec.md requires while two stages have no implementation
(ticket 14). The danger is not their existence - it is one being read as though
it were a result, because a contract-shaped table with plausible columns is what
real output looks like too.

Two things are checked here:

  - :func:`papipeline.standins.load_stage_table` refuses, by banner and by
    value, so a caller cannot consume a stand-in by accident;
  - *no other module reads those tables at all*. That is the structural half.
    Ticket 14 adds real consumers, and this assertion is written to fail when
    they appear, which is the point - see the note in its docstring.

**Which stages this covers is derived, not written down.**
:data:`TABLES` below is computed from
:data:`papipeline.standins.UNBUILT_WITH_STANDIN`, so narrowing that tuple
narrows every test here at once, and a test cannot quietly keep guarding a
stage that now has a real implementation. `variants` left the tuple in f5573f5;
its `variants.tsv` stand-in file is still on disk but is no longer guarded,
because the stage no longer has anything to stand in for. The loader tests
below still exercise a `variants.tsv`-shaped table as a "real" table, which is
the point of that test: the guard is specific to stand-ins, not to filenames.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from papipeline.errors import DataContractError
from papipeline.execution.contracts import STAGE_TABLES
from papipeline.standins import (
    STANDIN_BANNER,
    UNBUILT_WITH_STANDIN,
    StandInError,
    is_standin_table,
    load_stage_table,
)

REPO = Path(__file__).resolve().parents[2]
STANDIN_DIR = REPO / "test_data" / "standins" / "unbuilt_stages"
TABLES = tuple(STAGE_TABLES[s][0] for s in UNBUILT_WITH_STANDIN)


class TestTheFixturesDeclareThemselves:
    @pytest.mark.parametrize("filename", TABLES)
    def test_the_file_exists(self, filename):
        assert (STANDIN_DIR / filename).exists(), (
            f"{filename} is declared as an external input for an unbuilt stage; "
            "without it the TEST DAG cannot resolve"
        )

    @pytest.mark.parametrize("filename", TABLES)
    def test_the_header_says_stand_in(self, filename):
        lines = (STANDIN_DIR / filename).read_text(encoding="utf-8").splitlines()
        header = "\n".join(lines[:6])
        assert STANDIN_BANNER in header, (
            f"{filename} has no {STANDIN_BANNER!r} banner in its comment header. "
            "A reader who opens this file must not be able to mistake it for output."
        )

    @pytest.mark.parametrize("filename", TABLES)
    def test_the_header_says_delete_me(self, filename):
        """It should say what to do about it, not just what it is."""
        text = (STANDIN_DIR / filename).read_text(encoding="utf-8").lower()
        assert "ticket 14" in text, f"{filename} does not cite the ticket that removes it"

    @pytest.mark.parametrize("filename", TABLES)
    def test_it_satisfies_the_contract_minimum(self, filename):
        """One row, because `stage_spec` requires non_empty and min_rows(1).

        Zero rows would also satisfy Snakemake and break the observatory, which
        is how a 0-row stand-in would have passed the DAG and failed later.
        """
        rows = [
            line for line in (STANDIN_DIR / filename).read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.startswith("#")
        ]
        assert len(rows) == 2, (
            f"{filename} has {len(rows) - 1} data rows; the contract needs at "
            "least 1 and the file needs no more"
        )

    @pytest.mark.parametrize("filename", TABLES)
    def test_every_value_marks_itself(self, filename):
        for line in (STANDIN_DIR / filename).read_text(encoding="utf-8").splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            for cell in line.split("\t"):
                if cell.strip().lower() in {"standin", "standin_not_called", "standin_not_computed"}:
                    continue
                # headers and the synthetic sample id are not calls
                if line.split("\t")[0] == cell:
                    continue
                assert not cell.strip().isdigit(), (
                    f"{filename} contains a bare number in a result cell: {cell!r}. "
                    "A stand-in should not carry a value that reads as a measurement"
                )


class TestTheLoaderRefuses:
    @pytest.mark.parametrize("filename", TABLES)
    def test_it_detects_a_stand_in(self, filename):
        assert is_standin_table(STANDIN_DIR / filename) is True

    @pytest.mark.parametrize("filename", TABLES)
    @pytest.mark.parametrize("stage", UNBUILT_WITH_STANDIN)
    def test_loading_one_raises(self, filename, stage):
        with pytest.raises(StandInError) as excinfo:
            load_stage_table(STANDIN_DIR / filename, stage)
        message = str(excinfo.value)
        assert STANDIN_BANNER in message
        assert "ticket 14" in message

    def test_the_error_is_a_data_contract_error(self):
        """A caller catching PipelineError still stops."""
        assert issubclass(StandInError, DataContractError)
        for stage, filename in zip(UNBUILT_WITH_STANDIN, TABLES):
            with pytest.raises(DataContractError):
                load_stage_table(STANDIN_DIR / filename, stage)

    def test_a_real_table_loads(self, tmp_path):
        """The refusal must be specific, not a blanket failure.

        Named rather than indexed out of `UNBUILT_WITH_STANDIN`: that tuple is
        empty now that every stage is built, but `load_stage_table` is still the
        sanctioned reader and still has to accept a table that is not a
        stand-in. `variants` is named because it is the stage this class has
        always exercised the reader against, and its columns are the ones the
        assertion below writes.
        """
        stage = "variants"
        real = tmp_path / STAGE_TABLES[stage][0]
        real.write_text(
            "sample_id\tposition\treference\talternate\tvariant_type\tcall_status\n"
            "TEST_PA_001\t1234\tA\tG\tsnp\tconfirmed\n",
            encoding="utf-8",
        )
        assert is_standin_table(real) is False
        rows = load_stage_table(real, stage)
        assert len(rows) == 1
        assert rows[0]["call_status"] == "confirmed"

    def test_a_standin_with_edited_rows_is_still_refused(self, tmp_path):
        """The banner survives hand-editing, so the banner is checked too."""
        path = tmp_path / "variants.tsv"
        path.write_text(
            f"# {STANDIN_BANNER} - hand-edited\n"
            "sample_id\tposition\treference\talternate\tvariant_type\tcall_status\n"
            "TEST_PA_001\t1\tA\tG\tsnp\tconfirmed\n",
            encoding="utf-8",
        )
        assert is_standin_table(path) is True
        with pytest.raises(StandInError):
            load_stage_table(path, "variants")


class TestNoOtherModuleReadsTheseTables:
    """The structural half of the rule.

    A stand-in is only safe because nothing consumes it. When ticket 14 adds
    real consumers this test will fail, and that failure is the signal to
    replace `load_stage_table` with a real loader - not an obstacle to delete.
    """

    #: Modules allowed to name these files: the contract that declares them,
    #: this package, and the workflow that supplies them as external inputs.
    ALLOWED = {
        "papipeline.execution.contracts",
        "papipeline.standins",
    }

    def test_nothing_else_reads_them(self):
        offenders = []
        for path in sorted((REPO / "papipeline").rglob("*.py")):
            if "__pycache__" in path.parts:
                continue
            rel = path.relative_to(REPO).with_suffix("")
            module = ".".join(rel.parts)
            if module in self.ALLOWED:
                continue
            try:
                tree = ast.parse(path.read_text(encoding="utf-8"))
            except SyntaxError:
                continue
            for node in ast.walk(tree):
                literals = []
                if isinstance(node, ast.Constant) and isinstance(node.value, str):
                    literals = [node.value]
                for text in literals:
                    if any(text.endswith(f"/{t}") or text == t for t in TABLES):
                        offenders.append(f"{module}: {text!r}")
        assert not offenders, (
            "these modules read a table that is a stand-in fixture:\n  "
            + "\n  ".join(offenders)
            + "\n\nUse papipeline.standins.load_stage_table, which refuses a "
            "stand-in. When ticket 14 builds the real stages, replace that "
            "loader and revisit this test - see TASKS.md."
        )

    def test_the_workflow_declares_standins_as_inputs_only(self):
        """No rule wires a stand-in, because no stage is stood in for any more.

        `variants.tsv` is a legitimate output filename for the `variants` rule.
        What must never happen is a rule *producing* a stand-in, because then the
        stand-in would be indistinguishable from the stage's own output - and
        the stronger form of that, a rule *consuming* one for a stage that
        computes its own table, is what `dag-resolve` removed for
        `recombination`.
        """
        snakefile = (REPO / "workflow" / "Snakefile").read_text(encoding="utf-8")
        # Comments may still *mention* a retired name - that is the record of why
        # it went. What must not survive is a use, so comments are stripped
        # before anything is asserted about the names.
        code = "\n".join(
            line for line in snakefile.splitlines()
            if not line.lstrip().startswith("#")
        )

        # `STANDIN_DIR` and `RECOMBINATION_STANDIN` are gone. They existed so
        # TEST had *something* to read while stage 8 was unbuilt, and
        # `RECOMBINATION_STANDIN` sitting in `rule recombination`'s `test_only=`
        # is exactly what would have let a fabricated table keep satisfying the
        # edge from stage 9 after the stage began computing a real one. TEST now
        # reads committed gubbins-shaped output from
        # `test_data/intermediate/gubbins/` through the real parser.
        for retired in (
            "STANDIN_DIR",
            "RECOMBINATION_STANDIN",
            "VARIANTS_STANDIN",
            "COHORT_VARIANTS_STANDIN",
            "SIMILARITY_STANDIN",
        ):
            assert f"{retired} =" not in code, (
                f"{retired} is still declared as a module-level path. "
                "standins.UNBUILT_WITH_STANDIN is empty, so a stand-in path here "
                "is a fixture for a stage that computes its own output."
            )
            assert f"[{retired}]" not in code and f"{retired}," not in code, (
                f"{retired} is still wired into a rule's inputs"
            )

        # No rule may declare a stand-in as an output.
        import re

        for match in re.finditer(r"(?ms)^rule (\w+):\n.*?(?=^rule |\Z)", snakefile):
            block = match.group(0)
            outputs = re.search(r"(?m)^    output:\n((?:        .*\n)+)", block)
            if not outputs:
                continue
            assert "STANDIN" not in outputs.group(1), (
                f"rule {match.group(1)!r} declares a stand-in as an output; a "
                "stand-in must be an external input, never something a stage "
                "produces"
            )
