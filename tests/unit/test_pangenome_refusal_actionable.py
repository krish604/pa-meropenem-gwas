"""The panaroo refusals are actionable: command, flags, and both output paths.

**The defect.** The pipeline never runs panaroo - there is no `subprocess` call
in `papipeline/adapters/panaroo.py`, which is asserted from the AST below rather
than assumed. Its output is therefore *operator-provisioned*, and every refusal
that fires because it is absent has to say how to produce it.

None of them did. `build_pangenome_from_panaroo` said "panaroo did not write
gene_presence_absence.csv in <dir>" and stopped. That is a fact, not a remedy:
the reader is left to know that panaroo exists, that it is the tool that writes
that file, and which invocation produces it with the flags that matter. On
real Bakta data the obvious invocation *fails* - panaroo raises
`ValueError: Invalid gene sequence!` on any CDS whose length is not a multiple
of 3, is under 34 bp, or contains an internal stop, unless
`--remove-invalid-genes` is passed. That is 10 genes out of 63,359 on the
verified 10-isolate cohort, so an operator following the refusal without the flag
gets a crash and no indication why.

**What is asserted.** For each refusal that can fire on a missing or mismatched
panaroo output:

  * the command, including every flag in the invocation recorded in
    `docs/environment-arm64.md` section 3;
  * the *resolved* output directory the parser reads, not a placeholder;
  * **both** files the pipeline reads from it - `gene_presence_absence.csv` for
    stage 7 and `core_gene_alignment_filtered.aln` for stage 8 - because one
    command produces both and naming only one sends the operator round again;
  * that the pipeline does not invoke panaroo, so the premise of all of it holds.

**A rule was broken here and is now asserted against.** The first version of
stage 8's refusal said `panaroo build -i ... --remove-invalid-files`. That
subcommand word and that flag were written from memory and are wrong. AGENTS.md
rule 1: never invent tool flags. `test_the_documented_command_is_the_one_quoted`
pins the message to the documented invocation so it cannot drift back.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from papipeline.adapters import panaroo as panaroo_adapter
from papipeline.adapters.gubbins import CORE_ALIGNMENT_BASENAME
from papipeline.errors import DataContractError

REPO = Path(__file__).resolve().parents[2]
PANAROO_ADAPTER = REPO / "papipeline" / "adapters" / "panaroo.py"

#: The invocation `docs/environment-arm64.md` section 3 records as verified
#: end-to-end against real Bakta GFF3 (exit 0, 10,019 gene families).
VERIFIED_INVOCATION = (
    "panaroo -i *.gff3 -o out -t 4 --clean-mode moderate --remove-invalid-genes"
)


class TestThePipelineDoesNotRunPanaroo:
    """The premise of every message in this file."""

    def test_the_adapter_imports_no_subprocess(self):
        """Parsed, not grepped.

        A text search would match the word `subprocess` in a docstring - and the
        module docstring and these refusal messages both mention it, precisely
        because they explain that there is no such call. Only the AST can tell
        an import from a sentence.
        """
        tree = ast.parse(PANAROO_ADAPTER.read_text(encoding="utf-8"))
        imported = {
            alias.name.split(".")[0]
            for node in ast.walk(tree)
            if isinstance(node, (ast.Import, ast.ImportFrom))
            for alias in getattr(node, "names", [])
        }
        assert "subprocess" not in imported, (
            "adapters/panaroo.py now imports subprocess. panaroo's output is "
            "operator-provisioned and every refusal below tells the operator to "
            "run it; if the pipeline runs it, those messages are wrong."
        )

    def test_the_adapter_calls_nothing(self):
        """No `run`/`call`/`check_output`/`Popen` on any module object."""
        tree = ast.parse(PANAROO_ADAPTER.read_text(encoding="utf-8"))
        process_calls = {
            node.func.attr
            for node in ast.walk(tree)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
            and node.func.attr in {
                "run", "call", "check_call", "check_output", "Popen",
            }
        }
        assert not process_calls, (
            f"adapters/panaroo.py invokes {sorted(process_calls)}. Same "
            "consequence as the import: the operator is being told to do "
            "something the pipeline already does."
        )

    def test_the_stage_does_not_run_it_either(self):
        """`stages/pangenome.py` reads the output; it does not produce it."""
        tree = ast.parse(
            (REPO / "papipeline" / "stages" / "pangenome.py").read_text(
                encoding="utf-8"
            )
        )
        imported = {
            alias.name.split(".")[0]
            for node in ast.walk(tree)
            if isinstance(node, (ast.Import, ast.ImportFrom))
            for alias in getattr(node, "names", [])
        }
        assert "subprocess" not in imported, (
            "stages/pangenome.py imports subprocess. It must read "
            "operator-provisioned panaroo output, never produce it."
        )


class TestEveryPanarooRefusalIsActionable:
    """One place that checks the whole family, so a new refusal cannot be added
    without the same properties."""

    @pytest.fixture
    def panaroo_dir(self, tmp_path: Path) -> Path:
        directory = tmp_path / "panaroo"
        directory.mkdir()
        return directory

    def _assert_actionable(self, message: str, panaroo_dir: Path) -> None:
        assert "panaroo -i" in message, (
            f"the refusal must carry the invocation, not just the tool name:\n"
            f"{message}"
        )
        for flag in ("-i", "-o", "-t", "--clean-mode", "--remove-invalid-genes"):
            assert flag in message, (
                f"{flag} is in the verified invocation but not in the refusal. "
                "AGENTS.md rule 1."
            )
        assert "--remove-invalid-genes" in message, (
            "this flag is not cosmetic: panaroo RAISES on a CDS whose length is "
            "not a multiple of 3, is under 34 bp, or has an internal stop. On "
            "the verified 10-isolate cohort that is 10 genes out of 63,359, so "
            "an operator following a refusal without it gets a crash and no "
            "explanation."
        )
        assert str(panaroo_dir) in message, (
            "the output directory must be the RESOLVED path the parser reads. A "
            "placeholder leaves the operator to work out whether they are in the "
            "right run."
        )
        assert panaroo_adapter.PRESENCE_CSV in message, (
            "stage 7's own input must be named"
        )
        assert CORE_ALIGNMENT_BASENAME in message, (
            "stage 8's input must be named too. One command produces both, so "
            "naming only one sends the operator round again."
        )
        assert "does not run panaroo" in message, (
            "the refusal must say the tool is out of band, or it reads as a bug "
            "in this pipeline rather than a step the operator owns"
        )

    def test_the_missing_presence_table(self, panaroo_dir):
        with pytest.raises(DataContractError) as caught:
            panaroo_adapter.build_pangenome_from_panaroo(panaroo_dir, ["A", "B"])
        self._assert_actionable(str(caught.value), panaroo_dir)

    def test_the_cohort_mismatch(self, panaroo_dir):
        """A table from an earlier run is the likely cause, so say so."""
        (panaroo_dir / panaroo_adapter.PRESENCE_CSV).write_text(
            "Gene,Name,Annotation,OLD_ISOLATE\n"
            "group_1,group_1,,A\n"
            "group_2,group_2,,A\n",
            encoding="utf-8",
        )
        with pytest.raises(DataContractError) as caught:
            panaroo_adapter.build_pangenome_from_panaroo(panaroo_dir, ["A", "B"])
        message = str(caught.value)
        assert "different manifest" in message, (
            "the common cause of a mismatch is a table left by an earlier run; "
            "naming it saves the operator from assuming the parser is broken"
        )
        self._assert_actionable(message, panaroo_dir)

    def test_the_missing_toolchain(self):
        """The preflight, which runs before any path is resolved.

        So the directory is rendered as a placeholder rather than omitted: a
        refusal that names a missing tool and says nothing about where the output
        should go is only half a remedy.
        """
        message = panaroo_adapter.provision_instructions()
        assert "panaroo -i" in message
        assert "--remove-invalid-genes" in message
        assert "<this run's intermediate>/panaroo" in message, (
            "with no run in hand the directory is rendered the way the operator "
            "would recognise it, not dropped"
        )
        assert panaroo_adapter.PRESENCE_CSV in message
        assert CORE_ALIGNMENT_BASENAME in message

    def test_the_empty_pangenome_real_refusal_names_the_command(self, tmp_path):
        """`_refuse_an_empty_pangenome` in REAL is the third member of the family.

        A header-only presence table is panaroo failing, not a cohort with no
        genes, and the remedy is to re-run the command - which has to be quoted.
        """
        from papipeline.errors import StageError
        from papipeline.models import RunMode
        from papipeline.stages.pangenome import Pangenome, _refuse_an_empty_pangenome

        source_root = tmp_path / "intermediate"
        (source_root / "panaroo").mkdir(parents=True)
        with pytest.raises(StageError) as caught:
            _refuse_an_empty_pangenome(
                Pangenome(
                    sample_ids=("A", "B"),
                    genes=(),
                    presence={},
                    core=(),
                    accessory=(),
                ),
                RunMode.REAL,
                source_root,
            )
        message = str(caught.value)
        self._assert_actionable(
            message, source_root / "panaroo"
        )


class TestTheQuotedCommandIsTheDocumentedOne:
    """AGENTS.md rule 1, asserted rather than trusted."""

    def test_the_documented_invocation_still_exists(self):
        doc = (REPO / "docs" / "environment-arm64.md").read_text(encoding="utf-8")
        assert VERIFIED_INVOCATION in doc, (
            "docs/environment-arm64.md section 3 no longer records the verified "
            "invocation. PROVISION_COMMAND's provenance comment and this test "
            "have to be updated together, deliberately."
        )

    def test_every_flag_of_it_is_quoted(self):
        for flag in ("-i", "-o", "-t", "--clean-mode", "--remove-invalid-genes"):
            assert flag in panaroo_adapter.PROVISION_COMMAND, (
                f"{flag} is in the verified invocation but not in "
                "PROVISION_COMMAND"
            )

    def test_there_is_no_invented_subcommand(self):
        """The mistake this file documents, pinned so it cannot return."""
        assert "panaroo build" not in panaroo_adapter.PROVISION_COMMAND, (
            "`panaroo build` was written from memory and is not the invocation "
            "this project verified. panaroo's verified form here takes flags "
            "directly, with no subcommand word."
        )

    def test_the_output_directory_placeholder_is_the_intermediate_root(
        self,
    ):
        """`-o` points at the intermediate ROOT, not at the panaroo directory.

        `output_dir` is defined as `<intermediate_root>/panaroo`, and panaroo
        writes into `<-o>/panaroo`. Handing it the subdirectory leaves the
        pipeline looking one level away from where the files landed - a mistake
        the placeholder makes visible rather than hiding.
        """
        assert "-o <this run's intermediate>" in panaroo_adapter.PROVISION_COMMAND
        assert "-o <this run's intermediate>/panaroo" not in (
            panaroo_adapter.PROVISION_COMMAND
        ), (
            "this would put panaroo's output one directory too deep. "
            "adapters.panaroo.output_dir(intermediate_root) is "
            "<intermediate_root>/panaroo, and panaroo appends its own "
            "'panaroo' to -o."
        )