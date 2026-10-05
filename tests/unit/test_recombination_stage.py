"""Stage 8 (recombination): the contract, the three refusals, and the real run.

Five groups, in descending order of how much they would cost to be wrong.

1. **The three refusals**, each asserting its own *named* reason. One message
   covering two failures cannot tell an operator which one they have, so these
   assert on the distinguishing fragment, not merely on "it raised".
2. **The block-table contract**, including the `mean_branch_length` derivation,
   which is this stage's own arithmetic and the part most likely to be wrong.
3. **TEST-mode output shape**, and the `gubbins_dir` requirement.
4. **The binary-resolution machinery**, tested against a stand-in script rather
   than against gubbins, so the PATH-isolation logic is exercised on a machine
   where gubbins is not installed.
5. **One REAL-mode test**, `skipif`-gated on the source binary being present at
   its absolute path, which skips with a reason rather than failing.

**The internal-node label is the thing most likely to be silently wrong.**
gubbins writes ``(a:1,b:1)Node_1:5`` — the internal node's label sits *after* the
closing parenthesis, which is not standard Newick. A generic parser either drops
the label or attributes it to the preceding child, and both mistakes pair a
branch length with the wrong SNP count, which is the one number a reader
checks. That is asserted against a hand-checked tree below.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

from papipeline.adapters import gubbins as gubbins_adapter
from papipeline.errors import (
    DataContractError,
    ModeNotAllowedError,
    ToolExecutionError,
    ToolNotAvailableError,
)
from papipeline.models import RunMode
from papipeline.stages import recombination as stage

REPO = Path(__file__).resolve().parents[2]
SNAKEFILE = REPO / "workflow" / "Snakefile"

#: The source-built gubbins. `laptop.yaml` declares `tools/gubbins/bin` in
#: `tool_search_dirs` and `scripts/build_gubbins_from_source.sh` writes it there,
#: but that build is a gitignored artefact — so its absence must SKIP the real
#: test, not fail it. Located by absolute path rather than through `which`,
#: because `which` would also find the broken conda build.
SOURCE_GUBBINS = REPO / "tools" / "gubbins" / "bin" / "gubbins"

SKIP_NO_SOURCE_BUILD = pytest.mark.skipif(
    not SOURCE_GUBBINS.is_file(),
    reason=(
        f"the source-built gubbins is absent at {SOURCE_GUBBINS}. It is a "
        "gitignored build artefact; produce it with "
        "scripts/build_gubbins_from_source.sh. The conda osx-arm64 build cannot "
        "stand in for it — it segfaults with exit 139 (see "
        "tests/integration/test_gubbins_source_build.py)."
    ),
)


# --------------------------------------------------------------------------
# Fixtures. Deliberately gubbins-SHAPED, and checked against the real formats.
# --------------------------------------------------------------------------

PER_BRANCH_HEADER = "\t".join(
    [
        "Node", "Total SNPs", "Number of SNPs Inside Recombinations",
        "Number of SNPs Outside Recombinations", "Number of Recombination Blocks",
        "Bases in Recombinations Including Gaps",
        "Cumulative Bases in Recombinations Including Gaps",
        "Bases in Recombinations Excluding Gaps",
        "Cumulative Bases in Recombinations Excluding Gaps",
        "r/m", "rho/theta", "Genome Length", "Bases in Clonal Frame",
    ]
)


def write_statistics(directory: Path, rows) -> Path:
    """A `per_branch_statistics.csv` with the real 13-column header.

    The header is copied from a real gubbins run rather than invented, because
    the point of these tests is that the parser reads the file gubbins writes.
    A shortened header would let a position-based parser pass here and fail on
    real output.
    """
    path = directory / f"{gubbins_adapter.PREFIX}{gubbins_adapter.PER_BRANCH_STATISTICS_SUFFIX}"
    lines = [PER_BRANCH_HEADER]
    for node, snps, blocks in rows:
        lines.append("\t".join([node, str(snps), "0", str(snps), str(blocks),
                                "0", "0", "0", "0", "0.0", "0.0", "100", "100"]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_tree(directory: Path, newick: str) -> Path:
    path = directory / f"{gubbins_adapter.PREFIX}{gubbins_adapter.NODE_LABELLED_TREE_SUFFIX}"
    path.write_text(newick, encoding="utf-8")
    return path


@pytest.fixture
def gubbins_output(tmp_path):
    """A minimal but structurally faithful gubbins output directory.

    The tree is `(a:10,b:20)Node_1:100;` — one internal node whose label is after
    the parenthesis, which is the case a naive parser gets wrong.
    """
    directory = tmp_path / "gubbins"
    directory.mkdir()
    write_statistics(directory, [("a", 5, 1), ("b", 0, 0), ("Node_1", 3, 0)])
    write_tree(directory, "(a:10,b:20)Node_1:100;")
    (directory / f"{gubbins_adapter.PREFIX}{gubbins_adapter.SNP_ALIGNMENT_SUFFIX}").write_text(
        ">a\nACGT\n>b\nACGA\n", encoding="utf-8"
    )
    return directory


# --------------------------------------------------------------------------
# 1. The three refusals, each with its own named reason.
# --------------------------------------------------------------------------


class TestTheThreeRefusals:
    def test_a_shut_real_gate_names_the_flag_the_override_and_the_overlay(self):
        """Refusal (a): `runtime.allow_real_mode` is false.

        Asserted on the three things an operator needs in order to act: the
        config flag to change, the one-session environment override, and the
        overlay that currently holds it shut.
        """
        class _Config:
            machine_name = "laptop"
            runtime = {"allow_real_mode": False}

        with pytest.raises(ModeNotAllowedError) as excinfo:
            gubbins_adapter.require_real_gate(_Config())
        message = str(excinfo.value)
        assert "allow_real_mode" in message
        assert "PIPELINE_ALLOW_REAL_MODE=1" in message
        assert "laptop" in message

    def test_an_open_gate_lets_the_run_proceed(self):
        """A gate on the run, not a refusal of the stage.

        The taxonomy derives "refuses REAL" from the AST and treats a gate as the
        opposite; a gate that never opens would be a refusal wearing a gate's
        name.
        """
        class _Config:
            machine_name = "bigmachine"
            runtime = {"allow_real_mode": True}

        assert gubbins_adapter.require_real_gate(_Config()) is None

    def test_the_gate_is_checked_before_the_alignment(self):
        """Order is load-bearing.

        An operator whose gate is shut is not the operator whose alignment is
        missing, and telling them to supply an alignment would not unblock them.
        """
        class _Config:
            machine_name = "laptop"
            runtime = {"allow_real_mode": False}

        with pytest.raises(ModeNotAllowedError) as excinfo:
            gubbins_adapter.run_gubbins(
                Path("/nonexistent/core_gene_alignment_filtered.aln"),
                Path("/tmp/never-created"),
                runner=Path("/r"), binary=Path("/b"), threads=1, config=_Config(),
            )
        assert "allow_real_mode" in str(excinfo.value)
        assert "core_gene_alignment_filtered" not in str(excinfo.value)

    def test_a_missing_alignment_names_the_expected_path(self):
        """Refusal (b): the input is absent.

        Names the path AND says what the absence means, because "no
        recombination" and "did not look" are different findings and the message
        is the only place that distinction survives.
        """
        missing = Path("/nonexistent/intermediate/panaroo/core_gene_alignment_filtered.aln")
        with pytest.raises(DataContractError) as excinfo:
            gubbins_adapter.run_gubbins(
                missing,
                Path("/tmp/never-created"),
                runner=Path("/nonexistent/run_gubbins.py"),
                binary=Path("/nonexistent/gubbins"),
                threads=1,
            )
        message = str(excinfo.value)
        assert str(missing) in message
        assert "not an absence of" in message
        assert "panaroo" in message

    def test_a_segfaulting_binary_is_refused_by_name_as_exit_139(self):
        """Refusal (c): the binary dies with SIGSEGV.

        Asserted on `139`, the conda osx-arm64 build and the unresolved zlib
        symbols, because those are what make the diagnosis. "gubbins failed" is
        not actionable; "your binary is the mis-packaged one, here is the
        upstream issue" is.
        """
        with pytest.raises(ToolNotAvailableError) as excinfo:
            gubbins_adapter.refuse_segfaulting_binary(
                Path("/opt/conda/envs/pa-gubbins/bin/gubbins"),
                gubbins_adapter.CONDA_BUILD_VERSION,
                gubbins_adapter.SEGFAULT_EXIT,
                Path("/tmp/self_test"),
            )
        message = str(excinfo.value)
        assert "139" in message
        assert "SIGSEGV" in message
        assert "osx-arm64" in message
        assert "_gzopen" in message
        assert "452" in message
        assert gubbins_adapter.SOURCE_BUILD_VERSION in message

    def test_the_three_refusals_are_three_different_messages(self):
        """No message may serve two failures.

        If two causes shared a message, an operator holding it could not tell
        which one they had, and the whole point of naming the cause is lost.
        """
        messages = set()

        with pytest.raises(DataContractError) as missing:
            gubbins_adapter.run_gubbins(
                Path("/nonexistent/aln"), Path("/tmp/x"),
                runner=Path("/r"), binary=Path("/b"), threads=1,
            )
        messages.add(str(missing.value))

        with pytest.raises(ToolNotAvailableError) as segv:
            gubbins_adapter.refuse_segfaulting_binary(
                Path("/b"), "3.4.3", 139, Path("/tmp/x")
            )
        messages.add(str(segv.value))

        class _Config:
            machine_name = "laptop"
            runtime = {"allow_real_mode": False}

        with pytest.raises(ModeNotAllowedError) as gate:
            gubbins_adapter.require_real_gate(_Config())
        messages.add(str(gate.value))

        assert len(messages) == 3, "a refusal message is being reused for two causes"
        # And each is distinguishable on its own terms.
        assert "allow_real_mode" in str(gate.value)
        assert "no core gene alignment" in str(missing.value)
        assert "SIGSEGV" in str(segv.value)


# --------------------------------------------------------------------------
# 2. The block-table contract.
# --------------------------------------------------------------------------


class TestTheBlockTable:
    def test_the_columns_are_the_declared_contract(self):
        from papipeline.execution.contracts import STAGE_TABLES

        assert stage.BLOCK_COLUMNS == STAGE_TABLES["recombination"][1], (
            "the writer and the contract have drifted; one of them validates "
            "nothing while appearing to"
        )

    def test_it_is_not_a_per_sample_table(self):
        """`node` holds `Node_<n>` too, so the dense-per-sample property is false.

        Asserted because the table looks like a sample table — it has a leading
        identifier column and plausible counts — and adding it to
        `DENSE_PER_SAMPLE_STAGES` would assert something false about it.
        """
        from papipeline.execution.contracts import DENSE_PER_SAMPLE_STAGES

        assert "recombination" not in DENSE_PER_SAMPLE_STAGES

    def test_the_internal_node_label_is_read_from_after_the_parenthesis(self):
        """`(a:10,b:20)Node_1:100;` must yield a node labelled `Node_1`.

        gubbins is not writing standard Newick: the internal node's label sits
        *after* the closing parenthesis. A parser that reads the label before the
        children attributes it to `b` — so `b` would report a clade mean of
        (20+100)/2 = 60 and `Node_1` would not exist at all, which would fail the
        join. The assertions below pin both halves of that: the label exists, and
        no tip has absorbed it.
        """
        lengths = stage.mean_branch_lengths("(a:10,b:20)Node_1:100;")
        assert "Node_1" in lengths, "the post-parenthesis label was not read"
        assert lengths["b"] == pytest.approx(20.0), (
            "the internal label was attributed to the preceding child"
        )
        # Node_1's clade is {own edge 100, a:10, b:20}.
        assert lengths["Node_1"] == pytest.approx(130.0 / 3.0)

    def test_mean_branch_length_is_the_mean_over_the_whole_clade(self):
        """The root's clade mean is the mean of every edge, including its own.

        Hand-checked: edges are 10, 20 and 100, so the mean is 130/3.
        `Node_1`'s own clade is {its edge 100, a:10, b:20} = 130/3, and the root
        is `Node_1` here so both agree.
        """
        lengths = stage.mean_branch_lengths("(a:10,b:20)Node_1:100;")
        assert lengths["Node_1"] == pytest.approx(130.0 / 3.0)

    def test_a_tip_gets_its_own_edge_length_not_an_empty_cell(self):
        """The reason the clade definition was chosen.

        A tip has no children, so "mean of the immediate children" would leave
        every tip row's `mean_branch_length` empty — and an empty cell in a
        branch-length column is a hole a later reader fills with a zero, which
        then reads as "no evolution on this branch".
        """
        lengths = stage.mean_branch_lengths("(a:10,b:20)Node_1:100;")
        assert lengths["a"] == pytest.approx(10.0)
        assert lengths["b"] == pytest.approx(20.0)

    def test_a_zero_length_root_edge_is_kept_as_a_real_zero(self):
        """gubbins midpoint-roots, so the root edge can legitimately be 0.0.

        Filling it with a substitute would misreport the rooting; the zero is a
        fact about the tree.
        """
        lengths = stage.mean_branch_lengths("(a:10,b:20)Node_1:0.0;")
        assert lengths["Node_1"] == pytest.approx(30.0 / 3.0)

    def test_a_row_carries_all_four_declared_columns(self, gubbins_output):
        rows = stage.run(
            _config(), _manifest(), RunMode.TEST, gubbins_dir=gubbins_output
        )
        assert len(rows) == 3
        for row in rows:
            assert set(row) == set(stage.BLOCK_COLUMNS)

    def test_recombination_detected_is_a_zero_or_one_from_the_block_count(
        self, tmp_path
    ):
        """Read from `Number of Recombination Blocks`, not from the SNP count.

        A branch can carry many SNPs and no recombination (a long branch from
        substitution alone), so deriving the flag from `n_snps` would report
        recombination where there was none — the most consequential way this
        table could be wrong.
        """
        directory = tmp_path / "gb"
        directory.mkdir()
        write_statistics(
            directory,
            [
                ("many_snps_no_recombination", 900, 0),
                ("few_snps_recombination", 2, 1),
                ("Node_1", 902, 1),
            ],
        )
        write_tree(directory, "(many_snps_no_recombination:1,few_snps_recombination:1)Node_1:1;")
        rows = {r["node"]: r for r in stage.run(
            _config(), _manifest(), RunMode.TEST, gubbins_dir=directory
        )}
        assert rows["many_snps_no_recombination"]["recombination_detected"] == 0
        assert rows["few_snps_recombination"]["recombination_detected"] == 1

    def test_a_node_in_the_statistics_but_not_the_tree_is_refused(self, tmp_path):
        """Refuses rather than writing a partial join.

        A node in one file and not the other means the two describe different
        trees. Filling the gap with a blank branch length would produce a table
        that looks complete and pairs every later row with the wrong number.
        """
        directory = tmp_path / "gb"
        directory.mkdir()
        write_statistics(directory, [("a", 1, 0), ("b", 1, 0)])
        write_tree(directory, "(a:1)Node_9:1;")
        with pytest.raises(DataContractError) as excinfo:
            stage.run(_config(), _manifest(), RunMode.TEST, gubbins_dir=directory)
        message = str(excinfo.value)
        assert "b" in message
        assert "different trees" in message

    def test_a_labelled_node_the_statistics_omit_is_refused(self, tmp_path):
        """The inverse. A dropped node's SNP count would read as 'no recombination'."""
        directory = tmp_path / "gb"
        directory.mkdir()
        write_statistics(directory, [("a", 1, 0), ("b", 1, 0)])
        write_tree(directory, "(a:1,b:1)Node_1:1;")
        with pytest.raises(DataContractError) as excinfo:
            stage.run(_config(), _manifest(), RunMode.TEST, gubbins_dir=directory)
        assert "Node_1" in str(excinfo.value)

    def test_node_matching_is_by_equality_so_node_1_cannot_match_node_11(self):
        """`Node_1` must not match `Node_11` under a prefix or regex match."""
        # Both labels present in one tree; each keeps its OWN edge length.
        # Under a prefix match `Node_1` would absorb `Node_11`'s numbers, and
        # since the join in `block_rows` is a dict lookup on this map, that would
        # pair the wrong branch length with the wrong SNP count.
        lengths = stage.mean_branch_lengths("(Node_1:3,(Node_11:7)Node_12:11)Node_13:1;")
        assert lengths["Node_1"] == pytest.approx(3.0)
        assert lengths["Node_11"] == pytest.approx(7.0)
        assert set(lengths) == {"Node_1", "Node_11", "Node_12", "Node_13"}

        # And the join keeps them apart.
        statistics = [
            {"node": "Node_1", "n_snps": 1, "recombination_detected": 0},
            {"node": "Node_11", "n_snps": 2, "recombination_detected": 1},
            {"node": "Node_12", "n_snps": 3, "recombination_detected": 0},
            {"node": "Node_13", "n_snps": 4, "recombination_detected": 0},
        ]
        rows = {r["node"]: r for r in stage.block_rows(statistics, lengths)}
        assert rows["Node_1"]["n_snps"] == 1
        assert rows["Node_11"]["n_snps"] == 2
        assert rows["Node_11"]["mean_branch_length"] == pytest.approx(7.0)

    def test_the_table_is_written_with_the_declared_header(self, gubbins_output, tmp_path):
        out = tmp_path / "recombination.tsv"
        stage.run(_config(), _manifest(), RunMode.TEST, gubbins_dir=gubbins_output, out_path=out)
        lines = out.read_text(encoding="utf-8").splitlines()
        assert lines[0].split("\t") == list(stage.BLOCK_COLUMNS)


# --------------------------------------------------------------------------
# 3. TEST-mode shape and its required input.
# --------------------------------------------------------------------------


class TestTestModeShape:
    def test_test_mode_needs_an_explicit_gubbins_dir(self):
        """No default.

        Resolving one from configuration would let a TEST run report a table it
        measured from an input nobody chose — and the result would carry this
        stage's contract as though it had been selected.
        """
        with pytest.raises(DataContractError) as excinfo:
            stage.run(_config(), _manifest(), RunMode.TEST)
        message = str(excinfo.value)
        assert "gubbins_dir" in message
        assert "no default" in message

    def test_stub_mode_returns_the_header_and_no_rows(self):
        """Not shared with TEST.

        A stub table carrying one plausible fabricated row is the exact failure
        stub mode exists to prevent, so this asserts the row count is zero
        rather than merely that the call succeeded.
        """
        assert stage.run(_config(), _manifest(), RunMode.STUB) == []

    def test_a_missing_statistics_file_names_the_file(self, tmp_path):
        directory = tmp_path / "gb"
        directory.mkdir()
        write_tree(directory, "(a:1)Node_1:1;")
        with pytest.raises(DataContractError) as excinfo:
            stage.run(_config(), _manifest(), RunMode.TEST, gubbins_dir=directory)
        assert "per_branch_statistics.csv" in str(excinfo.value)
        assert "not that no recombination was found" in str(excinfo.value)

    def test_a_missing_tree_names_the_file_and_says_the_column_cannot_be_filled(
        self, tmp_path
    ):
        """The tree is the ONLY source of a branch length.

        `per_branch_statistics.csv` has no branch-length column, so without the
        tree the declared table cannot be filled at all. Saying so prevents a
        reader treating the omission as an oversight.
        """
        directory = tmp_path / "gb"
        directory.mkdir()
        write_statistics(directory, [("a", 1, 0)])
        with pytest.raises(DataContractError) as excinfo:
            stage.run(_config(), _manifest(), RunMode.TEST, gubbins_dir=directory)
        message = str(excinfo.value)
        assert "node_labelled" in message
        assert "no branch-length column" in message

    def test_the_masked_alignment_is_written_where_the_consumers_read_it(
        self, gubbins_output, tmp_path
    ):
        """`phylogeny_dir / core_snp_alignment.fasta`.

        That exact path is what `stages/phylogeny.py:368` reads and hands to
        IQ-TREE, and what `stages/gwas.py:502` records as its SNP source. An
        alignment written anywhere else is one no consumer reads.
        """
        phylogeny_dir = tmp_path / "phylo"
        stage.run(
            _config(), _manifest(), RunMode.TEST,
            gubbins_dir=gubbins_output, phylogeny_dir=phylogeny_dir,
        )
        written = phylogeny_dir / stage.SNP_ALIGNMENT_NAME
        assert written.is_file()
        assert written.read_text(encoding="utf-8").startswith(">a")

    def test_the_tree_is_kept_but_is_not_renamed_to_tree_nwk(self, gubbins_output, tmp_path):
        """Stage 9 owns `tree.nwk`.

        Two files claiming that path would leave a downstream reader unable to
        tell which inference its patristic distances were measured on.
        """
        source = gubbins_output / f"{gubbins_adapter.PREFIX}{gubbins_adapter.FINAL_TREE_SUFFIX}"
        source.write_text("(a:1,b:1)Node_1:1;", encoding="utf-8")
        phylogeny_dir = tmp_path / "phylo"
        written = stage.build_outputs(
            gubbins_dir=gubbins_output, phylogeny_dir=phylogeny_dir
        )
        assert "tree.nwk" not in written
        assert not (phylogeny_dir / "tree.nwk").exists()
        assert written["final_tree"].is_file()

    def test_a_missing_masked_alignment_is_refused_rather_than_reported_as_empty(
        self, tmp_path
    ):
        """Exit-0-with-no-alignment is a failure, not "no recombination"."""
        directory = tmp_path / "gb"
        directory.mkdir()
        write_statistics(directory, [("a", 1, 0)])
        write_tree(directory, "(a:1)Node_1:1;")
        with pytest.raises(DataContractError) as excinfo:
            stage.build_outputs(gubbins_dir=directory, phylogeny_dir=tmp_path / "p")
        assert "filtered_polymorphic_sites.fasta" in str(excinfo.value)


# --------------------------------------------------------------------------
# 4. Binary resolution and PATH isolation, without needing gubbins.
# --------------------------------------------------------------------------


def _fake_tool(directory: Path, name: str, *, body: str = "") -> Path:
    """A stand-in executable that behaves like the real tool.

    A shell script rather than a real binary so the PATH logic can be tested on
    a machine where gubbins is not installed — and so the *resolution* can be
    observed independently of whether the tool works, which is the whole reason
    this stage separates the two.
    """
    path = directory / name
    path.write_text(
        textwrap.dedent(
            f"""\
            #!/bin/sh
            echo "TOOL=$0"
            echo "PATH=$PATH"
            {body}
            """
        ),
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


class TestBinaryResolutionAndPathIsolation:
    def test_search_dirs_are_scanned_in_order_and_the_first_hit_wins(self, tmp_path):
        first = tmp_path / "a"
        second = tmp_path / "b"
        first.mkdir()
        second.mkdir()
        expected = _fake_tool(first, "gubbins")
        _fake_tool(second, "gubbins")
        assert gubbins_adapter.resolve_binary([str(first), str(second)]) == expected

    def test_a_directory_without_the_binary_is_skipped_not_fatal(self, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        real = tmp_path / "real"
        real.mkdir()
        _fake_tool(real, "gubbins")
        assert gubbins_adapter.resolve_binary([str(empty), str(real)]) == real / "gubbins"

    def test_resolution_reports_where_it_looked_when_nothing_is_found(self, tmp_path, monkeypatch):
        monkeypatch.setattr(gubbins_adapter.shutil, "which", lambda _: None)
        with pytest.raises(ToolNotAvailableError) as excinfo:
            gubbins_adapter.resolve_binary([str(tmp_path / "nowhere")])
        message = str(excinfo.value)
        assert "nowhere" in message
        assert "build_gubbins_from_source.sh" in message

    def test_the_isolated_directory_prepends_and_shadows_a_later_gubbins(self, tmp_path):
        """The core mechanism.

        `run_gubbins.py` has no flag naming a binary; it resolves `gubbins` by a
        left-to-right PATH scan and its only fallback *appends*, which can never
        beat an earlier entry. So the resolved binary has to be prepended.
        """
        source_dir = tmp_path / "src"
        conda_dir = tmp_path / "conda"
        source_dir.mkdir()
        conda_dir.mkdir()
        source = _fake_tool(source_dir, "gubbins")
        conda = _fake_tool(conda_dir, "gubbins")

        environment, isolated = gubbins_adapter.isolated_env(
            source, [str(conda_dir)], bin_dir=tmp_path / "iso"
        )
        entries = environment["PATH"].split(os.pathsep)
        assert entries[0] == str(isolated)
        assert entries.index(str(conda_dir)) > 0

        # And the shadowing resolves the way the tool's own scan would resolve it:
        # the FIRST directory containing an executable `gubbins` wins.
        winning_dir = next(
            entry for entry in entries if (Path(entry) / "gubbins").is_file()
        )
        assert winning_dir == str(isolated)
        assert (Path(winning_dir) / "gubbins").resolve() == source.resolve()
        assert (Path(winning_dir) / "gubbins").resolve() != conda.resolve()

    def test_the_isolated_directory_holds_exactly_one_entry(self, tmp_path):
        """Two entries would defeat the point.

        The purpose is to control one executable. A directory of convenience
        symlinks is a second place for a tool to come from unnoticed.
        """
        source_dir = tmp_path / "src"
        source_dir.mkdir()
        source = _fake_tool(source_dir, "gubbins")
        _fake_tool(source_dir, "raxml-ng")
        _, isolated = gubbins_adapter.isolated_env(source)
        assert [p.name for p in isolated.iterdir()] == ["gubbins"]

    def test_sysctl_is_retained_because_gubbins_needs_it(self, tmp_path):
        """Without `/usr/sbin`, gubbins' RAxML-NG builder raises UnboundLocalError.

        `utils.choose_executable_based_on_processor` probes CPU features by
        shelling out to `sysctl`, assigns a local only inside that branch, and
        reads it unconditionally on the next line. The workaround is therefore
        "keep a normal system PATH", not "make PATH more minimal" — so this
        asserts the tail is present rather than trimming it.
        """
        source_dir = tmp_path / "src"
        source_dir.mkdir()
        source = _fake_tool(source_dir, "gubbins")
        environment, _ = gubbins_adapter.isolated_env(source)
        entries = environment["PATH"].split(os.pathsep)
        assert "/usr/sbin" in entries
        assert "/usr/bin" in entries
        assert entries.index("/usr/sbin") > 0

    def test_the_probe_alignment_is_already_aligned_and_gap_free(self, tmp_path):
        """A gubbins input must be an alignment.

        Gaps would exercise the filtering path rather than the one under test,
        and a crash there would be reported as the packaging defect it is not.
        """
        path = gubbins_adapter.write_probe_alignment(tmp_path / "probe.aln")
        records, name, length = {}, None, set()
        for line in path.read_text(encoding="utf-8").splitlines():
            if line.startswith(">"):
                name = line[1:]
                records[name] = ""
            elif line.strip():
                records[name] += line.strip()
                length.add(len(line.strip()))
        assert len(records) == gubbins_adapter.PROBE_TAXA
        assert length == {60}, "sequence lines must be uniform so the file is aligned"
        assert gubbins_adapter.PROBE_SITES % 60 == 0, (
            "PROBE_SITES must divide by the line width or the probe is ragged"
        )
        assert set("".join(records.values())) <= set("ACGT")
        assert len(set(map(len, records.values()))) == 1

    def test_the_probe_alignment_is_deterministic(self, tmp_path):
        """A self-test that differs run to run cannot be compared to a baseline."""
        first = gubbins_adapter.write_probe_alignment(tmp_path / "one.aln").read_text()
        second = gubbins_adapter.write_probe_alignment(tmp_path / "two.aln").read_text()
        assert first == second

    def test_the_tree_builder_flag_is_the_spelling_run_gubbins_accepts(self):
        """`raxmlng`, no hyphen — read out of `run_gubbins.py --help`, not assumed.

        The hyphenated form is not accepted, so a plausible-looking value would
        fail only at the point where the stage is trying to do real work.
        """
        assert gubbins_adapter.TREE_BUILDER == "raxmlng"
        assert "-" not in gubbins_adapter.TREE_BUILDER

    def test_the_command_passes_the_alignment_by_absolute_path(self, tmp_path):
        """gubbins writes relative to its CWD, and is given a different CWD."""
        relative = Path("relative/aln.aln")
        command = gubbins_adapter.build_command(
            Path("/opt/env/bin/run_gubbins.py"), relative, tmp_path, 4
        )
        resolved = str(relative.resolve())
        assert command[-1] == resolved, (
            "gubbins writes relative to its CWD and is handed a different CWD, "
            "so a relative alignment path would be resolved against the wrong "
            "directory"
        )
        assert "/" in resolved and command[-1] != "relative/aln.aln"

    def test_the_command_carries_threads_and_prefix_from_arguments(self, tmp_path):
        """Hard rule 2: nothing environmental is hard-coded."""
        command = gubbins_adapter.build_command(
            Path("/r"), Path("/a.aln"), tmp_path, 7
        )
        assert "--threads" in command and "7" in command
        assert "--prefix" in command
        assert gubbins_adapter.PREFIX in command

    def test_a_non_139_failure_is_not_reported_as_the_segfault(self, tmp_path):
        """A usage error, a missing tree builder and a broken dylib differ.

        Telling an operator their tool is mis-packaged when it is merely
        mis-invoked sends them to rebuild a binary that is fine.
        """
        source_dir = tmp_path / "src"
        source_dir.mkdir()
        source = _fake_tool(source_dir, "gubbins", body="exit 3")
        exit_code, version = gubbins_adapter.self_test(source, tmp_path / "probe")
        assert exit_code == 3
        with pytest.raises(ToolExecutionError) as excinfo:
            gubbins_adapter.run_gubbins(
                _write_minimal_alignment(tmp_path),
                tmp_path / "work",
                runner=tmp_path / "runner",
                binary=source,
                threads=1,
                selftest=lambda: (exit_code, version),
            )
        message = str(excinfo.value)
        assert "not the exit-139 segfault" in message
        # It names the segfault only to rule it out, and does NOT assert the dylib
        # defect — that diagnosis belongs to the 139 path alone.
        assert "_gzopen" not in message
        assert "Run the binary by hand" in message

    def test_provenance_records_the_binary_version_and_not_the_banner(self, tmp_path):
        """The banner says 3.4.3 on BOTH builds; the binary says 3.4.2 on one.

        Recording the banner would report the broken conda build as the good
        source build, which is the exact inversion this module exists to avoid.
        """
        result = gubbins_adapter.GubbinsResult(
            work_dir=tmp_path, prefix="recombination",
            snp_alignment=tmp_path / "snp", final_tree=tmp_path / "tre",
            node_labelled_tree=tmp_path / "nl", per_branch_statistics=tmp_path / "pbs",
            log=tmp_path / "log", binary=tmp_path / "gubbins",
            binary_version=gubbins_adapter.SOURCE_BUILD_VERSION,
            isolated_bin_dir=tmp_path / "bin",
            argv=["run_gubbins.py", "--tree-builder", "raxmlng"],
            exit_code=0, seconds=31.64,
        )
        record = result.provenance()
        assert record["binary_version"] == gubbins_adapter.SOURCE_BUILD_VERSION
        assert record["exit_code"] == 0
        assert record["seconds"] == 31.64
        assert "raxmlng" in record["argv"]
        written = gubbins_adapter.write_provenance(result, tmp_path / "prov.json")
        assert json.loads(written.read_text())["binary_version"] == "3.4.2"

    def test_the_statistics_parser_refuses_a_shortened_header(self, tmp_path):
        """A missing column is refused rather than read positionally.

        `Total SNPs` and `Number of SNPs Inside Recombinations` have similar
        names and different meanings; a positional read swaps them and produces
        a number that looks right.
        """
        path = tmp_path / "pbs.csv"
        path.write_text("Node\tTotal SNPs\nA\t5\n", encoding="utf-8")
        with pytest.raises(DataContractError) as excinfo:
            gubbins_adapter.parse_per_branch_statistics(path)
        assert "Number of Recombination Blocks" in str(excinfo.value)

    def test_the_statistics_parser_refuses_a_truncated_row(self, tmp_path):
        """A short row means truncated, not "no SNPs"."""
        directory = tmp_path / "gb"
        directory.mkdir()
        path = directory / "pbs.csv"
        header = PER_BRANCH_HEADER.split("\t")
        path.write_text("\t".join(header) + "\nA\t5\n", encoding="utf-8")
        with pytest.raises(DataContractError) as excinfo:
            gubbins_adapter.parse_per_branch_statistics(path)
        assert "truncated" in str(excinfo.value)

    def test_the_statistics_parser_reads_columns_by_name_not_position(self, tmp_path):
        """Reordered columns must still join the right numbers."""
        path = tmp_path / "pbs.csv"
        reordered = [
            "Total SNPs", "Node", "Number of Recombination Blocks",
            "Genome Length", "Bases in Clonal Frame",
        ]
        # Row fields follow the reordered header: Total SNPs, Node, Blocks, ...
        path.write_text("\t".join(reordered) + "\n5\tA\t3\t100\t90\n", encoding="utf-8")
        rows = gubbins_adapter.parse_per_branch_statistics(path)
        assert rows == [{"node": "A", "n_snps": 5, "recombination_detected": 1}]


# --------------------------------------------------------------------------
# 5. One REAL-mode test, skipped cleanly when the source binary is absent.
# --------------------------------------------------------------------------


@SKIP_NO_SOURCE_BUILD
class TestAgainstTheSourceBuiltBinary:
    def test_the_source_build_is_the_one_that_does_not_segfault(self, tmp_path):
        """Run it for real, on a generated probe, and require a clean exit.

        The version string discriminates 3.4.2 (source) from 3.4.3 (conda), and
        the exit code confirms it independently: the conda build exits 139 here
        and the source build exits 0.
        """
        binary = gubbins_adapter.resolve_binary(
            [str(SOURCE_GUBBINS.parent)],
        )
        version = gubbins_adapter.binary_version(binary)
        exit_code, probe_version = gubbins_adapter.self_test(binary, tmp_path)
        assert probe_version == version
        assert version == gubbins_adapter.SOURCE_BUILD_VERSION, (
            f"{binary} reports {version}; the source build is expected to report "
            f"{gubbins_adapter.SOURCE_BUILD_VERSION}. A conda binary here would "
            "segfault (exit 139) rather than answer, so a 3.4.3 here means the "
            "wrong build was resolved."
        )
        assert exit_code == 0, f"the source build exited {exit_code}"

    def test_the_probe_run_writes_the_files_the_stage_reads(self, tmp_path):
        """The output names are verified against a real run, not remembered.

        `common.py:1543-1551` maps input names to output names by string
        substitution; a wrong constant is a stage that reports "no
        recombination" because it looked under the wrong name.
        """
        binary = gubbins_adapter.resolve_binary([str(SOURCE_GUBBINS.parent)])
        runner = _runner_or_skip()
        result = gubbins_adapter.run_gubbins(
            _write_minimal_alignment(tmp_path, sites=1200),
            tmp_path / "work",
            runner=runner,
            binary=binary,
            threads=1,
            search_dirs=[str(Path(binary).parent)],
        )
        for path in (
            result.snp_alignment,
            result.final_tree,
            result.node_labelled_tree,
            result.per_branch_statistics,
        ):
            assert path.is_file(), f"gubbins did not write {path.name}"
        assert result.exit_code == 0
        assert result.seconds > 0
        assert result.binary_version == gubbins_adapter.SOURCE_BUILD_VERSION

    def test_the_masked_alignment_is_variable_sites_only(self, tmp_path):
        """The name `core_snp_alignment.fasta` promises polymorphic sites only.

        stage 10 reports `snp_density_vs_core` from this file, and a file
        carrying invariant columns would make that ratio meaningless.
        """
        binary = gubbins_adapter.resolve_binary([str(SOURCE_GUBBINS.parent)])
        runner = _runner_or_skip()
        result = gubbins_adapter.run_gubbins(
            _write_minimal_alignment(tmp_path, sites=1200),
            tmp_path / "work",
            runner=runner,
            binary=binary,
            threads=1,
            search_dirs=[str(Path(binary).parent)],
        )
        records, name = {}, None
        for line in result.snp_alignment.read_text(encoding="utf-8").splitlines():
            if line.startswith(">"):
                name = line[1:]
                records[name] = ""
            elif line.strip():
                records[name] += line.strip()
        assert records
        widths = {len(v) for v in records.values()}
        assert len(widths) == 1, "the masked alignment must be aligned"
        assert set("".join(records.values())) <= set("ACGT"), (
            "gubbins' masked alignment should carry no gaps or ambiguity codes"
        )
        # Every column must vary, or it is not a polymorphic-sites alignment.
        first = next(iter(records.values()))
        for index in range(len(first)):
            assert len({seq[index] for seq in records.values()}) > 1, (
                f"column {index} is invariant, so this is not a "
                "polymorphic-sites-only alignment"
            )


# --------------------------------------------------------------------------
# Helpers.
# --------------------------------------------------------------------------


def _config():
    """The real laptop config, loaded through the real loader.

    Not a hand-rolled stub: `PipelineConfig` has fifteen required constructor
    arguments, and a stub would let this file pass against a config shape no
    overlay ever has.
    """
    from papipeline.config.loader import load_config

    return load_config(REPO / "config" / "science.yaml", machine="laptop")


def _manifest():
    from papipeline.manifest import SampleManifest

    return SampleManifest(samples=())


def _write_minimal_alignment(tmp_path: Path, sites: int = 420) -> Path:
    """A small already-aligned FASTA. `sites` must divide by 60 for uniformity."""
    import random

    rng = random.Random(1)
    base = [rng.choice("ACGT") for _ in range(sites)]
    path = tmp_path / "core_gene_alignment_filtered.aln"
    with path.open("w", encoding="utf-8") as handle:
        for index in range(5):
            seq = list(base)
            for _ in range(sites // 10):
                position = rng.randrange(sites)
                seq[position] = rng.choice([c for c in "ACGT" if c != seq[position]])
            handle.write(f">s{index + 1}\n")
            handle.write("\n".join("".join(seq[i:i + 60]) for i in range(0, sites, 60)))
            handle.write("\n")
    return path


def _runner_or_skip() -> Path:
    """`run_gubbins.py`, or skip with the reason it could not be found.

    Skipping rather than failing is correct here: the driver lives in the gubbins
    conda environment, which is a separate artefact from the source-built binary
    and is not installed on every machine that has the binary. A machine without
    the driver has not got a broken stage.
    """
    try:
        return gubbins_adapter.resolve_runner([str(SOURCE_GUBBINS.parent)])
    except gubbins_adapter.ToolNotAvailableError as exc:
        pytest.skip(str(exc))

# --------------------------------------------------------------------------
# 6. The Snakefile edge. Dry-run coverage of this is blocked by other rules.
# --------------------------------------------------------------------------
#
# `snakemake -n --config mode=REAL` cannot reach rule `recombination` on this
# machine: every earlier REAL rule (amr, mlst, ...) demands real tool output
# that is not on disk, so Snakemake raises MissingInputException for rule `amr`
# while building the DAG and never evaluates ours. That is honest — it is why the
# mode-splitting is asserted against the helper directly here.


def _stage_inputs_function():
    """`stage_inputs`, lifted out of the Snakefile as source and executed.

    A Snakefile is not valid Python - `rule all:` is Snakemake DSL - so the file
    cannot be imported and `ast.parse` raises on it. The function body *is*
    Python, so it is sliced out textually and executed with the two module-level
    names it closes over, re-bound per mode.

    This asserts against the real source rather than a copy: a copy would drift
    the moment someone edited the Snakefile, which is the whole failure mode
    these checks exist to catch.
    """
    source = SNAKEFILE.read_text(encoding="utf-8")
    start = re.search(r"^def stage_inputs\(", source, re.MULTILINE)
    assert start, "stage_inputs not found in the Snakefile; has it been renamed?"
    # `source[start.end():]`, not `rest[start.end():]` - `end` is an absolute
    # offset into `source`, so re-basing it onto an already-sliced string
    # silently lands hundreds of characters away and extracts the wrong region.
    rest = source[start.end():]
    # Slice to the NEXT top-level definition, not to the first `return`: the
    # function has an early return for STUB, and stopping there would exec a
    # truncated copy that only ever exercises one mode.
    following = re.search(r"^(def |CONFIG_FILES)", rest, re.MULTILINE)
    body = rest[: following.start()] if following else rest
    # `rest` begins *after* `def stage_inputs(`, so the signature has to be put
    # back for the slice to be a compilable function.
    definition = "def stage_inputs(" + body
    assert "real_only" in definition and "test_only" in definition, (
        "the extracted stage_inputs does not mention both mode-specific slots; "
        "the slice is wrong, not the Snakefile"
    )
    assert "CONFIG_FILES" in definition, (
        "stage_inputs no longer appends CONFIG_FILES, so a config edit would not "
        "invalidate a run"
    )

    def call(mode: str, produced, external=(), test_only=(), real_only=()):
        namespace = {
            "RunMode": RunMode,
            "RESOLVED_MODE": RunMode[mode],
            "CONFIG_FILES": ["config/science.yaml"],
        }
        exec(compile(definition, str(SNAKEFILE), "exec"), namespace)
        return namespace["stage_inputs"](
            list(produced), list(external),
            test_only=list(test_only), real_only=list(real_only),
        )

    return call


def _rule_body(name: str) -> str:
    """One rule's text with its docstring removed.

    Removing the docstring is not cosmetic. This rule's docstring *names*
    `pangenome_summary.tsv` while explaining why it no longer depends on it, so a
    search over the raw block matches its own explanation. The repository's own
    consistency test separates the two for the same reason.
    """
    source = SNAKEFILE.read_text(encoding="utf-8")
    start = re.search(rf"^rule {name}:\s*$", source, re.MULTILINE)
    assert start, f"rule {name} not found in the Snakefile"
    rest = source[start.end():]
    following = re.search(r"^rule [A-Za-z_][A-Za-z0-9_]*:\s*$", rest, re.MULTILINE)
    block = rest[: following.start()] if following else rest
    doc = re.search(r'"""(.*?)"""', block, re.DOTALL)
    return block.replace(doc.group(0), "") if doc else block


class TestTheSnakefileInputEdge:
    def test_the_alignment_edge_is_real_only(self):
        """REAL demands the alignment; TEST does not.

        Panaroo runs in REAL only - in TEST the pangenome partition is built from
        the annotation table - so no core gene alignment exists under the TEST
        intermediate root. Demanding one there would make `snakemake -n` fail on
        a file that by construction cannot be there, turning the dry run into a
        check that the TEST tree carries a panaroo product.
        """
        stage_inputs = _stage_inputs_function()
        alignment = ["intermediate/panaroo/core_gene_alignment_filtered.aln"]
        assert alignment[0] in stage_inputs("REAL", [], real_only=alignment)
        assert alignment[0] not in stage_inputs("TEST", [], real_only=alignment)

    def test_the_standin_edge_is_test_only_and_still_there(self):
        """The stand-in remains, because the stage is still in `UNBUILT_STAGES`.

        Demanding it in REAL would pass a fabricated table off as the stage's
        result, which is exactly what `StandInError` exists to refuse.
        """
        stage_inputs = _stage_inputs_function()
        standin = ["standins/unbuilt_stages/recombination.tsv"]
        assert standin[0] in stage_inputs("TEST", [], test_only=standin)
        assert standin[0] not in stage_inputs("REAL", [], test_only=standin)

    def test_stub_declares_no_file_inputs_at_all(self):
        """STUB fabricates everything, so external files must not be demanded.

        Otherwise the DAG cannot resolve, which is the one thing STUB exists to
        prove.
        """
        stage_inputs = _stage_inputs_function()
        assert stage_inputs(
            "STUB", [], external=["metadata.tsv"], test_only=["standin.tsv"],
            real_only=["core_aln"],
        ) == ["config/science.yaml"], (
            "STUB must drop `external` as well as the mode-specific slots: it "
            "fabricates every output and runs with no fixture on disk"
        )

    def test_config_files_are_declared_in_every_mode(self):
        """Otherwise editing a config file would not invalidate a run."""
        stage_inputs = _stage_inputs_function()
        for mode in ("STUB", "TEST", "REAL"):
            assert "config/science.yaml" in stage_inputs(mode, [])

    def test_the_rule_demands_the_alignment_and_not_a_metric_tally(self):
        """The defect this replaced.

        The rule consumed `pangenome_summary.tsv` - a `metric`/`value` table with
        no sequence in it. A gubbins run can never be satisfied by that, so the
        rule gated on a file that said nothing about the thing it gated, and a DAG
        that built would have looked correct while the stage had no input.
        """
        body = _rule_body("recombination")
        assert "pangenome_summary.tsv" not in body
        assert "OUT_PANGENOME" not in body
        assert "CORE_ALIGNMENT_IN" in body

    def test_the_alignment_basename_has_exactly_one_definition(self):
        """The name is the adapter's to declare, not the workflow's.

        A literal in the Snakefile would be a second definition of a location the
        code already owns - how `TREE_IN`/`TREE_META_IN` drifted once already.
        """
        source = SNAKEFILE.read_text(encoding="utf-8")
        assert "core_gene_alignment_filtered.aln" not in source
        assert gubbins_adapter.CORE_ALIGNMENT_BASENAME == (
            "core_gene_alignment_filtered.aln"
        )
        assert "CORE_ALIGNMENT_BASENAME" in source

    def test_the_pangenome_subdirectory_name_has_exactly_one_definition(self):
        """`"panaroo"` is `adapters.panaroo.OUTPUT_DIRNAME`'s to declare."""
        from papipeline.adapters.panaroo import OUTPUT_DIRNAME

        assert OUTPUT_DIRNAME == "panaroo"
        assert "PANGENOME_DIRNAME" in SNAKEFILE.read_text(encoding="utf-8")


# --------------------------------------------------------------------------
# Did gubbins converge?
# --------------------------------------------------------------------------
class TestConvergenceIsRecorded:
    """A tree from a non-converged run has to be identifiable as such.

    **Why this is not already covered by the exit code.** gubbins exits 0 both
    when it converges and when it hits its iteration ceiling, and it writes the
    final tree, the masked alignment and the per-branch statistics either way.
    A run that stopped at the ceiling therefore produces a complete-looking
    result set that a downstream consumer cannot distinguish from a converged
    one — which is precisely how a truncated analysis becomes a reported
    finding. The only place the verdict appears is the console output, so
    nothing downstream can see it unless it is parsed out and recorded.

    The real 10-isolate run hit this: ``Maximum number of iterations (5)
    reached.`` with ``exit_code: 0``.
    """

    CEILING = "Checking for convergence...\nMaximum number of iterations (5) reached.\n"
    CONVERGED = "Checking for convergence...\nConvergence after 3 iterations: Tree observed before.\n"

    def test_hitting_the_ceiling_is_recorded_as_not_converged(self):
        verdict = gubbins_adapter.parse_convergence(self.CEILING)
        assert verdict.converged is False
        assert verdict.iterations == 5
        assert verdict.iterations_max == 5
        assert verdict.recordable is True

    def test_convergence_is_recorded_as_converged_with_its_iteration(self):
        verdict = gubbins_adapter.parse_convergence(self.CONVERGED)
        assert verdict.converged is True
        assert verdict.iterations == 3
        assert verdict.recordable is True

    def test_a_log_that_says_neither_is_unknown_and_not_false(self):
        """A crashed run is not a run that converged.

        `False` is a finding and `None` is the absence of one. Reporting the
        absence of a verdict as "did not converge" would be defensible, but
        reporting it as `False` alongside a real ceiling result would make the
        two indistinguishable in the record.
        """
        for text in ("", "Running Gubbins to detect recombinations...\nTraceback (most recent call last):"):
            verdict = gubbins_adapter.parse_convergence(text)
            assert verdict.converged is None
            assert verdict.recordable is False
            assert verdict.iterations is None
            assert verdict.iterations_max is None

    def test_an_unparseable_log_does_not_raise(self):
        """Provenance is written on failure too, so it must not be a second failure.

        `run_gubbins` writes `provenance.json` on both paths (contract section
        5.4). A parser that raised on a truncated log would replace the run's
        real error with a parse error and destroy the question provenance
        exists to answer.
        """
        verdict = gubbins_adapter.parse_convergence("Maximum number of iterations")
        assert verdict.converged is None

    def test_the_ceiling_wins_when_both_lines_are_present(self):
        """A resumed run can carry an earlier attempt's convergence line.

        `--resume` continues from a previous run's tree, so one stdout can hold
        a convergence from the earlier attempt and a ceiling from this one. The
        weaker claim is the one about the run that actually produced the
        outputs, so it is the one that must survive.
        """
        text = self.CONVERGED + self.CEILING
        assert gubbins_adapter.parse_convergence(text).converged is False

    def test_provenance_carries_the_verdict(self, tmp_path):
        """The record is the whole point; a parser nobody reads is not a fix."""
        result = gubbins_adapter.GubbinsResult(
            work_dir=tmp_path, prefix="recombination",
            snp_alignment=tmp_path / "snp", final_tree=tmp_path / "tre",
            node_labelled_tree=tmp_path / "nl",
            per_branch_statistics=tmp_path / "pbs",
            log=tmp_path / "log", binary=tmp_path / "gubbins",
            binary_version=gubbins_adapter.SOURCE_BUILD_VERSION,
            isolated_bin_dir=tmp_path / "bin",
            argv=["run_gubbins.py"], exit_code=0, seconds=22.07,
            stdout=self.CEILING,
            convergence=gubbins_adapter.parse_convergence(self.CEILING),
        )
        record = result.provenance()
        assert record["converged"] is False
        assert record["iterations"] == 5
        assert record["iterations_max"] == 5
        assert "iteration ceiling" in record["convergence_note"]

        written = gubbins_adapter.write_provenance(result, tmp_path / "prov.json")
        on_disk = json.loads(written.read_text())
        assert on_disk["converged"] is False
        assert on_disk["iterations"] == 5

    def test_a_result_with_no_verdict_records_unknown_not_false(self, tmp_path):
        """`provenance()` must not invent a verdict from an empty stdout.

        A `GubbinsResult` built by hand — which the tests and any future caller
        both do — has no stdout, and the default has to say so rather than
        assert non-convergence it never observed.
        """
        result = gubbins_adapter.GubbinsResult(
            work_dir=tmp_path, prefix="recombination",
            snp_alignment=tmp_path / "snp", final_tree=tmp_path / "tre",
            node_labelled_tree=tmp_path / "nl",
            per_branch_statistics=tmp_path / "pbs",
            log=tmp_path / "log", binary=tmp_path / "gubbins",
            binary_version=gubbins_adapter.SOURCE_BUILD_VERSION,
            isolated_bin_dir=tmp_path / "bin",
        )
        record = result.provenance()
        assert record["converged"] is None
        assert record["iterations"] is None
        assert record["convergence_note"] == gubbins_adapter.UNKNOWN_CONVERGENCE.note

    def test_run_gubbins_reads_the_verdict_off_its_own_stdout(self, tmp_path):
        """End to end through the adapter, so the wiring is what is tested.

        The parser being correct and the adapter not calling it are different
        bugs, and only this asserts the second.
        """
        source_dir = tmp_path / "src"
        source_dir.mkdir()
        # A stand-in runner: prints gubbins' two terminal lines, writes the SNP
        # alignment the adapter insists on, exits 0.
        runner = source_dir / "run_gubbins.py"
        runner.write_text(
            textwrap.dedent(
                """\
                #!/usr/bin/env python3
                # gubbins writes relative to its CWD, which the adapter sets to
                # work_dir, so a bare relative name is the whole contract here.
                import pathlib
                print("Checking for convergence...")
                print("Maximum number of iterations (5) reached.")
                pathlib.Path("recombination.filtered_polymorphic_sites.fasta").write_text(">a\\nAC\\n")
                """
            ),
            encoding="utf-8",
        )
        runner.chmod(0o755)
        result = gubbins_adapter.run_gubbins(
            _write_minimal_alignment(tmp_path),
            tmp_path / "work",
            runner=runner,
            binary=_fake_tool(source_dir, "gubbins"),
            threads=1,
            selftest=lambda: (0, gubbins_adapter.SOURCE_BUILD_VERSION),
        )
        assert result.convergence.converged is False
        assert result.convergence.iterations == 5
        assert result.provenance()["converged"] is False
