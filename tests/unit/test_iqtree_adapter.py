"""Stage 10's REAL caller: IQ-TREE, from a core SNP alignment.

`stages.phylogeny.run` reads `tree.nwk` and a `tree_metadata.tsv` crosswalk and
had no producer for either, so a REAL run had no tree at all - the same shape of
gap stages 3 and 4 had. IQ-TREE 3.1.3 is installed, so this is missing code rather
than a missing dependency.

**The pangenome step above it cannot run here.** `phylogeny.aligner` is
`panaroo`, and panaroo is absent - `laptop.yaml` records it `available: false`
with the reason (no osx-arm64 build at any version). `gubbins`, which feeds the
recombination-filtered alignment, is absent too. So the core alignment cannot be
built on this machine and this adapter starts from an alignment that already
exists. That limit is stated rather than papered over: a caller that cannot run
end to end should say which step it needs, not imply the rest works.

**Every flag is one this IQ-TREE actually has**, taken from `iqtree --help` on
the installed binary rather than from IQ-TREE 2 habits:

* `-T` takes the thread count. `-nt` is IQ-TREE 1 and is not accepted;
* `-B 1000` is ultrafast bootstrap, and the help states `>= 1000`;
* `--alrt` is SH-aLRT, and `--seed` makes a stochastic search reproducible;
* `--prefix` sets every output name at once.

**The tip crosswalk is derived from the tree, not assumed.** IQ-TREE names tips
by the sequence identifier in the alignment, so in principle tip == sample_id.
Writing that as an assumption would be a crosswalk that is right until it is
wrong; `tree_metadata.tsv` exists precisely so a mismatch is *representable*,
and `stages.phylogeny.validate_tree_samples` requires the tip set to equal the
manifest exactly. So the labels are read out of the Newick and written as found.

**The stage's REAL branch needs stage 3's output too.** `build_tree_outputs`
writes the contract's five columns, and `lineage_label` is the MLST sequence
type (ruling R4), so it reads `mlst/mlst_results.tsv` under the run's
tool-output root and refuses, naming that path, when it is absent. `real_stage3`
supplies it: a REAL tree stage with no stage-3 table is refused, and three tests
in this file were relying on exactly that shape before it was.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from papipeline.adapters import iqtree as adapter
from papipeline.config.loader import PipelineConfig
from papipeline.models import RunMode

REPO = Path(__file__).resolve().parents[2]
HAS_IQTREE = shutil.which("iqtree") is not None


@pytest.fixture
def alignment(tmp_path):
    """A small but real DNA alignment, written by this test rather than
    committed - the committed fixture is 5 columns long, which is too few sites
    for the configured model to mean anything."""
    ids = [f"S{i:03d}" for i in range(6)]
    # A deliberately structured alignment: two groups differing at fixed sites,
    # so the tree is not arbitrary and a support value is not the only signal.
    rows = [
        "ACGTACGTAAGGCCTTACGGAT",
        "ACGTACGTAAGGCCTTACGGAT",
        "ACGTACGTAAGGCCTTACGGAT",
        "ACGTTCGTAAGGCCTTACGGAT",
        "ACGTTCGTAAGGCCTTACGGAT",
        "ACGTTCGTAAGGCCTTACGGAT",
    ]
    path = tmp_path / "core_snp_alignment.fasta"
    path.write_text(
        "\n".join(f">{i}\n{r}" for i, r in zip(ids, rows)) + "\n", encoding="utf-8"
    )
    return path, ids


#: The columns `stages/mlst.py` `_write_mlst_results` emits, and `load_mlst`
#: requires. Copied from `docs/data_contract.md` (`mlst/mlst_results.tsv`)
#: rather than invented: a fixture one column short fails in `read_tsv` with an
#: error about the FIXTURE, which says nothing about the stage-3 seam.
MLST_HEADER = "sample_id\tST\talleles\tMLST_status\tmlst_scheme\tallele_database"


def write_mlst_results(path, calls):
    """`path` as stage 3 writes it: one row per sample, `ST` or `.`."""
    lines = [MLST_HEADER]
    for sample_id, st in calls:
        alleles = ";".join(f"locus{i}:1" for i in range(7)) if st else "."
        lines.append(
            f"{sample_id}\t{st or '.'}\t{alleles}\t"
            f"{'typed' if st else 'no_call'}\tpaerMLST\tpubmlst.org"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


@pytest.fixture
def real_stage3(tmp_path, monkeypatch):
    """A REAL run whose stage-3 output exists, isolated under `tmp_path`.

    `PIPELINE_RESULTS_ROOT` is the redirect the loader reads in exactly one
    place, so setting it moves `tool_output_root(REAL)` without touching a
    tracked overlay - the same reason every other REAL-mode test in this repo
    uses it. Returns a writer taking `{sample_id: st}`.
    """
    from papipeline.config.loader import RESULTS_ROOT_ENV

    monkeypatch.setenv(RESULTS_ROOT_ENV, str(tmp_path / "results"))
    target = (
        tmp_path / "results" / RunMode.REAL.value.lower()
        / "intermediate" / "mlst" / "mlst_results.tsv"
    )

    def write(calls):
        return write_mlst_results(target, list(calls.items()))

    write.path = target
    return write


class TestTheCommandComesFromConfig:
    def test_the_model_comes_from_config(self, config: PipelineConfig):
        command = adapter.build_command(
            "iqtree", alignment=Path("a.fasta"), prefix=Path("/tmp/x"),
            model=config.raw["phylogeny"]["models"],
            threads=int(config.runtime["threads"]), seed=12345,
        )
        assert "-m" in command
        assert command[command.index("-m") + 1] == config.raw["phylogeny"]["models"] == "GTR+G+ASC", (
            "the configured model must reach -m verbatim. This asserted the "
            "substring 'GTR+G', which 'GTR+G+ASC' also contains - so the "
            "assertion kept passing while the model it was guarding changed."
        )

    def test_threads_use_dash_upper_t_not_dash_nt(self, config: PipelineConfig):
        """`-nt` is IQ-TREE 1 and is not accepted by the installed binary."""
        command = adapter.build_command(
            "iqtree", alignment=Path("a.fasta"), prefix=Path("/tmp/x"),
            model="GTR+G", threads=4, seed=1,
        )
        assert "-T" in command
        assert "-nt" not in command

    def test_support_replicates_are_at_least_the_documented_floor(self):
        """The help states `-B ... (>=1000)`, so a smaller value is invalid
        rather than merely weak."""
        command = adapter.build_command(
            "iqtree", alignment=Path("a.fasta"), prefix=Path("/tmp/x"),
            model="GTR+G", threads=1, seed=1, bootstrap=1000, alrt=1000,
        )
        assert command[command.index("-B") + 1] == "1000"
        assert command[command.index("--alrt") + 1] == "1000"

    def test_a_seed_is_always_passed(self, config: PipelineConfig):
        """A stochastic search without a seed is not reproducible, and
        reproducibility is the reason the seed exists."""
        command = adapter.build_command(
            "iqtree", alignment=Path("a.fasta"), prefix=Path("/tmp/x"),
            model="GTR+G", threads=1, seed=4242,
        )
        assert "--seed" in command
        assert "4242" in command

    def test_the_prefix_sets_every_output_name(self, config: PipelineConfig):
        command = adapter.build_command(
            "iqtree", alignment=Path("a.fasta"), prefix=Path("/tmp/out/tree"),
            model="GTR+G", threads=1, seed=1,
        )
        assert "--prefix" in command
        assert "/tmp/out/tree" in command

    def test_nothing_is_hard_coded(self):
        """The model must come from config, not from the module."""
        command = adapter.build_command(
            "iqtree", alignment=Path("a.fasta"), prefix=Path("/tmp/x"),
            model="TEST-MODEL", threads=1, seed=1,
        )
        assert "GTR+G" not in command


@pytest.mark.skipif(not HAS_IQTREE, reason="iqtree is not on PATH")
class TestAgainstTheRealTool:
    def test_it_produces_a_tree_whose_tips_are_the_alignment_ids(
        self, alignment, tmp_path
    ):
        from papipeline.stages.phylogeny import extract_tips

        path, ids = alignment
        result = adapter.build_tree(
            path, tmp_path, model="GTR+G", threads=2, seed=12345,
            bootstrap=1000, alrt=1000,
        )
        assert result.tree.is_file(), "IQ-TREE produced no .treefile"
        assert sorted(extract_tips(result.tree.read_text())) == sorted(ids)

    def test_the_crosswalk_is_derived_from_the_tree(self, alignment, tmp_path):
        """Not assumed. If the alignment renamed a sequence, the crosswalk must
        show the real tip, or `validate_tree_samples` fails on a file that
        claims to agree."""
        path, ids = alignment
        renamed = tmp_path / "renamed.fasta"
        # Every sequence, so the assertion below is about the whole crosswalk
        # rather than about one renamed tip among unchanged ones.
        renamed.write_text(
            path.read_text().replace(">S", ">tip_S"), encoding="utf-8"
        )
        result = adapter.build_tree(
            renamed, tmp_path, model="GTR+G", threads=2, seed=1,
            bootstrap=1000, alrt=1000,
        )
        crosswalk = dict(result.crosswalk)
        assert crosswalk, "no crosswalk rows were produced"
        assert all(tip.startswith("tip_") for tip in crosswalk.values()), (
            f"crosswalk claims tips are unchanged: {crosswalk}"
        )

    def test_a_missing_alignment_is_refused_before_launching_anything(
        self, tmp_path
    ):
        """Cheap and specific, rather than a tool error about a file."""
        with pytest.raises(Exception) as excinfo:
            adapter.build_tree(
                tmp_path / "absent.fasta", tmp_path, model="GTR+G",
                threads=1, seed=1, bootstrap=1000, alrt=1000,
            )
        assert "absent.fasta" in str(excinfo.value)

    def test_a_nonzero_exit_is_refused(self, tmp_path):
        from papipeline.errors import ToolExecutionError

        broken = tmp_path / "broken.fasta"
        broken.write_text("not a fasta at all\n", encoding="utf-8")
        with pytest.raises(ToolExecutionError):
            adapter.build_tree(
                broken, tmp_path, model="GTR+G", threads=1, seed=1,
                bootstrap=1000, alrt=1000,
            )


class TestTheAbsentUpstreamStepIsNamed:
    def test_the_aligner_is_named_and_its_absence_is_not_hidden(self):
        """`phylogeny.aligner` is panaroo, which is unavailable here.

        The adapter starts from an alignment, so on a machine without panaroo it
        cannot be the whole story. A caller that silently produced nothing would
        read as "stage 10 ran and found no tree"; refusing names the step.
        """
        assert adapter.REQUIRED_ALIGNMENT_KEY == "phylogeny.aligner"
        assert "panaroo" in adapter.REQUIRED_ALIGNMENT_TOOL

@pytest.mark.skipif(not HAS_IQTREE, reason="iqtree is not on PATH")
class TestTheStageRealBranchProducesBothFiles:
    """The adapter is only useful if the stage calls it and writes what it reads.

    `stages.phylogeny.run` reads `tree.nwk` and `tree_metadata.tsv`; an adapter
    that wrote them somewhere else would leave the stage exactly as short of a
    tree as before.
    """

    def test_it_writes_tree_and_crosswalk_under_the_names_the_stage_reads(
        self, alignment, tmp_path, config: PipelineConfig, real_stage3
    ):
        import copy

        from papipeline.manifest import SampleManifest
        from papipeline.models import Sample
        from papipeline.stages import phylogeny as stage

        path, ids = alignment
        # `lineage_label` is the MLST sequence type (R4), so the producer reads
        # stage 3's table; absent, it refuses rather than labelling every
        # sample "unknown".
        real_stage3({i: "235" for i in ids})
        phylo_dir = tmp_path / "phylogeny"
        phylo_dir.mkdir()
        shutil.copy(path, phylo_dir / "core_snp_alignment.fasta")

        prepared = copy.copy(config)
        raw = copy.deepcopy(dict(config.raw))
        raw["phylogeny"] = dict(raw.get("phylogeny") or {})
        raw["phylogeny"]["require_exact_sample_match"] = False
        # This fixture has ONE polymorphic column, so the configured +ASC model
        # and the filter cannot both apply to it: filtered it is a single site,
        # and IQ-TREE exits 2 on a single site (`Unknown sequence type`). The
        # configured pair is exercised by
        # `TestAscIsUsableOnceTheFilterHasRun` and
        # `TestTheStagePassesTheSwitchAndRecordsTheCount`; this test is about
        # the stage calling the adapter and writing both files.
        raw["phylogeny"]["models"] = "GTR+G"
        raw["phylogeny"]["asc_drop_partially_constant"] = False
        object.__setattr__(prepared, "raw", raw)

        manifest = SampleManifest(samples=[Sample(sample_id=i) for i in ids])
        stage.run(prepared, manifest, RunMode.REAL, phylo_dir)

        assert (phylo_dir / "tree.nwk").is_file()
        assert (phylo_dir / "tree_metadata.tsv").is_file()

    def test_it_refuses_when_the_alignment_is_absent(self, tmp_path, config):
        """Naming panaroo beats reporting an empty tree.

        With no alignment the stage would otherwise summarise a missing file and
        read as "stage 10 ran, no phylogeny" - an absence of signal where there
        was a missing input.
        """
        from papipeline.manifest import SampleManifest
        from papipeline.models import Sample
        from papipeline.stages import phylogeny as stage

        manifest = SampleManifest(samples=[Sample(sample_id="S000")])
        with pytest.raises(Exception) as excinfo:
            stage.run(config, manifest, RunMode.REAL, tmp_path)
        assert "panaroo" in str(excinfo.value)

    def test_an_existing_tree_is_not_rebuilt(self, alignment, tmp_path, config):
        """Re-running must not overwrite a tree that is already there.

        IQ-TREE is stochastic in its search even with a seed, so a rebuild would
        replace a result with a merely-similar one for no reason.
        """
        from papipeline.manifest import SampleManifest
        from papipeline.models import Sample
        from papipeline.stages import phylogeny as stage

        path, ids = alignment
        phylo_dir = tmp_path / "phylogeny"
        phylo_dir.mkdir()
        # Every tip, because `validate_tree_samples` requires the tip set to
        # equal the manifest - a short tree is refused before the rebuild check
        # is ever reached.
        nested = ",".join(f"{i}:0.0" for i in ids[1:])
        (phylo_dir / "tree.nwk").write_text(
            f"({ids[0]}:0.0,({nested}):0.0);\n", encoding="utf-8"
        )
        manifest = SampleManifest(samples=[Sample(sample_id=i) for i in ids])
        # No alignment present, so a rebuild attempt would raise. It does not.
        stage.run(config, manifest, RunMode.TEST, phylo_dir)
        assert (phylo_dir / "tree.nwk").read_text().count(";") == 1


# ---------------------------------------------------------------------------
# R8: +ASC, on an alignment whose constant-after-gaps columns have gone
# ---------------------------------------------------------------------------


class TestTheConstantAfterGapsFilter:
    """`adapters.iqtree.drop_partially_constant_columns`, which is what makes
    ``-m GTR+G+ASC`` runnable at all.

    The pinned binary refuses the model while any column is constant once gaps
    and N are ignored, and says so by name::

        ERROR: Invalid use of +ASC because of 20 invariant sites in the alignment

    so the filter is not an optimisation. It is the precondition, and the two
    halves of R8 are asserted together in
    `tests/unit/test_round11_config_keys.py::TestPhylogenyModelAndAscPreprocessing`.
    """

    def test_a_partially_constant_column_is_dropped(self, tmp_path):
        """`A A A A -` is invariant for the correction even though it carries a
        gap, so a test that only counted ungapped invariant columns would leave
        it in place and the tool would still refuse."""
        source = tmp_path / "in.fasta"
        source.write_text(
            ">S1\nACGT\n>S2\nACGT\n>S3\nACGT\n>S4\n-CG-\n>S5\nAANT\n",
            encoding="utf-8",
        )
        result = adapter.drop_partially_constant_columns(
            source, tmp_path / "out.fasta"
        )
        assert result.n_columns_in == 4
        assert result.n_columns_kept == 1, "only the polymorphic column may remain"
        assert result.n_columns_dropped == 3
        assert (tmp_path / "out.fasta").read_text() == (
            ">S1\nC\n>S2\nC\n>S3\nC\n>S4\nC\n>S5\nA\n"
        )

    def test_a_fully_gap_filled_column_is_dropped(self, tmp_path):
        """No base call at all is a stronger version of the same problem."""
        source = tmp_path / "in.fasta"
        source.write_text(">S1\nAC-\n>S2\nAC-\n>S3\nGT-\n", encoding="utf-8")
        result = adapter.drop_partially_constant_columns(
            source, tmp_path / "out.fasta"
        )
        assert (result.n_columns_in, result.n_columns_kept) == (3, 2)
        assert result.n_columns_dropped == 1, (
            "the all-gap column carries no information about base composition, "
            "so leaving it in is leaving in a column +ASC cannot see"
        )

    def test_the_sequence_order_and_ids_are_preserved(self, tmp_path):
        """The filter may not reorder or rename: the tip set is validated
        against the manifest downstream, so a reordering is invisible here and
        a rename is a hard failure there."""
        source = tmp_path / "in.fasta"
        source.write_text(
            ">first\nACGTA\n>second\nACGTA\n>third\nACGTT\n", encoding="utf-8"
        )
        result = adapter.drop_partially_constant_columns(
            source, tmp_path / "out.fasta"
        )
        assert result.n_columns_kept == 1
        assert (tmp_path / "out.fasta").read_text() == (
            ">first\nA\n>second\nA\n>third\nT\n"
        )

    def test_lowercase_and_question_marks_count_as_no_call(self, tmp_path):
        """`read_fasta` upper-cases residues; `?` is the other spelling of a
        missing call that appears in real alignments."""
        assert adapter.observed_states("Aa-?Nn") == {"A"}

    def test_an_alignment_with_nothing_variable_is_refused(self, tmp_path):
        """An empty filtered alignment would make the tool's failure name a file
        rather than the filter that emptied it, so the filter refuses first."""
        from papipeline.errors import DataContractError

        source = tmp_path / "in.fasta"
        source.write_text(">S1\nACG\n>S2\nACG\n", encoding="utf-8")
        with pytest.raises(DataContractError) as excinfo:
            adapter.drop_partially_constant_columns(source, tmp_path / "out.fasta")
        assert "constant" in str(excinfo.value)
        assert not (tmp_path / "out.fasta").exists(), (
            "a refused filter must not leave a truncated alignment behind"
        )

    def test_a_ragged_alignment_is_refused(self, tmp_path):
        from papipeline.errors import DataContractError

        source = tmp_path / "in.fasta"
        source.write_text(">S1\nACGT\n>S2\nACG\n", encoding="utf-8")
        with pytest.raises(DataContractError) as excinfo:
            adapter.drop_partially_constant_columns(source, tmp_path / "out.fasta")
        assert "rectangular" in str(excinfo.value)


@pytest.fixture
def asc_alignment(tmp_path):
    """An alignment with enough variable sites for the model to mean something,
    plus the two kinds of column the filter exists to remove.

    Twelve columns are genuinely polymorphic, one is constant once gaps are
    ignored (`A A A - A A`), and one carries no base call at all. The
    committed 20-column fixture cannot supply this: filtering it leaves a single
    site, and IQ-TREE exits on that (``It makes no sense to perform bootstrap
    with less than 4 sequences``).
    """
    ids = [f"S{i:03d}" for i in range(6)]
    columns = [
        "AAATTT", "CCCGGG", "GTAAAA", "TGGGGG", "AAAAAC", "CCTTTT",
        "GGGGGA", "TAAAAA", "ATTTCC", "CGGGGT", "GTTTTT", "TGACCA",
    ]
    partially_constant = "AAA-AA"
    uncalled = "N-N-N-"
    rows = [
        "".join(column[i] for column in columns)
        + partially_constant[i]
        + uncalled[i]
        for i in range(len(ids))
    ]
    path = tmp_path / "asc_alignment.fasta"
    path.write_text(
        "\n".join(f">{i}\n{r}" for i, r in zip(ids, rows)) + "\n", encoding="utf-8"
    )
    return path, ids


class TestAscIsUsableOnceTheFilterHasRun:
    """The behaviour R8 buys, against the real binary."""

    @pytest.mark.skipif(not HAS_IQTREE, reason="iqtree is not on PATH")
    def test_an_alignment_with_partially_constant_columns_does_not_crash(
        self, asc_alignment, tmp_path
    ):
        """The failure this prevents is `ERROR: Invalid use of +ASC because of
        20 invariant sites in the alignment` - exit 2, named, and total."""
        path, ids = asc_alignment
        result = adapter.build_tree(
            path, tmp_path, model="GTR+G+ASC", threads=2, seed=7,
            bootstrap=1000, alrt=1000, drop_partially_constant=True,
        )
        assert result.column_filter is not None
        assert result.column_filter.n_columns_dropped == 2, (
            "the two constant-after-gaps columns must go, or +ASC is refused"
        )
        assert result.tree.is_file()

    @pytest.mark.skipif(not HAS_IQTREE, reason="iqtree is not on PATH")
    def test_without_the_filter_the_same_alignment_is_refused(
        self, asc_alignment, tmp_path
    ):
        """The filter is what makes the model runnable, so the counterfactual is
        asserted rather than described: same alignment, same model, no filter."""
        from papipeline.errors import ToolExecutionError

        path, _ids = asc_alignment
        with pytest.raises(ToolExecutionError) as excinfo:
            adapter.build_tree(
                path, tmp_path, model="GTR+G+ASC", threads=2, seed=7,
                bootstrap=1000, alrt=1000, drop_partially_constant=False,
            )
        assert "ASC" in str(excinfo.value)

    @pytest.mark.skipif(not HAS_IQTREE, reason="iqtree is not on PATH")
    def test_the_search_reads_the_filtered_file_not_the_original(
        self, asc_alignment, tmp_path
    ):
        """A filter that computes the count and then points `-s` at the
        unfiltered alignment would pass every other test here."""
        path, _ids = asc_alignment
        result = adapter.build_tree(
            path, tmp_path, model="GTR+G+ASC", threads=2, seed=7,
            bootstrap=1000, alrt=1000, drop_partially_constant=True,
        )
        served = result.command[result.command.index("-s") + 1]
        assert served == str(result.column_filter.alignment)
        assert served != str(path)
        assert (tmp_path / adapter.FILTERED_ALIGNMENT_NAME).is_file()

    @pytest.mark.skipif(not HAS_IQTREE, reason="iqtree is not on PATH")
    def test_no_filter_means_no_filter_record(self, alignment, tmp_path):
        """`None` and "filtered, dropped nothing" must stay distinguishable, so
        an unfiltered run cannot report a count of zero it never measured."""
        path, _ids = alignment
        result = adapter.build_tree(
            path, tmp_path, model="GTR+G", threads=1, seed=1,
            bootstrap=1000, alrt=1000, drop_partially_constant=False,
        )
        assert result.column_filter is None
        assert not (tmp_path / adapter.FILTERED_ALIGNMENT_NAME).exists()


class TestTheStagePassesTheSwitchAndRecordsTheCount:
    """`phylogeny.asc_drop_partially_constant` must reach the adapter, and the
    count of dropped columns must be recorded rather than merely logged.

    The adapter is stubbed, so this asserts the wiring - which is the part that
    was inert - without needing a tree that a 6-tip fixture cannot support under
    a +ASC model.
    """

    def _prepared(self, config, **phylogeny_overrides):
        import copy

        raw = copy.deepcopy(dict(config.raw))
        raw["phylogeny"] = dict(raw.get("phylogeny") or {})
        raw["phylogeny"].update(phylogeny_overrides)
        prepared = copy.copy(config)
        object.__setattr__(prepared, "raw", raw)
        return prepared

    def _run_stage(self, monkeypatch, tmp_path, config, ids, real_stage3, **overrides):
        from papipeline.manifest import SampleManifest
        from papipeline.models import Sample
        from papipeline.stages import phylogeny as stage

        captured = {}

        def fake_build_tree(alignment, out_dir, **kwargs):
            captured["alignment"] = Path(alignment)
            captured.update(kwargs)
            out = Path(out_dir)
            out.mkdir(parents=True, exist_ok=True)
            (out / "iqtree.treefile").write_text(
                "(" + ",".join(f"{i}:0.0" for i in ids) + ");\n", encoding="utf-8"
            )
            return adapter.TreeResult(
                tree=out / "iqtree.treefile",
                report=None,
                log=None,
                tip_labels=list(ids),
                crosswalk={i: i for i in ids},
                command=["iqtree"],
                column_filter=adapter.ColumnFilterResult(
                    alignment=out / adapter.FILTERED_ALIGNMENT_NAME,
                    n_columns_in=70660,
                    n_columns_kept=70640,
                ),
            )

        monkeypatch.setattr(adapter, "build_tree", fake_build_tree)
        real_stage3({i: "235" for i in ids})
        phylo_dir = tmp_path / "phylogeny"
        phylo_dir.mkdir()
        (phylo_dir / "core_snp_alignment.fasta").write_text(
            ">x\nACGT\n", encoding="utf-8"
        )
        manifest = SampleManifest(samples=[Sample(sample_id=i) for i in ids])
        stage.run(self._prepared(config, **overrides), manifest, RunMode.REAL, phylo_dir)
        return captured, phylo_dir

    def test_the_config_key_reaches_the_adapter(
        self, monkeypatch, tmp_path, config, real_stage3
    ):
        captured, _ = self._run_stage(
            monkeypatch, tmp_path, config, ["S1", "S2"], real_stage3,
            require_exact_sample_match=False,
        )
        assert captured["drop_partially_constant"] is True, (
            "phylogeny.asc_drop_partially_constant is still inert: the stage did "
            "not pass it to the adapter"
        )
        assert captured["model"] == "GTR+G+ASC"

    def test_the_key_can_be_turned_off_and_then_it_is_off(
        self, monkeypatch, tmp_path, config, real_stage3
    ):
        captured, _ = self._run_stage(
            monkeypatch, tmp_path, config, ["S1", "S2"], real_stage3,
            require_exact_sample_match=False, asc_drop_partially_constant=False,
        )
        assert captured["drop_partially_constant"] is False

    def test_the_dropped_column_count_is_recorded_in_the_file(
        self, monkeypatch, tmp_path, config, real_stage3
    ):
        """A model claimed in a file with no count beside it has shown nothing."""
        _, phylo_dir = self._run_stage(
            monkeypatch, tmp_path, config, ["S1", "S2"], real_stage3,
            require_exact_sample_match=False,
        )
        header = (phylo_dir / "tree_metadata.tsv").read_text(encoding="utf-8")
        assert "alignment_columns_dropped=20" in header
        assert "alignment_columns_kept=70640" in header
