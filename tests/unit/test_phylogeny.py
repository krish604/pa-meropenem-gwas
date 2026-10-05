"""Tests for the phylogeny module (stage 10).

The tree/manifest correspondence check is the important one: a tree whose
tips do not match the manifest one-to-one silently corrupts every
lineage-aware stage downstream, so mismatches must be hard errors with
precise diagnostics.
"""

from __future__ import annotations

import pytest

from papipeline.errors import DataContractError, TreeSampleMismatchError
from papipeline.manifest import SampleManifest
from papipeline.models import Sample
from papipeline.stages import phylogeny as stage


class TestExtractTips:
    @pytest.mark.parametrize(
        "newick,expected",
        [
            ("(A,B);", ["A", "B"]),
            ("(A,B,(C,D));", ["A", "B", "C", "D"]),
            ("((A:0.1,B:0.2)0.9:0.3,(C:0.1,D:0.1)0.8:0.3)root;", ["A", "B", "C", "D"]),
            ("(TEST_PA_001,TEST_PA_002,TEST_PA_003);", ["TEST_PA_001", "TEST_PA_002", "TEST_PA_003"]),
        ],
    )
    def test_extracts_tips(self, newick, expected):
        assert stage.extract_tips(newick) == expected

    def test_internal_labels_excluded(self):
        tips = stage.extract_tips("((n1:0.1,A:0.2)n2:0.3,(B:0.1,n3:0.1)n4:0.3)n5;")
        assert sorted(tips) == ["A", "B"]

    def test_quoted_labels_preserved(self):
        tips = stage.extract_tips("('sample one':0.1,'sample two':0.2)n1;")
        assert sorted(tips) == ["sample one", "sample two"]

    def test_labels_with_spaces_unquoted(self):
        tips = stage.extract_tips("(sample one,sample two);")
        assert sorted(tips) == ["sample one", "sample two"]

    def test_no_branch_lengths(self):
        assert stage.extract_tips("(A,B,C);") == ["A", "B", "C"]

    def test_deeply_nested(self):
        assert sorted(stage.extract_tips("((((A,B),C),D),E);")) == list("ABCDE")

    def test_single_tip(self):
        assert stage.extract_tips("(A);") == ["A"]

    @pytest.mark.parametrize("bad", ["", "   ", "\n"])
    def test_empty_raises(self, bad):
        with pytest.raises(DataContractError, match="empty"):
            stage.extract_tips(bad)

    def test_unbalanced_raises(self):
        with pytest.raises(DataContractError):
            stage.extract_tips("((A,B);")

    def test_unterminated_quote_raises(self):
        with pytest.raises(DataContractError, match="Unterminated"):
            stage.extract_tips("('unclosed:0.1,B);")

    def test_root_label_after_clade_is_legal(self):
        """A trailing internal/root label is valid Newick, not an error."""
        assert stage.extract_tips("(A,B)@;") == ["A", "B"]

    def test_unexpected_character_inside_clade_raises(self):
        with pytest.raises(DataContractError, match="Unexpected character"):
            stage.extract_tips("(A,B];")


class TestValidateTreeSamples:
    def test_exact_match_passes(self, sample_manifest):
        stage.validate_tree_samples(sample_manifest.sample_ids, sample_manifest)

    def test_missing_tip_raises(self, sample_manifest):
        with pytest.raises(TreeSampleMismatchError) as exc:
            stage.validate_tree_samples(["TEST_A_01", "TEST_A_02"], sample_manifest)
        assert "TEST_A_03" in exc.value.context["missing_in_tree"]

    def test_extra_tip_raises(self, sample_manifest):
        tips = sample_manifest.sample_ids + ["TEST_A_99"]
        with pytest.raises(TreeSampleMismatchError) as exc:
            stage.validate_tree_samples(tips, sample_manifest)
        assert "TEST_A_99" in exc.value.context["not_in_manifest"]

    def test_duplicate_tip_raises(self, sample_manifest):
        tips = ["TEST_A_01", "TEST_A_01", "TEST_A_02", "TEST_A_03"]
        with pytest.raises(TreeSampleMismatchError, match="duplicate tip"):
            stage.validate_tree_samples(tips, sample_manifest)

    def test_both_missing_and_extra_reported(self, sample_manifest):
        with pytest.raises(TreeSampleMismatchError) as exc:
            stage.validate_tree_samples(
                ["TEST_A_01", "TEST_A_02", "GHOST"], sample_manifest
            )
        assert exc.value.context["n_missing"] == 1
        assert exc.value.context["n_extra"] == 1

    def test_order_does_not_matter(self, sample_manifest):
        reordered = list(reversed(sample_manifest.sample_ids))
        stage.validate_tree_samples(reordered, sample_manifest)


class TestSummariseTree:
    def test_validates_against_synthetic_manifest(self, config, manifest, phylogeny_dir):
        summary = stage.summarise_tree(
            config, phylogeny_dir / "tree.nwk", manifest
        )
        assert summary.n_tips == 20
        assert summary.has_branch_lengths

    def test_missing_tree_raises(self, config, manifest, tmp_path):
        with pytest.raises(DataContractError, match="not found"):
            stage.summarise_tree(config, tmp_path / "absent.nwk", manifest)

    def test_empty_tree_raises(self, config, manifest, tmp_path):
        path = tmp_path / "empty.nwk"
        path.write_text("")
        with pytest.raises(DataContractError, match="empty"):
            stage.summarise_tree(config, path, manifest)

    def test_mismatched_tree_raises(self, config, manifest, tmp_path):
        path = tmp_path / "wrong.nwk"
        path.write_text("(ONLY_ONE_TIP);\n")
        with pytest.raises(TreeSampleMismatchError):
            stage.summarise_tree(config, path, manifest)

    def test_can_be_disabled_by_config(self, config, manifest, tmp_path):
        path = tmp_path / "wrong.nwk"
        path.write_text("(ONLY_ONE_TIP);\n")
        relaxed = config.__class__(
            **{**config.__dict__, "phylogeny": config.phylogeny.__class__(
                **{**config.phylogeny.__dict__, "require_exact_sample_match": False}
            )}
        )
        summary = stage.summarise_tree(relaxed, path, manifest)
        assert summary.n_tips == 1

    def test_duplicate_tip_in_tree_file_raises(self, config, manifest, tmp_path):
        path = tmp_path / "dup.nwk"
        ids = manifest.sample_ids
        path.write_text("(" + ",".join(ids[:-1] + [ids[0]]) + ");\n")
        with pytest.raises(TreeSampleMismatchError, match="duplicate tip"):
            stage.summarise_tree(config, path, manifest)


class TestTreeMetadata:
    def test_loads_synthetic_metadata(self, phylogeny_dir):
        mapping = stage.load_tree_metadata(phylogeny_dir / "tree_metadata.tsv")
        assert len(mapping) == 20
        assert set(mapping.values()) == {"LINEAGE_A", "LINEAGE_B", "LINEAGE_C"}

    def test_missing_file_returns_empty(self, tmp_path):
        assert stage.load_tree_metadata(tmp_path / "absent.tsv") == {}


class TestAlignmentSummary:
    def test_summarises_both_alignments(self, phylogeny_dir):
        info = stage.summarise_alignments(
            phylogeny_dir / "core_alignment.fasta",
            phylogeny_dir / "core_snp_alignment.fasta",
        )
        assert info["core_n_sequences"] == 20
        assert info["snp_n_sites"] < info["core_n_sites"]

    def test_snp_density_is_a_fraction(self, phylogeny_dir):
        info = stage.summarise_alignments(
            phylogeny_dir / "core_alignment.fasta",
            phylogeny_dir / "core_snp_alignment.fasta",
        )
        assert 0 < info["snp_density_vs_core"] < 1

    def test_missing_alignment_is_tolerated(self, tmp_path):
        summary = stage.summarise_alignment(tmp_path / "absent.fasta")
        assert summary["n_sequences"] == 0


class TestRun:
    def test_returns_summary_and_alignment_info(self, config, manifest, phylogeny_dir):
        from papipeline.models import RunMode

        summary, info = stage.run(config, manifest, RunMode.TEST, phylogeny_dir)
        assert summary.n_tips == 20
        assert info["core_n_sequences"] == 20

    def test_alignment_info_empty_when_absent(self, config, manifest, tmp_path):
        from papipeline.models import RunMode

        newick = "(" + ",".join(manifest.sample_ids) + ");\n"
        (tmp_path / "tree.nwk").write_text(newick)
        _summary, info = stage.run(config, manifest, RunMode.TEST, tmp_path)
        assert info == {}


# ---------------------------------------------------------------------------
# R4: lineage_label is the MLST sequence type from stage 3
# ---------------------------------------------------------------------------


MLST_HEADER = "sample_id\tST\talleles\tMLST_status\tmlst_scheme\tallele_database"


@pytest.fixture
def stage3_table(tmp_path, monkeypatch):
    """Stage 3's contracted table, in an isolated run root.

    `PIPELINE_RESULTS_ROOT` is the redirect `loader.results_root` reads, in one
    place, so this moves `tool_output_root(REAL)` without touching an overlay.
    """
    from papipeline.config.loader import RESULTS_ROOT_ENV
    from papipeline.models import RunMode

    monkeypatch.setenv(RESULTS_ROOT_ENV, str(tmp_path / "results"))
    path = (
        tmp_path / "results" / RunMode.REAL.value.lower()
        / "intermediate" / "mlst" / "mlst_results.tsv"
    )

    def write(calls):
        path.parent.mkdir(parents=True, exist_ok=True)
        rows = [MLST_HEADER]
        for sample_id, st in calls.items():
            alleles = ";".join(f"locus{i}:{i + 1}" for i in range(7))
            rows.append(
                f"{sample_id}\t{st}\t{alleles}\ttyped\tpaerMLST\tpubmlst.org"
            )
        path.write_text("\n".join(rows) + "\n", encoding="utf-8")
        return path

    write.path = path
    return write


class TestLineageLabelsComeFromStageThree:
    def test_the_labels_are_the_sequence_types(
        self, config, manifest, stage3_table
    ):
        from papipeline.models import RunMode

        stage3_table({sid: str(100 + i) for i, sid in enumerate(manifest.sample_ids)})
        labels = stage.load_lineage_labels(config, RunMode.REAL, manifest)
        assert labels == {
            sid: str(100 + i) for i, sid in enumerate(manifest.sample_ids)
        }

    def test_it_reads_the_contracted_path_not_one_it_invents(
        self, config, manifest, stage3_table
    ):
        from papipeline.models import RunMode

        written = stage3_table({sid: "235" for sid in manifest.sample_ids})
        assert written == stage.mlst_results_path(config, RunMode.REAL), (
            "the producer must read the path docs/data_contract.md declares for "
            "stage 3; a second path is a second seam"
        )

    def test_a_missing_table_is_refused_naming_it(self, config, manifest, tmp_path,
                                                  monkeypatch):
        from papipeline.config.loader import RESULTS_ROOT_ENV
        from papipeline.errors import PipelineError
        from papipeline.models import RunMode

        monkeypatch.setenv(RESULTS_ROOT_ENV, str(tmp_path / "absent"))
        with pytest.raises(PipelineError) as excinfo:
            stage.load_lineage_labels(config, RunMode.REAL, manifest)
        assert str(stage.mlst_results_path(config, RunMode.REAL)) in str(excinfo.value)

    def test_an_untyped_sample_is_the_sentinel_not_a_dropped_row(
        self, config, manifest, stage3_table
    ):
        """The manifest decides cohort membership. A sample with no sequence type
        must still be present in the crosswalk, carrying a value the refusing
        stages can refuse on."""
        from papipeline.models import RunMode

        first = manifest.sample_ids[0]
        stage3_table({sid: "235" for sid in manifest.sample_ids if sid != first})
        labels = stage.load_lineage_labels(config, RunMode.REAL, manifest)
        assert set(labels) == set(manifest.sample_ids)
        assert labels[first] == stage.UNKNOWN_LINEAGE

    def test_the_method_is_st_and_nothing_else(self, config):
        assert config.lineage.method == "st"


class TestTheCrosswalkIsWrittenWithTheContractColumns:
    """`docs/data_contract.md` declares five columns for
    `phylogeny/tree_metadata.tsv`. The producer wrote two and let
    `extrasaction="ignore"` drop the rest, which is how a producer with no
    `lineage_label` in it looked like a producer that simply had nothing to
    say."""

    def _produce(self, monkeypatch, tmp_path, config, ids, **overrides):
        """Run the REAL branch with the adapter stubbed, so the file's COLUMNS
        can be asserted without a tree a fixture cannot support under +ASC."""
        import copy
        from pathlib import Path

        from papipeline.adapters import iqtree as adapter
        from papipeline.models import RunMode

        def fake_build_tree(alignment, out_dir, **_kwargs):
            out = Path(out_dir)
            out.mkdir(parents=True, exist_ok=True)
            (out / "iqtree.treefile").write_text(
                "(" + ",".join(f"{i}:0.0" for i in ids) + ");\n", encoding="utf-8"
            )
            return adapter.TreeResult(
                tree=out / "iqtree.treefile", report=None, log=None,
                tip_labels=list(ids), crosswalk={i: i for i in ids},
                command=["iqtree"],
            )

        monkeypatch.setattr(adapter, "build_tree", fake_build_tree)
        phylo_dir = tmp_path / "phylogeny"
        phylo_dir.mkdir()
        (phylo_dir / "core_snp_alignment.fasta").write_text(">x\nACGT\n")
        raw = copy.deepcopy(dict(config.raw))
        raw["phylogeny"] = dict(raw.get("phylogeny") or {})
        raw["phylogeny"].update(overrides)
        prepared = copy.copy(config)
        object.__setattr__(prepared, "raw", raw)
        stage.run(
            prepared,
            SampleManifest(samples=[Sample(sample_id=i) for i in ids]),
            RunMode.REAL,
            phylo_dir,
        )
        return phylo_dir / "tree_metadata.tsv"

    def test_all_five_contract_columns_are_present(
        self, monkeypatch, tmp_path, config, stage3_table
    ):
        from papipeline.io.tsv import read_tsv

        stage3_table({"S1": "235", "S2": "164"})
        path = self._produce(
            monkeypatch, tmp_path, config, ["S1", "S2"],
            require_exact_sample_match=False,
        )
        lines = path.read_text(encoding="utf-8").splitlines()
        header = next(l for l in lines if not l.startswith("#")).split("\t")
        assert header == list(stage.TREE_METADATA_COLUMNS), (
            f"the contract declares {stage.TREE_METADATA_COLUMNS}; the file has "
            f"{tuple(header)}"
        )
        rows = {r["sample_id"]: r for r in read_tsv(path)}
        assert rows["S1"]["lineage_label"] == "235"
        assert rows["S2"]["lineage_label"] == "164"

    def test_st_and_lineage_label_agree_because_r4_says_they_do(
        self, monkeypatch, tmp_path, config, stage3_table
    ):
        """Under R4 a lineage IS the sequence type, so the two columns are one
        fact recorded twice. Keeping both means tree-cut, when it exists, is a
        visible change in `st` rather than a silent reinterpretation."""
        from papipeline.io.tsv import read_tsv

        stage3_table({"S1": "235", "S2": "164"})
        path = self._produce(
            monkeypatch, tmp_path, config, ["S1", "S2"],
            require_exact_sample_match=False,
        )
        rows = {r["sample_id"]: r for r in read_tsv(path)}
        for row in rows.values():
            assert row["lineage_label"] == row["st"]

    def test_the_source_column_survives(self, monkeypatch, tmp_path, config,
                                        stage3_table):
        """It was dropped silently alongside lineage_label, which is the part of
        the original defect nobody had noticed."""
        from papipeline.io.tsv import read_tsv

        stage3_table({"S1": "235"})
        path = self._produce(
            monkeypatch, tmp_path, config, ["S1"],
            require_exact_sample_match=False,
        )
        assert read_tsv(path)[0]["source"] == "iqtree"

    def test_the_method_is_recorded_in_the_files_own_provenance(
        self, monkeypatch, tmp_path, config, stage3_table
    ):
        """`lineage.method` has to travel with the data. `read_tsv` skips `#`
        lines, so this costs the reader nothing and cannot drift from the file
        it describes."""
        from papipeline.models import RunMode

        stage3_table({"S1": "235"})
        path = self._produce(
            monkeypatch, tmp_path, config, ["S1"],
            require_exact_sample_match=False,
        )
        comments = [
            l for l in path.read_text(encoding="utf-8").splitlines()
            if l.startswith("#")
        ]
        assert comments, "the provenance header is missing entirely"
        head = "\n".join(comments)
        assert "lineage_method=st" in head
        assert "lineage_source=" in head
        assert str(stage.mlst_results_path(config, RunMode.REAL)) in head
        # And it is skipped by the reader, so the file is still a plain TSV.
        from papipeline.io.tsv import read_tsv
        assert read_tsv(path)[0]["sample_id"] == "S1"

    def test_a_second_lineage_method_is_refused_not_implemented(
        self, monkeypatch, tmp_path, config, stage3_table
    ):
        """`LineageConfig` already refuses an unaccepted method at load time;
        this asserts the stage would not quietly label things if it did not."""
        import copy

        from papipeline.models import RunMode

        stage3_table({"S1": "235"})
        prepared = copy.copy(config)
        object.__setattr__(prepared, "lineage",
                           type(config.lineage)(method="tree_cut"))
        with pytest.raises(Exception) as excinfo:
            stage.load_lineage_labels(
                prepared, RunMode.REAL,
                SampleManifest(samples=[Sample(sample_id="S1")]),
            )
        assert "tree_cut" in str(excinfo.value)
