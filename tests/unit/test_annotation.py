"""Tests for the annotation parser (stage 2).

The Bakta parsers are exercised against literal GFF/TSV text, so these tests
run with no Bakta install and no genome data.
"""

from __future__ import annotations

import pytest

import scripts.run_annotation_sweep as sweep

from papipeline.errors import DataContractError
from papipeline.stages import annotation as stage

GFF_HEADER = "##gff-version 3\n"

GFF_BODY = (
    "contig_1\tbakta\tCDS\t1\t300\t.\t+\t0\tID=1_1;locus_tag=oprD;product=porin;"
    "gene=oprD\n"
    "contig_1\tbakta\tCDS\t400\t800\t.\t-\t0\tID=1_2;locus_tag=mexR;"
    "product=transcriptional regulator\n"
    "contig_1\tbakta\ttRNA\t900\t1000\t.\t+\t.\tID=1_3;locus_tag=trna1\n"
    "contig_1\tbakta\tCDS\t1100\t1200\t.\t+\t0\tID=1_4;product=no name here\n"
)



def _without_annotation_db(config):
    """A copy of `config` with `annotation.bakta_db` removed.

    The un-configured state has to be built deliberately. It used to arrive for
    free, because `annotation.bakta_db` was simply absent from `science.yaml` -
    which meant the test below stopped testing its own subject the moment a
    database was pinned, and then invoked Bakta for real. Whether a key is
    *missing* from the repository's config is not something a test should depend
    on; that is a statement about today's repository, not about the code.
    """
    import copy

    unconfigured = copy.copy(config)
    raw = copy.deepcopy(dict(config.raw))
    raw.pop("annotation", None)
    object.__setattr__(unconfigured, "raw", raw)
    return unconfigured


class TestParseBaktaGff:
    def test_parses_named_cds_features(self):
        features = stage.parse_bakta_gff(GFF_HEADER + GFF_BODY)
        names = [f["gene_name"] for f in features]
        assert "oprD" in names
        assert "mexR" in names

    def test_skips_comments(self):
        features = stage.parse_bakta_gff(GFF_HEADER + GFF_BODY)
        assert all(not f["gene_name"].startswith("##") for f in features)

    def test_skips_unnamed_features(self):
        """A feature with no name cannot be matched across samples."""
        features = stage.parse_bakta_gff(GFF_HEADER + GFF_BODY)
        assert "no name here" not in [f["gene_name"] for f in features]

    def test_prefers_gene_over_locus_tag(self):
        features = stage.parse_bakta_gff(GFF_HEADER + GFF_BODY)
        oprd = [f for f in features if f["gene_name"] == "oprD"][0]
        assert oprd["gene_id"] == "oprD"
        assert oprd["product"] == "porin"

    def test_falls_back_to_locus_tag(self):
        gff = "contig_1\tbakta\tCDS\t1\t10\t.\t+\t0\tlocus_tag=mexZ;product=p\n"
        features = stage.parse_bakta_gff(gff)
        assert features[0]["gene_name"] == "mexZ"

    def test_captures_coordinates_and_strand(self):
        features = stage.parse_bakta_gff(GFF_HEADER + GFF_BODY)
        oprd = [f for f in features if f["gene_name"] == "oprD"][0]
        assert oprd["start"] == "1"
        assert oprd["end"] == "300"
        assert oprd["strand"] == "+"
        assert oprd["seqid"] == "contig_1"

    def test_trna_is_kept(self):
        features = stage.parse_bakta_gff(GFF_HEADER + GFF_BODY)
        assert any(f["gene_type"] == "tRNA" for f in features)

    def test_gene_feature_is_ignored(self):
        gff = "contig_1\tbakta\tgene\t1\t900\t.\t+\t.\tID=g1;locus_tag=oprD\n"
        assert stage.parse_bakta_gff(gff) == []

    def test_url_escapes_decoded(self):
        gff = "c\tb\tCDS\t1\t9\t.\t+\t0\tlocus_tag=a;product=x%2Cy%3Bz\n"
        features = stage.parse_bakta_gff(gff)
        assert features[0]["product"] == "x,y;z"

    def test_quoted_attribute_value(self):
        gff = 'c\tb\tCDS\t1\t9\t.\t+\t0\tlocus_tag=a;product="a, long product"\n'
        features = stage.parse_bakta_gff(gff)
        assert features[0]["product"] == "a, long product"

    def test_short_line_skipped(self):
        assert stage.parse_bakta_gff("c\tb\tCDS\t1\t9\n") == []

    def test_empty_input(self):
        assert stage.parse_bakta_gff("") == []


class TestParseBaktaTsv:
    TSV = (
        "Sequence Id,Type,Strand,Locus Tag,Start,Stop,Length,GC,Gene,Product,Contig\n"
        "contig_1,CDS,+,oprD,1,300,300,52.0,oprD,outer membrane porin,contig_1\n"
        "contig_1,CDS,-,mexR,400,800,401,48.0,,transcriptional regulator,contig_1\n"
    )

    def test_parses_rows(self):
        features = stage.parse_bakta_tsv(self.TSV)
        assert len(features) == 2

    def test_gene_column_preferred(self):
        features = stage.parse_bakta_tsv(self.TSV)
        assert features[0]["gene_name"] == "oprD"

    def test_falls_back_to_locus_tag(self):
        features = stage.parse_bakta_tsv(self.TSV)
        assert features[1]["gene_name"] == "mexR"

    def test_missing_required_column_raises(self):
        with pytest.raises(DataContractError, match="missing required columns"):
            stage.parse_bakta_tsv("Sequence Id,Type\nc1,CDS\n")

    def test_empty_tsv_raises(self):
        """An empty annotation file is a contract error, not zero genes.

        Returning [] here would silently hide a truncated or failed run.
        """
        with pytest.raises(DataContractError, match="no header"):
            stage.parse_bakta_tsv("")

    def test_comments_only_raises(self):
        with pytest.raises(DataContractError, match="no header"):
            stage.parse_bakta_tsv("# Annotated with Bakta\n# Software: v1.12.1\n")

    def test_header_only_tsv_yields_no_features(self):
        assert stage.parse_bakta_tsv("Sequence Id,Type,Strand,Locus Tag,Start,Stop\n") == []


REAL_BAKTA_TSV = (
    "# Annotated with Bakta\n"
    "# Software: v1.12.1\n"
    "# Database: v6.0, light\n"
    "# DOI: 10.1099/mgen.0.000685\n"
    "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n"
    "contig_1\tcds\t105\t1460\t-\toprD\toprD\touter membrane protein D\tUniRef:UniRef50_X\n"
    "contig_1\tcds\t2000\t2400\t+\tnalC\tnalC\tFused N-acetylmuramoyl-L-alanine\tamphyl_isoleucine\tUniRef:UniRef50_Y\n"
)


class TestParseRealBaktaTsv:
    """Real Bakta 1.12 output is TAB-delimited with a '#'-prefixed header.

    Both facts were discovered only by running against real annotation
    output; the original parser assumed comma-delimited rows with a plain
    header line and silently produced nothing useful.
    """

    def test_parses_real_format(self):
        features = stage.parse_bakta_tsv(REAL_BAKTA_TSV)
        assert len(features) == 2

    def test_header_is_read_from_the_comment_line(self):
        features = stage.parse_bakta_tsv(REAL_BAKTA_TSV)
        assert features[0]["gene_name"] == "oprD"
        assert features[0]["start"] == "105"
        assert features[0]["strand"] == "-"

    def test_human_preamble_is_not_mistaken_for_the_header(self):
        """'# Annotated with Bakta' must not become the column header."""
        features = stage.parse_bakta_tsv(REAL_BAKTA_TSV)
        assert all("Annotated" not in (f["gene_name"] or "") for f in features)

    def test_second_gene_parsed(self):
        features = stage.parse_bakta_tsv(REAL_BAKTA_TSV)
        assert features[1]["gene_name"] == "nalC"

    def test_standardise_integers(self):
        records = stage.standardise(stage.parse_bakta_tsv(REAL_BAKTA_TSV), "S1")
        assert records[0].start == 105
        assert records[0].end == 1460

    def test_gff_and_tsv_agree_on_gene_count(self):
        gff = stage.parse_bakta_gff(GFF_HEADER + GFF_BODY)
        tsv = stage.parse_bakta_tsv(REAL_BAKTA_TSV)
        assert len(gff) == 4 and len(tsv) == 2


class TestStandardise:
    def test_converts_to_records(self):
        features = stage.parse_bakta_gff(GFF_HEADER + GFF_BODY)
        records = stage.standardise(features, "TEST_PA_001")
        assert all(r.sample_id == "TEST_PA_001" for r in records)
        assert all(isinstance(r.start, int) for r in records)

    def test_records_carry_source(self):
        features = stage.parse_bakta_gff(GFF_HEADER + GFF_BODY)
        records = stage.standardise(features, "S1", source="bakta")
        assert all(r.annotation_source == "bakta" for r in records)

    def test_non_integer_coordinate_becomes_none(self):
        features = [{"seqid": "c", "gene_id": "g", "gene_name": "x", "start": "?", "end": "y"}]
        record = stage.standardise(features, "S1")[0]
        assert record.start is None
        assert record.end is None

    def test_null_gene_name_allowed_on_record(self):
        features = [{"seqid": "c", "gene_id": "g", "gene_name": None, "gene_type": "CDS"}]
        record = stage.standardise(features, "S1")[0]
        assert record.gene_name is None


class TestStandardisedRoundTrip:
    def test_write_then_read(self, tmp_path):
        from papipeline.io.tsv import write_tsv

        features = stage.parse_bakta_gff(GFF_HEADER + GFF_BODY)
        records = stage.standardise(features, "TEST_PA_001")
        path = tmp_path / "ann.tsv"
        write_tsv(path, [r.to_row() for r in records], list(stage.ANNOTATION_COLUMNS))
        reloaded = stage.read_standardised(path)
        assert len(reloaded) == len(records)
        assert reloaded[0].gene_name == records[0].gene_name
        assert reloaded[0].start == records[0].start


class TestLoadFromIntermediate:
    def test_loads_all_synthetic_samples(self, intermediate_root, manifest):
        annotations = stage.load_from_intermediate(intermediate_root, manifest)
        assert set(annotations) == set(manifest.sample_ids)
        assert all(len(v) > 0 for v in annotations.values())

    def test_sample_without_annotation_maps_to_empty_list(
        self, intermediate_root, manifest, tmp_path
    ):
        annotations = stage.load_from_intermediate(tmp_path, manifest)
        assert annotations[manifest.sample_ids[0]] == []

    def test_gene_names_helper(self, intermediate_root, manifest):
        annotations = stage.load_from_intermediate(intermediate_root, manifest)
        names = stage.gene_names(annotations["TEST_PA_001"])
        assert "oprD" in names or "mexR" in names
        assert names == sorted(set(names))

    def test_oprd_absent_from_at_least_one_sample(self, intermediate_root, manifest):
        """The fixture must exercise OprD absence."""
        annotations = stage.load_from_intermediate(intermediate_root, manifest)
        absent = [
            sid
            for sid, records in annotations.items()
            if "oprD" not in {r.gene_name for r in records}
        ]
        assert absent, "no synthetic sample lacks oprD"


class TestRunMode:
    def test_real_mode_demands_a_pinned_database(self, config, manifest, intermediate_root):
        """REAL mode is implemented, so it now fails on *configuration*.

        This assertion used to be ``NotImplementedError``, which recorded
        that the adapter did not exist. It does now, and the honest failure
        for a run with no pinned Bakta database is a configuration error -
        Bakta must never be pointed at a database it downloaded itself.
        """
        from papipeline.errors import StageError
        from papipeline.models import RunMode

        unconfigured = _without_annotation_db(config)
        with pytest.raises(StageError) as caught:
            stage.run(unconfigured, manifest, RunMode.REAL, intermediate_root)
        assert "pinned Bakta database" in str(caught.value)
        assert caught.value.context.get("stage") == "annotation"

    def test_real_mode_refuses_a_database_that_is_not_the_pinned_one(
        self, config, manifest, intermediate_root, tmp_path
    ):
        """The refusal the stage could not previously make.

        Pointing `annotation.bakta_db` at a *wrong* database used to be
        indistinguishable from pointing it at the right one: the stage resolved
        the path and ran. `docs/design/07-bakta-optimization.md` records tool /
        database disagreement as the lesson to encode as a preflight rather than
        rediscover mid-run, and a `StageError` about a missing setting cannot
        express it at all.
        """
        import copy
        import json

        from papipeline.errors import PipelineError
        from papipeline.models import RunMode

        impostor = tmp_path / "db-light"
        impostor.mkdir()
        (impostor / "version.json").write_text(
            json.dumps(
                {"date": "2025-09-01", "major": 6, "minor": 0, "type": "light"}
            ),
            encoding="utf-8",
        )

        redirected = copy.copy(config)
        raw = copy.deepcopy(dict(config.raw))
        raw["annotation"] = dict(raw.get("annotation") or {})
        raw["annotation"]["bakta_db"] = str(impostor)
        object.__setattr__(redirected, "raw", raw)

        with pytest.raises(PipelineError) as caught:
            stage.run(redirected, manifest, RunMode.REAL, intermediate_root)
        message = str(caught.value)
        assert "2025-09-01" in message and "2025-02-24" in message, (
            "the refusal must name both the database on disk and the pinned "
            "release, or the reader cannot tell which of the two is wrong"
        )

    def test_real_mode_is_still_gated_by_configuration(self, config, manifest):
        """Implementing the adapter did not unlock REAL mode.

        ``allow_real_mode: false`` is what keeps the analysis phase shut, and
        that gate is enforced in ``run_pipeline`` before any stage runs.
        """
        from papipeline.errors import ModeNotAllowedError
        from papipeline.models import RunMode
        from papipeline.run import resolve_mode

        assert config.runtime.get("allow_real_mode") is False
        with pytest.raises(ModeNotAllowedError):
            resolve_mode("REAL", config)
        assert RunMode.REAL.value == "REAL"


class TestBaktaOutputSelection:
    """Bakta writes three TSVs per genome and they are not interchangeable.

    ``bakta_outputs`` originally returned the ``.inference.tsv`` file. That
    table has no ``Gene`` or ``Product`` column, so every gene name was lost
    and gene-level analysis silently fell back to locus tags - all 11 stage-6
    loci appeared absent in all 65 genomes. The annotation table is the only
    one carrying gene names.
    """

    @staticmethod
    def _dir(tmp_path, name):
        d = tmp_path / name
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{name}.gff3").write_text("##gff-version 3\n")
        (d / f"{name}.tsv").write_text(
            "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n"
            "contig_1\tcds\t1\t100\t+\tLOCUS_1\toprD\tOprD family porin\t-\n"
        )
        (d / f"{name}.inference.tsv").write_text(
            "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tScore\tEvalue\t"
            "Query Cov\tSubject Cov\tId\tAccession\n"
        )
        (d / f"{name}.hypotheticals.tsv").write_text(
            "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tScore\tEvalue\t"
            "Query Cov\tSubject Cov\tId\tAccession\n"
        )
        return d

    def test_returns_the_annotation_table_not_the_inference_table(self, tmp_path):
        d = self._dir(tmp_path, "S1")
        _, table = sweep.bakta_outputs(d)
        assert table.name == "S1.tsv"

    def test_inference_table_is_available_separately(self, tmp_path):
        d = self._dir(tmp_path, "S1")
        assert sweep.bakta_inference_tsv(d).name == "S1.inference.tsv"

    def test_hypotheticals_table_is_never_returned(self, tmp_path):
        d = self._dir(tmp_path, "S1")
        _, table = sweep.bakta_outputs(d)
        assert "hypotheticals" not in table.name

    def test_gene_names_survive(self, tmp_path):
        """The end-to-end consequence: oprD must be recoverable."""
        d = self._dir(tmp_path, "S1")
        _, table = sweep.bakta_outputs(d)
        genes = {r["gene_name"] for r in stage.parse_bakta_tsv(table.read_text())}
        assert "oprD" in genes

    def test_completion_requires_the_annotation_table(self, tmp_path):
        d = self._dir(tmp_path, "S1")
        assert sweep.is_complete(d) is True
        (d / "S1.tsv").unlink()
        assert sweep.is_complete(d) is False

    def test_completion_is_false_without_a_gff(self, tmp_path):
        d = self._dir(tmp_path, "S1")
        (d / "S1.gff3").unlink()
        assert sweep.is_complete(d) is False

    def test_inference_table_alone_is_not_complete(self, tmp_path):
        """The exact false-positive this fix removes."""
        d = self._dir(tmp_path, "S1")
        (d / "S1.tsv").unlink()
        remaining = [p.name for p in d.glob("*.tsv")]
        assert any(n.endswith(".inference.tsv") for n in remaining)
        assert sweep.is_complete(d) is False


class TestOutputPruning:
    """Pruning must reclaim disk without ever invalidating a sample.

    At 200 genomes the difference is ~17 GB of Bakta output versus ~2 GB,
    and a prune that removed the GFF or the annotation TSV would silently
    force a full re-annotation.
    """

    @staticmethod
    def _sample(tmp_path, name="S1"):
        d = tmp_path / name
        d.mkdir(parents=True, exist_ok=True)
        stem = f"{name}_genomic"
        (d / f"{stem}.gff3").write_text("##gff-version 3\n")
        (d / f"{stem}.tsv").write_text("#Sequence Id\tGene\ncontig_1\toprD\n")
        (d / f"{stem}.inference.tsv").write_text("#Sequence Id\tScore\n")
        (d / f"{stem}.hypotheticals.tsv").write_text("#Sequence Id\tScore\n")
        (d / f"{stem}.log").write_text("bakta log\n")
        for ext in ("json", "embl", "gbff", "svg", "fna", "ffn", "faa", "png"):
            (d / f"{stem}.{ext}").write_text("x" * 4096)
        return d

    def test_keeps_the_files_the_pipeline_reads(self, tmp_path):
        d = self._sample(tmp_path)
        sweep.prune_outputs(d)
        kept = {p.name for p in d.iterdir()}
        assert any(n.endswith(".gff3") for n in kept)
        assert any(n.endswith(".tsv") and not n.endswith(".inference.tsv") for n in kept)
        assert any(n.endswith(".inference.tsv") for n in kept)
        assert any(n.endswith(".log") for n in kept)

    def test_removes_unused_formats(self, tmp_path):
        d = self._sample(tmp_path)
        sweep.prune_outputs(d)
        for ext in ("json", "embl", "gbff", "svg", "fna", "ffn", "faa", "png"):
            assert not list(d.glob(f"*.{ext}")), f"{ext} should have been pruned"

    def test_sample_stays_complete(self, tmp_path):
        d = self._sample(tmp_path)
        assert sweep.is_complete(d) is True
        sweep.prune_outputs(d)
        assert sweep.is_complete(d) is True

    def test_completion_still_uses_the_annotation_table(self, tmp_path):
        d = self._sample(tmp_path)
        sweep.prune_outputs(d)
        _, table = sweep.bakta_outputs(d)
        assert "inference" not in table.name
        assert "hypotheticals" not in table.name

    def test_is_idempotent(self, tmp_path):
        d = self._sample(tmp_path)
        sweep.prune_outputs(d)
        first = sorted(p.name for p in d.iterdir())
        removed, freed = sweep.prune_outputs(d)
        assert removed == 0 and freed == 0
        assert sorted(p.name for p in d.iterdir()) == first

    def test_reports_freed_bytes(self, tmp_path):
        d = self._sample(tmp_path)
        removed, freed = sweep.prune_outputs(d)
        assert removed == 8
        assert freed == 8 * 4096

    def test_missing_directory_is_not_an_error(self, tmp_path):
        assert sweep.prune_outputs(tmp_path / "nope") == (0, 0)

    def test_prune_all_walks_sample_directories(self, tmp_path):
        self._sample(tmp_path, "S1")
        self._sample(tmp_path, "S2")
        result = sweep.prune_all(tmp_path)
        assert result["files_removed"] == 16
        assert result["bytes_freed"] == 16 * 4096


class TestRunOneNeverRaises:
    """One bad genome must not abandon the whole queue.

    With 200 genomes an unhandled exception would cost every genome not yet
    started, so run_one converts every failure into a state record.
    """

    def test_missing_genome_is_a_record_not_an_exception(self, tmp_path):
        record = sweep.run_one(
            "PILOT100_001", tmp_path / "absent.fna", tmp_path / "out",
            tmp_path / "db", 1, None, 60,
        )
        assert record["status"] == "error"
        assert "not found" in record["error"]

    def test_empty_genome_is_a_record(self, tmp_path):
        genome = tmp_path / "empty.fna"
        genome.write_text("")
        record = sweep.run_one(
            "PILOT100_001", genome, tmp_path / "out", tmp_path / "db", 1, None, 60,
        )
        assert record["status"] == "error"
        assert "empty" in record["error"]

    def test_missing_bakta_binary_is_a_record(self, tmp_path, monkeypatch):
        genome = tmp_path / "g.fna"
        genome.write_text(">c\nACGT\n")
        def boom():
            raise RuntimeError("bakta is not on PATH")
        monkeypatch.setattr(sweep, "bakta_executable", boom)
        record = sweep.run_one(
            "PILOT100_001", genome, tmp_path / "out", tmp_path / "db", 1, None, 60,
        )
        assert record["status"] == "error"
        assert "unavailable" in record["error"]

    def test_unwritable_output_directory_is_a_record(self, tmp_path):
        genome = tmp_path / "g.fna"
        genome.write_text(">c\nACGT\n")
        blocker = tmp_path / "blocked"
        blocker.write_text("i am a file, not a directory")
        record = sweep.run_one(
            "PILOT100_001", genome, blocker / "out", tmp_path / "db", 1, None, 60,
        )
        assert record["status"] == "error"
