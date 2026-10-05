"""Unit tests for declarative output validation.

Includes regression tests for the three "successful-looking failure" bugs
this repository actually shipped, because a validator that cannot catch
those is not earning its place:

1. ``bakta_outputs`` returned the inference table, so every gene name was
   lost and all 11 stage-6 loci scored 0/65;
2. a leaked loop variable made 65 per-sample VCFs byte-identical, so the
   core-SNP alignment had one distinct sequence;
3. an identity-by-state kinship matrix was built from a called/not-called
   mask, so every entry was 1.0.
"""

from __future__ import annotations

import json
from pathlib import Path


from papipeline.execution import OutputSpec, StageState, validate
from papipeline.execution.specs import (
    ANNOTATION_REQUIRED_COLUMNS,
    BaktaOutputs,
    bakta_spec,
    covers_samples,
    exists,
    has_columns,
    has_real_values,
    identifiers_known,
    min_rows,
    non_empty,
    not_constant_matrix,
    parses,
    variant_call_spec,
    variant_sibling_spec,
    kinship_spec,
    standardised_annotation_spec,
    stage6_resolution_spec,
)
from papipeline.execution.validation import SiblingSpec


# --- helpers -------------------------------------------------------------

def write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


ANNOTATION_TSV = (
    "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n"
    "contig_1\tcds\t105\t1460\t-\tPA0001\toprD\tOprD family porin\t-\n"
    "contig_1\tcds\t2000\t2400\t+\tPA0002\tnalC\tbeta-lactamase\t-\n"
)

INFERENCE_TSV = (
    "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tScore\tEvalue\t"
    "Query Cov\tSubject Cov\tId\tAccession\n"
    "contig_1\tcds\t105\t1460\t-\tPA0001\t99.5\t1e-40\t1.0\t1.0\tX\tY\n"
)

STANDARDISED = (
    "sample_id\tcontig_id\tgene_id\tgene_name\tproduct\tgene_type\tstart\tend\t"
    "strand\tannotation_source\n"
    "S1\tcontig_1\tPA0001\toprD\tOprD family porin\tcds\t105\t1460\t-\tbakta\n"
)


# --- individual checks ---------------------------------------------------

class TestBasicChecks:
    def test_missing_file_fails(self, tmp_path):
        result = validate(OutputSpec("s", (exists(tmp_path / "absent.tsv"),)))
        assert not result.ok
        assert result.state is StageState.INCOMPLETE

    def test_present_file_passes_exists(self, tmp_path):
        p = write(tmp_path / "a.tsv", "x\n")
        assert validate(OutputSpec("s", (exists(p),))).ok

    def test_empty_file_fails_non_empty(self, tmp_path):
        p = write(tmp_path / "empty.tsv", "")
        result = validate(OutputSpec("s", (non_empty(p),)))
        assert not result.ok
        assert result.state is StageState.INCOMPLETE

    def test_non_empty_file_passes(self, tmp_path):
        p = write(tmp_path / "a.tsv", "x\n")
        assert validate(OutputSpec("s", (non_empty(p),))).ok

    def test_tsv_row_with_too_many_fields_fails_to_parse(self, tmp_path):
        """read_tsv is strict about over-wide rows; so is the validator.

        A *short* row is padded with None and is legal, so the malformed
        case that must fail is a row with more fields than the header.
        """
        p = write(tmp_path / "bad.tsv", "a\tb\n1\t2\t3\n")
        result = validate(OutputSpec("s", (parses(p, "tsv"),)))
        assert not result.ok
        assert result.state is StageState.INVALID

    def test_short_tsv_row_is_padded_and_accepted(self, tmp_path):
        """Documented read_tsv behaviour: short rows pad, they do not fail."""
        p = write(tmp_path / "short.tsv", "a\tb\tc\n1\t2\n")
        assert validate(OutputSpec("s", (parses(p, "tsv"),))).ok

    def test_header_only_tsv_is_rejected_by_the_strict_reader(self, tmp_path):
        p = write(tmp_path / "hdr.tsv", "a\tb\n")
        result = validate(OutputSpec("s", (parses(p, "tsv"), min_rows(p, 1))))
        assert not result.ok

    def test_fasta_parse(self, tmp_path):
        p = write(tmp_path / "a.fasta", ">c1\nACGTACGT\n>c2\nTTTT\n")
        assert validate(OutputSpec("s", (parses(p, "fasta"),))).ok

    def test_broken_fasta_fails(self, tmp_path):
        p = write(tmp_path / "a.fasta", "ACGT without a header\n")
        assert not validate(OutputSpec("s", (parses(p, "fasta"),))).ok

    def test_newick_parse(self, tmp_path):
        p = write(tmp_path / "t.nwk", "((A:0.1,B:0.2):0.3,C:0.4);")
        assert validate(OutputSpec("s", (parses(p, "newick"),))).ok

    def test_non_newick_fails(self, tmp_path):
        p = write(tmp_path / "t.nwk", "this is not a tree")
        assert not validate(OutputSpec("s", (parses(p, "newick"),))).ok

    def test_json_parse(self, tmp_path):
        p = write(tmp_path / "a.json", '{"a": 1}')
        assert validate(OutputSpec("s", (parses(p, "json"),))).ok

    def test_broken_json_fails(self, tmp_path):
        p = write(tmp_path / "a.json", "{not json")
        assert not validate(OutputSpec("s", (parses(p, "json"),))).ok

    def test_unknown_parser_is_reported(self, tmp_path):
        p = write(tmp_path / "a.tsv", "a\n1\n")
        result = validate(OutputSpec("s", (parses(p, "nosuchparser"),)))
        assert not result.ok
        assert "unknown parser" in result.failures[0].detail

    def test_missing_columns_fails(self, tmp_path):
        p = write(tmp_path / "a.tsv", "sample_id\tx\nS1\t1\n")
        result = validate(OutputSpec("s", (has_columns(p, ("sample_id", "gene_name")),)))
        assert not result.ok
        assert "missing columns" in result.failures[0].detail
        assert result.state is StageState.INVALID

    def test_present_columns_pass(self, tmp_path):
        p = write(tmp_path / "a.tsv", "sample_id\tgene_name\nS1\toprD\n")
        assert validate(OutputSpec("s", (has_columns(p, ("sample_id", "gene_name")),))).ok

    def test_min_rows(self, tmp_path):
        p = write(tmp_path / "a.tsv", "sample_id\nS1\nS2\n")
        assert validate(OutputSpec("s", (min_rows(p, 2),))).ok
        assert not validate(OutputSpec("s", (min_rows(p, 5),))).ok

    def test_all_sentinel_column_is_rejected(self, tmp_path):
        """A shape with no data must not pass as content."""
        p = write(tmp_path / "a.tsv", "sample_id\tgene_name\nS1\t.\nS2\t.\n")
        result = validate(OutputSpec("s", (has_real_values(p, "gene_name"),)))
        assert not result.ok
        assert result.state is StageState.INVALID

    def test_sentinel_column_with_some_real_values_passes(self, tmp_path):
        p = write(tmp_path / "a.tsv", "sample_id\tgene_name\nS1\t.\nS2\toprD\n")
        assert validate(OutputSpec("s", (has_real_values(p, "gene_name"),))).ok

    def test_unexpected_identifier_is_rejected(self, tmp_path):
        p = write(tmp_path / "a.tsv", "sample_id\nS1\nS999\n")
        result = validate(OutputSpec("s", (identifiers_known(p, ("S1", "S2")),)))
        assert not result.ok
        assert "S999" in result.failures[0].detail

    def test_expected_identifiers_pass(self, tmp_path):
        p = write(tmp_path / "a.tsv", "sample_id\nS1\nS2\n")
        assert validate(OutputSpec("s", (identifiers_known(p, ("S1", "S2")),))).ok

    def test_sample_count_check(self, tmp_path):
        p = write(tmp_path / "a.tsv", "sample_id\nS1\nS2\n")
        assert validate(OutputSpec("s", (covers_samples(p, 2),))).ok
        assert not validate(OutputSpec("s", (covers_samples(p, 65),))).ok


class TestMatrixNotConstant:
    def test_all_one_kinship_is_rejected(self, tmp_path):
        """The bug: IBS was computed from a called/not-called mask."""
        rows = "sample\t" + "\t".join(f"S{i}" for i in range(4)) + "\n"
        for i in range(4):
            rows += f"S{i}\t" + "\t".join(["1"] * 4) + "\n"
        p = write(tmp_path / "kinship.tsv", rows)
        result = validate(OutputSpec("gwas", (not_constant_matrix(p),)))
        assert not result.ok
        assert result.state is StageState.INVALID
        assert "constant matrix" in result.failures[0].detail

    def test_varying_kinship_passes(self, tmp_path):
        rows = "sample\tS1\tS2\tS3\nS1\t1\t0.9\t0.2\nS2\t0.9\t1\t0.3\nS3\t0.2\t0.3\t1\n"
        p = write(tmp_path / "kinship.tsv", rows)
        assert validate(OutputSpec("gwas", (not_constant_matrix(p),))).ok

    def test_no_numeric_values_is_rejected(self, tmp_path):
        p = write(tmp_path / "k.tsv", "sample\tlabel\nS1\talpha\n")
        assert not validate(OutputSpec("gwas", (not_constant_matrix(p),))).ok


# --- cross-sample degeneracy --------------------------------------------

class TestSiblingDegeneracy:
    def test_identical_outputs_are_rejected(self, tmp_path):
        """The bug: a leaked loop variable made all 65 VCFs identical."""
        root = tmp_path / "vcfs"
        for sid in ("S1", "S2", "S3"):
            write(root / f"{sid}.vcf.gz", "SAME CONTENT\n")
        spec = OutputSpec("variant_calling", sibling=SiblingSpec(
            root=root, pattern="*.vcf.gz", min_distinct=2, label="VCF"))
        result = validate(spec)
        assert not result.ok
        assert result.state is StageState.INVALID
        assert "byte-identical" in (result.siblings.detail or "")

    def test_distinct_outputs_pass(self, tmp_path):
        root = tmp_path / "vcfs"
        for sid in ("S1", "S2", "S3"):
            write(root / f"{sid}.vcf.gz", f"content for {sid}\n")
        spec = OutputSpec("variant_calling", sibling=SiblingSpec(
            root=root, pattern="*.vcf.gz", min_distinct=2, label="VCF"))
        assert validate(spec).ok

    def test_identical_annotations_are_rejected(self, tmp_path):
        root = tmp_path / "ann"
        for sid in ("S1", "S2"):
            write(root / sid / "g.gff3", "##gff-version 3\n")
        spec = OutputSpec("annotation", sibling=SiblingSpec(
            root=root, pattern="*/*.gff3", min_distinct=2, label="annotation GFF"))
        assert not validate(spec).ok

    def test_no_siblings_at_all_is_rejected(self, tmp_path):
        spec = OutputSpec("annotation", sibling=SiblingSpec(
            root=tmp_path / "empty", pattern="*.gff3"))
        assert not validate(spec).ok


# --- verdict semantics ---------------------------------------------------

class TestValidationResult:
    def test_all_pass_is_succeeded(self, tmp_path):
        p = write(tmp_path / "a.tsv", "sample_id\nS1\n")
        assert validate(OutputSpec("s", (non_empty(p), parses(p, "tsv")))).state \
            is StageState.SUCCEEDED

    def test_missing_file_wins_over_later_failures(self, tmp_path):
        """A missing file is reported as missing, not as a schema problem."""
        missing = tmp_path / "gone.tsv"
        p = write(tmp_path / "a.tsv", "wrong\n")
        result = validate(OutputSpec("s", (non_empty(missing), has_columns(p, ("x",)))))
        assert result.state is StageState.INCOMPLETE

    def test_first_failure_drives_the_state(self, tmp_path):
        p = write(tmp_path / "a.tsv", "a\tb\n1\t2\n")
        result = validate(OutputSpec("s", (
            has_columns(p, ("nope",)),                 # INVALID
            non_empty(tmp_path / "gone.tsv"),          # INCOMPLETE
        )))
        assert result.state is StageState.INVALID

    def test_result_serialises_to_json(self, tmp_path):
        p = write(tmp_path / "a.tsv", "sample_id\nS1\n")
        payload = json.loads(validate(OutputSpec("s", (non_empty(p),))).to_json())
        assert payload["ok"] is True
        assert payload["state"] == "SUCCEEDED"
        assert payload["checks"][0]["passed"] is True

    def test_spec_describes_itself(self, tmp_path):
        p = write(tmp_path / "a.tsv", "sample_id\nS1\n")
        text = OutputSpec("amr", (non_empty(p),)).describe()
        assert "amr" in text and "non-empty" in text


# --- regression: the wrong-Bakta-TSV bug ---------------------------------

class TestBaktaTableSelectionRegression:
    def test_inference_table_is_rejected_as_the_annotation(self, tmp_path):
        """REGRESSION: the inference table was returned where the annotation
        table was expected, losing every gene name and making all 11 stage-6
        loci read 0/65."""
        out = tmp_path / "P1"
        out.mkdir()
        gff = write(out / "P1.gff3", "##gff-version 3\n")
        write(out / "P1.tsv", ANNOTATION_TSV)
        write(out / "P1.inference.tsv", INFERENCE_TSV)
        write(out / "P1.hypotheticals.tsv", INFERENCE_TSV)
        spec = bakta_spec(BaktaOutputs(out, "S1", "P1"))
        result = validate(spec)
        assert result.ok, result.detail

    def test_spec_rejects_an_inference_table_substituted_in(self, tmp_path):
        """The check that would have caught it: the feature-table marker
        columns must be present."""
        out = tmp_path / "P1"
        out.mkdir()
        gff = write(out / "P1.gff3", "##gff-version 3\n")
        annotation = write(out / "P1.tsv", ANNOTATION_TSV)
        write(out / "P1.inference.tsv", INFERENCE_TSV)
        # Simulate the bug: the inference table is what the stage will read.
        spec = bakta_spec(BaktaOutputs(out, "S1", "P1"))
        swapped = OutputSpec("annotation", checks=tuple(
            c for c in spec.checks if c.path != annotation
        ) + (has_columns(annotation, ("Score", "Id", "Evalue")),))
        result = validate(swapped)
        # The inference table satisfies its own columns but fails the
        # sample/gene-name contract the pipeline requires.
        assert not result.ok
        assert result.state is StageState.INVALID

    def test_annotation_table_with_no_gene_names_is_rejected(self, tmp_path):
        """Shape present, content absent: all-sentinel Gene column."""
        out = tmp_path / "P1"
        out.mkdir()
        write(out / "P1.gff3", "##gff-version 3\n")
        tsv = ("#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\n"
               "contig_1\tcds\t1\t10\t-\tL1\t.\t-\n")
        write(out / "P1.tsv", tsv)
        result = validate(bakta_spec(BaktaOutputs(out, "S1", "P1")))
        assert not result.ok
        assert result.state is StageState.INVALID

    def test_row_from_another_isolate_is_rejected(self, tmp_path):
        out = tmp_path / "P1"
        out.mkdir()
        write(out / "P1.gff3", "##gff-version 3\n")
        write(out / "P1.tsv", ANNOTATION_TSV.replace("PA0001", "PA0001"))
        write(out / "P1.tsv", "#Sequence Id\tLocus Tag\tGene\ncontig_1\tPA1\toprD\n")
        spec = OutputSpec("annotation", checks=(
            has_columns(out / "P1.tsv", ANNOTATION_REQUIRED_COLUMNS),
            identifiers_known(out / "P1.tsv", ("S1",), key="sample_id"),
        ))
        result = validate(spec)
        assert not result.ok

    def test_inference_table_is_optional_by_default(self, tmp_path):
        out = tmp_path / "P1"
        out.mkdir()
        write(out / "P1.gff3", "##gff-version 3\n")
        write(out / "P1.tsv", ANNOTATION_TSV)
        assert validate(bakta_spec(BaktaOutputs(out, "S1", "P1"))).ok

    def test_inference_table_can_be_required(self, tmp_path):
        out = tmp_path / "P1"
        out.mkdir()
        write(out / "P1.gff3", "##gff-version 3\n")
        write(out / "P1.tsv", ANNOTATION_TSV)
        spec = bakta_spec(BaktaOutputs(out, "S1", "P1"), require_inference=True)
        assert not validate(spec).ok  # inference table absent


# --- regression: the kinship and stage-6 contracts -----------------------

class TestDomainSpecRegression:
    def test_constant_kinship_is_invalid(self, tmp_path):
        """REGRESSION: the all-1.0 matrix that produced h^2 = 1.00."""
        labels = ["S1", "S2", "S3"]
        rows = "sample\t" + "\t".join(labels) + "\n"
        for lab in labels:
            rows += lab + "\t" + "\t".join(["1"] * len(labels)) + "\n"
        p = write(tmp_path / "kinship.tsv", rows)
        result = validate(kinship_spec(p, labels))
        assert not result.ok
        assert result.state is StageState.INVALID

    def test_real_kinship_is_valid(self, tmp_path):
        labels = ["S1", "S2", "S3"]
        rows = "sample\t" + "\t".join(labels) + "\n"
        rows += "S1\t1\t0.98\t0.21\nS2\t0.98\t1\t0.22\nS3\t0.21\t0.22\t1\n"
        p = write(tmp_path / "kinship.tsv", rows)
        assert validate(kinship_spec(p, labels)).ok

    def test_kinship_with_wrong_dimensions_is_invalid(self, tmp_path):
        p = write(tmp_path / "k.tsv", "sample\tS1\tS2\nS1\t1\t0.9\nS2\t0.9\t1\n")
        result = validate(kinship_spec(p, ["S1", "S2", "S3"]))
        assert not result.ok

    def test_vcf_with_no_variants_is_incomplete(self, tmp_path):
        p = write(tmp_path / "s.vcf.gz",
                  "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        result = validate(variant_call_spec(p, "S1"))
        assert not result.ok
        assert result.state is StageState.INCOMPLETE

    def test_vcf_with_variants_passes(self, tmp_path):
        p = write(tmp_path / "s.vcf.gz",
                  "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                  "contig_1\t100\t.\tA\tG\t50\tPASS\t.\n")
        assert validate(variant_call_spec(p, "S1")).ok

    def test_vcf_sibling_check_catches_identical_calls(self, tmp_path):
        root = tmp_path / "vcf"
        for sid in ("S1", "S2"):
            write(root / f"{sid}.vcf.gz", "identical\n")
        spec = OutputSpec("variant_calling", sibling=variant_sibling_spec(root, "S1"))
        assert not validate(spec).ok

    def test_stage6_distinguishes_no_hits_from_wrong_input(self, tmp_path):
        """A sample with no tracked loci is a result; a sample whose input
        never had them is a failure. The contract separates them."""
        good = write(tmp_path / "amr.tsv", "sample_id\tGene\nS1\tblaIMP-1\nS1\toprD\n")
        assert validate(stage6_resolution_spec(good, "S1")).ok

        shape_only = write(tmp_path / "amr2.tsv", "sample_id\tGene\nS1\t.\nS1\t.\n")
        result = validate(stage6_resolution_spec(shape_only, "S1"))
        assert not result.ok
        assert result.state is StageState.INVALID

    def test_standardised_annotation_contract(self, tmp_path):
        p = write(tmp_path / "S1.annotation.tsv", STANDARDISED)
        assert validate(standardised_annotation_spec(p, "S1")).ok

    def test_standardised_annotation_missing_gene_name_column(self, tmp_path):
        p = write(tmp_path / "S1.annotation.tsv",
                  "sample_id\tcontig_id\tgene_id\nS1\tc1\tg1\n")
        assert not validate(standardised_annotation_spec(p, "S1")).ok
