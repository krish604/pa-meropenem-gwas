"""The interval lookup, against the real PAO1 GFF and a synthetic alignment.

The coordinate-frame bug this module exists to avoid is asserted here rather
than described in a docstring: a variant at a PAO1 position must resolve to the
PAO1 gene covering it, and must NOT resolve to whatever sits at the same
*offset* on some assembly contig.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.adapters import gff

GFF_HEADER = "##gff-version 3\n"


def _gene_row(tag, start, end, strand="+", name="", product=""):
    attrs = f"ID=gene-{tag};Name={name};locus_tag={tag}"
    if product:
        attrs += f";product={product}"
    return f"NC_002516.2\tBakta\tgene\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n"


def _cds_row(tag, start, end):
    return (
        f"NC_002516.2\tBakta\tCDS\t{start}\t{end}\t.\t+\t0\t"
        f"ID=cds-{tag};Parent=gene-{tag};locus_tag={tag}\n"
    )


def _write_gff(tmp_path, *rows):
    path = tmp_path / "ref.gff"
    path.write_text(GFF_HEADER + "".join(rows), encoding="utf-8")
    return path


class TestLoadingIntervals:
    def test_it_indexes_the_requested_tags(self, tmp_path):
        path = _write_gff(
            tmp_path,
            _gene_row("PA0958", 100, 400, name="oprD"),
            _gene_row("PA0424", 500, 700, name="mexR"),
        )
        found = gff.load_gene_intervals(path, ["PA0958", "PA0424"])
        assert set(found) == {"PA0958", "PA0424"}
        assert found["PA0958"].name == "oprD"
        assert found["PA0958"].start == 100

    def test_cds_rows_do_not_double_count_a_locus(self, tmp_path):
        """A gene and its CDS share a locus_tag.

        Counting both makes every locus report two features. The locus resolver
        hit exactly this and failed all eleven rows on it.
        """
        path = _write_gff(
            tmp_path,
            _gene_row("PA0958", 100, 400),
            _cds_row("PA0958", 100, 400),
        )
        found = gff.load_gene_intervals(path, ["PA0958"])
        assert found["PA0958"].start == 100
        assert len(found) == 1

    def test_a_missing_gff_is_an_error_not_an_empty_index(self, tmp_path):
        """An empty index reports every isolate as carrying no regulator
        variant - a fabricated negative on the central mechanisms."""
        with pytest.raises(gff.PipelineError) as excinfo:
            gff.load_gene_intervals(tmp_path / "absent.gff", ["PA0958"])
        assert "coordinate" in str(excinfo.value).lower()

    def test_a_tag_absent_from_the_gff_is_an_error(self, tmp_path):
        path = _write_gff(tmp_path, _gene_row("PA0958", 100, 400))
        with pytest.raises(gff.PipelineError) as excinfo:
            gff.load_gene_intervals(path, ["PA0958", "PA9999"])
        assert "PA9999" in str(excinfo.value)

    def test_no_tags_requested_is_not_an_error(self, tmp_path):
        """`mexS` has no PAO1 counterpart, so the set can legitimately be empty."""
        path = _write_gff(tmp_path, _gene_row("PA0958", 100, 400))
        assert gff.load_gene_intervals(path, []) == {}


class TestCoordinateFrame:
    def test_a_pao1_position_resolves_to_the_pao1_gene(self, tmp_path):
        path = _write_gff(
            tmp_path, _gene_row("PA0958", 1000, 2000, name="oprD")
        )
        interval = gff.load_gene_intervals(path, ["PA0958"])["PA0958"]
        hit = gff.gene_at([interval], "NC_002516.2", 1500)
        assert hit is not None and hit.locus_tag == "PA0958"

    def test_a_draft_assembly_contig_never_matches(self, tmp_path):
        """The bug this design exists to prevent.

        The cohort's own assemblies are annotated on scaffolds like
        `JFJU01000001.1`. A variant in PAO1 coordinates must not resolve against
        one, even though the offsets overlap numerically.
        """
        path = _write_gff(tmp_path, _gene_row("PA0958", 1000, 2000))
        interval = gff.load_gene_intervals(path, ["PA0958"])["PA0958"]
        assert gff.gene_at([interval], "JFJU01000001.1", 1500) is None

    def test_boundaries_are_inclusive(self, tmp_path):
        path = _write_gff(tmp_path, _gene_row("PA0958", 1000, 2000))
        interval = gff.load_gene_intervals(path, ["PA0958"])["PA0958"]
        assert gff.gene_at([interval], "NC_002516.2", 1000) is not None
        assert gff.gene_at([interval], "NC_002516.2", 2000) is not None
        assert gff.gene_at([interval], "NC_002516.2", 999) is None
        assert gff.gene_at([interval], "NC_002516.2", 2001) is None


class TestCodonReading:
    def _interval(self, strand="+"):
        return gff.GeneInterval(
            locus_tag="T", contig="C", start=100, end=109, strand=strand
        )

    def test_plus_strand(self):
        # reference bases 100..109 = "AAACCCGGGT" (1-based)
        seq = "A" * 99 + "AAACCCGGGT"
        assert gff.codon_at(seq, self._interval("+"), 101) == "AAA"
        assert gff.codon_at(seq, self._interval("+"), 102) == "AAA"
        # Codons are anchored at the gene start (position 100), so they run
        # 100-102, 103-105, 106-108 - position 103 opens the second.
        assert gff.codon_at(seq, self._interval("+"), 103) == "CCC"

    def test_minus_strand_reads_the_reverse_complement(self):
        """A stop is a stop on either strand; this is why the strand is read."""
        seq = "A" * 99 + "TTTGGGCCCA"
        # Plus-strand positions 100-102 are "TTT". Read on the minus strand
        # that codon is its reverse complement, "AAA" - which is how a stop on
        # one strand is recognised as a stop on either.
        assert gff.codon_at(seq, self._interval("-"), 102) == "AAA"
        assert gff.reverse_complement("TTT") == "AAA"

    def test_a_codon_running_off_the_end_is_none(self):
        seq = "A" * 99 + "AC"
        assert gff.codon_at(seq, self._interval("+"), 101) is None

    def test_a_position_outside_the_gene_is_none(self):
        assert gff.codon_at("A" * 200, self._interval("+"), 500) is None


class TestClassification:
    def _call(self, ref, alt, position=101):
        return gff.VariantCall(
            sample_id="S", contig="C", position=position,
            reference=ref, alternate=alt,
        )

    def test_a_substitution_without_reference_context_is_an_snv(self):
        """Understates rather than invents - the deliberate direction."""
        assert gff.classify(self._call("A", "G")) == "SNV"

    def test_a_substitution_creating_a_stop_is_a_premature_stop(self):
        interval = gff.GeneInterval(
            locus_tag="T", contig="C", start=100, end=109, strand="+"
        )
        # The gene starts at 100, so codons run 100-102 ("TAA"), 103-105, ...
        # Substituting position 102 A->G turns "TAA" into "TGA", a stop.
        seq = "A" * 99 + "TAAGGGCCCA"
        assert gff.codon_at(seq, interval, 101) == "TAA"
        assert gff.classify(
            self._call("A", "G", 102), seq, interval
        ) == "PREMATURE_STOP"

    def test_a_substitution_within_a_codon_is_still_an_snv(self):
        interval = gff.GeneInterval(
            locus_tag="T", contig="C", start=100, end=109, strand="+"
        )
        seq = "A" * 99 + "TAAGGGCCCA"
        # Position 101 is the codon's first base; T->C leaves "CAA", not a stop.
        assert gff.classify(self._call("T", "C", 101), seq, interval) == "SNV"

    def test_an_in_frame_indel_is_an_indel(self):
        assert gff.classify(self._call("A", "ATTT")) == "INDEL"

    def test_an_out_of_frame_indel_is_a_frameshift(self):
        assert gff.classify(self._call("A", "AT")) == "FRAMESHIFT"
        assert gff.classify(self._call("A", "ATTT")) == "INDEL"

    def test_a_deletion_is_classified_by_its_length_change(self):
        # Length change, not total length: a 3 bp deletion shifts the frame by
        # zero, so it is in-frame however long the alleles are.
        assert gff.classify(self._call("ATTT", "A")) == "INDEL"
        assert gff.classify(self._call("ATTTT", "A")) == "FRAMESHIFT"

    def test_an_identical_allele_is_unknown(self):
        assert gff.classify(self._call("A", "A")) == "UNKNOWN"

    def test_a_missing_allele_is_unknown(self):
        assert gff.classify(self._call("", "A")) == "UNKNOWN"


class TestAssignment:
    def test_an_overlapping_position_is_reported_for_both_genes(self, tmp_path):
        """Dropping one would under-report the study's central mechanisms."""
        a = gff.GeneInterval(locus_tag="A", contig="C", start=100, end=500,
                             strand="+", name="oprD")
        b = gff.GeneInterval(locus_tag="B", contig="C", start=400, end=900,
                             strand="+", name="ampC")
        call = gff.VariantCall("S", "C", 450, "A", "G")
        assert {name for name, _ in gff.assign(call, [a, b])} == {"oprD", "ampC"}

    def test_a_position_outside_every_gene_assigns_nothing(self):
        a = gff.GeneInterval(locus_tag="A", contig="C", start=100, end=500,
                             strand="+", name="oprD")
        call = gff.VariantCall("S", "C", 900, "A", "G")
        assert gff.assign(call, [a]) == []


class TestAgainstTheRealReference:
    """The pinned PAO1 GFF, not a fixture."""

    @pytest.fixture(scope="class")
    def config(self):
        from papipeline.config.loader import load_config

        return load_config("config/science.yaml", machine="config/machines/smoke.yaml")

    def test_every_screenable_locus_resolves(self, config):
        from papipeline.config.loader import load_regulators

        specs = load_regulators(Path("config/regulators.tsv"))
        screenable = {n: s for n, s in specs.items() if s.locus_tag}
        assert len(screenable) == 10, (
            f"expected 10 screenable loci, got {len(screenable)}: "
            f"{sorted(n for n, s in specs.items() if not s.locus_tag)} lack a tag"
        )
        gff_path = config.reference_fasta().with_suffix(".gff")
        intervals = gff.load_gene_intervals(
            gff_path, [s.locus_tag for s in screenable.values()]
        )
        assert len(intervals) == 10

    def test_oprd_resolves_to_its_verified_tag(self, config):
        from papipeline.config.loader import load_regulators

        specs = load_regulators(Path("config/regulators.tsv"))
        gff_path = config.reference_fasta().with_suffix(".gff")
        interval = gff.load_gene_intervals(
            gff_path, [specs["oprD"].locus_tag]
        )[specs["oprD"].locus_tag]
        assert interval.locus_tag == "PA0958"
        assert interval.name == "oprD"
        # A whole gene, not a fragment.
        assert 90 <= interval.end - interval.start + 1 <= 15000

    def test_mexz_sits_beside_the_mexxy_pump_it_regulates(self, config):
        """The adjacency the locus resolver relied on to identify it."""
        from papipeline.config.loader import load_regulators

        specs = load_regulators(Path("config/regulators.tsv"))
        gff_path = config.reference_fasta().with_suffix(".gff")
        interval = gff.load_gene_intervals(
            gff_path, [specs["mexZ"].locus_tag]
        )[specs["mexZ"].locus_tag]
        assert interval.locus_tag == "PA2020"
        mexxy = gff.load_gene_intervals(gff_path, ["PA2018", "PA2019"])
        gaps = [
            min(abs(i.start - f.end), abs(f.start - i.end))
            for f in mexxy.values()
            for i in [interval]
        ]
        assert min(gaps) < 2000, f"mexZ is {min(gaps)} bp from the mexXY genes"