"""Stage 2 must not reject a genome because its two Bakta tables differ in length.

**The bug this locks down.** `_run_by_calling_tool` compared the GFF and TSV
feature counts with `len(gff) != len(features)` and rejected on any difference.
Real data rejects **every** prepared genome:

    Stage 2: 967 of 967 cohort members were not annotated
             (957 no assembly, 10 annotated but unusable):
             PDT000034122.1=GFF/TSV disagree (5912 vs 5913)

So stage 2 ran Bakta on all ten prepared assemblies, got ~5,900 features from
each, and then discarded every result. `02_annotation_summary.tsv` recorded
`n_records = 0` for all 967 cohort members - a plausible-looking zero rather
than an error, which is why it survived a completed run unnoticed. Everything
downstream - regulators, mechanisms, GWAS - would have been built on nothing.

**Why the counts differ, and why that is correct.** Verified against the real
output: the 5,912 GFF feature IDs are a *perfect subset* of the TSV's 5,944.
Bakta's `inference.tsv` carries rows the GFF filter legitimately excludes -
derived features, and rows whose gene id is `-`. The pair is consistent; they are
simply not the same length.

**The invariant that actually catches a broken pair** is subset containment: if
the GFF and TSV came from different genomes, or one were truncated, some GFF ID
would have no TSV counterpart. Equality proves nothing that containment does not
prove better, and it rejects valid data besides.
"""

from __future__ import annotations

import pytest

from papipeline.stages import annotation as stage_annotation


GFF_HEADER = "##gff-version 3\n"
GFF_ROW = (
    "contig_1\tBakta\t{feature}\t{start}\t{stop}\t.\t+\t.\t"
    "ID={gene};Name={name}\n"
)
#: Bakta writes the header on a `#`-prefixed line, distinct from the `# `
#: human-readable preamble above it. Column order matches BAKTA_TSV_COLUMNS.
TSV_HEADER = (
    "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n"
)
TSV_ROW = "{contig}\t{feature}\t{start}\t{stop}\t-\t{gene}\t{gene}\tnp\t-\n"


def _gff(rows):
    return GFF_HEADER + "".join(rows)


def _tsv(rows):
    return TSV_HEADER + "".join(rows)


class TestTheCrossCheckIsContainmentNotEquality:
    def test_a_tsv_with_extra_derived_rows_is_accepted(self):
        """The real case: 5,912 GFF features inside 5,944 TSV rows.

        Rejecting this is what zeroed stage 2 for the whole cohort.
        """
        gff = parse = stage_annotation.parse_bakta_gff(
            _gff([GFF_ROW.format(feature="CDS", start=i, stop=i + 100,
                                 gene=f"ABLLKF_{i:05d}", name=f"gene{i}")
                  for i in range(1, 4)])
        )
        tsv = stage_annotation.parse_bakta_tsv(
            _tsv([TSV_ROW.format(contig="contig_1", feature="cds", start=i,
                                 stop=i + 100, gene=f"ABLLKF_{i:05d}")
                  for i in range(1, 4)]
                 + [TSV_ROW.format(contig="contig_1", feature="misc",
                                   start=900, stop=950, gene="-")])
        )
        assert len(gff) == 3 and len(tsv) == 4
        assert stage_annotation.gff_features_are_covered(gff, tsv) is True, (
            "every GFF feature has a TSV counterpart, so the pair is consistent "
            "even though the lengths differ"
        )

    def test_a_gff_feature_missing_from_the_tsv_is_rejected(self):
        """The failure the cross-check exists for.

        A GFF from a different genome, or a truncated table, leaves a GFF ID with
        no TSV counterpart. Equality also catches this, but so does containment
        - without rejecting valid data.
        """
        gff = stage_annotation.parse_bakta_gff(
            _gff([GFF_ROW.format(feature="CDS", start=i, stop=i + 100,
                                 gene=f"ABLLKF_{i:05d}", name=f"gene{i}")
                  for i in range(1, 4)])
        )
        tsv = stage_annotation.parse_bakta_tsv(
            _tsv([TSV_ROW.format(contig="contig_1", feature="cds", start=i,
                                 stop=i + 100, gene=f"ABLLKF_{i:05d}")
                  for i in (1, 2)])  # gene 3 missing
        )
        assert stage_annotation.gff_features_are_covered(gff, tsv) is False

    def test_rows_with_no_gene_id_do_not_cause_a_false_pass(self):
        """A `-` id on both sides must not make containment vacuously true."""
        gff = stage_annotation.parse_bakta_gff(
            _gff([GFF_ROW.format(feature="CDS", start=1, stop=100,
                                 gene="ABLLKF_00001", name="g1")])
        )
        tsv = stage_annotation.parse_bakta_tsv(
            _tsv([TSV_ROW.format(contig="contig_1", feature="misc", start=900,
                                 stop=950, gene="-")])
        )
        assert stage_annotation.gff_features_are_covered(gff, tsv) is False

    def test_real_counts_differ_and_are_still_consistent(self):
        """The specific numbers from the smoke cohort.

        Recorded so a future change cannot quietly re-tighten this into equality
        and silently re-empty stage 2.
        """
        gff = stage_annotation.parse_bakta_gff(
            _gff([GFF_ROW.format(feature="CDS", start=i, stop=i + 100,
                                 gene=f"ABLLKF_{i:05d}", name=f"g{i}")
                  for i in range(1, 5913)])
        )
        tsv = stage_annotation.parse_bakta_tsv(
            _tsv([TSV_ROW.format(contig="contig_1", feature="cds", start=i,
                                 stop=i + 100, gene=f"ABLLKF_{i:05d}")
                  for i in range(1, 5913)]
                 + [TSV_ROW.format(contig="contig_1", feature="misc",
                                   start=9, stop=9, gene="-")] * 32)
        )
        assert len(gff) == 5912 and len(tsv) == 5944
        assert stage_annotation.gff_features_are_covered(gff, tsv) is True


class TestStageTwoProducesRecords:
    def test_the_rejection_reason_is_containment(self, config, monkeypatch):
        """Regression on the message the run logs.

        "GFF/TSV disagree" is what the completed 11:54 run printed, and it reads
        like a tool failure. It was the pipeline rejecting its own valid input.
        """
        from pathlib import Path

        source = (
            Path(__file__).resolve().parents[2] / "papipeline" / "stages" / "annotation.py"
        ).read_text(encoding="utf-8")
        assert "GFF/TSV disagree" not in source, (
            "the stage still rejects on a count mismatch; the invariant is "
            "containment, and this message misattributes it to the tool"
        )