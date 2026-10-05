"""The GPA table must survive a write/read round trip. This is the test that
catches a reader and a writer disagreeing about the format.

**Why a round-trip test and not a shape test.** A test asserting the header is
`gene, sample_id, present` passes even if the reader cannot parse what the writer
produces. Asserting the columns is checking the declaration, not the data. The
failure this guards against is a writer and a reader that are each individually
reasonable and jointly wrong - which is exactly what happened when the reshape was
first attempted: the writer was moved to long format, the reader was left reading
wide, `GPA_COLUMNS` said one thing, and a narrow `-k` selector reported 78
passes because nothing asserted the contract.

A round trip is the only assertion that notices. It fails the moment the two ends
disagree, regardless of which shape either one believes in.

**Written before any format code changed, and watched to fail.** If this passes
on an unmodified tree, it is testing nothing - that is the first thing to check,
not the last.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.stages import pangenome as stage

SAMPLES = ("S1", "S2", "S3")

# Deliberately not symmetric: one gene in every sample (core), one in two
# (accessory), one in exactly one (genome-specific). A round trip that passes
# these cannot be collapsing the presence set.
PRESENCE = {
    "gene_core": {"S1", "S2", "S3"},
    "gene_accessory": {"S1", "S2"},
    "gene_specific": {"S3"},
}


@pytest.fixture()
def pangenome():
    return stage.partition(SAMPLES, PRESENCE)


class TestRoundTrip:
    def test_presence_survives_write_then_read(self, pangenome, tmp_path):
        paths = stage.write_outputs(pangenome, tmp_path)
        samples, presence = stage.read_gene_presence_absence(
            paths["gene_presence_absence"]
        )
        assert samples == SAMPLES, "sample column set changed across the round trip"
        assert presence == PRESENCE, (
            "the presence mapping did not survive the round trip; a reader and a "
            "writer disagree about the GPA format"
        )

    def test_the_partition_survives_the_round_trip(self, pangenome, tmp_path):
        """Core/accessory is our partition, derived from presence.

        Deriving it from a re-read table is a stronger check than comparing the
        object to itself: it proves the written table still carries enough
        information to rebuild the same partition.
        """
        paths = stage.write_outputs(pangenome, tmp_path)
        samples, presence = stage.read_gene_presence_absence(
            paths["gene_presence_absence"]
        )
        rebuilt = stage.partition(samples, presence)
        assert rebuilt.core == pangenome.core
        assert rebuilt.accessory == pangenome.accessory
        assert len(rebuilt.core) == 1, "one gene is in all three samples"
        assert len(rebuilt.accessory) == 2

    def test_a_genome_specific_gene_is_not_read_as_core(self, pangenome, tmp_path):
        """The specific failure mode of a lossy format.

        `gene_specific` is in one isolate of three. If the format cannot express
        *which* isolate, it collapses to `n_samples`/`frequency`, and a gene that
        is absent from two samples starts looking present in all three - the
        exact inversion D6 depends on being able to avoid.
        """
        paths = stage.write_outputs(pangenome, tmp_path)
        _, presence = stage.read_gene_presence_absence(paths["gene_presence_absence"])
        assert presence["gene_specific"] == {"S3"}
        assert "S1" not in presence["gene_specific"]


class TestDeclaredShape:
    """The contract, pinned so a silent reversion is visible."""

    def test_gpa_columns_are_the_long_shape(self):
        assert stage.GPA_COLUMNS == ("gene", "sample_id", "present")

    def test_the_written_header_matches_the_declared_columns(self, pangenome, tmp_path):
        paths = stage.write_outputs(pangenome, tmp_path)
        header = [
            line
            for line in paths["gene_presence_absence"]
            .read_text(encoding="utf-8")
            .splitlines()
            if line and not line.startswith("#")
        ][0]
        assert tuple(header.split("\t")) == stage.GPA_COLUMNS

    def test_the_reader_does_not_depend_on_column_position(self, pangenome, tmp_path):
        """Columns are located by name.

        pyseer and roary both emit a leading metadata column before the genome
        columns, and treating those as genomes is the easy way to compute a core
        count that is wrong but plausible.
        """
        paths = stage.write_outputs(pangenome, tmp_path)
        samples, presence = stage.read_gene_presence_absence(
            paths["gene_presence_absence"]
        )
        assert set(samples) == set(SAMPLES)
        assert "sample_id" not in samples and "gene" not in samples
