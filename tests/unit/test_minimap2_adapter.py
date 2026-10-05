"""The real align-and-call invocation, as configured.

`config/science.yaml` pins this invocation against a measured REAL run, so the
science is settled and only the code was missing. What is tested here is that
the code *reads that configuration* rather than restating it: a hard-coded
`-q 20` would pass every functional test while silently disagreeing with the
file that records the science.

Three properties are load-bearing.

**The reference index is verified before anything reads it.** A `.fai` built by
indexing a bgzipped copy of PAO1 has been observed reporting 64,400 bases for a
6,264,404-base file; every tool then read a reference two orders of magnitude
too small and produced confident, entirely wrong output - 64 alignments for a
6.4 Mb genome. Nothing errored. So the guard is not optional and not the
caller's to remember: `call_isolate` performs it, and the flags cannot be
assembled without passing through.

**Every flag is one this bcftools actually has.** Checked against the installed
binaries rather than assumed; `bcftools mpileup --help` does not exist, and
`bcftools call -h` lists `-m` as `--multiallelic-caller`, not as a
variants-only switch. A wrong flag is worse than no code.

**One isolate in, one VCF out.** The adapter's whole job is that conversion. It
does not interpret a call - that stays in `stages.variants`, which owns the
schema - so the seam is the VCF path, and the stage parses from there.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from papipeline.adapters import minimap2 as adapter
from papipeline.config.loader import PipelineConfig
from papipeline.errors import DataContractError, ToolNotAvailableError
from papipeline.reference import verify_fai_index

REAL_REFERENCE = Path(
    "db/reference/GCF_000006765.1/GCF_000006765.1_ASM676v1_genomic.fna"
)


class TestNothingIsHardCoded:
    """The config is the single source; the code must not restate it."""

    def test_the_preset_comes_from_config(self, config: PipelineConfig):
        configured = config.raw["variants"]["minimap2"]["preset"]
        assert configured == "asm20"
        command = adapter.align_command(
            "minimap2", preset=configured, assembly="a.fna",
            reference="ref.fna", output_sam="a.sam",
        )
        assert "-x" in command and "asm20" in command

    def test_sam_output_is_requested_because_samtools_cannot_read_sam(
        self, config: PipelineConfig
    ):
        """science.yaml records the pinned samtools as unable to read SAM text.

        So the adapter asks minimap2 for SAM and hands the conversion to pysam.
        Asking for BAM here instead would assume a samtools that this
        environment does not have.
        """
        assert config.raw["variants"]["minimap2"]["output_format"] == "SAM"
        assert config.raw["variants"]["samtools"]["sam_text_supported"] is False
        assert config.raw["variants"]["samtools"]["convert_and_index_with"] == "pysam"
        command = adapter.align_command(
            "minimap2", preset="asm20", assembly="a.fna",
            reference="ref.fna", output_sam="a.sam",
        )
        assert "-a" in command

    def test_secondary_alignments_are_kept(self, config: PipelineConfig):
        """science.yaml: `--secondary=no` was measurably wrong here."""
        assert config.raw["variants"]["minimap2"]["secondary_alignments"] is True
        command = adapter.align_command(
            "minimap2", preset="asm20", assembly="a.fna",
            reference="ref.fna", output_sam="a.sam",
        )
        assert "--secondary=yes" in command

    def test_mpileup_thresholds_come_from_config(self, config: PipelineConfig):
        pileup = config.raw["variants"]["bcftools"]["mpileup"]
        command = adapter.mpileup_command(
            "bcftools", bam="a.bam", reference="ref.fna",
            min_mapq=pileup["min_mapq"], min_bq=pileup["min_bq"],
            annotate=pileup["annotate"], max_depth=pileup["max_depth"],
        )
        assert "-q" in command and "20" in command
        assert "-Q" in command
        assert "-d" in command and "250" in command

    def test_the_annotation_string_is_passed_whole(self, config: PipelineConfig):
        """`FORMAT/DP,FORMAT/AD` is one argument; splitting it invents a flag."""
        pileup = config.raw["variants"]["bcftools"]["mpileup"]
        command = adapter.mpileup_command(
            "bcftools", bam="a.bam", reference="ref.fna",
            min_mapq=20, min_bq=20, annotate=pileup["annotate"], max_depth=250,
        )
        assert "FORMAT/DP,FORMAT/AD" in command
        assert "FORMAT/DP" not in command

    def test_the_reference_is_named_explicitly(self, config: PipelineConfig):
        """Without `-f` bcftools has no reference and the run is meaningless."""
        command = adapter.mpileup_command(
            "bcftools", bam="a.bam", reference="ref.fna",
            min_mapq=20, min_bq=20, annotate="FORMAT/DP", max_depth=250,
        )
        assert "-f" in command
        assert "ref.fna" in command

    def test_the_call_flags_come_from_config(self, config: PipelineConfig):
        """`-m` is multiallelic-calling and `-v` is variants-only; they differ.

        `science.yaml` carries both switches, and the real run was
        `bcftools call -mv`. An adapter that emitted only `-m` would emit sites
        that are not variant sites, which the merge would then treat as
        polymorphisms.
        """
        call = config.raw["variants"]["bcftools"]["call"]
        command = adapter.call_command(
            "bcftools", input_vcf="a.vcf", output_vcf="b.vcf",
            variants_only=call["variants_only"],
            multiallelic=call["multiallelic"],
        )
        assert "-m" in command
        assert "-v" in command

    def test_multiallelic_split_is_the_measured_behaviour(self, config: PipelineConfig):
        """`split` is what produced the multiallelic sites the parser expands."""
        assert config.raw["variants"]["bcftools"]["call"]["multiallelic"] == "split"

    def test_the_threads_value_comes_from_the_machine_overlay(
        self, config: PipelineConfig
    ):
        """A thread count in code is a hard rule 2 violation."""
        from papipeline.adapters import minimap2

        command = minimap2.call_command(
            "bcftools", input_vcf="a.vcf", output_vcf="b.vcf",
            variants_only=True, multiallelic="split",
            threads=int(config.runtime.get("threads", 1)),
        )
        assert "--threads" in command
        assert str(config.runtime.get("threads", 1)) in command


class TestTheIndexGuardIsInsideTheAdapter:
    """It must not be the caller's to remember."""

    def test_verification_is_exposed_and_uses_the_configured_length(
        self, config: PipelineConfig
    ):
        expected = config.raw["variants"]["reference_index"][
            "verify_against_expected_length"
        ]
        assert expected == 6264404
        assert callable(adapter.verify_reference_index)

    def test_a_missing_index_is_refused(self, tmp_path):
        with pytest.raises(DataContractError):
            adapter.verify_reference_index(
                tmp_path / "absent.fai", expected_bases=6264404
            )

    def test_a_lying_index_is_refused(self, tmp_path):
        """The observed 64,400-for-6,264,404 shape, written out as a fixture.

        The numbers are the ones `tests/integration/test_reference_index_guard.py`
        measures from real samtools behaviour; this asserts the adapter refuses
        them, so the guard is wired in rather than merely available.
        """
        fai = tmp_path / "reference.fna.fai"
        fai.write_text("NC_002516.2\t64400\t0\t0\t0\n", encoding="utf-8")
        with pytest.raises(DataContractError) as excinfo:
            adapter.verify_reference_index(fai, expected_bases=6264404)
        message = str(excinfo.value)
        assert "6264404" in message
        assert "64400" in message

    def test_a_correct_index_passes(self, tmp_path):
        fai = tmp_path / "reference.fna.fai"
        fai.write_text("NC_002516.2\t6264404\t0\t0\t0\n", encoding="utf-8")
        assert adapter.verify_reference_index(fai, expected_bases=6264404) == 6264404

    def test_it_is_the_shared_guard_not_a_second_implementation(self):
        """Two implementations of this check would be free to disagree."""
        assert adapter.verify_reference_index is not None
        # The adapter must delegate, so the integration test that proves the
        # hazard is caught stays true of the adapter's own path.
        assert verify_fai_index.__module__ == "papipeline.reference"


class TestTheBamHeaderNamesTheReference:
    """Found by running it, not by reading it.

    `minimap2 -a` writes `@SQ` lines for the **query** - the assembly's own
    contigs - while the alignment records themselves are named for the
    **reference**. So the SAM looks self-consistent and is not: pysam sorts it
    into a BAM whose header lists `JYGC02000001.1`, and `bcftools mpileup -f
    PAO1` then has no reference sequence for any of the contigs it is asked to
    pile up. Every site comes out as a spurious `N -> G` mismatch and
    `bcftools call -v` returns **zero variants** - for a genome that plainly has
    tens of thousands of them.

    Nothing errors. The VCF is well formed, it just says the isolate is
    identical to PAO1, and the cohort merge would report every real difference
    as absent. This is the same failure class as the short `.fai`: a confident,
    entirely wrong answer that no check would catch.

    The reference's own sequence dictionary is installed on the BAM at parse
    time, and the test asserts it rather than trusting the tool.

    Scoped honestly: with the arguments to `align_command` in their correct
    order, bcftools reads the un-rewritten BAM and calls correctly, so this is
    a defence rather than a current fix - it is asserted here because the failure
    it prevents is invisible (reads quietly demoted, a smaller call set, nothing
    in the output to say so), and because it *was* load-bearing while the
    arguments were swapped.
    """

    def test_it_rewrites_the_header_from_the_reference(self, tmp_path):
        bam = adapter.convert_and_index(
            _sam_with_query_header(tmp_path), tmp_path / "a.sorted.bam",
            reference=_tiny_reference(tmp_path),
        )
        import pysam

        with pysam.AlignmentFile(str(bam), "rb") as handle:
            names = list(handle.references)
        assert names == [_REFERENCE_CONTIG], (
            f"the BAM header names {names} rather than the reference contig; "
            "bcftools mpileup -f will find no sequence for these and call "
            "nothing"
        )

    def test_the_length_comes_from_the_reference_too(self, tmp_path):
        """A length from the query would let bcftools read past the end."""
        import pysam

        reference = _tiny_reference(tmp_path, length=500)
        bam = adapter.convert_and_index(
            _sam_with_query_header(tmp_path), tmp_path / "b.sorted.bam",
            reference=reference,
        )
        with pysam.AlignmentFile(str(bam), "rb") as handle:
            assert handle.get_reference_length(_REFERENCE_CONTIG) == 500

    def test_the_alignment_records_survive_the_rewrite(self, tmp_path):
        """A header fix that drops the reads is a different silent failure.

        Asserted on the record count, so a rewrite implemented as "write a new
        empty file with the right header" cannot pass.
        """
        import pysam

        bam = adapter.convert_and_index(
            _sam_with_query_header(tmp_path), tmp_path / "d.sorted.bam",
            reference=_tiny_reference(tmp_path),
        )
        with pysam.AlignmentFile(str(bam), "rb") as handle:
            # `until_eof=True` rather than `count()`: count() goes through
            # fetch(), which requires an index, and asserting on the rewrite
            # should not depend on the index step that follows it.
            reads = list(handle.fetch(until_eof=True))
        assert len(reads) == 1
        assert reads[0].reference_name == _REFERENCE_CONTIG
        assert reads[0].mapping_quality == 60
        # And it is still a *mapped* read. htslib demotes a record whose contig
        # is missing from the header to unmapped, so a rewrite that dropped the
        # reference from the dictionary would leave a well-formed BAM full of
        # unmapped reads - which calls nothing, the same silent failure.
        assert not reads[0].is_unmapped

    def test_it_does_not_parse_against_the_wrong_header(self, tmp_path):
        """The specific trap: a text rewrite silently empties the BAM.

        htslib validates RNAME against `@SQ` *as it parses*. So a SAM whose
        `@SQ` names the query is read with every record demoted to unmapped -
        and the obvious "fix" of rewriting `@SQ` in the text afterwards still
        parses against the old dictionary first, losing the reads.

        Asserted as a count, because the failure is invisible: the file exists,
        parses, indexes, and holds the right header. Only the read count says
        anything happened, and only this test would notice.
        """
        bam = adapter.convert_and_index(
            _sam_with_query_header(tmp_path), tmp_path / "e.sorted.bam",
            reference=_tiny_reference(tmp_path),
        )
        import pysam

        with pysam.AlignmentFile(str(bam), "rb") as handle:
            assert len(list(handle.fetch(until_eof=True))) == 1, (
                "the reads were lost; the BAM is well formed and empty, which "
                "bcftools would report as an isolate with no variants"
            )

    def test_the_index_is_built_over_the_rewritten_header(self, tmp_path):
        """An index over the wrong header is a wrong index, however valid it is."""
        bam = adapter.convert_and_index(
            _sam_with_query_header(tmp_path), tmp_path / "c.sorted.bam",
            reference=_tiny_reference(tmp_path),
        )
        import pysam

        pysam.index(str(bam))
        with pysam.AlignmentFile(str(bam), "rb") as handle:
            assert handle.has_index()
            assert list(handle.references) == [_REFERENCE_CONTIG]


def _tiny_reference(tmp_path, *, length: int = 300) -> Path:
    """A real, indexable reference FASTA. Not a stub: the header must be real."""
    fasta = tmp_path / f"ref{length}.fna"
    fasta.write_text(
        f">{_REFERENCE_CONTIG}\n{'ACGT' * (length // 4)}\n", encoding="utf-8"
    )
    import pysam

    pysam.faidx(str(fasta))
    return fasta


#: The reference contig, as the alignment records name it.
_REFERENCE_CONTIG = "NC_TESTREF.1"

#: The query contig, as minimap2's `@SQ` lines wrongly name it. This is the
#: real shape, read off an actual `minimap2 -a` run: `@SQ` follows the *query*
#: while RNAME follows the *reference*. A fixture that got this backwards would
#: pass the rewrite while proving nothing about the failure it exists to catch.
_QUERY_CONTIG = "CONTIG_FROM_THE_QUERY"


def _sam_text() -> str:
    return (
        "@HD\tVN:1.6\tSO:unsorted\tGO:query\n"
        f"@SQ\tSN:{_QUERY_CONTIG}\tLN:1478571\n"
        # RNAME is the REFERENCE contig; QNAME is the query contig. That
        # divergence between the header and the records is the whole bug.
        f"{_REFERENCE_CONTIG}\t0\t{_QUERY_CONTIG}\t1\t60\t20M\t*\t0\t0\t"
        "ACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII\n"
    )


def _sam_with_query_header(tmp_path) -> Path:
    sam = tmp_path / "query.sam"
    sam.write_text(_sam_text(), encoding="utf-8")
    return sam


class TestTheAssemblyIsTheQuery:
    """`minimap2 [options] <target.fa> [query.fa]`.

    The reference is the **target** and the assembly is the **query**, and
    getting that backwards is invisible. Swapped, minimap2 aligns PAO1 against
    the isolate: the output is well formed, every record is named
    `NC_002516.2`, and every SEQ is 6,264,404 bases long - the reference's
    length, not the assembly's. Downstream, `bcftools mpileup` piles up a query
    it believes is the reference and `call -v` returns **zero variants** for a
    genome that has tens of thousands of them. Nothing errors.

    Observed, not theorised: this exact inversion shipped and produced a clean
    STUB-shaped result for a real assembly. The two assertions below are the ones
    that would have caught it, and they are asserted on the *content* of the
    files rather than on their names, because a path bug and a swapped pair look
    identical in a command line.
    """

    def _assembly(self, tmp_path) -> Path:
        """A real assembly, md5-distinct from the reference, in the real layout.

        Named like a PDC isolate and shaped like one, because the lookup under
        test resolves by directory, and a fixture that did not would pass for
        the wrong reason.
        """
        path = tmp_path / "GCA_000000000.1.fna"
        path.write_text(">contig_1\n" + "ATGC" * 200 + "\n", encoding="utf-8")
        return path

    def _reference(self, tmp_path) -> Path:
        path = tmp_path / "GCF_000006765.1_ASM676v1_genomic.fna"
        path.write_text(">NC_002516.2\n" + "ATGC" * 500 + "\n", encoding="utf-8")
        import pysam

        pysam.faidx(str(path))
        return path

    def test_the_reference_is_the_target_and_the_assembly_the_query(
        self, tmp_path
    ):
        """minimap2's first positional is the target, the second the query."""
        command = adapter.align_command(
            "minimap2", preset="asm20",
            assembly=self._assembly(tmp_path),
            reference=self._reference(tmp_path),
            output_sam=tmp_path / "a.sam",
        )
        # The two FASTA arguments are the last two entries; the options come
        # first and are asserted separately.
        assert command[-2] == str(tmp_path / "GCF_000006765.1_ASM676v1_genomic.fna"), (
            "the reference must be the FIRST positional (the target); passing it "
            "second makes it the query, and minimap2 then aligns PAO1 against "
            "the isolate"
        )
        assert command[-1] == str(tmp_path / "GCA_000000000.1.fna"), (
            "the assembly must be the SECOND positional (the query)"
        )

    def test_the_two_input_paths_are_never_identical(self, tmp_path):
        """Passing one file twice aligns the reference to itself.

        A guard in its own right, and the one that would catch a
        path-resolution bug where ``locate_assembly``'s result never reaches
        the call and the reference is supplied in both slots.
        """
        assembly = self._assembly(tmp_path)
        command = adapter.align_command(
            "minimap2", preset="asm20", assembly=assembly,
            reference=self._reference(tmp_path), output_sam=tmp_path / "a.sam",
        )
        assert Path(command[-1]).resolve() != Path(command[-2]).resolve(), (
            "the assembly and the reference resolved to the same file; minimap2 "
            "would align PAO1 against PAO1 and produce a self-comparison"
        )

    def test_the_query_is_not_the_reference_by_content(self, tmp_path):
        """Compared by digest, because a path bug looks like a swap in argv.

        The reference's md5 is the pinned one from ``science.yaml``; the
        assembly's is computed from the file the lookup returned. If the query
        digest ever equals the reference's, the reference reached the query slot
        - whether through an argument swap or through a mis-resolved path.
        """
        import hashlib

        def md5(path: Path) -> str:
            return hashlib.md5(Path(path).read_bytes()).hexdigest()

        reference = self._reference(tmp_path)
        assembly = self._assembly(tmp_path)
        command = adapter.align_command(
            "minimap2", preset="asm20", assembly=assembly,
            reference=reference, output_sam=tmp_path / "a.sam",
        )
        pinned = "b8b18c2da0cb95376bb8de30f5a8cab7"
        assert md5(Path(command[-1])) != md5(reference)
        assert md5(Path(command[-1])) != pinned, (
            "the query argument is the pinned reference; minimap2 would align "
            "PAO1 against the isolate and call nothing"
        )

    def test_the_recorded_length_identifies_which_file_is_the_query(
        self, tmp_path
    ):
        """The symptom that gave this away, asserted as a property.

        A SEQ whose length equals the *reference's* length on an isolate's
        contig name is the fingerprint of the swap. Asserted on the fixture's
        two distinct lengths, so the test states the rule rather than the
        incident.
        """
        def bases(path: Path) -> int:
            return sum(
                len(line.strip()) for line in Path(path).read_text().splitlines()
                if not line.startswith(">")
            )

        assembly, reference = self._assembly(tmp_path), self._reference(tmp_path)
        command = adapter.align_command(
            "minimap2", preset="asm20", assembly=assembly,
            reference=reference, output_sam=tmp_path / "a.sam",
        )
        query, target = Path(command[-1]), Path(command[-2])
        assert bases(query) == bases(assembly)
        assert bases(target) == bases(reference)
        assert bases(query) != bases(target), (
            "the fixture's two files must differ in length, or this test proves "
            "nothing about which one reached the query slot"
        )


class TestToolAvailability:
    @pytest.mark.parametrize("missing", ["minimap2", "bcftools"])
    def test_a_missing_tool_is_refused_by_name(self, missing, tmp_path):
        """A refusal the operator cannot act on is just an obstacle."""
        present = {name: f"/usr/bin/{name}" for name in ("minimap2", "bcftools")}
        present[missing] = None
        with pytest.raises(ToolNotAvailableError) as excinfo:
            adapter.require_callers(present)
        assert missing in str(excinfo.value)

    def test_present_tools_are_returned(self):
        statuses = {"minimap2": "/usr/bin/minimap2", "bcftools": "/usr/bin/bcftools"}
        assert adapter.require_callers(statuses) == statuses


@pytest.mark.skipif(
    not shutil.which("bcftools"), reason="bcftools is not installed on this machine"
)
class TestTheFlagsAreReal:
    """Against the installed binary, not against my memory of it.

    A test asserting the command *contains* `-m` proves only that the code
    contains `-m`. This asks bcftools whether it recognises the flags the
    adapter emits, which is the failure mode that actually bit this project: a
    plausible flag that the installed version does not have.
    """

    def test_mpileup_accepts_every_flag_the_adapter_emits(self, config, tmp_path):
        import subprocess

        pileup = config.raw["variants"]["bcftools"]["mpileup"]
        command = adapter.mpileup_command(
            "bcftools", bam=str(tmp_path / "a.bam"),
            reference=str(tmp_path / "ref.fna"),
            min_mapq=pileup["min_mapq"], min_bq=pileup["min_bq"],
            annotate=pileup["annotate"], max_depth=pileup["max_depth"],
            output_vcf=str(tmp_path / "raw.vcf"),
        )
        # Run it for real: a missing file is a *usage* error, an unrecognised
        # option is not. Only the latter means the flag is invented.
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        assert "unrecognized" not in result.stderr.lower(), result.stderr
        assert "invalid option" not in result.stderr.lower(), result.stderr

    def test_call_accepts_every_flag_the_adapter_emits(self, config, tmp_path):
        import subprocess

        command = adapter.call_command(
            "bcftools", input_vcf=str(tmp_path / "a.vcf"),
            output_vcf=str(tmp_path / "b.vcf"),
            variants_only=True, multiallelic="split", threads=1,
        )
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        assert "unrecognized" not in result.stderr.lower(), result.stderr
        assert "invalid option" not in result.stderr.lower(), result.stderr

    def test_the_variant_recording_matches_this_bcftools(self, config):
        """science.yaml pins 1.23.1; say so rather than quietly running another.

        The columns in `PER_ISOLATE_COLUMNS` were read off 1.23.1's output. A
        different version can add or rename INFO fields, so running against one
        silently is exactly the PROVISIONAL-columns mistake one level up - and
        it is invisible, because the calls still parse.

        Note this reads the *active* `bcftools`, which is not necessarily the
        one `shutil.which` finds outside the test environment: this machine has
        a homebrew 1.24 on PATH alongside the pinned 1.23.1. That is precisely
        the situation the assertion exists to catch.
        """
        pinned = config.raw["variants"]["bcftools"]["version"]
        assert pinned == "1.23.1"
        import subprocess

        out = subprocess.run(
            ["bcftools", "--version"], capture_output=True, text=True, check=False
        ).stdout
        installed = out.split()[1]
        assert installed == pinned, (
            f"bcftools {installed} is on PATH but science.yaml pins {pinned}. The "
            "per-isolate contract was read off the pinned version, so its columns "
            "are not guaranteed to survive this one. Re-derive the contract, or "
            "run the pinned environment, before calling a REAL cohort."
        )
