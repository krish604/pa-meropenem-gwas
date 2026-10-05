"""Stage 6's mpileup cost is bounded, and the bound is visible on the command line.

The defect: in the round-12 REAL run over 10 isolates, `bcftools mpileup` was
SIGKILLed (`returncode=-9`) on `PDT000034122.1`, with the argv
`... -d 250 -m 1 -o .../PDT000034122.1.raw.vcf ...sorted.bam`. The isolate then
stayed a cohort member with no calls, so `variants.tsv` and `06_regulators.tsv`
carried 9 rows against a cohort of 10.

**What the investigation actually found, and it is not what the argv suggests.**
`-d` is bcftools' only memory lever, but it engages ONLY where true depth exceeds
it. Measured against the pinned 6,264,404 bp PAO1 (peak RSS, `/usr/bin/time -l`):

* depth 1 -- the geometry `minimap2 -a --secondary=yes` normally produces --
  94 MB at `-d 250`, and `-d 1` is indistinguishable. **The cap does not bind, so
  it is not what an ordinary isolate costs.**
* 60 x 500 kb contig records stacked on one locus -- true depth 60, which is what
  a repeat locus under `--secondary=yes` looks like -- 127 MB at `-d 250`,
  77 MB at `-d 20`, 52 MB at `-d 1`.

Every other candidate lever was measured and rejected: `-r`/`-t` region splitting
(63.2 MB whole vs 62.9 MB for a 20 kb region -- no reduction), `-L/--max-idepth`
(62.5 / 62.8 / 62.6 MB at 250 / 20 / 1 -- no effect), `-O z` and `-O b`
(127.4 / 127.9 MB vs 125.2 -- slightly WORSE, compression buffers), piping to
`bcftools call` instead of `-o` (126.8 MB -- no reduction), `--threads 2` (88.5 vs
94.5 MB -- within noise). **`--max-BP` / `-b` does not exist in bcftools 1.23.1**
(absent from the usage block; `--max-BP` is rejected as unrecognized; `grep -c
"max-BP"` on the binary is 0), so the flag the defect report reached for is not
available and is not used.

**So `max_depth` is not lowered here.** `-d` is a SUBSAMPLING cap: past it mpileup
keeps a subset of reads, moving `FORMAT/DP` and `FORMAT/AD` and therefore which
alleles `call -mv` emits. `QUAL` is constant at 30.4183 across every called site,
so DP/AD is not recoverable downstream and the call genuinely moves. That is a
SCIENCE decision about the callable variant set. This file's job is to make the
bound EXPLICIT, keep it out of a hidden literal, and pin the two things that
would otherwise be lost silently.

**What "the callable variant set is preserved" means here, precisely.** It means
the emitted argv differs from the pre-fix argv ONLY by the addition of `-O z`,
and the output is renamed to carry the `.gz` extension. Everything that decides
*which variants are called* is byte-identical: the reference (`-f`), the mapping
quality floor (`-q`), the base quality floor (`-Q`), the annotation set (`-a`),
the depth cap (`-d 250`, unchanged value) and `min_ireads` (`-m 1`, R14).
`-O` selects an on-disk encoding and is read by no downstream logic. That claim is
not asserted from the source; it is MEASURED, by
`test_the_callable_variant_set_is_unchanged_by_this_fix`, which runs the real
pinned bcftools over a synthetic BAM twice -- once with the pre-fix argv and once
with the post-fix argv -- and requires the `call -mv` output to be identical.
Compression that changed a single DP, AD or allele would change that output.

**R12.** Every test here either builds the argv through the real
`mpileup_command`/`call_settings` (so it exercises production code with
synthetic inputs) or drives it with an INJECTED FAKE `bcftools` runner on `PATH`.
No bioinformatics tool is executed except the pinned `bcftools` in the one test
that says so explicitly and skips when it is absent.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import List

import pytest

from papipeline.adapters import minimap2 as adapter
from papipeline.config.loader import (
    DEFAULT_MPILEUP_MAX_DEPTH,
    DEFAULT_MPILEUP_MIN_IREADS,
    DEFAULT_MPILEUP_OUTPUT_TYPE,
    PipelineConfig,
)
from papipeline.stages import variants as stage

#: The pinned binary's own usage block for the flags this fix relies on, kept as
#: a literal so the assertions are against the tool's words and not against this
#: file's memory of them. `bcftools mpileup --help` does NOT exist on 1.23.1 --
#: it exits 1 with "unrecognized option `--help'" -- so a bare `bcftools
#: mpileup` is the authoritative text.
MPILEUP_DEPTH_USAGE = (
    "-d, --max-depth INT     Max raw per-file depth; avoids excessive memory "
    "usage [250]"
)
MPILEUP_OUTPUT_TYPE_USAGE = (
    "  -O, --output-type TYPE  'b' compressed BCF; 'u' uncompressed BCF;\n"
    "                          'z' compressed VCF; 'v' uncompressed VCF; "
    "0-9 compression level [v]"
)

#: The flags the pinned binary documents. The regression test compares every
#: emitted flag against this, so a future flag that 1.23.1 does not have is
#: caught here rather than as a "call failed" in a REAL run.
MPILEUP_FLAGS_THAT_EXIST = frozenset({
    "-d", "--max-depth", "-f", "--fasta-ref", "-q", "--min-MQ", "-Q",
    "--min-BQ", "-a", "--annotate", "-m", "--min-ireads", "-L",
    "--max-idepth", "-o", "--output", "-O", "--output-type", "-r",
    "--regions", "-t", "--targets", "-R", "--regions-file", "-T",
    "--targets-file", "--threads", "--no-version", "--samples", "-s",
    "--no-BAQ", "-B", "--full-BAQ", "-D", "-E", "-C", "--adjust-MQ",
    "--max-BQ", "--delta-BQ", "--ignore-RG", "--ls", "--skip-all-set",
    "--ns", "--skip-any-set", "--lu", "--skip-all-unset", "--nu",
    "--skip-any-unset", "-S", "--samples-file", "-b", "--bam-list", "-g",
    "--gvcf", "-v", "--verbosity", "-W", "--write-index", "-X", "--config",
    "-e", "--ext-prob", "-F", "--gap-frac", "-h", "--tandem-qual", "-I",
    "--skip-indels", "-M", "--max-read-len", "-o", "--open-prob", "-p",
    "--per-sample-mF", "-P", "--platforms", "--ar", "--ambig-reads",
    "--indel-bias", "--del-bias", "--score-vs-ref", "--indel-size",
    "--indels-2.0", "--indels-cns", "--seqq-offset", "--no-indels-cns",
    "--poly-mqual", "--seed", "--max-BP", "--no-reference",
})


def _argv(config: PipelineConfig, tmp_path: Path) -> List[str]:
    """The argv the STAGE builds, end to end, through production code."""
    settings = stage.call_settings(config)
    pileup = settings["bcftools"]["mpileup"]
    return adapter.mpileup_command(
        "bcftools",
        bam=tmp_path / "s.bam",
        reference=tmp_path / "ref.fna",
        min_mapq=int(pileup["min_mapq"]),
        min_bq=int(pileup["min_bq"]),
        annotate=str(pileup["annotate"]),
        max_depth=int(pileup["max_depth"]),
        min_ireads=int(pileup["min_ireads"]),
        output_type=str(pileup["output_type"]),
        output_vcf=tmp_path / f"s.{adapter.RAW_VCF}",
    )


# ---------------------------------------------------------------------------
# (a) + (b): the bound is present in the emitted argv, and the pre-fix argv
# does not have it -- so the test genuinely fails before the fix.
# ---------------------------------------------------------------------------


class TestTheBoundIsPresentOnTheCommandLine:
    def test_the_output_encoding_is_bounded_and_named(
        self, config: PipelineConfig, tmp_path: Path
    ):
        """`-O z` is the fix's one argv addition. Asserted, not assumed."""
        argv = _argv(config, tmp_path)
        assert "-O" in argv, argv
        assert argv[argv.index("-O") + 1] == DEFAULT_MPILEUP_OUTPUT_TYPE == "z"

    def test_the_depth_bound_reaches_the_argv_explicitly(
        self, config: PipelineConfig, tmp_path: Path
    ):
        """`-d` was already emitted; this pins that it is emitted, and how.

        The pre-fix failure mode was not a missing `-d` but a depth bound whose
        value was decided in two places at once -- `science.yaml` and a bare
        `250` literal in the adapter. This asserts the number on the command
        line is the loader's one number.
        """
        argv = _argv(config, tmp_path)
        assert argv[argv.index("-d") + 1] == str(DEFAULT_MPILEUP_MAX_DEPTH) == "250"

    def test_the_adapter_no_longer_carries_its_own_depth_literal(self):
        """The hidden literal is gone; the default is named and imported.

        `pileup_cfg.get("max_depth", 250)` in `call_isolate` is what made the
        bound invisible: change `science.yaml` and the adapter still had a
        second, undiscoverable number. Asserted against the SOURCE, because the
        behaviour it guards is the absence of a literal.
        """
        source = Path(adapter.__file__).read_text(encoding="utf-8")
        assert 'pileup_cfg.get("max_depth", 250)' not in source
        assert "DEFAULT_MPILEUP_MAX_DEPTH" in source

    def test_the_raw_vcf_name_carries_the_gz_extension(self):
        """`-O z` only compresses when the FILENAME ends `.gz`.

        Measured: `-O z -o mismatch.vcf` wrote 69,488,493 bytes of UNCOMPRESSED
        data and exited 0 with no warning, against 1,861,798 bytes for the same
        command ending `.gz`. So the extension is half the fix, and a test that
        only checked the flag would pass against a no-op.
        """
        assert adapter.RAW_VCF.endswith(".gz")


# ---------------------------------------------------------------------------
# (c): R14 guard. `min_ireads` stays 1. Pinned so no future change to this
# file can quietly remove or raise it.
# ---------------------------------------------------------------------------


class TestR14MinIreadsIsStillOne:
    def test_min_ireads_is_still_1_on_the_argv(
        self, config: PipelineConfig, tmp_path: Path
    ):
        argv = _argv(config, tmp_path)
        assert argv[argv.index("-m") + 1] == str(DEFAULT_MPILEUP_MIN_IREADS) == "1"

    def test_min_ireads_is_unrelated_to_the_fix(self, tmp_path: Path):
        """`-m` must not be confused with the depth cap.

        This is the exact confusion the R14 docstring warns about: `-m` is
        `--min-ireads` (gapped reads supporting an indel), NOT a depth filter, and
        `-d` is the depth knob. So raising `-m` "to bound memory" would be
        wrong, and removing it would reinstate the pinned default of 2, which no
        assembly input can satisfy.
        """
        argv = adapter.mpileup_command(
            "bcftools", bam=tmp_path / "b.bam", reference=tmp_path / "r.fna",
            min_mapq=20, min_bq=20, annotate="FORMAT/DP", max_depth=250,
        )
        assert argv[argv.index("-m") + 1] == "1"
        assert argv[argv.index("-d") + 1] == "250"
        assert argv[argv.index("-m") + 1] != argv[argv.index("-d") + 1]

    def test_the_loader_default_is_still_one(self):
        assert DEFAULT_MPILEUP_MIN_IREADS == 1


# ---------------------------------------------------------------------------
# (d): no invented flags. Checked against the installed binary's own usage.
# ---------------------------------------------------------------------------


class TestNoFlagIsInvented:
    @pytest.mark.skipif(
        not shutil.which("bcftools"),
        reason="bcftools is not installed on this machine",
    )
    def test_the_usage_block_really_says_what_this_file_quotes(self):
        """Non-vacuity for the literals above.

        If the installed bcftools were not 1.23.1, the quoted usage lines would
        be decoration. This requires them to appear in the binary's own output.
        """
        out = subprocess.run(
            ["bcftools", "mpileup"], capture_output=True, text=True, check=False
        )
        usage = out.stdout + out.stderr
        assert MPILEUP_DEPTH_USAGE in usage
        assert MPILEUP_OUTPUT_TYPE_USAGE in usage

    @pytest.mark.skipif(
        not shutil.which("bcftools"),
        reason="bcftools is not installed on this machine",
    )
    def test_max_bp_does_not_exist_and_is_not_emitted(self, tmp_path: Path):
        """The flag the defect report reached for is absent from 1.23.1.

        `-b` IS a real flag but it is `--bam-list FILE`, nothing to do with
        depth, so emitting `-b <n>` would be worse than emitting nothing.
        """
        result = subprocess.run(
            ["bcftools", "mpileup", "--max-BP", "1000", str(tmp_path / "x.bam")],
            capture_output=True, text=True, check=False,
        )
        assert "unrecognized" in (result.stdout + result.stderr).lower()
        assert "-b" not in _argv(_config(), tmp_path)

    def test_every_emitted_flag_is_one_the_pinned_bcftools_documents(
        self, config: PipelineConfig, tmp_path: Path
    ):
        """No emitted token may be a flag 1.23.1 does not have.

        Positionals (the reference, the output path, the BAM) are not flags and
        are skipped; `-a`'s value contains no leading dash and is skipped too.
        """
        argv = _argv(config, tmp_path)
        value_of = {
            argv[i + 1] for i, token in enumerate(argv[:-1])
            if token.startswith("-") and token in {"-f", "-q", "-Q", "-d", "-m",
                                                   "-a", "-o", "-O"}
        }
        for token in argv:
            if not token.startswith("-"):
                continue
            assert token not in value_of, f"value mistaken for a flag: {token}"
            assert token in MPILEUP_FLAGS_THAT_EXIST, f"invented flag: {token}"


# ---------------------------------------------------------------------------
# (e): the callable variant set, MEASURED rather than asserted from the source.
# ---------------------------------------------------------------------------


#: The edits the synthetic assembly carries against its reference. Written as
#: (reference_start_0based, ref_allele, alt_allele) with ``None`` for a pure
#: insertion, so the CIGAR is DERIVED from them rather than asserted.
#:
#: A CIGAR is the whole point here. An earlier version of this fixture wrote a
#: single ``<len>M`` CIGAR and then asserted an indel had been called; mpileup
#: read that as a perfect match with no gaps, so there was no indel candidate
#: and the assertion failed. The deletion has to be IN THE CIGAR to be visible
#: to `--min-ireads`, which is the flag under test.
#: What each edit DOES, as (reference_start, ref_allele_length,
#: substitute_with_or_None_for_a_deletion). The reference bases themselves are
#: read from the generated reference rather than written down here, because a
#: hand-written allele has to be kept in step with the seed and drifting out of
#: step makes the fixture fail for a reason that has nothing to do with stage 6.
_SYNTHETIC_EDITS = (
    (2_000, 1, "C"),
    (5_000, 1, "G"),
    (9_001, 1, "A"),
    (14_500, 1, "C"),
    (17_500, 1, "G"),
    (18_000, 60, None),   # a 60 bp deletion: the indel that needs `-m 1`
)
_SYNTHETIC_REF_LEN = 20_000


def _write_synthetic_pileup_fixture(tmp_path: Path) -> Path:
    """A small reference plus one depth-1 assembly alignment over it.

    Built by hand rather than by running minimap2, so the test needs no aligner
    and the geometry is exactly the one the pipeline feeds mpileup: ONE record
    whose SEQ is a whole contig.
    """
    import random

    random.seed(20260905)
    ref = "".join(random.choice("ACGT") for _ in range(_SYNTHETIC_REF_LEN))
    ref_path = tmp_path / "ref.fna"
    ref_path.write_text(
        ">NC_SYN.1\n"
        + "\n".join(ref[i:i + 70] for i in range(0, len(ref), 70))
        + "\n",
        encoding="utf-8",
    )
    # `bcftools index` refuses a plain FASTA ("in a format that cannot be
    # usefully indexed", exit 255); `samtools faidx` is what writes the .fai
    # mpileup's `-f` needs. Same htslib, different subcommand.
    subprocess.run(
        ["samtools", "faidx", str(ref_path)],
        capture_output=True, text=True, check=True,
    )
    bam_path = tmp_path / "s.sorted.bam"
    _write_single_record_bam(bam_path, ref)
    return bam_path


def _query_and_cigar(ref: str) -> "tuple[str, str]":
    """Apply :data:`_SYNTHETIC_EDITS` to ``ref``, returning (query, CIGAR)."""
    query: list = []
    cigar: list = []
    cursor = 0

    def emit_match(run: int) -> None:
        if run:
            query.append(ref[cursor:cursor + run])
            cigar.append(f"{run}M")

    for start, length, alt_base in _SYNTHETIC_EDITS:
        emit_match(start - cursor)
        if alt_base is None:                         # deletion
            cigar.append(f"{length}D")
        else:                                        # substitution
            query.append(ref[start] if alt_base == ref[start]
                         else alt_base)
            cigar.append("1M")
        cursor = start + length
    emit_match(len(ref) - cursor)
    return "".join(query), "".join(cigar)


def _write_single_record_bam(path: Path, ref: str) -> None:
    """One whole-contig alignment record -- the shape `minimap2 -a` emits."""
    import pysam

    query, cigar = _query_and_cigar(ref)
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"LN": len(ref), "SN": "NC_SYN.1"}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        record = pysam.AlignedSegment()
        record.query_name = "contig_1"
        record.query_sequence = query
        record.flag = 0
        record.reference_id = 0
        record.reference_start = 0
        record.mapping_quality = 60
        record.cigarstring = cigar
        record.query_qualities = pysam.qualitystring_to_array("I" * len(query))
        record.next_reference_id = -1
        record.next_reference_start = -1
        record.template_length = 0
        out.write(record)
    pysam.index(str(path))


@pytest.mark.skipif(
    not (shutil.which("bcftools") and shutil.which("samtools")),
    reason="bcftools or samtools is not installed on this machine",
)
class TestTheCallableVariantSetIsUnchanged:
    """(e). Pre-fix argv vs post-fix argv, real bcftools, compared by OUTPUT.

    "Preserved" is asserted the only way worth trusting: not by reading the two
    commands and reasoning about them, but by requiring that the `call -mv`
    records are identical. A change to compression that moved one DP, one AD or
    one allele would show up here.
    """

    def _calls(self, tmp_path: Path, *, output_type: str, raw_name: str) -> str:
        bam = _write_synthetic_pileup_fixture(tmp_path)
        raw = tmp_path / raw_name
        called = tmp_path / "called.vcf"
        argv = adapter.mpileup_command(
            "bcftools", bam=bam, reference=tmp_path / "ref.fna",
            min_mapq=20, min_bq=20, annotate="FORMAT/DP,FORMAT/AD",
            max_depth=DEFAULT_MPILEUP_MAX_DEPTH,
            min_ireads=DEFAULT_MPILEUP_MIN_IREADS,
            output_type=output_type, output_vcf=raw,
        )
        pileup = subprocess.run(argv, capture_output=True, text=True, check=False)
        assert pileup.returncode == 0, pileup.stderr
        assert "unrecognized" not in pileup.stderr.lower()
        call = adapter.call_command(
            "bcftools", input_vcf=raw, output_vcf=called,
            variants_only=True, multiallelic="split",
        )
        called_run = subprocess.run(
            call, capture_output=True, text=True, check=False
        )
        assert called_run.returncode == 0, called_run.stderr
        # Drop the header entirely: it records the argv, which is the one thing
        # that is SUPPOSED to differ. What must match is the called records.
        return "\n".join(
            line for line in called.read_text(encoding="utf-8").splitlines()
            if not line.startswith("#")
        )

    def test_pre_fix_and_post_fix_argv_call_the_same_variants(
        self, config: PipelineConfig, tmp_path: Path
    ):
        before = self._calls(
            tmp_path, output_type="v", raw_name="pre.raw.vcf"
        )
        after = self._calls(
            tmp_path, output_type="z", raw_name=f"post.{adapter.RAW_VCF}"
        )
        assert before, "the fixture must call something, or the comparison is vacuous"
        assert before == after

    def test_the_fixture_really_exercises_snvs_and_indels(
        self, tmp_path: Path
    ):
        """Non-vacuity for the equality assertion above.

        Comparing two empty call sets would pass while proving nothing, so this
        requires the fixture to produce calls, and requires at least one of them
        to be an INDEL -- which only happens because `min_ireads` is 1. If a
        future change raised `-m`, this fails and says why.
        """
        calls = self._calls(
            tmp_path, output_type="z", raw_name=f"f.{adapter.RAW_VCF}"
        )
        rows = [row.split("\t") for row in calls.splitlines()]
        assert len(rows) >= 5, f"fixture called almost nothing: {calls}"
        # An indel is a row whose REF or ALT is longer than one base -- NOT
        # necessarily ALT: bcftools writes a deletion the other way round, as a
        # 60 bp REF against a 1 bp ALT. Checking ALT alone would pass on a
        # fixture containing no indels at all.
        assert any(len(r[3]) > 1 or len(r[4]) > 1 for r in rows), (
            "no indel was called, so this fixture is not exercising min_ireads: "
            f"{[(r[1], r[3], r[4]) for r in rows]}"
        )

    def test_the_indel_is_present_only_because_min_ireads_is_1(
        self, tmp_path: Path
    ):
        """R14's payload, measured: `-m 2` loses the indel this fix must not lose.

        The pinned bcftools default of 2 is a gapped-reads threshold an assembly
        can never satisfy -- there is exactly ONE gapped read per indel. This
        pair is the assertion that `-m 1` is load-bearing rather than decorative,
        run on this fixture's own data rather than quoted from the docstring.
        """
        def indels(min_ireads: str) -> int:
            bam = _write_synthetic_pileup_fixture(tmp_path)
            raw = tmp_path / f"m{min_ireads}.{adapter.RAW_VCF}"
            called = tmp_path / f"m{min_ireads}.calls.vcf"
            subprocess.run(
                adapter.mpileup_command(
                    "bcftools", bam=bam, reference=tmp_path / "ref.fna",
                    min_mapq=20, min_bq=20, annotate="FORMAT/DP,FORMAT/AD",
                    max_depth=DEFAULT_MPILEUP_MAX_DEPTH,
                    min_ireads=int(min_ireads), output_vcf=raw,
                ),
                capture_output=True, text=True, check=True,
            )
            subprocess.run(
                adapter.call_command(
                    "bcftools", input_vcf=raw, output_vcf=called,
                    variants_only=True, multiallelic="split",
                ),
                capture_output=True, text=True, check=True,
            )
            return sum(
                1 for line in called.read_text(encoding="utf-8").splitlines()
                if not line.startswith("#")
                and (len(line.split("\t")[3]) > 1 or len(line.split("\t")[4]) > 1)
            )

        assert indels("1") >= 1, "the fixture's indel vanished even at -m 1"
        assert indels("2") == 0, (
            "the fixture no longer distinguishes -m 1 from -m 2, so it cannot "
            "demonstrate that R14's threshold matters"
        )

    def test_compression_actually_happens(self, tmp_path: Path):
        """The fix is not a no-op: the compressed file is materially smaller.

        Without this, a build that dropped `-O z` would still pass the
        equality test above, because both sides would then be uncompressed and
        still equal. This is what makes (e) mean something.
        """
        import gzip

        bam = _write_synthetic_pileup_fixture(tmp_path)
        plain = tmp_path / "plain.vcf"
        packed = tmp_path / f"packed.{adapter.RAW_VCF}"
        common = dict(
            bam=bam, reference=tmp_path / "ref.fna", min_mapq=20, min_bq=20,
            annotate="FORMAT/DP,FORMAT/AD", max_depth=DEFAULT_MPILEUP_MAX_DEPTH,
            min_ireads=DEFAULT_MPILEUP_MIN_IREADS,
        )
        subprocess.run(
            adapter.mpileup_command(
                "bcftools", output_type="v", output_vcf=plain, **common
            ),
            capture_output=True, text=True, check=True,
        )
        subprocess.run(
            adapter.mpileup_command(
                "bcftools", output_type="z", output_vcf=packed, **common
            ),
            capture_output=True, text=True, check=True,
        )
        assert packed.stat().st_size < plain.stat().st_size / 2
        # And the bytes are the same records: bcftools call accepts the .gz.
        assert gzip.open(packed, "rt", encoding="utf-8").readline().startswith("##")


# ---------------------------------------------------------------------------
# N7: the config key, its consumer and these tests land together. Assert the
# whole chain rather than any one end of it.
# ---------------------------------------------------------------------------


class TestTheKeyAndItsConsumerAreWired:
    def test_science_yaml_states_the_key(self, config: PipelineConfig):
        pileup = config.raw["variants"]["bcftools"]["mpileup"]
        assert pileup["output_type"] == "z"
        assert pileup["max_depth"] == 250

    def test_the_stage_resolves_both_keys_from_the_loader(
        self, config: PipelineConfig
    ):
        settings = stage.call_settings(config)
        pileup = settings["bcftools"]["mpileup"]
        assert pileup["max_depth"] == DEFAULT_MPILEUP_MAX_DEPTH
        assert pileup["output_type"] == DEFAULT_MPILEUP_OUTPUT_TYPE

    def test_a_config_without_the_keys_still_resolves(
        self, config: PipelineConfig, tmp_path: Path
    ):
        """A missing section yields the DEFAULT, not a KeyError.

        The same defensive rule `variants_mpileup_min_ireads` follows: an
        annotation-less config is a legitimate STUB/TEST run.
        """
        import copy

        raw = copy.deepcopy(dict(config.raw))
        raw.get("variants", {}).get("bcftools", {}).pop("mpileup", None)
        stripped = PipelineConfig(
            root=config.root, raw=raw, antibiotics=config.antibiotics,
            allowed_phenotypes=config.allowed_phenotypes,
            analysis=dict(config.analysis), mechanisms=dict(config.mechanisms),
            regulators=dict(config.regulators),
            references=dict(config.references),
            antibiotic_specs=dict(config.antibiotic_specs),
            organism=dict(config.organism),
            qc=config.qc, gwas=config.gwas, phylogeny=config.phylogeny,
            convergence=config.convergence, paths=dict(config.paths),
            runtime=dict(config.runtime),
        )
        settings = stage.call_settings(stripped)
        pileup = settings["bcftools"]["mpileup"]
        assert pileup["max_depth"] == DEFAULT_MPILEUP_MAX_DEPTH
        assert pileup["output_type"] == DEFAULT_MPILEUP_OUTPUT_TYPE
        assert pileup["min_ireads"] == DEFAULT_MPILEUP_MIN_IREADS

    def test_call_settings_does_not_mutate_the_shared_config(
        self, config: PipelineConfig
    ):
        """`config.raw` is shared by every later stage in the run."""
        before = dict(config.raw["variants"]["bcftools"]["mpileup"])
        stage.call_settings(config)
        assert config.raw["variants"]["bcftools"]["mpileup"] == before


def _config() -> PipelineConfig:
    from papipeline.config import load_config

    return load_config()


# ---------------------------------------------------------------------------
# The INJECTED FAKE RUNNER (R12): stage 6 driven end to end with no
# bioinformatics tool installed and no fixture on disk.
# ---------------------------------------------------------------------------


class TestTheInjectedFakeBcftoolsRunnerRefusesAnUnboundedCommand:
    """The argv is checked by a runner that enforces the bound, not by a string.

    A fake `bcftools` first on `PATH` writes a minimal VCF at whatever `-o`
    names and refuses (exit 9, the SIGKILL shape) when the argv carries no
    compression bound. This drives the REAL adapter and the REAL stage path, so
    it fails against the pre-fix command and passes against this one.
    """

    @pytest.fixture
    def fake_bcftools(self, tmp_path: Path, monkeypatch):
        log = tmp_path / "argv.log"
        script = tmp_path / "bcftools"
        script.write_text(
            "#!/bin/sh\n"
            f'printf "%s\\n" "$*" >> "{log}"\n'
            # SIGKILL's exit status: what the round-12 run saw.
            'if [ "$1" = "mpileup" ]; then\n'
            '  bounded=no\n'
            '  prev=\n'
            '  for tok in "$@"; do\n'
            '    if [ "$prev" = "-O" ]; then bounded=yes; fi\n'
            '    prev=$tok\n'
            '  done\n'
            '  if [ "$bounded" = no ]; then\n'
            '    echo "[mpileup] unbounded output: no -O" >&2\n'
            "    exit 9\n"
            "  fi\n"
            '  out=\n'
            '  prev=\n'
            '  for tok in "$@"; do\n'
            '    if [ "$prev" = "-o" ]; then out=$tok; fi\n'
            '    prev=$tok\n'
            "  done\n"
            '  printf "##fileformat=VCFv4.2\\n" > "$out"\n'
            '  printf "##contig=<ID=NC_SYN.1,length=20000>\\n" >> "$out"\n'
            '  printf "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\tS1\\n" >> "$out"\n'
            '  printf "NC_SYN.1\\t2000\\t.\\tA\\tC\\t30.4183\\t.\\tDP=1\\tGT:AD\\t0/1:1,0\\n" >> "$out"\n'
            "  exit 0\n"
            "fi\n"
            "if [ \"$1\" = \"call\" ]; then\n"
            '  out=\n'
            '  prev=\n'
            '  for tok in "$@"; do\n'
            '    if [ "$prev" = "-o" ]; then out=$tok; fi\n'
            '    prev=$tok\n'
            "  done\n"
            '  printf "##fileformat=VCFv4.2\\n" > "$out"\n'
            '  printf "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\tS1\\n" >> "$out"\n'
            '  printf "NC_SYN.1\\t2000\\t.\\tA\\tC\\t30.4183\\t.\\tDP=1\\tGT:AD\\t0/1:1,0\\n" >> "$out"\n'
            "  exit 0\n"
            "fi\n"
            "exit 1\n",
            encoding="utf-8",
        )
        script.chmod(0o755)
        monkeypatch.setenv("PATH", f"{tmp_path}{':'}{__import__('os').environ['PATH']}")
        return log

    def test_the_emitted_argv_carries_the_bound_the_runner_enforces(
        self, config: PipelineConfig, tmp_path: Path, fake_bcftools: Path
    ):
        argv = _argv(config, tmp_path)
        result = subprocess.run(argv, capture_output=True, text=True, check=False)
        assert result.returncode != 9, (
            "the fake bcftools refused the argv: no -O was present, so the "
            "command this fix is supposed to change is still unbounded"
        )
        assert result.returncode == 0, result.stderr
        assert "-O z" in fake_bcftools.read_text(encoding="utf-8")

    def test_the_runner_really_would_refuse_an_unbounded_argv(
        self, config: PipelineConfig, tmp_path: Path, fake_bcftools: Path
    ):
        """Non-vacuity: the guard above is enforced by the runner, so prove the
        runner refuses. Without this, a fake that always exits 0 would make the
        previous test meaningless."""
        # The pre-fix shape: `-O` absent entirely, which is what `call_isolate`
        # emitted at f44fde8. Built by DELETING `-O` from the real argv rather
        # than by hand-writing a command, so this really is "the current
        # unbounded command" and cannot drift away from what the adapter emits.
        full = adapter.mpileup_command(
            "bcftools", bam=tmp_path / "b.bam", reference=tmp_path / "r.fna",
            min_mapq=20, min_bq=20, annotate="FORMAT/DP",
            max_depth=DEFAULT_MPILEUP_MAX_DEPTH,
            output_vcf=tmp_path / "u.vcf",
        )
        unbounded = [
            token for i, token in enumerate(full)
            if token != "-O" and full[i - 1] != "-O"
        ]
        assert "-O" not in unbounded
        result = subprocess.run(
            unbounded, capture_output=True, text=True, check=False
        )
        assert result.returncode == 9, result.stderr
        assert "unbounded" in result.stderr