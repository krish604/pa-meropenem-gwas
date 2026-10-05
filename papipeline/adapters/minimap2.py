"""The stage-6 adapter: one assembly in, one VCF of called variants out.

`config/science.yaml` already pins this invocation, against a measured REAL run
rather than an expectation of what the tools do (see the `variants:` block and
the measurements recorded beside it). So this module is deliberately thin: it
reads that configuration, assembles argument vectors from it, and refuses when a
tool or a reference index is missing.

Three decisions are load-bearing enough to state.

**The reference index is verified before anything reads it, here.** A `.fai`
built by indexing a *bgzipped* copy of PAO1 has been observed reporting 64,400
bases for a 6,264,404-base file. Every tool then read a reference two orders of
magnitude too small and produced confident, entirely wrong output — 64
alignments for a 6.4 Mb genome, and thousands of "variants" in a 1 Mb window,
all sharing one constant QUAL. Nothing errored, which is the whole danger.
:func:`verify_reference_index` is therefore called by :func:`call_isolate` and
cannot be skipped by a caller who forgets; it delegates to
:func:`papipeline.reference.verify_fai_index` so there is one implementation of
the check rather than two that can disagree.

**SAM text, not a samtools conversion.** The pinned samtools in this environment
(0.1.19) auto-detects input format by extension, refuses SAM text outright, and
has no `-O` flag; it only indexes FASTA. So minimap2 emits SAM and pysam does
the sort and index, which bcftools requires anyway since it refuses unsorted or
unindexed input. Both facts are recorded in `science.yaml` under
`variants.samtools`, and both are read from there rather than restated.

**No interpretation.** This module produces a VCF path. What the calls *mean* —
the columns, the INFO fields carried and dropped, the multiallelic expansion —
belongs to :mod:`papipeline.stages.variants`, which owns the schema. Keeping the
seam at the file is what lets the stage be tested against committed fixtures
with no aligner installed.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

from ..config.loader import (
    DEFAULT_MPILEUP_MAX_DEPTH,
    DEFAULT_MPILEUP_MIN_IREADS,
    DEFAULT_MPILEUP_OUTPUT_TYPE,
)
from ..errors import PipelineError, ToolExecutionError, ToolNotAvailableError
from ..logging_utils import get_logger
from ..reference import verify_fai_index

LOGGER = get_logger("adapters.minimap2")

#: The two tools this adapter drives. Named so a refusal can say which is absent.
CALLER_TOOLS: Sequence[str] = ("minimap2", "bcftools")

#: VCF suffix for the mpileup output. Kept as a constant because the `call`
#: step's input is a different file from its output and mixing them up produces
#: a VCF that parses and is empty.
#:
#: The `.gz` is load-bearing, not cosmetic. `-O z` is emitted with
#: :data:`DEFAULT_MPILEUP_OUTPUT_TYPE`, and bcftools decides compression from the
#: OUTPUT FILENAME's extension, not from the flag alone: with `-O z` and a name
#: ending in plain `.vcf` it wrote 69,488,493 bytes of UNCOMPRESSED data and
#: exited 0 without a warning, against 1,861,798 bytes for the same command
#: ending `.gz`. So renaming the file without setting the flag would have done
#: nothing, and setting the flag without renaming it would have done nothing
#: while looking like a fix.
RAW_VCF = "raw.vcf.gz"
CALLED_VCF = "calls.vcf"


def resolve_executable(name: str, override: Optional[str] = None) -> str:
    """Locate one executable, or refuse by name.

    The environment-variable override exists because these tools are usually
    installed inside a conda/micromamba environment that is not on the login
    ``PATH``; hardcoding a personal absolute path would make the repository
    unrunnable anywhere else.
    """
    if override:
        candidate = Path(override).expanduser()
        if not candidate.is_file():
            raise ToolNotAvailableError(
                f"{name} is configured at a path that is not a file",
                tool=name, configured=str(override),
            )
        return str(candidate)
    found = shutil.which(name)
    if not found:
        raise ToolNotAvailableError(
            f"{name} was not found on PATH. REAL-mode variant calling runs the "
            "real tools; it is never substituted or simulated. Put it on PATH, "
            "or set "
            f"{name.upper()}_BIN to its full path.",
            tool=name,
            hint=f"Install the pinned environment declared in science.yaml under "
                 f"variants.bcftools / variants.minimap2, or set {name.upper()}_BIN.",
        )
    return found


def require_callers(
    statuses: Mapping[str, Optional[str]],
) -> Dict[str, str]:
    """Every tool this adapter needs, or a refusal naming the absent one.

    Args:
        statuses: Tool name -> executable path, ``None`` where unavailable.

    Returns:
        The same mapping, which is the resolved executable per tool.

    Raises:
        ToolNotAvailableError: Any required tool is missing. The message names
            it, because a refusal naming no tool leaves the operator guessing.
    """
    resolved: Dict[str, str] = {}
    for name in CALLER_TOOLS:
        executable = statuses.get(name)
        if not executable:
            raise ToolNotAvailableError(
                f"Stage 6 needs {name}, which is not installed. Naming the one "
                "that is absent, so the fix is unambiguous.",
                tool=name,
                required=",".join(CALLER_TOOLS),
            )
        resolved[name] = str(executable)
    return resolved


def align_command(
    executable: str,
    *,
    preset: str,
    assembly: Path,
    reference: Path,
    output_sam: Path,
) -> List[str]:
    """``minimap2 -x <preset> -a --secondary=yes <reference> <assembly>``.

    **Order is load-bearing.** minimap2's usage is
    ``minimap2 [options] <target.fa> [query.fa]``: the **reference is the
    target** and the **assembly is the query**.

    Swapped, minimap2 aligns PAO1 against the isolate, and the result is
    silently wrong rather than wrong-looking. Observed here: every record came
    out named ``NC_002516.2`` with a SEQ of 6,264,404 bases - the *reference's*
    length - under an ``@SQ`` built from the assembly's 24 contigs. Nothing
    errored. ``bcftools mpileup`` then piled up a query it believed was the
    reference and ``call -v`` returned **zero variants** for a genome with tens
    of thousands of them, and the cohort merge would have reported every real
    difference as absent.

    ``-a`` because the conversion path is pysam's, which reads SAM text; see
    the module docstring. ``--secondary=yes`` because
    ``config/science.yaml`` records that ``--secondary=no`` was measurably wrong
    on this cohort.
    """
    return [
        str(executable),
        "-x", str(preset),
        "-a",
        "--secondary=yes",
        str(reference),
        str(assembly),
    ]


def mpileup_command(
    executable: str,
    *,
    bam: Path,
    reference: Path,
    min_mapq: int,
    min_bq: int,
    annotate: str,
    max_depth: int = DEFAULT_MPILEUP_MAX_DEPTH,
    min_ireads: int = DEFAULT_MPILEUP_MIN_IREADS,
    output_type: str = DEFAULT_MPILEUP_OUTPUT_TYPE,
    output_vcf: Optional[Path] = None,
) -> List[str]:
    """``bcftools mpileup -f <reference> -q -Q -a -d -m -O [-o] <bam>``.

    ``-f`` is not optional: without it bcftools has no reference, and the calls
    it produces would be about nothing in particular. ``annotate`` is passed as
    one argument, because ``FORMAT/DP,FORMAT/AD`` is a single value and
    splitting it would invent a flag that does not exist.

    ``-m`` is ``--min-ireads``, verified against the pinned bcftools 1.23.1
    (``bcftools mpileup`` usage line 53) as "Minimum number gapped reads for
    indel candidates [2]". It is NOT a read-depth filter -- ``-d`` above is the
    depth knob. It is emitted unconditionally because a threshold of 2 cannot be
    satisfied by an assembly (one gapped read per indel), and omitting the flag
    would silently reinstate that unusable default. See
    :data:`papipeline.config.loader.DEFAULT_MPILEUP_MIN_IREADS`.

    ``-d`` is bcftools' only memory lever and it is bounded here EXPLICITLY
    rather than left to the tool's default, so the value is visible in the argv
    and in `run_manifest.json`. It is bcftools' own inherited default of 250,
    not a threshold this project chose; see
    :data:`papipeline.config.loader.DEFAULT_MPILEUP_MAX_DEPTH` for the measured
    RSS curve and for why lowering it is a science decision about the callable
    variant set rather than an engineering one.

    ``-O`` is emitted unconditionally and defaults to ``'z'``. Both bcftools
    AND this function's caller rely on the two agreeing, and the failure when
    they disagree is silent: ``-O z`` writing to a filename that does not end
    ``.gz`` produces uncompressed data, exit code 0, no warning. That is why
    :data:`RAW_VCF` ends in ``.vcf.gz``. ``bcftools call`` reads it either way
    (verified), so this is safe, and it does not change what is called.
    """
    command = [
        str(executable), "mpileup",
        "-f", str(reference),
        "-q", str(int(min_mapq)),
        "-Q", str(int(min_bq)),
        "-a", str(annotate),
        "-d", str(int(max_depth)),
        "-m", str(int(min_ireads)),
        "-O", str(output_type),
    ]
    if output_vcf is not None:
        command += ["-o", str(output_vcf)]
    command.append(str(bam))
    return command


def call_command(
    executable: str,
    *,
    input_vcf: Path,
    output_vcf: Path,
    variants_only: bool,
    multiallelic: Any,
    threads: Optional[int] = None,
) -> List[str]:
    """``bcftools call -m -v -o <out> <in>``.

    ``-m`` is ``--multiallelic-caller`` and ``-v`` is ``--variants-only``; they
    are different switches and only together do they mean what the pinned
    ``call -mv`` run meant. Emitting ``-m`` alone would emit every site,
    including reference matches, and the merge would treat them as
    polymorphisms.
    """
    command = [str(executable), "call"]
    if multiallelic in ("split", True, "multiallelic"):
        command.append("-m")
    if variants_only:
        command.append("-v")
    if threads is not None:
        command += ["--threads", str(int(threads))]
    command += ["-o", str(output_vcf), str(input_vcf)]
    return command


def verify_reference_index(fai_path: Path, *, expected_bases: int) -> int:
    """Check the reference's ``.fai`` against the configured base count.

    A thin, named delegation to :func:`papipeline.reference.verify_fai_index`, so
    the adapter's path cannot drift from the guard that integration tests prove
    catches a truncating index.
    """
    return verify_fai_index(Path(fai_path), expected_bases=int(expected_bases))


@dataclass(frozen=True)
class CallResult:
    """What one isolate's calling produced."""

    sample_id: str
    vcf: Path
    commands: Sequence[Sequence[str]] = field(default_factory=tuple)


def _run(command: Sequence[str], *, cwd: Optional[Path] = None) -> None:
    """Run one tool invocation, refusing loudly on a non-zero exit.

    A non-zero exit from any of these is a real failure, not a warning: a
    truncated alignment still produces a parseable VCF, and an empty one looks
    exactly like an isolate that genuinely carries no variant.
    """
    LOGGER.info("exec: %s", " ".join(command))
    try:
        completed = subprocess.run(
            list(command), cwd=str(cwd) if cwd else None,
            capture_output=True, text=True, check=False,
        )
    except OSError as exc:
        raise ToolExecutionError(
            "Failed to launch external tool",
            command=" ".join(command), error=str(exc),
        ) from exc
    if completed.returncode != 0:
        raise ToolExecutionError(
            "External tool exited non-zero",
            command=" ".join(command),
            returncode=completed.returncode,
            stderr=(completed.stderr or "")[-2000:],
        )


def convert_and_index(
    sam_path: Path, bam_path: Path, *, reference: Optional[Path] = None
) -> Path:
    """Sort and index SAM text into a BAM bcftools will accept.

    pysam does this because the pinned samtools cannot read SAM text at all
    (``config/science.yaml`` -> ``variants.samtools``). bcftools refuses an
    unsorted or unindexed input, so neither step is optional.

    Args:
        sam_path: The SAM minimap2 wrote.
        bam_path: Where to write the sorted, indexed BAM.
        reference: The reference FASTA, whose sequence dictionary is installed
            on the BAM. Optional but strongly preferred - see below.

    **On the reference header.** minimap2 writes ``@SQ`` lines for the
    *reference* (its target) while the records are named for the *query* (the
    assembly's contigs). Left alone, that mismatch means htslib resolves records
    against a dictionary that does not contain them and demotes them to
    unmapped; observed as 7 of 94 reads dropped on a real assembly, with
    ``unrecognized reference name`` on stderr.

    With the argument order correct, bcftools still reads the un-rewritten BAM
    and calls correctly, so this rewrite is **defensive rather than currently
    load-bearing**. It is kept because the cost is one small header rewrite and
    the failure it prevents is invisible: reads quietly demoted, a
    correspondingly reduced call set, and nothing in the output to say so. It
    was demonstrably load-bearing while the arguments to :func:`align_command`
    were swapped, where the same mismatch dropped *every* read and the stage
    reported zero variants.
    """
    try:
        import pysam
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise PipelineError(
            "pysam is required to convert SAM to a sorted, indexed BAM, and it "
            "is not importable. The pinned samtools cannot read SAM text, so "
            "there is no fallback path.",
            hint="Install pysam in the analysis environment.",
        ) from exc

    if reference is None:
        pysam.sort("-o", str(bam_path), str(sam_path))
        pysam.index(str(bam_path))
        return bam_path

    header = _reference_header(Path(reference))
    with pysam.AlignmentFile(str(sam_path), "r", header=header) as source:
        unsorted = Path(bam_path).with_name(f"{Path(bam_path).stem}.unsorted.bam")
        with pysam.AlignmentFile(str(unsorted), "wb", header=header) as sink:
            for read in source.fetch(until_eof=True):
                sink.write(read)

    pysam.sort("-o", str(bam_path), str(unsorted))
    pysam.index(str(bam_path))
    return bam_path


def reference_dictionary(reference: Path) -> List[Tuple[str, int]]:
    """``(contig, length)`` for the reference, read from its own ``.fai``.

    Read from the index rather than by scanning the FASTA, so the lengths are
    the ones bcftools itself resolves and there is no second parser free to
    disagree with it.
    """
    fai = Path(f"{reference}.fai")
    if not fai.is_file():
        raise PipelineError(
            f"Reference index not found at {fai}. The alignment's sequence "
            "dictionary has to come from the reference's own index; without it "
            "the header keeps naming the assembly's contigs and bcftools calls "
            "nothing.",
            path=str(fai),
        )

    references: List[Tuple[str, int]] = []
    for line in fai.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        references.append((parts[0], int(parts[1])))
    if not references:
        raise PipelineError(
            f"Reference index at {fai} lists no contigs, so no sequence "
            "dictionary can be taken from it.",
            path=str(fai),
        )
    return references


def _reference_header(reference: Path) -> Any:
    """An htslib header naming the reference, for parsing SAM against.

    Supplied to pysam at parse time rather than patched into the SAM text
    afterwards. That distinction is not cosmetic: htslib validates RNAME against
    the header *as it parses*, so a SAM already carrying a mismatched ``@SQ`` is
    read with the affected records demoted to unmapped, and a text-level rewrite
    applied afterwards has already lost them.

    ``@HD`` sort order is left at ``unsorted``: the sort happens next, and the
    sort order in the header is what a reader trusts.
    """
    references = reference_dictionary(reference)
    return {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": name, "LN": length} for name, length in references],
    }


def call_isolate(
    sample_id: str,
    *,
    assembly: Path,
    reference: Path,
    workdir: Path,
    settings: Mapping[str, Any],
    tools: Mapping[str, str],
    threads: int,
) -> CallResult:
    """Align one assembly to the reference and call its variants.

    Args:
        sample_id: The isolate. Only used to name files and to label the result.
        assembly: This isolate's assembly FASTA.
        reference: The verified PAO1 FASTA.
        workdir: Scratch directory for this isolate's SAM and BAM.
        settings: The ``config.raw['variants']`` block. Read here rather than
            restated, so the code cannot disagree with the pinned science.
        tools: Executable per tool name, from :func:`require_callers`.
        threads: From the machine overlay, never from here.

    Returns:
        A :class:`CallResult` pointing at the called VCF.

    Raises:
        DataContractError: The reference index is missing or does not match the
            configured base count. This happens *before* minimap2 runs, because
            a short index produces confident wrong output rather than an error.
        ToolExecutionError: Any invocation exited non-zero.
    """
    minimap2_cfg = settings.get("minimap2") or {}
    samtools_cfg = settings.get("samtools") or {}
    bcftools_cfg = settings.get("bcftools") or {}
    pileup_cfg = bcftools_cfg.get("mpileup") or {}
    call_cfg = bcftools_cfg.get("call") or {}
    index_cfg = settings.get("reference_index") or {}

    reference = Path(reference)
    expected_bases = index_cfg.get("verify_against_expected_length")
    if expected_bases is None:
        raise PipelineError(
            "variants.reference_index.verify_against_expected_length is not set. "
            "Without it the reference index cannot be checked, and an index that "
            "under-reports the reference produces confident wrong coordinates.",
        )
    verify_reference_index(
        Path(f"{reference}.fai"), expected_bases=int(expected_bases)
    )

    converter = str(samtools_cfg.get("convert_and_index_with") or "")
    if converter != "pysam":
        # A new value in the config would otherwise be ignored and the SAM
        # handed to whatever the code happens to do.
        raise PipelineError(
            f"variants.samtools.convert_and_index_with is {converter!r}, which "
            "this adapter does not implement. Only 'pysam' is supported: the "
            "pinned samtools refuses SAM text and has no -O flag, so a new "
            "converter needs code, not a configuration value.",
            configured=converter, supported="pysam",
        )

    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    sam_path = workdir / f"{sample_id}.sam"
    bam_path = workdir / f"{sample_id}.sorted.bam"
    raw_vcf = workdir / f"{sample_id}.{RAW_VCF}"
    called_vcf = workdir / f"{sample_id}.{CALLED_VCF}"

    commands: List[Sequence[str]] = []
    align = align_command(
        tools["minimap2"],
        preset=str(minimap2_cfg.get("preset") or "asm20"),
        assembly=Path(assembly), reference=reference, output_sam=sam_path,
    )
    commands.append(align)
    with sam_path.open("w", encoding="utf-8") as handle:
        LOGGER.info("exec: %s > %s", " ".join(align), sam_path)
        try:
            completed = subprocess.run(
                align, stdout=handle, stderr=subprocess.PIPE, text=True, check=False
            )
        except OSError as exc:
            raise ToolExecutionError(
                "Failed to launch minimap2", command=" ".join(align), error=str(exc)
            ) from exc
        if completed.returncode != 0:
            raise ToolExecutionError(
                "minimap2 exited non-zero", command=" ".join(align),
                returncode=completed.returncode,
                stderr=(completed.stderr or "")[-2000:],
            )

    # The reference is passed because minimap2's SAM header names the *query*
    # contigs while the records name the reference; without the rewrite,
    # bcftools finds no sequence and silently calls nothing. See
    # `convert_and_index`.
    convert_and_index(sam_path, bam_path, reference=reference)

    pileup = mpileup_command(
        tools["bcftools"],
        bam=bam_path, reference=reference,
        min_mapq=int(pileup_cfg.get("min_mapq", 20)),
        min_bq=int(pileup_cfg.get("min_bq", 20)),
        annotate=str(pileup_cfg.get("annotate") or "FORMAT/DP,FORMAT/AD"),
        max_depth=int(pileup_cfg.get("max_depth", DEFAULT_MPILEUP_MAX_DEPTH)),
        min_ireads=int(
            pileup_cfg.get("min_ireads", DEFAULT_MPILEUP_MIN_IREADS)
        ),
        output_type=str(
            pileup_cfg.get("output_type", DEFAULT_MPILEUP_OUTPUT_TYPE)
        ),
        output_vcf=raw_vcf,
    )
    commands.append(pileup)
    _run(pileup)

    call = call_command(
        tools["bcftools"],
        input_vcf=raw_vcf, output_vcf=called_vcf,
        variants_only=bool(call_cfg.get("variants_only", True)),
        multiallelic=call_cfg.get("multiallelic", "split"),
        threads=threads,
    )
    commands.append(call)
    _run(call)

    if not called_vcf.is_file():
        raise ToolExecutionError(
            "bcftools call reported success but wrote no VCF",
            command=" ".join(call), expected=str(called_vcf),
        )

    return CallResult(sample_id=sample_id, vcf=called_vcf, commands=tuple(commands))


__all__ = [
    "CALLER_TOOLS",
    "CALLED_VCF",
    "CallResult",
    "RAW_VCF",
    "align_command",
    "call_command",
    "call_isolate",
    "convert_and_index",
    "mpileup_command",
    "require_callers",
    "resolve_executable",
    "verify_reference_index",
]
