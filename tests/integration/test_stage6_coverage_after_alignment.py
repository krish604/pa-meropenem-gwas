"""Stage 6 must measure locus coverage AFTER the caller wrote the alignments.

**The defect, as the REAL 10-isolate run reported it.** ``Stage 6 cannot
measure per-locus alignment coverage: 10 manifest sample(s) have no alignment on
disk`` (``PA_AMR_pipeline_real_run_10_isolates/PER_STAGE_OUTCOMES.md``). The
message blamed missing BAMs; the BAMs were missing because nothing had run yet.

**Where the wrong order was.** In the ``variants`` branch of ``run_pipeline``,
``derive_locus_coverage`` sat ABOVE ``stage_variants.run``. It reads
``<workdir>/<sample>/<sample>.sorted.bam`` and the only producer of those files
is ``adapters.minimap2.call_isolate``, reached from ``stage_variants.run``. On a
clean intermediate root the measurement therefore could never find an alignment
for any isolate, on any run, ever - it is an unconditional refusal, not a data
problem. It is an internal ordering bug in one branch's body; the Snakefile
prerequisite was never wrong.

**What this file does NOT do.** It does not stub the measurement. The coverage
producer under test is the REAL ``run.derive_locus_coverage``, reading REAL BAMs
written by an INJECTED fake caller (R12: synthetic inputs, injected runners).
Stubbing the measurement would have passed before the fix too, which is what
makes this a regression test rather than a description.

**Why the fake caller writes real BAMs.** ``derive_locus_coverage`` opens each
file with ``pysam`` and fetches the regulator intervals off it, so the file has
to be a real, indexed BAM over the pinned reference's contig, and it has to
actually span the loci - a synthetic row of zero coverage would measure 0.0 and
the honest threshold would then suppress every regulator record, which is
correct behaviour and would stop the run one line later for an unrelated reason.
The records are therefore built from the pinned reference's own GFF through the
production interval loader, and the query sequence is written to match the
reference span, so the measured fraction is a real 1.0 rather than a fabricated
one.
"""

from __future__ import annotations

import copy
import dataclasses
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pytest
import yaml

#: Spelled from the orchestrator's own constant rather than re-spelled here: the
#: measurement and the producer must agree on the filename, and a second
#: spelling of it in a test is how they stop agreeing silently.
from papipeline.run import SORTED_BAM_SUFFIX

REPO = Path(__file__).resolve().parents[2]
SCIENCE = REPO / "config" / "science.yaml"
SMOKE = REPO / "config" / "machines" / "smoke.yaml"

#: Two isolates: the cohort size is irrelevant to the ordering claim, and two is
#: the smallest that still exercises a per-sample loop in both halves.
ISOLATES: Tuple[str, ...] = ("GCA_000000001.1", "GCA_000000002.1")

#: A stop-free synthetic stretch. Shape only; nothing here is a biological claim.
_ASSEMBLY_SEQUENCE = "ATGGCTAGC" + "GCTGATCGATCGATCGTAGCTAGCTAGC" * 70

_DATABASE_STRING = "v6.0, light"

_FEATURE_TSV = (
    "# Annotated with Bakta\n# Software: v1.12.1\n"
    f"# Database: {_DATABASE_STRING}\n"
    "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n"
    "c1\tcds\t100\t400\t+\tTESTPA_00001\toprD\tOprD family porin\t-\n"
)
_GFF3 = (
    "##gff-version 3\n"
    "c1\tBakta\tCDS\t100\t400\t.\t+\t0\tID=TESTPA_00001;Name=TESTPA_00001\n"
)
_INFERENCE_TSV = (
    "# Software: v1.12.1\n"
    f"# Database: {_DATABASE_STRING}\n"
    "Sequence Id\tLocus Tag\tScore\tEvalue\nc1\tTESTPA_00001\t45.2\t1e-40\n"
)


@pytest.fixture()
def synthetic_real_run(tmp_path):
    """A REAL-mode run whose every input is synthetic and whose roots are tmp.

    Built as a patched copy of the COMMITTED smoke overlay rather than a
    hand-written config, so the mode under test is the one the repository ships.
    ``reuse_tool_output: require`` is asserted, not assumed: it is what keeps
    Bakta unreachable, since curated output is staged for both isolates.

    The real clinical PDC table is never read or written; the manifest comes from
    a synthetic one written here.
    """
    from papipeline.config.loader import load_config
    from papipeline.models import RunMode

    genomes = tmp_path / "smoke_genomes"
    genomes.mkdir(parents=True)
    for isolate in ISOLATES:
        (genomes / f"{isolate}.fna").write_text(
            f">c1\n{_ASSEMBLY_SEQUENCE}\n", encoding="utf-8"
        )

    pdc = tmp_path / "PDC_synthetic.tsv"
    pdc.write_text(
        "Isolate\tAssembly\tBioSample\n"
        + "".join(
            f"{isolate}\tGCA_{index:09d}.1\tSAMN{index:08d}\n"
            for index, isolate in enumerate(ISOLATES)
        ),
        encoding="utf-8",
    )
    subset = tmp_path / "isolates.txt"
    subset.write_text("\n".join(ISOLATES) + "\n", encoding="utf-8")

    database = tmp_path / "bakta-db" / "db-light"
    database.mkdir(parents=True)
    (database / "bakta.db").write_text("stub\n", encoding="utf-8")
    (database / "variant.json").write_text("stub\n", encoding="utf-8")
    (database / "version.json").write_text(
        '{"major": 6, "minor": 0, "type": "light", "date": "2025-02-24"}',
        encoding="utf-8",
    )

    results = tmp_path / "results"
    phenotype = tmp_path / "phenotypes"
    phenotype.mkdir()
    (phenotype / "phenotypes.tsv").write_text(
        "sample_id\tmeropenem\tMIC\n"
        + "".join(f"{isolate}\tR\t8\n" for isolate in ISOLATES),
        encoding="utf-8",
    )

    overlay = yaml.safe_load(SMOKE.read_text(encoding="utf-8"))
    assert overlay["annotation"]["reuse_tool_output"] == "require", (
        "this harness relies on `require` to keep Bakta unreachable; the committed "
        f"smoke overlay says {overlay['annotation']['reuse_tool_output']!r}"
    )
    overlay["paths"].update({
        "smoke_genome_dir": str(genomes),
        "pdc_table": str(pdc),
        "smoke_phenotype_dir": str(phenotype),
        "data_root": str(tmp_path / "data"),
        "results_root": str(results),
        "reports_root": str(results / "reports"),
        "intermediate_root": str(results / "intermediate"),
        "status_dir": str(tmp_path / "status"),
    })
    overlay["cohort"]["subset_file"] = str(subset)
    overlay["runtime"]["allow_real_mode"] = True
    overlay["runtime"]["observatory"] = {"enabled": False}
    overlay_path = tmp_path / "smoke-coverage-order.yaml"
    overlay_path.write_text(yaml.safe_dump(overlay), encoding="utf-8")

    config = load_config(SCIENCE, machine=overlay_path)
    raw = copy.deepcopy(dict(config.raw or {}))
    raw.setdefault("annotation", {})
    raw["annotation"]["bakta_db"] = str(database)
    config = dataclasses.replace(config, raw=raw)

    intermediate = Path(config.intermediate_root(RunMode.REAL))
    for isolate in ISOLATES:
        curated = intermediate / "bakta" / isolate
        curated.mkdir(parents=True, exist_ok=True)
        (curated / f"{isolate}.tsv").write_text(_FEATURE_TSV, encoding="utf-8")
        (curated / f"{isolate}.gff3").write_text(_GFF3, encoding="utf-8")
        (curated / f"{isolate}.inference.tsv").write_text(
            _INFERENCE_TSV, encoding="utf-8"
        )
        # Not required by `decide_reuse`; required by `run._oprd_locus_inputs`,
        # which names the curated `.faa` as one of its four declared inputs.
        (curated / f"{isolate}.faa").write_text(
            ">TESTPA_00001_oprD\nMKKIAVTQ\n", encoding="utf-8"
        )

    return {
        "config": config,
        "overlay_path": overlay_path,
        "intermediate": intermediate,
    }


@pytest.fixture()
def injected_amr(monkeypatch):
    """AMR's declared table with one row, because a synthetic genome has none.

    Stage 4 runs AMRFinderPlus for real; against a synthetic genome it correctly
    finds nothing, and the output contract then refuses a header-only table, so
    the run stops at stage 4. This writes the contract stage 4 itself declares -
    a stand-in for a tool result, not for anything under test here.
    """
    import papipeline.stages.amr as stage_amr
    from papipeline.execution.contracts import STAGE_TABLES
    from papipeline.io.tsv import write_tsv

    columns = list(STAGE_TABLES["amr"][1])
    row = {column: "-" for column in columns}
    row["sample_id"] = ISOLATES[0]
    row["gene"] = "blaOXA-1"
    row["source"] = "amrfinderplus"

    def run(config, manifest, mode, tool_output_root, *args, **kwargs):
        write_tsv(
            Path(config.intermediate_root(mode)) / "stages" / "04_amr.tsv",
            [row],
            columns,
        )
        return {}

    monkeypatch.setattr(stage_amr, "run", run)
    return row


@pytest.fixture()
def fake_blast_runner():
    """An injected command runner that writes what the blast tools would write.

    ``run_pipeline``'s ``oprd_structural_runner`` parameter exists so this can
    exist; ``None`` runs the real ``makeblastdb`` and ``tblastn``. The fake also
    writes ``<prefix>.nin``, which ``run_tblastn`` checks for, and returns an
    empty table - a VALID tblastn result meaning zero hits, which is what makes
    the absent/not-assessed verdicts reachable without fabricating an alignment.
    """
    from papipeline.adapters.external import CommandResult

    calls: List[List[str]] = []

    def runner(command):
        calls.append(list(command))
        if Path(command[0]).name == "makeblastdb":
            prefix = command[command.index("-out") + 1]
            Path(f"{prefix}.nin").write_text("synthetic index\n", encoding="utf-8")
        return CommandResult(
            command=list(command), returncode=0, stdout="", stderr="",
            duration_s=0.0,
        )

    runner.calls = calls
    return runner


@pytest.fixture()
def aligned_caller(monkeypatch, synthetic_real_run):
    """An injected ``stage_variants.run`` that writes REAL, INDEXED BAMs.

    This is the stand-in for ``adapters.minimap2.call_isolate``: the only thing
    it replaces is the aligner, which cannot produce a meaningful alignment
    between a synthetic genome and PAO1. It writes
    ``<workdir>/<sample>/<sample>.sorted.bam`` - the exact path
    ``run.derive_locus_coverage`` opens, built from ``run.SORTED_BAM_SUFFIX``
    rather than re-spelled - covering every screenable regulator locus in full,
    and indexes it, because the measurement fetches intervals off the index.

    It also returns one variant row placed INSIDE a screened locus, from the
    pinned reference's own GFF through the production interval loader and the
    pinned FASTA through pysam. Without a call inside a locus the run would stop
    at ``derive_regulator_table``'s own honest "ran and found nothing" refusal,
    which is correct and would mask the ordering claim.
    """
    import pysam

    import papipeline.stages.variants as stage_variants
    from papipeline.adapters import gff as gff_adapter

    config = synthetic_real_run["config"]
    specs = {
        gene: spec for gene, spec in dict(config.regulators).items()
        if spec.screenable
    }
    intervals = gff_adapter.load_gene_intervals(
        config.reference_gff(), [spec.locus_tag for spec in specs.values()]
    )
    loci = [
        (gene, intervals[spec.locus_tag])
        for gene, spec in sorted(specs.items())
        if spec.locus_tag in intervals
    ]
    assert loci, (
        "no regulator locus was found in the pinned reference GFF, so this "
        "harness cannot write a BAM that covers one"
    )

    reference = config.reference_fasta()
    with pysam.FastaFile(str(reference)) as fasta:
        contigs = dict(zip(fasta.references, fasta.lengths))

    written: List[Path] = []

    def write_bam(path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        header = pysam.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6", "SO": "unsorted"},
            "SQ": [{"SN": name, "LN": length} for name, length in contigs.items()],
        })
        # The reference's own bases over each span, read through pysam, so each
        # record is a faithful full-length alignment rather than a shape that
        # merely looks like one. Read inside the block: the handle the fixture
        # opened earlier is closed by the time the caller runs.
        with pysam.FastaFile(str(reference)) as source:
            # Keyed by the interval, not the contig: several loci share one
            # contig, and a per-contig key would hand every record the LAST
            # interval's sequence.
            sequences = {
                (interval.contig, interval.start, interval.end): source.fetch(
                    interval.contig, interval.start - 1, interval.end
                ).upper()
                for _gene, interval in loci
            }
        # Coordinate-sorted, in that order and not the order `loci` happens to be
        # in. `samtools index` refuses an unsorted BAM, and the refusal names the
        # offending POSITION rather than the cause, so a gene-ordered list fails
        # here with a message about a missing file.
        ordered = sorted(loci, key=lambda item: (item[1].contig, item[1].start))
        with pysam.AlignmentFile(str(path), "wb", header=header) as out:
            for _gene, interval in ordered:
                length = interval.end - interval.start + 1
                segment = pysam.AlignedSegment(header)
                segment.query_name = "synthetic"
                segment.reference_id = header.get_tid(interval.contig)
                segment.reference_start = interval.start - 1
                segment.mapping_quality = 60
                segment.cigarstring = f"{length}M"
                segment.query_sequence = sequences[
                    (interval.contig, interval.start, interval.end)
                ]
                segment.query_qualities = pysam.qualitystring_to_array(
                    "I" * length
                )
                segment.flag = 0
                out.write(segment)
        pysam.index(str(path))
        written.append(path)

    columns = list(stage_variants.PER_ISOLATE_COLUMNS)
    gene, interval = loci[0]
    position = interval.start + 1  # 1-based, strictly inside the CDS
    with pysam.FastaFile(str(reference)) as fasta:
        reference_base = fasta.fetch(interval.contig, position - 1, position).upper()
    row = {column: "-" for column in columns}
    row.update({
        "sample_id": ISOLATES[0],
        "chrom": interval.contig,
        "pos": str(position),
        "ref": reference_base,
        "alt": "ACGT"[("ACGT".index(reference_base) + 1) % 4],
    })
    assert set(row) == set(columns), (row, columns)

    def run(config, manifest, mode, *, data_root, workdir, statuses):
        assert isinstance(workdir, Path), (
            f"the caller was handed {type(workdir).__name__}, not the Path "
            "derive_locus_coverage's `variants_workdir` names"
        )
        # Out-param, filled by the real stage so stage 6a can tell "carried
        # nothing" from "never read". Both isolates here are genuinely called.
        assert isinstance(statuses, dict)
        for sample_id in manifest.sample_ids:
            write_bam(
                workdir / sample_id / f"{sample_id}{SORTED_BAM_SUFFIX}"
            )
            statuses[sample_id] = stage_variants.STATUS_CALLED
        # run.py rewrites `variants.tsv` from this RETURN value, so returning
        # rows is what makes the table non-empty; writing the file here would be
        # overwritten with a header one line later.
        return {ISOLATES[0]: [row], ISOLATES[1]: []}

    monkeypatch.setattr(stage_variants, "run", run)
    return {"bams": written, "loci": loci, "row": row}


@pytest.fixture()
def captured_coverage(monkeypatch):
    """Record the ``locus_coverage`` the regulator producer actually received.

    A wrapper around ``run.derive_regulator_table`` - not a replacement: the real
    producer still runs and still writes its table. What is captured is the
    argument the orchestrator computed, which is the thing under test.
    """
    import papipeline.run as run

    seen: Dict[str, Any] = {}
    original = run.derive_regulator_table

    def recording(config, manifest, **kwargs):
        seen["locus_coverage"] = kwargs.get("locus_coverage")
        seen["min_locus_coverage"] = kwargs.get("min_locus_coverage")
        return original(config, manifest, **kwargs)

    monkeypatch.setattr(run, "derive_regulator_table", recording)
    return seen


class TestStageSixMeasuresCoverageAfterTheAlignmentsAreWritten:
    """The ordering claim, driven through `run_pipeline` on synthetic input."""

    @staticmethod
    def _drive(env, runner):
        """Run the orchestrator; return ``(result_or_None, stage_error_or_None)``.

        The pair shape exists so a test can say WHICH failure occurred. Before
        the fix the drive raises the stage 6 coverage refusal, and that specific
        refusal is the regression signal - a `TypeError` or a contract refusal
        would be a different defect and must not be reported as this one.
        """
        from papipeline.errors import StageError
        from papipeline.run import run_pipeline

        try:
            result = run_pipeline(
                config_path=SCIENCE,
                config=env["config"],
                mode="REAL",
                only=["variants"],
                write_html=False,
                machine=env["overlay_path"],
                oprd_structural_runner=runner,
            )
        except StageError as exc:
            return None, exc
        return result, None

    def test_stage_six_completes_on_a_clean_intermediate_root(
        self, synthetic_real_run, fake_blast_runner, injected_amr, aligned_caller
    ):
        """THE REGRESSION. Before the fix this raised; now stage 6 completes.

        The intermediate root starts with NO BAM anywhere - the only writer is
        the injected caller, and it writes them when stage 6 reaches it. So this
        passes only if the measurement happens after that call, which is exactly
        the claim. Any run that measured earlier raised with `n_missing` equal to
        the cohort size.
        """
        variants_root = synthetic_real_run["intermediate"] / "variants"
        assert not variants_root.exists(), (
            f"{variants_root} already exists, so the harness is not starting "
            "from the clean intermediate root the ordering claim is about"
        )

        result, error = self._drive(synthetic_real_run, fake_blast_runner)
        assert error is None, (
            "stage 6 failed on a clean intermediate root:\n"
            f"{error}\n\nIf this is the `cannot measure per-locus alignment "
            "coverage ... have no alignment on disk` refusal, the measurement is "
            "running before the caller wrote the BAMs. Every isolate will be "
            "reported missing, because on a fresh tree NO isolate has an "
            "alignment yet."
        )
        assert result is not None
        assert result.stage_status.get("variants") == "completed", (
            f"stage 6 did not complete: {result.stage_status}"
        )

    def test_the_bams_the_measurement_read_were_written_by_the_caller(
        self, synthetic_real_run, fake_blast_runner, injected_amr, aligned_caller
    ):
        """The producer really ran, at the path the consumer opens.

        Guards against the fix being made by MOVING the measurement somewhere
        that no longer depends on the aligner: if the caller never ran, or wrote
        elsewhere, this fails even though stage 6 might still complete.
        """
        result, error = self._drive(synthetic_real_run, fake_blast_runner)
        assert error is None, f"stage 6 failed: {error}"

        expected = sorted(
            synthetic_real_run["intermediate"] / "variants" / sample_id
            / f"{sample_id}{SORTED_BAM_SUFFIX}"
            for sample_id in ISOLATES
        )
        assert sorted(aligned_caller["bams"]) == expected, (
            f"the caller wrote {sorted(aligned_caller['bams'])}, not {expected}"
        )
        for bam in expected:
            assert bam.is_file(), f"{bam} is not on disk"
            assert Path(f"{bam}.bai").is_file(), (
                f"{bam} was not indexed; the measurement fetches intervals off "
                "the index, so an unindexed BAM would measure nothing rather "
                "than fail"
            )

    def test_the_coverage_handed_to_the_screen_was_measured_not_stubbed(
        self, synthetic_real_run, fake_blast_runner, injected_amr,
        aligned_caller, captured_coverage
    ):
        """The screen received REAL measured coverage for every (isolate, locus).

        Non-vacuity in the direction that matters: the previous version of this
        harness had to INJECT the measurement, because a synthetic genome aligned
        against PAO1 measures 0.0 everywhere. This one does not, because the
        injected caller writes alignments that genuinely span the loci - so a
        coverage of 1.0 here is a measurement, and the measurement is exactly
        what the ordering defect prevented.

        The threshold is the module's own constant, read through the live call
        site rather than typed in, and it is asserted to be present - a screen
        handed coverage with no threshold is the gap this file's subject came
        from.
        """
        result, error = self._drive(synthetic_real_run, fake_blast_runner)
        assert error is None, f"stage 6 failed: {error}"

        coverage = captured_coverage.get("locus_coverage")
        assert coverage is not None, (
            "the regulator producer was handed no coverage at all; stage 6 then "
            "reports a clean screen over loci no aligner reached"
        )
        genes = [gene for gene, _interval in aligned_caller["loci"]]
        for sample_id in ISOLATES:
            for gene in genes:
                value = coverage.get((sample_id, gene))
                assert value is not None, (
                    f"no coverage was measured for ({sample_id}, {gene})"
                )
                assert value == 1.0, (
                    f"({sample_id}, {gene}) measured {value}, not 1.0, even "
                    "though the injected caller wrote an alignment spanning the "
                    "locus in full. A partial measurement here means the reader "
                    "and the writer disagree about the file or the interval."
                )
        assert captured_coverage.get("min_locus_coverage") is not None, (
            "the producer was not given a coverage threshold, so nothing could "
            "ever be suppressed"
        )

    def test_bakta_was_never_invoked_by_this_drive(
        self, synthetic_real_run, fake_blast_runner, injected_amr,
        aligned_caller, monkeypatch
    ):
        """The harness holds `require`, so no tool is executed for annotation.

        Reported because the other half of this round's work is about executing
        Bakta: a regression test that quietly ran the tool would be worse than
        one that failed.
        """
        import papipeline.adapters.bakta as bakta_adapter

        def forbidden(*args, **kwargs):
            raise AssertionError(
                "Bakta was invoked under reuse_tool_output=require; curated "
                "output is staged, so the tool must not be reached"
            )

        monkeypatch.setattr(bakta_adapter, "run_bakta", forbidden)
        self._drive(synthetic_real_run, fake_blast_runner)
        # No assertion: `forbidden` raising IS the failure.