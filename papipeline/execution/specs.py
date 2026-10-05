"""Declared output contracts for the stages that have a real one.

These are *data*, not code paths. A spec says what a stage must produce;
:mod:`papipeline.execution.validation` turns that into a verdict. Adding a
contract for a stage never changes that stage's behaviour, and a stage with
no contract still runs - it simply reports "no output spec declared" rather
than claiming a validated success.

Each spec encodes a contract that already exists in the repository. No
biological threshold is introduced here. Where a check mentions a specific
locus or column it is because an existing stage already requires it
(``config/gene_families.tsv``, ``papipeline/stages/annotation.py``,
``papipeline/stages/validation.py``).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from .validation import (
    Check,
    CheckKind,
    OutputSpec,
    SiblingSpec,
    covers_samples,
    exists,
    has_columns,
    has_real_values,
    identifiers_known,
    identifies_file,
    min_rows,
    non_empty,
    not_constant_matrix,
    parses,
)

#: The columns the pipeline's own annotation reader requires
#: (``annotation.read_standardised`` -> ``required_columns=("sample_id",
#: "gene_name")``). A file without them cannot be consumed downstream, so a
#: stage that produced one has not succeeded.
ANNOTATION_REQUIRED_COLUMNS = ("sample_id", "gene_name")

#: Columns that identify an annotation table as the *feature* table rather
#: than the expert-system inference table. This is the check for a bug that
#: shipped: ``bakta_outputs`` returned ``*.inference.tsv``, every gene name
#: was lost, and all 11 stage-6 loci scored 0/65 - a uniform false negative
#: that read exactly like a biological result.
BAKTA_ANNOTATION_MARKER_COLUMNS = ("Gene", "Locus Tag")
BAKTA_INFERENCE_MARKER_COLUMNS = ("Score", "Id", "Evalue")

#: Column names the pipeline's own Bakta reader produces. Checked alongside
#: the raw names because the normalised view is what stages 6 and 9 consume,
#: while the raw view is what distinguishes a feature table from an
#: inference table - both normalise to the same key set.
BAKTA_NORMALISED_COLUMNS = ("seqid", "gene_id", "gene_name")

#: Loci the pipeline already tracks in ``config/gene_families.tsv`` and
#: ``config/mechanisms.tsv``. Used only to prove the annotation table is
#: resolvable, never to assert a count.
STAGE6_LOCI = ("oprD", "nalC", "nalD", "nfxB", "mexT", "mexS", "mexR",
               "mexZ", "ampR", "ampD", "dacB")


@dataclass(frozen=True)
class BaktaOutputs:
    """Resolved paths for one genome's Bakta output directory."""

    out_dir: Path
    sample_id: str
    genome_stem: str

    @property
    def gff(self) -> Path:
        return self.out_dir / f"{self.genome_stem}.gff3"

    @property
    def annotation_tsv(self) -> Path:
        return self.out_dir / f"{self.genome_stem}.tsv"

    @property
    def inference_tsv(self) -> Path:
        return self.out_dir / f"{self.genome_stem}.inference.tsv"


def bakta_spec(
    outputs: BaktaOutputs,
    *,
    require_inference: bool = False,
    require_faa: bool = False,
) -> OutputSpec:
    """Contract for one genome's annotation.

    Checks, in the order they matter:

    1. the GFF and the annotation TSV exist and are non-empty;
    2. the TSV parses with the pipeline's strict reader;
    3. it carries the columns the pipeline requires;
    4. it carries the feature-table marker columns, which is what
       distinguishes it from the inference table;
    5. it contains at least one real gene name, so "shape but no data" is
       rejected;
    6. every row belongs to this sample.

    ``require_faa`` exists because the current contract does not consume
    protein FASTA; it defaults off and is available for a stage that does.
    No check here judges whether the annotation is *biologically* right.
    """
    checks = [
        non_empty(outputs.gff, description="annotation GFF is present and non-empty"),
        parses(outputs.gff, "gff", description="GFF is line-oriented text with a directive"),
        identifies_file(outputs.gff, outputs.genome_stem,
                        description="GFF is named after the genome that was submitted"),
        non_empty(outputs.annotation_tsv,
                  description="annotation table is present and non-empty"),
        parses(outputs.annotation_tsv, "bakta",
               description="annotation table parses with the Bakta reader"),
        has_columns(outputs.annotation_tsv, BAKTA_ANNOTATION_MARKER_COLUMNS,
                    dialect="bakta_raw",
                    description="is the feature table, not the inference table"),
        has_columns(outputs.annotation_tsv, BAKTA_NORMALISED_COLUMNS,
                    dialect="bakta",
                    description="normalises to the columns stages 6 and 9 consume"),
        has_real_values(outputs.annotation_tsv, "gene_name", dialect="bakta",
                        description="at least one row carries a resolvable gene name"),
        identifies_file(outputs.annotation_tsv, outputs.genome_stem,
                        description="table is named after the genome that was submitted"),
    ]
    if require_inference:
        checks.append(non_empty(outputs.inference_tsv,
                               description="inference table is present"))
    if require_faa:
        checks.append(non_empty(outputs.out_dir / f"{outputs.genome_stem}.faa",
                               description="protein FASTA is present"))
    return OutputSpec(
        stage="annotation",
        checks=tuple(checks),
        expectations={
            "genome_stem": outputs.genome_stem,
            "sample_id": outputs.sample_id,
        },
    )


def annotation_sibling_spec(annotation_root: Path, stem_of) -> SiblingSpec:
    """Cross-genome degeneracy check for a cohort's annotation.

    Catches the class of failure where every genome's annotation is the same
    artefact, which no per-genome check can see.
    """
    return SiblingSpec(
        root=annotation_root,
        pattern="*/*.gff3",
        min_distinct=2,
        label="annotation GFF",
        description="per-genome annotations are not all byte-identical",
    )


def variant_call_spec(
    vcf: Path,
    sample_id: str,
    *,
    min_variants: int = 1,
) -> OutputSpec:
    """Contract for one isolate's variant calls.

    A VCF with no variant records, or one that does not mention its own
    sample, has not done the job. ``min_variants`` defaults to 1: a call
    with zero records is structurally a failure of the calling step for a
    real assembly, not a biological result. It is a parameter rather than a
    constant so a caller can state its own floor explicitly.
    """
    return OutputSpec(
        stage="variant_calling",
        checks=(
            non_empty(vcf, description="VCF is present and non-empty"),
            parses(vcf, "vcf", description="VCF header and records are well formed"),
            Check(CheckKind.MIN_ROWS, path=vcf, minimum=min_variants, dialect="vcf",
                  description="VCF carries variant records"),
        ),
        expectations={"sample_id": sample_id, "min_variants": min_variants},
    )


def variant_sibling_spec(vcf_root: Path, sample_id: str) -> SiblingSpec:
    """Cross-isolate degeneracy check for per-isolate VCFs.

    This is the check for the bug where a leaked loop variable made every
    isolate's VCF byte-identical, so the core-SNP alignment had exactly one
    distinct sequence and the tree had no resolution.
    """
    return SiblingSpec(
        root=vcf_root,
        pattern=f"{sample_id}*",
        min_distinct=2,
        label="variant call",
        description="per-isolate VCFs are not all byte-identical",
    )


def kinship_spec(kinship: Path, expected_labels: Sequence[str]) -> OutputSpec:
    """Contract for a pyseer similarity/kinship matrix.

    The decisive check is :func:`not_constant_matrix`. A kinship matrix
    where every entry is 1.0 says "all samples are identical", which cannot
    be true of distinct isolates; that exact bug produced h^2 = 1.00 and a
    set of associations with no plausible mechanism. It is rejected here as
    ``INVALID`` rather than being carried into the association test.
    """
    labels = [str(label) for label in expected_labels]
    return OutputSpec(
        stage="gwas",
        checks=(
            non_empty(kinship, description="kinship matrix is present and non-empty"),
            parses(kinship, "tsv", description="matrix is readable tabular data"),
            not_constant_matrix(kinship,
                                description="matrix is not all-1.0/constant"),
            covers_samples(kinship, len(labels), key="sample",
                           description="one row per isolate"),
        ),
        expectations={"n_labels": len(labels), "labels": labels[:5]},
    )


def standardised_annotation_spec(
    path: Path, sample_id: str, *, min_records: int = 1
) -> OutputSpec:
    """Contract for the standardised annotation TSV stage 2 and stage 9 read."""
    return OutputSpec(
        stage="standardised_annotation",
        checks=(
            non_empty(path, description="standardised annotation is non-empty"),
            parses(path, "tsv", description="passes the strict TSV reader"),
            has_columns(path, ANNOTATION_REQUIRED_COLUMNS,
                        description="carries the columns stage 9 requires"),
            has_real_values(path, "gene_name",
                            description="at least one record carries a gene name"),
            min_rows(path, min_records, description="has at least one record"),
            identifiers_known(path, (sample_id,), key="sample_id",
                              description="every row belongs to this isolate"),
        ),
        expectations={"sample_id": sample_id},
    )


def stage6_resolution_spec(path: Path, sample_id: str) -> OutputSpec:
    """Contract for the table stage 6 reads to resolve AMR loci.

    This separates two situations that look identical in a results table
    but are not: a sample that genuinely carries none of the tracked loci,
    and a sample whose input never contained the loci because the wrong
    file was selected. The second is a pipeline failure and is reported as
    such here.
    """
    return OutputSpec(
        stage="amr_locus_resolution",
        checks=(
            non_empty(path, description="AMR determinant table is non-empty"),
            parses(path, "tsv", description="passes the strict TSV reader"),
            has_columns(path, ("sample_id", "Gene"),
                        description="carries sample and gene columns"),
            has_real_values(path, "Gene", description="at least one gene name is resolvable"),
            identifiers_known(path, (sample_id,), key="sample_id",
                              description="every row belongs to this isolate"),
        ),
        expectations={
            "sample_id": sample_id,
            "tracked_loci": list(STAGE6_LOCI),
            "note": "tracked loci are used to prove the input was resolvable, "
                    "never to assert how many should be present",
        },
    )
