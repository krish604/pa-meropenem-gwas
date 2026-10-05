"""Stage 5 - P. aeruginosa resistance mechanism interpretation.

This is the layer that turns a raw gene detection or a variant call into a
named mechanism, using ``config/mechanisms.tsv``.

The central guarantee: interpretation can only ever *lower* a claim to the
gene's configured ceiling, never raise it. A determinant that reaches this
stage as ``DETECTED`` leaves as ``DETECTED``. Only stages 12 and 13, which
have phenotype and phylogeny in hand, may move a claim to ``ASSOCIATED`` or
``SUPPORTED``.

The stage is also strict about OprD. Detecting the ``oprD`` gene is evidence
of an intact locus, which is reported as ``reduced_permeability`` /
``locus_present`` and is *not* evidence of reduced permeability. Loss of the
locus is what carries the mechanism, and that arrives from stage 6 or 7.

A variant **allele** is a different observation from the presence of the gene,
and this stage now says so. ``oprD_V359L`` used to be interpreted by gene name
alone, so it produced the same output as intact ``oprD``: mechanism
``reduced_permeability``, evidence ``DETECTED``, and a note claiming the locus
is intact - in one row, two contradictory statements, one of which
``config/mechanisms.tsv`` explicitly forbids (*"Presence of oprD indicates an
intact locus... Loss-of-function must come from stage 6 or stage 7"*). Since
stage 11 adopted ``--organism`` this is live on real data, not hypothetical.
:func:`interpret_determinant` therefore drops the gene's presence-based note for
a variant allele and substitutes one that says what the row actually observed,
and :func:`annotate_intact_locus` declines to certify a locus that carries a
variant. No ``ClaimStatus`` is substituted - see
:data:`_VARIANT_ALLELE_NOTE` for why ``PREDICTED`` cannot mean what it looks
like it means.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from ..config.loader import PipelineConfig
from ..knowledge import ceiling_for, mechanism_map_for_antibiotic
from ..logging_utils import get_logger
from ..manifest import SampleManifest
from ..models import (
    AmrDeterminant,
    ClaimStatus,
    MechanismCall,
    RegulatorVariant,
    RunMode,
    StructuralCallStatus,
    StructuralVariant,
)

LOGGER = get_logger("stages.mechanisms")

MECHANISM_COLUMNS: Tuple[str, ...] = (
    "sample_id",
    "antibiotic",
    "determinant",
    "mechanism",
    "evidence_level",
    "gene",
    "gene_type",
    "notes",
)

#: Genes whose *presence* does not constitute the mechanism. For these, a
#: detection is reported as evidence of an intact locus only.
LOSS_OF_FUNCTION_REQUIRED = frozenset({"oprD"})

#: Note prefix for a determinant that is a variant **allele** rather than a
#: detection of the gene.
#:
#: No :class:`ClaimStatus` is substituted here, and that is a deliberate
#: non-decision rather than an oversight. `ClaimStatus` is ordered
#: SUPPORTED > ASSOCIATED > PREDICTED > DETECTED > UNKNOWN
#: (`models.CLAIM_STATUS_PRECEDENCE`), so `PREDICTED` is a *stronger* claim than
#: `DETECTED` and cannot express a downgrade - `ceiling_for(config, "oprD",
#: PREDICTED)` hands back DETECTED, because DETECTED is oprD's ceiling. The
#: claim the tool actually made - an allele was detected - is true, so
#: `DETECTED` stands. What was wrong was never the strength of the claim; it was
#: that the row asserted the *presence* reading alongside it, in the knowledge
#: table's own words, for a locus the same row says is altered. That is a
#: contradiction, and it is removed by not making the claim rather than by
#: inventing a status to hold it down with.
_VARIANT_ALLELE_NOTE = (
    "variant allele called at this locus, not a detection of the gene; an "
    "altered allele is not a locus-integrity observation, and loss-of-function "
    "must be established by stage 6 or stage 7"
)


def _notes_for(
    gene: Optional[str], mechanism: str, detail: str, spec_notes: Optional[str]
) -> str:
    parts = [detail]
    if spec_notes:
        parts.append(spec_notes)
    return " | ".join(p for p in parts if p)


def _is_variant_allele(gene: str, determinant: AmrDeterminant) -> bool:
    """Whether this determinant is an allele call rather than a gene detection.

    Only meaningful for the genes in :data:`LOSS_OF_FUNCTION_REQUIRED`, which is
    why the caller checks that set first: for every other gene in the knowledge
    table a variant is the observation of interest (``nalC`` variants are what
    the table says it records) and its configured note is written for exactly
    that case. oprD's note is written for the opposite case, which is the whole
    contradiction.
    """
    return bool(determinant.variant)


def interpret_determinant(
    config: PipelineConfig,
    determinant: AmrDeterminant,
    antibiotic: str,
) -> Optional[MechanismCall]:
    """Interpret one acquired determinant.

    Returns ``None`` when the determinant is not in the knowledge table, so
    the caller can count it as unmapped rather than inventing a mechanism.

    A determinant carrying a ``variant`` at a
    :data:`LOSS_OF_FUNCTION_REQUIRED` gene is an allele call, not a detection of
    the gene. Reading only the gene name - as this did until the round-12 fix -
    made ``oprD_V359L`` byte-identical to intact ``oprD`` and emitted
    ``config/mechanisms.tsv``'s own presence claim (*"Presence of oprD indicates
    an intact locus"*) beside an altered locus, in one row. See
    :data:`_VARIANT_ALLELE_NOTE` for why no ``ClaimStatus`` is substituted.
    """
    gene = determinant.gene or determinant.determinant.split("~")[0]
    spec = mechanism_map_for_antibiotic(config, antibiotic).get(gene)
    if spec is None:
        LOGGER.debug(
            "Determinant %s has no mechanism mapping; leaving unmapped", gene
        )
        return None

    is_allele = gene in LOSS_OF_FUNCTION_REQUIRED and _is_variant_allele(
        gene, determinant
    )
    if is_allele:
        detail = _VARIANT_ALLELE_NOTE
        spec_notes = None
        LOGGER.debug(
            "Determinant %s carries variant %s at a loss-of-function locus; "
            "the gene's presence note does not apply to an altered allele",
            determinant.determinant,
            determinant.variant,
        )
    else:
        detail = (
            "acquired determinant detected; detection does not establish "
            "phenotypic resistance"
        )
        spec_notes = spec.notes

    evidence = ceiling_for(config, gene, determinant.claim_status)
    return MechanismCall(
        sample_id=determinant.sample_id,
        antibiotic=antibiotic,
        determinant=determinant.determinant,
        mechanism=spec.mechanism,
        evidence_level=evidence,
        gene=gene,
        gene_type=spec.gene_type,
        notes=_notes_for(gene, spec.mechanism, detail, spec_notes),
    )


def interpret_variant(
    config: PipelineConfig,
    variant: RegulatorVariant,
    antibiotic: str,
) -> Optional[MechanismCall]:
    """Interpret one chromosomal variant call.

    The ``GENE_ABSENCE`` / disruptive variant types carry the mechanism for
    loci such as oprD whose *loss* is the relevant observation.

    A variant at a :data:`LOSS_OF_FUNCTION_REQUIRED` gene is an allele call, so
    the gene's presence-based note is dropped exactly as in
    :func:`interpret_determinant`. This is not a hypothetical on real data: the
    stage-6 screen on the 10 smoke isolates emitted **769** oprD rows
    (`round12/cooc/o5_report.json`), every one of which carried
    ``config/mechanisms.tsv``'s *"Presence of oprD indicates an intact locus"*
    beside a description of a variant in oprD.
    """
    spec = mechanism_map_for_antibiotic(config, antibiotic).get(variant.gene)
    if spec is None:
        LOGGER.debug("No mechanism mapping for gene %s", variant.gene)
        return None

    if variant.gene in LOSS_OF_FUNCTION_REQUIRED:
        detail = _VARIANT_ALLELE_NOTE
        spec_notes = None
    else:
        detail = f"variant {variant.variant_type} detected in {variant.gene}"
        spec_notes = spec.notes

    evidence = ceiling_for(config, variant.gene, variant.call_status)
    return MechanismCall(
        sample_id=variant.sample_id,
        antibiotic=antibiotic,
        determinant=variant.variant,
        mechanism=spec.mechanism,
        evidence_level=evidence,
        gene=variant.gene,
        gene_type=spec.gene_type,
        notes=_notes_for(variant.gene, spec.mechanism, detail, spec_notes),
    )


def interpret_structural(
    config: PipelineConfig,
    sv: StructuralVariant,
    antibiotic: str,
) -> Optional[MechanismCall]:
    """Interpret a structural variant.

    A ``candidate`` call is downgraded to ``PREDICTED`` and a
    ``not_assessable`` call yields no mechanism call at all. Neither is
    ever reported at ``DETECTED``.
    """
    if not sv.affected_gene:
        return None
    spec = mechanism_map_for_antibiotic(config, antibiotic).get(sv.affected_gene)
    if spec is None:
        return None

    if sv.call_status is StructuralCallStatus.NOT_ASSESSABLE:
        LOGGER.debug(
            "SV %s is not assessable; no mechanism call emitted", sv.variant_id
        )
        return None

    if sv.call_status is StructuralCallStatus.CANDIDATE:
        evidence = ClaimStatus.PREDICTED
    else:
        evidence = ceiling_for(config, sv.affected_gene, ClaimStatus.DETECTED)

    return MechanismCall(
        sample_id=sv.sample_id,
        antibiotic=antibiotic,
        determinant=sv.variant_id,
        mechanism=spec.mechanism,
        evidence_level=evidence,
        gene=sv.affected_gene,
        gene_type=spec.gene_type,
        notes=_notes_for(
            sv.affected_gene,
            spec.mechanism,
            f"structural {sv.variant_type} ({sv.call_status.value}) affecting "
            f"{sv.affected_gene}",
            spec.notes,
        ),
    )


def annotate_intact_locus(
    config: PipelineConfig,
    sample_id: str,
    gene: str,
    antibiotic: str,
    source: str,
    *,
    variant_allele: Optional[str] = None,
) -> Optional[MechanismCall]:
    """Report the *presence* of a loss-of-function-relevant locus.

    This deliberately records ``locus_present`` and stays at ``DETECTED``.
    It exists so the master table can distinguish "gene seen" from "gene
    disrupted", which is the distinction rule 1 depends on.

    ``variant_allele`` is the ``V359L`` in ``oprD_V359L`` for the same sample
    and gene, when one exists. A locus carrying a variant allele is not an
    intact locus in the sense this call asserts, so the call is **declined**
    rather than emitted at a lower status: there is no status of "intact but
    altered", and emitting one of the two would state the thing the variant
    contradicts. ``run`` records the supersession in the variant row's own
    notes, so the omission is stated in the table rather than being silent.
    """
    if gene not in LOSS_OF_FUNCTION_REQUIRED:
        return None
    if variant_allele:
        LOGGER.debug(
            "%s carries variant %s at %s; declining to certify the locus as "
            "intact",
            sample_id,
            variant_allele,
            gene,
        )
        return None
    spec = mechanism_map_for_antibiotic(config, antibiotic).get(gene)
    if spec is None:
        return None
    return MechanismCall(
        sample_id=sample_id,
        antibiotic=antibiotic,
        determinant=gene,
        mechanism="locus_intact",
        evidence_level=ClaimStatus.DETECTED,
        gene=gene,
        gene_type=spec.gene_type,
        notes=(
            f"{gene} locus detected ({source}); presence indicates an intact "
            "locus and is not evidence of reduced permeability or of "
            "susceptibility"
        ),
    )


def run(
    config: PipelineConfig,
    manifest: SampleManifest,
    mode: RunMode,
    antibiotic: str,
    amr_calls: Dict[str, List[AmrDeterminant]],
    regulator_variants: Optional[Dict[str, List[RegulatorVariant]]] = None,
    structural_variants: Optional[Dict[str, List[StructuralVariant]]] = None,
    annotated_genes: Optional[Dict[str, Sequence[str]]] = None,
    *,
    intermediate_root: Optional[Path] = None,
) -> Dict[str, List[MechanismCall]]:
    """Stage 5 entry point.

    Args:
        amr_calls: Stage 4 output, sample -> determinants. **Required**, and
            deliberately so: `convergence` and `cooccurrence` make their
            equivalent optional, and copying that would let a caller reach this
            stage with no stage-4 result at all and receive an empty table that
            reads as "no mechanism found". Pass `{}` to mean "there is
            nothing"; do not make the argument omissible. In REAL these are
            **re-read from disk** and the argument is ignored, because a
            standalone or Snakemake invocation has no in-memory stage-4 result
            to inherit.
        regulator_variants: Stage 6 output. Re-read from disk in REAL.
        structural_variants: Stage 7 output. Re-read from disk in REAL.
        annotated_genes: Stage 2 output reduced to gene names per sample, used
            only to report intact OprD loci. Re-read from disk in REAL.
        intermediate_root: REAL only. Where the four inputs above are read
            from. Required in REAL; its absence is a refusal, not a fallback to
            the passed mappings.

    Raises:
        StageError: In REAL, when ``intermediate_root`` is absent or a declared
            input's table is absent. Nothing is written, so no partial table is
            left behind.
        NotImplementedError: In STUB, which fabricates this table instead.
    """
    if mode is not RunMode.TEST:
        if mode is RunMode.REAL:
            amr_calls, regulator_variants, structural_variants, annotated_genes = (
                _run_real_inputs(
                    config,
                    manifest,
                    mode,
                    antibiotic,
                    intermediate_root=intermediate_root,
                )
            )
        else:
            # STUB fabricates this table in `papipeline.stub`, so reaching the
            # real computation here would mean STUB stopped stubbing.
            raise NotImplementedError(
                "STUB-mode mechanisms is fabricated by papipeline.stub and has "
                f"no caller here. (mode={getattr(mode, 'value', mode)})"
            )

    config.require_antibiotic(antibiotic)
    mechanism_map = mechanism_map_for_antibiotic(config, antibiotic)
    calls: Dict[str, List[MechanismCall]] = {sid: [] for sid in manifest.sample_ids}

    #: (sample, gene) -> the variant allele seen there. A locus carrying one of
    #: these is not certified intact; see `annotate_intact_locus`.
    variant_alleles: Dict[Tuple[str, str], str] = {}

    unmapped = 0
    for sample_id, determinants in amr_calls.items():
        if sample_id not in calls:
            continue
        for determinant in determinants:
            gene = determinant.gene or determinant.determinant.split("~")[0]
            if determinant.variant and gene in LOSS_OF_FUNCTION_REQUIRED:
                variant_alleles.setdefault((sample_id, gene), str(determinant.variant))
            call = interpret_determinant(config, determinant, antibiotic)
            if call is None:
                unmapped += 1
                continue
            calls[sample_id].append(call)

    for sample_id, variants in (regulator_variants or {}).items():
        for variant in variants:
            call = interpret_variant(config, variant, antibiotic)
            if call is not None and sample_id in calls:
                calls[sample_id].append(call)

    for sample_id, variants in (structural_variants or {}).items():
        for sv in variants:
            call = interpret_structural(config, sv, antibiotic)
            if call is not None and sample_id in calls:
                calls[sample_id].append(call)

    for sample_id, genes in (annotated_genes or {}).items():
        if sample_id not in calls:
            continue
        for gene in genes:
            if gene not in LOSS_OF_FUNCTION_REQUIRED or gene not in mechanism_map:
                continue
            allele = variant_alleles.get((sample_id, gene))
            call = annotate_intact_locus(
                config,
                sample_id,
                gene,
                antibiotic,
                source="annotation",
                variant_allele=allele,
            )
            if call is None:
                # No locus_intact row for a locus carrying a variant allele. The
                # variant row's own `notes` already say the locus-integrity
                # reading does not apply, so the table states why the
                # certification is absent instead of dropping it silently.
                LOGGER.debug(
                    "%s: no intact-locus annotation for %s (variant %s)",
                    sample_id,
                    gene,
                    allele,
                )
                continue
            calls[sample_id].append(call)

    total = sum(len(v) for v in calls.values())
    by_mechanism: Dict[str, int] = {}
    for records in calls.values():
        for call in records:
            by_mechanism[call.mechanism] = by_mechanism.get(call.mechanism, 0) + 1

    LOGGER.info(
        "Stage 5: %d mechanism calls across %d mechanisms | unmapped "
        "determinants: %d",
        total,
        len(by_mechanism),
        unmapped,
    )
    if by_mechanism:
        LOGGER.info(
            "Stage 5 mechanism distribution: %s",
            ", ".join(f"{k}={v}" for k, v in sorted(by_mechanism.items())),
        )
    return calls


#: What holds each REAL input back, for the refusal message. Mirrors the shape
#: used by `convergence` and `cooccurrence` so a reader who has met one
#: recognises the other.
#:
#: Each names the *stage that produces the file*. Note what is NOT said: none of
#: these claims the table will be non-empty once the stage has run, because a
#: present-but-empty table is a result this stage accepts (see
#: `_run_real_inputs`).
_MECHANISM_BLOCKERS = {
    "amr_calls": (
        "stage 4 (amr), which runs AMRFinderPlus per isolate and writes "
        "`04_amr.tsv`"
    ),
    "regulator_variants": (
        "stage 6 (variants), whose `regulators/regulator_variants.tsv` is the "
        "only route to an oprD loss-of-function mechanism call. Its producer "
        "also writes `regulators/regulator_screen.json` beside it"
    ),
    "structural_variants": (
        "stage 7 (structural variants), which writes the folded "
        "`07_structural_variants.tsv`"
    ),
    "annotated_genes": (
        "stage 2 (annotation), whose per-isolate "
        "`<intermediate_root>/annotation/<sample_id>.annotation.tsv` is the "
        "only on-disk record carrying gene *names*. `02_annotation_summary.tsv` "
        "carries counts, so the intact-locus annotation cannot be recovered "
        "from it"
    ),
}


def real_input_paths(
    config: PipelineConfig,
    intermediate_root: Path,
    mode: RunMode,
) -> Dict[str, str]:
    """Where each REAL input is read from, for the refusal message.

    Stage-table locations come from ``contracts`` so this cannot name a stale
    path. ``annotated_genes`` is a *directory* of per-isolate tables, not a
    file, because that is the only shape annotation writes.

    Every path is derived from ``intermediate_root`` rather than from the
    config's own roots, including the annotation directory. In REAL those are
    the same directory, but a caller holding the root somewhere else - a test, a
    Snakemake rule with a scratch root - must have all four inputs move
    together. Mixing the two would read three inputs from the tree under test
    and the fourth from `results/real`, which is the arrangement that makes a
    REAL run pick up a previous run's annotations.
    """
    from ..execution.contracts import internal_table_path, table_path

    stage_dir = Path(intermediate_root) / "stages"
    return {
        "amr_calls": str(table_path(stage_dir, "amr")),
        "regulator_variants": str(
            Path(intermediate_root) / "regulators" / "regulator_variants.tsv"
        ),
        "structural_variants": str(internal_table_path(stage_dir, "structural_variants")),
        "annotated_genes": str(Path(intermediate_root) / "annotation"),
    }


def _present_but_empty(name: str, path: Path, exc: Exception) -> bool:
    """Whether ``exc`` means *present with no rows* rather than *malformed*.

    `io.tsv.read_tsv` refuses a header-only table, and
    `regulators.load_regulator_variants` refuses it too - deliberately, and for
    the regulator table with a message that says so. This stage needs the other
    reading: the producer ran and found nothing, which is a result (see
    `_run_real_inputs`), and the file being *there* is the evidence that it ran.

    Narrow on purpose. A genuinely malformed table - too few fields, a missing
    required column - must still propagate, because "the table is broken" and
    "the table is empty" are different failures and collapsing them is how a
    parse error becomes a cohort carrying no mechanism call.
    """
    message = str(exc)
    if "no data rows" in message or "That is a result" in message:
        LOGGER.info(
            "Stage 5 input %s at %s is present and carries no rows. The "
            "producing stage ran and found nothing; stage 5 continues.",
            name, path,
        )
        return True
    return False


def load_amr_calls(
    path: Path,
    antibiotic: str,
    mechanism_map: Mapping[str, object],
) -> Dict[str, List[AmrDeterminant]]:
    """Sample -> stage-4 determinants, read from ``04_amr.tsv``.

    A missing file yields ``{}``; ``_run_real_inputs`` refuses on that before
    this is reached. A present-but-empty file also yields ``{}``, because
    ``read_tsv`` raises on a header-only table and that is a result rather than a
    failure - see :func:`_present_but_empty`.
    """
    from ..errors import DataContractError
    from .amr import load_amr_table

    table = Path(path)
    if not table.is_file():
        return {}
    try:
        records = load_amr_table(
            table, antibiotic=antibiotic, mechanism_map=mechanism_map
        )
    except DataContractError as exc:
        if not _present_but_empty("amr_calls", table, exc):
            raise
        return {}
    by_sample: Dict[str, List[AmrDeterminant]] = {}
    for determinant in records:
        by_sample.setdefault(str(determinant.sample_id), []).append(determinant)
    return by_sample


def load_regulator_calls(
    config: PipelineConfig,
    path: Path,
    intermediate_root: Optional[Path] = None,
) -> Dict[str, List[RegulatorVariant]]:
    """Sample -> stage-6 regulator variants, from ``regulator_variants.tsv``.

    A missing file yields ``{}``. A present-but-empty file raises out of
    ``regulators.load_regulator_variants`` with a message that says the screen
    ran and found nothing; that is caught and read as ``{}``, keeping the
    loader's own sidecar-based distinction intact rather than flattening it.
    """
    from ..errors import DataContractError
    from .regulators import load_regulator_variants

    table = Path(path)
    if not table.is_file():
        return {}
    try:
        records = load_regulator_variants(config, table, intermediate_root)
    except DataContractError as exc:
        if not _present_but_empty("regulator_variants", table, exc):
            raise
        return {}
    by_sample: Dict[str, List[RegulatorVariant]] = {}
    for variant in records:
        by_sample.setdefault(str(variant.sample_id), []).append(variant)
    return by_sample


def load_structural_calls(
    config: PipelineConfig, path: Path
) -> Dict[str, List[StructuralVariant]]:
    """Sample -> structural variants, from ``07_structural_variants.tsv``.

    A missing file yields ``{}``; see :func:`load_amr_calls` for the rest.
    """
    from ..errors import DataContractError
    from .sv import load_structural_variants

    table = Path(path)
    if not table.is_file():
        return {}
    try:
        records = load_structural_variants(config, table)
    except DataContractError as exc:
        if not _present_but_empty("structural_variants", table, exc):
            raise
        return {}
    by_sample: Dict[str, List[StructuralVariant]] = {}
    for sv in records:
        by_sample.setdefault(str(sv.sample_id), []).append(sv)
    return by_sample


def load_annotated_gene_names(
    root: Path, manifest: SampleManifest
) -> Dict[str, List[str]]:
    """Sample -> gene names, from stage 2's per-isolate annotation tables.

    ``root`` is the *parent* of the ``annotation/`` directory, because that is
    where ``annotation.load_from_intermediate`` reads from and passing the
    directory itself would silently yield no records for every sample.
    """
    from .annotation import load_from_intermediate

    directory = Path(root)
    if not directory.is_dir():
        return {}
    records = load_from_intermediate(directory.parent, manifest)
    return {
        sample_id: sorted({r.gene_name for r in rows if r.gene_name})
        for sample_id, rows in records.items()
    }


def _run_real_inputs(
    config: PipelineConfig,
    manifest: SampleManifest,
    mode: RunMode,
    antibiotic: str,
    *,
    intermediate_root: Optional[Path],
):
    """Read all four inputs from disk and refuse if any carries nothing.

    Mirrors ``convergence._run_real`` and ``cooccurrence._run_real``: the REAL
    path re-reads rather than trusting the caller's mappings, because a
    standalone or Snakemake invocation has no in-memory stage result to
    inherit, and a stage that silently fell back to an empty argument would
    emit a table with no rows that reads as "no mechanism found".

    All four inputs are checked for **presence**, not for row count, and the
    distinction is the whole point. A missing table means the producing stage
    did not run, or its output was pruned, and a stage that computed over what
    happened to be there would emit a table whose emptiness reads as "no
    mechanism found". A table that is *present with no rows* means the producer
    ran and found nothing, which is a result, and refusing on it would throw
    away a real answer. `regulators.load_regulator_variants` draws the same line
    using the provenance sidecar; here there is no sidecar, so presence is what
    is available and it is what is used.

    Every input is still announced. A present-but-empty table is reported at
    DEBUG with its row count, because "the table is there and it is empty" is
    exactly the thing a reader of a sparse `05_mechanisms.tsv` cannot tell from
    "stage 5 never ran".
    """
    from .real_inputs import refuse_incomplete

    if intermediate_root is None:
        raise refuse_incomplete(
            "mechanisms",
            {
                "intermediate_root": (
                    "REAL mode reads its inputs from disk, and no intermediate "
                    "root was supplied, so there is nowhere to read them from"
                )
            },
            blocked_by={
                "intermediate_root": (
                    "a caller that passes intermediate_root - run.py, or the "
                    "Snakemake rule for this stage"
                )
            },
        )

    root = Path(intermediate_root)
    paths = real_input_paths(config, root, mode)
    mechanism_map = mechanism_map_for_antibiotic(config, antibiotic)

    absent = {
        name: path
        for name, path in paths.items()
        if not Path(path).exists()
    }
    if absent:
        raise refuse_incomplete(
            "mechanisms",
            {
                name: (
                    "no table at this path. An absent table means the "
                    "producing stage did not run or its output was pruned; a "
                    "table that is present with no rows is a result and is not "
                    "refused."
                )
                for name in absent
            },
            paths=paths,
            blocked_by=_MECHANISM_BLOCKERS,
        )

    amr_calls = load_amr_calls(Path(paths["amr_calls"]), antibiotic, mechanism_map)
    regulator_variants = load_regulator_calls(
        config, Path(paths["regulator_variants"]), root
    )
    structural_variants = load_structural_calls(
        config, Path(paths["structural_variants"])
    )
    annotated_genes = load_annotated_gene_names(
        Path(paths["annotated_genes"]), manifest
    )

    for name, loaded in (
        ("amr_calls", amr_calls),
        ("regulator_variants", regulator_variants),
        ("structural_variants", structural_variants),
        ("annotated_genes", annotated_genes),
    ):
        rows = sum(len(v) for v in loaded.values())
        if rows:
            LOGGER.debug("Stage 5 REAL input %s: %d rows", name, rows)
        else:
            LOGGER.info(
                "Stage 5 REAL input %s: the table is present and carries no "
                "rows. That is the producing stage having run and found "
                "nothing, not stage 5 having found nothing.",
                name,
            )

    return amr_calls, regulator_variants, structural_variants, annotated_genes


def dedupe(calls: Iterable[MechanismCall]) -> List[MechanismCall]:
    """Remove exact duplicate mechanism calls, preserving order."""
    seen = set()
    out: List[MechanismCall] = []
    for call in calls:
        key = (call.sample_id, call.antibiotic, call.determinant, call.mechanism)
        if key in seen:
            continue
        seen.add(key)
        out.append(call)
    return out
