"""Stage 13 - convergence analysis.

Question: does a resistance-associated determinant occur in **independent
phylogenetic backgrounds**, or is it simply a feature of one lineage?

Classification, evaluated in this order and documented so the categories are
mutually exclusive:

1. ``widespread_background`` - carried by more than
   ``convergence.widespread_fraction`` of the **analysed** samples. Too common
   to say anything about convergence. "Analysed" is
   :func:`analysed_sample_count`, never ``len(manifest)`` - see that function
   for why the two differ by two orders of magnitude.
2. ``rare_isolated`` - carried by at most ``convergence.rare_max_samples``
   samples. Too few to distinguish chance from convergence.
3. ``lineage_associated`` - all carriers fall in a single lineage. The
   determinant tracks the background, not the phenotype.
4. ``recurrent_convergent`` - present in at least
   ``convergence.min_independent_lineages`` independent lineages and
   therefore a candidate for convergent selection.
5. ``unknown`` - no lineage information available.

Every call is accompanied by the lineage distribution, so a reader can see
the basis for the category. Convergence is a *pattern in the data*, not a
demonstration of a shared selective cause, and the notes column says so.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

from ..config.loader import PipelineConfig
from ..io.tsv import write_tsv
from ..logging_utils import get_logger
from ..manifest import SampleManifest
from ..models import ConvergenceCall, ConvergenceCategory, GwasResult, RunMode

LOGGER = get_logger("stages.convergence")

CONVERGENCE_COLUMNS: Tuple[str, ...] = (
    "determinant",
    "independent_lineages",
    "branch_count",
    "distribution",
    "convergence_category",
)

UNKNOWN_LINEAGE = "unknown"


@dataclass(frozen=True)
class DeterminantCarriers:
    """A determinant and the samples carrying it."""

    determinant: str
    carriers: Tuple[str, ...]

    @property
    def n_carriers(self) -> int:
        return len(self.carriers)


def collect_carriers(
    determinants: Mapping[str, Sequence[str]],
) -> List[DeterminantCarriers]:
    """Build carrier records from sample -> determinants mappings."""
    inverted: Dict[str, List[str]] = {}
    for sample_id, names in determinants.items():
        for name in set(names):
            inverted.setdefault(name, []).append(sample_id)
    return [
        DeterminantCarriers(determinant=name, carriers=tuple(sorted(samples)))
        for name, samples in sorted(inverted.items())
    ]


def lineage_distribution(
    carriers: Sequence[str], lineages: Mapping[str, str]
) -> Dict[str, int]:
    """Count carriers per lineage."""
    counts: Dict[str, int] = {}
    for sample_id in carriers:
        lineage = lineages.get(sample_id, UNKNOWN_LINEAGE)
        counts[lineage] = counts.get(lineage, 0) + 1
    return dict(sorted(counts.items()))


def analysed_sample_count(lineages: Mapping[str, str]) -> int:
    """How many samples this stage actually analysed.

    **Not ``len(manifest)``, and the difference is not cosmetic.** The manifest
    is the whole PDC roster - 967 isolates - even on a run that analyses ten,
    because cohort membership is the manifest's job and an unprepared isolate
    must stay visible so it can refuse rather than vanish. ``widespread_fraction``
    is 0.75, so with 967 in the denominator a determinant has to be carried by
    726 isolates to read as background: on a ten-isolate run
    :func:`classify` could never return ``widespread_background`` at all, and
    the category that exists to stop an over-claim was unreachable. This is the
    same defect ``run._assemblies_analysed`` was written to fix for the
    report's N; ``tests/unit/seam_call_site_defaulted_allowlist.txt`` records it
    for this call site, and its entry is the one finding in that file marked
    fixed - in this function, not at the call site.

    The population is **the samples carrying a real lineage label**, which is
    what makes it the analysed set rather than a subset chosen to flatter the
    result:

    * a sample with no lineage is one whose determinants classify
      ``UNKNOWN``. Counting it in the denominator would divide by a population
      this stage cannot classify anything about, deflating every carrier
      fraction and pushing calls *away* from ``widespread_background`` - the
      same over-claim in the opposite direction.
    * a sample with a lineage and no determinant **is** counted. It was
      analysed and carries nothing, which is the observation that makes a
      determinant's frequency mean anything. Dropping those samples is how a
      determinant carried by every sequenced isolate reads as rare.

    Both directions therefore move toward the conservative category, which is
    the right direction for a stage whose output is a scientific claim.

    Returns:
        The number of samples with a non-sentinel lineage label. Zero when the
        crosswalk is empty, which makes the ``widespread_background`` test in
        :func:`classify` inoperative rather than dividing by zero - and in REAL
        that state is already a refusal (see :func:`_run_real`).
    """
    return len(
        {
            sample_id
            for sample_id, lineage in lineages.items()
            if lineage and lineage != UNKNOWN_LINEAGE
        }
    )


def classify(
    carriers: DeterminantCarriers,
    distribution: Mapping[str, int],
    config: PipelineConfig,
    n_samples: int,
) -> ConvergenceCategory:
    """Assign a convergence category. See the module docstring for ordering."""
    rules = config.convergence

    if n_samples and carriers.n_carriers / n_samples > rules.widespread_fraction:
        return ConvergenceCategory.WIDESPREAD_BACKGROUND
    if carriers.n_carriers <= rules.rare_max_samples:
        return ConvergenceCategory.RARE_ISOLATED

    known = {k: v for k, v in distribution.items() if k != UNKNOWN_LINEAGE}
    if not known:
        return ConvergenceCategory.UNKNOWN

    independent = len(known)
    if independent < rules.min_independent_lineages:
        return ConvergenceCategory.LINEAGE_ASSOCIATED
    return ConvergenceCategory.RECURRENT_CONVERGENT


def analyse(
    config: PipelineConfig,
    carriers: Sequence[DeterminantCarriers],
    lineages: Mapping[str, str],
    n_samples: int,
) -> List[ConvergenceCall]:
    """Classify every determinant."""
    calls: List[ConvergenceCall] = []
    for record in carriers:
        distribution = lineage_distribution(record.carriers, lineages)
        category = classify(record, distribution, config, n_samples)
        calls.append(
            ConvergenceCall(
                determinant=record.determinant,
                independent_lineages=len(
                    {k for k in distribution if k != UNKNOWN_LINEAGE}
                ),
                branch_count=sum(distribution.values()),
                distribution=distribution,
                convergence_category=category,
            )
        )

    by_category: Dict[str, int] = {}
    for call in calls:
        key = call.convergence_category.value
        by_category[key] = by_category.get(key, 0) + 1
    LOGGER.info(
        "Stage 13: %d determinants classified | %s",
        len(calls),
        ", ".join(f"{k}={v}" for k, v in sorted(by_category.items())) or "none",
    )
    return calls


def convergent_calls(
    calls: Sequence[ConvergenceCall],
) -> List[ConvergenceCall]:
    """Only the ``recurrent_convergent`` calls."""
    return [
        c for c in calls if c.convergence_category is ConvergenceCategory.RECURRENT_CONVERGENT
    ]


def supported_associations(
    calls: Sequence[ConvergenceCall],
    gwas_results: Sequence[GwasResult],
    config: PipelineConfig,
) -> List[ConvergenceCall]:
    """Determinants that are both GW-associated and phylogenetically robust.

    This is the only place in the pipeline that may set
    ``ClaimStatus.SUPPORTED``, and it requires three independent things to
    line up: a significant association, carriage in multiple lineages, and
    absence of the widespread-background explanation. Anything less stays at
    ``ASSOCIATED``.
    """
    threshold = config.gwas.significance_threshold
    hits = {
        r.feature
        for r in gwas_results
        if r.adjusted_p_value is not None and r.adjusted_p_value <= threshold
    }
    return [
        c
        for c in calls
        if c.convergence_category is ConvergenceCategory.RECURRENT_CONVERGENT
        and c.determinant in hits
    ]


def load_amr_names(
    path: Path, antibiotic: str, mechanism_map: Optional[Mapping[str, object]] = None
) -> Dict[str, List[str]]:
    """Sample -> acquired determinant names, read from stage 4's table.

    Reads the table from disk rather than taking stage 4's in-memory
    ``Dict[str, List[AmrDeterminant]]``. That coupling does not exist when this
    stage runs on its own or under Snakemake, where a rule depends on a *file*;
    and the pattern for reading a stage's own dependency is already
    ``stages.amr.load_amr_table``, which this calls rather than re-parsing, so
    the two spellings of the AMR schema cannot drift apart.

    Only ``AMR`` determinant types are carried: convergence is about acquired
    determinants appearing in independent backgrounds, and folding in
    chromosomal or intrinsic calls would answer a different question.
    """
    from ..execution.contracts import table_path
    from .amr import load_amr_table

    stage_dir = Path(path).parent
    amr_table = table_path(stage_dir, "amr")
    if not amr_table.is_file():
        # `load_amr_table` ends in `read_tsv`, which *raises* on a missing file.
        # Returning `{}` instead lets the caller's completeness check name
        # `amr_calls` as absent, which is the refusal a reader can act on; a
        # bare "TSV file not found" from three frames down says nothing about
        # which stage input is missing or what would produce it.
        return {}
    table = load_amr_table(
        amr_table,
        antibiotic=antibiotic,
        mechanism_map=mechanism_map,
    )
    names: Dict[str, List[str]] = {}
    for determinant in table:
        if str(determinant.determinant_type).upper() != "AMR":
            continue
        name = determinant.gene or determinant.determinant
        if not name:
            continue
        names.setdefault(str(determinant.sample_id), []).append(str(name))
    return {sample: sorted(set(values)) for sample, values in names.items()}


def load_variant_names(path: Path) -> Dict[str, List[str]]:
    """Sample -> regulator-locus variant names, read from stage 6's table.

    ``{gene}:{variant}`` is the spelling ``run.py`` uses for the in-memory path,
    so a REAL run and a TEST run classify the same strings.
    """
    from ..io.tsv import read_tsv

    table = Path(path)
    if not table.exists():
        return {}
    names: Dict[str, List[str]] = {}
    for row in read_tsv(table, required_columns=("sample_id",)):
        gene = row.get("gene")
        variant = row.get("variant")
        if not (gene and variant):
            continue
        names.setdefault(str(row["sample_id"]), []).append(f"{gene}:{variant}")
    return {sample: sorted(set(values)) for sample, values in names.items()}


#: Where each REAL input is read from, for the refusal message. The stage's
#: own table location comes from ``contracts`` so this cannot name a stale path.
def real_input_paths(
    config: PipelineConfig,
    intermediate_root: Path,
    mode: RunMode,
    phylogeny_dir: Optional[Path] = None,
) -> Dict[str, str]:
    """Where each REAL input is read from.

    ``phylogeny_dir`` defaults to the configured stage-9 location but is
    overridable, because stage 9 writes outside the intermediate root and a
    caller - or a test - may legitimately hold the tree elsewhere. Without the
    override the only way to point this stage at a tree is to write into the
    real results tree, which is how the first version of these tests came to
    leave a fixture behind in `results/`.
    """
    from ..execution.contracts import table_path

    stage_dir = Path(intermediate_root) / "stages"
    phylo = Path(phylogeny_dir) if phylogeny_dir else config.phylogeny_dir(mode)
    return {
        "amr_calls": str(table_path(stage_dir, "amr")),
        "regulator_variants": str(
            Path(intermediate_root) / "regulators" / "regulator_variants.tsv"
        ),
        "lineages": str(phylo / "tree_metadata.tsv"),
    }


def run(
    config: PipelineConfig,
    manifest: SampleManifest,
    mode: RunMode,
    amr_calls: Optional[Mapping[str, Sequence[str]]] = None,
    regulator_variants: Optional[Mapping[str, Sequence[str]]] = None,
    lineages: Optional[Mapping[str, str]] = None,
    n_samples: Optional[int] = None,
    *,
    intermediate_root: Optional[Path] = None,
    antibiotic: Optional[str] = None,
    phylogeny_dir: Optional[Path] = None,
) -> List[ConvergenceCall]:
    """Stage 13 entry point.

    **Where the inputs come from.** This module reads nothing from disk itself
    except through the ``load_*`` helpers, and only in REAL. In TEST the caller
    supplies in-memory mappings; in REAL every input is re-read from the
    contracted paths under ``intermediate_root`` (see :func:`real_input_paths`)
    and the three mapping arguments are **ignored**. That is deliberate: a
    standalone or Snakemake invocation of this stage has no upstream in-memory
    result to inherit, so the on-disk table is the only input there is.

    ``run.py`` may therefore pass ``amr_calls``, ``regulator_variants`` and
    ``lineages`` or omit all three; the outcome is identical in REAL. It must
    pass ``intermediate_root`` and ``antibiotic``, and it may pass
    ``phylogeny_dir`` and ``n_samples`` but need not - both have correct
    defaults (:func:`real_input_paths` and :func:`analysed_sample_count`).

    Args:
        manifest: cohort membership. **Not** the denominator for anything - see
            :func:`analysed_sample_count` - which is why it is accepted and
            unused rather than removed, since it is the orchestrator's stage
            signature.
        amr_calls: sample -> acquired determinant names. Ignored in REAL.
        regulator_variants: sample -> variant names in regulator loci.
            Ignored in REAL.
        lineages: sample -> lineage. Ignored in REAL.
        n_samples: the analysed sample count, when the caller knows it. Left
            out, it is derived from the crosswalk.
        intermediate_root: REAL only. Where the three inputs above are read
            from. Required in REAL; its absence is a refusal, not a fallback to
            the passed mappings.
        antibiotic: REAL only. Selects the AMR rows for the drug under
            analysis.
        phylogeny_dir: REAL only. Overrides ``config.phylogeny_dir(mode)``,
            because stage 9 writes outside the intermediate root.

    Returns:
        One :class:`~papipeline.models.ConvergenceCall` per determinant. The
        caller writes them with ``write_tsv`` to
        ``contracts.table_path(stage_dir, "convergence")``; this function
        writes nothing.

    Raises:
        StageError: In REAL, when any declared input carries no data. Nothing
            is written, so no partial table is left behind.
        NotImplementedError: In STUB, which fabricates this table instead.
    """
    if mode is not RunMode.TEST:
        if mode is RunMode.REAL:
            return _run_real(
                config,
                manifest,
                mode,
                amr_calls,
                regulator_variants,
                lineages,
                n_samples,
                intermediate_root=intermediate_root,
                antibiotic=antibiotic,
                phylogeny_dir=phylogeny_dir,
            )
        # STUB fabricates this table in `papipeline.stub`, so reaching the real
        # computation here would mean STUB stopped stubbing.
        raise NotImplementedError(
            "STUB-mode convergence is fabricated by papipeline.stub and has no "
            f"caller here. (mode={getattr(mode, 'value', mode)})"
        )

    combined: Dict[str, List[str]] = {}
    for source in (amr_calls or {}, regulator_variants or {}):
        for sample_id, names in source.items():
            bucket = combined.setdefault(sample_id, [])
            bucket.extend(names)

    carriers = collect_carriers(combined)
    return analyse(
        config,
        carriers,
        lineages or {},
        n_samples or analysed_sample_count(lineages or {}),
    )


def _run_real(
    config: PipelineConfig,
    manifest: SampleManifest,
    mode: RunMode,
    amr_calls: Optional[Mapping[str, Sequence[str]]],
    regulator_variants: Optional[Mapping[str, Sequence[str]]],
    lineages: Optional[Mapping[str, str]],
    n_samples: Optional[int],
    *,
    intermediate_root: Optional[Path],
    antibiotic: Optional[str],
    phylogeny_dir: Optional[Path] = None,
) -> List[ConvergenceCall]:
    """The REAL path: read from disk, then require every input to be present.

    Ordered deliberately - read, then check, then compute. Checking before
    reading would have nothing to check; computing before checking would leave
    a partial table behind when the check fires.
    """
    from .real_inputs import missing_required_inputs, refuse_incomplete

    if intermediate_root is None:
        raise refuse_incomplete(
            "convergence",
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
    paths = real_input_paths(config, root, mode, phylogeny_dir)
    drug = antibiotic or str(
        (getattr(config, "raw", None) or {}).get("project", {}).get(
            "primary_antibiotic"
        )
        or config.antibiotics[0]
    )

    amr_calls = load_amr_names(paths["amr_calls"], drug)
    regulator_variants = load_variant_names(Path(paths["regulator_variants"]))
    lineages = _load_lineages(Path(paths["lineages"]))

    missing = missing_required_inputs(
        {
            "amr_calls": amr_calls,
            "regulator_variants": regulator_variants,
            "lineages": lineages,
        },
        sentinels={"lineages": UNKNOWN_LINEAGE},
    )
    if missing:
        # The reason each is absent is not the same for all three, and the
        # reader needs to know which are merely unwritten and which are written
        # but degenerate.
        raise refuse_incomplete(
            "convergence",
            missing,
            paths=paths,
            blocked_by=_CONVERGENCE_BLOCKERS,
        )

    combined: Dict[str, List[str]] = {}
    for source in (amr_calls, regulator_variants):
        for sample_id, names in source.items():
            combined.setdefault(sample_id, []).extend(names)

    carriers = collect_carriers(combined)
    # The analysed count, never `len(manifest)`. A REAL manifest is the whole
    # PDC roster (967 isolates) whatever this run actually analysed, so the
    # `widespread_fraction` denominator would be the roster and
    # `classify` could never return WIDESPREAD_BACKGROUND - the category that
    # prevents an over-claim would be unreachable. See `analysed_sample_count`.
    analysed = n_samples or analysed_sample_count(lineages)
    LOGGER.info(
        "Stage 13: denominator is %d analysed sample(s)%s", analysed,
        "" if n_samples else " (derived; caller passed no n_samples)",
    )
    return analyse(config, carriers, lineages, analysed)


#: What holds each input back, for the refusal message. Kept beside the stage
#: that needs it so the two cannot disagree about why a run cannot proceed.
#:
#: These strings were rewritten in round 12 because they documented a state that
#: no longer exists. The previous text is quoted here because the *shape* of the
#: mistake is worth keeping visible: it named two real tool gaps and then went on
#: to assert three things that had since become false, so a reader who resolved
#: every gap the message offered would arrive at the same wall with no
#: explanation left.
_CONVERGENCE_BLOCKERS = {
    "amr_calls": "stage 4 (amr) writes it - stages/04_amr.tsv; run with "
    "--only amr first",
    # CORRECTED. The previous text said "stage 6 (variants) has no REAL
    # producer yet - ticket 14/15; only papipeline.testing.synthetic writes
    # this table". That was true when written and is false now:
    # `stages.regulators.produce_regulator_variants` exists and
    # `run.derive_regulator_table` (run.py:1808) calls it in REAL, writing
    # `<tool_output_root>/regulators/regulator_variants.tsv`. In REAL
    # `tool_output_root` IS `intermediate_root` (loader.py:880), so the path
    # this stage reads is the path the producer writes.
    "regulator_variants": (
        "stage 6 (variants) then stage 6a (regulators) - "
        "run.derive_regulator_table calls "
        "stages.regulators.produce_regulator_variants, which writes "
        "regulators/regulator_variants.tsv into this run's intermediate root"
    ),
    # CORRECTED. The previous text claimed, in order: that
    # `build_tree_outputs` "writes no lineage_label - it emits only sample_id,
    # tree_tip_label and source"; that "lineage_label has no REAL producer at
    # all; only papipeline.testing.synthetic writes it"; and that producing one
    # "needs a clade-assignment rule ... and config/science.yaml phylogeny:
    # defines ... no lineage rule at all". All three are now false. Ruling R4
    # defined a lineage as the MLST sequence type, `config/science.yaml` carries
    # `lineage.method: st`, `phylogeny.load_lineage_labels` reads stage 3's
    # table, and `build_tree_outputs` writes the five-column
    # `TREE_METADATA_COLUMNS` including `lineage_label`. Only gap 1 survives.
    "lineages": (
        "stage 9 (phylogeny), and specifically the core gene alignment it "
        "builds its tree from. Gap 1, the tools: `phylogeny.aligner` is "
        "panaroo, which the pipeline never invokes (adapters/panaroo.py has no "
        "subprocess), so `core_gene_alignment_filtered.aln` is "
        "operator-provisioned; the gubbins-masked copy that becomes "
        "core_snp_alignment.fasta is stage 8's job and the pinned gubbins "
        "osx-arm64 build segfaults. The earlier claim that this stage also had "
        "no lineage_label producer was false and has been removed: "
        "phylogeny.build_tree_outputs writes lineage_label from stage 3's MLST "
        "sequence type (lineage.method: st), and a crosswalk arriving without "
        "that column means an upstream stage failed rather than that this input "
        "is unobtainable"
    ),
}


def _load_lineages(path: Path) -> Dict[str, str]:
    """Stage 9's sample -> lineage crosswalk, read from disk.

    ``stages.phylogeny.load_tree_metadata`` is reused rather than re-parsed, so
    the two spellings of the tree-metadata schema cannot drift. Unlike that
    function it does **not** warn-and-return-empty on a missing file: here an
    absent crosswalk is a refusal (see `_run_real`), and silently returning
    ``{}`` would turn "no tree" into "every determinant is lineage-unknown".
    """
    if not Path(path).exists():
        return {}
    from .phylogeny import load_tree_metadata

    return load_tree_metadata(path)
