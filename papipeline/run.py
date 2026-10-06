"""Native orchestrator.

Executes the sixteen stages, writing each stage's outputs to disk and
recording provenance. This is the primary entry point and needs nothing
beyond the pinned Python dependencies.

The Snakemake workflow in ``workflow/`` is a thin wrapper that calls these
same stage functions, so the two execution paths cannot diverge in their
scientific behaviour.

Stage ordering
--------------
Stages are numbered as in the specification, but they are *executed* in
dependency order rather than numeric order. Stage 5 (mechanism
interpretation) consumes the variant calls from stage 6 and the structural
calls from stage 7, so the execution order is 1, 2, 3, 4, 6, 7, 5, 8, ...
The report lists stages in numeric order; the log lists them in execution
order. :data:`EXECUTION_ORDER` is the single source of truth for the
ordering, and the Snakefile mirrors it.

Mode gating
-----------
``TEST`` reads only from ``test_data/`` and writes results under
``results/test/``. ``REAL`` reads from ``data/`` and is refused unless
``runtime.allow_real_mode`` is explicitly set, which it is not during the
build phase. The two roots are different directories by construction, so
synthetic and real data cannot be mixed.
"""

from __future__ import annotations

import json
import platform
import shutil
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, FrozenSet, List, Mapping, Optional, Sequence, Set, Tuple

from . import __version__
from .adapters import UNPROBED_VERSION, detect_all, probe_on_demand
from .execution.contracts import internal_table_path, table_path
from .config.loader import PipelineConfig, load_config
from .errors import ModeNotAllowedError, PipelineError, StageError
from .io.tsv import read_tsv, write_tsv
from .knowledge import efflux_regulator_genes
from .logging_utils import configure_logging, get_logger
from . import reference, stub
from .adapters.bakta import genome_stem_for
from .adapters.panaroo import OUTPUT_DIRNAME as _PANAROO_OUTPUT_DIRNAME
from .adapters import tblastn as tblastn_adapter
from .adapters.oprd_locus import (
    LocusResolution,
    StructuralCall,
    reference_cds_nucleotides,
    reference_protein,
    resolve_isolate,
    structural_call_from_paths,
    write_tblastn_query,
)
from .assemblies import locate_assembly
from .execution.specs import BaktaOutputs
from .join import join_manifest_to_phenotype
from .manifest import (
    SampleManifest,
    discover_manifest,
    discover_pdc_manifest,
)
from .execution.semantics import resumable
from .observatory import wiring
from .models import RunMode
from .stages import amr as stage_amr
from .stages import annotation as stage_annotation
from .stages import cohort_variants as stage_cohort_variants
from .stages import cooccurrence as stage_cooccurrence
from .stages import convergence as stage_convergence
from .stages import gwas as stage_gwas
from .stages import gwas_features as stage_gwas_features
from .stages import integration as stage_integration
from .stages import mechanisms as stage_mechanisms
from .stages import mlst as stage_mlst
from .stages import pangenome as stage_pangenome
from .stages import phenotype as stage_phenotype
from .stages import phylogeny as stage_phylo
from .stages import similarity as stage_similarity
from .stages import regulators as stage_regulators
from .stages.regulators import MIN_LOCUS_COVERAGE_FRACTION as MIN_LOCUS_COVERAGE
from .stages import reporting as stage_reporting
from .stages.report_tables import load_run_tables
from .stages import sv as stage_sv
from .stages import validation as stage_validation
from .stages import variants as stage_variants
from .stages import virulence as stage_virulence
from .viz import (
    amr_heatmap_matrix,
    cooccurrence_network_table,
    convergence_table,
    determinant_distribution_table,
    efflux_regulator_by_phenotype_table,
    gene_by_phenotype_table,
    gwas_association_table,
    mechanism_by_phenotype_table,
    mechanism_distribution_table,
    mechanism_heatmap_matrix,
    oprd_by_phenotype_table,
    oprd_status_per_sample,
    phenotype_distribution_table,
    tree_annotation_table,
)

LOGGER = get_logger("run")

#: Stage name -> the config flag that enables it.
STAGE_FLAGS: Mapping[str, str] = {
    "mlst": "mlst",
    "amr": "amr",
    "virulence": "virulence",
    "pangenome": "pangenome",
    "phylogeny": "phylogeny",
    "gwas": "gwas",
    "convergence": "convergence",
    "cooccurrence": "cooccurrence",
    # Unbuilt (ticket 14). Off by default; asking for one is an error rather
    # than an empty table. See UNBUILT_STAGES.
    "variants": "variants",
    "cohort_variants": "cohort_variants",
    "recombination": "recombination",
    "similarity": "similarity",
}

#: Stages in spec.md D1 order, as reported.
#:
  #: Sixteen, per spec.md D1 (amended 2026-09-29 from fifteen).
#: - `mechanisms`, `regulators`, `structural_variants`, `integration` - are now
#: internal steps of the stage that owns them (spec.md:351), so they do not
#: appear here. Their code and tests are unchanged.
STAGE_ORDER: Tuple[str, ...] = (
    "validation",
    "annotation",
    "mlst",
    "amr",
    "virulence",
    "variants",
    "cohort_variants",
    "pangenome",
    "recombination",
    "phylogeny",
    "similarity",
    "phenotype",
    "gwas",
    "convergence",
    "cooccurrence",
    "reporting",
)

#: Stages in dependency order, as executed. Identical to STAGE_ORDER: the
#: internal steps folded into `amr`, `cooccurrence` and `reporting` run inside
#: their owner, in the order the owner already ran them.
EXECUTION_ORDER: Tuple[str, ...] = STAGE_ORDER

    #: Stages in the spec's sixteen that have no implementation.
#:
#: spec.md D1 lists them, so the DAG must list them too. But there is no code
#: behind them, and a stage that emitted a header-only table would be
#: indistinguishable from a stage that ran and found nothing - and "no sample
#: is similar to any other" is a claim, not an absence. So outside STUB they
#: refuse, naming the ticket that will build them.
#:
#: STUB is the exception, and necessarily so: fabricating a plausible output is
#: what that mode is for.
#:
#: `variants` and `cohort_variants` left this map in f5573f5, which built the
#: dispatch that finally called their parsers. What remains is genuinely
#: absent code, not a parser with no caller.
#: `similarity` left this mapping in 2223332, once it had a TEST path (1d54f31)
#: and a REAL caller (2223332). `variants` and `cohort_variants` left it earlier
#: for the same reason.
#: `recombination` left it in this reconciliation: the module was already
#: written (447 lines, a real refusal and a self-test) and had NO caller, which
#: is the half-built shape `tests/integration/test_stage_taxonomy_is_self_verifying.py`
#: was written to catch. Dropping the name is not the whole change - the
#: dispatch branch below is the other half, and `analysis.recombination` in
#: `config/science.yaml` is the third.
#:
#: **Empty on purpose.** Every stage in `STAGE_ORDER` is dispatched. A stage
#: added here without a dispatch branch is a stage nobody runs; a dispatch
#: branch added without removing the name here is a stage that runs while the
#: taxonomy says it is absent. `test_the_dispatch_covers_the_execution_order_exactly`
#: asserts the two agree, so neither half can be landed alone.
UNBUILT_STAGES: Mapping[str, str] = {}

#: Stages that are runnable in TEST and **refuse REAL on purpose**.
#:
#: This is a third state, and it is the one the taxonomy used to leave implicit.
#: `UNBUILT_STAGES` says a stage has no implementation. These have working TEST
#: implementations and a deliberate REAL refusal: without an injected engine
#: `gwas` would fall through to `ReferenceEngine` - Fisher's exact with BH and
#: no population-structure correction - which cannot detect the lineage
#: confounding a real cohort, so it refuses instead (f90d6c9, 46f5005). Same
#: reasoning for `convergence` and `cooccurrence` (46f5005).
#:
#: Reported separately because "runnable" cannot mean one thing and its opposite
#: at the same time. `tests/integration/test_stage_taxonomy_is_self_verifying.py`
#: asserts this set equals the set of stage modules whose `run()` refuses REAL -
#: whether the `NotImplementedError` is unconditional or guarded by the mode -
#: derived from their AST. So a stage cannot be added here without a refusal,
#: nor acquire one without being added here.
#:
#: Not to be confused with a stage whose `run()` raises `NotImplementedError`
#: because `runtime.allow_real_mode` is false. That is a gate on the run, not a
#: refusal of the stage: `similarity` does it and works once the gate is open.
REAL_REFUSING_STAGES: Mapping[str, str] = {
    "gwas": "no REAL engine: ReferenceEngine cannot handle lineage confounding",
    "convergence": "REAL convergence would be indistinguishable from TEST",
    "cooccurrence": "REAL co-occurrence cannot be tested on a synthetic cohort",
}

#: The ticket that will build them. Named in the refusal so the reader is told
#: what to do, not only what is missing.
UNBUILT_TICKET = "ticket 14"

#: What each stage needs in order to produce a meaningful result.
#:
#: Rewritten for the fold (spec.md:351): `regulators` is not a stage, so
#: `convergence` and `cooccurrence` no longer name it, and `mechanisms` is an
#: internal step of `cooccurrence` rather than a prerequisite of it.
PREREQUISITES: Mapping[str, FrozenSet[str]] = {
    "variants": frozenset({"amr"}),
    "cohort_variants": frozenset({"variants"}),
    "recombination": frozenset({"pangenome"}),
    "similarity": frozenset({"phylogeny"}),
    "gwas": frozenset({"phenotype", "variants", "pangenome"}),
    "convergence": frozenset({"amr", "phylogeny"}),
    "cooccurrence": frozenset({"amr", "variants", "annotation"}),
    "phylogeny": frozenset({"recombination"}),
    # What the report actually consumes, since it now also builds the master
    # table (spec.md:351). `--only reporting` has to bring all of this in, or
    # the run fails on a missing upstream stage.
    "reporting": frozenset({
        "validation", "phenotype", "amr", "cooccurrence", "gwas",
        "convergence", "phylogeny", "variants", "recombination", "similarity",
        # Listed so the stage is actually scheduled. Nothing analyses its output
        # yet - the contract is still open - but a stage nothing depends on never
        # runs at all, so a consumer is needed even a placeholder one. Reporting
        # is the honest one: a report that claims a complete run should have run
        # every stage, and this one refuses outside STUB rather than pretending.
        "cohort_variants",
    }),
}


def resolve_prerequisites(selected: Set[str]) -> Set[str]:
    """Expand a stage selection with its transitive prerequisites."""
    needed: Set[str] = set()
    frontier = set(selected)
    while frontier:
        stage = frontier.pop()
        for dependency in PREREQUISITES.get(stage, frozenset()):
            if dependency not in needed and dependency not in selected:
                needed.add(dependency)
                frontier.add(dependency)
    return needed


@dataclass
class RunResult:
    """Everything a completed run produced."""

    mode: RunMode
    antibiotic: str
    manifest: SampleManifest
    outputs: Dict[str, Path] = field(default_factory=dict)
    stage_status: Dict[str, str] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    master_table: List[Any] = field(default_factory=list)
    report_paths: Dict[str, Path] = field(default_factory=dict)
    figure_data: Dict[str, Any] = field(default_factory=dict)
    run_manifest: Optional[Path] = None
    #: What D5 verification found for the reference this run used, or None
    #: when the run is not against a reference. Recorded so the manifest
    #: states which genome the coordinates are about.
    reference: Optional[Mapping[str, Any]] = None

    @property
    def n_samples(self) -> int:
        return len(self.manifest)


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def resolve_mode(mode: str, config: PipelineConfig) -> RunMode:
    """Validate the requested run mode against configuration."""
    try:
        resolved = RunMode(str(mode).strip().upper())
    except ValueError:
        raise PipelineError(
            "Unknown run mode",
            mode=str(mode),
            allowed=",".join(m.value for m in RunMode),
        ) from None

    if resolved is RunMode.REAL and not bool(
        config.runtime.get("allow_real_mode", False)
    ):
        raise ModeNotAllowedError(
            "REAL mode is disabled: runtime.allow_real_mode is false in the "
            f"machine overlay for {config.machine.name!r}. To carry this cohort, "
            "run on the analysis machine instead - "
            "--machine bigmachine, whose overlay has no sample cap - or raise "
            "allow_real_mode there deliberately. Do not enable it on a machine "
            "that cannot hold the cohort: the cap exists because that machine's "
            "memory cannot. REAL is enabled only for the separate analysis "
            "phase, after this pipeline skeleton has been reviewed.",
            requested="REAL",
        )
    return resolved


def _run_one_stage(
    stage: str,
    execute: Any,
    *,
    obs: Any,
    stage_dir: Path,
    n_samples: int,
    tool_versions: Optional[Mapping[str, Any]] = None,
    prerequisites: Sequence[str] = (),
    classification: Any = None,
) -> None:
    """Run one stage, observed or not.

    With the observatory disabled this is a direct call inside the
    pipeline's own ``try``/``except``, exactly as the pipeline behaved
    before observability existed. With it enabled the same callable is
    handed to :func:`~papipeline.execution.run_task`, which executes it
    once and then persists state, validates the stage's declared outputs
    and publishes the transition.

    Either way the stage body runs exactly once, and a failure still stops
    the run. What the observatory adds is the record, not a new behaviour.
    """
    if not obs.enabled:
        try:
            execute(stage)
        except StageError:
            raise
        except Exception as exc:  # noqa: BLE001 - surfaced with stage context
            LOGGER.exception("Stage %s failed", stage)
            raise StageError(
                f"Stage failed: {stage}", stage=stage, error=str(exc)
            ) from exc
        return

    from .execution import RetryPolicy, StageState, TaskContext, run_task
    from .observatory.events import ExecutionEvent

    extra: List[Path] = []
    if stage == "phylogeny":
        companion = stage_dir / "10_alignment_summary.tsv"
        if companion.exists():
            extra.append(companion)

    version = None
    if tool_versions:
        status = tool_versions.get("bakta" if stage == "annotation" else stage)
        version = getattr(status, "version", None) if status is not None else None
        # A tool that is present but was never probed reports UNPROBED_VERSION,
        # and storing that string here would record it as a measured version.
        # `None` is what this code already wrote for an absent tool, and NULL in
        # a provenance row is visibly a gap rather than a plausible number. The
        # version that MATTERS - the one a per-genome Bakta task records - is
        # measured on demand by `_bakta_tool_version`, inside the stage that
        # decides to run the tool.
        if version == UNPROBED_VERSION:
            version = None

    inputs: List[str] = []
    if stage == "annotation":
        inputs = [str(p) for p in sorted((stage_dir.parent).glob("*.fna"))[:0]]

    ctx = TaskContext(
        run_key=obs.settings.run_key,
        stage=stage,
        subject="",
        spec=obs.spec_for(stage, stage_dir=stage_dir, n_samples=n_samples,
                           extra_paths=extra),
        tool_version=version,
        config_hash=obs.hash,
        input_ids=inputs,
        log_path=stage_dir.parent / "logs" / f"{stage}.log",
    )

    if prerequisites:
        obs.mark_unblocked(stage, prerequisites)

    result = run_task(
        ctx,
        lambda: execute(stage),
        store=obs.store,
        # One attempt: the pipeline has always stopped at the first stage
        # failure. Retrying a scientific stage here would be a behaviour
        # change dressed up as resilience.
        policy=RetryPolicy(max_attempts=1),
        event_sink=obs.event_sink,
    )

    if result.state is StageState.SUCCEEDED:
        return

    # The stage failed. Validation is tried anyway, because the reason a
    # stage fails is very often that it produced nothing, and saying
    # "INCOMPLETE: outputs absent" is more useful than "exit 1".
    detail = (
        result.validation.detail if result.validation is not None
        else (result.outcome.detail or f"state {result.state.value}")
    )
    LOGGER.error("Stage %s ended %s: %s", stage, result.state.value, detail)
    obs.bus.emit(ExecutionEvent(
        name="STAGE_FAILED", stage=stage, subject="",
        run_key=obs.settings.run_key, state=result.state.value,
        attempt=result.attempts, detail=detail,
        payload={"validation": json.loads(result.validation.to_json())
                  if result.validation else None},
    ))
    raise StageError(
        f"Stage failed: {stage}", stage=stage, error=detail,
    )


def discover_run_manifest(
    config: PipelineConfig,
    mode: RunMode,
    *,
    pdc_path: Optional[Path] = None,
) -> SampleManifest:
    """Discover the cohort for a run, honouring the smoke overlay.

    Two sources, chosen by the overlay rather than by the mode - both a smoke
    run and a full-cohort analysis are ``RunMode.REAL``, so the mode cannot
    distinguish them:

    * an overlay setting ``paths.smoke_genome_dir`` builds its cohort from
      ``PDC_essential.tsv`` and points the isolates it has at that directory.
      This is what makes the smoke overlay mean anything: it analyses the ten
      prepared assemblies, not whatever ``data/`` happens to hold;
    * every other overlay reads ``data/metadata/sample_metadata.tsv``, which is
      the behaviour it has always had and which this function does not change.

    ``PDC_essential.tsv`` is read directly and no ``sample_metadata.tsv`` is
    written. A populated copy would be a second statement of the cohort, free to
    drift from the first.

    Args:
        config: Loaded configuration; the overlay decides which source applies.
        mode: The run's mode.
        pdc_path: The PDC table. Defaults to the repository's copy.

    Returns:
        The manifest. In the smoke case every isolate in the PDC table is a
        member, and only those with a file in the smoke directory carry an
        ``assembly_path`` - so an isolate with no prepared assembly stays in the
        cohort and refuses per sample later, rather than silently disappearing.
    """
    if not config.is_smoke_overlay():
        return apply_cohort_subset(
            config,
            discover_manifest(
                config.metadata_dir(mode),
                pattern=config.sample_id_pattern,
                unique_attributes=config.sample_id_unique_attributes,
            ),
        )

    genome_dir = config.machine.smoke_genome_dir()
    if not genome_dir.is_dir():
        raise PipelineError(
            f"The smoke overlay points at {genome_dir}, which does not exist. "
            "There is no fall-back to `data/`: a bounded run that fell back "
            "would analyse the full cohort while appearing to be a 10-assembly "
            "one, which is the outcome this overlay exists to prevent.",
            genome_dir=str(genome_dir),
        )

    if pdc_path is None:
        pdc_path = pdc_table_path(config)
    pdc_path = Path(pdc_path)
    if not pdc_path.is_file():
        raise PipelineError(
            f"PDC isolate table not found at {pdc_path}. The smoke overlay "
            "builds its cohort from this file rather than from "
            f"`data/metadata/`, so it is required, and it is named as "
            f"`paths.pdc_table` in {pdc_table_config_source(config)} - an "
            "untracked clinical file (AGENTS.md rule 4: never committed), so a "
            "fresh checkout will not have one and every path has to be "
            "provisional until it is provisioned. This is the named input the "
            "Snakefile could not declare: under a smoke overlay no rule "
            "produces it and no other rule needs it, so `snakemake --dry-run` "
            "reported a complete DAG and then all sixteen stages died here.",
            pdc_path=str(pdc_path),
            config_key="paths.pdc_table",
        )

    assembly_paths = {
        fasta.stem: str(fasta)
        for fasta in sorted(genome_dir.glob("*.fna"))
    }
    manifest = apply_cohort_subset(config, discover_pdc_manifest(
        pdc_path, assembly_paths=assembly_paths
    ))
    prepared = sum(1 for s in manifest if s.assembly_path)
    LOGGER.info(
        "Smoke overlay: %d isolates in the cohort (after cohort.subset_file), "
        "%d with a prepared assembly in %s",
        len(manifest), prepared, genome_dir,
    )
    return manifest


#: The PDC isolate table's path, relative to the repository root.
#:
#: A constant and a default, not the only definition: the overlays name it under
#: ``paths.pdc_table`` and this is what an overlay that does not gets. It is a
#: *relative* name on purpose - the table is untracked clinical data that lives
#: beside the repository rather than inside it, so an absolute path here would
#: be correct only on the machine that wrote it.
PDC_TABLE_RELPATH = "PDC_essential.tsv"


def pdc_table_path(config: PipelineConfig) -> Path:
    """Where the PDC isolate table is, from configuration.

    **This is what makes the table a NAMED INPUT rather than a convention.**
    ``PDC_essential.tsv`` was required by ``discover_run_manifest`` and declared
    by nothing: not a rule, not an overlay key, not ``docs/data_contract.md``.
    Under a smoke overlay the cohort is built from this file and from nowhere
    else, so ``snakemake --dry-run`` reported a complete DAG and then every one
    of the sixteen stage invocations died at manifest load on the same missing
    file (E2E-DISCOVER D7). A required input with no name is an input no operator
    can provision, because nothing tells them it is one.

    So it is now a ``paths`` key with a default, read through the overlay like
    every other location, and the refusal above names the key. AGENTS.md rule 2
    forbids typing a path into code; this reads one from configuration and falls
    back to the historical relative name rather than to an absolute path.

    Args:
        config: The loaded configuration. ``machine`` may be ``None`` for a
            self-contained test config, in which case the repository root is
            used - the same root the overlays resolve against.

    Returns:
        The configured path. Not checked for existence: the refusal belongs to
        the caller, which knows whether this run needs the table at all.
    """
    declared = None
    if config.machine is not None:
        declared = config.machine.paths.get("pdc_table")
    relative = str(declared) if declared else PDC_TABLE_RELPATH
    path = Path(relative)
    return path if path.is_absolute() else (Path(config.root) / path)


def pdc_table_config_source(config: PipelineConfig) -> str:
    """Which overlay, if any, names the PDC table. For the refusal message."""
    if config.machine is None:
        return "the configuration"
    return str(getattr(config.machine, "source", "the machine overlay"))


def apply_cohort_subset(
    config: PipelineConfig, manifest: SampleManifest
) -> SampleManifest:
    """`cohort.subset_file`, applied: the cohort is only the ids it lists.

    **This is where the key stops being a comment.** `CohortConfig` parsed
    `subset_file`, resolved it against the repository root and read it - and
    nothing called `read_subset_file`, so the smoke overlay's ten isolates were
    *not* the cohort. `discover_run_manifest` built all 967 PDC isolates and the
    bounded run analysed every one of them, which is the exact outcome the
    overlay exists to prevent.

    Applied here, at manifest load, rather than in a stage: every stage reads
    `manifest`, so applying it once here bounds all of them, and applying it
    later would leave stage 1's validation counts describing a cohort the run
    did not analyse.

    Three decisions, each load-bearing:

    * **A listed id missing from the manifest is a refusal naming it.** Not a
      drop. A silently dropped isolate changes the cohort while the run reports
      the subset it was given, and the caller cannot tell which of the two
      happened.
    * **Manifest order is preserved; the list selects membership, not order.**
      Every downstream consumer - the phenotype join above all - iterates in
      manifest order, and reordering the cohort by an operator's file would make
      results depend on that file's line order for no stated reason.
    * **A duplicate in the list is a refusal.** `SampleManifest` would raise on
      the duplicate anyway, but naming it as a configuration error is more
      useful than "duplicate sample_id in the manifest".

    Returns the manifest unchanged when no subset is configured, which is the
    state of `laptop` and `bigmachine` - `science.yaml` says `null`.
    """
    members = config.cohort.read_subset_file()
    if members is None:
        return manifest

    path = config.cohort.resolve_subset_file()
    duplicated = sorted({m for m in members if members.count(m) > 1})
    if duplicated:
        raise PipelineError(
            f"cohort.subset_file lists {len(duplicated)} id(s) more than once: "
            f"{', '.join(duplicated[:10])}"
            + (" ..." if len(duplicated) > 10 else "")
            + f". Read from {path}. A cohort is a set of isolates; a repeated "
            "one is a configuration error, not a second cohort member.",
            subset_file=str(path),
            duplicated=",".join(duplicated[:10]),
        )

    index = manifest.index()
    missing = [member for member in members if member not in index]
    if missing:
        raise PipelineError(
            f"cohort.subset_file names {len(missing)} isolate(s) the manifest "
            f"does not contain: {', '.join(missing[:10])}"
            + (" ..." if len(missing) > 10 else "")
            + f". Read from {path}, which resolved to "
            f"{config.cohort.resolve_subset_file()}. The manifest decides "
            "cohort membership - "
            + (
                "for the smoke overlay that is PDC_essential.tsv"
                if config.is_smoke_overlay()
                else f"{config.metadata_dir(RunMode.REAL)}"
            )
            + " - so an isolate absent from it cannot be analysed. Dropping it "
            "instead would run a different cohort from the one the list names "
            "and still report the list, which is the substitution this refusal "
            "exists to prevent (AGENTS.md rule 5).",
            subset_file=str(path),
            missing=",".join(missing[:10]),
            n_missing=len(missing),
        )

    # Manifest order, not file order: see the docstring.
    return SampleManifest(
        [sample for sample in manifest if sample.sample_id in set(members)]
    )


def capped_sample_count(
    config: PipelineConfig, manifest: SampleManifest
) -> int:
    """The number `max_samples` is compared against, for this overlay.

    Two overlays, two quantities, and the difference is not cosmetic:

    * a **smoke** overlay declares ``max_samples`` in *assemblies* - it names
      the number of prepared genomes in its own directory. Its cohort is
      deliberately the whole of ``PDC_essential.tsv`` (all 967 isolates, every
      one a member), so comparing that against 10 would refuse every correctly
      wired smoke run before it started. The count is therefore the members
      that actually have a prepared sequence.
    * every other overlay declares it in *cohort members*, which is what it has
      always meant and always counted. A member with no assembly still counts,
      because the cap is about how much the machine can carry, and a member it
      cannot read is still work it was asked to do.

    Deliberately a *count*, not a second cap. ``enforce_sample_cap`` keeps its
    signature and its single raise site, so the smoke overlay header's claim
    that "there is no second cap mechanism" remains true - only the number fed
    to it differs.

    Note what this does **not** do: it does not shrink the cohort. Membership
    stays 967 and an isolate with no prepared assembly is still a member, which
    is what keeps ``merge_calls``'s denominator honest. What changes is which
    number the cap measures.
    """
    if config.is_smoke_overlay():
        return sum(1 for sample in manifest if sample.assembly_path)
    return len(manifest)


def run_pipeline(
    config_path: Optional[Path] = None,
    mode: str = "TEST",
    only: Optional[Sequence[str]] = None,
    skip: Optional[Sequence[str]] = None,
    write_html: bool = True,
    config: Optional[PipelineConfig] = None,
    observatory: Optional["wiring.Observatory"] = None,
    resume: bool = False,
    machine: "str | Path | None" = "laptop",
    oprd_structural_runner: Optional[Any] = None,
) -> RunResult:
    """Run the pipeline.

    Args:
        config_path: Path to ``config/science.yaml``. Defaults to
            ``<cwd>/config/science.yaml``.
        machine: Machine overlay - a name (``"laptop"``, ``"smoke"``) or a path
            to a ``machines/*.yaml`` file. Last in the signature so that adding
            it could not shift an existing positional argument. Ignored when
            ``config`` is passed, since the caller has already chosen.
        mode: ``TEST`` or ``REAL``.
        only: Run only these stages plus their prerequisites.
        skip: Skip these stages.
        write_html: Also emit the HTML report.
        config: Pre-loaded configuration, to avoid a second parse.
        observatory: Optional handle from
            :func:`papipeline.observatory.wiring.create`. When enabled, every
            stage runs *inside* :func:`~papipeline.execution.run_task`, so
            its state is persisted, its declared outputs are validated and
            its transitions are published to the observatory. When absent or
            disabled the stage is called directly and behaviour is identical
            to the pipeline before observability existed.
        resume: Continue an interrupted run rather than starting a new
            execution. With ``resume=False`` (the default) every stage
            executes, even if the run key already holds a SUCCEEDED record;
            that is a deliberate re-run, not a continuation. With
            ``resume=True``, stages whose recorded outputs still re-validate
            are skipped and only the rest are re-executed.
        oprd_structural_runner: Injected command runner for the two blast
            tools the oprD structural screen drives,
            ``runner(command=...) -> CommandResult``. ``None`` runs them for
            real. Present so a test can reach the REAL branch of this call site
            with synthetic inputs and faked subprocesses (R12) - which is the
            only way to test a call site that is otherwise reachable only under
            ``mode="REAL"``. Not a production parameter; nothing in
            ``scripts/`` passes it.

    Returns:
        A :class:`RunResult` with output paths and per-stage status.
    """
    if config is None:
        config = load_config(config_path, machine=machine)
    resolved = resolve_mode(mode, config)
    configure_logging(str(config.runtime.get("log_level", "INFO")))

    only_set: Optional[Set[str]] = set(only) if only else None
    skip_set: Set[str] = set(skip) if skip else set()

    unknown = (only_set or set()) | skip_set
    unknown -= set(STAGE_ORDER)
    if unknown:
        raise PipelineError(
            "Unknown stage name(s)",
            stages=",".join(sorted(unknown)),
            allowed=",".join(STAGE_ORDER),
        )

    if only_set is not None:
        effective = only_set | resolve_prerequisites(only_set)
    else:
        effective = set(STAGE_ORDER)

    def skip_reason(stage: str) -> Optional[str]:
        """Why this stage will not run, or ``None`` if it will.

        Split out from the decision so the reasons are *values*, not log lines.
        A reason that exists only as an ``LOGGER.info`` call is invisible to
        ``run_manifest.json``, which is the record a reader of a finished run has
        - and a stage that was asked for and did not run is then indistinguishable
        from a stage nobody asked for. That is exactly what happened to
        `recombination`: `analysis.recombination` was absent from science.yaml,
        so `stage_enabled` was False, this function returned False, and the run
        reported success with stage 9 building an unmasked tree.

        Returns ``None`` for a stage outside ``effective``: it was not asked for,
        which is not a skip and must not be recorded as one. ``--only gwas``
        leaving `similarity` alone is not an omission; `analysis.similarity`
        being false while `similarity` was requested is.
        """
        if stage not in effective:
            return None
        if stage in skip_set:
            return "named in --skip"
        if resolved is RunMode.STUB:
            # STUB fabricates every declared output, unbuilt stages included -
            # fabricating is the entire point of that mode - so neither the
            # unbuilt check nor the config flag may switch a stage off here.
            return None
        if stage in UNBUILT_STAGES:
            # Refuse *here*, before any work, and only when the stage would
            # actually run. Disabled and unrequested is a skip, like any
            # optional stage; asked for, or enabled, is an error. Returning
            # False instead would make `--only variants` look like it had run
            # and found nothing.
            #
            # Checked before the config flag, because a stage that is off by
            # default would otherwise skip silently and `--only variants` would
            # appear to have succeeded.
            wanted = only_set is not None and stage in only_set
            if wanted or config.stage_enabled(STAGE_FLAGS[stage]):
                raise NotImplementedError(
                    f"Stage {stage!r} is not implemented ({UNBUILT_TICKET}). It "
                    f"is one of the sixteen in spec.md D1, but there is no code "
                    f"behind it: {UNBUILT_STAGES[stage]}. Refusing rather than "
                    "writing an empty table, because an empty result is a "
                    "finding and this stage has not looked for anything. Run in "
                    "STUB mode to exercise the DAG."
                )
            return (
                f"unbuilt ({UNBUILT_TICKET}: {UNBUILT_STAGES[stage]}) and not "
                "requested"
            )
        flag = STAGE_FLAGS.get(stage)
        if flag and not config.stage_enabled(flag):
            # The key named, not just the fact. A reader who finds this in the
            # manifest needs to know which of the two gates closed, and this one
            # is the config flag rather than `runtime.allow_real_mode`.
            return f"disabled in configuration (analysis.{flag} is false)"
        return None

    #: Why each scheduled stage will not run. Computed up front, before any
    #: stage executes, so the refusal below can fire before the expensive ones
    #: have been paid for - the same ordering the sample cap uses.
    skipped_reasons: Dict[str, str] = {
        stage: reason
        for stage in EXECUTION_ORDER
        for reason in (skip_reason(stage),)
        if reason is not None
    }

    def enabled(stage: str) -> bool:
        if stage not in effective:
            return False
        reason = skipped_reasons.get(stage)
        if reason is None:
            return True
        # WARNING, not INFO. `--- stage X ---` is INFO, so a skip at INFO sits at
        # the same level as the lines it is meant to contrast with, and an
        # operator watching a long run at WARNING never sees it.
        LOGGER.warning("Stage %s skipped: %s", stage, reason)
        return False

    # A stage that something REQUESTED needs cannot be skipped. `--skip gwas` on
    # a full run used to produce a report that claimed a complete run without
    # ever running GWAS, which `PREREQUISITES["reporting"]` says in as many
    # words. Refusing names both stages, because naming only the skipped one
    # tells the reader nothing they did not already know.
    _requested = sorted(only_set) if only_set is not None else sorted(STAGE_ORDER)
    unmet: List[str] = []
    for stage in _requested:
        needed = sorted(
            resolve_prerequisites({stage}) | set(PREREQUISITES.get(stage, ()))
        )
        for upstream in needed:
            if upstream in skipped_reasons:
                unmet.append(
                    f"  - {stage!r} needs {upstream!r}: {skipped_reasons[upstream]}"
                )
    if unmet:
        raise PipelineError(
            "A stage this run was asked for depends on a stage that will not "
            f"run ({len(unmet)} unmet prerequisite(s)). Refusing rather than "
            "producing output that looks complete and is not:\n"
            + "\n".join(unmet)
            + "\n\nTo run it anyway, drop the --skip that causes it and let the "
            "upstream stage run; to omit the downstream stage too, name it in "
            "--skip as well, so the run says out loud that it is partial.",
            unmet="; ".join(unmet),
        )

    antibiotic = str(
        (config.raw.get("project") or {}).get("primary_antibiotic")
        or config.antibiotics[0]
    )
    config.require_antibiotic(antibiotic)

    # `assembly_root`, not `data_root`: under a smoke overlay a stage that
    # needs a sequence must read the overlay's bounded genome directory.
    # `data_root` put 602 unauthorised assemblies in front of a cohort of ten,
    # so stage 1 reported `pass=10` while stage 2 read 602.
    data_root = config.assembly_root(resolved)
    results_root = config.results_root(resolved)
    intermediate = config.intermediate_root(resolved)
    # Named after the accessor it calls, not after the concept. Several stages
    # declare this same parameter as `intermediate_root`, so a local called
    # `source_dir` made `grep source_dir` here and `grep intermediate_root` in
    # the stages look like a wiring disagreement when they were not - and stage
    # 7 was once handed `intermediate` instead of this, read its own outputs as
    # inputs, and reported a 0-gene pan-genome over a fixture declaring 10
    # (tests/unit/test_pangenome_caller_root.py).
    #
    # The two roots are NOT interchangeable. `intermediate_root` is where this
    # run writes; `tool_output_root` is where stage *input tables* are read from
    # (loader.py:537-540). In REAL they are the same path, which is what makes a
    # swap invisible there; in TEST they are `results/test/intermediate` and
    # `test_data/intermediate` respectively, so a swap reads a directory of this
    # run's own output. `tests/unit/test_stage_root_invariant.py` asserts the
    # mapping for every root-consuming stage so the swap cannot come back.
    tool_output_root = config.tool_output_root(resolved)
    reports_root = config.reports_root(resolved)

    reference_record: Optional[Mapping[str, Any]] = None
    if resolved is RunMode.REAL:
        # spec.md D5, and the only mode that calls variants against the
        # reference. Verified here, before any stage runs, because a truncated
        # or substituted genome shifts every coordinate silently - the run
        # would produce numbers and they would be about a different sequence.
        # STUB and TEST are excluded deliberately: neither calls variants, and
        # STUB must work with no reference on disk at all.
        reference_record = reference.verify_reference(config, resolved)

    if not data_root.exists():
        raise PipelineError(
            "Data root for this mode does not exist",
            mode=resolved.value,
            path=str(data_root),
        )

    for directory in (results_root, intermediate, reports_root):
        directory.mkdir(parents=True, exist_ok=True)

    stage_dir = intermediate / "stages"
    stage_dir.mkdir(parents=True, exist_ok=True)

    outputs: Dict[str, Path] = {}
    status: Dict[str, str] = {}
    warnings: List[str] = []
    tools = detect_all()

    unpinned = config.unpinned_references()
    if unpinned:
        warnings.append(
            f"{len(unpinned)} of {len(config.references)} references in "
            "config/references.tsv are unpinned. Every database version must be "
            "pinned and recorded before a real analysis (scientific rule 8)."
        )
    provenance = build_provenance(config, tools)

    if resolved is RunMode.STUB:
        # STUB fabricates every output and must not read a cohort, a fixture or
        # a manifest. Discovering one here would make a stub run depend on
        # exactly the files the mode exists to avoid needing.
        manifest = SampleManifest(samples=())
    else:
        # Overlay-driven, not mode-driven: the smoke overlay names its own
        # genome directory, and every other overlay reads `data/metadata/` as
        # it always has. See `discover_run_manifest`.
        manifest = discover_run_manifest(config, resolved)

    # The machine's cohort limit, checked as soon as the cohort is known and
    # before any stage runs. Enforcing this later - after annotation, say -
    # would be a refusal the user learns about hours too late, which is the
    # situation the cap exists to prevent.
    #
    # The count is overlay-dependent: a smoke overlay's `max_samples` is in
    # prepared assemblies, everyone else's is in cohort members. See
    # `capped_sample_count`; the cap itself is unchanged.
    config.enforce_sample_cap(capped_sample_count(config, manifest))

    result = RunResult(
        mode=resolved, antibiotic=antibiotic, manifest=manifest, outputs=outputs
    )
    result.reference = reference_record
    result.warnings = warnings

    LOGGER.info(
        "=== %s mode | %d samples | antibiotic=%s | data=%s ===",
        resolved.value,
        len(manifest),
        antibiotic,
        data_root,
    )

    def mark(name: str, state: str = "completed") -> None:
        status[name] = state

    def record(key: str, path: Path) -> None:
        outputs[key] = Path(path)

    # Mutable state shared between stages.
    qc: List[Any] = []
    qc_summary: Dict[str, Any] = {}
    annotations: Dict[str, List[Any]] = {}
    annotated_gene_names: Dict[str, List[str]] = {}
    mlst_calls: Dict[str, Optional[Any]] = {}
    amr_calls: Dict[str, List[Any]] = {}
    regulator_variants: Dict[str, List[Any]] = {}
    #: Per-isolate oprD locus verdicts from `resolve_oprd_loci`, or None
    #: outside REAL. None and {} are different here: `oprd_status_per_sample`
    #: reads None as "fall back to the annotation symbol" and any mapping
    #: as "these are the verdicts".
    oprd_locus_resolution: Optional[Dict[str, Any]] = None
    #: Per-isolate oprD reading-frame verdicts from `resolve_oprd_structure`,
    #: or None outside REAL. Kept apart from `oprd_locus_resolution` because
    #: `viz.oprd_status_per_sample` treats a structural verdict as authoritative,
    #: so merging a structural `not_assessed` into the carriage mapping would
    #: take every sample from `intact` to `not_assessed` and add nothing.
    oprd_structural_calls: Optional[Dict[str, StructuralCall]] = None
    #: Per-isolate calls from stage 6, keyed by isolate. Handed to
    #: `cohort_variants` (6a) so the merge consumes this run's calls rather than
    #: re-reading a file that could disagree with them.
    variant_calls: Dict[str, List[Any]] = {}
    #: Stage 6's per-isolate verdict, same handoff for the same reason. The
    #: calls alone cannot express it: an isolate that produced no variant and
    #: one that could not be read are both an empty list. Passed in memory so
    #: `an_calls` in the merge reflects THIS run, not a sidecar that could be
    #: stale.
    variant_statuses: Dict[str, str] = {}
    structural: Dict[str, List[Any]] = {}
    mechanism_calls: Dict[str, List[Any]] = {}
    virulence: Dict[str, List[Any]] = {}
    lineages: Dict[str, str] = {}
    tree_summary = None
    alignment_info: Dict[str, Any] = {}
    phenotype_calls: List[Any] = []
    phenotype_by_sample: Dict[str, Any] = {}
    gwas_results: List[Any] = []
    gwas_input: Optional[Any] = None
    convergence_calls: List[Any] = []
    cooccurrence: List[Any] = []
    master: List[Any] = []

    # Ownership: a handle passed in belongs to the caller, who may still be
    # reading from it after the run. Only a handle created here is closed
    # here.
    owns_observatory = observatory is None
    obs = observatory if observatory is not None else wiring.create(config)
    from .execution.semantics import Intent, RunKind, classify

    # What kind of execution is this? Decided from the store before any
    # stage runs, so a resume knows which stages it may skip.
    intent = Intent.RESUME if resume else Intent.EXECUTE
    verdict = classify(obs.store if obs.enabled else None,
                       obs.settings.run_key, intent=intent)
    obs.classification = verdict
    LOGGER.info("Execution kind: %s (%s)", verdict.kind.value, verdict.reason)
    executed: List[str] = []
    skipped: List[str] = []
    run_failed = False
    if obs.enabled:
        obs.announce_run(verdict, mode=resolved.value, stages=list(EXECUTION_ORDER))
    try:
        if obs.enabled:
            LOGGER.info(
                "Observatory enabled: run=%s store=%s",
                obs.settings.run_key, obs.settings.db_path,
            )

        def _execute_stage(stage: str) -> None:
            """Run one stage.

            The body below is the pipeline's own stage chain, unchanged and
            in the same order. It is a closure rather than inline code only
            so that :func:`papipeline.execution.run_task` can invoke it and
            the stage's state, validation and events are recorded around it.
            """
            nonlocal qc, qc_summary, annotations, annotated_gene_names, mlst_calls, amr_calls, regulator_variants, oprd_locus_resolution, oprd_structural_calls, variant_calls, variant_statuses, structural, mechanism_calls, virulence, lineages, tree_summary, alignment_info, phenotype_calls, phenotype_by_sample, gwas_results, gwas_input, convergence_calls, cooccurrence, master

            if resolved is RunMode.STUB:
                # STUB runs no parser, no classifier and no tool. It writes the
                # outputs the contract declares and stops, so the workflow's
                # shape, the event stream and the dashboard are all exercised
                # without anything real being read or computed.
                for key, path in stub.fabricate(stage, stage_dir).items():
                    record(key, path)
                if stage == "phylogeny":
                    # Stage 9 writes two files that are not its contract table:
                    # the tree and the tip crosswalk, into `phylogeny_dir` rather
                    # than `stage_dir`. `workflow/Snakefile` declares both as
                    # OUTPUTS of `rule phylogeny` - they used to be declared
                    # externals, which asked the DAG to resolve as inputs the
                    # files the rule itself produces - so a STUB run has to write
                    # them too or Snakemake stops on a missing output.
                    #
                    # Fabricated, not analysed: the tree below is a single
                    # unlabelled node and the crosswalk has a header and no
                    # rows. STUB exists to prove the DAG's shape with nothing
                    # real on disk, and a tip with a label would be a fabricated
                    # result wearing a result's shape.
                    for key, path in fabricate_phylogeny_tree_outputs(
                        config, resolved
                    ).items():
                        record(key, path)
                if stage == "reporting":
                    # Assigned to the result rather than `record`ed: the report
                    # is a pair of paths, not a stage table, and `record` maps
                    # one key to one Path.
                    result.report_paths = stub.fabricate_report(
                        reports_root, stub.stub_report_name(config, resolved)
                    )
                mark(stage)
                return

            if stage == "validation":
                qc = stage_validation.run(config, manifest, resolved)
                record(
                    "validation",
                    write_tsv(
                        table_path(stage_dir, "validation"),
                        stage_validation.rows(qc),
                        list(stage_validation.QC_COLUMNS),
                    ),
                )
                qc_summary = stage_validation.summarise(qc)
                if qc_summary.get("n_flagged"):
                    warnings.append(
                        f"Stage 1 flagged {qc_summary['n_flagged']} of "
                        f"{qc_summary['n_samples']} assemblies. Flagged samples are "
                        "retained in the analysis, not excluded."
                    )
                mark("validation")

            elif stage == "annotation":
                annotations = stage_annotation.run(
                    config, manifest, resolved, tool_output_root,
                    # Forwarded, not omitted. `tool_version` and
                    # `database_version` reach every per-GENOME Bakta
                    # `TaskContext` (`adapters/bakta.py:218-224`), and this
                    # function already resolves both for its own stage-level
                    # context at `_run_one_stage` - so omitting them recorded
                    # `NULL` where the value was in hand. A NULL tool version in
                    # the observatory is a provenance gap, not a style choice.
                    tool_version=_bakta_tool_version(tools, config),
                    database_version=_bakta_database_version(config),
                    # The observatory handle, so the per-genome Bakta tasks are
                    # persisted and published like every other tool's. Both are
                    # `None` when the observatory is disabled, which is what the
                    # stage already handles: `run_task` returns early from
                    # `announce` when there is no sink.
                    store=obs.store,
                    event_sink=obs.event_sink,
                    run_key=obs.settings.run_key,
                )
                annotated_gene_names = {
                    sid: sorted({r.gene_name for r in records if r.gene_name})
                    for sid, records in annotations.items()
                }
                record(
                    "annotation",
                    write_tsv(
                        table_path(stage_dir, "annotation"),
                        [
                            {
                                "sample_id": sid,
                                "n_records": len(records),
                                "n_named_genes": len(annotated_gene_names.get(sid, [])),
                            }
                            for sid, records in sorted(annotations.items())
                        ],
                        ["sample_id", "n_records", "n_named_genes"],
                    ),
                )
                mark("annotation")

            elif stage == "mlst":
                mlst_calls = stage_mlst.run(
                    config, manifest, resolved, tool_output_root,
                    data_root=data_root,
                )
                record(
                    "mlst",
                    write_tsv(
                        table_path(stage_dir, "mlst"),
                        [c.to_row() for c in mlst_calls.values() if c],
                        list(stage_mlst.MLST_COLUMNS),
                    ),
                )
                mark("mlst")

            elif stage == "amr":
                amr_calls = stage_amr.run(
                    config, manifest, resolved, tool_output_root, antibiotic,
                    data_root=data_root,
                )
                record(
                    "amr",
                    write_tsv(
                        table_path(stage_dir, "amr"),
                        [d.to_row() for calls in amr_calls.values() for d in calls],
                        list(stage_amr.AMR_COLUMNS),
                    ),
                )
                # spec.md:351 - `structural_variants` is an internal step of
                # `amr`: it normalises a tool's SV table, so it belongs with the
                # other screen. sv.py needs no alignment, so it is not blocked
                # on ticket 14.
                structural = stage_sv.run(config, manifest, resolved, tool_output_root)
                record(
                    "structural_variants",
                    write_tsv(
                        internal_table_path(stage_dir, "structural_variants"),
                        [s.to_row() for calls in structural.values() for s in calls],
                        list(stage_sv.SV_COLUMNS),
                    ),
                )
                candidates = sum(
                    1
                    for calls in structural.values()
                    for s in calls
                    if s.call_status.value == "candidate"
                )
                not_assessable = sum(
                    1
                    for calls in structural.values()
                    for s in calls
                    if s.call_status.value == "not_assessable"
                )
                if candidates:
                    warnings.append(
                        f"Stage 4 retained {candidates} candidate structural "
                        "variants. Candidates are reported as candidates and are "
                        "never promoted to confirmed."
                    )
                if not_assessable:
                    warnings.append(
                        f"Stage 4 could not assess {not_assessable} regions. These "
                        "are recorded as not_assessable, which is not the same as "
                        "'no variant present'."
                    )
                mark("amr")

            elif stage == "variants":
                # spec.md:351 - `regulators` is an internal step of `variants`,
                # so this stage owns the chromosomal screen as well as the
                # calling. Built and tested long before the calling existed, and
                # preserved unchanged.
                #
                # The calls are held in `variant_calls` rather than discarded,
                # because `cohort_variants` below merges *this run's* calls. A
                # re-read here would be a second read that could disagree.
                variant_calls = stage_variants.run(
                    config, manifest, resolved,
                    data_root=data_root,
                    workdir=intermediate / VARIANTS_WORK_DIRNAME,
                    statuses=variant_statuses,
                )
                record(
                    "variants",
                    write_tsv(
                        table_path(stage_dir, "variants"),
                        [row for rows in variant_calls.values() for row in rows],
                        list(stage_variants.PER_ISOLATE_COLUMNS),
                    ),
                )
                # AFTER the call above, not before it. `derive_locus_coverage`
                # reads `<workdir>/<sample>/<sample>.sorted.bam`, and the only
                # producer of those files is `stage_variants.run` ->
                # `adapters.minimap2.call_isolate`. Measured any earlier in this
                # branch it reads alignments the same run has not written yet,
                # and on a clean intermediate root the raise is unconditional:
                # that is the stage 6 failure a REAL 10-isolate run hit, every
                # isolate reported as `n_missing`.
                #
                # It stays here - after the calling, before both consumers -
                # because both consumers need it and neither produces it. The
                # ordering defect was in this branch's body, not in the DAG, so
                # no Snakefile prerequisite changed.
                #
                # REAL only: in TEST the committed gubbins-shaped and
                # variant-shaped fixtures are the stage's inputs and no BAM
                # exists, so measuring coverage there would refuse for want of
                # a file the mode does not produce.
                locus_coverage: Optional[Dict[Tuple[str, str], float]] = (
                    derive_locus_coverage(
                        config,
                        manifest,
                        variants_workdir=intermediate / VARIANTS_WORK_DIRNAME,
                    )
                    if resolved is RunMode.REAL
                    else None
                )
                # The regulator screen's REAL input is the ON-DISK stage input
                # table (docs/data_contract.md:114), so in REAL it is produced
                # here from the table just written and then read back by the
                # stage below rather than handed over in this scope. TEST reads
                # the committed fixture and nothing is produced. `variant_calls`
                # is deliberately not reused: the screen's input is the file, and
                # `cohort_variants` below owns the in-memory copy.
                if resolved is not RunMode.TEST:
                    # Measured from the alignments stage 6 just wrote, and
                    # handed to BOTH halves of the screen: the producer below
                    # applies the threshold as it writes the table, and the
                    # reader re-applies it to whatever it reads back. Passing it
                    # only to one half is how a parameter gets accepted and
                    # ignored - which is exactly what happened while nothing
                    # supplied it at all.
                    # Both coverage arguments, and the threshold is the module's
                    # own constant rather than a number typed here: it is a
                    # scientific threshold, so AGENTS.md rule 2 forbids
                    # re-stating it in the orchestrator, and spelling it makes
                    # the value the screen actually applies visible at the one
                    # place a reader looks for it.
                    derive_regulator_table(
                        config,
                        manifest,
                        calls_table=table_path(stage_dir, "variants"),
                        tool_output_root=tool_output_root,
                        locus_coverage=locus_coverage,
                        min_locus_coverage=MIN_LOCUS_COVERAGE,
                    )
                regulator_variants = stage_regulators.run(
                    config, manifest, resolved, tool_output_root,
                    locus_coverage=locus_coverage,
                    min_locus_coverage=MIN_LOCUS_COVERAGE,
                )
                record(
                    "regulators",
                    write_tsv(
                        internal_table_path(stage_dir, "regulators"),
                        [v.to_row() for calls in regulator_variants.values() for v in calls],
                        list(stage_regulators.REGULATOR_COLUMNS),
                    ),
                )
                # oprD carriage, resolved by protein identity rather than by
                # Bakta's annotation symbol. Computed here because this is the
                # stage that owns the oprD call screen; consumed by `reporting`,
                # which is where the figure is built. `None` outside REAL, and
                # `None` is passed on as `None` rather than as `{}` - see
                # `resolve_oprd_loci`.
                oprd_locus_resolution = resolve_oprd_loci(
                    config,
                    manifest,
                    resolved,
                    genomes_dir=data_root,
                    intermediate_root=intermediate,
                )
                # The structural half of the same screen: the same reference CDS,
                # searched against the isolate's ASSEMBLY rather than its
                # proteome, so the reading frame is judged rather than the
                # protein's identity. Deliberately a separate mapping and a
                # separate table - see `resolve_oprd_structure`'s docstring for
                # why the verdicts are not merged into `oprd_locus_resolution`.
                oprd_structural_calls = resolve_oprd_structure(
                    config,
                    manifest,
                    resolved,
                    genomes_dir=data_root,
                    intermediate_root=intermediate,
                    runner=oprd_structural_runner,
                )
                if oprd_structural_calls:
                    record(
                        "oprd_structural_calls",
                        write_oprd_structural_table(
                            oprd_structural_calls,
                            intermediate
                            / OPRD_LOCUS_WORK_DIRNAME
                            / OPRD_STRUCTURAL_TABLE_NAME,
                        ),
                    )
                mark("variants")

            elif stage == "cohort_variants":
                # Cohort-wide, from the calls the stage above just produced.
                # The filter itself - the minor-allele-frequency bounds - lives
                # in `stage_cohort_variants.merge_calls`; it is not restated
                # here, because two implementations of it would be free to
                # disagree and nothing would catch it.
                cohort_rows = stage_cohort_variants.run(
                    config, manifest, variant_calls,
                    statuses=variant_statuses,
                )
                record(
                    "cohort_variants",
                    write_tsv(
                        table_path(stage_dir, "cohort_variants"),
                        cohort_rows,
                        list(stage_cohort_variants.MERGE_COLUMNS),
                    ),
                )
                mark("cohort_variants")

            elif stage == "virulence":
                virulence = stage_virulence.run(
                    config, manifest, resolved, tool_output_root
                )
                record(
                    "virulence",
                    write_tsv(
                        table_path(stage_dir, "virulence"),
                        [v.to_row() for calls in virulence.values() for v in calls],
                        list(stage_virulence.VIRULENCE_COLUMNS),
                    ),
                )
                mark("virulence")

            elif stage == "pangenome":
                pangenome = stage_pangenome.run(
                    config, manifest, resolved, tool_output_root, annotations or None
                )
                paths = stage_pangenome.write_outputs(pangenome, stage_dir)
                record("pangenome", paths["pangenome_summary"])
                # Stage 12's feature matrix, built HERE rather than in the
                # `gwas` branch, and the reason is ordering rather than taste:
                # the producer reads three tables - stage 7's gene matrix,
                # stage 4's determinant table and stage 6's regulator table -
                # and every one of them is written by a stage that runs before
                # stage 7. `gwas.run` itself refuses REAL before it reads
                # anything (`stages/gwas.py:2065`), so calling the producer from
                # there would mean the refusal fires first and the producer is
                # never reached in a REAL run at all.
                #
                # Written to `tool_output_root`, which is the directory
                # `gwas.run` receives as its `intermediate_root` and therefore
                # the exact path it reads
                # (`Path(intermediate_root)/"gwas"/GWAS_FEATURE_TABLE_NAME`).
                #
                # REAL only, and the same gate `derive_regulator_table` uses
                # eleven lines above with its reason restated here because
                # `gwas_features`'s own docstring says the opposite: it claims
                # the producer does not branch on mode. It cannot, because in
                # TEST the three sources are NOT the same tables - stage 7
                # writes its gene matrix into `results/test/intermediate/stages`
                # while `tool_output_root(TEST)` is `test_data/intermediate`,
                # where no `stages/` directory exists at all. So in TEST the
                # committed `test_data/intermediate/gwas/gwas_features.tsv` is
                # the stage input, and producing there would overwrite a
                # committed fixture with a table built from a different root.
                # STUB is excluded for the reason it always was: it fabricates
                # every declared output and runs no parser.
                if resolved is RunMode.REAL:
                    record(
                        "gwas_features",
                        stage_gwas_features.produce_gwas_features(
                            config, manifest, resolved, tool_output_root
                        ),
                    )
                mark("pangenome")

            elif stage == "recombination":
                # Stage 8. Between 7 and 9 per spec.md D1, and the branch that
                # makes stage 9's `core_snp_alignment.fasta` input have a
                # producer: in REAL, `derive_recombination_tables` runs gubbins
                # and copies the masked alignment there. Before this branch
                # existed, `recombination` was in UNBUILT_STAGES, so
                # `stage_enabled("recombination")` was False, `enabled()`
                # returned False, the stage logged at INFO and was skipped -
                # and stage 9 then built its tree from whatever alignment it
                # found, silently unmasked. That is the silent-skip shape, and
                # `analysis.recombination` was absent from science.yaml the
                # whole time, which is why it was silent.
                record(
                    "recombination",
                    derive_recombination_tables(
                        config,
                        manifest,
                        resolved,
                        intermediate_root=intermediate,
                        phylogeny_dir=config.phylogeny_dir(resolved),
                        tool_output_root=tool_output_root,
                        stage_dir=stage_dir,
                    ),
                )
                mark("recombination")

            elif stage == "phylogeny":
                phylo_dir = config.phylogeny_dir(resolved)
                tree_summary, alignment_info = stage_phylo.run(
                    config, manifest, resolved, phylo_dir
                )
                lineages = stage_phylo.load_tree_metadata(
                    phylo_dir / "tree_metadata.tsv"
                )
                row = tree_summary.to_row()
                record(
                    "phylogeny_summary",
                    write_tsv(
                        table_path(stage_dir, "phylogeny"),
                        [row],
                        list(row.keys()),
                    ),
                )
                if alignment_info:
                    record(
                        "alignment_summary",
                        write_tsv(
                            stage_dir / "10_alignment_summary.tsv",
                            [{"metric": k, "value": v} for k, v in alignment_info.items()],
                            ["metric", "value"],
                        ),
                    )
                mark("phylogeny")

            elif stage == "similarity":
                # Reconciled once TEST (1d54f31) and REAL (2223332) both
                # existed. Added at the same time as the drop from
                # UNBUILT_STAGES, because a stage declared built but never
                # dispatched fails its own output contract: the validator looks
                # for `similarity.tsv`, nothing writes it, and the run stops
                # with "file does not exist" - which is how the reconciliation
                # first failed, in 22 tests at once.
                #
                # `variants` set the precedent: dropping the name from
                # UNBUILT_STAGES is not the whole change.
                similarity_rows = stage_similarity.run(
                    config,
                    manifest,
                    resolved,
                    tree_path=config.phylogeny_dir(resolved) / "tree.nwk",
                    out_path=stage_dir / "similarity.tsv",
                )
                record(
                    "similarity",
                    write_tsv(
                        table_path(stage_dir, "similarity"),
                        similarity_rows,
                        list(similarity_rows[0].keys()) if similarity_rows else ["sample_id"],
                    ),
                )
                mark("similarity")

            elif stage == "phenotype":
                phenotype_calls = stage_phenotype.load_phenotype(
                    config,
                    config.phenotype_dir(resolved),
                    antibiotic,
                    sample_ids=manifest.sample_ids,
                )
                # The directional join governs the manifest edge: the manifest
                # decides what exists, the phenotype table is evidence about it
                # and may be larger. A genome with no phenotype row is a hard
                # failure here, not a warning - before this the run continued
                # with the sample silently missing.
                # require_trait=False: the GWAS model is still binary R-vs-S,
                # so what the run needs is the categorical call, not a measured
                # MIC. Demanding a trait here would drop every sample without
                # one - 17 of the 20 committed fixtures - on the way to a model
                # that would not have used the trait anyway. This flips to True
                # when the continuous model lands (ticket 16).
                cohort = join_manifest_to_phenotype(
                    manifest.sample_ids, phenotype_calls, require_trait=False
                )
                phenotype_calls = [sample.call for sample in cohort.joined]
                phenotype_by_sample = {
                    sample.sample_id: sample.call for sample in cohort.joined
                }
                record(
                    "phenotype_exclusions",
                    cohort.write_exclusions(
                        stage_dir / "11_phenotype_exclusions.tsv"
                    ),
                )
                if cohort.excluded:
                    warnings.append(
                        f"{len(cohort.excluded)} phenotype rows were excluded "
                        f"for samples not in the manifest: {cohort.counts.summary()}"
                    )
                record(
                    "phenotype",
                    write_tsv(
                        table_path(stage_dir, "phenotype"),
                        [c.to_row() for c in phenotype_calls],
                        [
                            "sample_id",
                            "antibiotic",
                            "phenotype",
                            "MIC",
                            "MIC_unit",
                            "zone_diameter",
                            "zone_unit",
                            "source",
                        ],
                    ),
                )
                mark("phenotype")

            elif stage == "gwas":
                if phenotype_calls:
                    gwas_results, gwas_input = stage_gwas.run(
                        config,
                        manifest,
                        resolved,
                        tool_output_root,
                        phenotype_calls,
                        lineages=lineages,
                        engine=_build_gwas_engine(
                            config, resolved, intermediate
                        ),
                    )
                    record(
                        "gwas",
                        write_tsv(
                            table_path(stage_dir, "gwas"),
                            [r.to_row() for r in gwas_results],
                            list(stage_gwas.GWAS_COLUMNS),
                        ),
                    )
                    if gwas_input is not None and gwas_input.n_excluded:
                        warnings.append(
                            f"Stage 12 excluded {gwas_input.n_excluded} samples whose "
                            "category is outside the configured outcome groups. "
                            "Excluded samples are counted, not reassigned."
                        )
                    mark("gwas")
                else:
                    mark("gwas", "skipped_no_phenotype")
                    warnings.append(
                        "Stage 12 was skipped: no usable phenotype data. GWAS is "
                        "never run on an imputed or absent outcome."
                    )

            elif stage == "convergence":
                amr_names = {
                    sid: [d.gene or d.determinant for d in calls if (d.gene or d.determinant)]
                    for sid, calls in amr_calls.items()
                }
                variant_names = {
                    sid: [f"{v.gene}:{v.variant}" for v in calls]
                    for sid, calls in regulator_variants.items()
                }
                convergence_calls = stage_convergence.run(
                    config, manifest, resolved, amr_names, variant_names, lineages,
                    intermediate_root=intermediate, antibiotic=antibiotic,
                    # The ANALYSED sample count, stated here rather than left to
                    # the stage's own fallback. `widespread_fraction` is 0.75,
                    # so the denominator decides whether a ubiquitous
                    # determinant reads as background or as a convergent-
                    # selection candidate: with the 967-member PDC roster in the
                    # denominator, a determinant carried by 8 of 9 analysed
                    # isolates scores 0.008 and is published
                    # `recurrent_convergent` - the inversion CONV measured. The
                    # stage derives the same number from the crosswalk when this
                    # is None, so passing it is belt-and-braces; passing it makes
                    # the denominator a decision at the call site rather than a
                    # fact of a fallback nobody chose.
                    n_samples=_assemblies_analysed(manifest, qc_summary),
                    # Stated for the same reason stage 14 states it: the config
                    # accessor is the writer's own, so passing the same value
                    # makes the two callers visibly agree.
                    phylogeny_dir=config.phylogeny_dir(resolved),
                )
                record(
                    "convergence",
                    write_tsv(
                        table_path(stage_dir, "convergence"),
                        [c.to_row() for c in convergence_calls],
                        list(stage_convergence.CONVERGENCE_COLUMNS),
                    ),
                )
                mark("convergence")

            elif stage == "cooccurrence":
                # spec.md:351 - `mechanisms` is an internal step of
                # `cooccurrence`: the determinant classes whose co-occurrence
                # is tested are the mechanism classes. The report only renders
                # them, so they have to exist before this stage consumes them.
                mechanism_calls = stage_mechanisms.run(
                    config,
                    manifest,
                    resolved,
                    antibiotic,
                    amr_calls,
                    regulator_variants=regulator_variants,
                    structural_variants=structural,
                    annotated_genes=annotated_gene_names or None,
                    # Stage 14's own input root, given to stage 5 as well. In
                    # REAL `mechanisms` RE-READS all four of its inputs from this
                    # path and ignores the arguments above (`mechanisms.py`
                    # `_run_real_inputs`), so omitting it made a REAL run of
                    # stage 5 alone refuse with "no intermediate root was
                    # supplied" while the identical value sat one call site
                    # below.
                    intermediate_root=intermediate,
                )
                record(
                    "mechanisms",
                    write_tsv(
                        internal_table_path(stage_dir, "mechanisms"),
                        [c.to_row() for calls in mechanism_calls.values() for c in calls],
                        list(stage_mechanisms.MECHANISM_COLUMNS),
                    ),
                )
                amr_names = {
                    sid: [d.gene or d.determinant for d in calls if (d.gene or d.determinant)]
                    for sid, calls in amr_calls.items()
                }
                variant_names = {
                    sid: [f"{v.gene}:{v.variant}" for v in calls]
                    for sid, calls in regulator_variants.items()
                }
                mech_names = {
                    sid: sorted({c.mechanism for c in calls})
                    for sid, calls in mechanism_calls.items()
                }
                cooccurrence = stage_cooccurrence.run(
                    config,
                    manifest,
                    resolved,
                    amr_names,
                    variant_names,
                    mech_names,
                    lineages=lineages,
                    intermediate_root=intermediate,
                    antibiotic=antibiotic,
                    # The writer's own accessor, spelled at the call site for
                    # the same reason stage 13's is: `cooccurrence.real_input_paths`
                    # overrides it because stage 9 writes outside the
                    # intermediate root, and passing the same value here makes
                    # the reader and the writer visibly agree on where that is.
                    phylogeny_dir=config.phylogeny_dir(resolved),
                )
                record(
                    "cooccurrence",
                    write_tsv(
                        table_path(stage_dir, "cooccurrence"),
                        [c.to_row() for c in cooccurrence],
                        list(stage_cooccurrence.COOCCURRENCE_COLUMNS),
                    ),
                )
                mark("cooccurrence")

            elif stage == "reporting":
                # spec.md:351 - `integration` is an internal step of `report`:
                # the master table exists in order to be reported, so it is
                # built here and this stage is its owner.
                master = stage_integration.run(
                    config,
                    manifest,
                    resolved,
                    antibiotic,
                    phenotype_by_sample,
                    amr_calls,
                    mechanism_calls,
                    regulator_variants,
                    structural,
                    mlst_calls,
                    lineages,
                    virulence,
                    gwas_results,
                    convergence_calls,
                )
                master_path = internal_table_path(stage_dir, "master_table")
                stage_integration.write_master_table(master, master_path)
                record("master_table", master_path)
                result.master_table = master
                # Marked before the context is built so the report's own
                # stage table does not list itself as not_run.
                mark("reporting")
                figure_data = _prepare_figure_data(
                    config,
                    manifest,
                    phenotype_by_sample,
                    amr_calls,
                    mechanism_calls,
                    regulator_variants,
                    annotated_gene_names,
                    lineages,
                    mlst_calls,
                    gwas_results,
                    convergence_calls,
                    cooccurrence,
                    oprd_locus_resolution,
                )
                result.figure_data = figure_data
                figure_manifest = _figure_manifest(figure_data)
                record(
                    "figure_data",
                    write_tsv(
                        table_path(stage_dir, "reporting"),
                        [
                            {"figure": name, "kind": kind, "n_items": n}
                            for name, (kind, n) in figure_manifest.items()
                        ],
                        ["figure", "kind", "n_items"],
                    ),
                )
                _write_figure_data(figure_data, stage_dir / "figure_data")

                context = _build_report_context(
                    config=config,
                    mode=resolved,
                    antibiotic=antibiotic,
                    manifest=manifest,
                    status=status,
                    warnings=warnings,
                    provenance=provenance,
                    figure_data=figure_data,
                    qc_summary=qc_summary,
                    alignment_info=alignment_info,
                    gwas_input=gwas_input,
                    convergence_calls=convergence_calls,
                    cooccurrence=cooccurrence,
                    phenotype_calls=phenotype_calls,
                    stage_refusals=skipped_reasons,
                    oprd_calls=(
                        list(oprd_structural_calls.values())
                        if oprd_structural_calls
                        else (
                            list(oprd_locus_resolution.values())
                            if oprd_locus_resolution
                            else None
                        )
                    ),
                )
                result.report_paths = stage_reporting.write_report(
                    context, reports_root, write_html=write_html
                )

        for stage in EXECUTION_ORDER:
            if not enabled(stage):
                continue

            # A resume may skip a stage whose outputs still re-validate. A
            # re-run never skips, which is what makes it a re-run.
            if verdict.kind is RunKind.RESUME and obs.enabled:
                decision = resumable(
                    obs.store, obs.settings.run_key, stage, "",
                    obs.spec_for(stage, stage_dir=stage_dir,
                                 n_samples=len(manifest)),
                    config_hash=obs.hash,
                )
                if decision.skip:
                    LOGGER.info(
                        "--- stage %s skipped: %s ---", stage, decision.reason)
                    skipped.append(stage)
                    obs.announce_skipped(stage, reason=decision.reason)
                    mark(stage, "skipped_validated")
                    continue
                LOGGER.info("--- stage %s to re-execute: %s ---",
                            stage, decision.reason)

            LOGGER.info("--- stage %s ---", stage)
            executed.append(stage)
            _run_one_stage(
                stage,
                _execute_stage,
                obs=obs,
                stage_dir=stage_dir,
                n_samples=len(manifest),
                tool_versions=tools,
                prerequisites=sorted(PREREQUISITES.get(stage, ())),
                classification=verdict,
            )
    except StageError as exc:
        run_failed = True
        failed_stage = str(getattr(exc, "context", {}).get("stage", "") or "") or None
        if obs.enabled:
            obs.announce_run_finished(
                verdict, executed=executed, skipped=skipped,
                failed=failed_stage or "unknown",
                detail=str(exc))
        raise
    finally:
        if obs.enabled and not run_failed:
            obs.announce_run_finished(
                verdict, executed=executed, skipped=skipped,
                detail=(f"{len(executed)} stage(s) executed, "
                        f"{len(skipped)} skipped as already validated"),
            )
        if obs.enabled and owns_observatory:
            obs.close()

    result.stage_status = status
    result.warnings = warnings
    result.run_manifest = _write_run_manifest(
        results_root=results_root,
        mode=resolved,
        antibiotic=antibiotic,
        status=status,
        outputs=outputs,
        tools=tools,
        provenance=provenance,
        skipped=skipped_reasons,
    )

    LOGGER.info(
        "=== %s mode complete: %d stages executed, %d samples ===",
        resolved.value,
        len(status),
        len(manifest),
    )
    return result


def _bakta_tool_version(
    tools: Dict[str, Any], config: PipelineConfig
) -> Optional[str]:
    """The Bakta version THIS RUN will annotate with, or ``None``.

    Read from the same map every other stage's version comes from, and named
    rather than inlined so the annotation call site states what it forwards
    instead of where the number comes from. ``None`` when Bakta is absent, which
    is the honest value: the stage's own preflight refuses in that case, and
    recording a version for a tool that was never found would be a fabrication.

    **The version is probed on demand, and only when this run is about to
    execute Bakta.** ``detect_all()`` reports presence without running anything,
    so the sweep's value for Bakta is :data:`UNPROBED_VERSION` by design. This
    function replaces it with a measured version in exactly the one case where
    measuring is legitimate:

    * ``annotation.reuse_tool_output == "require"`` - Bakta is NEVER invoked
      (``stages/annotation.py`` refuses the isolate instead), so it is NEVER
      probed. Running ``bakta --version`` here would execute the tool on a run
      configured specifically so the tool is not executed, which is the
      violation a REAL tripwire shim logged (pids 93671, 93679).
    * any other mode - stage 2 will invoke Bakta for at least one genome, so the
      version is about to be asked for by a real run of the tool, and recording
      it is what provenance means.

    ``"prefer"`` is treated as "will run Bakta", which is what it means: it is
    the mode that falls through to the tool whenever no verified output exists.
    """
    status = tools.get("bakta") if tools else None
    if status is None or not getattr(status, "available", False):
        return None
    if getattr(status, "version", None) != UNPROBED_VERSION:
        return getattr(status, "version", None)
    if config.reuse_tool_output() == "require":
        # Deliberately left unprobed, and the value forwarded is None rather
        # than UNPROBED_VERSION: this argument reaches the per-genome Bakta
        # TaskContext, and under `require` there is no Bakta task for it to
        # reach. A version string here would be a tool version for an execution
        # that never happened.
        return None
    return probe_on_demand("bakta", tools).version


def _bakta_database_version(config: PipelineConfig) -> Optional[str]:
    """The pinned Bakta database version from ``config/references.tsv``.

    Read from the reference table rather than re-derived, so the version the
    per-genome tasks record is the same one ``build_provenance`` puts in
    ``run_manifest.json`` - two readers of one pin, which is the whole point of
    the table. ``None`` when no Bakta reference is pinned, which leaves the
    field NULL and visible rather than filled with a guess.
    """
    for spec in config.references.values():
        if spec.tool == "bakta":
            return spec.database_version
    return None


def build_provenance(
    config: PipelineConfig, tools: Mapping[str, Any]
) -> List[Dict[str, Any]]:
    """Merge the pinned reference table with the versions actually detected.

    The configured ``tool_version`` wins when it is pinned; otherwise the
    detected version is recorded and ``version_status`` stays ``unpinned`` so
    the gap is visible in the provenance record.
    """
    rows: List[Dict[str, Any]] = []
    for spec in config.references.values():
        detected = tools.get(spec.tool) if spec.tool else None
        tool_version = spec.tool_version
        # A detected-but-unprobed status carries UNPROBED_VERSION, which is the
        # absence of a measurement, not a measurement. Overwriting the pin's own
        # `UNPINNED` with it would replace "we have not pinned this" with a
        # string that reads like a version, and would hide the very gap
        # `version_status` exists to make visible.
        if (
            spec.version_status != "pinned"
            and detected is not None
            and detected.version != UNPROBED_VERSION
        ):
            tool_version = detected.version or tool_version
        rows.append(
            {
                "reference_id": spec.reference_id,
                "tool": spec.tool,
                "tool_version": tool_version,
                "database": spec.database,
                "database_version": spec.database_version,
                "version_status": spec.version_status,
            }
        )
    return rows


#: The PAO1 oprD locus tag the query protein is translated from.
#:
#: A constant, not configuration, for the reason
#: :func:`~papipeline.adapters.oprd_locus.reference_protein` gives: it derives the
#: query from the two pinned reference files, so a query from anywhere else is
#: not comparable to them and every identity number becomes unreproducible.
OPRD_LOCUS_TAG = "PA0958"

#: Where the per-isolate blastp scratch lives, under the run's intermediate root.
#:
#: Layout, not environment, exactly as
#: :data:`~papipeline.adapters.gubbins.WORK_DIRNAME` works: the root itself comes
#: from configuration.
OPRD_LOCUS_WORK_DIRNAME = "oprd_locus"

#: Stage 12's pyseer work directory, under the run's own intermediate root.
#:
#: Named here, not inside the ``gwas`` branch, because the engine needs it at
#: construction time and the branch is the only place a reader looks. It is
#: spelled by joining a name onto a root the configuration owns, so it cannot
#: land outside the run's tree.
GWAS_WORK_DIRNAME = "gwas_work"


def _build_gwas_engine(
    config: PipelineConfig,
    mode: RunMode,
    intermediate_root: Path,
) -> Optional[Any]:
    """The engine stage 12 is given, or ``None`` outside REAL.

    **``None`` outside REAL, and ``None`` is what makes the refusal correct.**
    ``stages.gwas.run`` refuses every non-TEST mode whose engine is ``None`` or
    a :class:`~papipeline.stages.gwas.ReferenceEngine`, so handing it nothing
    in REAL keeps the stage refusing rather than letting it fall through to
    Fisher's exact with BH and no population-structure correction - the outcome
    that refusal exists to prevent. In TEST the stage reads the committed
    fixtures and supplies its own engine, so injecting one here would change no
    number and only add a pyseer lookup to a mode that must not run tools
    (spec.md D8).

    **Why the engine at all.** ``PyseerEngine`` is a complete class - ``--lmm``,
    kinship wiring, and the two-pass step-12a threshold reduction - and had no
    caller anywhere in pipeline code, so the whole of it was reachable only from
    tests. ``gwas.py``'s own guard says the remaining blocker is that the stage
    *flags* a lineage-confounded feature rather than removing the confounding
    from the model; that limitation is unchanged by this wiring and is still
    stated in the stage's refusal message, which fires for the tool being absent.

    The executable is resolved by ``config.machine.tool_candidates`` - ``PATH``
    first, then the overlay's ``tool_search_dirs`` - rather than by
    ``shutil.which`` alone, so an overlay that keeps pyseer outside ``PATH``
    reaches it. **A pyseer that cannot be resolved returns ``None``, so the
    stage refuses with its own message naming pyseer**, rather than this function
    inventing an engine around a missing binary.

    The core-SNP alignment is passed when stage 8 wrote one, because pyseer
    estimates relatedness far better from core SNPs than from a few hundred gene
    presence/absence columns (``gwas.py:1044-1050`` documents the coarseness).
    Its absence is not fatal: the engine falls back to the Gram matrix, and says
    so on ``kinship_source``.
    """
    if mode is not RunMode.REAL:
        return None

    # `tool_candidates` already documents that it returns an EMPTY tuple when
    # nothing resolves, so the emptiness is the check rather than an index that
    # could raise. A raise here would replace the stage's refusal - which names
    # pyseer and what to do about it - with an IndexError from here.
    candidates = (
        config.machine.tool_candidates("pyseer")
        if config.machine is not None else ()
    )
    if not candidates:
        return None
    executable = candidates[0]

    # Imported here rather than at module scope for the reason
    # `derive_recombination_tables` imports it locally: `tests/unit/
    # test_r13_stage8_stage9_handoff.py` resolves the module a
    # `stage_recombination.X` call reaches from the import that binds the alias
    # *inside the function*, so a module-level binding would make that class of
    # check unresolvable rather than merely redundant. Two local imports of one
    # module is a smaller cost than an unresolvable call site.
    from .stages import recombination as stage_recombination

    workdir = Path(intermediate_root) / GWAS_WORK_DIRNAME
    snp_alignment = (
        Path(config.phylogeny_dir(mode))
        / stage_recombination.SNP_ALIGNMENT_NAME
    )
    return stage_gwas.PyseerEngine(
        executable,
        workdir,
        snp_alignment=snp_alignment if snp_alignment.is_file() else None,
    )


#: Stage 8's input directory, named by panaroo's owner rather than spelled here.
#:
#: `adapters.panaroo.OUTPUT_DIRNAME`, imported under an alias for the same reason
#: `workflow/Snakefile` aliases it: stage 7 reads the presence table from that
#: directory and stage 8 reads the core alignment from it, so a second literal
#: here is a second answer to a question the adapter already answers.
PANGENOME_OUTPUT_DIRNAME = _PANAROO_OUTPUT_DIRNAME


def fabricate_phylogeny_tree_outputs(
    config: PipelineConfig, mode: RunMode
) -> Dict[str, Path]:
    """STUB: the two files stage 9 writes outside its contract table.

    ``tree.nwk`` and ``tree_metadata.tsv`` live in ``phylogeny_dir``, not
    ``stage_dir``, because ``stages.similarity`` and ``stages.gwas`` read them
    from there by name. ``workflow/Snakefile`` declares them as outputs of
    ``rule phylogeny``; when they were declared externals instead, the rule
    demanded its own products as inputs and the REAL DAG could not resolve.

    Both files are fabricated and say so. The tree is one unlabelled node: a
    named tip would be a fabricated taxon, and ``tree_metadata.tsv`` is a header
    with no rows rather than a header-plus-crosswalk, because a crosswalk over
    zero samples cannot be distinguished from a fabricated one by shape.
    """
    phylogeny_dir = Path(config.phylogeny_dir(mode))
    phylogeny_dir.mkdir(parents=True, exist_ok=True)
    tree = phylogeny_dir / "tree.nwk"
    tree.write_text(
        "()STUB_NOT_A_TREE:0.0;\n",
        encoding="utf-8",
    )
    crosswalk = write_tsv(
        phylogeny_dir / "tree_metadata.tsv",
        [],
        list(stage_phylo.REQUIRED_TREE_METADATA),
    )
    return {"tree": tree, "tree_metadata": crosswalk}


def derive_regulator_table(
    config: PipelineConfig,
    manifest: SampleManifest,
    *,
    calls_table: Path,
    tool_output_root: Path,
    locus_coverage: Optional[Mapping[Tuple[str, str], float]] = None,
    min_locus_coverage: float = MIN_LOCUS_COVERAGE,
) -> Path:
    """REAL: produce the regulator screen's stage input table, or refuse by name.

    **The input is the file.** ``docs/data_contract.md:114`` declares
    ``regulators/regulator_variants.tsv`` a *stage input table*, and its section
    says that in REAL the tools write those into the run's intermediate
    directory. ``regulators`` had no producer at that call site, so every REAL
    run refused at stage 6 on ``missing='calls_by_isolate'``. This reads the
    per-isolate calls the ``variants`` branch just wrote
    (``<intermediate>/stages/variants.tsv``), rebuilds the ``sample_id -> rows``
    mapping from that file, and hands it to
    :func:`~papipeline.stages.regulators.produce_regulator_variants`, which
    writes the contracted table. ``stage_regulators.run`` then reads the table
    back, so the value the rest of the pipeline receives is the one on disk.

    Why the file rather than the ``variant_calls`` mapping already in scope: a
    second read is normally the thing to avoid, but here the file *is* the
    declared input, and ``cohort_variants`` downstream keeps the in-memory copy
    for the merge. Two readers of one table, one of which is the contract.

    **Every manifest sample is seeded, including the ones with no rows.** An
    absent key means "never screened" and an empty list means "called, found
    nothing"; collapsing them is how an unanalysable genome becomes a reported
    negative. See ``stages.variants.run`` and
    ``stages.regulators._run_real``.

    Args:
        calls_table: This run's per-isolate call table,
            ``table_path(stage_dir, "variants")``.
        tool_output_root: The stage-input root the regulator table is
            contracted to. In REAL this is also ``intermediate_root``
            (``config/loader.py:529-543``), so the table lands in this run's own
            intermediate directory, which is where the contract says a REAL
            stage input table lives.

    Returns:
        The path written.

    Raises:
        StageError: A declared input is absent - the call table, or the pinned
            reference the screen intersects against - or the call table names a
            sample the manifest does not. Each refusal names what is missing,
            where it was expected, and what produces it, so a reader deep in the
            DAG can act rather than infer.

    Args (continued):
        locus_coverage: ``(sample_id, gene) -> covered fraction``, forwarded to
            the producer so a locus the aligner only partly read is reported
            ``not_assessed`` rather than clean. See
            :func:`derive_locus_coverage`, which measures it.
        min_locus_coverage: The coverage threshold below which a regulator
            variant is not assessed, forwarded unchanged to
            :func:`~papipeline.stages.regulators.produce_regulator_variants`.

            **This parameter existed as a call-site keyword with nothing to
            receive it.** ``run.py:1262`` passed ``min_locus_coverage`` here and
            this function did not declare it, so every REAL run of stage 6 raised
            ``TypeError: derive_regulator_table() got an unexpected keyword
            argument 'min_locus_coverage'`` at that line - after stages 1-5 had
            run and before the stage produced anything. The companion
            ``stage_regulators.run`` call passes the same keyword to a callee
            that does declare it, which is why only this half was broken and why
            no targeted test saw it: every test of this function called it
            directly, without the keyword. Found by driving the call site
            through ``run_pipeline`` with a synthetic REAL cohort
            (``tests/integration/test_run_pipeline_call_sites.py``), which is the
            only way to reach a REAL-only call site; reverting this parameter
            makes that test fail with the exact ``TypeError`` above. The
            threshold is the regulators module's own constant, re-exported here
            as ``MIN_LOCUS_COVERAGE`` rather than restated as a number.
    """
    table = Path(calls_table)
    gff_path = config.reference_gff()
    fasta_path = config.reference_fasta()

    absent = [
        f"  - {what}: not found at {path} ({producer})"
        for what, path, producer in (
            (
                "stage-6 per-isolate calls",
                table,
                "written by this branch: stage_variants.run -> "
                "table_path(stage_dir, 'variants')",
            ),
            (
                "pinned reference annotation",
                gff_path,
                f"the pinned reference's own annotation, resolved from "
                f"config.reference_fasta() = {fasta_path}",
            ),
            (
                "pinned reference sequence",
                fasta_path,
                "declared by the machine overlay's reference section",
            ),
        )
        if not Path(path).is_file()
    ]
    if absent:
        raise StageError(
            "REAL-mode stage 6 refuses to run the regulator screen: "
            f"{len(absent)} declared input(s) are absent. The screen intersects "
            "this run's calls with the pinned reference's CDS features, and with "
            "an input missing it would either report a cohort-wide clean screen "
            "or attribute variants to unrelated genes - both fabricated "
            "negatives on the study's central mechanisms. Missing:\n"
            + "\n".join(absent),
            stage="variants",
            missing=str(len(absent)),
        )

    rows = read_tsv(table, required_columns=stage_variants.PER_ISOLATE_COLUMNS)
    calls: Dict[str, List[Dict[str, Any]]] = {sid: [] for sid in manifest.sample_ids}
    foreign: List[str] = []
    for row in rows:
        sample_id = str(row["sample_id"])
        if sample_id not in calls:
            foreign.append(sample_id)
            continue
        calls[sample_id].append(dict(row))
    if foreign:
        # Rule 5: a sample-ID mismatch is a hard failure with a clear message,
        # never a silent drop or a fuzzy join.
        raise StageError(
            f"{len(set(foreign))} row(s) of the stage-6 call table name sample(s) "
            f"the manifest does not: {sorted(set(foreign))[:10]}. The manifest "
            "decides cohort membership; a call from an isolate it does not name "
            "cannot be screened, because screening it would add a member the "
            "study does not contain.",
            stage="variants",
            foreign=",".join(sorted(set(foreign))[:10]),
        )

    LOGGER.info(
        "Stage 6: regulator screen reads %s (%d row(s) over %d isolate(s))",
        table, len(rows), len(calls),
    )
    target = stage_regulators.regulator_table_path(tool_output_root)
    stage_regulators.produce_regulator_variants(
        config, manifest, Path(tool_output_root), calls_by_isolate=calls,
        locus_coverage=locus_coverage,
    )
    return target


def derive_recombination_tables(
    config: PipelineConfig,
    manifest: SampleManifest,
    mode: RunMode,
    *,
    intermediate_root: Path,
    phylogeny_dir: Path,
    tool_output_root: Path,
    stage_dir: Path,
) -> Path:
    """Stage 8: mask recombination with gubbins and write `recombination.tsv`.

    **What this stage's input is, and who makes it.** The input is a *core gene
    alignment*, which panaroo writes as
    ``<intermediate>/panaroo/core_gene_alignment_filtered.aln``. **The pipeline
    never runs panaroo** - `papipeline/adapters/panaroo.py` imports no
    `subprocess` and contains no invocation - so the alignment is
    operator-provisioned, and its absence is a refusal naming the command to
    run and the file to produce. That is the same seam stage 9 reads
    (`core_snp_alignment.fasta`), so a REAL run without panaroo has already
    been told at stage 7.

    **The two modes do different work, and neither fabricates anything.**

    * ``REAL`` runs gubbins. The binary is resolved from the overlay's
      ``paths.tool_search_dirs`` and then *proved* by
      :func:`~papipeline.adapters.gubbins.self_test`, which runs it against a
      5-taxon probe alignment before the real run writes anything. A binary
      that segfaults (exit 139 - the broken conda ``osx-arm64`` build) is
      refused by name, because a real run would otherwise fail the same way
      after filling a scratch directory with intermediates. gubbins' masked
      alignment is then copied to where stages 9 and 12 read it.
    * ``TEST`` runs the real parser and the real writer over committed
      gubbins-shaped output under the TEST tool-output root. No external tool,
      exactly as ``stages.similarity.run(tree_path=...)`` reads a fixture. The
      masked alignment is **not** copied in TEST: `phylogeny_dir` under TEST is
      the committed ``test_data/phylogeny/``, and writing into it would replace
      a fixture with this run's output.

    Args:
        intermediate_root: This run's intermediate root. In REAL it is also the
            root the operator's panaroo output is read from.
        phylogeny_dir: Where the masked alignment is written in REAL.
        tool_output_root: Where TEST reads its gubbins-shaped output from.
        stage_dir: Where ``recombination.tsv`` is written, in both modes.

    Returns:
        The path of ``recombination.tsv``.

    Raises:
        ModeNotAllowedError: REAL with ``runtime.allow_real_mode`` false. Raised
            by the adapter before anything else, so an operator whose gate is
            shut is not also told their alignment is missing.
        ToolNotAvailableError: no gubbins binary, no ``run_gubbins.py``, or a
            binary that exits 139 on the self-test.
        StageError: the panaroo core alignment is absent, naming the exact
            command and output path that would produce it.
    """
    from .adapters import gubbins as gubbins_adapter
    from .adapters import panaroo as panaroo_adapter
    from .stages import recombination as stage_recombination

    table = table_path(stage_dir, "recombination")

    if mode is RunMode.TEST:
        gubbins_dir = Path(tool_output_root) / gubbins_adapter.WORK_DIRNAME
        rows = stage_recombination.run(
            config,
            manifest,
            mode,
            gubbins_dir=gubbins_dir,
            out_path=table,
        )
        LOGGER.info(
            "Stage 8: %d node row(s) from %s", len(rows), gubbins_dir,
        )
        return table

    # The gate FIRST, before the alignment is complained about. `run_gubbins`
    # checks it again, and deliberately so: it is the adapter's own contract and
    # it must hold for a caller that does not come through here. But this
    # function has an input to check before it reaches the adapter, and an
    # operator whose gate is shut does not need to be told their alignment is
    # missing - that is not their problem, and fixing it would not unblock them.
    gubbins_adapter.require_real_gate(config)

    alignment = Path(intermediate_root) / PANGENOME_OUTPUT_DIRNAME / (
        gubbins_adapter.CORE_ALIGNMENT_BASENAME
    )
    if not alignment.is_file():
        raise StageError(
            "REAL-mode stage 8 (recombination) has no core gene alignment. "
            f"Expected it at {alignment}. "
            + panaroo_adapter.provision_instructions(
                Path(intermediate_root) / PANGENOME_OUTPUT_DIRNAME
            )
            + " A missing alignment is a missing input, not an absence of "
            "recombination: this stage has not looked for anything.",
            stage="recombination",
            alignment=str(alignment),
        )

    search_dirs = config.machine.resolved_tool_search_dirs()
    work_dir = Path(intermediate_root) / gubbins_adapter.WORK_DIRNAME
    result = gubbins_adapter.run_gubbins(
        alignment,
        work_dir,
        runner=gubbins_adapter.resolve_runner(search_dirs),
        binary=gubbins_adapter.resolve_binary(search_dirs),
        threads=int(config.runtime.get("threads", 1) or 1),
        search_dirs=search_dirs,
        config=config,
    )
    gubbins_adapter.write_provenance(result, work_dir / "provenance.json")
    stage_recombination.build_outputs(
        gubbins_dir=work_dir, phylogeny_dir=Path(phylogeny_dir),
    )
    rows = stage_recombination.block_rows(
        gubbins_adapter.parse_per_branch_statistics(result.per_branch_statistics),
        stage_recombination.mean_branch_lengths(
            Path(result.node_labelled_tree).read_text(encoding="utf-8")
        ),
    )
    write_tsv(table, rows, list(stage_recombination.BLOCK_COLUMNS))
    LOGGER.info(
        "Stage 8: %d node row(s), %d with recombination detected",
        len(rows), sum(r["recombination_detected"] for r in rows),
    )
    return table


#: The stage-6 scratch directory, relative to the run's intermediate root.
#:
#: Named here for the same reason every other path constant in this file is: the
#: work directory is passed to ``stage_variants.run`` at one call site and read
#: back at another, and two spellings of one directory are how the reader ends up
#: looking somewhere the writer never wrote.
VARIANTS_WORK_DIRNAME = "variants"

#: The suffix ``adapters.minimap2.call_isolate`` gives a sorted, indexed BAM.
#:
#: Spelled from the adapter's own format rather than reconstructed: the adapter
#: builds ``f"{sample_id}.sorted.bam"``, and this reads the same stem.
SORTED_BAM_SUFFIX = ".sorted.bam"


def derive_locus_coverage(
    config: PipelineConfig,
    manifest: SampleManifest,
    *,
    variants_workdir: Path,
) -> Dict[Tuple[str, str], float]:
    """Per ``(isolate, locus)`` aligned fraction, from stage 6's own BAMs.

    **This is the producer the regulator screen's coverage machinery never had.**
    ``stages.regulators`` takes ``locus_coverage`` and ``min_locus_coverage`` and
    uses them to suppress records at a locus the aligner only partly read,
    reporting it ``not_assessed`` with the named reason
    ``locus_coverage_below_minimum``. Nothing in the repository ever supplied
    one, so ``not_assessed_reason(None, 0.9)`` returned ``None``, no record was
    ever suppressed, and ``locus_coverage_rows`` reported ``coverage=None`` for
    every pair: the threshold existed, had nothing to measure, and the stage
    reported a clean screen over loci no aligner had reached. That is a
    fabricated negative on the study's central mechanisms, and it was silent.

    **Where the evidence comes from: the alignments stage 6 already produced.**
    Stage 6 writes ``<variants_workdir>/<sample_id>/<sample_id>.sorted.bam`` for
    every isolate it managed to align (``adapters.minimap2.call_isolate``), and
    the covered fraction over a locus is computed by
    ``stages.regulators.locus_coverage_fraction`` from exactly those records -
    counting primary, secondary **and supplementary** alignments, because on
    this project's own data 8 of 10 isolates reach oprD only on a supplementary
    record and excluding them would report most of the cohort as unseen.

    **A missing alignment is a refusal naming the file, never a zero.** A pair
    with no BAM at all is not "0% covered"; it is "never looked for", and the
    two must not share a value. So an absent BAM raises rather than contributing
    ``0.0`` - a fabricated zero would suppress the whole locus and report
    ``not_assessed`` for an isolate the aligner never even attempted, which
    would replace one fabrication with another.

    Args:
        variants_workdir: This run's stage-6 scratch directory,
            ``intermediate/variants``.

    Returns:
        ``(sample_id, gene) -> covered fraction`` for every manifest sample and
        every screenable locus.

    Raises:
        StageError: ``pysam`` is not importable, the pinned GFF is absent or
            names none of the configured locus tags, or any manifest sample has
            no BAM. Each refusal names what is missing and what produces it.
    """
    from .adapters import gff as gff_adapter
    from .stages import regulators as regulators_stage

    try:
        import pysam
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise StageError(
            "Stage 6 cannot measure per-locus alignment coverage: `pysam` is "
            "not importable. Stage 6 already requires it to convert minimap2's "
            "SAM into the sorted, indexed BAM this function reads "
            "(variants.samtools.convert_and_index_with is pinned to pysam), so "
            "its absence means stage 6 never produced the alignments either. "
            "This is not a clean screen at a locus nobody aligned to - it is "
            "the absence of the measurement.",
            stage="variants", missing="pysam", error=str(exc),
        ) from exc

    gff_path = config.reference_gff()
    specs = {
        gene: spec for gene, spec in dict(config.regulators).items()
        if spec.screenable
    }
    if not gff_path.is_file():
        raise StageError(
            "Stage 6 cannot measure per-locus alignment coverage: the pinned "
            f"reference annotation was not found at {gff_path}. Coverage is "
            "measured against that file's CDS intervals; without it there is "
            "no interval to measure, and reporting coverage=None would be "
            "indistinguishable from a locus no aligner reached.",
            stage="variants", gff=str(gff_path),
        )
    if not specs:
        raise StageError(
            "Stage 6 cannot measure per-locus alignment coverage: no locus in "
            "config/regulators.tsv carries a PAO1 locus_tag, so no locus has "
            "an interval to measure. Resolve the tags with "
            "scripts/resolve_locus_tags.py.",
            stage="variants", n_loci=len(dict(config.regulators)),
        )

    # The adapter raises naming the file and the tags it could not find; both
    # refusals already exist and are tested, so they are not restated here.
    intervals_by_tag = gff_adapter.load_gene_intervals(
        gff_path, [spec.locus_tag for spec in specs.values()]
    )
    loci = sorted(
        (
            (gene, intervals_by_tag[spec.locus_tag])
            for gene, spec in specs.items()
            if spec.locus_tag in intervals_by_tag
        )
    )

    root = Path(variants_workdir)
    coverage: Dict[Tuple[str, str], float] = {}
    missing: List[str] = []
    for sample_id in manifest.sample_ids:
        bam = root / sample_id / f"{sample_id}{SORTED_BAM_SUFFIX}"
        if not bam.is_file():
            missing.append(
                f"  - {sample_id}: no alignment at {bam} (written by stage 6, "
                "adapters.minimap2.call_isolate)"
            )
            continue
        with pysam.AlignmentFile(str(bam), "rb") as handle:
            for gene, interval in loci:
                coverage[(sample_id, gene)] = (
                    regulators_stage.locus_coverage_fraction(
                        handle.fetch(interval.contig, interval.start - 1, interval.end),
                        interval.start,
                        interval.end,
                    )
                )
    if missing:
        raise StageError(
            f"Stage 6 cannot measure per-locus alignment coverage: {len(missing)} "
            "manifest sample(s) have no alignment on disk. Coverage is read off "
            "the BAMs stage 6 wrote; a missing one is an isolate the aligner "
            "never searched, which is not the same as a locus with zero "
            "coverage and must not be recorded as one. Missing:\n"
            + "\n".join(missing),
            stage="variants", n_missing=len(missing),
        )
    return coverage


def _oprd_locus_inputs(
    manifest: SampleManifest,
    genomes_dir: Path,
    intermediate_root: Path,
) -> List[Tuple[str, Path, Path, Path]]:
    """Every declared oprD-locus input, per sample, or a named refusal.

    Declared rather than discovered. Four inputs per isolate, each with a path
    that comes from somewhere else in this pipeline rather than from a guess:

    * the **assembly**, via :func:`~papipeline.assemblies.locate_assembly` -
      the same lookup annotation uses, so the resolver cannot analyse a
      sequence the annotation stage did not;
    * the isolate's `.gff3` and `.faa`, in the stage 2 output directory, named
      the way :func:`~papipeline.execution.specs.bakta_spec` names them, so a
      rename of Bakta's output cannot silently strand this.

    Returns:
        One ``(sample_id, assembly, gff3, faa)`` tuple per manifest sample, in
        manifest order.

    Raises:
        StageError: naming every missing input at once, with its expected path.
            Named rather than defaulted because each of these has exactly one
            remedy and a reader who is not told which one cannot act.
    """
    missing: List[str] = []
    resolved: List[Tuple[str, Path, Path, Path]] = []
    for sample in manifest:
        try:
            assembly = locate_assembly(sample, Path(genomes_dir))
        except PipelineError as exc:
            missing.append(f"  - assembly for {sample.sample_id}: {exc}")
            continue
        outputs = BaktaOutputs(
            out_dir=Path(intermediate_root) / "bakta" / sample.sample_id,
            sample_id=sample.sample_id,
            genome_stem=genome_stem_for(sample),
        )
        # Spelled as `bakta_spec(require_faa=True)` spells it rather than as a
        # new property, so this and the annotation contract cannot name different
        # files. `require_faa` is off there today - nothing consumed the
        # proteome - which is exactly why the absence needs naming here.
        faa = outputs.out_dir / f"{outputs.genome_stem}.faa"
        absent = [
            f"  - {sample.sample_id}: {what} not found at {path} "
            f"(written by stage 2, Bakta)"
            for what, path in (("gff3", outputs.gff), ("faa", faa))
            if not Path(path).is_file()
        ]
        if absent:
            missing.extend(absent)
            continue
        resolved.append(
            (sample.sample_id, assembly, Path(outputs.gff), Path(faa))
        )

    if missing:
        raise StageError(
            "REAL-mode oprD locus resolution refuses to run: "
            f"{len(missing)} declared input(s) are absent. The locus resolver is "
            "what stands between the annotation symbol and the study's central "
            "negative, and an unrun resolver leaves the symbol test in charge - "
            "so this stops rather than reporting `absent` for isolates whose "
            "locus was never looked for. Missing:\n" + "\n".join(missing),
            stage="variants",
        )
    return resolved


def resolve_oprd_loci(
    config: PipelineConfig,
    manifest: SampleManifest,
    mode: RunMode,
    *,
    genomes_dir: Path,
    intermediate_root: Path,
) -> Optional[Dict[str, LocusResolution]]:
    """Resolve every isolate's oprD locus, or refuse by name.

    **Why this is wired here at all.** :func:`viz.oprd_status_per_sample` decides
    carriage from Bakta's ``gene=`` symbol when it is given no
    ``locus_resolution``, and on real data that test reports ``absent`` where
    the gene is present: across ten isolates, 34 of the 36 ``gene=oprD`` CDS
    features are OprD/OprP/OprQ paralogs, and 8 of 10 carry the true locus
    unlabelled. See :mod:`papipeline.adapters.oprd_locus`.

    **The returned mapping must cover every manifest sample, and that is the
    invariant this function exists to hold.** ``oprd_status_per_sample`` falls
    through to the symbol test for any sample absent from the mapping, so a
    mapping that silently omits a sample reintroduces the fabricated negative
    for exactly that sample. A sample the resolver could not resolve therefore
    gets a `LocusResolution` carrying a **refusal verdict**, which
    ``oprd_status_per_sample`` maps to ``not_assessed`` - never ``absent``.

    **A missing input is a refusal, not a resolution.** Nothing here invents a
    verdict to keep going: an absent assembly, an absent ``.faa``/``.gff3`` or a
    missing ``blastp`` raises by name, because a resolver that returned
    "unresolved" for an isolate it never searched would report the same
    ``not_assessed`` as one it searched and found nothing - collapsing "not
    assessed" into "looked for, not there".

    Returns:
        ``None`` outside REAL. ``None`` is not an empty mapping and the
        distinction is load-bearing: ``oprd_status_per_sample`` treats ``None``
        as "run the symbol test" and any mapping as "these are the verdicts",
        so an empty dict would leave the symbol test running for every sample.
    """
    if mode is not RunMode.REAL:
        return None

    program = shutil.which("blastp")
    if program is None:
        raise StageError(
            "REAL-mode oprD locus resolution refuses to run: `blastp` is not on "
            "PATH. The resolver decides oprD carriage by protein identity "
            "against the pinned PAO1 PA0958 protein, and blastp is what performs "
            "that search; without it there is no evidence to decide from. "
            "Install the NCBI BLAST+ `blastp` (AGENTS.md rule 1: verify the "
            "version with `blastp -version` rather than assuming one). This is "
            "not the same as an isolate that lacks oprD.",
            stage="variants",
        )

    query = reference_protein(
        config.reference_gff(), config.reference_fasta(), OPRD_LOCUS_TAG
    )
    work_root = Path(intermediate_root) / OPRD_LOCUS_WORK_DIRNAME
    threads = int(config.runtime.get("threads", 1) or 1)

    resolutions: Dict[str, LocusResolution] = {}
    for sample_id, _assembly, gff3, faa in _oprd_locus_inputs(
        manifest, genomes_dir, intermediate_root
    ):
        verdict = resolve_isolate(
            sample_id,
            faa,
            gff3,
            query,
            workdir=work_root / sample_id,
            program=program,
            threads=threads,
        )
        LOGGER.info(
            "oprD locus %s: %s (coverage %.1f%%)",
            sample_id, verdict.verdict, verdict.coverage_pct,
        )
        resolutions[sample_id] = verdict
    return resolutions


#: The pinned oprD CDS, written once per run under the locus work root.
#:
#: A layout constant, not a filename spelled at the call site, because the
#: structural screen and the blastp resolver read the same directory and two
#: spellings of one filename are how the two drift apart.
OPRD_REFERENCE_CDS_NAME = "reference_oprd_cds.fna"

#: Where one isolate's tblastn artefacts live, under the locus work root.
#:
#: A subdirectory rather than the work root itself, and the reason is load
#: bearing: ``resolve_isolate`` writes its blastp query to
#: ``oprd_locus/<sample_id>/reference_protein.faa``, and the tblastn query has
#: the same filename with **different bytes** - the blastp query is the 443-residue
#: protein, the tblastn query is that plus the reference's own terminal stop
#: (444 residues, ``oprd_locus.reference_codon_count``). Writing both to one path
#: would mean the second silently destroys the auditable record of what blastp
#: searched. Separate directories, same filename.
OPRD_STRUCTURAL_DIRNAME = "tblastn"

#: The per-isolate structural verdicts table.
#:
#: **No contracted path exists for this**, which is a finding rather than a
#: choice: ``contracts.py`` declares no oprD table and the Snakefile declares no
#: oprD output, so the verdicts reach the report through
#: ``ReportContext.oprd_calls`` (see REPORT.md §2.5). It is written here anyway,
#: because a verdict that exists only in memory is not a record and this is the
#: study's central negative. Filed to BACKLOG with the contracted-path request.
OPRD_STRUCTURAL_TABLE_NAME = "structural_calls.tsv"

#: The columns of that table. Declared once, next to the writer.
OPRD_STRUCTURAL_COLUMNS: Tuple[str, ...] = (
    "sample_id",
    "structural_verdict",
    "structural_reason",
    "identity_pct",
    "coverage_pct",
    "orf_aa_length",
    "reference_aa_length",
    "orf_length_fraction",
    "lesion_type",
    "truncation_aa",
    "internal_stop_codon",
    "internal_stop_position_aa",
    "stop_is_reference_terminator",
    "frameshift",
    "split_abutting_fragments",
    "at_contig_edge",
    "contig",
    "strand",
    "footprint_start",
    "footprint_end",
    "n_hsp",
    "search_recorded",
    "search_command",
    "search_database",
    "search_parameters",
    "search_n_hits",
)


def resolve_oprd_structure(
    config: PipelineConfig,
    manifest: SampleManifest,
    mode: RunMode,
    *,
    genomes_dir: Path,
    intermediate_root: Path,
    runner: Optional[Any] = None,
) -> Optional[Dict[str, StructuralCall]]:
    """Search every isolate's assembly for the reference oprD CDS, and judge it.

    **This is the second half of the oprD locus screen, and it is the half that
    measures structure.** :func:`resolve_oprd_loci` decides *whether* the locus
    is present by protein identity against the isolate's proteome;
    :func:`~papipeline.adapters.oprd_locus.structural_call_from_paths` decides
    whether the reading frame is intact, from a translated search of the
    isolate's assembled nucleotides. Round 12's WIRE.md A1 established that the
    decision logic existed and was unreachable: no run ever produced the
    ``tblastn_table_path`` it wants, so every isolate came back ``not_assessed``
    with a reason naming a file that had never been written.

    **The query is built by
    :func:`~papipeline.adapters.oprd_locus.write_tblastn_query` and by nothing
    else.** It carries the reference's terminal stop, so its residue count is
    ``reference_codon_count(reference_cds)``. Handing this screen
    :func:`~papipeline.adapters.oprd_locus.reference_protein`'s 443-residue
    string instead is what made ``stop_is_reference_terminator`` unsatisfiable
    in production - ``qend <= qlen`` can never exceed the query's own length, so
    an alignment of a 443-residue query cannot reach reference codon 444 - and
    the visible symptom was every isolate reading ``disrupted``.

    **Every input is declared and every absence refuses by name**, in the same
    shape as the blastp resolver above: the pinned reference
    (:func:`reference_cds_nucleotides` names the GFF and the FASTA), ``tblastn``
    and ``makeblastdb`` (:func:`~papipeline.adapters.tblastn.run_tblastn` names
    whichever is absent), the assembly and the isolate's Bakta output
    (:func:`_oprd_locus_inputs` names all of them at once), and the search
    table this function writes and then reads back. Nothing here invents a
    verdict to keep going.

    **The verdicts are kept apart from the blastp resolutions on purpose.**
    ``viz.oprd_status_per_sample`` treats a structural verdict as authoritative
    and falls back to ``is_resolved`` otherwise, so merging a structural
    ``not_assessed`` into the carriage mapping would take every sample from
    ``intact`` to ``not_assessed`` - destroying the presence screen - in
    exchange for no new information. A refusal to measure structure is not
    evidence about presence. The two therefore reach the pipeline side by side:
    this mapping is what the report's oprD section reads, and
    ``oprd_locus_resolution`` keeps driving the carriage figure unchanged.

    Args:
        runner: Injected command runner for ``makeblastdb`` and ``tblastn``
            (``runner(command=...) -> CommandResult``). ``None`` runs the real
            tools. This is the seam R12 authorises: the production code path,
            over synthetic inputs, with the subprocesses faked.

    Returns:
        ``None`` outside REAL, for the same reason
        :func:`resolve_oprd_loci` returns ``None`` there. Inside REAL, one
        :class:`~papipeline.adapters.oprd_locus.StructuralCall` per manifest
        sample - the coverage invariant :func:`resolve_oprd_loci` documents
        applies here too.
    """
    if mode is not RunMode.REAL:
        return None

    reference_cds = reference_cds_nucleotides(
        config.reference_gff(), config.reference_fasta(), OPRD_LOCUS_TAG
    )
    work_root = Path(intermediate_root) / OPRD_LOCUS_WORK_DIRNAME
    reference_path = work_root / OPRD_REFERENCE_CDS_NAME
    reference_path.parent.mkdir(parents=True, exist_ok=True)
    reference_path.write_text(
        f">{OPRD_LOCUS_TAG}\n{reference_cds}\n", encoding="utf-8"
    )

    search_dirs = [
        Path(d) for d in config.machine.resolved_tool_search_dirs()
    ] if config.machine is not None else []
    threads = int(config.runtime.get("threads", 1) or 1)

    calls: Dict[str, StructuralCall] = {}
    for sample_id, assembly, _gff3, _faa in _oprd_locus_inputs(
        manifest, genomes_dir, intermediate_root
    ):
        isolate_root = work_root / sample_id / OPRD_STRUCTURAL_DIRNAME
        query = write_tblastn_query(
            isolate_root / "reference_protein.faa", reference_cds
        )
        table = isolate_root / "tblastn.tsv"
        info = tblastn_adapter.run_tblastn(
            query,
            assembly,
            table,
            tool_search_dirs=search_dirs,
            threads=threads,
            runner=runner,
        )
        verdict = structural_call_from_paths(
            sample_id,
            reference_cds_path=reference_path,
            tblastn_table_path=table,
            assembly_path=assembly,
            search=info.search,
        )
        LOGGER.info(
            "oprD structure %s: %s (%d HSP row(s))",
            sample_id, verdict.verdict, info.n_hits,
        )
        calls[sample_id] = verdict
    return calls


def write_oprd_structural_table(
    calls: Mapping[str, StructuralCall], path: Path
) -> Path:
    """Write the per-isolate structural verdicts, and return the path."""
    rows = []
    for sample_id in sorted(calls):
        row = calls[sample_id].as_row()
        rows.append({name: row.get(name) for name in OPRD_STRUCTURAL_COLUMNS})
    return write_tsv(Path(path), rows, list(OPRD_STRUCTURAL_COLUMNS))


def _prepare_figure_data(
    config: PipelineConfig,
    manifest: SampleManifest,
    phenotype: Mapping[str, Any],
    amr: Mapping[str, Sequence[Any]],
    mechanisms: Mapping[str, Sequence[Any]],
    variants: Mapping[str, Sequence[Any]],
    annotated_gene_names: Mapping[str, Sequence[str]],
    lineages: Mapping[str, str],
    mlst: Mapping[str, Optional[Any]],
    gwas_results: Sequence[Any],
    convergence: Sequence[Any],
    cooccurrence: Sequence[Any],
    oprd_locus_resolution: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    """Build every figure's data table. All thirteen, fully data-driven."""
    sample_ids = manifest.sample_ids

    oprd_status = oprd_status_per_sample(
        variants, annotated_gene_names or None, oprd_locus_resolution
    )
    efflux_regulators = efflux_regulator_genes(config)

    phenotypes, phenotype_counts = phenotype_distribution_table(phenotype)
    determinants, determinant_counts = determinant_distribution_table(amr)
    mechanism_names, mechanism_counts = mechanism_distribution_table(mechanisms)
    oprd_phenotypes, oprd_states, oprd_counts = oprd_by_phenotype_table(
        oprd_status, phenotype
    )
    efflux_phenotypes, efflux_genes, efflux_counts = (
        efflux_regulator_by_phenotype_table(
            variants, efflux_regulators, phenotype
        )
    )
    convergence_categories, convergence_counts = convergence_table(convergence)
    nodes, edges = cooccurrence_network_table(cooccurrence)

    return {
        "1_phenotype_distribution": {
            "kind": "bar",
            "title": f"{config.organism.get('name')} phenotype distribution",
            "categories": phenotypes,
            "counts": dict(phenotype_counts),
        },
        "2_determinant_distribution": {
            "kind": "bar",
            "title": "AMR determinant distribution",
            "categories": determinants,
            "counts": dict(determinant_counts),
        },
        "3_mechanism_distribution": {
            "kind": "bar",
            "title": "Resistance mechanism distribution",
            "categories": mechanism_names,
            "counts": dict(mechanism_counts),
        },
        "4_gene_by_phenotype": {
            "kind": "grouped_bar",
            "title": "AMR gene x phenotype",
            "data": gene_by_phenotype_table(amr, phenotype),
        },
        "5_mechanism_by_phenotype": {
            "kind": "grouped_bar",
            "title": "Mechanism x phenotype",
            "data": mechanism_by_phenotype_table(mechanisms, phenotype),
        },
        "6_oprd_by_phenotype": {
            "kind": "stacked_bar",
            "title": "OprD status x phenotype",
            "phenotypes": oprd_phenotypes,
            "states": oprd_states,
            "counts": {f"{p}|{s}": v for (p, s), v in oprd_counts.items()},
        },
        "7_efflux_regulator_by_phenotype": {
            "kind": "grouped_bar",
            "title": "Efflux-regulator alterations x phenotype",
            "phenotypes": efflux_phenotypes,
            "genes": efflux_genes,
            "counts": {f"{p}|{g}": v for (p, g), v in efflux_counts.items()},
        },
        "8_amr_heatmap": {
            "kind": "heatmap",
            "title": "AMR gene presence/absence",
            "matrix": amr_heatmap_matrix(amr, sample_ids),
        },
        "9_mechanism_heatmap": {
            "kind": "heatmap",
            "title": "Mechanism presence/absence",
            "matrix": mechanism_heatmap_matrix(mechanisms, sample_ids),
        },
        "10_tree_annotations": {
            "kind": "table",
            "title": "Phylogeny tip annotations",
            "rows": tree_annotation_table(
                sample_ids, amr, mechanisms, phenotype, lineages, mlst
            ),
        },
        "11_gwas_associations": {
            "kind": "manhattan",
            "title": "GWAS associations",
            "rows": gwas_association_table(gwas_results, config),
        },
        "12_convergence": {
            "kind": "bar",
            "title": "Convergence categories",
            "categories": convergence_categories,
            "counts": dict(convergence_counts),
        },
        "13_cooccurrence_network": {
            "kind": "network",
            "title": "AMR co-occurrence network (association only)",
            "nodes": nodes,
            "edges": edges,
        },
    }


def _figure_manifest(figure_data: Mapping[str, Any]) -> Dict[str, Tuple[str, int]]:
    """(kind, n_items) per figure, for the run summary."""
    out: Dict[str, Tuple[str, int]] = {}
    for name, payload in figure_data.items():
        kind = str(payload.get("kind", "unknown"))
        if kind == "bar":
            n = len(payload.get("categories", []))
        elif kind == "heatmap":
            n = len(payload["matrix"][0])
        elif kind in ("grouped_bar", "stacked_bar"):
            n = len(payload.get("data", {})) or len(payload.get("genes", []))
        elif kind in ("table", "manhattan"):
            n = len(payload.get("rows", []))
        elif kind == "network":
            n = len(payload.get("edges", []))
        else:
            n = 0
        out[name] = (kind, n)
    return out


def _write_figure_data(figure_data: Mapping[str, Any], out_dir: Path) -> None:
    """Persist each figure's prepared table as TSV.

    The plotting layer reads these files, so figure data is inspectable
    without running any plotting code.
    """
    import json

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, payload in figure_data.items():
        path = out_dir / f"{name}.json"
        path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")


def _assemblies_analysed(
    manifest: SampleManifest,
    qc_summary: Optional[Mapping[str, Any]],
) -> int:
    """How many assemblies the run actually read.

    Not ``len(manifest)``. A smoke run's manifest is deliberately the *whole*
    PDC roster - every isolate is a member, so the ten with no prepared assembly
    refuse per sample instead of silently vanishing - which made the report stamp
    "N=967 assemblies" for a run that analysed ten. The marker's own clause,
    "not the full cohort analysis", was thereby false of the number beside it.

    Stage 1 already counts this properly: ``qc_summary['n_samples']`` is the
    number of AssemblyQC records, one per assembly that was read. Preferring it
    also makes the marker agree with the report body, which prints the same
    count when it warns about flagged samples.

    Falls back to the manifest when stage 1 did not run. Reporting a number that
    may be too high beats refusing to build a report at all, and for a full
    cohort the two agree anyway.
    """
    if qc_summary and qc_summary.get("n_samples") is not None:
        return int(qc_summary["n_samples"])
    return len(manifest.samples)


def _build_report_context(
    config: PipelineConfig,
    mode: RunMode,
    antibiotic: str,
    manifest: SampleManifest,
    status: Mapping[str, str],
    warnings: Sequence[str],
    provenance: Sequence[Mapping[str, Any]],
    figure_data: Mapping[str, Any],
    qc_summary: Optional[Mapping[str, Any]] = None,
    alignment_info: Optional[Mapping[str, Any]] = None,
    gwas_input: Optional[Any] = None,
    convergence_calls: Sequence[Any] = (),
    cooccurrence: Sequence[Any] = (),
    phenotype_calls: Sequence[Any] = (),
    stage_refusals: Optional[Mapping[str, str]] = None,
    oprd_calls: Optional[Sequence[Any]] = None,
) -> stage_reporting.ReportContext:
    """Assemble the report context from run artefacts.

    Three of the arguments exist because the report **reads the run's own output
    tables back from disk** (`stages.report_tables.load_run_tables`, merged in
    Phase 2), and the report can only do that if the context carries the tables.
    ``run_tables`` is left to default ``None`` here and resolved by
    ``write_report`` itself from ``config.intermediate_root(mode)`` - the same
    directory either way, so the orchestrator passing it would add a parameter
    without adding a guarantee. The two the orchestrator alone knows are passed:

    * ``stage_refusals`` - stage name -> why it was refused or skipped. Only
      ``run_pipeline`` knows this; the report otherwise derives per-stage
      status from whether a table exists on disk, which cannot distinguish "ran
      and found nothing" from "was refused";
    * ``oprd_calls`` - the per-isolate oprD locus verdicts. **There is no
      contracted path for them** (``contracts.py`` declares no oprD table and the
      Snakefile declares no oprD output), so this is the only route to the
      section and omitting it is what makes the report print "not produced: no
      contracted path exists" about a verdict the run measured. The structural
      verdicts are preferred over the blastp resolutions because they carry
      reading-frame evidence; either shape is accepted by the section.
    """
    sections: List[Tuple[str, str]] = []
    tables: List[Tuple[str, List[str], List[Dict[str, Any]]]] = []

    sections.append(
        (
            "Stages executed",
            stage_reporting.markdown_table(
                ["stage", "status"],
                [
                    {"stage": stage, "status": status.get(stage, "not_run")}
                    for stage in STAGE_ORDER
                ],
            ),
        )
    )

    if qc_summary:
        tables.append(
            (
                "Stage 1 assembly QC summary",
                ["metric", "value"],
                [
                    {"metric": k, "value": v}
                    for k, v in qc_summary.items()
                    if k != "flagged_sample_ids"
                ],
            )
        )

    if alignment_info:
        tables.append(
            (
                "Stage 10 alignment summary",
                ["metric", "value"],
                [{"metric": k, "value": v} for k, v in alignment_info.items()],
            )
        )

    if gwas_input is not None:
        tables.append(
            (
                "Stage 12 GWAS input composition",
                ["metric", "value"],
                [
                    {"metric": "samples tested", "value": gwas_input.n_samples},
                    {"metric": "positive group (R)", "value": gwas_input.n_positive},
                    {"metric": "negative group (S)", "value": gwas_input.n_negative},
                    {"metric": "excluded", "value": gwas_input.n_excluded},
                    {"metric": "features constructed", "value": len(gwas_input.features)},
                    {
                        "metric": "features testable",
                        "value": len(gwas_input.testable_features()),
                    },
                ],
            )
        )

    gwas_rows = figure_data.get("11_gwas_associations", {}).get("rows", [])
    if gwas_rows:
        tables.append(
            (
                "Stage 12 associations (strongest adjusted p-value first)",
                [
                    "feature",
                    "feature_type",
                    "adjusted_p_value",
                    "effect",
                    "frequency",
                    "lineage_linked",
                    "passes_threshold",
                ],
                gwas_rows,
            )
        )
        lineage_linked = [r for r in gwas_rows if r.get("lineage_linked")]
        if lineage_linked:
            sections.append(
                (
                    "Lineage-linked features",
                    "The features below are carried almost entirely by samples from "
                    "a single lineage. An association test cannot separate them from "
                    "that lineage, so they are reported as lineage-linked rather than "
                    "as resistance associations (scientific rule 4).\n\n"
                    + stage_reporting.markdown_table(
                        ["feature", "dominant_lineage_share"],
                        [
                            {
                                "feature": r["feature"],
                                "dominant_lineage_share": r["dominant_lineage_share"],
                            }
                            for r in lineage_linked
                        ],
                    ),
                )
            )

    if convergence_calls:
        tables.append(
            (
                "Stage 13 convergence classification",
                ["determinant", "independent_lineages", "convergence_category", "distribution"],
                [c.to_row() for c in convergence_calls],
            )
        )

    threshold = float(
        (config.raw.get("cooccurrence") or {}).get("significance_threshold", 0.05)
    )
    significant_pairs = [
        c
        for c in cooccurrence
        if c.adjusted_p_value is not None and c.adjusted_p_value <= threshold
    ]
    if significant_pairs:
        tables.append(
            (
                "Stage 14 co-occurrence pairs (association only, not causal)",
                ["feature_a", "feature_b", "feature_type", "n_both", "statistic_value", "adjusted_p_value"],
                [c.to_row() for c in significant_pairs],
            )
        )

    tables.append(
        (
            "Figures prepared",
            ["figure", "kind", "n_items"],
            [
                {"figure": name, "kind": kind, "n_items": n}
                for name, (kind, n) in _figure_manifest(figure_data).items()
            ],
        )
    )

    return stage_reporting.ReportContext(
        mode=mode,
        config=config,
        generated_at=utc_now(),
        antibiotic=antibiotic,
        n_samples=_assemblies_analysed(manifest, qc_summary),
        sections=sections,
        tables=tables,
        warnings=list(warnings),
        provenance=list(provenance),
        # Reporting reads nothing from disk, so it cannot detect its own absent
        # inputs - it is told what it was handed. These are the three
        # finding-bearing aggregates rendered above, and each is dropped
        # silently when empty, so each is what the stage's REAL-mode guard
        # checks (`stages.reporting.REPORTING_REQUIRED_INPUTS`).
        declared_inputs={
            "convergence_calls": list(convergence_calls),
            "cooccurrence": list(cooccurrence),
            "gwas_associations": list(gwas_rows),
        },
        # Stage 11's calls, so `reporting.phenotype_report` - the production
        # caller of `phenotype.provenance_report` - can reach them. Without
        # this the report omitted the provenance section entirely in TEST and
        # STUB, and in REAL printed the "no phenotype calls were handed to this
        # report" absence, which is a claim about this run being false: stage 11
        # had loaded the calls and the orchestrator was holding them.
        #
        # `None` and not `[]` when there are none, deliberately.
        # `phenotype_report` returns `None` only for `None`, and its own
        # docstring rules out the alternative: a zero-row provenance report
        # "would read as" a cohort whose standard was wholly unrecorded, which
        # is the claim `NO_STANDARD_RECORDED` is careful to make only when it is
        # true. `None` reaches the REAL branch, which says the absence is an
        # absence in the run rather than a finding - the honest reading when
        # stage 11 produced nothing at all.
        phenotype_calls=list(phenotype_calls) if phenotype_calls else None,
        # The run's own output tables, read back from disk. Resolved HERE rather
        # than left to `write_report`, so a caller that inspects the context
        # before writing gets the tables rather than `None` - and so the read
        # happens at a point where `status` is already complete, which is what
        # the per-stage table compares the two against.
        run_tables=load_run_tables(config, mode),
        # The per-stage status the orchestrator recorded. Preferred over the
        # report's own disk-derived fallback because only the orchestrator knows
        # why a stage was refused; `stage_refusals` supplies the reason.
        stage_status=dict(status),
        stage_refusals=dict(stage_refusals or {}),
        oprd_calls=list(oprd_calls) if oprd_calls else None,
    )


def _read_previous_manifest(path: Path, mode: RunMode) -> Dict[str, Any]:
    """Whatever the last writer of ``path`` recorded, if it is this mode's.

    Returns ``{}`` for a missing, unreadable, malformed or **different-mode**
    manifest. The mode check is what stops a TEST run's record being carried
    into a REAL one: the two write to different trees in practice, but a
    redirected ``PIPELINE_RESULTS_ROOT`` can put them side by side, and a stage
    that "completed" under synthetic fixtures has completed nothing about a real
    cohort.

    A malformed file is treated as absent rather than raising. The manifest is a
    record, not an input: a run that cannot read the previous one still has
    every fact it needs to write its own, and refusing to run because a JSON file
    is corrupt would be a worse failure than overwriting it.
    """
    if not path.is_file():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    if not isinstance(payload, dict):
        return {}
    if payload.get("run_mode") != mode.value:
        return {}
    return payload


def _write_run_manifest(
    results_root: Path,
    mode: RunMode,
    antibiotic: str,
    status: Mapping[str, str],
    outputs: Mapping[str, Path],
    tools: Mapping[str, Any],
    provenance: Sequence[Mapping[str, Any]],
    skipped: Optional[Mapping[str, str]] = None,
) -> Path:
    """Write ``run_manifest.json``: the reproducibility record for the run.

    ``stages_skipped`` is a first-class key, not an omission. It carries each
    stage that was scheduled and did not run, mapped to the reason - which is
    what turned the silent `recombination` skip from invisible into visible.
    Without it a reader can only infer an omission from an absence, and an
    absence is indistinguishable from a pipeline that has no such stage.

    **A stage is recorded ``completed`` only because some invocation of this run
    completed it, and every stage this invocation completed is listed in
    ``stages_completed_this_run``.** Both halves are load-bearing and they fix
    D10 from two directions:

    * This function used to overwrite ``stages`` with the stages *this process*
      ran. Every Snakefile rule is a separate process (``workflow/Snakefile``'s
      documented KNOWN LIMITATION), so after five completed stages the manifest
      read ``{"virulence": "completed"}`` - stages 1-4 erased, and a reader of
      the reproducibility record had no way to know they had run. The union with
      the previous manifest keeps them.
    * The union alone is not enough: a cumulative map alone cannot say *which*
      stage the writer of this manifest actually ran, so a reader cannot tell a
      stage that completed in this run from one left over from a previous
      invocation of the same mode. ``stages_completed_this_run`` says it, and it
      is **exactly** the stages this process marked - no inference, no
      timestamps, and nothing a stale file could widen.

    The union is over stages **some invocation completed**, so it can never
    contain a stage that did not run. That is the invariant: the manifest does
    not invent a completion, it stops forgetting one.
    """
    path = Path(results_root) / "run_manifest.json"
    previous = _read_previous_manifest(path, mode)

    stages: Dict[str, str] = {
        str(name): str(state)
        for name, state in (previous.get("stages") or {}).items()
    }
    stages.update({str(name): str(state) for name, state in status.items()})

    skipped_map: Dict[str, str] = {
        str(name): str(reason)
        for name, reason in (previous.get("stages_skipped") or {}).items()
    }
    skipped_map.update({str(name): str(r) for name, r in (skipped or {}).items()})

    completed_now = sorted(
        str(name) for name, state in status.items() if str(state) == "completed"
    )

    outputs_map: Dict[str, str] = {
        str(name): str(value)
        for name, value in (previous.get("outputs") or {}).items()
    }
    outputs_map.update({str(k): str(v) for k, v in outputs.items()})

    payload = {
        "pipeline_version": __version__,
        "run_mode": mode.value,
        "generated_at_utc": utc_now(),
        "antibiotic": antibiotic,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "stages": stages,
        # Not a filter over `stages` and not a timestamp: the exact set this
        # invocation completed. Present so a reader can always answer "did THIS
        # run do it" without inferring from a cumulative map.
        "stages_completed_this_run": completed_now,
        "stages_skipped": skipped_map,
        "outputs": outputs_map,
        "tools_detected": {
            name: {
                "executable": tool.executable,
                "version": tool.version,
                "available": tool.available,
            }
            for name, tool in tools.items()
        },
        "references": list(provenance),
    }
    path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    LOGGER.debug("Run manifest written to %s", path)
    return path
