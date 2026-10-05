"""Typed accessors over the YAML configuration and the knowledge tables.

Design rules:

* configuration is read once into frozen-ish dataclasses, never re-parsed
  inside a stage;
* the knowledge tables (``mechanisms``/``regulators``/``references``) are
  editable data files, so lookups go through
  :mod:`papipeline.knowledge` and unknown entries raise a typed error
  instead of being silently ignored;
* no sample identifiers and no biological results live in configuration.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import os

import yaml

from ..errors import AntibioticError, ConfigError, SampleCapExceeded
from ..io.tsv import read_tsv
from ..logging_utils import get_logger
from ..models import ClaimStatus, Phenotype, RunMode

LOGGER = get_logger("config")

DEFAULT_CONFIG_RELPATH = Path("config/science.yaml")

#: The repository root, derived from this module's own location rather than from
#: the working directory, so a configured RELATIVE path means the same file no
#: matter where the command was invoked from. Used by
#: :meth:`CohortConfig.resolve_subset_file` and by :func:`_resolve_overlay`.
REPO_ROOT = Path(__file__).resolve().parents[2]

#: Environment variable that redirects a run's outputs. See
#: :meth:`PipelineConfig.results_root`.
RESULTS_ROOT_ENV = "PIPELINE_RESULTS_ROOT"

#: Environment variable that opens REAL mode for one session, without editing a
#: tracked overlay. See :func:`_allow_real_mode_override`.
#:
#: An environment variable rather than a local overlay file because the
#: committed overlay must keep REAL mode shut - two tests assert exactly that,
#: and they should keep guarding it. A run that had to edit the file to enable
#: itself made the working tree dirty for its whole duration, so nothing else
#: could be committed while it ran, and those two tests failed for as long as it
#: existed. The variable is inherited by every Snakemake worker, which a
#: gitignored overlay merged at load would not be.
ALLOW_REAL_MODE_ENV = "PIPELINE_ALLOW_REAL_MODE"

#: Spellings accepted as "yes". Anything else is an error rather than a guess:
#: a typo that happened to be truthy would enable a mode that shuts by design.
_TRUTHY = frozenset({"1", "true", "yes", "on"})
_FALSEY = frozenset({"", "0", "false", "no", "off"})
MACHINES_DIRNAME = "machines"
#: Names accepted by ``load_config(machine=...)``.
#:
#: ``smoke`` is here because ``config/machines/smoke.yaml`` is a committed,
#: validated overlay that exists to authorise a bounded REAL run, and the name
#: did not resolve (D12). Every invocation had to pass the *path*
#: ``--machine config/machines/smoke.yaml`` instead, which the Snakefile
#: interpolates verbatim into a shell command at twelve places - so the overlay
#: was resolved against the CWD rather than the worktree root, and a run from
#: anywhere but the repository root named a file that was not there. A committed
#: overlay that cannot be named is an overlay every caller has to spell by hand.
#:
#: Naming it here grants nothing: ``runtime.allow_real_mode`` is still ``false``
#: in that file and ``resolve_mode`` still refuses REAL. It makes the overlay
#: *addressable*, which is what a name is for.
KNOWN_MACHINES = ("laptop", "bigmachine", "smoke")

#: Top-level keys a machine overlay owns. Anything else must live in
#: science.yaml; see ``load_machine_config``.
#:
#: ``reference`` was once here, which meant each machine could name its own
#: reference genome and two machines could analyse against different assemblies
#: without complaint. Which reference the science uses is a scientific fact; only
#: where the file lives (``paths.reference_fasta``) is per-machine.
MACHINE_OWNED_KEYS = frozenset(
    {"machine", "machine_class", "paths", "runtime", "tools"}
)

#: The only keys a per-tool entry in an overlay may carry: whether the tool
#: exists on this machine, and where it is. A version or a database release here
#: would be a second source of truth for pinned versions, which
#: ``config/references.tsv`` already is.
TOOL_ENTRY_KEYS = frozenset({"available", "path"})

#: Multiple-testing corrections step 12a can compute.
#:
#: Exactly one, because pyseer 1.1.2 documents exactly one: its own
#: ``scripts/count_patterns.py``, vendored at ``scripts/gwas/count_patterns.py``
#: because the wheel does not package it. That script's whole logic is
#: ``sort -u | wc -l`` followed by ``alpha / count``. There is no flag for a
#: second correction, so accepting a second name in configuration would only
#: ever produce a threshold that some other code computed - a number the run
#: reported without having applied it.
SUPPORTED_CORRECTIONS = frozenset({"bonferroni"})


# --------------------------------------------------------------------------
# Knowledge table records
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class AntibioticSpec:
    """One row of ``config/antibiotics.tsv``."""

    antibiotic: str
    antibiotic_class: str
    mechanism_scope: str
    phenotype_source_standard: str
    notes: Optional[str] = None

    def applies_to(self, configured: Sequence[str]) -> bool:
        return self.antibiotic in configured


@dataclass(frozen=True)
class MechanismSpec:
    """One row of ``config/mechanisms.tsv``.

    ``claim_ceiling`` is the strongest :class:`ClaimStatus` that a raw
    detection of this gene alone may be reported as. It is enforced by stage 5
    so that a gene detection can never be escalated to an association.
    """

    mechanism: str
    gene: str
    gene_type: str
    antibiotic: str
    biological_role: str
    evidence_level: ClaimStatus
    mechanism_class: str
    gene_role: str
    claim_ceiling: ClaimStatus
    notes: Optional[str] = None

    def relevant_to(self, antibiotic: str) -> bool:
        return self.antibiotic in ("all", antibiotic, "all_beta_lactams") or (
            antibiotic in self.antibiotic.split(",")
        )


@dataclass(frozen=True)
class RegulatorSpec:
    """One row of ``config/regulators.tsv``."""

    gene: str
    mechanism: str
    gene_role: str
    screen_for: Tuple[str, ...]
    promoter_screen: bool
    reference_length: Optional[int]
    evidence_level: ClaimStatus
    notes: Optional[str] = None
    #: PAO1 RefSeq locus tag from `config/regulators.tsv`, resolved and verified
    #: by scripts/resolve_locus_tags.py. Optional, because a row may exist with
    #: no counterpart in PAO1 - `mexS` has none. A locus WITHOUT a tag cannot be
    #: screened against reference coordinates, so `screenable` reports that
    #: rather than the caller discovering it mid-run.
    locus_tag: Optional[str] = None

    @property
    def screenable(self) -> bool:
        """True when this locus can be resolved against the PAO1 reference."""
        return bool(self.locus_tag)

    @property
    def variant_classes(self) -> frozenset:
        return frozenset(self.screen_for)


def _allow_real_mode_override(configured: bool) -> bool:
    """Resolve `runtime.allow_real_mode`, honouring the session override.

    The environment can only *open* a gate the committed overlay left shut. It
    cannot close one that is open, because a stray `PIPELINE_ALLOW_REAL_MODE=0`
    in someone's shell would then silently re-shut a machine whose overlay
    deliberately enables REAL mode - and the resulting failure would look like
    the overlay's own doing.

    Raises:
        ConfigError: The variable is set to something unrecognised. Guessing
            would mean a mistyped value silently decided whether real genomes
            were analysed.
    """
    raw = os.environ.get(ALLOW_REAL_MODE_ENV)
    if raw is None:
        return configured
    value = raw.strip().lower()
    if value in _TRUTHY:
        return True
    if value in _FALSEY:
        return configured
    raise ConfigError(
        f"{ALLOW_REAL_MODE_ENV}={raw!r} is not a value this pipeline understands. "
        f"Use one of {sorted(_TRUTHY)} to enable REAL mode, or one of "
        f"{sorted(_FALSEY)} to leave the overlay's value alone. It is refused "
        "rather than guessed, because this variable decides whether real genomes "
        "are analysed.",
        env_var=ALLOW_REAL_MODE_ENV, value=raw,
    )


@dataclass(frozen=True)
class ReferenceSpec:
    """One row of ``config/references.tsv``."""

    reference_id: str
    tool: Optional[str]
    tool_version: Optional[str]
    database: Optional[str]
    database_version: Optional[str]
    version_status: str
    source: Optional[str] = None
    notes: Optional[str] = None

    @property
    def is_pinned(self) -> bool:
        return self.version_status == "pinned"


# --------------------------------------------------------------------------
# Sub-configuration views
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class QcConfig:
    min_assembly_size: Optional[int]
    max_assembly_size: Optional[int]
    min_n50: Optional[int]
    max_ambiguous_bases: Optional[int]
    hooks: Mapping[str, bool] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, raw: Mapping[str, Any]) -> "QcConfig":
        return cls(
            min_assembly_size=raw.get("min_assembly_size"),
            max_assembly_size=raw.get("max_assembly_size"),
            min_n50=raw.get("min_n50"),
            max_ambiguous_bases=raw.get("max_ambiguous_bases"),
            hooks=dict(raw.get("hooks") or {}),
        )


@dataclass(frozen=True)
class UniquePatternConfig:
    """Step 12a: the variant-matrix reduction that sets the testing threshold.

    ``alpha`` is the family-wise error rate; the threshold is
    ``alpha / n_unique_patterns``. ``method`` names the correction and is
    accepted only as ``bonferroni``, because that is the one pyseer documents
    and the one the vendored helper implements. A second name here would be a
    claim the pipeline cannot honour, so an unknown one is a configuration
    error rather than something to be ignored.

    Both keys are **required**, with no default. That is the same treatment the
    machine overlay's helper limits get, for the same reason: each has an
    upstream or conventional default (the helper's ``0.05``, the sole supported
    correction), and inheriting one would leave an unconfigured run reporting a
    threshold as though someone had chosen it. ``enabled`` was carried here for
    a while and read by nothing - a flag that changes nothing is worse than no
    flag - so it is gone rather than left as decoration.
    """

    method: str
    alpha: float

    @classmethod
    def from_dict(cls, raw: Mapping[str, Any]) -> "UniquePatternConfig":
        method = str(raw.get("method", "")).strip().lower()
        if not method:
            raise ConfigError(
                "gwas.unique_patterns.method is required. It is the multiple-"
                "testing correction step 12a applies, and the run reports the "
                "threshold that correction produced; leaving it unset would "
                "mean reporting a number without saying how it was derived."
            )
        if method not in SUPPORTED_CORRECTIONS:
            raise ConfigError(
                f"gwas.unique_patterns.method is {method!r}, which this pipeline "
                f"does not implement. Supported: "
                f"{', '.join(sorted(SUPPORTED_CORRECTIONS))}. Declaring a "
                f"correction the vendored helper cannot compute would leave the "
                f"multiple-testing threshold undefined while the run reported a "
                f"number."
            )
        if raw.get("alpha") is None:
            raise ConfigError(
                "gwas.unique_patterns.alpha is required. It is the family-wise "
                "error rate the Bonferroni correction is applied at, and it is "
                "a scientific constant: inheriting the helper's own default of "
                "0.05 would make an unconfigured run report a threshold as "
                "though someone had chosen it."
            )
        alpha = float(raw["alpha"])
        if not 0.0 < alpha <= 1.0:
            raise ConfigError(
                f"gwas.unique_patterns.alpha is {alpha}; a family-wise error "
                f"rate lies in (0, 1].",
                alpha=alpha,
            )
        return cls(method=method, alpha=alpha)


@dataclass(frozen=True)
class GwasConfig:
    enabled: bool
    tool: str
    feature_types: Tuple[str, ...]
    positive: Tuple[str, ...]
    negative: Tuple[str, ...]
    excluded: Tuple[str, ...]
    min_samples_per_group: int
    correction: str
    significance_threshold: float
    filter_model: str
    kinship: str
    lineage_confound_threshold: float
    unique_patterns: UniquePatternConfig

    @classmethod
    def from_dict(cls, raw: Mapping[str, Any]) -> "GwasConfig":
        outcome = raw.get("outcome") or {}
        return cls(
            enabled=bool(raw.get("enabled", False)),
            tool=str(raw.get("tool", "pyseer")),
            feature_types=tuple(raw.get("feature_types") or ()),
            positive=tuple(outcome.get("positive") or ()),
            negative=tuple(outcome.get("negative") or ()),
            excluded=tuple(outcome.get("excluded") or ()),
            min_samples_per_group=int(raw.get("min_samples_per_group", 3)),
            correction=str(raw.get("correction", "fdr_bh")),
            significance_threshold=float(raw.get("significance_threshold", 0.05)),
            filter_model=str(raw.get("filter", "mixed")),
            kinship=str(raw.get("kinship", "unitig")),
            lineage_confound_threshold=float(
                raw.get("lineage_confound_threshold", 0.9)
            ),
            unique_patterns=UniquePatternConfig.from_dict(
                raw.get("unique_patterns") or {}
            ),
        )


@dataclass(frozen=True)
class PhylogenyConfig:
    aligner: str
    snp_caller: str
    tree_builder: str
    models: str
    require_exact_sample_match: bool
    min_core_alignment_fraction: float
    asc_drop_partially_constant: bool = False

    @classmethod
    def from_dict(cls, raw: Mapping[str, Any]) -> "PhylogenyConfig":
        return cls(
            aligner=str(raw.get("aligner", "panaroo")),
            snp_caller=str(raw.get("snp_caller", "snp_sites")),
            tree_builder=str(raw.get("tree_builder", "iqtree2")),
            models=str(raw.get("models", "GTR+G")),
            require_exact_sample_match=bool(
                raw.get("require_exact_sample_match", True)
            ),
            min_core_alignment_fraction=float(
                raw.get("min_core_alignment_fraction", 0.5)
            ),
            # DEFAULT False, and the default is the conservative direction:
            # this key says whether columns constant after gaps/Ns are ignored
            # are dropped before a +ASC fit. A caller that has not chosen gets
            # the alignment exactly as committed, so the alignment a test
            # fixture pins is still the alignment that is analysed. It is not
            # the same switch as `models`: `models` reaches the argv as `-m`
            # (see papipeline/adapters/iqtree.py build_command), whereas this
            # key is a preprocessing decision that no argv carries yet.
            asc_drop_partially_constant=bool(
                raw.get("asc_drop_partially_constant", False)
            ),
        )


#: `variants.bcftools.mpileup.min_ireads`'s DEFAULT.
#:
#: Ruling R6: this pipeline's input is assembled genomes, and an assembly is a
#: single consensus sequence, so there is exactly ONE gapped read per indel.
#: bcftools 1.23.1 defaults this to 2 -- a READS threshold, verified from
#: `bcftools mpileup` usage line 53:
#:
#:   -m, --min-ireads INT   Minimum number gapped reads for indel
#:                           candidates [2]
#:
#: At the bcftools default no assembly indel can ever satisfy it, so every one
#: is discarded before being considered and stage 6 reports the locus as clean.
#: 1 is therefore not a loosened threshold chosen for sensitivity; it is the
#: only value consistent with a depth-1 input.
#:
#: It is NOT a read-depth filter. `-d`/`max_depth` is the depth knob and is
#: separate.
DEFAULT_MPILEUP_MIN_IREADS: int = 1

#: `variants.bcftools.mpileup.max_depth`'s DEFAULT.
#:
#: **INHERITED, NOT CHOSEN.** This is bcftools 1.23.1's own default, verified
#: from the usage block printed by a bare `bcftools mpileup`:
#:
#:   -d, --max-depth INT   Max raw per-file depth; avoids excessive memory
#:                         usage [250]
#:
#: It was written into `science.yaml` to record what the tool does, and no
#: measurement in this repository chose it. It is stated here, as one named
#: number, because the adapter used to carry its own bare `250` literal: two
#: places held the same fact and one of them was invisible.
#:
#: **What it does, measured.** `-d` is bcftools' only memory lever, and it is
#: monotone, but it engages ONLY where true depth exceeds it. Against the pinned
#: PAO1 reference (6,264,404 bp), peak RSS:
#:
#: * depth 1 (the ordinary assembly case, `minimap2 -a --secondary=yes`):
#:   94 MB at `-d 250`, and `-d 1` is indistinguishable. The cap never binds.
#: * 60 x 500 kb contig records stacked on one locus (true depth 60, which is
#:   what `--secondary=yes` yields at a repeat locus): 127 MB at `-d 250`,
#:   77 MB at `-d 20`, 52 MB at `-d 1`. So the cap is worth ~2.4x there.
#:
#: **Why it is still 250.** `-d` is a SUBSAMPLING cap: past it, mpileup keeps a
#: subset of reads, which moves `FORMAT/DP` and `FORMAT/AD` and therefore which
#: alleles `bcftools call -mv` emits. `science.yaml` records that `QUAL` is
#: constant at 30.4183 across every called site, so DP/AD cannot be recovered
#: downstream and the call genuinely moves. Lowering this number is a SCIENCE
#: decision about the callable variant set, not an engineering one, and it is
#: deliberately left alone until that decision is made.
DEFAULT_MPILEUP_MAX_DEPTH: int = 250

#: `variants.bcftools.mpileup.output_type`'s DEFAULT.
#:
#: `'z'` is bgzf-COMPRESSED VCF, verified from the same usage block:
#:
#:   -O, --output-type TYPE  'b' compressed BCF; 'u' uncompressed BCF;
#:                           'z' compressed VCF; 'v' uncompressed VCF;
#:                           0-9 compression level [v]
#:
#: Stage 6's dominant UNBOUNDED cost is not memory but output VOLUME. With
#: `-a FORMAT/DP,FORMAT/AD` across a whole 6.26 Mb reference, mpileup emits a
#: record for every covered position: measured 750 MB of raw VCF for ONE isolate,
#: so ~7.5 GB for the 10-isolate smoke cohort and ~625 GB for 835. `'z'` cuts
#: that 40x, measured 750 MB -> 19 MB.
#:
#: It buys no RSS: measured 127.4 MB with `-O z` against 125.2 MB with `-O v`.
#: It is here because the volume is real and unbounded, not because it is the
#: fix for a SIGKILL.
#:
#: **The extension is load-bearing and the failure is SILENT.** `-O z` with an
#: output filename that does not end `.gz` writes UNCOMPRESSED data and exits 0
#: with no warning: measured `-O z -o mismatch.vcf` produced 69,488,493 bytes,
#: byte-for-byte the size of the `-O v` output, while `-O z -o ok.vcf.gz`
#: produced 1,861,798. So `RAW_VCF` in the adapter is `raw.vcf.gz` and not
#: `raw.vcf`; changing the flag alone would have been a no-op dressed as a fix.
#:
#: Compressing does not change what is called: same BAM, same `-q/-Q/-a/-d/-m`,
#: `-O v` vs `-O z`, and the `bcftools call -mv` output is identical.
DEFAULT_MPILEUP_OUTPUT_TYPE: str = "z"

#: `lineage.method`'s DEFAULT, and the only method this project accepts today.
#:
#: Ruling R4: `lineage_label` is the MLST sequence type from stage 3. The method
#: is recorded so a result always says how its lineages were defined.
#:
#: The tree-cut alternative is DEFERRED, deliberately: it needs a support
#: threshold, a minimum clade size, and a rule for isolates the tree places
#: ambiguously, and none of those has been chosen. An accepted-values set is
#: the place that keeps a deferred method from being spelled into a config file
#: and then silently believed by a reader.
DEFAULT_LINEAGE_METHOD: str = "st"

#: Accepted values for `lineage.method`.
#:
#: Deliberately one entry. See :data:`DEFAULT_LINEAGE_METHOD`.
ACCEPTED_LINEAGE_METHODS: Tuple[str, ...] = (DEFAULT_LINEAGE_METHOD,)

#: `annotation.reuse_tool_output`'s DEFAULT, and the only value that changes no
#: behaviour at all.
#:
#: Reuse means "trust Bakta output already on disk instead of invoking Bakta",
#: so a default of anything else would make every future run of this pipeline
#: silently inherit a trust decision nobody made for it - including runs on
#: machines and cohorts that did not exist when the default was chosen.
DEFAULT_REUSE_TOOL_OUTPUT: str = "off"

#: Accepted values for `annotation.reuse_tool_output`.
#:
#:   ``off``     never reuse; every genome is annotated by the tool. The default.
#:   ``prefer``  reuse per-sample output that verifies, run the tool for the rest.
#:   ``require`` reuse per-sample output that verifies; record a refusal, by
#:               name, for the ones that do not. The tool is never invoked, so
#:               this is the mode that makes a Bakta-free run possible at all.
#:
#: Three and not two, because `prefer` cannot answer the question the refusal
#: exists for: when no output is verifiable, `prefer` quietly spends an hour of
#: CPU re-annotating the whole cohort, which on a cohort that was assembled
#: never to be annotated is not a fallback but a second way to get the wrong
#: answer. Both non-default modes are per-sample decisions; neither is a
#: cohort-wide verdict.
REUSE_TOOL_OUTPUT_MODES: Tuple[str, ...] = ("off", "prefer", "require")


@dataclass(frozen=True)
class LineageConfig:
    """How a sample's lineage label is defined."""

    method: str = DEFAULT_LINEAGE_METHOD

    @classmethod
    def from_dict(cls, raw: Mapping[str, Any]) -> "LineageConfig":
        method = str(raw.get("method", DEFAULT_LINEAGE_METHOD)).strip()
        if method not in ACCEPTED_LINEAGE_METHODS:
            # Refused rather than accepted-and-ignored. A lineage_label is the
            # grouping variable for convergence and co-occurrence, so an
            # unrecognised method would not fail loudly -- it would regroup the
            # cohort by an unknown rule and report the result as if it were the
            # rule this project uses.
            raise ConfigError(
                "lineage.method is not a method this project accepts",
                path="lineage.method",
                configured=method,
                accepted=",".join(ACCEPTED_LINEAGE_METHODS),
            )
        return cls(method=method)


#: `cohort.subset_file`'s DEFAULT: no subset, i.e. the whole manifest.
#:
#: `None` rather than `""` because the two mean different things and only one of
#: them is a refusal. `None` is "no subset was configured"; an empty string is
#: "a subset was configured and it is empty", which selects zero samples and is
#: a configuration error rather than a valid run.
DEFAULT_COHORT_SUBSET_FILE: Optional[str] = None


@dataclass(frozen=True)
class CohortConfig:
    """Which manifest members this run analyses."""

    #: Repository-relative path to a newline-delimited list of sample_ids, or
    #: ``None`` for "the whole roster".
    subset_file: Optional[str] = DEFAULT_COHORT_SUBSET_FILE

    @classmethod
    def from_dict(cls, raw: Mapping[str, Any]) -> "CohortConfig":
        value = raw.get("subset_file", DEFAULT_COHORT_SUBSET_FILE)
        if value is None:
            return cls(subset_file=None)
        text = str(value).strip()
        if not text:
            # Distinguished from None on purpose. An empty subset selects no
            # samples at all, and a run that analyses nothing while reporting
            # success is the failure this key must not be able to express.
            raise ConfigError(
                "cohort.subset_file is set but empty. Null means 'no subset, "
                "the whole manifest'; an empty path selects zero samples and "
                "is a configuration error, not a valid empty subset.",
                path="cohort.subset_file",
                configured=repr(value),
            )
        return cls(subset_file=text)

    def resolve_subset_file(self) -> Optional[Path]:
        """The absolute path of the member list, or ``None`` for no subset.

        A relative `subset_file` is resolved against :data:`REPO_ROOT`, NOT
        against the working directory. Two reasons, and both have bitten:

        * A Snakemake rule's shell command runs with whatever CWD the scheduler
          happened to choose, so a CWD-relative path names a different file in
          different jobs -- or no file at all.
        * `local/` (and `data/`, and `db/`) are UNTRACKED SYMLINKS into the
          canonical worktree in a dev worktree. Resolving through them would
          make the list this worktree reads a property of another worktree's
          disk, which is precisely what a per-worktree local file exists to
          avoid.

        Absolute paths are returned unchanged, so an operator who needs to point
        at a list outside the repository can still do so.
        """
        if self.subset_file is None:
            return None
        candidate = Path(self.subset_file)
        if candidate.is_absolute():
            return Path(os.path.normpath(candidate))
        # normpath, not resolve(): it tidies `repo/../x` without following
        # symlinks, and following them is exactly what must not happen here.
        return Path(os.path.normpath(REPO_ROOT / candidate))

    def read_subset_file(self) -> Optional[List[str]]:
        """The isolate ids this run analyses, or ``None`` for the whole roster.

        Reads through :meth:`resolve_subset_file`. Raises :class:`ConfigError` if
        the configured list does not exist, naming the resolved path -- an
        absent member list must never degrade silently into "no subset", because
        that reads as the full cohort while the operator asked for ten isolates.
        """
        path = self.resolve_subset_file()
        if path is None:
            return None
        if not path.is_file():
            raise ConfigError(
                "cohort.subset_file names a file that does not exist. Falling "
                "back to the whole roster would analyse the full cohort while "
                "appearing to be the bounded run that was asked for.",
                path="cohort.subset_file",
                configured=self.subset_file,
                resolved=str(path),
            )
        ids = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
        members = [i for i in ids if i and not i.startswith("#")]
        if not members:
            raise ConfigError(
                "cohort.subset_file is empty. Null means 'no subset, the whole "
                "manifest'; a file with no ids is a configuration error.",
                path="cohort.subset_file",
                resolved=str(path),
            )
        return members


@dataclass(frozen=True)
class ConvergenceConfig:
    min_independent_lineages: int
    widespread_fraction: float
    rare_max_samples: int
    minimum_branch_support: float

    @classmethod
    def from_dict(cls, raw: Mapping[str, Any]) -> "ConvergenceConfig":
        return cls(
            min_independent_lineages=int(raw.get("min_independent_lineages", 2)),
            widespread_fraction=float(raw.get("widespread_fraction", 0.75)),
            rare_max_samples=int(raw.get("rare_max_samples", 2)),
            minimum_branch_support=float(raw.get("minimum_branch_support", 0.0)),
        )


# --------------------------------------------------------------------------
# Top-level configuration
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class PipelineConfig:
    """Fully resolved configuration for one pipeline run."""

    root: Path
    raw: Mapping[str, Any]
    antibiotics: Tuple[str, ...]
    allowed_phenotypes: Tuple[Phenotype, ...]
    analysis: Mapping[str, bool]
    mechanisms: Mapping[str, MechanismSpec]
    regulators: Mapping[str, RegulatorSpec]
    references: Mapping[str, ReferenceSpec]
    antibiotic_specs: Mapping[str, AntibioticSpec]
    organism: Mapping[str, Any]
    qc: QcConfig
    gwas: GwasConfig
    phylogeny: PhylogenyConfig
    convergence: ConvergenceConfig
    paths: Mapping[str, Any]
    runtime: Mapping[str, Any]
    machine: Optional[MachineConfig] = None
    lineage: LineageConfig = LineageConfig()
    cohort: CohortConfig = CohortConfig()

    # -- machine-owned settings, read through the overlay -----------------
    @property
    def threads(self) -> Optional[int]:
        return (self.runtime or {}).get("threads")

    @property
    def memory_mb(self) -> Optional[int]:
        return (self.runtime or {}).get("memory_mb")

    @property
    def max_samples(self) -> Optional[int]:
        """The machine's sample cap, or ``None`` when uncapped.

        ``None`` is a real, meaningful answer here: the analysis machine is
        deliberately uncapped, which is not the same as not having looked.
        """
        return (self.runtime or {}).get("max_samples")

    # -- step 12a: the vendored unique-pattern helper's own limits --------
    #
    # These bound a *separate process*. The vendored helper shells out to
    # `sort -S <memory> -T <dir>`, and `sort` has an address space of its own,
    # so the pipeline's `memory_mb` is not a statement about it. Each accessor
    # raises rather than falling back to a number, because every one of these
    # has an upstream default that would then be silently inherited: a run
    # would appear to be configured and would actually be guessing.

    @property
    def unique_patterns_memory_mb(self) -> int:
        """Memory ceiling handed to the helper's `sort`, in Mb.

        Upstream defaults to 1024. Inherited silently, that default would look
        like a configuration decision nobody made.
        """
        value = (self.runtime or {}).get("unique_patterns_memory_mb")
        if value is None:
            raise ConfigError(
                "runtime.unique_patterns_memory_mb is not set. The vendored "
                "unique-pattern helper shells out to `sort`, which has an "
                "address space of its own, so this pipeline's runtime."
                "memory_mb is not a bound for it. Add the key to the machine "
                "overlay."
            )
        return int(value)

    @property
    def unique_patterns_cores(self) -> int:
        """Cores handed to the helper's `sort`.

        Upstream defaults to 1. Above 1 the helper adds `--parallel=<n>`, so
        this key decides whether the sort is single- or multi-threaded.
        """
        value = (self.runtime or {}).get("unique_patterns_cores")
        if value is None:
            raise ConfigError(
                "runtime.unique_patterns_cores is not set. It decides whether "
                "the helper's `sort` runs with --parallel; add the key to the "
                "machine overlay rather than letting the helper's own default "
                "of 1 stand."
            )
        return int(value)

    @property
    def unique_patterns_temp_dir(self) -> Path:
        """Scratch directory handed to the helper's `sort -T`.

        `sort` spills to this directory on a patterns file too large for
        memory, so it must be a real, writable directory - not a path someone
        assumed exists. Upstream defaults to ``/tmp``.

        Declared under ``runtime`` rather than ``paths`` because ``paths``
        values are repository-relative by an invariant this repository enforces
        (see ``tests/unit/test_reference_ownership.py``), and a scratch
        directory is by nature a location on one machine's filesystem.
        """
        value = (self.runtime or {}).get("unique_patterns_temp_dir")
        if not value:
            raise ConfigError(
                "runtime.unique_patterns_temp_dir is not set. It is where the "
                "helper's `sort -T` spills when the patterns file exceeds "
                "runtime.unique_patterns_memory_mb, so it has to be a real "
                "directory; add the key to the machine overlay rather than "
                "letting the helper's own /tmp default stand."
            )
        return Path(str(value)).expanduser()

    # -- sample identity ---------------------------------------------------
    @property
    def sample_id_pattern(self) -> Optional[str]:
        """The configured sample-id format rule, if one is declared.

        ``None`` means the structural checks apply and no pattern does. That is
        a real distinction: an undeclared pattern is a gap in the configuration,
        not a licence to accept anything.
        """
        section = (self.raw or {}).get("sample_id") or {}
        return section.get("pattern") or None

    @property
    def sample_id_unique_attributes(self) -> Tuple[str, ...]:
        """Provenance attributes that must identify at most one assembly."""
        section = (self.raw or {}).get("sample_id") or {}
        raw = section.get("unique_attributes")
        if not raw:
            return ("assembly", "isolate_id", "biosample")
        if isinstance(raw, str):
            return (raw,)
        return tuple(str(item) for item in raw)

    @property
    def machine_name(self) -> Optional[str]:
        return self.machine.name if self.machine else None

    def is_smoke_overlay(self) -> bool:
        """Whether this overlay declares itself the bounded smoke run.

        True when the overlay sets ``paths.smoke_genome_dir``. That key is the
        opt-in: it is what says "read the prepared assemblies from here", so an
        overlay without it is an ordinary one and is treated as such.

        Kept as a single predicate because the question is asked from two
        places - manifest discovery and report marking - and answering it
        twice invites them to disagree about which reports are smoke reports
        and which runs read the smoke genomes. A run that resolved 10
        assemblies but produced a report with no marker, or a report marked
        SMOKE over a full cohort, would both come from the two disagreeing.
        """
        if self.machine is None:
            return False
        return bool((self.machine.paths or {}).get("smoke_genome_dir"))

    def tool_available(self, tool: str) -> Optional[bool]:
        if self.machine is None:
            return None
        return self.machine.tool_available(tool)

    def tool_candidates(self, tool: str) -> Tuple[str, ...]:
        """Locations to try for ``tool``. Empty without a machine overlay."""
        if self.machine is None:
            return ()
        return self.machine.tool_candidates(tool)

    def reference_gff(self) -> Path:
        """The pinned reference's gene annotation, derived from its own FASTA.

        Stage 6's variants are in PAO1 coordinates, so resolving them to genes
        needs the annotation in that same frame. The path is derived from the
        FASTA's rather than configured separately, so the two cannot drift onto
        different assemblies - which is the failure a GFF from another release
        produces: every gene boundary silently shifted, nothing raised.
        """
        return self.reference_fasta().with_suffix(".gff")

    def reference_fasta(self) -> Path:
        """The reference this run is configured against, via the machine overlay.

        The *identity* is science, read from ``raw['reference']``; only the path
        is per-machine. Without this delegating wrapper the only way to reach
        the reference was to know the overlay carried it and reach through
        ``config.machine`` - so nothing called it, and D5 verification had no
        path to check. Same shape as :meth:`data_root_for`.
        """
        identity = (self.raw or {}).get("reference") or {}
        if self.machine is not None:
            return self.machine.reference_fasta(identity)
        name = identity.get("accession") or "reference"
        return self.root / "db" / "reference" / str(name) / f"{name}_genomic.fna"

    def data_root_for(self, mode: RunMode) -> Path:
        """Data root for ``mode``, via the machine overlay when there is one."""
        if self.machine is not None:
            return self.machine.data_root_for(mode)
        return self.data_root(mode)

    def enforce_sample_cap(self, sample_count: int) -> None:
        """Refuse a cohort this machine cannot carry. See MachineConfig."""
        if self.machine is None:
            return
        self.machine.enforce_sample_cap(sample_count)

    # -- derived paths ----------------------------------------------------
    @property
    def config_dir(self) -> Path:
        return self.root / "config"

    def stage_enabled(self, stage: str) -> bool:
        return bool(self.analysis.get(stage, False))

    def variants_mpileup_min_ireads(self) -> int:
        """`variants.bcftools.mpileup.min_ireads`, with its DEFAULT.

        DEFAULT is :data:`DEFAULT_MPILEUP_MIN_IREADS`, and it is stated HERE
        rather than inline at the call site so that the number the loader
        promises and the number the adapter falls back to are one number.

        Reads through the whole chain defensively: this project's input is
        assembled genomes, and a run that has not configured the mpileup block
        at all is a legitimate STUB/TEST run, not a malformed config. A missing
        section yields the DEFAULT, not a KeyError.
        """
        variants = self.raw.get("variants") or {}
        bcftools = variants.get("bcftools") or {}
        mpileup = bcftools.get("mpileup") or {}
        value = mpileup.get("min_ireads", DEFAULT_MPILEUP_MIN_IREADS)
        return int(value)

    def reuse_tool_output(self) -> str:
        """`annotation.reuse_tool_output`, validated against its accepted set.

        DEFAULT is :data:`DEFAULT_REUSE_TOOL_OUTPUT` (``off``), which is to say
        "no behaviour changes unless a run asks for it". See
        :data:`REUSE_TOOL_OUTPUT_MODES` for what each value does.

        A missing section yields the DEFAULT rather than a ``KeyError``, for the
        same reason `variants_mpileup_min_ireads` does: an annotation-less
        config is a legitimate STUB/TEST run, not a malformed one.

        An unrecognised value raises rather than falling back to ``off``. The
        failure it prevents is not the loud kind: a typo like ``reuse`` or
        ``reuse: true`` that resolved to ``off`` would produce a run that
        re-annotated everything and reported no reason why, which reads exactly
        like a correctly configured run.

        **The machine overlay may state it, for the overlays
        :data:`OVERLAY_SCIENCE_EXEMPT_KEYS` permits.** "Is there verified
        pre-computed annotation on this machine?" is a machine fact, exactly like
        ``runtime.allow_real_mode``; what the run does with the records is science
        and stays in science.yaml. The overlay's value WINS over the science
        file's when it states one, because that is the whole point of the
        exemption - and ``MachineConfig.reuse_tool_output`` is ``None`` rather
        than ``off`` precisely so "the overlay said nothing" is distinguishable
        from "the overlay chose off".

        Returns:
            One of :data:`REUSE_TOOL_OUTPUT_MODES`.
        """
        section = (self.raw or {}).get("annotation") or {}
        value = section.get("reuse_tool_output", DEFAULT_REUSE_TOOL_OUTPUT)
        if self.machine is not None and self.machine.reuse_tool_output is not None:
            value = self.machine.reuse_tool_output
        if value is None:
            return DEFAULT_REUSE_TOOL_OUTPUT
        if isinstance(value, bool):
            # Not a defensive branch: YAML 1.1 reads a bare `off`, `no`, `on`,
            # `yes` and `false` as booleans, so the natural spelling of the
            # default arrives here as `False`. The config file quotes it, but a
            # machine overlay or a hand edit need not, and the error they would
            # otherwise get - "False is not one of off, prefer, require" - reads
            # as though the pipeline had invented a rule.
            raise ConfigError(
                f"annotation.reuse_tool_output is the boolean {value!r}. YAML "
                "reads a bare off/no/on/yes/false as a boolean, so the default "
                'has to be written as \'off\' with quotes. Set it to "off", '
                '"prefer" or "require", or remove the key.',
                annotation_reuse_tool_output=value,
                accepted=list(REUSE_TOOL_OUTPUT_MODES),
            )
        value = str(value)
        if value not in REUSE_TOOL_OUTPUT_MODES:
            raise ConfigError(
                f"annotation.reuse_tool_output is {value!r}, which is not one of "
                f"{', '.join(REUSE_TOOL_OUTPUT_MODES)}. Set it to one of those "
                "three, or remove the key to take the default ('off'). A value "
                "this function does not recognise is refused rather than read as "
                "'off': a silent fallback would re-annotate the whole cohort and "
                "report nothing about why, which is indistinguishable from a run "
                "that was configured correctly.",
                annotation_reuse_tool_output=value,
                accepted=list(REUSE_TOOL_OUTPUT_MODES),
            )
        return value

    def data_root(self, mode: RunMode) -> Path:
        """Return the data root for ``mode``.

        TEST and REAL resolve to *different* directories by construction so
        synthetic and real data can never be mixed.
        """
        if mode is RunMode.TEST:
            return self.root / str(self.paths.get("test_data_root", "test_data"))
        return self.root / str(self.paths.get("data_root", "data"))

    def phylogeny_dir(self, mode: RunMode) -> Path:
        """Where stage 10 reads its inputs and writes its tree.

        Not `assembly_root / "phylogeny"`, which is what `run.py` used to build.
        That conflated an output location with the assembly root: under TEST it
        landed on the committed fixture by coincidence, and under REAL it wrote
        into the read-only `data/` tree - and once `data_root` became
        `assembly_root` to bound the sequence-reading stages, into
        `db/smoke_genomes/phylogeny`, among the curated FASTA files.

        So the two are asked separately. TEST keeps the committed fixture, which
        is the existing contract; REAL gets a directory under the run's
        intermediate root, which is writable and belongs to the run.
        """
        if mode is RunMode.TEST:
            return self.data_root(mode) / "phylogeny"
        return self.intermediate_root(mode) / "phylogeny"

    def assembly_root(self, mode: RunMode) -> Path:
        """Where *this run's* assemblies are read from.

        ``data_root`` for an ordinary run. For a smoke overlay it is the overlay's
        own genome directory instead - the same directory
        :func:`papipeline.run.discover_run_manifest` builds the cohort's
        ``assembly_path`` entries from, and the same one the sample cap is
        counted against.

        This exists because a stage that reached for ``data_root`` directly would
        resolve *every* isolate in the roster that happens to have a sequence in
        ``data/`` - 602 of them - rather than the ten the overlay authorises. The
        cohort is deliberately the whole roster so unprepared isolates can refuse
        per sample; that only works if the stages that need a sequence read from
        the bounded directory. Reading from ``data/`` does not widen the cohort,
        it widens the *work*, and a bounded run then annotates sixty times more
        genomes than it claims.

        So: one accessor, so the stage and the manifest builder cannot disagree
        about which directory the run is allowed to read.
        """
        if self.is_smoke_overlay():
            return self.machine.smoke_genome_dir()
        return self.data_root(mode)

    def results_root(self, mode: RunMode) -> Path:
        """Where a run's outputs go, honouring the redirect override.

        The override is an environment variable rather than a Snakemake
        ``--config`` key because every rule shells out to a *separate* process
        that recomputes its own paths from this config. A ``--config`` value
        would therefore be seen by the workflow and not by the stage runner, so
        the run would quietly write to the real tree while the DAG believed it
        was writing somewhere else - which is how a test ends up passing and the
        developer's results directory ends up full of stub output.

        Read in exactly one place so the workflow and the worker cannot
        disagree about where a file is going.
        """
        override = os.environ.get(RESULTS_ROOT_ENV)
        if override:
            return Path(override) / mode.value.lower()
        return self.root / str(self.paths.get("results_root", "results")) / mode.value.lower()

    def reports_root(self, mode: RunMode) -> Path:
        """Reports follow the results redirect when one is set."""
        if os.environ.get(RESULTS_ROOT_ENV):
            return self.results_root(mode) / "reports"
        return self.root / str(self.paths.get("reports_root", "reports"))

    def intermediate_root(self, mode: RunMode) -> Path:
        return self.results_root(mode) / "intermediate"

    def tool_output_root(self, mode: RunMode) -> Path:
        """Where stage input tables are read from.

        In TEST mode these are the checked-in synthetic fixtures under
        ``test_data/intermediate/``. In REAL mode the external tools
        (Bakta, AMRFinderPlus, mlst, ...) write into the run's intermediate
        directory and the stages read them back from there.

        Keeping this distinct from :meth:`intermediate_root` is what stops a
        run from reading its own outputs as if they were inputs, and stops
        synthetic fixtures from being written into ``results/``.
        """
        if mode is RunMode.TEST:
            return self.data_root(mode) / "intermediate"
        return self.intermediate_root(mode)

    def metadata_dir(self, mode: RunMode) -> Path:
        return self.data_root(mode) / "metadata"

    def phenotype_dir(self, mode: RunMode) -> Path:
        """Where phenotype tables live for this run.

        A smoke overlay has its own directory, because `data/phenotype/` is the
        full cohort's: it is empty, and stage 11 would otherwise either refuse
        for want of a table or read one describing 584 isolates while reporting
        on 10. The overlay already keeps its genomes separate for the same
        reason, and the phenotype half follows it.
        """
        if self.is_smoke_overlay():
            return self.machine.smoke_phenotype_dir()
        return self.data_root(mode) / "phenotype"

    # -- knowledge lookups -----------------------------------------------
    def mechanism_for_gene(self, gene: str) -> MechanismSpec:
        """Look up a gene in the mechanism table.

        Raises:
            UnknownGeneError: The gene is absent. Unknown genes are surfaced,
                never silently passed through, so the table can be extended
                deliberately.
        """
        from ..errors import UnknownGeneError

        key = gene.strip()
        if key in self.mechanisms:
            return self.mechanisms[key]
        # Determinants are often written as gene~allele or with a source tag.
        base = key.split("~")[0].strip()
        if base in self.mechanisms:
            return self.mechanisms[base]
        raise UnknownGeneError(
            "Gene is not present in config/mechanisms.tsv",
            gene=key,
            hint="Add a row to config/mechanisms.tsv before rerunning",
        )

    def has_gene(self, gene: str) -> bool:
        base = gene.split("~")[0].strip()
        return base in self.mechanisms

    def regulator(self, gene: str) -> Optional[RegulatorSpec]:
        return self.regulators.get(gene.strip())

    def require_antibiotic(self, antibiotic: str) -> str:
        """Validate an antibiotic against both science.yaml and antibiotics.tsv."""
        if antibiotic not in self.antibiotics:
            raise AntibioticError(
                "Antibiotic is not configured for this project",
                antibiotic=antibiotic,
                configured=",".join(self.antibiotics),
            )
        if antibiotic not in self.antibiotic_specs:
            raise AntibioticError(
                "Antibiotic is in science.yaml but missing from config/antibiotics.tsv",
                antibiotic=antibiotic,
            )
        return antibiotic

    def unpinned_references(self) -> List[ReferenceSpec]:
        return [r for r in self.references.values() if not r.is_pinned]

    def references_path(self) -> Path:
        """The pinned-reference table itself.

        Exposed because a stage preflight has to re-read the file from disk, not
        the specs already parsed from it: the preflight's job is to compare what
        a database reports about itself against what the contract says, and a
        cached parse would make the comparison partly circular. Without this,
        each preflight reconstructs `config/references.tsv` by hand and they can
        drift.
        """
        return self.root / "config" / "references.tsv"


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------


def _resolve_root(config_path: Optional[Path]) -> Path:
    if config_path is None:
        return Path.cwd()
    config_path = Path(config_path).resolve()
    if not config_path.exists():
        raise ConfigError("Config file not found", path=str(config_path))
    # The pipeline root is the parent of the config/ directory.
    return config_path.parent.parent


def load_antibiotics(path: Path) -> Dict[str, AntibioticSpec]:
    """Load and validate ``config/antibiotics.tsv``."""
    rows = read_tsv(
        path,
        required_columns=("antibiotic", "antibiotic_class"),
        unique_columns=("antibiotic",),
    )
    specs: Dict[str, AntibioticSpec] = {}
    for row in rows:
        name = row["antibiotic"]
        if not name:
            raise ConfigError(
                "Antibiotic row has an empty antibiotic name", path=str(path)
            )
        specs[name] = AntibioticSpec(
            antibiotic=name,
            antibiotic_class=row["antibiotic_class"] or "unknown",
            mechanism_scope=row.get("mechanism_scope") or "unknown",
            phenotype_source_standard=row.get("phenotype_source_standard") or "unknown",
            notes=row.get("notes"),
        )
    return specs


def load_mechanisms(path: Path) -> Dict[str, MechanismSpec]:
    """Load and validate ``config/mechanisms.tsv``."""
    required = (
        "mechanism",
        "gene",
        "gene_type",
        "antibiotic",
        "biological_role",
        "evidence_level",
        "mechanism_class",
        "gene_role",
        "claim_ceiling",
    )
    rows = read_tsv(path, required_columns=required)
    specs: Dict[str, MechanismSpec] = {}
    for row in rows:
        gene = row["gene"]
        if not gene:
            raise ConfigError("Mechanism row has an empty gene", path=str(path))
        if gene in specs:
            raise ConfigError(
                "Duplicate gene in mechanisms table", path=str(path), gene=gene
            )
        specs[gene] = MechanismSpec(
            mechanism=row["mechanism"] or "unknown",
            gene=gene,
            gene_type=row["gene_type"] or "unknown",
            antibiotic=row["antibiotic"] or "all",
            biological_role=row["biological_role"] or "",
            evidence_level=_parse_claim(row["evidence_level"], path, gene, "evidence_level"),
            mechanism_class=row["mechanism_class"] or "unknown",
            gene_role=row["gene_role"] or "unknown",
            claim_ceiling=_parse_claim(row["claim_ceiling"], path, gene, "claim_ceiling"),
            notes=row.get("notes"),
        )
    return specs


def load_regulators(path: Path) -> Dict[str, RegulatorSpec]:
    """Load and validate ``config/regulators.tsv``."""
    required = (
        "gene",
        "mechanism",
        "gene_role",
        "screen_for",
        "promoter_screen",
        "reference_length",
        "evidence_level",
    )
    rows = read_tsv(path, required_columns=required, unique_columns=("gene",))
    specs: Dict[str, RegulatorSpec] = {}
    for row in rows:
        gene = row["gene"]
        if not gene:
            raise ConfigError("Regulator row has an empty gene", path=str(path))
        screen = tuple(
            s.strip() for s in (row["screen_for"] or "").replace(" ", "").split(",") if s.strip()
        )
        length_raw = row.get("reference_length")
        try:
            length = int(length_raw) if length_raw else None
        except (TypeError, ValueError):
            raise ConfigError(
                "Regulator reference_length is not an integer",
                path=str(path),
                gene=gene,
                value=str(length_raw),
            ) from None
        specs[gene] = RegulatorSpec(
            gene=gene,
            mechanism=row["mechanism"] or "unknown",
            gene_role=row["gene_role"] or "unknown",
            screen_for=screen,
            promoter_screen=str(row["promoter_screen"]).strip().lower()
            in {"true", "yes", "1"},
            reference_length=length,
            evidence_level=_parse_claim(row["evidence_level"], path, gene, "evidence_level"),
            notes=row.get("notes"),
            locus_tag=(row.get("locus_tag") or "").strip() or None,
        )
    return specs


def load_references(path: Path) -> Dict[str, ReferenceSpec]:
    """Load ``config/references.tsv``."""
    required = ("reference_id", "tool", "tool_version", "database", "database_version", "version_status")
    rows = read_tsv(path, required_columns=required, unique_columns=("reference_id",))
    specs: Dict[str, ReferenceSpec] = {}
    for row in rows:
        rid = row["reference_id"]
        status = (row["version_status"] or "unpinned").strip().lower()
        if status not in {"pinned", "unpinned", "na"}:
            raise ConfigError(
                "Unknown version_status in references table",
                path=str(path),
                reference_id=rid,
                value=status,
            )
        specs[rid] = ReferenceSpec(
            reference_id=rid,
            tool=row.get("tool"),
            tool_version=row.get("tool_version"),
            database=row.get("database"),
            database_version=row.get("database_version"),
            version_status=status,
            source=row.get("source"),
            notes=row.get("notes"),
        )
    return specs


def _parse_claim(value: Optional[str], path: Path, gene: str, column: str) -> ClaimStatus:
    """Parse a claim status from a knowledge table, rejecting 'causal'."""
    text = (value or "").strip().upper()
    if text == "CAUSAL":
        raise ConfigError(
            "A causal claim is not a permitted evidence level in this pipeline",
            path=str(path),
            gene=gene,
            column=column,
            hint="Use DETECTED, PREDICTED, ASSOCIATED, SUPPORTED or UNKNOWN",
        )
    try:
        return ClaimStatus(text)
    except ValueError:
        raise ConfigError(
            "Unknown evidence level",
            path=str(path),
            gene=gene,
            column=column,
            value=str(value),
            allowed=",".join(s.value for s in ClaimStatus),
        ) from None



# --------------------------------------------------------------------------
# Machine overlays
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class MachineConfig:
    """One machine: where things live, how many resources, what can run.

    An overlay is deliberately thin. It may declare only the keys in
    :data:`MACHINE_OWNED_KEYS`; a scientific section appearing here is an
    error, which is what makes it impossible for two machines to disagree
    about the science rather than merely unlikely.
    """

    name: str
    machine_class: str
    root: Path
    paths: Mapping[str, Any]
    threads: int
    memory_mb: int
    allow_real_mode: bool
    log_level: str
    max_samples: Optional[int]
    observatory: Mapping[str, Any]
    unique_patterns_memory_mb: Optional[int]
    unique_patterns_cores: Optional[int]
    unique_patterns_temp_dir: Optional[str]
    tools: Mapping[str, Mapping[str, Any]]
    reference: Mapping[str, Any]
    source: Path

    #: The cohort bound THIS overlay states, for the overlays permitted to state
    #: one (see :data:`OVERLAY_SCIENCE_EXCEPTIONS`). Default -- the whole roster,
    #: i.e. whatever `config/science.yaml` says -- for every overlay that is not
    #: exempt. Kept as its own field rather than folded into `paths` so that a
    #: reader of the overlay can see it is a cohort claim and not a directory.
    cohort: CohortConfig = CohortConfig()

    #: The overlay's ``annotation.reuse_tool_output``, for the overlays allowed
    #: to declare it (see :data:`OVERLAY_SCIENCE_EXEMPT_KEYS`). ``None`` means
    #: "this overlay says nothing", which is NOT the same as ``off``: an overlay
    #: that says nothing inherits science.yaml's value, and an overlay that says
    #: ``off`` has chosen the default deliberately. Keeping them distinct is why
    #: this is ``None`` rather than :data:`DEFAULT_REUSE_TOOL_OUTPUT`.
    reuse_tool_output: Optional[str] = None

    # -- tool availability -------------------------------------------------
    def tool_available(self, tool: str) -> Optional[bool]:
        """Availability of ``tool`` on this machine.

        Returns ``None`` when the tool is not declared at all. An undeclared
        tool is *unknown*, which is deliberately distinct from both
        "available" and "missing": a stage that needs an undeclared tool
        should surface that the config does not know, rather than silently
        concluding the tool is present.
        """
        entry = self.tools.get(tool)
        if entry is None:
            return None
        value = entry.get("available")
        return None if value is None else bool(value)

    def tool_database(self, tool: str) -> Optional[str]:
        entry = self.tools.get(tool) or {}
        return entry.get("database")

    @property
    def tool_search_dirs(self) -> Tuple[str, ...]:
        """Extra directories to search, declared by this machine's overlay.

        Deliberately empty by default. A machine that needs a tool outside
        ``PATH`` - a hand-built Gubbins, a tool from a source build - says so
        here, so the location is a fact of the machine rather than a literal
        buried in a script. That is the difference between a pre-flight check
        that reports the truth and one that reports what a developer once
        typed.

        These are the values **as declared**, verbatim. A relative entry is
        resolved by `resolved_tool_search_dirs`, not here, so that a test can
        assert what the overlay says without also asserting where the
        repository happens to be checked out.
        """
        raw = self.paths.get("tool_search_dirs") or ()
        if isinstance(raw, str):
            raw = (raw,)
        return tuple(str(item) for item in raw)

    def resolved_tool_search_dirs(self) -> Tuple[str, ...]:
        """`tool_search_dirs` made absolute against the repository root.

        **Why this exists.** Joining a relative entry against the process
        working directory makes tool discovery depend on where you happen to be
        standing. A `tools/gubbins/bin` entry resolved against the CWD works
        under pytest, which runs from the repository root, and silently resolves
        to nothing from a Snakemake run directory, a cron job, or a git
        worktree checked out elsewhere - and a tool that resolves to nothing is
        reported as absent, not as misconfigured. That is the same
        "a missing tool reported as absent-by-omission" failure
        `tests/unit/test_tool_discovery.py` was written to prevent, and it is
        worse than a hard-coded path because it looks correct.

        Resolving against `self.root` is what makes a repository-relative entry
        mean the same thing in every checkout, which is also what lets the
        built artifact stop depending on the absolute prefix it was installed
        with. Absolute entries are passed through untouched: a Homebrew prefix
        like `/opt/homebrew/opt/mummer/bin` is a machine fact and must not be
        reinterpreted.
        """
        resolved: List[str] = []
        for directory in self.tool_search_dirs:
            path = Path(directory)
            resolved.append(str(path if path.is_absolute() else (self.root / path)))
        return tuple(resolved)

    def tool_candidates(self, tool: str) -> Tuple[str, ...]:
        """Every location to try for ``tool``, ``PATH`` first.

        ``PATH`` is authoritative and comes first. Configured search
        directories follow, made absolute by `resolved_tool_search_dirs`.
        Nothing else is consulted, and a directory that does not exist is
        skipped rather than returned as a dead path.
        """
        import shutil

        seen: List[str] = []
        resolved = shutil.which(tool)
        if resolved:
            seen.append(resolved)
        for directory in self.resolved_tool_search_dirs():
            candidate = str(Path(directory) / tool)
            if candidate in seen:
                continue
            if not Path(candidate).exists():
                continue
            seen.append(candidate)
        return tuple(seen)

    # -- derived paths -----------------------------------------------------
    def _path(self, key: str, default: str) -> Path:
        return self.root / str(self.paths.get(key, default))

    def data_root_for(self, mode: RunMode) -> Path:
        if mode is RunMode.TEST:
            return self._path("test_data_root", "test_data")
        return self._path("data_root", "data")

    def results_root_for(self, mode: RunMode) -> Path:
        return self._path("results_root", "results") / mode.value.lower()

    def reports_root_for(self, mode: RunMode) -> Path:
        return self._path("reports_root", "results/reports")

    def intermediate_root_for(self, mode: RunMode) -> Path:
        return self.results_root_for(mode) / "intermediate"

    def status_dir(self) -> Path:
        return self._path("status_dir", "status")

    def db_root(self) -> Path:
        return self._path("db_root", "db")

    def reference_fasta(self, science: Optional[Mapping[str, Any]] = None) -> Path:
        """Where the pinned reference lives on this machine.

        The *identity* of the reference comes from ``science.yaml``; only the
        path is per-machine. Taking it from the overlay alone would let two
        machines disagree about which genome the science was run against.
        """
        identity = (science or {}).get("reference") or {}
        name = identity.get("accession") or "reference"
        if not identity:
            LOGGER.warning(
                "No reference section in the science configuration; cannot "
                "confirm which reference genome this run uses"
            )
        return self._path(
            "reference_fasta", f"db/reference/{name}/{name}_genomic.fna"
        )

    # -- the safety rail ----------------------------------------------------
    def smoke_genome_dir(self) -> Path:
        """The smoke overlay's genome directory. Required, with no default.

        Deliberately *not* `data_root_for(RunMode.REAL)`. A bounded smoke run
        that fell back to `data/` when this key was absent would analyse the full
        cohort while appearing to be a 10-assembly run - the one outcome the
        overlay exists to prevent. So the key is required of any overlay that
        declares itself a smoke overlay, and the refusal names it.

        `load_machine_config` already refuses a smoke overlay missing this key,
        at load, so the guard below is a second line for a `MachineConfig` built
        directly rather than through the loader.

        Raises:
            ConfigError: The overlay is a smoke overlay without the key.
        """
        value = self.paths.get("smoke_genome_dir")
        if not value:
            raise ConfigError(
                "This overlay declares itself a smoke overlay but does not set "
                "`paths.smoke_genome_dir`. There is no default and no fall-back: "
                "pointing it at `data_root` would analyse the full cohort while "
                "appearing to be a bounded run. Add `smoke_genome_dir` under "
                f"`paths` in {self.source}.",
                path=str(self.source),
                key="smoke_genome_dir",
            )
        return self.root / str(value)

    def smoke_phenotype_dir(self) -> Path:
        """The smoke overlay's phenotype directory. Required, with no default.

        Same reasoning as :meth:`smoke_genome_dir`, applied to the other half of
        the input. ``data/phenotype/`` is the **full cohort's** directory and is
        empty; a bounded run's ten-isolate table written there would sit exactly
        where the real one belongs, and the two are indistinguishable by path.

        Without this, stage 11 either refuses (no phenotype anywhere) or reads
        the real cohort's table, which is the second of those two failures: a
        10-isolate run reporting on 584.

        Raises:
            ConfigError: A non-smoke overlay was asked - which is a programming
                error, not a configuration one - or a directly-built smoke
                overlay lacks the key.
        """
        if not self.paths.get("smoke_genome_dir"):
            raise ConfigError(
                "smoke_phenotype_dir is a smoke-overlay concept, but this "
                "overlay does not set `paths.smoke_genome_dir`, so it is not a "
                "smoke overlay. A non-smoke run reads `data/phenotype/` - use "
                "PipelineConfig.phenotype_dir(RunMode.REAL) instead.",
                path=str(self.source),
                keys=",".join(sorted(self.paths or {})),
            )
        value = self.paths.get("smoke_phenotype_dir")
        if not value:
            raise ConfigError(
                "This overlay declares itself a smoke overlay but does not set "
                "`paths.smoke_phenotype_dir`. There is no default and no "
                "fall-back to `data/phenotype/`: that is the full cohort's "
                "directory, and a bounded run's table written there would be "
                "indistinguishable from the real one. Add `smoke_phenotype_dir` "
                f"under `paths` in {self.source}.",
                path=str(self.source),
                key="smoke_phenotype_dir",
            )
        return self.root / str(value)

    def enforce_sample_cap(self, sample_count: int) -> None:
        """Refuse a cohort this machine cannot carry.

        Raises:
            SampleCapExceeded: ``sample_count`` exceeds ``max_samples``. The
                message names the cap, the observed count, and the file to
                change, because a refusal the user cannot act on is just an
                obstacle.
        """
        if self.max_samples is None:
            return
        if sample_count <= self.max_samples:
            return
        raise SampleCapExceeded(
            f"This machine is configured for at most {self.max_samples} samples, "
            f"but the run was given {sample_count}. "
            f"Raise or remove `runtime.max_samples` in {self.source} to run a "
            f"larger cohort here; for the full-scale run use the analysis "
            f"machine (config/machines/bigmachine.yaml), which has no cap.",
            max_samples=self.max_samples,
            requested=sample_count,
            config=str(self.source),
        )


def _overlay_root(overlay_path: Path) -> Path:
    """Where the paths a machine overlay declares are rooted.

    An overlay declares machine facts - how many samples this machine can hold,
    where its reference lives. It does not declare *where* it should be rooted,
    and rooting it at its own file location makes that an accident of where the
    file happens to sit.

    So: an overlay inside a real project layout (``<root>/config/machines/*.yaml``)
    is rooted at that project. Anything else is rooted at this repository, which
    is what a named overlay resolves to anyway. Without this, passing an overlay
    by path rooted it at ``/`` - ``path.parent.parent.parent`` of a file in a
    temporary directory - and every derived path became ``/data``,
    ``/results``. That was invisible while no production caller used the
    delegating accessors, and it is not a dormant problem: the bounded REAL
    smoke run is meant to use a custom overlay, which makes a path-loaded
    overlay the normal case.
    """
    path = Path(overlay_path)
    candidate = path.parent.parent.parent
    # Only trust the ancestry if it really is a project root, i.e. it has the
    # config/ directory this overlay came out of.
    if (candidate / "config" / MACHINES_DIRNAME).is_dir():
        return candidate
    return Path(__file__).resolve().parents[2]


#: The one deliberate exception to "an overlay may not declare a science-owned
#: section", keyed by the overlay's own ``machine:`` name and then by the single
#: section it may add.
#:
#: WHY IT EXISTS. `cohort.subset_file` is a science key, so
#: `config/science.yaml` states it as `null` -- the whole roster -- and the
#: laptop and the analysis machine inherit that. But the bounded REAL smoke run
#: is not a cohort-wide analysis, and "the whole roster" is the single outcome
#: that overlay exists to prevent. Its ten isolates have to be named somewhere
#: that does not also bound the laptop, which rules out science.yaml.
#:
#: WHY IT IS THIS NARROW. Three limits, each of them asserted:
#:
#:   * only the `smoke` overlay is exempt, named here rather than inferred from
#:     a path, so a second overlay cannot acquire the exemption by being added;
#:   * only `cohort` is exempted, and its contents are then validated by
#:     CohortConfig.from_dict, so the door is one section wide and not one key;
#:   * `laptop` and `bigmachine` are still refused, which is the invariant the
#:     guard exists for: those two machines must not be able to disagree about
#:     the science.
#:
#: The overlay already bounds the cohort in a machine-owned key
#: (`runtime.max_samples: 10`), so this adds a second machine-owned statement of
#: the same bound rather than a new kind of claim.
OVERLAY_SCIENCE_EXCEPTIONS: Mapping[str, frozenset] = {
    "smoke": frozenset({"cohort"}),
}

#: The one key of a science-owned section the `smoke` overlay may declare, and
#: why it is a machine fact wearing a science key's name.
#:
#: WHY IT EXISTS. `annotation.reuse_tool_output` decides whether stage 2 runs
#: Bakta or IMPORTS Bakta output already on disk. "Is there verified pre-computed
#: annotation on THIS machine?" is a machine fact - exactly like
#: `runtime.allow_real_mode` and `runtime.max_samples`, which the overlay already
#: owns beside it. What the run does with the resulting records is science and
#: stays in science.yaml; whether this machine may reuse instead of recompute is
#: not a scientific claim and must not become one machine's.
#:
#: WHY IT IS ONE KEY AND NOT A SECTION EXEMPTION. `cohort` above is exempted
#: wholesale because its validator (`CohortConfig.from_dict`) covers the whole
#: section and there is nothing else in it. `annotation` has other keys, and
#: exempting the section would let this overlay set any of them - which is the
#: hole :data:`OVERLAY_SCIENCE_EXCEPTIONS` exists to prevent. So the exemption is
#: keyed, and every other key in the block is still refused by name.
#:
#: `require` on the smoke overlay is what makes a Bakta-free REAL run possible at
#: all: it is the only mode in which the tool is NEVER invoked, so a run against
#: pre-computed annotation cannot silently re-spend 68 minutes per cohort
#: re-deriving it. `prefer` cannot answer that question - when nothing verifies it
#: quietly runs the tool over a cohort assembled never to be annotated, which is
#: not a fallback but a second way to get the wrong answer.
OVERLAY_SCIENCE_EXEMPT_KEYS: Mapping[str, Mapping[str, frozenset]] = {
    "smoke": {"annotation": frozenset({"reuse_tool_output"})},
}


def _validate_exempt_sections(
    machine_name: str, raw: Mapping[str, Any], path: Path
) -> None:
    """Type-check the science sections an overlay is permitted to declare.

    An exemption that skipped validation would be a hole with a comment on it:
    the overlay could then say ``subset_file: ""`` and reach the run without
    ever passing the refusal that the science file's own value would get. So the
    exempt block goes through exactly the class the science file uses.
    """
    validators = {"cohort": CohortConfig.from_dict}
    for section in sorted(set(OVERLAY_SCIENCE_EXCEPTIONS.get(machine_name, ()))):
        if section not in raw:
            continue
        block = raw[section]
        if not isinstance(block, Mapping):
            raise ConfigError(
                f"{machine_name!r} declares {section!r} as a section, but the "
                "value is not a mapping of its keys",
                path=str(path),
            )
        validators[section](block)


def load_machine_config(overlay: "str | Path") -> MachineConfig:
    """Load and validate one machine overlay.

    Args:
        overlay: A machine name (``"laptop"``) or a path to a ``*.yaml`` in
            ``config/machines/``.

    Raises:
        ConfigError: The file is missing or malformed, omits its machine name
            or its resources, or declares a section that belongs to science.
    """
    path = _resolve_overlay(overlay)

    with path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)
    if not isinstance(raw, dict):
        raise ConfigError("Machine overlay root must be a mapping", path=str(path))

    science_keys = _science_top_level_keys(path)
    declared = set(raw)
    machine_name = str(raw.get("machine") or "").strip()
    exempt = OVERLAY_SCIENCE_EXCEPTIONS.get(machine_name, frozenset())
    # Keyed exemptions widen `exempt` by one SECTION but only by one KEY inside
    # it - see OVERLAY_SCIENCE_EXEMPT_KEYS. The keys are checked here, before
    # the section is allowed at all, so a typo names itself rather than being
    # swallowed by a section-level exemption.
    exempt_keys = OVERLAY_SCIENCE_EXEMPT_KEYS.get(machine_name, {})
    keyed = {
        section
        for section, keys in exempt_keys.items()
        if section in declared and set(raw.get(section) or {}) & set(keys)
    }
    exempt = frozenset(set(exempt) | keyed)
    for section, keys in sorted(exempt_keys.items()):
        block = raw.get(section)
        if not isinstance(block, dict):
            continue
        extra = sorted(set(block) - set(keys))
        if extra:
            raise ConfigError(
                f"The {machine_name!r} machine overlay may declare only "
                f"{', '.join(sorted(keys))} under {section!r}; it also declares "
                f"{', '.join(extra)}, which is science and belongs in "
                f"{_science_path_for(path)} so that every machine shares one "
                "science and cannot disagree about it.",
                path=str(path),
                science_only_keys=f"{section}.{','.join(extra)}",
            )
    collisions = sorted((declared & science_keys) - set(exempt))
    if collisions:
        raise ConfigError(
            "A machine overlay may not declare sections that belong to the "
            f"science configuration ({', '.join(collisions)}). Move them to "
            f"{_science_path_for(path)} so that every machine shares one "
            f"science and cannot disagree about it.",
            path=str(path),
            science_only_keys=",".join(collisions),
        )
    if exempt & (declared & science_keys):
        # Validated through the same class the science file uses, so an exempted
        # section cannot carry a value the science file would have refused.
        _validate_exempt_sections(machine_name, raw, path)

    for tool, entry in (raw.get("tools") or {}).items():
        if not isinstance(entry, dict):
            raise ConfigError(
                f"Tool entry {tool!r} in a machine overlay must be a mapping",
                path=str(path),
            )
        extra = sorted(set(entry) - TOOL_ENTRY_KEYS)
        if extra:
            raise ConfigError(
                f"Tool {tool!r} in a machine overlay declares "
                f"{', '.join(extra)}. A machine overlay may only say whether a "
                f"tool is available here and where it is; tool versions and "
                f"database releases are pinned in config/references.tsv, and "
                f"duplicating them per machine is how two machines end up "
                f"reading different databases.",
                path=str(path),
                tool=tool,
                disallowed=",".join(extra),
            )

    allowed_here = MACHINE_OWNED_KEYS | set(exempt)
    unknown = sorted(declared - allowed_here)
    if unknown:
        raise ConfigError(
            f"Unrecognised key(s) in a machine overlay: {', '.join(unknown)}. "
            f"A machine overlay may only declare: "
            f"{', '.join(sorted(allowed_here))}.",
            path=str(path),
        )

    name = machine_name
    if not name:
        raise ConfigError(
            "A machine overlay must declare its own name under `machine:`",
            path=str(path),
        )

    # The one exempt science key this overlay may carry, if any. `None` rather
    # than a default, so "said nothing" and "said off" stay distinguishable -
    # see the field's docstring on MachineConfig.
    reuse_override = None
    annotation_block = raw.get("annotation")
    if isinstance(annotation_block, dict) and "reuse_tool_output" in annotation_block:
        reuse_override = annotation_block["reuse_tool_output"]

    runtime = raw.get("runtime")
    if not isinstance(runtime, dict):
        raise ConfigError(
            "A machine overlay must declare a `runtime:` mapping with at least "
            "`threads` and `memory_mb`",
            path=str(path),
        )
    if runtime.get("threads") is None:
        raise ConfigError(
            "A machine overlay must declare `runtime.threads`; there is no "
            "default, because a silent default is how a large-cohort run ends up "
            "single-threaded",
            path=str(path),
        )
    if runtime.get("memory_mb") is None:
        raise ConfigError(
            "A machine overlay must declare `runtime.memory_mb`",
            path=str(path),
        )

    # Step 12a's helper limits.
    #
    # Optional at load and required at use, which is the opposite of
    # `threads`/`memory_mb` above - and deliberately so. Those two bound
    # everything the pipeline does; these three bound one `sort`, in one stage,
    # and an overlay for a bounded smoke run that never reaches stage 12 should
    # not have to declare them to be valid. :class:`PipelineConfig` raises on
    # each access when they are absent, naming the key, so the run that
    # actually needs them is the one that stops; and
    # tests/unit/test_gwas_unique_patterns.py asserts every committed overlay
    # declares all three, so a deletion is still caught in CI. What is not
    # allowed is silently inheriting the helper's own defaults - 1024 Mb, 1
    # core, /tmp - which is what those accessors exist to prevent.
    max_samples = runtime.get("max_samples")
    # A smoke overlay is validated HERE, at load, not on first use. The
    # guarantee is that the overlay refuses to *start* without a genome
    # directory; a lazy check would let the overlay load and the run begin, and
    # only refuse once something asked for the path - at which point the failure
    # reads as a missing directory rather than a misconfigured overlay.
    paths_raw = dict(raw.get("paths") or {})
    is_smoke = str(raw.get("machine") or "") == "smoke"
    if is_smoke and not paths_raw.get("smoke_genome_dir"):
        raise ConfigError(
            "The smoke overlay must declare `paths.smoke_genome_dir`, and it is "
            "absent. There is no default and no fall-back to `data_root`: "
            "resolving genomes from `data/` would analyse the full cohort while "
            f"appearing to be a bounded run. Add `smoke_genome_dir` to {path}",
            path=str(path),
            key="smoke_genome_dir",
        )
    # The phenotype half of the same argument. `data/phenotype/` is the full
    # cohort's directory and is empty, so without its own a bounded run either
    # refuses at stage 11 for want of a table or reads the real cohort's and
    # reports on 584 isolates while claiming 10. Checked at load for the reason
    # given above: the overlay must not *start* misconfigured.
    if is_smoke and not paths_raw.get("smoke_phenotype_dir"):
        raise ConfigError(
            "The smoke overlay must declare `paths.smoke_phenotype_dir`, and it "
            "is absent. There is no default and no fall-back to "
            "`data/phenotype/`: that is the full cohort's directory, and a "
            "bounded run's table written there would be indistinguishable from "
            "the real one. Add `smoke_phenotype_dir` to "
            f"{path} (see scripts/smoke/build_smoke_phenotype.py)",
            path=str(path),
            key="smoke_phenotype_dir",
        )
    return MachineConfig(
        name=name,
        machine_class=str(raw.get("machine_class") or "unknown"),
        root=_overlay_root(path),
        paths=dict(raw.get("paths") or {}),
        threads=int(runtime["threads"]),
        memory_mb=int(runtime["memory_mb"]),
        allow_real_mode=bool(runtime.get("allow_real_mode", False)),
        log_level=str(runtime.get("log_level") or "INFO"),
        max_samples=None if max_samples is None else int(max_samples),
        reuse_tool_output=reuse_override,
        observatory=dict(runtime.get("observatory") or {}),
        unique_patterns_memory_mb=(
            None
            if runtime.get("unique_patterns_memory_mb") is None
            else int(runtime["unique_patterns_memory_mb"])
        ),
        unique_patterns_cores=(
            None
            if runtime.get("unique_patterns_cores") is None
            else int(runtime["unique_patterns_cores"])
        ),
        unique_patterns_temp_dir=(
            None
            if runtime.get("unique_patterns_temp_dir") is None
            else str(runtime["unique_patterns_temp_dir"])
        ),
        tools={
            str(k): dict(v or {}) for k, v in (raw.get("tools") or {}).items()
        },
        reference=dict(raw.get("reference") or {}),
        source=path,
        cohort=(
            CohortConfig.from_dict(raw["cohort"])
            if isinstance(raw.get("cohort"), Mapping)
            else CohortConfig()
        ),
    )


def _resolve_overlay(overlay: "str | Path") -> Path:
    """Turn a machine name or path into a readable overlay path."""
    if not isinstance(overlay, (str, Path)):
        raise ConfigError(f"Machine must be a name or a path, got {type(overlay).__name__}")
    candidate = Path(overlay)
    text = str(overlay)

    if text in KNOWN_MACHINES:
        # Derive the repository root from this module's own location rather
        # than from the working directory, so a machine name resolves the same
        # way regardless of where the command was invoked from.
        repo_root = REPO_ROOT
        for base in (repo_root, Path.cwd()):
            path = base / "config" / MACHINES_DIRNAME / f"{text}.yaml"
            if path.is_file():
                return path
        raise ConfigError(
            f"Machine {text!r} is known but {MACHINES_DIRNAME}/{text}.yaml was not "
            f"found under config/. Known machines: {', '.join(KNOWN_MACHINES)}",
            machine=text,
        )

    if candidate.suffix in (".yaml", ".yml") or candidate.exists():
        path = candidate if candidate.is_absolute() else (Path.cwd() / candidate)
        if not path.is_file():
            raise ConfigError("Machine overlay file not found", path=str(path))
        return path

    known = ", ".join(KNOWN_MACHINES)
    raise ConfigError(
        f"Unknown machine {candidate.name!r}. Known machines: {known}. "
        f"Pass one of those names, or a path to a {MACHINES_DIRNAME}/*.yaml file.",
        machine=str(overlay),
        known=known,
    )


def _science_path_for(overlay_path: Path) -> Path:
    """The science file an overlay is validated against.

    A sibling ``science.yaml`` wins, so a self-contained project directory is
    self-contained. Failing that, fall back to this repository's science file:
    an overlay copied elsewhere must still be rejected for redeclaring a
    scientific section, rather than becoming unchecked just because it moved.
    """
    sibling = overlay_path.parent.parent / "science.yaml"
    if sibling.is_file():
        return sibling
    return Path(__file__).resolve().parents[2] / "config" / "science.yaml"


def _science_top_level_keys(overlay_path: Path) -> set:
    """Top-level sections declared by the science configuration."""
    science = _science_path_for(overlay_path)
    if not science.is_file():
        return set()
    with science.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)
    return set(raw) if isinstance(raw, dict) else set()


def load_config(
    config_path: Optional[Path] = None,
    machine: "str | Path | None" = "laptop",
) -> PipelineConfig:
    """Load, validate and freeze the full pipeline configuration.

    Args:
        config_path: Path to the science configuration. Defaults to
            ``<cwd>/config/science.yaml``.
        machine: A machine name (``"laptop"``, ``"bigmachine"``) or a path to a
            machine overlay. Defaults to ``"laptop"``.

            The default is the *capped* machine on purpose. An unspecified
            machine should fail safe, and the laptop is the one that refuses a
            cohort it cannot carry. The analysis machine must be asked for
            explicitly.

    Raises:
        ConfigError: Any structural or cross-file inconsistency, or a machine
            overlay that redeclares a scientific section.
    """
    if config_path is None:
        config_path = Path.cwd() / DEFAULT_CONFIG_RELPATH
    config_path = Path(config_path)
    root = _resolve_root(config_path)

    with config_path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)
    if not isinstance(raw, dict):
        raise ConfigError("Config root must be a mapping", path=str(config_path))

    if machine is None:
        # Callers that build a self-contained config in a tmp dir (the test
        # suite does this) pass machine=None to mean "this file is everything".
        machine_paths, machine_runtime, machine_cfg = {}, {}, None
    else:
        inline = sorted(k for k in ("paths", "runtime") if raw.get(k) is not None)
        if inline:
            # A science file that also carries paths or resources is
            # self-contained. Silently layering a machine overlay on top of it
            # would discard the caller's own values, and silently ignoring them
            # would hide a mistake, so the ambiguity is made impossible.
            raise ConfigError(
                f"This configuration declares its own {', '.join(inline)}, so it is "
                f"self-contained; a machine overlay would override "
                f"{', '.join(inline)}. Pass machine=None if the file really is "
                f"self-contained, or remove {', '.join(inline)} from it and name a "
                f"machine (one of: {', '.join(KNOWN_MACHINES)}).",
                path=str(config_path),
                self_contained_keys=",".join(inline),
            )
        overlay = load_machine_config(machine)
        machine_cfg = overlay
        machine_paths, machine_runtime = dict(overlay.paths), {
            "threads": overlay.threads,
            "memory_mb": overlay.memory_mb,
            "allow_real_mode": _allow_real_mode_override(overlay.allow_real_mode),
            "log_level": overlay.log_level,
            "observatory": dict(overlay.observatory),
            # Step 12a's helper runs `sort` as its own process. Carried
            # explicitly rather than folded into `memory_mb`/`threads`, because
            # `sort` has an address space this pipeline does not share and a
            # run-wide limit says nothing about it. Omitted when the overlay
            # does not declare them, so the properties below still raise rather
            # than the overlay's silence becoming an inherited default.
            "unique_patterns_memory_mb": overlay.unique_patterns_memory_mb,
            "unique_patterns_cores": overlay.unique_patterns_cores,
            "unique_patterns_temp_dir": overlay.unique_patterns_temp_dir,
        }
        if overlay.max_samples is not None:
            machine_runtime["max_samples"] = overlay.max_samples

    antibiotics = tuple(raw.get("antibiotics") or ())
    if not antibiotics:
        raise ConfigError(
            "science.yaml must list at least one antibiotic", path=str(config_path)
        )
    if len(set(antibiotics)) != len(antibiotics):
        raise ConfigError(
            "science.yaml lists a duplicate antibiotic", path=str(config_path)
        )

    config_dir = root / "config"
    antibiotic_specs = load_antibiotics(config_dir / "antibiotics.tsv")
    for antibiotic in antibiotics:
        if antibiotic not in antibiotic_specs:
            raise ConfigError(
                "Antibiotic is listed in science.yaml but absent from antibiotics.tsv",
                antibiotic=antibiotic,
                path=str(config_dir / "antibiotics.tsv"),
            )
    for name in antibiotic_specs:
        if name not in antibiotics:
            LOGGER.warning(
                "Antibiotic %s is in antibiotics.tsv but not enabled in science.yaml", name
            )

    phenotype_cfg = raw.get("phenotype") or {}
    allowed_raw = phenotype_cfg.get("allowed_values") or [p.value for p in Phenotype]
    try:
        allowed = tuple(Phenotype(str(v).strip().upper()) for v in allowed_raw)
    except ValueError as exc:
        raise ConfigError(
            "Invalid phenotype value in science.yaml",
            path=str(config_path),
            error=str(exc),
            allowed=",".join(p.value for p in Phenotype),
        ) from None

    analysis = dict(raw.get("analysis") or {})

    return PipelineConfig(
        root=root,
        raw=raw,
        antibiotics=antibiotics,
        allowed_phenotypes=allowed,
        analysis=analysis,
        mechanisms=load_mechanisms(config_dir / "mechanisms.tsv"),
        regulators=load_regulators(config_dir / "regulators.tsv"),
        references=load_references(config_dir / "references.tsv"),
        antibiotic_specs=antibiotic_specs,
        organism=dict(raw.get("organism") or {}),
        qc=QcConfig.from_dict(raw.get("qc") or {}),
        gwas=GwasConfig.from_dict(raw.get("gwas") or {}),
        phylogeny=PhylogenyConfig.from_dict(raw.get("phylogeny") or {}),
        convergence=ConvergenceConfig.from_dict(raw.get("convergence") or {}),
        lineage=LineageConfig.from_dict(raw.get("lineage") or {}),
        # An overlay permitted to state a cohort bound OVERRIDES the science
        # file's `null`, and only for the overlays named in
        # OVERLAY_SCIENCE_EXCEPTIONS. Every other machine inherits science.yaml,
        # which is the invariant the guard above exists to keep.
        cohort=(
            machine_cfg.cohort
            if machine_cfg is not None
            and machine_cfg.cohort.subset_file is not None
            else CohortConfig.from_dict(raw.get("cohort") or {})
        ),
        paths=machine_paths or dict(raw.get("paths") or {}),
        runtime=machine_runtime or dict(raw.get("runtime") or {}),
        machine=machine_cfg,
    )
