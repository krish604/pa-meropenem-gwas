"""Stage 1 of the meropenem build: the cohort join and size gate.

Before any analysis stage runs, this module answers three questions with
numbers instead of assumptions:

1. **How many isolates were tested** for the antibiotic (phenotype rows read),
2. **How many of those have an assembly** (a manifest genome to analyse),
3. **How many are R, I and S** - the binary S/I/R categories the laboratory
   reported. Phenotype is binary S/I/R only; there are no MICs here, and no
   category is ever converted into one.

It then applies the configured **intermediate-isolate policy** and refuses the
run if too few resistant isolates remain to support an association scan.

The join is not reimplemented. :func:`papipeline.join.join_manifest_to_phenotype`
provides the directional contract (spec.md D3): a manifest genome with zero or
two phenotype rows is a hard failure naming the sample; a phenotype row with
no assembly is *excluded*, not failed, because the phenotype table
legitimately describes more isolates than were downloaded. The gate only adds
policy on top of that result:

* ``cohort_gate.intermediate_policy`` (default ``exclude``) decides whether
  isolates the laboratory called ``I`` *enter* the analysis cohort. They are
  always counted in the report either way, and they are never folded into
  ``R`` or ``S`` - the laboratory's call stands (spec.md D4, ruling R7).
* ``cohort_gate.min_resistant`` (default 100) is the size floor. Below it the
  gate raises :class:`CohortGateError`, whose message names the config key so
  the remedy travels with the failure. An association scan over fewer
  resistant isolates finds lineage structure, not resistance.

Both keys live under the ``cohort_gate`` section of the science configuration
(``config/science.yaml``). Nothing here reads a path, a thread count or a
threshold from anywhere else.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from .config.loader import PipelineConfig
from .errors import ConfigError, PipelineError
from .join import ExclusionReason, JoinedSample, join_manifest_to_phenotype
from .logging_utils import get_logger
from .models import Phenotype
from .stages.phenotype import interpret_mic, load_phenotype

LOGGER = get_logger("cohort_gate")

#: The section of the science configuration this gate reads.
CONFIG_SECTION = "cohort_gate"

#: The two keys, in full. The refusal message names these strings verbatim so
#: a reader can paste them into ``config/science.yaml``.
INTERMEDIATE_POLICY_CONFIG_KEY = "cohort_gate.intermediate_policy"
MIN_RESISTANT_CONFIG_KEY = "cohort_gate.min_resistant"

#: Code defaults, used when the section or a key is absent. The shipped
#: configuration states both explicitly; these exist so a stripped-down config
#: still gates instead of silently not gating.
DEFAULT_INTERMEDIATE_POLICY = "exclude"
DEFAULT_MIN_RESISTANT = 100

#: The only policies that exist. Anything else - including a value that would
#: fold ``I`` into ``R`` - is a configuration error naming the key.
INTERMEDIATE_POLICIES: Tuple[str, ...] = ("exclude", "include")


class CohortGateError(PipelineError):
    """The cohort is too small for the configured analysis.

    Carries the numbers as context (``antibiotic``, ``n_resistant``,
    ``min_resistant``, ``config_key``), which ``PipelineError`` renders into
    ``str(exc)`` - so a log line containing the error also contains the count,
    the threshold and the key to change.
    """


@dataclass(frozen=True)
class CohortGateReport:
    """What the gate counted, in one immutable record.

    ``n_resistant`` / ``n_intermediate`` / ``n_susceptible`` count rows that
    have an assembly (the rows that can enter the cohort); rows in other
    categories (``SDD``, ``ND``) are included in ``n_with_assembly`` only.
    The three category counts therefore need not sum to ``n_with_assembly``.
    """

    antibiotic: str
    intermediate_policy: str
    min_resistant: int
    n_tested: int
    n_with_assembly: int
    n_no_assembly: int
    n_resistant: int
    n_intermediate: int
    n_susceptible: int
    n_in_cohort: int
    cohort_sample_ids: Tuple[str, ...]

    def summary(self) -> str:
        """One line with the five headline counts, for logs and run reports."""
        return (
            f"Cohort gate {self.antibiotic}: tested={self.n_tested} "
            f"with_assembly={self.n_with_assembly} R={self.n_resistant} "
            f"I={self.n_intermediate} S={self.n_susceptible} "
            f"cohort={self.n_in_cohort} policy={self.intermediate_policy} "
            f"min_resistant={self.min_resistant}"
        )

    def as_record(self, *, enforced: bool) -> Dict[str, Any]:
        """The counts as a JSON-safe record for the run result and manifest.

        ``enforced`` is not derivable from the counts - it says which mode
        took this gate - and a reader of ``run_manifest.json`` needs it to
        tell a gated run from a soft one. Everything else is a straight read
        of the fields, so the manifest and the result cannot disagree.
        """
        return {
            "status": "evaluated",
            "antibiotic": self.antibiotic,
            "intermediate_policy": self.intermediate_policy,
            "min_resistant": self.min_resistant,
            "n_tested": self.n_tested,
            "n_with_assembly": self.n_with_assembly,
            "n_no_assembly": self.n_no_assembly,
            "n_resistant": self.n_resistant,
            "n_intermediate": self.n_intermediate,
            "n_susceptible": self.n_susceptible,
            "n_in_cohort": self.n_in_cohort,
            "cohort_sample_ids": list(self.cohort_sample_ids),
            "enforced": bool(enforced),
        }


def _read_gate_settings(
    config: PipelineConfig,
) -> Tuple[str, int]:
    """Read and validate the two gate keys from the science configuration.

    Raises:
        ConfigError: The section is not a mapping, the policy is not one of
            :data:`INTERMEDIATE_POLICIES`, or the minimum is not a
            non-negative integer. Every message names its config key.
    """
    section: Any = config.raw.get(CONFIG_SECTION)
    if section is None:
        section = {}
    if not isinstance(section, Mapping):
        raise ConfigError(
            f"`{CONFIG_SECTION}` must be a mapping of gate keys in the science "
            "configuration",
            config_key=CONFIG_SECTION,
        )

    policy_raw = section.get("intermediate_policy")
    policy = DEFAULT_INTERMEDIATE_POLICY if policy_raw is None else str(policy_raw)
    if policy not in INTERMEDIATE_POLICIES:
        raise ConfigError(
            f"Unknown intermediate-isolate policy {policy_raw!r}; "
            f"`{INTERMEDIATE_POLICY_CONFIG_KEY}` accepts "
            f"{', '.join(INTERMEDIATE_POLICIES)}. Intermediates are counted "
            "under either policy and are never folded into R or S.",
            config_key=INTERMEDIATE_POLICY_CONFIG_KEY,
        )

    minimum_raw = section.get("min_resistant")
    if minimum_raw is None:
        minimum = DEFAULT_MIN_RESISTANT
    elif isinstance(minimum_raw, bool) or not isinstance(minimum_raw, int):
        raise ConfigError(
            f"`{MIN_RESISTANT_CONFIG_KEY}` must be a non-negative integer, "
            f"got {minimum_raw!r}",
            config_key=MIN_RESISTANT_CONFIG_KEY,
        )
    else:
        minimum = minimum_raw
    if minimum < 0:
        raise ConfigError(
            f"`{MIN_RESISTANT_CONFIG_KEY}` must be a non-negative integer, "
            f"got {minimum}",
            config_key=MIN_RESISTANT_CONFIG_KEY,
        )
    return policy, minimum


def _apply_intermediate_policy(
    manifest_ids: Sequence[str],
    joined_by_id: Mapping[str, JoinedSample],
    calls_by_id: Mapping[str, Any],
    excluded: Mapping[str, Any],
    policy: str,
) -> List[JoinedSample]:
    """Decide cohort membership in manifest order, on top of the join's result.

    * ``exclude``: an ``I`` isolate the join kept (only possible with a
      measured MIC) is dropped; the join's own exclusion of an MIC-less ``I``
      row stands.
    * ``include``: an ``I`` isolate the join excluded as an unusable category
      is re-admitted - that exclusion encodes the *binary model's* rule
      (spec.md D4), and this key is the run's chance to say the cohort is
      wider than the model. ``ND`` is never re-admitted: not determined
      carries no measurement by definition, and ``SDD`` is outside this key's
      scope.

    Nothing here ever relabels a category; re-admitted isolates keep the call
    the laboratory made.
    """
    cohort: List[JoinedSample] = []
    for sample_id in manifest_ids:
        kept = joined_by_id.get(sample_id)
        if kept is not None:
            if policy == "exclude" and kept.call.phenotype is Phenotype.I:
                continue
            cohort.append(kept)
            continue
        if policy != "include":
            continue
        drop = excluded.get(sample_id)
        call = calls_by_id.get(sample_id)
        if (
            drop is not None
            and drop.reason is ExclusionReason.UNUSABLE_CATEGORY
            and call is not None
            and call.phenotype is Phenotype.I
        ):
            cohort.append(
                JoinedSample(
                    sample_id=sample_id,
                    call=call,
                    log2_mic=interpret_mic(call.mic),
                )
            )
    return cohort


def evaluate_cohort_gate(
    config: PipelineConfig,
    *,
    manifest_ids: Sequence[str],
    phenotype_dir: Path,
    antibiotic: str,
    enforce: bool = True,
) -> CohortGateReport:
    """Join assemblies to S/I/R calls, report the counts, enforce the gate.

    Args:
        config: Loaded pipeline configuration (science + machine).
        manifest_ids: Sample identifiers of the assemblies, in manifest order.
        phenotype_dir: Directory holding ``<antibiotic>_phenotype.tsv``.
        antibiotic: Which antibiotic to gate on. Must be configured.
        enforce: Whether ``cohort_gate.min_resistant`` refuses the run. The
            default, ``True``, is the gate doing its job: below the floor the
            scan is underpowered and the refusal names the key. The run
            entrypoint passes ``False`` for TEST, where the committed fixtures
            carry fewer resistant isolates than the floor by construction and
            their numbers are never reported as findings - so the counts are
            still taken, and only the refusal is withheld. Real-mode runs pass
            the default. ``enforce=False`` never suppresses a ``ConfigError``,
            a join failure or a missing table: those are about the data being
            wrong, not about the cohort being small, and no mode tolerates
            them.

    Returns:
        The five counts (n tested, n with assembly, n R, n I, n S) plus the
        policy, the threshold and the resulting cohort.

    Raises:
        AntibioticError: The antibiotic is not configured in the science
            configuration.
        ConfigError: A ``cohort_gate`` key is present but invalid; the message
            names the key.
        PhenotypePreflightError: The phenotype table is missing or empty.
        DataContractError: The directional join failed - a manifest genome
            with zero or two phenotype rows, or a duplicate row.
        CohortGateError: Resistant isolates fall below
            ``cohort_gate.min_resistant``; the message names the key. Raised
            only when ``enforce`` is true.
    """
    config.require_antibiotic(antibiotic)
    policy, min_resistant = _read_gate_settings(config)

    calls = load_phenotype(config, Path(phenotype_dir), antibiotic)
    manifest_ids = list(manifest_ids)
    manifest_set = set(manifest_ids)
    calls_by_id = {call.sample_id: call for call in calls}

    with_assembly = [c for c in calls if c.sample_id in manifest_set]
    n_tested = len(calls)
    n_with_assembly = len(with_assembly)
    n_resistant = sum(1 for c in with_assembly if c.phenotype is Phenotype.R)
    n_intermediate = sum(1 for c in with_assembly if c.phenotype is Phenotype.I)
    n_susceptible = sum(1 for c in with_assembly if c.phenotype is Phenotype.S)

    # The directional join, reused. Every manifest genome must have exactly
    # one phenotype row (hard failure otherwise); rows with no assembly are
    # excluded with a reason, not failed. `require_trait=False`: the phenotype
    # is categorical S/I/R and there is no continuous trait to require.
    joined = join_manifest_to_phenotype(manifest_ids, calls, require_trait=False)

    cohort = _apply_intermediate_policy(
        manifest_ids,
        {s.sample_id: s for s in joined.joined},
        calls_by_id,
        joined.excluded,
        policy,
    )

    if n_resistant < min_resistant:
        if enforce:
            raise CohortGateError(
                f"Cohort gate refused for {antibiotic}: {n_resistant} resistant "
                f"isolate(s) with an assembly, below the minimum of "
                f"{min_resistant}. The threshold is `{MIN_RESISTANT_CONFIG_KEY}` "
                "in the science configuration (config/science.yaml); raise it "
                "only for a deliberately small run, or enlarge the cohort.",
                antibiotic=antibiotic,
                n_resistant=n_resistant,
                min_resistant=min_resistant,
                config_key=MIN_RESISTANT_CONFIG_KEY,
            )
        # Not suppressed, not silent: the floor is still reported as unmet, so
        # a log of a non-enforcing run reads the same as the refusal it would
        # have been - just without stopping. The counts below carry both
        # numbers too, so the manifest says the same thing.
        LOGGER.warning(
            "Cohort gate below the minimum for %s: %d resistant isolate(s) "
            "with an assembly, below %d (%s) - not enforced in this run mode",
            antibiotic,
            n_resistant,
            min_resistant,
            MIN_RESISTANT_CONFIG_KEY,
        )

    report = CohortGateReport(
        antibiotic=antibiotic,
        intermediate_policy=policy,
        min_resistant=min_resistant,
        n_tested=n_tested,
        n_with_assembly=n_with_assembly,
        n_no_assembly=n_tested - n_with_assembly,
        n_resistant=n_resistant,
        n_intermediate=n_intermediate,
        n_susceptible=n_susceptible,
        n_in_cohort=len(cohort),
        cohort_sample_ids=tuple(s.sample_id for s in cohort),
    )
    LOGGER.info(report.summary())
    if report.n_no_assembly:
        LOGGER.info(
            "%d %s phenotype row(s) name isolates with no assembly; excluded "
            "from the cohort, not failed (the phenotype table legitimately "
            "describes more isolates than were downloaded)",
            report.n_no_assembly,
            antibiotic,
        )
    return report
