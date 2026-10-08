"""Stage 12 - imipenem GWAS.

Two engines, one interface:

* :class:`PyseerEngine` - the preferred engine. Writes pyseer's expected
  input files and invokes it.
* :class:`ReferenceEngine` - a pure-Python/SciPy implementation of a simple
  contingency-table test with Benjamini-Hochberg correction.

The reference engine exists so the *pipeline mechanics* (input construction,
binarisation, multiple-testing correction, lineage-confound flagging) are
testable on a machine without pyseer. It is **not** a substitute for
pyseer's population-structure model: it performs no kinship correction, so
its p-values are anti-conservative in a structured cohort. Every result it
produces is stamped with that caveat, and a lineage-confounded feature is
flagged rather than reported as a hit.

What this stage will not do (rules 2, 3 and 4 in docs/scientific_rules.md):

* call a feature causal;
* report an association without also reporting the lineage distribution,
  so a lineage artefact is visible rather than hidden in a p-value;
* treat ``I`` as ``R`` or ``S``. Categories in ``gwas.excluded`` are dropped
  from the analysis and counted, not silently reassigned.
"""

from __future__ import annotations

import re
import shutil
import sys

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

from ..config.loader import PipelineConfig
from ..errors import DataContractError, ToolNotAvailableError
from ..io.tsv import read_tsv, write_tsv
from ..logging_utils import get_logger
from ..manifest import SampleManifest
from ..models import GwasResult, PhenotypeCall, RunMode

LOGGER = get_logger("stages.gwas")

GWAS_COLUMNS: Tuple[str, ...] = (
    "feature",
    "feature_type",
    "effect",
    "p_value",
    "adjusted_p_value",
    "effect_size",
    "frequency",
    "lineage_distribution",
    "model",
)

FEATURE_COLUMNS: Tuple[str, ...] = ("sample_id",)

#: Minimum number of samples carrying a feature for it to be testable.
MIN_CARRIERS_DEFAULT = 2

#: Wrapper that applies the scipy/statsmodels aliases pyseer 1.1.2 expects.
_PYSEER_COMPAT = Path(__file__).resolve().parents[2] / "scripts" / "pyseer_compat.py"

#: A feature carried by more than this fraction of samples is background.
MAX_FEATURE_FREQUENCY = 0.95


@dataclass
class GwasInput:
    """The constructed GWAS inputs, ready to hand to an engine."""

    sample_ids: Tuple[str, ...]
    binary_outcome: Dict[str, int]
    features: Dict[str, Dict[str, int]]
    feature_types: Dict[str, str]
    lineages: Dict[str, str]
    n_positive: int
    n_negative: int
    n_excluded: int
    excluded_samples: Tuple[str, ...] = ()

    @property
    def n_samples(self) -> int:
        return len(self.sample_ids)

    def testable_features(self) -> List[str]:
        """Features with enough carriers and non-carriers to test."""
        out: List[str] = []
        for name, vector in self.features.items():
            carriers = sum(vector[s] for s in self.sample_ids)
            if MIN_CARRIERS_DEFAULT <= carriers <= MAX_FEATURE_FREQUENCY * self.n_samples:
                out.append(name)
            else:
                LOGGER.debug(
                    "Skipping feature %s: %d/%d carriers (outside testable range)",
                    name,
                    carriers,
                    self.n_samples,
                )
        return sorted(out)


def parse_feature_name(name: str) -> Tuple[str, str]:
    """Split ``type__label`` into its parts.

    A feature with no ``__`` separator is treated as an unspecified type.
    """
    if "__" in name:
        feature_type, _, label = name.partition("__")
        return feature_type.strip(), label.strip()
    return "unspecified", name.strip()


def build_input(
    config: PipelineConfig,
    manifest: SampleManifest,
    phenotype_calls: Sequence[PhenotypeCall],
    feature_rows: Sequence[Mapping[str, Optional[str]]],
    lineages: Optional[Mapping[str, str]] = None,
) -> GwasInput:
    """Construct the binarised outcome and the feature matrix.

    Args:
        phenotype_calls: Stage 11 output for the target antibiotic.
        feature_rows: Stage 12 feature table (one row per sample).
        lineages: Optional sample -> lineage map, used for the lineage
            distribution and the confound check.

    Raises:
        DataContractError: A manifest sample has no usable feature row, or
            the positive and negative groups are too small to test.
    """
    gwas = config.gwas
    # The project's antibiotic, by name. `config.antibiotics[0]` is a position
    # rather than a choice, so it silently names whichever drug happens to be
    # listed first - which disagrees with `run.py`, which reads the same
    # preference from `project.primary_antibiotic`. The value reaches the
    # "too few samples in one group" error, so a failed run could name the wrong
    # drug as the one under analysis. Mirrors run.py's resolution, including its
    # fallback, so the two cannot drift apart again.
    antibiotic = str(
        (getattr(config, "raw", None) or {}).get("project", {}).get(
            "primary_antibiotic"
        )
        or config.antibiotics[0]
    )

    by_sample: Dict[str, PhenotypeCall] = {
        call.sample_id: call for call in phenotype_calls
    }

    binary: Dict[str, int] = {}
    excluded: List[str] = []
    for sample_id in manifest.sample_ids:
        call = by_sample.get(sample_id)
        if call is None:
            excluded.append(sample_id)
            continue
        value = call.phenotype.value
        if value in gwas.positive:
            binary[sample_id] = 1
        elif value in gwas.negative:
            binary[sample_id] = 0
        else:
            excluded.append(sample_id)

    n_pos = sum(binary.values())
    n_neg = len(binary) - n_pos
    if n_pos < gwas.min_samples_per_group or n_neg < gwas.min_samples_per_group:
        raise DataContractError(
            "Too few samples in one outcome group to run GWAS",
            n_positive=n_pos,
            n_negative=n_neg,
            minimum=gwas.min_samples_per_group,
            # Named so the refusal says which key to change. A count and a
            # threshold with no key is a dead end for an operator, and this is
            # the refusal a real cohort hits first: `gwas.min_samples_per_group`
            # in config/science.yaml. It is also the *only* sample-size floor
            # this stage applies, so the smallest cohort that can be analysed at
            # all is twice it - 2 * 3 = 6 with the committed value. That floor
            # is not a statement about power: a mixed model on six samples runs
            # and means nothing, which is why the result carries the
            # underpowered flag rather than being suppressed.
            minimum_key="gwas.min_samples_per_group",
            minimum_cohort_n=2 * gwas.min_samples_per_group,
            antibiotic=antibiotic,
            n_excluded=len(excluded),
        )

    rows_by_sample: Dict[str, Mapping[str, Optional[str]]] = {}
    for row in feature_rows:
        sample_id = str(row.get("sample_id") or "")
        if not sample_id:
            raise DataContractError("GWAS feature row has no sample_id")
        if sample_id in rows_by_sample:
            raise DataContractError(
                "Duplicate sample in GWAS feature table", sample_id=sample_id
            )
        rows_by_sample[sample_id] = row

    missing = [s for s in binary if s not in rows_by_sample]
    if missing:
        raise DataContractError(
            "Samples with phenotype but no GWAS feature row",
            n_missing=len(missing),
            examples=",".join(missing[:10]),
        )

    features: Dict[str, Dict[str, int]] = {}
    feature_types: Dict[str, str] = {}
    allowed_types = set(gwas.feature_types)
    warned_types: Set[str] = set()

    for sample_id in sorted(binary):
        row = rows_by_sample[sample_id]
        for column, raw in row.items():
            if column == "sample_id" or raw is None:
                continue
            text = str(raw).strip().lower()
            if text in {"1", "1.0", "true", "present"}:
                value = 1
            elif text in {"0", "0.0", "false", "absent"}:
                value = 0
            else:
                LOGGER.debug(
                    "Non-binary value %r for feature %s in %s; treated as absent",
                    raw,
                    column,
                    sample_id,
                )
                value = 0
            features.setdefault(column, {})[sample_id] = value
            feature_type, _label = parse_feature_name(column)
            if feature_type not in allowed_types and feature_type not in warned_types:
                warned_types.add(feature_type)
                LOGGER.warning(
                    "Feature prefix %r (e.g. %s) is not in "
                    "config.gwas.feature_types %s. Such features are still "
                    "tested but flagged; add the prefix to the config to "
                    "declare it.",
                    feature_type,
                    column,
                    sorted(allowed_types),
                )
            feature_types[column] = feature_type

    LOGGER.info(
        "Stage 12 input: %d samples | R-equivalent n=%d | S-equivalent n=%d | "
        "excluded %d | %d features",
        len(binary),
        n_pos,
        n_neg,
        len(excluded),
        len(features),
    )

    return GwasInput(
        sample_ids=tuple(sorted(binary)),
        binary_outcome=binary,
        features=features,
        feature_types=feature_types,
        lineages=dict(lineages or {}),
        n_positive=n_pos,
        n_negative=n_neg,
        n_excluded=len(excluded),
        excluded_samples=tuple(excluded),
    )


# --------------------------------------------------------------------------
# Step 12a - unique variant patterns
# --------------------------------------------------------------------------
#
# Split out from the association step because the multiple-testing threshold is
# the one number in stage 12 that a reader cannot check from the association
# table. It is alpha divided by the number of *distinct* variant patterns, and
# nothing in `12_gwas.tsv` says how large that number was. So it is computed in
# its own step, its inputs and its result are written to disk, and the count is
# reported.
#
# The count comes from pyseer's own helper, `scripts/count_patterns.py`,
# vendored into this repository because pyseer 1.1.2's wheel does not package
# it (the wheel declares five console scripts and it is not among them; see the
# file's provenance header). Substituting a locally written correction would
# have been the alternative, and it would have been a silent change of the
# science for a number the report presents as pyseer's.

#: pyseer's helper, vendored. Resolved from this module's own location rather
#: than from the working directory or an environment variable, so it is the
#: same file however the run was invoked.
COUNT_PATTERNS_HELPER = Path(__file__).resolve().parents[2] / "scripts" / "gwas" / "count_patterns.py"

#: Characters a path may contain when handed to the helper.
#:
#: Not a style preference. The helper interpolates BOTH its ``patterns``
#: argument and its ``--temp`` value into one ``shell=True`` string with no
#: quoting, so a path containing a space silently changes the command it runs
#: and one containing a metacharacter runs that instead. The vendored file must
#: stay byte-identical to upstream, so the caller validates the paths instead of
#: fixing the quoting. Both are checked: the patterns path is per-run and comes
#: from the work directory, and the temp directory comes from the machine
#: overlay - lower risk, but the same defect, and leaving one of the two
#: unguarded would be an asymmetry with no defensible reason behind it.
#:
#: Deliberately narrower than the base64 alphabet. The *contents* of the
#: patterns file are base64 and may hold any of those characters; this governs
#: paths, which this pipeline builds from directory names and filenames. The
#: set covers every path it produces.
_SAFE_PATH_CHARS = frozenset(
    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-./"
)

#: pyseer's own tested-variant count, from the stderr line it writes after the
#: association loop: ``'%d tested variants\n' % tested``
#: (``pyseer/__main__.py``). ``tested`` is incremented immediately before each
#: ``patterns.write(...)``, so this is the line count of the patterns file.
#:
#: Matched rather than searched for by substring, and pinned against the
#: installed package by
#: ``test_the_tested_count_string_is_pyseers_own`` - pyseer writes four
#: adjacent counts to stderr and `loaded` is `tested + prefiltered`, so picking
#: the wrong one would report a different population.
#:
#: Public because that test needs it, and because a check that fails open - a
#: pattern that silently stops matching - should be one caller can see.
PYSEER_TESTED_COUNT_PATTERN = r"^(\d+) tested variants$"

_PYSEER_TESTED_COUNT_RE = re.compile(PYSEER_TESTED_COUNT_PATTERN)

#: ``count_patterns.py``'s ``mem_adjust``: it passes ``--memory - mem_adjust``
#: to ``sort``. Below this the subtraction goes negative and ``sort`` is handed
#: ``-S <negative>M``, which it rejects with a message about buffer size rather
#: than about the configuration key that was wrong.
_HELPER_MIN_MEMORY_MB = 11

#: The helper adds ``--parallel=<n>`` to ``sort`` only when ``--cores`` is above
#: 1. That flag is not in POSIX, so whether a machine's ``sort`` has it is a
#: fact about that machine and is probed at run time rather than assumed.
_HELPER_PARALLEL_FLAG = "--parallel"

#: pyseer's association-table column names -> the schema
#: :func:`parse_pyseer_output` reads. From `pyseer/__main__.py`'s `header`
#: list plus `pyseer/utils.py::format_output`'s value order.
#:
#: Module-level rather than a literal inside `_normalise_columns` so that a
#: caller can state what the stage expects pyseer to have emitted, and a test
#: can check that expectation without scraping the function's source.
#: The columns without which a pyseer output file is not an association table.
_PYSEER_HEADER_REQUIRED = frozenset({"variant", "lrt-pvalue"})

#: The columns ``parse_pyseer_output`` requires of a *normalised* table. Named
#: so the refusal that reports a missing column reports the same set the check
#: used, rather than a second literal that could drift from it.
_PARSE_REQUIRED_COLUMNS = frozenset({"feature"})

PYSEER_COLUMN_RENAME: Mapping[str, str] = {
    "variant": "feature",
    "lrt-pvalue": "pvalue",
    "filter-pvalue": "filter_pvalue",
    "beta-std-err": "beta_stderr",
    "variant_h2": "variant_h2",
    "af": "freq",
}

#: The engine's own provenance record. Deliberately not a stage output table:
#: it describes the tool invocation rather than the biology, and it belongs with
#: the 12a artefacts in the work directory (see docs/data_contract.md, "Step
#: 12a artefacts"), not in ``intermediate/stages/`` where a consumer would look
#: for a result.
PYSEER_PROVENANCE_FILENAME = "pyseer_provenance.tsv"

#: File names step 12a writes inside the work directory.
PATTERNS_FILENAME = "patterns.txt"
UNIQUE_PATTERNS_FILENAME = "unique_patterns.tsv"
UNIQUE_PATTERNS_SUMMARY_FILENAME = "unique_patterns_summary.tsv"

UNIQUE_PATTERNS_COLUMNS: Tuple[str, ...] = ("pattern", "n_variants")

UNIQUE_PATTERNS_SUMMARY_COLUMNS: Tuple[str, ...] = (
    "metric",
    "value",
)


@dataclass(frozen=True)
class UniquePatternReduction:
    """The outcome of step 12a.

    Attributes:
        n_variants_tested: Patterns in the file, i.e. the number of variants
            pyseer actually tested. This is the *input* to the reduction, not
            its divisor - see :attr:`n_unique_patterns`. Not the number of
            variants in the input matrix either: prefiltered and
            out-of-frequency variants never reach the file.
        n_unique_patterns: Distinct patterns among them. This, not
            ``n_variants_tested``, is the number of independent tests.
        n_redundant_variants: ``n_variants_tested - n_unique_patterns``.
            The variants whose presence/absence vector duplicates another
            variant's, and which therefore carry no independent information.
        alpha: Family-wise error rate, from ``config.gwas.unique_patterns``.
        threshold: ``alpha / n_unique_patterns``, the value handed to pyseer's
            ``--lrt-pvalue``.
        helper_threshold: The threshold the vendored helper itself printed.
            Kept to cross-check ``threshold`` against upstream's own
            arithmetic rather than only against this module's.
        patterns_file: pyseer's ``--output-patterns`` output, read but not
            written here - pyseer writes it.
        survivors_file: One row per surviving unique pattern.
        summary_file: The audit record: counts, alpha, method, threshold and
            where every input came from.
        method: The correction, ``bonferroni``.
    """

    n_variants_tested: int
    n_unique_patterns: int
    n_redundant_variants: int
    alpha: float
    threshold: float
    helper_threshold: float
    patterns_file: Path
    survivors_file: Path
    summary_file: Path
    method: str

    def as_summary_rows(self) -> List[Dict[str, object]]:
        """The audit record, as rows for :data:`UNIQUE_PATTERNS_SUMMARY_COLUMNS`.

        Every number a reader needs to recompute the threshold is here,
        including the two files to diff and the helper that produced the count.
        A threshold nobody can reconstruct is a threshold nobody should trust.
        """
        return [
            {"metric": "method", "value": self.method},
            {"metric": "alpha", "value": self.alpha},
            {
                "metric": "variants_tested",
                "value": self.n_variants_tested,
            },
            {
                "metric": "unique_patterns",
                "value": self.n_unique_patterns,
            },
            {
                "metric": "redundant_variants",
                "value": self.n_redundant_variants,
            },
            {"metric": "threshold", "value": self.threshold},
            {"metric": "helper_threshold", "value": self.helper_threshold},
            {"metric": "patterns_file", "value": str(self.patterns_file)},
            {"metric": "survivors_file", "value": str(self.survivors_file)},
            {
                "metric": "count_helper",
                "value": str(COUNT_PATTERNS_HELPER),
            },
        ]


def _require_helper_safe_path(
    path: Path, what: str = "patterns file", key: Optional[str] = None
) -> Path:
    """The absolute path, or a refusal explaining why it cannot be passed.

    The helper cannot be given a quoted path, so anything outside
    :data:`_SAFE_PATH_CHARS` is refused here. The check exists because the
    failure it prevents is invisible: ``sort`` reading the wrong file yields a
    wrong count, and a wrong count yields a wrong threshold, and a wrong
    threshold yields p-values that look fine.

    Args:
        path: The path to validate.
        what: What the path is, for the error message. The helper receives two
            paths and they are configured in different places.
        key: The configuration key to change, when there is one. The patterns
            path comes from the work directory rather than from configuration,
            so it has none to name; the scratch directory has exactly one, and
            an error that did not name it would leave the reader guessing which
            of the machine overlay's keys was wrong.
    """
    resolved = Path(path).expanduser()
    if not resolved.is_absolute():
        resolved = (Path.cwd() / resolved).resolve()
    text = str(resolved)
    unsafe = sorted(set(text) - _SAFE_PATH_CHARS)
    if unsafe:
        raise DataContractError(
            f"Refusing to hand the vendored unique-pattern helper a {what} "
            "containing characters it cannot safely receive. The helper "
            "interpolates its file and scratch-directory arguments into a "
            "shell command without quoting, so a space would change the "
            "command and a metacharacter would be executed. "
            + (
                f"Move the location somewhere without these characters rather "
                f"than relaxing this check; the key to change is {key}."
                if key
                else "Move the work directory somewhere without these "
                "characters rather than relaxing this check."
            ),
            path=text,
            unsafe_characters="".join(unsafe),
            **({"key": key} if key else {}),
        )
    return resolved


def _require_parallel_sort(cores: int) -> None:
    """Refuse ``--cores > 1`` unless this machine's ``sort`` supports it.

    The helper adds ``--parallel=<n>`` only above one core.
    ``--parallel`` originated as a GNU coreutils extension and is not in the
    POSIX spec, so its presence varies by implementation - checked here rather
    than assumed, because where it is absent the run otherwise fails with a
    ``sort`` usage error that says nothing about the configuration key that
    caused it.

    (macOS is not the counter-example it might seem: ``/usr/bin/sort`` 2.3-Apple
    does accept ``--parallel``. Whether a given machine's ``sort`` does is a
    fact about that machine, which is the whole reason this probes rather than
    infers.)

    Raises:
        DataContractError: ``sort`` here has no ``--parallel``.
    """
    import subprocess

    try:
        help_text = subprocess.run(
            ["sort", "--help"],
            capture_output=True,
            text=True,
            timeout=30,
        ).stdout
    except (OSError, subprocess.SubprocessError) as exc:
        raise DataContractError(
            f"runtime.unique_patterns_cores is {cores}, so the helper would "
            f"pass --parallel={cores} to `sort`, and this machine's `sort` "
            f"could not be asked whether it supports that ({exc}). Set the key "
            f"to 1, or install a `sort` that has --parallel.",
            cores=cores,
            key="runtime.unique_patterns_cores",
        ) from None
    if _HELPER_PARALLEL_FLAG not in help_text:
        raise DataContractError(
            f"runtime.unique_patterns_cores is {cores}, so the helper would "
            f"pass --parallel={cores} to `sort`, but this machine's `sort` "
            f"does not support {_HELPER_PARALLEL_FLAG}. It originated as a GNU "
            "coreutils extension and is not in POSIX, so not every "
            "implementation has it. Set runtime.unique_patterns_cores to 1, or "
            "install a `sort` that has the flag. Lowering it costs nothing "
            "here: the count is `sort -u` over a file of 24-byte lines, which "
            "is not a sort that benefits from threads.",
            cores=cores,
            key="runtime.unique_patterns_cores",
        )
    LOGGER.info(
        "`sort` here supports %s; the helper will use %d cores",
        _HELPER_PARALLEL_FLAG,
        cores,
    )


def parse_count_patterns_output(stdout: str) -> Tuple[int, float]:
    """Read the vendored helper's two-line report.

    Without ``--threshold`` the helper prints::

        Patterns:\t<n>
        Threshold:\t<alpha/n in %.2E>

    ``%.2E`` is three significant digits, so the printed threshold is a
    cross-check on this module's arithmetic, not a value to hand pyseer. The
    count is exact.

    Raises:
        DataContractError: The output is not in that shape. A helper that
            printed something else has changed, and guessing at which line is
            which would put an unverified number into the threshold.
    """
    fields: Dict[str, str] = {}
    for line in stdout.splitlines():
        if not line.strip():
            continue
        key, separator, value = line.partition(":")
        if not separator:
            raise DataContractError(
                "The vendored unique-pattern helper printed a line with no "
                "':' separator, so its output shape has changed. Expected "
                "'Patterns:\\t<n>' and 'Threshold:\\t<value>' as written by "
                f"pyseer's scripts/count_patterns.py. Got: {line!r}",
                line=line,
            )
        if key.strip() in fields:
            raise DataContractError(
                f"The vendored unique-pattern helper printed {key.strip()!r} "
                "twice, so its output shape has changed. Refused rather than "
                "taking one of the two.",
                line=key.strip(),
            )
        fields[key.strip()] = value.strip()
    # An *extra* keyed line is refused too, not ignored. Two keys are expected;
    # a third means the helper changed, and silently dropping it is how a
    # renamed field turns into a silently-defaulted number.
    missing = sorted({"Patterns", "Threshold"} - set(fields))
    extra = sorted(set(fields) - {"Patterns", "Threshold"})
    if missing or extra:
        raise DataContractError(
            "The vendored unique-pattern helper's output shape has changed. "
            f"Expected exactly 'Patterns' and 'Threshold'; "
            f"{'missing ' + ', '.join(missing) if missing else ''}"
            f"{'unexpected ' + ', '.join(extra) if extra else ''}. Re-read "
            "pyseer's scripts/count_patterns.py before parsing it again.",
            missing=",".join(missing),
            unexpected=",".join(extra),
            printed=",".join(sorted(fields)),
        )
    try:
        count = int(fields["Patterns"])
        threshold = float(fields["Threshold"])
    except ValueError as exc:
        raise DataContractError(
            "The vendored unique-pattern helper printed a non-numeric count or "
            f"threshold ({exc}).",
            patterns=fields["Patterns"],
            threshold=fields["Threshold"],
        ) from None
    if count < 1:
        raise DataContractError(
            "The vendored unique-pattern helper counted "
            f"{count} unique patterns. A Bonferroni threshold of "
            "alpha / 0 is undefined, so this is refused rather than divided.",
            unique_patterns=count,
            patterns_file=fields["Patterns"],
        )
    return count, threshold


def count_unique_patterns(
    config: PipelineConfig, patterns_file: Path
) -> Tuple[int, float]:
    """Run pyseer's vendored helper over a patterns file.

    The helper's ``--memory``, ``--cores`` and ``--temp`` are all passed from
    the machine overlay rather than left at their upstream defaults, because
    they are facts about this machine that the overlay is where this repository
    keeps them.

    Args:
        config: The resolved run configuration.
        patterns_file: pyseer's ``--output-patterns`` output.

    Returns:
        ``(n_unique_patterns, helper_threshold)`` as the helper printed them.

    Raises:
        DataContractError: The helper is missing, produced no output, failed,
            or produced output this module cannot verify.
    """
    import subprocess

    helper = COUNT_PATTERNS_HELPER
    if not helper.is_file():
        raise DataContractError(
            "pyseer's unique-pattern helper is not at the vendored path. "
            "pyseer 1.1.2 does not package scripts/count_patterns.py, so this "
            "repository vendors it; if the file has been moved or deleted, the "
            "multiple-testing threshold has no defined source.",
            helper=str(helper),
        )
    patterns = _require_helper_safe_path(patterns_file)
    if not patterns.is_file():
        raise DataContractError(
            "No patterns file to reduce. Step 12a reads the file pyseer writes "
            "for --output-patterns, so this is a missing intermediate rather "
            "than an empty cohort - a cohort with no testable variant would "
            "have produced a file with zero lines, which is caught below with a "
            "different message.",
            patterns_file=str(patterns),
        )

    temp_dir = config.unique_patterns_temp_dir
    if not temp_dir.is_dir():
        raise DataContractError(
            "The unique-pattern helper's scratch directory does not exist, so "
            "its `sort` has nowhere to spill. Fix "
            "runtime.unique_patterns_temp_dir in the machine overlay rather "
            "than letting `sort` fall back somewhere else.",
            temp_dir=str(temp_dir),
            key="runtime.unique_patterns_temp_dir",
        )
    # The helper interpolates this into the same unquoted shell string as the
    # patterns path, so it gets the same check. See _SAFE_PATH_CHARS.
    temp_dir = _require_helper_safe_path(
        temp_dir, what="scratch directory", key="runtime.unique_patterns_temp_dir"
    )

    memory_mb = config.unique_patterns_memory_mb
    if memory_mb < _HELPER_MIN_MEMORY_MB:
        raise DataContractError(
            f"runtime.unique_patterns_memory_mb is {memory_mb}, below the "
            f"minimum of {_HELPER_MIN_MEMORY_MB}. The helper reserves 10 Mb "
            "for itself and passes the remainder to `sort -S`, so a smaller "
            "value would reach `sort` as a negative buffer size - a complaint "
            "about a buffer, naming a configuration key nobody set.",
            value=memory_mb,
            minimum=_HELPER_MIN_MEMORY_MB,
            key="runtime.unique_patterns_memory_mb",
        )

    cores = config.unique_patterns_cores
    if cores > 1:
        _require_parallel_sort(cores)

    cmd = [
        sys.executable,
        str(helper),
        str(patterns),
        # Deliberately NOT --threshold: that form prints only the threshold, at
        # three significant digits, and the count is the thing 12a exists to
        # report.
        "--alpha",
        str(config.gwas.unique_patterns.alpha),
        "--memory",
        str(memory_mb),
        "--cores",
        str(cores),
        "--temp",
        str(temp_dir),
    ]
    LOGGER.info("Counting unique variant patterns: %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise DataContractError(
            "pyseer's unique-pattern helper failed with exit code "
            f"{proc.returncode}. Its stderr is below; the usual cause is a "
            "patterns file that is not where the run thinks it is, or a "
            "scratch directory it cannot write to.",
            helper=str(helper),
            returncode=proc.returncode,
            stderr=proc.stderr.strip()[:2000],
        )
    if not proc.stdout.strip():
        raise DataContractError(
            "pyseer's unique-pattern helper printed nothing, so the count is "
            "unknown. Refusing to derive a threshold from it.",
            helper=str(helper),
            stderr=proc.stderr.strip()[:2000],
        )
    count, threshold = parse_count_patterns_output(proc.stdout)
    LOGGER.info(
        "Unique variant patterns: %d (alpha=%s, threshold=%.3E)",
        count,
        config.gwas.unique_patterns.alpha,
        threshold,
    )
    return count, threshold


def tally_patterns(patterns_file: Path) -> Tuple[Dict[str, int], int]:
    """Count patterns per distinct value, independently of the helper.

    Returns ``(pattern -> how many variants share it, total patterns)``.

    This exists to be checked against :func:`count_unique_patterns`, not to
    replace it. The helper is upstream's arithmetic and its answer is the one
    that goes into the threshold; this is a second opinion on the number of
    distinct values, and a disagreement means one of the two is reading a
    different file than the run thinks - which is exactly the failure a
    multiple-testing threshold must never survive. The count of *lines* also
    comes from here, because pyseer does not write it anywhere.
    """
    tally: Dict[str, int] = {}
    total = 0
    with Path(patterns_file).open("r", encoding="utf-8") as handle:
        for line in handle:
            value = line.strip()
            if not value:
                continue
            tally[value] = tally.get(value, 0) + 1
            total += 1
    return tally, total


def reduce_unique_patterns(
    config: PipelineConfig,
    patterns_file: Path,
    workdir: Path,
) -> UniquePatternReduction:
    """Step 12a: reduce the variant matrix to unique patterns.

    Writes two files into ``workdir``: ``unique_patterns.tsv``, one row per
    surviving distinct pattern with the number of variants sharing it, and
    ``unique_patterns_summary.tsv``, the audit record. Both are readable TSV so
    a user can check the count against `sort -u <patterns> | wc -l` themselves.

    Args:
        config: The resolved run configuration.
        patterns_file: pyseer's ``--output-patterns`` output for the pass that
            tested every eligible variant.
        workdir: Directory the outputs are written to.

    Raises:
        DataContractError: The patterns file is missing or empty, the helper
            disagrees with the independent tally, or the helper's own
            arithmetic disagrees with ``alpha / n``.
    """
    patterns_file = Path(patterns_file)
    workdir = Path(workdir)
    # Existence is checked here, before the tally opens the file, so a missing
    # intermediate is reported as the missing intermediate it is rather than as
    # an `open()` traceback pointing at the reduction step.
    if not patterns_file.is_file():
        raise DataContractError(
            "No patterns file to reduce. Step 12a reads the file pyseer writes "
            "for --output-patterns, so this is a missing intermediate rather "
            "than an empty cohort - a cohort with no testable variant would have "
            "produced a file with zero lines, which is caught below with a "
            "different message.",
            patterns_file=str(patterns_file),
        )
    tally, n_variants = tally_patterns(patterns_file)
    if n_variants == 0:
        raise DataContractError(
            "The patterns file is empty, so pyseer tested no variant and there "
            "is no unique-pattern count to reduce to. This is a cohort with no "
            "testable feature, not a failed reduction; the association step "
            "cannot report a threshold derived from nothing.",
            patterns_file=str(patterns_file),
        )

    n_unique, helper_threshold = count_unique_patterns(config, patterns_file)
    if n_unique != len(tally):
        raise DataContractError(
            "The vendored helper counted "
            f"{n_unique} unique patterns but the patterns file holds {len(tally)} "
            f"distinct values across {n_variants} patterns. The two disagreed, "
            "so one of them is not reading the file this run wrote, and the "
            "threshold derived from either would be unfounded.",
            helper_count=n_unique,
            independent_count=len(tally),
            patterns_tested=n_variants,
            patterns_file=str(patterns_file),
        )

    alpha = config.gwas.unique_patterns.alpha
    threshold = alpha / n_unique
    # The helper prints `%.2E`, i.e. three significant digits, so agreement is
    # asserted at that resolution and no closer. Anything looser would let a
    # genuinely different alpha or count pass.
    if not _agrees_to_printed_precision(threshold, helper_threshold):
        raise DataContractError(
            "alpha / n_unique_patterns "
            f"({threshold:.6E}) does not agree with the threshold the vendored "
            f"helper computed ({helper_threshold:.6E}). The helper's own "
            "arithmetic is upstream's, so a disagreement means the config and "
            "the helper were not given the same alpha, and one threshold is "
            "not the threshold this stage claims to use.",
            alpha=alpha,
            n_unique_patterns=n_unique,
            computed=threshold,
            helper=helper_threshold,
        )

    workdir.mkdir(parents=True, exist_ok=True)
    survivors_file = workdir / UNIQUE_PATTERNS_FILENAME
    # `sort -u` is byte order under LC_ALL=C, which the helper also requests, so
    # this file is in the same order the count was taken in. A reader who sorts
    # the patterns file themselves and diffs gets a match rather than a puzzle.
    write_tsv(
        survivors_file,
        ({"pattern": pattern, "n_variants": tally[pattern]} for pattern in sorted(tally)),
        UNIQUE_PATTERNS_COLUMNS,
        header_comment=[
            "Step 12a: unique variant patterns surviving the reduction.",
            f"{n_unique} distinct patterns from {n_variants} tested variants; "
            f"{n_variants - n_unique} were redundant.",
            f"threshold = alpha / unique_patterns = {alpha} / {n_unique} = "
            f"{threshold:.6E} ({_method_label(config)})",
            f"Source: {patterns_file}",
        ],
    )

    reduction = UniquePatternReduction(
        n_variants_tested=n_variants,
        n_unique_patterns=n_unique,
        n_redundant_variants=n_variants - n_unique,
        alpha=alpha,
        threshold=threshold,
        helper_threshold=helper_threshold,
        patterns_file=patterns_file,
        survivors_file=survivors_file,
        summary_file=workdir / UNIQUE_PATTERNS_SUMMARY_FILENAME,
        method=config.gwas.unique_patterns.method,
    )
    write_tsv(
        reduction.summary_file,
        reduction.as_summary_rows(),
        UNIQUE_PATTERNS_SUMMARY_COLUMNS,
        header_comment=[
            "Step 12a: the multiple-testing threshold and how it was derived.",
            "Every value needed to recompute it is here, including the file "
            "the count came from and the helper that produced it.",
        ],
    )
    LOGGER.info(
        "Step 12a: %d tested variants reduced to %d unique patterns "
        "(%d redundant) | threshold %.3E = %s / %d",
        n_variants,
        n_unique,
        n_variants - n_unique,
        threshold,
        alpha,
        n_unique,
    )
    return reduction


def _method_label(config: PipelineConfig) -> str:
    return config.gwas.unique_patterns.method


def _agrees_to_printed_precision(value: float, printed: float) -> bool:
    """Whether two floats match at three significant digits.

    The helper prints ``%.2E``, so that is the resolution its answer carries.
    Comparing at full float precision would fail on a value that is correct,
    and comparing at two decimals would pass on one that is not.
    """
    return f"{value:.2E}" == f"{printed:.2E}"


# --------------------------------------------------------------------------
# Engines
# --------------------------------------------------------------------------


class GwasEngine:
    """Interface for a GWAS engine."""

    name = "abstract"

    def run(self, gwas_input: GwasInput, config: PipelineConfig) -> List[GwasResult]:
        raise NotImplementedError


class ReferenceEngine(GwasEngine):
    """Fisher's exact test with Benjamini-Hochberg correction.

    No population-structure correction. Use for pipeline mechanics testing
    only; results carry an explicit warning in the ``model`` field.
    """

    name = "reference_fisher"

    def run(self, gwas_input: GwasInput, config: PipelineConfig) -> List[GwasResult]:
        from scipy.stats import fisher_exact

        pvalues: Dict[str, float] = {}
        raw: Dict[str, Tuple[Optional[float], Optional[float], float, Dict[str, int]]] = {}

        for feature in gwas_input.testable_features():
            vector = gwas_input.features[feature]
            carriers_pos = sum(
                vector.get(s, 0) for s in gwas_input.sample_ids if gwas_input.binary_outcome[s] == 1
            )
            carriers_neg = sum(
                vector.get(s, 0) for s in gwas_input.sample_ids if gwas_input.binary_outcome[s] == 0
            )
            n_pos, n_neg = gwas_input.n_positive, gwas_input.n_negative
            table = [
                [carriers_pos, n_pos - carriers_pos],
                [carriers_neg, n_neg - carriers_neg],
            ]
            odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
            pvalues[feature] = float(p_value)

            n_carriers = carriers_pos + carriers_neg
            frequency = n_carriers / gwas_input.n_samples if gwas_input.n_samples else None
            effect = float(odds_ratio) if odds_ratio is not None and odds_ratio != float("inf") else None

            distribution: Dict[str, int] = {}
            for sample_id in gwas_input.sample_ids:
                if vector.get(sample_id, 0):
                    lineage = gwas_input.lineages.get(sample_id, "unknown")
                    distribution[lineage] = distribution.get(lineage, 0) + 1
            raw[feature] = (effect, p_value, frequency or 0.0, distribution)

        adjusted = benjamini_hochberg(pvalues)

        results: List[GwasResult] = []
        for feature in sorted(raw, key=lambda f: (adjusted.get(f, 1.0), raw[f][1])):
            effect, p_value, frequency, distribution = raw[feature]
            results.append(
                GwasResult(
                    feature=feature,
                    feature_type=gwas_input.feature_types.get(feature, "unspecified"),
                    effect=effect,
                    p_value=p_value,
                    adjusted_p_value=adjusted.get(feature),
                    effect_size=effect,
                    frequency=frequency,
                    lineage_distribution=distribution,
                    model=f"{self.name}:NO_KINSHIP_CORRECTION",
                )
            )
        return results


class PyseerEngine(GwasEngine):
    """pyseer backend. Writes pyseer's inputs and invokes the executable.

    Two passes over pyseer, which is what pyseer's own documentation requires
    and not a choice made here. The threshold is ``alpha / n_unique_patterns``,
    and ``n_unique_patterns`` is only knowable after pyseer has decided which
    variants it will test - pyseer writes the patterns file during the fit, and
    only for variants that survived the frequency and prefilter gates. So:

    * **12a** runs pyseer with ``--output-patterns`` and both p-value
      thresholds stated as open. Every variant pyseer tests has its pattern
      recorded, because pyseer writes patterns before it evaluates its own
      p-value gate - so those flags do not narrow the file either way. The
      association table 12a produces is not read; only the patterns file is
      kept.
    * **12b** runs pyseer again with ``--lrt-pvalue`` set to the threshold 12a
      derived, and that table is the association result.

    The cost is one extra full model fit. That is inherent to the correction,
    not overhead this implementation added: pyseer 1.1.2 has no way to report
    the number of distinct patterns without fitting the model.
    """

    name = "pyseer"

    def __init__(
        self,
        executable: str,
        workdir: Path,
        snp_alignment: Optional[Path] = None,
    ) -> None:
        """
        Args:
            executable: pyseer entry point (see ``scripts/pyseer_compat.py``).
            workdir: Where inputs, the kinship matrix and results are written.
            snp_alignment: Optional core-SNP alignment. When given, kinship is
                estimated from these SNPs, which is what pyseer documents and
                resolves relatedness far better than gene presence/absence.
                With only a few hundred core genes the feature matrix
                describes relatedness very coarsely, and a coarse kinship
                inflates the estimated heritability.
        """
        self.executable = executable
        self.workdir = Path(workdir)
        self.snp_alignment = Path(snp_alignment) if snp_alignment else None
        self.kinship_source = "gene_presence_absence"
        #: Set by :meth:`run` once step 12a has run. ``None`` before then,
        #: which is a real distinction: an engine that has run has an
        #: auditable threshold and one that has not does not.
        self.reduction: Optional[UniquePatternReduction] = None
        #: The version pyseer reports about itself, resolved once by
        #: :meth:`require_executable`. ``None`` until then. A result with no
        #: version is a result from an unknown build of the tool that produced
        #: every p-value in it.
        self.executable_version: Optional[str] = None
        #: The argv of each pyseer invocation, in call order, so the exact
        #: command behind a table is recoverable from the run rather than from
        #: a reader's guess at what the code would have assembled.
        self.commands: List[List[str]] = []
        #: Where :meth:`require_executable` found pyseer. ``None`` until then.
        self.resolved_executable: Optional[str] = None
        #: The cohort size the provenance record states. Set by :meth:`run`;
        #: ``None`` before then, which is a real distinction for the same
        #: reason :attr:`reduction` is.
        self._n_samples: Optional[int] = None

    def write_inputs(self, gwas_input: GwasInput) -> Dict[str, Path]:
        """Write the phenotype, feature and kinship files pyseer expects."""
        self.workdir.mkdir(parents=True, exist_ok=True)
        phenotype_path = self.workdir / "phenotype.tsv"
        feature_path = self.workdir / "features.tsv"
        lineage_path = self.workdir / "lineages.tsv"

        write_tsv(
            phenotype_path,
            [
                {"sample": s, "phenotype": gwas_input.binary_outcome[s]}
                for s in gwas_input.sample_ids
            ],
            ["sample", "phenotype"],
        )
        columns = ["sample"] + sorted(gwas_input.features)
        write_tsv(
            feature_path,
            [
                {"sample": s, **{c: gwas_input.features[c].get(s, 0) for c in columns[1:]}}
                for s in gwas_input.sample_ids
            ],
            columns,
        )
        write_tsv(
            lineage_path,
            [
                {"sample": s, "lineage": gwas_input.lineages.get(s, "unknown")}
                for s in gwas_input.sample_ids
            ],
            ["sample", "lineage"],
        )
        return {
            "phenotype": phenotype_path,
            "features": feature_path,
            "lineages": lineage_path,
        }

    def require_executable(self) -> str:
        """Resolve pyseer and record the version it reports about itself.

        Called from :meth:`_invoke` rather than from ``__init__`` so that
        constructing an engine costs nothing and a caller who never reaches an
        invocation is never refused for a tool it did not run. The resolution
        and the version are cached: both pyseer passes would otherwise repeat
        the same filesystem search and the same subprocess.

        Refuses **by name** when pyseer cannot be found, in the shape the other
        adapters use (:class:`~papipeline.errors.ToolNotAvailableError`), rather
        than letting a missing tool surface as a ``FileNotFoundError`` from
        deep inside a subprocess whose stderr this stage redirects to a log
        nobody is watching. A caller resolves the path the way every other
        adapter does - ``config.machine.tool_candidates("pyseer")`` - and hands
        the result in.

        Returns:
            The version string pyseer printed.
        """
        if self.executable_version is not None:
            return self.executable_version

        import subprocess

        candidate = shutil.which(self.executable) or (
            self.executable if Path(self.executable).is_file() else None
        )
        if candidate is None:
            raise ToolNotAvailableError(
                f"Stage 12 (gwas) cannot find {self.executable!r}. It looks on "
                "PATH only; a machine overlay that keeps the tool elsewhere has "
                "to name it in tool_search_dirs, and the caller resolves it with "
                "config.machine.tool_candidates(...) and passes the result in. "
                "Without pyseer there is no mixed model, and this stage refuses "
                "rather than falling back to ReferenceEngine, whose p-values "
                "carry no population-structure correction.",
                tool=self.executable,
                executable=str(self.executable),
            )

        try:
            proc = subprocess.run(
                [candidate, "--version"],
                capture_output=True,
                text=True,
                timeout=60,
                check=False,
            )
        except OSError as exc:
            raise ToolNotAvailableError(
                f"Stage 12 (gwas) found {self.executable!r} at {candidate} but "
                f"could not run it ({exc}). A file that exists is not a tool "
                "that works; this is a fact about the path, not about pyseer.",
                tool=self.executable,
                executable=candidate,
            ) from exc

        text = (proc.stdout or "") + (proc.stderr or "")
        match = re.search(r"pyseer\s+(\S+)", text)
        if proc.returncode != 0 or match is None:
            raise ToolNotAvailableError(
                f"Stage 12 (gwas) ran {self.executable!r} at {candidate} but it "
                f"did not report a version (exit {proc.returncode}). Every "
                "p-value this stage would report comes from that tool, so a "
                "build that cannot say which one it is is refused rather than "
                "used unlabelled.",
                tool=self.executable,
                executable=candidate,
                returncode=proc.returncode,
                output=text.strip()[:400],
            )

        self.executable_version = match.group(1)
        self.resolved_executable = candidate
        LOGGER.info("pyseer %s at %s", self.executable_version, candidate)
        return self.executable_version

    def _base_command(
        self,
        config: PipelineConfig,
        gwas_input: GwasInput,
        paths: Mapping[str, Path],
        kinship: Path,
    ) -> List[str]:
        """The pyseer invocation both passes share.

        The frequency bounds and the ``--lmm`` choice are the same in 12a and
        12b, and they decide *which* variants are tested - so they decide the
        size of the unique-pattern set the threshold is divided by. They live in
        one place so the two passes cannot drift apart and quietly disagree
        about how many tests were performed.

        The frequency floor is ``1 / n_samples``, not a configured constant: a
        feature carried by fewer than one isolate cannot be tested, so a fixed
        floor would either admit untestable features or exclude a genuine
        single-carrier signal, depending on cohort size. ``n_samples`` of zero
        is refused rather than given a fallback, because it means the outcome
        groups were empty - which `build_input` already refuses - and a fallback
        here would quietly put a cohort-sized floor in place of that refusal.
        """
        if not gwas_input.n_samples:
            raise DataContractError(
                "Refusing to run pyseer with an empty cohort: the frequency "
                "floor is derived from the sample count, and a run with no "
                "samples has no correct floor to derive.",
                n_samples=gwas_input.n_samples,
            )
        min_af = 1.0 / gwas_input.n_samples
        threads = config.threads
        if threads is None:
            raise DataContractError(
                "runtime.threads is not set, so pyseer has no process count to "
                "be given. It is a machine fact and belongs in the machine "
                "overlay; inferring one here would be a guess about hardware "
                "nobody on this machine described.",
                key="runtime.threads",
            )
        # The interpreter runs the compatibility shim, NOT `self.executable`
        # directly, and that is not an accident. Measured on this machine's
        # pinned environment: `pyseer --phenotypes ... --lmm` dies with
        # `KeyError: 'arange'` from `scipy/__init__.py`, because pyseer 1.1.2
        # (2021) calls names that scipy and statsmodels have since removed.
        # `scripts/pyseer_compat.py` rebinds each one to the identical function
        # it used to name and changes no fitted value, so every p-value is
        # still pyseer's own.
        #
        # `self.executable` is therefore the *availability declaration*: it is
        # resolved and version-checked by `require_executable`, so a missing
        # pyseer is refused by name instead of surfacing as a traceback in a
        # log this stage redirects. Invoking it is what the shim is for.
        command = [
            sys.executable, str(_PYSEER_COMPAT),
            "--phenotypes", str(paths["phenotype"]),
            "--phenotype-column", "phenotype",
            "--pres", str(self.workdir / "features.rtab"),
            "--similarity", str(kinship),
            "--lmm",
            "--min-af", f"{min_af:.6f}",
            "--max-af", f"{1 - min_af:.6f}",
            "--cpu", str(max(1, int(threads))),
        ]
        return command

    def _invoke(
        self,
        cmd: Sequence[str],
        table: Path,
        log: Path,
        phase: str,
    ) -> None:
        """Run pyseer, keeping its stdout table and stderr log apart.

        pyseer writes the association table to stdout and progress to stderr.
        Merging them loses the table inside the log, which is how a run ends up
        reporting zero associations while logging thousands of lines.
        """
        import subprocess

        # Before anything is launched: pyseer must exist and must be able to
        # say which version it is. Both passes need this, so it is resolved
        # once and cached on the engine.
        self.require_executable()

        # Recorded here rather than in `_base_command`, because the argv that
        # runs is not the one that method returns: each pass appends its own
        # flags to it. A provenance record holding the shared prefix of two
        # different invocations describes neither.
        self.commands.append(list(cmd))
        table.parent.mkdir(parents=True, exist_ok=True)
        if table.exists():
            table.unlink()
        LOGGER.info("pyseer %s: %s", phase, " ".join(cmd))
        with table.open("w") as out, log.open("w") as err:
            proc = subprocess.run(cmd, stdout=out, stderr=err)
        if proc.returncode != 0:
            raise RuntimeError(
                f"pyseer failed during {phase} with exit code "
                f"{proc.returncode}; see {log}"
            )

    def run(self, gwas_input: GwasInput, config: PipelineConfig) -> List[GwasResult]:
        """Run pyseer's mixed model and parse the association table.

        pyseer 1.1.2 takes the population structure matrix as *input*
        (``--similarity`` for the random/mixed model, ``--distances`` for
        fixed effects); it does not derive one from the features. The
        kinship matrix is therefore built here from the feature matrix as a
        normalised Gram matrix, which is positive semi-definite by
        construction - the property the LMM requires of its covariance
        structure.

        The mixed model is not optional. This is a near-clonal
        *P. aeruginosa* cohort, so an uncorrected test reports lineage as
        significance. ``--lmm`` makes pyseer fit the association within the
        kinship structure.

        Returns:
            The parsed associations from pass 12b, and - on
            ``self.reduction`` - the unique-pattern reduction the threshold came
            from.

        Raises:
            DataContractError: From step 12a; see :func:`reduce_unique_patterns`.
            RuntimeError: pyseer itself failed.
        """
        if not gwas_input.testable_features():
            # Refused here rather than left to pyseer, and this is not tidiness.
            # With no feature to test, `_write_kinship`'s fallback builds an
            # all-zero Gram matrix; pyseer then dies inside its own eigendecomposition
            # with `numpy.linalg.LinAlgError: Eigenvalues did not converge` and
            # exit 1 (observed against pyseer 1.1.2). That is a tool traceback in
            # a log this stage redirects, in place of a statement about the
            # cohort - and the cohort is the thing that is actually wrong.
            raise DataContractError(
                "No feature in this cohort is testable, so there is nothing to "
                "associate and no threshold to derive. A feature is testable "
                f"when between {MIN_CARRIERS_DEFAULT} and "
                f"{MAX_FEATURE_FREQUENCY:.0%} of samples carry it; the feature "
                "table supplied has "
                f"{len(gwas_input.features)} feature(s) across "
                f"{gwas_input.n_samples} sample(s). A cohort that small cannot "
                "support a mixed model, and a kinship matrix built from nothing "
                "is singular - pyseer exits non-zero on it rather than reporting "
                "a null result.",
                n_samples=gwas_input.n_samples,
                n_features=len(gwas_input.features),
                min_carriers=MIN_CARRIERS_DEFAULT,
                max_frequency=MAX_FEATURE_FREQUENCY,
            )
        self._n_samples = gwas_input.n_samples
        paths = self.write_inputs(gwas_input)
        rtab = self.workdir / "features.rtab"
        kinship = self.workdir / "kinship.tsv"
        self._write_rtab(gwas_input, rtab)
        self._write_kinship(gwas_input, kinship)
        LOGGER.info("Kinship estimated from: %s", self.kinship_source)

        cmd = self._base_command(config, gwas_input, paths, kinship)

        # --- 12a: derive the multiple-testing threshold from unique patterns --
        patterns_file = self.workdir / PATTERNS_FILENAME
        if patterns_file.exists():
            patterns_file.unlink()
        LOGGER.info(
            "Step 12a: first pyseer pass, to record a pattern for every "
            "eligible variant. Its association table is not read - only the "
            "patterns file is used - but it is left on disk as "
            "pyseer_pass1.tsv rather than deleted, so a run can be inspected."
        )
        self._invoke(
            [*cmd, "--output-patterns", str(patterns_file),
             # Both p-value gates stated explicitly as open, rather than left at
             # pyseer's defaults, so a future change to those defaults cannot
             # alter how many patterns the file ends up holding.
             #
             # They are NOT what bounds the pattern set. `pyseer/__main__.py`
             # writes a pattern for every non-prefilter variant *before* it
             # evaluates `x.filter`, so the p-value gates cannot remove a
             # variant from this file at all; and both already default to 1.
             # What bounds it is `--min-af`/`--max-af` and the prefilter gate,
             # which are in the shared command above. An earlier version of
             # this comment claimed the flags were the bound. They were not.
             "--filter-pvalue", "1",
             "--lrt-pvalue", "1"],
            table=self.workdir / "pyseer_pass1.tsv",
            log=self.workdir / "pyseer_pass1.log",
            phase="pass 12a (pattern collection)",
        )
        self.reduction = reduce_unique_patterns(config, patterns_file, self.workdir)
        self._cross_check_pyseer_tested_count(self.workdir / "pyseer_pass1.log")

        # --- 12b: associate, at the threshold 12a derived ---------------------
        output = self.workdir / "pyseer_raw.tsv"
        LOGGER.info(
            "Step 12b: second pyseer pass at --lrt-pvalue %.6E "
            "(alpha %s / %d unique patterns)",
            self.reduction.threshold,
            self.reduction.alpha,
            self.reduction.n_unique_patterns,
        )
        self._invoke(
            [*cmd, "--lrt-pvalue", repr(self.reduction.threshold)],
            table=output,
            log=self.workdir / "pyseer.log",
            phase="pass 12b (association)",
        )
        normalised = self._normalise_columns(output)
        model = "pyseer:mixed:bonferroni_unique_patterns"
        results = parse_pyseer_output(
            normalised,
            gwas_input,
            model=model,
            correction=config.gwas.correction,
        )
        self.write_provenance(model, results)
        return results

    def write_provenance(self, model: str, results: Sequence[GwasResult]) -> Path:
        """Record which tool produced these p-values, and with what arguments.

        Every number in ``12_gwas.tsv`` comes from pyseer, and pyseer is a
        Python package whose behaviour has already changed once between the
        version this repository pins and the library versions beside it - which
        is why ``scripts/pyseer_compat.py`` exists at all. A result that does
        not carry the version and the argv cannot be told apart from one
        produced by a different build, so both are written next to the
        association table they explain.

        Written even when ``results`` is empty. An empty association table is a
        result, and "no association, from pyseer X, with these arguments" is
        the statement a reader needs; a run that failed leaves no file at all,
        which is a different fact.
        """
        path = self.workdir / PYSEER_PROVENANCE_FILENAME
        rows = [
            {
                "key": "pyseer_version",
                "value": self.executable_version or "unresolved",
            },
            {
                "key": "pyseer_executable",
                "value": self.resolved_executable or self.executable,
            },
            {
                "key": "pyseer_invoked_as",
                "value": " ".join([sys.executable, str(_PYSEER_COMPAT)]),
            },
            {"key": "kinship_source", "value": self.kinship_source},
            {"key": "model", "value": model},
            {"key": "n_samples", "value": str(self._n_samples)},
            {
                "key": "n_features_reported",
                "value": str(len(results)),
            },
            {
                "key": "n_unique_patterns",
                "value": str(
                    self.reduction.n_unique_patterns if self.reduction else "unset"
                ),
            },
            {
                "key": "lrt_pvalue_threshold",
                "value": (
                    repr(self.reduction.threshold) if self.reduction else "unset"
                ),
            },
        ]
        for index, command in enumerate(self.commands, start=1):
            rows.append({"key": f"command_{index}", "value": " ".join(command)})
        write_tsv(
            path,
            rows,
            ("key", "value"),
            header_comment=[
                "What produced the association table: the tool version, the "
                "path it was resolved to, and the exact argv of each pass.",
                "A result whose tool version is unknown cannot be reproduced.",
            ],
        )
        return path

    @staticmethod
    def _write_rtab(gwas_input: GwasInput, path: Path) -> None:
        """Write a roary/piggy-style .rtab: features as rows, samples as columns.

        pyseer reads the header with ``str.split()`` and then drops the first
        token, so the top-left cell must hold a non-empty label. A leading
        tab is silently swallowed by ``str.split()``, which shifts every
        column by one and makes pyseer reject each row as a header/row
        length mismatch.
        """
        samples = list(gwas_input.sample_ids)
        columns = sorted(gwas_input.features)
        with path.open("w") as handle:
            handle.write("Gene\t" + "\t".join(samples) + "\n")
            for column in columns:
                vector = gwas_input.features[column]
                handle.write(
                    column + "\t"
                    + "\t".join(str(int(vector.get(s, 0))) for s in samples) + "\n"
                )

    def _write_kinship(self, gwas_input: GwasInput, path: Path) -> None:
        """Write the kinship matrix for the mixed model.

        Core SNPs are used when a core-SNP alignment is available: identity
        by state across tens of thousands of sites is a far sharper
        description of relatedness than a gene presence/absence matrix over
        a few hundred core genes. The fallback is retained so the stage
        still runs without a phylogeny, and the source used is recorded.
        """
        import numpy as np

        samples = list(gwas_input.sample_ids)
        matrix = None

        if self.snp_alignment is not None and self.snp_alignment.exists():
            try:
                matrix = self._identity_by_state(self.snp_alignment, samples)
                self.kinship_source = f"core_snp_alignment:{self.snp_alignment.name}"
            except Exception as exc:  # noqa: BLE001
                LOGGER.warning(
                    "Could not build SNP kinship from %s (%s); falling back to "
                    "the feature matrix",
                    self.snp_alignment, exc,
                )
                matrix = None

        if matrix is None:
            columns = sorted(gwas_input.features)
            features = np.array(
                [
                    [float(gwas_input.features[c].get(s, 0)) for c in columns]
                    for s in samples
                ]
            )
            if features.size == 0:
                features = np.zeros((len(samples), 1))
            matrix = (features @ features.T) / features.shape[1]
            self.kinship_source = "gene_presence_absence"

        with path.open("w") as handle:
            handle.write("sample\t" + "\t".join(samples) + "\n")
            for index, sample in enumerate(samples):
                handle.write(
                    sample + "\t"
                    + "\t".join(f"{float(v):.6f}" for v in matrix[index]) + "\n"
                )

    @staticmethod
    def _identity_by_state(alignment: Path, samples: Sequence[str]):
        """Identity-by-state proportion from a core-SNP alignment.

        Returns an (n, n) matrix in which entry (i, j) is the fraction of
        sites at which the two sequences carry the same allele, ignoring
        positions that are missing in either. This is positive
        semi-definite by construction, as the LMM requires.
        """
        import numpy as np

        wanted = set(samples)
        seqs: Dict[str, List[str]] = {}
        name: Optional[str] = None
        chunks: List[str] = []
        with Path(alignment).open() as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if name in wanted:
                        seqs[name] = list("".join(chunks))
                    name = line[1:].split()[0]
                    chunks = []
                else:
                    chunks.append(line)
            if name in wanted:
                seqs[name] = list("".join(chunks))

        ordered = [s for s in samples if s in seqs]
        if len(ordered) < 2:
            raise ValueError(
                f"only {len(ordered)} of {len(samples)} samples found in {alignment.name}"
            )
        width = min(len(seqs[s]) for s in ordered)
        # Encode the ALLELE, not merely whether the site is called. A
        # called/not-called mask makes every called site compare equal, so
        # identity is 1.0 for every pair and the matrix carries no
        # information about relatedness at all.
        codes = {"A": 0, "C": 1, "G": 2, "T": 3}
        arr = np.array(
            [[codes.get(c.upper(), -1) for c in seqs[s][:width]] for s in ordered],
            dtype=np.int8,
        )
        called = arr >= 0
        n = len(ordered)
        out = np.zeros((n, n))
        for i in range(n):
            for j in range(i, n):
                both = called[i] & called[j]
                total = int(both.sum())
                value = float((arr[i][both] == arr[j][both]).sum()) / total if total else 0.0
                out[i, j] = out[j, i] = value
        return out

    def _cross_check_pyseer_tested_count(self, log: Path) -> None:
        """Compare the patterns file against pyseer's own tally, from its log.

        pyseer writes ``%d tested variants`` to stderr
        (``pyseer/__main__.py``, alongside the loaded/filtered/printed counts).
        That is pyseer's own statement of how many variants reached the model,
        and it is the size of the population the reduction runs over - so when
        it disagrees with the patterns file, one of the two is misreading the
        run and every count derived from the file is untrustworthy.

        To be precise about the arithmetic, because the two numbers are easy to
        confuse: the **divisor** in ``alpha / n_unique_patterns`` is
        :attr:`UniquePatternReduction.n_unique_patterns` - the number of
        *distinct* patterns, which is what makes it a Bonferroni denominator.
        ``n_variants_tested`` is the input to the reduction, and this check is
        what establishes that the input population is the one pyseer intended.
        This method validates the input, not the divisor.

        A log without the line is not a disagreement: a caller who captured
        pyseer's stderr elsewhere, or an older pyseer that did not emit it,
        leaves nothing to compare, and inventing a second opinion is not better
        than recording that there was none.

        Raises:
            DataContractError: pyseer and the patterns file disagree.
        """
        if not Path(log).is_file():
            LOGGER.debug("No pyseer log at %s; no independent count to check", log)
            return
        stated: Optional[int] = None
        for line in Path(log).read_text(errors="replace").splitlines():
            match = _PYSEER_TESTED_COUNT_RE.match(line.strip())
            if match:
                stated = int(match.group(1))
        if stated is None:
            LOGGER.info(
                "pyseer's log at %s does not report a tested-variant count, so "
                "the pattern count could not be cross-checked against pyseer's "
                "own; the reduction stands on the helper and the local tally.",
                log,
            )
            return
        if stated != self.reduction.n_variants_tested:
            raise DataContractError(
                f"pyseer reported testing {stated} variants but the patterns "
                f"file holds {self.reduction.n_variants_tested}. The threshold "
                "is derived by reducing the patterns in that file to their "
                "distinct values, so a file that does not hold the variants "
                "pyseer tested means the unique-pattern count - and therefore "
                "the threshold - rests on a population neither pyseer nor "
                "this file supports. Refused rather than reported.",
                pyseer_tested=stated,
                patterns_in_file=self.reduction.n_variants_tested,
                log=str(log),
            )
        LOGGER.info(
            "Pattern count agrees with pyseer's own tested-variant count (%d)",
            stated,
        )

    @staticmethod
    def _normalise_columns(raw: Path) -> Path:
        """Map pyseer's column names onto the parser's schema.

        Three shapes are distinguished, and the distinction is the point:

        * **no lines at all** - a truncated or unwritten file, not a result.
          pyseer prints its header before the association loop, so an exit-0 run
          never produces this. Left empty for the parser to refuse.
        * **a header and no rows** - pyseer ran and no variant passed
          ``--lrt-pvalue``. After step 12a that is the ordinary outcome of a
          null result rather than a corner case, and it becomes an empty
          result.
        * **a header and rows** - normalised and parsed.

        A single line that is *not* a header is none of these. It is a
          truncated write or a stray diagnostic, and it must not be read as an
          empty association table: that would turn a broken run into a reported
          negative result. :func:`parse_pyseer_output` refuses it.

        Every branch rewrites the target, so a previous run's
        ``pyseer_results.tsv`` can never be read as this one's.

        Returns:
            The path of the normalised table.
        """
        target = raw.parent / "pyseer_results.tsv"
        with raw.open() as handle:
            lines = [line for line in handle if not line.startswith("#")]
        if not lines:
            # Truncated or unwritten. Left as an empty file, which
            # `parse_pyseer_output` refuses - distinct from the header-only
            # case below, which is a real "nothing passed" result.
            LOGGER.error(
                "pyseer output at %s is empty, so it is not an association "
                "table; it will be refused rather than read as no results.",
                raw,
            )
            target.write_text("")
            return target
        header = lines[0].rstrip("\n").split("\t")
        if not _PYSEER_HEADER_REQUIRED <= set(header):
            # Not pyseer's table. Passed through unchanged, deliberately: the
            # parser's own `feature`-column check then refuses it by name, with
            # the columns it did find. Normalising it here would mean guessing
            # which line was meant to be a header, which is the guess that
            # turns a broken run into an empty result.
            # The names come from the gate's own set, not from a second
            # literal: a message that named columns the gate no longer required
            # would go stale silently, and a refusal whose explanation no longer
            # matches the refusal is worse than no explanation.
            LOGGER.error(
                "pyseer output at %s is missing the columns %s, so it is not "
                "an association table. It is passed through un-normalised and "
                "will be refused.",
                raw,
                "/".join(sorted(_PYSEER_HEADER_REQUIRED)),
            )
            shutil.copyfile(raw, target)
            return target
        if len(lines) == 1:
            target.write_text("feature\tpvalue\n")
            return target
        out_header = [PYSEER_COLUMN_RENAME.get(name, name) for name in header]
        with target.open("w") as handle:
            handle.write("\t".join(out_header) + "\n")
            for line in lines[1:]:
                if line.strip():
                    handle.write(line if line.endswith("\n") else line + "\n")
        return target


#: The correction this parser applies to ``adjusted_p_value`` when pyseer did
#: not supply one. Named once so the refusal, the docstring and the default all
#: quote the same string rather than three literals that can drift apart.
SUPPORTED_ADJUSTED_P_CORRECTION = "fdr_bh"


def parse_pyseer_output(
    path: Path,
    gwas_input: GwasInput,
    model: str = "pyseer:mixed",
    correction: str = SUPPORTED_ADJUSTED_P_CORRECTION,
) -> List[GwasResult]:
    """Parse a pyseer association table into the standard schema.

    Recognised columns: ``feature``, ``pvalue``, ``adj.pvalue``, ``OR``,
    ``beta``, ``freq``. Missing adjusted p-values are computed locally.

    **Which correction is in force, and why the two are not the same thing.**

    Two corrections apply, at two different places, and conflating them is the
    error this docstring exists to prevent:

    * **Bonferroni, upstream, in step 12a.** ``alpha / n_unique_patterns`` is
      handed to pyseer as ``--lrt-pvalue``, and pyseer drops every variant at or
      above it (``pyseer/lmm.py``: ``if p_values[i] >= lrt_pvalue``). This is
      the gate that decides what appears in the table at all, and
      ``docs/scientific_rules.md`` section 3a pins it to pyseer's own
      ``count_patterns.py``. It is *not* what fills ``adjusted_p_value``.
    * **Benjamini-Hochberg, here.** ``adjusted_p_value`` is the column the
      contract declares and the one :func:`significant` filters on, at
      ``config.gwas.significance_threshold``. ``docs/scientific_rules.md``
      section 3 lists ``benjamini_hochberg()`` as the enforcement for that
      column.

    The ``adj.pvalue`` branch is **unreachable against pyseer 1.1.2**, which is
    worth stating rather than leaving as a silent fallback: that package
    contains no FDR computation anywhere (verified against the installed
    ``site-packages/pyseer`` - no ``adj.pvalue``, no ``fdr``, no
    ``multipletests``), and its header list ends at ``notes``. Reading that
    column would therefore be reading something no supported pyseer emits. It
    is kept because it costs one ``row.get`` and it is the right thing to do if
    a future pyseer grows the column, but the BH branch below is the one every
    real run takes - which is also what makes the rank bug it used to contain
    live rather than dormant.

    Args:
        correction: The correction to apply to ``adjusted_p_value`` where
            pyseer did not supply one. Read from ``config.gwas.correction``.
            Any value this function does not implement is refused rather than
            ignored, because silently applying BH under a key that says
            something else is how a run reports a science nobody chose.

    Raises:
        DataContractError: The file is not an association table, or
            ``correction`` is not one this function computes.
    """
    if correction != SUPPORTED_ADJUSTED_P_CORRECTION:
        raise DataContractError(
            f"config.gwas.correction is {correction!r} and this parser computes "
            f"only {SUPPORTED_ADJUSTED_P_CORRECTION!r}. Accepting a name and "
            "then applying a different correction would let a run report a "
            "multiple-testing procedure nobody chose, so the key is refused "
            "rather than honoured with the wrong arithmetic. Note this is a "
            "different correction from step 12a's, which is Bonferroni on "
            "unique patterns and is what gates which variants pyseer reports.",
            correction=correction,
            supported=SUPPORTED_ADJUSTED_P_CORRECTION,
        )
    path = Path(path)
    if not path.exists():
        LOGGER.warning("pyseer output not found: %s", path)
        return []

    # Order matters here, because three different situations all present as
    # "one line" or fewer and only two of them are results.
    #
    # A header carrying no rows is a legitimate result: pyseer emits exactly
    # that when no variant passes --lrt-pvalue, which after step 12a is the
    # ordinary outcome of a cohort with no association. `read_tsv` refuses a
    # header-only file - correctly, for a table that is supposed to hold
    # records - so that case is recognised here instead.
    #
    # A single line that is *not* a header is a truncated write or a stray
    # diagnostic. Reading it as "no associations" would convert a broken run
    # into a reported negative result, so it is refused with the columns that
    # were actually found.
    lines = path.read_text().splitlines()
    if not lines:
        # Not a result. With the real binary, exit 0 always implies at least a
        # header (`pyseer/__main__.py` prints one before the association loop),
        # so a zero-byte table is a truncated or unwritten file - the same class
        # of thing as the single-non-header-line case below, and refused for the
        # same reason: reporting it as "no associations" would convert a broken
        # run into a negative scientific result.
        raise DataContractError(
            "pyseer produced an empty output file. With the real binary an "
            "exit-0 run always writes a header before the association loop, so "
            "this is a truncated or unwritten file rather than a run that found "
            "nothing - and a run that found nothing writes the header alone. "
            "Refused rather than read as 'no associations'.",
            path=str(path),
        )
    columns = set(lines[0].split("\t"))
    if not _PARSE_REQUIRED_COLUMNS <= columns:
        raise DataContractError(
            "pyseer output carries no 'feature' column, so it is not an "
            "association table this stage can read. The columns present are "
            f"{','.join(sorted(columns)) or '(none)'}. A truncated write or a "
            "diagnostic printed where the table belongs both look like this, "
            "and are not the same as a run that found nothing.",
            path=str(path),
            columns=",".join(sorted(columns)),
            required_columns=",".join(sorted(_PARSE_REQUIRED_COLUMNS)),
        )
    if not [line for line in lines[1:] if line.strip()]:
        LOGGER.info(
            "pyseer reported no feature passing the threshold step 12a derived; "
            "there is no association to report. The threshold and how it was "
            "derived are in %s",
            path.parent / UNIQUE_PATTERNS_SUMMARY_FILENAME,
        )
        return []

    rows = read_tsv(path, required_columns=tuple(sorted(_PARSE_REQUIRED_COLUMNS)))
    parsed: List[Tuple[str, Optional[float], Optional[float], Optional[float], Optional[float]]] = []
    for row in rows:
        feature = str(row["feature"])
        parsed.append(
            (
                feature,
                _to_float(row.get("pvalue") or row.get("p_value")),
                _to_float(row.get("adj.pvalue") or row.get("adjusted_p_value")),
                _to_float(row.get("OR")),
                _to_float(row.get("beta")),
            )
        )

    # Index -> p-value, for the rows that need an adjustment, so the adjusted
    # value can be looked up by POSITION rather than by p-value.
    #
    # The previous lookup was `computed.get(str(need_adjustment.index(p)))`,
    # and `.index()` returns the first match. Two features sharing a raw
    # p-value therefore received the *same* adjusted value, which is wrong even
    # when the two ties are adjacent: BH is monotone in rank, so two equal
    # p-values at different ranks have different adjusted values, and `index()`
    # hands both the earlier rank's. It is not a rare shape here - pyseer
    # prints `%.2E`, so six significant digits collapse into far fewer distinct
    # values than there are features, and a null cohort ties constantly. The
    # keys of `computed` are positions, so they are what the lookup uses.
    need_adjustment = [
        p for _, p, adj, _, _ in parsed if p is not None and adj is None
    ]
    computed = benjamini_hochberg({str(i): p for i, p in enumerate(need_adjustment)})
    #: Position in ``parsed`` -> its adjusted p-value, or ``None`` when the row
    #: carried one from pyseer or has no usable p-value.
    adjusted_by_position: List[Optional[float]] = []
    cursor = 0
    for _feature, p_value, adjusted, _or, _beta in parsed:
        if adjusted is not None or p_value is None:
            adjusted_by_position.append(adjusted)
            continue
        adjusted_by_position.append(computed.get(str(cursor)))
        cursor += 1

    results: List[GwasResult] = []
    for index, (feature, p_value, adjusted, odds_ratio, beta) in enumerate(parsed):
        vector = gwas_input.features.get(feature)
        frequency: Optional[float] = None
        distribution: Dict[str, int] = {}
        if vector is not None:
            carriers = sum(vector.get(s, 0) for s in gwas_input.sample_ids)
            frequency = carriers / gwas_input.n_samples if gwas_input.n_samples else None
            for sample_id in gwas_input.sample_ids:
                if vector.get(sample_id, 0):
                    lineage = gwas_input.lineages.get(sample_id, "unknown")
                    distribution[lineage] = distribution.get(lineage, 0) + 1
        adjusted_value = adjusted_by_position[index]
        feature_type, _ = parse_feature_name(feature)
        results.append(
            GwasResult(
                feature=feature,
                feature_type=gwas_input.feature_types.get(feature, feature_type),
                effect=odds_ratio if odds_ratio is not None else beta,
                p_value=p_value,
                adjusted_p_value=adjusted_value,
                effect_size=odds_ratio if odds_ratio is not None else beta,
                frequency=frequency,
                lineage_distribution=distribution,
                model=model,
            )
        )
    return results


def _to_float(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


# --------------------------------------------------------------------------
# Statistics
# --------------------------------------------------------------------------


def benjamini_hochberg(pvalues: Mapping[str, float]) -> Dict[str, float]:
    """Benjamini-Hochberg FDR adjustment.

    Implemented directly so that multiple-testing correction is a tested
    property of this codebase rather than a dependency of the engine.
    """
    items = [(k, v) for k, v in pvalues.items() if v is not None]
    if not items:
        return {}
    items.sort(key=lambda kv: kv[1])
    n = len(items)
    adjusted: Dict[str, float] = {}
    previous = 1.0
    for rank in range(n, 0, -1):
        key, p_value = items[rank - 1]
        value = min(previous, p_value * n / rank)
        adjusted[key] = min(1.0, value)
        previous = adjusted[key]
    return adjusted


def is_lineage_confounded(
    result: GwasResult, gwas_input: GwasInput, threshold: float
) -> bool:
    """Whether a feature is essentially confined to a single lineage.

    A feature carried only by samples from one lineage cannot be separated
    from that lineage by any association test, so it is reported as
    lineage-linked rather than as a resistance association.
    """
    if not result.lineage_distribution:
        return False
    total = sum(result.lineage_distribution.values())
    if total == 0:
        return False
    top = max(result.lineage_distribution.values())
    return (top / total) >= threshold


def flag_lineage_linked(
    results: Sequence[GwasResult], gwas_input: GwasInput, config: PipelineConfig
) -> List[GwasResult]:
    """Return the subset of results that are lineage-confounded."""
    threshold = config.gwas.lineage_confound_threshold
    return [r for r in results if is_lineage_confounded(r, gwas_input, threshold)]


def significant(
    results: Sequence[GwasResult], config: PipelineConfig
) -> List[GwasResult]:
    """Results passing the configured adjusted-p threshold."""
    threshold = config.gwas.significance_threshold
    return [
        r
        for r in results
        if r.adjusted_p_value is not None and r.adjusted_p_value <= threshold
    ]


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------


#: The stage-12 feature table: the one input this stage reads from disk, and
#: the one whose producer this stage does not own. Named as a constant so the
#: refusal that reports it missing, and any caller that has to provision it,
#: quote one string.
GWAS_FEATURE_TABLE_NAME = "gwas_features.tsv"

#: Which stage's contract declares that table, and where the declaration lives.
#: The refusal names both, because "file not found" without a producer sends the
#: reader looking for a bug in this stage rather than for the missing edge.
GWAS_FEATURE_TABLE_CONTRACT = "docs/data_contract.md, 'gwas/gwas_features.tsv'"


def load_features(path: Path) -> List[Dict[str, Optional[str]]]:
    """Load the stage 12 feature table.

    Raises:
        DataContractError: The table is absent. The message names the file, the
            stage that consumes it and the stage whose contract declares it,
            because an absent input and an unreadable one are different faults
            and only one of them is fixed by writing a file somewhere.
    """
    path = Path(path)
    if not path.exists():
        raise DataContractError(
            f"Stage 12 (gwas) cannot find its feature table {path.name!r} at "
            f"{path}. This is the one input this stage reads from disk: one row "
            "per manifest sample, `sample_id` plus one binary column per feature "
            "named `<type>__<label>`. It is an input, not an output - this stage "
            "reads it and never writes it, so no amount of re-running stage 12 "
            "will produce it. Stage 12 (gwas) is the consumer; the schema is "
            f"declared by {GWAS_FEATURE_TABLE_CONTRACT}. Provision it there, or "
            "wire the producer that writes it before this stage runs.",
            path=str(path),
            input_name=GWAS_FEATURE_TABLE_NAME,
            contract=GWAS_FEATURE_TABLE_CONTRACT,
            consumer_stage="gwas",
        )
    return read_tsv(path, required_columns=("sample_id",), unique_columns=("sample_id",))


def run(
    config: PipelineConfig,
    manifest: SampleManifest,
    mode: RunMode,
    intermediate_root: Path,
    phenotype_calls: Sequence[PhenotypeCall],
    lineages: Optional[Mapping[str, str]] = None,
    engine: Optional[GwasEngine] = None,
) -> Tuple[List[GwasResult], GwasInput]:
    """Stage 12 entry point.

    Returns:
        ``(results, constructed_input)``. The input is returned so that
        later stages and the report can state exactly which samples and
        features were tested.

    Raises:
        NotImplementedError: In ``REAL`` mode without a real engine. See the
            guard below.
        DataContractError: The feature table is absent, or one outcome group is
            below the configured minimum.
    """
    if mode is not RunMode.TEST and (
        engine is None or isinstance(engine, ReferenceEngine)
    ):
        # Placed first, before the feature table is read, so refusing REAL does
        # not itself depend on intermediates being present.
        #
        # The engine condition is part of THIS guard rather than a nested `if`,
        # and that is deliberate. The taxonomy self-verification test
        # (`test_it_matches_the_stage_modules_that_refuse`) classifies a raise
        # by its INNERMOST enclosing guard, to tell `similarity`'s raise inside
        # an `allow_real_mode` gate (a gate, not a refusal) from a raise that
        # refuses REAL itself. A nested `if isinstance(engine, ...)` would make
        # the innermost guard the engine check, which mentions no mode, and the
        # stage would be classified as NOT refusing REAL - while `run.py` lists
        # it in REAL_REFUSING_STAGES, correctly, because an engine-less REAL
        # run does raise. Folding the condition into the mode guard keeps the
        # classifier's rule working and keeps the two descriptions agreeing.
        #
        # The refusal is now CONDITIONAL. `PyseerEngine` is a complete class
        # (`--lmm`, kinship wiring, two passes at the step-12a threshold) and a
        # REAL run that injects one must be allowed to use it; refusing it would
        # be refusing a mixed model because a guard could not tell one from a
        # contingency table.
        #
        # What is still refused is the engine-less REAL run, and the reason has
        # changed since the original message. That message said "no caller yet",
        # which stopped being true when `PyseerEngine` landed, so it explained a
        # live block with a dead reason. The two blockers that remain are:
        #
        #   1. No real engine was injected, so the run would reach
        #      `ReferenceEngine` - Fisher's exact with BH and no
        #      population-structure correction, documented for pipeline
        #      mechanics only. `run.py` DOES pass `engine=`: its stage-12 call
        #      site passes `engine=_build_gwas_engine(config, resolved,
        #      intermediate)`. What that helper returns in REAL is `None`
        #      unless pyseer resolves from the machine config (PATH first,
        #      then the overlay's `tool_search_dirs`), and it is that `None` -
        #      or a caller that explicitly injects a `ReferenceEngine` - that
        #      reaches this guard. So the live condition is "pyseer is missing
        #      or the wrong engine was injected", not "the call site was never
        #      wired".
        #   2. Even with `PyseerEngine`, this stage does not correct for
        #      *lineage* as a model term. pyseer corrects for relatedness
        #      (a random effect via `--similarity`); it does not take lineage as
        #      a fixed effect here, because doing so needs pyseer's
        #      `--lineage`, which in turn requires a `--distances` matrix
        #      (`pyseer/__main__.py` exits without one when `--lmm --lineage`
        #      is passed, and `pyseer/input.py`'s `load_structure` indexes the
        #      file as `m.loc[ids, ids]`, so the header must name every
        #      sample). Stage 10 leaves no such file: its own writer emits the
        #      square matrix, but `run.py`'s dispatch then overwrites that same
        #      path with the contract form (`STAGE_TABLES["similarity"]` -
        #      `sample_id` plus one packed `distances` column), which pyseer
        #      cannot index. Lineage is therefore handled by *flagging* a
        #      lineage-confounded feature (`is_lineage_confounded`), not by
        #      removing the confounding from the model. A flagged feature is
        #      reported as lineage-linked rather than as resistance-associated;
        #      it is not corrected.
        #
        # Both halves are stated rather than one, because they call for
        # different fixes and a reader told only one of them will conclude the
        # wrong thing about what is missing.
        raise NotImplementedError(
            "REAL-mode GWAS refuses to run without a real engine. Two "
            "blockers, both real. (1) No engine was injected, so this run "
            "would reach ReferenceEngine - Fisher's exact with BH and no "
            "population-structure correction - which is documented for "
            "pipeline mechanics testing only; construct PyseerEngine and "
            "pass it as engine=. (2) Even then, this stage corrects for "
            "relatedness via pyseer's --lmm random effect but does not "
            "correct for lineage as a model term: a feature confined to one "
            "lineage is FLAGGED as lineage-linked rather than deconfounded, "
            "because doing so needs pyseer's --lineage and a --distances "
            "matrix stage 10 does not write in usable form. The TEST path "
            f"reads the committed fixtures. (mode={getattr(mode, 'value', mode)})"
        )

    feature_path = Path(intermediate_root) / "gwas" / GWAS_FEATURE_TABLE_NAME
    features = load_features(feature_path)

    gwas_input = build_input(
        config, manifest, phenotype_calls, features, lineages=lineages
    )

    if engine is None:
        # TEST-only, deliberately. The guard above has already refused a REAL
        # run that arrives here without an engine, so this fallback can never be
        # the one a real cohort is analysed with.
        engine = ReferenceEngine()

    results = engine.run(gwas_input, config)
    confounded = flag_lineage_linked(results, gwas_input, config)
    if confounded:
        LOGGER.warning(
            "%d of %d features are confined to a single lineage and are "
            "lineage-linked rather than resistance-associated: %s",
            len(confounded),
            len(results),
            ", ".join(r.feature for r in confounded[:8]),
        )

    hits = significant(results, config)
    LOGGER.info(
        "Stage 12: %d features tested | %d pass adj.p<=%s | %d lineage-confounded "
        "| engine=%s",
        len(results),
        len(hits),
        config.gwas.significance_threshold,
        len(confounded),
        engine.name,
    )
    return results, gwas_input
