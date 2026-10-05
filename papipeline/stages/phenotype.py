"""Stage 11 - phenotype input interface.

Responsibilities
----------------
* validate the phenotype table against the data contract;
* restrict to the configured antibiotic(s);
* normalise R/I/S/SDD/ND;
* carry quantitative fields **only** when they were actually measured.

Hard rules enforced here (see ``docs/scientific_rules.md``):

* R/I/S are never converted into MIC or zone diameter;
* a zone diameter is only read from a real measurement column, never
  derived from a category;
* a sample with no phenotype record stays ``None`` downstream rather than
  being defaulted to S or R;
* an unknown antibiotic is a hard error, not a silent drop;
* the primary contrast is R vs S (ruling R7); ``I`` and ``SDD`` are excluded
  and counted, never merged into either arm.

Where the AST provenance lives
------------------------------

``ast_method`` / ``ast_standard`` / ``ast_edition`` are read here and stored on
:class:`~papipeline.models.PhenotypeCall`, and ``PhenotypeCall.to_row()``
exports them. They are **not** written to ``11_phenotype.tsv``: that table's
column list is declared in ``papipeline/execution/contracts.py``
(``STAGE_TABLES['phenotype']``) and omits all three, and ``io.tsv.write_tsv``
passes ``extrasaction="ignore"`` to ``csv.DictWriter``, so they are dropped at
the write. In memory is therefore the only layer that holds them, and any
consumer reading the table back from disk sees none. Verified end to end by
``tests/unit/test_phenotype_ast_provenance_layers.py``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from ..config.loader import PipelineConfig
from ..errors import PhenotypeError, PipelineError
from ..io.tsv import read_tsv
from ..logging_utils import get_logger
from ..models import Phenotype, PhenotypeCall

LOGGER = get_logger("stages.phenotype")

REQUIRED_COLUMNS: Tuple[str, ...] = ("sample_id", "antibiotic", "phenotype")
QUANTITATIVE_COLUMNS: Tuple[str, ...] = (
    "MIC",
    "MIC_unit",
    "zone_diameter",
    "zone_unit",
)

PHENOTYPE_FILENAME_TEMPLATE = "{antibiotic}_phenotype.tsv"

#: What `config/science.yaml` writes into `breakpoints.edition` when no edition
#: has been supplied. Read as a value rather than re-declared, so that changing
#: the sentinel in the config changes what this stage accepts with it.
UNSUPPLIED_EDITION = "UNSUPPLIED"

#: What `BreakpointStatus.from_config` substitutes when `breakpoints.standard`
#: is absent, and the value that makes "no standard was declared" detectable.
UNSPECIFIED_STANDARD = "unspecified"

#: The values that mean "this field was never declared". Two, not one, because
#: they arrive by different routes: a config that omits `edition` gets
#: ``UNSPECIFIED_STANDARD`` from ``from_config``'s own default, while a config
#: that sets ``edition: "UNSUPPLIED"`` states the omission explicitly. Both are
#: an undeclared edition and both must be refused. Checking only the first
#: would let every config that simply left the key out reinterpret MICs.
UNDECLARED_VALUES = frozenset({UNSUPPLIED_EDITION, UNSPECIFIED_STANDARD})

#: Every susceptibility category ``Phenotype`` defines, and therefore every
#: code this stage can recognise as a real AST result.
#:
#: A code outside this set is a different thing from a code the model excludes.
#: ``I`` is a genuine category that this model refuses to compare; ``RESISTANT``
#: is a spelling of nothing the pipeline knows, and deciding it meant ``R``
#: would be a guess about a clinical record. The two refusals are therefore
#: separate, and both name the offending code.
AST_CATEGORIES = frozenset(p.value for p in Phenotype)

#: The PRIMARY phenotype contrast (ruling R7): resistant versus susceptible.
#: These two, and only these two, enter the association model.
PRIMARY_POSITIVE = "R"
PRIMARY_NEGATIVE = "S"

#: Categories the primary contrast EXCLUDES, and which are therefore COUNTED
#: rather than used. ``I`` (susceptible with increased exposure) and ``SDD``
#: (susceptible, standard-dose dosing) are both real, both measured, and both
#: outside R-vs-S. Folding either into R or S is the failure R7 exists to
#: prevent: the intermediate would be reported as a resistance call the
#: laboratory never made, and the SDD would be reported as a plain susceptible.
EXCLUDED_CATEGORIES: Tuple[str, ...] = ("I", "SDD", "ND")


class PhenotypePreflightError(PipelineError):
    """An input stage 11 needs is missing. Named so the cause is unambiguous.

    Separate from :class:`PhenotypeError` on purpose. A ``PhenotypeError`` means
    a row is wrong and the fix is to correct that row; a
    ``PhenotypePreflightError`` means the inputs do not yet support what is
    being asked of this stage, and the fix is to provision something. Two
    failures, two different remedies, and a caller that has to tell them apart
    cannot do it if both arrive as the same type.
    """


def _optional_text(value: Optional[str]) -> Optional[str]:
    """A sentinel or blank provenance field is absent, not a literal string."""
    if value is None:
        return None
    text = str(value).strip()
    if not text or text in {".", "-", "NA", "na", "none", "None"}:
        return None
    return text


def _to_float(value: Optional[str], column: str, sample_id: str) -> Optional[float]:
    """Parse an optional quantitative field. A present-but-unparsable value
    is an error; an absent value stays ``None``."""
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        raise PhenotypeError(
            "Quantitative phenotype field is not a number",
            sample_id=sample_id,
            column=column,
            value=str(value),
        ) from None


def parse_phenotype_row(
    row: Mapping[str, Optional[str]],
    allowed: Sequence[Phenotype],
    antibiotics: Sequence[str],
    require_antibiotic_config: bool = True,
    *,
    log_base: int = 2,
) -> PhenotypeCall:
    """Validate and convert one phenotype row.

    Raises:
        PhenotypeError: Unknown category, unknown antibiotic, or a
            quantitative field that cannot be parsed.
    """
    raw_id = row.get("sample_id")
    if not raw_id:
        raise PhenotypeError("Phenotype row has no sample_id")
    sample_id = str(raw_id)

    raw_antibiotic = row.get("antibiotic")
    if not raw_antibiotic:
        raise PhenotypeError(
            "Phenotype row has no antibiotic", sample_id=sample_id
        )
    antibiotic = str(raw_antibiotic)
    if require_antibiotic_config and antibiotic not in antibiotics:
        raise PhenotypeError(
            "Phenotype row references an antibiotic that is not configured",
            sample_id=sample_id,
            antibiotic=antibiotic,
            configured=",".join(antibiotics),
        )

    raw_phenotype = row.get("phenotype")
    if raw_phenotype is None:
        raise PhenotypeError(
            "Phenotype row has no phenotype value", sample_id=sample_id
        )
    token = str(raw_phenotype).strip().upper()
    if token not in {p.value for p in allowed}:
        # Two refusals, because a code outside the model's contrast and a code
        # that is not an AST category at all are different mistakes with
        # different fixes. Both name the offending code in the message text as
        # well as in the context, so a reader sees it without reading kwargs.
        allowed_values = ",".join(p.value for p in allowed)
        if token not in AST_CATEGORIES:
            raise PhenotypeError(
                f"Phenotype value {token!r} is not in the allowed set "
                f"({allowed_values}) and is not a recognised AST category "
                f"either. This pipeline reads one of: "
                f"{','.join(sorted(AST_CATEGORIES))}. {token!r} is not a "
                "spelling of any of them, and mapping it onto R or S would be "
                "a guess about a clinical record rather than a reading of one. "
                "Correct the source row, or add the category to "
                "`phenotype.allowed_values` in config/science.yaml if the "
                "source genuinely uses it.",
                sample_id=sample_id,
                antibiotic=antibiotic,
                value=token,
                allowed=allowed_values,
                recognised_categories=",".join(sorted(AST_CATEGORIES)),
            )
        raise PhenotypeError(
            f"Phenotype value {token!r} is a real AST category but is not in "
            f"the allowed set for this run: {allowed_values}. It is refused, "
            "not folded: the primary contrast is "
            f"{PRIMARY_POSITIVE}-vs-{PRIMARY_NEGATIVE} (ruling R7), and "
            f"{', '.join(EXCLUDED_CATEGORIES)} are excluded and counted "
            "instead of being reassigned to either side.",
            sample_id=sample_id,
            antibiotic=antibiotic,
            value=token,
            allowed=allowed_values,
            excluded_categories=",".join(EXCLUDED_CATEGORIES),
        )
    phenotype = Phenotype(token)

    mic = _to_float(row.get("MIC"), "MIC", sample_id)
    mic_unit = row.get("MIC_unit")
    zone = _to_float(row.get("zone_diameter"), "zone_diameter", sample_id)
    zone_unit = row.get("zone_unit")

    if mic is not None and phenotype in (Phenotype.ND,):
        raise PhenotypeError(
            "An MIC is reported for a sample with no determinate category",
            sample_id=sample_id,
            phenotype=phenotype.value,
            hint="ND means not determined; it cannot carry an MIC",
        )

    return PhenotypeCall(
        sample_id=sample_id,
        antibiotic=antibiotic,
        log2_mic=interpret_mic(mic, base=log_base),
        ast_method=_optional_text(row.get("ast_method")),
        ast_standard=_optional_text(row.get("ast_standard")),
        ast_edition=_optional_text(row.get("ast_edition")),
        phenotype=phenotype,
        mic=mic,
        mic_unit=mic_unit,
        zone_diameter=zone,
        zone_unit=zone_unit,
        source=str(row.get("source") or "unknown"),
    )


def phenotype_path(
    phenotype_dir: Path, antibiotic: str, filename: Optional[str] = None
) -> Path:
    """Resolve the phenotype file for one antibiotic."""
    name = filename or PHENOTYPE_FILENAME_TEMPLATE.format(antibiotic=antibiotic)
    return Path(phenotype_dir) / name


def require_ast_table(path: Path, *, antibiotic: str) -> Path:
    """Refuse to proceed without a readable AST table, naming what builds one.

    The failure this prevents is a quiet one: without a table the stage has
    nothing to read, and a caller that treated "no records" as a result would
    report a cohort with no susceptibility data as though the cohort had been
    tested and found uniformly uninformative. A missing AST table is a missing
    input.

    Deliberately does not check the *contents*. Whether a cohort's calls may be
    pooled is a separate question, answered by
    :func:`require_reinterpretation_is_sourced` for the re-interpretation case
    and by ``provenance_report`` for comparability. Deciding here what a row is
    allowed to omit would mean inventing a rule that
    ``docs/scientific_rules.md`` does not state.
    """
    table = Path(path)
    if not table.is_file():
        raise PhenotypePreflightError(
            f"no {antibiotic} AST table at {table}. Stage 11 reads "
            "susceptibility calls; it does not obtain them, so there is no "
            "fallback and nothing to compute while you wait. For a bounded run "
            "the table is produced by `scripts/smoke/build_smoke_phenotype.py` "
            "from the `AST phenotypes` column of `PDC_essential.tsv`, written "
            "to the overlay's `paths.smoke_phenotype_dir`. For the full cohort "
            "it is `data/phenotype/<antibiotic>_phenotype.tsv` and must be "
            "supplied. This is a missing input, not an empty result.",
            path=str(table),
            antibiotic=antibiotic,
        )
    if table.stat().st_size == 0:
        raise PhenotypePreflightError(
            f"the {antibiotic} AST table at {table} is empty. An empty file is "
            "not a cohort with no resistance: it is a table that was created "
            "and never written, and reading it would silently produce a "
            "zero-sample result.",
            path=str(table),
            antibiotic=antibiotic,
        )
    return table


def require_reinterpretation_is_sourced(
    status: "BreakpointStatus",
) -> "BreakpointStatus":
    """Refuse to reinterpret an MIC against thresholds with no stated source.

    Two separate refusals, because they are two separate omissions and the
    remedy for each is different:

    1. **no standard.** A threshold set is meaningless without knowing which
       standard it came from, and for imipenem the choice is not cosmetic:
       CLSI and EUCAST do not place the S/I/R boundaries at the same MIC for
       *Pseudomonas aeruginosa*. Reinterpreting against the wrong one does not
       produce a slightly different answer, it produces a differently
       *categorised* isolate.
    2. **no edition.** Breakpoint tables are revised between editions. S/I/R
       values carrying no year cannot be checked against the edition the source
       laboratory used, so the result cannot be reproduced or audited.

    Neither refusal picks a standard. The repository names CLSI
    (``config/antibiotics.tsv`` -> ``phenotype_source_standard`` =
    ``CLSI_M100_table_2B``, and ``config/science.yaml`` -> ``breakpoints.
    standard: "CLSI M100"``), but which standard *this cohort's* laboratories
    actually applied is a fact about the source data that no file in this
    repository records. ``docs/scientific_rules.md`` is silent on whether a
    categorical call with no declared standard may still be reported, so this
    function refuses only the step that would *change* a call - re-interpretation
    - and leaves the reported source calls exactly as the laboratory made them.

    Returns ``status`` unchanged when no threshold is configured, which is the
    shipped state: with ``thresholds: {}`` nothing is reinterpreted and there
    is nothing to be wrong about.
    """
    if not status.is_configured:
        return status

    if not status.standard_declared:
        raise PhenotypePreflightError(
            "S/I/R thresholds are configured but no susceptibility standard is "
            "declared, so there is nothing to say which standard's breakpoints "
            "they are. Declaring `breakpoints.standard` in config/science.yaml "
            "is the fix; guessing one here is not, because CLSI and EUCAST "
            "disagree on the imipenem S/I/R boundaries for Pseudomonas "
            "aeruginosa and picking the wrong one re-categorises isolates "
            "rather than shifting a number. Note that config/antibiotics.tsv "
            "records the *intended* standard for this antibiotic, which is not "
            "the same thing as a declaration of what the source laboratory "
            "used.",
            thresholds=",".join(sorted(status.thresholds)),
            hint="Set breakpoints.standard in config/science.yaml",
        )

    if not status.edition_declared:
        raise PhenotypePreflightError(
            f"S/I/R thresholds are configured and attributed to "
            f"{status.standard}, but no breakpoint edition or year is declared "
            f"(breakpoints.edition is {status.edition!r}). Breakpoint tables "
            "are revised between editions, so unversioned thresholds cannot be "
            "checked against the edition the source laboratory used and the "
            "result cannot be reproduced. Set `breakpoints.edition` to the "
            "edition and year the S/I/R values were taken from.",
            standard=status.standard,
            edition=status.edition,
            hint="Set breakpoints.edition in config/science.yaml",
        )
    return status


def load_phenotype(
    config: PipelineConfig,
    phenotype_dir: Path,
    antibiotic: str,
    filename: Optional[str] = None,
    sample_ids: Optional[Sequence[str]] = None,
) -> List[PhenotypeCall]:
    """Load and validate phenotype records for one antibiotic.

    Args:
        config: Loaded pipeline configuration.
        phenotype_dir: Directory holding the phenotype TSV.
        antibiotic: Antibiotic to load. Must be configured.
        filename: Override for the file name.
        sample_ids: Optional manifest to check coverage against. Samples
            with no record are simply absent from the result; they are never
            given a default phenotype.

    Returns:
        Validated calls, ordered as in the file.

    Raises:
        PhenotypePreflightError: The AST table is missing or empty, or S/I/R
            thresholds are configured without a declared standard or edition.
        PhenotypeError: A row is malformed.
    """
    config.require_antibiotic(antibiotic)
    status = BreakpointStatus.from_config(config)
    # Before the read, not after: a table that cannot be interpreted must not
    # be interpreted, and finding that out after parsing every row would make
    # the same mistake it is refusing to make, just later.
    require_reinterpretation_is_sourced(status)
    log_base = status.log_base
    path = require_ast_table(
        phenotype_path(Path(phenotype_dir), antibiotic, filename),
        antibiotic=antibiotic,
    )
    rows = read_tsv(
        path,
        required_columns=REQUIRED_COLUMNS,
        unique_together=[("sample_id", "antibiotic")],
    )

    calls: List[PhenotypeCall] = []
    for row in rows:
        if row.get("antibiotic") != antibiotic:
            # A file may legitimately hold several antibiotics; skip the rest
            # rather than failing, but only for a configured antibiotic.
            other = str(row.get("antibiotic"))
            if other in config.antibiotics:
                continue
            raise PhenotypeError(
                "Phenotype file contains an unconfigured antibiotic",
                antibiotic=other,
                path=str(path),
            )
        calls.append(
            parse_phenotype_row(
                row, config.allowed_phenotypes, config.antibiotics, log_base=log_base
            )
        )

    if sample_ids is not None:
        present = {c.sample_id for c in calls}
        missing = [s for s in sample_ids if s not in present]
        if missing:
            LOGGER.warning(
                "%d of %d samples have no %s phenotype record: %s",
                len(missing),
                len(sample_ids),
                antibiotic,
                ",".join(missing[:10]) + ("..." if len(missing) > 10 else ""),
            )

    counts: Dict[str, int] = {}
    for call in calls:
        counts[call.phenotype.value] = counts.get(call.phenotype.value, 0) + 1
    LOGGER.info(
        "Loaded %d %s phenotype records: %s",
        len(calls),
        antibiotic,
        ", ".join(f"{k}={v}" for k, v in sorted(counts.items())),
    )
    # R7, in the log as well as in the report: the primary contrast is R vs S,
    # and anything held out of it has to be visible as a held-out count here.
    # A run log that lists I=6 beside R=7 and S=3 without saying so reads as
    # thirteen usable samples.
    report = provenance_report(calls)
    if report.excluded:
        LOGGER.info(
            "Primary contrast %s-vs-%s uses %d of %d %s records; %d excluded "
            "and counted, not merged into either arm: %s",
            PRIMARY_POSITIVE,
            PRIMARY_NEGATIVE,
            report.analysed,
            report.total,
            antibiotic,
            report.excluded,
            ", ".join(
                f"{k}={v}" for k, v in report.excluded_counts.items()
            ),
        )
    return calls


def phenotype_category_matrix(
    calls: Sequence[PhenotypeCall],
) -> Dict[str, Dict[str, str]]:
    """Long-format category matrix: sample -> antibiotic -> R/I/S/SDD/ND.

    This is the *reported* view. It is kept separate from
    :func:`phenotype_matrix`, which is the continuous trait the association
    model consumes, because the two answer different questions and conflating
    them is how a categorical label ends up standing in for a measurement.
    """
    matrix: Dict[str, Dict[str, str]] = {}
    for call in calls:
        matrix.setdefault(call.sample_id, {})[call.antibiotic] = call.phenotype.value
    return matrix


def phenotype_matrix(
    calls: Sequence[PhenotypeCall], *, log_base: int = 2
) -> Dict[str, float]:
    """The continuous trait pyseer consumes: sample_id -> log_base(MIC).

    A row with no measured MIC is **absent**, not zero. A zero would be a
    fabricated measurement sitting in the middle of the trait's range, and it
    would pull the association towards it. Absence is honest; zero is not.
    """
    return continuous_matrix(calls, log_base=log_base)



def binarise(
    call: PhenotypeCall, positive: Sequence[str], negative: Sequence[str]
) -> Optional[int]:
    """Map a category to 1/0 for GWAS, or ``None`` to exclude.

    A category that appears in neither ``positive`` nor ``negative`` is
    excluded rather than guessed. This is what keeps ``I`` from silently
    becoming either R or S.
    """
    value = call.phenotype.value
    if value in positive:
        return 1
    if value in negative:
        return 0
    return None


# --------------------------------------------------------------------------
# Continuous trait and breakpoint interpretation
# --------------------------------------------------------------------------


def interpret_mic(
    mic: Optional[float], base: int = 2
) -> Optional[float]:
    """The continuous trait pyseer consumes: log_base(MIC).

    Computed numerically, never by bucketing into a breakpoint. A MIC of
    0.2 mg/L is a measurement, and rounding it to "susceptible" discards the
    information the continuous model exists to use.

    Returns ``None`` for a value with no logarithm: zero, a negative number, or
    no measurement at all. A missing trait must stay missing - a zero here
    would be a fabricated measurement that quietly pulls the model.
    """
    if mic is None:
        return None
    value = float(mic)
    if value <= 0:
        return None
    if base <= 1:
        raise ValueError(f"log base must be greater than 1; got {base!r}")
    return math.log(value, base)


@dataclass(frozen=True)
class BreakpointStatus:
    """Whether a standard is actually available to interpret MICs against.

    The standard and edition are recorded *whether or not thresholds exist*,
    so a result always states which standard it was interpreted against -
    including the honest case where nothing was applied.
    """

    standard: str
    edition: str
    organism: str
    agent: str
    log_base: int
    thresholds: Mapping[str, float]
    exceptions: Mapping[str, str]
    source: str

    @property
    def is_configured(self) -> bool:
        """True only when a usable S/I/R triple is present."""
        return {"S", "I", "R"}.issubset(self.thresholds)

    @property
    def standard_declared(self) -> bool:
        """Whether a susceptibility standard was actually named.

        False both when the key is absent and when it holds the
        ``unspecified`` default, because those are the same fact: nobody said
        which standard.
        """
        return self.standard not in UNDECLARED_VALUES

    @property
    def edition_declared(self) -> bool:
        """Whether a breakpoint edition or year was actually named."""
        return self.edition not in UNDECLARED_VALUES

    @property
    def reason(self) -> str:
        if self.is_configured:
            return (
                f"MICs are interpreted against {self.standard} "
                f"{self.edition} (S/I/R thresholds configured)."
            )
        return (
            f"No S/I/R thresholds are configured, so no MIC was reinterpreted "
            f"into a category. Source calls stand unchanged. Intended standard: "
            f"{self.standard} ({self.edition}). CLSI M100 is a licensed "
            f"document; supply the {self.agent} S/I/R values to enable "
            f"re-interpretation."
        )

    @classmethod
    def from_config(cls, config, *, source: str = "science.yaml") -> "BreakpointStatus":
        section = dict((config.raw or {}).get("breakpoints") or {})
        thresholds = {
            str(k).upper(): float(v)
            for k, v in (section.get("thresholds") or {}).items()
            if v not in (None, "")
        }
        return cls(
            standard=str(section.get("standard") or "unspecified"),
            edition=str(section.get("edition") or "unspecified"),
            organism=str(section.get("organism") or "unspecified"),
            agent=str(section.get("agent") or ""),
            log_base=int(section.get("log_base") or 2),
            thresholds=thresholds,
            exceptions={
                str(k): str(v) for k, v in (section.get("exceptions") or {}).items()
            },
            source=source,
        )

    def category_for(self, mic: Optional[float]) -> Optional[str]:
        """S/I/R for a measured MIC, or ``None`` if no standard applies.

        ``None`` here means "not interpretable", which is different from any
        susceptibility category and is never collapsed into one.
        """
        if not self.is_configured or mic is None or float(mic) <= 0:
            return None
        value = float(mic)
        if value <= self.thresholds["S"]:
            return "S"
        if value <= self.thresholds["I"]:
            return "I"
        return "R"


@dataclass(frozen=True)
class ProvenanceReport:
    """How comparable this cohort's susceptibility calls actually are.

    Two things are recorded, and they are not the same question:

    * **category accounting** (``counts`` and the properties below) - how many
      rows carry each AST category. R and S are the primary contrast; ``I``,
      ``SDD`` and ``ND`` are excluded from it and counted (ruling R7). A cohort
      with 20 intermediate calls and a cohort with none are different cohorts,
      and only the count says which.
    * **provenance gaps** - how many rows lack a method, a standard or an
      edition. A gap here is not a defect in the row; it is a fact about the
      source, and it decides whether the calls may be pooled at all.
    """

    total: int
    with_mic: int
    lacking_method: int
    lacking_standard: int
    lacking_edition: int
    #: category value -> number of rows carrying it. Every category present is
    #: a key, including a zero-free absence: a category with no rows is simply
    #: not a key.
    counts: Mapping[str, int] = field(default_factory=dict)

    # -- R7: the primary contrast, and what is held out of it ---------------

    @property
    def positive(self) -> int:
        """Rows in the primary (resistant) arm."""
        return int(self.counts.get(PRIMARY_POSITIVE, 0))

    @property
    def negative(self) -> int:
        """Rows in the secondary (susceptible) arm."""
        return int(self.counts.get(PRIMARY_NEGATIVE, 0))

    @property
    def analysed(self) -> int:
        """Rows that actually enter the R-vs-S contrast."""
        return self.positive + self.negative

    @property
    def excluded_counts(self) -> Dict[str, int]:
        """Held-out categories and their counts, in a fixed order.

        Ordered by :data:`EXCLUDED_CATEGORIES`, not alphabetically, so the
        report reads I, SDD, ND - the order a reader expects and the order
        the reasoning goes in. Categories with no rows are omitted rather than
        shown as zero, because "no intermediate calls in this cohort" and "the
        intermediate count is zero" are the same fact and only one of them
        needs saying.
        """
        return {
            category: int(self.counts[category])
            for category in EXCLUDED_CATEGORIES
            if self.counts.get(category)
        }

    @property
    def excluded(self) -> int:
        """How many rows the primary contrast refused.

        Reported, never absorbed. A cohort whose excluded rows were folded
        into R or S would produce a larger apparent sample and a different
        answer, which is why this number is carried rather than recomputed
        downstream where nothing would notice it.
        """
        return sum(self.excluded_counts.values())

    @property
    def standard_recorded_everywhere(self) -> bool:
        """True when no row is missing its AST standard.

        False on an empty cohort as well as on a populated one, because
        "every row has a standard" is not a claim an empty cohort supports.
        """
        return self.total > 0 and self.lacking_standard == 0

    def summary(self) -> str:
        held = self.excluded_counts
        held_text = (
            ", ".join(f"{k}={v}" for k, v in held.items()) if held else "none"
        )
        return (
            f"{self.total} phenotype rows: {self.analysed} enter the "
            f"{PRIMARY_POSITIVE}-vs-{PRIMARY_NEGATIVE} contrast "
            f"({PRIMARY_POSITIVE}={self.positive}, {PRIMARY_NEGATIVE}="
            f"{self.negative}); {self.excluded} excluded and counted, not "
            f"merged into either arm ({held_text}). "
            f"{self.with_mic} carry a measured MIC. "
            f"AST method missing on {self.lacking_method}, standard missing on "
            f"{self.lacking_standard}, edition missing on {self.lacking_edition}. "
            f"Calls without comparable provenance may not be pooled confidently."
        )


def provenance_report(calls: Sequence[PhenotypeCall]) -> ProvenanceReport:
    """Count what the source did and did not tell us, per row.

    A missing column is treated as missing on every row, which is the honest
    reading: the file did not carry it.
    """
    counts: Dict[str, int] = {}
    for call in calls:
        value = call.phenotype.value
        counts[value] = counts.get(value, 0) + 1
    return ProvenanceReport(
        total=len(calls),
        with_mic=sum(1 for c in calls if c.mic is not None),
        lacking_method=sum(1 for c in calls if not c.ast_method),
        lacking_standard=sum(1 for c in calls if not c.ast_standard),
        lacking_edition=sum(1 for c in calls if not c.ast_edition),
        counts=counts,
    )


def continuous_matrix(
    calls: Sequence[PhenotypeCall], *, log_base: int = 2
) -> Dict[str, float]:
    """sample_id -> log_base(MIC), for the rows that have a measured MIC.

    Rows without a measurement are absent rather than zero.
    """
    matrix: Dict[str, float] = {}
    for call in calls:
        value = interpret_mic(call.mic, base=log_base)
        if value is None:
            continue
        matrix[call.sample_id] = value
    return matrix
