"""Stage 16 (part 2) - report generation.

Produces a Markdown (and optionally HTML) report describing what the pipeline
ran and what it produced.

**Every** non-real mode is explicitly labelled as not-real, in the title and in
the banner. TEST's numbers come from synthetic fixtures in `test_data/`, so the
report says so and leads with the pipeline's own health checks (unpinned
references, QC flags, sample counts) rather than with any finding - a reader
should not be able to mistake a smoke test for an analysis. STUB's outputs are
fabricated outright, from no input at all, and are labelled with the same
`STUB_BANNER` the workflow's own stub reporter uses.

The labelling is a lookup table keyed on every `RunMode` member rather than an
`if mode is TEST / else`, because that `else` is what once handed STUB the real
banner and the real title (`_BANNERS`). Both labels fail closed.
"""

from __future__ import annotations

import html as html_module
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from .. import __version__
from ..config.loader import PipelineConfig
from ..logging_utils import get_logger
from ..models import ClaimStatus, ConvergenceCategory, RunMode
from ..reporting.smoke import require_smoke_marker, smoke_marker
from ..stub import STUB_BANNER
from .phenotype import (
    EXCLUDED_CATEGORIES,
    PRIMARY_NEGATIVE,
    PRIMARY_POSITIVE,
    BreakpointStatus,
    ProvenanceReport,
    provenance_report,
)
from .report_tables import NOT_PRODUCED, RunTables, TableRead, load_run_tables

LOGGER = get_logger("stages.reporting")

TEST_BANNER = (
    "SYNTHETIC TEST DATA - NOT BIOLOGICAL RESULTS. Every number in this "
    "report was produced from generated fixtures in `test_data/`. Nothing "
    "here describes a real Pseudomonas aeruginosa genome, and no biological "
    "conclusion may be drawn from it."
)

REAL_BANNER = (
    "REAL DATASET REPORT. Interpret every association according to "
    "docs/scientific_rules.md: detection is not resistance, and association "
    "is not causation."
)

#: Provenance label per mode, keyed on **every** member of `RunMode`.
#:
#: This was `TEST_BANNER if mode is RunMode.TEST else REAL_BANNER`, and
#: `RunMode` has three members, so STUB took the REAL branch: a STUB report came
#: out titled `Pseudomonas aeruginosa Imipenem AMR Report` and bannered
#: `REAL DATASET REPORT`.
#:
#: STUB is the mode where that mislabelling is least recoverable. TEST's numbers
#: come from committed fixtures, so they are stable - wrong, but reproducibly
#: wrong, and a reader can eventually notice. STUB's are derived from nothing at
#: all: no genome, no tool, no parser, no fixture (`papipeline/stub.py`). There
#: is no input to re-run and diff against, so a STUB number can never afterwards
#: be shown to be fabricated. The banner is the only thing marking it.
#:
#: STUB's wording is imported from `papipeline.stub` rather than restated, so
#: that both producers of this artefact - `stub.fabricate_report`, which the
#: Snakemake rule calls, and `write_report` here - say the same thing. They were
#: two independent paths to one file, and only one of them was honest.
_BANNERS: Mapping[RunMode, str] = {
    RunMode.STUB: STUB_BANNER,
    RunMode.TEST: TEST_BANNER,
    RunMode.REAL: REAL_BANNER,
}

#: Report title per mode. Same reasoning as `_BANNERS`: a two-way
#: `is RunMode.TEST` test gave STUB the real title, which is the first line of
#: the document and the stem of the filename a reader is sent.
_TITLES: Mapping[RunMode, str] = {
    RunMode.STUB: "Pseudomonas aeruginosa Imipenem AMR - STUB Run Report",
    RunMode.TEST: "Pseudomonas aeruginosa Imipenem AMR - TEST Pipeline Report",
    RunMode.REAL: "Pseudomonas aeruginosa Imipenem AMR Report",
}


def banner_for(mode: RunMode) -> str:
    """The provenance banner for ``mode``.

    Falls back to :data:`~papipeline.stub.STUB_BANNER` for a mode the table does
    not know, so the miss under-claims rather than over-claims. Under-claiming
    provenance is recoverable; a fabricated report bannered REAL is the defect
    this whole table exists to prevent, and it must not be the default.
    """
    return _BANNERS.get(mode, STUB_BANNER)  # type: ignore[arg-type]


def title_for(mode: RunMode) -> str:
    """The report title for ``mode``, failing closed as :func:`banner_for` does."""
    return _TITLES.get(mode, _TITLES[RunMode.STUB])  # type: ignore[arg-type]


# --------------------------------------------------------------------------
# R16: the power flag
# --------------------------------------------------------------------------

#: Ruling R16's flag for a statistic computed on the ten-isolate cohort.
#:
#: Declared as a constant because it is an assertion a reader may rely on and
#: because it is the literal the ruling names. The template below generalises it
#: to whatever ``n`` a section actually used, so a report cannot claim a cohort
#: size it did not compute on - and at ``n == 10`` it renders exactly this
#: string, which is asserted in the tests.
R16_FLAG = "n=10, underpowered, not a finding"

#: The flag, parameterised by the cohort size the statistic was computed on.
POWER_FLAG_TEMPLATE = "n={n}, underpowered, not a finding"

#: What the flag means, stated once so no section paraphrases it.
POWER_FLAG_MEANING = (
    "Ruling R16. This is a description of the isolates in this run and nothing "
    "more. It is not evidence about *Pseudomonas aeruginosa* imipenem "
    "susceptibility in general, and no association, rate or difference below may "
    "be cited as a finding."
)


def power_flag(n: int) -> str:
    """The R16 flag for a statistic computed on ``n`` isolates."""
    return POWER_FLAG_TEMPLATE.format(n=n)


def power_flag_line(n: int) -> str:
    """The block-quoted flag a statistic-bearing section carries."""
    return f"> **{power_flag(n)}** - {POWER_FLAG_MEANING}"


@dataclass
class ReportContext:
    """Everything the report needs. Assembled by the orchestrator."""

    mode: RunMode
    config: PipelineConfig
    generated_at: str
    antibiotic: str
    n_samples: int
    sections: List[Tuple[str, str]] = field(default_factory=list)
    tables: List[Tuple[str, List[str], List[Dict[str, Any]]]] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    provenance: List[Dict[str, Any]] = field(default_factory=list)
    #: The finding-bearing aggregates this report was built from, by input name.
    #:
    #: Reporting reads nothing from disk - the orchestrator assembles the
    #: context - so unlike `convergence` and `cooccurrence` it cannot go looking
    #: for its inputs and notice they are absent. It has to be *told* what it
    #: was given, and that is what this field is. Empty means the caller
    #: declared nothing, which in REAL is an absence, not a waiver.
    declared_inputs: Mapping[str, Any] = field(default_factory=dict)
    #: Stage 11's calls, when the run has them.
    #:
    #: The provenance section is derived from these, and this is the field that
    #: makes it reachable: `phenotype_report` below is what calls
    #: `provenance_report`, so the only way the report states that the testing
    #: standard was not recorded is that the caller handed the calls over.
    #: ``None`` means "no phenotype was loaded", which the report states as an
    #: absence rather than passing over in silence - a section that appears
    #: only when there is something to say must say so when there is not.
    phenotype_calls: Optional[Sequence[Any]] = None
    #: The run's own output tables, read back from disk.
    #:
    #: ``None`` means the caller did not read them, and the report then says so
    #: in every section rather than falling back to the in-memory aggregates.
    #: That fallback is the defect this field exists to end: a report built from
    #: memory cannot tell "found nothing" from "never ran", and the round-12
    #: seam analysis found 24 output artefacts with no reader at all
    #: (`pa-artifacts/round12/seam-matrix.md` §e).
    run_tables: Optional[RunTables] = None
    #: Per-stage status, stage name -> state, as the orchestrator recorded it.
    #:
    #: Preferred over the disk-derived fallback in the per-stage status table,
    #: because only the orchestrator knows why a stage was refused or skipped.
    #: When empty the table falls back to whether each stage's table is on disk,
    #: which is a weaker claim and is labelled as one.
    stage_status: Mapping[str, str] = field(default_factory=dict)
    #: Stage name -> why it was refused or skipped, for the stages that have one.
    stage_refusals: Mapping[str, str] = field(default_factory=dict)
    #: Per-isolate oprD locus verdicts, when the run has them.
    #:
    #: ``None`` and ``[]`` mean different things, as they do for
    #: `resolve_oprd_loci`: ``None`` is "no locus search ran", and the report
    #: says the table is not produced and why. There is **no** contracted path
    #: for these verdicts - `papipeline/execution/contracts.py` declares no
    #: oprD table and `workflow/Snakefile` declares no oprD output - so they
    #: cannot be read from disk and this field is the only route to them. That
    #: is a finding, not a design: see the BACKLOG entry for a contracted oprD
    #: verdict table.
    oprd_calls: Optional[Sequence[Any]] = None


def markdown_table(
    columns: Sequence[str], rows: Sequence[Mapping[str, Any]]
) -> str:
    """Render a Markdown table. Empty input yields an explicit note."""
    if not rows:
        return "_No rows._\n"
    header = "| " + " | ".join(columns) + " |"
    divider = "| " + " | ".join("---" for _ in columns) + " |"
    lines = [header, divider]
    for row in rows:
        cells = []
        for column in columns:
            value = row.get(column)
            cells.append("-" if value is None else str(value))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def _kv_table(mapping: Mapping[str, Any]) -> str:
    return markdown_table(["field", "value"], [{"field": k, "value": v} for k, v in mapping.items()])


#: The exact wording the report uses when the source recorded no AST standard.
#:
#: Named, not spelled inline, because it is an assertion a reader may rely on:
#: the pipeline must say that the standard was not recorded rather than merely
#: omitting a standard column and leaving the reader to infer the gap. Changing
#: it is a content change to every report, so it is declared once here.
NO_STANDARD_RECORDED = "testing standard not recorded in source data"


def phenotype_report(context: "ReportContext") -> Optional[ProvenanceReport]:
    """Stage 11's provenance, or ``None`` when no phenotype was loaded.

    **This is the production caller of
    :func:`~papipeline.stages.phenotype.provenance_report`.** It is called from
    :func:`build_markdown`, so it runs on the reporting path of every run that
    carries phenotype calls - not only from a test.

    ``None`` rather than an empty report when there are no calls: a
    zero-row report and a report about a cohort whose provenance is wholly
    absent are different claims, and a table of zeros would read as the first.
    """
    calls = context.phenotype_calls
    if calls is None:
        return None
    return provenance_report(calls)


def phenotype_provenance_section(
    config: PipelineConfig, report: ProvenanceReport
) -> Tuple[str, str]:
    """Render stage 11's provenance: what the cohort is, and what it is not.

    Two things a reader cannot get from the association table.

    1. **Which categories are actually in the contrast.** R and S are the
       primary phenotype (ruling R7). ``I`` and ``SDD`` are real, measured
       categories that the model refuses, and they are counted here rather
       than merged - a report that showed only R and S would let a reader
       assume the cohort had no intermediate isolates.
    2. **Whether the calls are comparable.** A call whose laboratory recorded
       no standard cannot be checked against a standard, and this report says
       so in those words rather than by leaving the field blank.
    """
    held = report.excluded_counts

    lines: List[str] = []
    lines.append(
        f"Stage 11 loaded {report.total} phenotype record(s). "
        f"**{report.analysed} enter the {PRIMARY_POSITIVE}-vs-"
        f"{PRIMARY_NEGATIVE} contrast** ({PRIMARY_POSITIVE}={report.positive}, "
        f"{PRIMARY_NEGATIVE}={report.negative})."
    )
    lines.append("")
    if report.excluded:
        lines.append(
            f"{report.excluded} record(s) are excluded from that contrast and "
            f"counted here: "
            f"{', '.join(f'{k}={v}' for k, v in held.items())}. The excluded "
            f"categories are {', '.join(EXCLUDED_CATEGORIES)}, and ruling R7 "
            f"forbids merging any of them into {PRIMARY_POSITIVE} or "
            f"{PRIMARY_NEGATIVE}: an intermediate is not a resistance call and "
            "an SDD is not a plain susceptible. They are reported as measured "
            "and left out of the comparison."
        )
    else:
        lines.append(
            f"No record is excluded from that contrast: this cohort carries no "
            f"{', '.join(EXCLUDED_CATEGORIES)} call. Had it carried one it would "
            "be counted here and left out of the comparison, never merged into "
            f"{PRIMARY_POSITIVE} or {PRIMARY_NEGATIVE} (ruling R7)."
        )
    lines.append("")
    lines.append(_kv_table({
        "phenotype records loaded": report.total,
        f"{PRIMARY_POSITIVE} (primary, positive arm)": report.positive,
        f"{PRIMARY_NEGATIVE} (primary, negative arm)": report.negative,
        "entering the primary contrast": report.analysed,
        "excluded and counted, not merged": report.excluded,
        "carrying a measured MIC": report.with_mic,
    }))
    lines.append("")
    lines.append(_kv_table({
        "rows lacking AST method": report.lacking_method,
        "rows lacking AST standard": report.lacking_standard,
        "rows lacking AST edition": report.lacking_edition,
    }))

    if not report.standard_recorded_everywhere:
        if report.lacking_standard == report.total:
            lines.append("")
            lines.append(
                f"**The {NO_STANDARD_RECORDED}.** The phenotype table carries "
                "no `ast_standard` column (or carries it empty), so which "
                "breakpoint standard these calls were made against is unknown "
                "and cannot be recovered from this run. Calls made under "
                "different standards are not interchangeable - CLSI and EUCAST "
                "do not place the imipenem S/I/R boundaries for *Pseudomonas "
                "aeruginosa* at the same MIC - so these calls are reported as "
                "the source laboratory made them and are not re-interpreted "
                "against any threshold."
            )
        else:
            lines.append("")
            lines.append(
                f"{report.lacking_standard} of {report.total} records carry no "
                "AST standard, so **for those rows the testing standard not "
                "recorded in source data**. Calls made under different "
                "standards are not interchangeable, so these records are not "
                "assumed comparable with the rest of the cohort."
            )

    status = BreakpointStatus.from_config(config)
    lines.append("")
    lines.append(
        f"No S/I/R threshold was applied to these calls: `breakpoints."
        f"thresholds` is empty, so no MIC was reinterpreted into a category and "
        "the source calls stand unchanged. `config/science.yaml` names "
        f"{status.standard} ({status.edition}) as the standard this repository "
        f"*intends* the source laboratories to have used for {status.agent}. "
        "That is an **intention held by this repository, not a provenance fact "
        "about these isolates**: it is not present in the source data, and no "
        "file in this repository records which standard the source "
        f"laboratories actually applied. {status.standard} is a licensed "
        f"document; supply the {status.agent} S/I/R values to enable "
        "re-interpretation."
    )

    return ("Phenotype provenance", "\n".join(lines))


def is_smoke_report(config: PipelineConfig) -> bool:
    """Whether this report is from the bounded smoke run rather than a real analysis.

    Keyed on the **overlay**, not the mode. Both a smoke run and a full-cohort
    analysis are ``RunMode.REAL``, so mode cannot tell them apart - but only the
    smoke overlay sets ``smoke_genome_dir``. A full-cohort report stamped
    "SMOKE TEST" would be confidently wrong, which is the same failure class as
    omitting the marker from a smoke one and no less bad.

    Delegates to :meth:`PipelineConfig.is_smoke_overlay`, which manifest
    discovery also uses, so the two cannot disagree about which runs are
    bounded.
    """
    return config.is_smoke_overlay()


# --------------------------------------------------------------------------
# Sections read back from the run's own output tables
# --------------------------------------------------------------------------

#: Column of the stage-4 table that marks a point-mutated determinant row.
#:
#: **The stage-4 table carries no `point_mutation` flag, and nothing may
#: pretend otherwise.** `determinant_type` is AMRFinderPlus' `Element Type`
#: column, which is `AMR` for a point mutation and for an intact gene alike
#: (`papipeline/stages/amr.py:406`), so it does not discriminate. What does is
#: `variant`: `amr.point_variants_in_report` populates it from the report's own
#: `Subtype` column, accepting only `POINT` and `POINT_DISRUPT`
#: (`papipeline/stages/amr.py:69-73`). So "variant non-empty" means "the tool
#: told us this element is mutated", which is the strongest statement this table
#: supports.
#:
#: It does NOT mean disruptive. In this cohort `POINT` and `POINT_DISRUPT`
#: partition catalogued and uncatalogued alleles respectively, with the
#: frameshift `oprD_Y237TerfsTer0` among the uncatalogued ones
#: (`pa-artifacts/round12/AMRX.md` §4.3), and the table records neither. Any
#: claim about which substitutions disrupt oprD needs evidence this file does
#: not contain.
POINT_MUTATION_COLUMN = "variant"

#: The oprD gene symbol, as the AMR table spells it.
OPRD_GENE = "oprD"

#: How many rows a per-sample table may contribute to a summary before the
#: report says "the rest" rather than printing all of them.
TOP_N = 25


def _counts_body(
    counter: Mapping[str, int], *, label: str, total_label: str = "rows"
) -> str:
    """A two-column count table, largest first, with the total stated."""
    if not counter:
        return (
            f"No {label} were recorded. A table with no {label} is not the same "
            f"as a stage that did not run; this section is rendered because the "
            f"file behind it was read and held rows, and every value in it was "
            f"blank.\n"
        )
    ordered = sorted(counter.items(), key=lambda kv: (-kv[1], str(kv[0])))
    shown = ordered[:TOP_N]
    omitted = len(ordered) - len(shown)
    body = markdown_table(
        [label, total_label], [{label: k, total_label: v} for k, v in shown]
    )
    if omitted:
        body += (
            f"\n_{omitted} further {label}(s) present in the file and not shown "
            f"here; the file is the record._\n"
        )
    return body


def _absent_body(table: TableRead, *, what: str) -> str:
    """The body of a section whose table could not be read.

    Always says *what* was wanted and *why* it is not there. Never blank, never
    a zero-row table: a header-only file passed to ``markdown_table`` renders as
    `_No rows._`, which reads as a measurement of zero rather than as a gap in
    the run.
    """
    return (
        f"**{NOT_PRODUCED}.** {table.reason}\n\n"
        f"This section reports {what}. With the table absent there is nothing "
        f"to report, and the run is not making a claim about it either way.\n"
    )


def _flagged(body: str, n: int) -> str:
    """Prefix a statistic-bearing body with the R16 power flag."""
    return f"{power_flag_line(n)}\n\n{body}"


def _stage_status_rows(
    context: ReportContext, tables: RunTables
) -> Tuple[List[Dict[str, Any]], str]:
    """Per-stage status, and a sentence saying where the statuses came from.

    Two sources, and which one is in use is stated because they make different
    claims. The orchestrator's map knows about refusals, skips and unbuilt
    stages. The disk-derived fallback knows only whether a file is there, so a
    stage refused upstream and a stage that crashed look identical - which is why
    it is labelled rather than presented as the run's own record.
    """
    from .report_tables import STAGE_TABLE_KEYS

    rows: List[Dict[str, Any]] = []
    if context.stage_status:
        provenance_note = (
            "Status is the orchestrator's own record of this run "
            "(`ReportContext.stage_status`), which knows why a stage was "
            "refused or skipped."
        )
        for stage, state in context.stage_status.items():
            # `is not None`, not truthiness: `TableRead.__bool__` is `present`,
            # so testing the object would report "no table declared for this
            # stage" about every declared table that is simply absent from disk.
            table = tables.get(stage)
            rows.append(
                {
                    "stage": stage,
                    "status": state,
                    "reason": context.stage_refusals.get(stage)
                    or (
                        ""
                        if table is not None and table.present
                        else (
                            table.reason
                            if table is not None
                            else "no table declared for this stage"
                        )
                    ),
                }
            )
        return rows, provenance_note

    provenance_note = (
        "Status is **derived from disk**, not from the run's own record: a stage "
        "is called `completed` when its contracted output table is present and "
        "readable. That cannot distinguish a stage that was refused, one that "
        "was skipped and one that crashed after writing its table, so the "
        "orchestrator should hand over `ReportContext.stage_status` and "
        "`stage_refusals` for a run whose stage outcomes are in dispute."
    )
    for stage, _key in STAGE_TABLE_KEYS:
        table = tables[stage]
        rows.append(
            {
                "stage": stage,
                "status": "completed" if table.present else NOT_PRODUCED,
                "reason": "" if table.present else table.reason,
            }
        )
    return rows, provenance_note


def run_summary_section(context: ReportContext, tables: RunTables) -> str:
    """What ran, on how many isolates, and what each stage did."""
    rows, note = _stage_status_rows(context, tables)
    lines = [
        _kv_table(
            {
                "run mode": context.mode.value,
                "isolates in this run": context.n_samples,
                "antibiotic": context.antibiotic,
                "generated (UTC)": context.generated_at,
                "pipeline version": __version__,
                "stage tables read from disk": f"{tables.n_present} present, "
                f"{tables.n_absent} not produced",
                "table directory": str(tables.stage_dir),
            }
        ),
        "",
        f"### Per-stage status\n",
        f"{note}\n",
        markdown_table(["stage", "status", "reason"], rows),
    ]
    return "\n".join(lines)


def cohort_section(context: ReportContext, tables: RunTables) -> str:
    """Cohort size, sequence-type distribution and phenotype counts."""
    validation = tables["validation"]
    mlst = tables["mlst"]
    phenotype = tables["phenotype"]

    lines: List[str] = []
    lines.append(
        _kv_table(
            {
                "isolates in this run": context.n_samples,
                "isolates with a validation row": len(
                    validation.distinct("sample_id")
                ) if validation.present else NOT_PRODUCED,
                "assemblies screened (stage 2 rows)": (
                    len(annotation.distinct("sample_id"))
                    if (annotation := tables["annotation"]).present
                    else NOT_PRODUCED
                ),
            }
        )
    )

    lines.append("")
    if mlst.present:
        st_counts = mlst.count("ST")
        lines.append("### Sequence types (stage 3)\n")
        lines.append(_counts_body(st_counts, label="ST", total_label="isolates"))
    else:
        lines.append("### Sequence types (stage 3)\n")
        lines.append(_absent_body(mlst, what="the cohort's sequence-type distribution"))

    lines.append("")
    if phenotype.present:
        pheno_counts = phenotype.count("phenotype")
        lines.append("### Susceptibility categories (stage 11)\n")
        lines.append(
            _counts_body(pheno_counts, label="category", total_label="isolates")
        )
        lines.append("")
        lines.append(
            f"Source: `{phenotype.path}`. "
            f"{NO_STANDARD_RECORDED} is stated in the phenotype provenance "
            f"section above, because this table carries no AST standard column "
            f"and the categories are therefore the source laboratory's own "
            f"labels rather than calls re-interpreted against a breakpoint "
            f"standard.\n"
        )
    else:
        lines.append("### Susceptibility categories (stage 11)\n")
        lines.append(
            _absent_body(phenotype, what="the cohort's susceptibility categories")
        )

    return _flagged("\n".join(lines), context.n_samples)


def _amr_body(tables: RunTables) -> str:
    amr = tables["amr"]
    if not amr.present:
        return _absent_body(amr, what="the cohort's antimicrobial-resistance calls")

    lines: List[str] = []
    samples = set(amr.distinct("sample_id"))
    is_point = [
        bool(str(r.get(POINT_MUTATION_COLUMN) or "").strip()) for r in amr.rows
    ]
    point_rows = [r for r, point in zip(amr.rows, is_point) if point]
    intact_rows = [r for r, point in zip(amr.rows, is_point) if not point]
    genes = amr.distinct("gene")
    point_genes = sorted(
        {str(r.get("gene")) for r in point_rows if str(r.get("gene") or "").strip()}
    )

    lines.append(
        _kv_table(
            {
                "determinant rows": amr.n_rows,
                "isolates with at least one determinant": len(samples),
                "distinct genes": len(genes),
                "intact-gene rows": len(intact_rows),
                "point-mutation rows": len(point_rows),
                "genes carrying a point mutation": len(point_genes),
                "source": str(amr.path),
            }
        )
    )

    lines.append("")
    lines.append("### Determinants by gene\n")
    gene_counts = amr.count("gene")
    lines.append(
        _counts_body(gene_counts, label="gene", total_label="determinant rows")
    )

    lines.append("")
    lines.append(
        f"### Point mutations ({len(point_rows)} row(s))\n"
    )
    if point_rows:
        lines.append(
            markdown_table(
                [
                    "sample_id",
                    "gene",
                    "variant",
                    "determinant",
                    "determinant_type",
                    "identity_pct",
                    "coverage_pct",
                    "claim_status",
                ],
                point_rows,
            )
        )
        lines.append(
            f"\nA row is in this table because its `{POINT_MUTATION_COLUMN}` "
            f"column is non-empty, which `stages.amr` fills in only for "
            f"AMRFinderPlus `Subtype` values `POINT` or `POINT_DISRUPT`. "
            f"`determinant_type` is the tool's `Element Type` and is `AMR` for "
            f"both a point mutation and an intact gene, so it does not "
            f"discriminate and is shown for completeness only.\n\n"
            f"**These rows are not evidence that the mutation disrupts the "
            f"gene.** This table records no substitution class and no "
            f"reading-frame assessment. The structural reading-frame verdicts "
            f"are in the oprD verdict table below, and they are a separate "
            f"evidence source with its own contract.\n"
        )
    else:
        lines.append(
            f"No determinant row carries a `{POINT_MUTATION_COLUMN}` value, so "
            f"this table is empty of point mutations. That is what the run "
            f"recorded; whether the screen was configured to report point "
            f"mutations at all is a property of the tool invocation, and the "
            f"stage writes its organism flag into the file's `#` provenance "
            f"header for exactly this reason.\n"
        )

    oprd_rows = [
        r for r in amr.rows if str(r.get("gene") or "").strip() == OPRD_GENE
    ]
    lines.append("")
    lines.append(f"### {OPRD_GENE} calls from the AMR screen\n")
    if oprd_rows:
        lines.append(
            markdown_table(
                [
                    "sample_id",
                    "determinant",
                    "variant",
                    "identity_pct",
                    "coverage_pct",
                    "claim_status",
                ],
                oprd_rows,
            )
        )
    else:
        lines.append(
            f"No `{OPRD_GENE}` row is in `{amr.path}`. Absence from the AMR "
            f"table is **not** evidence that the locus is absent: the "
            f"structural verdict table below is the only source that may say "
            f"so, and only for an isolate it actually searched.\n"
        )
    return "\n".join(lines)


def amr_section(context: ReportContext, tables: RunTables) -> str:
    """Antimicrobial-resistance summary, point mutations included."""
    amr = tables["amr"]
    n = len(amr.distinct("sample_id")) if amr.present else context.n_samples
    return _flagged(_amr_body(tables), n)


def _virulence_body(tables: RunTables) -> str:
    vf = tables["virulence"]
    if not vf.present:
        return _absent_body(vf, what="the cohort's virulence-factor calls")
    lines = [
        _kv_table(
            {
                "virulence rows": vf.n_rows,
                "isolates with at least one factor": len(vf.distinct("sample_id")),
                "distinct virulence factors": len(vf.distinct("virulence_factor")),
                "distinct genes": len(vf.distinct("gene")),
                "source": str(vf.path),
            }
        ),
        "",
        "### Virulence factors detected\n",
        _counts_body(
            vf.count("virulence_factor"),
            label="virulence_factor",
            total_label="rows",
        ),
        "",
        f"A factor detected by homology is a factor **detected**. "
        f"`docs/scientific_rules.md` rule 1: gene presence is never reported "
        f"as phenotypic virulence, and this table records no phenotype.\n",
    ]
    return "\n".join(lines)


def virulence_section(context: ReportContext, tables: RunTables) -> str:
    vf = tables["virulence"]
    n = len(vf.distinct("sample_id")) if vf.present else context.n_samples
    return _flagged(_virulence_body(tables), n)


#: The oprD verdict table's columns, in the order a reader needs them.
OPRD_VERDICT_COLUMNS: Tuple[str, ...] = (
    "sample_id",
    "verdict",
    "lesion_type",
    "position",
    "truncation_aa",
    "identity_pct",
    "coverage_pct",
)


def oprd_verdict_rows(calls: Sequence[Any]) -> List[Dict[str, Any]]:
    """One row per isolate's oprD verdict, with the evidence behind it.

    Accepts both shapes `resolve_oprd_loci` can hand over:

    * a :class:`~papipeline.adapters.oprd_locus.LocusResolution` (blastp against
      the isolate proteome), which carries no reading-frame evidence, so its
      lesion columns are empty and say so;
    * a :class:`~papipeline.adapters.oprd_locus.StructuralCall` (tblastn against
      the assembly), which carries a
      :class:`~papipeline.adapters.oprd_locus.LocusEvidence` and therefore a
      lesion type, a stop position and a truncation in residues.

    A structural call is the one that judges the reading frame, and RULING 2
    defines "disrupted" as a frameshift **or** a premature stop - so an evidence
    cell that is empty is not a missing measurement to be filled in later, it is
    the reason the verdict cannot be a structural one.
    """
    rows: List[Dict[str, Any]] = []
    for index, call in enumerate(calls or ()):
        evidence = getattr(call, "evidence", None)
        sample_id = str(getattr(call, "sample_id", index))
        row: Dict[str, Any] = {
            "sample_id": sample_id,
            "verdict": getattr(call, "verdict", "unknown"),
            "lesion_type": getattr(evidence, "lesion_type", None),
            "position": (
                getattr(evidence, "internal_stop_position_aa", None)
                if evidence is not None
                else None
            ),
            "truncation_aa": (
                getattr(evidence, "truncation_aa", None)
                if evidence is not None
                else None
            ),
            "identity_pct": getattr(call, "identity_pct", None),
            "coverage_pct": getattr(call, "coverage_pct", None),
        }
        rows.append(row)
    return rows


#: Why the oprD verdict table cannot be read from disk.
#:
#: Stated rather than left implicit, because "not produced" without a reason is
#: indistinguishable from a stage that was never built.
OPRD_NOT_ON_DISK = (
    f"{NOT_PRODUCED}: no contracted path exists for the oprD locus verdict. "
    f"`papipeline/execution/contracts.py` declares no oprD table in "
    f"STAGE_TABLES or INTERNAL_TABLES and `workflow/Snakefile` declares no "
    f"oprD output, so `run.py`'s in-memory `oprd_locus_resolution` is the only "
    f"copy of these verdicts and nothing on disk holds them."
)


def _oprd_body(context: ReportContext, tables: RunTables) -> str:
    if context.oprd_calls is None:
        return (
            f"**{OPRD_NOT_ON_DISK}** The orchestrator also handed over no "
            f"`oprd_calls`, so the run holds no verdicts to report.\n\n"
            f"This is a gap in the run's durable record rather than a result. "
            f"The central negative of this study turns on oprD carriage, and a "
            f"verdict that exists only in memory is a verdict that cannot be "
            f"re-derived from the artefacts the run left behind.\n"
        )
    if not context.oprd_calls:
        return (
            f"**{NOT_PRODUCED}.** The run completed a locus search that "
            f"returned no isolate at all. An empty search and a search that "
            f"never ran are different claims and this states the first.\n"
        )

    rows = oprd_verdict_rows(context.oprd_calls)
    n = len(rows)
    lines = [
        _kv_table(
            {
                "isolates with an oprD verdict": n,
                "verdicts": len({r["verdict"] for r in rows}),
                "evidence source": (
                    "structural (blastp/tblastn against the pinned PAO1 PA0958 "
                    "sequence)"
                ),
            }
        ),
        "",
        "### oprD verdict and its evidence\n",
        markdown_table(list(OPRD_VERDICT_COLUMNS), rows),
        "",
        "### Verdicts by call\n",
        _counts_body(
            _verdict_counts(rows),
            label="verdict",
            total_label="isolates",
        ),
        "",
        f"**How to read the evidence columns.** `verdict` is the decision. "
        f"`lesion_type` is what broke the reading frame - RULING 2 defines "
        f"`disrupted` as a frameshift **or** a premature stop and nothing else, "
        f"decided by lesion type and never by ORF length. `position` is the "
        f"internal stop codon in amino-acid coordinates on the 444-residue PAO1 "
        f"reference; `truncation_aa` counts the reference residues the ORF does "
        f"not encode. Truncation is **evidence, not a decision**: a locus with a "
        f"large truncation and a native stop is `intact`.\n\n"
        f"An empty lesion cell means the verdict came from a presence search "
        f"that does not judge the reading frame, not that the locus was found "
        f"whole. Only a structural call may report `absent`, and only for an "
        f"isolate it actually searched.\n",
    ]
    return _flagged("\n".join(lines), n)


def _verdict_counts(rows: Sequence[Mapping[str, Any]]) -> Dict[str, int]:
    counter: Dict[str, int] = {}
    for row in rows:
        key = str(row.get("verdict"))
        counter[key] = counter.get(key, 0) + 1
    return counter


def oprd_section(context: ReportContext, tables: RunTables) -> str:
    return _oprd_body(context, tables)


def _variant_kind(ref: Any, alt: Any) -> str:
    """SNV, MNV, insertion, deletion, or `unknown` for an unusable row.

    Decided from the REF and ALT alleles alone, on the VCF convention: a single
    base on each side is a substitution, a single base against several is an
    insertion, several against a single base is a deletion, and equal lengths
    above one are a multi-nucleotide substitution.
    """
    r = str(ref or "")
    a = str(alt or "")
    if not r or not a:
        return "unknown"
    if len(r) == 1 and len(a) == 1:
        return "SNV"
    if len(r) == 1:
        return "insertion"
    if len(a) == 1:
        return "deletion"
    if len(r) == len(a):
        return "MNV"
    return "complex"


def variants_body(tables: RunTables) -> str:
    variants = tables["variants"]
    cohort = tables["cohort_variants"]
    if not variants.present:
        return _absent_body(variants, what="the cohort's variant calls")

    kinds: Dict[str, int] = {}
    per_sample: Dict[str, int] = {}
    for row in variants.rows:
        kind = _variant_kind(row.get("ref"), row.get("alt"))
        kinds[kind] = kinds.get(kind, 0) + 1
        sample = str(row.get("sample_id") or "")
        per_sample[sample] = per_sample.get(sample, 0) + 1

    snv = kinds.get("SNV", 0)
    non_snv = variants.n_rows - snv

    lines = [
        _kv_table(
            {
                "per-isolate variant rows": variants.n_rows,
                "isolates with at least one call": len(per_sample),
                "SNV (REF and ALT both one base)": snv,
                "non-SNV (every other kind)": non_snv,
                "source": str(variants.path),
            }
        ),
        "",
        "### Calls by variant kind\n",
        _counts_body(kinds, label="kind", total_label="rows"),
        "",
        f"**SNV versus non-SNV is a statement about the allele strings, not "
        f"about biology.** A row is an SNV because its REF and ALT are one base "
        f"each; the classification is made here, from the two columns the file "
        f"carries, and it says nothing about whether the call is correct, "
        f"whether the position is variable across the cohort, or whether an "
        f"indel-sized event is a caller artefact. Per-isolate calls are "
        f"confounded with species-wide fixed differences - "
        f"`contracts.py` says so on the `variants` entry - so a cohort-wide "
        f"merge is required before any of these are GWAS features.\n",
    ]

    lines.append("")
    lines.append("### Cohort-wide merged variants (stage 6a)\n")
    if cohort.present:
        ac = [r.get("ac") for r in cohort.rows]
        af = [r.get("af") for r in cohort.rows]
        lines.append(
            _kv_table(
                {
                    "cohort variant rows": cohort.n_rows,
                    "distinct chromosomes": len(cohort.distinct("chrom")),
                    "allele count range": (
                        f"{min(int(a) for a in ac if a is not None)}-"
                        f"{max(int(a) for a in ac if a is not None)}"
                        if any(a is not None for a in ac)
                        else "not recorded"
                    ),
                    "allele frequency range": (
                        f"{min(float(a) for a in af if a is not None)}-"
                        f"{max(float(a) for a in af if a is not None)}"
                        if any(a is not None for a in af)
                        else "not recorded"
                    ),
                    "source": str(cohort.path),
                }
            )
        )
    else:
        lines.append(_absent_body(cohort, what="the cohort-wide variant merge"))

    return "\n".join(lines)


def variants_section(context: ReportContext, tables: RunTables) -> str:
    variants = tables["variants"]
    n = len(variants.distinct("sample_id")) if variants.present else context.n_samples
    return _flagged(variants_body(tables), n)


def pangenome_section(context: ReportContext, tables: RunTables) -> str:
    summary = tables["pangenome"]
    if not summary.present:
        return _flagged(
            _absent_body(summary, what="the pangenome summary"), context.n_samples
        )
    lines = [
        f"Source: `{summary.path}`.\n",
        _kv_table(
            {
                str(r.get("metric")): r.get("value") for r in summary.rows
            }
        ),
    ]
    for key, label in (
        ("core_genes", "Core genes (stage 7)"),
        ("accessory_genes", "Accessory genes (stage 7)"),
        ("gene_presence_absence", "Gene presence/absence matrix"),
    ):
        table = tables[key]
        lines.append("")
        lines.append(f"### {label}\n")
        if table.present:
            lines.append(f"{table.n_rows} row(s) in `{table.path}`.")
        else:
            lines.append(_absent_body(table, what=label.lower()))
    return _flagged("\n".join(lines), context.n_samples)


def tree_similarity_section(context: ReportContext, tables: RunTables) -> str:
    """Tree and distance summary. The distances' UNITS are stated, not implied."""
    from .similarity import DISTANCE_UNIT, METHOD, NOT_A_SNP_COUNT

    tree = tables["phylogeny"]
    similarity = tables["similarity"]
    units = tables.similarity_units

    lines: List[str] = []
    if tree.present:
        lines.append("### Tree (stage 9)\n")
        lines.append(markdown_table(list(tree.rows[0].keys()), list(tree.rows)))
        lines.append(f"\nSource: `{tree.path}`.\n")
    else:
        lines.append("### Tree (stage 9)\n")
        lines.append(_absent_body(tree, what="the tree summary"))

    lines.append("")
    lines.append("### Alignment summary (stage 9)\n")
    alignment = tables["alignment_summary"]
    if alignment.present:
        lines.append(
            _kv_table({str(r.get("metric")): r.get("value") for r in alignment.rows})
        )
    else:
        lines.append(
            _absent_body(
                alignment,
                what=(
                    "the alignment summary. Stage 9 writes it only when "
                    "`alignment_info` is non-empty, which in turn requires both "
                    "`core_alignment.fasta` and `core_snp_alignment.fasta`; the "
                    "seam analysis recorded that neither has an in-repo producer "
                    "in REAL (`pa-artifacts/round12/seam-matrix.md` §G5)"
                ),
            )
        )

    lines.append("")
    lines.append("### Pairwise distances (stage 10)\n")
    if similarity.present:
        lines.append(
            _kv_table(
                {
                    "rows": similarity.n_rows,
                    "isolates": len(similarity.distinct("sample_id")),
                    "source": str(similarity.path),
                }
            )
        )
    else:
        lines.append(_absent_body(similarity, what="the pairwise distance matrix"))

    lines.append("")
    lines.append("#### Units of the distances\n")
    if units is not None and units.present:
        lines.append(
            _kv_table({str(k): v for k, v in units.rows[0].items()})
        )
        lines.append(f"\nSource: `{units.path}`.\n")
        lines.append(f"**Unit: {DISTANCE_UNIT}.** {NOT_A_SNP_COUNT}\n")
    else:
        reason = units.reason if units is not None else "no sidecar was read"
        lines.append(
            f"**Unit: {DISTANCE_UNIT}** - this is the unit stage 10 computes in "
            f"(`stages.similarity.DISTANCE_UNIT`), stated here because the "
            f"report must not leave a reader to guess.\n\n"
            f"**{NOT_PRODUCED}: the units sidecar.** {reason} The distances "
            f"above therefore carry their unit only in this sentence. Stage 10 "
            f"writes `similarity.units.json` beside the matrix for exactly this "
            f"purpose, and the sidecar's reader had no caller until this "
            f"section: `pa-artifacts/round12/seam-matrix.md` §1 recorded that "
            f"`similarity.read_units_sidecar` had zero callers, so the file was "
            f"written on every similarity run and read by nothing.\n\n"
            f"{NOT_A_SNP_COUNT}\n"
        )

    lines.append("")
    lines.append(
        f"The method recorded by stage 10 is `{METHOD}`, and the distance is "
        f"the sum of branch lengths along the tree - so it is only as good as "
        f"the tree above, and a distance matrix computed over a tree with no "
        f"rooting rationale will not be symmetric under any re-rooting that "
        f"matters for these numbers.\n"
    )
    return _flagged("\n".join(lines), _tree_n(context, tree, similarity))


def _tree_n(context: ReportContext, tree: TableRead, similarity: TableRead) -> int:
    """The cohort size these tree/distance numbers were computed on.

    Taken from the files when they carry one, so the R16 flag states the ``n``
    the numbers actually have rather than the ``n`` the orchestrator believed
    the run had. Those disagree whenever a stage covered a subset of the cohort,
    and a flag that overstates the cohort overstates the power.
    """
    if similarity.present:
        n = len(similarity.distinct("sample_id"))
        if n:
            return n
    if tree.present and tree.rows:
        tips = _as_float(tree.rows[0].get("n_tips"))
        if tips:
            return int(tips)
    return context.n_samples


def _as_float(value: Any) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _gwas_threshold(config: PipelineConfig) -> float:
    """The GWAS significance threshold, from the config accessor the stage uses.

    Read from ``config.gwas.significance_threshold`` - the same attribute
    `stages.gwas` reads - rather than from a raw dict lookup, so a report section
    cannot disagree with the stage that produced the numbers it is summarising.
    """
    return float(config.gwas.significance_threshold)


def _cooccurrence_threshold(config: PipelineConfig) -> float:
    """The co-occurrence threshold, from the key `stages.cooccurrence` reads."""
    return float(
        (config.raw.get("cooccurrence") or {}).get("significance_threshold", 0.05)
    )


def _association_body(context: ReportContext, tables: RunTables) -> str:
    gwas = tables["gwas"]
    convergence = tables["convergence"]
    cooccurrence = tables["cooccurrence"]

    lines: List[str] = []

    lines.append("### Stage 12 - GWAS associations\n")
    if gwas.present:
        adjusted = [
            r.get("adjusted_p_value")
            for r in gwas.rows
            if r.get("adjusted_p_value") is not None
        ]
        floats = [f for f in (_as_float(a) for a in adjusted) if f is not None]
        best = min(floats) if floats else None
        threshold = _gwas_threshold(context.config)
        below = [f for f in floats if f <= threshold]
        untested = gwas.n_rows - len(adjusted)
        lines.append(
            _kv_table(
                {
                    "association rows": gwas.n_rows,
                    "distinct features": len(gwas.distinct("feature")),
                    "rows carrying an adjusted p-value": len(adjusted),
                    f"rows at or below the configured threshold ({threshold:g})": len(
                        below
                    ),
                    "rows with no adjusted p-value (never testable)": untested,
                    "strongest adjusted p-value": (
                        f"{best:g}" if best is not None else "not recorded"
                    ),
                    "source": str(gwas.path),
                }
            )
        )
        if gwas.rows:
            lines.append("")
            lines.append(
                markdown_table(
                    [
                        "feature",
                        "feature_type",
                        "effect",
                        "p_value",
                        "adjusted_p_value",
                        "frequency",
                    ],
                    gwas.rows,
                )
            )
        lines.append(
            f"\n**Association is not causation.** `docs/scientific_rules.md` "
            f"rule 2, and the word 'causal' is not a finding anywhere in this "
            f"pipeline. Stage 12's `ReferenceEngine` applies no population-"
            f"structure correction, so a feature carried almost entirely by one "
            f"lineage will test significant for reasons that have nothing to do "
            f"with susceptibility; the `lineage_distribution` column is the "
            f"evidence for that, and it is in the file.\n"
        )
    else:
        lines.append(_absent_body(gwas, what="the cohort's association results"))

    lines.append("")
    lines.append("### Stage 13 - convergence\n")
    if convergence.present:
        lines.append(
            _kv_table(
                {
                    "determinant rows": convergence.n_rows,
                    "distinct determinants": len(convergence.distinct("determinant")),
                    "source": str(convergence.path),
                }
            )
        )
        lines.append("")
        lines.append(
            markdown_table(
                [
                    "determinant",
                    "independent_lineages",
                    "branch_count",
                    "distribution",
                    "convergence_category",
                ],
                convergence.rows,
            )
        )
        lines.append(
            f"\n`convergence_category` is the tool's own classification of how "
            f"widely a determinant is spread across the tree. A determinant "
            f"classified `UNKNOWN` was not classified - it is the sentinel this "
            f"stage's own REAL-mode guard treats as an absence "
            f"(`stages.reporting._REPORTING_SENTINELS`).\n"
        )
    else:
        lines.append(_absent_body(convergence, what="the convergence classification"))

    lines.append("")
    lines.append("### Stage 14 - co-occurrence\n")
    if cooccurrence.present:
        adjusted = [
            r.get("adjusted_p_value") for r in cooccurrence.rows
            if r.get("adjusted_p_value") is not None
        ]
        threshold = _cooccurrence_threshold(context.config)
        significant = [
            r for r in cooccurrence.rows
            if _as_float(r.get("adjusted_p_value")) is not None
            and float(_as_float(r.get("adjusted_p_value"))) <= threshold
        ]
        lines.append(
            _kv_table(
                {
                    "pairs tested": cooccurrence.n_rows,
                    "feature types": len(cooccurrence.distinct("feature_type")),
                    f"pairs at or below the configured threshold "
                    f"({threshold:g})": len(significant),
                    "pairs carrying an adjusted p-value": len(adjusted),
                    "source": str(cooccurrence.path),
                }
            )
        )
        if significant:
            lines.append("")
            lines.append(
                markdown_table(
                    [
                        "feature_a",
                        "feature_b",
                        "feature_type",
                        "n_both",
                        "statistic_value",
                        "adjusted_p_value",
                    ],
                    significant,
                )
            )
        lines.append(
            f"\n**Association only.** Two features observed together in a "
            f"cohort of this size is not interaction, and the "
            f"`interpretation_limit` column in the file is where the tool "
            f"records that.\n"
        )
    else:
        lines.append(
            _absent_body(cooccurrence, what="the co-occurrence results")
        )

    return "\n".join(lines)


def _as_float(value: Any) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def association_section(context: ReportContext, tables: RunTables) -> str:
    """GWAS, convergence and co-occurrence. Each carries the R16 flag."""
    return _flagged(_association_body(context, tables), context.n_samples)


def _provenance_body(context: ReportContext, tables: RunTables) -> str:
    """Tool versions, plus the phenotype AST provenance the table dropped.

    ``stages.phenotype.to_row`` emits thirteen keys and the stage-11 table is
    written with eight, so ``ast_method``, ``ast_standard`` and ``ast_edition``
    never reach disk (`pa-artifacts/round12/seam-matrix.md` §f). They are still
    in `ReportContext.phenotype_calls`, so the report can state what they were -
    or state that they were not recorded, which is the answer the committed
    fixture requires.
    """
    lines: List[str] = []
    pheno = phenotype_report(context)
    if pheno is None:
        lines.append(
            f"**{NOT_PRODUCED}: the phenotype testing-standard provenance.** No "
            f"stage 11 calls were handed to this report, so nothing at all is "
            f"known about which testing standard produced them. That is "
            f"*stronger* than an unrecorded standard - it says the run holds no "
            f"provenance, rather than that the source recorded none - and it is "
            f"why no S/I/R interpretation appears anywhere in this report.\n"
        )
    else:
        rows = [
            {
                "field": "records",
                "value": pheno.total,
                "note": "phenotype records loaded by stage 11",
            },
            {
                "field": "rows lacking ast_method",
                "value": pheno.lacking_method,
                "note": "how the laboratory made the call is unknown",
            },
            {
                "field": "rows lacking ast_standard",
                "value": pheno.lacking_standard,
                "note": "which breakpoint standard was applied is unknown",
            },
            {
                "field": "rows lacking ast_edition",
                "value": pheno.lacking_edition,
                "note": "which edition of that standard was applied is unknown",
            },
        ]
        lines.append(
            "### Phenotype testing-standard provenance\n\n"
            + markdown_table(["field", "value", "note"], rows)
        )
        if not pheno.standard_recorded_everywhere:
            lines.append(
                f"\n**{NO_STANDARD_RECORDED}.** These three fields "
                f"(`ast_method`, `ast_standard`, `ast_edition`) are dropped at "
                f"the stage-11 write: `PhenotypeCall.to_row` emits thirteen "
                f"keys and `run.py` hands `write_tsv` eight, and "
                f"`papipeline/io/tsv.py` discards the rest. The report recovers "
                f"them from the in-memory calls, which is the only place they "
                f"still exist.\n"
            )
        else:
            lines.append(
                "\nEvery phenotype record carries a testing standard, so the "
                "calls in this run are comparable with each other.\n"
            )
    lines.append("")
    lines.append("### Stage tables this report was read from\n")
    lines.append(
        markdown_table(
            ["stage", "path", "rows", "state"],
            [
                {
                    "stage": stage,
                    "path": str(table.path),
                    "rows": table.n_rows,
                    "state": "read" if table.present else NOT_PRODUCED,
                }
                for stage, table in sorted(tables.tables.items())
            ],
        )
    )
    return "\n".join(lines)


def provenance_section(context: ReportContext, tables: RunTables) -> str:
    return _provenance_body(context, tables)


#: Section builders, in the order the report renders them.
#:
#: Keyed on a name rather than inlined in `build_markdown` so a test can assert
#: the list is complete: R2 requires every one of these to be present in every
#: mode, and a section dropped from this tuple would silently disappear.
SECTION_BUILDERS: Tuple[Tuple[str, Any], ...] = (
    ("Run summary and per-stage status", run_summary_section),
    ("Cohort overview", cohort_section),
    ("AMR summary", amr_section),
    ("Virulence summary", virulence_section),
    ("oprD locus verdicts", oprd_section),
    ("Variants summary", variants_section),
    ("Pangenome summary", pangenome_section),
    ("Tree and similarity summary", tree_similarity_section),
    ("Associations: GWAS, convergence, co-occurrence", association_section),
    ("Provenance of the tables read", provenance_section),
)

SectionT = Tuple[str, str]


def build_run_sections(context: ReportContext) -> List[SectionT]:
    """Every section, read back from the run's own output tables.

    Reads the tables once, up front, so a report says what is on disk rather
    than what it was handed. A caller that hands over no tables at all still
    gets all ten sections, each stating that the tables were not read - which is
    the honest answer and is what R2 requires of a missing input.
    """
    tables = context.run_tables
    if tables is None:
        tables = load_run_tables(context.config, context.mode)
        context.run_tables = tables
    return [(heading, builder(context, tables)) for heading, builder in SECTION_BUILDERS]


# --------------------------------------------------------------------------
# Figures - reusing the plotting code that already exists
# --------------------------------------------------------------------------

#: One `_HeatmapMember` per isolate, so `figures.heatmap_figure` can be called.
#:
#: `heatmap_figure` reads `m.pilot_id` and a feature list, and it is written for
#: the pilot-100 cohort. Rather than fork it - and with it the rule that no
#: figure is drawn without sufficient data, and the association annotation on
#: every co-occurrence figure - the report supplies the two attributes it needs.
#: The field is named `pilot_id` because that is the attribute the shared
#: function reads, and renaming the function's contract to suit this caller
#: would be a larger change than adapting to it.
@dataclass(frozen=True)
class _HeatmapMember:
    pilot_id: str
    genes: Tuple[str, ...] = ()


def render_report_figures(
    context: ReportContext, tables: RunTables, out_dir: Path
) -> List[Tuple[str, Optional[Path], bool, str]]:
    """Draw the figures this run's tables support, and report the ones it does not.

    Returns one ``(name, path, created, reason)`` per attempt. A figure that was
    not drawn carries the reason it was not drawn, so the report can state it -
    the same rule `papipeline.pilot.figures` enforces itself ("no figure is drawn
    without sufficient data"): a blank plot is worse than no plot.

    **No new plotting stack.** Every figure here is produced by a function that
    already existed, in `papipeline/pilot/figures.py`, on data reshaped into the
    row shape that function documents. `matplotlib` is already a pinned
    dependency (`environment/environment.yml`, `matplotlib-base=3.11.2`) and is
    imported lazily by that module.
    """
    from ..pilot import figures as figs

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    label = context.antibiotic.capitalize()

    results: List[Tuple[str, Optional[Path], bool, str]] = []

    amr = tables["amr"]
    if amr.present:
        gene_counts: Dict[str, int] = {}
        gene_mechanism: Dict[str, str] = {}
        members: Dict[str, set] = {}
        for row in amr.rows:
            gene = str(row.get("gene") or row.get("determinant") or "").strip()
            if not gene:
                continue
            sample = str(row.get("sample_id") or "")
            gene_counts[gene] = gene_counts.get(gene, 0) + 1
            gene_mechanism.setdefault(gene, str(row.get("mechanism") or "unknown"))
            members.setdefault(sample, set()).add(gene)
        results.append(
            _figure(
                figs.gene_frequency_figure(
                    [
                        {
                            "Gene": gene,
                            "n_isolates_detected": count,
                            "Mechanism": gene_mechanism.get(gene, "unknown"),
                        }
                        for gene, count in gene_counts.items()
                    ],
                    f"{label} AMR determinant detection frequency",
                    out_dir / f"{context.mode.value.lower()}_amr_genes.png",
                ),
                f"the AMR table at {amr.path}",
            )
        )
        results.append(
            _figure(
                figs.heatmap_figure(
                    [_HeatmapMember(s, tuple(sorted(g))) for s, g in sorted(members.items())],
                    lambda m: list(m.genes),
                    f"{label} AMR gene carriage",
                    out_dir / f"{context.mode.value.lower()}_amr_heatmap.png",
                    row_label="Gene",
                ),
                f"the AMR table at {amr.path}",
            )
        )

    phenotype = tables["phenotype"]
    mechanisms = tables["mechanisms"]
    if mechanisms.present and phenotype.present:
        by_sample = {
            str(r.get("sample_id")): str(r.get("phenotype"))
            for r in phenotype.rows
        }
        grouped: Dict[str, Dict[str, int]] = {}
        for row in mechanisms.rows:
            mechanism = str(row.get("mechanism") or "unknown")
            bucket = grouped.setdefault(mechanism, {})
            call = by_sample.get(str(row.get("sample_id")))
            if call:
                bucket[call] = bucket.get(call, 0) + 1
        results.append(
            _figure(
                figs.mechanism_frequency_figure(
                    [
                        {"Mechanism": mechanism, **{f"n_{k}": v for k, v in calls.items()}}
                        for mechanism, calls in sorted(grouped.items())
                    ],
                    f"{label} resistance mechanism by susceptibility category",
                    out_dir / f"{context.mode.value.lower()}_mechanisms.png",
                ),
                f"the mechanism table at {mechanisms.path} with the phenotype "
                f"table at {phenotype.path}",
            )
        )

    cooccurrence = tables["cooccurrence"]
    if cooccurrence.present:
        results.append(
            _figure(
                figs.gene_pairs_figure(
                    [
                        {
                            "Gene_A": str(r.get("feature_a")),
                            "Gene_B": str(r.get("feature_b")),
                            "n_isolates_with_both": r.get("n_both"),
                        }
                        for r in cooccurrence.rows
                    ],
                    f"{label} most frequently co-occurring features",
                    out_dir / f"{context.mode.value.lower()}_cooccurrence_pairs.png",
                ),
                f"the co-occurrence table at {cooccurrence.path}",
            )
        )

    for name, path, created, reason in results:
        LOGGER.info(
            "reporting: figure %s %s%s",
            name,
            "written" if created else "not drawn",
            f" ({reason})" if reason else "",
        )
    return results


def _figure(result: Any, source: str) -> Tuple[str, Optional[Path], bool, str]:
    """Normalise a `figures.FigureResult`, naming the table it needed.

    The reason column carries two things, in this order: what the plotting
    module said, and which table the data would have come from. The second is
    what makes a skipped figure actionable - "only 2 distinct gene(s)" says the
    cohort is small, and "from <path>" says which file to go and look at.
    """
    created = bool(getattr(result, "created", False))
    said = str(getattr(result, "reason", "") or "").strip()
    if created:
        return (
            str(getattr(result, "name", "?")),
            getattr(result, "path", None),
            True,
            f"drawn from {source}",
        )
    reason = said or "the plotting module gave no reason"
    return (
        str(getattr(result, "name", "?")),
        None,
        False,
        f"{reason} - not drawn; the table it would have used is {source}",
    )


def figures_section(results: Sequence[Tuple[str, Optional[Path], bool, str]]) -> str:
    """The figures table, including every figure that was not drawn and why."""
    rows = [
        {
            "figure": name,
            "created": "yes" if created else "no",
            "path": str(path) if path else "-",
            "reason": reason,
        }
        for name, path, created, reason in results
    ]
    body = [
        "Figures are drawn by `papipeline.pilot.figures`, the plotting code this "
        "repository already had, on the run's own tables. That module's rule is "
        "enforced rather than assumed: **no figure is drawn without sufficient "
        "data**, and a figure that was skipped says below why. A blank plot is "
        "worse than no plot.\n",
        markdown_table(["figure", "created", "path", "reason"], rows),
    ]
    if not rows:
        body.append(
            f"\n**{NOT_PRODUCED}: no figure was attempted**, because none of the "
            f"tables this report reads was present.\n"
        )
    else:
        body.append(
            f"\nEvery figure above describes this cohort and nothing more. "
            f"Ruling R16 applies to them as it does to every number here: a "
            f"figure carries the `{POWER_FLAG_TEMPLATE}` flag with the `n` stated "
            f"in the section it summarises.\n"
        )
    return "\n".join(body)


def build_markdown(context: ReportContext) -> str:
    """Render the full Markdown report."""
    mode = context.mode
    banner = banner_for(mode)

    parts: List[str] = []
    title = title_for(mode)

    parts.append(f"# {title}\n")
    if is_smoke_report(context.config):
        # Attached here and verified in `write_report`. The count is the one the
        # run actually used, so the marker cannot claim a number the report does
        # not support.
        parts.append(
            f"> **{smoke_marker(context.n_samples)}**\n>\n"
            f"> {banner}\n"
        )
    else:
        parts.append(f"> **WARNING**\n>\n> {banner}\n")
    parts.append("## Run metadata\n")
    parts.append(
        _kv_table(
            {
                "pipeline_version": __version__,
                "run_mode": mode.value,
                "generated_at_utc": context.generated_at,
                "organism": context.config.organism.get("name"),
                "antibiotic": context.antibiotic,
                "n_samples": context.n_samples,
                "config_root": str(context.config.root),
            }
        )
    )

    parts.append("## Scientific rules in force\n")
    parts.append(
        "This report is generated under the rules in `docs/scientific_rules.md`.\n\n"
        "* Gene presence is never reported as phenotypic resistance.\n"
        "* Association is never reported as causation. The word 'causal' is not "
        "used as a finding anywhere in this pipeline.\n"
        "* R/I/S categories are never converted into MIC or zone diameter.\n"
        "* Candidate structural variants are never promoted to confirmed.\n"
        "* Database versions are always recorded; unpinned ones are reported below.\n"
        "* Missing values remain missing.\n"
    )

    if context.warnings:
        parts.append("## Pipeline health warnings\n")
        for warning in context.warnings:
            parts.append(f"* {warning}\n")
        parts.append("\n")

    # Before the caller's own sections: this is provenance about the cohort the
    # rest of the report is about, so it belongs above the findings rather than
    # in an appendix a reader reaches only after trusting the finding.
    pheno = phenotype_report(context)
    if pheno is not None:
        heading, body = phenotype_provenance_section(context.config, pheno)
        parts.append(f"## {heading}\n")
        parts.append(body.rstrip() + "\n\n")
    elif context.mode is RunMode.REAL:
        # An absence in REAL is a gap in the run, not a stylistic choice, and
        # passing over it silently is what let a report state an association
        # without saying anything about the phenotype it tested.
        parts.append("## Phenotype provenance\n")
        parts.append(
            "No stage 11 phenotype calls were handed to this report, so it "
            "states nothing about the cohort's susceptibility calls, their "
            "categories or their AST provenance. This is an absence in the run, "
            "not a finding.\n\n"
        )
    else:
        # Not REAL, and not silent either. The REAL wording above is a claim
        # about a real run and would be a false alarm on a STUB report, so the
        # non-REAL absence is stated in its own words - and in words that do not
        # reuse the REAL sentence, so a reader (or a test) cannot confuse the
        # two. It is stated at all because a section that appears only when there
        # is something to say must say so when there is not (R2: never a
        # silently absent section).
        parts.append("## Phenotype provenance\n")
        parts.append(
            f"**{NOT_PRODUCED}.** Stage 11 handed this {context.mode.value} "
            f"report no phenotype calls, so it says nothing about susceptibility "
            f"calls, their categories or their AST provenance. This "
            f"{context.mode.value} run does not analyse real isolates, so there "
            f"is no finding here to over- or under-claim.\n\n"
        )

    # Read back from disk, before the caller's own sections: a reader who has
    # just been told what mode this run was in needs the per-stage status and
    # the cohort it is about before any aggregate the orchestrator supplied.
    for heading, body in build_run_sections(context):
        parts.append(f"## {heading}\n")
        parts.append(body.rstrip() + "\n\n")

    for heading, body in context.sections:
        parts.append(f"## {heading}\n")
        parts.append(body.rstrip() + "\n\n")

    for caption, columns, rows in context.tables:
        parts.append(f"### {caption}\n")
        parts.append(markdown_table(columns, rows))

    if context.provenance:
        parts.append("## Provenance\n")
        parts.append(
            "Tool and database versions recorded for this run. An unpinned "
            "version means the reference must be pinned before any real "
            "analysis.\n\n"
        )
        parts.append(
            markdown_table(
                ["reference_id", "tool", "tool_version", "database", "database_version", "version_status"],
                context.provenance,
            )
        )

    parts.append(
        "\n---\n\nGenerated by `papipeline.stages.reporting`. "
        "No scientific conclusion about Pseudomonas aeruginosa imipenem "
        "susceptibility is asserted by this report.\n"
    )
    return "".join(parts)


def markdown_to_html(markdown_text: str) -> str:
    """Minimal, dependency-free Markdown subset renderer.

    Handles the constructs this report actually emits: headings, blockquote
    warnings, bullet lists, pipe tables, bold text and rules. It is not a
    general Markdown implementation.
    """
    out: List[str] = []
    lines = markdown_text.splitlines()
    index = 0

    while index < len(lines):
        line = lines[index]

        if line.startswith("#"):
            level = len(line) - len(line.lstrip("#"))
            out.append(f"<h{level}>{html_module.escape(line[level:].strip())}</h{level}>")
            index += 1
        elif line.startswith("|"):
            table_lines: List[str] = []
            while index < len(lines) and lines[index].startswith("|"):
                table_lines.append(lines[index])
                index += 1
            out.append(_render_html_table(table_lines))
        elif line.startswith("> "):
            quote_lines = []
            while index < len(lines) and lines[index].startswith(">"):
                quote_lines.append(lines[index].lstrip("> ").strip())
                index += 1
            inner = "<br>".join(html_module.escape(q) for q in quote_lines if q)
            out.append(f'<div class="warning"><strong>WARNING</strong><br>{inner}</div>')
        elif line.startswith("* "):
            items = []
            while index < len(lines) and lines[index].startswith("* "):
                items.append(f"<li>{_inline_html(lines[index][2:])}</li>")
                index += 1
            out.append("<ul>" + "".join(items) + "</ul>")
        elif line.strip() == "---":
            out.append("<hr>")
            index += 1
        elif not line.strip():
            index += 1
        else:
            out.append(f"<p>{_inline_html(line)}</p>")
            index += 1

    return (
        "<!DOCTYPE html><html><head><meta charset='utf-8'>"
        f"<title>{html_module.escape(_html_title(markdown_text))}</title>"
        "<style>"
        "body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;"
        "max-width:1100px;margin:2rem auto;padding:0 1rem;line-height:1.5}"
        "table{border-collapse:collapse;margin:1rem 0;font-size:.9rem}"
        "th,td{border:1px solid #ccc;padding:.35rem .6rem;text-align:left}"
        "th{background:#f4f4f4}"
        ".warning{border:2px solid #b00;background:#fff4f4;padding:.75rem;"
        "margin:1rem 0;border-radius:4px}"
        "h1{border-bottom:2px solid #333;padding-bottom:.3rem}"
        "</style></head><body>" + "\n".join(out) + "</body></html>"
    )


def _html_title(markdown_text: str) -> str:
    """The HTML ``<title>``, taken from the document's own H1.

    **This was the static string `Pa Imipenem AMR Pipeline Report`, which names
    no mode.** Everything else about the STUB mislabelling had been fixed - the
    Markdown banner, the Markdown title, the filename - and the browser tab and
    the window title still read the same on a STUB report, a TEST report and the
    real one. That is the same defect one layer out: a reader who has three
    report tabs open cannot tell them apart without opening each one.

    Derived from the H1 rather than from a ``mode`` parameter, so the HTML title
    cannot disagree with the Markdown title by construction: both come from
    ``_TITLES[mode]``, and this reads what the renderer actually emitted. A
    document with no H1 falls back to a mode-neutral string, never to a
    real-sounding one - an untitled report is recoverable, a mislabelled one is
    not.
    """
    for line in markdown_text.splitlines():
        if line.startswith("# "):
            return line[2:].strip() or _UNTITLED_HTML_TITLE
    return _UNTITLED_HTML_TITLE


#: Used only when the rendered document carries no H1 at all.
_UNTITLED_HTML_TITLE = "Pa Imipenem AMR Pipeline Report (untitled document)"


def _render_html_table(table_lines: Sequence[str]) -> str:
    parsed = [
        [cell.strip() for cell in line.strip().strip("|").split("|")]
        for line in table_lines
    ]
    if not parsed:
        return ""
    header = parsed[0]
    body = [r for r in parsed[2:] if not all(set(c) <= {"-", ":", " "} for c in r)]
    head_html = "".join(f"<th>{html_module.escape(c)}</th>" for c in header)
    body_html = "".join(
        "<tr>" + "".join(f"<td>{html_module.escape(c)}</td>" for c in row) + "</tr>"
        for row in body
    )
    return f"<table><thead><tr>{head_html}</tr></thead><tbody>{body_html}</tbody></table>"


def _inline_html(text: str) -> str:
    """Escape text, then restore ``**bold**`` and ``` `code` ``` spans."""
    escaped = html_module.escape(text)
    escaped = _CODE_RE.sub(r"<code>\1</code>", escaped)
    escaped = _BOLD_RE.sub(r"<strong>\1</strong>", escaped)
    return escaped


_CODE_RE = re.compile(r"`([^`]+)`")
_BOLD_RE = re.compile(r"\*\*([^*]+)\*\*")


# --------------------------------------------------------------------------
# REAL-mode refusal on incomplete input
# --------------------------------------------------------------------------

#: The aggregates a REAL report is built from.
#:
#: This list follows from what reporting actually renders, not from what would
#: be convenient to require. `_build_report_context` emits each of these three
#: **conditionally** - `run.py:1615` for gwas rows, `run.py:1653` for
#: convergence, `run.py:1670` for significant co-occurrence pairs - and drops the
#: table without comment when the collection is empty. So each is precisely one
#: whose absence the finished report would otherwise hide: the reader gets a
#: file titled `Pseudomonas aeruginosa Imipenem AMR Report`, bannered
#: `REAL DATASET REPORT`, silently missing whole stages.
#:
#: They are also exactly the three analysis stages `run.PREREQUISITES` already
#: makes reporting depend on (`run.py:213`), so the guard states in the artefact
#: what the DAG already asserts.
#:
#: `provenance` is deliberately **not** here: it is derived from config
#: references (`run.py:656`), not from a stage, so its absence says nothing
#: about whether an analysis ran.
REPORTING_REQUIRED_INPUTS: Tuple[str, ...] = (
    "convergence_calls",
    "cooccurrence",
    "gwas_associations",
)

#: Stands in for "this association exists but was never testable".
#:
#: A stage-12 row whose ``adjusted_p_value`` is ``None`` rendered into the
#: associations table reads as a finding with no p-value. A table where *every*
#: row is like that is a stage claiming to have run and having tested nothing -
#: which is an absence of data wearing the costume of a result. Normalising to a
#: string lets it be declared as a sentinel; `real_inputs` treats ``None`` as
#: "no sentinel declared", so the raw value cannot be.
GWAS_UNTESTED = "no_adjusted_p_value"

#: Per-input sentinel: the value that means "no information" rather than a real
#: label. Every value being the sentinel counts as empty.
_REPORTING_SENTINELS: Mapping[str, str] = {
    # A stage-13 table in which every determinant classifies UNKNOWN is a
    # classification that found nothing - the same failure as an all-`unknown`
    # lineage map, and the reason `is_effectively_empty` exists.
    "convergence_calls": ConvergenceCategory.UNKNOWN.value,
    "gwas_associations": GWAS_UNTESTED,
}

#: What holds each input back, for the refusal message.
_REPORTING_BLOCKERS: Mapping[str, str] = {
    "convergence_calls": (
        "stage 13 (convergence), which refuses in REAL on its own two inputs: "
        "intermediate/regulators/regulator_variants.tsv is read by every stage "
        "but written by none outside the TEST fixture generator, and "
        "lineage_label has no REAL producer at all (only "
        "papipeline.testing.synthetic writes it)"
    ),
    "cooccurrence": (
        "stage 14 (cooccurrence); its mutation namespace reads that same "
        "orphan regulator_variants.tsv path"
    ),
    "gwas_associations": (
        "stage 12 (gwas), which needs stage 6b (pangenome) and so panaroo, "
        "which has no installable osx-arm64 build"
    ),
}


def real_input_paths(config: PipelineConfig, mode: RunMode) -> Dict[str, str]:
    """Where each REAL input's table lives, for the refusal message.

    The locations come from ``contracts`` so this cannot name a stale path -
    the same reason ``convergence.real_input_paths`` does it that way. These are
    where `run.py` writes the tables; reporting itself reads them from memory,
    so the path names where the aggregate *should have been produced*, which is
    what makes the refusal actionable.
    """
    from ..execution.contracts import table_path

    stage_dir = config.intermediate_root(mode) / "stages"
    return {
        "convergence_calls": str(table_path(stage_dir, "convergence")),
        "cooccurrence": str(table_path(stage_dir, "cooccurrence")),
        "gwas_associations": str(table_path(stage_dir, "gwas")),
    }


def declared_input_view(
    declared: Mapping[str, Any]
) -> Dict[str, Dict[str, Any]]:
    """Reshape the declared aggregates into what the shared check can read.

    ``real_inputs.is_effectively_empty`` inspects a **mapping** of values, and
    skips any value that is itself a container - so handing it the raw lists
    would make every input look empty, or (with a sentinel) every input look
    full. ``cooccurrence`` avoids this the same way, by projecting each
    namespace down to its non-empty members before checking
    (`cooccurrence.py:462`).

    So each aggregate is projected to ``{row key: the value worth testing for a
    sentinel}``: the classification for stage 13, the pair key for stage 14, and
    a normalised p-value for stage 12.
    """
    view: Dict[str, Dict[str, Any]] = {}
    for name in REPORTING_REQUIRED_INPUTS:
        rows = list(declared.get(name) or ())
        if name == "convergence_calls":
            view[name] = {
                str(getattr(row, "determinant", index)): getattr(
                    getattr(row, "convergence_category", None), "value", None
                )
                for index, row in enumerate(rows)
            }
        elif name == "cooccurrence":
            # No sentinel: a stage-14 pair list has no "no information" value,
            # only a length. Checked for emptiness alone.
            view[name] = {
                f"{getattr(row, 'feature_a', index)}|"
                f"{getattr(row, 'feature_b', index)}": getattr(
                    row, "feature_type", None
                )
                for index, row in enumerate(rows)
            }
        else:
            view[name] = {
                str(row.get("feature", index)): (
                    row.get("adjusted_p_value")
                    if row.get("adjusted_p_value") is not None
                    else GWAS_UNTESTED
                )
                for index, row in enumerate(rows)
                if isinstance(row, Mapping)
            }
    return view


def refuse_if_incomplete(context: "ReportContext") -> None:
    """Refuse a REAL report whose declared inputs carry no usable content.

    Raises before anything reaches disk, so no partial report is left behind -
    a half-written file bannered `REAL DATASET REPORT` is the dangerous
    artefact, not the exception.
    """
    from .real_inputs import missing_required_inputs, refuse_incomplete

    missing = missing_required_inputs(
        declared_input_view(context.declared_inputs),
        sentinels=_REPORTING_SENTINELS,
    )
    if missing:
        raise refuse_incomplete(
            "reporting",
            missing,
            paths=real_input_paths(context.config, RunMode.REAL),
            blocked_by=_REPORTING_BLOCKERS,
        )


def write_report(
    context: ReportContext,
    out_dir: Path,
    write_html: bool = True,
    write_figures: bool = True,
) -> Dict[str, Path]:
    """Write the Markdown (and optionally HTML) report, and its figures.

    The figures are drawn **before** the Markdown is rendered and their results
    are appended to the context, because a report that lists a figure it did not
    draw - or draws one it never mentions - is worse than either. Drawing first
    also means a plotting failure cannot leave a report on disk claiming a figure
    exists.

    Args:
        context: The run's report context.
        out_dir: Where the report and its figures go. Figures land in
            ``out_dir / "figures"``.
        write_html: Also write the HTML rendering.
        write_figures: Draw the figures this run's tables support. ``False``
            still renders a Figures section, stating that none was attempted -
            never a silently absent one (R2).

    Returns:
        ``{"markdown": path, "html": path}`` plus one ``figure.<name>`` entry per
        figure actually drawn. Figures that were skipped are **not** returned:
        a key whose value is a path that does not exist would be a claim the
        caller could not check.

    Raises:
        StageError: In REAL, when a declared aggregate carries no usable
            content. Checked before ``out_dir`` is even created, so a refusal
            leaves no report and no directory behind.
    """
    if context.mode is RunMode.REAL:
        # Before the mkdir, and before the name is chosen: the check must cost
        # nothing on disk, or a refused run leaves the very artefact it exists
        # to prevent.
        refuse_if_incomplete(context)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    name = str(
        (context.config.raw.get("reporting", {}) or {}).get(
            "test_report_name", "test_pipeline_report"
        )
    )
    if context.mode is RunMode.REAL:
        name = "Pseudomonas_aeruginosa_Imipenem_AMR_Report"

    markdown_path = out_dir / f"{name}.md"

    figures: List[Tuple[str, Optional[Path], bool, str]] = []
    if write_figures:
        tables = context.run_tables
        if tables is None:
            tables = load_run_tables(context.config, context.mode)
            context.run_tables = tables
        try:
            figures = render_report_figures(context, tables, out_dir / "figures")
        except Exception as exc:  # pragma: no cover - plotting is best effort
            # A figure that cannot be drawn is reported as not drawn. Letting a
            # matplotlib failure abort the report would trade a missing picture
            # for a missing report, and the report is the artefact that carries
            # the numbers.
            LOGGER.warning("reporting: figures were not drawn (%s)", exc)
            figures = [("<figures>", None, False, f"plotting failed: {exc}")]
    else:
        figures = [
            (
                "<all>",
                None,
                False,
                "no figure was attempted: `write_report(write_figures=False)`",
            )
        ]
    context.sections.append(("Figures", figures_section(figures)))

    markdown_text = build_markdown(context)

    if is_smoke_report(context.config):
        # Verified here rather than assumed. `build_markdown` attaches the
        # marker, so this is normally a formality - which is the point: it is
        # the check that holds when attachment is not. A report whose body
        # somehow lacks the marker must not reach disk, because a smoke report
        # that reads like the full-cohort analysis is the exact thing this
        # overlay exists to prevent.
        #
        # Placed before the write so a failure leaves no file behind: a
        # half-written report carrying no marker is the dangerous artefact,
        # not the exception.
        require_smoke_marker(markdown_text, context.n_samples)

    markdown_path.write_text(markdown_text, encoding="utf-8")

    written: Dict[str, Path] = {"markdown": markdown_path}
    if write_html:
        html_path = out_dir / f"{name}.html"
        html_path.write_text(markdown_to_html(markdown_text), encoding="utf-8")
        written["html"] = html_path
    for name_, path, created, _reason in figures:
        if created and path is not None:
            written[f"figure.{name_}"] = Path(path)

    LOGGER.info("Report written to %s", markdown_path)
    return written
