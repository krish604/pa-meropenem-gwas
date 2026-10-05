"""Fabricate a stage's declared outputs without running its body.

STUB mode exists so the whole DAG, the event stream and the dashboard can be
exercised in seconds, on a machine with no bioinformatics tool installed, no
fixture on disk, and no real parser running (spec.md D8). Every rule still has
to write something at the path the next rule reads, or the DAG is not being
exercised at all - it is being short-circuited.

A fabricated output is therefore *shaped*, not *populated*: the stage's
contract-declared header, zero rows. That is deliberate on three counts:

  - the header comes from :data:`~papipeline.execution.contracts.STAGE_TABLES`,
    the same single source the real writer and the output contract validate
    against, so a stub cannot drift from the shape the pipeline promises;
  - zero rows keeps "no data" distinguishable from "data that happens to be
    empty", which a fabricated row would erase;
  - no fabricated values means nothing downstream can mistake a stub for a
    result. There is deliberately no sample_id here, not even a fake one.

What this module is not: a mock of the stages. It never calls into
``papipeline.stages``, so no parser, validator or classifier runs.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Sequence, Tuple

from .execution.contracts import INTERNAL_TABLES, required_columns, table_path
from .io.tsv import write_tsv

#: Internal steps folded into a stage per spec.md:351, as
#: ``owner -> (step names in INTERNAL_TABLES,)``. The owner fabricates them,
#: because the workflow declares the owner's rule as producing them.
FOLDED_INTO: Dict[str, Tuple[str, ...]] = {
    "amr": ("structural_variants",),
    # `regulators` belongs to stage 6, not stage 4. `run.py` writes
    # 06_regulators.tsv inside the `variants` branch, and rule variants now
    # declares it as a second output - so STUB must fabricate it there too.
    # It was absent here, and STUB failed with a MissingOutputException naming a
    # file its own owner was supposed to write.
    "variants": ("regulators",),
    "cooccurrence": ("mechanisms",),
    "reporting": ("master_table",),
}

#: Outputs a stage declares that are not its principal table. Keyed by stage so
#: that a per-stage invocation writes what that stage owns and nothing else -
#: the workflow calls one stage at a time.
EXTRA_TABLES: Dict[str, Tuple[Tuple[str, Tuple[str, ...]], ...]] = {
    "pangenome": (
        ("gene_presence_absence.tsv", ("sample_id", "gene", "present")),
        ("core_genes.tsv", ("gene",)),
        ("accessory_genes.tsv", ("gene",)),
    ),
    "phylogeny": (
        ("10_alignment_summary.tsv", ("metric", "value")),
    ),
}

#: Banner stamped into a fabricated report, so a stub report can never be
#: mistaken for a real one on disk.
STUB_BANNER = (
    "STUB run: every stage output below is fabricated. No genome, tool, "
    "parser or fixture was involved. Nothing here is a scientific result."
)


def fabricate(stage: str, stage_dir: Path) -> Dict[str, Path]:
    """Write the declared outputs for ``stage`` and return them by key.

    Args:
        stage: A stage name present in the output contract.
        stage_dir: Directory the stage's tables live in.

    Returns:
        Mapping of output key to the path written, matching the keys the real
        stage chain records so that downstream consumers - the run manifest,
        the report, the dashboard - cannot tell a stub from a real run by the
        shape of the mapping alone.
    """
    stage_dir = Path(stage_dir)
    stage_dir.mkdir(parents=True, exist_ok=True)

    written: Dict[str, Path] = {}
    written[stage] = write_tsv(
        table_path(stage_dir, stage), [], required_columns(stage)
    )

    for filename, columns in EXTRA_TABLES.get(stage, ()):  # type: ignore[arg-type]
        key = Path(filename).stem
        written[key] = write_tsv(stage_dir / filename, [], columns)

    # The internal steps spec.md:351 folded into this stage still declare an
    # output, and the workflow declares the owner's rule as producing it. A stub
    # run has to write them too, or the rule finishes having produced nothing
    # and Snakemake stops on a missing output.
    for name in FOLDED_INTO.get(stage, ()):
        filename, columns = INTERNAL_TABLES[name]
        written[name] = write_tsv(stage_dir / filename, [], columns)

    return written


def fabricate_report(reports_root: Path, name: str) -> Dict[str, Path]:
    """Write the stub Markdown and HTML report.

    Args:
        reports_root: Directory the report is written to.
        name: Report basename, without extension. Supplied by the caller
            because it is config-derived and mode-dependent - the same name the
            real reporter uses.

    Returns:
        Mapping of ``markdown``/``html`` to the path written.
    """
    reports_root = Path(reports_root)
    reports_root.mkdir(parents=True, exist_ok=True)

    markdown = f"# STUB run\n\n{STUB_BANNER}\n"
    written = {"markdown": reports_root / f"{name}.md"}
    written["markdown"].write_text(markdown, encoding="utf-8")

    html = reports_root / f"{name}.html"
    html.write_text(
        "<!doctype html>\n<html><head><meta charset='utf-8'>"
        f"<title>STUB run</title></head><body><h1>STUB run</h1>"
        f"<p>{STUB_BANNER}</p></body></html>\n",
        encoding="utf-8",
    )
    written["html"] = html
    return written


def stub_report_name(config, mode) -> str:
    """The report basename for ``mode``, from config, never invented.

    Mirrors the naming rule in :func:`papipeline.stages.reporting.write_report`
    so the workflow's declared output and the reporter agree without either
    duplicating the other.
    """
    if getattr(mode, "value", mode) == "REAL":
        return "Pseudomonas_aeruginosa_Imipenem_AMR_Report"
    configured = (config.raw.get("reporting", {}) or {}).get(
        "test_report_name", "test_pipeline_report"
    )
    return str(configured)


__all__ = [
    "EXTRA_TABLES",
    "STUB_BANNER",
    "fabricate",
    "fabricate_report",
    "stub_report_name",
]
