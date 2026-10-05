"""The report's figures come from the plotting code that already exists, and every skip says why.

**What this pins.** Ruling R4 for this stage: *use the existing visualization
code for figures where it already exists - do not write a new plotting stack.*
The code that already exists is `papipeline/pilot/figures.py`, and it carries two
rules worth more than anything a fresh plotting layer would have:

* **no figure is drawn without sufficient data** - a skipped figure returns the
  reason, because a blank plot is worse than no plot;
* every co-occurrence figure carries `ASSOCIATION_NOTE`, so a bar chart of pairs
  cannot be read as interaction.

Those rules are only worth anything if the report uses them, so this file
asserts the report calls them - by monkeypatching the four entry points and
checking they were reached - as well as the end-to-end behaviour: a cohort with
enough rows produces PNGs, and a cohort with too few produces a table of skips
that name their reasons rather than blank sections.

Synthetic inputs throughout (ruling R12). The tables are the pipeline's own
contracted files written by `write_tsv`-shaped helpers; only the data is fake.
"""

from __future__ import annotations

from pathlib import Path
from typing import List

import pytest

from papipeline.config.loader import PipelineConfig
from papipeline.models import RunMode
from papipeline.stages import reporting

from report_table_fixtures import (  # noqa: F401  (stage_dir is a fixture)
    context,
    populated_tables,
    stage_dir,
)


#: Enough isolates, genes, mechanisms and pairs for the shared plotting module's
#: own sufficiency rules to be satisfied: `MIN_ROWS = 3` distinct genes,
#: `MIN_CATEGORIES = 2` mechanisms, `MIN_ROWS = 3` pairs, and a heatmap needing
#: `MIN_CELLS = 10` populated cells (`papipeline/pilot/figures.py`).
RICH_N = 6
RICH_GENES = ("oprD", "blaKPC-2", "blaNDM-1", "mexR", "mexZ", "nalC", "gyrA", "parC")
RICH_MECHANISMS = ("loss_of_function", "enzymatic_inactivation", "target_mutation")


def rich_tables(config: PipelineConfig, root: Path) -> None:
    """A synthetic cohort large enough that every figure is drawable.

    Deliberately written so that the *shared* module's insufficiency rules all
    pass. If it drew nothing, the test would be measuring the fixture rather
    than the report - and a report that never draws anything cannot be
    distinguished from one that never tries.
    """
    from papipeline.stages.report_tables import read_table  # noqa: F401
    from papipeline.execution.contracts import INTERNAL_TABLES, internal_table_path

    from report_table_fixtures import write_contract_table, write_table

    def write_internal(root: Path, key: str, rows) -> Path:
        return write_table(
            internal_table_path(root, key), rows, INTERNAL_TABLES[key][1]
        )

    amr_rows = []
    for index in range(RICH_N):
        for offset, gene in enumerate(RICH_GENES):
            if (index + offset) % 3:
                amr_rows.append(
                    {
                        "sample_id": f"TEST_PA_{index + 1:03d}",
                        "antibiotic": "imipenem",
                        "determinant": gene,
                        "gene": gene,
                        "mechanism": RICH_MECHANISMS[offset % len(RICH_MECHANISMS)],
                        "evidence_source": "amrfinderplus",
                        "database": "AMRFinderPlus",
                        "database_version": "2026-08-07.1",
                        "claim_status": "DETECTED",
                    }
                )
    write_contract_table(root, "amr", amr_rows)
    write_internal(
        root,
        "mechanisms",
        [
            {
                "sample_id": f"TEST_PA_{index + 1:03d}",
                "antibiotic": "imipenem",
                "mechanism": RICH_MECHANISMS[index % len(RICH_MECHANISMS)],
            }
            for index in range(RICH_N)
        ],
    )
    write_contract_table(
        root,
        "phenotype",
        [
            {
                "sample_id": f"TEST_PA_{index + 1:03d}",
                "phenotype": "R" if index % 2 else "S",
            }
            for index in range(RICH_N)
        ],
    )
    write_contract_table(
        root,
        "cooccurrence",
        [
            {
                "feature_a": RICH_GENES[index],
                "feature_b": RICH_GENES[index + 1],
                "feature_type": "gene",
                "n_both": 2 + index,
                "adjusted_p_value": 0.5,
            }
            for index in range(5)
        ],
    )


def render_and_write(config: PipelineConfig, root: Path, out: Path) -> str:
    ctx = context(config, RunMode.TEST)
    ctx.run_tables = reporting.load_run_tables(config, RunMode.TEST, stage_dir=root)
    reporting.write_report(ctx, out, write_html=False, write_figures=True)
    return (out / "test_pipeline_report.md").read_text(encoding="utf-8")


class TestTheFiguresSectionIsAlwaysPresent:
    def test_it_is_rendered_even_with_no_figures(
        self, config: PipelineConfig, stage_dir: Path, tmp_path: Path
    ):
        """R2 applies to figures too: a skipped figure is named, never absent."""
        text = render_and_write(config, stage_dir, tmp_path / "out")
        assert "## Figures" in text
        assert "no figure was attempted" in text

    def test_it_is_rendered_when_figure_writing_is_switched_off(
        self, config: PipelineConfig, stage_dir: Path, tmp_path: Path
    ):
        populated_tables(config, stage_dir)
        ctx = context(config, RunMode.TEST)
        ctx.run_tables = reporting.load_run_tables(
            config, RunMode.TEST, stage_dir=stage_dir
        )
        reporting.write_report(
            ctx, tmp_path / "off", write_html=False, write_figures=False
        )
        text = (tmp_path / "off" / "test_pipeline_report.md").read_text(encoding="utf-8")
        assert "## Figures" in text
        assert "write_figures=False" in text


class TestTheExistingPlottingCodeIsUsed:
    """R4: no new plotting stack. Asserted by interception, not by inspection."""

    def test_every_figure_is_produced_by_pilot_figures(
        self, config: PipelineConfig, stage_dir: Path, tmp_path: Path, monkeypatch
    ):
        populated_tables(config, stage_dir)
        from papipeline.pilot import figures as figs

        called: List[str] = []
        for name in (
            "gene_frequency_figure",
            "mechanism_frequency_figure",
            "gene_pairs_figure",
            "heatmap_figure",
        ):
            original = getattr(figs, name)

            def spy(*args, _n=name, _f=original, **kwargs):
                called.append(_n)
                return _f(*args, **kwargs)

            monkeypatch.setattr(figs, name, spy)

        ctx = context(config, RunMode.TEST)
        ctx.run_tables = reporting.load_run_tables(
            config, RunMode.TEST, stage_dir=stage_dir
        )
        results = reporting.render_report_figures(
            ctx, ctx.run_tables, tmp_path / "figs"
        )
        assert called, (
            "render_report_figures drew nothing with papipeline.pilot.figures. "
            "This stage is required to reuse the plotting code that already "
            "exists rather than write a new one."
        )
        assert {name for name, _p, _c, _r in results} == {
            "test_amr_genes",
            "test_amr_heatmap",
            "test_mechanisms",
            "test_cooccurrence_pairs",
        }

    def test_a_heatmap_member_carries_the_two_attributes_the_shared_function_reads(
        self, config: PipelineConfig, stage_dir: Path, tmp_path: Path
    ):
        """`figures.heatmap_figure` reads `m.pilot_id`; the report supplies it.

        Pinned as an attribute-name assertion rather than a comment, because the
        field is called `pilot_id` on a non-pilot cohort and that is exactly the
        kind of thing a later tidy-up renames - which would break the shared
        function silently.
        """
        populated_tables(config, stage_dir)
        from papipeline.pilot import figures as figs

        member = reporting._HeatmapMember("TEST_PA_001", ("oprD",))
        assert member.pilot_id == "TEST_PA_001"
        assert member.genes == ("oprD",)
        # The shared function accepts it without any adaptation beyond the two
        # attributes, which is the point of not forking it.
        result = figs.heatmap_figure(
            [member], lambda m: list(m.genes), "t", tmp_path / "h.png"
        )
        assert result.created is False, (
            "one member cannot make a heatmap; the module's own insufficient-data "
            "rule should have refused it"
        )
        assert result.reason


class TestFiguresAreDrawnOnlyWhenThereIsEnoughData:
    def test_a_populated_cohort_draws_its_figures(
        self, config: PipelineConfig, stage_dir: Path, tmp_path: Path
    ):
        rich_tables(config, stage_dir)
        ctx = context(config, RunMode.TEST)
        ctx.run_tables = reporting.load_run_tables(
            config, RunMode.TEST, stage_dir=stage_dir
        )
        results = reporting.render_report_figures(ctx, ctx.run_tables, tmp_path / "f")
        drawn = {n: p for n, p, created, _r in results if created}
        assert drawn, "a cohort with three genes and three isolates drew no figure"
        for name, path in drawn.items():
            assert Path(path).is_file() and Path(path).stat().st_size > 0, (
                f"{name} was reported as drawn but {path} is not a readable file"
            )

    def test_a_skipped_figure_names_the_table_it_would_have_used(
        self, config: PipelineConfig, stage_dir: Path, tmp_path: Path
    ):
        """The populated fixture's heatmap is refused by the shared module's own
        sparse-cell rule, which is the case worth pinning: the skip must carry a
        reason a reader can act on.
        """
        populated_tables(config, stage_dir)
        ctx = context(config, RunMode.TEST)
        ctx.run_tables = reporting.load_run_tables(
            config, RunMode.TEST, stage_dir=stage_dir
        )
        results = reporting.render_report_figures(ctx, ctx.run_tables, tmp_path / "f")
        skipped = {n: r for n, _p, created, r in results if not created}
        assert skipped, "the sparse fixture should have been refused somewhere"
        for name, reason in skipped.items():
            assert "not drawn" in reason, f"{name} was skipped with no reason"
            assert ".tsv" in reason, (
                f"{name}'s skip reason does not name the table it would have used"
            )

    def test_the_report_table_matches_what_was_drawn(
        self, config: PipelineConfig, stage_dir: Path, tmp_path: Path
    ):
        """A report that lists a figure it did not draw is worse than either.

        Asserted by comparing the returned paths against the table, because the
        two are produced in the same call and could still disagree - the table is
        rendered from `figures_section(results)` and the return value is built
        from the same list, so a filter introduced in one and not the other would
        show up here.
        """
        rich_tables(config, stage_dir)
        out = tmp_path / "out"
        text = render_and_write(config, stage_dir, out)
        body = text[text.index("## Figures") :]
        assert "| figure | created | path | reason |" in body
        for name, path, created, _reason in reporting.render_report_figures(
            context(config, RunMode.TEST),
            reporting.load_run_tables(config, RunMode.TEST, stage_dir=stage_dir),
            out / "figures",
        ):
            row = next(
                line for line in body.splitlines() if line.startswith(f"| {name} |")
            )
            if created:
                assert Path(path).name in row, (
                    f"{name} was drawn to {path} but the report does not say so"
                )
            else:
                assert row.split("|")[2].strip() == "no", (
                    f"{name} was not drawn but the report claims it was"
                )

    def test_only_drawn_figures_are_returned_as_paths(
        self, config: PipelineConfig, stage_dir: Path, tmp_path: Path
    ):
        """A returned path that does not exist is a claim the caller cannot check."""
        rich_tables(config, stage_dir)
        written = reporting.write_report(
            context(config, RunMode.TEST), tmp_path / "out", write_html=False
        )
        for key, path in written.items():
            if key.startswith("figure."):
                assert Path(path).is_file()