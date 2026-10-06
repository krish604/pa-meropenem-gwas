from __future__ import annotations

from pathlib import Path
from typing import Dict, Sequence

import pytest

from papipeline.config.loader import PipelineConfig
from papipeline.execution.contracts import (
    INTERNAL_TABLES,
    internal_table_path,
    required_columns,
    table_path,
)
from papipeline.models import RunMode
from papipeline.stages import reporting
from papipeline.stages.report_tables import (
    NOT_PRODUCED,
    SNAKEFILE_ONLY_TABLES,
    STAGE_TABLE_KEYS,
    load_run_tables,
    read_table,
)

# Imported by name: `tests/unit/` is not a package, so pytest puts it on
# `sys.path`. `stage_dir` is re-exported so it is collected as this module's own
# fixture, and the rest is the one copy of these builders.
from report_table_fixtures import (  # noqa: F401
    COHORT_N,
    REQUIRED_SECTIONS,
    amr_row,
    context,
    header_only,
    populated_tables,
    real_guard_payload,
    render,
    stage_dir,
    variant_row,
    write_contract_table,
    write_table,
)


# --------------------------------------------------------------------------
# R1 - the paths are the ones the repository declares
# --------------------------------------------------------------------------


class TestEveryPathComesFromADeclaration:
    def test_each_stage_table_is_read_from_its_contracted_path(
        self, config: PipelineConfig, stage_dir: Path
    ):
        tables = load_run_tables(config, RunMode.REAL, stage_dir=stage_dir)
        for stage, key in STAGE_TABLE_KEYS:
            assert tables[stage].path == table_path(stage_dir, key), (
                f"stage {stage} was read from {tables[stage].path}, which is not "
                f"the path contracts.table_path declares for {key}"
            )

    def test_each_folded_step_is_read_from_its_contracted_path(
        self, config: PipelineConfig, stage_dir: Path
    ):
        tables = load_run_tables(config, RunMode.REAL, stage_dir=stage_dir)
        for key in INTERNAL_TABLES:
            assert tables[key].path == internal_table_path(stage_dir, key)

    def test_the_snakefile_only_tables_are_the_four_the_snakefile_declares(self):
        """Three pan-genome gene tables and the alignment summary, no others.

        Named explicitly so a fifth "we found this one too" entry cannot join
        the list without a decision about where it is declared.
        """
        assert [name for name, _f, _d in SNAKEFILE_ONLY_TABLES] == [
            "gene_presence_absence",
            "core_genes",
            "accessory_genes",
            "alignment_summary",
        ]

    def test_every_stage_in_stage_order_has_a_table_entry(self):
        """A dispatched stage with no table here would be a missing section.

        Asserted against `run.STAGE_ORDER`, the authoritative list, rather than
        against `STAGE_TABLE_KEYS` - a list checked against itself proves nothing.
        """
        from papipeline.run import STAGE_ORDER

        assert set(STAGE_ORDER) == {stage for stage, _ in STAGE_TABLE_KEYS}, (
            "run.STAGE_ORDER and report_tables.STAGE_TABLE_KEYS have drifted; a "
            "stage in one and not the other has no report section"
        )

    def test_the_units_sidecar_is_read_from_beside_the_matrix(
        self, config: PipelineConfig, stage_dir: Path
    ):
        from papipeline.stages.similarity import units_sidecar_path

        tables = load_run_tables(config, RunMode.REAL, stage_dir=stage_dir)
        assert tables.similarity_units is not None
        assert tables.similarity_units.path == units_sidecar_path(
            table_path(stage_dir, "similarity")
        )

    def test_a_table_that_violates_its_contract_is_not_produced_with_the_violation(
        self, config: PipelineConfig, stage_dir: Path
    ):
        """A wrong header is a broken run, and the report must name it.

        The alternative - reading the rows anyway - reports a table whose columns
        are not the ones every other stage agreed on, which is the class of
        defect `contracts.py` exists to catch.
        """
        write_table(table_path(stage_dir, "amr"), [{"not_a_column": 1}])
        tables = load_run_tables(config, RunMode.REAL, stage_dir=stage_dir)
        assert not tables["amr"].present
        assert "missing required columns" in tables["amr"].reason
        assert required_columns("amr")

    def test_cohort_variants_without_an_calls_is_refused(
        self, config: PipelineConfig, stage_dir: Path
    ):
        """`an_calls` is part of the contract now, so its absence is a defect.

        Added by fix 3 of `pa-artifacts/round12/FIX-MPILEUP.md`. A table written
        before it existed must not be read as though it were one that reports its
        denominator: without the column, "all ten isolates were read" and "nine
        were and the tenth never was" are byte-identical, which is exactly the
        indistinguishability the column was added to remove.
        """
        write_table(
            table_path(stage_dir, "cohort_variants"),
            [{"chrom": "NC_002516.2", "pos": "1", "ref": "A", "alt": "G",
              "ac": "1", "an": "2", "af": "0.5"}],
        )
        tables = load_run_tables(config, RunMode.REAL, stage_dir=stage_dir)
        assert not tables["cohort_variants"].present
        assert "missing required columns" in tables["cohort_variants"].reason
        assert "an_calls" in required_columns("cohort_variants")


# --------------------------------------------------------------------------
# R2 - a missing table is named, with its reason, in every mode
# --------------------------------------------------------------------------


class TestAMissingTableIsNamedNotBlanked:
    @pytest.mark.parametrize("mode", list(RunMode))
    def test_every_section_is_present_whatever_is_on_disk(self, config, mode, stage_dir):
        """An empty results tree produces every section, not a short report."""
        text = render(config, mode, stage_dir)
        for heading in REQUIRED_SECTIONS:
            assert f"## {heading}" in text, f"{mode.value} report has no {heading!r} section"

    def test_a_missing_table_says_not_produced_and_why(self, config, stage_dir):
        text = render(config, RunMode.REAL, stage_dir)
        assert NOT_PRODUCED in text
        assert "no file at" in text, (
            "the section says a table is missing without saying the file was "
            "not there, so a reader cannot tell a missing file from an "
            "unreadable one"
        )

    def test_no_section_renders_the_empty_tables_placeholder(self, config, stage_dir):
        """`_No rows._` is what `markdown_table` prints for zero rows.

        Rendered, it reads as a measurement of zero. With no table at all the
        honest rendering is the absent sentence, so the placeholder must not
        appear anywhere in a report built from an empty tree.
        """
        text = render(config, RunMode.REAL, stage_dir)
        assert "_No rows._" not in text

    def test_a_header_only_table_is_not_a_result_of_zero(self, config, stage_dir):
        """The distinction that matters most, and the easiest to get wrong.

        A header-only file is what an unwritten table looks like. Every one of
        these contracts requires at least one row (`contracts.stage_spec` adds
        `min_rows(table, 1)`), so a header-only file has already failed its
        contract - and reporting it as "0 determinants" would report a contract
        violation as science.
        """
        header_only(stage_dir, "amr")
        tables = load_run_tables(config, RunMode.REAL, stage_dir=stage_dir)
        assert not tables["amr"].present
        assert "header and no rows" in tables["amr"].reason
        text = reporting.build_markdown(
            context(config, RunMode.REAL, run_tables=tables)
        )
        assert "**not produced.**" in text
        assert "carries a header and no rows" in text

    def test_an_empty_file_is_distinguished_from_a_missing_one(self, stage_dir):
        table_path(stage_dir, "amr").write_text("", encoding="utf-8")
        read = read_table("amr", table_path(stage_dir, "amr"))
        assert not read.present
        assert "empty" in read.reason

    def test_a_stage_with_no_declared_output_says_so_rather_than_guessing(
        self, stage_dir
    ):
        read = read_table("nonesuch", stage_dir / "nonesuch.tsv")
        assert not read.present
        assert "not produced" in read.absent_sentence

    def test_the_run_tables_view_does_not_invent_a_reason_for_a_declared_table(
        self, config: PipelineConfig, stage_dir: Path
    ):
        """`RunTables.__getitem__` uses `present` as its truth value.

        `or`-chaining a `TableRead` against a fallback therefore substituted the
        fallback's reason - "no output is declared for X" - for the real one on
        every declared-but-absent table. It is pinned because the resulting
        message is confidently false.
        """
        tables = load_run_tables(config, RunMode.REAL, stage_dir=stage_dir)
        assert "recombination" in tables.tables, (
            "recombination is dispatched by run.py and must be a declared table"
        )
        assert tables["recombination"].path == table_path(stage_dir, "recombination")
        assert "no output is declared" not in tables["recombination"].reason


# --------------------------------------------------------------------------
# R1 - the sections carry what they are required to carry
# --------------------------------------------------------------------------


class TestTheReportCarriesTheRunTables:
    def test_the_amr_summary_counts_point_mutation_rows_separately(
        self, config, stage_dir
    ):
        populated_tables(config, stage_dir)
        text = render(config, RunMode.REAL, stage_dir)
        assert "point-mutation rows | 1" in text
        assert "intact-gene rows | 2" in text

    def test_the_amr_summary_says_the_variant_column_is_the_discriminator(
        self, config, stage_dir
    ):
        """`determinant_type` is AMRFinderPlus' `Element Type` and is `AMR`
        for a point mutation and an intact gene alike, so it cannot be the test.
        """
        populated_tables(config, stage_dir)
        text = render(config, RunMode.REAL, stage_dir)
        assert reporting.POINT_MUTATION_COLUMN == "variant"
        assert "`variant` column is non-empty" in text
        assert "does not\ndiscriminate" in text or "does not discriminate" in text

    def test_the_amr_summary_refuses_to_call_a_point_mutation_disruptive(
        self, config, stage_dir
    ):
        """AMRX.md §4.3: `POINT` and `POINT_DISRUPT` partition catalogued and
        uncatalogued alleles, and the frameshift `oprD_Y237TerfsTer0` is among
        the uncatalogued ones. The table records neither, so the report must not
        assert disruption from it.
        """
        populated_tables(config, stage_dir)
        text = render(config, RunMode.REAL, stage_dir)
        assert "not evidence that the mutation disrupts the" in text

    def test_the_cohort_overview_states_st_distribution_and_phenotype_counts(
        self, config, stage_dir
    ):
        populated_tables(config, stage_dir)
        text = render(config, RunMode.REAL, stage_dir)
        assert "### Sequence types (stage 3)" in text
        assert "| ST | isolates |" in text
        assert "### Susceptibility categories (stage 11)" in text
        assert "| category | isolates |" in text
        assert "| R | 1 |" in text and "| I | 1 |" in text

    def test_the_phenotype_section_does_not_merge_the_excluded_categories(
        self, config, stage_dir
    ):
        populated_tables(config, stage_dir)
        text = render(config, RunMode.REAL, stage_dir)
        assert "| I | 1 |" in text, (
            "the I count is missing; R7 forbids merging I or SDD into R or S, so "
            "they must be counted in their own right"
        )

    def test_the_variants_summary_splits_snv_from_non_snv(self, config, stage_dir):
        populated_tables(config, stage_dir)
        text = render(config, RunMode.REAL, stage_dir)
        assert "SNV (REF and ALT both one base) | 1" in text
        assert "non-SNV (every other kind) | 3" in text
        assert "insertion | 1" in text and "deletion | 1" in text and "MNV | 1" in text

    def test_the_variant_kind_test_is_the_allele_string_not_the_biology(
        self, config, stage_dir
    ):
        populated_tables(config, stage_dir)
        text = render(config, RunMode.REAL, stage_dir)
        assert "statement about the allele strings, not about biology" in text

    def test_the_pangenome_summary_comes_from_the_declared_table(
        self, config, stage_dir
    ):
        populated_tables(config, stage_dir)
        text = render(config, RunMode.REAL, stage_dir)
        assert "core_gene_families" in text and "4828" in text

    def test_the_tree_summary_is_read_from_the_phylogeny_table(
        self, config, stage_dir
    ):
        populated_tables(config, stage_dir)
        text = render(config, RunMode.REAL, stage_dir)
        assert "### Tree (stage 9)" in text
        assert "has_branch_lengths" in text

    def test_each_association_result_carries_its_own_flag(self, config, stage_dir):
        populated_tables(config, stage_dir)
        text = render(config, RunMode.REAL, stage_dir)
        assert "### Stage 12 - GWAS associations" in text
        assert "### Stage 13 - convergence" in text
        assert "### Stage 14 - co-occurrence" in text
        assert text.count(f"n={COHORT_N}, underpowered, not a finding") >= 1

    def test_provenance_lists_the_tables_it_read_and_the_ast_gaps(
        self, config, stage_dir
    ):
        populated_tables(config, stage_dir)
        ctx = context(config, RunMode.REAL)
        ctx.run_tables = load_run_tables(config, RunMode.REAL, stage_dir=stage_dir)
        ctx.phenotype_calls = []
        text = reporting.build_markdown(ctx)
        assert "### Phenotype testing-standard provenance" in text
        assert reporting.NO_STANDARD_RECORDED in text
        assert "rows lacking ast_standard" in text
        assert "### Stage tables this report was read from" in text

    def test_the_run_summary_lists_every_stage_with_a_status(
        self, config, stage_dir
    ):
        populated_tables(config, stage_dir)
        text = render(config, RunMode.REAL, stage_dir)
        assert "| stage | status | reason |" in text
        for stage, _key in STAGE_TABLE_KEYS:
            assert f"| {stage} |" in text, f"{stage} has no status row"

    def test_a_status_from_the_orchestrator_is_preferred_over_the_disk_fallback(
        self, config, stage_dir
    ):
        populated_tables(config, stage_dir)
        ctx = context(
            config,
            RunMode.REAL,
            stage_status={"amr": "refused", "gwas": "completed"},
            stage_refusals={"amr": "AMRFinderPlus database is not present at db/..."},
            run_tables=load_run_tables(config, RunMode.REAL, stage_dir=stage_dir),
        )
        text = reporting.build_markdown(ctx)
        assert "| amr | refused | AMRFinderPlus database is not present" in text
        assert "the orchestrator's own record" in text

    def test_an_orchestrator_status_for_an_absent_table_names_the_missing_file(
        self, config, stage_dir
    ):
        """The same truthiness trap as the table view, in the status table.

        `TableRead.__bool__` is `present`, so `if table:` substituted "no table
        declared for this stage" for the real reason on every declared-but-absent
        table - and that message is confidently false, naming a missing
        declaration where there is one.
        """
        ctx = context(
            config,
            RunMode.REAL,
            stage_status={"mlst": "completed", "similarity": "completed"},
            run_tables=load_run_tables(config, RunMode.REAL, stage_dir=stage_dir),
        )
        text = reporting.build_markdown(ctx)
        assert "no table declared for this stage" not in text
        assert text.count("no file at") >= 2

    def test_the_disk_fallback_says_where_its_statuses_came_from(
        self, config, stage_dir
    ):
        """A stage refused, one skipped and one that crashed look identical to
        a file-presence test, so the weaker claim must be labelled as one.
        """
        populated_tables(config, stage_dir)
        text = render(config, RunMode.REAL, stage_dir)
        assert "derived from disk" in text
        assert "ReportContext.stage_status" in text


# --------------------------------------------------------------------------
# R1 - the oprD verdict table, with its evidence
# --------------------------------------------------------------------------


class _Evidence:
    lesion_type = "premature_stop"
    internal_stop_position_aa = 237
    truncation_aa = 207


class _StructuralCall:
    def __init__(self, sample_id: str, verdict: str = "disrupted"):
        self.sample_id = sample_id
        self.verdict = verdict
        self.evidence = _Evidence()
        self.identity_pct = 97.96
        self.coverage_pct = 100.0
        self.reason = "premature stop at codon 237"


class _PresenceOnlyResolution:
    """What `resolve_oprd_loci` actually returns today: no reading-frame evidence."""

    def __init__(self, sample_id: str, verdict: str = "resolved"):
        self.sample_id = sample_id
        self.verdict = verdict
        self.identity_pct = 94.58
        self.coverage_pct = 100.0


class TestTheOprdVerdictTable:
    def test_it_carries_verdict_lesion_position_and_truncation(self, config, stage_dir):
        ctx = context(
            config,
            RunMode.REAL,
            oprd_calls=[_StructuralCall("TEST_PA_001")],
            run_tables=load_run_tables(config, RunMode.REAL, stage_dir=stage_dir),
        )
        text = reporting.build_markdown(ctx)
        for column in reporting.OPRD_VERDICT_COLUMNS:
            assert column in text, f"the oprD table has no {column} column"
        assert "disrupted" in text
        assert "premature_stop" in text
        assert "| 237 |" in text
        assert "| 207 |" in text

    def test_it_states_ruling_2_so_disrupted_cannot_be_read_as_short(
        self, config, stage_dir
    ):
        ctx = context(
            config,
            RunMode.REAL,
            oprd_calls=[_StructuralCall("TEST_PA_001")],
            run_tables=load_run_tables(config, RunMode.REAL, stage_dir=stage_dir),
        )
        text = reporting.build_markdown(ctx)
        assert "RULING 2" in text
        assert "never by ORF length" in text

    def test_truncation_is_evidence_not_the_decision(self, config, stage_dir):
        ctx = context(
            config,
            RunMode.REAL,
            oprd_calls=[_StructuralCall("TEST_PA_001")],
            run_tables=load_run_tables(config, RunMode.REAL, stage_dir=stage_dir),
        )
        text = reporting.build_markdown(ctx)
        assert "evidence, not a decision" in text

    def test_a_presence_only_verdict_leaves_the_evidence_empty_and_says_why(
        self, config, stage_dir
    ):
        """A `LocusResolution` judges presence, not the reading frame. Its lesion
        cells are empty because the evidence does not exist, and the report must
        say that rather than let a blank column read as "no lesion found".
        """
        rows = reporting.oprd_verdict_rows([_PresenceOnlyResolution("TEST_PA_001")])
        assert rows[0]["lesion_type"] is None
        assert rows[0]["truncation_aa"] is None
        ctx = context(
            config,
            RunMode.REAL,
            oprd_calls=[_PresenceOnlyResolution("TEST_PA_001")],
            run_tables=load_run_tables(config, RunMode.REAL, stage_dir=stage_dir),
        )
        text = reporting.build_markdown(ctx)
        assert "does not judge the reading frame" in text

    def test_an_absent_locus_search_is_named_because_there_is_no_path_to_read(
        self, config, stage_dir
    ):
        ctx = context(
            config, RunMode.REAL, oprd_calls=None,
            run_tables=load_run_tables(config, RunMode.REAL, stage_dir=stage_dir),
        )
        text = reporting.build_markdown(ctx)
        assert "no contracted path exists for the oprD locus verdict" in text
        assert "STAGE_TABLES" in text and "workflow/Snakefile" in text

    def test_an_empty_search_is_distinguished_from_no_search(self, config, stage_dir):
        ctx = context(
            config, RunMode.REAL, oprd_calls=[],
            run_tables=load_run_tables(config, RunMode.REAL, stage_dir=stage_dir),
        )
        text = reporting.build_markdown(ctx)
        assert "returned no isolate at all" in text

    def test_the_amr_table_saying_nothing_about_oprd_is_not_absence(self, config, stage_dir):
        """The fabricated negative this study already paid for once."""
        write_contract_table(
            stage_dir,
            "amr",
            [amr_row(sample_id="TEST_PA_001", gene="blaKPC-2", determinant="blaKPC-2")],
        )
        text = render(config, RunMode.REAL, stage_dir)
        assert "Absence from the AMR\ntable is **not** evidence" in text or (
            "Absence from the AMR table is **not** evidence" in text
        )


# --------------------------------------------------------------------------
# R1 - similarity states its units
# --------------------------------------------------------------------------


class TestSimilarityStatesItsUnits:
    def test_the_units_come_from_the_sidecar_when_it_is_there(
        self, config, stage_dir
    ):
        populated_tables(config, stage_dir)
        text = render(config, RunMode.REAL, stage_dir)
        assert "#### Units of the distances" in text
        assert "expected substitutions per site" in text
        assert "patristic" in text
        assert "n_pairs" in text

    def test_a_missing_sidecar_is_still_stated_as_a_unit_and_flagged(
        self, config, stage_dir
    ):
        """The unit is a property of the computation, so it is stated even when
        the file recording it is absent - and the absent file is named, because
        `read_units_sidecar` had zero callers until this section existed.
        """
        populated_tables(config, stage_dir)
        from papipeline.stages.similarity import units_sidecar_path

        units_sidecar_path(table_path(stage_dir, "similarity")).unlink()
        text = render(config, RunMode.REAL, stage_dir)
        assert "expected substitutions per site" in text
        assert "units sidecar" in text
        assert NOT_PRODUCED in text

    def test_the_report_says_the_numbers_are_not_snp_counts(self, config, stage_dir):
        """The single most likely misreading of a bare distance matrix."""
        populated_tables(config, stage_dir)
        text = render(config, RunMode.REAL, stage_dir)
        assert "NOT SNP counts" in text

    def test_an_absent_alignment_summary_names_why_it_is_absent(
        self, config, stage_dir
    ):
        populated_tables(config, stage_dir)
        (stage_dir / "10_alignment_summary.tsv").unlink()
        text = render(config, RunMode.REAL, stage_dir)
        assert "alignment summary" in text
        assert "core_alignment.fasta" in text


# --------------------------------------------------------------------------
# R16 - the power flag
# --------------------------------------------------------------------------


class TestThePowerFlag:
    def test_the_literal_r16_names_is_the_flag_at_n_ten(self):
        """R16 names `n=10, underpowered, not a finding` verbatim."""
        assert reporting.R16_FLAG == "n=10, underpowered, not a finding"
        assert reporting.power_flag(COHORT_N) == reporting.R16_FLAG

    def test_the_flag_names_the_cohort_actually_used(self):
        """A flag that hard-codes 10 while the numbers came from 20 states a
        cohort that was not analysed.
        """
        assert reporting.power_flag(20) == "n=20, underpowered, not a finding"

    @pytest.mark.parametrize(
        "heading",
        [
            "Cohort overview",
            "AMR summary",
            "Virulence summary",
            "Variants summary",
            "Pangenome summary",
            "Tree and similarity summary",
            "Associations: GWAS, convergence, co-occurrence",
        ],
    )
    def test_every_statistic_section_carries_it(self, config, stage_dir, heading):
        populated_tables(config, stage_dir)
        text = render(config, RunMode.REAL, stage_dir)
        start = text.index(f"## {heading}")
        rest = text[start:]
        end = rest.find("\n## ", 1)
        body = rest if end < 0 else rest[:end]
        assert "underpowered, not a finding" in body, f"{heading} has no R16 flag"

    def test_the_flag_says_what_it_forbids(self):
        line = reporting.power_flag_line(10)
        assert "not evidence about" in line
        assert "may" in line and "be cited as a finding" in line


# --------------------------------------------------------------------------
# R3 - every artefact is mode-aware
# --------------------------------------------------------------------------


class TestEveryArtefactIsModeAware:
    @pytest.mark.parametrize("mode", list(RunMode))
    def test_the_html_title_names_the_mode(self, config, mode, tmp_path):
        ctx = context(config, mode, **real_guard_payload(mode))
        written = reporting.write_report(
            ctx, tmp_path / mode.value, write_html=True, write_figures=False
        )
        html = written["html"].read_text(encoding="utf-8")
        assert f"<title>{reporting.title_for(mode)}</title>" in html, (
            "the HTML <title> was the static string "
            "'Pa Imipenem AMR Pipeline Report', which names no mode: a reader "
            "with three report tabs open could not tell them apart without "
            "opening each one"
        )

    def test_the_stub_html_title_is_not_the_real_one(self, config, tmp_path):
        written = reporting.write_report(
            context(config, RunMode.STUB), tmp_path, write_html=True, write_figures=False
        )
        html = written["html"].read_text(encoding="utf-8")
        assert reporting.title_for(RunMode.REAL) not in html

    def test_a_document_with_no_heading_gets_a_mode_neutral_title(self):
        assert "Real" not in reporting._html_title("no heading here")
        assert reporting._html_title("# Something\nrest") == "Something"

    def test_an_untitled_document_fails_closed_rather_than_to_the_real_title(self):
        assert reporting._UNTITLED_HTML_TITLE != reporting.title_for(RunMode.REAL)

    @pytest.mark.parametrize("mode", list(RunMode))
    def test_the_phenotype_absence_is_stated_in_every_mode(self, config, mode, tmp_path):
        """Not REAL means the absence wording differs, not that it is dropped.

        R2: a section that appears only when there is something to say must say
        so when there is not.
        """
        text = reporting.build_markdown(context(config, mode))
        assert "## Phenotype provenance" in text
        assert "phenotype calls" in text
        if mode is RunMode.REAL:
            assert "This is an absence in the run" in text
        else:
            assert "not produced" in text