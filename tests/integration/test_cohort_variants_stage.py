"""A sixteenth stage, cohort-wide, sitting between per-sample calls and their consumers.

Per-isolate variant calling produces calls against PAO1, and those conflate two
different things: a position where *every* isolate differs from the reference is
a species-wide fixed difference, not a polymorphism. Only a multi-isolate
comparison can separate them, so the merge is its own stage rather than a second
phase of `variants` — the per-sample/cohort distinction is load-bearing in this
DAG, and one rule with two input cardinalities would undercut it.

This pins the *taxonomy*, and after the gate-opening commit the behaviour too.
The stage was in `UNBUILT_STAGES` and refused outside STUB; it now runs - reading
the committed per-isolate VCFs in TEST - so what is pinned is that the DAG knows
about it, that it is cohort-wide rather than per-sample, that its position is
right, and that everything consuming cohort-level variant data reaches it through
`variants`. See :class:`TestItRunsRatherThanRefusing`.

Note the naming: `cohort_variants` is deliberately *not* `variants_merge`. The
suffix pattern already in the taxonomy (`cooccurrence`, `structural_variants`)
is what these read like, and a name that described the mechanism rather than the
artefact would date badly once the merge also handles indels and multiallelic
sites.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.execution.contracts import (
    DENSE_PER_SAMPLE_STAGES,
    PER_SAMPLE_STAGES,
    STAGE_TABLES,
    required_columns,
)
from papipeline.io.tsv import read_tsv
from papipeline.run import (
    EXECUTION_ORDER,
    PREREQUISITES,
    STAGE_ORDER,
    UNBUILT_STAGES,
    resolve_prerequisites,
)

NEW_STAGE = "cohort_variants"
PREDECESSOR = "variants"


@pytest.fixture(scope="module")
def completed_run_rows(tmp_path_factory, pipeline_root):
    """One real TEST run, read back as the merge's rows.

    Module-scoped because a full run is not cheap, and the assertions below are
    about the shape of the table rather than about the run that produced it.
    """
    from papipeline.config.loader import RESULTS_ROOT_ENV, load_config
    from papipeline.run import run_pipeline

    science = pipeline_root / "config" / "science.yaml"
    redirect = tmp_path_factory.mktemp("cohort-merge")
    import os

    previous = os.environ.get(RESULTS_ROOT_ENV)
    os.environ[RESULTS_ROOT_ENV] = str(redirect)
    try:
        result = run_pipeline(
            config_path=science, mode="TEST", config=load_config(science)
        )
    finally:
        if previous is None:
            os.environ.pop(RESULTS_ROOT_ENV, None)
        else:
            os.environ[RESULTS_ROOT_ENV] = previous
    return read_tsv(result.outputs[NEW_STAGE])


class TestTheStageExists:
    def test_it_is_in_the_dag(self):
        assert NEW_STAGE in STAGE_ORDER
        assert NEW_STAGE in EXECUTION_ORDER

    def test_the_dag_has_sixteen_stages(self):
        """spec.md D1 amended from fifteen to sixteen."""
        assert len(STAGE_ORDER) == 16
        assert len(EXECUTION_ORDER) == 16

    def test_it_is_unique_in_the_dag(self):
        assert STAGE_ORDER.count(NEW_STAGE) == 1

    def test_it_has_an_output_contract(self):
        assert NEW_STAGE in STAGE_TABLES
        assert STAGE_TABLES[NEW_STAGE][0] == "cohort_variants.tsv"

    def test_it_declares_a_header(self):
        assert required_columns(NEW_STAGE)


class TestItIsCohortWide:
    """The distinction that is the whole reason this is a separate stage."""

    def test_it_is_not_a_per_sample_stage(self):
        assert NEW_STAGE not in PER_SAMPLE_STAGES
        assert NEW_STAGE not in DENSE_PER_SAMPLE_STAGES

    def test_variants_remains_per_sample(self):
        """Folding the merge in would have removed this; it must stay."""
        assert PREDECESSOR in PER_SAMPLE_STAGES

    def test_the_two_differ(self):
        """Guards against the new stage being quietly added to the per-sample set."""
        assert (PREDECESSOR in PER_SAMPLE_STAGES) is not (
            NEW_STAGE in PER_SAMPLE_STAGES
        )


class TestPosition:
    def test_it_comes_after_variants(self):
        assert STAGE_ORDER.index(NEW_STAGE) > STAGE_ORDER.index(PREDECESSOR)

    def test_it_comes_before_phylogeny(self):
        """Recombination masking and phylogeny both want a cohort SNP set."""
        assert STAGE_ORDER.index(NEW_STAGE) < STAGE_ORDER.index("phylogeny")

    def test_it_needs_variants(self):
        assert PREDECESSOR in PREREQUISITES[NEW_STAGE]

    def test_only_reporting_consumes_it(self):
        """Reporting is its sole consumer, and only so the stage actually runs.

        A stage nothing depends on is never scheduled, so a dry run stays green
        and says nothing about it - the orphan the reachability check exists to
        catch. `reporting` is the honest consumer: nothing *analyses* the output
        yet, because the contract is still open, but a run that claims to have
        executed every stage must have executed this one.
        """
        consumers = [s for s, deps in PREREQUISITES.items() if NEW_STAGE in deps]
        assert consumers == ["reporting"]

    def test_no_analysis_consumes_it_yet(self):
        """`gwas` and `phylogeny` must not read it before its contract is fixed.

        Wiring either to a provisional table would settle the output schema by
        accident, from whatever `gwas.run()` happens to read today - which is
        the mistake the PROVISIONAL column set already was, one level down.
        """
        for analysis in ("gwas", "phylogeny", "convergence", "cooccurrence"):
            assert NEW_STAGE not in PREREQUISITES.get(analysis, frozenset())

    def test_selecting_it_pulls_in_the_per_sample_calls(self):
        """`--only cohort_variants` must not run with no calls to merge."""
        needed = resolve_prerequisites({NEW_STAGE})
        assert PREDECESSOR in needed

    def test_it_is_reachable_from_its_own_prerequisites(self):
        """`resolve_prerequisites` walks dependencies, not dependents.

        Asking for the new stage pulls `variants` in. It is not *pulled in* by
        resolving everything else, because nothing depends on it yet - that is
        the deliberate no-consumer state, pinned by the test above. Asserting it
        were reachable from the full stage set would contradict that.
        """
        assert PREDECESSOR in resolve_prerequisites({NEW_STAGE})


class TestItRunsRatherThanRefusing:
    """The gate is open, so this class asserts the stage works, not that it refuses.

    It was `TestItIsUnbuilt`, and every assertion in it was about the refusal:
    that the stage is listed in `UNBUILT_STAGES`, that the refusal names ticket
    14, that the two refusals read differently. The gate-opening commit made all
    three false by design, so the class was rewritten rather than deleted - the
    distinction between `variants` and `cohort_variants` it used to pin is still
    worth stating, just as a statement about what the two stages *are* now.

    What replaces it: the stage produces a real table in TEST, from the
    committed per-isolate VCFs, and REAL is still refused on every machine
    overlay. Mirrors
    `tests/unit/test_variant_stage_branches.py::TestVariantStagesReachableButRealGated`
    so the two files state the same invariant the same way.
    """

    def test_it_is_no_longer_declared_unbuilt(self):
        assert NEW_STAGE not in UNBUILT_STAGES, (
            f"{NEW_STAGE} is refused as unbuilt; the branch and the merge both "
            "exist, so the refusal is now false"
        )

    def test_its_predecessor_is_no_longer_declared_unbuilt(self):
        assert PREDECESSOR not in UNBUILT_STAGES

    def test_it_is_enabled_in_the_science_config(self, config):
        """The other half of the gate, and the one that fails quietly.

        A stage absent from `analysis` is disabled, and a disabled stage makes
        its Snakemake rule execute nothing and exit non-zero. Both stages were
        un-refused before either was configured, which is what broke the first
        attempt at the reconciliation.
        """
        assert config.stage_enabled(NEW_STAGE) is True
        assert config.stage_enabled(PREDECESSOR) is True

    def test_it_produces_a_table_in_test_mode(self, config, tmp_path, pipeline_root):
        """Runs for real, over the committed fixtures, and writes its contract.

        The table is the stage's principal output and the merge is its whole
        purpose, so "runs and wrote a merge" is the minimum that distinguishes
        this from a stage that refused.
        """
        from papipeline.config.loader import RESULTS_ROOT_ENV, load_config
        from papipeline.run import run_pipeline

        monkey = pytest.MonkeyPatch()
        monkey.setenv(RESULTS_ROOT_ENV, str(tmp_path / "test"))
        try:
            result = run_pipeline(
                config_path=pipeline_root / "config" / "science.yaml",
                mode="TEST", config=load_config(
                    pipeline_root / "config" / "science.yaml"
                ),
            )
        finally:
            monkey.undo()

        assert result.stage_status.get(NEW_STAGE) == "completed", (
            f"{NEW_STAGE} did not complete: "
            f"{result.stage_status.get(NEW_STAGE)!r}"
        )
        path = result.outputs.get(NEW_STAGE)
        assert path is not None, f"no {NEW_STAGE} output was recorded"
        header = Path(path).read_text(encoding="utf-8").splitlines()
        assert header[0].split("\t") == list(STAGE_TABLES[NEW_STAGE][1])
        assert len(header) > 1, (
            f"{NEW_STAGE} wrote a header and no rows; the committed fixtures "
            "merge to a non-empty table, so an empty one means the merge ran "
            "and found nothing, which would be a different finding"
        )

    def test_the_merge_is_cohort_wide_not_per_sample(self, completed_run_rows):
        """The distinction that made this a separate stage rather than a phase.

        No `sample_id` column: one row per variable site across the cohort. If a
        future change pivoted the merge wide, this table would start naming
        isolates per column and every downstream consumer would inherit that
        shape - the reason the long form was chosen (pyseer's `--pres` path).
        """
        assert "sample_id" not in STAGE_TABLES[NEW_STAGE][1]
        for row in completed_run_rows:
            assert set(row) == set(STAGE_TABLES[NEW_STAGE][1])

    @pytest.mark.parametrize(
        "overlay", ["laptop.yaml", "bigmachine.yaml", "smoke.yaml"]
    )
    def test_real_execution_is_still_gated_on_every_overlay(self, overlay):
        """`allow_real_mode` is false everywhere, so REAL is still refused.

        Read per-overlay rather than from the config the tests happen to load: a
        single read would pass while another overlay had been flipped, and this
        flag is the only thing between the pipeline and a REAL cohort being
        analysed on a machine that cannot hold it.
        """
        import yaml

        from tests.unit.test_variant_stage_branches import REPO

        data = yaml.safe_load(
            (REPO / "config" / "machines" / overlay).read_text(encoding="utf-8")
        )
        assert data["runtime"]["allow_real_mode"] is False, (
            f"{overlay} has allow_real_mode true. Opening a stage's gate is not "
            "authorisation to run REAL data; that is authorised separately."
        )

    def test_real_mode_is_still_refused_through_the_public_interface(self, config):
        """The behaviour, not just the flag."""
        from papipeline.errors import ModeNotAllowedError
        from papipeline.run import run_pipeline

        with pytest.raises(ModeNotAllowedError) as excinfo:
            run_pipeline(config=config, mode="REAL")
        assert "allow_real_mode" in str(excinfo.value)


class TestEveryStageIsAccountedFor:
    """The invariant the previous taxonomy work established, now at sixteen."""

    @pytest.mark.parametrize("stage", STAGE_ORDER)
    def test_it_has_an_output_contract(self, stage):
        assert stage in STAGE_TABLES, f"{stage} has no declared output contract"

    @pytest.mark.parametrize("stage", STAGE_ORDER)
    def test_its_prerequisites_exist_and_precede_it(self, stage):
        for dependency in PREREQUISITES.get(stage, frozenset()):
            assert dependency in STAGE_ORDER, f"{stage} names unknown {dependency}"
            assert STAGE_ORDER.index(dependency) < STAGE_ORDER.index(stage), (
                f"{stage} depends on {dependency}, which comes later"
            )

    def test_no_folded_stage_name_leaked_back_in(self):
        """The four modules D1 folded away must stay folded."""
        for folded in ("mechanisms", "regulators", "structural_variants", "integration"):
            assert folded not in STAGE_ORDER
