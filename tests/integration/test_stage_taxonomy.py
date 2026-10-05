"""The stage taxonomy is the spec's sixteen, and the unbuilt stages say so.

spec.md D1 fixes the canonical list - sixteen since the 2026-09-29 amendment
that added `cohort_variants` - and says (spec.md:349-353) that
the existing sixteen-stage code is *mapped onto* it, with `mechanisms`,
`regulators`, `structural_variants` and `integration` becoming internal steps
and the ordering constants rewritten.

The property that matters most here is the one about the unbuilt stages.
`recombination` and `similarity` have no implementation. An unbuilt
stage that emitted a header-only table would be indistinguishable from a stage
that ran and found nothing - and "no sample is similar to any other" is a
finding. So they fail loudly outside STUB, and fabricate only inside it.

`variants` and `cohort_variants` were in that group until the gate-opening
commit; they left it, and :data:`UNBUILT` below is what they left behind.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.execution.contracts import INTERNAL_TABLES, STAGE_TABLES
from papipeline.run import PREREQUISITES, STAGE_ORDER, UNBUILT_STAGES, run_pipeline

#: spec.md D1, in order. The code keeps its own spelling for stages that map
#: one-to-one, and adopts the spec's where the code has no equivalent.
SPEC_D1 = (
    "validate", "annotate", "mlst", "amr", "virulence", "variants",
    "cohort_variants",  # 6a, added by the 2026-09-29 amendment
    "pangenome", "recombination", "phylogeny", "similarity", "phenotype",
    "gwas", "convergence", "combination", "report",
)

CODE_NAME_FOR_SPEC = {
    "validate": "validation",
    "annotate": "annotation",
    "combination": "cooccurrence",
    "report": "reporting",
}

#: The four that spec.md:351 demotes to internal steps of another stage.
FOLDED = ("mechanisms", "regulators", "structural_variants", "integration")

#: Derived from the code, not restated. A hard-coded copy keeps asserting a stage
#: is absent after it has a caller - which is exactly what happened to
#: `similarity`, and what `variants` and `cohort_variants` before it.
UNBUILT = tuple(sorted(UNBUILT_STAGES))


class TestTheStageListIsSixteen:
    def test_there_are_sixteen(self):
        assert len(STAGE_ORDER) == 16, (
            f"expected the spec's sixteen, got {len(STAGE_ORDER)}: {STAGE_ORDER}"
        )

    def test_membership_matches_the_spec(self):
        expected = {CODE_NAME_FOR_SPEC.get(s, s) for s in SPEC_D1}
        assert set(STAGE_ORDER) == expected, (
            "the stage set does not match spec.md D1\n"
            f"  extra in code: {sorted(set(STAGE_ORDER) - expected)}\n"
            f"  missing:        {sorted(expected - set(STAGE_ORDER))}"
        )

    @pytest.mark.parametrize("name", FOLDED)
    def test_a_folded_module_is_no_longer_a_stage(self, name):
        assert name not in STAGE_ORDER, (
            f"{name!r} was folded into another stage per spec.md:351 but is "
            "still a stage of its own"
        )
        assert name not in STAGE_TABLES, (
            f"{name!r} is not a stage, so it must not own a stage's principal "
            "output; it belongs in INTERNAL_TABLES"
        )
        assert name in INTERNAL_TABLES or name == "integration", (
            f"{name!r} still writes a file, so its location needs declaring"
        )

    def test_every_stage_has_a_declared_output_table(self):
        missing = [s for s in STAGE_ORDER if s not in STAGE_TABLES]
        assert not missing, f"stages with no output contract: {missing}"

    def test_no_output_is_declared_as_both_a_stage_and_an_internal_table(self):
        overlap = set(STAGE_TABLES) & set(INTERNAL_TABLES)
        assert not overlap, (
            f"{overlap} cannot be both a stage's principal output and an "
            "internal step's; the duplication would let the two drift apart"
        )

    def test_prerequisites_only_name_real_stages(self):
        for stage, needs in PREREQUISITES.items():
            assert stage in STAGE_ORDER, f"PREREQUISITES names unknown stage {stage!r}"
            for need in needs:
                assert need in STAGE_ORDER, (
                    f"stage {stage!r} depends on {need!r}, which is not a stage"
                )

    def test_a_folded_module_is_still_a_real_module(self):
        """The code survives the fold; only the stage entry goes.

        Their content and tests are explicitly preserved by spec.md:351, so
        dropping the modules would be a regression dressed as a rename.
        """
        import importlib

        for module in ("mechanisms", "regulators", "sv", "integration"):
            assert importlib.import_module(f"papipeline.stages.{module}") is not None


class TestUnbuiltStagesAreHonest:
    @pytest.mark.parametrize("stage", UNBUILT)
    def test_declared_in_the_dag(self, stage):
        assert stage in STAGE_ORDER
        assert stage in STAGE_TABLES
        assert stage in UNBUILT_STAGES

    @pytest.mark.parametrize("stage", UNBUILT)
    def test_stub_fabricates_it(self, config, monkeypatch, tmp_path, stage):
        from pathlib import Path

        from papipeline.config.loader import RESULTS_ROOT_ENV

        monkeypatch.setenv(RESULTS_ROOT_ENV, str(tmp_path / "stub"))
        result = run_pipeline(config=config, mode="STUB")
        path = result.outputs.get(stage)
        assert path is not None, f"STUB produced no output for {stage!r}"

        lines = Path(path).read_text(encoding="utf-8").splitlines()
        assert lines[0].split("\t") == list(STAGE_TABLES[stage][1]), (
            f"{stage!r} did not fabricate its contract header"
        )
        assert len(lines) == 1, (
            f"{stage!r} fabricated {len(lines) - 1} rows; a stub must invent no "
            "data, or it becomes indistinguishable from a result"
        )

    @pytest.mark.parametrize("stage", UNBUILT)
    def test_test_mode_refuses_rather_than_faking(self, config, tmp_path, stage):
        from papipeline.config.loader import RESULTS_ROOT_ENV

        monkey = pytest.MonkeyPatch()
        monkey.setenv(RESULTS_ROOT_ENV, str(tmp_path / "refuse"))
        try:
            with pytest.raises(NotImplementedError) as excinfo:
                run_pipeline(config=config, mode="TEST", only=[stage])
        finally:
            monkey.undo()

        message = str(excinfo.value)
        assert "not implemented" in message.lower(), message
        assert "14" in message, (
            f"the refusal should cite the ticket that will build it: {message}"
        )

    @pytest.mark.parametrize("stage", UNBUILT)
    def test_no_output_is_written_when_it_refuses(self, config, monkeypatch, tmp_path, stage):
        """The dangerous outcome is a file, not an exception.

        A header-only table left behind by a refused run would be
        indistinguishable from a stage that ran and found nothing.
        """
        import glob
        from pathlib import Path

        from papipeline.config.loader import RESULTS_ROOT_ENV

        redirected = tmp_path / "refused"
        monkeypatch.setenv(RESULTS_ROOT_ENV, str(redirected))
        with pytest.raises(NotImplementedError):
            run_pipeline(config=config, mode="TEST", only=[stage])

        # Exact filename, not a substring: structural_variants folds into amr and
        # writes 07_structural_variants.tsv, which contains "variants".
        expected = STAGE_TABLES[stage][0]
        leaked = [p for p in glob.glob(str(redirected / "**" / "*.tsv"), recursive=True)
                  if Path(p).name == expected]
        assert not leaked, f"a refused stage left output behind: {leaked}"


class TestSpecEdgesAreDeclared:
    """The edges spec.md requires, asserted one by one.

    Reachability cannot police this: `variants` stays reachable through
    `reporting` even if the `gwas` edge is deleted, so a reachability-only check
    passes on a DAG that has quietly lost a dependency the spec calls required
    (D6: the SNP/indel feature family is one of three GWAS input families).

    So the edges are asserted directly, and deleting one is a red test.
    """

    #: stage -> (name, expected output variable) each must consume.
    SPEC_EDGES = {
        "gwas": ("OUT_VARIANTS",),              # D6: SNP/indel feature family
        "phylogeny": ("OUT_RECOMBINATION",),    # D1: stage 9 consumes stage 8
        "reporting": ("OUT_VARIANTS", "OUT_RECOMBINATION", "OUT_SIMILARITY",
                     "OUT_MECHANISMS", "OUT_MASTER" if False else "OUT_AMR"),
    }

    @staticmethod
    def _rule(snakefile: str, name: str) -> str:
        import re

        match = re.search(
            rf"(?ms)^rule {name}:\n.*?(?=^rule |\Z)", snakefile
        )
        assert match, f"no rule named {name!r}"
        return match.group(0)

    @pytest.mark.parametrize("stage,outputs", sorted(
        (k, v) for k, v in SPEC_EDGES.items()
    ))
    def test_the_edge_is_declared(self, stage, outputs):
        snakefile = (Path(__file__).resolve().parents[2] / "workflow" / "Snakefile").read_text()
        block = self._rule(snakefile, stage)
        for output in outputs:
            assert output in block, (
                f"rule {stage!r} does not declare {output} as an input. The spec "
                "requires the edge; reachability alone will not catch its loss."
            )

    def test_the_unbuilt_stages_are_rules_in_the_dag(self):
        snakefile = (Path(__file__).resolve().parents[2] / "workflow" / "Snakefile").read_text()
        for stage in UNBUILT:
            assert f"rule {stage}:" in snakefile, (
                f"{stage!r} is one of the spec's sixteen but has no rule"
            )

    def test_no_rule_exists_for_a_folded_stage(self):
        """The four are internal steps; a rule for them is the old taxonomy."""
        snakefile = (Path(__file__).resolve().parents[2] / "workflow" / "Snakefile").read_text()
        for folded in ("regulators", "structural_variants", "mechanisms", "integration"):
            assert f"rule {folded}:" not in snakefile, (
                f"{folded!r} was folded into another stage (spec.md:351) but "
                "still has a rule of its own"
            )


# --------------------------------------------------------------------------
# The cascade, asserted rather than checked by hand
# --------------------------------------------------------------------------

class TestTheCascadeThroughTheWorkflow:
    """The DAG consequences of the unbuilt stages, end to end.

    Verified by hand once and green in CI ever since, which is not the same
    thing: a hand check is not re-run when something else moves.
    """

    @staticmethod
    def _snakefile() -> str:
        return (Path(__file__).resolve().parents[2] / "workflow" / "Snakefile").read_text()

    def _dry_run(self, tmp_path, mode):
        import os
        import subprocess

        from papipeline.config.loader import RESULTS_ROOT_ENV

        return subprocess.run(
            ["snakemake", "--snakefile",
             str(Path(__file__).resolve().parents[2] / "workflow" / "Snakefile"),
             "--config", f"mode={mode}", "machine=laptop",
             "--cores", "1", "--dry-run"],
            cwd=Path(__file__).resolve().parents[2],
            env={**os.environ, RESULTS_ROOT_ENV: str(tmp_path / mode)},
            capture_output=True, text=True, timeout=600,
        )

    @pytest.mark.parametrize("mode", ["TEST", "STUB"])
    def test_the_dag_reaches_all_sixteen_stages(self, tmp_path, mode):
        """Every declared stage must be scheduled, or a rule never runs."""
        import re

        result = self._dry_run(tmp_path, mode)
        assert result.returncode == 0, (
            f"{mode} DAG did not resolve:\n{result.stdout[-2000:]}"
        )
        combined = result.stdout + result.stderr
        scheduled = {m.group(1) for m in re.finditer(r"(?m)^\s*rule (\w+):", combined)}
        missing = sorted(set(STAGE_ORDER) - scheduled)
        assert not missing, f"{mode}: these stages are never scheduled: {missing}"

    @pytest.mark.parametrize("stage", UNBUILT)
    def test_test_mode_cannot_run_an_unbuilt_stage(self, tmp_path, stage):
        """Refused, named, and nothing written.

        The last part is the one that matters: a header-only table left behind
        by a refused run would read as a stage that looked and found nothing.
        """
        import os
        import subprocess

        from papipeline.config.loader import RESULTS_ROOT_ENV

        results_root = tmp_path / f"refuse-{stage}"
        target = results_root / "test" / "intermediate" / "stages" / STAGE_TABLES[stage][0]
        result = subprocess.run(
            ["snakemake", "--snakefile",
             str(Path(__file__).resolve().parents[2] / "workflow" / "Snakefile"),
             "--config", "mode=TEST", "machine=laptop",
             f"status_log={tmp_path / 'events.jsonl'}",
             "--cores", "1", str(target)],
            cwd=Path(__file__).resolve().parents[2],
            env={**os.environ, RESULTS_ROOT_ENV: str(results_root)},
            capture_output=True, text=True, timeout=600,
        )
        combined = result.stdout + result.stderr
        assert result.returncode != 0, "an unbuilt stage ran in TEST"
        assert "not implemented" in combined.lower(), combined[-2000:]
        assert "14" in combined, "the refusal must cite the ticket that builds it"
        assert not target.exists(), (
            f"{stage} refused but left {target} behind; an empty table is a "
            "finding and this stage looked for nothing"
        )
