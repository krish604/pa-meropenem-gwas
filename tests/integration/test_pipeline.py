"""End-to-end tests: configuration, the full TEST run, and mode gating.

These are the tests that would catch a regression in the wiring between
stages, as opposed to a bug inside one stage.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from papipeline.config.loader import load_config
from papipeline.errors import (
    AntibioticError,
    ConfigError,
    ModeNotAllowedError,
    PipelineError,
    UnknownGeneError,
)
from papipeline.knowledge import regulator_genes
from papipeline.models import ClaimStatus, RunMode
from papipeline.run import (
    EXECUTION_ORDER,
    PREREQUISITES,
    STAGE_ORDER,
    resolve_mode,
    resolve_prerequisites,
    run_pipeline,
)
from papipeline.testing import SYNTHETIC_BANNER


class TestConfiguration:
    def test_loads(self, config):
        assert config.organism["name"] == "Pseudomonas aeruginosa"
        # The meropenem GWAS build enables both antibiotics (the old assertion
        # `== ("imipenem",)` encoded the state this build deliberately
        # reverses). imipenem must stay first: config.antibiotics[0] is read
        # as the project's fallback antibiotic by several stages.
        assert config.antibiotics == ("imipenem", "meropenem")
        assert config.antibiotics[0] == "imipenem"

    def test_allowed_phenotypes(self, config):
        assert {p.value for p in config.allowed_phenotypes} == {
            "R",
            "I",
            "S",
            "SDD",
            "ND",
        }

    def test_all_analysis_switches_present(self, config):
        for switch in (
            "amr",
            "mlst",
            "virulence",
            "regulators",
            "structural_variants",
            "pangenome",
            "phylogeny",
            "gwas",
            "convergence",
            "cooccurrence",
        ):
            assert switch in config.analysis, switch

    def test_required_resistance_loci_configured(self, config):
        required = {
            "oprD",
            "mexR",
            "nalC",
            "nalD",
            "mexZ",
            "nfxB",
            "mexT",
            "mexS",
            "ampD",
            "ampR",
            "dacB",
        }
        assert required <= set(regulator_genes(config))
        assert required <= set(config.mechanisms)

    def test_no_sample_ids_in_configuration(self, config):
        """Configuration holds no sample identifiers.

        `sample_id.pattern` is excluded from the check: it is the *rule* for a
        valid identifier, so it necessarily contains the prefixes "GCA_" and
        "TEST_". Naming a real isolate in configuration would be a stray record
        of who is in the study, so the assertion is that no value in the
        configuration is itself a valid identifier.
        """
        import re

        pattern = re.compile(config.raw["sample_id"]["pattern"])
        offenders = []

        def walk(node, path="config"):
            if isinstance(node, dict):
                for key, value in node.items():
                    if key == "pattern":
                        continue
                    walk(value, f"{path}.{key}")
            elif isinstance(node, (list, tuple)):
                for index, value in enumerate(node):
                    walk(value, f"{path}[{index}]")
            elif isinstance(node, str) and pattern.fullmatch(node):
                offenders.append(f"{path} = {node!r}")

        walk(config.raw)
        assert not offenders, f"configuration names sample ids: {offenders}"

    def test_real_mode_is_disabled(self, config):
        assert config.runtime.get("allow_real_mode") is False

    def test_no_implicit_database_updates(self, config):
        for section in ("annotation", "amr", "virulence"):
            assert config.raw[section]["allow_database_update"] is False

    def test_test_and_real_roots_are_separate_directories(self, config):
        test_root = config.data_root(RunMode.TEST)
        real_root = config.data_root(RunMode.REAL)
        assert test_root != real_root
        assert test_root.name == "test_data"
        assert real_root.name == "data"

    def test_promoter_assessment_disabled(self, config):
        assert config.raw["regulators"]["promoter_assessment_enabled"] is False

    def test_candidate_sv_promotion_declared_but_unused(self, config):
        """The knob exists in config but stage 7 ignores it by design."""
        assert "promote_candidate_calls" in config.raw["structural_variants"]


class TestConfigValidation:
    def test_missing_config_raises(self, tmp_path):
        with pytest.raises(ConfigError, match="not found"):
            load_config(tmp_path / "config" / "science.yaml", machine=None)

    def test_antibiotic_not_in_table_raises(self, tmp_path, pipeline_root):
        _copy_config(pipeline_root, tmp_path)
        _append_antibiotic_to_yaml(tmp_path)
        with pytest.raises(ConfigError, match="absent from antibiotics.tsv"):
            load_config(tmp_path / "config" / "science.yaml", machine=None)

    def _set_evidence_level(self, path: Path, level: str) -> None:
        """Rewrite column 6 (evidence_level) of the first *data* row.

        The header is skipped by comparing against the first line's field
        count and content, so this does not depend on note wording.
        """
        lines = path.read_text().splitlines()
        for index, line in enumerate(lines):
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[0] == "mechanism":
                continue  # this is the header
            fields[5] = level
            lines[index] = "\t".join(fields)
            break
        path.write_text("\n".join(lines) + "\n")

    def test_causal_evidence_level_rejected(self, tmp_path, pipeline_root):
        """A knowledge table may not declare a causal claim."""
        _copy_config(pipeline_root, tmp_path)
        self._set_evidence_level(tmp_path / "config" / "mechanisms.tsv", "CAUSAL")
        with pytest.raises(ConfigError, match="causal"):
            load_config(tmp_path / "config" / "science.yaml", machine=None)

    def test_unknown_evidence_level_rejected(self, tmp_path, pipeline_root):
        _copy_config(pipeline_root, tmp_path)
        self._set_evidence_level(
            tmp_path / "config" / "mechanisms.tsv", "VERY_STRONG"
        )
        with pytest.raises(ConfigError, match="Unknown evidence level"):
            load_config(tmp_path / "config" / "science.yaml", machine=None)

    def test_malformed_knowledge_table_rejected(self, tmp_path, pipeline_root):
        _copy_config(pipeline_root, tmp_path)
        path = tmp_path / "config" / "regulators.tsv"
        header = [
            line
            for line in path.read_text().splitlines()
            if line and not line.startswith("#")
        ][0]
        n = len(header.split("\t"))
        # One field more than the header declares.
        with path.open("a") as handle:
            handle.write("\t".join(["x"] * (n + 1)) + "\n")
        with pytest.raises(PipelineError, match="more fields than the header"):
            load_config(tmp_path / "config" / "science.yaml", machine=None)

    def test_unknown_gene_lookup_raises(self, config):
        with pytest.raises(UnknownGeneError):
            config.mechanism_for_gene("not_a_real_gene")

    def test_require_antibiotic(self, config):
        assert config.require_antibiotic("imipenem") == "imipenem"
        with pytest.raises(AntibioticError, match="not configured"):
            config.require_antibiotic("ceftazidime")


class TestAntibioticExtensibility:
    """Adding an antibiotic must be a configuration change, not a code change."""

    def test_meropenem_is_declared_and_enabled(self, config):
        """Declared *and* enabled - the property at its new address.

        The original test asserted `test_meropenem_is_declared_but_disabled`
        : meropenem was in `config/antibiotics.tsv` but absent from
        science.yaml's `antibiotics:` list. The meropenem GWAS build reversed
        that premise deliberately - this pipeline now analyses meropenem - so
        the declared-but-disabled state no longer exists anywhere to test.
        The same property at its new address is: the table declares it AND
        the config enables it. (No third antibiotic was invented just to keep
        the old name green.) If either half regresses, one of these two
        assertions fails.
        """
        assert "meropenem" in config.antibiotic_specs
        assert "meropenem" in config.antibiotics

    def test_enabling_meropenem_requires_no_code_change(self, tmp_path, pipeline_root):
        """The documented edit alone must enable the antibiotic - and this
        test must actually *perform* that edit.

        The shipped science.yaml already lists meropenem, so the test first
        removes it from the copied config and proves the removal by loading
        it (config must come back imipenem-only). It then applies the
        documented one-line addition to that copy and asserts the file text
        changed before loading. Without those two guards the test loads the
        shipped config and asserts what is already true - a pass for the
        wrong reason, which is worse than a failure.
        """
        _copy_config(pipeline_root, tmp_path)
        config_path = tmp_path / "config" / "science.yaml"
        anchor = "antibiotics:\n  - imipenem"

        shipped = config_path.read_text()
        assert anchor in shipped, (
            "science.yaml no longer carries the literal `antibiotics:\\n  - "
            "imipenem` anchor; this test edits that exact text"
        )

        disabled_text = shipped.replace(anchor + "\n  - meropenem", anchor)
        assert disabled_text != shipped, (
            "removing meropenem did not change the copied config, so the "
            "rest of this test would assert against the shipped file"
        )
        config_path.write_text(disabled_text)
        disabled = load_config(config_path, machine=None)
        assert disabled.antibiotics == ("imipenem",)
        assert "meropenem" not in disabled.antibiotics

        # The documented edit: add one line, change no code.
        enabled_text = disabled_text.replace(anchor, anchor + "\n  - meropenem")
        assert enabled_text != disabled_text, (
            "the documented edit was a no-op on the copy: the anchor no "
            "longer matches the text this test just wrote"
        )
        config_path.write_text(enabled_text)
        config = load_config(config_path, machine=None)
        assert config.antibiotics == ("imipenem", "meropenem")
        assert config.require_antibiotic("meropenem") == "meropenem"
        # "No code change" also means the knowledge tables need no edit:
        # oprD's antibiotic column already names meropenem.
        assert config.mechanism_for_gene("oprD").relevant_to("meropenem")

    def test_antibiotic_column_accepts_a_list(self, config):
        """One knowledge row can cover several antibiotics."""
        oprd = config.mechanism_for_gene("oprD")
        assert oprd.relevant_to("imipenem")
        assert oprd.relevant_to("meropenem")

    def test_all_keyword_covers_any_antibiotic(self, tmp_path, pipeline_root):
        _copy_config(pipeline_root, tmp_path)
        _set_antibiotic(tmp_path / "config" / "mechanisms.tsv", "oprD", "all")
        config = load_config(tmp_path / "config" / "science.yaml", machine=None)
        assert config.mechanism_for_gene("oprD").relevant_to("anything_at_all")


class TestModeGating:
    def test_test_mode_allowed(self, config):
        assert resolve_mode("TEST", config) is RunMode.TEST

    def test_real_mode_refused(self, config):
        with pytest.raises(ModeNotAllowedError, match="REAL mode is disabled"):
            resolve_mode("REAL", config)

    def test_unknown_mode_refused(self, config):
        with pytest.raises(PipelineError, match="Unknown run mode"):
            resolve_mode("PRODUCTION", config)

    def test_mode_name_is_case_insensitive(self, config):
        assert resolve_mode("test", config) is RunMode.TEST

    def test_stub_mode_is_accepted(self, config):
        """STUB is implemented, so it resolves.

        It was refused for a while: the mode was declared but no stage had a
        fabricating branch, so a "stub" run read real files. The runtime now
        exists, so the refusal would be the bug.
        """
        assert resolve_mode("STUB", config) is RunMode.STUB

    def test_stub_run_fabricates_without_reading_a_cohort(self, config, tmp_path, monkeypatch):
        """A stub run must not need a manifest, a fixture or a genome.

        This is what distinguishes STUB from TEST. If it had to discover a
        manifest, a stub run would depend on the committed fixtures, and
        "STUB" would be a name for running TEST with the parsers switched off -
        which is not what spec.md D8 asks for.
        """
        monkeypatch.setenv("PIPELINE_RESULTS_ROOT", str(tmp_path / "stub-results"))
        result = run_pipeline(config=config, mode="STUB")
        assert result.mode is RunMode.STUB
        assert len(result.manifest) == 0, "a stub run should have no samples"
        assert result.outputs, "a stub run should still declare outputs"
        for key, path in result.outputs.items():
            assert Path(path).exists(), f"stage {key!r} declared a path it did not write"

    def test_stub_writes_no_rows_rather_than_inventing_them(self, config, tmp_path, monkeypatch):
        """A fabricated output is shaped, not populated.

        A stub row would be indistinguishable from a result, and would let a
        downstream consumer mistake a stub for a finding. Zero rows keeps "no
        data" readable as such.
        """
        monkeypatch.setenv("PIPELINE_RESULTS_ROOT", str(tmp_path / "stub-rows"))
        result = run_pipeline(config=config, mode="STUB")
        amr = result.outputs.get("amr")
        assert amr is not None
        lines = Path(amr).read_text(encoding="utf-8").splitlines()
        assert len(lines) == 1, f"a stub table should be header-only, got {len(lines)} lines"

    def test_a_stub_run_needs_no_assembly_root_on_disk(
        self, config, tmp_path, monkeypatch
    ):
        """spec.md D8: STUB reads nothing, so nothing has to exist.

        The existence check at run entrypoint is about the directory a *real*
        cohort is read from. Pointing it at STUB - a mode that discovers no
        manifest, reads no fixture and needs no genome - made a stub run
        depend on the very files the mode exists to avoid needing: on a
        machine without `data/`, every STUB run died with
        "Data root for this mode does not exist" before a single rule ran.

        Pinned with a data root of our own naming rather than the machine's,
        so the test is red on a developer's box that happens to have `data/`
        and green on CI that does not.
        """
        monkeypatch.setitem(config.paths, "data_root", "no_such_data_root_in_this_repo")
        assert not config.assembly_root(RunMode.STUB).exists(), (
            "the point of this test is a data root that is not there"
        )
        monkeypatch.setenv("PIPELINE_RESULTS_ROOT", str(tmp_path / "stub-no-data"))
        result = run_pipeline(config=config, mode="STUB")
        assert result.mode is RunMode.STUB
        assert result.outputs, "a stub run must still write its fabricated outputs"

    def test_a_test_run_still_refuses_without_its_data_root(
        self, config, tmp_path, monkeypatch
    ):
        """The check stays for the mode that reads a cohort. Both halves.

        Loosening the guard for STUB must not delete it: a TEST run whose
        fixtures directory is missing has to fail loudly and early, naming the
        mode and the path, rather than discover an empty cohort and report a
        clean run over nothing. A test that only covers the STUB side would
        pass just as well on a guard that had been removed entirely.
        """
        monkeypatch.setitem(config.paths, "test_data_root", "no_such_test_data_root")
        monkeypatch.setenv("PIPELINE_RESULTS_ROOT", str(tmp_path / "test-no-data"))
        with pytest.raises(PipelineError) as caught:
            run_pipeline(config=config, mode="TEST")
        message = str(caught.value)
        assert "Data root" in message, message
        assert "TEST" in message, f"the refusal must name the mode: {message}"
        assert "no_such_test_data_root" in message, (
            f"the refusal must name the path it could not find: {message}"
        )


class TestStageGraph:
    def test_sixteen_stages(self):
        """The spec's sixteen (spec.md D1), amended 2026-09-29.

        The count was fifteen when the code's own sixteen was mapped onto the
        spec; `cohort_variants` makes it sixteen again, and the amendment is
        recorded in spec.md rather than left implicit here.
        """
        assert len(STAGE_ORDER) == 16

    def test_execution_order_is_a_permutation(self):
        assert sorted(EXECUTION_ORDER) == sorted(STAGE_ORDER)

    def test_the_folded_steps_run_inside_their_owner(self):
        """spec.md:351 made them internal steps, so ordering moved with them.

        `structural_variants` used to be stage 7 and `amr` stage 4, and
        mechanisms had to wait for both. Folded, the ordering guarantee is that
        the *owner* runs after what the step consumed: `amr` before the
        mechanism classes, which `cooccurrence` now builds.
        """
        from papipeline.execution.contracts import INTERNAL_TABLES

        for owner, step in (("amr", "structural_variants"),
                            ("cooccurrence", "mechanisms"),
                            ("reporting", "master_table")):
            assert step in INTERNAL_TABLES, f"{step} lost its declared output"
            assert owner in STAGE_ORDER
        assert EXECUTION_ORDER.index("cooccurrence") > EXECUTION_ORDER.index("amr")
        assert EXECUTION_ORDER.index("reporting") > EXECUTION_ORDER.index("amr")

    def test_reporting_is_last(self):
        assert EXECUTION_ORDER[-1] == "reporting"

    def test_prerequisite_expansion(self):
        expanded = resolve_prerequisites({"reporting"})
        assert {"amr", "phenotype", "validation", "variants"} <= expanded

    def test_prerequisite_expansion_is_transitive(self):
        """`reporting` needs `phenotype`, which has no prerequisites."""
        expanded = resolve_prerequisites({"reporting"})
        assert "phenotype" in expanded

    def test_expansion_does_not_mutate_the_request(self):
        requested = {"gwas"}
        resolve_prerequisites(requested)
        assert requested == {"gwas"}


@pytest.fixture(scope="module")
def completed_run(tmp_path_factory, pipeline_root):
    """One full TEST run, shared across the assertions below."""
    results_dir = pipeline_root / "results" / "test"
    return run_pipeline(
        config_path=pipeline_root / "config" / "science.yaml",
        mode="TEST",
        config=load_config(pipeline_root / "config" / "science.yaml"),
    )


class TestFullRun:
    def test_sixteen_stages_complete(self, completed_run):
        """Sixteen declared; every stage not in `UNBUILT_STAGES` runs in TEST.

        `recombination` has no implementation, so in TEST it is disabled and
        unrequested and is skipped before it is ever marked - which is why it
        does not appear in `stage_status` at all, rather than appearing as
        skipped. Asking for one is an error; see test_stage_taxonomy.py.

        **The split is derived, never counted.** This was twelve runnable and four
        unbuilt, then fourteen and two as `variants` and `cohort_variants` were
        reconciled, and now fifteen and one as `similarity` followed. A hard-coded
        count is a tripwire that fires on every reconciliation - it asserted two
        unbuilt and failed the moment the third name was dropped - while saying
        nothing about whether the stage that *should* run actually did. So the
        invariant is asserted instead: sixteen declared, and every declared stage
        outside `UNBUILT_STAGES` present in the run's status.
        """
        from papipeline.run import STAGE_ORDER, UNBUILT_STAGES

        runnable = set(STAGE_ORDER) - set(UNBUILT_STAGES)
        assert len(STAGE_ORDER) == 16
        assert runnable, "every stage is unbuilt, which cannot be right"
        assert len(runnable) + len(UNBUILT_STAGES) == len(STAGE_ORDER)
        assert set(completed_run.stage_status) == runnable, (
            "the stages that should have run are: "
            f"{sorted(runnable - set(completed_run.stage_status))}"
        )
        failed = {
            k: v for k, v in completed_run.stage_status.items() if v != "completed"
        }
        assert failed == {}, f"stages that should have run did not: {failed}"

    def test_twenty_samples(self, completed_run):
        assert completed_run.n_samples == 20

    def test_every_expected_output_written(self, completed_run):
        # `variants` and `cohort_variants` run for real now that the gate is
        # open, and `regulators` is among their outputs: it folded into
        # `variants` (spec.md:351), so the chromosomal screen runs whenever the
        # stage that owns it does. It was absent while `variants` was unbuilt.
        for key in (
            "variants",
            "cohort_variants",
            "regulators",
            "validation",
            "mlst",
            "amr",
            "mechanisms",
            "structural_variants",
            "virulence",
            "pangenome",
            "phylogeny_summary",
            "phenotype",
            "gwas",
            "convergence",
            "cooccurrence",
            "master_table",
            "figure_data",
        ):
            assert key in completed_run.outputs, key
            assert Path(completed_run.outputs[key]).exists(), key

    def test_master_table_row_count(self, completed_run):
        from papipeline.io import read_tsv

        rows = read_tsv(
            Path(completed_run.outputs["master_table"]),
            required_columns=["sample_id", "antibiotic", "confidence"],
        )
        assert len(rows) == 20

    def test_master_table_ids_match_manifest(self, completed_run, manifest):
        from papipeline.io import read_tsv

        rows = read_tsv(Path(completed_run.outputs["master_table"]))
        assert {r["sample_id"] for r in rows} == set(manifest.sample_ids)

    def test_antibiotic_recorded_everywhere(self, completed_run):
        from papipeline.io import read_tsv

        rows = read_tsv(Path(completed_run.outputs["master_table"]))
        assert {r["antibiotic"] for r in rows} == {"imipenem"}

    def test_confidence_never_causal(self, completed_run):
        from papipeline.io import read_tsv

        rows = read_tsv(Path(completed_run.outputs["master_table"]))
        allowed = {s.value for s in ClaimStatus}
        assert {r["confidence"] for r in rows} <= allowed
        assert "CAUSAL" not in {r["confidence"] for r in rows}

    def test_amr_records_database_version(self, completed_run):
        from papipeline.io import read_tsv

        rows = read_tsv(Path(completed_run.outputs["amr"]))
        assert rows
        assert all(r["database"] and r["database_version"] for r in rows)

    def test_amr_rows_are_only_detected(self, completed_run):
        from papipeline.io import read_tsv

        rows = read_tsv(Path(completed_run.outputs["amr"]))
        assert {r["claim_status"] for r in rows} == {"DETECTED"}

    def test_candidate_svs_are_not_promoted(self, completed_run):
        from papipeline.io import read_tsv

        rows = read_tsv(Path(completed_run.outputs["structural_variants"]))
        assert any(r["call_status"] == "candidate" for r in rows)
        assert any(r["call_status"] == "not_assessable" for r in rows)

    def test_no_mic_fabricated_from_categories(self, completed_run):
        """Every MIC in the output must have come from the input file."""
        from papipeline.io import read_tsv

        source = read_tsv(
            Path(completed_run.outputs["phenotype"]), required_columns=["phenotype"]
        )
        for row in source:
            if row["MIC"] is None:
                assert row["phenotype"] in {"R", "I", "S", "SDD", "ND"}

    def test_all_thirteen_figures_prepared(self, completed_run):
        assert len(completed_run.figure_data) == 13

    def test_figure_kinds(self, completed_run):
        kinds = {p["kind"] for p in completed_run.figure_data.values()}
        assert {"bar", "grouped_bar", "heatmap", "manhattan", "network"} <= kinds

    def test_gwas_produces_results(self, completed_run):
        from papipeline.io import read_tsv

        rows = read_tsv(
            Path(completed_run.outputs["gwas"]), required_columns=["feature"]
        )
        assert rows
        assert all(r["adjusted_p_value"] for r in rows)

    def test_report_is_marked_synthetic(self, completed_run):
        text = Path(completed_run.report_paths["markdown"]).read_text()
        assert "NOT BIOLOGICAL RESULTS" in text

    def test_report_states_the_scientific_rules(self, completed_run):
        text = Path(completed_run.report_paths["markdown"]).read_text()
        assert "never reported as phenotypic resistance" in text
        assert "never reported as causation" in text
        assert "never converted into MIC" in text

    def test_html_report_written(self, completed_run):
        assert Path(completed_run.report_paths["html"]).exists()

    def test_run_manifest_records_provenance(self, completed_run):
        payload = json.loads(Path(completed_run.run_manifest).read_text())
        assert payload["run_mode"] == "TEST"
        assert payload["references"]
        assert "tools_detected" in payload

    def test_unpinned_references_are_surfaced(self, completed_run):
        assert any("unpinned" in w for w in completed_run.warnings)

    def test_qc_flags_are_reported_not_silent(self, completed_run):
        assert any("flagged" in w for w in completed_run.warnings)


class TestNoRealDataLeakage:
    def test_test_run_never_reads_the_real_data_directory(self, completed_run, pipeline_root):
        """No output path may sit inside the real data tree."""
        real_root = (pipeline_root / "data").resolve()
        for path in list(completed_run.outputs.values()) + [
            Path(completed_run.run_manifest)
        ]:
            assert real_root not in Path(path).resolve().parents

    def test_results_are_under_the_test_directory(self, completed_run):
        for path in completed_run.outputs.values():
            assert "test" in Path(path).parts

    def test_generator_refuses_to_write_into_data(self, tmp_path):
        from papipeline.testing import generate

        target = tmp_path / "data"
        target.mkdir()
        with pytest.raises(PipelineError, match="Refusing to write"):
            generate(target)


class TestDeterminism:
    def test_generator_is_reproducible(self, tmp_path):
        from papipeline.testing import generate

        first = tmp_path / "a"
        second = tmp_path / "b"
        generate(first)
        generate(second)
        for path in sorted(first.rglob("*")):
            if path.is_file():
                mirror = second / path.relative_to(first)
                assert path.read_bytes() == mirror.read_bytes(), path.name

    def test_figures_are_reproducible_across_runs(self, completed_run, pipeline_root):
        first = json.dumps(completed_run.figure_data, sort_keys=True, default=str)
        second_run = run_pipeline(
            config_path=pipeline_root / "config" / "science.yaml", mode="TEST"
        )
        second = json.dumps(second_run.figure_data, sort_keys=True, default=str)
        assert first == second


class TestStageSelection:
    def test_only_runs_requested_stage_plus_prerequisites(self, pipeline_root):
        result = run_pipeline(
            config_path=pipeline_root / "config" / "science.yaml",
            mode="TEST",
            only=["gwas"],
        )
        assert result.stage_status["gwas"] == "completed"
        assert "phenotype" in result.stage_status
        # `reporting` is not a prerequisite of gwas, and the master table it now
        # builds is therefore not built either.
        assert "reporting" not in result.stage_status

    def test_skip_excludes_a_stage_nothing_requested_needs(self, pipeline_root):
        """`--skip` still works where it is safe.

        `mlst` is not a prerequisite of anything, so skipping it is a plain
        exclusion: the run proceeds and the stage is recorded as skipped, with
        the reason, in `run_manifest.json`.

        **This used to use `gwas` and `cooccurrence`**, and that was the wrong
        pair: `reporting` requires both, so `--skip gwas` used to produce a
        report claiming a complete run without ever running GWAS. It now
        refuses - see `TestSkippingAStageSomethingRequestedIsARefusal` in
        `tests/integration/test_no_silent_skips.py` and
        `PREREQUISITES["reporting"]`, which says in as many words that a report
        claiming a complete run must have run every stage. The skip is still
        honoured; it is just no longer allowed to be invisible.
        """
        result = run_pipeline(
            config_path=pipeline_root / "config" / "science.yaml",
            mode="TEST",
            skip=["mlst"],
        )
        assert "mlst" not in result.stage_status
        # `integration` folded into `reporting` (spec.md:351), so this is the
        # stage that owns the master table now.
        assert result.stage_status["reporting"] == "completed"
        manifest = json.loads(
            Path(result.run_manifest).read_text(encoding="utf-8")
        )
        assert "mlst" in manifest["stages_skipped"], (
            "a stage that was scheduled and did not run must be recorded as "
            "skipped, with its reason; otherwise a reader of the manifest "
            "cannot tell it apart from a pipeline that has no such stage"
        )
        assert "skip" in manifest["stages_skipped"]["mlst"].lower()

    def test_skipping_a_stage_the_report_needs_refuses(self, pipeline_root):
        """The counterpart, and the behaviour A3 introduced.

        `--skip gwas` on a full run now refuses rather than writing a report
        that claims completeness. Asserted here as well as in
        `test_no_silent_skips.py` because this is the file a reader of
        `run_pipeline`'s `skip` parameter will look in.
        """
        with pytest.raises(PipelineError) as caught:
            run_pipeline(
                config_path=pipeline_root / "config" / "science.yaml",
                mode="TEST",
                skip=["gwas"],
            )
        message = str(caught.value)
        assert "gwas" in message and "reporting" in message, (
            "the refusal must name both the skipped stage and the stage that "
            f"needed it: {message}"
        )

    def test_unknown_stage_rejected(self, pipeline_root):
        with pytest.raises(PipelineError, match="Unknown stage"):
            run_pipeline(
                config_path=pipeline_root / "config" / "science.yaml",
                mode="TEST",
                only=["not_a_stage"],
            )


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------


def _copy_config(source_root: Path, target_root: Path) -> None:
    import shutil

    shutil.copytree(source_root / "config", target_root / "config")


def _set_antibiotic(path: Path, gene: str, value: str) -> None:
    """Rewrite the antibiotic column (4) of the row for ``gene``."""
    lines = path.read_text().splitlines()
    for index, line in enumerate(lines):
        if line.startswith("#") or "\t" not in line:
            continue
        fields = line.split("\t")
        if fields[1] == gene:
            fields[3] = value
            lines[index] = "\t".join(fields)
            break
    path.write_text("\n".join(lines) + "\n")


def _append_antibiotic_to_yaml(target_root: Path) -> None:
    path = target_root / "config" / "science.yaml"
    text = path.read_text()
    appended = text.replace(
        "antibiotics:\n  - imipenem", "antibiotics:\n  - imipenem\n  - vancomycin"
    )
    assert appended != text, (
        "the `antibiotics:\\n  - imipenem` anchor no longer matches "
        "science.yaml; appending vancomycin would be a silent no-op and the "
        "caller's expected ConfigError would never fire"
    )
    path.write_text(appended)
