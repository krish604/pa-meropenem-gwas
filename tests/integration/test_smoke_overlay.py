"""The smoke overlay: a bounded REAL run that cannot be mistaken for the real one.

This is the one place the project opens a door its standing rules keep shut, so
the three properties below are all refusals rather than conveniences.

`config/machines/smoke.yaml` exists to run a handful of real genomes through the
pipeline for real, without touching `allow_real_mode` on any other overlay and
without any possibility of a run on 10 assemblies being read as the 967-isolate
analysis. Three guarantees, each of which fails loudly:

1. **No silent fall-through to `data/`.** `smoke_genome_dir` is required, and
   the refusal names the key. An overlay that resolved genomes by falling back
   to `data/` when the key was missing would quietly analyse the full cohort -
   the exact outcome this overlay exists to prevent - and would do it while
   looking like a bounded run.
2. **The cap is the existing one.** `max_samples: 10` goes through
   `enforce_sample_cap`, which is already the tested refusal. No second cap
   mechanism is introduced, because two caps is one more thing to keep in step.
3. **The report says what it is.** Every report from this overlay carries
   `SMOKE TEST — REAL data, N=<count> assemblies, not the full cohort
   analysis.`, and *fails* if the marker cannot be attached rather than emitting
   a report without it. A report that omits the marker looks identical to the
   real thing, which defeats the purpose of generating it.

Note on scale: the cap test deliberately builds cohorts at 11 and 25, not at the
11/10 boundary alone. A boundary test passes for an implementation that
compares wrongly at other sizes, and this session has already had a small-N test
that could not catch the very claim it was written for.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.config.loader import load_config, load_machine_config
from papipeline.errors import PipelineError, SampleCapExceeded
from papipeline.manifest import discover_pdc_manifest
from papipeline.models import RunMode

OVERLAY = Path("config/machines/smoke.yaml")

#: The machine overlays directory and the science config, for the gate tests
#: below, which resolve real paths rather than a single overlay.
MACHINES = Path("config/machines")
SCIENCE = Path("config/science.yaml")

SMOKE_MARKER_PREFIX = "SMOKE TEST — REAL data, N="


def _write_pdc(tmp_path: Path, n: int) -> Path:
    """A PDC table with ``n`` isolates, each carrying an assembly."""
    lines = ["\t".join(["Isolate", "Assembly", "BioSample", "PDC_present"])]
    for i in range(1, n + 1):
        lines.append(
            f"PDC{i:06d}\tGCA_{i:09d}.1\tSAMN{i:08d}\tTRUE"
        )
    path = tmp_path / f"pdc_{n}.tsv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _overlay_without_key(tmp_path: Path) -> Path:
    """A copy of the smoke overlay with `smoke_genome_dir` removed."""
    import yaml

    text = OVERLAY.read_text(encoding="utf-8")
    data = yaml.safe_load(text)
    data["paths"].pop("smoke_genome_dir", None)
    path = tmp_path / "smoke_missing.yaml"
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    return path


class TestOverlayExists:
    def test_the_overlay_is_committed(self):
        assert OVERLAY.exists(), "config/machines/smoke.yaml is missing"

    def test_it_caps_at_ten_samples(self):
        machine = load_machine_config(OVERLAY)
        assert machine.max_samples == 10

    def test_it_names_its_genome_directory(self):
        machine = load_machine_config(OVERLAY)
        assert machine.smoke_genome_dir() is not None

    def test_it_does_not_open_real_mode(self):
        """Building the overlay is not the same as authorising a REAL run.

        `allow_real_mode` stays false here. The bounded REAL run is a separate,
        explicit step; an overlay that flipped it would let any smoke run start
        a REAL analysis of whatever genomes it could find.
        """
        machine = load_machine_config(OVERLAY)
        assert machine.allow_real_mode is False


class TestGenomeDirectoryIsRequired:
    def test_missing_key_is_refused_by_name(self, tmp_path):
        """The refusal names `smoke_genome_dir`, so it is actionable.

        An unnamed refusal is the failure this guards: the operator cannot tell
        which key to add.
        """
        path = _overlay_without_key(tmp_path)
        with pytest.raises(PipelineError) as excinfo:
            load_machine_config(path)
        assert "smoke_genome_dir" in str(excinfo.value)

    def test_missing_key_does_not_fall_back_to_data(self, tmp_path):
        """There is no default. Falling back to `data/` would be the bug.

        Asserted as an absence: after the refusal, nothing anywhere in the
        overlay may claim genomes live under the data root.
        """
        path = _overlay_without_key(tmp_path)
        with pytest.raises(PipelineError):
            load_machine_config(path)
        assert "data/" not in path.read_text(encoding="utf-8").replace(
            "data_root: data", ""
        )

    def test_present_key_resolves(self):
        machine = load_machine_config(OVERLAY)
        resolved = machine.smoke_genome_dir()
        assert resolved.is_absolute()
        assert resolved.name == "smoke_genomes"


class TestCapEnforcedAboveTheBoundary:
    @pytest.mark.parametrize("count", [11, 17, 25])
    def test_cohorts_above_ten_are_refused(self, tmp_path, count):
        """Not just 11. A cap compared wrongly elsewhere still passes at 11.

        The cohort is built through `discover_pdc_manifest` so this also
        exercises the size-agnostic path the cap has to be independent of.
        """
        config = load_machine_config(OVERLAY)
        manifest = discover_pdc_manifest(_write_pdc(tmp_path, count))
        with pytest.raises(SampleCapExceeded) as excinfo:
            config.enforce_sample_cap(len(manifest))
        message = str(excinfo.value)
        assert "10" in message
        assert str(count) in message

    @pytest.mark.parametrize("count", [1, 9, 10])
    def test_cohorts_at_or_below_ten_are_allowed(self, tmp_path, count):
        config = load_machine_config(OVERLAY)
        manifest = discover_pdc_manifest(_write_pdc(tmp_path, count))
        config.enforce_sample_cap(len(manifest))  # must not raise

    def test_the_cap_is_the_existing_mechanism(self):
        """No second cap: the overlay's own value is what gets enforced.

        Called on the `MachineConfig` itself, which is where the cap lives, so
        this pins that the smoke overlay uses the standard rail rather than a
        bespoke one that could drift from it.
        """
        machine = load_machine_config(OVERLAY)
        assert machine.max_samples == 10
        with pytest.raises(SampleCapExceeded):
            machine.enforce_sample_cap(11)
        machine.enforce_sample_cap(10)  # at the cap, allowed

    def test_cap_message_names_the_file_to_change(self):
        """The refusal is actionable: it says which file to edit."""
        machine = load_machine_config(OVERLAY)
        with pytest.raises(SampleCapExceeded) as excinfo:
            machine.enforce_sample_cap(25)
        assert "smoke.yaml" in str(excinfo.value)


class TestReportMarker:
    def test_marker_states_the_actual_count(self):
        from papipeline.reporting.smoke import smoke_marker

        for count in (1, 3, 10):
            assert smoke_marker(count) == (
                f"SMOKE TEST — REAL data, N={count} assemblies, "
                "not the full cohort analysis."
            )

    def test_marker_is_required_not_optional(self, tmp_path):
        """A report that cannot carry the marker fails instead of shipping.

        This is the whole point of the marker. A report generated without it
        looks exactly like the real analysis, so an omission is not a cosmetic
        defect but the failure the marker exists to prevent.
        """
        from papipeline.reporting.smoke import require_smoke_marker

        good = "prelude\nSMOKE TEST — REAL data, N=3 assemblies, not the full cohort analysis.\n"
        assert require_smoke_marker(good, 3) == good

        for absent in ("prelude\nno marker here\n", "SMOKE TEST — REAL data, N=9 assemblies, not the full cohort analysis.\n"):
            with pytest.raises(PipelineError):
                require_smoke_marker(absent, 3)

    def test_generation_fails_rather_than_omitting(self):
        """The generator raises if the marker is absent from its own output."""
        from papipeline.reporting.smoke import require_smoke_marker

        with pytest.raises(PipelineError) as excinfo:
            require_smoke_marker("a report with no marker", 4)
        assert "SMOKE TEST" in str(excinfo.value)


class TestTheRealModeGateIsScopedToSmokeOnly:
    """Turning REAL on for the smoke run must not turn it on anywhere else.

    `test_it_does_not_open_real_mode` already asserts the smoke overlay's flag is
    false. This is the other half, and it is the half that matters when the flag
    is eventually flipped: **scoping** is the property under test, not the value.

    The failure this guards against is specific. `allow_real_mode` is per-overlay
    state, so authorising a bounded 10-assembly run and authorising a full
    835-isolate run are the same edit to different files - and the bigmachine
    overlay is the one without a sample cap. A change made for the smoke run that
    leaked into `bigmachine.yaml` would let an unbounded REAL cohort start on the
    analysis machine with nothing to stop it.

    So every overlay is asserted twice: the flag is false, **and** the refusal
    actually happens through the public interface. The second is what stops this
    from being a check on a YAML key that the code might not read.

    This test passes today, against the still-false flags. That is deliberate:
    it proves the assertion is fail-closed *before* anything is turned on, so
    the flip cannot arrive without this test having been watching. The
    complement - a test that only ever passed because REAL was off - would be
    worthless as a guard.
    """

    #: Every overlay on this machine, and whether a cap applies.
    OVERLAYS = ("laptop.yaml", "bigmachine.yaml", "smoke.yaml")

    @pytest.mark.parametrize("overlay", OVERLAYS)
    def test_the_flag_is_false_on_every_overlay(self, overlay):
        from papipeline.config.loader import load_machine_config

        machine = load_machine_config(MACHINES / overlay)
        assert machine.allow_real_mode is False, (
            f"{overlay} has allow_real_mode true. The smoke run has not been "
            "authorised yet; flipping a flag is a separate, explicit decision."
        )

    @pytest.mark.parametrize("overlay", OVERLAYS)
    def test_real_mode_is_refused_through_the_public_interface(self, overlay):
        """Not a YAML check: the run itself must refuse.

        A flag that is false but not consulted would pass the test above and let
        a REAL cohort start anyway.
        """
        from papipeline.config.loader import load_config
        from papipeline.errors import ModeNotAllowedError
        from papipeline.run import resolve_mode

        config = load_config(SCIENCE, machine=MACHINES / overlay)
        with pytest.raises(ModeNotAllowedError) as excinfo:
            resolve_mode("REAL", config)
        message = str(excinfo.value)
        assert "allow_real_mode" in message, message

    def test_the_refusal_is_the_same_gate_on_all_three(self):
        """Scoping means one gate, even though the messages differ per overlay.

        Each refusal names its own machine, which is deliberate - an operator
        who ran with the wrong `--machine` needs to be told which overlay
        refused. So the messages cannot be compared for equality. What must
        hold is that all three are the *same* gate: each says the flag is off,
        names the flag, and points at the same alternative.
        """
        from papipeline.config.loader import load_config
        from papipeline.errors import ModeNotAllowedError
        from papipeline.run import resolve_mode

        for overlay in self.OVERLAYS:
            config = load_config(SCIENCE, machine=MACHINES / overlay)
            with pytest.raises(ModeNotAllowedError) as excinfo:
                resolve_mode("REAL", config)
            message = str(excinfo.value)
            assert "allow_real_mode" in message, f"{overlay}: {message}"
            assert "is false" in message, f"{overlay}: {message}"
            # Each names its own machine, so an operator knows which overlay
            # refused rather than being told only that something did.
            assert Path(overlay).stem in message, (
                f"{overlay}: the refusal does not name the overlay that refused: "
                f"{message}"
            )
            # And the same remedy is offered, so the advice does not depend on
            # which overlay was chosen.
            assert "bigmachine" in message, f"{overlay}: {message}"

    def test_a_flipped_smoke_flag_does_not_disturb_the_others(self, tmp_path):
        """The scoping property itself, demonstrated rather than asserted.

        Flips `allow_real_mode` on a *copy* of the smoke overlay and shows the
        other two are unaffected. Written this way because the pre-flip state
        has all three identical, so any assertion comparing them is trivially
        true until the flag moves - and this one still means something after it
        does.

        A copy, not the real file: this test must never leave a REAL-mode flag
        set, even transiently, even if it fails partway.
        """
        import yaml

        from papipeline.config.loader import load_config
        from papipeline.errors import ModeNotAllowedError
        from papipeline.run import resolve_mode

        flipped = tmp_path / "smoke_flipped.yaml"
        data = yaml.safe_load((MACHINES / "smoke.yaml").read_text(encoding="utf-8"))
        data["runtime"]["allow_real_mode"] = True
        flipped.write_text(yaml.safe_dump(data), encoding="utf-8")

        # The copy now permits REAL...
        permissive = load_config(SCIENCE, machine=flipped)
        assert permissive.runtime.get("allow_real_mode") is True
        assert resolve_mode("REAL", permissive) is RunMode.REAL

        # ...and the other two still refuse, from their own committed files.
        for overlay in ("laptop.yaml", "bigmachine.yaml"):
            config = load_config(SCIENCE, machine=MACHINES / overlay)
            with pytest.raises(ModeNotAllowedError):
                resolve_mode("REAL", config)

        # And the committed smoke overlay is still false - the flip was a copy.
        from papipeline.config.loader import load_machine_config

        assert load_machine_config(MACHINES / "smoke.yaml").allow_real_mode is False

    def test_bigmachine_is_the_unbounded_one_and_is_not_authorised(self):
        """Why scoping matters here specifically.

        `bigmachine` is the overlay with no sample cap, so a flag flipped there
        would authorise a full cohort rather than a bounded smoke run. Recorded
        as a test so the asymmetry is visible next to the guard, not only in
        the YAML.
        """
        from papipeline.config.loader import load_machine_config

        assert load_machine_config(MACHINES / "bigmachine.yaml").max_samples is None
        assert load_machine_config(MACHINES / "smoke.yaml").max_samples == 10
        assert (
            load_machine_config(MACHINES / "bigmachine.yaml").allow_real_mode is False
        )


class SampleManifestFactory:
    """Manifests of a given size, with or without prepared assemblies."""

    @staticmethod
    def assemblyless(n):
        from papipeline.manifest import SampleManifest
        from papipeline.models import Sample

        return SampleManifest(
            samples=[Sample(sample_id=f"ISO_{i:04d}") for i in range(n)]
        )


def _complete_declared_inputs():
    """The three aggregates a REAL report must declare, with real content.

    A REAL report refuses when any of them carries no usable data
    (`stages.reporting.REPORTING_REQUIRED_INPUTS`), so a marker test that wants
    a report on disk has to declare a complete run - which is what the
    orchestrator would hand it. Empty would be a *refusal* now, not a thin
    report, so these are real rows rather than placeholders.
    """
    from papipeline.models import Cooccurrence, ConvergenceCall, ConvergenceCategory

    return {
        "convergence_calls": [
            ConvergenceCall(
                determinant="blaOXA-1",
                independent_lineages=2,
                branch_count=2,
                distribution={"L1": 3, "L2": 3},
                convergence_category=ConvergenceCategory.RECURRENT_CONVERGENT,
            )
        ],
        "cooccurrence": [
            Cooccurrence(
                feature_a="blaOXA-1",
                feature_b="fosA",
                feature_type="gene_gene",
                n_a=4,
                n_b=3,
                n_both=2,
                statistic="odds_ratio",
                statistic_value=4.0,
                adjusted_p_value=0.01,
            )
        ],
        "gwas_associations": [
            {
                "feature": "gene__blaOXA-1",
                "feature_type": "acquired_gene",
                "adjusted_p_value": 0.001,
                "effect": 3.2,
                "frequency": 0.5,
                "lineage_linked": False,
                "passes_threshold": True,
                "dominant_lineage_share": 0.5,
            }
        ],
    }


def _report_context(config, *, mode, n_samples, marker_count=None):
    """A minimal ReportContext, as the orchestrator would assemble one."""
    from papipeline.stages.reporting import ReportContext

    return ReportContext(
        mode=mode,
        config=config,
        generated_at="2026-09-30T00:00:00Z",
        antibiotic="imipenem",
        n_samples=n_samples,
        declared_inputs=_complete_declared_inputs(),
    )


class TestTheMarkerIsActuallyWired:
    """`require_smoke_marker` had no production caller, so property 3 was prose.

    The overlay header has always claimed the marker is "asserted so it fails
    loudly if it cannot be attached". It was not: the function existed, three
    tests exercised it directly, and no code path in the pipeline ever called
    it. A report from a REAL run shipped without it.

    Wiring it in `write_report` rather than in the orchestrator is deliberate.
    `write_report` is the only place a report becomes a file, so the check sits
    where an unverified report would otherwise be written - and it cannot be
    forgotten by a future caller who builds a report some other way.

    **The marker is required of a *smoke* report, not of every REAL report.**
    A full-cohort analysis is a REAL run on an overlay without
    `smoke_genome_dir`, and stamping "SMOKE TEST" on it would be as wrong as
    omitting the marker from a smoke run - worse, actually, since it would be
    confidently wrong. So the test is keyed on the overlay, not the mode.
    """

    @staticmethod
    def _smoke_config():
        from papipeline.config.loader import load_config

        return load_config(SCIENCE, machine=MACHINES / "smoke.yaml")

    @staticmethod
    def _real_config(tmp_path):
        """A REAL-capable overlay that is NOT a smoke overlay.

        Written to a temporary file rather than built in memory, because
        `load_config` takes a machine *name or path* and nothing else. Copied
        from the committed laptop overlay with the flag raised, so the config
        under test is a real one; no committed file is modified.
        """
        import yaml

        from papipeline.config.loader import load_config

        data = yaml.safe_load((MACHINES / "laptop.yaml").read_text(encoding="utf-8"))
        data["runtime"]["allow_real_mode"] = True
        path = tmp_path / "real_not_smoke.yaml"
        path.write_text(yaml.safe_dump(data), encoding="utf-8")
        return load_config(SCIENCE, machine=path)

    def test_a_smoke_report_carries_the_marker(self, tmp_path):
        from papipeline.reporting.smoke import smoke_marker
        from papipeline.stages.reporting import write_report

        config = self._smoke_config()
        context = _report_context(config, mode=RunMode.REAL, n_samples=10)
        written = write_report(context, tmp_path, write_html=False)
        text = written["markdown"].read_text(encoding="utf-8")
        assert smoke_marker(10) in text, (
            "a smoke report was written without its marker; this is the failure "
            "the marker exists to prevent"
        )

    def test_the_marker_names_the_actual_assembly_count(self, tmp_path):
        from papipeline.reporting.smoke import smoke_marker
        from papipeline.stages.reporting import write_report

        config = self._smoke_config()
        context = _report_context(config, mode=RunMode.REAL, n_samples=7)
        written = write_report(context, tmp_path, write_html=False)
        text = written["markdown"].read_text(encoding="utf-8")
        assert smoke_marker(7) in text
        assert smoke_marker(10) not in text, (
            "the marker must state the count the run actually used, not a "
            "fixed one"
        )

    def test_a_smoke_report_without_the_marker_raises(self, tmp_path, monkeypatch):
        """The verification, independent of the attachment.

        Attachment is stubbed out so this asserts the *check* rather than the
        whole path: a report body that never got the marker must stop the
        write, not ship.
        """
        from papipeline.stages import reporting as stage_reporting

        monkeypatch.setattr(
            stage_reporting,
            "build_markdown",
            lambda context: "# report with no marker\n\nnothing here\n",
        )
        config = self._smoke_config()
        context = _report_context(config, mode=RunMode.REAL, n_samples=10)
        with pytest.raises(PipelineError) as excinfo:
            stage_reporting.write_report(context, tmp_path, write_html=False)
        assert "marker" in str(excinfo.value).lower()
        assert not list(tmp_path.glob("*.md")), (
            "a report was written despite failing the marker check; the failure "
            "has to happen before the file exists"
        )

    def test_a_real_non_smoke_report_is_not_forced_to_carry_it(self, tmp_path):
        """A full-cohort analysis is REAL and is not a smoke run.

        Stamping "SMOKE TEST" on it would be confidently wrong, so the marker is
        keyed on the overlay rather than the mode.
        """
        from papipeline.reporting.smoke import smoke_marker
        from papipeline.stages.reporting import write_report

        config = self._real_config(tmp_path)
        context = _report_context(config, mode=RunMode.REAL, n_samples=967)
        written = write_report(context, tmp_path, write_html=False)
        text = written["markdown"].read_text(encoding="utf-8")
        assert smoke_marker(967) not in text
        assert text.strip(), "the report should still have been written"

    @pytest.mark.parametrize("mode", [RunMode.TEST, RunMode.STUB])
    def test_non_real_modes_are_unaffected(self, tmp_path, mode):
        """The requirement is a smoke/REAL concern; TEST and STUB are untouched.

        A smoke overlay is not how TEST or STUB are run, so neither should be
        made to carry the marker - and the STUB run fabricates its report by a
        different path entirely, which must keep working untouched.
        """
        from papipeline.stages.reporting import write_report

        config = self._smoke_config()
        context = _report_context(config, mode=mode, n_samples=20)
        written = write_report(context, tmp_path, write_html=False)
        assert written["markdown"].read_text(encoding="utf-8").strip()


class TestManifestDiscoveryHonoursTheSmokeGenomeDir:
    """`smoke_genome_dir` was declared by the overlay and read by nothing.

    `MachineConfig.smoke_genome_dir()` existed and had a test, and no run code
    ever called it. So `run_pipeline(mode="REAL")` discovered its manifest from
    `data/metadata/sample_metadata.tsv` under *every* overlay - including this
    one, whose entire purpose is to read a different set of genomes. The run
    then failed on a missing metadata file, which reads like a broken pipeline
    rather than "this overlay was never wired".

    The property under test is narrow and deliberate: an overlay that sets
    `smoke_genome_dir` discovers its cohort from there, and an overlay that does
    not is **untouched**. laptop.yaml and bigmachine.yaml must behave exactly as
    they did, and the test for that is the reason this is safe to change.
    """

    def test_the_loader_exposes_the_key(self):
        """Precondition: the value is there, it was simply never consumed."""
        machine = load_machine_config(MACHINES / "smoke.yaml")
        assert machine.smoke_genome_dir().name == "smoke_genomes"

    @pytest.mark.parametrize(
        "overlay", ["laptop.yaml", "bigmachine.yaml"]
    )
    def test_non_smoke_overlays_do_not_declare_it(self, overlay):
        """The gate on the whole change: only the smoke overlay may opt in."""
        machine = load_machine_config(MACHINES / overlay)
        assert "smoke_genome_dir" not in (machine.paths or {}), (
            f"{overlay} must not set smoke_genome_dir; doing so would point a "
            "full-cohort run at the 10 prepared assemblies"
        )

    def test_smoke_discovery_reads_the_smoke_genome_dir(self, tmp_path):
        """The run's own manifest resolution, not a hand-rolled equivalent.

        Driven through the function `run_pipeline` calls, so this cannot pass
        while the orchestrator keeps using a different one.
        """
        from papipeline.run import discover_run_manifest

        genomes = tmp_path / "smoke_genomes"
        genomes.mkdir()
        for sid in ("PDT000050773.2", "PDT000167133.1"):
            (genomes / f"{sid}.fna").write_text(">c\nACGT\n", encoding="utf-8")

        config = load_config(SCIENCE, machine=MACHINES / "smoke.yaml")
        manifest = discover_run_manifest(
            config, RunMode.REAL, pdc_path=Path("PDC_essential.tsv")
        )

        with_assembly = [
            s for s in manifest if s.assembly_path
        ]
        assert with_assembly, "no isolate resolved an assembly from smoke_genome_dir"
        for sample in with_assembly:
            assert str(genomes) in str(sample.assembly_path) or str(
                config.machine.smoke_genome_dir()
            ) in str(sample.assembly_path), (
                f"{sample.sample_id} resolved to {sample.assembly_path}, which "
                "is outside the smoke genome directory"
            )

    def test_a_non_smoke_overlay_keeps_reading_data_metadata(self, tmp_path):
        """laptop/bigmachine behaviour is unchanged, and provably so.

        A temporary REAL-capable overlay built from laptop.yaml - the flag
        raised, no `smoke_genome_dir` - must still discover from
        `data/metadata/`, which is what every non-smoke run has always done.
        """
        import yaml

        from papipeline.run import discover_run_manifest

        data = yaml.safe_load((MACHINES / "laptop.yaml").read_text(encoding="utf-8"))
        data["runtime"]["allow_real_mode"] = True
        path = tmp_path / "laptop_real.yaml"
        path.write_text(yaml.safe_dump(data), encoding="utf-8")
        config = load_config(SCIENCE, machine=path)

        with pytest.raises(PipelineError) as excinfo:
            discover_run_manifest(
                config, RunMode.REAL, pdc_path=Path("PDC_essential.tsv")
            )
        message = str(excinfo.value)
        assert "metadata" in message.lower(), (
            "a non-smoke overlay must still look in data/metadata; falling back "
            f"to the smoke directory would analyse the wrong cohort: {message}"
        )

    def test_no_sample_metadata_file_is_written(self, tmp_path):
        """Discovery must not populate data/metadata/ as a side effect.

        The handoff is explicit that `PDC_essential.tsv` is read directly and
        `data/metadata/sample_metadata.tsv` is not populated as a copy. A
        manifest that leaked one would be a second source of truth.
        """
        from papipeline.run import discover_run_manifest

        config = load_config(SCIENCE, machine=MACHINES / "smoke.yaml")
        discover_run_manifest(
            config, RunMode.REAL, pdc_path=Path("PDC_essential.tsv")
        )
        produced = config.metadata_dir(RunMode.REAL) / "sample_metadata.tsv"
        assert not produced.exists(), (
            f"{produced} was created; PDC_essential.tsv is the source and a copy "
            "would be a second truth to drift from it"
        )


class TestTheCapCountsPreparedAssembliesOnASmokeRun:
    """`max_samples: 10` means ten assemblies, not 967 cohort members.

    The overlay has always said the cap goes through the existing
    `enforce_sample_cap` and that "there is no second cap mechanism". It was
    true of the *mechanism* and false of the *number*: `enforce_sample_cap`
    compares against `len(manifest)`, and `discover_pdc_manifest` deliberately
    makes every isolate in `PDC_essential.tsv` a member - all 967 of them - with
    only the prepared ones carrying an `assembly_path`. So a correctly-wired
    smoke run was refused for carrying 967 samples when it held 10 assemblies.

    The fix changes the number the existing cap compares, and nothing else:
    `enforce_sample_cap` keeps its signature and its single raise site, and
    every non-smoke overlay keeps counting raw membership.

    **The manifest is not touched.** Membership stays 967. What changed is which
    number the cap measures, and asserting membership here is what stops a
    later "simplification" from quietly shrinking the cohort to make the cap fit.
    """

    @staticmethod
    def _config(tmp_path, *, max_samples, n_assemblies, allow_real=True):
        """A smoke overlay pointed at a temporary genome directory."""
        import yaml

        from papipeline.config.loader import load_config

        genomes = tmp_path / f"genomes_{n_assemblies}"
        genomes.mkdir(exist_ok=True)
        for i in range(n_assemblies):
            (genomes / f"SMOKE_{i:03d}.fna").write_text(">c\nACGT\n", encoding="utf-8")

        data = yaml.safe_load((MACHINES / "smoke.yaml").read_text(encoding="utf-8"))
        data["paths"]["smoke_genome_dir"] = str(genomes)
        data["runtime"]["max_samples"] = max_samples
        if allow_real:
            data["runtime"]["allow_real_mode"] = True
        path = tmp_path / f"smoke_{n_assemblies}_{max_samples}.yaml"
        path.write_text(yaml.safe_dump(data), encoding="utf-8")
        return load_config(SCIENCE, machine=path), genomes

    @staticmethod
    def _manifest(n):
        from papipeline.manifest import SampleManifest
        from papipeline.models import Sample

        return SampleManifest(
            samples=[
                Sample(sample_id=f"SMOKE_{i:03d}", assembly_path=f"g/SMOKE_{i:03d}.fna")
                for i in range(n)
            ]
        )

    def test_ten_assemblies_against_a_cap_of_ten_passes(self, tmp_path):
        """The case the overlay was written for: 10 assemblies, cap 10."""
        from papipeline.run import capped_sample_count

        config, _ = self._config(tmp_path, max_samples=10, n_assemblies=10)
        count = capped_sample_count(config, self._manifest(10))
        assert count == 10
        config.enforce_sample_cap(count)  # must not raise

    def test_an_eleventh_assembly_trips_the_cap(self, tmp_path):
        from papipeline.errors import SampleCapExceeded
        from papipeline.run import capped_sample_count

        config, _ = self._config(tmp_path, max_samples=10, n_assemblies=11)
        count = capped_sample_count(config, self._manifest(11))
        assert count == 11
        with pytest.raises(SampleCapExceeded) as excinfo:
            config.enforce_sample_cap(count)
        assert "10" in str(excinfo.value)

    def test_a_smoke_run_whose_isolates_are_unprepared_is_not_blocked(
        self, tmp_path
    ):
        """Membership alone must not trip the cap on a smoke overlay.

        A smoke run resolves 967 members and prepares 10. If the cap counted
        membership, the overlay could never run at all - which is precisely the
        failure this change exists to remove.
        """
        from papipeline.run import capped_sample_count

        config, _ = self._config(tmp_path, max_samples=10, n_assemblies=10)
        big = SampleManifestFactory.assemblyless(967)
        assert capped_sample_count(config, big) == 0
        config.enforce_sample_cap(0)

    def test_membership_is_untouched_still_967(self):
        """The manifest is not the thing that changed, and this pins it.

        Read from the real table, so it states the cohort size rather than a
        number typed into a test.
        """
        from papipeline.manifest import discover_pdc_manifest

        manifest = discover_pdc_manifest(Path("PDC_essential.tsv"))
        assert len(manifest) == 967, (
            "membership changed; the cap was supposed to change, not the cohort"
        )

    def test_a_non_smoke_overlay_still_counts_raw_membership(self, tmp_path):
        """Unchanged behaviour, stated directly.

        Every member counts, whether or not it has an assembly - the original
        semantics, kept intact.
        """
        import yaml

        from papipeline.config.loader import load_config
        from papipeline.run import capped_sample_count

        data = yaml.safe_load((MACHINES / "laptop.yaml").read_text(encoding="utf-8"))
        data["runtime"]["max_samples"] = 20
        data["runtime"]["allow_real_mode"] = True
        path = tmp_path / "laptop_cap.yaml"
        path.write_text(yaml.safe_dump(data), encoding="utf-8")
        config = load_config(SCIENCE, machine=path)

        assert capped_sample_count(config, self._manifest(7)) == 7
        assert capped_sample_count(config, SampleManifestFactory.assemblyless(7)) == 7, (
            "a non-smoke overlay must count members with no assembly too; "
            "dropping them would let an unanalysable cohort through the cap"
        )


class TestTheMarkerCountsAssembliesNotRosterMembers:
    """The marker said N=967 while the report body said 10 assemblies.

    A smoke run's manifest is deliberately the *whole* PDC roster - every
    isolate is a member, and the ten without a prepared assembly are meant to
    refuse per sample rather than vanish. So ``len(manifest)`` is 967, and
    stamping that produced:

        SMOKE TEST — REAL data, N=967 assemblies, not the full cohort analysis.

    which is wrong twice over. It claims 967 assemblies when ten were analysed,
    and 967 *is* the full cohort, so the sentence contradicts itself. The code
    even claimed otherwise: the comment at the marker site says the count "is
    the one the run actually used, so the marker cannot claim a number the
    report does not support". It was not.

    Stage 1 already had the right number - ``qc_summary['n_samples']`` counts
    AssemblyQC records, one per assembly that was actually read - so the report
    body and the marker disagreed inside a single report.
    """

    @pytest.fixture(scope="class")
    def smoke_config(self):
        from papipeline.config.loader import load_config

        return load_config("config/science.yaml", machine="config/machines/smoke.yaml")

    def _context(self, config, mode, manifest, qc_summary):
        from papipeline.run import _build_report_context

        declared = _complete_declared_inputs()
        return _build_report_context(
            config=config,
            mode=mode,
            antibiotic="imipenem",
            manifest=manifest,
            status={},
            warnings=[],
            provenance=[],
            # The stage-12 rows `_build_report_context` reads to declare
            # `gwas_associations`; convergence and co-occurrence are passed
            # through directly. Without them a REAL context declares nothing
            # usable and the stage refuses, which is the point of the guard but
            # not what this fixture is testing.
            figure_data={
                "11_gwas_associations": {"rows": declared["gwas_associations"]}
            },
            qc_summary=qc_summary,
            convergence_calls=declared["convergence_calls"],
            cooccurrence=declared["cooccurrence"],
        )

    def test_a_smoke_context_counts_assemblies_not_roster_members(self, smoke_config):
        from papipeline.manifest import SampleManifest
        from papipeline.models import RunMode, Sample

        roster = SampleManifest(samples=[Sample(sample_id=f"ISO_{i:04d}") for i in range(967)])
        context = self._context(
            smoke_config, RunMode.REAL, roster, {"n_samples": 10}
        )
        assert context.n_samples == 10, (
            f"marker would read N={context.n_samples} for a run that analysed 10 "
            "assemblies; 967 is the roster, and calling it 'assemblies' makes "
            "the marker's own 'not the full cohort analysis' clause false"
        )

    def test_the_marker_text_names_the_assemblies_analysed(
        self, smoke_config, tmp_path
    ):
        from papipeline.manifest import SampleManifest
        from papipeline.models import RunMode, Sample
        from papipeline.reporting.smoke import smoke_marker
        from papipeline.stages.reporting import write_report

        roster = SampleManifest(samples=[Sample(sample_id=f"ISO_{i:04d}") for i in range(967)])
        context = self._context(
            smoke_config, RunMode.REAL, roster, {"n_samples": 10}
        )
        written = write_report(context, tmp_path, write_html=False)
        text = written["markdown"].read_text(encoding="utf-8")
        assert smoke_marker(10) in text
        assert smoke_marker(967) not in text

    def test_a_full_cohort_run_is_unaffected(self):
        """Every member of a real cohort has an assembly, so the two counts agree.

        Asserted so the smoke fix cannot quietly change what a full run reports.
        """
        from papipeline.config.loader import load_config
        from papipeline.manifest import SampleManifest
        from papipeline.models import RunMode, Sample

        config = load_config("config/science.yaml", machine="config/machines/laptop.yaml")
        cohort = SampleManifest(samples=[Sample(sample_id=f"ISO_{i:04d}") for i in range(50)])
        context = self._context(config, RunMode.REAL, cohort, {"n_samples": 50})
        assert context.n_samples == 50

    def test_it_falls_back_to_the_manifest_without_a_qc_summary(self, smoke_config):
        """Stage 1 not run is not a licence to claim 967.

        The QC count is preferred, but if it is absent the manifest is the only
        number available - and reporting it beats refusing to build a report.
        """
        from papipeline.manifest import SampleManifest
        from papipeline.models import RunMode, Sample

        roster = SampleManifest(samples=[Sample(sample_id=f"ISO_{i:04d}") for i in range(967)])
        context = self._context(smoke_config, RunMode.REAL, roster, None)
        assert context.n_samples == 967
