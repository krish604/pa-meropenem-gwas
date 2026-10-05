"""GWAS must refuse REAL rather than fall through to a mechanics-testing engine.

`stages.gwas.run` accepted `mode` and never branched on it. With no engine
injected - which is what `run.py` does - it fell through to `ReferenceEngine`:
Fisher's exact with BH correction and **no population-structure correction**,
whose own docstring says "Use for pipeline mechanics testing only".

Its results are self-labelled: the `model` field records the limitation, so a
reader can see what produced them. That is better than silence, and it is why
this was left alone while `convergence` and `cooccurrence` were fixed. But
self-labelling is not the same as being permitted. A real analysis should not be
able to run at all on an engine built for pipeline mechanics, because the
lineage confounding it cannot detect is precisely what a real cohort has most of
- and the report's `model` field is a caveat, not a control.

So REAL now refuses, matching `similarity`, `convergence` and `cooccurrence`.
TEST is untouched, and so are both engines: this is about which one REAL may
reach, not about wiring `PyseerEngine` up. That is separate work, and it depends
on `similarity`'s REAL caller existing first.

`GwasEngine`'s own `NotImplementedError` is an abstract-method stub and was
never a REAL guard - which is why this went unnoticed for as long as it did.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.errors import DataContractError
from papipeline.models import RunMode
from papipeline.stages import gwas as stage_gwas


def _fixture_ids():
    """Sample ids taken from the committed fixture, not invented.

    They have to match: `build_input` refuses a manifest sample with no usable
    feature row, so a made-up `S0` would fail on the fixture rather than on
    anything this file is testing.
    """
    path = Path("test_data/intermediate/gwas/gwas_features.tsv")
    if not path.is_file():
        return None
    ids, seen = [], set()
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        sample_id = line.split("\t", 1)[0]
        if sample_id != "sample_id" and sample_id not in seen:
            seen.add(sample_id)
            ids.append(sample_id)
    return ids


def _manifest():
    from papipeline.manifest import SampleManifest
    from papipeline.models import Sample

    return SampleManifest(
        samples=[Sample(sample_id=i) for i in (_fixture_ids() or [])]
    )


def _phenotypes():
    """Both outcome groups, so `min_samples_per_group` is satisfied.

    Otherwise a TEST-mode test would be asserting the group floor rather than
    engine dispatch.
    """
    from papipeline.models import Phenotype, PhenotypeCall

    ids = _fixture_ids() or []
    return [
        PhenotypeCall(
            sample_id=i, antibiotic="imipenem",
            phenotype=Phenotype("R" if n % 2 == 0 else "S"), mic=None,
        )
        for n, i in enumerate(ids)
    ]


class TestGwasRefusesReal:
    def test_it_raises_before_reaching_any_engine(self, config, tmp_path):
        """The guard stands between an engine-less REAL run and any engine.

        The engine choice is the thing being restricted, so the refusal has to
        come first. With no engine injected there is nothing to reach, and
        reaching `ReferenceEngine` by default is the fault.

        The paired case - REAL *with* a real engine - is no longer a refusal
        and is covered by `TestRealRunsWhenARealEngineIsInjected` below.
        """
        reached = []

        class SpyEngine(stage_gwas.GwasEngine):
            name = "spy"

            def run(self, gwas_input, cfg):
                reached.append(True)
                return []

        # `engine=None`, which is what `run.py` passes today.
        with pytest.raises(NotImplementedError) as excinfo:
            stage_gwas.run(
                config, _manifest(), RunMode.REAL, tmp_path, _phenotypes(),
            )
        assert not reached, "an engine ran despite the REAL refusal"
        assert "gwas" in str(excinfo.value).lower()

    def test_an_explicitly_injected_reference_engine_is_still_refused(
        self, config, tmp_path
    ):
        """Naming the mechanics engine is not a way past the guard.

        `engine=None` and `engine=ReferenceEngine()` are the same request - "run
        this cohort on an engine with no population-structure correction" - so
        they get the same answer. Accepting the explicit form would make the
        refusal a check on the caller's discretion.
        """
        reached = []

        class SpyReference(stage_gwas.ReferenceEngine):
            name = "spy_reference"

            def run(self, gwas_input, cfg):
                reached.append(True)
                return []

        with pytest.raises(NotImplementedError):
            stage_gwas.run(
                config, _manifest(), RunMode.REAL, tmp_path, _phenotypes(),
                engine=SpyReference(),
            )
        assert not reached, "a ReferenceEngine subclass ran in REAL mode"

    def test_the_message_says_why_and_names_the_mode(self, config, tmp_path):
        with pytest.raises(NotImplementedError) as excinfo:
            stage_gwas.run(config, _manifest(), RunMode.REAL, tmp_path, _phenotypes())
        message = str(excinfo.value)
        assert "REAL" in message
        assert "population" in message.lower() or "confound" in message.lower(), (
            "the refusal must say what is missing, or a reader cannot tell "
            "whether it is a wiring gap or a deliberate limit"
        )

    def test_it_names_the_engine_it_refuses_to_fall_through_on(
        self, config, tmp_path
    ):
        with pytest.raises(NotImplementedError) as excinfo:
            stage_gwas.run(config, _manifest(), RunMode.REAL, tmp_path, _phenotypes())
        assert "reference" in str(excinfo.value).lower(), (
            "the refusal should name ReferenceEngine, since that is the engine "
            "a REAL run would otherwise have silently used"
        )

    def test_the_message_names_pyseer_as_the_engine_to_inject(
        self, config, tmp_path
    ):
        """Not merely "not ReferenceEngine", but which engine to use.

        `PyseerEngine` has existed as a complete class for some time. A refusal
        that only said what it would not accept left a reader with nothing to
        act on, which is the same defect as the original message's "no caller
        yet" - it named a problem and not the fix.
        """
        with pytest.raises(NotImplementedError) as excinfo:
            stage_gwas.run(config, _manifest(), RunMode.REAL, tmp_path, _phenotypes())
        assert "PyseerEngine" in str(excinfo.value)

    def test_the_message_states_that_lineage_is_flagged_not_corrected(
        self, config, tmp_path
    ):
        """The surviving limitation, stated as the limitation it is.

        The stage corrects for relatedness with pyseer's random effect but does
        not put lineage in the model, so a lineage-confounded feature is
        *flagged*. A message that let a reader assume the confounding was
        removed would be the more dangerous error of the two.
        """
        with pytest.raises(NotImplementedError) as excinfo:
            stage_gwas.run(config, _manifest(), RunMode.REAL, tmp_path, _phenotypes())
        message = str(excinfo.value).lower()
        assert "lineage" in message
        assert "flag" in message, (
            "the refusal must distinguish lineage being flagged from lineage "
            "being corrected for"
        )

    def test_the_message_no_longer_claims_there_is_no_caller(self, config, tmp_path):
        """The stale clause, pinned so it cannot come back.

        `PyseerEngine` is a complete, wired class with two pyseer passes and a
        kinship matrix. "has no caller yet" described the codebase before it
        landed; a refusal that still says so explains a live block with a dead
        reason, and a reader who trusts it will not go looking for the two
        blockers that are actually current.
        """
        with pytest.raises(NotImplementedError) as excinfo:
            stage_gwas.run(config, _manifest(), RunMode.REAL, tmp_path, _phenotypes())
        assert "no caller yet" not in str(excinfo.value)

    def test_test_mode_still_reaches_an_engine(self, config, tmp_path):
        """The guard is a branch, not a replacement. TEST must still compute."""
        import shutil

        source = Path("test_data/intermediate/gwas/gwas_features.tsv")
        if not source.is_file():
            pytest.skip("gwas fixture not generated in this checkout")
        target = tmp_path / "gwas"
        target.mkdir(parents=True)
        shutil.copy(source, target / "gwas_features.tsv")
        results, gwas_input = stage_gwas.run(
            config, _manifest(), RunMode.TEST, tmp_path, _phenotypes()
        )
        assert isinstance(results, list)
        assert gwas_input is not None


class TestRealRunsWhenARealEngineIsInjected:
    """The other half of the conditional, which the old guard made unreachable.

    `PyseerEngine` exists and is complete. A REAL run that injects one must be
    allowed to use it; refusing it would be refusing a kinship-corrected mixed
    model because a guard could not tell one from a contingency table. What that
    run then needs from disk is the ordinary input contract, and the refusal
    for *that* names the file and its producer.
    """

    def test_a_real_engine_is_not_refused_out_of_hand(self, config, tmp_path):
        source = Path("test_data/intermediate/gwas/gwas_features.tsv")
        if not source.is_file():
            pytest.skip("gwas fixture not generated in this checkout")
        (tmp_path / "gwas").mkdir(parents=True)
        (tmp_path / "gwas" / "gwas_features.tsv").write_text(
            source.read_text(encoding="utf-8"), encoding="utf-8"
        )
        seen = []

        class SpyEngine(stage_gwas.GwasEngine):
            name = "spy"

            def run(self, gwas_input, cfg):
                seen.append(gwas_input.n_samples)
                return []

        results, gwas_input = stage_gwas.run(
            config, _manifest(), RunMode.REAL, tmp_path, _phenotypes(),
            engine=SpyEngine(),
        )
        assert results == []
        assert seen, "the injected engine never ran in REAL mode"

    def test_the_missing_input_names_the_file_and_its_producer(self, config, tmp_path):
        """An absent input is reported as an absent input, not a broken stage.

        `gwas_features.tsv` is stage 12's only on-disk input and this stage never
        writes it, so re-running stage 12 cannot fix it. A refusal that said only
        "TSV file not found" sent the reader looking for a bug in the reader.
        """
        reached = []

        class SpyEngine(stage_gwas.GwasEngine):
            name = "spy"

            def run(self, gwas_input, cfg):
                reached.append(True)
                return []

        with pytest.raises(DataContractError) as excinfo:
            stage_gwas.run(
                config, _manifest(), RunMode.REAL, tmp_path, _phenotypes(),
                engine=SpyEngine(),
            )
        message = str(excinfo.value)
        assert "gwas_features.tsv" in message, (
            "the refusal must name the missing input by filename"
        )
        assert "data_contract" in message, (
            "the refusal must name where the input's schema is declared, so a "
            "reader can find the producer rather than the reader"
        )
        assert not reached, "the engine ran without its feature table"


class TestTheSampleFloorNamesItsKey:
    def test_too_few_samples_names_the_config_key_and_the_cohort_floor(
        self, config, tmp_path
    ):
        """A count and a threshold with no key is a dead end for an operator.

        This is the refusal a real cohort meets first, and the key is the only
        actionable thing in it.
        """
        from papipeline.models import Phenotype, PhenotypeCall

        ids = _fixture_ids() or []
        calls = [
            PhenotypeCall(
                sample_id=i, antibiotic="imipenem",
                phenotype=Phenotype("R" if n < 2 else "S"), mic=None,
            )
            for n, i in enumerate(ids)
        ]
        with pytest.raises(DataContractError) as excinfo:
            build_input_stubbed(config, _manifest(), calls)
        assert "gwas.min_samples_per_group" in str(excinfo.value)

    def test_the_stated_cohort_floor_is_twice_the_group_minimum(self, config):
        """The minimum n is derived, not chosen.

        It is `2 * gwas.min_samples_per_group` because both groups must clear
        that floor independently. Stating it as a derived number keeps a single
        configured scientific threshold as the only thing anybody chose.
        """
        assert config.gwas.min_samples_per_group == 3
        assert 2 * config.gwas.min_samples_per_group == 6


def build_input_stubbed(config, manifest, calls):
    """`build_input` over rows that cannot clear either outcome group."""
    from papipeline.stages.gwas import build_input

    return build_input(
        config, manifest, calls,
        [{"sample_id": s, "gene__x": 1} for s in manifest.sample_ids],
    )


class TestTheGuardIsStructural:
    def test_run_py_injects_only_a_real_engine_or_none(self):
        """The reason a REAL run reached `ReferenceEngine` at all, restated.

        **R10, round 12 Phase 3.** This was
        `test_run_py_injects_no_engine`, and it asserted that the text
        ``stage_gwas.run(`` is not followed by ``engine=`` within 320
        characters. That assertion was a *pin on a gap*: it existed so that adding
        a REAL caller "means changing this deliberately", and it is now changed
        deliberately - `run.py` injects `PyseerEngine` in REAL.

        What the pin was really protecting is restated rather than dropped,
        because the protection is still needed and the old text no longer
        expresses it. Two things could go wrong at this call site and both are
        the failure this file exists for:

        1. the injected engine is a **real** one (`PyseerEngine`), or `None` so
           that `gwas.run`'s own guard refuses with its message naming pyseer; and
        2. nothing anywhere makes `run.py` construct a `ReferenceEngine` for a
           REAL run - that class is the mechanics-testing contingency table, and
           a REAL run reaching it would publish
           `reference_fisher:NO_KINSHIP_CORRECTION`, which is precisely what the
           refusal exists to prevent.

        So the assertion is behavioural and stronger than the text match it
        replaces: it calls `run._build_gwas_engine` in REAL with a resolvable
        pyseer and asserts the result is a `PyseerEngine`, and asserts
        `ReferenceEngine` appears nowhere in `run.py` at all.
        """
        import papipeline.run as run_module
        from papipeline.stages.gwas import PyseerEngine, ReferenceEngine

        class _FakeMachine:
            def tool_candidates(self, tool):
                return ("/opt/pyseer/bin/pyseer",)

        class _FakeConfig:
            machine = _FakeMachine()
            runtime = {"threads": 2}

            def phylogeny_dir(self, mode):
                from pathlib import Path
                return Path("/nonexistent/phylogeny")

        engine = run_module._build_gwas_engine(
            _FakeConfig(), RunMode.REAL, Path("/nonexistent/intermediate")
        )
        assert isinstance(engine, PyseerEngine), (
            f"run.py's REAL gwas branch built {type(engine).__name__!r}, not a "
            "PyseerEngine. Anything other than a real engine or None means a "
            "REAL run would reach ReferenceEngine - Fisher's exact with BH and "
            "no population-structure correction - which is the outcome the "
            "stage's refusal exists to prevent."
        )
        # Not a text search for the name: `_build_gwas_engine`'s own docstring
        # names `ReferenceEngine` to explain why it is never constructed, and a
        # docstring cannot be evidence about the code. What matters is whether
        # the name is REACHABLE in run.py's namespace - an import or an
        # assignment would make constructing it one edit away.
        assert not hasattr(run_module, "ReferenceEngine"), (
            "run.py binds the name ReferenceEngine. Constructing it for a REAL "
            "run would publish reference_fisher results stamped "
            "NO_KINSHIP_CORRECTION as a real cohort's associations, which is "
            "what the stage's refusal exists to prevent."
        )

    def test_no_engine_means_the_stage_refuses_rather_than_falling_back(self):
        """`None` is the safe answer when pyseer cannot be resolved.

        The other half of the pin above, and the reason `_build_gwas_engine`
        returns `None` instead of raising: `gwas.run`'s guard refuses a non-TEST
        mode whose engine is `None`, and its message NAMES pyseer. A helper that
        raised instead would replace a refusal that tells the operator which tool
        to install with one that does not.
        """
        import papipeline.run as run_module

        class _NoPyseer:
            def tool_candidates(self, tool):
                return ()

        class _FakeConfig:
            machine = _NoPyseer()
            runtime = {"threads": 2}

        assert run_module._build_gwas_engine(
            _FakeConfig(), RunMode.REAL, Path("/nonexistent")
        ) is None
        assert run_module._build_gwas_engine(
            _FakeConfig(), RunMode.TEST, Path("/nonexistent")
        ) is None, (
            "TEST must not get an engine either: spec.md D8 says TEST invokes no "
            "real tools, and the stage reads committed fixtures and supplies "
            "ReferenceEngine itself."
        )
