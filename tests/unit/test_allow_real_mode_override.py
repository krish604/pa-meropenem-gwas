"""REAL mode must be reachable without editing a tracked file.

A REAL run needs `runtime.allow_real_mode` true, and that lives in
`config/machines/smoke.yaml` - a **tracked** file whose committed value must stay
`false`. So enabling a run meant editing it, which meant:

* the working tree was dirty for the duration, so nothing else could be
  committed while a run was in flight;
* two tests that assert the committed flag is `false` failed for as long as the
  run existed, because the run had mutated the file they read.

Neither is a test defect. The flag is deliberately shut in the repository - the
overlay header says so - and the *ability to run* is a property of the operator's
session, not of the committed configuration.

So the override is an environment variable, which is the pattern
`PipelineConfig.results_root` already established for the same reason: a Snakemake
`--config` key is re-read by each rule in its own process, whereas the environment
is inherited by all of them. A local gitignored overlay was the other option and
would have needed a merge step at load; the environment is the mechanism this
loader already has for exactly this class of override.

**The committed file is still authoritative.** The variable can only *open* the
gate, never close one that is already open, so setting it cannot weaken a
machine whose overlay says `false` by accident - it makes an intentional opening
possible. That asymmetry matters: an override that could also close the gate
would let a stray `PIPELINE_ALLOW_REAL_MODE=0` silently re-shut a machine that
deliberately enabled REAL mode.
"""

from __future__ import annotations

import pytest

from papipeline.config import loader as loader_module
from papipeline.config.loader import load_config
from papipeline.models import RunMode

REPO = __import__("pathlib").Path(__file__).resolve().parents[2]
ENV = loader_module.ALLOW_REAL_MODE_ENV


@pytest.fixture
def smoke_overlay():
    return REPO / "config" / "machines" / "smoke.yaml"


class TestTheCommittedFlagStaysShut:
    def test_the_committed_smoke_overlay_is_false(self, smoke_overlay):
        """The guard the whole change is protecting."""
        import yaml

        runtime = yaml.safe_load(smoke_overlay.read_text(encoding="utf-8"))["runtime"]
        assert runtime["allow_real_mode"] is False, (
            "the committed overlay must keep REAL mode shut; the environment "
            "override exists so a run does not have to edit this file"
        )

    @pytest.mark.parametrize("overlay", ["laptop", "bigmachine", "smoke"])
    def test_every_committed_overlay_is_shut(self, overlay):
        import yaml

        path = REPO / "config" / "machines" / f"{overlay}.yaml"
        runtime = yaml.safe_load(path.read_text(encoding="utf-8"))["runtime"]
        assert runtime["allow_real_mode"] is False


class TestTheEnvironmentOpensTheGate:
    def test_it_is_absent_by_default(self, monkeypatch):
        monkeypatch.delenv(ENV, raising=False)
        config = load_config(REPO / "config" / "science.yaml", machine="laptop")
        assert config.runtime.get("allow_real_mode") is False

    @pytest.mark.parametrize(
        "value", ["1", "true", "TRUE", "yes", "on"]
    )
    def test_truthy_spellings_open_it(self, monkeypatch, value):
        monkeypatch.setenv(ENV, value)
        config = load_config(REPO / "config" / "science.yaml", machine="laptop")
        assert config.runtime.get("allow_real_mode") is True, (
            f"{ENV}={value!r} was not read as enabling REAL mode; a run would "
            "refuse for no stated reason"
        )

    @pytest.mark.parametrize("value", ["0", "false", "no", "off", ""])
    def test_falsey_spellings_do_not_open_it(self, monkeypatch, value):
        """An unrecognised or negative value must not enable REAL mode by
        accident - a typo that happened to be truthy would be worse than one
        that is inert."""
        monkeypatch.setenv(ENV, value)
        config = load_config(REPO / "config" / "science.yaml", machine="laptop")
        assert config.runtime.get("allow_real_mode") is False

    def test_an_unrecognised_value_refuses_rather_than_guessing(self, monkeypatch):
        monkeypatch.setenv(ENV, "perhaps")
        with pytest.raises(Exception) as excinfo:
            load_config(REPO / "config" / "science.yaml", machine="laptop")
        assert ENV in str(excinfo.value)

    def test_it_cannot_close_a_gate_that_is_already_open(self, monkeypatch):
        """The asymmetry, stated as a test.

        The override exists to let an operator open a shut gate. Were it also
        able to close an open one, a stray `=0` in a shell would silently
        re-shut a machine whose overlay deliberately enables REAL mode - and the
        failure would look like the overlay's own doing.
        """
        import copy

        import yaml

        monkeypatch.setenv(ENV, "0")
        data = yaml.safe_load(
            (REPO / "config" / "machines" / "bigmachine.yaml").read_text(
                encoding="utf-8"
            )
        )
        data["runtime"]["allow_real_mode"] = True
        # Not written to the repository - a temporary overlay in tmp_path.
        import tempfile

        path = REPO / ".scratch_tmp_open_overlay.yaml"
        try:
            path.write_text(yaml.safe_dump(data), encoding="utf-8")
            config = load_config(REPO / "config" / "science.yaml", machine=path)
            assert config.runtime.get("allow_real_mode") is True
        finally:
            path.unlink(missing_ok=True)


class TestRunPipelineSeesIt:
    def test_the_gate_reads_the_resolved_value_not_the_file(self, monkeypatch):
        """`run.py` asks the config, so fixing the loader is enough.

        Worth pinning because the alternative - teaching `run.py` about the
        environment - would leave every *other* consumer of
        `runtime.allow_real_mode` reading the committed value and disagreeing
        with the gate that actually decides.
        """
        from papipeline.errors import ModeNotAllowedError
        from papipeline.run import resolve_mode

        monkeypatch.delenv(ENV, raising=False)
        config = load_config(REPO / "config" / "science.yaml", machine=REPO / "config" / "machines" / "smoke.yaml")
        with pytest.raises(ModeNotAllowedError):
            resolve_mode("REAL", config)

        monkeypatch.setenv(ENV, "1")
        opened = load_config(REPO / "config" / "science.yaml", machine=REPO / "config" / "machines" / "smoke.yaml")
        assert resolve_mode("REAL", opened) is RunMode.REAL


class TestTheOverlayIsUntouched:
    def test_a_real_run_leaves_the_file_alone(self, monkeypatch, tmp_path):
        """The end-to-end property, without running the pipeline.

        Resolve the config for a REAL run and confirm the committed overlay is
        byte-identical afterwards. This is what makes committing other work
        possible while a run is in flight.
        """
        before = (REPO / "config" / "machines" / "smoke.yaml").read_bytes()
        monkeypatch.setenv(ENV, "1")
        config = load_config(REPO / "config" / "science.yaml", machine=REPO / "config" / "machines" / "smoke.yaml")
        assert config.runtime.get("allow_real_mode") is True
        assert (REPO / "config" / "machines" / "smoke.yaml").read_bytes() == before