"""The `smoke` overlay may state one science key, and only that one.

**What this file is for.** ``annotation.reuse_tool_output`` decides whether
stage 2 runs Bakta or IMPORTS Bakta output already on disk, and stage 2 reads it
through ``config.reuse_tool_output()``. The overlay declares it so the mode is a
property of the run rather than of whatever ``science.yaml`` happens to say.

**What it asserts about the mode: that it is `require`.** `prefer` reuses
verified curated output and EXECUTES BAKTA when it is absent, so it leaves the
one tool this overlay exists not to run reachable. Round 12 set `require`,
reverted it to `prefer`, and recorded the reason as "`require` refuses the tool
when curated output is absent, so a smoke run refuses to annotate ALL 967 roster
members and produces nothing at all".

Both halves of that are wrong, and the reversion was the wrong fix for a real
cause. THE COHORT IS TEN: `config/machines/smoke.yaml:63` sets
`cohort.subset_file`, and `papipeline/run.py:501` applies it, so the manifest a
smoke run carries holds the ten isolates in `local/smoke_isolates.txt`. 967 is
`discover_pdc_manifest`'s unfiltered membership - the number
`papipeline/run.py:645` documents that `enforce_sample_cap` deliberately does
not measure on a smoke overlay. WHAT ACTUALLY BLOCKS `require` is that curated
Bakta output must be STAGED at `<intermediate_root>/bakta/<sample_id>/`, which
is where `run_bakta` would have written it
(`papipeline/stages/annotation.py:909`, `:690`). That is provisioning, and it
is fixed by provisioning.

What the reversion actually fixed was this file: the mode was asserted as a
LITERAL string, so `require` failed three tests that compared YAML text and
checked `!= "require"`. No stage-2 behaviour test failed at any point. So the
rest of this module pins the MECHANISM - the overlay may declare exactly this
key, the reader validates the value against the accepted set, and a bad value
is refused rather than coerced - and the one behavioural test is stated in both
directions in `test_smoke_require_never_reaches_bakta`.

**Why that needs an exemption at all.** ``annotation`` is a top-level section of
``config/science.yaml``, and ``load_machine_config`` refuses a machine overlay
that declares a science-owned section - the rule that makes it impossible for two
machines to disagree about the science. The exemption is
``OVERLAY_SCIENCE_EXEMPT_KEYS``, and it is **keyed, not section-wide**, so this
overlay can say ``reuse_tool_output`` and nothing else under ``annotation``.

The exemption is a hole in a guard, so it is tested from both sides: the one key
it permits works, and the keys it does not are still refused **by name**.
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path

import pytest
import yaml

from papipeline.config.loader import (
    OVERLAY_SCIENCE_EXEMPT_KEYS,
    REPO_ROOT,
    REUSE_TOOL_OUTPUT_MODES,
    ConfigError,
    load_config,
    load_machine_config,
)

SCIENCE = REPO_ROOT / "config" / "science.yaml"
SMOKE = REPO_ROOT / "config" / "machines" / "smoke.yaml"


def _overlay_with(tmp_path: Path, block: str) -> Path:
    """The committed smoke overlay, patched. Never edits the committed file."""
    source = SMOKE.read_text(encoding="utf-8")
    patched = source + "\n" + block
    assert patched != source
    out = tmp_path / "smoke-patched.yaml"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(patched, encoding="utf-8")
    return out


#: Matches the `reuse_tool_output` assignment in the committed overlay, whatever
#: its value is.
_REUSE_LINE = re.compile(
    r'^(?P<indent>\s*)reuse_tool_output:.*$', re.MULTILINE
)


def _with_reuse_mode(source: str, value: str) -> str:
    """The committed overlay text with `reuse_tool_output` set to `value`.

    Written as a PATTERN over the key rather than as a string replacement of the
    value, because the value is not what these tests are about. Round 12 wrote
    `.replace('reuse_tool_output: "prefer"', ...)` and it silently became a no-op
    the moment the committed mode changed: the intended typo was never injected,
    `pytest.raises` saw nothing, and the tests failed for a reason that had
    nothing to do with the reader's validation. A pattern over the key cannot
    rot that way - and it asserts the patch applied, which is the thing that
    would have caught it immediately.
    """
    patched, count = _REUSE_LINE.subn(
        lambda m: f'{m.group("indent")}reuse_tool_output: {value}', source
    )
    assert count == 1, (
        f"expected exactly one reuse_tool_output assignment in the committed "
        f"smoke overlay, found {count}"
    )
    return patched


class TestTheOneExemptKeyWorks:
    def test_the_committed_smoke_overlay_says_which_mode_it_is(self):
        """Read through the public accessor, not by text-searching the YAML.

        The accessor is what stage 2 reads (``annotation.py:906``), so asserting
        on it makes the claim about behaviour rather than about a YAML file. The
        assertion is that the overlay DECLARES a mode and the reader returns it
        verbatim - which is what "the overlay may state this key" means - not
        that the mode is any particular one.
        """
        raw = yaml.safe_load(SMOKE.read_text(encoding="utf-8"))
        declared = (raw.get("annotation") or {}).get("reuse_tool_output")
        assert declared is not None, (
            "the smoke overlay no longer declares annotation.reuse_tool_output, "
            "so stage 2 falls back to the science.yaml default and the overlay "
            "cannot control reuse at all"
        )

        config = load_config(SCIENCE, machine=SMOKE)
        assert config.reuse_tool_output() == declared, (
            f"the overlay declares {declared!r} but the accessor returns "
            f"{config.reuse_tool_output()!r}, so stage 2 is not reading what "
            "the overlay says"
        )

    def test_the_committed_mode_is_require_so_bakta_is_never_reachable(self):
        """The mode, stated directly, with the reason it is not `prefer`.

        Round 12 set this to `require`, reverted it to `prefer` on a diagnosis
        that turned out to be wrong, and the reason it was wrong is the point
        of this test: the reversion was made because three tests in this file
        compared the YAML against the LITERAL string `"prefer"` or asserted
        `!= "require"`. No test of stage 2's behaviour failed at any point. The
        claim that `require` "refuses ALL 967 roster members" does not describe
        this machine: the smoke cohort is the TEN isolates in
        `local/smoke_isolates.txt` (`config/machines/smoke.yaml:63`,
        applied at `papipeline/run.py:501`), and 967 is the unfiltered
        manifest membership that `papipeline/run.py:645` documents the cap does
        not measure.

        So the mode is `require`, and it is asserted on the value rather than
        on its spelling in the file, because the value is what stage 2 reads.
        """
        config = load_config(SCIENCE, machine=SMOKE)
        assert config.reuse_tool_output() == "require", (
            "the smoke overlay must set reuse_tool_output: require. `prefer` "
            "reuses verified curated output and EXECUTES BAKTA when it is "
            "absent, which leaves the one tool this overlay exists not to run "
            "reachable. `require` never invokes it: an isolate whose curated "
            "output is absent is recorded as failed with its refusal named "
            "(`papipeline/stages/annotation.py:942-950`). Curated output is "
            "staged at <intermediate_root>/bakta/<sample_id>/, which is where "
            "`run_bakta` would have written it."
        )

    def test_require_is_the_only_mode_that_cannot_reach_the_tool(self):
        """Why the assertion above is not a preference between two strings.

        `off` and `prefer` both reach `run_bakta` on some path; `require` is
        the only one of the three whose control flow has no edge to it. Asserted
        against the stage's own source rather than its docstring, because the
        docstring is a claim and the branch is the mechanism.
        """
        import inspect

        from papipeline.stages import annotation as stage

        source = inspect.getsource(stage.run)
        # The tool is reached from the one `run_bakta` call, and it sits BELOW
        # the `require` refusal, which `continue`s. So the requirement is that
        # the refusal appears before the call and returns control first.
        refusal = source.index('if reuse_mode == "require":')
        invocation = source.index("run_bakta(")
        assert refusal < invocation, (
            "the require refusal must come before the tool call, or a genome "
            "whose curated output is absent would reach the tool"
        )
        between = source[refusal:invocation]
        assert between.count("continue") >= 1, (
            "the require branch must `continue` before the run_bakta call is "
            "reached; found no control transfer between them"
        )

    def test_the_key_is_the_one_named_in_the_exemption(self):
        """The file and the code cannot disagree about which key it is."""
        raw = yaml.safe_load(SMOKE.read_text(encoding="utf-8"))
        declared = set(raw.get("annotation") or {})
        permitted = OVERLAY_SCIENCE_EXEMPT_KEYS["smoke"]["annotation"]
        assert declared <= permitted, (
            f"config/machines/smoke.yaml declares {sorted(declared)} under "
            f"`annotation`, which is more than the exemption permits "
            f"({sorted(permitted)}). load_machine_config would refuse this "
            "overlay at load; the exemption table and the file must move "
            "together."
        )
        assert declared, (
            "the smoke overlay declares nothing under `annotation`, so stage 2 "
            "would take the science.yaml default and the overlay would have no "
            "say over reuse at all"
        )

    def test_it_is_a_recognised_mode_not_a_typo(self):
        config = load_config(SCIENCE, machine=SMOKE)
        assert config.reuse_tool_output() in REUSE_TOOL_OUTPUT_MODES

    def test_only_the_smoke_overlay_has_the_exemption(self):
        """A second overlay must not acquire it by being added.

        `cohort.subset_file` has the same constraint and the same reason: the
        exemption names the overlay rather than inferring it from a path.
        """
        assert set(OVERLAY_SCIENCE_EXEMPT_KEYS) == {"smoke"}, (
            "an overlay other than `smoke` may now declare a science key. Each "
            "exemption is a deliberate decision about one machine; widen this "
            "only with the reason written down beside it."
        )


class TestTheExemptionIsKeyedNotSectionWide:
    def test_another_key_under_annotation_is_still_refused(self, tmp_path):
        """The guard the exemption must not have opened.

        ``annotation`` is a science section with more than one key, so
        exempting the section would let this overlay set any of them. The
        refusal names the offending key so a reader knows which one to move.
        """
        overlay = _overlay_with(tmp_path, "annotation:\n  something_else: 1\n")
        with pytest.raises(ConfigError) as excinfo:
            load_machine_config(overlay)
        message = str(excinfo.value)
        assert "annotation.something_else" in message, (
            f"the refusal does not name the offending key: {message}"
        )
        assert "science.yaml" in message, (
            "the refusal does not say where the key belongs"
        )

    def test_the_exemption_is_not_a_general_purpose_door(self, tmp_path):
        """A *different* science section is still refused outright."""
        overlay = _overlay_with(tmp_path, "gwas:\n  enabled: false\n")
        with pytest.raises(ConfigError) as excinfo:
            load_machine_config(overlay)
        assert "gwas" in str(excinfo.value)

    def test_a_bad_value_is_refused_by_the_reader_not_accepted(self, tmp_path):
        """A typo must not resolve to a silent default.

        ``reuse`` or ``require: true`` resolving to ``off`` would re-annotate the
        whole cohort and report nothing about why, which reads exactly like a run
        that was configured correctly. The loader refuses; this asserts the
        refusal reaches the caller through the overlay path too.

        The accessor is CALLED inside the ``raises`` block, not merely loaded:
        ``reuse_tool_output()`` validates at read time, so ``load_config`` alone
        would succeed and this test would pass for the wrong reason - proving
        nothing about the refusal.
        """
        source = SMOKE.read_text(encoding="utf-8")
        text = _with_reuse_mode(source, "reuse")
        assert text != source, (
            "could not patch reuse_tool_output out of the committed overlay; "
            "the spelling of the key in config/machines/smoke.yaml has changed"
        )
        out = tmp_path / "smoke-typo.yaml"
        out.write_text(text, encoding="utf-8")
        with pytest.raises(ConfigError) as excinfo:
            load_config(SCIENCE, machine=out).reuse_tool_output()
        assert "reuse_tool_output" in str(excinfo.value)

    def test_an_unquoted_off_is_refused_rather_than_read_as_false(self, tmp_path):
        """YAML 1.1 reads a bare `off` as the boolean ``false``.

        The loader refuses a boolean by name rather than guessing, so the natural
        spelling of the default produces an error that says what is wrong instead
        of one that reads as though the pipeline had invented a rule.
        """
        source = SMOKE.read_text(encoding="utf-8")
        text = _with_reuse_mode(source, "off")
        assert text != source, (
            "could not patch reuse_tool_output out of the committed overlay"
        )
        out = tmp_path / "smoke-bare-off.yaml"
        out.write_text(text, encoding="utf-8")
        with pytest.raises(ConfigError) as excinfo:
            load_config(SCIENCE, machine=out).reuse_tool_output()
        assert "boolean" in str(excinfo.value)


class TestSilenceIsNotAChoice:
    def test_the_other_overlays_inherit_the_science_default(self):
        """`None` on MachineConfig, not `off`.

        An overlay that says nothing must be distinguishable from one that chose
        ``off``: the first inherits science.yaml, the second overrode it. Had the
        field defaulted to ``off`` those would be the same value and the overlay's
        silence would become a decision nobody made.
        """
        for machine in ("laptop", "bigmachine"):
            overlay = load_machine_config(machine)
            assert overlay.reuse_tool_output is None, (
                f"{machine} states annotation.reuse_tool_output "
                f"({overlay.reuse_tool_output!r}) without declaring it; the field "
                "must be None so 'said nothing' is distinguishable from 'said off'"
            )
            config = load_config(SCIENCE, machine=machine)
            assert config.reuse_tool_output() == "off", (
                f"{machine} no longer inherits the science default"
            )


@pytest.mark.skipif(shutil.which("python3") is None, reason="no interpreter")
def test_the_overlay_is_committed_and_machine_is_known():
    """`--machine smoke` resolves by name.

    D12: the name did not resolve, so every invocation passed a PATH that the
    Snakefile interpolates verbatim into a shell command - resolved against the
    CWD rather than the worktree root.
    """
    from papipeline.config.loader import KNOWN_MACHINES

    assert "smoke" in KNOWN_MACHINES
    config = load_config(SCIENCE, machine="smoke")
    assert config.machine is not None
    assert config.machine.name == "smoke"