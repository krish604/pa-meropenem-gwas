"""An overlay must declare the tools the stages it runs actually invoke.

`tool_available` returns `None` for a tool the overlay never mentions, and
`None` is deliberately distinct from both `true` and `false`: it means the
config does not know. That distinction is the whole point of the block, and it
is what makes a missing entry a silent hole rather than a loud failure - a
stage whose tool is undeclared gets `None`, which reads as "not missing" to
anything testing it for falsiness.

The smoke overlay had exactly that hole. It declared `seqkit`, `bakta` and
`mlst`, and omitted `amrfinder` and `blast`, which is stage 4's entire
toolchain. Both are installed and both are declared `true` in the laptop and
bigmachine overlays, so the omission was not a statement about the machine - it
was a gap in the file, and it would have surfaced as an unexplained `None` at
the worst possible moment: mid-run, on the stage that blocks everything after
it.

So the invariant is asserted structurally. Every overlay declares stage 4's
toolchain, and declares it as a real boolean rather than leaving it implied.

A stronger check - that every `available: true` claim resolves to a real
binary - was written and removed, because it is unsound twice over. `blast` is
a suite name with no binary of that name (the binaries are `blastn`,
`makeblastdb`, and so on), so a config key is not a binary name. And
`bigmachine.yaml` describes a different machine, so verifying its claims
against this environment is meaningless. Neither gap is worth a test that
cannot be made true.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

MACHINES = Path("config/machines")

#: Stage 4 runs AMRFinderPlus, which in turn shells out to BLAST. Both are the
#: stage's toolchain: an AMRFinderPlus that cannot find blastn fails at
#: detection time, not at a stage boundary where the message can be useful.
STAGE4_TOOLCHAIN = ("amrfinder", "blast")

#: Stage 3's tool, asserted for the same reason as stage 4's.
STAGE3_TOOL = ("mlst",)


def _overlays() -> list:
    return sorted(p for p in MACHINES.glob("*.yaml"))


class TestEveryOverlayDeclaresTheToolsItsStagesInvoke:
    @pytest.mark.parametrize("path", _overlays(), ids=lambda p: p.stem)
    @pytest.mark.parametrize("tool", STAGE3_TOOL + STAGE4_TOOLCHAIN)
    def test_the_tool_is_declared(self, path: Path, tool: str):
        tools = yaml.safe_load(path.read_text(encoding="utf-8")).get("tools") or {}
        assert tool in tools, (
            f"{path.name} does not declare `{tool}`. An undeclared tool reads as "
            "unknown rather than missing, so a stage needing it gets None and "
            "fails mid-run instead of here."
        )

    @pytest.mark.parametrize("path", _overlays(), ids=lambda p: p.stem)
    @pytest.mark.parametrize("tool", STAGE3_TOOL + STAGE4_TOOLCHAIN)
    def test_the_declaration_is_a_boolean_not_a_bare_key(self, path: Path, tool: str):
        """`mlst:` with no value is None, which is the hole this file exists to
        close - it reads as unknown, not as false."""
        tools = yaml.safe_load(path.read_text(encoding="utf-8")).get("tools") or {}
        assert isinstance(tools[tool].get("available"), bool), (
            f"{path.name} declares `{tool}` without a boolean `available`, so it "
            "records nothing a stage can act on"
        )
