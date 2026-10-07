"""The ``unitig-caller`` adapter: stub in TEST/STUB, refusal in REAL.

AGENTS.md rule 1 says tool flags are verified with ``--help`` rather than
assumed, and ``unitig-caller`` is not installed on the machine this package was
written against, so no flag of that tool has ever been run here. The adapter
therefore records that honestly instead of inventing a command line:

* :data:`FLAGS_VERIFIED` is ``False`` and :data:`VERIFIED_FLAGS` is empty -
  nothing is claimed as checked;
* :func:`unitig_tool_status` reports presence without executing anything
  (``adapters.external.detect_tools`` runs no binary by default);
* TEST and STUB runs write a stub table carrying the isolate set, which proves
  the interface and the file contract without pretending unitigs were called;
* a REAL run raises :class:`UnverifiedFlagsError` naming the command that was
  never run, and writes nothing.

When ``unitig-caller --help`` has actually been run and its flags recorded,
this is the one module to change: set :data:`VERIFIED_FLAGS`, then give the
REAL branch a command built from those flags and from configuration.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence, Tuple

from ..adapters.external import ToolStatus, detect_tools
from ..errors import PipelineError
from ..io.tsv import write_tsv
from ..logging_utils import get_logger
from ..models import RunMode

LOGGER = get_logger("layers.unitigs")

#: The tool this adapter wraps, and where its output lives.
UNITIG_TOOL = "unitig-caller"
UNITIG_TABLE_FILENAME = "unitigs.tsv"

#: The verified flag set. Empty because the tool is absent, and an empty set
#: is the point: no flag is claimed as checked.
FLAGS_VERIFIED = False
VERIFIED_FLAGS: Tuple[str, ...] = ()

#: Why REAL refuses. Records exactly what was not run, so the next reader
#: knows which command to verify rather than having to guess.
UNVERIFIED_REASON = (
    f"{UNITIG_TOOL} flags are not verified: `{UNITIG_TOOL} --help` has not "
    "been run in this environment - the tool was absent from PATH when this "
    "adapter was written - so no flag of it is known, and AGENTS.md rule 1 "
    "forbids inventing one. REAL mode therefore refuses to invoke the tool."
)


class UnverifiedFlagsError(PipelineError):
    """Raised instead of running ``unitig-caller`` with guessed flags."""


def unitig_tool_status() -> ToolStatus:
    """Presence of ``unitig-caller`` on PATH, without executing it."""
    return detect_tools([UNITIG_TOOL])[UNITIG_TOOL]


def run_unitig_screen(
    *, mode: RunMode, sample_ids: Sequence[str], output_dir: Path
) -> Path:
    """Write the unitig table for this mode.

    TEST and STUB write a stub: one row per isolate, no unitig features, and a
    banner that says the tool was not invoked and which command was never run.
    REAL refuses, and refuses *before* writing anything, because an artefact
    left behind by a refusal reads as a screen that ran.

    Raises:
        UnverifiedFlagsError: In REAL mode, while the flags are unverified.
    """
    if mode is RunMode.REAL:
        status = unitig_tool_status()
        raise UnverifiedFlagsError(
            f"{UNVERIFIED_REASON} Detected on PATH right now: "
            f"{status.available}."
        )

    path = Path(output_dir) / UNITIG_TABLE_FILENAME
    banner = [
        f"tool: {UNITIG_TOOL}",
        "stub: true",
        f"mode: {mode.value}",
        f"note: {UNITIG_TOOL} was not invoked; `{UNITIG_TOOL} --help` has not "
        "been run in this environment, so no flag of it is verified and no "
        "unitig feature is produced",
    ]
    write_tsv(
        path,
        [{"sample_id": sample_id} for sample_id in sample_ids],
        ["sample_id"],
        header_comment=banner,
    )
    LOGGER.info(
        "layers: wrote stub %s for %d isolates; %s was not invoked",
        path,
        len(sample_ids),
        UNITIG_TOOL,
    )
    return path
