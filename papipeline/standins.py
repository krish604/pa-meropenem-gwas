"""Stand-in fixtures, and the only sanctioned way to read one.

A stand-in fixture is a one-row table under
``test_data/standins/unbuilt_stages/``, headed ``STAND-IN FIXTURE``, whose
values say plainly that nothing computed them. The DAG declares the edges the
spec requires, so a stage that is *declared* but not *built* needs something to
satisfy them in TEST, or the workflow cannot resolve.

**The tuple is empty.** Every stage in `STAGE_ORDER` is now built and
dispatched: `variants` and `cohort_variants` in f5573f5, `similarity` in
2223332, and `recombination` in the `dag-resolve` reconciliation. Each has a
TEST path that reads committed *real* output, so none of them needs a stand-in.
The mechanism is kept because it is the right answer for the next stage that is
declared and not built, and because `load_stage_table` is the sanctioned reader
whether or not anything calls it today.

The risk a stand-in carries is not that it exists. It is that it is read as
though it were output, because a header-only contract table with plausible
columns is exactly what a real result looks like. So:

* the fixtures declare themselves, in a comment header and in every value that
  could be mistaken for a call;
* :func:`load_stage_table` is the reader for these tables and it refuses;
* a test asserts no other module reads them.

The committed fixture *files* under ``test_data/standins/unbuilt_stages/`` are
deliberately left on disk. They are unreferenced, which is harmless, and the
stand-in test asserts a file exists for every stage the tuple lists - so
deleting them has to happen together with the tuple, not before it. See
TASKS.md.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .errors import DataContractError
from .io.tsv import read_tsv

#: Banner every stand-in fixture carries in its comment header. Checked first,
#: because it survives even if the rows are edited.
STANDIN_BANNER = "STAND-IN FIXTURE"

#: Values a stand-in puts in the cell a real row would put a call in. A
#: substring test, so `standin_not_computed_2` is caught too.
STANDIN_VALUE_PREFIX = "standin"

#: The stages whose output is stood in for. Named here so the loader can
#: refuse by stage name and not only by content.
#:
#: `variants` left this tuple in f5573f5, `similarity` in 2223332, and
#: `recombination` in this reconciliation - the dispatch branch that promoted it
#: is `run.derive_recombination_tables`, and TEST now reads committed
#: gubbins-shaped output from `test_data/intermediate/gubbins/`. Neither is
#: stood in for now: all three have a TEST path and a REAL caller, and
#: `test_data/` holds real output rather than a stand-in. A stage leaves this
#: tuple when it stops being absent - not when its stand-in merely stops being
#: read.
UNBUILT_WITH_STANDIN: Tuple[str, ...] = ()


class StandInError(DataContractError):
    """A stand-in fixture was offered where real output was required.

    A subclass of DataContractError because that is what it is: a file that
    does not satisfy the contract it is being read against. The run fails
    rather than proceeding, which is the point - a stand-in silently consumed
    is a fabricated result with a plausible header.
    """


def is_standin_table(path: Path) -> bool:
    """Whether ``path`` is a stand-in fixture rather than stage output.

    Two independent signals, either of which is enough:

    * the comment header carries :data:`STANDIN_BANNER`;
    * any cell starts with :data:`STANDIN_VALUE_PREFIX`.

    The header is checked separately from the values because a stand-in whose
    rows were edited by hand keeps its banner.
    """
    path = Path(path)
    if not path.exists():
        return False

    text = path.read_text(encoding="utf-8")
    for line in text.splitlines():
        if not line.startswith("#"):
            break
        if STANDIN_BANNER in line:
            return True

    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        for cell in line.split("\t"):
            if cell.strip().lower().startswith(STANDIN_VALUE_PREFIX):
                return True
    return False


def load_stage_table(
    path: Path,
    stage: str,
    *,
    required_columns: Sequence[str] = (),
) -> List[Dict[str, Optional[str]]]:
    """Read a table, refusing if it turns out to be a stand-in.

    Args:
        path: The table to read.
        stage: Stage that is supposed to have written it, for the error message.
        required_columns: Forwarded to :func:`read_tsv`.

    Returns:
        The rows, when the file is real output.

    Raises:
        StandInError: The file is a stand-in fixture. This is not a fallback -
            a caller that wanted stage output has been handed a placeholder, and
            continuing would fabricate a result.
    """
    path = Path(path)
    if is_standin_table(path):
        raise StandInError(
            f"Refusing to read a stand-in fixture as {stage!r} output",
            path=str(path),
            stage=stage,
            hint=(
                f"{path} carries the {STANDIN_BANNER} banner: it stands in for an "
                "unbuilt stage (ticket 14), contains no computed values, and "
                "must not be read as a result."
            ),
        )
    return read_tsv(path, required_columns=required_columns)
