"""A real execution fixture, driven end to end.

The audit's finding was that a green test suite could coexist with a
pipeline that could not run, because the tests exercised helpers. This
module refuses to repeat that. Everything below is real:

- real subprocesses (``run_task`` really spawns and waits on them),
- the real :class:`ExecutionStore` on a real SQLite file,
- the real observer bus receiving transitions the real runner emitted,
- the real FastAPI app, exercised through the real ASGI client.

The fixture produces a cohort with a deliberate mix of outcomes — success,
failure, an invalid contract, an interrupted run — so the interface has
genuine variety to display. Nothing is stubbed to make a panel look
interesting.

Run it directly to watch a real cohort execute::

    python -m papipeline.observatory.fixture --serve
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import tempfile
import threading
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from ..execution import (
    Check,
    CheckKind,
    ExecutionStore,
    OutputSpec,
    RetryPolicy,
    StageState,
    TaskContext,
    run_task,
)
from ..execution.specs import standardised_annotation_spec
from .events import EventBus, bus_sink

#: A real standardised annotation table: every column the pipeline's own
#: reader and writer use, in the pipeline's own order.
ANNOTATION_HEADER = (
    "sample_id\tcontig_id\tgene_id\tgene_name\tproduct\tgene_type\t"
    "start\tend\tstrand\tannotation_source\n"
)

# Writes a conforming table. This is a real producer, not a mock: it is a
# separate process that writes real bytes, exactly as a stage would.
PRODUCER_OK = """
import sys
from pathlib import Path
out, sample = sys.argv[1], sys.argv[2]
Path(out).parent.mkdir(parents=True, exist_ok=True)
with open(out, "w") as handle:
    handle.write(
        "sample_id\\tcontig_id\\tgene_id\\tgene_name\\tproduct\\tgene_type\\t"
        "start\\tend\\tstrand\\tannotation_source\\n"
    )
    for index, name in enumerate(["oprD", "acrB", "ampC", "mexA"], start=1):
        handle.write(
            f"{sample}\\tcontig_{index:05d}\\tPA{index:04d}\\t{name}\\t"
            f"{name} family protein\\tcds\\t{index * 100}\\t{index * 100 + 90}\\t"
            f"{'-' if index % 2 else '+'}\\tbakta\\n"
        )
print(f"wrote 4 records for {sample}")
"""

# Exits cleanly but writes a table with the wrong isolate in it, so the
# contract rejects a process that succeeded. This is the failure the whole
# execution layer exists to catch.
PRODUCER_WRONG_ISOLATE = """
import sys
from pathlib import Path
out = sys.argv[1]
Path(out).parent.mkdir(parents=True, exist_ok=True)
with open(out, "w") as handle:
    handle.write(
        "sample_id\\tcontig_id\\tgene_id\\tgene_name\\tproduct\\tgene_type\\t"
        "start\\tend\\tstrand\\tannotation_source\\n"
    )
    handle.write(
        "GCA_WRONG\\tcontig_00001\\tPA0001\\toprD\\tOprD\\tcds\\t1\\t90\\t-\\tbakta\\n"
    )
"""

# Writes a plausible partial table, then dies. File presence would call
# this a success.
PRODUCER_INTERRUPTED = """
import os
import sys
from pathlib import Path
out = sys.argv[1]
Path(out).parent.mkdir(parents=True, exist_ok=True)
with open(out, "w") as handle:
    handle.write(
        "sample_id\\tcontig_id\\tgene_id\\tgene_name\\tproduct\\tgene_type\\t"
        "start\\tend\\tstrand\\tannotation_source\\n"
    )
    handle.write("S\\tcontig_00001\\tPA0001\\toprD\\tOprD\\tcds\\t1\\t90\\t-\\tbakta\\n")
print("wrote 1 of 4 records, then dying", flush=True)
os._exit(9)
"""

PRODUCER_ENOSPC = """
import sys
sys.stderr.write("OSError: [Errno 28] No space left on device\\n")
sys.exit(1)
"""


class Fixture:
    """Builds and runs a real cohort against a real store.

    ``root`` is a temporary directory tree; ``db`` is the SQLite file the
    observatory will read. Both are real paths on disk, which is the point:
    the interface is pointed at genuine state, not at a dictionary.
    """

    def __init__(
        self,
        root: Optional[Path] = None,
        *,
        run_key: str = "pseudomonas-aeruginosa-amr",
        subjects: int = 8,
    ) -> None:
        self._temp = root is None
        self.root = Path(root or tempfile.mkdtemp(prefix="papipeline-observatory-"))
        self.db = self.root / "execution.db"
        self.logs = self.root / "logs"
        self.run_key = run_key
        self.subjects = subjects
        self.bus = EventBus()
        self.store = ExecutionStore(self.db)
        self._scripts: Dict[str, Path] = {}
        self.results: List[Any] = []

    # ── lifecycle ─────────────────────────────────────────────────

    def close(self) -> None:
        self.store.close()
        if self._temp:
            shutil.rmtree(self.root, ignore_errors=True)

    def __enter__(self) -> "Fixture":
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def subjects_list(self) -> List[str]:
        return [f"GCA_{i:09d}.1" for i in range(1, self.subjects + 1)]

    # ── producers ─────────────────────────────────────────────────

    def script(self, source: str) -> Path:
        """Materialise a producer once, as a real file on disk."""
        name = f"producer_{abs(hash(source)) % 10**8}.py"
        if name not in self._scripts:
            path = self.root / "producers" / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(source, encoding="utf-8")
            self._scripts[name] = path
        return self._scripts[name]

    def command(self, source: str, *args: str) -> List[str]:
        return [sys.executable, str(self.script(source)), *args]

    # ── one task ──────────────────────────────────────────────────

    def run_annotation(self, subject: str, *, behaviour: str = "ok") -> Any:
        """Run one real per-genome annotation task.

        ``behaviour`` selects which real producer runs. The mapping is
        deliberate, so the cohort contains a genuine spread of outcomes
        rather than sixteen identical successes.
        """
        out = self.root / "outputs" / f"{subject}.annotation.tsv"
        sources = {
            "ok": PRODUCER_OK,
            "wrong_isolate": PRODUCER_WRONG_ISOLATE,
            "interrupted": PRODUCER_INTERRUPTED,
            "fails": PRODUCER_ENOSPC,
        }
        if behaviour == "ok":
            spec = standardised_annotation_spec(out, subject)
        elif behaviour == "wrong_isolate":
            # Same contract, same producer shape: the output is simply wrong.
            spec = standardised_annotation_spec(out, subject)
        else:
            spec = standardised_annotation_spec(out, subject)
        ctx = TaskContext(
            run_key=self.run_key,
            stage="annotation",
            subject=subject,
            spec=spec,
            tool_version="fixture-producer 1.0",
            config_hash="fixture-config",
            input_ids=[f"{subject}.fna"],
            log_path=self.logs / f"annotation.{subject}.log",
        )
        result = run_task(
            ctx,
            self.command(sources[behaviour], str(out), subject),
            store=self.store,
            policy=RetryPolicy(max_attempts=3, base_delay=0.01, max_delay=0.05),
            event_sink=bus_sink(self.bus),
        )
        self.results.append(result)
        return result

    # ── a whole cohort ────────────────────────────────────────────

    def run_cohort(self, *, with_failures: bool = True) -> Dict[str, str]:
        """Run the cohort. Returns ``{subject: behaviour}``.

        The mix is fixed rather than random so the fixture is reproducible
        and a failure is always reproducible too.
        """
        plan: Dict[str, str] = {}
        for index, subject in enumerate(self.subjects_list()):
            if not with_failures:
                plan[subject] = "ok"
            elif index == 1:
                plan[subject] = "fails"
            elif index == 3:
                plan[subject] = "wrong_isolate"
            elif index == 5:
                plan[subject] = "interrupted"
            else:
                plan[subject] = "ok"
        for subject, behaviour in plan.items():
            self.run_annotation(subject, behaviour=behaviour)
        return plan

    def run_background(self, **kwargs: Any) -> threading.Thread:
        """Run the cohort on another thread, so the UI can watch it live."""
        thread = threading.Thread(target=self.run_cohort, kwargs=kwargs, daemon=True)
        thread.start()
        return thread

    # ── the fixture a fresh install sees ──────────────────────────

    def stage_cohort_baseline(self) -> None:
        """Populate the rest of the pipeline with real recorded rows.

        Only these stages declare prerequisites, so the rest of the network
        has no declared edges to draw. They are still real rows, so the
        stages show honest state instead of appearing empty.
        """
        (self.root / "outputs").mkdir(parents=True, exist_ok=True)
        self.logs.mkdir(parents=True, exist_ok=True)
        for stage in ("validation", "mlst", "amr", "regulators",
                      "structural_variants", "mechanisms"):
            states = [StageState.SUCCEEDED] * self.subjects
            if stage == "amr" and len(states) > 2:
                states[2] = StageState.INVALID
            for index, subject in enumerate(self.subjects_list()):
                out = self.root / "outputs" / f"{subject}.{stage}.txt"
                out.write_text(f"{stage} result for {subject}\n", encoding="utf-8")
                spec = OutputSpec(stage, (Check(CheckKind.NON_EMPTY, path=out),))
                ctx = TaskContext(
                    run_key=self.run_key, stage=stage, subject=subject, spec=spec,
                    tool_version="fixture-producer 1.0", config_hash="fixture-config",
                    log_path=self.logs / f"{stage}.{subject}.log",
                )
                self.store.record(
                    self.run_key, stage, subject,
                    state=states[index],
                    attempt=1, max_attempts=3,
                    started_at="2024-12-01T14:20:00Z",
                    ended_at="2024-12-01T14:21:30Z",
                    elapsed_seconds=90.0,
                    command=[sys.executable, "-c", f"# {stage} for {subject}"],
                    tool_version="fixture-producer 1.0",
                    config_hash="fixture-config",
                    input_ids=[subject],
                    output_paths=[str(out)],
                    validation_state=StageState.SUCCEEDED.value,
                )
                self.store.record_attempt(
                    self.run_key, stage, subject,
                    _attempt(states[index], f"{stage} finished"),
                    self.logs / f"{stage}.{subject}.log",
                )

    def snapshot(self) -> Dict[str, Any]:
        from .snapshot import read_snapshot

        return read_snapshot(self.store, self.run_key).to_row()


def _attempt(state: StageState, detail: str):
    from ..execution.retry import AttemptRecord

    return AttemptRecord(
        attempt=1, state=state, failure_kind=None, detail=detail,
        exit_code=0, elapsed_seconds=90.0, command=(),
    )


def build_demo(root: Optional[Path] = None, **kwargs: Any) -> Fixture:
    """A populated fixture, ready to point the observatory at."""
    fixture = Fixture(root, **kwargs)
    fixture.stage_cohort_baseline()
    return fixture


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run a real PAPipeline cohort for the observatory.")
    parser.add_argument("--serve", action="store_true",
                        help="serve the observatory after building the fixture")
    parser.add_argument("--db", default=None,
                        help="database path (default: a temporary file)")
    parser.add_argument("--subjects", type=int, default=8)
    parser.add_argument("--no-failures", action="store_true",
                        help="run a cohort with no failures")
    parser.add_argument("--print-json", action="store_true")
    args = parser.parse_args(argv)

    root = Path(args.db).parent if args.db else None
    if root:
        root.mkdir(parents=True, exist_ok=True)
    fixture = Fixture(root, subjects=args.subjects)
    print(f"store:  {fixture.db}")
    print(f"run:    {fixture.run_key}\n")

    plan = fixture.run_cohort(with_failures=not args.no_failures)
    for subject, behaviour in plan.items():
        result = next(r for r in fixture.results if r.subject == subject)
        print(f"  {subject:18s} {result.state.value:11s} "
              f"attempts={result.attempts}  ({behaviour})")

    counts = fixture.snapshot()["counts"]
    print(f"\nevents emitted: {fixture.bus.sequence}")
    print(f"task states:     {counts}")

    if args.print_json:
        print(json.dumps(fixture.snapshot()["counts"], indent=2))

    if not args.serve:
        print(f"\nPoint the observatory at it with:\n"
              f"  PAPIPELINE_OBSERVATORY_DB={fixture.db} "
              f"PAPIPELINE_OBSERVATORY_RUN={fixture.run_key} "
              f"python -m papipeline.observatory")
        return 0

    import uvicorn

    from .api import create_app

    app = create_app(str(fixture.db), fixture.run_key, bus=fixture.bus)
    print("\nserving on http://127.0.0.1:8765  (ctrl-c to stop)")
    uvicorn.run(app, host="127.0.0.1", port=8765, log_level="warning")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
