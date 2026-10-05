"""`python -m dashboard` — start the read-only results dashboard.

    python -m dashboard --results <path> [--port 8765] [--allow-launch]

Binds `127.0.0.1` by default. A non-loopback bind requires `PA_DASH_TOKEN`, and
`0.0.0.0` is refused outright (UI-D4.3): the dashboard is reached over an SSH
tunnel, and a wildcard bind on an analysis machine is the data leak spec.md D10
names.

The launcher is **disabled unless `--allow-launch`**, and even with the flag the
only action it offers is a dry run. A real run needs the typed phrase
`run real samples`, the `bigmachine` overlay, and `PIPELINE_ALLOW_REAL_MODE` on
the child process — and this build has no endpoint that starts one.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Optional, Sequence

from .server.app import build_state, create_app
from .server.launcher import EXPENSIVE_TOOLS, REAL_PHRASE
from .server.security import (
    BindRefused,
    DEFAULT_HOST,
    DEFAULT_PORT,
    TOKEN_ENV,
)

REPO_ROOT = Path(__file__).resolve().parent.parent


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m dashboard",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--results",
        "--results-root",
        dest="results",
        type=Path,
        default=None,
        help=(
            "The results root to read. Auto-detected when omitted: "
            "PA_DASH_RESULTS_ROOT, then a delivery bundle, then the live REAL, "
            "TEST and STUB roots from the machine overlay's paths.results_root."
        ),
    )
    parser.add_argument(
        "--bundle",
        dest="bundle",
        type=Path,
        default=None,
        help="A delivery bundle root, probed for 02_stage_outputs/ and the rest.",
    )
    parser.add_argument(
        "--state-dir",
        dest="state_dir",
        type=Path,
        default=None,
        help=(
            "Where the index and the dashboard's own log go. Defaults to "
            "$PA_DASH_STATE_DIR, then ~/.pa_dashboard. Never inside the results "
            "root (UI-D1)."
        ),
    )
    parser.add_argument(
        "--host",
        default=DEFAULT_HOST,
        help=f"Bind address. Loopback by default; {TOKEN_ENV} is required for anything else.",
    )
    parser.add_argument("--port", type=int, default=DEFAULT_PORT, help="Bind port.")
    parser.add_argument(
        "--machine",
        default=None,
        help="Machine overlay name for the config (laptop, bigmachine).",
    )
    parser.add_argument(
        "--allow-launch",
        action="store_true",
        help=(
            "Enable the launcher. Even with it, the only action is a dry run "
            "that starts nothing."
        ),
    )
    parser.add_argument(
        "--no-serve",
        action="store_true",
        help="Build the app, report what opened, and exit. Useful for a smoke check.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    # A state directory inside the opened results root would break UI-D1, and
    # the failure is far clearer here than at the first write.
    if args.state_dir is not None and args.results is not None:
        candidate = args.state_dir.expanduser().resolve()
        root = args.results.expanduser().resolve()
        if candidate == root or root in candidate.parents:
            print(
                f"dashboard: refusing --state-dir {candidate}: it is inside the "
                f"results root {root}. UI-D1 makes the state directory the only "
                f"write location, and it must be outside the results root by "
                f"construction.",
                file=sys.stderr,
            )
            return 2

    try:
        state = build_state(
            results_root=args.results,
            bundle=args.bundle,
            repo_root=REPO_ROOT,
            state_dir=args.state_dir,
            host=args.host,
            port=args.port,
            allow_launch=args.allow_launch,
            machine=args.machine,
        )
    except BindRefused as exc:
        print(f"dashboard: {exc}", file=sys.stderr)
        return 2

    app = create_app(state)

    if state.source is None:
        failure = state.failure
        print("dashboard: no results source could be opened.", file=sys.stderr)
        if failure is not None:
            print(failure.render(), file=sys.stderr)
        print(
            "The API is still running and will answer 503 with the same "
            "message, so the failure is inspectable over HTTP.",
            file=sys.stderr,
        )
    else:
        print(
            f"dashboard: opened a {state.kind} results root at "
            f"{state.source.root} (manifest writer: {state.manifest_writer})",
            file=sys.stderr,
        )

    print(
        f"dashboard: state directory {state.state_dir.root} "
        f"(the only write location; UI-D1)",
        file=sys.stderr,
    )
    print(
        "dashboard: launcher "
        + ("ENABLED (dry run only; nothing is started)" if state.allow_launch else "DISABLED (pass --allow-launch)")
        + f". A REAL run would need the phrase {REAL_PHRASE!r} typed exactly, "
        f"the bigmachine overlay, and PIPELINE_ALLOW_REAL_MODE on the child "
        f"process only. Expensive tools to expect: "
        f"{', '.join(EXPENSIVE_TOOLS)}.",
        file=sys.stderr,
    )

    if args.no_serve:
        return 0

    try:
        import uvicorn
    except ImportError:
        print(
            "dashboard: uvicorn is not installed in this environment, so the "
            "server cannot be started. The application object was built "
            "successfully; import it with "
            "`from dashboard.server.app import build_state, create_app` to use "
            "it under another ASGI server.",
            file=sys.stderr,
        )
        return 1

    if not state.bind.loopback and not os.environ.get(TOKEN_ENV):
        # Belt and braces: `resolve_bind` already refused this, so reaching
        # here means the check was bypassed rather than that it passed.
        print(
            f"dashboard: refusing to serve on {state.bind.host} without "
            f"{TOKEN_ENV}.",
            file=sys.stderr,
        )
        return 2

    uvicorn.run(app, host=state.bind.host, port=state.bind.port, log_level="info")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())