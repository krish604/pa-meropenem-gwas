"""PAPipeline Computational Observatory.

A read-only observability surface over the execution state in
:mod:`papipeline.execution`. It renders what the pipeline is actually
doing: which stage, for which genome, in which state, having run which
command, with which result.

It is not a simulator. No value in this package is generated, defaulted or
estimated for display; anything the engine does not record is reported as
unavailable. See :mod:`papipeline.observatory.snapshot` for the projection
and :mod:`papipeline.observatory.metrics` for what can and cannot honestly
be measured on a host.

The package has no effect on how the pipeline runs. The one addition to
the execution layer is an optional ``event_sink`` parameter on
``run_task``, which is ``None`` everywhere it is not explicitly passed.
"""

from __future__ import annotations

__all__ = [
    "EventBus",
    "ExecutionEvent",
    "HostMetrics",
    "Node",
    "Edge",
    "EdgeKind",
    "Snapshot",
    "Throughput",
    "adapt",
    "aggregate_state",
    "bus_sink",
    "build_graph",
    "create_app",
    "event_for_state",
    "graph_payload",
    "read_snapshot",
    "sample_host",
    "task_detail",
]

from .events import (  # noqa: E402
    BUS,
    EventBus,
    ExecutionEvent,
    adapt,
    bus_sink,
    event_for_state,
)
from .graph import (  # noqa: E402
    Edge,
    EdgeKind,
    Node,
    build_graph,
    dependency_edges,
    graph_payload,
    sequence_edges,
)
from .metrics import HostMetrics, psutil_available, sample_host  # noqa: E402
from .snapshot import (  # noqa: E402
    Snapshot,
    aggregate_state,
    read_snapshot,
    task_detail,
)
from .throughput import Throughput, median_duration  # noqa: E402


def create_app(*args, **kwargs):
    """Build the FastAPI app.

    Imported lazily so that reading state from this package does not
    require FastAPI to be installed. A pipeline that only wants to project
    its store into a snapshot should not be forced to pull in a web stack.
    """
    from .api import create_app as _create_app

    return _create_app(*args, **kwargs)
