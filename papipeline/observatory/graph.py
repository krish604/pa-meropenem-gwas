"""The task network, built from the pipeline's own declared structure.

Nothing here invents a dependency. Nodes are the stages in
:data:`papipeline.run.EXECUTION_ORDER`; edges are the prerequisites in
:data:`papipeline.run.PREREQUISITES`. Where the UI also shows the order
stages execute in, that is a *separate, separately-labelled* edge class
(:data:`EdgeKind.SEQUENCE`) so that a scheduling relationship is never
displayed as though it were a data dependency.

Layout is deterministic. The same structure always produces the same
coordinates, which is what makes the network diffable and testable; a
force-directed layout that re-seeds per request would make the picture
move for no reason and would be untestable.
"""

from __future__ import annotations

import enum
import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Dict, List, Optional, Sequence, Tuple

if TYPE_CHECKING:  # pragma: no cover
    from ..run import EXECUTION_ORDER, PREREQUISITES


def _pipeline_shape() -> tuple:
    """The pipeline's own order and prerequisites, resolved on first use.

    Imported lazily rather than at module scope because ``papipeline.run``
    imports the observatory to record its stages, so a top-level import
    here would be a cycle. Deferring it keeps the observatory a read-only
    observer with no import-time dependency on the orchestrator.
    """
    from ..run import EXECUTION_ORDER, PREREQUISITES

    return EXECUTION_ORDER, PREREQUISITES


def _order(order=None):
    return tuple(order) if order is not None else _pipeline_shape()[0]


def default_order() -> tuple:
    """The pipeline's execution order, as a tuple."""
    return _pipeline_shape()[0]


def _prereqs(prerequisites=None):
    return _pipeline_shape()[1] if prerequisites is None else prerequisites


class EdgeKind(str, enum.Enum):
    """Why two stages are connected.

    The distinction is load-bearing. ``DEPENDENCY`` means stage B reads
    stage A's output and will not run without it. ``SEQUENCE`` means only
    that A is scheduled before B; B does not necessarily read A. The UI
    draws them differently and the legend says so.
    """

    DEPENDENCY = "dependency"
    SEQUENCE = "sequence"


@dataclass(frozen=True)
class Node:
    """One stage in the network."""

    stage: str
    label: str
    depth: int
    index: int
    x: float = 0.0
    y: float = 0.0
    #: Populated from the execution store; never from a default guess.
    total: int = 0
    counts: Dict[str, int] = field(default_factory=dict)
    state: str = "PENDING"
    elapsed_seconds: float = 0.0

    def to_row(self) -> Dict[str, object]:
        return {
            "stage": self.stage,
            "label": self.label,
            "depth": self.depth,
            "index": self.index,
            "x": self.x,
            "y": self.y,
            "total": self.total,
            "counts": dict(self.counts),
            "state": self.state,
            "elapsed_seconds": self.elapsed_seconds,
        }


@dataclass(frozen=True)
class Edge:
    """A directed connection between two stages."""

    source: str
    target: str
    kind: str

    def to_row(self) -> Dict[str, object]:
        return {"source": self.source, "target": self.target, "kind": self.kind}


#: Display names. The stage keys are the pipeline's own; these are only
#: labels, and a stage missing from this table falls back to its key
#: title-cased rather than being dropped.
STAGE_LABELS: Dict[str, str] = {
    "validation": "Genome Validation",
    "annotation": "Bakta",
    "mlst": "MLST",
    "amr": "AMR",
    "regulators": "Regulator",
    "structural_variants": "Structural Variants",
    "mechanisms": "Mechanism",
    "virulence": "Virulence",
    "pangenome": "Pan-genome",
    "phylogeny": "Phylogeny",
    "phenotype": "Phenotype",
    "gwas": "GWAS",
    "convergence": "Convergence",
    "cooccurrence": "Co-occurrence",
    "integration": "Integrated Model",
    "reporting": "Reporting",
}


def label_for(stage: str) -> str:
    return STAGE_LABELS.get(stage, stage.replace("_", " ").title())


def dependency_edges(
    order: Optional[Sequence[str]] = None,
    prerequisites: Optional[Dict[str, Sequence[str]]] = None,
) -> List[Edge]:
    """Real data dependencies, from ``PREREQUISITES`` only.

    An edge whose source is not in ``order`` is dropped rather than
    rendered as a floating node: a dependency on something the pipeline
    does not execute cannot be shown as a real link.
    """
    prereqs = _prereqs(prerequisites)
    order = _order(order)
    known = set(order)
    out: List[Edge] = []
    for target in order:
        for source in sorted(prereqs.get(target, ())):
            if source in known:
                out.append(Edge(source, target, EdgeKind.DEPENDENCY.value))
    return out


def sequence_edges(order: Optional[Sequence[str]] = None) -> List[Edge]:
    """Scheduling order, as a distinct and separately-labelled edge class."""
    order = _order(order)
    return [
        Edge(a, b, EdgeKind.SEQUENCE.value)
        for a, b in zip(order, order[1:])
    ]


def depth_of(stage: str, order: Optional[Sequence[str]] = None) -> int:
    """Longest path from any root, over real dependency edges.

    Using longest-path rather than ``index`` means a stage sits below the
    deepest thing it waits for, so arrows point downward on screen.
    """
    order = _order(order)
    prereqs = {k: sorted(v) for k, v in _prereqs().items()}
    memo: Dict[str, int] = {}

    def walk(node: str, seen: Tuple[str, ...] = ()) -> int:
        if node in memo:
            return memo[node]
        if node in seen:  # a cycle; treat as a root rather than recursing
            return 0
        parents = [p for p in prereqs.get(node, ()) if p in set(order)]
        value = 1 + max((walk(p, seen + (node,)) for p in parents), default=-1)
        memo[node] = value
        return value

    return walk(stage)


def layout(
    order: Optional[Sequence[str]] = None,
    *,
    x_gap: float = 150.0,
    y_gap: float = 132.0,
    sweeps: int = 6,
) -> Dict[str, Tuple[float, float]]:
    """Deterministic two-dimensional coordinates.

    ``x`` comes from :data:`EXECUTION_ORDER`, because the order stages run
    in is real and declared. ``y`` is settled by a barycentric sweep: each
    stage drifts toward the vertical midpoint of the stages it actually
    depends on, collisions are then pushed apart, and the sweep repeats.

    It is deterministic on purpose - no random seeding, no dependence on
    iteration order, no floating-point accumulation across runs. A
    scientific instrument should draw the same picture twice, and a layout
    that re-seeds per request would both look unstable and be untestable.
    """
    order = _order(order)
    position = {stage: i for i, stage in enumerate(order)}
    parents: Dict[str, List[str]] = {s: [] for s in order}
    children: Dict[str, List[str]] = {s: [] for s in order}
    for edge in dependency_edges(order):
        parents[edge.target].append(edge.source)
        children[edge.source].append(edge.target)

    # Seed y by declared order so the sweep starts from something stable.
    y = {stage: float(i) for stage, i in position.items()}

    for _ in range(max(1, sweeps)):
        moved = False
        for stage in order:
            neighbours = parents[stage] + children[stage]
            if not neighbours:
                continue
            target = sum(y[n] for n in neighbours) / len(neighbours)
            if abs(target - y[stage]) > 1e-9:
                y[stage] = 0.5 * y[stage] + 0.5 * target
                moved = True
        # Resolve overlap: walk in order and push each node clear of the
        # previous one. Sorting by y first makes this independent of the
        # order stages happen to be visited in.
        placed: List[Tuple[float, str]] = []
        for stage in sorted(order, key=lambda s: (y[s], position[s])):
            value = y[stage]
            if placed and value - placed[-1][0] < 1.0:
                value = placed[-1][0] + 1.0
            placed.append((value, stage))
            y[stage] = value
        if not moved:
            break

    # Centre on the origin so pan/zoom starts balanced.
    values = list(y.values())
    offset = (max(values) + min(values)) / 2.0
    return {
        stage: (position[stage] * x_gap, (y[stage] - offset) * y_gap)
        for stage in order
    }



def build_graph(
    order: Optional[Sequence[str]] = None,
    prerequisites: Optional[Dict[str, Sequence[str]]] = None,
) -> Tuple[List[Node], List[Edge]]:
    order = _order(order)

    """Nodes and edges for the whole pipeline. State fields are left blank.

    Blanking matters: this function knows the *structure* and nothing about
    any particular run. :mod:`papipeline.observatory.snapshot` fills the
    state in from the execution store, so no value here can be a guess.
    """
    coords = layout(order)
    order = _order(order)
    nodes = [
        Node(
            stage=stage,
            label=label_for(stage),
            depth=depth_of(stage, order),
            index=order.index(stage),
            x=coords[stage][0],
            y=coords[stage][1],
        )
        for stage in order
    ]
    edges = dependency_edges(order, prerequisites) + sequence_edges(order)
    return nodes, edges


def graph_payload(
    order: Optional[Sequence[str]] = None,
    prerequisites: Optional[Dict[str, Sequence[str]]] = None,
) -> Dict[str, object]:
    order = _order(order)

    """The structure, as the client receives it."""
    nodes, edges = build_graph(order, prerequisites)
    return {
        "nodes": [n.to_row() for n in nodes],
        "edges": [e.to_row() for e in edges],
        "stages": list(order),
        "legend": {
            "dependency": "data dependency: the target reads the source's output",
            "sequence": "execution order only: no output is passed between them",
        },
    }
