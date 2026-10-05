"""Refuse a REAL run whose declared inputs are not all present.

Both `convergence` and `cooccurrence` had the same defect in opposite forms.
Each accepted `mode`, never branched on it, and computed from whatever
intermediates happened to exist - producing a result indistinguishable from
TEST mode's. `46f5005` replaced that with a flat refusal: no REAL caller at
all. That was correct as far as it went, but it left the stages unreachable
even once their inputs are real, which is a different problem with the same
cause.

The fix is not to remove the guard. It is to make the guard **check its inputs**
rather than refuse unconditionally, so the stage runs when - and only when -
every input it declares is actually present.

The distinction that matters, and why this is a module rather than two copies:

* An **empty** mapping and a mapping that is present-but-degenerate are
  different failures with the same symptom. `lineages={"S1": "unknown"}` is not
  an absence of lineage data; it is a *claim* that every sample is in a lineage
  called `unknown`, which would classify every determinant `UNKNOWN` and read
  as "no convergence found". So :func:`is_effectively_empty` treats
  all-sentinel as empty, and each input declares its own sentinel.
* A **refusal must name what is missing.** `require_tool` in
  `adapters/external.py` is the idiom this follows: raise a typed error whose
  message says what was asked for, what is absent, and what would fix it. A
  refusal naming only the stage is not actionable, and this one fires deep in a
  DAG where the reader has no other signal.

Deliberately **not** a warning-and-continue. A stage that logs "some inputs
missing" and writes a partial table produces exactly the artefact `46f5005`
refused to produce: a file that reads as a finding and is not one.
"""

from __future__ import annotations

from typing import Any, Dict, Mapping, Optional, Sequence

from ..errors import StageError


def is_effectively_empty(
    values: Optional[Mapping[str, Any]],
    *,
    sentinel: Optional[str] = None,
) -> bool:
    """Whether a stage input carries no usable content.

    Args:
        values: The input mapping, or ``None`` when it was never supplied.
        sentinel: A value that means "no information" rather than a real
            label. Every value being the sentinel is treated as empty, because
            a mapping full of sentinels is a stage *claiming* to have data and
            having none - which is worse than an absent mapping, since the
            absence at least announces itself.

    Returns:
        True when there is nothing to compute from.
    """
    if not values:
        return True
    if sentinel is None:
        return False
    return all(
        value == sentinel
        for value in values.values()
        if not isinstance(value, (list, tuple, set, frozenset, dict))
    )


def missing_required_inputs(
    inputs: Mapping[str, Mapping[str, Any]],
    *,
    sentinels: Optional[Mapping[str, str]] = None,
) -> Dict[str, str]:
    """Which declared inputs carry no usable content, and why.

    Args:
        inputs: input name -> that input's mapping.
        sentinels: input name -> the sentinel value that means "absent" for
            that input. Inputs not named here are checked for emptiness only.

    Returns:
        input name -> a one-line reason. Empty when everything is present, which
        is the only case in which a REAL run should proceed.
    """
    sentinels = sentinels or {}
    missing: Dict[str, str] = {}
    for name, values in inputs.items():
        sentinel = sentinels.get(name)
        if sentinel is None and is_effectively_empty(values):
            missing[name] = "no rows carry a value"
        elif sentinel is not None and is_effectively_empty(values, sentinel=sentinel):
            missing[name] = (
                f"every row is the {sentinel!r} sentinel, which is an absence "
                f"of the data rather than a value"
            )
    return missing


def refuse_incomplete(
    stage: str,
    missing: Mapping[str, str],
    *,
    paths: Optional[Mapping[str, str]] = None,
    blocked_by: Optional[Mapping[str, str]] = None,
) -> "StageError":
    """The error a REAL stage raises when a declared input is absent.

    Built rather than raised so a caller can attach context and re-raise, and
    so the message is constructed in exactly one place.

    Args:
        stage: The stage refusing, for the message.
        missing: input name -> reason, from :func:`missing_required_inputs`.
        paths: input name -> the path it would be read from. Included so the
            reader can see *where* the absent file was expected.
        blocked_by: input name -> what would produce it (a ticket, a tool, a
            blocked stage). This is what turns "missing" into "missing, and
            here is what is holding it".

    Returns:
        The error to raise. Nothing is written before it is raised, so no
        partial table is left behind.
    """
    paths = paths or {}
    blocked_by = blocked_by or {}
    detail = []
    for name in sorted(missing):
        line = f"  - {name}: {missing[name]}"
        if name in paths:
            line += f" (expected at {paths[name]})"
        if name in blocked_by:
            line += f" [blocked by {blocked_by[name]}]"
        detail.append(line)
    return StageError(
        f"REAL-mode {stage} refuses to run: "
        f"{len(missing)} of its declared inputs carry no data. Computing it "
        f"anyway would write a partial table that reads as a finding. "
        f"Missing: "
        + "\n".join(detail),
        stage=stage,
        missing=",".join(sorted(missing)),
    )