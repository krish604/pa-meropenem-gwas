"""The smoke-run report marker: asserted, never optional.

A report produced from the bounded smoke overlay must say what it is. A report
that omitted this looks identical to the real-cohort analysis, and the operator
reading it later has no way to tell 10 assemblies from 967. That is the whole
reason the marker exists, so an absent or wrong marker is raised rather than
tolerated: generating a clean-looking report is the failure mode, not the fix.
"""

from __future__ import annotations

from ..errors import PipelineError

#: Prefixed with the count and suffixed with the caveat, deliberately verbose:
#: it should be impossible to miss when skimming a report.
SMOKE_MARKER_PREFIX = "SMOKE TEST — REAL data, N="
SMOKE_MARKER_SUFFIX = " assemblies, not the full cohort analysis."


def smoke_marker(assembly_count: int) -> str:
    """The exact marker line for a run over ``assembly_count`` assemblies."""
    return f"{SMOKE_MARKER_PREFIX}{assembly_count}{SMOKE_MARKER_SUFFIX}"


def require_smoke_marker(text: str, assembly_count: int) -> str:
    """Return ``text`` if it carries the marker for ``assembly_count``.

    Args:
        text: A generated report.
        assembly_count: The number of assemblies the run actually used.

    Returns:
        The report unchanged, once verified.

    Raises:
        PipelineError: The marker is absent, or names a different count. A
            marker for another count is treated as absent, because a report
            claiming 9 assemblies when 10 were used is worse than one with no
            marker at all - it is confidently wrong.
    """
    marker = smoke_marker(assembly_count)
    if marker not in text:
        raise PipelineError(
            "Smoke-run report is missing its required marker, or states a "
            f"different assembly count. Expected exactly: {marker}",
            expected_count=assembly_count,
        )
    return text
