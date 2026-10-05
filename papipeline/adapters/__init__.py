"""External tool detection and adapter base classes."""

from __future__ import annotations

from .external import (
    UNPROBED_VERSION,
    VERSION_FLAGS,
    CommandResult,
    ToolAdapter,
    ToolStatus,
    detect_all,
    detect_tools,
    probe_on_demand,
    probe_version,
    require_optional_tool,
    require_tool,
)

__all__ = [
    "UNPROBED_VERSION",
    "VERSION_FLAGS",
    "CommandResult",
    "ToolAdapter",
    "ToolStatus",
    "detect_all",
    "detect_tools",
    "probe_on_demand",
    "probe_version",
    "require_optional_tool",
    "require_tool",
]
