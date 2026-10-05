"""Real host metrics, read from the operating system.

There is no per-task CPU or memory in this build, because the execution
engine does not record it: the store holds only a static host capacity
snapshot taken at write time, and ``run_task`` does not surface a child
process id. So this module reports **host** metrics, honestly labelled, and
the per-task fields are ``None`` rather than the host's numbers wearing a
task's name.

Every probe is guarded. A metric that cannot be read is reported as
unavailable; a metric that is not supported on this platform is omitted
rather than approximated.
"""

from __future__ import annotations

import os
import platform
import shutil
import socket
import time
from dataclasses import dataclass, field
from typing import Any, Dict, Optional

try:  # psutil is the only way to get real instantaneous CPU on macOS
    import psutil
except ImportError:  # pragma: no cover - exercised on hosts without psutil
    psutil = None  # type: ignore[assignment]


def psutil_available() -> bool:
    return psutil is not None


@dataclass
class HostMetrics:
    """One sample of the machine's real state.

    ``None`` always means "not measurable here", never "zero" and never a
    guess. The UI renders ``None`` as an explicit dash.
    """

    sampled_at: float
    host: str
    platform: str
    python: str
    cpu_count: int
    cpu_percent: Optional[float] = None
    load_average_1m: Optional[float] = None
    memory_total_mb: Optional[int] = None
    memory_used_mb: Optional[int] = None
    memory_percent: Optional[float] = None
    disk_read_bps: Optional[float] = None
    disk_write_bps: Optional[float] = None
    uptime_seconds: Optional[float] = None
    notes: Dict[str, str] = field(default_factory=dict)

    def to_row(self) -> Dict[str, Any]:
        return {
            "sampled_at": self.sampled_at,
            "host": self.host,
            "platform": self.platform,
            "python": self.python,
            "cpu_count": self.cpu_count,
            "cpu_percent": self.cpu_percent,
            "load_average_1m": self.load_average_1m,
            "memory_total_mb": self.memory_total_mb,
            "memory_used_mb": self.memory_used_mb,
            "memory_percent": self.memory_percent,
            "disk_read_bps": self.disk_read_bps,
            "disk_write_bps": self.disk_write_bps,
            "uptime_seconds": self.uptime_seconds,
            "notes": dict(self.notes),
        }


class _DiskRate:
    """Disk throughput needs two samples; a single read is not a rate."""

    def __init__(self) -> None:
        self._last: Optional[tuple] = None
        self._last_at: Optional[float] = None

    def sample(self) -> tuple:
        now = time.monotonic()
        if psutil is None:
            return (None, None)
        io = psutil.disk_io_counters()
        if io is None:
            return (None, None)
        read, write, at = float(io.read_bytes), float(io.write_bytes), now
        if self._last is None or at <= self._last_at:  # type: ignore[operator]
            self._last, self._last_at = (read, write), at
            return (None, None)
        dt = at - self._last_at  # type: ignore[operator]
        rates = ((read - self._last[0]) / dt, (write - self._last[1]) / dt)
        self._last, self._last_at = (read, write), at
        return (max(0.0, rates[0]), max(0.0, rates[1]))


_DISK = _DiskRate()


def sample_host(disk_root: str = "/") -> HostMetrics:
    """Read the real state of this machine once.

    ``cpu_percent`` is a delta since the previous call, so the first call
    legitimately returns ``None`` — there is nothing to compare against
    yet. Returning 0.0 there would be a fabricated reading.
    """
    notes: Dict[str, str] = {}
    metrics = HostMetrics(
        sampled_at=time.time(),
        host=socket.gethostname(),
        platform=platform.platform(),
        python=platform.python_version(),
        cpu_count=os.cpu_count() or 0,
    )

    if psutil is None:
        notes["psutil"] = (
            "psutil is not installed, so CPU, memory and disk I/O cannot be "
            "sampled. Install it; the UI will show these as unavailable "
            "until then."
        )
        return metrics

    try:
        # interval=None compares against the previous call, so this is a
        # real delta rather than a blocking one-second measurement.
        metrics.cpu_percent = float(psutil.cpu_percent(interval=None))
    except Exception:  # noqa: BLE001 - never let a probe break the UI
        notes["cpu_percent"] = "unavailable on this host"
    try:
        metrics.load_average_1m = os.getloadavg()[0]
    except (AttributeError, OSError):
        notes["load_average"] = "not available on this platform"
    try:
        vm = psutil.virtual_memory()
        metrics.memory_total_mb = int(vm.total / (1024 * 1024))
        metrics.memory_used_mb = int((vm.total - vm.available) / (1024 * 1024))
        metrics.memory_percent = float(vm.percent)
    except Exception:  # noqa: BLE001
        notes["memory"] = "unavailable on this host"
    try:
        read_bps, write_bps = _DISK.sample()
        metrics.disk_read_bps = read_bps
        metrics.disk_write_bps = write_bps
        if read_bps is None and write_bps is None:
            notes["disk_io"] = "first sample; a rate needs two readings"
    except Exception:  # noqa: BLE001
        notes["disk_io"] = "unavailable on this host"
    try:
        metrics.uptime_seconds = float(time.time() - psutil.boot_time())
    except Exception:  # noqa: BLE001
        metrics.uptime_seconds = None
    try:
        shutil.disk_usage(disk_root)
    except OSError:
        notes["disk_root"] = f"{disk_root} is not readable"
    return metrics
