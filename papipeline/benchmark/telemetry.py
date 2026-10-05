"""Process and system telemetry for the benchmark.

Sampling only, on an interval, in a thread per observed process. There is no
permanent monitoring loop: the sampler for a process starts when the process
starts and stops when it exits, and a sampler that cannot measure reports
``None`` rather than a plausible number.

``None`` is load-bearing throughout. A missing measurement written as ``0``
or as ``0.0`` is a fabricated measurement, and the benchmark's whole value
depends on being able to tell "we did not measure this" from "we measured
zero".
"""

from __future__ import annotations

import threading
import time
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional

try:
    import psutil
except ImportError:  # pragma: no cover
    psutil = None  # type: ignore[assignment]


def psutil_available() -> bool:
    return psutil is not None


@dataclass
class ProcessSample:
    """What was measured about one child process."""

    pid: Optional[int] = None
    started_at: Optional[float] = None
    ended_at: Optional[float] = None
    wall_seconds: Optional[float] = None
    cpu_seconds: Optional[float] = None
    peak_rss_mb: Optional[float] = None
    average_rss_mb: Optional[float] = None
    cpu_percent_mean: Optional[float] = None
    samples: int = 0
    threads: Optional[int] = None
    note: str = ""

    def finish(self) -> "ProcessSample":
        """Close the sample. Always yields a wall time if we have a start."""
        if self.ended_at is None:
            self.ended_at = time.time()
        if self.started_at is None:
            self.started_at = self.ended_at
        self.wall_seconds = self.ended_at - self.started_at
        return self

    def to_row(self) -> Dict[str, Any]:
        return {
            "pid": self.pid, "started_at": self.started_at,
            "ended_at": self.ended_at, "wall_seconds": self.wall_seconds,
            "cpu_seconds": self.cpu_seconds, "peak_rss_mb": self.peak_rss_mb,
            "average_rss_mb": self.average_rss_mb,
            "cpu_percent_mean": self.cpu_percent_mean,
            "samples": self.samples, "bakta_threads": self.threads,
            "note": self.note,
        }


class _ProcessSampler:
    """Samples one pid until it exits."""

    def __init__(self, pid: int, interval: float = 1.0) -> None:
        self.pid = pid
        self.interval = interval
        self.result = ProcessSample(pid=pid, started_at=time.time())
        self._stop = threading.Event()
        self._thread: Optional[threading.Thread] = None
        self._rss: List[float] = []
        self._cpu: List[float] = []
        self._lock = threading.Lock()

    def start(self) -> "ProcessSample":
        if psutil is None:
            self.result.note = "psutil is not installed; no process metrics"
            self.result.finish()
            return self.result
        try:
            handle = psutil.Process(self.pid)
        except Exception as exc:  # noqa: BLE001 - the process may already be gone
            self.result.note = f"process not observable: {exc}"
            self.result.finish()
            return self.result
        # Prime the CPU counter; the first call would otherwise return 0.0,
        # which is a measurement artefact rather than an idle process.
        try:
            handle.cpu_percent(None)
        except Exception:  # noqa: BLE001
            pass
        self._thread = threading.Thread(target=self._loop, args=(handle,),
                                        daemon=True)
        self._thread.start()
        return self.result

    def _loop(self, handle) -> None:
        while not self._stop.is_set():
            try:
                with handle.oneshot():
                    rss = handle.memory_info().rss / 2 ** 20
                    cpu = handle.cpu_percent(None)
                    times = handle.cpu_times()
                with self._lock:
                    self._rss.append(rss)
                    self._cpu.append(cpu)
                    # Read CPU time while the process is still alive. It is
                    # unavailable once the process has exited, so collecting
                    # it only at the end would always report None.
                    self.result.cpu_seconds = float(times.user + times.system)
                    self.result.samples += 1
            except Exception:  # noqa: BLE001 - the process exited; that is normal
                break
            self._stop.wait(self.interval)
        self._stop.set()
        self._collect_final(handle)
        self.result.finish()

    def _collect_final(self, handle) -> None:
        try:
            times = handle.cpu_times()
            self.result.cpu_seconds = float(times.user + times.system)
        except Exception:  # noqa: BLE001
            pass
        try:
            memory = handle.memory_info().rss / 2 ** 20
            with self._lock:
                self._rss.append(memory)
                self.result.samples += 1
        except Exception:  # noqa: BLE001
            pass
        with self._lock:
            rss = list(self._rss)
            cpu = list(self._cpu)
        if rss:
            self.result.peak_rss_mb = max(rss)
            self.result.average_rss_mb = sum(rss) / len(rss)
        if cpu:
            self.result.cpu_percent_mean = sum(cpu) / len(cpu)

    def stop(self) -> ProcessSample:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=10)
        return self.result


def sample_process(
    pid: int,
    *,
    interval: float = 1.0,
    threads: Optional[int] = None,
    sink: Optional[Callable[..., Any]] = None,
    run_key: str = "",
    stage: str = "",
    subject: str = "",
) -> ProcessSample:
    """Start sampling ``pid`` and return the accumulating sample.

    A small proxy object is returned so the caller can ``finish()`` it. The
    sampler thread stops on its own once the process exits, so a forgotten
    ``finish()`` cannot leak a thread indefinitely.
    """
    sampler = _ProcessSampler(pid, interval=interval)
    sampler.result.threads = threads
    sampler.start()
    return _ProcessHandle(sampler)


class _ProcessHandle:
    """What callers hold: a sample plus a way to finish it."""

    def __init__(self, sampler: _ProcessSampler) -> None:
        self._sampler = sampler

    def __getattr__(self, name: str) -> Any:
        return getattr(self._sampler.result, name)

    def finish(self) -> ProcessSample:
        return self._sampler.stop()


@dataclass
class SampleSystem:
    """System-wide measurements over one benchmark configuration."""

    started_at: float = 0.0
    ended_at: float = 0.0
    wall_seconds: float = 0.0
    peak_ram_gb: Optional[float] = None
    average_ram_gb: Optional[float] = None
    peak_ram_percent: Optional[float] = None
    cpu_percent_mean: Optional[float] = None
    cpu_percent_peak: Optional[float] = None
    load_average_1m: Optional[float] = None
    load_average_peak: Optional[float] = None
    disk_read_bps: Optional[float] = None
    disk_write_bps: Optional[float] = None
    samples: int = 0
    note: str = ""

    def to_row(self) -> Dict[str, Any]:
        return {
            "wall_seconds": round(self.wall_seconds, 2),
            "peak_ram_gb": self.peak_ram_gb,
            "average_ram_gb": self.average_ram_gb,
            "peak_ram_percent": self.peak_ram_percent,
            "cpu_percent_mean": self.cpu_percent_mean,
            "cpu_percent_peak": self.cpu_percent_peak,
            "load_average_1m": self.load_average_1m,
            "load_average_peak": self.load_average_peak,
            "disk_read_bps": self.disk_read_bps,
            "disk_write_bps": self.disk_write_bps,
            "samples": self.samples,
            "note": self.note,
        }


class SystemSampler:
    """Samples the machine for the duration of one configuration."""

    def __init__(self, interval: float = 1.0) -> None:
        self.interval = interval
        self.result = SampleSystem()
        self._stop = threading.Event()
        self._thread: Optional[threading.Thread] = None
        self._ram: List[float] = []
        self._cpu: List[float] = []
        self._load: List[float] = []

    def start(self) -> "SystemSampler":
        self.result.started_at = time.time()
        if psutil is None:
            self.result.note = "psutil is not installed; no system metrics"
            return self
        psutil.cpu_percent(None)
        self._thread = threading.Thread(target=self._loop, daemon=True)
        self._thread.start()
        return self

    def _loop(self) -> None:
        while not self._stop.is_set():
            try:
                ram = psutil.virtual_memory()
                self._ram.append(ram.used / 2 ** 30)
                self._cpu.append(psutil.cpu_percent(None))
                self._load.append(psutil.getloadavg()[0])
                self.result.peak_ram_percent = max(
                    self.result.peak_ram_percent or 0.0, ram.percent)
            except Exception:  # noqa: BLE001
                break
            self._stop.wait(self.interval)

    def stop(self) -> SampleSystem:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=10)
        self.result.ended_at = time.time()
        self.result.wall_seconds = self.result.ended_at - self.result.started_at
        self.result.samples = len(self._ram)
        if self._ram:
            self.result.peak_ram_gb = max(self._ram)
            self.result.average_ram_gb = sum(self._ram) / len(self._ram)
        if self._cpu:
            self.result.cpu_percent_mean = sum(self._cpu) / len(self._cpu)
            self.result.cpu_percent_peak = max(self._cpu)
        if self._load:
            self.result.load_average_1m = self._load[-1]
            self.result.load_average_peak = max(self._load)
        if psutil is not None:
            try:
                io = psutil.disk_io_counters()
                if io is not None:
                    self.result.disk_read_bps = float(io.read_bytes)
                    self.result.disk_write_bps = float(io.write_bytes)
            except Exception:  # noqa: BLE001
                pass
        return self.result
