"""§5.4 — provenance: what ran, from what, and what is not recorded.

Tool versions and availability (`tools_detected` or `stage_tool_requirements`,
whichever the manifest writer produced) · the `references` rows with
`version_status` and the `unpinned_references` list · the annotation-reuse
record and the literal `"Bakta executions: N"` · tripwire and watcher log tails
with a size cap · SHA-256 per stage artefact · git HEAD and
`git status --porcelain` of the opened root · the config snapshot · free disk
on the results volume and system load average.

**"Bakta executions: N" is reported as a count or not at all** (assumption A6).
Neither manifest writer carries `bakta_executions`, so when the manifest does
not record it the field is `None` and the UI renders `not reported`. It is
**never `0`**: `0` is a false claim of proof — it says Bakta was never invoked,
which is exactly the claim the reuse provenance exists to make checkable.
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from papipeline.stages.report_tables import NOT_PRODUCED, TableRead

#: The literal the UI shows, and the string the count is formatted into.
BAKTA_EXECUTIONS_TEMPLATE = "Bakta executions: {n}"
NOT_REPORTED = "not reported"

#: No log path is contracted anywhere (A7), so an absent log is `not produced`
#: with its probes named rather than an empty log that reads as "nothing was
#: reported".
LOG_PROBE_NAMES: Tuple[str, ...] = ("tripwire", "watcher")


def bakta_executions_label(count: Optional[int]) -> str:
    """`"Bakta executions: N"`, or `not reported`.

    The `None` case is the honest one: `bakta_executions` is in neither manifest
    writer, so a count that is not recorded must not become a zero.
    """
    if count is None:
        return NOT_REPORTED
    return BAKTA_EXECUTIONS_TEMPLATE.format(n=int(count))


def tools_from_manifest(manifest: Optional[Mapping[str, Any]]) -> Tuple[Dict[str, Any], str]:
    """The tool matrix, and which manifest key it came from.

    The pre-run writer produces `stage_tool_requirements` (stage -> tool ->
    info) and the run writer produces `tools_detected` (tool -> info). The
    field name says which, so a reader knows what shape they are looking at.
    """
    if not isinstance(manifest, Mapping):
        return {}, "none"
    detected = manifest.get("tools_detected")
    if isinstance(detected, Mapping):
        return dict(detected), "tools_detected"
    required = manifest.get("stage_tool_requirements")
    if isinstance(required, Mapping):
        # Flatten stage -> tool -> info into tool -> info, keeping the stages
        # that declared each tool. Both writers' rows have `available` and
        # `version`; the run writer adds `executable`.
        merged: Dict[str, Any] = {}
        for stage, tools in required.items():
            if not isinstance(tools, Mapping):
                continue
            for name, info in tools.items():
                if not isinstance(info, Mapping):
                    continue
                entry = merged.setdefault(
                    name, {"available": None, "version": None, "executable": None, "stages": []}
                )
                entry["stages"] = list(entry.get("stages", [])) + [str(stage)]
                if info.get("available") is not None:
                    entry["available"] = info.get("available")
                if info.get("version") is not None:
                    entry["version"] = info.get("version")
        return merged, "stage_tool_requirements"
    return {}, "none"


def references_from_manifest(manifest: Optional[Mapping[str, Any]]) -> List[Dict[str, Any]]:
    """The `references` rows, one shape under either writer.

    Both writers call `run.build_provenance`, so the provenance page has one
    shape to render.
    """
    if not isinstance(manifest, Mapping):
        return []
    rows = manifest.get("references")
    if not isinstance(rows, Sequence):
        return []
    return [dict(row) for row in rows if isinstance(row, Mapping)]


def unpinned_from_manifest(manifest: Optional[Mapping[str, Any]]) -> Optional[List[str]]:
    """`unpinned_references`, or None when the writer did not produce it.

    An unpinned reference is **reported as a failure, not a silent pass**
    (rule 8). Absent is None rather than `[]`, so "the run recorded none" is not
    mistaken for "the writer never said".
    """
    if not isinstance(manifest, Mapping):
        return None
    value = manifest.get("unpinned_references")
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return [str(v) for v in value]
    return None


def read_json(path: Optional[Path]) -> Tuple[Optional[Any], str]:
    """A JSON sidecar, or the sentence saying why there is not one."""
    if path is None:
        return None, f"{NOT_PRODUCED}: no such path is declared or present"
    path = Path(path)
    if not path.exists():
        return None, f"{NOT_PRODUCED}: no file at {path}"
    try:
        return json.loads(path.read_text(encoding="utf-8")), ""
    except (OSError, ValueError) as exc:
        return None, f"{NOT_PRODUCED}: {path} is unreadable: {exc}"


def bakta_from_reuse(read: Optional[TableRead]) -> Dict[str, Any]:
    """The annotation-reuse record, and the count it makes checkable.

    `annotation.write_reuse_provenance` writes one row per sample the stage
    decided about, reused or not, with `reused`, `reason`, `source_dir`,
    `sha256_feature_tsv`, `sha256_gff3` and `bakta_version`. A run that reused
    output invoked Bakta fewer times than the cohort size, and this record is
    what makes that number checkable.
    """
    if read is None or not read.present:
        return {
            "available": False,
            "reason": read.reason if read is not None else f"{NOT_PRODUCED}: no reuse record was probed for",
            "rows": [],
            "n_rows": None,
            "n_reused": None,
            "n_not_reused": None,
            "bakta_versions": [],
        }
    rows = [dict(r) for r in read.rows]
    reused = sum(1 for r in rows if str(r.get("reused")).strip().lower() in ("true", "1", "yes"))
    versions = sorted({str(r.get("bakta_version")) for r in rows if r.get("bakta_version")})
    sources = sorted({str(r.get("source_dir")) for r in rows if r.get("source_dir")})
    return {
        "available": True,
        "reason": "",
        "rows": rows,
        "n_rows": len(rows),
        "n_reused": reused,
        "n_not_reused": len(rows) - reused,
        "bakta_versions": versions,
        "reused_from": sources[0] if len(sources) == 1 else None,
    }


def bakta_executions(
    manifest: Optional[Mapping[str, Any]],
    reuse: Mapping[str, Any],
) -> Optional[int]:
    """How many times Bakta was invoked — or None, never 0.

    Read from the manifest if some future writer records it. Otherwise derived
    from the reuse record as `n_rows - n_reused`, **only when that record is
    complete enough to mean it**: a record holding one row per sample. When
    neither holds it, None, which renders as `not reported`.
    """
    if isinstance(manifest, Mapping):
        for key in ("bakta_executions", "bakta_execution_count"):
            if key in manifest:
                value = manifest[key]
                if isinstance(value, bool):
                    return None
                if isinstance(value, int):
                    return value
    if reuse.get("available") and reuse.get("n_rows") and reuse.get("n_not_reused") is not None:
        # Every non-reused sample ran the tool once.
        return int(reuse["n_not_reused"])
    return None


def digest(path: Path, *, chunk: int = 1 << 20) -> str:
    """SHA-256, streamed. A large artefact is not read into memory."""
    hasher = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(chunk)
            if not block:
                break
            hasher.update(block)
    return hasher.hexdigest()


def digests(paths: Mapping[str, Path]) -> List[Dict[str, Any]]:
    """One digest per readable artefact, with the size and mtime it was taken at.

    A digest is computed at read time, so it describes the file as it is now
    rather than as it was when a run wrote it. The mtime travels with it.
    """
    out: List[Dict[str, Any]] = []
    for key, path in sorted(paths.items()):
        path = Path(path)
        if not path.exists() or not path.is_file():
            continue
        try:
            stat = path.stat()
            out.append(
                {
                    "key": key,
                    "path": str(path.resolve()),
                    "sha256": digest(path),
                    "size": stat.st_size,
                    "mtime_ns": stat.st_mtime_ns,
                }
            )
        except OSError:
            continue
    return out


def git_facts(repo_root: Optional[Path]) -> Dict[str, Any]:
    """HEAD, branch, dirty flag and `git status --porcelain` of the opened root.

    Read through `subprocess` with an argument list, never a shell string, so a
    path with a space or a metacharacter cannot become a command. A directory
    that is not a git worktree yields nulls rather than a refusal: the
    dashboard reports a run, not a repository's cleanliness.
    """
    empty: Dict[str, Any] = {
        "head": None,
        "branch": None,
        "dirty": None,
        "status_lines": [],
        "reason": None,
    }
    if repo_root is None:
        empty["reason"] = "no repository root was opened"
        return empty
    repo_root = Path(repo_root)
    if not (repo_root / ".git").exists():
        empty["reason"] = f"{repo_root} is not a git worktree"
        return empty
    try:
        head = _git(repo_root, ["rev-parse", "HEAD"])
        branch = _git(repo_root, ["rev-parse", "--abbrev-ref", "HEAD"])
        status = _git(repo_root, ["status", "--porcelain"])
    except (OSError, subprocess.SubprocessError) as exc:
        empty["reason"] = f"git could not be read at {repo_root}: {exc}"
        return empty
    lines = [line for line in status.splitlines() if line.strip()]
    return {
        "head": head or None,
        "branch": branch or None,
        "dirty": bool(lines) if status is not None else None,
        "status_lines": lines[:500],
        "reason": None,
    }


def _git(repo_root: Path, args: Sequence[str]) -> str:
    result = subprocess.run(
        ["git", "-C", str(repo_root), *args],
        capture_output=True,
        text=True,
        timeout=15,
        check=False,
    )
    if result.returncode != 0:
        return ""
    return result.stdout.strip()


def config_snapshot(repo_root: Optional[Path], machine: Optional[str]) -> Dict[str, str]:
    """`config/science.yaml` plus the machine overlay actually loaded.

    Verbatim text, not a parsed re-serialisation: a config snapshot whose
    formatting or comments differ from the file it claims to snapshot is not a
    snapshot.
    """
    if repo_root is None:
        return {}
    root = Path(repo_root)
    out: Dict[str, str] = {}
    for key, path in (
        ("science.yaml", root / "config" / "science.yaml"),
        ("machine_overlay", root / "config" / "machines" / f"{machine}.yaml" if machine else None),
    ):
        if path is None or not path.exists():
            continue
        try:
            out[key] = path.read_text(encoding="utf-8")
        except OSError:
            continue
    return out


def host_facts(path: Optional[Path]) -> Dict[str, Any]:
    """Free disk on the results volume, and the system load average.

    `disk_free_bytes` is None when the volume cannot be stat'd — `0` bytes free
    would be a claim about a filesystem nobody measured.
    """
    facts: Dict[str, Any] = {
        "disk_free_bytes": None,
        "disk_total_bytes": None,
        "load_average": None,
    }
    if path is not None:
        try:
            usage = os.statvfs(str(Path(path)))
            facts["disk_free_bytes"] = int(usage.f_bavail * usage.f_frsize)
            facts["disk_total_bytes"] = int(usage.f_blocks * usage.f_frsize)
        except (OSError, AttributeError, ValueError):
            pass
    try:
        facts["load_average"] = [float(v) for v in os.getloadavg()]
    except (OSError, AttributeError):
        # Not every platform has getloadavg. `None` is "not available here",
        # which is not the same as a load of zero.
        facts["load_average"] = None
    return facts


def tail_text(path: Optional[Path], *, lines: int = 500, max_bytes: int = 1 << 20) -> Dict[str, Any]:
    """The last `lines` lines of a log, with a byte cap.

    Read backwards in chunks so a 10 GB log is streamed rather than buffered
    (DESIGN §11, "a huge file exhausting memory").
    """
    if path is None:
        return {"present": False, "reason": f"{NOT_PRODUCED}: no such log path is declared", "path": None, "n_lines": None, "lines": []}
    path = Path(path)
    if not path.exists():
        return {"present": False, "reason": f"{NOT_PRODUCED}: no file at {path}", "path": str(path), "n_lines": None, "lines": []}
    collected: List[str] = []
    total: Optional[int] = None
    try:
        size = path.stat().st_size
        total = size
        with open(path, "rb") as handle:
            handle.seek(max(0, size - max_bytes))
            block = handle.read(max_bytes)
        text_lines = block.decode("utf-8", errors="replace").splitlines()
        collected = text_lines[-lines:]
    except OSError as exc:
        return {"present": False, "reason": f"{NOT_PRODUCED}: {path} is unreadable: {exc}", "path": str(path), "n_lines": None, "lines": []}
    return {
        "present": True,
        "reason": "",
        "path": str(path),
        "n_lines": len(collected),
        "size_bytes": total,
        "lines": collected,
    }


def provenance(
    manifest: Optional[Mapping[str, Any]],
    reuse: Mapping[str, Any],
    *,
    repo_root: Optional[Path],
    results_root: Optional[Path],
    machine: Optional[str],
    writer: str,
) -> Dict[str, Any]:
    """The whole `/api/provenance` payload.

    Assembled once so the two endpoints cannot disagree about a count.
    """
    tools, tool_source = tools_from_manifest(manifest)
    return {
        "manifest_writer": writer,
        "tools": tools,
        "tool_source": tool_source,
        "references": references_from_manifest(manifest),
        "unpinned_references": unpinned_from_manifest(manifest),
        # None, not 0. See the module docstring.
        "bakta_executions": bakta_executions(manifest, reuse),
        "bakta_executions_label": bakta_executions_label(
            bakta_executions(manifest, reuse)
        ),
        "annotation_reuse": {
            "available": reuse.get("available", False),
            "reason": reuse.get("reason", ""),
            "n_rows": reuse.get("n_rows"),
            "n_reused": reuse.get("n_reused"),
        },
        "git": git_facts(repo_root),
        "config_snapshot": config_snapshot(repo_root, machine),
        "host": host_facts(results_root),
    }


__all__ = [
    "BAKTA_EXECUTIONS_TEMPLATE",
    "LOG_PROBE_NAMES",
    "NOT_REPORTED",
    "bakta_executions",
    "bakta_executions_label",
    "bakta_from_reuse",
    "config_snapshot",
    "digest",
    "digests",
    "git_facts",
    "host_facts",
    "provenance",
    "read_json",
    "references_from_manifest",
    "tail_text",
    "tools_from_manifest",
    "unpinned_from_manifest",
]