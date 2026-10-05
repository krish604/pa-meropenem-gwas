"""UI-D4 — containment for served files, and the bind/token rule.

Three independent checks, all required (DESIGN §9 UI-D4, §11):

1. **Resolve, then contain.** `candidate = (root / requested).resolve()`,
   which follows every symlink, and then
   `candidate == root or root in candidate.parents`. A traversal that escapes
   fails here even if it never wrote `..` — a symlink inside the root pointing
   at `/etc` resolves to `/etc/passwd` and is not under the root.
2. **Allowlist, then suffix.** Only the declared suffixes are served.
   `text/html` is served as an attachment with `nosniff`.
3. **Bind check.** A non-loopback host requires `PA_DASH_TOKEN`, and
   `0.0.0.0` is always refused. The token is compared with
   `secrets.compare_digest` and is only ever accepted as a header: a query
   string lands in shell history, in the address bar and in every proxy log.
"""

from __future__ import annotations

import os
import secrets
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence, Tuple

#: The token header. Never a query parameter — see the module docstring.
TOKEN_HEADER = "X-PA-Dash-Token"
TOKEN_ENV = "PA_DASH_TOKEN"

#: The default bind. Loopback only, reached over an SSH tunnel (spec.md D10).
DEFAULT_HOST = "127.0.0.1"

#: The documented default port (openapi.yaml `servers`).
DEFAULT_PORT = 8765

#: Loopback names. `::1` is loopback too and needs no token.
LOOPBACK_HOSTS = frozenset({"127.0.0.1", "::1", "localhost", "127.0.0.0/8"})

#: Suffixes that may be served (DESIGN §9 UI-D4.2). `md`, `html`, `json`,
#: `tsv`, `txt`, `csv`, `nwk`, `newick`, `tre`, `svg`, `png`, `pdf`, `log`,
#: `fa`, `fasta`, `faa`, `gff`, `gff3`, `vcf`.
ALLOWED_SUFFIXES = frozenset({
    "md", "html", "htm", "json", "tsv", "txt", "csv", "nwk", "newick", "tre",
    "svg", "png", "pdf", "log", "fa", "fasta", "faa", "gff", "gff3", "vcf",
})

#: Suffixes served as a download rather than inline, so a report's own markup
#: cannot execute in the dashboard's origin.
ATTACHMENT_SUFFIXES = frozenset({"html", "htm"})

#: MIME types by suffix. Only what the allowlist needs.
CONTENT_TYPES = {
    "md": "text/markdown; charset=utf-8",
    "html": "text/html; charset=utf-8",
    "htm": "text/html; charset=utf-8",
    "json": "application/json",
    "tsv": "text/tab-separated-values; charset=utf-8",
    "txt": "text/plain; charset=utf-8",
    "csv": "text/csv; charset=utf-8",
    "nwk": "text/plain; charset=utf-8",
    "newick": "text/plain; charset=utf-8",
    "tre": "text/plain; charset=utf-8",
    "svg": "image/svg+xml",
    "png": "image/png",
    "pdf": "application/pdf",
    "log": "text/plain; charset=utf-8",
    "fa": "text/plain; charset=utf-8",
    "fasta": "text/plain; charset=utf-8",
    "faa": "text/plain; charset=utf-8",
    "gff": "text/plain; charset=utf-8",
    "gff3": "text/plain; charset=utf-8",
    "vcf": "text/plain; charset=utf-8",
}

#: The one hard refusal that ignores the token entirely.
FORBIDDEN_HOSTS = frozenset({"0.0.0.0", "::", ""})


class ContainmentError(Exception):
    """The requested path is not servable. The message is the user-facing one."""

    def __init__(self, message: str, status_code: int = 400) -> None:
        super().__init__(message)
        self.message = message
        self.status_code = status_code


class BindRefused(Exception):
    """The server refuses to start on this host."""


@dataclass(frozen=True)
class BindPolicy:
    """The outcome of the bind check, decided once at startup."""

    host: str
    port: int
    loopback: bool
    token_required: bool

    @property
    def address(self) -> str:
        return f"{self.host}:{self.port}"


def is_loopback(host: str) -> bool:
    """Whether ``host`` names this machine only.

    `localhost` is accepted because it resolves to loopback, and a `127.0.0.0/8`
    address is loopback by definition. Anything else is a non-local bind and
    needs a token.
    """
    if not host:
        return False
    candidate = host.strip().strip("[]").lower()
    if candidate in LOOPBACK_HOSTS:
        return True
    if candidate.startswith("127."):
        return True
    return candidate in {"::1", "0:0:0:0:0:0:0:1"}


def resolve_bind(host: Optional[str], port: Optional[int]) -> BindPolicy:
    """Decide the bind, or refuse.

    `0.0.0.0` is refused even with a token: a wildcard bind on an analysis
    machine is the data leak `spec.md` D10 names, and a shared secret does not
    make it one.
    """
    chosen_host = (host or DEFAULT_HOST).strip()
    chosen_port = int(port) if port else DEFAULT_PORT

    if chosen_host in FORBIDDEN_HOSTS:
        raise BindRefused(
            f"refusing to bind {chosen_host!r}. A wildcard bind exposes every "
            f"result on this machine to the network, which is the data leak "
            f"spec.md D10 names. Bind {_DEFAULT_EXAMPLE} and reach the "
            f"dashboard over an SSH tunnel."
        )

    loopback = is_loopback(chosen_host)
    if loopback:
        return BindPolicy(
            host=chosen_host,
            port=chosen_port,
            loopback=True,
            token_required=bool(os.environ.get(TOKEN_ENV)),
        )

    if not os.environ.get(TOKEN_ENV):
        raise BindRefused(
            f"refusing to bind the non-loopback host {chosen_host!r} without "
            f"{TOKEN_ENV} set. Either bind {_DEFAULT_EXAMPLE} — the normal "
            f"case, reached over an SSH tunnel — or set {TOKEN_ENV} and send it "
            f"as {TOKEN_HEADER} on every /api request. The token is never "
            f"accepted as a query parameter: a query string lands in shell "
            f"history, in the address bar and in every proxy log."
        )
    return BindPolicy(
        host=chosen_host,
        port=chosen_port,
        loopback=False,
        token_required=True,
    )


_DEFAULT_EXAMPLE = "127.0.0.1"


def configured_token() -> Optional[str]:
    """The token, or None. Read from the environment, never from the client."""
    value = os.environ.get(TOKEN_ENV)
    return value if value else None


def check_token(supplied: Optional[str], *, required: bool) -> None:
    """Authorise a request, or raise `ContainmentError`.

    Compared with `secrets.compare_digest` so the comparison time does not
    depend on how many leading characters matched.
    """
    if not required:
        return
    expected = configured_token()
    if not expected:
        # The token disappeared after startup. Fail closed: an unauthenticated
        # non-local bind is the thing UI-D4.3 exists to prevent.
        raise ContainmentError(
            f"this server is bound to a non-loopback host and {TOKEN_ENV} is no "
            f"longer set, so no request can be authorised. Restart the "
            f"dashboard.",
            status_code=503,
        )
    if not supplied or not secrets.compare_digest(str(supplied), str(expected)):
        raise ContainmentError(
            f"missing or wrong {TOKEN_HEADER}. It is required because this "
            f"server is bound to a non-loopback host.",
            status_code=401,
        )


def resolve_contained(root: Path, requested: str) -> Path:
    """Resolve `requested` under `root`, or refuse.

    Args:
        root: The opened results root, already `resolve()`d.
        requested: A path relative to the root, as the client sent it. The
            handler decodes nothing itself; the framework has already decoded
            it, and this re-checks after `resolve()`.

    Raises:
        ContainmentError: The path is absolute, names `..`, escapes the root
            through a symlink, or does not exist.
    """
    root = Path(root).resolve()
    text = str(requested or "")
    if not text.strip():
        raise ContainmentError("no path was given.", status_code=400)

    # An absolute request is refused rather than reinterpreted. It is either a
    # bug in the client or an attempt, and neither has a legitimate reading
    # under a root-relative API.
    candidate = Path(text)
    if candidate.is_absolute() or text.startswith("~"):
        raise ContainmentError(
            f"{text!r} is absolute. This endpoint serves paths relative to the "
            f"opened results root.",
            status_code=400,
        )
    if any(part == ".." for part in candidate.parts):
        raise ContainmentError(
            f"{text!r} contains '..'. The path is refused rather than normalised, "
            f"so a traversal is visible in the log instead of silently resolved.",
            status_code=400,
        )

    resolved = (root / candidate).resolve()
    if resolved != root and root not in resolved.parents:
        raise ContainmentError(
            f"{text!r} resolves to {resolved}, which is outside the opened "
            f"results root {root}. A symlink inside the root that points "
            f"elsewhere fails here for the same reason a '..' does.",
            status_code=400,
        )
    if not resolved.exists():
        raise ContainmentError(
            f"no file at {resolved} under the opened results root.",
            status_code=404,
        )
    if not resolved.is_file():
        raise ContainmentError(
            f"{resolved} is not a regular file.",
            status_code=404,
        )
    return resolved


def check_suffix(path: Path, *, allowed: Iterable[str] = ALLOWED_SUFFIXES) -> str:
    """The suffix is on the allowlist, or refuse. Returns the suffix."""
    suffix = Path(path).suffix.lower().lstrip(".")
    if suffix not in set(allowed):
        raise ContainmentError(
            f"{Path(path).name} has suffix {('.' + suffix) if suffix else '(none)'}, "
            f"which is not on the serving allowlist. Only "
            f"{', '.join(sorted(set(allowed)))} are served.",
            status_code=403,
        )
    return suffix


def content_headers(suffix: str, filename: str) -> Sequence[Tuple[str, str]]:
    """Response headers for a served file, attachment or inline."""
    headers: list[Tuple[str, str]] = [
        # The report's own markup cannot execute in the dashboard's origin, and
        # the browser will not sniff a different type out of a text/html body.
        ("X-Content-Type-Options", "nosniff"),
        ("Cache-Control", "no-store"),
    ]
    if suffix in ATTACHMENT_SUFFIXES:
        headers.append(("Content-Disposition", f'attachment; filename="{filename}"'))
    return headers


__all__ = [
    "ALLOWED_SUFFIXES",
    "ATTACHMENT_SUFFIXES",
    "BindPolicy",
    "BindRefused",
    "CONTAINMENT_ALLOWED",
    "ContainmentError",
    "CONTENT_TYPES",
    "DEFAULT_HOST",
    "DEFAULT_PORT",
    "FORBIDDEN_HOSTS",
    "LOOPBACK_HOSTS",
    "TOKEN_ENV",
    "TOKEN_HEADER",
    "check_suffix",
    "check_token",
    "configured_token",
    "content_headers",
    "is_loopback",
    "resolve_bind",
    "resolve_contained",
]

#: Kept as a name because the openapi response for a rejected path is 403 and a
#: reader of the code should see that 403 is the suffix decision, not 400.
CONTAINMENT_ALLOWED = ALLOWED_SUFFIXES