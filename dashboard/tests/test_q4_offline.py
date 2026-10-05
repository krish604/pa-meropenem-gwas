"""Q4 — the dashboard is offline: no external host is ever contacted.

The scan covers every file under ``dashboard/static/`` including ``vendor/``.
Comments and documentation are excluded from the executable scan (the vendored
d3 banners and ``vendor/LICENSES.md`` carry non-fetching upstream URLs), but the
raw occurrences are asserted to be exactly those benign ones, so a real
network reference cannot hide in a comment.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import List, Tuple

from conftest import REPO_ROOT

STATIC = REPO_ROOT / "dashboard" / "static"

#: XML namespace URIs are identifiers, not fetches (DESIGN §11 / Q4).
ALLOWED_URL = re.compile(r"^https?://(www\.)?w3\.org/")

#: The forbidden third-party/analytics tokens, checked in executable content.
FORBIDDEN_TOKENS = ("//cdn", "unpkg", "jsdelivr", "googleapis", "fonts.google", "preconnect", "analytics")

#: Upstream URLs that are non-fetching: d3's banner comments and the vendored
#: LICENSES.md documentation.
BENIGN_URL_HOSTS = ("d3js.org", "unpkg.com", "opensource.org", "github.com", "creativecommons.org")


def _strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    text = re.sub(r"<!--.*?-->", "", text, flags=re.S)
    text = re.sub(r"(?m)//.*$", "", text)
    return text


def _executable_files() -> List[Path]:
    out = []
    for path in sorted(STATIC.rglob("*")):
        if not path.is_file():
            continue
        if path.suffix.lower() in (".md",):
            continue  # documentation, never executed by the browser
        out.append(path)
    return out


def test_no_external_host_in_executable_static_files():
    offenders: List[str] = []
    for path in _executable_files():
        text = _strip_comments(path.read_text(encoding="utf-8", errors="replace"))
        for match in re.finditer(r"https?://[^\s\"')]+", text):
            if not ALLOWED_URL.match(match.group(0)):
                offenders.append(f"{path.relative_to(REPO_ROOT)}: {match.group(0)}")
        for token in FORBIDDEN_TOKENS:
            if token in text:
                line = text[: text.index(token)].count("\n") + 1
                offenders.append(f"{path.relative_to(REPO_ROOT)}:{line}: {token!r}")
    assert not offenders, "external hosts referenced:\n" + "\n".join(offenders)


def test_no_fetch_or_xhr_to_an_absolute_url():
    offenders: List[str] = []
    for path in _executable_files():
        text = path.read_text(encoding="utf-8", errors="replace")
        for line_number, line in enumerate(text.splitlines(), start=1):
            if "fetch(" not in line and "XMLHttpRequest" not in line:
                continue
            if re.search(r"""(fetch|open)\s*\(\s*['"`]https?://""", line):
                offenders.append(f"{path.relative_to(REPO_ROOT)}:{line_number}: {line.strip()}")
    assert not offenders, "absolute fetch/XHR:\n" + "\n".join(offenders)


def test_raw_external_urls_are_only_benign_upstream_references():
    """Every raw `http(s)://` in the tree is a namespace or a non-fetching doc URL."""
    offenders: List[str] = []
    for path in sorted(STATIC.rglob("*")):
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        for match in re.finditer(r"https?://[^\s\"')>]+", text):
            url = match.group(0)
            if ALLOWED_URL.match(url):
                continue
            if any(host in url for host in BENIGN_URL_HOSTS):
                continue
            offenders.append(f"{path.relative_to(REPO_ROOT)}: {url}")
    assert not offenders, "unexpected external URL:\n" + "\n".join(offenders)
