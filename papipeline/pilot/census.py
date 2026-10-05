"""Rendering of the pilot census, shared by the report and the deck.

`exclusion_breakdown` (`pilot/cohort.py`) returns counts keyed by a reason-type
prefix -- `corrupt_binary_data`, `undersized`, `zero_byte_file`,
`file_not_found`, `not_fasta_no_header`, plus `other` for anything unmatched.
Those keys are machine identifiers and the report reads better with words, so
they go through the explicit table below.

Two rules, both deliberate:

- **No auto-derivation.** Splitting the snake_case identifier on underscores is
  not a substitute for knowing what the reason means, and it would silently
  produce plausible wrong English for a key added later. An unmapped key is
  rendered as the key itself.
- **Nothing is dropped.** A reason absent from the table still appears, with its
  count. A reason that silently vanishes reads as "none occurred", which is the
  one outcome worse than an ugly sentence.
"""

from __future__ import annotations

from typing import Iterable, Mapping

REASON_PHRASES: Mapping[str, str] = {
    "corrupt_binary_data": "corrupt",
    "zero_byte_file": "zero-byte",
    "undersized": "undersized",
    "oversized": "oversized",
    "file_not_found": "missing",
    "not_fasta_no_header": "not a valid FASTA",
    "other": "other",
}


def reason_phrase(key: str) -> str:
    """Readable label for a reason key; the key itself if unmapped."""
    return REASON_PHRASES.get(key, key)


def describe_reasons(reasons: Mapping[str, int]) -> str:
    """Every reason and count, as prose.

    Includes unmapped keys verbatim rather than dropping or aggregating them, so
    the rendered text always accounts for the full excluded total.
    """
    items: list[tuple[str, int]] = sorted(reasons.items(), key=lambda kv: (-kv[1], kv[0]))
    if not items:
        return "no assemblies were excluded"
    return ", ".join(f"{reason_phrase(key)}: {count:,}" for key, count in items)


def describe_selection(selected: int, discovered: object) -> str:
    """The 'N of M available' phrase, with the total from the run."""
    return f"{selected:,} of {discovered:,} available assemblies"
