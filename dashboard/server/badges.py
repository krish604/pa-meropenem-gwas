"""§4.2 — the six-badge taxonomy, and the one place a status sentence is made.

The badge is the **single source of truth for the wording** (UI-D2). No page
writes a status sentence from its own logic; it asks this module for the label,
the reason and the tone.

| badge | meaning |
|---|---|
| `completed` | ran, and its declared output is present and readable |
| `running` | announced started, no terminal event |
| `failed` | the process did not succeed |
| `refused` | it declined to run, by name |
| `not_assessed` | ran, and its output cannot answer the question |
| `not_run` | never ran |

Every non-`completed` badge carries a **non-empty reason**, sourced in this
order and no other:

1. `stages_skipped[s]`, verbatim.
2. `REAL_REFUSING_STAGES[s]`, verbatim.
3. `TableRead.reason`.
4. The event that explains it: the `fail` event's timestamp and subject.
5. `"absent from run_manifest.stages and no event explains it"` — the last
   resort, never the first.

`not_assessed` and `not_run` are kept apart on purpose: a stage that ran and
produced nothing answerable is not a stage that never ran, and the pipeline's
own reporting makes the same distinction at `reporting.py:465`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Mapping, Optional, Tuple

from papipeline.run import REAL_REFUSING_STAGES, STAGE_ORDER, UNBUILT_STAGES

#: The six badges. The literal strings the openapi enum lists.
BADGES: Tuple[str, ...] = (
    "completed",
    "running",
    "failed",
    "refused",
    "not_assessed",
    "not_run",
)

#: Human labels. One per badge, so two pages cannot disagree.
LABELS: Mapping[str, str] = {
    "completed": "completed",
    "running": "running",
    "failed": "failed",
    "refused": "refused",
    "not_assessed": "not assessed",
    "not_run": "not run",
}

#: A glyph beside every badge. Status is never carried by hue alone (§7.1).
GLYPHS: Mapping[str, str] = {
    "completed": "●",
    "running": "◐",
    "failed": "✕",
    "refused": "⊘",
    "not_assessed": "○",
    "not_run": "–",
}

#: A tone token the frontend maps to a colour. Never a colour value here —
#: the palette is SHELL's and is read from CSS custom properties (§13).
TONES: Mapping[str, str] = {
    "completed": "ok",
    "running": "busy",
    "failed": "bad",
    "refused": "warn",
    "not_assessed": "unknown",
    "not_run": "unknown",
}

#: The last-resort sentence (step 5 of the reason order).
LAST_RESORT = "absent from run_manifest.stages and no event explains it"

#: The sentence for a manifest written before the run by the provenance writer
#: (`scripts/common/write_provenance.py:101`). Different from "the run did not
#: record this stage", and the UI shows the first (§4).
PRE_RUN_MANIFEST_REASON = (
    "no stage record: this manifest was written by write_provenance.py before "
    "the run, which records no stage outcomes"
)

#: The manifest's own state vocabulary (DESIGN §4.1), as read from
#: `run_manifest.stages`.
STATE_COMPLETED = "completed"
STATE_PREFIX_SKIPPED = "skipped_"


@dataclass(frozen=True)
class StageBadge:
    """One stage's badge, with the reason that justifies it."""

    key: str
    label: str
    glyph: str
    tone: str
    reason: str

    def as_dict(self) -> Dict[str, Any]:
        return {
            "key": self.key,
            "label": self.label,
            "glyph": self.glyph,
            "tone": self.tone,
            "reason": self.reason,
        }


def badge_for(key: str) -> StageBadge:
    """The badge object. The only place a label, glyph or tone is chosen."""
    if key not in BADGES:
        raise KeyError(
            f"{key!r} is not one of the six UI badges ({', '.join(BADGES)}). "
            f"`papipeline/execution/state.py` StageState is the execution "
            f"layer's vocabulary and is deliberately NOT this one (§4.2)."
        )
    return StageBadge(
        key=key,
        label=LABELS[key],
        glyph=GLYPHS[key],
        tone=TONES[key],
        reason="",
    )


def is_skipped_state(state: str) -> bool:
    """Whether a manifest state is one of the `skipped_*` values."""
    return str(state).startswith(STATE_PREFIX_SKIPPED)


def classify(
    stage: str,
    *,
    manifest_state: Optional[str],
    table_present: bool,
    table_reason: str,
    stages_skipped: Mapping[str, str],
    run_mode: Optional[str],
    has_start_event: bool = False,
    has_terminal_event: bool = False,
    has_fail_event: bool = False,
    fail_event_sentence: str = "",
    manifest_records_stages: bool = True,
) -> StageBadge:
    """Apply §4.2 to one stage and return its badge with a non-empty reason.

    The six-row table in DESIGN §4.2, in its own order:

    - `completed`  — `stages[s] == "completed"` **and** the table is present.
    - `running`    — a `start` event with no later `done`/`fail`.
    - `failed`     — an event with `event == "fail"`, or `stages[s] == "failed"`.
    - `refused`    — `s` in `REAL_REFUSING_STAGES` and `run_mode == "REAL"`;
                     or `s` in `UNBUILT_STAGES`; or `stages[s]` starts with
                     `skipped_` and the run recorded a cause.
    - `not_assessed` — `stages[s] == "completed"` but the table is absent or
                     header-only.
    - `not_run`    — absent from `stages`, or in `stages_skipped`.

    `running` and `failed` are checked before `completed` on purpose: an event
    is a fact about the process and the manifest is a fact about the run's
    record of it, and a stage that is mid-run has no terminal record yet.
    """
    skipped_reason = stages_skipped.get(stage)

    # 1. running: announced started, no terminal event.
    if has_start_event and not has_terminal_event:
        return _with_reason(
            "running",
            "a start event with no later done or fail event; the stage is "
            "announced in progress and its outcome is not recorded yet",
        )

    # 2. failed: an event said so, or the manifest says so.
    if has_fail_event:
        return _with_reason(
            "failed",
            fail_event_sentence
            or "a fail event names this stage; the process did not succeed",
        )
    if manifest_state == "failed":
        return _with_reason(
            "failed",
            skipped_reason
            or "run_manifest.stages records this stage as failed",
        )

    # 3. refused: it declined to run, by name.
    refusing = REAL_REFUSING_STAGES.get(stage)
    if refusing and run_mode == "REAL":
        # Step 2 of the reason order, verbatim.
        return _with_reason("refused", refusing)
    if stage in UNBUILT_STAGES:
        return _with_reason(
            "refused",
            skipped_reason
            or (
                f"unbuilt ({UNBUILT_STAGES[stage]}) and not requested"
                if isinstance(UNBUILT_STAGES.get(stage), str)
                else "unbuilt and not requested"
            ),
        )
    if manifest_state and is_skipped_state(manifest_state):
        return _with_reason(
            "refused",
            skipped_reason
            or (
                f"run_manifest.stages records {manifest_state!r} for this "
                f"stage and no cause was recorded"
            ),
        )
    if skipped_reason and manifest_state is None and manifest_records_stages is False:
        return _with_reason("refused", skipped_reason)

    # 4. completed: ran AND its declared output is present and readable.
    if manifest_state == STATE_COMPLETED and table_present:
        return badge_for("completed")

    # 5. not_assessed: ran, and its output cannot answer the question.
    if manifest_state == STATE_COMPLETED and not table_present:
        # Step 3 of the reason order: the TableRead reason.
        return _with_reason("not_assessed", table_reason or LAST_RESORT)

    # 6. not_run: never ran. Distinct from not_assessed on purpose.
    if not manifest_records_stages:
        return _with_reason("not_run", PRE_RUN_MANIFEST_REASON)
    if skipped_reason:
        return _with_reason("not_run", skipped_reason)
    if manifest_state is None:
        return _with_reason("not_run", LAST_RESORT)
    return _with_reason(
        "not_run",
        f"run_manifest.stages records {manifest_state!r} for this stage, which "
        f"is none of completed, failed or skipped_*; there is no table either, "
        f"so nothing about this stage can be reported",
    )


def _with_reason(key: str, reason: str) -> StageBadge:
    """A badge with a non-empty reason.

    The last resort is used when nothing else supplied one, because a
    non-`completed` badge with an empty reason is indistinguishable from a bug.
    """
    base = badge_for(key)
    return StageBadge(
        key=base.key,
        label=base.label,
        glyph=base.glyph,
        tone=base.tone,
        reason=reason or LAST_RESORT,
    )


def badge_map() -> Dict[str, Dict[str, Any]]:
    """The whole vocabulary, for the frontend's `ctx.badges`.

    One object, so two pages cannot disagree about what `not_assessed` means
    and neither can disagree with the report (§7).
    """
    return {
        key: StageBadge(
            key=key, label=LABELS[key], glyph=GLYPHS[key], tone=TONES[key], reason=""
        ).as_dict()
        for key in BADGES
    }


def order() -> Tuple[str, ...]:
    """The display order, which is `STAGE_ORDER`."""
    return tuple(STAGE_ORDER)


__all__ = [
    "BADGES",
    "GLYPHS",
    "LABELS",
    "LAST_RESORT",
    "PRE_RUN_MANIFEST_REASON",
    "TONES",
    "StageBadge",
    "badge_for",
    "badge_map",
    "classify",
    "is_skipped_state",
    "order",
]