# Triage Labels

The skills speak in terms of five canonical triage roles. This file maps those roles to the actual label strings used in this repo's issue tracker.

Because the tracker is local markdown, a "label" is a value on the `Status:` line at the top of the issue file, written as a bare string (no `#` prefix).

| Label in mattpocock/skills | Label in our tracker | Meaning                                  |
| -------------------------- | -------------------- | ---------------------------------------- |
| `needs-triage`             | `needs-triage`       | Maintainer needs to evaluate this issue  |
| `needs-info`               | `needs-info`         | Waiting on reporter for more information |
| `ready-for-agent`          | `ready-for-agent`    | Fully specified, ready for an AFK agent  |
| `ready-for-human`          | `ready-for-human`    | Requires human implementation            |
| `wontfix`                  | `wontfix`            | Will not be actioned                     |

When a skill mentions a role (e.g. "apply the AFK-ready triage label"), use the corresponding label string from this table.

## Build state, recorded separately from triage

A ticket can be `ready-for-agent` or `wontfix` and still be only half finished.
"Done" is not one state, and collapsing it is how a unit that is finished,
well-tested and unreachable from a run comes to read as a delivered feature.
So a ticket's `Status:` line carries a build state after the triage label:

| Build state | Means |
| --- | --- |
| `built` | Implemented and unit-tested, but **nothing calls it**. Not reachable from a run. |
| `wired` | Implemented, tested, and reachable from the code path a run actually takes. |
| `verified` | Wired, **and** an end-to-end test exercises it on realistic input. |

Examples of the difference that matters: a sample-cap function with thirteen
tests and no caller is `built`; the same cap enforced before the first stage
and proven by a 21-sample run is `verified`. Report the weakest state that is
true — never the strongest.

Edit the right-hand column to match whatever vocabulary you actually use.
