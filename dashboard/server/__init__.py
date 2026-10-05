"""Read-only HTTP surface for a *Pseudomonas aeruginosa* imipenem AMR run.

Owned by BACKEND (Phase 1). See `dashboard/DESIGN.md` for the contract this
implements and `dashboard/openapi.yaml` for the wire schema.

Nothing in this package writes to the results root (UI-D1). The only write
location is the state directory, default `~/.pa_dashboard/` (UI-D1, §10).
"""

from __future__ import annotations

__all__ = ["__version__"]

__version__ = "0.1.0-phase1"