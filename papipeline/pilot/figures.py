"""Figure generation for the pilot-100 analysis.

All figures are 300 DPI, per the pilot specification. Every figure is built
from the analysis tables, so the underlying data can be inspected without
plotting anything.

Two rules are enforced rather than assumed:

* **No figure is drawn without sufficient data.** :func:`build_all` skips a
  figure and records the reason when the underlying table is empty or has too
  few populated categories to be informative. A blank or near-blank plot is
  worse than no plot.
* **No figure asserts a relationship.** Counts and frequencies are shown.
  Co-occurrence figures carry an "association only" annotation, and no figure
  labels a gene as causing a phenotype.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

from ..logging_utils import get_logger

LOGGER = get_logger("pilot.figures")

DPI = 300

#: A figure is skipped below these thresholds, and the reason is reported.
MIN_CATEGORIES = 2
MIN_ROWS = 3
MIN_CELLS = 10

#: Above this many rows a heatmap is unreadable at figure size; the most
#: frequent rows are kept and the omission is noted on the figure.
MAX_HEATMAP_ROWS = 60

#: Annotation placed on co-occurrence figures.
ASSOCIATION_NOTE = "Co-occurrence is association only; it does not imply interaction or causation."


@dataclass
class FigureResult:
    """Outcome of one figure attempt."""

    name: str
    path: Optional[Path]
    created: bool
    reason: Optional[str] = None

    def to_row(self) -> Dict[str, Any]:
        return {
            "figure": self.name,
            "created": "TRUE" if self.created else "FALSE",
            "path": str(self.path) if self.path else ".",
            "reason": self.reason or ".",
        }


def _pyplot():
    """Import pyplot with a non-interactive backend."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _save(fig, path: Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    return path


# --------------------------------------------------------------------------
# bar figures
# --------------------------------------------------------------------------


def gene_frequency_figure(
    rows: Sequence[Mapping[str, Any]], title: str, path: Path, top_n: int = 25
) -> FigureResult:
    """Horizontal bar chart of the most frequent genes."""
    name = path.stem
    if len(rows) < MIN_ROWS:
        return FigureResult(name, None, False, f"only {len(rows)} distinct gene(s)")

    ordered = sorted(rows, key=lambda r: -float(r.get("n_isolates_detected") or 0))
    shown = ordered[:top_n]
    omitted = len(ordered) - len(shown)
    if not shown:
        return FigureResult(name, None, False, "no gene had a non-zero count")

    plt = _pyplot()
    labels = [str(r["Gene"]) for r in shown][::-1]
    counts = [float(r.get("n_isolates_detected") or 0) for r in shown][::-1]
    categories = sorted({str(r.get("Mechanism")) for r in shown})

    height = max(4.0, 0.26 * len(shown))
    fig, ax = plt.subplots(figsize=(10, height))
    colors = plt.get_cmap("tab20")(  # type: ignore[attr-defined]
        [categories.index(str(r.get("Mechanism"))) % 20 for r in shown][::-1]
    )
    ax.barh(labels, counts, color=colors)
    ax.set_xlabel("Isolates with gene detected (n)")
    ax.set_title(title)
    ax.grid(axis="x", linestyle=":", alpha=0.4)
    for index, value in enumerate(counts):
        ax.text(value, index, f" {value:g}", va="center", fontsize=7)
    if omitted > 0:
        ax.text(
            0.99,
            0.01,
            f"{omitted} lower-frequency gene(s) not shown",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=7,
            style="italic",
        )
    fig.tight_layout()
    return FigureResult(name, _save(fig, path), True)


def mechanism_frequency_figure(
    rows: Sequence[Mapping[str, Any]], title: str, path: Path
) -> FigureResult:
    """Grouped bar chart of mechanism counts by phenotype category."""
    name = path.stem
    if len(rows) < MIN_CATEGORIES:
        return FigureResult(
            name, None, False, f"only {len(rows)} mechanism categor(y/ies) present"
        )

    phenotype_keys: List[str] = []
    for row in rows:
        for key in row:
            if key.startswith("n_") and key not in {
                "n_isolates",
                "n_isolates_with_phenotype",
                "n_genes_in_category",
            }:
                label = key[2:]
                if label not in phenotype_keys:
                    phenotype_keys.append(label)
    if not phenotype_keys:
        return FigureResult(name, None, False, "no phenotype split available")

    mechanisms = [str(r["Mechanism"]) for r in rows]
    plt = _pyplot()
    fig, ax = plt.subplots(figsize=(max(8, 1.6 * len(mechanisms)), 5.5))
    width = 0.8 / len(phenotype_keys)
    positions = range(len(mechanisms))
    for offset, key in enumerate(phenotype_keys):
        values = [float(r.get(f"n_{key}") or 0) for r in rows]
        ax.bar(
            [p + offset * width for p in positions],
            values,
            width=width,
            label=key,
        )
    ax.set_xticks([p + 0.4 - width / 2 for p in positions])
    ax.set_xticklabels(mechanisms, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Isolates (n)")
    ax.set_title(title)
    ax.legend(title="Phenotype", fontsize=8)
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    fig.tight_layout()
    return FigureResult(name, _save(fig, path), True)


# --------------------------------------------------------------------------
# heatmaps
# --------------------------------------------------------------------------


def _binary_matrix(
    members: Sequence[Any],
    feature_of: Callable[[Any], Sequence[str]],
    feature_names: Optional[Sequence[str]] = None,
) -> Tuple[List[str], List[str], List[List[int]]]:
    names = list(feature_names) if feature_names else []
    if not names:
        collected = {f for m in members for f in feature_of(m)}
        names = sorted(collected)
    order = sorted(
        names,
        key=lambda n: -sum(1 for m in members if n in feature_of(m)),
    )[:MAX_HEATMAP_ROWS]
    labels = [m.pilot_id for m in members]
    matrix = [
        [1 if feature in feature_of(m) else 0 for m in members] for feature in order
    ]
    return order, labels, matrix


def heatmap_figure(
    members: Sequence[Any],
    feature_of: Callable[[Any], Sequence[str]],
    title: str,
    path: Path,
    row_label: str = "Feature",
) -> FigureResult:
    """Binary presence/absence heatmap, features as rows and isolates as columns."""
    name = path.stem
    features, labels, matrix = _binary_matrix(members, feature_of)
    if len(features) < MIN_ROWS:
        return FigureResult(
            name, None, False, f"only {len(features)} distinct {row_label}(s)"
        )
    populated = sum(sum(row) for row in matrix)
    if populated < MIN_CELLS:
        return FigureResult(
            name, None, False, f"only {populated} populated cell(s); too sparse to plot"
        )

    plt = _pyplot()
    width = max(8.0, 0.12 * len(labels))
    height = max(4.0, 0.16 * len(features))
    fig, ax = plt.subplots(figsize=(width, height))
    ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap="Blues", vmin=0, vmax=1)
    ax.set_yticks(range(len(features)))
    ax.set_yticklabels(features, fontsize=6)
    step = max(1, len(labels) // 40)
    ax.set_xticks(range(0, len(labels), step))
    ax.set_xticklabels(labels[::step], rotation=90, fontsize=5)
    ax.set_xlabel("Isolate (Pilot_ID)")
    ax.set_ylabel(row_label)
    ax.set_title(title)
    truncated = len(matrix) < len(set(f for m in members for f in feature_of(m)))
    if truncated:
        ax.text(
            0.99,
            1.06,
            f"showing the {MAX_HEATMAP_ROWS} most frequent {row_label}(s)",
            transform=ax.transAxes,
            ha="right",
            fontsize=7,
            style="italic",
        )
    fig.tight_layout()
    return FigureResult(name, _save(fig, path), True)


# --------------------------------------------------------------------------
# co-occurrence figures
# --------------------------------------------------------------------------


def gene_pairs_figure(
    rows: Sequence[Mapping[str, Any]], title: str, path: Path, top_n: int = 25
) -> FigureResult:
    """Bar chart of the most frequent co-occurring gene pairs."""
    name = path.stem
    if len(rows) < MIN_ROWS:
        return FigureResult(name, None, False, f"only {len(rows)} gene pair(s) observed")

    ordered = sorted(rows, key=lambda r: -float(r.get("n_isolates_with_both") or 0))
    shown = ordered[:top_n]
    labels = [f"{r['Gene_A']} + {r['Gene_B']}" for r in shown][::-1]
    counts = [float(r.get("n_isolates_with_both") or 0) for r in shown][::-1]

    height = max(4.0, 0.3 * len(shown))
    plt = _pyplot()
    fig, ax = plt.subplots(figsize=(10, height))
    ax.barh(labels, counts, color="#4C72B0")
    ax.set_xlabel("Isolates carrying both genes (n)")
    ax.set_title(title)
    ax.grid(axis="x", linestyle=":", alpha=0.4)
    for index, value in enumerate(counts):
        ax.text(value, index, f" {value:g}", va="center", fontsize=7)
    fig.text(0.5, -0.01, ASSOCIATION_NOTE, ha="center", fontsize=7, style="italic")
    fig.tight_layout()
    return FigureResult(name, _save(fig, path), True)


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------


def build_all(
    out_dir: Path,
    antibiotic: str,
    gene_rows: Sequence[Mapping[str, Any]],
    gene_freq_rows: Sequence[Mapping[str, Any]],
    mechanism_rows: Sequence[Mapping[str, Any]],
    pair_rows: Sequence[Mapping[str, Any]],
    members: Sequence[Any],
) -> List[FigureResult]:
    """Build every figure for one antibiotic, plus the shared heatmaps.

    Heatmaps are antibiotic-independent (they show gene and mechanism
    carriage across the cohort) and are therefore built once.
    """
    out_dir = Path(out_dir)
    label = antibiotic.capitalize()
    results: List[FigureResult] = []

    results.append(
        gene_frequency_figure(
            gene_freq_rows,
            f"Pilot-100: {label} gene detection frequency",
            out_dir / f"pilot100_{antibiotic}_genes.png",
        )
    )
    results.append(
        mechanism_frequency_figure(
            mechanism_rows,
            f"Pilot-100: {label} resistance mechanism distribution",
            out_dir / f"pilot100_{antibiotic}_mechanisms.png",
        )
    )
    results.append(
        gene_pairs_figure(
            pair_rows,
            f"Pilot-100: {label} most frequent gene pairs",
            out_dir / f"pilot100_{antibiotic}_gene_pairs.png",
        )
    )
    return results


def build_shared(
    out_dir: Path, members: Sequence[Any], families: Mapping[str, str]
) -> List[FigureResult]:
    """Build the two antibiotic-independent heatmaps."""
    from . import mechanisms as mech_mod

    out_dir = Path(out_dir)
    return [
        heatmap_figure(
            members,
            lambda m: sorted(m.genes),
            "Pilot-100: AMR gene presence/absence (AMRFinderPlus)",
            out_dir / "pilot100_gene_heatmap.png",
            row_label="Gene",
        ),
        heatmap_figure(
            members,
            lambda m: sorted(
                mech_mod.functional_category(g, m.gene_annotations.get(g), families)
                for g in m.genes
            ),
            "Pilot-100: resistance mechanism presence/absence",
            out_dir / "pilot100_mechanism_heatmap.png",
            row_label="Mechanism",
        ),
    ]
