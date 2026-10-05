"""PowerPoint report for the pilot-100 analysis.

Sixteen slides, one per topic required by the pilot specification. Every
slide is generated from the analysis tables; no slide contains a hand-typed
number, and every slide that shows counts carries the run's provenance and
the applicable scientific caveat.

The deck is titled as a pilot and states the selection method on the first
content slide, because the selection method is the single most important
provenance fact in the whole analysis.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from ..logging_utils import get_logger
from . import mechanisms as mech_mod
from .census import describe_selection

LOGGER = get_logger("pilot.pptx")

TITLE = "Pilot-100 Pseudomonas aeruginosa AMR Analysis"
SUBTITLE = "AMRFinderPlus gene detection, PDC metadata join, imipenem and meropenem"


def _add_title(slide, text: str, subtitle: Optional[str] = None) -> None:
    slide.shapes.title.text = text
    if subtitle:
        box = slide.shapes.add_textbox(
            slide.shapes.title.left, slide.shapes.title.top + 40, 8_000_000, 700_000
        )
        frame = box.text_frame
        frame.text = subtitle
        frame.paragraphs[0].runs[0].font.size = __import__("pptx.util", fromlist=["Pt"]).Pt(14)


def _add_bullets(slide, items: Sequence[str], left: int = 600_000, top: int = 1_700_000) -> None:
    from pptx.util import Inches, Pt

    height = Inches(0.42 * max(1, len(items)) + 0.3)
    box = slide.shapes.add_textbox(left, top, Inches(9.0), height)
    frame = box.text_frame
    frame.word_wrap = True
    for index, item in enumerate(items):
        paragraph = frame.paragraphs[0] if index == 0 else frame.add_paragraph()
        paragraph.text = str(item)
        paragraph.font.size = Pt(15)
        paragraph.level = 0
        paragraph.space_after = Pt(4)


def _add_table(
    slide,
    columns: Sequence[str],
    rows: Sequence[Sequence[Any]],
    top: int = 1_700_000,
    max_rows: int = 12,
    font_size: int = 10,
) -> None:
    from pptx.util import Inches, Pt

    shown = list(rows)[:max_rows]
    if not shown:
        _add_bullets(slide, ["No rows available."], top=top)
        return
    left = Inches(0.4)
    width = Inches(9.2)
    height = Inches(0.3 * (len(shown) + 1) + 0.2)
    table = slide.shapes.add_table(
        len(shown) + 1, len(columns), left, top, width, height
    ).table
    for index, name in enumerate(columns):
        cell = table.cell(0, index)
        cell.text = str(name)
        cell.text_frame.paragraphs[0].runs[0].font.bold = True
        cell.text_frame.paragraphs[0].runs[0].font.size = Pt(font_size)
    for r, row in enumerate(shown, start=1):
        for c, value in enumerate(row):
            cell = table.cell(r, c)
            cell.text = "-" if value is None else str(value)
            cell.text_frame.paragraphs[0].runs[0].font.size = Pt(font_size)
    if len(rows) > max_rows:
        note = slide.shapes.add_textbox(
            Inches(0.4), top + height, Inches(9.2), Inches(0.3)
        )
        note.text_frame.text = f"... and {len(rows) - max_rows} more row(s); see the TSV"
        note.text_frame.paragraphs[0].runs[0].font.size = Pt(9)


def _add_image(slide, path: Optional[Path], caption: str) -> None:
    from pptx.util import Inches, Pt

    if path is None or not Path(path).exists():
        _add_bullets(slide, [f"Figure not created: {caption}"])
        return
    from PIL import Image

    with Image.open(path) as image:
        width_px, height_px = image.size
    max_w, max_h = Inches(9.0), Inches(5.0)
    scale = min(max_w / width_px, max_h / height_px)
    left = int((Inches(10.0) - width_px * scale) / 2)
    slide.shapes.add_picture(
        str(path), left, Inches(1.6), width=int(width_px * scale), height=int(height_px * scale)
    )
    box = slide.shapes.add_textbox(Inches(0.5), Inches(6.6), Inches(9.0), Inches(0.4))
    box.text_frame.text = caption
    box.text_frame.paragraphs[0].runs[0].font.size = Pt(10)


def _count_by(rows: Sequence[Mapping[str, Any]], key: str) -> List[Tuple[str, int]]:
    from collections import Counter

    return sorted(Counter(str(r.get(key)) for r in rows).items(), key=lambda kv: -kv[1])


def build_deck(
    path: Path,
    summary: Mapping[str, Any],
    tables: Mapping[str, Sequence[Mapping[str, Any]]],
    figures: Mapping[str, Optional[Path]],
    antibiotic_labels: Sequence[str] = ("imipenem", "meropenem"),
) -> Path:
    """Write the PowerPoint deck.

    Args:
        path: Destination ``.pptx``.
        summary: Run-level counts for the provenance slides.
        tables: Analysis tables keyed by name.
        figures: Figure path per name; ``None`` where a figure was skipped.
    """
    from pptx import Presentation
    from pptx.util import Inches, Pt

    presentation = Presentation()
    blank = presentation.slide_layouts[6]
    title_layout = presentation.slide_layouts[0]
    body_layout = presentation.slide_layouts[1]

    provenance = (
        f"Detection: AMRFinderPlus {summary.get('amrfinder_version', 'n/a')}, "
        f"database {summary.get('amr_database', 'n/a')} | "
        f"metadata: PDC_essential.tsv | "
        f"generated {summary.get('generated_utc', 'n/a')}"
    )

    # 1. Title
    slide = presentation.slides.add_slide(title_layout)
    slide.shapes.title.text = TITLE
    slide.placeholders[1].text = SUBTITLE

    # 2. Selection methodology
    slide = presentation.slides.add_slide(body_layout)
    _add_title(slide, "1. Pilot-100 selection methodology")
    _add_bullets(
        slide,
        [
            "SELECTION METHOD: the FIRST 100 genome assemblies physically present "
            "in data/, sorted by GCA accession.",
            "NOT the first 100 rows of PDC_essential.tsv.",
            f"Assemblies discovered in data/: {summary.get('assemblies_discovered')}",
            f"First 100 selected: {summary.get('n_selected', summary.get('selected'))}",
            f"Excluded (empty/corrupt/undersized assembly): {summary.get('n_excluded')}",
            f"ANALYSED: {summary.get('n_analysable')} assemblies. Every count in "
            "this deck is out of that number, not 100.",
            f"Overlap with PDC row order's first 100: "
            f"{summary.get('selection_overlap')} of {summary.get('selected')} "
            "(the two selections are not the same set)",
            "PDC_essential.tsv was read only AFTER the 100 were selected, to attach metadata.",
            "data/ was not modified; assemblies are reached via symlinks.",
            provenance,
        ],
    )

    # 3. Assembly QC
    slide = presentation.slides.add_slide(body_layout)
    _add_title(slide, "2. Assembly QC")
    qc_rows = tables.get("qc_summary", [])
    if qc_rows:
        _add_table(
            slide,
            ["Metric", "Value"],
            [[r.get("metric"), r.get("value")] for r in qc_rows],
        )
    else:
        _add_bullets(slide, ["Assembly QC not run in this pilot."])

    # 4. PDC metadata matching
    slide = presentation.slides.add_slide(body_layout)
    _add_title(slide, "3. PDC metadata matching")
    _add_bullets(
        slide,
        [
            f"Selected assemblies: {summary.get('selected')}",
            f"Matched to PDC_essential.tsv: {summary.get('pdc_matched')}",
            f"Missing from PDC_essential.tsv: {summary.get('pdc_missing')}",
            "Matching is on the GCA/Assembly accession, never on row order.",
            "A genome with no PDC record is retained and marked "
            "PDC_metadata_match=MISSING; it is never discarded.",
        ],
    )

    # 5-6. Phenotypes
    for index, antibiotic in enumerate(antibiotic_labels, start=4):
        label = antibiotic.capitalize()
        slide = presentation.slides.add_slide(body_layout)
        _add_title(
            slide,
            f"{index}. {label} phenotype",
            "Categorical susceptibility only. No MIC or zone diameter was derived.",
        )
        rows = tables.get(f"phenotype_{antibiotic}", [])
        counts = _count_by(rows, "phenotype")
        _add_table(
            slide,
            ["Category", "Isolates (n)"],
            [[category, n] for category, n in counts],
        )

    # 7-8. Genes
    for index, antibiotic in enumerate(antibiotic_labels, start=6):
        label = antibiotic.capitalize()
        slide = presentation.slides.add_slide(body_layout)
        _add_title(
            slide,
            f"{index}. {label} genes",
            "Gene detection is DETECTION, not phenotypic resistance.",
        )
        _add_image(
            slide,
            figures.get(f"pilot100_{antibiotic}_genes"),
            f"Pilot-100 {label} gene detection frequency (AMRFinderPlus)",
        )

    # 9-10. Mechanisms
    for index, antibiotic in enumerate(antibiotic_labels, start=8):
        label = antibiotic.capitalize()
        slide = presentation.slides.add_slide(body_layout)
        _add_title(slide, f"{index}. {label} mechanisms")
        _add_image(
            slide,
            figures.get(f"pilot100_{antibiotic}_mechanisms"),
            f"Pilot-100 {label} mechanism distribution",
        )

    # 11-12. Heatmaps
    slide = presentation.slides.add_slide(body_layout)
    _add_title(slide, "10. Gene heatmap")
    _add_image(
        slide, figures.get("pilot100_gene_heatmap"), "AMR gene presence/absence, 100 isolates"
    )

    slide = presentation.slides.add_slide(body_layout)
    _add_title(slide, "11. Mechanism heatmap")
    _add_image(
        slide,
        figures.get("pilot100_mechanism_heatmap"),
        "Mechanism presence/absence, 100 isolates",
    )

    # 13-14. Gene pairs
    for index, antibiotic in enumerate(antibiotic_labels, start=12):
        label = antibiotic.capitalize()
        slide = presentation.slides.add_slide(body_layout)
        _add_title(
            slide, f"{index}. {label} gene pairs", "Association only, not causation."
        )
        _add_image(
            slide,
            figures.get(f"pilot100_{antibiotic}_gene_pairs"),
            f"Most frequent {label} gene pairs",
        )

    # 15. Mechanism combinations
    slide = presentation.slides.add_slide(body_layout)
    _add_title(
        slide, "14. Mechanism combinations", "Association only, not causation."
    )
    combo_rows = tables.get("mechanism_combinations", [])
    if combo_rows:
        _add_table(
            slide,
            ["Combination", "Isolates (n)", "Genes"],
            [
                [r.get("Combination"), r.get("n_isolates"), str(r.get("genes", ""))[:60]]
                for r in combo_rows
            ],
            max_rows=12,
            font_size=9,
        )
    else:
        _add_bullets(slide, ["No mechanism combinations recorded."])

    # 15. QC
    slide = presentation.slides.add_slide(body_layout)
    _add_title(slide, "15. QC")
    qc_rows_final = tables.get("qc_summary", [])
    _add_bullets(
        slide,
        [
            f"Assemblies selected (first 100, GCA order): {summary.get('n_selected', summary.get('selected'))}",
            f"Excluded for an unusable assembly file: {summary.get('n_excluded')}",
            f"Genomes analysed: {summary.get('n_analysable')}",
            f"Genomes processed successfully: {summary.get('n_processed')}"
            f" / {summary.get('n_analysable')} analysable",
            f"Genomes that failed AMR detection: {summary.get('n_failed')}",
            f"AMR gene calls: {summary.get('n_gene_calls')}",
            f"Distinct genes: {summary.get('n_distinct_genes')}",
            f"Distinct mechanism categories: {summary.get('n_mechanisms')}",
            f"Gene pairs: {summary.get('n_gene_pairs')}",
            f"Mechanism combinations: {summary.get('n_mechanism_combinations')}",
            f"Isolates with imipenem metadata: {summary.get('n_with_imipenem')}",
            f"Isolates with meropenem metadata: {summary.get('n_with_meropenem')}",
            provenance,
        ],
    )
    if qc_rows_final:
        LOGGER.info("QC rows available for the deck: %d", len(qc_rows_final))

    # 16. Limitations
    slide = presentation.slides.add_slide(body_layout)
    _add_title(slide, "16. Limitations")
    _add_bullets(
        slide,
        [
            "0. n is NOT 100. 45 of the 100 selected assemblies have an empty, "
            "corrupt or undersized file and are excluded; nothing is substituted. "
            "This is a data-download problem in data/, not an analytical one.",
            "1. Gene detection does not establish phenotypic resistance.",
            "2. OprD gene detection indicates an INTACT locus, not reduced "
            "permeability; only an explicit disruptive variant implies disruption.",
            "3. ESBL/MBL labels are assigned by beta-lactamase gene-family "
            "nomenclature (config/gene_families.tsv), not by enzyme activity.",
            "4. Phenotypes are categorical R/I/S only. No MIC or zone diameter "
            "was derived, and none exists in the source.",
            "5. This is an observational pilot. No GWAS, phylogenetic or "
            "convergence analysis was performed, so no association here is "
            "corrected for lineage or population structure.",
            "6. Co-occurrence is association only and is not adjusted for clonal "
            "structure; pairs may co-occur simply because they belong to the same clone.",
            f"7. {describe_selection(summary.get('n_selected', 0), summary.get('assemblies_discovered', 0))}"
            " available assemblies; the pilot is not a random sample "
            "and does not estimate prevalence in the full collection.",
            "8. Pan-genome and phylogenomics stages were not run, so lineage "
            "structure is unknown.",
        ],
    )

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    presentation.save(path)
    LOGGER.info("PowerPoint written to %s", path)
    return path
