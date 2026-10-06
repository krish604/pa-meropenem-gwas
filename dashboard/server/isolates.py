"""The joined per-isolate table (DESIGN §5.2).

One row per `sample_id` in the manifest — or, where no manifest is readable, in
the union of the per-sample tables, with `membership_source` naming which.

**A sample present in one table and absent from another still gets a row.** The
cells with no source read `not assessed`, not `0` and not blank. This is rule 9
of `docs/scientific_rules.md` on the client, and the distinction is load-bearing
in three places:

| column | a present table with no row | a missing table |
|---|---|---|
| `amr_genes` | `[]` — a measurement of zero | `not assessed` |
| `n_virulence` | `0` | `not assessed` |
| `n_snv` / `n_indel` | `0` (the provenance file writes every key at zero) | `not assessed` |

The rule from the pipeline's own docstring, restated: **a table that exists and
holds no row for a sample is a measurement of zero; a table that does not exist
is not a measurement.**

`oprd_state` and `oprd_verdict` come from `oprd.py` and are `not_assessed`
unless a verdict source exists — never inferred from the master table's
`chromosomal_mutation`, which is the derivation that would manufacture the
study's central negative (DESIGN §4.3).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from papipeline.stages.report_tables import TableRead
from papipeline.viz import OPRD_STATE_ORDER

#: The literals this module renders. One place, so no page can invent another.
NOT_ASSESSED = "not assessed"
NOT_PRODUCED = "not produced"
NOT_REPORTED = "not reported"

#: `data_contract.md:259` — the sentinel the master-table builder refuses, so
#: one reaching here is a discrepancy worth surfacing rather than rendering.
SENTINEL_LINEAGE = "unknown"

#: The isolate view's columns, in order. Also the sortable/filterable set.
COLUMNS: Tuple[str, ...] = (
    "sample_id",
    "st",
    "amr_genes",
    "amr_point_mutations",
    "n_virulence",
    "oprd_state",
    "oprd_verdict",
    "phenotype_sir",
    "lineage",
    "n_snv",
    "n_indel",
)

#: Substring search runs over these.
SEARCHABLE: Tuple[str, ...] = (
    "sample_id",
    "st",
    "amr_genes",
    "amr_point_mutations",
    "oprd_verdict",
    "phenotype_sir",
    "lineage",
)


@dataclass
class _Index:
    """Per-table lookup structures, built once and reused per request.

    At 900 rows and six contributing tables a per-request scan is six passes
    per row; this is one pass per table and a dict lookup per row, which is
    what makes `/api/isolates` instant (UI-D5).
    """

    amr_genes: Dict[str, List[str]] = field(default_factory=dict)
    amr_mutations: Dict[str, List[str]] = field(default_factory=dict)
    amr_present: bool = False
    virulence: Dict[str, int] = field(default_factory=dict)
    virulence_present: bool = False
    mlst: Dict[str, Optional[str]] = field(default_factory=dict)
    mlst_status: Dict[str, Optional[str]] = field(default_factory=dict)
    mlst_present: bool = False
    phenotype: Dict[str, str] = field(default_factory=dict)
    phenotype_present: bool = False
    lineage: Dict[str, Optional[str]] = field(default_factory=dict)
    lineage_present: bool = False
    variants: Dict[str, Tuple[int, int]] = field(default_factory=dict)
    variants_present: bool = False


def _present(read: TableRead) -> bool:
    return bool(read is not None and read.present)


def build_index(
    *,
    mlst: TableRead,
    amr: TableRead,
    virulence: TableRead,
    phenotype: TableRead,
    master: TableRead,
    tree_metadata: TableRead,
    variants_provenance: Optional[Mapping[str, Any]],
) -> _Index:
    """One pass per contributing table.

    Every column read here comes from a table the pipeline declares. Nothing is
    derived from a file that does not exist.
    """
    index = _Index()

    # -- stage 3: MLST ----------------------------------------------------
    index.mlst_present = _present(mlst)
    if index.mlst_present:
        for row in mlst.rows:
            sample = row.get("sample_id")
            if not sample:
                continue
            status = row.get("MLST_status")
            st = row.get("ST")
            # `no_call` means the stage looked and could not call; the ST is
            # then not an answer. Distinct from a missing row.
            index.mlst_status[str(sample)] = status
            index.mlst[str(sample)] = (
                None if str(status).strip() == "no_call" else st
            )

    # -- stage 4: AMR -----------------------------------------------------
    index.amr_present = _present(amr)
    if index.amr_present:
        genes: Dict[str, List[str]] = {}
        mutations: Dict[str, List[str]] = {}
        for row in amr.rows:
            sample = row.get("sample_id")
            if not sample:
                continue
            sample = str(sample)
            gene = row.get("gene")
            determinant = row.get("determinant")
            variant = row.get("variant")
            # `variant` non-empty means the tool says this element is mutated.
            # It does NOT mean disruptive - `determinant_type` is `AMR` for a
            # point mutation and for an intact gene alike, so `variant` is
            # what discriminates. The column is labelled accordingly.
            if variant is not None and str(variant).strip():
                mutations.setdefault(sample, []).append(
                    f"{gene or determinant}:{variant}"
                )
            name = gene or determinant
            if name is not None and str(name).strip():
                bucket = genes.setdefault(sample, [])
                if name not in bucket:
                    bucket.append(str(name))
        index.amr_genes = genes
        index.amr_mutations = mutations

    # -- stage 5: virulence ----------------------------------------------
    index.virulence_present = _present(virulence)
    if index.virulence_present:
        counts: Dict[str, int] = {}
        for row in virulence.rows:
            sample = row.get("sample_id")
            if not sample:
                continue
            counts[str(sample)] = counts.get(str(sample), 0) + 1
        index.virulence = counts

    # -- stage 11: phenotype ---------------------------------------------
    index.phenotype_present = _present(phenotype)
    if index.phenotype_present:
        calls: Dict[str, str] = {}
        for row in phenotype.rows:
            sample = row.get("sample_id")
            value = row.get("phenotype")
            if not sample or value is None:
                continue
            calls[str(sample)] = str(value)
        index.phenotype = calls

    # -- lineage: master table first, else tree_metadata ------------------
    # The master table REFUSES a sentinel lineage
    # (`data_contract.md:259`), so one reaching here is surfaced as a
    # discrepancy rather than rendered.
    if _present(master):
        index.lineage_present = True
        for row in master.rows:
            sample = row.get("sample_id")
            if sample:
                value = row.get("lineage")
                index.lineage[str(sample)] = (
                    NOT_ASSESSED
                    if value is None or str(value).strip() in ("", SENTINEL_LINEAGE)
                    else str(value)
                )
    elif _present(tree_metadata):
        for row in tree_metadata.rows:
            sample = row.get("sample_id")
            if sample:
                value = row.get("lineage_label")
                index.lineage[str(sample)] = (
                    None if value is None or str(value).strip() == ""
                    else str(value)
                )

    # -- stage 6 accounting ----------------------------------------------
    # Every key is present at zero when the file exists, so `0` here is a real
    # measurement. Absent file -> not_assessed, because "called nothing" and
    # "never screened" must not look alike.
    if isinstance(variants_provenance, Mapping):
        per_isolate = variants_provenance.get("per_isolate")
        if isinstance(per_isolate, Sequence):
            index.variants_present = True
            for entry in per_isolate:
                if not isinstance(entry, Mapping):
                    continue
                sample = entry.get("sample_id")
                if not sample:
                    continue
                snv = entry.get("n_snv")
                indel = entry.get("n_indel")
                index.variants[str(sample)] = (
                    int(snv) if isinstance(snv, (int, float)) else None,
                    int(indel) if isinstance(indel, (int, float)) else None,
                )
    return index


def cohort(
    manifest_samples: Optional[Sequence[str]],
    tables: Mapping[str, TableRead],
) -> Tuple[List[str], str]:
    """The cohort, and which artefact defined it.

    `manifest` when the run's manifest lists the samples; `union_of_tables`
    when no manifest is readable and the cohort is the union of the per-sample
    tables. Never silently narrowed: a sample missing from one table still gets
    a row.
    """
    if manifest_samples:
        return sorted({str(s) for s in manifest_samples}), "manifest"
    seen: set = set()
    for read in tables.values():
        if not _present(read):
            continue
        for row in read.rows:
            sample = row.get("sample_id")
            if sample:
                seen.add(str(sample))
    return sorted(seen), "union_of_tables"


def row_for(sample: str, index: _Index, oprd: Optional[Mapping[str, Any]] = None) -> Dict[str, Any]:
    """One isolate's row, with every absent cell named rather than defaulted."""
    oprd_row = (oprd or {}).get(sample) or {}
    oprd_state = oprd_row.get("display_state")
    oprd_verdict = oprd_row.get("verdict")

    if sample not in index.mlst:
        st: Any = NOT_ASSESSED
    else:
        st = index.mlst[sample] if index.mlst[sample] is not None else NOT_ASSESSED

    amr_genes = index.amr_genes.get(sample)
    if amr_genes is None:
        genes: Any = NOT_ASSESSED
    else:
        genes = list(amr_genes)

    mutations = index.amr_mutations.get(sample)
    if mutations is None:
        point_mutations: Any = NOT_ASSESSED
    else:
        point_mutations = list(mutations)

    if not index.virulence_present:
        n_virulence: Any = NOT_ASSESSED
    else:
        n_virulence = index.virulence.get(sample, 0)

    # A missing phenotype table and a present table with no row for this sample
    # both read `not assessed`: the run recorded no SIR call, and `None` would
    # render as a blank cell, which is exactly the plausible absence UI-D2
    # refuses. (When the table is present and holds a row, the value is verbatim,
    # and `ND` is a recorded value.)
    if not index.phenotype_present:
        phenotype: Any = NOT_ASSESSED
    else:
        phenotype = index.phenotype.get(sample, NOT_ASSESSED)
    lineage = index.lineage.get(sample, NOT_ASSESSED) if index.lineage_present else NOT_ASSESSED

    if index.variants_present and sample in index.variants:
        snv, indel = index.variants[sample]
        n_snv: Any = snv if snv is not None else NOT_ASSESSED
        n_indel: Any = indel if indel is not None else NOT_ASSESSED
    else:
        n_snv = NOT_ASSESSED
        n_indel = NOT_ASSESSED

    return {
        "sample_id": sample,
        "st": st,
        "amr_genes": genes,
        "amr_point_mutations": point_mutations,
        "n_virulence": n_virulence,
        "oprd_state": oprd_state if oprd_state in OPRD_STATE_ORDER else "not_assessed",
        "oprd_verdict": oprd_verdict,
        "phenotype_sir": phenotype,
        "lineage": lineage,
        "n_snv": n_snv,
        "n_indel": n_indel,
        # The three derived filters DESIGN §5.2 promises. They are computed
        # here, server-side, from the same index the row is built from, so a
        # filter can never disagree with the cells beside it.
        "has_amr_gene": bool(isinstance(genes, list) and genes),
        "has_point_mutation": bool(isinstance(point_mutations, list) and point_mutations),
        "virulence_min": n_virulence if isinstance(n_virulence, int) else None,
    }


def contributing_sources(index: _Index) -> Dict[str, Any]:
    """Which tables actually answered, and what each one means for the cells.

    Served beside the rows so a client can say "no virulence row" and "no
    virulence table" apart without re-deriving it, and so `n_virulence: 0` is
    visibly a measurement rather than an absent source.
    """
    return {
        "mlst": index.mlst_present,
        "amr": index.amr_present,
        "virulence": index.virulence_present,
        "phenotype": index.phenotype_present,
        "lineage": index.lineage_present,
        "variants_provenance": index.variants_present,
    }


def absent_behaviour() -> Dict[str, str]:
    """Per column, what `not assessed` means and what a value means.

    The distinction is the whole point of the column set, so it is stated in
    the response rather than left to each page.
    """
    return {
        "amr_genes": (
            "[] is a measurement of zero from a present 04_amr.tsv; "
            "'not assessed' means the table is not there"
        ),
        "amr_point_mutations": (
            "rows whose `variant` column is non-empty, i.e. the tool reported "
            "the element as mutated - NOT that it is disruptive"
        ),
        "n_virulence": (
            "0 when 08_virulence.tsv is present and holds no row for this "
            "sample; 'not assessed' when the table is absent"
        ),
        "oprd_state": (
            "'not_assessed' unless a verdict source exists. Never inferred from "
            "the master table's chromosomal_mutation column."
        ),
        "phenotype_sir": (
            "R/I/S/SDD/ND for the run's antibiotic. 'ND' is a recorded value, "
            "distinct from an absent row. Never converted to an MIC."
        ),
        "n_snv": (
            "0 is a real measurement when variants_provenance.json exists - it "
            "writes every key at zero. Absent file reads 'not assessed', because "
            "'called nothing' and 'never screened' must not look alike."
        ),
        "lineage": (
            "from 15_master_table.tsv, else tree_metadata.tsv. The master table "
            "refuses a sentinel lineage, so a sentinel here is surfaced rather "
            "than rendered."
        ),
    }


def detail_tables(
    sample: str,
    *,
    reads: Mapping[str, TableRead],
) -> List[Dict[str, Any]]:
    """Every contributing table for one isolate, absent ones included.

    So the page can say which tables had nothing for this isolate and which
    tables were not there at all — a distinction `[]` cannot carry.
    """
    out: List[Dict[str, Any]] = []
    for key, read in reads.items():
        rows = [dict(r) for r in read.rows if str(r.get("sample_id")) == sample] if _present(read) else []
        out.append(
            {
                "key": key,
                "present": _present(read),
                "reason": "" if _present(read) else read.reason,
                "path": str(read.path),
                "rows": rows,
            }
        )
    return out


__all__ = [
    "COLUMNS",
    "NOT_ASSESSED",
    "NOT_PRODUCED",
    "NOT_REPORTED",
    "OPRD_STATE_ORDER",
    "SEARCHABLE",
    "SENTINEL_LINEAGE",
    "absent_behaviour",
    "build_index",
    "contributing_sources",
    "cohort",
    "detail_tables",
    "row_for",
]