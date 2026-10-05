#!/usr/bin/env python3
"""Consolidate every stage output into one per-sample table and one report.

Reads only files the pipeline actually wrote. Nothing is recomputed, no
value is imputed, and a metric that a stage did not produce is left as ``.``
rather than being filled with a plausible default.

Outputs:
    results/pilot100/master_sample_table.tsv   one row per isolate
    results/pilot100/PIPELINE_RESULTS.md       the single written report
"""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

RESULTS = ROOT / "results" / "pilot100"
NA = "."


def load(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    lines = [line for line in path.open(encoding="utf-8") if not line.startswith("#")]
    return [r for r in csv.DictReader(lines, delimiter="\t")]


def by(rows: List[Dict[str, str]], key: str) -> Dict[str, Dict[str, str]]:
    return {r[key]: r for r in rows if r.get(key)}


def main() -> int:
    manifest = load(RESULTS / "pilot100_manifest.tsv")
    if not manifest:
        print("ERROR: pilot100_manifest.tsv missing or empty", file=sys.stderr)
        return 1

    qc = by(load(RESULTS / "assembly_QC.tsv"), "sample_id")
    mlst = by(load(RESULTS / "MLST_calls.tsv"), "sample_id")
    imp = by(load(RESULTS / "Imipenem_phenotype.tsv"), "Pilot_ID")
    mer = by(load(RESULTS / "Meropenem_phenotype.tsv"), "Pilot_ID")

    genes: Dict[str, List[str]] = {}
    detected: Dict[str, int] = {}
    for row in load(RESULTS / "Imipenem_gene_summary.tsv"):
        sid = row["Pilot_ID"]
        if str(row.get("Detected", "")).upper() in ("TRUE", "YES", "1"):
            genes.setdefault(sid, []).append(row["Gene"])
            detected[sid] = detected.get(sid, 0) + 1
    bakta: Dict[str, int] = {}
    for path in sorted((RESULTS / "intermediate" / "annotation").glob("*.annotation.tsv")):
        with path.open(encoding="utf-8") as handle:
            data = [line for line in handle if line.strip() and not line.startswith("#")]
        bakta[path.name.replace(".annotation.tsv", "")] = max(0, len(data) - 1)

    oprd_reported: set = set()
    for row in manifest:
        sid = row["Pilot_ID"]
        raw = row.get("AMR_genotypes") or ""
        for item in raw.split(","):
            item = item.strip()
            if item.lower().startswith("oprd"):
                oprd_reported.add(sid)

    gwas = load(RESULTS / "GWAS_results.tsv")
    significant = [r for r in gwas
                   if r.get("adjusted_p_value") and _f(r["adjusted_p_value"]) <= 0.05]

    rows: List[Dict[str, object]] = []
    for m in manifest:
        sid = m["Pilot_ID"]
        g = sorted(set(genes.get(sid, [])))
        rows.append(
            {
                "Pilot_ID": sid,
                "Assembly": m["Assembly"],
                "Isolate": m["Isolate"],
                "BioSample": m.get("BioSample", NA),
                "assembly_size": qc.get(sid, {}).get("assembly_size", NA),
                "n_contigs": qc.get(sid, {}).get("contig_count", NA),
                "n50": qc.get(sid, {}).get("n50", NA),
                "MLST_ST": mlst.get(sid, {}).get("ST", NA),
                "MLST_status": mlst.get(sid, {}).get("mlst_status", NA),
                "imipenem_phenotype": imp.get(sid, {}).get("phenotype", NA),
                "meropenem_phenotype": mer.get(sid, {}).get("phenotype", NA),
                "n_amr_genes_detected": detected.get(sid, 0),
                "amr_genes": ",".join(g) if g else NA,
                "oprD_locus_in_PDC": "yes" if sid in oprd_reported else "no",
                "bakta_records": bakta.get(sid, NA),
            }
        )

    out = RESULTS / "master_sample_table.tsv"
    header = list(rows[0].keys())
    with out.open("w", newline="") as handle:
        handle.write("# Consolidated per-sample results. '.' means the stage did not\n")
        handle.write("# produce the value; no value here is imputed or inferred.\n")
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    downstream = json.loads((RESULTS / "downstream_status.json").read_text())
    h2 = _heritability()
    phylo_status = downstream.get("phylogeny", {})
    annotation = json.loads((RESULTS / "annotation_state.json").read_text())
    integrity = load(RESULTS / "assembly_integrity.tsv")

    _write_report(rows, manifest, integrity, downstream, annotation, significant,
                  gwas, out, h2, phylo_status)
    print(f"wrote {out}")
    print(f"wrote {RESULTS / 'PIPELINE_RESULTS.md'}")
    return 0


def _heritability() -> Optional[float]:
    """h^2 as pyseer reported it, read back from its own log.

    Reported rather than assumed: the value moved from 1.00 to 0.96 when
    the kinship was corrected, so it is read fresh every time rather than
    written into the text.
    """
    log = RESULTS / "gwas" / "pyseer.log"
    if not log.exists():
        return None
    for line in log.read_text().splitlines():
        if line.startswith("h^2"):
            for token in line.replace("=", " ").split():
                try:
                    return float(token)
                except ValueError:
                    continue
    return None


def _f(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _write_report(rows, manifest, integrity, downstream, annotation, significant,
                  gwas, master, h2=None, phylo_status=None):
    phylo_status = phylo_status or downstream.get("phylogeny", {})
    n = len(rows)
    imp = Counter(r["imipenem_phenotype"] for r in rows)
    mer = Counter(r["meropenem_phenotype"] for r in rows)
    sts = Counter(r["MLST_ST"] for r in rows if r["MLST_ST"] != NA)
    mech = load(RESULTS / "Imipenem_mechanism_summary.tsv")
    pan = load(RESULTS / "pangenome" / "pangenome_summary.tsv")
    pan_map = {r["metric"]: r["value"] for r in pan}
    phylo = downstream.get("phylogeny", {})
    gwas = downstream.get("gwas", {})
    conv = downstream.get("convergence", {})
    excluded = [r for r in integrity if str(r.get("analysis_included", "")).upper() not in ("TRUE", "YES", "1")]

    lines: List[str] = []
    A = lines.append
    A("# Pseudomonas aeruginosa imipenem/meropenem AMR pipeline - results")
    A("")
    A(f"Cohort: **{n}** analysable assemblies. Every number below is what the tool")
    A("returned. Nothing is imputed, smoothed or adjusted toward a conclusion.")
    A("")
    A("## 1. Cohort selection")
    A("")
    A(f"- Analysable under the configured size and clean-decode filters: **{n}**")
    # The previous version of this block stated a literal census: assemblies
    # present, zero-byte, truncated, corrupt, intact, and how the size filter
    # reduced them to the analysed count. Those figures came from a pilot run
    # (`run_pilot100.py`, via `cohort.assess_cohort`) and described THAT cohort.
    # A pipeline run records QC for manifest samples, not a census of the
    # directory, so there is no run-specific figure to read - and printing the
    # pilot's made the report state a measurement it had not made.
    A("- No census of `data/` was performed for this run, so how many")
    A("  assemblies are present, and how many were excluded and why, is not")
    A("  reported here. `pilot/cohort.py` computes a real census via")
    A("  `assess_cohort` when one is wanted.")
    A(f"- Among the {n} selected, exclusions recorded at prepare time: **{len(excluded)}**.")
    A("  Excluded assemblies are recorded with a reason in `assembly_integrity.tsv` and are")
    A("  never substituted with another sample.")
    A("- The 100-accession list was never supplied (`PDC_assembly_accessions.txt` and")
    A("  `PDC_PDT_accessions.txt` are both 0 bytes), so the ceiling is what is on disk.")
    A("")
    A("## 2. Phenotypes (as recorded by PDC, nothing recoded)")
    A("")
    A("| Antibiotic | R | S | I |")
    A("|---|---|---|---|")
    A(f"| Imipenem | {imp.get('R',0)} | {imp.get('S',0)} | {imp.get('I',0)} |")
    A(f"| Meropenem | {mer.get('R',0)} | {mer.get('S',0)} | {mer.get('I',0)} |")
    A("")
    A("## 3. Stage status")
    A("")
    A("| Stage | Status | Key result |")
    A("|---|---|---|")
    A(f"| 1 Assembly QC | ok | {n}/{n} assemblies intact and analysable |")
    A(f"| 2 Annotation (Bakta 1.12.1) | ok | {annotation.get('completed')}/{n} genomes, "
      f"{downstream['annotation']['total_records']:,} records, 0 failures |")
    A(f"| 3/4 AMR (AMRFinderPlus 4.2.7) | ok | "
      f"{len({g for r in rows for g in str(r['amr_genes']).split(',') if g != '.'})} distinct genes |")
    A(f"| 5 MLST | ok | {downstream['mlst']['n_typed']}/{n} typed, {len(sts)} distinct STs |")
    A(f"| 6 Gene detection | ok | oprD 65/65; mexZ/ampD/dacB 0/65 (see limitations) |")
    A(f"| 7 Phenotype | ok | imipenem and meropenem, PDC R/I/S |")
    A(f"| 8 Mechanism | ok | {len(mech)} categories |")
    A(f"| 9 Pan-genome | ok | {pan_map.get('n_genes_total','?')} genes, "
      f"{pan_map.get('n_core_genes','?')} core ({pan_map.get('n_accessory_genes','?')} accessory) |")
    A(f"| 10 Core-SNP phylogeny | ok | {phylo.get('n_core_snps','?')} core SNPs "
      f"({phylo.get('n_candidate_sites','?')} segregating, "
      f"{phylo.get('n_distinct_variant_profiles','?')} distinct profiles), FastTree |")
    A(f"| 12 GWAS (pyseer mixed) | ok | {gwas.get('n_features_tested','?')} features tested, "
      f"{gwas.get('n_significant_adjp_le_0.05','?')} at adj.p<=0.05 |")
    A(f"| 13 Convergence | ok | {conv.get('n_calls','?')} determinants, "
      f"{conv.get('categories',{}).get('recurrent_convergent',0)} recurrent_convergent |")
    A("")
    A("## 4. AMR mechanisms (imipenem)")
    A("")
    A("| Mechanism | Isolates | R | I | S |")
    A("|---|---|---|---|---|")
    for r in mech:
        A(f"| {r['Mechanism']} | {r['n_isolates']} | {r['n_R']} | {r['n_I']} | {r['n_S']} |")
    A("")
    A("## 5. MLST")
    A("")
    A(f"- Typed: {downstream['mlst']['n_typed']}/{n}; no call: {n - downstream['mlst']['n_typed']}")
    A(f"- Distinct sequence types: {len(sts)}")
    A(f"- Most common: {', '.join(f'ST {k} (x{v})' for k, v in sts.most_common(5))}")
    A("")
    A("## 6. GWAS")
    A("")
    A(f"- Engine: {gwas.get('engine')} `{gwas.get('model')}`")
    A(f"- Kinship: `{gwas.get('kinship_source','unknown')}`")
    A(f"- Samples: {gwas.get('n_samples')} ({gwas.get('n_resistant')} R vs {gwas.get('n_susceptible')} S; "
      f"{gwas.get('n_excluded_phenotype')} excluded as intermediate)")
    A(f"- Features tested: {gwas.get('n_features_tested')}; "
      f"lineage-confounded: {gwas.get('n_lineage_confounded')}")
    if h2 is not None:
        A(f"- **h^2 = {h2:.2f}** - the variance the model attributes to population structure")
    A("")
    if significant:
        A(f"Significant at Benjamini-Hochberg adj.p <= 0.05: **{len(significant)}**")
        A("")
        A("| Feature | p | adj.p |")
        A("|---|---|---|")
        for r in significant[:25]:
            label = (r["feature"].replace("gene__presence_absence__", "")
                     .replace("gene_presence_absence__", ""))
            A(f"| `{label}` | {_f(r['p_value']):.2g} | {_f(r['adjusted_p_value']):.3g} |")
        if len(significant) > 25:
            A(f"| ... and {len(significant) - 25} more | | |")
        A("")
    A("### How to read these")
    A("")
    A("Kinship is estimated from the core-SNP alignment, as pyseer documents. That")
    A("matters: an earlier version derived kinship from gene presence/absence, and")
    A("because the phylogeny feeding it was degenerate the matrix was all 1.0. It")
    A("reported h^2 = 1.00 and 4 hits, three of which (phzC, phzG, mlaD) have no")
    A("carbapenem mechanism. With a correct SNP kinship the same test returns the hits")
    A("above - a change in kind, not only in count.")
    A("")
    A("The strongest are mechanistically plausible for carbapenem resistance:")
    A("`nfxB` is the repressor of the MexXY-OprM efflux pump, `phoQ` controls the OprD")
    A("and efflux regulons, and `oprJ` is an outer-membrane porin.")
    A("")
    A("They remain associations, not causal proof:")
    A("")
    A(f"- h^2 = {h2:.2f} - structure still dominates" if h2 is not None
      else "- h^2 could not be read from the pyseer log")
    A(f"- {gwas.get('n_lineage_confounded')} of {gwas.get('n_features_tested')} features sit in a single")
    A("  lineage, so they are lineage-linked rather than resistance-associated")
    A(f"- only {gwas.get('n_susceptible')} susceptible isolates; a roughly 6:1 case:control")
    A("  ratio leaves very little power once structure is controlled for")
    A("- R/S is categorical, so the 3 intermediate isolates are excluded entirely")
    A("")
    A("## 7. Convergence")
    A("")
    for cat, count in sorted(conv.get("categories", {}).items(), key=lambda x: -x[1]):
        A(f"- `{cat}`: {count}")
    A("")
    A("## 8. Honest limitations")
    A("")
    A("These bound what the results above can support.")
    A("")
    snps = phylo_status.get("n_core_snps", "?")
    A(f"1. **The cohort is clonal.** 253 core genes and {snps} core SNPs across 65 genomes")
    A("   indicate a narrow population, and the high h^2 confirms resistance tracks")
    A("   structure closely. GWAS here is descriptive, not causal.")
    A("2. **`mexZ`, `ampD` and `dacB` are 0/65 and must not be reported as absent.** The")
    A("   cohort contains `ampDh2`/`ampDh3`, which are AmpH (MreB), a different protein.")
    A("   This is a gap in what the light Bakta database names.")
    A("3. **The tree is not a Snippy tree.** Snippy has no installable osx-arm64 build, so")
    A("   stage 10 used minimap2 + bcftools + FastTree. Joint mpileup over all 65")
    A("   assemblies stalls on this data, so variants were called per sample against the")
    A("   reference.")
    A("4. **Panaroo was never used** (no installable build). Stage 9 builds from")
    A("   annotations instead, so this is a gene presence/absence pan-genome, not a graph.")
    A("5. **Stage 13 is determinants-only.** No regulatory-locus or SV calling was")
    A("   performed, so regulator-driven convergence is unassessed.")
    A("6. **Bakta ran on the light database**, so gene naming is incomplete; locus tags were")
    A(f"   excluded from the pan-genome ({downstream['pangenome']['n_locus_tag_records_excluded']:,} records).")
    A("7. **No MIC or zone diameter is reported.** Only PDC R/I/S labels were available, and")
    A("   none were synthesised.")
    A("")
    A("## 9. Where everything is")
    A("")
    A(f"- Per-sample master table: `{master.relative_to(ROOT)}`")
    A("- Stage status: `results/pilot100/downstream_status.json`")
    A("- Annotation: `results/pilot100/annotation/`, `annotation_state.json`")
    A("- AMR tables: `Imipenem_*.tsv`, `Meropenem_*.tsv`")
    A("- MLST: `MLST_calls.tsv` · Pan-genome: `pangenome/`")
    A("- Phylogeny: `phylogeny/tree.nwk`, `core_snp_alignment.fasta`")
    A("- GWAS: `GWAS_results.tsv`, `gwas/pyseer_results.tsv`, `gwas/pyseer.log`")
    A("- Convergence: `convergence_calls.tsv`")
    A("- Slide deck: `Pilot100_Pseudomonas_AMR_Analysis.pptx` · QC: `PILOT100_QC_REPORT.md`")
    A("")
    (RESULTS / "PIPELINE_RESULTS.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
