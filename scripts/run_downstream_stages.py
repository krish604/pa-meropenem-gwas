#!/usr/bin/env python3
"""Run the pipeline stages that follow annotation, on the prepared cohort.

This driver does not reimplement any stage. It calls the existing stage
entry points in ``papipeline.stages`` in order and supplies the inputs each
one declares:

    stage 5   mlst        -> papipeline.stages.mlst.run
    stage 6   annotation  -> annotation intermediates from Bakta GFF
    stage 9   pangenome   -> papipeline.stages.pangenome.run
    stage 10  phylogeny   -> core-SNP alignment + tree, then
                             papipeline.stages.phylogeny.run
    stage 12  gwas        -> papipeline.stages.gwas.run with PyseerEngine
    stage 13  convergence -> papipeline.stages.convergence.run

Every stage writes its own outputs and reports what it could not do. A
stage that cannot run records the reason and the driver continues, so one
missing tool cannot silently truncate the rest of the run.

Usage:
    python scripts/run_downstream_stages.py --stage all
    python scripts/run_downstream_stages.py --stage phylogeny
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from papipeline.config.loader import load_config  # noqa: E402
from papipeline.io.tsv import write_tsv  # noqa: E402
from papipeline.logging_utils import configure_logging  # noqa: E402
from papipeline.manifest import manifest_from_rows  # noqa: E402
from papipeline.models import PhenotypeCall, Phenotype, RunMode  # noqa: E402

RESULTS = ROOT / "results" / "pilot100"
INTERMEDIATE = RESULTS / "intermediate"
ANNOTATION_SRC = RESULTS / "annotation"
MANIFEST_TSV = RESULTS / "pilot100_manifest.tsv"
STATUS_PATH = RESULTS / "downstream_status.json"

#: The two antibiotics this analysis reports on.
ANTIBIOTICS = ("imipenem", "meropenem")


# ---------------------------------------------------------------------------
# cohort
# ---------------------------------------------------------------------------


def read_rows(path: Path) -> List[Dict[str, str]]:
    """Read a TSV, skipping the leading ``#`` comment block."""
    with path.open(encoding="utf-8") as handle:
        lines = [line for line in handle if not line.startswith("#")]
    return [row for row in csv.DictReader(lines, delimiter="\t")]


def load_cohort() -> Tuple[List[Dict[str, str]], "object"]:
    rows = [r for r in read_rows(MANIFEST_TSV) if r.get("analysis_included", "").upper() == "TRUE"]
    config = load_config(ROOT / "config" / "science.yaml", machine=args.machine)
    manifest = manifest_from_rows(
        ({"sample_id": r["Pilot_ID"], "assembly_path": r["Genome_path"]} for r in rows),
        source="pilot100",
    )
    return rows, (config, manifest)


def phenotype_calls(rows: Sequence[Dict[str, str]]) -> List[PhenotypeCall]:
    """Phenotypes exactly as PDC recorded them.

    ``I`` is kept distinct from ``R`` and ``S``; nothing is recoded and no
    MIC is synthesised from an R/S/I label.
    """
    out: List[PhenotypeCall] = []
    for row in rows:
        label = (row.get("AST_phenotype") or "").strip()
        if not label or label == "MISSING":
            continue
        try:
            pheno = Phenotype(label)
        except ValueError:
            continue
        out.append(
            PhenotypeCall(
                sample_id=row["Pilot_ID"],
                antibiotic="imipenem",
                phenotype=pheno,
                source="PDC",
            )
        )
    return out


def phenotypes_by_sample(rows: Sequence[Dict[str, str]]) -> Dict[str, Dict[str, Phenotype]]:
    """Per-sample R/I/S for every antibiotic, read from the pilot outputs."""
    out: Dict[str, Dict[str, Phenotype]] = {}
    for antibiotic in ANTIBIOTICS:
        path = RESULTS / f"{antibiotic.capitalize()}_phenotype.tsv"
        if not path.exists():
            continue
        for row in read_rows(path):
            label = (row.get("phenotype") or "").strip()
            if label == "MISSING":
                continue
            try:
                out.setdefault(row["Pilot_ID"], {})[antibiotic] = Phenotype(label)
            except ValueError:
                continue
    return out


# ---------------------------------------------------------------------------
# stage 6: annotation intermediates
# ---------------------------------------------------------------------------


def write_annotation_intermediates(rows: Sequence[Dict[str, str]]) -> Dict[str, int]:
    """Convert each Bakta GFF into the pipeline's standardised annotation TSV.

    The GFF is used rather than the annotation TSV because the GFF carries
    both ``gene_id`` and ``gene_name`` in one record; the annotation TSV is
    the same content in tabular form. The main annotation table is preferred
    when both are present.
    """
    from papipeline.stages.annotation import parse_bakta_gff, standardise
    from scripts.run_annotation_sweep import bakta_gff

    out_dir = INTERMEDIATE / "annotation"
    out_dir.mkdir(parents=True, exist_ok=True)
    counts: Dict[str, int] = {}
    for row in rows:
        sample_id = row["Pilot_ID"]
        gff = bakta_gff(ANNOTATION_SRC / sample_id)
        if gff is None:
            counts[sample_id] = 0
            continue
        records = standardise(parse_bakta_gff(gff.read_text()), sample_id)
        write_tsv(out_dir / f"{sample_id}.annotation.tsv", [r.to_row() for r in records],
                  ("sample_id", "contig_id", "gene_id", "gene_name", "product",
                   "gene_type", "start", "end", "strand", "annotation_source"))
        counts[sample_id] = len(records)
    return counts


# ---------------------------------------------------------------------------
# stage 5: MLST
# ---------------------------------------------------------------------------


def _allele_token(raw: str) -> str:
    """Convert one mlst allele cell into the pipeline's ``locus:allele`` form.

    mlst writes ``acsA(1)``, ``acsA(3?)`` for a partial match, and ``acsA(-)``
    for no match. The pipeline's parse_alleles contract is ``locus:allele``,
    so the conversion happens at this boundary rather than by changing the
    parser the rest of the pipeline already depends on.
    """
    text = raw.strip()
    match = re.match(r"^([^(]+)\(([^)]*)\)$", text)
    if not match:
        return text.replace("(", "").replace(")", "")
    locus, allele = match.group(1).strip(), match.group(2).strip()
    partial = allele.endswith("?")
    allele = allele.rstrip("?")
    if allele in ("", "-"):
        return f"{locus}:-"
    return f"{locus}:{allele}" + ("?" if partial else "")


def stage_mlst(rows, config, manifest) -> Dict[str, object]:
    """Call MLST and hand the result to the existing stage module."""
    from papipeline.stages import mlst as mlst_mod

    raw_path = RESULTS / "mlst" / "mlst_raw.tsv"
    if not raw_path.exists():
        return {"status": "skipped", "reason": f"missing {raw_path}"}

    # Map by the manifest's own genome filename. Splitting the assembly name
    # on "_" does not work: "GCA_000710625.1" contains an underscore itself.
    by_filename: Dict[str, str] = {}
    by_assembly: Dict[str, str] = {}
    for r in rows:
        by_filename[Path(r["Genome_path"]).name] = r["Pilot_ID"]
        by_assembly[r["Assembly"]] = r["Pilot_ID"]

    out_rows: List[Dict[str, object]] = []
    unmatched: List[str] = []
    for line in raw_path.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        fname, scheme, st = parts[0], parts[1], parts[2]
        base = Path(fname).name
        sample_id = by_filename.get(base)
        if sample_id is None:
            m = re.match(r"(GC[AF]_\d+\.\d+)", base)
            sample_id = by_assembly.get(m.group(1)) if m else None
        if sample_id is None:
            unmatched.append(base)
            continue
        alleles = parts[3:]
        matched = len([a for a in alleles if "(" in a and not a.endswith("(?)")])
        null_st = st in ("-", "", "ND", "*")
        if null_st or matched == 0:
            status = "no_call"
        elif any(a.endswith("(?)") for a in alleles):
            status = "partial"
        else:
            status = "typed"
        out_rows.append(
            {
                "sample_id": sample_id,
                "mlst_scheme": scheme,
                "ST": "" if null_st else st,
                # mlst emits "acsA(1)"; the pipeline's parse_alleles contract
                # is "locus:allele". Converted here, at the boundary.
                "alleles": ",".join(_allele_token(a) for a in alleles),
                "MLST_status": status,
                "allele_database": "pubmlst",
            }
        )

    out_dir = INTERMEDIATE / "mlst"
    out_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(out_dir / "mlst_results.tsv", out_rows,
              ("sample_id", "mlst_scheme", "ST", "alleles", "MLST_status", "allele_database"))
    calls = mlst_mod.run(config, manifest, RunMode.REAL, INTERMEDIATE)
    typed = sum(1 for c in calls.values() if c.sequence_type)
    write_tsv(RESULTS / "MLST_calls.tsv",
              [{"sample_id": c.sample_id, "ST": c.sequence_type or ".",
                "mlst_status": c.mlst_status, "mlst_scheme": c.mlst_scheme,
                "alleles": ",".join(f"{k}:{v}" for k, v in sorted(c.alleles.items())),
                "allele_database": c.allele_database or "."} for c in
               sorted(calls.values(), key=lambda x: x.sample_id)],
              ("sample_id", "ST", "mlst_status", "mlst_scheme", "alleles", "allele_database"))
    return {"status": "ok", "n_calls": len(calls), "n_typed": typed,
            "n_unmatched_assemblies": len(unmatched)}


# ---------------------------------------------------------------------------
# stage 9: pangenome
# ---------------------------------------------------------------------------


def named_annotations(manifest) -> Dict[str, List]:
    """Standardised annotations restricted to records with a real gene name.

    ``standardise`` falls back to the Locus Tag when Bakta leaves the Gene
    column empty. That fallback is fine for a per-gene report but fatal for
    a pan-genome: Bakta mints run-specific tags (ABLLKF_00001,
    GGCMBK_01800), so every such tag is unique to one genome and each looks
    like a private gene. Including them inflated the pan-genome to 295,854
    "genes" with only 253 core, and would have handed pyseer ~295,000
    meaningless features.

    A record whose gene_name is simply its own Locus Tag therefore carries
    no cross-sample identity and is excluded. The count of excluded records
    is reported rather than silently dropped.
    """
    from papipeline.stages.annotation import load_from_intermediate

    raw = load_from_intermediate(INTERMEDIATE, manifest)
    kept: Dict[str, List] = {}
    dropped = 0
    for sample_id, records in raw.items():
        out = []
        for record in records:
            name = (record.gene_name or "").strip()
            if not name or name == "-" or name == (record.gene_id or "").strip():
                dropped += 1
                continue
            out.append(record)
        kept[sample_id] = out
    return kept, dropped


def stage_pangenome(rows, config, manifest) -> Dict[str, object]:
    from papipeline.stages import pangenome as pg

    try:
        annotations, dropped = named_annotations(manifest)
    except Exception as exc:  # noqa: BLE001
        return {"status": "failed", "reason": f"{type(exc).__name__}: {exc}"}

    pan = pg.build_from_annotations(annotations, manifest.sample_ids)
    paths = pg.write_outputs(pan, RESULTS / "pangenome")
    return {
        "status": "ok",
        "n_samples": pan.n_samples,
        "n_genes": len(pan.genes),
        "n_core": len(pan.core),
        "n_accessory": len(pan.accessory),
        "n_locus_tag_records_excluded": dropped,
        "gene_name_source": "Bakta Gene column only (locus-tag fallback excluded)",
        "outputs": {k: str(v) for k, v in paths.items()},
    }


# ---------------------------------------------------------------------------
# stage 10: phylogeny
# ---------------------------------------------------------------------------


def _run(cmd: Sequence[str], log: Path, stdout_to: Optional[Path] = None) -> Tuple[int, str]:
    """Run a command, appending its diagnostics to ``log``.

    ``stdout_to`` redirects standard output to a file, which minimap2 and
    FastTree need: both write their primary result to stdout, so capturing
    it into the log would bury the SAM/alignment and make it unusable.
    """
    log.parent.mkdir(parents=True, exist_ok=True)
    if stdout_to is not None:
        stdout_to.parent.mkdir(parents=True, exist_ok=True)
        with log.open("a") as err, stdout_to.open("w") as out:
            proc = subprocess.run(cmd, stdout=out, stderr=err)
    else:
        with log.open("a") as handle:
            proc = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT)
    return proc.returncode, log.read_text(errors="replace")[-2000:]


def _tool_candidates(name: str, config=None) -> Tuple[str, ...]:
    """Where to look for ``name``: PATH first, then the machine's own
    configured search directories.

    There is no built-in fallback location. An earlier revision of this file
    fell back to a literal macOS Homebrew path, which on the analysis machine
    does not exist - so a tool that was genuinely missing resolved to nothing
    and was reported as present. If this machine needs a tool outside PATH,
    the overlay says so under `paths.tool_search_dirs`.
    """
    if config is None or config.machine is None:
        return ()
    return config.machine.tool_candidates(name)


def _working_hts_tool(name: str, config=None) -> Optional[str]:
    """Locate a usable samtools/bcftools, skipping the 0.1.19 stubs.

    The pilot100 environment ships samtools and bcftools 0.1.19 (2014
    builds). samtools 0.1.19 cannot read modern SAM - `samtools sort` fails
    with "fail to open file" on its own output and `samtools view` reports
    "this is not a BAM file" - so it silently produces no variants. A
    version check is done rather than trusting PATH order.
    """
    candidates = list(_tool_candidates(name, config))
    for path in candidates:
        if not path or not Path(path).exists():
            continue
        try:
            proc = subprocess.run([path, "--version"], capture_output=True, text=True, timeout=30)
        except Exception:  # noqa: BLE001
            continue
        text = (proc.stdout or "") + (proc.stderr or "")
        match = re.search(r"(\d+)\.(\d+)", text)
        if match and (int(match.group(1)), int(match.group(2))) >= (1, 0):
            return path
    return None


def stage_phylogeny(rows, config, manifest) -> Dict[str, object]:
    """Build a core-SNP phylogeny.

    Snippy cannot be installed on this platform: no osx-arm64 build is
    solvable because perl-bioperl, tabixpp and vcflib are unavailable. The
    substitution is the equivalent route with tools that are present:

        minimap2  assembly-to-reference alignment  (bwa is also available)
        samtools  sort to BAM
        bcftools  mpileup -> call -> per-sample consensus
        FastTree  tree from the SNP sites present in every sample

    Assembly inputs must be aligned before mpileup: bcftools mpileup reads
    BAM/CRAM/VCF, not FASTA. The reference is recorded explicitly, and the
    method is recorded with the result so this tree is never mistaken for a
    Snippy tree.
    """
    minimap2 = next(iter(_tool_candidates("minimap2", config)), None)
    samtools = _working_hts_tool("samtools", config)
    bcftools = _working_hts_tool("bcftools", config)
    fasttree = shutil.which("FastTree") or shutil.which("fasttree")
    for name, path in (("minimap2", minimap2), ("samtools", samtools),
                       ("bcftools", bcftools), ("fasttree", fasttree)):
        if not path or not Path(path).exists():
            return {"status": "failed", "reason": f"{name} not available"}

    out_dir = RESULTS / "phylogeny"
    work = out_dir / "work"
    work.mkdir(parents=True, exist_ok=True)
    log = work / "commands.log"
    if log.exists():
        log.unlink()

    ref = rows[0]["Genome_path"]
    ref_name = Path(ref).stem
    ref_copy = work / "reference.fna"
    if not ref_copy.exists():
        shutil.copy(ref, ref_copy)
    if not (ref_copy.with_suffix(".fna.fai")).exists():
        _run([samtools, "faidx", str(ref_copy)], log)

    # 1. align each assembly to the reference, 2. sort to BAM
    bams = []
    for row in rows:
        sid = row["Pilot_ID"]
        bam = work / f"{sid}.sorted.bam"
        bams.append(str(bam))
        if bam.exists():
            continue
        sam = work / f"{sid}.aligned.sam"
        # --secondary=no: supplementary/secondary alignments of the same
        # region produce spurious pileup evidence and false variant calls.
        code, tail = _run([minimap2, "-x", "asm20", "-a", "--secondary=no",
                           "-t", "4", str(ref_copy), row["Genome_path"]],
                          log, stdout_to=sam)
        if code != 0:
            return {"status": "failed", "reason": f"minimap2 failed for {sid}: {tail[-200:]}"}
        code, tail = _run([samtools, "sort", "-o", str(bam), str(sam)], log)
        if code != 0:
            return {"status": "failed", "reason": f"samtools sort failed for {sid}: {tail[-200:]}"}
        sam.unlink(missing_ok=True)

    # 3. call variants per sample against the reference.
    #
    # Joint mpileup over all 65 assemblies was tried first and abandoned: it
    # stalls at a fixed ~618 MB of output on this data and never finishes.
    # Per-sample calling is fast and, for the purpose here - a core-SNP
    # alignment, one isolate per sample - gives the same sites: a site is
    # core if every sample covers it and at least one carries a different
    # allele. Joint calling matters for low-frequency variants within a
    # pooled population, which is not what this cohort is.
    for row in rows:
        sid = row["Pilot_ID"]
        vcf = work / f"{sid}.vcf.gz"
        # Recomputed from sid. Referring to a `bam` variable that leaked out
        # of the alignment loop above silently called every sample from the
        # last sample's alignment, producing 65 byte-identical VCFs and a
        # core-SNP alignment in which all sequences were identical.
        sample_bam = work / f"{sid}.sorted.bam"
        if not sample_bam.exists():
            return {"status": "failed", "reason": f"alignment missing for {sid}"}
        if not vcf.exists():
            code, tail = _run([bcftools, "mpileup", "-f", str(ref_copy), "-q", "20",
                               "-Q", "20", "-d", "1000", "-a", "FORMAT/DP,AD",
                               "-Ou", "-o", str(work / f"{sid}.pileup.bcf"),
                               str(sample_bam)], log)
            if code != 0:
                return {"status": "failed", "reason": f"mpileup failed for {sid}: {tail[-200:]}"}
            code, tail = _run([bcftools, "call", "-mv", "-Oz", "-o", str(vcf),
                               str(work / f"{sid}.pileup.bcf")], log)
            if code != 0:
                return {"status": "failed", "reason": f"call failed for {sid}: {tail[-200:]}"}
            (work / f"{sid}.pileup.bcf").unlink(missing_ok=True)

    n_raw = sum(_count_variants(work / f"{r['Pilot_ID']}.vcf.gz", bcftools) for r in rows)

    # Guard: per-sample calls must actually differ. Identical VCFs across
    # samples mean the calling step read one input repeatedly, which yields
    # an alignment whose sequences are all equal and a tree with no
    # resolution - a silent, entirely plausible-looking failure.
    fingerprints = {}
    for row in rows:
        sid = row["Pilot_ID"]
        out = subprocess.run(
            [bcftools, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n",
             str(work / f"{sid}.vcf.gz")], capture_output=True, text=True)
        fingerprints[sid] = hashlib.md5(
            "".join(sorted(out.stdout.splitlines())).encode()
        ).hexdigest()
    if len(set(fingerprints.values())) == 1 and len(fingerprints) > 1:
        return {"status": "failed",
                "reason": ("every sample produced an identical variant call set; "
                           "per-sample calling is not reading per-sample input"),
                "n_samples": len(fingerprints)}
    n_distinct = len(set(fingerprints.values()))


    # 4. per-sample consensus from that sample's own calls.
    # bcftools consensus refuses an unindexed VCF ("could not load index"),
    # so each is indexed first.
    cons_dir = work / "consensus"
    cons_dir.mkdir(exist_ok=True)
    for row in rows:
        sid = row["Pilot_ID"]
        out = cons_dir / f"{sid}.fa"
        vcf = work / f"{sid}.vcf.gz"
        if out.exists():
            continue
        _run([bcftools, "index", "-f", str(vcf)], log)
        code, tail = _run([bcftools, "consensus", "-f", str(ref_copy), str(vcf),
                           "-o", str(out)], log)
        if code != 0 or not out.exists():
            return {"status": "failed",
                    "reason": f"bcftools consensus failed for {sid}: {tail[-200:]}"}

    # 5. core SNPs: segregating in the cohort and covered by nearly all of it
    ids = [r["Pilot_ID"] for r in rows]
    candidates, alts = _segregating_sites(work, ids, bcftools)
    if not candidates:
        return {"status": "failed", "reason": "no segregating SNP sites found",
                "n_variants_called": n_raw}
    contig_lengths = _contig_lengths(ref_copy, samtools)
    coverage = _coverage_at_sites(work, ids, candidates, samtools, contig_lengths, log)
    sites = sorted(coverage)
    if not sites:
        return {"status": "failed",
                "reason": (f"none of {len(candidates)} segregating sites is covered by "
                           f">=95% of samples"),
                "n_variants_called": n_raw, "n_candidate_sites": len(candidates)}
    aln = work / "core_snp_alignment.fasta"
    _write_snp_alignment(cons_dir, ids, sites, alts, aln)
    tree = out_dir / "tree.nwk"
    code, tail = _run([fasttree, str(aln)], log, stdout_to=tree)
    if code != 0:
        return {"status": "failed", "reason": f"FastTree failed: {tail[-200:]}"}
    shutil.copy(aln, out_dir / "core_snp_alignment.fasta")
    shutil.copy(aln, out_dir / "core_alignment.fasta")
    return {"status": "ok", "n_variants_called": n_raw,
            "n_distinct_variant_profiles": n_distinct,
            "n_candidate_sites": len(candidates),
            "n_core_snps": len(sites),
            "n_dropped_uncovered_sites": len(candidates) - len(sites),
            "n_samples": len(rows),
            "reference": ref_name,
            "method": "minimap2 asm20 + samtools + bcftools mpileup/call (per sample) + FastTree",
            "method_note": ("snippy has no installable osx-arm64 build; joint mpileup over "
                            "65 assemblies stalls on this data, so variants are called per "
                            "sample against the reference. A site enters the alignment only "
                            "if >=95% of samples have a read covering it."),
            "tree": str(tree), "alignment": str(out_dir / "core_snp_alignment.fasta")}


def _contig_lengths(reference: Path, samtools: str) -> Dict[str, int]:
    out = subprocess.run([samtools, "faidx", str(reference)],
                         capture_output=True, text=True)
    lengths: Dict[str, int] = {}
    for line in out.stdout.splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) >= 2:
            try:
                lengths[parts[0]] = int(parts[1])
            except ValueError:
                continue
    return lengths


def _count_variants(vcf: Path, bcftools: str) -> int:
    """Number of variant records in a VCF."""
    out = subprocess.run([bcftools, "view", "-H", str(vcf)],
                         capture_output=True, text=True)
    return len([line for line in out.stdout.splitlines() if line.strip()])


def _segregating_sites(work: Path, ids: Sequence[str], bcftools: str):
    """Positions where at least one sample carries an alternate allele.

    Returns ``(sites, alts)``. ``sites`` is a list of (contig, position) in
    reference order; ``alts`` maps (contig, position) -> {sample: allele}.

    Sites are keyed by contig *and* position: the reference has 18 contigs,
    so a bare position is ambiguous and keying on position alone silently
    collapses the alignment onto whichever contig was written last. Only
    substitutions are kept, because an indel shifts every downstream
    coordinate.

    Coverage is NOT taken from here. A per-sample VCF holds variant records
    only, so a site absent from a sample's VCF means "not called", not
    "uncovered"; requiring a call in every sample would empty the set.
    Coverage is measured from the alignments by :func:`_coverage_at_sites`.
    """
    alts: Dict[tuple, Dict[str, str]] = {}
    for sid in ids:
        out = subprocess.run(
            [bcftools, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n",
             str(work / f"{sid}.vcf.gz")], capture_output=True, text=True)
        if out.returncode != 0:
            continue
        for line in out.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            ref, alt = parts[2], parts[3].split(",")[0]
            if len(ref) != 1 or len(alt) != 1:
                continue  # indel: skip
            alts.setdefault((parts[0], int(parts[1])), {})[sid] = alt
    return sorted(alts), alts


def _coverage_at_sites(work: Path, ids: Sequence[str], sites: Sequence[tuple],
                       samtools: str, contig_lengths: Dict[str, int],
                       log: Path, min_fraction: float = 0.95):
    """How many samples have an alignment covering each candidate site.

    A site is usable when at least ``min_fraction`` of samples have a
    read covering it, so the alignment is not built from a handful of
    samples with the rest silently filled in with reference bases.
    """
    if not sites:
        return {}
    bed = work / "candidate_sites.bed"
    with bed.open("w") as handle:
        for contig, pos in sites:
            handle.write(f"{contig}\t{pos - 1}\t{pos}\n")
    depth: Dict[tuple, int] = {site: 0 for site in sites}
    for sid in ids:
        bam = work / f"{sid}.sorted.bam"
        if not bam.exists():
            continue
        out = subprocess.run(
            [samtools, "depth", "-a", "-q", "0", "-l", str(bed), str(bam)],
            capture_output=True, text=True)
        if out.returncode != 0:
            LOGGER.warning("samtools depth failed for %s: %s", sid, out.stderr[-200:])
            continue
        for line in out.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            try:
                key = (parts[0], int(parts[1]))
                # samtools may report neighbours of a BED interval; only the
                # candidate sites are of interest.
                if key in depth and int(float(parts[2])) > 0:
                    depth[key] += 1
            except ValueError:
                continue
    needed = max(1, int(round(min_fraction * len(ids))))
    return {site: n for site, n in depth.items() if n >= needed}


def _write_snp_alignment(cons_dir: Path, samples: Sequence[str], sites: Sequence[tuple],
                         alts: Dict[tuple, Dict[str, str]], out: Path) -> None:
    """Emit a FASTA of the core SNP sites, one row per sample.

    ``sites`` are (contig, position) pairs in reference order. A sample's
    base is its own called alternate allele where it has one, otherwise the
    base its consensus carries at that contig and position.

    The consensus is indexed by (contig, position) rather than position
    alone, so the 18 reference contigs cannot overwrite one another.
    """
    sequences: Dict[str, str] = {}
    for sid in samples:
        path = cons_dir / f"{sid}.fa"
        lookup: Dict[tuple, str] = {}
        contig: Optional[str] = None
        coords: List[int] = []
        chunks: List[str] = []

        def flush() -> None:
            if contig is not None:
                for offset, base in enumerate("".join(chunks), start=1):
                    lookup[(contig, offset)] = base

        with path.open() as handle:
            for line in handle:
                if line.startswith(">"):
                    flush()
                    contig = line[1:].strip().split()[0]
                    coords, chunks = [], []
                else:
                    chunks.append(line.strip())
        flush()
        sequences[sid] = "".join(
            alts.get(site, {}).get(sid) or lookup.get(site, "N") for site in sites
        )
    with out.open("w") as handle:
        for sid in samples:
            handle.write(f">{sid}\n{sequences[sid]}\n")


# ---------------------------------------------------------------------------
# stage 12: GWAS
# ---------------------------------------------------------------------------


def stage_gwas(rows, config, manifest) -> Dict[str, object]:
    """Stage 12: write the feature table and let the stage run pyseer.

    Feature columns use the declared ``<type>__<label>`` naming so the
    stage can recover the feature type. The phenotype is passed as
    recorded: R is positive, S negative, and I/SDD/ND are excluded by
    config.gwas.outcome rather than folded into either group.
    """
    from papipeline.stages import gwas as gwas_mod
    from papipeline.stages import pangenome as pg

    pheno = phenotypes_by_sample(rows)
    calls: List[PhenotypeCall] = []
    for sid, per in sorted(pheno.items()):
        label = per.get(config.antibiotics[0])
        if label is None:
            continue
        calls.append(PhenotypeCall(sample_id=sid, antibiotic=config.antibiotics[0],
                                   phenotype=label, source="PDC"))
    if not calls:
        return {"status": "failed", "reason": "no phenotypes available"}

    annotations, _dropped = named_annotations(manifest)
    pan = pg.build_from_annotations(annotations, manifest.sample_ids)
    amr = _amr_calls(rows)

    presence: Dict[str, set] = {g: set(pan.presence.get(g, ())) for g in pan.genes}
    amr_genes = sorted({g for calls_ in amr.values() for g in calls_})
    columns: List[Tuple[str, str, set]] = []
    for gene in pan.genes:
        columns.append(("gene__presence_absence", gene, presence[gene]))
    for gene in amr_genes:
        if gene in pan.genes:
            continue
        carriers = {sid for sid, genes in amr.items() if gene in genes}
        columns.append(("gene_presence_absence", gene, carriers))

    out_dir = INTERMEDIATE / "gwas"
    out_dir.mkdir(parents=True, exist_ok=True)
    feature_rows = []
    for sample in manifest.sample_ids:
        row: Dict[str, object] = {"sample_id": sample}
        for ftype, label, carriers in columns:
            row[f"{ftype}__{label}"] = int(sample in carriers)
        feature_rows.append(row)
    header = ["sample_id"] + [f"{t}__{lab}" for t, lab, _ in columns]
    write_tsv(out_dir / "gwas_features.tsv", feature_rows, tuple(header))

    lineages = _lineages(rows)
    # Prefer a core-SNP alignment for kinship: pyseer documents SNP-based
    # kinship, and it resolves relatedness far better than the gene
    # presence/absence matrix. Falls back to the features if absent.
    snp_aln = RESULTS / "phylogeny" / "core_snp_alignment.fasta"
    engine = gwas_mod.PyseerEngine(
        shutil.which("pyseer") or "pyseer", RESULTS / "gwas",
        snp_alignment=snp_aln if snp_aln.exists() else None,
    )
    results, gwas_input = gwas_mod.run(config, manifest, RunMode.REAL, INTERMEDIATE,
                                        calls, lineages=lineages, engine=engine)
    hits = gwas_mod.significant(results, config)
    confounded = gwas_mod.flag_lineage_linked(results, gwas_input, config)
    write_tsv(RESULTS / "GWAS_results.tsv",
              [r.to_row() for r in sorted(results, key=lambda x: (x.adjusted_p_value or 1.0))],
              tuple(results[0].to_row().keys()) if results else ("feature",))
    return {
        "status": "ok",
        "antibiotic": config.antibiotics[0],
        "n_samples": gwas_input.n_samples,
        "n_resistant": gwas_input.n_positive,
        "n_susceptible": gwas_input.n_negative,
        "n_excluded_phenotype": gwas_input.n_excluded,
        "excluded_samples": list(gwas_input.excluded_samples),
        "n_features_written": len(columns),
        "n_features_tested": len(results),
        "n_significant_adjp_le_0.05": len(hits),
        "n_lineage_confounded": len(confounded),
        "engine": engine.name,
        "model": "pyseer:mixed",
        "kinship_source": engine.kinship_source,
    }


def _lineages(rows: Sequence[Dict[str, str]]) -> Dict[str, str]:
    path = RESULTS / "MLST_calls.tsv"
    if not path.exists():
        return {r["Pilot_ID"]: "unknown" for r in rows}
    return {r["sample_id"]: (r["ST"] if r["ST"] not in (".", "") else "unknown")
            for r in read_rows(path)}


# ---------------------------------------------------------------------------
# stage 13: convergence
# ---------------------------------------------------------------------------


def stage_convergence(rows, config, manifest) -> Dict[str, object]:
    from papipeline.stages import convergence as conv

    pheno = phenotypes_by_sample(rows)
    lineages = _lineages(rows)
    amr_calls = _amr_calls(rows)
    regulator_variants = _regulator_variants(rows)
    calls = conv.run(config, manifest, RunMode.REAL, amr_calls, regulator_variants,
                     lineages, n_samples=len(rows))
    categories: Dict[str, int] = {}
    for call in calls:
        categories[call.convergence_category] = categories.get(call.convergence_category, 0) + 1
    write_tsv(RESULTS / "convergence_calls.tsv", [c.to_row() for c in calls],
              tuple(calls[0].to_row().keys()) if calls else ("sample_id",))
    return {"status": "ok", "n_calls": len(calls),
            "categories": categories,
            "regulator_variants_supplied": 0,
            "limitation": ("convergence is scored on acquired AMR determinants only; "
                           "no regulatory-locus variant calling and no SV calling were "
                           "performed, so regulator-driven convergence is not assessed")}


def _amr_calls(rows: Sequence[Dict[str, str]]) -> Dict[str, List[str]]:
    """AMR gene names per sample, from the AMRFinderPlus reports.

    Uses the pilot's v4 report parser rather than the main pipeline's
    ``load_amr_table``: the files on disk are raw AMRFinderPlus output, and
    the main loader expects the pipeline's own normalised schema. Reading
    them with the wrong loader would either fail or silently mis-map every
    determinant to ``unmapped``.
    """
    from papipeline.pilot.amr_detect import parse_amrfinder_report

    db_dir = Path(_amr_database_dir())
    version = ""
    try:
        from papipeline.pilot.amr_detect import database_version
        version = database_version(db_dir)
    except Exception:  # noqa: BLE001
        version = "unknown"

    out: Dict[str, List[str]] = {}
    for row in rows:
        path = RESULTS / "amrfinder" / f"{row['Pilot_ID']}.amrfinder.tsv"
        if not path.exists():
            out[row["Pilot_ID"]] = []
            continue
        determinants = parse_amrfinder_report(
            path, row["Pilot_ID"], "imipenem", "amrfinderplus", version
        )
        out[row["Pilot_ID"]] = sorted({d.gene for d in determinants if d.gene})
    return out


def _amr_database_dir() -> str:
    from papipeline.pilot.orchestrator import _locate_database
    return str(_locate_database())


def _regulator_variants(rows: Sequence[Dict[str, str]]) -> Dict[str, List[str]]:
    """Regulator-locus variant names per sample: none, and that is the truth.

    Stage 13 merges regulator variants with acquired determinants when
    scoring convergence, so anything passed here is treated as a determinant
    name. The PDC ``AMR_genotypes`` field holds raw metadata encodings such
    as ``blaPDC-3=COMPLETE`` - that is a genotype call for a locus, not a
    variant observed in a regulator, and feeding it in produces
    "determinants" like ``aac(3)-Ic=COMPLETE``.

    No variant calling was performed at nfxB, mexR, ampR or any other
    regulatory locus in this run, so this returns nothing rather than
    something misleading. The gap is recorded in the stage status.
    """
    return {row["Pilot_ID"]: [] for row in rows}


# ---------------------------------------------------------------------------


STAGES = {
    "mlst": stage_mlst,
    "pangenome": stage_pangenome,
    "phylogeny": stage_phylogeny,
    "gwas": stage_gwas,
    "convergence": stage_convergence,
}
ORDER = ["mlst", "pangenome", "phylogeny", "gwas", "convergence"]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", default="all",
                        help="all, or one of: " + ", ".join(ORDER))
    parser.add_argument("--log-file", type=Path, default=None)
    args = parser.parse_args(argv)
    configure_logging("INFO", logfile=args.log_file)

    rows, (config, manifest) = load_cohort()
    print(f"cohort: {len(rows)} included assemblies", flush=True)

    status: Dict[str, object] = {}
    if STATUS_PATH.exists():
        status = json.loads(STATUS_PATH.read_text())

    # Every downstream stage reads standardised annotations, so materialise
    # them once from the Bakta GFFs.
    print("\n=== stage: annotation intermediates ===", flush=True)
    counts = write_annotation_intermediates(rows)
    total = sum(counts.values())
    empty = [k for k, v in counts.items() if v == 0]
    print(json.dumps({"samples": len(counts), "total_records": total,
                      "samples_with_no_records": empty}, indent=2), flush=True)
    status["annotation"] = {"status": "ok", "samples": len(counts),
                            "total_records": total, "empty_samples": empty}
    STATUS_PATH.write_text(json.dumps(status, indent=2))

    todo = ORDER if args.stage == "all" else [args.stage]
    for name in todo:
        print(f"\n=== stage: {name} ===", flush=True)
        started = time.time()
        try:
            result = STAGES[name](rows, config, manifest)
        except Exception as exc:  # noqa: BLE001
            result = {"status": "failed", "reason": f"{type(exc).__name__}: {exc}"}
        result["seconds"] = round(time.time() - started, 1)
        status[name] = result
        STATUS_PATH.write_text(json.dumps(status, indent=2))
        print(json.dumps(result, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
