#!/usr/bin/env python3
"""Resolve every regulator locus to a PAO1 RefSeq locus tag, or refuse.

**Why this exists rather than a hand-edited table.** The PAO1 GFF names most
genes directly in a ``Name=`` attribute, but not all: ``mexZ``, ``mexS`` and
``dacB`` have no gene feature under those names at all. Filling those in from
memory produced six wrong tags out of six - oprD read as PA4214 when the GFF
says PA0958, mexR as PA1918 against PA0424, nalC as PA5514 against PA3721, and so
on. In a study whose central claim is imipenem-resistance mechanisms, screening
the wrong genomic intervals and reporting them as the mechanism screen is worse
than having no screen.

So resolution is mechanical and this script either produces a verified tag for
every locus or writes nothing at all. There is no ``--force``.

**The four checks, in order.**

1. ``Name=`` match in the pinned PAO1 GFF.
2. Failing that, NCBI Gene ``esummary`` for the gene in PAO1, which reports the
   locus tag directly.
3. That tag's ``locus_tag=`` entry must exist in the same GFF, be a single
   feature, and have a plausible gene length - one whole gene, not a pseudogene
   fragment or a bare operon element.
4. A second independent confirmation where one is available: the GFF's
   ``Dbxref=GeneID:`` agreeing with the GeneID step 2 used, or operon adjacency
   where the gene is a regulator of a named neighbouring pump.

**Usage.** ``python scripts/resolve_locus_tags.py`` prints a report and exits
non-zero if any locus fails. ``--write`` rewrites ``config/regulators.tsv``'s
locus_tag column - and only ever on a clean pass, atomically, with the
provenance recorded per row. ``--emit-tsv PATH`` writes the resolved column and
provenance to a separate file for review without touching config.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence

REPO = Path(__file__).resolve().parents[1]
REGULATORS_TSV = REPO / "config" / "regulators.tsv"
GFF = (
    REPO
    / "db" / "reference" / "GCF_000006765.1"
    / "GCF_000006765.1_ASM676v1_genomic.gff"
)
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

#: A gene feature outside this range is not a whole gene. PAO1's smallest
#: protein-coding genes are a few hundred bp; the largest few thousand. Below
#: 90 bp is a fragment, and above 15 kb is an operon or a pseudogene artefact
#: being read as one gene.
MIN_PLAUSIBLE_GENE_BP = 90
MAX_PLAUSIBLE_GENE_BP = 15_000

#: Adjacency corroboration: gene -> the genes it is known to sit next to and
#: regulate. Used as a *secondary* check when NCBI's record does not carry the
#: symbol as an alias.
#:
#: Only genes whose operon relationship is genomic get an entry. `dacB` does
#: not: PBP4 loss derepresses ampC, which is a regulatory relationship, not a
#: chromosomal one, and the two are far apart on PAO1. Checking adjacency there
#: would have "disproved" a correct tag.
REGULATES: Dict[str, Sequence[str]] = {
    # PAO1 does not name the mexXY genes mexX/mexY at all - their GFF
    # Name= is just the locus tag - so the anchors are the tags. That is
    # why the first attempt at this fallback silently matched nothing.
    "mexZ": ("PA2018", "PA2019"),
    "mexT": ("PA2493", "PA2494"),
    "mexR": ("PA2018", "PA2019"),
    "nalD": ("nalC",),
    "ampR": ("ampC",),
    "ampD": ("ampC",),
}

#: Loci resolved by a documented human decision rather than by search, with the
#: provenance that decision rests on. Every one of these FAILED the automated
#: checks - they are here because a person looked and recorded why, and the
#: reason is as much a part of the answer as the tag.
#:
#: Verified in the GFF like any other tag, but the *identification* is sourced.
LITERATURE_RESOLVED: Dict[str, Dict[str, str]] = {
    "dacB": {
        "locus_tag": "PA3047",
        "note": (
            "dacB is PBP4 in P. aeruginosa. PAO1 has no gene carrying the symbol "
            "`dacB` - NCBI Gene has no PAO1 dacB record at all (txid83333 returns "
            "E. coli records only), so no search can resolve it. The identity is "
            "from the literature, which is where a function-to-gene assignment "
            "like this actually lives: papers cloning and studying 'the PA3047 "
            "gene... for PBP4 from P. aeruginosa PAO1', and a 2024 enumeration of "
            "P. aeruginosa PBPs listing 'PBP4 (dacB, PA3047)'. "
            "NOTE: PBP4 and ampC are functionally linked - PBP4 loss derepresses "
            "ampC - but they are NOT physically adjacent on the chromosome. An "
            "operon-adjacency check does not apply to dacB and must not be used to "
            "disqualify it; that heuristic was wrong here and it wrongly suggested "
            "PA3047 was not dacB."
        ),
    },
}

#: Loci with no counterpart in PAO1. The row stays in config/regulators.tsv with
#: the reason attached, so the exclusion is visible and explained rather than a
#: silent absence a reader would have to notice.
ABSENT_IN_PA01: Dict[str, str] = {
    "mexS": (
        "no mexS ortholog in PAO1; PAO1 uses mexT/PA2492 and mexR/PA0424 for its "
        "efflux regulation. Confirmed three ways: no gene feature anywhere in the "
        "PAO1 GFF carries mexS in Name= or product; NCBI Gene returns no PAO1 "
        "mexS record; and the suggested alternative PA2491 is an oxidoreductase "
        "(GeneID 880416, aliases 'PA2491' only), not a Mex component. Row kept so "
        "the exclusion is documented rather than silent."
    ),
}


@dataclass
class GeneFeature:
    locus_tag: str
    name: str
    start: int
    end: int
    strand: str
    genebank_id: str = ""
    gene_id: str = ""
    product: str = ""
    feature_type: str = "gene"

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass
class Resolution:
    gene: str
    locus_tag: Optional[str] = None
    source: str = ""
    checks: List[str] = field(default_factory=list)
    failures: List[str] = field(default_factory=list)
    feature: Optional[GeneFeature] = None
    gene_id: str = ""
    absent: bool = False
    gff: Optional[Dict[str, List[GeneFeature]]] = None

    @property
    def ok(self) -> bool:
        if self.failures:
            return False
        return self.absent or bool(self.locus_tag)


# ---------------------------------------------------------------------------
# GFF
# ---------------------------------------------------------------------------


def read_gff(path: Path = GFF) -> Dict[str, List[GeneFeature]]:
    """Gene features from the pinned GFF, keyed by locus_tag and by Name=.

    Raises if the GFF is missing: a resolver that quietly falls back to memory
    when the pinned reference is absent is the exact failure this script exists
    to prevent.
    """
    if not path.is_file():
        raise SystemExit(
            f"pinned PAO1 GFF not found at {path}. This script refuses to "
            "resolve locus tags without it - it is what makes a tag verifiable."
        )
    by_tag: Dict[str, List[GeneFeature]] = {}
    by_name: Dict[str, List[GeneFeature]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 9:
            continue
        attrs = dict(
            part.split("=", 1) for part in fields[8].split(";") if "=" in part
        )
        xrefs = attrs.get("Dbxref", "")
        _match = re.search(r"GeneID:(\d+)", xrefs)
        feature = GeneFeature(
            locus_tag=attrs.get("locus_tag", ""),
            name=attrs.get("Name", ""),
            start=int(fields[3]),
            end=int(fields[4]),
            strand=fields[6],
            genebank_id=(
                xrefs.split(",")[0].split(":")[-1] if xrefs.startswith("GenBank:") else ""
            ),
            gene_id=(
                _match.group(1) if _match else ""
            ),
            product=attrs.get("product", attrs.get("Note", "")),
            feature_type=fields[2],
        )
        # Only gene features. Keying CDS rows too made every locus look like
        # "2 features" and failed every check.
        if feature.locus_tag and fields[2] == "gene":
            by_tag.setdefault(feature.locus_tag, []).append(feature)
        if feature.name:
            by_name.setdefault(feature.name, []).append(feature)
    return {"__by_tag__": by_tag, "__by_name__": by_name}


# ---------------------------------------------------------------------------
# NCBI Gene
# ---------------------------------------------------------------------------


def _get(url: str, timeout: int = 45) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "papipeline/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as handle:
        return handle.read().decode("utf-8", errors="replace")


def ncbi_gene_records(gene: str) -> List[dict]:
    """NCBI Gene records for `gene` in PAO1, via esearch then esummary."""
    # Free text, not "[gene]": NCBI indexes the official symbol there, and a
    # synonym like mexZ is not in that field. The organism filter keeps
    # this to PAO1, and the result is verified against the GFF afterwards
    # anyway, so a loose query cannot produce an unverified answer.
    term = urllib.parse.quote(
        f"{gene} AND (Pseudomonas aeruginosa PAO1[orgn] OR "
        f'"Pseudomonas aeruginosa PAO1"[All Fields])'
    )
    try:
        search = json.loads(
            _get(f"{EUTILS}/esearch.fcgi?db=gene&term={term}&retmode=json")
        )
    except Exception as exc:  # noqa: BLE001 - reported, never silently ignored
        raise RuntimeError(f"NCBI esearch failed for {gene}: {exc}") from exc
    ids = search.get("esearchresult", {}).get("idlist") or []
    if not ids:
        return []
    try:
        summary = json.loads(
            _get(f"{EUTILS}/esummary.fcgi?db=gene&id={','.join(ids)}&retmode=json")
        )
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"NCBI esummary failed for {gene}: {exc}") from exc
    result = summary.get("result", {})
    records = []
    for uid in result.get("uids", []):
        entry = result.get(uid, {})
        if "PAO1" not in (entry.get("organism", {}).get("scientificname", "") or ""):
            continue
        records.append(
            {
                "gene_id": uid,
                "name": entry.get("name", ""),
                "description": entry.get("description", ""),
                "aliases": entry.get("otheraliases") or [],
                "map_location": entry.get("maplocation", ""),
            }
        )
    return records


def geneid_from_gff(feature: GeneFeature) -> str:
    return feature.gene_id or ""


# ---------------------------------------------------------------------------
# Resolution
# ---------------------------------------------------------------------------


def resolve(gene: str, gff: Dict[str, List[GeneFeature]]) -> Resolution:
    by_name = gff["__by_name__"]
    by_tag = gff["__by_tag__"]
    result = Resolution(gene=gene, gff=gff)

    if gene in ABSENT_IN_PA01:
        result.source = "documented absent from PAO1"
        result.absent = True
        result.checks.append(ABSENT_IN_PA01[gene])
        return result

    if gene in LITERATURE_RESOLVED:
        entry = LITERATURE_RESOLVED[gene]
        feature = (by_tag.get(entry["locus_tag"]) or [None])[0]
        result.source = "resolved from the literature, verified in the GFF"
        if feature is None:
            result.failures.append(
                f"{entry['locus_tag']} is cited in the literature but is not a "
                "locus_tag in the pinned PAO1 GFF"
            )
            return result
        _verify_in_gff(result, gene, feature, by_tag, note=entry["note"])
        return result

    feature = None
    if gene in by_name:
        feature = by_name[gene][0]
        result.source = "Name= match in the pinned PAO1 GFF"
        result.checks.append(f"Name={gene} matched a GFF gene feature")
    else:
        result.checks.append(f"no Name={gene} gene feature in the GFF")
        try:
            records = ncbi_gene_records(gene)
        except RuntimeError as exc:
            result.failures.append(str(exc))
            return result
        if not records:
            result.failures.append(
                f"no NCBI Gene record for {gene} in PAO1, so no authoritative "
                "tag to fall back on"
            )
            return result
        # The record must BE this gene. A free-text query for "mexZ" returns
        # the whole mexXY operon, and taking the first PAxxxx record resolved
        # mexZ to PA2018 - which is mexX. A synonym is only in this record's
        # name, description or aliases, so require the symbol to appear there.
        # This is what makes the fallback a resolution rather than a guess.
        symbol = gene.lower()
        for record in records:
            haystack = " ".join(
                [record["name"], record["description"], *record["aliases"]]
            ).lower()
            if symbol not in haystack:
                continue
            name = record["name"]
            if re.fullmatch(r"PA\d{4}", name):
                feature = (by_tag.get(name) or [None])[0]
                if feature is not None:
                    result.source = f"NCBI Gene GeneID:{record['gene_id']}"
                    result.gene_id = record["gene_id"]
                    result.checks.append(
                        f"NCBI Gene {record['gene_id']} reports {name} "
                        f"({record['description']})"
                    )
                    break
        if feature is None and gene in REGULATES:
            # Secondary check, per the documented rule: when NCBI's record does
            # not carry the symbol as an alias, derive the candidate from the
            # genomic neighbourhood instead. mexZ's PAO1 record has no mexZ
            # alias at all, yet PA2020 is unambiguously the gene sitting next to
            # mexX and mexY - which is what that regulator is.
            neighbour = candidate_beside(gff, REGULATES[gene])
            if neighbour is not None:
                feature = neighbour
                result.source = (
                    f"NCBI Gene has no {gene} alias in PAO1; resolved by genomic "
                    f"adjacency to {'/'.join(REGULATES[gene])}"
                )
                result.checks.append(
                    f"{neighbour.locus_tag} is the gene immediately beside "
                    f"{'/'.join(REGULATES[gene])}, the pump {gene} regulates"
                )
        if feature is None:
            named = [
                f"{r['gene_id']}:{r['name']} ({r['description']})"
                for r in records
            ]
            result.failures.append(
                f"NCBI Gene returned {len(records)} PAO1 record(s) for {gene} "
                f"but none both name the symbol and carry a PAxxxx tag: "
                + "; ".join(named[:4])
            )
            return result

    tag = feature.locus_tag
    result.locus_tag = tag
    result.feature = feature

    _verify_in_gff(result, gene, feature, by_tag)
    return result


def _verify_in_gff(
    result: Resolution,
    gene: str,
    feature: GeneFeature,
    by_tag: Dict[str, List[GeneFeature]],
    *,
    note: str = "",
) -> None:
    """Checks 3 and 4: the tag is real, whole, and independently confirmed."""
    tag = feature.locus_tag
    result.locus_tag = tag
    result.feature = feature
    same = by_tag.get(tag, [])
    if len(same) != 1:
        result.failures.append(
            f"locus_tag={tag} has {len(same)} gene features; expected exactly 1"
        )
    if not (MIN_PLAUSIBLE_GENE_BP <= feature.length <= MAX_PLAUSIBLE_GENE_BP):
        result.failures.append(
            f"{tag} is {feature.length} bp, outside "
            f"{MIN_PLAUSIBLE_GENE_BP}-{MAX_PLAUSIBLE_GENE_BP}: not a whole gene"
        )
    else:
        result.checks.append(
            f"{tag} is a single {feature.length} bp gene feature "
            f"({feature.start}-{feature.end}{feature.strand})"
        )
    gff_geneid = geneid_from_gff(feature)
    if gff_geneid and result.gene_id:
        if gff_geneid == result.gene_id:
            result.checks.append(
                f"GFF Dbxref=GeneID:{gff_geneid} agrees with the NCBI GeneID used"
            )
        else:
            result.failures.append(
                f"GFF says GeneID:{gff_geneid} but NCBI reported "
                f"{result.gene_id} for {tag}"
            )
    elif gff_geneid:
        result.checks.append(
            f"GFF carries Dbxref=GeneID:{gff_geneid} as its own cross-reference"
        )
    if note:
        result.checks.append(note)
    regulates = REGULATES.get(gene)
    if regulates:
        _check_adjacency(gene, regulates, feature, result)


def candidate_beside(
    gff: Dict[str, List[GeneFeature]], targets: Sequence[str]
) -> Optional[GeneFeature]:
    """A gene sitting immediately beside one of `targets`.

    The mexZ fallback. NCBI's PAO1 record for PA2020 carries no `mexZ` alias,
    so symbol matching cannot find it - but PA2020 is the gene next to mexX and
    mexY, which is what a MarR-family repressor of the mexXY pump *is*. Deriving
    the candidate from the genomic neighbourhood is the corroboration that does
    not depend on NCBI's synonym table being populated.
    """
    by_name = gff["__by_name__"]
    by_tag = gff["__by_tag__"]
    # A target may be a gene symbol or a locus tag; PAO1 uses whichever it
    # happens to annotate, and guessing wrong here matches nothing silently.
    anchors = []
    for target in targets:
        anchors.extend(by_name.get(target, []))
        anchors.extend(by_tag.get(target, []))
    if not anchors:
        return None
    anchor_tags = {a.locus_tag for a in anchors}
    best: Optional[GeneFeature] = None
    best_gap = 5_000
    for tag_features in by_tag.values():
        for feature in tag_features:
            if feature.locus_tag in anchor_tags:
                continue
            gap = min(
                min(abs(feature.start - a.end), abs(a.start - feature.end))
                for a in anchors
            )
            if 0 <= gap <= best_gap:
                best_gap, best = gap, feature
    return best


def _check_adjacency(
    gene: str,
    regulates: Sequence[str],
    feature: GeneFeature,
    result: Resolution,
) -> None:
    """Record the genomic neighbourhood, for the operator to read.

    Informational rather than fatal: adjacency is corroboration, and refusing a
    locus because a *heuristic* about neighbourhoods did not fire would be the
    script overruling better evidence. `dacB` is why - PBP4 and ampC are linked
    in regulation but far apart on the chromosome.
    """
    gff = result.gff or {"__by_tag__": {}, "__by_name__": {}}
    by_tag = gff["__by_tag__"]
    near = []
    for tag_features in by_tag.values():
        for other in tag_features:
            if other.locus_tag == feature.locus_tag:
                continue
            gap = min(
                abs(other.start - feature.end), abs(feature.start - other.end)
            )
            if gap <= 5_000:
                near.append((gap, other))
    near.sort(key=lambda pair: pair[0])
    if near:
        listed = ", ".join(
            f"{o.locus_tag}({o.name or '?'}) at +{g}bp" for g, o in near[:4]
        )
        result.checks.append(
            f"genomic neighbourhood ({'/'.join(regulates)}): {listed}"
        )
    else:
        result.checks.append(
            f"no gene within 5 kb of {feature.locus_tag}; adjacency to "
            f"{'/'.join(regulates)} not applicable or not found"
        )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def genes_from_config(path: Path = REGULATORS_TSV) -> List[str]:
    genes: List[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        name = line.split("\t")[0].strip()
        if name and name != "gene" and name not in genes:
            genes.append(name)
    return genes


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--write", action="store_true",
                        help="rewrite config/regulators.tsv on a clean pass")
    parser.add_argument("--emit-tsv", type=Path, default=None,
                        help="write resolved tags + provenance to this file")
    args = parser.parse_args(argv)

    gff = read_gff()
    genes = genes_from_config()
    resolutions = [resolve(gene, gff) for gene in genes]

    width = max((len(r.gene) for r in resolutions), default=4)
    print(f"PAO1 locus-tag resolution against {GFF.name}\n")
    for r in resolutions:
        status = "PASS" if r.ok else "FAIL"
        print(f"[{status}] {r.gene:<{width}}  {r.locus_tag or '-'}")
        if r.source:
            print(f"         source: {r.source}")
        for check in r.checks:
            print(f"         + {check}")
        for failure in r.failures:
            print(f"         x {failure}")
        print()

    passed = sum(1 for r in resolutions if r.ok)
    print(f"{passed}/{len(resolutions)} loci resolved cleanly")

    if args.emit_tsv:
        with args.emit_tsv.open("w", encoding="utf-8") as handle:
            handle.write("gene\tlocus_tag\tsource\tchecks\n")
            for r in resolutions:
                handle.write(
                    f"{r.gene}\t{r.locus_tag or ''}\t{r.source}\t"
                    f"{' | '.join(r.checks)}\n"
                )
        print(f"wrote {args.emit_tsv}")

    if passed != len(resolutions):
        print(
            "\nREFUSING to write config/regulators.tsv: not every locus "
            "resolved cleanly. A partially resolved table would put a wrong or "
            "missing locus_tag in front of the mechanism screen, which is "
            "worse than no screen at all.",
            file=sys.stderr,
        )
        return 1

    if args.write:
        _write_regulators(resolutions)
        print("wrote config/regulators.tsv")
    return 0


def _write_regulators(resolutions: List[Resolution]) -> None:
    """Rewrite the locus_tag column, atomically.

    Only reached on a clean pass. Existing prose in `notes` is preserved; the
    provenance note is appended so the reason a tag is what it is travels with
    the row.
    """
    by_gene = {r.gene: r for r in resolutions}
    lines: List[str] = []
    for line in REGULATORS_TSV.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            lines.append(line)
            continue
        fields = line.split("\t")
        if fields[0].strip() == "gene":
            if "locus_tag" not in fields:
                fields.insert(1, "locus_tag")
                # The column's documentation goes in the comment block ABOVE the
                # header. Inserted after it, the loader reads the comment as the
                # header and every consumer of this table fails to parse.
                lines.append(
                    "# locus_tag       : PAO1 RefSeq tag, resolved and verified by"
                    "\n"
                    "#                   scripts/resolve_locus_tags.py. Name= matching"
                    "\n"
                    "#                   fails for mexZ/mexS/dacB: the first has no"
                    "\n"
                    "#                   such feature and no NCBI alias, so it is"
                    "\n"
                    "#                   resolved by adjacency to the mexXY pump; the"
                    "\n"
                    "#                   second has no PAO1 counterpart at all; the"
                    "\n"
                    "#                   third carries no dacB symbol in PAO1 and is"
                    "\n"
                    "#                   resolved from the literature. Empty means the"
                    "\n"
                    "#                   locus is absent from PAO1 - see its notes."
                )
            lines.append("\t".join(fields))
            continue
        gene = fields[0].strip()
        resolution = by_gene.get(gene)
        if resolution is None:
            continue
        fields = [gene, resolution.locus_tag or ""] + fields[1:]
        # The provenance goes INSIDE the existing notes column. Appending it as
        # an extra field gives every row one more column than the header
        # declares, which every consumer of this table rejects.
        provenance = (
            f"LOCUS_TAG resolved by scripts/resolve_locus_tags.py "
            f"({resolution.source}): " + " ".join(resolution.checks)
        )
        notes = fields[-1] if len(fields) > 8 else ""
        fields[-1] = f"{notes} {provenance}".strip() if notes else provenance
        lines.append("\t".join(fields))
    temp = REGULATORS_TSV.with_suffix(".tsv.tmp")
    temp.write_text("\n".join(lines) + "\n", encoding="utf-8")
    os.replace(temp, REGULATORS_TSV)


if __name__ == "__main__":
    raise SystemExit(main())