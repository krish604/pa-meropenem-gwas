"""Stage 6a: the cohort-wide merge of per-isolate variant calls.

Per-isolate calls conflate two different things. A site every isolate differs
from PAO1 at is a **species-wide fixed difference**, not a polymorphism, and no
per-sample processing can tell those apart - only a multi-isolate comparison
can. This stage is that comparison.

Output is long/tidy: one row per `(chrom, pos, ref, alt)`, with per-site carrier
and cohort counts. That shape is settled by reading pyseer, which is wide only
on its `--pres` path, where samples come from the *header* row. A wide table
would carry one column per isolate - 835 columns - and every downstream consumer
would inherit it. The pivot to an Rtab matrix is a separate concern and has no
consumer yet, so it has no home.

The filter is on **minor-allele frequency**, `1/N <= maf <= max_maf`, not on a
raw carrier count. That is not a preference; the measurement is in
`docs/design/cohort-variant-merge.md`: across four real isolates, 68.3% of sites
are carried by exactly one, so "carried by at least two isolates" retains 31.7%
of the genome - the opposite of conservative - and a fixed count does not scale,
being 50% of a 4-isolate cohort and 0.24% of an 835-isolate one.

`max_maf` is not defaulted. It is unset because a minor-allele frequency needs a
real minor allele and four isolates cannot show a rare one; guessing it is the
PROVISIONAL-columns mistake one level up. Callers that need it must pass it.
"""

from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from ..errors import DataContractError
from ..logging_utils import get_logger

LOGGER = get_logger("stages.cohort_variants")

#: The merge output contract. Long/tidy; see the module docstring.
MERGE_COLUMNS: Tuple[str, ...] = ("chrom", "pos", "ref", "alt", "ac", "an", "af")

#: `(chrom, pos, ref, alt)` for one isolate's call.
Call = Tuple[str, str, str, str]


def merge_calls(
    calls_by_isolate: Mapping[str, Sequence[Call]],
    *,
    max_maf: Optional[float] = None,
    require_max_maf: bool = False,
    min_dissenters: Optional[int] = 2,
) -> List[Dict[str, str]]:
    """Merge per-isolate calls into cohort-polymorphic sites.

    A site is retained when its minor-allele frequency is at least ``1/N``, which
    makes "carried by at least two isolates" implicit rather than a separate
    rule: a lone carrier is a 1/N allele by definition.

    Args:
        calls_by_isolate: Each isolate's calls, keyed by isolate.
        max_maf: Optional explicit frequency ceiling. Usually left unset: the
            symmetric ``min_dissenters`` bound below is the operative one.
        min_dissenters: Isolates that must *not* carry a site for it to be
            retained. Defaults to 2 - the symmetric partner of the two-carrier
            floor, since a lone dissenter is as suspect as a lone carrier. Pass
            None to disable the ceiling.
        require_max_maf: Refuse when neither bound was supplied.

    Returns:
        One row per retained site, sorted by position then allele.

    Raises:
        DataContractError: No isolates, a malformed call, or a required upper
            bound that was not supplied.
    """
    if not calls_by_isolate:
        raise DataContractError(
            "Cohort merge received no isolates; a cohort of zero cannot be "
            "compared with itself",
            n_isolates=0,
        )
    cohort_size = len(calls_by_isolate)
    ceiling_af = max_maf
    if min_dissenters is not None:
        # `af > 1 - d/an` is exactly `an - ac < d`, so the symmetric rule and an
        # explicit frequency ceiling are the same bound; no separate path.
        symmetric = 1 - (min_dissenters / cohort_size)
        ceiling_af = symmetric if ceiling_af is None else min(ceiling_af, symmetric)

    if require_max_maf and ceiling_af is None:
        raise DataContractError(
            "max_maf is required here but was not supplied. An unset upper "
            "bound keeps every retained site as a minor allele, which holds at "
            "a handful of isolates and is wrong at full scale. Measure the "
            "allele-frequency spectrum and pass it explicitly - see "
            "docs/design/cohort-variant-merge.md.",
            n_isolates=len(calls_by_isolate),
        )

    carriers: Dict[Call, int] = {}
    for isolate, calls in calls_by_isolate.items():
        seen: set = set()
        for call in calls:
            chrom, pos, ref, alt = call
            if not (chrom and pos and ref and alt):
                raise DataContractError(
                    f"Variant call from {isolate!r} has no usable locus; it "
                    "cannot be merged or audited",
                    isolate=isolate,
                    call=str(call),
                )
            if call in seen:
                # One isolate cannot carry the same allele twice; a duplicate is
                # a caller bug, and counting it twice would inflate the allele
                # frequency for that site.
                continue
            seen.add(call)
            carriers[call] = carriers.get(call, 0) + 1

    rows: List[Dict[str, str]] = []
    for (chrom, pos, ref, alt), count in carriers.items():
        ac = count
        an = cohort_size
        # Allele frequency of THIS variant, not the minor allele of a biallelic
        # genotype. `min(ac, an - ac)` is wrong here: it reads the non-carriers
        # as the alternative allele, so a lone carrier at N=3 scores 1/3 - the
        # same as a genuine third of the cohort - and the floor stops filtering
        # anything. The carriers are the allele being measured; the rest simply
        # did not call it.
        af = ac / an
        if ac < 2:
            # A single carrier: a 1/N allele by definition, which is the
            # "at least two isolates" rule expressed without hard-coding a count
            # that would not scale from 4 to 835.
            continue
        if ceiling_af is not None and af > ceiling_af:
            # Too few isolates dissent. By the floor's argument - one carrier is
            # indistinguishable from one isolate's error - a site nearly every
            # isolate carries is as suspect as a site nearly none of them do.
            # This asks for support on both sides; it does not claim the site is
            # biologically fixed.
            continue
        rows.append(
            {
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "ac": str(ac),
                "an": str(an),
                "af": f"{af:.6g}",
            }
        )
    rows.sort(key=lambda r: (r["chrom"], int(r["pos"]) if r["pos"].isdigit() else 0, r["ref"], r["alt"]))
    return rows


def run(
    config: Any,
    manifest: Any,
    calls_by_isolate: Mapping[str, Sequence[Mapping[str, str]]],
) -> List[Dict[str, str]]:
    """Stage 6a entry point: merge the per-isolate calls cohort-wide.

    Args:
        config: Loaded pipeline configuration. Unused today - `max_maf` is
            deliberately not defaulted, so there is nothing to read yet - but
            it is taken so a ceiling decided in `science.yaml` can be threaded
            through here without changing the call sites.
        manifest: The cohort. Used for one thing only, and it is load-bearing:
            see below.
        calls_by_isolate: Per-isolate rows from the stage above, keyed by
            isolate.

    Returns:
        One row per retained site, in :data:`MERGE_COLUMNS` order.

    **The manifest, not the mapping, fixes the cohort size.** :func:`merge_calls`
    derives ``an = len(calls_by_isolate)``, so any isolate missing from that
    mapping is missing from the denominator - and its absence is invisible: no
    error, and a table that looks entirely well formed. Every manifest sample is
    therefore seeded, so an isolate whose assembly could not be aligned is
    counted as carrying no variant rather than vanishing.

    That is not tidiness. An isolate with no assembly is an ordinary member of
    the cohort, and a 20-isolate cohort in which 5 produced calls must report
    ``an=20``. Reporting ``an=5`` would inflate every allele frequency four-fold.
    """
    cohort: Dict[str, List[Call]] = {
        sample_id: [] for sample_id in manifest.sample_ids
    }

    for isolate, rows in calls_by_isolate.items():
        if isolate not in cohort:
            # Not a silent drop: a call from an isolate the manifest never named
            # means the two disagree about who is in the study, and this is
            # exactly the class of mismatch rule 5 makes a hard failure.
            raise DataContractError(
                f"Variant calls reference {isolate!r}, which is not in the "
                "manifest. The manifest decides cohort membership; a call from "
                "an isolate it does not name cannot be merged, because it would "
                "inflate the denominator with a sample the study does not contain.",
                isolate=isolate,
                n_manifest=len(cohort),
            )
        for row in rows:
            cohort[isolate].append(
                (str(row["chrom"]), str(row["pos"]), str(row["ref"]), str(row["alt"]))
            )

    # No `max_maf`, and no `require_max_maf`: the ceiling has now been measured
    # at N=4, 12 and 34 and remains unresolvable. The fixed-difference mode
    # *shrank* as N grew - impossible for a genuinely species-fixed set - and the
    # upper tail sits at 94%, not 100%, which is the signature of shared
    # assembly artefacts rather than of a biological difference. The symmetric
    # bound below is the operative filter. See
    # docs/design/cohort-variant-merge.md 5c.
    rows = merge_calls(cohort)
    LOGGER.info(
        "Stage 6a: %d polymorphic sites across %d isolates", len(rows), len(cohort)
    )
    return rows
