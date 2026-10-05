"""Stage 6 - per-isolate variant calls against the reference.

The contract here is what `bcftools mpileup | bcftools call -mv` actually emits,
measured on a real run rather than assumed. See `tests/unit/test_variants_contract.py`,
which records the run the columns came from.

Scope: **one isolate against the reference.** That is what per-isolate calling
means here, and it has a consequence worth stating. An assembly is a consensus
sequence, so depth is 1 at every locus and QUAL is a constant; more importantly,
these calls conflate species-wide fixed differences with cohort-informative
polymorphisms. A position where every isolate differs from PAO1 is not
polymorphic, and only a multi-isolate merge can tell the two apart. So this
stage's output is the *input to a merge*, not the GWAS's feature table.

Deliberately not carried from the VCF: `SGB`, a bcftools implementation detail
whose name would outlive its meaning. `QUAL` is kept for completeness, but no
threshold may be derived from it here - see the module docstring.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

from ..adapters import minimap2 as _aligner
from ..assemblies import locate_assembly
from ..config.loader import (
    DEFAULT_MPILEUP_MAX_DEPTH,
    DEFAULT_MPILEUP_OUTPUT_TYPE,
    PipelineConfig,
)
from ..errors import (
    DataContractError,
    PipelineError,
    ToolExecutionError,
)
from ..logging_utils import get_logger
from ..manifest import SampleManifest
from ..models import RunMode

#: Re-exported at module level so the stage and its tests refer to the caller by
#: one name. `run()` reaches it through this binding, which is what lets a test
#: substitute it - a `from ..adapters import minimap2` *inside* `_run_by_aligning`
#: would resolve at call time and ignore the substitution.
call_isolate = _aligner.call_isolate

LOGGER = get_logger("stages.variants")

#: The per-isolate call contract, in file order. Derived from a real run.
PER_ISOLATE_COLUMNS = (
    "sample_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "qual",
    "filter",
    "GT",
    "AC",
    "AN",
    "DP4",
    "MQ",
    "MQ0F",
)

#: INFO fields lifted from the VCF. `SGB` is excluded on purpose.
_INFO_FIELDS = ("AC", "AN", "DP4", "MQ", "MQ0F")

#: Sidecar written beside the stage's per-isolate calls. See
#: :func:`write_provenance`.
PROVENANCE_NAME = "variants_provenance.json"

#: Substring the zero-indel warning must contain. Named rather than inlined so
#: the test that greps the log and the code that emits it cannot drift.
NO_INDEL_WARNING = "suspicious for assembly input"


def _info(info_field: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for item in (info_field or "").split(";"):
        if not item or "=" not in item:
            continue
        key, _, value = item.partition("=")
        out[key] = value
    return out


def parse_vcf(text: str, sample_id: Optional[str]) -> List[Dict[str, str]]:
    """Parse per-isolate calls out of a VCF.

    A multiallelic site (``ALT=C,G``) becomes one row per ALT allele, because a
    comma-separated ALT is two calls sharing one line and keeping only the first
    would lose half of them silently.

    Args:
        text: VCF contents.
        sample_id: The isolate these calls belong to. Required: the isolate is
            the join key, and an unattributable call cannot be merged.

    Returns:
        One row per called allele, in file order.

    Raises:
        DataContractError: No sample id, or a record with no usable position.
    """
    if not sample_id:
        raise DataContractError(
            "Per-isolate variant calls must carry a sample_id; the isolate is "
            "the join key for the cohort merge"
        )

    rows: List[Dict[str, str]] = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 8:
            continue
        chrom, pos, _id, ref, alt, qual, filt, info = fields[:8]
        if pos in ("", "."):
            raise DataContractError(
                "Variant record has no position", chrom=chrom, alt=alt
            )
        parsed = _info(info)
        genotype = fields[9].split(":")[0] if len(fields) > 9 else ""
        for allele in (alt.split(",") if alt else []):
            if allele in ("", "."):
                continue
            row = {
                "sample_id": sample_id,
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": allele,
                "qual": qual,
                "filter": filt,
                "GT": genotype,
            }
            for key in _INFO_FIELDS:
                row[key] = parsed.get(key, "")
            rows.append(row)
    return rows


def load_precomputed(
    calls_dir: Path, sample_id: str
) -> List[Dict[str, str]]:
    """Parse one isolate's calls from an already-produced VCF.

    The TEST-mode path, and the ingestion route for externally generated calls.
    It runs the same :func:`parse_vcf` a real run does, so the awkward shapes
    real bcftools emits - a multiallelic ALT, ``DP4`` as a four-tuple - are
    exercised here rather than only in a parser unit test.

    Args:
        calls_dir: Directory holding ``<sample_id>.vcf``.
        sample_id: The isolate whose calls are wanted.

    Returns:
        The parsed rows, which is empty when the VCF carries no records.

    Raises:
        DataContractError: No VCF for this isolate. A missing file and an empty
            one are different findings, and only the second is a result.
    """
    path = Path(calls_dir) / f"{sample_id}.vcf"
    if not path.is_file():
        raise DataContractError(
            f"No variant calls found for {sample_id} at {path}. A missing VCF "
            "is not an empty result: it means the isolate was never called "
            "against the reference, which is a different finding from being "
            "called and carrying nothing.",
            sample_id=sample_id,
            path=str(path),
        )
    rows = parse_vcf(path.read_text(encoding="utf-8"), sample_id=sample_id)
    LOGGER.debug("Stage 6: %d precomputed calls for %s", len(rows), sample_id)
    return rows


# ---------------------------------------------------------------------------
# Provenance: how many SNVs and how many indels, per isolate.
# ---------------------------------------------------------------------------
#
# Why this exists. The stage's contract is a row per called ALLELE, and it
# carries no column saying what KIND of allele it was - that is the consumer's
# inference from `len(ref)`/`len(alt)`. Which means the stage can emit 2,447
# rows, 507 of them 1 nt, and say nothing whatsoever about the fact that it
# called ZERO indels for six consecutive rounds. A pre-round-11 REAL run over
# the 10-isolate smoke subset called 539,527 alleles of which exactly 28 were
# non-SNV, none inside oprD, and the table looked entirely well formed.
#
# The 28/539,527 ratio is the whole finding and it is invisible in the output.
# So it is recorded here, per isolate, next to the calls.

#: The provenance keys, in file order. `n_alleles` is `n_snv + n_indel`, stated
#: rather than summed by the reader so a mismatch is visible.
PROVENANCE_COLUMNS: tuple = ("sample_id", "n_calls", "n_snv", "n_indel", "n_alleles")


def is_indel(ref: str, alt: str) -> bool:
    """Whether one called allele is not a single-base substitution.

    A substitution is one reference base for one alternate base; anything else
    changes the length of the sequence, so it is an indel by the only
    definition available here. An MNV (``AC>GT``) is therefore an indel, and
    says so rather than being quietly bucketed with the substitutions - it is a
    0 nt length change that is not a substitution, and the two cases call for
    different downstream handling.

    Args:
        ref: The REF allele.
        alt: One ALT allele. Already split on commas by :func:`parse_vcf`.

    Returns:
        True when the allele changes length.
    """
    return len(ref) != 1 or len(alt) != 1


def count_alleles(rows: Sequence[Mapping[str, str]]) -> Dict[str, int]:
    """``n_snv`` / ``n_indel`` for one isolate's rows.

    Args:
        rows: The isolate's rows, in :data:`PER_ISOLATE_COLUMNS` form.

    Returns:
        Counts keyed ``n_snv``, ``n_indel``, ``n_alleles``. Every key is always
        present, including at zero: an absent count is indistinguishable from
        a count that was never taken.
    """
    n_snv = 0
    for row in rows:
        if not is_indel(str(row.get("ref") or ""), str(row.get("alt") or "")):
            n_snv += 1
    n_indel = len(rows) - n_snv
    return {"n_snv": n_snv, "n_indel": n_indel, "n_alleles": len(rows)}


def call_provenance(
    calls_by_isolate: Mapping[str, Sequence[Mapping[str, str]]],
) -> Dict[str, Any]:
    """The per-isolate SNV/indel accounting for a whole stage run.

    One entry per isolate **in the mapping**, in sorted order, so the record is
    byte-stable across runs of the same cohort. An isolate with no calls gets
    ``n_snv: 0, n_indel: 0`` rather than no entry, which is the same reasoning
    :func:`run` applies to its return value: "called and found nothing" and
    "never screened" must not look alike.

    Args:
        calls_by_isolate: ``sample_id -> rows``.

    Returns:
        ``{"per_isolate": [...], "totals": {...}}``. The totals are recomputed
        from the per-isolate entries rather than accumulated alongside them, so
        the two cannot disagree.
    """
    per_isolate: List[Dict[str, Any]] = []
    for sample_id in sorted(calls_by_isolate):
        counts = count_alleles(calls_by_isolate[sample_id])
        per_isolate.append(
            {
                "sample_id": sample_id,
                "n_calls": len(calls_by_isolate[sample_id]),
                **counts,
            }
        )

    totals = {
        "isolates": len(per_isolate),
        "n_calls": sum(e["n_calls"] for e in per_isolate),
        "n_snv": sum(e["n_snv"] for e in per_isolate),
        "n_indel": sum(e["n_indel"] for e in per_isolate),
    }
    totals["n_alleles"] = totals["n_snv"] + totals["n_indel"]
    return {"per_isolate": per_isolate, "totals": totals}


def write_provenance(
    calls_by_isolate: Mapping[str, Sequence[Mapping[str, str]]],
    path: Path,
) -> Dict[str, Any]:
    """Write the accounting to ``path`` as JSON, and return it.

    Best effort in exactly one direction: the sidecar is written, but a failure
    to write it is logged and does not fail the stage. The calls are the result;
    the provenance is the receipt, and losing the receipt must not lose the
    result. The reverse - logging a receipt that was never written - is why the
    write is attempted before the summary is logged.

    Args:
        calls_by_isolate: ``sample_id -> rows``.
        path: Where to write. Parent directories are created.

    Returns:
        The same payload :func:`call_provenance` built, so a caller that only
        wants the numbers never has to read the file back.
    """
    payload = call_provenance(calls_by_isolate)
    path = Path(path)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n",
                        encoding="utf-8")
    except OSError as exc:
        LOGGER.warning(
            "Stage 6: could not write the call provenance sidecar at %s (%s). "
            "The per-isolate SNV/indel counts are in the log line below, but "
            "nothing on disk records them.",
            path, exc,
        )
    return payload


def summarise_call_provenance(
    calls_by_isolate: Mapping[str, Sequence[Mapping[str, str]]],
    *,
    mode: RunMode,
    path: Optional[Path] = None,
) -> Dict[str, Any]:
    """Record the SNV/indel accounting, and warn when there are no indels at all.

    The warning exists because a REAL assembly run that calls thousands of
    substitutions and not one indel is far more likely to be a broken caller
    than a clean cohort. It is a WARNING and not a refusal: the run's result is
    legitimate, and this stage does not get to overrule it. But silence is not
    an option either - the failure this catches produced a clean-looking table
    for six rounds.

    The condition is deliberately the **cohort total**, not per isolate. A
    single isolate with no indels is ordinary; ten of ten with none, on input
    that is by construction a set of draft assemblies each carrying its own
    indels relative to PAO1, is the signature of a filter that cannot see them.

    Args:
        calls_by_isolate: ``sample_id -> rows``.
        mode: The run's mode. The warning fires for ``REAL`` only; TEST reads
            committed fixtures whose contents are a property of the fixture,
            not of this pipeline.
        path: Optional sidecar destination. See :func:`write_provenance`.

    Returns:
        The payload, so a caller can assert on the numbers without re-reading.
    """
    payload = (
        write_provenance(calls_by_isolate, path)
        if path is not None
        else call_provenance(calls_by_isolate)
    )
    totals = payload["totals"]

    LOGGER.info(
        "Stage 6: %d alleles across %d isolates (%d SNV, %d indel)",
        totals["n_alleles"], totals["isolates"],
        totals["n_snv"], totals["n_indel"],
    )

    if mode is RunMode.REAL and totals["n_indel"] == 0:
        LOGGER.warning(
            "Stage 6: ZERO indels called across all %d isolates, against %d "
            "SNVs. An assembly is a consensus sequence, so every isolate "
            "carries its own indels relative to the reference and at least "
            "some should be callable; a run with substitutions and no indels "
            "at all is %s, and the loci concerned will read as clean. Check "
            "bcftools mpileup --min-ireads (an assembly supplies ONE gapped "
            "read per indel, so the tool default of 2 discards every one of "
            "them) before treating this as a result.",
            totals["isolates"], totals["n_snv"], NO_INDEL_WARNING,
        )
    return payload


def call_settings(config: PipelineConfig) -> Dict[str, Any]:
    """The ``variants`` settings block, with every key the adapter reads resolved.

    Exists so the whole configuration-to-argv chain is one named, testable step
    rather than a dozen lines buried in :func:`_run_by_aligning`. Ruling R6's
    ``min_ireads`` is the key that has to survive it: parsed in
    ``science.yaml``, read by :meth:`PipelineConfig.variants_mpileup_min_ireads`,
    copied here, and finally emitted by
    :func:`papipeline.adapters.minimap2.mpileup_command`.

    Args:
        config: Loaded pipeline configuration.

    Returns:
        A deep-enough copy of ``config.raw['variants']`` in which
        ``bcftools.mpileup.min_ireads`` is present and equal to the loader's
        answer, even when the raw block omits it. Every other key is passed
        through untouched.

    The copy is deliberate. ``config.raw`` is shared by every later stage in the
    run, so filling a default into the nested dict in place would edit the
    configuration out from under the rest of the pipeline.
    """
    settings = dict(config.raw.get("variants") or {})
    bcftools_cfg = dict(settings.get("bcftools") or {})
    mpileup_cfg = dict(bcftools_cfg.get("mpileup") or {})
    mpileup_cfg.setdefault("min_ireads", config.variants_mpileup_min_ireads())
    mpileup_cfg.setdefault("max_depth", DEFAULT_MPILEUP_MAX_DEPTH)
    mpileup_cfg.setdefault("output_type", DEFAULT_MPILEUP_OUTPUT_TYPE)
    bcftools_cfg["mpileup"] = mpileup_cfg
    settings["bcftools"] = bcftools_cfg
    return settings


def run(
    config: PipelineConfig,
    manifest: SampleManifest,
    mode: RunMode,
    *,
    data_root: Path,
    workdir: Path,
) -> Dict[str, List[Dict[str, str]]]:
    """Stage 6 entry point: per-isolate calls against the reference.

    Args:
        config: Loaded pipeline configuration.
        manifest: The cohort. Every sample here is a member of it.
        mode: ``TEST`` reads committed fixtures; ``REAL`` aligns and calls.
        data_root: Mode's data root. ``TEST`` finds its VCFs under
            ``<data_root>/variants/``.
        workdir: Scratch directory for intermediate alignments.

    Returns:
        ``sample_id -> rows``, with a key for **every** manifest sample. An
        empty list means "called, found nothing"; an absent key would mean
        "never screened", and collapsing the two would make an unanalysable
        genome indistinguishable from a clean one.

        The empty entries are load-bearing downstream, not tidiness:
        :func:`~papipeline.stages.cohort_variants.merge_calls` derives its
        denominator ``an`` from ``len(calls_by_isolate)``, so dropping the
        isolates that produced no calls would inflate every allele frequency the
        cohort reports.
    """
    if mode is RunMode.TEST:
        return _run_from_fixtures(manifest, data_root)

    return _run_by_aligning(config, manifest, data_root, workdir)


def _run_from_fixtures(
    manifest: SampleManifest, data_root: Path
) -> Dict[str, List[Dict[str, str]]]:
    """TEST: parse the committed per-isolate VCFs.

    There is no reference in `test_data/` and REAL mode is closed, so TEST
    cannot align anything. It can still drive the real parser over output shaped
    the way bcftools shapes it, which is the point of the committed fixtures -
    a TEST path that skipped `parse_vcf` would leave the awkward shapes untested.
    """
    calls_dir = Path(data_root) / "variants"
    if not calls_dir.is_dir():
        raise DataContractError(
            f"No variant calls directory at {calls_dir}. TEST mode reads the "
            "committed per-isolate VCFs from here; without them the stage has "
            "nothing to parse, and returning an empty cohort instead would be "
            "indistinguishable from a cohort in which nothing was found.",
            path=str(calls_dir),
        )

    calls: Dict[str, List[Dict[str, str]]] = {}
    for sample_id in manifest.sample_ids:
        if (calls_dir / f"{sample_id}.vcf").is_file():
            calls[sample_id] = load_precomputed(calls_dir, sample_id)
        else:
            # Present, empty: the isolate is a cohort member that carries no
            # called variant. `an` must still count it.
            calls[sample_id] = []

    called = sum(1 for rows in calls.values() if rows)
    LOGGER.info(
        "Stage 6: %d/%d isolates produced calls", called, len(calls),
    )
    # TEST writes no sidecar: there is no scratch workdir on that path, and the
    # counts of a committed fixture are a property of the fixture. The summary
    # is still logged, with the zero-indel warning suppressed - see
    # `summarise_call_provenance`.
    summarise_call_provenance(calls, mode=RunMode.TEST)
    return calls


def _run_by_aligning(
    config: PipelineConfig,
    manifest: SampleManifest,
    data_root: Path,
    workdir: Path,
) -> Dict[str, List[Dict[str, str]]]:
    """REAL: align each isolate to the reference and call its variants.

    Two failure modes are recorded rather than raised, and both keep the isolate
    in the cohort with no calls.

    **No assembly.** :func:`~papipeline.assemblies.locate_assembly` refuses, and
    that refusal is correct - substituting another sample's sequence would shift
    every coordinate silently. The isolate stays a cohort member carrying
    nothing, because :func:`merge_calls` derives ``an`` from the cohort size and
    dropping the isolate would inflate every allele frequency the merge reports.

    **The call itself failed.** 835 isolates were attempted and 34 were usable
    (docs/design/cohort-variant-merge.md 5c) - 82% of the assemblies are
    corrupt. A stage that stopped at the first failure could therefore never
    complete, and one that ignored the failures would report a cohort smaller
    than the study without saying so. So each failure is logged with its
    isolate and reason, and the count is reported at the end: an unanalysable
    genome is visible rather than quietly absent.
    """
    from ..adapters import minimap2 as aligner

    tools = aligner.require_callers({
        name: _which(name) for name in aligner.CALLER_TOOLS
    })
    reference = config.reference_fasta()
    settings = call_settings(config)
    threads = int(config.runtime.get("threads", 1) or 1)
    workdir = Path(workdir)

    calls: Dict[str, List[Dict[str, str]]] = {}
    unusable: List[str] = []

    for sample_id in manifest.sample_ids:
        try:
            assembly = locate_assembly(manifest.require(sample_id), data_root)
        except PipelineError as exc:
            LOGGER.warning(
                "Stage 6: no assembly for %s, recorded as a cohort member with "
                "no calls (%s)", sample_id, exc,
            )
            calls[sample_id] = []
            unusable.append(sample_id)
            continue

        try:
            result = call_isolate(
                sample_id,
                assembly=assembly,
                reference=reference,
                workdir=workdir / sample_id,
                settings=settings,
                tools=tools,
                threads=threads,
            )
        except (PipelineError, ToolExecutionError) as exc:
            LOGGER.warning(
                "Stage 6: calling failed for %s, recorded as a cohort member "
                "with no calls (%s)", sample_id, exc,
            )
            calls[sample_id] = []
            unusable.append(sample_id)
            continue

        # The same parser TEST uses, so the two modes cannot disagree about how
        # a call is read.
        calls[sample_id] = parse_vcf(
            Path(result.vcf).read_text(encoding="utf-8"), sample_id=sample_id
        )

    if unusable:
        LOGGER.warning(
            "Stage 6: %d of %d isolates produced no calls (%s). They remain "
            "cohort members, so the merge's denominator still counts them.",
            len(unusable), len(calls), ",".join(unusable[:10]),
        )
    summarise_call_provenance(
        calls, mode=RunMode.REAL, path=workdir / PROVENANCE_NAME
    )
    return calls


def _which(name: str) -> Optional[str]:
    """Locate one tool, honouring the ``<TOOL>_BIN`` override.

    An environment override exists because these tools usually live inside a
    micromamba environment that is not on the login ``PATH``.
    """
    import os

    override = os.environ.get(f"{name.upper()}_BIN")
    if override:
        return override
    import shutil

    return shutil.which(name)
