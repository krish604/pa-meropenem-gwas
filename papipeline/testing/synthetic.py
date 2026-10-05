"""Deterministic synthetic test-data generator.

============================= IMPORTANT =============================
EVERYTHING THIS MODULE PRODUCES IS SYNTHETIC TEST DATA.

IT IS NOT BIOLOGICAL RESULTS. IT IS NOT DERIVED FROM ANY REAL GENOME,
PHENOTYPE OR METADATA FILE. The statistical relationships between the
generated genotype, mechanism and phenotype columns are ARTEFACTS OF THE
GENERATOR SEED. They exist so that stages 12-15 have a signal to process and
so that assertions in the test suite are stable. They must never be read as
evidence about Pseudomonas aeruginosa.
===================================================================

Everything is derived from a single integer seed, so two runs with the same
seed produce byte-identical files. This is what makes the pipeline's own
outputs testable.

Generator design
----------------
The generator builds a deliberately simple, fully documented world:

* three lineages ``LINEAGE_A``/``B``/``C``, assigned round-robin;
* a latent binary ``resistance_module`` per sample, drawn from the seed;
* a mechanism profile (OprD disruption, efflux-regulator variants, an
  acquired carbapenemase determinant) drawn conditional on that module, so
  that stage 12 has a recoverable association and stage 13 has convergent
  structure to classify;
* a phenotype drawn from the module with deliberate noise, so that some
  resistant samples have no detected determinant and some susceptible
  samples do. This is required: if the synthetic data made gene presence
  perfectly predict resistance, the test suite would fail to exercise the
  very distinction (rule 1 in ``docs/scientific_rules.md``) that the pipeline
  exists to preserve.

The generator writes into ``test_data/`` only. It refuses to write into
``data/``.
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from ..errors import PipelineError
from ..io.tsv import write_fasta, write_tsv
from ..logging_utils import get_logger
from ..models import (
    DISRUPTIVE_VARIANT_TYPES,
    ClaimStatus,
    Phenotype,
    VariantType,
)
from ..stages.pangenome import GPA_COLUMNS

LOGGER = get_logger("testing.synthetic")

DEFAULT_SEED = 20240617
DEFAULT_N_SAMPLES = 20
SAMPLE_PREFIX = "TEST_PA_"
SYNTHETIC_BANNER = "SYNTHETIC TEST DATA - NOT BIOLOGICAL RESULTS"

LINEAGES: Tuple[str, ...] = ("LINEAGE_A", "LINEAGE_B", "LINEAGE_C")

#: Determinants the generator may emit. All are knowledge-table genes.
ACQUIRED_DETERMINANTS: Tuple[str, ...] = (
    "blaOXA-48",
    "blaOXA-23",
    "blaNDM-1",
    "blaVIM-1",
    "blaKPC-2",
    "blaTEM-1",
)

#: Genes the generator may report as reduced/absent (never "resistant").
LOSS_OF_FUNCTION_GENES: Tuple[str, ...] = ("oprD",)

#: Fallback screened variant classes per gene, used when no regulators table
#: is supplied. Mirrors config/regulators.tsv; the table is preferred so the
#: knowledge file stays the single source of truth.
DEFAULT_VARIANT_CHOICES: Dict[str, Tuple[str, ...]] = {
    "oprD": ("SNV", "FRAMESHIFT", "PREMATURE_STOP", "GENE_ABSENCE", "GENE_DISRUPTION"),
    "mexR": ("SNV", "FRAMESHIFT", "PREMATURE_STOP", "GENE_ABSENCE"),
    "nalC": ("SNV", "FRAMESHIFT", "PREMATURE_STOP", "GENE_ABSENCE"),
    "nalD": ("SNV", "FRAMESHIFT", "PREMATURE_STOP", "GENE_ABSENCE"),
    "mexZ": ("SNV", "FRAMESHIFT", "PREMATURE_STOP", "GENE_ABSENCE"),
    "nfxB": ("SNV", "FRAMESHIFT", "PREMATURE_STOP", "GENE_ABSENCE"),
    "mexT": ("SNV", "INSERTION", "FRAMESHIFT", "GENE_DISRUPTION"),
    "mexS": ("SNV", "FRAMESHIFT", "PREMATURE_STOP"),
    "ampD": ("SNV", "FRAMESHIFT", "GENE_ABSENCE"),
    "ampR": ("SNV", "FRAMESHIFT", "PREMATURE_STOP"),
    "dacB": ("SNV", "FRAMESHIFT", "PREMATURE_STOP", "GENE_ABSENCE"),
}

#: Variant classes the generator will use for a "variant present" call.
#: GENE_ABSENCE is handled separately because it also drives OprD status.
_EMITTABLE_VARIANT_TYPES: Tuple[str, ...] = (
    "SNV",
    "FRAMESHIFT",
    "PREMATURE_STOP",
    "GENE_ABSENCE",
)


def variant_choices_from_table(path: Path) -> Dict[str, Tuple[str, ...]]:
    """Read the screened variant classes per gene from ``config/regulators.tsv``.

    Keeping the generator driven by the knowledge table means a fixture can
    never emit a variant class the screen does not declare.
    """
    from ..io.tsv import read_tsv

    if not Path(path).exists():
        return dict(DEFAULT_VARIANT_CHOICES)
    rows = read_tsv(path, required_columns=("gene", "screen_for"))
    choices: Dict[str, Tuple[str, ...]] = {}
    for row in rows:
        gene = row["gene"]
        screened = tuple(
            v.strip()
            for v in (row.get("screen_for") or "").split(",")
            if v.strip()
        )
        emittable = tuple(v for v in screened if v in _EMITTABLE_VARIANT_TYPES)
        if emittable:
            choices[gene] = emittable
    return choices or dict(DEFAULT_VARIANT_CHOICES)


def disruptive_variant_choices_from_table(path: Path) -> Dict[str, Tuple[str, ...]]:
    """Per-gene *disruptive* classes the screen actually declares.

    Reads the raw ``screen_for`` column rather than the
    ``_EMITTABLE_VARIANT_TYPES``-filtered choices, because that filter exists to
    drop classes with no ``VariantType`` member (``INSERTION``) and was never
    meant to decide what counts as disruptive.

    ``GENE_ABSENCE`` is excluded: it is loss of function too, but the generator
    reports it as its own OprD state, and folding it in here would collapse
    absence into disruption - the exact distinction the
    ``oprD_absent`` / ``oprD_LoF`` pair exists to keep.

    Sharing :data:`papipeline.models.DISRUPTIVE_VARIANT_TYPES` with the stage
    that consumes these rows is the point. The generator and the stage now agree
    by construction on what "disruptive" means, so a fixture can exercise the
    same code path real calls take.
    """
    from ..io.tsv import read_tsv

    if not Path(path).exists():
        return {g: _DEFAULT_DISRUPTIVE for g in DEFAULT_VARIANT_CHOICES}

    rows = read_tsv(path, required_columns=("gene", "screen_for"))
    choices: Dict[str, Tuple[str, ...]] = {}
    for row in rows:
        gene = row["gene"]
        screened = (
            v.strip() for v in (row.get("screen_for") or "").split(",") if v.strip()
        )
        disruptive = tuple(
            v
            for v in screened
            if v in DISRUPTIVE_VARIANT_TYPES and v != VariantType.GENE_ABSENCE.value
        )
        if disruptive:
            choices[gene] = disruptive
    return choices


#: Fallback for :func:`disruptive_variant_choices_from_table`, mirroring
#: ``config/regulators.tsv``. oprD's three disruptive classes, in the order the
#: generator cycles through them.
_DEFAULT_DISRUPTIVE: Tuple[str, ...] = (
    "FRAMESHIFT",
    "PREMATURE_STOP",
    "GENE_DISRUPTION",
)

#: Regulator loci the generator may place a variant in.
REGULATOR_GENES: Tuple[str, ...] = tuple(sorted(DEFAULT_VARIANT_CHOICES))

VIRULENCE_FACTORS: Tuple[Tuple[str, str, str], ...] = (    ("vfh", "VIR_FACT", "adhesion"),
    ("exoS", "VIR_FACT", "cytotoxin"),
    ("exoT", "VIR_FACT", "cytotoxin"),
    ("lasB", "VIR_FACT", "protease"),
    ("pyocyanin_synthesis", "VIR_FACT", "pigment"),
    ("pilA", "VIR_FACT", "adhesion"),
)

#: Miniature contigs. These are deliberately far smaller than a real
#: P. aeruginosa genome so the fixtures stay small and fast. Stage 1 QC is
#: expected to flag them; that is correct behaviour, not a failure.
N_CONTIGS = 5
CONTIG_LENGTH = 2000


def sample_id_for(index: int) -> str:
    """``1 -> TEST_PA_001``."""
    return f"{SAMPLE_PREFIX}{index:03d}"


@dataclass
class SyntheticSample:
    """The latent truth for one synthetic sample, plus what was emitted."""

    index: int
    sample_id: str
    lineage: str
    resistance_module: bool
    oprd_status: str
    regulator_variant_gene: Optional[str]
    regulator_variant_type: Optional[str]
    acquired_determinant: Optional[str]
    virulence_genes: List[str]
    st: int
    phenotype: Phenotype
    mic: Optional[float] = None
    mic_unit: Optional[str] = None

    def to_row(self) -> Dict[str, object]:
        return {
            "sample_id": self.sample_id,
            "lineage": self.lineage,
            "resistance_module": "true" if self.resistance_module else "false",
            "oprD_status": self.oprd_status,
            "regulator_variant_gene": self.regulator_variant_gene,
            "regulator_variant_type": self.regulator_variant_type,
            "acquired_determinant": self.acquired_determinant,
            "phenotype": self.phenotype.value,
        }


# --------------------------------------------------------------------------
# Latent world
# --------------------------------------------------------------------------


def build_samples(
    n_samples: int = DEFAULT_N_SAMPLES,
    seed: int = DEFAULT_SEED,
    variant_choices: Optional[Dict[str, Tuple[str, ...]]] = None,
) -> List[SyntheticSample]:
    """Construct the latent synthetic cohort deterministically.

    Args:
        n_samples: Number of synthetic samples.
        seed: RNG seed. Same seed gives byte-identical output.
        variant_choices: Per-gene screened variant classes, normally read
            from ``config/regulators.tsv``. A variant class the screen does
            not declare is never emitted.
    """
    rng = random.Random(seed)
    choices = dict(variant_choices or DEFAULT_VARIANT_CHOICES)
    samples: List[SyntheticSample] = []

    for index in range(1, n_samples + 1):
        sid = sample_id_for(index)
        lineage = LINEAGES[(index - 1) % len(LINEAGES)]
        module = rng.random() < 0.45

        # OprD: present in ~85% overall. Disruption correlates with, but does
        # not determine, the module.
        if module:
            oprd = rng.choices(
                ["disrupted", "absent", "intact"], weights=[0.60, 0.15, 0.25]
            )[0]
        else:
            oprd = rng.choices(
                ["disrupted", "absent", "intact"], weights=[0.12, 0.03, 0.85]
            )[0]

        # At most one regulator variant, and only a class that gene screens.
        reg_gene: Optional[str] = None
        reg_type: Optional[str] = None
        pool = [g for g in REGULATOR_GENES if g != "oprD" and choices.get(g)]
        if module and pool and rng.random() < 0.55:
            reg_gene = rng.choice(pool)
            reg_type = rng.choice(choices[reg_gene])
        elif (not module) and pool and rng.random() < 0.18:
            reg_gene = rng.choice(pool)
            reg_type = rng.choice(choices[reg_gene])

        acquired: Optional[str] = None
        if module and rng.random() < 0.30:
            acquired = rng.choice(ACQUIRED_DETERMINANTS)
        elif (not module) and rng.random() < 0.06:
            # Deliberate false-positive-ish determinants: present without
            # resistance. Exists so the pipeline must not equate the two.
            acquired = rng.choice(ACQUIRED_DETERMINANTS)

        virulence = [
            name
            for name, _gene, _cat in VIRULENCE_FACTORS
            if rng.random() < 0.55
        ]

        # Phenotype: module drives R, with noise in both directions, plus a
        # small number of I / SDD / ND so every allowed value is exercised.
        roll = rng.random()
        if index % 17 == 0:
            phenotype = Phenotype.ND
        elif index % 13 == 0:
            phenotype = Phenotype.SDD
        elif module:
            phenotype = (
                Phenotype.R if roll < 0.72 else rng.choice([Phenotype.I, Phenotype.S])
            )
        else:
            phenotype = (
                Phenotype.S if roll < 0.78 else rng.choice([Phenotype.I, Phenotype.R])
            )

        mic = mic_unit = None
        if phenotype in (Phenotype.R, Phenotype.I) and rng.random() < 0.30:
            # Only a minority get a measured MIC, and only for R/I. Never
            # derived from the category: drawn independently, then floored
            # into a plausible range for the category it accompanies.
            mic = round(rng.uniform(4.0, 16.0) if phenotype is Phenotype.R else rng.uniform(1.5, 4.0), 1)
            mic_unit = "mg/L"

        samples.append(
            SyntheticSample(
                index=index,
                sample_id=sid,
                lineage=lineage,
                resistance_module=module,
                oprd_status=oprd,
                regulator_variant_gene=reg_gene,
                regulator_variant_type=reg_type,
                acquired_determinant=acquired,
                virulence_genes=virulence,
                st=rng.randint(1, 12),
                phenotype=phenotype,
                mic=mic,
                mic_unit=mic_unit,
            )
        )

    _ensure_oprd_coverage(samples)
    return samples


def _ensure_oprd_coverage(samples: List[SyntheticSample]) -> None:
    """Guarantee the cohort exercises all three OprD states.

    The OprD *absent* state is the one that carries the reduced-permeability
    mechanism, so a seed that happened to produce none would leave the most
    important branch of stages 5, 6, 7, 13 and 15 untested. If a state is
    missing, the lowest-indexed eligible sample is set to it. The result is
    still fully deterministic; only the seed-to-cohort mapping is nudged.
    """
    present = {s.oprd_status for s in samples}
    for state in ("absent", "disrupted", "intact"):
        if state in present:
            continue
        for sample in samples:
            if sample.oprd_status == "intact":
                sample.oprd_status = state
                break


# --------------------------------------------------------------------------
# File writers
# --------------------------------------------------------------------------


def _banner(extra: Sequence[str] = ()) -> List[str]:
    return [SYNTHETIC_BANNER, f"generator_seed={DEFAULT_SEED}", *extra]


def _mini_fasta(sample: SyntheticSample, rng_seed: int) -> List[Tuple[str, str]]:
    """Build a miniature, clearly synthetic FASTA for one sample."""
    rng = random.Random(rng_seed)
    records: List[Tuple[str, str]] = []
    for contig in range(1, N_CONTIGS + 1):
        sequence = "".join(rng.choice("ACGT") for _ in range(CONTIG_LENGTH))
        header = (
            f"{sample.sample_id}_contig{contig} SYNTHETIC_TEST_DATA "
            f"length={CONTIG_LENGTH} topology=linear"
        )
        records.append((header, sequence))
    return records


def write_metadata(
    out_dir: Path, samples: Sequence[SyntheticSample]
) -> Dict[str, Path]:
    """Write ``test_data/metadata/sample_metadata.tsv``."""
    path = out_dir / "metadata" / "sample_metadata.tsv"
    rows = [
        {
            "sample_id": s.sample_id,
            "assembly_path": f"../genomes/{s.sample_id}.fna",
            "lineage": s.lineage,
            "data_class": "synthetic_test",
        }
        for s in samples
    ]
    write_tsv(
        path,
        rows,
        ["sample_id", "assembly_path", "lineage", "data_class"],
        header_comment=_banner(["column lineage is a generator label, not a result"]),
    )
    return {"metadata": path}


def write_genomes(
    out_dir: Path, samples: Sequence[SyntheticSample], seed: int = DEFAULT_SEED
) -> Dict[str, Path]:
    """Write one miniature FASTA per synthetic sample."""
    written: Dict[str, Path] = {}
    for sample in samples:
        path = out_dir / "genomes" / f"{sample.sample_id}.fna"
        write_fasta(path, _mini_fasta(sample, seed + sample.index))
        written[sample.sample_id] = path
    return written


def write_phenotype(
    out_dir: Path, samples: Sequence[SyntheticSample], antibiotic: str = "imipenem"
) -> Dict[str, Path]:
    """Write ``test_data/phenotype/imipenem_phenotype.tsv``."""
    path = out_dir / "phenotype" / f"{antibiotic}_phenotype.tsv"
    rows = [
        {
            "sample_id": s.sample_id,
            "antibiotic": antibiotic,
            "phenotype": s.phenotype.value,
            "MIC": s.mic,
            "MIC_unit": s.mic_unit,
            "source": "synthetic_generator",
        }
        for s in samples
    ]
    write_tsv(
        path,
        rows,
        ["sample_id", "antibiotic", "phenotype", "MIC", "MIC_unit", "source"],
        header_comment=_banner(
            [
                "phenotype values are randomly assigned by the generator",
                "MIC is only populated where the generator drew one independently",
                "R/I/S are never converted into MIC",
            ]
        ),
    )
    return {"phenotype": path}


def write_annotation(
    out_dir: Path, samples: Sequence[SyntheticSample]
) -> Dict[str, Path]:
    """Write a per-sample annotation TSV in the standardised schema."""
    rng = random.Random(DEFAULT_SEED + 7)
    paths: Dict[str, Path] = {}
    for sample in samples:
        path = (
            out_dir
            / "intermediate"
            / "annotation"
            / f"{sample.sample_id}.annotation.tsv"
        )
        genes = [
            ("oprD", "outer membrane protein D", "porin", 400, 1723, "+"),
            ("mexR", "transcriptional regulator MexR", "regulator", 200, 941, "-"),
            ("mexZ", "transcriptional regulator MexZ", "regulator", 100, 1153, "+"),
            ("nfxB", "transcriptional regulator NfxB", "regulator", 300, 1596, "-"),
            ("ampD", "D-alanyl-D-alanine carboxypeptidase AmpD", "enzyme", 50, 485, "+"),
            ("dacB", "D-alanyl-D-alanine carboxypeptidase DacB", "enzyme", 100, 1063, "-"),
        ]
        if sample.oprd_status == "absent":
            genes = [g for g in genes if g[0] != "oprD"]
        rows = []
        for index, (name, product, gene_type, start, end, strand) in enumerate(genes, 1):
            rows.append(
                {
                    "sample_id": sample.sample_id,
                    "contig_id": f"{sample.sample_id}_contig1",
                    "gene_id": f"{sample.sample_id}_g{index:03d}",
                    "gene_name": name,
                    "product": product,
                    "gene_type": gene_type,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "annotation_source": "synthetic_generator",
                }
            )
        if sample.acquired_determinant:
            rows.append(
                {
                    "sample_id": sample.sample_id,
                    "contig_id": f"{sample.sample_id}_contig2",
                    "gene_id": f"{sample.sample_id}_g900",
                    "gene_name": sample.acquired_determinant,
                    "product": f"synthetic determinant {sample.acquired_determinant}",
                    "gene_type": "acquired",
                    "start": 100,
                    "end": 900,
                    "strand": "+",
                    "annotation_source": "synthetic_generator",
                }
            )
        write_tsv(
            path,
            rows,
            [
                "sample_id",
                "contig_id",
                "gene_id",
                "gene_name",
                "product",
                "gene_type",
                "start",
                "end",
                "strand",
                "annotation_source",
            ],
            header_comment=_banner(),
        )
        paths[sample.sample_id] = path
    LOGGER.debug("rng consumed for annotation: %s", rng.random())
    return paths


def write_mlst(out_dir: Path, samples: Sequence[SyntheticSample]) -> Dict[str, Path]:
    """Write the MLST table, including one deliberately incomplete call."""
    loci = ["abcZ", "adk", "aroE", "guaA", "gyrB", "ircD", "toxA"]
    rows = []
    for sample in samples:
        # Every 7th sample gets an incomplete profile, to exercise
        # MLST_status=partial and to keep "no call" distinct from "bad call".
        complete = (sample.index % 7) != 0
        alleles = {
            locus: str((sample.st * (i + 3)) % 40 + 1)
            for i, locus in enumerate(loci if complete else loci[:4])
        }
        rows.append(
            {
                "sample_id": sample.sample_id,
                "ST": str(sample.st) if complete else "NA",
                "alleles": ";".join(f"{k}:{v}" for k, v in sorted(alleles.items())),
                "MLST_status": "typed" if complete else "partial",
                "mlst_scheme": "pseudomonas_aeruginosa",
                "allele_database": "synthetic_scheme_v0",
            }
        )
    path = out_dir / "intermediate" / "mlst" / "mlst_results.tsv"
    write_tsv(
        path,
        rows,
        ["sample_id", "ST", "alleles", "MLST_status", "mlst_scheme", "allele_database"],
        header_comment=_banner(),
    )
    return {"mlst": path}


def write_amr(out_dir: Path, samples: Sequence[SyntheticSample], antibiotic: str = "imipenem") -> Dict[str, Path]:
    """Write stage 4 AMR calls.

    Note that an acquired determinant is emitted for some susceptible
    samples and absent for some resistant ones. This is intentional: it is
    the fixture that proves the pipeline does not equate detection with
    resistance.
    """
    rows = []
    for sample in samples:
        if not sample.acquired_determinant:
            continue
        rows.append(
            {
                "sample_id": sample.sample_id,
                "antibiotic": antibiotic,
                "determinant": f"{sample.acquired_determinant}~SYNTHETIC_ALLELE",
                "gene": sample.acquired_determinant,
                "variant": "SYNTHETIC_ALLELE",
                "determinant_type": "acquired",
                "mechanism": "acquired_determinant",
                "evidence_source": "synthetic_amrfinder_emulator",
                "database": "SYNTHETIC_DB",
                "database_version": "v0-synthetic",
                "confidence": "99.0",
                "claim_status": ClaimStatus.DETECTED.value,
                "identity_pct": "99.5",
                "coverage_pct": "100.0",
            }
        )
    path = out_dir / "intermediate" / "amr" / "amr_determinants.tsv"
    write_tsv(
        path,
        rows,
        [
            "sample_id",
            "antibiotic",
            "determinant",
            "gene",
            "variant",
            "determinant_type",
            "mechanism",
            "evidence_source",
            "database",
            "database_version",
            "confidence",
            "claim_status",
            "identity_pct",
            "coverage_pct",
        ],
        header_comment=_banner(
            ["detection does not imply phenotypic resistance"]
        ),
    )
    return {"amr": path}


#: Free-text ``effect`` per disruptive class, so the fixture does not claim an
#: insertion sequence caused every disruption it reports.
_OPRD_DISRUPTION_EFFECT: Dict[str, str] = {
    "FRAMESHIFT": "frameshift_truncates_coding_sequence",
    "PREMATURE_STOP": "premature_stop_truncates_coding_sequence",
    "GENE_DISRUPTION": "disrupted_by_insertion",
}


def write_regulator_variants(
    out_dir: Path,
    samples: Sequence[SyntheticSample],
    oprd_disruptive: Sequence[str] = _DEFAULT_DISRUPTIVE,
) -> Dict[str, Path]:
    """Write stage 6 variant calls, including OprD loss-of-function.

    A disrupted OprD is emitted as one of the *disruptive* classes the screen
    declares for oprD - FRAMESHIFT, PREMATURE_STOP or GENE_DISRUPTION - cycling
    deterministically by sample index.

    It used to be emitted as GENE_DISRUPTION unconditionally, which made
    FRAMESHIFT and PREMATURE_STOP unreachable at oprD in TEST mode. Stage 6's
    ``gene_status()`` matched those two classes by literal string against
    GENE_ABSENCE and GENE_DISRUPTION only, so it labelled them "variant" and
    ``oprD_LoF`` came out 0. The emulator was emitting precisely the two class
    names the buggy branch understood, which is why the defect was invisible
    until real calls arrived.

    Cycling by index rather than drawing from ``rng`` is deliberate: consuming a
    draw here would shift the shared RNG stream and rewrite every downstream
    fixture, burying this change in unrelated churn.
    """
    classes = tuple(oprd_disruptive) or _DEFAULT_DISRUPTIVE
    rows: List[Dict[str, object]] = []
    for sample in samples:
        if sample.oprd_status in ("disrupted", "absent"):
            if sample.oprd_status == "disrupted":
                variant_type = classes[(sample.index - 1) % len(classes)]
                effect = _OPRD_DISRUPTION_EFFECT.get(
                    variant_type, "coding_sequence_disrupted"
                )
            else:
                variant_type = VariantType.GENE_ABSENCE.value
                effect = "locus_absent"
            rows.append(
                {
                    "sample_id": sample.sample_id,
                    "gene": "oprD",
                    "variant": f"SYNTHETIC_OPRD_{variant_type}_001",
                    "variant_type": variant_type,
                    "position": 412,
                    "reference": "SYN",
                    "alternate": "SYN",
                    "effect": effect,
                    "mechanism": "reduced_permeability",
                    "confidence": "SYNTHETIC",
                    "evidence_source": "synthetic_variant_emulator",
                    "call_status": ClaimStatus.DETECTED.value,
                }
            )
        if sample.regulator_variant_gene and sample.regulator_variant_type:
            rows.append(
                {
                    "sample_id": sample.sample_id,
                    "gene": sample.regulator_variant_gene,
                    "variant": f"SYNTHETIC_{sample.regulator_variant_gene}_V1",
                    "variant_type": sample.regulator_variant_type,
                    "position": 100 + (sample.index * 7) % 900,
                    "reference": "A",
                    "alternate": "T" if sample.regulator_variant_type == "SNV" else "SYN",
                    "effect": "synthetic_effect",
                    "mechanism": "efflux",
                    "confidence": "SYNTHETIC",
                    "evidence_source": "synthetic_variant_emulator",
                    "call_status": ClaimStatus.DETECTED.value,
                }
            )
    path = out_dir / "intermediate" / "regulators" / "regulator_variants.tsv"
    write_tsv(
        path,
        rows,
        [
            "sample_id",
            "gene",
            "variant",
            "variant_type",
            "position",
            "reference",
            "alternate",
            "effect",
            "mechanism",
            "confidence",
            "evidence_source",
            "call_status",
        ],
        header_comment=_banner(["a variant call is not a resistance call"]),
    )
    return {"regulators": path}


def write_structural_variants(out_dir: Path, samples: Sequence[SyntheticSample]) -> Dict[str, Path]:
    """Write stage 7 calls, deliberately mixing all three call statuses."""
    rows: List[Dict[str, object]] = []
    for sample in samples:
        if sample.oprd_status == "disrupted":
            rows.append(
                {
                    "sample_id": sample.sample_id,
                    "variant_id": f"{sample.sample_id}_SV_OPRD_INS",
                    "variant_type": "insertion_sequence",
                    "position": 415,
                    "affected_gene": "oprD",
                    "size": 1200,
                    "evidence": "synthetic_contig_breakpoint",
                    "confidence": "high",
                    "call_status": "confirmed",
                    "mge": "SYN_IS",
                }
            )
        elif sample.index % 4 == 0:
            # Candidate call: must stay `candidate` and must never be promoted.
            rows.append(
                {
                    "sample_id": sample.sample_id,
                    "variant_id": f"{sample.sample_id}_SV_CAND_001",
                    "variant_type": "deletion",
                    "position": 5000,
                    "affected_gene": "mexXY",
                    "size": 900,
                    "evidence": "synthetic_low_support",
                    "confidence": "low",
                    "call_status": "candidate",
                    "mge": None,
                }
            )
        elif sample.index % 6 == 0:
            rows.append(
                {
                    "sample_id": sample.sample_id,
                    "variant_id": f"{sample.sample_id}_SV_NA_001",
                    "variant_type": "deletion",
                    "position": None,
                    "affected_gene": "mexAB",
                    "size": None,
                    "evidence": "synthetic_contig_gap",
                    "confidence": None,
                    "call_status": "not_assessable",
                    "mge": None,
                }
            )
    path = out_dir / "intermediate" / "structural_variants" / "structural_variants.tsv"
    write_tsv(
        path,
        rows,
        [
            "sample_id",
            "variant_id",
            "variant_type",
            "position",
            "affected_gene",
            "size",
            "evidence",
            "confidence",
            "call_status",
            "mge",
        ],
        header_comment=_banner(
            ["call_status distinguishes confirmed, candidate and not_assessable"]
        ),
    )
    return {"structural_variants": path}


def write_virulence(out_dir: Path, samples: Sequence[SyntheticSample]) -> Dict[str, Path]:
    """Write stage 8 virulence detections (independent of AMR by design)."""
    rows: List[Dict[str, object]] = []
    for sample in samples:
        for name, _gene, category in VIRULENCE_FACTORS:
            if name in sample.virulence_genes:
                rows.append(
                    {
                        "sample_id": sample.sample_id,
                        "virulence_factor": name,
                        "gene": name,
                        "category": category,
                        "database": "SYNTHETIC_VFDB",
                        "database_version": "v0-synthetic",
                        "confidence": "95.0",
                        "identity_pct": "97.0",
                    }
                )
    path = out_dir / "intermediate" / "virulence" / "virulence_factors.tsv"
    write_tsv(
        path,
        rows,
        [
            "sample_id",
            "virulence_factor",
            "gene",
            "category",
            "database",
            "database_version",
            "confidence",
            "identity_pct",
        ],
        header_comment=_banner(["virulence stage is independent of the AMR stage"]),
    )
    return {"virulence": path}


def write_tree(out_dir: Path, samples: Sequence[SyntheticSample]) -> Dict[str, Path]:
    """Write a synthetic newick tree with one tip per sample, grouped by lineage.

    The tip labels are exactly the manifest sample IDs so that stage 10's
    strict tree/manifest match validation passes.
    """
    grouped: Dict[str, List[SyntheticSample]] = {lineage: [] for lineage in LINEAGES}
    for sample in samples:
        grouped[sample.lineage].append(sample)

    counter = {"n": 0}

    def next_label() -> str:
        counter["n"] += 1
        return f"n{counter['n']}"

    def clade(members: List[SyntheticSample]) -> str:
        """Build a balanced clade. A single member is a leaf."""
        if len(members) == 1:
            return f"{members[0].sample_id}:0.05"
        midpoint = len(members) // 2
        return f"({clade(members[:midpoint])},{clade(members[midpoint:])}){next_label()}:0.05"

    parts = [clade(sorted(grouped[lineage], key=lambda s: s.index)) for lineage in LINEAGES]
    newick = f"({','.join(parts)}){next_label()}:0.0;"
    path = out_dir / "phylogeny" / "tree.nwk"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(newick + "\n", encoding="utf-8")
    return {"tree": path}


def write_tree_metadata(out_dir: Path, samples: Sequence[SyntheticSample]) -> Dict[str, Path]:
    """Write the tree/manifest crosswalk."""
    rows = [
        {
            "sample_id": s.sample_id,
            "tree_tip_label": s.sample_id,
            "lineage_label": s.lineage,
            "st": s.st,
            "source": "synthetic_generator",
        }
        for s in samples
    ]
    path = out_dir / "phylogeny" / "tree_metadata.tsv"
    write_tsv(
        path,
        rows,
        ["sample_id", "tree_tip_label", "lineage_label", "st", "source"],
        header_comment=_banner(),
    )
    return {"tree_metadata": path}


def write_alignments(out_dir: Path, samples: Sequence[SyntheticSample], n_sites: int = 240) -> Dict[str, Path]:
    """Write miniature core-gene and core-SNP alignments.

    The SNP alignment contains **only variable sites**, so the
    ``snp_density_vs_core`` figure reported by stage 10 is meaningful
    (a handful of SNPs over the core length) rather than trivially 1.0.
    """
    rng = random.Random(DEFAULT_SEED + 11)
    core_path = out_dir / "phylogeny" / "core_alignment.fasta"
    records = []
    for sample in samples:
        seq = "".join(rng.choice("ACGT") for _ in range(n_sites))
        records.append((sample.sample_id, seq))
    write_fasta(core_path, records)

    # Variable sites: every 50th position of the core alignment.
    variable_positions = list(range(0, n_sites, 50))
    snp_path = out_dir / "phylogeny" / "core_snp_alignment.fasta"
    snp_records = []
    for index, sample in enumerate(samples):
        chars = [
            "ACGT"[(index + position // 50) % 4] for position in variable_positions
        ]
        snp_records.append((sample.sample_id, "".join(chars)))
    write_fasta(snp_path, snp_records)
    return {"core_alignment": core_path, "core_snp_alignment": snp_path}


def write_gwas_features(out_dir: Path, samples: Sequence[SyntheticSample]) -> Dict[str, Path]:
    """Write the stage 12 feature matrix.

    Columns are the binary features a GWAS engine would consume. The latent
    ``resistance_module`` is deliberately recoverable, so a stage-12 test can
    assert that the association machinery finds *something* and that
    multiple-testing correction is applied.
    """
    rows = []
    for sample in samples:
        rows.append(
            {
                "sample_id": sample.sample_id,
                "gene__oprD_absent": "1" if sample.oprd_status == "absent" else "0",
                "gene__oprD_LoF": "1" if sample.oprd_status == "disrupted" else "0",
                "gene__acquired_determinant": "1" if sample.acquired_determinant else "0",
                "snp__TESTPOS_00100": "1" if sample.resistance_module else "0",
                "snp__TESTPOS_00200": "1" if sample.index % 2 == 0 else "0",
                "unitig__utg000001": "1" if sample.lineage == "LINEAGE_A" else "0",
                "unitig__utg000002": "1" if sample.lineage == "LINEAGE_B" else "0",
                "kmer__k0000001": "1" if sample.lineage == "LINEAGE_C" else "0",
            }
        )
    path = out_dir / "intermediate" / "gwas" / "gwas_features.tsv"
    write_tsv(
        path,
        rows,
        ["sample_id"]
        + [
            "gene__oprD_absent",
            "gene__oprD_LoF",
            "gene__acquired_determinant",
            "snp__TESTPOS_00100",
            "snp__TESTPOS_00200",
            "unitig__utg000001",
            "unitig__utg000002",
            "kmer__k0000001",
        ],
        header_comment=_banner(
            ["feature names are synthetic labels; no biological feature is implied"]
        ),
    )
    return {"gwas_features": path}


def write_pangenome(out_dir: Path, samples: Sequence[SyntheticSample]) -> Dict[str, Path]:
    """Write pan-genome style outputs derived from the synthetic cohort."""
    from ..io.tsv import read_tsv

    annotation_dir = out_dir / "intermediate" / "annotation"
    gene_sample: Dict[str, set] = {}
    for sample in samples:
        path = annotation_dir / f"{sample.sample_id}.annotation.tsv"
        for row in read_tsv(path, required_columns=("gene_name",)):
            name = row.get("gene_name")
            if name:
                gene_sample.setdefault(name, set()).add(sample.sample_id)

    total = len(samples)
    # Long format, matching stages/pangenome.py's GPA_COLUMNS. This fixture is a
    # consumer of that contract, so it moves in the same commit as the contract -
    # leaving it on the old shape is how a reader and a writer end up
    # disagreeing while every test that reads only this file still passes.
    presence_path = out_dir / "intermediate" / "pangenome" / "gene_presence_absence.tsv"
    write_tsv(
        presence_path,
        [
            {"gene": gene, "sample_id": sample.sample_id, "present": int(gene in ids)}
            for gene, ids in sorted(gene_sample.items())
            for sample in samples
        ],
        list(GPA_COLUMNS),
        header_comment=_banner(),
    )

    core = [g for g, ids in gene_sample.items() if len(ids) == total]
    accessory = [g for g, ids in gene_sample.items() if len(ids) < total]

    core_path = out_dir / "intermediate" / "pangenome" / "core_genes.tsv"
    write_tsv(core_path, [{"gene": g} for g in sorted(core)], ["gene"], header_comment=_banner())
    acc_path = out_dir / "intermediate" / "pangenome" / "accessory_genes.tsv"
    write_tsv(
        acc_path,
        [{"gene": g, "n_samples": len(gene_sample[g])} for g in sorted(accessory)],
        ["gene", "n_samples"],
        header_comment=_banner(),
    )
    summary_path = out_dir / "intermediate" / "pangenome" / "pangenome_summary.tsv"
    write_tsv(
        summary_path,
        [
            {
                "metric": "n_samples",
                "value": str(total),
            },
            {"metric": "n_genes_total", "value": str(len(gene_sample))},
            {"metric": "n_core_genes", "value": str(len(core))},
            {"metric": "n_accessory_genes", "value": str(len(accessory))},
            {
                "metric": "core_fraction",
                "value": f"{(len(core)/len(gene_sample)):.4f}" if gene_sample else "0",
            },
        ],
        ["metric", "value"],
        header_comment=_banner(),
    )
    return {
        "gene_presence_absence": presence_path,
        "core_genes": core_path,
        "accessory_genes": acc_path,
        "pangenome_summary": summary_path,
    }


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------


def generate(
    out_dir: Path,
    n_samples: int = DEFAULT_N_SAMPLES,
    seed: int = DEFAULT_SEED,
    antibiotic: str = "imipenem",
    write_genome_files: bool = True,
    regulators_table: Optional[Path] = None,
) -> Dict[str, Path]:
    """Generate the full synthetic fixture set under ``out_dir``.

    Refuses to write into a directory named ``data`` so synthetic fixtures
    can never be mixed with the real data tree.
    """
    out_dir = Path(out_dir).resolve()
    if out_dir.name == "data":
        raise PipelineError(
            "Refusing to write synthetic test data into a 'data' directory",
            out_dir=str(out_dir),
            hint="Synthetic fixtures belong in test_data/",
        )

    if regulators_table is None:
        regulators_table = Path.cwd() / "config" / "regulators.tsv"
    choices = variant_choices_from_table(regulators_table)
    disruptive_choices = disruptive_variant_choices_from_table(regulators_table)

    for sub in (
        "genomes",
        "metadata",
        "phenotype",
        "phylogeny",
        "intermediate/annotation",
        "intermediate/mlst",
        "intermediate/amr",
        "intermediate/regulators",
        "intermediate/structural_variants",
        "intermediate/virulence",
        "intermediate/gwas",
        "intermediate/pangenome",
    ):
        (out_dir / sub).mkdir(parents=True, exist_ok=True)

    samples = build_samples(
        n_samples=n_samples, seed=seed, variant_choices=choices
    )

    written: Dict[str, Path] = {}
    written.update(write_metadata(out_dir, samples))
    if write_genome_files:
        written.update(write_genomes(out_dir, samples, seed=seed))
    written.update(write_phenotype(out_dir, samples, antibiotic=antibiotic))
    written.update(write_annotation(out_dir, samples))
    written.update(write_mlst(out_dir, samples))
    written.update(write_amr(out_dir, samples, antibiotic=antibiotic))
    written.update(
        write_regulator_variants(
            out_dir,
            samples,
            oprd_disruptive=disruptive_choices.get(
                "oprD", _DEFAULT_DISRUPTIVE
            ),
        )
    )
    written.update(write_structural_variants(out_dir, samples))
    written.update(write_virulence(out_dir, samples))
    written.update(write_tree(out_dir, samples))
    written.update(write_tree_metadata(out_dir, samples))
    written.update(write_alignments(out_dir, samples))
    written.update(write_gwas_features(out_dir, samples))
    written.update(write_pangenome(out_dir, samples))

    LOGGER.info(
        "%s Generated %d synthetic samples in %s",
        SYNTHETIC_BANNER,
        len(samples),
        out_dir,
    )
    return written
