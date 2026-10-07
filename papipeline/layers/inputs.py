"""Reading the source tables the layer builders consume.

Five inputs, four of them committed fixtures in TEST mode:

* ``amr/amr_determinants.tsv`` - AMRFinderPlus calls: L1 acquired genes, L2
  PDC alleles and their residue positions, L5 target-gene variants;
* ``regulators/regulator_variants.tsv`` - the chromosomal screen: L2/L3/L4
  variant classes, loaded through
  :func:`papipeline.stages.regulators.load_regulator_variants` so promoter
  gating and the provenance sidecar behave exactly as they do everywhere else;
* ``structural_variants/structural_variants.tsv`` - confirmed versus candidate
  structural calls at oprD (IS insertion, large deletion);
* ``oprd_locus/structural_calls.tsv`` - the REAL-mode tblastn verdicts
  (``absent`` / ``disrupted`` / ``intact`` / ``not_assessed``), read only when
  present: ``contracts.py`` declares no oprD table, so this one is optional by
  construction and TEST fixtures do not carry it;
* the raw regulator and AMR rows, kept alongside the typed records because
  ``partial_call_flag`` needs the optional ``call_state`` column, which the
  typed dataclasses do not model.

Sample identifiers are checked against the manifest on every table, and a row
naming a sample the manifest does not contain is a hard failure naming the
sample (AGENTS.md rule 5). A manifest sample *without* a row is not a failure:
these tables are presence-dependent, so no row means the feature was not seen,
which is a finding rather than a missing record.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

from ..config.loader import PipelineConfig
from ..errors import SampleIdError
from ..io.tsv import read_tsv
from ..logging_utils import get_logger
from ..manifest import SampleManifest
from ..models import RegulatorVariant
from ..stages import regulators as reg_stage
from ..stages.gwas_features import AMR_REQUIRED

LOGGER = get_logger("layers")

#: Where each input lives under the stage-input root.
AMR_RELPATH = Path("amr") / "amr_determinants.tsv"
REGULATOR_RELPATH = Path("regulators") / "regulator_variants.tsv"
SV_RELPATH = Path("structural_variants") / "structural_variants.tsv"
OPRD_STRUCTURAL_RELPATH = Path("oprd_locus") / "structural_calls.tsv"

#: Columns this package needs beyond what the owning stage requires.
SV_REQUIRED: Tuple[str, ...] = (
    "sample_id",
    "variant_id",
    "variant_type",
    "call_status",
    "affected_gene",
)
OPRD_STRUCTURAL_REQUIRED: Tuple[str, ...] = ("sample_id", "structural_verdict")

#: The optional per-row call state. Not every producer emits it; when no
#: source table does, ``partial_call_flag`` is written as missing rather than
#: as a fabricated 0 (``docs/scientific_rules.md``: an unmeasurable value is
#: missing, not zero).
CALL_STATE_COLUMN = "call_state"

#: The states the ticket names, cross-checked against
#: ``pilot.pdc_fields.KNOWN_CALL_STATES`` by the tests rather than trusted
#: here.
PARTIAL_CALL_STATES = frozenset({"PARTIAL", "MISTRANSLATION", "HMM"})


@dataclass
class Sources:
    """Everything the five builders read, filtered to one antibiotic."""

    antibiotic: str
    sample_ids: Tuple[str, ...]
    amr_rows: Tuple[Mapping[str, Any], ...] = ()
    regulator_records: Tuple[RegulatorVariant, ...] = ()
    regulator_rows: Tuple[Mapping[str, Any], ...] = ()
    sv_rows: Tuple[Mapping[str, Any], ...] = ()
    #: ``None`` when the optional table was not produced (TEST, and REAL
    #: before stage 6 runs its structural screen). Not the same as empty: an
    #: absent table is "not assessed", an empty one would be a contract fault
    #: and is read strictly.
    oprd_structural_rows: Optional[Tuple[Mapping[str, Any], ...]] = None

    def rows_for_state_check(self) -> Sequence[Mapping[str, Any]]:
        """Every raw row that could carry a ``call_state``."""
        rows: list[Mapping[str, Any]] = list(self.amr_rows)
        rows.extend(self.regulator_rows)
        rows.extend(self.sv_rows)
        if self.oprd_structural_rows:
            rows.extend(self.oprd_structural_rows)
        return rows


def _cell(row: Mapping[str, Any], column: str) -> str:
    value = row.get(column)
    if value is None:
        return ""
    return str(value).strip()


def _require_within_cohort(
    expected: Sequence[str], found: Sequence[str], path: Path, source: str
) -> None:
    """Refuse a row naming a sample the manifest does not contain.

    Presence-dependent semantics, as documented in
    :func:`papipeline.stages.gwas_features._check_cohort`: a manifest sample
    with no row is carrying none of the features, which is a finding.
    """
    unexpected = sorted(set(found) - set(expected))
    if not unexpected:
        return
    raise SampleIdError(
        f"{source} does not describe this run's cohort. "
        f"unexpected_in_table={unexpected}. Sample-ID mismatches fail loudly "
        "and are never dropped: a feature matrix built over a different "
        "cohort is wrong in ways no later stage can detect. The usual cause "
        "is a table left on disk by an earlier run over a different manifest.",
        path=str(path),
        source=source,
        unexpected=unexpected,
    )


def load_sources(
    config: PipelineConfig,
    manifest: SampleManifest,
    input_root: Path,
    antibiotic: str,
) -> Sources:
    """Read every layer source under ``input_root`` for one antibiotic.

    Args:
        config: Loaded pipeline configuration.
        manifest: The cohort every table must belong to.
        input_root: The stage-input root (``tool_output_root`` for the run,
            ``test_data/intermediate`` in TEST).
        antibiotic: The antibiotic the AMR rows must be filtered to by the
            builders; recorded here so the bundle carries it.

    Raises:
        DataContractError: A table is missing, malformed or empty, or an
            optional table is present but empty.
        SampleIdError: A table names a sample the manifest does not contain.
    """
    root = Path(input_root)
    sample_ids = tuple(manifest.sample_ids)

    amr_path = root / AMR_RELPATH
    amr_rows = tuple(read_tsv(amr_path, required_columns=AMR_REQUIRED))
    _require_within_cohort(
        sample_ids, [str(row["sample_id"]) for row in amr_rows], amr_path,
        "the AMR determinant table",
    )

    regulator_path = root / REGULATOR_RELPATH
    regulator_rows = tuple(read_tsv(regulator_path, required_columns=reg_stage.REQUIRED))
    _require_within_cohort(
        sample_ids,
        [str(row["sample_id"]) for row in regulator_rows],
        regulator_path,
        "the regulator variant table",
    )
    # Typed records for the builders: this applies promoter gating, the
    # knowledge-table checks and the empty-table provenance rule exactly as
    # the stage that owns the screen applies them.
    regulator_records = tuple(
        reg_stage.load_regulator_variants(config, regulator_path, root)
    )

    sv_path = root / SV_RELPATH
    sv_rows = tuple(read_tsv(sv_path, required_columns=SV_REQUIRED))
    _require_within_cohort(
        sample_ids, [str(row["sample_id"]) for row in sv_rows], sv_path,
        "the structural variant table",
    )

    structural_path = root / OPRD_STRUCTURAL_RELPATH
    oprd_structural_rows: Optional[Tuple[Mapping[str, Any], ...]] = None
    if structural_path.exists():
        rows = tuple(
            read_tsv(structural_path, required_columns=OPRD_STRUCTURAL_REQUIRED)
        )
        _require_within_cohort(
            sample_ids,
            [str(row["sample_id"]) for row in rows],
            structural_path,
            "the oprD structural call table",
        )
        oprd_structural_rows = rows

    sources = Sources(
        antibiotic=antibiotic,
        sample_ids=sample_ids,
        amr_rows=amr_rows,
        regulator_records=regulator_records,
        regulator_rows=regulator_rows,
        sv_rows=sv_rows,
        oprd_structural_rows=oprd_structural_rows,
    )
    LOGGER.info(
        "layers: loaded %d AMR rows, %d regulator variants (%d raw rows), "
        "%d structural variants, %s oprD structural calls for %s",
        len(amr_rows),
        len(regulator_records),
        len(regulator_rows),
        len(sv_rows),
        "an optional" if oprd_structural_rows is not None else "no",
        antibiotic,
    )
    return sources


def partial_call_flags(sources: Sources) -> Dict[str, Optional[int]]:
    """Per-isolate ``partial_call_flag``: 1, 0, or ``None`` for untracked.

    ``PARTIAL``, ``MISTRANSLATION`` and ``HMM`` set the flag. If *no* source
    table in the run carries a ``call_state`` column at all, every value is
    ``None`` - the writer emits ``.`` - because "we never looked" and "we
    looked and saw no partial call" are different findings and only the second
    may be written as 0.
    """
    rows = sources.rows_for_state_check()
    tracked = any(_cell(row, CALL_STATE_COLUMN) for row in rows)
    if not tracked:
        LOGGER.warning(
            "layers: no source table carries a %s column; partial_call_flag is "
            "written as missing for all %d isolates rather than as 0",
            CALL_STATE_COLUMN,
            len(sources.sample_ids),
        )
        return {sample_id: None for sample_id in sources.sample_ids}

    flags: Dict[str, Optional[int]] = {
        sample_id: 0 for sample_id in sources.sample_ids
    }
    for row in rows:
        state = _cell(row, CALL_STATE_COLUMN).upper()
        if state in PARTIAL_CALL_STATES:
            sample_id = _cell(row, "sample_id")
            if sample_id in flags:
                flags[sample_id] = 1
    return flags
