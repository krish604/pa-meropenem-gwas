"""The names the dashboard reuses rather than redeclares.

Every constant below is **imported**, never typed. A path or a stage list
written into dashboard code would be a second answer to a question
`contracts.py` or `run.py` already answers, which is the drift class those
modules' own docstrings say they exist to prevent.

Import-time side effects: none. All of these modules were verified to import
silently, with the working tree unchanged afterwards (DESIGN §0), so nothing
here needs a `PAPipeline` module for its side effects.
"""

from __future__ import annotations

from typing import Any, Dict, Mapping

from papipeline.execution.contracts import (
    INTERNAL_TABLES,
    STAGE_TABLES,
    required_columns,
)
from papipeline.io.tsv import MISSING_SENTINELS, iter_tsv_dicts, read_tsv
from papipeline.run import (
    EXECUTION_ORDER,
    PREREQUISITES,
    REAL_REFUSING_STAGES,
    STAGE_ORDER,
    UNBUILT_STAGES,
)
from papipeline.stages.report_tables import (
    INTERNAL_TABLE_KEYS,
    NOT_PRODUCED,
    SNAKEFILE_ONLY_TABLES,
    STAGE_TABLE_KEYS,
    TableRead,
    read_json_sidecar,
    read_table,
)
from papipeline.stages.reporting import (
    POWER_FLAG_MEANING,
    POWER_FLAG_TEMPLATE,
    R16_FLAG,
    power_flag,
    power_flag_line,
)
from papipeline.stages.similarity import (
    DISTANCE_UNIT,
    NOT_A_SNP_COUNT,
    UNITS_SIDECAR_SUFFIX,
    patristic_distances,
    units_sidecar_path,
)

#: spec.md D1's row for each code stage name. Four differ (§1 of DESIGN).
#:
#: spec.md's own numbering is fifteen rows carrying sixteen names:
#: `cohort_variants` is row `6a`. So the spec number is a string here and may
#: be `"6a"`.
SPEC_NAMES: Mapping[str, str] = {
    "validation": "validate",
    "annotation": "annotate",
    "mlst": "mlst",
    "amr": "amr",
    "virulence": "virulence",
    "variants": "variants",
    "cohort_variants": "cohort_variants",
    "pangenome": "pangenome",
    "recombination": "recombination",
    "phylogeny": "phylogeny",
    "similarity": "similarity",
    "phenotype": "phenotype",
    "gwas": "gwas",
    "convergence": "convergence",
    "cooccurrence": "combination",
    "reporting": "report",
}

#: spec.md D1's row number per stage. `cohort_variants` is `6a`.
SPEC_NUMBERS: Mapping[str, str] = {
    "validation": "1",
    "annotation": "2",
    "mlst": "3",
    "amr": "4",
    "virulence": "5",
    "variants": "6",
    "cohort_variants": "6a",
    "pangenome": "7",
    "recombination": "8",
    "phylogeny": "9",
    "similarity": "10",
    "phenotype": "11",
    "gwas": "12",
    "convergence": "13",
    "cooccurrence": "14",
    "reporting": "15",
}


def stage_table_path(stage_dir, stage: str):
    """Where a stage's principal table lives, per `contracts.table_path`."""
    from papipeline.execution.contracts import table_path

    return table_path(stage_dir, stage)


def config_snapshot_facts(config: Any) -> Dict[str, Any]:
    """The configured minimum sample size, for the D3 power flag.

    Read from `config.gwas.min_samples_per_group` rather than typed, because it
    is a scientific setting in `config/science.yaml` and a dashboard that
    hard-codes it would report a threshold the pipeline did not use.
    """
    gwas = getattr(config, "gwas", None)
    value = getattr(gwas, "min_samples_per_group", None)
    return {"min_samples_per_group": int(value) if value is not None else None}


__all__ = [
    "DISTANCE_UNIT",
    "EXECUTION_ORDER",
    "INTERNAL_TABLE_KEYS",
    "INTERNAL_TABLES",
    "MISSING_SENTINELS",
    "NOT_A_SNP_COUNT",
    "NOT_PRODUCED",
    "POWER_FLAG_MEANING",
    "POWER_FLAG_TEMPLATE",
    "PREREQUISITES",
    "R16_FLAG",
    "REAL_REFUSING_STAGES",
    "SNAKEFILE_ONLY_TABLES",
    "SPEC_NAMES",
    "SPEC_NUMBERS",
    "STAGE_ORDER",
    "STAGE_TABLE_KEYS",
    "STAGE_TABLES",
    "TableRead",
    "UNITS_SIDECAR_SUFFIX",
    "UNBUILT_STAGES",
    "config_snapshot_facts",
    "iter_tsv_dicts",
    "patristic_distances",
    "power_flag",
    "power_flag_line",
    "read_json_sidecar",
    "read_table",
    "read_tsv",
    "required_columns",
    "stage_table_path",
    "units_sidecar_path",
]