"""Layer-prefixed feature encoding for the meropenem GWAS build.

Five layers of features over one cohort, each column named with its layer
prefix so provenance is visible in the name:

======  ==================================================================
Layer   What it encodes
======  ==================================================================
``L1``  Acquired determinants: ``any_MBL`` (NDM/VIM/IMP) and the KPC, GES
        and OXA groups.
``L2``  PDC: allele columns, residue-position flags, ampD/ampR/dacB loss and
        missense, and the ``PDC_high`` composite (allele **and** a regulator
        loss).
``L3``  oprD: ``absent``, ``LoF_tier1`` (stop, frameshift, IS insertion,
        large deletion), ``LoF_tier2`` (probable loss), the ``off_any``
        composite, and every other call as an individual feature - a
        missense never enters the composite.
``L4``  Efflux regulators: loss and missense for mexR, nalC, nalD, nfxB and
        mexZ, plus pump-high composites whose groups are read from
        ``config/mechanisms.tsv``.
``L5``  Targets: gyrA, parC, ftsI, ``any_QRDR`` and ``QRDR_count``, using the
        PAO1 numbering AMRFinderPlus emits.
======  ==================================================================

Outputs: one TSV per layer, ``all_layers.tsv``, ``<antibiotic>_feature_matrix.tsv``
(phenotype merged) and ``feature_dictionary.tsv``. A feature present in fewer
than ``layers.min_carriers`` isolates or in more than
``layers.max_prevalence`` of them is dropped everywhere, and the drop list is
logged and bannered. Every output carries the manifest's isolates exactly.

The ``unitig-caller`` adapter lives here too, but is a separate entry point
(:func:`run_unitig_screen`): it stubs in TEST/STUB and refuses in REAL while
its flags are unverified.

Test mode only unless a caller says otherwise: nothing in this package reads
``data/``, ``db/`` or ``PDC_essential.tsv``.
"""

from __future__ import annotations

from .builders import (
    build_layer1,
    build_layer2,
    build_layer3,
    build_layer4,
    build_layer5,
    build_layers,
    pump_groups,
)
from .encode import encode_layers
from .features import Feature, LayerBundle
from .inputs import (
    CALL_STATE_COLUMN,
    PARTIAL_CALL_STATES,
    Sources,
    load_sources,
    partial_call_flags,
)
from .settings import (
    DEFAULT_MAX_PREVALENCE,
    DEFAULT_MIN_CARRIERS,
    MAX_PREVALENCE_KEY,
    MIN_CARRIERS_KEY,
    LayerSettings,
    layer_settings,
)
from .unitigs import (
    FLAGS_VERIFIED,
    UNITIG_TABLE_FILENAME,
    UNITIG_TOOL,
    UNVERIFIED_REASON,
    VERIFIED_FLAGS,
    UnverifiedFlagsError,
    run_unitig_screen,
    unitig_tool_status,
)
from .writers import (
    ALL_LAYERS_FILENAME,
    DICTIONARY_FILENAME,
    FLAG_COLUMN,
    LAYER_FILENAMES,
    PHENOTYPE_COLUMN,
    DroppedFeature,
    LayerResult,
    apply_drop,
    feature_matrix_filename,
    write_layers,
)

__all__ = [
    "ALL_LAYERS_FILENAME",
    "CALL_STATE_COLUMN",
    "DEFAULT_MAX_PREVALENCE",
    "DEFAULT_MIN_CARRIERS",
    "DICTIONARY_FILENAME",
    "DroppedFeature",
    "FLAGS_VERIFIED",
    "Feature",
    "FLAG_COLUMN",
    "LAYER_FILENAMES",
    "MAX_PREVALENCE_KEY",
    "MIN_CARRIERS_KEY",
    "LayerBundle",
    "LayerResult",
    "LayerSettings",
    "PARTIAL_CALL_STATES",
    "PHENOTYPE_COLUMN",
    "Sources",
    "UNVERIFIED_REASON",
    "UNITIG_TABLE_FILENAME",
    "UNITIG_TOOL",
    "VERIFIED_FLAGS",
    "UnverifiedFlagsError",
    "apply_drop",
    "build_layer1",
    "build_layer2",
    "build_layer3",
    "build_layer4",
    "build_layer5",
    "build_layers",
    "encode_layers",
    "feature_matrix_filename",
    "layer_settings",
    "load_sources",
    "partial_call_flags",
    "pump_groups",
    "run_unitig_screen",
    "unitig_tool_status",
    "write_layers",
]
