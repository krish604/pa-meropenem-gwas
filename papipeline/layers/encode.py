"""The one entry point that reads sources and writes all eight layer outputs.

Wiring, deliberately in one place: choose the antibiotic the way stage 12 does
(so the feature table and the analysis cannot disagree about which drug is
under study), read the five sources, build the layers, resolve the drop
thresholds from configuration, merge the phenotype, and write.

The unitig stub is *not* produced here. It is a separate screen with its own
refusal semantics (see :mod:`papipeline.layers.unitigs`), and folding a tool
that REAL mode refuses to invoke into an entry point that otherwise succeeds
would hide that refusal.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

from ..config.loader import PipelineConfig
from ..logging_utils import get_logger
from ..manifest import SampleManifest
from ..models import RunMode
from ..stages.gwas_features import target_antibiotic
from ..stages.phenotype import load_phenotype
from .builders import build_layers
from .inputs import load_sources, partial_call_flags
from .settings import layer_settings
from .writers import LayerResult, write_layers

LOGGER = get_logger("layers")


def encode_layers(
    config: PipelineConfig,
    manifest: SampleManifest,
    *,
    mode: RunMode,
    input_root: Path,
    output_dir: Path,
    phenotype_dir: Path,
) -> LayerResult:
    """Build and write the layer feature tables for one cohort.

    Args:
        config: Loaded pipeline configuration.
        manifest: The cohort; every output carries exactly these isolates.
        mode: STUB/TEST/REAL, recorded in each file's provenance banner.
        input_root: Where the stage-input tables are read from
            (``tool_output_root`` in a run, ``test_data/intermediate`` in
            TEST).
        output_dir: Where the eight outputs are written.
        phenotype_dir: Where the phenotype table for the target antibiotic
            lives; the matrix's ``phenotype`` column is merged from it.

    Returns:
        The paths written, the surviving features, and the drop list.
    """
    antibiotic = target_antibiotic(config)
    sources = load_sources(config, manifest, Path(input_root), antibiotic)
    bundle = build_layers(config, sources, tuple(manifest.sample_ids))
    settings = layer_settings(config)

    calls = load_phenotype(
        config,
        Path(phenotype_dir),
        antibiotic,
        sample_ids=manifest.sample_ids,
    )
    cohort = set(manifest.sample_ids)
    phenotype_by_sample: Dict[str, Optional[str]] = {}
    for call in calls:
        if call.sample_id not in cohort:
            # A phenotype row with no assembly is excluded, not joined (spec
            # D3); the same row must not appear in a matrix over assemblies.
            continue
        phenotype_by_sample[call.sample_id] = call.phenotype.value

    result = write_layers(
        bundle,
        manifest=manifest,
        phenotype_by_sample=phenotype_by_sample,
        partial_flags=partial_call_flags(sources),
        antibiotic=antibiotic,
        settings=settings,
        mode=mode,
        output_dir=Path(output_dir),
    )
    LOGGER.info(
        "layers: %s matrix written as %s with %d features (%d phenotype "
        "records merged)",
        antibiotic,
        result.paths["matrix"].name,
        len(result.features),
        len(phenotype_by_sample),
    )
    return result
