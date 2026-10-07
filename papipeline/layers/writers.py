"""Writing the layer outputs: eight tables, one isolate set, a logged drop list.

Outputs, all under one directory (``<stage-input root>/layers/`` in a run,
``tmp_path`` in a test):

===========================  =========================================
``layer1_acquired.tsv``      ``L1_`` features only
``layer2_pdc.tsv``           ``L2_`` features only
``layer3_oprD.tsv``          ``L3_`` features only
``layer4_efflux.tsv``        ``L4_`` features only
``layer5_targets.tsv``       ``L5_`` features only
``all_layers.tsv``           every surviving feature plus the covariate
``<antibiotic>_feature_matrix.tsv``  features + ``partial_call_flag`` + phenotype
``feature_dictionary.tsv``   feature, layer, source_tokens, rule
===========================  =========================================

Two rules run here rather than in the builders:

* **the drop rule** (:func:`apply_drop`): a feature present in fewer than
  ``layers.min_carriers`` isolates or in more than
  ``layers.max_prevalence`` of them is dropped from *every* output, so the
  layer files, the union and the matrix always describe the same feature set.
  The drop list is logged and recorded in the provenance banner of
  ``all_layers.tsv``.
* **the isolate-set guard**: every per-isolate output must carry exactly the
  manifest's samples, one row each. A bundle that disagrees is refused with a
  message naming the output file, because eight tables that quietly disagree
  about the cohort are wrong in ways no later stage can detect.

The matrix is drug-parameterised by name: it is
``meropenem_feature_matrix.tsv`` exactly when ``project.primary_antibiotic`` is
``meropenem``. A file named for one drug while carrying another drug's
phenotype would be a lie, so the name follows the configuration instead of the
ticket's example.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from ..errors import DataContractError
from ..io.tsv import write_tsv
from ..logging_utils import get_logger
from ..manifest import SampleManifest
from ..models import RunMode
from .features import Feature, LayerBundle
from .settings import (
    MAX_PREVALENCE_KEY,
    MIN_CARRIERS_KEY,
    LayerSettings,
)

LOGGER = get_logger("layers")

#: One file per layer, plus the three cross-cutting outputs.
LAYER_FILENAMES: Dict[int, str] = {
    1: "layer1_acquired.tsv",
    2: "layer2_pdc.tsv",
    3: "layer3_oprD.tsv",
    4: "layer4_efflux.tsv",
    5: "layer5_targets.tsv",
}
ALL_LAYERS_FILENAME = "all_layers.tsv"
DICTIONARY_FILENAME = "feature_dictionary.tsv"
FLAG_COLUMN = "partial_call_flag"
PHENOTYPE_COLUMN = "phenotype"

#: How the partial-call states are spelled in the banner, for a reader who
#: finds a ``1`` and wants to know which call set it.
FLAG_STATES = "PARTIAL,MISTRANSLATION,HMM"


def feature_matrix_filename(antibiotic: str) -> str:
    """``<antibiotic>_feature_matrix.tsv`` - the ticket's example name for the
    meropenem build appears once the project is configured for meropenem."""
    return f"{antibiotic}_feature_matrix.tsv"


@dataclass(frozen=True)
class DroppedFeature:
    """One feature the drop rule removed, and why."""

    name: str
    layer: int
    carriers: int
    reason: str


@dataclass(frozen=True)
class LayerResult:
    """What was written, and what was dropped on the way."""

    antibiotic: str
    output_dir: Path
    paths: Dict[str, Path]
    features: Tuple[Feature, ...]
    dropped: Tuple[DroppedFeature, ...]
    settings: LayerSettings


def _carriers(values: Mapping[str, int]) -> int:
    return sum(1 for value in values.values() if value > 0)


def apply_drop(
    bundle: LayerBundle, settings: LayerSettings
) -> Tuple[List[Feature], List[DroppedFeature]]:
    """Split the bundle's features into those to keep and those to drop.

    A carrier is any isolate with a value above 0, so the same rule applies to
    binary features and to the one count feature. ``partial_call_flag`` is not
    a feature and is never dropped; nor is the phenotype.
    """
    n_isolates = len(bundle.sample_ids)
    kept: List[Feature] = []
    dropped: List[DroppedFeature] = []
    for feature in bundle.features:
        column = bundle.values.get(feature.name)
        if column is None:
            raise DataContractError(
                f"feature {feature.name} has no value column",
                feature=feature.name,
            )
        carriers = _carriers(column)
        reason = settings.drop_reason(carriers, n_isolates)
        if reason is None:
            kept.append(feature)
            continue
        dropped.append(
            DroppedFeature(
                name=feature.name,
                layer=feature.layer,
                carriers=carriers,
                reason=reason,
            )
        )
        LOGGER.warning(
            "layers: dropped %s (layer %s): %s",
            feature.name,
            feature.layer,
            reason,
        )
    if dropped:
        LOGGER.warning(
            "layers: dropped %d of %d features (%s); full list in the %s banner",
            len(dropped),
            len(bundle.features),
            ", ".join(drop.name for drop in dropped),
            ALL_LAYERS_FILENAME,
        )
    return kept, dropped


def _assert_isolate_set(
    rows: Sequence[Mapping[str, Any]], expected: Sequence[str], path: Path
) -> None:
    """Refuse an output that does not carry exactly the manifest's isolates."""
    found = [row.get("sample_id") for row in rows]
    missing = sorted(set(expected) - set(found))
    unexpected = sorted(str(s) for s in set(found) - set(expected))
    duplicates = sorted(
        {sample_id for sample_id in found if found.count(sample_id) > 1}
    )
    if not (missing or unexpected or duplicates):
        return
    raise DataContractError(
        f"{path.name} does not carry this run's isolate set. "
        f"missing={missing} unexpected={unexpected} duplicates={duplicates}. "
        "Every layer output must cover exactly the manifest's isolates, one "
        "row each: a matrix that silently drops or invents a sample is wrong "
        "in ways no later stage can detect.",
        path=str(path),
        missing=missing,
        unexpected=unexpected,
        duplicates=duplicates,
    )


def _feature_rows(
    bundle: LayerBundle, features: Sequence[Feature], sample_ids: Sequence[str]
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for sample_id in sample_ids:
        row: Dict[str, Any] = {"sample_id": sample_id}
        for feature in features:
            column = bundle.values.get(feature.name)
            if column is None or sample_id not in column:
                raise DataContractError(
                    f"feature {feature.name} has no value for {sample_id}",
                    feature=feature.name,
                    sample_id=sample_id,
                )
            row[feature.name] = column[sample_id]
        rows.append(row)
    return rows


def _flag_banner(partial_flags: Mapping[str, Optional[int]]) -> str:
    tracked = any(flag is not None for flag in partial_flags.values())
    if tracked:
        return (
            f"{FLAG_COLUMN}: tracked (call_state states={FLAG_STATES}); "
            "1 for a PARTIAL, MISTRANSLATION or HMM call, else 0"
        )
    return (
        f"{FLAG_COLUMN}: not_tracked (no source table carries a call_state "
        "column; written as '.' rather than a fabricated 0)"
    )


def _drop_banner(dropped: Sequence[DroppedFeature]) -> str:
    if not dropped:
        return "drop_list: (none)"
    return "drop_list: " + "; ".join(
        f"{drop.name}(layer {drop.layer}, carriers={drop.carriers}, {drop.reason})"
        for drop in dropped
    )


def write_layers(
    bundle: LayerBundle,
    *,
    manifest: SampleManifest,
    phenotype_by_sample: Mapping[str, Optional[str]],
    partial_flags: Mapping[str, Optional[int]],
    antibiotic: str,
    settings: LayerSettings,
    mode: RunMode,
    output_dir: Path,
) -> LayerResult:
    """Apply the drop rule, then write all eight outputs.

    Raises:
        DataContractError: A feature has no value for an isolate, or a
            per-isolate output would not carry exactly the manifest's
            isolates. The message names the offending file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    expected = tuple(manifest.sample_ids)
    if tuple(bundle.sample_ids) != expected:
        # Not fatal by itself - the guard below reports the first output that
        # would disagree - but worth saying out loud, since it is the usual
        # cause of that refusal.
        LOGGER.warning(
            "layers: bundle sample_ids differ from the manifest; the "
            "per-output isolate-set guard will refuse the first mismatch"
        )

    kept, dropped = apply_drop(bundle, settings)
    by_layer: Dict[int, List[Feature]] = {
        layer: sorted(
            (feature for feature in kept if feature.layer == layer),
            key=lambda feature: feature.name,
        )
        for layer in LAYER_FILENAMES
    }
    ordered = [feature for layer in LAYER_FILENAMES for feature in by_layer[layer]]

    base_banner = [
        "source: papipeline/layers",
        f"mode: {mode.value}",
        f"antibiotic: {antibiotic}",
        f"drop_rule: carriers < {settings.min_carriers} ({MIN_CARRIERS_KEY}) or "
        f"carriers/n > {settings.max_prevalence} ({MAX_PREVALENCE_KEY})",
        f"features: {len(kept)} kept of {len(bundle.features)} built",
        _flag_banner(partial_flags),
    ]
    drop_banner = base_banner + [_drop_banner(dropped)]

    paths: Dict[str, Path] = {}

    def _write(
        key: str,
        filename: str,
        rows: List[Dict[str, Any]],
        columns: Sequence[str],
        banner: Sequence[str],
        *,
        check_isolates: bool = True,
    ) -> None:
        path = output_dir / filename
        if check_isolates:
            _assert_isolate_set(rows, expected, path)
        write_tsv(path, rows, list(columns), header_comment=list(banner))
        paths[key] = path

    for layer, filename in LAYER_FILENAMES.items():
        features = by_layer[layer]
        rows = _feature_rows(bundle, features, bundle.sample_ids)
        _write(
            f"layer{layer}",
            filename,
            rows,
            ["sample_id", *[feature.name for feature in features]],
            base_banner,
        )

    all_rows = _feature_rows(bundle, ordered, bundle.sample_ids)
    for row, sample_id in zip(all_rows, bundle.sample_ids):
        row[FLAG_COLUMN] = partial_flags.get(sample_id)
    _write(
        "all_layers",
        ALL_LAYERS_FILENAME,
        all_rows,
        ["sample_id", *[feature.name for feature in ordered], FLAG_COLUMN],
        drop_banner,
    )

    for row, sample_id in zip(all_rows, bundle.sample_ids):
        row[PHENOTYPE_COLUMN] = phenotype_by_sample.get(sample_id)
    _write(
        "matrix",
        feature_matrix_filename(antibiotic),
        all_rows,
        [
            "sample_id",
            *[feature.name for feature in ordered],
            FLAG_COLUMN,
            PHENOTYPE_COLUMN,
        ],
        drop_banner,
    )

    dictionary_rows: List[Dict[str, Any]] = [
        {
            "feature": feature.name,
            "layer": feature.layer,
            "source_tokens": feature.source_tokens,
            "rule": feature.rule,
        }
        for feature in ordered
    ]
    _write(
        "dictionary",
        DICTIONARY_FILENAME,
        dictionary_rows,
        ["feature", "layer", "source_tokens", "rule"],
        drop_banner,
        check_isolates=False,
    )

    LOGGER.info(
        "layers: wrote %d outputs to %s (%d features kept, %d dropped)",
        len(paths),
        output_dir,
        len(kept),
        len(dropped),
    )
    return LayerResult(
        antibiotic=antibiotic,
        output_dir=output_dir,
        paths=paths,
        features=tuple(kept),
        dropped=tuple(dropped),
        settings=settings,
    )
