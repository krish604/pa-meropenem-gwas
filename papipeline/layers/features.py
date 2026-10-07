"""The two containers every layer builder and writer shares.

A :class:`Feature` is one column of a layer matrix plus the provenance that
makes it auditable: which layer it belongs to, which source tables it was read
from, and the exact rule that turns those rows into a value. The ticket asks
for a ``feature_dictionary.tsv`` of *feature, layer, source tokens, rule*, so
those four things travel with the feature from the builder onward rather than
being reconstructed by the writer.

A :class:`LayerBundle` holds one value per feature per isolate. Values are
integers: ``0``/``1`` for binary features, and a count for the one feature the
ticket defines as a count (``L5_QRDR_count``).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple

#: How a feature's values should be read. ``binary`` means presence/absence;
#: ``count`` means a non-negative integer that is *not* a presence flag, which
#: is why the drop rule counts a carrier as ``value > 0`` for both kinds.
BINARY = "binary"
COUNT = "count"


@dataclass(frozen=True)
class Feature:
    """One column of a layer matrix, with its provenance."""

    name: str
    layer: int
    value_kind: str = BINARY
    source_tokens: Tuple[str, ...] = ()
    rule: str = ""


@dataclass
class LayerBundle:
    """All five layers' features and values over one cohort.

    ``values[feature][sample_id]`` is the cell. Every feature must have an
    entry for every sample, and ``sample_ids`` must match the manifest - the
    writers refuse anything else rather than emitting tables that disagree
    about the cohort.
    """

    sample_ids: Tuple[str, ...]
    features: List[Feature] = field(default_factory=list)
    values: Dict[str, Dict[str, int]] = field(default_factory=dict)

    def layer_features(self, layer: int) -> List[Feature]:
        """The bundle's features for one layer, in bundle order."""
        return [feature for feature in self.features if feature.layer == layer]

    def value(self, feature_name: str, sample_id: str) -> int:
        return self.values[feature_name][sample_id]
