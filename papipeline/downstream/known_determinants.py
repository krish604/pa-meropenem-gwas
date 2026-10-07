"""What "already known" means for the downstream scans.

``config/known_determinants.tsv`` is DERIVED from ``config/mechanisms.tsv``:
:func:`derive_known_determinants` regenerates the rows in memory and the
committed file has to equal them, so a mechanisms-table edit that is not
mirrored fails :mod:`tests.unit.test_downstream_config_tables` rather than
silently calling a now-known determinant novel.

Three kinds of label a scan column can resolve to:

``gene``
    a row of ``config/mechanisms.tsv`` - ``oprD``, ``mexR``, ``blaNDM-1``...
``composite``
    ``any_MBL``, union of the genes whose ``biological_role`` says
    metallo-beta-lactamase. Derived from roles, never invented: an MBL
    composite over one gene is not a composite, so it is emitted only when at
    least two such genes exist.
``layer``
    the two layer composites this build declares -
    :data:`LAYER_COMPOSITES`, spelled exactly as
    ``papipeline/layers/builders.py`` spells them (``L2_PDC_high``,
    ``L3_<OPRD_GENE>_off_any``). They are layer-encoder constructs, not
    mechanisms-table rows, which is why they are declared in code and do not
    appear in the derived TSV.

Matching accepts every spelling a scan emits: with or without the ``L1_``..
``L5_`` layer prefix, with or without the ``type__`` feature prefix, any case
(``L3_oprd_off_any`` and ``L3_oprD_off_any`` are the same column), the
bla-stripped family name the layer encoder uses (``L1_KPC`` for
``blaKPC-2``), and a gene stem (``oprD_burden``, ``oprD_absent`` -> ``oprD``).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Tuple

from ..errors import DataContractError
from ..io.tsv import read_tsv, write_tsv
from ..layers.builders import OPRD_GENE

#: ``kind`` column: a single gene from ``config/mechanisms.tsv``.
KIND_GENE = "gene"

#: ``kind`` column: a union of genes (``any_MBL``).
KIND_COMPOSITE = "composite"

#: The metallo-beta-lactamase composite. Its members are the mechanisms-table
#: rows whose ``biological_role`` contains "metallo-beta-lactamase".
MBL_COMPOSITE = "any_MBL"

#: The two layer composites this build declares as known endpoints. Keys are
#: the composite names WITHOUT the layer prefix; values are the genes they
#: concern, or an empty tuple when the composite concerns no single gene
#: (``PDC_high`` is an allele-plus-regulator rule).
#:
#: Exactly two, because the interaction-pair contract says endpoints must
#: resolve to a knowledge-table row or to one of these two. Adding a third is
#: a decision about what is already known, not a convenience.
LAYER_COMPOSITES: Mapping[str, Tuple[str, ...]] = {
    "PDC_high": (),
    "oprD_off_any": (OPRD_GENE,),
}

#: File name of the derived table, relative to ``config/``.
KNOWN_DETERMINANTS_FILENAME = "known_determinants.tsv"

#: Columns of the derived table, in write order.
COLUMNS: Tuple[str, ...] = (
    "determinant",
    "kind",
    "mechanism",
    "mechanism_class",
    "claim_ceiling",
    "members",
    "source",
)

#: Where the rows came from, recorded per row so a reader can trace one.
SOURCE = "config/mechanisms.tsv"

#: Role text that makes a gene a member of the MBL composite.
MBL_ROLE = "metallo-beta-lactamase"

_LAYER_PREFIX = re.compile(r"^L\d+_")

#: Trailing ``-2`` / ``-24/40`` version suffix, stripped when reducing a gene
#: to its family: ``blaOXA-24/40`` -> ``oxa``.
_VERSION_SUFFIX = re.compile(r"-\d+(/\d+)*$")


def normalise_label(feature: str) -> str:
    """Strip the layer prefix and the ``type__`` prefix off a column name.

    ``L3_oprd_off_any`` -> ``oprD_off_any``; ``gene__oprD_absent`` ->
    ``oprD_absent``; ``gene_presence_absence__blaNDM-1`` -> ``blaNDM-1``.
    The original spelling is preserved in what is left, so the caller decides
    whether to compare case-insensitively.
    """
    label = str(feature).strip()
    label = _LAYER_PREFIX.sub("", label)
    if "__" in label:
        label = label.split("__", 1)[1]
    return label.strip()


def _family_key(name: str) -> Optional[str]:
    """``blaKPC-2`` -> ``kpc``; ``None`` for anything that is not a bla gene."""
    key = name.casefold()
    if not key.startswith("bla"):
        return None
    key = _VERSION_SUFFIX.sub("", key[3:])
    return key or None


def known_determinants_path(config) -> Path:
    """Where the derived table lives, resolved from the repository root."""
    return Path(config.root) / "config" / KNOWN_DETERMINANTS_FILENAME


def derive_known_determinants(config) -> List[Dict[str, Optional[str]]]:
    """Regenerate the table from ``config/mechanisms.tsv``.

    One row per mechanisms-table gene, plus the MBL composite when at least
    two genes carry the metallo-beta-lactamase role. ``members`` is ``None``
    for a gene (not ``""``) because :func:`papipeline.io.tsv.read_tsv` reads
    the empty cell back as ``None`` and the two must compare equal.
    """
    rows: List[Dict[str, Optional[str]]] = []
    for gene, spec in config.mechanisms.items():
        rows.append(
            {
                "determinant": gene,
                "kind": KIND_GENE,
                "mechanism": spec.mechanism,
                "mechanism_class": spec.mechanism_class,
                "claim_ceiling": spec.claim_ceiling.value,
                "members": None,
                "source": SOURCE,
            }
        )

    mbl = [
        gene
        for gene, spec in config.mechanisms.items()
        if MBL_ROLE in spec.biological_role.casefold()
    ]
    if len(mbl) >= 2:
        first = config.mechanisms[mbl[0]]
        rows.append(
            {
                "determinant": MBL_COMPOSITE,
                "kind": KIND_COMPOSITE,
                "mechanism": first.mechanism,
                "mechanism_class": first.mechanism_class,
                "claim_ceiling": first.claim_ceiling.value,
                "members": ",".join(mbl),
                "source": SOURCE,
            }
        )
    return rows


def read_known_determinants(path: Path) -> List[Dict[str, Optional[str]]]:
    """Read the committed table. :func:`papipeline.io.tsv.read_tsv` refuses a
    missing column, a duplicate determinant or a short row on its own."""
    return read_tsv(
        path,
        required_columns=COLUMNS,
        unique_columns=("determinant",),
    )


def write_known_determinants(config, path: Optional[Path] = None) -> Path:
    """Write the derivation to disk - the regeneration command.

    ``python3 -c "from papipeline.config.loader import load_config as L; \
from papipeline.downstream.known_determinants import write_known_determinants as W; \
W(L('config/science.yaml', machine='laptop'))"``
    """
    destination = Path(path) if path is not None else known_determinants_path(config)
    return write_tsv(
        destination,
        derive_known_determinants(config),
        COLUMNS,
        header_comment=(
            "Derived from config/mechanisms.tsv - do not edit by hand.",
            "Regenerate with papipeline.downstream.known_determinants"
            ".write_known_determinants(config); tests/unit/"
            "test_downstream_config_tables.py fails if this file drifts.",
        ),
    )


@dataclass(frozen=True)
class KnownDeterminants:
    """The derived table plus the indexes :meth:`is_known` looks things up in."""

    rows: Tuple[Mapping[str, Optional[str]], ...]
    genes: Tuple[str, ...]
    composites: Mapping[str, Tuple[str, ...]]
    layers: Mapping[str, Tuple[str, ...]]
    families: Mapping[str, Tuple[str, ...]]

    def is_known(self, feature: str) -> bool:
        """Whether anything in the knowledge or layer tables declares this label."""
        return self._resolve(feature) is not None

    def match_kind(self, feature: str) -> Optional[str]:
        """HOW ``feature`` resolves: ``gene``, ``composite``, ``layer``,
        ``family`` or ``stem`` - or ``None`` when nothing declares it.

        The kind matters to callers that must tell a determinant apart from a
        state named after it: ``any_MBL`` resolves as a ``composite`` (detecting
        ``blaNDM-1`` IS detecting it) while ``oprD_burden`` resolves as a
        ``stem`` (detecting an intact ``oprD`` is NOT detecting its loss).
        """
        resolved = self._resolve(feature)
        return resolved[0] if resolved is not None else None

    def members_for(self, feature: str) -> Optional[List[str]]:
        """The constituent genes of ``feature``, or ``None`` when it has none.

        ``any_MBL`` -> the MBL genes; ``gene__oprD_absent`` -> ``["oprD"]``;
        ``L1_KPC`` -> the KPC family's genes; ``PDC_high`` -> ``None``, a
        declared composite with no single gene behind it.
        """
        resolved = self._resolve(feature)
        if resolved is None:
            return None
        members = resolved[1]
        return list(members) if members else None

    def _resolve(self, feature: str) -> Optional[Tuple[str, Tuple[str, ...]]]:
        label = normalise_label(feature)
        key = label.casefold()

        for gene in self.genes:
            if gene.casefold() == key:
                return KIND_GENE, (gene,)

        if key in self.composites:
            return KIND_COMPOSITE, self.composites[key]

        if key in self.layers:
            return "layer", self.layers[key]

        family = _family_key(label)
        if key in self.families:
            return "family", self.families[key]
        if family is not None and family in self.families:
            return "family", self.families[family]

        stem = key.split("_", 1)[0]
        for gene in self.genes:
            if gene.casefold() == stem:
                return "stem", (gene,)

        return None


def load_known_determinants(config) -> KnownDeterminants:
    """Load the committed table and build the lookup indexes."""
    path = known_determinants_path(config)
    if not path.is_file():
        # The filename is deliberately not spelled out in the message: this
        # package cannot add it to the input inventory (that list lives in
        # tests/integration/test_seam_real_input_closure.py, which this package
        # does not own), and a refusal naming an uninventoried file fails that
        # guard. The path still reaches the reader through `path=`, which is
        # where a machine-readable location belongs anyway.
        raise DataContractError(
            "The derived known-determinants table is missing; it is derived "
            "from config/mechanisms.tsv and must be committed alongside it",
            path=str(path),
            regenerate="papipeline.downstream.known_determinants"
            ".write_known_determinants(config)",
        )
    rows = read_known_determinants(path)

    genes: List[str] = []
    composites: Dict[str, Tuple[str, ...]] = {}
    for row in rows:
        determinant = str(row.get("determinant") or "")
        if not determinant:
            raise DataContractError(
                "A known-determinants row has no determinant", path=str(path)
            )
        if row.get("kind") == KIND_GENE:
            genes.append(determinant)
        elif row.get("kind") == KIND_COMPOSITE:
            members = tuple(
                m for m in str(row.get("members") or "").split(",") if m.strip()
            )
            composites[determinant.casefold()] = members
        else:
            raise DataContractError(
                "A known-determinants row has an unknown kind",
                path=str(path),
                determinant=determinant,
                kind=row.get("kind"),
                accepted=",".join((KIND_GENE, KIND_COMPOSITE)),
            )

    # Genes grouped by the family the layer encoder names them with:
    # ``blaOXA-48``/``blaOXA-23``/``blaOXA-24/40`` -> ``oxa``, so ``L1_OXA``
    # and a bare ``OXA`` both resolve. Only bla genes have a family.
    families: Dict[str, Tuple[str, ...]] = {}
    for gene in genes:
        family = _family_key(gene)
        if family:
            families[family] = families.get(family, ()) + (gene,)

    return KnownDeterminants(
        rows=tuple(rows),
        genes=tuple(genes),
        composites=composites,
        layers={name.casefold(): tuple(members) for name, members in LAYER_COMPOSITES.items()},
        families=families,
    )
