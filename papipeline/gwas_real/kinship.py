"""The stage-10 kinship matrix: patristic distances, transformed for pyseer.

Spec user story 41 wants the similarity (kinship) matrix ``--lmm`` consumes to
be built from *the same tree* stage 10 reports, so the two cannot disagree.
pyseer ships no tree-to-similarity helper - the installed 1.1.2 declares five
console scripts and ``python -m pyseer.similarity`` takes variant files, not a
tree - so the conversion is this pipeline's own, and it is the classical one:

    K = -1/2 * (D^2 - row_mean - col_mean + grand_mean)

over ALL n entries of the squared distance matrix, diagonal included. That
double centring - not a mere demeaning - is what makes every row sum to zero
and, because a tree metric is of negative type, what makes ``K`` positive
semi-definite: the property ``pyseer/fastlmm/lmm_cov.py`` documents for its
kernel ("positive semi-definite Kernel matrix"). Both facts are *checked* here
rather than assumed: a trace that is not positive would make pyseer compute
``factor = float(len(p)) / np.diag(K.values).sum()`` (a division by zero deep
inside its fit), and a negative eigenvalue would make ``la.eigh`` measure this
module's bug instead of the cohort's relatedness. Both are fail-closed.

**The distances are stage 10's, not a second traversal.** This module calls
``stages.similarity.patristic_distances`` itself - one test in
``tests/unit/test_gwas_real_kinship.py`` spies on that call - because two
parsers of the same Newick that agree today can drift apart the moment one is
optimised, and nothing downstream would notice.

**Extra tips are allowed; missing ones are not.** The stage-9 tree may cover
more isolates than the analysable cohort (a non-binary phenotype excludes
samples from the model but not from the tree). The matrix is therefore taken
over the requested ids, and centred over exactly those: the kernel pyseer
rescales is defined on the samples it analyses. A *cohort* sample that is not
a tip is the other case, and it is a mismatch - silently dropping it would hand
pyseer a kinship matrix describing a different cohort than the phenotype file.

The sidecar beside the matrix follows the convention ``stages/similarity.py``
set for units, extended with provenance: a bare matrix of small decimals is
exactly the ambiguous case (substitutions per site, squared, centred - not
SNPs), and the file says so.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence

import numpy as np

from ..config.loader import PipelineConfig
from ..errors import DataContractError
from ..manifest import SampleManifest
from ..models import RunMode
from ..stages import similarity as stage_similarity
from ..logging_utils import get_logger

LOGGER = get_logger("gwas_real.kinship")

#: How ``K`` was obtained from the distances. Recorded so a reader never has to
#: guess whether these numbers came from the tree, from a feature matrix, or
#: from somewhere else.
METHOD = "gower_centered_patristic"

#: Leading column of the written matrix, and the key of each returned row.
ID_COLUMN = stage_similarity.ID_COLUMN

#: Key of a row's own similarity vector in the rows :func:`run` returns.
SIMILARITIES = "similarities"

#: Sidecar suffix, derived from the matrix's name so the pair cannot drift.
SIDECAR_SUFFIX = ".meta.json"

#: Eigenvalues below ``-PSD_TOLERANCE * scale`` fail closed. ``scale`` is at
#: least 1 so a matrix whose entries are all ~1e-12 is not judged against a
#: tolerance smaller than its own float noise.
PSD_TOLERANCE = 1e-8


def similarity_from_tree(
    tree_text: str, ids: Sequence[str]
) -> Dict[str, Dict[str, float]]:
    """The Gower-centred kinship matrix for ``ids``, read off a Newick tree.

    Args:
        tree_text: The stage-9 Newick string.
        ids: The cohort, in the order its rows must come out.

    Returns:
        ``{sample: {sample: k}}``, ordered as ``ids``.

    Raises:
        DataContractError: A cohort sample is not a tip in the tree; the
            transform's trace is not positive; or the matrix is not positive
            semi-definite within float tolerance.
    """
    ids = [str(i) for i in ids]
    if not ids:
        raise DataContractError(
            "Cannot build a kinship matrix for an empty cohort: there is no "
            "relatedness to describe and pyseer would divide by a zero trace.",
            n_samples=0,
        )
    if len(set(ids)) != len(ids):
        raise DataContractError(
            "The cohort contains duplicate sample ids, so it cannot name the "
            "rows of a square matrix.",
            n_samples=len(ids),
            n_distinct=len(set(ids)),
        )

    # The stage-10 computation itself, exactly once - see the module docstring.
    distances = stage_similarity.patristic_distances(tree_text)

    missing = [s for s in ids if s not in distances]
    if missing:
        raise DataContractError(
            f"{len(missing)} cohort sample(s) are not tips in the stage-9 tree: "
            f"{', '.join(missing[:10])}"
            + (" ..." if len(missing) > 10 else "")
            + ". The kinship matrix has to describe the same cohort as the "
            "phenotype file, so a sample in one and not the other is a "
            "mismatch to stop on - not a sample to drop. Dropping it would "
            "make pyseer analyse a relatedness structure that does not "
            "belong to the cohort it reports.",
            missing=",".join(missing[:10]),
            n_missing=len(missing),
            n_tips=len(distances),
        )

    n = len(ids)
    squared = [
        [distances[a][b] ** 2 for b in ids] for a in ids
    ]
    row_means = [sum(row) / n for row in squared]
    grand = sum(sum(row) for row in squared) / (n * n)

    kinship: Dict[str, Dict[str, float]] = {a: {} for a in ids}
    for i, a in enumerate(ids):
        for j, b in enumerate(ids):
            kinship[a][b] = -0.5 * (
                squared[i][j] - row_means[i] - row_means[j] + grand
            )

    trace = float(sum(kinship[a][a] for a in ids))
    if not trace > 0.0:
        raise DataContractError(
            f"The kinship matrix built from the stage-9 tree has trace {trace!r}, "
            "so the tree measures no relatedness at all (every distance is "
            "zero, or the centring cancelled). pyseer divides by this trace - "
            "factor = len(p) / np.diag(K.values).sum() in pyseer/lmm.py - and a "
            "zero trace would surface as a tool traceback rather than as a "
            "statement about the tree.",
            trace=trace,
            n_samples=n,
            method=METHOD,
        )

    matrix = np.asarray([[kinship[a][b] for b in ids] for a in ids], dtype=float)
    eigenvalues = np.linalg.eigvalsh(matrix)
    scale = max(1.0, float(np.abs(eigenvalues).max()))
    if float(eigenvalues.min()) < -PSD_TOLERANCE * scale:
        raise DataContractError(
            "The kinship matrix built from the stage-9 tree is not positive "
            f"semi-definite (minimum eigenvalue {float(eigenvalues.min()):.3E} "
            f"against a tolerance of {-PSD_TOLERANCE * scale:.3E}). Gower "
            "centring of a tree metric is PSD by construction, so a negative "
            "eigenvalue this large means the distances are not a tree metric - "
            "and pyseer's eigendecomposition would be fitting this defect "
            "rather than the cohort's relatedness.",
            min_eigenvalue=float(eigenvalues.min()),
            tolerance=-PSD_TOLERANCE * scale,
            method=METHOD,
        )
    return kinship


def _format(value: float) -> str:
    return f"{float(value):.10g}"


def units_sidecar_path(matrix_path: Path) -> Path:
    """`kinship.tsv` -> `kinship.tsv.meta.json`, beside it.

    Named ``.meta.json`` rather than the stage's ``.units.json`` because this
    record carries both the units and the provenance of the transform - which
    tree, which method, and the note about pyseer's trace rescaling below.
    """
    matrix_path = Path(matrix_path)
    return matrix_path.with_name(matrix_path.name + SIDECAR_SUFFIX)


def write_matrix(
    path: Path, ids: Sequence[str], kinship: Mapping[str, Mapping[str, float]],
    source_tree: Path,
) -> Path:
    """Write the square similarity matrix and its sidecar.

    The shape is the one ``pyseer/lmm.py::initialise_lmm`` reads back:
    ``pd.read_table(K_in, index_col=0)`` with ``K.index.astype(str)``, rows in
    the cohort's own order. Six-ten significant digits rather than a rounded
    display value: the fit rescales by trace, so precision here is precision in
    the kernel itself.
    """
    ids = [str(i) for i in ids]
    if not ids:
        raise DataContractError("Cannot write an empty kinship matrix")
    path = Path(path)

    lines = ["\t".join([ID_COLUMN, *ids])]
    for a in ids:
        lines.append("\t".join([a] + [_format(kinship[a][b]) for b in ids]))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    payload: Dict[str, Any] = {
        "matrix": path.name,
        "method": METHOD,
        "quantity": "Gower-centred kinship (similarity) between two samples",
        "input_unit": stage_similarity.DISTANCE_UNIT,
        "is_a_snp_count": False,
        "source_tree": str(source_tree),
        "n_samples": len(ids),
        "note": (
            "Entries are the stage-10 patristic distances squared (so the raw "
            "unit is (expected substitutions per site)^2) and double-centred by "
            "the Gower transform; they are NOT SNP counts and NOT distances. "
            "pyseer rescales K so its trace equals the sample count "
            "(pyseer/lmm.py: factor = len(p) / np.diag(K.values).sum()), so the "
            "scale of these raw entries does not reach any fitted value - what "
            "must be true before that rescale is that the TRACE is positive and "
            "the matrix is positive semi-definite, and both are checked before "
            "this file is written."
        ),
    }
    sidecar = units_sidecar_path(path)
    sidecar.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    return path


def run(
    config: PipelineConfig,
    manifest: SampleManifest,
    mode: RunMode,
    *,
    tree_path: Path | None = None,
    out_path: Path | None = None,
) -> List[Dict[str, Any]]:
    """Stage-shaped entry point: build the kinship matrix for this cohort.

    Args:
        config: Loaded configuration. REAL reads the stage-9 tree location
            from it; TEST takes the tree from ``tree_path``.
        manifest: The cohort. Every sample must be a tip in the tree.
        mode: ``TEST`` reads the committed fixture tree. ``REAL`` is gated on
            ``runtime.allow_real_mode`` and then reads the configured tree.
            ``STUB`` has no matrix to produce and is refused.
        tree_path: The stage-9 Newick tree (TEST only; no default).
        out_path: Where to write the matrix and its sidecar. Omit to skip
            writing.

    Returns:
        One row per manifest sample: ``{sample_id, similarities}``.
    """
    if mode is RunMode.TEST:
        if tree_path is None:
            raise DataContractError(
                "TEST-mode kinship needs `tree_path`: it reads the committed "
                "stage-9 fixture tree. There is no default, because resolving "
                "one from configuration would let a TEST run measure a tree "
                "nobody chose."
            )
        tree = Path(tree_path)
    elif mode is RunMode.REAL:
        if not bool(config.runtime.get("allow_real_mode", False)):
            raise NotImplementedError(
                "REAL-mode kinship is gated: runtime.allow_real_mode is false. "
                "Set it in the machine overlay, or open it for one session with "
                "PIPELINE_ALLOW_REAL_MODE=1. The committed overlays keep it shut "
                "so a REAL run cannot begin by accident. "
                f"(mode={getattr(mode, 'value', mode)})"
            )
        tree = config.phylogeny_dir(mode) / "tree.nwk"
    else:
        raise NotImplementedError(
            "STUB mode produces contract headers and zero rows; it does not "
            "compute a kinship matrix, and an empty one would tell pyseer every "
            "isolate is identical. Use TEST on the committed fixtures, or a "
            "REAL run with runtime.allow_real_mode open. "
            f"(mode={getattr(mode, 'value', mode)})"
        )

    if not tree.is_file():
        raise DataContractError(
            f"Stage-9 tree not found at {tree} (expected file name: tree.nwk). "
            "The kinship/similarity matrix pyseer's --lmm consumes is read off "
            "that tree; without it there is nothing to measure. A missing tree "
            "is a missing intermediate, not an absence of similarity.",
            path=str(tree),
            mode=getattr(mode, "value", mode),
        )

    ids = [str(s) for s in manifest.sample_ids]
    newick = tree.read_text(encoding="utf-8")
    kinship = similarity_from_tree(newick, ids)

    rows: List[Dict[str, Any]] = [
        {ID_COLUMN: sample, SIMILARITIES: dict(kinship[sample])}
        for sample in ids
    ]
    LOGGER.info(
        "Stage 10 kinship: %d x %d %s matrix from %s", len(ids), len(ids), METHOD, tree
    )
    if out_path is not None:
        write_matrix(Path(out_path), ids, kinship, tree)
    return rows


__all__ = [
    "ID_COLUMN",
    "METHOD",
    "SIMILARITIES",
    "run",
    "similarity_from_tree",
    "units_sidecar_path",
    "write_matrix",
]
