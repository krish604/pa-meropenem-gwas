"""The REAL-path pyseer adapter: stage-10 kinship and stage-12 invocation.

Two members, both package-level so callers and tests name them the same way:

* :mod:`papipeline.gwas_real.kinship` - the similarity (kinship) matrix
  ``--lmm`` consumes, built from the stage-9 tree by reusing stage 10's own
  patristic-distance traversal (``stages.similarity.patristic_distances``)
  and the classical Gower transform. ``papipeline/stages/similarity.py`` is
  the stage; this module is what pyseer needs from it.
* :mod:`papipeline.gwas_real.adapter` - the two-pass pyseer invocation with
  its mode gates, its input-cohort checks and step 12a's threshold reused via
  ``stages.gwas.reduce_unique_patterns``.

Nothing here is wired into ``workflow/Snakefile`` or the ``STAGE_*``
registries; that wiring is a separate package and is recorded in
``.build/gwas-engine.registry-notes.md``.
"""

from . import adapter, kinship  # noqa: F401  (the package surface, by name)
from .adapter import (  # noqa: F401
    OPT_IN_ENV,
    RealGwasRun,
    build_command,
    run,
)
from .kinship import (  # noqa: F401
    METHOD,
    units_sidecar_path,
)

__all__ = [
    "METHOD",
    "OPT_IN_ENV",
    "RealGwasRun",
    "adapter",
    "build_command",
    "kinship",
    "run",
    "units_sidecar_path",
]
