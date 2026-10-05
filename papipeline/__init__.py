"""Pseudomonas aeruginosa antimicrobial resistance analysis pipeline.

Imipenem-focused, assembly-based, configuration-driven.

The package is organised as:

* :mod:`papipeline.models`      - typed records and controlled vocabularies
* :mod:`papipeline.io`          - strict TSV/FASTA readers and writers
* :mod:`papipeline.config`      - configuration and knowledge-table loading
* :mod:`papipeline.knowledge`   - knowledge-table query helpers
* :mod:`papipeline.manifest`    - the master sample manifest (join key source)
* :mod:`papipeline.stages`      - the sixteen analysis stages
* :mod:`papipeline.adapters`    - external tool detection and adapters
* :mod:`papipeline.viz`         - data preparation for figures
* :mod:`papipeline.run`         - the native orchestrator
"""

from __future__ import annotations

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "errors",
    "logging_utils",
    "models",
]
