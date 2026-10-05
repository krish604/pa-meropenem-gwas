"""Test-support utilities: the synthetic fixture generator.

Nothing in this package is used by the analysis path. It exists so the
pipeline can be exercised end to end without touching real data.
"""

from __future__ import annotations

from .synthetic import (
    DEFAULT_N_SAMPLES,
    DEFAULT_SEED,
    SAMPLE_PREFIX,
    SYNTHETIC_BANNER,
    SyntheticSample,
    build_samples,
    generate,
    sample_id_for,
)

__all__ = [
    "DEFAULT_N_SAMPLES",
    "DEFAULT_SEED",
    "SAMPLE_PREFIX",
    "SYNTHETIC_BANNER",
    "SyntheticSample",
    "build_samples",
    "generate",
    "sample_id_for",
]
