"""Strict I/O helpers: TSV and FASTA.

Public surface::

    from papipeline.io import read_tsv, write_tsv, read_fasta, assembly_stats
"""

from __future__ import annotations

from .fasta import (
    AMBIGUITY_CODES,
    AssemblySummary,
    FastaRecord,
    assembly_stats,
    read_fasta,
    summarise_alignment,
)
from .tsv import (
    MISSING_SENTINELS,
    iter_tsv_dicts,
    read_tsv,
    write_fasta,
    write_lines,
    write_tsv,
)

__all__ = [
    "AMBIGUITY_CODES",
    "AssemblySummary",
    "FastaRecord",
    "MISSING_SENTINELS",
    "assembly_stats",
    "iter_tsv_dicts",
    "read_fasta",
    "read_tsv",
    "summarise_alignment",
    "write_fasta",
    "write_lines",
    "write_tsv",
]
