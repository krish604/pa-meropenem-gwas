"""Typed exception hierarchy for the pipeline.

Every error carries a human-readable message plus optional context so that
failures are actionable without reading a traceback.
"""

from __future__ import annotations

from typing import Any


class PipelineError(Exception):
    """Base class for all pipeline errors."""

    def __init__(self, message: str, **context: Any) -> None:
        self.message = message
        self.context = context
        super().__init__(self._render())

    def _render(self) -> str:
        if not self.context:
            return self.message
        rendered = ", ".join(f"{k}={v!r}" for k, v in sorted(self.context.items()))
        return f"{self.message} [{rendered}]"


class ConfigError(PipelineError):
    """Configuration file is missing, malformed, or internally inconsistent."""


class SampleCapExceeded(PipelineError):
    """A run asks a machine to carry more samples than it is configured for.

    This is a refusal, not a crash: the machine overlay sets the cap and the
    message says which file to change or which machine to use instead.
    """


class DataContractError(PipelineError):
    """An input file violates the documented data contract."""


class SampleIdError(DataContractError):
    """A sample identifier is malformed."""


class DuplicateSampleError(DataContractError):
    """A sample identifier appears more than once in a single-source table."""


class PhenotypeError(DataContractError):
    """Phenotype value is not in the configured allowed set."""


class AntibioticError(DataContractError):
    """Antibiotic is not configured for this project."""


class KnowledgeBaseError(PipelineError):
    """A knowledge table (mechanisms/regulators/references) is invalid."""


class UnknownGeneError(PipelineError):
    """A gene has no entry in the knowledge tables.

    Unknown genes are *not* silently dropped. They are surfaced so the
    knowledge tables can be extended deliberately.
    """


class TreeSampleMismatchError(DataContractError):
    """Phylogenetic tree tip labels do not match the master sample manifest."""


class ToolNotAvailableError(PipelineError):
    """An external tool required by an enabled stage is not installed."""


class ToolExecutionError(PipelineError):
    """An external tool ran but exited non-zero."""


class ModeNotAllowedError(PipelineError):
    """An operation was requested that the current run mode forbids."""


class StageError(PipelineError):
    """A pipeline stage failed."""
