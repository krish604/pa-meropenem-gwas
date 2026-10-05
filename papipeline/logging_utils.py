"""Logging setup shared by the library, the CLI wrappers and the workflow.

A single ``get_logger`` factory is used everywhere so that log format and
level are configured in exactly one place.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Optional

LOGGER_NAME = "papipeline"
_DEFAULT_FORMAT = "%(asctime)s %(levelname)-7s %(name)s: %(message)s"
_CONFIGURED = False


def configure_logging(
    level: str = "INFO",
    logfile: Optional[Path] = None,
    stream: bool = True,
) -> logging.Logger:
    """Configure the package logger. Safe to call more than once."""
    global _CONFIGURED

    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(level.upper())

    if not _CONFIGURED:
        handler: logging.Handler
        if stream:
            handler = logging.StreamHandler(sys.stderr)
            handler.setFormatter(logging.Formatter(_DEFAULT_FORMAT))
            logger.addHandler(handler)
        _CONFIGURED = True

    if logfile is not None:
        logfile = Path(logfile)
        logfile.parent.mkdir(parents=True, exist_ok=True)
        existing = {
            getattr(h, "baseFilename", None) for h in logger.handlers
        }
        if str(logfile.resolve()) not in existing:
            file_handler = logging.FileHandler(logfile)
            file_handler.setFormatter(logging.Formatter(_DEFAULT_FORMAT))
            file_handler.setLevel(level.upper())
            logger.addHandler(file_handler)

    logger.setLevel(level.upper())
    return logger


def get_logger(name: Optional[str] = None) -> logging.Logger:
    """Return a child logger under the package namespace."""
    if name is None:
        return logging.getLogger(LOGGER_NAME)
    return logging.getLogger(f"{LOGGER_NAME}.{name}")
