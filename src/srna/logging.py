"""Logging setup: one small helper to configure loguru for the pipeline."""

from __future__ import annotations

import sys
from pathlib import Path

from loguru import logger


def setup_logging(log_path: Path | None = None, level: str = "INFO") -> None:
    """Configure loguru to write to stderr and, optionally, a rotating log file.

    Safe to call multiple times; existing handlers are removed first so logs
    are never duplicated.
    """
    logger.remove()
    logger.add(
        sys.stderr,
        level=level,
        format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | {message}",
    )
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        logger.add(
            str(log_path),
            level=level,
            format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
            rotation="100 MB",
            retention="30 days",
            enqueue=True,
        )
