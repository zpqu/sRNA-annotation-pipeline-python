"""Command-line interface: ``srna`` entry point.

Runs the full small-RNA annotation pipeline. Options include the genome
assembly, sense-overlap strategy, chromosome naming convention and more.
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

from loguru import logger

from srna.config import PipelineConfig, canon_strategy
from srna.logging import setup_logging
from srna.pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    """Return the ``srna`` argument parser."""
    parser = argparse.ArgumentParser(
        prog="srna",
        description="Small-RNA annotation pipeline.",
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path.cwd(),
        help="project root containing data/ and output/ (default: current directory)",
    )
    parser.add_argument(
        "--genome",
        default="mm39",
        help="reference genome assembly (default: mm39)",
    )
    parser.add_argument(
        "--strategy",
        default="fully-contained",
        help="sense-overlap rule: fully-contained | union | any | comparison "
        "(default: fully-contained)",
    )
    parser.add_argument(
        "--chr-style",
        default="chr",
        help="chromosome naming convention: chr | none (default: chr)",
    )
    parser.add_argument(
        "--force-rebuild-db",
        action="store_true",
        help="rebuild the feature DB even when caches are valid",
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=None,
        help="directory for pipeline logs (default: <root>/logs)",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    """Run the pipeline from the command line."""
    args = build_parser().parse_args(argv)
    config = PipelineConfig(
        genome=args.genome,
        strategy=canon_strategy(args.strategy),
        chr_style=args.chr_style,
        force_rebuild_db=args.force_rebuild_db,
    )
    log_dir = args.log_dir or (args.root / "logs")
    log_path = log_dir / f"pipeline_{datetime.now():%Y%m%d_%H%M%S}.log"
    setup_logging(log_path)

    logger.info(
        "pipeline start: genome={} strategy={} root={}",
        config.genome,
        config.strategy.value,
        args.root,
    )
    run_pipeline(config, args.root)
    logger.info("pipeline finished; log written to {}", log_path)


if __name__ == "__main__":
    main()
