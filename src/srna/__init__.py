"""srna-annotation: production-quality small-RNA annotation pipeline."""

from .config import PipelineConfig, Strategy, canon_strategy
from .paths import PathResolver, derive_sample_label

__all__ = ["PipelineConfig", "PathResolver", "Strategy", "canon_strategy", "derive_sample_label"]

__version__ = "0.1.0"
