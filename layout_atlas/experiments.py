"""Compatibility module for callers that expect ``experiments.MODELS``.

The offline experiment package lives at the repository root as ``experiments/``.
This file only keeps the model list that ``micromap_experiment`` and similar
legacy entry points import as ``experiments.MODELS``.
"""
from core import EXAMPLE_MODELS as MODELS  # noqa: F401
