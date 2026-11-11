"""Unified layout evaluation for Open Meta Maps.

Scores Escher maps with Project 22 layout-benchmark rules (crossings, overlaps, long
edges, participant consistency, compartments, aspect/compactness, reference similarity)
plus the existing final-geometry evaluator used by the quality gate.

Pathway-specific metric code is ported from
https://github.com/elixir2026-project22/layout-benchmark (MIT); see ATTRIBUTION.md.
gdMetriX (GPL-3.0) is an optional import only and is never vendored.
"""

from .score import score_escher_map, score_paths

__all__ = ["score_escher_map", "score_paths"]
