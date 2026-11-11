"""Single evaluation path: benchmark-aligned metrics + final-geometry (layout_quality)."""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

from . import escher_bridge, findings, metrics
from .pathway import read_sbml
from .stability import position_displacement, reproducibility

Positions = dict[str, tuple[float, float]]


def _read_json(path: Path):
    return json.loads(path.read_text())


def _optional_layout_quality(escher_map):
    """Reuse the existing final-geometry evaluator (quality gate); never duplicate it."""
    try:
        import layout_quality
    except ImportError:
        return None, None
    try:
        summary, issues = layout_quality.evaluate(escher_map)
        return summary, issues
    except Exception as ex:  # noqa: BLE001 — evaluation aid must not fail the whole score
        return dict(error=f"{type(ex).__name__}: {ex}"), None


def score_escher_map(
    gem_path: str | Path,
    map_path: str | Path | None = None,
    escher_map=None,
    reference: Positions | None = None,
    reference_path: str | Path | None = None,
    include_geometry: bool = True,
    normalise: bool = True,
    runtime_seconds: float | None = None,
) -> dict[str, Any]:
    """Score one Escher map against its SBML model.

    Returns a single report dict with:
      - ``benchmark``: Project 22 metrics (pathway-specific always; gdMetriX when available)
      - ``geometry``: existing ``layout_quality.evaluate`` summary (curve crossings, collisions, …)
      - ``examples``: reaction ids / regions for the weakest spots
      - ``coverage`` / ``runtime_seconds``
    """
    gem_path = Path(gem_path)
    pathway = read_sbml(gem_path)
    if escher_map is None:
        if map_path is None:
            raise ValueError("map_path or escher_map is required")
        escher_map = escher_bridge.load_escher(map_path)
    t0 = time.perf_counter()
    pos = escher_bridge.escher_positions(pathway, escher_map, y_up=True)
    widths = {n: float(pathway.graph.nodes[n]["width"]) for n in pos if n in pathway.graph}
    if normalise:
        pos_metric = escher_bridge.normalise_positions(pos, widths)
    else:
        pos_metric = pos

    ref_pos = reference
    if reference_path is not None and ref_pos is None:
        ref_pos = _load_reference(pathway, reference_path)

    # Restrict pathway to placed nodes for fair graph metrics
    placed = [n for n in pathway.graph if n in pos_metric]
    sub = pathway.graph.subgraph(placed).copy()
    from .pathway import Pathway

    sub_pathway = Pathway(name=pathway.name, graph=sub, compartments=pathway.compartments, source=pathway.source)

    benchmark = metrics.all_benchmark_metrics(sub_pathway, pos_metric, reference=ref_pos)
    geom_summary, geom_issues = (None, None)
    # Final-geometry (curve sampling) is O(segments²); skip automatically on large GEMs unless forced.
    n_rxn = len(pathway.reactions)
    if include_geometry and n_rxn > 800:
        geom_summary = dict(skipped=True, reason=f"{n_rxn} reactions: use layout_issues.json from compose, or --force-geometry")
    elif include_geometry:
        geom_summary, geom_issues = _optional_layout_quality(escher_map)

    examples = findings.weak_spots(sub_pathway, pos_metric, issues=geom_issues)
    cov = escher_bridge.coverage(pathway, pos)
    seconds = runtime_seconds
    if seconds is None and map_path is not None:
        # Prefer compose report timing when present next to the map
        report = Path(map_path).with_name("report.json")
        if report.exists():
            try:
                rj = _read_json(report)
                fg = ((rj.get("layout_quality") or {}).get("final_geometry") or {})
                seconds = (
                    rj.get("seconds")
                    or (rj.get("timings") or {}).get("total")
                    or fg.get("seconds")
                )
            except Exception:
                seconds = None

    return dict(
        model=pathway.name,
        gem=str(gem_path.resolve()),
        map=str(Path(map_path).resolve()) if map_path else None,
        coverage=cov,
        runtime_seconds=seconds,
        metric_seconds=round(time.perf_counter() - t0, 3),
        dependencies=dict(gdmetrix=metrics.gdmetrix_available(), shapely=metrics.shapely_available()),
        benchmark=benchmark,
        geometry=geom_summary,
        examples=examples,
        note=(
            "Development evaluation (not a fair tool comparison). "
            "Atlas sources are never excluded for evaluation overlap. "
            "Pathway metrics adapted from elixir2026-project22/layout-benchmark (MIT); "
            "gdMetriX is optional (GPL-3.0)."
        ),
    )


def score_paths(
    gem: str | Path,
    map_json: str | Path,
    reference: str | Path | None = None,
    compare_map: str | Path | None = None,
    include_geometry: bool = True,
) -> dict[str, Any]:
    """CLI-friendly entry: paths in, report dict out. Optional compare_map → stability block."""
    report = score_escher_map(
        gem,
        map_json,
        reference_path=reference,
        include_geometry=include_geometry,
    )
    if compare_map:
        pathway = read_sbml(gem)
        a = escher_bridge.escher_positions(pathway, escher_bridge.load_escher(map_json))
        b = escher_bridge.escher_positions(pathway, escher_bridge.load_escher(compare_map))
        report["stability_vs_compare"] = position_displacement(a, b)
        report["reproducibility"] = reproducibility([a, b])
    return report


def _load_reference(pathway, path: str | Path) -> Positions | None:
    path = Path(path)
    if path.suffix.lower() == ".json":
        data = _read_json(path)
        # Escher map or {node: [x,y]}
        if isinstance(data, list) or (isinstance(data, dict) and "nodes" in data):
            return escher_bridge.escher_positions(pathway, data)
        if isinstance(data, dict):
            out = {}
            for k, v in data.items():
                if isinstance(v, (list, tuple)) and len(v) >= 2:
                    out[k] = (float(v[0]), float(v[1]))
            return out or None
    if path.suffix.lower() == ".sbgn":
        try:
            from .reference_sbgn import gallery_reference

            return gallery_reference(pathway, path)
        except Exception:
            return None
    return None
