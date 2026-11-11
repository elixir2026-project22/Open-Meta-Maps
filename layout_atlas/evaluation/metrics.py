"""Layout quality metrics aligned with Project 22 layout-benchmark rules.

Pathway-specific metrics are adapted from layout-benchmark `metrics.py` (MIT).
General graph-drawing metrics come from gdMetriX when installed (GPL-3.0; optional
import only — never vendored).
"""

from __future__ import annotations

import math
import os
from typing import Any

import numpy as np
from scipy.spatial import ConvexHull, procrustes
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr

from .pathway import PRODUCT, REACTANT, SPECIES, Pathway

Positions = dict[str, tuple[float, float]]
SYMMETRY_MAX_NODES = 300
CROSSING_MAX_EDGES = 1200  # number_of_crossings is quadratic; skip on large GEMs

try:
    from shapely import STRtree, box
    from shapely.geometry import Point, Polygon

    _HAS_SHAPELY = True
except ImportError:  # pragma: no cover - optional for core installs
    _HAS_SHAPELY = False

try:
    import gdMetriX

    _HAS_GDMETRIX = True
except ImportError:  # pragma: no cover
    gdMetriX = None  # type: ignore
    _HAS_GDMETRIX = False


def gdmetrix_available() -> bool:
    return _HAS_GDMETRIX


def shapely_available() -> bool:
    return _HAS_SHAPELY


def general_metrics(pathway: Pathway, pos: Positions) -> dict[str, float | None]:
    """gdMetriX metrics on the undirected graph. Values are None when gdMetriX is missing."""
    if not _HAS_GDMETRIX:
        return {
            "crossings": None,
            "crossing_density": None,
            "crossing_angular_resolution": None,
            "angular_resolution": None,
            "edge_length_deviation": None,
            "edge_orthogonality": None,
            "stress": None,
            "aspect_ratio": None,
            "area_per_node": None,
            "homogeneity": None,
            "coherence_to_average_flow": None,
            "visual_symmetry": None,
            "gdmetrix": 0.0,
        }
    g = pathway.undirected().subgraph([n for n in pathway.undirected() if n in pos]).copy()
    # gdMetriX needs every graph node in pos
    missing = [n for n in g if n not in pos]
    g.remove_nodes_from(missing)
    n = g.number_of_nodes()
    if n < 2 or g.number_of_edges() < 1:
        return {k: None for k in (
            "crossings", "crossing_density", "crossing_angular_resolution", "angular_resolution",
            "edge_length_deviation", "edge_orthogonality", "stress", "aspect_ratio", "area_per_node",
            "homogeneity", "coherence_to_average_flow", "visual_symmetry",
        )} | {"gdmetrix": 1.0}
    directed = pathway.graph.subgraph([n for n in pathway.graph if n in pos]).copy()
    e = g.number_of_edges()
    do_cross = e <= CROSSING_MAX_EDGES

    def _f(fn, *a, **kw):
        try:
            return fn(*a, **kw)
        except Exception:
            return None

    m = {
        "crossings": _f(gdMetriX.number_of_crossings, g, pos=pos) if do_cross else None,
        "crossing_density": _f(gdMetriX.crossing_density, g, pos=pos) if do_cross else None,
        "crossing_angular_resolution": _f(gdMetriX.crossing_angular_resolution, g, pos=pos) if do_cross else None,
        "angular_resolution": _f(gdMetriX.angular_resolution, g, pos=pos),
        "edge_length_deviation": _f(gdMetriX.edge_length_deviation, g, pos=pos),
        "edge_orthogonality": _f(gdMetriX.edge_orthogonality, g, pos=pos),
        "stress": _f(gdMetriX.stress, g, pos=pos),
        "aspect_ratio": _f(gdMetriX.aspect_ratio, g, pos=pos),
        "area_per_node": (lambda a: (a / max(n, 1)) if a is not None else None)(_f(gdMetriX.area, g, pos=pos)),
        "homogeneity": _f(gdMetriX.homogeneity, g, pos=pos),
        "coherence_to_average_flow": _f(gdMetriX.coherence_to_average_flow, directed, pos=pos),
        # visual_symmetry opens a matplotlib figure; keep it opt-in via env OPM_EVAL_SYMMETRY=1
        "visual_symmetry": (
            _f(gdMetriX.visual_symmetry, g, pos=pos)
            if n <= SYMMETRY_MAX_NODES and os.environ.get("OPM_EVAL_SYMMETRY") == "1"
            else None
        ),
        "gdmetrix": 1.0,
        "gdmetrix_crossings_computed": 1.0 if do_cross else 0.0,
    }
    return {k: (float(v) if v is not None else None) for k, v in m.items()}


def _boxes(pathway: Pathway, pos: Positions):
    nodes = [n for n in pathway.graph if n in pos]
    boxes = []
    for n in nodes:
        d = pathway.graph.nodes[n]
        x, y = pos[n]
        w, h = d["width"] / 2, d["height"] / 2
        boxes.append(box(x - w, y - h, x + w, y + h))
    return nodes, boxes


def overlaps(pathway: Pathway, pos: Positions) -> dict[str, float | None]:
    """Pairs of node boxes that overlap (touching is no overlap), and overlap area share."""
    if not _HAS_SHAPELY:
        return {"overlapping_pairs": None, "overlap_area_share": None}
    nodes, boxes_ = _boxes(pathway, pos)
    if not boxes_:
        return {"overlapping_pairs": 0.0, "overlap_area_share": 0.0}
    tree = STRtree(boxes_)
    pairs, area = 0, 0.0
    for i, b in enumerate(boxes_):
        for j in tree.query(b):
            if j > i:
                inter = b.intersection(boxes_[j]).area
                if inter > 0:
                    pairs += 1
                    area += inter
    total = sum(b.area for b in boxes_)
    return {"overlapping_pairs": float(pairs), "overlap_area_share": area / total if total else 0.0}


def participant_consistency(pathway: Pathway, pos: Positions) -> float | None:
    """Share of (reactant, product) pairs on opposite sides of their reaction (>90°)."""
    g = pathway.graph
    opposite, pairs = 0, 0
    for r in pathway.reactions:
        if r not in pos:
            continue
        c = np.array(pos[r], dtype=float)
        reactants = [np.array(pos[u]) - c for u, _, d in g.in_edges(r, data=True) if d.get("role") == REACTANT and u in pos]
        products = [np.array(pos[v]) - c for _, v, d in g.out_edges(r, data=True) if d.get("role") == PRODUCT and v in pos]
        for a in reactants:
            for b in products:
                pairs += 1
                opposite += float(np.dot(a, b)) < 0
    return opposite / pairs if pairs else None


def participant_consistency_details(pathway: Pathway, pos: Positions, limit: int = 40) -> list[dict[str, Any]]:
    """Worst reactions by opposite-side share (for localisation in the report)."""
    g = pathway.graph
    rows = []
    for r in pathway.reactions:
        if r not in pos:
            continue
        c = np.array(pos[r], dtype=float)
        reactants = [(u, np.array(pos[u]) - c) for u, _, d in g.in_edges(r, data=True) if d.get("role") == REACTANT and u in pos]
        products = [(v, np.array(pos[v]) - c) for _, v, d in g.out_edges(r, data=True) if d.get("role") == PRODUCT and v in pos]
        if not reactants or not products:
            continue
        pairs = opposite = 0
        for _, a in reactants:
            for _, b in products:
                pairs += 1
                opposite += float(np.dot(a, b)) < 0
        share = opposite / pairs
        rows.append(
            dict(
                reaction=r,
                label=g.nodes[r].get("label") or r,
                opposite_share=round(share, 3),
                pairs=pairs,
                x=round(float(pos[r][0]), 1),
                y=round(float(pos[r][1]), 1),
            )
        )
    rows.sort(key=lambda r: (r["opposite_share"], -r["pairs"], r["reaction"]))
    return rows[:limit]


def compartment_metrics(pathway: Pathway, pos: Positions) -> dict[str, float | None]:
    """Share of species inside another compartment's hull (0 = well separated)."""
    if not _HAS_SHAPELY:
        return {"species_in_other_compartment": None}
    by_comp: dict[str, list[str]] = {}
    for n, d in pathway.graph.nodes(data=True):
        if d["kind"] == SPECIES and d.get("compartment") and n in pos:
            by_comp.setdefault(d["compartment"], []).append(n)
    if len(by_comp) < 2:
        return {"species_in_other_compartment": None}
    hulls = {}
    for c, members in by_comp.items():
        pts = np.array([pos[n] for n in members])
        if len(pts) >= 3:
            try:
                hulls[c] = Polygon(pts[ConvexHull(pts).vertices])
            except (ValueError, RuntimeError):
                continue
    species = [n for members in by_comp.values() for n in members]
    if not species:
        return {"species_in_other_compartment": None}
    intruders = sum(
        1
        for n in species
        if any(
            c != pathway.graph.nodes[n]["compartment"] and h.contains(Point(pos[n]))
            for c, h in hulls.items()
        )
    )
    return {"species_in_other_compartment": intruders / len(species)}


def reference_similarity(pos: Positions, reference: Positions, k: int = 5) -> dict[str, float | None]:
    """Similarity after Procrustes alignment: disparity, distance Spearman, kNN preservation."""
    common = sorted(set(pos) & set(reference))
    if len(common) < 3:
        return {"procrustes": None, "distance_correlation": None, "knn_preservation": None}
    a = np.array([pos[n] for n in common], dtype=float)
    b = np.array([reference[n] for n in common], dtype=float)
    _, _, disparity = procrustes(b, a)
    da, db = pdist(a), pdist(b)
    corr = spearmanr(da, db).statistic
    k = min(k, len(common) - 1)

    def knn(xy):
        d = np.linalg.norm(xy[:, None, :] - xy[None, :, :], axis=2)
        np.fill_diagonal(d, math.inf)
        return [set(np.argsort(row)[:k]) for row in d]

    ka, kb = knn(a), knn(b)
    preservation = float(np.mean([len(x & y) / k for x, y in zip(ka, kb)]))
    return {
        "procrustes": float(disparity),
        "distance_correlation": float(corr),
        "knn_preservation": preservation,
        "reference_nodes": float(len(common)),
    }


def long_edge_metrics(pathway: Pathway, pos: Positions, long_factor: float = 3.0) -> dict[str, float | None]:
    """Edge-length statistics on placed undirected edges. ``long`` = length > long_factor × median."""
    lengths = []
    for u, v in pathway.undirected().edges():
        if u in pos and v in pos:
            lengths.append(float(np.hypot(pos[u][0] - pos[v][0], pos[u][1] - pos[v][1])))
    if not lengths:
        return {
            "edge_length_mean": None,
            "edge_length_median": None,
            "edge_length_p90": None,
            "long_edge_share": None,
            "long_edge_factor": float(long_factor),
        }
    arr = np.array(lengths, float)
    med = float(np.median(arr))
    return {
        "edge_length_mean": float(arr.mean()),
        "edge_length_median": med,
        "edge_length_p90": float(np.percentile(arr, 90)),
        "long_edge_share": float(np.mean(arr > long_factor * max(med, 1e-9))),
        "long_edge_factor": float(long_factor),
        "edge_count": float(len(arr)),
    }


def long_edge_details(pathway: Pathway, pos: Positions, long_factor: float = 3.0, limit: int = 30) -> list[dict]:
    med_lengths = []
    edges = []
    for u, v in pathway.undirected().edges():
        if u in pos and v in pos:
            L = float(np.hypot(pos[u][0] - pos[v][0], pos[u][1] - pos[v][1]))
            med_lengths.append(L)
            edges.append((u, v, L))
    if not med_lengths:
        return []
    med = float(np.median(med_lengths))
    thr = long_factor * max(med, 1e-9)
    rows = []
    for u, v, L in sorted(edges, key=lambda t: -t[2]):
        if L < thr:
            break
        rows.append(
            dict(
                u=u,
                v=v,
                length=round(L, 1),
                relative_to_median=round(L / max(med, 1e-9), 2),
                mid_x=round(0.5 * (pos[u][0] + pos[v][0]), 1),
                mid_y=round(0.5 * (pos[u][1] + pos[v][1]), 1),
            )
        )
        if len(rows) >= limit:
            break
    return rows


def compactness_fallback(pos: Positions) -> dict[str, float | None]:
    """Aspect ratio and area per node without gdMetriX (axis-aligned bbox)."""
    if len(pos) < 2:
        return {"aspect_ratio_bbox": None, "area_per_node_bbox": None}
    xy = np.array(list(pos.values()), float)
    w = float(xy[:, 0].max() - xy[:, 0].min())
    h = float(xy[:, 1].max() - xy[:, 1].min())
    return {
        "aspect_ratio_bbox": (max(w, h) / max(min(w, h), 1e-9)),
        "area_per_node_bbox": (w * h) / len(pos),
    }


def all_benchmark_metrics(pathway: Pathway, pos: Positions, reference: Positions | None = None) -> dict[str, float | None]:
    m = general_metrics(pathway, pos)
    m.update(overlaps(pathway, pos))
    m["participant_consistency"] = participant_consistency(pathway, pos)
    m.update(compartment_metrics(pathway, pos))
    m.update(long_edge_metrics(pathway, pos))
    m.update(compactness_fallback(pos))
    if reference is not None:
        m.update(reference_similarity(pos, reference))
    return m
