"""Localise weak metric regions: reaction ids and map windows for the report."""

from __future__ import annotations

from collections import defaultdict
from typing import Any

import numpy as np

from .metrics import long_edge_details, participant_consistency_details
from .pathway import Pathway

Positions = dict[str, tuple[float, float]]


def crossing_hotspots_from_issues(issues: dict | None, limit: int = 8) -> list[dict[str, Any]]:
    """Use existing layout_quality issues (curve crossings) when available."""
    if not issues:
        return []
    regions = list(issues.get("worst_regions") or [])[:limit]
    crossings = list(issues.get("crossings") or [])
    # Aggregate reactions that appear most often in algorithm crossings
    counts = defaultdict(int)
    for c in crossings:
        if c.get("kind") in ("expert_expert",):
            continue
        for rid in c.get("reactions") or []:
            counts[rid] += 1
    top = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[:limit]
    return {
        "worst_regions": regions,
        "reactions_with_most_crossings": [dict(reaction=r, crossings=n) for r, n in top],
    }


def weak_spots(pathway: Pathway, pos: Positions, issues: dict | None = None, limit: int = 12) -> dict[str, Any]:
    """Concrete examples: inconsistent reactions, long edges, crossing hotspots."""
    inconsistent = [r for r in participant_consistency_details(pathway, pos, limit=limit) if r["opposite_share"] < 0.5]
    long_edges = long_edge_details(pathway, pos, limit=limit)
    # Attach reaction ids when an endpoint is a reaction
    for e in long_edges:
        ends = []
        for n in (e["u"], e["v"]):
            if n in pathway.graph and pathway.graph.nodes[n].get("kind") == "reaction":
                ends.append(n)
        e["reaction_ends"] = ends
    hot = crossing_hotspots_from_issues(issues, limit=limit)
    # Dense windows from node coordinates (400-unit style, independent of issues)
    windows = densest_windows(pos, pathway, limit=5)
    return dict(
        inconsistent_reactions=inconsistent,
        long_edges=long_edges,
        crossing_hotspots=hot,
        dense_windows=windows,
    )


def densest_windows(pos: Positions, pathway: Pathway, cell: float = 400.0, limit: int = 5) -> list[dict]:
    if len(pos) < 4:
        return []
    xy = np.array(list(pos.values()), float)
    keys = list(pos)
    x0, y0 = float(xy[:, 0].min()), float(xy[:, 1].min())
    cells = defaultdict(list)
    for n, (x, y) in pos.items():
        cells[(int((x - x0) // cell), int((y - y0) // cell))].append(n)
    ranked = sorted(cells.items(), key=lambda kv: -len(kv[1]))
    out = []
    for (gx, gy), members in ranked:
        reactions = [n for n in members if pathway.graph.nodes.get(n, {}).get("kind") == "reaction"]
        out.append(
            dict(
                x=round(x0 + gx * cell, 1),
                y=round(y0 + gy * cell, 1),
                w=cell,
                h=cell,
                nodes=len(members),
                reactions=sorted(reactions)[:12],
            )
        )
        if len(out) >= limit:
            break
    return out
