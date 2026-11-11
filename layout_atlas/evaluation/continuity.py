"""Structure continuity on the final Escher map.

For each reused expert structure (or, as a fallback, each connected component of
expert-derived reactions linked by shared non-currency metabolites), measure how
contiguous it stays: share of reaction pairs that share a non-currency metabolite
and are drawn within a few median reaction spacings of each other, plus how many
structures split into separate geometric pieces.

Scale-free: adjacency uses k × median nearest-neighbour reaction spacing of the
whole map (default k=3). No model-specific thresholds.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Any

import numpy as np
from scipy.spatial import cKDTree

try:
    from core import CURRENCY, chem
except ImportError:  # pragma: no cover
    CURRENCY, chem = set(), lambda mid: mid


ADJACENCY_FACTOR = 3.0


def _is_currency(mid: str) -> bool:
    return chem(mid) in CURRENCY


def _reaction_positions(body: dict) -> dict[str, tuple[float, float]]:
    nodes = body["nodes"]
    out = {}
    for rx in body["reactions"].values():
        mid = None
        for s in rx["segments"].values():
            for nid in (s["from_node_id"], s["to_node_id"]):
                if nodes.get(nid, {}).get("node_type") == "midmarker":
                    mid = nid
                    break
            if mid is not None:
                break
        if mid is not None:
            out[rx["bigg_id"]] = (float(nodes[mid]["x"]), float(nodes[mid]["y"]))
    return out


def _reaction_mets(body: dict) -> dict[str, set[str]]:
    out = {}
    for rx in body["reactions"].values():
        mets = set()
        for m in rx.get("metabolites") or []:
            if isinstance(m, dict) and m.get("bigg_id"):
                mets.add(m["bigg_id"])
        out[rx["bigg_id"]] = mets
    return out


def _spacing(pos: dict[str, tuple[float, float]]) -> float:
    if len(pos) < 2:
        return 60.0
    xy = np.array(list(pos.values()), float)
    d = cKDTree(xy).query(xy, k=2)[0][:, 1]
    d = d[d > 0]
    return float(np.median(d)) if len(d) else 60.0


def _union_find(nodes: list[str]):
    parent = {n: n for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        a, b = find(a), find(b)
        if a != b:
            parent[a] = b

    return find, union


def structures_from_map(data: list | dict, structures: list[dict] | dict | None = None) -> dict[str, list[str]]:
    """id -> reaction ids. Prefers explicit structures, then metadata, then expert CC fallback."""
    body = data[1] if isinstance(data, list) else data
    meta = (data[0].get("metadata") if isinstance(data, list) else {}) or {}
    pos = _reaction_positions(body)

    groups: dict[str, list[str]] = {}
    raw = structures if structures is not None else meta.get("expert_structures")
    if isinstance(raw, dict):
        for sid, rids in raw.items():
            placed = [r for r in rids if r in pos]
            if len(placed) >= 2:
                groups[str(sid)] = placed
    elif isinstance(raw, list):
        for row in raw:
            if not isinstance(row, dict):
                continue
            sid = row.get("id") or row.get("code") or row.get("module_id")
            rids = row.get("reactions") or row.get("reaction_ids") or []
            placed = [r for r in rids if r in pos]
            if sid and len(placed) >= 2:
                groups[str(sid)] = placed
    if groups:
        return groups

    # Fallback: connected components of expert-derived reactions via non-currency metabolites
    prov = meta.get("reaction_provenance") or {}
    if isinstance(prov, dict):
        expert = [r for r, kind in prov.items() if kind == "expert_derived" and r in pos]
    else:
        expert = []
    if len(expert) < 2:
        return {}
    mets = _reaction_mets(body)
    find, union = _union_find(expert)
    for i, a in enumerate(expert):
        for b in expert[i + 1 :]:
            shared = {m for m in mets.get(a, set()) & mets.get(b, set()) if not _is_currency(m)}
            if shared:
                union(a, b)
    comps: dict[str, list[str]] = defaultdict(list)
    for r in expert:
        comps[find(r)].append(r)
    return {f"expert_cc_{i}": rs for i, rs in enumerate(comps.values()) if len(rs) >= 2}


def structure_continuity(
    data: list | dict,
    structures: list[dict] | dict | None = None,
    adjacency_factor: float = ADJACENCY_FACTOR,
) -> dict[str, Any]:
    """Return continuity summary for ``layout_quality`` core / evaluation reports."""
    body = data[1] if isinstance(data, list) else data
    pos = _reaction_positions(body)
    mets = _reaction_mets(body)
    groups = structures_from_map(data, structures)
    spacing = _spacing(pos)
    thr = adjacency_factor * max(spacing, 1e-9)

    pairs = contiguous = split = 0
    details = []
    for sid, rids in groups.items():
        placed = [r for r in rids if r in pos]
        if len(placed) < 2:
            continue
        edges = []
        for i, a in enumerate(placed):
            for b in placed[i + 1 :]:
                shared = {m for m in mets.get(a, set()) & mets.get(b, set()) if not _is_currency(m)}
                if not shared:
                    continue
                dist = float(np.hypot(pos[a][0] - pos[b][0], pos[a][1] - pos[b][1]))
                edges.append((a, b, dist))
        n_pairs = len(edges)
        n_contig = sum(1 for _, _, d in edges if d <= thr)
        pairs += n_pairs
        contiguous += n_contig
        find, union = _union_find(placed)
        for a, b, d in edges:
            if d <= thr:
                union(a, b)
        components = len({find(r) for r in placed})
        if components > 1:
            split += 1
        details.append(
            dict(
                structure=sid,
                reactions=len(placed),
                pairs=n_pairs,
                contiguous_pairs=n_contig,
                components=components,
                contiguous_pair_share=round(n_contig / n_pairs, 3) if n_pairs else None,
            )
        )

    structures_n = len(details)
    return dict(
        contiguous_pair_share=round(contiguous / pairs, 3) if pairs else None,
        pairs=pairs,
        contiguous_pairs=contiguous,
        structures=structures_n,
        split_structures=split,
        split_structure_share=round(split / structures_n, 3) if structures_n else None,
        reaction_spacing_px=round(spacing, 1),
        adjacency_factor=float(adjacency_factor),
        note=(
            "Share of expert-structure reaction pairs that share a non-currency metabolite "
            "and lie within adjacency_factor × median reaction spacing; split_structures count "
            "structures whose contiguous pairs form more than one geometric piece."
        ),
        details=details[:40],
    )
