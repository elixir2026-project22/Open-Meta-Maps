"""Read an Escher map into pathway node positions (and optional clone-aware drawing).

Adapted from layout-benchmark `external.escher_drawing` (MIT).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
from scipy.spatial import cKDTree

from .pathway import REACTION, SPECIES, Pathway

Positions = dict[str, tuple[float, float]]
SPACING = 1.5


def _id_lookup(pathway: Pathway) -> dict[str, str]:
    lookup = {}
    for n in pathway.graph:
        lookup[n] = n
        for prefix in ("M_", "R_", "S_", "s_", "r_"):
            if n.startswith(prefix):
                lookup.setdefault(n[len(prefix) :], n)
        # BiGG / Escher often drops compartment suffixes inconsistently; keep bare id too
        if "_" in n and n.rsplit("_", 1)[-1] in ("c", "e", "p", "m", "x", "n", "g", "r", "v", "l"):
            lookup.setdefault(n.rsplit("_", 1)[0], n)
    return lookup


def _body(escher_map: list | dict) -> dict:
    return escher_map[1] if isinstance(escher_map, list) else escher_map


def load_escher(path: str | Path) -> list | dict:
    return json.loads(Path(path).read_text())


def escher_positions(pathway: Pathway, escher_map: list | dict, y_up: bool = True) -> Positions:
    """One position per pathway node present in the map (primary metabolite / midmarker).

    Escher y grows downward; with ``y_up`` the y axis is flipped so metrics match the
    layout-benchmark convention.
    """
    body = _body(escher_map)
    nodes = body["nodes"]
    lookup = _id_lookup(pathway)
    pos: Positions = {}

    def xy(node):
        x, y = float(node["x"]), float(node["y"])
        return (x, -y) if y_up else (x, y)

    # metabolites: prefer primary clones
    for node in nodes.values():
        if node.get("node_type") != "metabolite":
            continue
        sid = lookup.get(node.get("bigg_id"))
        if sid is None:
            continue
        if sid not in pos or node.get("node_is_primary"):
            pos[sid] = xy(node)

    for reaction in body["reactions"].values():
        rid = lookup.get(reaction.get("bigg_id"))
        if rid is None or pathway.graph.nodes.get(rid, {}).get("kind") != REACTION:
            continue
        seg = nx.Graph()
        for s in reaction["segments"].values():
            seg.add_edge(s["from_node_id"], s["to_node_id"])
        mid = next((n for n in seg if nodes.get(n, {}).get("node_type") == "midmarker"), None)
        if mid is None:
            continue
        pos[rid] = xy(nodes[mid])
    return pos


def escher_drawing_graph(pathway: Pathway, escher_map: list | dict, y_up: bool = True):
    """Clone-aware drawn graph + positions (for overlap boxes that follow Escher clones)."""
    body = _body(escher_map)
    nodes = body["nodes"]
    lookup = _id_lookup(pathway)
    g = nx.DiGraph()
    pos: Positions = {}

    def xy(node):
        x, y = float(node["x"]), float(node["y"])
        return (x, -y) if y_up else (x, y)

    drawn = {}
    for nid, node in nodes.items():
        if node.get("node_type") != "metabolite" or node.get("bigg_id") not in lookup:
            continue
        sid = lookup[node["bigg_id"]]
        primary = bool(node.get("node_is_primary"))
        d = pathway.graph.nodes[sid]
        scale = 1.0 if primary else 0.6
        did = f"{sid}@{nid}"
        g.add_node(
            did,
            of=sid,
            kind=SPECIES,
            label=d["label"],
            width=d["width"] * scale,
            height=d["height"] * scale,
            compartment=d.get("compartment"),
            primary=primary,
        )
        pos[did] = xy(node)
        drawn[nid] = did

    for reaction in body["reactions"].values():
        rid = lookup.get(reaction.get("bigg_id"))
        if rid is None or pathway.graph.nodes.get(rid, {}).get("kind") != REACTION:
            continue
        seg = nx.Graph()
        for s in reaction["segments"].values():
            seg.add_edge(s["from_node_id"], s["to_node_id"])
        mid = next((n for n in seg if nodes.get(n, {}).get("node_type") == "midmarker"), None)
        if mid is None:
            continue
        d = pathway.graph.nodes[rid]
        g.add_node(
            rid,
            of=rid,
            kind=REACTION,
            label=d["label"],
            width=d["width"],
            height=d["height"],
            compartment=d.get("compartment"),
        )
        pos[rid] = xy(nodes[mid])
        coefficient = {m["bigg_id"]: m["coefficient"] for m in reaction.get("metabolites", [])}
        for nid in seg:
            if nid not in drawn:
                continue
            bigg = nodes[nid]["bigg_id"]
            if coefficient.get(bigg, 0) > 0:
                g.add_edge(rid, drawn[nid], role="product")
            else:
                g.add_edge(drawn[nid], rid, role="reactant")
    return g, pos


def normalise_positions(pos: Positions, widths: dict[str, float] | None = None) -> Positions:
    """Uniform scale so median nearest-neighbour distance is SPACING median widths (shape-preserving)."""
    if len(pos) < 2:
        return dict(pos)
    nodes = list(pos)
    xy = np.array([pos[n] for n in nodes], dtype=float)
    centre = xy.mean(axis=0)
    distances = cKDTree(xy).query(xy, k=2)[0][:, 1]
    distances = distances[distances > 0]
    if widths:
        wvals = [widths[n] for n in nodes if widths.get(n)]
        width = float(np.median(wvals)) if wvals else 1.0
    else:
        width = 1.0
    if not len(distances):
        scale = 1.0
    else:
        scale = SPACING * width / float(np.median(distances))
    return {n: (float((pos[n][0] - centre[0]) * scale), float((pos[n][1] - centre[1]) * scale)) for n in nodes}


def coverage(pathway: Pathway, pos: Positions) -> dict[str, float]:
    placed = set(pos)
    return {
        "nodes_total": float(pathway.graph.number_of_nodes()),
        "nodes_placed": float(len(placed & set(pathway.graph))),
        "reactions_placed": float(len(placed & set(pathway.reactions))),
        "species_placed": float(len(placed & set(pathway.species))),
        "placement_share": float(len(placed & set(pathway.graph)) / max(pathway.graph.number_of_nodes(), 1)),
    }
