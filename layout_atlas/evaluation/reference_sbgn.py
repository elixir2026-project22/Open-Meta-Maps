"""Hand-drawn SBGN reference positions (metabolismregulation gallery scheme).

Adapted from layout-benchmark `reference.py` (MIT) — only glyph matching + centres.
"""

from __future__ import annotations

from pathlib import Path
from xml.etree import ElementTree as ET

from .pathway import REACTION, Pathway

Positions = dict[str, tuple[float, float]]


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _glyph_for(node_id: str, glyphs) -> str | None:
    glyph = node_id.split("_", 1)[1] if "_" in node_id else node_id
    for candidate in (glyph, glyph.replace("_", "__")):
        if candidate in glyphs:
            return candidate
    return None


def read_sbgn_positions(path: str | Path) -> dict[str, tuple[float, float]]:
    root = ET.parse(path).getroot()
    pos = {}
    for el in root.iter():
        if _local(el.tag) != "glyph" or "id" not in el.attrib:
            continue
        bbox = next((c for c in el if _local(c.tag) == "bbox"), None)
        if bbox is None:
            continue
        x, y, w, h = (float(bbox.get(k)) for k in ("x", "y", "w", "h"))
        pos[el.get("id")] = (x + w / 2, -(y + h / 2))
    return pos


def read_sbgn_glyphs(path: str | Path) -> dict[str, tuple[str, str]]:
    root = ET.parse(path).getroot()
    out = {}
    for el in root.iter():
        if _local(el.tag) != "glyph" or not el.get("id") or not el.get("class"):
            continue
        label = next((c.get("text") for c in el if _local(c.tag) == "label"), None)
        out[el.get("id")] = (el.get("class"), label or "")
    return out


def _sbgn_processes(path: str | Path) -> dict[str, set[str]]:
    root = ET.parse(path).getroot()
    owner, kind = {}, {}
    for el in root.iter():
        if _local(el.tag) == "glyph" and el.get("id"):
            kind[el.get("id")] = el.get("class")
            for c in el:
                if _local(c.tag) == "port" and c.get("id"):
                    owner[c.get("id")] = el.get("id")
    links: dict[str, set[str]] = {}
    for el in root.iter():
        if _local(el.tag) != "arc":
            continue
        a, b = (owner.get(el.get(k), el.get(k)) for k in ("source", "target"))
        for x, y in ((a, b), (b, a)):
            if kind.get(y) == "process":
                links.setdefault(x, set()).add(y)
    return links


def match_glyphs(pathway: Pathway, sbgn: str | Path) -> dict[str, str]:
    glyphs = read_sbgn_glyphs(sbgn)
    by_scheme = {}
    for n in pathway.graph:
        glyph = n if n in glyphs else _glyph_for(n, glyphs)
        if glyph:
            by_scheme[n] = glyph
    if len(by_scheme) >= 0.5 * pathway.graph.number_of_nodes():
        return by_scheme
    match = {n: n for n, d in pathway.graph.nodes(data=True) if d["kind"] == REACTION and n in glyphs}
    processes = _sbgn_processes(sbgn)
    by_label: dict[str, list[str]] = {}
    for gid, (cls, label) in glyphs.items():
        if cls not in ("process", "compartment") and label:
            by_label.setdefault(label.strip().lower(), []).append(gid)
    used = set()
    for n, d in pathway.graph.nodes(data=True):
        if d["kind"] == REACTION:
            continue
        candidates = [g for g in by_label.get(d["label"].strip().lower(), []) if g not in used]
        if not candidates:
            continue
        reactions = {r for r in (set(pathway.graph.predecessors(n)) | set(pathway.graph.successors(n))) if r in match}
        best = max(candidates, key=lambda g: len(processes.get(g, set()) & reactions))
        match[n] = best
        used.add(best)
    return match


def gallery_reference(pathway: Pathway, sbgn: str | Path) -> Positions:
    glyphs = read_sbgn_positions(sbgn)
    return {n: glyphs[g] for n, g in match_glyphs(pathway, sbgn).items() if g in glyphs}
