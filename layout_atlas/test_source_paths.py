"""Regression tests for curated-curve reuse (``source_paths``).

Run: ``python3 test_source_paths.py``
"""
from __future__ import annotations

import json
import sqlite3
import time
import zlib

from core import DATA, ROOT, model_graph
import database
import ec_channel as ec
import model_io
import retrieval
import semantic_transfer as st
import source_paths as sp
from tool_pipeline import load_bridge

ECO_MAP_ID = "7c6b94063f4f96699811411b"


def _block_from_proposal(p: dict) -> dict:
    return dict(
        map_id=p["map_id"],
        node_mapping=p["node_mapping"],
        reaction_mapping=p["reaction_mapping"],
        positions=p["positions"],
        source=p.get("source"),
        reaction_ids=p.get("reaction_ids"),
    )


def _identity_block_for_reaction(data: dict, map_id: str, rkey: str, *, glyph_key_fn):
    """Place one source reaction at absolute coordinates; ``glyph_key_fn(gid)`` sets mapping keys."""
    rx = data["reactions"][rkey]
    node_mapping = {rkey: "r:" + rx["id"]}
    positions = {"r:" + rx["id"]: list(rx["xy"])}
    for gid in rx["chains"]:
        g = data["glyphs"].get(gid)
        if not g:
            continue
        tkey = glyph_key_fn(gid)
        tid = "m:" + gid
        node_mapping[tkey] = tid
        positions[tid] = list(g["xy"])
    # Pad with other map glyphs so the similarity fit has ≥3 correspondences.
    for gid, g in data["glyphs"].items():
        if len(node_mapping) >= 4:
            break
        if gid in rx["chains"] or g.get("kind") != "metabolite":
            continue
        tkey = glyph_key_fn(gid)
        if tkey in node_mapping:
            continue
        tid = "m:" + gid
        node_mapping[tkey] = tid
        positions[tid] = list(g["xy"])
    return dict(
        map_id=map_id,
        node_mapping=node_mapping,
        reaction_mapping={rkey: rx["id"]},
        positions=positions,
    )


def test_mapping_accepts_bare_and_prefixed_glyph_keys():
    """Escher modules store bare glyph ids; ``g:`` prefixes must also resolve."""
    data = sp.load()["maps"][ECO_MAP_ID]
    rkey, rx = next((k, r) for k, r in data["reactions"].items() if len(r["chains"]) >= 2)
    for glyph_key_fn in (lambda gid: gid, lambda gid: "g:" + gid):
        block = _identity_block_for_reaction(data, ECO_MAP_ID, rkey, glyph_key_fn=glyph_key_fn)
        paths, info = sp.paths_for_block(block, tolerance=5.0)
        assert info.get("too_few_correspondences", 0) == 0, info
        assert info.get("preserved", 0) >= 1, (info, list(paths)[:3])
        sample_gid = next(iter(rx["chains"]))
        assert sp._mapping_get(block["node_mapping"], sample_gid) is not None


def test_every_escher_block_reuses_curves_on_e_coli_core():
    """Selected Escher structures for e_coli_core must each reuse ≥1 curated chain."""
    t0 = time.time()
    model = model_io.load_model(ROOT / "examples/e_coli_core.xml")
    g = model_graph(model)
    core = database.modules()
    library = retrieval.Library(core)
    ps, _audit = library.candidates(g, model=model, max_roots=120)
    mapper = st.Mapper(model, load_bridge())
    strong, _sinfo = mapper.proposals(core, "open_core")
    selected, _stages = ec.select_with_ec(ps + strong, [])
    escher = [p for p in selected if p.get("source") == "escher"]
    assert escher, "expected at least one Escher structure on e_coli_core"
    failures = []
    for p in escher:
        paths, info = sp.paths_for_block(_block_from_proposal(p))
        preserved = int(info.get("preserved") or 0)
        if preserved < 1:
            failures.append(dict(
                map_id=p.get("map_id"), reactions=len(p.get("reaction_ids") or []),
                preserved=preserved, info={k: info[k] for k in info if info[k]},
                sample_keys=list(p["node_mapping"])[:8],
            ))
    assert not failures, failures
    print(f"PASS e_coli_core Escher curve reuse: {len(escher)} blocks, {time.time() - t0:.1f}s")


def test_bare_glyph_keys_on_e_coli_modules_reuse():
    """DB modules for the curated e_coli_core map reuse chains when placed at module coords."""
    con = sqlite3.connect(str(DATA / "open_core.sqlite"))
    checked = reused = 0
    for _id, payload in con.execute("select id, payload from modules where source='escher'"):
        m = json.loads(zlib.decompress(payload))
        if m.get("map_id") != ECO_MAP_ID or m.get("size", 0) < 2:
            continue
        node_mapping, positions, rmap = {}, {}, {}
        for n in m["nodes"]:
            tid = ("r:" + str(n["entity"])) if n["role"] == "reaction" else ("m:" + str(n["key"]))
            node_mapping[n["key"]] = tid
            positions[tid] = n["xy"]
        for r in m["reactions"]:
            if r["key"] in node_mapping:
                rmap[r["key"]] = r["id"]
        if len(set(node_mapping.values())) < 3:
            continue
        checked += 1
        _paths, info = sp.paths_for_block(dict(
            map_id=m["map_id"], node_mapping=node_mapping,
            reaction_mapping=rmap, positions=positions,
        ))
        if info.get("too_few_correspondences"):
            continue
        assert info.get("preserved", 0) >= 1, (m["id"], info)
        reused += 1
    assert checked >= 50 and reused >= 50, (checked, reused)
    print(f"PASS e_coli module curve reuse: reused={reused} checked={checked}")


if __name__ == "__main__":
    test_mapping_accepts_bare_and_prefixed_glyph_keys()
    test_bare_glyph_keys_on_e_coli_modules_reuse()
    test_every_escher_block_reuses_curves_on_e_coli_core()
    print("ok")
