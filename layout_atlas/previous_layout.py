"""Build-on-previous layout: pin reaction/metabolite positions from a prior map.json.

When enough reactions are already placed (≥ ``PACK_SKIP_SHARE``), module packing is
skipped and only free nodes are completed around the pinned anchors.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from curated_axes import _absolute_from_escher

# Skip block packing when at least this share of model reactions have a previous pin.
PACK_SKIP_SHARE = 0.80


def load_previous_anchors(path, model=None):
    """Read reaction midmarkers and primary metabolites from an Escher ``map.json``.

    Returns ``{ 'r:<id>': [x, y], 'm:<id>': [x, y] }``.  When ``model`` is given,
    only identifiers present in that model are kept.
    """
    path = Path(path)
    payload = json.loads(path.read_text(encoding='utf-8'))
    body = payload[1] if isinstance(payload, list) and len(payload) >= 2 else payload
    if not isinstance(body, dict):
        raise ValueError(f'Not an Escher map: {path}')
    abs_map = _absolute_from_escher(body, str(path))
    anchors = {}
    for rid, xy in (abs_map.get('reactions') or {}).items():
        anchors['r:' + rid] = [float(xy[0]), float(xy[1])]
    for mid, xy in (abs_map.get('metabolites') or {}).items():
        anchors['m:' + mid] = [float(xy[0]), float(xy[1])]
    if model is not None:
        rx = {r['id'] for r in model.get('reactions') or []}
        mets = {m['id'] for m in model.get('metabolites') or []}
        anchors = {
            k: v for k, v in anchors.items()
            if (k.startswith('r:') and k[2:] in rx) or (k.startswith('m:') and k[2:] in mets)
        }
    return anchors


def reaction_anchor_share(anchors, model):
    """Share of model reactions that have a previous midmarker pin."""
    rx = model.get('reactions') or []
    if not rx:
        return 0.0
    hits = sum(1 for r in rx if ('r:' + r['id']) in (anchors or {}))
    return hits / len(rx)


def should_skip_packing(anchors, model, threshold=PACK_SKIP_SHARE):
    return reaction_anchor_share(anchors, model) >= float(threshold)


def apply_block_positions(blocks, anchors):
    """Rewrite curated block glyph positions from previous anchors where present."""
    if not anchors:
        return blocks, 0
    n = 0
    for b in blocks or []:
        pos = dict(b.get('positions') or {})
        for k in list(pos):
            if k in anchors:
                pos[k] = list(anchors[k])
                n += 1
        b['positions'] = pos
    return blocks, n


def pin_scene(scene, anchors):
    """Hard-pin previous anchors into a display scene before export.

    Pins ``r:`` / ``m:`` keys and curated ``copy:…:met`` glyphs to the previous
    primary metabolite position.  Currency ``side:`` stubs follow their reaction
    midmarker by the same delta so spokes stay attached.  Must run before
    ``drawing.export`` (which applies one uniform canvas translation); do not
    re-apply absolute previous coordinates onto the exported Escher JSON.
    """
    if not scene or not anchors:
        return scene, 0
    pos = scene.get('positions') or {}
    n = 0
    old_mids = {k: list(pos[k]) for k in pos if k.startswith('r:') and k in anchors}
    for k, xy in anchors.items():
        if k in pos:
            pos[k] = [float(xy[0]), float(xy[1])]
            n += 1
    for k in list(pos):
        if not k.startswith('copy:'):
            continue
        met = k.rsplit(':', 1)[-1]
        mk = 'm:' + met
        if mk in anchors:
            pos[k] = [float(anchors[mk][0]), float(anchors[mk][1])]
            n += 1
    # Keep side stubs attached to pinned midmarkers.
    for k in list(pos):
        if not k.startswith('side:'):
            continue
        # side:<rid>:<met>
        parts = k.split(':', 2)
        if len(parts) < 3:
            continue
        rk = 'r:' + parts[1]
        if rk not in anchors or rk not in old_mids:
            continue
        delta = np.asarray(anchors[rk], float) - np.asarray(old_mids[rk], float)
        pos[k] = (np.asarray(pos[k], float) + delta).tolist()
        n += 1
    scene['positions'] = pos
    return scene, n


def pin_escher(data, anchors):
    """No-op kept for call-site compatibility; scene pinning + export translation is enough.

    Re-applying absolute previous coordinates after ``drawing.export`` mixes the
    previous canvas frame with the new export translation and scrambles edges.
    """
    return data, 0


def filter_previous_map(path, model):
    """Keep the previous Escher drawing for reactions still in ``model``.

    Drops removed reactions and any nodes that no remaining segment references.
    New model reactions are listed in the report for later placement; when the
    model only lost reactions, the filtered map is the previous layout itself.
    """
    from core import uid
    path = Path(path)
    payload = json.loads(path.read_text(encoding='utf-8'))
    if not (isinstance(payload, list) and len(payload) >= 2 and isinstance(payload[1], dict)):
        raise ValueError(f'Not an Escher map: {path}')
    header = dict(payload[0] or {})
    body = dict(payload[1])
    keep_rx = {r['id'] for r in model.get('reactions') or []}
    reactions = body.get('reactions') or {}
    nodes = body.get('nodes') or {}
    kept, used_nodes = {}, set()
    for rid, rx in reactions.items():
        bigg = rx.get('bigg_id')
        if bigg not in keep_rx:
            continue
        kept[rid] = rx
        for seg in (rx.get('segments') or {}).values():
            for nid in (seg.get('from_node_id'), seg.get('to_node_id')):
                if nid is not None:
                    used_nodes.add(str(nid))
    kept_nodes = {nid: dict(nodes[nid]) for nid in used_nodes if nid in nodes}
    # Keep the previous compare-map representative glyph for each remaining metabolite
    # (same rule as evaluation.escher_bridge.escher_positions: last primary wins).  Do not
    # move other clones — that would stretch edges that still use them.
    ref_nid = {}
    for nid, n in nodes.items():
        if n.get('node_type') != 'metabolite' or not n.get('bigg_id'):
            continue
        mid = n['bigg_id']
        if mid not in ref_nid or n.get('node_is_primary'):
            ref_nid[mid] = nid
    present = {n['bigg_id'] for n in kept_nodes.values()
               if n.get('node_type') == 'metabolite' and n.get('bigg_id')}
    for mid, nid in ref_nid.items():
        if mid not in present:
            continue
        if nid not in kept_nodes and nid in nodes:
            kept_nodes[nid] = dict(nodes[nid])
    for nid, n in kept_nodes.items():
        if n.get('node_type') != 'metabolite' or not n.get('bigg_id'):
            continue
        n['node_is_primary'] = (ref_nid.get(n['bigg_id']) == nid)
    drawn = {rx.get('bigg_id') for rx in kept.values() if rx.get('bigg_id')}
    new_reactions = sorted(keep_rx - drawn)
    meta = dict(header.get('metadata') or {})
    meta['model_hash'] = uid(model)
    meta['previous_layout'] = dict(
        source=str(path), kept_reactions=len(kept), removed_reactions=len(reactions) - len(kept),
        new_reactions=new_reactions, nodes=len(kept_nodes))
    header['metadata'] = meta
    out = [header, dict(body, nodes=kept_nodes, reactions=kept)]
    report = dict(kept_reactions=len(kept), removed_reactions=len(reactions) - len(kept),
                  new_reactions=new_reactions, nodes=len(kept_nodes))
    return out, report


def missing_model_reactions(path, model):
    """Reaction ids in ``model`` that have no midmarker in the previous map."""
    _, report = filter_previous_map(path, model)
    return list(report.get('new_reactions') or [])
