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
# Keep a local relayout of free (new) nodes only when crossings fall by this much.
CROSSING_IMPROVE_SHARE = 0.25
CROSSING_IMPROVE_MIN = 2
# Offset for parking a leftover reaction that has no drawn participant anchors.
NEW_COMPONENT_GAP = 120.0


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
    # Keep source_reactions only for reactions still in the model (graft fills new ones later).
    src = dict(meta.get('source_reactions') or {})
    rx_by_id = {r['id']: r for r in model.get('reactions') or []}
    pruned = {}
    for rid in keep_rx:
        if rid in src:
            pruned[rid] = src[rid]
        elif rid in drawn:
            r = rx_by_id[rid]
            pruned[rid] = dict(
                stoichiometry=dict(r.get('metabolites') or {}),
                lower_bound=r.get('lower_bound', 0),
                upper_bound=r.get('upper_bound', 0))
    meta['source_reactions'] = pruned
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


def unanchored_reaction_ids(model, anchors):
    """Model reaction ids with no previous midmarker pin (matched by stable id)."""
    anchors = anchors or {}
    return [r['id'] for r in model.get('reactions') or [] if ('r:' + r['id']) not in anchors]


def keep_local_relayout(crossings_before, crossings_after,
                        share=CROSSING_IMPROVE_SHARE, vmin=CROSSING_IMPROVE_MIN):
    """True only when a local relax cuts crossings enough to justify moving free nodes."""
    if crossings_before is None or crossings_after is None:
        return False
    before, after = int(crossings_before), int(crossings_after)
    drop = before - after
    return drop >= int(vmin) and drop >= float(share) * max(before, 1)


def scene_segment_crossings(model, scene):
    """Straight-segment crossing count for the current display scene (cheap gate signal)."""
    import drawing
    pos = scene.get('positions') or {}
    parts = scene.get('participants') or {}
    segments = []
    for r in model.get('reactions') or []:
        rk = 'r:' + r['id']
        if rk not in pos:
            continue
        for mid in r.get('metabolites') or {}:
            key = parts.get((r['id'], mid))
            if key in pos:
                segments.append((rk, key))
    if len(segments) >= 1200:
        return None
    return int(drawing.crossings(pos, segments) or 0)


def graft_new_reactions(data, model, new_reaction_ids, min_embed_anchors=2):
    """Add missing model reactions onto a filtered previous Escher map (pin-old / place-new).

    Matched geometry is left untouched.  New reactions with enough drawn participants are
    embedded beside them; otherwise they are parked as a new component.  Never skips a
    leftover — coverage stays complete for ``assert_chemistry``.
    """
    from core import currency as model_currency, uid
    new_reaction_ids = list(new_reaction_ids or [])
    if not new_reaction_ids:
        return data, dict(placed=0, embedded=0, new_components=0, leftovers=[])
    header, body = data[0], data[1]
    nodes = dict(body.get('nodes') or {})
    reactions = dict(body.get('reactions') or {})
    cur = set(model_currency(model))
    rx_by_id = {r['id']: r for r in model.get('reactions') or []}
    # Metabolite bigg_id -> node id (prefer primary).
    met_nodes = {}
    for nid, n in nodes.items():
        if n.get('node_type') != 'metabolite' or not n.get('bigg_id'):
            continue
        mid = n['bigg_id']
        if mid not in met_nodes or n.get('node_is_primary'):
            met_nodes[mid] = str(nid)
    xs = [float(n['x']) for n in nodes.values() if 'x' in n]
    ys = [float(n['y']) for n in nodes.values() if 'y' in n]
    if xs and ys:
        lo = np.array([min(xs), min(ys)], float)
        hi = np.array([max(xs), max(ys)], float)
    else:
        lo = np.array([0.0, 0.0]); hi = np.array([0.0, 0.0])
    park = np.array([hi[0] + NEW_COMPONENT_GAP, lo[1]], float)
    next_id = 0
    def fresh(prefix='n'):
        nonlocal next_id
        while True:
            nid = f'{prefix}{next_id}'
            next_id += 1
            if nid not in nodes and nid not in reactions:
                return nid
    embedded = new_components = 0
    for i, rid in enumerate(new_reaction_ids):
        r = rx_by_id[rid]
        core = [m for m in r['metabolites'] if m not in cur]
        anchor_xy = []
        for mid in core:
            nid = met_nodes.get(mid)
            if nid and nid in nodes:
                anchor_xy.append((mid, nid, np.array([nodes[nid]['x'], nodes[nid]['y']], float)))
        if len(anchor_xy) >= max(1, int(min_embed_anchors)):
            target = np.mean([p for _, _, p in anchor_xy], axis=0)
            direction = np.array([0.0, 1.0])
            if len(anchor_xy) >= 2:
                delta = anchor_xy[0][2] - anchor_xy[1][2]
                nrm = float(np.hypot(*delta))
                if nrm > 1e-6:
                    direction = np.array([-delta[1], delta[0]], float) / nrm
            mid_xy = target + direction * 55.0
            embedded += 1
        elif anchor_xy:
            mid_xy = anchor_xy[0][2] + np.array([0.0, 70.0])
            new_components += 1
        else:
            mid_xy = park + np.array([(i % 6) * 90.0, (i // 6) * 90.0])
            new_components += 1
        mid_nid = fresh('r')
        nodes[mid_nid] = dict(node_type='midmarker', x=float(mid_xy[0]), y=float(mid_xy[1]))
        sign = -1 if r.get('upper_bound', 0) <= 0 and r.get('lower_bound', 0) < 0 else 1
        segs = {}
        sid = 0
        for j, (mid, coeff) in enumerate(r['metabolites'].items()):
            if mid in met_nodes and met_nodes[mid] in nodes:
                met_nid = met_nodes[mid]
            else:
                met_nid = fresh('m')
                offset = np.array([30.0 * (j % 3 - 1), 50.0 + 20.0 * (j // 3)])
                xy = mid_xy + offset
                nodes[met_nid] = dict(
                    node_type='metabolite', bigg_id=mid, name=mid,
                    node_is_primary=mid not in cur,
                    x=float(xy[0]), y=float(xy[1]),
                    label_x=float(xy[0] + 10), label_y=float(xy[1] - 10))
                if mid not in cur:
                    met_nodes[mid] = met_nid
            segs[str(sid)] = dict(from_node_id=mid_nid, to_node_id=met_nid, b1=None, b2=None)
            sid += 1
        reactions[fresh('rx')] = dict(
            bigg_id=rid, name=r.get('name') or rid,
            reversibility=bool(r.get('lower_bound', 0) < 0 < r.get('upper_bound', 0)),
            label_x=float(mid_xy[0] + 10), label_y=float(mid_xy[1] - 10),
            gene_reaction_rule=r.get('gene_reaction_rule', ''), genes=[],
            metabolites=[dict(bigg_id=m, coefficient=c * sign) for m, c in r['metabolites'].items()],
            segments=segs)
    # Expand canvas to include grafted nodes.
    canvas = dict(body.get('canvas') or {})
    xs = [float(n['x']) for n in nodes.values() if 'x' in n]
    ys = [float(n['y']) for n in nodes.values() if 'y' in n]
    if xs and ys:
        pad = 40.0
        canvas = dict(
            x=float(min(xs) - pad), y=float(min(ys) - pad),
            width=float(max(xs) - min(xs) + 2 * pad),
            height=float(max(ys) - min(ys) + 2 * pad))
    meta = dict(header.get('metadata') or {})
    src = dict(meta.get('source_reactions') or {})
    for rid in new_reaction_ids:
        r = rx_by_id[rid]
        src[rid] = dict(
            stoichiometry=r['metabolites'],
            lower_bound=r.get('lower_bound', 0),
            upper_bound=r.get('upper_bound', 0))
    # Drop source_reactions entries for reactions no longer in the model.
    keep_rx = {r['id'] for r in model.get('reactions') or []}
    src = {k: v for k, v in src.items() if k in keep_rx}
    # Ensure every kept reaction has source metadata (filter path may already have them).
    for r in model.get('reactions') or []:
        src.setdefault(r['id'], dict(
            stoichiometry=r['metabolites'],
            lower_bound=r.get('lower_bound', 0),
            upper_bound=r.get('upper_bound', 0)))
    meta['source_reactions'] = src
    meta['model_hash'] = uid(model)
    meta['previous_layout'] = dict(meta.get('previous_layout') or {}, grafted=new_reaction_ids)
    header = dict(header, metadata=meta)
    out = [header, dict(body, nodes=nodes, reactions=reactions, canvas=canvas)]
    return out, dict(
        placed=len(new_reaction_ids), embedded=embedded, new_components=new_components,
        leftovers=list(new_reaction_ids))


def place_leftover_reactions(model, scene, anchors=None, cur=None, min_embed_anchors=2):
    """Pin-old / place-new: park every unanchored reaction so coverage stays complete.

    Reactions with ≥ ``min_embed_anchors`` drawn participant metabolites are embedded
    beside those anchors.  Otherwise the reaction (and any private metabolites) form a
    new component near the canvas — never dropped for lack of anchors.
    """
    from core import currency as model_currency
    if scene is None:
        return scene, dict(placed=0, embedded=0, new_components=0, leftovers=[])
    anchors = anchors or {}
    cur = set(cur if cur is not None else model_currency(model))
    scene = dict(scene, positions={k: list(v) for k, v in (scene.get('positions') or {}).items()})
    pos = scene['positions']
    parts = scene.setdefault('participants', {})
    copies = scene.setdefault('copies', {})
    labels = scene.setdefault('labels', {})
    leftovers = unanchored_reaction_ids(model, anchors)
    if not leftovers:
        # Still ensure every model reaction midmarker exists (defensive).
        leftovers = [r['id'] for r in model.get('reactions') or [] if ('r:' + r['id']) not in pos]
    if not leftovers:
        return scene, dict(placed=0, embedded=0, new_components=0, leftovers=[])

    drawn_met_keys = {}
    for key, mid in copies.items():
        if mid in cur:
            continue
        if key.startswith('side:'):
            continue
        drawn_met_keys.setdefault(mid, []).append(key)

    pts = [v for k, v in pos.items() if not k.startswith('side:')]
    if pts:
        arr = np.asarray(pts, float)
        lo, hi = arr.min(0), arr.max(0)
    else:
        lo, hi = np.array([0.0, 0.0]), np.array([0.0, 0.0])
    park_origin = np.array([hi[0] + NEW_COMPONENT_GAP, lo[1]], float)
    embedded = new_components = 0

    rx_by_id = {r['id']: r for r in model.get('reactions') or []}
    for i, rid in enumerate(leftovers):
        r = rx_by_id[rid]
        rk = 'r:' + rid
        core = [m for m in r['metabolites'] if m not in cur]
        anchor_xy = []
        for mid in core:
            keys = [k for k in drawn_met_keys.get(mid, []) if k in pos]
            if not keys:
                mk = 'm:' + mid
                if mk in pos:
                    keys = [mk]
            if keys:
                key = min(keys, key=lambda k: (0 if k.startswith('m:') else 1, k))
                anchor_xy.append((mid, key, np.asarray(pos[key], float)))
        if len(anchor_xy) >= max(1, int(min_embed_anchors)):
            target = np.mean([p for _, _, p in anchor_xy], axis=0)
            direction = np.array([0.0, 1.0])
            if len(anchor_xy) >= 2:
                delta = anchor_xy[0][2] - anchor_xy[1][2]
                n = float(np.hypot(*delta))
                if n > 1e-6:
                    direction = np.array([-delta[1], delta[0]], float) / n
            goal = target + direction * 55.0
            embedded += 1
        elif anchor_xy:
            goal = anchor_xy[0][2] + np.array([0.0, 70.0])
            new_components += 1
        else:
            goal = park_origin + np.array([(i % 6) * 90.0, (i // 6) * 90.0])
            new_components += 1
        old = np.asarray(pos.get(rk, goal), float)
        pos[rk] = [float(goal[0]), float(goal[1])]
        labels[rk] = rid
        for j, mid in enumerate(r['metabolites']):
            key = parts.get((rid, mid))
            if key is None:
                if mid in cur:
                    key = 'side:' + rid + ':' + mid
                elif mid in drawn_met_keys and any(k in pos for k in drawn_met_keys[mid]):
                    key = next(k for k in drawn_met_keys[mid] if k in pos)
                else:
                    key = 'm:' + mid
                parts[rid, mid] = key
                copies[key] = mid
                labels.setdefault(key, mid)
            if key.startswith('side:'):
                pos[key] = [float(goal[0] + 20.0 * (j % 3 - 1)), float(goal[1] + 40.0)]
                continue
            if key in pos and mid in drawn_met_keys and key in drawn_met_keys.get(mid, []):
                # Already-drawn participant — attach only; do not move the pin.
                continue
            delta = np.asarray(pos.get(key, old + np.array([0.0, 80.0])), float) - old
            if float(np.hypot(*delta)) < 1e-6:
                delta = np.array([0.0, 80.0 + 15.0 * j])
            pos[key] = (goal + delta).tolist()
            drawn_met_keys.setdefault(mid, []).append(key)

    scene['positions'] = pos
    return scene, dict(
        placed=len(leftovers), embedded=embedded, new_components=new_components,
        leftovers=list(leftovers))
