"""V7 curated-style connections: Bezier during layout; optional final orthogonal pass.

Default path: every algorithm-placed connection becomes a Bezier curve, solved for its real end points.
After source-path stamp, ``apply(..., mode='orthogonal')`` writes axis-parallel polylines with real
multimarker bends so official ``edge_orthogonality`` (which ignores Bezier handles) rises.

`edge_routing` (V5) rerouted only the blocked edges and rounded their corners.  The default Bezier path
styles every algorithm-placed connection with `b1` / `b2` (no libavoid).  The shape is not copied from a curated curve: the
curated drawings (`reaction_templates`, 22 Escher maps) only supply CONSTRAINTS for a participant of that role -- the tangent angle at the reaction centre
and at the glyph relative to the chord, the handle lengths as shares of the chord -- and the control points are SOLVED for the end points this map has.
Core participants are nearly straight in the curated maps (median tangent angle 0, small handles), cofactors bow to their own side by ~35-45 degrees; that is
what comes out here, so a short unobstructed connection is a nearly straight Bezier (control points on the chord) and is not bent just to look curved.

Per connection (curve level: everything is checked on the FINAL Bezier, sampled):
  candidates   template constraints, straight-ish (controls at 1/3 and 2/3), the template mirrored, a stronger template bow, small bows to both sides
  choice       staged with a small tolerance: no obstacle (glyph / label of another node or reaction) -> fewest crossings with the connections already
               placed and the expert connections -> no hook / extreme curvature -> shortest detour
  corridor     when every candidate hits an obstacle: A* over the padded corners of the obstacles near the chord (edge_routing._astar) finds a corridor of waypoints
               and a smooth multi-segment Bezier is fitted through it (tangent at each waypoint from its neighbours); a fitted curve that still hits something
               gets shorter handles, then a wider corridor; if nothing works the best candidate stays a Bezier and is listed with `degraded` -- never a polyline
Nothing moves a node; expert connections (`preserved_edges`, `source_paths`, any segment that already has controls) are untouched.  Written into `map.json`
(multimarker nodes for corridor waypoints, `b1` / `b2` on every segment) and `metadata.routed_edges` (kind `bezier_straightish` / `bezier_curved`, template
signature and source maps, `degraded`), so `line_render` draws exactly these segments.  Deterministic.
"""
import math
from collections import defaultdict

import numpy as np

import edge_routing as er
import reaction_templates as rt
from core import CURRENCY
import re

SAMPLES = 24
HANDLE = 1.0 / 3.0
MAX_DETOUR = 1.6
JOINT_MAX_EDGES = 8
FAN_MIN_EDGES = 12      # reactions with more connections than this get side-bundled fan-out candidates and an incremental joint search
FAN_AXIS = 0.6          # share of the reaction axis in the start direction of a fan-out curve
JOINT_MAX_BIG = 90


def _is_cur(bigg):
    return re.sub(r'_[a-z][0-9]?$', '', bigg) in CURRENCY


def _bezier(p, c1, c2, q, n=SAMPLES):
    t = np.linspace(0, 1, n + 1)[:, None]
    return (1 - t) ** 3 * p + 3 * (1 - t) ** 2 * t * c1 + 3 * (1 - t) * t ** 2 * c2 + t ** 3 * q


def _rot(v, deg):
    a = math.radians(deg)
    return np.array([v[0] * math.cos(a) - v[1] * math.sin(a), v[0] * math.sin(a) + v[1] * math.cos(a)])


def _controls(p, q, start_deg, end_deg, h1, h2, side):
    """Control points from the chord and the constraints (angles are signed towards the participant's own side `side`)."""
    chord = q - p
    L = float(np.hypot(*chord))
    d = chord / max(L, 1e-9)
    t0 = _rot(d, start_deg * side)             # initial tangent
    t1 = _rot(d, end_deg * side)               # travel direction at the glyph
    return p + t0 * (h1 * L), q - t1 * (h2 * L)


def _turning(samples):
    v = np.diff(samples, axis=0)
    n = np.hypot(v[:, 0], v[:, 1])
    ok = n > 1e-9
    v = v[ok] / n[ok][:, None]
    if len(v) < 2:
        return 0.0, 0.0
    ang = np.degrees(np.arccos(np.clip((v[:-1] * v[1:]).sum(1), -1, 1)))
    cross = v[:-1, 0] * v[1:, 1] - v[:-1, 1] * v[1:, 0]
    signed = np.where(cross >= 0, ang, -ang)
    return float(np.max(ang)), float(np.abs(signed.sum()))


def _count_hits(samples, others):
    """Sample pieces of a curve that touch an obstacle."""
    if not len(others) or len(samples) < 2:
        return 0
    samples = np.asarray(samples, float)
    A, B = samples[:-1], samples[1:]
    hit = er._hits_many(A, B, others, 2.0)
    # an obstacle sitting on an end of the connection is not counted, but along the curve every touch is
    ends = [0, len(A) - 1]
    hit[ends] = er._blocked_many(A[ends], B[ends], others, 2.0)
    return int(np.count_nonzero(hit.any(axis=1)))


def _evaluate(samples, rects, mine, S_other, label_rects=None, label_mine=None):
    """(glyph hits, crossings, anomaly, detour, label hits) of a sampled curve.  Glyph hits are hard; label hits are soft (the labels are moved afterwards)."""
    hits = _count_hits(samples, rects[~mine])
    label_hits = _count_hits(samples, label_rects[~label_mine]) if label_rects is not None and len(label_rects) else 0
    cross = er._seg_cross(samples[:-1], samples[1:], S_other)
    worst, total = _turning(samples)
    chord = float(np.hypot(*(samples[-1] - samples[0])))
    length = float(np.sum(np.hypot(*np.diff(samples, axis=0).T)))
    anomaly = int(worst > 70.0 or total > 200.0)
    return hits, cross, anomaly, length / max(chord, 1e-9), label_hits


def _fit_through(waypoints, t_start, t_end, handle=0.35):
    """Smooth multi-segment Bezier through the waypoints: [(p, c1, c2, q), ...]."""
    W = [np.array(w, float) for w in waypoints]
    T = []
    for i, w in enumerate(W):
        if i == 0:
            T.append(t_start)
        elif i == len(W) - 1:
            T.append(t_end)
        else:
            d = W[i + 1] - W[i - 1]
            T.append(d / max(float(np.hypot(*d)), 1e-9))
    segs = []
    for i in range(len(W) - 1):
        L = float(np.hypot(*(W[i + 1] - W[i])))
        segs.append((W[i], W[i] + T[i] * handle * L, W[i + 1] - T[i + 1] * handle * L, W[i + 1]))
    return segs


def _corridor(p, q, rects, mine, step_cost=35.0):
    margin = max(140.0, 0.5 * float(np.hypot(*(q - p))))
    lo, hi = np.minimum(p, q) - margin, np.maximum(p, q) + margin
    others = rects[~mine]
    inside = (others[:, 2] >= lo[0]) & (others[:, 0] <= hi[0]) & (others[:, 3] >= lo[1]) & (others[:, 1] <= hi[1])
    window = others[inside]
    out = []
    for bend_cost, pad in er.CANDIDATES:
        grown = np.stack([window[:, 0] - pad, window[:, 1] - pad, window[:, 2] + pad, window[:, 3] + pad], axis=1)
        corners = np.concatenate([grown[:, [0, 1]], grown[:, [2, 1]], grown[:, [0, 3]], grown[:, [2, 3]]]) if len(window) else np.zeros((0, 2))
        keep = np.ones(len(corners), bool)
        for i, c in enumerate(corners):
            if ((window[:, 0] - er.CLEAR < c[0]) & (c[0] < window[:, 2] + er.CLEAR) & (window[:, 1] - er.CLEAR < c[1]) & (c[1] < window[:, 3] + er.CLEAR)).any():
                keep[i] = False
        corners = corners[keep]
        if len(corners) > er.MAX_POINTS:
            dq = q - p
            d = np.abs(dq[0] * (corners[:, 1] - p[1]) - dq[1] * (corners[:, 0] - p[0])) / max(float(np.hypot(*dq)), 1e-9)
            corners = corners[np.argsort(d)[:er.MAX_POINTS]]
        path = er._astar(p, q, corners, window, bend_cost, pad)
        if path:
            out.append(path)
    return out


def _curated_block_rects(data, pad=18.0):
    """Axis-aligned padded boxes of curated (expert) structures for corridor avoidance.

    Built from midmarkers of reactions listed in ``metadata.expert_structures``.
    Empty when metadata is missing (incremental restyle paths).
    """
    meta = (data[0] or {}).get('metadata') or {}
    structures = meta.get('expert_structures') or []
    if not structures:
        return np.zeros((0, 4), float), np.zeros(0, dtype=object)
    nodes = data[1]['nodes']
    reactions = data[1]['reactions']
    mid_of = {}
    for rx in reactions.values():
        for s in rx['segments'].values():
            for nid in (s['from_node_id'], s['to_node_id']):
                n = nodes.get(nid) or {}
                if n.get('node_type') == 'midmarker':
                    mid_of[rx['bigg_id']] = (float(n['x']), float(n['y']))
                    break
            if rx['bigg_id'] in mid_of:
                break
    rects, owners = [], []
    for st in structures:
        pts = [mid_of[r] for r in (st.get('reactions') or []) if r in mid_of]
        if len(pts) < 2:
            continue
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        # Expand by a fixed pad so corridors clear the block body, not just midmarkers.
        span = max(max(xs) - min(xs), max(ys) - min(ys), 40.0)
        grow = max(pad, 0.15 * span)
        rects.append([min(xs) - grow, min(ys) - grow, max(xs) + grow, max(ys) + grow])
        owners.append('block:' + str(st.get('id') or st.get('code') or len(owners)))
    if not rects:
        return np.zeros((0, 4), float), np.zeros(0, dtype=object)
    return np.asarray(rects, float), np.asarray(owners, dtype=object)


def apply(data, skip_edges=(), only=None, mode=None, gem_path=None, pathway=None):
    """Style algorithm-placed connections in place.

    Default / ``mode='bezier'``: Bezier segments from curated-map constraints
    (pipeline stages, hotspot-stable).  ``mode='orthogonal'``: final pass that
    writes axis-parallel polylines with real multimarker bends for official
    ``edge_orthogonality`` (see ``orthogonal_route``).

    When `only` is a set of reaction ids, only those reactions are restyled
    (incremental path used by relocation/peripheral).  Curated block AABBs are
    appended as obstacles so algorithm corridors route around packed modules.
    ``gem_path`` (orthogonal mode) filters the in-venv accept gate to GEM ids.
    Pass a preloaded ``pathway`` to skip re-parsing the SBML.
    """
    if mode == 'orthogonal':
        import orthogonal_route as ortho

        def _stub(rx, mid_id, key, seg, t, P, q, u):
            chord = q - P
            L = float(np.hypot(*chord))
            if L < 1e-6:
                return
            side = 1.0 if (u[0] * chord[1] - u[1] * chord[0]) >= 0 else -1.0
            c1, c2 = _controls(P, q, 18.0, 18.0, HANDLE, HANDLE, side)
            seg['b1'] = dict(x=round(float(c1[0]), 2), y=round(float(c1[1]), 2))
            seg['b2'] = dict(x=round(float(c2[0]), 2), y=round(float(c2[1]), 2))

        return ortho.apply_orthogonal(
            data, skip_edges=skip_edges, only=only, bezier_stub_fn=_stub,
            gem_path=gem_path, pathway=pathway)
    import line_render
    body = data[1]
    nodes, reactions = body['nodes'], body['reactions']
    import layout_quality as lq
    all_rects, all_owner, all_kind, _ = lq._final_obstacles(data)
    is_node = np.array([k == 'node' for k in all_kind])
    rects, owner_arr = all_rects[is_node], np.array(all_owner, dtype=object)[is_node]
    label_rects, label_owner = all_rects[~is_node], np.array(all_owner, dtype=object)[~is_node]
    block_rects, block_owners = _curated_block_rects(data)
    if len(block_rects):
        rects = np.vstack([rects, block_rects]) if len(rects) else block_rects
        owner_arr = np.concatenate([owner_arr, block_owners]) if len(owner_arr) else block_owners
    skip = {tuple(e) for e in skip_edges}
    only_set = set(only) if only is not None else None
    # all current segments, by reaction (the crossing term sees what is already drawn)
    seg_rows, seg_owner = [], []
    for rx in reactions.values():
        for path in line_render._terminals(rx, nodes):
            if nodes[path[-1]].get('node_type') != 'metabolite':
                continue
            pts = np.array(er._path_points(rx, path, nodes))
            for a, b in zip(pts[:-1], pts[1:]):
                seg_rows.append((a[0], a[1], b[0], b[1]))
                seg_owner.append(rx['bigg_id'])
    S_all = np.array(seg_rows, float).reshape(-1, 4)
    S_owner = np.array(seg_owner, dtype=object)
    # spatial indexes: a connection is only checked against what lies near it (iMM904 has thousands of connections)
    # Larger cells on huge canvases keep neighbour lists short
    canvas = body.get('canvas') or {}
    extent = max(float(canvas.get('width', 0) or 0), float(canvas.get('height', 0) or 0), 1.0)
    CELL = 150.0 if extent < 20000 else (250.0 if extent < 80000 else 400.0)

    def cells(lo, hi):
        for gx in range(int(lo[0] // CELL), int(hi[0] // CELL) + 1):
            for gy in range(int(lo[1] // CELL), int(hi[1] // CELL) + 1):
                yield gx, gy
    node_ids = {}
    node_owner_int = np.array([node_ids.setdefault(o, len(node_ids)) for o in owner_arr], dtype=int)
    label_owner_int = np.array([node_ids.setdefault(o, len(node_ids)) for o in label_owner], dtype=int)
    rect_grid, label_grid, seg_grid = defaultdict(list), defaultdict(list), defaultdict(list)
    for i, r in enumerate(rects):
        for c in cells(r[:2], r[2:]):
            rect_grid[c].append(i)
    for i, r in enumerate(label_rects):
        for c in cells(r[:2], r[2:]):
            label_grid[c].append(i)
    rid_ids = {}
    cap = len(S_all) + 60000
    seg_arr = np.zeros((cap, 4))
    seg_own = np.full(cap, -1, dtype=int)
    seg_arr[:len(S_all)] = S_all
    seg_own[:len(S_all)] = [rid_ids.setdefault(o, len(rid_ids)) for o in S_owner]
    seg_n = [len(S_all)]
    for k in range(len(S_all)):
        for c in cells(np.minimum(S_all[k, :2], S_all[k, 2:]), np.maximum(S_all[k, :2], S_all[k, 2:])):
            seg_grid[c].append(k)
    cur_edge = {}

    def _gather(grid, lo, hi):
        # Order does not matter for hit tests; skip sorted() on the hot path
        seen = set()
        out = []
        for c in cells(lo, hi):
            for i in grid.get(c, ()):
                if i not in seen:
                    seen.add(i)
                    out.append(i)
        return np.fromiter(out, dtype=int, count=len(out)) if out else np.zeros(0, dtype=int)

    def ev(smp):
        lo, hi = smp.min(axis=0) - 8, smp.max(axis=0) + 8
        ni = _gather(rect_grid, lo, hi)
        li = _gather(label_grid, lo, hi)
        si = _gather(seg_grid, lo, hi)
        own = cur_edge['own']
        r_sub, l_sub = (rects[ni] if len(ni) else rects[:0]), (label_rects[li] if len(li) else label_rects[:0])
        s_sub = seg_arr[si][seg_own[si] != cur_edge['rid']] if len(si) else np.zeros((0, 4))
        return _evaluate(smp, r_sub, np.isin(node_owner_int[ni], own) if len(ni) else np.zeros(0, bool), s_sub, l_sub,
                         np.isin(label_owner_int[li], own) if len(li) else np.zeros(0, bool))
    holders = defaultdict(set)
    for rx in reactions.values():
        for s in rx['segments'].values():
            holders[s['to_node_id']].add(rx['bigg_id'])
    rect_of = {o: i for i, o in enumerate(owner_arr)}
    lab_of = defaultdict(list)
    for i, o in enumerate(label_owner):
        lab_of[o].append(i)
    joint_store = {}
    own_curves = defaultdict(list)      # sampled curves already decided for each reaction: a later curve of the same reaction must not cross them away from the reaction point
    report = dict(relocated=0, edges=0, expert_or_curved=0, bezier=0, straightish=0, curved=0, corridor=0, degraded=0, no_template=0, hits_after=0, crossings_after=0)
    routed, counter = [], [0]

    def new_node(x, y):
        counter[0] += 1
        nid = f'es{counter[0]}'
        nodes[nid] = dict(node_type='multimarker', x=round(float(x), 2), y=round(float(y), 2))
        return nid
    for rx in sorted(reactions.values(), key=lambda r: r['bigg_id']):
        rid = rx['bigg_id']
        if only is not None and rid not in only:
            continue
        coeff = {}
        for m in rx.get('metabolites', []):
            coeff.setdefault(m['bigg_id'], m['coefficient'])
        mid_id = next((s['from_node_id'] for s in rx['segments'].values() if nodes[s['from_node_id']].get('node_type') == 'midmarker'), None)
        if mid_id is None:
            continue
        P = np.array([nodes[mid_id]['x'], nodes[mid_id]['y']])
        edges = []
        for key, seg in list(rx['segments'].items()):
            if seg['from_node_id'] != mid_id or nodes[seg['to_node_id']].get('node_type') != 'metabolite':
                continue
            t = nodes[seg['to_node_id']]
            edges.append((key, seg, t))
        if not edges:
            continue
        subs = [np.array([t['x'], t['y']]) - P for _, _, t in edges if coeff.get(t['bigg_id'], 0) < 0]
        prods = [np.array([t['x'], t['y']]) - P for _, _, t in edges if coeff.get(t['bigg_id'], 0) > 0]
        u = (np.mean(prods, axis=0) if prods else np.zeros(2)) - (np.mean(subs, axis=0) if subs else np.zeros(2))
        u = u / np.hypot(*u) if np.hypot(*u) > 1e-6 else np.array([1.0, 0.0])
        cs = sum(1 for _, _, t in edges if coeff.get(t['bigg_id'], 0) < 0 and not _is_cur(t['bigg_id']))
        cp = sum(1 for _, _, t in edges if coeff.get(t['bigg_id'], 0) > 0 and not _is_cur(t['bigg_id']))
        ks = sum(1 for _, _, t in edges if coeff.get(t['bigg_id'], 0) < 0 and _is_cur(t['bigg_id']))
        kp = sum(1 for _, _, t in edges if coeff.get(t['bigg_id'], 0) > 0 and _is_cur(t['bigg_id']))
        chems = defaultdict(set)
        for _, _, t in edges:
            chems[re.sub(r'_[a-z][0-9]?$', '', t['bigg_id'])].add(t['bigg_id'])
        transport = any(len(v) > 1 for v in chems.values())
        for key, seg, t in edges:
            report['edges'] += 1
            report_key = (rid, t['bigg_id'])
            if report_key in skip or seg.get('b1') or seg.get('b2') or str(key).startswith('spseg'):
                report['expert_or_curved'] += 1
                continue
            q = np.array([t['x'], t['y']])
            chord = q - P
            L = float(np.hypot(*chord))
            if L < 1e-6:
                continue
            role = 'substrate' if coeff.get(t['bigg_id'], 0) < 0 else 'product'
            cur = _is_cur(t['bigg_id'])
            side = 1.0 if (u[0] * chord[1] - u[1] * chord[0]) >= 0 else -1.0
            st = rt.stats((cs, cp, ks, kp), transport, role, cur, True)
            if st is None:
                report['no_template'] += 1
            mine = (owner_arr == mid_id) | (owner_arr == seg['to_node_id'])
            label_mine = (label_owner == mid_id) | (label_owner == seg['to_node_id'])
            cur_edge['own'] = [node_ids.setdefault(mid_id, len(node_ids)), node_ids.setdefault(seg['to_node_id'], len(node_ids))]
            cur_edge['rid'] = rid_ids.setdefault(rid, len(rid_ids))
            base = st or dict(start_rel=0.0, end_rel=0.0, h1=HANDLE, h2=HANDLE, dev=0.0)
            s0, e0 = float(np.clip(base['start_rel'], -50, 50)), float(np.clip(base['end_rel'], -50, 50))
            h1, h2 = float(np.clip(base['h1'], 0.2, 0.5)), float(np.clip(base['h2'], 0.2, 0.5))
            cands = [('template', s0, e0, h1, h2, 1.0), ('straightish', 0.0, 0.0, HANDLE, HANDLE, 1.0), ('template_mirror', s0, e0, h1, h2, -1.0),
                     ('template_strong', s0 * 1.5, e0 * 1.5, min(h1 * 1.3, 0.55), min(h2 * 1.3, 0.55), 1.0),
                     ('bow_a', 14.0, 14.0, HANDLE, HANDLE, 1.0), ('bow_b', 14.0, 14.0, HANDLE, HANDLE, -1.0)]
            if abs(s0) < 3 and abs(e0) < 3:
                cands = [cands[1]] + [c for c in cands if c[0] != 'straightish']  # a core participant that is straight in the expert maps stays nearly straight
            if len(edges) > FAN_MIN_EDGES:
                # a very large reaction (biomass): bundle by side.  Substrates leave towards -u, products towards +u (the reaction axis), then bend to their glyph;
                # curves of one side leave in nearly the same direction, so they fan out instead of radiating from the point in every direction
                d = chord / L
                want = (u if role == 'product' else -u) * FAN_AXIS + d * (1.0 - FAN_AXIS)
                ang = math.degrees(math.atan2(d[0] * want[1] - d[1] * want[0], float(d @ want)))
                fan = float(np.clip(ang * side, -75.0, 75.0))
                fans = [('fan', fan, 0.0, 0.40, 0.30, 1.0), ('fan_soft', fan * 0.5, 0.0, 0.36, 0.33, 1.0)]
                cands = fans + cands
            best = None
            for name, a0, a1, g1, g2, mirror in cands:
                c1, c2 = _controls(P, q, a0, a1, g1, g2, side * mirror)
                smp = _bezier(P, c1, c2, q)
                hits, cross, anomaly, detour, lhits = ev(smp)
                cross += sum(_pair_cross(smp[3:], o[3:]) for o in own_curves[rid])
                score = (hits > 0, cross, anomaly, detour > MAX_DETOUR, lhits, round(detour, 2) if detour > 1.15 else 0)
                if best is None or score < best[0]:
                    best = (score, name, c1, c2, smp, hits, cross)
                if score[:3] == (False, 0, 0) and lhits == 0:
                    break
            joint_store.setdefault(rid, []).append(dict(key=key, seg=seg, P=P.copy(), q=q.copy(), side=side, cands=cands, mid=mid_id, to=seg['to_node_id'], chosen=best[1],
                                                       own=list(cur_edge['own']), rid_i=cur_edge['rid'], tname=t['bigg_id']))
            chain = None
            degraded = False
            # Corridor search is expensive; skip it on incremental restyles when a light hit is acceptable, or when the best candidate is already clean
            do_corridor = best[5] > 0 and (only_set is None or best[5] >= 2)
            if do_corridor:   # every candidate touches an obstacle: find a corridor and fit a smooth multi-segment Bezier through it
                for path in _corridor(P, q, rects, mine):
                    t_start = _rot(chord / L, s0 * side)
                    t_end = _rot(chord / L, e0 * side)
                    for handle in (0.35, 0.22):
                        segs = _fit_through(path, t_start, t_end, handle)
                        smp = np.vstack([_bezier(a, b, c, d, 12) for a, b, c, d in segs])
                        hits, cross, anomaly, detour, _l = ev(smp)
                        own_x = sum(_pair_cross(smp[3:], o[3:]) for o in own_curves[rid])
                        if hits == 0 and not anomaly and own_x == 0:
                            chain = (segs, smp, hits, cross)
                            break
                    if chain:
                        break
                if chain is None:
                    # a glyph only this reaction uses (a cofactor drawn beside it) may move a little: it belongs to the algorithm's layout, not to an expert structure
                    nid = seg['to_node_id']
                    if len(holders[nid]) == 1 and cur:
                        base = math.atan2(q[1] - P[1], q[0] - P[0])
                        options = []
                        for radius in (55.0, 75.0, 95.0):
                            for k in range(12):
                                ang = base + ((k + 1) // 2) * (math.pi / 6) * (1 if k % 2 else -1)
                                q2 = P + radius * np.array([math.cos(ang), math.sin(ang)])
                                box = np.array([q2[0] - 7, q2[1] - 7, q2[0] + 7, q2[1] + 7])
                                others_n = rects[owner_arr != nid]
                                if len(others_n) and ((others_n[:, 0] < box[2] + 4) & (box[0] - 4 < others_n[:, 2]) & (others_n[:, 1] < box[3] + 4) & (box[1] - 4 < others_n[:, 3])).any():
                                    continue
                                c1b, c2b = _controls(P, q2, 0.0, 0.0, HANDLE, HANDLE, side)
                                smp2 = _bezier(P, c1b, c2b, q2)
                                h2, x2, a2, d2, l2 = ev(smp2)
                                if h2 == 0 and not a2:
                                    options.append(((x2, l2, radius, k), q2, c1b, c2b, smp2))
                        if options:
                            _, q2, c1b, c2b, smp2 = min(options, key=lambda o: o[0])
                            dq = q2 - q
                            t['x'], t['y'] = round(float(q2[0]), 2), round(float(q2[1]), 2)
                            if 'label_x' in t:
                                t['label_x'], t['label_y'] = round(float(t['label_x'] + dq[0]), 2), round(float(t['label_y'] + dq[1]), 2)
                            rects[rect_of[nid]] = [q2[0] - er.NODE_HALF, q2[1] - er.NODE_HALF, q2[0] + er.NODE_HALF, q2[1] + er.NODE_HALF]
                            for c in cells(rects[rect_of[nid]][:2], rects[rect_of[nid]][2:]):
                                rect_grid[c].append(rect_of[nid])
                            for li in lab_of.get(nid, []):
                                label_rects[li] = label_rects[li] + np.array([dq[0], dq[1], dq[0], dq[1]])
                                for c in cells(label_rects[li][:2], label_rects[li][2:]):
                                    label_grid[c].append(li)
                            q = q2
                            chord = q - P
                            L = float(np.hypot(*chord))
                            best = ((False, 0, 0, False, 0, 0), 'relocated_glyph', c1b, c2b, smp2, 0, 0)
                            report['relocated'] += 1
                            joint_store[rid][-1]['fixed'] = True
                    if chain is None and best[5] > 0:
                        degraded = True
            if chain:
                joint_store[rid][-1]['fixed'] = True
                segs, smp, hits, cross = chain
                report['corridor'] += 1
                del rx['segments'][key]
                prev = mid_id
                for i, (a, b, c, d) in enumerate(segs):
                    last = i == len(segs) - 1
                    to = seg['to_node_id'] if last else new_node(*d)
                    counter[0] += 1
                    rx['segments'][f'es_s{counter[0]}'] = dict(from_node_id=prev, to_node_id=to, b1=dict(x=round(float(b[0]), 2), y=round(float(b[1]), 2)),
                                                               b2=dict(x=round(float(c[0]), 2), y=round(float(c[1]), 2)))
                    prev = to
                kind, name = 'bezier_curved', 'corridor'
            else:
                _, name, c1, c2, smp, hits, cross = best
                seg['b1'], seg['b2'] = dict(x=round(float(c1[0]), 2), y=round(float(c1[1]), 2)), dict(x=round(float(c2[0]), 2), y=round(float(c2[1]), 2))
                dev = float(np.max(np.abs((smp[:, 0] - P[0]) * chord[1] - (smp[:, 1] - P[1]) * chord[0]) / L)) / L
                kind = 'bezier_straightish' if dev < 0.05 else 'bezier_curved'
            own_curves[rid].append(smp)
            report['bezier'] += 1
            report['straightish' if kind == 'bezier_straightish' else 'curved'] += 1
            if degraded:
                report['degraded'] += 1
            report['hits_after'] += int(hits > 0)
            report['crossings_after'] += int(cross)
            routed.append(dict(reaction=rid, metabolite=t['bigg_id'], kind=kind, candidate=name, template_signature='s%d,%d,%d,%d' % (cs, cp, ks, kp) if st else None,
                               template_sources=[list(x) for x in st['source']] if st else None, degraded=bool(degraded)))
            for a, b in zip(smp[:-1], smp[1:]):
                k = seg_n[0]
                if k >= len(seg_arr):
                    seg_arr = np.vstack([seg_arr, np.zeros((len(seg_arr), 4))])
                    seg_own = np.concatenate([seg_own, np.full(len(seg_own), -1, dtype=int)])
                seg_arr[k] = (a[0], a[1], b[0], b[1])
                seg_own[k] = cur_edge['rid']
                seg_n[0] += 1
                for c in cells(np.minimum(a, b), np.maximum(a, b)):
                    seg_grid[c].append(k)
    # Incremental restyles only joint-route the reactions that were redrawn
    report.update(joint_route(data, joint_store, ev, cur_edge, rects, owner_arr, label_rects, label_owner, skip, only=only_set))
    data[0].setdefault('metadata', {})['routed_edges'] = routed
    return report


def _pair_cross(sa, sb):
    """Number of crossings between two sampled polylines (vectorised all pairs of pieces)."""
    if len(sa) < 2 or len(sb) < 2:
        return 0
    p, r = sa[:-1][:, None, :], (sa[1:] - sa[:-1])[:, None, :]
    q, s = sb[:-1][None, :, :], (sb[1:] - sb[:-1])[None, :, :]
    den = r[..., 0] * s[..., 1] - r[..., 1] * s[..., 0]
    with np.errstate(divide='ignore', invalid='ignore'):
        tt = ((q[..., 0] - p[..., 0]) * s[..., 1] - (q[..., 1] - p[..., 1]) * s[..., 0]) / den
        uu = ((q[..., 0] - p[..., 0]) * r[..., 1] - (q[..., 1] - p[..., 1]) * r[..., 0]) / den
    return int(np.sum((np.abs(den) > 1e-9) & (tt > 1e-3) & (tt < 1 - 1e-3) & (uu > 1e-3) & (uu < 1 - 1e-3)))


def _own_crossings(curves, skip_head=3, movable=None):
    """Crossings between the curves of ONE reaction away from the shared reaction point (pairs of two fixed curves are not counted: nothing can change them)."""
    n = 0
    for a in range(len(curves)):
        for b in range(a + 1, len(curves)):
            if movable is not None and not (movable[a] or movable[b]):
                continue
            n += _pair_cross(curves[a][skip_head:], curves[b][skip_head:])
    return n


def joint_route(data, store, ev, cur_edge, rects, owner_arr, label_rects, label_owner, skip, only=None):
    """Reaction-level joint choice: the curves of one reaction are not independent.  For every reaction whose curves cross each other away from the reaction point,
    the candidate of each curve is re-chosen jointly (coordinate descent over the stored candidates, a few rounds), scored as
    (curves through glyphs, crossings among the reaction's own curves, crossings with other connections, anomalies, labels hit, detour).  Only a strictly better
    combination replaces the current one.  A crossing that no candidate combination can remove stays and is counted.

    When `only` is a set of reaction ids (incremental restyle), other reactions are left untouched.
    """
    import line_render
    rep = dict(joint_reactions_with_own_crossings=0, joint_fixed=0, joint_own_crossings_before=0, joint_own_crossings_after=0)
    body = data[1]
    by_id = {r['bigg_id']: r for r in body['reactions'].values()}
    for rid, items in store.items():
        if only is not None and rid not in only:
            continue
        movable = [it for it in items if not it.get('fixed') and not str(it['key']).startswith('spseg')]
        if not movable or len(movable) > JOINT_MAX_BIG:
            continue
        big = len(movable) > JOINT_MAX_EDGES     # very large reactions (biomass): fan-out candidates + incremental search below
        # curves of this reaction that cannot be changed here (expert curves, corridor chains): they are obstacles for the movable ones
        mov_to = {it['to'] for it in movable}
        fixed_curves = []
        for path in line_render._terminals(by_id[rid], body['nodes']):
            if path[-1] not in mov_to and body['nodes'][path[-1]].get('node_type') == 'metabolite':
                fixed_curves.append(np.array(er._path_points(by_id[rid], path, body['nodes']), float))
        if len(movable) + len(fixed_curves) < 2:
            continue
        curves = []
        for it in movable:
            seg = it['seg']
            c1 = np.array([seg['b1']['x'], seg['b1']['y']]) if seg.get('b1') else None
            c2 = np.array([seg['b2']['x'], seg['b2']['y']]) if seg.get('b2') else None
            if c1 is None:
                break
            curves.append(_bezier(it['P'], c1, c2, it['q']))
        if len(curves) != len(movable):
            continue
        mflags = [True] * len(curves) + [False] * len(fixed_curves)
        now = _own_crossings(curves + fixed_curves, movable=mflags)
        if now == 0:
            continue
        rep['joint_reactions_with_own_crossings'] += 1
        rep['joint_own_crossings_before'] += now
        table = []
        for it in movable:
            cur_edge['own'] = it['own']
            cur_edge['rid'] = it['rid_i']
            row = []
            for name, a0, a1, g1, g2, mirror in it['cands']:
                c1, c2 = _controls(it['P'], it['q'], a0, a1, g1, g2, it['side'] * mirror)
                smp = _bezier(it['P'], c1, c2, it['q'])
                hits, cross, anomaly, detour, lhits = ev(smp)
                row.append((name, c1, c2, smp, (hits > 0, cross, anomaly, detour > MAX_DETOUR, lhits), detour))
            table.append(row)

        pair_cache = {}

        def own_cost(choice):
            smps = [table[i][c][3] for i, c in enumerate(choice)]
            n = 0
            for i in range(len(smps)):
                for j in range(i + 1, len(smps)):
                    key = (i, choice[i], j, choice[j])
                    if key not in pair_cache:
                        pair_cache[key] = _pair_cross(smps[i][3:], smps[j][3:])
                    n += pair_cache[key]
                for k, fc in enumerate(fixed_curves):
                    key = (i, choice[i], 'f', k)
                    if key not in pair_cache:
                        pair_cache[key] = _pair_cross(smps[i][3:], fc[3:])
                    n += pair_cache[key]
            return n

        def cost(choice):
            ind = [table[i][c][4] for i, c in enumerate(choice)]
            return (sum(x[0] for x in ind), own_cost(choice), sum(x[1] for x in ind), sum(x[2] for x in ind), sum(x[4] for x in ind), round(sum(max(table[i][c][5] - 1.15, 0) for i, c in enumerate(choice)), 2))
        choice = [int(np.argmin([r[4] + (r[5],) for r in row])) for row in table]
        if big:
            # incremental coordinate descent: changing curve i only changes its own pairs, so the total is updated by the difference
            def contrib(i, c, ch):
                n = 0
                for j in range(len(table)):
                    if j == i:
                        continue
                    key = (min(i, j), c if i < j else ch[j], max(i, j), ch[j] if i < j else c)
                    if key not in pair_cache:
                        si, sj = table[i][c][3], table[j][ch[j]][3]
                        pair_cache[key] = _pair_cross((si if i < j else sj)[3:], (sj if i < j else si)[3:])
                    n += pair_cache[key]
                for k, fc in enumerate(fixed_curves):
                    key = (i, c, 'f', k)
                    if key not in pair_cache:
                        pair_cache[key] = _pair_cross(table[i][c][3][3:], fc[3:])
                    n += pair_cache[key]
                return n
            ind0 = lambda i, c: table[i][c][4]
            for _ in range(3):
                improved = False
                for i in range(len(table)):
                    cur_c = choice[i]
                    base_k = (ind0(i, cur_c)[0], contrib(i, cur_c, choice), ind0(i, cur_c)[1], ind0(i, cur_c)[2], ind0(i, cur_c)[4], round(max(table[i][cur_c][5] - 1.15, 0), 2))
                    for c in range(len(table[i])):
                        if c == cur_c:
                            continue
                        k = (ind0(i, c)[0], contrib(i, c, choice), ind0(i, c)[1], ind0(i, c)[2], ind0(i, c)[4], round(max(table[i][c][5] - 1.15, 0), 2))
                        if k < base_k:
                            choice[i], base_k, improved = c, k, True
                if not improved:
                    break
            best = (0, own_cost(choice))
        else:
            best = cost(choice)
        for _ in range(0 if big else 3):
            improved = False
            for i in range(len(table)):
                for c in range(len(table[i])):
                    if c == choice[i]:
                        continue
                    trial = choice[:i] + [c] + choice[i + 1:]
                    k = cost(trial)
                    if k < best:
                        choice, best, improved = trial, k, True
            if not improved:
                break
        rep['joint_own_crossings_after'] += best[1]
        if best[1] < now:
            rep['joint_fixed'] += 1
        for it, row, c in zip(movable, table, choice):
            name, c1, c2 = row[c][0], row[c][1], row[c][2]
            it['seg']['b1'], it['seg']['b2'] = dict(x=round(float(c1[0]), 2), y=round(float(c1[1]), 2)), dict(x=round(float(c2[0]), 2), y=round(float(c2[1]), 2))
    return rep
