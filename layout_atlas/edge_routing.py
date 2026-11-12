"""V5 hybrid edge routing: connections of the reactions the algorithm placed are routed only where they need to be, then smoothed.

Until now every connection of a small map went through libavoid's purely orthogonal router (`routing.py`, hashed): all edges, expert or not,
became staircases of right angles, long detours crossed empty space, and the algorithm's lines looked nothing like the expert curves next to them.
Here the map is first exported without routing (every connection a straight line from the reaction to its metabolite glyph) and this module
then works on the finished `map.json`, as an Escher map:

  expert connections   untouched: the structure's own lines and curves (`preserved_edges`, `source_paths`, any segment that already has
                       Bezier controls)
  straight, free       an algorithm-placed connection whose straight line crosses no other glyph or label stays a straight line
  blocked              a straight line that crosses another node, reaction point or label gets a path around it: A* over the corners of the
                       obstacles near the line (corners padded), bends cost, 2-3 candidates with different bend cost and padding, the one with
                       no obstacle hits, then fewest crossings with other connections, then fewest bends, then shortest wins
  smoothed             the right-angle-free corners of a detour are rounded (cubic Bezier through the corner, radius limited by the adjacent
                       segments); a rounded corner that would touch an obstacle stays a sharp corner

The result is written into `map.json` itself (multimarker nodes for the corners, `b1`/`b2` for the rounded corners) and listed in
`metadata.routed_edges`; `line_render` draws exactly these segments, so the SVG and the Escher JSON show the same line.
`measure(map_json)` counts, for the same set of edges in any map, bends, obstacle hits and crossings (used for the libavoid baseline too).
Nothing here moves a node; GEM content is untouched.
"""
import math
from collections import defaultdict

import numpy as np

NODE_HALF = 7.0      # half width of a metabolite / reaction glyph box
LABEL_PAD = 2.0
CORNER_PAD = 10.0    # distance of the detour corners from an obstacle box
CLEAR = 6.0          # how far a free straight line has to stay from an obstacle box
MAX_POINTS = 260
CANDIDATES = ((35.0, CORNER_PAD), (120.0, CORNER_PAD), (35.0, CORNER_PAD * 2.2))  # (bend cost per 90 degrees, corner padding)
END_PAD = 14.0       # an obstacle this close to an end of a connection is not counted against it
ROUND = 20.0         # largest corner radius of a smoothed corner


def _font():
    import drawing
    return drawing.FONT


def obstacles(data):
    """All boxes that a connection must not cross: glyphs and labels.  Returns (rects [n,4], owner [n] -> node id or reaction key, kind)."""
    body = data[1]
    font = _font()
    rects, owner = [], []
    for nid, n in body['nodes'].items():
        t = n.get('node_type')
        if t == 'metabolite':
            rects.append((n['x'] - NODE_HALF, n['y'] - NODE_HALF, n['x'] + NODE_HALF, n['y'] + NODE_HALF))
            owner.append(nid)
            if 'label_x' in n:
                w = font.getlength(n.get('bigg_id', ''))
                rects.append((n['label_x'] - LABEL_PAD, n['label_y'] - 14, n['label_x'] + w + LABEL_PAD, n['label_y'] + 4))
                owner.append(nid)
        elif t == 'midmarker':
            rects.append((n['x'] - NODE_HALF + 1, n['y'] - NODE_HALF + 1, n['x'] + NODE_HALF - 1, n['y'] + NODE_HALF - 1))
            owner.append(nid)
    for rx in body['reactions'].values():
        mid = next((s['from_node_id'] for s in rx['segments'].values() if body['nodes'][s['from_node_id']].get('node_type') == 'midmarker'), None)
        if mid is None:
            continue
        w = font.getlength(rx['bigg_id'])
        rects.append((rx['label_x'] - LABEL_PAD, rx['label_y'] - 14, rx['label_x'] + w + LABEL_PAD, rx['label_y'] + 4))
        owner.append(mid)
    return np.array(rects, float).reshape(-1, 4), owner


def _hits_many(P, Q, rects, pad=0.0):
    """Boolean [segment, rect]: does segment P[i]-Q[i] cross the rect (box grown by pad)?  Slab test, vectorised."""
    P, Q = np.asarray(P, float).reshape(-1, 2), np.asarray(Q, float).reshape(-1, 2)
    x0, y0, x1, y1 = rects[:, 0] - pad, rects[:, 1] - pad, rects[:, 2] + pad, rects[:, 3] + pad
    d = Q - P
    shape = (len(P), len(rects))
    tmin, tmax = np.zeros(shape), np.ones(shape)
    ok = np.ones(shape, bool)
    for ax, lo, hi in ((0, x0, x1), (1, y0, y1)):
        flat = (np.abs(d[:, ax]) < 1e-9)[:, None]
        p = P[:, ax, None]
        ok &= ~flat | ((p >= lo) & (p <= hi))
        with np.errstate(divide='ignore', invalid='ignore'):
            a, b = (lo - p) / d[:, ax, None], (hi - p) / d[:, ax, None]
        tmin = np.where(flat, tmin, np.maximum(tmin, np.minimum(a, b)))
        tmax = np.where(flat, tmax, np.minimum(tmax, np.maximum(a, b)))
    return ok & (tmin <= tmax)


def _hits(p, q, rects, pad=0.0):
    """Boolean per rect: does the segment p-q cross it (box grown by pad)?"""
    return _hits_many(p, q, rects, pad)[0]


def _near_many(P, rects, pad):
    """Boolean [point, rect]: is the point inside the box grown by pad (an obstacle right at the end of a connection cannot be avoided)?"""
    x, y = P[:, 0, None], P[:, 1, None]
    return (rects[:, 0] - pad <= x) & (x <= rects[:, 2] + pad) & (rects[:, 1] - pad <= y) & (y <= rects[:, 3] + pad)


def _blocked_many(P, Q, rects, pad):
    """Boolean [segment, rect]: does P[i]-Q[i] cross a box that is not sitting on one of its two ends?"""
    P, Q = np.asarray(P, float).reshape(-1, 2), np.asarray(Q, float).reshape(-1, 2)
    if not len(rects):
        return np.zeros((len(P), 0), bool)
    return _hits_many(P, Q, rects, pad) & ~(_near_many(P, rects, END_PAD) | _near_many(Q, rects, END_PAD))


def _blocked(p, q, rects, pad):
    """Does p-q cross a box that is not sitting on one of its two ends?"""
    return _blocked_many(p, q, rects, pad)[0]


def _seg_cross(a, b, S):
    """Number of proper crossings between the pieces a[i]-b[i] (one piece or an array of them) and the segments in S [k,4] (x1,y1,x2,y2)."""
    if not len(S):
        return 0
    p = np.asarray(a, float).reshape(-1, 2)
    r = np.asarray(b, float).reshape(-1, 2) - p
    q, s = S[:, :2], S[:, 2:] - S[:, :2]
    den = r[:, 0, None] * s[:, 1] - r[:, 1, None] * s[:, 0]
    with np.errstate(divide='ignore', invalid='ignore'):
        t = ((q[:, 0] - p[:, 0, None]) * s[:, 1] - (q[:, 1] - p[:, 1, None]) * s[:, 0]) / den
        u = ((q[:, 0] - p[:, 0, None]) * r[:, 1, None] - (q[:, 1] - p[:, 1, None]) * r[:, 0, None]) / den
    return int(np.sum((np.abs(den) > 1e-9) & (t > 1e-3) & (t < 1 - 1e-3) & (u > 1e-3) & (u < 1 - 1e-3)))


def _astar(start, goal, corners, rects, bend_cost, pad):
    """Shortest bend-aware path start -> goal through `corners` avoiding `rects` (grown by CLEAR).  Returns list of points or None."""
    pts = np.vstack([start, goal, corners]) if len(corners) else np.vstack([start, goal])
    n = len(pts)
    goal_i = 1

    def visible(i):
        out = ~_blocked_many(np.broadcast_to(pts[i], pts.shape), pts, rects, CLEAR).any(axis=1)
        out[i] = False
        return out
    import heapq
    start_state = (0, -1)
    gbest, parent = {start_state: 0.0}, {start_state: None}
    heap = [(float(np.hypot(*(pts[0] - pts[goal_i]))), 0.0, start_state)]
    vis_cache = {}
    while heap:
        f, g, state = heapq.heappop(heap)
        if g > gbest.get(state, float('inf')):
            continue
        i, prev = state
        if i == goal_i:
            path = []
            while state is not None:
                path.append(state[0])
                state = parent[state]
            return [pts[k].tolist() for k in reversed(path)]
        if i not in vis_cache:
            vis_cache[i] = visible(i)
        for j in np.flatnonzero(vis_cache[i]):
            step = float(np.hypot(*(pts[j] - pts[i])))
            turn = 0.0
            if prev >= 0:
                a, b = pts[i] - pts[prev], pts[j] - pts[i]
                na, nb = float(np.hypot(*a)), float(np.hypot(*b))
                if na > 1e-9 and nb > 1e-9:
                    turn = math.degrees(math.acos(max(-1.0, min(1.0, float(a @ b) / (na * nb)))))
            ng = g + step + bend_cost * turn / 90.0
            ns = (int(j), i)
            if ng < gbest.get(ns, float('inf')):
                gbest[ns] = ng
                parent[ns] = state
                heapq.heappush(heap, (ng + float(np.hypot(*(pts[j] - pts[goal_i]))), ng, ns))
    return None


def _round_corners(path, rects, owner_mask):
    """-> list of ('L'|'C', ...) steps along `path` with right-angle-free corners; a corner that touches an obstacle stays sharp.
    Each step is (point, b1, b2): b1/b2 None for a straight step."""
    P = [np.array(p, float) for p in path]
    steps = []
    cur = P[0]
    for k in range(1, len(P) - 1):
        v, a, b = P[k], P[k - 1], P[k + 1]
        la, lb = float(np.hypot(*(v - a))), float(np.hypot(*(b - v)))
        r = min(ROUND, 0.4 * la, 0.4 * lb)
        turn_vec = (v - a) / max(la, 1e-9), (b - v) / max(lb, 1e-9)
        if r < 4.0 or abs(float(turn_vec[0] @ turn_vec[1]) - 1.0) < 1e-3:
            steps.append((v.tolist(), None, None))
            continue
        p1 = v - turn_vec[0] * r
        p2 = v + turn_vec[1] * r
        c1, c2 = p1 + (v - p1) * .55, p2 + (v - p2) * .55
        # sample the rounded corner: it must stay clear of obstacles
        t = np.linspace(0, 1, 9)[:, None]
        curve = (1 - t) ** 3 * p1 + 3 * (1 - t) ** 2 * t * c1 + 3 * (1 - t) * t ** 2 * c2 + t ** 3 * p2
        bad = False
        for u, w in zip(curve[:-1], curve[1:]):
            if _blocked(u, w, rects[owner_mask], 2.0).any():
                bad = True
                break
        if bad:
            steps.append((v.tolist(), None, None))
        else:
            steps.append((p1.tolist(), None, None))
            steps.append((p2.tolist(), c1.tolist(), c2.tolist()))
    steps.append((P[-1].tolist(), None, None))
    return steps


def route(data, skip_edges=()):
    """Route the algorithm-placed connections of `data` ([header, body] Escher map, exported without routing) in place.
    `skip_edges`: {(reaction id, metabolite id)} of expert connections.  Returns a report dict."""
    body = data[1]
    nodes, reactions = body['nodes'], body['reactions']
    rects, owner = obstacles(data)
    owner_arr = np.array(owner, dtype=object)
    skip = {tuple(e) for e in skip_edges}
    # every current straight segment, for the crossing term
    segments = []
    for rx in reactions.values():
        for s in rx['segments'].values():
            a, b = nodes[s['from_node_id']], nodes[s['to_node_id']]
            segments.append((a['x'], a['y'], b['x'], b['y']))
    S_all = np.array(segments, float).reshape(-1, 4)
    report = dict(edges=0, expert_or_curved=0, free_straight=0, rerouted=0, smoothed_corners=0, unresolved=0, bends_added=0)
    counter = [0]
    routed = []

    def new_node(x, y):
        counter[0] += 1
        nid = f'er{counter[0]}'
        nodes[nid] = dict(node_type='multimarker', x=round(float(x), 2), y=round(float(y), 2))
        return nid
    for rx in reactions.values():
        rid = rx['bigg_id']
        for key, seg in list(rx['segments'].items()):
            a, t = nodes[seg['from_node_id']], nodes[seg['to_node_id']]
            if a.get('node_type') != 'midmarker' or t.get('node_type') != 'metabolite':
                continue
            report['edges'] += 1
            if (rid, t.get('bigg_id')) in skip or seg.get('b1') or seg.get('b2') or str(key).startswith('spseg'):
                report['expert_or_curved'] += 1
                continue
            p, q = np.array([a['x'], a['y']]), np.array([t['x'], t['y']])
            mine = (owner_arr == seg['from_node_id']) | (owner_arr == seg['to_node_id'])
            others = rects[~mine]
            if not len(others) or not _blocked(p, q, others, CLEAR).any():
                report['free_straight'] += 1
                continue
            # window of obstacles around the line
            margin = max(140.0, 0.5 * float(np.hypot(*(q - p))))
            lo, hi = np.minimum(p, q) - margin, np.maximum(p, q) + margin
            inside = (others[:, 2] >= lo[0]) & (others[:, 0] <= hi[0]) & (others[:, 3] >= lo[1]) & (others[:, 1] <= hi[1])
            window = others[inside]
            best = None
            for bend_cost, pad in CANDIDATES:
                grown = np.stack([window[:, 0] - pad, window[:, 1] - pad, window[:, 2] + pad, window[:, 3] + pad], axis=1)
                corners = np.concatenate([grown[:, [0, 1]], grown[:, [2, 1]], grown[:, [0, 3]], grown[:, [2, 3]]])
                # corners inside another box are useless
                keep = np.ones(len(corners), bool)
                for i, c in enumerate(corners):
                    if ((window[:, 0] - CLEAR < c[0]) & (c[0] < window[:, 2] + CLEAR) & (window[:, 1] - CLEAR < c[1]) & (c[1] < window[:, 3] + CLEAR)).any():
                        keep[i] = False
                corners = corners[keep]
                if len(corners) > MAX_POINTS:
                    dq = q - p
                    d = np.abs(dq[0] * (corners[:, 1] - p[1]) - dq[1] * (corners[:, 0] - p[0])) / max(float(np.hypot(*dq)), 1e-9)
                    corners = corners[np.argsort(d)[:MAX_POINTS]]
                path = _astar(p, q, corners, window, bend_cost, pad)
                if path is None:
                    continue
                crossings = _seg_cross(path[:-1], path[1:], S_all)
                length = sum(float(np.hypot(*(np.array(w) - u))) for u, w in zip(path[:-1], path[1:]))
                score = (crossings, len(path) - 2, length)
                if best is None or score < best[0]:
                    best = (score, path)
            if best is None:
                report['unresolved'] += 1
                continue
            steps = _round_corners(best[1], window, np.ones(len(window), bool))
            # write: midmarker -> corners -> metabolite
            del rx['segments'][key]
            prev = seg['from_node_id']
            for k, (pt, b1, b2) in enumerate(steps):
                last = k == len(steps) - 1
                to = seg['to_node_id'] if last else new_node(*pt)
                counter[0] += 1
                rx['segments'][f'er_s{counter[0]}'] = dict(from_node_id=prev, to_node_id=to, b1=dict(x=round(b1[0], 2), y=round(b1[1], 2)) if b1 else None,
                                                           b2=dict(x=round(b2[0], 2), y=round(b2[1], 2)) if b2 else None)
                prev = to
            corners_rounded = sum(1 for _, b1, _ in steps if b1)
            report['rerouted'] += 1
            report['smoothed_corners'] += corners_rounded
            report['bends_added'] += len(best[1]) - 2
            routed.append(dict(reaction=rid, metabolite=t.get('bigg_id'), kind='smoothed' if corners_rounded else 'avoided', bends=len(best[1]) - 2,
                               straight_px=round(float(np.hypot(*(q - p))), 1), routed_px=round(best[0][2], 1)))
    data[0].setdefault('metadata', {})['routed_edges'] = routed
    return report


# ------------------------------------------------------------------------------------------------ measurement (same code for any map)
def _path_points(rx, path, nodes):
    pts = [(nodes[path[0]]['x'], nodes[path[0]]['y'])]
    for a, b in zip(path, path[1:]):
        seg = next((s for s in rx['segments'].values() if (s['from_node_id'], s['to_node_id']) in ((a, b), (b, a))), None)
        if seg and seg.get('b1') and seg.get('b2'):
            fwd = (seg['from_node_id'], seg['to_node_id']) == (a, b)
            b1, b2 = (seg['b1'], seg['b2']) if fwd else (seg['b2'], seg['b1'])
            p0, p3 = np.array(pts[-1]), np.array([nodes[b]['x'], nodes[b]['y']])
            c1, c2 = np.array([b1['x'], b1['y']]), np.array([b2['x'], b2['y']])
            for t in np.linspace(0, 1, 7)[1:]:
                pts.append(tuple(((1 - t) ** 3 * p0 + 3 * (1 - t) ** 2 * t * c1 + 3 * (1 - t) * t ** 2 * c2 + t ** 3 * p3).tolist()))
        else:
            pts.append((nodes[b]['x'], nodes[b]['y']))
    return pts


def measure(data, skip_edges=None):
    """Bends, obstacle hits and crossings of the algorithm-placed connections of any map.json (expert connections are counted separately).
    A bend is a path node between the reaction and the metabolite (a corner point); an obstacle hit is a path segment crossing the box of a
    glyph or label that does not belong to its own two ends; a crossing is a pair of path segments from different reactions that cross."""
    import line_render
    meta = data[0].get('metadata', {})
    skip = {tuple(e) for e in (skip_edges if skip_edges is not None else meta.get('preserved_edges') or [])}
    skip |= {(p['reaction'], p['metabolite']) for p in (meta.get('source_paths') or [])}
    body = data[1]
    nodes, reactions = body['nodes'], body['reactions']
    rects, owner = obstacles(data)
    owner_arr = np.array(owner, dtype=object)
    own, algo_paths, all_segments = defaultdict(int), [], []
    for rx in reactions.values():
        for path in line_render._terminals(rx, nodes):
            if nodes[path[-1]].get('node_type') != 'metabolite':
                continue
            pts = _path_points(rx, path, nodes)
            is_expert = (rx['bigg_id'], nodes[path[-1]].get('bigg_id')) in skip
            segs = [(u[0], u[1], w[0], w[1]) for u, w in zip(pts[:-1], pts[1:])]
            for s in segs:
                all_segments.append((s, rx['bigg_id']))
            if is_expert:
                own['expert_edges'] += 1
                continue
            own['algorithm_edges'] += 1
            curved = sum(1 for u, w in zip(path[:-1], path[1:]) if any(s.get('b1') and s.get('b2') and {s['from_node_id'], s['to_node_id']} == {u, w} for s in rx['segments'].values()))
            own['bends'] += len(path) - 2 - curved  # a rounded corner is two path nodes joined by one curved segment: one bend
            own['rounded_corners'] += curved
            ends = (path[0], path[-1])
            mask = ~np.isin(owner_arr, ends)
            hit = 0
            for si, s in enumerate(segs):
                end_segment = si == 0 or si == len(segs) - 1
                f = _blocked if end_segment else (lambda p, q, r, pad: _hits(p, q, r, pad))
                hit += int(f((s[0], s[1]), (s[2], s[3]), rects[mask], 0.0).any())
            own['obstacle_hits'] += hit
            algo_paths.append((rx['bigg_id'], segs))
    # crossings: algorithm-placed segment against every segment of another reaction (each pair once)
    S = np.array([s for s, _ in all_segments], float).reshape(-1, 4)
    owners_s = np.array([r for _, r in all_segments], dtype=object)
    cross = 0
    for rid, segs in algo_paths:
        other = S[owners_s != rid]
        for s in segs:
            cross += _seg_cross((s[0], s[1]), (s[2], s[3]), other)
    own['crossings_with_other_reactions'] = cross
    return dict(own)
