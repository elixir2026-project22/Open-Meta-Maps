"""Label placement on the FINAL curves.

`drawing.export` (hashed) places every label once, before any curve exists, so after the V7 Bezier connections a label can sit on a line.  This module moves
the labels of the finished `map.json` (`label_x` / `label_y` of metabolites and reactions) among candidate positions around their glyph (right, left, above,
below, the diagonals, and where it is now) and keeps, for each label, the position with the lowest cost.  A glyph never moves; only text moves.

  cost = 10 x (sampled curve pieces of OTHER connections that cross the label box) + 3 x (pieces of its own connections)
       + 10 x (overlap with another label) + 8 x (overlap with another glyph)
Greedy, multiple passes, the worst labels first, ties keep the current position then prefer 'right'.  Deterministic.

At high local density (many labels in a region), low-importance labels (currency metabolites like ATP, NAD) can be hidden
to reduce clutter.  This is optional and controlled by the `hide_currency` parameter.
"""
import re
import numpy as np
from collections import defaultdict

PADX, PADY_UP, PADY_DOWN = 2.0, 14.0, 4.0
CELL = 150.0
DENSITY_THRESHOLD = 6  # if more than this many labels in a cell, consider hiding currency
# Labels further than this from their glyph are treated as detached (must re-place).
# Intentional side/diagonal candidates stay well below this; rigid-restore bugs were 300–1500 px.
MAX_ANCHOR_DIST = 100.0
CURRENCY_PATTERNS = {'atp', 'adp', 'amp', 'nad', 'nadh', 'nadp', 'nadph', 'coa', 'fad', 'fadh2',
                     'h2o', 'h', 'pi', 'ppi', 'co2', 'o2', 'nh4', 'gtp', 'gdp', 'gmp', 'ctp', 'cdp',
                     'cmp', 'utp', 'udp', 'ump', 'accoa', 'succoa', 'malacp', 'malcoa'}


def _rect(lx, ly, w):
    return (lx - PADX, ly - PADY_UP, lx + w + PADX, ly + PADY_DOWN)


def _cells(lo, hi):
    for gx in range(int(lo[0] // CELL), int(hi[0] // CELL) + 1):
        for gy in range(int(lo[1] // CELL), int(hi[1] // CELL) + 1):
            yield gx, gy


def _seg_rect(S, r):
    """Boolean per segment [k,4]: does it touch the rectangle r?"""
    x0, y0, x1, y1 = r
    ax, ay, bx, by = S[:, 0], S[:, 1], S[:, 2], S[:, 3]
    dx, dy = bx - ax, by - ay
    tmin, tmax = np.zeros(len(S)), np.ones(len(S))
    ok = np.ones(len(S), bool)
    for p, d, lo, hi in ((ax, dx, x0, x1), (ay, dy, y0, y1)):
        zero = np.abs(d) < 1e-9
        with np.errstate(divide='ignore', invalid='ignore'):
            t1, t2 = (lo - p) / d, (hi - p) / d
        ok &= np.where(zero, (p >= lo) & (p <= hi), True)
        tmin = np.where(zero, tmin, np.maximum(tmin, np.minimum(t1, t2)))
        tmax = np.where(zero, tmax, np.minimum(tmax, np.maximum(t1, t2)))
    return ok & (tmin <= tmax)


def _is_currency_label(bigg_id):
    """Check if a metabolite is a common currency/cofactor that can be hidden in dense regions."""
    if not bigg_id:
        return False
    base = re.sub(r'_[a-z][0-9]?$', '', bigg_id.lower())
    return base in CURRENCY_PATTERNS


def optimize(data, passes=10, hide_currency=True):
    import edge_routing as er
    import line_render
    import drawing
    font = drawing.FONT
    body = data[1]
    nodes, reactions = body['nodes'], body['reactions']
    # all connection pieces, with the node ids of their two ends
    pieces, piece_ends = [], []
    for rx in reactions.values():
        for path in line_render._terminals(rx, nodes):
            if nodes[path[-1]].get('node_type') != 'metabolite':
                continue
            pts = np.array(er._path_points(rx, path, nodes), float)
            for a, b in zip(pts[:-1], pts[1:]):
                pieces.append((a[0], a[1], b[0], b[1]))
                piece_ends.append((path[0], path[-1]))
    S = np.array(pieces, float).reshape(-1, 4)
    grid = defaultdict(list)
    for k, s in enumerate(S):
        for c in _cells(np.minimum(s[:2], s[2:]), np.maximum(s[:2], s[2:])):
            grid[c].append(k)
    # glyph boxes
    glyphs = []
    for nid, n in nodes.items():
        if n.get('node_type') in ('metabolite', 'midmarker'):
            glyphs.append((nid, (n['x'] - 7, n['y'] - 7, n['x'] + 7, n['y'] + 7)))
    glyph_grid = defaultdict(list)
    for i, (nid, r) in enumerate(glyphs):
        for c in _cells(r[:2], r[2:]):
            glyph_grid[c].append(i)
    # labels: (owner node id, container, widths, anchor)
    labels = []
    for nid, n in nodes.items():
        if n.get('node_type') == 'metabolite' and 'label_x' in n:
            bigg_id = n.get('bigg_id', '')
            labels.append(dict(owner=nid, box=n, w=font.getlength(bigg_id), cx=n['x'], cy=n['y'],
                               is_currency=_is_currency_label(bigg_id), bigg_id=bigg_id))
    for rx in reactions.values():
        mid = next((s['from_node_id'] for s in rx['segments'].values() if nodes[s['from_node_id']].get('node_type') == 'midmarker'), None)
        if mid is not None:
            labels.append(dict(owner=mid, box=rx, w=font.getlength(rx['bigg_id']), cx=nodes[mid]['x'], cy=nodes[mid]['y'],
                               is_currency=False, bigg_id=rx['bigg_id']))
    label_grid = defaultdict(set)

    def put(i, rect):
        for c in _cells(rect[:2], rect[2:]):
            label_grid[c].add(i)

    def drop(i, rect):
        for c in _cells(rect[:2], rect[2:]):
            label_grid[c].discard(i)
    current = []
    for i, L in enumerate(labels):
        r = _rect(L['box']['label_x'], L['box']['label_y'], L['w'])
        current.append(r)
        put(i, r)

    def near_candidates(L):
        x, y, w = L['cx'], L['cy'], L['w']
        # Positions beside the glyph (never the current far-away spot).
        return [
                # First ring: close to glyph
                (x + 10, y + 4), (x - 10 - w, y + 4), (x - w / 2, y - 9), (x - w / 2, y + 19),
                (x + 8, y - 8), (x + 8, y + 17), (x - 8 - w, y - 8), (x - 8 - w, y + 17),
                # Second ring: further out for crowded places
                (x + 26, y + 4), (x - 26 - w, y + 4), (x - w / 2, y - 24), (x - w / 2, y + 36),
                (x + 22, y - 22), (x + 22, y + 32), (x - 22 - w, y - 22), (x - 22 - w, y + 32),
                # Third ring: even further for very dense regions
                (x + 40, y + 4), (x - 40 - w, y + 4), (x - w / 2, y - 38), (x - w / 2, y + 52)]

    def candidates(L):
        return [(L['box']['label_x'], L['box']['label_y'])] + near_candidates(L)

    def _anchor_dist(L, lx, ly):
        return float(np.hypot(lx - L['cx'], ly - L['cy']))

    def cost(i, rect, lx=None, ly=None):
        L = labels[i]
        c = 0.0
        if lx is None:
            lx, ly = L['box']['label_x'], L['box']['label_y']
        # Detached labels (far from their glyph) are always costly so the greedy
        # pass will pull them back even when they sit in empty space.
        if _anchor_dist(L, lx, ly) > MAX_ANCHOR_DIST:
            c += 50.0
        seen = set()
        for cell in _cells(rect[:2], rect[2:]):
            seen.update(grid.get(cell, ()))
        if seen:
            idx = np.array(sorted(seen), int)
            hit = _seg_rect(S[idx], rect)
            for k in idx[hit]:
                c += 3.0 if L['owner'] in piece_ends[k] else 10.0
        near = set()
        for cell in _cells(rect[:2], rect[2:]):
            near.update(label_grid.get(cell, ()))
        for j in near:
            if j != i:
                r = current[j]
                if r[0] < rect[2] and rect[0] < r[2] and r[1] < rect[3] and rect[1] < r[3]:
                    c += 10.0
        gl = set()
        for cell in _cells(rect[:2], rect[2:]):
            gl.update(glyph_grid.get(cell, ()))
        for g in gl:
            nid, r = glyphs[g]
            if nid != L['owner'] and r[0] < rect[2] and rect[0] < r[2] and r[1] < rect[3] and rect[1] < r[3]:
                c += 8.0
        return c
    before = sum(cost(i, current[i]) > 0 for i in range(len(labels)))
    moved = 0
    hidden = 0
    reanchored = 0

    # Identify dense regions for potential currency hiding
    dense_cells = set()
    if hide_currency:
        cell_counts = defaultdict(int)
        for i, L in enumerate(labels):
            cx, cy = L['cx'], L['cy']
            cell = (int(cx // CELL), int(cy // CELL))
            cell_counts[cell] += 1
        dense_cells = {c for c, cnt in cell_counts.items() if cnt > DENSITY_THRESHOLD}

    # Safety: snap any already-detached label onto the best free near-glyph spot
    # before the conflict-driven passes (root-cause stages should prevent this).
    for i, L in enumerate(labels):
        if L['box'].get('label_hidden'):
            continue
        if _anchor_dist(L, L['box']['label_x'], L['box']['label_y']) <= MAX_ANCHOR_DIST:
            continue
        best = None
        for lx, ly in near_candidates(L):
            rect = _rect(lx, ly, L['w'])
            c = cost(i, rect, lx=lx, ly=ly)
            if best is None or c < best[0] - 1e-9:
                best = (c, lx, ly, rect)
        if best is not None:
            _, lx, ly, rect = best
            drop(i, current[i])
            current[i] = rect
            put(i, rect)
            L['box']['label_x'], L['box']['label_y'] = round(float(lx), 2), round(float(ly), 2)
            moved += 1
            reanchored += 1

    for p in range(passes):
        order = sorted(range(len(labels)), key=lambda i: -cost(i, current[i]))
        changed = 0
        for i in order:
            L = labels[i]
            c0 = cost(i, current[i])
            if c0 == 0:
                continue
            best = (c0, 0, None)
            for k, (lx, ly) in enumerate(candidates(L)[1:], start=1):
                rect = _rect(lx, ly, L['w'])
                c = cost(i, rect, lx=lx, ly=ly)
                if c < best[0] - 1e-9:
                    best = (c, k, (lx, ly, rect))
            if best[2] is not None:
                lx, ly, rect = best[2]
                drop(i, current[i])
                current[i] = rect
                put(i, rect)
                L['box']['label_x'], L['box']['label_y'] = round(float(lx), 2), round(float(ly), 2)
                changed += 1
                moved += 1
        if not changed:
            break

    # After placement optimization, hide currency labels in still-conflicted dense regions
    after = sum(cost(i, current[i]) > 0 for i in range(len(labels)))
    if hide_currency and after > 0:
        for i, L in enumerate(labels):
            if not L.get('is_currency'):
                continue
            c = cost(i, current[i])
            if c <= 0:
                continue
            # Check if in a dense region
            cx, cy = L['cx'], L['cy']
            cell = (int(cx // CELL), int(cy // CELL))
            if cell in dense_cells:
                # Park on the glyph (schema requires label_x/y). label_hidden is
                # stripped before Escher export; SVG skips co-located currency text.
                drop(i, current[i])
                L['box']['label_x'] = round(float(L['cx']), 2)
                L['box']['label_y'] = round(float(L['cy']), 2)
                L['box']['label_hidden'] = True
                current[i] = (0.0, 0.0, 0.0, 0.0)
                hidden += 1
        after = sum(
            1 for i in range(len(labels))
            if not labels[i]['box'].get('label_hidden') and cost(i, current[i]) > 0
        )

    return dict(labels=len(labels), with_conflict_before=before, with_conflict_after=after,
                moves=moved, reanchored=reanchored, hidden_currency=hidden, dense_cells=len(dense_cells))


def reanchor_distant(data, max_dist=MAX_ANCHOR_DIST):
    """Final safety pass: re-place any label whose offset from its glyph exceeds ``max_dist``.

    Prefer fixing the stage that moved the node without the label; call this after late
    geometry (rigid restore, overlap clear, squeeze) so a detached label cannot survive.
    Only detached labels move; already-anchored and hidden labels are left alone.
    """
    import drawing
    font = drawing.FONT
    body = data[1]
    nodes, reactions = body['nodes'], body['reactions']
    # Build a lightweight label list with anchors.
    items = []
    for nid, n in nodes.items():
        if n.get('node_type') != 'metabolite' or 'label_x' not in n or n.get('label_hidden'):
            continue
        items.append(dict(box=n, w=font.getlength(n.get('bigg_id', '')),
                          cx=n['x'], cy=n['y']))
    for rx in reactions.values():
        if 'label_x' not in rx or rx.get('label_hidden'):
            continue
        mid = next((s['from_node_id'] for s in rx['segments'].values()
                    if nodes[s['from_node_id']].get('node_type') == 'midmarker'), None)
        if mid is None:
            continue
        items.append(dict(box=rx, w=font.getlength(rx['bigg_id']),
                          cx=nodes[mid]['x'], cy=nodes[mid]['y']))
    # Occupancy of already-good labels (and glyphs) so re-placed text avoids them.
    occupied = []
    for L in items:
        if float(np.hypot(L['box']['label_x'] - L['cx'], L['box']['label_y'] - L['cy'])) <= max_dist:
            occupied.append(_rect(L['box']['label_x'], L['box']['label_y'], L['w']))
    for n in nodes.values():
        if n.get('node_type') in ('metabolite', 'midmarker'):
            occupied.append((n['x'] - 7, n['y'] - 7, n['x'] + 7, n['y'] + 7))

    def _overlap(a, b):
        return a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]

    def _near(L):
        x, y, w = L['cx'], L['cy'], L['w']
        return [
            (x + 10, y + 4), (x - 10 - w, y + 4), (x - w / 2, y - 9), (x - w / 2, y + 19),
            (x + 8, y - 8), (x + 8, y + 17), (x - 8 - w, y - 8), (x - 8 - w, y + 17),
            (x + 26, y + 4), (x - 26 - w, y + 4), (x - w / 2, y - 24), (x - w / 2, y + 36),
            (x + 22, y - 22), (x + 22, y + 32), (x - 22 - w, y - 22), (x - 22 - w, y + 32),
            (x + 40, y + 4), (x - 40 - w, y + 4), (x - w / 2, y - 38), (x - w / 2, y + 52),
        ]

    fixed = 0
    for L in items:
        dist = float(np.hypot(L['box']['label_x'] - L['cx'], L['box']['label_y'] - L['cy']))
        if dist <= max_dist:
            continue
        best = None
        for lx, ly in _near(L):
            rect = _rect(lx, ly, L['w'])
            hits = sum(1 for o in occupied if _overlap(rect, o))
            score = (hits, float(np.hypot(lx - L['cx'], ly - L['cy'])))
            if best is None or score < best[0]:
                best = (score, lx, ly, rect)
        if best is None:
            continue
        _, lx, ly, rect = best
        L['box']['label_x'] = round(float(lx), 2)
        L['box']['label_y'] = round(float(ly), 2)
        occupied.append(rect)
        fixed += 1
    return dict(checked=len(items), reanchored=fixed, max_dist=max_dist)
