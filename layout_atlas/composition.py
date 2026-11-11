"""Final composition pass: modules that the layout scattered far from everything are brought back and the truly independent ones are arranged in order.

In iIT341 the module arrangement left 164 visually separate components on a canvas of 64,350 x 95,435 px: 69 single reactions and 33 pairs sat in far corners
(5th-nearest-reaction distance up to 31,128 px, median 296), the map used 0.5 % of a coarse grid.  No curve fitting can fix that; modules have to move.

Everything is decided from the GEM and the scene, never from names, and whole modules move rigidly (their inner geometry, the expert curves and protected
shared nodes are untouched):

  component    reactions joined through shared non-currency glyphs of the final scene (a connected visual component); its keys: the reaction points, the
               metabolite glyphs and the per-reaction currency glyphs of its reactions
  main         the biggest components (>= MAIN_MIN reactions): never moved
  related      a smaller component that shares a non-currency METABOLITE (the GEM, not the glyph) with a reaction of another component: it is moved next to
               those anchors if it is far from them (more than FAR_FACTOR x its own size + FAR_MIN px); nothing is added to the chemistry, the glyphs stay separate
               unless the instance step already shared them
  independent  shares no non-currency metabolite with any other reaction: arranged on shelves (rows of rectangles, equal gaps) beside the main network,
               grouped by compartment and sorted by size; the side with the smallest aspect-ratio distortion is used; no connection is drawn
Free-space search keeps CLEAR px from everything already placed (labels and curve channels included by the margin), so the pass does not create crowded places.
The pass reports what it moved and why; `quality_gate` / the final evaluation decide on the finished map (A/B switch OPM_COMPOSE=off).
"""
from collections import defaultdict

import numpy as np
from scipy.spatial import cKDTree

from core import chem, comp

MAIN_MIN = 20
FAR_FACTOR = 0.5
FAR_MIN = 700.0        # legacy fixed tolerance (OPM_COMPOSE_LONG=off)
TAU = 3.5              # a related component is a relocation candidate when its distance to its anchors exceeds TAU x the typical reaction spacing + FAR_FACTOR x its own size
import os
LONG = os.environ.get('OPM_COMPOSE_LONG', 'on') == 'on'
CLEAR = 110.0
SHELF_GAP = 120.0
ROW_WIDTH_MIN = 2000.0  # floor only; actual row width follows main-network size and target aspect
TARGET_ASPECT = 1.4    # target width/height ratio for the independent shelf block (slightly wider than tall)


def components(model, scene, cur):
    """Connected visual components: reaction ids -> dict(reactions, keys)."""
    parts = scene['participants']
    parent = {}

    def find(a):
        while parent.setdefault(a, a) != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for r in model['reactions']:
        rid = r['id']
        find('r:' + rid)
        for mid in r['metabolites']:
            key = parts[rid, mid]
            if key.startswith('side:'):
                continue
            parent[find(key)] = find('r:' + rid)
    groups = defaultdict(lambda: dict(reactions=[], keys=set()))
    for r in model['reactions']:
        rid = r['id']
        g = groups[find('r:' + rid)]
        g['reactions'].append(rid)
        g['keys'].add('r:' + rid)
        for mid in r['metabolites']:
            g['keys'].add(parts[rid, mid])
    return list(groups.values())


def _free_spot(local, target, tree, clear=CLEAR, tries=260):
    if tree is None:
        return target
    for i in range(tries):
        r = 0.0 if i == 0 else 45.0 * np.sqrt(i)
        a = i * 2.399963
        d = target + r * np.array([np.cos(a), np.sin(a)])
        if float(tree.query(local + d)[0].min()) >= clear:
            return d
    return None


def _aligned_spot(local, target, anchors, size, tree, centre):
    """Free place for a component next to its anchors.  Candidates: free spots on the 8 compass sides of the anchor cluster (distance = half the component size + a gap)
    plus the plain nearest free spot; score = distance to the anchor cluster + a bonus for sharing a centre line (x or y) with the nearest anchor, so the module lines up
    with what it belongs to instead of sitting at a random angle."""
    nearest = anchors[int(np.argmin(np.hypot(*(anchors - centre).T)))]
    r = size / 2.0 + 160.0
    cands = []
    plain = _free_spot(local, target, tree)
    if plain is not None:
        cands.append(plain)
    for k in range(8):
        a = k * np.pi / 4
        d = _free_spot(local, target + r * np.array([np.cos(a), np.sin(a)]), tree)
        if d is not None:
            cands.append(d)
    if not cands:
        return None

    def score(c):
        dx, dy = abs(c[0] - nearest[0]), abs(c[1] - nearest[1])
        return float(np.hypot(*(c - target))) + 0.6 * min(dx, dy)
    return min(cands, key=score)


def compose(model, scene, blocks, cur):
    """-> (new scene, report).  Moves small components; `scene` is not modified."""
    scene = dict(scene, positions={k: list(v) for k, v in scene['positions'].items()})
    pos, parts = scene['positions'], scene['participants']
    comps = components(model, scene, cur)
    rx = {r['id']: r for r in model['reactions']}
    metab_users = defaultdict(set)         # core metabolite -> reaction ids
    for r in model['reactions']:
        for mid in r['metabolites']:
            if mid not in cur:
                metab_users[mid].add(r['id'])
    comp_of = {rid: i for i, c in enumerate(comps) for rid in c['reactions']}
    main = [i for i, c in enumerate(comps) if len(c['reactions']) >= MAIN_MIN]
    report = dict(components=len(comps), main=len(main), moved_related=0, arranged_independent=0, kept_near=0, failed=0, moves=[])
    if not main:
        return scene, report
    movable_idx = {i for i in range(len(comps)) if i not in set(main)}
    rpts = np.array([pos[k] for k in pos if k.startswith('r:')], float)
    spacing = float(np.median(cKDTree(rpts).query(rpts, k=2)[0][:, 1])) if len(rpts) > 2 else 100.0
    report['spacing_px'] = round(spacing, 1)

    def tree_without(skip_keys):
        pts = [pos[k] for k in pos if k not in skip_keys]
        return cKDTree(np.array(pts, float)) if pts else None
    cur_keys_of = {i: {k for k in comps[i]['keys'] if k in pos} for i in range(len(comps))}
    # classify the movable components: a component is "related" only when it shares a core metabolite with a SETTLED component (the main network, or a
    # component that has already been brought next to it); components whose only partners are other stray components are settled in later rounds, and those
    # that never reach the main network are independent
    settled = set(main)
    related, independent = [], []
    pending = sorted(movable_idx, key=lambda i: (-len(comps[i]['reactions']), comps[i]['reactions'][0]))
    anchors_of = {}
    for i in pending:
        anchors_of[i] = {(mid, u) for rid in comps[i]['reactions'] for mid in rx[rid]['metabolites'] if mid not in cur for u in metab_users[mid] if comp_of[u] != i}
    order_rounds = []
    changed = True
    while changed and pending:
        changed = False
        for i in list(pending):
            if any(comp_of[u] in settled for _, u in anchors_of[i]):
                order_rounds.append(i)
                settled.add(i)
                pending.remove(i)
                changed = True
    rank = {i: -1 for i in main}
    rank.update({i: n for n, i in enumerate(order_rounds)})
    related = [(i, {(m, u) for m, u in anchors_of[i] if comp_of[u] in rank and rank[comp_of[u]] < rank[i]}) for i in order_rounds]
    independent = [(i, anchors_of[i]) for i in pending]
    # one occupancy tree for everything that is currently placed; updated as components move
    placed_pts = {k: list(v) for k, v in pos.items() if not k.startswith('side:')}

    def build_tree(exclude):
        pts = np.array([v for k, v in placed_pts.items() if k not in exclude], float)
        return cKDTree(pts) if len(pts) else None

    def move(i, delta):
        for k in cur_keys_of[i]:
            for table in (pos, placed_pts):
                if k in table:
                    table[k] = [table[k][0] + float(delta[0]), table[k][1] + float(delta[1])]
    # ---- related components: next to their anchors
    for i, anchors in related:
        keys = cur_keys_of[i]
        local = np.array([pos[k] for k in keys], float)
        size = float(np.hypot(*(local.max(0) - local.min(0))))
        centre = local.mean(0)
        anchor_pts = np.array([pos['r:' + u] for _, u in anchors if 'r:' + u in pos], float)
        if not len(anchor_pts):
            continue
        # the anchors that are nearest to the component decide where it goes
        d = np.hypot(*(anchor_pts - centre).T)
        near = anchor_pts[np.argsort(d)[:max(1, min(6, len(anchor_pts)))]]
        target_c = near.mean(0)
        dist_now = float(np.hypot(*(centre - target_c)))
        limit = (FAR_FACTOR * size + TAU * spacing) if LONG else (1.5 * size + FAR_MIN)
        if dist_now <= limit:
            report['kept_near'] += 1
            continue
        tree = build_tree(keys)
        if LONG:
            spot = _aligned_spot(local - centre, target_c, near, size, tree, centre)
        else:
            spot = _free_spot(local - centre, target_c, tree)
        if spot is None or (LONG and float(np.hypot(*(spot - target_c))) > 0.7 * dist_now):
            report['failed'] += 1          # no free place clearly closer to the anchors
            continue
        delta = spot - centre
        move(i, delta)
        report['moved_related'] += 1
        if len(report['moves']) < 400:
            report['moves'].append(dict(kind='related', reactions=comps[i]['reactions'][:4], moved_px=round(float(np.hypot(*delta)), 1), anchors=sorted({u for _, u in anchors})[:4]))
    # ---- independent components: shelves beside the main network
    if independent:
        main_pts = np.array([placed_pts[k] for i in main for k in cur_keys_of[i] if k in placed_pts], float)
        lo, hi = main_pts.min(0), main_pts.max(0)
        W = max(float(hi[0] - lo[0]), ROW_WIDTH_MIN)
        order = sorted(independent, key=lambda t: (_compartment(comps[t[0]], rx), -len(comps[t[0]]['reactions']), comps[t[0]]['reactions'][0]))
        boxes = []
        labels = scene.get('labels', {})
        try:
            import drawing
            text_w = lambda s: float(drawing.FONT.getlength(s))
        except Exception:
            text_w = lambda s: 6.5 * len(s)
        pads = {}
        for i, _ in order:
            ks = [k for k in cur_keys_of[i] if k in pos]
            arr = np.array([pos[k] for k in ks], float)
            # the visible rectangle: glyph points plus the widest label (a label sits beside its glyph) and a margin for the curves
            lw = min(max([text_w(str(labels.get(k, ''))) for k in ks] + [0.0]), 220.0)
            pads[i] = lw
            boxes.append((i, arr.min(0), arr.max(0) - arr.min(0) + np.array([60.0 + lw, 60.0])))
        # row width: aim for a TARGET_ASPECT block (wider than tall), not limited by main network width
        total_area = sum(b[2][0] * b[2][1] for b in boxes)
        # solve for row_w such that area / row_w * row_w ~ TARGET_ASPECT -> row_w ~ sqrt(area * TARGET_ASPECT)
        row_w = max(np.sqrt(total_area * TARGET_ASPECT), max(b[2][0] for b in boxes) + SHELF_GAP, ROW_WIDTH_MIN * 0.5)
        x = y = 0.0
        row_h = 0.0
        layout = []
        prev_group = None
        for i, bmin, size in boxes:
            g = _compartment(comps[i], rx)
            if x > 0 and (x + size[0] > row_w or (prev_group is not None and g != prev_group and x > row_w * 0.5)):
                x, y, row_h = 0.0, y + row_h + SHELF_GAP, 0.0
            layout.append((i, bmin, np.array([x, y])))
            x += size[0] + SHELF_GAP
            row_h = max(row_h, size[1])
            prev_group = g
        block_w, block_h = row_w, y + row_h
        # side with combined aspect ratio closest to 1:1 (square-ish canvas)
        main_w, main_h = float(hi[0] - lo[0]), float(hi[1] - lo[1])
        # each option: (x_origin, y_origin, combined_aspect_ratio)
        options = {'below': (lo[0], hi[1] + 400.0, max(main_w, block_w) / max(main_h + block_h + 400.0, 1.0)),
                   'above': (lo[0], lo[1] - 400.0 - block_h, max(main_w, block_w) / max(main_h + block_h + 400.0, 1.0)),
                   'right': (hi[0] + 400.0, lo[1], (main_w + block_w + 400.0) / max(main_h, block_h, 1.0)),
                   'left': (lo[0] - 400.0 - block_w, lo[1], (main_w + block_w + 400.0) / max(main_h, block_h, 1.0))}
        # prefer aspect ratio closest to 1.0; favor horizontal placements (right/left) when main is already tall
        side = min(options, key=lambda s: abs(np.log(max(options[s][2], 1e-6))))
        ox, oy = options[side][0], options[side][1]
        for i, bmin, off in layout:
            ks = [k for k in cur_keys_of[i] if k in pos]
            arr = np.array([pos[k] for k in ks], float)
            delta = np.array([ox, oy]) + off - arr.min(0) + np.array([30.0 + pads.get(i, 0.0) / 2.0, 30.0])
            move(i, delta)
            report['arranged_independent'] += 1
        # after-the-fact check: no two shelf rectangles (and none of them with the main network's rectangle) may overlap
        rects = []
        for i, _, _ in layout:
            arr = np.array([pos[k] for k in cur_keys_of[i] if k in pos], float)
            rects.append((arr.min(0)[0] - 30, arr.min(0)[1] - 30, arr.max(0)[0] + pads.get(i, 0.0) + 30, arr.max(0)[1] + 30))
        overlaps = sum(1 for a in range(len(rects)) for b in range(a + 1, len(rects)) if rects[a][0] < rects[b][2] and rects[b][0] < rects[a][2] and rects[a][1] < rects[b][3] and rects[b][1] < rects[a][3])
        main_rect = (lo[0], lo[1], hi[0], hi[1])
        with_main = sum(1 for r in rects if r[0] < main_rect[2] and main_rect[0] < r[2] and r[1] < main_rect[3] and main_rect[1] < r[3])
        report['shelf_check'] = dict(rect_overlaps=overlaps, overlaps_with_main_network=with_main)
        report['shelf'] = dict(side=side, rows=int(round(block_h / max(row_h + SHELF_GAP, 1))) , width=round(block_w), height=round(block_h), components=len(layout))
    scene['composition'] = report
    return scene, report


def _compartment(c, rx):
    cnt = defaultdict(int)
    for rid in c['reactions']:
        for mid in rx[rid]['metabolites']:
            cnt[comp(mid)] += 1
    return max(cnt, key=cnt.get) if cnt else ''
