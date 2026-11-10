"""Pathway-level organisation (V7.3 P0): functional regions and an initial coarse arrangement of the rigid modules.

The module arrangement before this step comes from the block packing (ELK) and the completion layouts; the pieces that were not packed (generated reactions) and the
pieces that were (expert blocks) end up where each stage left them, with no common idea of which belong together.  This pass looks at the whole scene once:

  units     an expert block (rigid) or a group of algorithm-drawn reactions that share glyphs (rigid as laid out); each unit has a rectangle (final glyph positions + margin)
  links     two units are linked by every non-currency metabolite their reactions share (hub metabolites used by many units are damped; > MAX_UNITS_PER_METABOLITE are ignored)
  regions   communities of the unit graph (Louvain, seeded, so repeatable); a region is a functional area - a soft grouping: a unit belongs to one region here, but the
            region object only steers the arrangement (within-region links pull a little harder) and is reported; nothing forces regions to be exclusive in the map
  arrange   a few hundred damped force steps move whole units rigidly: linked units are pulled to a rest distance (half their sizes + a gap), overlapping rectangles are
            pushed apart, the biggest unit stays where it is; then linked units whose centres almost line up are snapped onto a common centre line (weak alignment)
Only translations of whole units; glyphs shared by two units belong to the unit that holds most of their reactions (the others' connections simply get shorter or longer
and are re-fitted as Bezier curves later).  GEM content, protected connections and the inside of every expert block are untouched.
The result is judged on the finished map by the usual evaluation (A/B switch OPM_REGIONS=on|off); `scene['regions']` carries the region objects for later steps.
"""
from collections import defaultdict

import numpy as np

MAX_UNITS_PER_METABOLITE = 6
MARGIN = 50.0
GAP = 140.0
CLEAR = 100.0
ITER = 260
ALIGN_FRACTION = 0.30


def _owner_reaction(key):
    p = key.split(':')
    if p[0] in ('r', ):
        return p[1]
    if p[0] in ('side', 'inst'):
        return p[1]
    return None


def build_units(model, scene, blocks):
    pos, parts = scene['positions'], scene['participants']
    rx = {r['id']: r for r in model['reactions']}
    block_unit = {}
    for b in blocks:
        for rid in b['reaction_ids']:
            block_unit[rid] = 'b:' + b['id']
    parent = {}

    def find(a):
        while parent.setdefault(a, a) != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    free = [rid for rid in rx if rid not in block_unit]
    for rid in free:
        find(rid)
        for mid in rx[rid]['metabolites']:
            k = parts[rid, mid]
            if k.startswith('side:') or k.startswith('inst:') or k.startswith('copy:'):
                continue
            for other, o in ((o, rx[o]) for o in free if o != rid and mid in rx[o]['metabolites'] and parts[o, mid] == k):
                parent[find(other)] = find(rid)
    unit_of = dict(block_unit)
    for rid in free:
        unit_of[rid] = 'g:' + find(rid)
    holders = defaultdict(list)
    for (rid, mid), k in parts.items():
        holders[k].append(rid)
    keys_of = defaultdict(set)
    for k in pos:
        u = None
        if k.startswith('copy:'):
            u = 'b:' + k.split(':')[1]
        else:
            rid = _owner_reaction(k)
            if rid is not None:
                u = unit_of.get(rid)
            elif k in holders:
                count = defaultdict(int)
                for r in holders[k]:
                    count[unit_of[r]] += 1
                u = sorted(count, key=lambda x: (-count[x], x))[0]
        if u is not None:
            keys_of[u].add(k)
    reactions_of = defaultdict(set)
    for rid, u in unit_of.items():
        reactions_of[u].add(rid)
    return unit_of, keys_of, reactions_of


def organize(model, scene, blocks, cur):
    import networkx as nx
    from networkx.algorithms.community import louvain_communities
    scene = dict(scene, positions={k: list(v) for k, v in scene['positions'].items()})
    pos = scene['positions']
    unit_of, keys_of, reactions_of = build_units(model, scene, blocks)
    uids = sorted(u for u in keys_of if keys_of[u])
    index = {u: i for i, u in enumerate(uids)}
    n = len(uids)
    report = dict(units=n, expert_units=sum(u.startswith('b:') for u in uids), generated_units=sum(u.startswith('g:') for u in uids))
    if n < 3:
        return scene, report
    pts = [np.array([pos[k] for k in keys_of[u]], float) for u in uids]
    c0 = np.array([(p.min(0) + p.max(0)) / 2 for p in pts])
    wh = np.array([(p.max(0) - p.min(0)) + 2 * MARGIN for p in pts])
    # links through shared non-currency metabolites
    users = defaultdict(set)
    for r in model['reactions']:
        for mid in r['metabolites']:
            if mid not in cur:
                users[mid].add(unit_of[r['id']])
    W = defaultdict(float)
    for mid, us in users.items():
        us = sorted(us)
        if not 2 <= len(us) <= MAX_UNITS_PER_METABOLITE:
            continue
        for a in range(len(us)):
            for b in range(a + 1, len(us)):
                W[index[us[a]], index[us[b]]] += 1.0 / (len(us) - 1)
    G = nx.Graph()
    G.add_nodes_from(range(n))
    for (a, b), w in W.items():
        G.add_edge(a, b, weight=w)
    comms = louvain_communities(G, weight='weight', seed=0) if G.number_of_edges() else [{i} for i in range(n)]
    region_of = {}
    regions = []
    for ri, cset in enumerate(sorted(comms, key=lambda s: (-len(s), min(s)))):
        for i in cset:
            region_of[i] = ri
        regions.append(dict(id=ri, units=[uids[i] for i in sorted(cset)], reactions=sum(len(reactions_of[uids[i]]) for i in cset)))
    links = [(a, b, w * (1.5 if region_of[a] == region_of[b] else 1.0)) for (a, b), w in W.items()]
    sizes = np.array([len(reactions_of[u]) for u in uids])
    fixed = int(np.argmax(sizes))
    c = c0.copy()
    diag = np.hypot(wh[:, 0], wh[:, 1])

    def link_length(cc):
        return float(sum(w * np.hypot(*(cc[a] - cc[b])) for a, b, w in links))
    before = link_length(c)
    wmax = max((w for _, _, w in links), default=1.0)
    for it in range(ITER):
        cool = 1.0 - 0.7 * it / ITER
        disp = np.zeros_like(c)
        for a, b, w in links:
            d = c[b] - c[a]
            dist = float(np.hypot(*d)) or 1e-9
            rest = 0.28 * (diag[a] + diag[b]) + GAP
            if dist > rest:
                f = 0.06 * (w / wmax) ** 0.5 * (dist - rest)
                disp[a] += f * d / dist
                disp[b] -= f * d / dist
        dx = c[:, None, 0] - c[None, :, 0]
        dy = c[:, None, 1] - c[None, :, 1]
        ox = (wh[:, None, 0] + wh[None, :, 0]) / 2 + CLEAR - np.abs(dx)
        oy = (wh[:, None, 1] + wh[None, :, 1]) / 2 + CLEAR - np.abs(dy)
        hit = (ox > 0) & (oy > 0)
        np.fill_diagonal(hit, False)
        for i, j in zip(*np.nonzero(np.triu(hit))):
            if ox[i, j] < oy[i, j]:
                s = np.sign(dx[i, j]) or 1.0
                push = np.array([s * ox[i, j] * 0.5, 0.0])
            else:
                s = np.sign(dy[i, j]) or 1.0
                push = np.array([0.0, s * oy[i, j] * 0.5])
            disp[i] += push
            disp[j] -= push
        disp[fixed] = 0
        step = np.clip(disp, -35 * cool - 5, 35 * cool + 5)
        c += step
    # weak alignment of linked units whose centres almost line up on an axis
    aligned = 0
    order = sorted(links, key=lambda t: -t[2])
    for a, b, w in order:
        if w < np.median([x[2] for x in links]):
            break
        mover = b if sizes[b] <= sizes[a] else a
        other = a if mover == b else b
        if mover == fixed:
            continue
        for axis in (0, 1):
            gap = c[other][axis] - c[mover][axis]
            if 8.0 < abs(gap) < ALIGN_FRACTION * max(wh[mover][axis], wh[other][axis]):
                trial = c.copy()
                trial[mover][axis] += gap
                ox2 = (wh[mover, 0] + wh[:, 0]) / 2 + CLEAR - np.abs(trial[:, 0] - trial[mover, 0])
                oy2 = (wh[mover, 1] + wh[:, 1]) / 2 + CLEAR - np.abs(trial[:, 1] - trial[mover, 1])
                bad = (ox2 > 0) & (oy2 > 0)
                bad[mover] = False
                if not bad.any():
                    c = trial
                    aligned += 1
    after = link_length(c)
    delta = c - c0
    moved = 0
    for u, i in index.items():
        d = delta[i]
        if float(np.hypot(*d)) < 1.0:
            continue
        moved += 1
        for k in keys_of[u]:
            pos[k] = [pos[k][0] + float(d[0]), pos[k][1] + float(d[1])]
    scene['regions'] = regions
    report.update(links=len(links), regions=len(regions), region_sizes=[len(r['units']) for r in regions][:12], moved_units=moved, aligned_pairs=aligned,
                  weighted_link_length=[round(before), round(after)], max_move_px=round(float(np.hypot(*delta.T).max())), fixed_unit=uids[fixed])
    return scene, report
