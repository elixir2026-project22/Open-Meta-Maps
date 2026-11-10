"""Pathway backbone layout for small groups of reactions (V7, first slice: the regions `context.py` builds).

A region is a bipartite graph: reactions and the shared glyphs (metabolite nodes) they have in common.  The old layout was a spring embedding; it knows
nothing about direction, so a chain became a tangle and three reactions sharing ATP and ADP crossed each other.  Here the region is drawn as a layered
pathway, like the expert maps draw a pathway:

  direction   a shared glyph is a substrate of a reaction (coefficient < 0) or a product (> 0); an edge runs substrate -> reaction -> product
              (cycles are cut by reversing the edge that closes them, so a reaction that runs the other way sits beside its partner, not after it)
  layers      longest-path layering along that direction -> the main chain runs along +x, products to the right
  order       inside a layer the reactions / glyphs are ordered by the barycentre of their neighbours (3 sweeps), which removes the crossings of
              small bipartite groups (e.g. ATP / ADP shared by ATPM, ADK1, ATPS4r become a ladder between two columns)
  branches    rows spread above and below the main chain; a graph without any shared glyph keeps one reaction per row
Only the coordinates of this region are computed; the caller places the region among the rest of the map.  Deterministic.
"""
import networkx as nx
import numpy as np

DX, DY = 150.0, 95.0


def layered(reactions, region, cur=None):
    """reactions: {id: reaction dict (metabolites {mid: coefficient})}; region: dict(reactions=[ids], shared={mid: [rids]}).
    Returns {'r:<id>': (x, y), 's:<mid>': (x, y)} in local coordinates (x along the direction of the pathway)."""
    members = list(region['reactions'])
    shared = {m: [r for r in rs if r in members] for m, rs in region['shared'].items() if len([r for r in rs if r in members]) >= 1}
    D = nx.DiGraph()
    for rid in members:
        D.add_node('r:' + rid)
    for mid, rids in shared.items():
        D.add_node('s:' + mid)
        for rid in rids:
            coeff = reactions[rid]['metabolites'][mid]
            if coeff < 0:
                D.add_edge('s:' + mid, 'r:' + rid)   # substrate -> reaction
            else:
                D.add_edge('r:' + rid, 's:' + mid)   # reaction -> product
    # cut cycles: reverse the closing edge of every cycle found by DFS (deterministic order)
    while True:
        try:
            cycle = nx.find_cycle(D)
        except nx.NetworkXNoCycle:
            break
        u, v = sorted(cycle, key=lambda e: (e[0], e[1]))[0]
        D.remove_edge(u, v)
        D.add_edge(v, u)
    layer = {}
    for n in nx.topological_sort(D):
        layer[n] = max([layer[p] + 1 for p in D.predecessors(n)], default=0)
    # keep reactions and shared glyphs on alternating columns (an edge always spans exactly one column)
    layers = {}
    for n, l in layer.items():
        layers.setdefault(l, []).append(n)
    order = {n: i for l in sorted(layers) for i, n in enumerate(sorted(layers[l]))}
    U = D.to_undirected()
    for sweep in range(3):
        for l in (sorted(layers) if sweep % 2 == 0 else sorted(layers, reverse=True)):
            nodes_l = layers[l]
            bary = {}
            for n in nodes_l:
                neigh = [order[m] for m in U.neighbors(n) if layer[m] != l]
                bary[n] = (sum(neigh) / len(neigh)) if neigh else order[n]
            for i, n in enumerate(sorted(nodes_l, key=lambda n: (bary[n], n))):
                order[n] = i
    pos = {}
    for l in sorted(layers):
        count = len(layers[l])
        for n in layers[l]:
            pos[n] = (l * DX, (order[n] - (count - 1) / 2.0) * DY)
    # isolated components of the same region side by side (a region without shared glyphs between some members)
    comps = [sorted(c) for c in nx.connected_components(U)]
    if len(comps) > 1:
        shift = 0.0
        for comp in sorted(comps, key=lambda c: c[0]):
            ys = [pos[n][1] for n in comp]
            span = max(ys) - min(ys)
            for n in comp:
                pos[n] = (pos[n][0], pos[n][1] + shift - min(ys))
            shift += span + DY * 1.4
    return pos


# ------------------------------------------------------------------------------------------------ main chain of an automatic module
CHAIN_MIN, CHAIN_MAX = 3, 16
HUB_REACTIONS = 3        # a metabolite used by more reactions than this inside the module does not link a chain (hubs / cofactor-like)
STEP = 230.0             # distance between consecutive reactions on the chain
ROW = 8                  # reactions per row; longer chains fold back (snake) so a module does not become a long line
ROW_GAP = 210.0


def _direction(rx):
    """+1 forward, -1 backward, 0 uncertain (reversible): from the SBML bounds, no FBA."""
    lo, hi = rx.get('lower_bound', 0), rx.get('upper_bound', 0)
    if lo >= 0 and hi > 0:
        return 1
    if hi <= 0 and lo < 0:
        return -1
    return 0


def find_chain(reactions, rids, cur):
    """The longest chain of reactions of one module linked by core (non-currency, non-hub) metabolites, 3..16 long; None if there is none.
    Returns dict(reactions=[ids in order], links=[metabolite id between consecutive reactions], direction=+1/-1/0).
    The order follows the chain; its direction is the majority of the reactions' allowed directions (bounds only): reversible reactions do not vote."""
    rids = list(rids)
    users = {}
    for rid in rids:
        for mid in reactions[rid]['metabolites']:
            if mid not in cur:
                users.setdefault(mid, []).append(rid)
    G = nx.Graph()
    G.add_nodes_from(rids)
    via = {}
    for mid, us in users.items():
        if 2 <= len(us) <= HUB_REACTIONS:
            for i in range(len(us)):
                for j in range(i + 1, len(us)):
                    if not G.has_edge(us[i], us[j]):
                        G.add_edge(us[i], us[j])
                        via[frozenset((us[i], us[j]))] = mid
    best = []
    for comp in nx.connected_components(G):
        if len(comp) < CHAIN_MIN:
            continue
        sub = G.subgraph(comp)
        start = sorted(comp)[0]
        far = max(nx.single_source_shortest_path_length(sub, start).items(), key=lambda kv: (kv[1], kv[0]))[0]
        dist = nx.single_source_shortest_path_length(sub, far)
        far2 = max(dist.items(), key=lambda kv: (kv[1], kv[0]))[0]
        path = nx.shortest_path(sub, far, far2)
        if len(path) > len(best):
            best = path
    if len(best) < CHAIN_MIN:
        return None
    best = best[:CHAIN_MAX]
    links = [via[frozenset((a, b))] for a, b in zip(best[:-1], best[1:])]
    votes = 0
    for k, rid in enumerate(best):
        d = _direction(reactions[rid])
        if d == 0:
            continue
        nxt = links[k] if k < len(links) else None
        prv = links[k - 1] if k > 0 else None
        coeff = reactions[rid]['metabolites']
        along = 0
        if nxt is not None:      # forward along the chain means: the link metabolite is a product of this reaction in its allowed direction
            along = 1 if coeff[nxt] * d > 0 else -1
        elif prv is not None:
            along = 1 if coeff[prv] * d < 0 else -1
        votes += along
    if votes < 0:
        best, links = best[::-1], links[::-1]
    return dict(reactions=best, links=links, direction=1 if votes > 0 else (-1 if votes < 0 else 0), votes=votes)


def chain_positions(chain):
    """Fixed coordinates of the chain: reactions every STEP along x, the link metabolite midway; folded into rows of ROW reactions (snake)."""
    pos = {}
    n = len(chain['reactions'])
    for k, rid in enumerate(chain['reactions']):
        row, col = divmod(k, ROW)
        x = col * STEP if row % 2 == 0 else (ROW - 1 - col) * STEP
        pos['r:' + rid] = (x, row * ROW_GAP)
    for k, mid in enumerate(chain['links']):
        a, b = pos['r:' + chain['reactions'][k]], pos['r:' + chain['reactions'][k + 1]]
        pos['m:' + mid] = ((a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0)
    return pos
