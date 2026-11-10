"""V2: several arrangements of the modules on the global canvas, one chosen by module-level metrics.

`global_layout.layout_blocks` placed every module (expert blocks and automatic modules, all rigid rectangles) with one spring layout.  For a
large model the result had one weak spot: a bundle of long vertical connections through the middle of the map.  Here the arrangement is
built three ways, each scored on what it will do to the final map, and the best one is the one drawn:

  spring         the V1 arrangement (the previous behaviour; always a candidate, so the choice can never be worse than V1 by these measures)
  stress         stress majorization on the module-connection graph: modules that are N links apart end up ~N module-widths apart
  stress_square  the same on a square canvas, with the attraction compressed (sqrt) so that highly shared metabolites do not pull everything
                 into the middle

The module-connection graph: nodes are modules, an edge joins two modules that share a non-currency metabolite; a metabolite shared by more
than `HUB` modules is ignored (ATP/H2O-like hubs must not decide where modules go) and a metabolite shared by k modules adds 1/(k-1) to each
pair.  `links` (unweighted, one per pair per metabolite) is what the metrics count.

Metrics per candidate (all estimated from module rectangles; no routing, no drawing):
  links_over_250 / links_over_1000   module links whose rectangles are further apart than that (an estimate of long cross-module edges)
  module_crossings                   module links whose centre-to-centre segment passes through a third module
  gap_weighted                       sum of link weight x rectangle gap (reported; NOT used alone to choose)
  canvas / compactness / aspect      bounding box, share of it covered by modules, width / height
The choice: candidates with an extreme aspect ratio or a much emptier canvas than the best one are discarded first; among the rest the one with
the fewest links over 1000, then over 250, then module crossings (each compared with a small tolerance) wins.  The reason is written down.
Expert geometry is untouched: only whole modules move.  Deterministic.
"""
import time
from collections import defaultdict

import numpy as np

ASPECT_RANGE = (0.55, 2.4)        # canvas width / height outside this is discarded
COMPACT_TOLERANCE = 0.75          # a candidate must cover at least this share of the best candidate's compactness
LONG_TOLERANCE = 0.05             # counts within 5 % are treated as equal and the next criterion decides


def module_links(groups, cur, hub):
    """groups: [{reactions, _rx}]  ->  (weights {(a, b): w}, links {(a, b): number of shared metabolites})."""
    owner = defaultdict(set)
    for gi, grp in enumerate(groups):
        for rid in grp['reactions']:
            for mid in grp['_rx'][rid]['metabolites']:
                if mid not in cur:
                    owner[mid].add(gi)
    weights, links = defaultdict(float), defaultdict(int)
    for mid, gs in owner.items():
        gs = sorted(gs)
        if 1 < len(gs) <= hub:
            w = 1.0 / (len(gs) - 1)
            for a in range(len(gs)):
                for c in range(a + 1, len(gs)):
                    weights[gs[a], gs[c]] += w
                    links[gs[a], gs[c]] += 1
    return weights, links


def stress_layout(sizes, weights, aspect=1.0, compress=1.0, iterations=300):
    """Stress majorization (SMACOF) of the module-connection graph.  Target distance between linked modules is the sum of their half
    diagonals plus a gap; unconnected components are placed one after another."""
    import networkx as nx
    from scipy.sparse.csgraph import shortest_path
    n = len(sizes)
    sizes = np.asarray(sizes, float)
    radius = np.hypot(sizes[:, 0], sizes[:, 1]) / 2
    W = np.zeros((n, n))
    for (a, b), w in weights.items():
        W[a, b] = W[b, a] = w ** compress if compress != 1.0 else w
    # a strongly shared pair is closer: length 1 / (1 + w)
    L = np.where(W > 0, 1.0 / (1.0 + W), 0.0)
    D = shortest_path(L, directed=False)
    finite = np.isfinite(D)
    big = D[finite].max() if finite.any() else 1.0
    D = np.where(finite, D, big * 1.5 + 1.0)
    unit = float(np.median(radius)) * 2.2 + 90.0
    target = D * unit + (radius[:, None] + radius[None, :]) * .8
    np.fill_diagonal(target, 0.0)
    # start from a deterministic spiral (no random numbers) stretched to the aspect ratio
    angle = np.arange(n) * 2.399963
    rad = np.sqrt(np.arange(n) + 1.0) * unit * .7
    xy = np.stack([np.cos(angle) * rad, np.sin(angle) * rad], axis=1)
    weight = 1.0 / np.maximum(target, 1.0) ** 2
    np.fill_diagonal(weight, 0.0)
    wsum = weight.sum(1)
    for _ in range(iterations):
        diff = xy[:, None, :] - xy[None, :, :]
        dist = np.maximum(np.linalg.norm(diff, axis=2), 1e-6)
        # Gauss-Seidel-free SMACOF update (weighted, majorizing)
        ratio = weight * target / dist
        np.fill_diagonal(ratio, 0.0)
        new = (weight @ xy + ratio.sum(1)[:, None] * xy - ratio @ xy) / wsum[:, None]
        xy = new
    xy[:, 0] *= aspect ** .5
    xy[:, 1] /= aspect ** .5
    return xy


def _rect_gap(a_lo, a_hi, b_lo, b_hi):
    gap = np.maximum(0.0, np.maximum(b_lo - a_hi, a_lo - b_hi))
    return np.hypot(gap[..., 0], gap[..., 1])


def metrics(centres, sizes, weights, links):
    """The module-level measures of one arrangement (see the module docstring)."""
    half = np.asarray(sizes, float) / 2
    lo, hi = centres - half, centres + half
    pairs = sorted(links)
    a = np.array([p[0] for p in pairs], int)
    b = np.array([p[1] for p in pairs], int)
    cnt = np.array([links[p] for p in pairs], float)
    w = np.array([weights[p] for p in pairs], float)
    gap = _rect_gap(lo[a], hi[a], lo[b], hi[b]) if len(pairs) else np.zeros(0)
    crossings = 0
    for k in range(len(pairs)):
        p, q = centres[a[k]], centres[b[k]]
        d = q - p
        t0, t1 = np.zeros(len(centres)), np.ones(len(centres))
        ok = np.ones(len(centres), bool)
        for ax in (0, 1):
            if abs(d[ax]) < 1e-9:
                ok &= (p[ax] >= lo[:, ax]) & (p[ax] <= hi[:, ax])
            else:
                ta, tb = (lo[:, ax] - p[ax]) / d[ax], (hi[:, ax] - p[ax]) / d[ax]
                t0 = np.maximum(t0, np.minimum(ta, tb))
                t1 = np.minimum(t1, np.maximum(ta, tb))
        ok &= t0 < t1
        ok[a[k]] = ok[b[k]] = False
        crossings += int(ok.any()) * int(cnt[k])
    width, height = (hi.max(0) - lo.min(0))
    covered = float(np.sum((hi - lo)[:, 0] * (hi - lo)[:, 1]))
    return dict(links=int(cnt.sum()), links_over_250=int(cnt[gap > 250].sum()), links_over_1000=int(cnt[gap > 1000].sum()),
                module_crossings=int(crossings), gap_weighted=round(float((w * gap).sum()), 1), gap_median=round(float(np.median(gap)), 1) if len(gap) else 0.0,
                canvas=[round(float(width)), round(float(height))], aspect=round(float(width / max(height, 1.0)), 2),
                compactness=round(covered / max(float(width * height), 1.0), 3))


def choose(rows):
    """rows: [{name, metrics}] -> (index, reason)."""
    ok = [i for i, r in enumerate(rows) if ASPECT_RANGE[0] <= r['metrics']['aspect'] <= ASPECT_RANGE[1]] or list(range(len(rows)))
    best_compact = max(rows[i]['metrics']['compactness'] for i in ok)
    ok = [i for i in ok if rows[i]['metrics']['compactness'] >= COMPACT_TOLERANCE * best_compact]
    pool = list(ok)
    steps = []
    for key in ('links_over_1000', 'links_over_250', 'module_crossings'):
        low = min(rows[i]['metrics'][key] for i in pool)
        keep = [i for i in pool if rows[i]['metrics'][key] <= low * (1 + LONG_TOLERANCE) + 0.5]
        steps.append(f'{key}: best {low}, kept {[rows[i]["name"] for i in keep]}')
        pool = keep
        if len(pool) == 1:
            break
    winner = pool[0]
    discarded = [rows[i]['name'] for i in range(len(rows)) if i not in ok]
    return winner, dict(rule='aspect and compactness screen, then links_over_1000 > links_over_250 > module_crossings (5% tolerance)', discarded=discarded, steps=steps)


def place(items, weights, links, layout_blocks, log=None):
    """Build the candidates, score them, return (centres of the winner, report)."""
    sizes = np.array([i['size'] for i in items], float)
    specs = [('spring', dict()),
             ('stress', dict(init='stress', aspect=1.4)),
             ('stress_square', dict(init='stress', aspect=1.0, compress=.5))]
    rows, arrangements = [], []
    for name, spec in specs:
        t = time.time()
        if spec.get('init') == 'stress':
            init = stress_layout(sizes, weights, aspect=spec['aspect'], compress=spec.get('compress', 1.0))
            centres = layout_blocks(sizes, weights, aspect=spec['aspect'], init=init)
        else:
            centres = layout_blocks(sizes, weights)
        m = metrics(centres, sizes + 90.0, weights, links)
        rows.append(dict(name=name, metrics=m, seconds=round(time.time() - t, 1)))
        arrangements.append(centres)
        if log:
            log(f'  candidate {name}: >250 {m["links_over_250"]}, >1000 {m["links_over_1000"]}, crossings {m["module_crossings"]}, aspect {m["aspect"]}, compactness {m["compactness"]} ({rows[-1]["seconds"]} s)')
    win, why = choose(rows)
    why['chosen'] = rows[win]['name']
    if log:
        log(f'  chosen arrangement: {rows[win]["name"]}')
    return arrangements[win], dict(candidates=rows, selection=why)
