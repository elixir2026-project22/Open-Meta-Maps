"""Completion of the layout for reactions that no expert structure covers, without crowding.

drawing.complete runs fCoSE with a fixed ideal edge length (70) and one repulsion value for everything, and the result is not checked
against the size of the things that have to be drawn.  Uncovered reactions end up on top of each other (nearest-neighbour distance of
generated reactions ~31 px against ~70 px for expert-derived ones in e_coli_core), so reaction labels, side metabolites and routes overlap.

`complete_spaced` keeps the classic result as the starting point and then pushes only the FREE nodes (nodes of the core graph that are not
anchored by an expert structure) apart until every free node keeps a minimum distance from its neighbours, anchored ones included.
Anchored nodes never move, so the expert layout is untouched.  Deterministic: pairs are processed in sorted order, no random numbers.

`drawing.py` is not modified (its hash is recorded in saved results); the tool chooses `classic` or `spaced` with `COMPLETION`.
"""
import math

import numpy as np
from scipy.spatial import cKDTree

import drawing

# minimum centre distance (px, physical scale of the tool map; reaction labels are ~60-80 px wide at 12 px Menlo)
MIN_RX_RX = 130.0     # free reaction - any reaction
MIN_MET_MET = 85.0    # free metabolite - any metabolite
MIN_RX_MET = 65.0     # free reaction - non-adjacent metabolite
ITERATIONS = 300


def _need(a, b, spacing=None):
    rr, mm, rm = spacing or (MIN_RX_RX, MIN_MET_MET, MIN_RX_MET)
    ra, rb = a.startswith('r:'), b.startswith('r:')
    if ra and rb:
        return rr
    if not ra and not rb:
        return mm
    return rm


def spread(xy, g, anchors, iterations=ITERATIONS, spacing=None):
    """Push free nodes apart in place; returns how many nodes moved and the number of unresolved pairs."""
    keys = sorted(xy)
    index = {k: i for i, k in enumerate(keys)}
    pos = np.array([xy[k] for k in keys], dtype=float)
    free = np.array([k not in anchors for k in keys])
    if not free.any():
        return 0, 0
    radius = max(spacing or (MIN_RX_RX, MIN_MET_MET, MIN_RX_MET))
    start = pos.copy()
    unresolved = 0
    for it in range(iterations):
        tree = cKDTree(pos)
        pairs = sorted(tree.query_pairs(radius))
        push = np.zeros_like(pos)
        unresolved = 0
        for i, j in pairs:
            if not (free[i] or free[j]):
                continue
            a, b = keys[i], keys[j]
            if g.has_edge(a, b):
                continue  # connected nodes are held together by the edge length, not pushed apart
            need = _need(a, b, spacing)
            d = pos[j] - pos[i]
            dist = float(np.hypot(*d))
            if dist >= need:
                continue
            unresolved += 1
            if dist < 1e-6:
                angle = (index[a] * 2.399963) % (2 * math.pi)  # deterministic jitter for coincident nodes
                d, dist = np.array([math.cos(angle), math.sin(angle)]), 1.0
            step = d / dist * (need - dist)
            if free[i] and free[j]:
                push[i] -= step / 2
                push[j] += step / 2
            elif free[i]:
                push[i] -= step
            else:
                push[j] += step
        if not unresolved:
            break
        pos += push * 0.6
    moved = int((np.hypot(*(pos - start).T) > 1.0).sum())
    for k in keys:
        if k not in anchors:
            xy[k] = pos[index[k]].tolist()
    return moved, unresolved


def complete_spaced(g, anchors, seed=None, iterations=350):
    kwargs = {} if seed is None else dict(seed=seed)
    xy = drawing.complete(g, anchors, iterations=iterations, **kwargs)
    xy = {k: list(v) for k, v in xy.items()}
    moved, unresolved = spread(xy, g, anchors)
    complete_spaced.last = dict(moved=moved, unresolved_pairs=unresolved)
    return xy
