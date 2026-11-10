"""Opt-in P0 physical geometry: one display graph and explicit anchor contracts.

No source chemistry or atlas coordinates are modified. Historical experiment
conditions remain reproducible; `geometry` runs the new controlled ablation.
"""
import copy
import time
import numpy as np
import networkx as nx
from scipy.optimize import minimize
from scipy.spatial.distance import pdist
from core import currency, SEEDS

CONFIG = dict(min_core_edge=55., min_node_spacing=16., max_scale=3000.,
              soft_motion=120., collision_weight=30., shape_weight=.35,
              max_iterations=350, font_size=12)


def module_scale(model, proposal):
    """Uniform scale from target incidence and physical glyph/label footprint."""
    import drawing
    positions = proposal['positions']
    lookup = {r['id']: r for r in model['reactions']}
    lengths = []
    for rid in proposal['reaction_ids']:
        for mid in lookup[rid]['metabolites']:
            a, b = 'r:' + rid, 'm:' + mid
            if a in positions and b in positions:
                v = np.linalg.norm(np.array(positions[a]) - positions[b])
                if v > 1e-7:
                    lengths.append(v)
    arr = np.array(list(positions.values()))
    distances = pdist(arr) if len(arr) > 1 else np.array([1.])
    positive = distances[distances > 1e-7]
    edge_scale = CONFIG['min_core_edge'] / np.median(lengths) if lengths else 100.
    # Each label is measured using the export font. Density is a lower bound,
    # not a claim that all labels have been packed without collisions.
    footprints = sum((drawing.FONT.getlength(k.split(':', 1)[-1]) + 20) * 32
                     for k in positions)
    extent = np.ptp(arr, axis=0)
    density_scale = np.sqrt(footprints / max(float(np.prod(extent)), .01))
    spacing_scale = CONFIG['min_node_spacing'] / np.percentile(positive, 2) if len(positive) else 100.
    return float(min(CONFIG['max_scale'], max(100., edge_scale, density_scale, spacing_scale)))


def scene_graph(model, scene):
    """Exact visible biochemical edges, including local metabolite copies."""
    g = nx.Graph()
    for key in scene['positions']:
        g.add_node(key, role='reaction' if key.startswith('r:') else 'metabolite')
    for r in model['reactions']:
        for mid, coefficient in r['metabolites'].items():
            g.add_edge('r:' + r['id'], scene['participants'][r['id'], mid],
                       coefficient=coefficient)
    return g


def project_spacing(xy, hard, low, high, distance=16., steps=160):
    """Resolve glyph collisions without silently moving hard/bounded anchors."""
    xy = xy.copy()
    for _ in range(steps):
        moved = 0
        for i in range(len(xy)):
            delta = xy[i] - xy[i+1:]
            ds = np.linalg.norm(delta, axis=1)
            for off in np.flatnonzero(ds < distance - 1e-6):
                j = i + 1 + int(off)
                if hard[i] and hard[j]:
                    continue
                direction = xy[i] - xy[j]
                norm = np.linalg.norm(direction)
                if norm < 1e-9:
                    # Stable symmetry breaking, independent of process hash.
                    angle = (i * 17 + j * 31) * 2.399963
                    direction = np.array([np.cos(angle), np.sin(angle)])
                else:
                    direction /= norm
                shift = direction * (distance - norm + .01)
                if not hard[i] and not hard[j]:
                    xy[i] += shift / 2; xy[j] -= shift / 2
                elif not hard[i]: xy[i] += shift
                else: xy[j] -= shift
                xy[:] = np.maximum(low, np.minimum(high, xy))
                moved += 1
        if not moved: break
    return xy


def relax(g, positions, hints, shape_edges=(), iterations=350):
    """Joint positions: hard anchors, bounded soft priors, shape and spacing.

    Optimizer failures and unresolved constraints are reported, never called a
    successful layout merely because coordinates were returned.
    """
    if set(g) != set(positions): raise ValueError('Scene nodes do not match graph')
    if set(hints) - set(positions): raise ValueError('Unknown anchor nodes')
    keys = sorted(positions); ids = {k: i for i, k in enumerate(keys)}
    initial = np.array([positions[k] for k in keys], float)
    hard = np.array([bool(hints.get(k, {}).get('hard', False)) for k in keys])
    weights = np.array([.008 + .04 * hints[k].get('confidence', .5) if k in hints else .002
                        for k in keys])
    low = np.full(initial.shape, -np.inf); high = np.full(initial.shape, np.inf)
    for k, hint in hints.items():
        if k not in ids: continue
        i = ids[k]; reference = np.array(hint['xy'], float)
        radius = 0. if hint.get('hard') else float(hint.get('allowed_motion', CONFIG['soft_motion']))
        if radius < 0 or not np.isfinite(radius): raise ValueError('Invalid allowed_motion')
        initial[i] = reference
        low[i] = reference - radius; high[i] = reference + radius
    ea = np.array([[ids[a], ids[b]] for a, b in sorted(tuple(sorted(e)) for e in g.edges())], dtype=int).reshape(-1, 2)
    se = np.array([[ids[a], ids[b]] for a, b in shape_edges if a in ids and b in ids], dtype=int).reshape(-1, 2)
    original_vectors = initial[se[:, 0]] - initial[se[:, 1]]
    pa, pb = np.triu_indices(len(keys), 1)

    def energy(flat):
        xy = flat.reshape(-1, 2); grad = 2 * weights[:, None] * (xy - initial)
        value = float(np.sum(weights[:, None] * (xy - initial) ** 2))
        if len(ea):
            delta = xy[ea[:, 0]] - xy[ea[:, 1]]
            lengths = np.maximum(np.linalg.norm(delta, axis=1), 1e-8)
            error = lengths - 65.
            value += .015 * np.sum(error ** 2)
            force = .03 * error[:, None] * delta / lengths[:, None]
            np.add.at(grad, ea[:, 0], force); np.add.at(grad, ea[:, 1], -force)
        if len(se):
            error = xy[se[:, 0]] - xy[se[:, 1]] - original_vectors
            value += CONFIG['shape_weight'] * np.sum(error ** 2)
            force = 2 * CONFIG['shape_weight'] * error
            np.add.at(grad, se[:, 0], force); np.add.at(grad, se[:, 1], -force)
        delta = xy[pa] - xy[pb]; lengths = np.maximum(np.linalg.norm(delta, axis=1), 1e-8)
        mask = lengths < CONFIG['min_node_spacing']
        error = lengths[mask] - CONFIG['min_node_spacing']
        value += CONFIG['collision_weight'] * np.sum(error ** 2)
        force = 2 * CONFIG['collision_weight'] * error[:, None] * delta[mask] / lengths[mask, None]
        np.add.at(grad, pa[mask], force); np.add.at(grad, pb[mask], -force)
        return float(value), grad.ravel()

    result = minimize(energy, initial.ravel(), jac=True, method='L-BFGS-B',
                      bounds=list(zip(low.ravel(), high.ravel())),
                      options=dict(maxiter=iterations, ftol=1e-9, maxls=30))
    final = project_spacing(result.x.reshape(-1, 2), hard, low, high,
                            CONFIG['min_node_spacing'])
    assert np.isfinite(final).all()
    assert np.allclose(final[hard], initial[hard], atol=1e-7)
    assert np.all(final >= low - 1e-7) and np.all(final <= high + 1e-7)
    distances = pdist(final)
    error = final[se[:, 0]] - final[se[:, 1]] - original_vectors
    return {k: final[i].tolist() for i, k in enumerate(keys)}, dict(
        optimizer_converged=bool(result.success), optimizer_message=str(result.message),
        iterations=int(result.get('nit',0)), hard_anchors=int(hard.sum()),
        hard_anchor_max_drift=float(np.max(np.linalg.norm(final[hard]-initial[hard], axis=1))) if hard.any() else 0.,
        remaining_spacing_violations=int(np.sum(distances < CONFIG['min_node_spacing'] - 1e-5)),
        node_overlap_pairs=int(np.sum(distances < 10)),
        shape_vector_rmse=float(np.sqrt(np.mean(error ** 2))) if len(error) else 0.,
        motion_bound='axis-aligned box in display pixels',
        hint_max_drift=float(max((np.linalg.norm(final[ids[k]]-h['xy']) for k,h in hints.items() if k in ids), default=0.)))


def complete_core(g, hints, seed=SEEDS[0], iterations=350):
    """Typed hard/soft contract; old coordinate-only callers keep old behavior."""
    import drawing
    hints = {k: (v if isinstance(v, dict) else dict(xy=v, hard=True)) for k,v in hints.items()}
    positions = drawing.initial(g, {k:v['xy'] for k,v in hints.items()}, seed)
    return relax(g, positions, hints, iterations=iterations)[0]


def complete_scene(model, scene, blocks, hints=None, iterations=350):
    import drawing
    scene = copy.deepcopy(scene); hints = copy.deepcopy(hints or {})
    g = scene_graph(model, scene); shape_edges = []
    for b in blocks:
        node_keys = set()
        for rid in b['reaction_ids']:
            rk = 'r:'+rid
            if rk not in scene['positions']: continue
            node_keys.add(rk)
            for r in model['reactions']:
                if r['id'] != rid: continue
                for mid in r['metabolites']:
                    if 'm:'+mid not in b['positions']: continue
                    mk = scene['participants'][rid,mid]
                    if not mk.startswith('side:'): node_keys.add(mk)
        for key in node_keys:
            hints.setdefault(key, dict(xy=scene['positions'][key], hard=False,
                                      confidence=b.get('anchor_confidence', .5),
                                      allowed_motion=CONFIG['soft_motion'], module_id=b['id'],
                                      origin='retrieved_template'))
        shape_edges.extend(sorted(tuple(sorted((a,b))) for a,b in g.subgraph(node_keys).edges()))
    final, audit = relax(g, scene['positions'], hints, shape_edges, iterations)
    scene['positions'] = final
    scene['anchors'] = hints
    return scene, audit
