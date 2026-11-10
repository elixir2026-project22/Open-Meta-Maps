"""Reaction orientation: substrates on one side of the reaction symbol, products on the other.

The project's layout rules (hard rule: "placement of reaction participants should be consistent in relation to the reaction/process node")
and the expert maps agree: in the Escher e_coli_core map 89 % of the reactions have their substrates and products on opposite sides of the
reaction point (median angle 156 degrees).  The tool's maps had 28-30 %, for two reasons:

  * `drawing.display_scene` puts every currency copy of a reaction directly below it, whatever its role;
  * reactions outside the expert structures are placed by force layouts that know nothing about substrate and product sides.

This module corrects both after the scene exists (the hashed drawing files are not modified):

  free_reactions(model, scene, blocks)   moves reactions that no expert structure placed between their substrates and their products
  sides(model, scene, blocks)            puts currency substrates on the substrate side of their reaction and currency products on the
                                         product side, spread across the direction of the reaction
  angles(...) / summary(...)             the metric: angle between the mean substrate direction and the mean product direction

Expert structures are never moved: their reaction points and their glyph copies (`copy:` keys) stay where the expert drew them.
The orientation uses the stoichiometric sign of the participants (substrate = negative coefficient), not the flux direction.
"""
import math

import numpy as np

from core import currency


def _expert_reactions(blocks):
    return {rid for b in blocks for rid in b['reaction_ids']}


def _split(reaction, participants, positions, cur, only_core=True):
    subs, prods = [], []
    for mid, coefficient in reaction['metabolites'].items():
        if only_core and mid in cur:
            continue
        key = participants[reaction['id'], mid]
        (subs if coefficient < 0 else prods).append(np.array(positions[key], float))
    return subs, prods


def free_reactions(model, scene, blocks, blend=1.0, max_shift=600.0):
    """Move every reaction outside the expert structures to the point between its (non-currency) substrates and products."""
    cur = currency(model)
    expert = _expert_reactions(blocks)
    scene = dict(scene, positions={k: list(v) for k, v in scene['positions'].items()})
    positions, moved = scene['positions'], 0
    for r in model['reactions']:
        if r['id'] in expert:
            continue
        subs, prods = _split(r, scene['participants'], positions, cur)
        if not subs or not prods:
            continue
        target = (np.mean(subs, axis=0) + np.mean(prods, axis=0)) / 2
        here = np.array(positions['r:' + r['id']], float)
        shift = target - here
        norm = float(np.hypot(*shift))
        if norm > max_shift:
            shift *= max_shift / norm
        if norm > 1.0:
            positions['r:' + r['id']] = (here + blend * shift).tolist()
            moved += 1
    scene['orientation_free_moved'] = moved
    return scene


def sides(model, scene, blocks, distance=50.0, spread=30.0):
    """Currency substrates on the substrate side, currency products on the product side (`side:` glyphs only; they are not expert geometry)."""
    cur = currency(model)
    scene = dict(scene, positions={k: list(v) for k, v in scene['positions'].items()})
    positions = scene['positions']
    moved = 0
    for r in model['reactions']:
        here = np.array(positions['r:' + r['id']], float)
        subs, prods = _split(r, scene['participants'], positions, cur)
        if subs and prods:
            u = np.mean(prods, axis=0) - np.mean(subs, axis=0)
        elif subs:
            u = here - np.mean(subs, axis=0)
        elif prods:
            u = np.mean(prods, axis=0) - here
        else:
            u = np.array([1.0, 0.0])
        n = float(np.hypot(*u))
        u = np.array([1.0, 0.0]) if n < 1e-6 else u / n
        perp = np.array([-u[1], u[0]])
        groups = {-1: [], 1: []}
        for mid, coefficient in r['metabolites'].items():
            if mid in cur:
                key = scene['participants'][r['id'], mid]
                if key.startswith('side:'):
                    groups[-1 if coefficient < 0 else 1].append(key)
        for sign, keys in groups.items():
            for i, key in enumerate(keys):
                offset = (i - (len(keys) - 1) / 2) * spread
                positions[key] = (here + sign * u * distance + perp * offset).tolist()
                moved += 1
    scene['orientation_sides_moved'] = moved
    return scene


def separate_free_reactions(model, scene, blocks, reaction_gap=95.0, other_gap=42.0, passes=60):
    """After moving reactions between their substrates and products: keep reactions that no structure placed from landing on each other
    (reaction labels are 60-80 px wide) or on glyphs they are not connected to.  Only these reactions move; structures and glyphs stay."""
    from scipy.spatial import cKDTree
    expert = _expert_reactions(blocks)
    scene = dict(scene, positions={k: list(v) for k, v in scene['positions'].items()})
    positions = scene['positions']
    keys = sorted(positions)
    index = {k: i for i, k in enumerate(keys)}
    xy = np.array([positions[k] for k in keys], float)
    movable = np.zeros(len(keys), bool)
    for r in model['reactions']:
        if r['id'] not in expert:
            movable[index['r:' + r['id']]] = True
    is_reaction = np.array([k.startswith('r:') for k in keys])
    linked = {}
    for r in model['reactions']:
        i = index['r:' + r['id']]
        linked[i] = {index[scene['participants'][r['id'], m]] for m in r['metabolites']}
    moved = set()
    for _ in range(passes):
        pairs = cKDTree(xy).query_pairs(reaction_gap, output_type='ndarray')
        push = np.zeros_like(xy)
        any_violation = False
        for i, j in pairs:
            if not (movable[i] or movable[j]):
                continue
            need = reaction_gap if (is_reaction[i] and is_reaction[j]) else other_gap
            if is_reaction[i] != is_reaction[j]:
                a, b = (i, j) if is_reaction[i] else (j, i)
                if b in linked.get(a, ()):
                    continue
            d = xy[j] - xy[i]
            dist = float(np.hypot(*d))
            if dist >= need:
                continue
            if dist < 1e-6:
                angle = (i * 17 + j * 31) * 2.399963
                d, dist = np.array([np.cos(angle), np.sin(angle)]), 1.0
            step = d / dist * (need - dist)
            any_violation = True
            if movable[i] and movable[j]:
                push[i] -= step / 2
                push[j] += step / 2
            elif movable[i]:
                push[i] -= step
            else:
                push[j] += step
        if not any_violation:
            break
        xy += push * .7
        moved |= set(np.flatnonzero(np.hypot(*push.T) > 0))
    for k in keys:
        positions[k] = xy[index[k]].tolist()
    scene['orientation_separated'] = len(moved)
    return scene


def apply(model, scene, blocks, mode, stage):
    """`mode` off / sides / full; stage 'before_relax' (free reactions) or 'after_relax' (currency sides, final)."""
    if mode == 'full' and stage == 'before_relax':
        return free_reactions(model, scene, blocks)
    if mode in ('sides', 'full') and stage == 'after_relax':
        if mode == 'full':
            scene = separate_free_reactions(model, scene, blocks)
        return sides(model, scene, blocks)
    return scene


# ------------------------------------------------------------------------------------------------ metric
def angle_of(subs, prods):
    if not subs or not prods:
        return None
    a, b = np.mean(subs, axis=0), np.mean(prods, axis=0)
    na, nb = float(np.hypot(*a)), float(np.hypot(*b))
    if na < 1e-9 or nb < 1e-9:
        return None
    return math.degrees(math.acos(max(-1.0, min(1.0, float(a @ b) / (na * nb)))))


def angles_from_map(map_json, skip=None):
    """{reaction id: angle in degrees between the mean substrate and product directions}, from an Escher map (all drawn glyphs, as in the
    comparison with the expert map)."""
    body = map_json[1]
    nodes = body['nodes']
    out = {}
    for rx in body['reactions'].values():
        mid = None
        for seg in rx['segments'].values():
            for k in (seg['from_node_id'], seg['to_node_id']):
                if nodes[k]['node_type'] == 'midmarker':
                    mid = nodes[k]
        if mid is None:
            continue
        coefficient = {m['bigg_id']: m['coefficient'] for m in rx['metabolites']}
        subs, prods = [], []
        for seg in rx['segments'].values():
            for k in (seg['from_node_id'], seg['to_node_id']):
                n = nodes[k]
                if n['node_type'] == 'metabolite' and not (skip and skip(n['bigg_id'])):
                    v = np.array([n['x'] - mid['x'], n['y'] - mid['y']])
                    (subs if coefficient.get(n['bigg_id'], 0) < 0 else prods).append(v)
        angle = angle_of(subs, prods)
        if angle is not None:
            out[rx['bigg_id']] = angle
    return out


def summary(map_json, kinds=None, threshold=120.0):
    """Share of reactions with substrates and products on opposite sides (angle > threshold).

    `all_participants` counts every drawn glyph (the measure used to compare with the expert map); `core_participants` ignores currency
    glyphs, so it shows whether the core metabolites themselves are laid out substrate-to-product.  `kinds` (reaction id -> provenance kind)
    splits the reactions an expert structure placed from the ones the algorithm placed."""
    import re
    from core import CURRENCY
    is_cur = lambda m: re.sub(r'_[a-z][0-9]?$', '', m) in CURRENCY
    out = {}
    for label, skip in (('all_participants', None), ('core_participants', is_cur)):
        angles = angles_from_map(map_json, skip)
        if not angles:
            out[label] = dict(reactions=0)
            continue
        vals = np.array(list(angles.values()))
        row = dict(reactions=len(vals), opposite_sides_share=round(float((vals > threshold).mean()), 3), angle_median=round(float(np.median(vals))))
        if kinds:
            expert = np.array([a for r, a in angles.items() if kinds.get(r) in ('expert_derived', 'topology_layout_suggestion', 'functional_layout_suggestion')])
            auto = np.array([a for r, a in angles.items() if kinds.get(r) in ('generated', 'bridge')])
            if len(expert):
                row['expert_structure_reactions'] = dict(n=len(expert), opposite_sides_share=round(float((expert > threshold).mean()), 3), angle_median=round(float(np.median(expert))))
            if len(auto):
                row['algorithm_placed_reactions'] = dict(n=len(auto), opposite_sides_share=round(float((auto > threshold).mean()), 3), angle_median=round(float(np.median(auto))))
        out[label] = row
    out['threshold_deg'] = threshold
    out['opposite_sides_share'] = out['all_participants'].get('opposite_sides_share')  # headline number
    out['angle_median'] = out['all_participants'].get('angle_median')
    return out
