"""Peripheral attachment (V7): small exchange + transport modules go to the metabolic region their intracellular metabolite belongs to.

A peripheral module is the end of the network: an EXCHANGE reaction (one participant, the extracellular metabolite x_e) and the SIMPLE transport reaction(s) that
move x_e to the same chemical in another compartment (x_c) with nothing but protons / water as co-substrate.  Such a module is tidy on its own, but the
layout used to put it wherever the module arrangement left room, so a transporter whose x_c is used by other reactions ended up far from them (an island
with 0 crossings and 100 % curves that still reads wrongly).

What is decided from the GEM itself (stoichiometry and compartments; never from metabolite names):
  identification  EX reaction with a single non-currency participant x_e + non-expert transport reactions T whose non-currency participants are exactly x_e and one
                  x_c of the same chemical (id without the compartment suffix) and whose currency participants are only H+ / H2O; complex or ATP-coupled transporters
                  (NIabc ...) are not peripheral modules: the normal skeleton / local layout handles them
  anchors         the reactions OUTSIDE the module that use x_c as a core (non-currency) participant, ranked by trust: an expert structure (3) > a reaction the
                  V6 placement or the backbone placed (2) > any other placed reaction (1); the best rank decides
  no anchor       the module is a true independent component: it stays where it is and no connection is invented
  move            the whole module (EX, T, their shared x_e glyph and the currency glyphs beside them) is translated rigidly so that T sits PERIPH_DISTANCE px from
                  the x_c glyph, on the side away from the anchors; the move is kept only if T's connection to x_c gets clearly shorter and the free-spot search finds
                  room (nothing else moves; expert glyphs are never touched)
Counts for the quality panel: attachable (an anchor exists), attached (moved next to it), not attached (anchor exists, no move was possible), independent.
"""
from collections import defaultdict

import numpy as np
from scipy.spatial import cKDTree

from core import chem, comp

PERIPH_DISTANCE = 120.0
MIN_GAIN = 100.0
CLEAR = 80.0
SIMPLE_CURRENCY = {'h', 'h2o'}
CLOSE = 260.0   # a module whose x_c glyph is already this close to its anchors counts as attached where it is


def find_modules(model, cur, expert, blocks=(), simple=True):
    """Exchange + simple transport modules.  A reaction that an expert (or topology-matched) structure placed counts only if that whole structure is the module
    itself (a rigid two-reaction block can be moved as a unit); a block that holds more than the module is left alone."""
    rx = {r['id']: r for r in model['reactions']}
    block_of = {rid: tuple(sorted(b['reaction_ids'])) for b in blocks for rid in b['reaction_ids']}
    users = defaultdict(set)
    for r in model['reactions']:
        for mid in r['metabolites']:
            users[mid].add(r['id'])
    modules = []
    for r in model['reactions']:
        core = [m for m in r['metabolites'] if m not in cur]
        if len(r['metabolites']) != 1 or len(core) != 1:
            continue
        xe = core[0]
        transports = []
        for tid in sorted(users[xe] - {r['id']}):
            t = rx[tid]
            tcore = [m for m in t['metabolites'] if m not in cur]
            tcur = [m for m in t['metabolites'] if m in cur]
            if xe not in tcore or (simple and len(tcore) != 2):
                continue
            if simple:
                xc = next(m for m in tcore if m != xe)
                if chem(xc) != chem(xe) or comp(xc) == comp(xe):
                    continue
            else:
                # coupled transport: the metabolite crosses with co-transported ions (Na+, K+ ...) that are not in the currency list; the pair (x_e, x_c) is still the
                # same chemical in two compartments, the other participants are allowed
                same = [m for m in tcore if m != xe and chem(m) == chem(xe) and comp(m) != comp(xe)]
                if len(same) != 1 or len(tcore) > 6:
                    continue
                xc = same[0]
            if simple and any(chem(m) not in SIMPLE_CURRENCY for m in tcur):
                continue     # simple=False (expert sub-module mode 'branch'): symporters / antiporters / ATP-driven transporters with currency co-substrates are leaf branches too
            transports.append((tid, xc))
        if not 1 <= len(transports) <= 3:
            continue
        members = {r['id']} | {t for t, _ in transports}
        ok = True
        for m in members:
            if m in block_of and not set(block_of[m]) <= members:
                ok = False  # the structure that holds this reaction is bigger than the module
        if ok:
            modules.append(dict(exchange=r['id'], xe=xe, transports=transports))
    return modules


def attach(model, scene, blocks, cur, apply=True):
    """-> (new scene, report).  Moves peripheral modules next to their anchors (rigid), reports the four counts."""
    scene = dict(scene, positions={k: list(v) for k, v in scene['positions'].items()})
    pos, parts = scene['positions'], scene['participants']
    expert = {rid for b in blocks for rid in b['reaction_ids']}
    placed = set(scene.get('completion_context', {}).get('tiers', {}))
    rx = {r['id']: r for r in model['reactions']}
    users = defaultdict(set)
    for r in model['reactions']:
        for mid in r['metabolites']:
            if mid not in cur:
                users[mid].add(r['id'])
    holders = defaultdict(set)
    for (rid, mid), key in parts.items():
        holders[key].add(rid)
    modules = find_modules(model, cur, expert, blocks)
    report = dict(modules=len(modules), attachable=0, attached=0, not_attached=0, independent=0, proposed=0, details=[])
    members_all = {r for m in modules for r in [m['exchange']] + [t for t, _ in m['transports']]}
    movers = {k for k in pos}
    for mod in modules:
        members = [mod['exchange']] + [t for t, _ in mod['transports']]
        inside = set(members)
        # anchors of every transporter's intracellular end
        best_rank, anchors, xcs = 0, [], []
        for tid, xc in mod['transports']:
            xcs.append(xc)
            for u in users[xc] - inside:
                if u not in rx or 'r:' + u not in pos:
                    continue
                if u in members_all:
                    continue
                rank = 3 if u in expert else (2 if u in placed else 1)
                if rank > best_rank:
                    best_rank, anchors = rank, []
                if rank == best_rank:
                    anchors.append((u, xc))
        info = dict(exchange=mod['exchange'], transports=[t for t, _ in mod['transports']], xe=mod['xe'], xc=xcs, anchors=sorted({u for u, _ in anchors})[:6], rank=best_rank)
        if not anchors:
            report['independent'] += 1
            info['state'] = 'independent'
            report['details'].append(info)
            continue
        report['attachable'] += 1
        # the x_c glyph the transporter attaches to (a shared node, or an expert copy)
        tid0, xc0 = mod['transports'][0]
        own_key = parts[tid0, xc0]                       # the module's own x_c glyph (a shared node, or its structure's copy)
        own_xy = np.array(pos[own_key], float)
        anchor_glyphs = np.array([pos[parts[u, xc]] for u, xc in anchors], float)   # where the anchors draw that metabolite
        centre = anchor_glyphs.mean(axis=0)
        near_glyph = anchor_glyphs[int(np.argmin(np.hypot(*(anchor_glyphs - own_xy).T)))]
        away = own_xy - centre
        n = float(np.hypot(*away))
        away = away / n if n > 1e-6 else np.array([1.0, 0.0])
        # everything of the module moves together: its reactions, the x_e glyph, glyphs nobody else uses, and the glyphs of a structure that IS the module
        own_blocks = [bk for bk in blocks if set(bk['reaction_ids']) <= inside]
        keys = {'r:' + r for r in members}
        for bk in own_blocks:
            keys.update(k for k in pos if k.startswith('copy:' + bk['id'] + ':'))
        for r in members:
            for mid in rx[r]['metabolites']:
                k = parts[r, mid]
                if holders[k] <= inside and (not k.startswith('copy:') or any(k.startswith('copy:' + bk['id'] + ':') for bk in own_blocks)):
                    keys.add(k)
        keys = {k for k in keys if k in pos}
        t_pos = np.array(pos['r:' + tid0], float)
        g = near_glyph
        ref = own_xy
        old_len = float(np.hypot(*(g - ref)))
        target = g + away * PERIPH_DISTANCE
        delta = target - ref
        local = np.array([pos[k] for k in keys], float)
        others = np.array([v for k, v in pos.items() if k not in keys and not k.startswith('side:')], float)
        tree = cKDTree(others) if len(others) else None
        found = None
        for i in range(0, 90):
            r = 0.0 if i == 0 else 38.0 * np.sqrt(i)
            a = i * 2.399963
            d = delta + r * np.array([np.cos(a), np.sin(a)])
            if tree is None or float(tree.query(local + d)[0].min()) >= CLEAR * 0.6:
                found = d
                break
        new_len = float(np.hypot(*(g - (ref + found)))) if found is not None else old_len
        info['connection_px'] = [round(old_len, 1), round(new_len if found is not None else old_len, 1)]
        if old_len <= CLOSE:
            report['attached'] += 1
            info['state'] = 'attached'
            info['in_place'] = True
        elif found is None or old_len - new_len < MIN_GAIN:
            report['not_attached'] += 1
            info['state'] = 'not_attached'
            info['reason'] = 'no free place within reach' if found is None else 'gain below the minimum'
            info['x'], info['y'] = round(float(t_pos[0]), 1), round(float(t_pos[1]), 1)
        else:
            if apply:
                for k in keys:
                    pos[k] = (np.array(pos[k], float) + found).tolist()
                report['attached'] += 1
                info['state'] = 'attached'
            else:                      # proposed only: the quality gate decides on the final curves
                info['state'] = 'proposed'
            info['moved_px'] = round(float(np.hypot(*found)), 1)
            info['delta'] = [float(found[0]), float(found[1])]
            info['keys'] = sorted(keys)
            info['reactions'] = list(members)
        report['details'].append(info)
    report['proposed'] = sum(d['state'] == 'proposed' for d in report['details'])
    scene['peripheral'] = report
    return scene, report
