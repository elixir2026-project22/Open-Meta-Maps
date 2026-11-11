"""Low-confidence module decomposition (V7.3).

A structure retrieved by topology alone (`match_level == 'geometry'`: no chemical identity was confirmed) is a weak reason to move two reactions as one rigid unit.
When its reactions have very different real neighbourhoods in the target GEM, the rigid unit cannot be near the real partners of both (e_coli_core: BIOMASS + TALA, a
two-reaction topology match; TALA belongs next to the pentose-phosphate reactions, BIOMASS touches many regions, so the unit sat far from both).

What this does (after `select`, before the skeleton; evidence, not layout, is what changes):
  trigger    the proposal is topology-only AND two of its reactions have different external neighbourhoods: the Jaccard overlap of the sets of reactions OUTSIDE the proposal that
             share a non-currency metabolite with them is below JACCARD_MAX.  The number of participants is NOT a trigger by itself (a hub reaction shows up through the overlap).
  split      the proposal is replaced by one proposal per reaction, each keeping that reaction's own source geometry (its position and the glyphs it uses, as in
             `retrieval.clip_proposals`); the evidence level is unchanged (`match_level` stays topology-only), the chemistry of the target GEM is untouched.
  everything else (exact / adapted structures, confirmed chemistry) is left exactly as selected.
The decision is judged on the finished map like every other change (A/B switch OPM_LOWCONF=on|off).
"""


import copy
from collections import defaultdict

from core import uid

JACCARD_MAX = 0.35
DEGREE_RATIO = 4.0   # a member whose outside neighbourhood is this many times larger than another member's is a hub riding on a small structure


def external_neighbours(model, cur, rids, members):
    users = defaultdict(set)
    for r in model['reactions']:
        for m in r['metabolites']:
            if m not in cur:
                users[m].add(r['id'])
    rx = {r['id']: r for r in model['reactions']}
    out = {}
    for rid in rids:
        s = set()
        for m in rx[rid]['metabolites']:
            if m not in cur:
                s |= users[m]
        out[rid] = s - set(members)
    return out


def split(model, cur, selected):
    """-> (new selected list, report)."""
    import retrieval as me  # clip_proposals lives in retrieval (production)
    g = me.model_graph(model)
    rids_all = {r['id'] for r in model['reactions']}
    report = dict(checked=0, split=[], kept=0)
    out = []
    for p in selected:
        rids = sorted(set(p['reaction_ids']) & rids_all)
        if p.get('match_level') != 'geometry' or len(rids) < 2:
            out.append(p)
            continue
        report['checked'] += 1
        nb = external_neighbours(model, cur, rids, rids)
        worst, incompatible = 1.0, False
        rx_m = {r['id']: {m for m in r['metabolites'] if m not in cur} for r in model['reactions'] if r['id'] in rids}
        for i in range(len(rids)):
            for j in range(i + 1, len(rids)):
                a, b = nb[rids[i]], nb[rids[j]]
                union = a | b
                jac = (len(a & b) / len(union)) if union else 1.0
                worst = min(worst, jac)
                ratio = max(len(a), len(b)) / max(min(len(a), len(b)), 1)
                no_link = False
                # different neighbourhoods are normal inside a real pathway chain (each member touches its own neighbours); a pair is incompatible only when
                # one of them is a hub riding on a small structure (degree ratio) or the two share no non-currency metabolite at all
                if jac < JACCARD_MAX and ratio >= DEGREE_RATIO:
                    incompatible = True
        # the structure must also be connected inside: reactions of one pathway chain are linked through shared non-currency metabolites (not every pair, but as a graph)
        comp = {r: r for r in rids}

        def root(x):
            while comp[x] != x:
                comp[x] = comp[comp[x]]
                x = comp[x]
            return x
        for i in range(len(rids)):
            for j in range(i + 1, len(rids)):
                if rx_m[rids[i]] & rx_m[rids[j]]:
                    comp[root(rids[i])] = root(rids[j])
        if len({root(r) for r in rids}) > 1:
            incompatible = True
        if not incompatible:
            out.append(p)
            report['kept'] += 1
            continue
        for r in rids:
            keys = {'r:' + r} | set(g['r:' + r])
            q = copy.deepcopy(p)
            q['reaction_ids'] = [r]
            q['positions'] = {k: v for k, v in p['positions'].items() if k in keys}
            q['node_mapping'] = {k: v for k, v in p['node_mapping'].items() if v in q['positions']}
            q['reaction_mapping'] = {k: v for k, v in p['reaction_mapping'].items() if v == r}
            q['mapping_evidence'] = [e for e in p.get('mapping_evidence', []) if e.get('target_reaction') in (r, None)]
            q['source_glyph_mapping'] = {k: v for k, v in p.get('source_glyph_mapping', {}).items() if 'm:' + v in q['positions']}
            q['decomposed_from'] = dict(proposal=p['id'], reactions=rids, jaccard=round(worst, 3))
            q['id'] = uid([p['id'], 'lowconf-split', r])
            out.append(q)
        report['split'].append(dict(proposal=p['id'][:12], reactions=rids, jaccard=round(worst, 3), source=p.get('source')))
    return out, report
