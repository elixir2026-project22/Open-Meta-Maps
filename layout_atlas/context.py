"""Reaction placement context: reactions whose participants are all "currency" metabolites no longer float alone.

`core.currency` marks O2, CO2, NH4, Pi, H2O, H+, quinones and the nucleotide / coenzyme cofactors as currency, and the Composer draws a currency
metabolite as a private glyph next to every reaction.  A reaction with NO other participant (in e_coli_core: ATPM, ADK1, ATPS4r, NADH16, CYTBD,
THD2, NADTRHD, O2t, EX_o2_e, CO2t, EX_co2_e, NH4t, EX_nh4_e, PIt2r, EX_pi_e, H2Ot, EX_h2o_e, EX_h_e) therefore has nothing that connects it to anything:
every one of them became an island.  The expert Escher e_coli_core map does the opposite for exactly these reactions: O2t and EX_o2_e share one o2_e
node, CYTBD and O2t share o2_c, NADH16 and CYTBD share the quinone nodes, ATPM / ADK1 / ATPS4r share atp_c and adp_c, THD2 and NADTRHD share the
NAD(P)(H) nodes, and the groups sit next to each other.

What this module does (outside the hashed drawing files; expert structures are not touched):

  groups      reactions without any non-currency participant and outside the expert structures are grouped when they share a metabolite that
              at most MAX_SHARERS of them use; for each such metabolite ONE visual node is shared by the group (like the expert map).  Protons
              and other metabolites that too many of them use keep their private glyphs.
  regions     groups that share an extracellular ultra-common metabolite (H+ outside) are held together in one region, so the oxidative-
              phosphorylation reactions, ATP synthase, transhydrogenase and the proton exchange end up next to each other.
  placement   every region is laid out by itself (compact, spaced) and put where its related reactions are: the weighted mean of the reactions
              outside the group that use the same metabolites (a metabolite used by few reactions counts more).  A region whose metabolites are
              used by too many reactions to point anywhere stays where it was.

GEM content is unchanged: every reaction keeps all its participants; only which visual node a participant attaches to changes.
`isolated_reactions(model, scene, cur)` counts reactions that share no visual node with any other reaction (before / after).
"""
from collections import defaultdict, Counter

import numpy as np

from core import chem

MAX_SHARERS = 4      # a metabolite shared by more than this many unplaced reactions is not given a shared node (H+, H2O in the cell, ...)
PARTNER_MIN = 0.6    # total partner weight needed to move a region next to its partner reactions
SPACING = (110.0, 75.0, 60.0)
LAYOUT_STYLE_BACKBONE = False  # set by tool_pipeline when --layout-style expert: regions are laid out as layered pathways (pathway_backbone) instead of a spring embedding


def _extracellular(model_met, mid):
    return model_met.get(mid, {}).get('compartment') in ('e', 'extracellular') or mid.endswith('_e')


def find(model, cur, covered):
    """-> (regions, stats).  A region is dict(reactions, shared {mid: [rids]}, links {mid: [rids]})."""
    met = {m['id']: m for m in model['metabolites']}
    zero = [r for r in model['reactions'] if r['id'] not in covered and r['metabolites'] and all(m in cur for m in r['metabolites'])]
    ids = {r['id'] for r in zero}
    users = defaultdict(set)
    for r in zero:
        for mid in r['metabolites']:
            if chem(mid) in ('h', 'h2o', 'pi') and not _extracellular(met, mid):
                continue  # protons, water and phosphate inside the cell: every reaction keeps its own glyph (the expert map does not share them either)
            users[mid].add(r['id'])
    shared = {mid: sorted(u) for mid, u in users.items() if 2 <= len(u) <= MAX_SHARERS}
    links = {mid: sorted(u) for mid, u in users.items() if len(u) > MAX_SHARERS and _extracellular(met, mid) and chem(mid) in ('h', 'h2o')}
    parent = {rid: rid for rid in ids}

    def root(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for table in (shared, links):
        for mid, rids in table.items():
            for other in rids[1:]:
                parent[root(other)] = root(rids[0])
    groups = defaultdict(list)
    for rid in sorted(ids):
        groups[root(rid)].append(rid)
    regions = []
    for members in sorted(groups.values(), key=lambda m: m[0]):
        ms = set(members)
        regions.append(dict(reactions=members, shared={m: [r for r in rs if r in ms] for m, rs in shared.items() if set(rs) & ms},
                            links={m: [r for r in rs if r in ms] for m, rs in links.items() if set(rs) & ms}))
    stats = dict(reactions_without_core_participants=len(zero), regions=len(regions), singleton_regions=sum(len(r['reactions']) == 1 for r in regions),
                 shared_nodes=sum(len(r['shared']) for r in regions))
    return regions, stats


def place(model, core, cur, regions, covered):
    """Lay each region out and move it to its partners.  Returns (new core, {(rid, mid): node key}, report)."""
    import networkx as nx
    import completion
    core = {k: list(v) for k, v in core.items()}
    reactions = {r['id']: r for r in model['reactions']}
    degree = Counter(mid for r in model['reactions'] for mid in r['metabolites'])
    users = defaultdict(list)
    for r in model['reactions']:
        for mid in r['metabolites']:
            users[mid].append(r['id'])
    assign, moved, pinned = {}, 0, 0
    from scipy.spatial import cKDTree
    in_regions = {rid for region in regions for rid in region['reactions'] if len(region['reactions']) > 1}
    fixed = [v for k, v in core.items() if not (k.startswith('r:') and k[2:] in in_regions)]
    for k, region in enumerate(regions):
        members = region['reactions']
        if len(members) < 2:
            continue
        G = nx.Graph()
        for rid in members:
            G.add_node('r:' + rid)
        for mid, rids in region['shared'].items():
            G.add_node('s:' + mid)
            for rid in rids:
                G.add_edge('r:' + rid, 's:' + mid)
        if LAYOUT_STYLE_BACKBONE:
            import pathway_backbone
            pos = pathway_backbone.layered(reactions, region)
            xy = {n: [float(p[0]), float(p[1])] for n, p in pos.items()}
        else:
            for mid, rids in region['links'].items():
                G.add_node('v:' + mid)
                for rid in rids:
                    G.add_edge('r:' + rid, 'v:' + mid)
            scale = max(120.0, 70.0 * np.sqrt(len(G)))
            pos = nx.spring_layout(G, seed=7, iterations=80, scale=scale)
            xy = {n: [float(p[0]), float(p[1])] for n, p in pos.items()}
            completion.spread(xy, G, {}, spacing=SPACING)
        arr = np.array(list(xy.values()))
        centre_local = arr.mean(0)
        current = np.mean([core['r:' + rid] for rid in members if 'r:' + rid in core], axis=0)
        # where do the related reactions sit?
        weight, acc = 0.0, np.zeros(2)
        for mid in {m for rid in members for m in reactions[rid]['metabolites']}:
            for other in users[mid]:
                if other not in members and 'r:' + other in core:
                    w = 1.0 / max(degree[mid], 1)
                    weight += w
                    acc += w * np.array(core['r:' + other])
        target = acc / weight if weight >= PARTNER_MIN else current
        if weight >= PARTNER_MIN:
            moved += 1
        else:
            pinned += 1
        local = np.array([np.array(p) - centre_local for p in xy.values()])
        target = _free_spot(local, np.array(target, float), fixed)
        fixed.extend((local + target).tolist())
        for n, p in xy.items():
            q = (np.array(p) - centre_local + target).tolist()
            if n.startswith('r:'):
                core[n] = q
            elif n.startswith('s:'):
                core['m:' + n[2:]] = q
        for mid, rids in region['shared'].items():
            for rid in rids:
                assign[rid, mid] = 'm:' + mid
        region['placed_near_partners'] = weight >= PARTNER_MIN
    return core, assign, dict(regions_with_shared_nodes=sum(len(r['reactions']) > 1 for r in regions), moved_to_partners=moved, kept_in_place=pinned)


def _free_spot(local, target, fixed, clearance=110.0):
    """The point nearest to `target` (on a spiral) where the region's glyphs keep `clearance` px from everything placed so far."""
    from scipy.spatial import cKDTree
    if not len(fixed):
        return target
    tree = cKDTree(np.array(fixed, float))
    radius = float(np.hypot(*(local.max(0) - local.min(0)))) / 2 + clearance
    for i in range(0, 400):
        r = 0.0 if i == 0 else 60.0 * np.sqrt(i)
        a = i * 2.399963
        c = target + r * np.array([np.cos(a), np.sin(a)])
        d, _ = tree.query(local + c)
        if float(d.min()) >= clearance:
            return c
    return target


def apply(scene, assign, core):
    """Attach the participants of grouped reactions to the shared nodes instead of private currency glyphs."""
    scene = dict(scene, positions=dict(scene['positions']), labels=dict(scene['labels']), copies=dict(scene['copies']), participants=dict(scene['participants']))
    reasons = dict(scene.get('instance_reasons', {}))
    for (rid, mid), key in assign.items():
        old = scene['participants'][rid, mid]
        if old.startswith('side:'):
            for t in (scene['positions'], scene['labels'], scene['copies']):
                t.pop(old, None)
        scene['participants'][rid, mid] = key
        scene['positions'][key] = list(core[key])
        scene['labels'][key] = mid
        scene['copies'][key] = mid
        reasons[key] = 'group_shared'
    scene['instance_reasons'] = reasons
    return scene


def isolated_reactions(model, scene):
    """Reactions that share no visual node with any other reaction (currency glyphs drawn per reaction do not count)."""
    holders = defaultdict(set)
    for (rid, mid), key in scene['participants'].items():
        holders[key].add(rid)
    out = []
    for r in model['reactions']:
        if not any(len(holders[scene['participants'][r['id'], m]]) > 1 for m in r['metabolites']):
            out.append(r['id'])
    return out
