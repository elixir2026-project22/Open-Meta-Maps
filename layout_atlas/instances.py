"""V3: visual instances of metabolites, decided from the global layout.

Until now a metabolite outside the expert structures had ONE visual node (`m:<id>`), however many modules used it.  In the iMM904 global map
66 % of the connections longer than 1,000 px ran to ~30 such metabolites (glutamate, acetyl-CoA, 2-oxoglutarate, pyruvate ...) and 29 to the
biomass reaction.  Expert maps draw such a metabolite once per place where it is used.  This module decides, from the positions that the
global layout produced, how many visual instances each metabolite gets and which reactions each instance serves.  One GEM metabolite stays one
identity (`copies[key]` is always its id); instances are only places where it is drawn.

Three operations on the scene made by `drawing.display_scene` (the hashed drawing files are not modified):

  summary reactions   a reaction outside the expert structures with >= SUMMARY_MIN non-currency participants (biomass) gets a private instance
                      of every participant on a ring around it, instead of one long line per participant to nodes all over the map.  Every
                      participant stays drawn and connected; nothing is hidden or merged.  Reason `summary_reaction`.
  reuse / merge       a reaction attaches to the nearest existing instance of its metabolite (an expert copy, a shared node or an earlier
                      instance) when that is within NEAR px.  A shared node that no longer serves any reaction disappears (a merge).
  split               reactions of one metabolite that are further than NEAR from every existing instance are grouped (single linkage,
                      LINK px); a group of >= MIN_CLUSTER reactions (or one reaction farther than FAR) gets a new instance at its centre.
                      Reason `long_edge_split`.  At most MAX_NEW new instances per metabolite.

Expert glyphs and reaction points are never moved here; only the participant each non-expert reaction attaches to is chosen.  Currency
metabolites keep their per-reaction `side:` glyphs.  Deterministic.  The report gives long-edge counts before / after (on the scene before
relaxation, participant edges of non-currency metabolites) and the reduction per added instance.
"""
from collections import Counter, defaultdict

import numpy as np
from scipy.spatial import cKDTree

NEAR = 250.0
LINK = 260.0
FAR = 1000.0
MIN_CLUSTER = 2
MAX_NEW = 8
SUMMARY_MIN = 12          # absolute floor for summary-reaction treatment
SUMMARY_DEGREE_PCT = 0.92 # or: reactions whose non-currency degree is at/above this percentile of the model
import os as _os
SUMMARY_LAYOUT = _os.environ.get('OPM_SUMMARY_LAYOUT', 'columns')  # columns fan locally; ring radiates in every direction
SUMMARY_SPACING = 34.0
SPLIT_PLACED = 380.0  # a single reaction that anchor-driven placement put beside a structure gets its own recorded instance of a metabolite that is further than this
RESCUE = 700.0      # an otherwise isolated reaction attaches to the nearest instance used by another reaction up to this far
# High-degree non-currency metabolites (shared across many reactions) get more aggressive local cloning
HUB_METAB_PCT = 0.95
HUB_NEAR = 180.0
HUB_MAX_NEW = 12


def _long(model, scene, cur, thresholds=(250.0, 1000.0)):
    pos, parts = scene['positions'], scene['participants']
    counts = [0] * len(thresholds)
    for r in model['reactions']:
        rp = np.array(pos['r:' + r['id']], float)
        for mid in r['metabolites']:
            if mid in cur:
                continue
            d = float(np.hypot(*(np.array(pos[parts[r['id'], mid]], float) - rp)))
            for i, t in enumerate(thresholds):
                counts[i] += d > t
    return counts


def _clusters(points, link):
    n = len(points)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    if n > 1:
        for i, j in cKDTree(points).query_pairs(link):
            parent[find(i)] = find(j)
    groups = defaultdict(list)
    for i in range(n):
        groups[find(i)].append(i)
    return sorted(groups.values(), key=lambda g: g[0])


def _degree_threshold(model, cur, pct, floor):
    """Data-driven degree cutoff: max(floor, percentile of non-currency participant counts)."""
    degrees = [sum(1 for m in r['metabolites'] if m not in cur) for r in model['reactions']]
    if not degrees:
        return floor
    return max(floor, int(np.percentile(degrees, pct * 100)))


def _hub_metabolites(model, cur, pct):
    """Metabolites used by many reactions (top percentile of usage counts), excluding currency."""
    uses = defaultdict(int)
    for r in model['reactions']:
        for mid in r['metabolites']:
            if mid not in cur:
                uses[mid] += 1
    if not uses:
        return set()
    cut = max(8, int(np.percentile(list(uses.values()), pct * 100)))
    return {m for m, n in uses.items() if n >= cut}


def optimize(model, scene, blocks, cur, summary=True, reuse=True, split=True):
    """Returns (new scene, report).  `scene` is not modified."""
    scene = dict(scene, positions={k: list(v) for k, v in scene['positions'].items()}, labels=dict(scene['labels']), copies=dict(scene['copies']),
                 participants=dict(scene['participants']))
    pos, labels, copies, parts = scene['positions'], scene['labels'], scene['copies'], scene['participants']
    reasons = {}
    expert = {rid for b in blocks for rid in b['reaction_ids']}
    before = _long(model, scene, cur)
    rx = {r['id']: r for r in model['reactions']}
    shared_before = {k for k in pos if k.startswith('m:')}
    summary_cut = _degree_threshold(model, cur, SUMMARY_DEGREE_PCT, SUMMARY_MIN)
    hub_mets = _hub_metabolites(model, cur, HUB_METAB_PCT)
    report = dict(summary_reactions=[], added=0, reattached=0, summary_degree_cut=summary_cut, hub_metabolites=len(hub_mets))

    # 1. summary reactions: local stubs/columns instead of map-wide spokes (degree-driven, not name-driven)
    if summary:
        for rid in sorted(rx):
            r = rx[rid]
            core = [m for m in sorted(r['metabolites']) if m not in cur and not parts[rid, m].startswith('copy:')]
            if rid in expert or len(core) < summary_cut:
                continue
            centre = np.array(pos['r:' + rid], float)
            radius = max(160.0, len(core) * 46 / (2 * np.pi))
            subs = [m for m in core if r['metabolites'][m] < 0]
            prods = [m for m in core if r['metabolites'][m] >= 0]
            if SUMMARY_LAYOUT == 'columns':
                # two columns: substrates in one column on the left of the reaction point, products in one on the right, each ordered by where the
                # metabolite is drawn elsewhere (top to bottom), so the connections form two fans instead of a ring (a ring sends lines in every direction)
                spacing = SUMMARY_SPACING
                half = max(len(subs), len(prods))
                dist = max(200.0, 90.0 + half * spacing * 0.35)

                def yorder(mid):
                    return (pos['m:' + mid][1] if 'm:' + mid in pos else centre[1], mid)
                for members, sign in ((sorted(subs, key=yorder), -1.0), (sorted(prods, key=yorder), 1.0)):
                    n = len(members)
                    for i, mid in enumerate(members):
                        key = f'inst:{rid}:{mid}'
                        pos[key] = [float(centre[0] + sign * dist), float(centre[1] + (i - (n - 1) / 2.0) * spacing)]
                        labels[key], copies[key] = mid, mid
                        parts[rid, mid] = key
                        reasons[key] = 'summary_reaction'
                report['summary_reactions'].append(dict(reaction=rid, participants=len(core), layout='columns'))
                report['added'] += len(core)
                continue
            for members, start in ((subs, 90.0), (prods, -90.0)):
                for i, mid in enumerate(members):
                    angle = np.radians(start + 180.0 * (i + .5) / max(len(members), 1))
                    key = f'inst:{rid}:{mid}'
                    pos[key] = (centre + radius * np.array([np.cos(angle), np.sin(angle)])).tolist()
                    labels[key], copies[key] = mid, mid
                    parts[rid, mid] = key
                    reasons[key] = 'summary_reaction'
            report['summary_reactions'].append(dict(reaction=rid, participants=len(core)))
            report['added'] += len(core)

    # 2. reuse and split, per metabolite
    uses = defaultdict(list)
    for (rid, mid), key in parts.items():
        if key == 'm:' + mid and mid not in cur:
            uses[mid].append(rid)
    by_mid = defaultdict(list)
    for key, mid in copies.items():
        if not key.startswith('side:') and reasons.get(key) != 'summary_reaction':
            by_mid[mid].append(key)
    per_metabolite = {}
    group_of = {}      # (reaction, metabolite) -> protected record (conditional groups on a shared 'm:' node of a core metabolite)
    for rec in scene.get('protected', []):
        if not rec.get('strict') and rec['visual_instance_id'].startswith('m:') and rec['canonical_metabolite_id'] not in cur:
            for r in rec['attached_reaction_ids']:
                group_of[r, rec['canonical_metabolite_id']] = rec
    moved_groups = Counter()
    for mid in sorted(uses):
        rids = sorted(uses[mid])
        targets = [k for k in by_mid[mid] if k in pos]
        if not targets:
            continue
        # Hub metabolites (top percentile of reaction usage) clone more aggressively so they are not map-wide spokes
        near = HUB_NEAR if mid in hub_mets else NEAR
        max_new = HUB_MAX_NEW if mid in hub_mets else MAX_NEW
        rp = np.array([pos['r:' + rid] for rid in rids], float)
        added = 0
        loose = []
        # protected groups first: identity (same metabolite) -> group compatible -> geometry acceptable -> the whole group reuses an instance or gets one new
        done, handled = set(), set()
        for rid in rids:
            rec = group_of.get((rid, mid))
            if rec is None or id(rec) in done:
                continue
            done.add(id(rec))
            members = [r for r in rec['attached_reaction_ids'] if parts.get((r, mid)) == 'm:' + mid]
            if len(members) < 2:
                continue
            handled.update(members)
            mp = np.array([pos['r:' + r] for r in members], float)
            here = np.array(pos['m:' + mid], float)
            cost_now = float(np.sum(np.hypot(*(mp - here).T)))
            options = []
            for t in targets:
                if t == 'm:' + mid:
                    continue
                if t.startswith('copy:') and rec.get('anchor_region') and t.split(':')[1] != rec['anchor_region']:
                    continue  # a glyph of another expert structure is not this group's anchor
                d = np.hypot(*(mp - np.array(pos[t], float)).T)
                if float(d.max()) <= near and float(d.sum()) <= cost_now:
                    options.append((float(d.sum()), t))
            if options:
                best = min(options)[1]
            elif cost_now > 0 and float(np.min(np.hypot(*(mp - here).T))) > SPLIT_PLACED and added < max_new:
                best = f'inst:{mid}:{added}'
                pos[best] = (mp.mean(axis=0) + np.array([0.0, 55.0])).tolist()
                labels[best], copies[best] = mid, mid
                reasons[best] = 'long_edge_split'
                targets.append(best)
                added += 1
                report['added'] += 1
            else:
                best = 'm:' + mid
                moved_groups['kept'] += 1
            if best != 'm:' + mid:
                for r in members:
                    parts[r, mid] = best
                moved_groups['moved_as_a_whole'] += 1
        rids = [r for r in rids if r not in handled]
        rp = np.array([pos['r:' + rid] for rid in rids], float) if rids else np.zeros((0, 2))
        for rid, p in zip(rids, rp):
            tp = np.array([pos[t] for t in targets], float)
            d = np.hypot(*(tp - p).T)
            j = int(np.argmin(d))
            if reuse and d[j] <= near:
                if targets[j] != 'm:' + mid:
                    report['reattached'] += 1
                parts[rid, mid] = targets[j]
            else:
                loose.append((rid, p))
        if loose and split:
            pts = np.array([p for _, p in loose], float)
            for cluster in _clusters(pts, LINK):
                members = [loose[i] for i in cluster]
                centre = np.mean([p for _, p in members], axis=0)
                tp = np.array([pos[t] for t in targets], float)
                mean_d = float(np.mean([np.min(np.hypot(*(tp - p).T)) for _, p in members]))
                placed_here = scene.get('completion_context', {}).get('tiers', {})
                if added < max_new and (len(members) >= MIN_CLUSTER or mean_d > FAR or (mean_d > SPLIT_PLACED and any(rid in placed_here for rid, _ in members))):
                    key = f'inst:{mid}:{added}'
                    pos[key] = (centre + np.array([0.0, 55.0])).tolist()
                    labels[key], copies[key] = mid, mid
                    reasons[key] = 'long_edge_split'
                    targets.append(key)
                    added += 1
                    for rid, _ in members:
                        parts[rid, mid] = key
                    report['added'] += 1
                else:
                    for rid, p in members:
                        d = np.hypot(*(np.array([pos[t] for t in targets], float) - p).T)
                        parts[rid, mid] = targets[int(np.argmin(d))]
        if added:
            per_metabolite[mid] = added

    # 2b. a reaction outside the expert structures whose every glyph is private (shared with no other reaction) is attached to the nearest
    # instance that another reaction uses, whatever the usual NEAR rule says: an island is worse than a longer line
    holders = defaultdict(set)
    for (rid, mid), key in parts.items():
        holders[key].add(rid)
    rescued = 0
    for rid in sorted(rx):
        if rid in expert:
            continue
        mine = [(mid, parts[rid, mid]) for mid in sorted(rx[rid]['metabolites']) if mid not in cur]
        if not mine or any(len(holders[k]) > 1 for _, k in mine):
            continue
        rp = np.array(pos['r:' + rid], float)
        done = False
        for mid, k in mine:
            options = [t for t, m2 in copies.items() if m2 == mid and not t.startswith('side:') and reasons.get(t) != 'summary_reaction' and t != k and holders[t] - {rid}]
            if not options:
                continue
            dist = [float(np.hypot(*(np.array(pos[t], float) - rp))) for t in options]
            j = int(np.argmin(dist))
            if dist[j] <= RESCUE:
                holders[k].discard(rid)
                parts[rid, mid] = options[j]
                holders[options[j]].add(rid)
                done = True
        rescued += done
    report['isolated_reactions_attached'] = rescued

    # 3. shared nodes that no reaction uses any more were merged into other instances
    used = set(parts.values())
    removed = 0
    for key in sorted(shared_before):
        if key not in used and key in pos:
            for table in (pos, labels, copies):
                table.pop(key, None)
            removed += 1
    scene['instance_reasons'] = reasons
    after = _long(model, scene, cur)
    gained = before[0] - after[0]
    report['protected_groups'] = dict(moved_groups)
    report.update(shared_nodes_merged=removed, reactions_reattached=report.pop('reattached'),
                  edges_over_250=dict(before=int(before[0]), after=int(after[0])), edges_over_1000=dict(before=int(before[1]), after=int(after[1])),
                  long_edge_reduction_per_added_instance=round(gained / max(report['added'], 1), 2),
                  metabolites_split=len(per_metabolite), top_split=sorted(per_metabolite.items(), key=lambda kv: (-kv[1], kv[0]))[:15],
                  thresholds=dict(near=NEAR, link=LINK, far=FAR, min_cluster=MIN_CLUSTER, max_new=MAX_NEW, summary_min=SUMMARY_MIN))
    return scene, report


def overlap_breakdown(scene, radius=10.0):
    """The glyph pairs closer than `radius` (the node-overlap count of the export), by kind: real conflicts between different glyphs; `same_metabolite` are two
    instances of one metabolite lying on each other (a merge candidate); a shared node is ONE glyph and never an overlap."""
    keys = sorted(scene['positions'])
    xy = np.array([scene['positions'][k] for k in keys], float)
    out = Counter()
    for a, b in cKDTree(xy).query_pairs(radius):
        ka, kb = keys[a], keys[b]
        kind = tuple(sorted(k.split(':')[0] for k in (ka, kb)))
        same = ka[0] != 'r' and kb[0] != 'r' and scene['copies'].get(ka) == scene['copies'].get(kb)
        out['|'.join(kind) + ('|same_metabolite' if same else '')] += 1
    return dict(out)


def separate(model, scene, blocks, metabolite_gap=88.0, reaction_gap=70.0, passes=80, side_gap=None):
    """Keep the instances this module placed (`inst:`) and the shared nodes (`m:`) from landing on other glyphs: labels are 60-80 px wide, so
    two metabolite glyphs closer than `metabolite_gap`, or a metabolite glyph and a reaction it is not connected to closer than `reaction_gap`,
    have overlapping labels.  Only `inst:` and `m:` glyphs move (never expert copies or reactions).

    Currency ``side:`` glyphs use the same adaptive local-clearance rule as ``orientation.sides`` unless ``side_gap`` is given explicitly.
    """
    import orientation as orientation_module
    scene = dict(scene, positions={k: list(v) for k, v in scene['positions'].items()})
    pos = scene['positions']
    keys = sorted(pos)
    index = {k: i for i, k in enumerate(keys)}
    xy = np.array([pos[k] for k in keys], float)
    movable = np.array([k.startswith(('inst:', 'm:', 'side:')) for k in keys])
    is_side = np.array([k.startswith('side:') for k in keys])
    is_rx = np.array([k.startswith('r:') for k in keys])
    if side_gap is None:
        params = scene.get('orientation_spoke_params') or orientation_module.reaction_spoke_params(model, scene)
        side_need = np.array([
            float(params.get(k.split(':', 2)[1], {}).get('gap', orientation_module.side_gap_for(orientation_module.SIDE_MAX)))
            if k.startswith('side:') else 0.0 for k in keys
        ])
    else:
        side_need = np.full(len(keys), float(side_gap))
    linked = defaultdict(set)
    for (rid, mid), key in scene['participants'].items():
        a, b = index['r:' + rid], index[key]
        linked[a].add(b)
        linked[b].add(a)
    start = xy.copy()
    query_r = max(metabolite_gap, reaction_gap, float(side_need.max()) if len(side_need) else 0.0)
    for _ in range(passes):
        pairs = cKDTree(xy).query_pairs(query_r, output_type='ndarray')
        push = np.zeros_like(xy)
        hit = False
        for i, j in pairs:
            if not (movable[i] or movable[j]) or (is_rx[i] and is_rx[j]):
                continue
            if is_side[i] or is_side[j]:
                need = min(side_need[i] or side_need[j], side_need[j] or side_need[i])
                if is_rx[i] != is_rx[j]:
                    a, b = (i, j) if is_rx[i] else (j, i)
                    if b in linked[a]:
                        continue
            elif is_rx[i] != is_rx[j]:
                a, b = (i, j) if is_rx[i] else (j, i)
                if b in linked[a]:
                    continue
                need = reaction_gap
            else:
                need = metabolite_gap
            d = xy[j] - xy[i]
            dist = float(np.hypot(*d))
            if dist >= need:
                continue
            if dist < 1e-6:
                ang = (i * 17 + j * 31) * 2.399963
                d, dist = np.array([np.cos(ang), np.sin(ang)]), 1.0
            step = d / dist * (need - dist)
            hit = True
            if movable[i] and movable[j]:
                push[i] -= step / 2
                push[j] += step / 2
            elif movable[i]:
                push[i] -= step
            else:
                push[j] += step
        if not hit:
            break
        xy += push * .7
    moved = int((np.hypot(*(xy - start).T) > 1.0).sum())
    for k in keys:
        pos[k] = xy[index[k]].tolist()
    scene['instances_separated'] = moved
    return scene
