"""V6: anchor-driven placement of the reactions that no expert structure covers.

Reactions outside the expert structures were placed by a force layout that sees only the graph (and by bridge interpolation); a reaction could
end up ~500 px from the metabolites it uses (ACKr from ac_c / actp_c, ENO from 2pg_c / pep_c, FORt from for_c / for_e) because nothing tied it to
the place where those metabolites are drawn.  The expert maps draw every reaction beside the glyphs of its own metabolites.

This module places each such reaction from its ANCHORS, after `display_scene` and before the metabolite-instance step (the expert glyphs, the
expert reaction points and the reactions that `context.py` grouped are never moved):

  direct anchor      a glyph of one of the reaction's core (non-currency) metabolites inside an expert structure (`copy:` key).  A metabolite drawn
                     in several structures has several candidates; the one nearest to the reaction's best-specified anchor (the metabolite with
                     the fewest copies) is taken, and a candidate further than ADJACENT px from it is not an anchor of this reaction.
  propagated anchor  the glyph of a core metabolite of a reaction that was already placed in an earlier round (the chain of uncovered reactions
                     is followed round by round: it goes on while the next reactions still have a close anchor, up to MAX_ROUNDS).
  weight             1 / (number of reactions that use the metabolite - 2), floored: a metabolite used by many reactions (ATP, NADH, H+) is a weak anchor.

  T1 (embed)   two or more anchors that lie within ADJACENT px of each other: the reaction goes to their weighted centre, shifted outward from the
               structure so that it sits beside it, not on top of it
  T2 (edge)    one anchor: the reaction goes beside that anchor, away from the structure; the reaction's other metabolites follow when the
               instance step finds / makes a close instance of them
  T3 (auto)    no usable anchor: left where the layout put it (automatic module)
Placement is a CANDIDATE: it is applied only if the reaction's core connections get shorter (sum of distances to the glyphs they attach to,
after the move, at most ACCEPT of before); otherwise the reaction stays (T3).  A reaction's private shared-with-nobody glyphs move with it.
Every reaction keeps exactly one glyph; GEM content is untouched.  Deterministic (sorted order, no random numbers).

`attachment_rate(model, scene, cur, expert)`: of the algorithm-placed reactions that have a core metabolite in common with another reaction, the
share that share a visual node with another reaction (computed after this step and again on the final scene, so a later step that pulls a shared
node apart shows up).
"""
from collections import Counter, defaultdict

import numpy as np
from scipy.spatial import cKDTree

ADJACENT = 400.0     # anchors of one reaction must lie within this of each other
OUT_T1, OUT_T2 = 70.0, 110.0
ACCEPT = 0.85        # the move must bring the core connections to at most this share of their old length
MAX_ROUNDS = 6
CLEAR = 75.0
COPY_REACH = 380.0   # an instance of a metabolite further than this from the reaction does not count as close: a local instance would be made
COPY_EQUIV = 90.0    # length counted for a connection to such a local instance
SAVE_PER_COPY = 250.0  # a move that needs a new local instance must save at least this much connection length per instance


def _weight(degree):
    return 1.0 / max(degree - 2, 1)


def _free_spot(target, others_tree, clearance=CLEAR):
    if others_tree is None:
        return target
    for i in range(0, 120):
        r = 0.0 if i == 0 else 40.0 * np.sqrt(i)
        a = i * 2.399963
        c = target + r * np.array([np.cos(a), np.sin(a)])
        d, _ = others_tree.query(c)
        if d >= clearance:
            return c
    return target


def attachment_rate(model, scene, cur, expert, skip=()):
    holders = defaultdict(set)
    for (rid, mid), key in scene['participants'].items():
        holders[key].add(rid)
    users = defaultdict(set)
    for r in model['reactions']:
        for mid in r['metabolites']:
            if mid not in cur:
                users[mid].add(r['id'])
    eligible = attached = 0
    for r in model['reactions']:
        rid = r['id']
        if rid in expert or rid in skip:
            continue
        core = [m for m in r['metabolites'] if m not in cur]
        if not core or not any(len(users[m]) > 1 for m in core):
            continue
        eligible += 1
        if any(len(holders[scene['participants'][rid, m]]) > 1 for m in core):
            attached += 1
    return dict(eligible=eligible, attached=attached, rate=round(attached / eligible, 3) if eligible else None)


def record_protected(model, scene, blocks, cur, tiers, regions_of, context_members=()):
    """Groups of reactions that share one visual instance of a metabolite and must keep sharing it (a group is protected, not a coordinate):
      strict       expert: reactions of one expert structure sharing a `copy:` glyph
      conditional  `completion_context`: reactions this module placed (or expert reactions) around one shared `m:` node of a core metabolite, each within
                   ADJACENT px of it; `context_group`: the shared nodes `context.py` made for reactions without core participants
    The instance step may move a conditional group to another instance as a whole, never split it."""
    parts, pos = scene['participants'], scene['positions']
    reasons = scene.get('instance_reasons', {})
    expert = {rid for b in blocks for rid in b['reaction_ids']}
    holders = defaultdict(set)
    for (rid, mid), key in parts.items():
        holders[key].add(rid)
    out = []
    for key in sorted(holders):
        rids = sorted(holders[key])
        if len(rids) < 2:
            continue
        mid = scene['copies'][key]
        if key.startswith('copy:'):
            out.append(dict(canonical_metabolite_id=mid, visual_instance_id=key, attached_reaction_ids=rids, anchor_region=key.split(':')[1], protection_reason='expert', strict=True))
        elif key.startswith('m:') and mid in cur and reasons.get(key) == 'group_shared':
            out.append(dict(canonical_metabolite_id=mid, visual_instance_id=key, attached_reaction_ids=rids, anchor_region=None, protection_reason='context_group', strict=False))
        elif key.startswith('m:') and mid not in cur:
            if not all((r in tiers or r in expert) for r in rids) or not any(r in tiers for r in rids):
                continue
            kp = np.array(pos[key], float)
            if any(float(np.hypot(*(np.array(pos['r:' + r], float) - kp))) > ADJACENT for r in rids):
                continue
            regs = Counter(regions_of.get(r) for r in rids if regions_of.get(r))
            out.append(dict(canonical_metabolite_id=mid, visual_instance_id=key, attached_reaction_ids=rids, anchor_region=regs.most_common(1)[0][0] if regs else None, protection_reason='completion_context', strict=False))
    return out


def verify_protected(scene):
    """Every protected group must still be attached to ONE visual instance in the final scene (the instance may have changed: a whole group may move).
    Updates `visual_instance_id` and returns the groups that were split."""
    broken = []
    for rec in scene.get('protected', []):
        keys = {scene['participants'][r, rec['canonical_metabolite_id']] for r in rec['attached_reaction_ids']}
        if len(keys) == 1:
            rec['visual_instance_id'] = next(iter(keys))
        else:
            broken.append(rec)
    return broken


def place(model, scene, blocks, cur, skip=()):
    """-> (new scene, report)."""
    scene = dict(scene, positions={k: list(v) for k, v in scene['positions'].items()})
    pos, parts, copies = scene['positions'], scene['participants'], scene['copies']
    expert = {rid for b in blocks for rid in b['reaction_ids']}
    rx = {r['id']: r for r in model['reactions']}
    degree = Counter(m for r in model['reactions'] for m in r['metabolites'])
    block_of = {}
    centroid = {}
    for b in blocks:
        pts = np.array([v for k, v in b['positions'].items()], float)
        centroid[b['id']] = pts.mean(0)
    copies_of = defaultdict(list)  # core metabolite -> [(key, block id)]
    for key, mid in copies.items():
        if key.startswith('copy:') and mid not in cur:
            copies_of[mid].append((key, key.split(':')[1]))
    private = defaultdict(list)    # reaction -> 'm:' keys only that reaction uses
    holders = defaultdict(set)
    for (rid, mid), key in parts.items():
        holders[key].add(rid)
    for key, rids in holders.items():
        if key.startswith('m:') and len(rids) == 1 and key[2:] not in copies_of:  # a glyph of a metabolite drawn in a structure is a fixed anchor, not the reaction's own
            private[next(iter(rids))].append(key)
    private_all = {k for ks in private.values() for k in ks}
    free = sorted(rid for rid in rx if rid not in expert and rid not in skip and any(m not in cur for m in rx[rid]['metabolites']))
    before = attachment_rate(model, scene, cur, expert, skip)

    instances_of = defaultdict(list)
    for key, mid in copies.items():
        if not key.startswith('side:'):
            instances_of[mid].append(key)

    need = []

    def edge_sum(rid, rpos, potential=False):
        """Length of the reaction's core connections.  `potential`: to the nearest instance of each metabolite (the instance step re-attaches a
        reaction to the nearest instance within reach), otherwise to the glyph the reaction is attached to now.  The reaction's own glyphs move
        with it and do not count."""
        total, own = 0.0, set(private.get(rid, []))
        for m in rx[rid]['metabolites']:
            if m in cur or parts[rid, m] in own:
                continue
            if potential:
                keys = [k for k in instances_of[m] if k in pos and (k not in private_all or k in own)]
                d = min(float(np.hypot(*(np.array(pos[k]) - rpos))) for k in keys)
                if d > COPY_REACH:  # no instance within reach: the instance step makes a recorded local instance beside the reaction
                    need.append(m)
                    d = COPY_EQUIV
                total += d
            else:
                total += float(np.hypot(*(np.array(pos[parts[rid, m]]) - rpos)))
        return total

    stats, details = Counter(), {}
    regions_of = {}
    placed, placed_nodes = set(), {}
    tiers = {}
    fixed_tree = None
    for rnd in range(1, MAX_ROUNDS + 1):
        moved_this = 0
        pending = [rid for rid in free if rid not in placed]
        # fixed obstacles: everything that will not move in this round
        movable = {f'r:{r}' for r in pending} | {k for r in pending for k in private.get(r, [])}
        others = np.array([v for k, v in pos.items() if k not in movable and not k.startswith('side:')], float)
        for rid in pending:
            r = rx[rid]
            core = sorted(m for m in r['metabolites'] if m not in cur)
            cands = []
            for m in core:
                options = [(pos[k], b) for k, b in copies_of.get(m, [])]
                cands.append((m, options))
            # round 1: structure glyphs; later rounds also glyphs of reactions placed before
            propagated = []
            if rnd > 1:
                for m in core:
                    for other in sorted(placed):
                        if m in rx[other]['metabolites']:
                            propagated.append((m, pos['r:' + other]))
            anchors = []
            ranked = sorted(((len(opts), m, opts) for m, opts in cands if opts), key=lambda t: (t[0], t[1]))
            if ranked:
                _, pm, popts = ranked[0]
                cur_pos = np.array(pos['r:' + rid], float)
                prim_opt = min(popts, key=lambda o: float(np.hypot(*(np.array(o[0]) - cur_pos))))
                prim = np.array(prim_opt[0], float)
                anchors.append((pm, prim, _weight(degree[pm]), prim_opt[1]))
                for _, m, opts in ranked[1:]:
                    best = min(opts, key=lambda o: float(np.hypot(*(np.array(o[0]) - prim))))
                    if float(np.hypot(*(np.array(best[0]) - prim))) <= ADJACENT:
                        anchors.append((m, np.array(best[0], float), _weight(degree[m]), best[1]))
            if not anchors:
                for m, p in propagated:
                    anchors.append((m, np.array(p, float), _weight(degree[m]) * .6, None))
                    break
            if not anchors:
                continue
            w = np.array([a[2] for a in anchors])
            target = (np.array([a[1] for a in anchors]) * w[:, None]).sum(0) / w.sum()
            bid = anchors[0][3]
            direction = target - centroid.get(bid, target) if bid in centroid else np.zeros(2)
            n = float(np.hypot(*direction))
            direction = direction / n if n > 1e-6 else np.array([0.0, 1.0])
            tier = 'T1' if len(anchors) >= 2 else 'T2'
            goal = target + direction * (OUT_T1 if tier == 'T1' else OUT_T2)
            tree = cKDTree(others) if len(others) else None
            goal = _free_spot(goal, tree)
            here = np.array(pos['r:' + rid], float)
            delta = goal - here
            old = edge_sum(rid, here)
            # candidate: move the reaction and its private glyphs
            saved = {k: list(pos[k]) for k in ['r:' + rid] + private.get(rid, [])}
            for k in saved:
                pos[k] = (np.array(saved[k]) + delta).tolist()
            need.clear()
            new = edge_sum(rid, goal, potential=True)
            copies_needed = len(need)
            details[rid] = dict(tier=tier, anchors=[a[0] for a in anchors], old=round(old), new=round(new), round=rnd)
            details[rid]['copies'] = copies_needed
            if old > 1e-6 and new <= ACCEPT * old and old - new >= SAVE_PER_COPY * copies_needed:
                placed.add(rid)
                tiers[rid] = tier
                regions_of[rid] = bid
                stats[tier] += 1
                moved_this += 1
                others = np.vstack([others, goal[None, :]]) if len(others) else goal[None, :]
            else:
                for k, v in saved.items():
                    pos[k] = v
                stats['kept_T3'] += 1 if rnd == 1 else 0
        stats[f'round{rnd}'] = moved_this
        if not moved_this:
            break
    after = attachment_rate(model, scene, cur, expert, skip)
    scene['completion_context'] = dict(tiers=tiers, regions=regions_of)
    scene['protected'] = record_protected(model, scene, blocks, cur, tiers, regions_of, skip)
    report = dict(candidates=len(free), placed=len(placed), tiers=dict(Counter(tiers.values())), rounds={k: v for k, v in stats.items() if k.startswith('round')},
                  attachment_rate_before=before, attachment_rate_after=after, details=details,
                  thresholds=dict(adjacent=ADJACENT, accept=ACCEPT, max_rounds=MAX_ROUNDS))
    return scene, report
