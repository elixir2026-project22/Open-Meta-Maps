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
ACCEPT_BOUNDARY = 0.95  # looser when placing against final curated AABB (after-pack pass)
MAX_ROUNDS = 6
CLEAR = 75.0
COPY_REACH = 380.0   # an instance of a metabolite further than this from the reaction does not count as close: a local instance would be made
COPY_EQUIV = 90.0    # length counted for a connection to such a local instance
SAVE_PER_COPY = 250.0  # a move that needs a new local instance must save at least this much connection length per instance
SAVE_PER_COPY_BOUNDARY = 120.0


def _weight(degree):
    return 1.0 / max(degree - 2, 1)


def _block_aabb(blocks, scene=None):
    """block id -> (lo, hi, centre) from final curated positions.

    Prefer live scene positions (composition may have translated tiles without
    rewriting ``blocks[].positions``); fall back to the block dict.
    """
    out = {}
    pos = (scene or {}).get('positions') or {}
    for b in blocks or []:
        bid = b.get('id')
        prefix = 'copy:%s:' % bid
        pts = []
        for rid in b.get('reaction_ids') or []:
            rk = 'r:' + rid
            if rk in pos:
                pts.append(pos[rk])
        for k, v in pos.items():
            if k.startswith(prefix):
                pts.append(v)
        if len(pts) < 2:
            pts = [v for _, v in (b.get('positions') or {}).items()]
        if len(pts) < 1:
            continue
        arr = np.array(pts, float)
        lo, hi = arr.min(0), arr.max(0)
        out[bid] = (lo, hi, arr.mean(0))
    return out


def _face_clearance(lo, hi, others, side):
    """How far the named face looks into empty space (min gap to other AABBs along that normal)."""
    # Infinite if no other block sits in that half-plane; else gap between faces.
    gaps = []
    if side == 'W':
        x = float(lo[0])
        for o_lo, o_hi, _ in others:
            if o_hi[1] > lo[1] and o_lo[1] < hi[1] and float(o_hi[0]) <= x + 1e-6:
                gaps.append(x - float(o_hi[0]))
    elif side == 'E':
        x = float(hi[0])
        for o_lo, o_hi, _ in others:
            if o_hi[1] > lo[1] and o_lo[1] < hi[1] and float(o_lo[0]) >= x - 1e-6:
                gaps.append(float(o_lo[0]) - x)
    elif side == 'N':
        y = float(lo[1])
        for o_lo, o_hi, _ in others:
            if o_hi[0] > lo[0] and o_lo[0] < hi[0] and float(o_hi[1]) <= y + 1e-6:
                gaps.append(y - float(o_hi[1]))
    else:
        y = float(hi[1])
        for o_lo, o_hi, _ in others:
            if o_hi[0] > lo[0] and o_lo[0] < hi[0] and float(o_lo[1]) >= y - 1e-6:
                gaps.append(float(o_lo[1]) - y)
    return min(gaps) if gaps else 1e9


def _outward_from_aabb(anchor_xy, lo, hi, prefer_ew=False, others=None, free_space=True):
    """Unit outward normal of a face of the AABB.

    Default: nearest face to ``anchor_xy``.  With ``free_space`` and neighbouring
    AABBs, only switch away from the nearest face when that face is crowded
    (clearance below a soft threshold) and another nearby face is clearly freer —
    so spokes leave into open space without picking a long detour face.
    """
    x, y = float(anchor_xy[0]), float(anchor_xy[1])
    dist = {
        'W': abs(x - float(lo[0])),
        'E': abs(float(hi[0]) - x),
        'N': abs(y - float(lo[1])),
        'S': abs(float(hi[1]) - y),
    }
    pool = ('W', 'E') if prefer_ew else ('W', 'E', 'N', 'S')
    nearest = min(pool, key=lambda s: dist[s])
    side = nearest
    if free_space and others:
        other_boxes = list(others)
        near_clear = _face_clearance(lo, hi, other_boxes, nearest)
        if near_clear < 120.0:
            band = dist[nearest] * 1.25 + 30.0
            cand = [s for s in pool if dist[s] <= band]
            best = max(cand, key=lambda s: (_face_clearance(lo, hi, other_boxes, s), -dist[s]))
            if _face_clearance(lo, hi, other_boxes, best) > near_clear + 40.0:
                side = best
    normals = {
        'W': np.array([-1.0, 0.0]),
        'E': np.array([1.0, 0.0]),
        'N': np.array([0.0, -1.0]),
        'S': np.array([0.0, 1.0]),
    }
    return normals[side], side


def _clamp_outside_aabb(goal, lo, hi, pad):
    """If ``goal`` sits inside the block AABB, push it onto the nearest outer face + pad."""
    g = np.asarray(goal, float).copy()
    if not (lo[0] - 1e-9 <= g[0] <= hi[0] + 1e-9 and lo[1] - 1e-9 <= g[1] <= hi[1] + 1e-9):
        return g
    dist = {
        'W': g[0] - lo[0],
        'E': hi[0] - g[0],
        'N': g[1] - lo[1],
        'S': hi[1] - g[1],
    }
    side = min(dist, key=dist.get)
    if side == 'W':
        g[0] = lo[0] - pad
    elif side == 'E':
        g[0] = hi[0] + pad
    elif side == 'N':
        g[1] = lo[1] - pad
    else:
        g[1] = hi[1] + pad
    return g


def _seg_crosses_aabb(x0, y0, x1, y1, lo, hi):
    """True when the open segment crosses the interior of the AABB."""
    t_enter, t_leave = 0.0, 1.0
    for a0, a1, lo_a, hi_a in ((x0, x1, float(lo[0]), float(hi[0])),
                               (y0, y1, float(lo[1]), float(hi[1]))):
        da = a1 - a0
        if abs(da) < 1e-12:
            if a0 < lo_a or a0 > hi_a:
                return False
            continue
        ta, tb = (lo_a - a0) / da, (hi_a - a0) / da
        if ta > tb:
            ta, tb = tb, ta
        t_enter = max(t_enter, ta)
        t_leave = min(t_leave, tb)
        if t_enter > t_leave:
            return False
    return t_leave > t_enter and t_leave > 1e-6 and t_enter < 1.0 - 1e-6


def _dodge_crossed_blocks(goal, target, aabb, host_bid, pad):
    """If goal→target cuts through a foreign curated AABB, sit on that block's facing side.

    Stops algorithm midmarkers from aiming through TCA/glycolysis bodies; the link then
    leaves from the facing face instead of a long diagonal through the block.
    """
    g = np.asarray(goal, float).copy()
    t = np.asarray(target, float)
    g0 = np.asarray(goal, float)
    for bid, (lo, hi, _c) in aabb.items():
        if bid == host_bid:
            continue
        if not _seg_crosses_aabb(g[0], g[1], t[0], t[1], lo, hi):
            continue
        candidates = [
            np.array([lo[0] - pad, float(np.clip(g0[1], lo[1], hi[1]))]),
            np.array([hi[0] + pad, float(np.clip(g0[1], lo[1], hi[1]))]),
            np.array([float(np.clip(g0[0], lo[0], hi[0])), lo[1] - pad]),
            np.array([float(np.clip(g0[0], lo[0], hi[0])), hi[1] + pad]),
        ]
        g = min(candidates, key=lambda p: float(np.hypot(*(p - g0))))
    return g


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


def scene_longest_links(model, scene, cur, expert, top=10):
    """Top-N reaction–metabolite centre distances from the scene (pre-export diagnosis)."""
    pos, parts = scene['positions'], scene['participants']
    rows = []
    for r in model['reactions']:
        rid = r['id']
        rk = 'r:' + rid
        if rk not in pos:
            continue
        rp = np.array(pos[rk], float)
        for mid in r['metabolites']:
            if mid in cur:
                continue
            key = parts.get((rid, mid))
            if key is None or key not in pos:
                continue
            d = float(np.hypot(*(rp - np.array(pos[key], float))))
            rows.append(dict(reaction=rid, metabolite=mid, length_px=round(d, 1),
                             curated=rid in expert, key=key))
    rows.sort(key=lambda r: -r['length_px'])
    return dict(top=rows[:top], n_edges=len(rows),
                p90=round(float(np.percentile([r['length_px'] for r in rows], 90)), 1) if rows else 0.0)


def place(model, scene, blocks, cur, skip=(), boundary_aware=False, prefer_ew=None, min_anchors=1):
    """Place free reactions beside drawn metabolite anchors.

    ``min_anchors``: require at least this many usable anchors (default 1 = T2 allowed).
    With a previous layout, pass 2 so a new reaction is embedded only between already-drawn
    participants; otherwise it stays a separate component (T3).
    """
    """-> (new scene, report).

    ``boundary_aware``: place against final curated AABB face normals (after pack), not
    radial rays from block centroids.  ``prefer_ew``: force E/W faces when the scene is
    already tall (defaults to auto from content aspect when ``boundary_aware``).
    """
    scene = dict(scene, positions={k: list(v) for k, v in scene['positions'].items()})
    pos, parts, copies = scene['positions'], scene['participants'], scene['copies']
    expert = {rid for b in blocks for rid in b['reaction_ids']}
    rx = {r['id']: r for r in model['reactions']}
    degree = Counter(m for r in model['reactions'] for m in r['metabolites'])
    block_of = {}
    centroid = {}
    aabb = _block_aabb(blocks, scene)
    for b in blocks:
        pts = np.array([v for k, v in b['positions'].items()], float)
        centroid[b['id']] = pts.mean(0)
    accept = ACCEPT_BOUNDARY if boundary_aware else ACCEPT
    save_per_copy = SAVE_PER_COPY_BOUNDARY if boundary_aware else SAVE_PER_COPY
    if prefer_ew is None and boundary_aware:
        pts = np.array([v for k, v in pos.items() if not k.startswith('side:')], float)
        if len(pts) >= 2:
            lo, hi = pts.min(0), pts.max(0)
            prefer_ew = bool((hi[0] - lo[0]) / max(hi[1] - lo[1], 1e-9) < 1.25)
        else:
            prefer_ew = False
    prefer_ew = bool(prefer_ew)
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
            if not anchors or len(anchors) < max(1, int(min_anchors)):
                continue
            w = np.array([a[2] for a in anchors])
            target = (np.array([a[1] for a in anchors]) * w[:, None]).sum(0) / w.sum()
            bid = anchors[0][3]
            tier = 'T1' if len(anchors) >= 2 else 'T2'
            out = OUT_T1 if tier == 'T1' else OUT_T2
            attach_side = None
            if boundary_aware and bid in aabb:
                lo, hi, _c = aabb[bid]
                neighbour = [box for b2, box in aabb.items() if b2 != bid]
                direction, attach_side = _outward_from_aabb(
                    target, lo, hi, prefer_ew=prefer_ew, others=neighbour, free_space=True)
                goal = target + direction * out
                goal = _clamp_outside_aabb(goal, lo, hi, out)
                # If the chord still cuts another curated block, sit on that block's face.
                goal = _dodge_crossed_blocks(goal, target, aabb, bid, out)
            else:
                direction = target - centroid.get(bid, target) if bid in centroid else np.zeros(2)
                n = float(np.hypot(*direction))
                direction = direction / n if n > 1e-6 else np.array([0.0, 1.0])
                goal = target + direction * out
            tree = cKDTree(others) if len(others) else None
            goal = _free_spot(goal, tree)
            if boundary_aware and bid in aabb:
                lo, hi, _c = aabb[bid]
                goal = _clamp_outside_aabb(goal, lo, hi, out * 0.6)
                goal = _dodge_crossed_blocks(goal, target, aabb, bid, out * 0.6)
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
            if attach_side is not None:
                details[rid]['attach_side'] = attach_side
            if old > 1e-6 and new <= accept * old and old - new >= save_per_copy * copies_needed:
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
    # After boundary-aware rounds: force any algorithm midmarker still inside a curated
    # AABB out onto the nearest outer face (cuts algorithm×curated crossings through blocks).
    ejected = 0
    if boundary_aware and aabb:
        others = np.array([v for k, v in pos.items() if not k.startswith('side:')], float)
        tree = cKDTree(others) if len(others) else None
        for rid in free:
            rk = 'r:' + rid
            if rk not in pos:
                continue
            here = np.array(pos[rk], float)
            host = None
            for bid, (lo, hi, _c) in aabb.items():
                if lo[0] - 1e-9 <= here[0] <= hi[0] + 1e-9 and lo[1] - 1e-9 <= here[1] <= hi[1] + 1e-9:
                    host = bid
                    break
            if host is None:
                continue
            lo, hi, _c = aabb[host]
            goal = _clamp_outside_aabb(here, lo, hi, OUT_T2)
            goal = _free_spot(goal, tree)
            goal = _clamp_outside_aabb(goal, lo, hi, OUT_T2 * 0.6)
            delta = goal - here
            for k in [rk] + private.get(rid, []):
                pos[k] = (np.array(pos[k], float) + delta).tolist()
            details.setdefault(rid, {})['ejected_from_block'] = host
            ejected += 1
            if tree is not None:
                others = np.vstack([others, goal[None, :]]) if len(others) else goal[None, :]
                tree = cKDTree(others)
    after = attachment_rate(model, scene, cur, expert, skip)
    scene['completion_context'] = dict(tiers=tiers, regions=regions_of)
    scene['protected'] = record_protected(model, scene, blocks, cur, tiers, regions_of, skip)
    report = dict(candidates=len(free), placed=len(placed), tiers=dict(Counter(tiers.values())), rounds={k: v for k, v in stats.items() if k.startswith('round')},
                  attachment_rate_before=before, attachment_rate_after=after, details=details,
                  boundary_aware=bool(boundary_aware), prefer_ew=bool(prefer_ew) if prefer_ew is not None else None,
                  ejected_from_blocks=ejected,
                  thresholds=dict(adjacent=ADJACENT, accept=accept, max_rounds=MAX_ROUNDS,
                                  save_per_copy=save_per_copy))
    return scene, report
