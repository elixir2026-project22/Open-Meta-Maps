"""Layout quality measures that feed layout decisions (grown version by version; see Proposal/14).

Also contains metric_details and assert_chemistry, moved from geometry_experiment to avoid
loading experiment/UI modules during compose.

V1 (reaction orientation) needs:
  orientation   share of reactions with substrates and products on opposite sides of the reaction point (orientation.summary)
  label overlap pairs / node overlaps (drawing.export / geometry_experiment.metric_details, already in the map metric)
  expert geometry   how far the layout moved the expert structures from where they were placed (similarity-aligned residual of the reaction
                    points of each structure: translation, rotation and uniform scale are not counted as distortion)
  integrity       the three integrity checks of the pipeline (unchanged)
"""
import numpy as np

import source_paths as sp


def expert_geometry(blocks, scene):
    """Distortion of expert structures between their planned positions (`block['positions']`) and the final scene."""
    rows = []
    for b in blocks:
        src, dst = [], []
        for rid in b['reaction_ids']:
            k = 'r:' + rid
            if k in b['positions'] and k in scene['positions']:
                src.append(b['positions'][k])
                dst.append(scene['positions'][k])
        if len(src) < 3:
            continue
        fit = sp.similarity(src, dst)
        if fit is None:
            continue
        shift = float(np.mean(np.hypot(*(np.array(dst) - np.array(src)).T)))
        rows.append((fit['rms'], fit['max'], shift, len(src)))
    if not rows:
        return dict(structures=0)
    rms = np.array([r[0] for r in rows])
    weights = np.array([r[3] for r in rows], float)
    return dict(structures=len(rows), distortion_rms_px_mean=round(float(np.average(rms, weights=weights)), 2), distortion_rms_px_max=round(float(rms.max()), 2),
                distortion_max_px=round(float(max(r[1] for r in rows)), 2), mean_displacement_px=round(float(np.average([r[2] for r in rows], weights=weights)), 1),
                note='residual after translation, rotation and uniform scaling of each structure')


# ---------------------------------------------------------------------------------------------- final-geometry evaluator (V7, minimal version)
# Everything below is computed from the FINAL map.json (the same segments, Bezier controls and label offsets that the SVG draws), never from the layout
# stage's straight lines or rectangles.  Crossings are counted on the sampled final curves, by source: algorithm x algorithm, algorithm x expert,
# expert x expert; two connections that share a node (same reaction, same glyph) do not cross each other legitimately.
from collections import defaultdict as _dd

CELL = 200.0


def _grid_add(grid, idx, lo, hi):
    for gx in range(int(lo[0] // CELL), int(hi[0] // CELL) + 1):
        for gy in range(int(lo[1] // CELL), int(hi[1] // CELL) + 1):
            grid[gx, gy].append(idx)


def _grid_query(grid, lo, hi):
    out = set()
    for gx in range(int(lo[0] // CELL), int(hi[0] // CELL) + 1):
        for gy in range(int(lo[1] // CELL), int(hi[1] // CELL) + 1):
            out.update(grid.get((gx, gy), ()))
    return out


def _final_obstacles(data):
    """Glyph and label boxes with kind ('node' / 'label') and the node id each belongs to."""
    import drawing
    font = drawing.FONT
    body = data[1]
    rects, owner, kind, label_of = [], [], [], []
    for nid, n in body['nodes'].items():
        t = n.get('node_type')
        if t == 'metabolite':
            rects.append((n['x'] - er_NODE, n['y'] - er_NODE, n['x'] + er_NODE, n['y'] + er_NODE))
            owner.append(nid); kind.append('node'); label_of.append(n.get('bigg_id'))
            if 'label_x' in n and not n.get('label_hidden'):
                w = font.getlength(n.get('bigg_id', ''))
                rects.append((n['label_x'] - 2, n['label_y'] - 14, n['label_x'] + w + 2, n['label_y'] + 4))
                owner.append(nid); kind.append('label'); label_of.append(n.get('bigg_id'))
        elif t == 'midmarker':
            rects.append((n['x'] - er_NODE + 1, n['y'] - er_NODE + 1, n['x'] + er_NODE - 1, n['y'] + er_NODE - 1))
            owner.append(nid); kind.append('node'); label_of.append(None)
    for rx in body['reactions'].values():
        mid = next((s['from_node_id'] for s in rx['segments'].values() if body['nodes'][s['from_node_id']].get('node_type') == 'midmarker'), None)
        if mid is None or 'label_x' not in rx or rx.get('label_hidden'):
            continue
        w = font.getlength(rx['bigg_id'])
        rects.append((rx['label_x'] - 2, rx['label_y'] - 14, rx['label_x'] + w + 2, rx['label_y'] + 4))
        owner.append(mid); kind.append('label'); label_of.append(rx['bigg_id'])
    return np.array(rects, float).reshape(-1, 4), owner, kind, label_of


er_NODE = 7.0


def _seg_hits(a, b, S):
    """Indices (into S [k,4]) of the segments properly crossed by a-b, with intersection points."""
    if not len(S):
        return [], np.zeros((0, 2))
    p, r = np.array(a, float), np.array(b, float) - np.array(a, float)
    q, s = S[:, :2], S[:, 2:] - S[:, :2]
    den = r[0] * s[:, 1] - r[1] * s[:, 0]
    with np.errstate(divide='ignore', invalid='ignore'):
        t = ((q[:, 0] - p[0]) * s[:, 1] - (q[:, 1] - p[1]) * s[:, 0]) / den
        u = ((q[:, 0] - p[0]) * r[1] - (q[:, 1] - p[1]) * r[0]) / den
    ok = (np.abs(den) > 1e-9) & (t > 1e-3) & (t < 1 - 1e-3) & (u > 1e-3) & (u < 1 - 1e-3)
    idx = np.flatnonzero(ok)
    return idx, p + np.outer(t[idx], r)


def evaluate(data, skip_edges=None, orientation_kinds=None, max_issues=4000):
    """-> (summary, issues).  `summary` holds the six core indicators, the details and the worst local regions; `issues` lists where the problems are."""
    import edge_routing as er
    import line_render
    import orientation as orientation_module
    meta = data[0].get('metadata', {})
    skip = {tuple(e) for e in (skip_edges if skip_edges is not None else meta.get('preserved_edges') or [])}
    skip |= {(p['reaction'], p['metabolite']) for p in (meta.get('source_paths') or [])}
    routed = {(p['reaction'], p['metabolite']): p for p in (meta.get('routed_edges') or [])}
    body = data[1]
    nodes, reactions = body['nodes'], body['reactions']
    rects, owner, kind, label_of = _final_obstacles(data)
    owner_arr = np.array(owner, dtype=object)
    grid_r = _dd(list)
    for i, r in enumerate(rects):
        _grid_add(grid_r, i, r[:2], r[2:])
    edges = []
    for rx in reactions.values():
        for path in line_render._terminals(rx, nodes):
            if nodes[path[-1]].get('node_type') != 'metabolite':
                continue
            pts = np.array(er._path_points(rx, path, nodes), float)
            segs = [s for s in rx['segments'].values() if s['from_node_id'] in path and s['to_node_id'] in path]
            curved = [bool(s.get('b1') and s.get('b2')) for s in segs]
            edges.append(dict(rid=rx['bigg_id'], met=nodes[path[-1]].get('bigg_id'), path=path, pts=pts, expert=(rx['bigg_id'], nodes[path[-1]].get('bigg_id')) in skip,
                              curved=curved, ends=(path[0], path[-1])))
    algo = [e for e in edges if not e['expert']]
    # ---- crossings on the sampled final curves
    seg_a, seg_b, seg_edge = [], [], []
    for i, e in enumerate(edges):
        for a, b in zip(e['pts'][:-1], e['pts'][1:]):
            seg_a.append(a); seg_b.append(b); seg_edge.append(i)
    S = np.array([[a[0], a[1], b[0], b[1]] for a, b in zip(seg_a, seg_b)], float).reshape(-1, 4)
    seg_edge = np.array(seg_edge)
    grid_s = _dd(list)
    for k, (a, b) in enumerate(zip(seg_a, seg_b)):
        _grid_add(grid_s, k, np.minimum(a, b), np.maximum(a, b))
    cross = dict(algorithm_algorithm=0, algorithm_expert=0, expert_expert=0, same_reaction=0, same_reaction_expert=0)
    cross_points, per_edge_cross = [], _dd(int)
    for k, (a, b) in enumerate(zip(seg_a, seg_b)):
        i = seg_edge[k]
        cand = np.array(sorted(c for c in _grid_query(grid_s, np.minimum(a, b), np.maximum(a, b)) if seg_edge[c] > i), dtype=int)
        if not len(cand):
            continue
        idx, pts = _seg_hits(a, b, S[cand])
        for c, pt in zip(cand[idx], pts):
            j = seg_edge[c]
            if edges[i]['rid'] == edges[j]['rid']:      # curves of one reaction meet at the reaction point; a crossing elsewhere is a defect
                centre = nodes[edges[i]['ends'][0]]
                if float(np.hypot(pt[0] - centre['x'], pt[1] - centre['y'])) > 14.0:
                    both_expert = edges[i]['expert'] and edges[j]['expert']
                    cross['same_reaction_expert' if both_expert else 'same_reaction'] += 1
                    per_edge_cross[i] += 1; per_edge_cross[j] += 1
                    if len(cross_points) < max_issues:
                        cross_points.append(dict(kind='same_reaction_expert' if both_expert else 'same_reaction', reactions=[edges[i]['rid']], metabolites=[edges[i]['met'], edges[j]['met']], x=round(float(pt[0]), 1), y=round(float(pt[1]), 1)))
                continue
            if set(edges[i]['ends']) & set(edges[j]['ends']):
                continue
            ei, ej = edges[i]['expert'], edges[j]['expert']
            kindc = 'expert_expert' if ei and ej else ('algorithm_algorithm' if not ei and not ej else 'algorithm_expert')
            cross[kindc] += 1
            per_edge_cross[i] += 1; per_edge_cross[j] += 1
            if len(cross_points) < max_issues:
                cross_points.append(dict(kind=kindc, reactions=[edges[i]['rid'], edges[j]['rid']], x=round(float(pt[0]), 1), y=round(float(pt[1]), 1)))
    # ---- obstacles on the final curves (algorithm connections only)
    node_hits, label_hits, hit_points = 0, 0, []
    for e in algo:
        mask_own = np.isin(owner_arr, e['ends'])
        found_node = found_label = False
        n = len(e['pts']) - 1
        for si, (a, b) in enumerate(zip(e['pts'][:-1], e['pts'][1:])):
            cand = np.array(sorted(_grid_query(grid_r, np.minimum(a, b) - 2, np.maximum(a, b) + 2)), dtype=int)
            cand = cand[~mask_own[cand]] if len(cand) else cand
            if not len(cand):
                continue
            hit = (er._blocked(a, b, rects[cand], 0.0) if (si == 0 or si == n - 1) else er._hits(a, b, rects[cand], 0.0))
            for c in cand[np.flatnonzero(hit)]:
                if kind[c] == 'node' and not found_node:
                    found_node = True
                    hit_points.append(dict(kind='node', reaction=e['rid'], metabolite=e['met'], x=round(float(a[0]), 1), y=round(float(a[1]), 1)))
                if kind[c] == 'label' and not found_label:
                    found_label = True
                    hit_points.append(dict(kind='label', reaction=e['rid'], metabolite=e['met'], label=label_of[c], x=round(float(a[0]), 1), y=round(float(a[1]), 1)))
        node_hits += int(found_node); label_hits += int(found_label)
    # ---- glyph overlaps and label collisions
    from scipy.spatial import cKDTree
    centres = np.array([[n['x'], n['y']] for n in nodes.values() if n.get('node_type') in ('metabolite', 'midmarker')], float)
    node_overlap = list(cKDTree(centres).query_pairs(10)) if len(centres) > 1 else []
    lab_idx = [i for i, k in enumerate(kind) if k == 'label']
    label_overlap = 0
    for ii, i in enumerate(lab_idx):
        for j in _grid_query(grid_r, rects[i][:2], rects[i][2:]):
            if j > i and kind[j] == 'label' and rects[j][0] < rects[i][2] and rects[i][0] < rects[j][2] and rects[j][1] < rects[i][3] and rects[i][1] < rects[j][3]:
                label_overlap += 1
    # ---- length, detour, Bezier representation, curve anomalies
    lengths, detours = [], []
    bez_edges = poly_edges = ortho_edges = straight_edges = hooks = sharp = selfx = 0
    anomalies = []
    for e in algo:
        d = np.diff(e['pts'], axis=0)
        L = float(np.sum(np.hypot(d[:, 0], d[:, 1])))
        chord = float(np.hypot(*(e['pts'][-1] - e['pts'][0])))
        lengths.append(L)
        if chord >= 40:
            detours.append(L / chord)
        if e['curved'] and all(e['curved']):
            bez_edges += 1
            worst, total = _turn_stats(e['pts'])
            if total > 200:
                hooks += 1; anomalies.append(dict(kind='hook', reaction=e['rid'], metabolite=e['met']))
            elif worst > 70:
                sharp += 1; anomalies.append(dict(kind='sharp_turn', reaction=e['rid'], metabolite=e['met']))
            if _self_cross(e['pts']):
                selfx += 1; anomalies.append(dict(kind='self_intersection', reaction=e['rid'], metabolite=e['met']))
        else:
            if len(e['curved']) == 1:
                straight_edges += 1
            else:
                poly_edges += 1
                if _orthogonal(e['pts']):
                    ortho_edges += 1
    all_len = np.array([float(np.sum(np.hypot(*np.diff(e['pts'], axis=0).T))) for e in edges]) if edges else np.zeros(1)
    mids = np.array([[n['x'], n['y']] for n in nodes.values() if n.get('node_type') == 'midmarker'], float)
    spacing = float(np.median(cKDTree(mids).query(mids, k=2)[0][:, 1])) if len(mids) > 2 else 60.0
    # ---- worst local regions (400 px windows, stride 200)
    allx = centres[:, 0]; ally = centres[:, 1]
    windows = []
    if len(centres):
        x0, y0 = float(allx.min()), float(ally.min())
        events = [(c['x'], c['y'], 'crossing', c['reactions']) for c in cross_points if c['kind'] != 'expert_expert'] + \
                 [(h['x'], h['y'], 'obstacle', [h['reaction']]) for h in hit_points] + \
                 [(float(centres[i][0]), float(centres[i][1]), 'node_overlap', []) for i, _ in node_overlap]
        mid_pts = [(float(np.mean(e['pts'][:, 0])), float(np.mean(e['pts'][:, 1]))) for e in edges]
        cells = _dd(lambda: dict(crossings=0, obstacles=0, node_overlaps=0, edges=0, reactions=set()))
        for x, y, what, rs in events:
            for ox in (0, 1):
                for oy in (0, 1):
                    key = (int((x - x0) // 400) - ox + 1, int((y - y0) // 400) - oy + 1)
                    c = cells[key]
                    c['crossings' if what == 'crossing' else 'obstacles' if what == 'obstacle' else 'node_overlaps'] += 1
                    c['reactions'].update(rs)
        for x, y in mid_pts:
            for ox in (0, 1):
                for oy in (0, 1):
                    cells[(int((x - x0) // 400) - ox + 1, int((y - y0) // 400) - oy + 1)]['edges'] += 1
        ranked = sorted(cells.items(), key=lambda kv: (-(kv[1]['crossings'] + kv[1]['obstacles'] + kv[1]['node_overlaps']), -kv[1]['edges'], kv[0]))
        taken = []
        for (gx, gy), c in ranked:
            cx, cy = x0 + gx * 400 - 200 + 200, y0 + gy * 400 - 200 + 200
            if c['crossings'] + c['obstacles'] + c['node_overlaps'] == 0 or any(abs(cx - tx) < 300 and abs(cy - ty) < 300 for tx, ty in taken):
                continue
            taken.append((cx, cy))
            windows.append(dict(x=round(cx - 200, 1), y=round(cy - 200, 1), w=400, h=400, crossings=c['crossings'], obstacle_hits=c['obstacles'], node_overlaps=c['node_overlaps'],
                                edges=c['edges'], reactions=sorted(c['reactions'])[:10]))
            if len(windows) == 5:
                break
    ori = orientation_module.summary(data, orientation_kinds)
    # ---- global composition: space use, module spacing anomalies, aspect, densest place (cells scale with the typical reaction spacing)
    mids_xy = mids if len(mids) else np.zeros((0, 2))
    comp = {}
    if len(mids_xy) > 8:
        cell = max(spacing * 4.0, 50.0)
        gx = ((mids_xy[:, 0] - mids_xy[:, 0].min()) // cell).astype(int)
        gy = ((mids_xy[:, 1] - mids_xy[:, 1].min()) // cell).astype(int)
        occupied = len(set(zip(gx.tolist(), gy.tolist())))
        total_cells = int((gx.max() + 1) * (gy.max() + 1))
        d5 = cKDTree(mids_xy).query(mids_xy, k=min(6, len(mids_xy)))[0][:, -1]
        med = float(np.median(d5))
        outliers = np.flatnonzero(d5 > 8.0 * max(med, 1.0))
        counts = _dd(int)
        for a_, b_ in zip(gx.tolist(), gy.tolist()):
            counts[a_, b_] += 1
        W_, H_ = (mids_xy.max(0) - mids_xy.min(0))
        comp = dict(canvas=[round(float(W_)), round(float(H_))], aspect=round(float(W_ / max(H_, 1.0)), 2), space_utilisation=round(occupied / max(total_cells, 1), 4),
                    cell_px=round(cell), neighbour5_median_px=round(med), neighbour5_p90_px=round(float(np.percentile(d5, 90))), neighbour5_max_px=round(float(d5.max())),
                    isolated_reactions=int(len(outliers)), isolation_threshold_px=round(8.0 * med), max_reactions_in_one_cell=int(max(counts.values())))
        outlier_rows = [dict(x=round(float(mids_xy[i][0]), 1), y=round(float(mids_xy[i][1]), 1), neighbour5_px=round(float(d5[i]))) for i in outliers[:200]]
    else:
        outlier_rows = []
    n_alg = max(len(algo), 1)
    n_all = max(len(edges), 1)
    total_cross = sum(cross.values())
    summary = dict(
        edges=dict(total=len(edges), algorithm=len(algo), expert=len(edges) - len(algo)),
        core=dict(
            crossings=dict(total=total_cross, same_reaction=cross['same_reaction'], same_reaction_expert=cross['same_reaction_expert'], algorithm_algorithm=cross['algorithm_algorithm'], algorithm_expert=cross['algorithm_expert'], expert_expert=cross['expert_expert'],
                           per_100_edges=round(100.0 * (cross['algorithm_algorithm'] + cross['algorithm_expert']) / n_alg, 2)),
            collisions=dict(node_overlap_pairs=len(node_overlap), label_overlap_pairs=label_overlap, curves_through_nodes=node_hits, curves_through_labels=label_hits),
            length_p90=dict(px=round(float(np.percentile(all_len, 90)), 1), in_reaction_spacings=round(float(np.percentile(all_len, 90)) / max(spacing, 1.0), 2), reaction_spacing_px=round(spacing, 1)),
            core_orientation=dict(opposite_sides_share=(ori.get('core_participants') or {}).get('opposite_sides_share'), reactions=(ori.get('core_participants') or {}).get('reactions')),
            structure_continuity=None,
            bezier_coverage=dict(share=round(bez_edges / n_alg, 3), algorithm_edges=len(algo), bezier=bez_edges, straight_lines=straight_edges, polylines_with_bends=poly_edges,
                                 orthogonal_polylines=ortho_edges)),
        details=dict(
            length=dict(p50=round(float(np.percentile(all_len, 50)), 1), p90=round(float(np.percentile(all_len, 90)), 1), p99=round(float(np.percentile(all_len, 99)), 1),
                        over_250=int((all_len > 250).sum()), over_1000=int((all_len > 1000).sum())),
            detour_ratio=dict(median=round(float(np.median(detours)), 3) if detours else None, p90=round(float(np.percentile(detours, 90)), 3) if detours else None,
                              over_1_6=int(sum(d > 1.6 for d in detours))),
            curve_quality=dict(hooks=hooks, sharp_turns=sharp, self_intersections=selfx),
            crossing_density=dict(worst_window_crossings=windows[0]['crossings'] if windows else 0),
            orientation=ori, composition=comp),
        worst_regions=windows,
        note='computed from the final map.json; crossings are counted on sampled final curves, source-split; structure continuity needs the V7 backbone layout and is not available yet')
    issues = dict(isolated_reactions=outlier_rows, crossings=cross_points, obstacle_hits=hit_points, curve_anomalies=anomalies, worst_regions=windows,
                  node_overlaps=[dict(x=round(float(centres[i][0]), 1), y=round(float(centres[i][1]), 1)) for i, _ in node_overlap][:max_issues])
    return summary, issues


def _turn_stats(pts):
    v = np.diff(pts, axis=0)
    n = np.hypot(v[:, 0], v[:, 1])
    ok = n > 1e-9
    v = v[ok] / n[ok][:, None]
    if len(v) < 2:
        return 0.0, 0.0
    ang = np.degrees(np.arccos(np.clip((v[:-1] * v[1:]).sum(1), -1, 1)))
    cross = v[:-1, 0] * v[1:, 1] - v[:-1, 1] * v[1:, 0]
    return float(ang.max()), float(abs(np.where(cross >= 0, ang, -ang).sum()))


def _self_cross(pts):
    for i in range(len(pts) - 1):
        for j in range(i + 2, len(pts) - 1):
            a, b, c, d = pts[i], pts[i + 1], pts[j], pts[j + 1]
            r, s = b - a, d - c
            den = r[0] * s[1] - r[1] * s[0]
            if abs(den) < 1e-9:
                continue
            t = ((c[0] - a[0]) * s[1] - (c[1] - a[1]) * s[0]) / den
            u = ((c[0] - a[0]) * r[1] - (c[1] - a[1]) * r[0]) / den
            if 1e-3 < t < 1 - 1e-3 and 1e-3 < u < 1 - 1e-3:
                return True
    return False


def _orthogonal(pts):
    """Straight pieces that turn by ~90 degrees at every corner (a right-angle polyline)."""
    v = np.diff(pts, axis=0)
    n = np.hypot(v[:, 0], v[:, 1])
    v = v[n > 1e-6] / n[n > 1e-6][:, None]
    corners = 0
    for a, b in zip(v[:-1], v[1:]):
        ang = float(np.degrees(np.arccos(np.clip(a @ b, -1, 1))))
        if ang > 20:
            if abs(ang - 90) > 15:
                return False
            corners += 1
    return corners >= 1


def metric_details(path):
    """Actual exported glyphs/labels, excluding routing control points (moved from geometry_experiment)."""
    import drawing
    from core import read
    data = read(path)
    content = data[1]
    nodes = {k: v for k, v in content['nodes'].items() if v['node_type'] != 'multimarker'}
    xy = {k: [v['x'], v['y']] for k, v in nodes.items()}
    boxes = {}
    for key, v in nodes.items():
        if v['node_type'] == 'metabolite':
            label = v['bigg_id']
            x = v['label_x']
            y = v['label_y']
        else:
            rx = next(r for r in content['reactions'].values() if any(s['from_node_id'] == key or s['to_node_id'] == key for s in r['segments'].values()))
            label = rx['bigg_id']
            x = rx['label_x']
            y = rx['label_y']
        b = drawing.FONT.getbbox(label)
        boxes[key] = [x, y - (b[3] - b[1]), x + drawing.FONT.getlength(label), y]
    node_label = sum(b[0] < v[0] + 4 and b[2] > v[0] - 4 and b[1] < v[1] + 4 and b[3] > v[1] - 4
                     for b in boxes.values() for v in xy.values())
    canvas = content['canvas']
    values = np.array(list(xy.values()))
    segments = []
    for r in content['reactions'].values():
        segments.extend((s['from_node_id'], s['to_node_id']) for s in r['segments'].values())
    intersections = 0
    for a, b in segments:
        av = np.array([content['nodes'][a]['x'], content['nodes'][a]['y']])
        bv = np.array([content['nodes'][b]['x'], content['nodes'][b]['y']])
        delta = bv - av
        for k, box in boxes.items():
            if k in (a, b):
                continue
            lower = 0.
            upper = 1.
            for axis in (0, 1):
                if abs(delta[axis]) < 1e-9:
                    if av[axis] <= box[axis] or av[axis] >= box[axis + 2]:
                        upper = -1
                        break
                else:
                    lo = (box[axis] - av[axis]) / delta[axis]
                    hi = (box[axis + 2] - av[axis]) / delta[axis]
                    lower = max(lower, min(lo, hi))
                    upper = min(upper, max(lo, hi))
            intersections += int(lower < upper - 1e-8)
    return dict(visible_nodes=len(nodes), label_node_intersections=int(node_label),
                segment_label_intersections=intersections,
                canvas_width=canvas['width'], canvas_height=canvas['height'],
                canvas_area_per_reaction=canvas['width'] * canvas['height'] / len(content['reactions']),
                node_coordinate_extent=np.ptp(values, axis=0).tolist())


def assert_chemistry(model, path):
    """Validate that the exported map matches the model chemistry (moved from geometry_experiment)."""
    from core import read, uid
    data = read(path)
    lookup = {r['bigg_id']: r for r in data[1]['reactions'].values()}
    assert set(lookup) == {r['id'] for r in model['reactions']}
    for r in model['reactions']:
        sign = -1 if r['upper_bound'] <= 0 and r['lower_bound'] < 0 else 1
        assert {m['bigg_id']: m['coefficient'] for m in lookup[r['id']]['metabolites']} == {k: sign * v for k, v in r['metabolites'].items()}
        assert data[0]['metadata']['source_reactions'][r['id']]['stoichiometry'] == r['metabolites']
    assert data[0]['metadata']['model_hash'] == uid(model)
