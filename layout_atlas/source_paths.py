"""Original expert connection geometry ("source paths") for the Atlas modules.

The Atlas module records keep node positions only.  The expert maps they come from also drew every reaction as lines and curves
(Escher: chains of cubic Bezier segments between a metabolite, optional multimarkers and the reaction's midmarker).  That geometry
is still in the stored source observations (data/source_inputs/map_observations.json.gz); this module extracts it, keeps it as a
separate queryable data asset (data/extensions/source_paths/paths.json.gz) and resolves it for the modules the Composer selects.

Identity layers are kept apart: a path is addressed by VISUAL instance ids (source map -> reaction instance key -> metabolite glyph id,
with the segment keys of the chain) and carries the biochemical identity (reaction id, metabolite entity) separately.  Chains are
assembled from the original segment connectivity (from_node_id / to_node_id), never from coordinate proximity, and always oriented
from the reaction's midmarker to the metabolite glyph.  Escher segment geometry: cubic Bezier from the `from` node through b1, b2 to the
`to` node (b1/b2 null = straight); a reversed traversal swaps b1 and b2.

Use:
  python3 source_paths.py build        # extract -> data/extensions/source_paths/
  python3 source_paths.py show <map_id|name part> [reaction id]
  paths_for_block(block)               # per selected structure: chains in the FINAL frame, relative to the reaction's midpoint
Only Escher ('escher') maps are extracted so far; WikiPathways / Reactome follow the same layout of the data file.
"""
import gzip
import sys

from core import *

EXT = DATA / 'extensions/source_paths'
FILE = EXT / 'paths.json.gz'
_cache = {}


# ---------------------------------------------------------------------------------------------- extraction (Escher)
def _chain(mid, target, adjacency):
    """Unique path midmarker -> target over the reaction's segments; returns [(segment, forward)] or None."""
    stack, seen = [(mid, [])], {mid}
    found = []
    while stack:
        node, path = stack.pop()
        if node == target:
            found.append(path)
            continue
        for nxt, seg, forward in adjacency.get(node, []):
            if nxt not in seen:
                seen.add(nxt)
                stack.append((nxt, path + [(seg, forward)]))
    return found[0] if len(found) == 1 else None


def extract_map(m):
    """One Escher source map observation -> {glyphs, reactions:{key: {...}}} with original coordinates."""
    glyphs = m['glyphs'] if isinstance(m['glyphs'], dict) else {g['key']: g for g in m['glyphs']}
    out = dict(name=m['name'], source=m['source'], format=m['format'], source_sha256=m.get('source_sha256'), reactions={})
    used = set()
    stats = Counter()
    for r in m['reactions']:
        gl = {k: glyphs[k] for k in r.get('glyphs', []) if k in glyphs}
        mids = [k for k, g in gl.items() if g['kind'] == 'midmarker']
        if len(mids) != 1:
            stats['reaction_without_single_midmarker'] += 1
            continue
        mid = mids[0]
        adjacency = defaultdict(list)
        for p in r.get('paths', []):
            a, b = p['from_node_id'], p['to_node_id']
            adjacency[a].append((b, p['key'], True))
            adjacency[b].append((a, p['key'], False))
        segs = {p['key']: p for p in r.get('paths', [])}
        stats['segments'] += len(segs)  # unique source segments (chains share the middle segments)
        stats['bezier_segments'] += sum(bool(p.get('b1') and p.get('b2')) for p in segs.values())
        chains = {}
        for gid, g in gl.items():
            if g['kind'] != 'metabolite':
                continue
            path = _chain(mid, gid, adjacency)
            if path is None:
                stats['chain_missing_or_ambiguous'] += 1
                continue
            steps, prev = [], mid
            for key, forward in path:
                p = segs[key]
                a, b = (p['from_node_id'], p['to_node_id']) if forward else (p['to_node_id'], p['from_node_id'])
                c1, c2 = (p.get('b1'), p.get('b2')) if forward else (p.get('b2'), p.get('b1'))
                steps.append(dict(segment=key, **{'from': a, 'to': b}, controls=[[c1['x'], c1['y']], [c2['x'], c2['y']]] if c1 and c2 else None))
                used.update((a, b))
            chains[gid] = dict(entity=g.get('entity'), steps=steps)
            stats['chains'] += 1
        out['reactions'][r['key']] = dict(id=r['id'], name=r.get('name'), mid=mid, xy=r['xy'], chains=chains)
    out['glyphs'] = {k: dict(xy=glyphs[k]['xy'], kind=glyphs[k]['kind'], entity=glyphs[k].get('entity')) for k in used | {c for r in out['reactions'].values() for c in r['chains']} | {r['mid'] for r in out['reactions'].values()} if k in glyphs}
    out['stats'] = dict(stats)
    return out


def _half_point(pts):
    """Point at half of the polyline's arc length (same rule database.wp_diagram uses for the reaction centre) and its segment index."""
    lens = [math.hypot(b[0] - a[0], b[1] - a[1]) for a, b in zip(pts, pts[1:])]
    half, acc = sum(lens) / 2, 0.0
    for i, v in enumerate(lens):
        if acc + v >= half:
            t = (half - acc) / max(v, 1e-8)
            return [pts[i][0] + (pts[i + 1][0] - pts[i][0]) * t, pts[i][1] + (pts[i + 1][1] - pts[i][1]) * t], i
        acc += v
    return list(pts[0]), 0


def extract_gpml_map(m):
    """WikiPathways (GPML) map: every glyph centre plus the edges that carry explicit geometry (3 or more points).

    A GPML interaction is one line from a source glyph to a target glyph; the Atlas reaction node sits at the middle of its length.  Edges
    with two points and no connector type are plain straight lines: nothing to store (the Composer's straight edge is exact).  Connectors
    typed Elbow/Curved with only two points store their end points and anchors but NOT their shape (PathVisio computes it), so they are
    not extracted either.  `exact` is True for polylines whose every vertex is stored (Segmented, or untyped); Elbow/Curved connectors
    with waypoints keep the waypoints only (their smoothing / orthogonal routing between them is not reproduced).
    """
    glyphs = {g['key']: dict(xy=g['xy'], kind=g['kind'], width=g.get('width'), height=g.get('height')) for g in (m['glyphs'] if isinstance(m['glyphs'], list) else m['glyphs'].values())}
    out = dict(name=m['name'], source=m['source'], format='gpml', source_sha256=m.get('source_sha256'), glyphs=glyphs, reactions={})
    stats = Counter()
    for e in m['edges']:
        pts = [[float(p['X']), float(p['Y'])] for p in e.get('points', [])]
        connector = (e.get('attributes') or {}).get('ConnectorType') or 'Straight'
        stats['edges'] += 1
        stats['connector_' + connector] += 1
        if len(pts) < 3:
            stats['two_point_edges'] += 1
            continue
        if e.get('source') not in glyphs or e.get('target') not in glyphs:
            stats['edge_end_glyph_missing'] += 1
            continue
        mid, i = _half_point(pts)
        to_source = [mid] + pts[:i + 1][::-1]
        to_target = [mid] + pts[i + 1:]
        out['reactions'][e['key']] = dict(mid=mid, connector=connector, exact=connector in ('Segmented', 'Straight'), arrow=(e['points'][-1].get('ArrowHead')), source=e['source'], target=e['target'],
                                          chains={e['source']: dict(polyline=to_source), e['target']: dict(polyline=to_target)})
        stats['explicit_edges'] += 1
        stats['explicit_exact'] += connector in ('Segmented', 'Straight')
    out['stats'] = dict(stats)
    return out


def build():
    obs = json.loads(gzip.open(INPUTS / 'map_observations.json.gz', 'rt').read())
    maps, total = {}, Counter()
    totals_by_format = defaultdict(Counter)
    for m in obs:
        if m['source'] == 'escher' and m['format'] == 'escher':
            e = extract_map(m)
        elif m['source'] == 'wikipathways' and m['format'] == 'gpml':
            e = extract_gpml_map(m)
        else:
            continue
        maps[m['id']] = e
        totals_by_format[e['format']].update(e['stats'])
        totals_by_format[e['format']]['reactions'] += len(e['reactions'])
        if e['format'] == 'escher':
            total.update(e['stats'])
            total['reactions'] += len(e['reactions'])
    EXT.mkdir(parents=True, exist_ok=True)
    with gzip.open(FILE, 'wt') as f:
        json.dump(dict(version=1, maps=maps), f, ensure_ascii=False, separators=(',', ':'))
    summary = dict(extension='source paths (original expert connection geometry)', formats=['escher', 'gpml'], maps=len(maps), totals=dict(total), totals_by_format={k: dict(v) for k, v in totals_by_format.items()},
                   bezier_share=round(total['bezier_segments'] / max(total['segments'], 1), 4),
                   input=dict(file='data/source_inputs/map_observations.json.gz', sha256=sha((INPUTS / 'map_observations.json.gz').read_bytes())),
                   notes=['Absolute source coordinates; chains run from the midmarker to the metabolite glyph.',
                          'Modules resolve paths through their own node keys (r:<reaction key>, g:<glyph id>) and map_id; no per-module copy is stored.',
                          'WikiPathways: only edges with >=3 stored points are kept (two-point edges are straight lines or typed connectors whose shape PathVisio computes and the file does not store). Reactome (segments) is not extracted yet.'],
                   code_sha256=sha(Path(__file__).read_bytes()))
    dump(EXT / 'SOURCE.json', summary)
    return summary


def load():
    if 'data' not in _cache:
        _cache['data'] = json.loads(gzip.open(FILE, 'rt').read()) if FILE.exists() else dict(version=0, maps={})
    return _cache['data']


# ---------------------------------------------------------------------------------------------- resolve for a selected structure
def similarity(src, dst):
    """Least-squares similarity (scale, rotation, translation; reflection allowed and reported) src -> dst."""
    src, dst = np.asarray(src, float), np.asarray(dst, float)
    best = None
    for reflect in (False, True):
        s = src * np.array([1, -1]) if reflect else src
        mu_s, mu_d = s.mean(0), dst.mean(0)
        a, b = s - mu_s, dst - mu_d
        var = (a ** 2).sum()
        if var < 1e-9:
            continue
        u, d, vt = np.linalg.svd(b.T @ a)
        sign = np.sign(np.linalg.det(u @ vt)) or 1.0
        rot = u @ np.diag([1, sign]) @ vt
        scale = (d * [1, sign]).sum() / var
        t = mu_d - scale * (rot @ mu_s)
        pred = scale * (s @ rot.T) + t
        res = np.hypot(*(pred - dst).T)
        if best is None or res.max() < best['max']:
            best = dict(scale=float(scale), rot=rot, t=t, reflect=reflect, max=float(res.max()), rms=float(np.sqrt((res ** 2).mean())))
    return best


def _apply(fit, p):
    p = np.asarray(p, float)
    if fit['reflect']:
        p = p * np.array([1, -1])
    return fit['scale'] * (fit['rot'] @ p) + fit['t']


def paths_for_block(block, tolerance=3.0):
    """Source paths of one selected structure in the final frame, relative to each reaction's own midpoint.

    Returns (paths, info).  paths[(target reaction id, target metabolite id)] = dict(status='preserved', segs=[['L', dx, dy] | ['C', ...]],
    source={map, reaction, glyph, segments}); `info` counts what was eligible and why a chain was not reused.
    """
    info = Counter()
    data = load()['maps'].get(block.get('map_id'))
    if not data:
        info['no_source_paths_for_map'] += 1
        return {}, dict(info)
    if data.get('format') == 'gpml':
        return _paths_for_gpml_block(block, data, tolerance)
    corr_src, corr_dst, keyed = [], [], {}
    mapping = block['node_mapping']
    targets = Counter(mapping.values())
    for tkey, target in mapping.items():
        if targets[target] != 1 or target not in block['positions']:
            continue
        if tkey.startswith('r:'):
            r = data['reactions'].get(tkey[2:])
            xy = r['xy'] if r else None
        elif tkey.startswith('g:'):
            g = data['glyphs'].get(tkey[2:])
            xy = g['xy'] if g else None
        else:
            xy = None
        if xy is not None:
            corr_src.append(xy)
            corr_dst.append(block['positions'][target])
    if len(corr_src) < 3:
        info['too_few_correspondences'] += 1
        return {}, dict(info)
    fit = similarity(corr_src, corr_dst)
    if fit is None or fit['max'] > tolerance:
        info['transform_not_similarity'] += 1
        info['fit_max_residual_px'] = round(fit['max'], 2) if fit else None
        return {}, dict(info)
    info['fit_max_residual_px'] = round(fit['max'], 3)
    info['reflected'] = int(fit['reflect'])
    out = {}
    for tkey, target_r in block['reaction_mapping'].items():
        target_r = target_r if str(target_r).startswith('r:') else 'r:' + str(target_r)  # reaction_mapping stores plain reaction ids
        r = data['reactions'].get(tkey[2:] if tkey.startswith('r:') else tkey)
        if not r:
            info['reaction_not_in_source'] += 1
            continue
        mid = _apply(fit, data['glyphs'][r['mid']]['xy'])
        for gid, chain in r['chains'].items():
            tm = mapping.get('g:' + gid)
            if tm is None:
                info['glyph_not_in_module'] += 1  # e.g. currency metabolites: drawn by the Composer, not part of the module
                continue
            info['chains'] += 1
            segs, ok = [], True
            for st in chain['steps']:
                end = _apply(fit, data['glyphs'][st['to']]['xy']) - mid
                if st['controls']:
                    c1 = _apply(fit, st['controls'][0]) - mid
                    c2 = _apply(fit, st['controls'][1]) - mid
                    segs.append(['C', *c1.round(2).tolist(), *c2.round(2).tolist(), *end.round(2).tolist()])
                else:
                    segs.append(['L', *end.round(2).tolist()])
            # the chain must finish on the metabolite glyph that the structure drew (same target glyph, same place)
            gx = np.array(block['positions'][tm]) - np.array(block['positions'][target_r])
            deviation = float(np.hypot(*(np.array(segs[-1][-2:]) - gx)))
            if deviation > tolerance:
                info['endpoint_mismatch'] += 1
                continue
            info['preserved'] += 1
            out[(target_r.split(':', 1)[1], tm.split(':', 1)[1])] = dict(status='preserved', segs=segs, deviation=round(deviation, 3),
                                                                      source=dict(map=block['map_id'], reaction=tkey[2:], glyph=gid, segments=[s['segment'] for s in chain['steps']]))
    return out, dict(info)


def _paths_for_gpml_block(block, data, tolerance):
    """WikiPathways structure: polylines from the reaction's midpoint (middle of the original interaction line) to each glyph centre."""
    info = Counter()
    mapping = block['node_mapping']
    targets = Counter(mapping.values())
    src, dst = [], []
    for tkey, target in mapping.items():
        if targets[target] != 1 or target not in block['positions']:
            continue
        if tkey.startswith('r:'):
            r = data['reactions'].get(tkey[2:])
            xy = r['mid'] if r else None
        else:
            g = data['glyphs'].get(tkey[2:] if tkey.startswith('g:') else tkey)
            xy = g['xy'] if g else None
        if xy is not None:
            src.append(xy)
            dst.append(block['positions'][target])
    if len(src) < 3:
        info['too_few_correspondences'] += 1
        return {}, dict(info)
    fit = similarity(src, dst)
    if fit is None or fit['max'] > tolerance:
        info['transform_not_similarity'] += 1
        return {}, dict(info)
    info['fit_max_residual_px'] = round(fit['max'], 3)
    out = {}
    for tkey, target_r in block['reaction_mapping'].items():
        target_r = target_r if str(target_r).startswith('r:') else 'r:' + str(target_r)
        r = data['reactions'].get(tkey[2:] if tkey.startswith('r:') else tkey)
        if r is None:
            info['straight_or_unstored_edge'] += 1  # drawn as the exact straight line by the Composer
            continue
        mid = _apply(fit, r['mid'])
        for gid, chain in r['chains'].items():
            tm = mapping.get(gid) or mapping.get('g:' + gid)
            if tm is None:
                info['glyph_not_in_module'] += 1
                continue
            info['chains'] += 1
            pts = [_apply(fit, p) - mid for p in chain['polyline'][1:]]
            center = _apply(fit, data['glyphs'][gid]['xy']) - mid
            if not pts or np.hypot(*(pts[-1] - center)) > 0.5:
                pts.append(center)  # the stored end point lies on the glyph's outline: continue to the glyph centre
            segs = [['L', *p.round(2).tolist()] for p in pts]
            gx = np.array(block['positions'][tm]) - np.array(block['positions'][target_r])
            deviation = float(np.hypot(*(np.array(segs[-1][-2:]) - gx)))
            if deviation > tolerance:
                info['endpoint_mismatch'] += 1
                continue
            info['preserved'] += 1
            info['exact_polyline' if r['exact'] else 'waypoints_only'] += 1
            out[(target_r.split(':', 1)[1], tm.split(':', 1)[1])] = dict(status='preserved', segs=segs, deviation=round(deviation, 3), geometry_status='original' if r['exact'] else 'waypoints_only',
                                                                      source=dict(map=block['map_id'], reaction=tkey[2:] if tkey.startswith('r:') else tkey, glyph=gid, connector=r['connector']))
    return out, dict(info)


def adjust_chain(segs, delta, max_shift=40.0, max_turn_deg=40.0):
    """Constrained local adjustment of a chain whose end glyph was moved by `delta` (dx, dy) after the structure was placed.

    Only the LAST segment changes: its end point and (for a Bezier) its second control point move with the glyph, so the end tangent and the
    start tangent keep their directions and the body of the curve stays.  Refused (None) when the glyph moved too far or the chord of the
    last segment would turn or stretch too much; the Composer then keeps its own routing for that edge.
    """
    dx, dy = float(delta[0]), float(delta[1])
    shift = math.hypot(dx, dy)
    if shift > max_shift:
        return None, 'shift_too_large'
    last = list(segs[-1])
    p0 = np.array(segs[-2][-2:], float) if len(segs) > 1 else np.zeros(2)
    old_end = np.array(last[-2:], float)
    new_end = old_end + [dx, dy]
    a, b = old_end - p0, new_end - p0
    la, lb = float(np.hypot(*a)), float(np.hypot(*b))
    if la < 1e-6 or lb < 1e-6 or not 0.5 <= lb / la <= 2.0:
        return None, 'segment_stretched'
    turn = math.degrees(math.acos(max(-1.0, min(1.0, float(a @ b) / (la * lb)))))
    if turn > max_turn_deg:
        return None, 'segment_turned'
    if last[0] == 'C':
        last = ['C', last[1], last[2], round(last[3] + dx, 2), round(last[4] + dy, 2), round(new_end[0], 2), round(new_end[1], 2)]
    else:
        last = ['L', round(new_end[0], 2), round(new_end[1], 2)]
    return segs[:-1] + [last], shift


def write_curves_to_escher(data, found):
    """Write the reused curves into the Escher JSON itself, as standard segments (b1/b2 control points), so a downloaded map.json keeps them.

    For every reused chain the routed segments from the reaction's midmarker to the metabolite node are replaced by the source chain: new
    `multimarker` nodes for the chain's inner points and one segment per source segment (cubic Bezier with b1/b2 or a straight segment with
    null controls), ending on the metabolite node the map drew.  Segment direction stays midmarker -> metabolite, as in the routed export.
    Returns the entries that were written (others keep their routed segments).
    """
    body = data[1]
    nodes, reactions = body['nodes'], body['reactions']
    by_id = {r['bigg_id']: r for r in reactions.values()}
    written, counter = [], 0
    for f in found:
        rx = by_id.get(f['reaction'])
        if rx is None:
            continue
        segs = rx['segments']
        incoming = {s['to_node_id']: k for k, s in segs.items()}
        first = next(iter(segs.values()))
        mid_id = first['from_node_id']
        mx, my = nodes[mid_id]['x'], nodes[mid_id]['y']
        ex, ey = mx + f['segs'][-1][-2], my + f['segs'][-1][-1]
        ends = [nid for nid, k in incoming.items() if nodes[nid].get('node_type') == 'metabolite' and nodes[nid].get('bigg_id') == f['metabolite']]
        if not ends:
            continue
        end_id = min(ends, key=lambda n: math.hypot(nodes[n]['x'] - ex, nodes[n]['y'] - ey))
        # walk back over the routed chain (midmarker -> bends -> metabolite) and drop it
        chain, cur = [], end_id
        while cur != mid_id:
            k = incoming.get(cur)
            if k is None:
                chain = None
                break
            chain.append(k)
            cur = segs[k]['from_node_id']
        if chain is None:
            continue
        bends = [segs[k]['from_node_id'] for k in chain if segs[k]['from_node_id'] != mid_id]
        for k in chain:
            del segs[k]
        for b in bends:
            if not any(b in (s['from_node_id'], s['to_node_id']) for r2 in reactions.values() for s in r2['segments'].values()):
                nodes.pop(b, None)
        prev = mid_id
        for i, seg in enumerate(f['segs']):
            last = i == len(f['segs']) - 1
            if last:
                to = end_id
            else:
                counter += 1
                to = f'sp{counter}'
                nodes[to] = dict(node_type='multimarker', x=round(mx + seg[-2], 2), y=round(my + seg[-1], 2))
            if seg[0] == 'C':
                b1, b2 = dict(x=round(mx + seg[1], 2), y=round(my + seg[2], 2)), dict(x=round(mx + seg[3], 2), y=round(my + seg[4], 2))
            else:
                b1 = b2 = None
            counter += 1
            segs[f'spseg{counter}'] = dict(from_node_id=prev, to_node_id=to, b1=b1, b2=b2)
            prev = to
        written.append(f)
    return written


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'build':
        print(json.dumps(build(), indent=1))
    elif len(sys.argv) > 2 and sys.argv[1] == 'show':
        for mid, m in load()['maps'].items():
            if sys.argv[2] in mid or sys.argv[2] in m['name']:
                print(mid, m['name'], m['stats'])
                if len(sys.argv) > 3:
                    for k, r in m['reactions'].items():
                        if r['id'] == sys.argv[3]:
                            print(json.dumps(r, indent=1)[:3000])
    else:
        print(__doc__)
