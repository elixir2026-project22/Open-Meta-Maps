"""Reaction connection-style templates, extracted from the expert Escher maps (Atlas extension `reaction_geometry_templates`).

The 22 Escher maps in `data/source_inputs/raw/community/` (CC0 / CC BY, as recorded in the source registry) draw every reaction as lines and curves between
the reaction's mid-marker and the metabolite glyphs.  This module keeps, for every reaction and every participant, only the CONSTRAINTS of that drawing,
never the curve itself, so that they can be re-solved for the real end points of another map:

  frame        origin at the reaction's mid-marker, +x along the reaction axis (mean substrate position -> mean product position)
  per participant:
    role           substrate / product (from the stoichiometric coefficient)    cur  currency / cofactor (core.CURRENCY) or not
    start_rel      angle between the curve's initial tangent and the straight chord to the glyph (degrees, signed so that + means "bows to the participant's own side")
    end_rel        the same at the glyph end (travel direction at the glyph versus the chord)
    h1, h2         first / second Bezier handle length as a share of the chord
    stub           length of the straight stub that leaves the reaction centre, as a share of the chord (0 = the curve starts bending right away)
    dev            largest deviation of the whole chain from the chord, as a share of the chord
  per reaction:    counts (core substrates, core products, currency substrates, currency products), transport flag, source map and reaction id

Index key `"s{core_sub},{core_prod},{cur_sub},{cur_prod}"`; `recall(signature)` returns the templates of the closest signature (exact first, then the nearest
counts).  What is transferred is geometry and visual structure; the chemistry (participants, stoichiometry) of the target GEM never changes.

  python3 reaction_templates.py build
"""
import gzip
import hashlib
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from core import BASE, CURRENCY

RAW = BASE / 'data/source_inputs/raw/community'
EXT = BASE / 'data/extensions/reaction_geometry_templates'
FILE = EXT / 'templates.json.gz'
_cache = {}


def _is_currency(bigg_id):
    return re.sub(r'_[a-z][0-9]?$', '', bigg_id) in CURRENCY


def _ang(v):
    return math.degrees(math.atan2(v[1], v[0]))


def _wrap(a):
    return (a + 180.0) % 360.0 - 180.0


def _bez(p0, c1, c2, p3, n=10):
    t = np.linspace(0, 1, n)[:, None]
    return (1 - t) ** 3 * p0 + 3 * (1 - t) ** 2 * t * c1 + 3 * (1 - t) * t ** 2 * c2 + t ** 3 * p3


def _reaction(rx, nodes):
    """Descriptors of one Escher reaction or None."""
    adj = defaultdict(list)
    mid = None
    for seg in rx['segments'].values():
        a, b = seg['from_node_id'], seg['to_node_id']
        adj[a].append((b, seg, True))
        adj[b].append((a, seg, False))
        for k in (a, b):
            if nodes[k]['node_type'] == 'midmarker':
                mid = k
    if mid is None:
        return None
    coeff = {m['bigg_id']: m['coefficient'] for m in rx['metabolites']}
    P0 = np.array([nodes[mid]['x'], nodes[mid]['y']])
    parts = []
    for target, n in nodes.items():
        if n.get('node_type') != 'metabolite' or n.get('bigg_id') not in coeff:
            continue
        # unique chain mid -> target
        stack, found = [(mid, [])], None
        seen = {mid}
        while stack and found is None:
            node, path = stack.pop()
            if node == target:
                found = path
                break
            for nxt, seg, fwd in adj.get(node, []):
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append((nxt, path + [(seg, fwd, node, nxt)]))
        if not found:
            continue
        pts, first_dir, last_dir, h1, h2, stub = [P0], None, None, None, None, 0.0
        for i, (seg, fwd, a, b) in enumerate(found):
            pa, pb = np.array([nodes[a]['x'], nodes[a]['y']]), np.array([nodes[b]['x'], nodes[b]['y']])
            c1 = seg.get('b1') if fwd else seg.get('b2')
            c2 = seg.get('b2') if fwd else seg.get('b1')
            seg_len = float(np.hypot(*(pb - pa)))
            if c1 and c2:
                q1, q2 = np.array([c1['x'], c1['y']]), np.array([c2['x'], c2['y']])
                curve = _bez(pa, q1, q2, pb)
                d_start, d_end = q1 - pa, pb - q2
                if seg_len > 1e-6:  # handles of the last Bezier segment of the chain, as a share of that segment
                    h1, h2 = float(np.hypot(*d_start)) / seg_len, float(np.hypot(*d_end)) / seg_len
                pts.extend(curve[1:].tolist())
            else:
                d_start = d_end = pb - pa
                if i == 0 and len(found) > 1:
                    stub = seg_len
                pts.append(pb.tolist())
            if first_dir is None and np.hypot(*d_start) > 1e-6:
                first_dir = d_start
            last_dir = d_end if np.hypot(*d_end) > 1e-6 else last_dir
        Q = np.array([n['x'], n['y']])
        chord = Q - P0
        L = float(np.hypot(*chord))
        if L < 25 or first_dir is None or last_dir is None:
            continue
        pts = np.array(pts)
        dev = float(np.max(np.abs((pts[:, 0] - P0[0]) * chord[1] - (pts[:, 1] - P0[1]) * chord[0]) / L) / L)
        parts.append(dict(bigg=n['bigg_id'], coeff=float(coeff[n['bigg_id']]), cur=_is_currency(n['bigg_id']), chord=chord, L=L,
                          start=_wrap(_ang(first_dir) - _ang(chord)), end=_wrap(_ang(last_dir) - _ang(chord)),
                          h1=h1, h2=h2, stub=stub / L, dev=dev))
    if len(parts) < 2:
        return None
    subs = [p for p in parts if p['coeff'] < 0]
    prods = [p for p in parts if p['coeff'] > 0]
    if not subs or not prods:
        return None
    u = np.mean([p['chord'] for p in prods], axis=0) - np.mean([p['chord'] for p in subs], axis=0)
    if np.hypot(*u) < 1e-6:
        return None
    u = u / np.hypot(*u)
    out = []
    for p in parts:
        side = float(np.sign(u[0] * p['chord'][1] - u[1] * p['chord'][0])) or 1.0  # which side of the axis the glyph lies on
        out.append(dict(role='substrate' if p['coeff'] < 0 else 'product', cur=bool(p['cur']), side=int(side),
                        chord_angle=round(_wrap(_ang(p['chord']) - _ang(u)), 1), L=round(p['L'], 1),
                        start_rel=round(p['start'] * side, 1), end_rel=round(p['end'] * side, 1),
                        h1=None if p['h1'] is None else round(p['h1'], 3), h2=None if p['h2'] is None else round(p['h2'], 3),
                        stub=round(p['stub'], 3), dev=round(p['dev'], 3)))
    chems = defaultdict(set)
    for p in parts:
        chems[re.sub(r'_[a-z][0-9]?$', '', p['bigg'])].add(p['bigg'])
    transport = any(len(v) > 1 for v in chems.values())
    counts = (sum(1 for p in out if p['role'] == 'substrate' and not p['cur']), sum(1 for p in out if p['role'] == 'product' and not p['cur']),
              sum(1 for p in out if p['role'] == 'substrate' and p['cur']), sum(1 for p in out if p['role'] == 'product' and p['cur']))
    return dict(counts=counts, transport=transport, parts=out)


def build():
    files = sorted(RAW.glob('*/*.json'))
    templates, digest = [], hashlib.sha256()
    for f in files:
        raw = f.read_bytes()
        digest.update(raw)
        data = json.loads(raw)
        nodes = data[1]['nodes']
        for rx in data[1]['reactions'].values():
            t = _reaction(rx, nodes)
            if t:
                t.update(map=f.stem, reaction=rx['bigg_id'])
                templates.append(t)
    index = defaultdict(list)
    for i, t in enumerate(templates):
        index['s%d,%d,%d,%d' % tuple(t['counts'])].append(i)
    EXT.mkdir(parents=True, exist_ok=True)
    with gzip.open(FILE, 'wt', encoding='utf-8') as fh:
        json.dump(dict(version=1, templates=templates, index=index), fh)
    meta = dict(builder='reaction_templates.py', maps=len(files), templates=len(templates), signatures=len(index), raw_files_sha256=digest.hexdigest(),
                licence='Expert Escher community maps (CC0 / CC BY, see data/source_inputs/source_registry.json); only drawing constraints are stored, not curves or chemistry',
                stored='per participant: role, currency flag, side, chord angle, start/end tangent angle relative to the chord, handle length ratios, stub ratio, deviation ratio')
    (EXT / 'SOURCE.json').write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta, indent=1))


def load():
    if 'data' not in _cache:
        if not FILE.exists():
            _cache['data'] = None
        else:
            with gzip.open(FILE, 'rt', encoding='utf-8') as fh:
                _cache['data'] = json.load(fh)
    return _cache['data']


def recall(counts, transport=None, k=60):
    """Templates whose signature is closest to `counts` = (core subs, core prods, currency subs, currency prods); exact bucket first."""
    d = load()
    if d is None:
        return []
    key = 's%d,%d,%d,%d' % tuple(counts)
    if d['index'].get(key):
        pool = d['index'][key]
    else:
        def dist(s):
            c = [int(x) for x in s[1:].split(',')]
            return sum(abs(a - b) * (3 if i < 2 else 1) for i, (a, b) in enumerate(zip(c, counts)))
        best = sorted(d['index'], key=lambda s: (dist(s), -len(d['index'][s])))[:3]
        pool = [i for s in best for i in d['index'][s]]
    ts = [d['templates'][i] for i in pool]
    if transport is not None:
        same = [t for t in ts if t['transport'] == transport]
        ts = same or ts
    return ts[:max(k, 1) * 5]


def stats(counts, transport, role, cur, rank_side_same):
    """Median constraints for a participant of the given role / currency flag from the nearest templates; None if there are no templates."""
    ts = recall(counts, transport)
    rows = [p for t in ts for p in t['parts'] if p['role'] == role and p['cur'] == cur and p['h1'] is not None and p['h2'] is not None]
    if len(rows) < 3:
        rows = [p for t in ts for p in t['parts'] if p['role'] == role and p['h1'] is not None and p['h2'] is not None]
    if len(rows) < 3:
        return None
    med = lambda key: float(np.median([p[key] for p in rows]))
    return dict(start_rel=med('start_rel'), end_rel=med('end_rel'), h1=med('h1'), h2=med('h2'), dev=med('dev'), n=len(rows),
                source=sorted({(t['map'], t['reaction']) for t in ts})[:3])


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'build':
        build()
    else:
        d = load()
        print(None if d is None else (len(d['templates']), len(d['index'])))
