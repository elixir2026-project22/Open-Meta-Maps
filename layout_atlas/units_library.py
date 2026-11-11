"""Non-overlapping 'structure unit' library: an independent extension recut from the same source maps as the open core.

A unit is a connected set of 3-12 reactions that forms one drawing gestalt (cycle, fork or chain) on the *core* graph
of ONE source map (reaction/metabolite bipartite graph, currency metabolites already removed exactly as in
`database.diagram_modules` / the Escher core projection).  Every reaction of a map belongs to exactly one unit (or is
an appendix of a unit, or an orphan: an isolated <3-reaction island that is stored only as a drawing statistic).

Algorithm (deterministic, fixed ordering):
  1. cycles: short bipartite cycles (<=12 reactions) from a minimum cycle basis, taken greedily (shortest first) when all
     their reactions are still free;
  2. forks: a branch point (core metabolite shared by >=3 free reactions, or a reaction with >=3 free neighbours) with
     its arms (one reaction per arm, extended up to a size budget of 12);
  3. chains: what remains is cut at branch points; every unbranched path of >=3 reactions is a chain (chains >12 are
     split at the midpoint recursively);
  4. fragments (<3 reactions) are merged as appendices into an adjacent unit; fragments adjacent only to other fragments
     are unioned and kept when they reach 3 reactions.  Anything left is an orphan.
1- and 2-reaction pieces are never modules; `drawing_conventions.json` keeps orientation/angle/spacing statistics.

Build:  python3 units_library.py          (re-executes itself with PYTHONHASHSEED=0 for reproducibility)
"""
import gzip
import itertools
import sqlite3
import sys
import zlib

from core import *
import database
from semantics import rhea_index, parse_gpml

EXT = DATA / 'extensions/units'
RESULTS = OUT / 'units_extension'
CONFIG = dict(min_size=3, max_size=12, max_appendix_total=16, ring_rule=True, ring_min_reactions=6, fuse_min_shared=2, fuse_max_new=4, fork_max_arms=8, fork_arm_depth=3, cycle_max_reactions=12,
              variant_threshold=.25, mapping_cap=3000, step_budget=400000, max_module_nodes=100)
SMALL_CURRENCY_LABELS = {'h2o', 'water', 'h+', 'h(+)', 'atp', 'adp', 'amp', 'nad+', 'nadh', 'nadp+', 'nadph', 'coa-sh', 'coa', 'phosphate', 'pi', 'ppi'}
SMALL_CURRENCY_CHEBI = {'chebi:15377', 'chebi:15378', 'chebi:30616', 'chebi:456216', 'chebi:57287', 'chebi:57288', 'chebi:57945', 'chebi:57925'}


# ---------------------------------------------------------------- source maps -> core graphs

def escher_maps():
    """Reconstruct whole Escher maps as the union of their stored windows (absolute source-pixel coordinates)."""
    by = defaultdict(list)
    for m in sorted(source_modules('semantic_transfer'), key=lambda m: m['id']):
        if m['source'] == 'escher':
            by[m['map_id']].append(m)
    out = []
    for map_id in sorted(by):
        wins = by[map_id]
        glyphs, rx, edges, rxy, drift = {}, {}, {}, {}, 0.
        for m in wins:
            c, s = np.array(m['coordinate_frame']['center']), m['coordinate_frame']['scale']
            records = {str(r['key']): r for r in m['reactions']}
            ent = {n['key']: n.get('entity') for n in m['core_nodes']}
            for n in m['core_nodes']:
                xy = (np.array(n['xy']) * s + c)
                store = rxy if n['key'].startswith('r:') else glyphs
                if n['key'] in store:
                    drift = max(drift, float(np.abs(np.array(store[n['key']]['xy'] if isinstance(store[n['key']], dict) else store[n['key']]) - xy).max()))
                if n['key'].startswith('r:'):
                    rxy[n['key']] = xy.tolist()
                else:
                    glyphs[n['key']] = dict(xy=xy.tolist(), entity=n.get('entity'), label=n.get('label', ''))
            for a, b in m['core_edges']:
                rk = a if a.startswith('r:') else b
                mk = b if rk == a else a
                coef = records[rk[2:]]['coefficients'].get(ent[mk], 0)
                edges[(rk, mk)] = coef
            for r in m['reactions']:
                rx.setdefault('r:' + str(r['key']), dict(r, key='r:' + str(r['key']), references={'bigg.reaction': [r['id']]}))
        d = dict(id=map_id, source='escher', species=wins[0]['species'], lineage=wins[0]['lineage'], name=wins[0]['name'], license=wins[0]['license'],
                 source_sha256=wins[0]['source_sha256'], original_instance='data/source_inputs/module_observations.json.gz', window_count=len(wins),
                 max_coordinate_drift=drift, reactions=rx, rxy=rxy, glyphs=glyphs, edges=edges,
                 rorigin={k: 'observed' for k in rxy})
        out.append(d)
    return out


def _core_visible(r, glyphs):
    return [k for k in r['visible'] if glyphs[k].get('label', '').lower() not in SMALL_CURRENCY_LABELS and str(glyphs[k]['entity']).split('@')[0] not in SMALL_CURRENCY_CHEBI]


def diagram_map(d):
    """Whole-map core graph of a WikiPathways/Reactome diagram (same core filter as database.diagram_modules)."""
    rx, rxy, edges, glyphs = {}, {}, {}, {}
    for r in d['events']:
        core = _core_visible(r, d['glyphs'])
        if not core:
            continue
        rx[r['key']] = r
        rxy[r['key']] = list(r['xy'])
        for k in core:
            glyphs[k] = dict(xy=list(d['glyphs'][k]['xy']), entity=d['glyphs'][k]['entity'], label=d['glyphs'][k].get('label', ''))
            edges[(r['key'], k)] = r['visible'][k]
    return dict(id=d['id'], source=d['source'], species=d['species'], lineage=d['lineage'], name=d['name'], license=d['license'], source_sha256=d['source_sha256'],
                original_instance=d.get('path'), reactions=rx, rxy=rxy, glyphs=glyphs, edges=edges, rorigin={k: rx[k]['position_origin'] for k in rx})


def diagram_maps():
    rhea = rhea_index()
    original_maps = source_maps()
    bysha = {d['source_sha256']: d for d in original_maps if d['source'] == 'wikipathways'}
    out = []
    for p in sorted((INPUTS / 'raw/wp').rglob('*.gpml')):
        d = parse_gpml(p, rhea)
        if d['source_sha256'] in bysha:
            previous = bysha[d['source_sha256']]
            d['id'] = previous['id']
            d['lineage'] = previous['lineage']
        if not d['bindings']:
            continue
        out.append(diagram_map(database.wp_diagram(d, rhea)))
    acq = read(DATA / 'acquisition.json')
    for rec in acq['reactome_diagrams']:
        p = BASE / rec['layout_path'] if rec.get('layout_path') else DATA / 'reactome' / (rec['id'] + '.json')
        try:
            out.append(diagram_map(database.reactome_diagram(p)))
        except Exception:
            continue
    return out


def geometry_maps(exclude_ids):
    """Geometry-only WikiPathways maps (metabolite glyphs + interaction arrows, no reaction identity): segmented by SHAPE only.
    Each drawn interaction between two metabolite glyphs is a pseudo-reaction (a derived route waypoint, never an observed reaction centre)."""
    out = []
    for d in sorted(source_maps(), key=lambda d: d['id']):
        if d['source'] != 'wikipathways' or d.get('quality') != 'geometry_only' or d['id'] in exclude_ids:
            continue
        glyphs = {k: dict(xy=list(g['xy']), entity=g.get('entity'), label=g.get('label', '')) for k, g in sorted(d['glyphs'].items()) if g.get('kind') == 'Metabolite'}
        rx, rxy, edges = {}, {}, {}
        for e in sorted(d['edges'], key=lambda e: str(e['key'])):
            a, b = e.get('source'), e.get('target')
            if a in glyphs and b in glyphs and a != b:
                rk = 'route:' + str(e['key'])
                rx[rk] = dict(key=rk, id=rk)
                rxy[rk] = ((np.array(glyphs[a]['xy']) + np.array(glyphs[b]['xy'])) / 2).tolist()
                edges[(rk, a)] = -1; edges[(rk, b)] = 1
        if rx:
            out.append(dict(id=d['id'], source='wikipathways', species=d['species'], lineage=d['lineage'], name=d['name'], license=d['license'], source_sha256=d['source_sha256'],
                            original_instance='data/source_inputs/map_observations.json.gz', reactions=rx, rxy=rxy, glyphs=glyphs, edges=edges,
                            rorigin={k: 'derived_route_waypoint' for k in rx}, geometry_only=True))
    return out


def source_core_maps():
    semantic = escher_maps() + diagram_maps()
    geometry = geometry_maps({d['id'] for d in semantic})
    return sorted(semantic + geometry, key=lambda d: (d['source'], d['id']))


# ---------------------------------------------------------------- decomposition

def _bip(rkeys, r2m):
    g = nx.Graph()
    for r in sorted(rkeys):
        g.add_node(r, role='reaction')
        for m in sorted(r2m[r]):
            g.add_node(m, role='metabolite')
            g.add_edge(r, m)
    return g


def short_cycles(g, max_reactions):
    """Shortest cycle through every edge (bounded BFS, fixed neighbour order), 3..max_reactions reactions.

    This is the candidate set a minimum cycle basis is chosen from; `nx.minimum_cycle_basis` was measured at 268 s for one
    421-reaction fatty-acid map, so the bounded per-edge search is used instead (every basis cycle is a shortest cycle
    through at least one of its edges in these graphs, and the greedy shortest-first choice below is a basis-like selection).
    """
    limit = 2 * max_reactions - 1  # max path length in edges between the two ends of the removed edge
    out = set()
    for u, v in sorted(tuple(sorted(e)) for e in g.edges()):
        prev = {u: None}
        frontier = [u]
        depth = 0
        found = False
        while frontier and depth < limit and not found:
            depth += 1
            nxt = []
            for a in frontier:
                for b in sorted(g[a]):
                    if (a == u and b == v) or b in prev:
                        continue
                    prev[b] = a
                    if b == v:
                        found = True
                        break
                    nxt.append(b)
                if found:
                    break
            frontier = nxt
        if found:
            path = [v]
            while prev[path[-1]] is not None:
                path.append(prev[path[-1]])
            rs = tuple(sorted(n for n in path if g.nodes[n]['role'] == 'reaction'))
            if 3 <= len(rs) <= max_reactions:
                out.add(rs)
    return sorted(out, key=lambda rs: (len(rs), rs))


def decompose(rkeys, r2m, m2r, transport=None):
    """Returns (units, appendices, orphans). units: list of dict(shape, reactions(list), origin)."""
    cfg = CONFIG
    free = set(rkeys)
    units = []
    g = _bip(rkeys, r2m)
    # 1. cycles
    cycles = short_cycles(g, cfg['cycle_max_reactions'])
    transport = transport or set()

    def good(c):
        return len(c) >= cfg['ring_min_reactions'] and not (set(c) & transport)

    def fuse(seed):
        ring = set(seed)
        pure = not (ring & transport)
        grew = True
        while grew:  # ladder / fused rings, restricted so unrelated branches are not absorbed
            grew = False
            for c in cycles:
                cs = set(c)
                if cs <= ring or not cs <= (free | ring):
                    continue
                if len(cs & ring) < cfg['fuse_min_shared'] or len(cs - ring) > cfg['fuse_max_new'] or len(ring | cs) > cfg['max_size']:
                    continue
                if pure and (cs & transport):
                    continue
                ring |= cs
                grew = True
                break
        return ring
    # ring rule: first seed from cycles with >=6 reactions and no transport reaction, then from all remaining cycles
    for pool in ([c for c in cycles if good(c)], cycles) if cfg['ring_rule'] else (cycles,):
        for rs in pool:
            if not all(r in free for r in rs):
                continue
            ring = fuse(rs)
            units.append(dict(shape='cycle', reactions=sorted(ring), fused_cycles=len([c for c in cycles if set(c) <= ring])))
            free -= ring

    def free_at(m):
        return sorted(r for r in m2r[m] if r in free)

    # 2. forks
    def arms_of(c):
        if c.startswith('r:') or c in r2m:  # reaction center
            arms = []
            for m in sorted(r2m[c]):
                others = [r for r in free_at(m) if r != c]
                if others:
                    arms.append((m, others[0]))
            return arms
        return [(c, r) for r in free_at(c)]
    while True:
        centers = []
        for m in sorted(m2r):
            n = len(free_at(m))
            if n >= 3:
                centers.append((-n, m))
        for r in sorted(free):
            n = len(arms_of(r))
            if n >= 3:
                centers.append((-n, r))
        if not centers:
            break
        centers.sort()
        c = centers[0][1]
        is_reaction = c in r2m
        arms = arms_of(c)[:cfg['fork_max_arms']]
        members = ([c] if is_reaction else []) + [r for _, r in arms]
        members = list(dict.fromkeys(members))
        tips = [r for _, r in arms]
        depth = 1 if len(arms) > 4 else cfg['fork_arm_depth']
        taken = set(members)
        for _ in range(depth - 1):
            new_tips = []
            for t in tips:
                if len(taken) >= cfg['max_size']:
                    break
                cand = []
                for m in sorted(r2m[t]):
                    if m == c:
                        continue
                    others = [r for r in free_at(m) if r != t and r not in taken]
                    if len(others) == 1:
                        cand.append(others[0])
                if cand and len(taken) < cfg['max_size']:
                    taken.add(cand[0]); members.append(cand[0]); new_tips.append(cand[0])
            tips = new_tips
        if len(members) < cfg['min_size']:
            # cannot happen for >=3 arms; guard to guarantee progress
            free -= set(members)
            continue
        units.append(dict(shape='fork', reactions=sorted(members), center=c))
        free -= set(members)
    # 3. chains
    rg = nx.Graph()
    rg.add_nodes_from(sorted(free))
    for m in sorted(m2r):
        f = free_at(m)
        if len(f) == 2:
            rg.add_edge(*f)
    fragments = []
    for comp in sorted(sorted(c) for c in nx.connected_components(rg)):
        sub = rg.subgraph(comp)
        order = None
        if max(dict(sub.degree()).values(), default=0) <= 2:
            ends = sorted(n for n in sub if sub.degree(n) <= 1)
            start = ends[0] if ends else comp[0]
            order = [start]
            seen = {start}
            while True:
                nxt = [n for n in sorted(sub[order[-1]]) if n not in seen]
                if not nxt:
                    break
                order.append(nxt[0]); seen.add(nxt[0])
        else:  # defensive: cannot occur after the fork stage; cut at branch nodes
            fragments.append(list(comp)); continue
        if len(order) < cfg['min_size']:
            fragments.append(order); continue
        pieces = [order]
        while any(len(p) > cfg['max_size'] for p in pieces):
            pieces = [q for p in pieces for q in ([p] if len(p) <= cfg['max_size'] else [p[:len(p) // 2], p[len(p) // 2:]])]
        for p in pieces:
            units.append(dict(shape='chain', reactions=list(p)))
    # 4. merge fragments
    appendices = defaultdict(list)
    owner = {r: i for i, u in enumerate(units) for r in u['reactions']}
    frag = sorted((sorted(f) for f in fragments), key=lambda f: f[0])
    changed = True
    while changed and frag:
        changed = False
        rest = []
        for f in frag:
            adj = defaultdict(int)
            for r in f:
                for m in r2m[r]:
                    for o in m2r[m]:
                        if o in owner and o not in f:
                            adj[owner[o]] += 1
            options = sorted((i for i in adj if len(units[i]['reactions']) + len(appendices[i]) + len(f) <= cfg['max_appendix_total']),
                             key=lambda i: (len(appendices[i]) + len(units[i]['reactions']), -adj[i], i))
            if options:
                i = options[0]
                appendices[i].extend(f)
                for r in f:
                    owner[r] = i
                changed = True
            else:
                rest.append(f)
        frag = rest
    # fragments adjacent only to other fragments: union by shared metabolites
    leftover = [r for f in frag for r in f]
    orphans = []
    if leftover:
        lg = nx.Graph()
        lg.add_nodes_from(leftover)
        for m in sorted(m2r):
            f = [r for r in m2r[m] if r in set(leftover)]
            for a, b in itertools.combinations(sorted(f), 2):
                lg.add_edge(a, b)
        for comp in sorted(sorted(c) for c in nx.connected_components(lg)):
            if len(comp) >= cfg['min_size'] and len(comp) <= cfg['max_size']:
                units.append(dict(shape='chain', reactions=list(comp), origin='fragment_union'))
            else:
                orphans.append(comp)
    return units, appendices, orphans


# ---------------------------------------------------------------- unit records

def build_unit(d, u, appendix, unit_id):
    rks = sorted(set(u['reactions']) | set(appendix))
    rx = d['reactions']
    keys = set(rks)
    mets = sorted({m for (r, m) in d['edges'] if r in keys})
    keys |= set(mets)
    original = {}
    nodes = []
    for k in sorted(keys):
        if k in rx:
            r = rx[k]
            original[k] = d['rxy'][k]
            nodes.append(dict(key=k, role='reaction', xy=None, entity=r['id'], observed=d['rorigin'][k] == 'observed', position_origin=d['rorigin'][k]))
        else:
            g = d['glyphs'][k]
            original[k] = g['xy']
            nodes.append(dict(key=k, role='metabolite', xy=None, entity=g['entity'], observed=True, position_origin='observed', label=g.get('label', '')))
    norm = normalize(original)
    for n in nodes:
        n['xy'] = norm[n['key']]
    edges = []
    for (r, m), c in sorted(d['edges'].items()):
        if r in keys and m in keys:
            edges.append(dict(a=r, b=m, flow='substrate' if c < 0 else 'product', coefficient=c))
    return dict(id=unit_id, source=d['source'], species=d['species'], lineage=d['lineage'], map_id=d['id'], name=d['name'], license=d['license'],
                quality='geometry_only' if d.get('geometry_only') else 'semantic_transfer', size=len(rks), nodes=nodes, edges=edges, reactions=[] if d.get('geometry_only') else [rx[k] for k in rks], attachments=[], ports=[],
                geometry_only_flag=bool(d.get('geometry_only')),
                parent_map=dict(map_id=d['id'], name=d['name'], source=d['source'], parent_label_kind='source map / subsystem name (a label, not a unit)'),
                source_sha256=d['source_sha256'], original_instance=d['original_instance'], original_coordinates=original,
                review_status='not_required', partition='open_core', admitted_to_open_core=False,
                shape_type=u['shape'], appendix_reactions=sorted(appendix), core_reaction_count=len(u['reactions']),
                unit_definition='non-overlapping structure unit (see units_library.py); coordinates are the source map coordinates, normalised',
                biochemical_visibility='complete reference reaction; core metabolites only (currency removed)')


def is_transport(record):
    comps = {comp(e) for e in record.get('coefficients', {})} - {'unknown', 'local'}
    return len(comps) > 1


def metabolite_adjacency(d):
    r2m, m2r = defaultdict(set), defaultdict(set)
    for (r, m) in d['edges']:
        r2m[r].add(m); m2r[m].add(r)
    return r2m, m2r


def convention_stats(d, acc):
    """Orientation / angle / spacing priors from 1- and 2-reaction pieces of a source map (never stored as modules)."""
    r2m, m2r = metabolite_adjacency(d)
    src = d['source']
    dists = []
    vectors = {}
    for r, ms in sorted(r2m.items()):
        subs = [m for m in ms if d['edges'][(r, m)] < 0]
        prods = [m for m in ms if d['edges'][(r, m)] > 0]
        if len(subs) == 1 and len(prods) == 1:
            a, b = np.array(d['glyphs'][subs[0]]['xy']), np.array(d['glyphs'][prods[0]]['xy'])
            v = b - a
            n = float(np.linalg.norm(v))
            if n > 1e-6:
                vectors[r] = v / n
                dists.append(n)
    scale = float(np.median(dists)) if dists else 1.
    s = acc[src]
    s['maps'] += 1
    for r, v in vectors.items():
        ang = math.degrees(math.atan2(-v[1], v[0])) % 360  # screen y points down; report mathematical orientation
        s['single_reaction_direction_octants'][int(((ang + 22.5) % 360) // 45)] += 1
        axis = min(abs(math.sin(math.radians(ang))), abs(math.cos(math.radians(ang))))
        s['single_reaction_axis_aligned_5deg'] += int(axis < math.sin(math.radians(5)))
        s['single_reactions'] += 1
        a = np.array(d['glyphs'][[m for m in r2m[r] if d['edges'][(r, m)] < 0][0]]['xy']); b = np.array(d['glyphs'][[m for m in r2m[r] if d['edges'][(r, m)] > 0][0]]['xy'])
        s['_spacing'].append(float(np.linalg.norm(b - a)) / scale)
    for m, rs in sorted(m2r.items()):
        rs = sorted(r for r in rs if r in vectors)
        if len(rs) == 2:  # two reactions sharing a metabolite = a 2-reaction piece
            a, b = vectors[rs[0]], vectors[rs[1]]
            # orient both along the pair's shared metabolite: consecutive if one produces and the other consumes it
            sa, sb = d['edges'][(rs[0], m)], d['edges'][(rs[1], m)]
            if sa * sb < 0:
                turn = math.degrees(math.acos(max(-1., min(1., float(np.dot(a, b)))))) if sa > 0 else math.degrees(math.acos(max(-1., min(1., float(np.dot(a, b))))))
                s['consecutive_pairs'] += 1
                s['pair_turn_bins_30deg'][min(5, int(turn // 30))] += 1
    s['_scale'].append(scale)


def drawing_conventions(maps):
    acc = defaultdict(lambda: dict(maps=0, single_reactions=0, single_reaction_axis_aligned_5deg=0, single_reaction_direction_octants=[0] * 8,
                                   consecutive_pairs=0, pair_turn_bins_30deg=[0] * 6, _spacing=[], _scale=[]))
    for d in maps:
        if not d.get('geometry_only'):
            convention_stats(d, acc)
    out = {}
    for src, s in sorted(acc.items()):
        sp = np.array(s.pop('_spacing') or [0.])
        s.pop('_scale')
        s['spacing_substrate_to_product_over_map_median'] = {q: round(float(np.percentile(sp, p)), 3) for q, p in (('p10', 10), ('p25', 25), ('p50', 50), ('p75', 75), ('p90', 90))}
        s['axis_aligned_fraction'] = round(s['single_reaction_axis_aligned_5deg'] / max(1, s['single_reactions']), 4)
        s['pair_turn_note'] = 'angle between the direction vectors of two consecutive reactions: bin i covers [30i, 30(i+1)) degrees, last bin 150+; 0 = straight'
        out[src] = s
    return out


# ---------------------------------------------------------------- build

class _Matcher(nx.isomorphism.GraphMatcher):
    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self.steps = 0

    def syntactic_feasibility(self, a, b):
        self.steps += 1
        if self.steps > CONFIG['step_budget']:
            raise TimeoutError
        return super().syntactic_feasibility(a, b)


def all_mappings(a, b, cap):
    """Role-preserving isomorphisms a->b, at most `cap`; returns (mappings, budget_hit)."""
    out = []
    matcher = _Matcher(a, b, node_match=lambda x, y: x['role'] == y['role'])
    try:
        for mp in matcher.isomorphisms_iter():
            out.append(mp)
            if len(out) >= cap:
                return out, 0 if cap <= 8 else 1
        return out, 0
    except TimeoutError:
        return out, 1


def gz_dump(path, obj):
    with open(path, 'wb') as raw:
        with gzip.GzipFile(fileobj=raw, mode='wb', mtime=0, filename='') as f:
            f.write(json.dumps(obj, ensure_ascii=False, sort_keys=True).encode())


def build():
    EXT.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    maps = source_core_maps()
    mods, per_map, orphan_rows = [], [], []
    for d in maps:
        r2m, m2r = metabolite_adjacency(d)
        rks = sorted(r2m)
        units, appendices, orphans = decompose(rks, r2m, m2r, transport={r for r in rks if not d.get('geometry_only') and is_transport(d['reactions'][r])})
        owner = {}
        built = []
        for i, u in enumerate(units):
            uid_ = uid(['structure-unit', d['id'], sorted(set(u['reactions']) | set(appendices.get(i, [])))])
            m = build_unit(d, u, appendices.get(i, []), uid_)
            if len(m['nodes']) > CONFIG['max_module_nodes'] or not nx.is_connected(database.graph(m)):
                orphans.append(sorted(set(u['reactions']) | set(appendices.get(i, []))))
                continue
            built.append(m)
            for n in m['nodes']:
                if n['role'] == 'reaction':
                    owner[n['key']] = m['id']
        for m in built:  # interface ports: core metabolites shared with other units
            keys = {n['key'] for n in m['nodes'] if n['role'] == 'reaction'}
            ports = []
            for n in m['nodes']:
                if n['role'] != 'metabolite':
                    continue
                outside = sorted(r for r in m2r[n['key']] if r not in keys)
                if outside:
                    ports.append(dict(key=n['key'], metabolite=n['entity'], outside_reactions=outside, neighbour_units=sorted({owner[r] for r in outside if r in owner}), scope='whole_source_map_core'))
            m['ports'] = ports
        ids = sorted(m['id'] for m in built)
        for m in built:
            m['sibling_units'] = [x for x in ids if x != m['id']]
        mods.extend(built)
        sizes = [m['size'] for m in built]
        per_map.append(dict(map=d['id'], source=d['source'], name=d['name'], core_reactions=len(rks), units=len(built), unit_reactions=sum(sizes), orphan_reactions=len({r for o in orphans for r in o}),
                            shapes=dict(Counter(m['shape_type'] for m in built)), appendix_reactions=sum(len(m['appendix_reactions']) for m in built)))
        orphan_rows.extend(dict(map=d['id'], reactions=o) for o in orphans)
        if d['source'] == 'escher':
            per_map[-1]['window_count'] = d['window_count']; per_map[-1]['max_coordinate_drift'] = round(d['max_coordinate_drift'], 6)
    mods.sort(key=lambda m: m['id'])
    # families and variants.  Family = role-preserving graph isomorphism.  Variant = same drawing: the BEST geometry error over
    # ALL isomorphisms (exhaustive up to a mapping/step budget; cap hits are recorded) is below the threshold.  The existing
    # library only looked at the first 8 isomorphisms; that rule is tracked in parallel so the difference can be reported.
    families = defaultdict(list)
    familyrows = []
    audit = Counter()
    db_path = EXT / 'atlas.sqlite'
    if db_path.exists():
        db_path.unlink()
    db = sqlite3.connect(db_path)
    db.executescript('''CREATE TABLE modules(id TEXT PRIMARY KEY,size INTEGER,core_family TEXT,variant TEXT,shape_type TEXT,payload BLOB);
    CREATE TABLE reactions(module_id TEXT,reaction_key TEXT,source_label TEXT,PRIMARY KEY(module_id,reaction_key));
    CREATE INDEX family_idx ON modules(core_family);CREATE INDEX reaction_idx ON reactions(source_label);''')
    for i, m in enumerate(mods):
        gm = database.graph(m)
        h = nx.weisfeiler_lehman_graph_hash(gm, node_attr='role')
        bucket = (len(gm), gm.number_of_edges(), h)
        family = None
        for f in families[bucket]:
            mp, capped = all_mappings(gm, f['graph'], 1)
            audit['family_search_budget_hits'] += capped
            if mp:
                family = f
                break
        if family is None:
            family = dict(id=uid(['units-core', bucket, m['id']]), graph=gm, variants=[], variants8=[], instances=[])
            families[bucket].append(family); familyrows.append(family)
        best_variant = None
        for v in family['variants']:
            mp, capped = all_mappings(gm, v['graph'], CONFIG['mapping_cap'])
            audit['variant_cap_hits'] += capped
            audit['variant_comparisons'] += 1
            if mp and min(database.geometry_error(m, v['representative'], x) for x in mp) < CONFIG['variant_threshold']:
                best_variant = v
                break
        if best_variant is None:
            best_variant = dict(id=uid([family['id'], m['id']]), graph=gm, representative=m, instances=[])
            family['variants'].append(best_variant)
        # the existing rule (first 8 isomorphisms), tracked only for the comparison
        placed8 = False
        for v in family['variants8']:
            mp, _ = all_mappings(gm, v['graph'], 8)
            if mp and min(database.geometry_error(m, v['representative'], x) for x in mp) < CONFIG['variant_threshold']:
                placed8 = True
                break
        if not placed8:
            family['variants8'].append(dict(graph=gm, representative=m))
        variant = best_variant
        family['instances'].append(m['id']); variant['instances'].append(m['id'])
        m['core_family'] = family['id']; m['variant'] = variant['id']
        arr = np.array([n['xy'] for n in m['nodes']])
        m['geometry_fingerprint'] = uid(np.round(np.sort(np.linalg.norm(arr[:, None] - arr[None, :], axis=2).ravel()), 3).tolist())
        db.execute('INSERT INTO modules VALUES(?,?,?,?,?,?)', (m['id'], m['size'], m['core_family'], m['variant'], m['shape_type'], database.pack(m)))
        for r in m['reactions']:
            db.execute('INSERT INTO reactions VALUES(?,?,?)', (m['id'], r['key'], r['id']))
        if i % 2000 == 0:
            print('units index', i, '/', len(mods), flush=True)
    db.executescript('CREATE TABLE maps(map_id TEXT PRIMARY KEY,name TEXT,source TEXT,unit_ids TEXT);')
    by_map = defaultdict(list)
    for m in mods:
        by_map[m['map_id']].append(m)
    for map_id, ms in sorted(by_map.items()):
        db.execute('INSERT INTO maps VALUES(?,?,?,?)', (map_id, ms[0]['name'], ms[0]['source'], json.dumps([m['id'] for m in ms])))
    db.commit(); db.close()
    gz_dump(EXT / 'modules.json.gz', mods)
    dump(EXT / 'families.json', dict(families=[dict(id=f['id'], instances=f['instances'], variants=[dict(id=v['id'], template=v['representative']['id'], instances=v['instances']) for v in f['variants']]) for f in familyrows]))
    licenses = Counter(d['license'] for d in maps)
    source = dict(extension='structure units (non-overlapping recut of the open-core source maps)', config=CONFIG, partition='open_core', admitted_to_open_core=False,
                  source_maps=len(maps), source_map_ids_by_source=dict(Counter(d['source'] for d in maps)), licenses=dict(licenses),
                  license_note='Each unit keeps the licence, source_sha256 and original_instance pointer of its source map; no data from outside the open core.',
                  changes='Source maps recut into non-overlapping 3-12 reaction structure units; coordinates normalised per unit, originals kept in original_coordinates.',
                  geometry_only='WikiPathways maps without reaction identity are segmented by shape only (pseudo-reactions = drawn interactions between metabolite glyphs) and flagged quality=geometry_only; maps that also have semantic units are not duplicated.')
    groups = defaultdict(list)
    for m in mods:
        if m['reactions']:
            groups[uid(sorted(r['id'] for r in m['reactions']))].append(m['id'])
    dump(RESULTS / 'unit_reaction_groups.json', dict(note='units with the same set of reaction identifiers (across source maps); distinct sets and the units sharing each',
                                                      groups=len(groups), units=sum(len(v) for v in groups.values()), multi_map_groups=sum(len(v) > 1 for v in groups.values()),
                                                      members={k: v for k, v in sorted(groups.items())}))
    dump(EXT / 'SOURCE.json', source)
    sizes = [m['size'] for m in mods]
    summary = dict(source_maps=len(maps), modules=len(mods), unit_reactions=sum(sizes), source_core_reactions=sum(p['core_reactions'] for p in per_map),
                   orphan_reactions=sum(p['orphan_reactions'] for p in per_map), sizes=dict(sorted(Counter(sizes).items())), mean_size=round(float(np.mean(sizes)), 3),
                   shapes=dict(Counter(m['shape_type'] for m in mods)), appendix_reactions=sum(len(m['appendix_reactions']) for m in mods),
                   by_source=dict(Counter(m['source'] for m in mods)), families=len(familyrows), variants=sum(len(f['variants']) for f in familyrows), variants_under_first8_rule=sum(len(f['variants8']) for f in familyrows), alignment_audit=dict(audit),
                   singleton_families=sum(len(f['instances']) == 1 for f in familyrows), max_unit_size=max(sizes), over_12=sum(s > 12 for s in sizes),
                   overlap_check_reactions_in_multiple_units=sum(c > 1 for c in Counter((m['map_id'], n['key']) for m in mods for n in m['nodes'] if n['role'] == 'reaction').values()),
                   by_quality=dict(Counter(m['quality'] for m in mods)), source_layouts_with_units=len({m['map_id'] for m in mods}),
                   config=CONFIG, hash_seed=os.environ.get('PYTHONHASHSEED', 'random'),
                   code_sha256=sha(Path(__file__).read_bytes()),
                   limits=['Whole-map core graphs are rebuilt from stored source records: Escher maps as the union of their stored windows, WikiPathways/Reactome from GPML/diagram records.',
                           'Orphan islands (<3 reactions with no adjacent unit) are not modules.', 'Shape type is the formation gestalt; appendices may enlarge a unit beyond 12 (<=16).'])
    dump(RESULTS / 'library_summary.json', summary)
    dump(RESULTS / 'per_map.json', per_map)
    dump(RESULTS / 'orphans.json', orphan_rows)
    dump(RESULTS / 'drawing_conventions.json', dict(note='Orientation / angle / spacing priors from 1- and 2-reaction pieces (substrate->product vector of one-substrate/one-product reactions, consecutive reaction pairs). Not modules.',
                                                    by_source=drawing_conventions(maps)))
    return mods, summary


def derive_windows(parents, sizes=(2, 4, 8, 12)):
    """Windows derived from stored (non-redundant) modules: BFS sub-structures inside ONE parent, coordinates cropped from the parent.

    Same traversal as `retrieval.windows` but exhaustive over roots and confined to the parent.  Returned in the module record format
    (so retrieval.Library / semantic_transfer.Mapper / ec_channel.ECMapper accept them); the parent is `derived_from`.
    """
    out = []
    for parent in sorted(parents, key=lambda m: m['id']):
        rnodes = [n['key'] for n in parent['nodes'] if n['role'] == 'reaction']
        if len(rnodes) < 3:
            continue
        g = database.graph(parent)
        seen = set()
        for root in sorted(rnodes):
            queue, visited, order = [root], {root}, []
            for k in queue:
                if g.nodes[k]['role'] == 'reaction':
                    order.append(k)
                if len(order) >= max(sizes):
                    break
                for nb in sorted(g[k]):
                    if nb not in visited:
                        visited.add(nb); queue.append(nb)
            for n in sizes:
                if n < len(rnodes) and len(order) >= n:
                    rs = tuple(sorted(order[:n]))
                    if rs not in seen:
                        seen.add(rs)
                        out.append(_crop(parent, g, rs))
    return out


def _crop(parent, g, rs):
    keys = set(rs)
    for r in rs:
        keys.update(g[r])
    nodes = [dict(n) for n in parent['nodes'] if n['key'] in keys]
    xy = normalize({n['key']: n['xy'] for n in nodes})
    for n in nodes:
        n['xy'] = xy[n['key']]
    edges = [e for e in parent['edges'] if e['a'] in keys and e['b'] in keys]
    recs = [r for r in parent['reactions'] if r['key'] in set(rs)]
    ports = []
    for n in nodes:
        if n['role'] == 'metabolite':
            outside = sorted(r for r in g[n['key']] if r not in set(rs))
            if outside:
                ports.append(dict(key=n['key'], metabolite=n.get('entity'), outside_reactions=outside, scope='cropped_from_parent'))
    arr = np.array([n['xy'] for n in nodes])
    gm = nx.Graph()
    for n in nodes:
        gm.add_node(n['key'], role=n['role'])
    gm.add_edges_from((e['a'], e['b']) for e in edges)
    bucket = (len(gm), gm.number_of_edges(), nx.weisfeiler_lehman_graph_hash(gm, node_attr='role'))
    m = {k: v for k, v in parent.items() if k not in ('nodes', 'edges', 'reactions', 'ports', 'original_coordinates', 'sibling_units', 'appendix_reactions')}
    m.update(id=uid(['derived-window', parent['id'], rs]), size=len(rs), nodes=nodes, edges=edges, reactions=recs, ports=ports, derived_from=parent['id'], shape_type='window',
             core_family=uid(['derived-core', bucket]), variant=uid(['derived-variant', bucket]),
             geometry_fingerprint=uid(np.round(np.sort(np.linalg.norm(arr[:, None] - arr[None, :], axis=2).ravel()), 3).tolist()))
    return m


def modules(kind='units'):
    """Library switch for consumers (retrieval.Library, semantic_transfer.Mapper, ec_channel.ECMapper, tool_pipeline):
    pass `units_library.modules(kind)` where they now pass `database.modules()`.
      windows       the existing v0.3 overlapping-window library (unchanged)
      maximal       only the windows not contained in a larger window of the same map (results/module_redundancy.json.gz)
      maximal_sub   maximal + sub-windows derived from them at load time
      units         stored non-overlapping structure units
      units_sub     units + sub-windows derived from them (the hybrid: de-duplicated storage, windows derived on demand)
    """
    if kind == 'windows':
        return database.modules()
    if kind in ('maximal', 'maximal_sub'):
        flags = json.load(gzip.open(OUT / 'module_redundancy.json.gz', 'rt'))
        base = [m for m in database.modules() if flags[m['id']]['maximal']]
    elif kind in ('units', 'units_sub'):
        base = load()
    else:
        raise ValueError(kind)
    return base + (derive_windows(base) if kind.endswith('_sub') else [])


def load():
    with gzip.open(EXT / 'modules.json.gz', 'rt') as f:
        return json.load(f)


if __name__ == '__main__':
    if os.environ.get('PYTHONHASHSEED') != '0':
        os.environ['PYTHONHASHSEED'] = '0'
        os.execv(sys.executable, [sys.executable, str(Path(__file__).resolve())] + sys.argv[1:])
    print(json.dumps(build()[1], indent=1))
