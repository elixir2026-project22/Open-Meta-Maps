"""Independent, provenance-preserving BRENDA pathway-map extension (trial ingestion).

Source: the public BRENDA combined pathway map (SVG, CC BY 4.0). It is stored under
data/references/brenda/ and decoded here into pathway-level modules in the same record
format as the MicroMap extension.  Status: *unverified layout origin* -- the node and edge
geometry shows features of an automatic layered layout (see `layout_audit`), so these
modules are NOT counted as expert designs and are NOT admitted to the open core until the
data provider confirms how positions were produced.
"""
import itertools
import math
import re
import sqlite3
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict

from core import *
import database

EXT = DATA / 'extensions/brenda'
RESULTS = OUT / 'brenda_extension'
SVG = DATA / 'references/brenda/BRENDA_map.svg'
NS = '{http://www.w3.org/2000/svg}'
CONFIG = dict(window_sizes=[4, 8, 12, 24], max_module_nodes=100, pair_degree_cap=8,
              link_tolerance=8.0, variant_threshold=.25, variant_representatives=12)
KINDS = {'nodeEN': 'reaction', 'nodeHM': 'metabolite', 'nodeSM': 'metabolite', 'nodeCM': 'cofactor',
         'nodePathwayTitle': 'title'}


def _tokens(g):
    return g.get('class', '').split()


def _node(g):
    toks = _tokens(g)
    kind = next((t for t in toks if t in KINDS), None)
    if kind is None:
        return None
    m = re.match(r'translate\(([-\d.]+),([-\d.]+)\)', g.get('transform', ''))
    rect = g.find(NS + 'rect')
    text = g.find(NS + 'text')
    anchor = text.find(NS + 'a') if text is not None else None
    title = g.find(NS + 'title')
    href = anchor.get('href', '') if anchor is not None else ''
    n = dict(id=g.get('id'), kind=kind, role=KINDS[kind], pathway=next((t[3:] for t in toks if t.startswith('pw_')), None),
             x=float(m[1]), y=float(m[2]), w=float(rect.get('width')), h=float(rect.get('height')),
             name=(title.text or '').strip() if title is not None else '',
             label=''.join(anchor.itertext()).strip() if anchor is not None else ''.join(text.itertext()).strip() if text is not None else '')
    for t in toks:
        if t.startswith('rID'):
            n['reaction_token'] = t[3:]
        elif re.fullmatch(r'lid\d+', t):
            n['ligand_id'] = t[3:]
        elif re.fullmatch(r'gid\d+', t):
            n['group_id'] = t[3:]
    ec = re.search(r'ecno=([^&"]+)', href)
    if ec:
        n['ec'] = ec[1]
    lig = re.search(r'brenda_ligand_id=(\d+)', href)
    if lig:
        n['ligand_id'] = lig[1]
    return n


def decode(path=SVG):
    """Decode the combined map into pathways, nodes and node-to-node links (no coordinates altered)."""
    raw = Path(path).read_bytes()
    root = ET.fromstring(raw)
    nodes = [n for n in (_node(g) for g in root.iter(NS + 'g') if 'node' in _tokens(g)) if n]
    by_id = {n['id']: n for n in nodes}
    assert len(by_id) == len(nodes), 'duplicate node ids'
    pathways = {}
    for poly in root.iter(NS + 'polygon'):
        pid = poly.get('id', '')
        if not pid.startswith('pw_'):
            continue
        pts = [tuple(map(float, p.split(','))) for p in poly.get('points').split()]
        xs, ys = [p[0] for p in pts], [p[1] for p in pts]
        pathways[pid[3:]] = dict(slug=pid[3:], box=[min(xs), min(ys), max(xs), max(ys)], polygon_points=len(pts))
    for n in nodes:
        slug = n['pathway']
        if slug and slug not in pathways:
            pathways[slug] = dict(slug=slug, box=None, polygon_points=0)
    for n in nodes:
        if n['kind'] == 'nodePathwayTitle' and n['pathway'] in pathways:
            pathways[n['pathway']]['title'] = n['label'] or n['name']
    # A polygon id can differ from its nodes' pathway class (curcuminoid vs curcumin): pair an empty polygon with a
    # polygon-less node class when that class's title node lies inside the box.
    used = {n['pathway'] for n in nodes if n['pathway']}
    titles = {n['pathway']: n for n in nodes if n['kind'] == 'nodePathwayTitle'}
    for slug in [s for s, p in pathways.items() if s not in used and p['box']]:
        box = pathways[slug]['box']
        for other, p in list(pathways.items()):
            t = titles.get(other)
            if other in used and p['box'] is None and t and box[0] <= t['x'] <= box[2] and box[1] <= t['y'] <= box[3]:
                p.update(box=box, polygon_points=pathways[slug]['polygon_points'], polygon_id_alias=slug)
                del pathways[slug]
                break
    # spatial index of node rectangles
    G = 200.0
    grid = defaultdict(list)
    for n in nodes:
        if n['kind'] == 'nodePathwayTitle':
            continue
        for gx in range(int((n['x'] - n['w'] / 2) // G) - 1, int((n['x'] + n['w'] / 2) // G) + 2):
            for gy in range(int((n['y'] - n['h'] / 2) // G) - 1, int((n['y'] + n['h'] / 2) // G) + 2):
                grid[(gx, gy)].append(n)

    def nearest(px, py):
        best = None
        for n in grid.get((int(px // G), int(py // G)), []):
            dx = max(abs(px - n['x']) - n['w'] / 2, 0)
            dy = max(abs(py - n['y']) - n['h'] / 2, 0)
            d = math.hypot(dx, dy)
            if best is None or d < best[0]:
                best = (d, n)
        return best

    links, unmatched = [], Counter()
    groups = {g.get('id'): g for g in root.iter(NS + 'g') if g.get('id') in ('all-link', 'all-linkCof')}
    for gid, group in groups.items():
        for p in group.iter(NS + 'path'):
            m = re.fullmatch(r'M ([-\d.]+) ([-\d.]+)((?: L [-\d.]+ [-\d.]+)+)', p.get('d', '').strip())
            if not m:
                unmatched[gid + ':non_simple_path'] += 1
                continue
            pts = [[float(m[1]), float(m[2])]] + [[float(a), float(b)] for a, b in re.findall(r'L ([-\d.]+) ([-\d.]+)', m[3])]
            (x1, y1), (x2, y2) = pts[0], pts[-1]
            a, b = nearest(x1, y1), nearest(x2, y2)
            if not a or not b or a[0] > CONFIG['link_tolerance'] or b[0] > CONFIG['link_tolerance'] or a[1]['id'] == b[1]['id']:
                unmatched[gid + ':endpoint_not_on_node'] += 1
                continue
            links.append(dict(a=a[1]['id'], b=b[1]['id'], group=gid, inter_pathway='inter-pw-link' in p.get('class', ''),
                              markers=[p.get('marker-start') is not None, p.get('marker-end') is not None],
                              geometry=pts))
    return dict(canvas=[float(root.get('width')), float(root.get('height'))], pathways=pathways, nodes=by_id,
                links=links, unmatched_links=dict(unmatched), source_sha256=sha(raw))


def layout_audit(d):
    """Measured geometric features. Descriptive only; this cannot prove how positions were made."""
    nodes = [n for n in d['nodes'].values() if n['kind'] != 'nodePathwayTitle']
    per = defaultdict(list)
    for n in nodes:
        per[n['pathway']].append(n)
    sharey, sharex, ylevels = [], [], []
    for ns in per.values():
        if len(ns) < 8:
            continue
        ys, xs = Counter(round(n['y'], 1) for n in ns), Counter(round(n['x'], 1) for n in ns)
        sharey.append(sum(v for v in ys.values() if v > 1) / len(ns))
        sharex.append(sum(v for v in xs.values() if v > 1) / len(ns))
        ylevels.append(len(ys) / len(ns))
    med = lambda v: float(sorted(v)[len(v) // 2]) if v else None
    angles = Counter()
    for l in d['links']:
        (x1, y1), (x2, y2) = l['geometry'][0], l['geometry'][-1]
        th = abs(math.degrees(math.atan2(y2 - y1, x2 - x1))) % 90
        angles['axis_aligned' if min(th, 90 - th) < 2 else 'oblique'] += 1
    return dict(
        nodes=len(nodes), coordinates_with_integer_x=float(np.mean([n['x'] == int(n['x']) for n in nodes])),
        pathways_ge8_nodes=len(sharey), median_share_nodes_sharing_y_with_another=med(sharey),
        median_share_nodes_sharing_x_with_another=med(sharex), median_distinct_y_levels_per_node=med(ylevels),
        link_angles=dict(angles), links_with_bends=sum(len(l['geometry']) > 2 for l in d['links']),
        pathway_boxes_rectangular=sum(p['polygon_points'] == 4 for p in d['pathways'].values()), pathway_boxes=len(d['pathways']),
        interpretation='Rows/columns of identical coordinates inside each pathway and rectangular pathway boxes are '
                       'consistent with an automatic layered layout packed onto one canvas. This is evidence, not proof; '
                       'the provider has not confirmed how positions were produced.')


def pathway_graph(d, slug):
    """Reaction nodes, core compound nodes (HM/SM) and cofactor attachments of one pathway."""
    nodes = {k: n for k, n in d['nodes'].items() if n['pathway'] == slug}
    core, side, other = [], [], Counter()
    for l in d['links']:
        a, b = d['nodes'][l['a']], d['nodes'][l['b']]
        if a['pathway'] != slug or b['pathway'] != slug:
            continue
        if l['inter_pathway']:
            other['inter_pathway_link'] += 1
            continue
        roles = {a['role'], b['role']}
        if roles == {'reaction', 'metabolite'}:
            r, m = (a, b) if a['role'] == 'reaction' else (b, a)
            core.append((r['id'], m['id'], l))
        elif roles == {'reaction', 'cofactor'}:
            r, c = (a, b) if a['role'] == 'reaction' else (b, a)
            side.append((r['id'], c['id'], l))
        else:
            other['unparsed_' + '-'.join(sorted(roles))] += 1
    return nodes, core, side, dict(other)


def reaction_record(n):
    refs = {}
    if n.get('ec'):
        refs['ec-code'] = [n['ec']]
    if n.get('reaction_token'):
        refs['brenda.reaction.svg_token'] = [n['reaction_token']]
    return dict(key='r:' + n['id'], id=n.get('reaction_token') or n['id'], source_reaction_id=n.get('reaction_token'),
                ec=n.get('ec'), enzyme_name=n['name'], references=refs, coefficients={}, coefficients_origin='not_in_source_svg', reversible=None,
                reversibility_origin='not_in_source_svg', xy=[n['x'], n['y']], position_origin='observed_node_center_from_svg')


def extract(d):
    modules, excluded = [], []
    cite = dict(title='BRENDA pathway map (combined SVG)', doi='10.1093/nar/gkaf1113', url='https://brenda-enzymes.org/pathways/BRENDA_map.svg',
                paper='BRENDA in 2026, Nucleic Acids Research', license=dict(name='CC BY 4.0', uri='https://creativecommons.org/licenses/by/4.0/'),
                retrieved='2026-11-09', source_sha256=d['source_sha256'],
                changes='Pathway-level fragments decoded from the combined SVG; coordinates normalized for storage; originals kept.')
    audit_note = dict(status='unverified', counted_as_expert_design=False,
                      reason='geometry features are consistent with an automatic layered layout; provider confirmation pending')
    for slug in sorted(d['pathways']):
        nodes, core, side, other = pathway_graph(d, slug)
        rx = {n['id']: n for n in nodes.values() if n['role'] == 'reaction'}
        g = nx.Graph()
        for rid in rx:
            g.add_node(rid, role='reaction')
        for r, m, l in core:
            g.add_node(m, role='metabolite')
            g.add_edge(r, m)
        cof = defaultdict(list)
        for r, c, l in side:
            cof[r].append(c)
        uses = defaultdict(set)
        for r, m, l in core:
            uses[m].add(r)
        subsets = {(r,) for r in rx}
        for m in g:
            if g.nodes[m]['role'] == 'metabolite' and len(g[m]) <= CONFIG['pair_degree_cap']:
                subsets.update(tuple(sorted(p)) for p in itertools.combinations(sorted(g[m]), 2))
        for root in sorted(rx):
            order = [k for k in nx.bfs_tree(g, root) if k in rx][:max(CONFIG['window_sizes'])]
            for n in CONFIG['window_sizes']:
                if len(order) >= n:
                    subsets.add(tuple(sorted(order[:n])))
        for comp in nx.connected_components(g):
            rs = tuple(sorted(k for k in comp if k in rx))
            if len(rs) >= 2:
                subsets.add(rs)
        title = d['pathways'][slug].get('title') or slug
        for rs in sorted(subsets):
            keys = set(rs)
            for r in rs:
                keys.update(g[r])
            if len(keys) > CONFIG['max_module_nodes']:
                excluded.append(dict(pathway=slug, reactions=list(rs), reason='node_budget'))
                continue
            sub = g.subgraph(keys)
            if len(keys) <= len(rs) or not nx.is_connected(sub):
                excluded.append(dict(pathway=slug, reactions=list(rs), reason='empty_or_disconnected_core'))
                continue
            mnodes = []
            for k in sorted(keys):
                n = nodes[k]
                if n['role'] == 'reaction':
                    mnodes.append(dict(key='r:' + k, role='reaction', xy=[n['x'], n['y']], entity=n.get('reaction_token') or k,
                                       observed=True, position_origin='observed_node_center_from_svg', label=n['name'] or n.get('ec', ''), ec=n.get('ec')))
                else:
                    mnodes.append(dict(key='m:' + k, role='metabolite', xy=[n['x'], n['y']], entity='lid' + n.get('ligand_id', '?'),
                                       observed=True, position_origin='observed_node_center_from_svg', label=n['name'] or n['label'],
                                       ligand_id=n.get('ligand_id')))
            original = {n['key']: n['xy'] for n in mnodes}
            normalized = normalize(original)
            for n in mnodes:
                n['xy'] = normalized[n['key']]
            edges = [dict(a='r:' + r, b='m:' + m, flow='unknown', coefficient=None) for r, m, _ in core if r in rs and m in keys]
            ports = [dict(key='m:' + k, outside_reactions=sorted('r:' + r for r in uses[k] - set(rs)), scope='same_pathway_core_links')
                     for k in sorted(keys - set(rs)) if uses[k] - set(rs)]
            attach = [dict(alias='m:' + c, reaction='r:' + r, kind='cofactor', name=nodes[c]['name'], ligand_id=nodes[c].get('ligand_id'),
                           xy=[nodes[c]['x'], nodes[c]['y']]) for r in rs for c in cof.get(r, [])]
            modules.append(dict(
                id=uid(['brenda-instance', slug, rs]), source='brenda', species='pathway-level; organism assignment not in source SVG',
                lineage='brenda:' + slug, map_id=uid(['brenda-pathway', slug]), name=title, pathway_slug=slug,
                license='CC-BY-4.0', source_citation=cite, partition='research_extension', admitted_to_open_core=False,
                layout_origin=audit_note, quality='annotated_topology', size=len(rs), nodes=mnodes, edges=edges,
                reactions=[reaction_record(rx[r]) | dict(core=['m:' + m for m in g[r]]) for r in rs], attachments=attach, ports=ports,
                source_sha256=d['source_sha256'], original_instance=f'BRENDA_map.svg#pw_{slug}', original_coordinates=original,
                review_status='not_required', compartment_evidence='none in source SVG',
                biochemical_visibility='reaction and ligand identifiers per node; no stoichiometry or direction in the SVG'))
    return modules, excluded


def build():
    EXT.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    src = read(str(SVG) + '.source.json')
    assert sha(SVG.read_bytes()) == src['sha256'], 'BRENDA SVG changed; update its .source.json deliberately'
    d = decode()
    audit = layout_audit(d)
    mods, excluded = extract(d)
    dump(EXT / 'SOURCE.json', dict(src, partition='research_extension'))
    families, family_rows = defaultdict(list), []
    db = sqlite3.connect(EXT / 'atlas.sqlite')
    db.executescript('''CREATE TABLE IF NOT EXISTS sources(id TEXT PRIMARY KEY,payload BLOB);
    CREATE TABLE IF NOT EXISTS modules(id TEXT PRIMARY KEY,size INTEGER,core_family TEXT,variant TEXT,payload BLOB);
    CREATE TABLE IF NOT EXISTS reactions(module_id TEXT,reaction_key TEXT,source_label TEXT,PRIMARY KEY(module_id,reaction_key));
    CREATE INDEX IF NOT EXISTS family_idx ON modules(core_family);
    CREATE INDEX IF NOT EXISTS reaction_idx ON reactions(source_label);
    DELETE FROM sources;DELETE FROM modules;DELETE FROM reactions;''')
    db.execute('INSERT INTO sources VALUES(?,?)', (src['sha256'][:24], database.pack(dict(src, pathways=len(d['pathways'])))))
    for i, m in enumerate(mods):
        gm = database.graph(m)
        bucket = (len(gm), gm.number_of_edges(), nx.weisfeiler_lehman_graph_hash(gm, node_attr='role'))
        family = None
        for f in families[bucket]:
            if database.mappings(gm, f['graph'], 1):
                family = f
                break
        if family is None:
            family = dict(id=uid(['brenda-core', bucket, m['id']]), graph=gm, variants=[], instances=[])
            families[bucket].append(family)
            family_rows.append(family)
        variant = None
        for v in family['variants'][:CONFIG['variant_representatives']]:
            maps = database.mappings(gm, v['graph'], 4)
            if maps and min(database.geometry_error(m, v['representative'], mp) for mp in maps) < CONFIG['variant_threshold']:
                variant = v
                break
        if variant is None:
            variant = dict(id=uid([family['id'], m['id']]), graph=gm, representative=m, instances=[])
            family['variants'].append(variant)
        family['instances'].append(m['id'])
        variant['instances'].append(m['id'])
        m['core_family'], m['variant'] = family['id'], variant['id']
        arr = np.array([n['xy'] for n in m['nodes']])
        m['geometry_fingerprint'] = uid(np.round(np.sort(np.linalg.norm(arr[:, None] - arr[None, :], axis=2).ravel()), 3).tolist())
        db.execute('INSERT INTO modules VALUES(?,?,?,?,?)', (m['id'], m['size'], m['core_family'], m['variant'], database.pack(m)))
        for r in m['reactions']:
            db.execute('INSERT INTO reactions VALUES(?,?,?)', (m['id'], r['key'], r['id']))
        if i % 2000 == 0:
            print('BRENDA index', i, '/', len(mods), flush=True)
    db.commit()
    db.close()
    with gzip.open(EXT / 'modules.json.gz', 'wt') as f:
        json.dump(mods, f, ensure_ascii=False)
    dump(EXT / 'families.json', dict(families=[dict(id=f['id'], instances=f['instances'], variants=[
        dict(id=v['id'], template=v['representative']['id'], instances=v['instances']) for v in f['variants']]) for f in family_rows]))
    kinds = Counter(n['kind'] for n in d['nodes'].values())
    link_kinds = Counter('-'.join(sorted((KINDS[d['nodes'][l['a']]['kind']], KINDS[d['nodes'][l['b']]['kind']]))) + (':inter' if l['inter_pathway'] else '')
                         for l in d['links'])
    pathway_stats = {s: Counter(n['kind'] for n in d['nodes'].values() if n['pathway'] == s) for s in d['pathways']}
    summary = dict(
        source_file='data/references/brenda/BRENDA_map.svg', source_sha256=d['source_sha256'], license='CC BY 4.0',
        partition='research_extension', admitted_to_open_core=False, independent_expert_designs=0,
        layout_origin_status='unverified (features consistent with automatic layered layout); provider confirmation pending',
        canvas=d['canvas'], pathways=len(d['pathways']), pathways_with_polygon=sum(p['box'] is not None for p in d['pathways'].values()),
        nodes=dict(kinds), links_matched=len(d['links']), link_kinds=dict(link_kinds), unmatched_links=d['unmatched_links'],
        enzyme_nodes_with_reaction_token=sum(1 for n in d['nodes'].values() if n['kind'] == 'nodeEN' and n.get('reaction_token')),
        enzyme_nodes_with_ec=sum(1 for n in d['nodes'].values() if n['kind'] == 'nodeEN' and n.get('ec')),
        median_reactions_per_pathway=float(sorted(c['nodeEN'] for c in pathway_stats.values())[len(pathway_stats) // 2]),
        modules=len(mods), sizes=dict(sorted(Counter(m['size'] for m in mods).items())), families=len(family_rows),
        variants=sum(len(f['variants']) for f in family_rows), singleton_families=sum(len(f['instances']) == 1 for f in family_rows),
        excluded_module_instances=len(excluded), excluded_reasons=dict(Counter(e['reason'] for e in excluded)),
        layout_audit=audit, config=CONFIG, hash_seed=os.environ.get('PYTHONHASHSEED', 'random'), code_sha256=sha(Path(__file__).read_bytes()),
        limits=['Combined map decoded; no stoichiometry, reaction direction, compartments or organism assignment in the SVG.',
                'Positions are node centres of rectangles; link geometry is straight segments (a few links have one bend).',
                'Cofactor nodes (nodeCM) are stored as attachments, not core participants.',
                'Compound-compound, reaction-reaction and inter-pathway links are counted but not used as module edges.',
                'Layout origin unverified: not counted as expert designs and not admitted to the open core.'])
    dump(RESULTS / 'library_summary.json', summary)
    dump(RESULTS / 'excluded.json', dict(modules=excluded))
    return mods, summary


def load():
    with gzip.open(EXT / 'modules.json.gz', 'rt') as f:
        return json.load(f)


if __name__ == '__main__':
    # Variant grouping takes the first 4 isomorphisms found by database.mappings; their order depends on set iteration,
    # so the build re-executes itself with a fixed hash seed to be reproducible.
    if os.environ.get('PYTHONHASHSEED') != '0':
        os.environ['PYTHONHASHSEED'] = '0'
        os.execv(sys.executable, [sys.executable, str(Path(__file__).resolve())] + sys.argv[1:])
    print(json.dumps(build()[1], indent=1)[:3000])
