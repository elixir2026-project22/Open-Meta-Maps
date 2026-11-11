"""In-repo cutter for the 22 community Escher maps.

The historical windows in ``module_observations.json.gz`` were produced by a
program that is not in this repository. This module regenerates reaction
windows from the raw Escher JSON with a documented, replayable algorithm:

* parse each map with ``database.escher_diagram``
* connect reactions that share a **primary** metabolite **entity** (BiGG id)
  via drawn segments (duplicate glyphs of the same entity merge)
* breadth-first windows from **every** reaction root
* sizes ``1, 2, 4, 8, 12, 24, 40``, plus the full reachable component when it
  is off-grid and at most 40 reactions
* keep modules with ≤40 reactions and ≤100 nodes and a connected visible core

Parity against the stored observations is checked by
``test_escher_cutter_parity.py``. Exact set equality is not required when the
new cut preserves reaction coverage and improves retrieval gates.
"""
from __future__ import annotations

from collections import Counter, deque
from pathlib import Path

import networkx as nx
import numpy as np

from atlas_registry import LAYOUT_ORIGIN_EXPERT, SCHEMA_VERSION
from core import INPUTS, normalize, sha, uid
from database import escher_diagram, graph

WINDOW_SIZES = (1, 2, 4, 8, 12, 24, 40)
MAX_REACTIONS = 40
MAX_NODES = 100
METHOD = 'escher_community_windows'
VERSION = 'atlas-extraction-1'
RAW_ESCHER = INPUTS / 'raw' / 'escher'


def extraction_block(parameters, counts, *, inferred=False):
    """Provenance block shared with Stage A schema (method / version / parameters / counts)."""
    out = dict(method=METHOD, version=VERSION, parameters=dict(parameters), counts=dict(counts))
    if inferred:
        out['inferred'] = True
    return out


def default_parameters():
    return dict(
        max_roots='all',
        window_sizes=list(WINDOW_SIZES),
        link='primary_segment_entity',
        plus_full_component=True,
        max_reactions=MAX_REACTIONS,
        max_nodes=MAX_NODES,
    )


def raw_escher_paths(root=None):
    root = Path(root or RAW_ESCHER)
    return sorted(p for p in root.rglob('*.json') if p.is_file())


def _reaction_sort_key(key):
    tail = key[2:] if key.startswith('r:') else key
    return (0, int(tail)) if str(tail).isdigit() else (1, str(tail))


def _primary_entity_graph(diagram):
    """Bipartite graph: reactions ↔ primary metabolite entities (merged by BiGG id)."""
    g = nx.Graph()
    rx = {r['key']: r for r in diagram['events']}
    for r in diagram['events']:
        g.add_node(r['key'], role='reaction')
        for mk in r.get('visible') or {}:
            glyph = diagram['glyphs'].get(mk) or {}
            if not glyph.get('primary', True):
                continue
            ent = str(glyph.get('entity') or '').split('@')[0]
            if not ent:
                continue
            ek = 'e:' + ent
            g.add_node(ek, role='metabolite', entity=ent)
            g.add_edge(r['key'], ek)
    return g, rx


def reaction_windows(diagram, *, plus_full_component=True, max_reactions=MAX_REACTIONS):
    """Return distinct ordered reaction-key tuples (BFS prefixes / full components)."""
    g, rx = _primary_entity_graph(diagram)
    subsets = set()
    roots = sorted(rx, key=_reaction_sort_key)
    for root in roots:
        seen = {root}
        q = deque([root])
        order = []
        while q:
            u = q.popleft()
            if u in rx:
                order.append(u)

            def nbr_key(v, _rx=rx):
                if v in _rx:
                    return _reaction_sort_key(v)
                return (2, str(v))

            for v in sorted((x for x in g[u] if x not in seen), key=nbr_key):
                seen.add(v)
                q.append(v)
        for n in WINDOW_SIZES:
            if len(order) >= n:
                subsets.add(tuple(sorted(order[:n], key=_reaction_sort_key)))
        if plus_full_component and 1 < len(order) <= max_reactions and len(order) not in WINDOW_SIZES:
            subsets.add(tuple(sorted(order, key=_reaction_sort_key)))
    return [rs for rs in subsets if len(rs) <= max_reactions]


def _module_from_window(diagram, rs, *, path, parameters, map_counts):
    rx = {r['key']: r for r in diagram['events']}
    selected = []
    for k in rs:
        r = dict(rx[k])
        r.setdefault('references', {'bigg.reaction': [r['id']]})
        selected.append(r)
    keys = set(rs)
    for r in selected:
        for mk in r.get('visible') or {}:
            glyph = diagram['glyphs'].get(mk) or {}
            if glyph.get('primary', True):
                keys.add(mk)
    nodes = []
    for k in sorted(keys):
        if k in rx:
            nodes.append(dict(
                key=k, role='reaction', xy=rx[k]['xy'], entity=rx[k]['id'],
                observed=True, position_origin=rx[k].get('position_origin', 'observed'),
            ))
        else:
            glyph = diagram['glyphs'][k]
            nodes.append(dict(
                key=k, role='metabolite', xy=glyph['xy'], entity=glyph['entity'],
                observed=True, position_origin='observed', label=glyph.get('label', ''),
            ))
    if len(nodes) > MAX_NODES:
        return None
    xy = normalize({n['key']: n['xy'] for n in nodes})
    for n in nodes:
        n['xy'] = xy[n['key']]
    edges = []
    for r in selected:
        for mk, coeff in (r.get('visible') or {}).items():
            if mk not in keys:
                continue
            edges.append(dict(
                a=r['key'], b=mk,
                flow='substrate' if coeff < 0 else 'product',
                coefficient=coeff,
            ))
    gm = nx.Graph()
    for n in nodes:
        gm.add_node(n['key'], role=n['role'])
    for e in edges:
        gm.add_edge(e['a'], e['b'])
    if not len(gm) or not nx.is_connected(gm):
        return None
    # Ports: primary metabolites with reactions outside the window (entity graph).
    ent_graph, _ = _primary_entity_graph(diagram)
    ports = []
    for k in keys - set(rs):
        glyph = diagram['glyphs'].get(k) or {}
        ent = str(glyph.get('entity') or '').split('@')[0]
        if not ent:
            continue
        outside = sorted(
            n for n in ent_graph['e:' + ent]
            if n in rx and n not in rs
        )
        if outside:
            ports.append(dict(key=k, outside_reactions=outside, scope='available_source_map'))
    return dict(
        id=uid([diagram['id'], rs]),
        source='escher',
        species=diagram.get('species', 'unspecified'),
        lineage=diagram['lineage'],
        map_id=diagram['id'],
        name=diagram['name'],
        license=diagram['license'],
        quality='semantic_transfer',
        size=len(rs),
        nodes=nodes,
        edges=edges,
        reactions=selected,
        attachments=[],
        ports=ports,
        source_sha256=diagram['source_sha256'],
        review_status='not_required',
        original_instance=str(path),
        biochemical_visibility='complete reference reaction; partial observed drawing',
        layout_origin=dict(LAYOUT_ORIGIN_EXPERT),
        atlas_schema_version=SCHEMA_VERSION,
        extraction=extraction_block(parameters, map_counts),
    )


def _species_for_path(path):
    """Best-effort species from map_observations, else from the model folder name."""
    path = Path(path)
    try:
        from core import source_maps
        digest = sha(path.read_bytes())
        for m in source_maps():
            if m.get('source') == 'escher' and m.get('source_sha256') == digest and m.get('species'):
                return m['species']
    except Exception:
        pass
    folder = path.parent.name
    return {
        'e_coli_core': 'Escherichia coli',
        'iAF1260': 'Escherichia coli',
        'iJO1366': 'Escherichia coli',
        'iMM904': 'Saccharomyces cerevisiae',
        'RECON1': 'Homo sapiens',
    }.get(folder, 'unspecified')


def cut_map(path, *, species=None, parameters=None):
    """Cut one Escher JSON map into Atlas modules."""
    path = Path(path)
    parameters = dict(parameters or default_parameters())
    species = species or _species_for_path(path)
    diagram = escher_diagram(path, source='escher', license='CC0-1.0', species=species)
    windows = reaction_windows(
        diagram,
        plus_full_component=bool(parameters.get('plus_full_component', True)),
        max_reactions=int(parameters.get('max_reactions', MAX_REACTIONS)),
    )
    source_events = len(diagram['events'])
    id_of = {r['key']: r['id'] for r in diagram['events']}
    covered = {id_of[k] for rs in windows for k in rs if k in id_of}
    event_ids = set(id_of.values())
    coverage = (len(covered & event_ids) / len(event_ids)) if event_ids else 0.0
    map_counts = dict(
        source_events=source_events,
        modules=0,  # filled after emit
        coverage=round(coverage, 6),
    )
    modules = []
    for rs in sorted(windows, key=lambda t: (len(t), t)):
        m = _module_from_window(diagram, rs, path=path, parameters=parameters, map_counts=map_counts)
        if m is not None:
            modules.append(m)
    map_counts['modules'] = len(modules)
    for m in modules:
        m['extraction']['counts'] = dict(map_counts)
    # Prefer stable lineage from map_observations when the same file hash is known.
    lineage = diagram['lineage']
    map_id = diagram['id']
    try:
        from core import source_maps
        for prev in source_maps():
            if prev.get('source') == 'escher' and prev.get('source_sha256') == diagram['source_sha256']:
                lineage = prev.get('lineage') or lineage
                map_id = prev.get('id') or map_id
                break
    except Exception:
        pass
    for m in modules:
        m['lineage'] = lineage
        m['map_id'] = map_id
    map_record = dict(
        id=map_id,
        source='escher',
        species=species,
        lineage=lineage,
        name=diagram['name'],
        license=diagram['license'],
        source_sha256=diagram['source_sha256'],
        source_path=str(path),
        map_id=diagram.get('map_id'),
        provenance_status='manual_history_unverified',
        layout_origin=dict(LAYOUT_ORIGIN_EXPERT),
        atlas_schema_version=SCHEMA_VERSION,
        extraction=extraction_block(parameters, map_counts),
    )
    return dict(diagram=diagram, map=map_record, modules=modules, windows=windows, counts=map_counts)


def cut_all(paths=None, *, parameters=None):
    """Cut every community Escher map under ``raw/escher``."""
    parameters = dict(parameters or default_parameters())
    paths = list(paths or raw_escher_paths())
    maps, modules = [], []
    per_map = []
    for path in paths:
        result = cut_map(path, parameters=parameters)
        maps.append(result['map'])
        modules.extend(result['modules'])
        per_map.append(dict(
            path=str(path),
            source_sha256=result['diagram']['source_sha256'],
            name=result['diagram']['name'],
            **result['counts'],
        ))
    return dict(maps=maps, modules=modules, per_map=per_map, parameters=parameters)


def stored_escher_reaction_sets(observations=None, *, max_reactions=MAX_REACTIONS, max_nodes=MAX_NODES):
    """Reaction-id frozensets from stored Escher observations that pass size caps."""
    from core import source_modules
    obs = observations if observations is not None else source_modules()
    by_sha = {}
    for o in obs:
        if o.get('source') != 'escher':
            continue
        n_nodes = len(o.get('core_nodes', o.get('nodes') or []))
        size = o['size'] if o.get('quality') == 'semantic_transfer' else len(o.get('edges') or [])
        if not (1 <= size <= max_reactions and n_nodes <= max_nodes):
            continue
        key = frozenset(r['id'] for r in o.get('reactions') or [])
        by_sha.setdefault(o['source_sha256'], set()).add(key)
    return by_sha


def parity_report(paths=None, *, parameters=None):
    """Compare cutter reaction sets to stored observations (per map and overall)."""
    parameters = dict(parameters or default_parameters())
    paths = list(paths or raw_escher_paths())
    stored = stored_escher_reaction_sets(
        max_reactions=int(parameters.get('max_reactions', MAX_REACTIONS)),
        max_nodes=int(parameters.get('max_nodes', MAX_NODES)),
    )
    rows = []
    total_inter = total_stored = total_new = 0
    for path in paths:
        result = cut_map(path, parameters=parameters)
        digest = result['diagram']['source_sha256']
        new_sets = {frozenset(r['id'] for r in m['reactions']) for m in result['modules']}
        old_sets = stored.get(digest, set())
        inter = len(new_sets & old_sets)
        only_old = len(old_sets - new_sets)
        only_new = len(new_sets - old_sets)
        total_inter += inter
        total_stored += len(old_sets)
        total_new += len(new_sets)
        # reaction coverage on the raw map
        raw = __import__('json').loads(Path(path).read_bytes())[1]
        all_ids = {r.get('bigg_id') or rid for rid, r in raw.get('reactions', {}).items()}
        cov_new = set()
        for s in new_sets:
            cov_new |= s
        cov_old = set()
        for s in old_sets:
            cov_old |= s
        size_ok = all(1 <= m['size'] <= MAX_REACTIONS and len(m['nodes']) <= MAX_NODES for m in result['modules'])
        rows.append(dict(
            path=str(path),
            name=result['diagram']['name'],
            stored=len(old_sets),
            generated=len(new_sets),
            intersection=inter,
            only_stored=only_old,
            only_generated=only_new,
            recall=(inter / len(old_sets)) if old_sets else 1.0,
            reaction_coverage_stored=(len(cov_old & all_ids) / len(all_ids)) if all_ids else 0.0,
            reaction_coverage_generated=(len(cov_new & all_ids) / len(all_ids)) if all_ids else 0.0,
            size_bounds_ok=size_ok,
            source_events=result['counts']['source_events'],
        ))
    summary = dict(
        maps=len(rows),
        stored_modules=total_stored,
        generated_modules=total_new,
        intersection=total_inter,
        recall=(total_inter / total_stored) if total_stored else 1.0,
        precision=(total_inter / total_new) if total_new else 1.0,
        min_generated_coverage=min((r['reaction_coverage_generated'] for r in rows), default=1.0),
        min_stored_coverage=min((r['reaction_coverage_stored'] for r in rows), default=1.0),
        all_size_bounds_ok=all(r['size_bounds_ok'] for r in rows),
        parameters=parameters,
        method=METHOD,
        version=VERSION,
    )
    return dict(summary=summary, maps=rows)
