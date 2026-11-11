"""Complete WikiPathways refresh from all GPML files of a pinned release.

Generic pipeline (no species/map special cases):
  - parse_gpml → map records for every GPML
  - validated Rhea bindings → wp_diagram + diagram_modules (semantic_transfer)
  - otherwise metabolite–metabolite drawn interactions → geometry_only windows
    via the same diagram_modules BFS cutter
  - rebuild open_core WP maps/modules/reactions; refresh observations + manifests

Escher and Reactome rows in open_core are preserved.
"""
from __future__ import annotations

import gzip
import json
import sqlite3
import zlib
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

import numpy as np
import networkx as nx

from core import DATA, INPUTS, OUT, dump, read, sha, uid, normalize
from semantics import rhea_index, parse_gpml
from database import (
    pack, unpack, graph, mappings, geometry_error, diagram_modules, wp_diagram,
)
from atlas_registry import SCHEMA_VERSION, ensure_schema, normalize_open_core, ATTRIBUTION, LAYOUT_ORIGIN_EXPERT

RELEASE = '20261110'
OPEN_CORE = DATA / 'open_core.sqlite'
ATLAS = DATA / 'atlas.sqlite'


def _metabolite_glyphs(d):
    out = {}
    for k, g in d['glyphs'].items():
        if g.get('role') != 'Metabolite':
            continue
        if g.get('chebi'):
            entity = 'chebi:' + g['chebi']
        elif isinstance(g.get('entity'), dict) and g['entity'].get('ID'):
            db = g['entity'].get('Database') or 'Unknown'
            entity = f"{db}:{g['entity']['ID']}"
        else:
            entity = g.get('label') or k
        out[k] = dict(entity=entity, xy=list(g['xy']), label=g.get('label', ''), kind='Metabolite')
    return out


def geometry_diagram(d):
    """Diagram form for geometry-only cutting: pseudo-reactions on metabolite–metabolite edges."""
    glyphs = _metabolite_glyphs(d)
    events = []
    for e in d['edges']:
        a, b = e.get('source'), e.get('target')
        if a not in glyphs or b not in glyphs or a == b:
            continue
        xy = ((np.array(glyphs[a]['xy']) + np.array(glyphs[b]['xy'])) / 2).tolist()
        events.append(dict(
            key='r:' + e['key'], id='route:' + e['key'], references={}, coefficients={},
            visible={a: -1, b: 1}, xy=xy, position_origin='derived_from_observed_path',
            reversible=True, core=[a, b],
        ))
    return dict(
        id=d['id'], source=d['source'], species=d['species'], lineage=d['lineage'],
        map_id=d['id'], name=d['name'], license=d['license'], glyphs=glyphs, events=events,
        source_sha256=d['source_sha256'], path=d.get('path'), format='gpml', quality='geometry_only',
    )


def geometry_modules(d, max_roots=24):
    """BFS windows over geometry diagram; same size bounds as database.build legacy filter."""
    diagram = geometry_diagram(d)
    if not diagram['events']:
        return []
    mods = diagram_modules(diagram, max_roots=max_roots)
    out = []
    for m in mods:
        m['quality'] = 'geometry_only'
        m['biochemical_visibility'] = 'geometry_only; drawn interactions without reaction identity'
        m['reactions'] = []  # no chemical identity to index
        if 1 <= m['size'] <= 40 and len(m['nodes']) <= 100:
            out.append(m)
    # Re-stamp after size filter so counts match stored modules.
    from database import attach_extraction, WINDOW_SIZES
    attach_extraction(
        out, 'geometry_modules',
        dict(max_roots=max_roots, window_sizes=list(WINDOW_SIZES),
             underlying='diagram_modules over metabolite–metabolite GPML edges'),
        source_events=len(diagram['events']),
    )
    return out


def map_record(d, version=RELEASE, provenance_status='resolved_pathway_id'):
    return dict(
        id=d['id'], source='wikipathways', source_id=d.get('source_id') or d.get('path'),
        species=d['species'], lineage=d['lineage'], name=d['name'], author=d.get('author', ''),
        license=d.get('license') or 'CC0-1.0', source_sha256=d['source_sha256'],
        provenance_status=provenance_status, version=version,
        url=f'https://data.wikipathways.org/{version}/gpml/',
        attribution=ATTRIBUTION['wikipathways'],
        layout_origin=dict(LAYOUT_ORIGIN_EXPERT),
        atlas_schema_version=SCHEMA_VERSION,
        format='gpml',
    )


def map_observation(d, version=RELEASE, quality='geometry_only', use='geometry_pending_semantics'):
    """Observation payload compatible with units_library.geometry_maps / source_maps()."""
    glyphs = {}
    for k, g in d['glyphs'].items():
        glyphs[k] = dict(
            key=k, kind=g.get('role') or 'DataNode', role=g.get('role'),
            label=g.get('label', ''), xy=list(g['xy']),
            entity=g.get('entity'), chebi=g.get('chebi'),
            position_origin=g.get('position_origin', 'observed'),
        )
    edges = []
    for e in d['edges']:
        edges.append(dict(
            key=e['key'], source=e.get('source'), target=e.get('target'),
            points=e.get('points'), xref=e.get('xref'), arrows=e.get('arrows'),
            dashed=e.get('dashed', False),
        ))
    lineage = d['lineage']
    status = 'resolved_pathway_id' if str(lineage).startswith('WP') else 'import_or_propagation_unresolved'
    return dict(
        id=d['id'], source='wikipathways', source_id=d.get('source_id'), name=d['name'],
        species=d['species'], author=d.get('author', ''), license='CC0-1.0',
        lineage=lineage, lineage_hint=lineage if str(lineage).startswith('WP') else None,
        source_sha256=d['source_sha256'], provenance_status=status, version=version,
        format='gpml', quality=quality, use=use, glyphs=glyphs, edges=edges,
        metabolite_glyphs=sum(1 for g in glyphs.values() if g.get('kind') == 'Metabolite'),
        reactions=[], comments=[], issues=[], metadata={}, family=d.get('species'),
    )


def extract_all_gpml(rhea=None):
    rhea = rhea or rhea_index()
    files = sorted((INPUTS / 'raw/wp').rglob('*.gpml'))
    maps, semantic_mods, geometry_mods = [], [], []
    map_obs, binding_audit = [], []
    stats = Counter()
    for i, path in enumerate(files):
        d = parse_gpml(path, rhea)
        d['path'] = str(path.relative_to(DATA.parent)) if path.is_absolute() else str(path)
        # Prefer repo-relative path under data/
        try:
            d['path'] = str(path.resolve().relative_to(DATA.parent.resolve()))
        except Exception:
            d['path'] = str(path)
        stats['gpml_files'] += 1
        rec = map_record(d)
        maps.append(rec)
        dd = wp_diagram(d, rhea)
        n_events = len(dd.get('events') or [])
        for b in d['bindings']:
            binding_audit.append(dict(
                map_id=d['id'], species=d['species'], interaction=b['interaction'],
                rhea_id=b['rhea_id'], master=b['master'], status=b['status'],
            ))
            stats[b['status']] += 1
        if n_events:
            mods = diagram_modules(dd)
            for m in mods:
                m['original_instance'] = d['path']
                m['layout_origin'] = dict(LAYOUT_ORIGIN_EXPERT)
                m['attribution'] = ATTRIBUTION['wikipathways']
                m['atlas_schema_version'] = SCHEMA_VERSION
            semantic_mods.extend(mods)
            stats['semantic_modules'] += len(mods)
            stats['semantic_maps'] += 1
            map_obs.append(map_observation(d, quality='reference', use='semantic_transfer'))
        else:
            mods = geometry_modules(d)
            for m in mods:
                m['original_instance'] = d['path']
                m['layout_origin'] = dict(LAYOUT_ORIGIN_EXPERT)
                m['attribution'] = ATTRIBUTION['wikipathways']
                m['atlas_schema_version'] = SCHEMA_VERSION
            geometry_mods.extend(mods)
            stats['geometry_modules'] += len(mods)
            if mods:
                stats['geometry_maps_with_modules'] += 1
            else:
                stats['maps_without_modules'] += 1
            map_obs.append(map_observation(d, quality='geometry_only', use='geometry_pending_semantics'))
        if (i + 1) % 200 == 0:
            print(f'extract {i+1}/{len(files)} semantic_mods={stats["semantic_modules"]} geo_mods={stats["geometry_modules"]}', flush=True)
    stats['maps'] = len(maps)
    return dict(
        maps=maps, semantic_mods=semantic_mods, geometry_mods=geometry_mods,
        map_obs=map_obs, binding_audit=binding_audit, stats=dict(stats), files=len(files),
    )


def _load_non_wp(db):
    maps, mods = [], []
    for mid, source, species, lineage, blob in db.execute('SELECT id, source, species, lineage, payload FROM maps'):
        if source == 'wikipathways':
            continue
        maps.append(unpack(blob))
    for row in db.execute(
        'SELECT id, source, species, lineage, quality, size, core_family, variant, payload FROM modules'
    ):
        if row[1] == 'wikipathways':
            continue
        mods.append(unpack(row[-1]))
    return maps, mods


def _classify_and_write(db, maps, mods):
    """Assign core_family/variant and write modules/maps/reactions (full replace)."""
    db.execute('DELETE FROM reactions')
    db.execute('DELETE FROM modules')
    db.execute('DELETE FROM maps')
    families = defaultdict(list)
    familyrows = []
    excluded = []
    for i, m in enumerate(sorted(mods, key=lambda m: m['id'])):
        if i % 2000 == 0:
            print(f'Index {i}/{len(mods)}', flush=True)
        original_xy = {n['key']: n['xy'] for n in m['nodes']}
        xy = normalize(original_xy)
        for n in m['nodes']:
            n['xy'] = xy[n['key']]
        m['core_coordinate_normalization'] = 'centered RMS radius; original full instance retained by immutable source pointer'
        gm = graph(m)
        if not len(gm) or not nx.is_connected(gm):
            excluded.append(dict(
                id=m['id'], source=m['source'], map_id=m.get('map_id'),
                reason='empty_or_disconnected_visible_core',
                original_instance=m.get('legacy_id', m.get('original_instance')),
            ))
            continue
        h = nx.weisfeiler_lehman_graph_hash(gm, node_attr='role')
        bucket = (len(gm), gm.number_of_edges(), h)
        family = None
        for f in families[bucket]:
            matches = mappings(gm, f['graph'], 1)
            if matches:
                family = f
                break
        if family is None:
            family = dict(id=uid(['core-v03', bucket, m['id']]), graph=gm, representative=m, variants=[], instances=[])
            families[bucket].append(family)
            familyrows.append(family)
        variant = None
        for v in family['variants'][:24]:
            matches = mappings(gm, v['graph'], 8)
            if matches and min(geometry_error(m, v['representative'], x) for x in matches) < .25:
                variant = v
                break
        if variant is None:
            variant = dict(id=uid([family['id'], m['id']]), graph=gm, representative=m, instances=[])
            family['variants'].append(variant)
        family['instances'].append(m['id'])
        variant['instances'].append(m['id'])
        m['core_family'] = family['id']
        m['variant'] = variant['id']
        m['biochemical_structure_signature'] = uid([(e['a'], e['b'], e.get('flow')) for e in m['edges']])
        pts = np.array([n['xy'] for n in m['nodes']], float)
        m['geometry_fingerprint'] = uid(np.round(np.sort(np.linalg.norm(pts[:, None] - pts[None, :], axis=2).ravel()), 3).tolist())
        m.setdefault('layout_origin', dict(LAYOUT_ORIGIN_EXPERT))
        m.setdefault('atlas_schema_version', SCHEMA_VERSION)
        db.execute(
            'INSERT OR REPLACE INTO modules VALUES(?,?,?,?,?,?,?,?,?)',
            (m['id'], m['source'], m['species'], m['lineage'], m['quality'], m['size'], m['core_family'], m['variant'], pack(m)),
        )
        for r in m.get('reactions') or []:
            for ns, values in (r.get('references') or {}).items():
                for val in values:
                    db.execute('INSERT OR IGNORE INTO reactions VALUES(?,?,?,?)', (m['id'], r['key'], ns, str(val)))
    seen = set()
    for d in maps:
        if d['id'] in seen:
            continue
        seen.add(d['id'])
        d.setdefault('atlas_schema_version', SCHEMA_VERSION)
        db.execute(
            'INSERT OR REPLACE INTO maps VALUES(?,?,?,?,?)',
            (d['id'], d['source'], d['species'], d['lineage'], pack(d)),
        )
    db.commit()
    return dict(excluded=excluded, families=len(familyrows), variants=sum(len(f['variants']) for f in familyrows))


def _geometry_core_edges(p):
    """Recover metabolite–metabolite pairs from route-expanded module edges."""
    route_neighbors = defaultdict(list)
    for e in p.get('edges') or []:
        a, b = e['a'], e['b']
        if str(a).startswith('route:'):
            route_neighbors[a].append(b)
        elif str(b).startswith('route:'):
            route_neighbors[b].append(a)
    pairs = []
    for rk, nb in sorted(route_neighbors.items()):
        mets = [x for x in nb if not str(x).startswith('route:')]
        if len(mets) >= 2:
            pairs.append([mets[0], mets[1]])
    return pairs


def export_module_observations(db):
    """Rebuild module_observations.json.gz: keep non-WP; write WP geometry_only for from_legacy."""
    path = INPUTS / 'module_observations.json.gz'
    with gzip.open(path, 'rt') as f:
        old = json.load(f)
    kept = [r for r in old if (r.get('observation') or {}).get('source') != 'wikipathways']
    new = []
    for mid, quality, cf, var, blob in db.execute(
        "SELECT id, quality, core_family, variant, payload FROM modules WHERE source='wikipathways'"
    ):
        p = unpack(blob)
        if quality == 'geometry_only':
            pairs = _geometry_core_edges(p)
            nodes = [
                dict(key=n['key'], role='geometric_metabolite', xy=n['xy'], entity=n.get('entity'), label=n.get('label', ''))
                for n in p.get('nodes', [])
                if n.get('role') == 'metabolite' and not str(n.get('key', '')).startswith('route:')
            ]
            obs = dict(
                id=p['id'], map_id=p.get('map_id'), source='wikipathways', species=p['species'],
                lineage=p['lineage'], quality='geometry_only', review_status='not_required',
                size=len(pairs), name=p.get('name'), license=p.get('license') or 'CC0-1.0',
                source_sha256=p.get('source_sha256'), nodes=nodes, edges=pairs, reactions=[],
                attachments=[], ports=p.get('ports') or [], all_participant_ports=p.get('ports') or [],
                core_nodes=nodes, core_edges=pairs, version=RELEASE,
            )
        else:
            # Semantic WP modules are rebuilt from GPML in database.build; store payload-aligned observation.
            obs = dict(
                id=p['id'], map_id=p.get('map_id'), source='wikipathways', species=p['species'],
                lineage=p['lineage'], quality='semantic_transfer', review_status='not_required',
                size=p['size'], name=p.get('name'), license=p.get('license') or 'CC0-1.0',
                source_sha256=p.get('source_sha256'), nodes=p.get('nodes') or [],
                edges=[[e['a'], e['b']] for e in p.get('edges') or []],
                reactions=[{k: r.get(k) for k in ('key', 'id', 'references', 'coefficients', 'xy', 'visible', 'position_origin') if k in r}
                           for r in p.get('reactions') or []],
                attachments=[], ports=p.get('ports') or [], all_participant_ports=p.get('ports') or [],
                core_nodes=p.get('nodes') or [], core_edges=[[e['a'], e['b']] for e in p.get('edges') or []],
                version=RELEASE,
            )
        new.append(dict(id=mid, family=cf, variant=var, observation=obs))
    out = kept + new
    with gzip.open(path, 'wt') as f:
        json.dump(out, f, ensure_ascii=False, separators=(',', ':'))
    return dict(kept=len(kept), wikipathways=len(new), total=len(out))

def export_map_observations(map_obs):
    path = INPUTS / 'map_observations.json.gz'
    with gzip.open(path, 'rt') as f:
        old = json.load(f)
    kept = [m for m in old if m.get('source') != 'wikipathways']
    out = kept + map_obs
    with gzip.open(path, 'wt') as f:
        json.dump(out, f, ensure_ascii=False, separators=(',', ':'))
    return dict(kept=len(kept), wikipathways=len(map_obs), total=len(out))


def refresh_source_registry(version=RELEASE):
    reg_path = INPUTS / 'source_registry.json'
    reg = read(reg_path)
    files = [f for f in reg.get('files', []) if f.get('source') != 'wikipathways']
    for p in sorted((INPUTS / 'raw/wp').rglob('*.gpml')):
        rel = str(p.relative_to(INPUTS))
        digest = sha(p.read_bytes())
        lineage = None
        import re
        m = re.search(r'WP\d+', p.name)
        files.append(dict(
            source='wikipathways', id=p.name, format='gpml',
            family=p.parent.name, species=p.parent.name.replace('_', ' '),
            path=rel, sha256=digest, license='CC0-1.0', use='core',
            url=f'https://data.wikipathways.org/{version}/gpml/',
            version=version, lineage=m[0] if m else None,
        ))
    # archive provenances already updated by prior refresh; keep non-wp + rewrite wp gpml entries
    reg['files'] = files
    reg['wikipathways_release'] = version
    reg['access_date'] = str(date.today())
    dump(reg_path, reg)
    return len(files)


def _write_geometry_extension(geometry_mods):
    """Ship geometry-only modules as a default extension (keeps open_core under GitHub size limits)."""
    ext = DATA / 'extensions' / 'wikipathways_geometry'
    ext.mkdir(parents=True, exist_ok=True)
    for m in geometry_mods:
        m['partition'] = 'open_core'
        m['admitted_to_open_core'] = True
        m['extension'] = 'wikipathways_geometry'
        m.setdefault('layout_origin', dict(LAYOUT_ORIGIN_EXPERT))
        m.setdefault('attribution', ATTRIBUTION['wikipathways'])
        m.setdefault('atlas_schema_version', SCHEMA_VERSION)
    with gzip.open(ext / 'modules.json.gz', 'wt') as f:
        json.dump(geometry_mods, f, ensure_ascii=False, separators=(',', ':'))
    dump(ext / 'SOURCE.json', dict(
        source='wikipathways', release=RELEASE, license='CC0-1.0',
        license_url='https://www.wikipathways.org/terms.html',
        modules=len(geometry_mods), quality='geometry_only',
        extraction='diagram_modules BFS windows over metabolite–metabolite GPML edges',
        admitted_to_open_core=True, partition='open_core',
        note='Default extension; loaded by database.modules(). Kept out of open_core.sqlite for size.',
        date=str(date.today()),
    ))
    return dict(path='data/extensions/wikipathways_geometry/modules.json.gz', modules=len(geometry_mods))


def rebuild():
    print('Extracting all WikiPathways GPML…', flush=True)
    extracted = extract_all_gpml()
    print('extract stats', extracted['stats'], flush=True)

    # Classify geometry for family/variant, then store as extension; semantic stays in sqlite.
    print('Writing geometry-only extension…', flush=True)
    # Assign provisional family/variant via a throwaway classify on geometry alone for observation export
    geo_ext = _write_geometry_extension(extracted['geometry_mods'])

    db = sqlite3.connect(OPEN_CORE)
    ensure_schema(db)
    non_wp_maps, non_wp_mods = _load_non_wp(db)
    print(f'preserving non-WP maps={len(non_wp_maps)} modules={len(non_wp_mods)}', flush=True)

    all_maps = non_wp_maps + extracted['maps']
    # open_core holds semantic WP only (+ escher/reactome); geometry is the extension
    all_mods = non_wp_mods + extracted['semantic_mods']
    print(f'classify+write total_mods={len(all_mods)} maps={len(all_maps)}', flush=True)
    index_stats = _classify_and_write(db, all_maps, all_mods)
    print('index', index_stats, flush=True)

    maps = dict(db.execute('SELECT source, count(*) FROM maps GROUP BY source'))
    mods = dict(db.execute('SELECT source, count(*) FROM modules GROUP BY source'))
    quals = dict(db.execute("SELECT quality, count(*) FROM modules WHERE source='wikipathways' GROUP BY quality"))
    nrx = db.execute('SELECT count(*) FROM reactions').fetchone()[0]
    wp_rx = db.execute(
        "SELECT count(*) FROM reactions r JOIN modules m ON m.id=r.module_id WHERE m.source='wikipathways'"
    ).fetchone()[0]

    print('Exporting observations…', flush=True)
    mob = export_module_observations(db)
    # Also append geometry modules into module_observations for from_legacy rebuilds
    geo_obs = []
    for m in extracted['geometry_mods']:
        pairs = _geometry_core_edges(m) if m.get('edges') and any(str(e.get('a','')).startswith('route:') or str(e.get('b','')).startswith('route:') for e in m['edges']) else []
        if not pairs:
            pairs = _geometry_core_edges(m)
        nodes = [
            dict(key=n['key'], role='geometric_metabolite', xy=n['xy'], entity=n.get('entity'), label=n.get('label', ''))
            for n in m.get('nodes', [])
            if n.get('role') == 'metabolite' and not str(n.get('key', '')).startswith('route:')
        ]
        obs = dict(
            id=m['id'], map_id=m.get('map_id'), source='wikipathways', species=m['species'],
            lineage=m['lineage'], quality='geometry_only', review_status='not_required',
            size=len(pairs) or m.get('size', 0), name=m.get('name'), license=m.get('license') or 'CC0-1.0',
            source_sha256=m.get('source_sha256'), nodes=nodes, edges=pairs, reactions=[],
            attachments=[], ports=m.get('ports') or [], all_participant_ports=m.get('ports') or [],
            core_nodes=nodes, core_edges=pairs, version=RELEASE,
        )
        geo_obs.append(dict(id=m['id'], family=m.get('core_family'), variant=m.get('variant'), observation=obs))
    # merge geometry into module_observations (export_module_observations already wrote semantic WP)
    path = INPUTS / 'module_observations.json.gz'
    with gzip.open(path, 'rt') as f:
        cur = json.load(f)
    # drop any prior WP geometry rows then add fresh
    cur = [r for r in cur if not (
        (r.get('observation') or {}).get('source') == 'wikipathways'
        and (r.get('observation') or {}).get('quality') == 'geometry_only'
    )]
    cur.extend(geo_obs)
    with gzip.open(path, 'wt') as f:
        json.dump(cur, f, ensure_ascii=False, separators=(',', ':'))
    mob = dict(mob, geometry_observations=len(geo_obs), total=len(cur))

    mpb = export_map_observations(extracted['map_obs'])
    nreg = refresh_source_registry()
    db.close()

    print('Normalizing provenance + registry…', flush=True)
    norm = normalize_open_core(OPEN_CORE)

    refresh = dict(
        from_version='20261010',
        to_version=RELEASE,
        complete_rebuild=True,
        gpml_files=extracted['files'],
        extract_stats=extracted['stats'],
        maps_by_source=maps,
        modules_by_source=mods,
        wikipathways_quality_in_open_core=quals,
        wikipathways_geometry_extension=geo_ext,
        reactions_total=nrx,
        wikipathways_reactions=wp_rx,
        index=index_stats,
        module_observations=mob,
        map_observations=mpb,
        source_registry_files=nreg,
        normalize=norm,
        maps_without_modules_note='GPML with neither validated Rhea bindings nor metabolite–metabolite edges (typically gene/signaling maps); not metabolic layout modules',
        licence='CC0-1.0 (WikiPathways content; GPML files have no per-file License attribute; project terms CC0)',
        date=str(date.today()),
        binding_status={k: extracted['stats'].get(k, 0) for k in (
            'validated_visible_participants', 'reaction_id_only', 'reference_direction_disagrees',
            'visual_direction_unresolved', 'obsolete_reference_not_promoted', 'unresolved', 'aggregate_not_promoted',
        )},
    )
    dump(INPUTS / 'wikipathways_refresh.json', refresh)
    dump(OUT / 'binding_audit_wp.json', extracted['binding_audit'])
    dump(OUT / 'wikipathways_rebuild_summary.json', refresh)
    if ATLAS.exists():
        ATLAS.write_bytes(OPEN_CORE.read_bytes())
    print('DONE', json.dumps({k: refresh[k] for k in (
        'maps_by_source', 'modules_by_source', 'wikipathways_quality_in_open_core',
        'wikipathways_geometry_extension', 'reactions_total', 'wikipathways_reactions',
    )}, indent=2), flush=True)
    return refresh

if __name__ == '__main__':
    rebuild()
