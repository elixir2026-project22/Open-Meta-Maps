"""Human-maps extension builder (official escher/ JSON only; CC BY 4.0).

Downloads missing maps via ``source_registry.json`` remote URLs (or the local
``data/source_inputs/archives/human_maps.zip`` / ``data/references/human_maps``
cache), cuts windows with Stage B ``database.diagram_modules``, classifies
families, remaps MAR → BiGG/MetaNetX/Rhea via ``canonical_ids``, and writes
``data/extensions/human_maps/``.
"""
from __future__ import annotations

import gzip
import json
import sqlite3
import time
import urllib.request
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

import networkx as nx
import numpy as np

import database
from atlas_registry import LAYOUT_ORIGIN_EXPERT, SCHEMA_VERSION, ATTRIBUTION
from canonical_ids import load_human_gem_crosswalk, remap_human_maps_module
from core import DATA, INPUTS, ROOT, dump, read, sha, uid

SOURCE = 'human_maps'
LICENSE = 'CC-BY-4.0'
LICENSE_URI = 'https://creativecommons.org/licenses/by/4.0/'
REPO = 'https://github.com/SysBioChalmers/Human-maps'
EXT = DATA / 'extensions' / SOURCE
RAW = DATA / 'references' / SOURCE
UA = 'OpenPathwayMaps-research/0.3'


def _registry_files():
    reg = read(INPUTS / 'source_registry.json')
    return [f for f in reg.get('files', []) if f.get('source') == SOURCE]


def _dest_for(rel_id: str) -> Path:
    # registry id looks like 'compartment/escher/cytosol_1.json'
    return RAW / Path(rel_id)


def ensure_local(files=None, timeout=40):
    """Materialise official Escher JSON under data/references/human_maps/."""
    files = files or _registry_files()
    archive = INPUTS / 'archives' / 'human_maps.zip'
    fetched, cached, failed = [], [], []
    # Prefer zip when present
    if archive.exists():
        import zipfile
        with zipfile.ZipFile(archive) as z:
            by_name = {Path(n).name: n for n in z.namelist() if '/escher/' in n and n.endswith('.json')}
            for f in files:
                rel = f['id']
                dest = _dest_for(rel)
                if dest.exists() and dest.stat().st_size > 0:
                    cached.append(rel)
                    continue
                name = Path(rel).name
                if name not in by_name:
                    failed.append(dict(id=rel, error='missing_in_archive'))
                    continue
                dest.parent.mkdir(parents=True, exist_ok=True)
                raw = z.read(by_name[name])
                dest.write_bytes(raw)
                fetched.append(rel)
        return dict(fetched=fetched, cached=cached, failed=failed, via='archive')

    for f in files:
        rel = f['id']
        dest = _dest_for(rel)
        if dest.exists() and dest.stat().st_size > 0:
            cached.append(rel)
            continue
        url = f.get('remote_url') or f'{REPO}/raw/master/{rel}'
        try:
            req = urllib.request.Request(url, headers={'User-Agent': UA})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                raw = r.read()
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(raw)
            fetched.append(rel)
        except Exception as e:
            failed.append(dict(id=rel, url=url, error=str(e)))
    return dict(fetched=fetched, cached=cached, failed=failed, via='remote')


def classify(mods):
    families, family_rows = defaultdict(list), []
    kept, excluded = [], []
    for i, m in enumerate(mods):
        gm = database.graph(m)
        if not len(gm) or not nx.is_connected(gm):
            excluded.append(dict(id=m['id'], reason='empty_or_disconnected_visible_core'))
            continue
        h = nx.weisfeiler_lehman_graph_hash(gm, node_attr='role')
        bucket = (len(gm), gm.number_of_edges(), h)
        family = None
        for f in families[bucket]:
            if database.mappings(gm, f['graph'], 1):
                family = f
                break
        if family is None:
            family = dict(id=uid(['human-maps', bucket, m['id']]), graph=gm, representative=m, variants=[], instances=[])
            families[bucket].append(family)
            family_rows.append(family)
        variant = None
        for v in family['variants'][:24]:
            matches = database.mappings(gm, v['graph'], 8)
            if matches and min(database.geometry_error(m, v['representative'], x) for x in matches) < .25:
                variant = v
                break
        if variant is None:
            variant = dict(id=uid([family['id'], m['id']]), graph=gm, representative=m, instances=[])
            family['variants'].append(variant)
        family['instances'].append(m['id'])
        variant['instances'].append(m['id'])
        m['core_family'] = family['id']
        m['variant'] = variant['id']
        arr = np.array([n['xy'] for n in m['nodes']])
        m['geometry_fingerprint'] = uid(
            np.round(np.sort(np.linalg.norm(arr[:, None] - arr[None, :], axis=2).ravel()), 3).tolist()
        )
        kept.append(m)
        if i and i % 500 == 0:
            print('human_maps index', i, '/', len(mods), flush=True)
    return kept, family_rows, excluded


def build(access_date=None, version=None):
    EXT.mkdir(parents=True, exist_ok=True)
    RAW.mkdir(parents=True, exist_ok=True)
    access_date = access_date or str(date.today())
    version = version or access_date
    t0 = time.time()
    files = _registry_files()
    assert files, 'no human_maps rows in source_registry.json'
    materialise = ensure_local(files)
    print('human_maps materialise', {k: (len(v) if isinstance(v, list) else v) for k, v in materialise.items()}, flush=True)

    cite = dict(
        source=SOURCE, url=REPO, license=LICENSE, license_uri=LICENSE_URI,
        attribution=ATTRIBUTION.get(SOURCE) or 'Human-maps (CC BY 4.0); https://github.com/SysBioChalmers/Human-maps',
        access_date=access_date, version=version,
    )
    modules, map_records, per_map, failures = [], [], [], []
    for f in files:
        rel = f['id']
        path = _dest_for(rel)
        if not path.exists():
            failures.append(dict(id=rel, error='file_missing'))
            continue
        try:
            d = database.escher_diagram(
                path,
                source=SOURCE,
                license=LICENSE,
                species='Homo sapiens',
                reaction_namespace='human.gem',
            )
            d['lineage'] = f'{SOURCE}:{rel}'
            d['url'] = f.get('url') or REPO
            d['source_id'] = rel
            try:
                d['path'] = str(path.resolve().relative_to(ROOT.resolve()))
            except Exception:
                d['path'] = str(path)
        except Exception as e:
            failures.append(dict(id=rel, error=str(e)))
            continue
        mods = database.diagram_modules(d)
        for m in mods:
            m['partition'] = 'open_core'
            m['admitted_to_open_core'] = True
            m['layout_origin'] = dict(LAYOUT_ORIGIN_EXPERT)
            m['source_citation'] = cite
            m['attribution'] = cite['attribution']
            m['provenance_status'] = 'human_maps_escher_json'
            m['version'] = version
            m['url'] = d.get('url')
            m['source_id'] = rel
            m['lineage'] = d['lineage']
            m['map_id'] = d['id']
            m['original_instance'] = d.get('path')
            m['atlas_schema_version'] = SCHEMA_VERSION
        modules.extend(mods)
        map_records.append(dict(
            id=d['id'], source=SOURCE, source_id=rel, species=d['species'], lineage=d['lineage'],
            name=d['name'], license=LICENSE, source_sha256=d['source_sha256'],
            provenance_status='human_maps_escher_json', attribution=cite['attribution'],
            layout_origin=dict(LAYOUT_ORIGIN_EXPERT), version=version, url=d.get('url'),
            path=d.get('path'), format='escher', events=len(d['events']), glyphs=len(d['glyphs']),
            modules=len(mods), audit=d.get('audit'),
            extraction=(mods[0].get('extraction') if mods else None),
        ))
        per_map.append(dict(id=rel, events=len(d['events']), modules=len(mods), sha256=d['source_sha256']))

    kept, family_rows, excluded = classify(modules)
    crosswalk = load_human_gem_crosswalk()
    for m in kept:
        remap_human_maps_module(m, crosswalk)

    db = sqlite3.connect(EXT / 'atlas.sqlite')
    db.executescript(
        '''CREATE TABLE IF NOT EXISTS sources(id TEXT PRIMARY KEY, payload BLOB);
        CREATE TABLE IF NOT EXISTS maps(id TEXT PRIMARY KEY, payload BLOB);
        CREATE TABLE IF NOT EXISTS modules(id TEXT PRIMARY KEY, size INTEGER, core_family TEXT, variant TEXT, payload BLOB);
        CREATE TABLE IF NOT EXISTS reactions(module_id TEXT, reaction_key TEXT, namespace TEXT, identifier TEXT,
          PRIMARY KEY(module_id, reaction_key, namespace, identifier));
        DELETE FROM sources; DELETE FROM maps; DELETE FROM modules; DELETE FROM reactions;'''
    )
    db.execute('INSERT INTO sources VALUES(?,?)', (version, database.pack(cite)))
    for rec in map_records:
        db.execute('INSERT INTO maps VALUES(?,?)', (rec['id'], database.pack(rec)))
    for m in kept:
        db.execute(
            'INSERT INTO modules VALUES(?,?,?,?,?)',
            (m['id'], m['size'], m['core_family'], m['variant'], database.pack(m)),
        )
        for r in m.get('reactions') or []:
            for ns, values in (r.get('references') or {}).items():
                for val in values:
                    db.execute(
                        'INSERT OR IGNORE INTO reactions VALUES(?,?,?,?)',
                        (m['id'], r['key'], ns, str(val)),
                    )
    db.commit()
    db.close()

    with gzip.open(EXT / 'modules.json.gz', 'wt') as f:
        json.dump(kept, f, ensure_ascii=False)
    dump(EXT / 'maps.json', map_records)
    dump(EXT / 'SOURCE.json', dict(
        source=SOURCE, url=REPO, license=LICENSE, license_uri=LICENSE_URI,
        maps=len(map_records), modules=len(kept), families=len(family_rows),
        admitted_to_open_core=True, partition='open_core',
        import_policy='official escher/ JSON only via database.escher_diagram + diagram_modules; SVG/BY-SA excluded',
        by_sa_note=(
            'Companion SVG exports under the Human-maps repository embed Creative Commons '
            'Attribution-ShareAlike (CC BY-SA) metadata and are NOT imported.'
        ),
        access_date=access_date, version=version,
        license_scope='official escher/ JSON maps only',
        reaction_namespace='human.gem',
        canonicalization='Human-GEM reactions.tsv crosswalk → bigg.reaction / metanetx.reaction / rhea when annotated; MAR kept under human.gem',
        extraction=dict(
            methods=dict(Counter((m.get('extraction') or {}).get('method') for m in kept)),
            modules=len(kept), atlas_schema_version=SCHEMA_VERSION,
            parameters=dict(independent_only=True, root_policy='scaled_max_24_ceil_N_4_or_all_upto_48'),
        ),
        materialise=materialise, failures=failures, excluded_modules=len(excluded),
        seconds=round(time.time() - t0, 1),
    ))
    dump(EXT / 'import_audit.json', dict(
        maps=len(map_records), modules=len(kept), per_map=per_map,
        failures=failures, excluded=excluded[:200], zero_module_maps=[
            p for p in per_map if not p['modules']
        ],
    ))
    summary = dict(
        maps=len(map_records), modules=len(kept), families=len(family_rows),
        failures=len(failures), excluded=len(excluded), seconds=round(time.time() - t0, 1),
    )
    print('human_maps DONE', summary, flush=True)
    return summary


if __name__ == '__main__':
    build()
