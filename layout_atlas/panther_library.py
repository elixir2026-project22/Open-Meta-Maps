"""PANTHER Pathway CellDesigner maps (CC BY 4.0) via the generic CellDesigner importer.

Source: https://data.pantherdb.org/ftp/pathway/3.4/CD4.1/
Terms: https://pantherdb.org/tou.jsp
Licence: CC BY 4.0

Plain species labels are resolved through ``identity_labels`` (BiGG/ChEBI bridge
from public example models). Observed glyph centres and non-branch centre-lines
come from ``celldesigner.decode``; modules via ``database.diagram_modules``.
"""
from __future__ import annotations

import gzip
import json
import os
import sqlite3
import sys
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

from core import *
import database
from atlas_registry import LAYOUT_ORIGIN_EXPERT, ATTRIBUTION
import celldesigner

EXT = DATA / 'extensions/panther'
RESULTS = OUT / 'panther_extension'
RAW = INPUTS / 'raw/panther'
FTP = 'https://data.pantherdb.org/ftp/pathway/3.4/CD4.1/'
TOU = 'https://pantherdb.org/tou.jsp'
LICENSE = 'CC-BY-4.0'
LICENSE_URI = 'https://creativecommons.org/licenses/by/4.0/'
SOURCE = 'panther'
VERSION = '3.4/CD4.1'


def citation(version=VERSION, access_date=None):
    return dict(
        title='PANTHER Pathway CellDesigner maps',
        url=FTP,
        terms=TOU,
        license=dict(name='CC BY 4.0', uri=LICENSE_URI),
        version=version,
        retrieved=access_date or str(date.today()),
        attribution=(
            'PANTHER Pathway (https://pantherdb.org/), Creative Commons Attribution 4.0 International; '
            f'CellDesigner maps from {FTP}.'
        ),
        changes=(
            'Decoded with layout_atlas/celldesigner.py; plain labels resolved through '
            'identity_labels (BiGG/ChEBI bridge); module windows via database.diagram_modules.'
        ),
    )


def list_local():
    return sorted(RAW.glob('*.xml'))


def decode_all(version=VERSION, access_date=None):
    cite = citation(version, access_date)
    diagrams, failures = [], []
    for path in list_local():
        try:
            d = celldesigner.decode(
                path,
                source=SOURCE,
                license=LICENSE,
                name=path.stem.replace('_', ' '),
                species='PANTHER Pathway (multi-organism schematic)',
                lineage=f'{SOURCE}:{path.stem}',
                partition='open_core',
                admitted_to_open_core=True,
                source_id=path.stem,
                source_citation=cite,
                attribution=ATTRIBUTION.get(SOURCE) or cite['attribution'],
                layout_origin=dict(LAYOUT_ORIGIN_EXPERT),
                provenance_status='panther_celldesigner_cd4',
                url=FTP + path.name,
                version=version,
                terms=TOU,
            )
            diagrams.append(d)
        except Exception as e:
            failures.append(dict(path=str(path), error=str(e)))
    return diagrams, failures, cite


def classify(mods):
    families, family_rows = defaultdict(list), []
    kept, excluded = [], []
    for i, m in enumerate(mods):
        gm = database.graph(m)
        if not len(gm) or not nx.is_connected(gm):
            excluded.append(dict(id=m['id'], reason='empty_or_disconnected_visible_core'))
            continue
        bucket = (len(gm), gm.number_of_edges(), nx.weisfeiler_lehman_graph_hash(gm, node_attr='role'))
        family = None
        for f in families[bucket]:
            if database.mappings(gm, f['graph'], 1):
                family = f
                break
        if family is None:
            family = dict(id=uid(['panther-core', bucket, m['id']]), graph=gm, variants=[], instances=[])
            families[bucket].append(family)
            family_rows.append(family)
        variant = None
        for v in family['variants'][:24]:
            maps = database.mappings(gm, v['graph'], 8)
            if maps and min(database.geometry_error(m, v['representative'], mp) for mp in maps) < 0.25:
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
        if i % 200 == 0:
            print('panther index', i, '/', len(mods), flush=True)
    return kept, family_rows, excluded


def build(version=VERSION, access_date=None):
    EXT.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    access_date = access_date or str(date.today())
    assert list_local(), f'no PANTHER XML under {RAW}'
    diagrams, failures, cite = decode_all(version, access_date)
    modules, map_records, per_map = [], [], []
    for d in diagrams:
        mods = database.diagram_modules(d)
        for m in mods:
            m['partition'] = 'open_core'
            m['admitted_to_open_core'] = True
            m['layout_origin'] = dict(LAYOUT_ORIGIN_EXPERT)
            m['source_citation'] = cite
            m['attribution'] = cite['attribution']
            m['provenance_status'] = 'panther_celldesigner_cd4'
            m['version'] = version
            m['url'] = d.get('url')
            m['source_id'] = d.get('source_id')
            m['lineage'] = d['lineage']
            m['map_id'] = d['id']
        modules.extend(mods)
        map_records.append(dict(
            id=d['id'], source=SOURCE, source_id=d.get('source_id'), species=d['species'],
            lineage=d['lineage'], name=d['name'], license=LICENSE, source_sha256=d['source_sha256'],
            provenance_status='panther_celldesigner_cd4', attribution=cite['attribution'],
            layout_origin=dict(LAYOUT_ORIGIN_EXPERT), version=version, url=d.get('url'),
            path=d.get('path'), format='celldesigner', events=len(d['events']), glyphs=len(d['glyphs']),
            modules=len(mods), audit=d.get('audit'),
            extraction=(mods[0].get('extraction') if mods else database.extraction_block(
                'diagram_modules', dict(max_roots=24, window_sizes=[1, 2, 4, 8, 12, 24]),
                dict(source_events=len(d['events']), modules=0, reactions_covered=0, coverage=None),
            )),
        ))
        per_map.append(dict(source_id=d.get('source_id'), name=d['name'], sha256=d['source_sha256'],
                            events=len(d['events']), glyphs=len(d['glyphs']), modules=len(mods),
                            excluded=len(d.get('excluded') or [])))

    kept, family_rows, excluded = classify(modules)
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
        db.execute('INSERT INTO modules VALUES(?,?,?,?,?)', (m['id'], m['size'], m['core_family'], m['variant'], database.pack(m)))
        for r in m['reactions']:
            for ns, values in (r.get('references') or {}).items():
                for val in values:
                    db.execute('INSERT OR IGNORE INTO reactions VALUES(?,?,?,?)', (m['id'], r['key'], ns, str(val)))
    db.commit()
    db.close()
    with gzip.open(EXT / 'modules.json.gz', 'wt') as f:
        json.dump(kept, f, ensure_ascii=False)
    dump(EXT / 'maps.json', map_records)
    dump(EXT / 'families.json', dict(families=[
        dict(id=f['id'], instances=f['instances'], variants=[
            dict(id=v['id'], template=v['representative']['id'], instances=v['instances']) for v in f['variants']
        ]) for f in family_rows]))

    file_hashes = []
    for p in list_local():
        raw = p.read_bytes()
        file_hashes.append(dict(
            id=p.stem, source=SOURCE, path=str(p.relative_to(ROOT)).replace('\\', '/'),
            sha256=sha(raw), bytes=len(raw), license=LICENSE, license_uri=LICENSE_URI,
            version=version, url=FTP + p.name, format='celldesigner-sbml', access_date=access_date, terms=TOU,
        ))
    terms = DATA / 'terms/panther_tou.html'
    terms_meta = dict(path=str(terms.relative_to(ROOT)).replace('\\', '/'), sha256=sha(terms.read_bytes()),
                      bytes=terms.stat().st_size, url=TOU) if terms.exists() else None
    source_doc = dict(
        source=SOURCE, url=FTP, terms=TOU, license=LICENSE, license_uri=LICENSE_URI,
        version=version, access_date=access_date, maps=len(map_records), modules=len(kept),
        families=len(family_rows), admitted_to_open_core=True, partition='open_core',
        import_policy='CellDesigner CD4.1 via celldesigner.decode + identity_labels + diagram_modules',
        attribution=cite['attribution'], layout_origin=dict(LAYOUT_ORIGIN_EXPERT),
        files=file_hashes, licence_evidence=terms_meta, decode_failures=failures,
        excluded_modules=len(excluded), identity_layer='layout_atlas/identity_labels.py',
    )
    dump(EXT / 'SOURCE.json', source_doc)
    dump(EXT / 'import_audit.json', dict(per_map=per_map, failures=failures, excluded=excluded[:200]))
    summary = dict(
        source=SOURCE, version=version, access_date=access_date, license=LICENSE,
        maps=len(map_records), modules=len(kept), families=len(family_rows),
        variants=sum(len(f['variants']) for f in family_rows),
        sizes=dict(sorted(Counter(m['size'] for m in kept).items())),
        events_total=sum(x['events'] for x in per_map), glyphs_total=sum(x['glyphs'] for x in per_map),
        decode_failures=len(failures), excluded_modules=len(excluded),
        code_sha256=sha(Path(__file__).read_bytes()),
        celldesigner_sha256=sha((BASE / 'celldesigner.py').read_bytes()),
        identity_labels_sha256=sha((BASE / 'identity_labels.py').read_bytes()),
    )
    dump(RESULTS / 'library_summary.json', summary)
    return kept, summary, source_doc


def register_in_source_registry(source_doc):
    reg_path = INPUTS / 'source_registry.json'
    reg = read(reg_path) if reg_path.exists() else dict(schema=1, files=[], license_evidence=[])
    files = [f for f in reg.get('files', []) if f.get('source') != SOURCE]
    files.extend(source_doc.get('files') or [])
    reg['files'] = files
    evidence = [e for e in reg.get('license_evidence', []) if e.get('source') != SOURCE]
    if source_doc.get('licence_evidence'):
        le = source_doc['licence_evidence']
        evidence.append(dict(
            source=SOURCE, kind='license_evidence', access_date=source_doc.get('access_date'),
            url=le.get('url'), sha256=le.get('sha256'), bytes=le.get('bytes'),
            local_path=le.get('path'), license=LICENSE,
        ))
    reg['license_evidence'] = evidence
    reg['access_date'] = source_doc.get('access_date') or reg.get('access_date')
    reg['panther_pathway_version'] = source_doc.get('version')
    dump(reg_path, reg)
    return reg


if __name__ == '__main__':
    if os.environ.get('PYTHONHASHSEED') != '0':
        os.environ['PYTHONHASHSEED'] = '0'
        os.execv(sys.executable, [sys.executable, str(Path(__file__).resolve())] + sys.argv[1:])
    kept, summary, src = build()
    register_in_source_registry(src)
    print(json.dumps(summary, indent=2))
