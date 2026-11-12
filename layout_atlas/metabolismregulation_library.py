"""Metabolism Regulation gallery maps (CC BY 4.0) via the generic SBGN-ML importer.

Source: https://metabolismregulation.github.io/gallery/
Licence: Creative Commons Attribution 4.0 International
  (https://metabolismregulation.github.io/about/).
Layouts are expert / hand-drawn SBGN-ML (``.sbgn``) files published under
``downloads/`` (gallery) and additional atlas-downloads stems when present.

Imported into the default Atlas as extension ``metabolismregulation`` (loaded by
``database.modules()``), with full provenance in SOURCE.json and source_registry.
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
import sbgn

EXT = DATA / 'extensions/metabolismregulation'
RESULTS = OUT / 'metabolismregulation_extension'
RAW = INPUTS / 'raw/metabolismregulation'
TERMS = DATA / 'terms/metabolismregulation_about.html'

GALLERY_REPO = 'metabolismregulation/metabolismregulation.github.io'
GALLERY_PAGE = 'https://metabolismregulation.github.io/gallery/'
ABOUT_PAGE = 'https://metabolismregulation.github.io/about/'
LICENSE = 'CC-BY-4.0'
LICENSE_URI = 'https://creativecommons.org/licenses/by/4.0/'
SOURCE = 'metabolismregulation'

# Human-readable titles for gallery stems (fallback: title-case the slug).
TITLES = {
    'F001-glycolysis': 'Glycolysis',
    'F001-glycolysis-cardiac': 'Glycolysis (cardiac)',
    'F002-eicosanoids': 'Eicosanoid production',
    'F003-mevalonate': 'Mevalonate biosynthesis',
    'F004-glycogen-liver': 'Glycogen metabolism (liver)',
    'F004-glycogen-muscle': 'Glycogen metabolism (muscle)',
    'F005-srebp': 'SREBP signalling',
    'F006-acly': 'Acetyl-CoA synthesis',
    'F007-inos': 'iNOS pathway',
    'F008-enos': 'eNOS pathway',
    'F009-phosphoinositides': 'Phosphoinositide',
    'F012-camp': 'cAMP hydrolysis',
    'F014-carnitine': 'Acyl-carnitine synthesis',
    'F015-tag': 'Triacylglycerol hydrolysis',
    'F018-estrogen': 'Estradiol biosynthesis',
    'F019-glutaminase': 'Glutaminolysis',
    'F020-pyrimidine': 'Pyrimidine biosynthesis',
    'F021-ureacycle': 'Urea cycle',
    'F022-bcaa': 'Branched-chain amino acids',
    'F023-pdh': 'Pyruvate dehydrogenase',
    'F025-nnos': 'nNOS pathway',
    'F100-aspirin': 'Aspirin / eicosanoid (drug view)',
}


def _title(stem):
    if stem in TITLES:
        return TITLES[stem]
    slug = stem.split('-', 1)[-1].replace('-', ' ')
    return slug[:1].upper() + slug[1:]


def citation(version, access_date=None):
    return dict(
        title='Metabolism Regulation gallery maps',
        url=GALLERY_PAGE,
        repository=f'https://github.com/{GALLERY_REPO}',
        license=dict(name='CC BY 4.0', uri=LICENSE_URI),
        licence_page=ABOUT_PAGE,
        version=version,
        retrieved=access_date or str(date.today()),
        attribution=(
            'Map diagrams and layouts: Metabolism Regulation project '
            f'({GALLERY_PAGE}), Creative Commons Attribution 4.0 International.'
        ),
        changes=(
            'SBGN-ML layouts decoded with layout_atlas/sbgn.py; module windows via '
            'database.diagram_modules; coordinates normalized for storage; originals kept.'
        ),
    )


def list_local_sbgn():
    return sorted(RAW.glob('*.sbgn'))


def _repo_paths():
    meta = RAW / 'FILES.json'
    if not meta.exists():
        return {}
    return {row['name']: row for row in read(meta)}


def decode_all(version, access_date=None):
    cite = citation(version, access_date)
    paths = _repo_paths()
    diagrams, failures = [], []
    for path in list_local_sbgn():
        stem = path.stem
        info = paths.get(path.name) or {}
        repo_path = info.get('repo_path') or f'downloads/{path.name}'
        blob_url = f'https://github.com/{GALLERY_REPO}/blob/{version}/{repo_path}'
        try:
            d = sbgn.decode(
                path,
                source=SOURCE,
                license=LICENSE,
                name=_title(stem),
                species='Homo sapiens',
                lineage=f'{SOURCE}:{stem}',
                partition='open_core',
                admitted_to_open_core=True,
                reaction_namespace='sbgn.process',
                source_id=stem,
                source_citation=cite,
                attribution=ATTRIBUTION.get(SOURCE) or cite['attribution'],
                layout_origin=dict(LAYOUT_ORIGIN_EXPERT),
                provenance_status='gallery_sbgn_ml',
                url=blob_url,
                source_path=repo_path,
                version=version,
                gallery_page=GALLERY_PAGE,
            )
            diagrams.append(d)
        except Exception as e:
            failures.append(dict(path=str(path), error=str(e)))
    return diagrams, failures, cite


def classify(mods):
    """Assign core_family / variant / geometry_fingerprint (same scheme as BRENDA / Human-maps)."""
    families, family_rows = defaultdict(list), []
    kept = []
    excluded = []
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
            family = dict(id=uid(['mreg-core', bucket, m['id']]), graph=gm, variants=[], instances=[])
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
        if i % 500 == 0:
            print('metabolismregulation index', i, '/', len(mods), flush=True)
    return kept, family_rows, excluded


def build(version=None, access_date=None):
    EXT.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    if version is None:
        # Prefer a pinned commit recorded next to the raw files when present
        pin = RAW / 'SOURCE_COMMIT.txt'
        version = pin.read_text().strip() if pin.exists() else 'HEAD'
    access_date = access_date or str(date.today())
    assert list_local_sbgn(), f'no .sbgn files under {RAW}'

    diagrams, failures, cite = decode_all(version, access_date)
    modules = []
    map_records = []
    per_map = []
    for d in diagrams:
        mods = database.diagram_modules(d)
        for m in mods:
            m['partition'] = 'open_core'
            m['admitted_to_open_core'] = True
            m['layout_origin'] = dict(LAYOUT_ORIGIN_EXPERT)
            m['source_citation'] = cite
            m['attribution'] = cite['attribution']
            m['provenance_status'] = 'gallery_sbgn_ml'
            m['version'] = version
            m['url'] = d.get('url')
            m['source_id'] = d.get('source_id')
            # diagram_modules sets map_id to diagram id; keep lineage from the map stem
            m['lineage'] = d['lineage']
            m['map_id'] = d['id']
        modules.extend(mods)
        map_records.append(dict(
            id=d['id'],
            source=SOURCE,
            source_id=d.get('source_id'),
            species=d['species'],
            lineage=d['lineage'],
            name=d['name'],
            license=LICENSE,
            source_sha256=d['source_sha256'],
            provenance_status='gallery_sbgn_ml',
            attribution=cite['attribution'],
            layout_origin=dict(LAYOUT_ORIGIN_EXPERT),
            version=version,
            url=d.get('url'),
            path=d.get('path'),
            format='sbgn-ml',
            language=d.get('language'),
            events=len(d['events']),
            glyphs=len(d['glyphs']),
            modules=len(mods),
            audit=d.get('audit'),
            extraction=(mods[0].get('extraction') if mods else database.extraction_block(
                'diagram_modules',
                dict(max_roots='scaled', window_sizes=[1, 2, 4, 8, 12, 24], independent_only=True,
                     root_policy='scaled_max_24_ceil_N_4_or_all_upto_48'),
                dict(source_events=len(d['events']), modules=0, reactions_covered=0, coverage=None,
                     independent_modules=0),
            )),
        ))
        per_map.append(dict(
            source_id=d.get('source_id'),
            name=d['name'],
            sha256=d['source_sha256'],
            events=len(d['events']),
            glyphs=len(d['glyphs']),
            modules=len(mods),
            excluded_processes=len(d.get('excluded') or []),
        ))

    kept, family_rows, excluded = classify(modules)

    # Extension sqlite + modules.json.gz
    db = sqlite3.connect(EXT / 'atlas.sqlite')
    db.executescript(
        '''CREATE TABLE IF NOT EXISTS sources(id TEXT PRIMARY KEY, payload BLOB);
        CREATE TABLE IF NOT EXISTS maps(id TEXT PRIMARY KEY, payload BLOB);
        CREATE TABLE IF NOT EXISTS modules(id TEXT PRIMARY KEY, size INTEGER, core_family TEXT, variant TEXT, payload BLOB);
        CREATE TABLE IF NOT EXISTS reactions(module_id TEXT, reaction_key TEXT, namespace TEXT, identifier TEXT,
          PRIMARY KEY(module_id, reaction_key, namespace, identifier));
        DELETE FROM sources; DELETE FROM maps; DELETE FROM modules; DELETE FROM reactions;'''
    )
    db.execute('INSERT INTO sources VALUES(?,?)', (version[:24], database.pack(cite)))
    for rec in map_records:
        db.execute('INSERT INTO maps VALUES(?,?)', (rec['id'], database.pack(rec)))
    for m in kept:
        db.execute(
            'INSERT INTO modules VALUES(?,?,?,?,?)',
            (m['id'], m['size'], m['core_family'], m['variant'], database.pack(m)),
        )
        for r in m['reactions']:
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
    dump(
        EXT / 'families.json',
        dict(
            families=[
                dict(
                    id=f['id'],
                    instances=f['instances'],
                    variants=[
                        dict(id=v['id'], template=v['representative']['id'], instances=v['instances'])
                        for v in f['variants']
                    ],
                )
                for f in family_rows
            ]
        ),
    )

    paths = _repo_paths()
    file_hashes = []
    for p in list_local_sbgn():
        raw = p.read_bytes()
        info = paths.get(p.name) or {}
        repo_path = info.get('repo_path') or f'downloads/{p.name}'
        file_hashes.append(dict(
            id=p.stem,
            source=SOURCE,
            path=str(p.relative_to(ROOT)).replace('\\', '/'),
            repo_path=repo_path,
            sha256=sha(raw),
            bytes=len(raw),
            license=LICENSE,
            license_uri=LICENSE_URI,
            version=version,
            url=f'https://raw.githubusercontent.com/{GALLERY_REPO}/{version}/{repo_path}',
            format='sbgn-ml',
            access_date=access_date,
        ))

    terms_meta = None
    if TERMS.exists():
        terms_meta = dict(
            path=str(TERMS.relative_to(ROOT)).replace('\\', '/'),
            sha256=sha(TERMS.read_bytes()),
            bytes=TERMS.stat().st_size,
            url=ABOUT_PAGE,
        )

    source_doc = dict(
        source=SOURCE,
        url=GALLERY_PAGE,
        repository=f'https://github.com/{GALLERY_REPO}',
        license=LICENSE,
        license_uri=LICENSE_URI,
        licence_page=ABOUT_PAGE,
        version=version,
        access_date=access_date,
        maps=len(map_records),
        modules=len(kept),
        families=len(family_rows),
        admitted_to_open_core=True,
        partition='open_core',
        import_policy='gallery + atlas-downloads SBGN-ML via layout_atlas/sbgn.py + database.diagram_modules',
        attribution=cite['attribution'],
        layout_origin=dict(LAYOUT_ORIGIN_EXPERT),
        files=file_hashes,
        licence_evidence=terms_meta,
        decode_failures=failures,
        excluded_modules=len(excluded),
    )
    dump(EXT / 'SOURCE.json', source_doc)
    dump(EXT / 'import_audit.json', dict(per_map=per_map, failures=failures, excluded=excluded[:200]))

    summary = dict(
        source=SOURCE,
        version=version,
        access_date=access_date,
        license=LICENSE,
        maps=len(map_records),
        modules=len(kept),
        families=len(family_rows),
        variants=sum(len(f['variants']) for f in family_rows),
        sizes=dict(sorted(Counter(m['size'] for m in kept).items())),
        events_total=sum(x['events'] for x in per_map),
        decode_failures=len(failures),
        excluded_modules=len(excluded),
        code_sha256=sha(Path(__file__).read_bytes()),
        sbgn_importer_sha256=sha((BASE / 'sbgn.py').read_bytes()),
    )
    dump(RESULTS / 'library_summary.json', summary)
    return kept, summary, source_doc


def register_in_source_registry(source_doc):
    """Append / replace metabolismregulation file entries in source_registry.json."""
    reg_path = INPUTS / 'source_registry.json'
    reg = read(reg_path) if reg_path.exists() else dict(schema=1, files=[], license_evidence=[])
    prior = {f.get('path'): f for f in reg.get('files', []) if f.get('source') == SOURCE}
    files = [f for f in reg.get('files', []) if f.get('source') != SOURCE]
    for row in source_doc.get('files') or []:
        old = prior.get(row.get('path')) or {}
        row = dict(row)
        row.setdefault('extraction_method', old.get('extraction_method') or 'sbgn.decode + diagram_modules')
        row.setdefault('retrieval_date', old.get('retrieval_date') or source_doc.get('access_date'))
        for k, v in old.items():
            if k not in row and v is not None:
                row[k] = v
        files.append(row)
    reg['files'] = files
    evidence = [e for e in reg.get('license_evidence', []) if e.get('source') != SOURCE]
    if source_doc.get('licence_evidence'):
        le = source_doc['licence_evidence']
        evidence.append(dict(
            source=SOURCE,
            kind='license_evidence',
            access_date=source_doc.get('access_date'),
            url=le.get('url'),
            sha256=le.get('sha256'),
            bytes=le.get('bytes'),
            local_path=le.get('path'),
            license=LICENSE,
        ))
    reg['license_evidence'] = evidence
    reg['access_date'] = source_doc.get('access_date') or reg.get('access_date')
    reg.setdefault('metabolismregulation_release', source_doc.get('version'))
    dump(reg_path, reg)
    return reg


if __name__ == '__main__':
    if os.environ.get('PYTHONHASHSEED') != '0':
        os.environ['PYTHONHASHSEED'] = '0'
        os.execv(sys.executable, [sys.executable, str(Path(__file__).resolve())] + sys.argv[1:])
    ver = (RAW / 'SOURCE_COMMIT.txt').read_text().strip() if (RAW / 'SOURCE_COMMIT.txt').exists() else None
    kept, summary, src = build(version=ver)
    register_in_source_registry(src)
    print(json.dumps(summary, indent=2))
