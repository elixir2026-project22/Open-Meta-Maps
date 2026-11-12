"""iCH360 custom Escher maps (MIT) via database.escher_diagram + diagram_modules.

Source: https://github.com/marco-corrao/iCH360/tree/main/Visualisation
Licence: MIT (repo-root LICENSE; Visualisation/ has no overriding terms).
Corroboration: Zenodo paper deposit of the same Visualisation tree is CC BY 4.0
(https://doi.org/10.5281/zenodo.15048160). Paper (Corrao et al., PLOS Comp Biol)
describes them as custom-built Escher maps; Data Availability points at GitHub.

Imported as default Atlas extension ``ich360``. compressed_map.json is excluded
(CP pseudo-reactions). Pathway maps keep independent curated geometry even when
their reaction sets are subsets of full_map (Atlas already keeps overlapping
drawings from one source; retrieval ranks candidates).
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
from atlas_registry import LAYOUT_ORIGIN_CURATED, ATTRIBUTION

EXT = DATA / 'extensions/ich360'
RESULTS = OUT / 'ich360_extension'
RAW = INPUTS / 'raw/ich360'
TERMS_LICENSE = DATA / 'terms/ich360_LICENSE.txt'
TERMS_ZENODO = DATA / 'terms/ich360_zenodo_15048160.html'

REPO = 'marco-corrao/iCH360'
REPO_URL = f'https://github.com/{REPO}'
VIS_URL = f'{REPO_URL}/tree/main/Visualisation'
LICENSE = 'MIT'
LICENSE_URI = 'https://opensource.org/licenses/MIT'
SOURCE = 'ich360'
SPECIES = 'Escherichia coli'

TITLES = {
    'full_map': 'iCH360 full metabolic map',
    'aminoacids_biosynthesis': 'Amino-acid biosynthesis',
    'c1_pool': 'C1 metabolism',
    'fatty_acids_biosynthesis': 'Fatty-acid biosynthesis',
    'nucleotides_biosynthesis': 'Nucleotide biosynthesis',
    'cofactors_de_novo_biosynthesis': 'Cofactor de novo biosynthesis',
    'ions_metals_uptake': 'Ion and metal uptake',
    'kdo_biosynthesis': 'KDO-lipid biosynthesis',
    'murein_biosynthesis': 'Murein biosynthesis',
    'phospholipids_biosynthesis': 'Phospholipid biosynthesis',
}


def citation(version, access_date=None):
    return dict(
        title='iCH360 custom Escher maps',
        url=VIS_URL,
        repository=REPO_URL,
        license=dict(name='MIT', uri=LICENSE_URI),
        version=version,
        retrieved=access_date or str(date.today()),
        attribution=(
            'iCH360 Visualisation Escher maps (https://github.com/marco-corrao/iCH360), '
            'MIT License; Corrao et al., PLOS Computational Biology '
            '(doi:10.1371/journal.pcbi.1013564).'
        ),
        changes=(
            'Official Escher JSON decoded with database.escher_diagram; module windows via '
            'database.diagram_modules; compressed_map.json excluded (CP pseudo-reactions).'
        ),
        paper='Corrao et al. 2025, A compact model of Escherichia coli core and biosynthetic metabolism, '
              'PLOS Computational Biology doi:10.1371/journal.pcbi.1013564',
        zenodo_corroboration='https://doi.org/10.5281/zenodo.15048160 (CC BY 4.0 deposit of the same Visualisation tree)',
    )


def _commit():
    pin = RAW / 'SOURCE_COMMIT.txt'
    return pin.read_text().strip() if pin.exists() else 'HEAD'


def list_local():
    return sorted(p for p in RAW.rglob('*.json') if p.name != 'FILES.json')


def decode_all(version, access_date=None):
    cite = citation(version, access_date)
    diagrams, failures = [], []
    for path in list_local():
        stem = path.stem
        rel = str(path.relative_to(RAW)).replace('\\', '/')
        blob_url = f'https://github.com/{REPO}/blob/{version}/Visualisation/{rel}'
        try:
            d = database.escher_diagram(
                path,
                source=SOURCE,
                license=LICENSE,
                species=SPECIES,
                reaction_namespace='bigg.reaction',
                name=TITLES.get(stem, stem.replace('_', ' ')),
                lineage=f'{SOURCE}:{stem}',
                partition='open_core',
                admitted_to_open_core=True,
                source_id=stem,
                source_citation=cite,
                attribution=ATTRIBUTION.get(SOURCE) or cite['attribution'],
                layout_origin=dict(LAYOUT_ORIGIN_CURATED),
                provenance_status='ich360_custom_escher',
                url=blob_url,
                source_path=f'Visualisation/{rel}',
                version=version,
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
            family = dict(id=uid(['ich360-core', bucket, m['id']]), graph=gm, variants=[], instances=[])
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
            print('ich360 index', i, '/', len(mods), flush=True)
    return kept, family_rows, excluded


def build(version=None, access_date=None):
    EXT.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    version = version or _commit()
    access_date = access_date or str(date.today())
    assert list_local(), f'no Escher JSON under {RAW}; run tools/fetch_ich360_maps.py'

    diagrams, failures, cite = decode_all(version, access_date)
    modules, map_records, per_map = [], [], []
    for d in diagrams:
        mods = database.diagram_modules(d)
        for m in mods:
            m['partition'] = 'open_core'
            m['admitted_to_open_core'] = True
            m['layout_origin'] = dict(LAYOUT_ORIGIN_CURATED)
            m['source_citation'] = cite
            m['attribution'] = cite['attribution']
            m['provenance_status'] = 'ich360_custom_escher'
            m['version'] = version
            m['url'] = d.get('url')
            m['source_id'] = d.get('source_id')
            m['lineage'] = d['lineage']
            m['map_id'] = d['id']
            m['counted_as_expert_design'] = True
        modules.extend(mods)
        map_records.append(dict(
            id=d['id'], source=SOURCE, source_id=d.get('source_id'), species=d['species'],
            lineage=d['lineage'], name=d['name'], license=LICENSE, source_sha256=d['source_sha256'],
            provenance_status='ich360_custom_escher', attribution=cite['attribution'],
            layout_origin=dict(LAYOUT_ORIGIN_CURATED), version=version, url=d.get('url'),
            path=d.get('path'), format='escher', events=len(d['events']), glyphs=len(d['glyphs']),
            modules=len(mods), audit=d.get('audit'),
        ))
        per_map.append(dict(
            source_id=d.get('source_id'), name=d['name'], sha256=d['source_sha256'],
            events=len(d['events']), glyphs=len(d['glyphs']), modules=len(mods),
            skipped=d.get('audit'),
        ))

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
    db.execute('INSERT INTO sources VALUES(?,?)', (version[:24], database.pack(cite)))
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
        rel = str(p.relative_to(RAW)).replace('\\', '/')
        file_hashes.append(dict(
            id=p.stem, source=SOURCE,
            path=str(p.relative_to(ROOT)).replace('\\', '/'),
            repo_path=f'Visualisation/{rel}',
            sha256=sha(raw), bytes=len(raw), license=LICENSE, license_uri=LICENSE_URI,
            version=version,
            url=f'https://raw.githubusercontent.com/{REPO}/{version}/Visualisation/{rel}',
            format='escher', access_date=access_date, use='core',
            family='ecoli', species=SPECIES,
        ))

    licence_evidence = []
    if TERMS_LICENSE.exists():
        licence_evidence.append(dict(
            path=str(TERMS_LICENSE.relative_to(ROOT)).replace('\\', '/'),
            sha256=sha(TERMS_LICENSE.read_bytes()), bytes=TERMS_LICENSE.stat().st_size,
            url=f'https://raw.githubusercontent.com/{REPO}/{version}/LICENSE',
            kind='repo_LICENSE_MIT',
        ))
    if TERMS_ZENODO.exists():
        licence_evidence.append(dict(
            path=str(TERMS_ZENODO.relative_to(ROOT)).replace('\\', '/'),
            sha256=sha(TERMS_ZENODO.read_bytes()), bytes=TERMS_ZENODO.stat().st_size,
            url='https://doi.org/10.5281/zenodo.15048160',
            kind='zenodo_CC-BY-4.0_corroboration',
        ))

    source_doc = dict(
        source=SOURCE, url=VIS_URL, repository=REPO_URL, license=LICENSE, license_uri=LICENSE_URI,
        version=version, access_date=access_date, maps=len(map_records), modules=len(kept),
        families=len(family_rows), admitted_to_open_core=True, partition='open_core',
        import_policy=(
            'custom Escher JSON via database.escher_diagram + diagram_modules; '
            'compressed_map.json excluded (CP pseudo-reactions)'
        ),
        dedupe_note=(
            'Include full_map + pathway_maps + pathways_not_in_model (10 maps). '
            'Exclude compressed_map (pseudo-reactions). Pathway maps are reaction-subsets of '
            'full_map but keep independent curated geometry; Atlas already retains overlapping '
            'drawings from one source (e.g. Human-maps compartment vs subsystem) and retrieval ranks them.'
        ),
        attribution=cite['attribution'], layout_origin=dict(LAYOUT_ORIGIN_CURATED),
        files=file_hashes, licence_evidence=licence_evidence,
        decode_failures=failures, excluded_modules=len(excluded),
        paper_doi='10.1371/journal.pcbi.1013564',
        zenodo_doi='10.5281/zenodo.15048160',
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
    )
    dump(RESULTS / 'library_summary.json', summary)
    return kept, summary, source_doc


def register_in_source_registry(source_doc):
    reg_path = INPUTS / 'source_registry.json'
    reg = read(reg_path) if reg_path.exists() else dict(schema=1, files=[], license_evidence=[])
    prior = {f.get('path'): f for f in reg.get('files', []) if f.get('source') == SOURCE}
    files = [f for f in reg.get('files', []) if f.get('source') != SOURCE]
    for row in source_doc.get('files') or []:
        old = prior.get(row.get('path')) or {}
        row = dict(row)
        row.setdefault('extraction_method', old.get('extraction_method') or 'escher_diagram + diagram_modules')
        row.setdefault('retrieval_date', old.get('retrieval_date') or source_doc.get('access_date'))
        for k, v in old.items():
            if k not in row and v is not None:
                row[k] = v
        files.append(row)
    reg['files'] = files
    evidence = [e for e in reg.get('license_evidence', []) if e.get('source') != SOURCE]
    for le in source_doc.get('licence_evidence') or []:
        evidence.append(dict(
            source=SOURCE, kind=le.get('kind') or 'license_evidence',
            access_date=source_doc.get('access_date'), url=le.get('url'),
            sha256=le.get('sha256'), bytes=le.get('bytes'), local_path=le.get('path'),
            license=LICENSE if 'MIT' in str(le.get('kind', '')) or 'LICENSE' in str(le.get('kind', '')) else 'CC-BY-4.0',
        ))
    reg['license_evidence'] = evidence
    reg['access_date'] = source_doc.get('access_date') or reg.get('access_date')
    reg['ich360_visualisation_commit'] = source_doc.get('version')
    dump(reg_path, reg)
    return reg


if __name__ == '__main__':
    if os.environ.get('PYTHONHASHSEED') != '0':
        os.environ['PYTHONHASHSEED'] = '0'
        os.execv(sys.executable, [sys.executable, str(Path(__file__).resolve())] + sys.argv[1:])
    ver = _commit()
    kept, summary, src = build(version=ver)
    register_in_source_registry(src)
    print(json.dumps(summary, indent=2))
