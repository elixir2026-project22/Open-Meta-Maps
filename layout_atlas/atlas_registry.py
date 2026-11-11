"""Atlas schema version, provenance normalization, registry manifest, and validation.

Generic for any GEM-aligned Atlas source. Does not change layout algorithms.
"""
from __future__ import annotations

import gzip
import json
import sqlite3
import zlib
from collections import Counter
from datetime import date, datetime, timezone
from pathlib import Path

from core import DATA, INPUTS, dump, read, sha

SCHEMA_VERSION = 1
OPEN_CORE = DATA / 'open_core.sqlite'
ATLAS = DATA / 'atlas.sqlite'
REGISTRY_PATH = DATA / 'atlas_registry.json'
MANIFEST_PATH = DATA / 'input_manifest.json'

LAYOUT_ORIGIN_EXPERT = {
    'status': 'expert_hand_drawn',
    'counted_as_expert_design': True,
    'reason': 'source_map_expert_layout',
}

ATTRIBUTION = {
    'escher': 'Escher community-maps (CC0); https://github.com/escher/community-maps',
    'wikipathways': 'WikiPathways (CC0); https://www.wikipathways.org',
    'reactome': 'Reactome pathway data (CC0); https://reactome.org',
    'brenda': 'BRENDA; Chang et al. NAR 2026 doi:10.1093/nar/gkaf1113',
    'human_maps': 'Human-maps (CC BY 4.0); https://github.com/SysBioChalmers/Human-maps',
    'metabolismregulation': (
        'Metabolism Regulation gallery (CC BY 4.0); '
        'https://metabolismregulation.github.io/gallery/'
    ),
    'panther': (
        'PANTHER Pathway CellDesigner maps (CC BY 4.0); '
        'https://data.pantherdb.org/ftp/pathway/3.4/CD4.1/'
    ),
}

DEFAULT_EXTENSIONS = (
    'brenda',
    'human_maps',
    'wikipathways_geometry',
    'metabolismregulation',
    'panther',
)

DDL = '''
CREATE TABLE IF NOT EXISTS modules(
  id TEXT PRIMARY KEY, source TEXT, species TEXT, lineage TEXT,
  quality TEXT, size INTEGER, core_family TEXT, variant TEXT, payload BLOB);
CREATE INDEX IF NOT EXISTS core_family_idx ON modules(core_family);
CREATE INDEX IF NOT EXISTS modules_source_idx ON modules(source);
CREATE INDEX IF NOT EXISTS modules_lineage_idx ON modules(lineage);
CREATE TABLE IF NOT EXISTS reactions(
  module_id TEXT, reaction_key TEXT, namespace TEXT, identifier TEXT,
  PRIMARY KEY(module_id, reaction_key, namespace, identifier));
CREATE INDEX IF NOT EXISTS biological_idx ON reactions(namespace, identifier);
CREATE TABLE IF NOT EXISTS maps(
  id TEXT PRIMARY KEY, source TEXT, species TEXT, lineage TEXT, payload BLOB);
CREATE INDEX IF NOT EXISTS maps_source_idx ON maps(source);
CREATE INDEX IF NOT EXISTS maps_lineage_idx ON maps(lineage);
CREATE TABLE IF NOT EXISTS reviews(
  module_id TEXT, reviewer TEXT, status TEXT, payload BLOB,
  PRIMARY KEY(module_id, reviewer));
CREATE TABLE IF NOT EXISTS retained_observations(
  id TEXT PRIMARY KEY, reason TEXT, source_pointer TEXT);
CREATE TABLE IF NOT EXISTS meta(
  key TEXT PRIMARY KEY, value TEXT NOT NULL);
'''


def pack(v):
    return zlib.compress(json.dumps(v, ensure_ascii=False).encode(), 6)


def unpack(v):
    return json.loads(zlib.decompress(v))


def connect_path(path=None):
    path = Path(path or ATLAS)
    path.parent.mkdir(parents=True, exist_ok=True)
    db = sqlite3.connect(path)
    ensure_schema(db)
    return db


def ensure_schema(db):
    db.executescript(DDL)
    row = db.execute("SELECT value FROM meta WHERE key='schema_version'").fetchone()
    if row is None:
        db.execute(
            "INSERT INTO meta(key, value) VALUES('schema_version', ?)",
            (str(SCHEMA_VERSION),),
        )
        db.execute(
            "INSERT OR REPLACE INTO meta(key, value) VALUES('schema_updated', ?)",
            (datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),),
        )
        db.commit()
    return int(db.execute("SELECT value FROM meta WHERE key='schema_version'").fetchone()[0])


def refresh_input_manifest():
    """Rewrite data/input_manifest.json from current source_inputs files."""
    files = {
        str(p.relative_to(INPUTS)): sha(p.read_bytes())
        for p in sorted(INPUTS.rglob('*'))
        if p.is_file()
    }
    result = dict(version='0.3-source-inputs', generated=str(date.today()), files=files)
    dump(MANIFEST_PATH, result)
    return result


def _registry_by_sha():
    reg = read(INPUTS / 'source_registry.json')
    return {f['sha256']: f for f in reg.get('files', []) if f.get('sha256')}, reg


def _reactome_urls():
    acq = read(DATA / 'acquisition.json')
    out = {}
    for rec in acq.get('reactome_diagrams', []):
        rid = rec.get('id')
        if not rid:
            continue
        out[rid] = dict(
            url=rec.get('url') or rec.get('layout_url') or f'https://reactome.org/PathwayBrowser/#{rid}',
            version=str(acq.get('reactome_release') or '97'),
            source_id=rid,
        )
    return out


def _fix_escher_registry_paths(reg):
    """Point Escher registry paths at data/source_inputs/raw/escher/... when present."""
    changed = 0
    for f in reg.get('files', []):
        if f.get('source') != 'escher':
            continue
        sid = f.get('id') or ''
        candidate = INPUTS / 'raw' / 'escher' / sid
        if candidate.is_file():
            new_path = str(candidate.relative_to(DATA.parent)).replace('\\', '/')
            # store repo-relative path
            new_path = str(Path('data/source_inputs/raw/escher') / sid)
            if f.get('path') != new_path:
                f['path'] = new_path
                changed += 1
    return changed


def normalize_open_core(path=None):
    """Fill consistent provenance fields on maps/modules; set schema meta; write registry."""
    path = Path(path or OPEN_CORE)
    db = sqlite3.connect(path)
    ensure_schema(db)
    by_sha, reg = _registry_by_sha()
    path_fixes = _fix_escher_registry_paths(reg)
    if path_fixes:
        dump(INPUTS / 'source_registry.json', reg)
        by_sha, reg = _registry_by_sha()
    reactome = _reactome_urls()

    map_stats = Counter()
    map_ids = {}
    for mid, source, species, lineage, blob in db.execute(
        'SELECT id, source, species, lineage, payload FROM maps'
    ):
        p = unpack(blob)
        map_ids[mid] = p
        src = p.get('source') or source
        p['source'] = src
        if not p.get('license') and not p.get('licence'):
            map_stats['map_missing_license'] += 1
        lic = p.get('license') or p.get('licence')
        if lic:
            p['license'] = lic
        p.setdefault('species', species)
        p.setdefault('lineage', lineage)
        p.setdefault('attribution', ATTRIBUTION.get(src))
        p.setdefault('layout_origin', dict(LAYOUT_ORIGIN_EXPERT))
        meta = by_sha.get(p.get('source_sha256') or '')
        if meta:
            p.setdefault('url', meta.get('url'))
            p.setdefault('version', meta.get('version'))
            p.setdefault('source_id', p.get('source_id') or meta.get('id'))
            if meta.get('path') and 'modular_layout' not in str(meta.get('path')):
                p.setdefault('source_path', meta.get('path'))
        if src == 'wikipathways' and str(p.get('lineage', '')).startswith('WP'):
            p['provenance_status'] = 'resolved_pathway_id'
        elif src == 'escher':
            p.setdefault('provenance_status', 'manual_history_unverified')
        elif src == 'reactome':
            rid = (p.get('lineage') or '').split(':', 1)[-1]
            info = reactome.get(rid) or reactome.get(p.get('map_id') or '')
            if info:
                p.setdefault('url', info.get('url'))
                p.setdefault('version', info.get('version'))
                p.setdefault('source_id', info.get('source_id'))
            p.setdefault('provenance_status', 'reactome_release_member')
        p['atlas_schema_version'] = SCHEMA_VERSION
        db.execute(
            'UPDATE maps SET source=?, species=?, lineage=?, payload=? WHERE id=?',
            (p['source'], p.get('species') or species, p.get('lineage') or lineage, pack(p), mid),
        )
        map_stats['maps_updated'] += 1

    mod_stats = Counter()
    dangling = 0
    for mid, source, species, lineage, quality, size, core_family, variant, blob in db.execute(
        'SELECT id, source, species, lineage, quality, size, core_family, variant, payload FROM modules'
    ):
        p = unpack(blob)
        src = p.get('source') or source
        p['source'] = src
        lic = p.get('license') or p.get('licence')
        if not lic:
            mod_stats['module_missing_license'] += 1
        else:
            p['license'] = lic
        p.setdefault('species', species)
        p.setdefault('lineage', lineage)
        p.setdefault('attribution', ATTRIBUTION.get(src))
        if 'layout_origin' not in p or not p['layout_origin']:
            p['layout_origin'] = dict(LAYOUT_ORIGIN_EXPERT)
            mod_stats['layout_origin_set'] += 1
        p['atlas_schema_version'] = SCHEMA_VERSION
        map_id = p.get('map_id')
        if map_id and map_id not in map_ids:
            dangling += 1
        db.execute(
            'UPDATE modules SET source=?, species=?, lineage=?, payload=? WHERE id=?',
            (p['source'], p.get('species') or species, p.get('lineage') or lineage, pack(p), mid),
        )
        mod_stats['modules_updated'] += 1

    db.execute(
        "INSERT OR REPLACE INTO meta(key, value) VALUES('schema_version', ?)",
        (str(SCHEMA_VERSION),),
    )
    db.execute(
        "INSERT OR REPLACE INTO meta(key, value) VALUES('normalized', ?)",
        (datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),),
    )
    db.execute(
        "INSERT OR REPLACE INTO meta(key, value) VALUES('registry', ?)",
        ('data/atlas_registry.json',),
    )
    db.commit()

    manifest = refresh_input_manifest()
    registry = write_registry(db, manifest=manifest, path=path)
    db.close()

    # Keep runtime working copy in sync when present
    if path.resolve() == OPEN_CORE.resolve() and ATLAS.exists():
        ATLAS.write_bytes(path.read_bytes())

    return dict(
        schema_version=SCHEMA_VERSION,
        map_stats=dict(map_stats),
        mod_stats=dict(mod_stats),
        dangling_map_id=dangling,
        escher_registry_paths_fixed=path_fixes,
        registry=str(REGISTRY_PATH),
        manifest_files=len(manifest['files']),
    )


def write_registry(db=None, manifest=None, path=None):
    close = False
    if db is None:
        db = sqlite3.connect(path or OPEN_CORE)
        ensure_schema(db)
        close = True
    if manifest is None and MANIFEST_PATH.exists():
        manifest = read(MANIFEST_PATH)
    maps = dict(db.execute('SELECT source, count(*) FROM maps GROUP BY source'))
    mods = dict(db.execute('SELECT source, count(*) FROM modules GROUP BY source'))
    nrx = db.execute('SELECT count(*) FROM reactions').fetchone()[0]
    licences = Counter()
    for blob, in db.execute('SELECT payload FROM modules'):
        p = unpack(blob)
        licences[p.get('license') or p.get('licence') or '?'] += 1
    for name in DEFAULT_EXTENSIONS:
        ext = DATA / 'extensions' / name / 'modules.json.gz'
        if ext.exists():
            with gzip.open(ext, 'rt') as f:
                items = json.load(f)
            for m in items:
                licences[m.get('license') or m.get('licence') or '?'] += 1
            mods[name] = mods.get(name, 0) + len(items)
    ext_meta = {}
    for name in DEFAULT_EXTENSIONS:
        src = DATA / 'extensions' / name / 'SOURCE.json'
        if src.exists():
            ext_meta[name] = read(src)
    core_path = Path(path or OPEN_CORE)
    registry = dict(
        schema_version=SCHEMA_VERSION,
        generated=str(date.today()),
        open_core=dict(
            path='data/open_core.sqlite',
            sha256=sha(core_path.read_bytes()) if core_path.exists() else None,
            bytes=core_path.stat().st_size if core_path.exists() else None,
            maps=maps,
            modules={k: v for k, v in mods.items() if k in maps or k in ('escher', 'reactome', 'wikipathways')},
            modules_with_extensions=mods,
            reactions=nrx,
            licences=dict(licences),
        ),
        extensions=ext_meta,
        input_manifest=dict(
            version=(manifest or {}).get('version'),
            files=len((manifest or {}).get('files') or {}),
            sha256=sha(MANIFEST_PATH.read_bytes()) if MANIFEST_PATH.exists() else None,
        ),
        source_registry=dict(
            path='data/source_inputs/source_registry.json',
            sha256=sha((INPUTS / 'source_registry.json').read_bytes())
            if (INPUTS / 'source_registry.json').exists()
            else None,
            wikipathways_release=read(INPUTS / 'source_registry.json').get('wikipathways_release')
            if (INPUTS / 'source_registry.json').exists()
            else None,
        ),
        extraction=dict(pipeline='database.build / diagram_modules / from_legacy', atlas_schema_version=SCHEMA_VERSION),
    )
    dump(REGISTRY_PATH, registry)
    if close:
        db.close()
    return registry


def validate(path=None, include_extensions=True):
    """Return a list of error strings; empty means OK."""
    path = Path(path or OPEN_CORE)
    errors = []
    if not path.exists():
        return [f'missing database: {path}']
    db = sqlite3.connect(path)
    try:
        ensure_schema(db)
        ver = db.execute("SELECT value FROM meta WHERE key='schema_version'").fetchone()
        if not ver or int(ver[0]) < 1:
            errors.append('meta.schema_version missing or < 1')
        map_ids = {r[0] for r in db.execute('SELECT id FROM maps')}
        for mid, source, blob in db.execute('SELECT id, source, payload FROM maps'):
            p = unpack(blob)
            if not (p.get('license') or p.get('licence')):
                errors.append(f'map {mid}: missing license')
            if not p.get('source'):
                errors.append(f'map {mid}: missing source')
            if not p.get('source_sha256'):
                errors.append(f'map {mid}: missing source_sha256')
            if not p.get('layout_origin'):
                errors.append(f'map {mid}: missing layout_origin')
            if not p.get('provenance_status'):
                errors.append(f'map {mid}: missing provenance_status')
        for mid, source, blob in db.execute('SELECT id, source, payload FROM modules'):
            p = unpack(blob)
            if not (p.get('license') or p.get('licence')):
                errors.append(f'module {mid}: missing license')
            if not p.get('source'):
                errors.append(f'module {mid}: missing source')
            if not p.get('source_sha256'):
                errors.append(f'module {mid}: missing source_sha256')
            if not p.get('layout_origin'):
                errors.append(f'module {mid}: missing layout_origin')
            map_id = p.get('map_id')
            if map_id and map_id not in map_ids:
                errors.append(f'module {mid}: dangling map_id {map_id}')
        # indexes present
        idx = {r[0] for r in db.execute("SELECT name FROM sqlite_master WHERE type='index'")}
        for name in ('modules_source_idx', 'modules_lineage_idx', 'maps_source_idx', 'maps_lineage_idx'):
            if name not in idx:
                errors.append(f'missing index {name}')
    finally:
        db.close()

    if include_extensions:
        for name in DEFAULT_EXTENSIONS:
            ext = DATA / 'extensions' / name / 'modules.json.gz'
            if not ext.exists():
                continue
            with gzip.open(ext, 'rt') as f:
                items = json.load(f)
            for i, m in enumerate(items):
                if not (m.get('license') or m.get('licence')):
                    errors.append(f'extension {name}[{i}]: missing license')
                if not m.get('layout_origin'):
                    errors.append(f'extension {name}[{i}]: missing layout_origin')
                if not m.get('source'):
                    errors.append(f'extension {name}[{i}]: missing source')
                if name == 'metabolismregulation':
                    if not m.get('source_sha256'):
                        errors.append(f'extension {name}[{i}]: missing source_sha256')
                    cite = m.get('source_citation') or {}
                    if not (cite.get('attribution') or m.get('attribution')):
                        errors.append(f'extension {name}[{i}]: missing attribution')

    if not REGISTRY_PATH.exists():
        errors.append('missing data/atlas_registry.json')
    return errors


if __name__ == '__main__':
    import pprint
    pprint.pp(normalize_open_core())
    errs = validate()
    print('validation_errors', len(errs))
    for e in errs[:20]:
        print(' ', e)
    raise SystemExit(1 if errs else 0)
