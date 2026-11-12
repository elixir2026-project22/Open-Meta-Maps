"""Atlas schema version, provenance normalization, registry manifest, and validation.

Generic for any GEM-aligned Atlas source. Does not change layout algorithms.
"""
from __future__ import annotations

import gzip
import json
import sqlite3
import zlib
from collections import Counter, defaultdict
from datetime import date, datetime, timezone
from pathlib import Path

from core import DATA, INPUTS, dump, read, sha

SCHEMA_VERSION = 1
OPEN_CORE = DATA / 'open_core.sqlite'
ATLAS = DATA / 'atlas.sqlite'
REGISTRY_PATH = DATA / 'atlas_registry.json'
MANIFEST_PATH = DATA / 'input_manifest.json'
EXTRACTION_VERSION = 'atlas-extraction-1'

LAYOUT_ORIGIN_EXPERT = {
    'status': 'expert_hand_drawn',
    'counted_as_expert_design': True,
    'reason': 'source_map_expert_layout',
}

# Prefer this for new curated sources (user-facing wording: curated, not "expert").
LAYOUT_ORIGIN_CURATED = {
    'status': 'curated_hand_drawn',
    'counted_as_expert_design': True,
    'reason': 'source_map_curated_layout',
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
    'ich360': (
        'iCH360 Visualisation Escher maps (MIT); '
        'https://github.com/marco-corrao/iCH360; '
        'Corrao et al. doi:10.1371/journal.pcbi.1013564'
    ),
}

DEFAULT_EXTENSIONS = (
    'brenda',
    'human_maps',
    'wikipathways_geometry',
    'metabolismregulation',
    'panther',
    'ich360',
)

# Default cutter parameters used when backfilling extraction onto existing records.
# Stage B: size-scaled roots + persist only independent (non-nested) windows.
_DIAGRAM_PARAMS = dict(
    max_roots='scaled',
    window_sizes=[1, 2, 4, 8, 12, 24],
    independent_only=True,
    nest_factor=2.0,
    root_policy='all_roots_upto_48_else_max_roots_24',
)
_EXTRACTION_BY_SOURCE = {
    'escher': (
        'escher_community_windows',
        dict(
            max_roots='all',
            window_sizes=[1, 2, 4, 8, 12, 24, 40],
            link='primary_segment_entity',
            plus_full_component=True,
            max_reactions=40,
            max_nodes=100,
        ),
        EXTRACTION_VERSION,
        False,
    ),
    'reactome': ('diagram_modules', dict(_DIAGRAM_PARAMS), EXTRACTION_VERSION, False),
    'wikipathways': ('diagram_modules', dict(_DIAGRAM_PARAMS), EXTRACTION_VERSION, False),
    'human_maps': ('diagram_modules', dict(_DIAGRAM_PARAMS), EXTRACTION_VERSION, False),
    'metabolismregulation': ('diagram_modules', dict(_DIAGRAM_PARAMS), EXTRACTION_VERSION, False),
    'panther': ('diagram_modules', dict(_DIAGRAM_PARAMS), EXTRACTION_VERSION, False),
    'ich360': ('diagram_modules', dict(_DIAGRAM_PARAMS), EXTRACTION_VERSION, False),
    'brenda': (
        'brenda_library.extract',
        dict(
            window_sizes=[4, 8, 12, 24],
            pair_degree_cap=8,
            max_module_nodes=100,
            include_singletons=True,
            include_components=True,
        ),
        EXTRACTION_VERSION,
        False,
    ),
}

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


def _extraction_spec(source, quality=None):
    """Return (method, parameters, version, inferred) for a source/quality pair."""
    if source == 'wikipathways' and quality == 'geometry_only':
        return (
            'geometry_modules',
            dict(_DIAGRAM_PARAMS, underlying='diagram_modules over metabolite–metabolite GPML edges'),
            EXTRACTION_VERSION,
            False,
        )
    if source in _EXTRACTION_BY_SOURCE:
        return _EXTRACTION_BY_SOURCE[source]
    return ('unknown', {}, EXTRACTION_VERSION, True)


def _module_reaction_keys(m):
    keys = set()
    for r in m.get('reactions') or []:
        keys.add(r.get('key') or r.get('id'))
    if not keys:
        for n in m.get('nodes') or []:
            if n.get('role') == 'reaction':
                keys.add(n['key'])
    return keys


def _counts_for_modules(mods, source_events=None):
    covered = set()
    for m in mods:
        covered |= _module_reaction_keys(m)
    n_cov = len(covered)
    n_events = source_events if source_events is not None else (n_cov or None)
    coverage = round(n_cov / n_events, 4) if n_events else None
    return dict(
        source_events=n_events,
        modules=len(mods),
        reactions_covered=n_cov,
        coverage=coverage,
        independent_modules=len(mods),
    )


def _extraction_for_group(mods):
    """Build one extraction block for all modules of the same map cut.

    Preserve an existing complete block when every module already agrees
    (e.g. escher_community_windows from escher_cutter).
    """
    existing = [(m.get('extraction') or {}) for m in mods]
    methods = {e.get('method') for e in existing if e.get('method')}
    if len(methods) == 1 and all(e.get('method') and e.get('parameters') is not None for e in existing):
        # Keep cutter-authored method/params/version; refresh aggregate counts.
        base = dict(existing[0])
        base['counts'] = _counts_for_modules(mods, source_events=(base.get('counts') or {}).get('source_events'))
        return base
    src = mods[0].get('source')
    quality = mods[0].get('quality')
    method, params, version, inferred = _extraction_spec(src, quality)
    block = dict(
        method=method,
        version=version,
        parameters=dict(params),
        counts=_counts_for_modules(mods),
    )
    if inferred:
        block['inferred'] = True
    return block


def backfill_extension_extraction(names=None):
    """Rewrite extension modules.json.gz (and SOURCE.json extraction) with provenance blocks."""
    names = list(names or DEFAULT_EXTENSIONS) + ['units']
    stats = {}
    for name in names:
        ext = DATA / 'extensions' / name
        path = ext / 'modules.json.gz'
        if not path.exists():
            continue
        with gzip.open(path, 'rt') as f:
            items = json.load(f)
        groups = defaultdict(list)
        for m in items:
            groups[(m.get('source'), m.get('map_id') or m.get('id'))].append(m)
        stamped = 0
        for mods in groups.values():
            block = _extraction_for_group(mods)
            for m in mods:
                m['extraction'] = block
                stamped += 1
        with gzip.open(path, 'wt') as f:
            json.dump(items, f, ensure_ascii=False)
        src_path = ext / 'SOURCE.json'
        if src_path.exists():
            src = read(src_path)
            methods = Counter((m.get('extraction') or {}).get('method') for m in items)
            src['extraction'] = dict(
                methods=dict(methods),
                modules=len(items),
                atlas_schema_version=SCHEMA_VERSION,
            )
            dump(src_path, src)
        stats[name] = stamped
    return stats


def per_source_extraction_summary(modules_by_source):
    """Build registry extraction.by_source from module payloads."""
    out = {}
    for source, mods in sorted(modules_by_source.items()):
        methods = Counter((m.get('extraction') or {}).get('method') or '?' for m in mods)
        maps = {m.get('map_id') for m in mods if m.get('map_id')}
        covered = set()
        for m in mods:
            covered |= _module_reaction_keys(m)
        top_method = methods.most_common(1)[0][0] if methods else '?'
        sample = next(
            (m.get('extraction') for m in mods if (m.get('extraction') or {}).get('method') == top_method),
            {},
        ) or {}
        out[source] = dict(
            method=top_method,
            methods=dict(methods),
            parameters=sample.get('parameters') or {},
            version=sample.get('version'),
            inferred=bool(sample.get('inferred')),
            maps=len(maps),
            modules=len(mods),
            reactions_covered=len(covered),
        )
    return out


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

    mod_rows = []
    for mid, source, species, lineage, quality, size, core_family, variant, blob in db.execute(
        'SELECT id, source, species, lineage, quality, size, core_family, variant, payload FROM modules'
    ):
        mod_rows.append((mid, source, species, lineage, quality, size, core_family, variant, unpack(blob)))
    by_map = defaultdict(list)
    for mid, source, species, lineage, quality, size, core_family, variant, p in mod_rows:
        by_map[(p.get('source') or source, p.get('map_id') or mid)].append(p)
    extraction_by_map = {key: _extraction_for_group(mods) for key, mods in by_map.items()}

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
        key = (src, mid)
        if key in extraction_by_map:
            p['extraction'] = extraction_by_map[key]
            map_stats['map_extraction_set'] += 1
        elif not (p.get('extraction') or {}).get('method'):
            method, params, version, inferred = _extraction_spec(src)
            block = dict(method=method, version=version, parameters=dict(params), counts=dict(modules=0))
            if inferred:
                block['inferred'] = True
            p['extraction'] = block
            map_stats['map_extraction_default'] += 1
        db.execute(
            'UPDATE maps SET source=?, species=?, lineage=?, payload=? WHERE id=?',
            (p['source'], p.get('species') or species, p.get('lineage') or lineage, pack(p), mid),
        )
        map_stats['maps_updated'] += 1

    mod_stats = Counter()
    dangling = 0
    for mid, source, species, lineage, quality, size, core_family, variant, p in mod_rows:
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
        key = (src, map_id or mid)
        p['extraction'] = extraction_by_map[key]
        mod_stats['module_extraction_set'] += 1
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

    ext_stats = backfill_extension_extraction()

    manifest = refresh_input_manifest()
    registry = write_registry(db, manifest=manifest, path=path)
    db.close()

    if path.resolve() == OPEN_CORE.resolve() and ATLAS.exists():
        ATLAS.write_bytes(path.read_bytes())

    return dict(
        schema_version=SCHEMA_VERSION,
        map_stats=dict(map_stats),
        mod_stats=dict(mod_stats),
        extension_extraction=ext_stats,
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
    modules_by_source = defaultdict(list)
    for blob, in db.execute('SELECT payload FROM modules'):
        p = unpack(blob)
        licences[p.get('license') or p.get('licence') or '?'] += 1
        modules_by_source[p.get('source') or '?'].append(p)
    for name in DEFAULT_EXTENSIONS:
        ext = DATA / 'extensions' / name / 'modules.json.gz'
        if ext.exists():
            with gzip.open(ext, 'rt') as f:
                items = json.load(f)
            for m in items:
                licences[m.get('license') or m.get('licence') or '?'] += 1
                modules_by_source[m.get('source') or name].append(m)
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
        extraction=dict(
            atlas_schema_version=SCHEMA_VERSION,
            pipeline='database.build / escher_cutter / diagram_modules / from_legacy',
            cutters={
                'escher_community_windows': 'escher_cutter (all roots; sizes 1–40)',
                'diagram_modules': 'database.diagram_modules BFS windows',
                'geometry_modules': 'diagram_modules over metabolite–metabolite GPML edges',
                'brenda_library.extract': 'BRENDA SVG cutter (singletons, pairs, windows, components)',
                'units_library.decompose': 'non-overlapping structure units (companion layer)',
                'legacy_observations': 'from_legacy (geometry-only observations; unreproducible)',
            },
            by_source=per_source_extraction_summary(modules_by_source),
        ),
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
            if not (p.get('extraction') or {}).get('method'):
                errors.append(f'module {mid}: missing extraction.method')
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
