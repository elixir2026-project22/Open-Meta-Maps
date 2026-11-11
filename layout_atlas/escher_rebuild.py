"""Replace open-core Escher modules with the in-repo community cutter.

Preserves WikiPathways and Reactome rows (same pattern as wikipathways_rebuild).
Writes ``data/open_core.sqlite``, refreshes ``atlas.sqlite`` copy, registry, and a
small parity/gate summary under ``layout_atlas/results/``.
"""
from __future__ import annotations

import shutil
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

import networkx as nx
import numpy as np

from atlas_registry import (
    LAYOUT_ORIGIN_EXPERT,
    OPEN_CORE,
    SCHEMA_VERSION,
    ensure_schema,
    normalize_open_core,
    pack,
    unpack,
)
from core import DATA, dump, normalize, uid
from database import geometry_error, graph, mappings
from escher_cutter import cut_all, default_parameters, parity_report

ATLAS = DATA / 'atlas.sqlite'
OUT = Path(__file__).resolve().parent / 'results'


def _load_non_escher(db):
    maps, mods = [], []
    for row in db.execute('SELECT id, source, species, lineage, payload FROM maps'):
        if row[1] == 'escher':
            continue
        maps.append(unpack(row[-1]))
    for row in db.execute(
        'SELECT id, source, species, lineage, quality, size, core_family, variant, payload FROM modules'
    ):
        if row[1] == 'escher':
            continue
        mods.append(unpack(row[-1]))
    return maps, mods


def _classify_and_write(db, maps, mods):
    db.execute('DELETE FROM reactions')
    db.execute('DELETE FROM modules')
    db.execute('DELETE FROM maps')
    families = defaultdict(list)
    familyrows = []
    excluded = []
    for i, m in enumerate(sorted(mods, key=lambda m: m['id'])):
        if i % 2000 == 0:
            print(f'Index {i}/{len(mods)}', flush=True)
        xy = normalize({n['key']: n['xy'] for n in m['nodes']})
        for n in m['nodes']:
            n['xy'] = xy[n['key']]
        m['core_coordinate_normalization'] = (
            'centered RMS radius; original full instance retained by immutable source pointer'
        )
        gm = graph(m)
        if not len(gm) or not nx.is_connected(gm):
            excluded.append(dict(
                id=m['id'], source=m['source'], map_id=m.get('map_id'),
                reason='empty_or_disconnected_visible_core',
                original_instance=m.get('original_instance'),
            ))
            continue
        h = nx.weisfeiler_lehman_graph_hash(gm, node_attr='role')
        bucket = (len(gm), gm.number_of_edges(), h)
        family = None
        for f in families[bucket]:
            if mappings(gm, f['graph'], 1):
                family = f
                break
        if family is None:
            family = dict(
                id=uid(['core-v03', bucket, m['id']]), graph=gm,
                representative=m, variants=[], instances=[],
            )
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
        m['geometry_fingerprint'] = uid(
            np.round(np.sort(np.linalg.norm(pts[:, None] - pts[None, :], axis=2).ravel()), 3).tolist()
        )
        m.setdefault('layout_origin', dict(LAYOUT_ORIGIN_EXPERT))
        m.setdefault('atlas_schema_version', SCHEMA_VERSION)
        db.execute(
            'INSERT OR REPLACE INTO modules VALUES(?,?,?,?,?,?,?,?,?)',
            (m['id'], m['source'], m['species'], m['lineage'], m['quality'], m['size'],
             m['core_family'], m['variant'], pack(m)),
        )
        for r in m.get('reactions') or []:
            for ns, values in (r.get('references') or {}).items():
                for val in values:
                    db.execute(
                        'INSERT OR IGNORE INTO reactions VALUES(?,?,?,?)',
                        (m['id'], r['key'], ns, str(val)),
                    )
    seen = set()
    for d in maps:
        if d['id'] in seen:
            continue
        seen.add(d['id'])
        d.setdefault('atlas_schema_version', SCHEMA_VERSION)
        d.setdefault('layout_origin', dict(LAYOUT_ORIGIN_EXPERT))
        db.execute(
            'INSERT OR REPLACE INTO maps VALUES(?,?,?,?,?)',
            (d['id'], d['source'], d['species'], d['lineage'], pack(d)),
        )
    db.commit()
    return dict(
        excluded=len(excluded),
        families=len(familyrows),
        variants=sum(len(f['variants']) for f in familyrows),
        excluded_rows=excluded,
    )


def rebuild(open_core=None, *, parameters=None, sync_atlas=True):
    open_core = Path(open_core or OPEN_CORE)
    parameters = dict(parameters or default_parameters())
    print('Cutting community Escher maps...', flush=True)
    cut = cut_all(parameters=parameters)
    print(f"  maps={len(cut['maps'])} modules={len(cut['modules'])}", flush=True)
    print('Parity vs stored observations...', flush=True)
    parity = parity_report(parameters=parameters)
    print(
        f"  recall={parity['summary']['recall']:.3f} "
        f"gen_cov_min={parity['summary']['min_generated_coverage']:.3f} "
        f"stored={parity['summary']['stored_modules']} "
        f"generated={parity['summary']['generated_modules']}",
        flush=True,
    )

    db = sqlite3.connect(open_core)
    ensure_schema(db)
    maps, mods = _load_non_escher(db)
    before = dict(Counter(m['source'] for m in mods))
    before['escher'] = db.execute("SELECT COUNT(*) FROM modules WHERE source='escher'").fetchone()[0]
    maps.extend(cut['maps'])
    mods.extend(cut['modules'])
    print(f'Classifying {len(mods)} modules (non-Escher preserved)...', flush=True)
    stats = _classify_and_write(db, maps, mods)
    after = {s: n for s, n in db.execute('SELECT source, COUNT(*) FROM modules GROUP BY source')}
    db.execute(
        "INSERT OR REPLACE INTO meta(key, value) VALUES('schema_version', ?)",
        (str(SCHEMA_VERSION),),
    )
    db.commit()
    db.close()

    norm = normalize_open_core(open_core)
    if sync_atlas:
        shutil.copyfile(open_core, ATLAS)
        # clear cached open-core modules in this process
        try:
            from database import _core_modules_cached
            _core_modules_cached.cache_clear()
        except Exception:
            pass

    summary = dict(
        parameters=parameters,
        before_modules_by_source=before,
        after_modules_by_source=after,
        cut_per_map=cut['per_map'],
        parity=parity['summary'],
        classify=dict(excluded=stats['excluded'], families=stats['families'], variants=stats['variants']),
        normalize=norm,
        note=(
            'In-repo escher_community_windows replaces from_legacy for the 22 community maps. '
            'Stored module_observations.json.gz remains the parity reference and still feeds '
            'WikiPathways geometry_only rows via from_legacy.'
        ),
    )
    OUT.mkdir(parents=True, exist_ok=True)
    dump(OUT / 'escher_cutter_rebuild.json', summary)
    dump(OUT / 'escher_cutter_parity.json', parity)
    print(summary['after_modules_by_source'], flush=True)
    return summary


if __name__ == '__main__':
    rebuild()
