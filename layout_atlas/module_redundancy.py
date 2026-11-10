"""Derived index: how redundant are the stored module instances?

Modules are overlapping windows cut from source maps.  Within ONE source map, a module is *contained* when all its
core nodes also occur in another (larger) module of the same map; it then adds no drawing information that the larger
module does not already hold.  This script marks every module as `maximal` (not contained in another module of its
map) or `contained` (with the id of its smallest container) and reports the redundancy factor per source.
Read-only on the library; output: results/module_redundancy.json.gz and results/module_redundancy_summary.json.
"""
import gzip
from core import *
import database


def build():
    con = database.connect()
    by_map = defaultdict(list)
    for mid, source, size, payload in con.execute('SELECT id,source,size,payload FROM modules'):
        m = database.unpack(payload)
        keys = frozenset(n['key'] for n in m['nodes'] if not n['key'].startswith('route:'))
        react = frozenset(r['key'] for r in m.get('reactions', [])) or keys
        by_map[m['map_id']].append((mid, source, size, keys, react))
    flags, per_source = {}, defaultdict(Counter)
    distinct_cover = defaultdict(int)
    for map_id, items in by_map.items():
        items.sort(key=lambda x: (-len(x[3]), x[0]))
        for i, (mid, source, size, keys, react) in enumerate(items):
            container = None
            for other in items[:i]:
                if len(other[3]) > len(keys) and keys <= other[3]:
                    container = other[0]
                    break
            flags[mid] = dict(maximal=container is None, container=container)
            per_source[source]['modules'] += 1
            per_source[source]['maximal'] += container is None
            per_source[source]['sum_sizes'] += size
        distinct_cover[items[0][1]] += len(set().union(*[x[4] for x in items]))
    summary = dict(modules=len(flags), maximal=sum(f['maximal'] for f in flags.values()),
                   contained=sum(not f['maximal'] for f in flags.values()), source_maps=len(by_map),
                   by_source={s: dict(c, redundancy_factor=round(c['sum_sizes'] / max(1, distinct_cover[s]), 2), distinct_reactions_or_steps=distinct_cover[s]) for s, c in per_source.items()},
                   definition='contained = all core node keys also occur in a larger module of the same source map',
                   code_sha256=sha(Path(__file__).read_bytes()))
    summary['share_contained'] = summary['contained'] / summary['modules']
    OUT.mkdir(parents=True, exist_ok=True)
    with gzip.open(OUT / 'module_redundancy.json.gz', 'wt') as f:
        json.dump(flags, f)
    dump(OUT / 'module_redundancy_summary.json', summary)
    return summary


if __name__ == '__main__':
    print(json.dumps(build(), indent=1))
