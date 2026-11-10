"""Audit of the stored layout variants: how many would merge if every alignment (not only the first 8) were tried?

database.build() puts a module into an existing variant when the best geometry error over the FIRST 8 isomorphisms of
its core graph is < 0.25.  Graphs with many automorphisms (symmetric trees, leaves that can be swapped) have far more
isomorphisms, so two nearly identical drawings can be split.  This read-only audit recomputes the best error over up to
`CAP` isomorphisms for pairs of variant representatives inside each family and merges variants whose best error is
below the same 0.25 threshold (union-find).  Nothing in the stored library is changed.
  python3 variant_audit.py [max_variants_per_family]
"""
import itertools
import sys
import time
from core import *
import database

CAP = 3000
THRESHOLD = .25


def best_error(a, b):
    ga, gb = database.graph(a), database.graph(b)
    matcher = nx.isomorphism.GraphMatcher(ga, gb, node_match=lambda x, y: x['role'] == y['role'])
    best, n = None, 0
    for mp in itertools.islice(matcher.isomorphisms_iter(), CAP):
        n += 1
        e = database.geometry_error(a, b, mp)
        if best is None or e < best:
            best = e
            if best < THRESHOLD:
                break
    return best, n


def run(max_variants=12):
    lib = read(OUT / 'library.json')
    con = database.connect()
    t0 = time.time()
    rows, families_checked, before, after, pairs, capped = [], 0, 0, 0, 0, 0
    for fam in lib['families']:
        variants = fam['variants']
        if not 2 <= len(variants) <= max_variants:
            continue
        reps = {}
        for v in variants:
            row = con.execute('SELECT payload FROM modules WHERE id=?', (v['template'],)).fetchone()
            reps[v['id']] = database.unpack(row[0])
        parent = {v['id']: v['id'] for v in variants}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        ids = [v['id'] for v in variants]
        for i, j in itertools.combinations(range(len(ids)), 2):
            if find(ids[i]) == find(ids[j]):
                continue
            pairs += 1
            err, n = best_error(reps[ids[i]], reps[ids[j]])
            if n >= CAP and (err is None or err >= THRESHOLD):
                capped += 1
            if err is not None and err < THRESHOLD:
                parent[find(ids[i])] = find(ids[j])
        roots = len({find(x) for x in ids})
        families_checked += 1
        before += len(ids)
        after += roots
        if roots < len(ids):
            rows.append(dict(family=fam['id'], variants_before=len(ids), variants_after=roots))
    summary = dict(families_checked=families_checked, max_variants_per_family=max_variants, variants_before=before, variants_after=after,
                   merged_variants=before - after, share_merged=(before - after) / max(1, before), families_with_merges=len(rows),
                   pair_checks=pairs, pair_checks_hitting_cap_without_result=capped, iso_cap=CAP, threshold=THRESHOLD,
                   seconds=round(time.time() - t0, 1), examples=rows[:15],
                   note='Variants are represented by their template module; merging uses the same 0.25 threshold as database.build() but tries up to %d alignments instead of 8.' % CAP)
    dump(OUT / 'variant_audit.json', summary)
    return summary


if __name__ == '__main__':
    print(json.dumps(run(int(sys.argv[1]) if len(sys.argv) > 1 else 12), indent=1)[:2500])
