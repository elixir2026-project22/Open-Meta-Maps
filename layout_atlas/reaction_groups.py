"""Derived index: the same set of reactions drawn in different source maps ("reaction groups").

The atlas has two different notions of "same":

* family  -- same CONNECTIVITY SHAPE (graph isomorphism of the reaction/metabolite core), reactions may differ;
* reaction group (this module) -- the SAME REACTIONS (equal sets of canonical reaction identities), drawn in one or
  several source maps / GEMs.  Modules in one group are alternative expert drawings of the same chemistry.

Canonical reaction identity, in priority order (the origin is stored per reaction):
  1. Rhea master id given by the source reaction (`rhea_direct`);
  2. a BiGG reaction id that the fixed public Ref models annotate with exactly one Rhea master id (`rhea_via_bigg_annotation`);
  3. the BiGG reaction id itself (`bigg_id`);
  4. the Reactome reaction id (`reactome_id`).
Identities from different namespaces are never equated, so groups are a LOWER bound (a BiGG-only module and a Rhea-only
module of the same chemistry stay in different groups).  Geometry-only modules carry no reaction identity and are not
grouped.  Nothing in the stored library is modified; output goes to results/reaction_groups.json.gz.
"""
import gzip
from core import *
import database
from semantics import rhea_index


def bigg_to_rhea():
    table = defaultdict(set)
    master = rhea_index()['master']
    for name in ('e_coli_core', 'iML1515', 'iMM904', 'Human-GEM'):
        model = load_model(REF / 'models' / (name + '.xml'))
        for r in model['reactions']:
            refs = xrefs(r.get('annotations', []))
            ids = set(refs.get('bigg.reaction', [])) | ({r['id']} if name != 'Human-GEM' else set())
            rhea = {master.get(str(v), str(v)) for v in refs.get('rhea', [])}
            for b in ids:
                table[b] |= rhea
    return table


def canonical(reaction, master, b2r):
    refs = reaction.get('references', {})
    rhea = sorted({master.get(str(v), str(v)) for v in refs.get('rhea', [])})
    if rhea:
        return 'rhea:' + rhea[0], 'rhea_direct'
    for v in sorted(map(str, refs.get('bigg.reaction', []))):
        found = b2r.get(v)
        if found and len(found) == 1:
            return 'rhea:' + next(iter(found)), 'rhea_via_bigg_annotation'
        return 'bigg:' + v, 'bigg_id'
    for v in sorted(map(str, refs.get('reactome', []))):
        return 'reactome:' + v, 'reactome_id'
    return None, None


def build():
    master, b2r = rhea_index()['master'], bigg_to_rhea()
    mods = database.modules()
    by_set, module_index, origins, skipped = defaultdict(list), {}, Counter(), Counter()
    for m in mods:
        if not m.get('reactions'):
            skipped['no_reaction_identity_(geometry_only)'] += 1
            continue
        ids = []
        for r in m['reactions']:
            cid, how = canonical(r, master, b2r)
            if cid is None:
                ids = None
                break
            ids.append(cid)
            origins[how] += 1
        if ids is None or len(set(ids)) != len(ids):
            skipped['missing_or_duplicate_identity'] += 1
            continue
        by_set[tuple(sorted(ids))].append(m)
    groups = []
    for key in sorted(by_set, key=lambda k: (len(k), k)):
        ms = by_set[key]
        maps, sources, lineages = {x['map_id'] for x in ms}, {x['source'] for x in ms}, {x['lineage'] for x in ms}
        group = dict(id=uid(['reaction-group', key]), reactions=list(key), size=len(key), modules=sorted(x['id'] for x in ms),
                     source_maps=len(maps), sources=sorted(sources), lineages=len(lineages),
                     families=len({x['core_family'] for x in ms}), variants=len({x['variant'] for x in ms}),
                     cross_map=len(maps) >= 2, cross_source=len(sources) >= 2, cross_lineage=len(lineages) >= 2)
        groups.append(group)
        for x in ms:
            module_index[x['id']] = group['id']
    multi = [g for g in groups if g['cross_map']]
    summary = dict(
        modules_total=len(mods), modules_grouped=len(module_index), groups=len(groups), skipped=dict(skipped),
        identity_origin_counts=dict(origins),
        groups_in_more_than_one_source_map=len(multi), modules_in_those_groups=sum(len(g['modules']) for g in multi),
        groups_across_sources=sum(g['cross_source'] for g in groups), groups_across_lineages=sum(g['cross_lineage'] for g in groups),
        multi_map_groups_by_size=dict(sorted(Counter(g['size'] for g in multi).items())),
        multi_map_groups_by_size_tier=dict(Counter('micro(1)' if g['size'] == 1 else 'small(2-3)' if g['size'] <= 3 else 'structure(>=4)' for g in multi)),
        groups_with_several_families=sum(g['families'] > 1 for g in groups),
        interpretation='Lower bound: identities from different namespaces are not equated; geometry-only modules have no reaction identity.',
        code_sha256=sha(Path(__file__).read_bytes()))
    OUT.mkdir(parents=True, exist_ok=True)
    with gzip.open(OUT / 'reaction_groups.json.gz', 'wt') as f:
        json.dump(dict(summary=summary, groups=groups, module_group=module_index), f, ensure_ascii=False)
    dump(OUT / 'reaction_groups_summary.json', summary)
    return summary


if __name__ == '__main__':
    print(json.dumps(build(), indent=1))
