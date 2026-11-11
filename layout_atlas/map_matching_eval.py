"""Evaluate model-reaction → Atlas / BRENDA node matching for development models.

  python3 layout_atlas/map_matching_eval.py [--no-vicinity] [--json OUT]

Reports match rates for e_coli_core and iIT341 (iMM904 / iML1515 skipped during development).
Open-data path only.
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import *
import brenda_library as bl
import ec_channel as ec
import identifier_bridge as ib


MODELS = ['e_coli_core', 'iIT341']


def atlas_ec_coverage(model):
    """Share of non-boundary reactions whose (current) EC appears on an open-core module reaction."""
    import database
    modules = database.modules()
    by_ec = set()
    for m in modules:
        for r in m.get('reactions', []):
            ecs, _ = ec.source_ecs(r)
            by_ec |= ecs
    target = ec.target_ecs(model)
    non_b = [r for r in model['reactions']
             if len(r['metabolites']) > 1 and not r['id'].startswith(('EX_', 'DM_', 'SK_', 'sink_'))]
    hit = sum(1 for r in non_b if target.get(r['id'], ()) & by_ec)
    return dict(non_boundary=len(non_b), with_ec_on_atlas=hit,
                share=round(hit / max(1, len(non_b)), 4))


def evaluate(vicinity='propagate'):
    rows = []
    for name in MODELS:
        model = load_model(EXAMPLES / f'{name}.xml')
        matches = bl.match_model(model, vicinity=vicinity)
        summary = bl.match_summary(matches)
        non_b = summary['reactions'] - summary['status'].get('boundary', 0)
        on = summary['status'].get('on map', 0)
        related = summary['status'].get('related', 0)
        atlas = atlas_ec_coverage(model)
        # hand-checkable accepted sample: top-scoring EC-linked matches
        sample = []
        for m in matches.values():
            for c in m['candidates']:
                if c['accepted'] and 'ec' in c['via']:
                    sample.append(dict(reaction=m['reaction_id'], ec=c.get('ec'), score=round(c['score'], 3),
                                       via=c['via'], enzyme=c.get('enzyme_name')))
        sample.sort(key=lambda x: -x['score'])
        rows.append(dict(
            model=name, reactions=summary['reactions'], non_boundary=non_b,
            brenda_on_map=on, brenda_share=round(on / max(1, non_b), 4),
            brenda_related=related, brenda_nodes=summary['nodes_with_model_reactions'],
            vicinity_only=summary['vicinity_only_reactions'],
            atlas_ec=atlas, sample=sample[:12], thresholds=ib.THRESHOLDS,
        ))
        print(f"{name}: BRENDA on map {on}/{non_b} ({100 * on / max(1, non_b):.1f}%), "
              f"related {related}, nodes {summary['nodes_with_model_reactions']}, "
              f"vicinity-only {summary['vicinity_only_reactions']}; "
              f"Atlas EC coverage {atlas['with_ec_on_atlas']}/{atlas['non_boundary']}")
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--no-vicinity', action='store_true')
    p.add_argument('--json')
    a = p.parse_args()
    rows = evaluate(vicinity=None if a.no_vicinity else 'propagate')
    if a.json:
        dump(a.json, dict(models=rows, method='open-data BRENDA/Atlas matching; ideas from ecmapper, independent rewrite'))


if __name__ == '__main__':
    main()
