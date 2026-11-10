"""How do EC candidates fare against the reaction equations? (retrieval-level, no drawing)  python3 ec_screening_audit.py [model ...]"""
import sys
from core import *
import database
import ec_channel as ec
import ec_provenance as epv
import tool_pipeline as tp

MODELS = {'e_coli_core': REF / 'models/e_coli_core.xml', 'iYO844': DATA / 'models/iYO844.xml', 'iMM904': REF / 'models/iMM904.xml',
          'iML1515': REF / 'models/iML1515.xml', 'AraCore_v2_1': DATA / 'models/AraCore_v2_1.xml'}


def run(names):
    bridge, core = tp.load_bridge(), database.modules()
    ec.enrich(core)
    templates = {m['id']: m for m in core}
    rows = []
    for name in names:
        model = load_model(MODELS[name])
        em = ec.ECMapper(model, bridge)
        proposals, info = em.proposals(core, 'open_core EC')
        # per-pair statuses of ALL candidate pairs (not only kept proposals)
        pairs = Counter()
        for p in proposals:
            src = {r['key']: r for r in templates[p['template_id']]['reactions']}
            for ev in p['mapping_evidence']:
                pairs[epv.classify_pair(src.get(ev['source_key']), em.reactions[ev['target_reaction']], em)[0]] += 1
        kept, stats = epv.screen(proposals, templates, em)
        rows.append(dict(model=name, reactions=len(model['reactions']), proposals=len(proposals), pair_checks=dict(pairs), screening=stats,
                         kept_reactions=len({r for p in kept for r in p['reaction_ids']})))
        print(rows[-1], flush=True)
    dump(OUT / 'ec_channel/screening_audit.json', rows)


if __name__ == '__main__':
    run(sys.argv[1:] or list(MODELS))
