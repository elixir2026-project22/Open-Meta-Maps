"""Make a GEM that no identity channel of the atlas can recognise (test input for the from-scratch path).

Takes a real model (default examples/e_coli_core.xml) and keeps ONLY its network structure: reaction and metabolite ids are
replaced by random tokens without a BiGG-style compartment suffix, names and every annotation (bigg/rhea/ec/chebi/kegg/...)
are removed, genes are dropped.  Stoichiometry, bounds and compartment labels are kept, so the model is chemically a valid
network.  Because the atlas also offers topology-only ('geometry') candidates, the test that uses this file asserts that the
identity channels (complete equation, reaction reference, EC) return nothing; a run with the library switched off checks the
fully from-scratch path.
  python3 make_nomatch_gem.py [source.xml] [out.json]
"""
import random
import sys
from core import *


def scramble(model, seed=20261109):
    rng = random.Random(seed)
    token = lambda prefix, used: next(t for t in (prefix + '%08x' % rng.getrandbits(32) for _ in iter(int, 1)) if t not in used)
    mids, rids = {}, {}
    for m in model['metabolites']:
        mids[m['id']] = token('X', set(mids.values()))
    for r in model['reactions']:
        rids[r['id']] = token('R', set(rids.values()))
    return dict(id='scrambled_network', name='',
                metabolites=[dict(id=mids[m['id']], name='', compartment=m.get('compartment', 'c')) for m in model['metabolites']],
                reactions=[dict(id=rids[r['id']], name='', metabolites={mids[k]: v for k, v in r['metabolites'].items()},
                                lower_bound=r['lower_bound'], upper_bound=r['upper_bound'], gene_reaction_rule='') for r in model['reactions']],
                genes=[])


if __name__ == '__main__':
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'examples/e_coli_core.xml'
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / 'examples/scrambled_network.json'
    dump(out, scramble(load_model(src)))
    print('wrote', out)
