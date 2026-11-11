"""Synthetic checks of the EC screening logic (python3 test_ec_provenance.py)."""
import copy
from core import *
import database
import ec_channel as ec
import ec_provenance as ep
import tool_pipeline as tp


def main():
    model = load_model(ROOT / 'examples/e_coli_core.xml')
    bridge = tp.load_bridge()
    em = ec.ECMapper(model, bridge)
    core = database.modules()
    # Open-core Escher windows are labelled `community` in the Atlas payload (registry still says escher).
    _ESCHER = {'escher', 'community'}
    source = next(r for m in core if m['source'] in _ESCHER for r in m['reactions'] if r['id'] == 'PGI' and r.get('coefficients'))
    target = {r['id']: r for r in model['reactions']}
    cases = []
    cases.append(('same reaction, equation maps', ep.classify_pair(source, em.reactions['PGI'], em)[0], 'verified_equivalent'))
    cases.append(('source without stoichiometry', ep.classify_pair(dict(source, coefficients={}), em.reactions['PGI'], em)[0], 'insufficient_evidence'))
    unknown = copy.deepcopy(source)
    first = next(iter(unknown['coefficients']))
    unknown['coefficients']['bigg:unknown_compound_zzz@c'] = unknown['coefficients'].pop(first)
    cases.append(('one participant of unknown identity', ep.classify_pair(unknown, em.reactions['PGI'], em)[0], 'insufficient_evidence'))
    cases.append(('identified participants, different reaction', ep.classify_pair(source, em.reactions['PFK'], em)[0], 'incompatible'))
    ok = True
    for name, got, want in cases:
        print(('PASS' if got == want else 'FAIL'), name, '->', got)
        ok &= got == want
    # provenance classification of a hand-made selection
    p_eq = dict(id='a', source='escher', match_level='reaction_reference', complete_chemical_equation_verified=True, reaction_ids=['PGI'], mapping_evidence=[dict(target_reaction='PGI', complete_equation=True)])
    p_geo = dict(id='b', source='wikipathways', match_level='geometry', reaction_ids=['PFK'], mapping_evidence=[])
    p_ec_ok = dict(id='c', source='escher', match_level='ec_topology', reaction_ids=['FBP'], ec_status={'FBP': dict(status='verified_equivalent', shared_ec=['3.1.3.11'], ambiguity=1)}, mapping_evidence=[])
    p_ec_weak = dict(id='d', source='brenda', match_level='ec_topology', reaction_ids=['PYK'], ec_status={'PYK': dict(status='insufficient_evidence', reason='no stoichiometry', shared_ec=['2.7.1.40'], ambiguity=2)}, mapping_evidence=[])
    prov = ep.reaction_provenance(model, [p_eq, p_geo, p_ec_ok, p_ec_weak], {'ENO'})
    expected = dict(PGI='expert_derived', PFK='topology_layout_suggestion', FBP='expert_derived', PYK='functional_layout_suggestion', ENO='bridge', PGK='generated')
    for rid, want in expected.items():
        got = prov[rid]['provenance']
        print(('PASS' if got == want else 'FAIL'), 'provenance', rid, '->', got)
        ok &= got == want
    # --- chemical compatibility of two equations (identity / compatibility / layout are separate judgements) ---------------------------------
    import copy as _copy
    mdh = next(r for m in core if m['source'] in _ESCHER for r in m['reactions'] if r['id'] == 'MDH' and r.get('coefficients'))
    def cmp(src, tid, mutate=None):
        tgt = _copy.deepcopy(em.reactions[tid])
        if mutate:
            mutate(tgt)
        return ep.compare_reactions(src, tgt, em)
    def swap(tgt, pairs):
        tgt['metabolites'] = {pairs.get(k, k): v for k, v in tgt['metabolites'].items()}
    cases2 = []
    r = cmp(source, 'PGI'); cases2.append(('same equation', (r['status'], 'exact_equation' in r['difference_types']), ('compatible', True)))
    r = cmp(dict(source, coefficients={k: -v for k, v in source['coefficients'].items()}), 'PGI'); cases2.append(('reversed recording', (r['status'], 'direction_difference' in r['difference_types']), ('compatible', True)))
    r = cmp(dict(source, coefficients={k: 2 * v for k, v in source['coefficients'].items()}), 'PGI'); cases2.append(('proportional scaling', (r['status'], 'proportional_scaling' in r['difference_types']), ('compatible', True)))
    r = cmp(mdh, 'MDH', lambda t: swap(t, {'nad_c': 'nadp_c', 'nadh_c': 'nadph_c'})); cases2.append(('NAD -> NADP', (r['status'], 'cofactor_substitution' in r['difference_types']), ('conflict', True)))
    r = cmp(mdh, 'MDH', lambda t: t['metabolites'].pop('h_c')); cases2.append(('proton missing', (r['status'], 'protonation_difference' in r['difference_types']), ('conflict', True)))
    r = cmp(source, 'PFK'); cases2.append(('different core substrates', (r['status'], 'core_substrate_difference' in r['difference_types']), ('conflict', True)))
    r = cmp(unknown, 'PGI'); cases2.append(('unknown participant is not a conflict', r['status'], 'unknown'))
    r = cmp(dict(source, coefficients={}), 'PGI'); cases2.append(('no stoichiometry is not a conflict', r['status'], 'unknown'))
    irr = dict(source, reversible=False)
    r = cmp(irr, 'PGI', lambda t: t.update(lower_bound=-1000.0, upper_bound=0.0)); cases2.append(('irreversible, opposite direction', (r['status'], 'irreversibility_conflict' in r['difference_types'], r['direction_compatible']), ('compatible', True, False)))
    for name, got, want in cases2:
        print(('PASS' if got == want else 'FAIL'), 'compat', name, '->', got)
        ok &= got == want
    # identifier match with a conflicting equation: identity revoked (not deleted), layout lent only with independent similarity
    pfk_ec = sorted(ec.target_ecs(model)['PFK'])
    templates = {'t1': dict(reactions=[dict(source, key='r:1')]), 't2': dict(reactions=[dict(source, key='r:2', references={'ec-code': pfk_ec})])}
    mk = lambda pid, tid, key: dict(id=pid, template_id=tid, map_id='m', source='escher', species='x', match_level='reaction_reference', reaction_ids=['PFK'], reaction_mapping={key: 'PFK'},
                                    mapping_evidence=[dict(target_reaction='PFK', reference=True, complete_equation=False)])
    prov = ep.reaction_provenance(model, [mk('c1', 't1', 'r:1')], set(), templates=templates, mapper=em)['PFK']
    got = (prov['provenance'], prov['identity_evidence']['revoked'], prov['identity_evidence']['kind'], prov['chemical_compatibility']['status'], prov['layout_provenance'])
    want = ('topology_layout_suggestion', True, 'identifier_supported', 'conflict', 'expert_topological')
    print(('PASS' if got == want else 'FAIL'), 'identifier match, conflicting equation, nothing else similar ->', got); ok &= got == want
    if pfk_ec:
        prov = ep.reaction_provenance(model, [mk('c2', 't2', 'r:2')], set(), templates=templates, mapper=em)['PFK']
        got = (prov['provenance'], prov['identity_evidence']['revoked'], prov['layout_provenance'], prov['transfer']['mode'])
        want = ('functional_layout_suggestion', True, 'expert_adapted', 'adapted')
        print(('PASS' if got == want else 'FAIL'), 'identifier match, conflicting equation, same full EC ->', got); ok &= got == want
    prov = ep.reaction_provenance(model, [mk('c3', 't1', 'r:1')], set(), templates=dict(t1=dict(reactions=[dict(source, key='r:1', coefficients=dict(mdh['coefficients']) if False else source['coefficients'])])), mapper=em)['PFK']
    print('INFO evidence summary keys:', sorted(ep.evidence_summary({'PFK': prov}).keys()))
    print('ALL PASS' if ok else 'SOME FAILED')
    raise SystemExit(0 if ok else 1)


if __name__ == '__main__':
    main()
