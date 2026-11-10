"""Independent persisted-artifact checks for the separate MicroMap extension."""
import gzip
from core import *
import database


def verify(check_export,core_lookup):
    folder=OUT/'micromap_extension';summary=read(folder/'summary.json');library=summary['library'];checks=Counter()
    for name,digest in summary['code_sha256'].items():assert sha((BASE/name).read_bytes())==digest,('stale extension results',name)
    assert library['independent_designs']==1 and not library['admitted_to_open_core']
    assert library['source_sha256']==sha((DATA/'references/micromap/MicroMap.xml').read_bytes())
    ext=DATA/'extensions/micromap'
    with gzip.open(ext/'observations.json.gz','rt') as f:observations=json.load(f)
    with gzip.open(ext/'modules.json.gz','rt') as f:modules=json.load(f)
    assert len(modules)==library['modules'];lookup={m['id']:m for m in modules}
    assert not set(lookup)&set(core_lookup)
    all_uses=defaultdict(set)
    for r in observations['raw_reactions']:
        for link in r['glyph_links']:
            if link['valid']:all_uses[link['alias']].add('r:'+r['id'])
    for m in modules:
        assert m['source']=='micromap' and not m['admitted_to_open_core'] and m['license']=='CC-BY-NC-SA-4.0'
        assert m['source_citation']['doi']=='10.7910/DVN/FZKMJ8' and m['source_citation']['authors']
        normalized=normalize(m['original_coordinates']);g=database.graph(m);assert nx.is_connected(g)
        selected={r['key'] for r in m['reactions']}
        for n in m['nodes']:
            assert np.allclose(n['xy'],normalized[n['key']])
            assert n['observed']==(n['role']=='metabolite')
            if n['observed']:assert np.allclose(m['original_coordinates'][n['key']],observations['glyphs'][n['key']]['xy'])
        ports={p['key']:p for p in m['ports']}
        for n in m['nodes']:
            if n['role']=='reaction':continue
            expected=all_uses[n['key']]-selected
            assert set(ports.get(n['key'],{}).get('outside_reactions',[]))==expected
            checks['micromap_source_ports_checked']+=1
        checks['micromap_module_records']+=1
    db=sqlite3.connect(ext/'atlas.sqlite');assert db.execute('SELECT COUNT(*) FROM modules').fetchone()[0]==len(modules)
    for mid,payload in db.execute('SELECT id,payload FROM modules'):assert database.unpack(payload)==lookup[mid]
    applicability=read(folder/'applicability_index.json');assert applicability['expert_master_designs']==1
    assert db.execute('SELECT COUNT(*) FROM model_applicability').fetchone()[0]==sum(r['candidates'] for r in applicability['by_target'].values())
    assert db.execute('SELECT COUNT(*) FROM target_models').fetchone()[0]==len(summary['models'])
    assert not db.execute('SELECT a.module_id FROM model_applicability a LEFT JOIN modules m ON a.module_id=m.id WHERE m.id IS NULL').fetchall()
    checks['extension_indexed_applicability']=db.execute('SELECT COUNT(*) FROM model_applicability').fetchone()[0];db.close()
    combined={**core_lookup,**lookup}
    for model in summary['models']:
        assert model['code_sha256']==summary['code_sha256'];target=load_model(ROOT/model['target_file']);rlookup={r['id']:r for r in target['reactions']};g=model_graph(target)
        assert uid(target)==model['target_model_hash'];assert sha((ROOT/model['target_file']).read_bytes())==model['target_sha256']
        result_folder=folder/'models'/model['model'];uses=defaultdict(set)
        for r in target['reactions']:
            for mid in r['metabolites']:uses[mid].add(r['id'])
        for cond in model['conditions']:
            with gzip.open(result_folder/(cond['condition']+'-selected.json.gz'),'rt') as f:selected=json.load(f)
            used=[]
            for p in selected:
                assert p['expert_query_coordinates_used'] is False;template=combined[p['template_id']];snodes={n['key']:n for n in template['nodes']}
                for sk,tk in p['node_mapping'].items():assert snodes[sk]['role']==g.nodes[tk]['role']
                for e in template['edges']:assert g.has_edge(p['node_mapping'][e['a']],p['node_mapping'][e['b']])
                if p['match_mode']=='direct_full_equation_partial_visible_core':
                    sr={r['key']:r for r in template['reactions']}
                    for e in p['mapping_evidence']:
                        coeff=sr[e['source_key']]['coefficients'];met=e['metabolites'];tr=rlookup[e['target_reaction']]
                        assert set(coeff)==set(met) and set(met.values())==set(tr['metabolites'])
                        assert len(met.values())==len(set(met.values()))
                        assert all(abs(v*e['factor']-tr['metabolites'][met[k]])<1e-8 for k,v in coeff.items())
                        checks['extension_equation_mappings']+=1
                if p['source']=='micromap':
                    assert not p['complete_equation_verified'];assert p['source_citation']['authors']
                used.extend(p['reaction_ids'])
            assert len(used)==len(set(used));assert len(used)==cond['selected_reactions']
            status=read(result_folder/(cond['condition']+'-reaction-status.json'))
            assert set(status['reaction_status'])==set(rlookup)
            assert {k for k,v in status['reaction_status'].items() if v=='selected_module_not_yet_drawn'}==set(used)
            checks['extension_full_target_accounting']+=1
        cases=[(r['condition'],r,None) for r in model.get('drawing',[])]+[(r['condition'],r,r['view']) for r in model.get('fixed_views',[])]
        for cond,row,view in cases:
            dest=result_folder/cond if view is None else result_folder/cond/view
            ids=set(rlookup) if view is None else set(row['reaction_ids']);check_export(dest/'map.json',ids,target)
            state=read(dest/'skeleton.json');assert set(state['reaction_status'])==ids;assert state['full_target_model_hash']==uid(target)
            for b in state['modules']:
                assert set(b['reaction_ids'])<=ids
                assert set(b['reaction_mapping'].values())<=ids
                for port in b['target_ports']:assert port['outside_reactions']==sorted(uses[port['metabolite']]-set(b['reaction_ids']));checks['extension_full_model_ports']+=1
            for k,hint in state['anchors'].items():
                delta=np.abs(np.array(state['display_positions'][k])-hint['xy']);assert np.max(delta)<1e-7 if hint['hard'] else np.max(delta)<=hint['allowed_motion']+1e-7
            data=read(dest/'map.json');assert 'layout_sources' in data[0]['metadata'];assert (dest/'source_licenses.json').exists()
            checks['extension_executed_views']+=1
        checks['extension_target_models']+=1
    return checks
