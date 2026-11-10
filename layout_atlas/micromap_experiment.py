"""Same-mapper source ablation and P0-compatible target-only skeleton exports."""
import copy
import gzip
import time
from core import *
import database
import retrieval
import drawing
import geometry
import semantic_transfer
import micromap_library as ml

CONDITIONS=['G','G_D','G_D_M']


class CanonicalGraph(nx.Graph):
    """Keep seeded fallback layouts independent of Python set iteration."""
    def subgraph(self, nodes):
        keys=set(nodes)&set(self);out=CanonicalGraph()
        out.add_nodes_from((k,dict(self.nodes[k])) for k in sorted(keys))
        edges=sorted(tuple(sorted((a,b))) for a,b in self.edges() if a in keys and b in keys)
        out.add_edges_from((a,b,dict(self.edges[a,b])) for a,b in edges)
        return out


def save_gz(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    with gzip.open(path,'wt') as f:json.dump(obj,f,ensure_ascii=False)


def selected_stats(proposals,allids):
    selected=retrieval.select(copy.deepcopy(proposals))
    covered={r for p in selected for r in p['reaction_ids']}
    chemical={r for p in selected if p.get('complete_chemical_equation_verified',p['complete_equation_verified']) for r in p['reaction_ids']}
    candidate_chemical={r for p in proposals if p.get('complete_chemical_equation_verified',p['complete_equation_verified']) for r in p['reaction_ids']}
    micro={r for p in selected if p['source']=='micromap' for r in p['reaction_ids']}
    return selected,dict(target_reactions=len(allids),candidate_reactions=len({r for p in proposals for r in p['reaction_ids']}),candidate_chemical_reactions=len(candidate_chemical),candidate_chemical_ids=sorted(candidate_chemical),selected_reactions=len(covered),selected_chemical_reactions=len(chemical),selected_chemical_ids=sorted(chemical),selected_micromap_reactions=len(micro),selected_modules=len(selected),pending_before_bridges=len(allids-covered),selected_sources=dict(Counter(p['source'] for p in selected)))


def clip_proposals(selected,sub):
    ids={r['id'] for r in sub['reactions']};g=model_graph(sub);out=[]
    for p in selected:
        rids=sorted(set(p['reaction_ids'])&ids)
        # Single leftover reactions are drawing units, not multi-reaction skeletons.
        if len(rids)<2:continue
        keys={'r:'+r for r in rids}
        for r in rids:keys.update(g['r:'+r])
        q=copy.deepcopy(p);q['reaction_ids']=rids;q['positions']={k:v for k,v in p['positions'].items() if k in keys}
        q['node_mapping']={k:v for k,v in p['node_mapping'].items() if v in q['positions']}
        q['reaction_mapping']={k:v for k,v in p['reaction_mapping'].items() if v in rids}
        q['mapping_evidence']=[e for e in p.get('mapping_evidence',[]) if e.get('target_reaction') in rids or e.get('target_reaction') is None]
        q['source_glyph_mapping']={k:v for k,v in p.get('source_glyph_mapping',{}).items() if 'm:'+v in q['positions']}
        q['view_clipping']=dict(original_proposal=p['id'],full_reaction_ids=p['reaction_ids'],visible_reaction_ids=rids)
        q['id']=uid([p['id'],'fixed-view',rids]);out.append(q)
    return out


def draw_case(model,selected,folder,full_model=None):
    import interfaces
    from geometry_experiment import metric_details,assert_chemistry
    full_model=full_model or model
    original=uid(model);g=CanonicalGraph(model_graph(model));bridges=retrieval.bridges(model,selected)
    anchors,blocks,containers,routes=drawing.pack_modules(model,copy.deepcopy(selected),bridges,physical_scale=True)
    proposal_by_id={p['id']:p for p in selected}
    for b in blocks:
        p=proposal_by_id[b['id']]
        for key in ('license','partition','mapping_evidence','source_glyph_mapping','label_compartments_verified','source_to_initial_frame','source_path_policy','source_citation','view_clipping','source_sha256','original_instance'):
            if key in p:b[key]=p[key]
        if p['source']=='micromap':b['anchor_confidence']=.75
    for bridge in bridges:
        a=anchors.get(bridge['nodes'][0]);b=anchors.get(bridge['nodes'][-1])
        if a is None or b is None:continue
        for i,rid in enumerate(bridge['reaction_ids']):anchors['r:'+rid]=(np.array(a)+(np.array(b)-a)*(i+1)/(len(bridge['reaction_ids'])+1)).tolist()
    core=drawing.complete(g,anchors,iterations=350);scene=drawing.display_scene(model,core,blocks);scene,audit=geometry.complete_scene(model,scene,blocks)
    metric=drawing.export(model,core,folder,'map',blocks=blocks,full_model=full_model,route=True,scene=scene,label_obstacles=True)
    assert metric['routing']['status']=='ok',metric['routing']
    sources=[{k:b.get(k) for k in ('id','template_id','map_id','source','lineage','license','partition','source_citation','source_sha256','original_instance')} for b in blocks]
    dump(folder/'source_licenses.json',dict(layout_sources=sources,policy='Original licenses preserved; MicroMap-derived artifacts remain in the NC-SA research partition'))
    data=read(folder/'map.json');data[0]['metadata']['layout_sources']=sources
    from jsonschema import Draft4Validator
    Draft4Validator(read(ESCHER_SCHEMA)).validate(data);dump(folder/'map.json',data)
    from xml.sax.saxutils import escape
    svg=folder/'map.svg';text=svg.read_text();offset=text.index('>')+1;text=text[:offset]+'<metadata>'+escape(json.dumps(dict(layout_sources=sources),ensure_ascii=False))+'</metadata>'+text[offset:];svg.write_text(text)
    metric.update(metric_details(folder/'map.json'));metric['constraints']=audit;assert_chemistry(model,folder/'map.json');assert uid(model)==original
    used={r for p in selected for r in p['reaction_ids']};bridgeids={r for b in bridges for r in b['reaction_ids']};allids={r['id'] for r in model['reactions']}
    state=dict(schema=4,model=model['id'],target_model_hash=original,modules=blocks,bridges=bridges,containers=containers,
        anchors=scene['anchors'],display_positions=scene['positions'],
        participants=[dict(reaction=r,metabolite=m,glyph=k) for (r,m),k in scene['participants'].items()],
        reaction_status={rid:'expert_module' if rid in used else 'algorithmic_bridge' if rid in bridgeids else 'algorithm_completed' for rid in sorted(allids)},
        pending_before_completion=sorted(allids-used-bridgeids),
        license_policy='Per-module source licenses preserved; any MicroMap-derived display is a separate NC-SA research artifact',
        source_path_policy='Template center-line controls indexed; final connections are algorithmic target-only paths')
    typed_anchors=copy.deepcopy(state['anchors']);state['anchors']={}
    state=interfaces.annotate(full_model,state);state['anchors']=typed_anchors
    state['display_scope_reaction_ids']=sorted(allids);state['full_target_model_hash']=uid(full_model)
    dump(folder/'skeleton.json',state);dump(folder/'metrics.json',metric)
    return metric


def run(model_name=None,draw=True):
    import experiments
    ml.RESULTS.mkdir(parents=True,exist_ok=True)
    summary_path=ml.RESULTS/'library_summary.json'
    if not summary_path.exists() or any(sha((BASE/n).read_bytes())!=h for n,h in read(summary_path)['code_sha256'].items()):micro,summary=ml.build()
    else:micro=ml.load();summary=read(summary_path)
    run_code={n:sha((BASE/n).read_bytes()) for n in ['semantic_transfer.py','micromap_experiment.py','celldesigner.py','micromap_library.py','drawing.py','geometry.py','routing.py']}
    core=database.modules();template_lookup={m['id']:m for m in core+micro};bridge=semantic_transfer.annotation_bridge(experiments.MODELS);dump(ml.RESULTS/'annotation_bridge.json',bridge)
    all_rows=[]
    for name,path,species in experiments.MODELS:
        if model_name and name!=model_name:continue
        start=time.perf_counter();folder=ml.RESULTS/'models'/name;folder.mkdir(parents=True,exist_ok=True)
        model=load_model(path);allids={r['id'] for r in model['reactions']};mapper=semantic_transfer.Mapper(model,bridge)
        print('Direct mapping',name,'core',flush=True);core_proposals,core_audit=mapper.proposals(core,'open_core')
        print('Direct mapping',name,'MicroMap',len(core_proposals),'core candidates',flush=True);micro_proposals,micro_audit=mapper.proposals(micro,'micromap')
        old_folder=OUT/'models'/name/'all_geometry'
        with gzip.open(old_folder/'candidates.json.gz','rt') as f:old=json.load(f)
        for proposal in old+core_proposals+micro_proposals:
            template=template_lookup[proposal['template_id']]
            for key in ('license','source_sha256','original_instance','source_citation'):
                if key in template:proposal[key]=template[key]
            proposal.setdefault('partition',template.get('partition','open_core'))
        pools=[old,old+core_proposals,old+core_proposals+micro_proposals];rows=[];selections={}
        for condition,pool in zip(CONDITIONS,pools):
            selected,stats=selected_stats(pool,allids);selections[condition]=selected
            covered={r for p in selected for r in p['reaction_ids']};candidate_ids={r for p in pool for r in p['reaction_ids']}
            dump(folder/(condition+'-reaction-status.json'),dict(target_model_hash=uid(model),scope='full target, combinatorial selection only',reaction_status={rid:'selected_module_not_yet_drawn' if rid in covered else 'candidate_not_selected' if rid in candidate_ids else 'no_candidate_in_bounded_search' for rid in sorted(allids)}))
            rows.append(dict(condition=condition,**stats));save_gz(folder/(condition+'-selected.json.gz'),selected)
        baseline=set(rows[1]['candidate_chemical_ids']);union=set(rows[2]['candidate_chemical_ids'])
        result=dict(code_sha256=run_code,model=name,target_file=str(path.relative_to(ROOT)),target_sha256=sha(path.read_bytes()),target_model_hash=uid(model),
            conditions=rows,new_candidate_chemical_ids=sorted(union-baseline),
            mapper_core=core_audit,mapper_micro=micro_audit,
            source_model_annotations=bridge['inputs'],core_library_reactions_unchanged=True,
            chemical_mapping_scope='annotation-supported full equations; specific source compartment evidence may be label-only or role-remapped',
            counts_not_independent_designs=True,seconds=time.perf_counter()-start)
        save_gz(folder/'equation_mapping_audit.json.gz',mapper.mapping_rows)
        save_gz(folder/'direct_core_candidates.json.gz',core_proposals);save_gz(folder/'micromap_candidates.json.gz',micro_proposals)
        dump(folder/'mapping_summary.json',result)
        if draw and name=='e_coli_core':
            result['drawing']=[]
            for condition in CONDITIONS:
                print('Draw',name,condition,flush=True);row=draw_case(model,selections[condition],folder/condition)
                result['drawing'].append(dict(condition=condition,**row))
        if draw and name in ('iML1515','iYO844'):
            # Scope fixed from existing target partitions, before inspecting quality.
            views=read(old_folder/'views.json');chosen=sorted({0,len(views)//2,len(views)-1});result['fixed_views']=[]
            for index in chosen:
                view=views[index];sub=drawing.subset(model,view['reaction_ids'])
                for condition in ('G_D','G_D_M'):
                    print('Draw fixed view',name,view['id'],condition,flush=True)
                    selected=clip_proposals(selections[condition],sub);row=draw_case(sub,selected,folder/condition/view['id'],full_model=model)
                    result['fixed_views'].append(dict(view=view['id'],reaction_ids=view['reaction_ids'],condition=condition,**row))
        result['seconds']=time.perf_counter()-start;dump(folder/'summary.json',result);all_rows.append(result)
        print('MicroMap result',name,[(r['condition'],r['selected_reactions'],r['selected_chemical_reactions']) for r in rows],'new candidate chemical',len(result['new_candidate_chemical_ids']),flush=True)
    # Preserve independently generated targets when re-running a single model.
    if model_name:all_rows=[read(p) for p in sorted((ml.RESULTS/'models').glob('*/summary.json')) if read(p).get('code_sha256')==run_code]
    summary=dict(library=summary,models=all_rows,mapper_config=semantic_transfer.CONFIG,
       code_sha256=run_code,
       seed=SEEDS[0],interpretation='Applications, not independent learning tests. G_D versus G_D_M shares mapper and geometry; G versus G_D measures mapping-engine effects.')
    dump(ml.RESULTS/'summary.json',summary)
    import micromap_report;micromap_report.generate(summary)
    return summary

if __name__=='__main__':run()
