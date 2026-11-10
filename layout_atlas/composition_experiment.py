"""Frozen-library bounded composition search and full target-only export checks."""
import copy
import gzip
import shutil
import time
from core import *
import composition as c
import drawing
import micromap_experiment as pipeline

OUTDIR=OUT/'composition'


def gzread(path):
    with gzip.open(path,'rt') as f:return json.load(f)


def adapt_metadata(folder,selected):
    adaptations={p['id']:dict(parent=p.get('composition_parent_id',p['id']),pose=p.get('algorithmic_pose'),
                    source_instance_unchanged=True) for p in selected}
    state=read(folder/'skeleton.json')
    for b in state['modules']:b['composition_adaptation']=adaptations[b['id']]
    dump(folder/'skeleton.json',state)
    data=read(folder/'map.json');data[0]['metadata']['composition_adaptations']=adaptations;dump(folder/'map.json',data)
    # SVG already carries per-source attribution; keep target adaptation with it.
    from xml.sax.saxutils import escape
    svg=folder/'map.svg';text=svg.read_text();marker='</metadata>'
    text=text.replace(marker,escape(json.dumps(dict(composition_adaptations=adaptations),ensure_ascii=False))+marker,1)
    # Use one valid JSON metadata object instead of concatenated JSON documents.
    start=text.index('<metadata>')+len('<metadata>');end=text.index('</metadata>')
    sources=read(folder/'source_licenses.json')['layout_sources']
    text=text[:start]+escape(json.dumps(dict(layout_sources=sources,composition_adaptations=adaptations),ensure_ascii=False))+text[end:];svg.write_text(text)


def run_case(name,view,model,full_model,seed,pool,reference,destination):
    start=time.perf_counter()
    if destination.exists():shutil.rmtree(destination)
    destination.mkdir(parents=True)
    baseline=pipeline.draw_case(model,seed,destination/'seed',full_model=full_model);adapt_metadata(destination/'seed',seed)
    for key in ('node_overlap_pairs','label_overlap_pairs','straight_crossings','visible_nodes','canvas_area_per_reaction'):
        assert abs(baseline[key]-reference[key])<1e-6,('P1 baseline replay mismatch',name,view,key,baseline[key],reference[key])
    seed_state=read(destination/'seed/skeleton.json');context=c.Context(model,full_model,seed_state)
    required,required_chemical=c.coverage(seed)
    raw=c.mutations(seed,pool);raw_generated=len(raw)
    for row in raw:row['approximate']=context.cost(context.approximate_blocks(row['selected']))
    raw=c.diverse_take(raw,c.CONFIG['raw_budget'],'approximate')
    macro_input=c.diverse_take(raw,c.CONFIG['macro_previews'],'approximate')
    for row in macro_input:row['macro']=context.preview(row['selected'])
    finalists=c.diverse_take(macro_input,c.CONFIG['full_layouts']-1,'macro')
    dump(destination/'preview_log.json',dict(raw_generated=raw_generated,retained=len(raw),seed_context=context.cost(seed_state['modules']),previews=[{k:v for k,v in row.items() if k!='selected'}|dict(proposal_ids=[p['id'] for p in row['selected']]) for row in macro_input]))
    baseline_score=c.quality(baseline);best_score=baseline_score;best_folder='seed';best_selection=seed;best_metric=baseline
    trials=[dict(id='seed',kind='seed',folder='seed',metric=baseline,quality=baseline_score,acceptable=True,status='ok')]
    print('Composition',name,view,'raw',len(raw),'macro',len(macro_input),'full',len(finalists)+1,flush=True)
    for index,row in enumerate(finalists,1):
        trial_folder='trials/'+str(index)+'-'+row['id'][:12];folder=destination/trial_folder
        pipeline.save_gz(folder/'selected.json.gz',row['selected'])
        try:
            metric=pipeline.draw_case(model,row['selected'],folder,full_model=full_model);adapt_metadata(folder,row['selected'])
            assert c.valid_selection(row['selected'],required,required_chemical)
            score=c.quality(metric);allowed=c.acceptable(metric,baseline)
            record=dict(id=row['id'],kind=row['kind'],changed_reactions=row['changed_reactions'],folder=trial_folder,
                        approximate=row['approximate'],macro=row['macro'],metric=metric,quality=score,acceptable=allowed,status='ok')
            pipeline.save_gz(folder/'selected.json.gz',row['selected'])
            if allowed and score<best_score-1e-8:best_score=score;best_folder=trial_folder;best_selection=row['selected'];best_metric=metric
            print(' trial',index,row['kind'],'cross',metric['straight_crossings'],'labels',metric['label_overlap_pairs'],'guard',allowed,flush=True)
        except Exception as error:
            record=dict(id=row['id'],kind=row['kind'],folder=trial_folder,status='failed',error=type(error).__name__+': '+str(error))
            print(' trial',index,'failed',record['error'],flush=True)
        trials.append(record);dump(destination/'trial_log.json',trials)
    shutil.copytree(destination/best_folder,destination/'best')
    pipeline.save_gz(destination/'best/selected.json.gz',best_selection)
    pipeline.save_gz(destination/'seed/selected.json.gz',seed)
    after,chem=c.coverage(best_selection)
    result=dict(model=name,view=view,target_model_hash=uid(model),full_target_model_hash=uid(full_model),
        display_reaction_ids=sorted(r['id'] for r in model['reactions']),
        source_selection='P1 G_D_M, frozen candidate libraries; no target expert coordinates',
        required_reactions=sorted(required),required_chemical_reactions=sorted(required_chemical),
        before_covered=len(required),after_covered=len(after),before_chemical=len(required_chemical),after_chemical=len(chem),
        baseline=baseline,best=best_metric,baseline_quality=baseline_score,best_quality=best_score,
        improved=best_folder!='seed',winner=best_folder,
        raw_generated=raw_generated,raw_candidates=len(raw),macro_previews=len(macro_input),full_layout_attempts=len(trials),
        trials=trials,seconds=time.perf_counter()-start)
    dump(destination/'summary.json',result)
    print('Composition result',name,view,baseline['straight_crossings'],'->',best_metric['straight_crossings'],
          baseline['label_overlap_pairs'],'->',best_metric['label_overlap_pairs'],flush=True)
    return result


def run(model_name=None):
    source=OUT/'micromap_extension/summary.json';p1=read(source)
    for filename,digest in p1['code_sha256'].items():assert sha((BASE/filename).read_bytes())==digest,('P1 stale',filename)
    input_hashes={}
    for record in p1['models']:
        base=OUT/'micromap_extension/models'/record['model']
        paths=[ROOT/record['target_file'],base/'G_D_M-selected.json.gz',base/'direct_core_candidates.json.gz',base/'micromap_candidates.json.gz',OUT/'models'/record['model']/'all_geometry/candidates.json.gz']
        for path in paths:input_hashes[str(path.relative_to(ROOT))]=sha(path.read_bytes())
    cases=[];provenance={}
    import sqlite3,database
    source_db=sqlite3.connect(DATA/'atlas.sqlite');extension_db=sqlite3.connect(DATA/'extensions/micromap/atlas.sqlite')
    for row in p1['models']:
        name=row['model']
        if model_name and model_name!=name:continue
        folder=OUT/'micromap_extension/models'/name;model=load_model(ROOT/row['target_file'])
        assert sha((ROOT/row['target_file']).read_bytes())==row['target_sha256']
        seed=gzread(folder/'G_D_M-selected.json.gz')
        with gzip.open(OUT/'models'/name/'all_geometry/candidates.json.gz','rt') as f:original=json.load(f)
        pool=original+gzread(folder/'direct_core_candidates.json.gz')+gzread(folder/'micromap_candidates.json.gz')
        for p in pool:
            if p['template_id'] not in provenance:
                db=extension_db if p['source']=='micromap' else source_db
                item=db.execute('SELECT payload FROM modules WHERE id=?',(p['template_id'],)).fetchone();assert item
                template=database.unpack(item[0]);provenance[p['template_id']]={k:template[k] for k in ('license','partition','source_sha256','original_instance','source_citation') if k in template}
            for key,value in provenance[p['template_id']].items():p.setdefault(key,value)
            p.setdefault('partition','research_extension' if p['source']=='micromap' else 'open_core')
        # Reuse enriched versions for seed/known candidates to preserve attribution.
        enriched={p['id']:p for p in gzread(folder/'G_D_M-selected.json.gz')}
        pool=[enriched.get(p['id'],p) for p in pool]
        if row.get('drawing'):
            ref=next(r for r in row['drawing'] if r['condition']=='G_D_M')
            cases.append(run_case(name,'full',model,model,seed,pool,ref,OUTDIR/name/'full'))
        for ref in row.get('fixed_views',[]):
            if ref['condition']!='G_D_M':continue
            sub=drawing.subset(model,ref['reaction_ids']);local_seed=pipeline.clip_proposals(seed,sub);local_pool=pipeline.clip_proposals(pool,sub)
            cases.append(run_case(name,ref['view'],sub,model,local_seed,local_pool,ref,OUTDIR/name/ref['view']))
    source_db.close();extension_db.close()
    if not cases:raise ValueError('No predefined bounded drawing cases for this model')
    code={f:sha((BASE/f).read_bytes()) for f in ('composition.py','composition_experiment.py')}
    signature=uid([code,c.CONFIG,sha(source.read_bytes()),input_hashes])
    for case in cases:case['experiment_signature']=signature;dump(OUTDIR/case['model']/case['view']/'summary.json',case)
    if model_name:cases=[read(p) for p in sorted(OUTDIR.glob('*/*/summary.json')) if read(p).get('experiment_signature')==signature]
    result=dict(cases=cases,candidate_input_sha256=input_hashes,config=c.CONFIG,code_sha256=code,p1_summary_sha256=sha(source.read_bytes()),
                experiment_signature=signature,seed=SEEDS[0],
                interpretation='Target-only bounded heuristic search, not learned preference or independent readability validation; one source design is not multiplied by its adaptations')
    dump(OUTDIR/'summary.json',result)
    import composition_report;composition_report.generate(result)
    return result

if __name__=='__main__':run()
