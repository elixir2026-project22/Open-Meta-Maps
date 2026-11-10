"""Source-isolated geometric ranking and same-range whole-model application ablations."""
import gzip,time,pickle,html
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from core import *
import database,retrieval,drawing,pages



def ranking(mods):
 bysource=defaultdict(list)
 for m in sorted(mods,key=lambda m:m['id']):
  if m['size'] in (2,4,8) and 5<=len(m['nodes'])<=35 and sum(n['observed'] for n in m['nodes'])>=3:bysource[m['source']].append(m)
 queries=[]
 for source,vals in sorted(bysource.items()):
  counts=Counter();chosen=[]
  for m in vals:
   if counts[m['lineage']]>=2:continue
   chosen.append(m);counts[m['lineage']]+=1
   if len(chosen)==40:break
  queries.extend(chosen)
 allrows=[];summary=[];rankers={}
 for held in sorted(bysource):
  lib=retrieval.Library([m for m in mods if m['source']!=held]);rows=[]
  for q in queries:
   g=database.graph(q);nr=sum(g.nodes[k]['role']=='reaction' for k in g);ps,a=lib.candidates(g,exclude_lineages=[q['lineage']],exclude_geometry=[retrieval.duplicate_signature(q)],max_roots=1,query_sizes=(nr,));truth={n['key']:n['xy'] for n in q['nodes'] if n['observed']};cs=[]
   for p in ps[:24]:
    keys=sorted(set(truth)&set(p['positions']))
    if len(keys)<3:continue
    frame=fit_frame([p['positions'][k] for k in keys],[truth[k] for k in keys]);fit=transform({k:p['positions'][k] for k in keys},frame);error=float(np.sqrt(np.mean([np.sum((np.array(fit[k])-truth[k])**2) for k in keys])))
    cs.append(dict(template_id=p['template_id'],lineage=p['lineage'],source=p['source'],features=p['features'],score=p['score'],error=error,match_level=p['match_level']))
   rows.append(dict(query=q['id'],map_id=q['map_id'],lineage=q['lineage'],source=q['source'],quality=q['quality'],held_source=held,split='test' if q['source']==held else 'train',candidates=cs))
  train=[r for r in rows if r['split']=='train'];test=[r for r in rows if r['split']=='test'];eligible=[r for r in train if len(r['candidates'])>=2 and max(c['error'] for c in r['candidates'])-min(c['error'] for c in r['candidates'])>=.03];test_multi=sum(len(r['candidates'])>=2 for r in test);ready=len(eligible)>=50 and test_multi>=20;reg=None;linear=None;scale=None;active=[]
  if ready:
   xx=[];yy=[];weights=[];pairs=[];labels=[]
   for q in eligible:
    center=np.mean([c['error'] for c in q['candidates']]);n=len(q['candidates'])
    for c in q['candidates']:xx.append(c['features']);yy.append(center-c['error']);weights.append(1/n)
    for i,a in enumerate(q['candidates']):
     for b in q['candidates'][i+1:]:
      if abs(a['error']-b['error'])<.03:continue
      diff=np.array(a['features'])-b['features'];label=int(a['error']<b['error']);pairs.extend([diff,-diff]);labels.extend([label,1-label])
   reg=HistGradientBoostingRegressor(max_iter=80,max_leaf_nodes=7,min_samples_leaf=8,l2_regularization=2,random_state=SEEDS[0]).fit(xx,yy,sample_weight=weights);rankers[held]=reg
   if pairs:
    active=np.where(np.std(pairs,axis=0)>1e-8)[0].tolist();scale=StandardScaler().fit(np.array(pairs)[:,active]);linear=LogisticRegression(C=1,fit_intercept=False,max_iter=500).fit(scale.transform(np.array(pairs)[:,active]),labels)
   dest=OUT/'rankers';dest.mkdir(parents=True,exist_ok=True)
   with (dest/(held+'.pkl')).open('wb') as f:pickle.dump(reg,f)
  outcomes=defaultdict(list);bymap=defaultdict(lambda:defaultdict(list));wrong=0
  for q in test:
   cs=q['candidates']
   if not cs:continue
   choices=dict(rule=int(np.argmax([c['score'] for c in cs])),oracle=int(np.argmin([c['error'] for c in cs])))
   if reg:choices['gradient_boosting']=int(np.argmax(reg.predict([c['features'] for c in cs])))
   if linear is not None:choices['logistic']=int(np.argmax(linear.decision_function(scale.transform(np.array([c['features'] for c in cs])[:,active]))))
   q['selected']={name:cs[i]['template_id'] for name,i in choices.items()};q['errors']={name:cs[i]['error'] for name,i in choices.items()}
   for name,i in choices.items():outcomes[name].append(cs[i]['error']);bymap[q['map_id']][name].append(cs[i]['error'])
  row=dict(held_source=held,train_queries=len(train),training_preference_queries=len(eligible),test_queries=len(test),test_covered=sum(bool(r['candidates']) for r in test),test_multivariant=test_multi,status='trained' if ready else 'insufficient_preferences',thresholds=dict(training=50,test_multivariant=20),map_macro_errors={name:float(np.mean([np.mean(v[name]) for v in bymap.values()])) for name in outcomes},active_logistic_features=[retrieval.FEATURES[i] for i in active],feature_names=retrieval.FEATURES,ground_truth='Observed source coordinates only; derived reaction centers never labels');summary.append(row);allrows.extend(rows);print('Rank',held,row['status'],len(eligible),'train prefs',test_multi,'test choices',flush=True)
 dump(OUT/'ranking_summary.json',summary)
 with gzip.open(OUT/'ranking_queries.json.gz','wt') as f:json.dump(allrows,f,ensure_ascii=False)
 return rankers

MODELS=[('e_coli_core',REF/'models/e_coli_core.xml','Escherichia coli'),('iML1515',REF/'models/iML1515.xml','Escherichia coli'),('iMM904',REF/'models/iMM904.xml','Saccharomyces cerevisiae'),('Human-GEM',REF/'models/Human-GEM.xml','Homo sapiens'),('iYO844',DATA/'models/iYO844.xml','Bacillus subtilis'),('AraCore_v2_1',DATA/'models/AraCore_v2_1.xml','Arabidopsis thaliana')]

def fixed_views(model,g,max_size=50):
 metas={m['id']:m for m in model['metabolites']};bycomp=defaultdict(list)
 for r in model['reactions']:
  roles=tuple(sorted({metas[m].get('compartment','unknown') for m in r['metabolites']}));bycomp[roles].append(r['id'])
 out=[]
 for roles,ids in sorted(bycomp.items()):
  allowed={'r:'+rid for rid in ids};order=[];visited=set()
  # Incidence traversal produces a fixed, method-independent display scope.
  for root in sorted(allowed):
   if root in visited:continue
   queue=[root];visited.add(root)
   for r in queue:
    order.append(r[2:])
    for mid in sorted(g[r]):
     for nr in sorted(g[mid]):
      if nr in allowed and nr not in visited:visited.add(nr);queue.append(nr)
  for i in range(0,len(order),max_size):out.append(dict(id='view-'+str(len(out)),compartments=list(roles),reaction_ids=order[i:i+max_size]))
 return out

def run_model(name,path,species,mods,condition='all_geometry',max_roots=120,ranker_source=None):
 start=time.perf_counter();model=load_model(path);original=uid(model);g=model_graph(model);library=retrieval.Library([m for m in mods if condition!='semantic_only' or m['quality']=='semantic_transfer']);ps,audit=library.candidates(g,model=model,max_roots=max_roots,query_sizes=(2,4,8,12));ranker=None
 if ranker_source:
  with (OUT/'rankers'/(ranker_source+'.pkl')).open('rb') as f:selector=pickle.load(f)
  def ranker(features):
   features=list(features);chemical=features[1];features[1]=0
   return .6*chemical+float(selector.predict([features])[0])-.1*features[2]-.2*features[3]
 selected=retrieval.select(ps,ranker=ranker);bridge_list=retrieval.bridges(model,selected);anchors,blocks,containers,macro=drawing.pack_modules(model,selected,bridge_list)
 folder=OUT/'models'/name/condition;folder.mkdir(parents=True,exist_ok=True);used={rid for p in selected for rid in p['reaction_ids']};strict={rid for p in selected if p['complete_equation_verified'] for rid in p['reaction_ids']};chemical={rid for p in selected if p.get('complete_chemical_equation_verified',p['complete_equation_verified']) for rid in p['reaction_ids']};refs={rid for p in selected if p['match_level']=='reaction_reference' for rid in p['reaction_ids']};bridgeids={r for b in bridge_list for r in b['reaction_ids']};allids={r['id'] for r in model['reactions']};states={rid:'expert_module' if rid in used else 'algorithmic_bridge' if rid in bridgeids else 'pending' for rid in sorted(allids)}
 # Place bridge centers along their actual network endpoints, never copy template chemistry.
 for b in bridge_list:
  startpos=anchors.get(b['nodes'][0]);endpos=anchors.get(b['nodes'][-1])
  if startpos is None or endpos is None:continue
  for i,rid in enumerate(b['reaction_ids']):anchors['r:'+rid]=(np.array(startpos)+(np.array(endpos)-startpos)*(i+1)/(len(b['reaction_ids'])+1)).tolist()
 manifest=dict(schema=3,model=name,model_id=model['id'],target_file=os.path.relpath(path,ROOT),target_sha256=sha(path.read_bytes()),target_model_hash=original,condition=condition,evaluation='application; all admissible library sources; not an independent generalization score',selector=ranker_source or 'rule',anchors={k:dict(xy=v,hard=False,confidence=.8 if k.startswith('r:') else .4,movable=True) for k,v in anchors.items()},modules=blocks,bridges=bridge_list,containers=containers,macro_routes=macro,reaction_status=states,pending=sorted(allids-used-bridgeids),unmatched_reason={rid:'candidate_conflict' if any(rid in p['reaction_ids'] for p in ps) else 'searched_no_compatible_core' if rid in set(audit['queried_reaction_ids']) else 'outside_current_search_budget' for rid in allids-used-bridgeids},node_copy_policy='module-local copies preserved; shared chemical gets explicit derived junction')
 import interfaces;manifest=interfaces.annotate(model,manifest)
 dump(folder/'skeleton.json',manifest);dump(folder/'retrieval_summary.json',audit)
 with gzip.open(folder/'candidates.json.gz','wt') as f:json.dump(ps,f,ensure_ascii=False)
 views=fixed_views(model,g);rows=[];fail=[]
 for view in views:
  sub=drawing.subset(model,view['reaction_ids']);vg=model_graph(sub);hints={k:v for k,v in anchors.items() if k in vg};xy=drawing.complete(vg,hints,iterations=180);row=drawing.export(sub,xy,folder,view['id'],blocks=blocks,full_model=model,route=False);rows.append(dict(**view,**row))
 # Full small-model controlled comparison includes identical participants and all 95 reactions.
 comparisons=[]
 if len(model['reactions'])<=150:
  for method,hints in [('generic',{}),('skeleton_completed',anchors)]:
   xy=drawing.complete(g,hints,iterations=350);comparisons.append(dict(method=method,**drawing.export(model,xy,folder,method,blocks=blocks if hints else [],full_model=model,route=True)))
 assert uid(model)==original;assert set(r for v in views for r in v['reaction_ids'])==allids
 dump(folder/'views.json',rows);dump(folder/'failures.json',fail)
 (folder/'index.html').write_text(pages.view_index(folder,name,condition,len(allids),len(strict),len(refs),len(used),len(bridgeids),rows))
 result=dict(model=name,condition=condition,reactions=len(allids),strict_equation_reactions=len(strict),chemical_equation_reactions=len(chemical),reference_compatible_reactions=len(refs),module_reactions=len(used),geometry_only_reactions=len(used-chemical-refs),selected_geometry_database_reactions=len({rid for p in selected if p['template_quality']=='geometry_only' for rid in p['reaction_ids']}),selected_core_subgraph_modules=sum(p['match_mode']=='core_subgraph' for p in selected),searched_reaction_union=audit['query_reaction_union'],bridge_reactions=len(bridgeids),modules=len(selected),pending=len(allids-used-bridgeids),views=len(views),candidate_queries=audit['query_windows'],candidate_modules=len(ps),full_small_comparisons=comparisons,seconds=time.perf_counter()-start,view_label_overlaps=sum(r['label_overlap_pairs'] for r in rows),view_node_overlaps=sum(r['node_overlap_pairs'] for r in rows),status='ok');dump(folder/'summary.json',result);print('Model',name,condition,len(used),'/',len(allids),'views',len(views),flush=True);return result

def model_experiments(mods):
 rows=[]
 for name,path,species in MODELS:
  if not path.exists():rows.append(dict(model=name,status='source_unavailable'));continue
  for condition in ('semantic_only','all_geometry'):
   try:rows.append(run_model(name,path,species,mods,condition))
   except Exception as e:rows.append(dict(model=name,condition=condition,status='failed',error=str(e)));print('Model failed',name,condition,str(e),flush=True)
 dump(OUT/'model_summary.json',rows);return rows
if __name__=='__main__':
 mods=database.modules();ranking(mods);model_experiments(mods)
