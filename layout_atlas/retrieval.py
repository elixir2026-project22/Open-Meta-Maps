"""Evidence-ranked semantic and geometric retrieval; real incidence, diverse candidates."""
import itertools
from functools import lru_cache
from core import *
import database
from semantics import rhea_index

def duplicate_signature(m):
 return uid([m['geometry_fingerprint'],sorted([sorted((chem(k),float(v)) for k,v in r.get('coefficients',{}).items()) for r in m['reactions']])])

def layout_priority(m):
 """Expert Escher/WikiPathways/Reactome first; possibly-automatic layouts (e.g. BRENDA) rank lower when equivalent."""
 lo=m.get('layout_origin') or {}
 status=lo.get('status') if isinstance(lo,dict) else None
 if status in ('possibly_automatic','unverified','automatic'): return -0.75
 if m.get('source')=='brenda': return -0.75
 return 0.0

FEATURES=['size','biochemical_level','port_mismatch','port_direction_cost','aspect','compactness','axis_fraction','source_support','initial_fit','external_degree_x_aspect','port_fraction_x_spread']

@lru_cache(maxsize=1)
def public_metabolite_links():
 """BiGG metabolite id -> ChEBI ids, from `data/identifier_bridge.json` (explicit annotations of the public BiGG / Human-GEM models, with their origin)."""
 import identifier_bridge
 links=defaultdict(set)
 for bigg,chebis in identifier_bridge.metabolite_chebi().items():links[bigg].update(chebis)
 return links

def target_info(model):
 ri=rhea_index();meta={m['id']:m for m in model['metabolites']};metrefs={}
 for mid,m in meta.items():
  refs=xrefs(m.get('annotations',[]));refs.setdefault('bigg.metabolite',[])
  # BiGG-style ids end with _<compartment>; always treat the chemical base as a BiGG alias (generic, not model-list gated)
  base=chem(mid)
  if base and base!=mid:refs['bigg.metabolite']=sorted(set(refs['bigg.metabolite'])|{base})
  cands=set(refs.get('chebi',[]))
  for alias in refs['bigg.metabolite']:cands.update(public_metabolite_links().get(alias,set()))
  refs['compartment']=m.get('compartment','unknown');refs['name']=m.get('name','')
  refs['chebi_candidate_evidence']='target annotations and explicit public model cross-references; ambiguous sets retained'
  refs['chebi_ph73']=sorted(cands|{ri['ph73'][v]['id'] for v in cands if v in ri['ph73']});metrefs[mid]=refs
 rxrefs={}
 for r in model['reactions']:
  refs=xrefs(r.get('annotations',[]));refs['rhea_master']=sorted({ri['master'].get(v,v) for v in refs.get('rhea',[])})
  if model['id'] in ('e_coli_core','iML1515','iMM904','iYO844'):refs['bigg.reaction']=sorted(set(refs.get('bigg.reaction',[]))|{r['id']})
  rxrefs[r['id']]=refs
 return meta,metrefs,rxrefs

def binding_candidates(entity,target_mids,metrefs,label=None):
 """Target metabolites that can bind a source entity (ChEBI / BiGG / normalised name).

 `label` is used for name-only Atlas nodes (BRENDA ligands, some CellDesigner species) when
 `entity` is not a ChEBI or BiGG id.
 """
 if not isinstance(entity,str) and not label:return []
 if isinstance(entity,str) and entity.startswith('chebi:'):
  c=entity.split('@')[0].split(':')[1];ri=rhea_index();candidates={c,ri['ph73'].get(c,{}).get('id',c)}
  return [mid for mid in target_mids if candidates & set(metrefs[mid].get('chebi_ph73',[]))]
 if isinstance(entity,str):
  if entity in target_mids and chem(entity) in metrefs[entity].get('bigg.metabolite',[]):return [entity]
  options=[mid for mid in target_mids if chem(entity) in metrefs[mid].get('bigg.metabolite',[])]
  matching_comp=[mid for mid in options if comp(entity) not in ('unknown','local') and metrefs[mid].get('compartment')==comp(entity)]
  if matching_comp or options:return matching_comp or options
 # name-only fallback (BRENDA / CellDesigner): compare normalised metabolite names
 name=label or (entity if isinstance(entity,str) and not str(entity).startswith(('lid','r:','m:','chebi:')) else None)
 if name:
  import identifier_bridge as ib
  want=ib.name_keys(name)
  if want:
   return [mid for mid in target_mids if ib.name_keys(metrefs[mid].get('name') or '') & want]
 return []

def windows(g,max_roots=120,sizes=(2,4,8)):
 rs=sorted(k for k in g if g.nodes[k]['role']=='reaction');seen=set();step=max(1,len(rs)//max_roots)
 for root in rs[::step]:
  # Traverse the actual bipartite graph; no artificial chain among neighbors.
  queue=[root];visited={root};order=[]
  for k in queue:
   if g.nodes[k]['role']=='reaction':order.append(k)
   if len(order)>=max(sizes):break
   for n in sorted(g[k]):
    if n not in visited:visited.add(n);queue.append(n)
  for n in sizes:
   if len(order)<n:continue
   selected=tuple(sorted(order[:n]))
   if selected in seen:continue
   seen.add(selected);keys=set(selected)
   for r in selected:keys.update(g[r])
   yield selected,g.subgraph(keys).copy()

class Library:
 def __init__(self,mods):
  self.modules={m['id']:m for m in mods};self.buckets=defaultdict(list);self.biological=defaultdict(set);self.support=defaultdict(set);self.rep={}
  for m in sorted(mods,key=lambda m:m['id']):
   if m['size']<2:continue
   m['duplicate_signature']=duplicate_signature(m)
   self.support[m['core_family']].add(m['lineage'])
   # Up to 24 real variants per topology, interleaving sources for diversity.
   sig=self.signature(database.graph(m));vk=(sig,m['geometry_fingerprint'])
   if vk not in self.rep:
    self.rep[vk]=m;self.buckets[sig].append(m)
   for r in m['reactions']:
    for ns,vals in r.get('references',{}).items():
     for v in vals:self.biological[ns,str(v)].add(m['id'])
  self.partial2=[]
  for sig,vals in self.buckets.items():
   bysource=defaultdict(list)
   for m in vals:bysource[m['source']].append(m)
   self.buckets[sig]=[x for row in itertools.zip_longest(*[bysource[k] for k in sorted(bysource)]) for x in row if x is not None][:48]
   if sig[2]==2 and 4<=sig[0]<=9:
    self.partial2.extend(self.buckets[sig][:3])
 @staticmethod
 def signature(g):return (len(g),g.number_of_edges(),sum(g.nodes[k]['role']=='reaction' for k in g),tuple(sorted(dict(g.degree()).values())))

 def eligible(self,m,exclude_lineages,exclude_species,exclude_geometry):return m['lineage'] not in exclude_lineages and m['species']!=exclude_species and m.get('duplicate_signature',duplicate_signature(m)) not in exclude_geometry

 def candidates(self,g,uses=None,model=None,allow_geometry=True,exclude_lineages=(),exclude_species=None,exclude_geometry=(),max_roots=120,query_sizes=(2,4,8)):
  proposals=[];seen=set();info=target_info(model) if model else None;rlookup={r['id']:r for r in model['reactions']} if model else {};window_count=0;queried_reactions=set()
  # Exact topology gives the common search space; biochemical evidence ranks, not excludes, variants.
  query_list=list(windows(g,max_roots=max_roots,sizes=query_sizes));seen_windows={tuple(rs) for rs,_ in query_list}
  if info:
   reverse=defaultdict(set)
   for rid,refs in info[2].items():
    for ns,vs in refs.items():
     for v in vs:reverse[('rhea' if ns=='rhea_master' else ns),str(v)].add(rid)
   added=0
   for m in sorted(self.modules.values(),key=lambda m:(-m['size'],m['id'])):
    if m['size']<2 or not m['reactions'] or not self.eligible(m,exclude_lineages,exclude_species,exclude_geometry):continue
    ids=[];valid=True
    for sr in m['reactions']:
     choices=set()
     for ns,vs in sr.get('references',{}).items():
      for v in vs:choices.update(reverse[ns,str(v)])
     if len(choices)!=1:valid=False;break
     ids.append('r:'+next(iter(choices)))
    rs=tuple(sorted(ids))
    if not valid or len(set(rs))!=len(rs) or rs in seen_windows:continue
    keys=set(rs)
    for rk in rs:keys.update(g[rk])
    query_list.append((rs,g.subgraph(keys).copy()));seen_windows.add(rs);added+=1
    if added>=600:break
  for rs,tg in query_list:
   window_count+=1;queried_reactions.update(r[2:] for r in rs);sig=self.signature(tg);generic=nx.spring_layout(tg,seed=SEEDS[0],iterations=35,scale=1);targetports={k:len(set(g[k])-set(rs)) for k in tg if g.nodes[k]['role']=='metabolite'};pool=[]
   source_templates=list(self.buckets.get(sig,[]));partial_fallback=bool(model and len(rs)==2 and not source_templates);biological_templates=[]
   if info:
    target_refs=info[2];bioids=set()
    for rk in rs:
     for ns,vals in target_refs[rk[2:]].items():
      for value in vals:bioids.update(self.biological.get(('rhea' if ns=='rhea_master' else ns,str(value)),set()))
    for mid in sorted(bioids):
     m=self.modules[mid]
     if len(m['reactions'])!=len(rs):continue
     if all(any(set(vs)&set(target_refs[rk[2:]].get('rhea_master' if ns=='rhea' else ns,[])) for rk in rs for ns,vs in sr.get('references',{}).items()) for sr in m['reactions']):biological_templates.append(m)
    biological_templates=biological_templates[:24]
   if partial_fallback:source_templates=sorted([m for m in self.partial2 if len(m['nodes'])<=len(tg)],key=lambda m:(-len(m['nodes']),m['id']))[:18]
   source_templates=biological_templates+[m for m in source_templates if m['id'] not in {b['id'] for b in biological_templates}]
   for m in source_templates:
    if not self.eligible(m,exclude_lineages,exclude_species,exclude_geometry) or not allow_geometry and m['quality']=='geometry_only':continue
    sg=database.graph(m);partial_mode=partial_fallback or len(m['nodes'])!=len(tg);direct=None
    if info and m in biological_templates:
     mm={};occupied=set();valid=True
     for sr in m['reactions']:
      choices=[rk for rk in rs if any(set(vs)&set(info[2][rk[2:]].get('rhea_master' if ns=='rhea' else ns,[])) for ns,vs in sr.get('references',{}).items())]
      if len(choices)!=1 or choices[0] in occupied:valid=False;break
      mm[sr['key']]=choices[0];occupied.add(choices[0])
     if valid:
      for node in m['nodes']:
       if node['role']=='reaction':continue
       choices=binding_candidates(node.get('entity'),[k[2:] for k in tg if tg.nodes[k]['role']=='metabolite'],info[1],label=node.get('label'))
       if len(choices)!=1 or 'm:'+choices[0] in occupied:valid=False;break
       mm[node['key']]='m:'+choices[0];occupied.add('m:'+choices[0])
     if valid and all(tg.has_edge(mm[e['a']],mm[e['b']]) for e in m['edges']):direct=mm
    if direct is not None:maps=[direct];partial_mode=True
    elif partial_mode:
     # Skip VF2 when template is much larger than the query window (no useful subgraph match expected)
     if len(sg)>len(tg)*2+2:maps=[]
     else:
      try:maps=[{v:k for k,v in mp.items()} for mp in itertools.islice(database.Matcher(tg,sg,node_match=lambda a,b:a['role']==b['role'],budget=database._iso_budget(tg,sg)).subgraph_isomorphisms_iter(),3)]
      except TimeoutError:maps=[]
    else:maps=database.mappings(sg,tg,3)
    for mp in maps:
     # Align to an algorithmic start only, never the query's expert coordinates.
     source={n['key']:n for n in m['nodes']};ks=sorted(mp);frame=fit_frame([source[k]['xy'] for k in ks],[generic[mp[k]] for k in ks]);positions=transform({mp[k]:source[k]['xy'] for k in ks},frame);rmap={r['key']:mp[r['key']][2:] for r in m['reactions'] if r['key'] in mp};biol=0;mapping_checks=[]
     if info and m['reactions']:
      meta,mrefs,rrefs=info;ok_refs=True;ok_equations=True;any_ref=False
      for sr in m['reactions']:
       tid=rmap.get(sr['key']);tr=rlookup.get(tid)
       if tr is None:ok_refs=False;ok_equations=False;continue
       reference=False
       for ns,vs in sr.get('references',{}).items():
        tns='rhea_master' if ns=='rhea' else ns
        if set(map(str,vs))&set(rrefs[tid].get(tns,[])):reference=True;break
       any_ref|=reference;ok_refs&=reference
       # Entire-equation check uses pairwise candidate identifiers, preserving ambiguity.
       mapping={};used=set()
       for sm in sr['coefficients']:
        options=binding_candidates(sm,tr['metabolites'],mrefs)
        if len(options)!=1 or options[0] in used:ok_equations=False;continue
        mapping[sm]=options[0];used.add(options[0])
       ratio=[]
       if set(mapping)==set(sr['coefficients']) and set(mapping.values())==set(tr['metabolites']):ratio=[tr['metabolites'][tm]/sv for sm,sv in sr['coefficients'].items() if sv for tm in [mapping[sm]]]
       complete=bool(ratio) and max(ratio)-min(ratio)<1e-8
       ok_equations&=complete;mapping_checks.append(dict(source_reaction=sr['id'],target_reaction=tid,reference=reference,complete_equation=complete,compartments_verified=complete and all(comp(sm) not in ('unknown','local') and comp(sm)==mrefs[tm]['compartment'] for sm,tm in mapping.items()),metabolite_candidates={sm:binding_candidates(sm,tr['metabolites'],mrefs) for sm in sr['coefficients']}))
      core_agrees=True
      for sk,tk in mp.items():
       if source[sk]['role']=='metabolite' and isinstance(source[sk].get('entity'),str):
        if tk[2:] not in binding_candidates(source[sk]['entity'],[tk[2:]],mrefs):core_agrees=False
      biol=3 if ok_refs and core_agrees else 2 if ok_equations and core_agrees else 1 if any_ref and core_agrees else 0
     ports=sum(bool(n) for n in targetports.values());srcports=sum(p.get('key') in {n['key'] for n in m['nodes']} for p in m['ports']);xy=np.array(list(positions.values()));extent=np.ptp(xy,axis=0);area=max(float(np.prod(extent)),.05);edgevec=[np.array(positions[a])-positions[b] for a,b in tg.edges() if a in positions and b in positions];axis=np.mean([min(abs(v[0]),abs(v[1]))/max(np.linalg.norm(v),1e-8)<.1 for v in edgevec]);initial=float(np.sqrt(np.mean([np.sum((positions[k]-generic[k])**2) for k in positions])))
     # Edge-role disagreement is geometric context, not a claim of biochemical equivalence.
     flow_cost=0;count=0
     for e in m['edges']:
      a,b=mp[e['a']],mp[e['b']]
      if e.get('flow') in ('substrate','product') and tg.edges[a,b].get('role'):
       flow_cost+=e['flow']!=tg.edges[a,b]['role'];count+=1
     spread=float(np.std(xy[:,0])+np.std(xy[:,1]));features=[len(rs),biol,abs(ports-srcports),flow_cost/max(1,count),float(extent[0]/max(extent[1],.05)),len(tg)/area,float(axis),len(self.support[m['core_family']]),initial,sum(targetports.values())*float(extent[0]/max(extent[1],.05)),ports/max(len(tg),1)*spread]
     p=dict(id=uid([m['id'],rs,sorted(mp.items())]),template_id=m['id'],map_id=m['map_id'],lineage=m['lineage'],source=m['source'],species=m['species'],layout_origin=m.get('layout_origin'),license=m.get('license'),source_citation=m.get('source_citation'),geometry_fingerprint=m['geometry_fingerprint'],reaction_ids=[r[2:] for r in rs],node_mapping=mp,reaction_mapping=rmap,positions={k:np.array(v).tolist() for k,v in positions.items()},position_origins={mp[k]:source[k]['position_origin'] for k in mp},match_level=['geometry','partial_reference','complete_equation','reaction_reference'][biol],complete_chemical_equation_verified=biol>=2 and bool(mapping_checks) and all(x['complete_equation'] for x in mapping_checks),complete_equation_verified=biol>=2 and bool(mapping_checks) and all(x['complete_equation'] and x['compartments_verified'] for x in mapping_checks),mapping_evidence=mapping_checks,features=features,ports=targetports,score=2*len(rs)+.6*biol-.1*features[2]-.2*features[3]-.2*initial+layout_priority(m),expert_query_coordinates_used=False,template_quality=m['quality'],match_mode='core_subgraph' if partial_mode else 'exact_core',unmapped_target_nodes=sorted(set(tg)-set(mp.values())))
     pool.append(p)
   # Distinct geometry, not just distinct source identifiers; bounded to 24 candidates/query.
   distinct={}
   for p in sorted(pool,key=lambda p:(-p['score'],p['id'])):
    ps=np.array([p['positions'][k] for k in sorted(p['positions'])]);key=uid(np.round(np.linalg.norm(ps[:,None]-ps[None,:],axis=2),2).tolist());distinct.setdefault(key,p)
   for p in list(distinct.values())[:24]:
    if p['id'] not in seen:proposals.append(p);seen.add(p['id'])
  return proposals,dict(query_windows=window_count,candidates=len(proposals),levels=dict(Counter(p['match_level'] for p in proposals)),max_candidates_per_query=24,root_budget=max_roots,queried_reaction_ids=sorted(queried_reactions),query_reaction_union=len(queried_reactions))

# Compactness is features[5] = matched_nodes / bbox_area (see candidate scoring above).
# Coverage-safe: among same-source proposals with identical reaction sets,
# keep the best selection_score and break remaining ties by higher compactness.
# Stage B:
# - size-1/2 geometry-only windows stay in the library but are demoted so they
#   cannot displace denser / biochemical covers;
# - complete-equation matches get a score bonus;
# - when many free-standing size≤2 CE crumbs compete with size≥4 non-CE covers,
#   all size≤2 CE crumbs are locked before geometry fill-in (structural pool
#   density only — no source- or model-id branches). Sparse pools keep pure
#   score order so larger CE packs are not fragmented.
GEOMETRY_CRUMB_DEMOTION = 3.0
COMPLETE_EQUATION_BONUS = 4.0
# Escalate to size≤2 CE lock once this many free-standing CE crumbs overlap a
# size≥4 non-CE proposal. Sparse pools (few such conflicts) stay on score order.
CE_CRUMB_LOCK_MIN_COMPETING = 8

def _is_complete_equation(p):
 return bool(p.get('complete_chemical_equation_verified') or p.get('match_level')=='complete_equation')

def _free_standing_ce_crumb(p, ce_sets):
 """Size≤2 CE crumb not nested in any larger CE proposal in the pool."""
 if not (_is_complete_equation(p) and len(p.get('reaction_ids') or [])<=2):
  return False
 s=frozenset(p['reaction_ids'])
 return not any(s<t for t in ce_sets if len(t)>len(s))

def _ce_crumb_overlaps_large_non_ce(p, pool):
 """True when a size≥4 non-CE proposal shares reactions with this crumb."""
 s=set(p['reaction_ids'])
 for q in pool:
  if q is p or _is_complete_equation(q):
   continue
  if len(q.get('reaction_ids') or [])>=4 and set(q['reaction_ids'])&s:
   return True
 return False

def _competing_free_standing_ce_crumbs(pool, ce_sets):
 """Free-standing CE crumbs that share reactions with a size≥4 non-CE cover."""
 return [p for p in pool if _free_standing_ce_crumb(p, ce_sets) and _ce_crumb_overlaps_large_non_ce(p, pool)]

def select(proposals,ranker=None):
 if not proposals: return []
 for p in proposals:
  utility=float(ranker(p['features'])) if ranker else p['score']-2*len(p['reaction_ids'])
  p['selection_score']=2*len(p['reaction_ids'])+utility+layout_priority({'source':p.get('source'),'layout_origin':p.get('layout_origin')})
  biol=float(p['features'][1])
  n_rxn=len(p.get('reaction_ids') or [])
  if biol<1 and n_rxn<=2:
   p['selection_score']-=GEOMETRY_CRUMB_DEMOTION
  if _is_complete_equation(p):
   p['selection_score']+=COMPLETE_EQUATION_BONUS
 best={}; drop=set()
 for p in proposals:
  key=(p.get('source'), frozenset(p['reaction_ids']))
  cur=best.get(key)
  def better(a,b):
   return (a['selection_score'], float(a['features'][5]), a['id'])>(b['selection_score'], float(b['features'][5]), b['id'])
  if cur is None or better(p,cur):
   if cur is not None: drop.add(id(cur))
   best[key]=p
  else:
   drop.add(id(p))
 pool=[p for p in proposals if id(p) not in drop]
 ce_sets=[frozenset(p['reaction_ids']) for p in pool if _is_complete_equation(p)]
 # Dense CE-crumb / large-geometry competition → lock all size≤2 CE covers.
 # Sparse pools → score order (avoids fragmenting larger CE packs).
 competing=_competing_free_standing_ce_crumbs(pool, ce_sets)
 if len(competing)>=CE_CRUMB_LOCK_MIN_COMPETING:
  lock=[p for p in pool if _is_complete_equation(p) and len(p.get('reaction_ids') or [])<=2]
 else:
  lock=[]
 rest=[p for p in pool if p not in lock]
 chosen=[];used=set()
 for p in sorted(lock,key=lambda p:(-p['selection_score'],p['id']))+sorted(rest,key=lambda p:(-p['selection_score'],-len(p['reaction_ids']),p['id'])):
  if set(p['reaction_ids'])&used:continue
  chosen.append(p);used.update(p['reaction_ids'])
 # One-for-up-to-two local exchange, preserving all internal reaction exclusions.
 # Never let geometry-only replace complete-equation covers in this exchange.
 for p in sorted(pool,key=lambda p:(-p['selection_score'],-len(p['reaction_ids']),p['id'])):
  conflicts=[c for c in chosen if set(c['reaction_ids'])&set(p['reaction_ids'])]
  if not (1<=len(conflicts)<=2 and p['selection_score']>sum(c['selection_score'] for c in conflicts)+.05):
   continue
  p_geo=float(p['features'][1])<1 and not _is_complete_equation(p)
  if p_geo and any(_is_complete_equation(c) for c in conflicts):
   continue
  chosen=[c for c in chosen if c not in conflicts]+[p]
 return chosen


def bridges(model,chosen,max_new=4):
 g=model_graph(model);owner={('r:'+rid):i for i,p in enumerate(chosen) for rid in p['reaction_ids']};paths=[];seen=set()
 for root in sorted(owner)[::max(1,len(owner)//100)]:
  q=[(root,[root],0)];visited={root}
  for k,path,nnew in q:
   if k in owner and owner[k]!=owner[root]:
    new=[v[2:] for v in path if v.startswith('r:') and v not in owner]
    if new and tuple(new) not in seen:paths.append(dict(reaction_ids=new,nodes=path,from_module=owner[root],to_module=owner[k],origin='target_network_bridge'));seen.add(tuple(new))
    break
   for n in sorted(g[k]):
    nn=nnew+int(n.startswith('r:') and n not in owner)
    if n not in visited and nn<=max_new:visited.add(n);q.append((n,path+[n],nn))
 out=[];uf=nx.utils.UnionFind(range(len(chosen)));occupied=set(owner);added=set()
 for p in sorted(paths,key=lambda p:(len(p['reaction_ids']),p['reaction_ids'])):
  a,b=p['from_module'],p['to_module']
  if uf[a]==uf[b] or set(p['reaction_ids'])&added:continue
  uf.union(a,b);added.update(p['reaction_ids']);out.append(p)
 return out


def clip_proposals(selected, sub):
    """Clip expert proposals to a reaction subset (fixed views). Moved from micromap_experiment."""
    import copy
    ids = {r['id'] for r in sub['reactions']}
    g = model_graph(sub)
    out = []
    for p in selected:
        rids = sorted(set(p['reaction_ids']) & ids)
        if len(rids) < 2:
            continue
        keys = {'r:' + r for r in rids}
        for r in rids:
            keys.update(g['r:' + r])
        q = copy.deepcopy(p)
        q['reaction_ids'] = rids
        q['positions'] = {k: v for k, v in p['positions'].items() if k in keys}
        q['node_mapping'] = {k: v for k, v in p['node_mapping'].items() if v in q['positions']}
        q['reaction_mapping'] = {k: v for k, v in p['reaction_mapping'].items() if v in rids}
        q['mapping_evidence'] = [e for e in p.get('mapping_evidence', []) if e.get('target_reaction') in rids or e.get('target_reaction') is None]
        q['source_glyph_mapping'] = {k: v for k, v in p.get('source_glyph_mapping', {}).items() if 'm:' + v in q['positions']}
        q['view_clipping'] = dict(original_proposal=p['id'], full_reaction_ids=p['reaction_ids'], visible_reaction_ids=rids)
        q['id'] = uid([p['id'], 'fixed-view', rids])
        out.append(q)
    return out
