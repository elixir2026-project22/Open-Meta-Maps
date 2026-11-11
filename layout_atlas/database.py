"""Dual indexes: biochemical observations and core-shape templates with provenance."""
import functools
import gzip
import zlib,itertools,xml.etree.ElementTree as ET
from core import *
from semantics import rhea_index,parse_gpml
from atlas_registry import SCHEMA_VERSION, ensure_schema, pack, unpack, connect_path

def connect():
 """Open the runtime Atlas DB (atlas.sqlite), ensuring schema v1 objects exist."""
 DATA.mkdir(exist_ok=True)
 return connect_path(DATA/'atlas.sqlite')

@functools.lru_cache(maxsize=1)
def _core_modules_cached():
 """Cached open-core modules (immutable during one process)."""
 db=connect();out=tuple(unpack(r[0]) for r in db.execute('SELECT payload FROM modules ORDER BY id'));db.close();return out

def modules(include_brenda=True, include_human_maps=True, include_wikipathways_geometry=True, include_metabolismregulation=True, include_panther=True):
 """Open-core modules from atlas.sqlite, plus default extensions when present.

 WikiPathways geometry-only windows ship as an extension (size) and load by default.
 Metabolism Regulation gallery SBGN-ML and PANTHER CellDesigner maps ship as extensions and load by default.
 """
 out=list(_core_modules_cached())
 for flag, rel in (
  (include_brenda,'extensions/brenda/modules.json.gz'),
  (include_human_maps,'extensions/human_maps/modules.json.gz'),
  (include_wikipathways_geometry,'extensions/wikipathways_geometry/modules.json.gz'),
  (include_metabolismregulation,'extensions/metabolismregulation/modules.json.gz'),
  (include_panther,'extensions/panther/modules.json.gz'),
 ):
  path=DATA/rel
  if flag and path.exists():
   with gzip.open(path,'rt') as f: out.extend(json.load(f))
 return out

def graph(m):
 g=nx.Graph()
 for n in m['nodes']:g.add_node(n['key'],role=n['role'])
 for e in m['edges']:g.add_edge(e['a'],e['b'],role=e.get('flow','unknown'),**{k:v for k,v in e.items() if k not in ('a','b','flow')})
 return g

class Matcher(nx.isomorphism.GraphMatcher):
 def __init__(self,*a,budget=5000,**k):
  super().__init__(*a,**k);self.steps=0;self.budget=budget
 def syntactic_feasibility(self,a,b):
  self.steps+=1
  if self.steps>self.budget:raise TimeoutError('isomorphism budget')
  return super().syntactic_feasibility(a,b)

def _iso_budget(a,b):
 """Scale VF2 step budget with graph size: small exact matches stay thorough; large subgraphs abort sooner."""
 n=max(len(a),len(b))
 if n<=8:return 5000
 if n<=16:return 2500
 return 1200

def mappings(a,b,limit=16):
 try:return list(itertools.islice(Matcher(a,b,node_match=lambda x,y:x['role']==y['role'],budget=_iso_budget(a,b)).isomorphisms_iter(),limit))
 except TimeoutError:return []

def geometry_error(a,b,mp):
 aa={n['key']:n['xy'] for n in a['nodes']};bb={n['key']:n['xy'] for n in b['nodes']};ks=sorted(mp);frame=fit_frame([aa[k] for k in ks],[bb[mp[k]] for k in ks]);f=transform({k:aa[k] for k in ks},frame);return float(np.sqrt(np.mean([np.sum((np.array(f[k])-bb[mp[k]])**2) for k in ks])))

def from_legacy(m):
 core=m.get('core_nodes',m['nodes']);keys={n['key'] for n in core};nodes=[];edges=[]
 for n in core:nodes.append(dict(key=n['key'],role='reaction' if n['key'].startswith('r:') else 'metabolite',xy=n['xy'],entity=n.get('entity'),position_origin='observed',observed=True))
 if m['quality']=='geometry_only':
  # Insert route waypoints for layout transfer, never label them observed reaction centers.
  for i,(a,b) in enumerate(m['edges']):
   rk='route:'+str(i);nodes.append(dict(key=rk,role='reaction',entity=None,xy=((np.array(next(n['xy'] for n in nodes if n['key']==a))+next(n['xy'] for n in nodes if n['key']==b))/2).tolist(),position_origin='derived_route_waypoint',observed=False));edges.extend([dict(a=a,b=rk,flow='unknown'),dict(a=rk,b=b,flow='unknown')])
  reactions=[];size=len(m['edges'])
 else:
  rx={str(r['key']):r for r in m['reactions']};entities={n['key']:n.get('entity') for n in core}
  for a,b in m.get('core_edges',[]):
   rk=a if a.startswith('r:') else b;mk=b if rk==a else a;c=rx[rk[2:]]['coefficients'].get(entities[mk],0);edges.append(dict(a=a,b=b,flow='substrate' if c<0 else 'product',coefficient=c))
  reactions=[dict(r,key='r:'+str(r['key']),references={'bigg.reaction':[r['id']]} if m['source']=='escher' else {'reactome':[r['id']]}) for r in m['reactions']];size=m['size']
 return dict(id=m['id'],source=m['source'],species=m['species'],lineage=m['lineage'],map_id=m['map_id'],name=m['name'],license=m['license'],quality=m['quality'],size=size,nodes=nodes,edges=edges,reactions=reactions,attachments=m['attachments'],ports=m.get('all_participant_ports',m['ports']),legacy_id=m['id'],source_sha256=m['source_sha256'],review_status='not_required',original_instance='data/source_inputs/module_observations.json.gz')

def diagram_modules(d,max_roots=24):
 # d has explicit visible glyphs and semantic or annotated event records.
 g=nx.Graph();rx={r['key']:r for r in d['events']}
 def core_visible(r):
  return [k for k in r['visible'] if d['glyphs'][k].get('label','').lower() not in {'h2o','water','h+','h(+)','atp','adp','amp','nad+','nadh','nadp+','nadph','coa-sh','coa','phosphate','pi','ppi'} and str(d['glyphs'][k]['entity']).split('@')[0] not in {'chebi:15377','chebi:15378','chebi:30616','chebi:456216','chebi:57287','chebi:57288','chebi:57945','chebi:57925'}]
 for k,r in rx.items():
  g.add_node(k,role='reaction')
  for mk in core_visible(r):g.add_node(mk,role='metabolite');g.add_edge(k,mk)
 subsets=set()
 roots=sorted(rx);step=max(1,len(roots)//max_roots)
 for root in roots[::step]:
  order=[k for k in nx.bfs_tree(g,root) if k in rx][:24]
  for n in (1,2,4,8,12,24):
   if len(order)>=n:subsets.add(tuple(sorted(order[:n])))
 out=[]
 for rs in sorted(subsets):
  selected=[rx[k] for k in rs];keys=set(rs)
  for r in selected:keys.update(core_visible(r))
  nodes=[];edges=[]
  for k in sorted(keys):
   if k in rx:nodes.append(dict(key=k,role='reaction',xy=rx[k]['xy'],entity=rx[k]['id'],observed=rx[k]['position_origin']=='observed',position_origin=rx[k]['position_origin']))
   else:
    n=d['glyphs'][k];nodes.append(dict(key=k,role='metabolite',xy=n['xy'],entity=n['entity'],observed=True,position_origin='observed',label=n.get('label','')))
  xy=normalize({n['key']:n['xy'] for n in nodes})
  for n in nodes:n['xy']=xy[n['key']]
  for r in selected:
   for k,v in r['visible'].items():
    if k not in keys:continue
    edges.append(dict(a=r['key'],b=k,flow='substrate' if v<0 else 'product',coefficient=v))
  ports=[]
  for k in keys-set(rs):
   outside=sorted(n for n in g[k] if n not in rs)
   if outside:ports.append(dict(key=k,outside_reactions=outside,scope='available_source_map'))
  out.append(dict(id=uid([d['id'],rs]),source=d['source'],species=d['species'],lineage=d['lineage'],map_id=d['id'],name=d['name'],license=d['license'],quality='semantic_transfer',size=len(rs),nodes=nodes,edges=edges,reactions=selected,attachments=[],ports=ports,source_sha256=d['source_sha256'],review_status='not_required',original_instance=d.get('path'),biochemical_visibility='complete reference reaction; partial observed drawing'))
 return out

def wp_diagram(d,rhea):
 events=[];glyphs={}
 for bind in d['bindings']:
  if bind['status']!='validated_visible_participants':continue
  edge=bind['layout_path'];pts=np.array([[float(p['X']),float(p['Y'])] for p in edge['points']]);lens=np.linalg.norm(np.diff(pts,axis=0),axis=1);half=lens.sum()/2;acc=0;center=pts[0]
  for i,v in enumerate(lens):
   if acc+v>=half:center=pts[i]+(pts[i+1]-pts[i])*((half-acc)/max(v,1e-8));break
   acc+=v
  visible={}
  for key,part in bind['visible_participants'].items():
   gn=d['glyphs'][key];entity='chebi:'+part['chebi']+'@'+part['location'];glyphs[key]=dict(gn,entity=entity);visible[key]=part['coefficient']*bind['orientation']
  rx=bind['reaction'];coeff={k:v*bind['orientation'] for k,v in rx['coefficients'].items()}
  events.append(dict(key='r:'+edge['key'],id='RHEA:'+bind['master'],references={'rhea':[bind['master']]},coefficients=coeff,visible=visible,xy=center.tolist(),position_origin='derived_from_observed_path',reversible=bool(bind['drawing_bidirectional']),direction_evidence=dict(drawing_bidirectional=bind['drawing_bidirectional'],reference_direction=bind['reference_direction'],flux_bounds='not supplied'),path=edge['points'],core=list(visible)))
 return dict(d,glyphs=glyphs,events=events)

def reactome_diagram(path):
 a=read(path);gg=read(path.with_name(path.stem+'.graph.json'));entities={n['dbId']:n for n in gg.get('nodes',[])};events_by_id={e['dbId']:e for e in gg.get('edges',[])};glyphs={};rawnodes={str(n['id']):n for n in a.get('nodes',[])}
 for k,n in rawnodes.items():
  ident=entities.get(n.get('reactomeId'),{}).get('standardIdentifier','');p=n.get('prop',{})
  if ident.startswith('chebi:'):glyphs[k]=dict(entity=ident+'@local',xy=[p.get('x',0)+p.get('width',0)/2,p.get('y',0)+p.get('height',0)/2],label=n.get('displayName',''))
 events=[]
 for e in a.get('edges',[]):
  event=events_by_id.get(e.get('reactomeId'),{});visible={};coeff={};ok=event.get('schemaClass')=='Reaction'
  for role,sign in [('inputs',-1),('outputs',1)]:
   found=[]
   for entry in e.get(role,[]):
    k=str(entry['id'])
    if k not in glyphs:ok=False;continue
    found.append(rawnodes[k]['reactomeId']);matches=[c for c in rawnodes[k].get('connectors',[]) if c['edgeId']==e['id'] and c.get('type')==('INPUT' if sign<0 else 'OUTPUT')];values=[c.get('stoichiometry',{}).get('value') for c in matches]
    if len(values)!=1 or values[0] is None or float(values[0])<=0:ok=False;continue
    visible[k]=sign*float(values[0]);mid=glyphs[k]['entity'];coeff[mid]=coeff.get(mid,0)+visible[k]
   if set(found)!=set(event.get(role,[])):ok=False
  if ok and coeff:
   p=e.get('position',{});events.append(dict(key='r:'+str(e['id']),id=event.get('stId',str(e['reactomeId'])),references={'reactome':[event.get('stId',str(e['reactomeId']))]},coefficients={k:v for k,v in coeff.items() if v},visible=visible,xy=[p.get('x',0),p.get('y',0)],position_origin='observed',reversible=False,paths=e.get('segments',[])))
 rid=path.stem
 return dict(id=uid([rid,sha(path.read_bytes())]),source='reactome',species='Homo sapiens',lineage='reactome:'+rid,map_id=rid,name=a.get('displayName',rid),license='CC0-1.0-data',glyphs=glyphs,events=events,source_sha256=sha(path.read_bytes()),path=str(path.relative_to(BASE)))


def escher_diagram(path, source='escher', license='CC0-1.0', species='unspecified', reaction_namespace='bigg.reaction', **meta):
 """Parse an official Escher JSON map ([metadata, content]) into the diagram form used by diagram_modules.

 Generic for any GEM-aligned Escher map: reaction and metabolite ids are taken from the map's bigg_id fields
 (which may hold BiGG, Human-GEM MAR ids, or other namespaces) and recorded under reaction_namespace.
 """
 path=Path(path); raw=json.loads(path.read_bytes())
 if not isinstance(raw,list) or len(raw)<2: raise ValueError('Escher map must be a [metadata, content] list: '+str(path))
 header,content=raw[0],raw[1]; nodes=content.get('nodes') or {}; reactions=content.get('reactions') or {}; glyphs={}
 for nid,n in nodes.items():
  kind=n.get('node_type')
  if kind=='metabolite':
   entity=n.get('bigg_id') or n.get('name') or str(nid)
   glyphs[str(nid)]=dict(entity=entity,xy=[n.get('x',0),n.get('y',0)],label=n.get('name') or entity,primary=bool(n.get('node_is_primary')),kind='metabolite')
  elif kind in ('midmarker','multimarker'):
   glyphs[str(nid)]=dict(entity=None,xy=[n.get('x',0),n.get('y',0)],label='',kind=kind)
 events=[]; skipped=Counter()
 for rid,r in reactions.items():
  mets=r.get('metabolites') or []
  if not mets: skipped['no_metabolites']+=1; continue
  coeff={m['bigg_id']:float(m['coefficient']) for m in mets if m.get('bigg_id') is not None}
  if not coeff: skipped['no_id_coeff']+=1; continue
  segs=r.get('segments') or {}; used={str(s.get('from_node_id')) for s in segs.values()}|{str(s.get('to_node_id')) for s in segs.values()}
  visible={}
  for nid in used:
   g=glyphs.get(nid)
   if g and g.get('kind')=='metabolite' and g['entity'] in coeff: visible[nid]=coeff[g['entity']]
  if not visible:
   for nid,g in glyphs.items():
    if g.get('kind')=='metabolite' and g['entity'] in coeff: visible[nid]=coeff[g['entity']]
  if not visible: skipped['no_visible_glyphs']+=1; continue
  mids=[nid for nid in used if glyphs.get(nid,{}).get('kind')=='midmarker']
  if mids: xy=glyphs[mids[0]]['xy']
  elif r.get('label_x') is not None: xy=[r['label_x'],r.get('label_y',0)]
  else:
   pts=[glyphs[n]['xy'] for n in visible]; xy=[sum(p[0] for p in pts)/len(pts),sum(p[1] for p in pts)/len(pts)]
  rid_s=str(rid); rx_id=r.get('bigg_id') or r.get('name') or rid_s
  events.append(dict(key='r:'+rid_s,id=rx_id,references={reaction_namespace:[rx_id]},coefficients=coeff,visible=visible,xy=xy,position_origin='observed',reversible=bool(r.get('reversibility')),core=list(visible)))
 digest=sha(path.read_bytes()); name=header.get('map_name') or path.stem
 out=dict(id=uid([source,digest]),source=source,species=species,lineage=source+':'+digest,map_id=header.get('map_id') or digest[:24],name=name,license=license,glyphs=glyphs,events=events,source_sha256=digest,path=str(path),format='escher',audit=dict(skipped))
 out.update(meta); return out


def build():
 rhea=rhea_index();mods=[];original_maps=source_maps();bysha={d['source_sha256']:d for d in original_maps if d['source']=='wikipathways'};maps=[{k:d.get(k) for k in ('id','source','source_id','species','lineage','name','license','source_sha256','provenance_status')} for d in original_maps];binding_audit=[]
 for m in source_modules():
  if m['source']=='reactome':continue  # use the fixed v97 archive, keep v0.2 snapshots only as historical metadata
  if len(m.get('core_nodes',m['nodes']))<=100 and 1<= (m['size'] if m['quality']=='semantic_transfer' else len(m['edges']))<=40:mods.append(from_legacy(m))
 for p in sorted((INPUTS/'raw/wp').rglob('*.gpml')):
  d=parse_gpml(p,rhea)
  if d['source_sha256'] in bysha:
   previous=bysha[d['source_sha256']];d['id']=previous['id'];d['lineage']=previous['lineage']
  if not d['bindings']:continue
  binding_audit.extend(dict(map_id=d['id'],species=d['species'],interaction=b['interaction'],rhea_id=b['rhea_id'],master=b['master'],status=b['status']) for b in d['bindings'])
  dd=wp_diagram(d,rhea);maps.append({k:d[k] for k in ['id','source','source_id','species','lineage','name','author','source_sha256','license']});mods.extend(diagram_modules(dd))
 acq=read(DATA/'acquisition.json')
 for rec in acq['reactome_diagrams']:
  p=BASE/rec['layout_path'] if rec.get('layout_path') else DATA/'reactome'/(rec['id']+'.json')
  try:
   d=reactome_diagram(p);maps.append({k:d[k] for k in ['id','source','species','lineage','name','source_sha256','license']});mods.extend(diagram_modules(d))
  except Exception as e:binding_audit.append(dict(source='reactome',map_id=rec['id'],status='parse_failure',error=str(e)))
 db=connect();db.execute('DELETE FROM reactions');db.execute('DELETE FROM modules');db.execute('DELETE FROM maps');families=defaultdict(list);familyrows=[];variants={};timeouts=0;excluded=[]
 for i,m in enumerate(sorted(mods,key=lambda m:m['id'])):
  if i%3000==0:print('Index',i,'/',len(mods),flush=True)
  original_xy={n['key']:n['xy'] for n in m['nodes']};xy=normalize(original_xy)
  for n in m['nodes']:n['xy']=xy[n['key']]
  m['core_coordinate_normalization']='centered RMS radius; original full instance retained by immutable source pointer'
  gm=graph(m)
  if not len(gm) or not nx.is_connected(gm):
   excluded.append(dict(id=m['id'],source=m['source'],map_id=m['map_id'],reason='empty_or_disconnected_visible_core',original_instance=m.get('legacy_id',m.get('original_instance'))));continue
  h=nx.weisfeiler_lehman_graph_hash(gm,node_attr='role');bucket=(len(gm),gm.number_of_edges(),h);family=None;mp=None
  for f in families[bucket]:
   matches=mappings(gm,f['graph'],1)
   if matches:family=f;mp=matches[0];break
  if family is None:
   family=dict(id=uid(['core-v03',bucket,m['id']]),graph=gm,representative=m,variants=[],instances=[]);families[bucket].append(family);familyrows.append(family)
  variant=None
  for v in family['variants'][:24]:
   matches=mappings(gm,v['graph'],8)
   if matches and min(geometry_error(m,v['representative'],x) for x in matches)<.25:variant=v;break
  if variant is None:variant=dict(id=uid([family['id'],m['id']]),graph=gm,representative=m,instances=[]);family['variants'].append(variant)
  family['instances'].append(m['id']);variant['instances'].append(m['id']);m['core_family']=family['id'];m['variant']=variant['id'];m['biochemical_structure_signature']=uid([(e['a'],e['b'],e.get('flow')) for e in m['edges']]);m['geometry_fingerprint']=uid(np.round(np.sort(np.linalg.norm(np.array([n['xy'] for n in m['nodes']])[:,None]-np.array([n['xy'] for n in m['nodes']])[None,:],axis=2).ravel()),3).tolist())
  db.execute('INSERT OR REPLACE INTO modules VALUES(?,?,?,?,?,?,?,?,?)',(m['id'],m['source'],m['species'],m['lineage'],m['quality'],m['size'],m['core_family'],m['variant'],pack(m)))
  for r in m['reactions']:
   for ns,values in r.get('references',{}).items():
    for val in values:db.execute('INSERT OR IGNORE INTO reactions VALUES(?,?,?,?)',(m['id'],r['key'],ns,str(val)))
 for d in maps:db.execute('INSERT OR REPLACE INTO maps VALUES(?,?,?,?,?)',(d['id'],d['source'],d['species'],d['lineage'],pack(d)))
 db.commit();summary=dict(version='0.3',module_count=db.execute('SELECT COUNT(*) FROM modules').fetchone()[0],by_source={s:n for s,n in db.execute('SELECT source,count(*) FROM modules GROUP BY source')},by_quality={s:n for s,n in db.execute('SELECT quality,count(*) FROM modules GROUP BY quality')},families=len(familyrows),singleton_families=sum(len(f['instances'])==1 for f in familyrows),variants=sum(len(f['variants']) for f in familyrows),rhea_bindings=dict(Counter(b['status'] for b in binding_audit)),rhea_occurrences=sum('rhea_id' in b for b in binding_audit),rhea_unique_ids=len({b['rhea_id'] for b in binding_audit if 'rhea_id' in b}),rhea_maps=len({b['map_id'] for b in binding_audit if 'rhea_id' in b}),rhea_species=len({b['species'] for b in binding_audit if 'rhea_id' in b}),registered_map_records=len({d['id'] for d in maps}),source_maps=db.execute('SELECT count(distinct json_map) FROM (SELECT id AS json_map FROM maps)').fetchone()[0],review_requirement='none',classification_limit='5000 backtracking steps; <=8 geometric mappings; 24 variant representatives; conservative splitting',core_definition='Observed primary core or explicit GPML endpoints; source attachments separately preserved')
 summary['used_source_layouts']=len({m['map_id'] for m in mods if m['id'] not in {x['id'] for x in excluded}})
 dump(OUT/'excluded_modules.json',excluded);summary['excluded_visible_core_instances']=len(excluded)
 dump(OUT/'database_summary.json',summary);dump(OUT/'binding_audit.json',binding_audit);dump(OUT/'library.json',dict(families=[dict(id=f['id'],template=f['representative']['id'],instances=f['instances'],variants=[dict(id=v['id'],template=v['representative']['id'],instances=v['instances']) for v in f['variants']]) for f in familyrows]));db.close();print(summary,flush=True)
 return summary
if __name__=='__main__':build()
