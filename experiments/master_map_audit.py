"""Reference-only CellDesigner/master-map audit; never changes the active atlas."""

import sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(__file__).resolve().parents[1]
if str(_ROOT / 'layout_atlas') not in _sys.path:
    _sys.path.insert(0, str(_ROOT / 'layout_atlas'))
if str(_ROOT) not in _sys.path:
    _sys.path.insert(0, str(_ROOT))

import xml.etree.ElementTree as ET,gzip,urllib.request,hashlib,tarfile,io,time
from core import *

MICRO=DATA/'references/micromap'
def canonical(label):
 m=re.fullmatch(r'(.+)\[([^\]]+)\]',label.strip())
 return m[1]+'@'+m[2] if m else None

def parse(path):
 root=ET.parse(path).getroot();species={};aliases={};reactions=[];issues=Counter()
 for e in root.iter():
  tag=local(e)
  if tag=='species' and e.tag.startswith('{http://www.sbml.org/sbml/'):
   species[e.get('id')]=dict(id=e.get('id'),name=e.get('name',''),identity=canonical(e.get('name','')),sbml_compartment=e.get('compartment'))
  elif tag=='speciesAlias':
   b=next((c for c in e if local(c)=='bounds'),None)
   if b is None:issues['alias_without_bounds']+=1;continue
   aliases[e.get('id')]=dict(species=e.get('species'),bounds={k:float(b.get(k,0)) for k in ('x','y','w','h')},observed=True)
 for e in root.iter():
  if local(e)!='reaction' or not e.tag.startswith('{http://www.sbml.org/sbml/'):continue
  coefficients=defaultdict(float);refs=[];ok=True
  for group in e:
   if local(group) not in ('listOfReactants','listOfProducts'):continue
   sign=-1 if local(group)=='listOfReactants' else 1
   for sr in group:
    sid=sr.get('species');raw=sr.get('stoichiometry','1')
    try:v=float(raw)
    except (TypeError,ValueError):ok=False;issues['symbolic_stoichiometry']+=1;continue
    if sid not in species or species[sid]['identity'] is None:ok=False;issues['unresolved_participant_name']+=1;continue
    coefficients[species[sid]['identity']]+=sign*v
  core=[]
  for link in e.iter():
   tag=local(link)
   if tag not in ('baseReactant','baseProduct','reactantLink','productLink'):continue
   sid=link.get('species',link.get('reactant',link.get('product')));aid=link.get('alias');valid=aid in aliases and sid in species and aliases[aid]['species']==sid
   refs.append(dict(alias=aid,species=sid,kind=tag,valid=valid))
   if not valid:issues['invalid_or_mismatched_alias_link']+=1
   elif tag in ('baseReactant','baseProduct'):core.append(aid)
  reversible=e.get('reversible');paths=[dict(attributes=c.attrib,points=c.text or '') for c in e.iter() if local(c)=='editPoints']
  reactions.append(dict(id=e.get('id'),name=e.get('name','').strip(),coefficients={k:v for k,v in coefficients.items() if v},equation_numeric_complete=ok,reversible=(reversible=='true') if reversible is not None else True,reversibility_origin='explicit' if reversible is not None else 'SBML_L2_default_true',glyph_links=refs,core_aliases=sorted(set(core)),paths=paths))
 return dict(species=species,aliases=aliases,reactions=reactions,issues=dict(issues))

def equation_matches(rx,target):
 coeff=target['metabolites'];tc=defaultdict(float)
 for mid,v in coeff.items():tc[chem(mid)+'@'+comp(mid)]+=float(v)
 tc={k:v for k,v in tc.items() if v};sc=rx['coefficients']
 if not rx['equation_numeric_complete'] or set(tc)!=set(sc):return False
 ratios=[tc[k]/sc[k] for k in sc if sc[k]]
 return bool(ratios) and max(ratios)-min(ratios)<1e-8

def acquire():
 MICRO.mkdir(parents=True,exist_ok=True)
 def fetch(url,name,md5=None):
  path=MICRO/name
  if path.exists():raw=path.read_bytes()
  else:
   with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0'}),timeout=90) as response:raw=response.read()
   path.write_bytes(raw)
  if md5:assert hashlib.md5(raw).hexdigest()==md5
  dump(str(path)+'.source.json',dict(url=url,sha256=sha(raw),bytes=len(raw),md5=hashlib.md5(raw).hexdigest(),access_date='2026-11-06'))
  return raw
 meta=json.loads(fetch('https://dataverse.harvard.edu/api/datasets/10740668/versions/4.0','dataset-v4.0.json'))['data']
 for entry in meta['files']:
  f=entry['dataFile']
  if f['filename'] in ('MicroMap.xml','MicroMap Subsystem Overview.xlsx'):fetch('https://dataverse.harvard.edu/api/access/datafile/'+str(f['id']),f['filename'],f['checksum']['value'])
 fetch('https://api.datacite.org/dois/10.7910/dvn/fzkmj8','datacite-master.json')
 snapshot=MICRO/'agora2-class-dataset.json'
 if not snapshot.exists():fetch('https://dataverse.harvard.edu/api/datasets/:persistentId/versions/2.0?persistentId=doi:10.7910/DVN/PAHQH4','agora2-class-dataset.json')
 return meta

def acquire_overlay():
 meta=read(MICRO/'agora2-class-dataset.json')['data'];item=next(f['dataFile'] for f in meta['files'] if f['dataFile']['filename'].endswith('.tar.gz'));path=MICRO/item['filename'];url='https://dataverse.harvard.edu/api/access/datafile/'+str(item['id'])
 if not path.exists():
  with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0'}),timeout=90) as r:raw=r.read()
  path.write_bytes(raw)
 else:raw=path.read_bytes()
 assert hashlib.md5(raw).hexdigest()==item['checksum']['value'];dump(str(path)+'.source.json',dict(url=url,sha256=sha(raw),md5=hashlib.md5(raw).hexdigest(),bytes=len(raw),access_date='2026-11-06'))
 with tarfile.open(path,'r:gz') as tar:
  for member in tar:
   if member.isfile() and member.name.endswith('.xml'):
    data=tar.extractfile(member).read();dest=MICRO/'overlay-example.xml';dest.write_bytes(data);dump(str(dest)+'.source.json',dict(url=url,archive_member=member.name,sha256=sha(data),bytes=len(data),access_date='2026-11-06',layout_parent='MicroMap.xml'));return dest
 raise ValueError('No XML overlay found')

def run(overlay_path=None):
 start=time.perf_counter();d=parse(MICRO/'MicroMap.xml');meta=read(MICRO/'dataset-v4.0.json')['data'];rs=d['reactions'];byname=defaultdict(list)
 for r in rs:byname[r['name']].append(r)
 graph=nx.Graph();reaction_lookup={r['id']:r for r in rs}
 for r in rs:
  rk='r:'+r['id'];graph.add_node(rk,role='reaction')
  for a in r['core_aliases']:graph.add_node('m:'+a,role='metabolite');graph.add_edge(rk,'m:'+a)
 samples=[]
 for rk in sorted(k for k in graph if k.startswith('r:'))[::100]:
  connected=[k for k in nx.bfs_tree(graph,rk) if k.startswith('r:')][:4]
  if len(connected)<2:continue
  keys=set(connected)
  for key in connected:keys.update(graph[key])
  rows=[reaction_lookup[k[2:]] for k in connected];nodes=[]
  for k in keys:
   if k.startswith('m:'):
    a=d['aliases'][k[2:]];b=a['bounds'];nodes.append(dict(key=k,xy=[b['x']+b['w']/2,b['y']+b['h']/2],position_origin='observed_glyph',species=d['species'][a['species']]))
   else:
    aids=reaction_lookup[k[2:]]['core_aliases'];xy=[[d['aliases'][a]['bounds']['x']+d['aliases'][a]['bounds']['w']/2,d['aliases'][a]['bounds']['y']+d['aliases'][a]['bounds']['h']/2] for a in aids];nodes.append(dict(key=k,xy=np.mean(xy,0).tolist(),position_origin='derived_endpoint_mean_not_expert_reaction_center'))
  samples.append(dict(id=uid(['micromap-research',connected]),reaction_ids=connected,nodes=nodes,edges=[list(e) for e in graph.subgraph(keys).edges()],reaction_records=rows,license=meta['license']['name'],lineage=sha((MICRO/'MicroMap.xml').read_bytes()),partition='research_reference_only',admitted_to_core=False))
  if len(samples)>=12:break
 import experiments
 matches=[]
 for name,path,_ in experiments.MODELS:
  model=load_model(path);identified=[];equations=[]
  for target in model['reactions']:
   sources=byname.get(target['id'],[])
   if sources:identified.append(target['id'])
   if any(equation_matches(r,target) for r in sources):equations.append(target['id'])
  state=read(OUT/'models'/name/'all_geometry/skeleton.json');current={rid for b in state['modules'] if b['complete_equation_verified'] for rid in b['reaction_ids']}
  matches.append(dict(model=name,reactions=len(model['reactions']),exact_reaction_label_candidates=len(identified),complete_label_compartment_equation_candidates=len(equations),beyond_current_selected_strict_skeleton=len(set(equations)-current),candidate_ids=identified,equation_ids=equations,not_a_measured_layout_transfer=True))
 overlay=None
 if overlay_path:
  o=parse(overlay_path);shared=set(d['aliases'])&set(o['aliases']);same=sum(d['aliases'][a]['bounds']==o['aliases'][a]['bounds'] and d['aliases'][a]['species']==o['aliases'][a]['species'] for a in shared);overlay=dict(file=str(overlay_path.relative_to(BASE)),member=read(str(overlay_path)+'.source.json').get('archive_member'),aliases=len(o['aliases']),shared_alias_ids=len(shared),same_species_and_bounds=same,shared_reaction_ids=len({r['id'] for r in rs}&{r['id'] for r in o['reactions']}),independent_expert_layout=False)
 names=[s['identity'] for s in d['species'].values() if s['identity']];summary=dict(dataset_doi='10.7910/DVN/FZKMJ8',dataset_version='4.0',license=meta['license'],partition='research_reference_only',admitted_to_core=False,reactions=len(rs),species_records=len(d['species']),unique_label_compartment_ids=len(set(names)),unique_chemical_labels=len({s.split('@')[0] for s in names}),glyph_aliases=len(d['aliases']),invalid_links=d['issues'],complete_numeric_equations=sum(r['equation_numeric_complete'] for r in rs),reaction_names_with_multiple_records=sum(len(v)>1 for v in byname.values()),reversibility_origins=dict(Counter(r['reversibility_origin'] for r in rs)),sbml_compartments=dict(Counter(s['sbml_compartment'] for s in d['species'].values())),label_compartments=dict(Counter(s.rsplit('@',1)[1] for s in names)),reaction_paths_retained=sum(bool(r['paths']) for r in rs),provisional_local_modules=len(samples),model_candidates=matches,overlay_comparison=overlay,seconds=time.perf_counter()-start,limits=['Exact VMH-like reaction labels and metabolite labels only; no forced cross-namespace equivalence.','Label compartments are not real membrane containers.','Derived reaction centers never expert supervision.','NC-SA source remains separate from open core; current atlas numbers unchanged.'])
 dump(OUT/'master_map_audit.json',summary)
 with gzip.open(MICRO/'provisional_modules.json.gz','wt') as f:json.dump(samples,f,ensure_ascii=False)
 print(json.dumps({k:v for k,v in summary.items() if k not in ('model_candidates',)},indent=2),flush=True);print([(x['model'],x['exact_reaction_label_candidates'],x['complete_label_compartment_equation_candidates'],x['beyond_current_selected_strict_skeleton']) for x in matches],flush=True);return summary
if __name__=='__main__':
 acquire();run(acquire_overlay())
