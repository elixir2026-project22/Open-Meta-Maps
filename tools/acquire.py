"""Fixed source snapshots and reproducible failure records; no silent release upgrades."""

import sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(__file__).resolve().parents[1]
if str(_ROOT / 'layout_atlas') not in _sys.path:
    _sys.path.insert(0, str(_ROOT / 'layout_atlas'))
if str(_ROOT) not in _sys.path:
    _sys.path.insert(0, str(_ROOT))

import urllib.request,concurrent.futures,zipfile,io,re
from core import *

def fetch(url,path):
 p=DATA/path;p.parent.mkdir(parents=True,exist_ok=True)
 if p.exists():
  data=p.read_bytes();prov=read(p.with_suffix(p.suffix+'.source.json'))
  if sha(data)!=prov['sha256']:raise ValueError('cached source modified: '+str(p))
  return data
 with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'OpenPathwayMaps-research/0.3'}),timeout=40) as r:data=r.read()
 p.write_bytes(data);dump(p.with_suffix(p.suffix+'.source.json'),dict(url=url,sha256=sha(data),bytes=len(data),access_date='2026-11-05'));return data

def acquire():
 freeze();fail=[]
 for name in ['rhea-release.properties','LICENSE.txt','rdf/rhea.rdf.gz','tsv/rhea-directions.tsv','tsv/rhea-obsoletes.tsv','tsv/chebi_pH7_3_mapping.tsv','tsv/rhea2xrefs.tsv']:
  try:
   if name!='LICENSE.txt' and name!='rhea-release.properties':
    version=(DATA/'rhea/rhea-release.properties').read_text()
    if not re.search(r'\b143\b',version):raise ValueError('Expected Rhea 143; change snapshot explicitly')
   data=fetch('https://ftp.expasy.org/databases/rhea/'+name,'rhea/'+Path(name).name);print('Rhea',name,len(data),flush=True)
  except Exception as e:fail.append(dict(source='rhea',file=name,error=str(e)))
 targets=[('iYO844','https://bigg.ucsd.edu/static/models/iYO844.xml'),('AraCore_v2_1','https://raw.githubusercontent.com/pwendering/ArabidopsisCoreModel/b8cf4869172f211393fd6c5b97e58d800da13352/AraCore_v2_1/AraCore_v2_1.xml')]
 for name,url in targets:
  try:fetch(url,'models/'+name+'.xml')
  except Exception as e:fail.append(dict(source=name,error=str(e)))
 try:fetch('https://raw.githubusercontent.com/pwendering/ArabidopsisCoreModel/b8cf4869172f211393fd6c5b97e58d800da13352/LICENSE','models/LICENSE.AraCore')
 except Exception as e:fail.append(dict(source='AraCore_license',error=str(e)))
 # Reuse the exact downloaded archives, preserving licensing separation.
 refs=[]
 for source in ('human_maps','yeast_maps'):
  path=INPUTS/'archives'/(source+'.zip')
  with zipfile.ZipFile(path) as z:
   for name in z.namelist():
    if '/escher/' not in name or not name.endswith('.json'):continue
    raw=z.read(name);dest=DATA/'references'/source/Path(name).name;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(raw)
    refs.append(dict(source=source,path=str(dest.relative_to(BASE)),sha256=sha(raw),archive_sha256=sha(path.read_bytes()),status='mixed_manual_scripted; file_license_unresolved; excluded_from_core_training'))
 try:
  raw=fetch('https://codeload.github.com/SysBioChalmers/Yeast-maps/zip/93b787c6986a834cac323de667aac4b93547a457','references/yeast_original.zip')
  with zipfile.ZipFile(io.BytesIO(raw)) as z:
   files=[n for n in z.namelist() if n.endswith(('SBMLfiles/Yeast8.xml','LICENSE'))]
   for n in files:
    dest=DATA/'references/yeast_original'/Path(n).name;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(z.read(n))
  refs.append(dict(source='yeast_original',files=files,archive_sha256=sha(raw),status='historical CellDesigner source; source lineage and license kept separate'))
 except Exception as e:fail.append(dict(source='yeast_original',error=str(e)))
 try:
  raw=fetch('https://api.github.com/repos/MetabolicAtlas/data-generation/commits?path=maps/RULES.md&per_page=1','references/rules_history.json');refs.append(dict(source='rules_history',commits=json.loads(raw),status='No accessible path history if list empty; README claim not treated as recovered rules'))
 except Exception as e:fail.append(dict(source='rules_history',error=str(e)))
 # Obtain a hierarchy; source bytes freeze the discovery before requesting diagrams.
 urls=['https://reactome.org/ContentService/data/eventsHierarchy/9606','https://reactome.org/download/current/ReactomePathwaysRelation.txt','https://zenodo.org/api/records?q=title%3AReactome&sort=mostrecent&size=1']
 hierarchy=None
 for i,url in enumerate(urls):
  try:
   raw=fetch(url,'reactome/discovery-'+str(i)+'.txt')
   if i==0:
    hierarchy=json.loads(raw)
    if not isinstance(hierarchy,list):raise ValueError('Unexpected eventsHierarchy payload')
  except Exception as e:fail.append(dict(source='reactome_discovery',url=url,error=str(e)))
 records=[]
 if hierarchy:
  ids=set()
  def walk(node,metabolic=False):
   flag=metabolic or node.get('stId')=='R-HSA-1430728'
   if flag and node.get('stId') and node.get('diagram',node.get('hasDiagram',False)):ids.add(node['stId'])
   for ch in node.get('children',[]):walk(ch,flag)
  for n in hierarchy:walk(n)
  dump(DATA/'reactome/diagram_ids.json',sorted(ids))
  def get(rid):
   try:
    a=fetch(f'https://reactome.org/download/current/diagram/{rid}.json','reactome/'+rid+'.json');b=fetch(f'https://reactome.org/download/current/diagram/{rid}.graph.json','reactome/'+rid+'.graph.json');json.loads(a);json.loads(b);return dict(id=rid,layout_sha256=sha(a),graph_sha256=sha(b)),None
   except Exception as e:return None,dict(source='reactome_diagram',id=rid,error=str(e))
  with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
   for ok,err in ex.map(get,sorted(ids)):
    if ok:records.append(ok)
    if err:fail.append(err)
 dump(DATA/'acquisition.json',dict(references=refs,reactome_diagrams=records,failures=fail,reactome_release_archive_status='Discovery snapshot and file hashes; quarterly archive not verified unless recorded explicitly',rhea_expected_release=143))
 print('Acquisition complete:',len(refs),'references',len(records),'Reactome diagrams',len(fail),'failures',flush=True)
 return fail
if __name__=='__main__':acquire()
