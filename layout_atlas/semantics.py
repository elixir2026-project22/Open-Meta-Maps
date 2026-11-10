"""Rhea 143 relationships and loss-aware GPML event binding."""
import gzip,csv,xml.etree.ElementTree as ET
from functools import lru_cache
from core import *
RDF='{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
RH='http://rdf.rhea-db.org/'

@lru_cache(maxsize=1)
def rhea_index():
 path=DATA/'rhea/index.json'
 if path.exists():
  cached=read(path)
  if 'obsolete_ids' not in cached:
   with (DATA/'rhea/rhea-obsoletes.tsv').open() as f:cached['obsolete_ids']=[row['RHEA_ID'] for row in csv.DictReader(f,delimiter='\t')]
   dump(path,cached)
  return cached
 desc=defaultdict(lambda:defaultdict(list))
 needed={'subClassOf','side','contains','compound','coefficient','equation','isTransport','location','curatedOrder'}
 with gzip.open(DATA/'rhea/rhea.rdf.gz') as f:
  for _,e in ET.iterparse(f,events=['end']):
   if e.tag!=RDF+'Description':continue
   subject=e.get(RDF+'about','')
   if subject.startswith(RH):
    key=subject[len(RH):]
    for ch in e:
     name=local(ch)
     if name in needed or re.fullmatch(r'contains\d+',name):
      v=ch.get(RDF+'resource',ch.text or '')
      if v.startswith(RH):v=v[len(RH):]
      desc[key][name].append(v)
   e.clear()
 directions={};dirtype={}
 with (DATA/'rhea/rhea-directions.tsv').open() as f:
  for row in csv.DictReader(f,delimiter='\t'):
   for kind,rid in row.items():directions[rid]=row['RHEA_ID_MASTER'];dirtype[rid]=kind
 ph={}
 with (DATA/'rhea/chebi_pH7_3_mapping.tsv').open() as f:
  for row in csv.DictReader(f,delimiter='\t'):ph[row['CHEBI']]=dict(id=row['CHEBI_PH7_3'],origin=row['ORIGIN'])
 reactions={}
 for rid in sorted(set(directions.values())):
  record=desc[rid];participants=[];coeff={};complete=True
  for side in record.get('side',[]):
   sign=-1 if side.endswith('_L') else 1
   for part in set(desc[side].get('contains',[])):
    pd=desc[part];compound=next(iter(pd.get('compound',[])),None);cd=desc[compound] if compound else {}
    classes=[x for x in cd.get('subClassOf',[]) if 'CHEBI_' in x];chebi=classes[0].rsplit('CHEBI_',1)[-1] if len(classes)==1 else None
    count=[int(n[8:]) for n,values in desc[side].items() if re.fullmatch(r'contains\d+',n) and part in values]
    if not count:
     try:count=[float(pd['coefficient'][0])]
     except (KeyError,ValueError,IndexError):count=[]
    locations=pd.get('location',cd.get('location',[]));location=locations[0] if locations else 'local'
    if chebi is None or len(count)!=1:complete=False
    participants.append(dict(chebi=chebi,coefficient=sign*count[0] if len(count)==1 else None,location=location,source_participant=part))
    if chebi and len(count)==1:
     mid='chebi:'+chebi+'@'+location;coeff[mid]=coeff.get(mid,0)+sign*count[0]
  if not participants:complete=False
  reactions[rid]=dict(id=rid,equation=(record.get('equation') or [''])[0],participants=participants,coefficients={k:v for k,v in coeff.items() if v},complete=complete,is_transport=record.get('isTransport')==['true'])
 out=dict(release=143,master=directions,direction_type=dirtype,ph73=ph,reactions=reactions);dump(path,out);return out

def child(e,name):return next((x for x in e if local(x)==name),None)
def get_chebi(entity):
 if not entity or 'chebi' not in entity.get('Database','').lower():return None
 m=re.search(r'(\d+)$',entity.get('ID',''));return m[1] if m else None

def parse_gpml(path,rhea):
 root=ET.parse(path).getroot();glyphs={};edges=[];anchors={};groups={};enzymes=[]
 for e in root:
  k=e.get('GraphId');tag=local(e);gr=child(e,'Graphics')
  if tag=='Group':groups[k]=dict(e.attrib)
  if tag=='DataNode' and gr is not None:
   xr=child(e,'Xref');entity=dict(xr.attrib) if xr is not None else None
   glyphs[k]=dict(key=k,role=e.get('Type','unknown'),label=e.get('TextLabel',''),xy=[float(gr.get('CenterX',0)),float(gr.get('CenterY',0))],entity=entity,chebi=get_chebi(entity),group=e.get('GroupRef'),position_origin='observed')
 for e in root:
  if local(e) not in ('Interaction','GraphicalLine'):continue
  gr=child(e,'Graphics')
  if gr is None:continue
  points=[dict(x.attrib) for x in gr if local(x)=='Point'];xr=child(e,'Xref');xref=dict(xr.attrib) if xr is not None else None;k=e.get('GraphId',uid(points))
  aa=[]
  for a in gr:
   if local(a)=='Anchor':anchors[a.get('GraphId')]=dict(a.attrib,interaction=k);aa.append(dict(a.attrib))
  if len(points)<2:continue
  edge=dict(key=k,source=points[0].get('GraphRef'),target=points[-1].get('GraphRef'),points=points,xref=xref,anchors=aa,arrows=[p.get('ArrowHead','') for p in (points[0],points[-1])],dashed=gr.get('LineStyle','').lower()=='dashed',graphics=dict(gr.attrib))
  edges.append(edge)
 for edge in edges:
  if edge['source'] in glyphs and glyphs[edge['source']]['role'] in ('GeneProduct','Protein') or edge['target'] in glyphs and glyphs[edge['target']]['role'] in ('GeneProduct','Protein'):enzymes.append(edge)
 bindings=[]
 for edge in edges:
  xr=edge['xref'];isr=xr and xr.get('Database','').lower()=='rhea'
  if not isr:continue
  rid=str(xr.get('ID','')).split(':')[-1];master=rhea['master'].get(rid);reaction=rhea['reactions'].get(master);a=glyphs.get(edge['source']);b=glyphs.get(edge['target']);status='unresolved';orientation=None;matches={}
  if edge['dashed']:status='aggregate_not_promoted'
  elif reaction and a and b and a['role']==b['role']=='Metabolite':
   for orient in (1,-1):
    current={}
    for g,sign in [(a,-orient),(b,orient)]:
     ids={g['chebi']} if g['chebi'] else set();ids|={rhea['ph73'][i]['id'] for i in list(ids) if i in rhea['ph73']}
     poss=[p for p in reaction['participants'] if p['chebi'] in ids and p['coefficient'] is not None and p['coefficient']*sign>0]
     if len(poss)==1:current[g['key']]=poss[0]
    if len(current)==2:orientation=orient;matches=current;break
   status='validated_visible_participants' if orientation and reaction['complete'] else 'reaction_id_only'
  elif master:status='reaction_id_only'
  arrowtypes={'Arrow','mim-conversion','mim-modification'};left=edge['arrows'][0] in arrowtypes;right=edge['arrows'][1] in arrowtypes;edge_direction=-1 if left and not right else 1 if right else None
  if status=='validated_visible_participants':
   if edge_direction is None:status='visual_direction_unresolved'
   else:
    orientation*=edge_direction;expected={'RHEA_ID_LR':1,'RHEA_ID_RL':-1}.get(rhea['direction_type'].get(rid))
    if expected and not(left and right) and orientation!=expected:status='reference_direction_disagrees'
  if rid in rhea.get('obsolete_ids',[]):status='obsolete_reference_not_promoted'
  bindings.append(dict(drawing_bidirectional=left and right,visual_direction=edge_direction,reference_direction=rhea['direction_type'].get(rid),interaction=edge['key'],rhea_id=rid,master=master,status=status,orientation=orientation,visible_participants=matches,reaction=reaction,layout_path=edge,derived_reaction_center=True))
 return dict(id=uid([str(Path('data')/path.relative_to(INPUTS)),sha(path.read_bytes())]),source_id=path.name,name=root.get('Name',path.stem),species=root.get('Organism',path.parent.name.replace('_',' ')),author=root.get('Author',''),lineage=re.search(r'WP\d+',path.name)[0] if re.search(r'WP\d+',path.name) else uid(path.name),glyphs=glyphs,edges=edges,anchors=anchors,groups=groups,enzyme_connections=enzymes,bindings=bindings,source_sha256=sha(path.read_bytes()),license='CC0-1.0',source='wikipathways')

if __name__=='__main__':
 index=rhea_index();print('Rhea',len(index['reactions']),'complete',sum(r['complete'] for r in index['reactions'].values()))
