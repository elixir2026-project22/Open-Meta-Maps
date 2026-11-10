"""Refine chemical-equation evidence into chemistry+known-compartment evidence.
Does not rerun or alter selection/positions; records the stricter reporting distinction.
"""
import gzip,re
from core import *
import interfaces

def run():
 rows=read(OUT/'model_summary.json');extra=[]
 for path in (OUT/'models').glob('*/learned*/summary.json'):extra.append(read(path))
 audits=[]
 for row in rows+extra:
  if row['status']!='ok':continue
  folder=OUT/'models'/row['model']/row['condition'];state=read(folder/'skeleton.json');model=load_model(ROOT/state['target_file']);meta={m['id']:m for m in model['metabolites']}
  with gzip.open(folder/'candidates.json.gz','rt') as f:ps=json.load(f)
  for p in ps:
   for e in p['mapping_evidence']:
    e['compartments_verified']=bool(e['complete_equation']) and all(len(v)==1 and comp(k) not in ('unknown','local') and comp(k)==meta[v[0]].get('compartment') for k,v in e['metabolite_candidates'].items())
   p['complete_chemical_equation_verified']=p['match_level'] in ('complete_equation','reaction_reference') and bool(p['mapping_evidence']) and all(e['complete_equation'] for e in p['mapping_evidence'])
   p['complete_equation_verified']=p['complete_chemical_equation_verified'] and all(e['compartments_verified'] for e in p['mapping_evidence'])
  byid={p['id']:p for p in ps};chem=set();strict=set()
  for b in state['modules']:
   p=byid[b['id']];b['complete_chemical_equation_verified']=p['complete_chemical_equation_verified'];b['complete_equation_verified']=p['complete_equation_verified'];b['anchor_confidence']=.9 if p['complete_equation_verified'] else .7 if p['match_level']=='reaction_reference' else .35
   if p['complete_chemical_equation_verified']:chem.update(b['reaction_ids'])
   if p['complete_equation_verified']:strict.update(b['reaction_ids'])
  refs={rid for b in state['modules'] if b['match_level']=='reaction_reference' for rid in b['reaction_ids']};used={rid for b in state['modules'] for rid in b['reaction_ids']};row['geometry_only_reactions']=len(used-chem-refs)
  row['chemical_equation_reactions']=len(chem);row['strict_equation_reactions']=len(strict);row['compartment_evidence_policy']='Strict requires complete chemical equation and explicit matching source/target compartments; local/unknown never certifies compartment.'
  with gzip.open(folder/'candidates.json.gz','wt') as f:json.dump(ps,f,ensure_ascii=False)
  dump(folder/'skeleton.json',interfaces.annotate(model,state));dump(folder/'summary.json',row)
  index=folder/'index.html'
  if index.exists():index.write_text(re.sub(r'Strict equation coverage \d+',f'Strict equation coverage {len(strict)}',index.read_text()))
  audits.append(dict(model=row['model'],condition=row['condition'],chemical_equation_reactions=len(chem),strict_compartment_equation_reactions=len(strict)))
 dump(OUT/'model_summary.json',rows);dump(OUT/'compartment_mapping_audit.json',audits);return audits
if __name__=='__main__':print(run())
