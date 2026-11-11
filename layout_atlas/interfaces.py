"""Complete target interfaces and confidence, independent of the matched core projection."""
from core import *

def annotate(model,skeleton):
 uses=defaultdict(set);lookup={r['id']:r for r in model['reactions']};metas={m['id']:m for m in model['metabolites']};cur=currency(model);node_sources=defaultdict(list)
 for r in model['reactions']:
  for mid in r['metabolites']:uses[mid].add(r['id'])
 for b in skeleton['modules']:
  rids=set(b['reaction_ids']);mids={mid for rid in rids for mid in lookup[rid]['metabolites']};ports=[]
  for mid in sorted(mids):
   outside=sorted(uses[mid]-rids)
   if outside:ports.append(dict(metabolite=mid,compartment=metas[mid].get('compartment','unknown'),outside_reactions=outside,kind='side_participant' if mid in cur else 'core_or_boundary',xy=b['positions'].get('m:'+mid),position_status='placed_copy' if 'm:'+mid in b['positions'] else 'pending_side_arrangement'))
  b['target_ports']=ports
  for key in b['positions']:node_sources[key].append(b)
 for key,a in skeleton['anchors'].items():
  bs=node_sources.get(key,[]);a['hard']=False;a['movable']=True
  if not bs:a.update(confidence=.25,origin='algorithmic_bridge')
  elif len(bs)>1 and key.startswith('m:'):a.update(confidence=.4,origin='derived_shared_identity_junction',copy_modules=[b['id'] for b in bs])
  else:a.update(confidence=max(b['anchor_confidence'] for b in bs),origin='template_arrangement')
 skeleton['target_interface_scope']='All source-model participants including side species; unmatched chemistry remains target-owned';return skeleton
