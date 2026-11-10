"""Bounded predeclared representative sensitivity, not a full-library rerun."""
from core import *
import database

def run():
 byfamily=defaultdict(list)
 for m in database.modules():byfamily[m['core_family']].append(m)
 selected=[sorted(v,key=lambda m:m['id'])[:8] for k,v in sorted(byfamily.items()) if len(v)>=3][:100];counts={str(t):0 for t in (.15,.25,.35)};invariance=[]
 for vals in selected:
  distances={}
  for i,a in enumerate(vals):
   b=json.loads(json.dumps(a));xy={n['key']:n['xy'] for n in a['nodes']};rot=np.array([[0,-1],[-1,0]])
   for n in b['nodes']:n['xy']=(np.array(n['xy'])@rot*3+[14,-21]).tolist()
   invxy=normalize({n['key']:n['xy'] for n in b['nodes']})
   for n in b['nodes']:n['xy']=invxy[n['key']]
   invariance.append(database.geometry_error(a,b,{k:k for k in xy}))
   for j,c in enumerate(vals[:i]):
    mp=database.mappings(database.graph(a),database.graph(c),8);distances[i,j]=min([database.geometry_error(a,c,m) for m in mp] or [float('inf')])
  for threshold in (.15,.25,.35):
   reps=[]
   for i in range(len(vals)):
    if not any(distances[i,j]<threshold for j in reps):reps.append(i)
   counts[str(threshold)]+=len(reps)
 result=dict(families=len(selected),instances=sum(map(len,selected)),selection='First100 family IDs with >=3 instances; first8 instance IDs each; <=8 isomorphisms',threshold_variant_counts=counts,rigid_scale_reflection_max_error=max(invariance),note='Bounded sensitivity sample, not replacement full-library counts; supports auditing threshold dependence, not human correctness.');dump(OUT/'variant_sensitivity.json',result);return result
if __name__=='__main__':print(run())
