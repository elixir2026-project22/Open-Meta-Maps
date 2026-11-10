"""Target-model B1 -> source-held-out EGNN deployment, without target coordinate labels."""
import time,torch,copy
from core import *
import drawing,completion,neural

def run(name='e_coli_core',condition='all_geometry',held='community',seed=SEEDS[0]):
 folder=OUT/'models'/name/condition;skeleton=read(folder/'skeleton.json');model=load_model(ROOT/skeleton['target_file']);metas={m['id']:m for m in model['metabolites']};net=neural.EGNN();checkpoint=OUT/'egnn'/f'{held}-{seed}.pt';net.load_state_dict(torch.load(checkpoint,map_location='cpu',weights_only=True));net.eval();rows=[];start=time.perf_counter()
 for view in read(folder/'views.json'):
  sub=drawing.subset(model,view['reaction_ids']);g=model_graph(sub);raw={k:v['xy'] for k,v in skeleton['anchors'].items() if k in g}
  if len(raw)<2:rows.append(dict(view=view['id'],status='insufficient_anchors'));continue
  frame=completion.coordinate_frame(raw);hints=completion.to_local(raw,frame);nodes=[dict(key=k,**v) for k,v in g.nodes(data=True)];edges=[dict(a=a,b=b,coefficient=v.get('coefficient',0),flow=v.get('role','unknown')) for a,b,v in g.edges(data=True)]
  # No labels, expert geometry, or hidden coordinates exist in the deployment input.
  task=dict(input=dict(nodes=nodes,edges=edges,hints=hints,hard_anchors={},confidence={k:skeleton['anchors'][k]['confidence'] for k in raw}))
  out=neural.predict(net,task);xy={k:((np.array(v)/100*frame['scale'])@np.array(frame['rotation']).T+frame['center']).tolist() for k,v in out.items()};blocks=copy.deepcopy(skeleton['modules'])
  for block in blocks:
   for k,pos in block['positions'].items():
    if k in xy and k in raw:block['positions'][k]=(np.array(pos)+np.array(xy[k])-raw[k]).tolist()
  dump(folder/(view['id']+'-egnn-input.json'),task)
  dump(folder/(view['id']+'-egnn-positions.json'),xy)
  result=drawing.export(sub,xy,folder,view['id']+'-egnn',blocks=blocks,full_model=model,route=False);rows.append(dict(view=view['id'],status='ok',**result))
 result=dict(model=name,condition=condition,checkpoint=str(checkpoint.relative_to(BASE)),checkpoint_sha256=sha(checkpoint.read_bytes()),seed=seed,status='ok' if any(r['status']=='ok' for r in rows) else 'no_eligible_view',views=rows,seconds=time.perf_counter()-start,interpretation='Deployment feasibility on a known model; not independent evidence of full-model learning gain. Remaining side participants/labels use deterministic geometry rules.')
 dump(folder/'egnn_deployment.json',result);return result
if __name__=='__main__':print(run())
