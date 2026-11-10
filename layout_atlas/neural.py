"""Small 2-D EGNN residual completion; prediction only receives graph and hints."""
import time,random
import torch
from torch import nn
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from core import *
import completion,drawing

torch.set_num_threads(1)

def tensors(task):
 ip=task['input'];g=completion.task_graph(task);keys=sorted(g);idx={k:i for i,k in enumerate(keys)};initial=drawing.initial(g,ip['hints'],SEEDS[0]);features=[]
 for k in keys:
  n=g.nodes[k];c=[0.]*4
  if n.get('compartment') not in ('unknown','local',None):c[int(uid(n['compartment'])[:8],16)%4]=1
  features.append([float(n['role']=='reaction'),float(n['role']=='metabolite'),float(n.get('reversible',False)),len(g[k])/4,ip['confidence'].get(k,0),float(k in ip['hard_anchors']),*c])
 edges=[];attrs=[]
 for a,b,data in g.edges(data=True):
  c=float(data.get('coefficient',0));feat=[float(data.get('flow')=='substrate'),float(data.get('flow')=='product'),math.log1p(abs(c))]
  for x,y in [(a,b),(b,a)]:edges.append([idx[x],idx[y]]);attrs.append(feat)
 h=torch.tensor(features,dtype=torch.float32);x=torch.tensor([initial[k] for k in keys],dtype=torch.float32);e=torch.tensor(edges,dtype=torch.long).T;ea=torch.tensor(attrs,dtype=torch.float32);hard=torch.tensor([k in ip['hard_anchors'] for k in keys],dtype=torch.bool)
 return keys,h,x,e,ea,hard

class Layer(nn.Module):
 def __init__(self,hidden=64):
  super().__init__();self.edge=nn.Sequential(nn.Linear(hidden*2+4,hidden),nn.SiLU(),nn.Linear(hidden,hidden),nn.SiLU());self.coord=nn.Sequential(nn.Linear(hidden,hidden),nn.SiLU(),nn.Linear(hidden,1),nn.Tanh());self.node=nn.Sequential(nn.Linear(hidden*2,hidden),nn.SiLU(),nn.Linear(hidden,hidden))
 def forward(self,h,x,e,ea,hard,fixed):
  a,b=e;delta=(x[a]-x[b])/100;d2=(delta*delta).sum(1,keepdim=True);m=self.edge(torch.cat([h[a],h[b],d2,ea],1));weight=self.coord(m);move=delta/(1+torch.sqrt(d2+1e-9))*weight*32
  dx=torch.zeros_like(x);dx.index_add_(0,a,move);degree=torch.zeros((len(h),1),device=h.device);degree.index_add_(0,a,torch.ones((len(a),1),device=h.device));agg=torch.zeros_like(h);agg.index_add_(0,a,m)
  x=torch.where(hard[:,None],fixed,x+dx/degree.clamp_min(1));h=h+self.node(torch.cat([h,agg/degree.clamp_min(1)],1));return h,x

class EGNN(nn.Module):
 def __init__(self):super().__init__();self.embedding=nn.Linear(10,64);self.layers=nn.ModuleList([Layer(64) for _ in range(4)])
 def forward(self,h,x,e,ea,hard):
  h=self.embedding(h);fixed=x.clone()
  for layer in self.layers:h,x=layer(h,x,e,ea,hard,fixed)
  return x

class OffsetRidge:
 def __init__(self):self.reg=make_pipeline(StandardScaler(),Ridge(alpha=10))
 def features(self,task):
  g=completion.task_graph(task);ip=task['input'];xy=drawing.initial(g,ip['hints']);keys=sorted(set(g)-set(ip['hard_anchors']));xx=[]
  for k in keys:
   ns=list(g[k]);v=np.mean([xy[n] for n in ns],0)-xy[k] if ns else np.zeros(2);hints=[ip['hints'][n] for n in ns if n in ip['hints']];d=np.mean(hints,0)-xy[k] if hints else np.zeros(2);xx.append([float(g.nodes[k]['role']=='reaction'),len(ns),len(hints),*v,*d,ip['confidence'].get(k,0)])
  return keys,np.array(xx),xy
 def fit(self,tasks):
  xx=[];yy=[]
  for t in tasks:
   keys,x,xy=self.features(t);labels=t['labels']['observed_xy']
   for k,v in zip(keys,x):
    if k in t['hidden']:xx.append(v);yy.append(np.array(labels[k])-xy[k])
  self.reg.fit(xx,yy);return self
 def predict(self,task):
  keys,x,xy=self.features(task)
  if len(keys):
   for k,v in zip(keys,self.reg.predict(x)):xy[k]=(np.array(xy[k])+v).tolist()
  return completion.constrain(completion.task_graph(task),xy,task['input']['hints'],task['input']['hard_anchors'],task['input']['confidence'])

def predict(model,task,postprocess=True):
 keys,h,x,e,ea,hard=tensors(task)
 with torch.no_grad():out=model(h,x,e,ea,hard).numpy()
 xy={k:v.tolist() for k,v in zip(keys,out)}
 return completion.constrain(completion.task_graph(task),xy,task['input']['hints'],task['input']['hard_anchors'],task['input']['confidence']) if postprocess else xy

def train_fold(tasks,held,seed):
 train=[t for t in tasks if t['held_source']==held and t['split']=='train'];test=[t for t in tasks if t['held_source']==held and t['split']=='test'];maps=sorted({t['lineage'] for t in train});valmaps=set(maps[::max(1,len(maps)//4)][:max(1,len(maps)//5)]);training=[t for t in train if t['lineage'] not in valmaps];validation=[t for t in train if t['lineage'] in valmaps]
 if len(training)<10 or not validation or not test:return None,dict(held=held,seed=seed,status='insufficient_independent_tasks'),None
 random.seed(seed);np.random.seed(seed);torch.manual_seed(seed);model=EGNN();optimizer=torch.optim.Adam(model.parameters(),lr=.001,weight_decay=1e-4);cache={t['id']:tensors(t) for t in training+validation};best=None;bestloss=float('inf');patience=0;history=[];start=time.perf_counter()
 def loss_batch(batch):
  hs=[];xs=[];es=[];eas=[];hardm=[];labels=[];indices=[];weights=[];hintidx=[];hintvals=[];offset=0
  for t in batch:
   keys,h,x,e,ea,hard=cache[t['id']];hs.append(h);xs.append(x);es.append(e+offset);eas.append(ea);hardm.append(hard)
   for k in t['hidden']:
    indices.append(keys.index(k)+offset);labels.append(t['labels']['observed_xy'][k]);weights.append(1/max(1,len(t['hidden'])))
   for k in t['input']['hints']:
    if k not in t['input']['hard_anchors']:hintidx.append(keys.index(k)+offset);hintvals.append(t['input']['hints'][k])
   offset+=len(keys)
  out=model(torch.cat(hs),torch.cat(xs),torch.cat(es,dim=1),torch.cat(eas),torch.cat(hardm));target=torch.tensor(labels,dtype=torch.float32);w=torch.tensor(weights,dtype=torch.float32);value=(((out[indices]-target)**2).mean(1)*w).sum()/len(batch)/10000
  if hintidx:value=value+.03*((out[hintidx]-torch.tensor(hintvals,dtype=torch.float32))**2).mean()/10000
  return value
 for epoch in range(60):
  model.train();order=list(training);random.shuffle(order);train_losses=[]
  for startidx in range(0,len(order),12):
   batch=order[startidx:startidx+12];optimizer.zero_grad();value=loss_batch(batch);value.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),2);optimizer.step();train_losses.append(float(value.detach()))
  model.eval()
  with torch.no_grad():v=float(loss_batch(validation))
  history.append(dict(epoch=epoch,train=float(np.mean(train_losses)),validation=v))
  if v<bestloss-1e-5:bestloss=v;best={k:x.detach().clone() for k,x in model.state_dict().items()};patience=0
  else:patience+=1
  if patience>=8:break
 model.load_state_dict(best);folder=OUT/'egnn';folder.mkdir(parents=True,exist_ok=True);torch.save(model.state_dict(),folder/(held+'-'+str(seed)+'.pt'));record=dict(held=held,seed=seed,status='trained',layers=4,hidden=64,coordinate_step_scale=32,epochs=len(history),training_maps=len({t['map_id'] for t in training}),validation_maps=len({t['map_id'] for t in validation}),test_maps=len({t['map_id'] for t in test}),train_lineages=sorted({t['lineage'] for t in training}),validation_lineages=sorted(valmaps),test_lineages=sorted({t['lineage'] for t in test}),seconds=time.perf_counter()-start,history=history,torch=torch.__version__,device='cpu',supervision='only withheld observed expert coordinates; derived centers masked')
 dump(folder/(held+'-'+str(seed)+'.json'),record);return model,record,OffsetRidge().fit(training)

def evaluate():
 tasks=read(OUT/'completion_tasks.json')['tasks'];rows=[];records=[]
 for held in sorted({t['held_source'] for t in tasks}):
  tests=[t for t in tasks if t['held_source']==held and t['split']=='test']
  for seed in SEEDS:
   model,record,ridge=train_fold(tasks,held,seed);records.append(record)
   for t in tests:
    for method in ('initial','fcose','ridge','egnn'):
     if method in ('ridge','egnn') and model is None:continue
     try:
      start=time.perf_counter();xy=completion.baseline(t,method,seed) if method in ('initial','fcose') else ridge.predict(t) if method=='ridge' else predict(model,t);metric=completion.measure(t,xy);rows.append(dict(task=t['id'],query=t['query'],map_id=t['map_id'],held_source=held,source=t['source'],kind=t['anchor_kind'],mask=t['mask'],seed=seed,method=method,status='ok',seconds=time.perf_counter()-start,**metric))
     except Exception as e:rows.append(dict(task=t['id'],seed=seed,method=method,status='failed',error=str(e)))
   print('Completion fold',held,seed,record['status'],len(tests),'tests',flush=True)
 dump(OUT/'completion_results.json',rows);dump(OUT/'egnn_training.json',records)
 return rows
if __name__=='__main__':evaluate()
