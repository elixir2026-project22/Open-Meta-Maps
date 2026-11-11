"""Frozen-candidate P0 ablations, opt-in and separate from saved v0.3 results."""

import sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(__file__).resolve().parents[1]
if str(_ROOT / 'layout_atlas') not in _sys.path:
    _sys.path.insert(0, str(_ROOT / 'layout_atlas'))
if str(_ROOT) not in _sys.path:
    _sys.path.insert(0, str(_ROOT))

import copy
import gzip
import html
import time
from collections import Counter
from pathlib import Path
from core import *
import drawing
import geometry
import pages




def metric_details(path):
    """Actual exported glyphs/labels, excluding routing control points."""
    data = read(path); content = data[1]
    nodes = {k:v for k,v in content['nodes'].items() if v['node_type']!='multimarker'}
    xy = {k:[v['x'],v['y']] for k,v in nodes.items()}; boxes={}
    for key,v in nodes.items():
        if v['node_type']=='metabolite':
            label=v['bigg_id'];x=v['label_x'];y=v['label_y']
        else:
            rx=next(r for r in content['reactions'].values() if any(s['from_node_id']==key or s['to_node_id']==key for s in r['segments'].values()))
            label=rx['bigg_id'];x=rx['label_x'];y=rx['label_y']
        b=drawing.FONT.getbbox(label)
        boxes[key]=[x,y-(b[3]-b[1]),x+drawing.FONT.getlength(label),y]
    node_label=sum(b[0]<v[0]+4 and b[2]>v[0]-4 and b[1]<v[1]+4 and b[3]>v[1]-4
                   for b in boxes.values() for v in xy.values())
    canvas=content['canvas']; values=np.array(list(xy.values()))
    # A line enters a rectangle iff the clipped open parameter interval is nonempty.
    segments=[]
    for r in content['reactions'].values():
        segments.extend((s['from_node_id'],s['to_node_id']) for s in r['segments'].values())
    intersections=0
    for a,b in segments:
        av=np.array([content['nodes'][a]['x'],content['nodes'][a]['y']]);bv=np.array([content['nodes'][b]['x'],content['nodes'][b]['y']]);delta=bv-av
        for k,box in boxes.items():
            if k in (a,b):continue
            lower=0.;upper=1.
            for axis in (0,1):
                if abs(delta[axis])<1e-9:
                    if av[axis]<=box[axis] or av[axis]>=box[axis+2]:upper=-1;break
                else:
                    lo=(box[axis]-av[axis])/delta[axis];hi=(box[axis+2]-av[axis])/delta[axis]
                    lower=max(lower,min(lo,hi));upper=min(upper,max(lo,hi))
            intersections+=int(lower<upper-1e-8)
    return dict(visible_nodes=len(nodes),label_node_intersections=int(node_label),
                segment_label_intersections=intersections,
                canvas_width=canvas['width'],canvas_height=canvas['height'],
                canvas_area_per_reaction=canvas['width']*canvas['height']/len(content['reactions']),
                node_coordinate_extent=np.ptp(values,axis=0).tolist())


def assert_chemistry(model, path):
    data=read(path); lookup={r['bigg_id']:r for r in data[1]['reactions'].values()}
    assert set(lookup)=={r['id'] for r in model['reactions']}
    for r in model['reactions']:
        sign=-1 if r['upper_bound']<=0 and r['lower_bound']<0 else 1
        assert {m['bigg_id']:m['coefficient'] for m in lookup[r['id']]['metabolites']}=={k:sign*v for k,v in r['metabolites'].items()}
        assert data[0]['metadata']['source_reactions'][r['id']]['stoichiometry']==r['metabolites']
    assert data[0]['metadata']['model_hash']==uid(model)


def selected_candidates(folder):
    manifest=read(folder/'skeleton.json')
    with gzip.open(folder/'candidates.json.gz','rt') as f:pool={p['id']:p for p in json.load(f)}
    return manifest,[copy.deepcopy(pool[b['id']]) for b in manifest['modules']]


def run(model_name='e_coli_core'):
    import experiments
    source=OUT/'models'/model_name/'all_geometry'
    manifest,selected=selected_candidates(source)
    target=next(x[1] for x in experiments.MODELS if x[0]==model_name)
    model=load_model(target); fingerprint=uid(model);g=model_graph(model)
    if len(model['reactions'])>150:
        raise ValueError('P0 full-model ablation is bounded to <=150 reactions; use fixed subviews for large targets')
    destination=OUT/'geometry_p0'/model_name;destination.mkdir(parents=True,exist_ok=True)
    conditions=['legacy_replay','scaled_hard','display_soft','label_aware']
    results=[]; solved=None; scene_ids=None; solved_views={}; view_rows=[]
    views=read(source/'views.json')
    for condition in conditions:
        start=time.perf_counter()
        if condition=='legacy_replay':
            blocks=copy.deepcopy(manifest['modules']);anchors={k:v['xy'] for k,v in manifest['anchors'].items()}
        else:
            anchors,blocks,containers,routes=drawing.pack_modules(model,copy.deepcopy(selected),manifest['bridges'],physical_scale=True)
            for bridge in manifest['bridges']:
                a=anchors.get(bridge['nodes'][0]);b=anchors.get(bridge['nodes'][-1])
                if a is None or b is None:continue
                for i,rid in enumerate(bridge['reaction_ids']):anchors['r:'+rid]=(np.array(a)+(np.array(b)-a)*(i+1)/(len(bridge['reaction_ids'])+1)).tolist()
        core=drawing.complete(g,anchors,iterations=350)
        scene=drawing.display_scene(model,core,blocks)
        identity=(set(scene['positions']),scene['participants'],scene['copies'])
        if scene_ids is None:scene_ids=identity
        assert scene_ids==identity,'Fixed display graph must not change across ablations'
        audit={}
        if condition in ('display_soft','label_aware'):
            if solved is None:solved=geometry.complete_scene(model,scene,blocks)
            scene,audit=copy.deepcopy(solved)
        folder=destination/condition
        row=drawing.export(model,core,folder,'map',blocks=blocks,full_model=model,route=True,
                           scene=scene,label_obstacles=condition=='label_aware')
        assert row['routing']['status']=='ok',row['routing']
        row.update(metric_details(folder/'map.json'));assert_chemistry(model,folder/'map.json')
        row.update(condition=condition,seconds=time.perf_counter()-start,constraints=audit,
                   module_scales=[dict(id=b['id'],reactions=len(b['reaction_ids']),scale=b.get('physical_scale',100)) for b in blocks])
        # Persist exact pre-export scene IDs and typed evidence, not only reaction sets.
        dump(folder/'scene.json',dict(positions=scene['positions'],labels=scene['labels'],copies=scene['copies'],
             participants=[dict(reaction=r,metabolite=m,glyph=k) for (r,m),k in scene['participants'].items()],
             module_positions={b['id']:b['positions'] for b in blocks},anchors=scene.get('anchors',{}),constraints=audit))
        # Fixed target partitions use identical glyphs/visible chemistry in every condition.
        # They are dependent views of one model, not extra independent expert maps.
        for view in views:
            sub=drawing.subset(model,view['reaction_ids']);rids=set(view['reaction_ids'])
            base_scene=drawing.display_scene(model,core,blocks)
            participants={(r,m):k for (r,m),k in base_scene['participants'].items() if r in rids}
            needed={'r:'+r for r in rids}|set(participants.values())
            subscene=dict(positions={k:v for k,v in base_scene['positions'].items() if k in needed},
                          labels={k:v for k,v in base_scene['labels'].items() if k in needed},
                          copies={k:v for k,v in base_scene['copies'].items() if k in needed},participants=participants)
            subaudit={}
            if condition in ('display_soft','label_aware'):
                if view['id'] not in solved_views:solved_views[view['id']]=geometry.complete_scene(sub,subscene,blocks)
                subscene,subaudit=copy.deepcopy(solved_views[view['id']])
            viewfolder=folder/view['id']
            vr=drawing.export(sub,core,viewfolder,'map',scene=subscene,full_model=model,
                              route=False,label_obstacles=condition=='label_aware')
            vr.update(metric_details(viewfolder/'map.json'));assert_chemistry(sub,viewfolder/'map.json')
            vr.update(view=view['id'],condition=condition,constraints=subaudit)
            dump(viewfolder/'metrics.json',vr);view_rows.append(vr)
        row['seconds']=time.perf_counter()-start
        dump(folder/'metrics.json',row);results.append(row)
        print(condition,json.dumps({k:row[k] for k in ('node_overlap_pairs','label_overlap_pairs','label_node_intersections','straight_crossings','canvas_area_per_reaction','seconds')}),flush=True)
    assert uid(model)==fingerprint
    summary=dict(model=model_name,model_hash=fingerprint,target_sha256=sha(target.read_bytes()),
                 seed=SEEDS[0],source_skeleton_sha256=sha((source/'skeleton.json').read_bytes()),
                 source_candidates_sha256=sha((source/'candidates.json.gz').read_bytes()),
                 code_sha256={**{name:sha((BASE/name).read_bytes()) for name in ['geometry.py','drawing.py','routing.py']},
                              'experiments/geometry_experiment.py':sha((_ROOT/'experiments'/'geometry_experiment.py').read_bytes())},
                 display_graph_hash=uid([sorted(scene_ids[0]),sorted((r,m,k) for (r,m),k in scene_ids[1].items()),scene_ids[2]]),
                 frozen_candidate_ids=[p['id'] for p in selected],
                 config=geometry.CONFIG,results=results,fixed_view_results=view_rows,
                 interpretation='Frozen-candidate geometry ablation, not new database coverage or learning evidence; fixed glyph policy across all four conditions')
    dump(destination/'summary.json',summary)
    cols=['condition','node_overlap_pairs','label_overlap_pairs','label_node_intersections','straight_crossings','canvas_area_per_reaction']
    (destination/'index.html').write_text(pages.geometry_index(summary['results'],destination))
    from experiments import geometry_report;geometry_report.generate(summary)
    return summary

if __name__=='__main__':run()
