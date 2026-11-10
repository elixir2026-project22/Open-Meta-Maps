"""Replay a fixed large-model view under a caller-selected Python hash seed."""
import os
import tempfile
import gzip
from core import *
import drawing
import micromap_experiment as experiment


def check():
    folder=OUT/'micromap_extension/models/iML1515';summary=read(folder/'summary.json')
    model=load_model(ROOT/summary['target_file'])
    views=read(OUT/'models/iML1515/all_geometry/views.json');view=views[len(views)//2]
    with gzip.open(folder/'G_D_M-selected.json.gz','rt') as f:selected=json.load(f)
    sub=drawing.subset(model,view['reaction_ids']);selected=experiment.clip_proposals(selected,sub)
    expected=read(folder/'G_D_M'/view['id']/'map.json')[1]['nodes']
    with tempfile.TemporaryDirectory() as tmp:
        experiment.draw_case(sub,selected,Path(tmp),full_model=model)
        actual=read(Path(tmp)/'map.json')[1]['nodes']
    assert set(actual)==set(expected),'Node or routing point identities differ'
    error=max(max(abs(actual[k][axis]-expected[k][axis]) for axis in ('x','y')) for k in actual)
    assert error<1e-8,error
    result=dict(model='iML1515',view=view['id'],condition='G_D_M',python_hash_seed=os.environ.get('PYTHONHASHSEED'),
                layout_seed=SEEDS[0],maximum_coordinate_difference=error,nodes_including_route_controls=len(actual),
                target_sha256=summary['target_sha256'],code_sha256=sha(Path(__file__).read_bytes()),
                interpretation='Same target, frozen choices and runtime; checks component order and drawing replay, not independent scientific replication')
    dump(OUT/'micromap_extension/reproducibility_check.json',result);return result

if __name__=='__main__':print(check())
