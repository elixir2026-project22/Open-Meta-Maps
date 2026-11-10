"""Replay the entire bounded E. coli composition search under another hash seed."""
import os
import subprocess
from core import *


def check():
    folder=OUT/'composition/e_coli_core/full';before=read(folder/'summary.json')
    old_nodes=read(folder/'best/map.json')[1]['nodes']
    env=dict(os.environ);env['PYTHONHASHSEED']='42'
    result=subprocess.run([sys.executable,str(ROOT/'project.py'),'compose','--model','e_coli_core'],cwd=ROOT,env=env,text=True,capture_output=True,timeout=240)
    if result.returncode:raise RuntimeError(result.stdout+'\n'+result.stderr)
    after=read(folder/'summary.json');nodes=read(folder/'best/map.json')[1]['nodes']
    assert before['winner']==after['winner'];assert len(before['trials'])==len(after['trials'])
    for a,b in zip(before['trials'],after['trials']):
        assert (a['id'],a['status'],a.get('acceptable'))==(b['id'],b['status'],b.get('acceptable'))
        if a['status']=='ok':assert abs(a['quality']-b['quality'])<1e-7
    assert set(old_nodes)==set(nodes)
    error=max(max(abs(nodes[k][axis]-old_nodes[k][axis]) for axis in ('x','y')) for k in nodes)
    assert error<1e-7,error
    audit=dict(model='e_coli_core',view='full',python_hash_seed='42',layout_seed=SEEDS[0],
               complete_layouts_replayed=len(after['trials']),same_winner=True,
               max_coordinate_difference=error,code_sha256=sha(Path(__file__).read_bytes()),
               interpretation='Same runtime and inputs; deterministic search replay, not independent scientific replication')
    dump(OUT/'composition/reproducibility_check.json',audit);return audit

if __name__=='__main__':print(check())
