"""Independent checks of composition choices, immutable parentage and guards."""
import gzip
from core import *
import composition as c
import micromap_experiment as pipeline


def gzread(path):
    with gzip.open(path,'rt') as f:return json.load(f)


def verify(check_export):
    root=OUT/'composition';summary=read(root/'summary.json');checks=Counter()
    assert summary['config']==c.CONFIG
    assert summary['p1_summary_sha256']==sha((OUT/'micromap_extension/summary.json').read_bytes())
    for name,digest in summary['code_sha256'].items():assert sha((BASE/name).read_bytes())==digest,('stale composition',name)
    for path,digest in summary['candidate_input_sha256'].items():assert sha((ROOT/path).read_bytes())==digest,('changed composition input',path)
    p1={m['model']:m for m in read(OUT/'micromap_extension/summary.json')['models']};last_name=None;pool={}
    for case in summary['cases']:
        name=case['model'];record=p1[name];target=load_model(ROOT/record['target_file'])
        sub=drawing_subset(target,case['display_reaction_ids']);assert uid(sub)==case['target_model_hash']
        assert uid(target)==case['full_target_model_hash']
        if name!=last_name:
            folder=OUT/'micromap_extension/models'/name
            candidates=gzread(OUT/'models'/name/'all_geometry/candidates.json.gz')+gzread(folder/'direct_core_candidates.json.gz')+gzread(folder/'micromap_candidates.json.gz')
            pool={p['id']:p for p in candidates};last_name=name
        folder=root/name/case['view'];required=set(case['required_reactions']);chemical=set(case['required_chemical_reactions'])
        assert case['full_layout_attempts']<=c.CONFIG['full_layouts'];assert case['macro_previews']<=c.CONFIG['macro_previews']
        assert case['raw_candidates']<=c.CONFIG['raw_budget']
        eligible=[]
        for trial in case['trials']:
            dest=folder/trial['folder'];assert (dest/'selected.json.gz').exists()
            if trial['status']=='failed':assert trial['error'];checks['composition_recorded_failures']+=1;continue
            selected=gzread(dest/'selected.json.gz');assert c.valid_selection(selected,required,chemical)
            for proposal in selected:
                clipping=proposal.get('view_clipping')
                parent_id=clipping['original_proposal'] if clipping else proposal.get('composition_parent_id',proposal['id'])
                parent=pool[parent_id]
                if clipping:
                    cropped=pipeline.clip_proposals([parent],sub);assert len(cropped)==1;parent=cropped[0]
                assert parent['id']==proposal.get('composition_parent_id',proposal['id'])
                for key in ('reaction_ids','reaction_mapping','node_mapping'):assert proposal[key]==parent[key]
                assert c.chemical(proposal)==c.chemical(parent)
                assert set(proposal['positions'])==set(parent['positions'])
                if proposal.get('algorithmic_pose'):
                    pose=proposal['algorithmic_pose'];assert pose['degrees'] in c.CONFIG['rotations'];assert pose['source_instance_unchanged']
                    center=np.array(pose['center']);matrix=np.array(pose['matrix']);assert np.allclose(matrix.T@matrix,np.eye(2),atol=1e-8)
                    assert np.linalg.det(matrix)>.999999
                    for k,xy in parent['positions'].items():assert np.allclose((np.array(xy)-center)@matrix+center,proposal['positions'][k],atol=1e-8)
                    checks['composition_recorded_rotations']+=1
                else:
                    for k,xy in parent['positions'].items():assert np.allclose(xy,proposal['positions'][k],atol=1e-8)
                checks['composition_source_parentages']+=1
            check_export(dest/'map.json',case['display_reaction_ids'],target)
            state=read(dest/'skeleton.json');assert set(state['reaction_status'])==set(case['display_reaction_ids'])
            assert set(state['display_scope_reaction_ids'])==set(case['display_reaction_ids'])
            for k,hint in state['anchors'].items():
                delta=np.abs(np.array(state['display_positions'][k])-hint['xy'])
                assert np.max(delta)<1e-7 if hint['hard'] else np.max(delta)<=hint['allowed_motion']+1e-7
            metric=trial['metric'];assert abs(c.quality(metric)-trial['quality'])<1e-7
            allowed=True if trial['kind']=='seed' else c.acceptable(metric,case['baseline'])
            assert allowed==trial['acceptable']
            if allowed:eligible.append(trial['quality'])
            checks['composition_checked_layouts']+=1
        assert c.acceptable(case['best'],case['baseline'])
        assert case['best_quality']<=min(eligible)+1e-7
        chosen=gzread(folder/'best/selected.json.gz');covered,cc=c.coverage(chosen)
        assert len(covered)==case['after_covered'] and len(cc)==case['after_chemical']
        assert c.valid_selection(chosen,required,chemical)
        check_export(folder/'best/map.json',case['display_reaction_ids'],target)
        assert read(folder/'best/metrics.json')==case['best']
        checks['composition_fixed_cases']+=1
    return checks


def drawing_subset(model,ids):
    import drawing
    return drawing.subset(model,ids)
