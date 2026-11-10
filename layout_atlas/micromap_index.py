"""Publish model applicability separately from immutable expert source instances."""
import gzip
import sqlite3
from core import *
import database


def publish(summary=None):
    summary=summary or read(OUT/'micromap_extension/summary.json')
    db=sqlite3.connect(DATA/'extensions/micromap/atlas.sqlite');db.execute('PRAGMA foreign_keys=ON')
    db.executescript('''CREATE TABLE IF NOT EXISTS target_models(id TEXT PRIMARY KEY,model_hash TEXT,payload BLOB);
    CREATE TABLE IF NOT EXISTS model_applicability(proposal_id TEXT,target_model TEXT,module_id TEXT,evidence TEXT,selected INTEGER,payload BLOB,PRIMARY KEY(proposal_id,target_model),FOREIGN KEY(module_id) REFERENCES modules(id));
    CREATE INDEX IF NOT EXISTS applicability_module_idx ON model_applicability(module_id);
    CREATE INDEX IF NOT EXISTS applicability_target_idx ON model_applicability(target_model);
    CREATE TABLE IF NOT EXISTS applicability_reactions(target_model TEXT,reaction_id TEXT,proposal_id TEXT,PRIMARY KEY(target_model,reaction_id,proposal_id));
    CREATE INDEX IF NOT EXISTS applicability_reaction_idx ON applicability_reactions(target_model,reaction_id);
    DELETE FROM applicability_reactions;DELETE FROM model_applicability;DELETE FROM target_models;''')
    counts={}
    for model in summary['models']:
        name=model['model'];folder=OUT/'micromap_extension/models'/name
        with gzip.open(folder/'G_D_M-selected.json.gz','rt') as f:selected={p['id'] for p in json.load(f) if p['source']=='micromap'}
        with gzip.open(folder/'micromap_candidates.json.gz','rt') as f:candidates=json.load(f)
        db.execute('INSERT INTO target_models VALUES(?,?,?)',(name,model['target_model_hash'],database.pack({k:model[k] for k in ('model','target_file','target_sha256','target_model_hash','conditions')})))
        for p in candidates:
            # Target initial coordinates are algorithmic, not new expert examples.
            payload={k:p[k] for k in ('id','template_id','map_id','lineage','reaction_ids','node_mapping','reaction_mapping','mapping_evidence','source_glyph_mapping','label_compartments_verified','license','source_citation')}
            payload['scope']='model applicability only; no target expert coordinates or new independent design'
            evidence='full_equation_source_labels_match' if p['label_compartments_verified'] else 'full_equation_compartment_role_remap'
            db.execute('INSERT INTO model_applicability VALUES(?,?,?,?,?,?)',(p['id'],name,p['template_id'],evidence,int(p['id'] in selected),database.pack(payload)))
            for rid in p['reaction_ids']:db.execute('INSERT INTO applicability_reactions VALUES(?,?,?)',(name,rid,p['id']))
        counts[name]=dict(candidates=len(candidates),selected_modules=len(selected),candidate_reactions=len({rid for p in candidates for rid in p['reaction_ids']}))
    db.commit();db.close()
    result=dict(by_target=counts,expert_master_designs=1,counts_not_new_expert_instances=True,code_sha256=sha(Path(__file__).read_bytes()))
    dump(OUT/'micromap_extension/applicability_index.json',result)
    return result

if __name__=='__main__':print(publish())
