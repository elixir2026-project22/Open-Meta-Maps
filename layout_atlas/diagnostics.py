"""Matched-corpus family comparison and explicit retained-but-unmatched observations."""
import sqlite3
from core import *
import database

def run():
 oldrows={r['id']:(r['family'],r['variant'],json.dumps(r['observation'])) for r in source_module_records()};db=database.connect();new={r[0]:(r[1],r[2]) for r in db.execute('SELECT id,core_family,variant FROM modules') if r[0] in oldrows};shared=set(new);missing=set(oldrows)-shared;reasons=Counter()
 db.execute('CREATE TABLE IF NOT EXISTS retained_observations(id TEXT PRIMARY KEY,reason TEXT,source_pointer TEXT)')
 for key in sorted(missing):
  m=json.loads(oldrows[key][2]);reason='reactome_replaced_by_fixed_v97' if m['source']=='reactome' else 'outside_v03_size_bound' if len(m.get('core_nodes',m['nodes']))>100 or (m['size'] if m['quality']=='semantic_transfer' else len(m['edges']))>40 else 'disconnected_visible_core';reasons[reason]+=1;db.execute('INSERT OR REPLACE INTO retained_observations VALUES(?,?,?)',(key,reason,'data/source_inputs/module_observations.json.gz#'+key))
 db.commit();summary=dict(legacy_total=len(oldrows),same_instances_compared=len(shared),legacy_full_display_families=len({oldrows[k][0] for k in shared}),v03_core_shape_families=len({new[k][0] for k in shared}),legacy_layout_variants=len({oldrows[k][1] for k in shared}),v03_core_variants=len({new[k][1] for k in shared}),retained_not_in_active_matcher=len(missing),reasons=dict(reasons),note='Same retained IDs compared; core role abstraction differs, so fewer families does not prove expert agreement or readability.')
 dump(OUT/'matched_corpus_comparison.json',summary);db.close();return summary
if __name__=='__main__':print(run())
