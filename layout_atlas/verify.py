"""Independent verification of current source inputs and persisted experiment artifacts."""
import gzip,subprocess,platform,importlib.metadata,re,zlib,hashlib
from jsonschema import Draft4Validator
from core import *
import database

def verify():
 checks=Counter();tests=subprocess.run([sys.executable,'-m','unittest','discover','-s',str(BASE),'-p','test_v03.py','-v'],capture_output=True,text=True);(OUT/'test_run.log').write_text(tests.stdout+tests.stderr);assert tests.returncode==0,tests.stderr;checks['unit_tests']=int(re.search(r'Ran (\d+) tests',tests.stderr)[1])
 for path,digest in freeze()['files'].items():
  assert sha((INPUTS/path).read_bytes())==digest,path;checks['source_input_files_hashed']+=1
 for meta in DATA.rglob('*.source.json'):
  record=read(meta);path=Path(str(meta)[:-len('.source.json')]);assert sha(path.read_bytes())==record['sha256'],str(path);checks['cached_sources_sha256']+=1
 acq=read(DATA/'acquisition.json')
 for entry in acq['reactome_diagrams']:
  assert entry['archive_record']==21383214
  for member in entry['members']:
   raw=(BASE/member['path']).read_bytes();assert sha(raw)==member['sha256'];assert zlib.crc32(raw)&0xffffffff==member['zip_crc32'];checks['fixed_archive_members']+=1
 mods=database.modules();lookup={m['id']:m for m in mods};db=database.connect();mapids={r[0] for r in db.execute('SELECT id FROM maps')};db.close()
 for m in mods:
  assert m['map_id'] in mapids;keys={n['key'] for n in m['nodes']};assert all(e['a'] in keys and e['b'] in keys for e in m['edges']);assert m['source_sha256'];assert all(n['observed']==(n['position_origin']=='observed') for n in m['nodes']);checks['module_records']+=1
 with gzip.open(OUT/'ranking_queries.json.gz','rt') as f:ranking=json.load(f)
 for q in ranking:
  for c in q['candidates']:
   assert c['source']!=q['held_source'] and c['lineage']!=q['lineage'];assert c['template_id'] in lookup;checks['isolated_ranking_candidates']+=1
 for row in read(OUT/'ranking_summary.json'):
  if row['status']=='trained':assert row['training_preference_queries']>=50 and row['test_multivariant']>=20
 tasks=read(OUT/'completion_tasks.json')['tasks'];pairs=defaultdict(list)
 for t in tasks:
  pairs[t['held_source'],t['query'],t['mask']].append(t);assert not(set(t['derived_unlabelled'])&set(t['hidden']));assert set(t['input']['hints'])==set(t['visible']);assert t['retrieved_source']!=t['held_source'] and t['retrieved_lineage']!=t['lineage'];assert not any('xy' in n for n in t['input']['nodes']);assert set(t['hidden']).isdisjoint(t['visible']);checks['typed_completion_tasks']+=1
 for ts in pairs.values():
  assert {t['anchor_kind'] for t in ts}=={'expert','perturbed','retrieved'};assert len({tuple(t['hidden']) for t in ts})==len({tuple(t['visible']) for t in ts})==1;checks['identical_mask_triplets']+=1
 for r in read(OUT/'egnn_training.json'):
  assert r['status']=='trained';assert not(set(r['train_lineages'])&set(r['validation_lineages']));assert not((set(r['train_lineages'])|set(r['validation_lineages']))&set(r['test_lineages']));assert r['layers']==4 and r['hidden']==64;checks['isolated_egnn_runs']+=1
 for r in read(OUT/'completion_results.json'):
  assert r['status']=='ok',r;assert r['hard_anchor_drift']<1e-5;checks['completion_evaluations']+=1
 validator=Draft4Validator(read(REF/'vendor/escher/map-schema-1-0-0.json'))
 def check_export(path,ids,target):
  data=read(path);validator.validate(data);expected={r['id']:r for r in target['reactions']};actual={r['bigg_id']:r for r in data[1]['reactions'].values()};assert set(actual)==set(ids),(str(path),'scope')
  nodes=data[1]['nodes']
  for rid,r in actual.items():
   original=expected[rid];sign=-1 if original['upper_bound']<=0 and original['lower_bound']<0 else 1
   assert {m['bigg_id']:m['coefficient']*sign for m in r['metabolites']}==original['metabolites'],(str(path),rid,'stoich')
   meta=data[0]['metadata']['source_reactions'][rid];assert meta['stoichiometry']==original['metabolites'];assert meta['lower_bound']==original['lower_bound'] and meta['upper_bound']==original['upper_bound']
   network=nx.Graph([(s['from_node_id'],s['to_node_id']) for s in r['segments'].values()]);seen_mets={nodes[k]['bigg_id'] for k in network if nodes[k]['node_type']=='metabolite'};assert seen_mets==set(original['metabolites']);assert nx.is_connected(network);checks['chemically_intact_exported_reactions']+=1
  checks['schema_valid_views']+=1
 model_rows=read(OUT/'model_summary.json')+[read(p) for p in (OUT/'models').glob('*/learned*/summary.json')]
 for row in model_rows:
  assert row['status']=='ok',row;folder=OUT/'models'/row['model']/row['condition'];state=read(folder/'skeleton.json');target=load_model(ROOT/state['target_file']);assert sha((ROOT/state['target_file']).read_bytes())==state['target_sha256'];assert uid(target)==state['target_model_hash'];rids={r['id'] for r in target['reactions']};assert set(state['reaction_status'])==rids
  uses=defaultdict(set)
  for r in target['reactions']:
   for mid in r['metabolites']:uses[mid].add(r['id'])
  covered=[];target_graph=model_graph(target)
  for b in state['modules']:
   assert len(b['reaction_ids'])>=2;covered.extend(b['reaction_ids']);assert b['template_id'] in lookup
   source=lookup[b['template_id']];mapping=b['node_mapping']
   for n in source['nodes']:assert target_graph.nodes[mapping[n['key']]]['role']==n['role']
   for e in source['edges']:assert target_graph.has_edge(mapping[e['a']],mapping[e['b']])
   for p in b['target_ports']:assert p['outside_reactions']==sorted(uses[p['metabolite']]-set(b['reaction_ids']));checks['true_target_ports']+=1
  assert row['strict_equation_reactions']==len({rid for b in state['modules'] if b['complete_equation_verified'] for rid in b['reaction_ids']})
  assert len(covered)==len(set(covered));assert set(covered)=={k for k,v in state['reaction_status'].items() if v=='expert_module'}
  views=read(folder/'views.json');ids=[r for v in views for r in v['reaction_ids']];assert len(ids)==len(set(ids)) and set(ids)==rids
  for v in views:check_export(folder/v['file'],v['reaction_ids'],target)
  for v in row['full_small_comparisons']:check_export(folder/v['file'],sorted(rids),target)
  dep=folder/'egnn_deployment.json'
  if dep.exists():
   for v in read(dep)['views']:
    if v['status']=='ok':check_export(folder/v['file'],next(x['reaction_ids'] for x in views if x['id']==v['view']),target)
  checks['complete_target_accounting_runs']+=1
 # New opt-in P0 artifacts share the exact same target and glyph policy.
 for summary_path in (OUT/'geometry_p0').glob('*/summary.json'):
  folder=summary_path.parent;summary=read(summary_path);source=OUT/'models'/summary['model']/'all_geometry';state=read(source/'skeleton.json');target=load_model(ROOT/state['target_file']);rids={r['id'] for r in target['reactions']}
  assert summary['target_sha256']==sha((ROOT/state['target_file']).read_bytes());assert summary['model_hash']==uid(target)
  assert summary['source_skeleton_sha256']==sha((source/'skeleton.json').read_bytes());assert summary['source_candidates_sha256']==sha((source/'candidates.json.gz').read_bytes())
  for name,digest in summary['code_sha256'].items():assert sha((BASE/name).read_bytes())==digest,('stale P0 result',name)
  assert summary['frozen_candidate_ids']==[b['id'] for b in state['modules']]
  glyph_graphs=[];scopes={v['id']:v['reaction_ids'] for v in read(source/'views.json')}
  for row in summary['results']:
   dest=folder/row['condition'];scene=read(dest/'scene.json');glyph_graphs.append(uid([sorted(scene['positions']),scene['participants'],scene['copies']]))
   check_export(dest/'map.json',rids,target);assert row['routing']['status']=='ok'
   for key,hint in scene['anchors'].items():
    delta=np.abs(np.array(scene['positions'][key])-hint['xy']);assert np.max(delta)<1e-7 if hint['hard'] else np.max(delta)<=hint['allowed_motion']+1e-7
    checks['p0_anchor_contracts']+=1
   for view,ids in scopes.items():check_export(dest/view/'map.json',ids,target)
   checks['p0_ablation_conditions']+=1
  assert len(set(glyph_graphs))==1;checks['p0_frozen_display_graphs']+=1
 if (OUT/'micromap_extension/summary.json').exists():
  import micromap_verify;checks.update(micromap_verify.verify(check_export,lookup))
 if (OUT/'composition/summary.json').exists():
  import composition_verify;checks.update(composition_verify.verify(check_export))
 source_db=database.connect();published=sqlite3.connect(DATA/'open_core.sqlite');source_db.backup(published)
 used_maps={m['map_id'] for m in mods}
 for (mid,) in published.execute('SELECT id FROM maps').fetchall():
  if mid not in used_maps:published.execute('DELETE FROM maps WHERE id=?',(mid,))
 published.commit();published.close();source_db.close()
 packages={p:importlib.metadata.version(p) for p in ('networkx','numpy','scipy','scikit-learn','Pillow','jsonschema','torch')};result=dict(status='passed',checks=dict(checks),python=platform.python_version(),label_font=dict(path=str(__import__('drawing').FONT_PATH),size=12),packages=packages,seeds=SEEDS,code_sha256={p.name:sha(p.read_bytes()) for p in BASE.glob('*.py')},database_sha256=sha((DATA/'atlas.sqlite').read_bytes()),review_policy='No required human or dual-review gate',notes=['Automated checks are not human certification.','Schema and chemistry validity do not establish readability.','Whole Reactome archive checksum not independently verified; selected members checked.']);dump(OUT/'verification.json',result);print(json.dumps(result['checks'],indent=2));return result
if __name__=='__main__':verify()
