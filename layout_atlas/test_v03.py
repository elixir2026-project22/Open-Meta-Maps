"""Behavioral tests for the v0.3 research contracts."""
import unittest,tempfile,copy
from pathlib import Path
from core import *
import semantics,database,retrieval,drawing,completion

class ContractTests(unittest.TestCase):
 def test_rhea_directions_and_cardinality(self):
  r=semantics.rhea_index();self.assertEqual(r['master']['10001'],'10000');self.assertTrue(r['reactions']['10000']['complete']);self.assertTrue(any(p['coefficient'] is not None and abs(p['coefficient'])>1 for p in r['reactions']['21040']['participants']))
 def test_ph_mapping_preserves_origin(self):
  r=semantics.rhea_index();item=next(iter(r['ph73'].values()));self.assertIn(item['origin'],('curation','computation'))
 def fixture(self,dashed=False):
  return f'''<Pathway xmlns="http://pathvisio.org/GPML/2013a" Name="fixture" Organism="test"><DataNode GraphId="a" Type="Metabolite"><Graphics CenterX="0" CenterY="0"/><Xref Database="ChEBI" ID="CHEBI:15377"/></DataNode><DataNode GraphId="b" Type="Metabolite"><Graphics CenterX="100" CenterY="0"/><Xref Database="ChEBI" ID="CHEBI:15378"/></DataNode><Interaction GraphId="e"><Graphics LineStyle="{'Dashed' if dashed else 'Solid'}"><Point X="0" Y="0" GraphRef="a"/><Point X="100" Y="0" GraphRef="b" ArrowHead="Arrow"/><Anchor GraphId="enzyme_anchor" Position="0.5"/></Graphics><Xref Database="Rhea" ID="10000"/></Interaction></Pathway>'''
 def parse_fixture(self,dashed=False):
  # Put the fixture below INPUTS so provenance path normalization is realistic.
  with tempfile.TemporaryDirectory(dir=INPUTS) as tmp:
   p=Path(tmp)/'WPTEST.gpml';p.write_text(self.fixture(dashed));return semantics.parse_gpml(p,semantics.rhea_index())
 def test_gpml_retains_interaction_xref_and_anchor(self):
  d=self.parse_fixture();self.assertEqual(d['edges'][0]['xref']['Database'],'Rhea');self.assertEqual(d['anchors']['enzyme_anchor']['interaction'],'e');self.assertEqual(d['edges'][0]['arrows'][1],'Arrow')
 def test_dashed_not_promoted(self):self.assertEqual(self.parse_fixture(True)['bindings'][0]['status'],'aggregate_not_promoted')
 def test_unknown_bindings_not_invented(self):
  d=self.parse_fixture();self.assertFalse(semantics.get_chebi({'Database':'HMDB','ID':'HMDB00001'}));self.assertNotEqual(d['bindings'][0]['status'],'validated_visible_participants')
 def test_shape_family_ignores_side_attachment(self):
  m=dict(nodes=[dict(key='r:x',role='reaction'),dict(key='m:a',role='metabolite')],edges=[dict(a='r:x',b='m:a')]);n=copy.deepcopy(m);n['attachments']=['atp','h'];self.assertTrue(database.mappings(database.graph(m),database.graph(n)))
 def test_reflection_alignment_preserves_roles(self):
  a={'a':[0,0],'b':[1,0],'c':[1,2]};b={k:[10-3*v[1],20+3*v[0]] for k,v in a.items()};f=fit_frame(list(a.values()),list(b.values()));got=transform(a,f);self.assertLess(max(np.linalg.norm(np.array(got[k])-b[k]) for k in a),1e-7)
 def test_graph_uses_real_incidence(self):
  model=dict(metabolites=[dict(id=x) for x in ('A_c','B_c','C_c','D_c')],reactions=[dict(id='z',lower_bound=0,upper_bound=1,metabolites={'A_c':-1,'B_c':1}),dict(id='a',lower_bound=0,upper_bound=1,metabolites={'B_c':-1,'C_c':1}),dict(id='q',lower_bound=0,upper_bound=1,metabolites={'B_c':-1,'D_c':1})]);g=model_graph(model);self.assertEqual(set(g['m:B_c']),{'r:z','r:a','r:q'});self.assertFalse(g.has_edge('r:z','r:a'))
 def test_structural_mapping_does_not_need_names(self):
  a=nx.path_graph(5);b=nx.relabel_nodes(a,{i:'new'+str(i) for i in a});nx.set_node_attributes(a,{i:'reaction' if i%2 else 'metabolite' for i in a},'role');nx.set_node_attributes(b,{str('new'+str(i)):'reaction' if i%2 else 'metabolite' for i in a},'role');self.assertTrue(database.mappings(a,b))
 def test_candidate_evidence_retains_ambiguous_aliases(self):
  refs={'a':{'chebi_ph73':['1','2']},'b':{'chebi_ph73':['1']}};self.assertEqual(set(retrieval.binding_candidates('chebi:1@local',['a','b'],refs)),{'a','b'})
 def test_cover_preserves_internal_exclusion(self):
  ps=[dict(id=k,reaction_ids=ids,features=[],score=sc) for k,ids,sc in [('a',['1','2'],4),('b',['2','3'],5),('c',['4','5'],4)]];sel=retrieval.select(ps);ids=[r for p in sel for r in p['reaction_ids']];self.assertEqual(len(ids),len(set(ids)))
 def test_anchor_frame_only_sees_provided_points(self):
  a={'a':[10,20],'b':[100,20],'c':[40,80]};f=completion.coordinate_frame(a);local=completion.to_local(a,f);shift={k:(np.array(v)+[300,-90]).tolist() for k,v in a.items()};self.assertTrue(all(np.allclose(local[k],completion.to_local(shift,completion.coordinate_frame(shift))[k]) for k in a))
 def test_soft_hard_constraints(self):
  g=nx.path_graph(['a','b','c']);xy={'a':[0,0],'b':[50,0],'c':[100,0]};hints={'a':[10,10],'c':[120,0]};out=completion.constrain(g,xy,hints,{'a':[10,10]},{'a':1,'c':.3});self.assertEqual(out['a'],[10.,10.]);self.assertNotEqual(out['c'],hints['c'])
 def test_backend_uses_preset(self):
  out=drawing.backend(dict(backend='fcose',elements=[dict(data=dict(id='x'),position=dict(x=25,y=45))],debug_initial=True));self.assertEqual(out['x'],{'x':25,'y':45})
 def test_egnn_rigid_equivariance(self):
  import torch,neural
  torch.manual_seed(1);m=neural.EGNN().eval();h=torch.randn(4,10);x=torch.tensor([[0.,0.],[50.,20.],[90.,60.],[50.,80.]]);e=torch.tensor([[0,1,1,2,2,3],[1,0,2,1,3,2]]);ea=torch.randn(6,3);hard=torch.tensor([True,False,False,False]);rot=torch.tensor([[0.,-1.],[1.,0.]])
  with torch.no_grad():a=m(h,x,e,ea,hard);b=m(h,x@rot+100,e,ea,hard)
  self.assertTrue(torch.allclose(a@rot+100,b,atol=2e-4));self.assertTrue(torch.equal(a[0],x[0]))
 def test_export_chemistry_copies_and_reverse(self):
  model=dict(id='toy',metabolites=[dict(id=k,compartment='c') for k in ['A_c','B_c','atp_c']],reactions=[dict(id='r',name='r',lower_bound=-10,upper_bound=0,metabolites={'A_c':-2,'B_c':1,'atp_c':-1}),dict(id='blocked',name='b',lower_bound=0,upper_bound=0,metabolites={'A_c':-1,'B_c':1})]);g=model_graph(model);xy=drawing.initial(g,{});original=copy.deepcopy(model)
  with tempfile.TemporaryDirectory() as tmp:
   path=Path(tmp);drawing.export(model,xy,path,'out');data=read(path/'out.json')
  self.assertEqual(model,original);self.assertEqual(data[0]['metadata']['source_reactions']['r']['stoichiometry']['A_c'],-2);self.assertEqual(len(data[1]['reactions']),2)
 def test_archive_is_fixed_official_record(self):
  r=read(OUT/'reactome_archive.json');self.assertEqual(r['record'],21383214);self.assertEqual(r['version'],'v97');self.assertEqual(r['status'],'range_indexed')
 def test_no_target_model_is_expert_source(self):
  db=database.connect();self.assertNotIn('aracore',{r[0].lower() for r in db.execute('SELECT DISTINCT source FROM modules')});db.close()
 def test_paired_masks_and_no_derived_labels(self):
  tasks=read(OUT/'completion_tasks.json')['tasks'];groups=defaultdict(list)
  for t in tasks:groups[t['held_source'],t['query'],t['mask']].append(t);self.assertFalse(set(t['hidden'])&set(t['derived_unlabelled']))
  self.assertTrue(groups)
  for ts in groups.values():self.assertEqual(len({tuple(t['hidden']) for t in ts}),1);self.assertEqual(len({tuple(t['visible']) for t in ts}),1)
 def test_wp_semantic_geometry_same_source_same_lineage(self):
  byhash=defaultdict(set)
  for m in database.modules():
   if m['source']=='wikipathways':byhash[m['source_sha256']].add(m['lineage'])
  self.assertTrue(all(len(v)==1 for v in byhash.values()))
 def test_explicit_region_does_not_move_hard_anchor(self):
  g=nx.path_graph(['a','b']);nx.set_node_attributes(g,'c','compartment');out=completion.constrain(g,{'a':[100,100],'b':[200,200]},{},{'a':[100,100]},{},regions={'c':[0,0,50,50]});self.assertEqual(out['a'],[100.,100.]);self.assertTrue(all(0<=x<=50 for x in out['b']))
 def test_prediction_never_needs_hidden_labels(self):
  import neural
  t=copy.deepcopy(read(OUT/'completion_tasks.json')['tasks'][0]);t.pop('labels');out=neural.predict(neural.EGNN().eval(),t);self.assertEqual(set(out),{n['key'] for n in t['input']['nodes']})
 def test_pinned_libavoid_orthogonal_route(self):
  import routing
  paths=routing.avoid({'a':[0,0],'b':[100,100],'obstacle':[50,50]},[('a','b')]);self.assertGreaterEqual(len(paths[0]),2)
  for a,b in zip(paths[0],paths[0][1:]):self.assertTrue(abs(a[0]-b[0])<1e-5 or abs(a[1]-b[1])<1e-5)
 def test_organized_project_roots_and_runtime(self):
  self.assertEqual(BASE,ROOT/'layout_atlas');self.assertTrue((INPUTS/'module_observations.json.gz').exists());self.assertTrue((BASE/'runtime/layout_backend.cjs').exists());self.assertTrue((REF/'model_input.py').exists())
 def test_physical_scale_responds_to_dense_modules(self):
  import geometry
  model=dict(reactions=[dict(id='r',metabolites={'A_c':-1,'B_c':1})])
  sparse=dict(reaction_ids=['r'],positions={'r:r':[0,0],'m:A_c':[-1,0],'m:B_c':[1,0]});dense=copy.deepcopy(sparse);dense['positions']={k:(np.array(v)*.01).tolist() for k,v in sparse['positions'].items()}
  self.assertGreater(geometry.module_scale(model,dense),geometry.module_scale(model,sparse));self.assertEqual(model['reactions'][0]['metabolites']['A_c'],-1)
 def test_typed_geometry_hard_and_soft_contract(self):
  import geometry
  g=nx.path_graph(['hard','soft','free']);positions={'hard':[0,0],'soft':[0,0],'free':[5,0]};hints={'hard':dict(xy=[0,0],hard=True),'soft':dict(xy=[0,0],hard=False,confidence=.5,allowed_motion=30)}
  out,audit=geometry.relax(g,positions,hints,iterations=60);self.assertEqual(out['hard'],[0.,0.]);self.assertGreater(np.linalg.norm(out['soft']),10);self.assertTrue(all(abs(v)<=30+1e-7 for v in out['soft']));self.assertEqual(audit['hard_anchor_max_drift'],0)
 def test_scene_export_consumes_solved_copies(self):
  import geometry
  model=dict(id='copies',metabolites=[dict(id=k,compartment='c') for k in ('A_c','B_c','C_c')],reactions=[dict(id='one',lower_bound=0,upper_bound=10,metabolites={'A_c':-1,'B_c':1}),dict(id='two',lower_bound=0,upper_bound=10,metabolites={'B_c':-1,'C_c':1})])
  core={'r:one':[0,0],'r:two':[100,0],'m:A_c':[-40,0],'m:B_c':[50,0],'m:C_c':[140,0]};blocks=[dict(id='first',reaction_ids=['one'],positions={'m:B_c':[30,0]})];scene=drawing.display_scene(model,core,blocks);key=scene['participants']['one','B_c'];self.assertNotEqual(key,scene['participants']['two','B_c']);scene['positions'][key]=[30,45]
  self.assertEqual(len(geometry.scene_graph(model,scene).edges()),4)
  with tempfile.TemporaryDirectory() as tmp:
   drawing.export(model,core,Path(tmp),'scene',blocks=blocks,scene=scene);data=read(Path(tmp)/'scene.json')
  positions=[v for v in data[1]['nodes'].values() if v.get('bigg_id')=='B_c'];self.assertEqual(len(positions),2);self.assertAlmostEqual(abs(positions[0]['y']-positions[1]['y']),45)
 def test_incompatible_hard_anchors_are_reported(self):
  import geometry
  g=nx.Graph();g.add_nodes_from(['a','b']);out,audit=geometry.relax(g,{'a':[0,0],'b':[0,0]},{k:dict(xy=[0,0],hard=True) for k in g},iterations=5)
  self.assertEqual(out,{'a':[0.,0.],'b':[0.,0.]});self.assertGreater(audit['remaining_spacing_violations'],0)
 def test_label_router_avoids_rectangular_obstacle(self):
  import routing
  box=[40,-20,70,20];paths=routing.avoid({'a':[0,0],'b':[110,0]},[('a','b')],obstacles=[box]);self.assertGreaterEqual(len(paths[0]),4)
  for a,b in zip(paths[0],paths[0][1:]):
   self.assertTrue(abs(a[0]-b[0])<1e-5 or abs(a[1]-b[1])<1e-5)
   self.assertFalse(abs(a[1]-b[1])<1e-5 and -20<a[1]<20 and max(a[0],b[0])>40 and min(a[0],b[0])<70)
 def test_celldesigner_local_control_frame(self):
  import celldesigner
  self.assertTrue(np.allclose(celldesigner.local_to_global([10,20],[10,120],[[.5,.2]])[0],[-10,70]))
  with self.assertRaises(ValueError):celldesigner.local_to_global([0,0],[0,0],[[.5,0]])
 def test_celldesigner_process_segment_and_unobserved_center(self):
  import celldesigner
  xml='<sbml xmlns="http://www.sbml.org/sbml/level2/version4" xmlns:cd="http://www.sbml.org/2001/ns/celldesigner"><model><annotation><cd:extension><cd:listOfSpeciesAliases><cd:speciesAlias id="a" species="s1"><cd:bounds x="0" y="0" w="20" h="20"/></cd:speciesAlias><cd:speciesAlias id="b" species="s2"><cd:bounds x="100" y="0" w="20" h="20"/></cd:speciesAlias></cd:listOfSpeciesAliases></cd:extension></annotation><listOfSpecies><species id="s1" name="a[c]" compartment="default"/><species id="s2" name="b[c]" compartment="default"/></listOfSpecies><listOfReactions><reaction id="r" name="R"><annotation><cd:extension><cd:baseReactants><cd:baseReactant species="s1" alias="a"/></cd:baseReactants><cd:baseProducts><cd:baseProduct species="s2" alias="b"/></cd:baseProducts><cd:connectScheme connectPolicy="direct" rectangleIndex="1"/><cd:editPoints>0.25,0.1 0.75,0.1</cd:editPoints></cd:extension></annotation><listOfReactants><speciesReference species="s1"/></listOfReactants><listOfProducts><speciesReference species="s2"/></listOfProducts></reaction></listOfReactions></model></sbml>'
  with tempfile.TemporaryDirectory(dir=DATA) as tmp:
   p=Path(tmp)/'fixture.xml';p.write_text(xml);d=celldesigner.decode(p)
  self.assertEqual(len(d['events']),1);r=d['events'][0];self.assertTrue(np.allclose(r['xy'],[60,20]));self.assertNotEqual(r['xy'],[60,10]);self.assertIn('not_observed',r['position_origin']);self.assertEqual(len(r['path_centerline']),4)
 def test_equation_matching_reverse_scale_and_role_constraints(self):
  import semantic_transfer as st
  result=st.unique_equation_mapping({'a@c':-2,'b@c':1},{'metabolites':{'x':4,'y':-2}},{'a@c':['x'],'b@c':['y']},{'x':'m','y':'m'})
  self.assertEqual(len(result),1);self.assertEqual(result[0]['factor'],-2);self.assertFalse(result[0]['label_compartments_exact']);self.assertEqual(result[0]['compartment_roles'],{'c':'m'})
  bad=st.unique_equation_mapping({'a@c':-1,'a@e':1},{'metabolites':{'x':-1,'y':1}},{'a@c':['x'],'a@e':['y']},{'x':'c','y':'c'});self.assertEqual(bad,[])
 def test_equation_matching_preserves_ambiguity(self):
  import semantic_transfer as st
  result=st.unique_equation_mapping({'a@c':-1,'b@c':-1},{'metabolites':{'x':-1,'y':-1}},{'a@c':['x','y'],'b@c':['x','y']},{'x':'c','y':'c'});self.assertEqual(len(result),2)
 def test_spelling_does_not_erase_stereochemistry(self):
  import semantic_transfer as st
  self.assertEqual(st.spelling('mal__L'),'mal_L');self.assertNotEqual(st.spelling('mal__L'),st.spelling('mal__D'));self.assertNotEqual(st.spelling('ATP'),st.spelling('atp'))
 def test_annotations_map_without_reaction_name_identity(self):
  import semantic_transfer as st
  qualifier='{http://biomodels.net/biology-qualifiers/}is'
  model=dict(id='foreign_namespace',metabolites=[dict(id=k,compartment='c',annotations=[dict(qualifier=qualifier,uri='https://identifiers.org/bigg.metabolite/'+a)]) for k,a in [('M1','aa'),('M2','bb')]],reactions=[dict(id='DIFFERENT',metabolites={'M1':-1,'M2':1})]);mapper=st.Mapper(model,dict(links={}))
  result=mapper.match_reaction(dict(id='SOURCE',coefficients={'aa@c':-1,'bb@c':1}));self.assertEqual(result[0]['target_reaction'],'DIFFERENT');self.assertEqual(result[0]['metabolites'],{'aa@c':'M1','bb@c':'M2'})
 def test_micromap_partition_not_in_open_core(self):
  self.assertNotIn('micromap',{m['source'] for m in database.modules()})
  path=OUT/'micromap_extension/library_summary.json'
  if path.exists():self.assertFalse(read(path)['admitted_to_open_core']);self.assertEqual(read(path)['independent_designs'],1)
 def test_repeated_source_glyphs_are_not_silently_merged(self):
  import semantic_transfer as st
  qualifier='{http://biomodels.net/biology-qualifiers/}is'
  model=dict(id='toy',metabolites=[dict(id=k,compartment='c',annotations=[dict(qualifier=qualifier,uri='https://identifiers.org/bigg.metabolite/'+a)]) for k,a in [('M1','aa'),('M2','bb')]],reactions=[dict(id='ONE',lower_bound=0,upper_bound=10,metabolites={'M1':-1,'M2':1}),dict(id='TWO',lower_bound=0,upper_bound=10,metabolites={'M2':-1,'M1':1})])
  module=dict(id='source',size=2,core_family='family',lineage='one',reactions=[dict(key='r:s1',coefficients={'aa@c':-1,'bb@c':1}),dict(key='r:s2',coefficients={'bb@c':-1,'aa@c':1})],nodes=[dict(key=k,role=role,entity=entity) for k,role,entity in [('r:s1','reaction',None),('r:s2','reaction',None),('a1','metabolite','aa@c'),('a2','metabolite','aa@c'),('b','metabolite','bb@c')]],edges=[dict(a=a,b=b) for a,b in [('r:s1','a1'),('r:s1','b'),('r:s2','b'),('r:s2','a2')]])
  out,audit=st.Mapper(model,dict(links={})).proposals([module],'test');self.assertEqual(out,[]);self.assertGreater(audit['skips']['repeated_core_copy_requires_graphic_matcher'],0)
 def test_extension_fallback_graph_has_canonical_component_order(self):
  from micromap_experiment import CanonicalGraph
  g=CanonicalGraph();g.add_edges_from([('z','a'),('a','m')]);g.edges['z','a']['coefficient']=-2
  sub=g.subgraph({'z','m','a'});self.assertEqual(list(sub),['a','m','z']);self.assertEqual(sub.edges['a','z']['coefficient'],-2)
 def test_composition_rotation_preserves_identity_and_distances(self):
  import composition as c
  p=dict(id='p',reaction_ids=['R1','R2'],positions={'r:R1':[0,0],'m:A':[2,0],'r:R2':[2,3]},reaction_mapping={'s':'R1'},source='toy')
  original=copy.deepcopy(p);q=c.rotate(p,90);self.assertEqual(p,original);self.assertEqual(q['reaction_mapping'],p['reaction_mapping']);self.assertEqual(q['reaction_ids'],p['reaction_ids']);self.assertTrue(q['algorithmic_pose']['source_instance_unchanged'])
  for a in p['positions']:
   for b in p['positions']:self.assertAlmostEqual(np.linalg.norm(np.array(p['positions'][a])-p['positions'][b]),np.linalg.norm(np.array(q['positions'][a])-q['positions'][b]))
 def test_composition_requires_exact_coverage_sets(self):
  import composition as c
  def proposal(ids,chemical=True):return dict(reaction_ids=ids,complete_chemical_equation_verified=chemical)
  self.assertTrue(c.valid_selection([proposal(['A','B'])],{'A','B'},{'A'}));self.assertFalse(c.valid_selection([proposal(['A','C'])],{'A','B'},{'A'}));self.assertFalse(c.valid_selection([proposal(['A','B'],False)],{'A','B'},{'A'}));self.assertFalse(c.valid_selection([proposal(['A','B']),proposal(['B','C'])],{'A','B'},{'A'}))
 def test_composition_guard_rejects_hidden_regression(self):
  import composition as c
  baseline=dict(node_overlap_pairs=0,label_overlap_pairs=2,straight_crossings=10,canvas_area_per_reaction=100,routing=dict(metrics=dict(proper_crossings=20)),constraints=dict(shape_vector_rmse=1))
  better=copy.deepcopy(baseline);better['straight_crossings']=5;self.assertTrue(c.acceptable(better,baseline));better['label_overlap_pairs']=3;self.assertFalse(c.acceptable(better,baseline));better=copy.deepcopy(baseline);better['canvas_area_per_reaction']=151;self.assertFalse(c.acceptable(better,baseline))
 def test_composition_context_ignores_currency_hub(self):
  import composition as c
  model=dict(metabolites=[dict(id=k,name=k,compartment='c') for k in ('A_c','B_c','atp_c')],reactions=[dict(id='R1',metabolites={'A_c':-1,'atp_c':-1}),dict(id='R2',metabolites={'B_c':1,'atp_c':1})])
  context=c.Context(model,model,dict(modules=[]));blocks=[dict(reaction_ids=['R1'],positions={'r:R1':[0,0],'m:A_c':[-50,0]}),dict(reaction_ids=['R2'],positions={'r:R2':[500,0],'m:B_c':[550,0]})]
  self.assertEqual(context.cost(blocks)['interface_links'],0)
 def master_fixture(self,wrong_alias=False):
  import master_map_audit
  xml=f'''<sbml xmlns="http://www.sbml.org/sbml/level2/version4" xmlns:cd="http://www.sbml.org/2001/ns/celldesigner"><model><annotation><cd:extension><cd:listOfSpeciesAliases><cd:speciesAlias id="a1" species="s1"><cd:bounds x="0" y="0" w="20" h="20"/></cd:speciesAlias><cd:speciesAlias id="a2" species="s1"><cd:bounds x="10" y="0" w="20" h="20"/></cd:speciesAlias><cd:speciesAlias id="b1" species="s2"><cd:bounds x="100" y="0" w="20" h="20"/></cd:speciesAlias></cd:listOfSpeciesAliases></cd:extension></annotation><listOfSpecies><species id="s1" name="a[c]" compartment="default"/><species id="s2" name="b[c]" compartment="default"/></listOfSpecies><listOfReactions><reaction id="r1" name="R"><annotation><cd:extension><cd:baseReactants><cd:baseReactant species="s1" alias="{'b1' if wrong_alias else 'a1'}"/></cd:baseReactants><cd:baseProducts><cd:baseProduct species="s2" alias="b1"/></cd:baseProducts></cd:extension></annotation><listOfReactants><speciesReference species="s1" stoichiometry="2"/></listOfReactants><listOfProducts><speciesReference species="s2"/></listOfProducts></reaction></listOfReactions></model></sbml>'''
  with tempfile.TemporaryDirectory() as tmp:
   path=Path(tmp)/'map.xml';path.write_text(xml);return master_map_audit.parse(path)
 def test_master_preserves_graphical_copies(self):
  d=self.master_fixture();self.assertEqual(len(d['species']),2);self.assertEqual(len(d['aliases']),3);self.assertNotEqual(d['aliases']['a1']['bounds'],d['aliases']['a2']['bounds'])
 def test_master_label_compartment_not_sbml_container(self):
  d=self.master_fixture();self.assertEqual(d['species']['s1']['identity'],'a@c');self.assertEqual(d['species']['s1']['sbml_compartment'],'default')
 def test_master_numeric_stoichiometry_and_sbml_default(self):
  r=self.master_fixture()['reactions'][0];self.assertEqual(r['coefficients'],{'a@c':-2.0,'b@c':1.0});self.assertTrue(r['reversible']);self.assertEqual(r['reversibility_origin'],'SBML_L2_default_true')
 def test_master_rejects_wrong_alias_identity(self):
  d=self.master_fixture(True);self.assertEqual(d['issues']['invalid_or_mismatched_alias_link'],1);self.assertNotIn('b1',[r['alias'] for r in d['reactions'][0]['glyph_links'] if r['kind']=='baseReactant' and r['valid']])
if __name__=='__main__':unittest.main()
