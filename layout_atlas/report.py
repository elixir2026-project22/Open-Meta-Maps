"""Report and diagnostic figures generated exclusively from saved results."""
import gzip,html
from core import *
import pages



def table(headers,rows):return '| '+' | '.join(headers)+' |\n|'+ '|'.join(['---']*len(headers))+'|\n'+'\n'.join('| '+' | '.join(str(x) for x in row)+' |' for row in rows)
def generate():
 db=read(OUT/'database_summary.json');matched=read(OUT/'matched_corpus_comparison.json');ranks=read(OUT/'ranking_summary.json');models=read(OUT/'model_summary.json');sensitivity=read(OUT/'variant_sensitivity.json');assert all(r['status']=='ok' for r in models);tasks=read(OUT/'completion_tasks.json');results=read(OUT/'completion_results.json');training=read(OUT/'egnn_training.json');ref=read(OUT/'reference_source_audit.json');verification=read(OUT/'verification.json')
 agg=defaultdict(lambda:defaultdict(list));reads=defaultdict(lambda:defaultdict(list));seedagg=defaultdict(lambda:defaultdict(list))
 for r in results:
  if r['status']!='ok':continue
  key=(r['held_source'],r['kind'],r['method']);agg[key][r['map_id']].append(r['hidden_rmse_local_px']);seedagg[key][r['seed']].append(r['hidden_rmse_local_px'])
  for metric in ('straight_crossings','node_overlaps','hint_drift'):reads[key+(metric,)][r['map_id']].append(r[metric])
 mean=lambda v:float(np.mean([np.mean(x) for x in v.values()]))
 cr=[]
 for source in ('community','reactome','wikipathways'):
  for kind in ('expert','perturbed','retrieved'):
   cr.append([source,kind,len(agg[source,kind,'egnn']),*[f'{mean(agg[source,kind,method]):.2f}' for method in ('initial','fcose','ridge','egnn')]])
 mr=[];diagnostics=[]
 for name in sorted({r['model'] for r in models}):
  a=next(r for r in models if r['model']==name and r['condition']=='semantic_only');b=next(r for r in models if r['model']==name and r['condition']=='all_geometry');mr.append([name,b['reactions'],b['strict_equation_reactions'],b.get('chemical_equation_reactions','not saved'),b['reference_compatible_reactions'],a['module_reactions'],b['module_reactions'],b['selected_geometry_database_reactions'],b['bridge_reactions'],b['pending'],b['searched_reaction_union']])
  folder=OUT/'models'/name/'all_geometry';state=read(folder/'skeleton.json');model=load_model(ROOT/state['target_file']);transport=0;exchange=0;metas={m['id']:m for m in model['metabolites']}
  for r in model['reactions']:
   transport+=len({metas[m].get('compartment','unknown') for m in r['metabolites']})>1;exchange+=len(r['metabolites'])==1
  with gzip.open(folder/'candidates.json.gz','rt') as f:ps=json.load(f)
  cumulative=[];seen=set()
  for src in ('community','reactome','wikipathways'):
   seen.update(rid for p in ps if p['source']==src for rid in p['reaction_ids']);cumulative.append(dict(added_source=src,candidate_reaction_union=len(seen),fraction=len(seen)/b['reactions']))
  cur=currency(model);owner=defaultdict(list);mg=nx.Graph();mg.add_nodes_from(range(len(state['modules'])))
  for i,block in enumerate(state['modules']):
   for r in model['reactions']:
    if r['id'] in block['reaction_ids']:
     for mid in r['metabolites']:
      if mid not in cur:owner[mid].append(i)
  for ids in owner.values():
   for i in ids[1:]:mg.add_edge(ids[0],i)
  fragments=nx.number_connected_components(mg)
  for bridge in state['bridges']:mg.add_edge(bridge['from_module'],bridge['to_module'])
  fragments_after=nx.number_connected_components(mg)
  diagnostics.append(dict(skeleton_fragments_before_bridges=fragments,skeleton_fragments_after_bridges=fragments_after,model=name,total=b['reactions'],transport_reactions=transport,single_participant_exchange_like=exchange,pending_reasons=dict(Counter(state['unmatched_reason'].values())),source_candidate_accumulation=cumulative,interpretation='Candidate unions under a common full-library query set; not a causal source-removal experiment. Transport/exchange categories overlap.'))
 dump(OUT/'coverage_diagnostics.json',diagnostics)
 ranktable=table(['Held-out source','Training preference queries','Multi-candidate test queries','Rule','Logistic regression','Gradient boosting','Candidate upper bound'],[[r['held_source'],r['training_preference_queries'],r['test_multivariant'],*[f'{r["map_macro_errors"][m]:.4f}' for m in ('rule','logistic','gradient_boosting','oracle')]] for r in ranks])
 modeltable=table(['Model','All reactions','Strict equation with compartments','Chemical equation','Reaction-identity support','Semantic-library modules only','Plus geometry-library modules','Q2 reactions actually selected','Algorithmic bridges','Pending layout','Searched reactions'],mr)
 comptable=table(['Held-out source','Anchor condition','Test maps','Neighbourhood start','fCoSE','Ridge','EGNN'],cr)
 readability=table(['Source (retrieval anchors)','fCoSE crossings','EGNN crossings','fCoSE node overlaps','EGNN node overlaps'],[[s,*[f'{mean(reads[s,"retrieved",m,metric]):.2f}' for metric in ('straight_crossings','node_overlaps') for m in ('fcose','egnn')]] for s in ('community','reactome','wikipathways')])
 small=next(r for r in models if r['model']=='e_coli_core' and r['condition']=='all_geometry');st=table(['Method','Straight crossings','Label overlaps','Node overlaps','Orthogonal-routing crossings'],[[r['method'],r['straight_crossings'],r['label_overlap_pairs'],r['node_overlap_pairs'],r['routing'].get('metrics',{}).get('proper_crossings','failed')] for r in small['full_small_comparisons']])
 strict_status=verification['status']
 text=f'''# Expert layout database v0.3: semantic recovery, multi-candidate skeletons and paired completion

Date: 2026-11-06. Sources were acquired across 5-6 November; this report is generated by `layout_atlas/report.py` from saved results. The goal stays **A database -> B1 initial skeleton -> B2 completion of the remaining known reactions in the target**. A runnable technical closed loop is now complete. Quality acceptance uses source tracing, automatic consistency checks and reproducible experiments, with no two-person review step; learning gain is not yet stable.

## 1. Version, evidence and delivery status

The current version runs under a standalone `layout_atlas/`; old code, old reports and historical outputs have been deleted. Necessary raw material and source observations are stored separately in `data/source_inputs/`, and this round verified {verification['checks']['source_input_files_hashed']:,} input files. Database instances, families, variants, source maps, independent drawing authors and the number of human acceptances are counted separately; a window is not an independent sample.

| Deliverable | Status and entry point |
|---|---|
| Dual-index database | `../layout_atlas/data/atlas.sqlite`; the `modules` layout index and the `reactions` identity-candidate index |
| Instance and variant catalog | [static catalog](../layout_atlas/results/catalog/index.html), with the full `modules.json.gz` exchange file |
| Source and mapping audit | `data/acquisition.json`, `results/binding_audit.json`, `reference_source_audit.json` |
| Six-model skeletons | `results/models/<model>/<condition>/`; all reaction states, interfaces, module mappings, hierarchical views, SVG / Escher JSON |
| Ranking and completion | LR / gradient boosting actually trained, three paired anchor tasks, Ridge / fCoSE / 4-layer 64-dimensional EGNN, 9 CPU checkpoints |
| Automatic acceptance | `{strict_status}`; {verification['checks']['unit_tests']} behaviour tests, plus source-isolation, original-model and export-integrity checks |

## 2. A: what the existing material actually recovered

The GPML Interaction Xref recovered **{db['rhea_occurrences']} Rhea references, {db['rhea_unique_ids']} normalized IDs, {db['rhea_maps']} maps and {db['rhea_species']} species**. The 942 from the earlier raw string scan and the 941 from this normalization are not two different numbers of biochemical reactions; direction IDs, duplicate references and string normalization must each be considered separately.

{table(['Binding status','Occurrences'], sorted(db['rhea_bindings'].items()))}

Only bindings whose visible participants are compatible with the reference reaction and whose direction is interpretable were promoted. The 1008 occurrences with only a reaction ID keep their evidence but are not counted as 1008 qualified single-step layouts. Dashed aggregate processes, obsolete IDs, unclear arrowheads and inconsistent directions were not forcibly split or used to overwrite labels.

Uses the RDF, direction table, obsolete-ID table, cross-references and pH 7.3 mapping of **Rhea release 143**; 19,005 reactions were parsed, of which 15,289 have fully resolvable participants and numeric stoichiometry. The pH correspondence keeps `curation/computation`. Multiple ChEBI candidates remain candidate sets; participant and equation checks take place on specific source-target pairs. The [official Rhea download notes](https://www.rhea-db.org/help/download) clearly distinguish the direction table, the obsolete-ID table and the source of the protonation mapping.

The real nodes, paths, anchors, enzyme links and groupings of GPML are kept; reaction centres inserted along a path are marked `derived_from_observed_path` and are not used as expert reaction-coordinate supervision. When the original map draws only the main participants, a "complete reaction + partially visible layout" is stored, and the remaining target participants are arranged algorithmically at export. Midpoints of pure-geometry links are likewise derived routing points. A multi-step corridor cannot pass as several expert single-step reactions; see the [WikiPathways indirect-step tutorial](https://academy.wikipathways.org/stages/met-1-indirect/index.html).

### Reactome pinned version and layered reference

The metabolic hierarchy was enumerated and **115/115 pairs of Reactome v97 layout JSON + graph-structure files** were obtained, pinned by the [official archive record 21383214](https://zenodo.org/records/21383214). ZIP range reads were used, with per-member length, CRC32 and SHA256 checks; it is not claimed that the published MD5 of the whole 4.7 GB archive was verified. Earlier server-side 429 failures are kept in the historical log, and the number of extraction failures in the last round was zero. Only human-curated maps are used; maps inferred for other species are not treated as independent designs.

The adapter promotes only records whose graph event is a Reaction with clear input/output correspondence and connection stoichiometry. Public Reactome->Rhea cross-references are added to the candidate index, but they are not equivalence assertions that need no checking.

The upstream Human/Yeast Escher maps total **{ref['escher_reference_files']}** files (Human {ref['by_source']['human_maps']}, Yeast {ref['by_source']['yeast_maps']}); identity and geometry indexes were extracted for reference. The CellDesigner file of the historical yeast commit `93b787c6986a834cac323de667aac4b93547a457` was obtained: {ref['historic_yeast']['reaction_elements']} reaction elements and {ref['historic_yeast']['species_aliases']} speciesAlias entries. Licence conflicts and the per-element "kept / moved / added / auto-generated" correspondence remain unresolved, and **everything stays in the reference layer and does not enter core training or the module library**. The historical query for the drawing-rule path was empty; it is still unrecovered evidence and cannot be cited as a rule that has been read.

### Dual index and over-splitting check

The currently active library has **{db['module_count']:,} instances, {db['families']:,} core structural families and {db['variants']:,} layout variants**; {db['singleton_families']:,} families are singletons. By source: community {db['by_source']['community']:,}, Reactome {db['by_source']['reactome']:,}, WikiPathways {db['by_source']['wikipathways']:,}. Of these, `geometry_only` is {db['by_quality']['geometry_only']:,} and `semantic_transfer` is {db['by_quality']['semantic_transfer']:,}. The latter means reaction semantics are retained, not that every coordinate was directly observed.

These instances use {db['used_source_layouts']} source-layout records; the registry's {db['registered_map_records']} entries include historical metadata and cannot be treated as an equal number of independent expert designs. Full author-level deduplication across websites and historical versions is not yet done.

On the same batch of **{matched['same_instances_compared']:,} retained instances**, the old full-display families number {matched['legacy_full_display_families']:,} and the new core-shape families {matched['v03_core_shape_families']:,}; variants go from {matched['legacy_layout_variants']:,} to {matched['v03_core_variants']:,}. This shows that separating side participants changes the granularity of classification, **not that the merging is necessarily correct**. The core-shape index ignores chemical identity; direction, compartment and equation are checked in the biochemical index and retrieval evidence. Pure topological compatibility does not count as biochemical equivalence. Rotation, translation, scale and mirroring are aligned on real instances; representatives are actual instances, and no "expert map" is synthesized by averaging.

Another {matched['retained_not_in_active_matcher']} old observations did not enter the active matcher (disconnected visible core 2676, pinned-version replacement 88, scale out of bounds 174); they are kept in `retained_observations` and were not silently deleted. The graph-isomorphism budget and the cap of 24 variant representatives lead to conservative splitting. A fixed draw of the first 100 families with at least 3 instances, the first 8 instances per family, {sensitivity['instances']} instances in all: thresholds 0.15/0.25/0.35 give {sensitivity["threshold_variant_counts"]["0.15"]}/{sensitivity["threshold_variant_counts"]["0.25"]}/{sensitivity["threshold_variant_counts"]["0.35"]} variants respectively; the maximum residual of the rigid, mirror and scale transforms is about 1.51e-15. This bounded sample shows that the variant count depends on the threshold; it cannot be taken as a library-wide sensitivity conclusion, nor does it replace manual review.

![Full-model coverage](../layout_atlas/results/figures/coverage.png)

## 3. B1: how multi-candidate retrieval becomes a skeleton

Retrieval keeps four levels of evidence: reaction identity, complete equation, role core and pure geometry. Generic windows plus target windows derived from source reaction identities prevent reliable biochemical correspondences from failing to enter the candidates at all because of different BFS window ranges. Each query keeps at most 24 distinct geometric arrangements with explicit node / reaction mappings. A source core may cover only the main participants of the target reaction, but the additional participants in the target are kept and marked as awaiting algorithmic arrangement.

Multi-reaction modules that share no internal reactions are selected through weighted coverage, interface cost and one local swap. Bridging searches the **real target bipartite graph**; a template never adds a source reaction to the target. ELK lays out the compartment composition containers, module ports and shared chemical-identity junctions, and the module-local copies are then instantiated. Cross-compartment transport uses composition containers, not real membrane geometry; full nested organelle boundary constraints are not implemented. The [official ELK capabilities](https://eclipse.dev/elk/reference/algorithms/org-eclipse-elk-layered.html) are backend capabilities and do not mean this prototype has solved every compartment problem.

{modeltable}

"Strict equation with compartments" requires explicit and consistent source compartments, and `local/unknown` does not count as verified; "chemical equation" checks only complete participants, stoichiometry and direction. These two and "reaction-identity support" may overlap and cannot be added; module coverage includes geometric arrangements and cannot be called strict biochemical coverage. The Q2 column counts only the `geometry_only` source templates actually selected. The semantic library already contains expert geometry, so adding the Q2 library does not guarantee more reaction coverage. Pending layout means reactions the B1 skeleton has not yet arranged; the export views fill them in algorithmically, which does not mean the reaction was not exported. The full-model denominator keeps all exchange, transport and complex reactions; retrieval is also limited by a nominal 120 generic roots (step sampling may give slightly more) and a 600 biochemical-window budget. "Searched" and "the library definitely has no drawing" are different concepts; per-reaction reasons for non-coverage are in `coverage_diagnostics.json` and the skeleton list.

All six-model results are **application coverage**, using qualified source libraries, and are not independent generalisation scores. The independent learning evidence comes from the leave-one-source-out experiments below. The source-accumulation diagnostic gives an upper bound of candidate coverage under a common query set and cannot be read as a causal per-source acquisition-gain curve. Views are partitioned in advance by target compartment and true connectivity, with at most 50 reactions per figure and every reaction appearing exactly once; an 80-reaction cut-out no longer stands in for the full model.

### Learning geometry selection

Queries provide no expert coordinates or main/side labels. Candidate-attachment features come from the generic initial graph, and labels are used only to compute similarity-transform error at observed points. Sources are held out, and queries and layout-derived copies are excluded from the training candidates. All three folds exceed the engineering thresholds of 50 training preference queries and 20 multi-candidate test queries, and logistic regression and gradient boosting were actually trained; this threshold is not proof of statistical sufficiency.

{ranktable}

The values are normalized shape errors macro-averaged over maps, lower is better; the candidate upper bound is for evaluation only. Improvement varies with source and model, and a stable advantage over the rule cannot be claimed. The local ranking test is mainly geometry selection and should not pass as a validated cross-namespace biochemical ranking model. The production CLI can load the source-fold selector, but overall coverage and composition remain heuristic and no composition objective was learned. In another E. coli production run a gradient-boosting model from the leave-community-out fold was loaded; its result is stored separately as `learned_geometry_community/`, only as an interface-runnability comparison, and is not mixed into the independent test scores above.

## 4. B2: fair tasks, actual models and results

Source isolation comes first, then masking and the retrieval library are constructed. Candidates exclude the query's own source lines and duplicate layouts. Of the 36 pre-selected queries, 53 query-fold combinations across the three folds have no qualified full-core candidate, forming **{len(tasks['tasks'])} tasks** that actually involve {len({t['map_id'] for t in tasks['tasks']})} maps; the source maps available for independent testing total 17. Failures with no candidate are kept, and a condition success rate alone must not be reported.

For each query / mask / source fold, expert, perturbed and retrieved use the **same visible/hidden nodes**. The three fixed tasks this round are random 20%, contiguous branch 50% and module-region 80% masking; the Cartesian product of all ratios and mask types is not yet covered. The graph keeps reaction and metabolite types, direction / reversibility, stoichiometry, participant roles, compartment fields and graphical-copy keys; bidirectional message passing is still used, and edge features keep the original chemical roles. Derived centres do not enter hidden-coordinate supervision.

The frame is determined by the visible anchors; the aligned expert hidden coordinates exist only in labels, and the prediction function runs without labels. Expert/perturbed anchors are hard-fixed, retrieval anchors are softly constrained by confidence; the output is then adjusted with neighbourhood, spacing and hint constraints. Explicitly provided compartment rectangles can limit movable nodes; the current local tasks have no reliable membrane boundary, so the compartment-violation metric is recorded as undefined and must not be reported as zero.

A 4-layer, hidden-64-dimensional 2D EGNN residual network was implemented, on CPU with three fixed seeds and early stopping on source-grouped validation; hard anchor constraints are preserved layer by layer. Rotation / translation / mirror equivariance has unit tests, but the deterministic ordering in graph initialization is not a guarantee of full end-to-end equivariance. The [original EGNN paper](https://proceedings.mlr.press/v139/satorras21a.html) did not validate metabolic-graph drawing, and this project's results cannot borrow the paper's performance as endorsement.

{comptable}

Units are coordinate units at an RMS scale of 100 in the anchor frame, not original canvas pixels. Tasks and seeds are averaged within each map first, then macro-averaged over maps. There are very few test maps, and overlapping windows or the three seeds are not treated as independent replicates. Derived reaction centres have no RMSE label; per-node-type results are in the raw JSON. EGNN does better on community, while on Reactome it is not consistently better than fCoSE, so a general gain cannot yet be shown.

{readability}

The table above measures core-graph positions only and cannot replace readability of complete labels, side participants and paths. The three full-scope views of E. coli have additionally been run through "actual B1 skeleton -> label-free EGNN prediction -> target chemistry export", with file names `view-*-egnn.*`, for deployment feasibility and not as an end-to-end generalisation conclusion. Side metabolites, labels and paths are still at the rule stage.

![Paired completion results](../layout_atlas/results/figures/completion.png)

## 5. Full case and routing: regressions explicitly kept

E. coli core, all 95 reactions, same target-participant scope:

{st}

After the skeleton, straight-line crossings and label overlaps fall, but **node overlap rises**. libavoid produces orthogonal paths without guaranteeing fewer crossings; avoiding node boxes is not the same as avoiding labels, and collinear coincidences are counted separately. The current geometry results still fall short of an "ideal database + automatic publication-grade drawing". The quality of nodes, labels and paths should not be offset by a single RMSE.

Fixed case-selection rule: the full small model uses all of E. coli, with no sub-region picked; the coverage gap uses the model with the highest share pending layout among full models; ranking disagreement takes the first query by lexicographic ID; completion improvement / regression each take the first case by lexicographic ID under the same condition. See `results/cases/`. This avoids showing only the best pictures; readability conclusions rest on the saved objective metrics, with no two-person blind-rating threshold.

## 6. Acceptance, remaining items and next steps

Passed source verification, completeness of model reaction states, fidelity of stoichiometry / direction / participants, actual ports, source-fold isolation, paired masking, no hard-anchor drift, and the Escher schema check. For full counts see `verification.json`. When the source coordinate frame is unclear, compartment migration is not declared verified merely because the chemical equations agree. Chemical completeness does not equal readability.

Under the current research scope, the two-person review, the 20+80 review tasks and the two-person blind-rating acceptance threshold are cancelled. Automatic check results are not labelled as human certification. Per-element provenance and file-level licences of Human / Yeast remain unresolved. Not completed: globally optimal arrangement of large models, real membrane-boundary constraints, label-avoiding routing, full verification of source derivation relations, a full sensitivity analysis of the variant threshold, and cross-species stress tests. These must not be written up as implemented capabilities.

Next steps proceed on auditable evidence:

1. Improve the automatic quality audit: first check GPML single-step semantics, pure-geometry corridors, graphical copies, stoichiometry and ports; classify entries by evidence completeness into transferable, geometry-reference and unresolved, with no manual sign-off threshold.
2. For high-gap models, distinguish "missing identity bridge / extra side metabolite or compartment incompatibility / no template for the core structure / insufficient search budget". Extend reliable mappings and high-value species pathways, recovering existing material first rather than blindly adding files.
3. Optimize coverage target, interface quality and block composition separately; add avoidance based on actual label rectangles, module copy-point constraints and real compartment boundaries. Eliminate the node-overlap regression first, then decide whether a learned composition is needed.
4. Re-run fCoSE / EGNN with the same masking on newly added independent source maps, add full-branch and module-joining tasks, and add leave-species-out stress tests; judge gains by map-paired confidence intervals, crossings, overlaps and connection completeness. Without a stable gain, keep using rules and solvers and retain the learning interface.
5. In the second completion sub-stage, study side participants, labels and path arrangement specifically; keep all reactions from the target model and never predict or invent biochemical reactions.

For rerun commands, source versions and implementation boundaries see the [v0.3 README](../layout_atlas/README.md). Database trustworthiness, skeleton usability and learning gain are still accepted separately.
'''
 (ROOT/'Proposal/06_curated_atlas_and_learning.md').write_text(text);dump(OUT/'report_statistics.json',dict(database=db,completion_map_macro={':'.join(k):mean(v) for k,v in agg.items()},model_rows=mr,readability={':'.join(k):mean(v) for k,v in reads.items()},source='Persisted JSON only; report.py'))
 figures(models,ranks,agg);cases(tasks,results)
 import navigation;navigation.generate()
 return dict(report='Proposal/06_curated_atlas_and_learning.md',status='generated')

def figures(models,ranks,agg):
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 dest=OUT/'figures';dest.mkdir(exist_ok=True);names=[r['model'] for r in models if r['condition']=='all_geometry'];ys=[r for r in models if r['condition']=='all_geometry'];fig,ax=plt.subplots(figsize=(10,4));x=np.arange(len(ys));ax.bar(x,[100*r['module_reactions']/r['reactions'] for r in ys],label='Any module arrangement');ax.bar(x,[100*r['strict_equation_reactions']/r['reactions'] for r in ys],label='Complete equation checked');ax.set_xticks(x,names,rotation=15);ax.set_ylabel('% of all target reactions');ax.set_ylim(0,100);ax.legend();fig.tight_layout();fig.savefig(dest/'coverage.png',dpi=170);plt.close(fig)
 fig,axes=plt.subplots(1,3,figsize=(12,3.6),sharey=True)
 for ax,source in zip(axes,('community','reactome','wikipathways')):
  for method in ('initial','fcose','ridge','egnn'):
   val=[np.mean([np.mean(v) for v in agg[source,k,method].values()]) for k in ('expert','perturbed','retrieved')];ax.plot(range(3),val,'o-',label=method)
  ax.set_title(source);ax.set_xticks(range(3),['expert','perturbed','retrieved'],rotation=15)
 axes[0].set_ylabel('Map-macro hidden RMSE (anchor frame)');axes[-1].legend();fig.tight_layout();fig.savefig(dest/'completion.png',dpi=170);plt.close(fig)

def cases(tasks,results):
 import torch,neural,completion,matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 dest=OUT/'cases';dest.mkdir(exist_ok=True);rows={ (r['task'],r['method']):r for r in results if r['status']=='ok' and r['seed']==SEEDS[0]};selected=[]
 for direction in ('improves','worsens'):
  for task in sorted(tasks['tasks'],key=lambda t:t['id']):
   if task['split']!='test' or task['anchor_kind']!='retrieved':continue
   a=rows.get((task['id'],'egnn'));b=rows.get((task['id'],'fcose'))
   if a and b and ((a['hidden_rmse_local_px']<b['hidden_rmse_local_px'])==(direction=='improves')):selected.append((direction,task));break
 cards=[]
 for direction,t in selected:
  model=neural.EGNN();model.load_state_dict(torch.load(OUT/'egnn'/f'{t["held_source"]}-{SEEDS[0]}.pt',weights_only=True,map_location='cpu'));model.eval();egnn=neural.predict(model,t);fcose=completion.baseline(t,'fcose');truth=t['labels']['observed_xy'];g=completion.task_graph(t);fig,axs=plt.subplots(1,3,figsize=(12,4))
  for ax,name,xy in zip(axs,('Observed expert labels','fCoSE','EGNN'),(truth,fcose,egnn)):
   for a,b in g.edges():
    if a in xy and b in xy:ax.plot([xy[a][0],xy[b][0]],[xy[a][1],xy[b][1]],color='#aab8bf',lw=.8,zorder=0)
   keys=sorted(xy);arr=np.array([xy[k] for k in keys]);ax.scatter(arr[:,0],arr[:,1],c=['#d07b31' if g.nodes[k]['role']=='reaction' else '#157f87' for k in keys],s=[60 if k in t['visible'] else 15 for k in keys]);ax.set_title(name);ax.set_aspect('equal');ax.invert_yaxis()
  fig.suptitle(direction+' / first task ID / '+t['held_source']);fig.tight_layout();fig.savefig(dest/(direction+'.png'),dpi=150);plt.close(fig);dump(dest/(direction+'.json'),dict(task=t,fcose=fcose,egnn=egnn,selection='First lexicographic task ID with specified error difference; not best/worst.'));cards.append(f'<h2>{direction}</h2><img style="width:100%" src="{direction}.png"><a href="{direction}.json">complete case</a>')
 with gzip.open(OUT/'ranking_queries.json.gz','rt') as f:rq=json.load(f)
 disagreement=next((r for r in sorted(rq,key=lambda r:(r['held_source'],r['query'])) if r.get('selected',{}).get('gradient_boosting') and r['selected']['rule']!=r['selected']['gradient_boosting']),None);dump(dest/'first_ranking_disagreement.json',disagreement)
 (dest/'index.html').write_text(pages.cases_page(dest,cards))
if __name__=='__main__':print(generate())
