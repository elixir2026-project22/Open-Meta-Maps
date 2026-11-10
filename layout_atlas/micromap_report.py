"""Research extension report, rendered from saved source and target evidence."""
from core import *
import html
import pages,ui
import xml.etree.ElementTree as ET


def source_svg_bounds(text):
    root=ET.fromstring(text);points=[]
    for element in root.iter():
        if local(element)=='polyline':points.extend([list(map(float,p.split(','))) for p in element.get('points','').split()])
        elif local(element)=='circle':points.append([float(element.get('cx')),float(element.get('cy'))])
    if points:
        arr=np.array(points);mn=arr.min(0)-30;size=np.ptp(arr,axis=0)+60
        text=re.sub(r'viewBox="[^"]*"','viewBox="'+' '.join(map(str,[*mn,*size]))+'"',text,count=1)
    return text




def generate(summary=None):
    s=summary or read(OUT/'micromap_extension/summary.json');lib=s['library'];models=s['models']
    table=[];selected=[];draw=[];views=[]
    for model in models:
        rows={r['condition']:r for r in model['conditions']};a=rows['G'];b=rows['G_D'];c=rows['G_D_M']
        table.append(f"| {model['model']} | {a['target_reactions']} | {a['candidate_chemical_reactions']} | {b['candidate_chemical_reactions']} | {c['candidate_chemical_reactions']} | {len(model['new_candidate_chemical_ids'])} |")
        selected.append(f"| {model['model']} | {a['selected_chemical_reactions']} / {a['selected_reactions']} | {b['selected_chemical_reactions']} / {b['selected_reactions']} | {c['selected_chemical_reactions']} / {c['selected_reactions']} | {c['selected_micromap_reactions']} |")
        for r in model.get('drawing',[]):draw.append(f"| {model['model']} | {r['condition']} | {r['reactions']} | {r['visible_nodes']} | {r['node_overlap_pairs']} | {r['label_overlap_pairs']} | {r['straight_crossings']} | {r['routing']['metrics']['proper_crossings']} | {r['canvas_area_per_reaction']:.0f} |")
        for r in model.get('fixed_views',[]):views.append(f"| {model['model']} / {r['view']} | {r['condition']} | {r['reactions']} | {r['visible_nodes']} | {r['node_overlap_pairs']} | {r['label_overlap_pairs']} | {r['straight_crossings']} |")
    table='\n'.join(table);selected='\n'.join(selected);draw='\n'.join(draw);views='\n'.join(views)
    text=f'''# Independent MicroMap extension library: semantic mapping, geometry integration and source-increment experiment

Date: 2026-11-06. Generated from saved results. Follows the [07 survey](07_master_map_reuse.md), the [08 plan](08_next_version_plan.md) and the [09 engineering experiment](09_p0_geometry_results.md). This round implements a first closed loop for P1; it does not claim to have solved branched CellDesigner replay, every cross-namespace problem, or any learning gain.

## 1. Integration status and licence

MicroMap has moved from survey-only raw files to an **independently searchable research extension package**: a dual-index SQLite database with a model-applicability table, JSON exchange files, source observations, a catalog, equation mapping, module selection and target-graphic export. The original open core remains 17,221 instances, 3,127 families and 6,094 variants; the new package wrote nothing into the original core SQLite or the existing training set.

The master map is pinned to dataset v4.0, DOI 10.7910/DVN/FZKMJ8, source SHA256 `{lib['source_sha256']}`. The extension layouts inherit CC BY-NC-SA 4.0; they are stored separately from the original core, and combined output keeps each module's licence, design source and mapping. A paper's licence does not replace the data licence. [Pinned source data](https://dataverse.harvard.edu/api/datasets/10740668/versions/4.0).

| Extension package metric | Count |
|---|---:|
| Independent human master-map designs | {lib['independent_designs']} |
| Raw reaction records | {lib['source_reactions']} |
| Reactions entering non-branched geometry decoding | {lib['decoded_reactions']} |
| Source graphical copies kept | {lib['glyphs']} |
| Module instances | {lib['modules']} |
| Structural families in the extension package | {lib['families']} |
| Layout variants in the extension package | {lib['variants']} |

Module size distribution: `{lib['sizes']}`. Geometric families are grouped by reaction/metabolite roles and connectivity, with orientation and chemical content stored separately; families/variants come from bounded isomorphism and geometric comparison, and cannot be added directly to the core family counts to claim the same number of new independent structures. Variant retrieval compares only the first 12 representatives and at most 4 role mappings, and may over-split; these are records of the current algorithm, not human-confirmed counts of unique drawings. All instances still come from one master map and cannot be counted as 12,517 new human-drawn maps or as new species.

## 2. Geometry recovery and limits

Source graphics, connection endpoint centres, local control points and the segment carrying the process node are read from CellDesigner non-branched direct connections. Local points are transformed to canvas coordinates using the endpoint vector and its perpendicular; the reaction centre is recovered from the midpoint of the source segment and marked as a **derived coordinate**, and does not enter supervision of observed expert centres. The original control-point and anchor fields are always kept.

Decoding audit: `{lib['parser_audit']}`. Branches, abnormal references, participant cancellation and similar problems have per-case exclusion records. Only decodable source centre lines are drawn; graphic boundary clipping, complex branches and full replay of side connections are not yet implemented and are not filled in with fake paths.

The geometry rules follow the [CellDesigner version 4 specification, sections 11.1.3 and 11.9](https://www.celldesigner.org/documents/CellDesigner4ExtensionTagSpecificationE.pdf) and the [cd2sbgnml author geometry implementation](https://github.com/sbgn/cd2sbgnml/blob/a35eff8f98b3c831f0d4b652e47bca155b3d05e5/src/main/java/fr/curie/cd2sbgnml/graphics/GeometryUtils.java), accessed 2026-11-06. That implementation serves as a cross-check source for the transformation conventions; none of its converter capabilities is claimed as implemented in this project.

Source-side participants, source ports and all graphical copies are kept in the observation layer. When the first version matches against the current target identity graph and several internal source copies would have to be collapsed into one target core node, the candidate is explicitly rejected and the reason recorded; the first copy is not silently picked. Source paths are stored and replayed, but final target connections are still re-routed on the target network; the source polylines are not claimed to be migrated unchanged.

The subsystem Excel measurement holds only a name list without reaction membership, so reactions are not automatically given named-pathway labels. Modules are formed from the graphical connections actually present in the source map, and include single-reaction units, adjacent reaction pairs and 4/8/12/24-reaction windows under a fixed root budget.

## 3. Identity bridging and fair comparison

This round builds candidate bridges from the BiGG, ChEBI and MetaNetX annotations of the pinned local models, keeping many-to-many candidate relations and the original file SHA. Underscore-escaped spellings such as `mal_L` and `mal__L` are unified as identifier candidates, but D/L, case and all chemical tokens are preserved; entities are not force-merged on name similarity.

For reaction-level candidates the complete participant set, numeric stoichiometry, consistent ratios and direction are checked; bounded bipartite matching keeps unique complete correspondences, and ambiguity is not forced. Source compartment roles may not collapse into the same compartment in the target, and compartment correspondences inside a module must be consistent. Source `c/e/p` come from labels rather than true membrane containers in the source SBML; exact label correspondences and role remapping are recorded separately, and MicroMap candidates are not counted as strict matches with verified true source compartments.

The official MNXref and Rhea resources were checked, but this round did not separately download the full MNXref release or update the existing Rhea 143. Gains in the report come from the actual annotation bridges of the pinned models and equation checks, and cannot be attributed to resources not yet integrated. The [MNXref rules](https://www.metanetx.org/mnxdoc/mnxref.html) state that its generic compartments and undirected representation still need checking against the original equations; [Rhea downloads](https://www.rhea-db.org/help/download) provide direction and obsolete-ID tables. Accessed 2026-11-06.

To distinguish an improved mapping from an added data source, three conditions are set up:

- **G**: the original candidate pool, using the P0-corrected physical scale, joint solving of graphical copies and label avoidance.
- **G_D**: G plus the candidates obtained by running the new direct equation mapper on the existing core library.
- **G_D_M**: G_D plus MicroMap; the source library and the extension library use the same mapping, composition cap and scoring rules.

At most 6,000 candidates per library and at most 16 reaction-correspondence combinations per module; whether the total candidate budget is hit is recorded; the per-module combination cap is 16, and the number of combinations not traversed has not yet been counted module by module. A bounded search finding nothing cannot be interpreted as the whole source definitely having no drawing. Conservative reaction reuse is allowed, outgoing ports of the fixed subgraph are computed against the complete target GEM, and unknown source boundaries are retained; all target participants and stoichiometry still come from the target model; no target expert coordinates or main/side labels are used.

## 4. Candidates and selection under full-model denominators

In the table below, "chemical candidate coverage" means that at least one qualifying multi-reaction module candidate exists; it does not equal final drawing quality or true compartment validation. The new column is the difference in candidate reaction sets between G_D_M and G_D, which removes false additions that were merely "not selected from the old library".

| Target | All reactions | G chemical candidates | G_D chemical candidates | G_D_M chemical candidates | New MicroMap candidate reactions |
|---|---:|---:|---:|---:|---:|
{table}

Selection uses the same greedy plus local-swap heuristic. Each cell below is "selected equation-compatible chemical reactions / all module-arranged reactions"; a geometric match is not passed off as a chemical correspondence.

| Target | G | G_D | G_D_M | Reactions finally arranged by MicroMap modules |
|---|---:|---:|---:|---:|
{selected}

These are application retrieval/composition statistics, not a cross-design learning test. The large models have not yet had all scopes drawn under the new package, so module-composition coverage cannot be called a finished whole-graph drawing; only the fixed display scopes below received new geometry export. The existing full views of the large models keep the original implementation's results. Every model also has a complete reaction-status list; entries whose geometry layout was not executed are marked as selected but not drawn, or pending, and are not silently omitted.

## 5. Actual drawing results: no hiding of no-gain or regression

E. coli core uses the full 95-reaction scope:

| Target | Condition | Reactions | Displayed nodes | Node proximity pairs | Label overlaps | Straight crossings | Orthogonal crossings | Area/reaction (px²) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
{draw}

The large models use the first, middle and last scopes of the existing fixed-view list, with the selection rule fixed before looking at geometry metrics. Each scope compares G_D with G_D_M, and the actual reaction counts are listed row by row; these subgraphs do not represent the complete GEM:

| Scope | Condition | Reactions | Displayed nodes | Node proximity pairs | Label overlaps | Straight crossings |
|---|---|---:|---:|---:|---:|---:|
{views}

The two conditions use the same copy-generation rules, but the selected modules differ, so the actual graphical nodes and connection graph may still differ; the table lists the displayed node counts explicitly. This is an end-to-end source/selection comparison, and crossing differences cannot be attributed to the solver alone.

If the added library increases candidates/biochemical evidence but makes module selection and geometry worse for some target, that is a result this round must keep. It indicates that candidate quality, composition and interface cost still need improvement; library size or a single coverage number cannot be used to claim an overall drawing improvement.

## 6. Files, rerun and follow-up

```sh
python3 project.py micromap
python3 project.py micromap --model e_coli_core
python3 project.py verify
```

- Independent extension SQLite and JSON: `layout_atlas/data/extensions/micromap/`. There, `modules/reactions` index the source layouts and `model_applicability/applicability_reactions` index the target correspondences; a target model can reuse the same instance without producing new expert-instance counts.
- [Source catalog](../layout_atlas/results/micromap_extension/catalog/index.html).
- [Full results and source/code versions](../layout_atlas/results/micromap_extension/summary.json).
- Each target's `equation_mapping_audit.json.gz`, per-library candidates, per-condition selections, reaction status, skeleton and graphical output are in `results/micromap_extension/models/<model>/`; branched and invalid entries are in the [exclusion record](../layout_atlas/results/micromap_extension/excluded.json).
- All derived colour maps of the master map share the same design-source ID and produce no independent training/test sources. This round trained no new ranker or completion model.

Follow-up priorities are decided from the saved failure reasons: if reaction/metabolite annotation is insufficient, add identity bridges; if matching succeeds but no multi-reaction block can be formed, improve scale and coverage selection; if added blocks cause interface conflicts, optimize composition rather than piling up more files. Complex branched geometry, graphical-copy matching and true compartment boundaries remain to be done. Only under a unified candidate and graphical contract should the learning comparison be rerun.
'''
    path=ROOT/'Proposal/10_micromap_extension_results.md';path.write_text(text)
    folder=OUT/'micromap_extension'
    (folder/'index.html').write_text(pages.micromap_index(models))
    citation=read(DATA/'extensions/micromap/SOURCE.json')
    credit='<p>Dataset creator: '+html.escape('; '.join(citation['authors']))+'; '+html.escape(citation['title'])+' v4.0. <a href="https://doi.org/10.7910/DVN/FZKMJ8">Source DOI</a>. Adapted source center-lines and algorithm-derived process centers. CC BY-NC-SA4.0.</p>'
    catalog=folder/'catalog/index.html';content=catalog.read_text();marker='<h1>MicroMap research layout extension</h1>'
    if credit not in content:content=content.replace(marker,marker+credit)
    content=re.sub(r'<svg\b.*?</svg>',lambda match:source_svg_bounds(match.group(0)),content,flags=re.S)
    import sqlite3,database
    db=sqlite3.connect(DATA/'extensions/micromap/atlas.sqlite')
    for instance in re.findall(r'<small>Instance ([0-9a-f]+)',content):
        row=db.execute('SELECT payload FROM modules WHERE id=?',(instance,)).fetchone()
        if row:
            module=database.unpack(row[0]);labels=', '.join(r['id'] for r in module['reactions']);caption='<p>Source reaction IDs: '+html.escape(labels)+'</p>'
            anchor='<small>Instance '+instance
            if caption+anchor not in content:content=content.replace(anchor,caption+anchor)
    db.close();catalog.write_text(content)
    from xml.sax.saxutils import escape
    for svg in (folder/'catalog').glob('*.svg'):
        text=source_svg_bounds(svg.read_text());svg.write_text(text)
        if '<metadata>' not in text:
            offset=text.index('>')+1;svg.write_text(text[:offset]+'<metadata>'+escape(json.dumps(citation,ensure_ascii=False))+'</metadata>'+text[offset:])
    ui.rewrap(catalog,'MicroMap research layout extension',section='micromap');ui.localize_micromap(catalog);return path

if __name__=='__main__':print(generate())
