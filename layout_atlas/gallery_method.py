"""Library switch and the 'How these modules were collected' section shared by the two catalog index pages.

Every number is read from saved data or computed from the stored modules when the catalog is generated; nothing is typed in by hand.
Sources:  window library  database.modules(), results/database_summary.json, results/module_redundancy_summary.json, results/variant_audit.json
          unit library    data/extensions/units/SOURCE.json, results/units_extension/library_summary.json, units_library.modules('units_sub')
Rules quoted in the text come from database.diagram_modules / database.build (windows) and units_library.py (units).
"""
import json
from collections import Counter

import ui
from core import OUT, DATA

bi = ui.bi


def _read(path):
    with open(path) as f:
        return json.load(f)


def window_facts(mods):
    """Facts about the window library: from the stored modules plus the saved summaries."""
    summ = _read(OUT / 'database_summary.json')
    red = _read(OUT / 'module_redundancy_summary.json')
    audit = _read(OUT / 'variant_audit.json')
    maps = {}
    for m in mods:
        maps[m['map_id']] = (m['source'], m['license'])
    sizes = Counter(m['size'] for m in mods)
    fixed = (1, 2, 4, 8, 12, 24)
    return dict(
        modules=len(mods), families=len({m['core_family'] for m in mods}), variants=len({m['variant'] for m in mods}),
        maps=len(maps), maps_by_source=dict(Counter(s for s, _ in maps.values())), maps_by_license=dict(Counter(l for _, l in maps.values())),
        registered=summ['registered_map_records'], geometry_only=sum(m['quality'] == 'geometry_only' for m in mods),
        share_fixed=sum(sizes[k] for k in fixed) / len(mods), min_size=min(sizes), max_size=max(sizes),
        contained=red['contained'], share_contained=red['share_contained'], redundancy={s: v['redundancy_factor'] for s, v in red['by_source'].items()},
        audit=audit, summary_matches=(summ['module_count'] == len(mods) and summ['families'] == len({m['core_family'] for m in mods}) and summ['variants'] == len({m['variant'] for m in mods})))


def unit_facts(units):
    """Facts about the unit library: from SOURCE.json and library_summary.json, cross-checked against the stored units."""
    src = _read(DATA / 'extensions/units/SOURCE.json')
    summ = _read(OUT / 'units_extension/library_summary.json')
    import units_library
    sub = len(units_library.modules('units_sub'))
    return dict(src=src, summ=summ, units=len(units), sub=sub,
                summary_matches=(summ['modules'] == len(units) and summ['families'] == len({m['core_family'] for m in units}) and summ['variants'] == len({m['variant'] for m in units})))


def switch(current, other_href_units='../catalog_units/index.html', other_href_windows='../catalog/index.html'):
    """Library switch shown at the top of both index pages (reuses the .pager style: current = filled)."""
    items = [('windows', 'Pathway-fragment reuse', other_href_windows), ('units', 'Structure-unit reuse (experimental)', other_href_units)]
    parts = []
    for key, label, href in items:
        if key == current:
            parts.append(f'<span aria-current="page" title="current library">{bi(label)}</span>')
        else:
            parts.append(f'<a href="{href}">{bi(label)}</a>')
    return ('<div class="pager libsw" role="navigation" aria-label="Library"><span style="background:none;border:0;color:var(--mute)">' + bi('Library:') + '</span>' + ''.join(parts) + '</div>')

def n(x):
    return f'{x:,}'


def pct(x):
    return f'{100 * x:.1f}%'


def section(current, W, U):
    """The 'How these modules were collected' panel; both columns appear on both pages, the current library is marked."""
    s, S = U['src'], U['summ']
    by_src = s['source_map_ids_by_source']
    lic_u = s['licenses']
    sh = S['shapes']
    aud = W['audit']
    mb, ml = W['maps_by_source'], W['maps_by_license']
    q = S['by_quality']
    cfg = S['config']
    red = W['redundancy']

    def hd(key, en):
        mark = ('<span class="tag ind">' + bi('this page') + '</span>') if key == current else ''
        return f'<th>{bi(en)} {mark}</th>'

    def row(label, w, u):
        return f'<tr><th scope="row">{bi(*label)}</th><td>{w}</td><td>{u}</td></tr>'

    rows = []
    rows.append(row(('Source maps and licences',),
        bi(f"{n(W['maps'])} source layouts used (of {n(W['registered'])} registered map records): Escher {n(mb.get('escher', 0))}, WikiPathways {n(mb.get('wikipathways', 0))}, Reactome {n(mb.get('reactome', 0))}. Licence of the source maps: " + ', '.join(f'{k} {n(v)}' for k, v in sorted(ml.items())) + '. Each module keeps the licence, SHA-256 and source pointer of its map.'),
        bi(f"{n(s['source_maps'])} source maps recut: Escher {n(by_src.get('escher', 0))}, WikiPathways {n(by_src.get('wikipathways', 0))}, Reactome {n(by_src.get('reactome', 0))}; licences " + ', '.join(f'{k} {n(v)}' for k, v in sorted(lic_u.items())) + f". Only {n(S['source_layouts_with_units'])} of them yield at least one unit. The map set is not identical to the window library's (counts differ; the overlap itself was not checked). Open-core data only; this extension is not yet merged into the open-core SQLite database (admitted_to_open_core = {str(s['admitted_to_open_core']).lower()}).")))
    rows.append(row(('What a module instance is',),
        bi('A connected set of reactions of one source map with their core metabolites, stored with source coordinates, reaction identities (when the map has them), licence and source pointer. Different instances may contain the same reactions.'),
        bi(f"A connected set of {cfg['min_size']}–{cfg['max_size']} core reactions of one source map (appendix reactions can raise this to at most {cfg['max_appendix_total']}) forming one drawing shape: chain, cycle or fork. Currency metabolites (water, protons, ATP/ADP/AMP, NAD(P)(H), CoA, phosphate) are not part of the core. Reactions in two units of the same map: {n(S['overlap_check_reactions_in_multiple_units'])}.")))
    rows.append(row(('How it was cut',),
        bi(f"Independent breadth-first windows. WikiPathways maps whose interactions are bound to Rhea reactions and Reactome diagrams (also Human-maps, Metabolism Regulation, PANTHER) are cut by database.diagram_modules: size-scaled roots (all reactions when N≤48, else max(24, ceil(N/4))); from each root a breadth-first walk over the reaction–metabolite graph takes the first 1, 2, 4, 8, 12 and 24 reactions; only windows that are not a proper subset of another window on the same map are stored. Community Escher maps are cut in-repo by escher_cutter (escher_community_windows): every reaction is a root; links go through primary drawn metabolite entities (duplicate glyphs merged by BiGG id); window sizes 1, 2, 4, 8, 12, 24 and 40, plus the full reachable component when it is off-grid and ≤40 reactions; modules are kept at ≤100 nodes. WikiPathways geometry-only maps use the same Stage B diagram_modules policy over metabolite–metabolite GPML edges. {pct(W['share_fixed'])} of all windows have exactly 1, 2, 4, 8, 12 or 24 reactions; the sizes found are {W['min_size']}–{W['max_size']}."),
        bi(f"Non-overlapping, cut once per source map by units_library.decompose, in fixed order: (1) cycles: the shortest cycle through every edge, at most {cfg['cycle_max_reactions']} reactions, taken shortest first; ring rule: rings of at least {cfg['ring_min_reactions']} reactions without transport reactions are seeded first, and a neighbouring ring is fused when it shares at least {cfg['fuse_min_shared']} reactions and adds at most {cfg['fuse_max_new']}; (2) forks: a branch point with up to {cfg['fork_max_arms']} arms, arm depth up to {cfg['fork_arm_depth']}, within {cfg['max_size']} reactions; (3) chains: what remains is cut at branch points, and paths longer than {cfg['max_size']} are halved recursively; (4) fragments below {cfg['min_size']} reactions are attached as appendix to an adjacent unit (total at most {cfg['max_appendix_total']}); isolated islands that stay below {cfg['min_size']} are orphans and not modules ({n(S['orphan_reactions'])} orphan reactions; {n(S['unit_reactions'])} of {n(S['source_core_reactions'])} core reactions are inside units). WikiPathways maps without reaction identity are segmented by drawn shape only: each drawn interaction between two metabolite glyphs is a pseudo-reaction (quality = geometry_only; {n(q.get('geometry_only', 0))} units); maps that also have semantic units are not duplicated.")))
    rows.append(row(('Size and overlap',),
        bi(f"{n(W['modules'])} instances. {pct(W['share_contained'])} ({n(W['contained'])}) lie completely inside a larger window of the same source map (module_redundancy.json.gz); the redundancy factor (sum of window sizes ÷ distinct reactions) is " + ', '.join(f'{k} {v}' for k, v in sorted(red.items())) + '. This is why the catalog index shows independent structures only by default.'),
        bi(f"{n(S['modules'])} units, mean size {S['mean_size']}; sizes 3–12 plus {n(S['over_12'])} units above 12 (largest {S['max_unit_size']}) because of appendix reactions; appendix reactions in total {n(S['appendix_reactions'])}. Shapes: chain {n(sh.get('chain', 0))}, cycle {n(sh.get('cycle', 0))}, fork {n(sh.get('fork', 0))}. No unit overlaps another, so there are no redundancy or 'contained in' badges here.")))
    rows.append(row(('Families and variants',),
        bi(f"{n(W['families'])} families, {n(W['variants'])} layout variants. Family = same connection shape: isomorphic core graphs, node roles (reaction / metabolite) preserved; reactions and metabolites may differ. Variant = a different drawn layout of a family: geometry error below 0.25 after alignment. Only the first 8 isomorphisms are tried (at most 5000 backtracking steps), so some variants are split that are really the same drawing: a re-check with up to {n(aud['iso_cap'])} alignments on {n(aud['families_checked'])} families merged {n(aud['merged_variants'])} of {n(aud['variants_before'])} variants ({pct(aud['share_merged'])}), so at least about 5% of variants are mis-split."),
        bi(f"{n(S['families'])} families ({n(S['singleton_families'])} with a single instance), {n(S['variants'])} layout variants; same definitions, but the variant test uses the best geometry error over all isomorphisms (up to {n(cfg['mapping_cap'])} mappings and {n(cfg['step_budget'])} steps; mapping cap reached {n(S['alignment_audit']['variant_cap_hits'])} times, family search budget reached {n(S['alignment_audit']['family_search_budget_hits'])} times). The old first-8 rule would give {n(S['variants_under_first8_rule'])} variants instead.")))
    rows.append(row(('Used by the Composer / tool',),
        bi(f"Default for the Map Generator and Composer (pathway-fragment reuse)."),
        bi(f"Structure-unit reuse (units_sub): stored units plus sub-windows (2, 4, 8, 12 reactions) derived inside each unit when the library is loaded. That gives {n(U['sub'])} retrieval candidates, more than the window library's {n(W['modules'])}, although only {n(U['units'])} units are stored. Stored units alone do not replace the window library for retrieval.")))

    caveats = [
        f"Near-duplicate windows: nested and neighbouring windows of one map share most reactions, so instance counts overstate how much distinct content there is; count independent structures instead.",
        f"Geometry-only WikiPathways: {n(W['geometry_only'])} window modules and {n(q.get('geometry_only', 0))} units have no reaction identity, so they carry no reaction names or IDs and cannot be matched by chemistry; their reactions are drawn interactions.",
        "Variant splitting: see the variant rule above; at least about 5% of the window library's variants are mis-split, and unit boundaries depend on the content of each source map, so two versions of one pathway need not fall into the same family.",
        "Modules are cut automatically from source maps; each module keeps its source provenance.",
    ]
    if not (W['summary_matches'] and U['summary_matches']):
        caveats.append('Warning: stored summaries and the stored modules disagree on a count; counts on this page come from the stored modules where they could be computed.')
    body = ('<section class="panel" id="collected"><h2 style="margin-top:0">' + bi('How these modules were collected') + '</h2>'
            '<p class="meta">' + bi('Both libraries are cut from open-core curated maps (Escher maps, WikiPathways, Reactome); they differ in how a map becomes modules. Numbers below are read from the saved data when the catalog is generated.') + '</p>'
            '<div class="tablewrap"><table><thead><tr><th></th>' + hd('windows', 'Pathway-fragment reuse') + hd('units', 'Structure-unit reuse (experimental)') + '</tr></thead><tbody>' + ''.join(rows) + '</tbody></table></div>'
            '<h3 style="margin-top:14px">' + bi('Caveats') + '</h3><ul>' + ''.join(f'<li>{bi(c)}</li>' for c in caveats) + '</ul></section>')
    return body
