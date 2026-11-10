"""Light-weight HTML index pages built from saved results (no heavy imports)."""
import html
from core import *
import ui


def composition_index(s):
    e=html.escape;rows=[];figs=''
    for case in s['cases']:
        a=case['baseline'];b=case['best'];prefix=case['model']+'/'+case['view']
        rows.append([f'<a href="#{e(case["model"])}-{e(case["view"])}">{e(prefix)}</a>',f"{case['before_chemical']} → {case['after_chemical']}",f"{case['before_covered']} → {case['after_covered']}",f"{a['label_overlap_pairs']} → {b['label_overlap_pairs']}",f"{a['straight_crossings']} → {b['straight_crossings']}",'yes' if case['improved'] else 'no (original kept)'])
        figs+=f'<article id="{e(case["model"])}-{e(case["view"])}"><h2>{e(prefix)}</h2><p><a href="{prefix}/summary.json">All trials and failures</a> · <a href="{prefix}/seed/map.svg">Original selection</a> · <a href="{prefix}/best/map.json">Final Escher JSON</a></p><a href="{prefix}/best/map.svg"><img src="{prefix}/best/map.svg" alt="{e(prefix)}" loading="lazy"></a></article>'
    body='<div class="note">Frozen P1 libraries. Preserve covered reaction sets and equation-compatible sets. Same target-only completion and routing. No learned preference or target expert coordinates.</div>'+ui.table(['Case','Chemical coverage','Module coverage','Label overlaps','Straight crossings','Accepted change'],rows,raw=True)+figs
    return ui.page(OUT/'composition/index.html','Module composition search',body,section='composition',lead='<span class="stage">P2</span>Module variant and pose selection under fixed coverage, before → after.')


def micromap_index(models):
    e=html.escape;rows=[];figs=''
    for m in models:
        r={x['condition']:x for x in m['conditions']};rows.append([m['model'],r['G']['target_reactions'],*[r[k]['candidate_chemical_reactions'] for k in ('G','G_D','G_D_M')],len(m['new_candidate_chemical_ids'])])
        for r in m.get('drawing',[]):figs+=f'<article><h2>{e(m["model"])} / {e(r["condition"])}</h2><a href="models/{m["model"]}/{r["condition"]}/map.svg"><img src="models/{m["model"]}/{r["condition"]}/map.svg" alt="{e(m["model"])} {e(r["condition"])}" loading="lazy"></a></article>'
        if m.get('fixed_views'):
            figs+=f'<h2><span class="dn">{e(m["model"])}</span> '+ui.bi('fixed views','固定视图')+'</h2><ul>'+''.join(f'<li><a href="models/{m["model"]}/{r["condition"]}/{r["view"]}/map.svg">{e(r["condition"])} / {e(r["view"])} (</a>'+ui.bi(f'{r["reactions"]} reactions',f'{r["reactions"]} 个反应')+')</li>' for r in m['fixed_views'])+'</ul>'
    body=('<div class="note">Licence: CC BY-NC-SA 4.0 (MicroMap master v4.0, DOI 10.7910/DVN/FZKMJ8). Open core unchanged. G: old pool + P0; G_D: same-mapper core; G_D_M: add MicroMap.</div>'
        '<p><a href="catalog/index.html">Source module catalog</a> · <a href="summary.json">Full saved results</a></p>'
        +ui.table(['Target','All reactions','G candidates','G_D candidates','G_D_M candidates','Added chemical candidates'],rows)+figs)
    return ui.page(OUT/'micromap_extension/index.html','MicroMap extension experiment',body,section='micromap',lead='<span class="stage">P1</span>One independent master design added as a separate research source; candidate counts per source condition.')


def geometry_index(results,destination):
    cols=['condition','node_overlap_pairs','label_overlap_pairs','label_node_intersections','straight_crossings','canvas_area_per_reaction']
    rows=[[round(r[k],2) if isinstance(r[k],float) else r[k] for k in cols] for r in results]
    figs=''.join(f'<article><h2>{html.escape(r["condition"])}</h2><a href="{r["condition"]}/map.svg"><img src="{r["condition"]}/map.svg" alt="{html.escape(r["condition"])}" loading="lazy"></a><p><a href="{r["condition"]}/map.json">Escher JSON</a> · <a href="{r["condition"]}/metrics.json">metrics</a></p></article>' for r in results)
    body='<div class="note">All 95 reactions; same modules, glyph identities and biochemical participant edges. Same 12 px font. Larger physical scale is reported separately from solver effects.</div>'+ui.table([c.replace('_',' ') for c in cols],rows)+figs
    return ui.page(destination/'index.html','Drawing-geometry experiment',body,section='geometry',lead='<span class="stage">P0</span>E. coli core, fixed-candidate geometry ablation: four drawing conditions on the same frozen candidates.')


def view_index(folder,name,condition,total,strict,refs,used,bridges,rows):
 e=html.escape;stats='<div class="stats">'+ui.stat(total,'reactions accounted')+ui.stat(strict,'strict equation coverage')+ui.stat(refs,'reference-compatible')+ui.stat(used,'geometry/module coverage')+ui.stat(bridges,'algorithmic bridges')+'</div>'
 lead=ui.bi(f'All {total} reactions accounted. Strict equation coverage {strict}, reference-compatible {refs}, geometry/module coverage {used}, algorithmic bridges {bridges}. These are application results.',f'全部 {total} 个反应均已核对。严格方程覆盖 {strict}，参考相容 {refs}，几何／模块覆盖 {used}，算法桥接 {bridges}。以上为应用结果。')
 table=ui.table(['View','Reactions','Compartments','Files'],[[e(v['id']),len(v['reaction_ids']),e(', '.join(map(str,v['compartments']))),f'<a href="{v["id"]}.svg">SVG</a> · <a href="{v["id"]}.json">Escher JSON</a>'] for v in rows],raw=True)
 return ui.page(folder/'index.html',f'{name} / {condition}',stats+table,lead=lead,crumbs=[('Home',ui.rel(folder/'index.html','index.html')),(name,None)])


def cases_page(dest,cards):
 body=''.join(f'<article>{c}</article>' for c in cards)+'<p><a href="first_ranking_disagreement.json">First ranking disagreement</a></p>'
 return ui.page(dest/'index.html','Deterministically selected cases',body,section='cases',lead='Full E. coli and the largest pending-fraction model are in ../models/. Local cases select the first ID improving/worsening RMSE; readability is never claimed from RMSE alone.')
