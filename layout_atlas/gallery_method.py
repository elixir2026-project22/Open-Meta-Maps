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
    items = [('windows', ('Pathway-fragment reuse', '通路片段复用'), other_href_windows), ('units', ('Structure-unit reuse (experimental)', '结构单元复用（实验性）'), other_href_units)]
    parts = []
    for key, (en, zh), href in items:
        if key == current:
            parts.append(f'<span aria-current="page" title="current library">{bi(en, zh)}</span>')
        else:
            parts.append(f'<a href="{href}">{bi(en, zh)}</a>')
    return ('<div class="pager libsw" role="navigation" aria-label="Library"><span style="background:none;border:0;color:var(--mute)">' + bi('Library:', '模块库：') + '</span>' + ''.join(parts) + '</div>')


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

    def hd(key, en, zh):
        mark = ('<span class="tag ind">' + bi('this page', '本页') + '</span>') if key == current else ''
        return f'<th>{bi(en, zh)} {mark}</th>'

    def row(label, w, u):
        return f'<tr><th scope="row">{bi(*label)}</th><td>{w}</td><td>{u}</td></tr>'

    rows = []
    rows.append(row(('Source maps and licences', '来源图与许可'),
        bi(f"{n(W['maps'])} source layouts used (of {n(W['registered'])} registered map records): Escher {n(mb.get('escher', 0))}, WikiPathways {n(mb.get('wikipathways', 0))}, Reactome {n(mb.get('reactome', 0))}. Licence of the source maps: " + ', '.join(f'{k} {n(v)}' for k, v in sorted(ml.items())) + '. Each module keeps the licence, SHA-256 and source pointer of its map.',
           f"使用了 {n(W['maps'])} 张来源布局（已登记的地图记录共 {n(W['registered'])} 条）：Escher 社区 {n(mb.get('escher', 0))}、WikiPathways {n(mb.get('wikipathways', 0))}、Reactome {n(mb.get('reactome', 0))}。来源图许可：" + '、'.join(f'{k} {n(v)}' for k, v in sorted(ml.items())) + '。每个模块都保留其来源图的许可、SHA-256 和来源指针。'),
        bi(f"{n(s['source_maps'])} source maps recut: Escher {n(by_src.get('escher', 0))}, WikiPathways {n(by_src.get('wikipathways', 0))}, Reactome {n(by_src.get('reactome', 0))}; licences " + ', '.join(f'{k} {n(v)}' for k, v in sorted(lic_u.items())) + f". Only {n(S['source_layouts_with_units'])} of them yield at least one unit. The map set is not identical to the window library's (counts differ; the overlap itself was not checked). Open-core data only; the extension is flagged admitted_to_open_core = {str(s['admitted_to_open_core']).lower()}.",
           f"重新切分了 {n(s['source_maps'])} 张来源图：Escher 社区 {n(by_src.get('escher', 0))}、WikiPathways {n(by_src.get('wikipathways', 0))}、Reactome {n(by_src.get('reactome', 0))}；许可：" + '、'.join(f'{k} {n(v)}' for k, v in sorted(lic_u.items())) + f"。其中只有 {n(S['source_layouts_with_units'])} 张能切出至少一个单元。两个库用的地图集合并不完全相同（数量不同；两者交集本身未核对）。只用开放核心数据；该扩展库的标记为 admitted_to_open_core = {str(s['admitted_to_open_core']).lower()}（尚未并入开放核心）。")))
    rows.append(row(('What a module instance is', '模块实例是什么'),
        bi('A connected set of reactions of one source map with their core metabolites, stored with source coordinates, reaction identities (when the map has them), licence and source pointer. Different instances may contain the same reactions.',
           '同一张来源图中相连的一组反应及其核心代谢物，连同来源坐标、反应标识（若该图有）、许可和来源指针一起保存。不同实例可以包含相同的反应。'),
        bi(f"A connected set of {cfg['min_size']}–{cfg['max_size']} core reactions of one source map (appendix reactions can raise this to at most {cfg['max_appendix_total']}) forming one drawing shape: chain, cycle or fork. Currency metabolites (water, protons, ATP/ADP/AMP, NAD(P)(H), CoA, phosphate) are not part of the core. Reactions in two units of the same map: {n(S['overlap_check_reactions_in_multiple_units'])}.",
           f"同一张来源图中相连的 {cfg['min_size']}–{cfg['max_size']} 个核心反应（附属反应可使总数最多到 {cfg['max_appendix_total']}），构成一个绘图形状：链、环或叉。货币代谢物（水、质子、ATP/ADP/AMP、NAD(P)(H)、辅酶A、磷酸）不属于核心。同一张图内出现在两个单元里的反应数：{n(S['overlap_check_reactions_in_multiple_units'])}。")))
    rows.append(row(('How it was cut', '切分方式'),
        bi(f"Overlapping breadth-first windows. WikiPathways maps whose interactions are bound to Rhea reactions and Reactome diagrams are cut by database.diagram_modules: evenly spaced root reactions per map (every k-th reaction with k = N÷24 rounded down and at least 1, N = reactions of the map: all reactions for small maps, roughly 24–47 roots for large ones); from each root a breadth-first walk over the reaction–metabolite graph takes the first 1, 2, 4, 8, 12 and 24 reactions, and every distinct reaction set is stored as a window. Escher maps and WikiPathways maps without reaction identity come from stored source-observation records (module_observations.json.gz), kept when they have at most 100 nodes and 1–40 reactions (geometry-only: drawn interactions); the program that first cut those records is not part of this code base (not verified). {pct(W['share_fixed'])} of all windows have exactly 1, 2, 4, 8, 12 or 24 reactions; the sizes found are {W['min_size']}–{W['max_size']}.",
           f"重叠的广度优先窗口。WikiPathways 中交互已绑定 Rhea 反应的图和 Reactome 图由 database.diagram_modules 切分：每张图取等距的起点反应（每隔 k 个取一个，k = N÷24 向下取整且至少为 1，N 为该图的反应数：小图取全部反应，大图约 24–47 个起点）；从每个起点沿反应—代谢物图做广度优先遍历，取前 1、2、4、8、12、24 个反应，每个不同的反应集合保存为一个窗口。Escher 社区图和无反应标识的 WikiPathways 图来自已保存的来源观测记录（module_observations.json.gz），保留不超过 100 个节点、含 1–40 个反应的记录（仅几何：按画出的交互计）；最初切分这些记录的程序不在本代码库中（未核实）。全部窗口中有 {pct(W['share_fixed'])} 恰好含 1、2、4、8、12 或 24 个反应；实际规模范围为 {W['min_size']}–{W['max_size']}。"),
        bi(f"Non-overlapping, cut once per source map by units_library.decompose, in fixed order: (1) cycles: the shortest cycle through every edge, at most {cfg['cycle_max_reactions']} reactions, taken shortest first; ring rule: rings of at least {cfg['ring_min_reactions']} reactions without transport reactions are seeded first, and a neighbouring ring is fused when it shares at least {cfg['fuse_min_shared']} reactions and adds at most {cfg['fuse_max_new']}; (2) forks: a branch point with up to {cfg['fork_max_arms']} arms, arm depth up to {cfg['fork_arm_depth']}, within {cfg['max_size']} reactions; (3) chains: what remains is cut at branch points, and paths longer than {cfg['max_size']} are halved recursively; (4) fragments below {cfg['min_size']} reactions are attached as appendix to an adjacent unit (total at most {cfg['max_appendix_total']}); isolated islands that stay below {cfg['min_size']} are orphans and not modules ({n(S['orphan_reactions'])} orphan reactions; {n(S['unit_reactions'])} of {n(S['source_core_reactions'])} core reactions are inside units). WikiPathways maps without reaction identity are segmented by drawn shape only: each drawn interaction between two metabolite glyphs is a pseudo-reaction (quality = geometry_only; {n(q.get('geometry_only', 0))} units); maps that also have semantic units are not duplicated.",
           f"不重叠，每张来源图只切一次，由 units_library.decompose 按固定顺序执行：（1）环：经过每条边的最短环，不超过 {cfg['cycle_max_reactions']} 个反应，从短到长选取；环规则：先以不含转运反应、至少 {cfg['ring_min_reactions']} 个反应的环为种子，相邻的环若共享至少 {cfg['fuse_min_shared']} 个反应且新增不超过 {cfg['fuse_max_new']} 个反应就融合；（2）叉：一个分支点，最多 {cfg['fork_max_arms']} 条臂，臂深最多 {cfg['fork_arm_depth']}，总数不超过 {cfg['max_size']} 个反应；（3）链：剩余部分在分支点处切开，超过 {cfg['max_size']} 个反应的路径递归对半切；（4）不足 {cfg['min_size']} 个反应的碎片作为附属并入相邻单元（总数不超过 {cfg['max_appendix_total']}）；仍不足 {cfg['min_size']} 个的孤立小岛是孤岛，不算模块（孤岛反应 {n(S['orphan_reactions'])} 个；{n(S['source_core_reactions'])} 个核心反应中有 {n(S['unit_reactions'])} 个在单元内）。没有反应标识的 WikiPathways 图只按画出的形状切分：两个代谢物图元之间每条画出的交互视为一个伪反应（quality = geometry_only；共 {n(q.get('geometry_only', 0))} 个单元）；已有语义单元的图不重复切分。")))
    rows.append(row(('Size and overlap', '规模与重叠'),
        bi(f"{n(W['modules'])} instances. {pct(W['share_contained'])} ({n(W['contained'])}) lie completely inside a larger window of the same source map (module_redundancy.json.gz); the redundancy factor (sum of window sizes ÷ distinct reactions) is " + ', '.join(f'{k} {v}' for k, v in sorted(red.items())) + '. This is why the catalog index shows independent structures only by default.',
           f"{n(W['modules'])} 个实例。其中 {pct(W['share_contained'])}（{n(W['contained'])} 个）完全落在同一来源图中更大的窗口之内（module_redundancy.json.gz）；冗余系数（窗口规模之和 ÷ 不同反应数）为 " + '、'.join(f'{k} {v}' for k, v in sorted(red.items())) + '。因此图录索引默认只显示独立结构。'),
        bi(f"{n(S['modules'])} units, mean size {S['mean_size']}; sizes 3–12 plus {n(S['over_12'])} units above 12 (largest {S['max_unit_size']}) because of appendix reactions; appendix reactions in total {n(S['appendix_reactions'])}. Shapes: chain {n(sh.get('chain', 0))}, cycle {n(sh.get('cycle', 0))}, fork {n(sh.get('fork', 0))}. No unit overlaps another, so there are no redundancy or 'contained in' badges here.",
           f"{n(S['modules'])} 个单元，平均规模 {S['mean_size']}；规模 3–12，另有 {n(S['over_12'])} 个因附属反应超过 12（最大 {S['max_unit_size']}）；附属反应共 {n(S['appendix_reactions'])} 个。形状：链 {n(sh.get('chain', 0))}、环 {n(sh.get('cycle', 0))}、叉 {n(sh.get('fork', 0))}。单元互不重叠，所以这里没有冗余或“包含于”标记。")))
    rows.append(row(('Families and variants', '家族与变体'),
        bi(f"{n(W['families'])} families, {n(W['variants'])} layout variants. Family = same connection shape: isomorphic core graphs, node roles (reaction / metabolite) preserved; reactions and metabolites may differ. Variant = a different drawn layout of a family: geometry error below 0.25 after alignment. Only the first 8 isomorphisms are tried (at most 5000 backtracking steps), so some variants are split that are really the same drawing: a re-check with up to {n(aud['iso_cap'])} alignments on {n(aud['families_checked'])} families merged {n(aud['merged_variants'])} of {n(aud['variants_before'])} variants ({pct(aud['share_merged'])}), so at least about 5% of variants are mis-split.",
           f"{n(W['families'])} 个家族、{n(W['variants'])} 个布局变体。家族 = 连接形状相同：核心图同构，保持节点角色（反应／代谢物）；反应和代谢物可以不同。变体 = 家族的不同画法：对齐后几何误差小于 0.25。只尝试前 8 个同构（回溯至多 5000 步），所以有些本来是同一画法的变体被拆开了：对 {n(aud['families_checked'])} 个家族用至多 {n(aud['iso_cap'])} 次对齐重新检查，{n(aud['variants_before'])} 个变体中合并了 {n(aud['merged_variants'])} 个（{pct(aud['share_merged'])}），因此至少约 5% 的变体被错误拆分。"),
        bi(f"{n(S['families'])} families ({n(S['singleton_families'])} with a single instance), {n(S['variants'])} layout variants; same definitions, but the variant test uses the best geometry error over all isomorphisms (up to {n(cfg['mapping_cap'])} mappings and {n(cfg['step_budget'])} steps; mapping cap reached {n(S['alignment_audit']['variant_cap_hits'])} times, family search budget reached {n(S['alignment_audit']['family_search_budget_hits'])} times). The old first-8 rule would give {n(S['variants_under_first8_rule'])} variants instead.",
           f"{n(S['families'])} 个家族（其中 {n(S['singleton_families'])} 个只有一个实例）、{n(S['variants'])} 个布局变体；定义相同，但变体判定取全部同构中的最小几何误差（至多 {n(cfg['mapping_cap'])} 个映射、{n(cfg['step_budget'])} 步；映射上限触顶 {n(S['alignment_audit']['variant_cap_hits'])} 次，家族搜索预算触顶 {n(S['alignment_audit']['family_search_budget_hits'])} 次）。旧的“前 8 个同构”规则会得到 {n(S['variants_under_first8_rule'])} 个变体。")))
    rows.append(row(('Used by the Composer / tool', '组合器／工具如何使用'),
        bi(f"Default (option windows).", '默认使用（选项“标准窗口”）。'),
        bi(f"Only through the option units_sub = stored units + sub-windows (2, 4, 8, 12 reactions) derived inside each unit when the library is loaded. That gives {n(U['sub'])} retrieval candidates, more than the window library's {n(W['modules'])}, although only {n(U['units'])} units are stored. Proposal/13_units_and_tool.md concludes that the stored units alone do not replace the window library for retrieval (not re-tested here).",
           f"只能通过选项 units_sub = 已存单元 + 加载时在每个单元内部派生的子窗口（2、4、8、12 个反应）。这样得到 {n(U['sub'])} 个检索候选，比窗口库的 {n(W['modules'])} 个更多，尽管实际只存了 {n(U['units'])} 个单元。Proposal/13_units_and_tool.md 的结论是：单独的存储单元不能替代窗口库用于检索（本页未重新测试）。")))

    caveats = [
        (f"Near-duplicate windows: nested and neighbouring windows of one map share most reactions, so instance counts overstate how much distinct content there is; count independent structures instead.",
         '近似重复的窗口：同一张图中嵌套或相邻的窗口共享大部分反应，所以实例数会夸大不同内容的多少；应改数独立结构。'),
        (f"Geometry-only WikiPathways: {n(W['geometry_only'])} window modules and {n(q.get('geometry_only', 0))} units have no reaction identity, so they carry no reaction names or IDs and cannot be matched by chemistry; their reactions are drawn interactions.",
         f"仅几何的 WikiPathways：{n(W['geometry_only'])} 个窗口模块和 {n(q.get('geometry_only', 0))} 个单元没有反应标识，所以没有反应名或ID，也不能按化学匹配；它们的“反应”只是画出的交互。"),
        ("Variant splitting: see the variant rule above; at least about 5% of the window library's variants are mis-split, and unit boundaries depend on the content of each source map, so two versions of one pathway need not fall into the same family.",
         '变体拆分：见上面的变体规则；窗口库至少约 5% 的变体被错误拆分；单元边界取决于每张来源图的内容，同一通路的两个版本不一定落入同一个家族。'),
        ("Independence: observed source layouts do not prove independent designers; automatic checks only, no human certification is claimed.",
         '独立性：观测到的来源布局不能证明是独立设计者；只有自动检查，不声称经过人工认证。'),
    ]
    if not (W['summary_matches'] and U['summary_matches']):
        caveats.append(('Warning: stored summaries and the stored modules disagree on a count; counts on this page come from the stored modules where they could be computed.', '警告：已保存的汇总与已存模块在某个计数上不一致；本页的计数在可计算处取自已存模块。'))
    body = ('<section class="panel" id="collected"><h2 style="margin-top:0">' + bi('How these modules were collected', '这些模块是怎么收集的') + '</h2>'
            '<p class="meta">' + bi('Both libraries are cut from open-core expert maps (Escher maps, WikiPathways, Reactome); they differ in how a map becomes modules. Numbers below are read from the saved data when the catalog is generated.',
                                    '两个库都切自开放核心的专家图（Escher 社区图、WikiPathways、Reactome），区别在于一张图如何变成模块。下面的数字在生成图录时从已保存的数据中读取。') + '</p>'
            '<div class="tablewrap"><table><thead><tr><th></th>' + hd('windows', 'Pathway-fragment reuse', '通路片段复用') + hd('units', 'Structure-unit reuse (experimental)', '结构单元复用（实验性）') + '</tr></thead><tbody>' + ''.join(rows) + '</tbody></table></div>'
            '<h3 style="margin-top:14px">' + bi('Caveats', '需要注意') + '</h3><ul>' + ''.join(f'<li>{bi(a, b)}</li>' for a, b in caveats) + '</ul></section>')
    return body
