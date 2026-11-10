"""Chinese strings for every static English text node on the generated pages (exact-match dictionary)."""
ZH_EXTRA={
# home
'Current approach and results':'当前方案与结果','Measured drawing improvements':'绘图改进的实测结果','Independent research source (CC BY-NC-SA)':'独立研究来源（CC BY-NC-SA）','Module choice and pose under fixed coverage':'覆盖保持下的模块与姿态选择','Four drawing conditions on E. coli core':'E. coli core 的四种绘图条件','17,221 instances / 6,094 layout variants':'17,221个实例 / 6,094个布局变体','Commands and layout of the code':'命令与代码结构',
'Current research report':'当前研究报告','Drawing-geometry experiment (P0)':'绘图几何实验（P0）','MicroMap extension experiment (P1)':'MicroMap扩展库实验（P1）','Module composition search (P2)':'模块组合搜索（P2）','Improvement and failure cases':'改善与失败案例','Run instructions':'运行说明','Deterministically selected cases':'确定性选取的案例',
'Expert E. coli core map':'专家绘制的 E. coli core 图','Current context-composition result':'当前上下文组合结果','Full E. coli map from the current method':'当前方法生成的完整 E. coli 图',
'Browse':'浏览','Expert maps and current results':'专家图与当前结果','Six-model full-scope results':'六个模型的完整范围结果','Data and verification':'数据与验证','Definitions':'术语定义',
'Model / all views':'模型／全部视图','Module coverage / reactions':'模块覆盖／总反应数','Strict equation coverage':'严格方程覆盖','Data':'数据','Skeleton and mapping':'骨架与映射','Data directory':'数据目录说明','Results directory':'结果目录说明','Automatic verification':'自动验收结果','Paired completion scores':'配对补全逐案例得分',
'Common command:':'常用命令：','Browsing existing results does not require rerunning the experiments.':'查看现有结果无需重跑实验。',
'sources':'个来源','Sources':'来源',
'source map':'来源图',
'Files organised on 2026-11-06. Old versions removed; required source inputs and current results are kept.':'文件整理于2026-11-06。旧版本已删除；保留必要的来源输入和当前结果。',
'Module instance':'模块实例','Family':'家族','Variant':'变体（画法）','Geometry only':'仅几何','Display codes':'显示编号',
'One module (a connected set of reactions) as drawn on one source map.':'一张来源图上画出的一个模块（相连的一组反应）。','All instances with the same connectivity shape: reactions and metabolites can differ; only the connection pattern is the same. Same colour in the catalog.':'连接形状相同的全部实例：反应和代谢物可以不同，只是连接模式相同；在图录中颜色相同。','A different drawn layout of one family; instances with identical layouts share a variant.':'同一家族的不同布局画法；布局相同的实例属于同一变体。','The source drawing has positions but no reaction identifiers, so only the layout is stored.':'来源图只有位置而没有反应标识符，因此只保存布局。',
# catalog
'Search the variants':'检索布局变体','Browse all pages':'浏览全部页面','Layout':'布局缩略图','Module and reactions':'模块与反应','Family · variant':'家族 · 变体','Source':'来源','Quality':'质量','Size':'规模','Instances':'实例数','Licence':'许可','Order':'排序','Reset':'重置','Search':'搜索','Species':'物种','Min size':'最小规模','Max size':'最大规模',
'Downloads:':'下载：','complete JSON exchange (modules.json.gz)':'完整JSON交换文件（modules.json.gz）','variant index (index.json)':'变体索引（index.json）','Overview':'总览','Catalog':'图录','Filter this page':'筛选本页','all':'全部','Source map':'来源图','Module code':'模块编号','Technical details (original IDs)':'技术细节（原始ID）',
'module instances':'模块实例','structural families':'结构家族','layout variants':'布局变体','catalog pages':'图录页数','species':'物种','singleton families':'单例家族',
'Module instances by source':'按来源统计的模块实例','Layout variants by source':'按来源统计的布局变体','Provenance and licence':'来源与许可',
'Group by family':'按家族分组','Size (large first)':'规模（大者在前）','Name A–Z':'名称 A–Z','Show all families':'显示全部家族','What do these terms mean?':'这些术语是什么意思？',
'name, species, family or reaction id':'名称、物种、家族或反应ID','name, species, family or instance id':'名称、物种、家族或实例ID','name, species, family id or reaction id':'名称、物种、家族或反应ID',
'Family = modules with the same connectivity shape (reactions and metabolites can differ; only the connection pattern is the same). Variant = a different drawn layout of that shape.':'家族 = 连接形状相同的模块（反应和代谢物可以不同，只是连接模式相同）；变体（画法）= 该形状的不同布局画法。',
'Same colour = same connectivity shape (family; reactions and metabolites can differ); variants = different layouts of that shape.':'相同颜色 = 相同连接形状（家族，反应和代谢物可以不同）；变体 = 该形状的不同画法。',
'Solid points are observed in the source layout; pale points are derived. Circles = metabolites, lines = reactions (as in the expert maps). Families ignore side attachments; chemical equivalence is a separate index.':'实心点来自来源布局中的观测，浅色点为推导所得。圆点 = 代谢物，线 = 反应（与专家图一致）。家族划分忽略侧边附属物；化学等价另有单独索引。',
'Escher map':'Escher 地图','geometry only, no reaction IDs':'仅几何，无反应ID','Click to show only this family':'点击仅显示该家族','Showing one family only':'仅显示一个家族',
'Variant {a} of {b}':'画法 {a} / {b}','Same connectivity shape: {a} instances · {b} layout variants':'同连接形状 {a} 个实例 · {b} 种画法','{a} of {b} layout variants match':'{a} / {b} 个布局变体匹配',' (first {n} shown)':'（显示前 {n} 个）','{a} of {b} cards on this page shown':'本页 {b} 张卡片中显示 {a} 张',
# composition
'Case':'案例','Chemical coverage':'化学覆盖','Module coverage':'模块覆盖','Label overlaps':'标签重叠','Straight crossings':'直线交叉','Accepted change':'是否采用','yes':'是','no (original kept)':'否（保留原选择）',
'All trials and failures':'全部试验与失败记录','Original selection':'原始选择','Final Escher JSON':'最终 Escher JSON',
'Frozen P1 libraries. Preserve covered reaction sets and equation-compatible sets. Same target-only completion and routing. No learned preference or target expert coordinates.':'冻结P1库；保留已覆盖反应集合与方程相容集合；相同的仅目标信息的补全与布线；不使用学习偏好，也不使用目标专家坐标。',
'Module variant and pose selection under fixed coverage, before → after.':'覆盖保持不变时的模块变体与姿态选择，前 → 后。','Module composition search':'模块组合搜索',
# micromap
'Target':'目标','All reactions':'全部反应','G candidates':'G 候选数','G_D candidates':'G_D 候选数','G_D_M candidates':'G_D_M 候选数','Added chemical candidates':'新增化学候选数','Source module catalog':'来源模块图录','Full saved results':'完整保存结果',
'Licence: CC BY-NC-SA 4.0 (MicroMap master v4.0, DOI 10.7910/DVN/FZKMJ8). Open core unchanged. G: old pool + P0; G_D: same-mapper core; G_D_M: add MicroMap.':'许可：CC BY-NC-SA 4.0（MicroMap master v4.0，DOI 10.7910/DVN/FZKMJ8）。开放核心库未改动。G：旧候选池 + P0；G_D：同映射器核心库；G_D_M：加入MicroMap。',
'One independent master design added as a separate research source; candidate counts per source condition.':'一个独立的总设计作为单独研究来源加入；按来源条件统计候选数。','MicroMap extension experiment':'MicroMap扩展库实验',
'Source DOI':'来源DOI','. Adapted source center-lines and algorithm-derived process centers. CC BY-NC-SA4.0.':'。已改编来源中心线与算法推导的过程中心。CC BY-NC-SA4.0。','MicroMap research layout extension':'MicroMap研究布局扩展库',
'CC BY-NC-SA4.0. One independent master design. Not in the open core. Deterministic first6 IDs per size; observed glyphs and decoded center-line controls, algorithm-derived process centers.':'CC BY-NC-SA4.0。一个独立的总设计，不在开放核心库中。每种规模按ID确定性取前6个；图元为观测所得，中心线控制点为解码所得，过程中心由算法推导。',
# geometry
'node overlap pairs':'节点重叠对数','label overlap pairs':'标签重叠对数','label node intersections':'标签与节点相交数','straight crossings':'直线交叉数','canvas area per reaction':'每反应画布面积','condition':'条件',
'E. coli core, fixed-candidate geometry ablation: four drawing conditions on the same frozen candidates.':'E. coli core 固定候选的几何消融：同一批冻结候选上的四种绘图条件。','Drawing-geometry experiment':'绘图几何实验',
'All 95 reactions; same modules, glyph identities and biochemical participant edges. Same 12 px font. Larger physical scale is reported separately from solver effects.':'全部95个反应；模块、图元身份和生化参与物连线均相同；字体同为12像素。更大的物理尺度与求解器效果分开报告。',
'Escher JSON':'Escher JSON','metrics':'指标',
# models
'reactions accounted':'已核对反应数','strict equation coverage':'严格方程覆盖','reference-compatible':'参考相容反应','geometry/module coverage':'几何／模块覆盖','algorithmic bridges':'算法桥接反应','View':'视图','Reactions':'反应数','Compartments':'区室','Files':'文件','Home':'首页',
# cases
'Full E. coli and the largest pending-fraction model are in ../models/. Local cases select the first ID improving/worsening RMSE; readability is never claimed from RMSE alone.':'完整 E. coli 和待布局比例最大的模型见 ../models/。局部案例取第一个使RMSE改善／变差的ID；绝不仅凭RMSE声称可读性提升。',
'Site':'站点导航','by source':'按来源','improves':'改善','worsens':'变差','First ranking disagreement':'首个排序分歧','complete case':'完整案例',
}
