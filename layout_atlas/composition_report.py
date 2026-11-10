"""Generate the composition-search report from saved trials, including failures."""
import html
from core import *
import pages




def generate(result=None):
    s=result or read(OUT/'composition/summary.json');rows=[]
    for case in s['cases']:
        a=case['baseline'];b=case['best']
        rows.append(f"| {case['model']} / {case['view']} | {case['before_chemical']}→{case['after_chemical']} | {case['before_covered']}→{case['after_covered']} | {a['node_overlap_pairs']}→{b['node_overlap_pairs']} | {a['label_overlap_pairs']}→{b['label_overlap_pairs']} | {a['straight_crossings']}→{b['straight_crossings']} | {a['routing']['metrics']['proper_crossings']}→{b['routing']['metrics']['proper_crossings']} | {case['full_layout_attempts']} |")
    p1=read(OUT/'micromap_extension/summary.json');reference=next(m for m in p1['models'] if m['model']=='e_coli_core');old=next(r for r in reference['drawing'] if r['condition']=='G_D');old_cover=next(r for r in reference['conditions'] if r['condition']=='G_D')['selected_chemical_reactions'];ecoli=next((case for case in s['cases'] if case['model']=='e_coli_core' and case['view']=='full'),None)
    ecoli_note=(f"E. coli core保留{ecoli['after_chemical']}条方程相容模块覆盖（未加入MicroMap时为{old_cover}）。新组合的直线交叉为{ecoli['best']['straight_crossings']}，仍高于未加入MicroMap时的{old['straight_crossings']}；标签重叠也仍高于旧条件。因此这是部分修复，不是彻底消除扩库退化。" if ecoli else '')
    table='\n'.join(rows);improved=sum(c['improved'] for c in s['cases']);failed=sum(t['status']=='failed' for c in s['cases'] for t in c['trials']);attempts=sum(c['full_layout_attempts'] for c in s['cases']);rejected=sum(t['status']=='ok' and not t['acceptable'] for c in s['cases'] for t in c['trials']);unconverged=sum(not c['best']['constraints']['optimizer_converged'] for c in s['cases'])
    text=f'''# P2实测：覆盖保持下的模块变体选择与接口组合

日期：2026-11-06。由`composition_report.py`读取保存的全部试验生成。承接[10扩库实验](10_micromap_extension_results.md)：MicroMap扩大了方程相容候选，但E. coli core替换模板后交叉从124升至204。本轮只改目标上下文中的组合与姿态选择，不增加新专家来源或训练模型。

## 1. 实现与输入

新命令读取P1固定候选池、目标模型、选中模块及其映射，不输入目标专家坐标、主侧标签或专家偏好。三个源库、身份映射、目标化学和P0补全/布线实现保持相同；不改写专家源坐标和模块库。

候选变化包括真实来源的等反应集合画法、90/180/270度旋转、一个模块替换最多两个冲突块，以及小模块的两块拆分。旋转只变换几何，代谢物身份、底物/产物及计量不交换。所有算法姿态记录父候选、旋转矩阵与中心，不计作新增专家实例或专家布局变体。

必须保留原选择的**具体反应集合**和其中方程相容反应集合，不仅维持总数。允许增加覆盖，不允许内部反应重复占用。

## 2. 目标上下文与预算

先以当前骨架为算法坐标框架，估计实际非通用参与物接口的距离、朝向、缺失端口和块体碰撞，再用相同ELK后端预览块体排布。端口外连依据完整目标GEM，固定子图之外没有坐标的区域不冒充已知位置。ATP、水等通用参与物不作为模块邻近性的证据，最终导出仍保留它们。

预算和权重在本轮试验前写入`composition.CONFIG`：预览前最多保留120个候选变化、16次宏布局预览、8次完整布局（包括原选择）。每个已有反应集合最多取4种不同几何代表；拆分子池最多24种反应子集。不同变化类型分层入选预算，避免全部预算落在同一类姿态上。原始变化的枚举量另存于`raw_generated`，120是保留上限而不是所有轻量评分调用的上限。它是有界启发式搜索，不是全局最优覆盖求解。

完整布局用同一物理尺度、图形副本共同求解、12像素字体、标签调整和libavoid路由。候选评分使用实际生成图的交叉、标签冲突、面积和模块形变：`{s['config']['quality_weights']}`。这些指标来自目标网络生成图，不是专家答案。

新结果若比原选择增加节点近接、标签重叠、直线交叉或正交交叉，就不部署；面积不能超过原选择1.5倍，形变也受限。没有满足条件的改进就保留原选择。它不能保证其他所有指标都单调改善，边段穿越标签、形变、面积和时间仍逐试验保存。

## 3. 固定范围结果

范围沿用P1：完整95反应E. coli core，以及iML1515/iYO844原有首、中、末三个固定视图。本轮不重新挑选效果好的范围。每个原选择先重放，与P1主要几何指标一致后才搜索。

| 范围 | 方程相容覆盖 | 模块覆盖 | 节点近接 | 标签重叠 | 直线交叉 | 正交交叉 | 完整布局次数 |
|---|---:|---:|---:|---:|---:|---:|---:|
{table}

{ecoli_note}

{len(s['cases'])}个固定范围中，{improved}个接受了新的组合/姿态；共尝试{attempts}次完整布局，执行失败{failed}次，另有{rejected}次因防退化规则拒绝部署。未改善、预算未覆盖及原选择回退均保留，不把规则回退计作学习成功。

上表覆盖以实际绘制范围为分母；大模型子图不代表全部GEM已得到新排布。它们来自三个目标模型，且每模型的窗口相互依赖，不能以试验次数或窗口数量扩大独立样本数。

## 4. 如何理解证据

这是目标指标驱动的组合搜索。由于生成过程中已经使用这些指标选择结果，指标改善本身不是独立的人类可读性验证，也不是学习模型的泛化成绩。权重未用测试专家坐标调节；后续应在新来源/目标上检验成本和效果。

不同来源代表和合并/拆分可能改变图形副本数量；所有反应、计量与目标参与物仍完整，节点数在逐案例指标中另列。因此该对照衡量端到端组合，而非单独归因于补全求解器。保留所有候选的失败与退化指标，最终优胜者不是只挑最好图片后删除其他结果。

最终采用结果中{unconverged}/{len(s['cases'])}没有达到优化器的收敛标志，虽然已通过本轮锚点/间距及化学检查，仍只能称为固定350步预算下的可行结果，不能声称已收敛。

算法仍受候选库缺口、预览近似、预算和现有ELK端口实现限制；本轮没有实现完整真实膜边界或全GEM最优组合。没有训练新的排序器或EGNN；数据库复用、组合改进和学习增益分开判断。

## 5. 重运行与文件

```sh
python3 project.py compose
python3 project.py compose --model e_coli_core
python3 layout_atlas/composition_replay.py
python3 project.py verify
```

- [全部案例与试验入口](../layout_atlas/results/composition/index.html)。
- [完整汇总](../layout_atlas/results/composition/summary.json)。
- 每个范围的`seed/`保留原选择重放，`trials/`保存完整布局试验，`best/`保存实际采用的结果；源许可与算法姿态一并保留。
- `preview_log.json`保存候选保留与宏布局代价；单案例`summary.json`保存覆盖集合、逐试验质量、接受/拒绝标记及失败原因；根汇总固定代码、P1输入版本、权重和随机种子。

下一步按失败记录推进：若只旋转就有效，应先改善接口的显式方向约束；若没有可用等价画法，则补充独立设计或改善多尺度覆盖；若预览与完整布局排序不一致，则改进上下文代价。学习器应在这套真实候选与组合任务上另做来源隔离评价。
'''
    path=ROOT/'Proposal/11_context_composition_results.md';path.write_text(text)
    (OUT/'composition/index.html').write_text(pages.composition_index(s))
    return path

if __name__=='__main__':print(generate())
