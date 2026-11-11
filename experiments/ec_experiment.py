"""Evaluate the EC channel: reactions it adds beyond complete-equation matching, and how often it agrees with it.

Retrieval-level only (no routing/learning): candidates are generated from the fixed libraries and counted.
  python3 ec_experiment.py [model ...]
"""

import sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(__file__).resolve().parents[1]
if str(_ROOT / 'layout_atlas') not in _sys.path:
    _sys.path.insert(0, str(_ROOT / 'layout_atlas'))
if str(_ROOT) not in _sys.path:
    _sys.path.insert(0, str(_ROOT))

import gzip
import sys
import time

from core import *
import database
import semantic_transfer as st
import ec_channel as ec
import brenda_library as bl

RESULTS = OUT / 'ec_channel'
# Same six targets as MODELS (copied so this retrieval-level study does not import the drawing stack).
MODELS = [('e_coli_core', EXAMPLES / 'e_coli_core.xml', 'Escherichia coli'), ('iML1515', EXAMPLES / 'iML1515.xml', 'Escherichia coli'),
          ('iMM904', EXAMPLES / 'iMM904.xml', 'Saccharomyces cerevisiae'), ('Human-GEM', HUMAN_GEM, 'Homo sapiens'),
          ('iYO844', DATA / 'models/iYO844.xml', 'Bacillus subtilis'), ('AraCore_v2_1', DATA / 'models/AraCore_v2_1.xml', 'Arabidopsis thaliana')]


def agreement(strong, other):
    """Compare EC proposals with complete-equation proposals on the same source module."""
    by_module = defaultdict(set)
    for p in strong:
        by_module[p['template_id']].add(frozenset(p['reaction_ids']))
    both = [p for p in other if p['template_id'] in by_module]
    equal = sum(frozenset(p['reaction_ids']) in by_module[p['template_id']] for p in both)
    ec_sets = defaultdict(set)
    for p in other:
        ec_sets[p['template_id']].add(frozenset(p['reaction_ids']))
    strong_sets = sum(len(v) for k, v in by_module.items() if k in ec_sets)
    recovered = sum(len(v & ec_sets[k]) for k, v in by_module.items() if k in ec_sets)
    return dict(modules_with_both=len({p['template_id'] for p in both}), ec_proposals_on_those_modules=len(both),
                ec_proposals_equal_to_equation_match=equal,
                precision_vs_equation=equal / len(both) if both else None,
                equation_sets_on_those_modules=strong_sets, equation_sets_recovered_by_ec=recovered,
                recall_vs_equation=recovered / strong_sets if strong_sets else None)


def examples(proposals, modules, k=12):
    by_id = {m['id']: m for m in modules}
    rows = []
    for p in sorted(proposals, key=lambda p: (-len(p['reaction_ids']), p['id']))[:k]:
        m = by_id[p['template_id']]
        rows.append(dict(source=p['source'], module=m['name'], size=len(p['reaction_ids']), target_reactions=p['reaction_ids'],
                         shared_ec=[c['shared_ec'] for c in p['mapping_evidence']], unmapped_source_nodes=len(p['unmapped_source_nodes'])))
    return rows


def select_conditions(name, model, allids, strong, core_ec, br_ec, core, brenda):
    """Selection conditions (retrieval-level; selections are written in the format the drawing stage consumes).

    G_D     existing pool (saved retrieval candidates + complete-equation core candidates) -- the MicroMap study's G_D
    G_D_E   + EC candidates from the open-core library, gap-filling only
    G_D_E_B + EC candidates from the BRENDA extension (layout origin unverified), gap-filling only
    """
    old_path = OUT / 'models' / name / 'all_geometry' / 'candidates.json.gz'
    with gzip.open(old_path, 'rt') as f:
        old = json.load(f)
    templates = {m['id']: m for m in core + brenda}
    for p in strong + core_ec + br_ec:
        t = templates[p['template_id']]
        for key in ('license', 'source_sha256', 'original_instance', 'source_citation'):
            if key in t:
                p[key] = t[key]
        p.setdefault('partition', t.get('partition', 'open_core'))
    base = old + strong
    folder = RESULTS / 'models' / name
    folder.mkdir(parents=True, exist_ok=True)
    rows = {}
    for condition, pools in (('G_D', []), ('G_D_E', [core_ec]), ('G_D_E_B', [core_ec, br_ec])):
        selected, stages = ec.select_with_ec(base, pools)
        covered = {r for p in selected for r in p['reaction_ids']}
        chem_cov = {r for p in selected if p.get('complete_chemical_equation_verified', p.get('complete_equation_verified')) for r in p['reaction_ids']}
        ec_sel = [p for p in selected if p['match_level'] == 'ec_topology']
        rows[condition] = dict(selected_modules=len(selected), selected_reactions=len(covered), selected_chemical_reactions=len(chem_cov),
                               ec_selected_modules=len(ec_sel), ec_selected_reactions=len({r for p in ec_sel for r in p['reaction_ids']}),
                               ec_selected_from_brenda=len([p for p in ec_sel if p['source'] == 'brenda']),
                               ec_selected_sizes=dict(sorted(Counter(len(p['reaction_ids']) for p in ec_sel).items())),
                               selected_reaction_share=len(covered) / len(allids), stage_modules=[len(s) for s in stages],
                               match_levels=dict(Counter(p['match_level'] for p in selected)))
        with gzip.open(folder / (condition + '-selected.json.gz'), 'wt') as f:
            json.dump(selected, f, ensure_ascii=False)
    rows['gain_E_over_baseline'] = rows['G_D_E']['selected_reactions'] - rows['G_D']['selected_reactions']
    rows['gain_B_over_E'] = rows['G_D_E_B']['selected_reactions'] - rows['G_D_E']['selected_reactions']
    rows['base_pool'] = dict(saved_retrieval_candidates=len(old), complete_equation_candidates=len(strong))
    return rows


def run(names=None):
    RESULTS.mkdir(parents=True, exist_ok=True)
    core = database.modules()
    brenda = bl.load()
    enrich_core, enrich_brenda = ec.enrich(core), ec.enrich(brenda)
    bridge = st.annotation_bridge(MODELS)
    summary = dict(source_ec_coverage=dict(core=enrich_core, brenda=enrich_brenda), models=[],
                   code_sha256={p: sha((BASE / p).read_bytes()) for p in ('ec_channel.py', 'ec_experiment.py', 'brenda_library.py')},
                   config=ec.CONFIG, brenda_status='unverified layout origin; research extension; not counted as expert designs')
    for name, path, species in MODELS:
        if names and name not in names:
            continue
        t0 = time.time()
        model = load_model(path)
        rids = {r['id'] for r in model['reactions']}
        tec = ec.target_ecs(model)
        strong, sinfo = st.Mapper(model, bridge).proposals(core, 'core complete-equation')
        core_ec, cinfo = ec.ECMapper(model, bridge).proposals(core, 'core EC')
        br_ec, binfo = ec.ECMapper(model, bridge).proposals(brenda, 'BRENDA EC')
        S, C, B = set(sinfo['candidate_reaction_ids']), set(cinfo['candidate_reaction_ids']), set(binfo['candidate_reaction_ids'])
        selection = select_conditions(name, model, rids, strong, core_ec, br_ec, core, brenda)
        row = dict(model=name, reactions=len(rids), target_reactions_with_full_ec=sum(bool(v) for v in tec.values()),
                   complete_equation=dict(proposals=sinfo['candidates'], reactions=len(S)),
                   core_ec=dict(proposals=cinfo['candidates'], reactions=len(C), skips=cinfo['skips'], budget_reached=cinfo['budget_reached']),
                   brenda_ec=dict(proposals=binfo['candidates'], reactions=len(B), skips=binfo['skips'], budget_reached=binfo['budget_reached']),
                   added_by_core_ec=len(C - S), added_by_brenda_ec=len(B - S), added_by_either=len((C | B) - S),
                   union_all=len(S | C | B), union_equation_only=len(S),
                   agreement_core_ec=agreement(strong, core_ec),
                   proposal_levels=dict(equation=sinfo['candidates'], ec=cinfo['candidates'] + binfo['candidates']),
                   examples_core_ec_only=examples([p for p in core_ec if not set(p['reaction_ids']) <= S], core),
                   examples_brenda_ec_only=examples([p for p in br_ec if not set(p['reaction_ids']) <= S], brenda),
                   selection=selection, seconds=round(time.time() - t0, 1))
        summary['models'].append(row)
        dump(RESULTS / 'summary.json', summary)
        print(name, {c: (row['selection'][c]['selected_reactions'], row['selection'][c]['ec_selected_reactions']) for c in ('G_D', 'G_D_E', 'G_D_E_B')}, flush=True)
        print(name, {k: row[k] for k in ('reactions', 'target_reactions_with_full_ec', 'union_equation_only', 'added_by_core_ec', 'added_by_brenda_ec', 'union_all')},
              row['agreement_core_ec'], row['seconds'], flush=True)
    return summary


if __name__ == '__main__':
    run(sys.argv[1:] or None)
