"""Draw the EC-channel selections with the same P0-compatible pipeline as the MicroMap study.

Needs the drawing stack (PIL, jsonschema, node with elkjs/fcose, libavoid for this architecture).
  PATH=$HOME/.local/node/bin:$PATH python3 ec_drawing.py [model ...]

Views (fixed before drawing, same rule as micromap_experiment): the first, middle and last 50-reaction views of
results/models/<model>/all_geometry/views.json, plus the view with the most EC-selected reactions (chosen by EC
presence in the G_D_E selection, not by how the drawing looks).  Conditions: G_D, G_D_E, G_D_E_B (see ec_experiment).
"""
import gzip
import sys
import time

from core import *
import database
import drawing
import micromap_experiment as me
import routing_arch

RESULTS = OUT / 'ec_channel'
CONDITIONS = ['G_D', 'G_D_E', 'G_D_E_B']
MODELS = {'iML1515': REF / 'models/iML1515.xml', 'iMM904': REF / 'models/iMM904.xml', 'iYO844': DATA / 'models/iYO844.xml',
          'AraCore_v2_1': DATA / 'models/AraCore_v2_1.xml'}
KEYS = ('node_overlap_pairs', 'label_overlap_pairs', 'straight_crossings', 'label_node_intersections', 'segment_label_intersections',
        'canvas_area_per_reaction', 'visible_nodes')


def load_selected(name, condition):
    with gzip.open(RESULTS / 'models' / name / (condition + '-selected.json.gz'), 'rt') as f:
        return json.load(f)


def choose_views(name, selections):
    views = read(OUT / 'models' / name / 'all_geometry' / 'views.json')
    fixed = sorted({0, len(views) // 2, len(views) - 1})
    ec_reactions = {r for p in selections['G_D_E'] if p['match_level'] == 'ec_topology' for r in p['reaction_ids']}
    richest = max(range(len(views)), key=lambda i: (len(ec_reactions & set(views[i]['reaction_ids'])), -i))
    chosen = [(i, 'fixed_first_middle_last') for i in fixed]
    if richest not in fixed:
        chosen.append((richest, 'most_ec_selected_reactions'))
    else:
        chosen = [(i, why + ('+most_ec_selected_reactions' if i == richest else '')) for i, why in chosen]
    return views, chosen, ec_reactions


def run(names=None):
    routing_arch.use_native_binaries()  # no-op on the arm64 machine the checked-in binaries were built on
    summary_path = RESULTS / 'drawing_summary.json'
    summary = read(summary_path) if summary_path.exists() else dict(models={})
    for name, path in MODELS.items():
        if names and name not in names:
            continue
        t0 = time.time()
        model = load_model(path)
        selections = {c: load_selected(name, c) for c in CONDITIONS}
        views, chosen, ec_reactions = choose_views(name, selections)
        rows = []
        for index, why in chosen:
            view = views[index]
            sub = drawing.subset(model, view['reaction_ids'])
            in_view_ec = len(ec_reactions & set(view['reaction_ids']))
            for condition in CONDITIONS:
                clipped = me.clip_proposals(selections[condition], sub)
                module_reactions = {r for p in clipped for r in p['reaction_ids']}
                ec_in_view = {r for p in clipped if p['match_level'] == 'ec_topology' for r in p['reaction_ids']}
                folder = RESULTS / 'drawings' / name / view['id'] / condition
                print('Draw', name, view['id'], condition, flush=True)
                try:
                    metric = me.draw_case(sub, clipped, folder, full_model=model)
                    row = dict(status='ok', **{k: metric.get(k) for k in KEYS}, routed_crossings=metric['routing']['metrics']['proper_crossings'],
                               bends=metric['routing']['metrics']['bends'])
                except Exception as e:  # keep failures as results, never silently drop a condition
                    row = dict(status='failed', error=repr(e)[:400])
                rows.append(dict(model=name, view=view['id'], view_reactions=len(view['reaction_ids']), view_choice=why,
                                 ec_selected_reactions_in_view_total=in_view_ec, condition=condition,
                                 clipped_modules=len(clipped), module_reactions=len(module_reactions), ec_module_reactions=len(ec_in_view),
                                 ec_modules=sum(p['match_level'] == 'ec_topology' for p in clipped),
                                 brenda_modules=sum(p['source'] == 'brenda' for p in clipped), **row))
        summary['models'][name] = dict(rows=rows, seconds=round(time.time() - t0, 1))
        dump(summary_path, summary)
    summary['code_sha256'] = {p: sha((BASE / p).read_bytes()) for p in ('ec_drawing.py', 'ec_channel.py', 'micromap_experiment.py', 'drawing.py', 'geometry.py', 'routing.py')}
    summary['interpretation'] = ('Applications on fixed 50-reaction views; same drawing pipeline for all conditions; EC modules are small and weakly identified; '
                                 'BRENDA modules have unverified layout origin. Not a learning test.')
    dump(summary_path, summary)
    return summary


if __name__ == '__main__':
    run(sys.argv[1:] or None)
