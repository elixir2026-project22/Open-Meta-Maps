"""Window library vs structure-unit library with exactly the same matching and selection code.

Libraries (all through identical code)
  W      current open-core window library (17,221 overlapping BFS windows); geometry candidates = the saved
         `results/models/<m>/all_geometry/candidates.json.gz` (output of `retrieval.Library(core).candidates`, as in ec_experiment)
  M      maximal-only: windows not contained in a larger window of the same map (results/module_redundancy.json.gz), stored as is
  M_sub  M plus sub-windows derived from the maximal modules (BFS inside one parent, coordinates cropped from the parent)
  U      structure units (units_library.py; non-overlapping, 3-12 reactions + appendices, geometry-only maps by shape)
  U_sub  U plus sub-windows derived from the units (hybrid: de-duplicated store, windows derived on demand)

Every library goes through identical code: `retrieval.Library.candidates` (geometry/reference candidates), complete-equation
`semantic_transfer.Mapper`, `ec_channel.ECMapper`, `ec_channel.select_with_ec` (-> `retrieval.select`), `retrieval.bridges`.
Conditions: D = geometry+equation pool, DE = D plus EC gap filling (as G_D / G_D_E in ec_experiment); S / SE = equation (+EC) pool only.

  python3 units_experiment.py retrieval [model ...]      # plain python3 is enough
  PATH=$HOME/.local/node/bin:$PATH /private/tmp/claude-501/venv2/bin/python units_experiment.py drawing [model ...]
  python3 units_experiment.py report
"""
import copy
import gzip
import sys
import time

from core import *
import database
import semantic_transfer as st
import ec_channel as ec
import retrieval
import units_library as ul

RESULTS = OUT / 'units_experiment'
MODELS = [('e_coli_core', REF / 'models/e_coli_core.xml', 'Escherichia coli'), ('iML1515', REF / 'models/iML1515.xml', 'Escherichia coli'),
          ('iMM904', REF / 'models/iMM904.xml', 'Saccharomyces cerevisiae'), ('Human-GEM', REF / 'models/Human-GEM.xml', 'Homo sapiens'),
          ('iYO844', DATA / 'models/iYO844.xml', 'Bacillus subtilis'), ('AraCore_v2_1', DATA / 'models/AraCore_v2_1.xml', 'Arabidopsis thaliana')]
DRAW_MODELS = {'iML1515': REF / 'models/iML1515.xml', 'iMM904': REF / 'models/iMM904.xml', 'iYO844': DATA / 'models/iYO844.xml', 'AraCore_v2_1': DATA / 'models/AraCore_v2_1.xml'}
LIBS = ['W', 'M', 'M_sub', 'U', 'U_sub']
CONDS = ['D', 'DE', 'S', 'SE']
QUERY = {lib: (2, 4, 8, 12) for lib in LIBS}
KEEP = ('license', 'source_sha256', 'original_instance', 'source_citation')


KIND = {'W': 'windows', 'M': 'maximal', 'M_sub': 'maximal_sub', 'U': 'units', 'U_sub': 'units_sub'}


def libraries():
    return {lib: ul.modules(kind) for lib, kind in KIND.items()}


def tag(pool, templates):
    for p in pool:
        t = templates[p['template_id']]
        for k in KEEP:
            if k in t:
                p[k] = t[k]
        p.setdefault('partition', t.get('partition', 'open_core'))
    return pool


def metrics(model, selected, allids):
    covered = {r for p in selected for r in p['reaction_ids']}
    chem_cov = {r for p in selected if p.get('complete_chemical_equation_verified', p.get('complete_equation_verified')) for r in p['reaction_ids']}
    sizes = [len(p['reaction_ids']) for p in selected]
    bridges = retrieval.bridges(model, selected)
    n = max(1, len(covered))
    return dict(selected_modules=len(selected), selected_reactions=len(covered), coverage=len(covered) / len(allids), chemical_coverage=len(chem_cov) / len(allids),
                modules_per_100_covered=100 * len(selected) / n, mean_unit_size=float(np.mean(sizes)) if sizes else 0.,
                size_histogram=dict(sorted(Counter(sizes).items())), modules_with_2_or_less=sum(s <= 2 for s in sizes), bridges=len(bridges),
                bridge_reactions=len({r for b in bridges for r in b['reaction_ids']}), distinct_source_maps=len({p['map_id'] for p in selected}),
                match_levels=dict(Counter(p['match_level'] for p in selected)), sources=dict(Counter(p['source'] for p in selected)))


def run_model(name, path, bridge, libs):
    t0 = time.time()
    model = load_model(path)
    allids = {r['id'] for r in model['reactions']}
    g = model_graph(model)
    folder = RESULTS / 'models' / name
    folder.mkdir(parents=True, exist_ok=True)
    mapper, ecmapper = st.Mapper(model, bridge), ec.ECMapper(model, bridge)
    with gzip.open(OUT / 'models' / name / 'all_geometry' / 'candidates.json.gz', 'rt') as f:
        saved = json.load(f)
    row = dict(model=name, reactions=len(allids), libraries={})
    for lib in LIBS:
        t1 = time.time()
        mods = copy.deepcopy(libs[lib])
        templates = {m['id']: m for m in mods}
        ec.enrich(mods)
        if lib == 'W':
            cands, audit = saved, dict(query_windows=None, candidates=len(saved), source='saved all_geometry candidates')
        else:
            library = retrieval.Library(mods)
            cands, audit = library.candidates(g, model=model, max_roots=120, query_sizes=QUERY[lib])
            audit = {k: v for k, v in audit.items() if k != 'queried_reaction_ids'}
            tag(cands, templates)
        strong, sinfo = mapper.proposals(mods, lib + ' complete-equation')
        ecp, einfo = ecmapper.proposals(mods, lib + ' EC')
        tag(strong, templates); tag(ecp, templates)
        pools = dict(D=(cands + strong, []), DE=(cands + strong, [ecp]), S=(strong, []), SE=(strong, [ecp]))
        info = dict(library_modules=len(mods), library_reactions=sum(len(m['reactions']) for m in mods), geometry_candidates=len(cands), retrieval=audit,
                    equation_candidates=sinfo['candidates'], equation_candidate_reactions=len(sinfo['candidate_reaction_ids']), equation_budget_reached=sinfo['budget_reached'],
                    ec_candidates=einfo['candidates'], ec_candidate_reactions=len(einfo['candidate_reaction_ids']), ec_budget_reached=einfo['budget_reached'], conditions={})
        for cond, (base, extra) in pools.items():
            selected, stages = ec.select_with_ec(base, extra)
            info['conditions'][cond] = dict(metrics(model, selected, allids), stage_modules=[len(s) for s in stages], ec_selected_modules=sum(p['match_level'] == 'ec_topology' for p in selected))
            with gzip.open(folder / f'{lib}-{cond}-selected.json.gz', 'wt') as f:
                json.dump(selected, f, ensure_ascii=False)
        info['seconds'] = round(time.time() - t1, 1)
        row['libraries'][lib] = info
        print(name, lib, {c: (v['selected_reactions'], v['selected_modules'], round(v['coverage'], 3), v['bridges']) for c, v in info['conditions'].items()}, info['seconds'], flush=True)
        dump(folder / 'summary.json', dict(row, seconds=round(time.time() - t0, 1)))
    row['seconds'] = round(time.time() - t0, 1)
    dump(folder / 'summary.json', row)
    return row


def retrieval_stage(names=None):
    RESULTS.mkdir(parents=True, exist_ok=True)
    bridge = st.annotation_bridge(MODELS)
    libs = libraries()
    for name, path, _ in MODELS:
        if names and name not in names:
            continue
        run_model(name, path, bridge, libs)


# ---------------------------------------------------------------- drawing (needs the drawing stack)
DRAW_KEYS = ('node_overlap_pairs', 'label_overlap_pairs', 'straight_crossings', 'label_node_intersections', 'segment_label_intersections', 'canvas_area_per_reaction', 'visible_nodes')
DRAW_LIBS = ['W', 'M_sub', 'U', 'U_sub']


def drawing_stage(names=None):
    import drawing
    import micromap_experiment as me
    import routing_arch
    routing_arch.use_native_binaries()
    for name, path in DRAW_MODELS.items():
        if names and name not in names:
            continue
        t0 = time.time()
        summary_path = RESULTS / 'drawing' / f'{name}.json'
        summary = dict(models={})
        model = load_model(path)
        sel = {}
        for lib in DRAW_LIBS:
            with gzip.open(RESULTS / 'models' / name / f'{lib}-DE-selected.json.gz', 'rt') as f:
                sel[lib] = json.load(f)
        views = read(OUT / 'models' / name / 'all_geometry' / 'views.json')
        chosen = sorted({0, len(views) // 2, len(views) - 1})  # same fixed first/middle/last rule as ec_drawing, chosen before drawing
        rows = []
        for index in chosen:
            view = views[index]
            sub = drawing.subset(model, view['reaction_ids'])
            for lib in DRAW_LIBS:
                clipped = me.clip_proposals(sel[lib], sub)
                folder = RESULTS / 'drawings' / name / view['id'] / lib
                print('Draw', name, view['id'], lib, flush=True)
                try:
                    metric = me.draw_case(sub, clipped, folder, full_model=model)
                    row = dict(status='ok', **{k: metric.get(k) for k in DRAW_KEYS}, routed_crossings=metric['routing']['metrics']['proper_crossings'], bends=metric['routing']['metrics']['bends'])
                except Exception as e:
                    row = dict(status='failed', error=repr(e)[:400])
                reacts = {r for p in clipped for r in p['reaction_ids']}
                rows.append(dict(model=name, view=view['id'], view_reactions=len(view['reaction_ids']), library=lib, clipped_modules=len(clipped), module_reactions=len(reacts),
                                 mean_clipped_size=float(np.mean([len(p['reaction_ids']) for p in clipped])) if clipped else 0., **row))
        summary['models'][name] = dict(rows=rows, seconds=round(time.time() - t0, 1))
        summary['code_sha256'] = {p: sha((BASE / p).read_bytes()) for p in ('units_experiment.py', 'units_library.py', 'micromap_experiment.py', 'drawing.py', 'geometry.py', 'routing.py')}
        summary['interpretation'] = 'DE selections of each library drawn on the same fixed first/middle/last 50-reaction views with the same pipeline. Application, not a learning test.'
        dump(summary_path, summary)


# ---------------------------------------------------------------- acceptance cases (stage 1)
CASE_MAPS = dict(reactome_pufa='d53a9ded377ab025303afc66', wp_purine_disorders='c1abc063f57574a7f98c0fc1', wp_purine='e672cafe9bb312138d5586e8',
                 ijo1366_central='e63b0b880daf3cfba5d13db3', recon_glycolysis_tca_ppp='7261af45735b7f5d90d5ccca')
TCA = {'CS', 'ACONT', 'ICDHyr', 'AKGDH', 'SUCOAS', 'SUCDi', 'FUM', 'MDH'}


def describe(m, by_id):
    return dict(unit=m['id'][:8], shape=m['shape_type'], size=m['size'], appendix=len(m['appendix_reactions']), family=m['core_family'][:6], variant=m['variant'][:6],
                reactions=[r['id'] for r in m['reactions']],
                ports=[dict(metabolite=p['metabolite'], outside=len(p['outside_reactions']), units=[by_id[u]['shape_type'] + str(by_id[u]['size']) for u in p['neighbour_units'] if u in by_id]) for p in m['ports']])


def cases():
    units = ul.load()
    by_id = {m['id']: m for m in units}
    window_counts = Counter(m['map_id'] for m in database.modules())
    out = {}
    for label, map_id in CASE_MAPS.items():
        mine = sorted((m for m in units if m['map_id'] == map_id), key=lambda m: (-m['size'], m['id']))
        out[label] = dict(map_id=map_id, name=mine[0]['name'] if mine else None, window_modules_in_current_library=window_counts[map_id], units=len(mine),
                          unit_reactions=sum(m['size'] for m in mine), shapes=dict(Counter(m['shape_type'] for m in mine)), sizes=dict(sorted(Counter(m['size'] for m in mine).items())),
                          detail=[describe(m, by_id) for m in mine] if label in ('reactome_pufa',) else None)
    # (b) units of the two purine maps with the same reaction set -> family / variant after exhaustive alignment
    a = [m for m in units if m['map_id'] == CASE_MAPS['wp_purine_disorders']]
    b = [m for m in units if m['map_id'] == CASE_MAPS['wp_purine']]
    pairs = []
    for x in a:
        for y in b:
            rx, ry = {r['id'] for r in x['reactions']}, {r['id'] for r in y['reactions']}
            if len(rx & ry) >= max(3, min(len(rx), len(ry)) // 2):
                gx, gy = database.graph(x), database.graph(y)
                ex, cx = ul.all_mappings(gx, gy, CONFIG_CAP)
                e8 = ex[:8]
                err = lambda ms: min(database.geometry_error(x, y, mp) for mp in ms) if ms else None
                pairs.append(dict(shape_a=x['shape_type'], shape_b=y['shape_type'], size_a=x['size'], size_b=y['size'], shared_reactions=len(rx & ry), isomorphisms_found=len(ex), cap_hit=bool(cx),
                                  best_error_first8=err(e8), best_error_exhaustive=err(ex), same_family=x['core_family'] == y['core_family'], same_variant=x['variant'] == y['variant']))
    out['purine_pairs'] = pairs
    # (c) TCA cycle and chains in iJO1366 central metabolism
    c = [m for m in units if m['map_id'] == CASE_MAPS['ijo1366_central']]
    out['tca_units'] = [describe(m, by_id) | dict(tca_reactions=sorted(TCA & {r['id'] for r in m['reactions']})) for m in c if TCA & {r['id'] for r in m['reactions']}]
    chains = sorted((m for m in units if m['shape_type'] == 'chain' and m['source'] == 'community'), key=lambda m: (-m['size'], m['id']))[:3]
    out['longest_chains_escher'] = [dict(describe(m, by_id), map=m['name']) for m in chains]
    out['cycle_examples_escher'] = [dict(describe(m, by_id), map=m['name']) for m in sorted((m for m in c if m['shape_type'] == 'cycle'), key=lambda m: (-m['size'], m['id']))[:4]]
    dump(RESULTS / 'acceptance_cases.json', out)
    return out


CONFIG_CAP = ul.CONFIG['mapping_cap']


# ---------------------------------------------------------------- coverage-loss diagnosis
def diagnose():
    out = {}
    for name, _, _ in MODELS:
        folder = RESULTS / 'models' / name
        if not (folder / 'U-DE-selected.json.gz').exists():
            continue
        def load_sel(lib):
            with gzip.open(folder / f'{lib}-DE-selected.json.gz', 'rt') as f:
                return json.load(f)
        sel = {lib: load_sel(lib) for lib in LIBS if (folder / f'{lib}-DE-selected.json.gz').exists()}
        owner = {lib: {r: p for p in ps for r in p['reaction_ids']} for lib, ps in sel.items()}
        lost = sorted(set(owner['W']) - set(owner['U']))
        row = dict(reactions_covered_W=len(owner['W']), reactions_covered_U=len(owner['U']), lost_by_U=len(lost), gained_by_U=len(set(owner['U']) - set(owner['W'])))
        row['lost_by_size_of_W_module'] = dict(sorted(Counter(len(owner['W'][r]['reaction_ids']) for r in lost).items()))
        row['lost_covered_only_by_W_modules_of_size_le_2'] = sum(len(owner['W'][r]['reaction_ids']) <= 2 for r in lost)
        row['lost_by_W_match_level'] = dict(Counter(owner['W'][r]['match_level'] for r in lost))
        row['lost_by_W_source'] = dict(Counter(owner['W'][r]['source'] for r in lost))
        for lib in ('U_sub', 'M', 'M_sub'):
            if lib in owner:
                row['recovered_by_' + lib] = sum(r in owner[lib] for r in lost)
        out[name] = row
    dump(RESULTS / 'coverage_diagnosis.json', out)
    return out


# ---------------------------------------------------------------- report
def report():
    rows = {}
    for name, _, _ in MODELS:
        p = RESULTS / 'models' / name / 'summary.json'
        if p.exists():
            rows[name] = read(p)
    out = dict(models=rows, decision_inputs=decide(rows))
    dump(RESULTS / 'summary.json', out)
    print(json.dumps(out['decision_inputs'], indent=1))


def decide(rows, base='W', new='U', cond='DE'):
    """Decision rule: coverage drop <= 2 pp, fragmentation (modules per 100 covered reactions) down >= 30 %, drawing not worse in most views."""
    out = {}
    for name, r in rows.items():
        a, b = r['libraries'][base]['conditions'][cond], r['libraries'][new]['conditions'][cond]
        out[name] = dict(coverage_change_pp=round(100 * (b['coverage'] - a['coverage']), 2), fragmentation_change_pct=round(100 * (b['modules_per_100_covered'] / a['modules_per_100_covered'] - 1), 1),
                         coverage_ok=(a['coverage'] - b['coverage']) <= .02, fragmentation_ok=b['modules_per_100_covered'] <= .7 * a['modules_per_100_covered'])
    return out


if __name__ == '__main__':
    if os.environ.get('PYTHONHASHSEED') != '0':
        os.environ['PYTHONHASHSEED'] = '0'
        os.execv(sys.executable, [sys.executable, str(Path(__file__).resolve())] + sys.argv[1:])
    stage, args = (sys.argv[1] if len(sys.argv) > 1 else 'retrieval'), sys.argv[2:] or None
    {'retrieval': retrieval_stage, 'drawing': drawing_stage, 'cases': lambda a: print(json.dumps(cases(), indent=1)[:12000]), 'diagnose': lambda a: print(json.dumps(diagnose(), indent=1))}.get(stage, lambda a: report())(args)
