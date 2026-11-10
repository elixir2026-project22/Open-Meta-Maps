"""EC-number evidence channel for module transfer (opt-in, additive).

A source module reaction is linked to a target-model reaction when both carry the same *fully specified*
EC number (a.b.c.d with four numeric levels).  EC is a weak identity: one EC can cover several reactions and many
reactions (transport, exchange, biomass, spontaneous) have no EC.  The channel therefore

* never replaces the strong channel (complete-equation matching in `semantic_transfer`);
* requires the module's reaction/metabolite incidence to fit the target network (graph check), not just EC equality;
* maps a core metabolite only when it is uniquely determined (identity candidates or structure) and leaves the rest
  unmapped instead of guessing;
* labels every proposal `match_level='ec_topology'`, biochemical level 1, with the shared EC numbers as evidence.

Existing retrieval, ranking and saved results are untouched; this module only adds candidates.
"""
import csv
import itertools
import re
from functools import lru_cache

from core import *
import database
import semantic_transfer as st

CONFIG = dict(module_combinations=16, proposal_budget=6000, min_module_size=2, max_unmapped_fraction=.5)
FULL_EC = re.compile(r'^\d+\.\d+\.\d+\.\d+$')


def full_ecs(values):
    return {v for v in values if isinstance(v, str) and FULL_EC.match(v)}


@lru_cache(maxsize=1)
def rhea_to_ec():
    table = defaultdict(set)
    with (DATA / 'rhea/rhea2xrefs.tsv').open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['DB'] == 'EC':
                for key in (row['RHEA_ID'], row['MASTER_ID']):
                    table[key].add(row['ID'])
    return table


@lru_cache(maxsize=1)
def bigg_to_ec():
    """BiGG reaction id -> fully specified EC numbers, from `data/identifier_bridge.json` (explicit annotations of the public BiGG / Human-GEM models)."""
    import identifier_bridge
    return identifier_bridge.bigg_to_ec()


def source_ecs(reaction):
    """Fully specified EC numbers of a source module reaction, with their origin."""
    origin = defaultdict(set)
    for v in reaction.get('references', {}).get('ec-code', []):
        origin['source_annotation'].add(str(v))
    for v in reaction.get('references', {}).get('rhea', []):
        origin['rhea2xrefs'] |= rhea_to_ec().get(str(v), set())
    for v in reaction.get('references', {}).get('bigg.reaction', []):
        origin['bigg_annotation'] |= bigg_to_ec().get(str(v), set())
    merged = full_ecs(set().union(*origin.values())) if origin else set()
    return merged, {k: sorted(full_ecs(v)) for k, v in origin.items() if full_ecs(v)}


def enrich(modules):
    """Attach `ec_numbers` / `ec_origin` to module reactions in memory (no stored data is rewritten)."""
    have = total = 0
    for m in modules:
        for r in m.get('reactions', []):
            ecs, origin = source_ecs(r)
            r['ec_numbers'], r['ec_origin'] = sorted(ecs), origin
            total += 1
            have += bool(ecs)
    return dict(reactions=total, with_full_ec=have)


def target_ecs(model):
    out = {}
    for r in model['reactions']:
        ann = r.get('annotations', [])
        out[r['id']] = full_ecs(re.search(r'ec-code/([^/?#]+)', a.get('uri', ''))[1] for a in ann if 'ec-code/' in a.get('uri', ''))
    return out


class ECMapper(st.Mapper):
    """Subclass of the equation mapper so metabolite identity candidates can be reused as a filter."""

    def __init__(self, model, bridge):
        super().__init__(model, bridge)
        self.ec = target_ecs(model)
        self.by_ec = defaultdict(set)
        for rid, ecs in self.ec.items():
            for e in ecs:
                self.by_ec[e].add(rid)
        self.graph = model_graph(model)
        self.audit_ec = Counter()

    def reaction_matches(self, reaction):
        ecs = set(reaction.get('ec_numbers') or ())
        out = defaultdict(set)
        for e in ecs:
            for rid in self.by_ec.get(e, ()):
                out[rid].add(e)
        return [dict(target_reaction=rid, shared_ec=sorted(v)) for rid, v in sorted(out.items())]

    def metabolite_options(self, node, parent_targets):
        """Target metabolites for a module metabolite: neighbours shared by all parent reactions, filtered by identity if known."""
        domain = None
        for rid in parent_targets:
            nb = {k for k in self.graph['r:' + rid]}
            domain = nb if domain is None else domain & nb
        if not domain:
            return [], 'no_common_target_metabolite'
        identity = {'m:' + t for t in self.candidates(node.get('entity'))} if isinstance(node.get('entity'), str) else set()
        if identity & domain:
            return sorted(identity & domain), 'identity_filtered'
        return sorted(domain), 'structure_only'

    def proposals(self, modules, source_label):
        g = self.graph
        skips, out, matched = Counter(), [], set()
        support = defaultdict(set)
        for m in modules:
            support[m['core_family']].add(m['lineage'])

        @lru_cache(maxsize=4096)
        def initial(rids):
            keys = {'r:' + r for r in rids}
            for r in rids:
                keys.update(g['r:' + r])
            sub = nx.Graph()
            sub.add_nodes_from((k, dict(g.nodes[k])) for k in sorted(keys))
            sub.add_edges_from((a, b, dict(g.edges[a, b])) for a, b in sorted(tuple(sorted(e)) for e in g.subgraph(keys).edges()))
            return sub, nx.spring_layout(sub, seed=SEEDS[0], iterations=35, scale=1.)

        examined = 0
        for m in sorted(modules, key=lambda m: (-m['size'], m['id'])):
            if m['size'] < CONFIG['min_module_size'] or not m.get('reactions') or len(m['reactions']) != m['size']:
                continue
            examined += 1
            choices = [self.reaction_matches(r) for r in m['reactions']]
            if any(not c for c in choices):
                skips['source_reaction_without_ec_match'] += 1
                continue
            keys = {r['key'] for r in m['reactions']}
            parents = defaultdict(set)
            for e in m['edges']:
                if e['a'] in keys:
                    parents[e['b']].add(e['a'])
                if e['b'] in keys:
                    parents[e['a']].add(e['b'])
            for combo in itertools.islice(itertools.product(*choices), CONFIG['module_combinations']):
                rids = tuple(sorted(c['target_reaction'] for c in combo))
                if len(set(rids)) != len(rids):
                    skips['duplicate_target_reaction'] += 1
                    continue
                rmap = {r['key']: c['target_reaction'] for r, c in zip(m['reactions'], combo)}
                node_map = {k: 'r:' + v for k, v in rmap.items()}
                evidence, ambiguous, used = [], [], set()
                # Fit check: every module reaction pair that shares a core metabolite must share a target metabolite.
                fit = True
                for n in m['nodes']:
                    if n['role'] == 'reaction':
                        continue
                    pk = [rmap[p] for p in parents[n['key']] if p in rmap]
                    options, how = self.metabolite_options(n, pk)
                    if how == 'no_common_target_metabolite':
                        fit = False
                        break
                    options = [o for o in options if o not in used]
                    if len(options) == 1:
                        node_map[n['key']] = options[0]
                        used.add(options[0])
                        evidence.append(dict(node=n['key'], target=options[0], how=how))
                    else:
                        ambiguous.append(n['key'])
                if not fit:
                    skips['core_incidence_does_not_fit_target'] += 1
                    continue
                if any(not g.has_edge(node_map[e['a']], node_map[e['b']]) for e in m['edges'] if e['a'] in node_map and e['b'] in node_map):
                    skips['mapped_edge_missing_in_target'] += 1
                    continue
                if len(node_map) < 3:
                    skips['too_few_mapped_nodes'] += 1
                    continue
                tg, generic = initial(rids)
                source = {n['key']: n for n in m['nodes']}
                ks = sorted(node_map)
                frame = fit_frame([source[k]['xy'] for k in ks], [generic[node_map[k]] for k in ks])
                positions = transform({node_map[k]: source[k]['xy'] for k in ks}, frame)
                targetports = {k: len(set(g[k]) - {'r:' + rid for rid in rids}) for k in tg if g.nodes[k]['role'] == 'metabolite'}
                srcports, ports = len(m.get('ports', [])), sum(bool(v) for v in targetports.values())
                initial_fit = float(np.sqrt(np.mean([np.sum((positions[k] - generic[k]) ** 2) for k in positions])))
                xy = np.array(list(positions.values()))
                extent = np.ptp(xy, axis=0)
                area = max(float(np.prod(extent)), .05)
                spread = float(np.std(xy[:, 0]) + np.std(xy[:, 1]))
                vec = [np.array(positions[node_map[e['a']]]) - positions[node_map[e['b']]] for e in m['edges'] if e['a'] in node_map and e['b'] in node_map]
                axis = float(np.mean([min(abs(v[0]), abs(v[1])) / max(np.linalg.norm(v), 1e-8) < .2 for v in vec])) if vec else 0.
                features = [len(rids), 1, abs(ports - srcports), 0., float(extent[0] / max(extent[1], .05)), len(tg) / area, axis,
                            len(support[m['core_family']]), initial_fit, sum(targetports.values()) * float(extent[0] / max(extent[1], .05)),
                            ports / max(len(tg), 1) * spread]
                checks = [dict(source_reaction=r['id'], source_key=r['key'], target_reaction=c['target_reaction'], shared_ec=c['shared_ec'],
                               source_ec_origin=r.get('ec_origin', {}), target_ec=sorted(self.ec[c['target_reaction']]))
                          for r, c in zip(m['reactions'], combo)]
                p = dict(id=uid(['ec-topology', m['id'], rids, sorted(node_map.items())]), template_id=m['id'], map_id=m['map_id'], lineage=m['lineage'],
                         source=m['source'], species=m['species'], license=m['license'], source_citation=m.get('source_citation'),
                         partition=m.get('partition', 'open_core'), geometry_fingerprint=m['geometry_fingerprint'], reaction_ids=list(rids),
                         node_mapping=node_map, reaction_mapping=rmap, positions={k: np.array(v).tolist() for k, v in positions.items()},
                         position_origins={node_map[k]: source[k]['position_origin'] for k in ks}, source_to_initial_frame=frame,
                         match_level='ec_topology', complete_chemical_equation_verified=False, complete_equation_verified=False,
                         label_compartments_verified=False, mapping_evidence=checks, metabolite_evidence=evidence,
                         unmapped_source_nodes=sorted(ambiguous), source_node_count=len(m['nodes']), features=features,
                         score=2 * len(rids) + .6 - .1 * features[2] - .2 * initial_fit - .3 * len(ambiguous) / max(1, len(m['nodes'])),
                         expert_query_coordinates_used=False, template_quality=m['quality'], match_mode='ec_number_plus_incidence_fit',
                         layout_origin=m.get('layout_origin'), counted_as_expert_design=m['source'] != 'brenda',
                         unmapped_target_nodes=sorted(set(tg) - set(node_map.values())),
                         ec_policy='fully specified EC equality; weak identity; strong equation channel outranks it')
                out.append(p)
                matched.update(rids)
                if len(out) >= CONFIG['proposal_budget']:
                    break
            if len(out) >= CONFIG['proposal_budget']:
                break
        return out, dict(library=source_label, examined_multi_modules=examined, candidates=len(out), proposal_budget=CONFIG['proposal_budget'],
                         budget_reached=len(out) >= CONFIG['proposal_budget'], skips=dict(skips), candidate_reaction_ids=sorted(matched),
                         candidate_reactions=len(matched), config=CONFIG,
                         interpretation='EC-supported candidates with an incidence-fit check; not a certification of chemical equivalence')


def usable(proposal):
    """EC proposals whose source module is mostly mapped; heavily unmapped modules carry little layout information."""
    n = max(1, proposal.get('source_node_count', 1))
    return len(proposal['unmapped_source_nodes']) / n <= CONFIG['max_unmapped_fraction']


def select_with_ec(strong_pool, ec_pools):
    """Strong-first selection.

    1. `retrieval.select` over the existing strong pool, unchanged.
    2. For each EC pool in order (expert-library EC first, unverified extensions last), select only among proposals
       that do not touch an already covered reaction, again with `retrieval.select`.
    EC therefore fills gaps; it never displaces or overlaps a complete-equation or geometry selection.
    Returns (all selected, per-stage selected lists).
    """
    import copy
    import retrieval
    stages = [retrieval.select(copy.deepcopy(strong_pool))]
    covered = {r for p in stages[0] for r in p['reaction_ids']}
    for pool in ec_pools:
        free = [p for p in pool if usable(p) and not set(p['reaction_ids']) & covered]
        chosen = retrieval.select(copy.deepcopy(free))
        stages.append(chosen)
        covered |= {r for p in chosen for r in p['reaction_ids']}
    return [p for stage in stages for p in stage], stages
