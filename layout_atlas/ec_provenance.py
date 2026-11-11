"""EC extended mode: screen EC candidates against the reaction equation and label every reaction with its provenance.

An EC number is only weak evidence: reactions that share an EC can still be different chemistry.  Every EC candidate pair
(source reaction in a library module, target reaction of the uploaded GEM) is therefore checked against the equation when
the evidence allows it:

  verified_equivalent    every participant of the source reaction has an identity in the target model and exactly one
                         stoichiometric assignment maps the whole equation onto the target reaction -> may count as an
                         identity-supported (expert-derived) match; the evidence kind is `ec_plus_equation`
  insufficient_evidence  the source has no stoichiometry, or some participant has no known identity in the target, or the
                         assignment is ambiguous -> the layout is only a `functional_layout_suggestion`; no GEM reaction
                         equivalence is claimed
  incompatible           every participant is identified in the target but the equation is different -> the mapping is
                         rejected (the whole candidate structure is dropped)

Provenance types of a drawn reaction (`reaction_provenance.json`, map legend):
  expert_derived                identity-supported structure (complete equation, reaction identifier, or EC + verified equation)
  topology_layout_suggestion    structure chosen by connection shape only (no identity evidence)
  functional_layout_suggestion  structure chosen through a shared EC number without equation verification
  bridge                        reaction on a network bridge between structures
  generated                     laid out algorithmically
Identity evidence, context compatibility and ambiguity are stored separately for every reaction.

Three judgements are kept apart for every reaction drawn from an expert structure (they are NOT one scale):
  identity_evidence      is the GEM reaction the same reaction as the source reaction?   stoichiometry_verified / identifier_supported / ec_only / topology_only / none
                         (+ `verified`: only a checked equation verifies; an identifier match with a conflicting equation is revoked, not deleted)
  chemical_compatibility are the two equations compatible?  compatible / conflict / unknown, with `difference_types`:
                         exact_equation, direction_difference, proportional_scaling, protonation_difference, cofactor_substitution,
                         core_substrate_difference, compartment_difference, irreversibility_conflict
  layout_provenance      can the expert layout be reused?  expert_exact / expert_adapted / expert_functional / expert_topological / de_novo
                         (`transfer.mode` exact / adapted / topology).  A conflict never reads as "same reaction"; it may still lend its layout (adapted)
                         when an independent similarity holds (same full EC number, or half of the core participants shared).
The five `provenance` colours above are kept and derived from these.
"""
import re
from html import escape, unescape

from core import *
import semantic_transfer as st

KINDS = ['expert_derived', 'topology_layout_suggestion', 'functional_layout_suggestion', 'bridge', 'generated']
COLORS = {'expert_derived': '#2e7d32', 'topology_layout_suggestion': '#8e6bbf', 'functional_layout_suggestion': '#c9a227', 'bridge': '#3d6fa3', 'generated': '#9aa9b2'}
# Display names for the tool-page provenance legend (and SVG <title> tooltips).
# Key ``expert_derived`` is kept for report.json / reaction_provenance.json compatibility.
LABELS = {'expert_derived': 'Curated', 'topology_layout_suggestion': 'Topology-only layout suggestion',
          'functional_layout_suggestion': 'Functional (EC) layout suggestion', 'bridge': 'Network bridge',
          'generated': 'Generated'}


def classify_pair(source_reaction, target_reaction, mapper):
    """Equation check of one EC candidate pair; `mapper` is a semantic_transfer.Mapper (or ECMapper) of the target model."""
    coefficients = (source_reaction or {}).get('coefficients') or {}
    if not coefficients:
        return 'insufficient_evidence', 'the source reaction has no stoichiometry'
    options = {e: mapper.candidates(e) for e in coefficients}
    if any(not v for v in options.values()):
        return 'insufficient_evidence', 'some participants have no known identity in the target model'
    local = {}
    for entity, values in options.items():
        token, _ = st.entity_parts(entity)
        direct = set(mapper.inverted.get(token, set())) & set(target_reaction['metabolites'])
        local[entity] = sorted(direct) if direct else values
    answers = st.unique_equation_mapping(coefficients, target_reaction, local, mapper.compartments)
    if len(answers) == 1:
        return 'verified_equivalent', 'the complete equation maps onto the target reaction'
    if len(answers) > 1:
        return 'insufficient_evidence', 'the equation assignment is ambiguous'
    return 'incompatible', 'all participants are identified but the equations differ'


def screen(proposals, templates, mapper):
    """Keep EC proposals without an incompatible pair; attach `ec_status` (per target reaction) to the kept ones."""
    kept, stats = [], Counter()
    for p in proposals:
        template = templates[p['template_id']]
        source = {r['key']: r for r in template['reactions']}
        statuses, reject = {}, False
        for ev in p['mapping_evidence']:
            sr, tid = source.get(ev['source_key']), ev['target_reaction']
            status, why = classify_pair(sr, mapper.reactions[tid], mapper)
            statuses[tid] = dict(status=status, reason=why, shared_ec=ev.get('shared_ec'), ambiguity=len(mapper.reaction_matches(sr)) if sr else None)
            reject |= status == 'incompatible'
        if reject:
            stats['proposals_rejected_incompatible'] += 1
            continue
        p['ec_status'] = statuses
        kept.append(p)
        stats['proposals_kept'] += 1
        for s in statuses.values():
            stats['reactions_' + s['status']] += 1
    return kept, dict(stats)


# ------------------------------------------------------------------------------------------------ chemical compatibility of two equations
PROTONIC = {'h', 'h2o', 'oh1', 'h2', 'chebi:15377', 'chebi:15378', 'chebi:13365', 'chebi:15379'}
COFACTOR_FAMILIES = [
    {'nad', 'nadh', 'nadp', 'nadph', 'chebi:57540', 'chebi:57945', 'chebi:58349', 'chebi:57783'},
    {'atp', 'adp', 'amp', 'gtp', 'gdp', 'gmp', 'ctp', 'cdp', 'cmp', 'utp', 'udp', 'ump', 'itp', 'idp', 'imp', 'ppi', 'pi', 'chebi:30616', 'chebi:456216', 'chebi:456215'},
    {'fad', 'fadh2', 'fmn', 'fmnh2'},
    {'q8', 'q8h2', 'mqn8', 'mql8', 'mqn7', 'mql7', 'ubiquinone', 'ubiquinol'},
    {'coa', 'chebi:57287', 'chebi:57288'},
    {'o2', 'co2', 'nh4', 'so4', 'h2s', 'chebi:15379', 'chebi:16526', 'chebi:28938'},
]


def _name(entity):
    token, _ = st.entity_parts(entity)
    return (token or '').replace('bigg:', '')


def _is_cofactor(name):
    return name in PROTONIC or any(name in f for f in COFACTOR_FAMILIES) or name in {c for c in CURRENCY}


def _family(name):
    for i, f in enumerate(COFACTOR_FAMILIES):
        if name in f:
            return i
    return None


def target_direction(target):
    """+1 forward-only, -1 reverse-only, 0 reversible or blocked (bounds of the target GEM reaction)."""
    lo, hi = target.get('lower_bound', 0), target.get('upper_bound', 1000)
    if lo >= 0 and hi > 0:
        return 1
    if hi <= 0 and lo < 0:
        return -1
    return 0


def compare_reactions(source, target, mapper):
    """Chemical compatibility of a source (library) reaction with a target GEM reaction.

    Returns dict(status = compatible | conflict | unknown, difference_types = [...], reason, factor, direction_compatible, core_overlap,
    matched / source_unmatched / target_unmatched).  `conflict` only when every participant of the source equation has an identity in the
    target model and the equations still differ; missing identity or stoichiometry is `unknown`, never a conflict.
    """
    coefficients = (source or {}).get('coefficients') or {}
    out = dict(status='unknown', difference_types=[], reason='', factor=None, direction_compatible=None, core_overlap=None)
    if not coefficients:
        out['reason'] = 'the source reaction has no stoichiometry'
        return out
    tmets = set(target['metabolites'])
    options = {e: mapper.candidates(e) for e in coefficients}
    missing = [e for e, v in options.items() if not v]
    if missing:
        out['reason'] = 'some participants have no known identity in the target model'
        return out
    local = {}
    for entity, values in options.items():
        token, _ = st.entity_parts(entity)
        direct = set(mapper.inverted.get(token, set())) & tmets
        local[entity] = sorted(direct) if direct else values
    direct_match = {e: sorted(set(local[e]) & tmets) for e in coefficients}
    answers = st.unique_equation_mapping(coefficients, target, local, mapper.compartments)
    core_src = [e for e in coefficients if not _is_cofactor(_name(e))]
    matched_core = [e for e in core_src if direct_match[e]]
    matched_targets = {m for v in direct_match.values() for m in v}
    core_tgt = [m for m in tmets if not _is_cofactor(chem(m))]
    union = len(core_src) + len(core_tgt) - len(matched_core)
    out['core_overlap'] = round(len(matched_core) / union, 3) if union else 0.0
    if len(answers) == 1:
        a = answers[0]
        f = a['factor']
        out.update(status='compatible', factor=f)
        types = []
        if abs(f - 1) < 1e-9:
            types.append('exact_equation')
        else:
            if f < 0:
                types.append('direction_difference')
            if abs(abs(f) - 1) > 1e-9:
                types.append('proportional_scaling')
        if not a['label_compartments_exact'] and any(comp(e) not in ('unknown', 'local') for e in coefficients):
            types.append('compartment_difference')
        out['difference_types'] = types
        # irreversible on both sides and pointing opposite ways: the equation is the same, the arrow is not
        source_irreversible = (source or {}).get('reversible') is False
        td = target_direction(target)
        if source_irreversible and td != 0:
            out['direction_compatible'] = (td > 0) == (f > 0)
            if not out['direction_compatible']:
                out['difference_types'].append('irreversibility_conflict')
        out['reason'] = 'the complete equation maps onto the target reaction'
        return out
    if len(answers) > 1:
        out['reason'] = 'the equation assignment is ambiguous'
        return out
    # every source participant is identified in the model but the equations differ: say how
    src_unmatched = [e for e in coefficients if not direct_match[e]]
    tgt_unmatched = [m for m in tmets if m not in matched_targets]
    types = []
    names = [_name(e) for e in src_unmatched] + [chem(m) for m in tgt_unmatched]
    if names and all(n in PROTONIC for n in names):
        types.append('protonation_difference')
    elif names and all(_is_cofactor(n) for n in names):
        fams = {_family(n) for n in names}
        types.append('cofactor_substitution' if (fams - {None}) and not all(n in PROTONIC for n in names) else 'protonation_difference')
    else:
        types.append('core_substrate_difference')
        if any(_is_cofactor(n) for n in names):
            types.append('cofactor_substitution')
    out.update(status='conflict', difference_types=types, reason='all participants are identified but the equations differ',
               source_unmatched=sorted(src_unmatched), target_unmatched=sorted(tgt_unmatched))
    return out


def _source_for(p, rid, template):
    """The library reaction that was placed on target reaction `rid` in proposal `p`."""
    key = next((k for k, v in (p.get('reaction_mapping') or {}).items() if v == rid), None)
    if key is None:
        key = next((e.get('source_key') for e in p.get('mapping_evidence', []) if e.get('target_reaction') == rid), None)
    if key is None or template is None:
        return None
    return next((r for r in template.get('reactions', []) if r.get('key') == key or r.get('key') == 'r:' + str(key)), None)


ADAPT_MIN_CORE_OVERLAP = 0.5   # independent similarity needed for a conflicting reaction to lend its layout ("adapted")


def reaction_provenance(model, selected, bridge_ids, codes=None, templates=None, mapper=None):
    """One record per GEM reaction: the five-colour `provenance` plus the three independent judgements (see the module docstring).

    `templates` ({id: library module}) and `mapper` (a semantic_transfer.Mapper of the target) enable the equation comparison; without them
    only the identifier / EC evidence recorded on the proposals is used (chemical compatibility stays `unknown`).
    """
    import ec_channel as ec
    target_ec = ec.target_ecs(model) if mapper is not None else {}
    reactions = {r['id']: r for r in model['reactions']}
    owner = {}
    for p in selected:
        for r in p['reaction_ids']:
            owner[r] = p
    out = {}
    for r in model['reactions']:
        rid = r['id']
        rec = dict(reaction=rid, provenance='generated', identity_evidence=None, chemical_compatibility=None, layout_provenance='de_novo', transfer=None,
                   context_compatibility=None, ambiguity=None, structure=None)
        p = owner.get(rid)
        if p is None:
            if rid in bridge_ids:
                rec['provenance'] = 'bridge'
            out[rid] = rec
            continue
        level = p['match_level']
        ev = {e.get('target_reaction'): e for e in p.get('mapping_evidence', [])}.get(rid, {})
        template = (templates or {}).get(p.get('template_id'))
        source = _source_for(p, rid, template)
        comp_info = compare_reactions(source, reactions[rid], mapper) if (source is not None and mapper is not None) else dict(status='unknown', difference_types=[], reason='no source stoichiometry available', factor=None, direction_compatible=None, core_overlap=None)
        src_ec = ec.source_ecs(source)[0] if source is not None else set()
        shared_ec = sorted(src_ec & target_ec.get(rid, set()))
        rec['structure'] = dict(structure_id=p['id'], code=(codes or {}).get(p['id']), source=p['source'], match_level=level, license=p.get('license'), partition=p.get('partition', 'open_core'))
        rec['context_compatibility'] = dict(incidence_fit=True, compartments_verified=ev.get('compartments_verified'))
        equation_ok = bool(p.get('complete_chemical_equation_verified') or ev.get('complete_equation'))
        status = comp_info['status']
        if equation_ok:
            status = 'compatible'
        # --- identity -----------------------------------------------------------------------------------------------------------
        identity = dict(kind='topology_only', verified=False, namespaces=[], shared_ec=shared_ec, revoked=False)
        kind = 'topology_layout_suggestion'
        if level == 'ec_topology':
            st_ = (p.get('ec_status') or {}).get(rid, {})
            rec['ambiguity'] = st_.get('ambiguity')
            shared_ec = sorted(set(shared_ec) | set(st_.get('shared_ec') or []))
            identity.update(shared_ec=shared_ec, namespaces=['ec'])
            if st_.get('status') == 'verified_equivalent':
                identity.update(kind='stoichiometry_verified', verified=True, basis='ec_plus_equation')
                kind = 'expert_derived'
            else:
                identity.update(kind='ec_only', reason=st_.get('reason'))
                kind = 'functional_layout_suggestion'
        elif equation_ok:
            identity.update(kind='stoichiometry_verified', verified=True, basis='complete_equation', namespaces=['equation'])
            kind = 'expert_derived'
        elif level != 'geometry' and ev.get('reference'):
            identity.update(kind='identifier_supported', verified=False, basis='reaction_identifier', namespaces=['reaction_identifier'])
            if status == 'compatible':
                identity.update(kind='stoichiometry_verified', verified=True, basis='identifier_plus_equation')
                kind = 'expert_derived'
            elif status == 'conflict':
                identity.update(revoked=True, reason='the identifier matches but the equations differ: identity is not accepted')
                kind = 'topology_layout_suggestion'
            else:
                kind = 'expert_derived'  # identifier match with nothing against it; not equation-verified
        elif shared_ec and status != 'conflict':
            identity.update(kind='ec_only', namespaces=['ec'], basis='shared_full_ec')
            kind = 'functional_layout_suggestion'
        # --- layout provenance / transfer mode --------------------------------------------------------------------------------------
        independent = bool(shared_ec) or (comp_info.get('core_overlap') or 0) >= ADAPT_MIN_CORE_OVERLAP
        if identity['verified'] or (identity['kind'] == 'identifier_supported' and not identity['revoked']):
            layout, mode = 'expert_exact', 'exact'
        elif identity['kind'] == 'ec_only' and not identity['revoked']:
            layout, mode = 'expert_functional', 'adapted' if comp_info['status'] == 'conflict' else 'topology'
        elif comp_info['status'] == 'conflict' and independent:
            layout, mode = 'expert_adapted', 'adapted'
            kind = 'functional_layout_suggestion' if shared_ec else 'topology_layout_suggestion'
        else:
            layout, mode = 'expert_topological', 'topology'
        if equation_ok and not comp_info['difference_types']:
            comp_info = dict(comp_info, status='compatible', difference_types=['exact_equation'])
        elif equation_ok:
            comp_info = dict(comp_info, status='compatible')
        rec['provenance'] = kind
        rec['identity_evidence'] = identity
        rec['chemical_compatibility'] = dict(status=status, difference_types=comp_info['difference_types'], reason=comp_info['reason'], factor=comp_info['factor'],
                                             direction_compatible=comp_info['direction_compatible'], core_overlap=comp_info['core_overlap'])
        rec['layout_provenance'] = layout
        rec['transfer'] = dict(mode=mode, source_structure=p.get('template_id'), source_map=p.get('map_id'), source_kind=p['source'], source_species=p.get('species'),
                               source_reaction=(source or {}).get('id') or ev.get('source_reaction'),
                               basis=dict(shared_ec=shared_ec, core_overlap=comp_info['core_overlap'], difference_types=comp_info['difference_types']), adjustments=[])
        out[rid] = rec
    return out


def evidence_summary(provenance):
    """Counts of every dimension over a reaction_provenance() result, for report.json."""
    c = dict(identity_evidence=Counter(), identity_verified=Counter(), chemical_compatibility=Counter(), difference_types=Counter(), layout_provenance=Counter(), transfer_mode=Counter())
    for r in provenance.values():
        c['layout_provenance'][r['layout_provenance']] += 1
        if r['identity_evidence'] is None:
            continue
        c['identity_evidence'][r['identity_evidence']['kind']] += 1
        c['identity_verified']['verified' if r['identity_evidence']['verified'] else 'not_verified'] += 1
        c['chemical_compatibility'][r['chemical_compatibility']['status']] += 1
        for t in r['chemical_compatibility']['difference_types']:
            c['difference_types'][t] += 1
        c['transfer_mode'][r['transfer']['mode']] += 1
    out = {k: dict(v) for k, v in c.items()}
    out['identity_revoked'] = sum(1 for r in provenance.values() if r['identity_evidence'] and r['identity_evidence'].get('revoked'))
    return out


def recolor_svg(text, kinds):
    """Colour the reaction markers of an exported map by provenance (reaction markers are the orange circles followed by their label).

    Provenance counts are shown in the Map Generator UI, not drawn into the SVG.
    """
    def repl(m):
        label = unescape(m.group(4))
        kind = kinds.get(label)
        if kind is None:
            return m.group(0)
        return (f'<circle cx="{m.group(1)}" cy="{m.group(2)}" r="5.5" fill="{COLORS[kind]}"><title>{escape(label)} - {LABELS[kind]}</title></circle>'
                f'<text {m.group(3)}>{m.group(4)}</text>')
    return re.sub(r'<circle cx="([^"]+)" cy="([^"]+)" r="4" fill="#d77d38"/><text ([^>]*)>([^<]*)</text>', repl, text)
