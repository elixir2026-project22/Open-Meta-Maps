"""Runtime access to `data/identifier_bridge.json` (built by `tools/build_identifier_bridge.py`) plus shared metabolite-name matching helpers.

The JSON holds BiGG metabolite -> ChEBI and BiGG reaction -> EC.  Name normalisation, weighted Jaccard
scoring and one-to-one metabolite matching are used when Atlas nodes only carry metabolite *names*
(BRENDA pathway map, some CellDesigner drawings).  Method ideas follow ecmapper (GPL-3.0; ideas only —
this file is an independent rewrite for OpenMetaMaps).  Label→identity for importers lives in
`identity_labels` (uses this bridge); do not duplicate that index here.
"""
import json
import re
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

FILE = Path(__file__).resolve().parent.parent / 'data/identifier_bridge.json'
FULL_EC = re.compile(r'^\d+\.\d+\.\d+\.\d+$')

# --- metabolite name keys (name-only map sources) -------------------------------------------------
_STEREO = re.compile(r'(?<![a-z0-9])(?:all-trans|\((?:[0-9]*[rsez][,0-9rsez]*)\)|dl|[dl]|cis|trans)-')
# anomeric alpha-/beta- only when followed by D-/L- (keep beta-alanine, alpha-ketoglutarate)
_ANOMER = re.compile(r'(?<![a-z0-9])(?:alpha|beta)-(?=(?:[nos]-[a-z]+-)?(?:dl|[dl])-)')
_SIDE = re.compile(r'\s*\[(?:side \d+|in|out|inside|outside)\]$')
_EITHER = re.compile(r'\(([a-z])\)', re.I)
_PROTON_NAMES = frozenset({'h+', 'h', 'proton', 'protons', 'hydron', 'hydrogenion'})
# ChEBI / BiGG ids that are protons (ignored in scoring, as in the BRENDA-map method)
_PROTON_IDS = frozenset({'chebi:15378', 'chebi:24636', 'chebi:10744', 'chebi:13357', 'bigg:h'})
# generic redox carriers on maps ("a quinone") vs specific model species
_GROUPS = {
    'quinone': re.compile(r'(?:a|an|\d*)(?:ubi|mena|demethylmena|plasto)?quinone\d*|coenzymeq\d*|oxidizedquinone'),
    'quinol': re.compile(r'(?:a|an|\d*)(?:ubi|mena|demethylmena|plasto)?quinol\d*|hydroquinone|reduced(?:ubi|mena)?quinone\d*'),
}

# Default acceptance thresholds (principled / cross-model; same values as the published BRENDA-map
# evaluation on eleven curated prokaryote models).  Not tuned to a single GEM.
THRESHOLDS = dict(
    min_score=0.6,                 # EC- or reaction-id-linked candidate
    min_score_unlinked=0.8,        # metabolites-only candidate
    max_gap=0.1,                   # keep near-best nodes (same reaction in several subsystems)
    min_score_vicinity=0.4,        # vicinity with shared EC / reaction id
    min_score_vicinity_unlinked=0.6,
    currency_weight=0.5,
    ec_min_score=0.6,              # metabolite scoring when selecting among candidate ECs
    ec_min_margin=0.1,
)

# BiGG / common currency ids used as half-weight keys (aligned with data/currency_metabolites.json roles)
_CURRENCY_BIGG = frozenset({
    'h2o', 'h', 'pi', 'ppi', 'coa', 'co2', 'o2', 'nh4', 'nh3', 'h2o2',
    'atp', 'adp', 'amp', 'gtp', 'gdp', 'gmp', 'ctp', 'cdp', 'cmp', 'utp', 'udp', 'ump',
    'nad', 'nadh', 'nadp', 'nadph', 'fad', 'fadh2', 'q8', 'q8h2', 'mqn8', 'mql8',
})


@lru_cache(maxsize=1)
def _bridge():
    return json.loads(FILE.read_text())


def metabolite_chebi():
    """BiGG metabolite id -> set of ChEBI ids (union over the models)."""
    return {k: set(v) for k, v in _bridge()['metabolite_chebi'].items()}


def bigg_to_ec():
    """BiGG reaction id -> set of fully specified EC numbers (partial numbers such as 1.1.1.- are ignored), union over the models."""
    table = defaultdict(set)
    for bigg, ecs in _bridge()['reaction_ec'].items():
        table[bigg] |= {e for e in ecs if FULL_EC.match(e)}
    return table


def sources():
    return _bridge()['sources']


@lru_cache(maxsize=1)
def bigg_names():
    """BiGG metabolite id -> set of display names from currency table and bridge source models.

    Used so BRENDA / map labels (e.g. \"NAD+\", \"1,3-diphosphoglycerate\") meet BiGG's longer
    chemical names without KEGG.  Built only from open files already in the repository.
    """
    out = defaultdict(set)
    currency_path = FILE.parent / 'currency_metabolites.json'
    if currency_path.exists():
        for e in json.loads(currency_path.read_text()).get('currency', []):
            if e.get('id') and e.get('name'):
                out[e['id']].add(e['name'])
            # short id itself as a name token (atp, nadh, …)
            if e.get('id'):
                out[e['id']].add(e['id'])
    # Common open abbreviations that maps use instead of BiGG's long names
    for bigg, aliases in {
        'nad': ('NAD+', 'NAD', 'nicotinamide adenine dinucleotide'),
        'nadh': ('NADH',),
        'nadp': ('NADP+', 'NADP'),
        'nadph': ('NADPH',),
        'coa': ('CoA', 'coenzyme A', 'CoA-SH'),
        'pi': ('phosphate', 'orthophosphate', 'Pi'),
        'ppi': ('diphosphate', 'pyrophosphate', 'PPi'),
        '13dpg': ('1,3-diphosphoglycerate', '1,3-bisphosphoglycerate', '3-phospho-D-glyceroyl phosphate'),
        '3pg': ('3-phosphoglycerate', '3-Phospho-D-glycerate'),
        '2pg': ('2-phosphoglycerate', '2-Phospho-D-glycerate'),
        'pep': ('phosphoenolpyruvate',),
        'g3p': ('D-glyceraldehyde 3-phosphate', 'glyceraldehyde 3-phosphate'),
        'dhap': ('dihydroxyacetone phosphate', 'glycerone phosphate'),
        'fdp': ('fructose-1,6-bisphosphate', 'D-fructose 1,6-bisphosphate'),
        'f6p': ('D-fructose 6-phosphate', 'fructose 6-phosphate'),
        'g6p': ('D-glucose 6-phosphate', 'glucose 6-phosphate'),
        'akg': ('2-oxoglutarate', 'alpha-ketoglutarate', '2-Ketoglutarate'),
        'oaa': ('oxaloacetate',),
        'accoa': ('acetyl-CoA', 'acetyl-CoA'),
        'icit': ('isocitrate',),
        'cit': ('citrate',),
        'succ': ('succinate',),
        'fum': ('fumarate',),
        'mal__L': ('(S)-malate', 'L-malate', 'malate'),
        'glu__L': ('L-glutamate', 'glutamate'),
        'gln__L': ('L-glutamine', 'glutamine'),
        'pyr': ('pyruvate',),
        'lac__D': ('D-lactate',),
        'lac__L': ('L-lactate',),
        'acald': ('acetaldehyde',),
        'etoh': ('ethanol',),
        'ac': ('acetate',),
        'actp': ('acetyl phosphate',),
        'q8': ('ubiquinone', 'ubiquinone-8', 'a quinone'),
        'q8h2': ('ubiquinol', 'ubiquinol-8', 'a quinol'),
    }.items():
        out[bigg].update(aliases)
    return {k: set(v) for k, v in out.items()}


def _normalize(name, stereo=False):
    n = _SIDE.sub('', (name or '').lower())
    if not stereo:
        n = _STEREO.sub('', _ANOMER.sub('', n))
    n = n.replace('bisphosph', 'diphosph').replace('ic acid', 'ate')
    return re.sub(r'[^a-z0-9+]', '', n)


def name_key(name):
    """Compare metabolite names ignoring case, punctuation, membrane sides and stereo/anomer descriptors."""
    return _normalize(name)


def name_keys(name):
    """Stereo-stripped key plus a key that keeps stereo descriptors (racemases stay two metabolites)."""
    return {_normalize(name), _normalize(name, stereo=True)}


def expand_names(names):
    """Expand optional letters in parentheses: NAD(P)H -> NADH and NADPH as well."""
    out = set()
    for name in names:
        if not name:
            continue
        out.add(name)
        if _EITHER.search(name):
            out |= {_EITHER.sub('', name), _EITHER.sub(r'\1', name)}
    return out


def metabolite_keys(names=(), compounds=(), bigg_ids=()):
    """Frozen set of comparison keys for one metabolite: compound ids, name keys, generic groups.

    Protons are omitted (empty frozenset): they do not participate in the score.
    """
    names = {n for n in expand_names(names) if n and name_key(n) not in _PROTON_NAMES}
    compounds = {c if str(c).startswith(('chebi:', 'bigg:')) else f'chebi:{c}' for c in compounds}
    compounds -= _PROTON_IDS
    bigg_ids = {b for b in bigg_ids if b and b != 'h'}
    if not names and not compounds and not bigg_ids:
        return frozenset()
    keys = set(compounds) | {f'bigg:{b}' for b in bigg_ids}
    for n in names:
        keys |= {f'name:{k}' for k in name_keys(n) if k and k not in _PROTON_NAMES}
        for g, pattern in _GROUPS.items():
            if pattern.fullmatch(name_key(n)):
                keys.add(f'group:{g}')
    keys -= _PROTON_IDS
    if not keys or keys <= {f'name:{k}' for k in _PROTON_NAMES}:
        return frozenset()
    return frozenset(keys)


def is_currency_keys(keys, currency_bigg=_CURRENCY_BIGG):
    if not keys:
        return True
    if keys & {f'bigg:{c}' for c in currency_bigg}:
        return True
    if keys & {f'group:{g}' for g in _GROUPS}:
        return True
    # common currency name keys
    currency_names = {
        'water', 'atp', 'adp', 'amp', 'nad', 'nadh', 'nadp', 'nadph', 'coa', 'coash',
        'phosphate', 'orthophosphate', 'diphosphate', 'pyrophosphate', 'co2', 'carbondioxide',
        'o2', 'oxygen', 'nh3', 'ammonia', 'nh4', 'ammonium', 'h2o2', 'hydrogenperoxide',
    }
    return bool({k[5:] for k in keys if k.startswith('name:')} & currency_names)


def weighted_jaccard(left, right, currency_weight=None):
    """Weighted Jaccard with one-to-one matching (maximum matching preferring non-currency).

    `left` / `right` are lists of frozensets (one per metabolite).  Currency metabolites
    count `currency_weight` (default half).  Each metabolite matches at most one partner.
    """
    if currency_weight is None:
        currency_weight = THRESHOLDS['currency_weight']
    if not left and not right:
        return 0.0
    weights = [currency_weight if is_currency_keys(k) else 1.0 for k in left]
    partner = {}

    def augment(i, visited):
        for j, other in enumerate(right):
            if j not in visited and left[i] & other:
                visited.add(j)
                if j not in partner or augment(partner[j], visited):
                    partner[j] = i
                    return True
        return False

    for i in sorted(range(len(left)), key=lambda i: -weights[i]):
        augment(i, set())
    matched = sum(weights[i] for i in partner.values())
    union = sum(weights) + sum(
        currency_weight if is_currency_keys(other) else 1.0
        for j, other in enumerate(right) if j not in partner
    )
    return matched / union if union else 0.0


def accept_score(score, linked, best=None, thresholds=None):
    """Whether a candidate score is accepted (linked = shared EC or reaction id)."""
    t = thresholds or THRESHOLDS
    floor = t['min_score'] if linked else t['min_score_unlinked']
    if score < floor:
        return False
    if best is not None and score < best - t['max_gap']:
        return False
    return True
