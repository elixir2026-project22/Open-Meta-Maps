"""Generic label → chemical identity resolution for layout importers.

Uses the existing identifier layer (`data/identifier_bridge.json` BiGG↔ChEBI and
public example-model metabolite names/annotations). Callers such as the
CellDesigner importer pass cleaned source labels; Label[compartment] forms are
honoured first when present.
"""
from __future__ import annotations

import re
from collections import defaultdict
from functools import lru_cache

from core import *
import identifier_bridge

# CellDesigner / PANTHER labels embed rudimentary markup instead of Unicode.
_MARKUP = (
    (re.compile(r'_br_', re.I), ' '),
    (re.compile(r'_super_(.+?)_endsuper_', re.I), r'\1'),
    (re.compile(r'_sub_(.+?)_endsub_', re.I), r'\1'),
    (re.compile(r'_italic_(.+?)_enditalic_', re.I), r'\1'),
    (re.compile(r'_plus_', re.I), '+'),
    (re.compile(r'\s+'), ' '),
)
_COMPARTMENT = re.compile(r'^(.+)\[([^\]]+)\]$')
# Leading multiplicity on a single-metabolite glyph ("2 ATP")
_MULT = re.compile(r'^(\d+)\s+(.+)$')


def clean_label(name: str) -> str:
    """Strip CellDesigner markup; keep chemistry tokens."""
    s = (name or '').strip()
    for pat, repl in _MARKUP:
        s = pat.sub(repl, s)
    return s.strip(' ,;')


def normalize_key(name: str) -> str:
    s = clean_label(name).lower()
    s = s.replace('α', 'alpha').replace('β', 'beta').replace('γ', 'gamma')
    s = re.sub(r'^(d|l|dl)-', '', s)
    s = s.replace('-', ' ').replace('_', ' ')
    s = re.sub(r'[^a-z0-9+ ]+', ' ', s)
    s = re.sub(r'\s+', ' ', s).strip()
    # common synonym collapses
    synonyms = {
        'h2o h2o': 'h2o',
        'h+': 'h+',
        'glucose': 'd glucose',
        'g6p': 'glucose 6 phosphate',
        'glucose 6 phosphate': 'glucose 6 phosphate',
        'fructose 6 phosphate': 'fructose 6 phosphate',
        'fructose 1 6 bisphosphate': 'fructose 1 6 bisphosphate',
        'fructose 1 6 phosphate': 'fructose 1 6 bisphosphate',
        'pep': 'phosphoenolpyruvate',
        'nad+': 'nad+',
        'nadh': 'nadh',
        'atp': 'atp',
        'adp': 'adp',
        'amp': 'amp',
        'pi': 'phosphate',
        'orthophosphate': 'phosphate',
        'co2': 'co2',
        'pyruvate': 'pyruvate',
    }
    return synonyms.get(s, s)


def _add(index, key, entry):
    if not key:
        return
    bucket = index[key]
    if entry not in bucket:
        bucket.append(entry)


@lru_cache(maxsize=1)
def label_index():
    """normalized label → list of dict(kind, id, bigg, chebi, models)."""
    index = defaultdict(list)
    chebi_of = identifier_bridge.metabolite_chebi()
    # Prefer curated BiGG display names from the public example models when present.
    model_paths = []
    for name in ('e_coli_core.xml', 'iIT341.xml', 'iMM904.xml', 'iML1515.xml'):
        p = EXAMPLES / name
        if p.exists():
            model_paths.append((name, p))
    for model_name, path in model_paths:
        try:
            model = load_model(path)
        except Exception:
            continue
        for m in model.get('metabolites', []):
            mid = m.get('id') or ''
            bigg = chem(mid)
            name = m.get('name') or ''
            raw = chebi_of.get(bigg) or set()
            chebis = sorted(raw.keys() if isinstance(raw, dict) else raw)
            entry = dict(bigg=bigg, chebi=chebis[:1], models=[model_name], compartment=comp(mid))
            _add(index, normalize_key(name), entry)
            _add(index, normalize_key(bigg), entry)
    # Also index BiGG ids themselves for direct hits
    for bigg, raw in (chebi_of or {}).items():
        chebis = sorted(raw.keys() if isinstance(raw, dict) else raw)
        entry = dict(bigg=bigg, chebi=chebis[:1], models=['identifier_bridge'], compartment=None)
        _add(index, normalize_key(bigg), entry)
    return dict(index)


def resolve_label(name: str, compartment: str | None = None, index=None):
    """Return (identity, meta) where identity is entity@comp or None.

    Preference: Label[compartment] → ChEBI@comp → BiGG@comp → cleaned label@comp.
    The last form keeps observed geometry usable when the bridge has no hit.
    """
    raw = (name or '').strip()
    if not raw:
        return None, dict(status='empty')
    m = _COMPARTMENT.fullmatch(raw)
    if m:
        return f'{m.group(1)}@{m.group(2)}', dict(status='label_compartment', label=m.group(1), compartment=m.group(2))

    cleaned = clean_label(raw)
    mult = None
    mm = _MULT.fullmatch(cleaned)
    if mm and ' + ' not in cleaned:
        mult = int(mm.group(1))
        cleaned = mm.group(2).strip()

    # Composite glyph names ("2 NADH + 2 H+") are not a single identity
    if ' + ' in cleaned:
        return None, dict(status='composite_label', label=cleaned)

    key = normalize_key(cleaned)
    idx = index if index is not None else label_index()
    hits = idx.get(key) or []
    comp_id = compartment or 'local'
    if hits:
        hit = hits[0]
        if hit.get('chebi'):
            cid = hit['chebi'][0]
            if not str(cid).lower().startswith('chebi:'):
                cid = 'chebi:' + str(cid)
            return f'{cid}@{comp_id}', dict(status='chebi_bridge', label=cleaned, bigg=hit.get('bigg'), multiplicity=mult, models=hit.get('models'))
        if hit.get('bigg'):
            return f'{hit["bigg"]}@{comp_id}', dict(status='bigg_bridge', label=cleaned, multiplicity=mult, models=hit.get('models'))

    # Observed label identity — enough for centre-line decode / diagram_modules
    if cleaned:
        return f'{cleaned}@{comp_id}', dict(status='observed_label', label=cleaned, multiplicity=mult, key=key)
    return None, dict(status='unresolved', label=raw)
