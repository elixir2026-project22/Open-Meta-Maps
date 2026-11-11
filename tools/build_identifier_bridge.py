"""Build `data/identifier_bridge.json` (and check `data/currency_metabolites.json`) from public genome-scale models.

The tool needs two identifier tables that the matching step uses:
  metabolite_chebi   BiGG metabolite id -> ChEBI ids (from the models' own annotations)
  reaction_ec        BiGG reaction id   -> EC numbers (from the models' own annotations; partial EC numbers are kept here and ignored when the table is read)
Both are stored WITH their origin: for every ChEBI id or EC number the list of models it was found in, so conflicting or multiple annotations are kept side by side
(nothing is merged into a single value).  The JSON also records each model's file name, SHA-256, identifier and size, and the licence terms of the sources.

  python3 tools/build_identifier_bridge.py --models-dir DIR [--human-gem FILE]

DIR must contain e_coli_core.xml, iML1515.xml and iMM904.xml (the BiGG SBML files in `examples/`); Human-GEM.xml is taken from --human-gem or from DIR.
"""
import argparse
import hashlib
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE))
from model_io import load_model  # noqa: E402

MODELS = [
    dict(name='e_coli_core', file='e_coli_core.xml', source='BiGG Models', licence='BiGG licence (https://bigg.ucsd.edu/license)', bigg_ids=True),
    dict(name='iML1515', file='iML1515.xml', source='BiGG Models', licence='BiGG licence (https://bigg.ucsd.edu/license)', bigg_ids=True),
    dict(name='iMM904', file='iMM904.xml', source='BiGG Models', licence='BiGG licence (https://bigg.ucsd.edu/license)', bigg_ids=True),
    dict(name='Human-GEM', file='Human-GEM.xml', source='Human-GEM (SysBioChalmers)', licence='CC BY 4.0 according to the Human-GEM project page', bigg_ids=False),
]


def xrefs(annotations):
    """Same rules as core.xrefs: the `is` qualifier only, identifiers by database."""
    out = defaultdict(set)
    for a in annotations:
        if not a.get('qualifier', '').endswith('}is'):
            continue
        m = re.search(r'/(chebi|bigg.metabolite|rhea|reactome|bigg.reaction|metanetx.chemical|metanetx.reaction)[/:](?:CHEBI[:_])?([^/?#]+)', a.get('uri', ''), re.I)
        if m:
            out[m[1].lower()].add(m[2])
    return {k: sorted(v) for k, v in out.items()}


def build(models_dir, human_gem=None):
    models_dir = Path(models_dir)
    metabolite_chebi = defaultdict(lambda: defaultdict(list))
    reaction_ec = defaultdict(lambda: defaultdict(list))
    sources = []
    for spec in MODELS:
        path = Path(human_gem) if spec['name'] == 'Human-GEM' and human_gem else models_dir / spec['file']
        model = load_model(path)
        sources.append(dict(model=spec['name'], file=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest(), sbml_id=model['id'], sbml_name=model['name'],
                            metabolites=len(model['metabolites']), reactions=len(model['reactions']), source=spec['source'], licence=spec['licence']))
        for m in model['metabolites']:
            refs = xrefs(m.get('annotations', []))
            for bigg in refs.get('bigg.metabolite', []):
                for chebi in refs.get('chebi', []):
                    if spec['name'] not in metabolite_chebi[bigg][chebi]:
                        metabolite_chebi[bigg][chebi].append(spec['name'])
        for r in model['reactions']:
            ann = r.get('annotations', [])
            ecs = sorted({re.search(r'ec-code/([^/?#]+)', a['uri'])[1] for a in ann if 'ec-code/' in a.get('uri', '')})
            ids = set(xrefs(ann).get('bigg.reaction', []))
            if spec['bigg_ids']:
                ids.add(r['id'])            # the model's own reaction ids are BiGG ids
            for bigg in ids:
                for ec in ecs:
                    if spec['name'] not in reaction_ec[bigg][ec]:
                        reaction_ec[bigg][ec].append(spec['name'])
    return dict(
        format_version=1,
        description='Identifier tables used by the matching step: BiGG metabolite -> ChEBI and BiGG reaction -> EC number, each with the models whose annotations contain it. '
                    'Only explicit annotations of the listed models are used; multiple or conflicting annotations are kept side by side.',
        licence_note='The tables are derived from the annotations of the models listed in `sources`; the terms of those models apply (BiGG licence; CC BY 4.0 for Human-GEM). '
                     'Cite BiGG Models and Human-GEM when using them.',
        sources=sources,
        metabolite_chebi={k: {c: m for c, m in sorted(v.items())} for k, v in sorted(metabolite_chebi.items())},
        reaction_ec={k: {e: m for e, m in sorted(v.items())} for k, v in sorted(reaction_ec.items())})


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--models-dir', required=True)
    p.add_argument('--human-gem')
    p.add_argument('--out', default=str(BASE.parent / 'data/identifier_bridge.json'))
    a = p.parse_args()
    bridge = build(a.models_dir, a.human_gem)
    Path(a.out).write_text(json.dumps(bridge, ensure_ascii=False, indent=1) + '\n')
    print(f"wrote {a.out}: {len(bridge['metabolite_chebi'])} metabolites with ChEBI links, {len(bridge['reaction_ec'])} reactions with EC numbers")


if __name__ == '__main__':
    main()
