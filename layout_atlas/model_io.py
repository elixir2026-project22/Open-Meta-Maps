"""Model input: SBML (read with libSBML) and COBRA / BiGG JSON -> one plain dictionary.

`load_model(path)` returns

  dict(id, name, annotations, notes,
       metabolites=[dict(id, name, compartment, annotations, notes[, chemicalFormula, charge])],
       reactions=[dict(id, name, annotations, notes, gene_reaction_rule, metabolites={metabolite id: coefficient}, lower_bound, upper_bound)],
       genes=[dict(id, name, label, annotations, notes)])

  annotations   list of dict(qualifier='{http://biomodels.net/biology-qualifiers/}is', uri='https://identifiers.org/...'): the element's own CV terms only
  coefficient   negative for reactants, positive for products, summed when a species appears twice, zero entries dropped

The input rules are the ones the layout code relies on and are strict on purpose: exactly one model, no rules / events / initial assignments, constant numeric
stoichiometry, finite non-negative coefficients, valid flux bounds, unique ids.  BiGG-style escaped ids (`M_`, `R_`, `G_` prefixes and `__45__` character codes) are
decoded.  JSON input (COBRA / BiGG JSON) is returned as it is.
"""
import json
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path

BIOL_NS = '{http://biomodels.net/biology-qualifiers/}'
MODEL_NS = '{http://biomodels.net/model-qualifiers/}'


def _decode(identifier, prefix):
    if identifier.startswith(prefix):
        identifier = identifier[len(prefix):]
    return re.sub(r'__(\d+)__', lambda m: chr(int(m[1])), identifier)


def _notes_text(element):
    raw = element.getNotesString() if element.isSetNotes() else ''
    if not raw:
        return ''
    try:
        return ' '.join(' '.join(ET.fromstring(raw).itertext()).split())
    except ET.ParseError:
        return ' '.join(re.sub(r'<[^>]+>', ' ', raw).split())


def _annotations(element):
    import libsbml
    out = []
    for i in range(element.getNumCVTerms()):
        term = element.getCVTerm(i)
        if term.getQualifierType() == libsbml.BIOLOGICAL_QUALIFIER:
            tag = BIOL_NS + libsbml.BiolQualifierType_toString(term.getBiologicalQualifierType())
        else:
            tag = MODEL_NS + libsbml.ModelQualifierType_toString(term.getModelQualifierType())
        for j in range(term.getNumResources()):
            out.append(dict(qualifier=tag, uri=term.getResourceURI(j)))
    return out


def _metadata(element):
    return dict(annotations=_annotations(element), notes=_notes_text(element))


def _gene_rule(association):
    if association is None:
        return ''
    if association.isGeneProductRef():
        return _decode(association.getGeneProduct(), 'G_')
    if association.isFbcAnd() or association.isFbcOr():
        joiner = ' and ' if association.isFbcAnd() else ' or '
        return '(' + joiner.join(_gene_rule(association.getAssociation(i)) for i in range(association.getNumAssociations())) + ')'
    return ''


def load_model(path):
    path = Path(path)
    if path.suffix.lower() == '.json':
        return json.loads(path.read_text(encoding='utf-8'))
    import libsbml
    doc = libsbml.readSBMLFromFile(str(path))
    log = doc.getErrorLog()
    fatal = [log.getError(i).getMessage() for i in range(log.getNumErrors()) if log.getError(i).getSeverity() >= libsbml.LIBSBML_SEV_FATAL]
    if fatal:
        raise ValueError('Invalid SBML XML: ' + ' '.join(fatal[0].split()))
    model = doc.getModel()
    if model is None:
        first = log.getError(0).getMessage() if log.getNumErrors() else ''
        raise ValueError('Invalid SBML XML or not an SBML document with one model' + (': ' + ' '.join(first.split())[:200] if first else ''))
    if model.getNumRules() or model.getNumEvents() or model.getNumInitialAssignments():
        raise ValueError('Dynamic SBML rules, events and initial assignments are unsupported')
    species = [model.getSpecies(i) for i in range(model.getNumSpecies())]
    id_map = {s.getId(): _decode(s.getId(), 'M_') for s in species}
    if len(id_map) != len(species) or len(set(id_map.values())) != len(species):
        raise ValueError('Duplicate metabolite IDs after SBML ID normalization')
    metabolites = []
    for s in species:
        record = dict(id=id_map[s.getId()], name=s.getName(), compartment=s.getCompartment(), **_metadata(s))
        plugin = s.getPlugin('fbc')
        if plugin is not None:
            if plugin.isSetChemicalFormula():
                record['chemicalFormula'] = plugin.getChemicalFormula()
            if plugin.isSetCharge():
                record['charge'] = str(plugin.getCharge())
        metabolites.append(record)
    genes = []
    model_fbc = model.getPlugin('fbc')
    if model_fbc is not None:
        for i in range(model_fbc.getNumGeneProducts()):
            g = model_fbc.getGeneProduct(i)
            genes.append(dict(id=_decode(g.getId(), 'G_'), name=g.getName(), label=g.getLabel(), **_metadata(g)))
    parameters = {p.getId(): p.getValue() for p in (model.getParameter(i) for i in range(model.getNumParameters())) if p.isSetValue()}
    reactions = []
    for k in range(model.getNumReactions()):
        reaction = model.getReaction(k)
        rid = _decode(reaction.getId(), 'R_')
        coefficients = {}
        for sign, count, get in ((-1, reaction.getNumReactants(), reaction.getReactant), (1, reaction.getNumProducts(), reaction.getProduct)):
            for i in range(count):
                ref = get(i)
                if ref.isSetStoichiometryMath() or (ref.isSetConstant() and not ref.getConstant()):
                    raise ValueError(f'Variable stoichiometry is unsupported: {rid}')
                raw = ref.getSpecies()
                if raw not in id_map:
                    raise ValueError(f'Unknown species {raw} in reaction {rid}')
                value = ref.getStoichiometry() if ref.isSetStoichiometry() else 1.0
                if not math.isfinite(value) or value < 0:
                    raise ValueError(f'Invalid stoichiometry in reaction {rid}')
                mid = id_map[raw]
                coefficients[mid] = coefficients.get(mid, 0) + sign * value
        # COBRA Level 2 stores bounds as kinetic-law parameters.
        legacy = {}
        law = reaction.getKineticLaw()
        if law is not None:
            for i in range(law.getNumParameters()):
                p = law.getParameter(i)
                if p.isSetValue():
                    legacy[p.getId()] = p.getValue()
        reversible = reaction.getReversible() if reaction.isSetReversible() else True
        bounds = [legacy.get('LOWER_BOUND', -1000 if reversible else 0), legacy.get('UPPER_BOUND', 1000)]
        plugin = reaction.getPlugin('fbc')
        rule = ''
        if plugin is not None:
            for i, (is_set, get) in enumerate(((plugin.isSetLowerFluxBound(), plugin.getLowerFluxBound), (plugin.isSetUpperFluxBound(), plugin.getUpperFluxBound))):
                if is_set:
                    ref = get()
                    if ref not in parameters:
                        raise ValueError(f'Unknown flux-bound parameter {ref} in {rid}')
                    bounds[i] = parameters[ref]
            if plugin.isSetGeneProductAssociation():
                rule = _gene_rule(plugin.getGeneProductAssociation().getAssociation())
        if any(math.isnan(b) for b in bounds) or bounds[0] > bounds[1]:
            raise ValueError(f'Invalid flux bounds in reaction {rid}')
        reactions.append(dict(id=rid, name=reaction.getName(), **_metadata(reaction), gene_reaction_rule=rule,
                              metabolites={m: c for m, c in coefficients.items() if c}, lower_bound=bounds[0], upper_bound=bounds[1]))
    if not metabolites or not reactions or len({r['id'] for r in reactions}) != len(reactions):
        raise ValueError('Model must contain metabolites and reactions with unique IDs')
    return dict(id=model.getId(), name=model.getName(), **_metadata(model), metabolites=metabolites, reactions=reactions, genes=genes)
