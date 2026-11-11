"""Generic SBGN-ML (process description) decoder for Atlas diagram_modules.

Parses glyphs and arcs into the same diagram form as CellDesigner / Escher
importers. Callers pass source / licence / name / species metadata. Layout
coordinates are taken from glyph bboxes (entity centres, process centres);
arc routing waypoints are retained on events when present.
"""
from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from urllib.parse import unquote, urlparse

from core import *

# Entity / EP glyphs that can be stoichiometric participants
ENTITY_CLASSES = {
    'simple chemical',
    'macromolecule',
    'macromolecule multimer',
    'complex',
    'complex multimer',
    'nucleic acid feature',
    'nucleic acid feature multimer',
    'unspecified entity',
    'perturbing agent',
    'empty set',
}
PROCESS_CLASSES = {
    'process',
    'omitted process',
    'uncertain process',
    'association',
    'dissociation',
}
MODIFIER_ARCS = {
    'catalysis',
    'stimulation',
    'inhibition',
    'necessary stimulation',
    'modulation',
}
SUBSTRATE_ARCS = {'consumption'}
PRODUCT_ARCS = {'production'}

_IDENT_RE = re.compile(
    r'(?:identifiers\.org/|purl\.obolibrary\.org/obo/|bioregistry\.io/)?'
    r'(CHEBI[:_]|chebi:|uniprot:|UNIPROT:|hgnc\.symbol:|ncbigene:|ec-code:|RHEA:|rhea:|kegg\.compound:|kegg\.reaction:)?'
    r'([A-Za-z0-9_.:\-]+)'
)


def _local(tag):
    return tag.rsplit('}', 1)[-1] if tag else tag


def _bbox_center(bbox):
    return [
        float(bbox.get('x', 0)) + float(bbox.get('w', 0)) / 2.0,
        float(bbox.get('y', 0)) + float(bbox.get('h', 0)) / 2.0,
    ]


def _label_of(glyph):
    for c in glyph:
        if _local(c.tag) == 'label' and c.get('text'):
            return c.get('text').strip()
    return ''


def _bbox_of(glyph):
    for c in glyph:
        if _local(c.tag) == 'bbox':
            return {k: float(c.get(k, 0)) for k in ('x', 'y', 'w', 'h')}
    return None


def _child_glyphs(glyph):
    return [c for c in glyph if _local(c.tag) == 'glyph']


def _resources(glyph):
    """Collect identifiers.org / RDF resource URIs under a glyph."""
    out = []
    for e in glyph.iter():
        for key, val in e.attrib.items():
            if key.endswith('resource') or key == 'resource' or _local(key) == 'resource':
                if val:
                    out.append(val)
        if e.text and 'identifiers.org' in e.text:
            out.append(e.text.strip())
    return out


def _normalize_ident(uri):
    """Map a resource URI or CURIE to a compact id string, or None."""
    if not uri:
        return None
    uri = unquote(uri.strip())
    if uri.startswith('#'):
        return None
    # identifiers.org/CHEBI:58165 or .../chebi/58165
    m = re.search(r'identifiers\.org/([^/\s]+)/([^/\s?#]+)', uri, re.I)
    if m:
        ns, val = m.group(1).lower(), m.group(2)
        if ns in ('chebi',):
            return f'chebi:{val.replace("CHEBI:", "").replace("chebi:", "")}'
        if ns in ('uniprot', 'uniprotkb'):
            return f'uniprot:{val}'
        if ns in ('hgnc.symbol', 'hgnc'):
            return f'hgnc.symbol:{val}'
        if ns in ('ncbigene', 'ncbi.gene'):
            return f'ncbigene:{val}'
        if ns in ('ec-code', 'ec'):
            return f'ec-code:{val}'
        if ns == 'rhea':
            return f'RHEA:{val}'
        if ns in ('kegg.compound', 'kegg.reaction'):
            return f'{ns}:{val}'
        return f'{ns}:{val}'
    m = re.search(r'obo/CHEBI[_:](\d+)', uri, re.I)
    if m:
        return f'chebi:{m.group(1)}'
    m = re.match(r'(chebi|uniprot|RHEA|rhea|ec-code|hgnc\.symbol|ncbigene):(.+)', uri, re.I)
    if m:
        ns = m.group(1)
        if ns.lower() == 'rhea':
            ns = 'RHEA'
        elif ns.lower() == 'chebi':
            ns = 'chebi'
        return f'{ns}:{m.group(2)}'
    return None


def _entity_id(glyph, label):
    idents = []
    for uri in _resources(glyph):
        ident = _normalize_ident(uri)
        if ident and ident not in idents:
            idents.append(ident)
    if idents:
        # Prefer ChEBI for chemicals, else first
        chebi = [i for i in idents if i.startswith('chebi:')]
        primary = chebi[0] if chebi else idents[0]
        return primary + '@local', idents
    if label:
        return label + '@local', []
    return glyph.get('id') + '@local', []


def _port_owner(port_or_glyph_id, glyph_ids, port_to_owner):
    if port_or_glyph_id in glyph_ids:
        return port_or_glyph_id
    if port_or_glyph_id in port_to_owner:
        return port_to_owner[port_or_glyph_id]
    # SBGN port ids are typically "<glyphId>.<n>"
    if '.' in port_or_glyph_id:
        base = port_or_glyph_id.rsplit('.', 1)[0]
        if base in glyph_ids:
            return base
    return port_or_glyph_id


def parse(path):
    """Parse an SBGN-ML document into glyphs, arcs and map metadata."""
    path = Path(path)
    root = ET.parse(path).getroot()
    issues = Counter()
    maps = [e for e in root.iter() if _local(e.tag) == 'map']
    map_el = maps[0] if maps else root
    language = map_el.get('language') or ''
    map_id = map_el.get('id') or path.stem

    glyphs = {}
    port_to_owner = {}
    # Only top-level map glyphs (not state variables nested inside)
    for g in map_el:
        if _local(g.tag) != 'glyph':
            continue
        gid = g.get('id')
        if not gid:
            issues['glyph_without_id'] += 1
            continue
        cls = (g.get('class') or '').strip()
        bbox = _bbox_of(g)
        if bbox is None:
            issues['glyph_without_bbox'] += 1
            continue
        label = _label_of(g)
        for pg in _child_glyphs(g):
            if (pg.get('class') or '') == 'port' or pg.get('id'):
                pass
        for c in g:
            if _local(c.tag) == 'port' and c.get('id'):
                port_to_owner[c.get('id')] = gid
        # Nested state / unit-of-information glyphs are not map nodes
        kind = None
        if cls in PROCESS_CLASSES:
            kind = 'process'
        elif cls in ENTITY_CLASSES:
            kind = 'entity'
        elif cls in ('compartment',):
            kind = 'compartment'
        elif cls in ('tag', 'annotation', 'submap'):
            kind = 'other'
        else:
            # nested-looking classes at top level
            if cls in ('state variable', 'unit of information', 'port'):
                kind = 'annotation'
            else:
                kind = 'other'
                issues['unhandled_glyph_class:' + (cls or '?')] += 1
        entity, idents = (None, [])
        if kind == 'entity':
            entity, idents = _entity_id(g, label)
        glyphs[gid] = dict(
            id=gid,
            class_=cls,
            kind=kind,
            label=label,
            bounds=bbox,
            xy=_bbox_center(bbox),
            entity=entity,
            identifiers=idents,
            ports=[c.get('id') for c in g if _local(c.tag) == 'port' and c.get('id')],
        )

    arcs = []
    for a in map_el:
        if _local(a.tag) != 'arc':
            continue
        cls = (a.get('class') or '').strip()
        src = a.get('source')
        tgt = a.get('target')
        if not src or not tgt:
            issues['arc_missing_endpoint'] += 1
            continue
        points = []
        for c in a:
            tag = _local(c.tag)
            if tag in ('start', 'end', 'next'):
                try:
                    points.append([float(c.get('x', 0)), float(c.get('y', 0))])
                except (TypeError, ValueError):
                    issues['arc_bad_point'] += 1
        arcs.append(dict(
            id=a.get('id'),
            class_=cls,
            source=src,
            target=tgt,
            points=points,
        ))
    return dict(
        map_id=map_id,
        language=language,
        glyphs=glyphs,
        arcs=arcs,
        port_to_owner=port_to_owner,
        issues=dict(issues),
        path=str(path),
    )


def decode(
    path,
    source='sbgn',
    license='CC-BY-4.0',
    name=None,
    species=None,
    lineage=None,
    partition='open_core',
    admitted_to_open_core=True,
    reaction_namespace='sbgn.process',
    **meta,
):
    """Decode SBGN-ML into glyphs + reaction events for database.diagram_modules."""
    path = Path(path)
    raw = parse(path)
    glyph_ids = set(raw['glyphs'])
    port_to_owner = raw['port_to_owner']
    audit = Counter()
    glyphs = {}
    for gid, g in raw['glyphs'].items():
        if g['kind'] != 'entity':
            continue
        glyphs[gid] = dict(
            key=gid,
            entity=g['entity'],
            label=g['label'] or g['entity'],
            xy=list(g['xy']),
            bounds=g['bounds'],
            sbgn_class=g['class_'],
            identifiers=g['identifiers'],
            observed=True,
            position_origin='observed',
        )

    # Group arcs by owning process
    by_process = defaultdict(lambda: dict(substrates=[], products=[], modifiers=[], other=[]))
    for arc in raw['arcs']:
        src = _port_owner(arc['source'], glyph_ids, port_to_owner)
        tgt = _port_owner(arc['target'], glyph_ids, port_to_owner)
        cls = arc['class_']
        if cls in SUBSTRATE_ARCS:
            # consumption: entity -> process
            proc, ent = tgt, src
            if raw['glyphs'].get(proc, {}).get('kind') != 'process':
                audit['consumption_target_not_process'] += 1
                continue
            by_process[proc]['substrates'].append((ent, arc))
        elif cls in PRODUCT_ARCS:
            proc, ent = src, tgt
            if raw['glyphs'].get(proc, {}).get('kind') != 'process':
                audit['production_source_not_process'] += 1
                continue
            by_process[proc]['products'].append((ent, arc))
        elif cls in MODIFIER_ARCS:
            # usually entity -> process
            if raw['glyphs'].get(tgt, {}).get('kind') == 'process':
                by_process[tgt]['modifiers'].append((src, arc))
            elif raw['glyphs'].get(src, {}).get('kind') == 'process':
                by_process[src]['modifiers'].append((tgt, arc))
            else:
                audit['modifier_without_process'] += 1
        else:
            audit['arc_unhandled:' + (cls or '?')] += 1

    events = []
    excluded = []
    for pid, g in sorted(raw['glyphs'].items()):
        if g['kind'] != 'process':
            continue
        bag = by_process.get(pid)
        if not bag or (not bag['substrates'] and not bag['products']):
            excluded.append(dict(id=pid, reason='process_without_consumption_or_production'))
            audit['process_without_stoichiometry_arcs'] += 1
            continue
        visible = {}
        coefficients = {}
        ok = True
        for ent, arc in bag['substrates']:
            if ent not in glyphs:
                ok = False
                audit['substrate_not_entity_glyph'] += 1
                continue
            entity = glyphs[ent]['entity']
            visible[ent] = visible.get(ent, 0) - 1
            coefficients[entity] = coefficients.get(entity, 0) - 1
        for ent, arc in bag['products']:
            if ent not in glyphs:
                ok = False
                audit['product_not_entity_glyph'] += 1
                continue
            entity = glyphs[ent]['entity']
            visible[ent] = visible.get(ent, 0) + 1
            coefficients[entity] = coefficients.get(entity, 0) + 1
        if not visible:
            excluded.append(dict(id=pid, reason='no_visible_entity_participants'))
            audit['process_no_visible_entities'] += 1
            continue
        # Drop cancelled entities
        coefficients = {k: v for k, v in coefficients.items() if v}
        visible = {k: v for k, v in visible.items() if v}
        if not visible:
            excluded.append(dict(id=pid, reason='participants_cancel'))
            audit['process_participants_cancel'] += 1
            continue
        rid = g['label'] or pid
        refs = {reaction_namespace: [rid]}
        # Promote Rhea / EC identifiers on the process glyph when present
        for ident in g.get('identifiers') or []:
            if ident.upper().startswith('RHEA:'):
                refs.setdefault('rhea', []).append(ident.split(':', 1)[-1])
            elif ident.lower().startswith('ec-code:'):
                refs.setdefault('ec-code', []).append(ident.split(':', 1)[-1])
        path_pts = []
        for ent, arc in bag['substrates'] + bag['products']:
            if arc.get('points'):
                path_pts.extend(arc['points'])
        events.append(dict(
            key='r:' + pid,
            id=rid,
            source_reaction_id=pid,
            references=refs,
            coefficients=coefficients,
            visible=visible,
            xy=list(g['xy']),
            position_origin='observed',
            reversible=False,
            reversibility_origin='sbgn_unspecified',
            core=list(visible),
            path_centerline=path_pts or None,
            modifiers=[dict(glyph=e, arc_class=a['class_'], arc_id=a.get('id')) for e, a in bag['modifiers']],
            sbgn_class=g['class_'],
            stoichiometry_origin='sbgn_arc_unit_coefficient',
        ))
        if ok:
            audit['decoded_process'] += 1
        else:
            audit['decoded_process_partial'] += 1

    digest = sha(path.read_bytes())
    out = dict(
        id=uid([source, digest]),
        source=source,
        species=species or 'unspecified',
        name=name or path.stem,
        lineage=lineage or (source + ':' + digest),
        source_sha256=digest,
        path=str(path.relative_to(BASE.parent)) if path.is_absolute() and BASE.parent in path.parents else str(path),
        license=license,
        partition=partition,
        admitted_to_open_core=admitted_to_open_core,
        format='sbgn-ml',
        language=raw['language'],
        map_id=raw['map_id'] or path.stem,
        glyphs=glyphs,
        events=events,
        audit=dict(audit),
        excluded=excluded,
        raw_glyph_count=len(raw['glyphs']),
        raw_arc_count=len(raw['arcs']),
        parse_issues=raw['issues'],
    )
    out.update(meta)
    return out
