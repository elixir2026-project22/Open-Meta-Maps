"""Export a composed Escher map.json to SBGN-ML, SBML Layout, and CellDesigner SBML.

Exports read the final map geometry and the input GEM. They do not change layout
coordinates or the model stoichiometry. Shared helpers live here so CLI compose
and the Map Generator write the same files.
"""
from __future__ import annotations

import math
import re
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
from xml.dom import minidom

import numpy as np

from core import DATA, comp, read
from model_io import load_model

SBGN_NS = 'http://sbgn.org/libsbgn/0.3'
SBGN_SCHEMA = DATA / 'sbgn_schema' / 'SBGN.xsd'
CD_NS = 'http://www.sbml.org/2001/ns/celldesigner'
SBML_L2_NS = 'http://www.sbml.org/sbml/level2/version4'

# Glyph sizes (SBGN Process Description). Centres match Escher node x/y.
PROCESS_W = PROCESS_H = 20.0
CHEM_PRIMARY = 40.0
CHEM_SECONDARY = 28.0
MACRO_W, MACRO_H = 80.0, 40.0
PORT_INSET = 0.0

# SBO → SBGN entity class. Default is simple chemical.
SBO_MACRO = {245, 252, 253, 254, 297, 336}  # macromolecule / protein family
SBO_SIMPLE = {247, 249, 250}  # simple chemical / ion / ribonucleic acid etc.


def _body(data):
    return data[1] if isinstance(data, list) else data


def _header(data):
    return data[0] if isinstance(data, list) else {}


def _xy(node):
    return float(node['x']), float(node['y'])


def _bezier_points(p, c1, c2, q, n=12):
    t = np.linspace(0.0, 1.0, n + 1)[:, None]
    pts = (1 - t) ** 3 * p + 3 * (1 - t) ** 2 * t * c1 + 3 * (1 - t) * t ** 2 * c2 + t ** 3 * q
    return pts.tolist()


def _segment_controls(reaction, a, b):
    for seg in reaction['segments'].values():
        if seg['from_node_id'] == a and seg['to_node_id'] == b:
            return seg.get('b1'), seg.get('b2')
        if seg['from_node_id'] == b and seg['to_node_id'] == a:
            return seg.get('b2'), seg.get('b1')
    return None, None


def _path_polyline(reaction, path, nodes, samples_per_curve=12):
    """Flatten Escher segments along node path into a polyline (includes endpoints)."""
    pts = []
    for a, b in zip(path, path[1:]):
        ax, ay = _xy(nodes[a])
        bx, by = _xy(nodes[b])
        b1, b2 = _segment_controls(reaction, a, b)
        if b1 and b2:
            curve = _bezier_points(
                np.array([ax, ay], float),
                np.array([b1['x'], b1['y']], float),
                np.array([b2['x'], b2['y']], float),
                np.array([bx, by], float),
                n=samples_per_curve,
            )
            pts.extend(curve if not pts else curve[1:])
        else:
            if not pts:
                pts.append([ax, ay])
            pts.append([bx, by])
    return pts


def _reaction_paths(reaction, nodes):
    """Polylines from midmarker to each metabolite, following stored segments."""
    out_edges = {}
    for seg in reaction['segments'].values():
        out_edges.setdefault(seg['from_node_id'], []).append(seg['to_node_id'])
        out_edges.setdefault(seg['to_node_id'], [])
    targets = {t for ts in out_edges.values() for t in ts}
    roots = [n for n in out_edges if n not in targets] or list(out_edges)[:1]
    lines = []

    def walk(node, path):
        nxt = out_edges.get(node, [])
        if not nxt:
            lines.append(path)
            return
        for n in nxt:
            if n in nodes:
                walk(n, path + [n])

    for root in roots:
        if root in nodes:
            walk(root, [root])
    return lines


def _midmarker(reaction, nodes):
    for seg in reaction['segments'].values():
        for nid in (seg['from_node_id'], seg['to_node_id']):
            if nodes.get(nid, {}).get('node_type') == 'midmarker':
                return nid
    return None


def _model_lookup(model):
    mets = {m['id']: m for m in model['metabolites']}
    rxns = {r['id']: r for r in model['reactions']}
    return mets, rxns


def _sbo_term(model_met, sbml_species=None):
    if sbml_species is not None and sbml_species.isSetSBOTerm():
        return sbml_species.getSBOTerm()
    for a in model_met.get('annotations') or []:
        uri = a.get('uri') or ''
        m = re.search(r'SBO[:_](\d+)', uri, re.I)
        if m:
            return int(m.group(1))
    return None


def _entity_class(model_met, sbml_species=None):
    sbo = _sbo_term(model_met, sbml_species)
    if sbo in SBO_MACRO:
        return 'macromolecule'
    if sbo in SBO_SIMPLE or sbo is None:
        return 'simple chemical'
    # Unknown SBO: prefer simple chemical for metabolites.
    return 'simple chemical'


def _glyph_size(entity_class, primary):
    if entity_class == 'macromolecule':
        return MACRO_W, MACRO_H
    side = CHEM_PRIMARY if primary else CHEM_SECONDARY
    return side, side


def _bbox(cx, cy, w, h):
    return dict(x=cx - w / 2.0, y=cy - h / 2.0, w=w, h=h)


def _safe_id(prefix, raw):
    """Build an XML NCName/ID. Prefix supplies the leading letter when raw starts with a digit."""
    text = re.sub(r'[^A-Za-z0-9_\-.]', '_', str(raw))
    full = f'{prefix}{text}'
    if not full or not re.match(r'^[A-Za-z_]', full):
        full = 'n' + full
    return full


def _sbml_id_map(doc):
    """Map decoded BiGG-style ids and raw SBML ids to SBML element ids."""
    model = doc.getModel()
    species = {}
    reactions = {}
    compartments = {}
    for i in range(model.getNumSpecies()):
        sp = model.getSpecies(i)
        sid = sp.getId()
        species[sid] = sid
        decoded = sid[2:] if sid.startswith('M_') else sid
        species.setdefault(decoded, sid)
    for i in range(model.getNumReactions()):
        rx = model.getReaction(i)
        rid = rx.getId()
        reactions[rid] = rid
        decoded = rid[2:] if rid.startswith('R_') else rid
        reactions.setdefault(decoded, rid)
    for i in range(model.getNumCompartments()):
        c = model.getCompartment(i)
        compartments[c.getId()] = c.getId()
    return species, reactions, compartments


def _register_namespaces():
    ET.register_namespace('', SBGN_NS)


def validate_sbgn_xsd(path):
    """Validate SBGN-ML against the vendored libSBGN 0.3 XSD. Returns (ok, errors)."""
    from lxml import etree
    if not SBGN_SCHEMA.is_file():
        return False, [f'missing schema: {SBGN_SCHEMA}']
    schema = etree.XMLSchema(etree.parse(str(SBGN_SCHEMA)))
    doc = etree.parse(str(path))
    ok = schema.validate(doc)
    errors = [str(e) for e in schema.error_log.filter_from_errors()]
    return ok, errors


def validate_sbml_libsbml(path):
    """Validate SBML with libSBML. Returns (ok, error_messages) — warnings are ignored."""
    import libsbml
    doc = libsbml.readSBMLFromFile(str(path))
    messages = []
    for i in range(doc.getNumErrors()):
        err = doc.getError(i)
        if err.getSeverity() >= libsbml.LIBSBML_SEV_ERROR:
            messages.append(err.getMessage())
    doc.checkConsistency()
    for i in range(doc.getNumErrors()):
        err = doc.getError(i)
        if err.getSeverity() >= libsbml.LIBSBML_SEV_ERROR:
            messages.append(err.getMessage())
    return not messages, messages


# ------------------------------------------------------------------------------------------------ SBGN-ML


def build_sbgn(data, model):
    """Build an SBGN-ML 0.3 Process Description document from Escher map + model."""
    _register_namespaces()
    body = _body(data)
    nodes = body['nodes']
    reactions = body['reactions']
    canvas = body.get('canvas') or {}
    mets, _ = _model_lookup(model)

    # Count clones per metabolite id
    clone_count = defaultdict(int)
    for n in nodes.values():
        if n.get('node_type') == 'metabolite' and n.get('bigg_id'):
            clone_count[n['bigg_id']] += 1

    root = ET.Element(f'{{{SBGN_NS}}}sbgn')
    map_el = ET.SubElement(root, f'{{{SBGN_NS}}}map', {
        'id': 'map1',
        'language': 'process description',
    })
    # Map bbox from canvas
    ET.SubElement(map_el, f'{{{SBGN_NS}}}bbox', {
        'x': f"{float(canvas.get('x', 0)):.4f}",
        'y': f"{float(canvas.get('y', 0)):.4f}",
        'w': f"{float(canvas.get('width', 1000)):.4f}",
        'h': f"{float(canvas.get('height', 1000)):.4f}",
    })

    # Compartment glyphs from metabolite extents
    by_comp = defaultdict(list)
    for n in nodes.values():
        if n.get('node_type') != 'metabolite':
            continue
        mid = n.get('bigg_id')
        c = (mets.get(mid) or {}).get('compartment') or comp(mid or '')
        by_comp[c].append(_xy(n))
    pad = 40.0
    for cid, pts in sorted(by_comp.items()):
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        x0, y0 = min(xs) - pad, min(ys) - pad
        w, h = max(xs) - min(xs) + 2 * pad, max(ys) - min(ys) + 2 * pad
        gid = _safe_id('comp_', cid)
        g = ET.SubElement(map_el, f'{{{SBGN_NS}}}glyph', {'id': gid, 'class': 'compartment', 'compartmentOrder': '1'})
        ET.SubElement(g, f'{{{SBGN_NS}}}label', {'text': cid})
        ET.SubElement(g, f'{{{SBGN_NS}}}bbox', {
            'x': f'{x0:.4f}', 'y': f'{y0:.4f}', 'w': f'{max(w, 80):.4f}', 'h': f'{max(h, 80):.4f}',
        })

    entity_glyph = {}  # escher node id -> sbgn glyph id
    drawn_mids = set()
    for nid, n in nodes.items():
        if n.get('node_type') != 'metabolite':
            continue
        mid = n.get('bigg_id') or nid
        drawn_mids.add(mid)
        met = mets.get(mid) or dict(id=mid, name=n.get('name') or mid, compartment=comp(mid))
        cls = _entity_class(met)
        primary = bool(n.get('node_is_primary', True))
        w, h = _glyph_size(cls, primary)
        cx, cy = _xy(n)
        bb = _bbox(cx, cy, w, h)
        gid = _safe_id('ep_', nid)
        entity_glyph[nid] = gid
        compartment = (met.get('compartment') or comp(mid) or '')
        attrs = {'id': gid, 'class': cls}
        if compartment and _safe_id('comp_', compartment) in {_safe_id('comp_', c) for c in by_comp}:
            attrs['compartmentRef'] = _safe_id('comp_', compartment)
        g = ET.SubElement(map_el, f'{{{SBGN_NS}}}glyph', attrs)
        ET.SubElement(g, f'{{{SBGN_NS}}}label', {'text': mid})
        if clone_count.get(mid, 0) > 1:
            ET.SubElement(g, f'{{{SBGN_NS}}}clone')
        ET.SubElement(g, f'{{{SBGN_NS}}}bbox', {
            'x': f"{bb['x']:.4f}", 'y': f"{bb['y']:.4f}", 'w': f"{bb['w']:.4f}", 'h': f"{bb['h']:.4f}",
        })

    # Every model species appears at least once (park undrawn species below the canvas).
    canvas_w = float(canvas.get('width', 1000))
    canvas_h = float(canvas.get('height', 1000))
    park_x, park_y = 40.0, canvas_h + 60.0
    parked = False
    for mid, met in mets.items():
        if mid in drawn_mids:
            continue
        parked = True
        cls = _entity_class(met)
        w, h = _glyph_size(cls, False)
        gid = _safe_id('ep_park_', mid)
        g = ET.SubElement(map_el, f'{{{SBGN_NS}}}glyph', {'id': gid, 'class': cls})
        ET.SubElement(g, f'{{{SBGN_NS}}}label', {'text': mid})
        ET.SubElement(g, f'{{{SBGN_NS}}}bbox', {
            'x': f'{park_x:.4f}', 'y': f'{park_y:.4f}', 'w': f'{w:.4f}', 'h': f'{h:.4f}',
        })
        park_x += w + 10
        if park_x > canvas_w:
            park_x = 40.0
            park_y += h + 10
    if parked:
        map_el.find(f'{{{SBGN_NS}}}bbox').set('h', f'{park_y + CHEM_PRIMARY + 20:.4f}')

    # Index metabolite nodes by bigg_id for attaching reaction participants
    by_bigg = defaultdict(list)
    for nid, n in nodes.items():
        if n.get('node_type') == 'metabolite' and n.get('bigg_id'):
            by_bigg[n['bigg_id']].append(nid)

    # XSD: all glyphs, then all arcs (never interleave).
    pending_arcs = []
    for rk, rx in reactions.items():
        rid = rx.get('bigg_id') or rk
        mid_id = _midmarker(rx, nodes)
        if mid_id is None:
            continue
        mx, my = _xy(nodes[mid_id])
        bb = _bbox(mx, my, PROCESS_W, PROCESS_H)
        gid = _safe_id('re_', rid)
        g = ET.SubElement(map_el, f'{{{SBGN_NS}}}glyph', {'id': gid, 'class': 'process'})
        ET.SubElement(g, f'{{{SBGN_NS}}}bbox', {
            'x': f"{bb['x']:.4f}", 'y': f"{bb['y']:.4f}", 'w': f"{bb['w']:.4f}", 'h': f"{bb['h']:.4f}",
        })
        port_in = f'{gid}.1'
        port_out = f'{gid}.2'
        ET.SubElement(g, f'{{{SBGN_NS}}}port', {
            'id': port_in,
            'x': f"{bb['x'] - PORT_INSET:.4f}",
            'y': f"{my:.4f}",
        })
        ET.SubElement(g, f'{{{SBGN_NS}}}port', {
            'id': port_out,
            'x': f"{bb['x'] + bb['w'] + PORT_INSET:.4f}",
            'y': f"{my:.4f}",
        })

        coeffs = {m['bigg_id']: float(m['coefficient']) for m in rx.get('metabolites') or []}
        used_nodes = set()
        for path in _reaction_paths(rx, nodes):
            end = nodes.get(path[-1])
            if not end or end.get('node_type') != 'metabolite':
                continue
            met_nid = path[-1]
            if met_nid in used_nodes:
                continue
            used_nodes.add(met_nid)
            mid = end.get('bigg_id')
            coef = coeffs.get(mid, 0.0)
            poly = _path_polyline(rx, path, nodes)
            if len(poly) < 2:
                continue
            eg = entity_glyph.get(met_nid)
            if not eg:
                continue
            if coef < 0:
                pts = list(reversed(poly))
                source, target, cls = eg, port_in, 'consumption'
            else:
                pts = poly
                source, target, cls = port_out, eg, 'production'
            pending_arcs.append((cls, source, target, pts))

        attached = {nodes[n].get('bigg_id') for n in used_nodes}
        for mid, coef in coeffs.items():
            if mid in attached:
                continue
            candidates = [n for n in by_bigg.get(mid, []) if n in entity_glyph]
            if not candidates:
                continue
            met_nid = candidates[0]
            eg = entity_glyph[met_nid]
            ex, ey = _xy(nodes[met_nid])
            if coef < 0:
                source, target, cls = eg, port_in, 'consumption'
                pts = [[ex, ey], [bb['x'], my]]
            else:
                source, target, cls = port_out, eg, 'production'
                pts = [[bb['x'] + bb['w'], my], [ex, ey]]
            pending_arcs.append((cls, source, target, pts))

    for arc_i, (cls, source, target, pts) in enumerate(pending_arcs):
        arc = ET.SubElement(map_el, f'{{{SBGN_NS}}}arc', {
            'id': f'arc{arc_i}',
            'class': cls,
            'source': source,
            'target': target,
        })
        ET.SubElement(arc, f'{{{SBGN_NS}}}start', {'x': f"{pts[0][0]:.4f}", 'y': f"{pts[0][1]:.4f}"})
        for p in pts[1:-1]:
            ET.SubElement(arc, f'{{{SBGN_NS}}}next', {'x': f"{p[0]:.4f}", 'y': f"{p[1]:.4f}"})
        ET.SubElement(arc, f'{{{SBGN_NS}}}end', {'x': f"{pts[-1][0]:.4f}", 'y': f"{pts[-1][1]:.4f}"})

    return root


def write_sbgn(path, data, model):
    root = build_sbgn(data, model)
    tree = ET.ElementTree(root)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tree.write(path, encoding='utf-8', xml_declaration=True)
    return path


# ------------------------------------------------------------------------------------------------ SBML Layout


def _set_bb(go, x, y, w, h, bb_id):
    bb = go.getBoundingBox()
    bb.setId(bb_id)
    bb.setX(float(x))
    bb.setY(float(y))
    bb.setWidth(float(w))
    bb.setHeight(float(h))


def build_sbml_layout(model_path, data, model=None):
    """Return an SBML document: input model unchanged + Layout (+ Render) package."""
    import libsbml
    model_path = Path(model_path)
    doc = libsbml.readSBMLFromFile(str(model_path))
    if doc.getNumErrors() and any(
        doc.getError(i).getSeverity() >= libsbml.LIBSBML_SEV_ERROR for i in range(doc.getNumErrors())
    ):
        raise ValueError(f'cannot read SBML model: {model_path}')
    if model is None:
        model = load_model(model_path)

    xmlns = libsbml.LayoutExtension.getXmlnsL3V1V1()
    doc.enablePackage(xmlns, 'layout', True)
    doc.setPackageRequired('layout', False)
    try:
        doc.enablePackage(libsbml.RenderExtension.getXmlnsL3V1V1(), 'render', True)
        doc.setPackageRequired('render', False)
        render_ok = True
    except Exception:
        render_ok = False

    sbml_model = doc.getModel()
    # Remove any previous layout package content so re-export is deterministic
    plugin = sbml_model.getPlugin('layout')
    while plugin.getNumLayouts() > 0:
        plugin.removeLayout(0)

    body = _body(data)
    nodes = body['nodes']
    reactions = body['reactions']
    canvas = body.get('canvas') or {}
    mets, _ = _model_lookup(model)
    sp_map, rx_map, comp_map = _sbml_id_map(doc)

    layout = plugin.createLayout()
    layout.setId('Layout_1')
    layout.setName('Open Meta Maps layout')
    dims = layout.getDimensions()
    dims.setWidth(float(canvas.get('width', 1000)))
    dims.setHeight(float(canvas.get('height', 1000)))

    # Compartment glyphs
    by_comp = defaultdict(list)
    for n in nodes.values():
        if n.get('node_type') != 'metabolite':
            continue
        mid = n.get('bigg_id')
        c = (mets.get(mid) or {}).get('compartment') or comp(mid or '')
        by_comp[c].append(_xy(n))
    pad = 40.0
    for cid, pts in sorted(by_comp.items()):
        sbml_cid = comp_map.get(cid)
        if not sbml_cid:
            continue
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        x0, y0 = min(xs) - pad, min(ys) - pad
        w, h = max(xs) - min(xs) + 2 * pad, max(ys) - min(ys) + 2 * pad
        cg = layout.createCompartmentGlyph()
        cg.setId(_safe_id('cg_', cid))
        cg.setCompartmentId(sbml_cid)
        _set_bb(cg, x0, y0, max(w, 80), max(h, 80), _safe_id('bb_cg_', cid))

    # Species glyphs: one per Escher metabolite node (clones → multiple glyphs)
    species_glyph_for_node = {}
    glyph_index = defaultdict(int)
    for nid, n in nodes.items():
        if n.get('node_type') != 'metabolite':
            continue
        mid = n.get('bigg_id')
        sbml_sid = sp_map.get(mid)
        if not sbml_sid:
            continue
        met = mets.get(mid) or {}
        cls = _entity_class(met)
        primary = bool(n.get('node_is_primary', True))
        w, h = _glyph_size(cls, primary)
        cx, cy = _xy(n)
        idx = glyph_index[sbml_sid]
        glyph_index[sbml_sid] += 1
        sgid = _safe_id('sg_', f'{mid}_{idx}')
        sg = layout.createSpeciesGlyph()
        sg.setId(sgid)
        sg.setSpeciesId(sbml_sid)
        _set_bb(sg, cx - w / 2, cy - h / 2, w, h, _safe_id('bb_', sgid))
        species_glyph_for_node[nid] = sgid
        # Text glyph for label
        tg = layout.createTextGlyph()
        tg.setId(_safe_id('tg_', sgid))
        tg.setGraphicalObjectId(sgid)
        label = n.get('name') or mid
        tg.setText(label)
        lx = float(n.get('label_x', cx))
        ly = float(n.get('label_y', cy - h / 2 - 8))
        _set_bb(tg, lx, ly - 6, max(len(label) * 6.0, 20), 12, _safe_id('bb_tg_', sgid))

    # Species without a drawn node still get a parked glyph so every species appears
    placed_species = {sp_map[n['bigg_id']] for n in nodes.values()
                      if n.get('node_type') == 'metabolite' and n.get('bigg_id') in sp_map}
    canvas_w = float(canvas.get('width', 1000))
    park_y = float(canvas.get('height', 1000)) + 40
    park_x = 40.0
    for mid, met in mets.items():
        sbml_sid = sp_map.get(mid)
        if not sbml_sid or sbml_sid in placed_species:
            continue
        sg = layout.createSpeciesGlyph()
        sgid = _safe_id('sg_park_', mid)
        sg.setId(sgid)
        sg.setSpeciesId(sbml_sid)
        _set_bb(sg, park_x, park_y, CHEM_SECONDARY, CHEM_SECONDARY, _safe_id('bb_', sgid))
        park_x += CHEM_SECONDARY + 10
        if park_x > canvas_w:
            park_x = 40.0
            park_y += CHEM_SECONDARY + 10

    by_bigg_nodes = defaultdict(list)
    for nid, n in nodes.items():
        if n.get('node_type') == 'metabolite' and n.get('bigg_id'):
            by_bigg_nodes[n['bigg_id']].append(nid)

    for rk, rx in reactions.items():
        rid = rx.get('bigg_id') or rk
        sbml_rid = rx_map.get(rid)
        if not sbml_rid:
            continue
        mid_id = _midmarker(rx, nodes)
        if mid_id is None:
            continue
        mx, my = _xy(nodes[mid_id])
        rg = layout.createReactionGlyph()
        rg.setId(_safe_id('rg_', rid))
        rg.setReactionId(sbml_rid)
        _set_bb(rg, mx - PROCESS_W / 2, my - PROCESS_H / 2, PROCESS_W, PROCESS_H, _safe_id('bb_rg_', rid))

        coeffs = {m['bigg_id']: float(m['coefficient']) for m in rx.get('metabolites') or []}
        used = set()
        srg_i = 0
        for path in _reaction_paths(rx, nodes):
            end = nodes.get(path[-1])
            if not end or end.get('node_type') != 'metabolite':
                continue
            met_nid = path[-1]
            if met_nid in used:
                continue
            used.add(met_nid)
            mid = end.get('bigg_id')
            sgid = species_glyph_for_node.get(met_nid)
            if not sgid:
                continue
            coef = coeffs.get(mid, 0.0)
            role = libsbml.SPECIES_ROLE_SUBSTRATE if coef < 0 else libsbml.SPECIES_ROLE_PRODUCT
            srg = rg.createSpeciesReferenceGlyph()
            srg.setId(_safe_id('srg_', f'{rid}_{srg_i}'))
            srg_i += 1
            srg.setSpeciesGlyphId(sgid)
            srg.setRole(role)
            poly = _path_polyline(rx, path, nodes)
            if coef < 0:
                poly = list(reversed(poly))
            if len(poly) >= 2:
                for a, b in zip(poly, poly[1:]):
                    ls = srg.createLineSegment()
                    ls.getStart().setX(float(a[0]))
                    ls.getStart().setY(float(a[1]))
                    ls.getEnd().setX(float(b[0]))
                    ls.getEnd().setY(float(b[1]))

        # Cover stoichiometric participants missing from the segment walk
        attached = {nodes[n].get('bigg_id') for n in used}
        for mid, coef in coeffs.items():
            if mid in attached:
                continue
            candidates = [n for n in by_bigg_nodes.get(mid, []) if n in species_glyph_for_node]
            if not candidates:
                continue
            met_nid = candidates[0]
            sgid = species_glyph_for_node[met_nid]
            ex, ey = _xy(nodes[met_nid])
            srg = rg.createSpeciesReferenceGlyph()
            srg.setId(_safe_id('srg_', f'{rid}_{srg_i}'))
            srg_i += 1
            srg.setSpeciesGlyphId(sgid)
            srg.setRole(libsbml.SPECIES_ROLE_SUBSTRATE if coef < 0 else libsbml.SPECIES_ROLE_PRODUCT)
            ls = srg.createLineSegment()
            if coef < 0:
                ls.getStart().setX(ex); ls.getStart().setY(ey)
                ls.getEnd().setX(mx); ls.getEnd().setY(my)
            else:
                ls.getStart().setX(mx); ls.getStart().setY(my)
                ls.getEnd().setX(ex); ls.getEnd().setY(ey)

    # Reactions with no midmarker still get a reaction glyph so coverage is complete
    drawn_rx = {rx_map[rx['bigg_id']] for rx in reactions.values() if rx.get('bigg_id') in rx_map}
    park_rx_x, park_rx_y = 40.0, park_y + 60
    for r in model['reactions']:
        sbml_rid = rx_map.get(r['id'])
        if not sbml_rid or sbml_rid in drawn_rx:
            continue
        rg = layout.createReactionGlyph()
        rg.setId(_safe_id('rg_park_', r['id']))
        rg.setReactionId(sbml_rid)
        _set_bb(rg, park_rx_x, park_rx_y, PROCESS_W, PROCESS_H, _safe_id('bb_rg_park_', r['id']))
        park_rx_x += PROCESS_W + 10
        if park_rx_x > canvas_w:
            park_rx_x = 40.0
            park_rx_y += PROCESS_H + 10

    if render_ok:
        _add_basic_render(layout)

    return doc


def _add_basic_render(layout):
    """Attach a minimal local render information to the layout when the render package is enabled."""
    import libsbml
    rplugin = layout.getPlugin('render')
    if rplugin is None:
        return
    info = rplugin.createLocalRenderInformation()
    info.setId('OpenMetaMaps_render')
    info.setName('Open Meta Maps')
    info.setProgramName('Open Meta Maps')
    try:
        style = info.createStyle('speciesStyle')
        style.getGroup().setStroke('black')
        style.getGroup().setStrokeWidth(1.0)
        style.getGroup().setFill('#dcebdf')
    except Exception:
        pass


def write_sbml_layout(path, model_path, data, model=None):
    import libsbml
    doc = build_sbml_layout(model_path, data, model=model)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ok = libsbml.writeSBMLToFile(doc, str(path))
    if not ok:
        raise RuntimeError(f'libsbml failed to write {path}')
    return path


# ------------------------------------------------------------------------------------------------ CellDesigner SBML


def _global_to_local(start, end, points):
    """Inverse of celldesigner.local_to_global for editPoints."""
    a = np.array(start, float)
    v = np.array(end, float) - a
    if float(np.linalg.norm(v)) < 1e-8:
        return []
    perp = np.array([-v[1], v[0]], float)
    # Solve [v | perp] @ [x, y] = p - a
    mat = np.column_stack([v, perp])
    out = []
    for p in points:
        try:
            xy = np.linalg.lstsq(mat, np.array(p, float) - a, rcond=None)[0]
            out.append((float(xy[0]), float(xy[1])))
        except Exception:
            continue
    return out


def _cd(tag):
    return f'{{{CD_NS}}}{tag}'


def _sbml2(tag):
    return f'{{{SBML_L2_NS}}}{tag}'


def build_celldesigner(data, model):
    """Build CellDesigner 4.x SBML (L2V4) with speciesAlias bounds and reaction editPoints."""
    ET.register_namespace('celldesigner', CD_NS)
    ET.register_namespace('', SBML_L2_NS)
    body = _body(data)
    nodes = body['nodes']
    reactions = body['reactions']
    canvas = body.get('canvas') or {}
    mets, rxns = _model_lookup(model)

    root = ET.Element(_sbml2('sbml'), {'level': '2', 'version': '4'})
    model_el = ET.SubElement(root, _sbml2('model'), {
        'id': _safe_id('model_', model.get('id') or 'map'),
        'name': model.get('name') or model.get('id') or 'map',
    })

    # Model annotation / CellDesigner extension
    ann = ET.SubElement(model_el, _sbml2('annotation'))
    ext = ET.SubElement(ann, _cd('extension'))
    ET.SubElement(ext, _cd('modelVersion')).text = '4.0'
    ET.SubElement(ext, _cd('modelDisplay'), {
        'sizeX': str(int(float(canvas.get('width', 1000)))),
        'sizeY': str(int(float(canvas.get('height', 1000)))),
    })

    # Compartments
    comps = sorted({m.get('compartment') or 'c' for m in model['metabolites']})
    loc = ET.SubElement(model_el, _sbml2('listOfCompartments'))
    list_ca = ET.SubElement(ext, _cd('listOfCompartmentAliases'))
    for i, cid in enumerate(comps):
        ET.SubElement(loc, _sbml2('compartment'), {'id': cid, 'name': cid, 'size': '1'})
        # bounds from metabolites
        pts = [_xy(n) for n in nodes.values()
               if n.get('node_type') == 'metabolite' and (mets.get(n.get('bigg_id')) or {}).get('compartment', comp(n.get('bigg_id') or '')) == cid]
        if not pts:
            pts = [(0.0, 0.0), (100.0, 100.0)]
        pad = 40.0
        xs, ys = [p[0] for p in pts], [p[1] for p in pts]
        x0, y0 = min(xs) - pad, min(ys) - pad
        w, h = max(xs) - min(xs) + 2 * pad, max(ys) - min(ys) + 2 * pad
        ca = ET.SubElement(list_ca, _cd('compartmentAlias'), {'id': f'ca{i}', 'compartment': cid})
        ET.SubElement(ca, _cd('bounds'), {
            'x': f'{x0:.2f}', 'y': f'{y0:.2f}', 'w': f'{max(w, 80):.2f}', 'h': f'{max(h, 80):.2f}',
        })

    ET.SubElement(ext, _cd('listOfComplexSpeciesAliases'))
    list_sa = ET.SubElement(ext, _cd('listOfSpeciesAliases'))
    ET.SubElement(ext, _cd('listOfGroups'))
    ET.SubElement(ext, _cd('listOfProteins'))
    ET.SubElement(ext, _cd('listOfGenes'))
    ET.SubElement(ext, _cd('listOfRNAs'))
    ET.SubElement(ext, _cd('listOfAntisenseRNAs'))
    ET.SubElement(ext, _cd('listOfLayers'))
    ET.SubElement(ext, _cd('listOfBlockDiagrams'))

    # Species + aliases (one alias per Escher metabolite node; ensure every model species appears)
    los = ET.SubElement(model_el, _sbml2('listOfSpecies'))
    species_el = {}
    for mid, met in mets.items():
        sid = _safe_id('s_', mid)
        species_el[mid] = sid
        ET.SubElement(los, _sbml2('species'), {
            'id': sid,
            'name': met.get('name') or mid,
            'compartment': met.get('compartment') or comp(mid) or comps[0],
            'initialAmount': '0',
        })

    alias_for_node = {}
    alias_by_mid = defaultdict(list)
    alias_i = 0
    for nid, n in nodes.items():
        if n.get('node_type') != 'metabolite':
            continue
        mid = n.get('bigg_id')
        if mid not in species_el:
            continue
        met = mets[mid]
        cls = _entity_class(met)
        primary = bool(n.get('node_is_primary', True))
        w, h = _glyph_size(cls, primary)
        cx, cy = _xy(n)
        aid = f'a{alias_i}'
        alias_i += 1
        alias_for_node[nid] = aid
        alias_by_mid[mid].append(aid)
        sa = ET.SubElement(list_sa, _cd('speciesAlias'), {'id': aid, 'species': species_el[mid]})
        ET.SubElement(sa, _cd('activity')).text = 'inactive'
        ET.SubElement(sa, _cd('bounds'), {
            'x': f'{cx - w / 2:.2f}', 'y': f'{cy - h / 2:.2f}', 'w': f'{w:.2f}', 'h': f'{h:.2f}',
        })
        ET.SubElement(sa, _cd('font'), {'size': '12'})
        ET.SubElement(sa, _cd('view'), {'state': 'usual'})
        usual = ET.SubElement(sa, _cd('usualView'))
        ET.SubElement(usual, _cd('innerPosition'), {'x': '0.0', 'y': '0.0'})
        ET.SubElement(usual, _cd('boxSize'), {'width': f'{w:.2f}', 'height': f'{h:.2f}'})
        ET.SubElement(usual, _cd('singleLine'), {'width': '1.0'})
        ET.SubElement(usual, _cd('paint'), {'color': 'ffccff66' if primary else 'ffffffff', 'scheme': 'Color'})

    # Park aliases for species never drawn
    park_x, park_y = 40.0, float(canvas.get('height', 1000)) + 40
    for mid in mets:
        if alias_by_mid.get(mid):
            continue
        aid = f'a{alias_i}'
        alias_i += 1
        alias_by_mid[mid].append(aid)
        sa = ET.SubElement(list_sa, _cd('speciesAlias'), {'id': aid, 'species': species_el[mid]})
        ET.SubElement(sa, _cd('activity')).text = 'inactive'
        ET.SubElement(sa, _cd('bounds'), {
            'x': f'{park_x:.2f}', 'y': f'{park_y:.2f}', 'w': f'{CHEM_SECONDARY:.2f}', 'h': f'{CHEM_SECONDARY:.2f}',
        })
        park_x += CHEM_SECONDARY + 10

    # Reactions
    lor = ET.SubElement(model_el, _sbml2('listOfReactions'))
    # Index drawn reactions
    drawn = {rx.get('bigg_id'): (rk, rx) for rk, rx in reactions.items()}

    def alias_xy(aid):
        for nid, a in alias_for_node.items():
            if a == aid:
                return _xy(nodes[nid])
        return (0.0, 0.0)

    for r in model['reactions']:
        rid = r['id']
        rx_el = ET.SubElement(lor, _sbml2('reaction'), {
            'id': _safe_id('r_', rid),
            'name': r.get('name') or rid,
            'reversible': 'true' if r['lower_bound'] < 0 < r['upper_bound'] else 'false',
        })
        reactants = [(mid, -c) for mid, c in r['metabolites'].items() if c < 0]
        products = [(mid, c) for mid, c in r['metabolites'].items() if c > 0]

        # Choose aliases: prefer those connected in the Escher drawing
        escher = drawn.get(rid)
        node_aliases = {}
        mid_xy = None
        if escher:
            _rk, erx = escher
            mid_nid = _midmarker(erx, nodes)
            if mid_nid:
                mid_xy = _xy(nodes[mid_nid])
            for path in _reaction_paths(erx, nodes):
                end = nodes.get(path[-1])
                if end and end.get('node_type') == 'metabolite' and path[-1] in alias_for_node:
                    node_aliases[end['bigg_id']] = alias_for_node[path[-1]]

        def pick_alias(mid):
            if mid in node_aliases:
                return node_aliases[mid]
            return (alias_by_mid.get(mid) or [None])[0]

        base_r = pick_alias(reactants[0][0]) if reactants else None
        base_p = pick_alias(products[0][0]) if products else None

        r_ann = ET.SubElement(rx_el, _sbml2('annotation'))
        r_ext = ET.SubElement(r_ann, _cd('extension'))
        ET.SubElement(r_ext, _cd('reactionType')).text = 'STATE_TRANSITION'
        brs = ET.SubElement(r_ext, _cd('baseReactants'))
        bps = ET.SubElement(r_ext, _cd('baseProducts'))
        if base_r and reactants:
            ET.SubElement(brs, _cd('baseReactant'), {'species': species_el[reactants[0][0]], 'alias': base_r})
        if base_p and products:
            ET.SubElement(bps, _cd('baseProduct'), {'species': species_el[products[0][0]], 'alias': base_p})

        # Side links for remaining participants
        if len(reactants) > 1:
            lrl = ET.SubElement(r_ext, _cd('listOfReactantLinks'))
            for mid, _sto in reactants[1:]:
                aid = pick_alias(mid)
                if not aid:
                    continue
                link = ET.SubElement(lrl, _cd('reactantLink'), {
                    'reactant': species_el[mid], 'alias': aid, 'targetLineIndex': '-1,0',
                })
                ET.SubElement(link, _cd('line'), {'width': '1.0', 'color': 'ff000000', 'type': 'Straight'})
        if len(products) > 1:
            lpl = ET.SubElement(r_ext, _cd('listOfProductLinks'))
            for mid, _sto in products[1:]:
                aid = pick_alias(mid)
                if not aid:
                    continue
                link = ET.SubElement(lpl, _cd('productLink'), {
                    'product': species_el[mid], 'alias': aid, 'targetLineIndex': '-1,1',
                })
                ET.SubElement(link, _cd('line'), {'width': '1.0', 'color': 'ff000000', 'type': 'Straight'})

        # editPoints in local frame between base reactant and base product
        edit_text = ''
        rect_index = '0'
        if base_r and base_p:
            start = alias_xy(base_r)
            end = alias_xy(base_p)
            interiors = []
            if mid_xy is not None:
                interiors = [mid_xy]
            locals_ = _global_to_local(start, end, interiors)
            if locals_:
                edit_text = ' '.join(f'{x:.8f},{y:.8f}' for x, y in locals_)
                rect_index = '0'
        scheme = ET.SubElement(r_ext, _cd('connectScheme'), {
            'connectPolicy': 'direct',
            'rectangleIndex': rect_index,
        })
        lod = ET.SubElement(scheme, _cd('listOfLineDirection'))
        ET.SubElement(lod, _cd('lineDirection'), {'index': '0', 'value': 'unknown'})
        if edit_text:
            ET.SubElement(lod, _cd('lineDirection'), {'index': '1', 'value': 'unknown'})
            ET.SubElement(r_ext, _cd('editPoints')).text = edit_text
        ET.SubElement(r_ext, _cd('line'), {'width': '1.0', 'color': 'ff000000'})

        # SBML reactants / products with alias annotations
        lorct = ET.SubElement(rx_el, _sbml2('listOfReactants'))
        for mid, sto in reactants:
            sr = ET.SubElement(lorct, _sbml2('speciesReference'), {
                'species': species_el[mid],
                'stoichiometry': str(sto),
            })
            aid = pick_alias(mid)
            if aid:
                s_ann = ET.SubElement(sr, _sbml2('annotation'))
                s_ext = ET.SubElement(s_ann, _cd('extension'))
                ET.SubElement(s_ext, _cd('alias')).text = aid
        lop = ET.SubElement(rx_el, _sbml2('listOfProducts'))
        for mid, sto in products:
            sr = ET.SubElement(lop, _sbml2('speciesReference'), {
                'species': species_el[mid],
                'stoichiometry': str(sto),
            })
            aid = pick_alias(mid)
            if aid:
                s_ann = ET.SubElement(sr, _sbml2('annotation'))
                s_ext = ET.SubElement(s_ann, _cd('extension'))
                ET.SubElement(s_ext, _cd('alias')).text = aid

    return root


def write_celldesigner(path, data, model):
    root = build_celldesigner(data, model)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Pretty-print for readability
    rough = ET.tostring(root, encoding='utf-8')
    try:
        pretty = minidom.parseString(rough).toprettyxml(indent='  ', encoding='utf-8')
        path.write_bytes(pretty)
    except Exception:
        ET.ElementTree(root).write(path, encoding='utf-8', xml_declaration=True)
    return path


# ------------------------------------------------------------------------------------------------ Orchestration


EXPORT_FILES = (
    ('map.sbgn', 'sbgn', 'Completed map, SBGN-ML'),
    ('map.sbml', 'sbml', 'Completed map, SBML with Layout'),
    ('map.celldesigner.xml', 'xml', 'Completed map, CellDesigner SBML'),
)


def write_all_exports(folder, model_path=None, map_path=None, model=None):
    """Write SBGN-ML, SBML Layout and CellDesigner SBML next to map.json.

    Returns a list of written filenames. Missing map.json yields an empty list.
    """
    folder = Path(folder)
    map_path = Path(map_path) if map_path else folder / 'map.json'
    if not map_path.is_file():
        return []
    data = read(map_path)
    if model is None:
        if model_path is None:
            upload = next(iter(sorted((folder / 'upload').glob('model.*'))), None)
            if upload is None:
                raise FileNotFoundError(f'no model for export in {folder}')
            model_path = upload
        model = load_model(model_path)
    model_path = Path(model_path) if model_path else None

    written = []
    write_sbgn(folder / 'map.sbgn', data, model)
    written.append('map.sbgn')

    if model_path is not None and model_path.suffix.lower() in ('.xml', '.sbml'):
        write_sbml_layout(folder / 'map.sbml', model_path, data, model=model)
        written.append('map.sbml')
    else:
        # JSON GEM: build a temporary SBML from CellDesigner export is not Layout;
        # write CellDesigner only plus a layout-less note via CellDesigner path.
        pass

    write_celldesigner(folder / 'map.celldesigner.xml', data, model)
    written.append('map.celldesigner.xml')
    return written


def collect_export_ids(data, model):
    """Return sets of reaction/species ids present in the Escher map body."""
    body = _body(data)
    rx = {r.get('bigg_id') for r in body['reactions'].values() if r.get('bigg_id')}
    sp = {n.get('bigg_id') for n in body['nodes'].values()
          if n.get('node_type') == 'metabolite' and n.get('bigg_id')}
    return rx, sp


def sbgn_ids(path):
    """Reaction ids (from process glyph ids) and species labels from SBGN-ML."""
    root = ET.parse(path).getroot()

    def local(tag):
        return tag.rsplit('}', 1)[-1]

    processes, entities = set(), set()
    for g in root.iter():
        if local(g.tag) != 'glyph':
            continue
        cls = (g.get('class') or '').strip()
        gid = g.get('id') or ''
        label = next((c.get('text') for c in g if local(c.tag) == 'label' and c.get('text')), None)
        if cls == 'process':
            rid = gid[3:] if gid.startswith('re_') else (label or gid)
            if rid:
                processes.add(rid)
        elif cls in ('simple chemical', 'macromolecule', 'complex', 'unspecified entity') and label:
            entities.add(label)
    return processes, entities


def layout_ids(path):
    """Reaction and species ids referenced by the SBML Layout package."""
    import libsbml
    doc = libsbml.readSBMLFromFile(str(path))
    model = doc.getModel()
    plugin = model.getPlugin('layout')
    layout = plugin.getLayout(0) if plugin and plugin.getNumLayouts() else None
    rx, sp = set(), set()
    if layout is None:
        return rx, sp
    for i in range(layout.getNumReactionGlyphs()):
        rg = layout.getReactionGlyph(i)
        rid = rg.getReactionId()
        if rid.startswith('R_'):
            rid = rid[2:]
        rx.add(rid)
    for i in range(layout.getNumSpeciesGlyphs()):
        sg = layout.getSpeciesGlyph(i)
        sid = sg.getSpeciesId()
        if sid.startswith('M_'):
            sid = sid[2:]
        sp.add(sid)
    return rx, sp


def celldesigner_ids(path):
    """Reaction and species ids from a CellDesigner SBML export."""
    root = ET.parse(path).getroot()

    def local(tag):
        return tag.rsplit('}', 1)[-1]

    rx, sp = set(), set()
    for e in root.iter():
        if local(e.tag) == 'reaction' and 'sbml' in e.tag:
            name = e.get('name') or ''
            rid = e.get('id') or ''
            # Our writer uses r_<id> and name=original id/name
            if name:
                rx.add(name)
            stripped = rid[2:] if rid.startswith('r_') else rid
            if stripped:
                rx.add(stripped)
        if local(e.tag) == 'species' and 'sbml' in e.tag:
            sid = e.get('id') or ''
            name = e.get('name') or ''
            stripped = sid[2:] if sid.startswith('s_') else sid
            if stripped:
                sp.add(stripped)
            if name:
                sp.add(name)
    return rx, sp


def coordinates_match(data, sbgn_path, tol=1.0):
    """True when every Escher metabolite/midmarker centre appears as a glyph centre within tol."""
    from evaluation.reference_sbgn import read_sbgn_positions
    # reference_sbgn flips y; compare absolute with our writer (no flip)
    root = ET.parse(sbgn_path).getroot()

    def local(tag):
        return tag.rsplit('}', 1)[-1]

    centres = []
    for g in root.iter():
        if local(g.tag) != 'glyph':
            continue
        cls = (g.get('class') or '').strip()
        if cls in ('clone marker', 'unit of information', 'state variable'):
            continue
        bbox = next((c for c in g if local(c.tag) == 'bbox'), None)
        if bbox is None:
            continue
        x = float(bbox.get('x')) + float(bbox.get('w')) / 2.0
        y = float(bbox.get('y')) + float(bbox.get('h')) / 2.0
        centres.append((x, y))

    body = _body(data)
    needed = []
    for n in body['nodes'].values():
        if n.get('node_type') in ('metabolite', 'midmarker'):
            needed.append(_xy(n))
    for x, y in needed:
        if not any(math.hypot(x - cx, y - cy) <= tol for cx, cy in centres):
            return False
    return True
