"""Line style rendering of an Escher map: a reaction is a LINE between its metabolites, not a node.

Expert maps (Escher, WikiPathways, Reactome) show metabolites as the nodes and reactions as curves/lines between them, with the reaction
name written beside the line.  The Escher mid-marker is only an anchor for the label and the branch point; it is not drawn here.
Exported reaction lines and arrowheads are a single dark stroke (#222); main vs currency width and opacity are unchanged.
The tool viewer may recolour those paths from the `data-kind` provenance attribute.

Input is the Escher JSON written by the drawing pipeline (map.json); nothing about the chemistry is changed, only how it is displayed.
The old marker style stays available (`reaction_style = "marker"`): the drawing code in drawing.py is not modified.
"""
import json
from html import escape

import ec_provenance as epv

FONT = 'Menlo,monospace'
PRIMARY_R, SECONDARY_R = 5.0, 3.2


def _terminals(reaction, nodes):
    """Polylines from the mid-marker to every metabolite of the reaction, following the stored route segments."""
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


def _segment(reaction, a, b):
    """Segment between nodes a and b in either direction; control points returned in the a -> b direction (or None)."""
    for seg in reaction['segments'].values():
        if seg['from_node_id'] == a and seg['to_node_id'] == b:
            return (seg.get('b1'), seg.get('b2')) if seg.get('b1') and seg.get('b2') else None
        if seg['from_node_id'] == b and seg['to_node_id'] == a:
            return (seg.get('b2'), seg.get('b1')) if seg.get('b1') and seg.get('b2') else None
    return None


def _curved_path(reaction, path, nodes):
    """SVG path data along `path` using the Escher segments' Bezier control points; returns (d, curved)."""
    d, curved = f"M{nodes[path[0]]['x']:.1f},{nodes[path[0]]['y']:.1f}", False
    for a, b in zip(path, path[1:]):
        ctl = _segment(reaction, a, b)
        if ctl:
            curved = True
            d += f" C{ctl[0]['x']:.1f},{ctl[0]['y']:.1f} {ctl[1]['x']:.1f},{ctl[1]['y']:.1f} {nodes[b]['x']:.1f},{nodes[b]['y']:.1f}"
        else:
            d += f" L{nodes[b]['x']:.1f},{nodes[b]['y']:.1f}"
    return d, curved


EDGE = '#222'  # exported stroke / arrowhead; the tool viewer may override via data-kind


def render(data, kinds=None):
    """SVG text of the map in line style.  `kinds` maps reaction id -> provenance kind.

    Lines and arrowheads are drawn in EDGE.  Each path carries data-kind when known so
    the Map Generator viewer can recolour by provenance; legend counts stay in the UI.
    """
    body = data[1]
    nodes, reactions = body['nodes'], body['reactions']
    canvas = body.get('canvas') or {}
    width, height = canvas.get('width', 1000), canvas.get('height', 1000)
    kinds = kinds or {}
    defs = ''.join(f'<marker id="ar-{k}" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="5" markerHeight="5" orient="auto-start-reverse"><path d="M0,0L8,4L0,8Z" fill="{EDGE}"/></marker>' for k in list(epv.COLORS) + ['none'])
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width:.0f} {height:.0f}" width="{width:.0f}" height="{height:.0f}"><rect width="100%" height="100%" fill="white"/><defs>{defs}</defs>']
    preserved = {(a, b) for a, b in (data[0].get('metadata', {}).get('preserved_edges') or [])}
    source_paths = {(p['reaction'], p['metabolite']): p for p in (data[0].get('metadata', {}).get('source_paths') or [])}
    routed_edges = {(p['reaction'], p['metabolite']) for p in (data[0].get('metadata', {}).get('routed_edges') or [])}
    curves_drawn = 0
    lines, top_lines, labels = [], [], []
    for rk, r in reactions.items():
        rid = r['bigg_id']
        kind = kinds.get(rid)
        mark = f'url(#ar-{kind if kind in epv.COLORS else "none"})'
        kind_attr = f' data-kind="{escape(kind)}"' if kind else ''
        reversible = bool(r.get('reversibility'))
        coefficient = {}
        for m in r.get('metabolites', []):
            coefficient.setdefault(m['bigg_id'], m['coefficient'])
        for path in _terminals(r, nodes):
            end = nodes[path[-1]]
            if end.get('node_type') != 'metabolite':
                continue
            keep = (rid, end.get('bigg_id')) in preserved
            product = coefficient.get(end.get('bigg_id'), 1) > 0
            d_curve, curved = _curved_path(r, path, nodes)
            entry = source_paths.get((rid, end.get('bigg_id')))
            if (curved or entry is not None) and (rid, end.get('bigg_id')) not in routed_edges:  # the expert's own line: Bezier segments or the polyline of an original connector, stored in the map as Escher segments
                curves_drawn += 1
                status = (entry or {}).get('status', 'preserved')
                title = f'<title>{escape(rid)} {escape(r.get("name") or "")} - {escape(epv.LABELS[kind]) if kind in epv.LABELS else "unclassified"} (expert {"curve" if curved else "line"}{", adjusted" if status == "adjusted" else ""})</title>'
                top_lines.append(f'<path data-rx="{escape(rid)}" data-met="{escape(end.get("bigg_id") or "")}"{kind_attr} d="{d_curve}" stroke="{EDGE}" stroke-width="2.4" stroke-linejoin="round" stroke-linecap="round" fill="none"{(" marker-end=" + chr(34) + mark + chr(34)) if (product or reversible) else ""}>{title}</path>')
                continue
            if (rid, end.get('bigg_id')) in routed_edges:  # an algorithm-routed detour with rounded corners: the same segments as in map.json
                arrows = f' marker-end="{mark}"' if (product or reversible) else ''
                title = f'<title>{escape(rid)} {escape(r.get("name") or "")} - {escape(epv.LABELS[kind]) if kind in epv.LABELS else "unclassified"} (routed around obstacles)</title>'
                lines.append(f'<path data-rx="{escape(rid)}" data-met="{escape(end.get("bigg_id") or "")}"{kind_attr} d="{d_curve}" stroke="{EDGE}" stroke-width="1.3" stroke-opacity=".7" stroke-linejoin="round" stroke-linecap="round" fill="none"{arrows}>{title}</path>')
                continue
            if keep:
                path = [path[0], path[-1]]  # edge of an expert structure: keep the structure's own straight connection, not the routed detour
            pts = ' L'.join(f"{nodes[k]['x']:.1f},{nodes[k]['y']:.1f}" for k in path)
            product = coefficient.get(end.get('bigg_id'), 1) > 0
            arrows = f' marker-end="{mark}"' if (product or reversible) else ''
            title = f'<title>{escape(rid)} {escape(r.get("name") or "")} - {escape(epv.LABELS[kind]) if kind in epv.LABELS else "unclassified"}</title>'
            (top_lines if keep else lines).append(f'<path data-rx="{escape(rid)}" data-met="{escape(end.get("bigg_id") or "")}"{kind_attr} d="M{pts}" stroke="{EDGE}" stroke-width="{2.4 if keep else 1.3}"{"" if keep else " stroke-opacity=\".7\""} stroke-linejoin="round" stroke-linecap="round" fill="none"{arrows}>{title}</path>')
        labels.append(f'<text data-rk="{escape(rk)}" data-rl="{escape(rid)}" x="{r["label_x"]:.1f}" y="{r["label_y"]:.1f}" font-family="{FONT}" font-size="12" fill="#17323d" font-weight="600" stroke="#fff" stroke-width="3" paint-order="stroke">{escape(rid)}</text>')
    svg += lines + top_lines  # routed connections first, the expert structures' own edges on top
    hidden = set((data[0].get('metadata') or {}).get('hidden_metabolite_labels') or [])
    for nid, n in nodes.items():
        if n.get('node_type') != 'metabolite':
            continue
        primary = n.get('node_is_primary', True)
        radius = PRIMARY_R if primary else SECONDARY_R
        fill, stroke = ('#dcebdf', '#5d7480') if primary else ('#ffffff', '#8a9ca5')
        svg.append(f'<circle data-n="{escape(nid)}" data-m="{escape(n.get("bigg_id", ""))}" cx="{n["x"]:.1f}" cy="{n["y"]:.1f}" r="{radius}" fill="{fill}" stroke="{stroke}" stroke-width="1.2"/>')
        if n.get('label_hidden') or nid in hidden:
            continue  # dense-region currency hide
        labels.append(f'<text data-nl="{escape(nid)}" x="{n["label_x"]:.1f}" y="{n["label_y"]:.1f}" font-family="{FONT}" font-size="{12 if primary else 10.5}" fill="{"#17323d" if primary else "#6b7f88"}" stroke="#fff" stroke-width="3" paint-order="stroke">{escape(n.get("bigg_id", ""))}</text>')
    svg += labels
    svg.append('</svg>')
    return ''.join(svg)


def render_file(json_path, svg_path, kinds=None):
    data = json.loads(open(json_path).read())
    text = render(data, kinds)
    meta = data[0].get('metadata', {}).get('layout_sources')
    if meta:
        offset = text.index('>') + 1
        text = text[:offset] + '<metadata>' + escape(json.dumps(dict(layout_sources=meta), ensure_ascii=False)) + '</metadata>' + text[offset:]
    open(svg_path, 'w').write(text)
