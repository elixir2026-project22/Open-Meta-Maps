#!/usr/bin/env python3
"""Reproducible Atlas figures and duplication stats for docs.

  PYTHONPATH=layout_atlas python3 tools/atlas_module_figures.py all

Writes under docs/images/atlas/:
  atlas_overview.png          — dense mosaic of all default-load modules
  atlas_examples.svg/.png     — n×m captioned example grid
  atlas_overlap.json          — per-source counts, unique-after-dedup, nesting, pairwise overlap
  worked_*.svg/.png           — worked-example panels for atlas_data_management.md

Figure rules: canvas sized to content (no crop), uniform cells, grid-aligned,
labels kept clear of strokes. Mosaic PNG kept under ~3 MB; SVGs only for small panels.
"""
from __future__ import annotations

import argparse
import html
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'layout_atlas'))

from core import DATA  # noqa: E402

OUT = ROOT / 'docs' / 'images' / 'atlas'
CANON_PREF = (
    'bigg.reaction', 'metanetx.reaction', 'rhea', 'reactome', 'ec-code',
    'human.gem', 'sbgn.process', 'celldesigner.reaction',
)

# Display order and subtle region tints (light fills behind thumbnails).
SOURCE_ORDER = [
    'escher', 'reactome', 'wikipathways', 'wikipathways_geometry',
    'human_maps', 'brenda', 'metabolismregulation', 'panther', 'ich360',
]
SOURCE_LABEL = {
    'escher': 'Escher',
    'reactome': 'Reactome',
    'wikipathways': 'WikiPathways (semantic)',
    'wikipathways_geometry': 'WikiPathways (geometry)',
    'human_maps': 'Human-maps',
    'brenda': 'BRENDA',
    'metabolismregulation': 'Metabolism Regulation',
    'panther': 'PANTHER',
    'ich360': 'iCH360',
}
SOURCE_TINT = {
    'escher': '#dcefe6',
    'reactome': '#e4eaf6',
    'wikipathways': '#e8f0dc',
    'wikipathways_geometry': '#f0f4e6',
    'human_maps': '#f5e8dc',
    'brenda': '#eee4f2',
    'metabolismregulation': '#e6f2f2',
    'panther': '#f2e6ea',
    'ich360': '#e8ecef',
}
SOURCE_STROKE = {
    'escher': '#3f8f6b',
    'reactome': '#4a6fa5',
    'wikipathways': '#6a8f3f',
    'wikipathways_geometry': '#7a9a55',
    'human_maps': '#b07a45',
    'brenda': '#8a6a9a',
    'metabolismregulation': '#3f8a8a',
    'panther': '#9a5a6a',
    'ich360': '#5a6670',
}

# Representative example picks: (source_key, prefer_size, name_substr or None)
EXAMPLE_SPEC = [
    ('escher', 8, 'glycolysis'),
    ('escher', 12, 'TCA'),
    ('escher', 4, None),
    ('reactome', 8, None),
    ('reactome', 4, None),
    ('wikipathways', 8, None),
    ('wikipathways', 4, None),
    ('wikipathways_geometry', 8, None),
    ('human_maps', 8, None),
    ('human_maps', 12, None),
    ('ich360', 8, None),
    ('ich360', 4, None),
    ('brenda', 8, None),
    ('brenda', 4, None),
    ('metabolismregulation', 8, None),
    ('panther', 4, None),
]


def source_key(m):
    src = m.get('source') or '?'
    if src == 'wikipathways' and m.get('quality') == 'geometry_only':
        return 'wikipathways_geometry'
    return src


def primary_reaction_id(r, src):
    refs = r.get('references') or {}
    for ns in CANON_PREF:
        vals = [v for v in (refs.get(ns) or []) if v]
        if vals:
            return f'{ns}:{sorted(vals)[0]}'
    rid = r.get('id')
    return f'native:{src}:{rid}' if rid is not None else f'native:{src}:?'


def reaction_set(m):
    src = m.get('source') or '?'
    rx = m.get('reactions') or []
    if rx:
        return frozenset(primary_reaction_id(r, src) for r in rx)
    # Geometry-only WikiPathways windows often store no chemical reaction ids;
    # fall back to the geometry fingerprint so dedup stays meaningful.
    fp = m.get('geometry_fingerprint')
    if fp:
        return frozenset([f'geometry:{fp}'])
    ents = sorted({n.get('entity') for n in (m.get('nodes') or []) if n.get('entity')})
    return frozenset([f'entities:{",".join(ents)}']) if ents else frozenset()


def load_modules():
    """Load default Atlas modules (open-core sqlite + extension gzips)."""
    atlas = DATA / 'atlas.sqlite'
    open_core = DATA / 'open_core.sqlite'
    if open_core.exists() and (not atlas.exists() or atlas.stat().st_size < open_core.stat().st_size // 2):
        import shutil
        shutil.copyfile(open_core, atlas)
    from database import modules as load_all
    return load_all()


def compute_overlap(mods):
    by_src = defaultdict(list)
    sets = {}
    for m in mods:
        key = source_key(m)
        rs = reaction_set(m)
        sets[m['id']] = rs
        by_src[key].append(m)

    per_source = {}
    global_set_to_ids = defaultdict(list)
    for key, items in by_src.items():
        local = Counter()
        for m in items:
            rs = sets[m['id']]
            local[rs] += 1
            global_set_to_ids[rs].append((key, m['id']))
        nested = 0
        sorted_sets = sorted({sets[m['id']] for m in items}, key=len)
        # Index larger sets by size buckets for nesting checks
        by_size = defaultdict(list)
        for s in sorted_sets:
            by_size[len(s)].append(s)
        sizes = sorted(by_size)
        for s in sorted_sets:
            if len(s) == 0:
                continue
            contained = False
            for L in sizes:
                if L <= len(s):
                    continue
                for t in by_size[L]:
                    if s < t:
                        contained = True
                        break
                if contained:
                    break
            if contained:
                # count modules with this set that are nested (each module once)
                nested += local[s]
        per_source[key] = dict(
            modules=len(items),
            unique_reaction_sets=len(local),
            exact_duplicate_extra=sum(c - 1 for c in local.values() if c > 1),
            nested_in_larger=nested,
        )

    # Global unique across all sources
    all_sets = set(sets.values())
    # Cross-source: identical reaction sets appearing in ≥2 sources
    cross_shared_sets = 0
    pairwise = {a: {b: 0 for b in SOURCE_ORDER} for a in SOURCE_ORDER}
    set_sources = defaultdict(set)
    for rs, owners in global_set_to_ids.items():
        srcs = {s for s, _ in owners}
        set_sources[rs] = srcs
        if len(srcs) >= 2:
            cross_shared_sets += 1
            src_list = sorted(srcs)
            for i, a in enumerate(src_list):
                for b in src_list[i + 1:]:
                    if a in pairwise and b in pairwise[a]:
                        pairwise[a][b] += 1
                        pairwise[b][a] += 1

    # Global nesting (any larger module anywhere)
    nested_global = 0
    uniq_sets_list = sorted(all_sets, key=len)
    by_size = defaultdict(list)
    for s in uniq_sets_list:
        by_size[len(s)].append(s)
    sizes = sorted(by_size)
    nested_set_flags = {}
    for s in uniq_sets_list:
        if len(s) == 0:
            nested_set_flags[s] = False
            continue
        contained = False
        for L in sizes:
            if L <= len(s):
                continue
            for t in by_size[L]:
                if s < t:
                    contained = True
                    break
            if contained:
                break
        nested_set_flags[s] = contained
    for m in mods:
        if nested_set_flags.get(sets[m['id']]):
            nested_global += 1

    summary = dict(
        total_modules=len(mods),
        unique_reaction_sets_global=len(all_sets),
        modules_nested_in_larger_global=nested_global,
        cross_source_identical_sets=cross_shared_sets,
        by_source={k: per_source.get(k, dict(modules=0, unique_reaction_sets=0, exact_duplicate_extra=0, nested_in_larger=0)) for k in SOURCE_ORDER},
        pairwise_identical_sets=pairwise,
    )
    return summary, by_src, sets


def _bbox(nodes):
    xs = [n['xy'][0] for n in nodes]
    ys = [n['xy'][1] for n in nodes]
    return min(xs), min(ys), max(xs), max(ys)


def draw_module_paths(m, x0, y0, cell, stroke='#374151', lw=0.9, node_r=1.4, reaction_fill='#d77d38', met_fill='#197d87'):
    """Return SVG fragments for a module fitted into a square cell (padding inside)."""
    nodes = m.get('nodes') or []
    edges = m.get('edges') or []
    if len(nodes) < 2:
        return []
    xmin, ymin, xmax, ymax = _bbox(nodes)
    w = max(xmax - xmin, 1e-6)
    h = max(ymax - ymin, 1e-6)
    pad = cell * 0.12
    scale = (cell - 2 * pad) / max(w, h)
    cx = x0 + cell / 2
    cy = y0 + cell / 2
    mx = (xmin + xmax) / 2
    my = (ymin + ymax) / 2
    pos = {}
    for n in nodes:
        pos[n['key']] = (cx + (n['xy'][0] - mx) * scale, cy + (n['xy'][1] - my) * scale)
    parts = []
    for e in edges:
        a, b = pos.get(e['a']), pos.get(e['b'])
        if not a or not b:
            continue
        parts.append(f'<line x1="{a[0]:.2f}" y1="{a[1]:.2f}" x2="{b[0]:.2f}" y2="{b[1]:.2f}" stroke="{stroke}" stroke-width="{lw}" stroke-linecap="round"/>')
    for n in nodes:
        x, y = pos[n['key']]
        fill = reaction_fill if n.get('role') == 'reaction' or str(n.get('key', '')).startswith('r:') else met_fill
        parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{node_r}" fill="{fill}"/>')
    return parts


def pick_examples(by_src):
    picked = []
    used = set()
    for src, size, substr in EXAMPLE_SPEC:
        items = by_src.get(src) or []
        cand = [m for m in items if m['id'] not in used and 2 <= (m.get('size') or 0) <= 24]
        if size:
            prefer = [m for m in cand if m.get('size') == size]
            if prefer:
                cand = prefer + [m for m in cand if m not in prefer]
        if substr:
            named = [m for m in cand if substr.lower() in (m.get('name') or '').lower()]
            if named:
                cand = named
        # Prefer denser / mid-size drawings
        cand = sorted(cand, key=lambda m: (-min(m.get('size') or 0, 12), m.get('name') or '', m['id']))
        if not cand:
            continue
        m = cand[0]
        used.add(m['id'])
        picked.append(m)
    return picked


def write_examples(picked, overlap_summary):
    cols = 4
    rows = int(math.ceil(len(picked) / cols))
    cell = 150
    cap_h = 44
    gap = 16
    left = 20
    top = 20
    W = left * 2 + cols * cell + (cols - 1) * gap
    H = top * 2 + rows * (cell + cap_h) + (rows - 1) * gap
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="-apple-system,Segoe UI,Helvetica,Arial,sans-serif">',
        f'<title>Example Atlas modules</title>',
        f'<rect width="{W}" height="{H}" fill="#ffffff"/>',
    ]
    for i, m in enumerate(picked):
        r, c = divmod(i, cols)
        x = left + c * (cell + gap)
        y = top + r * (cell + cap_h + gap)
        sk = source_key(m)
        parts.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" rx="8" fill="#fafafa" stroke="{SOURCE_STROKE.get(sk, "#888")}" stroke-width="1.2"/>')
        parts.extend(draw_module_paths(m, x, y, cell, stroke='#4b5563', lw=1.1, node_r=2.0))
        name = (m.get('name') or m.get('map_id') or '')[:42]
        cap1 = SOURCE_LABEL.get(sk, sk)
        cap2 = html.escape(name)
        cap3 = f"{m.get('size', len(m.get('reactions') or []))} reactions"
        tx = x + cell / 2
        parts.append(f'<text x="{tx:.1f}" y="{y + cell + 14}" text-anchor="middle" font-size="11" font-weight="600" fill="#1f2937">{html.escape(cap1)}</text>')
        parts.append(f'<text x="{tx:.1f}" y="{y + cell + 28}" text-anchor="middle" font-size="10" fill="#4b5563">{cap2}</text>')
        parts.append(f'<text x="{tx:.1f}" y="{y + cell + 40}" text-anchor="middle" font-size="10" fill="#6b7280">{cap3}</text>')
    parts.append('</svg>')
    svg = '\n'.join(parts)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'atlas_examples.svg').write_text(svg, encoding='utf-8')
    _svg_to_png(OUT / 'atlas_examples.svg', OUT / 'atlas_examples.png', width=900)
    return OUT / 'atlas_examples.svg', OUT / 'atlas_examples.png'


def _hex_rgb(h):
    h = h.lstrip('#')
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))


def write_overview_mosaic(by_src, summary):
    """Dense PNG mosaic via PIL (fast): contiguous source regions, tint, legend."""
    from PIL import Image, ImageDraw, ImageFont

    cell = 10
    gap = 1
    region_gap = 8
    pad = 14
    legend_line = 14
    legend_h = 36 + legend_line * (1 + sum(1 for sk in SOURCE_ORDER if (by_src.get(sk))))

    blocks = []
    for sk in SOURCE_ORDER:
        items = by_src.get(sk) or []
        if not items:
            continue
        cols = max(10, int(math.ceil(math.sqrt(len(items)) * 1.4)))
        rows = int(math.ceil(len(items) / cols))
        items = sorted(items, key=lambda m: (m.get('size') or 0, m['id']))
        blocks.append((sk, items, cols, rows))

    max_content_w = 2200
    layout = []
    x = pad
    y = pad
    row_h = 0
    canvas_w = pad
    for sk, items, cols, rows in blocks:
        bw = cols * (cell + gap)
        bh = rows * (cell + gap)
        if x > pad and x + bw > max_content_w + pad:
            x = pad
            y += row_h + region_gap
            row_h = 0
        layout.append((sk, items, cols, rows, x, y))
        row_h = max(row_h, bh)
        x += bw + region_gap
        canvas_w = max(canvas_w, x)
    content_bottom = max((y0 + rows * (cell + gap) for _, _, _, rows, _, y0 in layout), default=pad)
    canvas_h = content_bottom + pad + legend_h + pad
    canvas_w = max(canvas_w + pad, pad + 900)

    im = Image.new('RGB', (canvas_w, canvas_h), (255, 255, 255))
    draw = ImageDraw.Draw(im)
    try:
        font = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf', 12)
        font_sm = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf', 11)
        font_bold = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf', 13)
    except Exception:
        font = font_sm = font_bold = ImageFont.load_default()

    ink = (55, 65, 81)
    met = (25, 125, 135)
    rxn = (215, 125, 56)

    for sk, items, cols, rows, x0, y0 in layout:
        bw = cols * (cell + gap)
        bh = rows * (cell + gap)
        tint = _hex_rgb(SOURCE_TINT.get(sk, '#eeeeee'))
        stroke = _hex_rgb(SOURCE_STROKE.get(sk, '#888888'))
        draw.rectangle([x0 - 2, y0 - 2, x0 + bw + 1, y0 + bh + 1], fill=tint, outline=stroke, width=1)
        for i, m in enumerate(items):
            r, c = divmod(i, cols)
            cx = x0 + c * (cell + gap)
            cy = y0 + r * (cell + gap)
            _pil_draw_module(draw, m, cx, cy, cell, ink=ink, met=met, rxn=rxn)

    ly = content_bottom + pad
    draw.text((pad, ly), 'Default Atlas modules by source (cut windows; several can come from one map)', fill=(31, 41, 55), font=font_bold)
    ly += 18
    draw.text(
        (pad, ly),
        f"Total cut modules: {summary['total_modules']:,}   ·   Unique after dedup: {summary['unique_reaction_sets_global']:,}",
        fill=(75, 85, 99),
        font=font_sm,
    )
    ly += 16
    for sk in SOURCE_ORDER:
        info = summary['by_source'].get(sk) or {}
        if not info.get('modules'):
            continue
        tint = _hex_rgb(SOURCE_TINT.get(sk, '#eeeeee'))
        stroke = _hex_rgb(SOURCE_STROKE.get(sk, '#888888'))
        draw.rectangle([pad, ly, pad + 11, ly + 11], fill=tint, outline=stroke)
        label = (
            f"{SOURCE_LABEL.get(sk, sk)}: {info['modules']:,} modules"
            f" · {info['unique_reaction_sets']:,} unique after dedup"
        )
        draw.text((pad + 18, ly - 1), label, fill=(55, 65, 81), font=font_sm)
        ly += legend_line

    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / 'atlas_overview.png'
    im.save(path, format='PNG', optimize=True)
    _shrink_png(path, max_bytes=3_000_000)
    return path


def _pil_draw_module(draw, m, x0, y0, cell, ink=(55, 65, 81), met=(25, 125, 135), rxn=(215, 125, 56)):
    nodes = m.get('nodes') or []
    edges = m.get('edges') or []
    if len(nodes) < 2:
        return
    xmin, ymin, xmax, ymax = _bbox(nodes)
    w = max(xmax - xmin, 1e-6)
    h = max(ymax - ymin, 1e-6)
    pad = cell * 0.14
    scale = (cell - 2 * pad) / max(w, h)
    cx = x0 + cell / 2
    cy = y0 + cell / 2
    mx = (xmin + xmax) / 2
    my = (ymin + ymax) / 2
    pos = {}
    for n in nodes:
        pos[n['key']] = (
            cx + (n['xy'][0] - mx) * scale,
            cy + (n['xy'][1] - my) * scale,
        )
    for e in edges:
        a, b = pos.get(e['a']), pos.get(e['b'])
        if a and b:
            draw.line([a, b], fill=ink, width=1)
    for n in nodes:
        x, y = pos[n['key']]
        is_r = n.get('role') == 'reaction' or str(n.get('key', '')).startswith('r:')
        r = 1.2 if is_r else 1.0
        fill = rxn if is_r else met
        draw.ellipse([x - r, y - r, x + r, y + r], fill=fill)


def _shrink_png(path: Path, max_bytes: int):
    from PIL import Image
    if path.stat().st_size <= max_bytes:
        return
    im = Image.open(path).convert('RGB')
    # progressive downscale / quality
    for scale in (0.85, 0.75, 0.65, 0.55):
        w, h = im.size
        small = im.resize((max(1, int(w * scale)), max(1, int(h * scale))), Image.Resampling.LANCZOS)
        small.save(path, format='PNG', optimize=True)
        if path.stat().st_size <= max_bytes:
            return
    # last resort: palette
    Image.open(path).convert('P', palette=Image.ADAPTIVE, colors=64).save(path, optimize=True)


def _svg_to_png(svg_path: Path, png_path: Path, width: int = 800):
    try:
        import cairosvg
        cairosvg.svg2png(url=svg_path.resolve().as_uri(), write_to=str(png_path), output_width=width)
        return
    except Exception:
        pass
    # Fallback: rasterise via wand/PIL is unavailable for arbitrary SVG; leave SVG only.
    print(f'PNG preview skipped for {svg_path.name} (install cairosvg for PNG)', flush=True)


def write_worked_examples(by_src, mods, sets):
    """Small SVG panels for the worked-examples section."""
    OUT.mkdir(parents=True, exist_ok=True)
    paths = {}

    # (1) Cut example: pick an Escher map with several modules
    escher = by_src.get('escher') or []
    by_map = defaultdict(list)
    for m in escher:
        by_map[m.get('map_id')].append(m)
    # Prefer a map with 8–40 modules of varied sizes
    map_id, group = max(by_map.items(), key=lambda kv: (min(40, len(kv[1])), -abs(len(kv[1]) - 20)))
    group = sorted(group, key=lambda m: (m.get('size') or 0, m['id']))
    # Show whole-map hull + up to 6 module cells
    paths['cut'] = _worked_cut_figure(group[:1][0] if False else group, map_id)

    # (2) Canonicalisation: Human-maps reaction with bigg/mnx/rhea
    hm_cands = [
        m for m in (by_src.get('human_maps') or [])
        if sum(1 for r in (m.get('reactions') or []) if (r.get('references') or {}).get('bigg.reaction')) >= 3
    ]
    hm = max(hm_cands, key=lambda m: m.get('size') or 0) if hm_cands else None
    paths['ids'] = _worked_ids_figure(hm)

    # (3) Dedup across sources: find identical reaction set in two sources
    set_owners = defaultdict(list)
    for m in mods:
        rs = sets[m['id']]
        if len(rs) >= 3:
            set_owners[rs].append(m)
    pair = None
    for rs, owners in set_owners.items():
        srcs = {source_key(m) for m in owners}
        if len(srcs) >= 2:
            # pick one from each of two sources
            a = next(m for m in owners if source_key(m) == sorted(srcs)[0])
            b = next(m for m in owners if source_key(m) == sorted(srcs)[1])
            pair = (a, b, rs)
            break
    paths['dedup'] = _worked_dedup_figure(pair)

    # (4) Provenance record
    sample = next((m for m in escher if (m.get('size') or 0) >= 4 and m.get('extraction')), escher[0] if escher else None)
    paths['prov'] = _worked_provenance_figure(sample)

    # (5) Compact selection schematic (optional)
    paths['select'] = _worked_select_figure()

    return paths


def _worked_cut_figure(group, map_id):
    # Panel A: scatter all module bboxes on shared coords from first large module's map — simpler: show 6 modules in a row with nesting note
    show = group[:6] if len(group) >= 6 else group
    cell = 120
    gap = 14
    left = 20
    top = 36
    W = left * 2 + len(show) * cell + (len(show) - 1) * gap
    H = top + cell + 56
    name = html.escape((show[0].get('name') if show else map_id) or str(map_id))
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="-apple-system,Segoe UI,Helvetica,Arial,sans-serif">',
        f'<rect width="{W}" height="{H}" fill="#fff"/>',
        f'<text x="{left}" y="22" font-size="13" font-weight="700" fill="#1f2937">Cut windows from one Escher map</text>',
        f'<text x="{left}" y="{H - 14}" font-size="11" fill="#6b7280">{name} · {len(group)} modules stored from this map (showing {len(show)})</text>',
    ]
    for i, m in enumerate(show):
        x = left + i * (cell + gap)
        y = top
        nested = any(reaction_set(m) < reaction_set(o) for o in group if o['id'] != m['id'])
        stroke = '#b9852f' if nested else '#3f8f6b'
        parts.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" rx="8" fill="#fafafa" stroke="{stroke}" stroke-width="1.2"/>')
        parts.extend(draw_module_paths(m, x, y, cell, lw=1.0, node_r=1.8))
        parts.append(f'<text x="{x + cell/2:.1f}" y="{y + cell + 14}" text-anchor="middle" font-size="10" fill="#374151">{m.get("size", "?")} rx{" · nested" if nested else ""}</text>')
    parts.append('</svg>')
    path = OUT / 'worked_cut.svg'
    path.write_text('\n'.join(parts), encoding='utf-8')
    _svg_to_png(path, OUT / 'worked_cut.png', width=900)
    return path


def _worked_ids_figure(m):
    W, H = 720, 220
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="-apple-system,Segoe UI,Helvetica,Arial,sans-serif">',
        f'<rect width="{W}" height="{H}" fill="#fff"/>',
        f'<text x="20" y="28" font-size="13" font-weight="700" fill="#1f2937">Canonical reaction identifiers (Human-maps example)</text>',
    ]
    if not m:
        parts.append('<text x="20" y="60" font-size="12" fill="#6b7280">No example module found.</text></svg>')
        path = OUT / 'worked_ids.svg'
        path.write_text('\n'.join(parts), encoding='utf-8')
        return path
    rows = []
    for r in (m.get('reactions') or [])[:4]:
        refs = r.get('references') or {}
        rows.append((r.get('id'), refs))
    y = 56
    parts.append(f'<text x="20" y="48" font-size="11" fill="#6b7280">Map: {html.escape(m.get("name") or "")}</text>')
    # header
    headers = ['Native (human.gem)', 'bigg.reaction', 'metanetx.reaction', 'rhea']
    xs = [20, 200, 360, 540]
    for x, h in zip(xs, headers):
        parts.append(f'<text x="{x}" y="{y}" font-size="11" font-weight="600" fill="#1f2937">{h}</text>')
    y += 8
    parts.append(f'<line x1="20" y1="{y}" x2="{W-20}" y2="{y}" stroke="#d0d7de"/>')
    y += 18
    for rid, refs in rows:
        vals = [
            (refs.get('human.gem') or [rid])[0],
            (refs.get('bigg.reaction') or ['—'])[0],
            (refs.get('metanetx.reaction') or ['—'])[0],
            (refs.get('rhea') or ['—'])[0],
        ]
        for x, v in zip(xs, vals):
            parts.append(f'<text x="{x}" y="{y}" font-size="11" fill="#374151">{html.escape(str(v))}</text>')
        y += 22
    parts.append('</svg>')
    path = OUT / 'worked_ids.svg'
    path.write_text('\n'.join(parts), encoding='utf-8')
    _svg_to_png(path, OUT / 'worked_ids.png', width=720)
    return path


def _worked_dedup_figure(pair):
    W, H = 640, 260
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="-apple-system,Segoe UI,Helvetica,Arial,sans-serif">',
        f'<rect width="{W}" height="{H}" fill="#fff"/>',
        f'<text x="20" y="28" font-size="13" font-weight="700" fill="#1f2937">Same canonical reaction set in two sources</text>',
    ]
    if not pair:
        parts.append('<text x="20" y="60" font-size="12" fill="#6b7280">No cross-source identical set with ≥3 reactions in the default Atlas.</text></svg>')
        path = OUT / 'worked_dedup.svg'
        path.write_text('\n'.join(parts), encoding='utf-8')
        return path
    a, b, rs = pair
    cell = 150
    for i, m in enumerate((a, b)):
        x = 40 + i * (cell + 80)
        y = 50
        sk = source_key(m)
        parts.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" rx="8" fill="#fafafa" stroke="{SOURCE_STROKE.get(sk,"#888")}" stroke-width="1.2"/>')
        parts.extend(draw_module_paths(m, x, y, cell, lw=1.0, node_r=1.8))
        parts.append(f'<text x="{x+cell/2:.1f}" y="{y+cell+16}" text-anchor="middle" font-size="11" font-weight="600" fill="#1f2937">{html.escape(SOURCE_LABEL.get(sk, sk))}</text>')
        parts.append(f'<text x="{x+cell/2:.1f}" y="{y+cell+30}" text-anchor="middle" font-size="10" fill="#6b7280">{html.escape((m.get("name") or "")[:36])}</text>')
    parts.append('<text x="320" y="130" text-anchor="middle" font-size="18" fill="#6b7280">=</text>')
    parts.append(f'<text x="20" y="{H-16}" font-size="11" fill="#6b7280">{len(rs)} shared canonical reaction ids · both modules kept in the library; selection prefers the denser drawing at compose time</text>')
    parts.append('</svg>')
    path = OUT / 'worked_dedup.svg'
    path.write_text('\n'.join(parts), encoding='utf-8')
    _svg_to_png(path, OUT / 'worked_dedup.png', width=640)
    return path


def _worked_provenance_figure(m):
    W, H = 640, 240
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="ui-monospace,Menlo,monospace">',
        f'<rect width="{W}" height="{H}" fill="#fff"/>',
        f'<text x="20" y="28" font-size="13" font-weight="700" font-family="-apple-system,Segoe UI,Helvetica,Arial,sans-serif" fill="#1f2937">Provenance fields on one module</text>',
    ]
    if not m:
        parts.append('</svg>')
        path = OUT / 'worked_provenance.svg'
        path.write_text('\n'.join(parts), encoding='utf-8')
        return path
    ex = m.get('extraction') or {}
    lines = [
        f'source:   {m.get("source")}',
        f'map:      {(m.get("name") or m.get("map_id") or "")[:70]}',
        f'license:  {m.get("license")}',
        f'method:   {ex.get("method")}',
        f'version:  {ex.get("version")}',
        f'size:     {m.get("size")} reactions',
    ]
    y = 56
    for line in lines:
        parts.append(f'<text x="28" y="{y}" font-size="12" fill="#374151">{html.escape(line)}</text>')
        y += 22
    # tiny drawing
    parts.append(f'<rect x="{W-170}" y="48" width="140" height="140" rx="8" fill="#fafafa" stroke="#d0d7de"/>')
    parts.extend(draw_module_paths(m, W - 170, 48, 140, lw=1.0, node_r=1.8))
    parts.append('</svg>')
    path = OUT / 'worked_provenance.svg'
    path.write_text('\n'.join(parts), encoding='utf-8')
    _svg_to_png(path, OUT / 'worked_provenance.png', width=640)
    return path


def _worked_select_figure():
    W, H = 640, 200
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="-apple-system,Segoe UI,Helvetica,Arial,sans-serif">',
        f'<rect width="{W}" height="{H}" fill="#fff"/>',
        f'<text x="20" y="28" font-size="13" font-weight="700" fill="#1f2937">Compact selection among overlapping candidates</text>',
        # two boxes same reaction set
        f'<rect x="40" y="50" width="200" height="100" rx="8" fill="#f5f5f5" stroke="#9ca3af" stroke-width="1.2"/>',
        f'<text x="140" y="78" text-anchor="middle" font-size="12" fill="#6b7280">Sprawling whole-map window</text>',
        f'<text x="140" y="98" text-anchor="middle" font-size="11" fill="#9ca3af">same reaction set</text>',
        f'<text x="140" y="128" text-anchor="middle" font-size="11" fill="#9ca3af">lower density → dropped</text>',
        f'<text x="280" y="105" text-anchor="middle" font-size="20" fill="#6b7280">→</text>',
        f'<rect x="320" y="50" width="200" height="100" rx="8" fill="#e7f4ee" stroke="#3f8f6b" stroke-width="1.6"/>',
        f'<text x="420" y="78" text-anchor="middle" font-size="12" font-weight="600" fill="#1f2937">Compact pathway window</text>',
        f'<text x="420" y="98" text-anchor="middle" font-size="11" fill="#3f8f6b">same reaction set</text>',
        f'<text x="420" y="128" text-anchor="middle" font-size="11" fill="#3f8f6b">higher density → kept</text>',
        f'<text x="20" y="{H-18}" font-size="11" fill="#6b7280">Coverage still ranks first; compactness breaks ties for same-source biochemical matches with identical reaction sets.</text>',
        '</svg>',
    ]
    path = OUT / 'worked_select.svg'
    path.write_text('\n'.join(parts), encoding='utf-8')
    _svg_to_png(path, OUT / 'worked_select.png', width=640)
    return path


def run_all():
    print('Loading modules…', flush=True)
    mods = load_modules()
    print(f'Loaded {len(mods)} modules', flush=True)
    summary, by_src, sets = compute_overlap(mods)
    OUT.mkdir(parents=True, exist_ok=True)
    stats_path = OUT / 'atlas_overlap.json'
    stats_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print('Wrote', stats_path, flush=True)
    for sk in SOURCE_ORDER:
        info = summary['by_source'][sk]
        print(f"  {sk}: modules={info['modules']} unique={info['unique_reaction_sets']} nested={info['nested_in_larger']}", flush=True)
    print(f"global unique sets={summary['unique_reaction_sets_global']} nested_modules={summary['modules_nested_in_larger_global']} cross_sets={summary['cross_source_identical_sets']}", flush=True)

    print('Writing example grid…', flush=True)
    picked = pick_examples(by_src)
    ex = write_examples(picked, summary)
    print('Wrote', ex, flush=True)

    print('Writing overview mosaic (this may take a minute)…', flush=True)
    ov = write_overview_mosaic(by_src, summary)
    print('Wrote', ov, 'size', ov.stat().st_size, flush=True)

    print('Writing worked examples…', flush=True)
    wp = write_worked_examples(by_src, mods, sets)
    for k, p in wp.items():
        print(f'  {k}: {p}', flush=True)
    return dict(stats=stats_path, overview=ov, examples=ex, worked=wp, summary=summary)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('command', nargs='?', default='all', choices=['all', 'stats', 'overview', 'examples', 'worked'])
    args = p.parse_args()
    if args.command == 'all':
        run_all()
        return
    mods = load_modules()
    summary, by_src, sets = compute_overlap(mods)
    OUT.mkdir(parents=True, exist_ok=True)
    if args.command == 'stats':
        (OUT / 'atlas_overlap.json').write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')
        print(json.dumps(summary, indent=2))
    elif args.command == 'overview':
        write_overview_mosaic(by_src, summary)
    elif args.command == 'examples':
        write_examples(pick_examples(by_src), summary)
    elif args.command == 'worked':
        write_worked_examples(by_src, mods, sets)


if __name__ == '__main__':
    main()
