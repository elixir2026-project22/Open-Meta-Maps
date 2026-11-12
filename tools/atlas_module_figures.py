#!/usr/bin/env python3
"""Reproducible Atlas figures and duplication stats for docs.

  PYTHONPATH=layout_atlas python3 tools/atlas_module_figures.py all

Writes under docs/images/atlas/:
  atlas_overview.png          — 3×3 source panels; a legible sample of unique modules per source
  atlas_examples.svg/.png     — 4×4 captioned example grid (distinct maps)
  atlas_overlap.json          — per-source counts, unique-after-dedup, nesting, pairwise overlap
  worked_*.svg/.png           — worked-example panels

Figure rules: compact rows and columns, uniform cells within a panel,
grid-aligned, labels clear of strokes, no overlapping or cropped text,
readable at README width (~900 px).
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
from atlas_browser import reaction_set, source_key  # noqa: E402

OUT = ROOT / 'docs' / 'images' / 'atlas'
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
# Short labels for tight example-grid captions
SOURCE_LABEL_SHORT = {
    'escher': 'Escher',
    'reactome': 'Reactome',
    'wikipathways': 'WikiPathways',
    'wikipathways_geometry': 'WP geometry',
    'human_maps': 'Human-maps',
    'brenda': 'BRENDA',
    'metabolismregulation': 'Metabolism Reg.',
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


def load_modules():
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
        by_size = defaultdict(list)
        for s in sorted_sets:
            by_size[len(s)].append(s)
        sizes = sorted(by_size)
        for s in sorted_sets:
            if not s:
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
                nested += local[s]
        per_source[key] = dict(
            modules=len(items),
            unique_reaction_sets=len(local),
            exact_duplicate_extra=sum(c - 1 for c in local.values() if c > 1),
            nested_in_larger=nested,
        )

    all_sets = set(sets.values())
    cross_shared_sets = 0
    pairwise = {a: {b: 0 for b in SOURCE_ORDER} for a in SOURCE_ORDER}
    also_in_examples = []
    for rs, owners in global_set_to_ids.items():
        srcs = {s for s, _ in owners}
        if len(srcs) >= 2:
            cross_shared_sets += 1
            src_list = sorted(srcs)
            for i, a in enumerate(src_list):
                for b in src_list[i + 1:]:
                    pairwise[a][b] += 1
                    pairwise[b][a] += 1
            if len(rs) >= 4 and len(also_in_examples) < 20:
                also_in_examples.append(dict(
                    n_reactions=len(rs),
                    sources=src_list,
                    module_ids=[mid for _, mid in owners[:6]],
                ))

    nested_global = 0
    uniq_sets_list = sorted(all_sets, key=len)
    by_size = defaultdict(list)
    for s in uniq_sets_list:
        by_size[len(s)].append(s)
    sizes = sorted(by_size)
    nested_set_flags = {}
    for s in uniq_sets_list:
        contained = False
        if s:
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
        also_in_examples=also_in_examples,
        dedup_key='canonical reaction-id set; geometry_fingerprint for WikiPathways geometry-only',
    )
    return summary, by_src, sets


def _bbox(nodes):
    xs = [n['xy'][0] for n in nodes]
    ys = [n['xy'][1] for n in nodes]
    return min(xs), min(ys), max(xs), max(ys)


def module_density(m):
    nodes = m.get('nodes') or []
    if len(nodes) < 2:
        return 0.0
    xmin, ymin, xmax, ymax = _bbox(nodes)
    area = max((xmax - xmin) * (ymax - ymin), 1e-9)
    return len(nodes) / area


def draw_module_paths(m, x0, y0, cell, stroke='#111111', lw=1.0, node_r=1.6, reaction_fill='#111111', met_fill='#111111', met_open=True):
    """SVG fragments; curated style = black edges. Metabolites as open circles, reactions filled."""
    nodes = m.get('nodes') or []
    edges = m.get('edges') or []
    if len(nodes) < 2:
        return []
    xmin, ymin, xmax, ymax = _bbox(nodes)
    w = max(xmax - xmin, 1e-6)
    h = max(ymax - ymin, 1e-6)
    pad = cell * 0.14
    scale = (cell - 2 * pad) / max(w, h)
    cx = x0 + cell / 2
    cy = y0 + cell / 2
    mx = (xmin + xmax) / 2
    my = (ymin + ymax) / 2
    pos = {n['key']: (cx + (n['xy'][0] - mx) * scale, cy + (n['xy'][1] - my) * scale) for n in nodes}
    parts = []
    for e in edges:
        a, b = pos.get(e['a']), pos.get(e['b'])
        if a and b:
            parts.append(f'<line x1="{a[0]:.2f}" y1="{a[1]:.2f}" x2="{b[0]:.2f}" y2="{b[1]:.2f}" stroke="{stroke}" stroke-width="{lw}" stroke-linecap="round"/>')
    for n in nodes:
        x, y = pos[n['key']]
        is_r = n.get('role') == 'reaction' or str(n.get('key', '')).startswith('r:')
        if is_r:
            parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{node_r}" fill="{reaction_fill}"/>')
        elif met_open:
            parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{node_r}" fill="#ffffff" stroke="{met_fill}" stroke-width="{max(lw, 0.8)}"/>')
        else:
            parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{node_r}" fill="{met_fill}"/>')
    return parts


def _hex_rgb(h):
    h = h.lstrip('#')
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))



def unique_representatives(by_src, sets):
    """One module per unique reaction set within each source (prefer denser / larger)."""
    out = {}
    for sk, items in by_src.items():
        best = {}
        for m in items:
            rs = sets[m['id']]
            score = (module_density(m), m.get('size') or 0, len(m.get('nodes') or []))
            cur = best.get(rs)
            if cur is None or score > cur[0]:
                best[rs] = (score, m)
        out[sk] = [t[1] for t in best.values()]
    return out


def _drawable(m):
    nodes = m.get('nodes') or []
    edges = m.get('edges') or []
    return len(nodes) >= 2 and bool(edges)


def _module_features(m):
    """Size, node count, density and log aspect — used only for sampling distance."""
    nodes = m.get('nodes') or []
    size = float(m.get('size') or 0)
    n = float(len(nodes))
    dens = module_density(m)
    if len(nodes) >= 2:
        xmin, ymin, xmax, ymax = _bbox(nodes)
        aspect = math.log((xmax - xmin + 1e-9) / (ymax - ymin + 1e-9))
    else:
        aspect = 0.0
    return (size, n, dens, aspect)


def _allocate_even(k, n):
    """Spread k slots over n strata; leftover slots go to the later (larger-size) strata."""
    if n <= 0:
        return []
    base, rem = divmod(k, n)
    return [base + (1 if i >= n - rem else 0) for i in range(n)]


def _farthest_point(items, k):
    """Deterministic farthest-point sample in z-scored feature space."""
    if k <= 0:
        return []
    if len(items) <= k:
        return list(items)
    arr = np.array([_module_features(m) for m in items], dtype=float)
    std = arr.std(axis=0)
    std[std < 1e-12] = 1.0
    arr = (arr - arr.mean(axis=0)) / std
    order = sorted(range(len(items)), key=lambda i: (items[i].get('size') or 0, items[i]['id']))
    first = order[len(order) // 2]
    chosen = [first]
    remaining = set(range(len(items))) - {first}
    dmin = {i: float(np.linalg.norm(arr[i] - arr[first])) for i in remaining}
    while len(chosen) < k and remaining:
        nxt = max(remaining, key=lambda i: (dmin[i], items[i]['id']))
        chosen.append(nxt)
        remaining.remove(nxt)
        for i in remaining:
            d = float(np.linalg.norm(arr[i] - arr[nxt]))
            if d < dmin[i]:
                dmin[i] = d
    picked = [items[i] for i in chosen]
    return sorted(picked, key=lambda m: (-(m.get('size') or 0), m['id']))


def size_stratified_farthest_point(items, k):
    """Up to k drawable modules, equally from size tertiles, farthest-point inside each tertile."""
    items = [m for m in items if _drawable(m)]
    if len(items) <= k:
        return sorted(items, key=lambda m: (-(m.get('size') or 0), m['id']))
    ordered = sorted(items, key=lambda m: (m.get('size') or 0, m['id']))
    n = len(ordered)
    strata = []
    for i in range(3):
        a = int(round(i * n / 3))
        b = int(round((i + 1) * n / 3))
        if b > a:
            strata.append(ordered[a:b])
    selected = []
    for group, ck in zip(strata, _allocate_even(k, len(strata))):
        selected.extend(_farthest_point(group, ck))
    return sorted(selected, key=lambda m: (-(m.get('size') or 0), m['id']))


def _pil_draw_module(draw, m, x0, y0, cell, lw=None, node_r=None):
    """Draw one module into a square cell. Stroke and node size scale with the cell."""
    ink = (17, 17, 17)
    nodes = m.get('nodes') or []
    edges = m.get('edges') or []
    if len(nodes) < 2 or not edges:
        return
    xmin, ymin, xmax, ymax = _bbox(nodes)
    w = max(xmax - xmin, 1e-6)
    h = max(ymax - ymin, 1e-6)
    pad = cell * 0.12
    scale = (cell - 2 * pad) / max(w, h)
    cx = x0 + cell / 2
    cy = y0 + cell / 2
    mx = (xmin + xmax) / 2
    my = (ymin + ymax) / 2
    pos = {n['key']: (cx + (n['xy'][0] - mx) * scale, cy + (n['xy'][1] - my) * scale) for n in nodes}
    if lw is None:
        lw = max(2, int(round(cell * 0.015)))
    r_r = node_r if node_r is not None else max(2.2, cell * 0.032)
    r_m = max(2.0, cell * 0.026)
    for e in edges:
        a, b = pos.get(e['a']), pos.get(e['b'])
        if a and b:
            draw.line([a, b], fill=ink, width=lw)
    for n in nodes:
        x, y = pos[n['key']]
        is_r = n.get('role') == 'reaction' or str(n.get('key', '')).startswith('r:')
        if is_r:
            draw.ellipse([x - r_r, y - r_r, x + r_r, y + r_r], fill=ink)
        else:
            draw.ellipse([x - r_m, y - r_m, x + r_m, y + r_m], outline=ink, fill=(255, 255, 255), width=max(1, lw))


def write_overview_mosaic(by_src, summary, sets):
    """Equal 3×3 source panels with a size-stratified sample of unique modules, readable at README width."""
    from PIL import Image, ImageDraw, ImageFont

    final_w, final_h = 3600, 2700
    hi = 2  # draw at 2×, then downsample
    canvas_w, canvas_h = final_w * hi, final_h * hi
    pad = 40
    gap = 24
    title_h = 88
    legend_h = 96
    sample_k = 24
    cell_min = 160 * hi
    cols, rows = 6, 4
    cell_gap = 10

    try:
        font_title = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf', 44)
        font_legend = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf', 36)
    except Exception:
        font_title = font_legend = ImageFont.load_default()

    reps = unique_representatives(by_src, sets)
    blocks = []
    for sk in SOURCE_ORDER:
        items = reps.get(sk) or []
        info = summary['by_source'][sk]
        shown = size_stratified_farthest_point(items, sample_k)
        blocks.append(dict(sk=sk, items=shown, info=info, n_drawable=sum(1 for m in items if _drawable(m))))

    assert len(blocks) == 9, f'overview expects 9 sources, got {len(blocks)}'

    grid_w = canvas_w - 2 * pad
    grid_h = canvas_h - 2 * pad - legend_h
    panel_w = (grid_w - 2 * gap) // 3
    panel_h = (grid_h - 2 * gap) // 3
    inner_w = panel_w - 16
    inner_h = panel_h - title_h - 12
    cell = min((inner_w - (cols - 1) * cell_gap) // cols, (inner_h - (rows - 1) * cell_gap) // rows)
    if cell < cell_min:
        raise RuntimeError(f'overview cell {cell}px is below the {cell_min}px minimum at {hi}×')
    grid_pw = cols * cell + (cols - 1) * cell_gap
    grid_ph = rows * cell + (rows - 1) * cell_gap

    im = Image.new('RGB', (canvas_w, canvas_h), (255, 255, 255))
    draw = ImageDraw.Draw(im)

    for i, b in enumerate(blocks):
        pr, pc = divmod(i, 3)
        x0 = pad + pc * (panel_w + gap)
        y0 = pad + pr * (panel_h + gap)
        tint = _hex_rgb(SOURCE_TINT[b['sk']])
        stroke = _hex_rgb(SOURCE_STROKE[b['sk']])
        draw.rectangle([x0, y0, x0 + panel_w - 1, y0 + panel_h - 1], fill=tint, outline=stroke, width=4)
        info = b['info']
        label = (
            f"{SOURCE_LABEL.get(b['sk'], b['sk'])}  ·  "
            f"{info['modules']:,} cut · {info['unique_reaction_sets']:,} unique · "
            f"{len(b['items']):,} shown"
        )
        # Shrink the title font if the full counts would overflow the panel.
        title_font = font_title
        if draw.textlength(label, font=title_font) > panel_w - 28:
            try:
                title_font = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf', 36)
            except Exception:
                title_font = font_title
        if draw.textlength(label, font=title_font) > panel_w - 28:
            label = (
                f"{SOURCE_LABEL_SHORT.get(b['sk'], b['sk'])}  ·  "
                f"{info['modules']:,} cut · {info['unique_reaction_sets']:,} unique · "
                f"{len(b['items']):,} shown"
            )
        # Centre the title; keep it clear of the frame.
        tx = x0 + (panel_w - draw.textlength(label, font=title_font)) / 2
        draw.text((tx, y0 + 18), label, fill=(17, 17, 17), font=title_font)

        gx = x0 + (panel_w - grid_pw) // 2
        gy = y0 + title_h + (inner_h - grid_ph) // 2
        for j, m in enumerate(b['items']):
            r, c = divmod(j, cols)
            if r >= rows:
                break
            _pil_draw_module(draw, m, gx + c * (cell + cell_gap), gy + r * (cell + cell_gap), cell)

    legend = (
        f"A deterministic size-stratified sample of unique reaction sets in each source  ·  "
        f"cut {summary['total_modules']:,}  ·  unique {summary['unique_reaction_sets_global']:,}  ·  "
        f"up to {sample_k} modules shown per source"
    )
    lw = draw.textlength(legend, font=font_legend)
    draw.text(((canvas_w - lw) / 2, canvas_h - legend_h + 24), legend, fill=(31, 41, 55), font=font_legend)

    out = im.resize((final_w, final_h), Image.Resampling.LANCZOS)
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / 'atlas_overview.png'
    out.save(path, format='PNG', optimize=True)
    _shrink_png(path, max_bytes=5_000_000)
    return path


def _shrink_png(path: Path, max_bytes: int):
    from PIL import Image
    if path.stat().st_size <= max_bytes:
        return
    im = Image.open(path).convert('RGB')
    for scale in (0.9, 0.8, 0.7, 0.6):
        w, h = im.size
        small = im.resize((max(1, int(w * scale)), max(1, int(h * scale))), Image.Resampling.LANCZOS)
        small.save(path, format='PNG', optimize=True)
        if path.stat().st_size <= max_bytes:
            return


def _ellipsis(text, max_chars):
    text = text or ''
    if len(text) <= max_chars:
        return text
    return text[: max(1, max_chars - 1)].rstrip(' …') + '…'


def _wrap(text, max_chars, max_lines=2):
    """Word-wrap; truncated text always ends the last line with an ellipsis."""
    text = (text or '').strip()
    if not text:
        return []
    words = text.split()
    lines = []
    idx = 0
    while idx < len(words) and len(lines) < max_lines:
        cur = words[idx]
        idx += 1
        while idx < len(words):
            trial = f'{cur} {words[idx]}'
            if len(trial) <= max_chars:
                cur = trial
                idx += 1
            else:
                break
        if len(cur) > max_chars:
            cur = _ellipsis(cur, max_chars)
        lines.append(cur)
    if idx < len(words) and lines:
        # Force ellipsis even when the last line already fills max_chars.
        base = lines[-1].rstrip('…')
        lines[-1] = base[: max(1, max_chars - 1)].rstrip(' …') + '…'
    return lines


def pick_examples(by_src):
    """16 distinct maps (4×4), prefer size 6–12, one map per cell."""
    targets = [
        ('escher', ['core metabolism', 'glycolysis'], 8),
        ('escher', ['tca', 'citric', 'central carbon'], 8),
        ('escher', ['fatty', 'lipid'], 12),
        ('reactome', ['glycolysis', 'glucose'], 8),
        ('reactome', ['aspartate', 'nucleotide', 'amino'], 8),
        ('reactome', None, 8),
        ('wikipathways', ['glycolysis', 'tca', 'pentose'], 8),
        ('wikipathways', None, 6),
        ('wikipathways_geometry', None, 8),
        ('human_maps', ['glycolysis'], 8),
        ('human_maps', ['cholesterol', 'lipid', 'sphingolipid'], 8),
        ('ich360', ['amino', 'biosynthesis'], 8),
        ('ich360', ['nucleotide', 'cofactor'], 6),
        ('brenda', None, 8),
        ('metabolismregulation', ['glycolysis'], 6),
        ('panther', None, 4),
    ]
    picked = []
    used_maps = set()
    used_ids = set()
    for sk, substrs, size in targets:
        items = by_src.get(sk) or []
        cand = [
            m for m in items
            if m['id'] not in used_ids
            and (m.get('map_id'), m.get('name')) not in used_maps
            and 4 <= (m.get('size') or 0) <= 24
            and len(m.get('nodes') or []) >= 5
        ]
        if substrs:
            named = [m for m in cand if any(s in (m.get('name') or '').lower() for s in substrs)]
            if named:
                cand = named
        prefer = [m for m in cand if m.get('size') == size]
        pool = prefer or cand
        pool = sorted(pool, key=lambda m: (-(m.get('size') == size), -min(m.get('size') or 0, 12), m.get('name') or '', m['id']))
        if not pool:
            # relax size floor
            pool = sorted(
                [m for m in items if m['id'] not in used_ids and (m.get('map_id'), m.get('name')) not in used_maps and (m.get('size') or 0) >= 2],
                key=lambda m: (-(m.get('size') or 0), m['id']),
            )
        if not pool:
            continue
        m = pool[0]
        used_ids.add(m['id'])
        used_maps.add((m.get('map_id'), m.get('name')))
        picked.append(m)
    return picked[:16]


def write_examples(picked, summary):
    cols = 4
    rows = int(math.ceil(len(picked) / cols))
    cell = 148
    cap_h = 52
    gap = 14
    left = 18
    top = 18
    W = left * 2 + cols * cell + (cols - 1) * gap
    H = top * 2 + rows * (cell + cap_h) + (rows - 1) * gap
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="-apple-system,Segoe UI,Helvetica,Arial,sans-serif">',
        '<title>Example Atlas modules</title>',
        f'<rect width="{W}" height="{H}" fill="#ffffff"/>',
    ]
    for i, m in enumerate(picked):
        r, c = divmod(i, cols)
        x = left + c * (cell + gap)
        y = top + r * (cell + cap_h + gap)
        sk = source_key(m)
        parts.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" rx="8" fill="#fafafa" stroke="{SOURCE_STROKE.get(sk, "#888")}" stroke-width="1.2"/>')
        parts.extend(draw_module_paths(m, x, y, cell, stroke='#111', lw=1.05, node_r=2.0))
        tx = x + cell / 2
        max_chars = 22
        parts.append(f'<text x="{tx:.1f}" y="{y + cell + 14}" text-anchor="middle" font-size="11" font-weight="600" fill="#111">{html.escape(SOURCE_LABEL_SHORT.get(sk, sk))}</text>')
        name_lines = _wrap(m.get('name') or m.get('map_id') or '', max_chars, 2)
        for li, line in enumerate(name_lines):
            parts.append(f'<text x="{tx:.1f}" y="{y + cell + 28 + li * 12}" text-anchor="middle" font-size="10" fill="#374151">{html.escape(line)}</text>')
        parts.append(f'<text x="{tx:.1f}" y="{y + cell + 28 + len(name_lines) * 12}" text-anchor="middle" font-size="10" fill="#6b7280">{m.get("size", len(m.get("reactions") or []))} reactions</text>')
    parts.append('</svg>')
    OUT.mkdir(parents=True, exist_ok=True)
    svg_path = OUT / 'atlas_examples.svg'
    svg_path.write_text('\n'.join(parts), encoding='utf-8')
    png_path = OUT / 'atlas_examples.png'
    _svg_to_png(svg_path, png_path, width=900)
    return svg_path, png_path


def _svg_to_png(svg_path: Path, png_path: Path, width: int = 800):
    try:
        import cairosvg
        cairosvg.svg2png(url=svg_path.resolve().as_uri(), write_to=str(png_path), output_width=width)
    except Exception as exc:
        print(f'PNG preview skipped for {svg_path.name}: {exc}', flush=True)


def write_worked_examples(by_src, mods, sets):
    OUT.mkdir(parents=True, exist_ok=True)
    paths = {}
    paths['cut'] = _worked_cut_figure(by_src)
    paths['ids'] = _worked_ids_figure(by_src)
    paths['dedup'] = _worked_dedup_figure(by_src, mods, sets)
    paths['prov'] = _worked_provenance_figure(by_src, mods, sets)
    paths['select'] = _worked_select_figure(by_src, sets)
    return paths


def _worked_cut_figure(by_src):
    """Source-map region with cut-window outlines + 3–4 non-trivial modules."""
    from database import escher_diagram
    from escher_cutter import reaction_windows, raw_escher_paths

    # Prefer e_coli_core — small, readable
    paths = [p for p in raw_escher_paths() if 'e_coli_core' in str(p)]
    if not paths:
        paths = raw_escher_paths()[:1]
    path = paths[0]
    diagram = escher_diagram(path)
    windows = reaction_windows(diagram)
    # Build bboxes in original coordinates for size 4/8/12 windows
    rx = {r['key']: r for r in diagram['events']}
    glyphs = diagram['glyphs']

    def window_bbox(rs):
        pts = []
        for k in rs:
            if k in rx:
                pts.append(rx[k]['xy'])
            for mk in (rx[k].get('visible') or {}):
                g = glyphs.get(mk) or {}
                if g.get('primary', True) and 'xy' in g:
                    pts.append(g['xy'])
        if not pts:
            return None
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        return min(xs), min(ys), max(xs), max(ys)

    # Map bounds first (used for spatial separation of cut windows).
    all_xy = [r['xy'] for r in diagram['events']] + [g['xy'] for g in glyphs.values() if 'xy' in g]
    xmin, ymin = min(p[0] for p in all_xy), min(p[1] for p in all_xy)
    xmax, ymax = max(p[0] for p in all_xy), max(p[1] for p in all_xy)

    # Pick 3–4 non-trivial windows with distinct sizes and spatially separated bboxes.
    sized = defaultdict(list)
    for rs in windows:
        if 4 <= len(rs) <= 16:
            bb = window_bbox(rs)
            if bb and (bb[2] - bb[0]) > 1 and (bb[3] - bb[1]) > 1:
                sized[len(rs)].append((rs, bb))

    def bb_iou(a, b):
        ix0, iy0 = max(a[0], b[0]), max(a[1], b[1])
        ix1, iy1 = min(a[2], b[2]), min(a[3], b[3])
        iw, ih = max(0.0, ix1 - ix0), max(0.0, iy1 - iy0)
        inter = iw * ih
        if inter <= 0:
            return 0.0
        aa = max(a[2] - a[0], 1e-6) * max(a[3] - a[1], 1e-6)
        ba = max(b[2] - b[0], 1e-6) * max(b[3] - b[1], 1e-6)
        return inter / max(aa + ba - inter, 1e-6)

    def bb_center(bb):
        return ((bb[0] + bb[2]) / 2, (bb[1] + bb[3]) / 2)

    def distinct_enough(bb, chosen_bbs, max_iou=0.15, min_frac=0.14):
        if not chosen_bbs:
            return True
        cx, cy = bb_center(bb)
        span = max(xmax - xmin, ymax - ymin, 1e-6)
        for ob in chosen_bbs:
            if bb_iou(bb, ob) > max_iou:
                return False
            ox, oy = bb_center(ob)
            if math.hypot(cx - ox, cy - oy) / span < min_frac:
                return False
        return True

    chosen = []
    chosen_bbs = []
    for n in (4, 8, 12):
        cands = sorted(
            sized.get(n) or [],
            key=lambda t: -((t[1][2] - t[1][0]) * (t[1][3] - t[1][1])),
        )
        for item in cands:
            if distinct_enough(item[1], chosen_bbs):
                chosen.append(item)
                chosen_bbs.append(item[1])
                break
        if len(chosen) >= 3:
            break
    if len(chosen) < 3:
        # Relax IoU, still prefer separated centres.
        for max_iou, min_frac in ((0.25, 0.10), (0.40, 0.06), (0.55, 0.0)):
            for n in (4, 8, 12, 6, 10, 16):
                for item in sorted(sized.get(n) or [], key=lambda t: -((t[1][2] - t[1][0]) * (t[1][3] - t[1][1]))):
                    if item in chosen:
                        continue
                    if distinct_enough(item[1], chosen_bbs, max_iou=max_iou, min_frac=min_frac):
                        chosen.append(item)
                        chosen_bbs.append(item[1])
                    if len(chosen) >= 3:
                        break
                if len(chosen) >= 3:
                    break
            if len(chosen) >= 3:
                break
    chosen = chosen[:3]

    def window_as_module(rs):
        nodes, edges, seen = [], [], set()
        for k in rs:
            if k not in rx:
                continue
            r = rx[k]
            rkey = f"r:{r['id']}"
            nodes.append({'key': rkey, 'xy': r['xy'], 'role': 'reaction'})
            for mk in (r.get('visible') or {}):
                g = glyphs.get(mk) or {}
                if 'xy' not in g:
                    continue
                mkey = f'm:{mk}'
                if mkey not in seen:
                    nodes.append({'key': mkey, 'xy': g['xy'], 'role': 'metabolite'})
                    seen.add(mkey)
                edges.append({'a': rkey, 'b': mkey})
        return {
            'nodes': nodes,
            'edges': edges,
            'size': len(rs),
            'reactions': [{'id': rx[k]['id']} for k in rs if k in rx],
        }

    # Layout: title + legend (outside map) + map panel + matching modules.
    map_w, map_h = 360, 280
    mw = max(xmax - xmin, 1e-6)
    mh = max(ymax - ymin, 1e-6)
    scale = min((map_w - 24) / mw, (map_h - 24) / mh)
    legend_y = 48
    map_top = 72

    def tx(xy):
        return 20 + 12 + (xy[0] - xmin) * scale, map_top + 12 + (xy[1] - ymin) * scale

    cell = 118
    gap = 14
    right_x = 20 + map_w + 28
    W = right_x + max(len(chosen), 1) * (cell + gap) + 8
    H = map_top + max(map_h, cell + 40) + 36
    # Stroke outlines with very light tint (<=15% alpha) so the map stays visible.
    fill_colors = ['rgba(63,143,107,0.12)', 'rgba(185,133,47,0.12)', 'rgba(74,111,165,0.12)']
    stroke_colors = ['#3f8f6b', '#b9852f', '#4a6fa5']

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="-apple-system,Segoe UI,Helvetica,Arial,sans-serif">',
        f'<rect width="{W}" height="{H}" fill="#fff"/>',
        f'<text x="20" y="24" font-size="13" font-weight="700" fill="#111">Cut windows on a source map</text>',
        f'<text x="20" y="40" font-size="11" fill="#4b5563">{html.escape(path.stem)} · outlined windows to stored modules</text>',
    ]
    # Legend outside the map panel
    chip_x = 20
    for i, (rs, bb) in enumerate(chosen):
        parts.append(
            f'<rect x="{chip_x}" y="{legend_y}" width="54" height="16" rx="3" '
            f'fill="#fff" stroke="{stroke_colors[i % 3]}" stroke-width="1.4"/>'
        )
        parts.append(
            f'<text x="{chip_x + 27:.1f}" y="{legend_y + 12}" text-anchor="middle" '
            f'font-size="10" font-weight="600" fill="{stroke_colors[i % 3]}">{len(rs)} rx</text>'
        )
        chip_x += 62

    parts.append(f'<rect x="20" y="{map_top}" width="{map_w}" height="{map_h}" rx="8" fill="#fafafa" stroke="#d0d7de"/>')
    # faint full map
    for r in diagram['events']:
        for mk, coeff in (r.get('visible') or {}).items():
            g = glyphs.get(mk) or {}
            if 'xy' not in g:
                continue
            a, b = tx(r['xy']), tx(g['xy'])
            parts.append(f'<line x1="{a[0]:.1f}" y1="{a[1]:.1f}" x2="{b[0]:.1f}" y2="{b[1]:.1f}" stroke="#cbd5e1" stroke-width="0.7"/>')
    for r in diagram['events']:
        x, y = tx(r['xy'])
        parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2.0" fill="#111"/>')
    # Window outlines around actual module node bboxes (largest first).
    for i, (rs, bb) in sorted(enumerate(chosen), key=lambda t: -len(t[1][0])):
        x0, y0 = tx((bb[0], bb[1]))
        x1, y1 = tx((bb[2], bb[3]))
        pad = 6
        parts.append(
            f'<rect x="{min(x0,x1)-pad:.1f}" y="{min(y0,y1)-pad:.1f}" width="{abs(x1-x0)+2*pad:.1f}" height="{abs(y1-y0)+2*pad:.1f}" '
            f'fill="{fill_colors[i % 3]}" stroke="{stroke_colors[i % 3]}" stroke-width="2" rx="3"/>'
        )

    # Right panels draw the same window geometry as the outlines (exact match).
    for i, (rs, bb) in enumerate(chosen):
        match = window_as_module(rs)
        x = right_x + i * (cell + gap)
        y = map_top + 8
        parts.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" rx="8" fill="#fff" stroke="{stroke_colors[i % 3]}" stroke-width="1.3"/>')
        parts.extend(draw_module_paths(match, x, y, cell, stroke='#111', lw=1.0, node_r=1.8))
        parts.append(f'<text x="{x+cell/2:.1f}" y="{y+cell+14}" text-anchor="middle" font-size="10" fill="#111">{len(rs)} reactions</text>')

    parts.append(
        f'<text x="20" y="{H-12}" font-size="11" fill="#6b7280">'
        f'Left: source map with cut-window outlines. Right: resulting non-trivial modules '
        f'(sizes {", ".join(str(len(rs)) for rs, _ in chosen)}).</text>'
    )
    parts.append('</svg>')
    out = OUT / 'worked_cut.svg'
    out.write_text('\n'.join(parts), encoding='utf-8')
    _svg_to_png(out, OUT / 'worked_cut.png', width=980)
    return out


def _worked_ids_figure(by_src):
    hm_cands = [
        m for m in (by_src.get('human_maps') or [])
        if sum(1 for r in (m.get('reactions') or []) if (r.get('references') or {}).get('bigg.reaction')) >= 3
    ]
    m = max(hm_cands, key=lambda x: x.get('size') or 0) if hm_cands else None
    W = 740
    headers = ['Native (human.gem)', 'bigg.reaction', 'metanetx.reaction', 'rhea']
    xs = [20, 200, 380, 560]
    rows = []
    if m:
        for r in (m.get('reactions') or []):
            refs = r.get('references') or {}
            if refs.get('bigg.reaction') or refs.get('metanetx.reaction') or refs.get('rhea'):
                rows.append((r.get('id'), refs))
            if len(rows) >= 4:
                break
    body = [
        f'<text x="20" y="24" font-size="13" font-weight="700" fill="#111">Canonical reaction identifiers (Human-maps)</text>',
    ]
    if not m:
        body.append('<text x="20" y="48" font-size="12" fill="#6b7280">No example module found.</text>')
        H = 64
    else:
        body.append(f'<text x="20" y="42" font-size="11" fill="#6b7280">Map: {html.escape(_ellipsis(m.get("name") or "", 70))}</text>')
        y = 64
        for x, h in zip(xs, headers):
            body.append(f'<text x="{x}" y="{y}" font-size="11" font-weight="600" fill="#111">{h}</text>')
        y += 6
        body.append(f'<line x1="20" y1="{y}" x2="{W-20}" y2="{y}" stroke="#d0d7de"/>')
        y += 18
        for rid, refs in rows:
            vals = [
                (refs.get('human.gem') or [rid])[0],
                (refs.get('bigg.reaction') or ['—'])[0],
                (refs.get('metanetx.reaction') or ['—'])[0],
                (refs.get('rhea') or ['—'])[0],
            ]
            for x, v in zip(xs, vals):
                body.append(f'<text x="{x}" y="{y}" font-size="11" fill="#374151">{html.escape(str(v))}</text>')
            y += 20
        H = y + 4
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="-apple-system,Segoe UI,Helvetica,Arial,sans-serif">',
        f'<rect width="{W}" height="{H}" fill="#fff"/>',
        *body,
        '</svg>',
    ]
    path = OUT / 'worked_ids.svg'
    path.write_text('\n'.join(parts), encoding='utf-8')
    _svg_to_png(path, OUT / 'worked_ids.png', width=740)
    return path


def _worked_dedup_figure(by_src, mods, sets):
    """Exact cross-source duplicate (merged also_in record) + nesting case."""
    # Cross-source identical set with size >= 6
    owners_by_set = defaultdict(list)
    for m in mods:
        rs = sets[m['id']]
        if len(rs) >= 6:
            owners_by_set[rs].append(m)
    pair = None
    for rs, owners in owners_by_set.items():
        srcs = {source_key(m) for m in owners}
        if len(srcs) >= 2 and 'escher' in srcs and 'human_maps' in srcs:
            a = next(m for m in owners if source_key(m) == 'escher')
            b = next(m for m in owners if source_key(m) == 'human_maps')
            pair = (a, b, rs)
            break
    if pair is None:
        for rs, owners in owners_by_set.items():
            srcs = {source_key(m) for m in owners}
            if len(srcs) >= 2:
                aa = sorted(srcs)
                a = next(m for m in owners if source_key(m) == aa[0])
                b = next(m for m in owners if source_key(m) == aa[1])
                pair = (a, b, rs)
                break

    # Nesting within Escher same map
    nest = None
    escher = by_src.get('escher') or []
    by_map = defaultdict(list)
    for m in escher:
        by_map[m.get('map_id')].append(m)
    for mid, group in by_map.items():
        sets_local = {m['id']: sets[m['id']] for m in group}
        for child in group:
            cs = sets_local[child['id']]
            if not (4 <= len(cs) <= 8):
                continue
            for parent in group:
                if parent['id'] == child['id']:
                    continue
                ps = sets_local[parent['id']]
                if len(ps) >= len(cs) + 4 and cs < ps:
                    nest = (child, parent)
                    break
            if nest:
                break
        if nest:
            break

    cell = 120
    W = 920
    body = [
        f'<text x="20" y="24" font-size="13" font-weight="700" fill="#111">Dedup: identical sets merge; nested modules remain stored</text>',
        '<text x="20" y="48" font-size="12" font-weight="600" fill="#111">A. Same reaction set across sources (drawings may differ)</text>',
    ]
    # Left panel: exact duplicate → merged record
    if pair:
        a, b, rs = pair
        for i, m in enumerate((a, b)):
            x = 20 + i * (cell + 16)
            y = 60
            sk = source_key(m)
            body.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" rx="8" fill="#fafafa" stroke="{SOURCE_STROKE.get(sk,"#888")}" stroke-width="1.2"/>')
            body.extend(draw_module_paths(m, x, y, cell, stroke='#111', lw=1.0, node_r=1.7))
            body.append(f'<text x="{x+cell/2:.1f}" y="{y+cell+14}" text-anchor="middle" font-size="10" font-weight="600" fill="#111">{html.escape(_ellipsis(SOURCE_LABEL_SHORT.get(sk, sk), 16))}</text>')
            body.append(f'<text x="{x+cell/2:.1f}" y="{y+cell+28}" text-anchor="middle" font-size="9" fill="#6b7280">{html.escape(_ellipsis(m.get("name") or "", 18))}</text>')
        body.append('<path d="M280 120 L310 120 L310 114 L324 124 L310 134 L310 128 L280 128 Z" fill="#6b7280"/>')
        mx, my, mw, mh = 340, 70, 200, 110
        body.append(f'<rect x="{mx}" y="{my}" width="{mw}" height="{mh}" rx="8" fill="#e7f4ee" stroke="#3f8f6b" stroke-width="1.4"/>')
        body.append(f'<text x="{mx+10}" y="{my+22}" font-size="11" font-weight="700" fill="#111">Unique after dedup</text>')
        body.append(f'<text x="{mx+10}" y="{my+42}" font-size="10" fill="#374151">{len(rs)} canonical reactions</text>')
        body.append(f'<text x="{mx+10}" y="{my+60}" font-size="10" fill="#374151">also_in:</text>')
        body.append(f'<text x="{mx+10}" y="{my+76}" font-size="10" fill="#374151">{html.escape(SOURCE_LABEL_SHORT.get(source_key(a), source_key(a)))}</text>')
        body.append(f'<text x="{mx+10}" y="{my+92}" font-size="10" fill="#374151">{html.escape(SOURCE_LABEL_SHORT.get(source_key(b), source_key(b)))}</text>')
    # Right panel: nesting
    body.append('<text x="560" y="48" font-size="12" font-weight="600" fill="#111">B. Nested module</text>')
    if nest:
        child, parent = nest
        for i, m in enumerate((child, parent)):
            x = 560 + i * (cell + 16)
            y = 60
            sk = source_key(m)
            stroke = '#b9852f' if i == 0 else '#3f8f6b'
            body.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" rx="8" fill="#fafafa" stroke="{stroke}" stroke-width="1.2"/>')
            body.extend(draw_module_paths(m, x, y, cell, stroke='#111', lw=1.0, node_r=1.7))
            label = 'nested module' if i == 0 else 'larger window'
            body.append(f'<text x="{x+cell/2:.1f}" y="{y+cell+14}" text-anchor="middle" font-size="10" font-weight="600" fill="#111">{label}</text>')
            body.append(f'<text x="{x+cell/2:.1f}" y="{y+cell+28}" text-anchor="middle" font-size="9" fill="#6b7280">{m.get("size")} reactions</text>')
        body.append('<path d="M686 120 L716 120 L716 114 L730 124 L716 134 L716 128 L686 128 Z" fill="#6b7280"/>')

    cap = (
        'A: identical canonical reaction sets collapse to one unique entry with also_in provenance. '
        'B: smaller reaction sets nested in larger windows stay stored as multi-scale nested modules.'
    )
    cap_lines = _wrap(cap, 110, 3)
    y_cap = 60 + cell + 28 + 18
    for i, line in enumerate(cap_lines):
        body.append(f'<text x="20" y="{y_cap + i * 14}" font-size="11" fill="#4b5563">{html.escape(line)}</text>')
    H = y_cap + max(len(cap_lines) - 1, 0) * 14 + 12
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="-apple-system,Segoe UI,Helvetica,Arial,sans-serif">',
        f'<rect width="{W}" height="{H}" fill="#fff"/>',
        *body,
        '</svg>',
    ]
    path = OUT / 'worked_dedup.svg'
    path.write_text('\n'.join(parts), encoding='utf-8')
    _svg_to_png(path, OUT / 'worked_dedup.png', width=920)
    return path


def _worked_provenance_figure(by_src, mods, sets):
    escher = by_src.get('escher') or []
    # Prefer a module that also appears in another source
    sample = None
    also_in = []
    owners = defaultdict(list)
    for m in mods:
        owners[sets[m['id']]].append(m)
    for m in escher:
        if (m.get('size') or 0) < 6 or not m.get('extraction'):
            continue
        peers = owners.get(sets[m['id']]) or []
        srcs = sorted({source_key(p) for p in peers})
        if len(srcs) >= 2:
            sample = m
            also_in = srcs
            break
    if sample is None:
        sample = next((m for m in escher if (m.get('size') or 0) >= 6 and m.get('extraction')), escher[0] if escher else None)

    W = 700
    thumb = 140
    body = [
        f'<text x="20" y="24" font-size="13" font-weight="700" font-family="-apple-system,Segoe UI,Helvetica,Arial,sans-serif" fill="#111">Provenance fields on one module</text>',
    ]
    if not sample:
        H = 48
        parts = [
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="ui-monospace,Menlo,monospace">',
            f'<rect width="{W}" height="{H}" fill="#fff"/>',
            *body,
            '</svg>',
        ]
        path = OUT / 'worked_provenance.svg'
        path.write_text('\n'.join(parts), encoding='utf-8')
        return path
    ex = sample.get('extraction') or {}
    lines = [
        f'source:   {sample.get("source")}',
        f'map:      {_ellipsis(sample.get("name") or sample.get("map_id") or "", 60)}',
        f'license:  {sample.get("license")}',
        f'method:   {ex.get("method")}',
        f'version:  {ex.get("version")}',
        f'size:     {sample.get("size")} reactions',
    ]
    if also_in:
        lines.append('also_in:  ' + ', '.join(SOURCE_LABEL.get(s, s) for s in also_in))
    y = 48
    for line in lines:
        body.append(f'<text x="28" y="{y}" font-size="12" fill="#374151">{html.escape(line)}</text>')
        y += 18
    body.append(f'<rect x="{W-170}" y="40" width="{thumb}" height="{thumb}" rx="8" fill="#fafafa" stroke="#d0d7de"/>')
    body.extend(draw_module_paths(sample, W - 170, 40, thumb, stroke='#111', lw=1.0, node_r=1.8))
    H = max(y + 8, 40 + thumb + 8)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="ui-monospace,Menlo,monospace">',
        f'<rect width="{W}" height="{H}" fill="#fff"/>',
        *body,
        '</svg>',
    ]
    path = OUT / 'worked_provenance.svg'
    path.write_text('\n'.join(parts), encoding='utf-8')
    _svg_to_png(path, OUT / 'worked_provenance.png', width=700)
    return path


def _worked_select_figure(by_src, sets):
    """Real same-source same-set modules with different density."""
    best = None
    for sk in ('escher', 'human_maps', 'ich360'):
        groups = defaultdict(list)
        for m in by_src.get(sk) or []:
            rs = sets[m['id']]
            if len(rs) >= 6:
                groups[rs].append(m)
        for rs, items in groups.items():
            if len(items) < 2:
                continue
            scored = []
            for m in items:
                d = module_density(m)
                if d > 0:
                    scored.append((d, m))
            scored.sort(key=lambda t: -t[0])
            if len(scored) >= 2 and scored[0][0] > scored[-1][0] * 1.25:
                best = (sk, scored[0], scored[-1], rs)
                break
        if best:
            break

    cell = 160
    W, H = 720, 310
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="-apple-system,Segoe UI,Helvetica,Arial,sans-serif">',
        f'<rect width="{W}" height="{H}" fill="#fff"/>',
        f'<text x="20" y="24" font-size="13" font-weight="700" fill="#111">Compact selection among overlapping candidates</text>',
    ]
    if not best:
        parts.append('<text x="20" y="60" font-size="12" fill="#6b7280">No same-set density pair found.</text></svg>')
        path = OUT / 'worked_select.svg'
        path.write_text('\n'.join(parts), encoding='utf-8')
        return path
    sk, (d_hi, m_hi), (d_lo, m_lo), rs = best
    # sprawling left, compact right
    for i, (m, dens, title, stroke) in enumerate((
        (m_lo, d_lo, 'Lower density (dropped)', '#9ca3af'),
        (m_hi, d_hi, 'Higher density (kept)', '#3f8f6b'),
    )):
        x = 40 + i * (cell + 100)
        y = 44
        fill = '#f5f5f5' if i == 0 else '#e7f4ee'
        parts.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" rx="8" fill="{fill}" stroke="{stroke}" stroke-width="1.5"/>')
        parts.extend(draw_module_paths(m, x, y, cell, stroke='#111', lw=1.05, node_r=1.9))
        parts.append(f'<text x="{x+cell/2:.1f}" y="{y+cell+16}" text-anchor="middle" font-size="11" font-weight="600" fill="#111">{html.escape(title)}</text>')
        parts.append(f'<text x="{x+cell/2:.1f}" y="{y+cell+32}" text-anchor="middle" font-size="10" fill="#4b5563">{html.escape(_ellipsis(m.get("name") or "", 28))}</text>')
        parts.append(f'<text x="{x+cell/2:.1f}" y="{y+cell+46}" text-anchor="middle" font-size="10" fill="#6b7280">{len(rs)} reactions · density {dens:.3g}</text>')
    # arrow path between boxes (ASCII geometry, no unicode glyphs)
    parts.append('<path d="M250 124 L300 124 L300 116 L318 128 L300 140 L300 132 L250 132 Z" fill="#6b7280"/>')
    parts.append(
        f'<text x="20" y="{H - 16}" font-size="11" fill="#4b5563">'
        'Same source and identical reaction set; coverage still ranks first. '
        'Density = matched nodes / bounding-box area.</text>'
    )
    parts.append('</svg>')
    path = OUT / 'worked_select.svg'
    path.write_text('\n'.join(parts), encoding='utf-8')
    _svg_to_png(path, OUT / 'worked_select.png', width=720)
    return path


def run_all():
    print('Loading modules…', flush=True)
    mods = load_modules()
    print(f'Loaded {len(mods)} modules', flush=True)
    summary, by_src, sets = compute_overlap(mods)
    OUT.mkdir(parents=True, exist_ok=True)
    stats_path = OUT / 'atlas_overlap.json'
    # JSON-serialize frozensets already expanded
    stats_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print('Wrote', stats_path, flush=True)
    for sk in SOURCE_ORDER:
        info = summary['by_source'][sk]
        print(f"  {sk}: modules={info['modules']} unique={info['unique_reaction_sets']} nested={info['nested_in_larger']}", flush=True)
    print(
        f"global unique={summary['unique_reaction_sets_global']} "
        f"nested={summary['modules_nested_in_larger_global']} cross_sets={summary['cross_source_identical_sets']}",
        flush=True,
    )

    print('Writing example grid…', flush=True)
    picked = pick_examples(by_src)
    print('  picked maps:', [(source_key(m), m.get('name'), m.get('size')) for m in picked], flush=True)
    ex = write_examples(picked, summary)
    print('Wrote', ex, flush=True)

    print('Writing overview mosaic…', flush=True)
    ov = write_overview_mosaic(by_src, summary, sets)
    print('Wrote', ov, 'size', ov.stat().st_size, 'dims', end=' ', flush=True)
    try:
        from PIL import Image
        im = Image.open(ov)
        print(im.size, flush=True)
    except Exception:
        print(flush=True)

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
        write_overview_mosaic(by_src, summary, sets)
    elif args.command == 'examples':
        write_examples(pick_examples(by_src), summary)
    elif args.command == 'worked':
        write_worked_examples(by_src, mods, sets)


if __name__ == '__main__':
    main()
