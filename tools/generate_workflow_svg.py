#!/usr/bin/env python3
"""Generate docs/images/workflow.svg — single editable source for the README diagram.

  python3 tools/generate_workflow_svg.py
  python3 tools/generate_workflow_svg.py --check   # PNG at README width (~700 px)

Boxes sit on a fixed grid; arrow labels keep a ≥12 px gap from strokes and box edges.
"""
from __future__ import annotations

import argparse
import io
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'docs' / 'images' / 'workflow.svg'

CARD_W, CARD_H = 268, 150
GAP = 24
LEFT = 22
ATLAS_W, ATLAS_H = 300, 112
IN_W = OUT_W = 184  # side boxes share one size
IN_H = OUT_H = 124

ATLAS_Y = 30
ROW_Y = 178
COMP_PAD_TOP = 40
COMP_X = LEFT + IN_W + GAP
COMP_Y = ROW_Y - 14
CARD1_X = COMP_X + 18
CARD1_Y = COMP_Y + COMP_PAD_TOP
CARD2_X = CARD1_X + CARD_W + GAP
CARD2_Y = CARD1_Y
COMP_W = (CARD2_X + CARD_W + 18) - COMP_X
COMP_H = COMP_PAD_TOP + CARD_H + 18
OUT_X = COMP_X + COMP_W + GAP
IN_X = LEFT
IN_Y = CARD1_Y + (CARD_H - IN_H) // 2
OUT_Y = IN_Y
ATLAS_X = CARD1_X + (CARD_W - ATLAS_W) // 2
# Canvas is derived from the content so nothing can be clipped.
W = OUT_X + OUT_W + LEFT
H = COMP_Y + COMP_H + ATLAS_Y


def esc(s: str) -> str:
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def lines(cx, y0, texts, size=12, weight='400', fill='#374151', gap=16):
    out = []
    for i, t in enumerate(texts):
        fw = f' font-weight="{weight}"' if weight != '400' else ''
        out.append(
            f'<text x="{cx:.1f}" y="{y0 + i * gap}" text-anchor="middle" font-size="{size}"{fw} fill="{fill}">{esc(t)}</text>'
        )
    return '\n'.join(out)


def card(x, y, w, h, stroke, title, body, foot):
    # Title, then body with even spacing, foot near bottom with padding.
    body_y = y + 50
    foot_y = y + h - 16
    return '\n'.join([
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="10" fill="#ffffff" stroke="{stroke}" stroke-width="1.5"/>',
        lines(x + w / 2, y + 26, [title], size=14, weight='700', fill='#1f2937', gap=0),
        lines(x + w / 2, body_y, body, size=12, fill='#374151', gap=18),
        lines(x + w / 2, foot_y, [foot], size=11, fill='#6b7280', gap=0),
    ])


def build() -> str:
    atlas_cx = ATLAS_X + ATLAS_W / 2
    card1_cx = CARD1_X + CARD_W / 2
    card1_cy = CARD1_Y + CARD_H / 2
    card2_cy = CARD2_Y + CARD_H / 2
    in_right = IN_X + IN_W
    in_cy = IN_Y + IN_H / 2
    out_left = OUT_X
    out_cy = OUT_Y + OUT_H / 2

    v_x = atlas_cx
    v_y1 = ATLAS_Y + ATLAS_H + 2
    v_y2 = CARD1_Y
    # Keep the retrieval label well right of the vertical stroke and of the Composer legend.
    label_x = v_x + 22
    label_y = (v_y1 + v_y2) / 2 + 4

    # Title sits inside the dashed frame, left of the vertical arrow (arrow x ≈ card1 centre).
    title_x = COMP_X + 16
    assert title_x + 120 < v_x - 12, (title_x, v_x)

    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="-apple-system, 'Segoe UI', Helvetica, Arial, sans-serif">
  <title>Open Meta Maps workflow</title>
  <desc>A genome-scale model is matched against the Atlas of curated pathway structures; the Layout Composer keeps a curated skeleton, completes the rest (orientation, arrangement, spacing, Bézier connections, labels), and reports layout quality.</desc>
  <defs>
    <marker id="arrow" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="#4b5563"/></marker>
  </defs>
  <rect width="{W}" height="{H}" rx="10" fill="#ffffff" stroke="#d0d7de"/>

  <!-- Atlas -->
  <rect x="{ATLAS_X}" y="{ATLAS_Y}" width="{ATLAS_W}" height="{ATLAS_H}" rx="10" fill="#e7f4ee" stroke="#3f8f6b" stroke-width="1.5"/>
  {lines(atlas_cx, ATLAS_Y + 28, ['Atlas'], size=16, weight='700', fill='#1f2937', gap=0)}
  {lines(atlas_cx, ATLAS_Y + 54, [
      'curated pathway structures',
      'Escher · WikiPathways · Reactome',
      'identity · topology · geometry · licence',
  ], size=12, fill='#374151', gap=17)}

  <!-- Layout Composer frame; title inside, left of the vertical retrieval arrow -->
  <rect x="{COMP_X}" y="{COMP_Y}" width="{COMP_W}" height="{COMP_H}" rx="12" fill="#fbf6ea" stroke="#b9852f" stroke-width="1.5" stroke-dasharray="6 4"/>
  <text x="{title_x}" y="{COMP_Y + 22}" font-size="13" font-weight="700" fill="#1f2937">Layout Composer</text>

  <!-- Input -->
  <rect x="{IN_X}" y="{IN_Y}" width="{IN_W}" height="{IN_H}" rx="10" fill="#eef3f8" stroke="#6b7f95" stroke-width="1.5"/>
  {lines(IN_X + IN_W / 2, IN_Y + 28, ['Input model'], size=16, weight='700', fill='#1f2937', gap=0)}
  {lines(IN_X + IN_W / 2, IN_Y + 56, [
      'genome-scale',
      'metabolic model',
      'SBML · COBRA JSON',
  ], size=12, fill='#374151', gap=17)}

  <!-- Uniform stage cards -->
  {card(CARD1_X, CARD1_Y, CARD_W, CARD_H, '#3f8f6b', '1  Curated skeleton', [
      'match reactions to Atlas',
      'select non-overlapping structures',
      'curated geometry kept as drawn',
  ], foot='(translated only)')}
  {card(CARD2_X, CARD2_Y, CARD_W, CARD_H, '#b9852f', '2  Algorithmic layout', [
      'completion · orientation',
      'arrangement · spacing',
      'Bézier connections · labels',
  ], foot='every GEM reaction is drawn')}

  <!-- Output -->
  <rect x="{OUT_X}" y="{OUT_Y}" width="{OUT_W}" height="{OUT_H}" rx="10" fill="#eef3f8" stroke="#6b7f95" stroke-width="1.5"/>
  {lines(OUT_X + OUT_W / 2, OUT_Y + 28, ['Complete map'], size=16, weight='700', fill='#1f2937', gap=0)}
  {lines(OUT_X + OUT_W / 2, OUT_Y + 56, [
      'Escher JSON · SVG',
      'provenance per reaction',
      'layout quality evaluation',
  ], size=12, fill='#374151', gap=17)}

  <!-- Arrows -->
  <path d="M{in_right},{in_cy:.1f} L{CARD1_X},{card1_cy:.1f}" fill="none" stroke="#4b5563" stroke-width="1.8" marker-end="url(#arrow)"/>
  <path d="M{v_x:.1f},{v_y1} L{v_x:.1f},{v_y2}" fill="none" stroke="#4b5563" stroke-width="1.8" marker-end="url(#arrow)"/>
  <text x="{label_x:.1f}" y="{label_y:.1f}" font-size="11" fill="#4b5563">structure retrieval</text>
  <path d="M{CARD1_X + CARD_W},{card2_cy:.1f} L{CARD2_X},{card2_cy:.1f}" fill="none" stroke="#4b5563" stroke-width="1.8" marker-end="url(#arrow)"/>
  <path d="M{CARD2_X + CARD_W},{card2_cy:.1f} L{out_left},{out_cy:.1f}" fill="none" stroke="#4b5563" stroke-width="1.8" marker-end="url(#arrow)"/>
</svg>
'''
    return svg


def check_render(svg_path: Path, width: int = 700):
    import cairosvg
    from PIL import Image
    png = cairosvg.svg2png(url=svg_path.resolve().as_uri(), output_width=width)
    im = Image.open(io.BytesIO(png)).convert('RGB')
    shot = Path('/tmp/opm_workflow_check.png')
    im.save(shot)
    print(f'Render check wrote {shot} ({im.size[0]}×{im.size[1]})')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--check', action='store_true', help='also write a PNG at README display width')
    args = p.parse_args()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(build(), encoding='utf-8')
    print(f'Wrote {OUT}')
    if args.check:
        check_render(OUT)


if __name__ == '__main__':
    main()
