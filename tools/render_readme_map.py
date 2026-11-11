#!/usr/bin/env python3
"""Rasterise a composed map.svg to a compact PNG for the README example.

Example (default curated Bézier layout):

  python3 project.py compose --gem examples/e_coli_core.xml --out /tmp/e_coli_core_readme/ --source-paths
  python3 tools/render_readme_map.py /tmp/e_coli_core_readme/map.svg docs/images/e_coli_core_map.png

Needs the ``cairosvg`` package (``python3 -m pip install cairosvg``). Output is ~1200 px wide.
"""
from __future__ import annotations

import argparse
import io
import sys
from pathlib import Path


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('svg', type=Path, help='input map.svg from compose')
    p.add_argument('png', type=Path, help='output PNG path')
    p.add_argument('--width', type=int, default=1200, help='target width in CSS pixels (default 1200)')
    args = p.parse_args(argv)
    try:
        import cairosvg
        from PIL import Image
    except ImportError as e:
        raise SystemExit('Need cairosvg and Pillow: python3 -m pip install cairosvg Pillow') from e
    svg = args.svg.resolve()
    if not svg.is_file():
        raise SystemExit(f'SVG not found: {svg}')
    out = args.png.resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    w = max(200, int(args.width))
    png = cairosvg.svg2png(url=svg.as_uri(), output_width=w)
    im = Image.open(io.BytesIO(png)).convert('RGB')
    im.save(out, format='PNG', optimize=True)
    kb = out.stat().st_size / 1024
    print(f'Wrote {out} ({im.size[0]}×{im.size[1]}, {kb:.0f} KiB)')
    if kb > 500:
        print('warning: PNG is larger than ~500 KiB; consider a narrower --width', file=sys.stderr)


if __name__ == '__main__':
    main()
