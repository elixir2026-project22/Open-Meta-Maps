#!/usr/bin/env python3
"""Crop an Open Meta Maps Escher map.json to a preview PNG.

White background, content-tight crop, strokes and nodes scaled for a
target width (default 2800 px). Labels are kept only when their font
would stay readable at that width.

  PYTHONPATH=layout_atlas python3 tools/render_map_preview.py map.json out.png
"""
from __future__ import annotations

import argparse
import io
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'layout_atlas'))

from core import read  # noqa: E402
import line_render  # noqa: E402

# Target sizes in output pixels at --width. Nodes are capped so dense maps
# stay as networks, not solid blobs.
TARGET_STROKE = 1.2
MAX_NODE_PX = 2.2
MIN_LABEL_PX = 7.0
CROP_PAD_FRAC = 0.03
DEFAULT_WIDTH = 2800
MAX_BYTES = 3_000_000


def _content_bbox(data):
    body = data[1]
    nodes = body['nodes']
    xs, ys = [], []
    for n in nodes.values():
        xs.append(float(n['x']))
        ys.append(float(n['y']))
        if n.get('label_x') is not None:
            xs.append(float(n['label_x']))
        if n.get('label_y') is not None:
            ys.append(float(n['label_y']))
    for rx in body.get('reactions', {}).values():
        if rx.get('label_x') is not None:
            xs.append(float(rx['label_x']))
        if rx.get('label_y') is not None:
            ys.append(float(rx['label_y']))
        for seg in (rx.get('segments') or {}).values():
            for key in ('b1', 'b2'):
                pt = seg.get(key)
                if pt and 'x' in pt and 'y' in pt:
                    xs.append(float(pt['x']))
                    ys.append(float(pt['y']))
    if not xs:
        canvas = body.get('canvas') or {}
        return 0.0, 0.0, float(canvas.get('width') or 1000), float(canvas.get('height') or 1000)
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    pad = CROP_PAD_FRAC * max(xmax - xmin, ymax - ymin, 1.0)
    return xmin - pad, ymin - pad, xmax + pad, ymax + pad


def _scale_svg(svg: str, mul: float, labels: bool) -> str:
    """Multiply stroke, radius and font sizes so they hit the target pixel size."""

    def num(attr, x):
        return f'{attr}="{float(x.group(1)) * mul:.3f}"'

    svg = re.sub(r'stroke-width="([0-9.]+)"', lambda m: num('stroke-width', m), svg)
    svg = re.sub(r'\br="([0-9.]+)"', lambda m: num('r', m), svg)
    svg = re.sub(r'font-size="([0-9.]+)"', lambda m: num('font-size', m), svg)
    svg = re.sub(r'markerWidth="([0-9.]+)"', lambda m: num('markerWidth', m), svg)
    svg = re.sub(r'markerHeight="([0-9.]+)"', lambda m: num('markerHeight', m), svg)
    if not labels:
        svg = re.sub(r'<text\b[^>]*>.*?</text>', '', svg, flags=re.DOTALL)
    return svg


def _set_viewbox(svg: str, bbox, width: int, height: int) -> str:
    x0, y0, x1, y1 = bbox
    bw, bh = x1 - x0, y1 - y0
    new = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{x0:.2f} {y0:.2f} {bw:.2f} {bh:.2f}" width="{width}" height="{height}">'
    return re.sub(r'<svg\b[^>]*>', new, svg, count=1)


def render_preview(map_path: Path, out_path: Path, width: int = DEFAULT_WIDTH, max_bytes: int = MAX_BYTES) -> Path:
    from PIL import Image
    import cairosvg

    data = read(map_path)
    if not isinstance(data, list):
        import json
        data = json.loads(Path(map_path).read_text())
    svg = line_render.render(data)
    x0, y0, x1, y1 = _content_bbox(data)
    bw, bh = max(x1 - x0, 1.0), max(y1 - y0, 1.0)
    height = max(1, int(round(width * bh / bw)))
    scale = width / bw  # map units -> output pixels
    labels = (12.0 * scale) >= MIN_LABEL_PX
    # Keep map-unit sizes unless a stroke would drop below TARGET_STROKE px;
    # never let a primary node exceed MAX_NODE_PX in the output.
    mul = 1.0
    if 2.4 * scale < TARGET_STROKE:
        mul = TARGET_STROKE / (2.4 * scale)
    node_out = 5.0 * mul * scale  # line_render.PRIMARY_R
    if node_out > MAX_NODE_PX:
        mul = MAX_NODE_PX / (5.0 * scale)
    svg = _scale_svg(svg, mul, labels=labels)
    svg = _set_viewbox(svg, (x0, y0, x1, y1), width, height)
    # Force the background rect to the viewBox.
    svg = re.sub(
        r'<rect width="100%" height="100%" fill="white"/>',
        f'<rect x="{x0:.2f}" y="{y0:.2f}" width="{bw:.2f}" height="{bh:.2f}" fill="white"/>',
        svg,
        count=1,
    )
    png = cairosvg.svg2png(bytestring=svg.encode('utf-8'), output_width=width, output_height=height)
    im = Image.open(io.BytesIO(png)).convert('RGB')
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    im.save(out_path, format='PNG', optimize=True)
    if out_path.stat().st_size > max_bytes:
        for q_w in (2600, 2400, 2200, 2000):
            if q_w >= width:
                continue
            q_h = max(1, int(round(q_w * bh / bw)))
            small = im.resize((q_w, q_h), Image.Resampling.LANCZOS)
            small.save(out_path, format='PNG', optimize=True)
            if out_path.stat().st_size <= max_bytes:
                break
    return out_path


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('map', type=Path, help='Escher map.json')
    p.add_argument('out', type=Path, help='output PNG')
    p.add_argument('--width', type=int, default=DEFAULT_WIDTH)
    p.add_argument('--max-bytes', type=int, default=MAX_BYTES)
    args = p.parse_args()
    path = render_preview(args.map, args.out, width=args.width, max_bytes=args.max_bytes)
    from PIL import Image
    im = Image.open(path)
    print(f'{path} {im.size[0]}x{im.size[1]} {path.stat().st_size} bytes')


if __name__ == '__main__':
    main()
