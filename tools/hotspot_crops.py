#!/usr/bin/env python3
"""Render before/after crops of the top crossing hotspots from two map.json files."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'layout_atlas'))

from core import read  # noqa: E402
import layout_quality as lq  # noqa: E402


def _draw(ax, data, win, title):
    body = data[1]
    nodes, reactions = body['nodes'], body['reactions']
    x0, y0, w, h = win['x'], win['y'], win['w'], win['h']
    pad = 20
    ax.set_xlim(x0 - pad, x0 + w + pad)
    ax.set_ylim(y0 + h + pad, y0 - pad)  # SVG-style y-down
    ax.set_aspect('equal')
    ax.set_title(title, fontsize=10)
    ax.add_patch(plt.Rectangle((x0, y0), w, h, fill=False, ls='--', lw=1.2, ec='#c0392b'))
    # edges
    skip = lq._skip_set if False else None
    meta = data[0].get('metadata') or {}
    skip_edges = meta.get('preserved_edges') or []
    edges, _ = lq.build_edges(data, skip_edges)
    for e in edges:
        pts = e['pts']
        if not len(pts):
            continue
        # draw if any point near window
        if not ((pts[:, 0] >= x0 - pad) & (pts[:, 0] <= x0 + w + pad) &
                (pts[:, 1] >= y0 - pad) & (pts[:, 1] <= y0 + h + pad)).any():
            continue
        color = '#1e6b3a' if not e['expert'] else '#7a7a7a'
        ax.plot(pts[:, 0], pts[:, 1], color=color, lw=0.8, alpha=0.85)
    # nodes
    for n in nodes.values():
        t = n.get('node_type')
        if t not in ('metabolite', 'midmarker'):
            continue
        if not (x0 - pad <= n['x'] <= x0 + w + pad and y0 - pad <= n['y'] <= y0 + h + pad):
            continue
        if t == 'midmarker':
            ax.plot(n['x'], n['y'], 'o', ms=2.5, color='#222')
        else:
            ax.plot(n['x'], n['y'], 'o', ms=4, color='#d0d0d0', markeredgecolor='#444', markeredgewidth=0.4)
            bigg = n.get('bigg_id') or ''
            if bigg:
                ax.text(n['x'] + 4, n['y'] - 4, bigg, fontsize=5, color='#333')
    for rx in reactions.values():
        mid = next((s['from_node_id'] for s in rx['segments'].values()
                    if nodes.get(s['from_node_id'], {}).get('node_type') == 'midmarker'), None)
        if mid is None:
            continue
        m = nodes[mid]
        if x0 - pad <= m['x'] <= x0 + w + pad and y0 - pad <= m['y'] <= y0 + h + pad:
            ax.text(m['x'] + 6, m['y'] + 6, rx['bigg_id'], fontsize=6, fontweight='bold', color='#111')
    ax.set_xticks([]); ax.set_yticks([])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--before', required=True)
    ap.add_argument('--after', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--k', type=int, default=3)
    ap.add_argument('--prefix', default='hotspot')
    args = ap.parse_args()
    before = read(args.before)
    after = read(args.after)
    meta = before[0].get('metadata') or {}
    skip = meta.get('preserved_edges') or []
    summary, _ = lq.evaluate(before, skip)
    windows = summary.get('worst_regions') or []
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for i, win in enumerate(windows[:args.k]):
        fig, axes = plt.subplots(1, 2, figsize=(10, 5), dpi=140)
        _draw(axes[0], before, win, f"before  crossings={win['crossings']}")
        # recount crossings in the same box on the after map (shared detector events)
        _total, _cross, events = lq.algorithm_crossings(after, skip)
        after_n = sum(1 for e in events if e['kind'] in ('algorithm_algorithm', 'algorithm_expert')
                      and win['x'] <= e['x'] < win['x'] + win['w'] and win['y'] <= e['y'] < win['y'] + win['h'])
        _draw(axes[1], after, win, f"after  crossings≈{after_n}")
        fig.suptitle(f"Top hotspot #{i+1} @ ({win['x']:.0f},{win['y']:.0f})", fontsize=11)
        fig.tight_layout()
        dest = out / f'{args.prefix}_{i+1}.png'
        fig.savefig(dest, bbox_inches='tight')
        plt.close(fig)
        rows.append(dict(index=i + 1, path=str(dest), before_crossings=win['crossings'], after_crossings=after_n, region=win))
        print('wrote', dest)
    (out / f'{args.prefix}_summary.json').write_text(__import__('json').dumps(rows, indent=2) + '\n')


if __name__ == '__main__':
    main()
