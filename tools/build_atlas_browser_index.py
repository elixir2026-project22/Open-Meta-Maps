#!/usr/bin/env python3
"""Build the compact Atlas Browser search index (local cache).

  python3 tools/build_atlas_browser_index.py
  python3 tools/build_atlas_browser_index.py --page

Writes data/cache/atlas_browser/index.json.gz (gitignored). The tool server also
builds this on first use and rebuilds when data/atlas_registry.json changes.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'layout_atlas'))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, default=None, help='Output path (default: data/cache/atlas_browser/index.json.gz)')
    p.add_argument('--page', action='store_true', help='Also regenerate the Atlas Browser HTML page')
    args = p.parse_args(argv)
    import atlas_browser
    summary = atlas_browser.build_index(args.out)
    if args.page:
        path = atlas_browser.page_html()
        summary['page'] = str(path.relative_to(ROOT))
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
