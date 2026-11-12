#!/usr/bin/env python3
"""Build the compact Atlas Browser search index (reproducible).

  python3 tools/build_atlas_browser_index.py
  python3 tools/build_atlas_browser_index.py --out data/atlas_browser_index.json.gz

Reads open_core.sqlite + default extensions and writes a gzipped JSON index used by
`layout_atlas/atlas_browser.py` and the tool server `/api/atlas/*` endpoints.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'layout_atlas'))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, default=None, help='Output path (default: data/atlas_browser_index.json.gz)')
    p.add_argument('--page', action='store_true', help='Also regenerate the Atlas Browser HTML page')
    args = p.parse_args(argv)
    import atlas_browser
    summary = atlas_browser.build_index(args.out)
    if args.page:
        path = atlas_browser.page_html()
        summary['page'] = str(path.relative_to(ROOT))
    print(json_dumps(summary))
    return 0


def json_dumps(obj):
    import json
    return json.dumps(obj, ensure_ascii=False, indent=2)


if __name__ == '__main__':
    raise SystemExit(main())
