#!/usr/bin/env python3
"""Download optional evaluation models listed in manifest.json into downloaded/.

Usage:
  python3 examples/eval_models/fetch_eval_models.py
  python3 examples/eval_models/fetch_eval_models.py --id e_coli_iCH360
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "manifest.json"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch_one(entry: dict, dest_dir: Path, force: bool = False) -> Path:
    dest = dest_dir / entry["path"]
    dest.parent.mkdir(parents=True, exist_ok=True)
    expected = entry["sha256"].lower()
    if dest.exists() and not force:
        got = _sha256(dest)
        if got == expected:
            print(f"ok  {entry['path']} (cached)")
            return dest
        print(f"re-download {entry['path']}: sha256 mismatch ({got[:12]}…)", file=sys.stderr)
    print(f"get {entry['url']}")
    req = urllib.request.Request(entry["url"], headers={"User-Agent": "OpenMetaMaps-eval-fetch/1.0"})
    with urllib.request.urlopen(req, timeout=180) as r:
        data = r.read()
    got = hashlib.sha256(data).hexdigest()
    if got != expected:
        raise SystemExit(f"sha256 mismatch for {entry['path']}: got {got}, expected {expected}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    tmp.write_bytes(data)
    tmp.replace(dest)
    print(f"ok  {entry['path']}  licence={entry.get('licence')}  {got[:16]}…")
    return dest


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Fetch licensed evaluation models on demand.")
    p.add_argument("--manifest", type=Path, default=MANIFEST)
    p.add_argument("--id", action="append", help="Fetch only this manifest id (repeatable)")
    p.add_argument("--force", action="store_true", help="Re-download even when sha256 matches")
    p.add_argument("--list", action="store_true", help="Print manifest entries and exit")
    args = p.parse_args(argv)
    man = json.loads(args.manifest.read_text())
    files = list(man.get("files") or [])
    if args.list:
        for e in files:
            print(f"{e['id']}\t{e.get('licence')}\t{e['path']}\t{e.get('source')}")
        for s in man.get("skipped") or []:
            print(f"SKIP\t{s['id']}\t{s.get('reason')}")
        return 0
    if args.id:
        want = set(args.id)
        files = [e for e in files if e["id"] in want]
        missing = want - {e["id"] for e in files}
        if missing:
            raise SystemExit(f"unknown id(s): {', '.join(sorted(missing))}")
    dest_dir = ROOT / (man.get("download_dir") or "downloaded")
    for e in files:
        if not e.get("licence"):
            print(f"skip {e.get('id')}: no clear licence", file=sys.stderr)
            continue
        fetch_one(e, dest_dir, force=args.force)
    print(f"downloaded → {dest_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
