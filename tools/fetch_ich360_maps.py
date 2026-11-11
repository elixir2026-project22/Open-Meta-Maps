#!/usr/bin/env python3
"""Fetch pinned iCH360 custom Escher maps into data/source_inputs/raw/ich360/.

Source: https://github.com/marco-corrao/iCH360 (Visualisation/)
Pinned commit is written to SOURCE_COMMIT.txt; each file is SHA-256 checked
against the expected digests below (update digests when re-pinning).

Excluded by design: model_map/compressed_map.json (CP pseudo-reactions; not
BiGG-aligned; upstream readme says it cannot be used directly with the model).
"""
from __future__ import annotations

import hashlib
import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data/source_inputs/raw/ich360'
REPO = 'marco-corrao/iCH360'
# Tip of main when this pin was recorded (2026-11-11).
COMMIT = '08dd4388096dcf0c88f36d7229c0b06688582744'

# Maps admitted to the Atlas (10). compressed_map.json deliberately omitted.
FILES = {
    'model_map/full_map.json': '071232f8c77e0842b4ef27f40e3f607502558debbe7e38add79a14bd74b7945d',
    'pathway_maps/aminoacids_biosynthesis.json': 'd33812270f230d3519da63595ccdd423ec5a8ef04dfbdfec68f4ae1ecd715e56',
    'pathway_maps/c1_pool.json': '24112984da7a78944da8daf3f54fb0062ccafce05ae1f0ed0a1ab5be39688f59',
    'pathway_maps/fatty_acids_biosynthesis.json': '28babdb8351367f7d637e1da9d54bf4fe56c86f752abc43e85c4fe4670034052',
    'pathway_maps/nucleotides_biosynthesis.json': '6673791af9a5154282bd59b5b8239cae360285e0f7a218fff43e4078b68f1a64',
    'pathways_not_in_model/cofactors_de_novo_biosynthesis.json': 'c8788438fd930f423872bdff5f6bd63ae843c47afe1b409ddab8c411a766540b',
    'pathways_not_in_model/ions_metals_uptake.json': '90462aadb6f89385df5dd923546f476ac4d052bc55260ca400be59d4e7f2802b',
    'pathways_not_in_model/kdo_biosynthesis.json': '9419669cc7e51774d22d391fc3ff540d8472309a9139f91c1477196e1441dca0',
    'pathways_not_in_model/murein_biosynthesis.json': 'd11444b3cd53a7667621e83c6f6caf4736d4c7c45136410f2915d86b1c0e2145',
    'pathways_not_in_model/phospholipids_biosynthesis.json': 'cdb6fe8851aa5d006c065e4116db3df24e65b758617e32164ba48eb97696e252',
}

EXCLUDED = {
    'model_map/compressed_map.json': (
        'CP pseudo-reactions; upstream Visualisation/readme.md says the compressed '
        'map cannot be used directly with iCH360'
    ),
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def fetch(rel: str) -> bytes:
    url = f'https://raw.githubusercontent.com/{REPO}/{COMMIT}/Visualisation/{rel}'
    with urllib.request.urlopen(url, timeout=120) as resp:
        return resp.read()


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    dry = '--dry-run' in argv
    RAW.mkdir(parents=True, exist_ok=True)
    records = []
    for rel, expect in FILES.items():
        dest = RAW / rel
        if dry and dest.exists():
            raw = dest.read_bytes()
        else:
            raw = fetch(rel)
            if not dry:
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(raw)
        digest = sha256(raw)
        if digest != expect:
            raise SystemExit(f'checksum mismatch for {rel}: got {digest}, expected {expect}')
        records.append(dict(
            path=rel, sha256=digest, bytes=len(raw),
            url=f'https://raw.githubusercontent.com/{REPO}/{COMMIT}/Visualisation/{rel}',
        ))
        print(f'OK {rel} ({len(raw)} bytes)')
    if not dry:
        (RAW / 'SOURCE_COMMIT.txt').write_text(COMMIT + '\n')
        (RAW / 'FILES.json').write_text(json.dumps(dict(
            repository=f'https://github.com/{REPO}',
            commit=COMMIT,
            access_date='2026-11-11',
            included=records,
            excluded=[{'path': k, 'reason': v} for k, v in EXCLUDED.items()],
        ), indent=2) + '\n')
        # Keep Visualisation/readme.md for audit (no restrictive per-file notes).
        (RAW / 'readme.md').write_bytes(fetch('readme.md'))
    print(f'pinned {COMMIT}; {len(FILES)} maps under {RAW.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
