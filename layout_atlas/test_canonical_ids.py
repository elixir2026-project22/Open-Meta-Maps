"""Canonical reaction-id mapping checks (python3 test_canonical_ids.py)."""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from canonical_ids import (
    CANONICAL_NS,
    HUMAN_GEM_JSON,
    measure_canonical_share,
    _is_mar,
)
from sbgn import _normalize_ident


def _load(rel):
    with gzip.open(Path(__file__).resolve().parent.parent / rel, 'rt') as f:
        return json.load(f)


def main():
    assert _normalize_ident('https://identifiers.org/CHEBI:58165') == 'chebi:58165'
    assert _normalize_ident('http://identifiers.org/chebi/58165') == 'chebi:58165'

    assert HUMAN_GEM_JSON.exists(), HUMAN_GEM_JSON
    cross = json.loads(HUMAN_GEM_JSON.read_text())
    assert cross.get('license') == 'CC-BY-4.0'
    assert len(cross.get('reactions') or {}) >= 1000

    human = _load('data/extensions/human_maps/modules.json.gz')
    # MAR must not be the sole bigg.reaction identity after remap
    mar_as_bigg = 0
    native = 0
    for m in human:
        for r in m.get('reactions') or []:
            refs = r.get('references') or {}
            if refs.get('human.gem'):
                native += 1
            for v in refs.get('bigg.reaction') or []:
                if _is_mar(v):
                    mar_as_bigg += 1
    assert native > 0, 'human.gem native namespace missing'
    assert mar_as_bigg == 0, f'MAR still stored under bigg.reaction ({mar_as_bigg})'
    share = measure_canonical_share(human)
    assert share['share'] >= 0.80, share

    sbgn = _load('data/extensions/metabolismregulation/modules.json.gz')
    share_s = measure_canonical_share(sbgn)
    assert share_s['share'] >= 0.20, share_s
    assert share_s['by_namespace'].get('rhea', 0) >= 100

    panther = _load('data/extensions/panther/modules.json.gz')
    # Native ids retained; canonical share may remain near zero without inventing ids
    for m in panther[:5]:
        for r in m.get('reactions') or []:
            assert 'celldesigner.reaction' in (r.get('references') or {})

    print(
        f"PASS canonical ids: human_maps share={share['share']:.3f} "
        f"metabolismregulation share={share_s['share']:.3f} "
        f"namespaces={list(CANONICAL_NS)}"
    )


if __name__ == '__main__':
    main()
