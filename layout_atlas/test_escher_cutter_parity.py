"""Parity checks for the in-repo community Escher cutter.

Compares regenerated windows to ``module_observations.json.gz`` on:
  - reaction-set overlap (reported; exact equality is not required)
  - per-map reaction coverage (must stay 100% on both stored and generated)
  - size bounds (≤40 reactions, ≤100 nodes)

Run: ``python3 test_escher_cutter_parity.py``
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from escher_cutter import MAX_NODES, MAX_REACTIONS, cut_all, parity_report


def main():
    report = parity_report()
    summary = report['summary']
    print(json.dumps(summary, indent=2, sort_keys=True))
    assert summary['maps'] == 22, summary['maps']
    assert summary['all_size_bounds_ok'], 'generated modules exceed size caps'
    assert summary['min_generated_coverage'] >= 0.999, (
        f"generated reaction coverage dropped: {summary['min_generated_coverage']}"
    )
    assert summary['min_stored_coverage'] >= 0.999, (
        f"stored reaction coverage unexpected: {summary['min_stored_coverage']}"
    )
    # Soft set parity against size-capped stored observations. Exact equality is
    # not required: the historical cutter is absent, and retrieval gates showed
    # the in-repo cut improves coverage on e_coli_core / iIT341.
    assert summary['recall'] >= 0.65, f"set recall too low: {summary['recall']}"
    cut = cut_all()
    assert len(cut['modules']) >= 5000
    for m in cut['modules']:
        assert 1 <= m['size'] <= MAX_REACTIONS
        assert len(m['nodes']) <= MAX_NODES
        assert m.get('extraction', {}).get('method') == 'escher_community_windows'
        assert 'parameters' in m['extraction'] and 'counts' in m['extraction']
        for r in m['reactions']:
            assert 'bigg.reaction' in (r.get('references') or {})
    print(
        f"PASS escher cutter parity: recall={summary['recall']:.3f} "
        f"generated={summary['generated_modules']} stored={summary['stored_modules']} "
        f"coverage={summary['min_generated_coverage']:.3f}"
    )


if __name__ == '__main__':
    main()
