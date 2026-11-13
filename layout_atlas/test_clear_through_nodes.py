"""Gated clear-through-nodes restyle (e_coli_core example_results)."""
from pathlib import Path

from evaluation.escher_bridge import load_escher
import edge_style
import layout_quality as lq

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "example_results" / "e_coli_core" / "OpenMetaMaps_e_coli_core_result.json"


def test_curves_through_nodes_count_matches_evaluate():
    data = load_escher(EXAMPLE)
    skip = (data[0].get("metadata") or {}).get("preserved_edges") or []
    n = lq.curves_through_nodes_count(data, skip)
    summary, _issues = lq.evaluate(data, skip)
    assert n == summary["core"]["overlaps"]["curves_through_nodes"]
    assert n >= 0


def test_clear_through_nodes_no_regression_on_e_coli_core():
    data = load_escher(EXAMPLE)
    skip = (data[0].get("metadata") or {}).get("preserved_edges") or []
    before = lq.curves_through_nodes_count(data, skip)
    before_cross = lq.CrossingIndex(data, skip).headline()
    report = edge_style.clear_through_nodes(data, skip)
    after = lq.curves_through_nodes_count(data, skip)
    after_cross = lq.CrossingIndex(data, skip).headline()
    assert report["before"] == before
    assert after == report["after"]
    assert after_cross == report["after_crossings"]
    if report["accepted"]:
        assert after < before
        assert after_cross <= before_cross
    else:
        # Reverted: live map unchanged.
        assert after == before
        assert after_cross == before_cross
