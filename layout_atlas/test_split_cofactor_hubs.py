"""Gated hub split on committed example_results maps (e_coli_core / iIT341)."""
import json
from pathlib import Path

from cofactor_hubs import cofactor_hubs
from core import load_model
from evaluation.escher_bridge import load_escher
from hub_drawing import pathway_from_model
from layout_gate import score_map
from layout_quality_evaluator import HUB_DIRS, METRIC_EPS, metrics_not_worse
from split_cofactor_hubs import apply_gated

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"
EXAMPLE_RESULTS = ROOT / "example_results"


def _curated_reactions(name: str) -> set[str]:
    """Reactions whose layout comes from the curated atlas (not de-novo).

    Shipped example_results still label atlas reuse as expert_*; treat those the
    same as curated_* provenance strings.
    """
    rows = json.loads((EXAMPLE_RESULTS / name / "reaction_provenance.json").read_text())
    curated = set()
    for row in rows:
        for key in ("layout_provenance", "provenance"):
            value = str(row.get(key) or "")
            if value.startswith("curated") or value.startswith("expert"):
                curated.add(row["reaction"])
                break
    return curated


def test_gate_columns_cover_directed_metrics():
    assert "angular_resolution" in HUB_DIRS
    assert "crossing_angular_resolution" in HUB_DIRS
    assert "aspect_ratio" in HUB_DIRS
    assert "homogeneity" in HUB_DIRS
    assert "species_drawn" not in HUB_DIRS
    assert "reactions_drawn" not in HUB_DIRS
    assert len(HUB_DIRS) == 11


def test_angle_candidates_include_fallback_rotations():
    from split_cofactor_hubs import _FALLBACK_ROTATIONS, _angle_candidates

    # Empty movers → still expose fallbacks so stubborn hubs can retry small angles.
    angles = _angle_candidates([{}, {"nodes": {}, "reactions": {}}], movers=[], glyph="1", rid_to_key={})
    for ang in _FALLBACK_ROTATIONS:
        assert ang in angles


def test_split_e_coli_core_no_regression():
    raw = load_escher(EXAMPLE_RESULTS / "e_coli_core" / "OpenMetaMaps_e_coli_core_result.json")
    model = load_model(EXAMPLES / "e_coli_core.xml")
    pathway = pathway_from_model(model)
    before = score_map(pathway, raw)
    hubs0 = cofactor_hubs(raw, model)["hub_count"]
    assert hubs0 == 6
    out, report = apply_gated(
        raw, model, pathway=pathway, curated=_curated_reactions("e_coli_core")
    )
    assert report["accepted"]
    assert report["hubs_before"] == 6
    # Fallback angles + multi-pass clear more exchange hubs toward cofactors-as-clones.
    assert report["hubs_after"] == 2
    assert cofactor_hubs(out, model)["hub_count"] == 2
    assert report["hubs_after"] < report["hubs_before"]
    after = score_map(pathway, out)
    assert metrics_not_worse(before, after, HUB_DIRS, METRIC_EPS)
    clones = (out[0].get("metadata") or {}).get("hub_split_clones") or []
    assert clones
    # Pin reader keeps the original hub glyph for split metabolites.
    from curated_axes import _absolute_from_escher

    before_xy = _absolute_from_escher(raw[1], "t")["metabolites"]
    after_xy = _absolute_from_escher(out[1], "t", skip_nids=clones)["metabolites"]
    for mid in before_xy:
        if mid in after_xy:
            assert after_xy[mid] == before_xy[mid]


def test_split_iIT341_no_regression():
    raw = load_escher(EXAMPLE_RESULTS / "iIT341" / "OpenMetaMaps_iIT341_result.json")
    model = load_model(EXAMPLES / "iIT341.xml")
    pathway = pathway_from_model(model)
    before = score_map(pathway, raw)
    hubs0 = cofactor_hubs(raw, model)["hub_count"]
    assert hubs0 == 9
    out, report = apply_gated(
        raw, model, pathway=pathway, curated=_curated_reactions("iIT341")
    )
    assert report["accepted"]
    assert report["hubs_before"] == 9
    assert report["hubs_after"] == 5
    assert cofactor_hubs(out, model)["hub_count"] == 5
    assert report["hubs_after"] < report["hubs_before"]
    after = score_map(pathway, out)
    assert metrics_not_worse(before, after, HUB_DIRS, METRIC_EPS)
