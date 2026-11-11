"""Unit tests for final-geometry layout quality (python3 test_layout_quality.py)."""
import layout_quality as lq


def _rx(bigg, mid_id, met_id, met_bigg, label=(0, 0)):
    return {
        "bigg_id": bigg,
        "name": bigg,
        "label_x": label[0],
        "label_y": label[1],
        "metabolites": [{"bigg_id": met_bigg, "coefficient": -1}],
        "segments": {"0": {"from_node_id": mid_id, "to_node_id": met_id}},
    }


def _crossing_map():
    """Two algorithm×algorithm crossings: one near (150,150), one near (750,750)."""
    nodes = {
        "m1": {"node_type": "midmarker", "x": 100.0, "y": 100.0},
        "a1": {"node_type": "metabolite", "x": 200.0, "y": 200.0, "bigg_id": "A_c", "node_is_primary": True,
               "label_x": 210.0, "label_y": 210.0},
        "m2": {"node_type": "midmarker", "x": 100.0, "y": 200.0},
        "a2": {"node_type": "metabolite", "x": 200.0, "y": 100.0, "bigg_id": "B_c", "node_is_primary": True,
               "label_x": 210.0, "label_y": 90.0},
        "m3": {"node_type": "midmarker", "x": 700.0, "y": 700.0},
        "a3": {"node_type": "metabolite", "x": 800.0, "y": 800.0, "bigg_id": "C_c", "node_is_primary": True,
               "label_x": 810.0, "label_y": 810.0},
        "m4": {"node_type": "midmarker", "x": 700.0, "y": 800.0},
        "a4": {"node_type": "metabolite", "x": 800.0, "y": 700.0, "bigg_id": "D_c", "node_is_primary": True,
               "label_x": 810.0, "label_y": 690.0},
        # Anchor the window origin so region boxes are deterministic.
        "pad": {"node_type": "midmarker", "x": 0.0, "y": 0.0},
    }
    reactions = {
        "1": _rx("R1", "m1", "a1", "A_c", (100, 90)),
        "2": _rx("R2", "m2", "a2", "B_c", (100, 210)),
        "3": _rx("R3", "m3", "a3", "C_c", (700, 690)),
        "4": _rx("R4", "m4", "a4", "D_c", (700, 810)),
    }
    return [
        {"map_name": "toy-cross", "map_id": "toy-cross", "schema": "1.0.0", "metadata": {}},
        {"nodes": nodes, "reactions": reactions, "canvas": {"x": 0, "y": 0, "width": 900, "height": 900}},
    ]


def test_region_window_keys_membership():
    # Event at x0+250 belongs to windows starting at x0+50 (k=floor(250/200)-1=0) and x0+200 (k=1)?
    # floor(250/200)=1 → k in {0,1}; windows [0,400) and [200,600).
    keys = lq._region_window_keys(250.0, 250.0, 0.0, 0.0)
    assert set(keys) == {(0, 0), (0, 1), (1, 0), (1, 1)}
    for kx, ky in keys:
        assert 0 + kx * 200 <= 250 < 0 + kx * 200 + 400
        assert 0 + ky * 200 <= 250 < 0 + ky * 200 + 400


def test_worst_region_count_matches_drawn_box():
    data = _crossing_map()
    summary, issues = lq.evaluate(data, max_issues=10)
    crossings = [c for c in issues["crossings"] if c["kind"] in ("algorithm_algorithm", "algorithm_expert")]
    assert len(crossings) == 2, crossings
    # Uncapped list used for windows even when max_issues is tiny.
    assert summary["core"]["crossings"]["algorithm_algorithm"] == 2

    near = [c for c in crossings if abs(c["x"] - 150) < 1 and abs(c["y"] - 150) < 1]
    far = [c for c in crossings if abs(c["x"] - 750) < 1 and abs(c["y"] - 750) < 1]
    assert len(near) == 1 and len(far) == 1

    windows = summary["worst_regions"]
    assert windows, "expected at least one worst region"
    for w in windows:
        inside = sum(
            1 for c in crossings
            if w["x"] <= c["x"] < w["x"] + w["w"] and w["y"] <= c["y"] < w["y"] + w["h"]
        )
        assert w["crossings"] == inside, (w, inside, crossings)
        # Reactions listed must correspond to crossings inside the drawn box.
        for rid in w["reactions"]:
            assert any(
                rid in c["reactions"] and w["x"] <= c["x"] < w["x"] + w["w"] and w["y"] <= c["y"] < w["y"] + w["h"]
                for c in crossings
            ), (rid, w, crossings)

    # The near and far crossings must not share a single 400×400 window count of 2.
    assert all(w["crossings"] <= 1 for w in windows), windows
    assert max(w["crossings"] for w in windows) == 1


def test_window_uses_uncapped_crossings():
    data = _crossing_map()
    summary, _ = lq.evaluate(data, max_issues=0)
    assert summary["core"]["crossings"]["algorithm_algorithm"] == 2
    assert summary["details"]["crossing_density"]["worst_window_crossings"] == 1
    assert any(w["crossings"] == 1 for w in summary["worst_regions"])


if __name__ == "__main__":
    test_region_window_keys_membership()
    test_worst_region_count_matches_drawn_box()
    test_window_uses_uncapped_crossings()
    print("ok")
