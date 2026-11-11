"""Smoke tests for the unified layout evaluation path."""
from evaluation.pathway import PRODUCT, REACTANT, REACTION, SPECIES, Pathway
from evaluation import metrics
import networkx as nx


def _toy():
    g = nx.DiGraph()
    g.add_node("A", kind=SPECIES, label="A", compartment="c", width=60, height=24)
    g.add_node("B", kind=SPECIES, label="B", compartment="c", width=60, height=24)
    g.add_node("C", kind=SPECIES, label="C", compartment="e", width=60, height=24)
    g.add_node("r1", kind=REACTION, label="r1", compartment=None, width=10, height=10)
    g.add_edge("A", "r1", role=REACTANT)
    g.add_edge("r1", "B", role=PRODUCT)
    g.add_edge("r1", "C", role=PRODUCT)
    return Pathway(name="toy", graph=g)


def test_participant_consistency_opposite():
    p = _toy()
    pos = {"r1": (0.0, 0.0), "A": (-10.0, 0.0), "B": (10.0, 1.0), "C": (10.0, -1.0)}
    assert metrics.participant_consistency(p, pos) == 1.0


def test_participant_consistency_same_side():
    p = _toy()
    pos = {"r1": (0.0, 0.0), "A": (-10.0, 0.0), "B": (-10.0, 2.0), "C": (-8.0, -2.0)}
    assert metrics.participant_consistency(p, pos) == 0.0


def test_long_edge_share():
    p = _toy()
    pos = {"r1": (0.0, 0.0), "A": (-1.0, 0.0), "B": (1.0, 0.0), "C": (100.0, 0.0)}
    m = metrics.long_edge_metrics(p, pos, long_factor=3.0)
    assert m["long_edge_share"] is not None
    assert 0.0 < m["long_edge_share"] <= 1.0


def test_overlaps_no_touch():
    if not metrics.shapely_available():
        return
    p = _toy()
    pos = {"r1": (0.0, 0.0), "A": (-100.0, 0.0), "B": (100.0, 0.0), "C": (0.0, 100.0)}
    m = metrics.overlaps(p, pos)
    assert m["overlapping_pairs"] == 0.0


if __name__ == "__main__":
    test_participant_consistency_opposite()
    test_participant_consistency_same_side()
    test_long_edge_share()
    test_overlaps_no_touch()
    print("ok")
