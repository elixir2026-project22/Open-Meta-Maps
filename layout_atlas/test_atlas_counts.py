import importlib.util
from pathlib import Path

from atlas_browser import reaction_set

_spec = importlib.util.spec_from_file_location(
    'atlas_module_figures', Path(__file__).resolve().parents[1] / 'tools' / 'atlas_module_figures.py')
figures = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(figures)


def _semantic(mid, source, map_id, rxns):
    return dict(id=mid, source=source, map_id=map_id, quality='semantic_transfer', size=len(rxns),
                reactions=[dict(id=r, references={'bigg.reaction': [r]}) for r in rxns], nodes=[], edges=[])


def _geometry(mid, map_id, keys, fingerprint):
    return dict(id=mid, source='wikipathways', map_id=map_id, quality='geometry_only', size=len(keys),
                reactions=[], nodes=[dict(key=k) for k in keys], edges=[], geometry_fingerprint=fingerprint)


def test_distinct_structures_drop_nested_small_and_cross_source_duplicates():
    mods = [
        _semantic('a1', 'escher', 'm1', ['PGI']),
        _semantic('a2', 'escher', 'm1', ['PGI', 'PFK', 'FBA', 'TPI']),
        _semantic('a3', 'escher', 'm1', ['PGI', 'PFK', 'FBA', 'TPI', 'GAPD', 'PGK', 'PGM', 'ENO']),
        _semantic('a4', 'escher', 'm2', ['CS', 'ACONTa', 'ACONTb']),
        _semantic('a5', 'escher', 'm2', ['SUCDi', 'FUM']),
        _semantic('b1', 'human_maps', 'h1', ['CS', 'ACONTa', 'ACONTb']),
        _semantic('b2', 'human_maps', 'h1', ['CS', 'ACONTa']),
        _geometry('g1', 'w1', ['n1', 'n2', 'n3'], 'shapeA'),
        _geometry('g2', 'w1', ['n1', 'n2', 'n3', 'n4', 'n5', 'n6'], 'shapeB'),
        _geometry('g3', 'w2', ['p1', 'p2', 'p3', 'p4', 'p5', 'p6'], 'shapeB'),
    ]
    summary, _, _ = figures.compute_overlap(mods)
    by = summary['by_source']
    assert summary['total_source_maps'] == 5
    assert summary['total_modules'] == 10
    assert by['escher']['source_maps'] == 2
    assert by['escher']['distinct_structures'] == 2
    assert by['human_maps']['distinct_structures'] == 1
    assert by['wikipathways_geometry']['distinct_structures'] == 1
    assert summary['distinct_structures_global'] == 3


def test_maximal_windows_compare_within_one_map_only():
    mods = [
        _semantic('a1', 'escher', 'm1', ['PGI', 'PFK', 'FBA']),
        _semantic('a2', 'escher', 'm2', ['PGI', 'PFK', 'FBA', 'TPI']),
    ]
    sets = {m['id']: reaction_set(m) for m in mods}
    assert {m['id'] for m in figures.maximal_windows(mods, sets)} == {'a1', 'a2'}
