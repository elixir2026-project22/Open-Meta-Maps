"""Stage B cutter helpers: size-scaled roots and independent-window dedupe."""
from database import scaled_max_roots, independent_reaction_sets
from retrieval import (
    select,
    GEOMETRY_CRUMB_DEMOTION,
    COMPLETE_EQUATION_BONUS,
    CE_CRUMB_LOCK_MIN_COMPETING,
    _free_standing_ce_crumb,
    _ce_crumb_overlaps_large_non_ce,
    _competing_free_standing_ce_crumbs,
    _is_complete_equation,
)


def test_scaled_max_roots():
    assert scaled_max_roots(10) == 10
    assert scaled_max_roots(48) == 48
    assert scaled_max_roots(49) == 24
    assert scaled_max_roots(100) == 24
    assert scaled_max_roots(200) == 24
    assert scaled_max_roots(0) == 24


def test_independent_reaction_sets():
    big12 = tuple(f'r{i}' for i in range(12))
    small8 = tuple(f'r{i}' for i in range(8))
    kept = independent_reaction_sets([
        ('a', 'b'), ('a', 'b', 'c', 'd'), small8, big12,
        ('x', 'y'), ('x', 'y') + tuple(f'z{i}' for i in range(10)),
    ])
    kept_sets = [set(s) for s in kept]
    assert {'a', 'b', 'c', 'd'} in kept_sets
    assert {'a', 'b'} not in kept_sets
    assert set(big12) in kept_sets
    assert set(small8) not in kept_sets
    assert {'x', 'y'} in kept_sets


def _feat(biol, compact=1.0):
    return [0.0, biol, 0.0, 0.0, 0.0, compact]


def test_select_demotes_geometry_crumbs():
    crumb = dict(id='crumb', source='src_a', reaction_ids=['A', 'B'], score=0.5,
                 features=_feat(0.0), layout_origin={})
    dense = dict(id='dense', source='src_b', reaction_ids=['A', 'B', 'C', 'D'], score=0.4,
                 features=_feat(2.0), layout_origin={},
                 match_level='complete_equation', complete_chemical_equation_verified=True)
    assert [c['id'] for c in select([crumb, dense])] == ['dense']
    alone = select([crumb])
    assert alone and alone[0]['selection_score'] == (0.5 - GEOMETRY_CRUMB_DEMOTION)


def test_complete_equation_bonus_constant():
    assert COMPLETE_EQUATION_BONUS == 4.0
    assert COMPLETE_EQUATION_BONUS > 2.2
    assert CE_CRUMB_LOCK_MIN_COMPETING == 8
    assert CE_CRUMB_LOCK_MIN_COMPETING > 3


def test_free_standing_ce_crumb_helper():
    ce_big = dict(id='ce_big', reaction_ids=['A0', 'A1', 'A2', 'A3'],
                  match_level='complete_equation', complete_chemical_equation_verified=True)
    ce_nested = dict(id='ce_nested', reaction_ids=['A0', 'A1'],
                     match_level='complete_equation', complete_chemical_equation_verified=True)
    ce_free = dict(id='ce_free', reaction_ids=['B0', 'B1'],
                   match_level='complete_equation', complete_chemical_equation_verified=True)
    ce_sets = [frozenset(p['reaction_ids']) for p in (ce_big, ce_nested, ce_free)]
    assert _free_standing_ce_crumb(ce_free, ce_sets)
    assert not _free_standing_ce_crumb(ce_nested, ce_sets)
    assert not _free_standing_ce_crumb(ce_big, ce_sets)


def test_ce_crumb_overlaps_large_non_ce_helper():
    crumb = dict(id='ce', reaction_ids=['R1', 'R2'],
                 match_level='complete_equation', complete_chemical_equation_verified=True)
    geo4 = dict(id='geo4', reaction_ids=['R1', 'R2', 'R3', 'R4'], match_level='geometry')
    geo2 = dict(id='geo2', reaction_ids=['R1', 'R2'], match_level='geometry')
    assert _ce_crumb_overlaps_large_non_ce(crumb, [crumb, geo4])
    assert not _ce_crumb_overlaps_large_non_ce(crumb, [crumb, geo2])


def _ce_crumb(i, r0, r1, score=5.66):
    return dict(id=f'ce{i}', source='src_b', reaction_ids=[r0, r1], score=score,
                features=_feat(2.0), layout_origin={},
                match_level='complete_equation', complete_chemical_equation_verified=True)


def _geo4(i, rids, score=7.85):
    return dict(id=f'geo{i}', source='src_c', reaction_ids=list(rids), score=score,
                features=_feat(0.0), layout_origin={}, match_level='geometry')


def test_sparse_pool_keeps_score_order():
    """Fewer than CE_CRUMB_LOCK_MIN_COMPETING conflicts → no CE-crumb lock."""
    # 3 competing free-standing CE crumbs (below threshold 8).
    crumbs = [_ce_crumb(i, f'R{2*i}', f'R{2*i+1}') for i in range(3)]
    geos = [_geo4(i, [f'R{2*i}', f'R{2*i+1}', f'G{i}a', f'G{i}b'], score=20.0) for i in range(3)]
    # Disjoint denser CE pack — wins on score without a density lock.
    ce_pack = dict(id='ce_pack', source='src_b',
                   reaction_ids=['X0', 'X1', 'X2', 'X3'], score=12.0,
                   features=_feat(2.0), layout_origin={},
                   match_level='complete_equation', complete_chemical_equation_verified=True)
    pool = crumbs + geos + [ce_pack]
    ce_sets = [frozenset(p['reaction_ids']) for p in pool if _is_complete_equation(p)]
    competing = _competing_free_standing_ce_crumbs(pool, ce_sets)
    assert len(competing) == 3
    assert len(competing) < CE_CRUMB_LOCK_MIN_COMPETING
    chosen = select(pool)
    chosen_ids = {c['id'] for c in chosen}
    assert 'ce_pack' in chosen_ids
    # No lock: higher-scoring size-4 geometry beats the overlapping CE crumbs.
    assert all(f'geo{i}' in chosen_ids for i in range(3))
    assert all(f'ce{i}' not in chosen_ids for i in range(3))


def test_dense_competition_locks_all_ce_crumbs():
    """≥ CE_CRUMB_LOCK_MIN_COMPETING competing free-standing crumbs → lock all size≤2 CE."""
    n = CE_CRUMB_LOCK_MIN_COMPETING
    crumbs = [_ce_crumb(i, f'R{2*i}', f'R{2*i+1}') for i in range(n)]
    geos = [_geo4(i, [f'R{2*i}', f'R{2*i+1}', f'G{i}a', f'G{i}b'], score=20.0) for i in range(n)]
    # Nested CE crumb (subset of a larger CE) is also locked under density rule.
    ce_big = dict(id='ce_big', source='src_b',
                  reaction_ids=[f'N{i}' for i in range(4)], score=15.0,
                  features=_feat(2.0), layout_origin={},
                  match_level='complete_equation', complete_chemical_equation_verified=True)
    ce_nested = dict(id='ce_nested', source='src_b', reaction_ids=['N0', 'N1'], score=5.5,
                     features=_feat(2.0), layout_origin={},
                     match_level='complete_equation', complete_chemical_equation_verified=True)
    pool = crumbs + geos + [ce_big, ce_nested]
    ce_sets = [frozenset(p['reaction_ids']) for p in pool if _is_complete_equation(p)]
    competing = _competing_free_standing_ce_crumbs(pool, ce_sets)
    assert len(competing) == n
    chosen = select(pool)
    chosen_ids = {c['id'] for c in chosen}
    # Locked crumbs beat higher-scoring overlapping geometry.
    assert all(f'ce{i}' in chosen_ids for i in range(n))
    assert all(f'geo{i}' not in chosen_ids for i in range(n))
    # Local exchange still prefers the larger CE pack over its nested crumb.
    assert 'ce_big' in chosen_ids
    assert 'ce_nested' not in chosen_ids


def test_select_does_not_fragment_larger_ce_in_sparse_pool():
    ce_big = dict(id='ce_big', source='src_b',
                  reaction_ids=[f'A{i}' for i in range(8)], score=20.0,
                  features=_feat(2.0), layout_origin={},
                  match_level='complete_equation', complete_chemical_equation_verified=True)
    ce_nested = dict(id='ce_nested', source='src_b', reaction_ids=['A0', 'A1'], score=5.5,
                     features=_feat(2.0), layout_origin={},
                     match_level='complete_equation', complete_chemical_equation_verified=True)
    assert [c['id'] for c in select([ce_nested, ce_big])] == ['ce_big']


def test_lock_rule_has_no_source_or_model_branches():
    import inspect
    from retrieval import select as select_fn
    src = inspect.getsource(select_fn)
    for needle in ("model_id", "e_coli", "iIT341", "n_target", "source =="):
        assert needle not in src


if __name__ == '__main__':
    test_scaled_max_roots()
    test_independent_reaction_sets()
    test_select_demotes_geometry_crumbs()
    test_complete_equation_bonus_constant()
    test_free_standing_ce_crumb_helper()
    test_ce_crumb_overlaps_large_non_ce_helper()
    test_sparse_pool_keeps_score_order()
    test_dense_competition_locks_all_ce_crumbs()
    test_select_does_not_fragment_larger_ce_in_sparse_pool()
    test_lock_rule_has_no_source_or_model_branches()
    print('PASS')
