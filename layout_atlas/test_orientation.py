"""Synthetic checks of the orientation step."""
import numpy as np
from core import *
import orientation as o


def _one_reaction_scene():
    model = dict(id='t', metabolites=[dict(id=m, compartment='c') for m in ('a_c', 'b_c', 'atp_c', 'adp_c')],
                 reactions=[dict(id='R1', metabolites={'a_c': -1, 'atp_c': -1, 'b_c': 1, 'adp_c': 1}, lower_bound=0, upper_bound=1000)])
    scene = dict(positions={'r:R1': [100, 100], 'm:a_c': [20, 100], 'm:b_c': [180, 100], 'side:R1:atp_c': [100, 145], 'side:R1:adp_c': [100, 173]},
                 participants={('R1', 'a_c'): 'm:a_c', ('R1', 'b_c'): 'm:b_c', ('R1', 'atp_c'): 'side:R1:atp_c', ('R1', 'adp_c'): 'side:R1:adp_c'}, copies={}, labels={})
    return model, scene


def _two_reaction_model():
    return dict(id='t2', metabolites=[dict(id=m, compartment='c') for m in ('a_c', 'b_c', 'c_c', 'd_c', 'atp_c', 'adp_c', 'nad_c', 'nadh_c')],
                reactions=[dict(id='R1', metabolites={'a_c': -1, 'atp_c': -1, 'b_c': 1, 'adp_c': 1}, lower_bound=0, upper_bound=1000),
                           dict(id='R2', metabolites={'c_c': -1, 'nad_c': -1, 'd_c': 1, 'nadh_c': 1}, lower_bound=0, upper_bound=1000)])


def _participants():
    return {('R1', 'a_c'): 'm:a_c', ('R1', 'b_c'): 'm:b_c', ('R1', 'atp_c'): 'side:R1:atp_c', ('R1', 'adp_c'): 'side:R1:adp_c',
            ('R2', 'c_c'): 'm:c_c', ('R2', 'd_c'): 'm:d_c', ('R2', 'nad_c'): 'side:R2:nad_c', ('R2', 'nadh_c'): 'side:R2:nadh_c'}


def test_sides_places_currency_on_opposite_sides():
    model, scene = _one_reaction_scene()
    P = o.sides(model, scene, [], distance=50.0, spread=30.0)['positions']
    assert P['side:R1:atp_c'][0] < 100 and P['side:R1:adp_c'][0] > 100


def test_sides_improves_substrate_product_angle():
    model, scene = _one_reaction_scene()
    def angle(pos):
        subs = [np.array(pos[k]) - pos['r:R1'] for k in ('m:a_c', 'side:R1:atp_c')]
        prods = [np.array(pos[k]) - pos['r:R1'] for k in ('m:b_c', 'side:R1:adp_c')]
        return o.angle_of(subs, prods)
    P = o.sides(model, scene, [], distance=50.0, spread=30.0)['positions']
    assert angle(scene['positions']) < 120 and angle(P) > 150


def test_free_reactions_moves_algorithm_only():
    model, scene = _one_reaction_scene()
    scene2 = dict(scene, positions=dict(scene['positions'], **{'r:R1': [100, 400]}))
    moved = o.free_reactions(model, scene2, [])['positions']['r:R1']
    assert abs(moved[0] - 100) < 1 and abs(moved[1] - 100) < 1


def test_adaptive_spoke_length_from_local_clearance():
    model2 = _two_reaction_model(); parts = _participants()
    tight = dict(positions={'r:R1': [100, 100], 'm:a_c': [40, 100], 'm:b_c': [160, 100], 'side:R1:atp_c': [100, 100], 'side:R1:adp_c': [100, 100],
                            'r:R2': [160, 100], 'm:c_c': [130, 100], 'm:d_c': [220, 100], 'side:R2:nad_c': [160, 100], 'side:R2:nadh_c': [160, 100]},
                 participants=parts, copies={}, labels={})
    open_scene = dict(positions={'r:R1': [100, 100], 'm:a_c': [40, 100], 'm:b_c': [160, 100], 'side:R1:atp_c': [100, 100], 'side:R1:adp_c': [100, 100],
                                 'r:R2': [400, 100], 'm:c_c': [340, 100], 'm:d_c': [460, 100], 'side:R2:nad_c': [400, 100], 'side:R2:nadh_c': [400, 100]},
                      participants=parts, copies={}, labels={})
    assert o.reaction_spoke_params(model2, tight)['R1']['distance'] < o.reaction_spoke_params(model2, open_scene)['R1']['distance']
    assert abs(o.reaction_spoke_params(model2, open_scene)['R1']['distance'] - o.SIDE_MAX) < 1e-6


def test_map_spoke_floor_follows_p10():
    assert o.map_spoke_floor([37.0] * 10) <= o.SIDE_MIN + 2.5
    assert abs(o.map_spoke_floor([74.0] * 10) - o.SIDE_MAX) < 1e-6


def test_spokes_prefer_freer_side_of_axis():
    model2 = _two_reaction_model()
    crowded = dict(positions={'r:R1': [100, 100], 'm:a_c': [40, 100], 'm:b_c': [160, 100], 'side:R1:atp_c': [100, 100], 'side:R1:nad_c': [100, 100],
                              'r:R2': [100, 40], 'm:c_c': [70, 40], 'm:d_c': [130, 40]},
                   participants={('R1', 'a_c'): 'm:a_c', ('R1', 'b_c'): 'm:b_c', ('R1', 'atp_c'): 'side:R1:atp_c', ('R1', 'nad_c'): 'side:R1:nad_c', ('R2', 'c_c'): 'm:c_c', ('R2', 'd_c'): 'm:d_c'}, copies={}, labels={})
    model2['reactions'] = [dict(id='R1', metabolites={'a_c': -1, 'atp_c': -1, 'nad_c': -1, 'b_c': 1}, lower_bound=0, upper_bound=1000),
                           dict(id='R2', metabolites={'c_c': -1, 'd_c': 1}, lower_bound=0, upper_bound=1000)]
    assert o._freer_perp(np.array([100., 100.]), np.array([1., 0.]), 'R1', o._unrelated_segments(model2, crowded, currency(model2)), 80.) == 1.0
