"""Synthetic checks of the orientation step (python3 test_orientation.py)."""
import numpy as np
from core import *
import orientation as o


def main():
    model = dict(id='t', metabolites=[dict(id=m, compartment='c') for m in ('a_c', 'b_c', 'atp_c', 'adp_c')],
                 reactions=[dict(id='R1', metabolites={'a_c': -1, 'atp_c': -1, 'b_c': 1, 'adp_c': 1}, lower_bound=0, upper_bound=1000)])
    scene = dict(positions={'r:R1': [100, 100], 'm:a_c': [20, 100], 'm:b_c': [180, 100], 'side:R1:atp_c': [100, 145], 'side:R1:adp_c': [100, 173]},
                 participants={('R1', 'a_c'): 'm:a_c', ('R1', 'b_c'): 'm:b_c', ('R1', 'atp_c'): 'side:R1:atp_c', ('R1', 'adp_c'): 'side:R1:adp_c'},
                 copies={'m:a_c': 'a_c', 'm:b_c': 'b_c', 'side:R1:atp_c': 'atp_c', 'side:R1:adp_c': 'adp_c'}, labels={})
    ok = True
    out = o.sides(model, scene, [])
    P = out['positions']
    # currency substrate on the substrate side (left, x < 100), currency product on the product side (right, x > 100); expert glyphs untouched
    checks = [('currency substrate left of the reaction', P['side:R1:atp_c'][0] < 100), ('currency product right of the reaction', P['side:R1:adp_c'][0] > 100),
              ('core glyphs untouched', P['m:a_c'] == [20, 100] and P['m:b_c'] == [180, 100])]
    # before: both currency copies below the reaction -> angle far below 120; after: substrates and products on opposite sides
    def angle(pos):
        subs = [np.array(pos[k]) - pos['r:R1'] for k in ('m:a_c', 'side:R1:atp_c')]
        prods = [np.array(pos[k]) - pos['r:R1'] for k in ('m:b_c', 'side:R1:adp_c')]
        return o.angle_of(subs, prods)
    checks.append(('angle before < 120 and after > 150', (angle(scene['positions']) < 120 and angle(P) > 150)))
    # a reaction no structure placed is moved between its substrate and product
    scene2 = dict(scene, positions=dict(scene['positions'], **{'r:R1': [100, 400]}))
    moved = o.free_reactions(model, scene2, [])['positions']['r:R1']
    checks.append(('free reaction moves between substrate and product', abs(moved[0] - 100) < 1 and abs(moved[1] - 100) < 1))
    # an expert reaction is not moved
    blocks = [dict(reaction_ids=['R1'])]
    checks.append(('expert reaction stays', o.free_reactions(model, scene2, blocks)['positions']['r:R1'] == [100, 400]))
    for name, good in checks:
        print(('PASS' if good else 'FAIL'), name)
        ok &= bool(good)
    print('ALL PASS' if ok else 'SOME FAILED')
    raise SystemExit(0 if ok else 1)


if __name__ == '__main__':
    main()
