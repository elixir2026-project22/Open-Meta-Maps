"""Unit tests for metabolite-name matching, EC currentisation and BRENDA map acceptance rules."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import identifier_bridge as ib
import ec_channel as ec


class NameKeys(unittest.TestCase):
    def test_stereo_and_anomer(self):
        self.assertEqual(ib.name_key('D-Glucose'), ib.name_key('glucose'))
        self.assertEqual(ib.name_key('alpha-D-glucose'), ib.name_key('D-glucose'))
        # positional alpha/beta kept
        self.assertNotEqual(ib.name_key('beta-alanine'), ib.name_key('alanine'))

    def test_bis_vs_di_phosphate(self):
        self.assertEqual(ib.name_key('fructose-1,6-bisphosphate'), ib.name_key('D-Fructose 1,6-diphosphate'))

    def test_nadph_expansion(self):
        keys = set().union(*(ib.name_keys(n) for n in ib.expand_names(['NAD(P)H'])))
        self.assertTrue(ib.name_keys('NADH') & keys)
        self.assertTrue(ib.name_keys('NADPH') & keys)

    def test_membrane_side(self):
        self.assertEqual(ib.name_key('H+[side 1]'), ib.name_key('H+'))

    def test_racemase_two_keys(self):
        # L- and D- share the stereo-stripped key but keep a stereo key each
        a, b = ib.name_keys('L-alanine'), ib.name_keys('D-alanine')
        self.assertTrue(a & b)
        self.assertNotEqual(a, b)


class WeightedJaccard(unittest.TestCase):
    def test_one_to_one(self):
        # anomerase-like: two stereo forms must not both match a single alanine
        left = [ib.metabolite_keys(names=['L-alanine']), ib.metabolite_keys(names=['D-alanine'])]
        right = [ib.metabolite_keys(names=['alanine'])]
        score = ib.weighted_jaccard(left, right)
        self.assertLess(score, 1.0)
        self.assertGreater(score, 0.0)

    def test_currency_half_weight(self):
        left = [ib.metabolite_keys(names=['ATP']), ib.metabolite_keys(names=['glucose'])]
        right = [ib.metabolite_keys(names=['ATP']), ib.metabolite_keys(names=['glucose'])]
        self.assertAlmostEqual(ib.weighted_jaccard(left, right), 1.0)

    def test_accept_rules(self):
        self.assertTrue(ib.accept_score(0.6, linked=True))
        self.assertFalse(ib.accept_score(0.7, linked=False))  # needs 0.8 unlinked
        self.assertTrue(ib.accept_score(0.85, linked=False))
        self.assertFalse(ib.accept_score(0.7, linked=True, best=0.9))  # gap > 0.1


class EnzymeCurrent(unittest.TestCase):
    def test_transferred(self):
        # classic transfer used in models: 1.6.5.3 -> 7.1.1.2
        cur = set(ec.current_ecs(['1.6.5.3']))
        self.assertIn('7.1.1.2', cur)
        self.assertNotIn('1.6.5.3', cur)

    def test_deleted_dropped(self):
        table = ec.enzyme_table()
        deleted = next(ec_ for ec_, e in table.items() if e['deleted'])
        self.assertEqual(set(ec.current_ecs([deleted])), set())

    def test_drop_covered(self):
        self.assertEqual(ec.drop_covered(['2.7.1.-', '2.7.1.11']), ['2.7.1.11'])


class HandCheckSample(unittest.TestCase):
    """A few e_coli_core reactions whose BRENDA nodes are checkable by EC + metabolites."""

    @classmethod
    def setUpClass(cls):
        from model_io import load_model
        import brenda_library as bl
        if not (Path(__file__).resolve().parent.parent / 'data/extensions/brenda/modules.json.gz').exists():
            raise unittest.SkipTest('BRENDA extension not installed')
        cls.model = load_model(Path(__file__).resolve().parent.parent / 'examples/e_coli_core.xml')
        cls.matches = bl.match_model(cls.model, vicinity='propagate')

    def test_pgi_on_map(self):
        m = self.matches['PGI']
        self.assertEqual(m['status'], 'on map')
        ecs = {c['ec'] for c in m['candidates'] if c['accepted']}
        self.assertTrue(any(e and e.startswith('5.3.1.9') for e in ecs) or any('5.3.1.9' in (c.get('ec') or '') for c in m['candidates'] if c['accepted']))

    def test_pfk_on_map(self):
        self.assertEqual(self.matches['PFK']['status'], 'on map')

    def test_exchange_boundary(self):
        self.assertEqual(self.matches['EX_glc__D_e']['status'], 'boundary')


if __name__ == '__main__':
    unittest.main()
