"""Tests for Atlas Browser search filters and download manifest.

Run: python3 layout_atlas/test_atlas_browser.py
"""
from __future__ import annotations

import csv
import io
import json
import tempfile
import unittest
import zipfile
from pathlib import Path

import atlas_browser as ab


def _mod(mid, *, source='escher', name='Demo.Map', map_id='m1', license='CC0-1.0',
         attribution='Demo attribution', rx_id='PFK', refs=None, ec=None, mets=None, size=2):
    refs = dict(refs or {'bigg.reaction': [rx_id]})
    if ec:
        refs['ec-code'] = list(ec)
    mets = list(mets or ['g6p_c', 'f6p_c'])
    nodes = [
        dict(key='r:1', role='reaction', xy=[0.0, 0.0], entity=rx_id, observed=True, position_origin='observed'),
        dict(key='n1', role='metabolite', xy=[1.0, 0.0], entity=mets[0], label=mets[0].split('_')[0].upper(), observed=True, position_origin='observed'),
        dict(key='n2', role='metabolite', xy=[0.0, 1.0], entity=mets[-1], label=mets[-1].split('_')[0].upper(), observed=True, position_origin='observed'),
    ]
    edges = [
        dict(a='r:1', b='n1', flow='substrate', coefficient=-1.0),
        dict(a='r:1', b='n2', flow='product', coefficient=1.0),
    ]
    reactions = [dict(key='r:1', id=rx_id, references=refs, coefficients={nodes[1]['entity']: -1.0, nodes[2]['entity']: 1.0}, reversible=True)]
    return dict(
        id=mid, source=source, name=name, map_id=map_id, lineage=f'{source}:{map_id}',
        species='Escherichia coli', license=license, attribution=attribution,
        size=size, quality='semantic_transfer', nodes=nodes, edges=edges, reactions=reactions,
        source_sha256='abc', original_instance='test',
    )


class SearchIndex(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / 'index.json.gz'
        mods = [
            _mod('mod_pfk', rx_id='TESTKX', ec=['9.9.9.9'], mets=['g6p_c', 'f6p_c'], name='e_coli_core.Glycolysis', map_id='eco1'),
            _mod('mod_gapd', rx_id='GAPDX', ec=['1.2.1.12'], mets=['g3p_c', '13dpg_c'], name='e_coli_core.Glycolysis', map_id='eco1', source='escher'),
            _mod('mod_human', rx_id='MAR0001', refs={'human.gem': ['MAR0001'], 'bigg.reaction': ['HEX1']},
                 ec=['2.7.1.1'], mets=['glc_D_c', 'g6p_c'], source='human_maps', license='CC-BY-4.0',
                 attribution='Human-maps (CC BY 4.0)', name='Human-GEM.Glycolysis', map_id='hm1'),
            _mod('mod_bad', license='CC-BY-SA-4.0'),  # excluded
        ]
        ab._index = None
        # Avoid identifier_bridge EC bleed-through in unit tests.
        self._bridge = ab._reaction_ec_bridge
        ab._reaction_ec_bridge = lambda: {}
        self.summary = ab.build_index(self.path, modules=[('core', m) for m in mods])
        self.idx = ab.load_index(self.path, force=True)

    def tearDown(self):
        ab._reaction_ec_bridge = self._bridge
        ab._index = None
        self.tmp.cleanup()

    def test_skips_non_cc0_ccby_mit(self):
        self.assertEqual(self.summary['modules'], 3)
        self.assertEqual(self.summary['skipped_license'], 1)
        self.assertFalse(ab.license_ok('CC-BY-SA-4.0'))
        self.assertTrue(ab.license_ok('CC0-1.0-data'))

    def test_search_by_ec(self):
        r = ab.search({'ec': '9.9.9.9'}, index=self.idx)
        self.assertEqual(r['total'], 1)
        self.assertEqual(r['results'][0]['id'], 'mod_pfk')

    def test_search_by_reaction_canonical_and_source(self):
        r = ab.search({'reaction': 'HEX1'}, index=self.idx)
        self.assertEqual(r['total'], 1)
        self.assertEqual(r['results'][0]['id'], 'mod_human')
        r2 = ab.search({'reaction': 'MAR0001'}, index=self.idx)
        self.assertEqual(r2['total'], 1)

    def test_search_by_metabolite_and_map(self):
        r = ab.search({'metabolite': 'g6p'}, index=self.idx)
        self.assertEqual(r['total'], 2)
        r2 = ab.search({'map': 'human-gem'}, index=self.idx)
        self.assertEqual(r2['total'], 1)
        r3 = ab.search({'source': 'human_maps', 'license': 'CC-BY-4.0'}, index=self.idx)
        self.assertEqual(r3['total'], 1)

    def test_pagination(self):
        r = ab.search({}, page=1, per_page=2, index=self.idx)
        self.assertEqual(r['total'], 3)
        self.assertEqual(r['pages'], 2)
        self.assertEqual(len(r['results']), 2)


class DownloadManifest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / 'index.json.gz'
        self.mods = {
            'a': _mod('a', license='CC-BY-4.0', attribution='Human-maps (CC BY 4.0); https://example.test',
                      name='Map A', source='human_maps'),
            'b': _mod('b', license='CC0-1.0', attribution='Escher community-maps (CC0)', name='Map B'),
            'c': _mod('c', license='MIT', attribution='iCH360 (MIT)', name='Map C', source='ich360'),
        }
        ab._index = None
        ab.build_index(self.path, modules=[('core', m) for m in self.mods.values()])
        self.idx = ab.load_index(self.path, force=True)
        # Point load_module at in-memory fixtures
        self._orig = ab.load_module
        ab.load_module = lambda mid, index=None: self.mods.get(mid)

    def tearDown(self):
        ab.load_module = self._orig
        ab._index = None
        self.tmp.cleanup()

    def test_manifest_json_has_cc_by_attribution_text(self):
        text, rows, missing = ab.download_manifest(['a', 'b'], index=self.idx, fmt='json')
        self.assertFalse(missing)
        self.assertEqual(len(rows), 2)
        by = {r['module_id']: r for r in rows}
        self.assertIn('Attribution required', by['a']['attribution_text'])
        self.assertIn('CC-BY-4.0', by['a']['attribution_text'])
        self.assertIn('public domain', by['b']['attribution_text'].lower())
        payload = json.loads(text)
        self.assertEqual(len(payload['modules']), 2)

    def test_manifest_csv_columns(self):
        text, rows, missing = ab.download_manifest(['a'], index=self.idx, fmt='csv')
        self.assertFalse(missing)
        reader = csv.DictReader(io.StringIO(text))
        row = next(reader)
        self.assertEqual(row['module_id'], 'a')
        self.assertEqual(row['license'], 'CC-BY-4.0')
        self.assertIn('Attribution required', row['attribution_text'])

    def test_zip_contains_escher_and_manifests(self):
        blob, rows, missing = ab.build_download_zip(['a', 'c'], index=self.idx)
        self.assertEqual(len(rows), 2)
        with zipfile.ZipFile(io.BytesIO(blob)) as zf:
            names = set(zf.namelist())
            self.assertIn('manifest.json', names)
            self.assertIn('manifest.csv', names)
            self.assertTrue(any(n.startswith('modules/') and n.endswith('.json') for n in names))
            data = json.loads(zf.read('modules/a.json'))
            self.assertEqual(len(data), 2)
            self.assertEqual(data[0]['metadata']['license'], 'CC-BY-4.0')
            # No mirror: metabolite y order preserved after uniform scale+translate
            nodes = data[1]['nodes']
            ys = sorted(n['y'] for n in nodes.values())
            self.assertGreater(ys[-1], ys[0])


class EscherExport(unittest.TestCase):
    def test_no_rotation_or_mirror(self):
        m = _mod('geo')
        # Place three points that would change chirality if mirrored
        m['nodes'] = [
            dict(key='r:1', role='reaction', xy=[0.0, 0.0], entity='PFK', observed=True, position_origin='observed'),
            dict(key='n1', role='metabolite', xy=[2.0, 0.0], entity='a_c', label='A', observed=True, position_origin='observed'),
            dict(key='n2', role='metabolite', xy=[0.0, 1.0], entity='b_c', label='B', observed=True, position_origin='observed'),
        ]
        m['edges'] = [
            dict(a='r:1', b='n1', flow='substrate', coefficient=-1),
            dict(a='r:1', b='n2', flow='product', coefficient=1),
        ]
        data = ab.module_to_escher(m)
        nodes = data[1]['nodes']
        by_key = {}
        # Map back via bigg_id / midmarker
        for n in nodes.values():
            if n['node_type'] == 'metabolite':
                by_key[n['bigg_id']] = n
            else:
                by_key['rx'] = n
        # Cross product (n1-rx) x (n2-rx) should keep the same sign as (2,0)x(0,1)=+2
        ax = by_key['a_c']['x'] - by_key['rx']['x']
        ay = by_key['a_c']['y'] - by_key['rx']['y']
        bx = by_key['b_c']['x'] - by_key['rx']['x']
        by = by_key['b_c']['y'] - by_key['rx']['y']
        self.assertGreater(ax * by - ay * bx, 0)


if __name__ == '__main__':
    unittest.main()
