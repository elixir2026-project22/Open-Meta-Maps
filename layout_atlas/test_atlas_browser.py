"""Tests for Atlas Browser search, dedup and download manifest (pytest)."""
from __future__ import annotations

import csv
import io
import json
import tempfile
import zipfile
from pathlib import Path

import atlas_browser as ab


def _mod(mid, *, source='escher', name='Demo.Map', map_id='m1', license='CC0-1.0',
         attribution='Demo attribution', rx_id='PFK', refs=None, ec=None, mets=None, size=2,
         quality='semantic_transfer', extra_rx=None):
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
    reactions = [dict(key='r:1', id=rx_id, references=refs, coefficients={mets[0]: -1.0, mets[-1]: 1.0}, reversible=True)]
    for i, rid in enumerate(extra_rx or [], start=2):
        rk = f'r:{i}'
        nodes.append(dict(key=rk, role='reaction', xy=[float(i), 0.0], entity=rid, observed=True, position_origin='observed'))
        reactions.append(dict(key=rk, id=rid, references={'bigg.reaction': [rid]}, coefficients={}, reversible=False))
    return dict(
        id=mid, source=source, name=name, map_id=map_id, lineage=f'{source}:{map_id}',
        species='Escherichia coli', license=license, attribution=attribution,
        size=size if extra_rx is None else 1 + len(extra_rx or []),
        quality=quality, nodes=nodes, edges=edges, reactions=reactions,
        source_sha256='abc', original_instance='test',
    )


class TestSearchIndex:
    def setup_method(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / 'index.json.gz'
        mods = [
            _mod('mod_pfk', rx_id='TESTKX', ec=['9.9.9.9'], mets=['g6p_c', 'f6p_c'], name='e_coli_core.Glycolysis', map_id='e_coli_core'),
            _mod('mod_gapd', rx_id='GAPDX', ec=['1.2.1.12'], mets=['g3p_c', '13dpg_c'], name='e_coli_core.Glycolysis', map_id='e_coli_core', source='escher'),
            _mod('mod_human', rx_id='MAR0001', refs={'human.gem': ['MAR0001'], 'bigg.reaction': ['HEX1']},
                 ec=['2.7.1.1'], mets=['glc_D_c', 'g6p_c'], source='human_maps', license='CC-BY-4.0',
                 attribution='Human-maps (CC BY 4.0)', name='Human-GEM.Glycolysis', map_id='hm1'),
            # Same canonical reaction set as mod_pfk, different map/source → also_in after collapse.
            _mod('mod_pfk_alt', rx_id='TESTKX', refs={'bigg.reaction': ['TESTKX']}, ec=['9.9.9.9'],
                 mets=['g6p_c', 'f6p_c'], name='iJO1366.Glycolysis', map_id='iJO1366', source='escher',
                 attribution='Escher community-maps (CC0)'),
            # Second window on the same source map as mod_pfk (same reactions, larger size).
            _mod('mod_pfk_win2', rx_id='TESTKX', ec=['9.9.9.9'], mets=['g6p_c', 'f6p_c'],
                 name='e_coli_core.Glycolysis', map_id='e_coli_core', size=5),
            # Larger module ranks first.
            _mod('mod_big', rx_id='A', extra_rx=['B', 'C', 'D', 'E', 'F', 'G'], mets=['x_c', 'y_c'], name='Big.Map', map_id='big1'),
            _mod('mod_bad', license='CC-BY-SA-4.0'),
        ]
        ab._index = None
        self._bridge = ab._reaction_ec_bridge
        ab._reaction_ec_bridge = lambda: {}
        self.summary = ab.build_index(self.path, modules=[('core', m) for m in mods])
        self.idx = ab.load_index(self.path, force=True)

    def teardown_method(self):
        ab._reaction_ec_bridge = self._bridge
        ab._index = None
        self.tmp.cleanup()

    def test_skips_non_cc0_ccby_mit(self):
        assert self.summary['modules'] == 6
        assert self.summary['skipped_license'] == 1
        assert not ab.license_ok('CC-BY-SA-4.0')
        assert ab.license_ok('CC0-1.0-data')

    def test_search_by_ec_groups_source_maps(self):
        r = ab.search({'q': '9.9.9.9'}, index=self.idx)
        assert r['q_kind'] == 'ec'
        assert r['module_total'] == 3  # pfk, pfk_alt, pfk_win2
        assert r['source_map_total'] == 2  # e_coli_core.Glycolysis + iJO1366.Glycolysis
        # Cross-source identical reaction set collapses to one card with also_in.
        assert r['total'] == 1
        card = r['results'][0]
        assert card['windows_count'] >= 2  # pfk + win2 on e_coli_core (or alt alone)
        # Representative is the largest window on that map (win2 size 5 > pfk size 2).
        assert card['id'] == 'mod_pfk_win2'
        assert card['size'] == 5
        assert card['size_label'] == '5 reactions'
        assert card['also_in']['maps'] == 2
        assert card['also_in']['label'] == 'also in: 1 source / 2 maps'
        assert '9.9.9.9' in ' '.join(card['matched'])
        assert card['subtitle'].startswith(ab.source_label(card['source']))
        # Expanded windows list is size-desc and labels size once.
        assert [w['size'] for w in card['windows']] == sorted((w['size'] for w in card['windows']), reverse=True)
        assert card['windows'][0]['size_label'] == '5 reactions'

    def test_show_all_windows_expands(self):
        r = ab.search({'q': '9.9.9.9', 'variants': '1'}, index=self.idx)
        assert r['module_total'] == 3
        assert r['total'] == 3
        assert r['deduped'] is False

    def test_ranks_larger_first(self):
        r = ab.search({}, index=self.idx)
        assert r['results'][0]['id'] == 'mod_big'
        assert r['results'][0]['size'] == 7

    def test_search_by_reaction_canonical_and_source(self):
        r = ab.search({'q': 'HEX1'}, index=self.idx)
        assert r['q_kind'] == 'reaction'
        assert r['total'] == 1
        assert r['results'][0]['id'] == 'mod_human'
        assert 'HEX1' in r['results'][0]['matched'] or 'MAR0001' in r['results'][0]['matched']

    def test_search_by_metabolite_and_map(self):
        r = ab.search({'q': 'g6p_c'}, index=self.idx)
        assert r['q_kind'] == 'metabolite'
        assert r['module_total'] >= 2
        r2 = ab.search({'q': 'Human-GEM.Glycolysis'}, index=self.idx)
        assert r2['total'] >= 1
        r3 = ab.search({'source': 'human_maps', 'license': 'CC-BY-4.0'}, index=self.idx)
        assert r3['module_total'] == 1

    def test_pagination_and_count_label(self):
        r = ab.search({}, page=1, per_page=2, index=self.idx)
        assert r['module_total'] == 6
        assert r['source_map_total'] >= 1
        assert 'source maps' in r['count_label']
        assert 'modules' in r['count_label']
        assert len(r['results']) == 2

    def test_each_result_card_has_thumbnail(self):
        r = ab.search({}, page=1, per_page=24, index=self.idx)
        assert r['results']
        for card in r['results']:
            assert card.get('thumbnail'), card['id']
            assert card['thumbnail'].endswith('/thumbnail.svg')
            assert f"/api/atlas/modules/{card['id']}/thumbnail.svg" == card['thumbnail']

    def test_free_text_q_auto_detects_ec_reaction_metabolite_map(self):
        r_ec = ab.search({'q': '9.9.9.9'}, index=self.idx)
        assert r_ec['q_kind'] == 'ec'
        assert r_ec['q_kind_label'] == 'Matched as EC number'
        r_rx = ab.search({'q': 'HEX1'}, index=self.idx)
        assert r_rx['q_kind'] == 'reaction'
        r_met = ab.search({'q': 'g6p_c'}, index=self.idx)
        assert r_met['q_kind'] == 'metabolite'
        r_map = ab.search({'q': 'e_coli_core'}, index=self.idx)
        assert r_map['q_kind'] == 'map'
        assert r_map['module_total'] >= 1
        r_txt = ab.search({'q': 'Glycolysis'}, index=self.idx)
        assert r_txt['q_kind'] == 'text'
        assert r_txt['module_total'] >= 1


class TestQueryAndSpeciesHelpers:
    def test_classify_query(self):
        assert ab.classify_query('1.1.1.100') == 'ec'
        assert ab.classify_query('1.1.1.-') == 'ec'
        assert ab.classify_query('HEX1') == 'reaction'
        assert ab.classify_query('rhea:12345') == 'reaction'
        assert ab.classify_query('g6p_c') == 'metabolite'
        assert ab.classify_query('chebi:17925') == 'metabolite'
        assert ab.classify_query('e_coli_core') == 'map'
        assert ab.classify_query('iAF1260') == 'map'
        assert ab.classify_query('iAF1260.Fatty_Acid_Biosynthesis') == 'map'
        assert ab.classify_query('central carbon metabolism') == 'text'
        assert ab.classify_query('') == 'text'

    def test_species_display_hides_placeholders(self):
        assert ab.species_display('Escherichia coli') == 'Escherichia coli'
        assert ab.species_display('pathway-level; organism assignment not in source SVG') is None
        assert ab.species_display('unspecified') is None
        assert ab.species_display('') is None
        assert ab.species_display(None) is None

    def test_public_card_omits_unknown_species(self):
        m = _mod('brenda_x', source='brenda', name='Pathway')
        m['species'] = 'pathway-level; organism assignment not in source SVG'
        card = ab.public_card(m)
        assert card['species'] is None
        html = (ab.ROOT / 'layout_atlas' / 'results' / 'atlas_browser' / 'index.html').read_text(encoding='utf-8')
        assert "d.species?'" in html or 'var spRow=d.species?' in html

    def test_thumbnail_uses_atlas_examples_style(self):
        m = _mod('thumb')
        svg = ab.thumbnail_svg(m)
        assert 'stroke="#111"' in svg
        assert 'fill="#111111"' in svg
        assert 'fill="#ffffff"' in svg
        assert 'vector-effect="non-scaling-stroke"' in svg
        assert 'preserveAspectRatio="xMidYMid meet"' in svg
        # Frame grey lives on .ab-th so meet letterboxing has no visible side bands.
        html = (ab.ROOT / 'layout_atlas' / 'results' / 'atlas_browser' / 'index.html').read_text(encoding='utf-8')
        assert '.ab-th{' in html and 'background:#fafafa' in html


class TestDownloadManifest:
    def setup_method(self):
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
        self._orig = ab.load_module
        ab.load_module = lambda mid, index=None: self.mods.get(mid)

    def teardown_method(self):
        ab.load_module = self._orig
        ab._index = None
        self.tmp.cleanup()

    def test_manifest_json_has_cc_by_attribution_text(self):
        text, rows, missing = ab.download_manifest(['a', 'b'], index=self.idx, fmt='json')
        assert not missing
        assert len(rows) == 2
        by = {r['module_id']: r for r in rows}
        assert 'Attribution required' in by['a']['attribution_text']
        assert 'CC-BY-4.0' in by['a']['attribution_text']
        assert 'public domain' in by['b']['attribution_text'].lower()
        payload = json.loads(text)
        assert len(payload['modules']) == 2

    def test_manifest_csv_columns(self):
        text, rows, missing = ab.download_manifest(['a'], index=self.idx, fmt='csv')
        assert not missing
        row = next(csv.DictReader(io.StringIO(text)))
        assert row['module_id'] == 'a'
        assert row['license'] == 'CC-BY-4.0'
        assert 'Attribution required' in row['attribution_text']

    def test_zip_contains_escher_and_manifests(self):
        blob, rows, missing = ab.build_download_zip(['a', 'c'], index=self.idx)
        assert len(rows) == 2
        with zipfile.ZipFile(io.BytesIO(blob)) as zf:
            names = set(zf.namelist())
            assert 'manifest.json' in names
            assert 'manifest.csv' in names
            assert any(n.startswith('modules/') and n.endswith('.json') for n in names)
            data = json.loads(zf.read('modules/a.json'))
            assert len(data) == 2
            assert data[0]['metadata']['license'] == 'CC-BY-4.0'


class TestEscherExportAndThumbnail:
    def test_no_rotation_or_mirror(self):
        m = _mod('geo')
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
        for n in nodes.values():
            if n['node_type'] == 'metabolite':
                by_key[n['bigg_id']] = n
            else:
                by_key['rx'] = n
        ax = by_key['a_c']['x'] - by_key['rx']['x']
        ay = by_key['a_c']['y'] - by_key['rx']['y']
        bx = by_key['b_c']['x'] - by_key['rx']['x']
        by = by_key['b_c']['y'] - by_key['rx']['y']
        assert ax * by - ay * bx > 0

    def test_thumbnail_has_circles_and_lines(self):
        m = _mod('thumb')
        svg = ab.thumbnail_svg(m)
        assert '<circle' in svg
        assert '<line' in svg
        assert 'viewBox="' in svg
        assert 'preserveAspectRatio="xMidYMid meet"' in svg
        assert 'vector-effect="non-scaling-stroke"' in svg
        assert 'stroke-width="1.2"' in svg
        assert 'stroke="#111"' in svg

    def test_thumbnail_viewbox_pads_content(self):
        m = _mod('thumb')
        # Known coords: reaction (0,0), mets (1,0) and (0,1) → bbox 0..1,0..1
        svg = ab.thumbnail_svg(m)
        head = svg.split('>', 1)[0]
        assert 'preserveAspectRatio="xMidYMid meet"' in head
        import re
        vb = list(map(float, re.search(r'viewBox="([^"]+)"', head).group(1).split()))
        # At least 6% pad, plus room for node radius so circles are not clipped.
        assert vb[0] <= -0.06 and vb[1] <= -0.06
        assert vb[2] >= 1.12 and vb[3] >= 1.12

    def test_count_label_pluralisation(self):
        assert ab._count_label(1, 'reaction') == '1 reaction'
        assert ab._count_label(2, 'reaction') == '2 reactions'
        assert ab._count_label(1, 'source') == '1 source'
        assert ab._count_label(4, 'map') == '4 maps'

    def test_quality_plain_language(self):
        q = ab.quality_label('semantic_transfer')
        assert q['label'] == 'Curated reaction layout'
        assert 'curated' in q['tip'].lower()

    def test_attribution_dict_normalised(self):
        text = ab.attribution_string(
            {'title': 'BRENDA pathway map', 'doi': '10.1/x', 'url': 'https://example.test'},
            'brenda',
        )
        assert 'BRENDA pathway map' in text
        assert 'doi:10.1/x' in text
        assert '[object Object]' not in text

    def test_page_html_simplified_form(self):
        ab.page_html()
        html = (ab.ROOT / 'layout_atlas' / 'results' / 'atlas_browser' / 'index.html').read_text(encoding='utf-8')
        assert 'id="ec"' not in html
        assert 'id="reaction"' not in html
        assert 'id="metabolite"' not in html
        assert 'id="map"' not in html
        assert '1.1.1.100, HEX1, g6p_c, iJO1366' in html
        assert 'qDetect' in html
        assert 'windows' in html
        assert 'max-height:4.8em' in html
        assert 'aspect-ratio:3/2' in html
        assert 'object-fit:contain' in html
        assert 'reactionLabel' in html
        assert 'size_label' in html
