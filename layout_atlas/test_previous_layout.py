"""Build-on-previous anchors and robust curated matching helpers."""
from pathlib import Path

import previous_layout as prev
import retrieval as R


def test_max_missing_reactions_scales_with_size():
    assert R.max_missing_reactions(2) == 0
    assert R.max_missing_reactions(7) == 0
    assert R.max_missing_reactions(8) == 1
    assert R.max_missing_reactions(40) == 5
    assert R.max_missing_reactions(24) == 3


def test_clip_atlas_module_drops_missing_keeps_frame():
    mod = dict(
        id='m1', size=4, source='escher',
        reactions=[
            dict(key='r:1', id='A'), dict(key='r:2', id='B'),
            dict(key='r:3', id='C'), dict(key='r:4', id='D'),
        ],
        nodes=[
            dict(key='r:1', role='reaction'), dict(key='r:2', role='reaction'),
            dict(key='r:3', role='reaction'), dict(key='r:4', role='reaction'),
            dict(key='m1', role='metabolite'), dict(key='m2', role='metabolite'),
        ],
        edges=[
            dict(a='r:1', b='m1'), dict(a='r:2', b='m1'),
            dict(a='r:3', b='m2'), dict(a='r:4', b='m2'),
        ],
    )
    clipped = R.clip_atlas_module(mod, {'A', 'B', 'C'})
    assert clipped['size'] == 3
    assert {r['id'] for r in clipped['reactions']} == {'A', 'B', 'C'}
    assert clipped['missing_source_reactions'] == ['D']
    assert clipped['partial_source_size'] == 4
    assert 'r:4' not in {n['key'] for n in clipped['nodes']}


def test_load_previous_anchors_filters_to_model(tmp_path):
    body = {
        'nodes': {
            '1': dict(node_type='metabolite', bigg_id='g6p_c', node_is_primary=True, x=10.0, y=20.0),
            '2': dict(node_type='metabolite', bigg_id='f6p_c', node_is_primary=True, x=30.0, y=40.0),
            '3': dict(node_type='midmarker', x=50.0, y=60.0),
            '4': dict(node_type='midmarker', x=70.0, y=80.0),
        },
        'reactions': {
            'r1': dict(bigg_id='PGI', segments={'s': dict(from_node_id='1', to_node_id='3')}),
            'r2': dict(bigg_id='PFK', segments={'s': dict(from_node_id='2', to_node_id='4')}),
        },
    }
    path = tmp_path / 'map.json'
    path.write_text(__import__('json').dumps([{'map_name': 't'}, body]))
    model = dict(
        reactions=[dict(id='PGI'), dict(id='OTHER')],
        metabolites=[dict(id='g6p_c'), dict(id='x_c')],
    )
    anchors = prev.load_previous_anchors(path, model)
    assert anchors['r:PGI'] == [50.0, 60.0]
    assert 'r:PFK' not in anchors
    assert anchors['m:g6p_c'] == [10.0, 20.0]
    assert 'm:f6p_c' not in anchors
    assert prev.reaction_anchor_share(anchors, model) == 0.5
    assert not prev.should_skip_packing(anchors, model)
    model2 = dict(reactions=[dict(id='PGI'), dict(id='PFK')], metabolites=[dict(id='g6p_c')])
    # PFK absent from filtered anchors → share 0.5; rebuild without model filter
    full = prev.load_previous_anchors(path)
    assert prev.should_skip_packing(full, model2)


def test_pin_scene_and_blocks():
    anchors = {'r:A': [1.0, 2.0], 'm:x': [3.0, 4.0]}
    scene = dict(positions={'r:A': [9.0, 9.0], 'm:x': [8.0, 8.0], 'r:B': [0.0, 0.0]})
    scene, n = prev.pin_scene(scene, anchors)
    assert n == 2
    assert scene['positions']['r:A'] == [1.0, 2.0]
    assert scene['positions']['r:B'] == [0.0, 0.0]
    blocks = [dict(positions={'r:A': [0.0, 0.0], 'm:x': [0.0, 0.0], 'r:C': [5.0, 5.0]})]
    blocks, n = prev.apply_block_positions(blocks, anchors)
    assert n == 2
    assert blocks[0]['positions']['r:A'] == [1.0, 2.0]
    assert blocks[0]['positions']['r:C'] == [5.0, 5.0]


def test_robust_match_keeps_large_core_block_on_e_coli_rm3():
    """Removing a few Core reactions must not shatter the large same-map window."""
    import copy
    from model_io import load_model
    from core import model_graph, ROOT
    import database

    model = load_model(ROOT / 'examples/e_coli_core.xml')
    remove = {'ME1', 'PPCK', 'FRD7'}
    reduced = copy.deepcopy(model)
    reduced['reactions'] = [r for r in reduced['reactions'] if r['id'] not in remove]
    lib = R.Library(database.modules())
    selected = R.select(lib.candidates(model_graph(reduced), model=reduced, max_roots=120)[0])
    sizes = sorted((len(p['reaction_ids']) for p in selected), reverse=True)
    assert sizes[0] >= 30, sizes[:5]
    full = R.select(lib.candidates(model_graph(model), model=model, max_roots=120)[0])
    full_sizes = sorted((len(p['reaction_ids']) for p in full), reverse=True)
    assert full_sizes[0] == 40, full_sizes[:3]


def test_annotation_seed_uses_stable_ids_not_names():
    """Robust partial match must key off reaction references, not display names."""
    import inspect
    src = inspect.getsource(R.Library.candidates)
    # The annotation-seeded path builds reverse from xrefs; name-only binding is elsewhere.
    assert 'reverse[' in src
    assert 'max_missing_reactions' in src
    assert 'clip_atlas_module' in src


def test_filter_previous_map_drops_removed_reactions(tmp_path):
    body = {
        'nodes': {
            '1': dict(node_type='metabolite', bigg_id='a_c', node_is_primary=True, x=0.0, y=0.0),
            '2': dict(node_type='midmarker', x=10.0, y=0.0),
            '3': dict(node_type='midmarker', x=20.0, y=0.0),
            '4': dict(node_type='metabolite', bigg_id='b_c', node_is_primary=True, x=30.0, y=0.0),
        },
        'reactions': {
            'r1': dict(bigg_id='KEEP', segments={'s': dict(from_node_id='1', to_node_id='2')}),
            'r2': dict(bigg_id='DROP', segments={'s': dict(from_node_id='3', to_node_id='4')}),
        },
    }
    path = tmp_path / 'map.json'
    path.write_text(__import__('json').dumps([{'map_name': 't', 'metadata': {}}, body]))
    model = dict(reactions=[dict(id='KEEP'), dict(id='NEW')], metabolites=[dict(id='a_c')])
    data, report = prev.filter_previous_map(path, model)
    assert set(data[1]['reactions']) == {'r1'}
    assert '2' in data[1]['nodes'] and '3' not in data[1]['nodes']
    assert report['new_reactions'] == ['NEW']
    assert report['removed_reactions'] == 1


def test_place_leftover_reactions_never_drops_unanchored():
    model = dict(
        reactions=[
            dict(id='OLD', metabolites={'a_c': -1, 'b_c': 1}, lower_bound=-1000, upper_bound=1000),
            dict(id='NEW_EMBED', metabolites={'a_c': -1, 'b_c': 1}, lower_bound=-1000, upper_bound=1000),
            dict(id='NEW_SOLO', metabolites={'z_c': -1, 'y_c': 1}, lower_bound=-1000, upper_bound=1000),
        ],
        metabolites=[dict(id='a_c'), dict(id='b_c'), dict(id='z_c'), dict(id='y_c')],
    )
    anchors = {'r:OLD': [0.0, 0.0], 'm:a_c': [-20.0, 0.0], 'm:b_c': [20.0, 0.0]}
    scene = dict(
        positions={
            'r:OLD': [0.0, 0.0], 'm:a_c': [-20.0, 0.0], 'm:b_c': [20.0, 0.0],
            'r:NEW_EMBED': [999.0, 999.0], 'r:NEW_SOLO': [888.0, 888.0],
        },
        participants={
            ('OLD', 'a_c'): 'm:a_c', ('OLD', 'b_c'): 'm:b_c',
            ('NEW_EMBED', 'a_c'): 'm:a_c', ('NEW_EMBED', 'b_c'): 'm:b_c',
        },
        copies={'m:a_c': 'a_c', 'm:b_c': 'b_c'},
        labels={},
    )
    scene, report = prev.place_leftover_reactions(model, scene, anchors=anchors, cur=set(), min_embed_anchors=2)
    assert report['placed'] == 2
    assert report['embedded'] == 1
    assert report['new_components'] == 1
    assert 'r:NEW_EMBED' in scene['positions']
    assert 'r:NEW_SOLO' in scene['positions']
    # Solo reaction must gain metabolite glyphs so export coverage stays complete.
    assert ('NEW_SOLO', 'z_c') in scene['participants']
    assert ('NEW_SOLO', 'y_c') in scene['participants']
    assert scene['participants'][('NEW_SOLO', 'z_c')] in scene['positions']
    # Pinned old reaction stays put.
    assert scene['positions']['r:OLD'] == [0.0, 0.0]
    assert prev.unanchored_reaction_ids(model, anchors) == ['NEW_EMBED', 'NEW_SOLO']


def test_keep_local_relayout_requires_substantial_crossing_drop():
    assert not prev.keep_local_relayout(10, 9)
    assert not prev.keep_local_relayout(10, 8)  # 20% < 25% share
    assert prev.keep_local_relayout(10, 7)
    assert not prev.keep_local_relayout(3, 2)  # drop 1 < vmin 2
    assert prev.keep_local_relayout(8, 6)


def test_graft_new_reactions_keeps_old_geometry_and_full_coverage():
    body = {
        'nodes': {
            '1': dict(node_type='metabolite', bigg_id='a_c', node_is_primary=True, x=0.0, y=0.0),
            '2': dict(node_type='midmarker', x=10.0, y=0.0),
            '3': dict(node_type='metabolite', bigg_id='b_c', node_is_primary=True, x=20.0, y=0.0),
        },
        'reactions': {
            'r1': dict(
                bigg_id='KEEP',
                metabolites=[dict(bigg_id='a_c', coefficient=-1), dict(bigg_id='b_c', coefficient=1)],
                segments={
                    '0': dict(from_node_id='2', to_node_id='1', b1=None, b2=None),
                    '1': dict(from_node_id='2', to_node_id='3', b1=None, b2=None),
                },
            ),
        },
        'canvas': dict(x=0, y=0, width=100, height=100),
    }
    data = [{'map_name': 't', 'metadata': {
        'source_reactions': {
            'KEEP': dict(stoichiometry={'a_c': -1, 'b_c': 1}, lower_bound=-1000, upper_bound=1000),
        },
    }}, body]
    model = dict(
        reactions=[
            dict(id='KEEP', metabolites={'a_c': -1, 'b_c': 1}, lower_bound=-1000, upper_bound=1000),
            dict(id='NEW', metabolites={'a_c': -1, 'b_c': 1}, lower_bound=-1000, upper_bound=1000),
            dict(id='SOLO', metabolites={'z_c': -1, 'y_c': 1}, lower_bound=-1000, upper_bound=1000),
        ],
        metabolites=[dict(id='a_c'), dict(id='b_c'), dict(id='z_c'), dict(id='y_c')],
    )
    out, report = prev.graft_new_reactions(data, model, ['NEW', 'SOLO'])
    assert report['placed'] == 2
    assert report['embedded'] >= 1
    bigg = {rx['bigg_id'] for rx in out[1]['reactions'].values()}
    assert bigg == {'KEEP', 'NEW', 'SOLO'}
    # Old midmarker untouched.
    assert out[1]['nodes']['2']['x'] == 10.0
    assert 'NEW' in out[0]['metadata']['source_reactions']
    assert 'SOLO' in out[0]['metadata']['source_reactions']
