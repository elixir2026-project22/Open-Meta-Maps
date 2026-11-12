"""Translation-only curated axes + whole-block relocate_sub."""
import numpy as np

import curated_axes as CA
import core


def test_fit_frame_default_no_rotation():
    src = [[0, 0], [2, 0], [0, 2], [2, 2]]
    # Target is rotated 90° and scaled; translation-only must keep identity rotation.
    tgt = [[0, 0], [0, 4], [-4, 0], [-4, 4]]
    frame = core.fit_frame(src, tgt)
    assert frame['rotate'] is False
    assert frame['reflect'] is False
    rot = np.asarray(frame['rotation'], float)
    assert abs(np.linalg.det(rot) - 1.0) < 1e-9
    assert abs(rot[0, 0] - 1.0) < 1e-9 and abs(rot[1, 1] - 1.0) < 1e-9


def test_fit_frame_rotate_forbids_reflection():
    src = [[0, 0], [1, 0], [0, 1]]
    tgt = [[0, 0], [-1, 0], [0, 1]]  # mirror
    frame = core.fit_frame(src, tgt, rotate=True)
    assert frame['reflect'] is False
    assert np.linalg.det(frame['rotation']) > 0


def test_reapply_same_map_common_scale_and_offsets():
    # Two modules from one map; independent fits had different scales / absolute poses.
    selected = [
        dict(id='a', map_id='M',
             source_xy={'r:1': [0.0, 0.0], 'r:2': [10.0, 0.0]},
             positions={'r:1': [100.0, 50.0], 'r:2': [120.0, 50.0]},
             source_to_initial_frame=dict(scale=2.0, rotation=np.eye(2).tolist())),
        dict(id='b', map_id='M',
             source_xy={'r:3': [0.0, 20.0], 'r:4': [10.0, 20.0]},
             positions={'r:3': [0.0, 0.0], 'r:4': [30.0, 0.0]},
             source_to_initial_frame=dict(scale=3.0, rotation=np.eye(2).tolist())),
    ]
    out, report = CA.reapply_same_map_layout(selected)
    assert report['maps'] == 1 and report['modules'] == 2
    assert abs(report['common_scales']['M'] - 2.5) < 1e-9
    # Relative centroid offset must match source (× common scale).
    ca = np.mean([out[0]['positions'][k] for k in ('r:1', 'r:2')], axis=0)
    cb = np.mean([out[1]['positions'][k] for k in ('r:3', 'r:4')], axis=0)
    # Source centroids: (5,0) and (5,20) → delta (0, 20) × 2.5 = (0, 50)
    assert abs((cb - ca)[0]) < 1e-6
    assert abs((cb - ca)[1] - 50.0) < 1e-6
    for p in out:
        assert p['source_to_initial_frame']['rotate'] is False
        assert p['source_to_initial_frame']['reflect'] is False


def test_axes_report_kept_when_identity():
    blocks = [dict(id='x', map_id='M',
                   source_to_initial_frame=dict(rotation=np.eye(2).tolist(), reflect=False, scale=1.0))]
    rep = CA.axes_report(blocks)
    assert rep['blocks'] == 1 and rep['axes_kept'] == 1 and rep['mirrored'] == 0


def test_bake_same_map_keeps_internal_and_rows():
    # Two modules; internal spans must not change, centroid delta must follow source direction.
    blocks = [
        dict(id='a', map_id='M', physical_scale=100.0, reaction_ids=['1', '2'],
             source_xy={'r:1': [0.0, 0.0], 'r:2': [10.0, 0.0]},
             positions={'r:1': [0.0, 0.0], 'r:2': [100.0, 0.0]},
             source_to_initial_frame=dict(scale=1.0, rotation=np.eye(2).tolist())),
        dict(id='b', map_id='M', physical_scale=100.0, reaction_ids=['3', '4'],
             source_xy={'r:3': [0.0, 20.0], 'r:4': [10.0, 20.0]},
             positions={'r:3': [500.0, 500.0], 'r:4': [600.0, 500.0]},
             source_to_initial_frame=dict(scale=1.0, rotation=np.eye(2).tolist())),
        dict(id='c', map_id='M', physical_scale=100.0, reaction_ids=['5', '6'],
             source_xy={'r:5': [0.0, 40.0], 'r:6': [10.0, 40.0]},
             positions={'r:5': [800.0, 800.0], 'r:6': [900.0, 800.0]},
             source_to_initial_frame=dict(scale=1.0, rotation=np.eye(2).tolist())),
    ]
    before_span = np.ptp([blocks[0]['positions']['r:1'], blocks[0]['positions']['r:2']], axis=0)
    out, rep = CA.bake_same_map_positions(blocks, model_n_reactions=6)
    assert rep['groups'] == 1
    after_span = np.ptp([out[0]['positions']['r:1'], out[0]['positions']['r:2']], axis=0)
    assert np.allclose(before_span, after_span)
    ca = np.mean([out[0]['positions'][k] for k in ('r:1', 'r:2')], axis=0)
    cb = np.mean([out[1]['positions'][k] for k in ('r:3', 'r:4')], axis=0)
    # Source centroid delta is (0, 20) — purely vertical; common scale ≈ 10.
    assert abs((cb - ca)[0]) < 1e-6
    assert abs((cb - ca)[1] - 200.0) < 1e-4
    assert all(b.get('same_map_rigid') == 'M' for b in out)


def test_compress_padding_shortens_empty_runs_keeps_order():
    # Two columns with a tall empty band between y=0 and y=1000; dense spacing elsewhere.
    cloud = {
        'r:A': [0.0, 0.0], 'r:B': [100.0, 0.0],
        'r:C': [0.0, 100.0], 'r:D': [100.0, 100.0],
        'r:E': [0.0, 1000.0], 'r:F': [100.0, 1000.0],
        'r:G': [0.0, 1100.0], 'r:H': [100.0, 1100.0],
    }
    before = {k: list(v) for k, v in cloud.items()}
    pad = CA.compress_padding_xy(cloud, edge_lengths=[100.0] * 20)
    assert pad['applied'] is True
    assert pad['gaps_y'] >= 1
    assert pad['removed_y'] > 0
    # Order and x alignment preserved; empty band collapsed toward median spacing.
    assert cloud['r:A'][1] < cloud['r:C'][1] < cloud['r:E'][1] < cloud['r:G'][1]
    assert abs(cloud['r:A'][0] - before['r:A'][0]) < 1e-9
    assert (cloud['r:E'][1] - cloud['r:C'][1]) < (before['r:E'][1] - before['r:C'][1]) - 100
    # Replay remap matches in-place compression.
    assert CA.apply_padding_xy(before['r:E'], pad) == cloud['r:E']


def test_bake_same_map_absolute_rigid_from_sha(tmp_path, monkeypatch):
    # Absolute source coords place two fragments with original relative offset.
    sha = 'abc' * 21 + 'ab'  # 65 chars? need 64 hex — use real-looking
    sha = 'a' * 64
    # edge_lengths large enough that the 200px row gap is not treated as padding.
    abs_map = dict(reactions={'R1': [0.0, 0.0], 'R2': [100.0, 0.0],
                              'R3': [0.0, 200.0], 'R4': [100.0, 200.0]},
                   metabolites={'a': [50.0, 0.0], 'b': [50.0, 200.0]}, path='test',
                   edge_lengths=[400.0] * 30)
    monkeypatch.setattr(CA, 'absolute_source_coords', lambda s: abs_map if s == sha else None)
    blocks = [
        dict(id='a', map_id='M', source_sha256=sha, reaction_ids=['R1', 'R2'],
             source_xy={'r:R1': [0.0, 0.0], 'r:R2': [1.0, 0.0]},  # collapsed local
             positions={'r:R1': [10.0, 10.0], 'r:R2': [110.0, 10.0], 'm:a': [60.0, 10.0]},
             source_to_initial_frame=dict(scale=1.0, rotation=np.eye(2).tolist())),
        dict(id='b', map_id='M', source_sha256=sha, reaction_ids=['R3', 'R4'],
             source_xy={'r:R3': [0.0, 0.0], 'r:R4': [1.0, 0.0]},
             positions={'r:R3': [900.0, 900.0], 'r:R4': [1000.0, 900.0], 'm:b': [950.0, 900.0]},
             source_to_initial_frame=dict(scale=1.0, rotation=np.eye(2).tolist())),
        dict(id='c', map_id='M', source_sha256=sha, reaction_ids=['R1'],  # tiny third for share gate
             source_xy={'r:R1': [0.0, 0.0]},
             positions={'r:R1': [0.0, 0.0]},
             source_to_initial_frame=dict(scale=1.0, rotation=np.eye(2).tolist())),
    ]
    # Fix third block to have its own reactions so absolute rebuild has 3 blocks
    blocks[2] = dict(id='c', map_id='M', source_sha256=sha, reaction_ids=['R1', 'R2'],
                     source_xy={'r:R1': [0.0, 0.0], 'r:R2': [1.0, 0.0]},
                     positions={'r:R1': [0.0, 0.0], 'r:R2': [50.0, 0.0]},
                     source_to_initial_frame=dict(scale=1.0, rotation=np.eye(2).tolist()))
    # Actually blocks a and c share reactions — use distinct R5/R6 for c
    abs_map['reactions']['R5'] = [200.0, 0.0]
    abs_map['reactions']['R6'] = [300.0, 0.0]
    blocks[2] = dict(id='c', map_id='M', source_sha256=sha, reaction_ids=['R5', 'R6'],
                     source_xy={'r:R5': [0.0, 0.0], 'r:R6': [1.0, 0.0]},
                     positions={'r:R5': [0.0, 0.0], 'r:R6': [50.0, 0.0]},
                     source_to_initial_frame=dict(scale=1.0, rotation=np.eye(2).tolist()))
    out, rep = CA.bake_same_map_positions(blocks, model_n_reactions=6)
    assert rep['groups'] == 1
    assert 'absolute_rigid' in rep['mode'].values()
    assert abs(list(rep['scales'].values())[0] - 1.0) < 1e-9
    ca = np.mean([out[0]['positions']['r:R1'], out[0]['positions']['r:R2']], axis=0)
    cb = np.mean([out[1]['positions']['r:R3'], out[1]['positions']['r:R4']], axis=0)
    # Source centroid delta (0, 200) at scale 1:1 (no padding collapse with large median).
    assert abs((cb - ca)[0]) < 1e-6
    assert abs((cb - ca)[1] - 200.0) < 1e-6
    assert out[0].get('same_map_snapshot')
    assert 'padding_compression' in rep
    # Restore undoes relative drift.
    scene = dict(positions={k: list(v) for b in out for k, v in b['positions'].items()}, participants={})
    scene['positions']['r:R3'][0] += 400
    scene['positions']['r:R4'][0] += 400
    scene2, rr = CA.restore_rigid_group(scene, out)
    assert rr['restored'] > 0
    ca2 = np.mean([scene2['positions']['r:R1'], scene2['positions']['r:R2']], axis=0)
    cb2 = np.mean([scene2['positions']['r:R3'], scene2['positions']['r:R4']], axis=0)
    assert abs((cb2 - ca2)[0]) < 1.0
    assert abs((cb2 - ca2)[1] - 200.0) < 1.0


def test_clear_rigid_node_overlaps_nudges_non_rigid():
    blocks = [dict(id='b', same_map_rigid='M', reaction_ids=['R1'],
                   same_map_snapshot={'r:R1': [0.0, 0.0], 'm:met': [10.0, 0.0]},
                   positions={'r:R1': [0.0, 0.0], 'm:met': [10.0, 0.0]})]
    scene = dict(
        positions={'r:R1': [0.0, 0.0], 'm:met': [10.0, 0.0], 'm:other': [10.5, 0.2]},
        participants={('R1', 'met'): 'm:met', ('ALG', 'other'): 'm:other'},
    )
    scene2, rep = CA.clear_rigid_node_overlaps(scene, blocks, min_dist=12.0)
    assert rep['moved'] >= 1
    d = np.hypot(*(np.asarray(scene2['positions']['m:other']) - np.asarray(scene2['positions']['m:met'])))
    assert d >= 12.0 - 1e-6
    # Rigid glyph stays put.
    assert scene2['positions']['m:met'] == [10.0, 0.0]


def test_spatial_reaction_pieces_split_distant_clusters():
    # Two compact clusters far apart → two continuity pieces.
    xy = {
        'A': (0.0, 0.0), 'B': (10.0, 0.0), 'C': (0.0, 10.0),
        'D': (500.0, 500.0), 'E': (510.0, 500.0), 'F': (500.0, 510.0),
    }
    pieces = CA.spatial_reaction_pieces(list(xy), xy, adjacency_factor=1.5)
    assert len(pieces) == 2
    sizes = sorted(len(p) for p in pieces)
    assert sizes == [3, 3]


def test_expert_structure_entries_emit_spatial_pieces():
    blocks = [
        dict(id='mega', reaction_ids=['A', 'B', 'C', 'D', 'E', 'F'],
             positions={
                 'r:A': [0.0, 0.0], 'r:B': [10.0, 0.0], 'r:C': [0.0, 10.0],
                 'r:D': [500.0, 500.0], 'r:E': [510.0, 500.0], 'r:F': [500.0, 510.0],
             }),
        dict(id='tiny', reaction_ids=['X', 'Y'],
             positions={'r:X': [0.0, 0.0], 'r:Y': [5.0, 0.0]}),
    ]
    rows = CA.expert_structure_entries(blocks, codes={'mega': 'U1', 'tiny': 'U2'})
    assert any(r['id'].startswith('mega#cc') for r in rows)
    assert any(r['id'] == 'tiny' for r in rows)
    assert all(len(r['reactions']) >= 2 for r in rows)


def test_whole_block_unit_covers_all_reactions():
    import block_relocation as br
    model = dict(reactions=[
        dict(id='EX_mal_e', metabolites={'mal__L_e': -1}),
        dict(id='MALt2_2', metabolites={'mal__L_e': 1, 'mal__L_c': -1}),
        dict(id='FUM', metabolites={'mal__L_c': 1, 'fum_c': -1}),
        dict(id='OTHER', metabolites={'x_c': -1, 'y_c': 1}),
    ], metabolites=[])
    scene = dict(
        positions={
            'r:EX_mal_e': [0.0, 0.0], 'r:MALt2_2': [50.0, 0.0], 'r:FUM': [500.0, 0.0], 'r:OTHER': [800.0, 0.0],
            'm:mal__L_e': [25.0, 10.0], 'm:mal__L_c': [200.0, 10.0], 'm:fum_c': [520.0, 10.0],
            'm:x_c': [780.0, 10.0], 'm:y_c': [820.0, 10.0],
        },
        participants={
            ('EX_mal_e', 'mal__L_e'): 'm:mal__L_e',
            ('MALt2_2', 'mal__L_e'): 'm:mal__L_e',
            ('MALt2_2', 'mal__L_c'): 'm:mal__L_c',
            ('FUM', 'mal__L_c'): 'm:mal__L_c',
            ('FUM', 'fum_c'): 'm:fum_c',
            ('OTHER', 'x_c'): 'm:x_c',
            ('OTHER', 'y_c'): 'm:y_c',
        },
    )
    block = dict(id='blk', reaction_ids=['EX_mal_e', 'MALt2_2', 'FUM'])
    long_hint = [('EX_mal_e', 'mal__L_c', 450.0, 'm:mal__L_c')]
    unit = br._whole_block_unit(model, scene, block, cur=set(), long_hint=long_hint, spacing=100.0)
    assert unit is not None
    assert set(unit['reactions']) == {'EX_mal_e', 'MALt2_2', 'FUM'}
    assert 'r:EX_mal_e' in unit['keys'] and 'r:FUM' in unit['keys']
    assert unit['whole_block'] is True


def test_compress_empty_axis_edges_shortens_empty_corridor():
    # Vertical mid↔met with empty corridor; a second column occupies the same Y band.
    cloud = {
        'r:V': [0.0, 0.0], 'm:top': [0.0, 800.0],
        'r:S': [400.0, 0.0], 'r:S2': [400.0, 400.0], 'r:S3': [400.0, 800.0],
    }
    edges = [('r:V', 'm:top')]
    out = CA.compress_empty_axis_edges(cloud, edges, max_gap=150.0)
    assert out['applied'] and out['shortened'] == 1
    assert abs(cloud['m:top'][1] - cloud['r:V'][1]) <= 150.0 + 1e-6
    # Side column untouched.
    assert cloud['r:S3'][1] == 800.0


def test_compress_padding_joint_exposes_breaks():
    cloud = {
        'r:A': [0.0, 0.0], 'r:B': [0.0, 100.0],
        'r:C': [0.0, 1000.0], 'r:D': [0.0, 1100.0],
        'm:x': [0.0, 0.0],
    }
    edges = [('r:A', 'r:C')]  # keys must exist; use r-r as synthetic
    # Use proper r/m pair
    cloud = {
        'r:A': [0.0, 0.0], 'm:t': [0.0, 1000.0],
        'r:B': [300.0, 0.0], 'r:C': [300.0, 500.0], 'r:D': [300.0, 1000.0],
    }
    edges = [('r:A', 'm:t')]
    rep = CA.compress_padding_joint(cloud, edge_lengths=[100.0]*20, edges=edges, rounds=3)
    assert rep['applied']
    assert 'breaks_x' in rep and 'breaks_y' in rep
