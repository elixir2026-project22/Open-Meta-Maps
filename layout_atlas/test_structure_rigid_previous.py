"""Curated structures stay rigid under ``--previous`` pins.

For every curated structure, every member node's coordinates relative to the
structure's reference midmarker must match the curated source within 1e-6
(translation only).  Curve control points are midmarker-relative in
``source_paths`` and therefore inherit the same rigid translation when the
structure's midmarkers and metabolite ends move together.  Checked on
e_coli_core and iIT341 with and without ``--previous``.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

import previous_layout as prev

ROOT = Path(__file__).resolve().parents[1]
MODELS = ('e_coli_core', 'iIT341')
TOL = 1e-6


def _read_map(path: Path):
    return json.loads(path.read_text(encoding='utf-8'))


def _make_partial_previous(data, out_path: Path, keep_fraction: float = 0.5) -> Path:
    """Drop later curated structures so ``--previous`` is a partial earlier map."""
    header, body = data[0], data[1]
    es = list((header.get('metadata') or {}).get('expert_structures') or [])
    assert es, 'expected curated structures in map metadata'
    n_keep = max(1, int(len(es) * keep_fraction))
    drop = set()
    for s in es[n_keep:]:
        drop.update(s.get('reactions') or [])
    assert drop, 'partial previous must drop at least one structure'
    reactions = {
        rid: rx for rid, rx in (body.get('reactions') or {}).items()
        if rx.get('bigg_id') not in drop
    }
    used = set()
    for rx in reactions.values():
        for seg in (rx.get('segments') or {}).values():
            for nid in (seg.get('from_node_id'), seg.get('to_node_id')):
                if nid is not None:
                    used.add(str(nid))
    nodes = {nid: n for nid, n in (body.get('nodes') or {}).items() if nid in used}
    meta = dict(header.get('metadata') or {})
    meta['expert_structures'] = es[:n_keep]
    out = [dict(header, metadata=meta), dict(body, nodes=nodes, reactions=reactions)]
    out_path.write_text(json.dumps(out), encoding='utf-8')
    return out_path


def _scene_from_result(data):
    """Rebuild a display-like scene from an exported map + glyph_identity."""
    meta, body = data[0]['metadata'], data[1]
    gi = meta.get('glyph_identity') or {}
    nodes, rxs = body['nodes'], body['reactions']
    es = meta.get('expert_structures') or []
    rid_to_module = {}
    for s in es:
        mod = (s.get('parent') or s['id']).split('#')[0]
        for rid in s.get('reactions') or []:
            rid_to_module[rid] = mod
    copy_keys = {k for k in gi if k.startswith('copy:') and k.count(':') >= 2}
    pos, parts = {}, {}
    blocks_by_id = {}
    for rx in rxs.values():
        rid = rx.get('bigg_id')
        if not rid:
            continue
        mid_nid = None
        for seg in (rx.get('segments') or {}).values():
            for nid in (seg.get('from_node_id'), seg.get('to_node_id')):
                if (nodes.get(nid) or {}).get('node_type') == 'midmarker':
                    mid_nid = nid
                    break
            if mid_nid is not None:
                break
        if mid_nid is None:
            continue
        pos['r:' + rid] = [float(nodes[mid_nid]['x']), float(nodes[mid_nid]['y'])]
        mod = rid_to_module.get(rid)
        if mod is not None:
            blocks_by_id.setdefault(mod, dict(id=mod, reaction_ids=[]))
            if rid not in blocks_by_id[mod]['reaction_ids']:
                blocks_by_id[mod]['reaction_ids'].append(rid)
        for seg in (rx.get('segments') or {}).values():
            for nid in (seg.get('from_node_id'), seg.get('to_node_id')):
                if nid == mid_nid:
                    continue
                n = nodes.get(nid) or {}
                if n.get('node_type') != 'metabolite' or not n.get('bigg_id'):
                    continue
                met = n['bigg_id']
                key = None
                if mod is not None:
                    ck = f'copy:{mod}:{met}'
                    if ck in copy_keys:
                        key = ck
                if key is None:
                    sk = f'side:{rid}:{met}'
                    key = sk if sk in gi else 'm:' + met
                pos[key] = [float(n['x']), float(n['y'])]
                parts[(rid, met)] = key
    blocks = [dict(id=bid, reaction_ids=list(b['reaction_ids']))
              for bid, b in sorted(blocks_by_id.items())]
    return dict(positions=pos, participants=parts, copies={}, labels={}), blocks


def _structure_relatives_from_scene(scene, blocks=None):
    """Relative coords of every structure member to that structure's reference midmarker."""
    pos = scene['positions']
    membership = prev._structure_membership(scene, pos, blocks=blocks)
    out = {}
    for bid, keys in membership.items():
        keys = sorted(k for k in keys if k in pos)
        if not keys:
            continue
        ref = next((k for k in keys if k.startswith('r:')), keys[0])
        ref_xy = np.asarray(pos[ref], float)
        out[bid] = {k: (np.asarray(pos[k], float) - ref_xy).tolist() for k in keys}
    return out


def _assert_relatives_equal(a, b, label):
    assert a, f'{label}: no curated structures'
    md = 0.0
    n = 0
    for bid, ra in a.items():
        rb = b.get(bid) or {}
        for k, va in ra.items():
            if k not in rb:
                continue
            d = float(np.hypot(*(np.asarray(va, float) - np.asarray(rb[k], float))))
            md = max(md, d)
            n += 1
    assert n > 0, f'{label}: no shared structure points'
    assert md <= TOL, f'{label}: max relative deviation {md} exceeds {TOL} ({n} points)'


@pytest.mark.parametrize('model_name', MODELS)
def test_pin_scene_preserves_curated_relatives(model_name, tmp_path):
    """Curated source relatives survive ``pin_scene`` with and without ``--previous``."""
    example = ROOT / 'example_results' / model_name / f'OpenMetaMaps_{model_name}_result.json'
    if not example.is_file():
        pytest.skip(f'missing {example}')
    data = _read_map(example)
    scene, blocks = _scene_from_result(data)
    curated = _structure_relatives_from_scene(scene, blocks=blocks)
    assert curated, model_name

    # Without --previous: pin is a no-op; relatives stay identical to the curated source.
    scene_noprev = dict(
        positions={k: list(v) for k, v in scene['positions'].items()},
        participants=dict(scene['participants']),
        copies={}, labels={},
    )
    prev.pin_scene(scene_noprev, anchors={}, blocks=blocks)
    _assert_relatives_equal(
        curated, _structure_relatives_from_scene(scene_noprev, blocks=blocks),
        f'{model_name} without --previous')

    # With --previous: one structure-level translation; relatives still match curated source.
    partial = _make_partial_previous(data, tmp_path / f'{model_name}_partial.json')
    anchors = prev.load_previous_anchors(partial)
    scene_prev = dict(
        positions={k: list(v) for k, v in scene['positions'].items()},
        participants=dict(scene['participants']),
        copies={}, labels={},
    )
    prev.pin_scene(scene_prev, anchors, blocks=blocks)
    _assert_relatives_equal(
        curated, _structure_relatives_from_scene(scene_prev, blocks=blocks),
        f'{model_name} with --previous')
