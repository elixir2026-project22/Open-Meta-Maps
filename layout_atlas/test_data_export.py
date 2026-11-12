"""Tests for `python3 project.py data export` (pytest)."""
from __future__ import annotations

import json
import tarfile

import pytest

import data_packages as dp


@pytest.fixture(scope='module')
def panther(tmp_path_factory):
    tmp = tmp_path_factory.mktemp('export')
    runs = {name: tmp / f'{name}.tar.gz' for name in ('unique', 'again', 'windows')}
    dp.export_atlas(runs['unique'], sources=['panther'])
    dp.export_atlas(runs['again'], sources=['panther'])
    dp.export_atlas(runs['windows'], sources=['panther'], all_windows=True)
    return runs


def _read(path):
    with tarfile.open(path, 'r:gz') as tf:
        return {m.name: tf.extractfile(m).read() for m in tf.getmembers()}


def test_archive_layout(panther):
    files = _read(panther['unique'])
    names = sorted(files)
    assert {'README.md', 'THIRD_PARTY_NOTICES.md', 'manifest.json'} <= set(names)
    modules = [n for n in names if n.startswith('modules/')]
    assert modules and all(n.startswith('modules/panther/') and n.endswith('.json') for n in modules)
    header, body = json.loads(files[modules[0]])
    assert header['metadata']['source'] == 'panther'
    assert body['nodes'] and body['reactions']


def test_manifest_records(panther):
    files = _read(panther['unique'])
    manifest = json.loads(files['manifest.json'])
    records = manifest['modules']
    assert manifest['module_count'] == len(records) == len([n for n in files if n.startswith('modules/')])
    assert manifest['selection']['sources'] == ['panther']
    for rec in records:
        assert rec['file'] in files
        assert rec['source'] == 'panther' and rec['license'] == 'CC-BY-4.0'
        assert rec['reaction_count'] > 0 and rec['canonical_reaction_ids']
        assert rec['attribution_text'] and rec['provenance']['extraction']['method']
        assert isinstance(rec['also_in'], list) and all(p['source'] != 'panther' for p in rec['also_in'])
    sets = [tuple(r['canonical_reaction_ids']) for r in records]
    assert len(sets) == len(set(sets))


def test_all_windows_and_deterministic(panther):
    unique = json.loads(_read(panther['unique'])['manifest.json'])['module_count']
    windows = json.loads(_read(panther['windows'])['manifest.json'])['module_count']
    registry = json.loads(dp.REGISTRY.read_text())
    assert unique < windows == registry['open_core']['modules_with_extensions']['panther']
    assert panther['unique'].read_bytes() == panther['again'].read_bytes()

