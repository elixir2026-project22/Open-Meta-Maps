"""Round-trip export tests: Escher map.json → SBGN-ML, SBML Layout, CellDesigner SBML.

Uses e_coli_core and iIT341 only. Every model reaction and species must appear in
each export; SBGN-ML validates against the libSBGN 0.3 XSD; SBML validates with
libSBML (no errors); coordinates match map.json; the layout-benchmark SBGN reader
parses the file.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))

from core import EXAMPLES, read
from model_io import load_model
import map_exports
from evaluation.reference_sbgn import read_sbgn_glyphs, read_sbgn_positions

MODELS = [
    ('e_coli_core', EXAMPLES / 'e_coli_core.xml'),
    ('iIT341', EXAMPLES / 'iIT341.xml'),
]
CACHE = ROOT / 'layout_atlas' / 'results' / 'export_test_maps'


def _compose_if_needed(name: str, gem: Path) -> Path:
    """Return a folder with map.json (compose once, reuse across tests)."""
    out = CACHE / name
    if (out / 'map.json').is_file() and (out / 'map.sbgn').is_file():
        return out
    out.mkdir(parents=True, exist_ok=True)
    # Drop stale partial results
    for p in out.iterdir():
        if p.is_file():
            p.unlink()
        elif p.is_dir() and p.name in ('upload', 'full', 'views'):
            shutil.rmtree(p, ignore_errors=True)
    cmd = [
        sys.executable, str(ROOT / 'project.py'), 'compose',
        '--gem', str(gem),
        '--out', str(out),
        '--source-paths',
    ]
    proc = subprocess.run(cmd, cwd=str(ROOT), capture_output=True, text=True)
    if proc.returncode != 0:
        raise AssertionError(
            f'compose failed for {name} (code {proc.returncode})\n'
            f'stdout tail:\n{(proc.stdout or "")[-2000:]}\n'
            f'stderr tail:\n{(proc.stderr or "")[-2000:]}'
        )
    assert (out / 'map.json').is_file(), f'missing map.json in {out}'
    # Exports should already be written by the pipeline; rewrite if needed
    if not (out / 'map.sbgn').is_file():
        upload = next(out.glob('upload/model.*'))
        map_exports.write_all_exports(out, model_path=upload)
    return out


@pytest.fixture(scope='module', params=MODELS, ids=[m[0] for m in MODELS])
def composed(request):
    name, gem = request.param
    assert gem.is_file(), gem
    folder = _compose_if_needed(name, gem)
    model = load_model(gem)
    data = read(folder / 'map.json')
    return name, gem, folder, model, data


def test_export_files_exist(composed):
    name, gem, folder, model, data = composed
    for fname in ('map.json', 'map.sbgn', 'map.sbml', 'map.celldesigner.xml'):
        assert (folder / fname).is_file(), f'{name}: missing {fname}'


def test_sbgn_xsd_and_reader(composed):
    name, gem, folder, model, data = composed
    ok, errors = map_exports.validate_sbgn_xsd(folder / 'map.sbgn')
    assert ok, f'{name}: SBGN XSD errors: {errors[:10]}'
    glyphs = read_sbgn_glyphs(folder / 'map.sbgn')
    positions = read_sbgn_positions(folder / 'map.sbgn')
    assert glyphs and positions, f'{name}: layout-benchmark reader got no glyphs'
    assert any(cls == 'process' for cls, _ in glyphs.values())
    assert any(cls in ('simple chemical', 'macromolecule') for cls, _ in glyphs.values())


def test_sbml_layout_validates(composed):
    name, gem, folder, model, data = composed
    ok, errors = map_exports.validate_sbml_libsbml(folder / 'map.sbml')
    assert ok, f'{name}: SBML errors: {errors[:10]}'
    # Model chemistry unchanged: same species/reaction counts as input
    import libsbml
    original = libsbml.readSBMLFromFile(str(gem))
    exported = libsbml.readSBMLFromFile(str(folder / 'map.sbml'))
    assert original.getModel().getNumSpecies() == exported.getModel().getNumSpecies()
    assert original.getModel().getNumReactions() == exported.getModel().getNumReactions()
    plugin = exported.getModel().getPlugin('layout')
    assert plugin is not None and plugin.getNumLayouts() >= 1


def test_every_reaction_and_species_in_exports(composed):
    name, gem, folder, model, data = composed
    model_rx = {r['id'] for r in model['reactions']}
    model_sp = {m['id'] for m in model['metabolites']}

    map_rx, map_sp = map_exports.collect_export_ids(data, model)
    assert model_rx <= map_rx, f'{name}: map.json missing reactions {sorted(model_rx - map_rx)[:10]}'
    # Map nodes cover reaction participants; orphan species (if any) are added only in exports.
    assert map_sp, f'{name}: map.json has no metabolite nodes'

    sbgn_rx, sbgn_sp = map_exports.sbgn_ids(folder / 'map.sbgn')
    assert model_rx <= sbgn_rx, f'{name}: SBGN missing reactions {sorted(model_rx - sbgn_rx)[:10]}'
    assert model_sp <= sbgn_sp, f'{name}: SBGN missing species {sorted(model_sp - sbgn_sp)[:10]}'

    lay_rx, lay_sp = map_exports.layout_ids(folder / 'map.sbml')
    assert model_rx <= lay_rx, f'{name}: SBML Layout missing reactions {sorted(model_rx - lay_rx)[:10]}'
    assert model_sp <= lay_sp, f'{name}: SBML Layout missing species {sorted(model_sp - lay_sp)[:10]}'

    cd_rx, cd_sp = map_exports.celldesigner_ids(folder / 'map.celldesigner.xml')
    assert model_rx <= cd_rx, f'{name}: CellDesigner missing reactions {sorted(model_rx - cd_rx)[:10]}'
    assert model_sp <= cd_sp, f'{name}: CellDesigner missing species {sorted(model_sp - cd_sp)[:10]}'


def test_coordinates_match_map_json(composed):
    name, gem, folder, model, data = composed
    assert map_exports.coordinates_match(data, folder / 'map.sbgn', tol=1.5), (
        f'{name}: SBGN glyph centres do not match map.json'
    )


def test_exports_do_not_alter_map_json(composed):
    """Re-export must leave map.json bytes unchanged (layout metrics stay the same)."""
    name, gem, folder, model, data = composed
    before = (folder / 'map.json').read_bytes()
    map_exports.write_all_exports(folder, model_path=gem, model=model)
    after = (folder / 'map.json').read_bytes()
    assert before == after, f'{name}: map.json changed during export'


def test_synthetic_bezier_and_clones():
    """Unit-level: curved segments flatten to next-points; clones get a clone marker."""
    model = {
        'id': 'mini', 'name': 'mini', 'genes': [], 'annotations': [], 'notes': '',
        'metabolites': [
            dict(id='A_c', name='A', compartment='c', annotations=[]),
            dict(id='B_c', name='B', compartment='c', annotations=[]),
            dict(id='h_c', name='H+', compartment='c', annotations=[]),
        ],
        'reactions': [
            dict(id='R1', name='R1', metabolites={'A_c': -1, 'B_c': 1, 'h_c': 1},
                 lower_bound=-1000, upper_bound=1000, gene_reaction_rule='', annotations=[], notes=''),
        ],
    }
    data = [
        dict(schema='https://escher.github.io/escher/jsonschema/1-0-0#', map_name='mini',
             map_id='mini', homepage='', map_description='', metadata={}),
        dict(
            nodes={
                '0': dict(node_type='metabolite', x=0, y=100, bigg_id='A_c', name='A',
                          node_is_primary=True, label_x=0, label_y=80),
                '1': dict(node_type='metabolite', x=200, y=100, bigg_id='B_c', name='B',
                          node_is_primary=True, label_x=200, label_y=80),
                '2': dict(node_type='metabolite', x=200, y=40, bigg_id='h_c', name='H+',
                          node_is_primary=False, label_x=200, label_y=20),
                '3': dict(node_type='metabolite', x=50, y=40, bigg_id='h_c', name='H+',
                          node_is_primary=False, label_x=50, label_y=20),
                'm': dict(node_type='midmarker', x=100, y=100),
            },
            reactions={
                '0': dict(
                    bigg_id='R1', name='R1', reversibility=True, label_x=100, label_y=120,
                    gene_reaction_rule='', genes=[],
                    metabolites=[
                        dict(bigg_id='A_c', coefficient=-1),
                        dict(bigg_id='B_c', coefficient=1),
                        dict(bigg_id='h_c', coefficient=1),
                    ],
                    segments={
                        '0': dict(from_node_id='m', to_node_id='0', b1=None, b2=None),
                        '1': dict(from_node_id='m', to_node_id='1',
                                  b1=dict(x=130, y=80), b2=dict(x=170, y=80)),
                        '2': dict(from_node_id='m', to_node_id='2', b1=None, b2=None),
                    },
                ),
            },
            text_labels={},
            canvas=dict(x=0, y=0, width=300, height=200),
        ),
    ]
    import tempfile
    import xml.etree.ElementTree as ET
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / 'map.sbgn'
        map_exports.write_sbgn(path, data, model)
        ok, errors = map_exports.validate_sbgn_xsd(path)
        assert ok, errors
        root = ET.parse(path).getroot()
        ns = {'s': map_exports.SBGN_NS}
        clones = root.findall('.//{http://sbgn.org/libsbgn/0.3}clone')
        assert len(clones) >= 2
        nexts = root.findall('.//{http://sbgn.org/libsbgn/0.3}next')
        assert nexts, 'expected flattened Bézier next points'
