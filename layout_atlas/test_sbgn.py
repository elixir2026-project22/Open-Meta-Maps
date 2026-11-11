"""Smoke tests for the generic SBGN-ML importer (python3 test_sbgn.py)."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))

from core import DATA
import sbgn
import database


def main():
    raw = DATA / 'source_inputs/raw/metabolismregulation'
    files = sorted(raw.glob('*.sbgn'))
    assert files, f'no SBGN files under {raw}'
    # Small map
    path = raw / 'F020-pyrimidine.sbgn'
    assert path.exists(), path
    d = sbgn.decode(
        path,
        source='metabolismregulation',
        license='CC-BY-4.0',
        name='Pyrimidine biosynthesis',
        species='Homo sapiens',
    )
    assert d['format'] == 'sbgn-ml'
    assert d['source_sha256']
    assert d['glyphs'], 'expected entity glyphs'
    assert d['events'], 'expected process events'
    for e in d['events']:
        assert e['visible'], e
        assert e['position_origin'] == 'observed'
        assert e['xy']
    mods = database.diagram_modules(d)
    assert mods, 'diagram_modules returned nothing'
    assert all(m['license'] == 'CC-BY-4.0' for m in mods)
    # All gallery files decode without raising
    n_events = 0
    for p in files:
        dd = sbgn.decode(p, source='metabolismregulation', license='CC-BY-4.0')
        n_events += len(dd['events'])
    assert n_events >= 50, f'expected many decoded processes, got {n_events}'
    print(f'PASS sbgn importer: {len(files)} maps, {n_events} processes, sample modules={len(mods)}')


if __name__ == '__main__':
    main()
