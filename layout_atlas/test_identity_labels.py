"""Smoke tests for label→identity bridging and PANTHER CellDesigner decode."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))

from core import DATA
import identity_labels
import celldesigner
import database


def main():
    ident, meta = identity_labels.resolve_label('Glucose')
    assert ident and ident.startswith('chebi:'), ident
    assert meta['status'] == 'chebi_bridge'
    ident2, meta2 = identity_labels.resolve_label('2 ATP')
    assert ident2 and 'atp' in ident2.lower() or ident2.startswith('chebi:'), ident2
    assert meta2.get('multiplicity') == 2
    path = DATA / 'source_inputs/raw/panther/Glycolysis.xml'
    assert path.exists(), path
    d = celldesigner.decode(path, source='panther', license='CC-BY-4.0', name='Glycolysis')
    assert len(d['glyphs']) >= 20, len(d['glyphs'])
    assert len(d['events']) >= 8, (len(d['events']), d['audit'])
    mods = database.diagram_modules(d)
    assert mods, 'no modules from Glycolysis'
    ext = DATA / 'extensions/panther/modules.json.gz'
    assert ext.exists(), ext
    print(f"PASS identity/panther: glycolysis glyphs={len(d['glyphs'])} events={len(d['events'])} modules={len(mods)}")


if __name__ == '__main__':
    main()
