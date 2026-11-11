"""Shim: run the EC experiment from experiments/ec_experiment.py."""
import runpy
import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parent.parent
if __name__ == '__main__':
    sys.path.insert(0, str(_ROOT / 'layout_atlas'))
    sys.path.insert(0, str(_ROOT))
    sys.argv[0] = str(_ROOT / 'experiments' / 'ec_experiment.py')
    runpy.run_path(str(_ROOT / 'experiments' / 'ec_experiment.py'), run_name='__main__')
