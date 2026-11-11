"""Shim: run the geometry report from experiments/geometry_report.py."""
import runpy
import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parent.parent
if __name__ == '__main__':
    sys.path.insert(0, str(_ROOT / 'layout_atlas'))
    sys.path.insert(0, str(_ROOT))
    sys.argv[0] = str(_ROOT / 'experiments' / 'geometry_report.py')
    runpy.run_path(str(_ROOT / 'experiments' / 'geometry_report.py'), run_name='__main__')
