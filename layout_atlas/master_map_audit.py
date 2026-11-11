"""Shim: master map audit lives in experiments.master_map_audit."""
import runpy
import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
from experiments.master_map_audit import *  # noqa: F401,F403
if __name__ == '__main__':
    sys.argv[0] = str(_ROOT / 'experiments' / 'master_map_audit.py')
    runpy.run_path(str(_ROOT / 'experiments' / 'master_map_audit.py'), run_name='__main__')
