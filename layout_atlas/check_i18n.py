"""Shim: i18n check lives in tools.check_i18n."""
import runpy
import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
if __name__ == '__main__':
    sys.argv[0] = str(_ROOT / 'tools' / 'check_i18n.py')
    runpy.run_path(str(_ROOT / 'tools' / 'check_i18n.py'), run_name='__main__')
