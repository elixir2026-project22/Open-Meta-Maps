"""Offline/experiment entry: re-exports the production region_placement helper."""
import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parents[2]
if str(_ROOT / 'layout_atlas') not in sys.path:
    sys.path.insert(0, str(_ROOT / 'layout_atlas'))
from region_placement import *  # noqa: F401,F403
