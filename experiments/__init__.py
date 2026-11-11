"""Offline experiments and reports (not loaded by default compose)."""
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
_LA = _ROOT / 'layout_atlas'
if str(_LA) not in sys.path:
    sys.path.insert(0, str(_LA))

from core import EXAMPLE_MODELS as MODELS  # noqa: F401
