"""Shim: geometry report lives in experiments.geometry_report."""
import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
from experiments.geometry_report import *  # noqa: F401,F403
