"""Shim: low-confidence decomposition lives in experiments.optional.lowconf."""
import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
from experiments.optional.lowconf import *  # noqa: F401,F403
