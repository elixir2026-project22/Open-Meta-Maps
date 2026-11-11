"""Shim: acquisition helpers live in tools.acquire."""
import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
from tools.acquire import *  # noqa: F401,F403
