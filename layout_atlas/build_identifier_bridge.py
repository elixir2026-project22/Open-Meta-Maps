"""Shim: identifier bridge builder lives in tools.build_identifier_bridge."""
import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
from tools.build_identifier_bridge import *  # noqa: F401,F403
if __name__ == '__main__':
    from tools.build_identifier_bridge import main
    main()
