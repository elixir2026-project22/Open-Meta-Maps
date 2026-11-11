"""Validation checks for Atlas provenance and schema (python3 test_atlas_provenance.py)."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_registry import SCHEMA_VERSION, validate, ensure_schema, connect_path
from core import DATA


def main():
    db_path = DATA / 'open_core.sqlite'
    assert db_path.exists(), f'missing {db_path}'
    db = connect_path(db_path)
    ver = int(db.execute("SELECT value FROM meta WHERE key='schema_version'").fetchone()[0])
    db.close()
    assert ver >= SCHEMA_VERSION, f'schema_version {ver} < {SCHEMA_VERSION}'
    errors = validate(db_path, include_extensions=True)
    if errors:
        print(f'FAIL {len(errors)} provenance/schema issues')
        for e in errors[:50]:
            print(' ', e)
        if len(errors) > 50:
            print(f'  ... and {len(errors) - 50} more')
        raise SystemExit(1)
    print(f'PASS atlas provenance and schema_version={ver}')


if __name__ == '__main__':
    main()
