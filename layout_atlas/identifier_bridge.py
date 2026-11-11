"""Runtime access to `data/identifier_bridge.json` (built by `tools/build_identifier_bridge.py`): BiGG metabolite -> ChEBI and BiGG reaction -> EC number."""
import json
import re
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

FILE = Path(__file__).resolve().parent.parent / 'data/identifier_bridge.json'
FULL_EC = re.compile(r'^\d+\.\d+\.\d+\.\d+$')


@lru_cache(maxsize=1)
def _bridge():
    return json.loads(FILE.read_text())


def metabolite_chebi():
    """BiGG metabolite id -> set of ChEBI ids (union over the models)."""
    return {k: set(v) for k, v in _bridge()['metabolite_chebi'].items()}


def bigg_to_ec():
    """BiGG reaction id -> set of fully specified EC numbers (partial numbers such as 1.1.1.- are ignored), union over the models."""
    table = defaultdict(set)
    for bigg, ecs in _bridge()['reaction_ec'].items():
        table[bigg] |= {e for e in ecs if FULL_EC.match(e)}
    return table


def sources():
    return _bridge()['sources']
