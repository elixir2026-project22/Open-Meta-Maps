# Installation

Validated on a fresh clone (macOS, Python 3.12, Node.js, Apple clang). Linux should work with `clang++` installed; it has not been tested.

## Requirements

- Python 3.11+ (packages pinned in `layout_atlas/requirements.txt`)
- Node.js and npm (layout backend, pinned in `layout_atlas/runtime/package-lock.json`)
- `clang++` and network access for the one-time router build

## Steps

```sh
python3 -m pip install -r layout_atlas/requirements.txt
npm ci --prefix layout_atlas/runtime
python3 project.py setup
```

`setup` downloads the pinned libavoid source and builds the two edge routers into `layout_atlas/build/`. It is needed once.

On first run, if `data/atlas.sqlite` is missing it is created as a copy of `data/open_core.sqlite`.

## Run

```sh
python3 project.py tool                 # web tool, 127.0.0.1:8765
python3 project.py compose --gem examples/e_coli_core.xml --out results/e_coli_core/ --source-paths
```

The `e_coli_core` example takes roughly two minutes on a laptop (the structure-selection step dominates).

## Dependencies and bundled files

SBML files are read with [libSBML](https://sbml.org/software/libsbml/) (`python-libsbml`, installed by `pip install -r layout_atlas/requirements.txt`); COBRA / BiGG JSON is read directly (`layout_atlas/model_io.py`). The official Escher map schema is kept in `data/escher_schema/` (MIT, with its licence); the currency / cofactor list is `data/currency_metabolites.json` and the identifier tables (BiGG to ChEBI and EC) are `data/identifier_bridge.json`, both with their provenance inside the file.
