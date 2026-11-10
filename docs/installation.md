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
python3 project.py compose --gem examples/e_coli_core.xml --out results/e_coli_core/ --source-paths --layout-style expert
```

The `e_coli_core` example takes roughly two minutes on a laptop (the structure-selection step dominates).

## Bundled third-party files

A minimal subset of the `sbml2escher` reference implementation (SBML reader, currency-metabolite roles, Escher map schema, example models) is vendored in `third_party/sbml2escher/`; see `NOTICE.md` there. If a `Ref/sbml2escher-main/` checkout exists next to the code it is used instead.
