# Vendored minimal subset of sbml2escher

These files are the small subset of the `sbml2escher` reference implementation that Open Meta Maps needs at runtime:

- `models/*.xml`: symlinks to the identical example models in `examples/`
- `model_input.py` (SBML / COBRA JSON reader)
- `data/bacteria/central_reference.json` (currency-metabolite roles)
- `vendor/escher/map-schema-1-0-0.json` and `vendor/escher/LICENSE` (Escher map schema, MIT license)

`layout_atlas/core.py` uses `Ref/sbml2escher-main/` when present (the developer checkout) and falls back to this directory otherwise.

License of `model_input.py` and `central_reference.json`: **to be confirmed** (upstream license not yet verified).
