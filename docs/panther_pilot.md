# PANTHER Pathway import (CellDesigner)

Status: **imported** into the default Atlas as extension `panther` (2026-11-11).

## Source

- FTP CellDesigner CD4.1: https://data.pantherdb.org/ftp/pathway/3.4/CD4.1/
- Terms: https://pantherdb.org/tou.jsp · Licence **CC BY 4.0**
- Version: **3.4 / CD4.1**
- Raw XML: `data/source_inputs/raw/panther/` (8 metabolic pathways in this first tranche)

Note: `www.pantherdb.org` HTML endpoints still return Cloudflare 403 from some environments; the **data.pantherdb.org** FTP tree is reachable and is the acquisition path used here.

## Identity mapping

PANTHER species names are plain labels (e.g. `Glucose`), not `Label[compartment]`. The generic CellDesigner importer resolves them through `layout_atlas/identity_labels.py`, which uses the existing BiGG↔ChEBI tables in `data/identifier_bridge.json` plus metabolite display names from the public example models. Composite glyph labels (e.g. `2 NADH + 2 H+`) stay unresolved; those reactions are excluded rather than invented.

## Decode totals (8 maps)

| | Count |
|---|---:|
| Glyph aliases with centres | **362** |
| SBML reactions | 113 |
| Decoded non-branch centre-lines | **92** |
| Explicitly excluded | 21 |
| Atlas modules | **238** |

Build: `python3 layout_atlas/panther_library.py`. Modules load by default via `database.modules()`.

## Provenance

Per-file SHA-256, version, licence and attribution: `data/extensions/panther/SOURCE.json` and `data/source_inputs/source_registry.json`.
