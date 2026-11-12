# Metabolism Regulation gallery import

Status: imported 2026-11-11 into the default Atlas as extension `metabolismregulation`.

## Source

- Gallery: https://metabolismregulation.github.io/gallery/
- About / licence: https://metabolismregulation.github.io/about/ — **CC BY 4.0**
- Repository: https://github.com/metabolismregulation/metabolismregulation.github.io
- Pinned commit: see `data/source_inputs/raw/metabolismregulation/SOURCE_COMMIT.txt`
- Files: hand-drawn **SBGN-ML** (`.sbgn`) under upstream `downloads/` (gallery publication set) plus unique top-level `atlas-downloads/*.sbgn` stems (e.g. cardiac glycolysis). Converter variants (`converter/sbgn/*-advanced|straight`) are not imported.

## Importer

Generic SBGN-ML decoder: `layout_atlas/sbgn.py` → `database.diagram_modules`. Build / register: `layout_atlas/metabolismregulation_library.py`.

Provenance per file: source, CC BY 4.0, attribution, upstream commit version, SHA-256, and local path in `data/extensions/metabolismregulation/SOURCE.json` and `data/source_inputs/source_registry.json`.

## elixir2026-project22/data_shared

Checked 2026-11-11: https://github.com/elixir2026-project22/data_shared is **not publicly available** (GitHub 404 / `git ls-remote` repository not found). The layout-benchmark README points at that repo for SBML models (BRENDA pathways, E. coli core, iCH360, Reactome, Human2, Mycoplasma, Stress Knowledge Map) and names the Metabolism Regulation gallery as the hand-drawn reference layouts (already imported here). No additional CC0/CC BY hand-drawn layouts could be fetched from `data_shared`. Public org repos (`layout-benchmark`, `sbml2map`, `ecmapper`) do not ship other hand-drawn SBGN/Escher reference packs beyond the gallery fetch helper.

Gallery maps are imported into the default Atlas.
