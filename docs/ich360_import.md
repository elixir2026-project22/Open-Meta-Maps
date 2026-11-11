# iCH360 custom Escher maps import

Status: imported 2026-11-11 into the default Atlas as extension `ich360`.

## Source

- Repository: https://github.com/marco-corrao/iCH360
- Path: `Visualisation/` (custom-built Escher maps; Corrao et al., PLOS Comp Biol doi:10.1371/journal.pcbi.1013564)
- Pinned commit: see `data/source_inputs/raw/ich360/SOURCE_COMMIT.txt`
- Licence: **MIT** (repo-root `LICENSE`; `Visualisation/` has no overriding terms). Zenodo paper deposit of the same Visualisation tree is **CC BY 4.0** (doi:10.5281/zenodo.15048160) and corroborates open redistribution of the map JSON.
- Fetch: `python3 tools/fetch_ich360_maps.py` (checksum-pinned)

## Maps included (10)

| Path | Role |
|---|---|
| `model_map/full_map.json` | Full curated model map |
| `pathway_maps/*.json` (4) | Subsystem maps (independent geometry; reaction sets ⊂ full_map) |
| `pathways_not_in_model/*.json` (5) | Pathways outside iCH360 used for the biomass equivalent |

## Excluded

| Path | Reason |
|---|---|
| `model_map/compressed_map.json` | Uses CP pseudo-reactions; upstream `Visualisation/readme.md` says it cannot be used directly with the model |

## Dedupe rationale

Atlas already keeps overlapping drawings from one source (e.g. Human-maps compartment vs subsystem maps) and ranks retrieval candidates. Pathway maps share reaction ids with `full_map` but were measured to have substantially different curated coordinates, so they are kept as alternate drawings. `pathways_not_in_model` add reaction content absent from `full_map`.

## Importer

`layout_atlas/ich360_library.py` → `database.escher_diagram` + `diagram_modules`. `layout_origin.status` = `curated_hand_drawn`.

Provenance: `data/extensions/ich360/SOURCE.json`, `data/source_inputs/source_registry.json`, licence snapshots under `data/terms/ich360_*`.

## Compose metrics (before → after; `--source-paths`)

| | e_coli_core before | e_coli_core after | iIT341 before | iIT341 after |
|---|---:|---:|---:|---:|
| curated_coverage_share | 0.6947 | 0.6947 | 0.2329 | 0.2347 |
| crossings (counted) | 47 | 47 | 337 | 208 |
| blank_space_share | 0.5459 | 0.5459 | 0.8412 | 0.9912 |
| aspect_ratio (W/H) | 0.6771 | 0.6771 | 0.4767 | 0.8443 |
| iCH360 modules used | — | 0 | — | 4 (`full_map`) |
