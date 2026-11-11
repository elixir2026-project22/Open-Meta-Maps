# iCH360 custom Escher maps import

Status: imported 2026-11-11 into the default Atlas as extension `ich360`; rebased onto Stage A (`atlas-extraction-1` + compactness-aware select) and regenerated.

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

`layout_atlas/ich360_library.py` → `database.escher_diagram` + `diagram_modules` (BFS windows, sizes 1–24, `max_roots=24`). Each module carries Stage A `extraction` provenance (`method=diagram_modules`, `version=atlas-extraction-1`). `layout_origin.status` = `curated_hand_drawn`.

Provenance: `data/extensions/ich360/SOURCE.json`, `data/source_inputs/source_registry.json`, licence snapshots under `data/terms/ich360_*`. Default load includes **858** iCH360 modules (total default Atlas **69,132**).

## Compose metrics (main → ich360-maps; `--source-paths`, Stage A baseline)

| | e_coli_core main | e_coli_core +iCH360 | iIT341 main | iIT341 +iCH360 |
|---|---:|---:|---:|---:|
| curated_coverage_share | 0.7368 | 0.7368 | 0.3069 | 0.2924 |
| crossings (counted) | 46 | 46 | 299 | 239 |
| overlapping_pairs | 5 | 5 | 18 | 16 |
| unrelated_overlaps_total | 5 | 5 | 17 | 8 |
| blank_space_share | 0.5826 | 0.5826 | 0.9401 | 0.8605 |
| edge_len_max_rel | 12.187 | 12.187 | 27.975 | 27.871 |
| aspect_ratio (W/H) | 1.1726 | 1.1726 | 0.5394 | 0.8079 |
| runtime_s | 77.9 | 74.7 | 457.3 | 454.2 |
| iCH360 modules used | — | 0 | — | 7 (4×`full_map` + 3×`ions_metals_uptake`, all size-2 geometry) |

### iIT341 iCH360 selection note

`full_map` size-2 windows still win for exchange/transport motifs (e.g. native `GLNtex`+`EX_gln__L_e` geometry-matched onto iIT341 `AAtr`+`EX_aa_e`). Tighter pathway maps (`aminoacids_biosynthesis`, etc.) do not contain that periplasm/exchange pair, so they never enter the proposal pool for those targets. Stage A compactness twinning only applies when `biochemical_level ≥ 1`; these hits are geometry-only (`biol=0`), so compactness does not prefer pathway twins. The full_map window beats the previous Escher geometry cover mainly via lower `initial_fit` (−0.2×Δfit ≈ +0.011 selection_score).
