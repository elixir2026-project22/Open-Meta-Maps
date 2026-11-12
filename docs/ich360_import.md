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

`layout_atlas/ich360_library.py` → `database.escher_diagram` + Stage B `diagram_modules` (size-scaled roots, `independent_only`, `nest_factor=1.5`, large-map min size 4) with an extra keep-`size ≥ 4` filter. Each module carries `extraction` provenance (`method=diagram_modules`, `version=atlas-extraction-1`). `layout_origin.status` = `curated_hand_drawn`.

Provenance: `data/extensions/ich360/SOURCE.json`, `data/source_inputs/source_registry.json`, licence snapshots under `data/terms/ich360_*`. Default load includes **503** iCH360 modules (sizes 4/8/12/24 only; total default Atlas **68,777**).

## Compose metrics (`--source-paths`)

Nested iCH360 cutter (main / #24) vs Stage B ich360 recut (size≥4; other Atlas sources unchanged from main):

| metric | e_coli nested | e_coli Stage B | iIT341 nested | iIT341 Stage B |
|---|---:|---:|---:|---:|
| curated_coverage_share | 0.7368 | 0.7368 | 0.2852 | 0.2852 |
| crossings (counted) | 46 | 46 | 239 | 230 |
| overlapping_pairs | 5 | 5 | 16 | 17 |
| node_overlap_pairs | 1 | 1 | 2 | 1 |
| blank_space_share | 0.5826 | 0.5826 | 0.8605 | 0.9459 |
| runtime_s | 67.3 | 69.3 | 426.9 | 432.6 |
| iCH360 modules used | 0 | 0 | 7 (all size-2 geometry) | 0 |

Pre-iCH360 Stage A docs baseline for iIT341 was `curated_coverage_share=0.3069` (compose dirs without iCH360 measured ≈0.2996). Nested #24 dipped to 0.2924 in the #24 write-up / 0.2852 when re-measured on this branch’s select path.

### iIT341 iCH360 selection note

Under the nested cutter, `full_map` / `ions_metals_uptake` **size-2 geometry-only** windows won for exchange/transport motifs (e.g. native `GLNtex`+`EX_gln__L_e` geometry-matched onto iIT341 `AAtr`+`EX_aa_e`). Pathway maps lack that periplasm/exchange pair, so they never entered the proposal pool for those targets; Stage A compactness twinning did not apply at `biol=0`.

Stage B removes those crumbs from the library (`large_map_min_size=4` + ich360 keep-`size≥4`) and demotes any leftover size≤2 geometry proposals in `retrieval.select` (−3.0). On iIT341 that yields **zero** iCH360 selections and slightly better crossings/overlaps, but **`curated_coverage_share` stays 0.2852** — the size-2 hits were geometry placements, and dropping them does not restore pre-iCH360 complete-equation coverage. Recut/dedupe alone does not undo the 0.307→0.292 dip.
