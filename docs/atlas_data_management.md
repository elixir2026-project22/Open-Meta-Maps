# Atlas data management audit

Status: 2026-11-11. Engineering review of how extracted layouts are stored and registered. Not legal advice.

## Inventory

| Store | Role |
|---|---|
| `data/open_core.sqlite` | Redistributable Atlas (tracked). Tables: `maps`, `modules`, `reactions`, `reviews`, `retained_observations`. |
| `data/atlas.sqlite` | Runtime working copy (gitignored). Copied from `open_core.sqlite` on first run (`core.py`). |
| `data/extensions/*/modules.json.gz` | Default extensions (BRENDA, Human-maps, WikiPathways geometry, Metabolism Regulation, PANTHER) loaded by `database.modules()`. |
| `data/source_inputs/` | Observations, raw inputs, `source_registry.json`, terms snapshots. |
| `data/input_manifest.json` | SHA-256 of files under `source_inputs/` (version label `0.3-source-inputs`). |
| `data/acquisition.json` | Acquisition ledger (Reactome diagrams, reference paths, failures). |
| `data/packages.json` | Optional release packages (atlas-full). |

`database.build()` rebuilds `atlas.sqlite` from observations + raw GPML/Reactome; it does not currently write a schema version row or a machine-readable registry of the resulting open core.

## Provenance completeness (open core)

### Maps (2,323)

| Field | escher (22) | reactome (86) | wikipathways (2,215) |
|---|---|---|---|
| `source` / `license` / `source_sha256` / `species` / `lineage` | present | present | present |
| `source_id` | present | **missing** | present |
| `url` / `version` / `date` / `taxon_id` / `attribution` | **missing** | **missing** | **missing** |
| `provenance_status` | `manual_history_unverified` (22) | **missing** | **missing** |
| `layout_origin` | **missing** | **missing** | **missing** |

Registry rows in `source_registry.json` already hold `url`, `version`, `sha256`, `license` for most files; map payloads do not consistently copy them.

### Modules (17,199)

Payloads always duplicate the indexed columns (`id`, `source`, `species`, `lineage`, `quality`, `size`, `core_family`, `variant`) inside the zlib JSON blob. Every module has `license` and `source_sha256`. None of the open-core modules carry `layout_origin` or `provenance_status` (extensions do: BRENDA `possibly_automatic`, Human-maps / Metabolism Regulation `expert_hand_drawn`).

### `provenance_status = import_or_propagation_unresolved`

Found on **2,215** WikiPathways entries in `map_observations.json.gz` (not on sqlite map rows). Meaning:

- At observation import time the stable public pathway id was not written into `lineage` (opaque 24-char hash instead).
- `lineage_hint` already holds the WikiPathways id (`WP####`) for every such row.
- The open-core sqlite maps **already use** `lineage=WP####` (resolved at GPML rebuild). Observation rows still say unresolved and still record snapshot version `20261010` after the 20261110 refresh — stale relative to the DB.

Escher maps use `manual_history_unverified` (22): honest; upstream community-maps history was not independently re-verified map-by-map.

### Reactions / reviews / retained_observations

- `reactions`: identifier index only (no licence). Healthy.
- `reviews`: schema present, **0 rows** (unused).
- `retained_observations` (2,938): exclusion ledger (`disconnected_visible_core`, `outside_v03_size_bound`, `reactome_replaced_by_fixed_v97`). Table exists in the DB file but was omitted from `database.connect()` DDL.

## Duplicates and copies

- **Lineage ≠ unique map key.** 2,323 maps → 2,320 distinct `lineage` strings. The three collisions are related Escher pathway pairs (e.g. iAF1260 vs iJO1366 “Central metabolism”, two RECON1 carbohydrate spellings). Modules intentionally share lineage within a family (up to ~1.1k modules per lineage).
- **`atlas.sqlite` vs `open_core.sqlite`:** identical when both present after a fresh copy; `atlas.sqlite` is gitignored and may drift if `database.build()` is run locally.
- **Column/payload duplication:** every module column is repeated inside `payload`. Acceptable for self-contained exports; costly for migrations.
- **Near-duplicate modules:** many `core_family`/`variant` buckets hold multiple instances by design (geometry variants / multi-map coverage), not accidental clones.

## Registry / reproducibility gaps

- No **schema version** inside sqlite; `database.build()` summary uses a hard-coded `version: 0.3` string only in a JSON side file under `layout_atlas/results/`.
- Indexes: `core_family_idx`, `biological_idx` only. Missing `source` / `lineage` indexes on `maps` and `modules`.
- `source_registry.json`: 22 Escher paths still point at stale `../modular_layout/data/raw/...` instead of `data/source_inputs/raw/escher/...`.
- `input_manifest.json`: missing `wikipathways_refresh.json` after the WP refresh (1 file on disk not in manifest).
- `freeze()` refuses to refresh if the manifest already exists — rebuilds are not driven from a single up-to-date manifest.
- Extensions are outside sqlite; only `SOURCE.json` records their provenance (good for BRENDA/Human-maps; inconsistent field names `license_url` vs `license_uri`).

## WikiPathways 20261110 refresh verification

Compared `open_core.sqlite` at the parent of commit `e7e3a4f` vs that commit.

| | Before | After |
|---|---:|---:|
| maps escher / reactome / wikipathways | 22 / 86 / **934** | 22 / 86 / **2215** |
| modules escher / reactome / wikipathways | 9098 / 3405 / **4718** | 9098 / 3405 / **4696** |
| modules total | 17,221 | 17,199 (−22) |
| reactions total | 63,274 | 63,274 |
| WP reactions (join modules) | 4,879 | 4,879 |
| WP quality | geometry_only 3473, semantic 1245 | geometry_only 3451, semantic 1245 |

What actually happened:

1. **Archives were refreshed.** All 41 species zips were re-downloaded for `20261110` (every archive SHA changed). `raw/wp` holds 2,215 GPML files.
2. **Maps table was rebuilt for WikiPathways.** Opaque-lineage rows (934) were replaced with one row per GPML (`lineage=WP####`, 2,215). File-hash overlap with the previous maps table: **929** unchanged, **5** old hashes gone, **1,286** new hashes.
3. **Semantic modules were re-cut, not expanded.** All 1,245 `semantic_transfer` module ids were replaced, but they still come from the **same 65** `source_sha256` values as before. Unique WP reaction identifiers (`namespace`,`identifier`) are identical (338). Row count in `reactions` for WP stayed 4,879 because the same semantic content was reindexed under new module ids (4,879 rows deleted and 4,879 inserted; the full reactions table is **not** byte-identical).
4. **Geometry-only modules were mostly retained, not regenerated.** 3,451 kept by matching `source_sha256` still present in the new snapshot; **22 dropped** (Apis cerana urea-cycle + four yeast pathways) because those five GPML contents changed hash and were not re-extracted. **No new geometry-only modules** were created for the 1,286 newly registered maps.
5. **`map_observations.json.gz` was not updated** (still `version=20261010` for all 2,215 WP rows; still `import_or_propagation_unresolved`).
6. **Consistency bug:** all 3,451 retained geometry-only modules still pointed at **old** `map_id` values that no longer exist in `maps` (0 hits). `source_sha256` still matched a new map row.

So the first refresh ingested new GPML into `maps` and re-derived semantic modules from unchanged bindable files, but it was a **partial** Atlas rebuild. That was corrected by a **complete rebuild** (below).

## WikiPathways 20261110 complete rebuild

Pipeline: `layout_atlas/wikipathways_rebuild.py` over **all 2,215** GPML files (generic; no species/map special cases).

| | Partial refresh (prior) | Complete rebuild |
|---|---:|---:|
| WP maps | 2215 | 2215 |
| WP semantic modules (open_core) | 1245 | **1245** (re-cut from GPML) |
| WP geometry_only modules | 3451 retained / 22 dropped | **29,673** (re-cut; default extension) |
| Maps with ≥1 module | ~929 | **934** (65 semantic + 869 geometry) |
| Maps without modules | ~1286 | **1281** (no Rhea-validated events and no metabolite–metabolite edges — typically gene/signaling drawings) |
| WP reactions indexed | 4879 | **4879** (geometry_only has no chemical ids) |
| Reactions total | 63274 | **63274** |
| `map_observations` version | 20261010 / unresolved | **20261110** / `resolved_pathway_id` |
| Dangling `map_id` | 3451 | **0** |

Changed GPMLs (Apis urea cycle; yeast NAD / sulfate / ubiquinol / pyrimidine salvage) are re-cut into geometry_only modules. Licence: WikiPathways content is **CC0**; GPML files have no per-file `License` attribute; Atlas stamps `CC0-1.0`.

Geometry-only windows are stored in `data/extensions/wikipathways_geometry/` (loaded by default via `database.modules()`) so `open_core.sqlite` stays under GitHub’s size limit; semantic WP modules remain in sqlite.

## Reasonable fixes (this PR)

1. Document findings (this file), including partial-refresh verification and the complete rebuild.
2. Remap (interim) then **fully re-extract** WikiPathways modules from all 20261110 GPML; regenerate observations / registry / manifest.
3. Schema v1: `meta` table, source/lineage indexes; `ensure_schema()` on connect.
4. Normalize open-core provenance fields; write `data/atlas_registry.json`.
5. Validation test: fail on missing licence / source / `layout_origin` / schema version / dangling module `map_id`.

Out of scope (kept lean): dropping payload column duplication; deleting the empty `reviews` table; taxon-id coverage for every species string; treating gene-only WikiPathways drawings as metabolic geometry modules.
