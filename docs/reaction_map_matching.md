# Reaction–map matching evaluation

Open-data matching of GEM reactions to Atlas nodes (especially BRENDA name-only metabolites). Method ideas from [ecmapper](https://github.com/elixir2026-project22/ecmapper) (GPL-3.0; independent rewrite — no code copied). Thresholds: `identifier_bridge.THRESHOLDS` / [methodology.md](methodology.md).

Development regression models: **e_coli_core** and **iIT341** (iMM904 / iML1515 skipped). Numbers below are after rebasing onto the open-licence Atlas (#3) and import-shim cleanup (#5); BRENDA ships in the default library via `data/extensions`.

## BRENDA node match rates (non-boundary reactions)

| Model | Before (exact full EC ∩ BRENDA) | After (scored + current EC + names, no vicinity) | After + vicinity `propagate` |
|---|---:|---:|---:|
| e_coli_core | 49 / 75 (65.3%) | 49 / 75 (65.3%) | 49 / 75 (65.3%) |
| iIT341 | 323 / 477 (67.7%) | 178 / 477 (37.3%) | 205 / 477 (43.0%) |

The drop from the raw EC intersection to the scored “on map” rate is intentional: many shared EC numbers fail the metabolite Jaccard floor (`min_score=0.6`) and stay **related** (iIT341 153 without vicinity, 126 with). Vicinity recovers further nodes from unambiguous EC anchors (iIT341 +27 reactions).

Hand check (e_coli_core): PGI, PFK, GAPD, CS, ENO, FBA, ICDHyr land on the expected EC nodes with scores ≥ 0.75; exchanges are `boundary`; only PGL and G6PDH2r remain related.

## Compose coverage (default library includes BRENDA; EC optional)

| Model | Default (EC off) | + EC (`--ec`) | BRENDA structures selected | Generated reactions |
|---|---|---|---:|---:|
| e_coli_core | 11 structures; curated 66; generated 25 | unchanged (0 BRENDA / 0 ec_topology) | 0 | 25 |
| iIT341 | 86 structures; curated 129; generated 274; 1 BRENDA | 98 structures (stages 86+12); curated 130; functional 23; generated 257; 12 BRENDA | 12 | 257 (−17) |

Curated-module reaction coverage: iIT341 264 → 288. No regression on e_coli_core. Importer label→identity stays in `identity_labels`; matching keys/scoring stay in `identifier_bridge` (no duplicated label index).

## What changed in the code

Extended `identifier_bridge.py`, `ec_channel.py`, `brenda_library.py`, `retrieval.py`, `semantic_transfer.py` (no parallel matcher package). Removed the old “exact EC string only” behaviour by folding currentisation and scored name matching into the existing EC channel. Runtime cache for `enzyme.dat` is gitignored under `layout_atlas/runtime/`. BRENDA matching filters `database.modules()` to `source == 'brenda'` so default extension loading does not mix other maps into enzyme-node scoring.

```sh
python3 layout_atlas/map_matching_eval.py
python3 layout_atlas/test_reaction_map_matching.py
```
