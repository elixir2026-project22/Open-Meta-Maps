# Reaction–map matching evaluation

Open-data matching of GEM reactions to Atlas nodes (especially BRENDA name-only metabolites). Method ideas from [ecmapper](https://github.com/elixir2026-project22/ecmapper) (GPL-3.0; independent rewrite — no code copied). Thresholds: `identifier_bridge.THRESHOLDS` / [methodology.md](methodology.md).

## BRENDA node match rates (non-boundary reactions)

| Model | Before (exact full EC ∩ BRENDA, no score) | After (scored + current EC + names, no vicinity) | After + vicinity `propagate` |
|---|---:|---:|---:|
| e_coli_core | 48 / 75 (64.0%) | 49 / 75 (65.3%) | 49 / 75 (65.3%) |
| iIT341 | 315 / 477 (66.0%) | 178 / 477 (37.3%) | 205 / 477 (43.0%) |
| iMM904 | 686 / 1413 (48.5%) | 307 / 1413 (21.7%) | 376 / 1413 (26.6%) |

The drop from the raw EC intersection to the scored “on map” rate is intentional: many shared EC numbers fail the metabolite Jaccard floor (`min_score=0.6`) and stay **related** (iIT341 126, iMM904 375). Vicinity recovers further nodes from unambiguous EC anchors (iIT341 +27 reactions, iMM904 +69).

Hand check (e_coli_core): PGI, PFK, GAPD, CS, ENO, FBA, ICDHyr land on the expected EC nodes with scores ≥ 0.75; exchanges are `boundary`; only PGL and G6PDH2r remain related.

## Compose coverage (core Atlas vs core + EC + BRENDA)

| Model | Core only | + EC + BRENDA | BRENDA structures selected | Generated reactions |
|---|---|---|---:|---:|
| e_coli_core | 11 structures; expert_derived 66; generated 25 | same (core already covers the enzymatic set) | 0 | 25 |
| iIT341 | 88 structures; expert 134; generated 280 | 100 structures (stages 88+1+11); expert 135; functional 23; generated 257 | 11 | 257 (−23) |
| iMM904 | 118 structures; expert 226; generated 1195 | 153 structures (118+10+25); expert 235; functional 61; generated 1128 | 25 | 1128 (−67) |

Expert-module reaction coverage: iIT341 262 → 286; iMM904 364 → 434. No regression on e_coli_core.

## What changed in the code

Extended `identifier_bridge.py`, `ec_channel.py`, `brenda_library.py`, `retrieval.py`, `semantic_transfer.py` (no parallel matcher package). Removed the old “exact EC string only” behaviour by folding currentisation and scored name matching into the existing EC channel. Runtime cache for `enzyme.dat` is gitignored under `layout_atlas/runtime/`.

```sh
python3 layout_atlas/map_matching_eval.py
python3 layout_atlas/test_reaction_map_matching.py
```
