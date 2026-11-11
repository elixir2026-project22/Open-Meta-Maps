# Layout evaluation (development aid)

This is a **development** evaluation of Open Meta Maps layouts against the Project 22
[layout-benchmark](https://github.com/elixir2026-project22/layout-benchmark) rules. It is
**not** a fair comparison with other tools. Atlas data sources are never excluded because of
evaluation overlap.

CLI (optional deps: `pip install -r layout_atlas/requirements-eval.txt`):

```sh
python3 project.py compose --gem examples/e_coli_core.xml --out results/eval_runs/e_coli_core --source-paths --layout-style expert
python3 project.py evaluate --gem examples/e_coli_core.xml --map results/eval_runs/e_coli_core/map.json --out results/evaluation/e_coli_core.json
```

The evaluator merges:

1. **Benchmark-aligned metrics** (`layout_atlas/evaluation/`): pathway metrics ported from
   layout-benchmark (MIT; see `ATTRIBUTION.md`) plus optional gdMetriX (GPL-3.0, import only).
2. **Final-geometry** from the existing `layout_quality.evaluate` (curve crossings, label/node
   collisions) — the same block the quality gate uses — so there is one scoring path for
   development, without a second parallel quality module.

gdMetriX straight-line `crossings` are skipped when the undirected graph has more than 1200
edges (quadratic cost). For maps with more than 800 reactions the curve geometry block is
skipped in `evaluate` and should be read from compose’s `layout_issues.json` / `report.json`.

## Metric table

Maps composed with `--source-paths --layout-style expert` (2026-11-11). `geom_*` are from
`layout_quality` (Bezier-aware). Straight-line gdMetriX crossings are omitted for larger GEMs.

| model | runtime_s | participant_consistency | overlapping_pairs | overlap_area_share | long_edge_share | species_in_other_compartment | aspect_ratio | geom_crossings | geom_label_overlaps | stress | procrustes | knn_preservation |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| e_coli_core | 61 | 0.401 | 191 | 0.285 | 0.194 | 0.583 | 0.726 | 76 | 28 | 5.3e3 | — | — |
| iIT341 | 239 | 0.434 | 1249 | 0.278 | 0.288 | 0.847 | 0.131 | 547 | 52 | 2.6e5 | — | — |
| iMM904 | 640 | 0.444 | 3826 | 0.319 | 0.409 | 0.965 | 0.520 | 2322† | 300† | 1.9e6 | — | — |
| e_coli_iCH360 | 103 | 0.387 | 826 | 0.302 | 0.126 | 0.520 | 0.833 | 384 | 48 | 6.7e4 | — | — |
| Mycoplasma_iPS189 | 103 | 0.431 | 1183 | 0.508 | 0.177 | 0.971 | 0.669 | 231 | 27 | 1.5e5 | — | — |
| R-HSA-77289 | 20 | 0.484 | 2 | 0.001 | 0.341 | 0.025 | 0.602 | 137 | 1 | 2.7e3 | — | — |
| R-HSA-70171 | 12 | 0.536 | 0 | 0.000 | 0.275 | 0.000 | 0.291 | 15 | 2 | 5.5e2 | — | — |
| F001-glycolysis | 8.5 | 0.667 | 0 | 0.000 | 0.000 | — | 0.561 | 1 | 0 | 63 | 0.305 | 0.564 |

† iMM904 curve geometry from compose `report.json` (`layout_quality.final_geometry`), not from
`evaluate` (auto-skipped above 800 reactions).

**data_shared note:** `https://github.com/elixir2026-project22/data_shared` returned 404 here.
Public stand-ins are under `examples/eval_models/` (iCH360, Mycoplasma iPS189, Reactome pathways,
metabolismregulation glycolysis + SBGN reference). Human2 / full Human-GEM was not run (size).
Private BRENDA `brenda_*.xml` files were not available; re-run with a local clone when you have it.
`iML1515` was skipped by design.

Machine-readable copy: [`docs/evaluation/summary.csv`](evaluation/summary.csv).

## Structure continuity

Computed from the final `map.json` for each reused expert structure (≥2 reactions): share of
reaction pairs that share a non-currency metabolite and lie within **3× median reaction spacing**,
plus how many structures split into more than one geometric piece. Scale-free; no model-specific
thresholds. Shown on the tool Layout quality card (EN + zh).

| model | contiguous_pair_share | structures | split_structures | pairs |
|---|---:|---:|---:|---:|
| e_coli_core | 0.856 | 11 | 0 | 146 |
| iIT341 | 0.906 | 85 | 4 | 235 |
| iMM904 | 0.893 | 131 | 7 | 484 |

Expert structures stay mostly contiguous on these maps; residual splits are sparse (≈5% of
structures on the larger GEMs) and are listed under `layout_issues.json` → `split_structures`.

## Where we are weakest (with examples)

Ranked by how consistently they hurt maps across models (not by one model’s absolute scores).

### 1. Participant consistency (reactants / products not on opposite sides)

Benchmark share stays ~0.39–0.44 on GEMs (gallery glycolysis reaches 0.67). Many reactions read
“sideways” or fold substrates and products onto the same half-plane.

- **e_coli_core:** `R_ACKr`, `R_GLCpts`, `R_PIt2r`, `R_TKT2`, `R_ADK1`, `R_FBA` (opposite_share 0).
- **iIT341:** `R_ARGSS`, `R_DBTS`, `R_FTHFLi`, `R_GCALDD` (far modules; e.g. `R_ARGSS` near
  normalised ≈(157, 4667)).
- **e_coli_iCH360:** transporters `R_GLNabcpp`, `R_SULabcpp`, `R_ATPS4rpp`, `R_NADH17pp`.
- **Reactome R-HSA-77289:** `reaction_5690066`, `reaction_77263`, `reaction_77277`.

### 2. Compartment mixing

`species_in_other_compartment` is high on multi-compartment GEMs (0.58 core → 0.85–0.97 on
iIT341 / iMM904 / Mycoplasma). Hulls of different compartments heavily interpenetrate; transport
and periplasmic reactions are not organised as a boundary.

### 3. Long edges and high stress (especially currency / biomass)

Share of edges longer than 3× median length is 0.13–0.41 on GEMs. Worst edges are almost always
currency or biomass spokes.

- **e_coli_core:** `R_BIOMASS_Ecoli_core_w_GAM` ↔ `M_nadp_c` / `M_nadph_c` / `M_nadh_c` / `M_nad_c`
  (≈5.6–6.0× median). Crossing hotspots also list `BIOMASS_Ecoli_core_w_GAM`, `GLUN`, `PPS`.
- **iIT341:** `M_h_c` ↔ `R_SERt2r` / `R_PSD_HP` / `R_PSSA_HP` (≈16×); densest crossing window
  includes `BIOMASS_HP_published` (104 crossings involving that reaction).
- **iCH360:** `M_nadph_c` ↔ `R_EAR181y` / `R_3OAR181`; densest window dominated by `Biomass`
  (142 crossings).
- **Mycoplasma:** `M_pi_c` / `M_atp_c` ↔ `R_Coabc` / `R_SPMDabc` (≈22–25×); hotspot around
  `Biomass`, `PTRC*`, `TYRTRS`.

### 4. Box overlaps and label collisions

Normalised node-box overlapping pairs grow with model size (191 → 1249 → 3826). Escher-native
label overlaps stay high on GEMs (28–52; iMM904 report: 300 label pairs, 307 curves through
nodes). Pathway SBML maps without BiGG identity (Reactome, gallery) are much cleaner on overlaps
— showing the weakness is mainly in **GEM completion density**, not in the metric itself.

### 5. Global aspect / packing (iIT341)

iIT341 aspect ratio ≈0.13 (extremely elongated) with huge stress: modules sit in a long strip /
shelf rather than a compact canvas. Same composition pattern appears as high
`area_per_node` on iMM904.

### 6. Reference similarity (only where a hand layout exists)

F001-glycolysis vs its SBGN reference: Procrustes 0.31, distance Spearman 0.70, kNN preservation
0.56 — usable but not tight; local neighbourhoods are only partly preserved after alignment.

## Ranked, generalizable improvement ideas

No model-specific special cases; no thresholds tuned to one map. Implement on a layout-fixes
branch, not here.

1. **Reaction-axis orientation as a hard local constraint**  
   After placement (and after expert blocks are fixed), enforce opposite half-planes for
   non-currency reactants vs products for every reaction, including transporters. Extend the
   existing orientation pass to multi-substrate reactions and to edges that currently keep
   expert currency copies on the wrong side.

2. **Compartment-aware global packing**  
   Treat compartments as soft regions (or nested group layouts): place all cytosolic cores
   first, then membrane / periplasm / extracellular as bands or rings. Score and gate moves
   with `species_in_other_compartment`, not only crossings.

3. **Currency and biomass as clones / stubs, not long spokes**  
   Always draw high-degree currency metabolites as local clones (or stub glyphs) next to each
   reaction; keep a single optional “hub” legend. Cap biomass fan-out with a dedicated
   radial/stub layout so `BIOMASS_*` stops dominating crossing and length stats.

4. **Crossing-aware joint routing before accepting composition**  
   Run a cheap crossing estimate on algorithm edges when relocating blocks / shelves; reject
   shelves that create dense 400 px windows (as in e_coli_core ≈(250,1050) and iIT341 biomass
   region). Prefer rewiring attachment sides over accepting long detours.

5. **Compact multi-component packing**  
   Replace one-sided shelf stacking with 2D bin packing / squarified layout of independent
   components so aspect ratio stays near 1 and stress drops without model-specific spacing
   constants (scale by median reaction spacing only).

6. **Identity-poor pathway models**  
   For Reactome / non-BiGG SBML, strengthen topology + annotation bridges (ChEBI/EC) so expert
   geometry can attach; today these maps lean on topology-only modules and still show long
   edges between loosely placed islands.

## Reproducibility / stability

`evaluate --compare-map` compares two maps with Procrustes residual displacement (seed /
leave-k workflows). Compose does not yet expose a seed flag; wire that when measuring
reproducibility across seeds on the layout-fixes branch.
