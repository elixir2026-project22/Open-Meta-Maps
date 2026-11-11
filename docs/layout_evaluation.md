# Layout evaluation (development aid)

This is a **development** evaluation of Open Meta Maps layouts against the Project 22
[layout-benchmark](https://github.com/elixir2026-project22/layout-benchmark) rules and the
[layout-transfer-check](https://github.com/elixir2026-project22/layout-transfer-check) field
names. It is **not** a fair comparison with other tools. Atlas data sources are never excluded
because of evaluation overlap.

CLI (optional deps: `pip install -r layout_atlas/requirements-eval.txt`):

```sh
python3 project.py compose --gem examples/e_coli_core.xml --out results/eval_runs/e_coli_core --source-paths
python3 project.py evaluate --gem examples/e_coli_core.xml --map results/eval_runs/e_coli_core/map.json --out results/evaluation/e_coli_core.json
```

One computation path:

1. **`layout_quality.evaluate`** — final-geometry block shared by the quality gate, the tool
   **Layout quality** panel, and `evaluate` (curve crossings with main/currency split, overlaps,
   long lines, participant placement ≥120°, cofactor clones, curated-structure continuity, plus
   a curated **Detailed metrics** table aligned with layout goals).
2. **`layout_atlas/evaluation/`** — pathway metrics aligned with layout-benchmark (MIT; see
   `ATTRIBUTION.md`) plus optional gdMetriX (GPL-3.0, **import only**, never vendored; CLI extras).

User-facing panel strings use **curated** / **人工绘制** (not “expert”). The 120° opposite-side
threshold is a Project 22 first guess (transfer-check). Benchmark `participant_consistency` (>90°)
stays on the evaluate CLI only (not duplicated on the panel).

gdMetriX straight-line `crossings` are skipped when the undirected graph has more than 1200
edges. For maps with more than 800 reactions the curve geometry block is skipped in `evaluate`
and should be read from compose’s `layout_issues.json` / `report.json`.

## Metric table (e_coli_core · iIT341)

Maps composed with `--source-paths` (Bézier / curated layout default). Geometry columns come from
`layout_quality` on the final `map.json`.

Stage A main vs this branch (`--source-paths`):

| metric | e_coli_core main | e_coli_core +iCH360 | iIT341 main | iIT341 +iCH360 |
|---|---:|---:|---:|---:|
| curated_coverage_share | 0.7368 | 0.7368 | 0.3069 | 0.2924 |
| crossings (counted) | 46 | 46 | 299 | 239 |
| overlapping_pairs | 5 | 5 | 18 | 16 |
| unrelated_overlaps_total | 5 | 5 | 17 | 8 |
| blank_space_share | 0.5826 | 0.5826 | 0.9401 | 0.8605 |
| edge_len_max_rel | 12.187 | 12.187 | 27.975 | 27.871 |
| aspect_ratio (W/H) | 1.1726 | 1.1726 | 0.5394 | 0.8079 |
| runtime_s | 77.9 | 74.7 | 457.3 | 454.2 |
| iCH360 modules used | — | 0 | — | 7 |

Refresh notes (2026-11-11): rebased iCH360 onto Stage A (extraction provenance + compactness-aware select) and re-composed e_coli_core / iIT341. e_coli_core unchanged with iCH360 loaded. iIT341 selects 4×`full_map` + 3×`ions_metals_uptake` size-2 geometry windows; crossings/blank/aspect improve vs Stage A main while coverage dips slightly. Older benchmark-aligned rows remain in [`docs/evaluation/summary.csv`](evaluation/summary.csv).

Machine-readable copy: [`docs/evaluation/summary.csv`](evaluation/summary.csv).

**Test scope for this refresh:** only `e_coli_core` and `iIT341` (not iMM904 / iML1515).

## Panel cards

| Card | Definition |
|---|---|
| Connection crossings | algorithm×algorithm + algorithm×curated; tooltip splits main-metabolite vs currency-involved; curated×curated not counted |
| Overlaps | node + label overlap pairs; curves through nodes/labels in the subtitle |
| Long lines | length p90 in median reaction spacings; `edge_len_max_rel` = max/median |
| Participant placement | share of reactant/product pairs ≥120° apart (all and core-only) |
| Cofactor clones | `currency_hub_violations` (currency instance linking reactions farther than 4× median spacing); `currency_long_arc_share` (currency lines >3× spacing) |
| Curated-structure continuity | contiguous pair share of reused curated structures (3× median reaction spacing) |

Bézier coverage is a **warning** when below 100% (still a quality-gate hard check), not a card.

## Detailed metrics (panel)

Only metrics the layout optimises or should optimise (layout **targets**, including backlog).
Non-target diagnostics are hidden under Debug. Each row has a plain-language English label,
a one-sentence tip (higher/lower + target), and the responsible pipeline stage
(or **not yet optimised** / **待优化**). Values that clearly miss the target use a subtle
amber highlight.

| Group | Rows | Stage |
|---|---|---|
| Connections | Local crossing hotspot; Unrelated overlaps | Bézier routing · quality gate; spacing / relaxation · Bézier routing |
| Line shape | Detour; Hooks / sharp turns / self-intersections; Bends per connection | Bézier routing |
| Structure & flow | Exchange reactions on the periphery; Axis-aligned main segments; Axis-aligned reaction axes; Grid-aligned metabolites; Median chain turn angle; Cycle circularity / direction / order; Species in another compartment | composition & peripherals; **not yet optimised** (orderliness / octilinear snap; ring/cycle placement; compartment packing) |
| Space & shape | Map size; Aspect ratio vs target; Blank space / space utilisation; Largest empty rectangle; Isolated reactions; Most reactions in one cell | composition & packing |
| Reproducibility | Procrustes / mean displacement (when `evaluate --compare-map`) | **not yet optimised** (detail only, not a main card) |
| Debug diagnostics (hidden by default) | Main flow-direction share | debug only — not an optimisation target (would fight curated orientations) |
| Pipeline diagnostics (collapsed) | Composition-pass counts; Peripheral-module counts | composition & peripherals |

**Orderliness (待优化):** main-segment and reaction-axis shares within ±10° of horizontal/vertical
(cofactors excluded); metabolite grid share (tolerance 0.05 × median main-edge length); median
turn angle at degree-2 main metabolites along linear chains.

Not added from layout-transfer-check: raw density, straight-line crossing counts, catalyst-beside, arrow-weight.

**One bounding box:** map size, aspect ratio (W/H) and blank / utilisation all use the
content-node bbox (metabolites + midmarkers). `space_utilisation = 1 − blank_space_share`
on the same blank-space raster.

**Removed from the panel** (duplicates or pure statistics with no layout goal): connection
length percentiles / absolute px thresholds; `edge_len_max_rel` (already on Long lines);
`participant_consistency` >90° (duplicate of the Participant placement card); `area_per_node`;
subsystem conventions; 5th-neighbour distance percentiles; the second, inconsistent aspect /
space-utilisation row.

## Structure continuity

| model | contiguous_pair_share | structures | split_structures | pairs |
|---|---:|---:|---:|---:|
| e_coli_core | 0.856 | 11 | 0 | 146 |
| iIT341 | 0.950 | 90 | 3 | 238 |

## Reproducibility / stability

`evaluate --compare-map` compares two maps with Procrustes residual displacement. Compose does
not yet expose a seed flag.

## blank_space and unrelated_overlaps (MetaCarto ideas)

Clean-room reimplementations of ideas from
[MetaCarto](https://github.com/forxhunter/MetaCarto) `src/layout/metrics.py` (CC BY 4.0).
No MetaCarto code is copied. Single implementation: `evaluation.metrics.blank_space`
(also exposed as `blank_space_from_positions` for component packing).

- **blank_space** — rasterise nodes, segments and labels at cell size ≈ 0.6 × median
  reaction spacing; a cell is blank if no ink lies within ≈ 1 × spacing. Report
  `blank_space_share`, `space_utilisation = 1 − blank`, and `largest_empty_rect_share`.
- **unrelated_overlaps** — node-on-node, edge-through-node, and collinear edge-on-edge
  overlaps between reactions that share no non-currency metabolite.

## Licence note

Pathway-specific metric ports: MIT (layout-benchmark). MetaCarto metric *ideas*
credited above (CC BY 4.0). gdMetriX remains optional GPL-3.0 — imported when
installed, never vendored into this repository.
