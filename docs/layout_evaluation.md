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

| metric | e_coli_core | iIT341 |
|---|---:|---:|
| runtime_s | 66.9 | 121.0 |
| geom_crossings (algo×algo + algo×curated) | 47 | 337 |
| … main-metabolite / currency-involved | 5 / 42 | 194 / 143 |
| overlapping_pairs (geom nodes+labels) | 7 | 22 |
| long lines p90 (reaction spacings) | see report | see report |
| edge_len_max_rel | 15.91 | 15.86 |
| inout_opposite_share (≥120°, all / core) | 0.813 / 0.794 | 0.820 / 0.691 |
| participant_consistency (>90°, CLI) | 0.422 | 0.434 |
| currency_hub_violations | 0 | 3 |
| currency_long_arc_share | 0.000 | 0.002 |
| structure_continuity | 0.856 | 0.947 |
| split_structures | 0 | 3 |
| cycle_circularity | 0.289 | 0.377 |
| cycle_direction_consistency | 0.575 | 0.600 |
| cycle_order_consistency | 1.000 | 1.000 |
| flow_main_dir_share | 0.302 | 0.273 |
| bends_per_arc | 0.089 | 0.172 |
| exchange_periphery_share | 0.417 | 0.633 |
| blank_space_share | 0.546 | 0.842 |
| space_utilisation (= 1 − blank) | 0.454 | 0.158 |
| largest_empty_rect_share | 0.240 | 0.260 |
| unrelated_overlaps_total | 4 | 17 |
| aspect_ratio (W/H, content bbox) | 0.677 | 0.477 |
| species_in_other_compartment | 0.583 | 0.862 |
| overlapping_pairs (benchmark boxes) | 198 | 1381 |
| overlap_area_share | 0.322 | 0.296 |
| long_edge_share (>3× median) | 0.211 | 0.232 |

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

Only metrics the layout optimises or should optimise. Each row has a plain-language label
(EN/ZH), a one-sentence tip (higher/lower + target), and the responsible pipeline stage
(or **not yet optimised** / **待优化**). Values that clearly miss the target use a subtle
amber highlight.

| Group | Rows | Stage |
|---|---|---|
| Connections | Local crossing hotspot; Unrelated overlaps | Bézier routing · quality gate; spacing / relaxation · Bézier routing |
| Line shape | Detour; Hooks / sharp turns / self-intersections; Bends per connection | Bézier routing |
| Structure & flow | Exchange reactions on the periphery; Cycle circularity / direction / order; Species in another compartment | composition & peripherals; **not yet optimised** (ring/cycle placement for algorithm-placed cycles; compartment packing) |
| Space & shape | Map size; Aspect ratio vs target; Blank space / space utilisation; Largest empty rectangle; Isolated reactions; Most reactions in one cell | composition & packing |
| Reproducibility | Procrustes / mean displacement (when `evaluate --compare-map`) | **not yet optimised** (detail only, not a main card) |
| Debug diagnostics (hidden by default) | Main flow-direction share | debug only — not an optimisation target (would fight curated orientations) |
| Pipeline diagnostics (collapsed) | Composition-pass counts; Peripheral-module counts | composition & peripherals |

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
| iIT341 | 0.947 | (see report) | 3 | (see report) |

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
