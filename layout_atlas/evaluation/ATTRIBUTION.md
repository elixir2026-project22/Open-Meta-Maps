# Attribution

Pathway-specific layout metrics in this package (`metrics.py` overlaps, participant
consistency, compartment metrics, reference similarity; `pathway.py` SBML graph
reading; `escher_bridge.py` Escher map → drawing positions) are adapted from
[layout-benchmark](https://github.com/elixir2026-project22/layout-benchmark)
(`src/layout_benchmark/metrics.py`, `graph.py`, `external.py` `escher_drawing`,
`reference.py`), Copyright 2026 the layout-benchmark copyright holders, licensed
under the MIT License (see that repository's `LICENSE.md`).

Transfer-check field names (`currency_hub_violations`, `currency_long_arc_share`,
`inout_opposite_share`, `edge_len_max_rel`, `bends_per_arc`, `cycle_circularity`,
`cycle_direction`, `cycle_order_consistency`, `flow_main_dir_share`, …) follow
[layout-transfer-check](https://github.com/elixir2026-project22/layout-transfer-check)
and are reimplemented in `map_rules.py` / `layout_quality.py` (no GPL code copied).

`blank_space` and `unrelated_overlaps` are clean-room reimplementations of ideas from
[MetaCarto](https://github.com/forxhunter/MetaCarto) (`src/layout/metrics.py`,
CC BY 4.0). MetaCarto code is not copied; see `docs/layout_evaluation.md`.

gdMetriX is **not** included here. When installed separately it may be imported at
runtime for general graph-drawing metrics; it remains GPL-3.0 and optional.
