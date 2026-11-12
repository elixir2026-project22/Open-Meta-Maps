# Component packing diagnosis (e_coli_core, iIT341)

## Baseline

Gates and PNGs below compare this branch to **`54342af`** (layout baseline; main is
now **`7bc5c22`** after #35 docs + #36 Atlas Browser — no layout algorithm change vs
54342af). Absolute same-map rigid reuse at 1:1 plus padding compression is preserved.

## Padding compression

Curated BiGG excerpts often stretch straight H/V runs only to align columns/rows with the
canvas. On same-map rigid bake we:

1. Take absolute midmarker/metabolite coords for the shared-source group.
2. Collapse empty axis-aligned gaps longer than `1.25 × median mid↔met segment length`
   (monotone remap — order, direction, sibling alignment kept; no rotation/mirroring).
3. Place the compressed cloud at scale 1:1 (translation only) and snapshot for restore.
4. Replay the same remap when reusing source-path control points.

### e_coli_core (8 blocks / Core metabolism)

| Axis | Gaps compressed | Length removed |
|---|---:|---:|
| X | 3 | 164.93 px |
| Y | 2 | 448.86 px |
| **Total** | **5 segments** | **613.79 px** |

`max_gap` = 241.14 px (median BiGG mid↔met ≈ 192.9 × 1.25).

iIT341 has no dominant same-map rigid group → padding compress does not fire.

## Gates vs `54342af`

`seg` / `rx` / `turn` / `grid` = orderliness (main-segment axis share, reaction axis share,
chain-turn median°, metabolite grid share). Counted crossings = algorithm×algorithm +
algorithm×curated.

| Model | run | cross | aa | ac | cc | ov | long p90 | cont | aspect | blank | empty | sec | seg | rx | turn° | grid |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| e_coli_core | 54342af | 18 | 6 | 12 | 10 | 0 | 118.0 | 0.845 | 1.147 | 0.605 | 0.144 | 88.3 | 0.154 | 0.206 | 12.4 | 0.169 |
| e_coli_core | branch | **7** | **2** | **5** | 19 | **0** | 348.5 | **0.921** | **1.371** | **0.573** | **0.128** | **55.2** | **0.485** | **0.619** | **11.6** | **0.557** |
| iIT341 | 54342af | 238 | 139 | 99 | 14 | 2 | 130.6 | 0.938 | 0.808 | 0.860 | 0.221 | 417 | 0.257 | 0.218 | 20.2 | 0.611 |
| iIT341 | branch | **229** | 148 | **81** | **12** | **0** | 154.6 | **0.994** | **1.335** | **0.623** | 0.234 | **485** | **0.414** | **0.477** | **18.4** | **0.628** |

Canvas: e_coli 54342af ~2616×2268 → branch ~6798×5809 (rigid BiGG excerpt + algo tiles).
iIT 54342af landscape after packing; branch aspect 1.335 after partial X squeeze.

### Checklist

| Gate | e_coli_core | iIT341 |
|---|---|---|
| orderliness (seg/rx/grid) vs 54342af | **pass** (large lift) | **pass** (lift) |
| crossings ≤ +5% | **pass** (−61%) | **pass** (−3.8%) |
| blank ≤ 54342af | **pass** | **pass** |
| empty rect ≤ 54342af | **pass** | fail (+5.9%; partial squeeze keeps 4:3 band) |
| aspect 4:3–16:9 | **pass** | **pass** |
| runtime ≤ +20% | **pass** (−37%) | **pass** (+16%) |
| continuity ≥ 54342af | **pass** | **pass** |
| node overlaps ≤ 54342af | **pass** | **pass** |

### Fixes in this iteration

- Spatial continuity pieces after absolute restore (3× median NN).
- KD-tree escher overlap clear; peri gate skipped for `n_rx ≥ 200`.
- Hotspot after restore/squeeze; CrossingIndex uses evaluate skip set and rebuilds
  after every accept (no phantom deltas); larger try/time budgets.
- Small-map post-restore `edge_style` restyle.
- Partial empty-band squeeze: when full collapse would leave 4:3–16:9, collapse only
  to the nearer band edge.

### PNGs

- `*_before.png` = `54342af` baseline
- `*_after.png` = this branch
- `*_side_by_side.png` = before | after at **equal absolute map units** (OMM not independently shrunk)
- `e_coli_core_curated_bigg.png` = BiGG Core metabolism (when present)
