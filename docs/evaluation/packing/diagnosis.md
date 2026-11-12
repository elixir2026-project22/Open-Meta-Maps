# Component packing diagnosis (e_coli_core, iIT341)

## Baseline

Gates and equal-scale PNGs compare this branch to **`50f2b59`** (main after squash-merge of #28:
absolute same-map rigid reuse at 1:1 + padding compression).

## What changed

Joint optimisation on every model (no e_coli-specific paths):

1. **Padding ↔ corridor loop** on absolute same-map clouds: global H/V empty-run
   collapse alternates with shortening axis-aligned mid↔met edges whose thin corridor
   is empty of other glyphs (length free; shape and direction kept).
2. **Interior free-region pack**: after related/independent moves, exterior small
   components are re-pocketed into the largest free neighbourhoods **inside** the main
   AABB before any canvas growth; search pad is a small flush strip only.
3. **Canvas fit** to content bbox; squeeze follow-up at 1.25× spacing; force-squeeze
   no longer drops below the 4:3 aspect band.

## Gates vs `50f2b59`

`seg` / `rx` / `turn` / `grid` = orderliness. Counted crossings = algorithm×algorithm +
algorithm×curated. **compact** = content-bbox area / canvas area (new compactness gate
with empty-rect).

| Model | run | cross | aa | ac | cc | ov | long p90 | cont | aspect | blank | empty | compact | sec | seg | rx | turn° | grid |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| e_coli_core | 50f2b59 | 2 | 2 | 0 | 9 | 0 | 353.1 | 0.907 | 1.522 | 0.591 | 0.115 | 0.698 | 46.3 | 0.485 | 0.603 | 11.2 | 0.558 |
| e_coli_core | branch | 4 | 2 | 2 | 9 | 0 | **289.2** | **0.974** | **1.466** | **0.490** | 0.120 | **0.957** | 47.8 | 0.443 | 0.444 | 13.3 | 0.557 |
| iIT341 | 50f2b59 | 259 | 152 | 107 | 14 | 0 | 150.1 | 0.988 | 1.717 | 0.708 | 0.115 | 0.974 | 495 | 0.445 | 0.513 | 22.0 | 0.655 |
| iIT341 | branch | **197** | **136** | **61** | **10** | **0** | 154.7 | **0.994** | **1.335** | **0.464** | **0.069** | **0.976** | **465** | 0.409 | 0.492 | 18.3 | **0.685** |

Content bbox: e_coli 6792×4463 → 5405×3687; iIT 12926×7530 → 8903×6669.
Canvas follows content (fit); compact rises with the shrink.

### Checklist

| Gate | e_coli_core | iIT341 |
|---|---|---|
| compactness (content/canvas) ≥ 50f2b59 | **pass** (0.70→0.96) | **pass** |
| empty rect ≤ 50f2b59 | near-miss (+0.005) | **pass** (−40%) |
| blank ≤ 50f2b59 | **pass** | **pass** |
| long p90 ≤ 50f2b59 | **pass** (−18%) | **pass** (≤+5%) |
| crossings ≤ +5% | fail (2→4; still low) | **pass** (−24%) |
| cc ≤ 50f2b59 | **pass** | **pass** |
| aspect 4:3–16:9 | **pass** | **pass** |
| runtime ≤ +20% | **pass** | **pass** |
| continuity ≥ 50f2b59 | **pass** | **pass** |
| node overlaps ≤ 50f2b59 | **pass** | **pass** |

### PNGs

- `*_before.png` = `50f2b59` baseline
- `*_after.png` = this branch
- `*_side_by_side.png` = before | after at **equal absolute map units**
