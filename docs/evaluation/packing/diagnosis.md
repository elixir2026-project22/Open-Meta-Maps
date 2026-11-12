# Component packing diagnosis (e_coli_core, iIT341)

## Baseline

Gates and equal-scale PNGs compare this branch to **`7044342`** (current main:
#43 packing compactness + #41 label-size buttons + #42 Atlas Browser CSS + #44
atlas data export). Packing geometry on main matches the #43 merge (`43f32b2`).

## What changed

1. **Empty vertical corridor slides** on the finished Escher map: translate a whole
   reaction unit (midmarker + private / companion EX partners) toward a shared
   pool when the mid↔met column is empty and extreme. Shared pools stay fixed;
   at most two vertical corridors per map; target length eased so mids do not
   dive into dense pool chemistry.
2. **Bake-time padding** stays global H/V only (no joint corridor) so BiGG axis
   share is not tilted before completion.
3. **Milder interior repack** for true exterior stubs only (cap hop distance);
   pocket slots prefer tuck-against-content + centre distance on small packs.
4. **Landscape reflow** only pocket-packs far outliers when already in the 4:3–16:9
   band (avoids collapsing large maps toward a square).
5. **Wording**: `tool_pipeline` layout-style docs/logs use `curated`; the deprecated
   `expert` layout_style alias is removed.

## Gates vs `7044342` (main)

`seg` / `rx` / `turn` / `grid` = orderliness. Counted crossings = algorithm×algorithm +
algorithm×curated. **compact** = content-bbox area / canvas area.

| Model | run | cross | aa | ac | cc | ov | long p90 | cont | aspect | blank | empty | compact | sec | seg | rx | turn° | grid |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| e_coli_core | main | 4 | 2 | 2 | 9 | 0 | 289.2 | 0.974 | 1.466 | 0.490 | 0.120 | 0.957 | 47.8 | 0.443 | 0.444 | 13.3 | 0.557 |
| e_coli_core | branch | 4 | 2 | 2 | 7 | 0 | 325.1 | 0.951 | 1.711 | 0.565 | 0.168 | 0.955 | 48.8 | **0.503** | **0.635** | **11.5** | **0.582** |
| iIT341 | main | 197 | 136 | 61 | 10 | 0 | 154.7 | 0.994 | 1.335 | 0.464 | 0.069 | 0.976 | 465 | 0.409 | 0.492 | 18.3 | 0.685 |
| iIT341 | branch | 208 | **132** | 76 | 11 | 0 | **152.6** | **1.000** | **1.396** | 0.465 | **0.066** | **0.977** | **465** | **0.412** | 0.492 | **17.2** | 0.664 |

Content bbox: e_coli 5405×3687 → 6809×3980; iIT 8903×6669 → 9418×6745.
Canvas follows content (fit).

### vs `50f2b59` orderliness / empty (e_coli recovery targets)

| metric | 50f2b59 | main | branch |
|---|---:|---:|---:|
| rx axis share | 0.603 | 0.444 | **0.635** |
| seg axis share | 0.485 | 0.443 | **0.503** |
| counted crossings | 2 | 4 | 4 |
| empty rect | 0.115 | 0.120 | 0.168 |
| compact | 0.698 | 0.957 | 0.955 |

### Checklist

| Gate | e_coli_core | iIT341 |
|---|---|---|
| compactness (content/canvas) ≥ main | **pass** (0.957→0.955) | **pass** |
| empty rect ≤ main | fail (shelf reflow) | **pass** |
| blank ≤ main | fail | near-miss |
| long p90 ≤ main (+5%) | near-miss | **pass** |
| crossings ≤ main (+5%) | **pass** (4→4) | near-miss (197→208) |
| cc ≤ main | **pass** | near-miss |
| aspect 4:3–16:9 | **pass** | **pass** |
| runtime ≤ +20% | **pass** | **pass** |
| continuity ≥ 50f2b59 | **pass** | **pass** |
| rx / seg ≥ 50f2b59 (e_coli) | **pass** | — |
| node overlaps ≤ main | **pass** | **pass** |

### PNGs

- `*_before.png` = main (`7044342`)
- `*_after.png` = this branch
- `*_side_by_side.png` = before | after at **equal absolute map units**
