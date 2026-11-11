# Should `libs['ec']` be always on?

Measured on **e_coli_core**, **iIT341**, **iMM904** (iML1515 skipped), with the improved open-data matching, **EC off vs EC on** (open-core library only; BRENDA not enabled). Defaults were **not** changed.

Reproduce: `python3 layout_atlas/ec_extension_eval.py --summarize results/ec_extension_eval` (after compose jobs in that folder).

## Numbers

| Model | Setting | Structures | Expert-module rxns | exact / adapted / topological / functional | EC proposals kept / rejected | Crossings | Label overlaps | Node overlaps | Length p90 (px) | Runtime (s) |
|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|
| e_coli_core | off | 11 | 68 | 66 / 0 / 2 / 0 | — | 74 | 6 | 1 | 129.7 | 27.3 |
| e_coli_core | on | 11 | 68 | 66 / 0 / 2 / 0 | 198 / 80 | 74 | 6 | 1 | 129.7 | 29.0 |
| iIT341 | off | 88 | 262 | 134 / 1 / 127 / 0 | — | 455 | 14 | 1 | 133.8 | 109.8 |
| iIT341 | on | 89 | 264 | 135 / 1 / 127 / 1 | 311 / 169 | **475** | **19** | **4** | **139.2** | 104.5 |
| iMM904 | off | 118 | 364 | 226 / 0 / 138 / 0 | — | 2684 | 89 | 14 | 164.1 | 501.8 |
| iMM904 | on | 128 | 384 | 235 / 0 / 138 / 11 | 481 / 397 | 2583 | 90 | 15 | 159.9 | 503.4 |

Notes:

- `exact` / `adapted` / `topological` / `functional` are `evidence.layout_provenance` counts (`expert_exact`, `expert_adapted`, `expert_topological`, `expert_functional`).
- EC proposals are after `ec_provenance.screen` (incompatible dropped). Many kept proposals are `insufficient_evidence` and may still be selected as **functional** layout suggestions.
- e_coli_core: selection stages `[11, 0]` — EC adds nothing once the strong channel has covered the map.
- iIT341: stages `[88, 1]` — one EC structure (2 reactions). Layout quality slightly **worse**.
- iMM904: stages `[118, 10]` — ten EC structures (~20 reactions in expert modules, 11 still labelled functional). Crossings slightly better; label/node overlaps essentially flat.

## Sample of kept `ec_topology` placements

### iIT341 (1 structure)

| Target | Source | Shared EC | Verdict |
|---|---|---|---|
| FOLD3 (dihydropteroate synthase) | DHPS2 | 2.5.1.15 | **Correct** (same reaction) |
| DHFS | DHFS | 6.3.2.12 / 6.3.2.17 | **Correct** (equation-verified → expert_exact) |

### iMM904 (10 structures) — notable cases

| Targets | Shared EC | Verdict |
|---|---|---|
| HMGCOAR + MEVK1 | 1.1.1.34, 2.7.1.36 | OK (mevalonate path) |
| ASNN + ASNS1 | 3.5.1.1 / 3.5.1.38, 6.3.5.4 | OK (equation-verified) |
| MTHFC + MTHFD2i; MTHFCm + MTHFDm | folate ECs | OK (cytosolic / mitochondrial folate cycle) |
| PPM + RBK; G6PI3 + G6PDA | pentose / hexose-P | Plausible pathway neighbours |
| IGPS + TRPS1←TRPS3 | 4.1.1.48, 4.2.1.20 | Pathway OK; TRPS1 vs TRPS3 is a variant, not identity |
| **DGGH←MLTG1 + MALT** | **3.2.1.20** | **Wrong / promiscuous EC** (digalactoside hydrolase vs maltodextrin/maltose glucosidase; same broad glycosidase number) |
| **ECOAH11p←ECOAH1m + HACD10p←HACD1m** | **4.2.1.17, 1.1.1.35** | **Doubtful** (C26 peroxisomal fatty-acid steps placed on short-chain mitochondrial drawings) |
| 34HPLFM + TYRTAim | 1.1.1.110 / 1.1.1.237, 2.6.1.5 | Related aromatic metabolism; not equation-verified |

No transport/exchange reactions appeared among selected EC structures in this sample. The failures are **promiscuous ECs** and **chain-length / isoform reuse**, not transporters.

## Recommendation

**Keep the EC extension optional and off by default** (leave the UI checkbox). Do **not** make it always-on yet.

Reasons:

1. **Coverage gain is small on the prokaryote GEMs** (e_coli_core: zero; iIT341: +2 reactions) while screening still costs work and can worsen crossings/overlaps.
2. **Yeast gains coverage** (+20 expert-module reactions) but still admits clear false/doubtful functional placements when equations are not verified.
3. The channel already warns that EC is weak identity; turning it on silently would present functional suggestions as part of the default map without the user opting in.

### Fix before considering default-on

1. **Gate selection**: only admit EC proposals with `verified_equivalent` (or a high metabolite Jaccard on name-only sources) into the skeleton; keep `insufficient_evidence` out of default selection (or behind an explicit “allow weak EC” flag).
2. **Downrank promiscuous ECs** when unverified (at least glycoside hydrolase 3.2.1.20-class and generic acyl-CoA hydratase/dehydrogenase families used across chain lengths).
3. **Re-measure** the three models after (1)–(2): coverage, layout quality, and a fresh false-positive sample must not regress vs EC-off.

BRENDA-assisted EC (`--ec --brenda`) is a separate opt-in: it adds more functional coverage (see `docs/reaction_map_matching.md`) but carries the same unverified-equation risk and an unverified layout-origin warning.
