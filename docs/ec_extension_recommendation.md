# Should `libs['ec']` be always on?

Measured on **e_coli_core** and **iIT341** after the matching improvements and rebase onto open-licence Atlas + default BRENDA extension (iMM904 / iML1515 skipped for development). **EC off vs EC on**; defaults were **not** changed.

## Summary table

| Model | EC | Modules | Curated-module rxns | layout_prov exact / adapted / topo / func | EC kept / rejected | Crossings | Labels | Node ov. | p90 (px) | Runtime (s) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| e_coli_core | off | 11 | 68 | 66 / 0 / 2 / 0 | — | 75 | 7 | 0 | 121.9 | 36.4 |
| e_coli_core | on | 11 | 68 | 66 / 0 / 2 / 0 | 301 / 80 | 75 | 7 | 0 | 121.9 | 50.9 |
| iIT341 | off | 86 | 264 | 129 / 0 / 135 / 0 | — | 512 | 21 | 3 | 142.8 | 118.1 |
| iIT341 | on | 98 | 288 | 130 / 0 / 135 / 23 | 465 / 169 | 354 | 17 | 2 | 133.6 | 131.5 |

## Coverage and quality

- e_coli_core: EC on changes nothing placed; screening still runs (~15 s extra).
- iIT341: +12 structures (all `ec_topology`, mostly BRENDA), +24 curated-module reactions, −17 generated; crossings/labels/p90 improve. Functional provenance: 23 reactions.

## Sample of kept `ec_topology` placements (iIT341)

24 reactions with `match_level=ec_topology`. Notable:

| Reaction | Shared EC | Chemical | Source | Note |
|---|---|---|---|---|
| DHFS | 6.3.2.12 / 6.3.2.17 | compatible (`exact_equation`) | community | verified_equivalent — good |
| FOLD3 | 2.5.1.15 | unknown | community:DHPS2 | `ec_only` — unchecked equation |
| ADSK, HSK, CDPMEK, … | various | unknown | brenda | BRENDA SVG has no stoichiometry → `insufficient_evidence` / unknown chem |
| AHCYSNS | 3.2.2.9 | unknown | brenda | glycosidase-adjacent; name-only participants |
| No EX_/DM_/transport ids in the sample | | | | |

Most BRENDA-driven EC hits cannot verify the equation, so they remain functional layout suggestions. That is useful for coverage but unsafe as silent default identity.

## Recommendation

**Keep `libs['ec']` optional (UI checkbox); default stays off.**

Before considering default-on:

1. Admit only `verified_equivalent` (or high metabolite Jaccard) into the skeleton by default; keep unverified BRENDA EC as opt-in.
2. Downrank promiscuous / name-only EC matches when chemical status is `unknown`.
3. Re-measure e_coli_core and iIT341 (and later iMM904) after that gate.
