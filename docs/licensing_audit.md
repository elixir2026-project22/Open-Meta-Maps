# Licensing audit

Status: audit of what the repository currently redistributes, dated 2026-11-11. Evidence is either local (term snapshots and provenance files in `data/source_inputs/`) or a public web page checked on that date. **This is an engineering review, not legal advice.** Items marked *open* need a decision or confirmation before a public release.

Policy: the project redistributes **only** openly redistributable maps (CC0 / CC BY by default; MIT/permissive when it clearly covers the data), all included by default. No NC or GPL map packages remain. Module counts below match [`data/atlas_registry.json`](../data/atlas_registry.json) (default load = open-core SQLite + default extensions).

## 1. Atlas content (`open_core.sqlite` + default extensions)

| Source | Modules | Licence | Evidence | Result |
|---|---|---|---|---|
| WikiPathways (semantic in SQLite) | 1,245 | CC0 | `source_inputs/terms/wikipathways.html` (CC0 waiver for content); snapshot **20261110** | OK |
| WikiPathways (geometry-only extension) | 29,673 | CC0 | same; stored under `data/extensions/wikipathways_geometry/` for size | OK |
| Reactome | 3,405 | CC0 for data and derived files; CC BY 4.0 for pathway illustrations and artwork | `source_inputs/terms/reactome.html` | OK for layout data; if Reactome illustrations are ever added, CC BY 4.0 attribution is required |
| Escher maps | 13,047 | CC0-1.0 | Escher community-maps repository `LICENSE` is CC0 | OK |
| BRENDA pathway maps | 11,703 | CC BY 4.0 | https://www.brenda.dsmz.de/license.php; doi:10.1093/nar/gkaf1113 | OK (layout possibly automatic; retrieval deprioritises vs curated community maps) |
| Human-maps | 8,540 | CC BY 4.0 (JSON) | official `escher/` JSON only | OK — companion SVGs are CC BY-SA and are **not** imported |
| Metabolism Regulation gallery | 423 | CC BY 4.0 | https://metabolismregulation.github.io/about/; SBGN-ML under `source_inputs/raw/metabolismregulation/` | OK (hand-drawn; attribution required) |
| PANTHER Pathway | 238 | CC BY 4.0 | `data.pantherdb.org` CD4.1 CellDesigner; `source_inputs/raw/panther/` | OK (labels bridged via BiGG/ChEBI identity layer) |
| iCH360 Visualisation | 503 | MIT | Repo-root `LICENSE` (MIT) covers `Visualisation/` (no carve-out; `Visualisation/readme.md` has no restrictive notes). Paper Data Availability points at GitHub. Zenodo deposit doi:10.5281/zenodo.15048160 of the same Visualisation tree is CC BY 4.0 (corroboration). Local: `data/terms/ich360_LICENSE.txt`, `data/terms/ich360_zenodo_15048160.html`. Stage B keeps size ≥4 only. | OK (permissive; attribution / copyright notice retained) |
| **Total (default load)** | **68,777** | | | |

### Escher maps

The 22 maps come from the Escher community-maps repository (BiGG-based maps for e_coli_core, iAF1260, iJO1366, iMM904, RECON1). The upstream repository LICENSE is CC0; Atlas modules record `CC0-1.0`. The Atlas source label is `escher`.

## 2. Third-party files tracked in the repository

| Path | Licence | Result |
|---|---|---|
| `data/rhea/` | CC BY 4.0 (`rhea/LICENSE.txt`: copyright statement must be reproduced with each copy) | OK, keep `LICENSE.txt`; credit Rhea in the README |
| `data/source_inputs/archives/wp-*.zip` and `raw/wp/` | CC0 (WikiPathways); snapshot **20261110** | OK |
| `data/reactome/` | CC0 data (Reactome v97) | OK |
| `data/models/AraCore_v2_1.xml` | MIT (`LICENSE.AraCore`, Philipp Wendering) | OK |
| `data/models/iYO844.xml` | from BiGG (BiGG License) | **Open**: non-commercial/educational terms |
| `examples/` (e_coli_core, iML1515, iMM904, iIT341) | BiGG models (BiGG License) | **Open**: same as above; iIT341 source not recorded |
| `data/identifier_bridge.json` | annotations extracted from BiGG models and Human-GEM (CC-BY-4.0 according to its GitHub page); the Human-GEM file itself is no longer in the repository | cite BiGG and Human-GEM; terms of the source models apply |
| `data/escher_schema/map-schema-1-0-0.json` | Escher, MIT (`LICENSE` copied) | OK |
| `data/routing/source.zip` (adaptagrams, libavoid) | adaptagrams is LGPL; exact terms not re-checked | **Open**: confirm and add the LGPL notice |
| Human-maps (official `escher/` JSON) | CC BY 4.0 | **Included** (8,540 modules via `database.escher_diagram`). Companion SVG exports are **CC BY-SA** and are **not** imported |
| Metabolism Regulation gallery SBGN-ML | CC BY 4.0 | **Included** via `layout_atlas/sbgn.py` (`data/extensions/metabolismregulation/`). See [metabolismregulation_import.md](metabolismregulation_import.md). |
| PANTHER Pathway CellDesigner | CC BY 4.0 | **Included** (8 metabolic CD4.1 maps → 238 modules) via `celldesigner.py` + `identity_labels.py`. See [panther_pilot.md](panther_pilot.md). |
| iCH360 Visualisation Escher JSON | MIT | **Included** (10 maps → 503 Stage B modules, size ≥4; `compressed_map.json` excluded). See [ich360_import.md](ich360_import.md). |
| `data/extensions/{units,source_paths,reaction_geometry_templates}` | derived from the open-core maps | inherits the Escher-map status above |
| `data/identifier_bridge.json` (+ derived Mapper crosswalk) | BiGG↔ChEBI / EC from public model annotations; Mapper view built at load time | inherits the model licences above |

Excluded (not redistributed): Yeast-maps / original yeast map (GPL-3.0); MicroMap (NC-SA). `celldesigner.py` remains a generic CellDesigner importer.

## 3. Default extensions and caches (not optional packages)

- BRENDA pathway modules: CC BY 4.0; included by default under `data/extensions/brenda/` with per-module attribution to BRENDA + NAR 2026 doi:10.1093/nar/gkaf1113; layout provenance marked possibly automatic.
- PANTHER Pathway (CC BY 4.0): included from `data.pantherdb.org` FTP CD4.1. See [panther_pilot.md](panther_pilot.md).
- iCH360 Visualisation (MIT): included from pinned GitHub commit under `data/source_inputs/raw/ich360/`. See [ich360_import.md](ich360_import.md).
- `data/references/`, `Ref/`, `Presentation/`, `results/`: not tracked.

## 4. Open actions

1. Decide how to handle the BiGG-licensed example models: obtain written confirmation, replace them with permissively licensed equivalents, or state the non-commercial restriction in the README. (Escher community-maps layouts themselves are CC0 per upstream LICENSE.)
2. Confirm the exact LGPL terms for adaptagrams/libavoid and keep the notice next to `data/routing/` (MIT for the project code is compatible with linking a locally built LGPL binary).
3. Add Rhea and Human-GEM citations; keep every upstream `LICENSE` next to its data.

Done (kept as current facts, not a backlog): Yeast-maps excluded; Human-maps admitted under CC BY 4.0 via official Escher JSON; former `sbml2escher` dependency removed (`layout_atlas/model_io.py` + `data/currency_metabolites.json`); code is MIT and Atlas data / generated maps mix CC0, CC BY 4.0 and MIT (iCH360) with source attribution; PANTHER and iCH360 imported as default extensions.
