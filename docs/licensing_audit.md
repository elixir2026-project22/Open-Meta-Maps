# Licensing audit

Status: audit of what the repository currently redistributes, dated 2026-11-11. Evidence is either local (term snapshots and provenance files in `data/source_inputs/`) or a public web page checked on that date. **This is an engineering review, not legal advice.** Items marked *open* need a decision or confirmation before a public release.

Policy: the project redistributes **only** openly redistributable maps (CC0 / CC BY), all included by default. No NC or GPL map packages remain.

## 1. Atlas content (`open_core.sqlite` + default extensions)

Core SQLite modules (checked by reading the module payloads), plus default extensions loaded with the Atlas:

| Source | Modules | Licence | Evidence | Result |
|---|---|---|---|---|
| WikiPathways | 4,696 | CC0 | `source_inputs/terms/wikipathways.html` (CC0 waiver for content); snapshot **20261110** | OK |
| Reactome | 3,405 | CC0 for data and derived files; CC BY 4.0 for pathway illustrations and artwork | `source_inputs/terms/reactome.html` | OK for layout data; if Reactome illustrations are ever added, CC BY 4.0 attribution is required |
| Escher maps | 9,098 | CC0-1.0 | Escher community-maps repository `LICENSE` is CC0 | OK |
| BRENDA pathway maps | 11,703 | CC BY 4.0 | https://www.brenda.dsmz.de/license.php; doi:10.1093/nar/gkaf1113 | OK (layout possibly automatic; retrieval deprioritises vs expert maps) |
| Human-maps | 8,540 | CC BY 4.0 (JSON) | official `escher/` JSON only | OK — companion SVGs are CC BY-SA and are **not** imported |
| Metabolism Regulation gallery | 423 | CC BY 4.0 | https://metabolismregulation.github.io/about/; SBGN-ML under `source_inputs/raw/metabolismregulation/` | OK (hand-drawn; attribution required) |
| PANTHER Pathway | 238 | CC BY 4.0 | `data.pantherdb.org` CD4.1 CellDesigner; `source_inputs/raw/panther/` | OK (labels bridged via BiGG/ChEBI identity layer) |

### Escher maps

The 22 maps come from the Escher community-maps repository (BiGG-based maps for e_coli_core, iAF1260, iJO1366, iMM904, RECON1). The upstream repository LICENSE is CC0; Atlas modules record `CC0-1.0`. The Atlas source label is `escher` (renamed from the former `community` label).

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
| ~~`model_input.py`, `central_reference.json` (sbml2escher)~~ | removed 2026-11-10: the SBML reader is now libSBML-based (`layout_atlas/model_io.py`), the currency list is `data/currency_metabolites.json` | **Resolved** |
| `data/routing/source.zip` (adaptagrams, libavoid) | adaptagrams is LGPL; exact terms not re-checked | **Open**: confirm and add the LGPL notice |
| Yeast-maps / original yeast map | GPL-3.0 (Yeast-maps) | **Removed** from packages, acquisition and docs. Not redistributed. |
| Human-maps (official `escher/` JSON) | CC BY 4.0 | **Included** (8,540 modules via `database.escher_diagram`). Companion SVG exports are **CC BY-SA** and are **not** imported; provenance/`SOURCE.json` record the JSON-only CC BY 4.0 scope. |
| Metabolism Regulation gallery SBGN-ML | CC BY 4.0 | **Included** via generic `layout_atlas/sbgn.py` (`data/extensions/metabolismregulation/`). See [metabolismregulation_import.md](metabolismregulation_import.md). |
| PANTHER Pathway CellDesigner | CC BY 4.0 | **Included** (8 metabolic CD4.1 maps → 238 modules) via generic `celldesigner.py` + `identity_labels.py`. See [panther_pilot.md](panther_pilot.md). |
| `data/extensions/{units,source_paths,reaction_geometry_templates}` | derived from the open-core maps | inherits the Escher-map status above |
| `data/identifier_bridge.json` (+ derived Mapper crosswalk) | BiGG↔ChEBI / EC from public model annotations; Mapper view built at load time | inherits the model licences above |

## 3. Not in the repository (verified by `.gitignore` and tracked-file listing)

- BRENDA pathway modules: CC BY 4.0 (https://www.brenda.dsmz.de/license.php); included by default under `data/extensions/brenda/` with per-module attribution to BRENDA + NAR 2026 doi:10.1093/nar/gkaf1113; layout provenance marked possibly automatic. The former optional `brenda` data package entry was removed as redundant.
- PANTHER Pathway (CC BY 4.0): **included** from `data.pantherdb.org` FTP CD4.1 (www HTML still Cloudflare-blocked). See [panther_pilot.md](panther_pilot.md).
- MicroMap (NC-SA): removed entirely; `celldesigner.py` remains as a generic CellDesigner importer.
- `data/references/`, `Ref/`, `Presentation/`, `results/`: not tracked.

## 4. Recommended actions, in order

1. ~~Remove Yeast-maps / yeast_original~~ (done). Human-maps re-admitted under CC BY 4.0 via official Escher JSON.
2. Decide how to handle the BiGG-licensed example models: obtain written confirmation, replace them with permissively licensed equivalents, or state the non-commercial restriction in the README. (Escher community-maps layouts themselves are CC0 per upstream LICENSE.)
3. ~~Record the upstream licence of the `sbml2escher` files~~ (resolved: dependency removed); record the licence of adaptagrams/libavoid; `THIRD_PARTY_NOTICES.md` exists.
4. Add Rhea and Human-GEM citations; keep every upstream `LICENSE` next to its data.
5. Choose the software licence (not chosen yet) after step 2, since it determines what can be combined.
6. ~~Re-run the PANTHER CellDesigner pilot~~ (done via data.pantherdb.org; identity bridge + default extension).
