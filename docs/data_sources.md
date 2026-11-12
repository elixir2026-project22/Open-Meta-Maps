# Data sources and licensing

The default Atlas is `data/open_core.sqlite` plus the in-tree extensions under `data/extensions/`. Module counts: [`data/atlas_registry.json`](../data/atlas_registry.json). How maps are cut into modules (Stage B `diagram_modules` recut), deduplicated, given canonical reaction ids, stamped with extraction provenance, and ranked for compactness: [atlas_data_management.md](atlas_data_management.md). Layouts come from:

| Source | Licence as recorded in this repository | Notes |
|---|---|---|
| Escher maps | CC0-1.0 (Escher community-maps `LICENSE`) | Atlas source label `escher`; in-repo `escher_community_windows` cutter; BiGG-based map ids. See [licensing_audit.md](licensing_audit.md) |
| WikiPathways | CC0 (terms snapshot in `data/source_inputs/terms/`) | pinned GPML in `source_inputs/raw/wp/` (snapshot 20261110); semantic + geometry-only windows (Stage B recut) |
| Human-maps | CC BY 4.0 (JSON) | official `escher/` JSON only (`data/extensions/human_maps/`); MAR ids remapped via Human-GEM crosswalk. Companion SVG exports are CC BY-SA and are outside the import policy (JSON only) |
| BRENDA pathway maps | CC BY 4.0 | `data/extensions/brenda/`; layout possibly automatic |
| Metabolism Regulation gallery | CC BY 4.0 | hand-drawn SBGN-ML via `layout_atlas/sbgn.py` (`data/extensions/metabolismregulation/`); Rhea when ChEBI uniquely matches; see [metabolismregulation_import.md](metabolismregulation_import.md) |
| PANTHER Pathway | CC BY 4.0 | CellDesigner CD4.1 via `celldesigner.py` + `identity_labels.py` (`data/extensions/panther/`); native ids kept; see [panther_pilot.md](panther_pilot.md) |
| iCH360 Visualisation | MIT | custom Escher JSON via `database.escher_diagram` + `diagram_modules` (`data/extensions/ich360/`); import uses `full_map` and pathway maps; see [ich360_import.md](ich360_import.md) |
| Reactome | CC0 data; some illustrations CC BY 4.0 (terms snapshot in `source_inputs/terms/`) | pinned release archive members; Stage B `diagram_modules` recut |

Source versions, SHA-256 hashes, retrieval dates and extraction methods: `data/source_inputs/source_registry.json` and `data/input_manifest.json`. Every Atlas structure carries its source, licence and an `extraction` block; each generated map lists the structures it used in `licenses.json`.

Other bundled files: `data/rhea/` (Rhea reaction data); `data/identifier_bridge.json` (BiGG↔ChEBI and BiGG↔EC from e_coli_core, iML1515, iMM904 and Human-GEM — the single identifier store; Mapper crosswalk derived at load time, optional local `annotation_bridge.json.gz` cache); and `data/currency_metabolites.json`.

Licence exclusions and redistribution decisions: [licensing_audit.md](licensing_audit.md) and [THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md).

Providing a source reference does not by itself grant the right to redistribute third-party data; each source's terms apply.
