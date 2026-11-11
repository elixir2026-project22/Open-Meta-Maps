# Data sources and licensing

The default Atlas is `data/open_core.sqlite` plus the in-tree extensions under `data/extensions/`. Module counts: [`data/atlas_registry.json`](../data/atlas_registry.json). Layouts come only from:

| Source | Licence as recorded in this repository | Notes |
|---|---|---|
| Escher maps | CC0-1.0 (Escher community-maps `LICENSE`) | Atlas source label `escher`; BiGG-based map ids. See [licensing_audit.md](licensing_audit.md) |
| WikiPathways | CC0 (terms snapshot in `data/source_inputs/terms/`) | pinned GPML in `source_inputs/raw/wp/` (snapshot 20261110) |
| Human-maps | CC BY 4.0 (JSON) | official `escher/` JSON only (`data/extensions/human_maps/`); companion SVGs are CC BY-SA and are not imported |
| BRENDA pathway maps | CC BY 4.0 | `data/extensions/brenda/`; layout possibly automatic |
| Metabolism Regulation gallery | CC BY 4.0 | hand-drawn SBGN-ML via `layout_atlas/sbgn.py` (`data/extensions/metabolismregulation/`); see [metabolismregulation_import.md](metabolismregulation_import.md) |
| PANTHER Pathway | CC BY 4.0 | CellDesigner CD4.1 via `celldesigner.py` + `identity_labels.py` (`data/extensions/panther/`); see [panther_pilot.md](panther_pilot.md) |
| Reactome | CC0 data; some illustrations CC BY 4.0 (terms snapshot in `source_inputs/terms/`) | pinned release archive members |

Source versions, SHA-256 hashes, retrieval dates and extraction methods: `data/source_inputs/source_registry.json` and `data/input_manifest.json`. Every Atlas structure carries its source and licence; each generated map lists the structures it used in `licenses.json`.

Other bundled files: `data/rhea/` (Rhea reaction data); `data/identifier_bridge.json` (BiGG↔ChEBI and BiGG↔EC from e_coli_core, iML1515, iMM904 and Human-GEM — the single identifier store; Mapper crosswalk derived at load time, optional local `annotation_bridge.json.gz` cache); and `data/currency_metabolites.json`. The licences of the example models are to be confirmed.

## Not included

| Material | Reason |
|---|---|
| Human-maps companion SVG | SVG exports are CC BY-SA; only official Escher JSON (CC BY 4.0) is imported |
| Yeast-maps (GPL-3.0) / original yeast map | Excluded |
| MicroMap (NC-SA) | Excluded |
| `data/references/` | Local acquisition cache; not redistributed with the checkout |

MetaCarto is cited only as related work; none of its code or maps is included.

Full results and open items: [licensing_audit.md](licensing_audit.md).

Providing a source reference does not by itself grant the right to redistribute third-party data; each source's terms apply.
