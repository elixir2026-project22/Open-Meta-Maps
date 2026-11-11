# Data sources and licensing

The default Atlas is `data/open_core.sqlite`. Its layouts come only from:

| Source | Licence as recorded in this repository | Notes |
|---|---|---|
| Escher maps | CC0-1.0 (Escher community-maps `LICENSE`) | Atlas source label `escher`; BiGG-based map ids. See [licensing_audit.md](licensing_audit.md) |
| WikiPathways | CC0 (terms snapshot in `data/source_inputs/terms/`) | pinned GPML in `source_inputs/raw/wp/` (snapshot 20261110) |
| Human-maps | CC BY 4.0 | official `escher/` JSON only (`data/extensions/human_maps/`); SVGs with embedded BY-SA are not imported |
| BRENDA pathway maps | CC BY 4.0 | `data/extensions/brenda/`; layout possibly automatic |
| Reactome | CC0 data; some illustrations CC BY 4.0 (terms snapshot in `source_inputs/terms/`) | pinned release archive members |

Source versions, SHA-256 hashes and access dates: `data/source_inputs/source_registry.json`, `origin.json` and `data/input_manifest.json`. Every Atlas structure carries its source and licence; each generated map lists the structures it used in `licenses.json`.

Other bundled files: `data/rhea/` (Rhea reaction data), the identifier bridge `data/annotation_bridge.json` (built from the annotations of the public example models), and the example models in `examples/` the identifier tables `data/identifier_bridge.json` (BiGG metabolite to ChEBI and BiGG reaction to EC, extracted from the annotations of e_coli_core, iML1515, iMM904 and Human-GEM, with the model each entry comes from) and `data/currency_metabolites.json` (BiGG ids of currency / cofactor metabolites). The licences of the example models are to be confirmed.

## Not included

| Material | Reason |
|---|---|
| Human-maps companion SVG | SVG exports embed BY-SA; only official Escher JSON is imported |
| Yeast-maps (GPL-3.0) / original yeast map | Removed |
| MicroMap | Removed (NC-SA) |
| PANTHER Pathway SBML | Pilot acquisition blocked (Cloudflare 403); not imported — see [panther_pilot.md](panther_pilot.md) |
| `data/references/` | Local acquisition cache; not redistributed with the checkout |

MetaCarto is cited only as related work; none of its code or maps is included.

Full results and open items: [licensing_audit.md](licensing_audit.md).

Providing a source reference does not by itself grant the right to redistribute third-party data; each source's terms apply.
