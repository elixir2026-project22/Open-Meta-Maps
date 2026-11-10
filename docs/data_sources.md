# Data sources and licensing

The default Atlas is `data/open_core.sqlite`. Its layouts come only from:

| Source | Licence as recorded in this repository | Notes |
|---|---|---|
| Escher community maps | recorded as CC0-1.0 in the Atlas metadata, **not verified** | no licence text found for the maps; BiGG-based, and the BiGG licence restricts commercial use. See [licensing_audit.md](licensing_audit.md) |
| WikiPathways | CC0 (terms snapshot in `data/source_inputs/terms/`) | pinned GPML in `source_inputs/raw/wp/` |
| Reactome | CC0 data; some illustrations CC BY 4.0 (terms snapshot in `source_inputs/terms/`) | pinned release archive members |

Source versions, SHA-256 hashes and access dates: `data/source_inputs/source_registry.json`, `origin.json` and `data/input_manifest.json`. Every Atlas structure carries its source and licence; each generated map lists the structures it used in `licenses.json`.

Other bundled files: `data/rhea/` (Rhea reaction data), the identifier bridge `layout_atlas/results/micromap_extension/annotation_bridge.json` (built from the annotations of the public example models), and the example models in `examples/` the identifier tables `data/identifier_bridge.json` (BiGG metabolite to ChEBI and BiGG reaction to EC, extracted from the annotations of e_coli_core, iML1515, iMM904 and Human-GEM, with the model each entry comes from) and `data/currency_metabolites.json` (BiGG ids of currency / cofactor metabolites). The licences of the example models are to be confirmed.

## Not included

| Material | Reason |
|---|---|
| MicroMap extension | CC BY-NC-SA 4.0 (non-commercial, share-alike). It could be released as a separate data package under those terms, but it must not be merged into the core Atlas. |
| BRENDA-derived layouts | Layout rights not yet verified; only the rebuild scripts are included. |
| `data/references/` | Mixed licences; to be audited per file. |

MetaCarto is cited only as related work; none of its code or maps is included.

Full results and open items: [licensing_audit.md](licensing_audit.md).

Providing a source reference does not by itself grant the right to redistribute third-party data; each source's terms apply.
