# Licensing audit

Status: audit of what the repository currently redistributes, dated 2026-11-10. Evidence is either local (term snapshots and provenance files in `data/source_inputs/`) or a public web page checked on that date. **This is an engineering review, not legal advice.** Items marked *open* need a decision or confirmation before a public release.

## 1. Atlas content (`open_core.sqlite`)

All 17,221 modules in `open_core.sqlite` come from three sources (checked by reading the module payloads):

| Source | Modules | Licence | Evidence | Result |
|---|---|---|---|---|
| WikiPathways | 4,718 | CC0 | `source_inputs/terms/wikipathways.html` (CC0 waiver for content) | OK |
| Reactome | 3,405 | CC0 for data and derived files; CC BY 4.0 for pathway illustrations and artwork | `source_inputs/terms/reactome.html` | OK for layout data; if Reactome illustrations are ever added, CC BY 4.0 attribution is required |
| Escher community maps | 9,098 | recorded in the Atlas as CC0-1.0 | none found | **Open**: see below |

### Escher community maps (open)

The 22 maps come from BiGG-based Escher maps (e_coli_core, iAF1260, iJO1366, iMM904, RECON1). The `CC0-1.0` label in the Atlas metadata is not backed by any licence text in the raw files or in the provenance records. Checked externally: the Escher software is MIT-licensed, but the Escher project page says nothing about the licence of published maps. The BiGG site states its data is available under the "BiGG License" (<https://bigg.ucsd.edu/license>), which allows copying and redistribution for educational, research and non-profit purposes and says commercial use is not covered. Whether the Escher maps fall under it is not stated.

Action: confirm the licence of each of the 22 maps with their authors or the Escher/BiGG maintainers; until then treat these maps as "redistribution terms unconfirmed" and do not describe the Atlas as fully CC0.

## 2. Third-party files tracked in the repository

| Path | Licence | Result |
|---|---|---|
| `data/rhea/` | CC BY 4.0 (`rhea/LICENSE.txt`: copyright statement must be reproduced with each copy) | OK, keep `LICENSE.txt`; credit Rhea in the README |
| `data/source_inputs/archives/wp-*.zip` and `raw/wp/` | CC0 (WikiPathways) | OK |
| `data/reactome/` | CC0 data (Reactome v97) | OK |
| `data/models/AraCore_v2_1.xml` | MIT (`LICENSE.AraCore`, Philipp Wendering) | OK |
| `data/models/iYO844.xml` | from BiGG (BiGG License) | **Open**: non-commercial/educational terms |
| `examples/` (e_coli_core, iML1515, iMM904, iIT341) | BiGG models (BiGG License) | **Open**: same as above; iIT341 source not recorded |
| `data/identifier_bridge.json` | annotations extracted from BiGG models and Human-GEM (CC-BY-4.0 according to its GitHub page); the Human-GEM file itself is no longer in the repository | cite BiGG and Human-GEM; terms of the source models apply |
| `data/escher_schema/map-schema-1-0-0.json` | Escher, MIT (`LICENSE` copied) | OK |
| ~~`model_input.py`, `central_reference.json` (sbml2escher)~~ | removed 2026-11-10: the SBML reader is now libSBML-based (`layout_atlas/model_io.py`), the currency list is `data/currency_metabolites.json` | **Resolved** |
| `data/routing/source.zip` (adaptagrams, libavoid) | adaptagrams is LGPL; exact terms not re-checked | **Open**: confirm and add the LGPL notice |
| `human_maps.zip`, `yeast_maps.zip` (Human-maps, Yeast-maps) | recorded as "repository CC-BY-4.0 / GPL-3.0; embedded SVG BY-SA-4.0 conflict" | **Removed** (2026-11-10), together with their 251 map observations; they contributed 0 modules to the open core. `input_manifest.json`, `source_registry.json` and `extensions/source_paths/SOURCE.json` were updated to match. `acquire.py` still contains the download code for them (not run by default); they are not redistributed |
| `data/extensions/{units,source_paths,reaction_geometry_templates}` | derived from the open-core maps | inherits the Escher-map status above |
| `layout_atlas/results/micromap_extension/annotation_bridge.json` | derived from model annotations (identifier links) | inherits the model licences above; name is historical, contains no MicroMap layout |

## 3. Not in the repository (verified by `.gitignore` and tracked-file listing)

- MicroMap: CC BY-NC-SA 4.0 (non-commercial, share-alike). Not tracked. May be released as a separate data package under those terms only.
- BRENDA-derived layouts: the source note says the BRENDA online version is CC BY 4.0 but the pathway-map-specific terms are not stated; not tracked.
- `data/references/`, `Ref/`, `Presentation/`, `results/`: not tracked.
- Code files `micromap_*.py`, `celldesigner.py` are tracked; they contain no MicroMap data.

## 4. Recommended actions, in order

1. ~~Remove `human_maps.zip`, `yeast_maps.zip` and their derived observations~~ (done; repository history rewritten).
2. Decide how to handle the BiGG-licensed example models and the Escher community maps: obtain written confirmation, replace them with permissively licensed equivalents, or state the non-commercial restriction in the README.
3. ~~Record the upstream licence of the `sbml2escher` files~~ (resolved: dependency removed); record the licence of adaptagrams/libavoid; `THIRD_PARTY_NOTICES.md` exists.
4. Add Rhea and Human-GEM citations; keep every upstream `LICENSE` next to its data.
5. Choose the software licence (not chosen yet) after step 2, since it determines what can be combined.
