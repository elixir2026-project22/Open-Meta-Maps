# Third-party sources and attribution

Open Meta Maps publishes **processed structures**: the module graphs, reaction identities, geometric layouts and layout constraints that were extracted from existing expert-drawn pathway maps and pathway databases. Each source below is credited, and the licence as recorded in this project is stated. Where the recorded licence could not be confirmed, that is said explicitly. This file records provenance; it is not legal advice. See [docs/licensing_audit.md](docs/licensing_audit.md) for the evidence behind each entry.

## Expert layout sources (the Atlas)

| Source | What is used | Where it comes from | Licence as recorded |
|---|---|---|---|
| **Escher community maps** (RECON1, e_coli_core, iAF1260, iJO1366, iMM904; 22 maps) | module structures and layouts extracted from the maps; the original map files are kept as input in `data/source_inputs/raw/community/` | <https://github.com/escher/community-maps> (commit `6aa92abded6d0f9c4502ff681e14e40c5c00d61f`) | recorded as CC0-1.0; **to be confirmed** with the map authors (no licence text found in the map files) |
| **WikiPathways** | GPML pathway maps and structures extracted from them | <https://www.wikipathways.org>; terms snapshot `data/source_inputs/terms/wikipathways.html` (accessed 2026-11-05) | CC0 |
| **Reactome** (release 97) | pathway diagram and graph data, structures extracted from them | <https://reactome.org>; terms snapshot `data/source_inputs/terms/reactome.html` (accessed 2026-11-05) | CC0 for data; pathway illustrations are CC BY 4.0 and are not used |

## Reference data

| Source | What is used | Licence as recorded |
|---|---|---|
| **Rhea** (release 143) | reaction identities, directions and mappings (`data/rhea/`) | CC BY 4.0; the Rhea copyright statement is kept in `data/rhea/LICENSE.txt` |
| **BiGG Models** | example genome-scale models in `examples/` (e_coli_core, iML1515, iMM904) and `data/models/iYO844.xml` | BiGG licence; **non-commercial / educational terms apply and are to be reviewed**; the source of `iIT341.xml` is **not recorded** |
| **AraCore** (Arabidopsis core model v2.1, P. Wendering et al.) | `data/models/AraCore_v2_1.xml` | MIT (`data/models/LICENSE.AraCore`) |
| **Human-GEM** (SysBioChalmers) | `third_party/sbml2escher/models/Human-GEM.xml` | CC BY 4.0 according to the project page; please cite Human-GEM when using it |

## Software and schemas

| Component | Licence as recorded |
|---|---|
| **Escher** map schema (`third_party/sbml2escher/vendor/escher/`) | MIT |
| **sbml2escher** (model reader and currency-metabolite roles, vendored subset) | upstream licence **to be confirmed** (see `third_party/sbml2escher/NOTICE.md`) |
| **adaptagrams / libavoid** (edge routing, `data/routing/`) | LGPL; exact terms **to be confirmed**, notice to be added |

## Not included in this repository

- **MicroMap** (CC BY-NC-SA 4.0) and layouts derived from **BRENDA** are used only in a local research configuration and are not redistributed here.
- Human-maps and Yeast-maps archives were removed because of conflicting embedded licences.

## Related work

MetaCarto (<https://github.com/forxhunter/MetaCarto>) informed the design of some layout steps; its published methods were re-implemented independently, and no MetaCarto code or maps are included.
