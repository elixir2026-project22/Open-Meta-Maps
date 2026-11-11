# Extra models for layout evaluation

Public stand-ins used when `elixir2026-project22/data_shared` is not reachable from this environment.
They are for generality checks of the evaluation tool only; they are **not** added to the Atlas, and
Atlas sources are never excluded because of evaluation overlap.

| File | Role | Source |
|---|---|---|
| `e_coli_iCH360.xml` | Compact *E. coli* GEM | [marco-corrao/iCH360](https://github.com/marco-corrao/iCH360) |
| `Mycoplasma_iPS189.xml` | Small bacterial GEM (*M. genitalium*) | BioModels `MODEL1507180052` (iPS189) |
| `R-HSA-77289.sbml` | Reactome pathway SBML (benchmark manifest sample) | Reactome Content Service |
| `R-HSA-70171.sbml` | Human glycolysis pathway (Reactome) | Reactome Content Service |
| `F001-glycolysis.xml` + `.sbgn` | Pathway SBML with hand-drawn reference | [metabolismregulation gallery](https://metabolismregulation.github.io/gallery/) (CC BY 4.0) |

**Not run here:** `iML1515` (skipped by design); Human2 / full Human-GEM (size not feasible in this pass);
private `data_shared` BRENDA pathway XMLs (repository returned 404). Re-run with a local
`data_shared` clone when available: point `--gem` at those files; do not drop Atlas sources.
