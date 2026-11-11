# Extra models for layout evaluation

Public stand-ins used when `elixir2026-project22/data_shared` is not reachable.
They are for generality checks of the evaluation tool only; they are **not** added
to the Atlas, and Atlas sources are never excluded because of evaluation overlap.

Model files are **not** committed. Download on demand into the gitignored
`downloaded/` folder:

```bash
python3 examples/eval_models/fetch_eval_models.py
python3 examples/eval_models/fetch_eval_models.py --list
```

Licence for every entry was checked before listing; unclear licences are omitted
(see `skipped` in `manifest.json`).

| id | Role | Licence | Source |
|---|---|---|---|
| `e_coli_iCH360` | Compact *E. coli* GEM | MIT | [marco-corrao/iCH360](https://github.com/marco-corrao/iCH360) |
| `Mycoplasma_iPS189` | Small bacterial GEM (*M. genitalium*) | CC0-1.0 | BioModels `MODEL1507180052` |
| `R-HSA-77289` / `R-HSA-70171` | Reactome pathway SBML | CC0-1.0 (Reactome data) | Reactome Content Service |
| `F001-glycolysis` (+ `.sbgn`) | Pathway SBML + hand-drawn reference | CC-BY-4.0 | [metabolismregulation gallery](https://metabolismregulation.github.io/gallery/) |

**Not fetched:** `iML1515` (skipped by design); Human-GEM (size); private
`data_shared` BRENDA pathway XMLs (repository returned 404). Point `--gem` at a
local `data_shared` clone when available.
