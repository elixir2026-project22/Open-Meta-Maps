# Extra models for layout evaluation

Public stand-ins used when `elixir2026-project22/data_shared` is not reachable.
They are for generality checks of the evaluation scoring CLI only.

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

**Not fetched:** `iML1515` (excluded from this sample set); Human-GEM (too large for routine download); private `data_shared` BRENDA pathway XMLs. Point `--gem` at a local `data_shared` clone when available.
