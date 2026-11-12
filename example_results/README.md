# Example results

Maps composed from the bundled example models with the default settings of `python3 project.py compose ... --source-paths`. Every reaction and metabolite of the model is drawn.

| Model | Folder | Reactions | Metabolites |
|---|---|---:|---:|
| *Escherichia coli* core | [`e_coli_core/`](e_coli_core/) | 95 | 72 |
| *Helicobacter pylori* iIT341 | [`iIT341/`](iIT341/) | 554 | 485 |

## Files

Each model folder holds the final map and its reports.

| File | Content |
|---|---|
| `map.json` | Map in Escher JSON format; opens in [Escher](https://escher.github.io) and in the Map Generator |
| `map.svg` | Vector drawing of the map |
| `map.png` | Preview cropped to content, 2800 px wide; labels kept only when readable at that width |
| `licenses.json` | Licence and source of every reused curated structure |
| `modules.json` | Reused curated structures with their source map and reactions |
| `reaction_provenance.json` | Layout provenance of every reaction (curated structure or generated) |
| `layout_issues.json` | Crossings, overlaps and other layout issues with coordinates |
| `report.json` | Run settings, Atlas coverage, integrity checks and layout quality |

## Regenerate

```sh
python3 project.py compose --gem examples/e_coli_core.xml --out results/e_coli_core/ --source-paths
python3 project.py compose --gem examples/iIT341.xml --out results/iIT341/ --source-paths
PYTHONPATH=layout_atlas python3 tools/render_map_preview.py results/e_coli_core/map.json results/e_coli_core/map.png
```

## Licences

Reaction and metabolite content comes from the input models in [`examples/`](../examples/) and follows the [BiGG licence](http://bigg.ucsd.edu/license). Reused curated structures keep their source licences (CC0, CC BY 4.0, MIT), listed per structure in `licenses.json`; CC BY 4.0 sources require the attribution given there. See [THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md).
