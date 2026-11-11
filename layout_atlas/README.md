# layout_atlas: the Composer and the Atlas

This directory is the Python package behind Open Meta Maps. It builds a metabolic map from a genome-scale model by matching its reactions to expert-drawn structures (the Atlas), assembling those structures into a skeleton, laying out the remaining reactions algorithmically and drawing the result with Bézier connections. Installation and the quick start are in the [root README](../README.md); method notes are in [`docs/`](../docs/).

## Entry points

| Entry | What it does |
|---|---|
| `python3 project.py tool` | local web tool (`tool_server.py`, page in `results/tool/index.html`): upload a model, run the pipeline, preview, download |
| `python3 project.py compose --gem MODEL.xml --out DIR --source-paths --layout-style expert` | the same pipeline from the command line (`compose_cli.py` → `tool_pipeline.py`) |
| `python3 project.py data list` / `fetch` | optional data packages, see [docs/data_packages.md](../docs/data_packages.md) |
| `python3 project.py setup` | one-time build of the edge routers (needs `g++` or `clang++`, see below) |
| `python3 project.py catalog` | rebuild the module catalog pages from the Atlas |

A run writes `map.json` (Escher), `map.svg`, `report.json`, `reaction_provenance.json`, `layout_issues.json`, `licenses.json`, `skeleton.json` and `modules.json` into the output folder. `reaction_provenance.json` labels every reaction as expert-derived, topology-only suggestion, functional (EC) suggestion, bridge or generated; these are reported separately and never merged into one reuse rate.

## Pipeline and main modules

| Stage | Modules |
|---|---|
| Model input and identifiers | `model_io.py` (SBML via libSBML, COBRA / BiGG JSON), `identifier_bridge.py` (BiGG → ChEBI / EC tables from `data/identifier_bridge.json`), `../tools/build_identifier_bridge.py` (rebuilds them from the public models) |
| Matching and selection | `retrieval.py`, `semantic_transfer.py`, `ec_channel.py`, `ec_provenance.py`, `database.py` (reads `data/atlas.sqlite`) |
| Skeleton and global layout | `drawing.py`, `global_layout.py`, `pathway_backbone.py`, `completion.py` |
| Context and anchors for unmatched reactions | `context.py`, `completion_context.py`, `instances.py` (metabolite instances and protected shared connections), `orientation.py` |
| Whole-map composition | `composition.py` (recovers stray components, arranges independent modules in rows), `peripheral.py` (exchange and transport modules), `block_relocation.py` (moves expert blocks that sit far from the reactions they connect to) |
| Connections and labels | `edge_style.py` (Bézier connections from expert templates, joint routing per reaction), `reaction_templates.py`, `label_placement.py`, `line_render.py`, `routing.py` (libavoid) |
| Quality | `layout_quality.py` (crossings, collisions, long connections, composition metrics, `layout_issues.json`), `quality_gate.py` (accepts a change only if the final geometry is not worse) |
| Pages | `ui.py`, `navigation.py`, `toolpage.py`, `gallery*.py` |

Files whose hashes are recorded by saved results (`routing.py`, `drawing.py`, `geometry.py`, `retrieval.py`, `semantic_transfer.py`, `celldesigner.py`) should not be edited casually.

## Switches

The tool and `compose` expose the layout style (`--layout-style legacy|expert`, expert is what the web tool uses). The environment variables below switch individual steps; defaults are what a normal run uses.

| Variable | Default | Meaning |
|---|---|---|
| `OPM_COMPOSE`, `OPM_COMPOSE_LONG` | on | whole-map composition; relocation threshold normalised by reaction spacing |
| `OPM_PERIPHERAL` | on | attach exchange + transport modules to the pathway they belong to |
| `OPM_BLOCK_RELOC`, `OPM_RELOC_ACCEPT` | on, `balanced` | relocate expert blocks with a long connection; `strict` uses the strict quality gate |
| `OPM_EXPERT_SUBMODULE` | `branch` | move exchange + transport modules (including coupled transporters) that sit far from the rest of a large expert block next to it; the connections that cross the module boundary are re-drawn and recorded as `regenerated_expert_edges`. `on` = simple transports only, `off` = disabled |
| `OPM_RELOC_CLEAR` | off | experimental local space-making for relocations |
| `OPM_LOWCONF` | off | experimental: split topology-only structures whose reactions have unrelated neighbourhoods |
| `OPM_REGIONS`, `OPM_REGION_PLACE` | off | experimental region-level arrangement (did not improve the example maps) |
| `OPM_DROP_TOPO_BIG`, `OPM_SUMMARY_LAYOUT` | off, `ring` | experimental handling of very large reactions (biomass) |

## System Requirements

### C++ Compiler for libavoid Router

The edge router (`routing.py`) compiles a native libavoid binary. A C++11-capable compiler is required:

- **Linux**: Install `g++` (usually via `build-essential` on Debian/Ubuntu: `sudo apt-get install build-essential`)
- **macOS**: Install Xcode command-line tools: `xcode-select --install` (provides `clang++`)
- **Windows**: Not tested; consider WSL with Linux tools

The build prefers `g++` when available, falling back to `clang++`.

## Tests

```sh
cd layout_atlas
python3 test_orientation.py
python3 test_ec_provenance.py
python3 test_reaction_map_matching.py
python3 map_matching_eval.py          # BRENDA / Atlas match rates for e_coli_core, iIT341
```


## Data

See [`data/README.md`](../data/README.md). Sources, licences and attribution: [`THIRD_PARTY_NOTICES.md`](../THIRD_PARTY_NOTICES.md) and [`docs/licensing_audit.md`](../docs/licensing_audit.md).
