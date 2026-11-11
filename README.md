# Open Meta Maps

Open Meta Maps turns a genome-scale metabolic model (GEM) into a complete, editable metabolic map. It collects public community-drawn pathway layouts in an **Atlas**, reuses matching curated structures on the input model, and completes everything else algorithmically with Bézier connections.

ELIXIR BioHackathon Europe 2026 **Project 22** (Open Pathway Maps). Layout benchmarks and the evaluation CLI are development aids, not the project goal.

## Workflow

<img src="docs/images/workflow.svg" alt="Workflow: an input model is matched against the Atlas of curated pathway structures; the Layout Composer keeps a curated skeleton, completes the rest with Bézier connections, and reports layout quality." width="700">

*Diagram source: `python3 tools/generate_workflow_svg.py` (optional `--check` writes a 700 px PNG preview).*

Reactions of the input model are matched to curated structures in the Atlas. Matched geometry is kept as drawn; every other reaction is placed algorithmically. Nothing is dropped: all reactions and metabolites of the model appear in the map.

## Key features

- **Complete GEM coverage** — every reaction and metabolite is drawn; stoichiometry and compartments are unchanged.
- **Curated pathway layouts** — structures and geometry reused from public community maps in the Atlas (see sources below).
- **Bézier layout completion** — unmatched reactions are arranged algorithmically (`layout_style=curated` by default; `--layout-style legacy` is CLI-only for orthogonal libavoid).
- **Transparent provenance** — Escher JSON and SVG, plus a per-reaction record (curated, topology-only, EC suggestion, bridge, or generated).

## Quick start

Requirements: Python 3.11+, Node.js with npm, a C++ compiler (`clang++`) and network access for the first setup. SBML via libSBML; COBRA / BiGG JSON also supported.

```sh
git clone https://github.com/yangziwei96/OpenMetaMaps.git
cd OpenMetaMaps

python3 -m pip install -r layout_atlas/requirements.txt
npm ci --prefix layout_atlas/runtime
python3 project.py setup        # one-time: build the edge routers
```

### Map Generator

```sh
python3 project.py tool
# or: ./start_tool.command
```

Open <http://127.0.0.1:8765/tool/index.html>. Upload an SBML or COBRA JSON model, or pick an example (`e_coli_core`, `iIT341`, `iMM904`, `iML1515`). The viewer supports pan, wheel/pinch zoom and a minimap. The server exits when the last browser tab closes (~8 s grace); pass `--keep-running` to disable that.

### Command line

```sh
python3 project.py compose --gem examples/e_coli_core.xml --out results/e_coli_core/ --source-paths
```

Layout evaluation (development aid): `pip install -r layout_atlas/requirements-eval.txt`, then `python3 project.py evaluate --gem MODEL.xml --map DIR/map.json`. See [docs/layout_evaluation.md](docs/layout_evaluation.md).

## How it works

- **Atlas** — reusable pathway structures from community-drawn maps (topology, reaction identity, geometry, licence).
- **Matching** — EC currentisation, Rhea metabolite scoring, name normalisation, weighted Jaccard, vicinity propagation; EC similarity extension optional and **off by default** (`--ec`). See [docs/reaction_map_matching.md](docs/reaction_map_matching.md).
- **Layout Composer** — selects non-overlapping curated structures as a skeleton, places the rest, draws Bézier connections, reports layout quality. Curated reuse, topology-only references and EC suggestions are never merged into one reuse rate. See [docs/methodology.md](docs/methodology.md).

## Layout after curated reuse

Once matching has chosen Atlas structures, `tool_pipeline.draw_map` / `draw_global_map` run the stages below. Small models (≤150 reactions) take the full path; mid-size models (≤3 000) use a global plan (`global_layout.plan_global`) that starts with straight edges and then applies the same curated Bézier pass; larger models keep a skeleton and draw ~50-reaction views on demand.

| Stage | What it does | Main code |
|---|---|---|
| **Rigid curated skeleton** | Scale each matched structure and pack blocks as rigid groups with shared metabolites as junctions (`drawing.pack_modules`, `geometry.module_scale`). Network bridges get anchors between blocks. | `drawing.py`, `geometry.py` |
| **Complete uncovered reactions** | Place free nodes around fixed anchors (`drawing.complete`, `completion.complete_spaced`). Optional context sharing for currency-only reactions (`context.find` / `place`). Anchor uncovered reactions beside curated glyphs when possible (`completion_context.place`). | `completion.py`, `context.py`, `completion_context.py` |
| **Metabolite instances** | Decide how many visual copies a non-currency metabolite gets to shorten long edges (`instances.optimize`). Currency/cofactors stay private `side:` glyphs per reaction (`drawing.display_scene`, `data/currency_metabolites.json`). | `instances.py`, `drawing.py` |
| **Orientation** | Put free reactions between their core substrates and products; place currency copies on the substrate/product side (`orientation.apply`, default `full`). Curated reactions and `copy:` glyphs are not moved. | `orientation.py` |
| **Composition & peripherals** | Pull related stray components toward the main network and shelf independents (`composition.compose`, `OPM_COMPOSE=on`). Attach exchange/simple transport modules to intracellular anchors (`peripheral.attach`). Optional relocation of exchange/transport *inside* large curated blocks (`block_relocation.relocate_sub`, `OPM_EXPERT_SUBMODULE=branch`). Whole-block relocation is off unless `OPM_BLOCK_RELOC=on`. | `composition.py`, `peripheral.py`, `block_relocation.py` |
| **Spacing / relaxation** | Soft-preserve curated shape, shorten edges, and project minimum node spacing (`geometry.complete_scene` / `relax`; sparse path for large maps). | `geometry.py`, `global_layout.py` |
| **Bézier connections** | Default `layout_style=curated`: fit algorithm-placed edges to cubic Béziers from reaction templates, scoring obstacles, crossings, curvature and detour (`edge_style.apply`). Optional `--source-paths` reuses original Escher curves on matched edges (`source_paths`). `--layout-style legacy` keeps orthogonal libavoid. | `edge_style.py`, `source_paths.py`, `routing.py` |
| **Labels** | After curves exist, move labels among candidates to cut collisions (`label_placement.optimize`). | `label_placement.py` |
| **Integrity** | Model hash, every reaction drawn, visual stoichiometry (`tool_pipeline.check_integrity`); protected shared glyphs must not split. | `tool_pipeline.py`, `completion_context.py` |

**Quality evaluation** (development aid): compose writes `layout_quality.evaluate` into `layout_issues.json` / `report.json` (curve crossings by class, collisions, lengths, orientation share, structure continuity). `python3 project.py evaluate` adds Project 22 layout-benchmark metrics via `evaluation/score.py`. See [docs/layout_evaluation.md](docs/layout_evaluation.md).

**Determinism:** fixed layout seed `20261105` (`core.SEEDS`), ELK/fCoSE seeds, and `PYTHONHASHSEED=0` on pipeline re-exec. There is no compose `--seed` flag yet; `--compare-map` on evaluate measures displacement only.

### Layout rules (status)

| Rule | Status | Notes |
|---|---|---|
| Avoid crossings | Partial | Minimised in `edge_style` / quality gate; GEMs still show many algo crossings |
| Avoid overlaps | Partial | Spacing + label pass + gate; overlaps grow with map size |
| Avoid long lines | Partial | Instances, composition, peripherals help; currency/biomass spokes remain |
| Place participants consistently | Partial | `orientation` on by default; measured share still well below curated reference maps |
| Respect compartments | Not yet | Compartment used in packing/shelves; cross-compartment species metric still high |
| Aspect / compactness | Partial | Shelf aspect target; packing of independent components is still linear shelves |
| Reproducibility | Partial | Fixed seeds; no compose seed CLI; evaluate `--compare-map` only |
| Build on previous layout | Not yet | Displacement metric exists; not a compose objective |
| Scale to GEMs | Partial | Global ≤3 000 + tiled views; large maps start routing-degraded |

**Known gaps** (also [docs/limitations.md](docs/limitations.md), [docs/layout_evaluation.md](docs/layout_evaluation.md)): 2D component packing (in progress separately), harder orientation (incl. transporters), compartment-aware packing, more aggressive currency/biomass handling, crossing-aware composition acceptance, compose seed flag. Experimental off by default: `OPM_REGIONS`, `OPM_LOWCONF`, whole-block `OPM_BLOCK_RELOC`.

**Detailed-metrics backlog (not yet optimised / 待优化):** ring/cycle placement for algorithm-placed cycles (circularity / direction / order); species in another compartment (compartment packing); reproducibility / build-on-previous-layout (detail row when `evaluate --compare-map` is used; not a main card). Global flow-direction share is a **debug-only** diagnostic (hidden by default) — not an optimisation target, because forcing one reading direction would fight curated orientations. Not added from layout-transfer-check: raw density, straight-line crossing counts, catalyst-beside, arrow-weight.

## Outputs

| Output | Description |
|---|---|
| `map.json` | Editable Escher map |
| `map.svg` | Vector visualisation |
| `report.json` | Coverage and layout-quality measures |
| `reaction_provenance.json` | Per-reaction provenance and matching evidence |
| `layout_issues.json` | Crossings, collisions and other layout problems |
| `licenses.json` | Licensing metadata of Atlas structures used |

## Data sources and licensing

Default Atlas (`data/open_core.sqlite` plus default extensions). Counts from [`data/atlas_registry.json`](data/atlas_registry.json) (generated 2026-11-11; About page uses the same file):

| Source | Modules | Licence |
|---|---:|---|
| Escher (22 community maps) | 13,047 | CC0-1.0 |
| Reactome (release 97, 86 diagrams) | 3,405 | CC0 (data) |
| WikiPathways (20261110; 1,245 semantic + 29,673 geometry-only) | 30,918 | CC0-1.0 |
| Human-maps (130 official Escher JSON) | 8,540 | CC BY 4.0 |
| BRENDA (194 pathway maps; layout possibly automatic) | 11,703 | CC BY 4.0 |
| Metabolism Regulation (22 SBGN-ML maps) | 423 | CC BY 4.0 |
| PANTHER (3.4/CD4.1, 8 maps) | 238 | CC BY 4.0 |
| iCH360 (10 custom Escher maps; MIT) | 858 | MIT |
| **Total (default load)** | **69,132** | |

Not included: MicroMap, Yeast-maps, KEGG; iCH360 `compressed_map.json` (CP pseudo-reactions). Attribution: [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). Audit notes: [docs/licensing_audit.md](docs/licensing_audit.md).

**Code** is [MIT](LICENSE). **Atlas data and generated maps** mix CC0, CC BY 4.0 and MIT sources; keep the attributions in THIRD_PARTY_NOTICES.md.

## Documentation

- [Installation](docs/installation.md)
- [Data sources](docs/data_sources.md) · [iCH360 import](docs/ich360_import.md) · [Optional packages](docs/data_packages.md) · [Atlas layout](data/README.md)
- [Reaction–map matching](docs/reaction_map_matching.md) · [EC extension](docs/ec_extension_recommendation.md)
- [Methodology](docs/methodology.md) · [Limitations](docs/limitations.md)
- [Layout evaluation](docs/layout_evaluation.md)
- [Code overview](layout_atlas/README.md) · [About / Provenance](layout_atlas/results/about/index.html) (regenerate with `python3 layout_atlas/navigation.py`)

## Project status

Under active development. Focus: global pathway organisation, less visual overlap, consistent layouts across GEMs.

## Citation and acknowledgements

A citation will be added with the first release. Until then, cite this repository and the sources in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

Developed during ELIXIR BioHackathon Europe 2026, building on the Escher, WikiPathways, Reactome and wider metabolic-modelling communities.
