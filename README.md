# Open Meta Maps

Open Meta Maps turns a genome-scale metabolic model (GEM) into a complete, editable metabolic map. It has two parts: a curated-structure **Atlas** built from public community maps, and a **layout optimisation** pipeline that reuses matching Atlas geometry and completes the rest with Bézier connections.

ELIXIR BioHackathon Europe 2026 **[Project 22](https://github.com/elixir2026-project22)** (Open Pathway Maps).

## Workflow

<img src="docs/images/workflow.svg" alt="Workflow: an input model is matched against the Atlas of curated pathway structures; the Layout Composer keeps a curated skeleton, completes the rest with Bézier connections, and reports layout quality." width="700">

*Diagram source: `python3 tools/generate_workflow_svg.py` (optional `--check` writes a 700 px PNG preview).*

Every reaction and metabolite of the model appears in the map. Stoichiometry and compartments are unchanged.

## Quick start

Requirements: Python 3.11+, Node.js with npm, a C++ compiler (`clang++`) and network access for the first setup. SBML via libSBML; COBRA / BiGG JSON also supported.

```sh
git clone https://github.com/elixir2026-project22/Open-Meta-Maps.git
cd Open-Meta-Maps

python3 -m pip install -r layout_atlas/requirements.txt
npm ci --prefix layout_atlas/runtime
python3 project.py setup        # one-time: build the edge routers
```

### Map Generator

```sh
python3 project.py tool
# or: ./start_tool.command
```

Open <http://127.0.0.1:8765/tool/index.html>. Upload an SBML or COBRA JSON model, or pick an example (`e_coli_core`, `iIT341`, `iMM904`, `iML1515`). The viewer supports pan, wheel/pinch zoom and a minimap. The server exits when the last browser tab closes (~8 s grace); pass `--keep-running` to disable that. The same server also serves **Atlas Browser** (`/atlas_browser/index.html`) to search curated modules and download Escher JSON with a provenance/licence manifest (search index built on first use under `data/cache/`).

### Command line

```sh
python3 project.py compose --gem examples/e_coli_core.xml --out results/e_coli_core/ --source-paths
```

Rebuild after a small model edit while keeping the previous layout: pass `--previous` with the prior `map.json`. Matched reaction midmarkers and primary metabolites are pinned by stable id; when at least 80% of reactions are already placed, module packing is skipped and the previous Escher geometry is reused. New reactions are grafted on the chord between drawn substrate/product participants (opposite sides); leftovers without anchors are still parked so coverage stays complete. A local nudge of grafted midmarkers is kept when crossings do not worsen.

```sh
python3 project.py compose --gem examples/e_coli_core_edited.xml --out results/e_coli_core_edited/ --source-paths --previous results/e_coli_core/map.json
```

Layout evaluation: `pip install -r layout_atlas/requirements-eval.txt`, then `python3 project.py evaluate --gem MODEL.xml --map DIR/map.json`. See [docs/layout_evaluation.md](docs/layout_evaluation.md).

## 1. Data construction — the curated-structure Atlas

The Atlas stores reusable pathway modules (topology, reaction identity, geometry, licence) cut from public community maps. Default load: `data/open_core.sqlite` plus in-tree extensions. Counts from [`data/atlas_registry.json`](data/atlas_registry.json) (generated 2026-11-12; About uses the same file).

<img src="docs/images/atlas/atlas_overview.png" alt="Overview of the Atlas: one panel per source, sized by its unique reaction sets and filled with every unique module, with the largest distinct shapes drawn as framed tiles." width="700">

*One panel per source; panel area scales with its unique reaction sets, and every unique set is one tile. Framed tiles show the largest distinct shapes. Titles give the cut and unique-after-dedup counts. Regenerate: `PYTHONPATH=layout_atlas python3 tools/atlas_module_figures.py overview`.*

### Sources and licences

| Source | Modules | Unique after dedup | Licence |
|---|---:|---:|---|
| Escher (22 community maps) | 13,047 | 10,626 | CC0-1.0 |
| Reactome (release 97, 115 diagrams) | 1,811 | 1,800 | CC0 (data) |
| WikiPathways semantic (20261110) | 669 | 565 | CC0-1.0 |
| WikiPathways geometry-only | 15,611 | 8,489 | CC0-1.0 |
| Human-maps (130 official Escher JSON) | 4,690 | 4,551 | CC BY 4.0 |
| BRENDA (194 pathway maps; layout possibly automatic) | 11,703 | 9,913 | CC BY 4.0 |
| Metabolism Regulation (22 SBGN-ML maps) | 229 | 172 | CC BY 4.0 |
| PANTHER (3.4/CD4.1, 8 maps) | 138 | 116 | CC BY 4.0 |
| iCH360 (10 custom Escher maps) | 436 | 426 | MIT |
| **Total (default load)** | **48,334** | **36,103** | |

Totals count **cut modules** (reaction windows), not independent drawings: several modules can come from the same source map, and many windows nest inside larger ones from the same map. Dedup uses canonical reaction-id sets (geometry fingerprint for WikiPathways geometry-only windows). Stats: [`docs/images/atlas/atlas_overlap.json`](docs/images/atlas/atlas_overlap.json) via `tools/atlas_module_figures.py`.

Attribution: [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). Source notes: [docs/data_sources.md](docs/data_sources.md).

### Example modules

<img src="docs/images/atlas/atlas_examples.svg" alt="Grid of example Atlas modules with source, map name and reaction count captions." width="700">

*Two–three representative modules per source. Regenerate: `PYTHONPATH=layout_atlas python3 tools/atlas_module_figures.py examples`.*

### Extraction, identity and selection

- **Cut** — documented extractors: in-repo Escher community cutter (`escher_community_windows`, parity-tested); Stage B `diagram_modules` / `geometry_modules` (size-scaled roots, near-nested dedupe) for Reactome, WikiPathways, Human-maps, Metabolism Regulation, PANTHER and iCH360; BRENDA’s own extractor.
- **Deduplicate** — near-duplicate windows grouped into core families and geometric variants; Stage B also drops near-nested modules (`nest_factor=2`).
- **Canonical ids** — unique crosswalks map reaction references to BiGG, MetaNetX and/or Rhea where available (Human-maps via Human-GEM; Metabolism Regulation via participant ChEBI; Escher / iCH360 already BiGG).
- **Provenance** — each module records source, map, licence and an `extraction` block (method, version, parameters, counts).
- **Selection** — coverage ranks first; among same-source biochemical matches with the same reaction set, the denser drawing is kept.

Worked examples (cut boundaries, ID remap, cross-source dedup, provenance, compact selection): [docs/atlas_data_management.md](docs/atlas_data_management.md#worked-examples). On-disk layout: [data/README.md](data/README.md).

### Download the Atlas

```sh
python3 project.py data export                                       # whole Atlas, unique after dedup -> dist/openmetamaps-atlas-<date>.tar.gz
python3 project.py data export --source panther --licence CC-BY-4.0  # one source, one licence (both repeatable)
python3 project.py data export --all-windows --out atlas_all.tar.gz  # every cut window instead of unique reaction sets
```

The archive holds one Escher JSON per module (`modules/<source>/<module_id>.json`), `manifest.json` (source map, licence, attribution text, reaction count, canonical reaction ids, extraction provenance and `also_in` for identical reaction sets in other sources), a short README and `THIRD_PARTY_NOTICES.md`. The default export has 36,103 modules; re-running it gives a byte-identical archive. To pick individual modules, search the Atlas Browser and download the selection as the same Escher JSON with the same manifest fields. Details: [docs/data_packages.md](docs/data_packages.md).

---

## 2. Layout optimisation algorithm

### Matching and curated reuse

Reactions of the input model are matched to Atlas modules (EC currentisation, Rhea metabolite scoring, name normalisation, weighted Jaccard, vicinity propagation). Report curated, topology-only and EC suggestions as separate rates in `reaction_provenance.json`. EC similarity extension is **off by default** (`--ec`). See [docs/reaction_map_matching.md](docs/reaction_map_matching.md) · [docs/methodology.md](docs/methodology.md).

Matched geometry is kept as drawn (translated only). Everything else is placed by the Composer (`layout_style=curated` by default; `--layout-style legacy` is CLI-only for orthogonal libavoid).

### Layout stages

`tool_pipeline.draw_map` / `draw_global_map`. Small models (≤150 reactions) take the full path; mid-size (≤3 000) use a global plan then the same curated Bézier pass; larger models keep a skeleton and draw ~50-reaction views on demand.

| Stage | What it does | Main code |
|---|---|---|
| **Rigid curated skeleton** | Scale matched structures; pack rigid blocks with shared metabolites as junctions; network bridges as anchors | `drawing.py`, `geometry.py` |
| **Complete uncovered reactions** | Place free nodes around fixed anchors; optional context for currency-only reactions | `completion.py`, `context.py`, `completion_context.py` |
| **Metabolite instances** | Extra visual copies of non-currency metabolites to shorten long edges; currency stays private `side:` glyphs | `instances.py`, `drawing.py` |
| **Orientation** | Put free reactions between substrates and products; currency copies on the correct side | `orientation.py` |
| **Composition & peripherals** | Pull related components in; shelf independents; attach exchange/transport; optional relocation inside large curated blocks | `composition.py`, `peripheral.py`, `block_relocation.py` |
| **Spacing / relaxation** | Soft-preserve curated shape; shorten edges; minimum node spacing | `geometry.py`, `global_layout.py` |
| **Bézier connections** | Cubic Béziers from reaction templates; optional `--source-paths` for original Escher curves | `edge_style.py`, `source_paths.py`, `routing.py` |
| **Labels** | Move labels among candidates to cut collisions | `label_placement.py` |
| **Integrity** | Model hash, every reaction drawn, visual stoichiometry; protected shared glyphs must not split | `tool_pipeline.py` |

**Determinism:** fixed layout seed `20261105` (`core.SEEDS`), ELK/fCoSE seeds, `PYTHONHASHSEED=0` on re-exec. Displacement measured by `evaluate --compare-map` only.

### Example output

Maps composed from the bundled example models. Every reaction and metabolite is drawn.

<p align="center">
<img src="docs/images/examples/e_coli_core.png" alt="Escherichia coli core map (e_coli_core, 95 reactions)" width="700"><br>
<em>Escherichia coli</em> core (<code>e_coli_core</code>, 95 reactions)
</p>

<p align="center">
<img src="docs/images/examples/iIT341.png" alt="Helicobacter pylori iIT341 map (554 reactions)" width="700"><br>
<em>Helicobacter pylori</em> iIT341 (554 reactions)
</p>

*Previews are cropped to content, 2800 px wide. Labels are kept only when they stay readable at that width. Regenerate: `PYTHONPATH=layout_atlas python3 tools/render_map_preview.py DIR/map.json docs/images/examples/NAME.png`.*

---

## 3. Evaluation criteria

Compose writes `layout_quality.evaluate` into `layout_issues.json` / `report.json`. `python3 project.py evaluate` adds layout-rules metrics via `evaluation/score.py`. See [docs/layout_evaluation.md](docs/layout_evaluation.md).

### Metrics

| Group | Metrics | Code |
|---|---|---|
| **Layout-benchmark metrics** | Crossings and crossing density, crossing angular resolution, overlapping pairs and overlap area share, participant consistency, edge length deviation, edge orthogonality, angular resolution, aspect ratio, area per node, homogeneity, stress, coverage | `metrics.py`, `escher_bridge.py` |
| **Pathway-scorer rules** | Longest line over median edge length, share of lines longer than 2.5× the median, bends per line, cofactor clones (currency hubs and long currency arcs) | `map_rules.py` |
| **Map rules** | Cycle circularity, direction and order; axis-aligned main segments and reaction axes; metabolites on a grid; chain turns; main flow direction; exchange reactions on the periphery; blank space and largest empty rectangle; unrelated overlaps (node overlaps, lines through nodes, collinear lines); label overlaps; species outside their compartment | `map_rules.py` |
| **Curated-structure integrity** | Structure continuity and split structures: curated structures stay whole and are only translated | `continuity.py` |
| **Stability** | Reproducibility across runs; displacement and Procrustes disparity against a previous map (`--previous`) | `stability.py` |
| **Reference similarity** | Procrustes, distance correlation and kNN preservation against a reference SBGN map | `reference_sbgn.py` |

### Acceptance criteria for changes

A change to the layout pipeline is merged only when, on `e_coli_core` and `iIT341`:

- no layout-benchmark metric gets worse, measured with the official layout-benchmark scorer;
- every reaction and metabolite is drawn;
- curated structures keep their internal geometry (translation only) and `--previous` pinned nodes do not move;
- the canvas does not grow and no label is cropped;
- compose time does not increase by more than 20%.

### Optimisation targets

| Target | Optimising stage(s) |
|---|---|
| Connection crossings; unrelated overlaps | Bézier routing · quality gate; spacing / relaxation |
| Detour, hooks, bends | Bézier routing |
| Long lines; cofactor clones | Instances · composition · peripherals |
| Participant placement | Orientation |
| Curated-structure continuity | Matching / rigid skeleton |
| Exchange on periphery; aspect; blank space | Composition & packing |
| Axis-aligned segments / axes; grid metabolites; chain turns; cycle placement; compartment packing | **Not yet optimised** |
| Build-on-previous / Procrustes displacement | `--previous` pin/graft path; measure with `evaluate --compare-map` / layout-benchmark |

Layout targets (including backlog) stay on Detailed metrics; non-targets under Debug (including main flow-direction share).

**Known gaps** ([docs/limitations.md](docs/limitations.md)): 2D component packing, harder orientation (incl. transporters), compartment-aware packing, stronger currency/biomass handling, crossing-aware composition acceptance, compose seed flag. Experimental off by default: `OPM_REGIONS`, `OPM_LOWCONF`, whole-block `OPM_BLOCK_RELOC`.

### Layout rules (status)

| Rule | Status | Notes |
|---|---|---|
| Avoid crossings | Partial | Minimised in Bézier routing / quality gate; GEMs still show many algo crossings |
| Avoid overlaps | Partial | Spacing + label pass + gate; overlaps grow with map size |
| Avoid long lines | Partial | Instances, composition, peripherals help; currency/biomass spokes remain |
| Place participants consistently | Partial | Orientation on by default; share still below curated reference maps |
| Respect compartments | Not yet | Used in packing/shelves; cross-compartment species metric still high |
| Aspect / compactness | Partial | Shelf aspect target; independent components still linear shelves |
| Reproducibility | Partial | Fixed seeds; no compose seed CLI |
| Build on previous layout | Partial | `--previous map.json` reuses matched geometry and grafts leftovers |
| Scale to GEMs | Partial | Global ≤3 000 + tiled views; large maps start routing-degraded |

---

## Outputs

| Output | Description |
|---|---|
| `map.json` | Editable Escher map |
| `map.svg` | Vector visualisation |
| `report.json` | Coverage and layout-quality measures |
| `reaction_provenance.json` | Per-reaction provenance and matching evidence |
| `layout_issues.json` | Crossings, collisions and other layout problems |
| `licenses.json` | Licensing metadata of Atlas structures used |

## Licensing

**Code** is [MIT](LICENSE). **Atlas data and generated maps** mix CC0, CC BY 4.0 and MIT sources; keep the attributions in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Documentation

- [Installation](docs/installation.md)
- [Data sources](docs/data_sources.md) · [Atlas data management](docs/atlas_data_management.md) · [Atlas layout](data/README.md) · [iCH360 import](docs/ich360_import.md) · [Optional packages](docs/data_packages.md)
- [Reaction–map matching](docs/reaction_map_matching.md) · [EC extension](docs/ec_extension_recommendation.md)
- [Methodology](docs/methodology.md) · [Limitations](docs/limitations.md) · [Layout evaluation](docs/layout_evaluation.md)
- [Code overview](layout_atlas/README.md) · [About / Provenance](layout_atlas/results/about/index.html) (regenerate with `python3 layout_atlas/navigation.py`)

## Project status

Under active development. Focus: global pathway organisation, less visual overlap, consistent layouts across GEMs.

## Citation and acknowledgements

A citation will be added with the first release. Until then, cite this repository and the sources in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

Developed during ELIXIR BioHackathon Europe 2026, building on the Escher, WikiPathways, Reactome and wider metabolic-modelling communities.
