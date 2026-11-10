# Layout atlas v0.3 — Traceable Database, Skeleton, Multiple Candidates and Conditional Completion

The current version is completely independent of the old directories; original maps, source observations and licenses are in `data/source_inputs/`, and the Node backend is in `runtime/`. The old version has been deleted.

Unified entry: use `python3 project.py --help` from the project root.

## Directly viewable

- `results/catalog/index.html`: a 77-page catalog of representative figures of real instances; `modules.json.gz` is the complete JSON exchange file.
- `data/atlas.sqlite`: the `modules/core_family/variant` layout index and the `reactions(namespace,identifier)` biochemical candidate index; compressed payloads are read with `database.unpack`.
- `results/models/<model>/<condition>/index.html`: fixed-range subgraphs of all reactions; `skeleton.json` records module mappings, sources, coordinates, confidence, actual target ports, bridges and reactions pending layout.
- `results/ranking_queries.json.gz`, `completion_tasks.json`, `completion_results.json`: per-query, per-task and per-seed results.
- `results/verification.json`, `test_run.log`: automatic integrity checks; these checks do not represent manual acceptance or qualified readability.

Current acceptance checks the SHA256 of the required source files through `data/input_manifest.json`.

## Environment and commands

Python 3.11+, CPU, Node.js, clang++; exact Python versions are in `requirements.txt`. The Node version and the pinned dependencies of `elkjs/cytoscape/fcose` follow `runtime/package-lock.json`. First-time installation:

```sh
python3 -m pip install -r layout_atlas/requirements.txt
npm ci --prefix layout_atlas/runtime
```

Run from the project root. The following covers all stages in order and can also be rerun individually; when network sources are cached, SHA256 is checked and nothing is silently updated. Quality acceptance uses automatic checks, source traceability and reproducible experiments, with no two-person review or manual sign-off gate.

```sh
python3 layout_atlas/cli.py all
```

Step-by-step run:

```sh
python3 layout_atlas/cli.py freeze
python3 layout_atlas/cli.py acquire
python3 layout_atlas/cli.py archive
python3 layout_atlas/cli.py routing
python3 layout_atlas/cli.py build
python3 layout_atlas/cli.py catalog
python3 layout_atlas/cli.py ranking
python3 layout_atlas/cli.py tasks
python3 layout_atlas/cli.py learning
python3 layout_atlas/cli.py models
python3 layout_atlas/cli.py deploy
python3 layout_atlas/cli.py verify
python3 layout_atlas/cli.py report
```

Six target models: `e_coli_core`, `iML1515`, `iMM904`, `Human-GEM`, `iYO844`, `AraCore_v2_1`. Examples of single-model / learned-ranking application:

```sh
python3 layout_atlas/cli.py models --model e_coli_core --condition all_geometry
python3 layout_atlas/cli.py models --model e_coli_core --condition learned_geometry_community --ranker-source community
```

The learned-ranking parameter specifies the **held-out source fold**; the file is a trusted pickle produced by local training, and arbitrary pickles supplied by others should not be loaded. Production model application has no independent expert answers and cannot be used as a held-out learning score. Ranking learning is first evaluated on a geometric candidate pool that does not use target-expert coordinates; biochemical evidence levels are retained separately.

## Data and interface boundaries

`modules` contains complete sources, original instance pointers, core nodes, appendages, observed/derived coordinates, paths, reactions and ports; geometric families ignore concrete chemical identities and side-metabolite configuration, so chemical equivalence cannot be inferred from them. Rhea/ChEBI identities are stored as multiple candidates with evidence and are not compressed into globally unique entities. Derived GPML reaction centers and pathway-corridor midpoints are not used as expert coordinate labels.

Each reaction in the skeleton has the state `expert_module/algorithmic_bridge/pending`; for the first, `match_level` must still be consulted to tell whether it is a strict equation, reaction-identifier support, or geometric transfer. Export completes all target participants, and source templates do not add reactions that the target lacks. Bridges are algorithmic layout. Shared boundary metabolites keep module-local copies and derived connection points. Compartment blocks are grouped by the compartment combination of the target reactions, and cross-compartment transport is treated as a composite container; true membrane boundaries have not yet been generated.

The completion task `input` contains only the typed graph, hint coordinates, hard anchors, confidence and the anchor coordinate frame; `labels` are for training/evaluation only. The three anchor types for the same query and masking use exactly the same visible/hidden nodes. The EGNN prediction function can directly accept tasks without `labels`. Anchors in expert/perturbation tasks are hard-fixed; retrieval anchors are soft constraints. Compartment geometric boundaries that are not given cannot be filled in from test answers.

libavoid uses the pinned-commit source and original license, and performs orthogonal routing with node obstacles only; label rectangles are not treated as obstacles, so it is not claimed to solve all label problems. Its build log and commands are in `results/routing_build.json`.

`config.json` records the fixed parameters actually used this time; random seeds 20261105/06/07. Controlled budgets include a nominal 120 generic query roots (sampled by stride, actual count may be slightly more), at most 600 additional windows derived from reaction identifiers, and at most 24 drawings per query; unsearched regions are marked separately and cannot be interpreted as absolute absence from the database.



## MicroMap Master-Map Extension Survey

The pinned-version original files and license are in `data/references/micromap/`, independent of the open core; the current core library, skeleton and learning results do not mix in these data.

```sh
python3 layout_atlas/master_map_audit.py
```

Outputs `results/master_map_audit.json` and 12 research fragments. Complete equation candidates are not verified transferable coverage; colored maps derived from the master map do not count as independent expert drawings.

## P0: Drawing Geometry Improvement Experiment

```sh
python3 project.py geometry --model e_coli_core
```

With the existing 11 modules and 95 target reactions fixed, runs in order the original-implementation replay, physical scale, joint solving of graphical copies, and label-rectangle avoidance; saves per-condition results for the full graph and three fixed ranges. Output is in `results/geometry_p0/e_coli_core/`. The new interface explicitly distinguishes soft and hard anchors and unifies graphical coordinates. Hard fixes are kept, and soft hints use a finite per-axis movement range; label routing uses a separately compiled backend from the same pinned-commit libavoid source.

The current command limits the complete target to no more than 150 reactions, validating E. coli core first; the original `models` and existing training results were not replaced. This experiment adds no new templates or learning evidence. Canvas area, remaining label conflicts, the solver convergence flag and routing degradation are all retained.

## P1: Independent MicroMap Research Extension Library

```sh
python3 project.py micromap
python3 project.py micromap --model e_coli_core
```

The command builds an independent extension pack, then compares the current library with the library plus MicroMap using the same annotation/equation mapper. It records the reaction denominators of all six models, candidate-set differences, per-condition selections and undrawn status; E. coli core is exported in full, and iML1515/iYO844 use pre-fixed first/middle/last ranges for the geometry control. These subgraphs of the large models do not mean that new drawing has been run on the whole model.

The extension data are in `data/extensions/micromap/` and explicitly inherit CC BY-NC-SA 4.0; the original open core and the training set are not rewritten. `SOURCE.json` stores the dataset authors, DOI, version, license and change notes; the composed outputs' JSON/SVG/sidecar manifests keep each module's source and license.

[Results and catalog](results/micromap_extension/index.html). The first version decodes only non-branched centerlines and the derived centers on source edge segments; complex branches and full boundary clipping are not implemented; when multiple copies must be collapsed into the same target core node, the candidate is explicitly rejected rather than silently merged.

## P2: Module Composition Search in Context

```sh
python3 project.py compose
python3 project.py compose --model e_coli_core
python3 layout_atlas/composition_replay.py
```

With the P1 source libraries and candidate pool fixed, preserves the original specific coverage set and equation-compatible set, and compares real drawings, allowed rotations, and bounded coverage swaps/splits. Uses the actual interfaces of the target network, block previews and at most 8 full-layout selections; no target-expert coordinates are input, and no new learner is claimed to have been trained.

Results are in `results/composition/<model>/<range>/`: `seed/` is the original-selection replay, `trials/` retains accepts/rejects and failures, and `best/` is the result actually adopted. `preview_log.json` and `summary.json` record the budget, parent candidates, costs and coverage; source licenses and algorithmic poses go into the export and do not increase the expert-instance count.

[Entry to all ranges](results/composition/index.html). Large models use the three sub-ranges pre-fixed in P1, and their metrics cannot represent the whole graph. After rebuilding P1 candidates or changing the relevant code, the composition experiment must be rerun; verification rejects old results that used changed inputs.

`python3 project.py all` now refreshes the P0, MicroMap and P2 experiments in turn after rebuilding the core results, then runs unified acceptance, so that old experiment snapshots are not used after inputs are updated. This command includes data acquisition and all experiments, and takes significantly longer than viewing saved results.
