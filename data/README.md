# Data

Everything in this directory is read by the tool at run time. The Atlas (the structure database) is built from existing expert-drawn maps; sources, versions and licences are in [`THIRD_PARTY_NOTICES.md`](../THIRD_PARTY_NOTICES.md) and [`docs/licensing_audit.md`](../docs/licensing_audit.md).

## In this repository

| Path | Content |
|---|---|
| `open_core.sqlite` | the open-core Atlas: 17,221 modules and 63,274 reactions (WikiPathways, Reactome, Escher community maps). Modules store core nodes, appendages, coordinates, paths, reactions and ports as compressed JSON (read with `database.unpack`). On first run it is copied to `atlas.sqlite`, which is the file the tool reads |
| `source_inputs/` | source data used to build the Atlas: `raw/community/` (original Escher community map JSON), `raw/wp/` (pinned WikiPathways GPML), `module_observations.json.gz` and `map_observations.json.gz` (observations with original identifiers, families, geometry and provenance chains), `archives/` and `terms/` (source packages and licence snapshots), `origin.json` and `source_registry.json` (versions and hashes) |
| `input_manifest.json`, `acquisition.json` | SHA-256 of the required inputs and how they were obtained |
| `rhea/` | Rhea release 143 reactions, directions and identity mappings (CC BY 4.0, copyright statement kept) |
| `reactome/` | Reactome v97 pathway diagrams and graphs |
| `models/` | additional target models (AraCore, iYO844) with their source records |
| `extensions/` | data derived from the open-core maps: `units/` (structure units), `source_paths/` (expert connection geometry) and `reaction_geometry_templates/` (connection-style templates used for the Bézier routing) |
| `routing/` | pinned libavoid source, commit and licence |
| `bibliography/` | evidence list used by the literature check |
| `packages.json` | manifest of the optional data packages |

The number of module instances or source files is not the number of independent expert drawings: several modules can come from the same map.

## Not in this repository: optional data packages

The full database used in development also has the MicroMap and BRENDA extensions, reference maps and a full `atlas.sqlite` (same modules and reactions; its `maps` table also lists the Human-maps and Yeast-maps sources). Their licences differ, so they are separate downloads:

```sh
python3 project.py data list
python3 project.py data fetch micromap     # or: brenda, reference-maps, atlas-full, all
```

Details and terms: [`docs/data_packages.md`](../docs/data_packages.md). A normal run needs none of them.
