# Current Data

- `atlas.sqlite`: current modules, layouts, biochemical identity candidates and source index; payloads are compressed JSON, read with `database.unpack`.
- `open_core.sqlite`: the current core release copy exported by acceptance; rebuildable.
- `source_inputs/`: source data required by the current run, now independent of the old version.
  - `raw/wp/`: pinned GPML original maps; `raw/community/`: original JSON of community manually drawn maps.
  - `module_observations.json.gz`, `map_observations.json.gz`: source observations with original identifiers, families, geometry and provenance chains. These are data inputs only, not backups of old-version code or historical outputs.
  - `archives/`, `terms/`: pinned source packages and license evidence.
  - `origin.json`, `source_registry.json`: source versions, hashes and migration information.
- `input_manifest.json`: SHA256 of the required source inputs; `cleanup.json` records deletion of old versions and data migration.
- `rhea/`: pinned release143 reactions, directions and identity mappings.
- `reactome/`: pinned v97 archive members; `archive-v97/` is used in the end.
- `models/`: newly added target models; the rest are read from `Ref/sbml2escher-main/models/`.
- `references/`: reference sources that keep their original licenses; ineligible entries do not enter training.
- `bibliography/`: evidence list cited by the current literature check.
- `routing/`: pinned libavoid source, commit and license.

Do not interpret the instance count or the number of source files directly as the number of independent expert drawings. The whole current pipeline no longer reads the old directories.

## MicroMap Independent Extension

`extensions/micromap/` is a research-extension partition (CC BY-NC-SA 4.0); it is not in `open_core.sqlite` or the current training set. It contains source observations, module exchange JSON, an independent SQLite index and the `SOURCE.json` dataset attribution. The source XML, release checksums and canonical references remain in `references/micromap/`. The master maps and the derived colored maps share the same design origin, so file counts cannot be used to inflate the number of independent experts/species.
