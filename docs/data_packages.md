# Optional data packages

The repository contains the **open-core Atlas** (`layout_atlas/data/open_core.sqlite`: 17,221 modules and 63,274 reactions from WikiPathways, Reactome and the Escher community maps). The full structure database used during development contains a few more parts. They have their own licences and are therefore published as **separate downloads**, each with a notice inside:

| Package | Content | Terms |
|---|---|---|
| `atlas-full` | the full `atlas.sqlite` (same modules and reactions as the open core; the `maps` table additionally lists the Human-maps and Yeast-maps sources) | see the package notice |
| `micromap` | MicroMap extension and its reference files | CC BY-NC-SA 4.0 (non-commercial, share-alike), MicroMap v4.0, I. Thiele, doi:10.7910/DVN/FZKMJ8 |
| `brenda` | BRENDA pathway-map extension (trial) | layout origin not confirmed: research use only |
| `reference-maps` | Human-maps, Yeast-maps, original yeast map, rule history, with provenance records | per-source terms differ and are not all confirmed; see `docs/licensing_audit.md` |

```sh
python3 project.py data list                 # what exists, size, terms, installed or not
python3 project.py data fetch micromap       # download, verify the SHA-256, install into layout_atlas/data/
python3 project.py data fetch all
python3 project.py data fetch all --from DIR # install from a local folder or another base URL
```

- Each archive is checked against `layout_atlas/data/packages.json` before anything is installed. An existing file is kept as `<name>.bak`; symlinks are never touched.
- With every package installed, `layout_atlas/data/` has exactly the same files as the developer database, and the example map generated from a fresh clone is identical.
- The default tool run (no MicroMap / BRENDA options) does not need any package.
- Respect each package's terms (see `package_notices/<name>.md` after installing).
