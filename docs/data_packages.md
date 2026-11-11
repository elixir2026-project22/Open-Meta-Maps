# Optional data packages

The repository contains the **open-core Atlas** (`data/open_core.sqlite`: 17,221 modules and 63,274 reactions from WikiPathways, Reactome and the Escher maps). Optional data packages (separate downloads, each with a notice):

| Package | Content | Terms |
|---|---|---|
| `atlas-full` | the full `atlas.sqlite` (open-core Atlas snapshot) | see the package notice |
| `brenda` | BRENDA pathway-map extension (trial) | layout origin not confirmed: research use only |
| ~~`reference-maps`~~ | withdrawn: Yeast-maps / original yeast map dropped (GPL-3.0); Human-maps move into the open-core Atlas via official Escher JSON | — |

```sh
python3 project.py data list                 # what exists, size, terms, installed or not
python3 project.py data fetch brenda         # download, verify the SHA-256, install into data/
python3 project.py data fetch all
python3 project.py data fetch all --from DIR # install from a local folder or another base URL
```

- Each archive is checked against `data/packages.json` before anything is installed. An existing file is kept as `<name>.bak`; symlinks are never touched.
- With every package installed, `data/` has exactly the same files as the developer database, and the example map generated from a fresh clone is identical.
- Respect each package's terms (see `package_notices/<name>.md` after installing).
