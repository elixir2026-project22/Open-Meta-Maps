# PANTHER Pathway pilot (CellDesigner)

Status: **acquisition blocked** from this environment (2026-11-11). No bulk import. No PANTHER modules added to the Atlas.

## Goal

Pilot 5–10 metabolic PANTHER pathways through the generic CellDesigner importer (`layout_atlas/celldesigner.py`) and report whether observed coordinates and reaction centre-line mapping work. Licence target: CC BY 4.0.

## Acquisition attempts

All live `pantherdb.org` / `archive.pantherdb.org` / `ftp.pantherdb.org` endpoints returned **HTTP 403** with a Cloudflare challenge page (`Just a moment…`, ~5–45 KB HTML), including:

- Pathway list and diagram pages
- Guessed SBML URLs such as `/pathway/resource/P00051.sbml`, `/pathway/sbml.jsp?acc=P00051`
- Pathway zip mirrors under archive / data hosts

Automated browsers (Chrome headless, Playwright with system Chrome, curl_cffi browser impersonation, cloudscraper) did not clear the challenge (only `__cf_bm` cookie; body remained the interstitial). Wayback Machine had no archived SBML payloads for the probed resource URLs. BioModels / EuropePMC / Zenodo / OSF / Pathway Commons searches did not yield redistributable PANTHER CellDesigner SBML files.

No authentic PANTHER SBML was saved under `/tmp/panther_pilot/` (only challenge HTML).

## Decode / mapping findings

**Not exercised on PANTHER files** — none were obtained. The generic importer remains in place for a future re-run:

- `celldesigner.parse` / `celldesigner.decode` read CellDesigner SBML species aliases (bounds → glyph centres) and non-branch reaction extensions (`baseReactant` / `baseProduct` + `editPoints`).
- Callers supply `source` / `license` / attribution metadata; there is no MicroMap- or PANTHER-specific branch.

## Outcome

| Item | Result |
|---|---|
| Files downloaded | 0 |
| Modules imported into Atlas | 0 (pilot only; no bulk) |
| Coordinates / reaction mapping on PANTHER | Unknown until SBML can be fetched |
| Recommended follow-up | Re-run the pilot from a network that can pass the pantherdb.org Cloudflare check; decode 5–10 metabolic pathways; decide on default Atlas inclusion only after that report |

## Licence note

PANTHER Pathway content is advertised as CC BY 4.0. It is **not** redistributed here until acquisition and decode are confirmed.
