# Methodology

Product overview and figures: [README](../README.md) §1–2 and [About / Provenance](../layout_atlas/results/about/index.html). Atlas cut / provenance detail: [atlas_data_management.md](atlas_data_management.md).

1. **Atlas construction**: pathway maps from Escher, WikiPathways, Reactome, Human-maps, BRENDA, Metabolism Regulation, PANTHER and iCH360 are cut into structural units (modules) that keep reaction identity (BiGG, Rhea, EC, ChEBI cross-references where available), geometry, source and licence.
2. **Matching**: reactions of the input GEM are matched to Atlas structures through identifier and annotation bridges.
3. **Skeleton composition**: a set of non-conflicting curated structures is selected and assembled into a skeleton.
4. **Completion**: remaining reactions are placed with computational layout, with bridges to the skeleton.
5. **Routing and evaluation**: algorithm-placed connections use Bézier curves by default (`layout_style=curated`); `--layout-style legacy` keeps orthogonal libavoid. Layout issues are reported.

Every reaction is labelled in `reaction_provenance.json` as curated (community layout), topology-only suggestion, functional (EC) suggestion, bridge or generated. Coverage numbers measure geometry reuse from Atlas structures.

## Reaction–map matching (EC and name-only nodes)

Matching of model reactions to Atlas enzyme nodes — especially BRENDA (metabolites as names only) and CellDesigner drawings with sparse identifiers — uses Rhea, ExPASy ENZYME, BiGG / identifier bridge, and ChEBI.

| Step | Rule | Default thresholds |
|---|---|---|
| EC currentisation | Transferred EC numbers → successors; deleted dropped (ExPASy `enzyme.dat`) | — |
| EC disambiguation | Among several candidates, score metabolites against Rhea reactions of each EC (weighted Jaccard); keep the clear best | `ec_min_score=0.6`, `ec_min_margin=0.1` |
| Name normalisation | Case, punctuation, membrane sides; stereo / anomer descriptors ignored for matching but kept as a second key; `NAD(P)H` expanded; generic quinones grouped | — |
| Node scoring | Weighted Jaccard of metabolite keys; cofactors / currency at half weight; one-to-one matching | `currency_weight=0.5` |
| Acceptance | EC-linked vs metabolites-only floors; keep nodes within a gap of the best | `min_score=0.6`, `min_score_unlinked=0.8`, `max_gap=0.1` |
| Vicinity | From unambiguous EC matches, propagate through shared non-currency main metabolites | `min_score_vicinity=0.4`, `min_score_vicinity_unlinked=0.6` |

Thresholds are the cross-model defaults validated on eleven curated prokaryote BiGG models in the ecmapper BRENDA-map evaluation; they are not tuned to a single OpenMetaMaps example.

**Method attribution:** the matching *ideas* (EC currentisation, metabolite-weighted EC selection, name keys, weighted Jaccard with half-weight cofactors and one-to-one matching, score/margin acceptance, vicinity propagation) follow [ecmapper](https://github.com/elixir2026-project22/ecmapper) (GPL-3.0). OpenMetaMaps reimplements them independently under its own licence; no ecmapper code is copied, vendored or depended on.
