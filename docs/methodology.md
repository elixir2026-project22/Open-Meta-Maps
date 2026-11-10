# Methodology

1. **Atlas construction**: pathway maps from Escher community maps, WikiPathways and Reactome are cut into structural units (modules) that keep reaction identity (BiGG, Rhea, EC, ChEBI cross-references), geometry, source and licence.
2. **Matching**: reactions of the input GEM are matched to Atlas structures through identifier and annotation bridges.
3. **Skeleton composition**: a set of non-conflicting expert structures is selected and assembled into a skeleton.
4. **Completion**: remaining reactions are placed with computational layout, with bridges to the skeleton.
5. **Routing and evaluation**: edges are routed with libavoid and layout issues are reported.

Every reaction is labelled in `reaction_provenance.json` as expert-derived, topology-only suggestion, functional (EC) suggestion, bridge or generated. Coverage numbers describe geometry reuse, not readability.
