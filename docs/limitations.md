# Limitations

- Learning gains over the baseline are not yet consistent across models.
- Full maps can still show weak global pathway organisation and overlaps.
- Large models may use tiled or sub-range drawing; those metrics do not represent a single whole-graph layout.
- The default Atlas holds curated community maps from Escher, Human-maps, WikiPathways, Reactome, BRENDA, Metabolism Regulation, PANTHER and iCH360; coverage for organisms or pathways absent from these sources is mostly algorithmic.
- The first setup needs network access and a C++ compiler (`g++` or `clang++`). Validated on macOS; Linux is expected to work but not CI-tested; Windows is untested (WSL recommended).
