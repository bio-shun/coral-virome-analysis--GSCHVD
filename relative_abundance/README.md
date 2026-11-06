# Taxonomic abundance heatmaps

This module contains independent DNA virome, RNA virome and bacterial MAG family heatmaps. Each script reads the bundled inputs and exports the exact plotted source-data matrix and a PDF. Values are **ln(1 + mean TPM)**, not percentage relative abundances. The original aggregation, filtering, ordering, labels and colours are preserved.

Use the versions in `../environment/README.md`. Run from this directory:

```sh
Rscript dna/run_heatmap.R dna/inputs reproduced/dna
Rscript rna/run_heatmap.R rna/inputs reproduced/rna
Rscript bacterial_family/run_heatmap.R bacterial_family/inputs reproduced/bacterial_family
```

Both arguments are explicit input/output directories. Use another output directory for an independent rerun. No machine-specific location is required.

For DNA and RNA, the first four TSV columns are `Row_Order`, `Row_ID`, `Phylum` and `Label`; remaining columns use complete group identifiers. Numeric zeros are retained, missing values are `NA`, and values are exported without analytical rounding. The RNA script retains the project-specific Avila grouping used in the supplied analysis. See the bacterial submodule README for its distinct averaging definition.

The DNA input contains the supplied 24,083-MAG TPM matrix. RNA and bacterial inputs are separate datasets; their inclusion does not imply that they were recalculated from that DNA matrix. Bundled PDFs preserve the original large-format geometry and have not been redesigned for a journal's final print dimensions.
