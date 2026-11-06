# Bacterial MAG family heatmap

Run from this directory:

```sh
Rscript run_heatmap.R inputs results
```

Dependencies: `tidyverse`, `pheatmap`, `RColorBrewer`, and the R-supplied `grid` package. Inputs are `bacterial_mag_tpm.tsv`, `bacterial_mag_taxonomy.tsv` and `sample_metadata.tsv`. A different output directory can be supplied for independent reproduction.

The script retains 51 named families and combines remaining assignments into `Others (incl. Unknown)`. For each family and metadata group, it takes `mean(TPM, na.rm = TRUE)` across all MAG-by-sample cells. This is a mean of those cells, not a family sum or percentage relative abundance. No additional zero-total-sample filter is applied.

Families with a maximum group mean of at least 1 are retained; Others is always retained. Values are transformed with `log1p`. Metadata groups use `Project||Species||Ocean_Basin||Sample_Type||Broad_Status||Detail_Status`; plot labels display `Project_Species`. The fixed colour scale spans 0–15. Original ordering, status merging and colours are retained.

The bundled dataset produces 52 rows by 76 groups from 492 metadata-matched samples. The input has 75 additional TPM sample columns without matching metadata; these are excluded by the specified inner join. The output `Bacteria_Family_Heatmap_Log1p_SourceData.tsv` contains Phylum, Family and the 76 complete group labels, exporting the exact plotted matrix without analytical rounding. `Bacteria_Family_Project_Ocean_Heatmap.pdf` preserves the 24-by-16-inch geometry.
