# Auxiliary metabolic gene profiles

AMG abundance summaries and heatmap. Use the environment described in `../environment/README.md`.

## Files

- `inputs/`: DNA TPM matrix, sample metadata and AMG annotations.
- `scripts/run_amg_heatmap.R`: calculation and plotting.
- `results/`: plotted `amg_heatmap_log1p.tsv` and gene/group annotations.
- `figures/`: PDF heatmap and PNG preview.

## Run

From this folder:

```sh
Rscript scripts/run_amg_heatmap.R inputs reproduced
```

The output directory receives `results/` and `figures/`. Positive genome TPM records are joined to AMG annotations and summed within sample/category/gene; multiple annotation rows on one genome can contribute separately. This measure uses annotated genome TPM, rather than independently quantified gene expression.

Group means use observed positive sample/gene records; samples without a positive record do not enter the mean. Unobserved group/gene combinations are filled with zero, then displayed as `ln(1 + mean TPM)`. Gene/group order is fixed by the documented annotations and metadata; no clustering is applied.
