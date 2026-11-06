# Alpha diversity

Count-based DNA/RNA Shannon diversity, Chao1 richness, tests and figures. Use the environment described in `../environment/README.md`.

## Files

- `inputs/`: DNA/RNA count matrices, sample metadata and the RNA rarefaction RNG state.
- `scripts/`: calculation and plotting scripts.
- `results/`: source tables, statistical tests, QC and session information.
- `figures/`: PDF/SVG figures and PNG previews.

## Run

From this folder, using a separate output directory:

```sh
Rscript scripts/run_alpha_diversity.R inputs reproduced/results
Rscript scripts/plot_alpha_diversity.R reproduced/results reproduced/figures
```

Samples are matched to metadata, restricted to recognized statuses and filtered at 1,000 total counts before rarefaction. DNA uses seed 123; RNA restores the supplied RNG state to reproduce its original random draw. Both assays are calculated from counts; the RNG file is not a precomputed diversity table.

Within each project/basin/compartment/coral-species stratum, each eligible status requires at least three samples and Healthy must be present. Wilcoxon tests compare Healthy with eligible statuses; BH correction is applied within each stratum and diversity index. See script comments and QC outputs for filtering and grouping details.
