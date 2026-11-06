# LEfSe analysis

Project-aware DNA/RNA comparisons and the corresponding bubble figure.
`run_lefse.R` recalculates the statistics from TPM, taxonomy and metadata.

## Files

| Path | Contents |
| --- | --- |
| `inputs/` | Six DNA/RNA TPM, taxonomy and sample-metadata tables |
| `run_lefse.R` | Statistical analysis, Source Data export and figure generation |
| `results/source_data.tsv` | Complete comparison-by-feature results and figure data |
| `results/audit_tables/` | Comparison definitions, sample membership and detailed calculations |
| `figures/` | DNA/RNA bubble figure in PDF, SVG, PNG and TIFF formats |

## Run

Use R 4.4.1 and the dependencies in [../environment/README.md](../environment/README.md),
including lefser 1.16.2. Arial and Cairo PDF support are required for the supplied figure layout.
From this directory, run:

```sh
Rscript run_lefse.R inputs reproduced
```

Tables and figures are written to `reproduced/results/` and `reproduced/figures/`.
Use a new output folder to retain the supplied results; existing files with matching names are overwritten.
`results/source_data.tsv` is a generated output, not an input to the statistical calculation.

## Analysis notes

Comparisons are stratified by analysis project, ocean basin, coral host and sample compartment,
with at least three sample accessions per group. DNA and RNA are analysed separately.
For DSS, only the original labels `Healthy/Normal` and `Disease(Dark Spot Syndrome)`
(after trimming whitespace) are included; nutrient-enrichment variants are not included in that comparison.

Taxonomic features include formal families and higher-rank unclassified bins.
The historical `CPM` columns contain TPM-derived compositions rescaled to 1,000,000 per sample.
Each comparison uses seed 123, unadjusted KW P < 0.05 and |LDA| ≥ 2 for marker selection;
the figure further requires |LDA| ≥ 2.5. Within-comparison BH-adjusted values are reported separately
and do not determine selection or display. Missing statistics remain missing, rather than being replaced by zero.
`Display_Flag` identifies the displayed results; `Healthy_Samples` and `Condition_Samples` record their samples.
