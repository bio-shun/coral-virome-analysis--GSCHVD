# Virus–host signal ratio

Linkage-informed virus-to-host signal ratios and project-stratified status comparisons.

## Files

| Path | Contents |
| --- | --- |
| `inputs/virus_abundance.tsv` | Viral abundance matrix |
| `inputs/host_abundance.tsv` | Host abundance matrix |
| `inputs/sample_metadata.tsv` | Sample, project, host, compartment and status metadata |
| `inputs/virus_host_links.tsv` | Predicted virus–host pairs |
| `run_vhr.R` | Ratio calculation, statistical tests and figure export |
| `results/source_data.tsv` | Sample values, comparisons and group summaries |
| `results/source_data_dictionary.tsv` | Source Data field definitions |
| `figures/` | PDF, SVG, PNG and TIFF figures |
| `provenance/` | Input manifest, run summary and software environment |

## Run

Use R 4.4.1 and the dependencies in [../environment/README.md](../environment/README.md).
From this directory, run:

```sh
Rscript run_vhr.R inputs reproduced
```

Outputs are written to `reproduced/results/`, `reproduced/figures/` and `reproduced/provenance/`.
Use a new output folder to retain the supplied results; existing files with matching names are overwritten.

## Analysis notes

The ratio is the summed signal of unique linked viruses divided by that of unique linked hosts,
counting each entity once even when it has multiple predicted partners.
It is an aggregate signal ratio; interpreting it as an absolute biological virus-to-host ratio
would require comparable upstream quantitative calibration of the two abundance tables.

Only eligible metagenomic samples with positive viral and host signals enter the comparisons.
Each status is compared with Healthy/Normal within the same project, ocean basin, coral species
and sample compartment, with at least three samples per group.
Tests are two-sided Wilcoxon rank-sum tests on log10 ratios.
Primary BH adjustment is within each comparison context; global BH values are supplementary.
Seed 123 fixes point positioning. In `source_data.tsv`, `record_type` distinguishes sample,
comparison and group-summary records; fields that do not apply are `NA`.
