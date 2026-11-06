# Coral virome analysis

Analysis scripts, processed input data, result tables and figures supporting *Global coral viromes reveal hidden diversity and context-dependent ecological restructuring*.

## Folders

| Folder | Contents |
| --- | --- |
| [alpha_diversity](alpha_diversity/) | DNA/RNA virus Shannon diversity, Chao1 richness, statistical tests and figures |
| [beta_diversity](beta_diversity/) | PCoA, ecological-factor comparisons and permutation tests |
| [relative_abundance](relative_abundance/) | DNA/RNA viral and bacterial-family abundance heatmaps |
| [amg_profiles](amg_profiles/) | Auxiliary metabolic gene abundance summaries and heatmap |
| [lyso_analysis](lyso_analysis/) | Phage lifestyle composition and figure |
| [lefse_analysis](lefse_analysis/) | DNA/RNA virus LEfSe analysis, result tables and figures |
| [viral_network](viral_network/) | Viral association networks, network metrics and panel source tables |
| [virus_host_ratio](virus_host_ratio/) | Virus–host abundance ratios, statistical tests and figures |
| [source_data](source_data/) | Combined results, statistical summary and figure-specific Source Data |
| [environment](environment/) | R/package version records and optional dependency restoration |
| [scripts](scripts/) | Script for assembling the combined Source Data from module results |

## Use

Use **R 4.4.1** and the package versions in [environment/package_versions.tsv](environment/package_versions.tsv); [renv.lock](renv.lock) records the pinned dependencies. The recorded environment is Windows with Arial and Cairo support. Setup details are in [environment/README.md](environment/README.md).

Run the analyses individually using each folder's README. Scripts include comments and explicit input/output paths. The supplied processed inputs are the starting point; corresponding tables and figures are included. Use a separate output directory to retain the bundled results. Run LEfSe before the network analysis and supply its newly generated source table as described in [viral_network/README.md](viral_network/README.md).

## Source Data

In `source_data/`, `master_source_data.tsv` combines the source values, `statistical_summary.tsv` indexes reported statistics, and `tables/` contains figure-specific source tables. `source_table_index.tsv` and `field_dictionary.tsv` describe their provenance and fields. Preserve `raw_value` as text when importing the combined table.

To assemble the summary from the bundled module results, run from the repository root:

```sh
Rscript scripts/build_source_data.R . source_data
```

## Manuscript figures

The authors assemble panels and adjust layout, spacing, lettering and dimensions for the manuscript. These presentation changes preserve the underlying data, plotted values, group assignments, statistical results and scientific interpretation. Accordingly, final manuscript layouts may differ from direct script exports; the scripts and Source Data provide the quantitative basis of the figures.
