# Phage lifestyle figure specification

The stacked bars display group means of sample-wise Lysogenic/Unknown TPM
percentages among genomes represented in the lifestyle annotation table.
Groups are defined by species, compartment, ocean, health block and status.
The current 24,083-feature input supplies abundance values.

DSS samples, unrecognised compartments and zero-total samples are excluded by
the retained analysis rules. Lifestyle classes, facet order, palette and plot
dimensions follow the supplied workflow. No statistical tests are added.

The figure source table is `results/phage_lifestyle_stats.csv`. The R script
exports a vector PDF, a screen PNG preview and an R-native audit object.
