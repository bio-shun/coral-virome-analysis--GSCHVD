# Source Data: global four-factor evaluation

These files contain saved results and the numeric values used by Figure_Global_Four_Factor_Associations.
The exporter reads completed analysis tables only. It does not fit models, run permutations, correct P values or draw figures.

- Panel_A_global_marginal: all saved marginal PERMANOVA terms and 100*R2 bar heights.
- Panel_B_project_conditioned: all saved Project-conditioned partial dbRDA terms, raw adjusted R2 and display values/labels.
- Panel_C_within_project_coverage: Project eligibility counts, denominators and displayed fractions.
- Panel_D_project_specific: every eligible saved Project-specific result, including nonsignificant results; 100*R2 are the values summarized by the boxes and points.

CSV and TSV are parallel exports of the same table. Numeric P/Q and effect-size columns are not rounded; scientific notation denotes a number, not a threshold. Missing values are written as NA. An intentionally empty Plotted_label (Panel D) means that the renderer has no point-specific text label; it is not a missing statistic.
Panel B preserves negative Adj_R2 and separate Estimable/Model_estimable flags. As in the renderer, the bar height is zero for nonestimable/missing values and is clamped at zero for a negative adjusted R2. The saved raw result is never replaced by that display value.
Panel D display jitter is decorative and is not included in the exported R2-derived values; the source statistic remains exact. Boxplots are generated from these point statistics using the unchanged renderer.
Source_table and Source_row link each exported result to its original saved table. Original statistical columns are retained; the new Figure, Panel, Factor_label, Factor_order, Colour and Plotted_* columns document figure mapping.
Factor_order is top-to-bottom in horizontal Panels A/B and left-to-right in Panels C/D. Internal factor levels are reversed in Panels A/B before coord_flip, matching the renderer. Colour records the exact factor palette.

Saved Benjamini-Hochberg families: Panel A Q was adjusted across the four ecological factors separately within DNA and RNA. Supporting global PERMDISP Dispersion_Q was adjusted across the four factors separately within each viral dataset. Panel D Q and Dispersion_Q were adjusted across eligible Project-specific tests separately within each Virus x Factor combination, with PERMANOVA and PERMDISP corrected as separate families. The partial dbRDA source tables do not contain a Q column; no adjusted value is inferred or calculated by this export.

Supporting/ contains global partial dbRDA, global PERMDISP and sample totals when available. Execution audit and context records are included when available; absolute checkpoint paths are reduced to basenames. Checkpoint payloads and original input matrices are not copied.
An explicit saved test-context map is used to select global audit records. If this map is absent in a legacy bundle, the complete saved audit is retained with an explicit scope note; no ambiguous mapping is inferred.
Missing supporting tables: 
None.

Data dictionary: data_dictionary.csv / data_dictionary.tsv.
SHA256_manifest.csv / .tsv list every generated file except the manifests themselves, with relative paths and SHA-256 hashes.
Raw table provenance is in source_table_checksums.csv / .tsv; it records original filenames and hashes without private directory paths.

Project is the study/blocking variable. The four ecological factors remain Sampling region, Coral taxa, Compartments and Status.
The exporter preserves the existing replication, exchangeability and modelling assumptions; it does not validate biological independence or remove study confounding.
