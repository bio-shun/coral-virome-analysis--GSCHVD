# Figure Source Data

master_source_data.tsv is the combined result table: one original source-table cell per row. raw_value preserves the source text (including numerical precision); numeric_value is populated only for numerical literals. Missing values retain the original source notation. source_row is 1-based and excludes the header.

Combined TSV files use standard doubled-quote escaping. A standard CSV/TSV reader recovers embedded quotation marks in status annotations. Import identifier and raw_value columns as text to preserve their original spelling and precision in spreadsheet software.

source_table_index.tsv maps each module/panel to its original table, dimensions, units/measurement context and SHA-256. tables/ contains complete figure-specific source tables. field_dictionary.tsv lists all original fields; module READMEs and the Beta dictionaries define statistical methods and exclusions.

statistical_summary.tsv gathers reported P-value fields with module/panel/row provenance. It is an index of reported tests, not a new family of hypothesis tests; repeated values across figure summaries and detailed panel tables are intentional. Adjustments apply only to the families defined by each analysis. Selection flags in the original LEfSe/VHR/network tables must be respected.

No across-module mean, extra correction or unit conversion is applied. Figure values, statistics and annotations keep their original scientific meanings. Manuscript figure numbers should be assigned after the final figure order is set.

Regenerate with: Rscript scripts/build_source_data.R <analysis_root> <source_data_directory>
