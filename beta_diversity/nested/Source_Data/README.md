# PCoA figure source data

Each factor folder contains points and statistics in equivalent UTF-8 CSV and TSV formats. Choose either format; they contain the same values.
Keys are Figure + Virus + Factor + Stratum (and Sample for points). Panel order/tags restart within each composite and include displayed panels only.
Coordinates, centroids, axis percentages, colours and all numerical statistics are read from completed outputs; neither ordinations nor tests are rerun.
All eligible comparisons are retained. Status panels use nominal pooled PERMANOVA P < 0.05; Figure_panel_included is a display flag, not evidence of BH significance.
Numeric P/Q and R-squared values are not rounded for presentation. R-squared values are fractions; axis percentages are percentages. NA denotes unavailable/not estimable, never zero.
BH families are Virus x Factor, separately for pooled/blocked PERMANOVA and PERMDISP. Blocking restricts permutations within Project and does not condition Project sums of squares.
Group envelopes are 80% t ellipses for n >= 4 or convex hulls for n = 3, with spokes to centroids; they are descriptive geometry, not significance regions.
Source_file_checksums identifies the read-only archived inputs. Supporting holds available per-test permutation audit information. No abundance matrices or native checkpoint binaries are duplicated.

中文：每个因素文件夹含坐标与完整统计表，CSV/TSV内容相同。Status未入图的比较也予以保留；本次仅导出，不重算PCoA或置换检验。
