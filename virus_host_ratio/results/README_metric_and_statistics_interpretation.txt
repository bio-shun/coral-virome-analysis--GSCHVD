Metric:
Aggregate linkage-informed virus-to-host abundance ratio =
sum abundance of all viruses represented in the virus-host linkage table /
sum abundance of all hosts represented in the virus-host linkage table.

This is NOT a pair-specific virus-to-host ratio.

Statistical design:
Comparisons are performed only within the same Analysis_Project, Ocean_Basin,
Species and Sample_Type. Each retained status group contains at least 3
nonzero samples, including a Healthy/Normal reference.

Primary statistical test:
two-sided Wilcoxon rank-sum test on log10-transformed aggregate ratios.

Primary multiple-testing correction:
Benjamini-Hochberg adjustment across Healthy/Normal-vs-other-status
comparisons within each Analysis_Project x Ocean_Basin x Species x Sample_Type
context.

A global BH-adjusted P value is also exported as a supplementary reference.

Abundance-scale caveat:
If virus_abundance.tsv and host_abundance.tsv were independently normalized to separate
compositional totals, interpret this metric as a relative viral-to-host signal
ratio rather than an absolute virus-to-host abundance ratio.
