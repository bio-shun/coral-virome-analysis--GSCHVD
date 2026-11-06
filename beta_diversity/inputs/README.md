# Beta-diversity inputs

These are the fixed processed inputs used by both beta-diversity workflows. The DNA feature set contains 24,083 viral features. The scripts start from TPM matrices; they do not align reads, estimate counts or recompute TPM normalization.

| File | Contents and orientation |
| --- | --- |
| `24083_metagenome_TPM.tsv` | DNA TPM; one viral feature per row and sample libraries in columns; the first column contains feature identifiers. |
| `rna_virus_tpm.tsv` | RNA-virus TPM in the same orientation. |
| `sample_info_metagenome.tsv` | DNA sample metadata; one sample library per row. |
| `sample_info_metatranscriptome.tsv` | RNA sample metadata; one sample library per row. |
| `permutation_plan.csv` | Fixed publication permutation seeds, RNG settings and assignment hashes for 292 distinct tests. |

Required metadata columns are `Sample`, `Project`, `Ocean_Basin`, `Coral_Host_Species`, `Sample_Type` and `Status`. They supply the sample key, study identifier, sampling region, coral taxon, compartment and health status. Sample identifiers must be unique and match the matrix column names. `Project` is a study/blocking variable; a sample library is not automatically an independent biological replicate.

The workflows validate finite nonnegative TPM, remove features with zero abundance across the input samples, match sample IDs, remove zero-total samples and filter missing/unknown study or ecological labels. Label-normalization rules are explicit in `clean_species()`, `clean_compartment()` and `clean_status()`. They then take the square root of TPM and calculate Bray-Curtis distances. The resulting DNA analysis has 487 sample libraries, 15,240 retained features and 22 Projects; the RNA analysis has 744 sample libraries, 460 retained features and 44 Projects.

The permutation plan contains no stored test statistics. `Stable_key` identifies the exact data/model/permutation design; `Original_test_ID` links it to the archived analysis record. `Virus`, `Type` and `N` describe the test. `Requested_permutations` and `Actual_permutations` distinguish the 9,999-permutation request from finite-space enumeration. `Seed`, `RNG_kind`, `Normal_kind` and `Sample_kind` reconstruct the random assignments; `Permutation_SHA256` verifies the regenerated assignment matrix before model fitting. There are 262 tests with 9,999 assignments and 30 with 719 complete nonidentity assignments. The plan is specific to this fixed input snapshot. A new sample set or model requires a reviewed new plan.

Module checksums cover all five input files. Keep these files together and provide their directory through `--input-dir`; no personal filesystem path is required.
