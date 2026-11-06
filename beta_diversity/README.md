# Beta diversity and PCoA

DNA/RNA virome beta diversity, PCoA and ecological-factor comparisons.
The analysis uses square-root-transformed TPM and Bray–Curtis distances.

## Files

| Path | Contents |
| --- | --- |
| `inputs/` | DNA/RNA TPM matrices, sample metadata and `permutation_plan.csv` |
| `nested/run_pcoa.R` | Nested comparisons and PCoA figures |
| `factor_evaluation/run_factor_evaluation.R` | Global, project-conditioned and project-specific factor comparisons |
| `nested/Data/`, `factor_evaluation/Data/` | Statistical results and analysis objects |
| `nested/Figures/`, `factor_evaluation/Figures/` | Exported figures |
| `nested/Source_Data/`, `factor_evaluation/Source_Data/` | Figure coordinates, statistics and source-table descriptions |

## Run

Use R 4.4.1 and the dependencies in [../environment/README.md](../environment/README.md).
Run from this directory, choosing new or empty output folders:

```sh
Rscript nested/run_pcoa.R --input-dir inputs --output-dir reproduced/nested --mode all
Rscript factor_evaluation/run_factor_evaluation.R --input-dir inputs --output-dir reproduced/factor_evaluation --mode all
```

Each command recalculates statistics, draws figures and exports Source Data.
Outputs use the `Data/`, `Figures/` and `Source_Data/` structure beneath the selected output folder.
Both scripts default to four workers; add `--workers 1` for serial execution.
Use `--help` to see the additional analysis, figure-only and source-data modes.

The supplied permutation plan records seeds and assignment checksums, not precomputed test results.
Each test rebuilds its assignments and calculates its statistics: normally 9,999 permutations,
or all 719 non-identity assignments for the relevant finite restricted designs.
Project is a study-level block or conditioning variable, not an additional ecological factor.
The result and Source Data tables retain the test context, sample counts, P values and adjusted P values.
