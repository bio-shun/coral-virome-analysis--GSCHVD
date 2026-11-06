# Viral network analysis

Viral co-occurrence networks, bootstrap stability, topology comparisons and the composite network figure.

## Files

| Path | Contents |
| --- | --- |
| `input_data/` | TPM, taxonomy, sample metadata and network display coordinates |
| `run_network.R` | Network calculation, statistical tests and figure export |
| `results/network_analysis_master_results.tsv` | Consolidated network estimates and statistical results |
| `figures/` | Composite figure in PDF, SVG, PNG and TIFF formats |
| `figure_source_data/` | Panel-level node, edge, bootstrap and statistical tables |
| `process_data/` | Intermediate network calculation tables |
| `provenance/` | LEfSe annotation source and analysis records |

## Run

Use R 4.4.1 and the dependencies in [../environment/README.md](../environment/README.md).
Run from this directory. First calculate LEfSe, then explicitly supply that run's source table:

```sh
Rscript ../lefse_analysis/run_lefse.R ../lefse_analysis/inputs reproduced/lefse
Rscript run_network.R full input_data reproduced/network reproduced/lefse/results/source_data.tsv
```

Use new output folders to retain the supplied results.
Network outputs are written under `reproduced/network/` with the folder structure listed above.
The Microviridae LDA annotation is taken from the matching DNA/SCTLD LEfSe comparison
and checked against the network's sample membership; it is not used to construct network edges.
The source LEfSe table and selected annotation row are saved under the new output's `provenance/lefse/`.

## Analysis notes

The supplied cohort and display coordinates are retained. Seed 20260813 controls the analysis,
which uses 500 site-balanced bootstrap resamples and 999 within-location health-label permutations.
Topology P values are empirical and two-sided, with BH adjustment across eight metrics.
Panel c shows edge-retention percentages and SCTLD-minus-Healthy differences in percentage points;
its dots indicate retention ≥ 50%, not statistical significance.
The figure's displayed families are a subset of the families used for the network calculation.
