# Phage lifestyle composition

Phage lifestyle composition table and figure. Use the environment described in `../environment/README.md`.

## Files

- `inputs/24083_metagenome_TPM.tsv`: viral TPM matrix.
- `inputs/lyso.tsv`: genome lifestyle annotations.
- `inputs/sample_info_metagenome.tsv`: sample metadata.
- `run_lyso.R`: calculation and plotting.
- `results/` and `figures/`: composition table, audit object, PDF and preview.

## Run

From this folder:

```sh
Rscript run_lyso.R inputs reproduced
```

Genomes are restricted to the lifestyle table. Annotations containing `Lysogenic` are labelled Lysogenic; remaining annotations are labelled Unknown. DSS samples and unrecognized compartments are excluded. Each sample is normalized by total TPM of the target genomes; displayed groups are arithmetic means of sample-level percentages. Project is not a grouping variable.
