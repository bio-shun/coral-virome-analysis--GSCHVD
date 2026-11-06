# Validation report

This report records the completed rerun of the two standalone workflows from the
portable release inputs. It is an audit of computational equivalence, permutation
execution and figure rendering;
it is not a biological validation of the underlying study design.

## Statistical equivalence

The refactored workflows were run from the released inputs and compared with the
completed reference analysis. The comparison covered the global and nested result
tables (tables 01-08, including 07b), filtered distance/metadata objects and
the stored execution audit.

- DNA input: 487 retained samples, 15,240 retained viral features and 22 Projects.
- RNA input: 744 retained samples, 460 retained viral features and 44 Projects.
- Result rows compared: 137 across the primary statistical tables.
- Numeric differences: 0; maximum numeric difference: 0. Reported raw P and BH Q
  values matched exactly in the locked environment.
- Distinct planned design signatures verified: 292. All 305 computed audit
  records across both workflows matched model specifications, original seeds,
  actual permutation counts and assignment hashes.
- Nested PCoA coordinates: exact match across all 54 eligible strata.

All 292 distinct planned designs were computed with a request for 9,999
permutations. The two workflows contain 305 actual executions (198 nested and
107 factor-evaluation); 13 designs are computed independently in both workflows.
Of the executions, 273 use 9,999 permutations and 32 enumerate all 719 available
non-identity restricted permutations. The recorded per-test audit confirms `Cache_origin = Computed`,
the requested count and the actual count. Small strata retain complete finite
permutation spaces when fewer distinct assignments exist. The fixed seeds and
assignment hashes in `inputs/permutation_plan.csv` make the execution traceable.
Matching results were verified with the released inputs, recorded seed plan and
locked R 4.4.1 environment; numerical or binary figure identity across different
software versions/platforms was not established.

## Engine and failure-mode checks

Independent synthetic tests passed for deterministic permutation construction,
blocked and unblocked PERMANOVA, PERMDISP, partial dbRDA, serial versus parallel
execution and corrupted-cache rejection. A corrupted or signature-mismatched
checkpoint fails closed rather than silently falling back to a different test.

## Figure checks

The R figure code was rendered from the completed analysis outputs. The four
nested composite figures and the global four-factor summary rendered
successfully. Panel QA confirmed equal DNA/RNA panel geometry, four-sided panel
borders, shared factor palettes, italicised coral-taxa labels, direct Status
labels and no rendered glyph below the 5-point floor at the native exported
PDF dimensions. This does not guarantee the same text size after reducing a
large composite to manuscript column width; final-size inspection remains
necessary.

## Organised source-data export

The export-only CLI was run against the completed bundle with compute and
render functions trapped. Neither `vegan` nor `ape` was loaded, the archived
input hashes were unchanged, and no statistical test, BH correction or PCoA
was entered. The export contains 54 nested comparisons and 2,037 coordinate
rows (Sampling region 2/530, Coral taxa 8/517, Compartments 8/111, Status
36/879). Status retains all 36 comparisons, with 8 displayed and 28 omitted
from the composite by the existing nominal pooled-P display rule. Global
source tables contain 8 rows each for panels A-C and 41 Project-specific rows
for panel D. CSV/TSV round-trip checks, original-statistic joins, coordinate
sample-set checks and generated-file SHA-256 manifests all passed.

## Scope and limitations

Project is retained as a study-level heterogeneity/blocking variable, not as a
fifth ecological factor. Project-blocked tests only use strata with sufficient
within-Project replication; a non-estimable blocked test is not evidence of no
effect. The workflows preserve the original exchangeability assumptions and do
not establish that all samples are independent biological replicates or that
blocking removes every source of study confounding.
