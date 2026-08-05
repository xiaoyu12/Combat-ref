# Validation plan for ComBat-refQL 0.4.0

This plan was fixed before reviewing the 0.4.0 factorial results.

## Design

The factorial simulation uses seeds 104729 and 130363, 120 genes, two batches,
and two biological groups. It crosses no, sparse, and widespread mean batch
effects; source/reference dispersion ratios 0.5, 1, and 2; 2, 3, 5, and 10
replicates per biological stratum; and balanced or unbalanced groups. Alternating
rows introduce partial biological overlap and a deliberately weak reference.
Every scenario includes heterogeneous library depths, low- and high-expression
genes, correlated modules, and a strong outlier gene.

The four formal ComBat-refQL configurations are compared:

| Mean | Dispersion | Label |
|---|---|---|
| none | shared | 0.3.1-style comparator |
| ql_eb | shared | A-only |
| none | batch_scaled | B-only |
| ql_eb | batch_scaled | full 0.4.0 |

Raw-count edgeR design modelling supplies the inference baseline. Full runs also
execute ComBat-ref at commit `f0c6d3f313e9f1f29de2c8c8f12b98113443be55`
and ComBat-seq at commit `2673dd8f23052ea4354ffa5e7600fabcd3f08269`.
Their count-alignment results are labelled separately because their returned
objects do not expose directly comparable batch-contrast uncertainty.

## Metrics and acceptance criteria

Scenario rows record bias and RMSE for batch contrasts and log-dispersion shifts,
estimated dispersion ratios, biology-conditional mean and variance alignment,
biological log-fold-change RMSE, mean absolute transformation, exact reference
invariance, runtime, and edgeR-design FDR and power.

Acceptance requires deterministic reruns, exact reference invariance, no errors
for estimable factorial scenarios, zero transformation for reference samples,
finite diagnostic estimates, and explicit typed errors for rank-deficient or
zero-reference-information designs. Statistical comparisons are descriptive;
no superiority claim follows from this simulation alone.

## Commands

```sh
Rscript validation/simulation.R fast
Rscript validation/simulation.R full
```
