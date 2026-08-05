# ComBat-refQL 0.4.0 validation plan

## Simulation design

The factorial simulation fixes seeds 104729 and 130363, 120 genes, two batches,
and two biological groups. It crosses:

- no, sparse, and widespread mean batch effects;
- source/reference dispersion ratios 0.5, 1, and 2;
- 2, 3, 5, and 10 replicates per biological stratum;
- balanced and unbalanced groups.

Alternating scenarios introduce partial biological overlap and a deliberately
weak reference. Every scenario includes heterogeneous library depths, low- and
high-expression genes, heterogeneous gene dispersion, correlated modules, and
a strong outlier gene.

The study compares the single ComBat-refQL 0.4.0 method, raw-count edgeR design
modelling, pinned ComBat-ref, and pinned ComBat-seq. Comparator code remains in
validation and is not part of the production API.

## Metrics

Known-truth metrics are bias and RMSE of mean batch contrasts and dispersion
shifts, dispersion-ratio recovery, biological log-fold-change RMSE,
transformation magnitude under null and non-null effects, conditional mean and
variance alignment, exact reference invariance, runtime, and raw-count edgeR
FDR and power. Adjusted counts are not used as an inference replacement.

Acceptance requires deterministic completion without mapping errors, exact
reference preservation, lower transformation under the no-batch-effect
scenario than under widespread effects, finite recovery metrics, and no
unsupported performance claim.

## Real-data design

The GFRN and NASA GeneLab workflows reproduce the upstream sample selections
and filtering. Quantitative endpoints are conditional batch association in the
leading PCs, mean alignment, residual-dispersion alignment, and reference
invariance. Real data characterize empirical behavior and do not establish
truth-based recovery.
