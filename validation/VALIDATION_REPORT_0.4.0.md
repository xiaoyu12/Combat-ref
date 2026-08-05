# Validation report for ComBat-refQL 0.4.0

## Completed study

The fast subset produced 10 scenario-method rows. The full factorial produced
1,008 rows: 144 scenarios across four ComBat-refQL configurations, raw-count
edgeR, pinned ComBat-ref, and ComBat-seq. All scenarios completed. No tested
reference batch changed; reference invariance is not applicable to ComBat-seq,
which does not retain a reference batch.

| Method | Mean delta RMSE | Mean d RMSE | Mean biology LFC RMSE | Mean absolute change |
|---|---:|---:|---:|---:|
| none + shared | 0.383 | 0.462 | 0.573 | 16.55 |
| ql_eb + shared | 0.200 | 0.462 | 0.569 | 6.84 |
| none + batch_scaled | 0.382 | 0.200 | 0.567 | 17.32 |
| ql_eb + batch_scaled | 0.200 | 0.200 | 0.566 | 8.52 |
| ComBat-ref | — | — | 0.517 | 12.25 |
| ComBat-seq | — | — | 0.580 | 12.85 |
| raw edgeR design input | — | — | 0.575 | 0.00 |

The batch-scaled pathways had mean absolute log-dispersion error 0.332, 0.161,
and 0.107 when the true source/reference ratios were 0.5, 1, and 2. Shared
dispersion fixes the shift at zero and therefore has absolute errors 0.693, 0,
and 0.693 in those settings. QL-EB reduced contrast RMSE in this simulation and
reduced transformation magnitude. These are simulation results, not general
superiority claims.

Raw-count edgeR design modelling was the only workflow used for inferential FDR
and power. Its scenario values are retained in the tidy result files; adjusted
counts were not presented as an inference replacement.

## Unit and numerical validation

Deterministic unit tests cover EB uncertainty limits and complete pooling,
factor-order invariance, QL covariance agreement with an edgeR one-df test,
dispersion identities, residual-information totals and saturation, increasing
precision, exact zero mid-P probability, equal-distribution identity,
dispersion-sensitive transport, all four ablations, exact reference invariance,
chunk equivalence, rollback, and fixed-input determinism.

`devtools::document()` completed without warnings. `devtools::test()` passed
114 assertions with zero failures, warnings, or skips. `devtools::check(
document = FALSE, manual = FALSE, error_on = "never")` completed with zero
errors, warnings, or notes. `git diff --check` passed.

## Limitations

The factorial study is synthetic and uses two batches. The hierarchical model
has one batch-wide multiplier and cannot represent gene-by-batch dispersion
interactions. ComBat-ref and ComBat-seq do not expose contrast and dispersion
objects directly comparable to the new diagnostics, so their delta and shift
RMSE entries remain unavailable. Automatic-reference behavior warrants further
study with independently characterized technical references.

## Reproducibility

```sh
Rscript validation/simulation.R fast
Rscript validation/simulation.R full
```

Outputs are `validation/results-fast.tsv` and
`validation/results-full.tsv`. The upstream comparators are pinned to the
commits listed in the validation plan.

```text
R version 4.6.1 (2026-06-24)
Platform: aarch64-apple-darwin25.4.0
Running under: macOS Tahoe 26.6
BLAS: OpenBLAS 0.3.34
LAPACK: 3.12.1
Locale: C.UTF-8
Time zone: Europe/Rome
edgeR: 4.10.1
S7: 0.2.2
```
