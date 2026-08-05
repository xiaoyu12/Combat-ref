# ComBat-refQL 0.4.0 implementation report

## Algorithms and edgeR quantities

The implementation runs `normLibSizes(method = "TMMwsp")`, constructs a
full-rank design, and calls `glmQLFit(..., robust = TRUE, legacy = FALSE,
prior.count = 0, keep.unit.mat = TRUE)` for preliminary fitting. Installed
edgeR 4.10.1 accepts a dispersion matrix with the same dimensions as the count
matrix; 0.4.0 uses that supported path for the final fit.

`estimateDisp(..., robust = TRUE, trend.method = "locfit")` supplies the
moderated gene-level NB baseline. A bounded Brent optimization estimates each
batch log multiplier from the fixed-mean NB likelihood weighted by inverse
trigamma effective information. Bounds constrain multipliers to `[0.1, 10]`;
NB dispersions are kept finite and positive. The reference-centred shift is
exactly zero.

The coefficient source is `unshrunk.coefficients` when present and
`coefficients` otherwise. The contrast covariance uses the final design and
working weights

$$w_{gj}=a_{gj}\mu_{gj}/\{1+(\phi_{gj}/\bar q)\mu_{gj}\},$$

then multiplies $c^T(X^TW_gX)^{-1}c$ by `s2.post[g]`. A deterministic test
compares the resulting one-df Wald statistic with `glmQLFTest`; their simulated
correlation exceeds 0.999 and the 90th-percentile relative difference is below
0.03. No p-value is inverted to obtain an error.

The empirical-Bayes prior is fitted by profile marginal likelihood with an
explicit zero-variance boundary. Transport calls base R NB functions with
distinct source and target dispersions. It clips extreme probabilities,
validates integer-valued finite output, and rolls failed genes back atomically.

## Convergence and safeguards

Each batch multiplier uses one deterministic one-dimensional optimization and
records the optimizer convergence code and function-evaluation count. A source
with zero precision receives a neutral shift and `not_estimable`. A zero-
information reference raises a typed statistical error. Singular gene-level
working information contributes zero rather than stopping unrelated genes.

## File map

- `R/combat-ref-ql.R`: ordered orchestration and public controls.
- `R/fit-ql.R`: edgeR fits and QL-scaled contrast covariance.
- `R/build-design.R`: exact contrasts and effective information.
- `R/estimate-batch-dispersion.R`: hierarchical NB dispersion.
- `R/shrink-batch-effects.R`: normal-normal empirical Bayes.
- `R/transport-counts.R`: source-to-reference mid-P transport.
- `R/build-result.R` and S7 files: validated model and tidy diagnostics.

## Verification environment and commands

Development used R 4.6.1 on arm64 macOS, edgeR 4.10.1, and S7 0.2.2. The
recorded commands are:

```r
devtools::document()
devtools::test()
devtools::check(document = FALSE, manual = FALSE, error_on = "never")
```

```sh
Rscript validation/simulation.R fast
Rscript validation/simulation.R full
git diff --check
```

The final check results are recorded after execution in the validation report.

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
