# ComBat-refQL 0.1.0 implementation report

## Executive summary

The checkout started at `bbce315` on branch `ci/fix-pages-environment` with
version `0.0.0.9000`; 88 tests passed. Version 0.1.0 separates frozen legacy
behavior from a repaired NB baseline, defaults to TMMwsp, uses counterfactual
reference prediction and stable mid-P mapping, exposes unsupported genes, and
adds diagnostic-only edgeR QL support. No validation or superiority claim is made.

## Repository state and inventory

The pre-existing untracked `.serena/` directory was preserved. Material files:

| File | Change |
|---|---|
| `R/combat_ref.R` | New public baseline, legacy wrapper, and QL alias |
| `R/combat_ref_baseline.R` | Normalization, dispersion, fitting, prediction, mapping |
| `R/build_design.R` | Canonical named design and counterfactual design |
| `R/quantile_mapping.R` | Mid-P/randomized methods and safe tails |
| `R/validate_inputs.R` | Covariates, one-level group, and exact alignment |
| `R/result_class.R` | Method metadata and gene status |
| `inst/scripts/combat-ref` | Resolved methods, seed, and diagnostics TSV |
| `tests/testthat/` | Legacy, design, offset, mapping, QL, status, and CLI checks |
| Documentation | README, NEWS, changelog, migration, vignettes, validation report |

## Mathematical change

Historically, the null-group code joined two intercept-bearing designs and
deleted every all-one column. The reference batch then had an all-zero design
row and no gene baseline. The repaired design begins with `~ 0 + batch`, adds
only nonredundant biological terms, and requires `qr(design)$rank == ncol(design)`.

For each gene/sample, the fitted mean is
`exp(offset + X beta)`. The reference mean uses the same offset, group, and
covariates but replaces the batch indicators in `X` with the requested reference.
This removes coefficient-position assumptions and the discontinuous historical
`gamma < -0.2` cap.

## Normalization and offsets

The initial DGEList is created before algorithmic gene filtering. Raw library
sizes are set once, TMMwsp factors are calculated once with edgeR
`normLibSizes()`, and subsetting retains those sizes/factors. TMMwsp is the
default because its singleton pairing is intended for zero-heavy asymmetric
count matrices; TMM and none remain explicit alternatives.

## Quantile mapping

The default uses `F_old(y-1) + 0.5 P_old(Y=y)` and an NB inverse CDF under the
counterfactual mean/reference dispersion. Probabilities are clamped to finite
machine bounds. Randomized mapping is opt-in and requires a seed. QL dispersion
is never used as NB `size`.

## QL/NB distinction

`glmQLFit(..., robust=TRUE)` supplies moderated `s2.post` diagnostics. The model
stores these separately from edgeR tagwise NB dispersion. Coefficient shrinkage
was not implemented because no derived/tested coefficient-variance estimator and
prior-variance validation gate were established.

## Reference and unsupported genes

Automatic selection still uses comparable robust per-batch common-dispersion
point scores. Frequencies are explicitly unresolved (`reference_stable = NA`);
an explicit reference is recommended. Genes all zero within any batch remain
unchanged to preserve dimensions and are labeled `unsupported_zero_batch`.

## Tests and numerical evidence

- Baseline before edits: 88 pass, 0 fail/warn/skip.
- Implemented suite: 106 pass, 0 fail/warn/skip.
- Legacy frozen fixtures remain exact through `combat_ref_legacy()`.
- Smoke simulation: raw B/A ratio 1.7112; adjusted 1.6447. This is not a validation gate.

## API and migration

The default is `combat_ref(..., normalization="TMMwsp", fit_method="baseline",
quantile_method="midp")`. Historical output requires `combat_ref_legacy()`.
See `MIGRATION.md` for user-visible numerical changes.

## Risks and unresolved work

- Licensing and upstream redistribution terms remain unresolved.
- Automatic-reference bootstrap stability is not implemented.
- The broad preregistered simulation matrix, real-data validation, and peak-memory benchmark remain for 0.2.0.
- Corrected counts do not propagate correction uncertainty into DE inference;
  direct edgeR QL analysis of raw counts with batch and group remains preferred.

## Commands and environment

Executed:

```r
devtools::document()
devtools::test()
devtools::check(document = FALSE, manual = FALSE, error_on = "never")
pkgdown::build_site()
source("validation/simulation.R")
sessionInfo()
```

Final results: 106 tests passed with no failures, warnings, or skips; R CMD
check completed with 0 errors, 0 warnings, and 0 notes; the methodology vignette
and pkgdown site built successfully. R 4.6.1, edgeR 4.10.1, macOS Tahoe 26.5.2,
arm64. Starting commit: `bbce315`; implementation commit: `3a66f1c`.
