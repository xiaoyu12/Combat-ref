# ComBat-refQL 0.1.0 implementation report

## Executive summary

Version 0.1.0 provides one supported negative-binomial correction pathway. It
defaults to TMMwsp, uses a full-rank named design and counterfactual reference
prediction, exposes unsupported genes, and offers diagnostic-only edgeR QL
support. No validation or superiority claim is made.

## File inventory

| File | Responsibility |
|---|---|
| `R/combat_ref.R` | Public baseline and experimental QL alias |
| `R/combat_ref_baseline.R` | Normalization, dispersion, fitting, prediction, mapping |
| `R/build_design.R` | Canonical named design and counterfactual design |
| `R/quantile_mapping.R` | Mid-P, upper-CDF, randomized mapping and safe tails |
| `R/validate_inputs.R` | Input, alignment, and identifiability checks |
| `R/result_class.R` | Result metadata and gene status |
| `inst/scripts/combat-ref` | CLI and diagnostics TSV |
| `tests/testthat/` | Design, offset, mapping, QL, status, API, and CLI checks |

## Statistical model

The design begins with `~ 0 + batch`, adds nonredundant biological terms, and
requires `qr(design)$rank == ncol(design)`. For each gene and sample, the fitted
mean is `exp(offset + X beta)`. The counterfactual reference mean keeps the same
offset, group, and covariates while replacing only the batch indicators.

The initial DGEList is created before gene filtering. TMMwsp factors are
calculated once with edgeR `normLibSizes()`, and subsetting retains the same
library sizes and factors. TMM and no normalization remain explicit options.

The default discrete mapping uses
`F_old(y-1) + 0.5 P_old(Y=y)` and an NB inverse CDF under the reference mean and
dispersion. Probabilities are clamped to finite machine bounds. Randomized
mapping is opt-in and requires a seed.

## QL/NB distinction

`glmQLFit(..., robust=TRUE)` supplies moderated `s2.post` diagnostics. These are
stored separately from edgeR tagwise NB dispersion. QL dispersion is never used
as an NB distribution parameter. Coefficient shrinkage is not implemented
without a derived and validated uncertainty model.

## Reference and unsupported genes

Automatic selection uses robust per-batch common-dispersion point scores.
Stability is explicitly unresolved (`reference_stable = NA`), so an explicit
reference remains recommended. Genes all zero within a batch remain unchanged
to preserve dimensions and are labeled `unsupported_zero_batch`.

## Verification

Commands run:

```r
devtools::document()
devtools::test()
devtools::check(document = FALSE, manual = FALSE, error_on = "never")
pkgdown::build_site()
source("validation/simulation.R")
sessionInfo()
```

The final suite has 78 passing tests with no failures, warnings, or skips. R CMD
check completed with zero errors, warnings, or notes, and pkgdown built cleanly.
The smoke simulation moved the raw B/A ratio from 1.7112 to 1.6447; this is not
a validation gate. Environment: R 4.6.1, edgeR 4.10.1, macOS Tahoe 26.5.2,
arm64.

## Remaining work

- Bootstrap automatic-reference stability.
- Run the preregistered multi-scenario simulation and real-data study.
- Benchmark peak memory and runtime.
- Validate any future QL-informed coefficient shrinkage before enabling it.
- Continue to recommend direct edgeR QL analysis of raw counts for DE inference.
