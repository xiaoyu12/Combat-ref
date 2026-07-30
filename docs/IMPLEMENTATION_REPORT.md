# ComBat-refQL 0.1.0 implementation report

## Executive summary

Version 0.1.0 provides one supported negative-binomial correction pathway. It
defaults to TMMwsp, uses a full-rank named design and counterfactual reference
prediction, exposes unsupported genes, and offers diagnostic-only edgeR QL
support. No validation or superiority claim is made.

Work started from `719f2c8` on `main`. No `0.1.0` Git tag or GitHub release was
present, so these corrections retain version 0.1.0 rather than claiming 0.1.1.
The ending commit is recorded in the final handoff because this report is part
of that commit.

## File inventory

| File | Responsibility |
|---|---|
| `R/combat_ref.R` | Public baseline and experimental QL alias |
| `R/combat_ref_baseline.R` | Normalization, dispersion, fitting, prediction, mapping |
| `R/build_design.R` | Canonical named design and counterfactual design |
| `R/quantile_mapping.R` | Mid-P, upper-CDF, randomized mapping and safe tails |
| `R/validate_inputs.R` | Input, alignment, and identifiability checks |
| `R/prepare_inputs.R` | Safe coercion, explicit lossy policies, and audit log |
| `R/result_class.R` | Result metadata and gene status |
| `inst/scripts/combat-ref` | Structured CLI, file parsing, diagnostics, and exit codes |
| `tests/testthat/` | Design, offset, mapping, QL, status, API, and CLI checks |

## Statistical model

The design begins with `~ 0 + batch`, adds nonredundant biological terms, and
requires `qr(design)$rank == ncol(design)`. For each gene and sample, the fitted
mean is `exp(offset + X beta)`. The counterfactual reference mean keeps the same
offset, group, and covariates while replacing only the batch indicators.

Reference scoring no longer derives a local model by deleting global batch
columns. Each batch gets a fresh `~ 1 + group + covariates` design after local
factor levels are dropped. Non-estimable trailing terms are removed with an
audited warning; the final fallback is `~ 1`. Diagnostics retain the design,
dropped columns, missing levels, rank, residual degrees of freedom, and reason.

The initial DGEList is created before gene filtering. TMMwsp factors are
calculated once with edgeR `normLibSizes()`, and subsetting retains the same
library sizes and factors. TMM and no normalization remain explicit options.

The default discrete mapping uses
`F_old(y-1) + 0.5 P_old(Y=y)` and an NB inverse CDF under the reference mean and
dispersion. Probabilities are clamped to finite machine bounds. Randomized
mapping is opt-in and requires a seed.

## QL/NB distinction

Under the pinned edgeR >= 4.2.0 API, `glmQLFit(..., dispersion = NULL,
robust=TRUE, legacy=FALSE)` estimates a global NB dispersion for QL fitting.
Baseline tagwise NB dispersion, this QL global NB dispersion, and moderated
`s2.post` are separate diagnostic fields. QL dispersion is never used as an NB
distribution parameter. Coefficient shrinkage is not implemented without a
derived and validated uncertainty model.

## Reference and unsupported genes

Automatic selection uses robust per-batch common-dispersion point scores.
The selected batch and method are reported, but selection frequency is `NULL`
and stability is `NA` until resampling exists. Genes all zero within a batch
remain unchanged to preserve dimensions and are labeled
`unsupported_zero_batch`. Fit convergence is `NA` unless edgeR supplies a
gene-level failure indicator; detected numerical failures are `FALSE`.

## Input policy and CLI

The preparation layer losslessly handles data frames, numeric character
columns, a unique first gene-ID column, integer storage, unambiguous transpose,
factor conversion, one-level groups, method spelling, and exact metadata
reordering. Every action is printed in normal CLI mode and stored in
`diagnostics$input_actions`. Fractional rounding and duplicate-gene summation
require explicit options. Negative, non-finite, normalized-scale, unparseable,
or ambiguous inputs remain errors.

The CLI prints Input, Configuration, Design, Normalization, Reference,
Correction, Output, warning recap, and timing sections by default. `--verbose`
adds detail, `--quiet` is pipeline-safe, and `--validate-only` stops before the
fit. JSON diagnostics and gene-status TSV are optional. Exit codes are 2 for
input, 3 for design, 4 for fitting/mapping, and 5 for output failures.

```text
ComBat-refQL 0.1.0
Input          200 genes x 8 samples
Configuration  baseline; TMMwsp; mid-P; automatic reference
Design         full rank; residual degrees of freedom reported
Reference      selected batch and auditable local score table
Correction     supported, unchanged, and numerical-tail counts
Output         adjusted counts and requested diagnostics paths
```

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

The final unit run has 146 passing tests with no failures, warnings, or skips.
R CMD check completed with zero errors, warnings, or notes, and pkgdown built
successfully. Eight manual CLI scenarios passed: clean input, data-frame and
gene-ID handling, metadata reordering, explicit and automatic references,
QL-assisted fitting, rejected and explicitly rounded fractional counts, and an
unambiguous transpose. The standalone simulation command requires an installed
package; from the checkout it was run after `devtools::load_all(quiet=TRUE)`.
The smoke simulation moved the raw B/A ratio from 1.7112 to 1.6579; this is not
a validation gate. Environment: R 4.6.1, edgeR 4.10.1, macOS Tahoe 26.5.2,
arm64.

## Remaining work

- Bootstrap automatic-reference stability.
- Run the preregistered multi-scenario simulation and real-data study.
- Benchmark peak memory and runtime.
- Add broad real-file testing beyond the deterministic CLI fixtures.
- Validate any future QL-informed coefficient shrinkage before enabling it.
- Continue to recommend direct edgeR QL analysis of raw counts for DE inference.
