# ComBat-refQL

[![R CMD check](https://github.com/Thokas99/ComBat-refQL/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Thokas99/ComBat-refQL/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/Thokas99/ComBat-refQL/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/Thokas99/ComBat-refQL/actions/workflows/pkgdown.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

**Reference-batch correction for raw RNA-seq counts.**

ComBat-refQL combines TMMwsp normalization, edgeR v4 quasi-likelihood fitting,
counterfactual reference-batch prediction, and deterministic negative-binomial
mid-P transport. It preserves count-matrix dimensions and reports unsupported
genes and lossy input actions.

> [!IMPORTANT]
> QL fitting supplies counterfactual means; it does not define a count CDF.
> Negative-binomial and QL dispersions remain strictly separate.

## Table of contents

- [Why ComBat-refQL?](#why-combat-refql)
- [Installation](#installation)
- [Quick start](#quick-start)
- [What the result contains](#what-the-result-contains)
- [Input safety](#input-safety)
- [Statistical approach](#statistical-approach)
- [Limitations](#limitations)
- [Documentation](#documentation)
- [Citation](#citation)
- [License](#license)

## Why ComBat-refQL?

- **Count-aware:** models raw RNA-seq counts with negative-binomial means and
  dispersions.
- **Reference-based:** maps non-reference batches toward a selected reference
  while preserving biological-group and covariate terms.
- **Auditable:** reports designs, offsets, normalization factors, reference
  scores, unsupported genes, fallbacks, and input conversions.
- **Safe by default:** audits fractional-count rounding and rejects invalid or
  normalized expression values.

## Installation

Install the development version directly from GitHub with
[`pak`](https://pak.r-lib.org/):

```r
install.packages("pak")
pak::pak("Thokas99/ComBat-refQL")
```

Then load the package:

```r
library(combatrefql)
```

ComBat-refQL requires R 4.3 or newer and edgeR 4.2 or newer.

## Quick start

`counts` must contain genes in rows and samples in columns. The batch and group
vectors must describe those same samples.

```r
fit <- combat_ref_ql(
  counts = counts,
  batch = metadata$batch,
  group = metadata$condition,
  reference = "batch_A"
)

adjusted_counts <- fit@counts
```

If `reference` is omitted, ComBat-refQL selects one automatically and returns
the batch scores used in that decision. An explicit
reference is recommended when study design or external QC identifies the best
batch.

The defaults are:

| Setting | Default |
|---|---|
| Normalization | `TMMwsp` |
| Fit | robust edgeR v4 QL |
| Transport | deterministic NB mid-P |
| Reference | automatic unless specified |
| Messages | normal, quiet, or verbose |

## What the result contains

```r
fit@counts                          # non-negative, integer-valued doubles
fit@specification                   # selected reference and method
fit@gene_status                     # adjusted, unsupported, or failed genes
fit@diagnostics$input_actions       # audited lossy input actions
summary(fit)
```

## Input safety

Fractional values are rounded by default, announced in normal and verbose
output, and stored in `fit@diagnostics$input_actions`. Strict rejection is
available with:

```r
combat_ref_ql(..., fractional_counts = "error")
```

Negative, missing, infinite, unparseable, ambiguous, or normalized inputs remain
errors. TPM, CPM, logCPM, RPKM, and FPKM are not raw counts and are rejected.

## Statistical approach

For gene \(g\) and sample \(i\), ComBat-refQL fits

\[
\log(\mu_{gi}) = o_i + X_i\beta_g,
\]

where the offset uses the original library size and its normalization factor.
Counterfactual prediction changes only the batch indicators to those of the
reference batch, preserving group, covariates, and sample-specific offsets.

Automatic-reference scores use fresh intercept-bearing designs within each
batch. Missing factor levels and non-estimable terms are handled explicitly and
recorded. edgeR v4 QL fitting estimates a global underlying NB dispersion and
stores it separately from prior and posterior QL dispersions. Only the NB
dispersion enters the mid-P probability and quantile functions.

See the [statistical methodology vignette](vignettes/statistical-methodology.Rmd)
for equations, assumptions, and implementation details.

## Limitations

- The QL-first pathway has not passed broad simulation or external
  real-data validation and must not be considered production-ready.
- Automatic-reference stability is not yet bootstrapped; selection frequency
  is therefore not reported.
- Genes with an all-zero batch lack enough information for correction and are
  returned unchanged with an explicit status.
- Corrected counts do not propagate correction uncertainty into downstream
  differential-expression tests.
- Complete batch/group confounding cannot be repaired and is rejected.

## Documentation

- [Changelog](CHANGELOG.md)
- [Release notes](NEWS.md)
- [Implementation report](IMPLEMENTATION_REPORT.md)
- [Validation report](validation/VALIDATION_REPORT.md)

## Citation

The original ComBat-ref method is described in:

> Zhang X. Highly effective batch effect correction method for RNA-seq count
> data. *Computational and Structural Biotechnology Journal*. 2025.
> doi:[10.1016/j.csbj.2024.12.010](https://doi.org/10.1016/j.csbj.2024.12.010).
> PMID: [39802213](https://pubmed.ncbi.nlm.nih.gov/39802213/).

In R, use `citation("combatrefql")` for the package citation record.

## License

The newly written ComBat-refQL package implementation is available under the
[MIT License](LICENSE). This license does not invent or alter licensing terms
for upstream code, publications, or data.
