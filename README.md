# ComBat-refQL

**Reference-batch adjustment of RNA-seq count matrices with edgeR
quasi-likelihood modelling and deterministic count transport.**

[![R CMD check](https://github.com/Thokas99/ComBat-refQL/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Thokas99/ComBat-refQL/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/Thokas99/ComBat-refQL/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/Thokas99/ComBat-refQL/actions/workflows/pkgdown.yaml)
[![Version](https://img.shields.io/badge/version-0.4.0-315B6D.svg)](DESCRIPTION)
[![License: MIT](https://img.shields.io/badge/license-MIT-315B6D.svg)](LICENSE)

ComBat-refQL fits RNA-seq counts with edgeR quasi-likelihood methods, moderates
gene-specific batch contrasts using their estimation uncertainty, estimates
batch-specific negative-binomial dispersion relative to an observed reference,
and transports non-reference counts toward that reference. Observed reference
counts remain unchanged.

ComBat-refQL is independently implemented from the statistical principles
described in ComBat-seq and ComBat-ref.

## Why ComBat-refQL?

- One fitting function: `combat_ref_ql()`.
- Explicit or auditable automatic reference selection.
- Biological-group and sample-covariate modelling.
- TMMwsp normalization and edgeR v4 model fitting.
- QL-aware empirical-Bayes moderation of batch contrasts.
- Hierarchical batch-specific NB dispersion estimation.
- Deterministic, chunked NB mid-P count transport.
- Typed conditions, validated S7 objects, and tidy diagnostics.

## Installation

```r
install.packages("pak")
pak::pak("Thokas99/ComBat-refQL")
```

ComBat-refQL requires R 4.3 or later, edgeR 4.2 or later, and S7 0.2.0 or
later.

## Quick start

```r
library(combatrefql)

fit <- combat_ref_ql(
  counts = counts,
  batch = metadata$batch,
  group = metadata$condition,
  reference = "reference_batch"
)

adjusted_counts <- fit@counts
```

Use `reference = NULL` for automatic selection. The diagnostics record every
candidate score, the selected and second-best candidates, and their score
margin. Plausible fractional estimated counts are rounded by default and the
action is recorded; use `fractional_counts = "error"` to require integers at
entry.

Sample-aligned numeric or factor covariates enter the same model:

```r
fit <- combat_ref_ql(
  counts = counts,
  batch = metadata$batch,
  group = metadata$condition,
  covariates = metadata[colnames(counts), c("sex", "age"), drop = FALSE],
  reference = "reference_batch"
)
```

## How it works

ComBat-refQL builds a full-rank biological and batch design and fits a
preliminary edgeR v4 QL model. It uses GLM leverage and working weights to
estimate continuous residual information, then estimates a moderated
gene-dispersion baseline and reference-centred batch multipliers.

The final edgeR fit provides each source-versus-reference contrast and its
QL-scaled variance. A normal empirical-Bayes model moderates these contrasts,
which define reference-counterfactual means. Deterministic mid-P transport maps
each supported source count from its fitted source NB distribution to the
reference NB distribution. The same rule applies to zero and positive counts.

## Output and diagnostics

```r
fit
fit_summary <- summary(fit)

fit_summary@batches
fit_summary@shrinkage
fit_summary@dispersion
fit_summary@gene_outcomes
fit@diagnostics$mapping
fit@gene_status
```

Printing stays compact:

```text
ComBat-refQL fit

Data
  18,432 genes x 48 samples
  3 batches - 2 biological groups

Model
  Model fitting:   edgeR v4
  Batch contrasts: QL-moderated
  Dispersion:      hierarchical NB
  Normalization:   TMMwsp
  Transport:       negative-binomial mid-P
  Reference:       batch_A - explicit

Result
  17,904 adjusted
  528 unchanged
  0 failed

Runtime: 18.7 seconds
```

Adjusted counts support batch-aligned visualization, exploratory analysis, and
workflows that require a count-scale matrix. Differential-expression inference
uses the original counts with batch and biological variables represented in the
statistical design.

## Validation

Deterministic tests cover QL shrinkage limits, contrast covariance, residual
information, batch-specific dispersion identities, mid-P mapping, numerical
rollback, chunk invariance, and exact reference preservation. A factorial
simulation uses known mean and dispersion effects. Reproducible GFRN and NASA
GeneLab analyses evaluate empirical batch and biological structure against
uncorrected counts, ComBat-ref, and ComBat-seq.

See the [validation overview](validation/README.md) and
[validation report](validation/VALIDATION_REPORT.md). Real-data scripts fetch
pinned public inputs into an ignored local cache.

## Documentation

- [Statistical methodology](vignettes/statistical-methodology.Rmd)
- [Validation article](vignettes/validation.Rmd)
- [Function and object reference](https://thokas99.github.io/ComBat-refQL/reference/)
- [News](NEWS.md)

## Citation

- Zhang X. Highly effective batch effect correction method for RNA-seq count
  data. *Computational and Structural Biotechnology Journal*. 2025.
  [doi:10.1016/j.csbj.2024.12.010](https://doi.org/10.1016/j.csbj.2024.12.010).
- Zhang Y, Parmigiani G, Johnson WE. ComBat-seq: batch effect adjustment for
  RNA-seq count data. *NAR Genomics and Bioinformatics*. 2020;2(3):lqaa078.
  [doi:10.1093/nargab/lqaa078](https://doi.org/10.1093/nargab/lqaa078).
- Run `citation("edgeR")` for references appropriate to the installed edgeR
  version and `citation("combatrefql")` for package citation metadata.

## Acknowledgements

ComBat-refQL extends the reference-batch framework introduced in ComBat-ref.
The original ComBat-ref implementation is available under the MIT License and
is attributed to Xiaoyu Zhang. The package license is in [LICENSE](LICENSE);
the preserved upstream notice is in
[LICENSES/ComBat-ref-MIT.txt](LICENSES/ComBat-ref-MIT.txt).
