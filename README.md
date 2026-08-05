# ComBat-refQL

**ComBat-refQL provides a clean R interface for reference-batch adjustment of
RNA-seq count matrices.**

It combines current edgeR model-fitting infrastructure with deterministic
negative-binomial count transport, explicit reference selection, validated S7
objects, and structured diagnostics.

[![R CMD check](https://github.com/Thokas99/ComBat-refQL/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Thokas99/ComBat-refQL/actions/workflows/R-CMD-check.yaml)
[![Version](https://img.shields.io/badge/version-0.3.2-4C78A8.svg)](DESCRIPTION)
[![R](https://img.shields.io/badge/R-%E2%89%A5%204.3-276DC3.svg)](https://www.r-project.org/)
[![edgeR](https://img.shields.io/badge/edgeR-%E2%89%A5%204.2-87B13F.svg)](https://bioconductor.org/packages/edgeR/)
[![License: MIT](https://img.shields.io/badge/license-MIT-2E8B57.svg)](LICENSE)

ComBat-refQL is an independent package inspired by the reference-batch strategy
described in ComBat-ref and the count-modelling foundations of ComBat-seq. It
centres that workflow in one focused function with auditable inputs, fitting,
transport, and results.

## Core capabilities

- one fitting function, `combat_ref_ql()`;
- automatic or explicit reference selection;
- biological-group and sample-covariate modelling;
- validated alignment of counts and sample metadata;
- audited handling of fractional estimated counts;
- TMMwsp normalization and edgeR v4 model fitting;
- model-based reference-counterfactual means;
- deterministic, chunked negative-binomial mid-P transport;
- typed conditions and informative progress;
- validated S7 fit and summary objects;
- concise printing and tidy diagnostic components.

## Installation

Install the development release from GitHub:

```r
install.packages("pak")
pak::pak("Thokas99/ComBat-refQL")
```

ComBat-refQL requires R 4.3 or later, edgeR 4.2 or later, and S7 0.2.0 or
later.

## Quick start

Prepare a numeric count matrix with genes in rows and samples in columns. Its
column order must match the sample metadata or its identifiers must support
unambiguous alignment.

```r
library(combatrefql)

fit_ql <- combat_ref_ql(
  counts = counts,
  batch = metadata$batch,
  group = metadata$condition
)

fit_ql
summary(fit_ql)
adjusted_counts <- fit_ql@counts
```

With `reference = NULL`, the default, ComBat-refQL selects the reference from
biologically adjusted within-batch models. Fixed inputs produce a deterministic
selection and transport result.

Plausible fractional estimated counts are rounded by default. The action and
number of affected values are recorded in
`fit_ql@diagnostics$input_actions`. Set `fractional_counts = "error"` when the
input workflow requires integer-valued counts at entry.

## Explicit reference

Supply a reference batch when study design or external quality-control evidence
identifies the appropriate target:

```r
fit_ql <- combat_ref_ql(
  counts = counts,
  batch = metadata$batch,
  group = metadata$condition,
  reference = "batch_A"
)
```

## Sample covariates

Sample-aligned numeric or factor covariates enter the model with the biological
group. Here the metadata row names are the sample identifiers:

```r
fit_ql <- combat_ref_ql(
  counts = counts,
  batch = metadata$batch,
  group = metadata$condition,
  covariates = metadata[colnames(counts), c("sex", "age"), drop = FALSE]
)
```

## Printed fit

Printing a fit reports its data dimensions, model specification, selected
reference, gene outcomes, and elapsed time:

```text
── ComBat-refQL fit ─────────────────────────────────────────────

Data
  18,432 genes x 48 samples
  3 batches - 2 biological groups

Model
  Model fitting:   edgeR v4
  Normalization:  TMMwsp
  Transport:      negative-binomial mid-P
  Reference:      batch_A - explicit

Result
  ✔ 17,904 adjusted
  ℹ 528 unchanged
  ✖ 0 failed

Runtime: 18.7 seconds
Use summary(x) for diagnostics.
```

## Summary and adjusted counts

The fit and summary use S7 property access. The adjusted matrix retains the
validated gene and sample identifiers:

```r
fit_summary <- summary(fit_ql)
adjusted_counts <- fit_ql@counts

fit_summary@batches
fit_summary@gene_outcomes
fit_summary@input_actions
fit_ql@gene_status
fit_ql@diagnostics
```

`fit_summary@batches` records reference scores, the selected batch, and local
design details. `fit_ql@gene_status` records the outcome for every gene.

## Implementation overview

ComBat-refQL validates and aligns inputs, applies TMMwsp normalization, builds a
full-rank design, and fits the count matrix with the edgeR v4 quasi-likelihood
pipeline. For each non-reference sample, the design substitutes the reference
batch indicators while retaining its biological group and covariates. This
produces model-based reference-counterfactual means.

Observed counts are transported to those means with deterministic
negative-binomial mid-P probabilities. Transport runs by non-reference batch
and gene chunk; `chunk_size` controls the temporary transport allocation. The
implementation records input actions, reference scoring, mapping outcomes,
model diagnostics, warnings, and timings in tidy components.

The [statistical methodology vignette](vignettes/statistical-methodology.Rmd)
documents the model, reference scoring, and transport equations.

## Recommended analysis use

Adjusted counts support batch-aligned visualization, exploratory analysis, and
workflows that require a count-scale matrix. Differential-expression inference
uses the original counts with batch and biological variables represented in the
statistical design.

An explicit reference supported by study design or external quality control
provides the clearest interpretation. Automatic selection supplies a
deterministic reference score and records the locally estimable biological
terms used for each candidate.

## Citations and attribution

ComBat-refQL is independently implemented from the statistical principles
described in ComBat-seq and ComBat-ref.

- Zhang X. Highly effective batch effect correction method for RNA-seq count
  data. *Computational and Structural Biotechnology Journal*. 2025.
  [doi:10.1016/j.csbj.2024.12.010](https://doi.org/10.1016/j.csbj.2024.12.010).
- Zhang Y, Parmigiani G, Johnson WE. ComBat-seq: batch effect adjustment for
  RNA-seq count data. *NAR Genomics and Bioinformatics*. 2020;2(3):lqaa078.
  [doi:10.1093/nargab/lqaa078](https://doi.org/10.1093/nargab/lqaa078).
- Run `citation("edgeR")` for references appropriate to the installed edgeR
  version.

Run `citation("combatrefql")` for package citation metadata.

ComBat-refQL extends the reference-batch framework introduced in ComBat-ref.
The original ComBat-ref implementation is available under the MIT License and
is attributed to Xiaoyu Zhang. The package license is in [LICENSE](LICENSE),
and the preserved upstream notice is in
[LICENSES/ComBat-ref-MIT.txt](LICENSES/ComBat-ref-MIT.txt).
