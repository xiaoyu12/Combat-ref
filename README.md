# ComBat-refQL

**A clean R interface for reference-batch adjustment of RNA-seq count data.**

[![R CMD check](https://github.com/Thokas99/ComBat-refQL/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Thokas99/ComBat-refQL/actions/workflows/R-CMD-check.yaml)
[![R](https://img.shields.io/badge/R-%E2%89%A5%204.3-276DC3.svg)](https://www.r-project.org/)
[![edgeR](https://img.shields.io/badge/edgeR-%E2%89%A5%204.2-87B13F.svg)](https://bioconductor.org/packages/edgeR/)

ComBat-refQL provides a compact workflow for adjusting non-reference RNA-seq batches toward a selected reference batch. It combines current edgeR model-fitting infrastructure with deterministic negative-binomial count transport and returns a structured S7 result.

The package focuses on a practical analysis experience:

- one main fitting function;
- automatic or explicit reference selection;
- support for biological groups and sample covariates;
- validated inputs and informative errors;
- concise progress reporting;
- deterministic, chunked count transport;
- structured `print()` and `summary()` output.

ComBat-refQL is an independent package inspired by the reference-batch strategy described in ComBat-ref and the count-modelling foundations of ComBat-seq.

## Installation

Install the development version from GitHub:

```r
install.packages("pak")
pak::pak("Thokas99/ComBat-refQL")
```

Load the package:

```r
library(combatrefql)
```

ComBat-refQL requires R 4.3 or later, edgeR 4.2 or later, and S7 0.2.0 or later.

## Quick start

Provide a raw count matrix with genes in rows and samples in columns. Sample metadata must follow the same order as the matrix columns or be named for automatic alignment.

```r
fit_ql <- combat_ref_ql(
  counts = counts,
  batch = metadata$batch,
  group = metadata$condition
)
```

When `reference` is omitted, the package evaluates the available batches and selects the reference with the smallest valid dispersion score.

Select a reference explicitly when it is established by the study design or external quality-control information:

```r
fit_ql <- combat_ref_ql(
  counts = counts,
  batch = metadata$batch,
  group = metadata$condition,
  reference = "batch_A"
)
```

Additional sample covariates can be included in the model:

```r
fit_ql <- combat_ref_ql(
  counts = counts,
  batch = metadata$batch,
  group = metadata$condition,
  covariates = metadata[c("sex", "age")],
  reference = "batch_A"
)
```

Fractional estimated counts are rounded by default and recorded in the fit diagnostics:

```r
fit_ql <- combat_ref_ql(
  counts = counts,
  batch = metadata$batch,
  group = metadata$condition,
  fractional_counts = "round"
)
```

Use quiet output inside automated workflows:

```r
fit_ql <- combat_ref_ql(
  counts = counts,
  batch = metadata$batch,
  group = metadata$condition,
  verbosity = "quiet"
)
```

## Results

Printing the fitted object gives a concise overview of the data, model, selected reference, gene outcomes, and runtime:

```r
fit_ql
```

```text
── ComBat-refQL fit ─────────────────────────────────────────────

Data
  18,432 genes x 48 samples
  3 batches - 2 biological groups

Model
  QL GLM:         edgeR v4
  Normalization:  TMMwsp
  Transport:      negative-binomial mid-P
  Reference:      batch_A - explicitly selected

Result
  ✔ 17,904 adjusted
  ℹ 528 unchanged
  ✔ 0 failed

Runtime: 18.7 seconds
Use summary(x) for diagnostics.
```

Use `summary()` for the design, reference scores, model diagnostics, gene outcomes, recorded input actions, and timings:

```r
summary(fit_ql)
```

The adjusted count matrix and tidy diagnostic components remain available from the S7 object:

```r
adjusted_counts <- fit_ql@counts

fit_ql@specification
fit_ql@samples
fit_ql@gene_status
fit_ql@diagnostics
```

## Implementation

ComBat-refQL performs the workflow in four stages:

1. validates and aligns the count matrix and sample metadata;
2. calculates TMMwsp normalization factors and fits the edgeR model;
3. selects or validates the reference batch and calculates reference-counterfactual means;
4. transports supported non-reference counts through deterministic negative-binomial mid-P quantiles.

Gene chunks control temporary memory allocations during count transport:

```r
fit_ql <- combat_ref_ql(
  counts = counts,
  batch = metadata$batch,
  group = metadata$condition,
  chunk_size = 5000L
)
```

Adjusted counts are suitable for batch-aligned visualization, exploratory analysis, and downstream workflows requiring a count-scale matrix. Differential-expression inference is performed on the original counts with batch and biological variables represented in the statistical design.

## Citation

If you use ComBat-refQL, cite the package and the methods that motivated its implementation:

```r
citation("combatrefql")
```

### ComBat-ref

Zhang X. Highly effective batch effect correction method for RNA-seq count data. *Computational and Structural Biotechnology Journal*. 2025. [doi:10.1016/j.csbj.2024.12.010](https://doi.org/10.1016/j.csbj.2024.12.010)

### ComBat-seq

Zhang Y, Parmigiani G, Johnson WE. ComBat-seq: batch effect adjustment for RNA-seq count data. *NAR Genomics and Bioinformatics*. 2020;2(3):lqaa078. [doi:10.1093/nargab/lqaa078](https://doi.org/10.1093/nargab/lqaa078)

### edgeR

Run the following command to obtain the citations appropriate for the installed edgeR version:

```r
citation("edgeR")
```
