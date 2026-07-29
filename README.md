# ComBat-refQL

ComBat-refQL is a development R package for the baseline ComBat-ref RNA-seq
reference-batch correction algorithm. The package is being polished around the
existing method; quasi-likelihood correction is future work and is not
implemented.

## Table of contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Input requirements](#input-requirements)
- [Command-line usage](#command-line-usage)
- [Returned value](#returned-value)
- [Development status](#development-status)
- [Citation](#citation)
- [License](#license)

## Overview

`combatrefql` wraps the current ComBat-ref implementation as a small R library.
The exported function is `combat_ref()`, which returns corrected counts and a
compact result object.

The statistical correction remains the baseline ComBat-ref algorithm: negative
binomial fitting, dispersion estimation, automatic minimum-dispersion reference
selection, and negative-binomial quantile mapping are preserved. The package
does not run `glmQLFit()` and does not implement quasi-likelihood count
correction.

## Installation

```r
install.packages("pak")
pak::pak("Thokas99/ComBat-refQL")

library(combatrefql)
```

The package is not released on CRAN or Bioconductor.

## Quick start

```r
library(combatrefql)

counts <- matrix(
  c(10, 12, 30, 33, 8, 9, 28, 31,
    4, 5, 14, 16, 3, 4, 15, 17,
    50, 55, 120, 130, 48, 52, 125, 135,
    1, 2, 4, 5, 1, 2, 5, 6),
  nrow = 4,
  byrow = TRUE,
  dimnames = list(paste0("gene", 1:4), paste0("sample", 1:8))
)

metadata <- data.frame(
  batch = rep(c("batch_A", "batch_B"), each = 2, times = 2),
  condition = rep(c("control", "treated"), each = 4),
  row.names = colnames(counts)
)

fit <- combat_ref(
  counts = counts,
  batch = metadata$batch,
  group = metadata$condition,
  verbose = TRUE
)

adjusted_counts <- fit$adjusted_counts

fit
summary(fit)
```

Use an explicit reference batch when the study design requires one:

```r
fit <- combat_ref(
  counts = counts,
  batch = metadata$batch,
  group = metadata$condition,
  reference_batch = "batch_A",
  verbose = TRUE
)
```

## Input requirements

- `counts` must contain raw, non-negative integer RNA-seq counts.
- Genes are rows and samples are columns.
- `batch` must have one value per sample and align with `colnames(counts)`.
- `group` is optional and must also align with `colnames(counts)` when supplied.
- TPM, CPM, logCPM and other normalized expression values are not accepted.
- Batch and group must not be perfectly confounded.
- Each batch must contain enough samples for the baseline model fit.

## Command-line usage

```bash
combat-ref \
  --counts counts.tsv \
  --metadata metadata.csv \
  --batch-column batch \
  --group-column condition \
  --output adjusted_counts.tsv
```

The metadata CSV must use sample identifiers as row names in its first column.
The CLI reports metadata reordering when sample IDs match but order differs.

## Returned value

`combat_ref()` returns a `combatref_result` object with stable fields:

- `adjusted_counts`: corrected count matrix;
- `reference_batch`: selected or requested reference batch;
- `batch`: validated batch factor;
- `group`: validated group factor or `NULL`;
- `removed_genes`: genes that were all zero within at least one batch;
- `diagnostics`: compact counts, batch sizes, group sizes and reference details;
- `call`: matched function call.

Set `keep_model = TRUE` to include selected internal model quantities under
`model`.

## Development status

Current scope:

- baseline ComBat-ref package API;
- regression tests protecting legacy corrected counts;
- CLI wrapper around `combat_ref()`;
- restored source-only validation resources under `validation/`.

Deferred scope:

- quasi-likelihood correction;
- performance claims;
- public release until licensing is clarified.

## Citation

This package is derived from the original ComBat-ref research implementation by
Xiaoyu Zhang and collaborators. The original method and results should be cited
separately:

- Zhang X. and collaborators. ComBat-ref for RNA-seq batch correction. PubMed:
  <https://pubmed.ncbi.nlm.nih.gov/39802213/>

## License

Licensing is provisional. The package wrapper code and upstream-derived
statistical implementation are documented in `LICENSE` and `LICENSE_NOTES.md`.
No explicit upstream license was found during inspection, so redistribution
terms require manual clarification before public release.
