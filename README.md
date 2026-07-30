# ComBat-refQL

Reference-batch correction for raw RNA-seq counts. Version 0.1.0 provides a
repaired negative-binomial ComBat-ref baseline and an explicitly experimental
edgeR quasi-likelihood diagnostic layer. QL dispersion is never used as an NB
distribution parameter.

## Contents

- [Install](#install)
- [Quick start](#quick-start)
- [Methods](#methods)
- [Output](#output)
- [Limitations](#limitations)
- [License](#license)

## Install

```r
remotes::install_github("Thokas99/ComBat-refQL")
```

## Quick start

```r
library(combatrefql)

fit <- combat_ref(
  counts,
  batch = metadata$batch,
  group = metadata$condition,
  reference_batch = "batch_A"
)
fit$adjusted_counts
fit$gene_status
```

The defaults are `normalization = "TMMwsp"`, `fit_method = "baseline"`, and
`quantile_method = "midp"`. Automatic reference selection is available when
`reference_batch = NULL`, but an explicit reference is more reproducible.

Use `combat_ref_ql()` to add experimental QL diagnostics while retaining NB
count mapping.

## Methods

The repaired model uses named batch coefficients plus biological covariates,
one set of effective library sizes calculated before filtering, and
counterfactual prediction with each sample's batch set to the reference batch.
Unsupported genes remain unchanged and are identified in `gene_status`.

For differential expression, analyze raw counts directly with an edgeR QL model
including batch and biological group. Adjusted counts do not propagate batch
correction uncertainty.

## Output

`combat_ref()` returns adjusted counts, selected reference, resolved methods,
gene-level status, compact diagnostics, and the call. `keep_model = TRUE` also
returns coefficients, observed/reference means, design matrices, and separate
NB/QL dispersions.

The CLI accepts the same method choices and can write a diagnostics TSV; run
`combat-ref --help` for options.

## Limitations

- QL support is experimental and diagnostic-only; it is not a QL count distribution.
- Automatic reference scores are dispersion-based point estimates in 0.1.0;
  stability probabilities require future resampling validation.
- Genes all zero in any batch cannot be fitted and are returned unchanged.
- Corrected-count performance has not been externally validated.

## License

The package is distributed under the MIT license in [LICENSE](LICENSE).
