# ComBat-refQL

[![R CMD check](https://github.com/Thokas99/ComBat-refQL/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Thokas99/ComBat-refQL/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/Thokas99/ComBat-refQL/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/Thokas99/ComBat-refQL/actions/workflows/pkgdown.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

**Reference-batch correction for raw RNA-seq counts.**

ComBat-refQL combines TMMwsp normalization, a negative-binomial model,
counterfactual reference-batch prediction, and deterministic mid-P quantile
mapping. It preserves count-matrix dimensions, reports unsupported genes, and
records every automatic input repair.

> [!IMPORTANT]
> The QL-assisted pathway is experimental and diagnostic. For differential
> expression, analyze raw counts directly with an edgeR QL model containing
> batch and biological group.

## Table of contents

- [Why ComBat-refQL?](#why-combat-refql)
- [Installation](#installation)
- [Quick start](#quick-start)
- [What the result contains](#what-the-result-contains)
- [Command-line interface](#command-line-interface)
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
- **Safe by default:** never silently rounds fractional counts, aggregates
  duplicate genes, or accepts normalized expression values.
- **Practical:** offers both an R API and an informative command-line tool.

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

ComBat-refQL requires R 4.1 or newer and edgeR 4.2 or newer.

## Quick start

`counts` must contain genes in rows and samples in columns. The batch and group
vectors must describe those same samples.

```r
fit <- combat_ref(
  counts = counts,
  batch = metadata$batch,
  group = metadata$condition,
  reference_batch = "batch_A"
)

adjusted_counts <- fit$adjusted_counts
```

If `reference_batch` is omitted, ComBat-refQL selects one automatically and
returns the batch scores and local designs used in that decision. An explicit
reference is recommended when study design or external QC identifies the best
batch.

The defaults are:

| Setting | Default |
|---|---|
| Normalization | `TMMwsp` |
| Fit | negative-binomial baseline |
| Quantile mapping | deterministic mid-P |
| Reference | automatic unless specified |
| Messages | normal interactively, quiet in scripts |

### Experimental QL diagnostics

```r
fit_ql <- combat_ref_ql(
  counts,
  batch = metadata$batch,
  group = metadata$condition,
  reference_batch = "batch_A"
)
```

This adds edgeR quasi-likelihood diagnostics; corrected counts are still mapped
with a negative-binomial distribution. QL dispersion is never substituted into
`pnbinom()` or `qnbinom()`.

## What the result contains

```r
fit$adjusted_counts                 # non-negative, integer-valued doubles
fit$reference_batch                 # selected or requested reference
fit$gene_status                     # adjusted, unsupported, or failed genes
fit$diagnostics$input_actions       # every safe conversion performed
fit$diagnostics$reference_local_designs
summary(fit)
```

Set `keep_model = TRUE` only when fitted coefficients, means, and dispersions
are needed for detailed auditing; the default avoids retaining these large
objects.

## Command-line interface

```bash
combat-ref \
  --counts counts.tsv \
  --metadata metadata.csv \
  --batch-column batch \
  --group-column condition \
  --output adjusted_counts.tsv \
  --diagnostics-json diagnostics.json \
  --gene-status gene_status.tsv
```

A normal run prints compact Input, Configuration, Design, Normalization,
Reference selection, Correction, Warning, and Output sections. Useful modes:

| Option | Purpose |
|---|---|
| `--validate-only` | Check files, alignment, design, and normalization without fitting |
| `--verbose` | Add model details and stage timings |
| `--quiet` | Print errors only for pipeline use |
| `--no-color` | Produce plain-text logs |
| `--log-file PATH` | Save a plain-text run log |
| `--diagnostics-json PATH` | Write machine-readable diagnostics |
| `--gene-status PATH` | Write the per-gene status table |

Supported files are TSV, CSV, clearly delimited TXT, `.tsv.gz`, and `.csv.gz`.
Exit codes are `0` success, `2` input, `3` design, `4` fit/mapping, and `5`
output writing.

## Input safety

Safe, unambiguous conversions are announced and stored in
`fit$diagnostics$input_actions`:

- data frames and numeric character columns to a numeric matrix;
- a unique first gene-ID column to row names;
- integer storage to numeric doubles;
- batch and group vectors to factors;
- a one-level group to `NULL`;
- metadata reordering when sample identifiers match exactly;
- transposition only when orientation is unambiguous;
- case-insensitive method names to their canonical spelling.

Potentially consequential conversions require explicit opt-in:

```r
combat_ref(..., fractional_counts = "round")
combat_ref(..., duplicate_genes = "sum")
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
recorded. The baseline fit uses tagwise NB dispersion. With edgeR 4, the
experimental QL fit uses a global NB dispersion and stores it separately from
both baseline tagwise and moderated QL dispersions.

See the [statistical methodology vignette](vignettes/statistical-methodology.Rmd)
for equations, assumptions, and implementation details.

## Limitations

- The QL-assisted pathway has not passed broad simulation or external
  real-data validation and must not be considered production-ready.
- Automatic-reference stability is not yet bootstrapped; selection frequency
  is therefore not reported.
- Genes with an all-zero batch lack enough information for correction and are
  returned unchanged with an explicit status.
- Corrected counts do not propagate correction uncertainty into downstream
  differential-expression tests.
- Complete batch/group confounding cannot be repaired and is rejected.

## Documentation

- [Migration guide](MIGRATION.md)
- [Changelog](CHANGELOG.md)
- [Release notes](NEWS.md)
- [Implementation report](docs/IMPLEMENTATION_REPORT.md)
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
