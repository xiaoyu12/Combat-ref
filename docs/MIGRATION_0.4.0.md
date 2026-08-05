# Migrating to ComBat-refQL 0.4.0

## Changed defaults

Version 0.4.0 defaults to `mean_shrinkage = "ql_eb"`,
`dispersion_model = "batch_scaled"`, and `zero_policy = "model"`. The fitted
mean correction therefore uses moderated source-versus-reference contrasts,
and transport can use distinct source and reference NB dispersions.

Adjusted counts can differ from 0.3.1 because QL uncertainty now controls the
strength of each gene correction, estimated batch multipliers alter the source
and target count distributions, and zeros follow the model-based mid-P rule.

## Compatibility configuration

Use the statistical 0.3.1-style comparator with:

```r
combat_ref_ql(
  counts, batch, group,
  mean_shrinkage = "none",
  dispersion_model = "shared",
  zero_policy = "preserve"
)
```

The A-only and B-only ablations set only one of `mean_shrinkage = "ql_eb"` or
`dispersion_model = "batch_scaled"`. These settings support validation and
method comparison; the combined pathway is the default.

## Diagnostics

`fit@diagnostics$batch_contrasts` contains raw contrasts, QL standard errors,
prior parameters, weights, posterior contrasts, and status. The
`dispersion` table contains baseline/reference/source NB dispersions,
batch shifts, multipliers, design-only and weighted residual information,
precision, convergence, and estimability. Mapping diagnostics add source and
target dispersion ranges, zero transitions, tail events, failures, and the
reference-invariance check. Reference diagnostics include the runner-up and
selection margin.

Adjusted counts support batch-aligned visualization and exploratory workflows.
Differential-expression inference uses original counts with batch and biology
represented in the statistical design.
