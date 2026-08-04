# combatrefql 0.2.0

This release intentionally breaks compatibility to make ComBat-refQL a focused
QL-first package.

- `combat_ref_ql()` is the sole fitting API and edgeR v4 QL is the sole model.
- Counterfactual means from `glmQLFit(..., robust = TRUE, legacy = FALSE)` are
  transported with one deterministic negative-binomial mid-P method.
- `CombatRefQLFit` and `CombatRefQLSummary` replace all old S3 result classes.
- CLI output is staged and non-redundant, with quiet, normal, and verbose modes.
- Fractional counts are rounded by default with a structured lossy-input audit.
- Typed conditions and tidy diagnostics make failures and results inspectable.
- Counterfactual prediction is batch-local and transport is gene-chunked.
- Baseline and legacy fits, method-selection arguments, randomized mapping,
  seeds, plotting, JSON diagnostics, the standalone CLI, deprecated aliases,
  compatibility wrappers, and migration machinery were removed.

Negative-binomial dispersion and QL dispersion are separate quantities. No QL
or posterior QL dispersion is passed to an NB probability function.
