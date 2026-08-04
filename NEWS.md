# combatrefql 0.3.1

- Remove the interactive progress bar while retaining concise stage messages
  and durable verbose mapping updates.
- Use the darkly Bootstrap theme for the pkgdown website.

# combatrefql 0.3.0

- Automatic reference scoring now adjusts for locally estimable biological
  groups and covariates.
- Local reference designs record dropped terms, missing levels, fallbacks, and
  residual degrees of freedom.
- The main fitting orchestrator and result construction were simplified without
  changing the public API.
- The README was shortened and detailed scientific limitations moved to
  `SCIENTIFIC_STATUS.md`.

# combatrefql 0.2.1

- Calculate observed and counterfactual means inside each gene chunk during
  transport, reducing peak temporary memory.
- Detect numerical mapping failures per gene, retain each failed gene's complete
  original row, and report typed structured failure diagnostics.
- Derive outcome totals from final gene statuses and validate their consistency
  in `CombatRefQLFit`.
- Refine normalized-expression detection so suspicious dense or equal-library
  matrices are rejected while plausible sparse fractional estimated counts can
  be audited and rounded.
- Streamline the README and validate every local link.

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
