# combatrefql 0.4.0

- Integrate QL-aware empirical-Bayes moderation of exact
  source-versus-reference batch contrasts into the single `combat_ref_ql()`
  method.
- Estimate a moderated gene-dispersion baseline and reference-centred
  batch-specific NB dispersion multipliers using continuous effective residual
  information.
- Refit edgeR v4 with gene-by-sample NB dispersions while keeping moderated QL
  uncertainty separate from NB count-distribution parameters.
- Transport non-reference observations from their fitted source mean and
  dispersion to reference-counterfactual mean and dispersion with deterministic
  NB mid-P mapping. Apply the same model rule to zero and positive counts.
- Preserve reference counts exactly and record saturated or weak-information
  dispersion strata explicitly.
- Record contrast moderation, reference-selection margins, batch dispersion,
  residual information, transport events, gene outcomes, and convergence in
  validated tidy diagnostics.
- Expand mathematical unit tests and factorial simulation validation. Add
  reproducible GFRN and NASA GeneLab workflows with pinned ComBat-ref and
  ComBat-seq comparisons.
- Consolidate the public interface around one statistical pathway and refresh
  the README, methodology vignette, validation report, and pkgdown site.

# combatrefql 0.3.2

- Reframe the package around its R interface for reference-batch adjustment.
- Align README, package metadata, generated documentation, citation metadata,
  terminology, progress messages, printing, and summaries.
- Add explicit ComBat-ref attribution and preserve its upstream MIT notice
  separately from the package license.
- Use runtime package-version lookup for startup output.

# combatrefql 0.3.1

- Remove the interactive progress bar while retaining concise stage messages
  and durable verbose mapping updates.
- Add a Bootstrap theme for the pkgdown website.

# combatrefql 0.3.0

- Add biologically adjusted within-batch reference scoring.
- Record locally retained and omitted biological terms, missing factor levels,
  fallbacks, residual degrees of freedom, and reference scores.
- Strengthen mapping rollback, input auditing, and package documentation.

# combatrefql 0.2.1

- Validate NB and QL dispersion separation throughout model fitting and
  transport.
- Improve chunked counterfactual prediction, deterministic transport tests,
  CI checks, and validation smoke testing.

# combatrefql 0.2.0

- Establish `combat_ref_ql()` as the focused fitting interface.
- Use edgeR v4 quasi-likelihood fitting, model-based
  reference-counterfactual means, and deterministic negative-binomial mid-P
  transport.
- Introduce validated S7 `CombatRefQLFit` and `CombatRefQLSummary` objects,
  typed conditions, audited fractional-count handling, concise progress, and
  tidy diagnostics.
- Separate NB dispersion used by count probabilities from QL dispersion used
  for model uncertainty.
- Add TMMwsp normalization, biological-group and sample-covariate modelling,
  automatic or explicit reference selection, and chunked transport.

# combatrefql 0.1.0

- Establish the independently maintained ComBat-refQL package history.
- Provide an initial reference-batch count-adjustment implementation with
  edgeR-based count modelling and counterfactual prediction.
