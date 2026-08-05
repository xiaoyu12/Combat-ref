# combatrefql 0.4.0

- Add QL-aware empirical-Bayes shrinkage of exact source-versus-reference
  batch contrasts and make it the default mean pathway.
- Add hierarchical batch-scaled negative-binomial dispersion estimation and
  refit edgeR v4 with the resulting gene-by-sample NB dispersion matrix.
- Keep moderated QL uncertainty separate from the NB dispersion parameters
  used for count transport.
- Transport from source to reference NB distributions with distinct
  dispersions. `zero_policy = "model"` now applies the mid-P rule to zeros;
  `zero_policy = "preserve"` provides the earlier behavior.
- Add `mean_shrinkage` and `dispersion_model` controls for the four formal
  comparator and ablation configurations.
- Record contrast moderation, effective residual information, batch dispersion
  shifts, transport events, and automatic-reference margins in validated tidy
  diagnostics.
- Preserve reference counts exactly. Saturated source strata receive a neutral
  dispersion shift with `not_estimable` status; a zero-information reference
  raises a typed statistical error.
- Add deterministic statistical tests and a factorial validation study with
  pinned ComBat-ref and ComBat-seq comparators.
- Add the 0.4.0 statistical specification, migration guide, implementation
  report, validation plan, validation report, and updated methodology vignette.
- Restore pkgdown configuration with the `site/` destination so authored files
  under `docs/` remain separate from generated website output.

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
- Add the darkly Bootstrap theme for the pkgdown website.

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
- Use edgeR v4 quasi-likelihood fitting, model-based reference-counterfactual
  means, and deterministic negative-binomial mid-P transport.
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
