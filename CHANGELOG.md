# Changelog

## 0.1.0 — auditable development release

| Area | Old behavior | New behavior | Consequence and migration | Code/tests |
|---|---|---|---|---|
| Null design | Reference rows had no gene baseline | `~ 0 + batch` foundation | Finite identifiable means; use legacy wrapper for old output | `R/build_design.R`, design tests |
| Depth | Recomputed/inconsistent library sizes | One pre-filter TMMwsp calculation | Common offsets across all stages | baseline and normalization tests |
| Prediction | Positional coefficient subtraction | Named counterfactual design | Preserves group/covariates and factor semantics | `R/combat_ref_baseline.R` |
| Stabilization | Discontinuous `-0.2` cap | Removed | No arbitrary threshold | counterfactual tests |
| Mapping | Upper endpoint with unchanged tails | Mid-P with clamping | Stable deterministic tails | quantile tests |
| QL | Not implemented | Diagnostic `glmQLFit()` layer | NB mapping still uses NB dispersion | QL separation test |
| Unsupported genes | Returned unchanged silently | Returned unchanged with status | Users can exclude or report them | gene-status test |
| Compatibility | Main API was historical | `combat_ref_legacy()` | Frozen fixtures remain unchanged | legacy tests |

No performance superiority or validation claim is made.
