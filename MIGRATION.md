# Migrating to 0.1.0

| Before | 0.1.0 | Action |
|---|---|---|
| No normalization | TMMwsp default | Set `normalization = "none"` for an unnormalized repaired fit. |
| Historical implementation | Repaired baseline | Call `combat_ref_legacy()` for exact historical output. |
| Broken `group = NULL` mean model | Named batch baselines | Expect corrected fitted means and numerical results. |
| Upper-CDF mapping | Mid-P mapping | Set `quantile_method = "upper"` for the alternative repaired mapping. |
| Positional batch coefficients and hard cap | Counterfactual prediction | Expect continuous, design-invariant changes. |
| Silent unsupported genes | `gene_status` | Inspect `unsupported_zero_batch` before interpreting output. |

QL diagnostics are opt-in with `combat_ref_ql()` or
`fit_method = "ql_assisted"`. They do not define a probability distribution.
