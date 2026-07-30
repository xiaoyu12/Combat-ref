# combatrefql 0.1.0

## Defaults

- The default normalization is TMMwsp and the default mapping is deterministic mid-P.

## Statistical corrections

- Fixed the missing gene baseline when `group = NULL` with one named coefficient per batch.
- Validation and fitting now use the same full-rank design.
- Replaced positional coefficient subtraction with counterfactual prediction.
- Removed the discontinuous `gamma < -0.2` fitted-mean cap.

## New methods

- Added experimental `fit_method = "ql_assisted"` and `combat_ref_ql()`.
- QL dispersion is diagnostic only and remains separate from NB dispersion.

## Normalization and offsets

- Added TMMwsp, TMM, and none; library sizes and factors are calculated once before filtering.

## Reference selection

- Reference scores and an explicit unresolved-stability field are returned.

## Quantile mapping

- Added mid-P, upper-CDF, and explicitly seeded randomized mapping with safe probability bounds.

## Diagnostics

- Added gene status, design, offset, dispersion, version, and numerical-tail diagnostics.

## Performance

- `keep_model = FALSE` no longer retains fitted-mean, coefficient, or design-copy model payloads.

## Tests and validation

- Added direct tests for TMMwsp, null-group design, counterfactual prediction,
  numerical mapping, QL/NB separation, and unsupported genes.

## Documentation

- Added migration, methodology, validation, changelog, and implementation reports.

## Known limitations

- Reference resampling and validated QL coefficient shrinkage are deferred.
