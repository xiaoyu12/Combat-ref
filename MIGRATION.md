# Using 0.1.0

The supported API is `combat_ref()` with these defaults:

```r
combat_ref(
  counts,
  batch,
  group = NULL,
  normalization = "TMMwsp",
  fit_method = "baseline",
  quantile_method = "midp"
)
```

Use `normalization = "none"` for unnormalized effective library sizes,
`quantile_method = "upper"` for upper-CDF mapping, and `combat_ref_ql()` for
experimental QL diagnostics. Inspect `gene_status` before interpreting genes
that are all zero within a batch.
