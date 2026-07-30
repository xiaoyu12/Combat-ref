# Using the corrected 0.1.0 pre-release

Defaults remain TMMwsp, baseline NB fitting, and mid-P mapping. QL users will
see a changed global NB dispersion because tagwise dispersion is no longer
passed to `glmQLFit()`.

Safe input coercions now occur automatically and are available at:

```r
fit$diagnostics$input_actions
```

Use `verbosity = "quiet"` to suppress messages without losing the audit.
Fractional counts and duplicate genes remain errors unless explicitly enabled
with `fractional_counts = "round"` or `duplicate_genes = "sum"`.

CLI users should replace `--diagnostics` with `--diagnostics-json`; optional
`--gene-status`, `--log-file`, `--validate-only`, `--no-color`, and
`--input-scale` are available. Failure codes are documented in `--help`.
