# ComBat-refQL implementation history

The current implementation exposes `combat_ref_ql()` as its fitting interface
and returns validated S7 fit and summary objects. The package uses TMMwsp
normalization, edgeR v4 model fitting, model-based reference-counterfactual
means, and deterministic negative-binomial mid-P transport.

The design contains named batch coefficients and locally estimable biological
terms. Automatic reference selection builds a fresh within-batch design for
each candidate and records the resulting score, retained terms, dropped terms,
missing levels, fallback, and residual degrees of freedom.

Transport runs by non-reference batch and gene chunk. Invalid mapped rows are
restored from the corresponding original rows and recorded in typed diagnostics.
The fit object retains adjusted counts, gene outcomes, model details, reference
diagnostics, input actions, warnings, and timings.

The 0.2.0 redesign established the current focused interface and replaced the
earlier multi-path API. Subsequent releases refined chunked transport, input
auditing, local reference designs, concise progress, and documentation.

Current verification commands are:

```r
devtools::document()
devtools::test()
devtools::check(args = "--no-manual")
```
