# ComBat-refQL 0.3.0 implementation report

The package has one orchestrator, one edgeR v4 QL fit, and one deterministic NB
mid-P kernel. Input preparation, normalization, support assessment, design
construction, biologically adjusted reference scoring, fitting, transport,
diagnostics, result construction, S7 validation, and reporting are separated by
responsibility. Numerical kernels do not print.

Each automatic-reference candidate receives a fresh local design. Locally
constant, redundant, or residual-df-exhausting biological terms are recorded
and omitted; intercept-only fallback is explicit in the batch diagnostics.

Peak-memory work is structural: observed and reference means are calculated
inside each `chunk_size` gene block for the current non-reference batch;
offsets, design, and normalization are reused; zero support is computed once;
and observed/reference mean matrices are omitted from the final object.

`model` retains coefficients, one NB dispersion, QL diagnostic vectors, and
the design needed for inspection. It does not retain input counts, fitted means,
counterfactual means, or transport intermediates.

Transport validates every mapped row. Any gene with invalid counterfactual
means or adjusted counts is marked failed, restored from the complete original
input row, and included in tidy failure and warning diagnostics. Outcome totals
are derived from these final gene statuses.
