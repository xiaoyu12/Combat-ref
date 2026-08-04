# ComBat-refQL 0.2.0 implementation report

The package has one orchestrator, one edgeR v4 QL fit, and one deterministic NB
mid-P kernel. Input preparation, design construction, reference scoring,
fitting, transport, diagnostics, S7 validation, and printing are separated by
responsibility. Numerical kernels do not print.

Peak-memory work is structural: reference means are calculated only for the
current non-reference batch; genes are transported in `chunk_size` blocks;
offsets, design, and normalization are reused; zero support is computed once;
and observed/reference mean matrices are omitted from the final object.

`model` retains coefficients, one NB dispersion, QL diagnostic vectors, and
the design needed for inspection. It does not retain input counts, fitted means,
counterfactual means, or transport intermediates.

