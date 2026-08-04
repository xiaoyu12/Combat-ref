# Scientific status

ComBat-refQL is experimental. Its software invariants are tested, but its
statistical operating characteristics have not been established broadly. This
document separates implemented behaviour from scientific claims that still
require evidence.

## Tested behaviour

- validated S7 fit and summary construction;
- input alignment, audited fractional rounding, and typed failures;
- robust edgeR v4 QL fitting with `legacy = FALSE`;
- strict separation of global NB dispersion from QL dispersions;
- deterministic NB mid-P transport and observed-zero preservation;
- chunk-equivalent mapping output and complete-row rollback after failures;
- biologically adjusted local reference designs, recorded dropped terms, and
  deterministic selection for fixed inputs;
- Linux and macOS package checks and clean installation.

These are implementation properties, not evidence of biological validity.

## Appropriate uses

- exploratory reference-batch harmonization of raw RNA-seq counts;
- method development where the chosen reference is justified independently;
- reproducible evaluation of deterministic count transport;
- inspecting how locally estimable biological terms affect automatic reference
  scores.

## Inappropriate uses

- treating corrected counts as uncertainty-free differential-expression data;
- using automatic selection as proof of batch quality;
- presenting the package as validated, production-ready, or superior;
- supplying TPM, CPM, log-expression, or other normalized matrices;
- correcting a design with complete batch–biology confounding.

For differential expression, analyse the original raw counts with batch and
the biological variables of interest in the edgeR QL design.

## Statistical limitations

- correction uncertainty is not propagated into downstream inference;
- genes without positive count support in every batch remain unchanged;
- zero preservation is a deliberate release constraint and changes the ideal
  quantile map at zero;
- one global underlying NB dispersion supplies the transport distribution;
- input-scale detection is necessarily heuristic;
- strong sparsity, composition shifts, outliers, and small batches need broader
  operating-characteristic evaluation.

## Automatic-reference limitations

Automatic scoring fits a fresh model within each candidate batch. It includes
the biological group and covariates only when they are locally variable,
linearly estimable, and leave positive residual degrees of freedom. Otherwise,
terms are dropped and an intercept-only fallback is used when necessary.

The smallest valid global NB dispersion is a heuristic score, not a formal
measure of technical quality. Candidate batches may differ in sample size,
composition, missing levels, or estimable covariates. Selection stability is
not bootstrapped, and the scoring rule is not claimed to be optimal or
validated. Prefer an explicit reference supported by study design and external
QC when available.

## Validation still required

- predefined simulations spanning mean, dispersion, balance, confounding,
  sparsity, composition bias, and batch size;
- type-I error, FDR, power, and biological-signal preservation assessments;
- external real-data comparisons with documented technical controls;
- reference-selection stability analysis, including bootstrap procedures;
- runtime and peak-memory profiling across realistic scales;
- comparisons with ComBat-seq and direct model-based alternatives without
  selective reporting.

The current [validation report](validation/VALIDATION_REPORT.md) contains only
an implementation smoke simulation and its explicit limitations.

## Requirements for production-readiness claims

A future production-readiness claim would require independent real-data
validation, broad preregistered simulations, documented failure boundaries,
stable reference selection, reproducible performance measurements, and clear
guidance for downstream inference. None of those claims is made in version
0.3.0.
