# ComBat-refQL Code Audit

Date: 2026-07-29

Scope: read-only audit of the current ComBat-refQL checkout. No scientific implementation was modified.

## 1. Repository Map

| File or path | Purpose | Type | Future package role | Recommendation |
|---|---|---|---|---|
| `README.md` | Project overview and current development status | Documentation | Keep | Rewritten for ComBat-refQL status and planned package direction. |
| `Combat_ref.R` | Current ComBat-ref count-correction implementation | Core code | Keep, modularize | Preserve until baseline tests exist; later split into model fitting, dispersion, reference-batch, and adjustment helpers. |
| `ComBat_seq.R` | Upstream-like ComBat-Seq implementation used for comparison | Core/reference code | Keep as comparison/internal baseline | Retain with upstream attribution; do not mix into new API without tests. |
| `helper_seq.R` | Shared helper functions, especially negative-binomial quantile mapping | Core code | Keep, modularize | Move reusable helpers into `R/` after tests; document rounding and zero behavior. |
| `batch_correct.R` | Command-line wrapper for count correction | CLI script | Retain initially, then move to `inst/scripts/` or deprecate | Use as baseline workflow fixture before replacing with R API. |
| `gfrn.mat` | Example count matrix | Data | Possibly `inst/extdata/` if small and licensed | Confirm provenance/license; use only small fixtures in package. |
| `gfrn_corrected.mat` | Generated corrected output | Generated output | Test fixture only if small and intentional | Prefer regenerate in tests; avoid committing generated outputs unless needed as baseline fixture. |
| `gfrn_samples.csv`, `gfrn_samples_new.csv` | Sample metadata | Data | Possibly `inst/extdata/` | Confirm which is canonical and document expected schema. |
| Root `batchFC*.png` | Simulation result figures | Generated output | Do not include in package source | Move to manuscript artifacts or remove from package branch. |
| `ComBat-seq.Rproj` | RStudio project file with old project name | Metadata | Optional | Rename only during package setup; do not change now. |
| `simulations/sim_DEpipe.R` | Main simulation pipeline | Simulation code | Move to `simulations/` or `inst/scripts/` outside package core | Preserve for validation; separate from package build. |
| `simulations/sim_DEpipe_helpers.R` | Simulation helpers and DE pipeline wrappers | Simulation/analysis code | Split; only generic utilities may move into package | Keep simulation-specific code outside `R/`. |
| `simulations/mcmc_seq.R` | MCMC/diagnostic helper using `rethinking` | Simulation/analysis code | Exclude from package core | Archive or document dependency if retained. |
| `simulations/plot_sim.R`, `plot_tpr_fpr.R`, `visualize.R`, `visualize_helpers.R` | Plotting and result aggregation | Generated-analysis code | Exclude from package core | Keep under `simulations/` or `analysis/`; remove hard-coded paths before reuse. |
| `simulations/test_sim_DEpipe.sh` | Shell simulation test runner | Test/script | Possible integration test template | Replace with small `testthat` fixtures; full simulation belongs outside routine checks. |
| `real_data_application/gfrn_application.R` | GFRN real-data batch-correction analysis | Real-data analysis | Exclude from package core | Move to `real-data/` or vignette only if data licensing permits. |
| `real_data_application/gfrn_DE.R` | GFRN downstream DE/GSEA analysis | Real-data analysis | Exclude from package core | Keep as validation/analysis, not exported package logic. |
| `real_data_application/gfrn_helpers.R` | DE method wrappers and analysis helpers | Analysis code | Mostly exclude | Some DE wrappers may inform validation tests, but should not be package API. |
| `real_data_application/*.png`, `*.csv`, `signature_data.rds` | Generated figures, pathway/signature data | Data/generated output | Exclude unless licensed and small | Confirm provenance; keep generated figures out of package source. |
| `nasa_data/nasa_application.R` | NASA real-data analysis | Real-data analysis | Exclude from package core | Keep as validation analysis if data permissions are clear. |
| `nasa_data/nasa.RData`, `*.png` | Real data and generated figure | Data/generated output | Exclude from package core | Confirm source/license; do not use as unit-test fixture. |
| `.DS_Store` | macOS metadata | Generated local artifact | None | Remove. |
| `.git/` | Git history | Metadata | Keep | Current branch is `main`, tracking `origin/main`. |

No `AGENTS.md`, `DESCRIPTION`, `NAMESPACE`, `LICENSE`, `CITATION`, `tests/`, `man/`, `vignettes/`, or package metadata files were found in the repository root.

## 2. Execution Flow

Current command-line workflow:

```text
input counts and metadata
-> validation and preprocessing
-> design-matrix construction
-> normalization
-> dispersion estimation
-> model fitting
-> reference-batch handling
-> count adjustment
-> output generation
```

1. Input counts and metadata
   - Entry point: `batch_correct.R`.
   - Important lines: `batch_correct.R:13` reads command-line arguments; `batch_correct.R:24` reads the count matrix; `batch_correct.R:27` reads sample metadata.
   - Inputs are expected as a tab-separated count matrix and a CSV-like sample file. The wrapper looks for a batch column by name matching.

2. Validation and preprocessing
   - File existence is checked in `batch_correct.R:20-21`.
   - Metadata column discovery is handled by local `find_column()` at `batch_correct.R:30`.
   - Core filtering in `Combat_ref.R` removes genes that are all zero within at least one batch (`Combat_ref.R:33`).
   - One-sample batches are rejected (`Combat_ref.R:28`).

3. Design-matrix construction
   - Core function: `ComBat_ref()` at `Combat_ref.R:19`.
   - Batch indices and counts are built around `Combat_ref.R:45-47`.
   - `build_design()` is defined inside `ComBat_ref()` at `Combat_ref.R:61`.
   - The design combines batch indicators with optional biological group design.
   - Confounding is checked by QR rank around `Combat_ref.R:73-74`.

4. Normalization
   - `Combat_ref.R:41` creates an edgeR `DGEList(counts=counts)`.
   - No `calcNormFactors()` call was found in `Combat_ref.R`; core correction appears to rely on edgeR offsets/library sizes from the `DGEList`.
   - `calcNormFactors()` appears only in downstream DE helpers (`simulations/sim_DEpipe_helpers.R:112`, `real_data_application/gfrn_helpers.R:43`).

5. Dispersion estimation
   - `Combat_ref.R` estimates batch-specific common/genewise dispersion through edgeR GLM dispersion helpers.
   - Reference batch selection is based on minimum common dispersion: `Combat_ref.R:95-97`.
   - `estimateDisp()` is not called in core correction. It is called in DE helper pipelines at `simulations/sim_DEpipe_helpers.R:125` and `real_data_application/gfrn_helpers.R:68`.

6. Model fitting
   - Core fitting uses `glmFit()` at `Combat_ref.R:139`.
   - Upstream-like comparison code uses `glmFit()` and `glmFit.default()` at `ComBat_seq.R:137` and `ComBat_seq.R:141`.
   - `glmQLFit()` is not called in core correction. It appears in DE helper code at `simulations/sim_DEpipe_helpers.R:126` and `real_data_application/gfrn_helpers.R:69`.

7. Reference-batch handling
   - `ComBat_ref.R` chooses the reference batch as the batch with the smallest estimated common dispersion (`Combat_ref.R:95-97`).
   - The reference batch is preserved; non-reference batches are mapped toward reference-batch mean/dispersion behavior.

8. Count adjustment
   - `Combat_ref.R:141-142` extracts batch coefficients and fitted means from the fitted model.
   - `Combat_ref.R:151-162` prepares old/new mean and dispersion quantities and calls `match_quantiles()`.
   - `helper_seq.R:47` defines `match_quantiles()`, which uses negative-binomial CDF/quantile mapping with `pnbinom()` and `qnbinom()`.
   - `Combat_ref.R:160-162` calls `match_quantiles(..., keep_zero=FALSE)`.

9. Output generation
   - `batch_correct.R:55` writes the corrected matrix with `write.table()`.
   - The current default output naming is script-controlled, not function-return controlled.

## 3. Statistical Audit

### Negative-Binomial Parameterization

The core algorithm uses edgeR negative-binomial GLM machinery. Helper code maps counts with base R negative-binomial CDF and quantile functions using `mu` and `size = 1 / phi`, where `phi` is dispersion. This implies variance `mu + phi * mu^2`.

### Normalization-Factor Handling

Core correction constructs a `DGEList` but does not call `calcNormFactors()`. The fitted GLM therefore depends on edgeR's default library-size offset behavior unless other offsets are constructed internally. DE analysis helpers do call TMM normalization before inference. This split must be documented because correction and downstream DE inference currently use different normalization handling.

### Design-Matrix Construction

The design matrix is built inside `ComBat_ref()` using batch indicators and optional biological group covariates. Rank checks detect some confounded designs. The design construction is embedded in the correction function, which makes it hard to test independently.

### Biological-Group Preservation

`group` is accepted by `ComBat_ref()` and used in design construction. The stated intent is to preserve biological group effects while adjusting batch effects. The exact behavior under missing groups, single-level groups, or batch/group confounding needs tests.

### Batch Coefficient Estimation

Batch coefficients are extracted from the fitted GLM after `glmFit()` (`Combat_ref.R:141`). The coefficient interpretation depends on the exact design coding. This needs regression tests before changing design construction.

### Dispersion Estimation and Pooling

The method selects a reference batch by minimum estimated common dispersion (`Combat_ref.R:95-97`) and uses common/genewise dispersion estimates for adjustment. The code comments describe adjustment toward the reference batch with minimum estimated dispersion.

### Reference-Batch Selection and Preservation

Reference batch is currently automatic, not user-specified. The planned API includes `reference_batch`; adding it will be result-changing unless implemented as an explicit override with tests.

### Fitted Means and Variance Assumptions

Fitted means come from edgeR GLM fitted values. Quantile mapping then uses old and new means/dispersion parameters. The validity of corrected counts depends on the mapping assumptions, not just on the GLM fit.

### Adjusted-Count Calculation

Adjustment is performed by negative-binomial quantile matching in `helper_seq.R:47-71`. `Combat_ref.R:160-162` calls this helper for non-reference batches. `keep_zero=FALSE` means zero counts may be remapped rather than preserved.

### Rounding, Truncation, and Integer Preservation

`qnbinom()` returns integer-valued quantiles. The output is count-like, but explicit tests should verify integer storage mode, non-negativity, and preservation of dimensions/dimnames.

### Zeros and Low-Count Genes

Genes all-zero in any batch are removed before fitting and later need restoration or documented omission behavior. The adjustment code uses `pmax(..., 1e-4)` for numerical stability around fitted means. Low-count and all-zero behavior is scientifically important and under-tested.

### Confounding and Balanced Designs

The code uses QR-rank checks for confounding. Additional tests are needed for perfectly confounded batch/group, single-group designs, insufficient replication, and unbalanced designs.

### Randomness and Reproducibility

Core correction appears deterministic. Simulation and real-data scripts set seeds in some places (`simulations/sim_DEpipe.R:12`, `real_data_application/gfrn_application.R:15`, `real_data_application/gfrn_DE.R:18`) but also contain random sampling in helpers. Reproducibility should be fixture-based.

### Current edgeR Use

Core correction:

- `DGEList()`: `Combat_ref.R:41`, `ComBat_seq.R:44`.
- `glmFit()`: `Combat_ref.R:139`, `ComBat_seq.R:137`.
- `glmFit.default()`: `ComBat_seq.R:141`.
- `estimateGLMCommonDisp()` / `estimateGLMTagwiseDisp()`: used in `Combat_ref.R` and `ComBat_seq.R`.

DE helper code:

- `DGEList()`: `simulations/sim_DEpipe_helpers.R:110`, `real_data_application/gfrn_helpers.R:42`.
- `calcNormFactors()`: `simulations/sim_DEpipe_helpers.R:112`, `real_data_application/gfrn_helpers.R:43`.
- `estimateDisp()`: `simulations/sim_DEpipe_helpers.R:125`, `real_data_application/gfrn_helpers.R:68`.
- `glmQLFit()`: `simulations/sim_DEpipe_helpers.R:126`, `real_data_application/gfrn_helpers.R:69`.

### Where QL Could and Could Not Be Used

QL could be evaluated for uncertainty-aware dispersion/model diagnostics, robust inference comparisons, and possibly as an internal estimation pathway if the correction algorithm formally defines how QL quantities alter fitted means, dispersion targets, or mapping uncertainty.

QL cannot be substituted mechanically for `glmFit()` to produce adjusted counts. `glmQLFit()` is an inference fit; it does not define the old-to-new negative-binomial quantile mapping used to generate corrected counts.

## 4. Software-Engineering Audit

| Severity | Location | Issue | Consequence | Recommended action | Could alter scientific results? |
|---|---|---|---|---|---|
| Critical | `Combat_ref.R:95-97` | Reference batch is chosen automatically as minimum common dispersion with no user override. | Users cannot reproduce a specified reference-batch analysis; API behavior may differ from study design. | Preserve current automatic mode as baseline; add tested `reference_batch` override later. | Yes. |
| Critical | `Combat_ref.R:139-162`, `helper_seq.R:47-71` | Count adjustment depends on fitted means, dispersion, and quantile mapping but is not isolated as a tested unit. | A small GLM/design change can silently alter corrected counts. | Extract adjustment only after baseline tests; test old/new NB parameter fixtures. | Yes. |
| High | `Combat_ref.R`, `ComBat_seq.R` | Core statistical code mixes validation, design construction, dispersion estimation, fitting, reference selection, and adjustment. | Hard to test or review individual statistical steps. | Modularize after baseline fixtures. | Yes. |
| High | `Combat_ref.R:73-74`, `ComBat_seq.R:84-85` | Confounding checks are embedded and not covered by tests. | Invalid designs may produce misleading corrections or unclear failures. | Add testthat cases for confounded and rank-deficient designs. | Yes. |
| High | `Combat_ref.R:33`, `ComBat_seq.R:36` | Genes all-zero in any batch are filtered internally. | Output dimensions or gene retention behavior may surprise users. | Define and test filtering/restoration policy. | Yes. |
| High | `batch_correct.R:24,27,55` | CLI performs direct file I/O and returns no structured object. | Cannot call library function directly or inspect diagnostics. | Keep wrapper as baseline; create R function returning structured result. | No, if wrapper delegates unchanged. |
| High | Repository root | No package metadata (`DESCRIPTION`, `NAMESPACE`), no tests, no documentation. | Not installable; no package checks or API contract. | Add minimal package skeleton after baseline tests. | No. |
| Medium | `batch_correct.R:30-38` | Metadata column discovery uses name matching. | Ambiguous column names can fail or select unexpected metadata. | Require explicit columns in future API; keep wrapper behavior tested. | Possibly. |
| Medium | `Combat_ref.R:61` | `build_design()` is nested inside `ComBat_ref()`. | Cannot unit-test design construction directly. | Extract to internal helper after baseline tests. | Possibly. |
| Medium | `Combat_ref.R:41` | Core correction lacks explicit normalization option. | Users may assume TMM normalization because DE helpers use it elsewhere. | Document current offset behavior; add explicit `normalization` argument later. | Yes. |
| Medium | `simulations/visualize.R:3-4` | Hard-coded absolute paths refer to old ComBat-seq locations. | Script is not portable. | Replace with project-relative paths when simulations are refactored. | No. |
| Medium | `nasa_data/nasa_application.R`, `real_data_application/*.R` | Analysis scripts source root files and write figures directly. | Repeated file I/O and generated output churn. | Move to `real-data/` workflows and keep generated output ignored. | No. |
| Medium | `simulations/sim_DEpipe.R` | Full simulation pipeline depends on many heavy packages and generated files. | Unsuitable for routine package tests. | Use small simulated fixtures for tests; keep full simulation as validation workflow. | No. |
| Medium | `simulations/sim_DEpipe_helpers.R:108-126`, `real_data_application/gfrn_helpers.R:41-69` | QL inference code is present only in DE wrappers, not correction. | Easy to confuse downstream QL DE testing with QL batch correction. | Keep terminology explicit; do not reuse as correction implementation without design. | Yes, if misused. |
| Low | `ComBat-seq.Rproj` | Project file still uses old ComBat-seq name. | Confusing project identity. | Rename during package skeleton stage. | No. |
| Low | `.DS_Store` | Local OS metadata is versioned/present. | Repository noise. | Remove and ignore. | No. |
| Low | Root PNGs and analysis PNGs | Generated figures are mixed with code. | Package source bloat and review noise. | Move generated outputs out of package source. | No. |
| Low | Many scripts | Naming alternates among `Combat`, `ComBat`, `ComBat_seq`, `Combat_ref`. | API and documentation inconsistency. | Standardize new public API as `combat_ref_ql()`. | No. |

## 5. Upstream Attribution and Licensing Audit

- Current remote: [Thokas99/ComBat-refQL](https://github.com/Thokas99/ComBat-refQL.git).
- Recent git history includes an upstream pointer mentioning [xiaoyu12/Combat-ref](https://github.com/xiaoyu12/Combat-ref), indicating this is likely a fork or copied derivative of the original ComBat-ref implementation.
- The old README credited Xiaoyu Zhang and linked the ComBat-ref publication.
- No `LICENSE` file was found. Do not invent licensing terms.
- No `CITATION` file was found. A future package should add `inst/CITATION` or `CITATION.cff` after confirming upstream citation requirements.
- Derivative files that copy or modify upstream code should carry explicit attribution comments or documentation.
- The package name `combatrefql` and method name `ComBat-refQL` distinguish this project from official ComBat-ref, but documentation must keep saying it is derivative until the relationship is legally and scientifically clarified.

Manual review needed:

- Upstream license and redistribution permission.
- Whether current data files can be redistributed in an R package.
- Whether figures/results are generated from licensed datasets.
- Exact author list and citation details for the original publication.

## 6. Proposed R-Package Architecture

Recommended target structure, adapted to this repository:

```text
ComBat-refQL/
|-- DESCRIPTION
|-- NAMESPACE
|-- LICENSE
|-- README.md
|-- R/
|   |-- combat_ref_ql.R
|   |-- combat_ref_original.R
|   |-- validate_inputs.R
|   |-- design_matrix.R
|   |-- estimate_dispersion.R
|   |-- fit_model.R
|   |-- adjust_counts.R
|   `-- utilities.R
|-- man/
|-- tests/
|   |-- testthat.R
|   `-- testthat/
|-- vignettes/
|-- inst/
|   |-- extdata/
|   `-- scripts/
|-- data-raw/
|-- simulations/
`-- real-data/
```

Package-role separation:

- Exported user functions: `combat_ref_ql()`, possibly `combat_ref_original()` for baseline comparison.
- Internal statistical functions: validation, design construction, dispersion estimation, model fitting, reference-batch selection, count adjustment.
- Simulations: keep under `simulations/`; not sourced by package.
- Example datasets: only small, licensed fixtures under `inst/extdata/`.
- Benchmarking/manuscript analyses: move real-data scripts and generated figures outside package core.
- CLI compatibility: keep `batch_correct.R` under `inst/scripts/` only after package API exists.

## 7. Proposed API

Candidate interface:

```r
combat_ref_ql(
  counts,
  batch,
  group = NULL,
  covariates = NULL,
  reference_batch = NULL,
  lib_size = NULL,
  normalization = c("TMM", "none"),
  dispersion = c("QL", "original"),
  return_model = FALSE
)
```

Assessment:

- `counts`: required; matrix-like non-negative integer counts, genes by samples.
- `batch`: required; length must equal `ncol(counts)`; at least two batches; sufficient replication.
- `group`: optional but scientifically important for preserving biological condition.
- `covariates`: useful, but must be checked for rank deficiency and alignment.
- `reference_batch`: should default to current automatic minimum-dispersion behavior only if `NULL`; explicit value must be validated.
- `lib_size`: optional; default should match edgeR library-size behavior. User-provided values need length and positivity checks.
- `normalization`: default is scientifically consequential. The safest first package release may expose `normalization = "none"` for exact baseline and add TMM only as a tested option.
- `dispersion`: should not default to `"QL"` until QL correction is specified and validated. Use `"original"` as the baseline path first.
- `return_model`: useful for diagnostics and reproducibility.

Proposed return class:

- A list with class `combatrefql_result`.
- Required fields: `counts`, `batch`, `group`, `reference_batch`, `removed_genes`, `normalization`, `dispersion_method`, `design`, `diagnostics`.
- Optional fields when `return_model = TRUE`: edgeR objects, fitted coefficients, fitted means, dispersion estimates, adjustment parameters.

Corrected counts should remain non-negative count-like values. Whether they must be stored as integer type should be tested and documented.

Original ComBat-ref behavior should be retained through `dispersion = "original"` or a separate `combat_ref_original()` wrapper until the QL method is validated.

## 8. Testing and Validation Plan

Small fixture tests before any algorithm change:

- invalid dimensions between `counts`, `batch`, `group`, and covariates;
- missing values in counts or metadata;
- duplicated sample names and duplicated gene IDs;
- single-batch input;
- one-sample batch and insufficient replication;
- zero-count genes and genes all-zero within one batch;
- all-zero samples;
- non-integer, negative, or non-numeric counts;
- perfectly confounded batch and group;
- invalid `reference_batch`;
- preservation of matrix dimensions and dimnames;
- deterministic repeated results;
- exact equivalence with current `batch_correct.R` output on a small fixture;
- exact equivalence with current `ComBat_ref()` output for the same fixture;
- simulation with known biological effect and known batch effect;
- type-I error, power, bias, false-positive rate, and batch-removal metrics;
- comparison between original and proposed QL pathways after QL is formally specified.

Regression fixtures should be small simulated matrices committed under `tests/testthat/fixtures/` or generated by deterministic test helpers. Large real datasets should not be used in routine unit tests.

## 9. Refactoring Roadmap

### Stage 1. Establish baseline behavior and regression tests

- Files affected: `tests/testthat.R`, `tests/testthat/`, small fixtures.
- Expected behavior: current `ComBat_ref()` and `batch_correct.R` outputs are captured without changing code.
- Tests required: exact corrected-count fixture, dimensions, dimnames, deterministic output.
- Scientific risks: low if no algorithm code changes.
- Completion criteria: tests fail if current output changes.

### Stage 2. Create minimal R-package skeleton

- Files affected: `DESCRIPTION`, `NAMESPACE`, `R/`, `.Rbuildignore`.
- Expected behavior: package installs and loads, but delegates to existing implementation.
- Tests required: `R CMD check`, namespace load, fixture test still passing.
- Scientific risks: low if function delegation is exact.
- Completion criteria: package can be installed locally and tests pass.

### Stage 3. Extract and document original algorithm

- Files affected: `R/combat_ref_original.R`, `R/adjust_counts.R`, `R/estimate_dispersion.R`, `R/fit_model.R`.
- Expected behavior: same corrected counts as current `Combat_ref.R`.
- Tests required: regression equivalence on fixtures and edge cases.
- Scientific risks: high because refactoring can alter design matrices or offsets.
- Completion criteria: byte- or value-equivalent outputs within documented tolerance.

### Stage 4. Replace script-level global state

- Files affected: `batch_correct.R`, `inst/scripts/`, `R/validate_inputs.R`.
- Expected behavior: CLI wrapper delegates to package function; R API returns structured result.
- Tests required: CLI integration test and direct function test produce same output.
- Scientific risks: medium if parsing or sample order changes.
- Completion criteria: sample/gene order preserved.

### Stage 5. Add validation and structured outputs

- Files affected: `R/validate_inputs.R`, `R/combat_ref_ql.R`.
- Expected behavior: invalid inputs fail early with clear messages.
- Tests required: invalid dimensions, missing values, confounding, zero counts.
- Scientific risks: medium because previous permissive inputs may now fail.
- Completion criteria: all planned input checks covered.

### Stage 6. Implement QL as separate internal pathway

- Files affected: `R/fit_model.R`, `R/estimate_dispersion.R`, new QL-specific helper.
- Expected behavior: original pathway unchanged; QL pathway opt-in.
- Tests required: original regression tests plus QL smoke tests and diagnostics.
- Scientific risks: critical until method is specified and validated.
- Completion criteria: QL quantities and count-adjustment role documented.

### Stage 7. Validate original and QL side by side

- Files affected: `simulations/`, `tests/testthat/`, validation reports.
- Expected behavior: comparison metrics are reproducible and interpretable.
- Tests required: type-I error, power, FPR, bias, batch-removal metrics.
- Scientific risks: critical for any performance claims.
- Completion criteria: validation report supports or rejects QL pathway claims.

### Stage 8. Optimize performance after correctness

- Files affected: core `R/` functions.
- Expected behavior: same outputs faster or with lower memory use.
- Tests required: regression equivalence and benchmark scripts.
- Scientific risks: medium if vectorization changes numerical behavior.
- Completion criteria: documented speed/memory improvement with unchanged results.

### Stage 9. Add vignettes, CI, package checks, documentation

- Files affected: `vignettes/`, `.github/workflows/`, `man/`, `README.md`.
- Expected behavior: examples run and documentation matches implementation.
- Tests required: `R CMD check`, vignette build, CI matrix.
- Scientific risks: low if docs do not make unvalidated claims.
- Completion criteria: clean package check.

### Stage 10. Prepare first development release

- Files affected: `NEWS.md`, `DESCRIPTION`, citation/license files, release docs.
- Expected behavior: installable development release with explicit limitations.
- Tests required: full package checks and fixture validation.
- Scientific risks: medium if release implies validated QL behavior too early.
- Completion criteria: release notes clearly separate original baseline from experimental QL.
