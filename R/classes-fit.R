#' ComBat-refQL fitted result
#'
#' A validated S7 result from [combat_ref_ql()]. Properties use `@` access.
#' @param counts Adjusted count matrix.
#' @param specification One-row method specification.
#' @param samples Tidy sample metadata.
#' @param gene_status Tidy gene outcome metadata.
#' @param diagnostics Structured diagnostic tables.
#' @param call Matched fitting call.
#' @param model Compact model information.
#' @export
CombatRefQLFit <- S7::new_class("CombatRefQLFit", properties = list(
  counts = S7::class_any, specification = S7::class_data.frame,
  samples = S7::class_data.frame, gene_status = S7::class_data.frame,
  diagnostics = S7::class_list, call = S7::class_call, model = S7::class_list
), constructor = function(counts, specification, samples, gene_status, diagnostics, call, model) {
  S7::new_object(S7::S7_object(), counts = counts, specification = specification,
    samples = samples, gene_status = gene_status, diagnostics = diagnostics,
    call = call, model = model)
}, validator = function(self) {
  errors <- character(); counts <- self@counts
  if (!is.matrix(counts) || !is.numeric(counts)) errors <- c(errors, "@counts must be a numeric matrix")
  if (is.matrix(counts)) {
    if (ncol(counts) != nrow(self@samples)) errors <- c(errors, "count columns must match sample rows")
    if (nrow(counts) != nrow(self@gene_status)) errors <- c(errors, "count rows must match gene-status rows")
    if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0)) errors <- c(errors, "adjusted counts must be finite and non-negative")
    if (any(abs(counts - round(counts)) > sqrt(.Machine$double.eps)) || !is.double(counts)) errors <- c(errors, "adjusted counts must be integer-valued doubles")
    if (anyDuplicated(colnames(counts)) || anyDuplicated(rownames(counts))) errors <- c(errors, "sample and gene identifiers must be unique")
    if (!identical(colnames(counts), self@samples$sample)) errors <- c(errors, "sample identifiers are not aligned")
    if (!identical(rownames(counts), self@gene_status$gene)) errors <- c(errors, "gene identifiers are not aligned")
  }
  if (!all(c("gene", "status", "reason", "adjusted", "zero_batches", "fit_valid", "mapping_valid") %in% names(self@gene_status))) errors <- c(errors, "@gene_status is missing required columns")
  if (!all(c("batches", "mapping", "mapping_failures", "timing",
    "input_actions", "warnings", "outcomes", "ql", "batch_contrasts",
    "dispersion") %in% names(self@diagnostics)))
    errors <- c(errors, "@diagnostics is missing required tables")
  batch_columns <- c("batch", "samples", "reference_score", "selected",
    "selection_method", "local_formula", "local_columns", "dropped_columns",
    "missing_levels", "fallback", "fallback_reason", "residual_df",
    "effective_residual_information", "second_best", "score_margin")
  if (is.data.frame(self@diagnostics$batches) &&
      !all(batch_columns %in% names(self@diagnostics$batches)))
    errors <- c(errors, "reference diagnostics are missing required columns")
  contrast_columns <- c("gene", "source_batch", "reference_batch",
    "raw_delta", "ql_variance", "ql_se", "prior_mean", "prior_variance",
    "shrinkage_weight", "posterior_delta", "valid", "status", "converged")
  if (!all(contrast_columns %in% names(self@diagnostics$batch_contrasts)))
    errors <- c(errors, "batch-contrast diagnostics are missing required columns")
  dispersion_columns <- c("gene", "batch", "reference_batch",
    "baseline_dispersion", "reference_dispersion", "source_dispersion",
    "d_batch", "multiplier", "design_residual_information",
    "weighted_residual_information", "precision_weight", "precision_formula",
    "estimable", "status", "converged", "iterations")
  if (!all(dispersion_columns %in% names(self@diagnostics$dispersion)))
    errors <- c(errors, "dispersion diagnostics are missing required columns")
  mapping_columns <- c("source_batch", "reference_batch", "genes_attempted",
    "genes_adjusted", "genes_failed", "source_dispersion_min",
    "source_dispersion_max", "target_dispersion_min", "target_dispersion_max",
    "zero_to_positive", "positive_to_zero", "tail_events",
    "reference_invariant", "seconds")
  if (!all(mapping_columns %in% names(self@diagnostics$mapping)))
    errors <- c(errors, "mapping diagnostics are missing required columns")
  if (length(self@specification$reference_batch) != 1L || !self@specification$reference_batch %in% self@samples$batch) errors <- c(errors, "reference batch is invalid")
  if (!all(c("method", "normalization", "transport", "zero_handling") %in%
      names(self@specification))) errors <- c(errors,
    "@specification is missing the integrated 0.4.0 method description")
  if (!all(c("coefficients", "nb_dispersion", "reference_dispersion",
      "batch_dispersion_shifts", "batch_dispersion_multipliers",
      "ql_dispersion", "posterior_ql_dispersion", "design",
      "design_matrix") %in% names(self@model)))
    errors <- c(errors, "@model is missing required 0.4.0 quantities")
  outcomes <- self@diagnostics$outcomes
  if (is.data.frame(outcomes) && nrow(outcomes)) {
    actual <- table(factor(self@gene_status$status, levels = outcomes$outcome))
    if (!identical(as.integer(actual), as.integer(outcomes$genes)))
      errors <- c(errors, "gene outcome totals do not match actual statuses")
  }
  if (length(errors)) errors
})
