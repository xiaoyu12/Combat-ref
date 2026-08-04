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
  if (!all(c("batches", "mapping", "mapping_failures", "timing", "input_actions", "warnings", "outcomes", "ql") %in% names(self@diagnostics))) errors <- c(errors, "@diagnostics is missing required tables")
  if (length(self@specification$reference_batch) != 1L || !self@specification$reference_batch %in% self@samples$batch) errors <- c(errors, "reference batch is invalid")
  outcomes <- self@diagnostics$outcomes
  if (is.data.frame(outcomes) && nrow(outcomes)) {
    actual <- table(factor(self@gene_status$status, levels = outcomes$outcome))
    if (!identical(as.integer(actual), as.integer(outcomes$genes)))
      errors <- c(errors, "gene outcome totals do not match actual statuses")
  }
  if (length(errors)) errors
})
