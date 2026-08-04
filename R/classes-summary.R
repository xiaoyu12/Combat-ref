#' ComBat-refQL summary
#'
#' A lightweight validated S7 summary returned by `summary(fit)`.
#' @param overview One-row size overview.
#' @param specification One-row method specification.
#' @param design One-row design diagnostics.
#' @param batches Reference scores by batch.
#' @param gene_outcomes Gene outcome counts and proportions.
#' @param ql QL and NB dispersion diagnostics.
#' @param input_actions Audited input transformations.
#' @param warnings Stored warning table.
#' @param timings Stage timing table.
#' @export
CombatRefQLSummary <- S7::new_class("CombatRefQLSummary", properties = list(
  overview = S7::class_data.frame, specification = S7::class_data.frame,
  design = S7::class_data.frame, batches = S7::class_data.frame,
  gene_outcomes = S7::class_data.frame, ql = S7::class_data.frame,
  input_actions = S7::class_data.frame, warnings = S7::class_data.frame,
  timings = S7::class_data.frame
), constructor = function(overview, specification, design, batches, gene_outcomes,
                          ql, input_actions, warnings, timings) {
  S7::new_object(S7::S7_object(), overview = overview, specification = specification,
    design = design, batches = batches, gene_outcomes = gene_outcomes, ql = ql,
    input_actions = input_actions, warnings = warnings, timings = timings)
}, validator = function(self) {
  errors <- character()
  if (nrow(self@overview) != 1L) errors <- c(errors, "@overview must have one row")
  if (!all(c("batch", "samples", "reference_score", "selected") %in% names(self@batches))) errors <- c(errors, "@batches is missing required columns")
  if (!all(c("outcome", "genes", "proportion") %in% names(self@gene_outcomes))) errors <- c(errors, "@gene_outcomes is missing required columns")
  if (length(errors)) errors
})
