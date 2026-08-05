#' @export
S7::method(summary, CombatRefQLFit) <- function(object, ...) {
  contrast <- object@diagnostics$batch_contrasts
  shrinkage <- do.call(rbind, lapply(split(contrast, contrast$source_batch),
    function(x) data.frame(source_batch = x$source_batch[1L],
      prior_mean = x$prior_mean[1L], prior_variance = x$prior_variance[1L],
      median_weight = stats::median(x$shrinkage_weight, na.rm = TRUE),
      valid_genes = sum(x$valid), status = x$status[1L])))
  dispersion <- object@diagnostics$dispersion
  dispersion <- do.call(rbind, lapply(split(dispersion, dispersion$batch),
    function(x) data.frame(batch = x$batch[1L], d_batch = x$d_batch[1L],
      multiplier = x$multiplier[1L],
      design_residual_information = x$design_residual_information[1L],
      median_weighted_information = stats::median(x$weighted_residual_information),
      total_precision = sum(x$precision_weight), status = x$status[1L],
      converged = x$converged[1L], iterations = x$iterations[1L])))
  CombatRefQLSummary(
    overview = data.frame(genes = nrow(object@counts), samples = ncol(object@counts), batches = length(unique(object@samples$batch)), groups = length(unique(stats::na.omit(object@samples$group))), warnings = nrow(object@diagnostics$warnings)),
    specification = object@specification, design = object@model$design,
    batches = object@diagnostics$batches, gene_outcomes = object@diagnostics$outcomes,
    ql = object@diagnostics$ql, shrinkage = shrinkage,
    dispersion = dispersion, input_actions = object@diagnostics$input_actions,
    warnings = object@diagnostics$warnings, timings = object@diagnostics$timing)
}
