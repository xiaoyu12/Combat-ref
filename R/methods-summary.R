#' @export
S7::method(summary, CombatRefQLFit) <- function(object, ...) {
  CombatRefQLSummary(
    overview = data.frame(genes = nrow(object@counts), samples = ncol(object@counts), batches = length(unique(object@samples$batch)), groups = length(unique(stats::na.omit(object@samples$group))), warnings = nrow(object@diagnostics$warnings)),
    specification = object@specification, design = object@model$design,
    batches = object@diagnostics$batches, gene_outcomes = object@diagnostics$outcomes,
    ql = object@diagnostics$ql, input_actions = object@diagnostics$input_actions,
    warnings = object@diagnostics$warnings, timings = object@diagnostics$timing)
}
